import geoTools
import marshal
import math
import numpy as np
import os.path
import pygplates as pgp
import sys
import time

from optapm import ModelSetup as ms, ProcessResults as pr
from functools import partial
import itertools
from datetime import datetime, timedelta
from no_net_rotation_model import NoNetRotationModel
from optimised_rotation_updater import OptimisedRotationUpdater
from plate_velocity_partitioner import PlateVelocityPartitioner
from continent_fragmentation import ContinentFragmentation
from trench_resolver import TrenchResolver

# All the config parameters are now in a separate module 'Optimised_config' that also
# gets imported into the pre-processing modules.
from Optimised_config import *


# Check the required pygplates version.
#
# PyGPlates version 0.19 can close the gaps in resolved topologies in the *deforming* model (along deforming lines).
# PyGPlates version 0.22 can handle topological lines (can get their sub-sub-segment plate IDs).
# PyGPlates version 0.25 greatly improved the speed of pygplates.RotationModel (workflow time reduced by about a factor of 8).
# PyGPlates version 1.0 supports calculating net rotation, and pickling (eg, can pickle a pygplates.RotationModel).
#
PYGPLATES_VERSION_REQUIRED = pgp.Version(1, 0)
# Check the imported pygplates version.
if not hasattr(pgp, 'Version') or pgp.Version.get_imported_version() < PYGPLATES_VERSION_REQUIRED:
    raise RuntimeError('{0}: Error - imported pygplates version {1} but version {2} or greater is required'.format(
            os.path.basename(__file__), pgp.Version.get_imported_version(), PYGPLATES_VERSION_REQUIRED))


if __name__ == '__main__':

    try:

        if use_parallel == IPYPARALLEL:

            import ipyparallel
            
            # Launch ipyparallel client.
            #
            # Can start the engines using:
            #
            #    ipcluster start -n 4
            #
            # ...in this case 4 engines/cores.
            #
            # Alternatively can start cluster in "IPython Clusters" tab of Jupyter notebook and
            # then call 
            rc = ipyparallel.Client(profile='default')
            print("Cores started: ", len(rc.ids))

            dview = rc[:]
            dview.block = True

            # UPDATE: This is now handled by importing inside the 'run_optimisation()' function.
            #
            # Make sure remote engines also import these modules.
            # This is because these modules are referenced in the 'run_optimisation()' function
            # which is executed on the remote engines (see dview.map() below), but
            # the remote engines don't import any import statements (outside 'run_optimisation()').
            #
            # with dview.sync_imports():
            #     from objective_function import ObjectiveFunction
            #     import nlopt
        
        elif use_parallel == MPI4PY:
        
            from mpi4py import MPI
            
            # It seems that if one process/rank raises an exception then we need to manually
            # kill the other MPI processes according to:
            #
            #   https://groups.google.com/forum/#!topic/mpi4py/RovYzJ8qkbc
            #
            # ...otherwise MPI Finalize (in the process that raised exception) will block waiting for
            # the other processes to finish, but they're waiting for input (gather) from the rank=0 process
            # resulting in a deadlock.
            #
            # This code was obtained from:
            #
            #   https://groups.google.com/forum/#!topic/mpi4py/ktAZWIfx8zI
            #
            # ...and is the easiest way to do this if we don't care about properly cleaning up the processes.
            #
            _excepthook = sys.excepthook
            def excepthook(t,v,tb):
                _excepthook(t,v,tb)
                if (not MPI.Is_finalized()
                    and MPI.Is_initialized()):
                    MPI.COMM_WORLD.Abort(1)
            sys.excepthook = excepthook
            
            mpi_comm = MPI.COMM_WORLD
            mpi_size = mpi_comm.Get_size()
            mpi_rank = mpi_comm.Get_rank()
        
        # else serial

        # Create 'optimisation' sub-directory of data model directory (if it doesn't already exist).
        optimisation_sub_dir = os.path.join(datadir, data_model, 'optimisation')
        if not os.path.exists(optimisation_sub_dir):
            print('{} does not exist, creating now... '.format(optimisation_sub_dir))
            os.mkdir(optimisation_sub_dir)

        # The age range to optimise over.
        #
        # Note: 'end_age' will equal 'actual_end_age' unless we are continuing an interrupted run.
        age_range = range(end_age + interval, start_age + interval, interval)

        # When using mpi4py we only print and collect/process results in one process (the one with rank/ID 0).
        if use_parallel != MPI4PY or mpi_rank == 0:

            # Load the topology features. They can take a long time to load (especially for a deforming model) so we
            # do it once instead of three times (once each for no-net-rotation, trench resolving and plate velocities).
            topology_features = []
            for topology_filename in topology_filenames:
                topology_features.extend(pgp.FeatureCollection(os.path.join(datadir, topology_filename)))
            
            # Manages updates to the rotation model due to optimisation.
            #
            # The creation/construction of this OptimisedRotationUpdater object also:
            #   Creates a single optimised rotation file by combining all unoptimised (input) rotations.
            #   The 005-000 rotation feature is inserted (or replaced if already existing in input) and
            #   defined such that the rotation of reference plate (obtained for each time using 'get_reference_params')
            #   relative to 000 is zero for each time from 'start_age' to 'end_age + interval'
            #   in 'interval' steps.
            optimised_rotation_updater = OptimisedRotationUpdater(
                    datadir,
                    original_rotation_filenames,
                    start_age,
                    end_age,
                    actual_end_age,
                    interval,
                    get_reference_params,
                    data_model,
                    model_name)
            
            # The filename of the single optimised rotation file just created
            # (relative to the 'data/' directory).
            rotfile = optimised_rotation_updater.get_optimised_rotation_filename()
            
            # Creates the no-net-rotation model.
            no_net_rotation_model = NoNetRotationModel(
                    datadir,
                    original_rotation_filenames,
                    topology_features,
                    start_age,
                    end_age,
                    actual_end_age,
                    data_model,
                    model_name)
            
            # The filename of single rotation file (containing entire rotation model) with
            # no net rotation.
            nnr_rotfile = no_net_rotation_model.get_no_net_rotation_filename()
            
            # Generates resolved trench features at each reconstruction time.
            trench_resolver = TrenchResolver(
                    datadir,
                    original_rotation_filenames,
                    topology_features,
                    data_model)
            
            # The filename used to store the trench features at the reconstruction time.
            # The same filename is used for all reconstruction times (it just gets overwritten at each time).
            tm_file = trench_resolver.get_trench_migration_filename()
            
            if plate_velocity_continental_polygons_file:
                plate_velocity_plate_features = list(
                    pgp.FeatureCollection(os.path.join(datadir, plate_velocity_continental_polygons_file)))
                plate_velocity_features_are_topologies = False

                # Continental fragmentation (global perimeter-to-area ratio) will be used to adjust the plate velocities weight.
                plate_velocity_fragmentation = ContinentFragmentation(
                        datadir,
                        original_rotation_filenames,
                        plate_velocity_plate_features,
                        plate_velocity_continental_fragmentation_point_spacing_degrees,
                        plate_velocity_continental_fragmentation_area_threshold_steradians,
                        plate_velocity_continental_fragmentation_gap_threshold_radians,
                        age_range)
            else:
                plate_velocity_plate_features = topology_features
                plate_velocity_features_are_topologies = True
                plate_velocity_fragmentation = None
            
            # Generates points and associated plate IDs at each reconstruction time (for plate velocities later on).
            # Use either topologies (continental+ocean) or just continental.
            plate_velocity_partitioner = PlateVelocityPartitioner(
                    datadir,
                    original_rotation_filenames,
                    plate_velocity_plate_features,
                    plate_velocity_features_are_topologies,
                    plate_velocity_fragmentation,
                    data_model,
                    plate_velocity_grid_spacing)
            # The filename used to store the (plate velocity) points and associated plate IDs at the reconstruction time.
            # The same filename is used for all reconstruction times (it just gets overwritten at each time).
            pv_file = plate_velocity_partitioner.get_plate_velocity_filename()

            
            print("Rotation file to be used: ", rotfile)
            print("TM data:", tm_data_type)
            print("TM method:", tm_method)
            print("Age range for model:", age_range)
            print("-------------------------------------------------------------------")
            print("")
            print(model_name)
            print("")

            print("Search type:", search)
            print("Search radius:", search_radius)
            print("")
            
            # Flush the print statements (for parallel code).
            sys.stdout.flush()



            min_results = []
            mean_results = []

            costs = []

            # Start timer over all time steps.
            main_start = time.time()
        
        
        # # This is probably not needed but make sure the rotation file has been written
        # # by the rank 0 process above before other rank processes continue.
        # # It's probably not needed because the first part of each iteration of time loop below does
        # # a scatter/broadcast which also synchronises processes before rotation file is read.
        # if use_parallel == MPI4PY:
        #     mpi_comm.barrier()
        
        
        #
        # Loop through all times.
        #
        
        for i in range(0, len(age_range)):
            
            ref_rotation_start_age = age_range[i]
            ref_rotation_end_age = ref_rotation_start_age - interval
            #ref_rotation_end_age = 0.
            
            
            # When using mpi4py we only prepare the data in one process (the one with rank/ID 0).
            if use_parallel != MPI4PY or mpi_rank == 0:
                
                print("Start age:", ref_rotation_start_age, "Ma")
                print("")
                
                # Incrementally build the no-net-rotation model as we go.
                # The results are updated to the file 'nnr_rotfile'.
                # NOTE: This does nothing if the entire no-net-rotation model was created in 'no_net_rotation_model.__init__()'.
                no_net_rotation_model.update_no_net_rotation(ref_rotation_start_age)
                
                # Generate the resolved trenches at time 'ref_rotation_start_age'.
                # The results are saved to the file 'tm_file'.
                # Note: The file only contains resolved trenches at time 'ref_rotation_start_age'.
                trench_resolver.generate_resolved_trenches(ref_rotation_start_age)
                
                # Generate the (plate velocity) points and associated plate IDs at time 'ref_rotation_start_age'.
                # The results are saved to the file 'pv_file'.
                # Note: The file only contains points and plate IDs partitioned at time 'ref_rotation_start_age'.
                plate_velocity_partitioner.generate_points_and_plate_ids(ref_rotation_start_age)
                
                # Determine reference plate ID (which could vary over time) and reference rotation file.
                ref_rotation_plate_id, ref_rotation_file = get_reference_params(ref_rotation_start_age)
                
                #
                # Testing getting reference rotation from no-net-rotation versus the previous optimised interval versus rotation file.
                #
                if ref_rotation_file == USE_NNR_REFERENCE_FRAME:
                    # No-net-rotation...
                    # If a reference rotation file is not provided then default to using no-net-rotation model.
                    ref_rotation_file = no_net_rotation_model.get_no_net_rotation_filename()
                elif ref_rotation_file == USE_OPTIMISED_REFERENCE_FRAME:
                    # Previous optimised interval...
                    # If a reference rotation file is not provided then default to using reference plate rotation from previous optimisation interval.
                    ref_rotation_file = rotfile

                # Ensure the optimised rotation file has valid rotations from start to end of current interval by
                # re-using the absolute optimisation from start of previous interval (end of current interval).
                # Once we've optimised the current interval we'll overwrite it, but it can get used before then
                # so it should have a reasonable value.
                #
                #   R(0->ts,000->ref_plate) = R(0->te,000->005) * R(0->ts,005->ref_plate)
                #
                _rotation_model = pgp.RotationModel(os.path.join(datadir, rotfile))
                plate_rotation_005_rel_000 = _rotation_model.get_rotation(
                        ref_rotation_end_age, 5, fixed_plate_id=0)
                plate_rotation_ref_plate_rel_005 = _rotation_model.get_rotation(
                        ref_rotation_start_age, ref_rotation_plate_id, fixed_plate_id=5)
                plate_rotation_ref_plate_rel_000 = plate_rotation_005_rel_000 * plate_rotation_ref_plate_rel_005
                optimised_rotation_updater.update_optimised_rotation(
                        plate_rotation_ref_plate_rel_000,
                        ref_rotation_plate_id,
                        ref_rotation_start_age)
                
                current_search_radius = search_radius
                current_models = models
                if expand_search_radius_on_ref_plate_switches and i > 0 and search == 'Initial':
                    # If the reference plate ID used in this iteration differs from the last iteration then temporarily
                    # expand the search diameter to 180 degrees since the two reference plate poles might differ a lot.
                    last_ref_rotation_plate_id, _ = get_reference_params(age_range[i-1])
                    if ref_rotation_plate_id != last_ref_rotation_plate_id:
                        current_search_radius = 180
                        # Expand number of models by the increase in area of small circle search radius 2*PI*(1 - cos(small_circle_radius)).
                        current_models = int(
                            (1.0 - math.cos(math.radians(current_search_radius))) /
                            (1.0 - math.cos(math.radians(search_radius)))
                            * models + 0.5)
                        print("Temporarily expanding search diameter to {0} from {1} at {2}Ma due to change in reference plate.".format(
                            current_search_radius, search_radius, ref_rotation_start_age))
                        print("Also proportionately expanding number of models to {0} from {1}.".format(current_models, models))
                        print("")

                # --------------------------------------------------------------------
                
                # Determine which components are enabled and their weightings (which could vary over time).
                enable_fracture_zones, fracture_zone_weight, fracture_zone_cost_func, fracture_zone_bounds = get_fracture_zone_params(ref_rotation_start_age)
                enable_net_rotation, net_rotation_weight, net_rotation_cost_func, net_rotation_bounds = get_net_rotation_params(ref_rotation_start_age)
                enable_trench_migration, trench_migration_weight, trench_migration_cost_func, trench_migration_bounds = get_trench_migration_params(ref_rotation_start_age)
                enable_hotspot_trails, hotspot_trails_weight, hotspot_trails_cost_func, hotspot_trails_bounds = get_hotspot_trail_params(ref_rotation_start_age)
                enable_plate_velocity, plate_velocity_weight, plate_velocity_cost_func, plate_velocity_bounds = get_plate_velocity_params(ref_rotation_start_age)

                # Gather parameters
                params = [current_search_radius, rotation_uncertainty, search_type, current_models, model_stop_condition, max_iter,
                        ref_rotation_plate_id, ref_rotation_start_age, ref_rotation_end_age, interpolation_resolution, rotation_age_of_interest,
                        enable_fracture_zones, enable_net_rotation, enable_trench_migration, enable_hotspot_trails, enable_plate_velocity,
                        fracture_zone_weight, net_rotation_weight, trench_migration_weight, hotspot_trails_weight, plate_velocity_weight,
                        fracture_zone_cost_func, net_rotation_cost_func, trench_migration_cost_func, hotspot_trails_cost_func, plate_velocity_cost_func,
                        fracture_zone_bounds, net_rotation_bounds, trench_migration_bounds, hotspot_trails_bounds, plate_velocity_bounds,
                        no_auto_ref_rot_longitude, no_auto_ref_rot_latitude, no_auto_ref_rot_angle, auto_calc_ref_pole, search, 
                        include_chains, interpolated_hotspot_trails, tm_method]

                # --------------------------------------------------------------------

                # Load all data
                data = ms.dataLoader(datadir, rotfile, ref_rotation_file, tm_file=tm_file, pv_file=pv_file, nnr_rotfile=nnr_rotfile, 
                                    ridge_file=ridge_file, isochron_file=isochron_file, isocob_file=isocob_file, 
                                    hst_file=hst_file, hs_file=hs_file, interpolated_hotspots=interpolated_hotspots)

                # Calculate starting conditions
                startingConditions = ms.modelStartConditions(params, data, plot)

                # Marshal each cost function into a code string so we can pass it over the network.
                cost_func_array = startingConditions[20]
                cost_func_code_string_array = [marshal.dumps(cost_func.__code__ if sys.version_info[0] >= 3 else cost_func.func_code)  # Python 2 vs 3.
                    for cost_func in cost_func_array]
                startingConditions[20] = cost_func_code_string_array
            
            
            if use_parallel == MPI4PY:
                
                # This is probably not needed but make sure the file (eg, optimised rotation, no-net rotations, trench migration
                # and plate velocity grid) have been written by the rank 0 process above before other rank processes continue.
                #
                # It's probably not needed because the scatter/broadcast just below should also synchronise all processes.
                #
                # But we did get an error opening the rotation file for reading by one of the processes in the objective function.
                mpi_comm.barrier()
                
                if mpi_rank == 0:
                    
                    # print 'all startingConditions[0]', startingConditions[0]
                    
                    # Divide the starting condition into two variables since we'll send them differently (to other processes).
                    xStartingCondition = startingConditions[0]  # this is a list of x
                    constantStartingConditions = startingConditions[1:]
                    
                    # If there are fewer x values than processes then some processes will get an empty list of x values.
                    if len(xStartingCondition) < mpi_size:
                        # Each process expects a list of x values.
                        xStartingCondition = [[x_item] for x_item in xStartingCondition]
                        # The last few processes get empty lists.
                        xStartingCondition.extend([[]] * (mpi_size - len(xStartingCondition)))
                    else:
                        # Divide the 'x' list among the processes.
                        num_x_per_rank = len(xStartingCondition) // mpi_size
                        new_x_list = []
                        for mpi_index in range(mpi_size):
                            # Each process gets the next 'num_x_per_rank' x values.
                            x_index = mpi_index * num_x_per_rank
                            new_x_list.append(xStartingCondition[x_index : x_index + num_x_per_rank])
                        # Distribute any remaining x values (if any) across the first few processes.
                        for x_index in range(mpi_size * num_x_per_rank, len(xStartingCondition)):
                            new_x_list[x_index - mpi_size * num_x_per_rank].append(xStartingCondition[x_index])
                        
                        xStartingCondition = new_x_list
                    
                else:
                    xStartingCondition = None
                    constantStartingConditions = None
                
                # These starting conditions *vary* across all processes so *scatter* them across all processes (from root process).
                xStartingCondition = mpi_comm.scatter(xStartingCondition, root=0)
                
                # These starting conditions are *constant* across all processes so just need to *broadcast* (from root process).
                constantStartingConditions = mpi_comm.bcast(constantStartingConditions, root=0)
                
                # Join 'x' values for the current process with the constant values back into a single list.
                startingConditions = []
                startingConditions.append(xStartingCondition)
                startingConditions.extend(constantStartingConditions)
            
            
            # Extract variables from starting conditions.
            (x, opt_n, N, lb, ub,  model_stop_condition, max_iter,
                rotation_file,
                ref_rotation_start_age, ref_rotation_end_age,
                ref_rotation_plate_id,
                Lats, Lons,
                spreading_directions, spreading_asymmetries, seafloor_ages,
                PID, CPID,
                data_array, weights_array, cost_func_code_string_array, bounds_array,
                trench_migration_file, plate_velocity_file, no_net_rotation_file, reformArray, trail_data,
                start_seeds, rotation_age_of_interest_age, data_array_labels_short,
                ref_rot_longitude, ref_rot_latitude, ref_rot_angle,
                seed_lons, seed_lats) = startingConditions[:35]

            if auto_calc_ref_pole == False:

                ref_rot_longitude = no_auto_ref_rot_longitude
                ref_rot_latitude = no_auto_ref_rot_latitude
                ref_rot_angle = no_auto_ref_rot_angle


            # When using mpi4py we only print in one process (the one with rank/ID 0).
            if use_parallel != MPI4PY or mpi_rank == 0:
                
                #print "Number of start seeds generated:", len(start_seeds)
                print("Optimised models to be run:", len(start_seeds))
                print(" ")
                
                # Flush the print statements (for parallel code).
                sys.stdout.flush()


            # --------------------------------------------------------------------
            # --------------------------------------------------------------------
            # Function to run optimisation routine

            def run_optimisation(x, opt_n, N, lb, ub, model_stop_condition, max_iter, interval, rotation_file, 
                                no_net_rotation_file, ref_rotation_start_age, Lats, Lons, spreading_directions, 
                                spreading_asymmetries, seafloor_ages, PID, CPID,
                                data_array, weights_array, cost_func_code_string_array, bounds_array,
                                trench_migration_file, plate_velocity_file, ref_rotation_end_age, ref_rotation_plate_id,
                                reformArray, trail_data, use_trail_age_uncertainty, trail_age_uncertainty_ellipse, tm_method):

                # Make sure remote nodes/cores also import these modules (when running code in parallel).
                #
                # Since this function we're in (ie, 'run_optimisation()') is executed on remote nodes/cores
                # (when running code in parallel), some parallelisation techniques (eg, ipyparallel) do not
                # process any import statements outside this function on the remote cores. Thus if we had
                # instead placed these import statements at the top of this file we could get 'ImportError's.
                #
                # We only need to import those modules explicitly referenced in this function.
                # For example, the 'objective_function' module will in turn import what it needs (so we don't have to).
                from objective_function import ObjectiveFunction
                import marshal
                import nlopt
                import sys
                import types

                # Turn cost function code strings back into functions.
                cost_func_array = [types.FunctionType(marshal.loads(cost_func_code_string), globals(), 'cost_func')
                                for cost_func_code_string in cost_func_code_string_array]

                # Load up the object function object once (eg, load rotation files).
                # NLopt will then call it multiple times.
                # NLopt will call this as 'obj_f(x, grad)' because 'obj_f' has a '__call__' method.
                obj_f = ObjectiveFunction(
                        interval, rotation_file, no_net_rotation_file, ref_rotation_start_age, Lats, Lons, spreading_directions,
                        spreading_asymmetries, seafloor_ages, PID, CPID, data_array, weights_array, cost_func_array, bounds_array,
                        trench_migration_file, plate_velocity_file, ref_rotation_end_age, ref_rotation_plate_id, reformArray, trail_data,
                        use_trail_age_uncertainty, trail_age_uncertainty_ellipse, tm_method)
                
                opt = nlopt.opt(nlopt.LN_COBYLA, opt_n)
                opt.set_min_objective(obj_f)
                opt.set_lower_bounds(lb)
                opt.set_upper_bounds(ub)

                # Select model stop condition
                if model_stop_condition != 'threshold':

                    opt.set_maxeval(max_iter)

                else:

                    opt.set_ftol_rel(1e-6)
                    opt.set_xtol_rel(1e-8)

                xopt = opt.optimize(x)
                minf = opt.last_optimum_value()
                
                # Debug print number of iterations needed to converge.
                #print(opt.get_numevals(), 'iterations performed')
                #sys.stdout.flush()

                # To debug the weighted cost functions (net rotation, trench migration, etc).
                #print('Min objective function costs', np.min(obj_f.debug_data_array, axis=0))
                #print('Max objective function costs', np.max(obj_f.debug_data_array, axis=0))
                #print('Mean objective function costs', np.mean(obj_f.debug_data_array, axis=0))
                #print('Std dev objective function costs', np.std(obj_f.debug_data_array, axis=0))
                #print('Median objective function costs', np.median(obj_f.debug_data_array, axis=0))
                #print('Median abs dev objective function costs', np.median(np.absolute(obj_f.debug_data_array - np.median(obj_f.debug_data_array, axis=0)), axis=0))
                #print('Optimal cost', minf)
                #sys.stdout.flush()

                return xopt, minf


            # --------------------------------------------------------------------
            # --------------------------------------------------------------------
            # Start optimisation

            # Wrap 'run_optimisation()' by passing all the constant parameters (ie, everything except 'x').
            runopt = partial(run_optimisation, opt_n=opt_n, N=N, lb=lb, ub=ub, 
                            model_stop_condition=model_stop_condition, max_iter=max_iter, interval=interval, rotation_file=rotation_file,
                            no_net_rotation_file=no_net_rotation_file, ref_rotation_start_age=ref_rotation_start_age, 
                            Lats=Lats, Lons=Lons, spreading_directions=spreading_directions, spreading_asymmetries=spreading_asymmetries,
                            seafloor_ages=seafloor_ages, PID=PID, CPID=CPID,
                            data_array=data_array, weights_array=weights_array, cost_func_code_string_array=cost_func_code_string_array, bounds_array=bounds_array,
                            trench_migration_file=trench_migration_file, plate_velocity_file=plate_velocity_file,
                            ref_rotation_end_age=ref_rotation_end_age, ref_rotation_plate_id=ref_rotation_plate_id,
                            reformArray=reformArray, trail_data=trail_data, use_trail_age_uncertainty=use_trail_age_uncertainty,
                            trail_age_uncertainty_ellipse=trail_age_uncertainty_ellipse, tm_method=tm_method)

            # Start timer for current time step.
            #start = time.time()

            #
            # Run optimisation in parallel or serial.
            #
            if use_parallel == IPYPARALLEL:
            
                # 'x' is a list, so distribute the elements across the processes.
                xopt = dview.map(runopt, x)
                
            elif use_parallel == MPI4PY:
                
                # print '%d:' % mpi_rank, 'x', x
                
                # Current process runs an optimisation on each element the sub-list it received from the root process.
                #
                # If there's too many processes (ie, not enough tasks to go around) then some processes
                # will have an empty list and hence have nothing to do here.
                xopt = [runopt(x_item) for x_item in x]
                
                # print '%d:' % mpi_rank, 'xopt', xopt
                
                # Gather results from all processes into the root (0) process.
                # Gathers a small list from each process, so root process will end up with a list of lists.
                xopt = mpi_comm.gather(xopt, root=0)
                if mpi_rank == 0:
                    # Flatten a list of lists into a single list.
                    # [[x1, x2], [x3, x4]] -> [x1, x2, x3, x4].
                    # Note that if some processes had no work to do then some lists will be empty, as in...
                    # [[x1], [x2], [x3], [x4], []] -> [[x1], [x2], [x3], [x4]]
                    # ...where there were 5 processes but only 4 'x' values to process.
                    xopt = list(itertools.chain.from_iterable(xopt))
                    
                    # print('all xopt', xopt)
                    # sys.stdout.flush()
                    
            else:
                
                # Calculate serially.
                xopt = [runopt(x_item) for x_item in x]

            # except Exception as e:
            
            #     text_file = open("Output.txt", "w")
            #     text_file.write("Model error: " + str(e))
            #     text_file.close()


            # When using mpi4py we only collect and process results in one process (the one with rank/ID 0).
            if use_parallel != MPI4PY or mpi_rank == 0:
                
                # Find minimum result from all models
                results = []

                for i in range(0, len(xopt)):

                    results.append(xopt[i][1])

                min_result_index = np.where(results == np.min(results))[0][0]
                min_result = xopt[min_result_index]

                print(" ")
                print("Optimisation complete.")
                print("Models produced:", len(xopt))
                #print('Minimum optimal cost:', np.min(results))
                print(" ")


                # Save results to pickle file located as '/model_output/
                output_file = pr.saveResultsToPickle(data_array, data_array_labels_short, ref_rotation_start_age, 
                                                    ref_rotation_end_age, current_search_radius, xopt, current_models, model_name)


                # Plot results
                rmin, rmean = pr.sortAndPlot(output_file, ref_rotation_start_age, ref_rotation_end_age, 
                                            rotation_age_of_interest_age, xopt, rotation_file, ref_rot_longitude,
                                            ref_rot_latitude, ref_rot_angle, seed_lons, seed_lats, 
                                            ref_rotation_plate_id, model_name, current_models, data_array_labels_short, 
                                            data_array, current_search_radius,
                                            plot)


                for j in range(0, len(xopt)):

                    costs.append(xopt[j][1])


                #print('rmin:', rmin)
                min_results.append(np.array(rmin[['Age', 'Ang', 'Lat', 'Lon', 'Minimum', 'Model']])[0])
                mean_results.append(np.array(rmean[['Age', 'Ang', 'Lat', 'Lon', 'Minimum', 'Model']])[0])

                ang = min_results[-1][1]
                plat, plon = geoTools.checkLatLon(min_results[-1][2], min_results[-1][3])
                #print("plat/plon/ang:", plat, plon, ang)


                # end_time = round(time.time() - start, 2)
                # sec = timedelta(seconds = float(end_time))
                # dt = datetime(1,1,1) + sec
                
                # print "Timestep completed in:"
                # print str(dt.day-1) + "d, " + str(dt.hour) + "h, " + str(dt.minute) + "m, " + str(dt.second) + "s."


                # --------------------------------------------------------------------
                # --------------------------------------------------------------------
                # Update the optimised rotation file with result

                optimised_rotation_ref_plate_rel_000 = pgp.FiniteRotation(
                        (plat, plon),
                        math.radians(ang))
                
                optimised_rotation_updater.update_optimised_rotation(
                        optimised_rotation_ref_plate_rel_000,
                        ref_rotation_plate_id,
                        ref_rotation_start_age)
                
                # Flush the print statements (for parallel code).
                sys.stdout.flush()


        # When using mpi4py we only collect and process results in one process (the one with rank/ID 0).
        if use_parallel != MPI4PY or mpi_rank == 0:
            
            # Save the final optimised model back to the original rotation files (or copies of them).
            optimised_rotation_updater.save_to_rotation_files()
            
            main_end_time = round(time.time() - main_start, 10)
            main_sec = timedelta(seconds = float(main_end_time))
            main_dt = datetime(1,1,1) + main_sec

            print("")
            print("")
            print("Model completed in:")
            print(str(main_dt.day-1) + "d, " + str(main_dt.hour) + "h, " + str(main_dt.minute) + "m, " + str(main_dt.second) + "s.")


            # Scaling (mean of 0-50Ma - 20 models)
            # 
            #     NR: 7:574, 3:465
            # 
            #     TM: 334
            # 
            #     HS: 398

            # Display result arrays
            print(np.mean(costs))

            # print "Mean of 20 models (0-50Ma)"
            # print ""
            # print "tm_eval =", 47 * 3
            # print "nr_eval =", 143


            # import pickle

            # with open('model_output/optAPM175_10-0Ma_10models_NR_TM_60.pkl', 'rb') as f:
            #     data = pickle.load(f)
                
            # print data
        
        
        # # This is probably not needed but we're getting garbage written to the rotation output file
        # # for some reason (even though only the rank 0 process writes to the rotation file).
        # #
        # # UPDATE: Problem was caused by pyGPlates Plates4 rotation writer having trouble with Unicode chars.
        # if use_parallel == MPI4PY:
        #     mpi_comm.barrier()

    except Exception as e:

        print("")
        print("! Caught exception: ", e)
        raise
