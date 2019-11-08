import sys
import os.path
import numpy as np
import math
import time
import warnings
import pygplates as pgp
import geoTools

from optapm import ModelSetup as ms, ProcessResults as pr
from functools import partial
import itertools
from datetime import datetime, timedelta
from optimised_rotation_updater import OptimisedRotationUpdater

# All the config parameters are now in a separate module 'Optimised_config' that also
# gets imported into the pre-processing modules.
from Optimised_config import *


def warning_format(message, category, filename, lineno, file=None, line=None):
    # return '{0}:{1}: {1}:{1}\n'.format(filename, lineno, category.__name__, message)
    return '{0}: {1}\n'.format(category.__name__, message)
# Print the warnings without the filename and line number.
# Users are not going to want to see that.
warnings.formatwarning = warning_format
# Always print warnings (not just the first time encountered at a particular location).
warnings.simplefilter("always")


# Check the required pygplates version.
#
# PyGPlates version 19 can close the gaps in resolved topologies in the *deforming* model (along deforming lines).
# PyGPlates version 22 can handle topological lines (can get their sub-sub-segment plate IDs).
# PyGPlates version 25 greatly improved the speed of pygplates.RotationModel.
#
# We really do need the speed afforded by pyGPlates version 25, so we'll make that a requirement,
# otherwise it can take about 8 times longer (eg, 16 hours instead of 2 hours).
PYGPLATES_VERSION_REQUIRED = pgp.Version(25)
# Check the imported pygplates version.
if not hasattr(pgp, 'Version') or pgp.Version.get_imported_version() < PYGPLATES_VERSION_REQUIRED:
    raise RuntimeError('{0}: Error - imported pygplates version {1} but version {2} or greater is required'.format(
            os.path.basename(__file__), pgp.Version.get_imported_version(), PYGPLATES_VERSION_REQUIRED))


if __name__ == '__main__':

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
        try:

            rc = ipyparallel.Client(profile='default')
            print "Cores started: ", len(rc.ids)

        except Exception as e:

            print ""
            print "! Caught exception: ", e
            sys.exit(1)

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


    age_range = np.arange(end_age + interval, start_age + interval, interval)

    # The main data directory is the 'data' sub-directory of the directory containing this source file.
    datadir = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'data', '')

    # When using mpi4py we only print and collect/process results in one process (the one with rank/ID 0).
    if use_parallel != MPI4PY or mpi_rank == 0:
        
        # Manages updates to the rotation model due to optimisation.
        #
        # The creation/construction of this OptimisedRotationUpdater object also:
        #   Creates a single optimised rotation file by combining all unoptimised (input) rotations.
        #   The 005-000 rotation feature is inserted (or replaced if already existing in input) and
        #   defined such that the rotation of reference plate (obtained for each time using
        #   'get_reference_params') relative to 000 is zero for each time in 'age_range'.
        optimised_rotation_updater = OptimisedRotationUpdater(
                datadir,
                unoptimised_rotation_filenames,
                age_range,
                get_reference_params,
                data_model,
                model_name)
        
        # The filename of the single optimised rotation file just created
        # (relative to the 'data/' directory).
        rotfile = optimised_rotation_updater.get_optimised_rotation_filename()
        
        
        print "Rotation file to be used: ", rotfile
        print "TM data:", tm_data_type
        print "TM method:", tm_method
        print "Age range for model:", age_range

        print "-------------------------------------------------------------------"
        print ""
        print model_name
        print ""

        print "Search type:", search
        print "Search radius:", search_radius
        print ""
        
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
    
    for i in xrange(0, len(age_range)):
        
        ref_rotation_start_age = age_range[i]
        ref_rotation_end_age = ref_rotation_start_age - interval
        #ref_rotation_end_age = 0.
        
        
        # Determine which components are enabled and their weightings (which could vary over time).
        fracture_zones, fracture_zone_weight = get_fracture_zone_params(ref_rotation_start_age)
        net_rotation, net_rotation_weight = get_net_rotation_params(ref_rotation_start_age)
        trench_migration, trench_migration_weight = get_trench_migration_params(ref_rotation_start_age)
        hotspot_trails, hotspot_trails_weight = get_hotspot_trail_params(ref_rotation_start_age)
        
        # Determine reference plate ID and PMAG rotation file (which could vary over time).
        ref_rotation_plate_id, pmag_rotfile = get_reference_params(ref_rotation_start_age)
        
        # When using mpi4py we only prepare the data in one process (the one with rank/ID 0).
        if use_parallel != MPI4PY or mpi_rank == 0:
            
            print "Start age:", ref_rotation_start_age, "Ma"
            print ""
            
            current_search_radius = search_radius
            current_models = models
            if expand_search_radius_on_ref_plate_switches and i > 0 and search == 'Initial':
                # If the reference plate ID used in this iteration differs from the last iteration then temporarily
                # expand the search diameter to 90 degrees since the two reference plate poles might differ a lot.
                # Is 90 as high as we can go?
                last_ref_rotation_plate_id, _ = get_reference_params(age_range[i-1])
                if ref_rotation_plate_id != last_ref_rotation_plate_id:
                    current_search_radius = 90
                    # Expand number of models by the increase in area of small circle search radius 2*PI*(1 - cos(small_circle_radius)).
                    # NOTE: The values given in search radius appear to be diameters (so halving them).
                    current_models = int(
                        (1.0 - math.cos(0.5 * math.radians(current_search_radius))) /
                        (1.0 - math.cos(0.5 * math.radians(search_radius)))
                        * models + 0.5)
                    print "Temporarily expanding search diameter to {0} from {1} at {2}Ma due to change in reference plate.".format(
                        current_search_radius, search_radius, ref_rotation_start_age)
                    print "Also proportionately expanding number of models to {0} from {1}.".format(current_models, models)
                    print ""

            # --------------------------------------------------------------------

            # Gather parameters
            params = [current_search_radius, rotation_uncertainty, search_type, current_models, model_stop_condition, max_iter,
                      ref_rotation_plate_id, ref_rotation_start_age, ref_rotation_end_age, interpolation_resolution, 
                      rotation_age_of_interest, fracture_zones, net_rotation, trench_migration, hotspot_trails,
                      no_auto_ref_rot_longitude, no_auto_ref_rot_latitude, no_auto_ref_rot_angle, auto_calc_ref_pole, search, 
                      fracture_zone_weight, net_rotation_weight, trench_migration_weight, hotspot_trails_weight,
                      include_chains, interpolated_hotspot_trails, tm_method]

            # --------------------------------------------------------------------

            # Load all data
            data = ms.dataLoader(datadir, rotfile, pmag_rotfile, nnr_rotfile=nnr_rotfile, nnr_relative_datadir=nnr_relative_datadir, 
                                 ridge_file=ridge_file, isochron_file=isochron_file, isocob_file=isocob_file, 
                                 hst_file=hst_file, hs_file=hs_file, interpolated_hotspots=interpolated_hotspots)

            # Calculate starting conditions
            startingConditions = ms.modelStartConditions(params, data, plot)
        
        
        if use_parallel == MPI4PY:
            
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
                    for mpi_index in xrange(mpi_size):
                        # Each process gets the next 'num_x_per_rank' x values.
                        x_index = mpi_index * num_x_per_rank
                        new_x_list.append(xStartingCondition[x_index : x_index + num_x_per_rank])
                    # Distribute any remaining x values (if any) across the first few processes.
                    for x_index in xrange(mpi_size * num_x_per_rank, len(xStartingCondition)):
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
            data_array,
            nnr_datadir, no_net_rotation_file, reformArray, trail_data,
            start_seeds, rotation_age_of_interest_age, data_array_labels_short,
            ref_rot_longitude, ref_rot_latitude, ref_rot_angle,
            seed_lons, seed_lats) = startingConditions[:31]

        if auto_calc_ref_pole == False:

            ref_rot_longitude = no_auto_ref_rot_longitude
            ref_rot_latitude = no_auto_ref_rot_latitude
            ref_rot_angle = no_auto_ref_rot_angle


        # When using mpi4py we only print in one process (the one with rank/ID 0).
        if use_parallel != MPI4PY or mpi_rank == 0:
            
            #print "Number of start seeds generated:", len(start_seeds)
            print "Optimised models to be run:", len(start_seeds)
            print " "
            
            # Flush the print statements (for parallel code).
            sys.stdout.flush()


        # --------------------------------------------------------------------
        # --------------------------------------------------------------------
        # Function to run optimisation routine

        def run_optimisation(x, opt_n, N, lb, ub, model_stop_condition, max_iter, interval, rotation_file, 
                             no_net_rotation_file, ref_rotation_start_age, Lats, Lons, spreading_directions, 
                             spreading_asymmetries, seafloor_ages, PID, CPID, data_array, nnr_datadir, 
                             ref_rotation_end_age, ref_rotation_plate_id, reformArray, trail_data,
                             fracture_zone_weight, net_rotation_weight, trench_migration_weight, hotspot_trails_weight,
                             use_trail_age_uncertainty, trail_age_uncertainty_ellipse, tm_method):

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
            import nlopt
            
            # Load up the object function object once (eg, load rotation files).
            # NLopt will then call it multiple times.
            # NLopt will call this as 'obj_f(x, grad)' because 'obj_f' has a '__call__' method.
            obj_f = ObjectiveFunction(
                    interval, rotation_file, no_net_rotation_file, ref_rotation_start_age, Lats, Lons, spreading_directions,
                    spreading_asymmetries, seafloor_ages, PID, CPID, data_array, nnr_datadir,
                    ref_rotation_end_age, ref_rotation_plate_id, reformArray, trail_data,
                    fracture_zone_weight, net_rotation_weight, trench_migration_weight, hotspot_trails_weight,
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
            print opt.get_numevals()
            sys.stdout.flush()

            return xopt, minf


        # --------------------------------------------------------------------
        # --------------------------------------------------------------------
        # Start optimisation

        # Wrap 'run_optimisation()' by passing all the constant parameters (ie, everything except 'x').
        runopt = partial(run_optimisation, opt_n=opt_n, N=N, lb=lb, ub=ub, 
                          model_stop_condition=model_stop_condition, max_iter=max_iter, interval=interval, rotation_file=rotation_file,
                          no_net_rotation_file=no_net_rotation_file, ref_rotation_start_age=ref_rotation_start_age, 
                          Lats=Lats, Lons=Lons, spreading_directions=spreading_directions, 
                          spreading_asymmetries=spreading_asymmetries, 
                          seafloor_ages=seafloor_ages, PID=PID, CPID=CPID, data_array=data_array, nnr_datadir=nnr_datadir,
                          ref_rotation_end_age=ref_rotation_end_age, ref_rotation_plate_id=ref_rotation_plate_id,
                          reformArray=reformArray, trail_data=trail_data, fracture_zone_weight=fracture_zone_weight,
                          net_rotation_weight=net_rotation_weight, trench_migration_weight=trench_migration_weight,
                          hotspot_trails_weight=hotspot_trails_weight, use_trail_age_uncertainty=use_trail_age_uncertainty,
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
                
                # print 'all xopt', xopt
                
        else:
            
            # Calculate serially.
            xopt = [runopt(x_item) for x_item in x]

        # print 'Mean opt_ind_data', np.mean(objective_function.opt_ind_data, axis=0)
        # print 'Median opt_ind_data', np.median(objective_function.opt_ind_data, axis=0)

        # except Exception as e:
        
        #     text_file = open("Output.txt", "w")
        #     text_file.write("Model error: " + str(e))
        #     text_file.close()


        # When using mpi4py we only collect and process results in one process (the one with rank/ID 0).
        if use_parallel != MPI4PY or mpi_rank == 0:
            
            # Find minimum result from all models
            results = []

            for i in xrange(0, len(xopt)):

                results.append(xopt[i][1])

            min_result_index = np.where(results == np.min(results))[0][0]
            min_result = xopt[min_result_index]

            print " "
            print "Optimisation complete."
            print "Models produced:", len(xopt)
            print " "


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


            for j in xrange(0, len(xopt)):

                costs.append(xopt[j][1])


            rmin = np.array(rmin)
            rmean = np.array(rmean)

            min_results.append(rmin[0])
            mean_results.append(rmean[0])

            plat, plon = geoTools.checkLatLon(min_results[-1][2], min_results[-1][3])


            # end_time = round(time.time() - start, 2)
            # sec = timedelta(seconds = float(end_time))
            # dt = datetime(1,1,1) + sec
            
            # print "Timestep completed in:"
            # print str(dt.day-1) + "d, " + str(dt.hour) + "h, " + str(dt.minute) + "m, " + str(dt.second) + "s."


            # --------------------------------------------------------------------
            # --------------------------------------------------------------------
            # Update the optimised rotation file with result

            optimised_rotation_ref_plate_rel_000 = pgp.FiniteRotation(
                    (np.double(round(plat, 2)), np.double(round(plon, 2))),
                    np.radians(np.double(round(min_results[-1][1], 2))))
            
            optimised_rotation_updater.update_optimised_rotation(
                    optimised_rotation_ref_plate_rel_000,
                    ref_rotation_plate_id,
                    ref_rotation_start_age)
            
            # Flush the print statements (for parallel code).
            sys.stdout.flush()


    # When using mpi4py we only collect and process results in one process (the one with rank/ID 0).
    if use_parallel != MPI4PY or mpi_rank == 0:
        
        main_end_time = round(time.time() - main_start, 10)
        main_sec = timedelta(seconds = float(main_end_time))
        main_dt = datetime(1,1,1) + main_sec

        print ""
        print ""
        print "Model completed in:"
        print str(main_dt.day-1) + "d, " + str(main_dt.hour) + "h, " + str(main_dt.minute) + "m, " + str(main_dt.second) + "s."


        # Scaling (mean of 0-50Ma - 20 models)
        # 
        #     NR: 7:574, 3:465
        # 
        #     TM: 334
        # 
        #     HS: 398

        # Display result arrays
        print np.mean(costs)

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
