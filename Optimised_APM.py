import sys
import os.path
import numpy as np
import time
import pygplates as pgp
import geoTools

from optapm import ModelSetup as ms, ProcessResults as pr
from functools import partial
import itertools
from datetime import datetime, timedelta


#
# Supported parallelisation methods (or None to disable parallelisation, eg, for testing).
#
MPI4PY = 0
IPYPARALLEL = 1

# Choose parallelisation method (or None to disable parallelisation, eg, for testing).
use_parallel = MPI4PY


#
# Set model parameters, load data, and calculate starting conditions
# 
# Sets all user-selected parameters for the mode run
# 
# Arguments:
# 
#     geographical_uncertainty : Number that approximately represents geographical uncertainty - 95% confidence limit around ref pole location
#     rotation_uncertainty : Number that represents the upper and lower bounds of the optimisation's angle variation
#     sample_space : Selects the sampling method to be used to generate start seeds. "Fisher" (spherical distribution) by default.
#     models : The total number of models to be produced. 1 model = 1 complete optimisation from 1 starting location
#     model_stop_condition : Type of condition to be used to terminate optimisation. "threshold" or "max_iter". Threshold by default.
#     max_iter : IF "max_iter" selected, sets maximum iterations regardless of successful convergence.
#     ref_rotation_plate_id : Plate to be used as fixed reference. 701 (Africa) by default.
#     ref_rotation_start_age : Rotation begin age.
#     ref_rotation_end_age : Rotation end age.
#     interpolation_resolution : Resolution in degrees used in the Fracture Zone calulcations
#     rotation_age_of_interest : The result we are interested in. Allows for 'windowed mean' approach. It is the midpoint between the start and end ages by default.
# 
# Data to be included in optimisation: True by default.
# 
#     fracture_zones : Boolean
#     net_rotation : Boolean
#     trench_migration : Boolean
#     hotspot_reconstruction : Boolean
#     hotspot_dispersion : Boolean
# 

# model_name = "optAPM176"
model_name = "optAPM9"

start_age = 410
end_age = 0
interval = 10

search = "Initial"
search_radius = 60
rotation_uncertainty = 30
auto_calc_ref_pole = True
models = 100

model_stop_condition = 'threshold'
max_iter = 5  # Only applies if model_stop_condition != 'threshold'

# Trench migration parameters
tm_method = 'pygplates' # 'pygplates' for new method OR 'convergence' for old method
# tm_data_type = 'muller2016' # 'muller2016' or 'shephard2013' or 'Global_Model_WD_Internal_Release_2016_v3'
tm_data_type = 'Global_Model_WD_Internal_Release_2016_v3' # 'muller2016' or 'shephard2013' or 'Global_Model_WD_Internal_Release_2016_v3'

# Hotspot parameters:
interpolated_hotspot_trails = True
use_trail_age_uncertainty = True

# Millions of years - e.g. 2 million years @ 50mm per year = 100km radius uncertainty ellipse
trail_age_uncertainty_ellipse = 1

include_chains = ['Louisville', 'Tristan', 'Reunion', 'St_Helena', 'Foundation', 'Cobb', 'Samoa', 'Tasmantid', 
                  'Hawaii']
#include_chains = ['Louisville', 'Tristan', 'Reunion', 'Hawaii', 'St_Helena', 'Tasmantid']


# Rotation file with existing APM rotations removed from 0-250Ma to be used:
if tm_data_type == 'Global_Model_WD_Internal_Release_2016_v3':

    original_rotfile = 'Global_Model_WD_Internal_Release_2016_v3/optimisation/all_rotations.rot'
    rotfile = 'Global_Model_WD_Internal_Release_2016_v3/optimisation/all_rotations_' + model_name + '.rot'

elif tm_data_type == 'muller2016':

    original_rotfile = 'Global_EarthByte_230-0Ma_GK07_AREPS.rot'
    rotfile = 'Global_EarthByte_230-0Ma_GK07_AREPS_' + model_name + '.rot'

elif tm_data_type == 'shephard2013':

    original_rotfile = 'Shephard_etal_ESR2013_Global_EarthByte_2013.rot'
    rotfile = 'Shephard_etal_ESR2013_Global_EarthByte_2013_' + model_name + '.rot'


# Large area grid search to find minima
if search == 'Initial':

    search_type = 'Random'

# Uses grid search minima as seed for targeted secondary search (optional)
elif search == 'Secondary':

    search_type = 'Uniform'
    search_radius = 15
    rotation_uncertainty = 30

    models = 60
    auto_calc_ref_pole = False

# Used when auto_calc_ref_pole is False.
no_auto_ref_rot_longitude = -53.5
no_auto_ref_rot_latitude = 56.6
no_auto_ref_rot_angle = -2.28

ref_rotation_plate_id = 701

interpolation_resolution = 5
rotation_age_of_interest = True

pmag_rotfile = 'Palaeomagnetic_Africa_S.rot'

if tm_method == 'convergence':

    if tm_data_type == 'muller2016':

        nnr_relative_datadir = 'TMData/'
        nnr_rotfile = 'Global_EarthByte_230-0Ma_GK07_AREPS_NNR.rot'


elif tm_method == 'pygplates':

    if tm_data_type == 'Global_Model_WD_Internal_Release_2016_v3':

        nnr_relative_datadir = 'TMData/Global_Model_WD_Internal_Release_2016_v3/'
        nnr_rotfile = 'Global_Model_WD_Internal_Release_2016_v3/optimisation/all_rotations_NNR.rot'

    elif tm_data_type == 'muller2016':

        nnr_relative_datadir = 'TMData/Muller_2016/'
        nnr_rotfile = 'Global_EarthByte_230-0Ma_GK07_AREPS_NNR.rot'

    elif tm_data_type == 'shephard2013':

        nnr_relative_datadir = 'TMData/Shephard_2013/'
        nnr_rotfile = 'Shephard_etal_ESR2013_Global_EarthByte_NNR_ORIGINAL.rot'


if tm_data_type == 'Global_Model_WD_Internal_Release_2016_v3':

    ridge_file = 'Global_Model_WD_Internal_Release_2016_v3/StaticGeometries/AgeGridInput/Global_EarthByte_GeeK07_Ridges_2016_v3.gpml'
    isochron_file = 'Global_Model_WD_Internal_Release_2016_v3/StaticGeometries/AgeGridInput/Global_EarthByte_GeeK07_Isochrons_2016_v3.gpml'
    isocob_file = 'Global_Model_WD_Internal_Release_2016_v3/StaticGeometries/AgeGridInput/Global_EarthByte_GeeK07_IsoCOB_2016_v3.gpml'

else:

    ridge_file = 'Global_EarthByte_230-0Ma_GK07_AREPS_Ridges.gpml'
    isochron_file = 'Global_EarthByte_230-0Ma_GK07_AREPS_Isochrons.gpmlz'
    isocob_file = 'Global_EarthByte_230-0Ma_GK07_AREPS_IsoCOB.gpml'

hst_file = 'HotspotTrails.geojson'
hs_file = 'HotspotCatalogue2.geojson'
interpolated_hotspots = 'interpolated_hotspot_chains_5Myr.xlsx'


# Don't plot in this workflow.
# This is so it can be run on an HPC cluster with no visualisation node.
plot = False


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

    rotation_file = datadir + rotfile


    # When using mpi4py we only print and collect/process results in one process (the one with rank/ID 0).
    if use_parallel != MPI4PY or mpi_rank == 0:
        
        #
        # Copy the original rotation file into the optimised version, and zero out the rotations
        # that we will soon replace with optimised rotations.
        #
        original_rotation_file = datadir + original_rotfile
        rotation_features = list(pgp.FeatureCollection(original_rotation_file))
        original_rotation_model = pgp.RotationModel(rotation_features)
        
        # for rotation_feature in rotation_features:
        #     
        #     total_reconstruction_pole = rotation_feature.get_total_reconstruction_pole()
        #     if total_reconstruction_pole:
        #         
        #         fixed_plate_id, moving_plate_id, rotation_sequence = total_reconstruction_pole
        #         if moving_plate_id == ref_rotation_plate_id:
        #             
        #             zero_rotation_time_samples_701_rel_fixed = []
        #             
        #             # Start with identity rotation at time 0Ma.
        #             zero_rotation_time_samples_701_rel_fixed.append(
        #                 pgp.GpmlTimeSample(pgp.GpmlFiniteRotation(pgp.FiniteRotation()), 0.0, 'optAPM'))
        #             
        #             for ref_rotation_start_age in age_range:
        #                 
        #                 # We want our 701 to 001 rotation to be zero.
        #                 # However the 701 sequence might have a fixed plate ID that is not 001.
        #                 # So convert zero 701 rel 001 rotation to the 701 rel 'fixed_plate_id' rotation to store in rotation feature.
        #                 #
        #                 #                         R(0->t,001->701) = R(0->t,001->fixed) * R(0->t,fixed->701)
        #                 #                                 Identity = R(0->t,001->fixed) * R(0->t,fixed->701)
        #                 #   inverse(R(0->t,001->fixed)) * Identity = R(0->t,fixed->701)
        #                 #                       R(0->t,fixed->701) = inverse(R(0->t,001->fixed))
        #                 #
        #                 zero_rotation_701_rel_fixed = original_rotation_model.get_rotation(ref_rotation_start_age, fixed_plate_id, fixed_plate_id=1).get_inverse()
        # 
        #                 zero_rotation_time_samples_701_rel_fixed.append(
        #                     pgp.GpmlTimeSample(pgp.GpmlFiniteRotation(zero_rotation_701_rel_fixed), ref_rotation_start_age, 'optAPM'))
        #             
        #             # Replace the 701 rotation sequence.
        #             rotation_feature.set_total_reconstruction_pole(
        #                 fixed_plate_id,
        #                 ref_rotation_plate_id,
        #                 pgp.GpmlIrregularSampling(zero_rotation_time_samples_701_rel_fixed))
        #             
        #             # Write the rotation file with zero 701-to-anchor rotations.
        #             pgp.FeatureCollection(rotation_features).write(rotation_file)
        #             break
        
        for rotation_feature in rotation_features:
            
            total_reconstruction_pole = rotation_feature.get_total_reconstruction_pole()
            if total_reconstruction_pole:
                
                fixed_plate_id, moving_plate_id, rotation_sequence = total_reconstruction_pole
                if moving_plate_id == 5 and fixed_plate_id == 0:
                    
                    zero_rotation_time_samples_005_rel_000 = []
                    
                    # Start with identity rotation at time 0Ma.
                    zero_rotation_time_samples_005_rel_000.append(
                        pgp.GpmlTimeSample(pgp.GpmlFiniteRotation(pgp.FiniteRotation()), 0.0, 'optAPM'))
                    
                    for ref_rotation_start_age in age_range:
                        
                        # We want our 701 to 000 rotation to be zero.
                        #
                        #                         R(0->t,000->701) = R(0->t,000->005) * R(0->t,005->701)
                        #                                 Identity = R(0->t,000->005) * R(0->t,005->701)
                        #   inverse(R(0->t,000->005)) * Identity = R(0->t,005->701)
                        #                       R(0->t,000->005) = inverse(R(0->t,005->701))
                        #
                        zero_rotation_005_rel_000 = original_rotation_model.get_rotation(ref_rotation_start_age, ref_rotation_plate_id, fixed_plate_id=0).get_inverse()

                        zero_rotation_time_samples_005_rel_000.append(
                            pgp.GpmlTimeSample(pgp.GpmlFiniteRotation(zero_rotation_005_rel_000), ref_rotation_start_age, 'optAPM'))
                    
                    # Replace the 005/000 rotation sequence.
                    rotation_feature.set_total_reconstruction_pole(
                        0,
                        5,
                        pgp.GpmlIrregularSampling(zero_rotation_time_samples_005_rel_000))
                    
                    # Write the rotation file with zero 701-to-anchor rotations.
                    pgp.FeatureCollection(rotation_features).write(rotation_file)
                    break
        
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
        
        
        fracture_zones   = False
        net_rotation     = True
        trench_migration = True
        # # Only use hotspot trails for 0-80Ma.
        if ref_rotation_start_age <= 80:
            hotspot_trails   = True
        else:
            hotspot_trails   = False
        
        # # sigma (i.e. cost / sigma = weight)
        fracture_zone_weight    = 1.0
        trench_migration_weight = 1.0
        hotspot_trails_weight   = 1.0
        # Reduce net rotation weight to half for times older than 80Ma
        # due to large number of synthetic (uncertain) plates.
        if ref_rotation_start_age <= 80:
            net_rotation_weight = 1.0
        else:
            # NOTE: These are inverse weights (ie, the constraint costs are *multiplied* by "1.0 / weight").
            net_rotation_weight = 2.0  # Gives a *multiplicative* weight of 0.5
        
        
        # When using mpi4py we only prepare the data in one process (the one with rank/ID 0).
        if use_parallel != MPI4PY or mpi_rank == 0:
            
            print "Start age:", ref_rotation_start_age, "Ma"
            print ""

            # --------------------------------------------------------------------

            # Gather parameters
            params = [search_radius, rotation_uncertainty, search_type, models, model_stop_condition, max_iter,
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

            rotation_model = data[0]

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
            _,
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
        
        
        # Debugging.
        # import objective_function
        # objective_function.opt_ind_data = []


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
                                                 ref_rotation_end_age, search_radius, xopt, models, model_name)


            # Plot results
            rmin, rmean = pr.sortAndPlot(output_file, ref_rotation_start_age, ref_rotation_end_age, 
                                         rotation_age_of_interest_age, xopt, rotation_file, ref_rot_longitude,
                                         ref_rot_latitude, ref_rot_angle, seed_lons, seed_lats, 
                                         ref_rotation_plate_id, model_name, models, data_array_labels_short, 
                                         data_array, search_radius,
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
            # Update rotation file with result

            rotation_model_tmp = pgp.FeatureCollection(rotation_file)

            # Find existing rotations for Africa
            opt_rotation_feature = None
            for rotation_feature in rotation_model_tmp:

                total_reconstruction_pole = rotation_feature.get_total_reconstruction_pole()

                if total_reconstruction_pole:

                    fixed_plate_id, moving_plate_id, rotation_sequence = total_reconstruction_pole

                    # if moving_plate_id == ref_rotation_plate_id:
                    #     opt_rotation_feature = rotation_feature
                    #     break

                    if moving_plate_id == 5 and fixed_plate_id == 0:
                        opt_rotation_feature = rotation_feature
                        break


            # Update existing rotation in the model with result
            if opt_rotation_feature:

                # adjustment_time = pgp.GeoTimeInstant(ref_rotation_start_age)

                for finite_rotation_sample in rotation_sequence.get_enabled_time_samples():

                    finite_rotation_time = finite_rotation_sample.get_time()

                    if finite_rotation_time == ref_rotation_start_age:

                        # new_rotation = pgp.FiniteRotation((np.double(round(min_results[-1][2], 2)), 
                        #                                    np.double(round(min_results[-1][3], 2))), 
                        #                                    np.radians(np.double(round(min_results[-1][1], 2))))

                        new_rotation = pgp.FiniteRotation((np.double(round(plat, 2)), 
                                                           np.double(round(plon, 2))), 
                                                           np.radians(np.double(round(min_results[-1][1], 2))))

                        # # Our new rotation is from 701 to 001 so remove the 'fixed_plate_id' to 001 part to get the
                        # # 701 to 'fixed_plate_id' part that get stored in this rotation feature.
                        # fixed_plate_rotation = rotation_model.get_rotation(ref_rotation_start_age, fixed_plate_id, fixed_plate_id=1)
                        # new_rotation = fixed_plate_rotation.get_inverse() * new_rotation

                        # Our new rotation is from 701 to 000 so remove the 701 to 005 part to get the
                        # 005 to 000 part that get stored in this rotation feature.
                        #
                        #                               R(0->t,000->701) = R(0->t,000->005) * R(0->t,005->701)
                        #   R(0->t,000->701) * inverse(R(0->t,005->701)) = R(0->t,000->005)
                        #
                        plate_rotation_701_rel_005 = rotation_model.get_rotation(
                                ref_rotation_start_age, ref_rotation_plate_id, fixed_plate_id=5)
                        new_rotation = new_rotation * plate_rotation_701_rel_005.get_inverse()
                        
                        finite_rotation_sample.get_value().set_finite_rotation(new_rotation)


            # Add result rotation pole to rotation file
            rotation_model_tmp.write(rotation_file)


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
