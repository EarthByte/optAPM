import numpy as np
import pandas as pd
import time
import pygplates as pgp
import pmagpy.pmag as pmag
import geoTools
import nlopt
import ipyparallel

from optapm import ModelSetup as ms, ProcessResults as pr
from functools import partial
from datetime import datetime, timedelta


# Launch parallel client
try:

    rc = ipyparallel.Client(profile='default')
    print "Cores started: ", len(rc.ids)

except Exception as e:

    print ""
    print "! Caught exception: ", e

print "Cores started: ", len(rc.ids)

dview = rc[:]
dview.block = True  # All calls synchronous (eg, dview.map).


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

model_name = "optAPM175"

start_age = 80
end_age = 0
interval = 10

search = "Initial"
search_radius = 60
rotation_uncertainty = 30
auto_calc_ref_pole = True
models = 2

model_stop_condition = 'threshold'
max_iter = 5  # 250

fracture_zones   = False
net_rotation     = True
trench_migration = True
hotspot_trails   = False

# # sigma (i.e. cost / sigma = weight)
fracture_zone_weight    = 1.
net_rotation_weight     = 1.
trench_migration_weight = 1.
hotspot_trails_weight   = 1.

# Trench migration parameters
tm_method = 'pygplates' # 'pygplates' for new method OR 'convergence' for old method
tm_data_type = 'muller2016' # 'muller2016' or 'shephard2013'

# Hotspot parameters:
interpolated_hotspot_trails = True
use_trail_age_uncertainty = True

# Millions of years - e.g. 2 million years @ 50mm per year = 100km radius uncertainty ellipse
trail_age_uncertainty_ellipse = 1

include_chains = ['Louisville', 'Tristan', 'Reunion', 'St_Helena', 'Foundation', 'Cobb', 'Samoa', 'Tasmantid', 
                  'Hawaii']
#include_chains = ['Louisville', 'Tristan', 'Reunion', 'Hawaii', 'St_Helena', 'Tasmantid']



min_results = []
mean_results = []

costs = []

age_range = np.arange(end_age + interval, start_age + interval, interval)

# Rotation file with existing APM rotations removed from 0-230Ma to be used:
if tm_data_type == 'muller2016':

    rotfile = 'Global_EarthByte_230-0Ma_GK07_AREPS_' + model_name + '.rot'

elif tm_data_type == 'shephard2013':

    rotfile = 'Shephard_etal_ESR2013_Global_EarthByte_2013_' + model_name + '.rot'

print "Rotation file to be used: ", rotfile
print "TM data:", tm_data_type
print "TM method:", tm_method
print "Age range for model:", age_range

print "-------------------------------------------------------------------"
print ""
print model_name
print ""

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

ref_rot_longitude = -53.5
ref_rot_latitude = 56.6
ref_rot_angle = -2.28

ref_rotation_plate_id = 701

interpolation_resolution = 5
rotation_age_of_interest = True

datadir = '/Users/John/Development/Usyd/source_code/other/Artemis/optAPM/data/'

pmag_rotfile = 'Palaeomagnetic_Africa_S.rot'

if tm_method == 'convergence':

    if tm_data_type == 'muller2016':

        nnr_datadir = 'TMData/'
        nnr_rotfile = 'Global_EarthByte_230-0Ma_GK07_AREPS_NNR.rot'


elif tm_method == 'pygplates':

    if tm_data_type == 'muller2016':

        nnr_datadir = 'TMData/Muller_2016/'
        nnr_rotfile = 'Global_EarthByte_230-0Ma_GK07_AREPS_NNR.rot'

    elif tm_data_type == 'shephard2013':

        nnr_datadir = 'TMData/Shephard_2013/'
        nnr_rotfile = 'Shephard_etal_ESR2013_Global_EarthByte_NNR_ORIGINAL.rot'

ridge_file = 'Global_EarthByte_230-0Ma_GK07_AREPS_Ridges.gpml'
isochron_file = 'Global_EarthByte_230-0Ma_GK07_AREPS_Isochrons.gpmlz'
isocob_file = 'Global_EarthByte_230-0Ma_GK07_AREPS_IsoCOB.gpml'
hst_file = 'HotspotTrails.geojson'
hs_file = 'HotspotCatalogue2.geojson'
interpolated_hotspots = 'interpolated_hotspot_chains_5Myr.xlsx'


print "Search type:", search
print "Search radius:", search_radius
print ""


# Loop through all times
for i in xrange(0, len(age_range)):

    # if fracture_zones == True:
    #     if age_range[i] <= 40:
    #         fracture_zones = True
    #     else:
    #         fracture_zones = False

    ref_rotation_start_age = age_range[i]
    ref_rotation_end_age = ref_rotation_start_age - interval
    #ref_rotation_end_age = 0.

    print "Start age:", ref_rotation_start_age, "Ma"
    print ""


    # --------------------------------------------------------------------

    # Gather parameters
    params = [search_radius, rotation_uncertainty, search_type, models, model_stop_condition, max_iter,
              ref_rotation_plate_id, ref_rotation_start_age, ref_rotation_end_age, interpolation_resolution, 
              rotation_age_of_interest, fracture_zones, net_rotation, trench_migration, hotspot_trails,
              ref_rot_longitude, ref_rot_latitude, ref_rot_angle, auto_calc_ref_pole, search, 
              fracture_zone_weight, net_rotation_weight, trench_migration_weight, hotspot_trails_weight,
              include_chains, interpolated_hotspot_trails, tm_method]

    # --------------------------------------------------------------------

    # Load all data
    data = ms.dataLoader(datadir, rotfile, pmag_rotfile, nnr_rotfile=nnr_rotfile, nnr_datadir=nnr_datadir, 
                         ridge_file=ridge_file, isochron_file=isochron_file, isocob_file=isocob_file, 
                         hst_file=hst_file, hs_file=hs_file, interpolated_hotspots=interpolated_hotspots)


    # Calculate starting conditions
    startingConditions = ms.modelStartConditions(params, data)


    # --------------------------------------------------------------------
    # --------------------------------------------------------------------
    # Objective function

    def obj_f(x, grad):

        import numpy as np
        import pygplates as pgp
        import optimisation_methods
        import obj_func_convergence
        import geoTools
        import pmagpy.pmag as pmag
        import subduction_convergence_for_absolute_plate_motion as scap

        from optapm import ObjectiveFunctions
        from scipy import stats

        print "Data array", data_array


        # Prepare rotation model for updates during optimisation - keeps rotations in memory
        rotation_model_tmp = pgp.FeatureCollection(rotation_file)



        #### -----------------------------------------------------------------------------------------
        #### 1. Calculate reconstructed data point locations

        tmp_opt_rlon = []
        tmp_opt_rlat = []
        opt_stats = []

        #import geoTools
        # Check incoming Africa finite rotation pole values
        lat_, lon_ = geoTools.checkLatLon(x[1], x[0])
        ang_ = x[2]


        #### -----------------------------------------------------------------------------------------
        #### 2. Find and update Africa rotation


        # Find existing rotation for Africa at correct time
        opt_rotation_feature = None
        for rotation_feature in rotation_model_tmp:

            total_reconstruction_pole = rotation_feature.get_total_reconstruction_pole()

            if total_reconstruction_pole:

                fixed_plate_id, moving_plate_id, rotation_sequence = total_reconstruction_pole

                if fixed_plate_id == 001 and moving_plate_id == 701:

                    opt_rotation_feature = rotation_feature
                    break


        # Update rotation file with proposed Africa rotation
        if opt_rotation_feature:

            adjustment_time = pgp.GeoTimeInstant(ref_rotation_start_age)

            for finite_rotation_samples in rotation_sequence.get_enabled_time_samples():

                finite_rotation_time = finite_rotation_samples.get_time()

                if finite_rotation_time == ref_rotation_start_age:

                    finite_rotation = finite_rotation_samples.get_value().get_finite_rotation()

                    new_rotation = pgp.FiniteRotation((np.double(lat_), np.double(lon_)), 
                                                      np.radians(np.double(ang_)))
                    finite_rotation_samples.get_value().set_finite_rotation(new_rotation)

        rotation_model_updated = pgp.RotationModel(rotation_model_tmp)




        #### -----------------------------------------------------------------------------------------
        #### 3. Calculate data fits


        #
        # Fracture zone orientation
        if data_array[0] == True:

            # Get skew values
            fz = optimisation_methods.Calc_Median(rotation_model_updated, PID, 
                                                  seafloor_ages, Lats, Lons, 
                                                  spreading_directions)


            tmp_fz_eval = fz[0] + fz[1]

            fz_eval = tmp_fz_eval / fracture_zone_weight



        #
        # Net rotation
        if data_array[1] == True:

            # Prepare no net rotation model for updates during optimisation - keeps rotations in memory
            nn_rotation_model = pgp.RotationModel(no_net_rotation_file)

            nr_timesteps = np.arange(ref_rotation_end_age, ref_rotation_start_age + 1, 2)

            PTLong1, PTLat1, PTangle1, SPLong, SPLat, SPangle, SPLong_NNR, SPLat_NNR, SPangle_NNR = \
            optimisation_methods.ApproximateNR_from_features(rotation_model_updated, nn_rotation_model, 
                                                             nr_timesteps, ref_rotation_plate_id)

            tmp_nr_eval = (np.sum(np.abs(PTangle1)) + np.mean(np.abs(PTangle1))) / 2

            nr_eval = tmp_nr_eval / net_rotation_weight



        #
        # Trench migration

        # Old method
        if data_array[2] == True and tm_method == 'convergence':

            kinArray = obj_func_convergence.kinloop(ref_rotation_end_age, ref_rotation_start_age, reformArray, 
                                                    rotation_model_tmp)

            cA = obj_func_convergence.kinstats(kinArray)
            cA = np.array(cA)

            trench_vel = -cA[:,6]
            trench_vel_SD = np.std(trench_vel)
            trench_numRetreating = len(np.where(trench_vel > 0)[0])
            trench_numAdvancing = len(trench_vel) - trench_numRetreating
            trench_numOver30 = len(np.where(trench_vel > 30)[0])
            trench_numLessNeg30 = len(np.where(trench_vel < -30)[0])
            trench_numTotal = len(trench_vel)
            trench_sumAbsVel_n = np.sum(np.abs(trench_vel)) / len(trench_vel)

            trench_percent_retreat = round((np.float(trench_numRetreating) / np.float(trench_numTotal)) * 100, 2)
            trench_percent_advance = 100. - trench_percent_retreat

            # Calculate cost
            #tm_eval_1 = (trench_percent_advance * 10) / trench_migration_weight
            #tm_eval_2 = (trench_sumAbsVel_n * 15) / trench_migration_weight

            # 1. trench percent advance + trench abs vel mean
            #tm_eval = (tm_eval_1 + tm_eval_2) / 2

            # 2. trench_abs_vel_mean
            #tm_eval_2 = (np.sum(np.abs(trench_vel)) / len(trench_vel)) / trench_migration_weight

            # 3. number of trenches in advance
            #tm_eval_3 = (trench_numAdvancing * 2) / trench_migration_weight

            # 4. abs median
            #tm_eval_4 = np.median(abs(trench_vel)) / trench_migration_weight

            # 5. standard deviation
            #tm_eval_5 = np.std(trench_vel) / trench_migration_weight

            # 6. variance
            #tm_stats = stats.describe(trench_vel)
            #tm_eval = tm_stats.variance / trench_migration_weight

            # 7. trench absolute motion abs vel mean
            #tm_eval_7 = ((np.sum(np.abs(trench_vel)) / len(trench_vel)) * 15) / trench_migration_weight

            #tm_eval = tm_eval_5



            #---- old ones
            # Minimise trench advance
            # tm_eval_1 = ((trench_percent_advance * 10) / trench_migration_weight)**2
            #tm_eval_1 = (trench_percent_advance * 10) / trench_migration_weight

            # Minimise trench velocities
            # tm_eval_2 = ((trench_sumAbsVel_n * 15) / trench_migration_weight)**2
            #tm_eval_2 = (trench_sumAbsVel_n * 15) / trench_migration_weight

            # Minimise trenches moving very fast (< or > 30)
            #tm_eval_3 = (trench_numOver30 + trench_numLessNeg30) * trench_migration_weight

            # # V1 (Original)
            # tmp_tm_eval = ((trench_vel_SD * (trench_numRetreating * trench_sumAbsVel_n)) / \
            #                (trench_numTotal - (trench_numOver30 + trench_numLessNeg30)))

            # tm_eval = tmp_tm_eval * trench_migration_weight



        # New method
        elif data_array[2] == True and tm_method == 'pygplates':

            tm_data = pgp.FeatureCollection(nnr_datadir + 'TMData_%sMa.gpml' % (int(ref_rotation_start_age)))

            # Calculate trench segment stats
            interval = 10.
            tm_stats = scap.subduction_absolute_motion(rotation_model_updated,
                                                       tm_data,
                                                       np.radians(1.),
                                                       ref_rotation_start_age - interval)

            # Process tm_stats to extract values for use in cost function
            trench_vel = []
            trench_obl = []

            for i in xrange(0, len(tm_stats)):

                trench_vel.append(tm_stats[i][2])
                trench_obl.append(tm_stats[i][3])

            trench_vel = np.array(trench_vel)
            trench_obl = np.array(trench_obl)

            # Scale velocities from cm to mm
            trench_vel = trench_vel * 10

            # Calculate trench orthogonal velocity
            tm_vel_orth = np.abs(trench_vel) * -np.cos(np.radians(trench_obl)) 


            trench_numTotal = len(tm_vel_orth)
            trench_numRetreating = len(np.where(tm_vel_orth > 0)[0])
            trench_numAdvancing = len(tm_vel_orth) - trench_numRetreating
            trench_percent_retreat = round((np.float(trench_numRetreating) / np.float(trench_numTotal)) * 100, 2)
            trench_percent_advance = 100. - trench_percent_retreat
            trench_sumAbsVel_n = np.sum(np.abs(tm_vel_orth)) / len(tm_vel_orth)
            trench_numOver30 = len(np.where(tm_vel_orth > 30)[0])
            trench_numLessNeg30 = len(np.where(tm_vel_orth < -30)[0])

            # Calculate cost
            #tm_eval_1 = (trench_percent_advance * 10) / trench_migration_weight
            #tm_eval_2 = (trench_sumAbsVel_n * 15) / trench_migration_weight

            # 1. trench percent advance + trench abs vel mean
            #tm_eval = (tm_eval_1 + tm_eval_2) / 2

            # 2. trench_abs_vel_mean orthogonal
            tm_eval_2 = (np.sum(np.abs(tm_vel_orth)) / len(tm_vel_orth)) / trench_migration_weight

            # 3. number of trenches in advance
            #tm_eval_3 = (trench_numAdvancing * 2) / trench_migration_weight

            # 4. abs median
            #tm_eval = np.median(abs(np.array(tm_vel_orth))) / trench_migration_weight

            # 5. standard deviation
            tm_eval_5 = (np.std(tm_vel_orth) / trench_migration_weight)

            # 6. variance
            #tm_stats = stats.describe(tm_vel_orth)
            #tm_eval_6 = tm_stats.variance / trench_migration_weight

            # 7. trench absolute motion abs vel mean
            #tm_eval_7 = ((np.sum(np.abs(trench_vel)) / len(trench_vel)) * 15) / trench_migration_weight

            tm_eval = ((tm_eval_2 + tm_eval_5) * 3) / trench_migration_weight
            
            # Original equation
            #tm_eval = ((tm_eval_5 * (trench_numRetreating * trench_sumAbsVel_n)) / \
            #           (trench_numTotal - (trench_numOver30 + trench_numLessNeg30)))
            
            #tm_eval = ((tm_eval_2 + tm_eval_5 + trench_numAdvancing) / trench_numRetreating) / trench_migration_weight



        # Hotspot trail distance misfit
        if data_array[3] == True:

            # returns: [point_distance_misfit, trail_distance_misfit, uncertainty, trail_name]
            hs = ObjectiveFunctions.hotspot_trail_misfit(trail_data, ref_rotation_start_age, 
                                                         rotation_model_updated, use_trail_age_uncertainty,
                                                         trail_age_uncertainty_ellipse)

            if use_trail_age_uncertainty == False:

                tmp_distance_median = np.median(hs[0])
                tmp_distance_sd = np.std(hs[0])

                hs_dist_eval = (tmp_distance_median + tmp_distance_sd) / hotspot_trails_weight


            elif use_trail_age_uncertainty == True:

                weighted_dist = []

                # Positively weight modelled distances that are less than uncertainty limit
                for i in xrange(0, len(hs[0])):

                    if hs[0][i] < hs[2][i]:

                        weighted_dist.append(hs[0][i] / 2)

                    else:

                        weighted_dist.append(hs[0][i] * 2)


                tmp_distance_median = np.median(weighted_dist)
                tmp_distance_sd = np.std(weighted_dist)

                hs_dist_eval = (tmp_distance_median + tmp_distance_sd) / hotspot_trails_weight
                #hs_dist_eval = tmp_distance_median / hotspot_trails_weight
                #hs_dist_eval = tmp_distance_sd / hotspot_trails_weight





        #### -----------------------------------------------------------------------------------------
        #### 3. Calculate evaluation return number

        # Scaling values
        alpha = 10
        beta = 100
        gamma = 1000


        opt_eval = 0

        # Fracture zones
        try:
            if fz_eval:
                opt_eval = opt_eval + (fz_eval * alpha)
        except:
            pass


        # Net rotation
        try:
            if nr_eval:
                opt_eval = opt_eval + (nr_eval * gamma)
        except:
            pass


        # Trench migration
        try:
            if tm_eval:
                #opt_eval = opt_eval + (tm_eval / alpha)
                opt_eval = opt_eval + tm_eval
        except:
            pass


        # Hotspot reconstruction distance + spherical dispersion statistics
        try:
            if hs_dist_eval and data_array[3] == True:

                # Distance only
                #opt_eval = opt_eval + hs_dist_eval

                # Kappa only
                #opt_eval = opt_eval + (hs_kappa_eval * 1e6)

                # Distance + Kappa
                #opt_eval = opt_eval + (((hs_kappa_eval * 1e6) + hs_dist_eval) / 1.5)

                # Distance misfit
                opt_eval = opt_eval + (hs_dist_eval / 2)

        except:
            pass




        #### ---------------------------------------------------------------------------------------------
        #### Return all calculated quantities     
        try:
            opt_eval_data.append(opt_eval)
        except:
            pass


        return opt_eval


    # --------------------------------------------------------------------
    # --------------------------------------------------------------------
    # Function to run optimisation routine

    def run_optimisation(x, opt_n, N, obj_f, lb, ub, model_stop_condition, max_iter, rotation_file, 
                         no_net_rotation_file, ref_rotation_start_age, Lats, Lons, spreading_directions, 
                         spreading_asymmetries, seafloor_ages, PID, CPID, data_array, nnr_datadir, 
                         ref_rotation_end_age, ref_rotation_plate_id, reformArray, trail_data,
                         fracture_zone_weight, net_rotation_weight, trench_migration_weight, hotspot_trails_weight,
                         use_trail_age_uncertainty, trail_age_uncertainty_ellipse, tm_method):

        import nlopt

        print tm_method

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


    # Map variables for use locally and by number of cores selected
    run_optimisation = dview['run_optimisation'] = run_optimisation
    opt_n = dview['opt_n'] = startingConditions[1]
    N = dview['N'] = startingConditions[2]
    obj_f = dview['obj_f'] = obj_f
    lb = dview['lb'] = startingConditions[3]
    ub = dview['ub'] = startingConditions[4]
    model_stop_condition = dview['model_stop_condition'] = startingConditions[5]
    max_iter = dview['max_iter'] = startingConditions[6]
    rotation_file = dview['rotation_file'] = data[1]
    ref_rotation_start_age = dview['ref_rotation_start_age'] = startingConditions[8]
    ref_rotation_end_age = dview['ref_rotation_end_age'] = startingConditions[9]
    ref_rotation_plate_id = dview['ref_rotation_plate_id'] = startingConditions[10]
    Lats = dview['Lats'] = startingConditions[11]
    Lons = dview['Lons'] = startingConditions[12]
    spreading_directions = dview['spreading_directions'] = startingConditions[13]
    spreading_asymmetries = dview['spreading_asymmetries'] = startingConditions[14]
    seafloor_ages = dview['seafloor_ages'] = startingConditions[15]
    PID = dview['PID'] = startingConditions[16]
    CPID = dview['CPID'] = startingConditions[17]
    data_array = dview['data_array'] = startingConditions[18]
    nnr_datadir = dview['nnr_datadir'] = startingConditions[19]
    no_net_rotation_file = dview['no_net_rotation_file'] = startingConditions[20]
    reformArray = dview['reformArray'] = startingConditions[21]
    trail_data = dview['trail_data'] = startingConditions[22]

    dview['fracture_zone_weight'] = fracture_zone_weight
    dview['net_rotation_weight'] = net_rotation_weight
    dview['trench_migration_weight'] = trench_migration_weight
    dview['hotspot_trails_weight'] = hotspot_trails_weight
    dview['use_trail_age_uncertainty'] = use_trail_age_uncertainty
    dview['trail_age_uncertainty_ellipse'] = trail_age_uncertainty_ellipse
    dview['tm_method'] = tm_method

    start_seeds = startingConditions[23]
    rotation_age_of_interest_age = startingConditions[24]
    data_array_labels_short = startingConditions[25]

    if auto_calc_ref_pole == True:

        ref_rot_longitude = startingConditions[26]
        ref_rot_latitude = startingConditions[27]
        ref_rot_angle = startingConditions[28]

    elif auto_calc_ref_pole == False:

        ref_rot_longitude = ref_rot_longitude
        ref_rot_latitude = ref_rot_latitude
        ref_rot_angle = ref_rot_angle

    seed_lons = startingConditions[29]
    seed_lats = startingConditions[30]

    x = startingConditions[0]

    minf = []

    #print "Number of start seeds generated:", len(start_seeds)
    print "Optimised models to be run:", len(start_seeds)
    print " "


    # --------------------------------------------------------------------
    # --------------------------------------------------------------------
    # Start optimisation

    # Start timer
    #start = time.time()
    main_start = time.time()

    # Run optimisation algorithm in parallel
    #try:

    prunopt = partial(run_optimisation, opt_n=opt_n, N=N, obj_f=obj_f, lb=lb, ub=ub, 
                      model_stop_condition=model_stop_condition, max_iter=max_iter, rotation_file=rotation_file,
                      no_net_rotation_file=no_net_rotation_file, ref_rotation_start_age=ref_rotation_start_age, 
                      Lats=Lats, Lons=Lons, spreading_directions=spreading_directions, 
                      spreading_asymmetries=spreading_asymmetries, 
                      seafloor_ages=seafloor_ages, PID=PID, CPID=CPID, data_array=data_array, nnr_datadir=nnr_datadir,
                      ref_rotation_end_age=ref_rotation_end_age, ref_rotation_plate_id=ref_rotation_plate_id,
                      reformArray=reformArray, trail_data=trail_data, fracture_zone_weight=fracture_zone_weight,
                      net_rotation_weight=net_rotation_weight, trench_migration_weight=trench_migration_weight,
                      hotspot_trails_weight=hotspot_trails_weight, use_trail_age_uncertainty=use_trail_age_uncertainty,
                      trail_age_uncertainty_ellipse=trail_age_uncertainty_ellipse, tm_method=tm_method)

    xopt = dview.map(prunopt, x)


    # except Exception as e:
    
    #     text_file = open("Output.txt", "w")
    #     text_file.write("Model error: " + str(e))
    #     text_file.close()


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
                                 plot=False)


    for j in xrange(0, len(xopt)):

        costs.append(xopt[j][1])


    rmin = np.array(rmin)
    rmean = np.array(rmean)

    min_results.append(rmin[0])
    mean_results.append(rmean[0])

    #import geoTools
    plat, plon = geoTools.checkLatLon(min_results[-1][2], min_results[-1][3])


    # end_time = round(time.time() - start, 2)
    # sec = timedelta(seconds = float(end_time))
    # dt = datetime(1,1,1) + sec
    
    # print "Timestep completed in:"
    # print str(dt.day-1) + "d, " + str(dt.hour) + "h, " + str(dt.minute) + "m, " + str(dt.second) + "s."


    # --------------------------------------------------------------------
    # --------------------------------------------------------------------
    # Update rotation file with result

    #rotation_file = datadir + 'Global_EarthByte_230-0Ma_GK07_AREPS_optAPM003.rot'
    rotation_file = datadir + rotfile

    rotation_model_tmp = pgp.FeatureCollection(rotation_file)

    # Find existing rotations for Africa
    opt_rotation_feature = None
    for rotation_feature in rotation_model_tmp:

        total_reconstruction_pole = rotation_feature.get_total_reconstruction_pole()

        if total_reconstruction_pole:

            fixed_plate_id, moving_plate_id, rotation_sequence = total_reconstruction_pole

            if fixed_plate_id == 001 and moving_plate_id == 701:
                opt_rotation_feature = rotation_feature
                break


    # Update existing rotation in the model with result
    if opt_rotation_feature:

        adjustment_time = pgp.GeoTimeInstant(ref_rotation_start_age)

        for finite_rotation_samples in rotation_sequence.get_enabled_time_samples():

            finite_rotation_time = finite_rotation_samples.get_time()

            if finite_rotation_time == ref_rotation_start_age:

                finite_rotation = finite_rotation_samples.get_value().get_finite_rotation()

                # new_rotation = pgp.FiniteRotation((np.double(round(min_results[-1][2], 2)), 
                #                                    np.double(round(min_results[-1][3], 2))), 
                #                                    np.radians(np.double(round(min_results[-1][1], 2))))

                new_rotation = pgp.FiniteRotation((np.double(round(plat, 2)), 
                                                   np.double(round(plon, 2))), 
                                                   np.radians(np.double(round(min_results[-1][1], 2))))

                finite_rotation_samples.get_value().set_finite_rotation(new_rotation)


    # Add result rotation pole to rotation file
    rotation_model_tmp.write(rotation_file)


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

print "Mean of 20 models (0-50Ma)"
print ""
print "tm_eval =", 47 * 3
print "nr_eval =", 143


import pickle

with open('model_output/optAPM175_10-0Ma_10models_NR_TM_60.pkl', 'rb') as f:
    data = pickle.load(f)
    
print data

