import sys
import numpy as np
import time
import pygplates as pgp
import geoTools
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
    sys.exit(1)

dview = rc[:]
dview.block = True

# Make sure remote engines also import these modules.
with dview.sync_imports():
    from objective_function import ObjectiveFunction
    import nlopt


def optimise_APM():

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
    models = 10

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
        # Function to run optimisation routine

        def run_optimisation(x, opt_n, N, lb, ub, model_stop_condition, max_iter, interval, rotation_file, 
                             no_net_rotation_file, ref_rotation_start_age, Lats, Lons, spreading_directions, 
                             spreading_asymmetries, seafloor_ages, PID, CPID, data_array, nnr_datadir, 
                             ref_rotation_end_age, ref_rotation_plate_id, reformArray, trail_data,
                             fracture_zone_weight, net_rotation_weight, trench_migration_weight, hotspot_trails_weight,
                             use_trail_age_uncertainty, trail_age_uncertainty_ellipse, tm_method):

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


        # Map variables for use locally and by number of cores selected
        run_optimisation = dview['run_optimisation'] = run_optimisation
        opt_n = dview['opt_n'] = startingConditions[1]
        N = dview['N'] = startingConditions[2]
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

        prunopt = partial(run_optimisation, opt_n=opt_n, N=N, lb=lb, ub=ub, 
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


if __name__ == '__main__':

    optimise_APM()