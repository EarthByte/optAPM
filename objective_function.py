import numpy as np
import pygplates as pgp
import optimisation_methods
import obj_func_convergence
import geoTools
import subduction_convergence_for_absolute_plate_motion as scap

from optapm import ObjectiveFunctions
from scipy import stats

# --------------------------------------------------------------------
# --------------------------------------------------------------------
# Objective function


class ObjectiveFunction(object):
    """
    Class containing objective function for optimisation.
    
    obj_f = ObjectiveFunction(...)
    result = obj_f(x, grad)
    """
    
    def __init__(
            self,
            interval,
            rotation_file,
            no_net_rotation_file,
            ref_rotation_start_age,
            Lats,
            Lons,
            spreading_directions,
            spreading_asymmetries,
            seafloor_ages,
            PID,
            CPID,
            data_array,
            nnr_datadir,
            ref_rotation_end_age,
            ref_rotation_plate_id,
            reformArray,
            trail_data,
            fracture_zone_weight,
            net_rotation_weight,
            trench_migration_weight,
            hotspot_trails_weight,
            use_trail_age_uncertainty,
            trail_age_uncertainty_ellipse,
            tm_method):
        
        #
        # Load/parse the feature collection files up front so we don't have to repeatedly do it in each objective function call.
        #
        # Prepare rotation model for updates during optimisation - keeps rotations in memory
        self.rotation_model_tmp = pgp.FeatureCollection(rotation_file)
        # Net rotation
        if data_array[1]:
            # Prepare no net rotation model for updates during optimisation - keeps rotations in memory
            self.nn_rotation_model = pgp.RotationModel(no_net_rotation_file)
        # Trench migration using pyGPlates.
        if data_array[2] and tm_method == 'pygplates':
            self.tm_data = pgp.FeatureCollection(nnr_datadir + 'TMData_%sMa.gpml' % (int(ref_rotation_start_age)))
        
        self.interval = interval
        self.rotation_file = rotation_file
        self.no_net_rotation_file = no_net_rotation_file
        self.ref_rotation_start_age = ref_rotation_start_age
        self.Lats = Lats
        self.Lons = Lons
        self.spreading_directions = spreading_directions
        self.spreading_asymmetries = spreading_asymmetries
        self.seafloor_ages = seafloor_ages
        self.PID = PID
        self.CPID = CPID
        self.data_array = data_array
        self.nnr_datadir = nnr_datadir
        self.ref_rotation_end_age = ref_rotation_end_age
        self.ref_rotation_plate_id = ref_rotation_plate_id
        self.reformArray = reformArray
        self.trail_data = trail_data
        self.fracture_zone_weight = fracture_zone_weight
        self.net_rotation_weight = net_rotation_weight
        self.trench_migration_weight = trench_migration_weight
        self.hotspot_trails_weight = hotspot_trails_weight
        self.use_trail_age_uncertainty = use_trail_age_uncertainty
        self.trail_age_uncertainty_ellipse = trail_age_uncertainty_ellipse
        self.tm_method = tm_method


    def __call__(self, x, grad):

        #### -----------------------------------------------------------------------------------------
        #### 1. Calculate reconstructed data point locations

        tmp_opt_rlon = []
        tmp_opt_rlat = []
        opt_stats = []

        # Check incoming Africa finite rotation pole values
        lat_, lon_ = geoTools.checkLatLon(x[1], x[0])
        ang_ = x[2]


        #### -----------------------------------------------------------------------------------------
        #### 2. Find and update Africa rotation


        # Find existing rotation for Africa at correct time
        opt_rotation_feature = None
        for rotation_feature in self.rotation_model_tmp:

            total_reconstruction_pole = rotation_feature.get_total_reconstruction_pole()

            if total_reconstruction_pole:

                fixed_plate_id, moving_plate_id, rotation_sequence = total_reconstruction_pole

                if fixed_plate_id == 001 and moving_plate_id == 701:

                    opt_rotation_feature = rotation_feature
                    break


        # Update rotation file with proposed Africa rotation
        if opt_rotation_feature:

            adjustment_time = pgp.GeoTimeInstant(self.ref_rotation_start_age)

            for finite_rotation_samples in rotation_sequence.get_enabled_time_samples():

                finite_rotation_time = finite_rotation_samples.get_time()

                if finite_rotation_time == self.ref_rotation_start_age:

                    finite_rotation = finite_rotation_samples.get_value().get_finite_rotation()

                    new_rotation = pgp.FiniteRotation((np.double(lat_), np.double(lon_)), 
                                                      np.radians(np.double(ang_)))
                    finite_rotation_samples.get_value().set_finite_rotation(new_rotation)

        rotation_model_updated = pgp.RotationModel(
                self.rotation_model_tmp,
                # OPTIMIZATION: We need to be careful setting this to False - we should ensure
                # that we'll never modify the rotation features 'self.rotation_model_tmp' while
                # 'rotation_model_updated' is being used (ie, calling one of its methods).
                clone_rotation_features=False)




        #### -----------------------------------------------------------------------------------------
        #### 3. Calculate data fits


        #
        # Fracture zone orientation
        if self.data_array[0] == True:

            # Get skew values
            fz = optimisation_methods.Calc_Median(rotation_model_updated, self.PID, 
                                                  self.seafloor_ages, self.Lats, self.Lons, 
                                                  self.spreading_directions)


            tmp_fz_eval = fz[0] + fz[1]

            fz_eval = tmp_fz_eval / self.fracture_zone_weight



        #
        # Net rotation
        if self.data_array[1] == True:

            nr_timesteps = np.arange(self.ref_rotation_end_age, self.ref_rotation_start_age + 1, 2)

            PTLong1, PTLat1, PTangle1, SPLong, SPLat, SPangle, SPLong_NNR, SPLat_NNR, SPangle_NNR = \
            optimisation_methods.ApproximateNR_from_features(rotation_model_updated, self.nn_rotation_model, 
                                                             nr_timesteps, self.ref_rotation_plate_id)

            tmp_nr_eval = (np.sum(np.abs(PTangle1)) + np.mean(np.abs(PTangle1))) / 2

            nr_eval = tmp_nr_eval / self.net_rotation_weight



        #
        # Trench migration

        # Old method
        if self.data_array[2] == True and self.tm_method == 'convergence':

            kinArray = obj_func_convergence.kinloop(
                    self.ref_rotation_end_age,
                    self.ref_rotation_start_age,
                    self.reformArray, 
                    self.rotation_model_tmp)

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
            #tm_eval_1 = (trench_percent_advance * 10) / self.trench_migration_weight
            #tm_eval_2 = (trench_sumAbsVel_n * 15) / self.trench_migration_weight

            # 1. trench percent advance + trench abs vel mean
            #tm_eval = (tm_eval_1 + tm_eval_2) / 2

            # 2. trench_abs_vel_mean
            #tm_eval_2 = (np.sum(np.abs(trench_vel)) / len(trench_vel)) / self.trench_migration_weight

            # 3. number of trenches in advance
            #tm_eval_3 = (trench_numAdvancing * 2) / self.trench_migration_weight

            # 4. abs median
            #tm_eval_4 = np.median(abs(trench_vel)) / self.trench_migration_weight

            # 5. standard deviation
            #tm_eval_5 = np.std(trench_vel) / self.trench_migration_weight

            # 6. variance
            #tm_stats = stats.describe(trench_vel)
            #tm_eval = tm_stats.variance / self.trench_migration_weight

            # 7. trench absolute motion abs vel mean
            #tm_eval_7 = ((np.sum(np.abs(trench_vel)) / len(trench_vel)) * 15) / self.trench_migration_weight

            #tm_eval = tm_eval_5



            #---- old ones
            # Minimise trench advance
            # tm_eval_1 = ((trench_percent_advance * 10) / self.trench_migration_weight)**2
            #tm_eval_1 = (trench_percent_advance * 10) / self.trench_migration_weight

            # Minimise trench velocities
            # tm_eval_2 = ((trench_sumAbsVel_n * 15) / self.trench_migration_weight)**2
            #tm_eval_2 = (trench_sumAbsVel_n * 15) / self.trench_migration_weight

            # Minimise trenches moving very fast (< or > 30)
            #tm_eval_3 = (trench_numOver30 + trench_numLessNeg30) * self.trench_migration_weight

            # # V1 (Original)
            # tmp_tm_eval = ((trench_vel_SD * (trench_numRetreating * trench_sumAbsVel_n)) / \
            #                (trench_numTotal - (trench_numOver30 + trench_numLessNeg30)))

            # tm_eval = tmp_tm_eval * self.trench_migration_weight



        # New method
        elif self.data_array[2] == True and self.tm_method == 'pygplates':

            # Calculate trench segment stats
            tm_stats = scap.subduction_absolute_motion(rotation_model_updated,
                                                       self.tm_data,
                                                       np.radians(1.),
                                                       self.ref_rotation_start_age - self.interval)

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
            #tm_eval_1 = (trench_percent_advance * 10) / self.trench_migration_weight
            #tm_eval_2 = (trench_sumAbsVel_n * 15) / self.trench_migration_weight

            # 1. trench percent advance + trench abs vel mean
            #tm_eval = (tm_eval_1 + tm_eval_2) / 2

            # 2. trench_abs_vel_mean orthogonal
            tm_eval_2 = (np.sum(np.abs(tm_vel_orth)) / len(tm_vel_orth)) / self.trench_migration_weight

            # 3. number of trenches in advance
            #tm_eval_3 = (trench_numAdvancing * 2) / self.trench_migration_weight

            # 4. abs median
            #tm_eval = np.median(abs(np.array(tm_vel_orth))) / self.trench_migration_weight

            # 5. standard deviation
            tm_eval_5 = (np.std(tm_vel_orth) / self.trench_migration_weight)

            # 6. variance
            #tm_stats = stats.describe(tm_vel_orth)
            #tm_eval_6 = tm_stats.variance / self.trench_migration_weight

            # 7. trench absolute motion abs vel mean
            #tm_eval_7 = ((np.sum(np.abs(trench_vel)) / len(trench_vel)) * 15) / self.trench_migration_weight

            tm_eval = ((tm_eval_2 + tm_eval_5) * 3) / self.trench_migration_weight
            
            # Original equation
            #tm_eval = ((tm_eval_5 * (trench_numRetreating * trench_sumAbsVel_n)) / \
            #           (trench_numTotal - (trench_numOver30 + trench_numLessNeg30)))
            
            #tm_eval = ((tm_eval_2 + tm_eval_5 + trench_numAdvancing) / trench_numRetreating) / self.trench_migration_weight



        # Hotspot trail distance misfit
        if self.data_array[3] == True:

            # returns: [point_distance_misfit, trail_distance_misfit, uncertainty, trail_name]
            hs = ObjectiveFunctions.hotspot_trail_misfit(self.trail_data, self.ref_rotation_start_age, 
                                                         rotation_model_updated, self.use_trail_age_uncertainty,
                                                         self.trail_age_uncertainty_ellipse)

            if self.use_trail_age_uncertainty == False:

                tmp_distance_median = np.median(hs[0])
                tmp_distance_sd = np.std(hs[0])

                hs_dist_eval = (tmp_distance_median + tmp_distance_sd) / self.hotspot_trails_weight


            elif self.use_trail_age_uncertainty == True:

                weighted_dist = []

                # Positively weight modelled distances that are less than uncertainty limit
                for i in xrange(0, len(hs[0])):

                    if hs[0][i] < hs[2][i]:

                        weighted_dist.append(hs[0][i] / 2)

                    else:

                        weighted_dist.append(hs[0][i] * 2)


                tmp_distance_median = np.median(weighted_dist)
                tmp_distance_sd = np.std(weighted_dist)

                hs_dist_eval = (tmp_distance_median + tmp_distance_sd) / self.hotspot_trails_weight
                #hs_dist_eval = tmp_distance_median / self.hotspot_trails_weight
                #hs_dist_eval = tmp_distance_sd / self.hotspot_trails_weight





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
            if hs_dist_eval and self.data_array[3] == True:

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
