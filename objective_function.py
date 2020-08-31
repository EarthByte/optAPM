import numpy as np
import pygplates as pgp
import optimisation_methods
import geoTools
import subduction_convergence_for_absolute_plate_motion as scap
import sys

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
            weights_array,
            cost_func_array,
            bounds_array,
            trench_migration_file,
            plate_velocity_file,
            ref_rotation_end_age,
            ref_rotation_plate_id,
            reformArray,
            trail_data,
            use_trail_age_uncertainty,
            trail_age_uncertainty_ellipse,
            tm_method):
        
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
        self.weights_array = weights_array
        self.bounds_array = bounds_array
        self.cost_func_array = cost_func_array
        self.trench_migration_file = trench_migration_file
        self.plate_velocity_file = plate_velocity_file
        self.ref_rotation_end_age = ref_rotation_end_age
        self.ref_rotation_plate_id = ref_rotation_plate_id
        self.reformArray = reformArray
        self.trail_data = trail_data
        self.use_trail_age_uncertainty = use_trail_age_uncertainty
        self.trail_age_uncertainty_ellipse = trail_age_uncertainty_ellipse
        self.tm_method = tm_method
        
        
        #
        # Load/parse the feature collection files up front so we don't have to repeatedly do it in each objective function call.
        #
        # Prepare rotation model for updates during optimisation - keeps rotations in memory
        rotation_features = pgp.FeatureCollection(rotation_file)
        self.rotation_features_updated = rotation_features
        # Also keep original rotation model to help when inserting rotation updates.
        self.rotation_model_original = pgp.RotationModel(rotation_features)
        # Net rotation
        if data_array[1]:
            # Prepare no net rotation model for updates during optimisation - keeps rotations in memory
            self.nn_rotation_model = pgp.RotationModel(no_net_rotation_file)
        # Trench migration using pyGPlates.
        if data_array[2] and tm_method == 'pygplates':
            self.tm_data = pgp.FeatureCollection(trench_migration_file)
        # Plate velocities.
        # Contains multi-points (and each multi-point has a plate ID).
        # Each multi-point represents those grid points that fall within a resolved plate
        # (or continental polygon) with a particular plate ID. Velocities can then be calculated using
        # each candidate rotation model and the plate IDs (and positions).
        if data_array[4]:
            self.pv_data = pgp.FeatureCollection(plate_velocity_file)


        #
        # Find the 005-000 rotation at correct time.
        #
        # Store the correct time sample for the 005-000 rotation as 'self.opt_finite_rotation_sample'
        # so we don't have to do it every time this objective function is called.

        for rotation_feature in rotation_features:
            total_reconstruction_pole = rotation_feature.get_total_reconstruction_pole()
            if total_reconstruction_pole:
                fixed_plate_id, moving_plate_id, rotation_sequence = total_reconstruction_pole
                if moving_plate_id == 5 and fixed_plate_id == 0:
                    for finite_rotation_sample in rotation_sequence.get_enabled_time_samples():
                        if finite_rotation_sample.get_time() == ref_rotation_start_age:
                            self.opt_finite_rotation_sample = finite_rotation_sample
                            break
                    break
        
        self.debug_count = 0
        
        # To debug the weighted cost functions (net rotation, trench migration, etc).
        self.debug_data_array = []


    def __call__(self, x, grad):

        # print self.debug_count
        # sys.stdout.flush()
        self.debug_count += 1
        
        #### -----------------------------------------------------------------------------------------
        #### 1. Calculate reconstructed data point locations

        tmp_opt_rlon = []
        tmp_opt_rlat = []
        opt_stats = []

        # Check incoming reference plate finite rotation pole values
        lat_, lon_ = geoTools.checkLatLon(x[1], x[0])
        ang_ = x[2]


        #### -----------------------------------------------------------------------------------------
        #### 2. Update reference plate rotation


        # Our new rotation is the 'ref_rotation_plate_id' relative to the optimisation root plate (000).
        new_rotation_ref_plate_rel_000 = pgp.FiniteRotation((np.double(lat_), np.double(lon_)), np.radians(np.double(ang_)))
        
        # Remove the 'ref_rotation_plate_id' relative to 005 part (of our new rotation) to get the
        # 005 relative to 000 (optimisation root plate) part that gets stored in the 005-000 rotation feature.
        #
        #                                     R(0->t,000->ref_plate) = R(0->t,000->005) * R(0->t,005->ref_plate)
        #   R(0->t,000->ref_plate) * inverse(R(0->t,005->ref_plate)) = R(0->t,000->005)
        #
        plate_rotation_ref_plate_rel_005 = self.rotation_model_original.get_rotation(
                self.ref_rotation_start_age,
                self.ref_rotation_plate_id,
                fixed_plate_id=5)
        new_rotation_005_rel_000 = new_rotation_ref_plate_rel_000 * plate_rotation_ref_plate_rel_005.get_inverse()
        
        # Update the 005-000 rotation.
        # Note that this modifies the state of 'self.rotation_features_updated' - in other words,
        # we're modifying a time sample of one of the rotation features in that list of features.
        self.opt_finite_rotation_sample.get_value().set_finite_rotation(new_rotation_005_rel_000)
        
        rotation_model_updated = pgp.RotationModel(
                self.rotation_features_updated,
                # OPTIMIZATION: We need to be careful setting this to False - we should ensure
                # that we'll never modify the rotation features 'self.rotation_features_updated' while
                # 'rotation_model_updated' is being used (ie, calling one of its methods).
                #
                # UPDATE: This optimisation is ignored for pygplates revision >= 25 since the
                #         'clone_rotation_features' argument has been deprecated (ie, subsequent modifications
                #         to rotation features no longer affect a RotationModel after it's created).
                #         Also the internals of RotationModel have been significantly optimised in revision 25.
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
            
            # Delegate cost evaluation to cost function.
            fz_eval = self.cost_func_array[0](fz)
            
            # Penalise out-of-bound cost values (if requested).
            if self.bounds_array[0]:
                fz_lower_bound, fz_upper_bound = self.bounds_array[0]
                # TODO: Might need to compare 'fz[1]' (ie, the  mean) instead of 'fz_eval'.
                if fz_eval < fz_lower_bound or fz_eval > fz_upper_bound:
                    # Arbitrary penalty on cost function (might need some tuning)
                    fz_eval += 10000.0

            fz_eval /= self.weights_array[0]



        #
        # Net rotation
        if self.data_array[1] == True:

            # List of times, not including 'self.ref_rotation_start_age' (ie, last value is
            # 'self.ref_rotation_start_age - 1') since 'optimisation_methods.ApproximateNR_from_features()'
            # calculates NR from 't+1' to 't' (and so 't+1' will be 'self.ref_rotation_start_age').
            #
            # This is important because sampling beyond 'self.ref_rotation_start_age' will create problems
            # since 'rotation_model_updated' has only been updated back to 'self.ref_rotation_start_age'.
            nr_timesteps = range(self.ref_rotation_end_age, self.ref_rotation_start_age, 1)

            PTLong1, PTLat1, PTangle1, SPLong, SPLat, SPangle, SPLong_NNR, SPLat_NNR, SPangle_NNR = \
            optimisation_methods.ApproximateNR_from_features(rotation_model_updated, self.nn_rotation_model, 
                                                             nr_timesteps, self.ref_rotation_plate_id)

            # Sum the net rotation absolute angles over 'nr_timesteps' to get full 'interval'.
            nr_over_interval = np.sum(np.abs(PTangle1))
            
            # Delegate cost evaluation to cost function.
            nr_eval = self.cost_func_array[1](PTLong1, PTLat1, PTangle1, SPLong, SPLat, SPangle, SPLong_NNR, SPLat_NNR, SPangle_NNR, nr_over_interval)
            
            # Penalise out-of-bound cost values (if requested).
            if self.bounds_array[1]:
                nr_lower_bound, nr_upper_bound = self.bounds_array[1]
                # 'nr_over_interval' is over whole 'interval', but the bounded values are in deg/Myr,
                # so convert our NR to deg/Myr (which is the units of the lower/upper bounds).
                nr_deg_per_myr = nr_over_interval / self.interval
                if nr_deg_per_myr < nr_lower_bound or nr_deg_per_myr > nr_upper_bound:
                    # Arbitrary penalty on cost function (might need some tuning)
                    nr_eval += 10000.0

            nr_eval /= self.weights_array[1]



        #
        # Trench migration

        # Old method
        if self.data_array[2] == True and self.tm_method == 'convergence':
            
            # No longer using old path (should use new pygplates path instead).
            # Have removed 'import obj_func_convergence'.
            raise NotImplementedError(
                'Deprecated old convergence path, use new pygplates trench migration path instead.')
            
            kinArray = obj_func_convergence.kinloop(
                    self.ref_rotation_end_age,
                    self.ref_rotation_start_age,
                    self.reformArray, 
                    self.rotation_features_updated)

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
            #tm_eval_1 = (trench_percent_advance * 10) / self.weights_array[2]
            #tm_eval_2 = (trench_sumAbsVel_n * 15) / self.weights_array[2]

            # 1. trench percent advance + trench abs vel mean
            #tm_eval = (tm_eval_1 + tm_eval_2) / 2

            # 2. trench_abs_vel_mean
            #tm_eval_2 = (np.sum(np.abs(trench_vel)) / len(trench_vel)) / self.weights_array[2]

            # 3. number of trenches in advance
            #tm_eval_3 = (trench_numAdvancing * 2) / self.weights_array[2]

            # 4. abs median
            #tm_eval_4 = np.median(abs(trench_vel)) / self.weights_array[2]

            # 5. standard deviation
            #tm_eval_5 = np.std(trench_vel) / self.weights_array[2]

            # 6. variance
            #tm_stats = stats.describe(trench_vel)
            #tm_eval = tm_stats.variance / self.weights_array[2]

            # 7. trench absolute motion abs vel mean
            #tm_eval_7 = ((np.sum(np.abs(trench_vel)) / len(trench_vel)) * 15) / self.weights_array[2]

            #tm_eval = tm_eval_5



            #---- old ones
            # Minimise trench advance
            # tm_eval_1 = ((trench_percent_advance * 10) / self.weights_array[2])**2
            #tm_eval_1 = (trench_percent_advance * 10) / self.weights_array[2]

            # Minimise trench velocities
            # tm_eval_2 = ((trench_sumAbsVel_n * 15) / self.weights_array[2])**2
            #tm_eval_2 = (trench_sumAbsVel_n * 15) / self.weights_array[2]

            # Minimise trenches moving very fast (< or > 30)
            #tm_eval_3 = (trench_numOver30 + trench_numLessNeg30) * self.weights_array[2]

            # # V1 (Original)
            # tmp_tm_eval = ((trench_vel_SD * (trench_numRetreating * trench_sumAbsVel_n)) / \
            #                (trench_numTotal - (trench_numOver30 + trench_numLessNeg30)))

            # tm_eval = tmp_tm_eval * self.weights_array[2]

            raise NotImplementedError("Trench migration cost no longer implemented using old convergence script - use new script instead")



        # New method
        elif self.data_array[2] == True and self.tm_method == 'pygplates':

            # Calculate trench segment stats.
            #
            # Note that migration velocities are calculated from 'self.ref_rotation_start_age' to
            # 'self.ref_rotation_start_age - self.interval' (which is 'self.ref_rotation_end_age').
            tm_stats = scap.subduction_absolute_motion(rotation_model_updated,
                                                       self.tm_data,
                                                       np.radians(1.),
                                                       self.ref_rotation_start_age,
                                                       velocity_delta_time=self.interval)

            # Process tm_stats to extract values for use in cost function
            trench_vel = np.array([tm_stat[2] for tm_stat in tm_stats])
            trench_obl = np.array([tm_stat[3] for tm_stat in tm_stats])

            # Scale velocities from cm/yr to mm/yr.
            # Note that mm/yr is same as km/Myr.
            trench_vel = trench_vel * 10

            # Calculate trench orthogonal velocity
            tm_vel_orth = np.abs(trench_vel) * -np.cos(np.radians(trench_obl)) 

            # Mean of trench orthogonal velocity.
            tm_mean_vel_orth = np.sum(tm_vel_orth) / len(tm_vel_orth)

            # Mean of absolute trench orthogonal velocity.
            tm_mean_abs_vel_orth = np.sum(np.abs(tm_vel_orth)) / len(tm_vel_orth)
            
            # Delegate cost evaluation to cost function.
            tm_eval = self.cost_func_array[2](trench_vel, trench_obl, tm_vel_orth, tm_mean_vel_orth, tm_mean_abs_vel_orth)
            
            # Penalise out-of-bound cost values (if requested).
            if self.bounds_array[2]:
                tm_lower_bound, tm_upper_bound = self.bounds_array[2]
                # Note that we compare the mean of trench orthogonal velocity (not 'tm_eval').
                # Note that we're not comparing 'absolute' velocity, so negative/positive values are trench advance/retreat.
                # Also note that mean and bounds velocities are in same units of mm/yr (equivalent to km/Myr).
                if tm_mean_vel_orth < tm_lower_bound or tm_mean_vel_orth > tm_upper_bound:
                    # Arbitrary penalty on cost function (might need some tuning)
                    tm_eval += 10000.0
            
            tm_eval /= self.weights_array[2]


        # Hotspot trail distance misfit
        if self.data_array[3] == True:

            # returns: [point_distance_misfit, trail_distance_misfit, uncertainty, trail_name]
            hs = ObjectiveFunctions.hotspot_trail_misfit(self.trail_data, self.ref_rotation_start_age, 
                                                         rotation_model_updated, self.use_trail_age_uncertainty,
                                                         self.trail_age_uncertainty_ellipse)
            
            if self.use_trail_age_uncertainty == False:

                distance_median = np.median(hs[0])
                distance_sd = np.std(hs[0])

            else:
                weighted_dist = []

                # Positively weight modelled distances that are less than uncertainty limit
                for i in range(0, len(hs[0])):

                    if hs[0][i] < hs[2][i]:

                        weighted_dist.append(hs[0][i] / 2)

                    else:

                        weighted_dist.append(hs[0][i] * 2)


                distance_median = np.median(weighted_dist)
                distance_sd = np.std(weighted_dist)
            
            # Delegate cost evaluation to cost function.
            hs_dist_eval = self.cost_func_array[3](hs, distance_median, distance_sd)
            
            # Penalise out-of-bound cost values (if requested).
            if self.bounds_array[3]:
                hs_lower_bound, hs_upper_bound = self.bounds_array[3]
                # TODO: Might need to compare (possibly weighted) mean instead of 'hs_dist_eval'.
                if hs_dist_eval < hs_lower_bound or hs_dist_eval > hs_upper_bound:
                    # Arbitrary penalty on cost function (might need some tuning)
                    hs_dist_eval += 10000.0
            
            hs_dist_eval /= self.weights_array[3]
 


        #
        # Plate velocities
        #
        if self.data_array[4] == True:
            
            # Calculate velocity vectors at all pre-calculated grid points ('self.pv_data').
            velocity_vectors = []
            
            # 'self.pv_data' contains multi-points (and each multi-point has a plate ID).
            # Each multi-point represents those grid points that fall within a resolved plate
            # (or continental polygon) with a particular plate ID. Velocities can then be calculated
            # using the current updated rotation model and the plate IDs (and positions).
            for multi_point_feature in self.pv_data:
                plate_id = multi_point_feature.get_reconstruction_plate_id()
                
                # Get equivalent stage rotation from 'self.ref_rotation_start_age' to
                # 'self.ref_rotation_start_age - self.interval'
                equivalent_stage_rotation = rotation_model_updated.get_rotation(
                        self.ref_rotation_start_age - self.interval,
                        plate_id,
                        self.ref_rotation_start_age)
                
                # There is only one (multi-point) geometry per feature but assume there could be more.
                for multi_point in multi_point_feature.get_geometries():
                    multi_point_velocity_vectors = pgp.calculate_velocities(
                            multi_point,
                            equivalent_stage_rotation,
                            self.interval,
                            # Units of km/Myr (equivalent to mm/yr)...
                            velocity_units=pgp.VelocityUnits.kms_per_my)
                    velocity_vectors.extend(multi_point_velocity_vectors)

            # Velocity vectors at all pre-calculated grid points.
            velocity_magnitudes = [velocity_vector.get_magnitude() for velocity_vector in velocity_vectors]
            
            median_velocity = np.median(velocity_magnitudes)
            
            # Delegate cost evaluation to cost function.
            pv_eval = self.cost_func_array[4](velocity_vectors, velocity_magnitudes, median_velocity)
            
            # Penalise out-of-bound cost values (if requested).
            if self.bounds_array[4]:
                pv_lower_bound, pv_upper_bound = self.bounds_array[4]
                # Note that we compare the median velocity (not 'pv_eval').
                # Currently they're the same, but might not be in future.
                # Also note that median and bounds velocities are in same units of mm/yr (equivalent to km/Myr).
                if median_velocity < pv_lower_bound or median_velocity > pv_upper_bound:
                    # Arbitrary penalty on cost function (might need some tuning)
                    pv_eval += 10000.0
            
            pv_eval /= self.weights_array[4]





        #### -----------------------------------------------------------------------------------------
        #### 3. Calculate evaluation return number

        # Scaling values
        scale_fracture_zones = 10
        scale_net_rotation = 1000.0 / 8.0
        scale_trench_migration = 1.0
        scale_hot_spots = 1.0 / 8.0
        scale_plate_velocity = 10.0

        # To debug the weighted cost functions (net rotation, trench migration, etc).
        # Stores data for the current iteration (which then gets appended to 'self.debug_data_array').
        debug_data = []

        opt_eval = 0

        # Fracture zones
        try:
            if self.data_array[0] == True:
                opt_eval = opt_eval + (fz_eval * scale_fracture_zones)
                debug_data.append(fz_eval * scale_fracture_zones)
        except:
            pass


        # Net rotation
        try:
            if self.data_array[1] == True:
                opt_eval = opt_eval + (nr_eval * scale_net_rotation)
                debug_data.append(nr_eval * scale_net_rotation)
        except:
            pass


        # Trench migration
        try:
            if self.data_array[2] == True:
                opt_eval = opt_eval + (tm_eval * scale_trench_migration)
                debug_data.append(tm_eval * scale_trench_migration)
        except:
            pass


        # Hotspot reconstruction distance + spherical dispersion statistics
        try:
            if self.data_array[3] == True:

                # Distance only
                #opt_eval = opt_eval + hs_dist_eval

                # Kappa only
                #opt_eval = opt_eval + (hs_kappa_eval * 1e6)

                # Distance + Kappa
                #opt_eval = opt_eval + (((hs_kappa_eval * 1e6) + hs_dist_eval) / 1.5)

                # Distance misfit
                opt_eval = opt_eval + (hs_dist_eval * scale_hot_spots)
                debug_data.append(hs_dist_eval * scale_hot_spots)

        except:
            pass


        # Plate velocities
        try:
            if self.data_array[4] == True:
                opt_eval = opt_eval + (pv_eval * scale_plate_velocity)
                debug_data.append(pv_eval * scale_plate_velocity)
        except:
            pass


        #### ---------------------------------------------------------------------------------------------
        #### Return all calculated quantities     
        # try:
        #     opt_eval_data.append(opt_eval)
        # except:
        #     pass
        
        # To debug the weighted cost functions (net rotation, trench migration, etc).
        self.debug_data_array.append(debug_data)

        return opt_eval
