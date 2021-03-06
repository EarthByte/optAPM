import math
import net_rotation
import numpy as np
import os.path
import pygplates
import sys
import warnings


class NoNetRotationModel(object):
    """
    Class to calculate no-net-rotation model.
    
    The original input rotation files are combined into a single no-net-rotation file,
    returned by 'get_no_net_rotation_filename()' (relative to the 'data/' directory).
    """
    
    # Whether to create the entire no-net-rotation model in constructor ('__init__').
    # If False then rotation model is incrementally built during optimisation iteration over reconstruction time.
    CREATE_NO_NET_ROTATION_MODEL_AT_INIT = False
    
    
    def __init__(
            self,
            data_dir,
            original_rotation_filenames,  # Relative to the 'data/' directory
            topology_features,
            start_age,
            end_age,
            data_model,
            model_name,
            # Temporary: Allow input of GPlates exported net rotation file.
            # TODO: Remove when we can calculate net rotation in pygplates for a deforming model.
            #       This class currently only calculates net rotation for non-deforming models.
            gplates_net_rotation_filename=None):  # Relative to the 'data/' directory
        """
        Create a single no-net-rotation file by combining all original (input) rotations and
        removing net rotation to produce a no-net-rotation model.
        """
        self.topology_features = topology_features
        self.data_dir = data_dir
        
        # The single combined no-net-rotation filename (relative to the 'data/' directory).
        self.no_net_rotation_filename = os.path.join(
                data_model, 'optimisation', 'no_net_rotation_model_' + model_name + '.rot')
        
        # Temporary: Allow input of GPlates exported net rotation file.
        # TODO: Remove when we can calculate net rotation in pygplates for a deforming model.
        #       This class currently only calculates net rotation for non-deforming models.
        if gplates_net_rotation_filename:
            self.gplates_net_rotations = self._read_gplates_net_rotations(
                    os.path.join(self.data_dir, gplates_net_rotation_filename))
        else:
            self.gplates_net_rotations = None
        
        #
        # Combine the original (input) rotation files into a single no-net-rotation file.
        #
        
        # Load all the original rotation feature collections.
        rotation_features_except_005_rel_000 = []
        for rotation_filename in original_rotation_filenames:
            # Read the current rotation file.
            rotation_feature_collection = pygplates.FeatureCollection(
                    os.path.join(self.data_dir, rotation_filename))
            rotation_features_except_005_rel_000.extend(rotation_feature_collection)
        
        # Remove any existing 005-000 rotation features and change fixed plate 000 to 005.
        rotation_feature_index = 0
        while rotation_feature_index < len(rotation_features_except_005_rel_000):
            rotation_feature = rotation_features_except_005_rel_000[rotation_feature_index]
            total_reconstruction_pole = rotation_feature.get_total_reconstruction_pole()
            if total_reconstruction_pole:
                fixed_plate_id, moving_plate_id, rotation_sequence = total_reconstruction_pole
                
                # Change fixed plate ID 000 to 005 (unless it's the 005-000 plate pair).
                # We want everything that references 000 to now reference 005.
                # Later we'll add a 005-000 sequence to store no-net-rotation adjustments.
                if fixed_plate_id == 0:
                    if moving_plate_id == 5:
                        # Exclude 005-000 rotation features (we'll add a new one later).
                        del rotation_features_except_005_rel_000[rotation_feature_index]
                        rotation_feature_index -= 1
                    else:
                        # Change the fixed plate ID from 000 to 005.
                        rotation_feature.set_total_reconstruction_pole(5, moving_plate_id, rotation_sequence)
            
            rotation_feature_index += 1
        
        # Make 005-000 the identity rotation and use it (and the other features) to calculate the net rotation.
        rotation_feature_zero_005_rel_000 = pygplates.Feature.create_total_reconstruction_sequence(
            0,
            5,
            pygplates.GpmlIrregularSampling([
                pygplates.GpmlTimeSample(pygplates.GpmlFiniteRotation(pygplates.FiniteRotation()), 0.0, 'NNR'),
                pygplates.GpmlTimeSample(pygplates.GpmlFiniteRotation(pygplates.FiniteRotation()), start_age, 'NNR')]))
        
        # Rotation model for zero 005-000.
        # Once we've calculated net rotation we'll remove net rotation (which will make 005-000 non-zero).
        rotation_model_zero_005_rel_000 = pygplates.RotationModel(
                [rotation_features_except_005_rel_000, rotation_feature_zero_005_rel_000])
        
        #
        # Create a new 005-000 rotation feature such that the total net rotation is zero.
        #
        
        self.no_net_rotation_time_samples_005_rel_000 = []
        
        # If we're not starting at 0Ma then attempt to re-use the existing partial no-net-rotation file.
        # This can save a lot of time if we need to re-start an interrupted optimisation run.
        if end_age != 0:
            try:
                partial_no_net_rotation_features = pygplates.FeatureCollection(
                        os.path.join(self.data_dir, self.no_net_rotation_filename))
            except pygplates.OpenFileForReadingError:
                warnings.warn('Attempted to re-use partial no-net-rotation file {0} starting at {1} '
                    'but could not open for reading, so starting at 0Ma instead'.format(
                        self.no_net_rotation_filename, end_age))
                end_age = 0
        
        # Get the initial 005-000 samples.
        # If we've started a new optimisation run then this will be the identity rotation at 0Ma.
        # If we are continuing a previous partial optimisation run then this will be all previous
        # 005-00 no-net-rotations computed so far.
        if end_age != 0:
            print('Re-using existing partial no-net-rotation file from 0Ma to {0}Ma: {1}'.format(
                    end_age, self.no_net_rotation_filename))
            sys.stdout.flush()
            
            # Look for existing 005-000 partial no-net-rotations in existing no-net-rotation file.
            # We only collect those samples with times '<= end_age'.
            for rotation_feature in partial_no_net_rotation_features:
                total_reconstruction_pole = rotation_feature.get_total_reconstruction_pole()
                if total_reconstruction_pole:
                    fixed_plate_id, moving_plate_id, rotation_sequence = total_reconstruction_pole
                    if moving_plate_id == 5 and fixed_plate_id == 0:
                        for finite_rotation_sample in rotation_sequence.get_enabled_time_samples():
                            if finite_rotation_sample.get_time() <= pygplates.GeoTimeInstant(end_age):
                                self.no_net_rotation_time_samples_005_rel_000.append(finite_rotation_sample)
                        break
            if not self.no_net_rotation_time_samples_005_rel_000:
                raise RuntimeError('Expected 005-000 in existing partial no-net-rotation file {0} from 0-{1}Ma'.format(
                        self.no_net_rotation_filename, end_age))
            
            # Get the no-net-rotation at time 'end_age'.
            no_net_rotation_at_end_age = self.no_net_rotation_time_samples_005_rel_000[-1].get_value().get_finite_rotation()
            # The no-net-rotation was obtained from the inverse of the net rotation.
            # So to get the net rotation at time 'end_age' we take the inverse again (to undo the other inverse).
            self.net_total_rotation = no_net_rotation_at_end_age.get_inverse()
            
            # Keep track of how update-to-date we are.
            self.last_update_time = end_age
            
        else:
            # Start with identity rotation at time 0Ma.
            self.no_net_rotation_time_samples_005_rel_000.append(
                pygplates.GpmlTimeSample(pygplates.GpmlFiniteRotation(pygplates.FiniteRotation()), 0.0, 'NNR'))
            
            # Start with identity net rotation at time 0Ma.
            self.net_total_rotation = pygplates.FiniteRotation()
            
            # Keep track of how update-to-date we are.
            self.last_update_time = 0
        
        
        if self.CREATE_NO_NET_ROTATION_MODEL_AT_INIT:
            print('Calculating no-net-rotation model for {0}-{1}Ma...'.format(end_age, start_age))
            sys.stdout.flush()
            
            # Calculate no-net-rotation at 1My increments from 'end_age+1' Ma to 'start_age' Ma.
            # We shouldn't need to sample older than 'start_age'.
            for time in range(end_age + 1, start_age + 1):
                
                #print time
                #sys.stdout.flush()
                
                if self.gplates_net_rotations:
                    net_stage_lat, net_stage_lon, net_stage_angle_per_my = self.gplates_net_rotations[time]
                    net_stage_rotation = pygplates.FiniteRotation(
                            pygplates.PointOnSphere(net_stage_lat, net_stage_lon),
                            math.radians(net_stage_angle_per_my))
                else:
                    # Calculate net stage rotation from 'time' to 'time-1'.
                    net_stage_rotation = net_rotation.calculate_net_rotation_internal_gplates(
                            rotation_model_zero_005_rel_000,
                            topology_features,
                            time,
                            velocity_method = net_rotation.VelocityMethod.T_TO_T_MINUS_DT,
                            velocity_delta_time = 1.0)
                
                # Accumulate net stage rotations going backward in time (hence the inverse stage rotation)
                # since finite rotations go backward in time in the rotation file.
                self.net_total_rotation = net_stage_rotation.get_inverse() * self.net_total_rotation
                
                # Remove net total rotation at current time from 005 rel 000 rotation.
                # Note that removing from 005 affects all other plates (since every plate circuit goes through 005).
                #
                #   R_net_rotation * R_no_net_rotation(0->t,000->005) = Identity
                #                    R_no_net_rotation(0->t,000->005) = inverse(R_net_rotation)
                #
                no_net_rotation_005_rel_000 = self.net_total_rotation.get_inverse()
                
                self.no_net_rotation_time_samples_005_rel_000.append(
                    pygplates.GpmlTimeSample(
                            pygplates.GpmlFiniteRotation(no_net_rotation_005_rel_000),
                            time,
                            'NNR'))
            
            # Keep track of how update-to-date we are.
            self.last_update_time = start_age
            
            # Create a new no-net-rotation sequence.
            no_net_rotation_feature_005_rel_000 = pygplates.Feature.create_total_reconstruction_sequence(
                0,
                5,
                pygplates.GpmlIrregularSampling(self.no_net_rotation_time_samples_005_rel_000))
            
            # All no-net-rotation features. This includes the 005-000 no-net-rotation feature.
            all_no_net_rotation_features = list(rotation_features_except_005_rel_000)
            all_no_net_rotation_features.append(no_net_rotation_feature_005_rel_000)
            
            # Write the no-net-rotation file.
            pygplates.FeatureCollection(all_no_net_rotation_features).write(
                    os.path.join(self.data_dir, self.no_net_rotation_filename))
            
            print('...finished calculating no-net-rotation model.')
            sys.stdout.flush()
        
        # Extra data to keep track of for use in later update
        # (if haven't created entire no-net-rotation file already above).
        self.rotation_model_zero_005_rel_000 = rotation_model_zero_005_rel_000
        self.no_net_rotation_features_except_005_rel_000 = rotation_features_except_005_rel_000
    
    
    def get_no_net_rotation_filename(self):
        """
        Return the filename (relative to the 'data/' directory) of the single combined rotation file containing no net rotation.
        """
        return self.no_net_rotation_filename
    
    
    def update_no_net_rotation(
            self,
            ref_rotation_start_age):
        """
        Update the no-net-rotation 005-000 sequence from last update to 'ref_rotation_start_age+1'.
        Note that it's 'ref_rotation_start_age+1' instead of 'ref_rotation_start_age' because that's the oldest time
        required by 'optimisation_methods.ApproximateNR_from_features()' for the current time interval.
        
        The no-net-rotation model is then written to the file 'get_no_net_rotation_filename()'.
        
        NOTE: This is just temporary. Ideally we should just calculate the entire no-net-rotation model
              at startup (ie, for the entire time range). But we're dividing it up into time sections
              for now so that we don't have to wait a long time at startup before seeing things happen.
        """
        
        # If we're already update-to-date then nothing to do.
        #
        # This can happen if we've already created the entire no-net-rotation model in '__init__()'.
        #
        # This can also happen if we read a partial no-net-rotation file in '__init__()'
        # (because we're continuing a previously interrupted run, ie, 'end_age > 0').
        if ref_rotation_start_age <= self.last_update_time:
            return
        
        print('Calculating no-net-rotation model for {0}-{1}Ma...'.format(self.last_update_time, ref_rotation_start_age))
        sys.stdout.flush()
        
        # Calculate no-net-rotation at 1My increments from where we left off until 'ref_rotation_start_age' Ma.
        #
        # We shouldn't need to sample older than 'ref_rotation_start_age' because at each iteration the
        # candidate optimised rotation models (that are compared to the no-net-rotation model) are
        # only updated to 'ref_rotation_start_age' and hence sampling older than that would
        # create problems (ie, sample unknown rotations).
        for time in range(self.last_update_time + 1, ref_rotation_start_age + 1):
            
            #print time
            #sys.stdout.flush()
            
            if self.gplates_net_rotations:
                net_stage_lat, net_stage_lon, net_stage_angle_per_my = self.gplates_net_rotations[time]
                net_stage_rotation = pygplates.FiniteRotation(
                        pygplates.PointOnSphere(net_stage_lat, net_stage_lon),
                        math.radians(net_stage_angle_per_my))
            else:
                # Calculate net stage rotation from 'time' to 'time-1'.
                net_stage_rotation = net_rotation.calculate_net_rotation_internal_gplates(
                        self.rotation_model_zero_005_rel_000,
                        self.topology_features,
                        time,
                        velocity_method = net_rotation.VelocityMethod.T_TO_T_MINUS_DT,
                        velocity_delta_time = 1.0)
            
            # Accumulate net stage rotations going backward in time (hence the inverse stage rotation)
            # since finite rotations go backward in time in the rotation file.
            self.net_total_rotation = net_stage_rotation.get_inverse() * self.net_total_rotation
            
            # Remove net total rotation at current time from 005 rel 000 rotation.
            # Note that removing from 005 affects all other plates (since every plate circuit goes through 005).
            #
            #   R_net_rotation * R_no_net_rotation(0->t,000->005) = Identity
            #                    R_no_net_rotation(0->t,000->005) = inverse(R_net_rotation)
            #
            no_net_rotation_005_rel_000 = self.net_total_rotation.get_inverse()
            
            self.no_net_rotation_time_samples_005_rel_000.append(
                pygplates.GpmlTimeSample(
                        pygplates.GpmlFiniteRotation(no_net_rotation_005_rel_000),
                        time,
                        'NNR'))
        
        self.last_update_time = ref_rotation_start_age
        
        # Create a new no-net-rotation sequence.
        no_net_rotation_feature_005_rel_000 = pygplates.Feature.create_total_reconstruction_sequence(
            0,
            5,
            pygplates.GpmlIrregularSampling(self.no_net_rotation_time_samples_005_rel_000))
        
        # All no-net-rotation features. This includes the 005-000 no-net-rotation feature.
        all_no_net_rotation_features = list(self.no_net_rotation_features_except_005_rel_000)
        all_no_net_rotation_features.append(no_net_rotation_feature_005_rel_000)
        
        # Write the no-net-rotation file.
        pygplates.FeatureCollection(all_no_net_rotation_features).write(
                os.path.join(self.data_dir, self.no_net_rotation_filename))
        
        print('...finished calculating no-net-rotation model.')
        sys.stdout.flush()
    
    
    def _read_gplates_net_rotations(
            self,
            gplates_net_rotations_filename):
        """
        Temporary function to read the net rotations file exported by GPlates.
        
        TODO: Remove this when we can calculate net rotation for deforming models using pygplates.
        """
        
        net_rotations = {}
        with open(gplates_net_rotations_filename, 'r') as net_rotations_file:
            for line_number, line in enumerate(net_rotations_file):

                # Make line number 1-based instead of 0-based.
                line_number = line_number + 1

                line = line.strip()
                
                # Skip empty lines.
                if not line:
                    continue
            
                # Split the line into strings (separated by commas).
                line_string_list = line.split(',')

                # Expecting 4 numbers
                if len(line_string_list) != 4:
                    continue

                # Attempt to convert each string into a floating-point number.
                try:
                    time = float(line_string_list[0])
                    net_stage_lat = float(line_string_list[1])
                    net_stage_lon = float(line_string_list[2])
                    net_stage_angle_per_my = float(line_string_list[3])
                except ValueError:
                    continue

                net_rotations[time] = net_stage_lat, net_stage_lon, net_stage_angle_per_my
        
        return net_rotations
