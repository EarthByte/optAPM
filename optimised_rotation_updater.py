import numpy as np
import os.path
import pygplates
import sys
import warnings

from ptt import remove_plate_rotations


class OptimisedRotationUpdater(object):
    """
    Class to manage updates to the rotation model due to optimisation.
    
    The original input rotation files are combined into an intermediate single optimised rotation file,
    returned by 'get_optimised_rotation_filename()' (relative to the 'data/' directory).
    
    Each optimised rotation 005-000 at each time is updated using 'update_optimised_rotation()'
    which writes the updates back out to the single optimised rotation file so that parallel
    processes can read it in for optimisation over the next time interval (pygplates,
    being a C++ extension, does not support pickling its objects for sending to parallel processes,
    hence the pygplates rotation features need to be saved and loaded from files).
    
    Finally the original rotation files are modified to reflect the optimisation
    using 'save_to_rotation_files()'.
    """
    
    def __init__(
            self,
            data_dir,
            original_rotation_filenames,  # Relative to the 'data/' directory.
            start_age,
            end_age,
            interval,
            reference_params_function,
            data_model,
            model_name):
        """
        Create a single optimised rotation file by combining all original (input) rotations.
        The 005-000 rotation feature is inserted (or replaced if already existing in input) and
        defined such that the rotation of reference plate (obtained for each time using
        'reference_params_function') relative to 000 is zero for each time from 'start_age' to
        'end_age + interval' in 'interval' steps.
        """
        
        self.data_dir = data_dir
        self.reference_params_function = reference_params_function
        self.model_name = model_name
        
        # First check that Africa has a zero present-day rotation.
        original_rotation_model = pygplates.RotationModel([os.path.join(self.data_dir, original_rotation_filename)
            for original_rotation_filename in original_rotation_filenames])
        if not original_rotation_model.get_rotation(0.0, 701).represents_identity_rotation():
            # If an exception is raised here then the original rotation model needs to be fixed before running the optimisation workflow.
            raise RuntimeWarning('Original rotation model has a non-zero finite rotation for Africa at present day.')
        
        # The single combined optimised rotation filename (relative to the 'data/' directory).
        self.optimised_rotation_filename = os.path.join(
                data_model, 'optimisation', 'optimised_rotation_model_' + self.model_name + '.rot')
        
        #
        # Combine the original (input) rotation files into a single optimised version, and
        # zero out the 005-000 rotations that we will soon replace with optimised rotations.
        #
        
        # Load all the original rotation feature collections.
        self.rotation_filenames = []
        self.rotation_feature_collections = []
        for rotation_filename in original_rotation_filenames:
            # Read the current rotation file.
            rotation_feature_collection = pygplates.FeatureCollection(
                    os.path.join(self.data_dir, rotation_filename))
            # Convert from a FeatureCollection to a list of Feature.
            rotation_feature_collection = list(rotation_feature_collection)
            # Keep track of original feature collections (these will get modified).
            self.rotation_feature_collections.append(rotation_feature_collection)
            # Also keep track of the original rotation filenames so
            # we can write back to them when we're finished optimising.
            self.rotation_filenames.append(rotation_filename)
        
        # Remove any existing 005-000 rotation features and change fixed plate 000 to 005.
        for rotation_feature_collection in self.rotation_feature_collections:
            rotation_feature_index = 0
            while rotation_feature_index < len(rotation_feature_collection):
                rotation_feature = rotation_feature_collection[rotation_feature_index]
                total_reconstruction_pole = rotation_feature.get_total_reconstruction_pole()
                if total_reconstruction_pole:
                    fixed_plate_id, moving_plate_id, rotation_sequence = total_reconstruction_pole
                    
                    # Change fixed plate ID 000 to 005 (unless it's the 005-000 plate pair).
                    # We want everything that references 000 to now reference 005.
                    # Later we'll add a 005-000 sequence to store optimised rotation adjustments.
                    if fixed_plate_id == 0:
                        if moving_plate_id == 5:
                            # Exclude 005-000 rotation features (we'll add a new one later).
                            del rotation_feature_collection[rotation_feature_index]
                            rotation_feature_index -= 1
                        else:
                            # Change the fixed plate ID from 000 to 005.
                            rotation_feature.set_total_reconstruction_pole(5, moving_plate_id, rotation_sequence)
                
                rotation_feature_index += 1
        
        # Rotation model with root plate 005.
        #
        # Excludes 005-000 rotation features (if any) and root plate is now 005 (instead of 000).
        original_rotation_model_anchor_005 = pygplates.RotationModel(self.rotation_feature_collections)
        
        #
        # Create a new 005-000 rotation feature such that 'ref_rotation_plate_id' rel 000 is zero (from 'start_age' to 'end_age').
        #
        
        rotation_time_samples_005_rel_000 = []
        
        # If we're not starting at 0Ma then attempt to re-use the existing partially optimised rotation file.
        # This can save a lot of time if we need to re-start an interrupted optimisation run.
        if end_age != 0:
            try:
                partially_optimised_rotation_features = pygplates.FeatureCollection(
                        os.path.join(self.data_dir, self.optimised_rotation_filename))
            except pygplates.OpenFileForReadingError:
                warnings.warn('Attempted to re-use partially optimised rotation file {0} starting at {1} '
                    'but could not open for reading, so starting at 0Ma instead'.format(
                        self.optimised_rotation_filename, end_age))
                end_age = 0
        
        # Get the initial 005-000 samples.
        # If we've started a new optimisation run then this will only be the identity rotation at 0Ma.
        # If we are continuing a previous partial optimisation run then this will be all previous
        # 005-00 rotations computed so far.
        if end_age != 0:
            print 'Re-using existing partially optimised rotation file from 0Ma to {0}Ma'.format(end_age)
            sys.stdout.flush()
            
            # Look for existing 005-000 partially optimised rotation in existing optimised rotation file.
            # We only collect those samples with times '<= end_age'.
            for rotation_feature in partially_optimised_rotation_features:
                total_reconstruction_pole = rotation_feature.get_total_reconstruction_pole()
                if total_reconstruction_pole:
                    fixed_plate_id, moving_plate_id, rotation_sequence = total_reconstruction_pole
                    if moving_plate_id == 5 and fixed_plate_id == 0:
                        for finite_rotation_sample in rotation_sequence.get_enabled_time_samples():
                            if finite_rotation_sample.get_time() <= pygplates.GeoTimeInstant(end_age):
                                rotation_time_samples_005_rel_000.append(finite_rotation_sample)
                        break
            if not rotation_time_samples_005_rel_000:
                raise RuntimeError('Expected 005-000 in existing partially optimised file {0} from 0-{1}Ma'.format(
                        self.optimised_rotation_filename, end_age))
        else:
            # Start with identity rotation at time 0Ma.
            rotation_time_samples_005_rel_000.append(
                pygplates.GpmlTimeSample(pygplates.GpmlFiniteRotation(pygplates.FiniteRotation()), 0.0, 'optAPM @absage'))
        
        # Add a rotation at each age in 'end_age + interval' to 'start_age' (inclusive) at 'interval' steps.
        for ref_rotation_start_age in np.arange(end_age + interval, start_age + interval, interval):
            
            # Get the reference plate ID (which could vary over time).
            ref_rotation_plate_id, _ = self.reference_params_function(ref_rotation_start_age)
            
            # We want our 'ref_rotation_plate_id' to 000 rotation to be zero.
            #
            #                 R(0->t,000->ref_plate) = R(0->t,000->005) * R(0->t,005->ref_plate)
            #                               Identity = R(0->t,000->005) * R(0->t,005->ref_plate)
            #   inverse(R(0->t,000->005)) * Identity = R(0->t,005->ref_plate)
            #                       R(0->t,000->005) = inverse(R(0->t,005->ref_plate))
            zero_rotation_005_rel_000 = original_rotation_model_anchor_005.get_rotation(
                    ref_rotation_start_age, ref_rotation_plate_id, anchor_plate_id=5).get_inverse()

            rotation_time_samples_005_rel_000.append(
                pygplates.GpmlTimeSample(
                        pygplates.GpmlFiniteRotation(zero_rotation_005_rel_000),
                        ref_rotation_start_age,
                        'optAPM @absage'))
        
        # Create a new 005/000 rotation sequence.
        # And keep track of 005-000 feature for later so we can update it as we optimise through time.
        self.optimised_rotation_feature = pygplates.Feature.create_total_reconstruction_sequence(
            0,
            5,
            pygplates.GpmlIrregularSampling(rotation_time_samples_005_rel_000))
        
        # Keep track of all optimised rotation features for later so we can write them to file
        # as we optimise through time. This includes the optimised 005-000 rotation feature.
        all_optimised_rotation_features = []
        for rotation_feature_collection in self.rotation_feature_collections:
            all_optimised_rotation_features.extend(rotation_feature_collection)
        all_optimised_rotation_features.append(self.optimised_rotation_feature)
        self.optimised_rotation_feature_collection = pygplates.FeatureCollection(all_optimised_rotation_features)
        
        # Write the rotation file with zero reference_plate-to-anchor rotations.
        self.optimised_rotation_feature_collection.write(
                os.path.join(self.data_dir, self.optimised_rotation_filename))
    
    
    def get_optimised_rotation_filename(self):
        """
        Return the filename (relative to the 'data/' directory) of the single combined optimised rotation file.
        """
        return self.optimised_rotation_filename
    
    
    def update_optimised_rotation(
            self,
            optimised_rotation_ref_plate_rel_000,
            ref_rotation_plate_id,
            ref_rotation_start_age):
        """
        Update the optimised 005-000 sequence at time 'ref_rotation_start_age', and
        write the updated rotation model to 'get_optimised_rotation_filename()'.
        """
        
        # The rotation model prior to updating the optimised rotation at the specified time.
        rotation_model = pygplates.RotationModel(self.optimised_rotation_feature_collection)
        
        # The optimised rotation *sequence* into which we'll insert the optimised rotation *pole* at specified time.
        _, _, optimised_rotation_sequence = self.optimised_rotation_feature.get_total_reconstruction_pole()
        
        # Update the optimised rotation sequence with optimised rotation pole at specified time.
        for optimised_finite_rotation_sample in optimised_rotation_sequence.get_enabled_time_samples():
            optimised_finite_rotation_time = optimised_finite_rotation_sample.get_time()
            if optimised_finite_rotation_time == ref_rotation_start_age:
                # Our optimised rotation is from 'ref_rotation_plate_id' to 000 so remove the
                # 'ref_rotation_plate_id' to 005 part to get the 005 to 000 part that gets
                # stored in the 005-000 rotation feature.
                #
                #                                     R(0->t,000->ref_plate) = R(0->t,000->005) * R(0->t,005->ref_plate)
                #   R(0->t,000->ref_plate) * inverse(R(0->t,005->ref_plate)) = R(0->t,000->005)
                #
                plate_rotation_ref_plate_rel_005 = rotation_model.get_rotation(
                        ref_rotation_start_age, ref_rotation_plate_id, fixed_plate_id=5)
                optimised_rotation_005_rel_000 = optimised_rotation_ref_plate_rel_000 * plate_rotation_ref_plate_rel_005.get_inverse()
                
                optimised_finite_rotation_sample.get_value().set_finite_rotation(optimised_rotation_005_rel_000)
        
        # Write all the rotation features (including the optimised 005-000 features) to rotation file.
        # This will write out the changes to the optimised rotation *sequence* we just made above.
        self.optimised_rotation_feature_collection.write(
                os.path.join(self.data_dir, self.optimised_rotation_filename))
    
    
    def save_to_rotation_files(self):
        """
        The original rotation files are modified to reflect the optimisation.
        
        Any rotations that referenced plate 000 are modified to include the optimised absolute plate motion.
        All other rotations remain unmodified.
        """
        
        # Remove the optimised rotation 005-000 and merge it into any rotations with fixed plate 005.
        #
        # First add the optimised rotation 005-000 to any existing feature collection so we don't have
        # to add one extra feature collection just to store the optimised rotation 005-000.
        # (it'll then get removed by 'remove_plates()').
        self.rotation_feature_collections[0].append(self.optimised_rotation_feature)
        self.rotation_feature_collections = remove_plate_rotations.remove_plates(
                self.rotation_feature_collections,
                [5])
        
        output_filename_suffix = '_' + self.model_name
        
        # Write the modified rotation feature collections back to disk.
        for rotation_feature_collection_index in range(len(self.rotation_feature_collections)):
            output_rotation_feature_collection = self.rotation_feature_collections[rotation_feature_collection_index]
            
            # Each output filename is the input filename with an optional suffix appended (before the extension).
            input_rotation_filename = self.rotation_filenames[rotation_feature_collection_index]
            if output_filename_suffix:
                dir, file_basename = os.path.split(input_rotation_filename)
                file_basename, file_ext = os.path.splitext(file_basename)
                output_rotation_filename = os.path.join(dir, 'optimisation', '{0}{1}{2}'.format(file_basename, output_filename_suffix, file_ext))
            else:
                output_rotation_filename = input_rotation_filename
            
            output_rotation_filename = os.path.join(self.data_dir, output_rotation_filename)
            
            output_rotation_feature_collection.write(output_rotation_filename)
