import os.path
import pygplates


class OptimisedRotationUpdater(object):
    """
    Class to manage updates to the rotation model due to optimisation.
    
    The unoptimised input rotation files are combined into a single optimised rotation file,
    returned by 'get_optimised_rotation_filename()' (relative to the 'data/' directory).
    
    Each optimised rotation 005-000 at each time is updated using 'update_optimised_rotation()'
    which writes the updates back out to the single optimised rotation file so that parallel
    processes can read it in for optimisation over the next time interval (pygplates,
    being a C++ extension, does not support pickling its objects for sending to parallel processes,
    hence the pygplates rotation features need to be saved and loaded from files).
    """
    
    def __init__(
            self,
            data_dir,
            unoptimised_rotation_filenames,  # Relative to the 'data/' directory.
            age_range,
            reference_params_function,
            data_model,
            model_name):
        """
        Create a single optimised rotation file by combining all unoptimised (input) rotations.
        The 005-000 rotation feature is inserted (or replaced if already existing in input) and
        defined such that the rotation of reference plate (obtained for each time using
        'reference_params_function') relative to 000 is zero for each time in 'age_range'.
        """
        
        self.data_dir = data_dir
        self.reference_params_function = reference_params_function
        
        # The single combined optimised rotation filename (relative to the 'data/' directory).
        self.optimised_rotation_filename = os.path.join(
                data_model, 'optimisation', 'optimised_rotations_' + model_name + '.rot')
        
        #
        # Combine the unoptimised (input) rotation files into a single optimised version, and
        # zero out the 005-000 rotations that we will soon replace with optimised rotations.
        #
        
        # Load all the unoptimised rotation features.
        rotation_features = []
        for unoptimised_rotation_filename in unoptimised_rotation_filenames:
            rotation_feature_collection = pygplates.FeatureCollection(
                    os.path.join(self.data_dir, unoptimised_rotation_filename))
            rotation_features.extend(rotation_feature_collection)
        
        # Remove any existing 005-000 rotation features and change fixed plate 000 to 005.
        rotation_features_except_005_000 = []
        for rotation_feature in rotation_features:
            total_reconstruction_pole = rotation_feature.get_total_reconstruction_pole()
            if total_reconstruction_pole:
                fixed_plate_id, moving_plate_id, rotation_sequence = total_reconstruction_pole
                
                # Change fixed plate ID 000 to 005 (unless it's the 005-000 plate pair).
                # We want everything that references 000 to now reference 005.
                # Later we'll add a 005-000 sequence to store optimised rotation adjustments.
                if fixed_plate_id == 0:
                    if moving_plate_id == 5:
                        # Exclude 005-000 rotation features (we'll add a new one below).
                        continue
                    else:
                        # Change the fixed plate ID from 000 to 005.
                        rotation_feature.set_total_reconstruction_pole(5, moving_plate_id, rotation_sequence)
            
            # All existing rotation features get added except 005-000.
            rotation_features_except_005_000.append(rotation_feature)
        
        # Rotation features now exclude any old 005-000 rotation features.
        rotation_features = rotation_features_except_005_000
        
        # Unoptimised rotation model.
        # Excludes 005-000 rotation features (if any) and root plate is now 005 (instead of 000).
        unoptimised_rotation_model_anchor_005 = pygplates.RotationModel(rotation_features)
        
        #
        # Create a new 005-000 rotation feature such that 'ref_rotation_plate_id' rel 000 is zero.
        #
        
        zero_rotation_time_samples_005_rel_000 = []
        
        # Start with identity rotation at time 0Ma.
        zero_rotation_time_samples_005_rel_000.append(
            pygplates.GpmlTimeSample(pygplates.GpmlFiniteRotation(pygplates.FiniteRotation()), 0.0, 'optAPM'))
        
        for ref_rotation_start_age in age_range:
            
            # Get the reference plate ID (which could vary over time).
            ref_rotation_plate_id, _ = self.reference_params_function(ref_rotation_start_age)
            
            # We want our 'ref_rotation_plate_id' to 000 rotation to be zero.
            #
            #                 R(0->t,000->ref_plate) = R(0->t,000->005) * R(0->t,005->ref_plate)
            #                               Identity = R(0->t,000->005) * R(0->t,005->ref_plate)
            #   inverse(R(0->t,000->005)) * Identity = R(0->t,005->ref_plate)
            #                       R(0->t,000->005) = inverse(R(0->t,005->ref_plate))
            zero_rotation_005_rel_000 = unoptimised_rotation_model_anchor_005.get_rotation(
                    ref_rotation_start_age, ref_rotation_plate_id, anchor_plate_id=5).get_inverse()

            zero_rotation_time_samples_005_rel_000.append(
                pygplates.GpmlTimeSample(
                        pygplates.GpmlFiniteRotation(zero_rotation_005_rel_000),
                        ref_rotation_start_age,
                        'optAPM'))
        
        # Create a new 005/000 rotation sequence.
        rotation_feature_005_000 = pygplates.Feature.create_total_reconstruction_sequence(
            0,
            5,
            pygplates.GpmlIrregularSampling(zero_rotation_time_samples_005_rel_000))
        # Keep track of 005-000 feature for later so we can update it as we optimise through time.
        self.optimised_rotation_feature = rotation_feature_005_000
        
        # Add the 005-000 rotation feature.
        rotation_features.append(rotation_feature_005_000)
        
        # Keep track of all rotation features for later so we can write them to file as we optimise through time.
        # This includes the optimised 005-000 rotation feature.
        self.optimised_rotation_feature_collection = pygplates.FeatureCollection(rotation_features)
        
        # Write the rotation file with zero reference_plate-to-anchor rotations.
        self.optimised_rotation_feature_collection.write(
                os.path.join(self.data_dir, self.optimised_rotation_filename))
    
    
    def get_optimised_rotation_filename(self):
        """
        Return the filename of the single combined optimised rotation file (relative to the 'data/' directory).
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
