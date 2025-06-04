import glob
import os.path
import pygplates


#### Input parameters ####

# Original (un-optimised) rotation files.
original_rotation_filenames = [
    '1000_0_rotfile.rot',
    '1800_1000_rotfile.rot',
]
#original_rotation_filenames = glob.glob('*.rot')

# This file contains ALL rotations, including the 005-000 optimised reference frame.
optimised_rotation_file = 'optimisation/optimised_rotation_model_20240725.rot'

# This file contains ALL rotations, including the 005-000 no-net-rotation reference frame.
no_net_rotation_file = 'optimisation/no_net_rotation_model_20240725.rot'

# Output rotation file to contain the reference frames being added (optimised, NNR, optimised-minus-pmag).
output_ref_frames_file = 'optimisation/reference_frames_20240725.rot'

# Plate IDs of reference frames to store in output rotation file.
optimised_ref_frame_plate_id = 6
no_net_rotation_ref_frame_plate_id = 7
optimised_minus_pmag_ref_frame_plate_id = 8

##########################


optimised_rotation_feature_collection = pygplates.FeatureCollection(optimised_rotation_file)
no_net_rotation_feature_collection = pygplates.FeatureCollection(no_net_rotation_file)

# Original (un-optimised) rotation model.
original_rotation_model = pygplates.RotationModel(original_rotation_filenames)

def find_005_000_rotation_sequence(rotation_feature_collection):
    rotation_feature_005_000 = None

    for rotation_feature in rotation_feature_collection:
        total_reconstruction_pole = rotation_feature.get_total_reconstruction_pole()
        if total_reconstruction_pole:
            fixed_plate_id, moving_plate_id, _ = total_reconstruction_pole
            
            if fixed_plate_id == 0 and moving_plate_id == 5:
                if rotation_feature_005_000 is not None:
                    raise ValueError("Found multiple rotation sequences with moving plate 5 and fixed plate 0")
                rotation_feature_005_000 = rotation_feature

    if rotation_feature_005_000 is None:
        raise ValueError("Found zero rotation sequences with moving plate 5 and fixed plate 0")

    return rotation_feature_005_000

# Find the optimised (005-000) rotation sequence.
optimised_rotation_feature = find_005_000_rotation_sequence(optimised_rotation_feature_collection)
optimised_rotation_feature = optimised_rotation_feature.clone()
# Change its moving plate ID (from 005 to <optimised_ref_frame_plate_id>).
optimised_rotation_feature.set(pygplates.PropertyName.gpml_moving_reference_frame, pygplates.GpmlPlateId(optimised_ref_frame_plate_id))
# Reverse its rotations since it'll get activated by setting the anchor plate to <optimised_ref_frame_plate_id> instead of 000.
for rotation_sample in optimised_rotation_feature.get_total_reconstruction_pole()[2].get_enabled_time_samples():
    finite_rotation = rotation_sample.get_value().get_finite_rotation()
    rotation_sample.get_value().set_finite_rotation(finite_rotation.get_inverse())
    rotation_sample.set_description('optAPM @absage')

# Find the no-net rotation (005-000) rotation sequence.
no_net_rotation_feature = find_005_000_rotation_sequence(no_net_rotation_feature_collection)
no_net_rotation_feature = no_net_rotation_feature.clone()
# Change its moving plate ID (from 005 to <no_net_rotation_ref_frame_plate_id>).
no_net_rotation_feature.set(pygplates.PropertyName.gpml_moving_reference_frame, pygplates.GpmlPlateId(no_net_rotation_ref_frame_plate_id))
# Reverse its rotations since it'll get activated by setting the anchor plate to <no_net_rotation_ref_frame_plate_id> instead of 000.
for rotation_sample in no_net_rotation_feature.get_total_reconstruction_pole()[2].get_enabled_time_samples():
    finite_rotation = rotation_sample.get_value().get_finite_rotation()
    rotation_sample.get_value().set_finite_rotation(finite_rotation.get_inverse())
    rotation_sample.set_description('NNR')

# Calculate a new optimised-minus-pmag sequence that combines:
# - a reversal of 701-000 (to remove PMAG frame), and
# - a reversal of existing optimised rotation (005-000)
#   - since it'll get activated by setting the anchor plate to <optimised_minus_pmag_ref_frame_plate_id> instead of 000.
rotation_time_samples_optimised_minus_pmag = []
for optimised_rotation_sample in optimised_rotation_feature.get_total_reconstruction_pole()[2].get_enabled_time_samples():
    time = optimised_rotation_sample.get_time()
    optimised_finite_rotation = optimised_rotation_sample.get_value().get_finite_rotation()

    # Get the 701-000 rotation.
    finite_rotation_701_000 = original_rotation_model.get_rotation(time, 701)

    # Remove pmag from optimised reference frame:
    #
    #  R(opt_minus_pmag->000)  = R(opt->000) * inverse[R(000->701)]
    #                          = inverse[R(000->opt)] * inverse[R(000->701)]
    #  R(000->opt_minus_pmag)) = inverse[R(opt_minus_pmag->000)]
    #                          = inverse[inverse[R(000->opt)] * inverse[R(000->701)]]
    finite_rotation_optimised_minus_pmag = (optimised_finite_rotation.get_inverse() * finite_rotation_701_000.get_inverse()).get_inverse()

    rotation_time_samples_optimised_minus_pmag.append(
        pygplates.GpmlTimeSample(
                pygplates.GpmlFiniteRotation(finite_rotation_optimised_minus_pmag),
                time,
                'optAPM without pmag'))

# Create a new optimised-minus-pmag rotation sequence.
optimised_minus_pmag_rotation_feature = pygplates.Feature.create_total_reconstruction_sequence(
    0,
    optimised_minus_pmag_ref_frame_plate_id,
    pygplates.GpmlIrregularSampling(rotation_time_samples_optimised_minus_pmag))

# Add reference frames to feature collection.
output_ref_frames_feature_collection = pygplates.FeatureCollection()
output_ref_frames_feature_collection.add(optimised_rotation_feature)
output_ref_frames_feature_collection.add(no_net_rotation_feature)
output_ref_frames_feature_collection.add(optimised_minus_pmag_rotation_feature)

# Write reference frames to output rotation file.
output_ref_frames_feature_collection.write(output_ref_frames_file)
