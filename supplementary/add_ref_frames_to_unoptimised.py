import glob
import os.path
import pygplates


#
# This script reads the original un-optimised rotation files and the optimised and
# the no-net rotation reference frames, and generates a new set of rotation files
# with various reference frames added.
#
# Currently the added reference frames are the original paleomag (PMAG), the optimised mantle,
# the no-net rotation, and an approximation of true polar wander (TPW) calculated by replacing
# the PMAG reference frame with the optimised mantle reference frame.
#


#### Input parameters ####

# Original (un-optimised) rotation files.
input_rotation_filenames = [
    '1000_0_rotfile.rot',
    '1800_1000_rotfile.rot',
]
#input_rotation_filenames = glob.glob('*.rot')

# This file contains the 005-000 optimised reference frame.
input_optimised_rotation_file = 'optimisation/optimised_rotation_model_20240725.rot'

# This file contains the 005-000 no-net-rotation reference frame.
input_no_net_rotation_file = 'optimisation/no_net_rotation_model_20240725.rot'

# The suffix to append to the basenames of the input rotation filenames to get the output rotation filenames.
#
# Note: The output rotation files also end up in the 'optimisation/' sub-directory.
output_rotation_filenames_suffix = '_20240725'

# Plate IDs of reference frames.
# 
# These are PMAG, optimised, no-net rotation and true polar wander (which is optimised with PMAG removed).
pmag_ref_frame_plate_id = 14
optimised_ref_frame_plate_id = 15
no_net_rotation_ref_frame_plate_id = 16
true_polar_wander_ref_frame_plate_id = 17

# The default reference frame that plate 000 should correspond to.
#
# It should be one of the above reference frame plate IDs.
default_reference_frame_plate_id = optimised_ref_frame_plate_id

#
# Note that the final rotation hierarchy will look something like...
#
#       015  016  017
#        |    |    |
#        -----------
#             |
#            000
#             |
#            014
#             |
#           -----
#           |   |
#          101 701
#           |   |
#
# ...where if 000 represents:
#   1. PMAG – then 014-000 is zero, and 000-015 (and 014-015) represents what 005-000 previously represented.
#   2. Optimised – then 000-015 is zero, and 014-000 (and 014-015) represents what 005-000 previously represented.
#
# That can be one version of the rotation model.
# And another can be a flattened version that removes the reference frame plate IDs (eg, 014, 015, 016, 017).
# The flattened version can be what's included with GPlately.
#

##########################


# Each input rotation file will get a corresponding output rotation file.
output_rotation_feature_collections = [pygplates.FeatureCollection(filename)
                                       for filename in input_rotation_filenames]

# The output rotation filenames associated with the input rotation files.
output_rotation_filenames = []
for filename in input_rotation_filenames:
    # Insert suffix into output rotation filenames.
    basename, ext = os.path.splitext(filename)
    output_filename = basename + output_rotation_filenames_suffix + ext
    # Place in 'optimisation/' sub-directory.
    output_filename = os.path.join('optimisation', output_filename)
    output_rotation_filenames.append(output_filename)

# The output rotation filename to store the optimised, NNR and true polar wander reference frames.
output_reference_frames_filename = 'optimisation/reference_frames{}.rot'.format(output_rotation_filenames_suffix)

##########################


# Change fixed plate 000 to <pmag_ref_frame_plate_id>
# (and remove any existing leftover 005-000 rotation features).
max_time_of_sequences_with_fixed_000 = 0.0
for rotation_feature_collection in output_rotation_feature_collections:

    # Exclude any leftover 005-000 rotation features.
    def is_feature_005_000(feature):
        total_reconstruction_pole = feature.get_total_reconstruction_pole()
        if total_reconstruction_pole:
            fixed_plate_id, moving_plate_id, _ = total_reconstruction_pole
            return fixed_plate_id == 0 and moving_plate_id == 5
        return False
    rotation_feature_collection.remove(is_feature_005_000)
    
    # Change fixed plate 000 to <pmag_ref_frame_plate_id>.
    for rotation_feature in rotation_feature_collection:
        total_reconstruction_pole = rotation_feature.get_total_reconstruction_pole()
        if total_reconstruction_pole:
            fixed_plate_id, moving_plate_id, rotation_sequence = total_reconstruction_pole
            
            # Change fixed plate ID 000 to <pmag_ref_frame_plate_id>.
            # We want everything that references 000 to now reference <pmag_ref_frame_plate_id>.
            if fixed_plate_id == 0:
                # Change the fixed plate ID from 000 to <pmag_ref_frame_plate_id>.
                rotation_feature.set_total_reconstruction_pole(
                    pmag_ref_frame_plate_id, moving_plate_id, rotation_sequence)
                # The oldest time in this sequence.
                max_time = max(time_sample.get_time() for time_sample in rotation_sequence.get_time_samples())
                # The oldest time in all sequences (with fixed plate 000).
                max_time_of_sequences_with_fixed_000 = max(max_time_of_sequences_with_fixed_000, max_time)

# Original (un-optimised) rotation model.
#
# Note: We use the modified input rotation features (with fixed plate 000 changed to 014).
#
# Note: It has no 000 plate ID. Instead it's anchor plate is the PMAG reference frame plate ID.
original_rotation_model_anchored_014 = pygplates.RotationModel(
    output_rotation_feature_collections,
    default_anchor_plate_id=pmag_ref_frame_plate_id)


def find_005_000_rotation_sequence(rotation_filename):
    rotation_sequence_005_000 = None

    for rotation_feature in pygplates.FeatureCollection(rotation_filename):
        total_reconstruction_pole = rotation_feature.get_total_reconstruction_pole()
        if total_reconstruction_pole:
            fixed_plate_id, moving_plate_id, rotation_sequence = total_reconstruction_pole
            
            if fixed_plate_id == 0 and moving_plate_id == 5:
                if rotation_sequence_005_000 is not None:
                    raise ValueError("Found multiple rotation sequences with moving plate 5 and fixed plate 0")
                rotation_sequence_005_000 = rotation_sequence

    if rotation_sequence_005_000 is None:
        raise ValueError("Found zero rotation sequences with moving plate 5 and fixed plate 0")

    return rotation_sequence_005_000

# Find the optimised (005-000) rotation sequence.
optimised_rotation_sequence = find_005_000_rotation_sequence(input_optimised_rotation_file)

# Find the no-net rotation (005-000) rotation sequence.
no_net_rotation_sequence = find_005_000_rotation_sequence(input_no_net_rotation_file)


def create_identity_rotation_sequence(rotation_description, begin_time, end_time=0.0):
    rotation_time_samples = []
    identity_rotation = pygplates.FiniteRotation()
    # Only need two time samples at start and end times.
    for time in (end_time, begin_time):
        rotation_time_samples.append(
            pygplates.GpmlTimeSample(
                    pygplates.GpmlFiniteRotation(identity_rotation),
                    time,
                    rotation_description))
    return rotation_time_samples

def clone_rotation_sequence(rotation_sequence, rotation_description):
    rotation_time_samples = []
    for rotation_sample in rotation_sequence:
        rotation_time_samples.append(
            pygplates.GpmlTimeSample(
                    pygplates.GpmlFiniteRotation(rotation_sample.get_value().get_finite_rotation()),
                    rotation_sample.get_time(),
                    rotation_description,
                    rotation_sample.is_enabled()))
    return rotation_time_samples

def calculate_true_polar_wander_rotation_sequence(optimised_rotation_sequence, rotation_description):
    rotation_time_samples = []
    for optimised_rotation_sample in optimised_rotation_sequence:
        time = optimised_rotation_sample.get_time()
        # R(000->005).
        optimised_finite_rotation = optimised_rotation_sample.get_value().get_finite_rotation()
        # R(014->701).
        finite_rotation_701_014 = original_rotation_model_anchored_014.get_rotation(time, 701)
        #
        # PMAG contain a mantle reference frame and true polar wander (TPW).
        # Whereas OPT only contains a mantle reference frame.
        #
        #   TPW + OPT = PMAG
        #   TPW       = PMAG - OPT
        #
        # So looking at Africa this becomes...
        #
        #   R(PMAG->701) - R(OPT->701) = R(PMAG->701) * R(701->OPT)
        #                              = R(PMAG->OPT)
        #
        # ...and, using the old optimised rotation R(000->005) (where 000 was OPT and 005 was PMAG)...
        #
        #   R(PMAG->OPT) = R(005->000)
        #                = inverse[R(000->005)]
        #
        # ...therefore...
        #
        #   R(PMAG->701) - R(OPT->701) = R(PMAG->OPT)
        #                              = inverse[R(000->005)]
        #
        # So we want 017->701 to be equivalent to the *inverse* of the optimised rotation 000->005:
        #
        #   TPW         = PMAG - OPT
        #
        #   R(017->701) = R(PMAG->701) - R(OPT->701)
        #               = R(PMAG->OPT)
        #               = inverse[R(000->005)]
        #
        # And expressing this relative to 014 becomes...
        #
        #   R(017->701)               = inverse[R(000->005)]
        #   R(017->014) * R(014->701) = inverse[R(000->005)]
        #   R(017->014)               = inverse[R(000->005)] * inverse[R(014->701)]
        #
        finite_rotation_true_polar_wander = optimised_finite_rotation.get_inverse() * finite_rotation_701_014.get_inverse()
        rotation_time_samples.append(
            pygplates.GpmlTimeSample(
                    pygplates.GpmlFiniteRotation(finite_rotation_true_polar_wander),
                    time,
                    rotation_description,
                    optimised_rotation_sample.is_enabled()))
    return rotation_time_samples

# Generate the rotation sequence for the PMAG reference frame.
#
# This depends on which reference frame 000 is assigned to.
if default_reference_frame_plate_id == pmag_ref_frame_plate_id:
    # Both 000 and 014 are PMAG.
    #
    # So 000->014 is zero (identity).
    rotation_time_samples_pmag = create_identity_rotation_sequence(
        f' Reference frames: paleomag ({pmag_ref_frame_plate_id:03d} and 000)',
        max_time_of_sequences_with_fixed_000)
elif default_reference_frame_plate_id == optimised_ref_frame_plate_id:
    # 000 is optimised (and 014 is PMAG).
    #
    # So 000->014 is the existing optimised rotation (000->005).
    rotation_time_samples_pmag = clone_rotation_sequence(
        optimised_rotation_sequence,
        f' Reference frames: paleomag ({pmag_ref_frame_plate_id:03d}) and optimised mantle (000)')
elif default_reference_frame_plate_id == no_net_rotation_ref_frame_plate_id:
    # 000 is no-net rotation (and 014 is PMAG).
    #
    # So 000->014 is the existing no-net rotation (000->005).
    rotation_time_samples_pmag = clone_rotation_sequence(
        no_net_rotation_sequence,
        f' Reference frames: paleomag ({pmag_ref_frame_plate_id:03d}) and no-net rotation (000)')
elif default_reference_frame_plate_id == true_polar_wander_ref_frame_plate_id:
    # 000 is TPW (and 014 is PMAG).
    #
    # So 000->014 is the calculated TPW rotation sequence.
    rotation_time_samples_pmag = calculate_true_polar_wander_rotation_sequence(
        optimised_rotation_sequence,
        f' Reference frames: paleomag ({pmag_ref_frame_plate_id:03d}) and approx true polar wander (000)')
else:
    raise ValueError('Default reference frame plate ID is not one of the reference frame plate IDs')

# Create a new PMAG rotation sequence R(000->014).
pmag_rotation_feature = pygplates.Feature.create_total_reconstruction_sequence(
    fixed_plate_id=0,
    moving_plate_id=pmag_ref_frame_plate_id,
    total_reconstruction_pole=pygplates.GpmlIrregularSampling(rotation_time_samples_pmag))

# Default rotation model (anchor plate 000).
default_rotation_model = pygplates.RotationModel(output_rotation_feature_collections + [pmag_rotation_feature])


def remove_pmag_from_rotation_sequence(rotation_sequence, rotation_description):
    rotation_time_samples = []
    for rotation_sample in rotation_sequence:
        time = rotation_sample.get_time()
        finite_rotation = rotation_sample.get_value().get_finite_rotation()
        # R(000->014).
        finite_rotation_014_000 = default_rotation_model.get_rotation(time, pmag_ref_frame_plate_id)
        # pid->014 = pid->000 * 000->014
        # pid->000 = pid->014 * inverse[000->014]
        finite_rotation = finite_rotation * finite_rotation_014_000.get_inverse()
        rotation_time_samples.append(
            pygplates.GpmlTimeSample(
                    pygplates.GpmlFiniteRotation(finite_rotation),
                    time,
                    rotation_description,
                    rotation_sample.is_enabled()))
    return rotation_time_samples

# Generate the rotation sequences for the remaining reference frames (opt, NNR, TPW).
#
# This depends on which reference frame 000 is assigned to.
if default_reference_frame_plate_id == pmag_ref_frame_plate_id:
    # Both 000 and 014 are PMAG.
    #
    # 015->000 is the existing optimised rotation (000->005).
    rotation_time_samples_optimised = clone_rotation_sequence(
        optimised_rotation_sequence,
        f' Reference frames: paleomag (000) and optimised mantle ({optimised_ref_frame_plate_id:03d})')
    # 016->000 is the existing no-net rotation (000->005).
    rotation_time_samples_no_net_rotation = clone_rotation_sequence(
        no_net_rotation_sequence,
        f' Reference frames: paleomag (000) and no-net rotation ({no_net_rotation_ref_frame_plate_id:03d})')
    # 017->000 is the calculated TPW rotation sequence.
    rotation_time_samples_true_polar_wander = calculate_true_polar_wander_rotation_sequence(
        optimised_rotation_sequence,
        f' Reference frames: paleomag (000) and approx true polar wander ({true_polar_wander_ref_frame_plate_id:03d})')
elif default_reference_frame_plate_id == optimised_ref_frame_plate_id:
    # 000 is optimised (and 014 is PMAG).
    #
    # 015->000 is zero (identity).
    rotation_time_samples_optimised = create_identity_rotation_sequence(
        f' Reference frames: optimised mantle (000 and {optimised_ref_frame_plate_id:03d})',
        max_time_of_sequences_with_fixed_000)
    # 016->014 is the existing no-net rotation (000->005).
    # 016->014 = 016->000 * 000->014
    # 016->000 = 016->014 * inverse[000->014]
    rotation_time_samples_no_net_rotation = remove_pmag_from_rotation_sequence(
        no_net_rotation_sequence,
        f' Reference frames: optimised mantle (000) and no-net rotation ({no_net_rotation_ref_frame_plate_id:03d})')
    # 017->014 is the calculated TPW rotation sequence.
    # 017->014 = 017->000 * 000->014
    # 017->000 = 017->014 * inverse[000->014]
    rotation_time_samples_true_polar_wander = calculate_true_polar_wander_rotation_sequence(optimised_rotation_sequence, '')
    rotation_time_samples_true_polar_wander = remove_pmag_from_rotation_sequence(
        rotation_time_samples_true_polar_wander,
        f' Reference frames: optimised mantle (000) and approx true polar wander ({true_polar_wander_ref_frame_plate_id:03d})')
elif default_reference_frame_plate_id == no_net_rotation_ref_frame_plate_id:
    # 000 is no-net rotation (and 014 is PMAG).
    #
    # 015->014 is the existing optimised rotation (000->005).
    # 015->014 = 015->000 * 000->014
    # 015->000 = 015->014 * inverse[000->014]
    rotation_time_samples_optimised = remove_pmag_from_rotation_sequence(
        optimised_rotation_sequence,
        f' Reference frames: no-net rotation (000) and optimised mantle ({optimised_ref_frame_plate_id:03d})')
    # 016->000 is zero (identity).
    rotation_time_samples_no_net_rotation = create_identity_rotation_sequence(
        f' Reference frames: no-net rotation (000 and {no_net_rotation_ref_frame_plate_id:03d})',
        max_time_of_sequences_with_fixed_000)
    # 017->014 is the calculated TPW rotation sequence.
    # 017->014 = 017->000 * 000->014
    # 017->000 = 017->014 * inverse[000->014]
    rotation_time_samples_true_polar_wander = calculate_true_polar_wander_rotation_sequence(optimised_rotation_sequence, '')
    rotation_time_samples_true_polar_wander = remove_pmag_from_rotation_sequence(
        rotation_time_samples_true_polar_wander,
        f' Reference frames: no-net rotation (000) and approx true polar wander ({true_polar_wander_ref_frame_plate_id:03d})')
elif default_reference_frame_plate_id == true_polar_wander_ref_frame_plate_id:
    # 000 is TPW (and 014 is PMAG).
    #
    # 015->014 is the existing optimised rotation (000->005).
    # 015->014 = 015->000 * 000->014
    # 015->000 = 015->014 * inverse[000->014]
    rotation_time_samples_optimised = remove_pmag_from_rotation_sequence(
        optimised_rotation_sequence,
        f' Reference frames: approx true polar wander (000) and optimised mantle ({optimised_ref_frame_plate_id:03d})')
    # 016->014 is the existing no-net rotation (000->005).
    # 016->014 = 016->000 * 000->014
    # 016->000 = 016->014 * inverse[000->014]
    rotation_time_samples_no_net_rotation = remove_pmag_from_rotation_sequence(
        no_net_rotation_sequence,
        f' Reference frames: approx true polar wander (000) and no-net rotation ({no_net_rotation_ref_frame_plate_id:03d})')
    # 017->000 is zero (identity).
    rotation_time_samples_true_polar_wander = create_identity_rotation_sequence(
        f' Reference frames: approx true polar wander (000 and {true_polar_wander_ref_frame_plate_id:03d})',
        max_time_of_sequences_with_fixed_000)
else:
    raise ValueError('Default reference frame plate ID is not one of the reference frame plate IDs')

# Create a new optimised rotation sequence R(015->000).
optimised_rotation_feature = pygplates.Feature.create_total_reconstruction_sequence(
    fixed_plate_id=optimised_ref_frame_plate_id,
    moving_plate_id=0,
    total_reconstruction_pole=pygplates.GpmlIrregularSampling(rotation_time_samples_optimised))

# Create a new no-net rotation sequence R(016->000).
no_net_rotation_feature = pygplates.Feature.create_total_reconstruction_sequence(
    fixed_plate_id=no_net_rotation_ref_frame_plate_id,
    moving_plate_id=0,
    total_reconstruction_pole=pygplates.GpmlIrregularSampling(rotation_time_samples_no_net_rotation))

# Create a new true polar wander rotation sequence R(017->000).
true_polar_wander_rotation_feature = pygplates.Feature.create_total_reconstruction_sequence(
    fixed_plate_id=true_polar_wander_ref_frame_plate_id,
    moving_plate_id=0,
    total_reconstruction_pole=pygplates.GpmlIrregularSampling(rotation_time_samples_true_polar_wander))

# Add reference frames to a feature collection.
output_reference_frames_feature_collection = pygplates.FeatureCollection()
output_reference_frames_feature_collection.add(pmag_rotation_feature)
output_reference_frames_feature_collection.add(optimised_rotation_feature)
output_reference_frames_feature_collection.add(no_net_rotation_feature)
output_reference_frames_feature_collection.add(true_polar_wander_rotation_feature)

# Write the reference frames to their associated output rotation file.
output_reference_frames_feature_collection.write(output_reference_frames_filename)

# Write the output rotation filenames associated with the input rotation files.
for file_index in range(len(output_rotation_feature_collections)):
    output_rotation_feature_collection = output_rotation_feature_collections[file_index]
    output_rotation_filename = output_rotation_filenames[file_index]
    output_rotation_feature_collection.write(output_rotation_filename)
