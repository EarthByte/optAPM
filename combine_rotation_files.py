import glob
import os.path
import pandas
import pygplates
import Optimised_config

#########################################################################
# Script to concatenate all rotation ('.rot') files in a directory, and #
# reduce the set of rotations to only those plate circuit paths         #
# supporting the plate IDs required by the optimization process         #
# (eg, trenches and hotspots).                                          #
# Also changes fixed plate references to 000 to reference 005 instead   #
# and adds in a zero rotation for the plate pair 005/000 to complete    #
# the plate circuit to 000. The optimisation workflow will store        #
# optimised rotations in 005/000 so they affect all plates (not just    #
# those that go through Africa 701).                                    #
#########################################################################

# Currently only used for pre-resolved subduction zones.
times_start = Optimised_config.start_age

# Whether to fix Africa (701) to plate 000 (ie, zero rotations for 701 relative to 000).
# This is only used as a test input for the optimization workflow, to compare with (non-fixed)
# workflow runs that change various constraint weights (eg, trench migration vs net rotation).
fix_701_to_000 = False

# The main data directory is the directory containing this source file.
base_dir = os.path.abspath(os.path.dirname(__file__))

# Directory containing the input rotation files.
rotation_data_dir = os.path.join(base_dir, 'data', Optimised_config.data_model)

# The combined output rotation file.
output_rotation_filename = os.path.join(rotation_data_dir, 'optimisation', 'all_rotations.rot')

# When fixing 701 to 000 we also need a rotation file containing *all* rotations since no-net-rotation
# is calculated in GPlates by first loading all rotation files and topologies.
# Also we can't just adjust the 701 moving plate sequences because this does not shift any plates that
# do *not* pass through 701 in their plate circuit to 000 (such as 901 earlier than 80Ma).
# For the optimisation workflow we didn't need to worry about that because not using Pacific hotspots earlier than 80Ma,
# and net rotation is pre-calculated and trenches all pass through 701.
# But for calculating NNR we need proper topologies that don't break earlier than 80Ma.
if fix_701_to_000:
    output_rotation_filename_fixed_701_for_NNR = os.path.join(rotation_data_dir, 'optimisation', 'all_rotations_fixed_701_for_NNR.rot')

# Gather all '.rot' files in the input directory.
input_rotation_filenames = glob.glob(os.path.join(rotation_data_dir, '*.rot'))


# Main trench migration data directory.
tm_data_dir = os.path.join(base_dir, 'data', 'TMData', Optimised_config.data_model)

# Interpolated hotspots (used when 'interpolated_hotspot_trails' is True in main script).
interpolated_hotspots_filename = os.path.join(base_dir, 'data', Optimised_config.interpolated_hotspots)


#
# Get all plate IDs used by the optimisation process.
#

required_plate_ids = set()

# Add plate IDs are pre-resolved subduction zones.
for time in xrange(0, times_start+1):
    
    # print time
    
    # Add reference plate IDs (used by Net Rotation objective function).
    ref_rotation_plate_id, _ = Optimised_config.get_reference_params(time)
    required_plate_ids.add(ref_rotation_plate_id)
    
    # Read TM gpml file.
    subduction_features = pygplates.FeatureCollection.read(os.path.join(tm_data_dir, 'TMData_' + str(time) + 'Ma.gpml'))
    
    # Add subduction plate IDs to the set.
    for subduction_feature in subduction_features:
        required_plate_ids.add(subduction_feature.get_reconstruction_plate_id())

# Add plate IDs used for interpolated hotspots.
interpolated_hotspot_data = pandas.read_excel(interpolated_hotspots_filename)

for plate_id in interpolated_hotspot_data['PlateID']:
    required_plate_ids.add(int(plate_id))


print 'Required plate IDs:', sorted(required_plate_ids)
# print len(required_plate_ids)


# Read all the input rotation files.
rotation_features = []
for input_rotation_filename in input_rotation_filenames:
    rotation_features.extend(
        pygplates.FeatureCollection.read(input_rotation_filename))

plate_circuits = {}

# Create the plate circuits paths (dictionary of moving->fixed plate edges of plate circuit graph).
for rotation_feature_index, rotation_feature in enumerate(rotation_features):
    # Get the rotation feature information.
    total_reconstruction_pole = rotation_feature.get_total_reconstruction_pole()
    if not total_reconstruction_pole:
        # Not a rotation feature.
        continue

    fixed_plate_id, moving_plate_id, _ = total_reconstruction_pole
    
    # Can add both direction of edge in plate circuit graph.
    # But currently only adding the usual forward direction since rotation files currently assume that
    # and we're moving towards the anchor (0) plate which is at the root of rotation files.
    plate_circuits.setdefault(moving_plate_id, []).append((fixed_plate_id, rotation_feature_index))
    # plate_circuits.setdefault(fixed_plate_id, []).append((moving_plate_id, rotation_feature_index))


rotation_required = [False] * len(rotation_features)

visited_plate_ids = set()
def visit_plate_circuit(moving_plate_id):
    if moving_plate_id in visited_plate_ids:
        return
    visited_plate_ids.add(moving_plate_id)
    
    moving_plate_circuits = plate_circuits.get(moving_plate_id)
    if moving_plate_circuits:
        for fixed_plate_id, rotation_feature_index in moving_plate_circuits:
            rotation_required[rotation_feature_index] = True
            # TODO: Should be terminate (ie, not visit) when fixed_plate_id==0 ?
            visit_plate_circuit(fixed_plate_id)

# Find those fixed/moving rotation pairs required to support plate circuit from require plate IDs to anchor 000.
for required_plate_id in required_plate_ids:
    visit_plate_circuit(required_plate_id)

required_rotation_features = []
for rotation_feature_index, rotation_feature in enumerate(rotation_features):
    if rotation_required[rotation_feature_index]:
        required_rotation_features.append(rotation_feature)


# Change fixed plate IDs from 000 to 005.
# But only if necessary (it's possible the input rotations came from a previous optimised run and already reference 005).
required_rotations_include_005_000 = False
for rotation_feature in required_rotation_features:
    # Get the rotation feature information.
    total_reconstruction_pole = rotation_feature.get_total_reconstruction_pole()
    if not total_reconstruction_pole:
        # Not a rotation feature.
        continue

    fixed_plate_id, moving_plate_id, rotation_sequence = total_reconstruction_pole
    
    # Change fixed plate ID 000 to 005 (unless it's the 005-000 plate pair).
    # We want everything that references 000 to now reference 005.
    # Later we'll add a 005-000 sequence to store optimised rotation adjustments (if there isn't already one).
    if fixed_plate_id == 0:
        if moving_plate_id == 5:
            required_rotations_include_005_000 = True
        else:
            rotation_feature.set_total_reconstruction_pole(5, moving_plate_id, rotation_sequence)

# If the rotation features don't already contain a 005-000 plate pair then add one
# so that we have a path to 000 since we changed 000 to 005 above.
# We'll add a zero rotation so that we don't modify the absolute rotations (of the various plate IDs)
# because this file will be used as input when subtracting out net rotations to produce a NNR rotation file.
if not required_rotations_include_005_000:
    rotation_feature_005_rel_000 = pygplates.Feature.create_total_reconstruction_sequence(
        0,
        5,
        pygplates.GpmlIrregularSampling([
            pygplates.GpmlTimeSample(pygplates.GpmlFiniteRotation(pygplates.FiniteRotation()), time, 'Absolute Reference Frame Optimisation')
                for time in (0.0, times_start)]))
    required_rotation_features.append(rotation_feature_005_rel_000)

pygplates.FeatureCollection(required_rotation_features).write(output_rotation_filename)


# Fix Africa (701) to plate 000 if requested (for normal optimization runs we don't do this).
if fix_701_to_000:
    
    unfixed_required_rotation_model = pygplates.RotationModel(required_rotation_features)
    for rotation_feature in required_rotation_features:
        total_reconstruction_pole = rotation_feature.get_total_reconstruction_pole()
        if total_reconstruction_pole:
            fixed_plate_id, moving_plate_id, rotation_sequence = total_reconstruction_pole
            if moving_plate_id == 701:
                rotation_samples = rotation_sequence.get_time_samples()
                for rotation_sample in rotation_samples:
                    # We want our 701 to 000 rotation to be zero.
                    # However the 701 sequence might have a fixed plate ID that is not 000.
                    # So convert zero 701 rel 000 rotation to the 701 rel 'fixed_plate_id' rotation to store in rotation feature.
                    #
                    #                         R(0->t,000->701) = R(0->t,000->fixed) * R(0->t,fixed->701)
                    #                                 Identity = R(0->t,000->fixed) * R(0->t,fixed->701)
                    #   inverse(R(0->t,000->fixed)) * Identity = R(0->t,fixed->701)
                    #                       R(0->t,fixed->701) = inverse(R(0->t,000->fixed))
                    #
                    zero_rotation_701_rel_conjugate = unfixed_required_rotation_model.get_rotation(
                        rotation_sample.get_time(),
                        fixed_plate_id,
                        fixed_plate_id=0).get_inverse()
                    rotation_sample.get_value().set_finite_rotation(zero_rotation_701_rel_conjugate)
    
    #
    # Create the special no-net-rotation file containing *all* rotations (not just required rotations)
    # that are adjusted such that Africa (701) is fixed to 000.
    #
    
    # Read all the input rotation files again.
    nnr_rotation_features = []
    for input_rotation_filename in input_rotation_filenames:
        nnr_rotation_features.extend(
            pygplates.FeatureCollection.read(input_rotation_filename))
    
    # Change fixed plate IDs from 000 to 005.
    # But only if necessary (it's possible the input rotations came from a previous optimised run and already reference 005).
    nnr_rotation_features_tmp = []
    for rotation_feature in nnr_rotation_features:
        total_reconstruction_pole = rotation_feature.get_total_reconstruction_pole()
        if total_reconstruction_pole:
            fixed_plate_id, moving_plate_id, rotation_sequence = total_reconstruction_pole
            # Change fixed plate ID 000 to 005 (unless it's the 005-000 plate pair).
            # We want everything that references 000 to now reference 005.
            # Later we'll add a 005-000 sequence to store rotation adjustments.
            if fixed_plate_id == 0 and  moving_plate_id != 5:
                rotation_feature.set_total_reconstruction_pole(5, moving_plate_id, rotation_sequence)
            if not (fixed_plate_id == 0 and moving_plate_id == 5):
                # Add all rotation features except 005-000. We'll add our own 005-000 later.
                nnr_rotation_features_tmp.append(rotation_feature)
    
    # NNR rotation features now exclude moving plate 005, but we'll add below.
    nnr_rotation_features = nnr_rotation_features_tmp
    nnr_rotation_model = pygplates.RotationModel(nnr_rotation_features)
    
    zero_rotation_time_samples_005_rel_000 = []
    # Start with identity rotation at time 0Ma.
    zero_rotation_time_samples_005_rel_000.append(
        pygplates.GpmlTimeSample(pygplates.GpmlFiniteRotation(pygplates.FiniteRotation()), 0.0, '701 fixed to 000 for NNR calculation'))
    for time in xrange(1, times_start+1):
        
        # We want our 701 to 000 rotation to be zero.
        #
        #                         R(0->t,000->701) = R(0->t,000->005) * R(0->t,005->701)
        #                                 Identity = R(0->t,000->005) * R(0->t,005->701)
        #   inverse(R(0->t,000->005)) * Identity = R(0->t,005->701)
        #                       R(0->t,000->005) = inverse(R(0->t,005->701))
        #
        # NOTE: We're setting *anchor* plate to 5 (instead of *fixed* plate) since we don't yet have a path to 000.
        zero_rotation_005_rel_000 = nnr_rotation_model.get_rotation(time, 701, anchor_plate_id=5).get_inverse()
        
        zero_rotation_time_samples_005_rel_000.append(
            pygplates.GpmlTimeSample(pygplates.GpmlFiniteRotation(zero_rotation_005_rel_000), time, '701 fixed to 000 for NNR calculation'))
    
    # Create a new 005/000 rotation sequence.
    rotation_feature_005_000 = pygplates.Feature.create_total_reconstruction_sequence(
        0,
        5,
        pygplates.GpmlIrregularSampling(zero_rotation_time_samples_005_rel_000))
    
    nnr_rotation_features.append(rotation_feature_005_000)
    
    pygplates.FeatureCollection(nnr_rotation_features).write(output_rotation_filename_fixed_701_for_NNR)
