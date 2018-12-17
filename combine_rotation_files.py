import glob
import os.path
import pandas
import pygplates

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

# Optimising for 0-410Ma.
# Currently only used for pre-resolved subduction zones.
times_start = 410

# The main data directory is the directory containing this source file.
base_dir = os.path.abspath(os.path.dirname(__file__))

# Directory containing the input rotation files.
rotation_data_dir = os.path.join(base_dir, 'data', 'Global_Model_WD_Internal_Release_2016_v3')

# The combined output rotation file.
output_rotation_filename = os.path.join(rotation_data_dir, 'optimisation', 'all_rotations.rot')

# Gather all '.rot' files in the input directory.
input_rotation_filenames = glob.glob(os.path.join(rotation_data_dir, '*.rot'))


# Main trench migration data directory.
tm_data_dir = os.path.join(base_dir, 'data', 'TMData', 'Global_Model_WD_Internal_Release_2016_v3')

# Interpolated hotspots (used when 'interpolated_hotspot_trails' is True in main script).
interpolated_hotspots_filename = os.path.join(base_dir, 'data', 'interpolated_hotspot_chains_5Myr.xlsx')


#
# Get all plate IDs used by the optimisation process.
#

required_plate_ids = set()

# Add Africa (used by Net Rotation objective function).
required_plate_ids.add(701)

# Add plate IDs are pre-resolved subduction zones.
for time in xrange(0, times_start+1):
    
    # print time
    
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


# Read all the rotation files.
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
                for time in (0.0, 600.0)]))
    required_rotation_features.append(rotation_feature_005_rel_000)

pygplates.FeatureCollection(required_rotation_features).write(output_rotation_filename)