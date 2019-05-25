from __future__ import print_function
import sys
import os.path
import math
import pygplates


# # Check the required pygplates version.
# # Need the bug fix for crash writing out rotation features to '.rot' files.
# # UPDATE: Workaround is to just not write out a feature name or description.
#
# PYGPLATES_VERSION_REQUIRED = pygplates.Version(20)
# # Check the imported pygplates version.
# if not hasattr(pygplates, 'Version') or pygplates.Version.get_imported_version() < PYGPLATES_VERSION_REQUIRED:
#     raise RuntimeError('{0}: Error - imported pygplates version {1} but version {2} or greater is required'.format(
#             os.path.basename(__file__), pygplates.Version.get_imported_version(), PYGPLATES_VERSION_REQUIRED))


#########################################################
# Script to remove net rotation from a rotation file.   #
#                                                       #
# Note that this removes net rotation from Africa 701.  #
# This has the undesired effect of bypassing plate      #
# circuits that don't go through Africa, but this is OK #
# for the optimisation workflow because, for the net    #
# rotation part, it only looks at Africa.               #
#########################################################


# The main directory is the directory containing this source file.
base_dir = os.path.abspath(os.path.dirname(__file__))

# Data sub-directory.
data_dir = os.path.join(base_dir, 'data')

rotation_filename = os.path.join(data_dir, 'Global_Model_WD_Internal_Release_2016_v3', 'optimisation', 'all_rotations.rot')
nnr_rotation_filename = os.path.join(data_dir, 'Global_Model_WD_Internal_Release_2016_v3', 'optimisation', 'all_rotations_NNR.rot')

net_rotations_filename = os.path.join(data_dir, 'Global_Model_WD_Internal_Release_2016_v3', 'optimisation', 'total-net-rotations.csv')


def read_net_rotations(net_rotations_filename):
    
    net_rotations = []
    with open(net_rotations_filename, 'r') as net_rotations_file:
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
                net_rotation = [float(line_string) for line_string in line_string_list]
            except ValueError:
                continue

            net_rotations.append(net_rotation)
    
    return net_rotations


# Read the total net rotation file exported by GPlates.
# NOTE: Must be exported by a version of GPlates that supports deforming networks in Net Rotation export.
net_stage_pole_data = read_net_rotations(net_rotations_filename)

if len(net_stage_pole_data) < 2:
    print("Need at least two net rotations in {0}.".format(net_rotations_filename), file=sys.stderr)
    sys.exit(1)

# Interval is time difference between first two net rotation times (in first column of first two rows).
net_stage_rotation_interval = net_stage_pole_data[0][0] - net_stage_pole_data[1][0]

# print('Net rotation interval:', net_stage_rotation_interval)
# print('Net rotations:', net_stage_pole_data)


# Load the rotation features from rotation files.
rotation_features = list(pygplates.FeatureCollection(rotation_filename))

# A rotation model using the rotation features before they are modified.
rotation_model = pygplates.RotationModel(rotation_features)


# Remove any moving plate 701 rotation features (we'll add an NNR 701 rotation feature later).
rotation_features_tmp = []
for rotation_feature in rotation_features:

    # Get the rotation feature information.
    total_reconstruction_pole = rotation_feature.get_total_reconstruction_pole()
    if not total_reconstruction_pole:
        # Not a rotation feature.
        continue

    _, moving_plate_id, _ = total_reconstruction_pole
    if moving_plate_id != 701:
        # Add all rotation features except moving plate 701.
        rotation_features_tmp.append(rotation_feature)

# Rotation features now exclude moving plate 701, but we'll add below.
rotation_features = rotation_features_tmp


pole_time_samples_701_rel_001 = []

# Start with identity rotation at time 0Ma.
pole_time_samples_701_rel_001.append(
    pygplates.GpmlTimeSample(pygplates.GpmlFiniteRotation(pygplates.FiniteRotation()), 0.0, 'NNR'))

# Start with identity net rotation at time 0Ma.
net_total_rotation = pygplates.FiniteRotation()

for time, net_stage_lat, net_stage_lon, net_stage_angle_per_my in reversed(net_stage_pole_data):
    
    # Accumulate net stage rotations going backward in time (hence the inverse stage rotation)
    # since finite rotations go backward in time in the rotation file.
    net_stage_rotation = pygplates.FiniteRotation(
            pygplates.PointOnSphere(net_stage_lat, net_stage_lon),
            math.radians(net_stage_angle_per_my * net_stage_rotation_interval))
    net_total_rotation = net_stage_rotation.get_inverse() * net_total_rotation
    
    # Remove net total rotation at current time from 701 rel 001 rotation.
    no_net_rotation_701_rel_001 = net_total_rotation.get_inverse() * rotation_model.get_rotation(time, 701, fixed_plate_id=1)
    
    pole_time_samples_701_rel_001.append(
        pygplates.GpmlTimeSample(pygplates.GpmlFiniteRotation(no_net_rotation_701_rel_001), time, 'NNR'))

# The time samples need to be wrapped into an irregular sampling property value.
total_reconstruction_pole_701_rel_001 = pygplates.GpmlIrregularSampling(pole_time_samples_701_rel_001)

# Create the total reconstruction sequence (rotation) feature.
rotation_feature_701_rel_001 = pygplates.Feature.create_total_reconstruction_sequence(
    1,
    701,
    total_reconstruction_pole_701_rel_001)

# Add the new NNR 701 sequence.
rotation_features.append(rotation_feature_701_rel_001)

# Write NNR rotation file.
pygplates.FeatureCollection(rotation_features).write(nnr_rotation_filename)
