from __future__ import print_function
import sys
import os.path
import math
import pygplates
import Optimised_APM  # To query reference plate ID over time.


# # Check the required pygplates version.
# # Need the bug fix for crash writing out rotation features to '.rot' files.
# # UPDATE: Workaround is to just not write out a feature name or description.
#
# PYGPLATES_VERSION_REQUIRED = pygplates.Version(20)
# # Check the imported pygplates version.
# if not hasattr(pygplates, 'Version') or pygplates.Version.get_imported_version() < PYGPLATES_VERSION_REQUIRED:
#     raise RuntimeError('{0}: Error - imported pygplates version {1} but version {2} or greater is required'.format(
#             os.path.basename(__file__), pygplates.Version.get_imported_version(), PYGPLATES_VERSION_REQUIRED))


##########################################################
# Script to remove net rotation from a rotation file.    #
#                                                        #
# Note that this removes net rotation from the reference #
# plate (which can change over time). Typically this is  #
# Africa 701. Normally this has the undesired effect of  #
# bypassing plate circuits that don't go through the     #
# reference plate, but this is OK for the optimisation   #
# workflow because, for the net rotation part, it only   #
# looks at the reference plate (at a particular time).   #
##########################################################


# The main directory is the directory containing this source file.
base_dir = os.path.abspath(os.path.dirname(__file__))

# Data sub-directory.
data_dir = os.path.join(base_dir, 'data')

rotation_filename = os.path.join(data_dir, 'Global_1000-0_Model_2017', 'optimisation', 'all_rotations.rot')
nnr_rotation_filename = os.path.join(data_dir, 'Global_1000-0_Model_2017', 'optimisation', 'no_net_rotations.rot')

net_rotations_filename = os.path.join(data_dir, 'Global_1000-0_Model_2017', 'optimisation', 'total-net-rotations.csv')



# The rotation model before net rotation is removed.
rotation_model = pygplates.RotationModel(rotation_filename)

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


# Each item is a 2-tuple (ref_rotation_plate_id, list of time samples).
nnr_ref_plate_sequences = []

# Start with identity rotation at time 0Ma for the first reference plate.
first_ref_rotation_plate_id, _ = Optimised_APM.get_reference_params(0)
nnr_ref_plate_sequences.append((
    first_ref_rotation_plate_id,
    [pygplates.GpmlTimeSample(pygplates.GpmlFiniteRotation(pygplates.FiniteRotation()), 0.0, 'NNR')]))

# Start with identity net rotation at time 0Ma.
net_total_rotation = pygplates.FiniteRotation()

for time, net_stage_lat, net_stage_lon, net_stage_angle_per_my in reversed(net_stage_pole_data):
    
    # Accumulate net stage rotations going backward in time (hence the inverse stage rotation)
    # since finite rotations go backward in time in the rotation file.
    net_stage_rotation = pygplates.FiniteRotation(
            pygplates.PointOnSphere(net_stage_lat, net_stage_lon),
            math.radians(net_stage_angle_per_my * net_stage_rotation_interval))
    net_total_rotation = net_stage_rotation.get_inverse() * net_total_rotation
    
    # The reference plate for the current time obtained from main optimisation workflow.
    ref_rotation_plate_id, _ = Optimised_APM.get_reference_params(time)
    
    # Remove net total rotation at current time from 'ref_rotation_plate_id' rel 000 rotation.
    #
    #   R_net_rotation * R_no_net_rotation(0->t,000->ref_plate) = R(0->t,000->ref_plate)
    #                    R_no_net_rotation(0->t,000->ref_plate) = inverse(R_net_rotation) * R(0->t,000->ref_plate)
    #
    no_net_rotation = net_total_rotation.get_inverse() * rotation_model.get_rotation(time, ref_rotation_plate_id)
    
    # If the reference plate changes then start a new rotation sequence.
    last_ref_rotation_plate_id = nnr_ref_plate_sequences[-1][0]
    if ref_rotation_plate_id != last_ref_rotation_plate_id:
        nnr_ref_plate_sequences.append((ref_rotation_plate_id, []))
    
    # Add current no-net-rotation to the current reference plate sequence.
    nnr_ref_plate_sequences[-1][1].append(
        pygplates.GpmlTimeSample(pygplates.GpmlFiniteRotation(no_net_rotation), time, 'NNR'))

# There is one feature for each contiguous reference plate sequence.
# This is all that's needed for the main optimisation workflow to compare rotations against no-net-rotation.
#
# For example, if there is only reference plate and it's 701 then there will be only one
# 701 no-net-rotation feature even though 701 might have been spread across 2 rotation files
# (as two separate rotation features).
# Another example, if the reference plate switches once from 701 to 101 then there'll be 2 NNR features.
nnr_rotation_features = []
for ref_rotation_plate_id, nnr_pole_time_samples in nnr_ref_plate_sequences:
    # Create the no-net-rotation total reconstruction sequence (rotation) feature.
    nnr_rotation_feature = pygplates.Feature.create_total_reconstruction_sequence(
        0,
        ref_rotation_plate_id,
        # The time samples need to be wrapped into an irregular sampling property value...
        pygplates.GpmlIrregularSampling(nnr_pole_time_samples))
    
    nnr_rotation_features.append(nnr_rotation_feature)

# Write NNR rotation file.
pygplates.FeatureCollection(nnr_rotation_features).write(nnr_rotation_filename)
