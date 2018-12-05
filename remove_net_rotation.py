import sys
import os.path
import math
import pygplates

#
# Script to remove net rotation from a rotation file.
#

# The main directory is the directory containing this source file.
base_dir = os.path.abspath(os.path.dirname(__file__))

# Data sub-directory.
data_dir = os.path.join(base_dir, 'data')

rotation_filename = os.path.join(data_dir, 'Global_EarthByte_230-0Ma_GK07_AREPS.rot')
nnr_rotation_filename = os.path.join(data_dir, 'Global_EarthByte_230-0Ma_GK07_AREPS_NNR2.rot')

net_rotation_pole_data = [
    (230, -12.876, 161.064, 0.106018     ),
    (220, 6.44506, -99.2005, 0.210344    ),
    (210, 7.47756, -101.467, 0.263705    ),
    (200, 3.91095, -107.254, 0.31531     ),
    (190, 4.60476, -177.216, 0.304054    ),
    (180, -8.86242, 110.454, 0.0742011   ),
    (170, -13.751, 59.5681, 0.339836     ),
    (160, -14.3206, 110.991, 0.17593     ),
    (150, 12.6728, 133.79, 0.206295      ),
    (140, -29.6835, -22.2692, 0.257134   ),
    (130, -38.2763, -14.5914, 0.197295   ),
    (120, -79.9333, 44.9839, 0.210802    ),
    (110, -61.8724, 157.771, 0.205893    ),
    (100, -14.4632, 170.454, 0.197388    ),
    (90, -9.50488, 127.643, 0.192914     ),
    (80, -15.1639, 114.874, 0.202884     ),
    (70, -20.3742, 46.1005, 0.0629387    ),
    (60, 4.39409, 80.8206, 0.232849      ),
    (50, -33.3891, 110.667, 0.147159     ),
    (40, -14.8551, 80.3278, 0.0727454    ),
    (30, -44.1182, 50.3136, 0.142597     ),
    (20, -42.3818, 13.4825, 0.164551     ),
    (10, -36.4544, 30.9531, 0.195274     )
]



# Load the rotation features from rotation files.
rotation_features = list(pygplates.FeatureCollection(rotation_filename))

# A rotation model using the rotation features before they are modified.
rotation_model = pygplates.RotationModel(rotation_features)

total_reconstruction_pole_701 = None
for rotation_feature_index, rotation_feature in enumerate(rotation_features):

    # Get the rotation feature information.
    total_reconstruction_pole = rotation_feature.get_total_reconstruction_pole()
    if not total_reconstruction_pole:
        # Not a rotation feature.
        continue

    fixed_plate_id, moving_plate_id, rotation_sequence = total_reconstruction_pole
    # We're only interested in rotation features with moving plate ID 701.
    if moving_plate_id != 701:
        continue

    total_reconstruction_pole_701 = total_reconstruction_pole
    break

if not total_reconstruction_pole_701:
    print "Rotation 701 not found."
    sys.exit(1)


pole_time_samples_701_rel_fixed = []

for time, lat, lon, angle in reversed(net_rotation_pole_data):
    net_rotation = pygplates.FiniteRotation(pygplates.PointOnSphere(lat, lon), math.radians(angle))
    
    no_net_rotation_701_rel_001 = net_rotation.get_inverse() * rotation_model.get_rotation(time, 701, fixed_plate_id=1)
    
    no_net_rotation_701_rel_fixed = rotation_model.get_rotation(time, fixed_plate_id, fixed_plate_id=1).get_inverse() * no_net_rotation_701_rel_001

    pole_time_samples_701_rel_fixed.append(
        pygplates.GpmlTimeSample(pygplates.GpmlFiniteRotation(no_net_rotation_701_rel_fixed), time, 'NNR'))

# The time samples need to be wrapped into an irregular sampling property value.
total_reconstruction_pole_701_rel_fixed = pygplates.GpmlIrregularSampling(pole_time_samples_701_rel_fixed)

# Create the total reconstruction sequence (rotation) feature.
# rotation_feature_701_rel_fixed = pygplates.Feature(pygplates.FeatureType.gpml_total_reconstruction_sequence)
# rotation_feature_701_rel_fixed.set_name('INA-AUS Muller et.al 2000')
# rotation_feature_701_rel_fixed.set_total_reconstruction_pole(701, fixed_plate_id, total_reconstruction_pole_701_rel_fixed)
rotation_feature_701_rel_fixed = pygplates.Feature.create_total_reconstruction_sequence(
    701,
    fixed_plate_id,
    total_reconstruction_pole_701_rel_fixed,
    name='NNR')

# Replace the original 701 sequence with the new NNR 701 sequence.
rotation_features[rotation_feature_index] = rotation_feature_701_rel_fixed

# Write NNR rotation file.
pygplates.FeatureCollection(rotation_features).write(nnr_rotation_filename)
