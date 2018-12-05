import glob
import os.path

#
# Script to concatenate all rotation ('.rot') files in a directory.
#

# The main data directory is the directory containing this source file.
data_dir = os.path.abspath(os.path.dirname(__file__))

# Directory containing the input rotation files.
input_data_dir = os.path.join(data_dir, 'data', 'Global_Model_WD_Internal_Release_2016_v3')

# The concatenated output rotation file.
output_rotation_filename = os.path.join(input_data_dir, 'optimisation', 'all_rotations.rot')

# Gather all '.rot' files in the input directory.
input_rotation_filenames = glob.glob(os.path.join(input_data_dir, '*.rot'))

# Concatenate the input rotation files in the output rotation file.
with open(output_rotation_filename, 'w') as output_rotation_file:
    for input_rotation_filename in input_rotation_filenames:
        with open(input_rotation_filename) as input_rotation_file:
            output_rotation_file.write(input_rotation_file.read())
