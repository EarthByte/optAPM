import glob
import os.path
import sys
import pygplates as pgp
import subduction_convergence_for_absolute_plate_motion as scap


# Check the required pygplates version.
# Require a version that closes the gaps in resolved topologies in the *deforming* model (along deforming lines).
PYGPLATES_VERSION_REQUIRED = pgp.Version(19)
# Check the imported pygplates version.
if not hasattr(pgp, 'Version') or pgp.Version.get_imported_version() < PYGPLATES_VERSION_REQUIRED:
    raise RuntimeError('{0}: Error - imported pygplates version {1} but version {2} or greater is required'.format(
            os.path.basename(__file__), pgp.Version.get_imported_version(), PYGPLATES_VERSION_REQUIRED))


####################################################################################
# Generates trench migration data required for optimisation routine.               #
#                                                                                  #
# The subducting trenches are reconstructed from 0-250Ma (in 1My intervals)        #
# and stored in files ready for quick retrieval during the optimization stage.     #
# This is done as a pre-process because resolving the boundaries of                #
# continuously-closing plates and deforming networks is relatively time consuming. #
####################################################################################

timeStart = 410

# If using optimised rotation model to reconstruct/resolve topologies.
# model = '1'


# The main data directory is the 'data' sub-directory of the directory containing this source file.
data_dir = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'data')

# Main trench migration data directory.
tm_data_dir = os.path.join(data_dir, 'TMData')


# Muller 2016
# starting_rotfile = '/Users/Mike/PhD/Plate_Models/Muller_etal_2016/Global_EarthByte_230-0Ma_GK07_AREPS.rot'
# output_data_dir = '/Users/Mike/PhD/Plate_Models/Muller_etal_2016/Resolved_subduction_zones/'
# topologies = '/Users/Mike/PhD/Plate_Models/Muller_etal_2016/Global_EarthByte_230-0Ma_GK07_AREPS_PlateBoundaries.gpml'

# Shephard 2016
# starting_rotfile = '/Users/Mike/PhD/Plate_Models/Shephard_etal_2013_ESR/GPlates/Shephard_etal_ESR2013_Global_EarthByte_2013.rot'
# output_data_dir = '/Users/Mike/PhD/Plate_Models/Shephard_etal_2013_ESR/GPlates/Resolved_subduction_zones/'
# topologies = '/Users/Mike/PhD/Plate_Models/Shephard_etal_2013_ESR/GPlates/Shephard_etal_ESR2013_Global_EarthByte_2013.gpml'

# Muller 2016
# starting_rotfiles = [os.path.join(data_dir, 'Global_EarthByte_230-0Ma_GK07_AREPS_optAPM' + model + '.rot')]

# Shephard 2013    
# starting_rotfile = '/Users/Mike/Projects/optAPM/model_output_0-80Ma/optAPM' + model + \
#                    '/Shephard_etal_ESR2013_Global_EarthByte_2013_optAPM' + model + '.rot'

# Muller 2016 alternative reference frames
# starting_rotfile = '/Users/Mike/PhD/Plate_Models/Muller_etal_2016/Alternative_reference_frames/' + \
#                    'Global_EarthByte_230-0Ma_GK07_AREPS_' + model + '.rot'

# output_data_dir = '/Users/Mike/Projects/optAPM/model_output_0-80Ma/optAPM' + model + \
#                   '/Subduction_zone_params/Resolved_topologies/'
# output_data_dir = '/Users/Mike/PhD/Plate_Models/Muller_etal_2016/Alternative_reference_frames/Resolved_topologies/' + model + '/'
# output_data_dir = os.path.join(data_dir, 'TMData', 'Muller_2016')

# topology_files = [
#     os.path.join(data_dir, 'Global_EarthByte_230-0Ma_GK07_AREPS_PlateBoundaries.gpml'),
#     os.path.join(data_dir, 'Global_EarthByte_230-0Ma_GK07_AREPS_Topology_BuildingBlocks.gpml')]
#topologies = '/Users/Mike/PhD/Plate_Models/Shephard_etal_2013_ESR/GPlates/Shephard_etal_ESR2013_Global_EarthByte_2013.gpml'

#
# Global_Model_WD_Internal_Release_2016_v3
#
output_data_dir = os.path.join(tm_data_dir, 'Global_Model_WD_Internal_Release_2016_v3')
input_data_dir = os.path.join(data_dir, 'Global_Model_WD_Internal_Release_2016_v3')
topology_files = glob.glob(os.path.join(input_data_dir, '*.gpml'))
# Gather all '.rot' files in the input directory.
starting_rotfiles = glob.glob(os.path.join(input_data_dir, '*.rot'))
# ...or use the optimised rotation file from the most recent model run.
# starting_rotfiles = os.path.join(data_dir, 'Global_Model_WD_Internal_Release_2016_v3', 'optimisation', 'all_rotations_optAPM' + model + '.rot')


# Load/parse the rotation and topologies files once (instead of repeating in each time loop iteration).
topologies = [pgp.FeatureCollection(topology_file) for topology_file in topology_files]
starting_rotation_model = pgp.RotationModel(starting_rotfiles)

for time in xrange(0, timeStart+1):
    
    print time
    
    # Resolve subduction zone features
    resolved_subduction_zone_features = scap.resolve_subduction_zones(starting_rotation_model, topologies, time)

    # Write to gpml file
    pgp.FeatureCollection(resolved_subduction_zone_features).write(os.path.join(output_data_dir, 'TMData_' + str(time) + 'Ma.gpml'))
