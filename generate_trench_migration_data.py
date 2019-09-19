import glob
import os.path
import sys
import warnings
import pygplates as pgp
import subduction_convergence_for_absolute_plate_motion as scap
import Optimised_config


# Check the required pygplates version.
# PyGPlates version 19 can close the gaps in resolved topologies in the *deforming* model (along deforming lines).
# PyGPlates version 22 can handle topological lines (can get their sub-sub-segment plate IDs).
PYGPLATES_VERSION_REQUIRED = pgp.Version(22)
# Check the imported pygplates version.
if not hasattr(pgp, 'Version') or pgp.Version.get_imported_version() < PYGPLATES_VERSION_REQUIRED:
    raise RuntimeError('{0}: Error - imported pygplates version {1} but version {2} or greater is required'.format(
            os.path.basename(__file__), pgp.Version.get_imported_version(), PYGPLATES_VERSION_REQUIRED))


def warning_format(message, category, filename, lineno, file=None, line=None):
    # return '{0}:{1}: {1}:{1}\n'.format(filename, lineno, category.__name__, message)
    return '{0}: {1}\n'.format(category.__name__, message)
# Print the warnings without the filename and line number.
# Users are not going to want to see that.
warnings.formatwarning = warning_format
# Always print warnings (not just the first time encountered at a particular location).
warnings.simplefilter("always")


####################################################################################
# Generates trench migration data required for optimisation routine.               #
#                                                                                  #
# The subducting trenches are reconstructed from 0-250Ma (in 1My intervals)        #
# and stored in files ready for quick retrieval during the optimization stage.     #
# This is done as a pre-process because resolving the boundaries of                #
# continuously-closing plates and deforming networks is relatively time consuming. #
####################################################################################

timeStart = Optimised_config.start_age

# If using optimised rotation model to reconstruct/resolve topologies.
# model = '1'


# The main data directory is the 'data' sub-directory of the directory containing this source file.
data_dir = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'data')

# Main trench migration data directory.
tm_data_dir = os.path.join(data_dir, 'TMData')


output_data_dir = os.path.join(tm_data_dir, Optimised_config.data_model)
input_data_dir = os.path.join(data_dir, Optimised_config.data_model)
topology_files = glob.glob(os.path.join(input_data_dir, '*.gpml'))
starting_rotfiles = glob.glob(os.path.join(input_data_dir, '*.rot'))


# Load/parse the rotation and topologies files once (instead of repeating in each time loop iteration).
topologies = [pgp.FeatureCollection(topology_file) for topology_file in topology_files]
starting_rotation_model = pgp.RotationModel(starting_rotfiles)

for time in xrange(0, timeStart+1):
    
    print time
    
    # Resolve subduction zone features
    resolved_subduction_zone_features = scap.resolve_subduction_zones(starting_rotation_model, topologies, time)

    # Write to gpml file
    pgp.FeatureCollection(resolved_subduction_zone_features).write(os.path.join(output_data_dir, 'TMData_' + str(time) + 'Ma.gpml'))
