import glob
import os.path


#########################################
# Optimisation configuration parameters #
#########################################


##########################################################################################
# Supported parallelisation methods (or None to disable parallelisation, eg, for testing).
#
MPI4PY = 0
IPYPARALLEL = 1

# Choose parallelisation method (or None to disable parallelisation, eg, for testing).
#use_parallel = MPI4PY  # For example, to use with 'mpiexec -n <cores> python Optimised_APM.py'.
use_parallel = None
##########################################################################################
# The root input data directory ('data/').
# This is the 'data/' sub-directory of the directory containing this source file.
datadir = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'data', '')


# The model name is suffixed to various output filenames.
model_name = "run_test"

start_age = 410
end_age = 0
interval = 10

models = 100


# The data model to run the optimisation on.
# This should be the name of the sub-directory in 'data/' and also the sub-directory in 'data/TMData/'.
data_model = 'Global_Model_WD_Internal_Release_2019_v2'

# The original rotation files (relative to the 'data/<data_model>/' directory).
#
# Can either:
#   1) use glob to automatically find all the '.rot' files (you don't need to do anything), or
#   2) explicitly list all the rotation files (you need to list the filenames).
#
# 1) Automatically gather all '.rot' files (and make filenames relative to the 'data/' directory).
original_rotation_filenames = [os.path.relpath(abs_path, datadir) for abs_path in
        glob.glob(os.path.join(datadir, data_model, '*.rot'))]
# 2) Explicitly list all the input rotation files (must be relative to the 'data/' directory).
#original_rotation_filenames = [
#  'Global_Model_WD_Internal_Release_2019_v2/rotation_file1.rot',
#  'Global_Model_WD_Internal_Release_2019_v2/rotation_file2.rot',
#]

# The topology files (relative to the 'data/<data_model>/' directory).
#
# Can either:
#   1) use glob to automatically find all the '.gpml' files (you don't need to do anything), or
#   2) explicitly list all the topology files (you need to list the filenames).
#
# 1) Automatically gather all '.gpml' and '.gpmlz' files (and make filenames relative to the 'data/' directory).
topology_filenames = [os.path.relpath(abs_path, datadir) for abs_path in
        glob.glob(os.path.join(datadir, data_model, '*.gpml')) + glob.glob(os.path.join(datadir, data_model, '*.gpmlz'))]
# 2) Explicitly list all the topology files (must be relative to the 'data/' directory).
#topology_filenames = [
#  'Global_Model_WD_Internal_Release_2019_v2/topology_file1.gpml',
#  'Global_Model_WD_Internal_Release_2019_v2/topology_file2.gpml',
#]

ridge_file = data_model + '/StaticGeometries/AgeGridInput/Global_EarthByte_GeeK07_Ridges_2019_v1.gpml'
isochron_file = data_model + '/StaticGeometries/AgeGridInput/Global_EarthByte_GeeK07_Isochrons_2019_v1.gpml'
isocob_file = data_model + '/StaticGeometries/AgeGridInput/Global_EarthByte_GeeK07_IsoCOB_2019_v1.gpml'
#
# For (data_model == 'Global_1000-0_Model_2017') or (data_model == 'Muller++_2015_AREPS_CORRECTED') ...
#
##################################################################################################################################
#
# There are no static geometries (besides coastlines) for this data model.
#
# NOTE: SO USING SAME FILES AS 'Global_Model_WD_Internal_Release_2019_v2'.
#       THIS IS OK IF WE'RE NOT INCLUDING FRACTURE ZONES (BECAUSE THEN THESE FILES ARE NOT USED FOR FINAL OPTIMISED ROTATIONS).
#
##################################################################################################################################
# ridge_file = 'Global_Model_WD_Internal_Release_2019_v2/StaticGeometries/AgeGridInput/Global_EarthByte_GeeK07_Ridges_2019_v1.gpml'
# isochron_file = 'Global_Model_WD_Internal_Release_2019_v2/StaticGeometries/AgeGridInput/Global_EarthByte_GeeK07_Isochrons_2019_v1.gpml'
# isocob_file = 'Global_Model_WD_Internal_Release_2019_v2/StaticGeometries/AgeGridInput/Global_EarthByte_GeeK07_IsoCOB_2019_v1.gpml'
#
# Original files used in original optimisation script...
#
# ridge_file = 'Global_EarthByte_230-0Ma_GK07_AREPS_Ridges.gpml'
# isochron_file = 'Global_EarthByte_230-0Ma_GK07_AREPS_Isochrons.gpmlz'
# isocob_file = 'Global_EarthByte_230-0Ma_GK07_AREPS_IsoCOB.gpml'


#
# Which components are enabled and their weightings.
#
# NOTE: The weights are inverse weights (ie, the constraint costs are *multiplied* by "1.0 / weight").
#
def get_fracture_zone_params(age):
    return False, 1.0  # Disable fracture zones.

def get_net_rotation_params(age):
    if age <= 80:
        return True, 1.0
    elif age <= 170:
        # NOTE: These are inverse weights (ie, the constraint costs are *multiplied* by "1.0 / weight").
        return True, 2.0 # Gives a *multiplicative* weight of 0.5
    else:
        # NOTE: These are inverse weights (ie, the constraint costs are *multiplied* by "1.0 / weight").
        return True, 5.0 # Gives a *multiplicative* weight of 0.2

def get_trench_migration_params(age):
    return True, 1.0

def get_hotspot_trail_params(age):
    # Only use hotspot trails for 0-80Ma.
    if age <= 80:
        return True, 1.0
    else:
        return False, 1.0

def get_plate_velocity_params(age):
    return True, 1.0

# The grid spacing (in degrees) between points in the grid used for plate velocity calculations.
plate_velocity_grid_spacing = 2.0


#
# Which reference plate ID and PMAG rotation file to use at which age.
#
def get_reference_params(age):
    if data_model == 'Global_1000-0_Model_2017':
        if age <= 550:
            ref_rotation_plate_id = 701
            pmag_rotfile = 'Global_1000-0_Model_2017/pmag/550_0_Palaeomagnetic_Africa_S.rot'
        else:
            ref_rotation_plate_id = 101
            pmag_rotfile = 'Global_1000-0_Model_2017/pmag/1000_550_Laurentia_pmag_reference.rot'
    else:
        ref_rotation_plate_id = 701
        pmag_rotfile = 'Palaeomagnetic_Africa_S.rot'
    
    return ref_rotation_plate_id, pmag_rotfile


search = "Initial"
search_radius = 60
# If True then temporarily expand search radius to 90 whenever the reference plate changes.
# Normally the reference plate stays constant at Africa (701), but does switch to 101 for the 1Ga model.
# It's off by default since it doesn't appear to change the results, and may sometimes cause job to fail
# on Artemis (presumably since 'models' is increased by a factor of 2.5) - although problem manifested
# as failure to read the rotation file being optimised, so it was probably something else.
expand_search_radius_on_ref_plate_switches = False
rotation_uncertainty = 30
auto_calc_ref_pole = True

model_stop_condition = 'threshold'
max_iter = 5  # Only applies if model_stop_condition != 'threshold'


# Trench migration parameters
tm_method = 'pygplates' # 'pygplates' for new method OR 'convergence' for old method
tm_data_type = data_model


# Hotspot parameters:
interpolated_hotspot_trails = True
use_trail_age_uncertainty = True

# Millions of years - e.g. 2 million years @ 50mm per year = 100km radius uncertainty ellipse
trail_age_uncertainty_ellipse = 1

include_chains = ['Louisville', 'Tristan', 'Reunion', 'St_Helena', 'Foundation', 'Cobb', 'Samoa', 'Tasmantid', 
                  'Hawaii']
#include_chains = ['Louisville', 'Tristan', 'Reunion', 'Hawaii', 'St_Helena', 'Tasmantid']


# Large area grid search to find minima
if search == 'Initial':

    search_type = 'Random'

# Uses grid search minima as seed for targeted secondary search (optional)
elif search == 'Secondary':

    search_type = 'Uniform'
    search_radius = 15
    rotation_uncertainty = 30

    models = 60
    auto_calc_ref_pole = False

# Used when auto_calc_ref_pole is False.
no_auto_ref_rot_longitude = -53.5
no_auto_ref_rot_latitude = 56.6
no_auto_ref_rot_angle = -2.28

interpolation_resolution = 5
rotation_age_of_interest = True

hst_file = 'HotspotTrails.geojson'
hs_file = 'HotspotCatalogue2.geojson'
interpolated_hotspots = 'interpolated_hotspot_chains_5Myr.xlsx'


# Don't plot in this workflow.
# This is so it can be run on an HPC cluster with no visualisation node.
plot = False
