#
# Supported parallelisation methods (or None to disable parallelisation, eg, for testing).
#
MPI4PY = 0
IPYPARALLEL = 1

# Choose parallelisation method (or None to disable parallelisation, eg, for testing).
use_parallel = MPI4PY


#
# Set model parameters, load data, and calculate starting conditions
# 
# Sets all user-selected parameters for the mode run
# 
# Arguments:
# 
#     geographical_uncertainty : Number that approximately represents geographical uncertainty - 95% confidence limit around ref pole location
#     rotation_uncertainty : Number that represents the upper and lower bounds of the optimisation's angle variation
#     sample_space : Selects the sampling method to be used to generate start seeds. "Fisher" (spherical distribution) by default.
#     models : The total number of models to be produced. 1 model = 1 complete optimisation from 1 starting location
#     model_stop_condition : Type of condition to be used to terminate optimisation. "threshold" or "max_iter". Threshold by default.
#     max_iter : IF "max_iter" selected, sets maximum iterations regardless of successful convergence.
#     ref_rotation_plate_id : Plate to be used as fixed reference. 701 (Africa) by default.
#     ref_rotation_start_age : Rotation begin age.
#     ref_rotation_end_age : Rotation end age.
#     interpolation_resolution : Resolution in degrees used in the Fracture Zone calulcations
#     rotation_age_of_interest : The result we are interested in. Allows for 'windowed mean' approach. It is the midpoint between the start and end ages by default.
# 
# Data to be included in optimisation: True by default.
# 
#     fracture_zones : Boolean
#     net_rotation : Boolean
#     trench_migration : Boolean
#     hotspot_reconstruction : Boolean
#     hotspot_dispersion : Boolean
# 

# The 'r' number is the Subversion revision number of deforming model 2016_v3.
model_name = "optAPM_run1"

# The data model to run the optimisation on.
# This should be the name of the sub-directory in 'data/' and also the sub-directory in 'data/TMData/'.
data_model = 'Muller++_2015_AREPS_CORRECTED'

start_age = 230
end_age = 0
interval = 10

search = "Initial"
search_radius = 60
# If True then temporarily expand search radius to 90 whenever the reference plate changes.
# Normally the reference plate stays constant at Africa (701), but does switch to 101 for the 1Ga model.
# It's off by default since it doesn't appear to change the results, and may sometimes cause job to fail
# on Artemis (presumably since 'models' is increased by a factor of 2.5) - although problem manifested
# as failure to read "all_rotations_<...>.rot" file, so it was probably something else.
expand_search_radius_on_ref_plate_switches = False
rotation_uncertainty = 30
auto_calc_ref_pole = True
models = 100

model_stop_condition = 'threshold'
max_iter = 5  # Only applies if model_stop_condition != 'threshold'


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
    else:
        # NOTE: These are inverse weights (ie, the constraint costs are *multiplied* by "1.0 / weight").
        return True, 2.0  # Gives a *multiplicative* weight of 0.5

def get_trench_migration_params(age):
    return True, 1.0

def get_hotspot_trail_params(age):
    # Only use hotspot trails for 0-80Ma.
    if age <= 80:
        return True, 1.0
    else:
        return False, 1.0


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


# Rotation file with existing APM rotations removed from 0-250Ma to be used:
if tm_data_type == 'Global_Model_WD_Internal_Release_2019_v1':

    original_rotfile = 'Global_Model_WD_Internal_Release_2019_v1/optimisation/all_rotations.rot'
    rotfile = 'Global_Model_WD_Internal_Release_2019_v1/optimisation/all_rotations_' + model_name + '.rot'

elif tm_data_type == 'Global_1000-0_Model_2017':

    original_rotfile = 'Global_1000-0_Model_2017/optimisation/all_rotations.rot'
    rotfile = 'Global_1000-0_Model_2017/optimisation/all_rotations_' + model_name + '.rot'

elif tm_data_type == 'Muller++_2015_AREPS_CORRECTED':

    original_rotfile = 'Muller++_2015_AREPS_CORRECTED/optimisation/all_rotations.rot'
    rotfile = 'Muller++_2015_AREPS_CORRECTED/optimisation/all_rotations_' + model_name + '.rot'

elif tm_data_type == 'muller2016':

    original_rotfile = 'Global_EarthByte_230-0Ma_GK07_AREPS.rot'
    rotfile = 'Global_EarthByte_230-0Ma_GK07_AREPS_' + model_name + '.rot'

elif tm_data_type == 'shephard2013':

    original_rotfile = 'Shephard_etal_ESR2013_Global_EarthByte_2013.rot'
    rotfile = 'Shephard_etal_ESR2013_Global_EarthByte_2013_' + model_name + '.rot'


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

#
# Which reference plate ID and PMAG rotation file to use at which age.
#
def get_reference_params(age):
    if tm_data_type == 'Global_1000-0_Model_2017':
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


if tm_method == 'convergence':

    if tm_data_type == 'muller2016':

        nnr_relative_datadir = 'TMData/'
        nnr_rotfile = 'Global_EarthByte_230-0Ma_GK07_AREPS_NNR.rot'


elif tm_method == 'pygplates':

    if tm_data_type == 'Global_Model_WD_Internal_Release_2019_v1':

        nnr_relative_datadir = 'TMData/Global_Model_WD_Internal_Release_2019_v1/'
        nnr_rotfile = 'Global_Model_WD_Internal_Release_2019_v1/optimisation/no_net_rotations.rot'

    elif tm_data_type == 'Global_1000-0_Model_2017':

        nnr_relative_datadir = 'TMData/Global_1000-0_Model_2017/'
        nnr_rotfile = 'Global_1000-0_Model_2017/optimisation/no_net_rotations.rot'

    elif tm_data_type == 'Muller++_2015_AREPS_CORRECTED':

        nnr_relative_datadir = 'TMData/Muller++_2015_AREPS_CORRECTED/'
        nnr_rotfile = 'Muller++_2015_AREPS_CORRECTED/optimisation/no_net_rotations.rot'

    elif tm_data_type == 'muller2016':

        nnr_relative_datadir = 'TMData/Muller_2016/'
        nnr_rotfile = 'Global_EarthByte_230-0Ma_GK07_AREPS_NNR.rot'

    elif tm_data_type == 'shephard2013':

        nnr_relative_datadir = 'TMData/Shephard_2013/'
        nnr_rotfile = 'Shephard_etal_ESR2013_Global_EarthByte_NNR_ORIGINAL.rot'


if tm_data_type == 'Global_Model_WD_Internal_Release_2019_v1':

    ridge_file = 'Global_Model_WD_Internal_Release_2019_v1/StaticGeometries/AgeGridInput/Global_EarthByte_GeeK07_Ridges_2019_v1.gpml'
    isochron_file = 'Global_Model_WD_Internal_Release_2019_v1/StaticGeometries/AgeGridInput/Global_EarthByte_GeeK07_Isochrons_2019_v1.gpml'
    isocob_file = 'Global_Model_WD_Internal_Release_2019_v1/StaticGeometries/AgeGridInput/Global_EarthByte_GeeK07_IsoCOB_2019_v1.gpml'

elif (tm_data_type == 'Global_1000-0_Model_2017' or
      tm_data_type == 'Muller++_2015_AREPS_CORRECTED'):

    ##################################################################################################################################
    #
    # There are no static geometries (besides coastlines) for this data model.
    #
    # NOTE: SO USING SAME FILES AS 'Global_Model_WD_Internal_Release_2019_v1'.
    #       THIS IS OK IF WE'RE NOT INCLUDING FRACTURE ZONES (BECAUSE THEN THESE FILES ARE NOT USED FOR FINAL OPTIMISED ROTATIONS).
    #
    ##################################################################################################################################
    ridge_file = 'Global_Model_WD_Internal_Release_2019_v1/StaticGeometries/AgeGridInput/Global_EarthByte_GeeK07_Ridges_2019_v1.gpml'
    isochron_file = 'Global_Model_WD_Internal_Release_2019_v1/StaticGeometries/AgeGridInput/Global_EarthByte_GeeK07_Isochrons_2019_v1.gpml'
    isocob_file = 'Global_Model_WD_Internal_Release_2019_v1/StaticGeometries/AgeGridInput/Global_EarthByte_GeeK07_IsoCOB_2019_v1.gpml'

else:

    ridge_file = 'Global_EarthByte_230-0Ma_GK07_AREPS_Ridges.gpml'
    isochron_file = 'Global_EarthByte_230-0Ma_GK07_AREPS_Isochrons.gpmlz'
    isocob_file = 'Global_EarthByte_230-0Ma_GK07_AREPS_IsoCOB.gpml'

hst_file = 'HotspotTrails.geojson'
hs_file = 'HotspotCatalogue2.geojson'
interpolated_hotspots = 'interpolated_hotspot_chains_5Myr.xlsx'


# Don't plot in this workflow.
# This is so it can be run on an HPC cluster with no visualisation node.
plot = False
