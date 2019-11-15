# Copyright (C) 2016  Michael G. Tetley
# EarthByte Group, University of Sydney / Seismological Laboratory, California Institute of Technology
# Contact email: michael.tetley@sydney.edu.au / mtetley@caltech.edu
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.


# Required libraries
import pygplates as pgp
import numpy as np
import pandas as pd

import pickle
import random
import textwrap
import geoTools
import pmagpy.ipmag as ipmag
import pmagpy.pmag as pmag
import isopolate
import math
import json

from HotSpotLoader import GetHotSpotTrailsFromGeoJSON, GetHotSpotLocationsFromGeoJSON

# The optimization workflow doesn't actually need to plot so we won't require user to install these modules.
# If the user plots then we'll get an AttributeError, in which case the try/except part should be removed.
try:
    import matplotlib.pyplot as plt
    from mpl_toolkits.basemap import Basemap
    from IPython.display import display, HTML
except ImportError:
    pass


class ModelSetup():


    @staticmethod
    def dataLoaderHelp():

        print textwrap.dedent("""
            Data Loader

            Loads all required data for optimisation routine.

            Arguments:
            ----------
            datadir : Absolute path to data directory
            rot_file : GPlates rotation *.rot file

            Optional arguments:

                nnr_rotfile : GPlates no net-rotation *.rot file. Required for net rotation calculations.
                tm_file : Trench migration data for the current reconstruction time. Required for trench migration calculations.
                pv_file : Plate velocity data (points and associated plate IDs) for the current reconstruction time.
                          Required for plate velocity calculations.

                ridge_file : GPlates *.gpml ridge file. Required for fracture zone calculations.
                isochron_file : GPlates *.gpml isochron file. Required for fracture zone calculations.
                isocob_file : GPlates *.gpml isocob file. Required for fracture zone calculations.

                hst_file : GeoJSON hotspot trail data. Relative path to datadir.
                hs_file : GeoJSON hotspot location data. Relative path to datadir.

            """)


    @staticmethod
    def dataLoader(datadir, rot_file, pmag_rot_file=None, tm_file=None, pv_file=None, nnr_rotfile=None,
                   ridge_file=None, isochron_file=None, isocob_file=None, hst_file=None, hs_file=None,
                   interpolated_hotspots=None):

        # Create rotation model
        rotation_file = datadir + rot_file
        rotation_model = pgp.RotationModel(rotation_file)

        # Create pmag rotation model
        if pmag_rot_file:

            pmag_rotation_file = datadir + pmag_rot_file
            pmag_rotation_model = pgp.RotationModel(pmag_rotation_file)

        # Check for and load optional arguments
        if tm_file:

            trench_migration_file = datadir + tm_file

        if pv_file:

            plate_velocity_file = datadir + pv_file

        if nnr_rotfile:

            no_net_rotation_file = datadir + nnr_rotfile

        if ridge_file:

            features_ri = pgp.FeatureCollection(datadir + ridge_file)
            RidgeFile_subset = pgp.FeatureCollection()

        if isochron_file:

            features_iso = pgp.FeatureCollection(datadir + isochron_file)
            IsochronFile_subset = pgp.FeatureCollection()

        if isocob_file:

            features_isocob = pgp.FeatureCollection(datadir + isocob_file)
            IsoCOBFile_subset = pgp.FeatureCollection()

        if hst_file:

            hst_datafile = hst_file
            hs_trails = GetHotSpotTrailsFromGeoJSON(datadir + hst_datafile)

        if hs_file:

            hs_datafile = hs_file
            hotspots = GetHotSpotLocationsFromGeoJSON(datadir + hs_datafile)

        if interpolated_hotspots:

            interpolated_hotspot_data = pd.read_excel(
                datadir + interpolated_hotspots,
                # pyGPlates (via boost.python) cannot convert numpy.int64 to regular Python float.
                # Not sure if happens only in 32-bit versions of pyGPlates (ie, on Windows).
                # So convert here instead...
                converters = {'Hotspot_lat' : float, 'Hotspot_lon' : float})


        # Return all loaded data
        data = [rotation_model, rotation_file, trench_migration_file, plate_velocity_file, no_net_rotation_file,
                features_ri, RidgeFile_subset, features_iso, IsochronFile_subset, features_isocob, IsoCOBFile_subset,
                hs_trails, hotspots, pmag_rotation_model, pmag_rotation_file, interpolated_hotspot_data]

        print "- Data loaded"
        print " "

        return data




    @staticmethod
    def modelParametersHelp():

        print textwrap.dedent("""
            Model Parameters

            Sets all user-selected parameters for the mode run

            Arguments:
            ----------
            geographical_uncertainty : Number that approximately represents geographical uncertainty - 95% confidence limit around ref pole location
            rotation_uncertainty : Number that represents the upper and lower bounds of the optimisation's angle variation
            sample_space : Selects the sampling method to be used to generate start seeds. "Fisher" (spherical distribution) by default.
            models : The total number of models to be produced. 1 model = 1 complete optimisation from 1 starting location
            model_stop_condition : Type of condition to be used to terminate optimisation. "threshold" or "max_iter". Threshold by default.
            max_iter : IF "max_iter" selected, sets maximum iterations regardless of successful convergence.
            ref_rotation_plate_id : Plate to be used as fixed reference. 701 (Africa) by default.
            ref_rotation_start_age : Rotation begin age.
            ref_rotation_end_age : Rotation end age.
            rotation_age_of_interest : The result we are interested in. Allows for 'windowed mean' approach. It is the midpoint between the start and end ages by default.

            Data to be included in optimisation: True by default.

            enable_fracture_zones : Boolean
            enable_net_rotation : Boolean
            enable_trench_migration : Boolean
            enable_hotspot_trails : Boolean
            enable_plate_velocity : Boolean

            """)


    @staticmethod
    def modelStartConditions(params, data, plot=True):

        # Translate data array
        (rotation_model,
         rotation_file,
         trench_migration_file,
         plate_velocity_file,
         no_net_rotation_file,
         features_ri,
         RidgeFile_subset,
         features_iso,
         IsochronFile_subset,
         features_isocob,
         IsoCOBFile_subset,
         hs_trails,
         hotspots,
         pmag_rotation_model,
         pmag_rotation_file,
         interpolated_hotspot_data) = data[:16]

        # Translate params array
        (geographical_uncertainty,
            rotation_uncertainty,
            search_type,
            models,
            model_stop_condition,
            max_iter,
            ref_rotation_plate_id,
            ref_rotation_start_age,
            ref_rotation_end_age,
            interpolation_resolution,
            rotation_age_of_interest,
            enable_fracture_zones, enable_net_rotation, enable_trench_migration, enable_hotspot_trails, enable_plate_velocity,
            fracture_zone_bounds, net_rotation_bounds, trench_migration_bounds, hotspot_trails_bounds, plate_velocity_bounds,
            ref_rot_longitude,
            ref_rot_latitude,
            ref_rot_angle,
            auto_calc_ref_pole,
            search,
            fz_weight,
            nr_weight,
            tm_weight,
            hs_weight,
            pv_weight,
            include_chains,
            interpolated_hotspot_trails,
            tm_method) = params[:34]

        #print fz_weight, nr_weight, tm_weight, hs_weight, pv_weight


        # Set rotation age of interest
        rotation_age_of_interest_age = ref_rotation_start_age - (0.5 * (ref_rotation_start_age - ref_rotation_end_age))

        # Optimisation data array
        data_array_labels = ['Fracture zone orientation', 'Net rotation', 'Trench migration',
                             'Hotspot trails', 'Plate velocity']

        weights_array = [fz_weight, nr_weight, tm_weight, hs_weight, pv_weight]

        data_array_labels_short = ['FZ', 'NR', 'TM', 'HS', 'PV']
        data_array = [enable_fracture_zones, enable_net_rotation, enable_trench_migration, enable_hotspot_trails, enable_plate_velocity]
        
        data_bounds = [fracture_zone_bounds, net_rotation_bounds, trench_migration_bounds, hotspot_trails_bounds, plate_velocity_bounds]

        # Array containing the name of all chains to be included in optimisation
        if ref_rotation_start_age <= 80:

            pass

            # chains including the pacific
            #include_chains = ['Reunion', 'Louisville', 'Tristan', 'St_Helena', 'Foundation', 'Samoa', 'Cobb', 'Caroline', 'Tasmantid']
            #print "Hotspot chains used:", str(include_chains)

        #else:

            # chains excluding the pacific
            #include_chains = ['Reunion', 'Tristan', 'St_Helena', 'Tasmantid']
            #print "Hotspot chains used:", str(include_chains)

        print "Optimisation parameters:"
        print " "
        print "Data constraints:"

        for i in xrange(0, len(data_array)):

            if data_array[i] == True:

                print "- " + data_array_labels[i] + ": weight(" + str(weights_array[i]) + ")"
                
                if data_bounds[i]:
                    print "  - bounds" + str(data_bounds[i])

        print " "
        print "Termination:"

        if model_stop_condition == "max_iter":

            print "- Maximum iteration: " + str(max_iter)

        else:

            print "- Threshold"


        print " "
        print "Sampling method:"
        print "- " + search_type
        print " "
        print "Reference rotation type:"


        # Prepare rotation model for updates during optimisation - keeps rotations in memory
        # rotation_model_tmp = pgp.FeatureCollection(rotation_file)


        # Calculate initial reference rotation for the reference plate (eg, Africa) from selected rotation model
        if auto_calc_ref_pole == True:

            print "- Auto-calc palaeomagnetic"

            ref_rot = pmag_rotation_model.get_rotation(np.double(ref_rotation_start_age), ref_rotation_plate_id, 0)
            ref_rot_pole, ref_rot_angle = ref_rot.get_euler_pole_and_angle()

            ref_rot_of_interest = pmag_rotation_model.get_rotation(np.double(rotation_age_of_interest_age), ref_rotation_plate_id, 0)
            ref_rot_of_interest_pole, ref_rot_of_interest_angle = ref_rot_of_interest.get_euler_pole_and_angle()

            # Convert finite rotations to lat, lon and degrees
            ref_rot_pole = pgp.convert_point_on_sphere_to_lat_lon_point(ref_rot_pole)

            ref_rot_longitude = ref_rot_pole.get_longitude()
            ref_rot_latitude = ref_rot_pole.get_latitude()
            ref_rot_angle = np.rad2deg(ref_rot_angle)

        elif auto_calc_ref_pole == False:

            print "- User reference"

            ref_rot_longitude = ref_rot_longitude
            ref_rot_latitude = ref_rot_latitude
            ref_rot_angle = ref_rot_angle

            ref_rot_of_interest = pmag_rotation_model.get_rotation(np.double(rotation_age_of_interest_age), ref_rotation_plate_id, 0)
            ref_rot_of_interest_pole, ref_rot_of_interest_angle = ref_rot_of_interest.get_euler_pole_and_angle()


        # else:

        #     print "Reference rotation type: EarthByte model"

        #     ref_rot = rotation_model.get_rotation(np.double(ref_rotation_start_age), ref_rotation_plate_id, 0)
        #     ref_rot_pole, ref_rot_angle = ref_rot.get_euler_pole_and_angle()

        #     ref_rot_of_interest = rotation_model.get_rotation(np.double(rotation_age_of_interest_age), ref_rotation_plate_id, 0)
        #     ref_rot_of_interest_pole, ref_rot_of_interest_angle = ref_rot_of_interest.get_euler_pole_and_angle()


        # Convert finite rotations to lat, lon and degrees
        # ref_rot_pole = pgp.convert_point_on_sphere_to_lat_lon_point(ref_rot_pole)

        # ref_rot_longitude = ref_rot_pole.get_longitude()
        # ref_rot_latitude = ref_rot_pole.get_latitude()
        # ref_rot_angle = np.rad2deg(ref_rot_angle)

        print " "
        print "Reference finite rotation pole for reference plate", ref_rotation_plate_id, "at", ref_rotation_start_age, "Ma:"
        print "- Lon:", ref_rot_longitude
        print "- Lat:", ref_rot_latitude
        print "- Angle:", ref_rot_angle
        print " "


        ref_rot_of_interest_pole = pgp.convert_point_on_sphere_to_lat_lon_point(ref_rot_of_interest_pole)

        ref_rot_longitude_of_interest = ref_rot_of_interest_pole.get_longitude()
        ref_rot_latitude_of_interest = ref_rot_of_interest_pole.get_latitude()
        ref_rot_angle_of_interest = np.rad2deg(ref_rot_of_interest_angle)

        # print "Reference finite rotation pole for reference plate", ref_rotation_plate_id, "at", rotation_age_of_interest_age, "Ma:"
        # print "- Lon:", ref_rot_longitude_of_interest
        # print "- Lat:", ref_rot_latitude_of_interest
        # print "- Angle:", ref_rot_angle_of_interest


        # Calculate start seeds (reference + (models - 1) random starting rotations within uncertainty limits)
        start_seeds_rotated = []
        start_seeds = []
        seed_history = []
        seed_lons = []
        seed_lats = []
        seed_angs = []


        # Generate uniform distribution of start seeds
        if search_type == 'Random':

            if search == "Initial":

                num_points = models * 5
                #num_points = models

            lons = []
            lats = []

            for i in xrange(0, num_points):

                theta = 2 * np.pi * np.random.random()
                phi = np.arccos(2 * np.random.random() - 1.0)

                x = np.cos(theta) * np.sin(phi)
                y = np.sin(theta) * np.sin(phi)
                z = np.cos(phi)

                point = pgp.convert_point_on_sphere_to_lat_lon_point((x,y,z))
                lats.append(point.get_latitude())
                lons.append(point.get_longitude())

            # Extract points within latitudinal zone of interest
            sample_lats = []
            sample_lons = []

            for i in xrange(0, len(lats)):

                if lats[i] > 90 - geographical_uncertainty:

                    sample_lats.append(lats[i])
                    sample_lons.append(lons[i])

            #print len(sample_lats)

            # Rotate points from pole to reference start seed
            for i in xrange(0, len(sample_lats)):

                start_seeds_rotated.append(pmag.dodirot(sample_lons[i], sample_lats[i], ref_rot_longitude, ref_rot_latitude))

            # Trim start seed array to match number of requested models
            if len(start_seeds_rotated) > models:

                trim = len(start_seeds_rotated) - models
                start_seeds_rotated = start_seeds_rotated[:-trim]

            # Sample gaussian array of angles using standard deviation of ref angle and uncertainty limits
            ang_array = [(ref_rot_angle - rotation_uncertainty), ref_rot_angle, (ref_rot_angle + rotation_uncertainty)]
            ang_array_sd = np.std(ang_array)

            ang_gaussian_array = []
            ang_gaussian_array.append(np.random.normal(ref_rot_angle, ang_array_sd, models))

            # Create start seeds array
            for i in xrange(0, len(start_seeds_rotated)):

                seed = [[start_seeds_rotated[i][0]], [start_seeds_rotated[i][1]], ref_rot_angle]

                seed_lons.append(start_seeds_rotated[i][0])
                seed_lats.append(start_seeds_rotated[i][1])
                seed_angs.append(ref_rot_angle)

                if seed not in seed_history:

                    start_seeds.append(seed)
                    seed_history.append(seed)

            #print len(start_seeds_rotated)
            #print len(start_seeds)
            #print len(seed_history)


            # For secondary minimisation, add inital min result to start seeds
            if auto_calc_ref_pole == False:

                seed = [[ref_rot_longitude], [ref_rot_latitude], ref_rot_angle]
                seed_lons.append(ref_rot_longitude)
                seed_lats.append(ref_rot_latitude)
                seed_angs.append(ref_rot_angle)

                if seed not in seed_history:

                    start_seeds.append(seed)
                    seed_history.append(seed)



            # (Lon, Lat, Ang)
            start_seeds = np.array(start_seeds)
            #print start_seeds
            #print len(start_seeds)


            # Set N to number of models to be run (optimisations)
            N = len(start_seeds[0])


            # Set lon, lat and angle error bounds.
            lb = []
            lb.append(np.min(seed_lons) - geographical_uncertainty)
            lb.append(np.min(seed_lats) - geographical_uncertainty)
            lb.append(np.min(seed_angs) - rotation_uncertainty)

            ub = []
            ub.append(np.max(seed_lons) + geographical_uncertainty)
            ub.append(np.max(seed_lats) + geographical_uncertainty)
            ub.append(np.max(seed_angs) + rotation_uncertainty)

            # print 'max_lons', np.max(seed_lons)
            # print 'max_lats', np.max(seed_lats)
            # print 'min_lons', np.min(seed_lons)
            # print 'min_lats', np.min(seed_lats)
            # print " "
            # print seed_lons
            # print " "



        if search_type == 'Uniform':

            if search == "Initial":

                # Works for search radius of 60
                #num_points = models * 4
                num_points = models * 5

            elif search == "Secondary":

                if geographical_uncertainty == 30:

                    num_points = models * 15

                elif geographical_uncertainty == 15:

                    num_points = models * 60

            # angle = np.pi * (3 - np.sqrt(5))
            # theta = angle * np.arange(num_points)
            # z = np.linspace(1 - 1.0 / num_points, 1.0 / num_points - 1, num_points)
            # radius = np.sqrt(1 - z * z)
            #
            # points = np.zeros((num_points, 3))
            # points[:,0] = radius * np.cos(theta)
            # points[:,1] = radius * np.sin(theta)
            # points[:,2] = z
            #
            # lons = []
            # lats = []
            #
            # for i in xrange(0, len(points)):
            #
            #     point = pgp.convert_point_on_sphere_to_lat_lon_point((points[i][0], points[i][1], points[i][2]))
            #     lats.append(point.get_latitude())
            #     lons.append(point.get_longitude())

            lons = []
            lats = []

            for i in xrange(0, num_points):

                theta = 2 * np.pi * np.random.random()
                phi = np.arccos(2 * np.random.random() - 1.0)

                x = np.cos(theta) * np.sin(phi)
                y = np.sin(theta) * np.sin(phi)
                z = np.cos(phi)

                point = pgp.convert_point_on_sphere_to_lat_lon_point((x,y,z))
                lats.append(point.get_latitude())
                lons.append(point.get_longitude())

            # Extract points within latitudinal zone of interest
            sample_lats = []
            sample_lons = []

            for i in xrange(0, len(lats)):

                if lats[i] > 90 - geographical_uncertainty:

                    sample_lats.append(lats[i])
                    sample_lons.append(lons[i])

            #print len(sample_lats)

            # Rotate points from pole to reference start seed
            for i in xrange(0, len(sample_lats)):

                start_seeds_rotated.append(pmag.dodirot(sample_lons[i], sample_lats[i], ref_rot_longitude, ref_rot_latitude))

            # Sample gaussian array of angles using standard deviation of ref angle and uncertainty limits
            ang_array = [(ref_rot_angle - rotation_uncertainty), ref_rot_angle, (ref_rot_angle + rotation_uncertainty)]
            ang_array_sd = np.std(ang_array)

            ang_gaussian_array = []
            ang_gaussian_array.append(np.random.normal(ref_rot_angle, ang_array_sd, models))

            # Create start seeds array
            for i in xrange(0, len(start_seeds_rotated)):

                seed = [[start_seeds_rotated[i][0]], [start_seeds_rotated[i][1]], ref_rot_angle]
                #seed = [[start_seeds_rotated[i][0]], [start_seeds_rotated[i][1]], ang_gaussian_array[0][i]]
                seed_lons.append(start_seeds_rotated[i][0])
                seed_lats.append(start_seeds_rotated[i][1])
                seed_angs.append(ref_rot_angle)
                #seed_angs.append(ang_gaussian_array[0][i])

                if seed not in seed_history:

                    start_seeds.append(seed)
                    seed_history.append(seed)

            #print len(start_seeds_rotated)
            #print len(start_seeds)
            #print len(seed_history)


            # For secondary minimisation, add inital min result to start seeds
            if auto_calc_ref_pole == False:

                seed = [[ref_rot_longitude], [ref_rot_latitude], ref_rot_angle]
                seed_lons.append(ref_rot_longitude)
                seed_lats.append(ref_rot_latitude)
                seed_angs.append(ref_rot_angle)

                if seed not in seed_history:

                    start_seeds.append(seed)
                    seed_history.append(seed)



            # (Lon, Lat, Ang)
            start_seeds = np.array(start_seeds)
            #print start_seeds
            #print len(start_seeds)


            # Set N to number of models to be run (optimisations)
            N = len(start_seeds[0])


            # Set lon, lat and angle error bounds.
            lb = []
            lb.append(np.min(seed_lons) - geographical_uncertainty)
            lb.append(np.min(seed_lats) - geographical_uncertainty)
            lb.append(np.min(seed_angs) - rotation_uncertainty)

            ub = []
            ub.append(np.max(seed_lons) + geographical_uncertainty)
            ub.append(np.max(seed_lats) + geographical_uncertainty)
            ub.append(np.max(seed_angs) + rotation_uncertainty)

            # print 'max_lons', np.max(seed_lons)
            # print 'max_lats', np.max(seed_lats)
            # print 'min_lons', np.min(seed_lons)
            # print 'min_lats', np.min(seed_lats)
            # print " "
            # print seed_lons
            # print " "

        # Generate a gaussian distribution of start seeds
        elif search_type == 'Fisher':

            # Calculate Fisher kappa distribution on sphere to represent geographical uncertainty
            if models == 100:

                k = geoTools.calcKfromA95(geographical_uncertainty / 9, models)

            elif models == 1000:

                 k = geoTools.calcKfromA95(geographical_uncertainty / 22, models)

            else:

                k = geoTools.calcKfromA95(geographical_uncertainty / 6, models)

            # Draw fisher distributed start seeds [lon lat] around pole
            if auto_calc_ref_pole == True:

                fdist = ipmag.fishrot(k, models, 0)

            elif auto_calc_ref_pole == False:

                # Subtract array with models - 1 as initial optimised minimum result is appended
                fdist = ipmag.fishrot(k, models - 1, 0)

            # Rotate start seeds to center on reference location
            for i in xrange(0, len(fdist)):

                start_seeds_rotated.append(pmag.dodirot(fdist[i][0], fdist[i][1], ref_rot_longitude, ref_rot_latitude))

            # Sample gaussian array of angles using standard deviation of ref angle and uncertainty limits
            ang_array = [(ref_rot_angle - rotation_uncertainty), ref_rot_angle, (ref_rot_angle + rotation_uncertainty)]
            ang_array_sd = np.std(ang_array)

            ang_gaussian_array = []
            ang_gaussian_array.append(np.random.normal(ref_rot_angle, ang_array_sd, models))


            # Create start seeds array
            for i in xrange(0, len(start_seeds_rotated)):

                if search == "Secondary":

                    seed = [[start_seeds_rotated[i][0]], [start_seeds_rotated[i][1]], ref_rot_angle]

                else:

                    seed = [[start_seeds_rotated[i][0]], [start_seeds_rotated[i][1]], ang_gaussian_array[0][i]]

                seed_lons.append(start_seeds_rotated[i][0])
                seed_lats.append(start_seeds_rotated[i][1])
                seed_angs.append(ang_gaussian_array[0][i])

                if seed not in seed_history:

                    start_seeds.append(seed)
                    seed_history.append(seed)

            # For secondary minimisation, add inital min result to start seeds
            if auto_calc_ref_pole == False:

                seed = [[ref_rot_longitude], [ref_rot_latitude], ref_rot_angle]
                seed_lons.append(ref_rot_longitude)
                seed_lats.append(ref_rot_latitude)
                seed_angs.append(ref_rot_angle)

                if seed not in seed_history:

                    start_seeds.append(seed)
                    seed_history.append(seed)


            # (Lon, Lat, Ang)
            start_seeds = np.array(start_seeds)
            #print start_seeds
            #print len(start_seeds)


            # Set N to number of models to be run (optimisations)
            N = len(start_seeds[0])


            # Set lon, lat and angle error bounds.
            lb = []
            lb.append(np.min(seed_lons) - geographical_uncertainty)
            lb.append(np.min(seed_lats) - geographical_uncertainty)
            lb.append(np.min(seed_angs) - rotation_uncertainty)

            ub = []
            ub.append(np.max(seed_lons) + geographical_uncertainty)
            ub.append(np.max(seed_lats) + geographical_uncertainty)
            ub.append(np.max(seed_angs) + rotation_uncertainty)



        elif search_type == '2D_cartesian':

            # Using a while loop here ensures total start seeds = model iterations even if duplicates are found
            while len(start_seeds) != models:

                tmp_seed = []

                # Set first seed to known reference plate rotation
                if len(start_seeds) < 1:

                    lon_rand = ref_rot_longitude
                    lat_rand = ref_rot_latitude
                    angle_rand = ref_rot_angle

                    seed = [[lon_rand], [lat_rand], [angle_rand]]

                    start_seeds.append(seed)
                    seed_history.append(seed)

                elif len(start_seeds) >= 1:

                    # Generate random start seeds from geographic uncertainty ranges for initial longitude and latitude
                    lon_rand = random.sample(range((int(ref_rot_longitude) - int(geographical_uncertainty)),
                                                   (int(ref_rot_longitude) + int(geographical_uncertainty))), 1)

                    lat_rand = random.sample(range((int(ref_rot_latitude) - int(geographical_uncertainty)),
                                                   (int(ref_rot_latitude) + int(geographical_uncertainty))), 1)

                    angle_rand = random.sample(range((int(ref_rot_angle) - int(rotation_uncertainty)),
                                                     (int(ref_rot_angle) + int(rotation_uncertainty))), 1)


                    seed = [lon_rand, lat_rand, angle_rand]


                    # Check to make sure latitude seed is > 0 and < 90
                    if seed[1][0] >= 90:

                        seed[1][0] = 89

                    elif seed[1][0] <= -90:

                        seed[1][0] = 1


                    # Check all seeds are within upper and lower bounds
                    if seed[0][0] >= ref_rot_longitude - geographical_uncertainty and \
                       seed[0][0] <= ref_rot_longitude + geographical_uncertainty and \
                       seed[1][0] >= ref_rot_latitude - geographical_uncertainty and \
                       seed[1][0] <= ref_rot_latitude + geographical_uncertainty and \
                       seed[2][0] >= ref_rot_angle - rotation_uncertainty and \
                       seed[2][0] <= ref_rot_angle + rotation_uncertainty:

                        if seed not in seed_history:

                            start_seeds.append(seed)
                            seed_history.append(seed)

                    else:

                        continue



            # (Lon, Lat, Ang)
            start_seeds = np.array(start_seeds)


            # Set N to number of models to be run (optimisations)
            N = len(start_seeds[0])


            # Set lon, lat and angle error bounds.
            lb = []
            lb.append(ref_rot_longitude - geographical_uncertainty)
            lb.append(ref_rot_latitude - geographical_uncertainty)
            lb.append(ref_rot_angle - rotation_uncertainty)

            ub = []
            ub.append(ref_rot_longitude + geographical_uncertainty)
            ub.append(ref_rot_latitude + geographical_uncertainty)
            ub.append(ref_rot_angle + rotation_uncertainty)




        # Append lon, lat and angle to make x = 3N elements for all start seeds (optimisation vector input)
        x = []
        x_ = []

        for i in xrange(0, len(start_seeds)):

            x_ = []
            x_.append(np.append(start_seeds[i][0], start_seeds[i][1]))
            x.append(np.append(x_, start_seeds[i][2]))


        opt_n = len(x[0])

        hs_eval_data = []




        # Build feature subset lists
        for feature in features_iso:
            if feature.get_valid_time()[1] <= ref_rotation_start_age:
                IsochronFile_subset.add(feature)

        for feature in features_ri:
            if feature.get_valid_time()[1] <= ref_rotation_start_age:
                RidgeFile_subset.add(feature)

        for feature in features_isocob:
            if feature.get_valid_time()[1] <= ref_rotation_start_age:
                IsoCOBFile_subset.add(feature)



        # This is where intertec is actually called - note that these one line interpolates the isochrons,
        # but does not reconstruct them
        recon_interpolated_isochrons = []

        output_features = isopolate.interpolate_isochrons(rotation_model,
                          [RidgeFile_subset, IsochronFile_subset, IsoCOBFile_subset],
                          interpolate=np.arange(ref_rotation_end_age, ref_rotation_start_age + 1, interpolation_resolution),
                          #interpolate=0.01,
                          tessellate_threshold_radians=math.radians(0.5),
                          output_scalar_spreading_direction=True,
                          output_scalar_spreading_rate=True,
                          output_scalar_spreading_asymmetry=True,
                          output_scalar_age=True,
                          print_debug_output=0)

        # Here we do the last step that the old intertec did, reconstructing the interpolated
        # isochrons to the selected time.
        pgp.reconstruct(output_features, rotation_file, recon_interpolated_isochrons, ref_rotation_start_age, 0)

        ## Step2
        # Take the coverage that intertec produced, and use it to derive arrays of data more easily analysed using numpy commands

        Lats = []
        Lons = []
        spreading_directions = []
        spreading_rates = []
        spreading_asymmetries = []
        seafloor_ages = []
        PID = []
        CPID = []

        for recon_interpolated_isochron in recon_interpolated_isochrons:

            tmp = recon_interpolated_isochron.get_feature()
            tmp2 = tmp.get_geometry(coverage_return=pgp.CoverageReturn.geometry_and_scalars)

            if tmp2:

                coverage_geometry, coverage_scalars = tmp2
                coverage_points = coverage_geometry.get_points()
                #tmp3 = coverage_points.get_points()

                for scalar in coverage_scalars.get(pgp.ScalarType.create_gpml('SpreadingDirection')):

                    spreading_directions.append(scalar)
                    #spreading_directions.append(coverage_scalars.get(pgp.ScalarType.create_gpml('SpreadingDirection')))

                for scalar in coverage_scalars.get(pgp.ScalarType.create_gpml('SpreadingRate')):

                    spreading_rates.append(scalar)

                for scalar in coverage_scalars.get(pgp.ScalarType.create_gpml('SpreadingAsymmetry')):

                    spreading_asymmetries.append(scalar)

                for scalar in coverage_scalars.get(pgp.ScalarType.create_gpml('Age')):

                    seafloor_ages.append(scalar)

                for point in coverage_points.get_points().to_lat_lon_array():

                    Lats.append(point[0])
                    Lons.append(point[1])
                    PID.append(tmp.get_reconstruction_plate_id())
                    CPID.append(tmp.get_conjugate_plate_id())

        # TRENCH MIGRATION
        # Load pre-computed migration data if needed (net rotation or trench migration)

        if tm_method == 'convergence':
            
            # DEPRECATED: We currently use "tm_method == 'convergence'" so this code is essentially deprecated.
            #             However, if it's brought then we need access to 'datadir' and 'data_model' (via function arguments).
            raise DeprecationWarning("Using the old 'convergence' of trench statistics is deprecated - use 'pygplates' instead.")
            nnr_datadir = datadir + 'TMData/' + data_model + '/'
            
            FNAME = nnr_datadir + 'data_%s_%s.txt' % (int(ref_rotation_start_age), int(ref_rotation_end_age))
            with open(FNAME, 'r') as f:
                reformArray = json.load(f)

        elif tm_method == 'pygplates':

            reformArray = []


        # HOTSPOT CHAINS
        # Get hotspot chain data into Dataframe
        trails_list = []

        for item in hs_trails:

            trails_list.append({ 'Chain': str(item) , 'Lon': hs_trails[item]['lon'], 'Lat': hs_trails[item]['lat'],
                                 'Age': hs_trails[item]['age'], 'Age_error': hs_trails[item]['age_error'],
                                 'PlateID': hs_trails[item]['plateid'], 'Geo_error': geographical_uncertainty,
                                 'Rot_error': rotation_uncertainty})

        hs_trailsDF = pd.DataFrame(trails_list)


        # Add geographical location of each hotspot to existing chains list
        for i, item in enumerate(hs_trailsDF['Chain']):

            for j, hotspot in enumerate(hotspots):

                if item == hotspot:

                    trails_list[i].update({'Hotspot': (hotspots[hotspot]['lat'], hotspots[hotspot]['lon'])})

        hs_trailsDF = pd.DataFrame(trails_list)


        # Build trail data subset
        if interpolated_hotspot_trails == False:

            print '- Using raw hotspot trail data'

            trail_data = []

            for i, item in enumerate(hs_trailsDF.Chain):

                if item in include_chains:

                    #print item

                    trail_data.append({'Chain': hs_trailsDF.Chain[i], 'Lon': hs_trailsDF.Lon[i], 'Lat': hs_trailsDF.Lat[i],
                                       'Age': hs_trailsDF.Age[i], 'Age_err': hs_trailsDF.Age_error[i], 'PlateID': hs_trailsDF.PlateID[i],
                                       'Hotspot_Location': hs_trailsDF.Hotspot[i], 'Geo_err': hs_trailsDF.Geo_error[i],
                                       'Rot_err': hs_trailsDF.Rot_error[i]})

        if interpolated_hotspot_trails == True:

            print '- Using interpolated hotspot trail data'

            trail_data = []

            for i, item in enumerate(include_chains):

                trail_lons = []
                trail_lats = []
                trail_ages = []
                trail_plateID = []

                for j in xrange(0, len(interpolated_hotspot_data)):

                    if interpolated_hotspot_data['Chain'][j] == item:

                        trail_lons.append(interpolated_hotspot_data['Lon'][j])
                        trail_lats.append(interpolated_hotspot_data['Lat'][j])
                        trail_ages.append(interpolated_hotspot_data['Age'][j])
                        trail_plateID.append(int(interpolated_hotspot_data['PlateID'][j]))
                        trail_hotspot = [interpolated_hotspot_data['Hotspot_lat'][j], interpolated_hotspot_data['Hotspot_lon'][j]]

                trail_data.append({'Chain': include_chains[i], 'Lon': tuple(trail_lons), 'Lat': tuple(trail_lats), 'Age': tuple(trail_ages),
                                   'PlateID': tuple(trail_plateID), 'Hotspot_Location': tuple(trail_hotspot)})


        print "- Created optimisation input vector"
        print "- Upper and lower bounds calculated"
        # print "- Precomputed isopolate output"
        # print "- Generated hotspot and trail data sets"
        print " "
        print "Start seeds generated:", len(seed_lons)



        if plot:
            # Plot start seeds
            m = Basemap(projection='robin',lat_0=ref_rot_latitude,lon_0=ref_rot_longitude,resolution='c',area_thresh=50000)
            plt.figure(figsize=(7, 7))
            #plt.title("Start seed distribution")
            m.drawcoastlines(linewidth=0.25)
            m.fillcontinents(color='bisque',zorder=1)
            m.drawmeridians(np.arange(0,360,30))
            m.drawparallels(np.arange(-90,90,30))

            ipmag.plot_vgp(m, seed_lons, seed_lats)
            ipmag.plot_pole(m, ref_rot_longitude, ref_rot_latitude, geographical_uncertainty, color='red')


        startingConditions = [x, opt_n, N, lb, ub, model_stop_condition, max_iter,
                              rotation_file, ref_rotation_start_age, ref_rotation_end_age, ref_rotation_plate_id,
                              Lats, Lons, spreading_directions, spreading_asymmetries, seafloor_ages, PID, CPID,
                              data_array, data_bounds, trench_migration_file, plate_velocity_file, no_net_rotation_file,
                              reformArray, trail_data, start_seeds, rotation_age_of_interest_age, data_array_labels_short,
                              ref_rot_longitude, ref_rot_latitude, ref_rot_angle,
                              seed_lons, seed_lats, ang_gaussian_array]


        return startingConditions



class ProcessResults():

    @staticmethod
    def saveResultsToPickle(data_array, data_array_labels_short, ref_rotation_start_age, ref_rotation_end_age, geographical_uncertainty, xopt, models, model_name):

        constraints = ''

        for i in xrange(0, len(data_array)):

            if data_array[i] == True:

                constraints = constraints + '_' + str(data_array_labels_short[i])


        # Save pickle file
        output_file = "model_output/" + model_name + "_" + str(int(ref_rotation_start_age)) + "-" + str(int(ref_rotation_end_age)) + "Ma_" \
        + str(models) + "models" + constraints + "_" + str(geographical_uncertainty) + ".pkl"

        pickle.dump(xopt, open(output_file, 'wb'))

        print "Results saved to:", output_file
        print " "

        return output_file



    @staticmethod
    def sortAndPlot(output_file, ref_rotation_start_age, ref_rotation_end_age, rotation_age_of_interest_age, xopt, rotation_file, ref_rot_longitude,
                    ref_rot_latitude, ref_rot_angle, seed_lons, seed_lats, ref_rotation_plate_id, model_name, models, data_array_labels_short, data_array,
                    geographical_uncertainty,
                    plot=True):

        # Load saved results
        xopt_ = pickle.load(open(output_file, 'rb'))

        results_dict = []
        stage_results_dict = []

        for i in xrange(0, len(xopt_)):

            # xopt returns as tuple (immutable) - convert to list for wrangling
            xopt_[i] = list(xopt_[i])

            # Assign a model number (zero-indexed) - the check is to make sure this happens only once if cell is re-run
            if len(xopt[i]) == 2:

                xopt_[i].append(i)

            # if xopt_[i][0][0] > 180:
            #
            #     xopt_lon = -360 + xopt_[i][0][0]
            #
            # if xopt_[i][0][0] < -180:
            #
            #     #xopt_lon = -180 + xopt_[i][0][0]
            #     xopt_lon = xopt_[i][0][0] + 360
            #
            # if xopt_[i][0][0] < 180 and xopt_[i][0][0] > -180:

            xopt_lon = xopt_[i][0][0]

            # Create dictionary (within a list) for dataframe
            results_dict.append({'Lon': xopt_lon, 'Lat': xopt_[i][0][1], 'Ang': xopt_[i][0][2], 'Minimum': xopt_[i][1],
                                 'Model': xopt_[i][2], 'Age': ref_rotation_start_age})

            # Create dictionary for stage rotation calculated at rotation_age_of_interest_age for dataframe
            stage_results_dict.append({'Lon': xopt_lon, 'Lat': xopt_[i][0][1],
                                       'Ang': xopt_[i][0][2] / ref_rotation_start_age * \
                                             (ref_rotation_start_age - rotation_age_of_interest_age),
                                       'Minimum': xopt_[i][1], 'Model': xopt_[i][2], 'Age': rotation_age_of_interest_age})


        # Create pandas dataframe for all results to sort ascending by 'Minimum' column
        results_dataframe = pd.DataFrame(results_dict)
        stage_results_dataframe = pd.DataFrame(stage_results_dict)

        results_dataframe = results_dataframe.sort_values(by='Minimum')
        stage_results_dataframe = stage_results_dataframe.sort_values('Minimum')

        # Get top 10% minimum results
        top_ten = results_dataframe.head(int(0.1 * len(results_dataframe)))
        stage_top_ten = stage_results_dataframe.head(int(0.1 * len(stage_results_dataframe)))

        # Get mean of top 10%
        top_ten_mean = top_ten.mean()
        top_ten_mean = pd.Series.to_frame(top_ten_mean)
        top_ten_mean = top_ten_mean.transpose()

        stage_top_ten_mean = stage_top_ten.mean()
        stage_top_ten_mean = pd.Series.to_frame(stage_top_ten_mean)
        stage_top_ten_mean = stage_top_ten_mean.transpose()

        # Get absolute minimum results
        min_result = results_dataframe.head(1)
        stage_min_result = results_dataframe.head(1)

        # Get reference rotation from EarthByte model
        ref_result = results_dataframe.loc[results_dataframe['Model'] == 0]


        if plot:
            
            # Plot results
            fig = plt.figure(figsize=(30,22),dpi=150)

            pmap = Basemap(projection='robin', lat_0=0, lon_0=0, resolution='l')
            pmap.drawmapboundary(fill_color='white')
            pmap.drawmeridians(np.arange(-180, 180, 30), labels=[0,0,0,1], fontsize=12)
            pmap.drawparallels(np.arange(-90, 90, 30), labels=[1,0,0,0], fontsize=12)

            pmap.fillcontinents(color='lightgray')

            plt.title('Optimised Africa rotations compared to EarthByte reference rotation for ' \
                      + str(ref_rotation_start_age) + ' - ' + str(ref_rotation_end_age) + ' Ma')

            # Convert results to basmap coords
            xall, yall = pmap(np.array(results_dataframe['Lon']), np.array(results_dataframe['Lat']))
            xten, yten = pmap(np.array(top_ten['Lon']), np.array(top_ten['Lat']))
            xmean, ymean = pmap(np.float(top_ten_mean['Lon']), np.float(top_ten_mean['Lat']))
            xmin, ymin = pmap(np.float(min_result['Lon']), np.float(min_result['Lat']))
            xref, yref = pmap(ref_rot_longitude, ref_rot_latitude)
            slon, slat = pmap(seed_lons, seed_lats)

            point_multiplier = 200


            # Plot all results
            for i in xrange(0, len(results_dataframe)):

                if results_dataframe.iloc[i]['Ang'] < 0:

                    pmap.scatter(xall[i], yall[i], marker='o', s=abs(results_dataframe.iloc[i]['Ang']) * point_multiplier,
                                 edgecolor='gray', zorder=2, label='All results' if i == 0 else "", facecolor='none',
                                 linewidth=1)

                    pmap.scatter(xall[i], yall[i], marker='.', s=5, color='gray', zorder=2)

                elif results_dataframe.iloc[i]['Ang'] > 0:

                    pmap.scatter(xall[i], yall[i], marker='o', s=abs(results_dataframe.iloc[i]['Ang']) * point_multiplier, c='gray',
                                 edgecolor='darkgray', zorder=2, label='All results' if i == 0 else "")

                    pmap.scatter(xall[i], yall[i], marker='.', s=5, color='white', zorder=2)

            # Plot top 10%
            for i in xrange(0, len(top_ten)):

                if top_ten.iloc[i]['Ang'] < 0:

                    pmap.scatter(xten[i], yten[i], marker='o', s=abs(top_ten.iloc[i]['Ang']) * point_multiplier, c='green', zorder=2,
                                 label='10%' if i == 0 else "", edgecolor='seagreen', facecolor='none', linewidth=2)

                    pmap.scatter(xten[i], yten[i], marker='.', s=5, c='green', zorder=2)


                elif top_ten.iloc[i]['Ang'] > 0:

                    pmap.scatter(xten[i], yten[i], marker='o', s=abs(top_ten.iloc[i]['Ang']) * point_multiplier, c='green', zorder=2,
                                 label='10%' if i == 0 else "")

                    pmap.scatter(xten[i], yten[i], marker='.', s=5, c='k', zorder=2)


            # Plot top 10% mean
            if np.float(top_ten_mean['Ang']) < 0:

                pmap.scatter(xmean, ymean, marker='o', s=abs(top_ten_mean['Ang']) * point_multiplier, c='orange', zorder=2,
                             label='10% mean', facecolor='none', edgecolor='orange', linewidth=2)

                pmap.scatter(xmean, ymean, marker='.', s=5, c='orange', zorder=2)

            elif np.float(top_ten_mean['Minimum']) > 0:

                pmap.scatter(xmean, ymean, marker='o', s=abs(top_ten_mean['Ang']) * point_multiplier, c='orange', zorder=2,
                             label='10% mean')

                pmap.scatter(xmean, ymean, marker='.', s=5, c='k', zorder=2)


            # Plot min result
            if np.float(min_result['Ang']) < 0:

                pmap.scatter(xmin, ymin, marker='o', s=abs(np.float(min_result['Ang'])) * point_multiplier, c='blue', zorder=2,
                             label='Absolute minimum', facecolor='none', edgecolor='royalblue', linewidth=2)

                pmap.scatter(xmin, ymin, marker='.', s=5, c='blue', zorder=2)

            elif np.float(min_result['Ang']) > 0:

                pmap.scatter(xmin, ymin, marker='o', s=abs(np.float(min_result['Ang'])) * point_multiplier, c='blue', zorder=2,
                             label='Absolute minimum')

                pmap.scatter(xmin, ymin, marker='.', s=5, c='k', zorder=2)


            # Plot EarthByte ref result
            if np.float(ref_rot_angle) < 0:

                pmap.scatter(xref, yref, marker='o', s=abs(ref_rot_angle) * point_multiplier, c='crimson', zorder=2,
                             label='Reference rotation', facecolor='none', edgecolor='crimson', linewidth=2)

                pmap.scatter(xref, yref, marker='.', s=5, c='crimson', zorder=2)

            elif np.float(ref_result['Ang']) > 0:

                pmap.scatter(xref, yref, marker='o', s=abs(ref_rot_angle) * point_multiplier, c='crimson', zorder=2,
                             label='Reference rotation')

                pmap.scatter(xref, yref, marker='.', s=5, c='k', zorder=2)


            # Plot seed points
            pmap.scatter(slon, slat, marker='.', s=10, c='k', zorder=2)


            plt.legend(scatterpoints=1, labelspacing=2, borderpad=1)


            constraints = ''

            for i in xrange(0, len(data_array)):

                if data_array[i] == True:

                    constraints = constraints + '_' + str(data_array_labels_short[i])

            plt.savefig('model_output/plots/' + model_name + "_" + str(int(ref_rotation_start_age)) + "-" + str(int(ref_rotation_end_age)) + "Ma_" \
            + str(models) + "models" + constraints + "_" + str(geographical_uncertainty) + '.png', format='png', dpi=300)

            plt.show()

            print " "
            print "Reference rotation from " + str(ref_rotation_start_age) + " - " + str(ref_rotation_end_age) + " for plate ID " + str(ref_rotation_plate_id) + "."
            print "Plot produced from " + str(len(xopt)) + " models. Size of circles proportional to angle magnitude."
            print "EarthByte reference model: ", rotation_file
            print " "
            print "Absolute minimum rotation:"
            print display(HTML(min_result.to_html(index=False, columns=['Lon', 'Lat', 'Ang'])))
            print " "
            print "10% mean rotation:"
            print display(HTML(top_ten_mean.to_html(index=False, columns=['Lon', 'Lat', 'Ang'])))
            print " "
            print "Starting ref rotation:"
            print "Lon:", ref_rot_longitude, " - Lat:", ref_rot_latitude, " - Ang:", ref_rot_angle
            # print " "
            # print "Minimised value:", min_result('Minimum')

        return min_result, top_ten_mean




class ObjectiveFunctions():

    # """
    # Hotspot data
    # """
    @staticmethod
    def hotspot_trail(trail_data, ref_rotation_start_age, rotation_model_updated):

        #  Reconstruct data points along chains
        tmp_dist_results = []
        tmp_kappa_results = []
        opt_stats = []
        trail_num = []
        trail_num = []

        for i in xrange(0, len(trail_data)):

            opt_results = []

            # Loop through individual trail point data
            for j in xrange(0, len(trail_data[i]['Lon'])):

                if trail_data[i]['Age'][j] <= ref_rotation_start_age:

                    # Reconstruct each observation data point back to 0Ma relative to proposed Africa location
                    tmp_latLonPoint = pgp.PointOnSphere(trail_data[i]['Lat'][j], trail_data[i]['Lon'][j])

                    # Get rotation for data point and reconstruct to 0Ma
                    point_rotation = rotation_model_updated.get_rotation(np.double(trail_data[i]['Age'][j]), trail_data[i]['PlateID'][j], 0., 1)

                    reconstructed_point = point_rotation * tmp_latLonPoint
                    reconstructed_point_degrees = reconstructed_point.to_lat_lon()

                    opt_results.append({'Chain': trail_data[i]['Chain'], 'Lon': reconstructed_point_degrees[1],
                        'Lat': reconstructed_point_degrees[0], 'Age': trail_data[i]['Age'][j]})


            # Calculate spherical distribution statistics
            if len(opt_results) > 1:

                opt_list = []
                trail_num.append(i)

                for k in xrange(0, len(opt_results)):

                    opt_list.append([opt_results[k]['Lon'], opt_results[k]['Lat'], 1.])

                opt_stats.append(pmag.fisher_mean(opt_list))

            else:

                pass


        #Calculate great circle distance from spherical distribution mean to known hotspot location
        tmp_dist = []
        tmp_kappa = []

        for l in xrange(0, len(opt_stats)):

            tmp_dist.append(geoTools.haversine(opt_stats[l]['dec'], opt_stats[l]['inc'],
                                              trail_data[trail_num[l]]['Hotspot_Location'][1],
                                              trail_data[trail_num[l]]['Hotspot_Location'][0]))

            tmp_kappa.append(opt_stats[l]['k'])


        hs_dist_eval = np.mean(tmp_dist) + np.median(tmp_dist)
        #hs_dist_eval = np.median(tmp_dist)
        hs_kappa_eval = (1. / np.median(tmp_kappa))


        return hs_dist_eval, hs_kappa_eval, tmp_dist, tmp_kappa, opt_stats, trail_num, np.mean(tmp_dist), np.median(tmp_dist), np.std(tmp_dist),\
               np.mean(tmp_kappa), np.median(tmp_kappa), np.std(tmp_kappa)




    # """
    # Fracture zones
    # """
    def Get_FZ_Directions(X1, Y1, X2, Y2):

        long1 = np.radians(X1)
        long2 = np.radians(X2)
        lat1 = np.radians(Y1)
        lat2 = np.radians(Y2)

        bearing = np.arctan2(np.sin(long2 - long1) * np.cos(lat2),
                             np.cos(lat1) * np.sin(lat2) - np.sin(lat1) * np.cos(lat2) *\
                             np.cos(long2 - long1))

        bearing = np.degrees(bearing)
        bearing = (bearing + 360) % 360

        return bearing


    @staticmethod
    def fracture_zones(rotation_model, PID, seafloor_ages, Lats, Lons, spreading_directions):

        # Create empty arrays in which to place the rotated coordinates
        BG = np.zeros(np.asarray(PID).shape) * np.nan
        MG = np.zeros(np.asarray(PID).shape) * np.nan

        # The Delta controls the time-averaging used for the velocities.
        # Note that velocities are always computed assuming that the given time is the end time of the Delta period,
        # for example if the reconstruction time is 45 and Delta = 5, the velocity is for the 50-45 Ma.
        Delta = 5.

        # Loop over each Grid Point
        for j in range (0, len(Lats)):

            plateid = PID[j]
            Age = seafloor_ages[j]

            if np.isfinite(plateid):

                if np.isnan(Age):

                    BG[j] = np.nan

                else:

                    #print Age, plateid
                    # Get the stage pole of the plate wrt mantle at the time of the age grid pixel
                    stage_rotation = rotation_model.get_rotation(Age + Delta, int(plateid), Age, 1, int(plateid))

                    # Find the angle between the
                    bearing1 = Get_FZ_Directions(Lons[j],
                                                 Lats[j],
                                                 stage_rotation.get_lat_lon_euler_pole_and_angle_degrees()[1],
                                                 stage_rotation.get_lat_lon_euler_pole_and_angle_degrees()[0])

                    BG[j] = bearing1
                    MG[j] = stage_rotation.get_rotation_distance((Lats[j],Lons[j])) / Delta


        A = np.fmod(BG + 90., 180.)
        B = np.fmod(spreading_directions, 180.)
        Bearing_Diff = np.remainder(A - B, 180.)
        Bearing_Diff = 90. - np.abs(90. - Bearing_Diff)

        median_skew = np.median(Bearing_Diff)
        mean_skew = np.mean(Bearing_Diff)
        skew_sd = np.std(Bearing_Diff)
        num_segments = len(Bearing_Diff)

        return median_skew, mean_skew, skew_sd, num_segments


    @staticmethod
    def hotspot_trail_misfit(trail_data, ref_rotation_start_age, rotation_model, use_trail_age_uncertainty, trail_age_uncertainty):

        global_point_distance_misfit = []
        global_trail_distance_misfit = []
        global_trail_uncertainty = []
        trail_name = []

        for i in xrange(0, len(trail_data)):

            if trail_data[i]['Age'][-1] >= ref_rotation_start_age - 10:

                opt_results = []

                # Loop through individual trail point data
                for j in xrange(0, len(trail_data[i]['Lon'])):

                    if trail_data[i]['Age'][j] <= ref_rotation_start_age and trail_data[i]['Age'][j] >= ref_rotation_start_age - 10:

                        opt_results.append({'Chain': trail_data[i]['Chain'], 'Lon': trail_data[i]['Lon'][j], 'Lat': trail_data[i]['Lat'][j],
                                            'Age': trail_data[i]['Age'][j], 'Hotspot_Location': trail_data[i]['Hotspot_Location'],
                                            'PlateID': trail_data[i]['PlateID'][j]})

                # If we have no results then skip the current hotspot.
                # This happens if none of the current hotspot's sample ages happen to fall with our 10My interval.
                if not opt_results:
                    continue

                chain = trail_data[i]['Chain']

                tmp_lat = []
                tmp_lon = []
                tmp_age = []
                tmp_hs = []
                tmp_pid = []

                for i in xrange(0, len(opt_results)):

                    tmp_lat.append(opt_results[i]['Lat'])
                    tmp_lon.append(opt_results[i]['Lon'])
                    tmp_age.append(opt_results[i]['Age'])
                    tmp_pid.append(opt_results[i]['PlateID'])
                    tmp_hs.append(opt_results[i]['Hotspot_Location'])

            #print tmp_age, tmp_lat, tmp_lon, tmp_hs, tmp_pid

                # Create synthetic hotspot track
                relativePlate = tmp_pid[0]
                times = np.array(tmp_age)
                times = np.insert(times, 0, 0)
                reconstruction_time = times[0]
                terrane_plateID = 0

                # print 'Chain', chain
                # print 'relative plate', relativePlate
                # print 'times', times
                # print 'hs', tmp_hs[0]
                # print ""

                seed_points_at_digitisation_time = pgp.PointOnSphere(tmp_hs[0])

                motion_path_feature = pgp.Feature.create_motion_path(seed_points_at_digitisation_time,
                                                                     times,
                                                                     valid_time = (times[-1], times[0]),
                                                                     relative_plate = relativePlate,
                                                                     reconstruction_plate_id = terrane_plateID)

                # Reconstruct motion paths
                motion_path_lons, motion_path_lats = [], []

                reconstructed_motion_paths = []

                pgp.reconstruct(motion_path_feature, rotation_model, reconstructed_motion_paths,
                                reconstruction_time, reconstruct_type = pgp.ReconstructType.motion_path)

                for reconstructed_motion_path in reconstructed_motion_paths:

                    path = np.empty([0,0])

                    for point in reconstructed_motion_path.get_motion_path():

                        path = np.append(path, point.to_lat_lon_array())

                    path_lon = np.flipud(path[1::2])
                    path_lat = np.flipud(path[::2])

                    if chain == 'Hawaii' or chain == 'Samoa' or chain == 'Louisville':

                        lon = np.array(path_lon)

                        lon_shift = np.asarray(lon)
                        lon_shift = np.where(lon_shift > 180+180, lon_shift-360 ,lon_shift)
                        lon_shift = np.where(lon_shift < 180-180, lon_shift+360 ,lon_shift)
                        itemindex = lon.shape[0]-np.where(lon_shift[0:-1]-lon_shift[1:] >= 180)[0]

                        motion_path_lons.append(lon_shift - 360)
                        motion_path_lats.append(path_lat)

                    else:

                        motion_path_lons.append(path_lon)
                        motion_path_lats.append(path_lat)

                # Append hotspot location to interpolated data
                tmp_lon = np.insert(np.array(tmp_lon), 0, tmp_hs[0][1])
                tmp_lat = np.insert(np.array(tmp_lat), 0, tmp_hs[0][0])

                # print tmp_lat
                # print tmp_lon

                point_distance = []
                ref_trail_distance = []
                ref_trail_azimuth = []
                new_trail_distance = []
                new_trail_azimuth = []
                trail_distance_misfit = []
                point_distance_misfit = []
                uncertainty = []

                trail_name.append(chain)

                # Iterate along motion paths and extract distance metrics
                if use_trail_age_uncertainty == False:

                    for i in xrange(1, len(tmp_lon)):

                        point_distance_tmp = geoTools.haversine(tmp_lon[i], tmp_lat[i], motion_path_lons[0][i], motion_path_lats[0][i])
                        point_distance.append(point_distance_tmp[0])

                        ref_trail_distance_tmp = geoTools.haversine(tmp_lon[i], tmp_lat[i], tmp_lon[i - 1], tmp_lat[i - 1])
                        ref_trail_distance.append(ref_trail_distance_tmp[0])
                        ref_trail_azimuth.append(ref_trail_distance_tmp[1])

                        new_trail_distance_tmp = geoTools.haversine(motion_path_lons[0][i], motion_path_lats[0][i], tmp_lon[i - 1], tmp_lat[i - 1])
                        new_trail_distance.append(new_trail_distance_tmp[0])
                        new_trail_azimuth.append(new_trail_distance_tmp[1])

                        trail_distance_misfit.append(abs(ref_trail_distance_tmp[0] - new_trail_distance_tmp[0]))
                        point_distance_misfit.append(abs(point_distance_tmp[0]))

                    global_point_distance_misfit.append(np.sum(point_distance_misfit))
                    global_trail_distance_misfit.append(np.sum(trail_distance_misfit))

                elif use_trail_age_uncertainty == True:

                    for i in xrange(1, len(tmp_lon)):

                        point_distance_tmp = geoTools.haversine(tmp_lon[i], tmp_lat[i], motion_path_lons[0][i], motion_path_lats[0][i])
                        point_distance.append(point_distance_tmp[0])

                        ref_trail_distance_tmp = geoTools.haversine(tmp_lon[i], tmp_lat[i], tmp_lon[i - 1], tmp_lat[i - 1])
                        ref_trail_distance.append(ref_trail_distance_tmp[0])
                        ref_trail_azimuth.append(ref_trail_distance_tmp[1])

                        seg_age = times[i] - times[i - 1]
                        seg_uncert = (ref_trail_distance_tmp[0] / seg_age) * trail_age_uncertainty
                        uncertainty.append(seg_uncert)

                        new_trail_distance_tmp = geoTools.haversine(motion_path_lons[0][i], motion_path_lats[0][i], tmp_lon[i - 1], tmp_lat[i - 1])
                        new_trail_distance.append(new_trail_distance_tmp[0])
                        new_trail_azimuth.append(new_trail_distance_tmp[1])

                        trail_distance_misfit.append(abs(ref_trail_distance_tmp[0] - new_trail_distance_tmp[0]))
                        point_distance_misfit.append(abs(point_distance_tmp[0]))

                    global_point_distance_misfit.append(np.sum(point_distance_misfit))
                    global_trail_distance_misfit.append(np.sum(trail_distance_misfit))
                    global_trail_uncertainty.append(np.mean(uncertainty))


        return global_point_distance_misfit, global_trail_distance_misfit, global_trail_uncertainty, trail_name
