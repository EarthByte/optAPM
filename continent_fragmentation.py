import math
import os
import os.path
import points_in_polygons
import proximity_query
import pygplates
from skimage import measure
import sys
import numpy as np


class ContinentFragmentation(object):
    """
    Class to calculate continental fragmentation (global perimeter-to-area ratio).
    """
    
    
    def __init__(
            self,
            data_dir,
            original_rotation_filenames,  # Relative to the 'data/' directory.
            continent_features,
            # Point spacing used when creating continent contours...
            continent_contouring_point_spacing_degrees,
            # Area threshold (in square radians) when creating continent contours...
            continent_contouring_area_threshold_steradians,
            # Distance threshold to ensure small gaps between continents are ignored during contouring...
            continent_contouring_gap_threshold_radians,
            age_range):
        """
        Load the continent features and use *original* rotation model to calculate fragment through all time to find normalisation factor.
        """
        
        self.data_dir = data_dir
        self.continent_features = continent_features

        # A point grid to calculate contour polygons representing the boundary of reconstructed static polygons that overlap each other.
        self.contouring_area_threshold_steradians = continent_contouring_area_threshold_steradians
        self.continent_contouring_gap_threshold_radians = continent_contouring_gap_threshold_radians
        self.contouring_point_spacing_degrees = continent_contouring_point_spacing_degrees
        lons = np.arange(-180.0, 180.001, self.contouring_point_spacing_degrees)
        lats = np.arange(-90.0, 90.001, self.contouring_point_spacing_degrees)
        self.contouring_grid_dimensions = len(lats), len(lons)

        self.age_range = age_range

        contouring_longitude_array, contouring_latitude_array = np.meshgrid(lons, lats)
        self.contouring_points = pygplates.MultiPointOnSphere(
                zip(contouring_latitude_array.flatten(), contouring_longitude_array.flatten()))
        
        self.debug_contour_polygons = False
        if self.debug_contour_polygons:
            self.debug_contour_polygon_features = []
            self.debug_time_interval = age_range[1] - age_range[0]
        
        # Load all the original rotation feature collections.
        rotation_features = []
        for rotation_filename in original_rotation_filenames:
            # Read the current rotation file.
            rotation_feature_collection = pygplates.FeatureCollection(
                    os.path.join(self.data_dir, rotation_filename))
            rotation_features.extend(rotation_feature_collection)
        
        # Load all original rotations into a rotation model.
        self.rotation_model = pygplates.RotationModel(rotation_features)
        
        # A dict mapping age to global fragmentation index.
        # This gets populated as needed.
        self.fragmentations = {}

        # We'll calculate this later if a *normalized* fragmentation index is ever requested since
        # that'll require calculating fragmentation indices over the entire age range.
        self.max_fragmentation = None

        # Debug output contour polygons to GPML.
        if self.debug_contour_polygons:
            # Easiest way to force contour polygons to be generated for all times is to get a *normalized* fragmentation index.
            self.get_fragmentation(age_range[0], normalize=True)
            pygplates.FeatureCollection(self.debug_contour_polygon_features).write('contour_polygons.gpmlz')
    
    
    def get_fragmentation(
            self,
            age,
            normalize=True):
        """
        Calculate the continental fragmentation index at the specified time.
        
        If *normalize* is True then the normalisation factor will be such that approx 1.0 will represent the
        maximum fragmentation over the age range specified in constructor.
        """

        # Calculate and cache result in dict (if not already cached).
        if age not in self.fragmentations:
            self.fragmentations[age] = self._calculate_fragmentation(age)
        
        fragmentation = self.fragmentations[age]

        if normalize:
            # First calculate the maximum fragmentation index (if we haven't already).
            if not self.max_fragmentation:
                print('Calculating continental fragmentation for {0}-{1}Ma...'.format(self.age_range[0], self.age_range[-1]))
                sys.stdout.flush()

                # Find the maximum fragmentation over the age range.
                # This also caches the fragmentation values over the age range.
                self.max_fragmentation = max(self.get_fragmentation(age, normalize=False) for age in self.age_range)

                #print('  min:', np.min(list(self.fragmentations.values())) / self.max_fragmentation)
                #print('  mean:', np.mean(list(self.fragmentations.values())) / self.max_fragmentation)
                #print('  dev:', np.std(list(self.fragmentations.values())) / self.max_fragmentation)
                #sys.stdout.flush()
            
            # A *normalised* version of the pre-calculated fragmentation at age.
            fragmentation /= self.max_fragmentation

        return fragmentation
    
    
    def _calculate_fragmentation(
            self,
            age):

        total_perimeter = 0.0
        total_area = 0.0
        
        # Reconstruct static continental polygons.
        reconstructed_feature_geometries = []
        pygplates.reconstruct(self.continent_features, self.rotation_model, reconstructed_feature_geometries, age)
        
        # Get a list of polygons.
        #
        # We should have polygons (not polylines) but turn into a polygon if happens to be a polyline
        # (but that actually only works if the polyline is a closed loop and not just part of a polygon's boundary).
        reconstructed_polygons = [pygplates.PolygonOnSphere(reconstructed_feature_geometry.get_reconstructed_geometry())
                for reconstructed_feature_geometry in reconstructed_feature_geometries]
        
        # Calculate contour polygons representing the boundary(s) of the reconstructed static polygons that overlap each other.
        reconstructed_contour_polygons = self.get_contour_polygons(reconstructed_polygons)

        # Update total perimeter and area.
        for contour_polygon in reconstructed_contour_polygons:
            total_perimeter += contour_polygon.get_arc_length()
            total_area += contour_polygon.get_area()

        # Debug output contour polygons.
        if self.debug_contour_polygons:
            self.debug_contour_polygon_features.extend(
                    pygplates.Feature.create_reconstructable_feature(
                            pygplates.FeatureType.gpml_unclassified_feature,
                            reconstructed_contour_polygon,
                            valid_time=(age + 0.5 * self.debug_time_interval, age - 0.5 * self.debug_time_interval))
                    for reconstructed_contour_polygon in reconstructed_contour_polygons)

        #print('age:', age, 'frag_index (1/km):', total_perimeter / total_area / 6371.0); sys.stdout.flush()
        return total_perimeter / total_area
    
    
    def get_contour_polygons(
            self,
            continent_polygons):
        """
        Find the boundaries of the specified (potentially overlapping/abutting) continent polygons as contour polygons.
        """

        contour_polygons = self._calculate_contour_polygons(continent_polygons)
        
        # Contour polygons smaller than this will be excluded.
        min_area = self.contouring_area_threshold_steradians
        # A contour polygon's area should not be more than half the global area.
        # It seems this can happen with pygplates revisions prior to 31 when there's a sliver polygon along the dateline
        # (that gets an area that's a multiple of PI, instead of zero).
        max_area = 2 * math.pi - 1e-4

        contour_polygons_above_area_threshold = []
        for contour_polygon in contour_polygons:
            # Exclude contour polygon if smaller than the threshold.
            contour_polygon_area = contour_polygon.get_area()
            if (contour_polygon_area > min_area and
                contour_polygon_area < max_area):
                contour_polygons_above_area_threshold.append(contour_polygon)

        return contour_polygons_above_area_threshold

    
    def _calculate_contour_polygons(
            self,
            continent_polygons):
        """
        Find the boundaries of the specified (potentially overlapping/abutting) continent polygons as contour polygons.

        Note that this function is based on code written by Andrew Merdith and Simon Williams:
        See https://github.com/amer7632/pyGPlates_examples/blob/master/Merdith_2019_GPC/Perimeter-to-area-ratio.ipynb
        """

        # Find the reconstructed continental polygon (if any) containing each point.
        continent_polygons_containing_points = points_in_polygons.find_polygons(
                self.contouring_points,
                continent_polygons)

        # Find the reconstructed continental polygon (if any) near each point.
        # If a point is outside all polygon but close enough to the outline of a polygon then it's considered inside a contour.
        # This ensures small gaps between continents are ignored during contouring.
        continent_polygons_near_points = proximity_query.find_closest_geometries_to_points(
                self.contouring_points,
                continent_polygons,
                distance_threshold_radians=self.continent_contouring_gap_threshold_radians)

        zval = []
        for contouring_point_index in range(len(self.contouring_points)):
            continent_polygon_containing_point = continent_polygons_containing_points[contouring_point_index]
            continent_polygon_near_point = continent_polygons_near_points[contouring_point_index]

            # If the current point is either inside a continent polygon or close enough to its outline
            # then mark the point as inside a contour.
            if (continent_polygon_containing_point is not None or
                continent_polygon_near_point is not None):
                zval.append(1)
            else:
                zval.append(0)
            
        bi = np.array(zval).reshape(self.contouring_grid_dimensions[0], self.contouring_grid_dimensions[1])
    
        # To handle edge effects, pad grid before making contour polygons.
        pad_hor = np.zeros((1, bi.shape[1]))
        pad_ver = np.zeros((bi.shape[0]+1, 1))
        pad1 = np.vstack((bi, pad_hor))
        pad2 = np.hstack((pad_ver, pad1))
        pad3 = np.hstack((pad2, pad_ver))
        contours = measure.find_contours(pad3, 0.5, fully_connected='low')

        contour_polygons = []
        for contour in contours:
            # To handle edge effects again - strip off parts of polygon
            # due to padding, and adjust from image coordinates to long/lat
            contour[:,1] = (contour[:,1] * self.contouring_point_spacing_degrees) - 1
            contour[:,0] = (contour[:,0] * self.contouring_point_spacing_degrees) - 1
            contour[np.where(contour[:,0] < 0.0), 0] = 0
            contour[np.where(contour[:,0] > 180.0), 0] = 180
            contour[np.where(contour[:,1] < 0.0), 1] = 0
            contour[np.where(contour[:,1] > 360.0), 1] = 360
            
            contour_polygon = pygplates.PolygonOnSphere(zip(contour[:,0] - 90, contour[:,1] - 180))

            contour_polygons.append(contour_polygon)

        return contour_polygons
