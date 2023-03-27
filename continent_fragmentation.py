from  collections import deque
import math
import os
import os.path
from ptt.continent_contours import ContouredContinent, ContinentContouring
import pygplates
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
            # Function (accepting time in Ma) to return area threshold (in square radians) when creating continent contours...
            continent_contouring_area_threshold_steradians_function,
            # Function (accepting time in Ma) to return distance threshold to ensure small gaps between continents are ignored during contouring...
            continent_contouring_gap_threshold_radians_function,
            age_range):
        """
        Load the continent features and use *original* rotation model to calculate fragment through all time to find normalisation factor.
        """
        
        # Load all the original rotation feature collections.
        rotation_features = []
        for rotation_filename in original_rotation_filenames:
            # Read the current rotation file.
            rotation_feature_collection = pygplates.FeatureCollection(
                    os.path.join(data_dir, rotation_filename))
            rotation_features.extend(rotation_feature_collection)
        
        self.continent_contouring = ContinentContouring(
                rotation_features,
                continent_features,
                continent_contouring_point_spacing_degrees,
                continent_contouring_area_threshold_steradians_function,
                continent_contouring_gap_threshold_radians_function)

        self.age_range = age_range
        
        # A dict mapping age to global fragmentation index.
        # This gets populated as needed.
        self.fragmentations = {}

        # We'll calculate this later if a *normalized* fragmentation index is ever requested since
        # that'll require calculating fragmentation indices over the entire age range.
        self.max_fragmentation = None
    
    
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
            self.fragmentations[age] = self.continent_contouring.get_fragmentation(age)
        
        fragmentation = self.fragmentations[age]

        if normalize:
            # First calculate the maximum fragmentation index (if we haven't already).
            if not self.max_fragmentation:
                print('Calculating continental fragmentation for {0}-{1}Ma...'.format(self.age_range[0], self.age_range[-1]))
                sys.stdout.flush()

                # Find the maximum fragmentation over the age range.
                # This also caches the fragmentation values over the age range.
                self.max_fragmentation = max(self.continent_contouring.get_fragmentation(age) for age in self.age_range)

                #print('  min:', np.min(list(self.fragmentations.values())) / self.max_fragmentation)
                #print('  mean:', np.mean(list(self.fragmentations.values())) / self.max_fragmentation)
                #print('  dev:', np.std(list(self.fragmentations.values())) / self.max_fragmentation)
                #sys.stdout.flush()
            
            # A *normalised* version of the pre-calculated fragmentation at age.
            fragmentation /= self.max_fragmentation

        return fragmentation
    
    
    def get_contoured_continents(
            self,
            age):
        """
        Reconstruct the continents and then find the boundaries of the specified (potentially overlapping/abutting) continent polygons.
        """

        # Calculate contoured continents representing the boundary(s) of the reconstructed continent polygons that overlap each other.
        return self.continent_contouring.get_contoured_continents(age)

    
    def calculate_contoured_continents(
            self,
            continent_polygons,
            age):
        """
        Find the boundaries of the specified (potentially overlapping/abutting) continent polygons as contour polygons.

        The 'age' is only used to look up the time-dependent contouring thresholds (using functions passed into constructor).
        """

        return self.continent_contouring.calculate_contoured_continents(continent_polygons, age)
