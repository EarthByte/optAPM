import os
import os.path
import pygplates
import sys


class ContinentFragmentation(object):
    """
    Class to calculate continental fragmentation (global perimeter-to-area ratio).
    """
    
    
    def __init__(
            self,
            data_dir,
            original_rotation_filenames,  # Relative to the 'data/' directory.
            topology_features,
            age_range):
        """
        Load the topology features and use *original* rotation model to calculate fragment through all time to find normalisation factor.
        
        The normalisation factor should be such that approx 1.0 will represent the maximum fragmentation when later using optimised
        rotation model to calculate fragmentation at each time interval.
        """
        
        self.data_dir = data_dir
        self.topology_features = topology_features
        
        # Load all the original rotation feature collections.
        rotation_features = []
        for rotation_filename in original_rotation_filenames:
            # Read the current rotation file.
            rotation_feature_collection = pygplates.FeatureCollection(
                    os.path.join(self.data_dir, rotation_filename))
            rotation_features.extend(rotation_feature_collection)
        
        # Load all original rotations into a rotation model.
        self.rotation_model = pygplates.RotationModel(rotation_features)
        
        # Find the maximum fragmentation over the age range.
        print 'Calculating continental fragmentation for {0}-{1}Ma...'.format(age_range[0], age_range[-1])
        sys.stdout.flush()
        self.max_fragmentation = max(self._calculate_fragmentation(age) for age in age_range)
    
    
    def get_fragmentation(
            self,
            age):
        """
        Calculate the normalised continental fragmentation index at the specified time and save them to the trench migration file.
        """

        # Return the *normalised* fragmentation.
        return self._calculate_fragmentation(age) / self.max_fragmentation
    
    
    def _calculate_fragmentation(
            self,
            age):
        
        # Resolve topologies.
        resolved_topologies = []
        pygplates.resolve_topologies(
                self.topology_features,
                self.rotation_model,
                resolved_topologies,
                age)

        total_perimeter = 0.0
        total_area = 0.0

        # Iterate over the resolved topologies.
        for resolved_topology in resolved_topologies:
            total_perimeter += resolved_topology.get_resolved_boundary().get_arc_length()
            total_area += resolved_topology.get_resolved_boundary().get_area()

        return total_perimeter / total_area