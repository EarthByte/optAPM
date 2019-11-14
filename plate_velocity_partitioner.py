import math
import os
import os.path
import points_in_polygons
import pygplates
import sys


class PlateVelocityPartitioner(object):
    """
    Class to partition a grid of points into resolved topologies or continental polygons at each reconstruction time.
    This is the bulk of the computation and enables the objective function to efficiently calculate velocities
    (from the positions and associated plate IDs we calculate) in its many iterations over candidate rotation models.
    """
    
    
    def __init__(
            self,
            data_dir,
            original_rotation_filenames,  # Relative to the 'data/' directory.
            topology_features,
            data_model,
            grid_spacing_degrees):
        """
        Load the topology features and original rotation model, and find plate IDs at a grid of point positions.
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
        
        # Load all the rotations into a rotation model.
        self.rotation_model = pygplates.RotationModel(rotation_features)
        
        # The plate velocity filename (relative to the 'data/' directory).
        self.plate_velocity_filename = os.path.join(
                data_model, 'optimisation', 'temp_plate_velocities.gpml')
        
        #
        # Create a multi-point uniformly lat-lon spaced.
        #
        self.points = []
        
        num_latitudes = int(math.floor(180.0 / grid_spacing_degrees))
        num_longitudes = int(math.floor(360.0 / grid_spacing_degrees))
        for lat_index in range(num_latitudes):
            # The 0.5 puts the point in the centre of the grid pixel.
            # This also avoids sampling right on the poles.
            lat = -90 + (lat_index + 0.5) * grid_spacing_degrees
            
            for lon_index in range(num_longitudes):
                # The 0.5 puts the point in the centre of the grid pixel.
                # This also avoids sampling right on the dateline.
                lon = -180 + (lon_index + 0.5) * grid_spacing_degrees
                
                self.points.append((lat, lon))
        
        # Convert a list of (lat,lon) to a pygplates.MultiPointOnSphere
        # (which behaves like an iterable sequence of pygplates.PointOnSphere).
        self.points = pygplates.MultiPointOnSphere(self.points)
    
    
    def __del__(self):
        # Remove temporary plate velocity file.
        try:
            os.remove(
                os.path.join(self.data_dir, self.plate_velocity_filename))
        except AttributeError:
            # 'self.plate_velocity_filename' might not exist if exception raised inside '__init__()'.
            pass
    
    
    def get_plate_velocity_filename(self):
        """
        Return the filename (relative to the 'data/' directory) of the file containing a grid of points
        masked by resolved topologies (or continental polygons) at the reconstruction time specified in
        most recent call to 'generate_points_and_plate_ids()'.
        """
        return self.plate_velocity_filename
    
    
    def generate_points_and_plate_ids(
            self,
            ref_rotation_start_age):
        """
        Generate the resolved trench features at the specified time and save them to the trench migration file.
        """
        
        # Resolve the plate polygons for the current time.
        resolved_topologies = []
        pygplates.resolve_topologies(self.topology_features, self.rotation_model, resolved_topologies, ref_rotation_start_age,
            # TODO: We are only looking at *non-deforming* boundaries, but need to also support deforming regions
            #       Requires extensions to pyGPlates for that.
            resolve_topology_types=pygplates.ResolveTopologyType.boundary)
        
        #
        # Assuming the resolved polygons are *non-overlapping*, which they should be (except for minor overap/gap errors),
        # find the single polygon (resolved boundary) containing each point.
        #
        point_plate_ids = points_in_polygons.find_polygons(
                self.points,
                # The resolved plate polygon geometries...
                [resolved_topology.get_resolved_boundary() for resolved_topology in resolved_topologies],
                # The plate ID of each resolved topology (this is what is returned by 'find_polygons')...
                [resolved_topology.get_feature().get_reconstruction_plate_id() for resolved_topology in resolved_topologies])
        
        # Dictionary of point lists indexed by plate ID (ie, each plate ID has a list of associated points).
        points_by_plate_id_dict = {}
        
        # Each point is contained by one resolved plate boundary.
        for point_index, point_plate_id in enumerate(point_plate_ids):
            # If point is not in any resolved boundary then it either fell in a tiny crack/gap
            # or the topologies don't have global coverage.
            # In either case we ignore it.
            if point_plate_id is not None:
                # Create dictionary entry (an empty list) if first time encountered plate ID.
                if point_plate_id not in points_by_plate_id_dict:
                    points_by_plate_id_dict[point_plate_id] = []
                # Add point to dictionary entry.
                points_by_plate_id_dict[point_plate_id].append(self.points[point_index])
        
        # Generate a multi-point feature for each plate ID (containing all points associated with each plate ID).
        partitioned_multi_point_features = []
        for plate_id in points_by_plate_id_dict.keys():
            # Create the multi-point geometry.
            points_with_plate_id = points_by_plate_id_dict[plate_id]
            multi_point = pygplates.MultiPointOnSphere(points_with_plate_id)
            
            # Create the multi-point feature.
            multi_point_feature = pygplates.Feature()
            multi_point_feature.set_geometry(multi_point)
            multi_point_feature.set_reconstruction_plate_id(plate_id)
            
            partitioned_multi_point_features.append(multi_point_feature)
        
        # Write points and their plate IDs to the plate velocities file.
        pygplates.FeatureCollection(partitioned_multi_point_features).write(
                os.path.join(self.data_dir, self.plate_velocity_filename))
