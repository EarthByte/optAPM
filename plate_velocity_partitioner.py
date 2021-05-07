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
            plate_features,
            plate_features_are_topologies,
            # Used to generate contour polygons when aggregrating blocks of continental polygons.
            # Only applies when plate features are continents (should be None if plate_features_are_topologies is True).
            continent_fragmentation,
            data_model,
            grid_spacing_degrees):
        """
        Load the topology features and original rotation model, and find plate IDs at a grid of point positions.
        """
        
        self.data_dir = data_dir
        self.plate_features = plate_features
        self.plate_features_are_topologies = plate_features_are_topologies
        self.continent_fragmentation = continent_fragmentation
        
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
            file_to_remove = os.path.join(self.data_dir, self.plate_velocity_filename)
            if os.path.exists(file_to_remove):
                os.remove(file_to_remove)
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
        Generate a multipoint for each plate (and associated plate ID) and save them to a file.
        """
        
        output_features = []

        #
        # Resolve topological plate polygons or reconstruct continental polygons.
        #
        if self.plate_features_are_topologies:
            # Resolve the topological plate polygons for the current time.
            resolved_topologies = []
            pygplates.resolve_topologies(self.plate_features, self.rotation_model, resolved_topologies, ref_rotation_start_age,
                # TODO: We are only looking at *non-deforming* boundaries, but need to also support deforming regions
                #       Requires extensions to pyGPlates for that.
                resolve_topology_types=pygplates.ResolveTopologyType.boundary)
            
            # Get a list of resolved polygons and a list of their plate IDs.
            polygons = []
            polygon_plate_ids = []
            for resolved_topology in resolved_topologies:
                polygons.append(resolved_topology.get_resolved_boundary())
                polygon_plate_ids.append(resolved_topology.get_feature().get_reconstruction_plate_id())
            
            # Find the resolved plate polygon (if any) containing each point.
            point_polygon_indices = points_in_polygons.find_polygons(
                    self.points,
                    # The resolved plate polygon geometries...
                    polygons,
                    # The index of each resolved plate polygon (this is what is returned by 'find_polygons')...
                    list(range(len(polygons))))
            
            # Dictionary of point lists indexed by plate ID (ie, each plate ID has a list of associated points).
            points_by_plate_id = {}
            
            # Each point is contained by one resolved plate boundary.
            for point_index, polygon_index in enumerate(point_polygon_indices):
                # If point is not in any resolved boundary then it either fell in a tiny crack/gap or the topologies don't have global coverage.
                # In which case we ignore the point.
                if polygon_index is not None:
                    point_plate_id = polygon_plate_ids[polygon_index]
                    # Create dictionary entry (an empty list) if first time encountered plate ID.
                    if point_plate_id not in points_by_plate_id:
                        points_by_plate_id[point_plate_id] = (polygon_index, [])
                    # Add point to dictionary entry.
                    points_by_plate_id[point_plate_id][1].append(self.points[point_index])
            
            # Generate a multi-point feature for each plate ID (containing all points associated with each plate ID).
            for plate_id, (polygon_index, points_with_plate_id) in points_by_plate_id.items():
                # Create the multi-point geometry.
                multi_point = pygplates.MultiPointOnSphere(points_with_plate_id)
                
                # Create the multi-point feature.
                multi_point_feature = pygplates.Feature()
                multi_point_feature.set_geometry(multi_point)
                multi_point_feature.set_reconstruction_plate_id(plate_id)
                
                # Set some contour parameters for later (using shapefile attributes).
                polygon = polygons[polygon_index]
                multi_point_feature.set_shapefile_attributes({
                        'contour_perimeter' : polygon.get_arc_length(),
                        'contour_area' : polygon.get_area()})
                
                output_features.append(multi_point_feature)
            
        else:  # Continental polygons...
            # Reconstruct the continental polygons.
            reconstructed_feature_geometries = []
            pygplates.reconstruct(self.plate_features, self.rotation_model, reconstructed_feature_geometries, ref_rotation_start_age)
            
            # Get a list of reconstructed continent polygons and a list of their plate IDs.
            reconstructed_continent_polygons = []
            continent_polygon_plate_ids = []
            for reconstructed_feature_geometry in reconstructed_feature_geometries:
                reconstructed_continent_polygons.append(reconstructed_feature_geometry.get_reconstructed_geometry())
                continent_polygon_plate_ids.append(reconstructed_feature_geometry.get_feature().get_reconstruction_plate_id())
            
            # Find the reconstructed continental polygon (if any) containing each point.
            continent_point_plate_ids = points_in_polygons.find_polygons(
                    self.points,
                    # The reconstructed continental polygons...
                    reconstructed_continent_polygons,
                    # The plate ID of each reconstructed continent (this is what is returned by 'find_polygons')...
                    continent_polygon_plate_ids)
            
            # Each point is contained by one reconstructed polygon (or none).
            continent_points = []
            continent_point_plate_ids = []
            for point_index, point_plate_id in enumerate(continent_point_plate_ids):
                # If point is not in any continental polygon then it is on oceanic crust and we ignore it.
                # If it's in two continental polygons then they overlap (which can happen when continental polygons are reconstructed),
                # in which case we asked for only one of the polygons to be kept by 'points_in_polygons.find_polygons()'.
                if point_plate_id is not None:
                    continent_points.append(self.points[point_index])
                    continent_point_plate_ids.append(point_plate_id)
            
            # Calculate contour polygons representing the boundary(s) of the reconstructed continental polygons that overlap each other.
            reconstructed_contour_polygons = self.continent_fragmentation.get_contour_polygons(reconstructed_continent_polygons)
            
            # Find the reconstructed contour polygon (if any) containing each continent point.
            continent_point_contour_polygon_indices = points_in_polygons.find_polygons(continent_points,
                    # The reconstructed contour polygons...
                    reconstructed_contour_polygons,
                    # The index of each reconstructed contour polygon (this is what is returned by 'find_polygons')...
                    list(range(len(reconstructed_contour_polygons))))
            
            # Dictionary of continent points indexed by contour polygon (use a contour polygon *index* since pygplates.PolygonOnSphere is not hashable).
            # Each group of continent points is itself a dict mapping continent plate ID to a list of continent points with that plate ID
            # (because a single contour polygon could encompass continental polygons with different plate IDs).
            continent_points_by_contour_polygon = {}
            
            # Each continent point is contained by one reconstructed contour polygon (or none).
            for continent_point_index, contour_polygon_index in enumerate(continent_point_contour_polygon_indices):
                # If point is not in any contour polygons then ignore it.
                #
                # This mainly happens when the area of a contour polygon is below the minimum and hence excludes
                # the continental polygons that that contour encompassed. Aside from that, if any other continent
                # points are outside any contour polygons then it means a contour polygon did not fully encompass
                # a set of continental polygons. But there shouldn't be many of those points.
                if contour_polygon_index is None:
                    continue
                
                # Create dictionary entry (an empty dict) if first time encountered contour polygon.
                if contour_polygon_index not in continent_points_by_contour_polygon:
                    continent_points_by_contour_polygon[contour_polygon_index] = {}
                continent_points_in_contour_polygon = continent_points_by_contour_polygon[contour_polygon_index]
                
                # The current continent point and its plate ID.
                continent_point = continent_points[continent_point_index]
                continent_point_plate_id = continent_point_plate_ids[continent_point_index]
                
                # Create dictionary entry (an empty list) if first time encountered continent plate ID.
                if continent_point_plate_id not in continent_points_in_contour_polygon:
                    continent_points_in_contour_polygon[continent_point_plate_id] = []

                # Add current continent point and its plate ID to dictionary entry.
                continent_points_in_contour_polygon[continent_point_plate_id].append(continent_point)
            
            # Scalar type used to store plate IDs in each contour coverage geometry.
            plate_id_scalar_type = pygplates.ScalarType.create_gpml('PlateID')

            # Generate a multi-point feature for each contour polygon (containing all points associated with contour polygon).
            for contour_polygon_index, continent_points_in_contour_polygon in continent_points_by_contour_polygon.items():
                
                # Create the multi-point feature.
                coverage_multi_point_feature = pygplates.Feature()

                # Accumulate the coverage points and associated plate IDs contained within the current contour polygon.
                # Note that points are grouped by plate ID (ie, first N points have plate_id_1, next M points have plate_id_2, etc).
                coverage_points = []
                coverage_plate_ids = []
                for continent_plate_id, points_in_continent in continent_points_in_contour_polygon.items():
                    coverage_points.extend(points_in_continent)
                    # Each point in the current continent has the same plate ID.
                    coverage_plate_ids.extend([continent_plate_id] * len(points_in_continent))

                # Create the multi-point coverage (geometry + plate IDs).
                coverage_multipoint = pygplates.MultiPointOnSphere(coverage_points)
                coverage_scalars = { plate_id_scalar_type : coverage_plate_ids }
                coverage_multi_point_feature.set_geometry((coverage_multipoint, coverage_scalars))
                
                # Set some contour parameters for later (using shapefile attributes).
                contour_polygon = reconstructed_contour_polygons[contour_polygon_index]
                coverage_multi_point_feature.set_shapefile_attributes({
                        # Mark the feature as a continent aggregate (for later in objective function).
                        # Any value will do (the presence of the attribute is all that's checked)...
                        'is_in_continent_contour' : 1,
                        'contour_perimeter' : contour_polygon.get_arc_length(),
                        'contour_area' : contour_polygon.get_area()})

                output_features.append(coverage_multi_point_feature)
        
        # Write points and their plate IDs to the plate velocities file.
        pygplates.FeatureCollection(output_features).write(
                os.path.join(self.data_dir, self.plate_velocity_filename))
