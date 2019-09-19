
"""
    Copyright (C) 2017 The University of Sydney, Australia
    
    This program is free software; you can redistribute it and/or modify it under
    the terms of the GNU General Public License, version 2, as published by
    the Free Software Foundation.
    
    This program is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.
    
    You should have received a copy of the GNU General Public License along
    with this program; if not, write to Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
"""


##################################################
# Find the convergence rate of subduction zones. #
#                                                #
# NOTE: This is a modified version of            #
# "subduction_convergence.py" useful for         #
# resolving topologies once (slow) and then      #
# calculating velocities using a rotation model  #
# with various absolute rotation adjustments     #
# applied.                                       #
##################################################


from __future__ import print_function
import argparse
import math
import os.path
import sys
import pygplates
import warnings


# PyGPlates version 22 can handle topological lines (can get their sub-sub-segment plate IDs).
PYGPLATES_VERSION_REQUIRED = pygplates.Version(22)
# Check the imported pygplates version.
if not hasattr(pygplates, 'Version') or pygplates.Version.get_imported_version() < PYGPLATES_VERSION_REQUIRED:
    raise ImportError('{0}: Error - imported pygplates version {1} but version {2} or greater is required'.format(
            os.path.basename(__file__), pygplates.Version.get_imported_version(), PYGPLATES_VERSION_REQUIRED))


# Resolves subduction zone topologies at 'time' and returns resolved geometries as a list of features.
#
# Note that features normally store present-day geometries but the returned features are reconstructed/resolved.

def resolve_subduction_zones(
        # Rotation model or feature collection(s), or list of features, or filename(s)...
        rotation_features_or_model,
        # Topology feature collection(s), or list of features, or filename(s) or any combination of those...
        topology_features,
        time,
        anchor_plate_id = 0):
    
    # Turn rotation data into a RotationModel (if not already).
    rotation_model = pygplates.RotationModel(rotation_features_or_model)
    
    # Turn topology data into a list of features (if not already).
    topology_features = pygplates.FeaturesFunctionArgument(topology_features)
    
    # Resolve our topological plate polygons (and deforming networks) to the current 'time'.
    # We generate both the resolved topology boundaries and the boundary sections between them.
    resolved_topologies = []
    shared_boundary_sections = []
    pygplates.resolve_topologies(topology_features.get_features(), rotation_model, resolved_topologies, time, shared_boundary_sections, anchor_plate_id)
    
    # List of resolved subduction zone features from shared subsegments for the current 'time'.
    resolved_subduction_features = []
    
    # Iterate over the shared boundary sections of all resolved topologies.
    for shared_boundary_section in shared_boundary_sections:
    
        # Skip sections that are not subduction zones.
        if shared_boundary_section.get_feature().get_feature_type() != pygplates.FeatureType.gpml_subduction_zone:
            continue
        
        # Iterate over the shared sub-segments of the current subducting line.
        # These are the parts of the subducting line that actually contribute to topological boundaries.
        for shared_sub_segment in shared_boundary_section.get_shared_sub_segments():
            
            # The plate ID of the trench line (as opposed to the subducting plate).
            #
            # Update: The plate IDs of the trench line and overriding plate can differ
            # even in a non-deforming model due to smaller plates, not modelled by topologies, moving
            # differently than the larger topological plate being modelled - and the trench line
            # having plate IDs of the smaller plates near them. For that reason we use the plate ID
            # of the trench line whenever we can.
            #
            # If the current shared sub-segment is part of a topological line then we obtain its sub-sub-segments
            # (since we require pyGPlates version 22 or above). This is because trench lines that are
            # topological lines might actually be deforming (or intended to be deforming) and hence their
            # plate ID is not meaningful or at least we can't be sure whether it will be zero or the
            # overriding plate (or something else). In this case we look at the plate IDs of the sub-sub-segments.
            #
            sub_segments_of_topological_line_sub_segment = shared_sub_segment.get_sub_segments()
            if sub_segments_of_topological_line_sub_segment:
                
                # Subduction polarity - default to 'Unknown' if subduction zone doesn't have one.
                subduction_polarity = shared_sub_segment.get_feature().get_enumeration(
                    pygplates.PropertyName.gpml_subduction_polarity,
                    'Unknown')
                
                # Iterate over the sub-sub-segments associated with the topological line.
                for sub_sub_segment in sub_segments_of_topological_line_sub_segment:
                    
                    # Create the resolved subduction zone feature.
                    resolved_subduction_feature = sub_sub_segment.get_resolved_feature()
                    
                    # Transfer the subduction polarity from the subduction zone feature to its current sub-segment feature.
                    # They are different features and hence only the subduction zone itself will have a subduction polarity.
                    resolved_subduction_feature.set_enumeration(
                        pygplates.PropertyName.gpml_subduction_polarity,
                        subduction_polarity,
                        # Just in case the sub-sub-segment feature type does not support a subduction polarity...
                        verify_information_model=pygplates.VerifyInformationModel.no)
                    
                    resolved_subduction_features.append(resolved_subduction_feature)
                    
            else: # It's not a topological line...
                
                # Create the resolved subduction zone feature.
                resolved_subduction_feature = shared_sub_segment.get_resolved_feature()
                
                resolved_subduction_features.append(resolved_subduction_feature)
                
    return resolved_subduction_features


# Tessellates all resolved subduction zones to within 'threshold_sampling_distance_radians' radians and
# returns the following parameters at each tessellates point:
#
# - point longitude
# - point latitude
# - subduction zone absolute (relative to anchor plate) velocity magnitude (in cm/yr)
# - subduction zone absolute velocity obliquity angle (angle between subduction zone normal vector and absolute velocity vector)
# - length of arc segment (in degrees) that current point is on
# - subducting arc normal azimuth angle (clockwise starting at North, ie, 0 to 360 degrees) at current point
#
# Also first applies the rotation adjustment 'absolute_rotation_adjustment' to the resolved geometries
# (this assumes the topologies were resolved using a rotation model that did not have this adjustment).
# Although note that the argument 'rotation_features_or_model' must include the rotation adjustment
# since it is used to calculate the velocities.
#
# The reason for this is an optimisation whereby the topologies only need be resolved once but the
# absolute motions determined many times with different absolute rotation adjustments.
#
# NOTE: The topologies should have already been resolved for time 'time' and their resolved geometries stored
# in features (which are passed in the 'resolved_topology_features' argument). The subduction zone plate ID
# and polarities should also be stored in the features.
#
# Note that the absolute velocity magnitudes are negative if the subduction zones absolute motion is heading
# in the direction of the overriding plates (if absolute obliquity angle is less than 90 and greater than -90).

def subduction_absolute_motion(
        # Rotation model or feature collection(s), or list of features, or filename(s)...
        rotation_features_or_model,
        # Resolved topology feature collection(s), or list of features, or filename(s) or any combination of those...
        resolved_topology_features,
        # Threshold sampling distance along subduction zones (in radians)...
        threshold_sampling_distance_radians,
        time,
        # Absolute rotation to apply to 'resolved_topology_features' (default to identity/zero rotation).
        absolute_rotation_adjustment = pygplates.FiniteRotation(),
        velocity_delta_time = 1.0,
        anchor_plate_id = 0):
    
    # Turn rotation data into a RotationModel (if not already).
    rotation_model = pygplates.RotationModel(rotation_features_or_model)
    
    # Turn resolved topology data into a list of features (if not already).
    resolved_topology_features = pygplates.FeaturesFunctionArgument(resolved_topology_features)
    
    # List of tesselated subduction zone shared subsegment points and associated convergence parameters
    # for the current 'time'.
    output_data = []
    
    # Iterate over the shared boundary sections of all resolved topologies.
    for shared_sub_segment_feature in resolved_topology_features.get_features():
        
        subduction_zone_plate_id = shared_sub_segment_feature.get_reconstruction_plate_id()
        
        # Get the rotation of the subduction zone relative to the anchor plate
        # from 'time + velocity_delta_time' to 'time'.
        #
        # Note: We don't need to convert to and from the stage rotation reference frame
        # like the above convergence because this stage rotation is relative to the anchor plate
        # and so the above to/from stage rotation frame conversion "R(0->t2,A->F)" is the
        # identity rotation since the fixed plate (F) is the anchor plate (A).
        subduction_zone_equivalent_stage_rotation = rotation_model.get_rotation(
                time,
                subduction_zone_plate_id,
                time + velocity_delta_time,
                anchor_plate_id=anchor_plate_id)
        
        # We need to reverse the subducting_normal vector direction if overriding plate is to
        # the right of the subducting line since great circle arc normal is always to the left.
        subduction_polarity = shared_sub_segment_feature.get_enumeration(pygplates.PropertyName.gpml_subduction_polarity)
        if (not subduction_polarity) or (subduction_polarity == 'Unknown'):
            warnings.warn('Unknown subduction polarity in feature "{0}".\n'.format(shared_sub_segment_feature.get_name()), RuntimeWarning)
            continue
        elif subduction_polarity == 'Left':
            subducting_normal_reversal = 1
        else:
            subducting_normal_reversal = -1
        
        shared_sub_segment_geometry = shared_sub_segment_feature.get_geometry()
        if not shared_sub_segment_geometry:
            warnings.warn('Missing subduction geometry in feature "{0}".\n'.format(shared_sub_segment_feature.get_name()), RuntimeWarning)
            continue
        
        # Apply absolute rotation adjustment to pre-resolved topology.
        shared_sub_segment_geometry = absolute_rotation_adjustment * shared_sub_segment_geometry
        
        # Ensure the shared sub-segment is tessellated to within the threshold sampling distance.
        tessellated_shared_sub_segment_polyline = shared_sub_segment_geometry.to_tessellated(threshold_sampling_distance_radians)
        
        # Iterate over the great circle arcs of the tessellated polyline to get the
        # arc midpoints, lengths and subducting normals.
        # There is an arc between each adjacent pair of points in the polyline.
        arc_midpoints = []
        arc_lengths = []
        subducting_arc_normals = []
        for arc in tessellated_shared_sub_segment_polyline.get_segments():
            if not arc.is_zero_length():
                arc_midpoints.append(arc.get_arc_point(0.5))
                arc_lengths.append(arc.get_arc_length())
                # The normal to the subduction zone in the direction of subduction (towards overriding plate).
                subducting_arc_normals.append(subducting_normal_reversal * arc.get_great_circle_normal())
        
        # Shouldn't happen, but just in case the shared sub-segment polyline coincides with a point.
        if not arc_midpoints:
            warnings.warn('Sub-segment polyline coincides with a point in feature "{0}".\n'.format(shared_sub_segment_feature.get_name()), RuntimeWarning)
            continue
        
        # The subducting arc normals relative to North (azimuth).
        # Convert global 3D normal vectors to local (magnitude, azimuth, inclination) tuples (one tuple per point).
        subducting_arc_local_normals = pygplates.LocalCartesian.convert_from_geocentric_to_magnitude_azimuth_inclination(
                arc_midpoints, subducting_arc_normals)
        
        # Calculate the absolute velocities at the arc midpoints.
        absolute_velocity_vectors = pygplates.calculate_velocities(
                arc_midpoints, subduction_zone_equivalent_stage_rotation,
                velocity_delta_time, pygplates.VelocityUnits.cms_per_yr)
        
        for arc_index in range(len(arc_midpoints)):
            arc_midpoint = arc_midpoints[arc_index]
            arc_length = arc_lengths[arc_index]
            subducting_arc_normal = subducting_arc_normals[arc_index]
            subducting_arc_normal_azimuth = subducting_arc_local_normals[arc_index][1]
            lat, lon = arc_midpoint.to_lat_lon()
            
            # The direction towards which we rotate from the subducting normal in a clockwise fashion.
            clockwise_direction = pygplates.Vector3D.cross(subducting_arc_normal, arc_midpoint.to_xyz())
            
            # Calculate the absolute rate parameters.
            absolute_velocity_vector = absolute_velocity_vectors[arc_index]
            if absolute_velocity_vector.is_zero_magnitude():
                absolute_velocity_magnitude = 0
                absolute_obliquity_degrees = 0
            else:
                absolute_velocity_magnitude = absolute_velocity_vector.get_magnitude()
                absolute_obliquity_degrees = math.degrees(pygplates.Vector3D.angle_between(
                        absolute_velocity_vector, subducting_arc_normal))
                # Anti-clockwise direction has range (0, -180) instead of (0, 180).
                if pygplates.Vector3D.dot(absolute_velocity_vector, clockwise_direction) < 0:
                    absolute_obliquity_degrees = -absolute_obliquity_degrees
                
                # See if the subduction zone absolute motion is heading in the direction of the
                # overriding plate. If it is then make the velocity magnitude negative to
                # indicate this. This could be inferred from the obliquity but it seems this
                # is the standard way to output convergence rate.
                #
                # Note that we are not calculating the motion of the subduction zone
                # relative to the overriding plate - they are usually attached to each other
                # and hence wouldn't move relative to each other.
                if math.fabs(absolute_obliquity_degrees) < 90:
                    absolute_velocity_magnitude = -absolute_velocity_magnitude
            
            # The data will be output in GMT format (ie, lon first, then lat, etc).
            output_data.append((
                    lon,
                    lat,
                    absolute_velocity_magnitude,
                    absolute_obliquity_degrees,
                    math.degrees(arc_length),
                    math.degrees(subducting_arc_normal_azimuth)))
    
    # Return data sorted since it's easier to compare results (when at least lon/lat is sorted).
    return sorted(output_data)


if __name__ == '__main__':
    
    def adjust_rotation_features(rotation_features, absolute_rotation_adjustment):
        for rotation_feature in rotation_features:
            
            total_reconstruction_pole = rotation_feature.get_total_reconstruction_pole()
            if not total_reconstruction_pole:
                # Not a rotation feature.
                continue

            fixed_plate_id, moving_plate_id, rotation_sequence = total_reconstruction_pole
            # If found the absolute plate motion rotation feature then reset it to the adjustment rotation.
            if moving_plate_id == 1 and fixed_plate_id == 0:
                # Set two finite rotations to the adjustment rotation at 0Ma and 600Ma.
                rotation_sequence[:] = [
                    pygplates.GpmlTimeSample(pygplates.GpmlFiniteRotation(absolute_rotation_adjustment), time)
                        for time in (0, 600)]
                
                return
    
    
    time = 1
    threshold_sampling_distance_radians = math.radians(1.0)
    absolute_rotation_adjustment = pygplates.FiniteRotation(pygplates.PointOnSphere.north_pole, math.radians(5))
    
    rotation_features = pygplates.FeatureCollection('../../../../sample_data/2.0/SampleData/FeatureCollections/Rotations/Matthews_etal_GPC_2016_410-0Ma_GK07.rot')
    topology_features = [
        pygplates.FeatureCollection('../../../../sample_data/2.0/SampleData/FeatureCollections/DynamicPolygons/Matthews_etal_GPC_2016_MesozoicCenozoic_PlateTopologies.gpmlz'),
        pygplates.FeatureCollection('../../../../sample_data/2.0/SampleData/FeatureCollections/DynamicPolygons/Matthews_etal_GPC_2016_Paleozoic_PlateTopologies.gpmlz'),
        pygplates.FeatureCollection('../../../../sample_data/2.0/SampleData/FeatureCollections/DynamicPolygons/Matthews_etal_GPC_2016_TopologyBuildingBlocks.gpmlz')]
    
    # One rotation model without adjustment and one with.
    original_rotation_model = pygplates.RotationModel(rotation_features)
    adjust_rotation_features(rotation_features, absolute_rotation_adjustment)
    adjusted_rotation_model = pygplates.RotationModel(rotation_features)
    
    # Resolved topologies with and without adjustment.
    resolved_subduction_zone_features = resolve_subduction_zones(original_rotation_model, topology_features, time)
    adjusted_resolved_subduction_zone_features = resolve_subduction_zones(adjusted_rotation_model, topology_features, time)
    if not resolved_subduction_zone_features or not adjusted_resolved_subduction_zone_features:
        sys.exit(1)
    
    # Can write intermediate results to a GPML file if desired.
    #
    #resolved_topologies_filename = 'tmp.gpml'
    #pygplates.FeatureCollection(resolved_subduction_zone_features).write(resolved_topologies_filename)
    
    # Apply adjustment to topologies that were resolved *without* adjustment.
    data1 = subduction_absolute_motion(
            adjusted_rotation_model,
            resolved_subduction_zone_features,
            threshold_sampling_distance_radians,
            time,
            absolute_rotation_adjustment)
    # Don't apply adjustment to topologies that were resolved *with* adjustment.
    data2 = subduction_absolute_motion(
            adjusted_rotation_model,
            adjusted_resolved_subduction_zone_features,
            threshold_sampling_distance_radians,
            time)
    
    # Both datasets should be the same within numerical tolerance.
    # Can visually compare the following text files.
    
    with open('data1.txt', 'w') as output_file:
        for line in data1:
            output_file.write(' '.join(str(item) for item in line) + '\n')
    
    with open('data2.txt', 'w') as output_file:
        for line in data2:
            output_file.write(' '.join(str(item) for item in line) + '\n')
