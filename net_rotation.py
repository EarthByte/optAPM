
"""
    Copyright (C) 2018 The University of Sydney, Australia
    
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

############################################################################################
#
# Python implementation of the Net Rotation export currently in GPlates (versions <= 2.1).
#
# An effort has been made to reproduce the implementation in GPlates as closely as possible.
# This does not preclude writing a better net rotation implementation (such as using uniformly
# spaced points on the globe, instead of the uniform lat-lon spacing used here).
#
# NOTE: Like GPlates, this implementation only works with *non-deforming* topologies.
#
############################################################################################


import math
import points_in_polygons
import pygplates

#
# Python 2 and 3 compatibility.
#
# Iterating over a dict.
try:
    dict.iteritems
except AttributeError:
    # Python 3
    def itervalues(d):
        return iter(d.values())
    def iteritems(d):
        return iter(d.items())
    def listvalues(d):
        return list(d.values())
    def listitems(d):
        return list(d.items())
else:
    # Python 2
    def itervalues(d):
        return d.itervalues()
    def iteritems(d):
        return d.iteritems()
    def listvalues(d):
        return d.values()
    def listitems(d):
        return d.items()


class VelocityMethod(object):
    """
    An enumeration (class with static integers) of velocity methods.
    
    This is essentially the GPlates enumeration GPlatesQtWidgets::VelocityMethodWidget::VelocityMethod.
    """
    T_TO_T_MINUS_DT = 0       # Velocity calculated over time interval [time, time - delta].
    T_PLUS_DT_TO_T = 1        # Velocity calculated over time interval [time + delta, time].
    T_PLUS_MINUS_HALF_DT = 2  # Velocity calculated over time interval [time + delta/2, time - delta/2].

# Default velocity method used in GPlates...
DEFAULT_VELOCITY_METHOD = VelocityMethod.T_TO_T_MINUS_DT

# Default velocity delta time used in GPlates...
DEFAULT_VELOCITY_DELTA_TIME = 10.0


def generate_net_rotation_grid_points_internal_gplates():
    """
    Creates a multi-point that matches those points used in the Net Rotation export in GPlates (versions <= 2.1).
    
    These are uniformly lat-lon spaced at 1 degree.
    """
    
    net_rotation_points = []
    
    for lat in range(-90, 91, 1):
        for lon in range(-180, 181, 1):
            net_rotation_points.append((lat, lon))
    
    return pygplates.MultiPointOnSphere(net_rotation_points)

# Default uniform lat-lon points at 1 degree resolution used in the Net Rotation export in GPlates.
DEFAULT_POINTS = generate_net_rotation_grid_points_internal_gplates()


def calculate_net_rotation_internal_gplates(
        rotation_features_or_model,
        topology_features,
        reconstruction_time,
        points = DEFAULT_POINTS,
        velocity_method = DEFAULT_VELOCITY_METHOD,
        velocity_delta_time = DEFAULT_VELOCITY_DELTA_TIME):
    """
    Calculate the net rotation over all *non-deforming* topological plates matching that of the
    Net Rotation export in GPlates (versions <= 2.1).
    
    Returns a pygplates.FiniteRotation.
    
    NOTE: Currently the code assumes a uniform lat-lon distribution of 'points' (in the cos(latitude) area weighting).
    
    NOTE: This only works with *non-deforming* topologies.
    """
    
    # Get the net rotations results for each individual plate as a key/value dictionary of plate_id/NetRotationResult.
    plate_net_rotations = calculate_plate_net_rotation_results_internal_gplates(
            rotation_features_or_model,
            topology_features,
            reconstruction_time,
            points,
            velocity_method,
            velocity_delta_time)
    
    # Iterate over the plates and accumulate unweighted rotation component of each plate and total weight across all plates.
    total_net_rotation_result = NetRotationResult.create_zero()
    for plate_id, plate_net_rotation_result in iteritems(plate_net_rotations):
        #print 'Plate ID:', plate_id
        #plate_net_rotation_result.debug_print()
        #print 'Net rotation:', plate_net_rotation_result.get_net_rotation()
        total_net_rotation_result.accumulate(plate_net_rotation_result)
    
    return total_net_rotation_result.get_net_rotation()


class NetRotationResult(object):
    """
    This is essentially the GPlates class GPlatesAppLogic::NetRotationUtils::NetRotationResult.
    
    According to GPlates:
    
      This is used for storing intermediate results during point-by-point net-rotation calculations.
      Each point used in the net-rotation calculation has its results stored here, which are later summed.
    """
    
    @staticmethod
    def create_zero():
        return NetRotationResult(pygplates.Vector3D.zero, 0.0, 0.0, 0.0)
    
    def __init__(
            self,
            rotation_component,
            weighting_factor,
            plate_area_component,
            plate_angular_velocity):
        
        self.rotation_component = rotation_component
        self.weighting_factor = weighting_factor
        self.plate_area_component = plate_area_component
        self.plate_angular_velocity = plate_angular_velocity
    
    def get_net_rotation(self):
        """
        Calculate weighted net rotation.
        
        Returns a pygplates.FiniteRotation.
        """
        # If zero weighting factor then return identity rotation.
        if self.weighting_factor == 0:
            return pygplates.FiniteRotation()
        
        weighted_rotation_component = (1.0 / self.weighting_factor) * self.rotation_component
        
        # If zero then return identity rotation.
        if weighted_rotation_component.is_zero_magnitude():
            return pygplates.FiniteRotation()
        
        # Extract total pole/angle from rotation rate vector.
        weighted_rotation_pole, weighted_rotation_angle_degrees_per_my = convert_net_rotation_xyz_to_pole(weighted_rotation_component)
        
        # Return total pole/angle as a FiniteRotation.
        return pygplates.FiniteRotation(
                weighted_rotation_pole,
                math.radians(weighted_rotation_angle_degrees_per_my))
    
    def accumulate(self, other_net_rotation_result):
        self.rotation_component = self.rotation_component + other_net_rotation_result.rotation_component
        self.weighting_factor += other_net_rotation_result.weighting_factor
        self.plate_area_component += other_net_rotation_result.plate_area_component
        # Note that accumulating 'plate_angular_velocity' is also absent from the GPlates function
        # GPlatesAppLogic::NetRotationUtils::sum_net_rotations().
        # It appears that 'plate_angular_velocity' is essentially unused in GPlates.
    
    def debug_print(self):
        print 'weighting_factor:', self.weighting_factor
        print 'Omega pre weight:', self.rotation_component
        if self.weighting_factor != 0:
            omega_post_weight = (1.0 / self.weighting_factor) * self.rotation_component
            print 'Omega post weight:', omega_post_weight
        net_rotation = self.get_net_rotation()
        if not net_rotation.represents_identity_rotation():
            omega_pole_latitude, omega_pole_longitude, omega_angle_degrees_per_my = net_rotation.get_lat_lon_euler_pole_and_angle_degrees()
            print 'Omega (pole, angle):', omega_pole_latitude, omega_pole_longitude, omega_angle_degrees_per_my


def calculate_plate_net_rotation_results_internal_gplates(
        rotation_features_or_model,
        topology_features,
        reconstruction_time,
        points = DEFAULT_POINTS,
        velocity_method = DEFAULT_VELOCITY_METHOD,
        velocity_delta_time = DEFAULT_VELOCITY_DELTA_TIME):
    """
    Calculate net rotation results for *non-deforming* topological plates matching that of the
    Net Rotation export in GPlates (versions <= 2.1).
    
    Returns a dictionary indexed by integer plate ID (keys) with values of NetRotationResult that
    are summed over those of 'points' that lie within each plate.
    
    NOTE: Currently the code assumes a uniform lat-lon distribution of 'points' (in the cos(latitude) area weighting).
    
    NOTE: This only works with *non-deforming* topologies.
    """
    
    if velocity_delta_time <= 0:
        raise ValueError("'velocity_delta_time' should be positive")
    
    # Turn rotation data into a RotationModel (if not already).
    rotation_model = pygplates.RotationModel(rotation_features_or_model)
    
    # Turn topology data into a list of features (if not already).
    topology_features = pygplates.FeaturesFunctionArgument(topology_features).get_features()
    
    # Resolve the plate polygons for the current time.
    resolved_topologies = []
    pygplates.resolve_topologies(topology_features, rotation_model, resolved_topologies, reconstruction_time,
        # NOTE: We are only looking at *non-deforming* boundaries...
        resolve_topology_types=pygplates.ResolveTopologyType.boundary)
    
    from_time, to_time = get_from_and_to_times(reconstruction_time, velocity_method, velocity_delta_time)
    
    # Dictionary of stage poles (in xyz form) indexed by plate ID.
    stage_pole_xyz_dict = {}
    
    # Each resolved boundary is a polygon.
    resolved_topology_polygons = []
    
    # Fill stage pole dictionary and extract resolved boundary polygons.
    for resolved_topology in resolved_topologies:
        resolved_topology_polygons.append(resolved_topology.get_resolved_boundary())
        
        plate_id = resolved_topology.get_feature().get_reconstruction_plate_id()
        if plate_id not in stage_pole_xyz_dict:
            stage_pole = get_stage_pole_internal_gplates(rotation_model, to_time, plate_id, from_time, 0)
            
            if stage_pole.represents_identity_rotation():
                # Store as zero vector.
                stage_pole_xyz_dict[plate_id] = pygplates.Vector3D.zero
            else:
                pole, angle_radians = stage_pole.get_euler_pole_and_angle()
                
                # Convert angle to degrees per My.
                angle_degrees_per_my = math.degrees(angle_radians) / velocity_delta_time
                
                # Store stage pole in xyz form.
                stage_pole_xyz_dict[plate_id] = convert_net_rotation_pole_to_xyz(pole, angle_degrees_per_my)
    
    # Dictionary of net rotations indexed by plate ID.
    plate_net_rotations = {}
    
    #
    # Assuming the resolved polygons are *non-overlapping*, which they should be (except for minor overap/gap errors),
    # find the single polygon (resolved boundary) containing each point.
    #
    point_plate_ids = points_in_polygons.find_polygons(
            points,
            resolved_topology_polygons,
            # The plate ID of each resolved topology (this is what is returned by 'find_polygons'...
            [resolved_topology.get_feature().get_reconstruction_plate_id() for resolved_topology in resolved_topologies])
    
    # Each point is contained by one resolved plate boundary.
    # Iterate over them and sum the point's contribution to that plate's net rotation.
    for point_index, point_plate_id in enumerate(point_plate_ids):
        if point_plate_id is None:
            # Point is not in any resolved boundary.
            # It either fell in a tiny crack/gap or the topologies don't have global coverage.
            continue
        
        point = points[point_index]
        
        # Get the contribution of current point in the plate containing it.
        net_rotation_contribution = calc_net_rotation_contribution(
                point,
                stage_pole_xyz_dict[point_plate_id])
        
        # Each plate starts at zero contribution.
        if point_plate_id not in plate_net_rotations:
            plate_net_rotations[point_plate_id] = NetRotationResult.create_zero()
        
        # Accumulate the contribution of the current point into the plate containing it.
        plate_net_rotations[point_plate_id].accumulate(net_rotation_contribution)
    
    return plate_net_rotations


def get_from_and_to_times(
        reconstruction_time,
        velocity_method,
        velocity_delta_time):
    """This is essentially the GPlates function get_older_and_younger_times() in 'ExportNetRotationAnimationStrategy.cc'."""
    
    if velocity_method == VelocityMethod.T_TO_T_MINUS_DT:
        return reconstruction_time, reconstruction_time - velocity_delta_time
    elif velocity_method == VelocityMethod.T_PLUS_DT_TO_T:
        return reconstruction_time + velocity_delta_time, reconstruction_time
    elif velocity_method == VelocityMethod.T_PLUS_MINUS_HALF_DT:
        return reconstruction_time + 0.5 * velocity_delta_time, reconstruction_time - 0.5 * velocity_delta_time
    
    raise ValueError("Unexpected value for 'velocity_method'")


def calc_net_rotation_contribution(
        point,
        stage_pole_xyz):
    """
    This is essentially the GPlates function GPlatesAppLogic::NetRotationUtils::calc_net_rotation_contribution().
    
    Note: As an optimisation (that differs from the GPlates source code), part of this function has been moved into the
    calling code (specifically the part that converts stage pole into xyz form, since only need to do once per stage pole).
    """
    
    point_x, point_y, point_z = point.to_xyz()
    point_xyz = pygplates.Vector3D(point_x, point_y, point_z)
    
    # cos(latitude)^2 is same as (x*x + y*y).
    cos_latitude_squared = point_x * point_x + point_y * point_y
    
    # This is stored in the returned NetRotationResult as the plate area component.
    cos_latitude = math.sqrt(cos_latitude_squared)
    
    # This is the vector triple cross-product "R x (W x R)" mentioned in the Torsvik 2010 paper:
    #   "Plate tectonics and net lithosphere rotation over the past 150 My".
    omega = pygplates.Vector3D.cross(
            point_xyz,
            pygplates.Vector3D.cross(stage_pole_xyz, point_xyz))
    
    
    # Omega is multiplied by 'cos_latitude' to account for the unequal area lat/lon spacing near the poles.
    # In other words, the surface integral over the unit sphere in spherical coordinates is:
    #
    #       / /                             /2pi       /pi
    #       | |                             |          |
    #   I = | | omega(cotheta, phi) |dA| =  | d(phi)   | omega(cotheta, phi) * sin(cotheta) d(cotheta)
    #       | |                             |          |
    #      / /                             / 0        / 0
    #
    #                                       /2pi       /pi/2
    #                                       |          |
    #                                    =  | d(phi)   | omega(theta, phi) * cos(theta) d(theta)
    #                                       |          |
    #                                      / 0        /-pi/2
    #
    # ...where sin(colatitude) = cos(latitude).
    # So that accounts for the cos(latitude) that we are using.
    omega = omega * cos_latitude
    
    # Later on we sum up the weighting factors over a plate (or over entire globe) and also sum up the cosine-weighted
    # omega contributions (noted above) and then divide the two to get the global or plate net rotation.
    # This is "sum(omega * cos_latitude) / sum(weighting_factor)", with the weighting factor in the denominator being
    # "cos_latitude * (x^2 + y^2)" which is the same as "cos(latitude) * cos(latitude)^2" which is the same as "cos(latitude)^3".
    # The reason for this weighting (normalization) factor is as follows...
    #
    # In this example, noted by Bernhard Steinberger, we have a single plate with a rotation 'Wz' about the global z-axis.
    # We want the net rotation integral (noted above in spherical coordinates) to equal 'Wz', so we have a normalisation factor 'k'
    # such that the right hand side equals the left hand side in the following
    # (where W = <0,0,Wz> is 3D rotation vector, and R = <x,y,z> is 3D position on unit-sphere) ...
    #
    #          /2pi       /pi/2
    #          |          |
    #   W =  k | d(phi)   | R x (W x R) cos(theta) d(theta)
    #          |          |
    #         / 0        /-pi/2
    #
    #          /2pi       /pi/2
    #          |          |
    #     =  k | d(phi)   | [W(R.R) - R(R.W)] cos(theta) d(theta)
    #          |          |
    #         / 0        /-pi/2
    #
    # ...where "R x (W x R) = [W(R.R) - R(R.W)]" using vector triple product expansion - https://en.wikipedia.org/wiki/Triple_product#Vector_triple_product
    #
    #                  /2pi       /pi/2
    #                  |          |
    #   <0,0,Wz>  =  k | d(phi)   | [<0,0,Wz> - <x,y,z> Wz sin(theta)] cos(theta) d(theta)
    #                  |          |
    #                 / 0        /-pi/2
    #
    # ...where we used "R.R = 1" since dot product of a unit vector with itself is the scalar one, and "R.W = <x,y,z>.<0,0,Wz> = Wz z = Wz sin(theta)"
    #
    #                  /2pi       /pi/2
    #                  |          |
    #   <0,0,Wz>  =  k | d(phi)   | [<0,0,Wz> - <cos(phi) cos(theta), sin(phi) cos(theta), sin(theta)> Wz sin(theta)] cos(theta) d(theta)
    #                  |          |
    #                 / 0        /-pi/2
    #
    # ...and noting that...
    #
    #    /2pi                 /2pi     
    #    |                    |        
    #    | d(phi) cos(phi) =  | d(phi) sin(phi) =  0
    #    |                    |        
    #   / 0                  / 0       
    #
    # ...we get...
    #
    #                  /2pi       /pi/2
    #                  |          |
    #   <0,0,Wz>  =  k | d(phi)   | [<0,0,Wz> - <0,0,sin(theta)> Wz sin(theta)] cos(theta) d(theta)
    #                  |          |
    #                 / 0        /-pi/2
    #
    # ...where only the z-component of vector on both sides is non-zero, which we write as...
    #
    #                  /2pi       /pi/2
    #                  |          |
    #         Wz  =  k | d(phi)   | [Wz - Wz sin(theta)^2] cos(theta) d(theta)
    #                  |          |
    #                 / 0        /-pi/2
    #
    #                             /pi/2
    #                             |
    #             =  k (2 pi) Wz  | [1 - sin(theta)^2] cos(theta) d(theta)
    #                             |
    #                            /-pi/2
    #
    #                             /pi/2
    #                             |
    #             =  k (2 pi) Wz  | cos(theta)^3 d(theta)
    #                             |
    #                            /-pi/2
    #
    # ...which results in the normalization factor 'k' being...
    #
    #                               1
    #         k =   ------------------------------------
    #                            /pi/2
    #                            |
    #                       2 pi | cos(theta)^3 d(theta)
    #                            |
    #                           /-pi/2
    #
    # ...hence the weighting factor (the denominator of 'k') is the integral/summation of "cos(latitude)^3".
    #
    # Also note that this integral for 'k' evaluates to (using "integral[cos(theta)^3] = sin(theta) - (1/3) sin(theta)^3"):
    #
    #         k =   3 / (8 pi)
    #
    # ...which is normalization factor mentioned in the Torsvik 2010 paper:
    #   "Plate tectonics and net lithosphere rotation over the past 150 My".
    #
    weighting_factor = cos_latitude_squared * cos_latitude  # cos(latitude)^3
    
    # We do, however, store 'cos_latitude' as the plate area (but plate area does not seemed to be used anywhere).
    return NetRotationResult(omega, weighting_factor, cos_latitude, 0.0)


def convert_net_rotation_xyz_to_pole(xyz):
    """This is essentially the GPlates function GPlatesAppLogic::NetRotationUtils::convert_net_rotation_xyz_to_pole()"""
    
    pole = pygplates.PointOnSphere(xyz.to_normalised().to_xyz())
    angle_degrees_per_my = xyz.get_magnitude()
    
    return (pole, angle_degrees_per_my)


def convert_net_rotation_pole_to_xyz(pole, angle_degrees_per_my):
    """This is essentially the GPlates function GPlatesAppLogic::NetRotationUtils::convert_net_rotation_pole_to_xyz()"""
    
    # This looks to be the 'rotation rate vector' mentioned in the Torsvik 2010 paper:
    #   "Plate tectonics and net lithosphere rotation over the past 150 My".
    return angle_degrees_per_my * pygplates.Vector3D(pole.to_xyz());


def get_stage_pole_internal_gplates(
        rotation_model,
        to_time,
        moving_plate_id,
        from_time,
        fixed_plate_id):
    """
    Calculate a stage pole the way it's done internally in GPlates.
    
    This is similar to just calling:
    
        rotation_model.get_rotation(to_time, moving_plate_id, from_time, fixed_plate_id)
    
    ...but GPlates will go ahead and use an identity rotation for each equivalent total rotation
    (in the stage pole calculation) when a plate ID is missing (for the particular time), thus the
    final stage pole can be the result of a mixture of identity and non-identity rotations.
    Whereas pygplates.RotationModel.get_rotation() will return an identity rotation if *any* plate IDs
    are missing at any to/from times.
    
    An attempt was made to change this in GPlates but it caused issues with certain datasets
    (such as topologies) that rely on the existing behaviour. For example, topologies broke when this
    was changed in GPlates, so the change was reverted.
    """
    
    #
    # The following code is essentially copied from the GPlates function
    # GPlatesAppLogic::RotationUtils::get_stage_pole().
    #
    # ...where 't1' is 'from_time', 't2' is 'to_time', 'M' is 'moving_plate_id' and 'F' is 'fixed_plate_id'.
    #
    
    # For t1, get the rotation for plate M w.r.t. anchor
    finite_rot_0_to_t1_M = rotation_model.get_rotation(from_time, moving_plate_id)

    # For t1, get the rotation for plate F w.r.t. anchor
    finite_rot_0_to_t1_F = rotation_model.get_rotation(from_time, fixed_plate_id)


    # For t2, get the rotation for plate M w.r.t. anchor
    finite_rot_0_to_t2_M = rotation_model.get_rotation(to_time, moving_plate_id)

    # For t2, get the rotation for plate F w.r.t. anchor
    finite_rot_0_to_t2_F = rotation_model.get_rotation(to_time, fixed_plate_id)

    # Compose these rotations so that we get
    # the stage pole from time t1 to time t2 for plate M w.r.t. plate F.

    finite_rot_t1 = finite_rot_0_to_t1_F.get_inverse() * finite_rot_0_to_t1_M

    finite_rot_t2 = finite_rot_0_to_t2_F.get_inverse() * finite_rot_0_to_t2_M

    stage_pole = finite_rot_t2 * finite_rot_t1.get_inverse()
    
    return stage_pole


if __name__ == '__main__':
    
    #
    # Some test code that finds the total net rotation at 50Ma of the dynamic plates in the 2.1 sample data.
    #
    
    import os.path
    import time
    
    sample_data_path = 'C:/Users/John/Development/Usyd/gplates/sample_data/2.1/SampleData'
    
    rotation_features = pygplates.FeatureCollection(os.path.join(sample_data_path, 'FeatureCollections/Rotations/Matthews_etal_GPC_2016_410-0Ma_GK07.rot'))
    
    topology_features = [
        pygplates.FeatureCollection(os.path.join(sample_data_path, 'FeatureCollections/DynamicPolygons/Matthews_etal_GPC_2016_MesozoicCenozoic_PlateTopologies.gpmlz')),
        pygplates.FeatureCollection(os.path.join(sample_data_path, 'FeatureCollections/DynamicPolygons/Matthews_etal_GPC_2016_Paleozoic_PlateTopologies.gpmlz')),
        pygplates.FeatureCollection(os.path.join(sample_data_path, 'FeatureCollections/DynamicPolygons/Matthews_etal_GPC_2016_TopologyBuildingBlocks.gpmlz'))]
    
    start_time = time.clock()
    
    net_rotation = calculate_net_rotation_internal_gplates(
            rotation_features,
            topology_features,
            reconstruction_time = 50.0,
            #velocity_method = VelocityMethod.T_TO_T_MINUS_DT,
            #velocity_method = VelocityMethod.T_PLUS_DT_TO_T,
            #velocity_method = VelocityMethod.T_PLUS_MINUS_HALF_DT,
            #velocity_delta_time = 1.0
            )
    
    end_time = time.clock()
    
    print 'Net rotation:', net_rotation
    
    print 'Net rotation calculation took', (end_time - start_time), 'seconds'
