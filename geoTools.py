"""
    Copyright (C) 2016  Michael G. Tetley
    EarthByte Group, University of Sydney
    Geological and Planetary Sciences, California Institute of Technology
    Contact email: michael.tetley@sydney.edu.au // mtetley@caltech.edu

    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation; either version 2
    of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.


    ## Geoscience Toolbox ##

    Set of tools to calculate useful things for geoscience and palaeomagnetics research
"""

# Import required libraries
import numpy as np
import pygplates as pgp
import pmagpy.ipmag as ipmag

import pylab

# The optimization workflow doesn't actually need to plot so we won't require user to install these modules.
# If the user plots then we'll get an AttributeError, in which case the try/except part should be removed.
try:
    import matplotlib.pyplot as plt
    from mpl_toolkits.basemap import Basemap
except ImportError:
    pass


""" GEOSCIENCE """
"""
    haversine

    Module to calculate the great circle distance and bearing between two points on a sphere.
    Returns distance in kilometers and bearing in degrees.
"""
def haversine(lon1, lat1, lon2, lat2):

    # convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])

    dlon = lon2 - lon1
    dlat = lat2 - lat1

    a = np.sin(dlat / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2) ** 2
    c = 2 * np.arcsin(np.sqrt(a))

    # Radius of the Earth
    km = 6371 * c

    bearing = np.arctan2(np.sin(lon2-lon1) * np.cos(lat2), np.cos(lat1) * np.sin(lat2) - np.sin(lat1) * np.cos(lat2) * np.cos(lon2-lon1))
    bearing = np.degrees(bearing)
    bearing = (bearing + 360) % 360

    return km, np.rad2deg(c), bearing


"""
    global_points_rand

    Module to generate a random distribution of lat / lon points on the Earth
    Returns length = samples list of latitudes and longitudes
"""
def global_points_rand(samples):

    lats = []
    lons = []

    for i in xrange(0, samples):

        theta = 2 * np.pi * np.random.random()
        phi = np.arccos(2 * np.random.random() - 1.0)

        x = np.cos(theta) * np.sin(phi)
        y = np.sin(theta) * np.sin(phi)
        z = np.cos(phi)

        point = pgp.convert_point_on_sphere_to_lat_lon_point((x,y,z))
        lats.append(point.get_latitude())
        lons.append(point.get_longitude())

    return lats, lons


"""
    global_points_uniform

    Module to generate a uniform (even) distribution of lat / lon points on the Earth
    Returns length = samples list of latitudes and longitudes
"""
def global_points_uniform(samples, plotResult=False, projection='robin'):

    lats = []
    lons = []

    angle = np.pi * (3 - np.sqrt(5))
    theta = angle * np.arange(samples)
    z = np.linspace(1 - 1.0 / samples, 1.0 / samples - 1, samples)
    radius = np.sqrt(1 - z * z)

    points = np.zeros((samples, 3))
    points[:,0] = radius * np.cos(theta)
    points[:,1] = radius * np.sin(theta)
    points[:,2] = z

    for i in xrange(0, len(points)):

        point = pgp.convert_point_on_sphere_to_lat_lon_point((points[i][0], points[i][1], points[i][2]))
        lats.append(point.get_latitude())
        lons.append(point.get_longitude())

    if plotResult == True:

        # Plot start seeds
        m = Basemap(projection=projection,lat_0=0,lon_0=0,resolution='c',area_thresh=50000)
        plt.figure(figsize=(7, 7))
        plt.title("n = " + str(samples) + " uniformly distributed global 'seed' locations")
        m.drawcoastlines(linewidth=0.25)
        m.fillcontinents(color='bisque',zorder=1)
        m.drawmeridians(np.arange(0,360,30))
        m.drawparallels(np.arange(-90,90,30))

        ipmag.plot_vgp(m, lons, lats)

        plt.show()


    return lats, lons


"""
    checkLatLon

    Simple function to check incoming lat and lon are within correct ranges.
    
    If the longitude is outside [-180, 180] then it is wrapped by 360 degrees to bring it into that range.
    If the latitude is outside [-90, 90] then it is moved to the other side of the North (90) or South (-90) pole.
"""
def checkLatLon(lat, lon):

    if lon > 180:

        lon_corrected = -360 + lon

    elif lon < -180:

        lon_corrected = lon + 360

    else:

        lon_corrected = lon


    if lat > 90:

        # Instead of wrapping from North to South pole we move to the other side of the North pole
        # (which also involves shifting the longitude by 180 degrees).
        
        lat_corrected = 180 - lat
        
        lon_corrected = 180 + lon_corrected
        if lon_corrected > 180:
            lon_corrected = -360 + lon_corrected

    elif lat < -90:

        # Instead of wrapping from South to North pole we move to the other side of the South pole
        # (which also involves shifting the longitude by 180 degrees).
        
        lat_corrected = -180 - lat
        
        lon_corrected = 180 + lon_corrected
        if lon_corrected > 180:
            lon_corrected = -360 + lon_corrected

    else:

        lat_corrected = lat

    return lat_corrected, lon_corrected


"""

    featureScaling

    Normalise dataset using feature scaling
"""
def featureScaling(data):

    data = np.array(data)
    data = data.astype(float)

    norm_data = []
    min_data = np.nanmin(data)
    max_data = np.nanmax(data)

    for i in xrange(0, len(data)):

        norm_data.append( (data[i] - min_data) / (max_data - min_data) )

    return np.array(norm_data)




""" PALAEOMAGNETICS """
"""
    calcKfromA95

    Calculate koenigsberger ratio (k) from alpha 95 confidence limit (A95) and number of samples (n).
    Returns k
"""
def calcKfromA95(alpha95, n):

    alpha95 = np.radians(alpha95)
    fac = 20.**(1. / (n - 1))
    r2 = 1. / (fac - np.cos(alpha95))
    r2 = r2 * n * (fac - 1.)
    k = (n - 1.) / (n - r2)

    return k
