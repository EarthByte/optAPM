""" Library containing all the required functions to calculate objective function data evaluations
    for optimisation of Africa reference frame

    Author: Simon Williams / Michael Tetley
            Earthbyte

    Date:   02/11/2015

"""

import math
import pygplates
import numpy as np


# Fracture zone orientation functions
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



def Calc_Median(rotation_model, PID, seafloor_ages, Lats, Lons, spreading_directions):

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



# Net rotation functions
def DiffPoles(poleA_lon_degrees,poleA_lat_degrees,angleA_degrees,poleB_lon_degrees,poleB_lat_degrees,angleB_degrees):

    angleA_radians = angleA_degrees*(math.pi/180)
    poleA_pos = pygplates.convert_lat_lon_point_to_point_on_sphere(pygplates.LatLonPoint(poleA_lat_degrees, poleA_lon_degrees))
    stage_poleA = pygplates.FiniteRotation(poleA_pos, angleA_radians)

    angleB_radians = angleB_degrees*(math.pi/180)
    poleB_pos = pygplates.convert_lat_lon_point_to_point_on_sphere(pygplates.LatLonPoint(poleB_lat_degrees, poleB_lon_degrees))
    stage_poleB = pygplates.FiniteRotation(poleB_pos, angleB_radians)

    diff_rotation = stage_poleA.get_inverse() * stage_poleB

    diff_lat_degrees, diff_lon_degrees, diff_angle_degrees = diff_rotation.get_lat_lon_euler_pole_and_angle_degrees()

    # threshold=0.001
    # if math.fabs(diff_angle_degrees) < threshold:
    #     print 'Stage poles almost the same'

    return diff_lon_degrees, diff_lat_degrees, diff_angle_degrees



def GetStagePoleSequenceForPlate(rotation_model,MovPlate,AnchorPlate,TIMELIST):
    LongList=[]
    LatList=[]
    angleList=[]

    count=0
    for TIME in TIMELIST:
        stageData = rotation_model.get_rotation(TIME, MovPlate, TIME+1, anchor_plate_id=AnchorPlate)
        pole,angle = stageData.get_euler_pole_and_angle()
        pole = pygplates.convert_point_on_sphere_to_lat_lon_point(pole)

        LongList.append(pole.get_longitude())
        LatList.append(pole.get_latitude())
        angleList.append(angle*180/math.pi)

        count = count+1

    return LongList,LatList,angleList


def ApproximateNR_from_features(rotation_model,rotation_model_NNR,TIMELIST,MovPlate):

    # For both rotation models, get the Stage Poles for the Moving plate relative to the spin axis
    SPLong_NNR,SPLat_NNR,SPangle_NNR = GetStagePoleSequenceForPlate(rotation_model_NNR,MovPlate,0,TIMELIST)
    SPLong,SPLat,SPangle = GetStagePoleSequenceForPlate(rotation_model,MovPlate,0,TIMELIST)

    NRLong = []
    NRLat = []
    NRangle = []

    # Iterate over each time step, get the Euler pole that is the difference between
    # the two rotation model stage poles
    for index in range(0,len(SPLong)):

        tmp1,tmp2,tmp3 = DiffPoles(SPLong[index],SPLat[index],SPangle[index],SPLong_NNR[index],SPLat_NNR[index],SPangle_NNR[index])
        NRLong.append(tmp1)
        NRLat.append(tmp2)
        NRangle.append(tmp3)

    return NRLong,NRLat,NRangle,SPLong,SPLat,SPangle,SPLong_NNR,SPLat_NNR,SPangle_NNR



def ApproximateNR(input_rotation_filename,input_rotation_filename_NNR,TIMELIST,MovPlate):

    # If the inputs are rotation files, call this function.
    # Alternatively, if the two rotation models (the model to be tested and the NNR model)
    # are already opened within memory within the calling python script, the function
    # 'ApproximateNR_from_features' should be called directly

    # NB in majority of cases, MovPlate will be 701

    # Get rotation models for both the file containing the reference frame to test, and the
    # 'baseline' rotation model with NNR
    rotation_model_NNR = pygplates.RotationModel(input_rotation_filename_NNR)
    rotation_model = pygplates.RotationModel(input_rotation_filename)

    # Call the function that gets the net rotation angle
    NRLong,NRLat,NRangle,SPLong,SPLat,SPangle,SPLong_NNR,SPLat_NNR,SPangle_NNR = \
        ApproximateNR_from_features(rotation_model,rotation_model_NNR,TIMELIST,MovPlate)


    return NRLong,NRLat,NRangle,SPLong,SPLat,SPangle,SPLong_NNR,SPLat_NNR,SPangle_NNR
