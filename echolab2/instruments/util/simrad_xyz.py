# coding=utf-8

#     National Oceanic and Atmospheric Administration (NOAA)
#     Alaskan Fisheries Science Center (AFSC)
#     Resource Assessment and Conservation Engineering (RACE)
#     Midwater Assessment and Conservation Engineering (MACE)

# THIS SOFTWARE AND ITS DOCUMENTATION ARE CONSIDERED TO BE IN THE PUBLIC DOMAIN
# AND THUS ARE AVAILABLE FOR UNRESTRICTED PUBLIC USE. THEY ARE FURNISHED "AS
# IS. THE AUTHORS, THE UNITED STATES GOVERNMENT, ITS INSTRUMENTALITIES,
# OFFICERS, EMPLOYEES, AND AGENTS MAKE NO WARRANTY, EXPRESS OR IMPLIED,
# AS TO THE USEFULNESS OF THE SOFTWARE AND DOCUMENTATION FOR ANY PURPOSE.
# THEY ASSUME NO RESPONSIBILITY (1) FOR THE USE OF THE SOFTWARE AND
# DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL SUPPORT TO USERS.

"""

    :synopsis:  Reads and writes Simrad/Kongsberg XYZ bottom detection data
                files.


| Developed by:  Rick Towler   <rick.towler@noaa.gov>
| National Oceanic and Atmospheric Administration (NOAA)
| Alaska Fisheries Science Center (AFSC)
| Midwater Assessment and Conservation Engineering Group (MACE)
|
| Author:
|       Rick Towler   <rick.towler@noaa.gov>
| Maintained by:
|       Rick Towler   <rick.towler@noaa.gov>

"""


import os
from datetime import datetime
import numpy as np



def read_xyz(xyz_filename, as_range=False):
    '''read_xyz will read a Simrad EK80 .xyz file and return the data in a
    dictionary.

    xyz_filename (string): The full path to the .xyz file to read
    as_range (bool): Set to True to return the detected bottom value as range
        to bottom. The recorded transducer draft will be subtracted from
        the bottom detection value.

    '''

    def convert_float(val):
        try:
            val = float(val)
        except:
            val = np.nan
        return val


    bottom_data = {}

    # Normalize filename and read the file
    xyz_filename = os.path.normpath(xyz_filename)
    with open(xyz_filename, 'r') as infile:
        xyz_data = infile.readlines()

    #  determine the number of line vertices
    n_points = len(xyz_data)

    # Simrad .xyz files contain Lat, Lon, depth, date, time, and transducer draft
    bottom_data['detected_bottom'] = np.empty((n_points), dtype=np.float32)
    bottom_data['latitude'] = np.empty((n_points), dtype=np.float32)
    bottom_data['longitude'] = np.empty((n_points), dtype=np.float32)
    bottom_data['transducer_draft'] = np.empty((n_points), dtype=np.float32)
    bottom_data['ping_time'] = np.empty((n_points), dtype='datetime64[ms]')

    # Loop thru the rows of data, parsing each line
    for idx, row in enumerate(xyz_data):
        #  split the row
        parts = row.split()
        n_parts = len(parts)

        if n_parts == 8:
            #  this is the XYZ format introduced in EK80 21.15.x with hemisphere
            (lat, lat_h, lon, lon_h, depth, date, time, draft) = parts

            #  convert lat/lon to floats
            lat = convert_float(lat)
            lon = convert_float(lon)

            #  add the sign to the lat/lon
            lat *= 1 if lat_h == 'N' else -1
            lon *= 1 if lon_h == 'E' else -1

        elif n_parts == 6:
            #  this is the OG XYZ with signed lat/lon
            (lat, lon, depth, date, time, draft) = parts

            #  convert lat/lon to floats
            lat = convert_float(lat)
            lon = convert_float(lon)

        else:
            #  this XYZ file is not what we were expecting
            raise Exception('Unknown XYZ format with %d fields' % n_parts)

        # Convert the time elements to datetime64
        bottom_data['ping_time'][idx] = np.datetime64(datetime.strptime(date + time, "%d%m%Y%H%M%S.%f"))

        # Convert depth and draft to floats
        bottom_data['transducer_draft'][idx] = convert_float(draft)
        if as_range:
            bottom_data['detected_bottom'][idx] = convert_float(depth) - bottom_data['transducer_draft'][idx]
        else:
            bottom_data['detected_bottom'][idx] = convert_float(depth)

        # Store lat and lon
        bottom_data['latitude'][idx] = lat
        bottom_data['longitude'][idx] = lon

    return bottom_data



def write_xyz(xyz_data, xyz_filename):
    '''write_xyz will write a Simrad EK80 .xyz file using the data provided in
    xyz_data dictionary.

    xyz_filename (string): The full path to the .xyz file to read
    xyz_data (dict): A dictionary containing the bottom detections, ping times,
        latitude, longitude and transducer draft data to write to the file. The
        dictionary must contain the following fields:

            xyz_data['detected_bottom'] - bottom detections as float
            xyz_data['latitude'] - latitudes as float
            xyz_data['longitude'] - longitudes as float
            xyz_data['transducer_draft'] - transducer drafts as float
            xyz_data['ping_time'] - ping times as datetime64[ms]

        All arrays must be the same length

    '''
    #  write_xyz has not been implemented
    pass

