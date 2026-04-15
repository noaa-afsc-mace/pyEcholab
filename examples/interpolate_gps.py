# -*- coding: utf-8 -*-
"""An example of extracting NMEA data and interpolating it.

The Kongsberg raw file format provides a couple of methods for storing ancillary
sensor data like GPS, speed, heading, pitch, roll, and heave. ASCII NMEA strings
are stored as NMEA datagrams and binary position and motion data can be stored
as MRU datagrams. Regardless of how this data is stored, it is collected
asyncronously so the GPS fixes, heading values etc. are not synced to the ping
data.

When working with transformed data (Sv, Sp, power) it is convienient to interpolate
these data to the ping times and the echosounder.read* functions will do this
for you. But in some cases you may want to interpolate this data for a different
use. This script demonstrates accessing NMEA and MRU data directly and interpolating
data on an arbitrary time axis.

If you simply want Sv data with NMEA/MRU data interpolated to the ping times, you
can use the echosounder module which will do this for you.

"""

import numpy as np
from echolab2.instruments import echosounder


'''
The first example shows how to read a raw file and interpolate NMEA based data.
'''

#  read in a raw file where the position data is stored as NMEA datagrams. To
#  keep things simple, we will only read the 38 kHz data.
rawfiles = ['C:/EK Test Data/EK80/DY2104/cw/DY2104-D20210602-T101508.raw']
print('Reading data...')
ek_data = echosounder.read(rawfiles, frequencies=[38000])
print(ek_data)

'''
<class 'echolab2.instruments.EK80.EK80'> at 0x1f68d2f2f10
    EK80 object contains data from 1 channel:
        WBT 978217-15 ES38-7_ES :: power/angle (352, 33875)
    data start time: 2021-06-02T10:15:08.535
      data end time: 2021-06-02T10:21:03.089
    number of pings: 352
'''

#  NMEA data is asynchronous to the ping data and it applies to all raw data contained
#  within the EK60/EK80 object. It is stored within an instance of the instruments.util.nmea_data
#  class in the EK object's nmea_data attribute. Get the reference to this object and print
#  out some info about it. 
nmea_data = ek_data.nmea_data
print(nmea_data)

'''
<class 'echolab2.instruments.util.nmea_data.nmea_data'> at 0x1f68cb4e110
      NMEA data start time: 2021-06-02T10:15:08.535
        NMEA data end time: 2021-06-02T10:21:03.089
  number of NMEA datagrams: 1961
         unique talker IDs: SD,GP
        unique message IDs: VLW,DTM,GGA,HBT,RMC,VTG,ZDA
'''

#  If you just want to interpolate NMEA data to ping times, then the easiest thing to do
#  is to use the echosounder module's get* methods which will do this for you. But if for
#  some reason you don't want to transform the sample data but still want NMEA data 
#  interpolated to the ping times, you can do it manually.

#  first you need to get the ping times. To do that, you need to get the channel data.
raw_data = ek_data.get_channel_data(frequencies=[38000])

#  The get_channel_data method will return a dict, keyed by the specified attribute, in this
#  case frequency. The value(s) will be a list of raw data objects associated with that
#  attribute. The values are *always* a list, even if you only have a single channel associated
#  with the specified attribute. This may seem overly complex and cumbersome, but there
#  are reasons why this is. For convienience, we'll unpack our raw data object:
raw_data = raw_data[38000][0]


#  Now that we have the raw data which contains the ping_time attribute, use the interpolate
#  method of the nmea_data class to return locations interpolated to the ping times. Since
#  NMEA data can contain a wide variety of data types and formats, you also need to pass the
#  message_type to interpolate. This can be the NMEA sentence type, like "GGA", "HDT", or
#  "VLW" or it can be a metatype as defined in the nmea_data class. The "position" metatype
#  will extract lat/lon data from "GGA", "GLL", or "RMC" datagrams.
ping_interp_fields, ping_interp_data = nmea_data.interpolate(raw_data.ping_time, 'position')
print("Position interpolated to ping times: " + str(len(ping_interp_data['latitude'])) +
        ' lat/lon pairs')
        
#  The interpolate method will return a list containing the field names of the data that
#  were returned, as well as a dictionary keyed by those filed names containing the data
#  for that field. The dictionary also contains the ping_time field which contains the
#  interpolated time vector. Print out a few values.
for i in range(min(5, len(ping_interp_data['ping_time']))):
    print("Time: {0} Lat: {1:.5f} Lon: {2:.5f}".format(ping_interp_data['ping_time'][i],
            ping_interp_data['latitude'][i], ping_interp_data['longitude'][i]))
print()
'''
Position interpolated to ping times: 352 lat/lon pairs
Time: 2021-06-02T10:15:08.535 Lat: 57.35914 Lon: -152.19629
Time: 2021-06-02T10:15:09.533 Lat: 57.35910 Lon: -152.19632
Time: 2021-06-02T10:15:10.534 Lat: 57.35906 Lon: -152.19635
Time: 2021-06-02T10:15:11.533 Lat: 57.35901 Lon: -152.19638
Time: 2021-06-02T10:15:12.534 Lat: 57.35896 Lon: -152.19642
'''

#  For some applications you may simply want to interpolate to an arbitrary time grid.
#  Here we'll create an array of lat/lon pairs every minute over the timespan of the data
arbitrary_times = np.arange(ek_data.start_time, ek_data.end_time,
        np.timedelta64(1, 'm'))
arb_interp_fields, arb_interp_data = nmea_data.interpolate(arbitrary_times, 'position')
print("Position interpolated to 1 minute intervals: " + str(len(arb_interp_data['latitude'])) +
        ' lat/lon pairs')
        
#  print out a few values
for i in range(min(5, len(arb_interp_data['ping_time']))):
    print("Time: {0} Lat: {1:.5f} Lon: {2:.5f}".format(arb_interp_data['ping_time'][i],
            arb_interp_data['latitude'][i], arb_interp_data['longitude'][i]))
print()
'''
Position interpolated to 1 minute intervals: 6 lat/lon pairs
Time: 2021-06-02T10:15:08.535 Lat: 57.35914 Lon: -152.19629
Time: 2021-06-02T10:16:08.535 Lat: 57.35612 Lon: -152.19813
Time: 2021-06-02T10:17:08.535 Lat: 57.35308 Lon: -152.19998
Time: 2021-06-02T10:18:08.535 Lat: 57.35006 Lon: -152.20176
Time: 2021-06-02T10:19:08.535 Lat: 57.34700 Lon: -152.20356
'''

#  You can interpolate any NMEA field that the nmea_data class knows how to parse.
#  Currently that is GGA, GLL, RMC, HDT, and VTG datagrams as well as the metatypes
#  'speed' (requires VTG datagrams), 'attitude' (requires SHR datagrams), 'distance'
#  (requires VLW datagrams), and 'position' (requries GGA and/or RMC and/or GLL)
#  The example data I selected contains VLW datagrams so we'll interpolate distance
#  at 30 second intervals

arbitrary_times = np.arange(ek_data.start_time, ek_data.end_time,
        np.timedelta64(30, 's'))
vlw_interp_fields, vlw_interp_data = nmea_data.interpolate(arbitrary_times, 'distance')
print("Vessel log (distance in nmi) interpolated to 30 second intervals: " +
        str(len(vlw_interp_data['trip_distance_nmi'])) + ' distance values')
        
#  print out a few values
for i in range(min(5, len(arb_interp_data['ping_time']))):
    print("Time: {0} Distance (VLW): {1:.4f}".format(vlw_interp_data['ping_time'][i],
            vlw_interp_data['trip_distance_nmi'][i]))
print()
'''
Vessel log (distance in nmi) interpolated to 30 second intervals: 12 distance values
Time: 2021-06-02T10:15:08.535 Distance (VLW): 46.7910
Time: 2021-06-02T10:15:38.535 Distance (VLW): 46.8876
Time: 2021-06-02T10:16:08.535 Distance (VLW): 46.9835
Time: 2021-06-02T10:16:38.535 Distance (VLW): 47.0798
Time: 2021-06-02T10:17:08.535 Distance (VLW): 47.1763
'''

print()




'''
In the following example, we will use lat/lon data contained within the extended motion
data (MRU1 datagram). This is a less common configuration so you may not have example
data to work with but it is included here to illustrate a
'''

#  read in a raw file where the position data was only stored in the motion datagrams.
print('Reading DriX data...')
rawfiles = ['C:/EK Test Data/EK80/Drix12/Drix12-D20230125-T160206.raw']
drix_data = echosounder.read(rawfiles)
print(drix_data)

'''
<class 'echolab2.instruments.EK80.EK80'> at 0x24d292ab9d0
    EK80 object contains data from 4 channels:
        WBT Tube 279893-15 ES120-7C_ES :: power/angle (4610, 14431)
        WBT Tube 279898-15 ES70-7C_ES :: power/angle (4610, 14431)
        WBT Tube 277104-7 ES38-18|200-18CR_ES :: power/angle (4610, 17317)
        WBT Tube 277104-8 ES38-18|200-18CR_ES :: power (4610, 17317)
    data start time: 2023-01-25T16:02:06.754
      data end time: 2023-01-25T16:48:16.919
    number of pings: 4610
'''

#  there is NMEA data in the file, but it is only the SDVLW datagrams generated by EK80.
#  (and there is only 1 of those in the file.)
drix_nmea_data = drix_data.nmea_data
print(drix_nmea_data)

'''
<class 'echolab2.instruments.util.nmea_data.nmea_data'> at 0x2076a5af310
      NMEA data start time: 2023-01-25T16:02:06.754
        NMEA data end time: 2023-01-25T16:02:06.754
  number of NMEA datagrams: 1
         unique talker IDs: SD
        unique message IDs: VLW
'''


#  now get a reference to the motion data. Like NMEA data, motion data is asynchronous to
#  the ping data and it applies to all pings. You can see below that this system was
#  configured to record motion data at a very high rate compared to the ping rate
#  as there are 96880 motion datagrams for 4610 pings.
drix_motion_data = drix_data.motion_data
print(drix_motion_data)

'''
<class 'echolab2.instruments.util.motion_data.motion_data'> at 0x24d2bcfa010
       MRU data start time: 2023-01-25T16:02:06.754
         MRU data end time: 2023-01-25T16:48:17.508
       Number of datagrams: 96880
  Has extended motion data: True
'''

#  Similar to what we did above, we'll interpolate values on a 1 minute interval. The
#  motion_data class has the same interpolate method as the nmea_data class, though it
#  is not as developed. There are no metatypes so you must specify the fields to interpret
#  and it does not return time.
arbitrary_times = np.arange(drix_data.start_time, drix_data.end_time,
        np.timedelta64(1, 'm'))
arb_interp_fields, arb_interp_data = drix_motion_data.interpolate(arbitrary_times,
        attributes=['latitude', 'longitude'])
print("Drix Position interpolated to 1 minute intervals: " + str(len(arb_interp_data['latitude'])) +
        ' lat/lon pairs')

#  print out a few values
for i in range(min(5, len(arbitrary_times))):
    print("Time: {0} Lat: {1:.5f} Lon: {2:.5f}".format(arbitrary_times[i],
            arb_interp_data['latitude'][i], arb_interp_data['longitude'][i]))
print()
'''
Drix Position interpolated to 1 minute intervals: 47 lat/lon pairs
Time: 2023-01-25T16:02:06.754 Lat: 47.71860 Lon: -122.40583
Time: 2023-01-25T16:03:06.754 Lat: 47.71858 Lon: -122.40584
Time: 2023-01-25T16:04:06.754 Lat: 47.71859 Lon: -122.40584
Time: 2023-01-25T16:05:06.754 Lat: 47.71859 Lon: -122.40583
Time: 2023-01-25T16:06:06.754 Lat: 47.71860 Lon: -122.40582
'''

print()
