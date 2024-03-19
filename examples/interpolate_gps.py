# -*- coding: utf-8 -*-
"""An example of extracting NMEA data and interpolating it.

This script demonstrates extracting NMEA lat/lon data and then interpolating it
a couple of different ways

THIS EXAMPLE IS A WORK IN PROGRESS AND IS INCOMPLETE


"""

import numpy as np
from echolab2.instruments import echosounder


'''
show examples of interpolating when working with processed_data, and when not working with processed_data

'''


#  read in a raw file where the position data is stored as NMEA datagrams.
print('Reading Dyson data...')
rawfiles = ['//nmfs.local/akc-race/MACE_Acoustic2/DY2308/EK80/DY2308_EK80-D20230717-T012621.raw']
dyson_data = echosounder.read(rawfiles)
print('Done.')

#  we just read a single EK80 file, so we know there will be a single EK80 object in
#  the list returned by echosounder.read() so we'll unpack that here and print some
#  basic info about the data in the EK80 object.
dyson_data = dyson_data[0]
print(dyson_data)

#  now get a reference to the NMEA data and print out some info about it. NMEA data
#  is asynchronous to the ping data and it applies to all pings so it is an attribute of
#  the EK80 object.
dyson_nmea_data = dyson_data.nmea_data
print(dyson_nmea_data)


#  get the 38 kHz Dyson data
dyson_raw_data = dyson_data.get_channel_data(frequencies=[38000])
dyson_raw_data = dyson_raw_data[38000][0]

#  If you're working with converted sample data (e.g. Sv) and you want to interpolate GPS
#  positions to ping times, the easiest thing to do is to get your processed_data object
#  and then call the apply_navigation() method, passing the nmea_data object. This will
#  interpolate the navigation data (including lat/lon) to the processed data ping times
#  adding the attributes to your processed data object.


#  But you may not always want/need to work with the sample data. In these cases, you can
#  call the nmea_data.interpolate() method directly, passing the time vector you want to
#  interpolate your nmea data to. Since NMEA data can contain a wide variety of data types
#  and formats, you also need to pass the message_type to interpolate. This can be the NMEA
#  sentence type, like "GGA", "HDT", or "VLW" or it can be a metatype as defined in the
#  nmea_data class. The "position" metatype will extract lat/lon data from "GGA", "GLL", or
#  "RMC" datagrams.
dy_ping_interp_fields, dy_ping_interp_data = dyson_nmea_data.interpolate(dyson_raw_data.ping_time, 'position')
print("Dyson position interpolated to ping times: " + str(len(dy_ping_interp_data['latitude'])) +
        ' lat/lon pairs')

#  For some applications you may simply want to interpolate to an arbitrary time grid.
#  Here we'll create an array of lat/lon pairs every minute over the timespan of the data
arbitrary_times = np.arange(dyson_data.start_time, dyson_data.end_time,
        np.timedelta64(1, 'm'))
arb_interp_fields, dy_arb_interp_data = dyson_nmea_data.interpolate(arbitrary_times, 'position')
print("Dyson Position interpolated to 1 minute intervals: " + str(len(dy_arb_interp_data['latitude'])) +
        ' lat/lon pairs')


'''
In the following example, we will use lat/lon data contained within the extended motion
data (MRU1 datagram).
'''



#  read in a raw file where the position data was only stored in the motion datagrams.
#  this file contains no NMEA data, but the extended motion datagram (MRU1 datagram)
#  contains lat/lon data we can use.
print('Reading DriX data...')
rawfiles = ['//nmfs.local/akc-race/MACE_Acoustic2/DY2308/DY2308_DriX_backup/Data_from_DRiX/Drix_survey_PC/EK80/DriX12_DY2308-D20230717-T012842.raw',
           '//nmfs.local/akc-race/MACE_Acoustic2/DY2308/DY2308_DriX_backup/Data_from_DRiX/Drix_survey_PC/EK80/DriX12_DY2308-D20230717-T013847.raw',
           '//nmfs.local/akc-race/MACE_Acoustic2/DY2308/DY2308_DriX_backup/Data_from_DRiX/Drix_survey_PC/EK80/DriX12_DY2308-D20230717-T014852.raw']
drix_data = echosounder.read(rawfiles)
print('Done.')

#  unpack our EK80 object and print some basic info about it
drix_data = drix_data[0]
print(drix_data)

#  now get a reference to the motion data. Like NMEA data, motion data is asynchronous to
#  the ping data and it applies to all pings.
drix_motion_data = drix_data.motion_data
print(drix_motion_data)


#  get the 38 kHz DriX data
drix_raw_data = drix_data.get_channel_data(frequencies=[38000])
drix_raw_data = drix_raw_data[38000][0]


#  and then call the motion_data interpolate method. We will specify that it only
#  interpolates and returns lat/lon (by default it will also include heave, pitch,
#  and roll.) The interpolate method will return a list containing the key names
#  of the fields that were interpolated as well as a dict containing the data.
ping_interp_fields, drix_ping_interp_data = drix_motion_data.interpolate(drix_raw_data.ping_time,
        attributes=['latitude', 'longitude'])
print("Drix Position interpolated to ping times: " + str(len(drix_ping_interp_data['latitude'])) +
        ' lat/lon pairs')

#  For some applications you may simply want to interpolate to an arbitrary time grid.
#  Here we'll create an array of lat/lon pairs every minute over the timespan of the data
arbitrary_times = np.arange(drix_data.start_time, drix_data.end_time,
        np.timedelta64(1, 'm'))
arb_interp_fields, drix_arb_interp_data = drix_motion_data.interpolate(arbitrary_times,
        attributes=['latitude', 'longitude'])
print("Drix Position interpolated to 1 minute intervals: " + str(len(drix_arb_interp_data['latitude'])) +
        ' lat/lon pairs')


pass
