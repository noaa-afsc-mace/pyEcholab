# -*- coding: utf-8 -*-
"""

"""
import os
import datetime
import numpy
from echolab2.instruments import EK80



use_file_date_for_zda = True

raw_dir = '//akc0ss-n086/RACE_GF_Acoustic/2025/GOA_2025_Alaska_Provider/'
rawfiles = ['akp-D20250531-T232902.raw',
            'akp-D20250601-T030347.raw',
            'akp-D20250604-T144109.raw',
            'akp-D20250604-T174412.raw']

raw_dir = '//akc0ss-n086/RACE_GF_Acoustic/2025/EBS_2025_NW_Explorer/'
rawfiles = ['NWE-D20250528-T075810.raw',
            'NWE-D20250529-T011153.raw',
            'NWE-D20250607-T012410.raw',
            'NWE-D20250610-T142534.raw']


for rawfile in rawfiles:

    ek80 = EK80.EK80()
    print("Reading " + rawfile)
    ek80.read_raw(raw_dir + rawfile)

    print(ek80)
    print()

    #  get the files modification time (LOCAL time last data was written)
    file_mod_time = os.path.getmtime(raw_dir + rawfile)
    file_timestamp = numpy.datetime64(datetime.datetime.fromtimestamp(file_mod_time))

    #  get the file's datetime - local time converted to UTC based on the PC clock and time zone
    file_datetime = rawfile[-21:-4]
    file_date = file_datetime[1:9]
    file_time = file_datetime[11:]
    file_dt64 = numpy.datetime64(datetime.datetime.strptime(file_date + " " + file_time, "%Y%m%d %H%M%S"))
    print("File time: ", file_dt64)
    print("File modification time: ", file_timestamp)

    #  get the raw NMEA datagrams
    raw_datagrams = ek80.nmea_data.get_datagrams(['ZDA', 'GGA'], return_raw=True)

    #  get the Simrad timestamp for the first ZDA datagram
    first_zda_datagram = raw_datagrams['ZDA']['time'][0]

    #  now get a datetime64 object representing the first ZDA time - this will always be in GMT
    first_zda = raw_datagrams['ZDA']['raw_string'][0]
    first_zda_parsed = first_zda.split(',')
    if use_file_date_for_zda:
        zda_date = file_date
    else:
        zda_date = first_zda_parsed[4] + first_zda_parsed[3] + first_zda_parsed[2]
    first_zda_time_str = (zda_date + " " + first_zda_parsed[1][:-3])
    first_zda_time = numpy.datetime64(datetime.datetime.strptime(first_zda_time_str, "%Y%m%d %H%M%S"))
    last_zda_datagram = raw_datagrams['ZDA']['time'][-1]


    print("First ZDA ES80 time: ", first_zda_datagram)
    print("First ZDA NMEA datagram (time is always in GMT): ", first_zda_time)
    data_dt = first_zda_time - first_zda_datagram
    #  print the difference between ACTUAL UTC and what the ES80 thinks in UTC based on the
    #  local timezone and time
    print("Time difference (hours) between UTC and ES80 data: ", data_dt / numpy.timedelta64(1,'h'))

    print("Time difference (hours) between UTC and PC: ", (file_timestamp-last_zda_datagram)/ numpy.timedelta64(1,'h'))

    print()
