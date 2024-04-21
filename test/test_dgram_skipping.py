# -*- coding: utf-8 -*-
"""

"""
import numpy as np
from echolab2.instruments import echosounder

# example of a MACE sequence file
#rawfile='//nmfs.local/akc-race/MACE_Acoustic2/DY2403/ek80/DY2403-D20240306-T164955.raw'

# EK60 example with bot
rawfile='C:/EK Test Data/EK60/DY1807/DY1807_EK60-D20180607-T223821.raw'

#  read a "normal" Dyson EK80 file
#rawfile='//nmfs.local/akc-race/MACE_Acoustic2/DY2308/EK80/DY2308_EK80-D20230702-T064241.raw'

#  read a group of Drix files - No NMEA, no bottom, Nav data in MRU1
#rawfile=['//nmfs.local/akc-race/MACE_Acoustic2/DY2308/DY2308_DriX_backup/Data_from_DRiX/Drix_survey_PC/EK80/DriX12_DY2308-D20230707-T214546.raw',
#         '//nmfs.local/akc-race/MACE_Acoustic2/DY2308/DY2308_DriX_backup/Data_from_DRiX/Drix_survey_PC/EK80/DriX12_DY2308-D20230707-T215551.raw',
#         '//nmfs.local/akc-race/MACE_Acoustic2/DY2308/DY2308_DriX_backup/Data_from_DRiX/Drix_survey_PC/EK80/DriX12_DY2308-D20230707-T220556.raw']



rawfile='C:/EK Test Data/EK80/CW/DY2207/DY2207_EK80-D20220610-T055930.raw'

print("Reading .raw data: " + rawfile)

#  read all 38 and 120 kHz data and skip the NMEA data
ek_data = echosounder.read(rawfile, frequencies=[38000,120000], nmea=False)
print("Raw data object info:")
print(ek_data)

##  read a single channel by channel ID
#ek_data = echosounder.read(rawfile, channel_ids=['GPT  38 kHz 009072033fa2 2-1 ES38B'])
#print("Raw data object info:")
#print(ek_data)

#  read a subset of the data by ping range
ek_data = echosounder.read(rawfile, start_ping=100, end_ping=199)
print("Raw data object info:")
print(ek_data)

##  read a subset of the data by time range
#start_time = np.datetime64('2018-06-07T22:40')
#end_time = np.datetime64('2018-06-07T22:44')
#ek_data = echosounder.read(rawfile, start_time=start_time, end_time=end_time,
#        frequencies=[18000])
#print("Raw data object info:")
#print(ek_data)


cal_objects = echosounder.get_calibration_from_raw(ek_data)
Sv_data = echosounder.get_Sv(ek_data, calibration=cal_objects, return_depth=False)


print()


