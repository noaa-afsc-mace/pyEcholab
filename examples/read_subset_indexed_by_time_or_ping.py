# -*- coding: utf-8 -*-
"""Simple example showing how to read subsets of data using the start/end ping
and time arguments as well as reading a random sample of pings from a file.

This example also shows how to use the echosounder.get_rawfile_info() function
which tries to use a .raw file's .idx file to quickly gather some basic information
about what is contained within the .raw file. This can be useful when building
up metadata about a raw file dataset.

"""

import random
import numpy as np
from echolab2.instruments import echosounder

# specify a raw file. While this example should work with any .raw file, the times
# and ping numbers are based on the example data file DY2104-D20210602-T102103.raw.
# If you want to run this with your own file, you will have to adjust the times
# and ping numbers appropriately.
raw_file = './data/EK80/cw/DY2104/raw/DY2104-D20210602-T102103.raw'


# first read a time range.
#
# times in echolab are represented by datetime64 objects with millisecond units.
# datetime64 objects with other units will be cast to millisecond units.
start_time = np.datetime64('2021-06-02T10:22:00', 'ms')
end_time = np.datetime64('2021-06-02T10:25:00', 'ms')

# Read a time span - if the time span is outside the span of your data,
# no data will be read. The resulting EK*0 object will contain ancillary
# data, but the sample data and associated axes data will be empty.
print(f"Reading a range of pings from time {start_time} to time {end_time} inclusive...")
data = echosounder.read(raw_file, start_time=start_time, end_time=end_time)
print(data)


# next, read data based on a range of ping numbers. Ping numbers are relative to
# the data file you are reading and always start at 1. Start and end ping numbers
# are inclusive. In this example we're reading 10 pings, nummbers 150-159.
# Same deal as with time, if your ping range is outside the range of your data
# the EK*0 object will be mostly void of data.
start_ping = 150
end_ping = 159

# Read a ping span
print(f"Reading a range of pings from ping {start_ping} to ping {end_ping} inclusive...")
data = echosounder.read(raw_file, start_ping=start_ping, end_ping=end_ping)
print(data)


# you can also read specific pings by passing arrays of ping numbers or ping times.
# Just like with passing ranges, ping numbers are relative to the file and start
# at 1. Here we will call echosounder.get_rawfile_info() to get some basic information
# including the total number of pings in the file. We'll use the total ping number
# to create the upper limit of the random list of pings to extract. The 
# get_rawfile_info() method returns a lot more info about raw file, and can be
# 
file_info = echosounder.get_rawfile_info(raw_file)

# Create an array of n_pings_to_read length of random pings to read.
n_pings_to_read = 15
# set the possible range of random pings to 1 thru the total pings in the file
ping_range = (1, file_info['n_pings'])
# Create the list of random pings - pings do not need to be in order and
# the pings will be read in the order they were written to the file.
read_pings = [random.randint(*ping_range) for _ in range(n_pings_to_read)]

# Finally, read random pings using the read_pings keyword
print("Reading the following random pings:",read_pings )
data = echosounder.read(raw_file, read_pings=read_pings)
print(data)


"""
Reading a range of pings from time 2021-06-02T10:22:00.000 to time 2021-06-02T10:25:00.000 inclusive...
<class 'echolab2.instruments.EK80.EK80'> at 0x1d29251d5d0
    EK80 object contains data from 5 channels:
        WBT 998500-15 ES18_ES :: power/angle (171, 24197)
        WBT 978217-15 ES38-7_ES :: power/angle (171, 33875)
        WBT 978213-15 ES70-7C_ES :: power/angle (171, 28229)
        WBT 976714-15 ES120-7C_ES :: power/angle (171, 28229)
        WBT 978208-15 ES200-7C_ES :: power/angle (171, 33875)
    data start time: 2021-06-02T10:22:00.535
      data end time: 2021-06-02T10:24:59.960
    number of pings: 171

Reading a range of pings from ping 150 to ping 159 inclusive...
<class 'echolab2.instruments.EK80.EK80'> at 0x1d2925c4950
    EK80 object contains data from 5 channels:
        WBT 998500-15 ES18_ES :: power/angle (10, 24197)
        WBT 978217-15 ES38-7_ES :: power/angle (10, 33875)
        WBT 978213-15 ES70-7C_ES :: power/angle (10, 28229)
        WBT 976714-15 ES120-7C_ES :: power/angle (10, 28229)
        WBT 978208-15 ES200-7C_ES :: power/angle (10, 33875)
    data start time: 2021-06-02T10:23:43.533
      data end time: 2021-06-02T10:23:53.875
    number of pings: 10

Reading the following random pings: [248, 306, 148, 207, 169, 62, 298, 270, 285, 121, 180, 1, 276, 218, 126]
<class 'echolab2.instruments.EK80.EK80'> at 0x1d295638210
    EK80 object contains data from 5 channels:
        WBT 998500-15 ES18_ES :: power/angle (15, 24197)
        WBT 978217-15 ES38-7_ES :: power/angle (15, 33875)
        WBT 978213-15 ES70-7C_ES :: power/angle (15, 28229)
        WBT 976714-15 ES120-7C_ES :: power/angle (15, 28229)
        WBT 978208-15 ES200-7C_ES :: power/angle (15, 33875)
    data start time: 2021-06-02T10:21:03.547
      data end time: 2021-06-02T10:26:20.218
    number of pings: 15
"""
