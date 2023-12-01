# -*- coding: utf-8 -*-
"""
@author: rick.towler

This example script plots a single ping as Sv for every channel in the
specified raw file.

"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
from echolab2.instruments import EK80


# Specify the ping number to plot.
ping_number = 24

# Define the path to the data file.
raw_filename = 'C:/EK Test Data/EK80/CW/complex/DY1802_EK80-D20180301-T185940.raw'

# Create an instance of our EK80 object. The EK80 object will contain
# all of the data elements read from the raw data files. This includes
# system configuration data, raw acoustic data, environmental data,
# and NMEA data (GPS, etc.)
ek80 = EK80.EK80()

# Read in the .raw data file. you can pass in a string defining the path
# to the data file, or a list of strings. If you provide a list of strings,
# each list element will be read and added to the EK80 object.
print('Reading raw file %s' % (raw_filename))
ek80.read_raw(raw_filename)

# Print some info about the state of our EK80 object.
print(ek80)

'''
<class 'echolab2.instruments.EK80.EK80'> at 0x1f28aa108b0
    EK80 object contains data from 5 channels:
        WBT 545612-15 ES18 :: complex-CW (88, 48723, 4)
        WBT 549762-15 ES38B :: complex-CW (88, 34106, 4)
        WBT 582207-15 ES70-7C :: complex-CW (88, 28422, 4)
        WBT 582214-15 ES120-7C :: complex-CW (88, 34106, 4)
        WBT 582215-15 ES200-7C :: complex-CW (88, 42633, 4)
    data start time: 2018-03-01T18:59:40.754
      data end time: 2018-03-01T19:06:39.186
    number of pings: 88
'''

# You can see that there is data from 5 different channels and
# the channel id, data type, and raw data array size is shown
# for each channel. For example:
#
#   WBT 545612-15 ES18 :: complex-CW (88, 48723, 4)
# 
# In this case, channel id "WBT 545612-15 ES18" has 88 pings of 
# complex CW data. Each ping has 48723 samples and since this
# is complex data, the array has a 3rd dimension where each
# sample has 4 sectors.
#
# The type of data stored in the file depends on how the EK80
# system was configured when the data were recorded. Data types
# can be complex-CW, complex-FM, or the "reduced" types power+angle,
# power, and angle.

# The biggest reason for the steep learning curve is the lack of
# proper documentation. You're going to have to get comfortable
# reading the method headers and using dir() to print out class
# attributes. Note that I have edited out all of the "private"
# attributes and methods (anything with a "_" in the name.)

'''
>>> dir(ek80)
[ 'annotations', 'append_raw', 'channel_ids', 'channel_number_map',
'end_ping', 'end_time', 'frequency_map', 'get_channel_data',
'motion_data', 'n_channels', 'n_files', 'n_pings', 'nmea_data',
'raw_data', 'raw_data_width', 'raw_tx', 'read_bot', 'read_channel_ids',
'read_end_ping', 'read_end_sample', 'read_end_time', 'read_frequencies',
'read_idx', 'read_max_sample_count', 'read_raw', 'read_start_ping',
'read_start_sample', 'read_start_time', 'start_ping', 'start_time',
'store_angles', 'store_complex', 'store_power', 'write_raw']
>>>
'''

# You can't always tell what is an attribute and what is a method but
# you can just print the attribute/method to figure that out:

'''
>>> print(ek80.n_channels)
5
>>> print(ek80.write_raw)
<bound method EK80.write_raw of <echolab2.instruments.EK80.EK80 object at 0x000001F28AA108B0>>
>>> 
'''

# n_channels is an attribute with a value of 5 where write_raw is a
# method. The most commonly used attributes are typically:
#
# channel_ids - a list containing the channel id's of the data read
# raw_data - the raw acoustic data, stored by channel id and data type.
# motion_data - vessel motion data stored in the raw file.
# nmea_data - the NMEA data stored in the raw file data.


# Set up the colormap for plotting.
color = iter(cm.rainbow(np.linspace(0, 1, ek80.n_channels)))

# Create a matplotlib figure to plot our echograms on.
fig = plt.figure(figsize=(7, 7))


# As noted above, the raw_data attribute stores the raw acoustic
# data by channel id. This attribute is a dict keyed by channel
# id. Here we will iterate thru the channel ids and plot a ping
# from each channel.

for channel_id in ek80.channel_ids:

    # Get a color for this channel.
    c = next(color)

    # raw data is stored in the EK80.raw_data attribute as a dict
    # keyed by channel id. Get the raw data for this channel id.
    raw_data = ek80.raw_data[channel_id]

    # raw data for a channel is stored in a list by data type.
    # If you read multiple files that contain different data types
    # there will be multiple elements in this list. If you read
    # only one file, or files with only one data type, there will
    # be only one element in this list.
    
    # In this case we have read a single file so we know we have one
    # raw_data object in this list so we grab the first element.
    raw_data = raw_data[0]

    # we can print some basic info about our raw data
    print(raw_data)

    # This of course will be different for each channel, but the
    # first channel in this file is the 18 kHz:
    '''
    <class 'echolab2.instruments.EK80.raw_data'> at 0x25bd2391ca0
                       channel: WBT 545612-15 ES18
        frequency (first ping): 18000.0
       pulse length (first ping): 0.001024
                 data start time: 2018-03-01T18:59:40.754
                   data end time: 2018-03-01T19:06:33.559
                 number of pings: 88
                      data type: CW complex
       complex array dimensions: (88,48723,4)

    >>> 
    '''

    # Just like with the ek80 object, you can use dir() to get more
    # info about the attributes and methods of the raw_data object.
    #dir(raw_data)

    # What one would do next depends on what they are trying to do,
    # but most people will want to convert the raw acoustic data to
    # Sv or Sp (aka TS). This process uses the system's calibration
    # parameters and raw data and converts it to a more usable form.
    
    # Simrad .raw data files store the system calibration parameters
    # configured at the time of recording and we can extract these
    # from the raw data.
    calibration = raw_data.get_calibration()

    # you can print out info about the calibration params. Doing so
    # will show each of the object's attributes and a sample value.
    # I encourage you to look at this but am commenting it out
    # here since it prints out a lot of stuff.
    #print(calibration)
    
    # Depending on your organization's data collection and processing
    # procedures, you may need to alter the calibration parameters
    # before converting the data. Some organizations will load 
    # calibration data into the EK80 before collecting data while others
    # always collect data with the system default values and apply
    # the cal values in post processing. If you are working with data
    # that doesn't have the calibration parameters loaded, you would
    # want to change the appropriate attributes of the calibration
    # object here.
    #
    # In this example we're not doing anything quantitative so we really
    # don't care and just take whatever params are in the .raw file.

    # next we will convert the raw acoustic data in our raw_data object
    # to Sv by calling the get_Sv() method, passing the calibration object.
    Sv_data = raw_data.get_Sv(calibration=calibration)


    # get_Sv returns data converted to Sv in a processed_data object. A
    # processed_data object stores the converted/processed/calibrated
    # data in a numpy array along with the axes data and optionally
    # other data associated with those axis. You can use print to get
    # some basic info about the data in the processed_data object
    print(Sv_data)
    '''
    <class 'echolab2.processing.processed_data.processed_data'> at 0x1eedcf58e80
                    channel(s): [WBT 545612-15 ES18]
                frequency (Hz): 18000.0
               data start time: 2018-03-01T18:59:40.754
                 data end time: 2018-03-01T19:06:33.559
               number of pings: 88
               data attributes: ping_time (88)
                                data (88,48723)
                                range (48723)
                                transducer_offset (88)
    '''

    # And now I think we are finally getting to closer to your answer?
    
    # We now have an object that contains Sv data from the 18 kHz that
    # is in an array 88 pings by 48723 samples. The ping_time attribute
    # (think x axis) stores the time (in GMT) of that ping and the range
    # attribute (think y axis) stores the range in meters from the
    # transducer face for each sample.
    
    # The attributes shown are the literal attribute names. For example,
    # to access the individual ping times, you reference the ping_time
    # attribute:
    '''
    >>> Sv_data.ping_time
    array(['2018-03-01T18:59:40.754', '2018-03-01T18:59:45.499',
           '2018-03-01T18:59:50.243', '2018-03-01T18:59:54.988',
           '2018-03-01T18:59:59.734', '2018-03-01T19:00:04.479',
           '2018-03-01T19:00:09.224', '2018-03-01T19:00:13.968',
           '2018-03-01T19:00:18.714', '2018-03-01T19:00:23.459',
           ...
           '2018-03-01T19:06:00.344', '2018-03-01T19:06:05.089',
           '2018-03-01T19:06:09.834', '2018-03-01T19:06:14.579',
           '2018-03-01T19:06:19.324', '2018-03-01T19:06:24.069',
           '2018-03-01T19:06:28.814', '2018-03-01T19:06:33.559'],
          dtype='datetime64[ms]')
    '''
    
    # range, contains the range in meters.
    '''
    >>>Sv_data.range
    array([1.00000000e-20, 2.05240003e-02, 4.10480006e-02, ...,
           9.99929294e+02, 9.99949818e+02, 9.99970342e+02])
    '''
    # We see here the first sample is effectively at zero meters
    # and the last sample is basically 1000 meters.
    
    # The data attribute contains the Sv data. To print the data
    # from the ->third<- ping
    '''
    >>> Sv_data.data[2,:]
    array([-422.17894215, -416.13466114, -410.98253101, ...,  -87.5875078 ,
            -87.64481624,  -87.7775565 ])
    '''
    
    # Note that you can access the same data by slicing the
    # processed_data object directly:
    '''
    >>> Sv_data[2,:]
    array([-422.17894215, -416.13466114, -410.98253101, ...,  -87.5875078 ,
            -87.64481624,  -87.7775565 ])
    '''

    # Again, dir() is your friend here. Look at the attributes and methods of
    # the processed_data object and don't hesitate to open the processed_data
    # class and read the docs in the method headers.

    # Plot Sv of a single ping for this channel.
    plt.plot(Sv_data[ping_number,:], Sv_data.range, color=c, label=channel_id)
    
# Label the figure and set other display properties.
plt.gca().invert_yaxis()
plt.ylabel('Range (m)')
plt.xlabel('Sv (dB)')
title = 'Ping %i' % (ping_number)
fig.suptitle(title, fontsize=14)
plt.legend()

# Display plot.
plt.show()
