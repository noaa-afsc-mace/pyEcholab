# -*- coding: utf-8 -*-
"""This example script demonstrates the use of numeric and boolean operators
on processed_data and Mask objects. It also provides an example of using the
processed_data.zeros_like() method to get an processed_data array we can use
to fill with the results of our simple analysis.  Lastly, it shows how to use
the processed_data.view() method to plot a subset of the data.

Note that this example is not intended to be an example of how to really do
frequency differencing, just the techniques needed to implement it.
"""

from matplotlib.pyplot import figure, show, subplots_adjust
from echolab2.instruments import echosounder
from echolab2.plotting.matplotlib import echogram
from echolab2.processing import mask, line
import numpy as np


# create a list of the files we want to read. The bottom detection files
# are assumed to be co-located with the data files.
#
# NOTE! - This example assumes that the data files will contain data at
#         18, 38, and 120 kHz and that bottom data will be available.
#
rawfiles = ['C:/EK Test Data/EK60/DY1807/raw/DY1807_EK60-D20180609-T123907.raw',
            'C:/EK Test Data/EK60/DY1807/raw/DY1807_EK60-D20180609-T124557.raw']

# read the data - this will return a dictionary keyed by channel ID,
# containing the raw data objects with data for that channel. Bottom data
# files will also be read if available.
print("Reading .raw data...")
raw_data = echosounder.read(rawfiles, frequencies=[18000, 38000, 120000])
print(raw_data)

# create calibration objects from the raw data. This is not strictly
# necessary in this case because if you do not provide calibration objects
# when calling the get_* methods, the values from the raw file will be
# used by default. If you were doing this for real, you would most likely
# read cal data from EK80 calibration XML files or an Echoview .ecs file
# using the appropriate echosounder class methods. 
cal_objects = echosounder.get_calibration_from_raw(raw_data)

# transform the raw data into Sv.
Sv_data = echosounder.get_Sv(raw_data, calibration=cal_objects)

# Now create a mask for each channel and apply surface and bottom lines
# to these masks such that we mask out samples near the surface and below the
# bottom. In this example, we will mask everything above 10m and below 0.5m
# above the detected bottom.

masks = {}
freq_chan_map = {}
for chan in Sv_data.keys():

    # Create an empty mask for this channel that matches our data array shape
    masks[chan] = mask.mask(like=Sv_data[chan])

    # Next create a new line based on this channel's bottom line that is
    # 0.5m shallower. This is our bottom exclusion line.
    bot_exclusion_line = Sv_data[chan].bottom_line - 0.5

    # Now create a surface exclusion line at 10m RANGE. We use the line.like()
    # constructor to return a line with ping times that match this channel
    surf_exclusion_line = line.like(Sv_data[chan], data=10)

    # Now apply our lines to our mask. When applying a line to a mask, you
    # will either set all samples that are above the line by setting
    # apply_above=True, or you will set all samples that are below the line
    # by setting apply_above=False. If you do not specify the value, the
    # mask samples will be set to True by default.
    #
    # We are creating an EXCLUSION mask, so samples we want to ignore will
    # be set to True.

    # Set all samples below our bottom exclusion line to True
    masks[chan].apply_line(bot_exclusion_line, apply_above=False)

    # And set all sample above our surface exclusion line to True
    masks[chan].apply_line(surf_exclusion_line, apply_above=True)

    # Now use this mask to set the masked sample data to NaN
    Sv_data[chan][masks[chan]] = np.nan

    # While the channel ID is great for uniquely keying data, in this example
    # we need to key the data by frequency. Since we're already looping
    # thru the channels building the masks, we will create a dictionary to 
    # map channel ID to frequency to make the differencing step below easier.
    freq_chan_map[Sv_data[chan].frequency] = chan

# Now lets compute some differences - the processed_data class implements the
# basic Python arithmetic operators so we can simply subtract processed_data
# objects like numeric objects.  Both, "regular" (+, -, *, /)  and 
# "in-place" (+=, -=, *=, /=) operators are implemented.  Regular operators
# return a new processed_data object with the same general properties containing
# the results of your operation.  The in-place operators will alter the data of
# the left hand side argument.

# 18 - 38
Sv_18m38 = Sv_data[freq_chan_map[18000]] - Sv_data[freq_chan_map[38000]]

# 120 - 38
Sv_120m38 = Sv_data[freq_chan_map[120000]] - Sv_data[freq_chan_map[38000]]


# Now we'll generate some masks identifying samples that fall within various
# ranges.
#
# The processed_data object also implements the Python comparison operators.
# These operators do an element by element comparison and will return a Mask
# object with samples set to the result of the comparison.

# For example, this operation will return a mask object where samples in the
# Sv_18m38 "channel" with a value greater than 6 will be set to True.  All
# other samples will be False.
jellies = Sv_18m38 > 6

# Now we're going to do an in-place and-ing where we'll take the results of
# our first operation and AND them with the results of this operation where
# we're setting all samples in the Sv_120m38 channel to True if they are less
# than -1.
jellies &= Sv_120m38 < -1

# Here we'll get crazy and do two comparisons and AND the results.  Masks
# support boolean operations, both in-place and regular.  Just make sure you
# group your expressions since the boolean operators have a higher precedence
# than the comparison operators.
euphausiids = (Sv_120m38 > 9) & (Sv_18m38 < -5)

# Do another comparison.
myctophids = (Sv_18m38 < -9) & (Sv_120m38 < -8)

# Do another comparison.  This one will comprise the results of 4 comparisons.
fish = (Sv_18m38 < 2) & (Sv_18m38 > -4)
fish &= (Sv_120m38 < 0) & (Sv_120m38 > -6)


# Now lets create a ProcessedData object the same shape as our other data
# arrays but with the data array set to zeros.
diff_results = Sv_18m38.zeros_like()

# Also, we'll use the masks to set the various samples to values that
# represent what we think they are.
diff_results[jellies] = 4
diff_results[euphausiids] = 7
diff_results[myctophids] = 15
diff_results[fish] = 18


# Create a matplotlib figure to plot our echograms on.
fig = figure()
subplots_adjust(left=0.1, bottom=.1, right=0.98, top=.90, wspace=None,
                hspace=1.5)

# Plot the original data.
ax = fig.add_subplot(4, 1, 1)
# Use the view method to return a processed_data object that is a view into
# our original data. We will plot all pings and samples 0-2000.
v_data = Sv_data[freq_chan_map[18000]].view((None, None, None),(0, 2000, None))
eg = echogram.Echogram(ax, v_data, threshold=[-70, -34])
ax.set_title("Original 18 kHz Sv Data")

ax = fig.add_subplot(4, 1, 2)
v_data = Sv_data[freq_chan_map[38000]].view((None, None, None), (0, 2000, None))
eg = echogram.Echogram(ax, v_data, threshold=[-70,-34])
ax.set_title("Original 38 kHz Sv Data")

ax = fig.add_subplot(4, 1, 3)
v_data = Sv_data[freq_chan_map[120000]].view((None, None, None), (0, 2000, None))
eg = echogram.Echogram(ax, v_data, threshold=[-70, -34])
ax.set_title("Original 120 kHz Sv Data")


# Plot our differencing data.
ax = fig.add_subplot(4, 1, 4)
v_results = diff_results.view((None, None, None), (0, 2000, None))
# note that we set the threshold to something that will work with the values
# we assigned to our results.
eg = echogram.Echogram(ax, v_results, threshold=[0, 20])
ax.set_title('Differencing results')

# Display the results.
show()
