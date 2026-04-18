# -*- coding: utf-8 -*-
"""Reads raw and bot/out files, detects bottom using super simple detector
and plots data.

This example implements a very simple bottom detector. It has been written to
demonstrate and explore how and where bottom detection would fit into the
current pyEcholab2 framework.

The afsc_bot_detector class used in this example is not intended to be used
for any real science. This is simply something to get development started and
to work out details around implementing additional processing algorithms.
"""

from matplotlib.pyplot import figure, show
from echolab2.instruments import echosounder
from echolab2.plotting.matplotlib import echogram
from echolab2.processing import afsc_bot_detector


# Create a list of .raw files to read
rawfiles = ['./data/EK60/DY1603/raw/DY1603_EK60-D20160308-T115724.raw']

# read the data - this will return a dictionary keyed by channel ID,
# containing the raw data objects with data for that channel. Bottom data
# files will also be read if available.
print("Reading .raw data...")
raw_data = echosounder.read(rawfiles, frequencies=[38000, 120000])
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

# Sometimes it's convienient to reference your data by frequency but the 
# echosounder.get_* methods return data keyed by channel ID. The 
# echosounder.get_data_by_frequency() method will return a reference given
# a frequency.
Sv_38 = echosounder.get_data_by_frequency(Sv_data, 38000.)
Sv_120 = echosounder.get_data_by_frequency(Sv_data, 120000.)

#  create an instance of our bottom detector and set the minimum detection distance.
#  since we'll be passing our bottom detector data on a depth grid, set this to the
#  minimum DEPTH in meters to search for the bottom. We'll also set the backstep.
#  Since the data are Sv, the backstep will be in dB.
bot_detector = afsc_bot_detector.afsc_bot_detector(search_min=15, backstep=32)


#  now use our simple detector to pick a bottom line for the 38. First we set the
#  backstep to a value appropriate for 38 kHz. The bottom detector class returns a
# pyEcholab2 Line object representing the bottom.

#  detect the bottom on the 38 kHz data
bottom_38_detected = bot_detector.detect(Sv_38)

#  Set this line's color to pink
bottom_38_detected.color = [1,0.4,0.8]

#  and pick the 120 bottom and make it pink. We'll first shorten the smoothing
#  window to something that works a bit better at higher frequencies.
bot_detector.window_len = 5
bottom_120_detected = bot_detector.detect(Sv_120)
bottom_120_detected.color = [1,0.4,0.8]

# Create matplotlib figures and display the results.
fig_38 = figure()
eg = echogram.Echogram(fig_38, Sv_38, threshold=[-70, -34])
fig_38.suptitle("Heave Corrected - 38kHz", fontsize=14)
eg.axes.set_title("With sounder detected (purple) and software detected " +
        "(pink) bottom lines", fontsize=10)
eg.add_colorbar(fig_38)
eg.plot_line(Sv_38.bottom_line, linewidth=2)
eg.plot_line(bottom_38_detected, linewidth=2)

fig_120 = figure()
eg = echogram.Echogram(fig_120, Sv_120, threshold=[-70, -34])
fig_120.suptitle("Heave Corrected - 120kHz", fontsize=14)
eg.axes.set_title("With sounder detected (purple) and software detected " +
        "(pink) bottom lines", fontsize=10)
eg.add_colorbar(fig_120)
eg.plot_line(Sv_120.bottom_line, linewidth=2)
eg.plot_line(bottom_120_detected, linewidth=2)

show()
