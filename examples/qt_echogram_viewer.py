# -*- coding: utf-8 -*-
"""This an example showing how to use the PyQt based echogram_viewer
which is a high level interface to the AFSC QImageViewer library.
While it can be used on its own, echogram_viewer is primarily an example
of how to embed the QEchogramViewer in your own GUI application

This example reads a couple of ES60 .raw files, reads the .out file,
gets Sv, and plots the echogram and the bottom line. It also shows
how

The QImageViewer library Requires PyQt6.

QEchogramViewer implements panning and zooming. When the application
has focus and the mouse is within the echogram window, you can pan and
zoom using the keyboard and mouse. The controls are:

    <CTRL> + Mouse Wheel zooms in/out
    <ALT> + Click and drag will pan

"""

import sys
from PyQt6 import QtWidgets
from echolab2.instruments import EK60
from echolab2.plotting.qt import echogram_viewer
from echolab2.processing import afsc_bot_detector


def read_write_callback(filename, cumulative_pct, cumulative_bytes, userref):
    '''
    read_write_callback is a simple example of using the progress_callback
    functionality of the EK60.read_raw and EK60.write_raw methods.
    '''
    if cumulative_pct > 100:
        return
    if cumulative_pct == 0:
        sys.stdout.write(filename)
    if cumulative_pct % 4:
        sys.stdout.write('.')
    if cumulative_pct == 100:
        sys.stdout.write('  done!\n')


# Specify a couple of raw file.
rawfiles = ['./data/ES60/OEX2017/raw/L0189-D20170624-T122804-ES60.raw',
            './data/ES60/OEX2017/raw/L0189-D20170624-T134748-ES60.raw']

# Specify the associated out file.
outfile = './data/ES60/OEX2017/raw/L0189-D20170624-T122804-ES60.out'

# For this example, we are going to do things The Hard Way(tm). Since reading
# ES60 .out files is currently not supported by the echosounder.read() function,
# we can't enjoy its convieniences and will read this data using the underlying
# classes directly.

# Create an instance of the EK60 instrument.
ek60 = EK60.EK60()

# Read the data. The progress_callback is not required and only shown here
# to be fancy. It is useful in GUI applications when you want to present
# the user with a progress bar.
print('Reading raw file:')
ek60.read_raw(rawfiles, progress_callback=read_write_callback,
        frequencies=38000)

# Next, read in the bottom data. For ES60 systems, bottom detections
# were recorded in .out files so we call the EK60.read_out() method.
# This should be done AFTER you read the raw file as .out files contain
# data for multiple .raw files and detections are matched to pings that
# are currently stored in your EK60 object.
ek60.read_out(outfile, progress_callback=read_write_callback)

# Get the 38 kHz raw data.
print('Getting Sv...')
raw_data = ek60.get_channel_data(frequencies=38000)

# Remember that when doing it The Hard Way (tm), raw data returned from
# get_channel_data is, for reasons, always in a list. Since I know that
# there is only one data type contained within the files we read, I know
# there will only be a single item in my list so I just grab the first one.
raw_data_38 = raw_data[38000][0]

# Print out some info about the data we just read in
print(raw_data_38)

# Get a calibration object - This returns a cal object that is populated
# with data from the .raw file.
cal_obj = raw_data_38.get_calibration()

# You can change cal parameters as needed. For example, say that
# at the time of recording an incorrect transducer draft and sound
# speed were entered into the system. We can change that here:
cal_obj.transducer_depth = 7.5
cal_obj.sound_speed = 1500

# Get Sv data. Since we will be plotting this with out bottom detections
# we will set return_depth to True to get data on a depth grid.
Sv_38 = raw_data_38.get_Sv(calibration=cal_obj, return_depth=True)

# Get the bottom detections as an echolab2 line. Bottom detections are
# *always* stored as depth with heave compensation applied (if applicable.)
# We must pass the calibration object to the raw_data.get_bottom method to 
# ensure that the depths are corrected for any changes between the recorded
# sound speed and transducer draft and the commanded sound speed and
# transducer draft values in the calibration object.
bottom_line = raw_data_38.get_bottom(calibration=cal_obj)

# The default line thickness of 1 is quite thin for display using
# QEchogram viewer. Change the thickness so the bottom line is a bit
# more visible
bottom_line.linewidth = 3.0


print('Plotting...')

# Create an application instance.
app = QtWidgets.QApplication([])

# Create the main application window and show it
eg_viewer = echogram_viewer.echogram_viewer()

#  show the application window
eg_viewer.show()

# Set the echogram data
eg_viewer.update_echogram(Sv_38)

# Add our bottom line - You will note that when you zoom in and look
# at the bottom line, it is correctly placed in the water column, even
# though we used a different sound speed and transducer draft than
# what was used when the data were recorded. This is the desired and
# expected result since we passed our calibration object with the
# modified values to the raw_data.get_bottom() method.
eg_viewer.add_line(bottom_line)

#  save the echogram at full resolution. Curently this does not
#  include the horizontal scaling that is appled to the sample data
#  to reduce the vertical exaggeration of the echograms.
#eg_viewer.save_image('test.png')

# Start event processing.
app.exec()
