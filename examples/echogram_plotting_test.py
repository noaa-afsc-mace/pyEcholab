# -*- coding: utf-8 -*-
"""echogram_plotting_test.py demonstrates plotting echograms using fake data.

It is primarily intended to test the processed_data and echogram classes. As
configured, this will plot a 60 pings x 100 samples array with samples that
are 0.5 m thick. Ping times will be jittered to simulate real data. A grid
will then be plotted on top of the samples.

If everything is working, the samples should be centered perfectly on the
grid lines.

"""

import numpy as np
from matplotlib.pyplot import figure, show
from echolab2.processing import processed_data, line
from echolab2.plotting.matplotlib import echogram


# Set adjust_x to False to see the effect of imshow's regular grid on
# plotting. When this is False, samples in the middle of the echogram will
# not be aligned properly on the grid.
adjust_x = True


# Create a processed_data object with test data. In this case, a 60 pings
# by 101 samples 1x1 checkerboard array with values alternating between 10
# and 0 and a sample thickness of 0.5 m. This gives us a sample array with a
# range of 0-50m
#
# We will set jitter_time to True to enable a non-uniform x axis values which
# mimicks real echosunder data. Since matplotlib's imshow displays the sample
# data on a regular grid, using non-uniform ping intervals tests the echogram's
# method for adjusting for the alongship shifts this creates.
test_array_size = (60,101)
data_block_shape = (1,1)
data_block_values = (10,0)
sample_thickness_m = 0.5

fake_Sv = processed_data.create_test_data(test_array_size[0], test_array_size[1],
        block_shape=data_block_shape, block_values=data_block_values, data_type='Sv', 
        jitter_time=True, sample_thickness=sample_thickness_m, range_start=0)


# Create a matplotlib figure to plot our echograms on and then create the echogram.
# Set show_grid to False since we'll plot our own simple grid that we center on
# our samples.
fig_1 = figure()
eg = echogram.Echogram(fig_1, fake_Sv, threshold=[0, 15], show_grid=False,
                adjust_x=adjust_x)
eg.axes.set_title("Echogram Plot Test")

# Create a horizontal line that we'll draw every 10 m at the vertical center of each sample.
h_line = line.line(ping_time=np.array([fake_Sv.ping_time[0], fake_Sv.ping_time[-1]]),
        data=np.array([0, 0]), color=[0.3,0,0.3])
n_lines = int((fake_Sv.range[-1] / (sample_thickness_m * 10)) + 0.5) + 1
print(n_lines)

# Now draw the horizontal lines
for i in range(n_lines):
    eg.plot_line(h_line)
    h_line = h_line + (sample_thickness_m * 10)
    print(h_line.data)

# Create a vertical line that we'll draw at the alongship center of each sample.
v_line = line.line(ping_time=np.array([fake_Sv.ping_time[0], fake_Sv.ping_time[0]]),
        data=np.array([fake_Sv.range[0], fake_Sv.range[-1]]), color=[0.3,0,0.3])

# Now draw the vertical lines
for t in fake_Sv.ping_time:
    v_line.ping_time = np.array([t,t])
    eg.plot_line(v_line)
    
# Display figure.
show()
