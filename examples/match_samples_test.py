# -*- coding: utf-8 -*-
"""match_samples_test.py tests the match_samples method by way of the insert
method. The append/insert methods will match samples of the "other" object
automagically if needed. This script tests this and then plots the results.

Two processed data objects are created with synthetic data. One has samples
with values alternating between 5 and 0, and the other has samples half the
thickness of the first, with alternating values of 10 and 0.

When the higher resolution data is inserted into the lower resolution data,
pairs of (10,0) samples from the higher resolution data will be averaged, 
resulting in a sample with the same thickness as the lower resultion samples, 
with a value of 5 (matching the lower resultion samples.)

When the lower resolution data is inserted into the higher resolution data,
the thicker samples are split into thinner samples. While you would assume the
result would be 2 samples with a value of 5 for each original sample with value
5, the target resampling grid doesn't quite line up that way. Since the ranges
both start at esentially 0, but the target grid's samples are 1/2 the thickness,
the target grid is shifted 1/4 of a sample relative to the source data. This
results in a pattern where sample values cycle from 5, 2.5, 0, 2.5 and back to
5.
"""
import numpy as np
from matplotlib.pyplot import figure, show
from echolab2.processing import processed_data, line
from echolab2.plotting.matplotlib import echogram

# create some synthetic data

# The first data object will be a 10 ping x 1000 sample Sv object. It will be populated
# with a 1 sample x 1 sample checkerboard of values 0 and 5. The sample thickness will
# be 0.1 meters, creating a range vector of 100 values from 1e-20 to 9.99 m.
#
# Note that we will set jitter_time to False and accept the default start_time so
# both data objects will have the same ping_times.
test_array_size = (10,100)
data_block_shape = (1,1)
data_block_values = (5,0)
sample_thickness_m = 0.1

fake_Sv_1 = processed_data.create_test_data(test_array_size[0], test_array_size[1],
        block_shape=data_block_shape, block_values=data_block_values, data_type='Sv', 
        jitter_time=False, sample_thickness=sample_thickness_m, range_start=0)


# The second data object will also be a 10 ping x 100 sample Sv object populated
# with a 1 sample x 1 sample checkerboard but the sample thickness will be 0.05m which
# is half of the first object's 0.1 m thickness. The checkerboard sample values will
# be 0 and 10, double the of values the first object. The range vector will be 1000
# values from 1e-20 to 4.995 m.
test_array_size = (10,100)
data_block_shape = (1,1)
data_block_values = (10,0)
sample_thickness_m = 0.05

fake_Sv_2 = processed_data.create_test_data(test_array_size[0], test_array_size[1],
        block_shape=data_block_shape, block_values=data_block_values, data_type='Sv', 
        jitter_time=False, sample_thickness=sample_thickness_m, range_start=0)


# Since the insert method alters the existing object, create a copy of fake_Sv_1
# to insert data into so we can plot both the original fake_Sv_1 and the copy with
# data from fake_Sv_2 inserted into it.
fake_Sv_3 = fake_Sv_1.copy()

# Now insert data from fake_Sv_2 into fake_Sv_3. Since these objects have different
# sample thicknesses, data from fake_Sv_2 will be resampled prior to insertion. We
# set time_order=True to insert the data in time order. Since both objects share
# the same ping times, the data will be interleaved.
#
# The expected result is an array 20 pings by 100 samples with a range from 1e-20
# to 9.99 m. The 0.05 m thick samples of fake_Sv_2 will be averaged into 0.1 m samples.
# Since the fake_Sv_2 sample values alternate between 10 and 0, the value of the
# averaged sample will be 5 which matches the data in fake_Sv_1.
fake_Sv_3.insert(fake_Sv_2, time_order=True)


# Plot the original data
fig = figure()
eg = echogram.Echogram(fig, fake_Sv_1, threshold=[0, 15], show_grid=True)
eg.axes.set_title("Fake Sv 1 with 0.1m sample thickness")

fig = figure()
eg = echogram.Echogram(fig, fake_Sv_2, threshold=[0, 15], show_grid=False)
eg.axes.set_title("Fake Sv 2 with 0.05m sample thickness")

# Draw the horizontal lines representing the target vertical resampling grid.
# This helps visualize how the thinner samples of fake_Sv_2 will be combined into
# the thicker samples of fake_Sv_1
h_line = line.line(ping_time=np.array([fake_Sv_3.ping_time[0], fake_Sv_3.ping_time[-1]]),
        data=np.array([0, 0]), color=[0.5,0,0.5])
target_grid_edges = np.append(fake_Sv_3.range - fake_Sv_3.sample_thickness / 2, 
        fake_Sv_3.range[-1] + fake_Sv_3.sample_thickness / 2)
for grid_line in target_grid_edges:
    h_line.data = np.array([grid_line,grid_line])
    eg.plot_line(h_line)


# And plot the merged data.
fig = figure()
eg = echogram.Echogram(fig, fake_Sv_3, threshold=[0, 15], show_grid=True)
eg.axes.set_title("Fake Sv 2 inserted into Fake Sv 1")

# Display figure.
show()


# Now do this the other way. Make a copy of fake_Sv_2
fake_Sv_3 = fake_Sv_2.copy()

# Now insert data from fake_Sv_1 into fake_Sv_3. The expected result is an array of
# 20 pings by 100 samples with ea range from 1e-20 to 4.995 m. Samples from
# fake_Sv_1 with ranges > 49.95 are dropped. The fake_Sv_1 samples will be split
fake_Sv_3.insert(fake_Sv_1, time_order=True)


# Plot the original data
fig = figure()
eg = echogram.Echogram(fig, fake_Sv_1, threshold=[0, 15], show_grid=False)
eg.axes.set_title("Fake Sv 1 with 0.1m sample thickness")

# Now draw the horizontal lines representing the target vertical resampling grid
# This helps visualize how the thicker samples of fake_Sv_1 will be split up into
# the thinner samples of fake_Sv_2.
h_line = line.line(ping_time=np.array([fake_Sv_3.ping_time[0], fake_Sv_3.ping_time[-1]]),
        data=np.array([0, 0]), color=[0.5,0,0.5])
target_grid_edges = np.append(fake_Sv_3.range - fake_Sv_3.sample_thickness / 2, 
        fake_Sv_3.range[-1] + fake_Sv_3.sample_thickness / 2)
for grid_line in target_grid_edges:
    h_line.data = np.array([grid_line,grid_line])
    eg.plot_line(h_line)


fig = figure()
eg = echogram.Echogram(fig, fake_Sv_2, threshold=[0, 15], show_grid=True)
eg.axes.set_title("Fake Sv 2 with 0.05m sample thickness")


# And plot the merged data.
fig = figure()
eg = echogram.Echogram(fig, fake_Sv_3, threshold=[0, 15], show_grid=True)
eg.axes.set_title("Fake Sv 1 inserted into Fake Sv 2")

# Display figure.
show()


print()
