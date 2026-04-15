

from glob import glob
import matplotlib.pyplot as plt
from echolab2.instruments import echosounder
from echolab2.processing import integration, grid
from echolab2.plotting.matplotlib import echogram


#  create a callback that prints a progress bar when reading
def print_progress(filename, cumulative_pct, cumulative_bytes,
                                callback_ref):
    if (cumulative_pct % 5):
        print('.', end='')
        
    if cumulative_pct >= 100:
        print()


#  set the path to the raw file 
raw_file =  'C:/EK Test Data/EK80/DY2600/cwfm/DY2600-D20260121-T213843.raw'

#  get a list of all available calibration files. The calibrations will be
#  matched to the data when they are read. Calibrations for channels not in
#  the data will be ignored.
cal_files = (glob('C:/EK Test Data/EK80/DY2600/cwfm//calibration/FM/*.xml') +
        glob('C:/EK Test Data/EK80/DY2600/cwfm//calibration/CW/*.xml'))

#  specify the channels to read. Set this to None or [] to read all channels
#  in the data file. In the example below, we read in 4 38 kHz channels.
#  1 is CW active, 2 is CW passive, 3 is FM active and 4 is FM passive. The
#  matching passive channels will be used for noise correction.
channel_ids=['WBT 978217-15 ES38-7_1','WBT 978217-15 ES38-7_2',
             'WBT 978217-15 ES38-7_3','WBT 978217-15 ES38-7_4']

#  read the raw file(s)
print('Reading' + raw_file, end='')
if channel_ids:
    ek_data = echosounder.read([raw_file], channel_ids=channel_ids,
            progress_callback=print_progress)
else:
    ek_data = echosounder.read([raw_file], progress_callback=print_progress)

#  now read in the calibration data. It is strongly encouraged that you
#  provide calibration files especially with FM data.
print("Reading calibration data...")
cal_data = echosounder.get_calibration_from_xml(ek_data, cal_files)

#  the get_calibration method will match the calibrations to the channels
#  in the data. We'll just print out what we found.
print('Read cal data for channels: ' + str(cal_data.keys()))

print("Getting Sv...")
# calculate Sv using the calibration object, and the default frequency domain method for FM
Sv = echosounder.get_Sv(ek_data, calibration=cal_data, frequency_resolution=1000)


Sv_noise_correct = echosounder.noise_correct(Sv, SNR_threshold=10)

#  apply 
Sv_bottom_clean = echosounder.apply_boundary_exclusions(Sv_noise_correct, exclude_above_line=5, bottom_offset=0.5)

integrated_data = {}

for channel in Sv_bottom_clean.keys():    # Perform integration, with or without bottom line exclusion
    # Initialize integrator without minimum threshold
    integrator = integration.integrator(min_threshold_applied=False)
    
    # Initialize grid for integration
    int_grid = grid.grid(interval_length=50, interval_axis='ping_number',
            layer_axis='range', layer_thickness=5,
            data=Sv_bottom_clean[channel])
    
    integrated_data[channel] = integrator.integrate(Sv_bottom_clean[channel], int_grid) 
    
    if Sv_bottom_clean[channel].frequency.size > 1:
        for f_idx, f in enumerate(Sv_bottom_clean[channel].frequency):
            #  create a figure to plot our echogram in
            fig = plt.figure()
            eg = echogram.Echogram(fig, Sv_bottom_clean[channel], threshold=[-70,-34],
                    frequency=f, cmap='viridis', grid=int_grid)
            eg.plot_integration_results(integrated_data[channel])
            if hasattr(Sv[channel], 'bottom_line'):
                eg.plot_line(Sv[channel].bottom_line, color=[0.9,0,0], linewidth=1.0)
            eg.add_colorbar(fig)
            eg.axes.set_title("Sv " + channel + " " + str(f) + " Hz")
    else:
        #  create a figure to plot our echogram in
        f = Sv_bottom_clean[channel].frequency
        fig = plt.figure()
        eg = echogram.Echogram(fig, Sv_bottom_clean[channel], threshold=[-70,-34],
                cmap='viridis', grid=int_grid)
        eg.plot_integration_results(integrated_data[channel])
        if hasattr(Sv[channel], 'bottom_line'):
            eg.plot_line(Sv[channel].bottom_line, color=[0.9,0,0], linewidth=1.0)
        eg.add_colorbar(fig)
        eg.axes.set_title("Sv " + channel + " " + str(f) + " Hz")


# Show our figures
plt.show()

print('done.')
