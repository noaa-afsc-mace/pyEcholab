# -*- coding: utf-8 -*-
"""cw_fm_splitter.py

This is a simple example showing how to split cw channels from a combined 
cw/fm file. into a cw only file. The concepts can easily be extended to write a
fm only file, or drop specific channels from a file.

"""

import os
from echolab2.instruments import echosounder



#  create a callback that prints a progress bar when reading and writing.
#  This is not required and is purely a convienience shown here  as an
#  example of how one would use this in a console application.
def print_progress(filename, cumulative_pct, cumulative_bytes, callback_ref):
    
    if (cumulative_pct % 5):
        print('.', end='')
        
    if cumulative_pct >= 100:
        print()



# set the input and output filenames
cwfm_in_file = 'C:/EK Test Data/cwfm/DY2600-D20260122-T004923.raw'
cw_out_file =  'C:/EK Test Data/cwfm/DY2600-D20260122-T004923-cw.raw'

#  parse the input filename and path and set output filename. The write_raw
#  method expects a dictionary with infile:outfile key value pairs.
infile_parts = os.path.split(cwfm_in_file)
ek80_out_file = {infile_parts[1]:cw_out_file}

#  read the input raw file
print('reading ' + cwfm_in_file)
ek_data = echosounder.read([cwfm_in_file], progress_callback=print_progress)

#  Iterate thru all channels and create a list of the cw channels to write.
cw_channels = []
for channel in ek_data.channel_ids:
    #  check if this channel is cw. The raw_data class also has the is_fm,
    #  and is_passive methods which could be used to include/exclude channels
    #  based on those properties too.
    if ek_data.get_channel_data()[channel][0].is_cw():
        print('found cw channel: ' + channel)
        cw_channels.append(channel)
n_channels = len(cw_channels)

#  in this example we are writing cw only files so if we don't have
#  any cw channels, we have nothing to do and skip the writing.
if n_channels > 0:

    #  we have at least one cw channel so write the output file
    print('Writing ' + cw_out_file)
    ek_data.write_raw(ek80_out_file, channel_ids=cw_channels, overwrite=True,
            progress_callback=print_progress)
            
print('done.')

