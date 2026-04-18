#!/usr/bin/env python
'''get_rawfile_info.py reads the raw file passed as a command line argument
and prints out a bunch of information about the raw file.

While useful on its own, this example shows how to access the file and
channel metadata within a raw file.
'''

import os
import argparse
import numpy as np
from echolab2.instruments import echosounder, EK80


#  create the argument parser. Set the application description.
parser = argparse.ArgumentParser(description='get_rawfile_info')

#  specify the positional arguments: filename
parser.add_argument("raw_file", help="The full path to the raw file to extract info from.")

#  parse our arguments
args = parser.parse_args()

#  normalize the path
raw_file = os.path.normpath(args.raw_file)

#  get the file size
size_mib = os.path.getsize(raw_file)  / (1024 * 1024)

#  get the base file name
filename = os.path.basename(raw_file)

#  read the raw file
ek_data = echosounder.read(raw_file)

#  compute the time span and data rate
time_span_minutes = (ek_data.end_time-ek_data.start_time) / np.timedelta64(1,'m')
rate_mib_per_hour = size_mib / (time_span_minutes / 60.)

#  get some info about the motion data
try:
    n_motion = len(ek_data.motion_data.times)
    has_pitch = np.any(ek_data.motion_data.pitch != 0.0)
    has_roll = np.any(ek_data.motion_data.roll != 0.0)
    has_heading = np.any(ek_data.motion_data.heading != 0.0)
    has_heave = np.any(ek_data.motion_data.heave != 0.0)
    has_extended = ek_data.motion_data.has_MRU1
    motion_start = ek_data.motion_data.times[0]
    motion_end = ek_data.motion_data.times[-1]
except:
    n_motion = 0
    has_pitch = False
    has_roll = False
    has_heading = False
    has_heave = False
    has_extended = False
    motion_start = None
    motion_end = None

#  get some info about the nmea data
try:
    n_nmea = ek_data.nmea_data.n_raw
    nmea_start = ek_data.nmea_data.nmea_times[0]
    nmea_end = ek_data.nmea_data.nmea_times[-1]
    nmea_talkers = ','.join(ek_data.nmea_data.talker_ids)
    nmea_message_ids= ','.join(ek_data.nmea_data.message_ids)
except:
    n_nmea = 0
    nmea_start = None
    nmea_end = None
    nmea_talkers = ''
    nmea_talkers = ''

#  get a reference to the file configuration header. When a file is read,
#  the header is parsed by channel ID and the global and channel specific
#  parameters are stored and a reference to these data are assigned to each
#  ping in the raw_data.configuration attribute. Since we are only reading
#  a single file, all of these references will point to the same data and
#  we can just grab the first one.
configuration = ek_data.raw_data[ek_data.channel_ids[0]][0].configuration[0]

# The EK60 and EK80 objects represent two general eras of raw file formats
# that containin different metadata so we have to handle them differently
if isinstance(ek_data, EK80.EK80):

    #  print out data about this EK80 like data

    print()
    print(f"get_rawfile_info.py")
    print()
    print(f"                  Raw data file: {filename}")
    print(f"   Acquisition application name: {configuration['application_name']}")
    print(f"Acquisition application version: {configuration['application_version']}")
    print(f"        Raw file format version: {configuration['file_format_version']}")
    print(f"                      Ping Mode: {configuration['ping_mode']}")
    print()
    print(f"                Number of pings: {ek_data.n_pings}")
    print(f"                First ping time: {ek_data.start_time}")
    print(f"                 Last ping time: {ek_data.end_time}")
    print(f"                 Data time span: {time_span_minutes:.2f} minutes")
    print(f"                      File size: {size_mib:.2f} MiB")
    print(f"           Data collection rate: {rate_mib_per_hour:.3f} MiB per hour")
    print()
    print(f"     Number of Motion datagrams: {n_motion}")
    print(f"     Motion datagram start time: {motion_start}")
    print(f"       Motion datagram end time: {motion_end}")
    print(f"                    Motion data: pitch:{has_pitch} roll:{has_roll} heading:{has_heading} heave:{has_heave}")
    print(f"       Has extended motion data: {has_extended}")
    print()
    print(f"       Number of NMEA datagrams: {n_nmea}")
    print(f"       NMEA datagram start time: {nmea_start}")
    print(f"         NMEA datagram end time: {nmea_end}")
    print(f"                NMEA talker IDs: {nmea_talkers}")
    print(f"               NMEA message IDs: {nmea_message_ids}")
    print()
    print(f"             Number of channels: {ek_data.n_channels}")

    for chan in ek_data.channel_ids:
        #  since we're only reading a single file, we can assume there will only
        #  be one data type associated with this channel
        chan_raw = ek_data.raw_data[chan][0]

        #  get a reference to this channel's configuration - again since we're
        #  reading a single file, all of these references will point to the same
        #  data so we just grab the first one.
        configuration = ek_data.raw_data[chan][0].configuration[0]

        #  compute the ping intervals
        ping_intervals_ms = np.diff(chan_raw.ping_time)

        #  get the FPGA firmware versions
        xcvr_version_info = configuration['transceiver_version'].split('\r\n')
        try:
            fpga_tx_firmware_ver = xcvr_version_info[8].split(":")[1].strip()
        except:
            fpga_tx_firmware_ver = 'unable to parse'
        try:
            fpga_rx_firmware_ver = xcvr_version_info[9].split(":")[1].strip()
        except:
            fpga_rx_firmware_ver = 'unable to parse'

        #  is it FM?
        if chan_raw.pulse_form[0] > 0:
            pulse_form = 'FM'
        else:
            pulse_form = 'CW'

        #  active or passive?
        if np.all(chan_raw.channel_mode):
            channel_mode = 'Passive'
        elif np.any(chan_raw.channel_mode):
            #  I don't know if these files can be created, but JIC
            channel_mode = 'Mixed active/passive'
        else:
            channel_mode = 'Active'

        print("")
        print(f"                         Channel ID: {chan_raw.channel_id}")
        print(f"                 Transceiver number: {configuration['transceiver_number']}")
        print(f"                   Transceiver Type: {configuration['transceiver_type']}")
        print(f"          Embedded software version: {configuration['transceiver_software_version']}")
        print(f"           FPGA TX firmware version: {fpga_tx_firmware_ver}")
        print(f"           FPGA RX firmware version: {fpga_rx_firmware_ver}")
        print(f"             Transceiver IP address: {configuration['ip_address']}")
        print(f"     Hardware channel configuration: {configuration['hw_channel_configuration']}")
        print(f"                    Transducer name: {configuration['transducer_name']}")
        print(f"           Transducer serial number: {configuration['transducer_serial_number']}")
        print(f"               Transducer beam type: {configuration['transducer_beam_type']}")
        print(f"               Transducer frequency: {configuration['transducer_frequency']}")
        print(f"       Transducer frequency minimum: {configuration['transducer_frequency_minimum']}")
        print(f"       Transducer frequency maximum: {configuration['transducer_frequency_maximum']}")
        print()
        print(f"                         Pulse form: {pulse_form}")
        print(f"                       Channel mode: {channel_mode}")
        print(f"                   Sample data type: {chan_raw.data_type}")
        print(f"                    Number of pings: {chan_raw.n_pings}")
        print(f"                  Number of samples: {chan_raw.n_samples}")
        print(f"                    First ping time: {chan_raw.ping_time[0]}")
        print(f"                     Last ping time: {chan_raw.ping_time[-1]}")
        print(f"              Minimum ping interval: {np.min(ping_intervals_ms)}")
        print(f"              Maximum ping interval: {np.max(ping_intervals_ms)}")
        print(f"                 Mean ping interval: {np.mean(ping_intervals_ms)}")

else:

    #  print out info about this EK60 like data

    print()
    print(f"get_rawfile_info.py")
    print()
    print(f"                  Raw data file: {filename}")
    print(f"   Acquisition application name: {configuration['sounder_name']}")
    print(f"Acquisition application version: {configuration['version']}")
    print(f"                    Survey name: {configuration['survey_name']}")
    print(f"                  Transect name: {configuration['transect_name']}")
    print()
    print(f"                Number of pings: {ek_data.n_pings}")
    print(f"                First ping time: {ek_data.start_time}")
    print(f"                 Last ping time: {ek_data.end_time}")
    print(f"                 Data time span: {time_span_minutes:.2f} minutes")
    print(f"                      File size: {size_mib:.2f} MiB")
    print(f"           Data collection rate: {rate_mib_per_hour:.3f} MiB per hour")
    print()
    print(f"     Number of Motion datagrams: {n_motion}")
    print(f"     Motion datagram start time: {motion_start}")
    print(f"       Motion datagram end time: {motion_end}")
    print(f"                    Motion data: pitch:{has_pitch} roll:{has_roll} heading:{has_heading} heave:{has_heave}")
    print(f"       Has extended motion data: {has_extended}")
    print()
    print(f"       Number of NMEA datagrams: {n_nmea}")
    print(f"       NMEA datagram start time: {nmea_start}")
    print(f"         NMEA datagram end time: {nmea_end}")
    print(f"                NMEA talker IDs: {nmea_talkers}")
    print(f"               NMEA message IDs: {nmea_message_ids}")
    print()
    print(f"             Number of channels: {ek_data.n_channels}")

    for chan in ek_data.channel_ids:
        #  since we're only reading a single file, we can assume there will only
        #  be one data type associated with this channel
        chan_raw = ek_data.raw_data[chan][0]

        #  get a reference to this channel's configuration - again since we're
        #  reading a single file, all of these references will point to the same
        #  data so we just grab the first one.
        configuration = ek_data.raw_data[chan][0].configuration[0]

        #  compute the ping intervals
        ping_intervals_ms = np.diff(chan_raw.ping_time)

        pulse_form = 'CW'

        #  active or passive?
        if np.all(chan_raw.transmit_mode):
            channel_mode = 'Passive'
        elif np.any(chan_raw.transmit_mode):
            #  I don't know if these files can be created, but JIC
            channel_mode = 'Mixed active/passive'
        else:
            channel_mode = 'Active'

        print("")
        print(f"                         Channel ID: {chan_raw.channel_id}")
        print(f"                 Transceiver number: {configuration['channel_number']}")
        print(f"                   Transceiver Type: {chan_raw.transceiver_type}")
        print(f"          Embedded software version: {configuration['gpt_software_version']}")

        print(f"               Transducer beam type: {configuration['beam_type']}")
        print(f"               Transducer frequency: {configuration['frequency']}")
        print()
        print(f"                         Pulse form: {pulse_form}")
        print(f"                       Channel mode: {channel_mode}")
        print(f"                   Sample data type: {chan_raw.data_type}")
        print(f"                    Number of pings: {chan_raw.n_pings}")
        print(f"                  Number of samples: {chan_raw.n_samples}")
        print(f"                    First ping time: {chan_raw.ping_time[0]}")
        print(f"                     Last ping time: {chan_raw.ping_time[-1]}")
        print(f"              Minimum ping interval: {np.min(ping_intervals_ms)}")
        print(f"              Maximum ping interval: {np.max(ping_intervals_ms)}")
        print(f"                 Mean ping interval: {np.mean(ping_intervals_ms)}")

print()
