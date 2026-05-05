#!/usr/bin/env python
'''echosounder_get_rawfile_info.py is a simple example showing how the
echosounder.get_rawfile_info() function works. 

echosounder.get_rawfile_info attempts to quickly extract the bulk of the
available raw file metadata. It does this by reading the configuration header
and, if available, the .idx file, to determine ping span, time span, channels,
and a host of transceiver and transducer parameters organized by channel ID.
If an .idx file is not available, it will resort to reading the whole file
which will be slower.

The amount of metadata between EK60 and EK80 files varies a lot. The EK80
file formats contain a lot more and the format has evolved to add more so
many fields may be empty depending on the raw file verion you read.

Since it (usually) doesn't read the entire raw file, it does not return 
information about NMEA or motion datagrams. If you need this, See the example
get_rawfile_info.py for an, er, example of extracting this information.

'''

from echolab2.instruments import echosounder




raw_file = './data/EK80/cwfm/DY2602/raw/DY2602-D20260310-T191652.raw'



#  read the raw file
file_info = echosounder.get_rawfile_info(raw_file)

size_mib = file_info['raw_file_size_bytes'] / 1024. / 1024.
n_channels = len(file_info['channels'])

print()
print(f"echosounder_get_rawfile_info.py")
print()
print(f"                  Raw data file: {file_info['raw_file_name']}")
print(f"                 .idx data file: {file_info['idx_file_name']}")
print(f"   Acquisition application name: {file_info['application_name']}")
print(f"Acquisition application version: {file_info['application_version']}")
print(f"        Raw file format version: {file_info['file_format_version']}")
print(f"                      Ping Mode: {file_info['ping_mode']}")
print()
print(f"    Starting global ping number: {file_info['start_ping']}")
print(f"      Ending global ping number: {file_info['end_ping']}")
print(f"                Number of pings: {file_info['n_pings']}")
print(f"                First ping time: {file_info['start_time']}")
print(f"                 Last ping time: {file_info['end_time']}")
print(f"                      File size: {size_mib:.2f} MiB")
print()
print(f"             Number of channels: {n_channels}")
print()

for chan in file_info['channels']:

    print("")
    print(f"                         Channel ID: {chan}")
    print(f"                 Transceiver number: {file_info['configuration'][chan]['transceiver_number']}")
    print(f"                   Transceiver type: {file_info['configuration'][chan]['transceiver_type']}")
    print(f"                   Transceiver name: {file_info['configuration'][chan]['transceiver_name']}")
    print(f"          Transceiver serial number: {file_info['configuration'][chan]['serial_number']}")
    print(f"                     Market segment: {file_info['configuration'][chan]['market_segment']}")
    print(f"               Multiplexing enabled: {file_info['configuration'][chan]['multiplexing']}")
    print(f"             Transceiver IP address: {file_info['configuration'][chan]['ip_address']}")
    print(f"                   Ethernet address: {file_info['configuration'][chan]['ethernet_address']}")
    print(f"                         Pulse form: {file_info['configuration'][chan]['pulse_form']}")
    print(f"                       Channel mode: {file_info['configuration'][chan]['channel_mode']}")
    print(f"              Transceiver impedance: {file_info['configuration'][chan]['impedance']}")
    print(f"           Rx sample frequency (Hz): {file_info['configuration'][chan]['rx_sample_frequency']}")
    print(f"     Hardware channel configuration: {file_info['configuration'][chan]['hw_channel_configuration']}")
    print(f"   Maximum transceiver Tx power (W): {file_info['configuration'][chan]['max_tx_power_transceiver']}")
    print(f"          Embedded software version: {file_info['configuration'][chan]['transceiver_software_version']}")
    print(f"           FPGA TX firmware version: {file_info['configuration'][chan]['fpga_tx_firmware_ver']}")
    print(f"           FPGA RX firmware version: {file_info['configuration'][chan]['fpga_rx_firmware_ver']}")
    print(f"                    Transducer name: {file_info['configuration'][chan]['transducer']['transducer_name']}")
    print(f"             Transducer custom name: {file_info['configuration'][chan]['transducer']['transducer_custom_name']}") 
    print(f"           Transducer serial number: {file_info['configuration'][chan]['transducer']['transducer_serial_number']}")
    print(f"               Transducer beam type: {file_info['configuration'][chan]['transducer']['transducer_beam_type']}")
    print(f"               Transducer frequency: {file_info['configuration'][chan]['transducer']['transducer_frequency']}")
    print(f"       Transducer frequency minimum: {file_info['configuration'][chan]['transducer']['transducer_frequency_minimum']}")
    print(f"       Transducer frequency maximum: {file_info['configuration'][chan]['transducer']['transducer_frequency_maximum']}")
    print(f"    Maximum transducer Tx power (W): {file_info['configuration'][chan]['transducer']['max_tx_power_transducer']}")
    print(f"                Transducer mounting: {file_info['configuration'][chan]['transducer']['transducer_mounting']}")
    print(f"             Transducer orientation: {file_info['configuration'][chan]['transducer']['transducer_orientation']}")
    print(f"                Transducer offset X: {file_info['configuration'][chan]['transducer']['transducer_offset_x']}")
    print(f"                Transducer offset Y: {file_info['configuration'][chan]['transducer']['transducer_offset_y']}")
    print(f"                Transducer offset Z: {file_info['configuration'][chan]['transducer']['transducer_offset_z']}")
    print(f"                 Transducer alpha X: {file_info['configuration'][chan]['transducer']['transducer_alpha_x']}")
    print(f"                 Transducer alpha Y: {file_info['configuration'][chan]['transducer']['transducer_alpha_y']}")
    print(f"                 Transducer alpha Z: {file_info['configuration'][chan]['transducer']['transducer_alpha_z']}")
    print()

print()
