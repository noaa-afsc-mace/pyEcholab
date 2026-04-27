# -*- coding: utf-8 -*-
"""This script demonstrates extracting Lat/Lon data, interpolating it,
and plotting the ship track using Folium.

Individual .raw files tend to cover relatively short time spans thus
plotted tracks are not very impressive unless you throw a lot of files
at this.

Folium generates an interactive map in the browser. Because even simple
maps can exceed the maximum amount of data browsers will accept in an
URL, this example will create an HTML file on disk. The script will try
to automatically load this file using your system's default browser.
If this doesn't work, you will have to load the file yourself, usually
by just clicking on it.

"""

from glob import glob
import numpy as np
from echolab2.instruments import echosounder
from echolab2.plotting.folium_map import folium_map

track_color = '#9402BDFF'
map_file_name = './echolab2_folium_trackline.html'


#  raw data files - define a list of files to read. The lat/lon data will be extracted
#  from these files and plotted. Each file will be plotted as a separate line segment.
#  in this example we use glob to grab all .raw files in the specified directory.
raw_files = glob('./data/ES60/OEX2017/raw/*.raw')

#  create a list to store the lat/lon data from each file
position_data = []

#  create a folium map object and add a layer for our trackline
map = folium_map()
track_layer = map.add_layer("trackline")

#  loop thru the files, extract GPS, interpolate, and add the line segment
#  for the file to the map.
for file in raw_files:

    #  read the raw file and get a reference to the NMEA data
    print("Reading " + file + "...")
    raw_data = echosounder.read(file)
    print(raw_data)
    nmea_data = raw_data.nmea_data

    #  GPS data is recorded in the raw data at the rate that it is received
    #  from the GPS unit, this is typically between 1-10 Hz. We don't want
    #  to plot that many data points so we will interpolate to 10 second
    #  intervals.

    #  first create a vector of times at 10 second intervals that span the
    #  data time range
    new_gps_times = np.arange(nmea_data.nmea_times[0],nmea_data.nmea_times[-1],
            np.timedelta64(10, 's'))

    #  next, use the interpolate method of the nmea_data class to interpolate
    #  the "position" data. "position" is a NMEA data metatype that will extract
    #  lat/lon data from GGA, RMC, and GLL NMEA sentences. See the
    #  echolab2.instruments.util.nmea_data class for more info.
    interp_fields, interp_data = nmea_data.interpolate(new_gps_times,
        'position')

    #  clean up nans since folium doesn't appreciate them
    not_nan = ~np.isnan(interp_data['latitude'])
    interp_data['latitude'] = interp_data['latitude'][not_nan]
    interp_data['longitude'] = interp_data['longitude'][not_nan]

    #  zip up the lat/lon data into a list of lat/lon tuples for plotting
    track_coords = zip(interp_data['latitude'], interp_data['longitude'])

    #  add this track section to the map
    map.add_polyline(track_coords, track_color, layer=track_layer)

#  render the map to file, then try to open it using the platform's default brower.
map.render_map(map_file_name, open=True)
print()
