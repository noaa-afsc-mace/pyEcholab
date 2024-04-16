# -*- coding: utf-8 -*-
"""This is a simple interface for plotting basic maps using folium.

This is a

"""


import os
import numpy as np
#from geographiclib.geodesic import Geodesic
#import branca.colormap
#from branca.element import Element
import webbrowser
import folium
from folium.features import DivIcon
from folium.plugins import GroupedLayerControl, BeautifyIcon, MousePosition

class folium_map(object):

    def __init__(self, tile_server=None, tile_server_attr=None, min_zoom=None,
            max_zoom=None, initial_zoom=None, title=None):
        super(folium_map, self).__init__()

        self.layers = []

        if tile_server is None:
            tile_server = ('https://server.arcgisonline.com/ArcGIS/rest/services/' +
                    'Ocean/World_Ocean_Base/MapServer/tile/{z}/{y}/{x}')
        if tile_server_attr is None:
            tile_server_attr = ('Tiles &copy; Esri &mdash; Sources: GEBCO, NOAA, ' +
                    'CHS, OSU, UNH, CSUMB, National Geographic, DeLorme, NAVTEQ, and Esri')
        if min_zoom is None:
            min_zoom = 5
        if max_zoom is None:
            max_zoom = 13
        if initial_zoom is None:
            initial_zoom = 8

        #  create the base map
        self.map = folium.Map(tiles=tile_server, attr=tile_server_attr,
            zoom_start=initial_zoom, max_zoom=max_zoom, min_zoom=min_zoom)

        if title is not None:
            #  add the title
            #title = ('<h2 style="position:absolute;z-index:100000;left:34vw" ><b>My Folium Map</b></h2>')
            self.map.get_root().html.add_child(folium.Element(title))


    def add_layer(self, layer_name):

        try:
            fg_idx = self.layers.index(layer_name)
            fg = self.layers[fg_idx]
        except:
            #  layer doesn't exist, add it
            fg = folium.FeatureGroup(layer_name)
            self.layers.append(fg)
            fg.add_to(self.map)

        return fg


    def add_polyline(self, line_verts, color, line_weight=1.5, layer=None, tooltip=None):


        if layer is None:
            layer = self.map

        if tooltip is None:
            folium.PolyLine(locations=line_verts, color=color,
                weight=line_weight).add_to(layer)
        else:
            folium.PolyLine(locations=line_verts, color=color,
                weight=line_weight, tooltip=tooltip).add_to(layer)


    def render_map(self, file_name, open=False):


        mapBounds = [[999,999],[-999,-999]]
        for fg in self.layers:
            thisBounds = fg.get_bounds()
            try:
                if mapBounds[0][0] >= thisBounds[0][0]:
                    mapBounds[0][0] = thisBounds[0][0]
                if mapBounds[0][1] >= thisBounds[0][1]:
                    mapBounds[0][1] = thisBounds[0][1]
                if mapBounds[1][0] <= thisBounds[1][0]:
                    mapBounds[1][0] = thisBounds[1][0]
                if mapBounds[1][1] <= thisBounds[1][1]:
                    mapBounds[1][1] = thisBounds[1][1]
            except:
                pass
        if min(min(mapBounds)) > -999 and max(max(mapBounds)) < 999:
            self.map.fit_bounds(mapBounds)

        #  add the position indicator
        MousePosition().add_to(self.map)

        #  save the map HTML file
        file_name = os.path.normpath(file_name)
        self.map.save(file_name)


        if open:
            webbrowser.open(file_name)
