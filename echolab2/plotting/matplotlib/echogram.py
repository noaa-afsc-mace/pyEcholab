# coding=utf-8

#     National Oceanic and Atmospheric Administration (NOAA)
#     Alaskan Fisheries Science Center (AFSC)
#     Resource Assessment and Conservation Engineering (RACE)
#     Midwater Assessment and Conservation Engineering (MACE)

#  THIS SOFTWARE AND ITS DOCUMENTATION ARE CONSIDERED TO BE IN THE PUBLIC DOMAIN
#  AND THUS ARE AVAILABLE FOR UNRESTRICTED PUBLIC USE. THEY ARE FURNISHED "AS
#  IS."  THE AUTHORS, THE UNITED STATES GOVERNMENT, ITS INSTRUMENTALITIES,
#  OFFICERS, EMPLOYEES, AND AGENTS MAKE NO WARRANTY, EXPRESS OR IMPLIED,
#  AS TO THE USEFULNESS OF THE SOFTWARE AND DOCUMENTATION FOR ANY PURPOSE.
#  THEY ASSUME NO RESPONSIBILITY (1) FOR THE USE OF THE SOFTWARE AND
#  DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL SUPPORT TO USERS.

"""

| Developed by:  Rick Towler   <rick.towler@noaa.gov>
| National Oceanic and Atmospheric Administration (NOAA)
| Alaska Fisheries Science Center (AFSC)
| Midwater Assessment and Conservation Engineering Group (MACE)
|
| Author:
|       Rick Towler   <rick.towler@noaa.gov>
| Maintained by:
|       Rick Towler   <rick.towler@noaa.gov>

"""

import numpy as np
from matplotlib import figure, axes
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.colors import LinearSegmentedColormap, Colormap
import echolab2.processing.grid


class Echogram(object):
    """This class generates echogram plots.

    The Echogram class provides basic plotting functions to display
    echolab2 data objects using matplotlib.

    In matplotlib there are 3 methods for displaying 2d bitmap data: imshow,
    pcolor, and pcolormesh. For echograms of any size, imshow is the only real
    option as the other two methods are quite slow. A major limitation of
    imshow for this application is that it assumes data are on a regular grid 
    which is only partially true for echosounder data. Y axis data are regular,
    but x axis data are not.
    
    Methods that plot data on top of an echogram like plot_line, will adjust
    the x axis of the data to align with the underlying sample data as displayed
    by imshow so visually things will look correct, but understand that the
    x values used for plotting will be slightly different than your recorded
    ping times.

    If you want to plot items directly on the echogram canvas, you will need
    to call the adjust_xdata method to adjust your x values prior to plotting.

    Attributes:
        figure: A figure for displaying echogram plots.
        axes: An axis of the plot.
        data_object: A PingData or ProcessedData object.
        data_attribute: An attribute of the PingData object.
        threshold: Plot upper and lower display thresholds.
        _simrad_color_table: A list of tuples based on the default SIMRAD
            EK500 color table, plus grey for no data.
        _simrad_cmap: Defines a simrad colormap for the plot.
        cmap (str or Colormap instance):  A string or Colormap instance for
            the plot.
        y_label_attribute: A string of the vertical axis attribute.  Either
            "range" or "depth".
    """

    def __init__(self, mpl_container, data_object, data_attribute=None,
                 threshold=None, cmap=None, frequency=None, grid=None,
                 y_axis=None, x_axis=None):
        """Initializes the Echogram class.

        Args:
            mpl_container: A figure or axis (or None if creating a figure
                           from scratch).
            data_object (object): A RawData or ProcessedData object.
            data_attribute (str): Optional. Specify the data attribute name of the
                                  attribute you want to plot. This is only required
                                  if you are passing something other than a
                                  ProcessedData or RawData object as the data_object.
            threshold (list): lower and upper thresholds in dB. Default is to scale
                              display thresholds to the min/max data values.
            cmap (str or Colormap instance):  A colormap for the plot. The default is
                               the Simrad EK500 color table.
            frequency (float): When plotting Sv(f) data objects, you can specify the
                               specific frequency to plot. If you are plotting Sv(f)
                               and do not specify a frequency, power will be averaged
                               across the band
            grid (object): Optional. An instance of echolab2.processing.grid object
                           that will be plotted. If no grid object is provided, a
                           default grid will be applied.
            y_axis (str): Optional. Specify the Y axis units to use. Set to 'range'
                          to use range as the Y axis, 'depth' for depth, and 'sample'
                          for sample number. When None is specified, we'll look for
                          range first and use that if available, then depth. If
                          neither are available then sample number will be used.
                          Default is None.
            x_axis (str): Optional. Specify the X axis units to use. Set to 'ping_time'
                          to use ping time as the X axis or 'ping_number' to use the
                          ping number. Note that ping number is relative to the data
                          being plotted (first ping is ping 1.) When None is specified,
                          we'll look for  ping_time first and use that if available. If
                          ping_time is not available, ping_number will be used.

        Raises:
            ValueError: a figure or subplot wasn't passed in.
        """

        # Determine what matplotlib container we have.  We must be passed
        # a figure or axes or None.  If we're passed a figure, we assume we're
        # rendering to the "active" axes.  If passed None, we create a figure.
        if mpl_container is None:
            self.figure = figure()
            self.axes = self.figure.gca()
        elif isinstance(mpl_container, figure.Figure):
            self.figure = mpl_container
            self.axes = mpl_container.gca()
        elif isinstance(mpl_container, axes.Subplot):
            self.figure = None
            self.axes = mpl_container
        else:
            raise ValueError("You must pass either a matplotlib figure or " +
                             "subplot specifying where the echogram will be "
                             "rendered.")

        # store some attributes
        self.data_object = data_object
        self.frequency = frequency
        self.threshold = threshold
        self.imshow_x_axis = None
        
        #  store the grid, if provided
        if isinstance(grid, echolab2.processing.grid.grid):
            self.grid = grid
        else:
            self.grid = None

        # The data attribute is only required when plotting ping_data objects.
        # If not provided, we assume we're plotting a processed_data object
        # whose data attribute is "data".
        if data_attribute:
            if not hasattr(self.data_object, data_attribute):
                raise ValueError("The data_attribute : " + data_attribute +
                                 "does not exist in the data_object provided.")
            else:
                data_attribute = getattr(self.data_object, data_attribute)
                self.echogram_data = data_attribute
        else:
            self.echogram_data = data_object.data
        
        # If we're working with FM data, power average across the band if no frequency is provided
        if len(np.unique(self.data_object.frequency)) > 1: 
            if self.frequency is not None:
                self.echogram_data = self.echogram_data[(np.abs(self.data_object.frequency -
                        self.frequency)).argmin()]
            else:
                self.echogram_data = 10*np.log10(np.nanmean(10**(self.echogram_data/10),axis=0))

        # Determine the vertical and horizontal axis attributes.

        # If the y axis is not specified, we look for an appropriate axes,
        # starting with range, then depth, then assume sample.
        if y_axis is None:
            if hasattr(self.data_object, 'range'):
                self.y_label_attribute = 'range'

            elif hasattr(self.data_object, 'depth'):
                self.y_label_attribute = 'depth'
            else:
                # We don't have range or depth attributes?, assume sample.
                self.y_label_attribute = 'sample'
        else:
            # If the y axis is specified, make sure we have that attribute
            if hasattr(self.data_object, y_axis):
                self.y_label_attribute = y_axis
            else:
                # Either sample was specified, or range/depth was specified but
                # we don't have that attribute so default to sample.
                self.y_label_attribute = 'sample'

        # Same deal with the X axis. If nothing is specified, default to ping
        # time if available. If not, use ping number.
        if x_axis is None:
            if hasattr(self.data_object, 'ping_time'):
                self.x_label_attribute = 'ping_time'
            else:
                self.x_label_attribute = 'ping_number'
        else:
            # The X axis has been specified
            if hasattr(self.data_object, x_axis):
                    self.x_label_attribute = x_axis
            else:
                # Either ping number was specified or ping time was specified
                # but is not available. Default to ping number.
                self.x_label_attribute = 'ping_number'

        # Set the default SIMRAD EK500 color table plus grey for NoData.
        self._simrad_color_table = [(1, 1, 1),
                                    (0.6235, 0.6235, 0.6235),
                                    (0.3725, 0.3725, 0.3725),
                                    (0, 0, 1),
                                    (0, 0, 0.5),
                                    (0, 0.7490, 0),
                                    (0, 0.5, 0),
                                    (1, 1, 0),
                                    (1, 0.5, 0),
                                    (1, 0, 0.7490),
                                    (1, 0, 0),
                                    (0.6509, 0.3255, 0.2353),
                                    (0.4705, 0.2353, 0.1568)]
        self._simrad_cmap = (LinearSegmentedColormap.from_list
                             ('Simrad', self._simrad_color_table))
        self._simrad_cmap.set_bad(color='grey')

        if cmap is None:
            self.cmap = self._simrad_cmap
        else:
            if type(cmap) is str:
                self.cmap = plt.get_cmap(cmap).copy()
            else:
                self.cmap = cmap.copy()
        self.cmap.set_bad(color='grey')

        self.update()


    def set_grid(self, grid_obj, update=True):
        """Sets a grid. The grid must be an echolab2.processing.grid object. Grid
        line locations, color, thickness, and line type are obtained from the
        grid object

        Args:
            grid_obj (obj): An echolab2.processing.grid object
            update (bool): Set to True to update the plot.
        """

        if isinstance(grid_obj, echolab2.processing.grid.grid):
            self.grid=grid_obj

        if update:
            self.update()
            

    def set_colormap(self, colormap, bad_data='grey', update=True):
        """Sets a colormap.

        Args:
            colormap (str): A color scheme for plotting.
            bad_data (str): The color to use for bad data values.  Needs to
                be a matplotlib color name.
            update (bool): Set to True to update the plot.
        """

        if isinstance(colormap, str):
            colormap = Colormap(colormap)
        self.cmap = colormap
        if bad_data:
            self.cmap.set_bad(color=bad_data)
        if update:
            self.update()


    def set_threshold(self, threshold=[-70,-34], update=True):
        """Sets a upper and lower threshold.

        Args:
            threshold (list): A two element list of upper and lower thresholds.
            update (bool): Set to True to update the plot.
        """

        if threshold:
            self.threshold = threshold
        else:
            self.threshold = None

        if update:
            self.update()


    def adjust_xdata(self, x_data):
        """Shifts x axis data when displaying data with imshow.

        For various reasons, echosounder ping times are not on a regular grid.
        Matplotlib's imshow forces data onto a regular grid which is OK for
        displaying sample data, but we run into issues when plotting data on top
        of an imshow echogram since the x axes no longer lines up with the samples.

        This method computes new x axis values that shifts the x axis data such
        that the values are centered on the sample pixel's x axis. Methods such
        as plot_line will apply this correction for you, but if you want to plot
        items diurectly onto the echogram canvas, you will need to adjust the x
        values using this method to get them to line up.

        This correction is required when plotting data using time for the
        alongtrack axis.

        Note that corrections do not need to be applied to the Y axis as samples
        are always on a regular range/depth grid.

        Args:
            x_data (array): A numpy array of datetime64 values that represent the
                            x axis data to be shifted.

        Returns:
            An array, adj_x, of shifted data.
        """

        if self.imshow_x_axis is not None and self.x_label_attribute == 'ping_time':
            xticks = np.interp(x_data.astype(np.float64),
                    self.data_object.ping_time.astype(np.float64), self.imshow_x_axis)
        else:
            # if the x axis is not ping time it must be ping number which
            # doesn't require correction.
            xticks = x_data

        return xticks


    def add_colorbar(self, fig, units='dB'):
        """
        add_colorbar adds a colorbar on the right side of the echogram in the provided
        figure.
        """
        cb = fig.colorbar(self.axes_image)
        cb.set_label(units)


    def plot_integration_results(self, int_obj, color=[0.7,0,0], size=10):
        """Plots nasc values output from echolab2.processing.integration class

        Args:
            int_obj (object): An instance of a PyEcholab integration results object.
            color (list): Defines the color of the text
            size (float): Defines the size of the text
        """
        
        self.grid = int_obj.grid
        if self.grid.layer_axis == 'depth':
            vert_vals =  self.grid.depth_middle
        else:
            vert_vals =  self.grid.range_middle

        # If we're working with FM data, average across the results if no frequency is provided
        # This probably doesn't make a lot of sense, so you shouldn't be doing this.
        if len(np.unique(self.data_object.frequency)) > 1: 
            if self.frequency is not None:
                #  frequency provided, extract nasc to plot
                freq_idx = (np.abs(self.data_object.frequency - self.frequency)).argmin()
                int_data = int_obj.nasc[freq_idx,:,:]
            else:
                #  no frequency specified, average nasc across band
                int_data = 10*np.log10(np.nanmean(10**(int_obj.nasc/10), axis=0))
        else:
            # results from CW data can be plotted directly
            int_data = int_obj.nasc
        
        #  adjust the x axis values if needed
        grid_mid_x = self.adjust_xdata(self.grid.time_middle.astype(np.float64))
        
        #  now loop thru the grid cells and plot the text
        for i, int_middle in enumerate(grid_mid_x):
                for v, layer_middle in enumerate(vert_vals):
                    self.axes.text(int_middle, layer_middle, f"{int_data[i,v]:.2f}", color=color,
                            fontsize=size, verticalalignment='center', horizontalalignment='center')


    def plot_line(self, line_obj, color=None, linestyle=None, linewidth=None):
        """Plots an echolab2 Line object on the echogram.
        
        NOTE! - Since imshow forces data onto a regular grid but ping times
        are usually not regular, this method will adjust the line's x axis
        values so they align visually with the other data on the 
        

        Args:
            line_obj (Line obj.): An instance of a PyEcholab Line object.
            color (list): Defines the color of the plotted line.
            linestyle (str): Defines the style of the plotted line.
            linewidth (float): Defines the width of the plotted line.
        """

        if color is None:
            color = line_obj.color
        if linestyle is None:
            linestyle = line_obj.linestyle
        if linewidth is None:
            linewidth = line_obj.linewidth

        # Get the line's horizontal axis. Adjust x values if needed
        if self.x_label_attribute == 'ping_time':
            # when plotting ping times, they must be adjusted
            xticks = self.adjust_xdata(line_obj.ping_time.astype(np.float64))
        else:
            xticks = np.interp(line_obj.ping_time.astype(np.float64),
                    self.data_object.ping_time.astype(np.float64), self.xticks)
    
        # Plot the data.
        self.axes.plot(xticks, line_obj.data, color=color, linestyle=linestyle,
                       linewidth=linewidth, label=line_obj.name)


    def update(self):
        """Updates the plot."""

        # This is a custom tick formatter for datetime64 values as floats.
        def format_datetime(x, pos=None):
            try:
                dt = x.astype('datetime64[ms]').astype('object')
                tick_label = dt.strftime("%H:%M:%S")
            except:
                tick_label = ''

            return tick_label

        # Get the thresholds if we have been given one.
        if self.threshold:
            threshold = self.threshold
        else:
            # Just use the min/max if we don't have thresholds.
            threshold = [np.nanmin(self.echogram_data),
                         np.nanmax(self.echogram_data)]

        # Transform the data so it looks right with imshow.
        echogram_data = np.flipud(np.rot90(self.echogram_data, 1))

        # Determine the vertical extent of the data and the y label.
        if self.y_label_attribute == 'sample':
            # Use sample number for the Y axis
            yticks = np.arange(echogram_data.shape[0]) + 1
            y_label = 'sample'
        elif hasattr(self.data_object, self.y_label_attribute):
            # Use range/depth as the Y axis
            yticks = getattr(self.data_object, self.y_label_attribute)
            y_label = self.y_label_attribute + ' (m)'

        if self.x_label_attribute == 'ping_number':
            # x ticks are ping number - generate ping number vector
            xticks = np.arange(self.data_object.shape[0], dtype='float32') + 1
            self.xticks = xticks
        else:
            # The x ticks are the pings times as serial time.
            xticks = self.data_object.ping_time.astype(np.float64)
        
        # Plot the sample data.

        # When using imshow with the extents keyword, you must pad the
        # extents a half pixel width on either side to ensure that the pixel
        # centers align with our underlying grid. Compute the half pixel widths
        dx = ((xticks[-1] - xticks[0]) / (echogram_data.shape[1] - 1)) / 2.0
        dy = ((yticks[0] - yticks[-1]) / (echogram_data.shape[0] - 1)) / 2.0
        
        # and compute the extents
        self.extents = [xticks[0] - dx, xticks[-1] + dx, yticks[-1]- dy, yticks[0] + dy]
        
        # imshow assumes data are on a regular grid but ping intervals are
        # almost never regular. Generate the x axis imshow will impose on the sample
        # data so we can use it to adjust elements plotted on the echogram
        self.imshow_x_axis = np.linspace(xticks[0], xticks[-1], self.echogram_data.shape[0])
        
        #  finally, plot the echogram using imshow
        self.axes_image = self.axes.imshow(
            echogram_data, cmap=self.cmap, vmin=threshold[0], vmax=threshold[1],
            aspect='auto', interpolation='none',  origin='upper',extent=self.extents)

        #  if there is a echolab grid object associated with this echogram render it
        if self.grid:
            if self.x_label_attribute == 'ping_time':
                grid_xdata = np.append(self.grid.time_start,
                        self.grid.time_end[-1]).astype(np.float64)
            else:
                grid_xdata = np.append(self.grid.ping_start,
                        self.grid.ping_end[-1]).astype(np.float64)
                
            grid_xdata = self.adjust_xdata(grid_xdata)

            if self.y_label_attribute == 'range':
                grid_ydata = np.append(self.grid.range_start,
                        self.grid.range_end[-1]).astype(np.float64)
            elif self.y_label_attribute == 'depth':
                grid_ydata = np.append(self.grid.range_start,
                        self.grid.range_end[-1]).astype(np.float64)
            else:
                grid_ydata = np.append(self.grid.sample_start,
                        self.grid.sample_end[-1]).astype(np.float64)
                    
            #  set the grid axes locations
            self.axes.set_xticks(grid_xdata)
            self.axes.set_yticks(grid_ydata)
            
            #  set the grid properties based on the grid object attributes
            self.axes.grid(True, color=self.grid.color, 
                    linewidth=self.grid.linewidth, linestyle=self.grid.linestyle)
        else:
            # No grid. Apply the defaule grid.
            self.axes.grid(True, color=[0,0,0])

        if self.x_label_attribute == 'ping_time':
            # Set our custom x-axis formatter.
            self.axes.xaxis.set_major_formatter(ticker.FuncFormatter(
                format_datetime))

            # Set the x axis label to the month-day-year of the first.
            # datetime64 we have in the data. This will fail if there are no
            # valid times so we'll not label the axis if that happens.
            try:
                x = self.axes.get_xticks()[0]
                dt = x.astype('datetime64[ms]').astype('object')

                x_label = dt.strftime("%m-%d-%Y")
            except:
                x_label = 'Ping Time'
        else:
            x_label = 'Ping Number'

        #  set the axes labels
        self.axes.set_xlabel(x_label)
        self.axes.set_ylabel(y_label)

        # There seems to be a bug in matplotlib where extra white space is
        # added to the figure around the echogram if other elements are
        # plotted on the echogram even if their data ranges fall *within* the
        # echogram.  This doesn't happen when we use ping number as the x axis,
        # but we don't want to use ping number because everything we plot on
        # the echogram would have to be mapped to ping number.  Weirdly, what
        # solves the issue is calling set_xlim/set_ylim without arguments
        # which should only return the current limits, but seems to fix this
        # issue.
        self.axes.set_xlim()
        self.axes.set_ylim()

