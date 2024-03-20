# coding=utf-8

#     National Oceanic and Atmospheric Administration (NOAA)
#     Alaskan Fisheries Science Center (AFSC)
#     Resource Assessment and Conservation Engineering (RACE)
#     Midwater Assessment and Conservation Engineering (MACE)

#  THIS SOFTWARE AND ITS DOCUMENTATION ARE CONSIDERED TO BE IN THE PUBLIC DOMAIN
#  AND THUS ARE AVAILABLE FOR UNRESTRICTED PUBLIC USE. THEY ARE FURNISHED "AS IS."
#  THE AUTHORS, THE UNITED STATES GOVERNMENT, ITS INSTRUMENTALITIES, OFFICERS,
#  EMPLOYEES, AND AGENTS MAKE NO WARRANTY, EXPRESS OR IMPLIED, AS TO THE USEFULNESS
#  OF THE SOFTWARE AND DOCUMENTATION FOR ANY PURPOSE. THEY ASSUME NO RESPONSIBILITY
#  (1) FOR THE USE OF THE SOFTWARE AND DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL
#  SUPPORT TO USERS.
"""
| Developed by:  Rick Towler   <rick.towler@noaa.gov>
| National Oceanic and Atmospheric Administration (NOAA)
| Alaska Fisheries Science Center (AFSC)
| Midwater Assesment and Conservation Engineering Group (MACE)
|
| Author:
|       Rick Towler   <rick.towler@noaa.gov>
| Maintained by:
|       Rick Towler   <rick.towler@noaa.gov>
"""

import numpy as np
from ...ping_data import ping_data


class motion_data(object):
    '''
    The motion_data class stores pitch, roll, heave, and heading data
    and provides a method to interpolate the data to your ping times.

    3/2024 - extended to support MRU1 datagrams based on the KM Binary sensor
             datagram, version 1. More info here:

             https://www3.mbari.org/products/mbsystem/formatdoc/KongsbergKmall/EMdgmFormat_RevH/html/kmBinary.html

    '''

    CHUNK_SIZE = 500

    def __init__(self):

        # Create a counter to keep track of the number of datagrams.
        self.n_raw = 0

        #  default to an MRU0 style object
        self.has_MRU1 = False

        # Create arrays to store MRU0 data
        self.times = np.empty(motion_data.CHUNK_SIZE, dtype='datetime64[ms]')
        self.heave = np.empty(motion_data.CHUNK_SIZE, dtype='f')
        self.pitch = np.empty(motion_data.CHUNK_SIZE, dtype='f')
        self.roll = np.empty(motion_data.CHUNK_SIZE, dtype='f')
        self.heading = np.empty(motion_data.CHUNK_SIZE, dtype='f')


    def add_datagram(self, datagram):
        """
        Add motion datagram data. This method will extract the values from the provided motion
        datagram and store it within this object.  This method will automatically convert from
        storing MRU0 data to MRU1 if passed an MRU1 datagram.

        Args:
            time (datetime64)
            datagram (dict)
        """

        # Check if this datagram has the same time as the previous datagram.
        # This simply filters replicate data when used with the EK60 class.
        if self.n_raw > 0 and self.times[self.n_raw-1] ==  datagram['timestamp']:
            # We already have this motion datagram stored.
            return

        # Increment datagram counter.
        self.n_raw += 1

        # Check if we need to resize our arrays.
        if self.n_raw > self.times.shape[0]:
            self._resize_arrays(self.times.shape[0] + motion_data.CHUNK_SIZE)

        # Add the fields common to MRU0 and MRU1
        self.times[self.n_raw-1] = datagram['timestamp']
        self.heave[self.n_raw-1] = datagram['heave']
        self.pitch[self.n_raw-1] = datagram['pitch']
        self.roll[self.n_raw-1] = datagram['roll']
        self.heading[self.n_raw-1] = datagram['heading']

        #  check if we've been passed a MRU1 datagram
        if 'status_word' in datagram:

            #  yes, this is an MRU1 datagram, check if we have attributes to store
            #  this data, if not, add the attributes to store this data
            if not self.has_MRU1:
                self._convert_to_MRU1()

            #  add the remainder of the MRU1 fields
            self.status[self.n_raw-1] = datagram['status_word']
            self.latitude[self.n_raw-1] = datagram['latitude']
            self.longitude[self.n_raw-1] = datagram['longitude']
            self.ellipsoid_height[self.n_raw-1] = datagram['ellipsoid_height']
            self.roll_rate[self.n_raw-1] = datagram['roll_rate']
            self.pitch_rate[self.n_raw-1] = datagram['pitch_rate']
            self.yaw_rate[self.n_raw-1] = datagram['yaw_rate']
            self.north_velocity[self.n_raw-1] = datagram['velocity_north']
            self.east_velocity[self.n_raw-1] = datagram['velocity_east']
            self.down_velocity[self.n_raw-1] = datagram['velocity_down']
            self.latitude_error[self.n_raw-1] = datagram['latitude_error']
            self.longitude_error[self.n_raw-1] = datagram['longitude_error']
            self.height_error[self.n_raw-1] = datagram['height_error']
            self.roll_error[self.n_raw-1] = datagram['roll_error']
            self.pitch_error[self.n_raw-1] = datagram['pitch_error']
            self.heading_error[self.n_raw-1] = datagram['heading_error']
            self.heave_error[self.n_raw-1] = datagram['heave_error']
            self.north_acceleration[self.n_raw-1] = datagram['accel_north']
            self.east_acceleration[self.n_raw-1] = datagram['accel_east']
            self.down_acceleration[self.n_raw-1] = datagram['accel_down']
            self.delayed_heave_utc_second[self.n_raw-1] = datagram['heave_delay_secs']
            self.delayed_heave_utc_nanoseconds[self.n_raw-1] = datagram['heave_delay_usecs']
            self.delayed_heave_m[self.n_raw-1] = datagram['heave_delay_m']


    def interpolate(self, p_data, attributes=None):
        """
        interpolate returns the requested motion data interpolated to the ping times
        that are present in the provided ping_data object.

            p_data is a ping_data object that contains the ping_time vector
                    to interpolate to. You can also simply pass an array of datetime64
                    objects to interpolate to if you are not working with processed_data
                    objects.
            attributes is a string or list of strings specifying the motion attribute(s)
                    to interpolate and return. If this argument is omitted or None, the
                    following attributes will be returned:

                        MRU0 data: 'heave', 'pitch', 'roll', 'heading'
                        MRU1 data:'latitude', 'longitude', 'heave', 'pitch', 'roll', 'heading'

        Returns a dictionary of numpy arrays keyed by attribute name that contain the
        interpolated data for that attribute.
        """
        # Create the dictionary to return
        out_data = {}

        # Check if we're been given specific attributes to interpolate. If not, we
        # interpolate what I am defining as
        if attributes is None:
            if self.has_MRU1:
                attributes = ['latitude', 'longitude', 'heave', 'pitch', 'roll', 'heading']
            else:
                attributes = ['heave', 'pitch', 'roll', 'heading']
        elif isinstance(attributes, str):
            # We have a string, put it in a list
            attributes = [attributes]

        #  check if the times are to be grabbed from a ping_data object
        if isinstance(p_data, ping_data):
            new_times = p_data.ping_time
        else:
            # If not a ping_data object, assume we've just been given times
            new_times = p_data

        # Work through the attributes and interpolate
        for attribute in attributes:
            try:
                # Interpolate this attribute using the time vector in the
                # provided ping_data object
                out_data[attribute] = np.interp(new_times.astype('d'),
                        self.times.astype('d'), getattr(self, attribute),
                        left=0, right=0)
            except:
                # Provided attribute doesn't exist
                out_data[attribute] = None

        return (attributes, out_data)


    def get_indices(self, start_time=None, end_time=None, time_order=True):
        """
        Return index of data contained in specified time range.

        get_indices returns an index array containing the indices contained
        in the range defined by the times provided. By default the indexes
        are in time order.

        Args:
            start_time is a datetime or datetime64 object defining the starting
                time of the data to return. If None, the start time is the
                earliest time.
            end_time is a datetime or datetime64 object defining the ending time
                of the data to return. If None, the end time is the latest time.
            time_order (bool): Control whether if indexes are returned in time
                order (True) or not.

        Returns: Index array containing indices of data to return.

        """
        #  Ensure that we have times to work with.
        if start_time is None:
            start_time = np.min(self.times)
        if end_time is None:
            end_time = np.max(self.times)

        # Sort time index if returning time ordered indexes.
        if time_order:
            primary_index = self.times.argsort()
        else:
            primary_index = self.times

        # Determine the indices of the data that fall within the time span
        # provided.
        mask = self.times[primary_index] >= start_time
        mask = np.logical_and(mask, self.times[primary_index] <= end_time)

        #  and return the indices that are included in the specified range
        return primary_index[mask]


    def trim(self):
        """
        Trim arrays to proper size after all data are added.

        trim is called when one is done adding data to the object. It
        removes empty elements of the data arrays.
        """

        self._resize_arrays(self.n_raw)


    def _resize_arrays(self, new_size):
        """
        Resize arrays if needed to hold more data.

        _resize_arrays expands our data arrays and is called when said arrays
        are filled with data and more data need to be added.

        Args:
            new_size (int): New size for arrays, Since these are all 1d
            arrays the value is simply an integer.

        """

        #  extend the common fields
        self.times = np.resize(self.times,(new_size))
        self.pitch = np.resize(self.pitch,(new_size))
        self.roll = np.resize(self.roll,(new_size))
        self.heading = np.resize(self.heading,(new_size))
        self.heave = np.resize(self.heave,(new_size))

        #  if we have MRU1 data, extend the MRU1 datagrams
        if self.has_MRU1:
            self.status = np.resize(self.status,(new_size))
            self.latitude = np.resize(self.latitude,(new_size))
            self.longitude = np.resize(self.longitude,(new_size))
            self.ellipsoid_height = np.resize(self.ellipsoid_height,(new_size))
            self.roll_rate = np.resize(self.roll_rate,(new_size))
            self.pitch_rate = np.resize(self.pitch_rate,(new_size))
            self.yaw_rate = np.resize(self.yaw_rate,(new_size))
            self.north_velocity = np.resize(self.north_velocity,(new_size))
            self.east_velocity = np.resize(self.east_velocity,(new_size))
            self.down_velocity = np.resize(self.down_velocity,(new_size))
            self.latitude_error = np.resize(self.latitude_error,(new_size))
            self.longitude_error = np.resize(self.longitude_error,(new_size))
            self.height_error = np.resize(self.height_error,(new_size))
            self.roll_error = np.resize(self.roll_error,(new_size))
            self.pitch_error = np.resize(self.pitch_error,(new_size))
            self.heading_error = np.resize(self.heading_error,(new_size))
            self.heave_error = np.resize(self.heave_error,(new_size))
            self.north_acceleration = np.resize(self.north_acceleration,(new_size))
            self.east_acceleration = np.resize(self.east_acceleration,(new_size))
            self.down_acceleration = np.resize(self.down_acceleration,(new_size))
            self.delayed_heave_utc_second = np.resize(self.delayed_heave_utc_second,(new_size))
            self.delayed_heave_utc_nanoseconds = np.resize(self.delayed_heave_utc_nanoseconds,(new_size))
            self.delayed_heave_m = np.resize(self.delayed_heave_m,(new_size))


    def _convert_to_MRU1(self):
        """
        _convert_to_MRU1 will extend the class to store all of the fields in the MRU1
        motion datagram introduced around EK80 version 23.x. The MRU1 datagram is based on
        and follows the format of the KM Binary sensor datagram, version 1. More info here:

        https://www3.mbari.org/products/mbsystem/formatdoc/KongsbergKmall/EMdgmFormat_RevH/html/kmBinary.html

        """

        #  when we extend the class for MRU1 data, create arrays filled with NaNs in case the object is
        #  extended after already reading MRU0 style data.
        self.status = np.full(self.times.size, 0, dtype='i4')
        self.latitude = np.full(self.times.size, np.nan, dtype='f8')
        self.longitude = np.full(self.times.size, np.nan, dtype='f8')
        self.ellipsoid_height = np.full(self.times.size, np.nan, dtype='f')
        self.roll_rate = np.full(self.times.size, np.nan, dtype='f')
        self.pitch_rate = np.full(self.times.size, np.nan, dtype='f')
        self.yaw_rate = np.full(self.times.size, np.nan, dtype='f')
        self.north_velocity = np.full(self.times.size, np.nan, dtype='f')
        self.east_velocity = np.full(self.times.size, np.nan, dtype='f')
        self.down_velocity = np.full(self.times.size, np.nan, dtype='f')
        self.latitude_error = np.full(self.times.size, np.nan, dtype='f')
        self.longitude_error = np.full(self.times.size, np.nan, dtype='f')
        self.height_error = np.full(self.times.size, np.nan, dtype='f')
        self.roll_error = np.full(self.times.size, np.nan, dtype='f')
        self.pitch_error = np.full(self.times.size, np.nan, dtype='f')
        self.heading_error = np.full(self.times.size, np.nan, dtype='f')
        self.heave_error = np.full(self.times.size, np.nan, dtype='f')
        self.north_acceleration = np.full(self.times.size, np.nan, dtype='f')
        self.east_acceleration = np.full(self.times.size, np.nan, dtype='f')
        self.down_acceleration = np.full(self.times.size, np.nan, dtype='f')
        self.delayed_heave_utc_second = np.full(self.times.size, 0, dtype='i4')
        self.delayed_heave_utc_nanoseconds = np.full(self.times.size, 0, dtype='i4')
        self.delayed_heave_m = np.full(self.times.size, np.nan, dtype='f')

        #  set the has_MRU1 attribute since we now have MRU1 data
        self.has_MRU1 = True


    def __str__(self):
        """
        Reimplemented string method that provides some basic info about the
        motion_data object.

        """

        #  print the class and address
        msg = str(self.__class__) + " at " + str(hex(id(self))) + "\n"

        #  print some more info about the motion_data instance
        if (self.n_raw > 0):
            msg = "{0}       MRU data start time: {1}\n".format(msg, self.times[0])
            msg = "{0}         MRU data end time: {1}\n".format(msg, self.times[self.n_raw-1])
            msg = "{0}       Number of datagrams: {1}\n".format(msg, self.n_raw)
            msg = "{0}  Has extended motion data: {1}\n".format(msg, self.has_MRU1)
        else:
            msg = msg + ("  motion_data object contains no data\n")

        return msg
