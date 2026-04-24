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
.. module:: echolab2.ping_data

    :synopsis: Base class for containers use to store data collected
               from fisheries sonar systems.

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

import copy
import warnings
import numpy as np
from scipy.interpolate import interp1d
from scipy.sparse import csr_matrix


class ping_data(object):
    """echolab2.ping_data is the base class for all classes that store time
    based data from fisheries sonar systems.

    This class is not intended to be instantiated by the user. It is a base
    class that defines the common data attributes and methods that the user
    facing classes share.

    Derived classes will add various attributes to this class that store
    either scalar values on a per-ping basis like sound velocity, transmit power
    transducer depth, etc. or attributes that store vector data such as
    sample power and sample angle data.

    One major assumption is that all data stored within an instance of our
    derived classes must exist on the same time grid. It is assumed that the
    index of a specific ping time should map to other attributes collected at
    that time. As a result, all attributes should share the same primary
    dimension.
    """

    def __init__(self):
        """Initializes ping_data class object.

        Creates and sets several internal properties used to store information
        about data and control operation of data processing of ping_data
        object instance. Code is heavily commented to facilitate use.
        """

        # Stores the total number of pings contained in our container.
        self.n_pings = -1

        # The number of samples in the 2d/3d sample arrays.
        self.n_samples = -1

        # A tuple describing the shape of the sample data array in the form
        # (n_pings, n_samples) or (n_pings, n_samples, s_sectors) for complex
        # data. If shape is None, the data arrays have not been allocated.
        self.shape = None

        # Allows the user to specify a dtype for the sample data.  This should
        # be set before any attributes are added.
        self.sample_dtype = 'float32'

        # Data_attributes is an internal list that contains the names of all
        # the class's "data attributes". The echolab2 package uses this
        # attribute list to generalize various functions that manipulate these
        # data.
        #
        # "data attributes" are attributes that store data by ping. They can
        # be 1d, such as ping_time, sample_interval, and transducer_depth. Or
        # they can be 2d, such as power, Sv, angle data, etc.  Echolab2
        # supports attributes that are 1d, 2d, and 3d numpy arrays.  When
        # subclassing, you must extend this list in your __init__ method to
        # contain all of the data attributes of that class that you want to
        # exist at instantiation (attributes can also be added later).

        # For the base class, we only define ping_time which is the only
        # required attribute that all data objects must have.
        self._data_attributes = ['ping_time']

        # Attributes are added using the add_data_attribute method. You can add
        # them manually by appending the name of the new attribute to the
        # _data_attributes dictionary and then setting the attribute using
        # setattr().

        # Object attributes are similar to data attributes except they are not
        # linked to a data axis (time, sample_number, range/depth). Another
        # important difference is that they are not resized when the data arrays
        # are resized. Object attributes can be set to any value or object and
        # do not have to be numpy arrays.

        # Similar to the data attributes you can add and remove object attributes
        # using the add_object_attribute and remove_object_attribute methods.
        self._object_attributes = ['sample_dtype']

        # When writing methods that operate on these data, we will not assume
        # that they exist.  An attribute should only exist if it contains data.


    def add_data_attribute(self, name, data):
        """Adds a "data attribute" to the class.

        Data attributes are attributes that are linked to one of the data
        axes. Data attributes are resized when the data arrays are resized.

        This method first checks if the new attribute shares the same
        dimensions as existing attributes, and if so, appends the attribute
        name to the internal list of data attributes and then creates an
        attribute by that name and sets it to the provided reference.

        Args:
            name (str): The attribute name to be added to the class.
            data (array): A numpy array containing the data you are adding
                as the attribute.

        Raises:
            ValueError: The attribute has a different number of samples than
                the other attributes.
            TypeError: The attribute is not a numpy array.
            ValueError: The attribute has a different number of pings than
                the other attributes.
        """
        # Get the new data's dimensions.
        data_height = -1
        if isinstance(data, np.ndarray):
            if data.ndim == 1:
                data_width = data.shape[0]
            if data.ndim > 1:
                data_width = data.shape[0]
                data_height = data.shape[1]
                # Check if n_samples has been set yet.  If not, set it.
                # Otherwise, check that the dimensions match.
                if self.n_samples < 0:
                    self.n_samples = data_height
                elif self.n_samples != data_height:
                    raise ValueError('Cannot add attribute. New attribute has ' +
                        'a different number of samples than the existing attributes.')
            if data.ndim == 3:
                # Add or update ourself.
                setattr(self, 'n_sectors', data.shape[2])
                self._object_attributes += ['n_sectors']
        else:
            # We only allow numpy arrays as data attributes.
            raise TypeError('Invalid data attribute type. Data attributes must ' +
                'be numpy arrays.')

        # Check if n_pings has been set yet.  If not, set it.  Otherwise,
        # check that the dimensions match.  When checking if dimensions
        # match, we allow a match on the number of pings OR the number of
        # samples since a 1d data attribute can be on either axis.
        if self.n_pings < 0:
            self.n_pings = data_width
        elif self.n_pings != data_width and self.n_samples != data_width:
            raise ValueError('Cannot add attribute. The new attribute has '
                    'a different number of pings or samples than the other attributes.')

        # Add the name to our list of attributes if it doesn't already exist.
        if name not in self._data_attributes:
            self._data_attributes.append(name)

        # Add or update ourself.
        setattr(self, name, data)

        #  update the shape attribute
        self.shape = self._shape()


    def add_object_attribute(self, name, data):
        """Adds a "object attribute" to the class.

        Object attributes are attributes that are not linked to a data axis.
        Object attributes are primarily used by the processed_data class to
        describe general attributes about the object or data like 'data_type'
        or 'is_log'. Object attributes can be any data type and are not
        resized or altered in any way when the data arrays are resized.

        Since the data is not linked to an axis there is no checking
        of dimensions before adding.

        Args:
            name (str): The attribute name to be added to the class.
            data (object): An python object/value

        """

        # Add the name to our list of attributes if it doesn't already exist.
        if name not in self._object_attributes:
            self._object_attributes.append(name)

        # Add or update ourself.
        setattr(self, name, data)


    def remove_object_attribute(self, name):
        """Removes a object attribute from the object.

        Args:
            name (str): The attribute name to be removed from the class.
        """

        # Try to remove the attribute given the name.  Silently fail if the
        # name is not in our list.
        try:
            self._object_attributes.remove(name)
            delattr(self, name)
        except:
            pass


    def remove_data_attribute(self, name):
        """Removes a data attribute from the object.

        Args:
            name (str): The attribute name to be removed from the class.
        """

        # Try to remove the attribute given the name.  Silently fail if the
        # name is not in our list.
        try:
            self._data_attributes.remove(name)
            delattr(self, name)
        except:
            pass


    def replace(self, obj_to_insert, ping_number=None, ping_time=None,
                index_array=None, _ignore_vertical_axes=False, force=False):
        """Replaces the data in this object with the data provided in the
        object to "insert".

        This method inserts data without shifting the existing data, resulting
        in the existing data being overwritten.  You must specify a ping
        number, ping time or provide an index array.  The number of pings
        replaced will be equal to the number of pings in the object you are
        adding.

        Args:
            obj_to_insert (ping_data): An instance of ping_data containing the
                replacement data to insert.
            ping_number (int): The ping number specifying the first ping to
                replace.
            ping_time (datetime): The ping time specifying the first ping to
                replace.
            index_array (array): A numpy array containing the indices of the
                pings you want to replace. Unlike when using a ping number or
                ping time, the pings do not have to be consecutive. When this
                keyword is present, the ping_number and ping_time keywords
                are ignored.
            _ignore_vertical_axes (bool): Controls whether to ignore vertical
                axes, range or depth, when resampling.

        Raises:
            ValueError: ping_number, ping_time or index array not provided.
            TypeError: The object provided isn't an instance of the ping_data
                class.
            TypeError: The frequency of the replacement data does not match
                the frequency of the data to be replaced.
            TypeError: The index array is not a numpy array.
            IndexError: The length of the index_array does not match the
                number of pings in the object providing the replacement data.
        """

        # Check that we have been given an insertion point or index array.
        if ping_number is None and ping_time is None and index_array is None:
            raise ValueError('Either ping_number or ping_time needs to be ' +
                             'defined or an index array needs to be provided ' +
                             'to specify a replacement point.')

        # Make sure that obj_to_insert class matches "this" class.
        if not isinstance(self, obj_to_insert.__class__):
            raise TypeError('The object provided as a source of replacement ' +
                            'pings must be an instance of ' +
                            str(self.__class__))

        # Make sure the data types are the same
        if self.data_type != obj_to_insert.data_type and not force:
            raise TypeError('The object you are inserting  does not have the same data ' +
                    'type as this object. This data type: ' + self.data_type +
                    ' object to insert data type: ' + obj_to_insert.data_type)

        # Make sure that the frequencies match.
        if not force:
            freq_match = self._check_frequencies(obj_to_insert)
            if not freq_match:
                raise ValueError('The object that you are inserting does not have the ' +
                        'same frequency attributes as this object. You cannot insert '+
                        'data with different frequency attributes.')

        # Get information about the shape of the data we're working with.
        my_pings = self.n_pings
        new_pings = obj_to_insert.n_pings
        my_samples = self.n_samples
        new_samples = obj_to_insert.n_samples

        # Determine the index of the replacement point or indices of the
        # pings we're inserting.
        if index_array is None:
            # Determine the index of the replacement point.
            replace_index = self.get_indices(start_time=ping_time,
                    end_time=ping_time, start_ping=ping_number,
                    end_ping=ping_number)[0]

            # Create an index array.
            replace_index = np.arange(new_pings) + replace_index

            # Clamp the index to the length of our existing data.
            replace_index = replace_index[replace_index < my_pings]

        else:
            # An explicit array is provided.  These will be a vector of
            # locations to replace.
            replace_index = index_array

            # Make sure the index array is a numpy array.
            if not isinstance(replace_index, np.ndarray):
                raise TypeError('index_array must be a numpy.ndarray.')

            # If we replace the index with a user provided index, make sure the
            # dimensions of the index and the object to insert match.
            if replace_index.shape[0] != new_pings:
                raise IndexError('The length of the index_array does not '
                        'match the number of pings in the object with ' +
                        'the replacement data.  These dimensions must match.')

        # Check if we need to vertically resize one of the arrays.  If so, we
        # resize the smaller array to the size of the larger array.  It will
        # automatically be padded with NaNs.
        if my_samples < new_samples:
            # Resize our data arrays.  Check if we have a limit on the max
            # number of samples.
            if hasattr(self, 'max_sample_number') and self.max_sample_number:
                # We have the attribute and a value is set.  Check if the new
                # array exceeds our max_sample_count.
                if new_samples > self.max_sample_number:
                    # We have to change our new_samples.
                    new_samples = self.max_sample_number
                    # Vertically trim the array we're inserting.
                    obj_to_insert.resize(new_pings, new_samples)
            # Because we're not inserting, we can set the new vertical sample
            # size here.
            self.resize(my_pings, new_samples)
        elif my_samples > new_samples:
            # Resize the object we're inserting.
            obj_to_insert.resize(new_pings, my_samples)

        # Work through our data properties inserting the data from
        # obj_to_insert.
        for attribute in self._data_attributes:

            # Check if we have data for this attribute.
            if not hasattr(self, attribute):
                # data_obj does not have this attribute, move along.
                continue

            # Get a reference to our data object's attribute.
            data = getattr(self, attribute)

            # Check if we're ignoring vertical axes such as range or depth.
            # We do this when operating on processed data objects since the
            # assumption with them is that 1d arrays have the dimension equal
            # to the number of pings (horizontal axis).
            if _ignore_vertical_axes and data.shape[0] == my_samples:
                continue

            # Check if the obj_to_insert shares this attribute.
            if hasattr(obj_to_insert, attribute):
                # Get a reference to our obj_to_insert's attribute.
                data_to_insert = getattr(obj_to_insert, attribute)

                # check if the replacement data has the same number of dimensions
                if data.ndim != data_to_insert.ndim:
                    raise ValueError('The replacement data for the' + attribute + ' does not have ' +
                        'the same number of dimensions as the data you are trying to replace. Replacement ' +
                        'data must have the same number of dimensions.')

                # insert the replacement data on top of the existing data.
                if data.ndim == 1:
                    data[replace_index] = data_to_insert[:]
                elif data.ndim == 2:
                    data[replace_index, :] = data_to_insert[:,:]
                elif data.ndim == 3:
                    data[replace_index, :, :] = data_to_insert[:,:,:]
            else:
                raise TypeError('The object providing the replacement data does not have ' +
                        'the ' + attribute + ' attribute. Objects containing replacement data ' +
                        'must at least have the same data attributes as the object whose ' +
                        'data you are replacing.')

        # Now update our global properties.
        if hasattr(obj_to_insert, 'channel_id'):
            if obj_to_insert.channel_id not in self.channel_id:
                self.channel_id += " :: " + obj_to_insert.channel_id

        # Update the size/shape attributes.
        self.n_samples = my_samples
        self.shape = self._shape()


    def delete(self, start_ping=None, end_ping=None, start_time=None,
               end_time=None, remove=True, index_array=None):
        """Deletes data from an echolab2 data object.

        This method deletes data by ping over the range defined by the start
        and end ping times. If remove is True, the data arrays are shrunk.
        If False, the arrays stay the same size and the data values are set
        to NaNs (or an appropriate value based on type).

        Args:
            start_ping (int): The starting ping of the range of pings to delete.
            end_ping (int): The ending ping of the range of pings to delete.
            start_time (datetime64): The starting time of the range of pings to
                delete.
            end_time (datetime64): The ending time of the range of pings to
                delete.

                You should set only one start and end point.

            remove (bool): Set to True to remove the specified pings and
                shrink the data arrays. When set to False, data in the
                deleted pings are set to Nan (or appropriate value for the
                data type).
            index_array (array): A numpy array containing the indices of the
                pings you want to delete. Unlike when using starts/ends to
                define the ping range to delete, the pings do not have to be
                consecutive. When this keyword is present, the start/end
                keywords are ignored.
        """
        # Determine the indices of the pings we're deleting.
        if index_array is None:
            # We haven't been provided an explicit array, so create one based
            # on provided ranges.
            del_idx = self.get_indices(start_time=start_time, end_time=end_time,
                    start_ping=start_ping, end_ping=end_ping)
        else:
            # Explicit array provided.
            del_idx = index_array

        # Determine the indices of the pings we're keeping.
        keep_idx = np.delete(np.arange(self.ping_time.shape[0]), del_idx)

        # Determine the number of pings we're keeping.
        new_n_pings = keep_idx.shape[0]

        # Work through the attributes to delete the data.  If we're removing
        # the pings, we first copy the data we're keeping to a contiguous
        # block before we resize all of the arrays (which will shrink them).
        # If we're not removing the pings, we simply set the values of the
        # various attributes we're deleting to NaNs.
        for attr_name in self._data_attributes:
            attr = getattr(self, attr_name)
            if isinstance(attr, np.ndarray) and attr.ndim > 1:
                if attr.ndim == 2:
                    # 2d array
                    if remove:
                        attr[0:new_n_pings, :] = attr[keep_idx, :]
                    else:
                        attr[del_idx, :] = np.nan
                elif attr.ndim == 3:
                    # 3d array
                    if remove:
                        attr[0:new_n_pings, :, :] = attr[keep_idx, :, :]
                    else:
                        attr[del_idx, :, :] = np.nan
            else:
                if remove:
                    try:
                        # Copy the data we're keeping into a contiguous block.
                        attr[0:new_n_pings] = attr[keep_idx]
                    except:
                        # Range/Depth will fail and that's OK
                        pass
                else:
                    # Set the data to NaN or appropriate value.
                    if np.issubdtype(attr.dtype, np.integer):
                        # Set all other integers to zero
                        attr[del_idx] = 0
                    else:
                        # This is a float like object so set to NaN or NaT
                        attr[del_idx] = np.nan

        # If we're removing the pings, shrink the arrays.
        if remove:
            self.resize(new_n_pings, self.n_samples)
            self.n_pings = new_n_pings


    def append(self, obj_to_append, force=False, time_order=False):
        """Appends another echolab2 data object to this one.

        The objects must be instances of the same class and share the same
        frequency to append one to the other.

        Args:
            time_order (bool): Set to True to force the result of the append
                to be ordered by time. Append implies adding data to the end
                of your data array, but depending on the times within the
                data you are appending, setting time_order=True can result
                in data being inserted throughout "this" object.
        """

        # Simply inserts a data object at the end of our internal array.
        self.insert(obj_to_append, ping_number=self.n_pings, force=force,
                    time_order=time_order)


    def insert(self, obj_to_insert, ping_number=None, ping_time=None,
               insert_after=True, index_array=None, force=False,
               time_order=False):
        """Inserts data from the provided echolab2 data object into
        this object.

        The insertion point is specified by ping number or time (you must
        specify a ping number or ping time). Existing data from the insertion
        point onward will be shifted after the inserted data.  After
        inserting data, the ping_number property is updated and the ping
        numbers will be re-numbered accordingly.

        Args:
            ping_number (int): The ping number specifying the insertion point
            ping_time (datetime): The ping time specifying the insertion point
            insert_after (bool): Set to True to insert *after* the specified
                ping time or ping number. Set to False to insert *at* the
                specified time or ping number.
            index_array (array): A numpy array containing the indices of the
                pings you want to insert. Unlike when using a ping number or
                ping time, the pings do not have to be consecutive. When this
                keyword is present, the ping_number, ping_time and
                insert_after keywords are ignored.
            force (bool): Set to true to disable all checks and force the
                insert even if the frequency, channel ID, etc are different.
            time_order (bool): Set to True to force the result of the append
                to be ordered by time.
                    Default: False

        Raises:
            ValueError: Insertion point not specified.
            TypeError: The object is not an instance of ping_data class.
            TypeError: The frequency of the object to be inserted
                doesn't match the frequency of this object.
            TypeError: Index array isn't a numpy array.
            IndexError: The length of the index array does not match the
                number of pings in the object to be inserted.
        """

        # if time_order is specified, it will override ping_number and ping_time
        # keywords. Set ping_number to n_pings to add data to the end of the arrays
        # before reordering by time.
        if time_order:
            ping_number = self.n_pings

        # Check that we have been given an insertion point or index array.
        if ping_number is None and ping_time is None and index_array is None:
            raise ValueError('Either ping_number or ping_time needs to be ' +
                             'defined or an index array needs to be provided ' +
                             'to specify an insertion point.')

        # Make sure that obj_to_insert class matches "this" class.
        if not isinstance(self, obj_to_insert.__class__):
            raise TypeError('The object you are inserting/appending must ' +
                            'be an instance of ' + str(self.__class__))

        # Make sure that the frequencies match.
        if not force:
            freq_match = self._check_frequencies(obj_to_insert)
            if not freq_match:
                raise ValueError('The object that you are inserting does not have the ' +
                        'same frequency attributes as this object. You cannot insert '+
                        'data with different frequency attributes.')

        # Get some info about the shape of the data we're working with.
        my_pings = self.n_pings
        new_pings = obj_to_insert.n_pings
        my_samples = self.n_samples
        new_samples = obj_to_insert.n_samples

        # Determine the index of the insertion point or indices of the pings
        # we're inserting.
        if index_array is None:
            # Determine the index of the insertion point.
            insert_index = self.get_indices(start_time=ping_time,
                    end_time=ping_time, start_ping=ping_number,
                    end_ping=ping_number)[0]

            # Check if we're inserting before or after the provided insert
            # point and adjust as necessary.
            if insert_after:
                # We're inserting *after*.
                insert_index += 1

            # Create an index array.
            insert_index = np.arange(new_pings) + insert_index

            # Generate the index used to move the existing pings.
            move_index = np.arange(my_pings)
            idx = move_index >= insert_index[0]
            move_index[idx] = move_index[idx] + new_pings

        else:
            # Explicit array provided.  These will be a vector of locations
            # to insert.
            insert_index = index_array

            # Make sure the index array is a numpy array.
            if (not isinstance(insert_index, np.ndarray)):
                raise TypeError('index_array must be a numpy.ndarray.')

            # If we are inserting with a user provided index, make sure the
            # dimensions of the index and the object to insert match.
            if insert_index.shape[0] != new_pings:
                raise IndexError('The length of the index_array does not ' +
                        'match the number of pings in the object you are inserting.')

            # Generate the index used to move the existing pings.
            move_index = np.arange(my_pings)
            for i in insert_index:
                idx = move_index >= i
                move_index[idx] = move_index[idx] + 1

        # Check if we need to vertically resize one of the arrays.  We
        # resize the smaller array to the size of the larger array.  It will
        # automatically be padded with NaNs.
        if my_samples < new_samples:
            # Resize our data arrays and check if we have a limit on the
            # max number of samples.
            if hasattr(self, 'max_sample_number') and self.max_sample_number:
                # We have the attribute and a value is set.  Check if the
                # new array exceeds our max_sample_count.
                if new_samples > self.max_sample_number:
                    # We have to change our new_samples.
                    new_samples = self.max_sample_number
                    # Vertically trim the array we're inserting.
                    obj_to_insert.resize(new_pings, new_samples)
            # Set the new vertical sample size.
            my_samples = new_samples
        elif my_samples > new_samples:
            # Resize the object we're inserting.
            obj_to_insert.resize(new_pings, my_samples)

        # Update the number of pings in the object we're inserting into
        # and then resize it.
        my_pings = my_pings + new_pings
        self.resize(my_pings, my_samples)

        # Generate the move index.
        move_idx = np.arange(move_index.shape[0])

        # Update ping time first, in case we're returning in time order
        self.ping_time[move_index[::-1],] = self.ping_time[move_idx[::-1]]
        # Insert the new data.
        self.ping_time[insert_index] = obj_to_insert.ping_time[:]
        if time_order:
            to_index = self.ping_time.argsort(kind='stable')
            self.ping_time[:] = self.ping_time[to_index]

        # Work through our data properties, inserting the data from obj_to_insert.
        for attribute in self._data_attributes:

            #  skip ping time since we've handled that
            if attribute == 'ping_time':
                continue

            # Check if we have data for this attribute.
            if not hasattr(self, attribute):
                # data_obj does not have this attribute, move along.
                continue

            # Get a reference to our data_obj's attribute.
            data = getattr(self, attribute)

            # Check if the obj_to_insert shares this attribute.
            if hasattr(obj_to_insert, attribute):
                # Get a reference to our obj_to_insert's attribute.
                data_to_insert = getattr(obj_to_insert, attribute)

                # check if the sata we're inserting data has the same number of dimensions
                if data.ndim != data_to_insert.ndim:
                    raise ValueError('The ' + attribute + ' data you are inserting/appending does not have ' +
                        'the same number of dimensions as the data you are inserting into or appending to.')

                # We have to handle the 2d and 1d differently.
                if data.ndim == 1 and data.shape[0] != my_samples:
                    # Skip vertical axis attributes, but move the other 1d data
                    # move right to left to avoid overwriting data before
                    # moving.
                    data[move_index[::-1],] = data[move_idx[::-1]]
                    # Insert the new data.
                    data[insert_index] = data_to_insert[:]
                    if time_order:
                        data[:] = data[to_index]
                elif data.ndim == 2:
                    # Move the existing data from right to left to avoid
                    # overwriting data yet to be moved.
                    data[move_index[::-1], :] = data[move_idx[::-1], :]
                    # Insert the new data.
                    data[insert_index, :] = data_to_insert[:, :]
                    if time_order:
                        data[:,:] = data[to_index,:]
                elif data.ndim == 3:
                    # Move 3d data
                    data[move_index[::-1], :, :] = data[move_idx[::-1], :, :]
                    # Insert the new data.
                    data[insert_index, :, :] = data_to_insert[:, :, :]
                    if time_order:
                        data[:,:,:] = data[to_index,:,:]
            else:
                raise TypeError('The object you are inserting/appending does not have ' +
                        'the ' + attribute + ' attribute. Objects that you insert or ' +
                        'append must at least have the same data attributes as the object you are ' +
                        'inserting into or appending to.')

        # Now update our global properties.
        if hasattr(obj_to_insert, 'channel_id'):
            if obj_to_insert.channel_id not in self.channel_id:
                self.channel_id += " :: " + obj_to_insert.channel_id

        # Update the size/shape attributes.
        self.n_pings = self.ping_time.shape[0]
        self.n_samples = my_samples
        self.shape = self._shape()


    def match_pings(self, other_data, match_to='cs', insert_only=False,
            remove_only=False):
        """Matches the ping times in this object to the ping times in the object
        provided. It does this by inserting and/or deleting pings. Ping times in
        the other object that aren't in this object are inserted. Ping times in
        this object that aren't in the other object are deleted from this object.
        If the time axes do not intersect at all, all of the data in this object
        will be deleted and replaced with empty pings for the ping times in the
        other object.

        Args:
            other_data (ping_data): A ping_data type object that this object
            will be matched to.

            match_to (str): Set to a string defining the precision of the match.
            A lower precision allows matching in cases where there is a slight
            time difference between the two data sources.

                cs : Match to a 100th of a second
                ds : Match to a 10th of a second
                s  : Match to the second

                Default: 'cs'

            insert_only (bool): Set to True to only insert pings into this object
            that aren't in the other object. Pings in this object that aren't in
            the other object will not be removed.

                Default: False

            remove_only (bool): Set to True to only remove pings in this object
            that aren't in the other object. Pings in the other object that aren't
            in this object will not be inserted. This is useful when you want to
            re-write raw data after matching pings.

                Default: False

        Returns:
            A dictionary with the keys 'inserted' and 'removed' containing the
            indices of the pings inserted and removed.

        """
        # Create a dict to store info on which pings were inserted/removed
        results = {'inserted':[], 'removed':[]}

        if match_to == 'cs':
            round_amt = np.uint64(5)
            truncate_to = -2
        elif match_to == 'ds':
            round_amt = np.uint64(50)
            truncate_to = -3
        elif match_to == 's':
            round_amt = np.uint64(500)
            truncate_to = -4

        # don't allow a recursive match
        if other_data is not self:

            # round our times to allow for a loose match window
            this_time = np.around(self.ping_time.astype('uint64') + round_amt,
                    truncate_to)
            other_time = np.around(other_data.ping_time.astype('uint64') + round_amt,
                    truncate_to)

            if not insert_only:
                #  remove any "extra" pings this object may have
                idx_out = np.isin(this_time, other_time, invert=True)
                idx_out = np.nonzero(idx_out)[0]
                if idx_out.size > 0:
                    # We have some extra pings, delete them
                    results['removed'] = idx_out
                    self.delete(index_array=idx_out)

            if not remove_only:
                # Insert any pings that this object is missing
                idx_in = np.isin(other_time, this_time, invert=True)
                idx_in = np.nonzero(idx_in)[0]
                if idx_in.size > 0:
                    # There were missing pings, we'll insert "empty" pings in their place
                    results['inserted'] = idx_in
                    self.insert(self.empty_like(len(idx_in)),
                            index_array=idx_in, force=True)

                    # Lastly, update the times
                    self.ping_time[idx_in] = other_data.ping_time[idx_in]

        return results


    def trim(self, n_pings=None):
        """Trims pings from an echolab2 data object to a given length.

        This method deletes pings from a data object to a length defined by
        n_pings.

        Args:
            n_pings (int): Number of pings (horizontal axis).
        """

        if not n_pings:
            n_pings = self.n_pings

        # Resize keeping the sample number the same.
        self.resize(n_pings, self.n_samples)


    def roll(self, roll_pings):
        """Rolls our data array elements along the ping axis.


        THIS METHOD IS UNTESTED AND PROBABLY DOESN'T WORK

        Elements that roll beyond the last position are re-introduced at the
        first position.

        Args:
            roll_pings ():
        """

        #TODO: Test these inline rolling functions
        #      Need to profile this code to see which methods are faster.
        #      Currently all rolling is implemented using np.roll which makes
        #      a copy of the data.
        #TODO: verify rolling direction
        #      Verify the correct rolling direction for both the np.roll
        #      calls and the 2 inline functions. I *think* the calls to
        #      np.roll are correct and the inline functions roll the wrong way.

        def roll_1d(data, n):
            # Rolls a 1d mostly in place.  Based on code found here:
            #    https://stackoverflow.com/questions/35916201/alternative-to
            #    -numpy-roll-without-copying-array
            # THESE HAVE NOT BEEN TESTED
            temp_view = data[:-n]
            temp_copy = data[-n]
            data[n:] = temp_view
            data[0] = temp_copy

        def roll_2d(data, n):
            # Rolls a 2d mostly in place.
            temp_view = data[:-n,:]
            temp_copy = data[-n,:]
            data[n:, :] = temp_view
            data[0, :] = temp_copy

        # Work through our list of attributes.
        for attr_name in self._data_attributes:

            # Get a reference to this attribute.
            attr = getattr(self, attr_name)

            # Resize the arrays using a technique dependent on the array
            # dimension.
            if attr.ndim == 1:
                attr = np.roll(attr, roll_pings)
                # attr[:] = roll_1d(attr, roll_pings)
            elif attr.ndim == 2:
                attr = np.roll(attr, roll_pings, axis=0)
                # attr[:] = roll_2d(attr, roll_pings)

            # Update the attribute.
            setattr(self, attr_name, attr)


    def resize(self, new_ping_dim, new_sample_dim):
        """Iterates through the provided list of attributes and resizes them.

        The size of the attributes in the instance of the provided object
        is resized given the new array dimensions.

        Args:
            new_ping_dim (int): Ping dimension gives the width of the array (
                horizontal axis).
            new_sample_dim (int): Sample dimension gives the height of the
                array (vertical axis).
        """

        def _resize2d(data, ping_dim, sample_dim, fill_value=np.nan):
            """
            _resize2d returns a new array of the specified dimensions with the
            data from the provided array copied into it. This function is
            used when we need to resize 2d arrays along the minor axis as
            ndarray.resize and numpy.resize don't maintain the order of the
            data in these cases.
            """
            # Create the new array.
            new_array = np.empty((ping_dim, sample_dim), dtype=self.sample_dtype)

            #  determine the copy bounds and fill the edges if needed
            n_pings = np.min((data.shape[0], ping_dim))
            if n_pings == data.shape[0]:
                new_array[n_pings:,:] = fill_value
            n_samps = np.min((data.shape[1], sample_dim))
            if n_samps == data.shape[1]:
                new_array[:, n_samps:] = fill_value
            new_array[n_pings:,:] = fill_value

            # Copy the data into our new array and return it.
            new_array[0:n_pings, 0:n_samps] = data[0:n_pings, 0:n_samps]
            return new_array


        def _resize3d(data, ping_dim, sample_dim, sector_dim, fill_value=np.nan):
            """
            _resize3d returns a new array of the specified dimensions with the
            data from the provided array copied into it. Same reasoning as above.
            """
            # Create the new array.
            new_array = np.empty((ping_dim, sample_dim, sector_dim), dtype=self.sample_dtype)

            #  determine the copy bounds and fill the edges if needed
            n_pings = np.min((data.shape[0], ping_dim))
            if n_pings == data.shape[0]:
                new_array[n_pings:,:,:] = fill_value
            n_samps = np.min((data.shape[1], sample_dim))
            if n_samps == data.shape[1]:
                new_array[:, n_samps:] = fill_value
            new_array[n_pings:,:,:] = fill_value

            # Copy the data into our new array and return it.
            new_array[0:n_pings, 0:n_samps, :] = data[0:n_pings, 0:n_samps, :]
            return new_array

        # Store the old sizes.
        old_sample_dim = self.n_samples
        old_ping_dim = self.ping_time.shape[0]

        # Ensure our values are integers.  Some platforms/versions don't
        # automatically coerce floats to integers when used as integer
        # arguments.
        new_ping_dim = int(new_ping_dim)
        new_sample_dim = int(new_sample_dim)

        # Work through our list of attributes.
        for attr_name in self._data_attributes:

            # Get a reference to this attribute.
            attr = getattr(self, attr_name)

            if isinstance(attr, np.ndarray):

                # Resize the arrays using a technique dependent on the array dimension.
                if attr.ndim == 1:
                    # 1d arrays can be an axis or be on the ping axes or sample axes and have
                    # to be handled differently.
                    if attr.shape[0] == old_sample_dim != new_sample_dim:
                        # Resize this sample axes attribute.
                        attr = np.resize(attr,(new_sample_dim))
                    elif attr.shape[0] == old_ping_dim != new_ping_dim:
                        # Resize this ping axes attribute.
                        attr = np.resize(attr,(new_ping_dim))
                elif attr.ndim == 2:
                    # Resize this 2d sample data array.
                    if new_sample_dim == old_sample_dim:
                        # If the minor axes isn't changing, we can use
                        # np.resize() function.
                        attr = np.resize(attr,(new_ping_dim, new_sample_dim))
                    else:
                        # If the minor axes is changing, we need to use our
                        # resize2d function.
                        attr = _resize2d(attr, new_ping_dim, new_sample_dim)
                elif attr.ndim == 3:
                    # Resize this 3d sample data array.
                    if new_sample_dim == old_sample_dim:
                        # If the minor axes isn't changing, we can use
                        # np.resize() function.
                        attr = np.resize(attr,(new_ping_dim, new_sample_dim, self.n_complex))
                    else:
                        # If the minor axes is changing, we need to use our
                        # resize2d function.
                        attr = _resize3d(attr, new_ping_dim, new_sample_dim, self.n_complex)

                #  Update the attribute.
                setattr(self, attr_name, attr)

        # Set the new shape attributes
        self.n_samples = new_sample_dim
        self.shape = self._shape()


    def get_indices(self, start_ping=None, end_ping=None, start_time=None,
                    end_time=None, time_order=True, **_):
        """Returns a index array containing where the indices in the
        range defined by the times and/or ping numbers provided are True.

        By default, the indices are returned in time order. If time_order is set
        to False, the data will be returned in the order they occur in the data
        arrays.

        Note that pings with "empty" times (ping time == NaT) will be sorted
        to the beginning of the index array for numpy versions < 1.18 and the
        END for versions >= 1.18

        Args:
            start_ping (int): The starting ping of the range of pings specified.
            end_ping (int): The ending ping of the range of pings specified.
            start_time (datetime64): The starting time of the range of pings
                specified.
            end_time (datetime64): The ending time of the range of pings
                specified.
            time_order (bool): Controls the order the indices will return.  If
                set to True, the indices will be in time order.  If False,
                the data will return in the order they occur in the data arrays.

        Returns:
            The indices that are included in the specified range.
        """

        # Generate the ping number vector.  We start counting pings at 1.
        ping_number = np.arange(self.n_pings) + 1

        # If starts and/or ends are omitted, assume first and last respectively.
        if start_ping is None and start_time is None:
            start_ping = ping_number[0]
        if end_ping is None and end_time is None:
            end_ping = ping_number[-1]

        # Get the primary index.
        if time_order:
            # Return indices in time order.  Note that empty ping times will be
            # sorted to the front for numpy versions < 1.18 and at the end for
            # versions >= 1.18
            primary_index = self.ping_time.argsort(kind='stable')
        else:
            # Return indices in ping order.
            primary_index = ping_number - 1

        # Generate a boolean mask of the values to return.
        if start_time is not None:
            mask = self.ping_time[primary_index] >= start_time
        elif start_ping is not None:
            mask = ping_number[primary_index] >= start_ping
        if end_time is not None:
            mask = np.logical_and(mask, self.ping_time[primary_index] <= end_time)
        elif end_ping is not None:
            mask = np.logical_and(mask, ping_number[primary_index] <= end_ping)

        # Return the indices that are included in the specified range.
        return primary_index[mask]


    def resample_by_axes(self, v_height, h_length, v_axis='range', h_axis='ping_number',
            method='area_weighted', skip_attributes=[]):
        """resample_by_axes resamples pings in both vertical and/or along track dimensions.
       
        This is similar to the Echoview "Resample by..." operators except that you have finer
        control over the vertical axes when resampling.
       
        Args:
            v_height (int, float, None): Specify the vertical height of the new samples in 
            the specified axis units. If depth or range are specified, the units are meters.
            If sample is specified, the units are sample numbers. If set to None, the
            vertical axis is not resampled.

            h_length (int, float, timestamp64, None): Specify the horizontal extent of the 
            new samples in the specified axis units. If ping_number is specified, the units
            are pings. If ping_time is specified, the units are time in datetime64 format.
            If h_axis is specified as trip_distance_m, the units are meters and for
            trip_distance_nmi the units are nautical miles. If set to None, the horizontal
            axis is not resampled.

            v_axis (str): Set to a string specifying the vertical axis to use when
            resampling.
            
                'range':  
                'depth':
                'sample':

                Default: range

            h_axis (str): Set to a string specifying the horizontal axis to use when
            resampling.
            
                'ping_number':  
                'ping_time':
                'trip_distance_m':
                'trip_distance_nmi'

                Default: ping_number

            method (str): Set to a string specifying the method to use when resampling. Options
            are:
                'area_weighted': Samples that are only partially covered by the new grid will be
                weighted according to the fraction of the sample that is covered. This is also
                known as conservative regridding and is the default method.

                'center_weighted': Samples that are only partially covered by the new grid will be
                weighted as 1 if the center of the sample is covered and 0 if it is not.

                Default: 'area_weighted'

            skip_attributes (list): Pass a list of data attributes specifying attributes that
            should be skipped when resampling data stored within this object. Normally
            this method will resample all of the objects data attributes (like lat/lon, SOG,
            bottom depth, etc.) so they are the same size and remain aligned with the object's
            new axes. This may not be suitable for all data and you can specify attributes to
            skip here.

        """

        warnings.warn("The resample_by_axes() method has not been fully tested and vetted. Please "
                "check all results for sanity.")

        # check if we have been given any new grid dimensions. If not, just return.
        if h_length is None and v_height is None:
            return

        # Transform ping time to a float. For datetime64[ms] objects this is the number of ms since the epoch.
        ping_time = self.ping_time.astype('float64')

        # Only regrid alongtrack if we have been given a new horizontal/alongtrack length
        if h_length is not None:

            # Get the horizontal axis data - perform any transformations as needed
            if h_axis == 'ping_number':
                # ping number is not an innate attribute and is generated
                h_axis_data = np.arange(self.n_pings, dtype='float32') + 1

            elif h_axis == 'ping_time':
                # Ping time needs to be converted to a float for gridding and we can use the already converted values
                h_axis_data = ping_time

                # Transform the horizontal length to a float. Since the interval length can be specified
                # in arbitrary time units we must first get it in ms then get that as a float64
                h_length = h_length.astype('datetime64[ms]').astype('float64')

            elif h_axis == 'trip_distance_m':
                # Trip distance in m a computed axis - we take vessel log (trip_distance_nmi) and convert to meters
                if  hasattr(self, 'trip_distance_nmi'):
                    # get a copy of the interval axis data
                    h_axis_data = getattr(self, h_axis).copy().astype('float32')
                    h_axis_data *= 1852
                else:
                    raise AttributeError("This object lacks the horizontal axis attribute 'trip_distance_nmi' " +
                            "which is required for the specified 'trip_distance_m' horizontal axis.")
            else:
                # This axis type doesn't need any special treatment
                if  hasattr(self, h_axis):
                    # get a copy of the interval axis data
                    h_axis_data = getattr(self, h_axis).astype('float32')
                else:
                    raise AttributeError("This object lacks the specified " +
                            "horizontal axis attribute '" + h_axis + "'.")

            # Generate the new horizontal grid attributes
            n_intervals, interval_edges, interval_centers, interval_pings, \
                    ping_interval_map = self._grid_axis(h_axis_data, h_length)

            # Compute new ping times for our new horizontal grid
            new_ping_times = np.empty((n_intervals), dtype='float64')
            for i in range(n_intervals):
                new_ping_times[i] = np.nansum(ping_time[np.ix_(ping_interval_map==i)]) / interval_pings[i]
            # Convert our new times to datetime64[ms]...
            new_ping_times = new_ping_times.astype('datetime64[ms]')
            # and update
            setattr(self, 'ping_time', new_ping_times)

        if v_height is not None:

            # Get the vertical axis data
            if v_axis.lower() == 'sample':
                # like ping number, sample number is not an innate attribute and is generated
                v_axis_data = np.arange(self.n_samples, dtype='float32') + 1
            else:
                if  hasattr(self, v_axis):
                    v_axis_data = getattr(self, v_axis)
                else:
                    raise AttributeError("This object lacks the specified " +
                            "vertical axis attribute '" + v_axis + "'.")

            # Generate the new vertical axis grid attributes
            n_layers, layer_edges, new_v_axis, layer_samples, \
                    sample_layer_map = self._grid_axis(v_axis_data, v_height)

            # update the vertical axis - if the the resampling axis is sample based,
            # we need to innterpolate the new range/depth axes values
            if v_axis.lower() == 'sample':
                if hasattr(self, 'range'):
                    new_axis = np.interp(new_v_axis, v_axis_data, self.range)
                    setattr(self, 'range', new_axis)
                if hasattr(self, 'depth'):
                    new_axis = np.interp(new_v_axis, v_axis_data, self.depth)
                    setattr(self, 'depth', new_axis)
            else:
                setattr(self, v_axis, new_v_axis)

        #  create the NaN mask of the data
        mask = np.isnan(self.data)

        # zero out NaNs so they don't contribute to the weighted average.
        self.data[mask] = 0

        # invert the mask to represent the valid data. It will be transformed, along
        # with the data, and then we'll divide it into the data to compute the mean.
        np.logical_not(mask, out=mask)

        # check if we have a horizontal axis to resample
        if h_length is not None:
            # area weighted resampling requires both the source and desination grids to
            # be defined the edges of the samples but our source grid is defined by the
            # centers. Since the horizontal axis is not regularly spaced, we just add
            # a final edge to the end of the horizontal axis data that is the same distance
            # from the last point as the last two points are from each other. 
            h_axis_thickness = h_axis_data[-1] - h_axis_data[-2]
            h_axis_edges = np.append(h_axis_data, h_axis_data[-1] + h_axis_thickness)

            # NOTE! These are for debugging - remove when finished
            self.resample_dest_interval_edges = interval_edges
            self.resample_source_interval_edges = h_axis_edges

            # compute the weights for the horizontal axis resampling.
            if method == 'area_weighted':
                h_weights = self._get_area_weighted_weights(h_axis_edges, interval_edges)
            elif method == 'center_weighted':
                h_weights = self._get_center_weighted_weights(h_axis_edges, interval_edges)
            else:
                raise ValueError("Invalid resampling method specified. Valid options are"
                        " 'area_weighted' and 'center_weighted'.")

            # resample the data along the horizontal axis. Handle the data and sample
            # masks separately.
            resampled_num = h_weights @ self.data
            resampled_den = h_weights @ mask.astype('float')
        else:
            # since we're not resampling the horizontal axis, just pass the data
            # through for the vertical resampling step.
            resampled_num = self.data
            resampled_den = mask.astype('float')

        #  check if we're resampling the vertical axis.
        if v_height is not None:
            # Like the horizontal axis, we need to define the edges of the source 
            # vertical axis. Since the vertical axis is regularly spaced, we can
            # just subtract and add half of the sample thickness to get the edges.
            v_axis_edges = np.append(v_axis_data - self.sample_thickness / 2, 
                    v_axis_data[-1] + self.sample_thickness / 2)
            
            # NOTE! These are for debugging - remove when finished
            self.resample_dest_layer_edges = layer_edges
            self.resample_source_layer_edges = v_axis_edges

            # compute the weights for the vertical axis resampling.
            if method == 'area_weighted':
                v_weights = self._get_area_weighted_weights(v_axis_edges, layer_edges)
            elif method == 'center_weighted':
                v_weights = self._get_center_weighted_weights(v_axis_edges, layer_edges)
            else:
                raise ValueError("Invalid resampling method specified. Valid options are"
                        " 'area_weighted' and 'center_weighted'.")

            # resample the data along the vertical axis. Again, we handle the data
            # and sample nan mask separately.
            resampled_num = resampled_num @ v_weights.T
            resampled_den = resampled_den @ v_weights.T

        # Compute the means by dividing the resampled data values by
        # the resampled "good samples" mask.
        resampled_data = np.divide(resampled_num, resampled_den, 
                        out=np.full_like(resampled_num, np.nan), 
                        where=resampled_den > 0)

        #  update our data attribute with our newly resampled data
        setattr(self, 'data', resampled_data)

        #  update data attribues
        self.n_samples = n_layers
        self.n_pings = n_intervals
        old_shape = self.shape
        self.shape = self.data.shape
        self.sample_thickness = layer_edges[1] - layer_edges[0]

        # The final step is to average the other data attributes so they remain aligned
        # with the new data axes. Here in the base class, we only handle numpy arrays.
        # Attributes of other types (e.g. pyEcholab line object) must be handled in 
        # the appropiate child classes.

        # Average the other data attributes - skip attributes we have already handled
        # or been explicitly told to ignore.
        skip_attributes.extend(['ping_time', 'data', 'range', 'depth'])
        for attr_name in self._data_attributes:
            if attr_name in skip_attributes:
                # This attribute has either been handled already or we have been
                # explicitly told to skip it. Move on. 
                continue

            # Get the attribute data
            attr_data = getattr(self, attr_name)
            
            # Handle attributes based on type. We can average numpy arrays directly if their
            # dimensions match the horizontal or vertical axis.  If they don't match, we will
            # emit a warning and just copy them as-is.

            if isinstance(attr_data, np.ndarray):
                #  we can average numpy data directly. 
                if attr_data.shape[0] == old_shape[0] and h_length is not None:
                    #  this is a horizontal (ping based) data attribute
                    new_attr_data = np.empty((n_intervals), dtype=attr_data.dtype)
                    for i in range(n_intervals):
                        this_data = attr_data[np.ix_(ping_interval_map==i)]
                        good_vals = interval_pings[i] - np.count_nonzero(np.isnan(this_data))
                        if good_vals > 0:
                            #  compute the mean of the non-NaN values in this interval
                            new_attr_data[i] = np.nansum(this_data) / good_vals
                        else:
                            #  all elements are NaN so the result is NaN
                            new_attr_data[i] = np.nan

                    #  update the attribute
                    setattr(self, attr_name, new_attr_data)

                elif attr_data.shape[0] == old_shape[1] and v_height is not None:
                    #  this is a vertical (sample based) data attribute
                    new_attr_data = np.empty((n_layers), dtype=attr_data.dtype)
                    for i in range(n_layers):
                        this_data = attr_data[np.ix_(sample_layer_map==i)]
                        good_vals = layer_samples[i] - np.count_nonzero(np.isnan(this_data))
                        if good_vals > 0:
                            #  compute the mean of the non-NaN values in this layer
                            new_attr_data[i] = np.nansum(this_data) / good_vals
                        else:
                            #  all elements are NaN so the result is NaN
                            new_attr_data[i] = np.nan

                    #  update the attribute
                    setattr(self, attr_name, new_attr_data)

                else:
                    # the size of this attribute doesn't match either axes so we can't resample it.
                    warnings.warn("The dimensions of data attribute " + attr_name + " do not match " +
                            "either axis. This attribute will be ignored.")


    def _get_area_weighted_weights(self, src_b, dst_b):
        """
        src_b: (N+1,) array of source boundaries
        dst_b: (M+1,) array of destination boundaries
        Returns: (M, N) sparse weight matrix
        """

        n_src, n_dst = len(src_b) - 1, len(dst_b) - 1
        src_l, src_h = src_b[:-1], src_b[1:]
        dst_l, dst_h = dst_b[:-1], dst_b[1:]

        weights = np.zeros((n_dst, n_src))
        for i in range(n_dst):
            # Calculate overlap of dst cell 'i' with all src cells
            overlap = np.maximum(0, np.minimum(src_h, dst_h[i]) - np.maximum(src_l, dst_l[i]))
            # Normalize by destination cell width
            weights[i, :] = overlap / (dst_h[i] - dst_l[i])
        
        return csr_matrix(weights)
    

    def _get_center_weighted_weights(self, src_b, dst_b):
        """
        src_b: (N+1,) array of source boundaries
        dst_b: (M+1,) array of destination boundaries
        Returns: (M, N) sparse weight matrix
        """

        # Calculate centers of source cells
        src_centers = (src_b[:-1] + src_b[1:]) / 2.0
        
        # Extract destination boundaries
        dst_l = dst_b[:-1, np.newaxis]
        dst_h = dst_b[1:, np.newaxis]
        
        # Compare all centers against all boundaries
        mask = (src_centers >= dst_l) & (src_centers < dst_h)
        
        # 4. Convert to float and then to sparse matrix
        return csr_matrix(mask.astype(float))


    def brute_force_resample(self, v_height, h_length, v_axis='range', h_axis='ping_number'):
        """
        This method loops thru the intervals/cells to compute the mean using a brute force
        method. It only computes center weighted averages.

        brute_force_resample only exists as a check against the resampling methods implemented
        in resample_by_axes.
        
        This method will be removed when testing is complete.
        """

        # check if we have been given any new grid dimensions. If not, just return.
        if h_length is None or v_height is None:
            raise AttributeError("Brute force resampling requires both horizontal and vertical grid dimensions.")

        # Transform ping time to a float. For datetime64[ms] objects this is the number of ms since the epoch.
        ping_time = self.ping_time.astype('float64')

        # Get the horizontal axis data - perform any transformations as needed
        if h_axis == 'ping_number':
            # ping number is not an innate attribute and is generated
            h_axis_data = np.arange(self.n_pings, dtype='float32') + 1

        elif h_axis == 'ping_time':
            # Ping time needs to be converted to a float for gridding and we can use the already converted values
            h_axis_data = ping_time

            # Transform the horizontal length to a float. Since the interval length can be specified
            # in arbitrary time units we must first get it in ms then get that as a float64
            h_length = h_length.astype('datetime64[ms]').astype('float64')

        elif h_axis == 'trip_distance_m':
            # Trip distance in m a computed axis - we take vessel log (trip_distance_nmi) and convert to meters
            if  hasattr(self, 'trip_distance_nmi'):
                # get a copy of the interval axis data
                h_axis_data = getattr(self, h_axis).copy().astype('float32')
                h_axis_data *= 1852
            else:
                raise AttributeError("This object lacks the horizontal axis attribute 'trip_distance_nmi' " +
                        "which is required for the specified 'trip_distance_m' horizontal axis.")
        else:
            # This axis type doesn't need any special treatment
            if  hasattr(self, h_axis):
                # get a copy of the interval axis data
                h_axis_data = getattr(self, h_axis).astype('float32')
            else:
                raise AttributeError("This object lacks the specified " +
                        "horizontal axis attribute '" + h_axis + "'.")

        # Generate the new horizontal grid attributes
        n_intervals, interval_edges, interval_centers, interval_pings, \
                ping_interval_map = self._grid_axis(h_axis_data, h_length)

        # Compute new ping times for our new horizontal grid
        new_ping_times = np.empty((n_intervals), dtype='float64')
        for i in range(n_intervals):
            new_ping_times[i] = np.nansum(ping_time[np.ix_(ping_interval_map==i)]) / interval_pings[i]
        # Convert our new times to datetime64[ms]...
        new_ping_times = new_ping_times.astype('datetime64[ms]')
        # and update
        setattr(self, 'ping_time', new_ping_times)

        # Get the vertical axis data
        if v_axis.lower() == 'sample':
            # like ping number, sample number is not an innate attribute and is generated
            v_axis_data = np.arange(self.n_samples, dtype='float32') + 1
        else:
            if  hasattr(self, v_axis):
                v_axis_data = getattr(self, v_axis)
            else:
                raise AttributeError("This object lacks the specified " +
                        "vertical axis attribute '" + v_axis + "'.")

        # Generate the new vertical axis grid attributes
        n_layers, layer_edges, new_v_axis, layer_samples, \
                sample_layer_map = self._grid_axis(v_axis_data, v_height)

        # update the vertical axis - if the the axis is sample based, we need
        # to innterpolate the new range/depth axes values
        if v_axis.lower() == 'sample':
            if hasattr(self, 'range'):
                new_axis = np.interp(new_v_axis, v_axis_data, self.range)
                setattr(self, 'range', new_axis)
            if hasattr(self, 'depth'):
                new_axis = np.interp(new_v_axis, v_axis_data, self.depth)
                setattr(self, 'depth', new_axis)
        else:
            setattr(self, v_axis, new_v_axis)

        #  finally resample the sample data
        resampled_data = self._brute_force(n_intervals, n_layers, interval_pings, ping_interval_map,
            layer_samples, sample_layer_map)
        setattr(self, 'data', resampled_data)

        #  update data attribues
        self.n_samples = n_layers
        self.n_pings = n_intervals
        old_shape = self.shape
        self.shape = self.data.shape
        self.sample_thickness = layer_edges[1] - layer_edges[0]

        # Average the other data attributes - skip attributes we handle spe
        skip_attrs = ['ping_time', 'data', 'range', 'depth']
        for attr_name in self._data_attributes:
            if attr_name in skip_attrs:
                # This attribute has been handled already - move on
                continue

            # Get the attribute data
            attr_data = getattr(self, attr_name)
            
            # Handle attributes based on type. We can average numpy arrays directly if their
            # dimensions match the horizontal or vertical axis.  If they don't match, we will
            # emit a warning and just copy them as-is. We can handle specific types as needed.
            # Most other types will be ignored. 

            if isinstance(attr_data, np.ndarray):
                #  we can average numpy data directly. 
                if attr_data.shape[0] == old_shape[0] and h_length is not None:
                    #  this is a horizontal (ping based) data attribute
                    new_attr_data = np.empty((n_intervals), dtype=attr_data.dtype)
                    for i in range(n_intervals):
                        this_data = attr_data[np.ix_(ping_interval_map==i)]
                        good_vals = interval_pings[i] - np.count_nonzero(np.isnan(this_data))
                        if good_vals > 0:
                            #  compute the mean of the non-NaN values in this interval
                            new_attr_data[i] = np.nansum(this_data) / good_vals
                        else:
                            #  all elements are NaN so the result is NaN
                            new_attr_data[i] = np.nan

                    #  update the attribute
                    setattr(self, attr_name, new_attr_data)

                elif attr_data.shape[0] == old_shape[1] and v_height is not None:
                    #  this is a vertical (sample based) data attribute
                    new_attr_data = np.empty((n_layers), dtype=attr_data.dtype)
                    for i in range(n_layers):
                        this_data = attr_data[np.ix_(sample_layer_map==i)]
                        good_vals = layer_samples[i] - np.count_nonzero(np.isnan(this_data))
                        if good_vals > 0:
                            #  compute the mean of the non-NaN values in this layer
                            new_attr_data[i] = np.nansum(this_data) / good_vals
                        else:
                            #  all elements are NaN so the result is NaN
                            new_attr_data[i] = np.nan

                    #  update the attribute
                    setattr(self, attr_name, new_attr_data)

                else:
                    warnings.warn("The dimensions of data attribute " + attr_name + " do not match " +
                            "either axis. This attribute will be ignored.")


    def _brute_force(self, n_intervals, n_layers, interval_pings, ping_interval_map,
            layer_samples, sample_layer_map):

        resamp_data = np.empty((n_intervals,n_layers), dtype=self.data.dtype)

        for i in range(n_intervals):
            ping_idx = ping_interval_map==i
            for j in range(n_layers):
                cell_data = self.data[np.ix_(ping_idx,sample_layer_map==j)]
                n_nans_in_cell = np.count_nonzero(np.isnan(cell_data))
                good_samps = (interval_pings[i] * layer_samples[j]) - n_nans_in_cell
                if good_samps > 0:
                    resamp_data[i,j] = np.nansum(cell_data) / good_samps
                else:
                    resamp_data[i,j] = np.nan

        return resamp_data


    def _check_frequencies(self, other_obj):
        '''_check_frequencies compares the frequency attributes of this object against
        the other_obj to see if they match. We require frequencies to match when performing
        certain operations like insert/append.

        Since the ping_data class is inherited by both the raw_data and processed_data
        classes, we have to handle checking both child clases. 

        raw_data objects for reduced or complex-CW store frequency by ping as a vector
            n-pings long and complex-FM objects store frequency as frequency_start and
            frequency_end vectors n-pings long.
        processed_data objects store frequency as a single value for non Sv(f) and
            as a vector n-frequencies long for Sv(f)/TS(f)

        
        '''
        # First, check if either object is complex-FM. raw_data objects containing
        # FM data will not have the frequency attribute as they have frequency_start
        # and frequency_end attributes instead.
        other_is_raw_fm = not hasattr(other_obj, 'frequency')
        self_is_raw_fm = not hasattr(self, 'frequency')

        # Next, regardless of the type, make sure they are the same
        if other_is_raw_fm != self_is_raw_fm:
            # they are not
            return False

        if other_is_raw_fm:
            # The frequency_start and frequency_end attributes of raw_data objects 
            # containing FM data will always be arrays. First check the frequency_start
                other_obj_freq = other_obj.frequency_start[~np.isnan(other_obj.frequency_start)]

                # If there is at least 1 frequency that isn't NaN, we check it
                if other_obj_freq.shape[0] > 0:

                    # check that the first element of other matches all in self, and that 
                    # the first element in other matches all in other. If both are true, then
                    # our frequencies match.
                    other_matches_self = np.all(self.frequency_start[~np.isnan(self.frequency_start)] == other_obj_freq[0])
                    other_matches_other = np.all(other_obj_freq == other_obj_freq[0])

                    if not (other_matches_self and other_matches_other):
                        # The frequency starts do not match
                        return False
                    
                # Now check the frequency end
                other_obj_freq = other_obj.frequency_end[~np.isnan(other_obj.frequency_end)]

                # If there is at least 1 frequency that isn't NaN, we check it
                if other_obj_freq.shape[0] > 0:

                    # check that the first element of other matches all in self, and that 
                    # the first element in other matches all in other. If both are true, then
                    # our frequencies match.
                    other_matches_self = np.all(self.frequency_end[~np.isnan(self.frequency_end)] == other_obj_freq[0])
                    other_matches_other = np.all(other_obj_freq == other_obj_freq[0])

                    if not (other_matches_self and other_matches_other):
                        # The frequency ends do not match
                        return False

        else:
            # Neither object contains complex-FM data, so we perform all checks on the
            # frequency attribute

            # Check if the other object's frequency attribute is an array
            if isinstance(other_obj.frequency, np.ndarray):
                # It is, which means the other object is a raw_data object or a processed_data 
                # object containing Sv(f) or TS(f)

                # Check that this object's frequency attribute is also an array
                if not isinstance(self.frequency, np.ndarray):
                    # it is not so this is not the driod we're looking for
                    return False

                #  get the other object's frequency info - ignoring NaNs
                other_obj_freqs = other_obj.frequency[~np.isnan(other_obj.frequency)]

                # If there is at least 1 frequency that isn't NaN, we check it
                if other_obj_freqs.shape[0] > 0:

                    # check that the first element of other matches all in self, and that 
                    # the first element in other matches all in other. If both are true, then
                    # our frequencies match.
                    other_matches_self = np.all(self.frequency[~np.isnan(self.frequency)] == other_obj_freqs[0])
                    other_matches_other = np.all(other_obj_freqs == other_obj_freqs[0])

                    if not (other_matches_self and other_matches_other):
                        # No match
                        return False

                else:
                    # The other object must have no data (frequency values are all NaN) so we
                    # match because you can insert/append "empty" data under the assumption that
                    # you know what you are doing and will populate it accordingly.
                    return True

            else:
                # the other object's frequency attribute is a scalar - this means that it is
                # a processed_data object containing a single frequency.

                # Check that this object's frequency attribute is also a scalar
                if isinstance(self.frequency, np.ndarray):
                    # it is not so this is not the driod we're looking for
                    return False
                
                # If we're here, we have two objects with scalar frequency attributes
                if self.frequency != other_obj.frequency:
                    return False

        # If we're here, we're good
        return True


    def _grid_axis(self, axis_data, axis_size):

        '''
        _grid_axis is an internal method that generates the grid parameters for the
        provided axis_data and size (horizontal length or vertical height)
        '''

        # Get the span of the axis data
        span = float(axis_data[-1]) - float(axis_data[0])

        # Compute the number of intervals/cells
        n_units = np.ceil(span / axis_size).astype('uint32')

        # compute the interval/cell edges, including the rightmost/bottommost edge
        axis_edges = (np.arange(n_units + 1) * float(axis_size)) + axis_data[0]
        #axis_edges[-1] = axis_data[-1]

        # create the axis mapping array - we include all pings/samples in an interval/cell
        # that are >= to the interval/cell start and < the interval/cell end EXCEPT FOR
        # THE LAST INTERVAL where we include the last ping/sample if it <= the interval/cell
        # end. intervals/samples that do not map to the grid are assigned a map value of -1.
        axis_map = np.full(axis_data.shape, -1, dtype='int32')
        n_els = np.full((n_units), 0, dtype='uint32')
        axis_centers = np.empty((n_units), dtype='float64')
        for b in range(n_units):
            if b < (n_units - 1):
                # for all intervals up to the last interval/layer - include >= start and < end
                mask = np.logical_and(axis_data >= axis_edges[b], axis_data < axis_edges[b+1])
            else:
                # for the last interval/layer, include >= start and <= end to ensure last ping/sample is captured
                mask = np.logical_and(axis_data >= axis_edges[b], axis_data <= axis_edges[b+1])

            #  store the number of elements, the map, and compute the center
            n_els[b] = mask.sum()
            axis_map[mask] = b
            axis_centers[b] = np.nansum(axis_data[np.ix_(axis_map==b)]) / n_els[b]

        return n_units, axis_edges, axis_centers, n_els, axis_map


    def _vertical_resample(self, data, sample_intervals,
                           unique_sample_intervals, resample_interval,
                           sample_offsets, min_sample_offset, is_power=True):
        """Internal method that vertically resamples sample data given a target
        sample interval.

        If the resampling factor is a whole number, samples will be replicated
        when expanding vertically and averaged when reduced. If the resampling
        factor is a fractional number, the samples will be interpolated.
        Interpolation will result in aliasing. Aliasing while upsampling is
        relatively minor, but it can be significant when downsampling. Care
        should be taken when downsampling. The default behavior is to upsample.

        This method also shifts samples vertically based on their sample
        offset so they are positioned correctly relative to each other. The
        first sample in the resulting array will have an offset that is the
        minimum of all offsets in the data.

        Args:
            data:
            sample_intervals:
            unique_sample_intervals:
            resample_interval:
            sample_offsets:
            min_sample_offset:
            is_power:

        Returns:
            The resampled data and the sampling interval used.
        """

        def isclose(a, b, rel_tol=1e-9, abs_tol=0.0):
            return abs(a - b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

        # Determine the number of pings in the new array.
        n_pings = data.shape[0]

        # Check if we need to substitute our resample_interval value.
        if resample_interval == 0:
            # Resample to the shortest sample interval in our data.
            resample_interval = min(unique_sample_intervals)
        elif resample_interval >= 1:
            # Resample to the longest sample interval in our data.
            resample_interval = max(unique_sample_intervals)

        # Generate a vector of sample counts.  The generalized method works
        # with both raw_data and processed_data classes and finds the first
        # non-NaN value searching from the "bottom up".
        if data.ndim == 3:
            # just use the 1st element for complex data types
            sample_counts = data.shape[1] - np.argmax(~np.isnan(np.fliplr(data[:,:,0])), axis=1)
        else:
            sample_counts = data.shape[1] - np.argmax(~np.isnan(np.fliplr(data)), axis=1)

        # Create a couple of dictionaries to store resampling parameters by
        # sample interval.  They will be used when we fill the output array
        # with the resampled data.
        resample_factor = {}
        pings_this_interval = {}
        sample_offsets_this_interval = {}

        # Determine the number of samples in the output array.  To do this,
        # we must loop through the sample intervals, determine the resampling
        # factor, then find the maximum sample count at that sample interval
        # (taking into account the sample's offset) and multiply by the
        # resampling factor to determine the max number of samples for that
        # sample interval.
        new_sample_dims = 0
        for sample_interval in unique_sample_intervals:
            # Set the resampling factor for pings with this sample interval
            resample_factor[sample_interval] = sample_interval / resample_interval

            # Determine the rows in this subset with this sample interval.
            pings_this_interval[sample_interval] = np.where(sample_intervals == sample_interval)[0]

            # Determine the net vertical shift for the samples with this sample interval.
            sample_offsets_this_interval[sample_interval] = \
                    sample_offsets[pings_this_interval[sample_interval]] - min_sample_offset

            # Also, determine the maximum number of samples for this sample interval.  This
            # has to be done on a row-by-row basis since sample number can change between
            # pings. We include the sample offset to ensure we have room to shift our
            # samples vertically by the offset.
            max_samples_this_sample_int = max(sample_counts[pings_this_interval[sample_interval]] +
                sample_offsets_this_interval[sample_interval])

            # Now compute the number of samples for this sample interval and
            # store the largest value over all of the intervals.
            max_dim_this_sample_int = int(round(max_samples_this_sample_int *
                    resample_factor[sample_interval]))
            if max_dim_this_sample_int > new_sample_dims:
                new_sample_dims = max_dim_this_sample_int

        # Now that we know the dimensions of the output array, create it and fill with NaNs.
        if data.ndim == 3:
            resampled_data = np.empty((n_pings, new_sample_dims, data.shape[2]), order='C')
        else:
            resampled_data = np.empty((n_pings, new_sample_dims),order='C')
        resampled_data.fill(np.nan)

        # Now fill the array with data. We loop through the sample intervals  and within an
        # interval, extract slices of data that share the same number of samples. We then
        # determine if we're expanding or shrinking the number of samples and if the resample
        # factor is a whole number or float.  If it is a whole number and we are expanding
        # we replicate existing sample data to fill out the expanded array. If reducing, we
        # take the mean of the samples. If the resample factor is not a whole number we interpolate.
        for sample_interval in unique_sample_intervals:
            # get an index of the pings with this sample interval
            pings = pings_this_interval[sample_interval]

            # Determine the unique sample_counts for this sample interval.
            unique_sample_counts = np.unique(sample_counts[pings])

            # check if the resample factor is an whole number. When it is a whole nummber, we
            # can reduce or expand the array, if it is a float we have to interpolate.
            if isclose(resample_factor[sample_interval], round(resample_factor[sample_interval])):
                # it is, we will replicate samples when expanding and take the mean when reducing
                use_interp = False
            else:
                # the resample factor is a float - we need to interpolate
                use_interp = True

            for count in unique_sample_counts:
                # get an index into the pings with this sample count
                p_n_samples = sample_counts[pings] == count

                # Combine the index array for this sample interval/sample count chunk of data.
                pings_this_interval_count = pings[p_n_samples]

                # Determine if we're reducing, expanding, or keeping the same number of samples.
                if resample_interval > sample_interval:
                    # We're reducing the number of samples.

                    # If we're resampling power, convert power to linear units.
                    if is_power:
                        this_data = 10.0 ** (data[pings_this_interval_count] / 10.0)
                    else:
                        this_data = data[pings_this_interval_count]

                    if use_interp:
                        # We can't reduce by averaging a whole number of samples so we have to interpolate.
                        xp = sample_interval * np.arange(count)
                        rsf = int(count * resample_factor[sample_interval])
                        yp = resample_interval * np.arange(rsf)
                        interp_f = interp1d(xp, this_data[:,0:count], kind='previous', axis=1,
                                bounds_error=False, fill_value=np.nan, assume_sorted=True)
                        this_data = interp_f(yp)
                    else:
                        # Reduce the number of samples by taking the mean.
                        n_mean = int(resample_factor[sample_interval])
                        this_data = np.mean(this_data.reshape(-1, n_mean), axis=1)

                    if is_power:
                        # Convert power back to log units.
                        this_data = 10.0 * np.log10(this_data)

                elif resample_interval < sample_interval:
                    # We're increasing the number of samples.

                    if use_interp:
                        # If we're resampling power, convert power to linear units.
                        if is_power:
                            this_data = 10.0 ** (data[pings_this_interval_count] / 10.0)
                        else:
                            this_data = data[pings_this_interval_count]

                        # We can't expand by replicating a whole number of samples so we have to interpolate.
                        xp = sample_interval * np.arange(count)
                        rsf = int(count * resample_factor[sample_interval])
                        yp = resample_interval * np.arange(rsf)
                        interp_f = interp1d(xp, this_data[:,0:count], kind='previous', axis=1,
                                bounds_error=False, fill_value=np.nan, assume_sorted=True)
                        this_data = interp_f(yp)

                        if is_power:
                            # Convert power back to log units.
                            this_data = 10.0 * np.log10(this_data)

                    else:
                        # Replicate the values to fill out the higher resolution array.
                        n_repeat = int(resample_factor[sample_interval])
                        this_data = data[pings_this_interval_count]
                        if data.ndim == 3:
                            this_data = np.repeat(this_data[:, 0:count,:], n_repeat, axis=1)
                        else:
                            this_data = np.repeat(this_data[:, 0:count], n_repeat, axis=1)

                else:
                    # The data exists on the resample_interval grid - no change
                    this_data = data[pings_this_interval_count, 0:count]

                # Assign new values to output array.  At the same time, we will shift the data by sample offset.
                unique_sample_offsets = np.unique(sample_offsets_this_interval[sample_interval]).astype('int')
                for offset in unique_sample_offsets:
                    if this_data.ndim == 3:
                        resampled_data[pings_this_interval_count, offset:offset + this_data.shape[1],:] = this_data
                    else:
                        resampled_data[pings_this_interval_count, offset:offset + this_data.shape[1]] = this_data

        # Return the resampled data and the sampling interval used.
        return resampled_data, resample_interval


    def _vertical_shift(self, data, sample_offsets, unique_sample_offsets,
                        min_sample_offset):
        """Adjusts the output array size and pads the top of the samples
        array to vertically shift the positions of the sample data in
        the output array.

        Pings with offsets greater than the minimum will be padded on the
        top, shifting them into their correct location relative to the other
        pings.  The result is an output array with samples that are properly
        aligned vertically relative to each other with a sample offset that is
        constant and equal to the minimum of the original sample offsets.

        This method is only called if our data has a constant sample interval,
        but varying sample offsets. If the data has multiple sample intervals
        the offset adjustment is done in vertical_resample.

        Args:
            data (array): A numpy array of data to be shifted.
            sample_offsets (array): A numpy array with the sample offset for
                each ping.
            unique_sample_offsets (list): The lis tof unique sample offset
                values.
            min_sample_offset (int):

        Returns:
            The shifted data array.
        """

        # Determine the new array size.
        new_sample_dims = (data.shape[1] + max(sample_offsets) -
                min_sample_offset)

        # Create the new array.
        shifted_data = np.empty((data.shape[0], new_sample_dims),
                dtype=self.sample_dtype, order='C')
        shifted_data.fill(np.nan)

        # Fill the array, looping over the different sample offsets.
        for offset in unique_sample_offsets:
            rows_this_offset = np.where(sample_offsets == offset)[0]
            start_index = offset - min_sample_offset
            end_index = start_index + data.shape[1]
            shifted_data[rows_this_offset, start_index:end_index] = \
                    data[rows_this_offset, 0:data.shape[1]]

        return shifted_data


    def _copy(self, obj):
        """Copies attributes.

        This is an internal helper method that is called by child "copy"
        methods to copy the data and object attributes.

        Args:
            obj (ping_data): The object to copy attributes to.

        Returns:
            The copy of the object.
        """

        # Copy the common attributes.
        obj.sample_dtype = self.sample_dtype
        obj.n_samples = self.n_samples
        obj.n_pings = self.n_pings
        obj.shape = self.shape
        obj._data_attributes = list(self._data_attributes)
        obj._object_attributes  = list(self._object_attributes)


        # Copy object attributes
        for attr_name in self._object_attributes:
            attr = getattr(self, attr_name)
            # check if this attribute is a numpy array
            if isinstance(attr, np.ndarray):
                # it is - use ndarray's copy method
                setattr(obj, attr_name, attr.copy())
            else:
                # it's not - use the copy module
                setattr(obj, attr_name, copy.deepcopy(attr))

        # Copy the data attributes
        for attr_name in obj._data_attributes:
            attr = getattr(self, attr_name)
            # data attributes are always numpy arrays so use ndarray's copy method
            setattr(obj, attr_name, attr.copy())

        # Return the copy.
        return obj


    def _shape(self):
        '''Internal method used to update the shape attribute
        '''
        shape = None
        if hasattr(self, 'power'):
            shape = self.power.shape
        elif hasattr(self, 'angles_alongship_e'):
            shape = self.angles_alongship_e.shape
        elif hasattr(self, 'complex'):
            shape = self.complex.shape
        elif hasattr(self, 'data'):
            shape = self.data.shape
        return shape


    def _like(self, obj, n_pings, value, empty_times=False, no_data=False):
        """Copies ping_data attributes and creates data arrays filled with the
        specified value.

        This is an internal helper method that is called by "empty_like" and
        "zeros_like" methods of child classes which copy the ping_data
        attributes into the provided ping_data based object as well as
        create "data" arrays that are filled with the specified value. All
        vertical axes will be copied without modification.

        If empty_times is False, the ping_time vector of this instance is copied
        to the new object. If it is True, the new ping_time vector is filled
        with NaT (not a time) values. If n_pings != self.n_pings THIS
        ARGUMENT IS IGNORED AND THE NEW PING VECTOR IS FILLED WITH NaT.

        The result should be a new object where horizontal axes (excepting
        ping_time) and sample data arrays are empty (NaN or NaT). The
        contents of the ping_time vector will depend on the state of the
        empty_times keyword. The new object's shape will be (n_pings,
        self.n_samples).

        Args:
            obj (ping_data): An empty object to copy attributes to.
            n_pings (int): Number of pings (horizontal axis)
            value (int,float): A scalar value to fill the array with.
            empty_times (bool): Controls whether ping_time data is copied
                over to the new object (TRUE) or if it will be filled with NaT
                values (FALSE).
            no_data (bool): Set to True to to set 2d and 3d data attributes
                to None, rather than creating numpy arrays. When False,
                numpy arrays are created. This allows you to avoid allocating
                the data arrays if you are planning on replacing them.
                This is primarily used internally. Default: False

        Returns:
            The object copy, obj.
        """
        # If n_pings is None, we create an empty array with the same number
        # of pings.
        if n_pings is None:
            n_pings = self.n_pings

        # Copy the common attributes.
        obj.sample_dtype = self.sample_dtype
        obj.n_samples = self.n_samples
        obj.n_pings = n_pings
        obj._data_attributes = list(self._data_attributes)
        obj._object_attributes  = list(self._object_attributes)

        # Copy object attributes - this is simple as there are no
        # size or type checks.
        for attr_name in self._object_attributes:
            attr = getattr(self, attr_name)
            # check if attribute is a numpy array
            if isinstance(attr, np.ndarray):
                # it is - use ndarray's copy method
                setattr(obj, attr_name, attr.copy())
            else:
                # it's not - use the copy module
                setattr(obj, attr_name, copy.deepcopy(attr))

        # Check if n_pings != self.n_pings.  If the new object's horizontal
        # axis is a different shape than this object's we can't copy
        # ping_time data since there isn't a direct mapping and we don't know
        # what the user wants here. This can/should be handled in the child
        # method if needed.
        if n_pings != self.n_pings:
            # We have to force an empty ping_time vector since the axes differ.
            empty_times = True

        # Create the dynamic attributes.
        for attr_name in self._data_attributes:

            # Get the attribute.
            attr = getattr(self, attr_name)
            
            if isinstance(attr, np.ndarray):

                if attr.shape[0] == self.n_samples:
                    # Copy all vertical axes w/o changing them.
                    data = attr.copy()
                else:
                    # Create an array with the appropriate shape filled with the
                    # specified value.
                    if attr.ndim == 1:
                        # Create an array with the same shape filled with the
                        # specified value.
                        data = np.empty(n_pings, dtype=attr.dtype)

                        # Check if this is the ping_time attribute and if we
                        # should copy this instance's ping_time data or create an
                        # empty ping_time vector
                        if attr_name == 'ping_time':
                            if empty_times:
                                data[:] = np.datetime64('NaT')
                            else:
                                data[:] = attr.copy()
                        elif data.dtype == 'datetime64[ms]':
                            data[:] = np.datetime64('NaT')
                        elif np.issubdtype(data.dtype, np.integer):
                            data[:] = 0
                        else:
                            data[:] = value
                    else:
                        # Check if we're supposed to create the sample data arrays
                        if no_data:
                            # No - we'll set them to None assuming the user will set them
                            data = None
                        else:
                            # Yes, create the data arrays
                            if attr.ndim == 2:
                                # Create the 2d array(s).
                                data = np.empty((n_pings, self.n_samples), dtype=attr.dtype)
                                data[:, :] = value
                            elif attr.ndim == 3:
                                #  must be a 3d attribute
                                data = np.empty((n_pings, self.n_samples, self.n_complex),
                                    dtype=attr.dtype)
                                data[:, :, :] = value
            else:
                #  For now we don't duplicate attributes that aren't numpy arrays
                data = None

            # Add the attribute to our empty object.  We can skip using
            # add_data_attribute here because we shouldn't need to check
            # dimensions and we've already handled the low level stuff like
            # copying the _data_attributes list, etc.
            setattr(obj, attr_name, data)

        return obj




    # def _center_weighted_resample(self, n_intervals, n_layers, interval_pings, ping_interval_map,
    #         layer_samples, sample_layer_map):

    #     resamp_data = np.empty((n_intervals,n_layers), dtype=self.data.dtype)

    #     for i in range(n_intervals):
    #         ping_idx = ping_interval_map==i
    #         for j in range(n_layers):
    #             cell_data = self.data[np.ix_(ping_idx,sample_layer_map==j)]
    #             n_nans_in_cell = np.count_nonzero(np.isnan(cell_data))
    #             good_samps = (interval_pings[i] * layer_samples[j]) - n_nans_in_cell
    #             if good_samps > 0:
    #                 resamp_data[i,j] = np.nansum(cell_data) / good_samps
    #             else:
    #                 resamp_data[i,j] = np.nan

    #     return resamp_data



    # def _area_weighted_resample(self, v_axis_data, h_axis_data, interval_edges, layer_edges):

    #     v_axis_edges = np.append(v_axis_data - self.sample_thickness / 2, 
    #             v_axis_data[-1] + self.sample_thickness / 2)
    #     h_axis_thickness = h_axis_data[-1] - h_axis_data[-2]
    #     h_axis_edges = np.append(h_axis_data, h_axis_data[-1] + h_axis_thickness)

    #     h_weights = self._get_weights_1d(h_axis_edges, interval_edges)
    #     v_weights = self._get_weights_1d(v_axis_edges, layer_edges)

    #     #  finally resample the sample data
    #     self.to_linear()

    #     #  create the NaN mask and zero out NaNs
    #     mask = np.isnan(self.data)
    #     self.data[mask] = 0
    #     #  invert the mask
    #     np.logical_not(mask, out=mask)

    #     resampled_num = h_weights @ self.data @ v_weights.T
    #     resampled_den = h_weights @ mask.astype('float') @ v_weights.T
    #     resampled_data = np.divide(resampled_num, resampled_den, 
    #                        out=np.full_like(resampled_num, np.nan), 
    #                        where=resampled_den > 0)


    #     return resampled_data

    '''
    def view_as_blocks(arr_in, block_shape):

        from numpy.lib.stride_tricks import as_strided

        """Block view of the input n-dimensional array (using re-striding).

        Blocks are non-overlapping views of the input array.

        Parameters
        ----------
        arr_in : ndarray, shape (M[, ...])
            Input array.
        block_shape : tuple
            The shape of the block. Each dimension must divide evenly into the
            corresponding dimensions of `arr_in`.

        Returns
        -------
        arr_out : ndarray
            Block view of the input array.

        Examples
        --------
        >>> import numpy as np
        >>> from skimage.util.shape import view_as_blocks
        >>> A = np.arange(4*4).reshape(4,4)
        >>> A
        array([[ 0,  1,  2,  3],
               [ 4,  5,  6,  7],
               [ 8,  9, 10, 11],
               [12, 13, 14, 15]])
        >>> B = view_as_blocks(A, block_shape=(2, 2))
        >>> B[0, 0]
        array([[0, 1],
               [4, 5]])
        >>> B[0, 1]
        array([[2, 3],
               [6, 7]])
        >>> B[1, 0, 1, 1]
        13

        >>> A = np.arange(4*4*6).reshape(4,4,6)
        >>> A  # doctest: +NORMALIZE_WHITESPACE
        array([[[ 0,  1,  2,  3,  4,  5],
                [ 6,  7,  8,  9, 10, 11],
                [12, 13, 14, 15, 16, 17],
                [18, 19, 20, 21, 22, 23]],
               [[24, 25, 26, 27, 28, 29],
                [30, 31, 32, 33, 34, 35],
                [36, 37, 38, 39, 40, 41],
                [42, 43, 44, 45, 46, 47]],
               [[48, 49, 50, 51, 52, 53],
                [54, 55, 56, 57, 58, 59],
                [60, 61, 62, 63, 64, 65],
                [66, 67, 68, 69, 70, 71]],
               [[72, 73, 74, 75, 76, 77],
                [78, 79, 80, 81, 82, 83],
                [84, 85, 86, 87, 88, 89],
                [90, 91, 92, 93, 94, 95]]])
        >>> B = view_as_blocks(A, block_shape=(1, 2, 2))
        >>> B.shape
        (4, 2, 3, 1, 2, 2)
        >>> B[2:, 0, 2]  # doctest: +NORMALIZE_WHITESPACE
        array([[[[52, 53],
                 [58, 59]]],
               [[[76, 77],
                 [82, 83]]]])
        """
        if not isinstance(block_shape, tuple):
            raise TypeError('block needs to be a tuple')

        block_shape = np.array(block_shape)
        if (block_shape <= 0).any():
            raise ValueError("'block_shape' elements must be strictly positive")

        if block_shape.size != arr_in.ndim:
            raise ValueError("'block_shape' must have the same length " "as 'arr_in.shape'")

        arr_shape = np.array(arr_in.shape)
        if (arr_shape % block_shape).sum() != 0:
            raise ValueError("'block_shape' is not compatible with 'arr_in'")

        # -- restride the array to build the block view
        new_shape = tuple(arr_shape // block_shape) + tuple(block_shape)
        new_strides = tuple(arr_in.strides * block_shape) + arr_in.strides

        arr_out = as_strided(arr_in, shape=new_shape, strides=new_strides)

        return arr_out
    '''