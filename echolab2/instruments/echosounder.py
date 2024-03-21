# coding=utf-8

#    National Oceanic and Atmospheric Administration
#    Alaskan Fisheries Science Center
#    Resource Assessment and Conservation Engineering
#    Midwater Assessment and Conservation Engineering

# THIS SOFTWARE AND ITS DOCUMENTATION ARE CONSIDERED TO BE IN THE PUBLIC DOMAIN
# AND THUS ARE AVAILABLE FOR UNRESTRICTED PUBLIC USE. THEY ARE FURNISHED "AS
# IS." THE AUTHORS, THE UNITED STATES GOVERNMENT, ITS INSTRUMENTALITIES,
# OFFICERS, EMPLOYEES, AND AGENTS MAKE NO WARRANTY, EXPRESS OR IMPLIED, AS TO
# THE USEFULNESS OF THE SOFTWARE AND DOCUMENTATION FOR ANY PURPOSE.  THEY
# ASSUME NO RESPONSIBILITY (1) FOR THE USE OF THE SOFTWARE AND DOCUMENTATION;
# OR (2) TO PROVIDE TECHNICAL SUPPORT TO USERS.

'''
.. module:: echolab2.instruments.echosounder

    :synopsis:  The top level interface for reading data files collected
                using scientific echosunders commonly used in fisheries
                stock assessment and research.


| Developed by:  Rick Towler   <rick.towler@noaa.gov>
| National Oceanic and Atmospheric Administration (NOAA)
| Alaska Fisheries Science Center (AFSC)
| Midwater Assesment and Conservation Engineering Group (MACE)
|
| Authors:
|       Rick Towler   <rick.towler@noaa.gov>
| Maintained by:
|       Rick Towler   <rick.towler@noaa.gov>

$Id$
'''

import os
from . import EK80
from . import EK60
from .util import simrad_utils


SIMRAD_EK60 = 0
SIMRAD_EK80 = 1


def read(files, ignore_errors=False, **kwargs):
    """

    """

    # if we're given a string, wrap it in a list
    if not isinstance(files, list):
            files = [files]

    data_object = None

    # Work through the list of input files
    for index, item in enumerate(files):
        filename = os.path.normpath(item)

        # Determine what kind of data file we have
        if data_object is None:
            data_type = _check_filetype(filename)

            if data_type == SIMRAD_EK60:
                # This is an EK60 file, we're going to assume all files are EK60
                data_object = EK60.EK60()

            elif data_type == SIMRAD_EK80:
                # This is an EK80 file, we're going to assume all files are EK80
                data_object = EK80.EK80()
            else:
                # We don't know what this is
                raise UnknownFormatError("Unknown file type encountered: " + filename)

        # Read the data file using the object for this data type
        data_object.append_raw(filename, **kwargs)

        #  now try to read the bottom data. We are assuming the .bot or .XYZ files are colocated
        #  with the .raw file and that they follow the normal Simrad naming convention.
        if data_type == SIMRAD_EK80:
            #  Simrad EK80 systems can generate both XYZ and.or .bot files
            bot_type, bot_files = simrad_utils.get_simrad_bottom_files(filename, data_object)
            if type == 'BOT':
                data_object.read_bot(bot_files)
            elif type == 'XYZ':
                for channel_id in bot_files:
                    data_object.read_xyz(channel_id, bot_files[channel_id])

        elif data_type == SIMRAD_EK60:
            #  Simrad EK60 systems do not generate XYZ files so we only check for .bot files
            _, bot_files = simrad_utils.get_simrad_bottom_files(filename, data_object,
                skip_xyz=True)
            if bot_files:
                data_object.read_bot(bot_files)


    return data_object



def _check_filetype(filename):

    # Read in the file header, this value can change if additional
    # instruments are introduced with larger headers.
    fh = open(filename, "rb")
    header = fh.read(8)
    fh.close()

    # Return the file type based on the header
    if header[4:8]==b'CON0':
        # Simrad EK60 style raw files start  with 4 bytes for the size then C O N 0
        return SIMRAD_EK60

    elif header[4:8]==b'XML0':
        # Simrad EK80 style raw files start with 4 bytes for the size then X M L 0
        return SIMRAD_EK80

    else:
       return -1



class UnknownFormatError(Exception):
    pass
