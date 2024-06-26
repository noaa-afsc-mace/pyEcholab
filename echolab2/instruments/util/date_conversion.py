﻿# coding=utf-8

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

'''
.. module:: echolab.instruments.util.unit_conversion


useful functions:

    nt_to_unix
    unix_to_nt

    datetime_to_unix
    unix_to_datetime


| Developed by:  Zac Berkowitz <zac.berkowitz@gmail.com> under contract for
| National Oceanic and Atmospheric Administration (NOAA)
| Alaska Fisheries Science Center (AFSC)
| Midwater Assesment and Conservation Engineering Group (MACE)
|
| Author:
|       Zac Berkowitz <zac.berkowitz@gmail.com>
| Maintained by:
|       Rick Towler   <rick.towler@noaa.gov>

$Id$
'''

import datetime
import numpy as np
from pytz import utc as pytz_utc
import logging


#NT epoch is Jan 1st 1601
UTC_NT_EPOCH = datetime.datetime(1601, 1, 1, 0, 0, 0, tzinfo=pytz_utc)
#Unix epoch is Jan 1st 1970
UTC_UNIX_EPOCH = datetime.datetime(1970, 1, 1, 0, 0, 0, tzinfo=pytz_utc)

EPOCH_DELTA_SECONDS = (UTC_UNIX_EPOCH - UTC_NT_EPOCH).total_seconds()

__all__ = ['nt_to_unix', 'unix_to_nt']

log = logging.getLogger(__name__)


def dt64_to_datetime(dt64):
    '''
    :param dt64: Numpy datetime64 object to convert to datetime
    :type dt64: datetime64

    Returns a datetime.datetime object representing the same time as the
    provided datetime64 object.

    source:
    https://stackoverflow.com/questions/13703720/converting-between-datetime-timestamp-and-datetime64

    '''

    ts = (dt64 - np.datetime64('1970-01-01T00:00:00')) / np.timedelta64(1, 's')

    return datetime.datetime.fromtimestamp(ts, tz=pytz_utc)


def dt64_to_nt(dt64):
    '''
    :param dt64: Numpy datetime64 object to convert to NT time
    :type dt64: datetime64

    Returns a tuple containing the NT time.

    This method was changed 8/2023 to use integer math to compute the NT time from
    the provided datetime64 object to reduce errors introduced when using FP math.

    '''

    if np.isnat(dt64):
        # NaT values are assumed to be NULL times
        lowDateTime = 0
        highDateTime = 0
    else:
        # Not NaT so we convert to NT time
        ts = (dt64 - np.datetime64('1970-01-01T00:00:00')) / np.timedelta64(1, 's')
        unix_datetime = datetime.datetime.fromtimestamp(ts, tz=pytz_utc)
        time_past_nt_epoch = (unix_datetime - UTC_NT_EPOCH)

        onehundred_ns_intervals = (((time_past_nt_epoch.days * 86400) + time_past_nt_epoch.seconds) *
                10000000) + (time_past_nt_epoch.microseconds * 10)

        lowDateTime = onehundred_ns_intervals & 0xFFFFFFFF
        highDateTime = onehundred_ns_intervals >> 32

    return lowDateTime, highDateTime


def nt_to_unix(nt_timestamp_tuple, return_datetime=True):
    '''
    :param nt_timestamp_tuple: Tuple of two longs representing the NT date
    :type nt_timestamp_tuple: (long, long)

    :param return_datetime:  Return a datetime object instead of float
    :type return_datetime: bool


    Returns a datetime.datetime object w/ UTC timezone
    calculated from the nt time tuple

    lowDateTime, highDateTime = nt_timestamp_tuple

    The timestamp is a 64bit count of 100ns intervals since the NT epoch
    broken into two 32bit longs, least significant first:

    This method was changed 8/2023 to use integer math to compute the NT time from
    the provided datetime64 object to reduce errors introduced when using FP math.
    This change only applies when returning a datetime object.

    '''

    lowDateTime, highDateTime = nt_timestamp_tuple

    if return_datetime:
        us_past_nt_epoch = ((highDateTime << 32) + lowDateTime) // 10
        return UTC_NT_EPOCH + datetime.timedelta(microseconds=us_past_nt_epoch)
    else:
        sec_past_nt_epoch = ((highDateTime << 32) + lowDateTime) * 1.0e-7
        sec_past_unix_epoch = sec_past_nt_epoch - EPOCH_DELTA_SECONDS
        return sec_past_unix_epoch


def unix_to_nt(unix_timestamp):
    '''
    Given a date, return the 2-element tuple used for timekeeping with SIMRAD echosounders

    #converting back may not yield the exact original date,
    #but will be within the datetime's precision

    '''

    if isinstance(unix_timestamp, datetime.datetime):
        if unix_timestamp.tzinfo is None:
            unix_datetime = pytz_utc.localize(unix_timestamp)

        elif unix_timestamp.tzinfo == pytz_utc:
            unix_datetime = unix_timestamp

        else:
            unix_datetime = pytz_utc.normalize(unix_timestamp.astimezone(pytz_utc))
    elif isinstance(unix_timestamp, np.datetime64):
        ts = (unix_timestamp - np.datetime64('1970-01-01T00:00:00Z')) / np.timedelta64(1, 's')
        unix_datetime = unix_to_datetime(ts)
    else:
        unix_datetime = unix_to_datetime(unix_timestamp)

    sec_past_nt_epoch = (unix_datetime - UTC_NT_EPOCH).total_seconds()

    onehundred_ns_intervals = int(sec_past_nt_epoch * 1e7)
    lowDateTime = onehundred_ns_intervals & 0xFFFFFFFF
    highDateTime = onehundred_ns_intervals >> 32

    return lowDateTime, highDateTime


def unix_to_datetime(unix_timestamp, tz=pytz_utc):
    '''
    :param unix_timestamp: Number of seconds since unix epoch (1/1/1970)
    :type unix_timestamp: float

    :param tz: timezone to use for conversion (default None = UTC)
    :type tz: None or tzinfo object (see datetime docs)

    :returns: datetime object
    :raises: ValueError if unix_timestamp is not of type float or datetime

    Returns a datetime object from a unix timestamp.  Simple wrapper for
    :func:`datetime.datetime.fromtimestamp`

    '''

    if isinstance(unix_timestamp, datetime.datetime):
        if unix_timestamp.tzinfo is None:
            unix_datetime = pytz_utc.localize(unix_timestamp)

        elif unix_timestamp.tzinfo == pytz_utc:
            unix_datetime = unix_timestamp
        else:
            unix_datetime = pytz_utc.normalize(unix_timestamp.astimezone(pytz_utc))

    elif isinstance(unix_timestamp, float):
        unix_datetime = pytz_utc.localize(datetime.datetime.fromtimestamp(unix_timestamp))

    else:
        errstr = 'Looking for a timestamp of type datetime.datetime or # of sec past unix epoch.\n'
        errstr += 'Supplied timestamp \'%s\' of type %s.' % \
            (str(unix_timestamp), type(unix_timestamp))
        raise ValueError(errstr)

    return unix_datetime


def datetime_to_unix(datetime_obj):
    '''
    :param datetime_obj: datetime object to convert
    :type datetime_obj: :class:`datetime.datetime`

    :param tz: Timezone to use for converted time -- if None, uses timezone
                information contained within datetime_obj
    :type tz: :class:datetime.tzinfo
    '''

    timestamp = (datetime_obj - UTC_UNIX_EPOCH).total_seconds()

    return timestamp


