# pyEcholab

## UNDER CONSTRUCTION

pyEcholab is a python package for reading, writing, processing, and plotting data from Simrad/Kongsberg sonar systems.

pyEcholab reads .raw data files created by the Simrad/Kongsberg EK/ES60, ES70, and EK/ES80 echjosounders, as well as the ME70 multibeam sonar, and many other Kongsberg sonar systems that output .raw files. The library provides direct access to the all data contained within the file allowing you to 



### Current Classes/Functionalities

&nbsp;&nbsp;**_Data Classes_**  
* PingData: Stores "ping" based data from fisheries sonar systems. Parent class for the ProcessedData and RawData classes.
* ProcessedData: Stores and manipulates a 2d sample data array along with 1d arrays that define the horizontal and vertical axes of that data. The horizontal axis is defined as 'ping_time' and the vertical axis is 'range' or 'depth'
* RawData: Contains a single channel’s data extracted from a Simrad raw file. This includes all data contained within the .raw file.

&nbsp;&nbsp;**_Reader Classes_**
*echosounder: The echosounder class provides an simple interface for reading .raw files.
* EK60: Reads in one or more EK60 files and generates a RawData class instance for each unique channel ID in the files. The class also contains numerous methods used to extract raw sample data or create ProcessedData objects containing transformed data such as Sv.
*  EK80: Reads in one or more EK8 files and generates a RawData class instance for each unique channel ID in the files. The class also contains numerous methods used to extract raw sample data or create ProcessedData objects containing transformed data such as Sv, TS, Sv(f)

&nbsp;&nbsp;**_Processing & Analysis Classes_**
* Mask: Creates, manipulates and applies masks to ProcessedData objects
* Line: Implements lines based on ping_time and depth/range values. It provides methods manipulating these values in various ways
* AlignPings: Operates on grouped ProcessedData or RawData objects aligning ping times across objects by either adding blank pings or deleting pings so all objects match the object with the longest or shortest pingtime array

&nbsp;&nbsp;**_Plotting and Export Classes_** 
* EchogramViewer: Uses PyQt to create interactive plot display of data objects. 
* Echogram: Provides basic plotting functions to display data objects using matplotlib

## Getting Started


pyEcholab and its core library echolab2 are currently hosted in this GitHub
repository (https://github.com/noaa-afsc-mace/pyEcholab) and are not currently 
available from any of the Python package management systems.
Download to your machine by cloning the git repository.

```
mkdir ~/<where_you_want_the_code>
cd ~/<where_you_want_the_code>
git clone https://github.com/noaa-afsc-mace/pyEcholab.git
```

### Prerequisites

pyEcholab2 requires  Python 3.6 or greater. For those not familiar with building Python environments, we recommend using an [Anaconda](https://www.anaconda.com/download/) Python distribution which includes all needed packages except basemap. The following packages are needed If you use your own Python installation: 

&nbsp;&nbsp;**_Required_**
* [matplotlib](https://matplotlib.org/) for plotting echograms.
* [numpy](http://www.numpy.org/) for large, multi-dimensional arrays.
* [pytz](http://pytz.sourceforge.net/) for cross platform timezone calculations.
* [future](https://pypi.org/project/future/) compatibility layer between Python 2 and Python 3.

&nbsp;&nbsp;**_Optional_**
* [PyQT4](https://wiki.python.org/moin/PyQt4) for GUI applications (see [example](https://github.com/CI-CMG/PyEcholab2/blob/master/examples/qt_echogram_viewer.py)).
* [basemap](https://matplotlib.org/basemap/) for plotting on maps (only used in [nmea example](https://github.com/CI-CMG/PyEcholab2/blob/master/examples/nmea_example.py) and currently only works with matplotlib 1.5.0rc3, basemap 1.0.8, and pyproj 1.9.5.1).
* [cartopy](https://scitools.org.uk/cartopy/docs/v0.15/installing.html#installing) for plotting on maps (replacing basemap).

### Installation

Once you have pyEcholab cloned to your machine, you need to tell Python where to find the package. Your Python development environment may provide a way to do this, if not, the following can be used at a bash prompt or added to
your shell initialisation file:

```
export PYTHONPATH=$PYTHONPATH:~/<where_you_want_the_code>
```

### Data
Raw files used by the example scripts need to be downloaded seperately to the examples/data folder before running. They can be downloaded from the NCEI FTP server: *ftp://ftp.ngdc.noaa.gov/pub/outgoing/mgg/wcd/pyEcholab_data/examples/*  
  
On Unix based systems, this can easily achieved as follows:

```
cd ~/<where_you_want_the_code>/examples/data
make
```

### Examples
Numerous example files are provided to introduce users to pyEcholab. These examples exist as both Python scripts and Jupyter Notebooks. The examples are heavily commented to explain each step.


## Authors

* **Rick Towler** - _NOAA/NMFS, Alaska Fisheries Science Center_
* **Zac Berkowitz** -_NOAA/NMFS, Alaska Fisheries Science Center_
* **Robert Levine** - _NOAA/NMFS, Alaska Fisheries Science Center_


## License

This software is licensed under the MIT License - see the
[LICENSE](LICENSE) file for details

