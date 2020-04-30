**This is an info file about the medicanes tracking algorithm <br>
developed in the MAR group of the University of Murcia in 2019. <br>
<br>
Authors: Enrique Pravia-Sarabia, Juan Pedro Montávez, Juan José Gómez-Navarro and Pedro Jiménez-Guerrero. <br>
<br>
Release date: December 2019. <br>
Last updated: March 2020.** <br>

### PART I: TRACKING ALGORITHM

To run the tracking algorithm, you should go to the *Code* folder and type the following command in the console:

```console
user@machine:dir$ bash findmedicanes.sh TRACKINGINPUT NCORES 
```

**TRACKINGINPUT** represents the path to a WRF output nc file which contains all the fields, timesteps and vertical levels, or a folder with the requested fields if they have been previously interpolated to pressure fields (see paragraph Requested fields below). Path to file or folder is to be introduced with no quotation marks. If path to WRF output file, it can be relative. If path to folder, an absolute path is mandatory. **NCORES** is the number of cores to be used in the parallelized calculation of the medicanes track. It is an optional arguments that defaults to the machine total number of cores if argument is not passed, saving one core to prevent overloading. 

Algorithm parameters can be changed in the *FindMedicanes.namelist* file inside the TrackingAlgorithm folder. This file must be in that same location when running the findmedicanes bash script. The *FindMedicanes.namelist.README* file contains detailed information about the meaning of each parameter, its possible values and its recommended/default value. 

#### Requested fields
In case of executing the algorithm over data different from WRF output, the input provided to the algorithm must be a folder  with four files called 'outputfile-slp.nc', 'outputfile-z.nc', 'outputfile-uvmet10-U.nc' and 'outputfile-uvmet10-V.nc', containing SLP, geopotential height in pressure levels, and 10-m wind U and V components. The SLP can be either in Pa or hPa, but the value of SLP threshold must be set accordingly in the TrackingAlgorithm/FindMedicanes.namelist file. Geopotential height must be in meters, and 10-m winds in km.h{-1}. If pressure levels for the geopotential height are in Pa, remember to change values of pressure levels for the calculation of Hart parameters in the TrackingAlgorithm/FindMedicanes.namelist file. The variables and dimensions names of the netCDF, as well as the spatial and temporal resolutions, need to be accordingly set in the TrackingAlgorithm/FindMedicanes.namelist file. 

**Additional notes**: <br>
<br>
The variables and dimensions names of the netCDF are contained in the netCDF metadata. The command ncdump -h of ncdump tool is recommended for this particular issue. In R, the ncdf4::nc_open() command can be useful. It will load the netCDF as a S4 object, with dimensions and variables saved in S4object$dim and S4object$var elements. 

Climate Data Operators, most known as CDO [https://code.mpimet.mpg.de/projects/cdo/embedded/cdo.pdf], are a useful tool for editing a netCDF. It will allow easy and quick modifications of units, spatial and temporal cropping on the netCDF before passing the algorithm. Nevertheless, the metadata of the netCDF is changed by CDO. Thus, if the WRF-python based pinterpy preprocessing tool is to be used (Code/pinterpy folder), it is highly recommended to use the NCO [http://nco.sourceforge.net/nco.pdf] operators instead. 

The name given to the output folder once the algorithm has terminated its execution is 'completed'. If another 'completed' folders are present in the running directory, then a number will be added to avoid overwriting. Output of this tool is a list object of R -*.RData*- with the found centers, which will be located inside the 'completedX' folder. 

#### System requirements

**Python version**: Python 3.6.8 (a correct working of pinterpy interpolation tool with previous python versions is expected but not guaranteed). <br>
**Python requested libraries**: <br>
   netCDF4 <br>
   wrf -> https://wrf-python.readthedocs.io/en/latest/ <br>
   ast <br>
   numpy <br>
   xarray <br>

**R version**: R 3.6.1 (a correct installation of requested packages in previous R versions is not guaranteed). <br>
**R requested packages**: <br>
   ncdf4 <br>
   oce <br>
     > install.packages('oce', dependencies=TRUE, repos='http://cran.rstudio.com/') <br>
       or <br>
     > install_github('dankelley/oce', ref='develop') <br>
   foreach <br>
   doParallel <br>


### PART II: TRACK POSTPROCESSING.

There are two postprocessing tools provided in this package, one for producing a plot with the track, and one for obtaining further data on the track parameters.

#### Plot of calculated track

To plot the calculated track -SEE PART I-, you should go to the *Code* folder and introduce the following command in the console:

```console
user@machine:dir$ bash plotmedicanestrack.sh PATHTOFOLDER ADJUST COMPLETE CONNECT SLPEXPANDN DTHRESHOLD DTTHRESHOLD
```

**PATHTOFOLDER** must be the path to the folder created by the PART I. Inside that folder, it should be present the track file -output of findmedicanes.sh-, along with a folder called 'output' containing the vertically interpolated files -slp, uvmet10 and z- in netcdf format. PATHTOFOLDER is to be introduced with no quotation marks. **ADJUST** should be one of TRUE/FALSE, and determines if the output map extent should be adjusted to the track, or the complete domain should be plotted instead. **COMPLETE** should be one of TRUE/FALSE, and determines if the track should be complemented with the SLP minimum path. **CONNECT** refers to the possibility of connecting track points trough a continuous line. SLPEXPANDN, RESOLUTION, DTHRESHOLD or DTTHRESHOLD should also be introduced as numbers. **SLPEXPANDN** is the number of track points that should be taken to expand the track at its beginning and ending, and is not necessary if complete is set to FALSE. **DTHRESHOLD** refers to the maximum allowed distance between two track points to be connected -in km-, and a value of DTHRESHOLD below 200 km is highly recommended for physical reasons if timestep is 1 hour, 200km*timestep for other timestep values; **DTTHRESHOLD** is the maximum allowed time interval -in timesteps- between two points to be connected. A value of 1 or 2 is recommended for this last parameter. For example, with a value of 1 only consecutive points in time are connected. RESOLUTION, DTHRESHOLD and DTTHRESHOLD are optional arguments. If not provided, the plot won't have connected track dots, and only sparse points will be shown.   

Output of this tool is a pdf called trackplot.pdf, located inside the provided folder. 

#### Obtain further data along calculated track

Go to the *Code* folder and type in the console the following command:

```console
user@machine:dir$ bash getmedicanestrackdata.sh PATHTOFOLDER TYPE
```

**PATHTOFOLDER** should be the path to the folder created in the PART I. Inside that folder, it should be present the *-track.RData file -output of findmedicanes.sh-, along with a folder called 'output' containing the vertically interpolated files -slp, uvmet10 and z- in netcdf format. PATHTOFOLDER is to be introduced with no quotation marks. **TYPE** should be 'reduced' or 'complete' -the quotation marks are allowed but not needed'. 

Output of this tool is a csv called trackingdf.csv, located inside the provided folder.

The 'reduced' TYPE value gives the following track parametersS:

**Date**:                              date and time in the format YY-MM-DD hh:mm:ss <br>
**Timestep**:                          corresponding timestep to the date in the simulation <br>
**x**:                                 coordinate of the center position in the longitudinal dimension <br>
**Lon**:                               longitudinal coordinate of the center position -in degrees- <br>
**y**:                                 coordinate of the center position in the latitudinal dimension <br>
**Lat**:                               latitudinal coordinate of the center position -in degrees- <br>
**center.slp.hPa**:                    SLP value at the center position -in hPa- <br>

The 'complete' TYPE value gives the following track parameters:

**Date**:                              date and time in the format YY-MM-DD hh:mm:ss <br>
**Timestep**:                          corresponding timestep to the date in the simulation <br>
**x**:                                 coordinate of the center position in the longitudinal dimension <br>
**Lon**:                               longitudinal coordinate of the center position -in degrees- <br>
**y**:                                 coordinate of the center position in the latitudinal dimension <br>
**Lat**:                               latitudinal coordinate of the center position -in degrees- <br>
**center.slp.hPa**:                    SLP value at the center position -in hPa- <br>
**inner.radius.km**:                   Medicane inner radius calculated as the distance from the center to the eyewall -max wind surface- <br>
**outer.radius.km**:                   Medicane outer radius calculated as the mean distance from the center to the zero vorticity closed line points <br>
**B.m**:                               B Hart thermal symmetry parameter value at the center position -in meters- <br>
**LTW**:                               LTW Hart lower thermal wind parameter value at the center position <br>
**UTW**:                               UTW Hart upper thermal wind parameter value at the center position <br>
**max.inZVradius.wind.kmh**:           maximum wind speed value inside a circle with the zero vorticity radius around the center position -in km.h{-1}- <br>
**center.wind.kmh**:                   wind speed value at the center position -in km.h{-1}- <br>
**min.inZVradius.wind.kmh**:           minimum wind speed value inside a circle with the zero vorticity radius around the center position -in km.h{-1}- <br>
**min.inZVradius.slp.x.position**:     x position of the minimum SLP value inside a circle with the zero vorticity radius around the center position -in km.h{-1}- <br>
**min.inZVradius.slp.Lon**:            longitude of the minimum SLP value inside a circle with the zero vorticity radius around the center position -in km.h{-1}- <br>
**min.inZVradius.slp.y.position**:     y position of the minimum SLP value inside a circle with the zero vorticity radius around the center position -in km.h{-1}- <br>
**min.inZVradius.slp.Lat**:            latitude of the minimum SLP value inside a circle with the zero vorticity radius around the center position -in km.h{-1}- <br>
**min.inZVradius.slp.hPa**:            minimum SLP value inside a circle with the zero vorticity radius around the center position -in hPa- <br>
**minSLP.center.distance.km**:         distance between the min.inZVradius.slp.hPa position and the center position -in km- <br>
**Saffir.Simpson.scale.category**:     storm category according to the Saffir-Simpson scale using the max.100km.wind.kmh wind speed




