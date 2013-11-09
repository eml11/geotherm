geotherm
========

1d time dependent geothermal model for netcdf input

Clone github repo with https://github.com/eml11/geotherm.git
Github repo at https://github.com/eml11/geotherm

Dependecies
  - gfortran compiler
  - netcdf-fortran

Tested on mac os x 10.9

To intall:
  - edit LDFLAGS in Makefile to include directory for netcdff and
    netcdf libaries if these are in a non-standard location
  - edit FDFLAGS in Makefile to include directory for netcdf.mod
    if this is in a non-standard location 
  - type make


Source code located in src/
Modules and Objects located in f90/
executible located in bin/

Supported Netcdf formats:
  - xyz  where x is time, y is depth, z is dependent variable
  - currently input only in depth possitive format


Author: Elliot Lynch
