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
  - edit LDFLAGS in Makefile to include directory for netcdff libary
  - edit FDFLAGS in Makefile to include directory for netcdf.mod 
  - type make

Source code located in src/
Modules and Objects located in f90/
executible located in bin/

Author:Elliot Lynch
