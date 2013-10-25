SRC=$(PWD)/src
OBJECTDIR=$(PWD)/f90
BIN=$(PWD)/bin
LDFLAGS=-L/distributions/netcdf-fortran-4.2/fortran/.libs -L/opt/local/lib -L/distributions/netcdf-fortran-4.2/f90 -lnetcdf -lnetcdff
FDFLAGS=-I/distributions/netcdf-fortran-4.2/f90 -I$(OBJECTDIR)
OUT=geotherm
OBJECTS=$(OBJECTDIR)/mathmodule.o $(OBJECTDIR)/helpermodule.o $(OBJECTDIR)/modelfilemodule.o $(OBJECTDIR)/model_helper.o

all: geotherm

geotherm: $(SRC)/GeoTherm.f90 $(OBJECTDIR)/mathmodule.mod $(OBJECTDIR)/helpermodule.mod $(OBJECTDIR)/module_modelfile.mod $(OBJECTDIR)/model_helper.mod
	#mkdir -p bin
	gfortran $(LDFLAGS) $(FDFLAGS) -o $(OUT) $(SRC)/GeoTherm.f90 $(OBJECTS)
	mv $(OUT) $(BIN)

$(OBJECTDIR)/mathmodule.mod: $(SRC)/MathModule.f90
	gfortran -c $(SRC)/MathModule.f90
	mv mathmodule.mod mathmodule.o $(OBJECTDIR)

$(OBJECTDIR)/helpermodule.mod: $(SRC)/HelperModule.f90
	gfortran -c $(FDFLAGS) $(SRC)/HelperModule.f90
	mv helpermodule.mod helpermodule.o $(OBJECTDIR)

$(OBJECTDIR)/module_modelfile.mod: $(SRC)/ModelFileModule.f90
	gfortran -c $(SRC)/ModelFileModule.f90
	mv module_modelfile.mod modelfilemodule.o $(OBJECTDIR)

$(OBJECTDIR)/model_helper.mod: $(SRC)/model_helper.f90
	gfortran -c $(SRC)/model_helper.f90
	mv model_helper.mod model_helper.o $(OBJECTDIR)

#$(BIN)/: $(SRC)/MathModule.f90 $(SRC)/HelperModule.f90 $(SRC)/ModelFileModule.f90 $(SRC)/model_helper.f90
	#mkdir -p $(OBJECTDIR)

clean:
	rm -f $(OBJECTDIR)/mathmodule.o $(OBJECTDIR)/helpermodule.o $(OBJECTDIR)/modelfilemodule.o $(OBJECTDIR)/model_helper.o
	rm -f $(OBJECTDIR)/mathmodule.mod $(OBJECTDIR)/helpermodule.mod $(OBJECTDIR)/module_modelfile.mod $(OBJECTDIR)/model_helper.mod
	rm -f $(BIN)/geotherm
