SRC=$(PWD)/src
OBJECTDIR=$(PWD)/f90
BIN=$(PWD)/bin
TESTS=$(PWD)/tests
FC=gfortran
LDFLAGS=-L/distributions/netcdf-fortran-4.2/fortran/.libs -L/opt/local/lib -lnetcdf -lnetcdff
FDFLAGS=-I/distributions/netcdf-fortran-4.2/f90 -I$(OBJECTDIR)
OUT=geotherm
OBJECTS=$(OBJECTDIR)/mathmodule.o $(OBJECTDIR)/helpermodule.o $(OBJECTDIR)/modelfilemodule.o $(OBJECTDIR)/model_helper.o

all: geotherm

geotherm: $(SRC)/GeoTherm.f90 $(OBJECTDIR)/mathmodule.mod $(OBJECTDIR)/helpermodule.mod $(OBJECTDIR)/module_modelfile.mod $(OBJECTDIR)/model_helper.mod
	$(FC) $(LDFLAGS) $(FDFLAGS) -o $(OUT) $(SRC)/GeoTherm.f90 $(OBJECTS)
	mkdir -p $(BIN); mv $(OUT) $_

$(OBJECTDIR)/mathmodule.mod: $(SRC)/MathModule.f90
	$(FC) -c $(SRC)/MathModule.f90
	mkdir -p $(OBJECTDIR); mv mathmodule.mod mathmodule.o $_

$(OBJECTDIR)/helpermodule.mod: $(SRC)/HelperModule.f90
	$(FC) -c $(FDFLAGS) $(SRC)/HelperModule.f90
	mkdir -p $(OBJECTDIR); mv helpermodule.mod helpermodule.o $_ 

$(OBJECTDIR)/module_modelfile.mod: $(SRC)/ModelFileModule.f90
	$(FC) -c $(SRC)/ModelFileModule.f90
	mkdir -p $(OBJECTDIR); mv module_modelfile.mod modelfilemodule.o $_

$(OBJECTDIR)/model_helper.mod: $(SRC)/model_helper.f90
	$(FC) -c $(SRC)/model_helper.f90
	mkdir -p $(OBJECTDIR); mv model_helper.mod model_helper.o $_

test:
	cd $(TESTS); make test

clean:
	rm -f $(OBJECTDIR)/mathmodule.o $(OBJECTDIR)/helpermodule.o $(OBJECTDIR)/modelfilemodule.o $(OBJECTDIR)/model_helper.o
	rm -f $(OBJECTDIR)/mathmodule.mod $(OBJECTDIR)/helpermodule.mod $(OBJECTDIR)/module_modelfile.mod $(OBJECTDIR)/model_helper.mod
	rm -f $(BIN)/geotherm
	rmdir $(OBJECTDIR)
	rmdir $(BIN)
