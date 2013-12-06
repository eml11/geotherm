SRC=$(PWD)/src
OBJECTDIR=$(PWD)/f90
BIN=$(PWD)/bin
TESTS=$(PWD)/tests
FC=gfortran
LDFLAGS=-L/usr/local/lib -lnetcdff -L/opt/local/lib -lnetcdf
FDFLAGS=-I/usr/local/include -I$(OBJECTDIR)
OUT=geotherm
OBJECTS=$(OBJECTDIR)/mathmodule.o $(OBJECTDIR)/helpermodule.o $(OBJECTDIR)/modelfilemodule.o $(OBJECTDIR)/EquationParts.o $(OBJECTDIR)/pressuresolver.o $(OBJECTDIR)/GeoChemSupprt.o $(OBJECTDIR)/modeldomainmodule.o $(OBJECTDIR)/modeloutput.o $(OBJECTDIR)/modelregionmodule.o $(OBJECTDIR)/modellogfile.o
MODULES=$(OBJECTDIR)/equationpartsmodule.mod $(OBJECTDIR)/mathmodule.mod $(OBJECTDIR)/helpermodule.mod $(OBJECTDIR)/module_modelfile.mod $(OBJECTDIR)/pressuresolver.mod $(OBJECTDIR)/geochem.mod $(OBJECTDIR)/modeldomainmodule.mod $(OBJECTDIR)/modeloutput.mod $(OBJECTDIR)/modelregionmodule.mod $(OBJECTDIR)/modellogfile.mod

all: geotherm domaingen unittests isostatics slicenetcdf

geotherm: $(SRC)/GeoTherm.f90 $(MODULES)
	$(FC) $(LDFLAGS) $(FDFLAGS) -o $(OUT) $(SRC)/GeoTherm.f90 $(OBJECTS)
	mkdir -p $(BIN); mv $(OUT) $(BIN)

unittests: $(SRC)/unit_tests.f90 $(MODULES)
	$(FC) $(LDFLAGS) $(FDFLAGS) -o unittests $(SRC)/unit_tests.f90 $(OBJECTS)
	mkdir -p $(TESTS)/bin; mv unittests $(TESTS)/bin

isostatics: $(SRC)/Isostatics.f90 $(MODULES) $(OBJECTDIR)/isostatichelper.mod
	$(FC) $(LDFLAGS) $(FDFLAGS) -o isostatics $(SRC)/Isostatics.f90 $(OBJECTS) $(OBJECTDIR)/isostatichelper.o
	mkdir -p $(BIN); mv isostatics $(BIN)

slicenetcdf: $(SRC)/SliceNetcdf.f90 $(OBJECTDIR)/helpermodule.mod $(OBJECTDIR)/mathmodule.mod
	$(FC) $(LDFLAGS) $(FDFLAGS) -o slicenetcdf $(SRC)/SliceNetcdf.f90 $(OBJECTDIR)/helpermodule.o $(OBJECTDIR)/mathmodule.o
	mkdir -p $(BIN); mv slicenetcdf $(BIN)

$(OBJECTDIR)/isostatichelper.mod: $(SRC)/isostatichelper.f90 $(OBJECTDIR)/mathmodule.mod $(OBJECTDIR)/modeldomainmodule.mod
	$(FC) -c $(LDFLAGS) $(FDFLAGS) $(SRC)/isostatichelper.f90
	mkdir -p $(OBJECTDIR); mv isostatichelper.mod isostatichelper.o $(OBJECTDIR)

domaingen: $(SRC)/generate_domain.f90
	$(FC) $(LDFLAGS) $(FDFLAGS) -o domaingen $(SRC)/generate_domain.f90
	mkdir -p $(BIN); mv domaingen $(BIN)

$(OBJECTDIR)/helpermodule.mod: $(SRC)/HelperModule.f90
	$(FC) -c $(FDFLAGS) -L$(OBJECTDIR) $(SRC)/HelperModule.f90
	mkdir -p $(OBJECTDIR); mv helpermodule.mod helpermodule.o $(OBJECTDIR)

$(OBJECTDIR)/module_modelfile.mod: $(SRC)/ModelFileModule.f90 $(OBJECTDIR)/pressuresolver.mod $(OBJECTDIR)/geochem.mod $(OBJECTDIR)/modeldomainmodule.mod $(OBJECTDIR)/helpermodule.mod $(OBJECTDIR)/modellogfile.mod
	$(FC) -c $(FDFLAGS) -L$(OBJECTDIR) $(SRC)/ModelFileModule.f90 $(OBJECTDIR)/pressuresolver.o $(OBJECTDIR)/GeoChemSupprt.o $(OBJECTDIR)/modeldomainmodule.o $(OBJECTDIR)/helpermodule.o $(OBJECTDIR)/modellogfile.o
	mkdir -p $(OBJECTDIR); mv module_modelfile.mod modelfilemodule.o $(OBJECTDIR) 

$(OBJECTDIR)/equationpartsmodule.mod: $(SRC)/EquationParts.f90 $(OBJECTDIR)/mathmodule.mod $(OBJECTDIR)/pressuresolver.mod $(OBJECTDIR)/modeldomainmodule.mod
	$(FC) -c -I$(OBJECTDIR) -L$(OBJECTDIR) $(SRC)/EquationParts.f90 $(OBJECTDIR)/mathmodule.o $(OBJECTDIR)/pressuresolver.o $(OBJECTDIR)/modeldomainmodule.o
	mkdir -p $(OBJECTDIR); mv equationpartsmodule.mod EquationParts.o $(OBJECTDIR)

$(OBJECTDIR)/pressuresolver.mod: $(SRC)/PressureSolver.f90 $(OBJECTDIR)/mathmodule.mod $(OBJECTDIR)/modeldomainmodule.mod
	$(FC) -c $(FDFLAGS) $(SRC)/PressureSolver.f90 $(OBJECTDIR)/mathmodule.o $(OBJECTDIR)/modeldomainmodule.o
	mkdir -p $(OBJECTDIR); mv pressuresolver.mod pressuresolver.o $(OBJECTDIR)

$(OBJECTDIR)/geochem.mod: $(SRC)/GeoChemSupprt.f90
	$(FC) -c $(FDFLAGS) $(SRC)/GeoChemSupprt.f90
	mkdir -p $(OBJECTDIR); mv geochem.mod GeoChemSupprt.o $(OBJECTDIR)

$(OBJECTDIR)/modeldomainmodule.mod: $(SRC)/modeldomainmodule.f90 $(OBJECTDIR)/geochem.mod $(OBJECTDIR)/modelregionmodule.mod
	$(FC) -c $(FDFLAGS) -L$(OBJECTDIR) $(SRC)/modeldomainmodule.f90 $(OBJECTDIR)/GeoChemSupprt.o $(OBJECTDIR)/modelregionmodule.o
	mkdir -p $(OBJECTDIR); mv modeldomainmodule.o modeldomainmodule.mod $(OBJECTDIR)

$(OBJECTDIR)/modelregionmodule.mod: $(SRC)/modelregionmodule.f90 $(OBJECTDIR)/geochem.mod 
	$(FC) -c $(FDFLAGS) -L$(OBJECTDIR) $(SRC)/modelregionmodule.f90 $(OBJECTDIR)/GeoChemSupprt.o
	mkdir -p $(OBJECTDIR); mv modelregionmodule.o modelregionmodule.mod $(OBJECTDIR)

$(OBJECTDIR)/modeloutput.mod: $(SRC)/modeloutput.f90  $(OBJECTDIR)/module_modelfile.mod $(OBJECTDIR)/pressuresolver.mod $(OBJECTDIR)/geochem.mod $(OBJECTDIR)/equationpartsmodule.mod
	$(FC) -c $(FDFLAGS) -L$(OBJECTDIR) $(SRC)/modeloutput.f90 $(OBJECTDIR)/modelfilemodule.o $(OBJECTDIR)/EquationParts.o $(OBJECTDIR)/pressuresolver.o $(OBJECTDIR)/GeoChemSupprt.o
	mkdir -p $(OBJECTDIR); mv modeloutput.mod modeloutput.o $(OBJECTDIR)

$(OBJECTDIR)/modellogfile.mod: $(SRC)/modellogfile.f90
	$(FC) -c $(SRC)/modellogfile.f90
	mkdir -p $(OBJECTDIR); mv modellogfile.mod modellogfile.o $(OBJECTDIR)

$(OBJECTDIR)/mathmodule.mod: $(SRC)/MathModule.f90
	$(FC) -c $(SRC)/MathModule.f90
	mkdir -p $(OBJECTDIR); mv mathmodule.mod mathmodule.o $(OBJECTDIR)

test: geotherm
	cd $(TESTS); make clean
	cd $(TESTS); make test

clean:
	rm -f $(OBJECTS) $(OBJECTDIR)/isostatichelper.mod $(OBJECTDIR)/isostatichelper.o
	rm -f $(MODULES)
	rm -f $(BIN)/geotherm $(BIN)/domaingen
	rmdir $(OBJECTDIR)
	rmdir $(BIN)
	cd $(TESTS); make clean
