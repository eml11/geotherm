SRC=$(PWD)/src
OBJECTDIR=$(PWD)/f90
BIN=$(PWD)/bin
TESTS=$(PWD)/tests
FC=gfortran
LDFLAGS=-L/usr/local/lib -lnetcdff -L/opt/local/lib -lnetcdf
FDFLAGS=-I/usr/local/include -I$(OBJECTDIR)
OUT=geotherm
OBJECTS=$(OBJECTDIR)/mathmodule.o $(OBJECTDIR)/helpermodule.o $(OBJECTDIR)/modelfilemodule.o $(OBJECTDIR)/EquationParts.o $(OBJECTDIR)/pressuresolver.o $(OBJECTDIR)/GeoChemSupprt.o
MODULES=$(OBJECTDIR)/equationpartsmodule.mod $(OBJECTDIR)/mathmodule.mod $(OBJECTDIR)/helpermodule.mod $(OBJECTDIR)/module_modelfile.mod $(OBJECTDIR)/pressuresolver.mod $(OBJECTDIR)/geochem.mod

all: geotherm domaingen unittests

geotherm: $(SRC)/GeoTherm.f90 $(MODULES)
	$(FC) $(LDFLAGS) $(FDFLAGS) -o $(OUT) $(SRC)/GeoTherm.f90 $(OBJECTS)
	mkdir -p $(BIN); mv $(OUT) $(BIN)

unittests: $(SRC)/unit_tests.f90 $(OBJECTDIR)/mathmodule.mod $(OBJECTDIR)/equationpartsmodule.mod $(OBJECTDIR)/module_modelfile.mod
	$(FC) -o unittests -I$(OBJECTDIR) $(SRC)/unit_tests.f90 $(OBJECTDIR)/mathmodule.o $(OBJECTDIR)/EquationParts.o $(OBJECTDIR)/modelfilemodule.o
	mkdir -p $(TESTS)/bin; mv unittests $(TESTS)/bin

domaingen: $(SRC)/generate_domain.f90
	$(FC) $(LDFLAGS) $(FDFLAGS) -o domaingen $(SRC)/generate_domain.f90
	mkdir -p $(BIN); mv domaingen $(BIN)

$(OBJECTDIR)/mathmodule.mod: $(SRC)/MathModule.f90 
	$(FC) -c $(SRC)/MathModule.f90
	mkdir -p $(OBJECTDIR); mv mathmodule.mod mathmodule.o $(OBJECTDIR)

$(OBJECTDIR)/helpermodule.mod: $(SRC)/HelperModule.f90 $(OBJECTDIR)/module_modelfile.mod $(OBJECTDIR)/pressuresolver.mod $(OBJECTDIR)/geochem.mod $(OBJECTDIR)/equationpartsmodule.mod
	$(FC) -c $(FDFLAGS) $(SRC)/HelperModule.f90 $(OBJECTDIR)/modelfilemodule.o $(OBJECTDIR)/EquationParts.o $(OBJECTDIR)/pressuresolver.o $(OBJECTDIR)/GeoChemSupprt.o
	mkdir -p $(OBJECTDIR); mv helpermodule.mod helpermodule.o $(OBJECTDIR)

$(OBJECTDIR)/module_modelfile.mod: $(SRC)/ModelFileModule.f90
	$(FC) -c $(SRC)/ModelFileModule.f90
	mkdir -p $(OBJECTDIR); mv module_modelfile.mod modelfilemodule.o $(OBJECTDIR)

$(OBJECTDIR)/equationpartsmodule.mod: $(SRC)/EquationParts.f90 $(OBJECTDIR)/mathmodule.mod $(OBJECTDIR)/module_modelfile.mod $(OBJECTDIR)/pressuresolver.mod
	$(FC) -c -I$(OBJECTDIR) -L$(OBJECTDIR) $(SRC)/EquationParts.f90 $(OBJECTDIR)/mathmodule.o $(OBJECTDIR)/modelfilemodule.o $(OBJECTDIR)/pressuresolver.o
	mkdir -p $(OBJECTDIR); mv equationpartsmodule.mod EquationParts.o $(OBJECTDIR)

$(OBJECTDIR)/pressuresolver.mod: $(SRC)/PressureSolver.f90 $(OBJECTDIR)/equationpartsmodule.mod
	$(FC) -c $(FDFLAGS) $(SRC)/PressureSolver.f90 $(OBJECTDIR)/mathmodule.o
	mkdir -p $(OBJECTDIR); mv pressuresolver.mod pressuresolver.o $(OBJECTDIR)

$(OBJECTDIR)/geochem.mod: $(SRC)/GeoChemSupprt.f90 $(OBJECTDIR)/module_modelfile.mod $(OBJECTDIR)/pressuresolver.mod $(OBJECTDIR)/equationpartsmodule.mod
	$(FC) -c $(FDFLAGS) $(SRC)/GeoChemSupprt.f90 $(OBJECTDIR)/modelfilemodule.o $(OBJECTDIR)/EquationParts.o $(OBJECTDIR)/pressuresolver.o
	mkdir -p $(OBJECTDIR); mv geochem.mod GeoChemSupprt.o $(OBJECTDIR)

test: clean geotherm
	cd $(TESTS); make clean
	cd $(TESTS); make test

clean:
	rm -f $(OBJECTS)
	rm -f $(MODULES)
	rm -f $(BIN)/geotherm $(BIN)/domaingen
	rmdir $(OBJECTDIR)
	rmdir $(BIN)
	cd $(TESTS); make clean
