SRCNETCDF=$(PWD)/src_netcdf
SRCBNDFL=$(PWD)/src_boundryfiles
RESULTS=$(PWD)/results
GEOTHERM=$(PWD)/../bin/geotherm
UNITTEST=$(PWD)/bin/unittests
BIN=$(PWD)/bin

test: test1 test2 test3 unittest testvelocity

unittest:
	$(unittest)

test1: $(GEOTHERM) constant_heatproduct_tindipendent.txt
	$(GEOTHERM) constant_heatproduct_tindipendent.txt 

test2: $(GEOTHERM) constant_heatproducthigh_tindipendent.txt
	$(GEOTHERM) constant_heatproducthigh_tindipendent.txt

test3: $(GEOTHERM) exponent_heatproduct_tindipendent.txt
	$(GEOTHERM) exponent_heatproduct_tindipendent.txt

testvelocity: constant_heatproduct_withvelocity.txt
	$(GEOTHERM) constant_heatproduct_withvelocity.txt

clean:
	rm -f $(RESULTS)/*.nc
	#rm -f $(UNITTEST)
	#rmdir $(BIN)
	#rmdir results
