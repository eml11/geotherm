SRCNETCDF=$(PWD)/src_netcdf
SRCBNDFL=$(PWD)/src_boundryfiles
RESULTS=$(PWD)/results
GEOTHERM=$(PWD)/../bin/geotherm
UNITTEST=$(PWD)/bin/unittests

test: test1 unittest

unittest:
	$(unittest)

test1: $(GEOTHERM) constant_heatproduct_tindipendent.txt
	$(GEOTHERM) constant_heatproduct_tindipendent.txt 

clean:
	rm -f $(RESULTS)/*.nc
	#rmdir results
