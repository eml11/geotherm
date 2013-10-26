SRCNETCDF=$(PWD)/src_netcdf
SRCBNDFL=$(PWD)/src_boundryfiles
GEOTHERM=$(PWD)/../bin/geotherm

test:

test1: $(GEOTHERM) constant_heatproduct_tindipendent.txt
	$(GEOTHERM) constant_heatproduct_tindipendent.txt 

clean:
	rmdir results
