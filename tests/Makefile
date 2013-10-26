SRCNETCDF=$(PWD)/src_netcdf
SRCBNDFL=$(PWD)/src_boundryfiles
RESULTS=$(PWD)/results
GEOTHERM=$(PWD)/../bin/geotherm

test: test1

test1: $(GEOTHERM) constant_heatproduct_tindipendent.txt
	$(GEOTHERM) constant_heatproduct_tindipendent.txt 

clean:
	rm -f $(RESULTS)/*.nc
	#rmdir results
