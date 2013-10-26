      program GEOTHERM
      use helpermodule
      use mathmodule
      use module_modelfile
      implicit none

      type (modelfile) modelfile_inst
      character (len = 256) :: filename

      call GETARG(1,filename)

      call READMDLF(modelfile_inst,filename)
      call compute_geotherm &
      &(modelfile_inst,modelfile_inst%ydim,modelfile_inst%tdim)

      end program

      subroutine compute_geotherm &
      &(modelfile_inst,n,m)
      use equationpartsmodule
      use module_modelfile
      use mathmodule
      use helpermodule
      type (modelfile), intent(in) :: modelfile_inst
      double precision :: input_tdata_ar(m)
      double precision :: input_qdata_ar(m)
      double precision :: tdata_artrs(m,n), tdata_ar(n,m)
      double precision :: qdata_artrs(m,n), qdata_ar(n,m)
      double precision :: bdash_ar(n,m)
      double precision :: velocity_ar(n,m)
      double precision :: density_ar(n,m)
      double precision :: thermlconduct_ar(n,m)
      double precision :: heatproduct_ar(n,m)
      double precision :: heatcapc_ar(n,m)
      double precision :: exintegral_ar(n,m)
      double precision :: inerintegral_ar(n,m)
      double precision :: outerintegral_ar(n,m)
      double precision :: incriment(2)
      double precision :: iner_dbl
      double precision :: outr_dbl

      call get_netcdf1d(modelfile_inst%gtempfile,input_tdata_ar)
      call get_netcdf1d(modelfile_inst%gqfluxfile,input_qdata_ar)
      call get_netcdf_wincrmnt &
      &(modelfile_inst%velocitynetcdf,velocity_ar,incriment,n,m)
      call get_netcdf(modelfile_inst%densitynetcdf,density_ar)
      call get_netcdf(modelfile_inst%heatproductnetcdf,heatproduct_ar)
      call get_netcdf(modelfile_inst%heatcapcnetcdf,heatcapc_ar)
      call get_netcdf &
      &(modelfile_inst%thermlconductnetcdf,thermlconduct_ar)

      call extend_ardimension(input_tdata_ar,tdata_artrs,n)
      call extend_ardimension(input_qdata_ar,qdata_artrs,n)
      tdata_ar = RESHAPE(tdata_artrs,(/n,m/))
      qdata_ar = RESHAPE(qdata_artrs,(/n,m/))

      call compute_bdashval &
      &(tdata_ar,qdata_ar,thermlconduct_ar, &
      &bdash_ar,incriment,n,m)
      call compute_exponentintegral &
      &(bdash_ar,velocity_ar,thermlconduct_ar/(density_ar*heatcapc_ar),&
      &exintegral_ar,incriment,n,m)
      call compute_init_inerintegral &
      &(exintegral_ar,bdash_ar,(density_ar*heatproduct_ar)/heatcapc_ar,&
      &thermlconduct_ar/(density_ar*heatcapc_ar), &
      &inerintegral_ar,incriment,n,m)
      call compute_inerintegralconstant &
      &(inerintegral_ar,exintegral_ar,qdata_ar, &
      &thermlconduct_ar,iner_dbl,n,m)
      call compute_init_outerintegral &
      &(inerintegral_ar,exintegral_ar,iner_dbl, &
      &outerintegral_ar,incriment,n,m)
      call compute_outerintegralconstant &
      &(outerintegral_ar,tdata_ar,outr_dbl,n,m)

      call write_netcdf &
      &(modelfile_inst%outfile,outerintegral_ar+outr_dbl,n,m)

      end subroutine

