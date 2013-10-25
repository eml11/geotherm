      program GEOTHERM
      use helpermodule
      use mathmodule
      use module_modelfile
      use model_helper
      implicit none
!234567
!read out n and m first then pass to a subroutine to read in the arrays
      type (physical_constants) physical_constants_inst
      type (model_constants) model_constants_inst
      type (modelfile) modelfile_inst
      character (len = 256) :: filename
      integer :: n = 1
      integer :: m = 1!temp until file reading is complete

      call READMDLF(modelfile_inst,filename)
      call compute_geotherm(modelfile_inst,n,m)

      end program

      subroutine compute_geotherm &
      &(modelfile_inst,n,m)
      use module_modelfile
      use model_helper
      use mathmodule
      use helpermodule
      type (modelfile), intent(in) :: modelfile_inst
      type (physical_constants) :: physical_constants_inst
      double precision :: input_tdata_ar(m)
      double precision :: input_qdata_ar(m)
      double precision :: tdata_artrs(m,n), tdata_ar(n,m)
      double precision :: qdata_artrs(m,n), qdata_ar(n,m)
      double precision :: bdash_ar(n,m)
      double precision :: velocity_ar(n,m)
      double precision :: kappa_ar(n,m)
      double precision :: heatproduct_ar(n,m)
      double precision :: heatcapc_ar(n,m)
      double precision :: exintegral_ar(n,m)
      double precision :: inerintegral_ar(n,m)
      double precision :: outerintegral_ar(n,m)
      double precision :: iner_dbl
      double precision :: outr_dbl

      call get_file(modelfile_inst%gtempfile,input_tdata_ar)
      call get_file(modelfile_inst%gqfluxfile,input_qdata_ar)
      call get_netcdf(modelfile_inst%velocitynetcdf,velocity_ar)
      call get_netcdf(modelfile_inst%kappanetcdf,kappa_ar)
      call get_netcdf(modelfile_inst%heatproductnetcdf,heatproduct_ar)
      call get_netcdf(modelfile_inst%heatcapcnetcdf,heatcapc_ar)

      call extend_ardimension(input_tdata_ar,tdata_artrs,n)
      call extend_ardimension(input_qdata_ar,qdata_artrs,n)
      tdata_ar = RESHAPE(tdata_artrs,(/n,m/))
      qdata_ar = RESHAPE(qdata_artrs,(/n,m/))

      call compute_bdashval &
      &(tdata_ar,qdata_ar,physical_constants_inst%kconstant, &
      &bdash_ar,n,m)
      call compute_exponentintegral &
      &(bdash_ar,velocity_ar,kappa_ar,exintegral_ar,n,m)
      call compute_init_inerintegral &
      &(exintegral_ar,bdash_ar,heatproduct_ar/heatcapc_ar,kappa_ar, &
      &inerintegral_ar,n,m)
      call compute_inerintegralconstant &
      &(inerintegral_ar,exintegral_ar,qdata_ar, &
      &physical_constants_inst%kconstant,iner_dbl,n,m)
      call compute_init_outerintegral &
      &(inerintegral_ar,exintegral_ar,iner_dbl,outerintegral_ar,n,m)
      call compute_outerintegralconstant &
      &(outerintegral_ar,tdata_ar,outr_dbl,n,m)

      end subroutine

      subroutine compute_bdashval &
      &(tdata_ar,qdata_ar,kconstant,retrn_ar,n,m)
      use mathmodule
      double precision :: tdata_ar(n,m)
      double precision :: qdata_ar(n,m)
      double precision :: kconstant
      double precision :: retrn_ar(n,m)
      integer :: n, m
      double precision, dimension(m,n) :: difftdata_ar

      
      call array_diff2d(RESHAPE(tdata_ar,(/m,n/)),difftdata_ar,n,m)
      retrn_ar = kconstant * &
      &(RESHAPE(difftdata_ar,(/n,m/))/qdata_ar)

      end subroutine

      subroutine compute_exponentintegral &
      &(bdash_ar,velocity_ar,kappa_ar,retrn_ar,n,m)
      use mathmodule
      double precision :: bdash_ar(n,m)
      double precision :: velocity_ar(n,m)
      double precision :: kappa_ar(n,m)
      double precision :: retrn_ar(n,m)
      integer n,m
      double precision, dimension(m,n) :: v_tintegral, b_tintegral
      double precision, dimension(n,m) :: v_yintegral, b_yintegral

      call array_integral2d & 
      &(RESHAPE((bdash_ar*velocity_ar)/kappa_ar,(/m,n/)),v_tintegral)
      call array_integral2d &
      (RESHAPE((bdash_ar*bdash_ar)/kappa_ar,(/m,n/)),b_tintegral)
      call array_integral2d(velocity_ar/kappa_ar,v_yintegral)
      call array_integral2d(1/kappa_ar,b_yintegral)
 
      retrn_ar = RESHAPE(v_tintegral,(/n,m/)) + &
      &RESHAPE(b_tintegral,(/n,m/)) + &
      &v_yintegral - bdash_ar*b_yintegral

      end subroutine

      subroutine compute_init_inerintegral &
      &(exintegral_ar,bdash_ar,theta_ar,kappa_ar,retrn_ar,n,m)
      use mathmodule
      double precision :: exintegral_ar(n,m)
      double precision :: bdash_ar(n,m)
      double precision :: retrn_ar(n,m)
      double precision :: theta_ar(n,m)
      double precision :: kappa_ar(n,m)
      double precision, dimension(n,m) :: integral_term, y_integral
      double precision, dimension(m,n) :: t_integral
     
      integral_term = theta_ar*kappa_ar*DEXP(-1*exintegral_ar)
      call array_integral2d(integral_term,y_integral)
      call array_integral2d &
      &(RESHAPE(integral_term*bdash_ar,(/m,n/)),t_integral)

      retrn_ar = y_integral + RESHAPE(t_integral,(/n,m/))

      end subroutine

      subroutine compute_inerintegralconstant &
      &(inerintegral_ar,exintegral_ar,qdata_ar,kconstant,retrn_dbl,n,m)
      double precision :: inerintegral_ar(n,m)
      double precision :: exintegral_ar(n,m)
      double precision :: qdata_ar(n,m)
      double precision :: kconstant
      double precision :: retrn_dbl

      retrn_dbl = &
      &((-1*qdata_ar(1,1)*qdata_ar(1,1))/kconstant) * &
      &EXP(exintegral_ar(1,1)) - inerintegral_ar(1,1)

      end subroutine

      subroutine compute_init_outerintegral &
      &(inerintegral_ar,exintegral_ar,iner_dbl,retrn_ar,n,m)
      use mathmodule
      double precision :: inerintegral_ar(n,m)
      double precision :: exintegral_ar(n,m)
      double precision :: retrn_ar(n,m)
      double precision :: iner_dbl
      double precision, dimension(n,m) :: integral_term, y_integral
      double precision, dimension(m,n) :: t_integral

      integral_term = EXP(exintegral_ar)*(inerintegral_ar + iner_dbl)
      call array_integral2d(integral_term,y_integral)
      call array_integral2d &
      &(RESHAPE(integral_term*bdash_ar,(/m,n/)),t_integral)

      retrn_ar = y_integral + RESHAPE(t_integral,(/n,m/))

      end subroutine

      subroutine compute_outerintegralconstant &
      &(outerintegral_ar,tdata_ar,retrn_dbl,n,m)
      double precision :: outerintegral_ar(n,m)
      double precision :: tdata_ar(n,m)
      double precision retrn_dbl

      retrn_dbl = tdata_ar(1,1) - outerintegral_ar(1,1)

      end subroutine
