      program reg_tests
      implicit none
      use equationpartsmodule
      use mathmodule 

      integer, parameter :: n = 100, m = 120
      double precision, dimension(n,m) :: t_ar,y_tar
      integer i,j
      double precision, dimension(n,m) :: wrk_ar1,wrk_ar2,wrk_ar3 &
      &wrk_ar4,wrk_ar5,wrk_ar6
      double precision :: incriment(2) = (/1,1/)

      do i=1,n
        do j=1,m
          t_ar(i,j) = j
          y_ar(i,j) = i
        end do
      end do

      wrk_ar1 = gaussian_nnorm(t_ar,y_ar*0,2,n,m)
      wrk_ar2 = gaussian_nnorm(t_ar,y_ar*0,3,n,m)
      wrk_ar3 = gaussian_nnorm(t_ar,y_ar*0,5,n,m)
      wrk_ar4 = (t_ar*wrk_ar1*wrk_ar3)/(wrk_ar2*2)
    
      call compute_bdashval &
      &(wrk_ar1,wrk_ar2,wrk_ar3,wrk_ar5,incriment,n,m)

      wrk_ar2 = gaussian_nnorm(t_ar,y_ar,3,n,m)
      wrk_ar3 = gaussian_nnorm(t_ar,y_ar,5,n,m)
      wrk_ar4 = (1/(t_ar*(-0.5-1/3.0+0.2)))*(wrk_ar1*wrk_ar2)/wrk_ar3 +&
      &(1/(y_ar*(-1/3.0+0.2)))*(wrk_ar2/wrk_ar3) + &
      & (1/(t_ar*(-1+0.2)))*(wrk_ar2*wrk_ar2)/wrk_ar3 + &
      &(1/(y_ar*(0.2)))*wrk_ar2/wrk_ar3
 
      call compute_exponentintegral &
      &(wrk_ar1,wrk_ar2,wrk_ar3,wrk_ar5,incriment,n,m)

      

      end program

      module equationpartsmodule

      contains

      subroutine compute_bdashval &
      &(tdata_ar,qdata_ar,kconstant,retrn_ar,incriment,n,m)
      use mathmodule
      double precision :: tdata_ar(n,m)
      double precision :: qdata_ar(n,m)
      double precision :: kconstant(n,m)
      double precision :: retrn_ar(n,m)
      double precision :: incriment(2)
      integer :: n, m
      double precision, dimension(m,n) :: difftdata_ar


      call array_diff2d(RESHAPE(tdata_ar,(/m,n/)), &
      &difftdata_ar,incriment(1),n,m)
      retrn_ar = kconstant * &
      &(RESHAPE(difftdata_ar,(/n,m/))/qdata_ar)

      end subroutine

      subroutine compute_exponentintegral &
      &(bdash_ar,velocity_ar,kappa_ar,retrn_ar,incriment,n,m)
      use mathmodule
      double precision :: bdash_ar(n,m)
      double precision :: velocity_ar(n,m)
      double precision :: kappa_ar(n,m)
      double precision :: retrn_ar(n,m)
      double precision :: incriment(2)
      integer n,m
      double precision, dimension(m,n) :: v_tintegral, b_tintegral
      double precision, dimension(n,m) :: v_yintegral, b_yintegral

      call array_integral2d &
      &(RESHAPE((bdash_ar*velocity_ar)/kappa_ar,(/m,n/)), &
      &v_tintegral,incriment(1),n,m)
      call array_integral2d &
      (RESHAPE((bdash_ar*bdash_ar)/kappa_ar,(/m,n/)), &
      &b_tintegral,incriment(1),n,m)
      call array_integral2d(velocity_ar/kappa_ar, &
      &v_yintegral,incriment(2),n,m)
      call array_integral2d(1/kappa_ar,b_yintegral,incriment(2),n,m)

      retrn_ar = RESHAPE(v_tintegral,(/n,m/)) + &
      &RESHAPE(b_tintegral,(/n,m/)) + &
      &v_yintegral - bdash_ar*b_yintegral

      end subroutine

      subroutine compute_init_inerintegral &
      &(exintegral_ar,bdash_ar,theta_ar,kappa_ar,retrn_ar,incriment,n,m)
      use mathmodule
      double precision :: exintegral_ar(n,m)
      double precision :: bdash_ar(n,m)
      double precision :: retrn_ar(n,m)
      double precision :: theta_ar(n,m)
      double precision :: kappa_ar(n,m)
      double precision :: incriment(2)
      double precision, dimension(n,m) :: integral_term, y_integral
      double precision, dimension(m,n) :: t_integral

      integral_term = theta_ar*kappa_ar*DEXP(-1*exintegral_ar)
      call array_integral2d(integral_term,y_integral,incriment(2),n,m)
      call array_integral2d &
      &(RESHAPE(integral_term*bdash_ar,(/m,n/)), &
      &t_integral,incriment(1),n,m)

      retrn_ar = y_integral + RESHAPE(t_integral,(/n,m/))

      end subroutine

      subroutine compute_inerintegralconstant &
      &(inerintegral_ar,exintegral_ar,qdata_ar,kconstant,retrn_dbl,n,m)
      double precision :: inerintegral_ar(n,m)
      double precision :: exintegral_ar(n,m)
      double precision :: qdata_ar(n,m)
      double precision :: kconstant(n,m)
      double precision :: retrn_dbl

      retrn_dbl = &
      &((-1*qdata_ar(1,1)*qdata_ar(1,1))/kconstant(1,1)) * &
      &EXP(exintegral_ar(1,1)) - inerintegral_ar(1,1)

      end subroutine

      subroutine compute_init_outerintegral &
      &(inerintegral_ar,exintegral_ar,iner_dbl,retrn_ar,incriment,n,m)
      use mathmodule
      double precision :: inerintegral_ar(n,m)
      double precision :: exintegral_ar(n,m)
      double precision :: retrn_ar(n,m)
      double precision :: iner_dbl
      double precision :: incriment(2)
      double precision, dimension(n,m) :: integral_term, y_integral
      double precision, dimension(m,n) :: t_integral

      integral_term = EXP(exintegral_ar)*(inerintegral_ar + iner_dbl)
      call array_integral2d(integral_term,y_integral,incriment(2),n,m)
      call array_integral2d &
      &(RESHAPE(integral_term*bdash_ar,(/m,n/)), &
      &t_integral,incriment(1),n,m)

      retrn_ar = y_integral + RESHAPE(t_integral,(/n,m/))

      end subroutine

      subroutine compute_outerintegralconstant &
      &(outerintegral_ar,tdata_ar,retrn_dbl,n,m)
      double precision :: outerintegral_ar(n,m)
      double precision :: tdata_ar(n,m)
      double precision retrn_dbl

      retrn_dbl = tdata_ar(1,1) - outerintegral_ar(1,1)

      end subroutine

      end module
