      module pressuresolver
      use mathmodule
      implicit none
      contains

      subroutine compute_pressure &
      &(density_ar,incompresibility_ar,retrn_ar,incriment,n,m)
      double precision :: density_ar(n,m)
      double precision :: incompresibility_ar(n,m)
      double precision :: retrn_ar(n,m)
      double precision :: gravity_const
      double precision :: y_integral(n,m)
      double precision :: incriment(2)
      integer n,m

      call array_integral2dydim(density_ar*gravity_const, & 
      &y_integral,incriment(2),n,m)

      retrn_ar = -incompresibility_ar * &
      &DLOG(1-y_integral*incompresibility_ar)
      end subroutine

      subroutine compute_pddensity &
      &(density_ar,incompresibility_ar,retrn_ar,incriment,n,m)
      double precision :: density_ar(n,m)
      double precision :: incompresibility_ar(n,m)
      double precision :: retrn_ar(n,m)
      double precision :: gravity_const
      double precision :: y_integral(n,m)
      double precision :: incriment(2)
      integer n,m

      call array_integral2dydim(density_ar*gravity_const, & 
      &y_integral,incriment(2),n,m)

      retrn_ar = (density_ar)/(1-y_integral*incompresibility_ar)

      end subroutine

      end module
