      program unit_tests
      use equationpartsmodule
      use mathmodule
      implicit none

      integer, parameter :: n = 100, m = 120
      double precision, dimension(n,m) :: t_ar,y_ar
      integer i,j
      double precision, dimension(n,m) :: wrk_ar1
      double precision, dimension(n,m) :: wrk_ar2
      double precision, dimension(n,m) :: wrk_ar3
      double precision, dimension(n,m) :: wrk_ar4
      double precision, dimension(n,m) :: wrk_ar5
      double precision, dimension(n,m) :: wrk_ar6

      double precision :: incriment(2) = (/1,1/)

      do i=1,n
        do j=1,m
          t_ar(i,j) = j
          y_ar(i,j) = i
        end do
      end do

      wrk_ar1 = gaussian_nnorm(t_ar,y_ar*0,2d0,n,m)
      wrk_ar2 = gaussian_nnorm(t_ar,y_ar*0,3d0,n,m)
      wrk_ar3 = gaussian_nnorm(t_ar,y_ar*0,5d0,n,m)
      wrk_ar4 = (t_ar*wrk_ar1*wrk_ar3)/(wrk_ar2*2)

      call compute_bdashval &
      &(wrk_ar1,wrk_ar2,wrk_ar3,wrk_ar5,incriment,n,m)

      print *,
      print *,
      print *, wrk_ar4-wrk_ar5

      wrk_ar2 = gaussian_nnorm(t_ar,y_ar,3d0,n,m)
      wrk_ar3 = gaussian_nnorm(t_ar,y_ar,5d0,n,m)
      wrk_ar4 = (1/(t_ar*(-0.5-1/3.0+0.2)))*(wrk_ar1*wrk_ar2)/wrk_ar3 +&
      &(1/(y_ar*(-1/3.0+0.2)))*(wrk_ar2/wrk_ar3) + &
      & (1/(t_ar*(-1+0.2)))*(wrk_ar2*wrk_ar2)/wrk_ar3 + &
      &(1/(y_ar*(0.2)))*wrk_ar2/wrk_ar3

      call compute_exponentintegral &
      &(wrk_ar1,wrk_ar2,wrk_ar3,wrk_ar5,incriment,n,m)

      print *,
      print *,
      print *, wrk_ar4-wrk_ar5

      wrk_ar4 = gaussian_nnorm(t_ar,y_ar,7d0,n,m)
      wrk_ar5 = (wrk_ar2*wrk_ar3*wrk_ar4)/(y_ar*(-1/3.0-0.2-1/7.0)) + &
      &(wrk_ar1*wrk_ar2*wrk_ar3*wrk_ar4)/(t_ar*(-1/3.0-0.2-1/7.0-0.5))

      call compute_init_inerintegral &
      &(wrk_ar2,wrk_ar1,wrk_ar3,wrk_ar4,wrk_ar6,incriment,n,m)

      print *,
      print *,
      print *, wrk_ar5-wrk_ar6

      end program
