      program unit_tests
      use equationpartsmodule
      use mathmodule
      implicit none

      integer, parameter :: n = 1000, m = 1200
      double precision, dimension(n,m) :: t_ar,y_ar
      integer i,j
      double precision, dimension(n,m) :: wrk_ar1
      double precision, dimension(n,m) :: wrk_ar2
      double precision, dimension(n,m) :: wrk_ar3
      double precision, dimension(n,m) :: wrk_ar4
      double precision, dimension(n,m) :: wrk_ar5
      double precision, dimension(n,m) :: wrk_ar6

      double precision :: incriment(2) = (/0.001,0.001/)

      do i=1,n
        do j=1,m
          t_ar(i,j) = i/1000.
          y_ar(i,j) = j/1000.
        end do
      end do
    
      wrk_ar1 = gaussian_nnorm(t_ar,y_ar,2d0,n,m)
      wrk_ar2 = -2*wrk_ar1/t_ar

      call array_integral2d(wrk_ar1,wrk_ar3,incriment(1),n,m)

      print *,
      print *,
      print *, wrk_ar3(3,1)+wrk_ar2(1,1)-wrk_ar3(1,1)," ", wrk_ar2(3,1)
      print *, MAXVAL( &
      &ABS((wrk_ar2(2:n,:)-wrk_ar3(2:n,:)+wrk_ar2(1,1)-wrk_ar3(1,1))/ &
      &wrk_ar2(2:n,:)))
      print *, MINVAL( &
      &ABS((wrk_ar2(2:n,:)-wrk_ar3(2:n,:)+wrk_ar2(1,1)-wrk_ar3(1,1))/ &
      &wrk_ar2(2:n,:)))

      wrk_ar2 = -t_ar*wrk_ar1/2

      call array_diff2d(wrk_ar1,wrk_ar3,incriment(1),n,m)

      print *,
      print *,
      print *, wrk_ar3(3,1)," ", wrk_ar2(3,1)
      print *, MAXVAL( &
      &ABS((wrk_ar2(2:n,:)-wrk_ar3(2:n,:))/wrk_ar2(2:n,:)))
      print *, MINVAL( &
      &ABS((wrk_ar2(2:n,:)-wrk_ar3(2:n,:))/wrk_ar2(2:n,:)))

      wrk_ar1 = gaussian_nnorm(t_ar,y_ar*0,2d0,n,m)
      wrk_ar2 = gaussian_nnorm(t_ar,y_ar*0,3d0,n,m)
      wrk_ar3 = gaussian_nnorm(t_ar,y_ar*0,5d0,n,m)
      wrk_ar4 = (t_ar*wrk_ar1*wrk_ar3)/(wrk_ar2*2)

      call compute_bdashval &
      &(wrk_ar1,wrk_ar2,wrk_ar3,wrk_ar5,incriment,n,m)

      print *,
      print *,
      print *, wrk_ar5(1,1)," ", wrk_ar4(1,1)
      print *, MAXVAL( &
      &ABS((wrk_ar4(2:n,:)-wrk_ar5(2:n,:))/wrk_ar4(2:n,:)))
      print *, MINVAL( &
      &ABS((wrk_ar4(2:n,:)-wrk_ar5(2:n,:))/wrk_ar4(2:n,:)))

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
      print *, wrk_ar5(1,1)," ", wrk_ar4(1,1)
      print *, MAXVAL( &
      &ABS((wrk_ar4(2:n,:)-wrk_ar5(2:n,:))/wrk_ar4(2:n,:)))
      print *, MINVAL( &
      &ABS((wrk_ar4(2:n,:)-wrk_ar5(2:n,:))/wrk_ar4(2:n,:)))
 
      wrk_ar4 = gaussian_nnorm(t_ar,y_ar,7d0,n,m)
      wrk_ar5 = (wrk_ar2*wrk_ar3*wrk_ar4)/(y_ar*(-1/3.0-0.2-1/7.0)) + &
      &(wrk_ar1*wrk_ar2*wrk_ar3*wrk_ar4)/(t_ar*(-1/3.0-0.2-1/7.0-0.5))

      call compute_init_inerintegral &
      &(wrk_ar2,wrk_ar1,wrk_ar3,wrk_ar4,wrk_ar6,incriment,n,m)

      print *,
      print *,
      print *, wrk_ar6(1,1)," ", wrk_ar5(1,1)
      print *, MAXVAL( &
      &ABS((wrk_ar5(2:n,:)-wrk_ar6(2:n,:))/wrk_ar5(2:n,:)))
      print *, MINVAL( &
      &ABS((wrk_ar5(2:n,:)-wrk_ar6(2:n,:))/wrk_ar5(2:n,:)))

      end program
