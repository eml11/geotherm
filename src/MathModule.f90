      module MathModule 
      implicit none
!234567
      contains

      subroutine array_integral2d(data_ar,retrn_ar,incriment,n,m) 
      double precision :: data_ar(:,:)
      double precision :: retrn_ar(:,:)
      double precision :: trap(n,m)
      double precision :: midpnt(n,m) 
      double precision, dimension(n,m) :: fowrd_ar
      double precision, dimension(n,m) :: back_ar
      double precision :: incriment, sumvl
      integer :: i,j,n,m

      !using simsons!

      fowrd_ar(n,:) = data_ar(n,:)
      fowrd_ar(:n,:) = data_ar(2:,:)
  
      back_ar(1,:) = data_ar(1,:)
      back_ar(2:,:) = data_ar(:n,:)

      !this appears to work - but realy souldn't
      !check with other function at some point
      midpnt = -1*data_ar*incriment
      trap = -1*incriment*(fowrd_ar+back_ar)/2

      do j=1,m
        sumvl = 0d0
        do i=1,n
          sumvl = sumvl + (2d0/3d0)*midpnt(i,j) + (1d0/3d0)*trap(i,j)
          retrn_ar(i,j) = sumvl
        end do
      end do

      end subroutine

      subroutine array_diff2d(data_ar,retrn_ar,incriment,n,m)
      double precision :: data_ar(:,:)
      double precision :: retrn_ar(:,:)
      double precision :: incriment
      integer :: n,m
      double precision, dimension(n,m) :: fowrd_ar
      double precision, dimension(n,m) :: back_ar

      fowrd_ar(n,:) = data_ar(n,:)
      fowrd_ar(:n,:) = data_ar(2:,:)

      back_ar(1,:) = data_ar(1,:)
      back_ar(2:,:) = data_ar(:n,:)

      !is sign correct? seems strange
      retrn_ar = (back_ar-fowrd_ar)/(2.0*incriment)
      retrn_ar(1,:) = retrn_ar(2,:)
      retrn_ar(n,:) = retrn_ar(n-1,:)

      end subroutine

      subroutine array_generalintegral2d &
      &(data_ar,param_ar,retrn_ar,incriment,n,m)
      double precision :: data_ar(:,:)
      double precision :: param_ar(:,:)
      double precision :: retrn_ar(:,:)
      double precision :: incriment
      integer :: n,m
      double precision, dimension(n,m) :: param_dif
      double precision, dimension(n,m) :: integral_term  

      call array_diff2d(param_ar,param_dif,incriment,n,m)
      integral_term = data_ar*param_dif
      call array_integral2d(integral_term,retrn_ar,incriment,n,m)

      end subroutine

      subroutine array_generaldiff2d &
      &(data_ar,param_ar,retrn_ar,incriment,n,m)
      double precision :: data_ar(:,:)
      double precision :: param_ar(:,:)
      double precision :: retrn_ar(:,:)
      double precision :: incriment
      integer :: n,m
      double precision, dimension(n,m) :: data_dif
      double precision, dimension(n,m) :: param_dif
 
      call array_diff2d(param_ar,param_dif,incriment,n,m)
      call array_diff2d(data_ar,data_dif,incriment,n,m)
      retrn_ar = data_dif/param_dif

      end subroutine

      subroutine extend_ardimension(data_ar,retrn_ar,m)
      double precision :: data_ar(:)
      double precision :: retrn_ar(:,:)
      integer :: m, i

      do i=1,m
        retrn_ar(:,i) = data_ar
      end do

      end subroutine

      function gaussian_nnorm(t_ar,y_ar,sigma,n,m)
      double precision :: t_ar(n,m)
      double precision :: y_ar(n,m)
      double precision :: gaussian_nnorm(n,m)
      double precision :: sigma
      integer n,m

      gaussian_nnorm = DEXP((t_ar*t_ar+y_ar*y_ar)/(2*sigma))      

      end function

      end module
