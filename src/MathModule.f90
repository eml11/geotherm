      module MathModule 
      implicit none
!234567
      contains

      subroutine array_integral2d(data_ar,retrn_ar,incriment,n,m) 
      double precision :: data_ar(:,:)
      double precision :: retrn_ar(:,:)
      double precision :: incriment
      integer :: i,n,m

      do i=1,n
        retrn_ar(i,:) = SUM(data_ar(:i,:),1)!may be 0  
      end do
      
      retrn_ar = retrn_ar/incriment

      end subroutine

      subroutine array_diff2d(data_ar,retrn_ar,incriment,n,m)
      double precision :: data_ar(:,:)
      double precision :: retrn_ar(:,:)
      double precision :: incriment
      integer :: n,m
      double precision, dimension(n,m) :: fowrd_ar
      double precision, dimension(n,m) :: back_ar

      fowrd_ar(1,:) = data_ar(1,:)
      fowrd_ar(1:,:) = data_ar(:n,:)

      back_ar(n,:) = data_ar(n,:)!may not work
      back_ar(:-1,:) = data_ar(1:,:)

      retrn_ar = (fowrd_ar-back_ar)/(2.0*incriment)

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

      end module
