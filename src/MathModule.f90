!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  
!  geotherm.
!  
!  Copyright (C) 2013 Elliot Lynch.
!  
!  Elliot M. Lynch, eml11@ic.ac.uk
!  Department of Earth Science and Engineering
!  Imperial College London
!  
!  This library is free software; you can redistribute it and/or
!  modify it under the terms of the GNU Lesser General Public
!  License as published by the Free Software Foundation,
!  version 2.1 of the License.
!  
!  This library is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!  Lesser General Public License for more details.
!  
!  You should have received a copy of the GNU Lesser General Public
!  License along with this library; if not, write to the Free Software
!  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
!  USA
!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      module MathModule 
      implicit none

      contains

      !> integrates data_ar along dimension 1
      !! using simsons method
      !! @param data_ar 2d input array
      !! @return retrn_ar integration of data_ar
      subroutine array_integral2d(data_ar,retrn_ar,incriment,n,m) 
      double precision :: data_ar(:,:)
      double precision :: retrn_ar(:,:)
      double precision :: trap(n,m)
      double precision :: midpnt(n,m) 
      double precision, dimension(n,m) :: fowrd_ar
      double precision, dimension(n,m) :: back_ar
      double precision :: incriment
      double precision :: sumvl(m)
      integer :: i,j,n,m      

      fowrd_ar(n,:) = data_ar(n,:)
      fowrd_ar(:n,:) = data_ar(2:,:)
  
      back_ar(1,:) = data_ar(1,:)
      back_ar(2:,:) = data_ar(:n,:)

      midpnt = data_ar*incriment
      trap = incriment*(fowrd_ar+back_ar)/2

      sumvl = 0d0*midpnt(1,:)
      
      do i=1,n
        sumvl = sumvl + (2d0/3d0)*midpnt(i,:) + (1d0/3d0)*trap(i,:)
        retrn_ar(i,:) = sumvl
      end do
      
      end subroutine

      !> integrates data_ar along dimension 2
      !! using simsons method
      !! @param data_ar 2d input array
      !! @return retrn_ar integration of data_ar
      subroutine array_integral2dydim(data_ar,retrn_ar,incriment,n,m)
      double precision :: data_ar(:,:)
      double precision :: retrn_ar(:,:)
      double precision :: trap(n,m)
      double precision :: midpnt(n,m)
      double precision, dimension(n,m) :: fowrd_ar
      double precision, dimension(n,m) :: back_ar
      double precision :: incriment
      double precision :: sumvl(n)
      integer :: i,j,n,m

      !using simsons!
      fowrd_ar(:,m) = data_ar(:,m)
      fowrd_ar(:,:m) = data_ar(:,2:)

      back_ar(:,1) = data_ar(:,1)
      back_ar(:,2:) = data_ar(:,:m)

      midpnt = data_ar*incriment
      trap = incriment*(fowrd_ar+back_ar)/2

      sumvl = 0d0*midpnt(:,1)
      do i=1,m
        sumvl = sumvl + (2d0/3d0)*midpnt(:,i) + (1d0/3d0)*trap(:,i)
        retrn_ar(:,i) = sumvl
      end do

      end subroutine

      !> differentiates data_ar along dimension 1
      !! use a center difference
      !! @param data_ar 2d input array
      !! @return retrn_ar integration of data_ar
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

      retrn_ar = (fowrd_ar-back_ar)/(2.0*incriment)
      retrn_ar(1,:) = retrn_ar(2,:)
      retrn_ar(n,:) = retrn_ar(n-1,:)

      end subroutine
      
      !> differentiates data_ar along dimension 2
      !! use a center difference
      !! @param data_ar 2d input array
      !! @return retrn_ar integration of data_ar
      subroutine array_diff2dydim(data_ar,retrn_ar,incriment,n,m)
      double precision :: data_ar(:,:)
      double precision :: retrn_ar(:,:)
      double precision :: incriment
      integer :: n,m
      double precision, dimension(n,m) :: fowrd_ar
      double precision, dimension(n,m) :: back_ar

      fowrd_ar(:,m) = data_ar(:,m)
      fowrd_ar(:,:m) = data_ar(:,2:)

      back_ar(:,1) = data_ar(:,1)
      back_ar(:,2:) = data_ar(:,:m)

      retrn_ar = (fowrd_ar-back_ar)/(2.0*incriment)
      retrn_ar(:,1) = retrn_ar(:,2)
      retrn_ar(:,m) = retrn_ar(:,m-1)

      end subroutine

      !> integral of an array wrt a second array
      !! unused for accuracy reasons
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

      !> differtial of an array wrt a second array
      !! unused for accuracy reasons
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

      !> extends 1d array along second dimesion by
      !! copying it m times
      !! @param data_ar 1d input array
      !! @return retrn_ar 2d output array
      subroutine extend_ardimension(data_ar,retrn_ar,m)
      double precision :: data_ar(:)
      double precision :: retrn_ar(:,:)
      integer :: m, i

      do i=1,m
        retrn_ar(:,i) = data_ar
      end do

      end subroutine

      !> unormalised gaussian used for testing
      !! @param t_ar array giving time variable
      !! @param y_ar array giving depth variable
      !! @param sigma constant specifying gaussian
      !! shape
      !! @return gaussian_nnorm gaussian array
      function gaussian_nnorm(t_ar,y_ar,sigma,n,m)
      double precision :: t_ar(n,m)
      double precision :: y_ar(n,m)
      double precision :: gaussian_nnorm(n,m)
      double precision :: sigma
      integer n,m

      gaussian_nnorm = DEXP(-1d0*((t_ar*t_ar+y_ar*y_ar)/(2*sigma)))    

      end function

      !> unormalised gaussian integral used for 
      !! testing
      !! @param integratedvar_ar coordinate integrated
      !! over
      !! @param normalvar_ar second coordinate
      !! @param alpha constant specifying result
      !! shape (note not the same as sigma in above)
      !! @return gaussian_ideal_integral result array
      function gaussian_ideal_integral &
      &(integratedvar_ar,normalvar_ar,alpha,n,m)
      double precision :: t_ar(n,m)
      double precision :: y_ar(n,m)
      double precision :: integratedvar_ar(n,m)
      double precision :: normalvar_ar(n,m)
      double precision :: gaussian_ideal_integral(n,m)
      double precision :: alpha
      double precision :: PI = 4*DATAN(1d0)
      integer n,m

      gaussian_ideal_integral = DEXP(-alpha*normalvar_ar**2d0) * &
      & DSQRT(PI/alpha)*DERF(DSQRT(alpha)*integratedvar_ar)/2d0

      end function

      end module
