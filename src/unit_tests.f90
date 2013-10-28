
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
      double precision, dimension(m,n) :: reshaped_ar
      double precision :: ingrl_cnst
 
      double precision :: incriment(2) = (/0.001,0.001/)

      do i=1,n
        do j=1,m
          t_ar(i,j) = i/1000.
          y_ar(i,j) = j/1000.
        end do
      end do
    
      wrk_ar1 = gaussian_nnorm(t_ar,y_ar,2d0,n,m)
      wrk_ar2 = -t_ar*wrk_ar1/2

      call array_diff2d(wrk_ar1,wrk_ar3,incriment(1),n,m)

      call test_output(wrk_ar2,wrk_ar3,n,m)

      call array_integral2d(wrk_ar1,wrk_ar2,incriment(1),n,m)
      call array_diff2d(wrk_ar2,wrk_ar3,incriment(1),n,m)

      call test_output(wrk_ar1,wrk_ar3,n,m)

      wrk_ar1 = gaussian_nnorm(t_ar,y_ar*0,2d0,n,m)
      wrk_ar2 = gaussian_nnorm(t_ar,y_ar*0,3d0,n,m)
      wrk_ar3 = gaussian_nnorm(t_ar,y_ar*0,5d0,n,m)
      wrk_ar4 = (-t_ar*wrk_ar1*wrk_ar3)/(wrk_ar2*2)

      call compute_bdashval &
      &(wrk_ar1,wrk_ar2,wrk_ar3,wrk_ar5,incriment,n,m)

      !assume fine for now issue with periodic section of ones
      !which must be explained, possibly down to RESHAPE

      call test_output(wrk_ar4,wrk_ar5,n,m)

      wrk_ar2 = gaussian_nnorm(t_ar,y_ar,3d0,n,m)
      wrk_ar3 = gaussian_nnorm(t_ar,y_ar,5d0,n,m)
      wrk_ar4 = (wrk_ar2-wrk_ar1)/wrk_ar3

      call compute_exponentintegral &
      &(wrk_ar1,wrk_ar2,wrk_ar3,wrk_ar5,incriment,n,m)

      call array_diff2d(wrk_ar5,wrk_ar2,incriment(1),n,m)
      call array_diff2d(RESHAPE(wrk_ar5,(/m,n/)), &
      &reshaped_ar,incriment(1),m,n)
      wrk_ar6 = wrk_ar2-RESHAPE(reshaped_ar,(/n,m/))/wrk_ar1

      call test_output(wrk_ar4,wrk_ar6,n,m)

      !print *, ABS( &
      !&(wrk_ar4(2:n,:)-wrk_ar6(2:n,:)+ingrl_cnst)/wrk_ar4(2:n,:))
 
      wrk_ar4 = gaussian_nnorm(t_ar,y_ar,7d0,n,m)
      wrk_ar5 = (wrk_ar2*wrk_ar3*wrk_ar4)/(y_ar*(-1/3.0-0.2-1/7.0)) + &
      &(wrk_ar1*wrk_ar2*wrk_ar3*wrk_ar4)/(t_ar*(-1/3.0-0.2-1/7.0-0.5))

      call compute_init_inerintegral &
      &(wrk_ar2,wrk_ar1,wrk_ar3,wrk_ar4,wrk_ar6,incriment,n,m)

      call test_output(wrk_ar5,wrk_ar6,n,m)

      end program

      subroutine test_output(ideal_ar,computed_ar,n,m)
      integer :: n,m
      double precision :: ideal_ar(n,m)
      double precision :: computed_ar(n,m)

      print *,
      print *,
      print *, computed_ar(4,4)," ", ideal_ar(4,4)
      print *, MAXVAL( &
      &ABS((ideal_ar(2:n,:)-computed_ar(2:n,:))/ideal_ar(2:n,:)))
      print *, MINVAL( &
      &ABS((ideal_ar(2:n,:)-computed_ar(2:n,:))/ideal_ar(2:n,:)))

      end subroutine
