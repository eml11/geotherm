         
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
