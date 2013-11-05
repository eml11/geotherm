
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

      module geochem
      implicit none

      contains
  
      function eclogitephase(pressure,n,m)
      double precision :: pressure(n,m)
      double precision :: eclogitephase(n,m)
      integer n,m

      eclogitephase = 0d13*pressure - 97d5

      end function 

      subroutine compute_reactionprogress &
      &(t_ar,temperature,diffusion_coeficient,grain_size,free_energy, &
      &retrn_ar,n,m)
      double precision :: temperature(n,m)
      double precision :: diffusion_coeficient(n,m)
      double precision :: grain_size(n,m)
      double precision :: free_energy(n,m)
      double precision :: retrn_ar(n,m)
      double precision :: t_ar(n,m)
      double precision :: caracteristic_time_ar(n,m)
      double precision, parameter :: gas_const 
      integer n,m

      caracteristic_time_ar = (grain_size(n,m)**2) * &
      &DEXP(free_energy/(gas_const*temperature))/diffusion_coeficient

      retrn_ar = 1d0 - DEXP(t_ar/caracteristic_time_ar)

      end subroutine

      subroutine compute_eclogite_content &
      &(t_ar,temperature,pressure,diffusion_coeficient, &
      &grain_size,free_energy,retrn_ar,n,m)
      double precision :: temperature(n,m)
      double precision :: pressure(n,m)
      double precision :: diffusion_coeficient(n,m)
      double precision :: grain_size(n,m)
      double precision :: free_energy(n,m)
      double precision :: retrn_ar(n,m)
      double precision :: t_ar(n,m)
      double precision :: caracteristic_time_ar(n,m)
      double precision, parameter :: gas_const
      integer n,m

      call compute_reactionprogress &
      &(t_ar,temperature,diffusion_coeficient,grain_size,free_energy, &
      &retrn_ar,n,m)

      where (temperature.GE.eclogitephase(pressure,n,m))
        retrn_ar = retrn_ar
      elsewhere
        retrn_ar = 0d0
      end where

      end module
