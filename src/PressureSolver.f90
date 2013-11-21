         
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
      use modeldomainmodule, NEWDOMAIN => NEW, &
      & DELETEDOMAIN => DELETE
      implicit none

      type pressurefield
        integer n,m
        double precision, allocatable :: pressure(:,:)
        double precision, allocatable :: density(:,:)
      endtype


      contains

      subroutine new(this,n,m)
      type (pressurefield) this
      integer n,m

      this%n = n
      this%m = m
      allocate ( this%pressure(n,m) )
      allocate ( this%density(n,m) )

      end subroutine

      !> computes depth dependent Pressure 
      !! @param density_ar Pressure indipendent
      !! density from input
      !! @param incompresibility_ar bulk modulus
      !! of substance
      !! @return retrn_ar depth dependent Pressure
      subroutine compute_pressure(this,domain)
      type (pressurefield) this
      type (modeldomain) domain
      double precision, parameter :: gravity_const = 9.81
      double precision :: y_integral(this%n,this%m)

      call array_integral2dydim(domain%density*gravity_const, & 
      &y_integral,domain%incriment(2),this%n,this%m)

      this%pressure = -domain%bulkmodulus * &
      &DLOG(1+y_integral/domain%bulkmodulus)

      end subroutine

      !> computes pressure dependent density 
      !! @param density_ar Pressure indipendent
      !! density from input
      !! @param incompresibility_ar bulk modulus
      !! of substance
      !! @return retrn_ar pressure dependent density
      subroutine compute_pddensity(this,domain)
      type (pressurefield) this
      type (modeldomain) domain
      double precision, parameter :: gravity_const = 9.81
      double precision :: y_integral(this%n,this%m)

      call array_integral2dydim(domain%density*gravity_const, & 
      &y_integral,domain%incriment(2),this%n,this%m)

      this%density = (domain%density)/(1+y_integral/domain%bulkmodulus)

      end subroutine
      
      subroutine delete(this)
      type (pressurefield) this
      
      deallocate ( this%pressure )
      deallocate ( this%density )

      end subroutine

      end module
