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

      type mineralphase
      
        integer n,m
        integer ID
        character(len=20) :: mineralname
        integer :: parent = 0
        double precision, allocatable :: mineralpart(:,:)
        double precision :: free_energy = 0d0
        double precision :: diffusion_coefficient = 1d0      
        double precision, allocatable :: density(:,:)
        double precision, allocatable :: heatproduction(:,:) 
        double precision, allocatable :: heatcapcity(:,:)
        double precision, allocatable :: thermalconductivity(:,:)
        double precision, allocatable :: bulkmodulus(:,:)
        double precision, allocatable :: grainsize(:,:)
        double precision, allocatable :: phaseline(:,:)
        double precision, allocatable :: velocity(:,:)
        double precision :: lowtemperaturephase(2) = (/0,0/)
        double precision :: lowpressurephase(2) = (/0,0/)

      end type

      contains
      
      !> allocates memory for a mineral
      !! phase
      !! @param this mineral phase instance
      subroutine NEW(this,n,m)
      integer n,m
      type (mineralphase) this
      
      this%n = n
      this%m = m
      allocate( this%mineralpart(n,m) )
      allocate( this%density(n,m) )      
      allocate( this%heatproduction(n,m) )
      allocate( this%heatcapcity(n,m) )
      allocate( this%thermalconductivity(n,m) )
      allocate( this%bulkmodulus(n,m) )
      allocate( this%grainsize(n,m) )
      allocate( this%velocity(n,m) )
      
      this%mineralpart = 0D0
      end subroutine
  
      !> deallocates memory for a mineral
      !! phase
      !! @param this mineral phase instance
      subroutine DELETE(this)
      type (mineralphase) this
      
      deallocate( this%mineralpart )
      deallocate( this%density )  
      deallocate( this%heatproduction )
      deallocate( this%heatcapcity )
      deallocate( this%thermalconductivity )
      deallocate( this%bulkmodulus )
      deallocate( this%grainsize )
      deallocate( this%velocity )      

      end subroutine
  
      !> Returns Temperature of eclogite Clapeyron
      !! slope for a given pressure
      !! @param pressure input pressure array
      !! @return eclogitephase temperature of Slope
      function eclogitephase(pressure,n,m)
      double precision :: pressure(n,m)
      double precision :: eclogitephase(n,m)
      integer n,m

      eclogitephase = 0.13d0*pressure - 97.5d0

      end function 

      !> Returns low temperature condition
      !! for stability given a pressure
      !! @param pressure input pressure array
      !! @return ltphaseline temperature of Slope
      function ltphaseline(this,pressure)
      type (mineralphase) this
      double precision pressure(this%n,this%m)
      double precision :: ltphaseline(this%n,this%m)

      ltphaseline = this%lowtemperaturephase(1)*pressure + &
      &this%lowtemperaturephase(2)

      end function
     
      !> Returns low pressure condition
      !! for stability given a temperature
      !! @param pressure input temperature array
      !! @return lpphaseline pressure of Slope
      function lpphaseline(this,temperature)
      type (mineralphase) this
      double precision temperature(this%n,this%m)
      double precision :: lpphaseline(this%n,this%m)
           
      lpphaseline = this%lowpressurephase(1)*temperature + &
      &this%lowpressurephase(2)

      end function

      !> Gives fraction of mineral phase produced
      !! @param temperature computed temperature
      !! array
      !! @param this mineral phase instance
      !! @param parent parent mineral phase instance
      !! of this
      !! @param t_ar time variable
      subroutine compute_reactionprogress &
      &(this,parent,t_ar,temperature)
      type (mineralphase) this, parent
      double precision :: temperature(:,:)
      double precision :: t_ar(this%n,this%m)
      double precision :: caracteristic_time_ar(this%n,this%m)
      double precision, parameter :: gas_const = 8.3144621 
      
      caracteristic_time_ar = (parent%grainsize**2) * &
      &DEXP(this%free_energy/(gas_const*temperature)) / &
      &parent%diffusion_coefficient

      this%mineralpart = this%mineralpart + &
      &(parent%mineralpart-parent%mineralpart* &
      &DEXP(-t_ar/caracteristic_time_ar))

      parent%mineralpart = parent%mineralpart* &
      &DEXP(-t_ar/caracteristic_time_ar)

      end subroutine

      !> combines reaction progress with
      !! clapeyron slope
      !! @param temperature computed temperature
      !! array
      !! @param pressure input pressure array
      !! @param this mineral phase instance
      !! @param parent parent mineral phase instance
      !! of this
      subroutine compute_part &
      &(this,parent,t_ar,temperature,pressure)
      type (mineralphase) this, parent
      double precision :: temperature(:,:)
      double precision :: pressure(:,:)

      double precision :: t_ar(this%n,this%m)
      double precision, parameter :: gas_const = 8.3144621 !tempory
      
      if (this%parent.NE.0) then
        call compute_reactionprogress(this,parent,t_ar,temperature)
      end if
      
      where (temperature .GE. &
      &ltphaseline(this,pressure))
        where (pressure .GE. &
        & lpphaseline(this,temperature))
          this%mineralpart = this%mineralpart
        elsewhere
          this%mineralpart = 0D0
        end where
      elsewhere
        this%mineralpart = 0D0
      end where

      end subroutine 

      end module
