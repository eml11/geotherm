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
      double precision, allocatable :: mineralpart(:,:)
      double precision :: free_energy
      double precision :: diffusion_coefficient      

      end type

      contains
      
      subroutine NEW(this,n,m)
      integer n,m
      type (mineralphase) this
      
      this%n = n
      this%m = m
      allocate( this%mineralpart(n,m) )
      
      end subroutine
  
      subroutine DELETE(this)
      type (mineralphase) this
      
      deallocate( this%mineralpart )
      
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

      !> Gives fraction of eclogite produced
      !! @param temperature computed temperature
      !! array
      !! @param diffusion_coeficient diffusion
      !! coefficient of limiting chemical species
      !! @param grain_size grain size of parent rock
      !! @param free_energy Gibbs free energy of
      !! reaction
      !! @param t_ar time variable
      !! @param caracteristic_time_ar caracteristic
      !! time for reaction to complete
      !! @return retrn_ar fraction of eclogite against
      !! total rock
      subroutine compute_reactionprogress &
      &(this,t_ar,tfield,model)
      use equationpartsmodule
      use module_modelfile

      type (mineralphase) this
      type (modelfile) model
      type (temperaturefield) tfield

      double precision :: t_ar(this%n,this%m)
      double precision :: caracteristic_time_ar(this%n,this%m)
      double precision, parameter :: gas_const = 8.3144621 

      caracteristic_time_ar = (model%grainsize**2) * &
      &DEXP(this%free_energy/(gas_const*tfield%temperature)) / &
      &this%diffusion_coefficient

      this%mineralpart = 1d0 - DEXP(-t_ar/caracteristic_time_ar)

      end subroutine

      !> combines reaction progress with
      !! clapeyron slope
      !! @param temperature computed temperature
      !! array
      !! @param pressure input pressure array
      !! @param grain_size grain size of parent rock
      !! @param t_ar time variable
      !! @return retrn_ar fraction of eclogite against
      !! total rock
      subroutine compute_eclogite_content &
      &(this,t_ar,tfield,pfield,model)
      use equationpartsmodule
      use module_modelfile
      use pressuresolver
      
      type (mineralphase) this
      type (modelfile) model
      type (pressurefield) pfield
      type (temperaturefield) tfield

      double precision :: t_ar(this%n,this%m)
      double precision, parameter :: gas_const = 8.3144621 !tempory

      call compute_reactionprogress(this,t_ar,tfield,model)

      where (tfield%temperature .LE. &
      &eclogitephase(pfield%pressure,this%n,this%m))
        this%mineralpart = this%mineralpart
      elsewhere
        this%mineralpart = 0d0*this%mineralpart
      end where

      end subroutine 

      end module
