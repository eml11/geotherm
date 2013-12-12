      
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

      module equationpartsmodule
      use modeldomainmodule, DOMAINNEW => NEW, &
      & DOMAINDELETE => DELETE
      use mathmodule
      use pressuresolver, PRESSURENEW => NEW, &
      PRESSUREDELETE => DELETE

      implicit none

      type temperaturefield

      integer n,m
      double precision, allocatable :: temperature(:,:)

      double precision, private, allocatable :: bderivative(:,:)
      double precision, private, allocatable :: expterm(:,:)
      double precision, private, allocatable :: innerintegral(:,:)
      double precision, private, allocatable :: outerintegral(:,:)
      double precision, private, allocatable :: innerconstant(:,:)
      double precision, private, allocatable :: outerconstant(:,:)

      end type

      contains


      !> allocates memory for a new
      !! temperature field
      !! @param this temperature field instance
      subroutine new(this,n,m)
      type (temperaturefield) this
      integer n,m
      this%n = n
      this%m = m
      allocate ( this%expterm(n,m) )
      allocate ( this%innerintegral(n,m) )
      allocate ( this%outerintegral(n,m) )
      allocate ( this%innerconstant(n,m) )
      allocate ( this%outerconstant(n,m) )
      allocate ( this%bderivative(n,m) )

      end subroutine

      !> sets retrn_ar to an array corrisponding to the
      !! time differential of b(t)
      !! @param tdata_ar ground temperature array
      !! @param qdata_ar ground heat flux array
      !! @return retrn_ar time derivative of b(t)
      subroutine compute_bdashval(this,domain)
      type (modeldomain) domain
      type (temperaturefield) this
      double precision, dimension(this%n,this%m) :: difftdata_ar

      call array_diff2d(domain%gtemp, &
      &difftdata_ar,domain%incriment(1),this%n,this%m)
      this%bderivative = domain%thermalconductivity * &
      &(difftdata_ar/domain%gqflux)

      end subroutine

      !> unused routine
      subroutine compute_bdashval_fromqx &
      &(tdata_ar,qdata_ar,qxdata_ar,kconstant,retrn_ar,incriment,n,m)
      double precision :: tdata_ar(n,m)
      double precision :: tfromqx_ar(n,m)
      double precision :: qdata_ar(n,m)
      double precision :: qxdata_ar(n,m)
      double precision :: kconstant(n,m)
      double precision :: retrn_ar(n,m)
      double precision :: incriment(2)
      integer :: n, m
      double precision, dimension(n,m) :: difftdata_ar

      call array_integral2d(qxdata_ar/kconstant,tfromqx_ar, &
      &incriment(2),n,m)
      call array_diff2d(-tfromqx_ar,difftdata_ar, &
      &incriment(1),n,m) 
      retrn_ar = kconstant * &
      &(difftdata_ar/qdata_ar)

      end subroutine

      !> sets retrn_ar to the exponented function
      !! in the solution to the homogeneous form
      !! @param bdash_ar time derivative of b(t)
      !! @param velocity_ar velocity of motion
      !! perpendicular to y direction
      !! @param kappa_ar thermal diffusivity
      !! @return retrn_ar function in exponent
      subroutine compute_exponentintegral(this,domain,pfield)
      type (modeldomain) domain
      type (temperaturefield) this
      type (pressurefield) pfield
      double precision :: kappa_ar(this%n,this%m)
      double precision, dimension(this%n,this%m) :: v_tintegral
      double precision, dimension(this%n,this%m) :: b_tintegral
      double precision, dimension(this%n,this%m) :: v_yintegral
      double precision, dimension(this%n,this%m) :: b_yintegral     
      double precision, dimension(this%n,this%m) :: argument

      kappa_ar = &
      &domain%thermalconductivity/ &
      &(domain%heatcapcity*pfield%density)

      argument = (this%bderivative*domain%velocity)/kappa_ar
      call array_integral2d(argument, &
      &v_tintegral,domain%incriment(1),this%n,this%m)
      
      argument = (this%bderivative*this%bderivative)/kappa_ar
      
      call array_integral2d(argument, &
      &b_tintegral,domain%incriment(1),this%n,this%m)
      
      argument = domain%velocity/kappa_ar
      
      call array_integral2dydim(argument, &
      &v_yintegral,domain%incriment(2),this%n,this%m)
      
      argument = this%bderivative/kappa_ar
      call array_integral2dydim(argument, &
      &b_yintegral,domain%incriment(2),this%n,this%m) 

      this%expterm = -v_tintegral + b_tintegral + &
      &v_yintegral - b_yintegral
      end subroutine

      !> sets retrn_ar to the first integral wrt
      !! eta (without specific constant)
      !! @param exintegral_ar function in exponent
      !! @param bdash_ar time derivative of b(t)
      !! @param A_ar 
      !! @param thermal_ar thermal conductivity
      !! @return retrn_ar first integration result
      subroutine compute_init_inerintegral(this,domain)
      type (modeldomain) domain
      type (temperaturefield) this
      double precision, dimension(this%n,this%m) :: integral_term
      double precision, dimension(this%n,this%m) :: t_integral
      double precision, dimension(this%n,this%m) :: y_integral
      
      integral_term = -domain%heatproduction*domain%density * &
      &DEXP(-1*this%expterm)/domain%thermalconductivity
      call array_integral2dydim &
      &(integral_term,y_integral,domain%incriment(2),this%n,this%m)
      
      integral_term = -domain%heatproduction*domain%density * &
      &this%bderivative * &
      &DEXP(-1*this%expterm)/domain%thermalconductivity
      call array_integral2d &
      &(integral_term, &
      &t_integral,domain%incriment(1),this%n,this%m)

      this%innerintegral = y_integral - t_integral
      end subroutine

      !> computes constant of integration for integral
      !! returned by subroutine compute_init_inerintegral
      !! sets ground flux to given boundry condition
      !! @param inerintegral_ar first integration result
      !! @param exintegral_ar function in exponent
      !! @param tdata_ar ground temperature array
      !! @param qdata_ar ground heat flux array
      !! @param bdash_ar time derivative of b(t)
      !! @param thermal_ar thermal conductivity
      !! @return retrn_dbl first constant of integration
      subroutine compute_inerintegralconstant(this,domain)
      type (modeldomain) domain
      type (temperaturefield) this
      double precision :: tdiff_ar(this%n,this%m)
      double precision :: tdata_diff_ar(this%n,this%m)
      double precision :: innerconstant(this%n)

      call array_diff2d(domain%gtemp,tdata_diff_ar, &
      &domain%incriment(1),this%n,this%m)

      where (this%bderivative.EQ.0d0)
        tdiff_ar = 0d0
      elsewhere
        tdiff_ar = tdata_diff_ar/this%bderivative
      end where

      innerconstant = &
      &(-(domain%gqflux(:,3)/domain%thermalconductivity(:,3)) + &
      &tdiff_ar(:,3)) * &
      &EXP(-this%expterm(:,3)) - this%innerintegral(:,3)

      call extend_ardimension(innerconstant,this%innerconstant,this%m)

      end subroutine

      !> sets retrn_ar to the second integral wrt
      !! eta (without specific constant)
      !! @param first integration result
      !! @param exintegral_ar function in exponent
      !! @param iner_dbl first constant of integration
      !! @return second integration result
      subroutine compute_init_outerintegral(this,domain)
      type (modeldomain) domain
      type (temperaturefield) this
      double precision, dimension(this%n,this%m) :: integral_term
      double precision, dimension(this%n,this%m) :: t_integral
      double precision, dimension(this%n,this%m) :: y_integral

      integral_term = &
      &EXP(this%expterm)*(this%innerintegral + this%innerconstant)

      call array_integral2dydim &
      &(integral_term,y_integral,domain%incriment(2),this%n,this%m)

      integral_term = this%bderivative * &
      &EXP(this%expterm)*(this%innerintegral + this%innerconstant)
      call array_integral2d(integral_term, &
      &t_integral,domain%incriment(1),this%n,this%m)

      this%outerintegral = y_integral - t_integral

      end subroutine

      !> computes constant of integration for integral
      !! returned by subroutine compute_init_outerintegral
      !! sets ground temperature to given boundry
      !! condition
      !! @param outerintegral_ar second integration result
      !! @param tdata_ar ground temperature array
      !! @return retrn_dbl second constant of integration
      subroutine compute_outerintegralconstant(this,domain)
      type (modeldomain) domain
      type (temperaturefield) this
      double precision :: outerintegral_ar(this%n,this%m)
      double precision :: tdata_ar(this%n,this%m)
      double precision outerconstant(this%n)

      outerconstant = domain%gtemp(:,1) - this%outerintegral(:,1)

      call extend_ardimension(outerconstant,this%outerconstant,this%m)

      end subroutine

      !> combines above routines to compute
      !! temperature
      !! @param this temperature field instance
      !! @param domain domain instance of the model
      !! @param pfield pressure field instance
      subroutine compute_temp(this,domain,pfield)
      type (modeldomain) domain
      type (temperaturefield) this
      type (pressurefield) pfield
   
      write(2,*) "subroutine compute_bdashval"
      call compute_bdashval(this,domain)
      write(2,*) "subroutine compute_exponentintegral" 
      call compute_exponentintegral(this,domain,pfield)
      write(2,*) "subroutine compute_init_inerintegral"
      call compute_init_inerintegral(this,domain)
      write(2,*) "subroutine compute_inerintegralconstant"
      call compute_inerintegralconstant(this,domain)
      write(2,*) "subroutine compute_init_outerintegral"
      call compute_init_outerintegral(this,domain)
      write(2,*) "subroutine compute_outerintegralconstant"
      call compute_outerintegralconstant(this,domain)

      this%temperature = this%outerintegral + this%outerconstant

      write(2,*) "finished computing temperature field"

      end subroutine
      
      !> deallocates memory for a new
      !! temperature field
      !! @param this temperature field instance
      subroutine delete(this)
      type (temperaturefield) this
      
      deallocate ( this%expterm )
      deallocate ( this%innerintegral )
      deallocate ( this%outerintegral )
      deallocate ( this%innerconstant )
      deallocate ( this%outerconstant )
      deallocate ( this%bderivative )
      
      end subroutine
      
      end module
