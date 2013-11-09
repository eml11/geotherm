      
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

      contains

      !> sets retrn_ar to an array corrisponding to the
      !! time differential of b(t)
      !! @param tdata_ar ground temperature array
      !! @param qdata_ar ground heat flux array
      !! @return retrn_ar time derivative of b(t)
      subroutine compute_bdashval &
      &(tdata_ar,qdata_ar,kconstant,retrn_ar,incriment,n,m)
      use mathmodule
      double precision :: tdata_ar(n,m)
      double precision :: qdata_ar(n,m)
      double precision :: kconstant(n,m)
      double precision :: retrn_ar(n,m)
      double precision :: incriment(2)
      integer :: n, m
      double precision, dimension(n,m) :: difftdata_ar

      call array_diff2d(tdata_ar, &
      &difftdata_ar,incriment(1),n,m)
      retrn_ar = kconstant * &
      &(difftdata_ar/qdata_ar)

      end subroutine

      subroutine compute_bdashval_fromqx &
      &(tdata_ar,qdata_ar,qxdata_ar,kconstant,retrn_ar,incriment,n,m)
      use mathmodule
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
      call array_diff2d(-tfromqx_ar,difftdata_ar,incriment(1),n,m) 
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
      subroutine compute_exponentintegral &
      &(bdash_ar,velocity_ar,kappa_ar,retrn_ar,incriment,n,m)
      use mathmodule
      double precision :: bdash_ar(n,m)
      double precision :: velocity_ar(n,m)
      double precision :: kappa_ar(n,m)
      double precision :: retrn_ar(n,m)
      double precision :: incriment(2)
      integer n,m
      double precision, dimension(n,m) :: v_tintegral, b_tintegral
      double precision, dimension(n,m) :: v_yintegral, b_yintegral
      
      call array_integral2d &
      &((bdash_ar*velocity_ar)/kappa_ar, &
      &v_tintegral,incriment(1),n,m)
      call array_integral2d &
      ((bdash_ar*bdash_ar)/kappa_ar, &
      &b_tintegral,incriment(1),n,m)
      call array_integral2dydim(velocity_ar/kappa_ar, &
      &v_yintegral,incriment(2),n,m)
      call array_integral2dydim(bdash_ar/kappa_ar, &
      &b_yintegral,incriment(2),n,m) 

      retrn_ar = -v_tintegral + b_tintegral + &
      &v_yintegral - b_yintegral

      !print *, retrn_ar(100,:)

      end subroutine

      !> sets retrn_ar to the first integral wrt
      !! eta (without specific constant)
      !! @param exintegral_ar function in exponent
      !! @param bdash_ar time derivative of b(t)
      !! @param A_ar 
      !! @param thermal_ar thermal conductivity
      !! @return retrn_ar first integration result
      subroutine compute_init_inerintegral &
      &(exintegral_ar,bdash_ar,A_ar,thermal_ar,retrn_ar,incriment,n,m)
      use mathmodule
      double precision :: exintegral_ar(n,m)
      double precision :: bdash_ar(n,m)
      double precision :: retrn_ar(n,m)
      double precision :: A_ar(n,m)
      double precision :: thermal_ar(n,m)
      double precision :: incriment(2)
      double precision, dimension(n,m) :: integral_term, t_integral
      double precision, dimension(n,m) :: y_integral
      !A_ar and thermal_ar incorrect, hacked for now
      integral_term = -A_ar*DEXP(-1*exintegral_ar)/thermal_ar
      call array_integral2dydim &
      &(integral_term,y_integral,incriment(2),n,m)
      call array_integral2d &
      &(integral_term*bdash_ar, &
      &t_integral,incriment(1),n,m)

      !print *, y_integral(100,:)

      retrn_ar = y_integral - t_integral

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
      subroutine compute_inerintegralconstant &
      &(inerintegral_ar,exintegral_ar,tdata_ar,qdata_ar,bdash_ar, &
      &thermal_ar,retrn_dbl,incriment,n,m)
      use mathmodule
      double precision :: inerintegral_ar(n,m)
      double precision :: exintegral_ar(n,m)
      double precision :: qdata_ar(n,m)
      double precision :: tdata_ar(n,m)
      double precision :: bdash_ar(n,m)
      double precision :: tdiff_ar(n,m)
      double precision :: thermal_ar(n,m)
      double precision :: retrn_dbl(n)      
      double precision :: tdata_diff_ar(n,m)
      double precision :: incriment(2)
      integer :: n,m

      call array_diff2d(tdata_ar,tdata_diff_ar,incriment(1),n,m)

      where (bdash_ar.EQ.0d0)
        tdiff_ar = 0d0
      elsewhere
        tdiff_ar = tdata_diff_ar/bdash_ar
      end where

      retrn_dbl = &
      &(-(qdata_ar(:,3)/thermal_ar(:,3))+tdiff_ar(:,3)) * &
      &EXP(-exintegral_ar(:,3)) - inerintegral_ar(:,3)

      !print *, thermal_ar(100,:)*EXP(exintegral_ar(100,:)) * &
      !&(inerintegral_ar(100,:) + retrn_dbl(100))

      end subroutine

      !> sets retrn_ar to the second integral wrt
      !! eta (without specific constant)
      !! @param first integration result
      !! @param exintegral_ar function in exponent
      !! @param iner_dbl first constant of integration
      !! @return second integration result
      subroutine compute_init_outerintegral &
      &(inerintegral_ar,exintegral_ar,bdash_ar,iner_dbl, &
      &retrn_ar,incriment,n,m)
      use mathmodule
      double precision :: inerintegral_ar(n,m)
      double precision :: exintegral_ar(n,m)
      double precision :: bdash_ar(n,m)
      double precision :: retrn_ar(n,m)
      double precision :: iner_dbl(n,m)
      double precision :: incriment(2)
      double precision, dimension(n,m) :: integral_term, t_integral
      double precision, dimension(n,m) :: y_integral

      integral_term = EXP(exintegral_ar)*(inerintegral_ar + iner_dbl)
      
      !print *, iner_dbl(100,:)

      call array_integral2dydim &
      &(integral_term,y_integral,incriment(2),n,m)
      call array_integral2d(integral_term*bdash_ar, &
      &t_integral,incriment(1),n,m)

      retrn_ar = y_integral - t_integral

      !print *, y_integral(100,:)

      end subroutine

      !> computes constant of integration for integral
      !! returned by subroutine compute_init_outerintegral
      !! sets ground temperature to given boundry
      !! condition
      !! @param outerintegral_ar second integration result
      !! @param tdata_ar ground temperature array
      !! @return retrn_dbl second constant of integration
      subroutine compute_outerintegralconstant &
      &(outerintegral_ar,tdata_ar,retrn_dbl,n,m)
      double precision :: outerintegral_ar(n,m)
      double precision :: tdata_ar(n,m)
      double precision retrn_dbl(n)

      retrn_dbl = tdata_ar(:,1) - outerintegral_ar(:,1)

      !print *, outerintegral_ar(500,1) + retrn_dbl

      end subroutine

      end module
