      

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


      !call array_diff2d(RESHAPE(tdata_ar,(/m,n/)), &
      !&difftdata_ar,incriment(1),m,n)
      !retrn_ar = kconstant * &
      !&(RESHAPE(difftdata_ar,(/n,m/))/qdata_ar)

      call array_diff2d(tdata_ar, &
      &difftdata_ar,incriment(1),n,m)
      retrn_ar = kconstant * &
      &(difftdata_ar/qdata_ar)

      end subroutine

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
      double precision, dimension(m,n) :: v_yintegral, b_yintegral

      !this is correct reshape for y

      call array_integral2d &
      &((bdash_ar*velocity_ar)/kappa_ar, &
      &v_tintegral,incriment(1),n,m)
      call array_integral2d &
      ((bdash_ar*bdash_ar)/kappa_ar, &
      &b_tintegral,incriment(1),n,m)
      call array_integral2d(RESHAPE(velocity_ar/kappa_ar,(/m,n/)), &
      &v_yintegral,incriment(2),m,n)
      call array_integral2d( &
      &RESHAPE(1/kappa_ar,(/m,n/)),b_yintegral,incriment(2),m,n)


      retrn_ar = -v_tintegral+ &
      &b_tintegral + &
      &RESHAPE(v_yintegral,(/n,m/)) -  &
      &RESHAPE(bdash_ar*b_yintegral,(/n,m/))

      end subroutine

      subroutine compute_init_inerintegral &
      &(exintegral_ar,bdash_ar,theta_ar,kappa_ar,retrn_ar,incriment,n,m)
      use mathmodule
      double precision :: exintegral_ar(n,m)
      double precision :: bdash_ar(n,m)
      double precision :: retrn_ar(n,m)
      double precision :: theta_ar(n,m)
      double precision :: kappa_ar(n,m)
      double precision :: incriment(2)
      double precision, dimension(n,m) :: integral_term, y_integral
      double precision, dimension(m,n) :: t_integral

      integral_term = theta_ar*kappa_ar*DEXP(-1*exintegral_ar)
      call array_integral2d(integral_term,y_integral,incriment(2),n,m)
      call array_integral2d &
      &(RESHAPE(integral_term*bdash_ar,(/m,n/)), &
      &t_integral,incriment(1),m,n)

      retrn_ar = y_integral + RESHAPE(t_integral,(/n,m/))

      end subroutine

      subroutine compute_inerintegralconstant &
      &(inerintegral_ar,exintegral_ar,qdata_ar,kconstant,retrn_dbl,n,m)
      double precision :: inerintegral_ar(n,m)
      double precision :: exintegral_ar(n,m)
      double precision :: qdata_ar(n,m)
      double precision :: kconstant(n,m)
      double precision :: retrn_dbl

      retrn_dbl = &
      &((-1*qdata_ar(1,1)*qdata_ar(1,1))/kconstant(1,1)) * &
      &EXP(exintegral_ar(1,1)) - inerintegral_ar(1,1)

      end subroutine

      subroutine compute_init_outerintegral &
      &(inerintegral_ar,exintegral_ar,iner_dbl,retrn_ar,incriment,n,m)
      use mathmodule
      double precision :: inerintegral_ar(n,m)
      double precision :: exintegral_ar(n,m)
      double precision :: retrn_ar(n,m)
      double precision :: iner_dbl
      double precision :: incriment(2)
      double precision, dimension(n,m) :: integral_term, y_integral
      double precision, dimension(m,n) :: t_integral

      integral_term = EXP(exintegral_ar)*(inerintegral_ar + iner_dbl)
      call array_integral2d(integral_term,y_integral,incriment(2),n,m)
      call array_integral2d &
      &(RESHAPE(integral_term*bdash_ar,(/m,n/)), &
      &t_integral,incriment(1),m,n)

      retrn_ar = y_integral + RESHAPE(t_integral,(/n,m/))

      end subroutine

      subroutine compute_outerintegralconstant &
      &(outerintegral_ar,tdata_ar,retrn_dbl,n,m)
      double precision :: outerintegral_ar(n,m)
      double precision :: tdata_ar(n,m)
      double precision retrn_dbl

      retrn_dbl = tdata_ar(1,1) - outerintegral_ar(1,1)

      end subroutine

      end module

