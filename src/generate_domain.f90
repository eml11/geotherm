
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

      program generatedomain

      implicit none
      integer :: n,m
      double precision :: incriment(2) 
      double precision :: subduction_angle
      double precision :: PI = 4*DATAN(1d0)
      double precision :: velocity      
      character (len = 256) :: filename

      filename = "output.nc"

      n = 2000
      m = 1200
      incriment(2) = 100d0
      subduction_angle = PI/6d0  
      velocity = 1.58*(10**(-9.0))    

      incriment(1)=(n*incriment(2))/(m*velocity*DSIN(subduction_angle))

      print *, incriment
      print *, incriment(1)*n,incriment(2)*m

      call work_routine(filename,subduction_angle,incriment,n,m)

      end program

      subroutine work_routine(filename,subduction_angle,incriment,n,m)
      use netcdf
      character (len = *) :: filename
      double precision :: subduction_angle
      double precision :: incriment(2)
      double precision :: region_ar(n,m)  
      double precision :: t_ar(n,m),y_ar(n,m)
      double precision :: subductLine(n,m)
      double precision :: oceanlith_thickness = 5000d0
      double precision :: landlith_thickness = 30000d0

      integer :: ncid, xvarid, yvarid, zvarid, dimids(2)
      integer :: x_dimid, y_dimid
      integer :: n,m,retval
      integer :: i,j

      do i=1,n
        do j=1,m
          t_ar(i,j) = i*incriment(1)
          y_ar(i,j) = j*incriment(2)
        end do
      end do

      subductLine = DTAN(subduction_angle)*t_ar * &
      &(incriment(2)/incriment(1))

      !use this to get around problem with gmt at some point
      !should really be ints for ids but concern as to
      !how gmt will process this
      where (y_ar.LE.subductLine)
        where(y_ar.LE.landlith_thickness)
          region_ar = 1000d0
        elsewhere
          region_ar = 3000d0
        endwhere
      elsewhere (y_ar.LE. &
      &(subductLine + (oceanlith_thickness/DCOS(subduction_angle))))
        region_ar = 2000d0
      elsewhere
        region_ar = 3000d0
      end where 

      call check( nf90_create(filename, NF90_CLOBBER, ncid) )
      
      call check(  nf90_def_dim(ncid, "x", n, x_dimid) )
      call check(  nf90_def_dim(ncid, "y", m, y_dimid) )
      
      dimids =  (/ x_dimid, y_dimid /)

      !call check( nf90_def_var &
      !&(ncid, "x", NF90_DOUBLE, x_dimid, xvarid)) 

      !call check( nf90_def_var &
      !&(ncid, "y", NF90_DOUBLE, y_dimid, yvarid))

      call check( nf90_def_var &
      &(ncid, "z", NF90_DOUBLE, dimids, zvarid) )

      call check(  nf90_enddef(ncid) )

      !call check(  nf90_put_var(ncid, xvarid, t_ar(:,1)) )
      !call check(  nf90_put_var(ncid, yvarid, y_ar(1,:)) )
      call check(  nf90_put_var(ncid, zvarid, region_ar) )

      end subroutine

      subroutine check(status)
      use netcdf
      integer, intent ( in) :: status
    
      if(status.NE.nf90_noerr) then 
        print *, trim(nf90_strerror(status))
        stop "Stopped"
      end if
      end subroutine
