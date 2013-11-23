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

      module HelperModule
      use netcdf
      use mathmodule
      implicit none

      contains

      !> basic netcdf reading subrutine
      !! @param filename name of netcdf file
      !! @return data_ar z data of netcdf
      subroutine get_netcdf(filename,data_ar)
      character (len = *) :: filename
      double precision :: data_ar(:, :)
            
      integer :: ncid, varid

      integer :: retval
      
      call check( nf90_open(filename, NF90_NOWRITE, ncid) )
      call check( nf90_inq_varid(ncid, "z", varid) )
      call check( nf90_get_var(ncid, varid, data_ar) )
      
      call check( nf90_close(ncid) )

      end subroutine

      !> ID netcdf reading subrutine
      !! @param filename name of netcdf file
      !! @return data_ar z data of netcdf
      subroutine get_idnetcdf(filename,data_ar,incriment,n,m)

      character (len = *) :: filename
      integer :: data_ar(:, :)
      double precision :: incriment(2)
      double precision, dimension(m) :: ydata
      double precision, dimension(n) :: tdata

      integer :: n,m
      integer :: ncid, zvarid, xvarid, yvarid

      integer :: retval

      call check( nf90_open(filename, NF90_NOWRITE, ncid) )
      call check( nf90_inq_varid(ncid, "z", zvarid) )
      call check( nf90_inq_varid(ncid, "x", xvarid) )
      call check( nf90_inq_varid(ncid, "y", yvarid) )
      call check( nf90_get_var(ncid, zvarid, data_ar) )

      call check( nf90_get_var(ncid, yvarid, ydata) )
      call check( nf90_get_var(ncid, xvarid, tdata) )

      incriment(1) = tdata(2) - tdata(1)
      incriment(2) = ydata(2) - ydata(1)

      !hack to fix netcdfs from gmt
      incriment(2) = -incriment(2)

      call check( nf90_close(ncid) )

      end subroutine


      !> 1d netcdf reading subroutine used for
      !! boundry conditions
      !! @param filename name of netcdf file
      !! @return data_ar z data of netcdf      
      subroutine get_netcdf1d(filename,retrn_ar,n,m)
      character (len = *) :: filename
      double precision :: data_ar(n)
      double precision :: retrn_ar(:,:)
      
      integer :: ncid, varid
      integer :: n,m
      integer :: retval

      call check( nf90_open(filename, NF90_NOWRITE, ncid) )
      call check( nf90_inq_varid(ncid, "z", varid) )
      call check( nf90_get_var(ncid, varid, data_ar) )

      call check( nf90_close(ncid) )

      call extend_ardimension(data_ar,retrn_ar,m)

      end subroutine

      !> netcdf reading subroutine called initially
      !! in order to get incriment of y/t
      !! @param filename name of netcdf file
      !! @return data_ar z data of netcdf 
      !! @return incriment spacing of y/t variables
      subroutine get_netcdf_wincrmnt(filename,data_ar,incriment,n,m)
  
      character (len = *) :: filename
      double precision :: data_ar(:, :)
      double precision :: incriment(2)
      double precision, dimension(m) :: ydata
      double precision, dimension(n) :: tdata

      integer :: n,m
      integer :: ncid, zvarid, xvarid, yvarid

      integer :: retval
      
      call check( nf90_open(filename, NF90_NOWRITE, ncid) )
      call check( nf90_inq_varid(ncid, "z", zvarid) )
      call check( nf90_inq_varid(ncid, "x", xvarid) )
      call check( nf90_inq_varid(ncid, "y", yvarid) )
      call check( nf90_get_var(ncid, zvarid, data_ar) )

      call check( nf90_get_var(ncid, yvarid, ydata) )
      call check( nf90_get_var(ncid, xvarid, tdata) )

      incriment(1) = tdata(2) - tdata(1)
      incriment(2) = ydata(2) - ydata(1)
      
      !hack to fix netcdfs from gmt
      incriment(2) = -incriment(2)

      call check( nf90_close(ncid) )

      end subroutine

      !> obsolete boundry condition file
      !! @param filename name of boundry file
      !! @return data_ar data from file    
      subroutine get_file(filename,data_ar)
      character (len = *) :: filename
      double precision :: data_ar(:)

      OPEN(1, file = filename)
      READ(1) data_ar

      end subroutine 
    
      !> netcdf error checking subroutine
      !! @param status return value of function
      !! from netcdf-fortran
      subroutine check(status)
      integer, intent ( in) :: status
    
      write(2,*) trim(nf90_strerror(status))
      if(status.NE.nf90_noerr) then
        stop "Stopped"
      end if
      end subroutine

      end module
