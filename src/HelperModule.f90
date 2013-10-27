
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
      implicit none

      contains

      subroutine get_netcdf(filename,data_ar)
      character (len = *) :: filename
      double precision :: data_ar(:, :)
      
      integer :: ncid, varid

      integer :: retval
      !retval is an error checking variable, should be nf_noerr (presumably 0)
      retval = nf90_open(filename, NF90_NOWRITE, ncid)
      retval = nf90_inq_varid(ncid, "z", varid) !change data, use z or something standard
      retval = nf90_get_var(ncid, varid, data_ar)
      
      retval = nf90_close(ncid) 

      end subroutine

      subroutine get_netcdf1d(filename,data_ar)
      character (len = *) :: filename
      double precision :: data_ar(:)

      integer :: ncid, varid

      integer :: retval

      !retval is an error checking variable, should be nf_noerr (presumably 0)
      retval = nf90_open(filename, NF90_NOWRITE, ncid)
      retval = nf90_inq_varid(ncid, "z", varid) !change data, use z or something standard
      retval = nf90_get_var(ncid, varid, data_ar)

      retval = nf90_close(ncid)

      end subroutine

      subroutine get_netcdf_wincrmnt(filename,data_ar,incriment,n,m)
  
      character (len = *) :: filename
      double precision :: data_ar(:, :)
      double precision :: incriment(2)
      double precision, dimension(n) :: ydata
      double precision, dimension(m) :: tdata

      integer :: n,m
      integer :: ncid, zvarid, xvarid, yvarid

      integer :: retval

      retval = nf90_open(filename, NF90_NOWRITE, ncid)
      retval = nf90_inq_varid(ncid, "z", zvarid)
      retval = nf90_inq_varid(ncid, "x", xvarid)
      retval = nf90_inq_varid(ncid, "y", yvarid)
      retval = nf90_get_var(ncid, zvarid, data_ar)

      retval = nf90_get_var(ncid, yvarid, ydata)
      retval = nf90_get_var(ncid, xvarid, tdata)

      incriment(1) = tdata(2) - tdata(1)
      incriment(2) = ydata(2) - ydata(1)

      retval = nf90_close(ncid)

      end subroutine

      subroutine get_file(filename,data_ar)
      character (len = *) :: filename
      double precision :: data_ar(:)

      OPEN(1, file = filename)
      READ(1) data_ar

      end subroutine     

      subroutine write_netcdf(filename,data_ar,n,m)
      character (len = *) :: filename
      double precision :: data_ar(n,m)

      integer :: ncid, varid, dimids(2)
      integer :: x_dimid, y_dimid
      integer :: n,m,retval
      
      retval = nf90_create(filename, NF90_CLOBBER, ncid)

      retval = nf90_def_dim(ncid, "x", m, x_dimid)
      retval = nf90_def_dim(ncid, "y", n, y_dimid)
      
      dimids =  (/ y_dimid, x_dimid /)

      !this wont be nf90_int
      retval = nf90_def_var(ncid, "z", NF90_DOUBLE, dimids, varid)
      retval = nf90_enddef(ncid)

      retval = nf90_put_var(ncid, varid, data_ar)
      retval = nf90_close(ncid)

      end subroutine

      end module
