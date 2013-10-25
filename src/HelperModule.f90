      module HelperModule
      use netcdf
      implicit none
!234567
      contains

      subroutine get_netcdf(filename,data_ar)
      character (len = *) :: filename
      double precision :: data_ar(:, :)
      
      integer :: ncid, varid

      integer :: x, y, retval

      !retval is an error checking variable, should be nf_noerr (presumably 0)
      retval = nf90_open(filename, NF90_NOWRITE, ncid)
      retval = nf90_inq_varid(ncid, "z", varid) !change data, use z or something standard
      retval = nf90_get_var(ncid, varid, data_ar)

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
