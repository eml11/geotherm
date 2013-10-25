      module HelperModule
      use netcdf
      implicit none
!234567
      contains

      subroutine get_netcdf(filename,data_ar)
      character (len = *) :: filename
      integer :: data_ar(:, :)
      
      integer :: ncid, varid

      integer :: x, y, retval

      !retval is an error checking variable, should be nf_noerr (presumably 0)
      retval = nf90_open(filename, NF90_NOWRITE, ncid)
      retval = nf90_inq_varid(ncid, "data", varid) !change data, use z or something standard
      retval = nf90_get_var(ncid, varid, data_ar)

      retval = nf90_close(ncid)

      end subroutine

      subroutine get_file(filename,data_ar)
      character (len = *) :: filename
      integer :: data_ar(:, :)

      OPEN(1, file = filename)
      READ(1) data_ar

      end subroutine     

      end module
