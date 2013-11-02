      program generatedomain
!234567
      implicit none
      integer :: n,m
      double precision :: incriment(2) 
      double precision :: subduction_angle
      double precision :: PI = 4*DATAN(1d0)
      double precision :: velocity      
      character (len = 256) :: filename

      filename = "output.nc"

      n = 2000
      m = 2000
      incriment(2) = 100d0
      subduction_angle = PI/6d0  
      velocity = 1d58*10**(-9)    

      incriment(1)=(n*incriment(2))/(m*velocity*DSIN(subduction_angle))

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

      integer :: ncid, varid, dimids(2)
      integer :: x_dimid, y_dimid
      integer :: n,m,retval
      integer :: i,j

      do i=1,n
        do j=1,m
          t_ar(i,j) = i/incriment(1)
          y_ar(i,j) = j/incriment(2)
        end do
      end do

      subductLine = -DTAN(subduction_angle)*(t_ar-n*incriment(1))

      print max(y_ar),min(subductLine)

      !should really be ints for ids but concern as to
      !how gmt will process this
      where (y_ar.LE.subductLine)
        where(y_ar.LE.landlith_thickness)
          region_ar = 1000d0
        elsewhere
          region_ar = 3000d0
        endwhere
      elsewhere (y_ar.LE. &
      &(subductLine - (oceanlith_thickness/DCOS(subduction_angle))))
        region_ar = 2000d0
      elsewhere
        region_ar = 3000d0
      end where 

      retval = nf90_create(filename, NF90_CLOBBER, ncid) 

      retval =  nf90_def_dim(ncid, "x", n, x_dimid)
      retval =  nf90_def_dim(ncid, "y", m, y_dimid)
      
      dimids =  (/ x_dimid, y_dimid /)

      retval = nf90_def_var &
      &(ncid, "z", NF90_DOUBLE, dimids, varid)

      retval =  nf90_enddef(ncid)

      retval =  nf90_put_var(ncid, varid, region_ar)

      end subroutine
