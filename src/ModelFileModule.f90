      module module_MODELFILE
      !use netcdf
      implicit none
!234567
      type modelfile
        character (len = 256) :: velocitynetcdf
        character (len = 256) :: kappanetcdf
        character (len = 256) :: heatproductnetcdf
        character (len = 256) :: heatcapcnetcdf
        character (len = 256) :: gtempfile
        character (len = 256) :: gqfluxfile
        integer :: ydim, tdim 
      end type

      contains

      subroutine READMDLF(this,filename)
      character (len = *) ::filename
      type (modelfile) this
      integer :: bool = 1      
      character (len = 1) :: typinput
      character (len = 256) :: modelfinput

      OPEN(1, file = filename)
      READ(1,*) this%ydim, this%tdim

      do while (bool .EQ. 1)
        READ(1,*) typinput, modelfinput
        if (typinput.EQ."E") then
          bool = 0
        else if (typinput.EQ."V") then
          this%velocitynetcdf = modelfinput
        else if (typinput.EQ."K") then
          this%kappanetcdf = modelfinput
        else if (typinput.EQ."H") then
          this%heatproductnetcdf = modelfinput
        else if (typinput.EQ."C") then
          this%heatcapcnetcdf = modelfinput
        else if (typinput.EQ."T") then
          this%gtempfile = modelfinput
        else if (typinput.EQ."Q") then
          this%gqfluxfile = modelfinput
        end if
      enddo

      end subroutine

      subroutine set_model_input(this)
      type (modelfile) this

      end subroutine      


      end module
