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

      module module_MODELFILE
      !use netcdf
      implicit none

      type modelfile
        character (len = 256) :: velocitynetcdf
        character (len = 256) :: densitynetcdf
        character (len = 256) :: heatproductnetcdf
        character (len = 256) :: heatcapcnetcdf
        character (len = 256) :: gtempfile
        character (len = 256) :: gqfluxfile
        character (len = 256) :: thermlconductnetcdf
        character (len = 256) :: outfile
        character (len = 256) :: incompresibilitynetcdf
        integer :: negativedown
        integer :: ydim, tdim 
      end type

      contains 

      !> subroutine to read file specifying
      !! model parameters
      !! @param this derived type storing
      !! file output
      !! @param filename name of file 
      !! specifying model parameters
      subroutine READMDLF(this,filename)
      character (len = *) ::filename
      type (modelfile) this
      integer :: bool = 1      
      character (len = 1) :: typinput
      character (len = 256) :: modelfinput

      OPEN(1, file = filename)
      READ(1,*) this%ydim, this%tdim
      this%negativedown = 0

      do while (bool .EQ. 1)
        READ(1,*) typinput, modelfinput
        if (typinput.EQ."E") then
          bool = 0
        else if (typinput.EQ."V") then
          this%velocitynetcdf = modelfinput
        else if (typinput.EQ."D") then
          this%densitynetcdf = modelfinput
        else if (typinput.EQ."K") then
          this%thermlconductnetcdf = modelfinput
        else if (typinput.EQ."H") then
          this%heatproductnetcdf = modelfinput
        else if (typinput.EQ."C") then
          this%heatcapcnetcdf = modelfinput
        else if (typinput.EQ."T") then
          this%gtempfile = modelfinput
        else if (typinput.EQ."Q") then
          this%gqfluxfile = modelfinput
        else if (typinput.EQ."O") then
          this%outfile = modelfinput
        else if (typinput.EQ."B") then
          this%incompresibilitynetcdf = modelfinput
        else if (typinput.EQ."R") then
          this%negativedown = 1
        end if
      enddo

      end subroutine

      end module
