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
        character (len = 256) :: gqxfluxfile
        character (len = 256) :: thermlconductnetcdf
        character (len = 256) :: outfile
        character (len = 256) :: incompresibilitynetcdf
        character (len = 256) :: grainsizenetcdf
        integer :: negativedown
        integer :: ydim, tdim
        double precision :: incriment(2)
        double precision, allocatable :: velocity(:,:)!Domain
        double precision, allocatable :: density(:,:)!mineral
        double precision, allocatable :: heatproduction(:,:)!mineral
        double precision, allocatable :: heatcapcity(:,:)!mineral
        double precision, allocatable :: gtemp(:,:)!Domain
        double precision, allocatable :: gqflux(:,:)!Domain
        double precision, allocatable :: thermalconductivity(:,:)!mineral
        double precision, allocatable :: bulkmodulus(:,:)!mineral
        double precision, allocatable :: grainsize(:,:) !mineral
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
      integer :: booltwo = 1
      integer :: boolthree = 1      
      character (len = 1) :: typinput
      character (len = 256) :: modelfinput
      integer ID, num, i
      double precision part

      OPEN(1, file = filename)
      READ(1,*) this%ydim, this%tdim
      this%negativedown = 0

      do while (bool.EQ.1)
        booltwo = 1
        READ(1,*) modelfinput
        if (modelfinput.EQ."Domain") then
          do while (booltwo.EQ.1)
            READ(1,*) modelfinput
            if (modelfinput.EQ."File") then
              READ(1,*) modelfinput
            if (modelfinput.EQ."Velocity") then
              READ(1,*) modelfinput
            if (modelfinput.EQ."Boundry") then
              boolthree = 1
              do while (boolthree.EQ.1)
                READ(1,*) modelfinput
                if (modelfinput.EQ."Temperature") then
                  READ(1,*) modelfinput
                else if (modelfinput.EQ."HeatFlux") then
                  READ(1,*) modelfinput
                else if (modelfinput.EQ."EndBoundry") then
                  boolthree = 0
                end if
              enddo
            else if (modelfinput.EQ."Region") then
              do while (boolthree.EQ.1)
                READ(1,*) num
                do i=1,num
                  READ(1,*) ID,part
                enddo
              enddo
            else if (modelfinput.EQ."EndDomain") then
              booltwo = 0
            end if
          enddo
        else if (modelfinput.EQ."Mineral") then
          do while (booltwo.EQ.1) then
            READ(1,*) modelfinput
            if (modelfinput.EQ."ID") then
              READ(1,*) ID
            else if (modelfinput.EQ."Density") then
              READ(1,*) typinput, modelfinput
              if (typinput.EQ."D") then
                READ(modelfinput,*) part
              end if
            else if (modelfinput.EQ."HeatProduction") then
              READ(1,*) typinput, modelfinput
              if (typinput.EQ."D") then
                READ(modelfinput,*) part
              end if
            else if (modelfinput.EQ."HeatCapcity") then
              READ(1,*) typinput, modelfinput
              if (typinput.EQ."D") then
                READ(modelfinput,*) part
              end if
            else if (modelfinput.EQ."ThermalConductivity") then
              READ(1,*) typinput, modelfinput
              if (typinput.EQ."D") then
                READ(modelfinput,*) part
              end if
            else if (modelfinput.EQ."BulkModulus") then
              READ(1,*) typinput, modelfinput
              if (typinput.EQ."D") then
                READ(modelfinput,*) part
              end if
            else if (modelfinput.EQ."GrainSize") then
              READ(1,*) typinput, modelfinput
              if (typinput.EQ."D") then
                READ(modelfinput,*) part
              end if
            else if (modelfinput.EQ."EndMineral") then
              booltwo = 0
            end if
          enddo
        else if (modelfinput.EQ."Output") then
          do while (booltwo.EQ.1)
            READ(1,*) modelfinput
            if (modelfinput.EQ."File") then
              READ(1,*) this%outfile
            else if (modelfinput.EQ."NegativeDown") then 
              this%negativedown = 1
            else if (modelfinput.EQ."EndOutput") then
              booltwo = 0
            end if
          enddo
        else if (modelfinput(1).EQ."!") then
          continue
        else if (modelfinput.EQ."End") then
          bool = 0
        end if
      enddo

      allocate ( this%density(this%tdim, this%ydim) )
      allocate ( this%heatproduction(this%tdim, this%ydim) )
      allocate ( this%heatcapcity(this%tdim, this%ydim) )
      allocate ( this%gtemp(this%tdim, this%ydim) )
      allocate ( this%gqflux(this%tdim, this%ydim) )
      allocate ( this%thermalconductivity(this%tdim, this%ydim) )
      allocate ( this%bulkmodulus(this%tdim, this%ydim) )
      allocate ( this%grainsize(this%tdim, this%ydim) )

      end subroutine
     
      subroutine new(this,n,m)
      type (modelfile) this
      integer n,m

      allocate ( this%velocity(n,m) )
      allocate ( this%density(n,m) )
      allocate ( this%heatproduction(n,m) )
      allocate ( this%heatcapcity(n,m) )
      allocate ( this%gtemp(n,m) )
      allocate ( this%gqflux(n,m) )
      allocate ( this%thermalconductivity(n,m) )
      allocate ( this%bulkmodulus(n,m) )
      allocate ( this%grainsize(n,m) )

      end subroutine
 
      subroutine delete(this) 
      type (modelfile) this
      
      deallocate ( this%velocity )
      deallocate ( this%density )
      deallocate ( this%heatproduction )
      deallocate ( this%heatcapcity )
      deallocate ( this%gtemp )
      deallocate ( this%gqflux )
      deallocate ( this%thermalconductivity )
      deallocate ( this%bulkmodulus )
      deallocate ( this%grainsize )

      end subroutine
     
      end module
