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
      use helpermodule
      use geochem, MINERALDELETE => DELETE, &
      & NEWMINERAL => NEW 
      use modeldomainmodule, NEWDOMAIN => NEW, &
      & DELETEDOMAIN => DELETE
      use modellogfile, LOGNEW => NEW

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
      end type

      contains 

      !> subroutine to read file specifying
      !! model parameters
      !! @param this derived type storing
      !! file output
      !! @param filename name of file 
      !! specifying model parameters
      subroutine READMDLF(this,filename,domain)
      
      character (len = *) ::filename
      type (modelfile) this
      integer :: bool = 1
      integer :: booltwo = 1
      integer :: boolthree = 1      
      character (len = 1) :: typinput
      character (len = 256) :: modelfinput
      integer ID, num, i, nominerals
      double precision part
      type (mineralphase), allocatable :: minerals(:)
      type (mineralphase), allocatable :: mineralstempar(:)
      type (modeldomain) domain
      type (logfile) lfile
      character(len=10) :: netcdfx, netcdfy, netcdfz

      call LOGNEW(lfile,filename)

      call writelog(lfile,"geotherm logfile")
      call writelog(lfile,"reading file: " // filename)
      OPEN(1, file = filename)
      READ(1,*) this%tdim, this%ydim
      this%negativedown = 0      

      nominerals = 0

      allocate( minerals(nominerals) )

      do while (bool.EQ.1)
        booltwo = 1
        READ(1,*) modelfinput
        if (modelfinput.EQ."Domain") then
          READ(1,*) ID
          WRITE(modelfinput,*) ID
          call writelog(lfile,"creating domain: " // modelfinput)
          call NEWDOMAIN(domain,ID,this%tdim,this%ydim)
          do while (booltwo.EQ.1)
            READ(1,*) modelfinput
            if (modelfinput.EQ."FileNonXYZ") then
              READ(1,*) netcdfx, netcdfy, netcdfz
              READ(1,*) modelfinput
              call get_varnetcdf(modelfinput, &
              &domain%geometry, (/netcdfx, netcdfy, netcdfz/), &
              &this%incriment,this%tdim,this%ydim)              
            else if (modelfinput.EQ."File") then
              READ(1,*) modelfinput
              call writelog &
              &(lfile,"reading geometry netcdf: " // modelfinput)
              call get_idnetcdf(modelfinput, &
              &domain%geometry,this%incriment,this%tdim,this%ydim)
              domain%incriment = this%incriment
            else if (modelfinput.EQ."Velocity") then
              READ(1,*) modelfinput
              call writelog &
              &(lfile,"reading velocity netcdf: " // modelfinput)
              call get_netcdf(modelfinput, &
              &domain%velocity)
            else if (modelfinput.EQ."RefranceFrame") then
              READ(1,*) domain%frameofrefrance
            else if (modelfinput.EQ."VelocityUnits") then
              READ(1,*) domain%unitsvelocity
            else if (modelfinput.EQ."Boundry") then
              boolthree = 1
              do while (boolthree.EQ.1)
                READ(1,*) modelfinput
                if (modelfinput.EQ."Temperature") then
                  READ(1,*) typinput,modelfinput
                  if (typinput.EQ."D") then
                    READ(modelfinput,*) part
                    WRITE(modelfinput,*) part
                    call writelog &
                    &(lfile,"setting gtemp to: " // modelfinput)
                    domain%gtemp = part
                  else
                    call writelog &
                    &(lfile,"reading gtemp netcdf: " // modelfinput)
                    call get_netcdf1d(modelfinput, &
                    &domain%gtemp,this%tdim, this%ydim)
                  end if
                else if (modelfinput.EQ."HeatFlux") then
                  READ(1,*) typinput,modelfinput
                  if (typinput.EQ."D") then
                    READ(modelfinput,*) part
                    WRITE(modelfinput,*) part
                    call writelog &
                    &(lfile,"setting gqflux to: " // modelfinput)
                    domain%gqflux = part
                  else
                    call writelog &
                    &(lfile,"reading gqflux netcdf: " // modelfinput)
                    call get_netcdf1d(modelfinput, &
                    &domain%gqflux,this%tdim, this%ydim)
                  end if
                else if (modelfinput.EQ."EndBoundry") then
                  boolthree = 0
                end if
              enddo
            else if (modelfinput.EQ."Region") then
              READ(1,*) ID
              READ(1,*) num
              WRITE(modelfinput,*) ID 
              call writelog(lfile,"creating region: " // modelfinput)
              call addregion(domain,ID,num)
              do i=1,num
                READ(1,*) ID,part
                call addmineral(domain,ID,part)
              enddo
            else if (modelfinput.EQ."EndDomain") then
              booltwo = 0
            end if
          enddo
        else if (modelfinput.EQ."Mineral") then
          nominerals = nominerals + 1
          allocate( mineralstempar(nominerals)  )         
          mineralstempar(:nominerals-1) = minerals
          deallocate( minerals )
          allocate( minerals(nominerals)  )
          minerals = mineralstempar
          deallocate(mineralstempar)
          call writelog(lfile,"creating new mineral")
          call NEWMINERAL(minerals(nominerals),this%tdim,this%ydim)
          do while (booltwo.EQ.1)
            READ(1,*) modelfinput
            if (modelfinput.EQ."ID") then
              READ(1,*) ID
              WRITE(modelfinput,*) ID
              call writelog(lfile,"ID: " // modelfinput)
              minerals(nominerals)%ID = ID
            else if (modelfinput.EQ."Density") then
              READ(1,*) typinput, modelfinput
              if (typinput.EQ."D") then
                READ(modelfinput,*) part
                call writelog &
                &(lfile,"setting density to: " // modelfinput)
                minerals(nominerals)%density = part
              else
                call writelog &
                &(lfile,"reading density netcdf: " // modelfinput)
                call get_netcdf(modelfinput, &
                &minerals(nominerals)%density)
              end if
            else if (modelfinput.EQ."HeatProduction") then
              READ(1,*) typinput, modelfinput
              if (typinput.EQ."D") then
                READ(modelfinput,*) part
                call writelog &
                &(lfile,"setting heatproduction to: " // modelfinput)
                minerals(nominerals)%heatproduction = part
              else
                call writelog &
               &(lfile,"reading heatproduction netcdf: " // modelfinput)
                call get_netcdf(modelfinput, &
                &minerals(nominerals)%heatproduction)
              end if
            else if (modelfinput.EQ."HeatCapcity") then
              READ(1,*) typinput, modelfinput
              if (typinput.EQ."D") then
                READ(modelfinput,*) part
                call writelog &
                &(lfile,"setting heatcapcity to: " // modelfinput)
                minerals(nominerals)%heatcapcity = part
              else
                call writelog &
                &(lfile,"reading heatcapcity netcdf: " // modelfinput)
                call get_netcdf(modelfinput, &
                &minerals(nominerals)%heatcapcity)
              end if
            else if (modelfinput.EQ."ThermalConductivity") then
              READ(1,*) typinput, modelfinput
              if (typinput.EQ."D") then
                READ(modelfinput,*) part
                call writelog &
                &(lfile,"setting thermalconductivity to: "//modelfinput)
                minerals(nominerals)%thermalconductivity = part
              else
                call writelog &
          &(lfile,"reading thermalconductivity netcdf: " // modelfinput)
                call get_netcdf(modelfinput, &
                &minerals(nominerals)%thermalconductivity)
              end if
            else if (modelfinput.EQ."BulkModulus") then
              READ(1,*) typinput, modelfinput
              if (typinput.EQ."D") then
                READ(modelfinput,*) part
                call writelog &
                &(lfile,"setting bulkmodulus to: " // modelfinput)
                minerals(nominerals)%bulkmodulus = part
              else
                call writelog &
                &(lfile,"reading bulkmodulus netcdf: " // modelfinput)
                call get_netcdf(modelfinput, &
                &minerals(nominerals)%bulkmodulus)
              end if
            else if (modelfinput.EQ."GrainSize") then
              READ(1,*) typinput, modelfinput
              if (typinput.EQ."D") then
                READ(modelfinput,*) part
                call writelog &
                &(lfile,"setting grainsize to: " // modelfinput)
                minerals(nominerals)%grainsize = part
              else
                call writelog &
                &(lfile,"reading grainsize netcdf: " // modelfinput)
                call get_netcdf(modelfinput, &
                &minerals(nominerals)%grainsize)
              end if
            else if (modelfinput.EQ."EndMineral") then
              booltwo = 0
            end if
          enddo
        else if (modelfinput.EQ."Output") then
          call writelog(lfile,"setting output")
          do while (booltwo.EQ.1)
            READ(1,*) modelfinput
            if (modelfinput.EQ."File") then
              READ(1,*) this%outfile
              call writelog(lfile,"ouput to: " // this%outfile)
            else if (modelfinput.EQ."NegativeDown") then 
              call writelog(lfile,"negative down true")
              this%negativedown = 1
            else if (modelfinput.EQ."EndOutput") then
              booltwo = 0
            end if
          enddo
        else if (modelfinput(1:2).EQ."!") then
          continue
        else if (modelfinput.EQ."End") then
          bool = 0
        end if
      enddo

      call writelog(lfile,"associating minerals to domain")
      call setminerals(domain,minerals)
      WRITE(2,*) "finished parsing: "      

      CLOSE(1) 

      end subroutine

      end module
