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

      program GEOTHERM
      use module_modelfile
      use modeloutput
      use equationpartsmodule
      use mathmodule
      use helpermodule
      use pressuresolver, PRESSUREFIELDNEW => NEW, &
      &PRESSUREFIELDDELETE => DELETE
      use geochem, MINERALDELETE => DELETE, &
      & MINERALNEW => NEW
      implicit none

      type (modelfile) modelfile_inst
      type (modeldomain) domain
      character (len = 256) :: filename
      type (pressurefield) :: pressurefield_inst
      type (temperaturefield) :: temperaturefield_inst
      double precision, allocatable :: t_ar(:,:)
      integer n,m,i,j

      !gets the name of the file use to specify
      !model parameters, set as a system argument
      call GETARG(1,filename)

      !reads the file used to specify model parameters
      !and stores these in the instace modelfile_inst
      ! of derived type modelfile
      call READMDLF(modelfile_inst,filename,domain)      

      n = modelfile_inst%tdim
      m = modelfile_inst%ydim

      allocate( t_ar(n,m) )

      !allocates memory for derived types
      write(2,*) "setting up pressure field"
      call PRESSUREFIELDNEW( pressurefield_inst,n,m )
      write(2,*) "setting up temperature field"
      call NEW( temperaturefield_inst,n,m )
      write(2,*) "offsetting geometry"
      !this is buggered
      call offsetgeometry( domain,n,m )
      write(2,*) "updating domain"
      call UPDATE( domain )
      write(2,*) "rescaling coordinates"
      call rescale( domain )
      print *,12
      call updateminerals( domain )

      !computes fields
      write(2,*) "computing pressure"
      call compute_pressure(pressurefield_inst,domain)
      write(2,*) "computing density"

      call compute_pddensity(pressurefield_inst,domain)
      write(2,*) "computing temperature"

      call compute_temp(temperaturefield_inst,domain, &
      &pressurefield_inst)
 
      !set up time array
      do j=1,m
        t_ar(:,j) = (/(i,i=1,n)/)*modelfile_inst%incriment(1)
      enddo
      
      !computes mineral reactions
      call computemineralparts( domain,t_ar, &
      &temperaturefield_inst%temperature,pressurefield_inst%pressure )
      
      !writes output netcdf
      write(2,*) "writing to output"
      call write_netcdf &
      &(modelfile_inst,temperaturefield_inst, &
      &pressurefield_inst,domain)
      
      !deallocating memory for derived types
      write(2,*) "deleting temperature field"
      call DELETE(temperaturefield_inst)
      write(2,*) "deleting pressure field"
      call PRESSUREFIELDDELETE(pressurefield_inst)
      write(2,*) "deleting domain"
      call DELETEDOMAIN( domain )
      write(2,*) "done"
      
      CLOSE(2)

      end program
