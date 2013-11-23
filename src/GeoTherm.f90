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
      type (mineralphase) eclogite
      double precision, allocatable :: eclogite_content(:,:)
      double precision, allocatable :: t_ar(:,:)
      integer n,m,i,j

      !gets the name of the file use to specify
      !model parameters, set as a sytem argument
      call GETARG(1,filename)

      !reads the file used to specify model parameters
      !and stores these in the instace modelfile_inst
      ! of derived type modelfile
      call READMDLF(modelfile_inst,filename,domain)
      
      !main work subroutine
      !call compute_geotherm &
      !&(modelfile_inst,domain,modelfile_inst%tdim,modelfile_inst%ydim)

      n = modelfile_inst%tdim
      m = modelfile_inst%ydim

      allocate( eclogite_content(n,m) )
      allocate( t_ar(n,m) )

      !call MODELDELETE(modelfile_inst)

      !end program

      !subroutine compute_geotherm(modelfile_inst,domain,n,m)
      !use equationpartsmodule
      !use module_modelfile, MODELDELETE => DELETE, &
      !MODELNEW => NEW 
      !use mathmodule
      !use helpermodule
      !use pressuresolver, PRESSUREFIELDNEW => NEW, &
      !&PRESSUREFIELDDELETE => DELETE
      !use geochem, MINERALDELETE => DELETE, &
      !& MINERALNEW => NEW
      !type (modelfile), intent(in) :: modelfile_inst
      !type (pressurefield) :: pressurefield_inst
      !type (temperaturefield) :: temperaturefield_inst
      !type (mineralphase) eclogite
      !type (modeldomain) domain
      !double precision :: eclogite_content(n,m)
      !double precision :: t_ar(n,m)
      !integer n,m,i,j

      call MINERALNEW( eclogite,n,m )
      write(2,*) "setting up pressure field"
      call PRESSUREFIELDNEW( pressurefield_inst,n,m )
      write(2,*) "setting up temperature field"
      call NEW( temperaturefield_inst,n,m )
      write(2,*) "updating domain"
      call UPDATE( domain )

      eclogite%free_energy = 326352.0
      eclogite%diffusion_coefficient = 0.00002

      !reading netcdfs from files sepcified in the modelfile
      !call get_netcdf1d(modelfile_inst%gtempfile, &
      !&modelfile_inst%gtemp,n,m)
      !call get_netcdf1d(modelfile_inst%gqfluxfile, &
      !&modelfile_inst%gqflux,n,m)
      !call get_netcdf_wincrmnt(modelfile_inst%velocitynetcdf, &
      !&modelfile_inst%velocity,modelfile_inst%incriment,n,m)
      !call get_netcdf(modelfile_inst%densitynetcdf, &
      !&modelfile_inst%density)
      !call get_netcdf(modelfile_inst%heatproductnetcdf, &
      !&modelfile_inst%heatproduction)
      !call get_netcdf(modelfile_inst%heatcapcnetcdf, &
      !&modelfile_inst%heatcapcity)
      !call get_netcdf &
      !&(modelfile_inst%thermlconductnetcdf, &
      !&modelfile_inst%thermalconductivity)
      !call get_netcdf &
      !&(modelfile_inst%incompresibilitynetcdf, &
      !&modelfile_inst%bulkmodulus)
      !call get_netcdf(modelfile_inst%grainsizenetcdf, &
      !&modelfile_inst%grainsize)

      !tempory
      !modelfile_inst%velocity = modelfile_inst%velocity/10.0
      
      write(2,*) "computing pressure"
      call compute_pressure(pressurefield_inst,domain)
      write(2,*) "computing density"
      call compute_pddensity(pressurefield_inst,domain)
      write(2,*) "computing temperature"
      call compute_temp(temperaturefield_inst,domain, &
      &pressurefield_inst)

      do j=1,m
        t_ar(:,j) = (/(i,i=1,n)/)*modelfile_inst%incriment(1)
      enddo
      
      call compute_eclogite_content(eclogite,t_ar, &
      &temperaturefield_inst%temperature, &
      &pressurefield_inst%pressure)
      
      !writes output netcdf
      write(2,*) "writing to output"
      call write_netcdf &
      &(modelfile_inst,temperaturefield_inst, &
      &pressurefield_inst,eclogite)

      write(2,*) "deleting temperature field"
      call DELETE(temperaturefield_inst)
      write(2,*) "deleting pressure field"
      call PRESSUREFIELDDELETE(pressurefield_inst)
      call MINERALDELETE(eclogite)
      write(2,*) "deleting domain"
      call DELETEDOMAIN( domain )
      write(2,*) "done"
      !end subroutine

      !call MODELDELETE(modelfile_inst)

      end program
