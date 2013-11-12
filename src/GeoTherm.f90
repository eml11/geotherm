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
      use helpermodule
      use mathmodule
      use module_modelfile
      implicit none

      type (modelfile) modelfile_inst
      character (len = 256) :: filename

      !gets the name of the file use to specify
      !model parameters, set as a sytem argument
      call GETARG(1,filename)

      !reads the file used to specify model parameters
      !and stores these in the instace modelfile_inst
      ! of derived type modelfile
      call READMDLF(modelfile_inst,filename)
      
      !main work subroutine
      call compute_geotherm &
      &(modelfile_inst,modelfile_inst%tdim,modelfile_inst%ydim)

      call DELETE(modelfile_inst)

      end program

      subroutine compute_geotherm(modelfile_inst,n,m)
      use equationpartsmodule
      use module_modelfile
      use mathmodule
      use helpermodule
      use pressuresolver, PRESSUREFIELDNEW => NEW, &
      &PRESSUREFIELDDELETE => DELETE
      use geochem
      type (modelfile), intent(in) :: modelfile_inst
      type (pressurefield) :: pressurefield_inst
      type (temperaturefield) :: temperaturefield_inst
      double precision :: eclogite_content(n,m)
      double precision :: t_ar(n,m)
      integer n,m,i,j

      call PRESSUREFIELDNEW( pressurefield_inst,n,m )
      call NEW( temperaturefield_inst,n,m )

      !reading netcdfs from files sepcified in the modelfile
      call get_netcdf1d(modelfile_inst%gtempfile, &
      &modelfile_inst%gtemp,n,m)
      call get_netcdf1d(modelfile_inst%gqfluxfile, &
      &modelfile_inst%gqflux,n,m)
      call get_netcdf_wincrmnt(modelfile_inst%velocitynetcdf, &
      &modelfile_inst%velocity,modelfile_inst%incriment,n,m)
      call get_netcdf(modelfile_inst%densitynetcdf, &
      &modelfile_inst%density)
      call get_netcdf(modelfile_inst%heatproductnetcdf, &
      &modelfile_inst%heatproduction)
      call get_netcdf(modelfile_inst%heatcapcnetcdf, &
      &modelfile_inst%heatcapcity)
      call get_netcdf &
      &(modelfile_inst%thermlconductnetcdf, &
      &modelfile_inst%thermalconductivity)
      call get_netcdf &
      &(modelfile_inst%incompresibilitynetcdf, &
      &modelfile_inst%bulkmodulus)
      call get_netcdf(modelfile_inst%grainsizenetcdf, &
      &modelfile_inst%grainsize)

      !tempory
      !velocity_ar = velocity_ar*0
      
      call compute_pressure(pressurefield_inst,modelfile_inst)
      
      call compute_pddensity(pressurefield_inst,modelfile_inst)

      call compute_temp(temperaturefield_inst,modelfile_inst,pressurefield_inst)

      do j=1,m
        t_ar(:,j) = (/(i,i=1,n)/)*modelfile_inst%incriment(1)
      enddo
      
      call compute_eclogite_content &
      &(t_ar,temperaturefield_inst%temperature, &
      &pressurefield_inst%pressure, &
      &modelfile_inst%grainsize,eclogite_content,n,m)
      
      !writes output netcdf
      call write_netcdf &
      &(modelfile_inst%outfile,temperaturefield_inst%temperature, &
      &pressurefield_inst%density,pressurefield_inst%pressure, &
      & eclogite_content,modelfile_inst%negativedown,n,m)

      DELETE(temperaturefield_inst)
      PRESSUREFIELDDELETE(pressurefield_inst)

      end subroutine
