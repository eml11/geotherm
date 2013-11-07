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

      end program

      subroutine compute_geotherm &
      &(modelfile_inst,n,m)
      use equationpartsmodule
      use module_modelfile
      use mathmodule
      use helpermodule
      use pressuresolver
      use geochem
      type (modelfile), intent(in) :: modelfile_inst
      double precision :: input_tdata_ar(n)
      double precision :: input_qdata_ar(n)
      double precision :: input_qxdata(n)
      double precision :: tdata_ar(n,m)
      double precision :: qdata_ar(n,m), qxdata_ar(n,m)
      double precision :: bdash_ar(n,m)
      double precision :: velocity_ar(n,m)
      double precision :: density_ar(n,m)
      double precision :: thermlconduct_ar(n,m)
      double precision :: heatproduct_ar(n,m)
      double precision :: heatcapc_ar(n,m)
      double precision :: kappa_ar(n,m)
      double precision :: pddensity(n,m)
      double precision :: pressure(n,m)
      double precision :: exintegral_ar(n,m)
      double precision :: inerintegral_ar(n,m)
      double precision :: outerintegral_ar(n,m)
      double precision :: incompresibility_ar(n,m)
      double precision :: eclogite_content(n,m)
      double precision :: temperature(n,m)
      double precision :: diffusion_coeficient(n,m)
      double precision :: grain_size(n,m)
      double precision :: incriment(2)
      double precision :: iner_dbl
      double precision :: outr_dbl
      double precision :: t_ar(n,m)
      integer n,m,i,j


      !reading netcdfs from files sepcified in the modelfile
      call get_netcdf1d(modelfile_inst%gtempfile,input_tdata_ar)
      call get_netcdf1d(modelfile_inst%gqfluxfile,input_qdata_ar)
      call get_netcdf1d(modelfile_inst%gqxfluxfile,input_qxdata)
      call get_netcdf_wincrmnt &
      &(modelfile_inst%velocitynetcdf,velocity_ar,incriment,n,m)
      call get_netcdf(modelfile_inst%densitynetcdf,density_ar)
      call get_netcdf(modelfile_inst%heatproductnetcdf,heatproduct_ar)
      call get_netcdf(modelfile_inst%heatcapcnetcdf,heatcapc_ar)
      call get_netcdf &
      &(modelfile_inst%thermlconductnetcdf,thermlconduct_ar)
      call get_netcdf &
      &(modelfile_inst%incompresibilitynetcdf,incompresibility_ar)
      call get_netcdf(modelfile_inst%grainsizenetcdf,grain_size)

      !this is a hack to fix issue with netcdfs created by gmt 
      incriment(2) = -incriment(2)

      !extend 1d netcdfs into 2d by copying along 2d dimension
      call extend_ardimension(input_tdata_ar,tdata_ar,m)
      call extend_ardimension(input_qdata_ar,qdata_ar,m)

      call compute_bdashval_fromqx &
      &(tdata_ar,qdata_ar,qxdata_ar,thermlconduct_ar, &
      &bdash_ar,incriment,n,m)

      !creates array corrisponding to b(t)'s time differential
      !call compute_bdashval &
      !&(tdata_ar,qdata_ar,thermlconduct_ar, &
      !&bdash_ar,incriment,n,m)
      
      call compute_pressure &
      &(density_ar,incompresibility_ar,pressure,incriment,n,m)
      
      call compute_pddensity &
      &(density_ar,incompresibility_ar,pddensity,incriment,n,m)

      kappa_ar = thermlconduct_ar/(pddensity*heatcapc_ar) 

      !creates array corrisponding to the function which appears
      !in the exponent
      call compute_exponentintegral &
      &(bdash_ar,velocity_ar,kappa_ar,exintegral_ar,incriment,n,m)
      
      !creates array corrisponding to the first integral with respect
      !to eta
      call compute_init_inerintegral &
      &(exintegral_ar,bdash_ar,density_ar*heatproduct_ar,&
      &thermlconduct_ar,inerintegral_ar,incriment,n,m)
 
      !creates a double which is the integration constant of the
      !integral calculated in the above subroutine required to set
      !the surface heatflux equal to the qdata_ar
      call compute_inerintegralconstant &
      &(inerintegral_ar,exintegral_ar,tdata_ar,qdata_ar, &
      &bdash_ar,thermlconduct_ar,iner_dbl,incriment,n,m)

      !creates array corrisponding to the second integral with respect
      !to eta
      call compute_init_outerintegral &
      &(inerintegral_ar,exintegral_ar,bdash_ar,iner_dbl, &
      &outerintegral_ar,incriment,n,m)

      !creates double which is the integration constant of the
      !integral calcultaed in the above subroutine required to set
      !the surface Temperature to tdata_ar
      call compute_outerintegralconstant &
      &(outerintegral_ar,tdata_ar,outr_dbl,n,m)
      
      temperature = outerintegral_ar+outr_dbl

      !diffusion_coeficient?
      do i=1,n
        do j=1,m
          t_ar = i*incriment(1)
        enddo
      enddo

      call compute_eclogite_content &
      &(t_ar,temperature,pressure, &
      &grain_size,eclogite_content,n,m)

      !writes output netcdf
      call write_netcdf &
      &(modelfile_inst%outfile,temperature, &
      &pddensity,pressure,eclogite_content, &
      &modelfile_inst%negativedown,n,m)

      end subroutine
