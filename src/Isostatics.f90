
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

      program ISOSTATICS
      use module_modelfile
      use modeloutput
      use equationpartsmodule
      use mathmodule
      use helpermodule
      use pressuresolver, PRESSUREFIELDNEW => NEW, &
      &PRESSUREFIELDDELETE => DELETE
      use geochem, MINERALDELETE => DELETE, &
      & MINERALNEW => NEW
      use isostatichelper
      implicit none

      type (modelfile) modelfile_inst
      type (modeldomain) domain
      character (len = 256) :: filename
      character (len = 256) :: outfile
      type (pressurefield) :: pressurefield_inst
      type (temperaturefield) :: temperaturefield_inst
      double precision, allocatable :: density(:,:)
      double precision, allocatable :: deltay(:)
      double precision, allocatable :: integral1(:,:)
      double precision, allocatable :: integral2(:,:)
      double precision, allocatable :: start_density(:,:)
      integer n,m,i,j
      integer retval,ncid,varid

      call GETARG(1,filename)
      call GETARG(2,outfile)

      call READMDLF(modelfile_inst,filename,domain)      

      n = modelfile_inst%tdim
      m = modelfile_inst%ydim

      call PRESSUREFIELDNEW(pressurefield_inst,n,m)

      allocate( density(n,m) )
      allocate( deltay(n) )
      allocate( integral1(n,m) )
      allocate( integral2(n,m) )
      allocate( start_density(n,m) )

      write(2,*) "offsetting geometry"
      call offsetgeometry( domain,n,m )
      write(2,*) "rescaling coordinates"
      call rescale( domain )
      write(2,*) "updating domain"
      call UPDATE( domain )
      call updateminerals( domain ) 
 
      call get_basedensity(modelfile_inst%outfile,density,domain)

      retval = nf90_open(modelfile_inst%outfile, NF90_NOWRITE, ncid) 
      retval = nf90_inq_varid(ncid, "Density", varid)
      retval = nf90_get_var(ncid, varid, start_density)
      retval = nf90_close(ncid)

      !print *, density 

      !density = (start_density/domain%density)*density
      domain%density = density 
      call compute_pddensity(pressurefield_inst,domain)
      density = pressurefield_inst%density 
      !print *, density
      !print *, density

      call domaindensityintegral(start_density,domain,integral1)
      call mineraldensityintegral(density,domain,integral2)

      call compute_deltay(domain, &
      &integral1(:,m),integral2(:,m),deltay)

      call write_output(outfile,deltay)
      
      write(2,*) "deleting domain"
      call DELETEDOMAIN( domain )
      call PRESSUREFIELDDELETE( pressurefield_inst )
      write(2,*) "done"
      CLOSE(2)

      end program  
