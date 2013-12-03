
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

      module modeloutput
      use netcdf
      use mathmodule
      use module_modelfile
      use equationpartsmodule
      use pressuresolver, PRESSURENEW => NEW, &
      & PRESSUREDELETE => DELETE
      use geochem, MINERALNEW => NEW, &
      & MINERALDELETE => DELETE
      implicit none

      contains

      !> netcdf writing subroutine
      !! @param filename name of netcdf file
      !! @param data_ar z data of netcdf 
      subroutine write_netcdf(model,tfield,pfield,domain)
      !&(filename,temp_ar,density_ar,pressure_ar,ecologite_ar, &
      !&negativedown,n,m)
      !type (modelfile) model
      type (modelfile) model
      type (modeldomain) domain
      type (temperaturefield) tfield
      type (pressurefield) pfield

      integer :: ncid, vvarid,tvarid, pvarid, dvarid, gvarid, dimids(2)
      integer :: mineralvarid(domain%minerals)
      integer :: x_dimid, y_dimid
      integer i      

      print *,2

      call wcheck( nf90_create(model%outfile, NF90_CLOBBER, ncid) )

      call wcheck( nf90_def_dim(ncid, "x", tfield%n, x_dimid) )
      call wcheck( nf90_def_dim(ncid, "y", tfield%m, y_dimid) )
      
      dimids =  (/ x_dimid, y_dimid /)
      print *,3
      call wcheck( nf90_def_var &
      &(ncid, "Temperature", NF90_DOUBLE, dimids, tvarid) )

      call wcheck( nf90_def_var &
      &(ncid, "Pressure", NF90_DOUBLE, dimids, pvarid) )

      call wcheck( nf90_def_var &
      &(ncid, "Density", NF90_DOUBLE, dimids, dvarid) )

      call wcheck( nf90_def_var &
      &(ncid, "Velocity", NF90_DOUBLE, dimids, vvarid) )

      !call wcheck ( nf90_def_var &
      !&(ncid, "EclogitePart", NF90_DOUBLE, dimids, evarid) )

      call wcheck ( nf90_def_var &
      &(ncid, "Geometry", NF90_DOUBLE, dimids, gvarid) )
      print *,4
      do i=1,domain%minerals
        call wcheck ( nf90_def_var &
        &(ncid, domain%mineralarray(i)%mineralname, NF90_DOUBLE, &
        &dimids, mineralvarid(i)) )
      enddo
      print *,5
      call wcheck( nf90_enddef(ncid) )

      if (model%negativedown.EQ.0) then
        call wcheck( nf90_put_var(ncid, dvarid, pfield%density) )
        call wcheck( nf90_put_var(ncid, tvarid, tfield%temperature) )
        call wcheck( nf90_put_var(ncid, pvarid, pfield%pressure) )
        !call wcheck( nf90_put_var(ncid, evarid, minphase%mineralpart) )
        call wcheck( nf90_put_var(ncid, gvarid, domain%geometry) )
        call wcheck( nf90_put_var(ncid, vvarid, domain%velocity) )
        do i=1,domain%minerals
          call wcheck( nf90_put_var(ncid, mineralvarid(i), &
          &domain%mineralarray(i)%mineralpart) )
        enddo
      else 
        call wcheck( &
        &nf90_put_var(ncid, dvarid, pfield%density(:,tfield%m:1:-1)) )
        call wcheck( &
        &nf90_put_var(ncid, tvarid,tfield%temperature(:,tfield%m:1:-1)))
        call wcheck( &
        &nf90_put_var(ncid, pvarid, pfield%pressure(:,tfield%m:1:-1)) )
        !call wcheck( &
        !&nf90_put_var(ncid,evarid, &
        !&minphase%mineralpart(:,tfield%m:1:-1)))
        call wcheck( &
        &nf90_put_var(ncid, gvarid, domain%geometry(:,tfield%m:1:-1)) )
        call wcheck( &
        &nf90_put_var(ncid, vvarid, domain%velocity(:,tfield%m:1:-1)) )
        do i=1,domain%minerals
          call wcheck( nf90_put_var(ncid, mineralvarid(i), &
          &domain%mineralarray(i)%mineralpart(:,tfield%m:1:-1)) )
        enddo
      end if
      print *,6
      end subroutine
      
      !> netcdf error checking subroutine
      !! @param status return value of function
      !! from netcdf-fortran
      subroutine wcheck(status)
      integer, intent ( in) :: status
    
      write(2,*) trim(nf90_strerror(status))
      if(status.NE.nf90_noerr) then
        stop "Stopped"
      end if
      end subroutine

      end module
