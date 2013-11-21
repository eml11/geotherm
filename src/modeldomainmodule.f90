
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

      module modeldomainmodule
      !use netcdf
      use geochem, MINERALDELETE => DELETE, &
      & NEWMINERAL => NEW
      use modelregionmodule, NEWREGION => NEW, &
      & DELETEREGION => DELETE, regionaddmineral => addmineral
      implicit none

      type modeldomain
        integer :: negativedown
        integer :: n,m
        integer :: regions
        integer, allocatable :: geometry(:,:)

        double precision :: incriment(2)

        double precision, allocatable :: gtemp(:,:)
        double precision, allocatable :: gqflux(:,:)
        double precision, allocatable :: velocity(:,:)
        double precision, allocatable :: density(:,:)
        double precision, allocatable :: heatproduction(:,:)
        double precision, allocatable :: heatcapcity(:,:)
        double precision, allocatable :: thermalconductivity(:,:)
        double precision, allocatable :: bulkmodulus(:,:)
        double precision, allocatable :: grainsize(:,:)
        !double precision, allocatable, target :: IDS(:)
        !double precision, pointer :: regionpointers(:,:)
        type (modelregion), allocatable :: regionarray(:)      
        type (mineralphase), allocatable :: mineralarray(:) 

        integer, private :: indxrg = 1

      end type

      contains
      
      subroutine NEW(this,regions,incriment,n,m)
      type (modeldomain) this
      integer regions
      integer n,m      
      double precision :: incriment(2)

      this%incriment = incriment
      this%regions = regions
      allocate( this%geometry(n,m) )      

      allocate( this%gtemp(n,m) ) 
      allocate( this%gqflux(n,m) )

      allocate( this%velocity(n,m) )
      allocate( this%density(n,m) )
      allocate( this%heatproduction(n,m) )
      allocate( this%heatcapcity(n,m) )
      allocate( this%thermalconductivity(n,m) )
      allocate( this%bulkmodulus(n,m) )
      allocate( this%grainsize(n,m) )

      allocate( this%regionarray(regions) )

      end subroutine

      subroutine DELETE(this)
      type (modeldomain) this
      integer i

      deallocate( this%gtemp )
      deallocate( this%gqflux )

      deallocate( this%velocity )
      deallocate( this%density )
      deallocate( this%heatproduction )
      deallocate( this%heatcapcity )
      deallocate( this%thermalconductivity )
      deallocate( this%bulkmodulus )
      deallocate( this%grainsize )

      do i=1,this%regions
        call DELETEREGION( this%regionarray(i) )
      enddo

      deallocate( this%regionarray )

      end subroutine

      subroutine addregion(this,inid,arraysize)
      type (modeldomain) this
      type (modelregion) region
      integer, target :: inid
      integer arraysize
      
      call NEWREGION(region,arraysize)
      region%ID = inid
      this%regionarray(this%indxrg) = region

      this%indxrg=this%indxrg+1

      end subroutine

      subroutine addmineral(this,mineralid,part)
      type (modeldomain) this
      integer mineralid
      double precision part
      
      call &
      &regionaddmineral(this%regionarray(this%indxrg-1),mineralid,part)

      end subroutine

      subroutine setminerals(this,in_mineralarray)
      type (modeldomain) this
      type (mineralphase) in_mineralarray(:)
      integer minerals
       
      minerals = SIZE(in_mineralarray)
      allocate( this%mineralarray(minerals) )
      this%mineralarray = in_mineralarray

      end subroutine

      subroutine UPDATE(this)
      type (modeldomain) this
      integer i,j   
      integer minerals   
      double precision, allocatable :: part(:)

      minerals = SIZE(this%mineralarray)
      allocate( part(minerals) )

      this%density = this%geometry*0d0
      this%heatproduction = this%geometry*0d0
      this%heatcapcity = this%geometry*0d0
      this%thermalconductivity = this%geometry*0d0
      this%bulkmodulus = this%geometry*0d0
      this%grainsize = this%geometry*0d0

      !should really be done with pointers

      do i=1,this%regions
          do j=1,minerals
            where &
          &(this%mineralarray(j)%ID.EQ.this%regionarray(i)%mineralids)
              part = this%regionarray(i)%mineralparts
            elsewhere
              part = 0d0
            end where

            where (this%geometry.EQ.this%regionarray(i)%ID)

              this%density = this%density + &
              &this%mineralarray(j)%density*part(1) 
              this%heatproduction = this%density + &
              &this%mineralarray(j)%density*part(1)
              this%heatcapcity = this%density + &
              &this%mineralarray(j)%density*part(1)
              this%thermalconductivity = this%density + &
              &this%mineralarray(j)%density*part(1)
              this%bulkmodulus = this%density + &
              &this%mineralarray(j)%density*part(1)
              this%grainsize = this%density + &
              &this%mineralarray(j)%density*part(1)

            end where
          enddo
      enddo

      end subroutine

      end module
