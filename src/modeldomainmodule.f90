
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
        integer minerals
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
        character(len=4) :: unitsvelocity = "MPS"
        double precision :: frameofrefrance = 0

        integer, private :: indxrg = 1

      end type

      contains
      
      subroutine NEW(this,regions,n,m)
      type (modeldomain) this
      integer regions
      integer n,m      
     
      this%n = n
      this%m = m 
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

      do i=1,this%minerals
        call MINERALDELETE( this%mineralarray(i) )
      enddo

      deallocate( this%mineralarray )

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
     
      !this is not getting called imediatly after region assignment so will have no meaning 
      call &
      &regionaddmineral(this%regionarray(this%indxrg-1),mineralid,part)

      end subroutine

      subroutine setminerals(this,in_mineralarray)
      type (modeldomain) this
      type (mineralphase) in_mineralarray(:)
 
      this%minerals = SIZE(in_mineralarray)
      allocate( this%mineralarray(this%minerals) )
      this%mineralarray = in_mineralarray

      end subroutine

      subroutine computemineralparts( this,t_ar,temperature,pressure )
      type (modeldomain) this
      integer i,j,minerals
      double precision :: t_ar(:,:)
      double precision :: temperature(:,:)
      double precision :: pressure(:,:)      

      !need to add in section on parent

      do i=1,this%minerals
       if (this%mineralarray(i)%parent.EQ.0) then
         call compute_part(this%mineralarray(i),this%mineralarray(i), &
         &t_ar,temperature,pressure)
       else
       do j=1,this%minerals
         if (j.EQ.i) then         
         continue
         end if
         if (this%mineralarray(j)%ID.EQ.this%mineralarray(i)%parent) then
         call compute_part(this%mineralarray(i),this%mineralarray(j), &
         &t_ar,temperature,pressure)
         exit
         end if
       enddo
       end if
      enddo
      
      end subroutine

      subroutine UPDATE(this)
      type (modeldomain) this
      integer i,j     
      double precision :: part

      this%density = this%geometry*0d0
      this%heatproduction = this%geometry*0d0
      this%heatcapcity = this%geometry*0d0
      this%thermalconductivity = this%geometry*0d0
      this%bulkmodulus = this%geometry*0d0
      this%grainsize = this%geometry*0d0

      !should really be done with pointers

      do i=1,this%regions
          do j=1,this%minerals
            !where &
          !&(this%mineralarray(j)%ID.EQ.this%regionarray(i)%mineralids)
              !part = this%regionarray(i)%mineralparts
            !elsewhere
            !  part = 0d0
            !end where
            if (ANY(this%mineralarray(j)%ID .EQ. &
            &this%regionarray(i)%mineralids)) then
              
            !else

            !endif
            
            part = MAXVAL(this%regionarray(i)%mineralparts,MASK = &
            &this%mineralarray(j)%ID .EQ. &
            &this%regionarray(i)%mineralids   )
            else
            part = 0
            endif


            !print *, "part ",part

            where (this%geometry.EQ.this%regionarray(i)%ID)

              this%density = this%density + &
              &this%mineralarray(j)%density*part
              this%heatproduction = this%heatproduction + &
              &this%mineralarray(j)%heatproduction*part
              this%heatcapcity = this%heatcapcity + &
              &this%mineralarray(j)%heatcapcity*part
              this%thermalconductivity = this%thermalconductivity + &
              &this%mineralarray(j)%thermalconductivity*part
              this%bulkmodulus = this%bulkmodulus + &
              &this%mineralarray(j)%bulkmodulus*part
              this%grainsize = this%grainsize + &
              &this%mineralarray(j)%grainsize*part

            end where
          enddo
      enddo

      end subroutine

      subroutine updateminerals(this)
      type (modeldomain) this
      integer i,j      
      double precision :: part

      !allocate( part(this%minerals) )
      do j=1,this%minerals
        this%mineralarray(j)%mineralpart = 0D0
      enddo

      !should really be done with pointers
      !probably gets overwritten by GeoChem subroutines
      do i=1,this%regions
          do j=1,this%minerals
            !where &
          !&(this%mineralarray(j)%ID.EQ.this%regionarray(i)%mineralids)
            !  part = this%regionarray(i)%mineralparts
            !elsewhere
            !  part = 0d0
            !end where
            !part = MAXVAL(this%regionarray(i)%mineralids,MASK = &
            !&this%mineralarray(j)%ID .EQ. &
            !&this%regionarray(i)%mineralids   )

            if (ANY(this%mineralarray(j)%ID .EQ. &
            &this%regionarray(i)%mineralids)) then

            !else

            !endif
            
            part = MAXVAL(this%regionarray(i)%mineralparts,MASK = &
            &this%mineralarray(j)%ID .EQ. &
            &this%regionarray(i)%mineralids   )
            else
            part = 0
            endif
            !print *, i"part ",part
            
            where (this%geometry.EQ.this%regionarray(i)%ID)
              this%mineralarray(j)%mineralpart = &
              &this%mineralarray(j)%mineralpart + part
            !elsewhere
            !  this%mineralarray(j)%mineralpart = 0D0
            end where
          enddo
      enddo    


      end subroutine

      subroutine offsetgeometry( this,n,m )
      type (modeldomain) this
      integer mask(n,m)      
      integer shift(n)
      integer baseid(n,m)
      double precision :: basevelo(n,m)
      integer i,n,m      
      !note domain should probably have it's own m,n
      where (this%geometry.EQ.0)
        mask = 1
      elsewhere
        mask = 0
      end where
       
      shift = SUM(mask,2)
     
      do i=1,m
        baseid(:,i) = this%geometry(:,1)
        basevelo(:,i) = this%velocity(:,1)
      enddo

      this%geometry = CSHIFT(this%geometry,-shift,2)
      this%velocity = CSHIFT(this%velocity,-shift,2)
      
      mask = 1D0
      
      where (this%geometry.EQ.0)
        this%geometry = baseid
        this%velocity = basevelo
      elsewhere
        this%geometry = this%geometry
        this%velocity = this%velocity
      end where
      
      this%geometry = this%geometry(:,m:1:-1)
      this%velocity = this%velocity(:,m:1:-1)

      end subroutine

      subroutine rescale( this )
      type (modeldomain) this
      double precision :: conversion = 0.01/(365.24*24*3600)  

      if (this%unitsvelocity.EQ."CMPA") then
        this%velocity = this%velocity*conversion
        this%frameofrefrance = this%frameofrefrance*conversion
      end if  

      if (this%frameofrefrance.NE.0D0) then
        this%incriment(1) = this%incriment(1)/this%frameofrefrance
      end if
      
      end subroutine

      end module
