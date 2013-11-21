
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

      module modelregionmodule
      !use netcdf
      use geochem, MINERALDELETE => DELETE, &
      & NEWMINERAL => NEW
      implicit none

      type modelregion

        integer ID
        integer minerals
        integer, allocatable :: mineralids(:)
        double precision, allocatable :: mineralparts(:)

        integer, private :: indxmn = 1

      endtype

      contains

      subroutine NEW(this,minerals)
      type (modelregion) this
      integer minerals      

      this%minerals = minerals
      allocate( this%mineralids(minerals) )
      allocate( this%mineralparts(minerals) )

      end subroutine

      subroutine DELETE(this)
      type (modelregion) this

      deallocate( this%mineralids )
      deallocate( this%mineralparts )

      end subroutine

      subroutine addmineral(this,mineralid,part)
      type (modelregion) this
      integer mineralid
      double precision part

      this%mineralids(this%indxmn) = mineralid
      this%mineralparts(this%indxmn) = part

      this%indxmn = &
      this%indxmn + 1

      end subroutine

      end module



