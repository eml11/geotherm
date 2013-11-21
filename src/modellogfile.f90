
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

      module modellogfile
      implicit none

      type logfile

      character(len=264) :: filename
      integer :: tofile = 1

      end type

      contains

      subroutine NEW(this,filename)
      type (logfile) this
      character(len=264) :: filename
      integer :: n
      
      !for now ignoring logfilename
      n = len(filename)
      !print *, len(TRIM(filename))
      !print *, filename(:len(TRIM(filename))-4) // "_logfile.log"
      !this%filename = filename(:len(TRIM(filename))-4) // "_logfile.log"
      this%filename = "./logfile.log"
      !print *, this%filename
      OPEN(2,FILE=this%filename)

      end subroutine

      subroutine writelog(this,string)
      type (logfile) this
      character(len=*) :: string

      if (this%tofile.EQ.1) then

        !OPEN(1, file = this%filename)
        WRITE(2,*) string

      end if
      end subroutine

      end module
