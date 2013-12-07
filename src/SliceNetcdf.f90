
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

      program SLICENETCDF
      use helpermodule
      implicit none

      character (len = 256) :: filename
      character (len = 256) :: outfile
      character (len = 10) :: variablename
      character (len = 8) :: char_n,char_m
      character (len = 1) :: char_ardim
      character (len = 8) :: char_possition
      double precision, allocatable :: data_input(:,:) 
      double precision, allocatable :: data_ar(:)
      
      integer :: n,m
      integer :: ardim, possition
      
      OPEN(2,FILE = "logfile.log")

      call GETARG(1,char_n)
      call GETARG(2,char_m)
      call GETARG(3,variablename)
      call GETARG(4,char_ardim)
      call GETARG(5,char_possition)
      call GETARG(6,filename)
      call GETARG(7,outfile)
      
      READ(char_n,*) n
      READ(char_m,*) m
      READ(char_ardim,*) ardim
      READ(char_possition,*) possition
      
      allocate( data_input(n,m) )
      
      if (ardim.EQ.1) then
        allocate( data_ar(m) )
      else if (ardim.EQ.2) then
       allocate( data_ar(n) )
      end if
      
      call get_vardoublenetcdf &
      &(filename,data_input,variablename)
      
      if (ardim.EQ.1) then
        data_ar = data_input(possition,:)
      else if (ardim.EQ.2) then
        data_ar = data_input(:,possition)
      end if
      
      call write_1dslice(outfile,data_ar)

      deallocate( data_input )
      deallocate( data_ar ) 

      end program  
