! Copyright (c) 2018 
! Jose Luis Casals Sainz <jluiscasalssainz@gmail.com>,
! Ángel Martín Pendás <angel@fluor.quimica.uniovi.es> and 
! Evelio Francisco Miguelez <evelio@uniovi.es>. 
!
! prejmol is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or (at
! your option) any later version.
! 
! prejmol is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
! This program reads an AIMPAC-like wavefunction file, prepares 
! the input of the graphical program Jmol in a file named 
! having the format of MOLDEN
!
program prejmol

  use mod_io, only: mline, ioinit, stdargs
  use mod_param, only: init_param
  use mod_wfn, only: init_wfn, maxorb, rdwfn
  implicit none

  character(len=:), allocatable :: optv !< command-line arguments 
  character(len=:), allocatable :: filedat !< file prefix
  integer(kind=4) :: ichange(maxorb+maxorb)

  call ioinit ()
  call stdargs (optv, filedat)

  call init_param ()
  call init_wfn ()

  ichange(1:maxorb+maxorb) = +1
  call rdwfn (ichange,trim(filedat))

end program prejmol
