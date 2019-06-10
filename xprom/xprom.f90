! Copyright (c) 2018 
! Jose Luis Casals Sainz <jluiscasalssainz@gmail.com>,
! Ángel Martín Pendás <angel@fluor.quimica.uniovi.es> and 
! Evelio Francisco Miguelez <evelio@uniovi.es>. 
!
! xprom is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or (at
! your option) any later version.
! 
! xprom is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
! This program reads an output file of xprom (POUT) after doing 
! a TES calculation of a given molecule and groups together different
! energetic quantities of this molecule to obtain several useful
! IQA energy properties of a group of atoms of this molecule. The
! atoms that define the group are read in interactively. TO RUN this
! program it is mandatory that the linux script 'xprom.csh' is in 
! the same directory than the output file of xprom. Optionally,
! the total energies of the reference state can be used.
!
program xprom

  use mod_io, only: ioinit, stdargs, string, equal, isinteger, &
                    isreal, isword, getline, lgetword, ferror, &
                    faterr, uin, uout, mline
  use mod_xprom, only: init_xprom, end_xprom, nodef, t2, loadpmd, &
                       full, resume
  implicit none                            

  character(len=:), allocatable :: optv !< command-line arguments 
  character(len=:), allocatable :: fileroot !< file prefix
  
  logical :: ok, file_exists
  integer(kind=4) :: lp
  character(len=:), allocatable :: line, subline, word

  character(len=mline) :: vchar
  
  call ioinit ()
  call stdargs (optv, fileroot)
  call init_xprom ()

  write (uout,'(1x,a)') string('#  ===================================================')
  write (uout,'(1x,a)') string('# |         XPROM: Parsing promolden outputs          |')
  write (uout,'(1x,a)') string('# |  (c) E. Francisco, University of Oviedo, 2018     |')
  write (uout,'(1x,a)') string('# |  (c) A. Martin Pendas, University of Oviedo, 2018 |')
  write (uout,'(1x,a)') string('# |  (c) J. L. Casals Sainz                           |')
  write (uout,'(1x,a)') string('#  ===================================================')
  write (uout,'(1x,a)') string('#')

  do while (getline(uin,line))
    lp = 1
    word = lgetword(line,lp)
    subline = line(lp:)

    if (equal(word,'#')) then
      continue
    
    else if (equal(word,'group')) then
      call pair ()
    
    else if (equal(word,'all')) then
      call full ()
    
    else if (equal(word,'t2')) then
      t2 = .true.
      write (uout,'(1x,a)') string('# Split correlation ON')

    else if (equal(word,'refstates')) then
      ok = isword(vchar, line, lp)
      if (.not.ok) call ferror('xprom', 'wrong refstates line', faterr) 
      !call loadref(vchar)
      nodef = .false.

    else if (equal(word,'load')) then
      ok = isword(vchar, line, lp)
      if (.not.ok) call ferror('xprom', 'wrong load line', faterr) 
      if (t2) then
        inquire (file="./xprom-t2.csh", exist=file_exists) 
        if (.not.file_exists) then
          call ferror ('xprom', 'file xprom-t2.csh does not exists', faterr)
        end if
        call system ('./xprom-t2.csh '//trim(vchar))
      else 
        inquire (file="./xprom.csh", exist=file_exists) 
        if (.not.file_exists) then
          call ferror ('xprom', 'file xprom.csh does not exists', faterr)
        end if
        call system ('./xprom.csh '//trim(vchar))
      end if
      call loadpmd ()
      call resume ()

    ! End of input
    else if (equal(word,'end')) then
      exit

    ! Exit main driver
    else if (equal(word,'exit') .or. equal(word,'quit')) then
      goto 1

    else
      call ferror('xprom', 'unknown option', faterr)

    end if
  end do

1 call end_xprom ()

end program
