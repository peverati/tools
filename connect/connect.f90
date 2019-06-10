! Copyright (c) 2018 
! Jose Luis Casals Sainz <jluiscasalssainz@gmail.com>,
! Ángel Martín Pendás <angel@fluor.quimica.uniovi.es> and 
! Evelio Francisco Miguelez <evelio@uniovi.es>. 
!
! connect is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or (at
! your option) any later version.
! 
! connect is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
program connect

  use mod_io, only: ioinit, stdargs, string, equal, isinteger, &
                    isreal, isword, getline, lgetword, ferror, &
                    faterr, uin, uout, mline
  use mod_geom, only: init_geom, readxyz, deallocate_space, cmcq, &
                      get_connect, covx
  implicit none                            

  character(len=:), allocatable :: optv !< command-line arguments 
  
  logical :: ok
  integer(kind=4) :: lp
  character(len=:), allocatable :: line, subline, word

  character(len=mline) :: vchar
  
  call ioinit ()
  call stdargs (optv)
  call init_geom ()

  write (uout,'(1x,a)') string('#  ===================================================')
  write (uout,'(1x,a)') string('# |     CONNECT: Analysis of molecular geometries     |')
  write (uout,'(1x,a)') string('# |  (c) E. Francisco, University of Oviedo, 2018     |')
  write (uout,'(1x,a)') string('# |  (c) A. Martin Pendas, University of Oviedo, 2018 |')
  write (uout,'(1x,a)') string('# |  (c) J. L. Casals Sainz                           |')
  write (uout,'(1x,a)') string('#  ===================================================')
  write (uout,'(1x,a)') string('#')

  write (uout,'(1x,a)') string('# +++ Begin to read input')
  do while (getline(uin,line))
    lp = 1
    word = lgetword(line,lp)
    subline = line(lp:)

    if (equal(word,'#')) then
      continue
    
    else if (equal(word,'covx')) then
      ok = isreal(covx, line, lp)
      if (.not.ok) call ferror('connect', 'wrong covx line', faterr) 
      write (uout,'(1x,a,1x,f6.3)') string('Covx changed to'), covx
    
    else if (equal(word,'load')) then
      ok = isword(vchar, line, lp)
      if (.not.ok) call ferror('connect', 'wrong load line', faterr) 
      call readxyz (vchar)
      exit 

    ! End of input
    else if (equal(word,'end')) then
      exit

    else
      call ferror('connect', 'unknown option', faterr)

    end if
  end do

  call deallocate_space ()

  write (uout,'(1x,a)') string('# +++ End of read input')

end program
