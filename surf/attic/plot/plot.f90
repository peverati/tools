! Copyright (c) 2018 
! Jose Luis Casals Sainz <jluiscasalssainz@gmail.com>,
! Ángel Martín Pendás <angel@fluor.quimica.uniovi.es> and 
! Evelio Francisco Miguelez <evelio@uniovi.es>. 
!
! plot is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or (at
! your option) any later version.
! 
! plot is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
program plot

  !$ use omp_lib
  use mod_prec, only: rp, ip
  use mod_param, only: verbose, debug, init_param
  use mod_plot, only: init_plot, parse_plot
  use mod_io, only: ferror, faterr, getdate, stdargs, equal, &
                    isinteger, isreal, isword, ioinit, getline,  &
                    uin, mline, string, flush_unit, uout, &
                    nwarns, ncomms, warning, lgetword
  implicit none
 
  character(len=:), allocatable :: optv !< command-line arguments 
  character(len=:), allocatable :: fileroot !< file prefix
 
  integer(kind=ip) :: lp
  real(kind=rp) :: tiempo1, tiempo2
  character(len=mline) :: wdate  
  character(len=:), allocatable :: line, subline, word

  ! Begin program
  call ioinit ()
  call stdargs (optv, fileroot)
  call getdate (wdate)
  write (uout,'(1x,a)') string('#  =====================================================')
  write (uout,'(1x,a)') string('# |       plot: procesing tool for surface files        |')
  write (uout,'(1x,a)') string('#  =====================================================')
  write (uout,'(1x,a)') string('#')
  write (uout,'(1x,a)') string('# Calculation starts at '//wdate)
  write (uout,'(1x,a)') string('# *** Begin to initialize data')
  call init_param ()
  call init_plot ()
  write (uout,'(1x,a)') string('# *** End of initialize data')
  if (index(optv,"d") /= 0) then
    debug = .true.
    verbose = .true.
    write (uout,'(1x,a)') string('# *** Debug mode enabled ***')
  end if
  if (index(optv,"h") /= 0) then
    write (uout,'(1x,a)') string('# +++ Help not yet available')
    goto 1
  end if

  call cpu_time (tiempo1)
  !$ tiempo1 = omp_get_wtime()
  ! Read the input file
  call flush_unit (uout)
  write (uout,'(1x,a)') string('# +++ Begin to read input')
  do while (getline(uin,line))
    lp = 1
    word = lgetword(line,lp)
    subline = line(lp:)

    if (equal(word,'#')) then
      continue

    ! surf file 
    else if (equal(word,'plot')) then
      call parse_plot ()

    ! Param options
    else if (equal(word,'verbose')) then
      verbose = .true.
      write (uout,'(1x,a)') string('# *** Verbose mode enabled ***')

    ! End of input
    else if (equal(word,'end')) then
      exit

    ! Exit main driver
    else if (equal(word,'exit') .or. equal(word,'quit')) then
      goto 1

    else
      call ferror ('plot', 'unknown option', faterr)

    end if
  end do
  write (uout,'(1x,a)') string('# +++ End of read input')
  write (uout,'(1x,a)') string('#')
  call flush_unit (uout)

1 call cpu_time (tiempo2)
  !$ tiempo2 = omp_get_wtime()
  call getdate (wdate)
  write (uout,'(1x,a,f16.6)') string('# Total elapsed time = '), tiempo2-tiempo1
  write (uout,'(" # Check : (",A," WARNINGS, ",A," COMMENTS)")') string(nwarns), string(ncomms)
  write (uout,'(1x,a)') string('# Calculation ends at '//wdate)

end program
