! Copyright (c) 2018 
! Jose Luis Casals Sainz <jluiscasalssainz@gmail.com>,
! Ángel Martín Pendás <angel@fluor.quimica.uniovi.es> and 
! Evelio Francisco Miguelez <evelio@uniovi.es>. 
!
! integ is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or (at
! your option) any later version.
! 
! integ is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
program integ

  !$ use omp_lib
  use mod_prec, only: rp, ip
  use mod_param, only: verbose, debug, init_param, isdata
  use mod_wfn, only: init_wfn, end_wfn, rdwfn
  use mod_gto, only: integwfn
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

  logical :: ok
  integer(kind=ip) :: nthreads
  character(len=mline) :: vchar

  ! Begin program
  call ioinit ()
  call stdargs (optv, fileroot)
  call getdate (wdate)
  write (uout,'(1x,a)') string('#  =====================================================')
  write (uout,'(1x,a)') string('# |      integ: compute energy given the 1/2-RDM        |')
  write (uout,'(1x,a)') string('#  =====================================================')
  write (uout,'(1x,a)') string('#')
  write (uout,'(1x,a)') string('# Calculation starts at '//wdate)
  write (uout,'(1x,a)') string('# *** Begin to initialize data')
  call init_param ()
  call init_wfn ()
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

    ! Parallel options
    else if (equal(word,'threads')) then
      ok = isinteger(nthreads, line, lp)
      ok = ok .and. nthreads.ne.0
      if (.not.ok) call ferror ('integ', 'wrong threads line', faterr) 
      !$ nthreads = abs(nthreads)
      !$ write (uout,'(1x,a,1x,i0)') string('# *** Number of threads changed to :'), nthreads
      !$ call omp_set_num_threads(nthreads)

    ! WFN 
    else if (equal(word,'load')) then
      if (isdata) then
        call end_wfn ()
        isdata = .false.
      end if
      ok = isword(vchar, line, lp)
      if (.not.ok) call ferror ('integ', 'wrong load line', faterr) 
      write (uout,'(1x,a,1x,a)') string('# *** Reading data from'), trim(vchar)
      call rdwfn (vchar)
      isdata = .true.

    ! Orders
    else if (equal(word,'integ')) then
      if (.not.isdata) call ferror ('integ', 'load wfn before integ order', faterr) 
      write (uout,'(1x,a)') string('# *** Computing integrals and energy')
      call integwfn ()

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
      call ferror ('integ', 'unknown option', faterr)

    end if
  end do
  write (uout,'(1x,a)') string('# +++ End of read input')
  write (uout,'(1x,a)') string('#')
  call flush_unit (uout)

  if (isdata) call end_wfn ()

1 call cpu_time (tiempo2)
  !$ tiempo2 = omp_get_wtime()
  call getdate (wdate)
  write (uout,'(1x,a,f16.6)') string('# Total elapsed time = '), tiempo2-tiempo1
  write (uout,'(" # Check : (",A," WARNINGS, ",A," COMMENTS)")') string(nwarns), string(ncomms)
  write (uout,'(1x,a)') string('# Calculation ends at '//wdate)

end program
