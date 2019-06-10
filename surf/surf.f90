! Copyright (c) 2018 
! Jose Luis Casals Sainz <jluiscasalssainz@gmail.com>,
! Ángel Martín Pendás <angel@fluor.quimica.uniovi.es> and 
! Evelio Francisco Miguelez <evelio@uniovi.es>. 
!
! surf is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or (at
! your option) any later version.
! 
! surf is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
program surf

  !$ use omp_lib
  use mod_prec, only: rp, ip
  use mod_param, only: verbose, debug, init_param, isdata
  use mod_surf, only: init_surf, parse_surface
  use mod_wfn, only: init_wfn, end_wfn, rdwfn, epsocc, rmaxatom, ncent, charge
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
  real(kind=rp) :: rmax
  integer(kind=ip) :: nthreads, intz, ic
  character(len=mline) :: vchar

  ! Begin program
  call ioinit ()
  call stdargs (optv, fileroot)
  call getdate (wdate)
  write (uout,'(1x,a)') string('#  =====================================================')
  write (uout,'(1x,a)') string('# |            surf: Compute QTAIM surfaces             |')
  write (uout,'(1x,a)') string('#  =====================================================')
  write (uout,'(1x,a)') string('#')
  write (uout,'(1x,a)') string('# Calculation starts at '//wdate)
  write (uout,'(1x,a)') string('# *** Begin to initialize data')
  call init_param ()
  call init_wfn ()
  call init_surf ()
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
      if (.not.ok) call ferror ('surf', 'wrong threads line', faterr) 
      !$ nthreads = abs(nthreads)
      !$ write (uout,'(1x,a,1x,i0)') string('# *** Number of threads changed to :'), nthreads
      !$ call omp_set_num_threads(nthreads)

    ! Surf orders
    else if (equal(word,'surface')) then
      if (.not.isdata) call ferror ('surf', 'load wfn before surface order', faterr) 
      call parse_surface (vchar)

    ! WFN 
    else if (equal(word,'load')) then
      if (isdata) then
        call end_wfn ()
        isdata = .false.
      end if
      ok = isword(vchar, line, lp)
      if (.not.ok) call ferror ('surf', 'wrong load line', faterr) 
      write (uout,'(1x,a,1x,a)') string('# *** Reading data from'), trim(vchar)
      call rdwfn (vchar)
      isdata = .true.

    else if (equal(word,'epsocc')) then
      ok = isreal(epsocc, line, lp)
      if (.not.ok) call ferror ('surf', 'wrong epsocc line', faterr) 
      epsocc = abs(epsocc)
      write (uout,'(1x,a,1x,e13.6)') string('# *** Variable epsocc changed to :'), epsocc

    else if (equal(word,'rmaxatom')) then
      ok = isreal(rmax, line, lp)
      ok = ok .and. rmax.ne.0.0_rp
      if (.not.ok) call ferror ('surf', 'wrong rmaxatom line', faterr) 
      rmaxatom = abs(rmax)
      write (uout,'(1x,a,1x,e13.6)') string('# *** Variable rmaxatom changed to :'), rmax

    else if (equal(word,'rmaxatom_atom')) then
      ok = isinteger(intz, line, lp)
      ok = ok .and. isreal(rmax, line, lp)
      ok = ok .and. rmax.ne.0.0_rp .and. intz.ne.0_ip .and. intz.le.ncent
      if (.not.ok) call ferror ('surf', 'wrong rmaxatom_atom line', faterr) 
      intz = abs(intz)
      rmax = abs(rmax)
      rmaxatom(intz) = rmax

    else if (equal(word,'rmaxatom_atom_z')) then
      ok = isinteger(intz, line, lp)
      ok = ok .and. isreal(rmax, line, lp)
      ok = ok .and. rmax.ne.0.0_rp .and. intz.ne.0_ip
      if (.not.ok) call ferror ('surf', 'wrong rmaxatom_atom_z line', faterr) 
      intz = abs(intz)
      rmax = abs(rmax)
      do ic = 1,ncent
        if (int(charge(ic),ip).eq.intz) then
          rmaxatom(ic) = rmax
        endif
      end do

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
      call ferror ('surf', 'unknown option', faterr)

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
