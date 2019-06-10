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
module mod_surf

  use mod_prec, only: rp, ip
  implicit none
  private

  integer(kind=ip), parameter, public :: minter = 10 
  integer(kind=ip), parameter, public :: maxtrial = 13
  ! The MAXSTART parameter has to do with the RSEARCH order, 
  ! RSEARCH nsearch nrstart, (rstart(i),i=1,nrstart)
  ! This order can be used when one wants to avoid using the
  ! default starting values of the radial coordinate in
  ! the initial search of interatomic surfaces. If nsearch is 0
  ! the starting values (rstart(i),i=1,nrstart) are used for all
  ! the atoms. If nsearch is different from 0, the rstart(i)
  ! (i=1,nrstart) values are used only for the atom 'nsearch'
  ! and all its symmetry equivalent atoms.
  integer(kind=ip), parameter, public :: maxstart = 40

  integer(kind=ip), public :: inuc
  real(kind=rp), public :: xnuc(3)
  real(kind=rp), allocatable, dimension(:,:), public :: xyzrho
  real(kind=rp), allocatable, dimension(:), public :: rmaxsurf
  integer(kind=ip), allocatable, dimension(:,:), public :: nangleb
  integer(kind=ip), public :: neqsurf
  integer(kind=ip), allocatable, dimension(:), public :: insurf
  real(kind=rp), allocatable, dimension(:,:), public :: rstart
  integer(kind=ip), allocatable, dimension(:), public :: nrsearch
  logical, allocatable, dimension(:), public :: lstart
  real(kind=rp), allocatable, dimension(:,:), public :: rlimsurf
  integer(kind=ip), allocatable, dimension(:), public :: nlimsurf

  ! options
  logical, public :: rotgrid
  real(kind=rp), public :: angx, angy, angz
  ! Precision in interatomic surfaces. 
  real(kind=rp), public :: epsilon
  ! If a non nuclear maximum is found whose distance to any nucleus 
  ! is larger than EPSISCP, promolden STOPS. This control of non nuclear
  ! maxima is performed in the 'iscp.f' routine. When Hydrogen atoms
  ! are involved it is very convenient to change the default value
  ! for EPSISCP (0.08) to a smaller value by using the 'EPSISCP value'
  ! order in the input file.
  real(kind=rp), public :: epsiscp
  ! NTRIAL  = number of sub-intervals in the search of the surface.
  ! RPRIMER = First point in the in the search of the surface.
  ! GEOFAC  = ((rmaxsurf-0.1d0)/rprimer)**(1d0/(ntrial-1))
  !           (RMAXSURF is defined in passed in integ.inc)
  integer(kind=ip), public :: ntrial
  real(kind=rp), public :: rprimer
  real(kind=rp), public :: geofac

  public :: init_surf, parse_surface

contains

  subroutine parse_surface (filename)

    use mod_param, only: pi
    use mod_wfn, only: ncent, charge
    use mod_io, only: ferror, faterr, equal, isinteger, isreal, warning, &
                      uin, string, uout, lgetword, getline, mline
    implicit none
    character(len=*), intent(in) :: filename

    integer(kind=ip) :: lp
    character(len=:), allocatable :: line, subline, word

    logical :: ok
    real(kind=rp) :: rmax
    integer(kind=ip) :: intz, ic, ii, iqudt, nphi, ntheta, npang
    character(len=mline), dimension(5) :: rqudstr

    call allocate_space_for_integ (ncent)

    ! Read the input file
    write (uout,'(1x,a)') string('# *** Reading surface options ***')
    do while (getline(uin,line))
      lp = 1
      word = lgetword(line,lp)
      subline = line(lp:)

      if (equal(word,'#')) then
        continue

      else if (equal(word,'agrid')) then
        ok = isinteger(iqudt, line, lp)
        ok = ok .and. isinteger(ntheta, line, lp)
        ok = ok .and. isinteger(nphi, line, lp)
        ok = ok .and. ntheta.ne.0_ip .and. nphi.ne.0_ip .and. iqudt.ne.0
        if (.not.ok) call ferror ('mod_surf', 'wrong agrid line', faterr) 
        iqudt = abs(iqudt)
        ntheta = abs(ntheta)
        nphi = abs(nphi)
        nangleb(:,1) = 0
        nangleb(:,2) = ntheta
        nangleb(:,3) = nphi 
        nangleb(:,4) = iqudt
        write (uout,'(1x,a,3(1x,i0))') string('# *** Surface agrid (iqudt,ntheta,nphi) :'), iqudt, ntheta, nphi

      else if (equal(word,'agrid_atom')) then
        ok = isinteger(intz, line, lp)
        ok = ok .and. isinteger(iqudt, line, lp)
        ok = ok .and. isinteger(ntheta, line, lp)
        ok = ok .and. isinteger(nphi, line, lp)
        ok = ok .and. ntheta.ne.0_ip .and. nphi.ne.0_ip .and. intz.ne.0_ip .and. intz.le.ncent .and. iqudt.ne.0
        if (.not.ok) call ferror ('mod_surf', 'wrong agrid_atom line', faterr) 
        iqudt = abs(iqudt)
        intz = abs(intz)
        ntheta = abs(ntheta)
        nphi = abs(nphi)
        nangleb(intz,1) = 0
        nangleb(intz,2) = ntheta
        nangleb(intz,3) = nphi
        nangleb(intz,4) = iqudt

      else if (equal(word,'agrid_atom_z')) then
        ok = isinteger(intz, line, lp)
        ok = ok .and. isinteger(iqudt, line, lp)
        ok = ok .and. isinteger(ntheta, line, lp)
        ok = ok .and. isinteger(nphi, line, lp)
        ok = ok .and. ntheta.ne.0_ip .and. nphi.ne.0_ip .and. intz.ne.0_ip .and. iqudt.ne.0
        if (.not.ok) call ferror ('mod_surf', 'wrong agrid_atom_z line', faterr) 
        iqudt = abs(iqudt)
        intz = abs(intz)
        ntheta = abs(ntheta)
        nphi = abs(nphi)
        do ic = 1,ncent
          if (int(charge(ic),ip).eq.intz) then
            nangleb(ic,1) = 0
            nangleb(ic,2) = ntheta
            nangleb(ic,3) = nphi
            nangleb(ic,4) = iqudt
          endif
        end do

      else if (equal(word,'lebedev')) then
        ok = isinteger(npang, line, lp)
        ok = ok .and. npang.ne.0_ip
        if (.not.ok) call ferror ('mod_surf', 'wrong lebedev line', faterr) 
        npang = abs(npang)
        call good_lebedev (npang)
        nangleb(:,1) = 1
        nangleb(:,2) = npang
        nangleb(:,3) = 0
        nangleb(:,4) = 0
        write (uout,'(1x,a,1x,i0)') string('# *** Surface lebedev changed to :'), npang

      else if (equal(word,'lebedev_atom')) then
        ok = isinteger(intz, line, lp)
        ok = ok .and. isinteger(npang, line, lp)
        ok = ok .and. npang.ne.0_ip .and. intz.ne.0_ip .and. intz.le.ncent
        if (.not.ok) call ferror ('mod_surf', 'wrong lebedev_atom line', faterr) 
        intz = abs(intz)
        npang = abs(npang)
        call good_lebedev (npang)
        nangleb(intz,1) = 0
        nangleb(intz,2) = npang
        nangleb(intz,3) = 0
        nangleb(intz,4) = 0

      else if (equal(word,'lebedev_atom_z')) then
        ok = isinteger(intz, line, lp)
        ok = ok .and. isinteger(npang, line, lp)
        ok = ok .and. npang.ne.0_ip .and. intz.ne.0_ip
        if (.not.ok) call ferror ('mod_surf', 'wrong lebedev_atom_z line', faterr) 
        intz = abs(intz)
        npang = abs(npang)
        call good_lebedev (npang)
        do ic = 1,ncent
          if (int(charge(ic),ip).eq.intz) then
            nangleb(ic,1) = 0
            nangleb(ic,2) = npang
            nangleb(ic,3) = 0
            nangleb(ic,4) = 0
          endif
        end do

      else if (equal(word,'epsiscp')) then
        ok = isreal(epsiscp, line, lp)
        ok = ok .and. epsiscp.ne.0.0_rp
        if (.not.ok) call ferror ('mod_surf', 'wrong epsiscp line', faterr) 
        epsiscp = abs(epsiscp)
        write (uout,'(1x,a,1x,e13.6)') string('# *** Variable epsiscp changed to :'), epsiscp

      else if (equal(word,'epsilon')) then
        ok = isreal(epsilon, line, lp)
        ok = ok .and. epsilon.ne.0.0_rp
        if (.not.ok) call ferror ('mod_surf', 'wrong epsilon line', faterr) 
        epsilon = abs(epsilon)
        write (uout,'(1x,a,1x,e13.6)') string('# *** Variable epsilon changed to :'), epsilon

      else if (equal(word,'ntrial')) then
        ok = isinteger(ntrial, line, lp)
        ok = ok .and. ntrial.ne.0_ip
        if (.not.ok) call ferror ('mod_surf', 'wrong ntrial line', faterr) 
        ntrial = abs(ntrial)
        if (mod(ntrial,2).eq.0) ntrial = ntrial + 1_ip
        write (uout,'(1x,a,1x,i0)') string('# *** Variable ntrial changed to :'), ntrial

      else if (equal(word,'rprimer')) then
        ok = isreal(rprimer, line, lp)
        if (.not.ok) call ferror ('mod_surf', 'wrong rprimer line', faterr) 
        rprimer = abs(rprimer)
        write (uout,'(1x,a,1x,e13.6)') string('# *** Variable rprimer changed to :'), rprimer

      else if (equal(word,'rmaxsurf')) then
        ok = isreal(rmax, line, lp)
        ok = ok .and. rmax.ne.0.0_rp
        if (.not.ok) call ferror ('mod_surf', 'wrong rmaxsurf line', faterr) 
        rmaxsurf = abs(rmax)
        write (uout,'(1x,a,1x,e13.6)') string('# *** Variable rmaxsurf changed to :'), rmax

      else if (equal(word,'rmaxsurf_atom')) then
        ok = isinteger(intz, line, lp)
        ok = ok .and. isreal(rmax, line, lp)
        ok = ok .and. rmax.ne.0.0_rp .and. intz.ne.0_ip .and. intz.le.ncent
        if (.not.ok) call ferror ('mod_surf', 'wrong rmaxsurf_atom line', faterr) 
        intz = abs(intz)
        rmax = abs(rmax)
        rmaxsurf(intz) = rmax

      else if (equal(word,'rmaxsurf_atom_z')) then
        ok = isinteger(intz, line, lp)
        ok = ok .and. isreal(rmax, line, lp)
        ok = ok .and. rmax.ne.0.0_rp .and. intz.ne.0_ip
        if (.not.ok) call ferror ('mod_surf', 'wrong rmaxsurf_atom_z line', faterr) 
        intz = abs(intz)
        rmax = abs(rmax)
        do ic = 1,ncent
          if (int(charge(ic),ip).eq.intz) then
            rmaxsurf(ic) = rmax
          endif
        end do

      else if (equal(word,'rotate')) then
        ok = isreal(angx,line,lp)
        ok = ok .and. isreal(angy,line,lp)
        ok = ok .and. isreal(angz,line,lp)
        if (.not.ok) call ferror('dosurf', 'wrong rotate line', faterr) 
        rotgrid = .true.
        if (abs(angx).lt.1d-2 .and. abs(angy).lt.1d-2 .and. abs(angz).lt.1d-2) rotgrid = .false.
        if (rotgrid) then
          angx = angx*pi/180.0_rp
          angy = angy*pi/180.0_rp
          angz = angz*pi/180.0_rp
          write (uout,'(1x,a,1x,3f6.2)') string('# *** Grid will be rotated (x,y,z) :'), angx,angy,angz
        else
          call ferror ('mod_surf', 'grid will be not rotated chech angs', warning)
        end if
      
      ! End of input
      else if (equal(word,'endsurface')) then
        exit

      else
        call ferror ('mod_surface', 'unknown option', faterr)

      end if
    end do

    ! Info
    rqudstr(1) = 'Gauss-Legendre'
    rqudstr(2) = 'Clenshaw-Curtis'
    rqudstr(3) = 'Gauss-Chebychev 1st kind'
    rqudstr(4) = 'Gauss-Chebychev 2nd kind'
    rqudstr(5) = 'Perez-Jorda (Gauss-Chebychev) 2nd kind'
    write (uout,'(1x,a)') string('# *** End reading surface options ***')
    write (uout,'(1x,a)') string('# *** Computing surfaces')
    write (uout,'(1x,a)') string('# Will be computed for the following atoms')
    write (uout,'(1x,a,1x,i0)') string('# Number of surfaces :'), neqsurf
    do ic = 1,neqsurf
      write (uout,'(1x,a,1x,i0,1x,i0)') string('#'), ic, insurf(ic)
    end do
    write (uout,'(1x,a)') string('# Rmaxsurf for the atoms') 
    do ic = 1,ncent
      write (uout,'(1x,a,1x,i0,1x,f6.2)') string('#'), ic, rmaxsurf(ic)
    end do
    write (uout,'(1x,a,1x,e13.6)') string('# Surface precision :'), epsilon
    write (uout,'(1x,a,1x,e13.6)') string('# EPSISCP parameter :'), epsiscp
    write (uout,'(1x,a,1x,i0)') string('# Ntrial :'), ntrial
    write (uout,'(1x,a,1x,e13.6)') string('# Rprimer :'), rprimer
    write (uout,'(1x,a)') string('# Follow integration parameters')
    do ic = 1,neqsurf
      ii = insurf(ic)
      if (nangleb(ii,1).eq.1) then
        write (uout, '(1x,a,1x,i0,1x,a,1x,i0)') &
        string('# Atom'), ii, string('lebedev points'), nangleb(ii,2)
      else if (nangleb(ii,1).eq.0) then
        write (uout, '(1x,a)') string('# Phi quadrature is always trapezoidal')
        write (uout, '(1x,a,1x,i0,1x,a,1x,i0,1x,i0,1x,a)') &
        string('# Atom'), ii, string('(ntheta,nphi,iqudt'), nangleb(ii,2), nangleb(ii,3), &
        string(rqudstr(nangleb(ii,4)))
      end if
    end do

    ! Do the job
    call surface (filename)
    call deallocate_space_for_integ ()

  end subroutine parse_surface

! Determine the surface for all atoms or specified atoms
  subroutine surface (filename)

    !$ use omp_lib, only: omp_get_wtime
    use mod_prec, only: i4=>ip
    use mod_memory, only: alloc, free
    use mod_param, only: pi
    use mod_io, only: fourchar, uout, mline, string, flush_unit
    use mod_quad, only: weightheta
    implicit none
 
    character(len=*), intent(in) :: filename
 
    integer(kind=i4) :: ntheta, nphi, npang, iqudt
    real(kind=rp), allocatable, dimension(:) :: ct, st, cp, sp, angw
    real(kind=rp), allocatable, dimension(:) :: tp, tw
    integer(kind=i4) :: i, ii, ip, it, j, k, nsurf
    integer(kind=i4) :: lsu, ltxt
    real(kind=rp) :: thang, time1, time2, phi
    real(kind=rp) :: rsurf(minter,2), delphi, rmaxs, rmins
    character(len=mline) :: files
    character(len=4) :: d4

    ! Init
    ltxt = 999
    lsu = 998
    files = trim(filename)//".surf"
  
    ! Find nucleus
    call findnuc ()

    ! Begin
    call flush_unit (uout)
    do ii = 1,neqsurf
      inuc = insurf(ii)
      d4 = fourchar(inuc)
      xnuc(:) = xyzrho(inuc,:)
      call flush_unit (uout)
      write (uout,'(1x,a,1x,i0)') string('# Computing SURFACE for atom'), inuc
      if (nangleb(inuc,1).eq.1_i4) then
        npang = nangleb(inuc,2)
        call good_lebedev (npang)
        call alloc ('surface', 'ct', ct, npang)
        call alloc ('surface', 'st', st, npang)
        call alloc ('surface', 'cp', cp, npang)
        call alloc ('surface', 'sp', sp, npang)
        call alloc ('surface', 'angw', angw, npang)
        call lebgrid (ct,st,cp,sp,angw,npang)
      else
        ntheta = nangleb(inuc,2)
        nphi = nangleb(inuc,3)
        iqudt = nangleb(inuc,4)
        npang = ntheta*nphi  
        call alloc ('surface', 'tp', tp, ntheta)
        call alloc ('surface', 'tw', tw, ntheta)
        call alloc ('surface', 'ct', ct, npang)
        call alloc ('surface', 'st', st, npang)
        call alloc ('surface', 'cp', cp, npang)
        call alloc ('surface', 'sp', sp, npang)
        call alloc ('surface', 'angw', angw, npang)
        call weightheta (iqudt,tp,tw,ntheta)
        delphi = 2.0_rp*pi/nphi
        i = 0_i4
        do ip = 0,nphi-1
          phi = ip*delphi
          do it = 1,ntheta
            i = i + 1_i4
            thang = tp(it)
            ct(i) = thang
            st(i) = sqrt(1.0_rp-thang*thang)
            cp(i) = cos(phi)
            sp(i) = sin(phi)
            angw(i) = tw(it)*delphi
          end do
        end do
        call free ('surface', 'tp', tp)
        call free ('surface', 'tw', tw)
      end if
      call allocate_space_for_surface (npang,minter)
      call cpu_time (time1)
      !$ time1 = omp_get_wtime()
      !$omp parallel default(none) &
      !$omp private(j,nsurf,rsurf) &
      !$omp shared(npang,ct,st,cp,sp,epsilon,rlimsurf,nlimsurf)
      !$omp do schedule(dynamic)
      do j = 1,npang
        call surf (ct(j),st(j),cp(j),sp(j),rsurf,nsurf)
        do k = 1,nsurf
          rlimsurf(j,k) = rsurf(k,2)
        end do
        nlimsurf(j) = nsurf
      end do
      !$omp end do nowait
      !$omp end parallel
      call cpu_time (time2)
      !$ time2 = omp_get_wtime()
      write (uout,'(1x,a,1x,f12.5)') string('# Elapsed seconds :'), time2-time1
      open (ltxt,file=trim(files)//"-txt"//d4)
      open (lsu,file=trim(files)//d4,form='unformatted')
      write (ltxt,1111) npang, inuc
      write (lsu) npang, inuc
      write (ltxt,3333) (nlimsurf(j),j=1,npang)
      write (ltxt,1090)
      write (lsu) (nlimsurf(j),j=1,npang)
      rmins = 1000_rp
      rmaxs = 0.0_rp
      do j = 1,npang
        nsurf = nlimsurf(j)
        rmins = min(rmins,rlimsurf(j,1))
        rmaxs = max(rmaxs,rlimsurf(j,nsurf))
        write (ltxt,2222) ct(j),st(j),cp(j),sp(j),angw(j),(rlimsurf(j,k),k=1,nsurf)
        write (lsu) ct(j),st(j),cp(j),sp(j),angw(j),(rlimsurf(j,k),k=1,nsurf)
      end do
      write (ltxt,2222) rmins,rmaxs
      write (lsu) rmins,rmaxs
      write (lsu) xnuc(1), xnuc(2), xnuc(3)
      write (ltxt,2223) xnuc
      close (ltxt)
      close (lsu)
      call deallocate_space_for_surface ()
      call free ('surface', 'ct', ct)
      call free ('surface', 'ct', st)
      call free ('surface', 'ct', cp)
      call free ('surface', 'ct', sp)
      call free ('surface', 'angw', angw)
    end do

1090 format (9x,'cos(theta)',13x,'sin(theta)',13x,'cos(phi)',15x,'sin(phi)',15x,'weight')
1111 format (2(1x,i5),' <--- (Angular points & Atom)')
3333 format (20(1x,i2),4x,'(Surface intersections)')
2222 format (15(1x,F22.15))
2223 format (3(1x,F22.15))
 
  end subroutine

  ! Determine the limit of the zero flux surface from xpoint along
  ! the theta phi direction. 
  subroutine surf (cot,sit,cop,sip,rsurf,nsurf)

    use mod_io, only: ferror, faterr
    use mod_wfn, only: ncent
    implicit none
    real(kind=rp), parameter :: half = 0.5_rp
 
    integer(kind=ip), intent(out) :: nsurf
    real(kind=rp), intent(in) :: cot, sit, cop, sip
    real(kind=rp), intent(out) :: rsurf(minter,2) 
 
    logical :: inf, good
    integer(kind=ip) :: nintersec, nt
    integer(kind=ip) :: isurf(minter,2), i, ia, ib, im, ncount
    real(kind=rp) :: sintcosp, sintsinp
    real(kind=rp) :: xin(3), xfin(3), xpoint(3), xmed(3), ra, ract
    real(kind=rp) :: xdeltain(3), xsurf(0:minter,3), rb, rm, cost
 
    cost = cot
    sintcosp = sit*cop
    sintsinp = sit*sip
 
    if (ncent.eq.1) then
      nsurf = 1
      rsurf(1,1) = 0.0_rp
      rsurf(1,2) = rmaxsurf(1)
      return
    end if
 
    inf = .false.
    ncount = 0_ip
    nintersec = 0_ip
    ia = inuc
    ra = 0.0_rp
    geofac = ((rmaxsurf(inuc)-0.1_rp)/rprimer)**(1.0_rp/(ntrial-1))
    if (.not.lstart(ia)) then
      do i = 1,ntrial
        ract = rprimer*geofac**(i-1)
        xdeltain(1) = ract*sintcosp
        xdeltain(2) = ract*sintsinp
        xdeltain(3) = ract*cost    
        xpoint(:) = xnuc(:) + xdeltain(:)
        call odeint (xpoint,0.1_rp,1.0_rp,inf,epsilon,xnuc)
        good = iscp(xpoint,ib)
        rb = ract
        if (ib.ne.ia .and. (ia.eq.inuc .or. ib.eq.inuc)) then
          if (ia.ne.inuc .or. ib.ne.-1) then
            nintersec = nintersec + 1_ip
            if (nintersec.gt.minter) then
              call ferror ('surf', 'increase minter in mod_surf', faterr)
            end if
            xsurf(nintersec,1) = ra
            xsurf(nintersec,2) = rb
            isurf(nintersec,1) = ia
            isurf(nintersec,2) = ib
          end if
        end if
        ia = ib
        ra = rb
      end do        ! looking for intersections
    else
      do i = 1,nrsearch(inuc)
        ract = rstart(inuc,i)
        xdeltain(1) = ract*sintcosp
        xdeltain(2) = ract*sintsinp
        xdeltain(3) = ract*cost
        xpoint(:) = xnuc(:) + xdeltain(:)
        call odeint (xpoint,0.1_rp,1.0_rp,inf,epsilon,xnuc)
        good = iscp(xpoint,ib)
        rb = ract
        if (ib.ne.ia .and. (ia.eq.inuc .or. ib.eq.inuc)) then
          if (ia.ne.inuc .or. ib.ne.-1) then
            nintersec = nintersec + 1_ip
            if (nintersec.gt.minter) then
              call ferror ('surf', 'increase minter in mod_surf', faterr)
            end if
            xsurf(nintersec,1) = ra
            xsurf(nintersec,2) = rb
            isurf(nintersec,1) = ia
            isurf(nintersec,2) = ib
          end if
        end if
        ia = ib
        ra = rb
      end do        ! looking for intersections
    end if

    ! We have now a consistent set of trial points. 
    ! Let us refine by bipartition, consistency check added. 
    ! No other nuclei basins can be found other than those explicit in isurf
    do i = 1,nintersec
      ia = isurf(i,1)
      ib = isurf(i,2)
      ra = xsurf(i,1)
      rb = xsurf(i,2)
      xin(1) = xnuc(1) + ra*sintcosp
      xin(2) = xnuc(2) + ra*sintsinp
      xin(3) = xnuc(3) + ra*cost
      xfin(1) = xnuc(1) + rb*sintcosp
      xfin(2) = xnuc(2) + rb*sintsinp
      xfin(3) = xnuc(3) + rb*cost
      do while (abs(ra-rb).gt.epsilon)
        ! Mean Value 
        xmed(:) = half*(xfin(:)+xin(:))    
        rm = (ra+rb)*half
        xpoint(:) = xmed(:)
        call odeint (xpoint,0.1_rp,1.0_rp,inf,epsilon,xnuc)
        good = iscp(xpoint,im)
        ! bipartition
        if (im.eq.ia) then
          xin(:) = xmed(:)
          ra = rm
        else if (im.eq.ib) then
          xfin(:) = xmed(:)
          rb = rm
        else
          if (ia.eq.inuc) then
            xfin(:) = xmed(:)
            rb = rm
          else
            xin(:) = xmed(:)
            ra = rm
          end if
        end if
      end do
      ! found. Mean value
      xpoint(:) = half*(xfin(:)+xin(:))    
      xsurf(i,3) = (ra+rb)*half
    end do

    ! organize pairs
    nsurf = 0_ip
    xsurf(0,3) = 0.0_rp
    ia = inuc
    nsurf = nintersec
    do i = 1,nsurf
      rsurf(i,2) = xsurf(i,3)
    end do
    nt = mod(nintersec,2)
    if (nt.eq.0) then
      nsurf = nsurf + 1_ip
      rsurf(nsurf,2) = rmaxsurf(inuc)
    end if
 
  end subroutine

  subroutine odeint (ystart,h1,iup,inf,eps,xnuc)
 
    use mod_io, only: faterr, ferror, string, warning
    use mod_fields, only: pointshells
    implicit none
    integer(kind=ip), parameter :: maxstp = 350
    real(kind=rp), parameter :: tiny = 1d-40
    real(kind=rp), parameter :: epsg = 1d-10
    real(kind=rp), parameter :: epsg1 = 1d-10

    ! Arguments
    logical, intent(inout) :: inf
    real(kind=rp), intent(inout) :: ystart(3)
    real(kind=rp), intent(in) :: h1
    real(kind=rp), intent(in) :: iup
    real(kind=rp), intent(in) :: eps
    real(kind=rp), intent(in) :: xnuc(3)

    ! Local vars
    integer(kind=ip) :: nstp, nuc
    real(kind=rp), parameter :: hmin = 0.0_rp ! minimal step size
    real(kind=rp) :: x1 ! intial point
    real(kind=rp) :: x2 ! final point
    real(kind=rp) :: p(3) ! initial solution point
    real(kind=rp) :: h ! initial step size
    real(kind=rp) :: x ! update point
    real(kind=rp) :: hnext ! next steep size
    real(kind=rp) :: dydx(3), y(3), yscal(3)
    real(kind=rp) :: a1, a2, a3
    real(kind=rp) :: rho, grad(3), gradmod
 
    inf = .false.
    p(1) = ystart(1)
    p(2) = ystart(2)
    p(3) = ystart(3)
    call pointshells (p,rho,grad,gradmod,inuc)
    if (gradmod.lt.epsg .and. rho.lt.epsg1) then
      inf = .true.
      return
    end if
 
    x1 = 0.0_rp
    x2 = 1d40*iup
    x = x1 ! initial point
    !h = sign(h1,x2-x1) ! initial steep size
    h = min(h1,x2-x1) ! initial steep size
    y(:) = ystart(:) ! initial point 
 
    do nstp = 1,maxstp
      call pointshells (y,rho,grad,gradmod,inuc)
      !write (*,*) y,rho,grad,gradmod,inuc
      dydx(:) = grad(:)
      yscal(:) = max(abs(y(:))+abs(h*dydx(:))+tiny,eps)
      if ((x+h-x2)*(x+h-x1).gt.0.0_rp) h = x2 - x
      call rkqs (y,dydx,x,h,eps,yscal,hnext)
      if ((x-x2)*(x2-x1).ge.0.0_rp .or. iscp(y,nuc)) then
        ystart(:) = y(:)
        return
      end if
      if (abs(hnext).lt.hmin) then
        call ferror ('mod_odeint/odeint', 'stepsize small than minimum', faterr)
      end if
      if (nstp.eq.maxstp) then
        call ferror ('mod_odeint/odeint', 'reached maxstp', warning)
      end if 
      h = hnext
    end do

    ! Test if the point is far from RMAXSURF from current atom. 
    a1 = y(1) - xnuc(1)
    a2 = y(2) - xnuc(2)
    a3 = y(3) - xnuc(3)
    if ((a1*a1+a2*a2+a3*a3).ge.5.0*5.0) then
      inf = .true.
      return
    else
      call ferror ('mod_odeint/odeint', 'Non nuclear maximum at : ' &
                                         //string(y(1),'e')//' '    &  
                                         //string(y(2),'e')//' '    &  
                                         //string(y(3),'e'), faterr) 
    end if
 
  end subroutine 

  subroutine findnuc ()

    use mod_io, only: string, faterr, uout, ferror
    use mod_wfn, only: xyz, charge, ncent
    use mod_fields, only: pointshells
    implicit none

    integer(kind=ip) :: i
    real(kind=rp) :: rho, grad(3), gradmod, p(3)
    logical :: inf

    do i = 1,ncent
      inuc = i
      p(:) = xyz(i,:)
      call gradrho (p,0.05_rp,1,inf)
      call pointshells (p,rho,grad,gradmod,inuc)
      if (gradmod.gt.1d-4) then
        if (charge(i).gt.2.0_rp) then
          write (uout,'(1x,a,i0,a)') '# Assuming nuclei ', i, ' position: check!'
          xyzrho(i,:) = xyz(i,:)
        else
          call ferror('findnuc', 'failed finding nucleus '//string(i), faterr)
        end if
      else
        write (uout,'(1x,a,i0,a)') '# Assuming nuclei ', i, ' position: check!'
        xyzrho(i,:) = xyz(i,:)
      end if
    end do
 
  end subroutine

  ! Integration of a trajectory in the vector field of the
  ! electron density.
  ! Input data:
  ! xpoint() .... starting point of the trajectory
  ! step ........ integration step. Enter negative value if default
  !            value is wanted.
  ! iup ......... +1 tells the routine to climb up in the gradient field
  !               (usually going to end in a nucleus).
  !               -1 tells the routine to go down in the field.
  ! Output data:
  ! xpoint() .... end point of the trajectory 
  subroutine gradrho (xpoint,hi,iup,inf)

    use mod_fields, only: pointshells
    implicit none
    real(kind=rp), parameter :: epsg = 1d-10
    real(kind=rp), parameter :: eps = 1d-8
    real(kind=rp), parameter :: epsg1 = 1d-10
    real(kind=rp), parameter :: hminimal = 1d-40
    integer(kind=ip), parameter :: mstep = 500
 
    logical, intent(inout) :: inf
    integer(kind=ip), intent(in) :: iup
    real(kind=rp), intent(in) :: hi
    real(kind=rp), intent(inout) :: xpoint(3)
 
    real(kind=rp) :: gradmod, rho, grad(3)
    integer(kind=ip) :: ier, niter, npoints, i
    real(kind=rp) :: xtemp(3), grdt(3), h0, escalar, grdmodule
 
    h0 = hi
    npoints = 1_ip
    inf = .false.
    call pointshells (xpoint,rho,grad,gradmod,inuc)
    if (gradmod.lt.epsg .and. rho.lt.epsg1) then
      inf = .true.
      return
    end if
    grdt(1) = grad(1)/gradmod
    grdt(2) = grad(2)/gradmod
    grdt(3) = grad(3)/gradmod
    grdmodule = gradmod
 
    escalar = 1_ip
    niter = 1_ip
    do while (grdmodule.gt.eps .and. niter.lt.mstep)
      niter = niter + 1_ip
      ier = 1_ip
      do while (ier.ne.0)
        xtemp(1) = xpoint(1) + h0*iup*grdt(1)
        xtemp(2) = xpoint(2) + h0*iup*grdt(2)
        xtemp(3) = xpoint(3) + h0*iup*grdt(3)
        call pointshells (xtemp,rho,grad,gradmod,inuc)
        escalar = 0.0_rp
        do i = 1,3
          escalar = escalar + grdt(i)*(grad(i)/(gradmod+1d-80))
        end do
        ! Good direction
        if (escalar.lt.0.707_rp) then
          ! It should'nt happen that h0 goes to zero, except if there
          ! are problems with the gradient, for instance the gradient
          ! has discontinuities or large rounding errors. Anyway, we
          ! introduce a safety check to avoid nasty infinite loops if
          ! h0 = 0.
          if (h0.ge.hminimal) then
            h0 = h0/2.0_rp
            ier = 1_ip
          else
            ier = 0_ip
          end if
        else
          if (escalar.gt.0.9_rp) h0 = min(hi,h0*1.6_rp)
          ier = 0_ip
          ! Yes
          do i = 1,3
            xpoint(i) = xtemp(i)
            grdt(i) = grad(i)/gradmod
          enddo
          grdmodule = gradmod
        end if
      end do
    end do
 
  end subroutine

  logical function iscp (p,nuc)
 
    use mod_fields, only: pointshells
    use mod_wfn, only: ncent
    implicit none
    integer(kind=ip), intent(out) :: nuc
    real(kind=rp), intent(in) :: p(3)
 
    real(kind=rp) :: gradmod, grad(3), rho
    real(kind=rp) :: x(3)
    integer(kind=ip) :: i

    iscp = .false.
    nuc = 0
    x(:) = p(:)
    call pointshells (x,rho,grad,gradmod,inuc)
    do i = 1,ncent
      if (abs(p(1)-xyzrho(i,1)).lt.epsiscp .and.  &
          abs(p(2)-xyzrho(i,2)).lt.epsiscp .and.  &
          abs(p(3)-xyzrho(i,3)).lt.epsiscp) then
        iscp = .true.
        nuc = i
      end if
    end do
 
    if (gradmod.le.1d-10) then
      iscp = .true.
      if (rho.le.1d-10) nuc = -1
    end if
 
  end function

  subroutine rkqs (y,dydx,x,htry,eps,yscal,hnext)
      
    use mod_io, only: ferror, faterr
    implicit none
    real(kind=rp), parameter :: safety = 0.9_rp
    real(kind=rp), parameter :: pgrow = -0.2_rp
    real(kind=rp), parameter :: pshrnk = -0.25_rp
    real(kind=rp), parameter :: errcon = 1.89d-4
    
    ! Arguments
    real(kind=rp), dimension(3), intent(in) :: dydx
    real(kind=rp), dimension(3), intent(inout) :: y
    real(kind=rp), dimension(3), intent(in) :: yscal
    real(kind=rp), intent(inout) :: x
    real(kind=rp), intent(in) :: eps
    real(kind=rp), intent(in) :: htry
    real(kind=rp), intent(out) :: hnext

    ! Local vars
    integer(kind=ip) :: i
    real(kind=rp), dimension(3) :: yerr, ytemp
    real(kind=rp) :: h, errmax, htemp, xnew
      
    h = htry
    hnext = 0.0_rp
    errmax = 0.0_rp
    yerr = 0.0_rp

    do
      call rkck (y,dydx,h,ytemp,yerr)
      ! adaptive and error estimation 
      errmax = 0.0_rp
      do i = 1,3
        errmax = max(errmax,abs(yerr(i)/yscal(i)))
      end do
      errmax = errmax/eps
      if (errmax.gt.1.0_rp) then
        htemp = safety*h*(errmax**pshrnk)
        !h = sign(max(abs(htemp),0.1_rp*abs(h)),h)
        h = min(max(abs(htemp),0.1_rp*abs(h)),h)
        xnew = x + h
        if (xnew.eq.x) then
          call ferror ('mod_odeint/rkqs', 'stepsize underflow', faterr)
          return
        end if
        cycle
      else
        if (errmax.gt.errcon) then
          hnext = safety*h*(errmax**pgrow)
        else
          hnext = 5.0_rp*h
        end if
        x = x + h
        y = ytemp
        return
      end if
    end do
 
  end subroutine rkqs

  ! Runge-Kutta-Cash-Karp embedded 4(5)-order, with local extrapolation.
  ! Runge-Kutta embedded 4th order, Cash-Karp parametrization.
  ! ##### 6 stages, 5th order
  ! This scheme is due to Cash and Karp, see [1].
  !    
  !   0    | 0
  !   1/5	 | 1/5
  !   3/10 | 3/40	         9/40
  !   3/5	 | 3/10	         -9/10	      6/5
  !   1	   | -11/54	       5/2	        -70/27	    35/27
  !   7/8	 | 1631/55296    175/512      575/13824   44275/110592     253/4096     0
  !  ----------------------------------------------------------------------------------------
  !        | 37/378        0           250/621      125/594          0            512/1771
  !        | 2825/27648    0           18575/48384  13525/55296      277/14336    1/4
  ! 
  ! [1] *A variable order Runge-Kutta method for initial value problems with rapidly varying right-hand sides*, J. R. Cash,
  ! A. H. Karp, ACM Transactions on Mathematical Software, vol. 16,  pp. 201--222, 1990, doi:10.1145/79505.79507.
  subroutine rkck (y,dydx,h,yout,yerr)
      
    use mod_fields, only: pointshells
    implicit none
  
    ! Arguments
    real(kind=rp), intent(in) :: dydx(3)
    real(kind=rp), intent(in) :: y(3)
    real(kind=rp), intent(out) :: yerr(3)
    real(kind=rp), intent(out) :: yout(3)
    real(kind=rp), intent(in) :: h
  
    ! Local vars
    real(kind=rp) :: rho, grad(3), gradmod
    real(kind=rp) :: ak2(3), ak3(3), ak4(3), ak5(3), ak6(3)
    real(kind=rp), parameter :: b21=0.2d0,                           &
                                b31=3.0d0/40.0d0,                    &
                                b32=9.0d0/40.0d0,                    &
                                b41=0.3d0,b42=-0.9d0,b43=1.2d0,      &
                                b51=-11.0d0/54.0d0,b52=2.5d0,        &
                                b53=-70.0d0/27.0d0,b54=35.d0/27.0d0, &
                                b61=1631.0d0/55296.0d0,              &
                                b62=175.0d0/512.0d0,                 &
                                b63=575.0d0/13824.0d0,               &
                                b64=44275.0d0/110592.0d0,            &
                                b65=253.0d0/4096.0d0                  
    real(kind=rp), parameter :: c1=37.0d0/378.0d0,c3=250.0d0/621.0d0,&
                                c4=125.0d0/594.0d0,                  &
                                c6=512.0d0/1771.0d0                   
    real(kind=rp), parameter :: dc1=c1-2825.0d0/27648.0d0,           &
                                dc3=c3-18575.0d0/48384.0d0,          &
                                dc4=c4-13525.0d0/55296.0d0,          &
                                dc5=-277.0d0/14336.0d0,              &
                                dc6=c6-0.25d0                         
   
    yout = y + b21*h*dydx
   
    call pointshells (yout,rho,grad,gradmod,inuc)
    ak2(:) = grad(:)
    yout = y + h*(b31*dydx+b32*ak2)

    call pointshells (yout,rho,grad,gradmod,inuc)
    ak3(:) = grad(:)
    yout = y + h*(b41*dydx+b42*ak2+b43*ak3)
  
    call pointshells (yout,rho,grad,gradmod,inuc)
    ak4(:) = grad(:)
    yout = y + h*(b51*dydx+b52*ak2+b53*ak3+b54*ak4)
  
    call pointshells (yout,rho,grad,gradmod,inuc)
    ak5(:) = grad(:)
    yout = y + h*(b61*dydx+b62*ak2+b63*ak3+b64*ak4+b65*ak5)

    call pointshells (yout,rho,grad,gradmod,inuc)
    ak6(:) = grad(:)
    yout = y + h*(c1*dydx+c3*ak3+c4*ak4+c6*ak6)
   
    yerr = h*(dc1*dydx+dc3*ak3+dc4*ak4+dc5*ak5+dc6*ak6)
 
  end subroutine

  subroutine init_surf ()

    implicit none

    ntrial = 11_ip
    rprimer = 0.4_rp
    inuc = 0_ip
    epsilon = 1d-5 
    epsiscp = 0.08_rp

  end subroutine
                                                                        
  subroutine allocate_space_for_integ (ncent)

    use mod_memory, only: alloc
    implicit none

    integer(kind=ip) :: i
    integer(kind=ip), intent(in) :: ncent

    neqsurf = ncent
    call alloc ('mod_surf', 'insurf', insurf, ncent)
    forall (i=1:neqsurf) insurf(i) = i
    call alloc ('mod_surf', 'xyzrho', xyzrho, ncent, 3)
    call alloc ('mod_surf', 'nangleb', nangleb, ncent, 4)
    nangleb(:,1) = 1
    nangleb(:,2) = 434
    nangleb(:,3) = 0
    nangleb(:,4) = 0
    call alloc ('mod_surf', 'rmaxsurf', rmaxsurf, ncent)
    rmaxsurf = 10.0_rp
    call alloc ('mod_surf', 'rstart', rstart, ncent, maxstart)
    call alloc ('mod_surf', 'nrsearch', nrsearch, ncent)
    ! By default no explicit values of the radial coordinate are given to 
    ! the atoms in the initial search of their interatomic surfaces.
    call alloc ('mod_surf', 'lstart', lstart, ncent)
    lstart = .false.

  end subroutine allocate_space_for_integ

  subroutine deallocate_space_for_integ ()

    use mod_memory, only: free
    implicit none

    call free ('mod_surf', 'xyzrho', xyzrho)
    call free ('mod_surf', 'rmaxsurf', rmaxsurf)
    call free ('mod_surf', 'nangleb', nangleb)
    call free ('mod_surf', 'rstart', rstart)
    call free ('mod_surf', 'nrsearch', nrsearch)
    call free ('mod_surf', 'lstart', lstart)
    call free ('mod_surf', 'insurf', insurf)

  end subroutine deallocate_space_for_integ

  subroutine allocate_space_for_surface (nangular,ncutsurf)

    use mod_memory, only: alloc
    implicit none
    integer(kind=ip), intent(in) :: nangular, ncutsurf

    call alloc ('mod_surf', 'rlimsurf', rlimsurf, nangular, ncutsurf)
    call alloc ('mod_surf', 'nlimsurf', nlimsurf, nangular)

  end subroutine allocate_space_for_surface

  subroutine deallocate_space_for_surface

    use mod_memory, only: free
    implicit none

    call free ('mod_surf', 'rlimsurf', rlimsurf)
    call free ('mod_surf', 'nlimsurf', nlimsurf)

  end subroutine deallocate_space_for_surface

end module mod_surf
