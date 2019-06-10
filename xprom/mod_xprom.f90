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
module mod_xprom

  implicit none
  public
      
  integer(kind=4), parameter :: mc = 120
  integer(kind=4), parameter :: mcmc = mc*(mc+1)/2

  logical :: nodef
  logical :: t2
  logical :: nosym

  integer(kind=4) :: neq, ncent
  character(len=2) :: symbol(mc)
  real(kind=8) :: xyz(mc,3), charge(mc), elec(mc), qq(mc)
  integer(kind=4) :: ig(mc,2)
  real(kind=8) :: nrepenergy

  integer(kind=4), dimension(mc) :: ineq, mult
  integer(kind=4), dimension(mc,2) :: idx
  integer(kind=4), dimension(mcmc,2) :: ipairc
  integer(kind=4), dimension(mcmc) :: icur,jcur

  integer(kind=4) :: ngroup(2)
  real(kind=8) :: def(mc)

  real(kind=8), dimension(mc) :: kin, pot, ee, eec, eexc, corr
  real(kind=8), dimension(mc) :: own, net, inter, eadd, eeff
  real(kind=8), dimension(mc,mc) :: nnab, enab, neab, eeab, corrab
  real(kind=8), dimension(mc,mc) :: intab, coulab, xcab, xvclab

  real(kind=8), dimension(mc) :: fkin, fpot, fee, feec, feexc, fcorr
  real(kind=8), dimension(mc) :: fown, fnet, finter, feadd, feeff
  real(kind=8), dimension(mc,mc) :: fnnab, fenab, fneab, feeab
  real(kind=8), dimension(mc,mc) :: fintab, fcoulab, fxcab, fxvclab, fcorrab

contains

  subroutine full ()

    implicit none
  
  end subroutine full

  subroutine init_xprom ()
  
    implicit none

    nodef = .true.
    t2 = .false.
    nosym = .true.

  end subroutine init_xprom

  subroutine end_xprom ()
  
    implicit none

  end subroutine end_xprom

  subroutine loadpmd ()      

    use mod_io, only: uout, faterr, ferror, mline, &
                      isreal, isword, equal, string
    use mod_datatm, only: maxzat, namatm
    implicit none

    integer(kind=4) :: lread, i, j, lp, ncompute
    integer(kind=4) :: iol1, iol2, n, icn, jcn
    character(len=mline) :: pmdout, pg, line
    logical :: ok

    lread = 600
    open (lread,file='data-from-promolden')
    read (lread,'(a)',end=999) pmdout
    read (lread,'(a)',end=999) pg
    write (uout,'(1x,a,1x,a)') string('# Output file :'), string(pmdout)
    write (uout,'(1x,a,1x,a)') string('# Point group :'), string(pg)
    read (lread,*,end=999) neq
    read (lread,*,end=999) ncent
    if (neq.gt.mc .or. ncent.gt.mc) then 
      call ferror('loadpmd','increase the value of the mc parameter', faterr)
    end if
    if (neq.ne.ncent) then 
      call ferror('loadpmd','symmetry not yet allowed', faterr)
    end if
    write (uout,'(1x,a,1x,2(1x,i0))') string('# ncent & neq : '), ncent, neq
    do i = 1,ncent
      lp = 1
      read (lread,'(a)',end=999) line
      do j = 1,3
        ok = isreal(xyz(i,j),line,lp)
      end do
      ok = isword(symbol(i),line,lp)
      do j = 1,maxzat
        if (equal(symbol(i),namatm(j))) charge(i) = real(j,8) 
      end do
      write (uout,'(1x,a,1x,a,3x,3(1x,f10.6))') '#', symbol(i), xyz(i,:)
    end do

    do i = 1,ncent
      ineq(i) = i
      mult(i) = 1
      idx(i,1) = i
    end do
    ! For all possible pairs of atoms, determine the indices of the
    ! equivalent pairs that are actually computed.
    ncompute = 0
    do i = 1,ncent
      do j = i+1,ncent
       ncompute = ncompute + 1
       ipairc(ncompute,1) = i
       ipairc(ncompute,2) = j
      end do
    end do

    ! By now no symmetry, this is redundant now
    do i = 1,ncompute
      iol1 = ipairc(i,1)
      iol2 = ipairc(i,2)
      icur(i) = ipairc(i,1)
      jcur(i) = ipairc(i,2)
    end do

    ! Read monocentric terms
    do i = 1,neq
      read (lread,*,end=999) kin(ineq(i))
    end do
    do i = 1,neq
      read (lread,*,end=999) pot(ineq(i))
    end do
    do i = 1,neq
      read (lread,*,end=999) ee(ineq(i))
    end do
    do i = 1,neq
      read (lread,*,end=999) eec(ineq(i))
    end do
    do i = 1,neq
      read (lread,*,end=999) eexc(ineq(i))
    end do
    do i = 1,neq
    read (lread,*,end=999) own(ineq(i))
    end do
    do i = 1,neq
      read (lread,*,end=999) net(ineq(i))
    end do
    do i = 1,neq
      read (lread,*,end=999) inter(ineq(i))
    end do
    do i = 1,neq
      read (lread,*,end=999) eadd(ineq(i))
    end do
    do i = 1,neq
      read (lread,*,end=999) eeff(ineq(i))
    end do
    if (t2) then
      do i = 1,neq
        read (lread,*,end=999) corr(ineq(i))
      end do
    end if

    ! Read and setup bicentric terms
    do i = 1,neq
      do j = 1,ncent
        if (j.ne.ineq(i)) then
         read (lread,*,end=999) nnab(ineq(i),j),enab(ineq(i),j), &
               neab(ineq(i),j),eeab(ineq(i),j),intab(ineq(i),j)
        end if
      end do
    end do
    do i = 1,neq
      do j = 1,ncent
        if (j.ne.ineq(i)) then
          read (lread,*,end=999) coulab(ineq(i),j),xcab(ineq(i),j)
          xvclab(ineq(i),j) = intab(ineq(i),j) - xcab(ineq(i),j)
        end if
      end do
    end do
    if (t2) then
      do i = 1,neq
        do j = 1,ncent
          if (j.ne.ineq(i)) then
            read (lread,*,end=999) corrab(ineq(i),j)
          end if
        end do
      end do
    end if
    do j = 1,ncent
      read (lread,*,end=999) qq(j)
      elec(j) = charge(j) - qq(j)
    end do

    n = 0
    do i = 1,ncent
      do j = i+1,ncent
        n = n + 1
        icn=icur(n)
        jcn=jcur(n)
        fnnab(i,j) = nnab(icn,jcn)
        feeab(i,j) = eeab(icn,jcn)
        fintab(i,j) = intab(icn,jcn)
        fcoulab(i,j) = coulab(icn,jcn)
        fxcab(i,j) = xcab(icn,jcn)
        fxvclab(i,j) = xvclab(icn,jcn)
        fneab(i,j) = neab(icn,jcn)
        fenab(i,j) = enab(icn,jcn)
        if (t2) then
          fcorrab(i,j) = corrab(icn,jcn)
        end if
        ! Symmetric
        fnnab(j,i) = nnab(icn,jcn)
        feeab(j,i) = eeab(icn,jcn)
        fintab(j,i) = intab(icn,jcn)
        fcoulab(j,i) = coulab(icn,jcn)
        fxcab(j,i) = xcab(icn,jcn)
        fxvclab(j,i) = xvclab(icn,jcn)
        fneab(j,i) = neab(icn,jcn)
        fenab(j,i) = enab(icn,jcn)
        if (t2) then
          fcorrab(j,i) = corrab(icn,jcn)
        end if
      end do
    end do

!999 call ferror('loadpmd', 'wrong file format', faterr)
999 continue

  end subroutine

  ! Print a resume of the output file
  subroutine resume ()

    use mod_io, only: uout
    implicit none

    integer(kind=4) :: i, j, k

    write (uout,*)
    write (uout,*) 'Monocentric energy components'
    write (uout,*) '-----------------------------'
    write (uout,11)
    do i = 1,ncent
      k = idx(i,1)
      fkin(i) = kin(ineq(k))
      fpot(i) = pot(ineq(k))
      fee(i) = ee(ineq(k))
      feec(i) = eec(ineq(k))
      feexc(i) = eexc(ineq(k))
      write (uout,13) i,fkin(i),fpot(i),fee(i),feec(i),feexc(i)
    end do
    write (uout,*)
    write (uout,12)
    do i = 1,ncent
      k = idx(i,1)
      fnet(i) = net(ineq(k))
      finter(i)= inter(ineq(k))
      feadd(i) = eadd(ineq(k))
      feeff(i) = eeff(ineq(k))
      fown(i) = own(ineq(k))
      write (uout,13) i,fnet(i),finter(i),feadd(i),feeff(i),fown(i)
    end do
    if (t2) then
      write (uout,*)
      write (uout,92)
      do i = 1,ncent
        k = idx(i,1)
        fcorr(i) = corr(ineq(k))
        write (uout,13) i,fcorr(i),qq(i)
      end do
    end if
    write (uout,*)
    write (uout,'(1x,a)') "Resume of some sumed terms"
    write (uout,'(1x,a)') "--------------------------------------------"
    write (uout,'(1x,a,f16.6)') "Sum of charges       : ", sum(qq) 
    write (uout,'(1x,a,f16.6)') "Sum of Kin intra     : ", sum(fkin) 
    write (uout,'(1x,a,f16.6)') "Sum of Nuclear own   : ", sum(fown) 
    write (uout,'(1x,a,f16.6)') "Sum of Pot intra     : ", sum(fpot) 
    write (uout,'(1x,a,f16.6)') "Sum of Coulomb intra : ", sum(feec) 
    write (uout,'(1x,a,f16.6)') "Sum of XC intra      : ", sum(feexc) 
    if (t2) write (uout,'(1x,a,f16.6)') "Sum of E_corr intra  : ", sum(fcorr) 
    write (uout,'(1x,a,f16.6)') "Sum of Enet intra    : ", sum(fnet) 
    write (uout,'(1x,a,f16.6)') "Sum of Einter intra  : ", sum(finter)*0.5 

    write (uout,*)
    write (uout,*) 'Bicentric energy components'
    write (uout,*) '---------------------------'
    write (uout,501)
    do i = 1,ncent
      write (uout,601) (i,j,feeab(i,j),j=1,ncent)
    end do
    write (uout,502)
    do i = 1,ncent
      write (uout,601) (i,j,fintab(i,j),j=1,ncent)
    end do
    write (uout,503)
    do i = 1,ncent
      write (uout,601) (i,j,fcoulab(i,j),j=1,ncent)
    end do
    write (uout,504)
    do i = 1,ncent
      write (uout,601) (i,j,fxcab(i,j),j=1,ncent)
    end do
    write (uout,505)
    do i = 1,ncent
      write (uout,601) (i,j,fxvclab(i,j),j=1,ncent)
    end do
    if (t2) then
      write (uout,506)
      do i = 1,ncent
        write (uout,601) (i,j,fcorrab(i,j),j=1,ncent)
      end do
      write (uout,509)
    end if
    write (uout,*)
    write (uout,'(1x,a)') "Resume of some sumed terms"
    write (uout,'(1x,a)') "--------------------------------------------"
    write (uout,'(1x,a,f16.6)') "Sum of NE+EN inter   : ", sum(fneab)*0.5+&
                                                           sum(fenab)*0.5
    write (uout,'(1x,a,f16.6)') "Sum of Coulonb inter : ", sum(fcoulab)*0.5
    write (uout,'(1x,a,f16.6)') "Sum of XC inter      : ", sum(fxcab)*0.5
    write (uout,'(1x,a,f16.6)') "Sum of Classic inter : ", sum(fxvclab)*0.5
    write (uout,'(1x,a,f16.6)') "Sum of XC + Classic  : ", sum(fintab)*0.5
    if (t2) then
      write (uout,'(1x,a,f16.6)') "Sum of E_corr inter  : ", sum(fcorrab)*0.5
    end if
  
    write (uout,*)
    write (uout,'(1x,a)') "Resume of total energy terms"
    write (uout,'(1x,a)') "--------------------------------------------"
    write (uout,'(1x,a,f16.6)') "*** Total Nuclear rep. : ", sum(fnnab)*0.5
    write (uout,'(1x,a,f16.6)') "*** Total Nuc-Electron : ", sum(fown)+&
                                                             sum(fneab)*0.5+&
                                                             sum(fenab)*0.5
    write (uout,'(1x,a,f16.6)') "*** Total Kinetic      : ", sum(fkin)
    write (uout,'(1x,a,f16.6)') "*** Total Coulomb      : ", sum(feec)+&
                                                             sum(fcoulab)*0.5
    write (uout,'(1x,a,f16.6)') "*** Total Exchange     : ", (sum(feexc)+&
                                                             sum(fxcab)*0.5)
    if (t2) then
      write (uout,'(1x,a,f16.6)') "*** Total Correlation  : ", sum(fcorr)+&
                                                               sum(fcorrab)*0.5
      write (uout,'(1x,a,f16.6)') "*** Total ENERGY       : ", sum(fintab)*0.5+&
                                                               sum(fnet)
    end if
    write (uout,*)

501 format (1x,'Interatomic EE repulsions')
502 format (1x,'Interatomic total interactions')
503 format (1x,'Interatomic Coulomb interactions')
504 format (1x,'Interatomic XC interactions')
505 format (1x,'Interatomic Classic interactions')
506 format (1x,'Interatomic Correlation interactions')
509 format (1x,'Total Delocalization Indexes')
601 format (4(1x,'(',2I4,' ) = ',F13.6))
13  format (I4,11(1x,F13.6))
11  format (' Atom',7x,'Kin',11x,'Pot',11x,'EErep',9x,&
            'EECoul',8x,'EExc',/1x,73('-'))
12  format (' Atom',7x,'Enet',10x,'Einter',8x,'Eadd',10x,&
            'Eeff',6x,'el-own-nuc'/1x,73('-'))
92  format (' Atom',6x,'E_corr',8x,'Charge',/1x,73('-'))

  end subroutine

end module mod_xprom
