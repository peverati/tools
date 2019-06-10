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
module mod_wfn
  
  use mod_prec, only: ip, rp
  implicit none
  public

  integer(kind=ip), parameter :: mgrp = 500_ip
  integer(kind=ip), parameter :: ngtoh = 21_ip
  integer(kind=ip), parameter :: maxtype = 56_ip

  ! nlm keeps the nlm values of x^n y^l z^m gaussian
  integer(kind=ip) :: nlm(maxtype,3)

  integer(kind=ip) :: nprims
  integer(kind=ip) :: nmo
  integer(kind=ip) :: ncent
  integer(kind=ip) :: maxgrp
  integer(kind=ip) :: numshells

  real(kind=rp), allocatable, dimension(:,:) :: rcutte
  real(kind=rp), allocatable, dimension(:,:) :: coeff
  integer(kind=ip), allocatable, dimension(:) :: npc
  integer(kind=ip), allocatable, dimension(:) :: ngroup
  integer(kind=ip), allocatable, dimension(:,:) :: icenat
  integer(kind=ip), allocatable, dimension(:,:) :: nzexp
  integer(kind=ip), allocatable, dimension(:,:,:) :: nuexp
  real(kind=rp), allocatable, dimension(:,:) :: coords
  real(kind=rp), allocatable, dimension(:) :: oexp
  real(kind=rp), allocatable, dimension(:,:) :: rint
  real(kind=rp), allocatable, dimension(:) :: occ
  real(kind=rp), allocatable, dimension(:) :: eorb
  real(kind=rp), allocatable, dimension(:) :: charge
  character(len=8), allocatable, dimension(:) :: atnam
  integer(kind=ip), allocatable, dimension(:) :: icen
  integer(kind=ip), allocatable, dimension(:) :: ityp

  ! OPTIONS:
  ! Any gaussian primitive and gaussian derivative will be assumed
  ! to be zero if it is smaller than cuttz.
  real(kind=rp) :: cuttz

contains

  subroutine allocate_space_for_wfn ()
  
    use mod_memory, only: alloc
    use mod_io, only: faterr, ferror
    implicit none

    integer(kind=ip) :: ier

    call alloc ('mod_wfn', 'coeff', coeff, nmo, nprims)
    call alloc ('mod_wfn', 'npc', npc, ncent)
    call alloc ('mod_wfn', 'ngroup', ngroup, ncent)
    call alloc ('mod_wfn', 'icenat', icenat, nprims, ncent)
    call alloc ('mod_wfn', 'nzexp', nzexp, ncent, mgrp)
    call alloc ('mod_wfn', 'nuexp', nuexp, ncent, mgrp, ngtoH)
    call alloc ('mod_wfn', 'coords', coords, ncent, 3)
    call alloc ('mod_wfn', 'oexp', oexp, nprims)
    call alloc ('mod_wfn', 'rint', rint, ncent, ncent)
    call alloc ('mod_wfn', 'occ', occ, nmo)
    call alloc ('mod_wfn', 'eorb', eorb, nmo)
    call alloc ('mod_wfn', 'charge', charge, ncent)
    call alloc ('mod_wfn', 'icen', icen, nprims)
    call alloc ('mod_wfn', 'ityp', ityp, nprims)
    call alloc ('mod_wfn', 'rcutte', rcutte, ncent, mgrp)
    if (.not.allocated(atnam)) then
      allocate (atnam(ncent),stat=ier) 
      if (ier.ne.0) then
        call ferror('mod_wfn', 'cannot allocate atnam', faterr)
      end if
    end if
 
  end subroutine allocate_space_for_wfn

  subroutine deallocate_space_for_wfn ()

    use mod_io, only: faterr, ferror
    use mod_memory, only: free
    implicit none

    integer(kind=ip) :: ier

    call free ('mod_wfn', 'coeff', coeff)
    call free ('mod_wfn', 'npc', npc)
    call free ('mod_wfn', 'ngroup', ngroup)
    call free ('mod_wfn', 'icenat', icenat)
    call free ('mod_wfn', 'nzexp', nzexp)
    call free ('mod_wfn', 'nuexp', nuexp)
    call free ('mod_wfn', 'coords', coords)
    call free ('mod_wfn', 'oexp', oexp)
    call free ('mod_wfn', 'rint', rint)
    call free ('mod_wfn', 'occ', occ)
    call free ('mod_wfn', 'eorb', eorb)
    call free ('mod_wfn', 'charge', charge)
    call free ('mod_wfn', 'icen', icen)
    call free ('mod_wfn', 'ityp', ityp)
    call free ('mod_wfn', 'rcutte', rcutte)
    if (allocated(atnam)) then
      deallocate (atnam,stat=ier) 
      if (ier.ne.0) then
        call ferror('mod_wfn', 'cannot deallocate atnam', faterr)
      end if
    end if

  end subroutine deallocate_space_for_wfn
                                                                        
  subroutine init_wfn()
                
    implicit none

    cuttz = 1d-14

    !.p's
    nlm(2,1)=1      !px
    nlm(2,2)=0      !px
    nlm(2,3)=0      !px
    nlm(3,1)=0      !py
    nlm(3,2)=1      !py
    nlm(3,3)=0      !py
    nlm(4,1)=0      !pz
    nlm(4,2)=0      !pz
    nlm(4,3)=1      !pz
    !.d's
    nlm(5,1)=2      !xx
    nlm(6,2)=2      !yy
    nlm(7,3)=2      !zz
    nlm(8,1)=1      !xy
    nlm(8,2)=1
    nlm(9,1)=1      !xz
    nlm(9,3)=1
    nlm(10,2)=1     !yz
    nlm(10,3)=1
    !.f's
    nlm(11,1)=3     !xxx
    nlm(12,2)=3     !yyy
    nlm(13,3)=3     !zzz    
    nlm(14,1)=2     !xxy
    nlm(14,2)=1
    nlm(15,1)=2     !xxz
    nlm(15,3)=1
    nlm(16,2)=2     !yyz
    nlm(16,3)=1 
    nlm(17,1)=1     !xyy
    nlm(17,2)=2 
    nlm(18,1)=1     !xzz
    nlm(18,3)=2     
    nlm(19,2)=1     !yzz
    nlm(19,3)=2
    nlm(20,1)=1     !xyz
    nlm(20,2)=1 
    nlm(20,3)=1 
    !.g's
    nlm(21,1)=4     !xxxx
    nlm(22,2)=4     !yyyy
    nlm(23,3)=4     !zzzz
    nlm(24,1)=3     !xxxy
    nlm(24,2)=1
    nlm(25,1)=3     !xxxz
    nlm(25,3)=1 
    nlm(26,1)=1     !xyyy
    nlm(26,2)=3 
    nlm(27,2)=3     !yyyz
    nlm(27,3)=1 
    nlm(28,1)=1     !xzzz 
    nlm(28,3)=3 
    nlm(29,2)=1     !yzzz
    nlm(29,3)=3
    nlm(30,1)=2     !xxyy
    nlm(30,2)=2 
    nlm(31,1)=2     !xxzz
    nlm(31,3)=2
    nlm(32,2)=2     !yyzz 
    nlm(32,3)=2 
    nlm(33,1)=2     !xxyz 
    nlm(33,2)=1
    nlm(33,3)=1
    nlm(34,1)=1     !xyyz
    nlm(34,2)=2 
    nlm(34,3)=1
    nlm(35,1)=1     !xyzz
    nlm(35,2)=1
    nlm(35,3)=2
    !.h's
    nlm(36,1)=0
    nlm(36,2)=0
    nlm(36,3)=5
    nlm(37,1)=0
    nlm(37,2)=1
    nlm(37,3)=4
    nlm(38,1)=0
    nlm(38,2)=2
    nlm(38,3)=3
    nlm(39,1)=0
    nlm(39,2)=3
    nlm(39,3)=2
    nlm(40,1)=0
    nlm(40,2)=4
    nlm(40,3)=1
    nlm(41,1)=0
    nlm(41,2)=5
    nlm(41,3)=0
    nlm(42,1)=1
    nlm(42,2)=0
    nlm(42,3)=4
    nlm(43,1)=1
    nlm(43,2)=1
    nlm(43,3)=3
    nlm(44,1)=1
    nlm(44,2)=2
    nlm(44,3)=2
    nlm(45,1)=1
    nlm(45,2)=3
    nlm(45,3)=1
    nlm(46,1)=1
    nlm(46,2)=4
    nlm(46,3)=0
    nlm(47,1)=2
    nlm(47,2)=0
    nlm(47,3)=3
    nlm(48,1)=2
    nlm(48,2)=1
    nlm(48,3)=2
    nlm(49,1)=2
    nlm(49,2)=2
    nlm(49,3)=1
    nlm(50,1)=2
    nlm(50,2)=3
    nlm(50,3)=0
    nlm(51,1)=3
    nlm(51,2)=0
    nlm(51,3)=2
    nlm(52,1)=3
    nlm(52,2)=1
    nlm(52,3)=1
    nlm(53,1)=3
    nlm(53,2)=2
    nlm(53,3)=0
    nlm(54,1)=4
    nlm(54,2)=0
    nlm(54,3)=1
    nlm(55,1)=4
    nlm(55,2)=1
    nlm(55,3)=0
    nlm(56,1)=5
    nlm(56,2)=0
    nlm(56,3)=0
  
  end subroutine init_wfn

  subroutine end_wfn()

    implicit none

    call deallocate_space_for_wfn ()

  end subroutine end_wfn

  subroutine rdwfn (wfnfile)

    use mod_memory, only: alloc, free
    use mod_io, only: ferror, faterr, mline, udat, string, warning, uout
    implicit none
 
    ! Arguments
    character(len=*), intent(in) :: wfnfile
 
    ! Local vars
    integer(kind=ip) :: i, iwfn, j, k 
    real(kind=rp) :: tote, x1, x2, y1, y2, z1, z2, dis, gamma
    character(len=80) :: wfnttl
    character(len=4) :: mode
    character(len=17) :: label
    character(len=8) :: check

    ! Init data
    open (udat,file=wfnfile,status='old') 
    iwfn = udat
    read (iwfn,101) wfnttl
    read (iwfn,102) mode, nmo, nprims, ncent
    write (uout,'(1x,a,1x,i0)') string('# Number of centers :'), ncent
    call allocate_space_for_wfn ()
    do i = 1,ncent
      read (iwfn,103) atnam(i),j,(coords(j,k),k=1,3),charge(j)
      atnam = adjustl(atnam)
      write (uout,'(1x,a,1x,i0,1x,a,1x,f4.1,1x,3f12.6)') '#', i, &
                    atnam(i)(1:4), charge(i), coords(i,:)
    end do

    ! Evaluate internuclear distances
    do i = 1,ncent
      x1 = coords(i,1)
      y1 = coords(i,2)
      z1 = coords(i,3)
      do j = 1,i
        x2 = coords(j,1)
        y2 = coords(j,2)
        z2 = coords(j,3)
        dis = sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
        rint(i,j) = dis
        rint(j,i) = dis
      end do
    end do

    ! Read 
    read (iwfn,104) (icen(i),i=1,nprims)
    read (iwfn,104) (ityp(i),i=1,nprims)
    read (iwfn,105) (oexp(i),i=1,nprims)
    do i = 1,nprims
      if (ityp(i).gt.maxtype) then
        call ferror('rdwfn', 'cannot work with , i- or higher primitives', faterr)
      end if
    end do
    do i = 1,nmo
      read (iwfn,106) occ(i), eorb(i) 
      read (iwfn,107) (coeff(i,j),j=1,nprims)
    end do
    read (iwfn,108) check
    if (check .ne. 'END DATA') then
      call ferror('rdwfn', 'end card 1 not found', faterr)
    endif

    ! Reduce and set info
    call setupwfn ()
    call filtergto ()

    ! Special cases
    read (iwfn,109) label,tote,gamma

    close (iwfn)

    ! formats
101 format (a80)
102 format (4x,a4,10x,3(i5,15x))
103 format (a8,11x,i3,2x,3f12.8,10x,f5.1)
104 format (20x,20i3)
105 format (10x,5e14.7)
106 format (35x,f12.8,15x,f12.8)
107 format (5e16.8)
108 format (a8)
109 format (a17,f20.12,18x,f13.8)
 
  end subroutine

  subroutine setupwfn ()

    use mod_prec, only: rp, ip
    use mod_memory, only: alloc, free
    use mod_io, only: ferror, faterr, uout, string
    implicit none

    real(kind=rp) :: alph
    integer(kind=ip) :: npa, i, ic, icentro, inda, isud
    integer(kind=ip) :: isuk, itd, itip, itk, j, k, m, npcant
    real(kind=rp), allocatable, dimension(:) :: oexpa
    real(kind=rp), allocatable, dimension(:,:) :: coefa
    integer(kind=ip), allocatable, dimension(:) :: icena
    integer(kind=ip), allocatable, dimension(:) :: itypa

    ! Init
    call alloc ('setupwfn', 'oexpa', oexpa, nprims)
    call alloc ('setupwfn', 'coefa', coefa, nmo, nprims)
    call alloc ('setupwfn', 'icena', icena, nprims)
    call alloc ('setupwfn', 'itypa', itypa, nprims)
    npa = 1_ip
    icena(1) = icen(1)
    oexpa(1) = oexp(1)
    itypa(1) = ityp(1)
    coefa(1:nmo,1) = coeff(1:nmo,1)
    cyclej: do j= 2,nprims
      do m = 1,npa
        if (icen(j).eq.icena(m) .and. &
            ityp(j).eq.itypa(m) .and. &
            abs(oexp(j)-oexpa(m)).le.1d-10) then
          coefa(1:nmo,m) = coefa(1:nmo,m)+coeff(1:nmo,j)
          cycle cyclej
        end if
      end do
      npa = npa + 1_ip
      icena(npa) = icen(j)
      oexpa(npa) = oexp(j)
      itypa(npa) = ityp(j)
      coefa(1:nmo,npa) = coeff(1:nmo,j)
    end do cyclej

    ! Recompute the original variables
    write (uout,'(1x,a,1x,i0)') string('# Number of molecular orbitals'), nmo
    write (uout,'(1x,a,1x,i0,1x,a,1x,i0)') string('# Input number of primitives'), &
                                           nprims, 'reduced to', npa

    nprims = npa
    do j = 1,nprims
      icen(j) = icena(j)
      oexp(j) = oexpa(j)
      ityp(j) = itypa(j)
      coeff(1:nmo,j) = coefa(1:nmo,j)
    end do

    ! Determine primitives corresponding to each center.
    do ic = 1,ncent
      npc(ic) = 0_ip
    end do
    do j = 1,nprims
      ic = icen(j)
      npc(ic) = npc(ic) + 1_ip
      inda = npc(ic)
      icenat(inda,ic) = j
    end do

    ! Classify primitives in each center by types and similar exponents.
    do ic = 1,ncent
      scyclej: do j = 1,npc(ic)
        k = icenat(j,ic)
        itk = ityp(k)
        isuk = nlm(itk,1) + nlm(itk,2) + nlm(itk,3)
        if (j.eq.1) then
          ngroup(ic) = 1_ip
          nzexp(ic,1) = 1_ip
          nuexp(ic,1,1) = k
        else
          do m = 1,ngroup(ic)
            inda = nuexp(ic,m,1)
            itd = ityp(inda)
            isud = nlm(itd,1) + nlm(itd,2) + nlm(itd,3)
            if (abs(oexp(k)-oexp(inda)).lt.1d-8) then
              if (itk.eq.1.and.itd.eq.1) then
                call ferror ('setupwfn', 'two s primitives with equal exponents', faterr)
              else
                if (isuk.eq.isud) then
                  nzexp(ic,m) = nzexp(ic,m) + 1_ip
                  nuexp(ic,m,nzexp(ic,m)) = k
                  cycle scyclej
                end if
              end if
            end if
          end do
          ngroup(ic) = ngroup(ic) + 1_ip
          if (ngroup(ic).gt.mgrp) then
            call ferror ('setupwfn', 'increase mgrp in mod_wfn file', faterr)
          end if
          nzexp(ic,ngroup(ic)) = 1_ip
          nuexp(ic,ngroup(ic),1) = k
        end if   
      end do scyclej
    end do
  
    ! Reconstruct the values of ityp(),icen(),oexp(), and coef()
    i = 0_ip
    do ic = 1,ncent
      do m = 1,ngroup(ic)
        do k = 1,nzexp(ic,m)
          j = nuexp(ic,m,k)
          alph = oexp(j)
          itip = ityp(j)
          icentro = icen(j)
          i = i + 1_ip
          itypa(i) = itip
          oexpa(i) = alph
          icena(i) = icentro
          coefa(1:nmo,i) = coeff(1:nmo,j)
        end do
      end do
    end do
    
    call free ('setupwfn', 'coeff', coeff)
    call alloc ('setupwfn', 'coeff', coeff, nmo, nprims)
    do i = 1,nprims
      ityp(i) = itypa(i)
      oexp(i) = oexpa(i)
      icen(i) = icena(i)
      coeff(1:nmo,i) = coefa(1:nmo,i)
    end do
   
    npcant = 0_ip
    do ic = 1,ncent
      do k = 1,npc(ic)
        icenat(k,ic) = k + npcant
      end do
      npcant = npcant + npc(ic)
    end do
  
    ! Reconstruct the values of nuexp(). Now, they are ordered.
    ! Determine also which is the maximum value of ngroup(ic)
    ! Determine the total number of shells.
    i = 0_ip
    numshells = 0_ip
    maxgrp = 0_ip
    do ic = 1,ncent
      numshells = numshells + ngroup(ic)
      if (ngroup(ic).gt.maxgrp) maxgrp = ngroup(ic)
      do m = 1,ngroup(ic)
        do k = 1,nzexp(ic,m)
          i = i + 1_ip
          nuexp(ic,m,k) = i
        end do
      end do
    end do

    ! Deallocate arrays.
    call free ('setupwfn', 'oexpa', oexpa)
    call free ('setupwfn', 'coefa', coefa)
    call free ('setupwfn', 'icena', icena)
    call free ('setupwfn', 'itypa', itypa)
  
  end subroutine

  subroutine filtergto ()

    use mod_io, only: uout
    use mod_param, only: verbose
    implicit none
 
    real(kind=rp) :: zz, x1
    integer(kind=ip) :: ic, i, k, lsum, m
    character(len=1) :: lb(0:5), lbl
    data (lb(i),i=0,5) /'S','P','D','F','G','H'/

    ! Maximum distance at which it is necessary to compute a shell.
    write (uout,'(1x,a,1x,e17.10)') '# Cutoff fot GTOS, eps =', cuttz
    do ic = 1,ncent
      if (verbose) write (uout,'(1x,a,1x,i0)') '# Center', ic
      do m = 1,ngroup(ic)
        i = nuexp(ic,m,1)
        lsum = nlm(ityp(i),1)+nlm(ityp(i),2)+nlm(ityp(i),3)
        zz = oexp(i)
        x1 = 0.1_rp
        do 
          if (x1**lsum*exp(-zz*x1*x1).le.abs(cuttz)) exit
          x1 = x1 + 0.1_rp
        end do
        rcutte(ic,m) = x1
        if (verbose) then
          lbl = lb(lsum)
          write (uout,613) lbl,zz,x1,(nuexp(ic,m,k),k=1,nzexp(ic,m))
        end if
      end do
    end do
 
    do ic = 1,ncent
      do m = 1,ngroup(ic)
        rcutte(ic,m) = rcutte(ic,m)*rcutte(ic,m)
      end do
    end do

613 format (1x,'# ',a,' Shell Exp = ',e16.8,4x,'Cutoff = ',f13.6,4x,'Primitives : ',21(1x,i0))

  end subroutine

end module mod_wfn
