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
module mod_wfn

  implicit none
  public
      
  integer(kind=4), parameter :: maxat = 1000
  integer(kind=4), parameter :: maxprim = 8000
  integer(kind=4), parameter :: maxorb = 1000
  integer(kind=4), parameter :: maxtype = 35
  integer(kind=4), parameter :: mcent=maxat
  integer(kind=4), parameter :: mgrp = 500
      
  integer(kind=4) :: nprims, nmo, ncent, npor

  integer(kind=4) :: npc(mcent) 
  integer(kind=4) :: icenat(maxprim,mcent)
  integer(kind=4) :: ngroup(mcent)
  integer(kind=4) :: nzexp(mcent,mgrp)
  integer(kind=4) :: nuexp(mcent,mgrp,15)
     
  real(kind=8) :: xnorm(maxprim)
  real(kind=8) :: xyz(maxat,3),oexp(maxprim), charge(maxat)  
  real(kind=8) :: coef(maxorb,maxprim), occ(maxprim),eorb(maxorb) 
  integer(kind=4) :: icen(maxprim), ityp(maxprim), nlm(35,3) 

  character(len=8) :: atnam(maxat)
  character(len=4) :: lgto(35)

contains

  subroutine rdwfn (ichange,wfnfile)

    use mod_io, only: faterr, ferror, mline, uin
    use mod_param, only: pi, facd
    implicit none

    integer(kind=4), intent(in) :: ichange(maxorb+maxorb)
    character(len=*), intent(in) :: wfnfile

    real(kind=8) :: aux, coefmod, zz, dblisu2
    real(kind=8) :: oexpa(maxprim), coefa(maxorb,maxprim), tote, gamma
    integer(kind=4) :: nsum, nne, icpr, ik, imold, incar, it1, intip
    integer(kind=4) :: it1d, it2, it2d, it3, it3d, itip, kmas, lui
    integer(kind=4) :: icena(maxprim), itypa(maxprim), itr(maxprim)
    integer(kind=4) :: nmotot, m, npa, ic, inda, isuk, isud, numorbitals
    integer(kind=4) :: p2m(35), m2p(35), ngtos(0:4,2), i, j, k, iwfn
    character(len=1) :: iotipe(0:4)
    character(len=mline) :: moldeninp
    data  iotipe(0) /'s'/, iotipe(1) /'p'/, iotipe(2) /'d'/, &
          iotipe(3) /'f'/, iotipe(4) /'g'/
    logical :: icorr, equalexp, eqtype

    character(len=80) :: wfnttl
    character(len=4) :: mode
    character(len=17) :: label
    character(len=8) :: check
 
    ! Type "i" in molden corresponds to type "p2m(i)" in promolden
!   data (p2m(i),i=1,35)             ! These 5 lines must substitute
!   & /1, 2, 3, 4, 5, 6, 7, 8, 9,10, ! the following 5 non-commented
!   & 11,12,13,17,14,15,18,19,16,20, ! lines in case that the WFN
!   & 35,25,21,34,33,29,24,26,22,32, ! file comes from 'gamess02'.
!   & 30,23,31,28,27/                !                                     
    data (p2m(i),i=1,35)             &
      /1, 2, 3, 4, 5, 6, 7, 8, 9,10, &
      11,12,13,17,14,15,18,19,16,20, &
      21,22,23,24,25,26,27,28,29,30, &
      31,32,33,34,35/
    ! Type "i" in promolden corresponds to type "m2p(i)" in molden
!   data (m2p(i),i=1,35)             ! These 5 lines must substitute
!   & /1, 2, 3, 4, 5, 6, 7, 8, 9,10, ! the following 5 non-commented   
!   & 11,12,13,15,16,19,14,17,18,20, ! lines in case that the WFN      
!   & 23,29,32,27,22,28,35,34,26,31, ! file comes from 'gamess02'.     
!   & 33,30,25,24,21/                !                                
    data (m2p(i),i=1,35)             &
      /1, 2, 3, 4, 5, 6, 7, 8, 9,10, &
      11,12,13,15,16,19,14,17,18,20, &
      21,22,23,24,25,26,27,28,29,30, &
      31,32,33,34,35/ 
 
    icorr = .false.
    iwfn = uin
    rewind (iwfn)

    read (iwfn,101) wfnttl
    read (iwfn,102) mode, nmo, nprims, ncent
    nmotot = nmo+nmo
    npor = nprims
    if (nmo .gt. maxorb) then
      call ferror ('rdwfn', 'max num. orbitals exceeded', faterr)
    end if
    if (nprims .gt. maxprim) then
      call ferror ('rdwfn', 'max num. primitives exceeded', faterr)
    end if
    if (ncent .gt. maxat) then
      call ferror ('rdwfn', 'max num. atoms exceeded', faterr)
    end if
 
    do i = 1,ncent
      read (iwfn,103) atnam(i),j,(xyz(j,k),k=1,3),charge(j)
    end do
    read (iwfn,104) (icen(i),i=1,nprims)
    read (iwfn,104) (ityp(i),i=1,nprims)
    read (iwfn,105) (oexp(i),i=1,nprims)
    do i = 1,nprims
      if (ityp(i).gt.maxtype) then
        call ferror ('rdwfn', 'cannot work with , h- or higher primitives', faterr)
      end if
    end do
    do i = 1,nmo
      occ(i)=0d0
      read (iwfn,106) occ(i),eorb(i) 
      read (iwfn,107) (coef(i,j),j=1,nprims)
    end do
    read (iwfn,108) check
    if (check .ne. 'END DATA') then
      call ferror ('rdwn', 'card not found', faterr)
    end if
    read (iwfn,109) label,tote,gamma

    if (label(1:5).eq.'MCSCF' .or. label(1:5).eq.'ALDET' &
                              .or. label(1:5).eq.'GENCI') then
      read (iwfn,*) label
      icorr = .true.
      do i = nmo+1,nmo+nmo
        occ(i) = 0d0
        read (iwfn,*) label
        read (iwfn,107) (coef(i,j),j=1,nprims)
      end do
      read (iwfn,108) check
      if (check .ne. 'end data') then
        call ferror ('rdwn', 'card not found', faterr)
      end if
    end if
 
    ! Since in WFN files each MO is expressed in terms of unnormalized
    ! cartesian gaussians, it could be that a given gaussian is
    ! repeated in the expansion. Here, we elliminate it adding
    ! up the coefficients.
    npa = 0
    do j = 1,nprims
      if (j.eq.1) then
        npa = npa + 1
        icena(npa) = icen(j)
        oexpa(npa) = oexp(j)
        itypa(npa) = ityp(j)
        do i = 1,nmotot
          coefa(i,npa) = coef(i,j)
        end do
      else
        do m = 1,npa
          if (icen(j).eq.icena(m) .and. ityp(j).eq.itypa(m) .and. &
                                   abs(oexp(j)-oexpa(m)).le.1d-10) then
            do i=1,nmotot
              coefa(i,m) = coefa(i,m) + coef(i,j)
            end do
            goto 10
          end if
        end do
        npa = npa + 1
        icena(npa) = icen(j)
        oexpa(npa) = oexp(j)
        itypa(npa) = ityp(j)
        do i = 1,nmotot
          coefa(i,npa) = coef(i,j)
        end do
      end if
 10 end do
 
    ! Recompute the original variables
    nprims = npa
    do j = 1,nprims
      icen(j) = icena(j)
      oexp(j) = oexpa(j)
      ityp(j) = itypa(j)
      do i = 1,nmotot
        coef(i,j) = coefa(i,j)
      end do
    end do

    ! Primitives corresponding to each center.
    do ic = 1,ncent
      npc(ic) = 0
    end do
    do j = 1,nprims
      ic = icen(j)
      npc(ic) = npc(ic) + 1
      inda = npc(ic)
      icenat(inda,ic) = j
    end do
 
    ! Classify primitives in each center by types and similar exponents.
    do ic = 1,ncent
      do j = 1,npc(ic)
        k = icenat(j,ic)
        isuk = nlm(ityp(k),1)+nlm(ityp(k),2)+nlm(ityp(k),3)
        if (j.eq.1) then
          ngroup(ic) = 1
          nzexp(ic,1) = 1
          nuexp(ic,1,1) = k
        else
          do m = 1,ngroup(ic)
            inda = nuexp(ic,m,1)
            isud = nlm(ityp(inda),1)+nlm(ityp(inda),2)+nlm(ityp(inda),3)
            equalexp = abs(oexp(k)-oexp(inda)).lt.1d-8
            if (equalexp) then
              eqtype = ityp(k).eq.1.and.ityp(inda).eq.1
              if (eqtype) then
                call ferror('rdwfn','two s gtos with equal exponents', faterr)
              else
                if (isuk.eq.isud) then
                  nzexp(ic,m)=nzexp(ic,m)+1
                  nuexp(ic,m,nzexp(ic,m))=k
                  goto 33
                end if
              end if
            end if
          end do
          ngroup(ic) = ngroup(ic) + 1
          nzexp(ic,ngroup(ic)) = 1
          nuexp(ic,ngroup(ic),1) = k
        end if   
33    end do
    end do
 
    ! First and last type in each symmetry
    ngtos(0,1) = 1
    ngtos(0,2) = 1
    ngtos(1,1) = 2
    ngtos(1,2) = 4
    ngtos(2,1) = 5
    ngtos(2,2) = 10
    ngtos(3,1) = 11
    ngtos(3,2) = 20
    ngtos(4,1) = 21
    ngtos(4,2) = 35
               
    ! Mapp the types between PROMOLDEN and MOLDEN criteria.
    ! Normalization constants of primitite GTOs are also computed.
    do ic = 1,ncent
      do m = 1,ngroup(ic)
        zz = oexp(nuexp(ic,m,1))
        do k = 1,nzexp(ic,m)
          nne = nuexp(ic,m,k)
          itip = ityp(nne)
          it1 = nlm(itip,1)
          it2 = nlm(itip,2)
          it3 = nlm(itip,3)
          it1d = it1+it1-1
          it2d = it2+it2-1
          it3d = it3+it3-1
          nsum = it1+it2+it3
          dblisu2 = (nsum+nsum+3)*0.25d0
          aux = (2d0/pi)**0.75d0*2d0**nsum*zz**dblisu2
          xnorm(nne) = aux/sqrt(facd(it1d)*facd(it2d)*facd(it3d)) 
          ! Find correspondence between PROMOLDEN and MOLDEN types.
          kmas = 0
          if (nsum.gt.0) kmas = ngtos(nsum-1,2)
          do ik = 1,nzexp(ic,m)
            itip = ityp(nuexp(ic,m,ik))
            imold = m2p(itip)
            if (imold.eq.k+kmas) itr(nne) = nuexp(ic,m,ik)
          end do
        end do
      end do
    end do
 
    ! Prepare input of MOLDEN
    moldeninp = wfnfile(1:len_trim(wfnfile))//".molden"
    lui = 55
    open (unit=lui,file=moldeninp)
    write (lui,400)
    do ic = 1,ncent
      incar = int(charge(ic))
      write (lui,401) atnam(ic)(3:4),ic,incar,(xyz(ic,k),k=1,3)
    end do
    write (lui,402)
    do ic = 1,ncent
      write (lui,403) ic
      do m = 1,ngroup(ic)
        i = ityp(nuexp(ic,m,1))
        write (lui,404) iotipe(nlm(i,1)+nlm(i,2)+nlm(i,3))
        write (lui,405) oexp(nuexp(ic,m,1))
      end do
      write (lui,406) '  '
    end do
    write (lui,406) '  '
    write (lui,407)
    numorbitals = nmo
    if (icorr) numorbitals = nmo + nmo
    do i = 1,numorbitals
      write (lui,408) eorb(i)
      if (occ(i).ge.0d0) write (lui,409) 
      if (occ(i).lt.0d0) write (lui,4091) 
      write (lui,410) occ(i)
      icpr = 0
      do ic = 1,ncent
        do m = 1,ngroup(ic)
          do k = 1,nzexp(ic,m)
            icpr = icpr+1
            nne = nuexp(ic,m,k)
            itip = ityp(nne)
            coefmod = ichange(i)*coef(i,itr(nne))/xnorm(itr(nne))
            it1 = nlm(itip,1)
            it2 = nlm(itip,2)
            it3 = nlm(itip,3)
            nsum= it1+it2+it3
            kmas = 0
            if (nsum.gt.0) kmas = ngtos(nsum-1,2)
            intip = p2m(k+kmas)
            write (lui,412) icpr,coefmod
          end do
        end do
      end do
    end do
         
400  format ('[Molden Format]'/'[Atoms] AU')
401  format (1x,a,I5,I3,3F13.6)
402  format ('[GTO]')
403  format (I3,1x,'0')
404  format (1x,a,'   1 1.00')
405  format (F18.10,'  1.00000000')
406  format (a)
407  format ('[MO]')
408  format (1x,'Ene=',F11.4)
409  format (1x,'Spin= Alpha')
4091 format (1x,'Spin= Beta ')
410  format (1x,'Occup=',F11.6)
412  format (1x,I6,F16.8)
101  format (A80)
102  format (4X,A4,10X,3(I5,15X))
103  format (A8,11X,I3,2X,3F12.8,10X,F5.1)
104  format (20X,20I3)
105  format (10X,5E14.7)
106  format (35X,F12.8,15X,F12.8)
107  format (5E16.8)
108  format (A8)
109  format (a17,F20.12,18X,F13.8)

  end subroutine

  subroutine init_wfn()

    implicit none
    integer(kind=4) :: i,j
    
    nlm = 0
    do j = 1,4
      do i = 1,35
        lgto(i)(j:j)=' '
      end do
    end do
    lgto(1)(1:1)='s'

    nlm(2,1)=1      !px
    nlm(3,2)=1      !py
    nlm(4,3)=1      !pz
    lgto(2)(1:1)='x'
    lgto(3)(1:1)='y'
    lgto(4)(1:1)='z'

    nlm(5,1)=2      !xx
    nlm(6,2)=2      !yy
    nlm(7,3)=2      !zz
    nlm(8,1)=1      !xy
    nlm(8,2)=1
    nlm(9,1)=1      !xz
    nlm(9,3)=1
    nlm(10,2)=1     !yz
    nlm(10,3)=1
    lgto(5 )(1:2)='xx'
    lgto(6 )(1:2)='yy'
    lgto(7 )(1:2)='zz'
    lgto(8 )(1:2)='xy'
    lgto(9 )(1:2)='xz'
    lgto(10)(1:2)='yz'

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
    lgto(11)(1:3)='xxx'
    lgto(12)(1:3)='yyy'
    lgto(13)(1:3)='zzz'
    lgto(14)(1:3)='xxy'
    lgto(15)(1:3)='xxz'
    lgto(16)(1:3)='yyz'
    lgto(17)(1:3)='xyy'
    lgto(18)(1:3)='xzz'
    lgto(19)(1:3)='yzz'
    lgto(20)(1:3)='xyz'

    nlm(21,3)=4     !zzzz
    nlm(22,2)=1     !yzzz
    nlm(22,3)=3
    nlm(23,2)=2     !yyzz 
    nlm(23,3)=2 
    nlm(24,2)=3     !yyyz
    nlm(24,3)=1 
    nlm(25,2)=4     !yyyy
    nlm(26,1)=1     !xzzz 
    nlm(26,3)=3 
    nlm(27,1)=1     !xyzz
    nlm(27,2)=1
    nlm(27,3)=2
    nlm(28,1)=1     !xyyz
    nlm(28,2)=2 
    nlm(28,3)=1
    nlm(29,1)=1     !xyyy
    nlm(29,2)=3 
    nlm(30,1)=2     !xxzz
    nlm(30,3)=2
    nlm(31,1)=2     !xxyz 
    nlm(31,2)=1
    nlm(31,3)=1
    nlm(32,1)=2     !xxyy
    nlm(32,2)=2 
    nlm(33,1)=3     !xxxz
    nlm(33,3)=1 
    nlm(34,1)=3     !xxxy
    nlm(34,2)=1
    nlm(35,1)=4     !xxxx
    lgto(21)(1:4)='zzzz'
    lgto(22)(1:4)='yzzz'
    lgto(23)(1:4)='yyzz'
    lgto(24)(1:4)='yyyz'
    lgto(25)(1:4)='yyyy'
    lgto(26)(1:4)='xzzz'
    lgto(27)(1:4)='xyzz'
    lgto(28)(1:4)='xyyz'
    lgto(29)(1:4)='xyyy'
    lgto(30)(1:4)='xxzz'
    lgto(31)(1:4)='xxyz'
    lgto(32)(1:4)='xxyy'
    lgto(33)(1:4)='xxxz'
    lgto(34)(1:4)='xxxy'
    lgto(35)(1:4)='xxxx'

  end subroutine

end module mod_wfn
