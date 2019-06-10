! Copyright (c) 2018 
! Jose Luis Casals Sainz <jluiscasalssainz@gmail.com>,
! Ángel Martín Pendás <angel@fluor.quimica.uniovi.es> and 
! Evelio Francisco Miguelez <evelio@uniovi.es>. 
!
! ortho is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or (at
! your option) any later version.
! 
! ortho is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
module mod_gto
  
  use mod_prec, only: ip, rp, i8=>size_t
  implicit none
  private

  real(kind=rp), public :: epsortho

  ! libcint data
  integer(kind=ip) :: CINTnatm
  integer(kind=ip) :: CINTnbas
  integer(kind=ip), allocatable :: CINTatm(:,:)
  integer(kind=ip), allocatable :: CINTbas(:,:)
  real(kind=rp), allocatable :: CINTenv(:)
  real(kind=rp), allocatable :: CINTbuff(:,:)
  real(kind=rp), allocatable :: CINTcoef(:,:)
  integer(kind=ip), dimension(2) :: CINTshls
  integer(kind=i8) :: CINTopt

  !
  ! cint:  s   x y z   xx xy xz yy yz zz   xxx xxy xxz xyy xyz xzz yyy yyz yzz zzz
  !        1   2 3 4    5  6  7  8  9 10    11  12  13  14  15  16  17  18  19  20
  ! wfn :  s   x y z   xx yy zz xy xz yz   xxx yyy zzz xxy xxz yyz xyy xzz yzz xyz
  !        1   2 3 4    5  8 10  6  7  9    11  17  20  12  13  18  14  16  19  15
  !      -----------------------------------------------------------------------------------
  !        0   0 0 0    0  2  3 -2 -2 -1     0   5   7  -2  -2   2  -3  -2   0  -5
  ! cint:  xxxx xxxy xxxz xxyy xxyz xxzz xyyy xyyz xyzz xzzz yyyy yyyz yyzz yzzz zzzz
  !         21   22   23   24   25   26   27   28   29   30   31   32   33   34   35
  ! wfn :  xxxx yyyy zzzz xxxy xxxz xyyy yyyz xzzz yzzz xxyy xxzz yyzz xxyz xyyz xyzz
  !         21   31   35   22   23   27   32   30   34   24   26   33   25   28   29
  !      -----------------------------------------------------------------------------------
  !          0    9   12   -2   -2    1    5    2    5   -6   -5    1   -8   -6   -6
  ! cint:  xxxxx xxxxy xxxxz xxxyy xxxyz xxxzz xxyyy xxyyz xxyzz xxzzz xyyyy xyyyz xyyzz
  !          36    37    38    39    40    41    42    43    44    45    46    47    48
  ! wfn :  zzzzz yzzzz yyzzz yyyzz yyyyz yyyyy xzzzz xyzzz xyyzz xyyyz xyyyy xxzzz xxyzz 
  !          56    55    54    53    52    51    50    49    48    47    46    45    44
  !      -----------------------------------------------------------------------------------
  !            
  ! cint:  xyzzz xzzzz yyyyy yyyyz yyyzz yyzzz yzzzz zzzzz
  !          49    50    51    52    53    54    55    56
  ! wfn :  xxyyz xxyyy xxxzz xxxyz xxxyy xxxxz xxxxy xxxxx 
  !          43    42    41    40    39    38    37    36
  !      -----------------------------------------------------------------------------------
  !   
  !
  integer(kind=ip), parameter, dimension(56) :: CINTmask = (/&
    0, 0,0,0, 0,2,3,-2,-2,-1, 0,5,7,-2,-2,2,-3,-2,0,-5, &
    0,9,12,-2,-2,1,5,2,5,-6,-5,1,-8,-6,-6, &
    20,18,16,14,12,10,8,6,4,2,0,-2,-4,-6,-8,-10,-12,-14,-16,-18,-20/)

  ! atm data 
  integer(kind=ip), parameter :: CHARGE_OF  = 1
  integer(kind=ip), parameter :: PTR_COORD  = 2
  integer(kind=ip), parameter :: NUC_MOD_OF = 3
  integer(kind=ip), parameter :: PTR_ZETA   = 4
  integer(kind=ip), parameter :: ATM_SLOTS  = 6
  ! bas data
  integer(kind=ip), parameter :: ATOM_OF    = 1
  integer(kind=ip), parameter :: ANG_OF     = 2
  integer(kind=ip), parameter :: NPRIM_OF   = 3
  integer(kind=ip), parameter :: NCTR_OF    = 4
  integer(kind=ip), parameter :: KAPPA_OF   = 5
  integer(kind=ip), parameter :: PTR_EXP    = 6
  integer(kind=ip), parameter :: PTR_COEFF  = 7
  integer(kind=ip), parameter :: BAS_SLOTS  = 8
  ! some info
  integer(kind=ip), parameter :: PTR_ENV_START = 20
  integer(kind=ip), parameter :: POINT_NUC     = 1

  public :: orthowfn, init_gto

contains

  subroutine orthowfn ()

    use mod_io, only: uout, string, warning, ferror
    use mod_wfn, only: nprims, nmo, coef
    use mod_param, only: debug, verbose
    implicit none

    logical :: ortho
    real(kind=rp) :: solap
    real(kind=rp), dimension(nprims, nprims) :: sprim
    real(kind=rp), dimension(nmo, nmo) :: smo
    integer(kind=ip), external :: CINTcgto_cart
    integer(kind=ip) :: p, q, di, dj, i, j, idx, jdx

    write (uout,'(1x,a,1x,f12.8)') string('# Tolerance for checks, epsortho : '), epsortho
    write (uout,'(1x,a)') string('# *** Computing overlap matrix in AO basis')
    ortho = .true.
    sprim = 0.0_rp
    call gtogto (sprim)
    if (debug) then
      write (uout,'(1x,a)') string('# *** S_ij in AO basis')
      do i = 1,nprims
        write (uout,'(2x,100(1x,f12.6))') (sprim(i,j),j=1,i) 
      end do
    end if
    write (uout,'(1x,a)') string('# *** Rotating AO integrals to MO basis')
    smo = matmul(coef,matmul(sprim,transpose(coef)))
    write (uout,'(1x,a)') string('# Checking orthogonality of canonical MOs')
    do i = 1,nmo
      do j = 1,i
        solap = smo(i,j)
        if (i.eq.j) then
          if (abs(abs(solap)-1.0) > epsortho) then
            ortho = .false.
            write (uout,'(1x,a,1x,i0,1x,a,1x,f12.6)') '# Canonical MO number', i, &
                                                      'not exactly normalized : ', solap
          end if
        else
          if (abs(solap) > epsortho) then
            ortho = .false.
            write (uout,'(1x,a,1x,i0,1x,i0,1x,a,1x,f12.6)') '# Canonical MO number', i, j, &
                                                            'not exactly orthogonal : ', solap
          end if
        end if
      end do
    end do
    if (ortho) then
      write (uout,'(1x,a)') string('# *** WFN is ok !!!!')
    else 
      call ferror ('ortho', 'check wfn file maybe wrong', warning)
    end if
    if (verbose) then
      write (uout,'(1x,a)') string('# *** S_ij in MO basis')
      do i =1,nmo
        write (uout,'(2x,100(1x,f12.6))') (smo(i,j),j=1,i)
      end do
    end if

    write (uout,'(1x,a)') string('# *** Computing overlap matrix in AO basis libcint version')
    ortho = .true.
    sprim = 0.0_rp
    call init_cint ()
    !if (debug) then
    !  write (uout,'(1x,a)') string('# *** CINTbuff in AO basis')
    !end if
    idx = 0
    !when parallel in omp CINTshls and CINTbuff should be threadprivate
    do p = 1,CINTnbas
      CINTshls(1) = p-1 
      di = CINTcgto_cart(CINTshls(1), CINTbas)
      ! parallel omp here
      jdx = 0
      do q = 1,p
        CINTshls(2) = q-1
        dj = CINTcgto_cart(CINTshls(2), CINTbas)
        call allocate_space_for_cint_buff (di,dj)
        CINTbuff = 0.0_rp
        call cint1e_ovlp_cart (CINTbuff, CINTshls, CINTatm, CINTnatm, CINTbas, CINTnbas, CINTenv)
        !if (debug) then
        !  do i = 1,di
        !    write (uout,'(2x,a,100(1x,f12.6))') "*** ", (CINTbuff(i,j),j=1,dj)
        !  end do
        !end if
        sprim(idx+1:idx+di,jdx+1:jdx+dj) = CINTbuff(:,:)
        sprim(jdx+1:jdx+dj,idx+1:idx+di) = transpose(CINTbuff(:,:))
        call deallocate_space_for_cint_buff ()
        jdx = jdx + dj
      end do
      idx = idx + di
    end do
    if (debug) then
      write (uout,'(1x,a)') string('# *** S_ij in AO basis')
      do i = 1,nprims
        write (uout,'(2x,100(1x,f12.6))') (sprim(i,j),j=1,i) 
      end do
    end if
    write (uout,'(1x,a)') string('# *** Rotating AO integrals to MO basis')
    smo = matmul(CINTcoef,matmul(sprim,transpose(CINTcoef)))
    write (uout,'(1x,a)') string('# Checking orthogonality of canonical MOs')
    do i = 1,nmo
      do j = 1,i
        solap = smo(i,j)
        if (i.eq.j) then
          if (abs(abs(solap)-1.0) > epsortho) then
            ortho = .false.
            write (uout,'(1x,a,1x,i0,1x,a,1x,f12.6)') '# Canonical MO number', i, &
                                                      'not exactly normalized : ', solap
          end if
        else
          if (abs(solap) > epsortho) then
            ortho = .false.
            write (uout,'(1x,a,1x,i0,1x,i0,1x,a,1x,f12.6)') '# Canonical MO number', i, j, &
                                                            'not exactly orthogonal : ', solap
          end if
        end if
      end do
    end do
    if (ortho) then
      write (uout,'(1x,a)') string('# *** WFN is ok !!!!')
    else 
      call ferror ('ortho', 'check wfn file maybe wrong', warning)
    end if
    if (verbose) then
      write (uout,'(1x,a)') string('# *** S_ij in MO basis')
      do i =1,nmo
        write (uout,'(2x,100(1x,f12.6))') (smo(i,j),j=1,i)
      end do
    end if
    call end_cint ()

  end subroutine orthowfn

  subroutine init_gto ()

    implicit none

    epsortho = 1e-6

  end subroutine init_gto

  !> Allocate pointers and data
  subroutine init_cint ()
 
    use mod_futils, only: mergesort
    use mod_memory, only: alloc, free
    use mod_io, only: uout, ferror, faterr
    use mod_param, only: pi, debug 
    use mod_wfn, only: ncent, charge, xyz, ityp, oexp, ngroup, nuexp, &
                       nmo, nprims, coef, numshells, nzexp
    implicit none
  
    integer(kind=ip) :: i, jj, iat, bas_id, di, idx, ish, l, nua, itipa
    integer(kind=ip) :: ic, m, k, itip, itypa(nprims)
    real(kind=rp), allocatable :: envbuf(:)
    real(kind=rp) :: factor, coefa(nmo,nprims)
    integer(kind=ip) :: offset
    data offset/0/
    
    integer(kind=ip), external :: CINTcgto_cart
    integer(kind=ip) :: iord(100)

    CINTnatm = ncent
    CINTnbas = numshells

    call alloc ('mod_cint', 'CINTatm', CINTatm, ATM_SLOTS, CINTnatm)
    call alloc ('mod_cint', 'CINTbas', CINTbas, BAS_SLOTS, CINTnbas)
    call alloc ('mod_cint', 'CINTcoef', CINTcoef, nmo, nprims)
    call alloc ('mod_cint', 'envbuf', envbuf, CINTnatm*4+CINTnbas*50*20)

    !> To avoid complex reorder in integrals, just make a reorder of the mo_coeff
    !> First move ityp to ordered positions then apply mask
    do ic = 1,ncent
      do m = 1,ngroup(ic)
        k = nuexp(ic,m,1)
        do jj = 1,nzexp(ic,m)
          iord(jj) = jj
        end do
        jj = nzexp(ic,m)
        call mergesort (ityp(k:k+jj-1),iord(1:jj),1,jj)
        do jj = 1,nzexp(ic,m)
          itypa(k+jj-1) = ityp(k+iord(jj)-1)
          coefa(:,k+jj-1) = coef(:,k+iord(jj)-1)
          itip = itypa(k+jj-1)
          idx = CINTmask(itip)
          CINTcoef(:,(k+jj-1)+idx) = coefa(:,(k+jj-1))
        end do
      end do
    end do

    !> Print basis info
    write (uout,'(1x,a)') "# *** Begin initialization of libcint data"
    write (uout,'(1x,a,i5)') "# Total number of shells : ", CINTnbas

    bas_id = 0
    do iat = 1,ncent
      ! atomic info
      CINTatm(CHARGE_OF,iat) = int(charge(iat),ip) ! atomic nuclear charge
      CINTatm(PTR_COORD,iat) = offset ! pointer to coordinates
      envbuf(offset+1) = xyz(iat,1) ! x (Bohr) 
      envbuf(offset+2) = xyz(iat,2) ! y (Bohr) 
      envbuf(offset+3) = xyz(iat,3) ! z (Bohr) 
      offset = offset + 3 ! advance pointer
      CINTatm(PTR_ZETA,iat) = POINT_NUC ! point nuclear model
      ! print atomic info
      if (debug) then
        write (uout,'(1x,a,i3)') "# Info for atom : ", iat
        write (uout,'(1x,6i8)') CINTatm(:,iat)
      end if
      ! basis info
      do ish = 1,ngroup(iat)
        bas_id = bas_id + 1 ! bas identifier
        nua = nuexp(iat,ish,1) ! pointer to shell exp
        itipa = ityp(nua) ! pointer to angular momentum
        CINTbas(ATOM_OF,bas_id) = iat - 1 ! atom_id, 0-based
        if (itipa.eq.1) then
          l = 0
          factor = 1d0
          factor = sqrt(4*pi)  !sp factor
        else if (itipa.eq.2) then
          l = 1
          factor = 1d0
          factor = sqrt((4*pi)/3) !sp factor
        else if (itipa.eq.5) then
          l = 2
          factor = 1d0
        else if (itipa.eq.11) then
          l = 3
          factor = 1d0
        else if (itipa.eq.21) then
          l = 4
          factor = 1d0
        else if (itipa.eq.36 .or. itipa.eq.56) then
          l = 5
          factor = 1d0
        else 
          call ferror('mod_cint' , 'not supported i or higher primitives', faterr)
        end if
        CINTbas(ANG_OF,bas_id) = l ! angular momentum
        CINTbas(NPRIM_OF,bas_id) = 1 ! primitives
        CINTbas(NCTR_OF,bas_id) = 1 ! contradted
        CINTbas(KAPPA_OF,bas_id) = 0 ! kappa for spinor basis
        ! exponents
        CINTbas(PTR_EXP,bas_id) = offset ! pointer to exp of primitive GTO
        envbuf(offset+1) = oexp(nua) ! orbital exponent
        offset = offset + 1 ! advance pointer due to exponent
        ! coefficients
        CINTbas(PTR_COEFF,bas_id) = offset ! pointer to contraction coefficients
        envbuf(offset+1) = 1d0*factor !factor due to sp contraction
        offset = offset + 1 ! advance pointer due to primitives
        ! print basis info
        if (debug) write (uout,'(1x,8i8)') CINTbas(:,bas_id)
      end do
    end do

    ! Finish storing info
    call alloc ('mod_cint', 'CINTenv', CINTenv, offset)
    do i = 1, offset
      CINTenv(i) = envbuf(i)
    end do
    call free ('mod_cint', 'envbuf', envbuf)
    if (debug) then
      write (uout,'(1x,a)') "# CINTEnv info "
      write (uout,'(19(1x,f12.6))') (CINTenv(i),i=1,offset)
    end if
    CINTopt = 0
    call cint2e_cart_optimizer(CINTopt, CINTatm, CINTnatm, CINTbas, CINTnbas, CINTenv)

    if (debug) then
      do i = 1,CINTnbas
        di = CINTcgto_cart(i-1, CINTbas)
        write (uout,'(1x,a,1x,i0,1x,a,1x,i0)') "# CINT Primitives in shell", i, ":", di
      end do
    end if
    write (uout,'(1x,a)') "# *** End initialization of libcint data"

  end subroutine init_cint

  ! Overlap matrix between primitive Cartesian Gaussian Functions
  subroutine gtogto (sprim)

    use mod_param, only: pi
    use mod_wfn, only: ncent, xyz, nlm, ityp, oexp, nprims, ngroup, nuexp, nzexp
    implicit none
    integer(kind=ip), parameter :: lamx = 12
 
    real(kind=rp), intent(out) :: sprim(nprims,nprims)
 
    real(kind=rp), dimension(ncent,ncent) :: ab2
    real(kind=rp) :: ax(1:3), bx(1:3), za, zb, p, pioverp, abaux
    real(kind=rp) :: prefactor, prod, xmu
    integer(kind=ip) :: i, j, k, l, m, ica, icb, nua, nub, n
    integer(kind=ip) :: itipa, itipb, la, lb, ka, kb, ma, mb

    ! ceabx() are the coefficients, except for the factor 
    ! EXP(-XMU*R_AB^2), where XMU=a*b/(a+b), that result from the 
    ! expansion of the product of two primitive cartesian Gaussian
    real(kind=rp) :: ceabx(-1:2*lamx,-1:lamx,-1:lamx,3)
 
    do ica = 1,ncent
      do icb = 1,ica
        abaux = 0.0_rp
        do j = 1,3
          abaux = abaux + (xyz(ica,j)-xyz(icb,j))**2
        end do
        ab2(ica,icb) = abaux
        ab2(icb,ica) = abaux
      end do
    end do
   
    ! Compute the electronic molecular electrostatic potential.
    do ica = 1,ncent
      ax(1:3) = xyz(ica,1:3)
      do ma = 1,ngroup(ica)   
        nua = nuexp(ica,ma,1)
        itipa = ityp(nua)
        la = nlm(itipa,1) + nlm(itipa,2) + nlm(itipa,3)
        za = oexp(nua)
        do icb = 1,ica
          bx(1:3) = xyz(icb,1:3)
          do mb = 1,ngroup(icb)
            nub = nuexp(icb,mb,1)
            itipb = ityp(nub)
            lb = nlm(itipb,1) + nlm(itipb,2) + nlm(itipb,3)
            zb = oexp(nub)
            p = za + zb
            xmu = za*zb/p
            prefactor = exp(-xmu*ab2(ica,icb))
            pioverp = pi/p
            pioverp = sqrt(pioverp*pioverp*pioverp)
            do j = 1,3
              call etijcalc (j,lamx,la,lb,ceabx,za,zb,ax(j),bx(j))
            end do
            ! Compute the target functions for all the products of Gaussians
            do ka = 1,nzexp(ica,ma)
              nua = nuexp(ica,ma,ka)
              itipa = ityp(nua)
              i = nlm(itipa,1)
              k = nlm(itipa,2)
              m = nlm(itipa,3)
              do kb = 1,nzexp(icb,mb)
                nub = nuexp(icb,mb,kb)
                if (nua.ge.nub) then
                  itipb = ityp(nub)
                  j = nlm(itipb,1)
                  l = nlm(itipb,2)
                  n = nlm(itipb,3)
                  prod = ceabx(0,i,j,1)*ceabx(0,k,l,2)*ceabx(0,m,n,3)
                  sprim(nua,nub) = prod*pioverp*prefactor
                  sprim(nub,nua) = sprim(nua,nub)
                end if
              end do
            end do
          end do
        end do
      end do
    end do
   
  end subroutine

  subroutine etijcalc (m,lamx,la,lb,ce,a,b,ax,bx)
 
    use mod_io, only: faterr, ferror

    implicit none
    integer(kind=ip) :: la,lb,lab,i,j,t,i1,j1,t1,m,lamx
    real(kind=rp) :: ce(-1:2*lamx,-1:lamx,-1:lamx,3)
    real(kind=rp) :: a,b,ax,bx,p,ab,pa,pb,tp

    if (la.lt.0) call ferror ('etijcalc', 'fatal error, la < 0', faterr)
    if (lb.lt.0) call ferror ('etijcalc', 'fatal error, lb < 0', faterr)
    if (la.gt.lamx) call ferror ('etijcalc', 'fatal error, la > lamx', faterr)
    if (lb.gt.lamx) call ferror ('etijcalc', 'fatal error, lb > lamx', faterr)
  
    lab = la + lb
    ce(-1:lab,-1:la,-1:lb,m) = 0.0_rp
    ce(0,0,0,m) = 1.0_rp
    if (lab.eq.0) return
    p  = a + b
    ab = ax - bx
    pa = -b*ab/p
    pb = +a*ab/p
    tp = 1.0_rp/(2.0_rp*p)
    do i = 0,la
      i1 = i-1
      do j = 0,lb
        j1 = j-1
        do t = 1,i+j
          t1 = t-1
          ce(t,i,j,m) = tp*(i*ce(t1,i1,j,m) + j*ce(t1,i,j1,m))/real(t,rp)
        end do
        if (i.lt.la) ce(0,i+1,j,m) = pa*ce(0,i,j,m) + ce(1,i,j,m)
        if (j.lt.lb) ce(0,i,j+1,m) = pb*ce(0,i,j,m) + ce(1,i,j,m)
      end do
    end do 
 
  end subroutine

  !> Destroy pointers
  subroutine end_cint ()
    use mod_memory, only: free
    implicit none
    call free ('mod_cint', 'CINTatm', CINTatm)
    call free ('mod_cint', 'CINTbas', CINTbas)
    call free ('mod_cint', 'CINTenv', CINTenv)
    call free ('mod_cint', 'CINTcoef', CINTcoef)
    call CINTdel_optimizer (CINTopt)
  end subroutine end_cint

  ! Temporal way, efficent should use stack memory
  subroutine allocate_space_for_cint_buff (di,dj)
    use mod_io, only: ferror, faterr
    implicit none
    integer(kind=ip), intent(in) :: di, dj
    integer(kind=ip) :: ier
    allocate (CINTbuff(di,dj),stat=ier)
    if (ier.ne.0) call ferror ("mod_cint", "error in alloc CINTbuff", faterr)
    CINTbuff = 0.0_rp
  end subroutine allocate_space_for_cint_buff
 
  subroutine deallocate_space_for_cint_buff ()
    use mod_io, only: ferror, faterr
    implicit none
    integer(kind=ip) :: ier
    deallocate (CINTbuff,stat=ier)
    if (ier.ne.0) call ferror ("mod_cint", "error in dealloc CINTbuff", faterr)
  end subroutine deallocate_space_for_cint_buff

end module mod_gto
