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
module mod_gto
  
  use mod_prec, only: ip, rp, i8=>size_t
  implicit none
  private

  ! libcint data
  integer(kind=ip) :: CINTnatm
  integer(kind=ip) :: CINTnbas
  integer(kind=ip), allocatable :: CINTatm(:,:)
  integer(kind=ip), allocatable :: CINTbas(:,:)
  real(kind=rp), allocatable :: CINTenv(:)
  real(kind=rp), allocatable :: CINTcoef(:,:)
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
  INTEGER(KIND=IP), PARAMETER :: PTR_RINV_ORIG = 4
  integer(kind=ip), parameter :: POINT_NUC     = 1

  public :: integwfn

contains

  subroutine integwfn ()

    use mod_memory, only: alloc, free
    use mod_io, only: uout, string
    use mod_wfn, only: nprims, nmo, rdm1, rdm2
    implicit none

    real(kind=rp) :: e1, e2, enuc
    integer(kind=ip) :: p, q, r, s, di, dj, dk, dl, t, u
    integer(kind=ip) :: idx, jdx, kdx, ldx, i, j, pq, rs
    integer(kind=ip), dimension(2) :: CINTshls1
    integer(kind=ip), dimension(4) :: CINTshls2
    real(kind=rp), allocatable, dimension(:,:) :: buf1, x, y, tmp
    real(kind=rp), allocatable, dimension(:,:) :: h1e, h1e_ao
    real(kind=rp), allocatable, dimension(:,:,:,:) :: buf2
    real(kind=rp), allocatable, dimension(:,:,:,:) :: h2e
    real(kind=rp), allocatable, dimension(:,:,:,:) :: h2e_ao
    real(kind=rp), allocatable, dimension(:,:,:,:) :: temp

    integer(kind=ip), external :: CINTcgto_cart

    call init_cint ()

    enuc = nnrep ()
    write (uout,'(1x,a,1x,f16.8)') string('# *** NN repulsion energy'), enuc

    write (uout,'(1x,a)') string('# *** Computing 1-electron integrals')
    call alloc ('integwfn', 'buf1', buf1, 21, 21)
    call alloc ('integwfn', 'h1e_ao', h1e_ao, nprims, nprims)
    call alloc ('integwfn', 'h1e', h1e, nmo, nmo)
    idx = 0
    do p = 1,CINTnbas
      CINTshls1(1) = p-1 
      di = CINTcgto_cart(CINTshls1(1), CINTbas)
      jdx = 0
      do q = 1,p
        CINTshls1(2) = q-1
        dj = CINTcgto_cart(CINTshls1(2), CINTbas)
        call cint1e_kin_cart (buf1(1:di,1:dj), CINTshls1, CINTatm, CINTnatm, CINTbas, CINTnbas, CINTenv)
        h1e_ao(idx+1:idx+di,jdx+1:jdx+dj) = buf1(1:di,1:dj)
        call cint1e_nuc_cart (buf1(1:di,1:dj), CINTshls1, CINTatm, CINTnatm, CINTbas, CINTnbas, CINTenv)
        h1e_ao(idx+1:idx+di,jdx+1:jdx+dj) = h1e_ao(idx+1:idx+di,jdx+1:jdx+dj) + buf1(1:di,1:dj)
        h1e_ao(jdx+1:jdx+dj,idx+1:idx+di) = transpose(h1e_ao(idx+1:idx+di,jdx+1:jdx+dj))
        jdx = jdx + dj
      end do
      idx = idx + di
    end do
    write (uout,'(1x,a)') string('# *** Rotating 1-electron integrals to MO basis')
    h1e = matmul(CINTcoef,matmul(h1e_ao,transpose(CINTcoef)))
    e1 = sum(h1e*rdm1)
    write (uout,'(1x,a,1x,f16.8)') string('# *** 1-electron energy'), e1
    call free ('integwfn', 'buf1', buf1)
    call free ('integwfn', 'h1e_ao', h1e_ao)
    call free ('integwfn', 'h1e', h1e)

    write (uout,'(1x,a)') string('# *** Computing 2-electron integrals')
    call alloc ('integwfn', 'buf2', buf2, 21, 21, 21, 21)
    call alloc ('integwfn', 'h2e_ao', h2e_ao, nprims, nprims, nprims, nprims)
    call alloc ('integwfn', 'h2e', h2e, nmo, nmo, nmo, nmo)
    call alloc ('integwfn', 'temp', temp, nprims, nprims, nmo, nmo)
    call alloc ('integwfn', 'x', x, nprims, nprims)
    call alloc ('integwfn', 'y', y, nmo, nmo)
    call alloc ('integwfn', 'tmp', tmp, nprims, nprims)
    idx = 0
    pq = 0
    !TODO: only four fold symmetry, add 8 eight fold symmetry
    do p = 1,CINTnbas
      CINTshls2(1) = p-1 
      di = CINTcgto_cart(CINTshls2(1), CINTbas)
      jdx = 0
      do q = 1,p
        pq = pq + 1
        CINTshls2(2) = q-1 
        dj = CINTcgto_cart(CINTshls2(2), CINTbas)
        kdx = 0
        rs = 0
        do r = 1,CINTnbas
          CINTshls2(3) = r-1 
          dk = CINTcgto_cart(CINTshls2(3), CINTbas)
          ldx = 0
          do s = 1,r
            rs = rs + 1
            CINTshls2(4) = s-1 
            dl = CINTcgto_cart(CINTshls2(4), CINTbas)
            call cint2e_cart(buf2(1:di,1:dj,1:dk,1:dl), CINTshls2, CINTatm, CINTnatm, CINTbas, CINTnbas, CINTenv, CINTopt)
            h2e_ao(idx+1:idx+di,jdx+1:jdx+dj,kdx+1:kdx+dk,ldx+1:ldx+dl) = buf2(1:di,1:dj,1:dk,1:dl)
            do t = kdx+1,kdx+dk
              do u = ldx+1,ldx+dl
                tmp(1:di,1:dj) = h2e_ao(idx+1:idx+di,jdx+1:jdx+dj,t,u)
                h2e_ao(jdx+1:jdx+dj,idx+1:idx+di,t,u) = transpose(tmp(1:di,1:dj))
              end do
            end do
            do t = idx+1,idx+di
              do u = jdx+1,jdx+dj
                tmp(1:dk,1:dl) = h2e_ao(t,u,kdx+1:kdx+dk,ldx+1:ldx+dl)
                h2e_ao(t,u,ldx+1:ldx+dl,kdx+1:kdx+dk) = transpose(tmp(1:dk,1:dl))
                h2e_ao(u,t,ldx+1:ldx+dl,kdx+1:kdx+dk) = h2e_ao(t,u,ldx+1:ldx+dl,kdx+1:kdx+dk) 
              end do
            end do
            ldx = ldx + dl
          end do
          kdx = kdx + dk
        end do
        jdx = jdx + dj
      end do
      idx = idx + di
    end do
    write (uout,'(1x,a)') string('# *** Rotating 2-electron integrals to MO basis')
    write (uout,'(1x,a)') string('# *** First half transform')
    do p = 1,nprims 
      do q = 1,p
        x = h2e_ao(p,q,:,:)
        y = matmul(CINTcoef,matmul(x,transpose(CINTcoef))) 
        temp(p,q,:,:) = y
        temp(q,p,:,:) = transpose(y)
      end do
    end do
    write (uout,'(1x,a)') string('# *** Second half transform')
    do i = 1,nmo
      do j = 1,i
        x = temp(:,:,i,j)
        y = matmul(CINTcoef,matmul(x,transpose(CINTcoef))) 
        h2e(:,:,i,j) = y
        h2e(:,:,j,i) = transpose(y)
      end do
    end do
    write (uout,'(1x,a)') string('# *** Rotated 2-electron integrals')
    e2 = sum(h2e*rdm2)
    write (uout,'(1x,a,1x,f16.8)') string('# *** 2-electron energy'), e2*0.5
    write (uout,'(1x,a,1x,f16.8)') string('# *** Total energy'), e2*0.5 + e1 + enuc
    call free ('integwfn', 'h2e', h2e)
    call free ('integwfn', 'h2e_ao', h2e_ao)
    call free ('integwfn', 'buf2', buf2)
    call free ('integwfn', 'temp', temp)
    call free ('integwfn', 'tmp', tmp)
    call free ('integwfn', 'x', x)
    call free ('integwfn', 'y', y)

    call end_cint ()

  end subroutine integwfn

  real(kind=rp) function nnrep ()

    use mod_wfn, only: ncent, charge, rint
    implicit none

    integer(kind=ip) :: i, j

    nnrep = 0.0_rp
    do i = 1,ncent
      do j = i+1,ncent
        nnrep = nnrep + charge(i)*charge(j)/rint(i,j)
      end do
    end do
    return

  end function

  !> Allocate pointers and data
  subroutine init_cint ()
 
    use mod_futils, only: mergesort
    use mod_memory, only: alloc, free
    use mod_io, only: uout, ferror, faterr
    use mod_param, only: pi, debug 
    use mod_wfn, only: ncent, charge, xyz, ityp, oexp, ngroup, nuexp, &
                       nmo, nprims, coef, numshells, nzexp
    implicit none
  
    integer(kind=ip) :: i, jj, iat, bas_id, di, idx, ish, nua, itipa, l
    integer(kind=ip) :: ic, m, k, itip, itypa(nprims)
    real(kind=rp) :: factor, coefa(nmo,nprims)
    integer(kind=ip) :: offset
    data offset/PTR_ENV_START/
    
    integer(kind=ip), external :: CINTcgto_cart
    integer(kind=ip) :: iord(1000)

    CINTnatm = ncent
    CINTnbas = numshells

    call alloc ('mod_cint', 'CINTatm', CINTatm, ATM_SLOTS, CINTnatm)
    call alloc ('mod_cint', 'CINTbas', CINTbas, BAS_SLOTS, CINTnbas)
    call alloc ('mod_cint', 'CINTcoef', CINTcoef, nmo, nprims)
    call alloc ('mod_cint', 'CINTenv', CINTenv, PTR_ENV_START+CINTnatm*4+CINTnbas*2)

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
    write (uout,'(1x,a,i5)') "# +++ Total number of shells : ", CINTnbas

    bas_id = 0
    do iat = 1,ncent
      ! atomic info
      CINTatm(CHARGE_OF,iat) = int(charge(iat),ip) ! atomic nuclear charge
      CINTatm(PTR_COORD,iat) = offset ! pointer to coordinates
      CINTenv(offset+1) = xyz(iat,1) ! x (Bohr) 
      CINTenv(offset+2) = xyz(iat,2) ! y (Bohr) 
      CINTenv(offset+3) = xyz(iat,3) ! z (Bohr) 
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
          factor = 1.0_rp
          factor = sqrt(4.0_rp*pi)  !sp factor
        else if (itipa.eq.2) then
          l = 1
          factor = 1.0_rp
          factor = sqrt((4.0_rp*pi)/3) !sp factor
        else if (itipa.eq.5) then
          l = 2
          factor = 1.0_rp
        else if (itipa.eq.11) then
          l = 3
          factor = 1.0_rp
        else if (itipa.eq.21) then
          l = 4
          factor = 1.0_rp
        else if (itipa.eq.36 .or. itipa.eq.56) then
          l = 5
          factor = 1.0_rp
        else 
          call ferror ('mod_cint' , 'not supported i or higher primitives', faterr)
        end if
        CINTbas(ANG_OF,bas_id) = l ! angular momentum
        CINTbas(NPRIM_OF,bas_id) = 1 ! primitives
        CINTbas(NCTR_OF,bas_id) = 1 ! contradted
        CINTbas(KAPPA_OF,bas_id) = 0 ! kappa for spinor basis
        ! exponents
        CINTbas(PTR_EXP,bas_id) = offset ! pointer to exp of primitive GTO
        CINTenv(offset+1) = oexp(nua) ! orbital exponent
        offset = offset + 1 ! advance pointer due to exponent
        ! coefficients
        CINTbas(PTR_COEFF,bas_id) = offset ! pointer to contraction coefficients
        CINTenv(offset+1) = 1.0_rp*factor !factor due to sp contraction
        offset = offset + 1 ! advance pointer due to primitives
        ! print basis info
        if (debug) write (uout,'(1x,8i8)') CINTbas(:,bas_id)
      end do
    end do

    ! Finish storing info
    if (debug) then
      write (uout,'(1x,a)') "# CINTEnv info "
      write (uout,'(19(1x,f12.6))') (CINTenv(i),i=1,offset)
    end if
    CINTopt = 0
    call cint2e_cart_optimizer (CINTopt, CINTatm, CINTnatm, CINTbas, CINTnbas, CINTenv)

    if (debug) then
      do i = 1,CINTnbas
        di = CINTcgto_cart(i-1, CINTbas)
        write (uout,'(1x,a,1x,i0,1x,a,1x,i0)') "# CINT Primitives in shell", i, ":", di
      end do
    end if
    write (uout,'(1x,a)') "# *** End initialization of libcint data"

  end subroutine init_cint

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

end module mod_gto
