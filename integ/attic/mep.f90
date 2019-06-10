
  subroutine mepwfn ()

    use mod_memory, only: alloc, free
    use mod_io, only: uout, string, warning, ferror
    use mod_wfn, only: nprims, nmo, rdm1
    implicit none

    real(kind=rp) :: e1
    integer(kind=ip), external :: CINTcgto_cart
    real(kind=rp), dimension(nprims, nprims) :: mep_ao
    real(kind=rp), dimension(nmo, nmo) :: mep_mo
    integer(kind=ip) :: p, q, di, dj, idx, jdx

    mep_ao = 0.0_rp
    e1 = 0.0_rp
    CINTxpoint = 0.0_rp
    write (uout,'(1x,a,1x,3(f12.8,1x))') string('# *** Computing MEP at point'), CINTxpoint
    call init_cint ()
    CINTenv(PTR_RINV_ORIG+1:PTR_RINV_ORIG+3) = CINTxpoint
    idx = 0
    do p = 1,CINTnbas
      CINTshls1(1) = p-1 
      di = CINTcgto_cart(CINTshls1(1), CINTbas)
      jdx = 0
      do q = 1,p
        CINTshls1(2) = q-1
        dj = CINTcgto_cart(CINTshls1(2), CINTbas)
        call allocate_space_for_cint_buff1 (di,dj)
        call cint1e_rinv_cart (CINTbuff1, CINTshls1, CINTatm, CINTnatm, CINTbas, CINTnbas, CINTenv)
        mep_ao(idx+1:idx+di,jdx+1:jdx+dj) = CINTbuff1(:,:)
        mep_ao(jdx+1:jdx+dj,idx+1:idx+di) = transpose(mep_ao(idx+1:idx+di,jdx+1:jdx+dj))
        call deallocate_space_for_cint_buff1 ()
        jdx = jdx + dj
      end do
      idx = idx + di
    end do
    mep_mo = matmul(CINTcoef,matmul(mep_ao,transpose(CINTcoef)))
    e1 = sum(mep_mo*rdm1)
    write (uout,'(1x,a,1x,f16.8)') string('# *** MEP value:'), e1

  end subroutine

