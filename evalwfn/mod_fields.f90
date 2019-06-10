module mod_fields

  use mod_prec, only: rp, ip
  implicit none
  private

  public :: density, denmos

contains
  
  subroutine density (p,rho)

    use mod_wfn, only: xyz=>coords, ncent, nlm, ityp, occ, nmo, oexp, coeff, &
                       ngroup, nuexp, nzexp, rcutte
    implicit none
 
    real(kind=rp), intent(in) :: p(3)
    real(kind=rp), intent(out) :: rho

    integer(kind=ip) :: it(3), i, ic, itip, j, jj, k, m, n
    real(kind=rp) :: aexp, cfj, dis2, x2, f123, ori, x
    real(kind=rp) :: xcoor(3), fun(3), dp2
    real(kind=rp), dimension(nmo) :: gun
 
    rho = 0.0_rp 
    gun = 0.0_rp

    do ic = 1,ncent
      xcoor(:) = p(:)-xyz(ic,:)
      dis2 = xcoor(1)*xcoor(1)+xcoor(2)*xcoor(2)+xcoor(3)*xcoor(3)
      do m = 1,ngroup(ic)
        k = nuexp(ic,m,1)
        if (dis2.gt.rcutte(ic,m)) goto 2
        ori = -oexp(k)
        dp2 = ori+ori
        aexp = exp(ori*dis2)
        do jj = 1,nzexp(ic,m)
          i = nuexp(ic,m,jj)
          itip = ityp(i)
          it(:) = nlm(itip,:)
          do j=1,3
            n = it(j)
            x = xcoor(j)
            if (n.eq.0) then
              fun(j) = 1.0_rp
            else if (n.eq.1) then
              fun(j)=x
            else if (n.eq.2) then
              x2 = x*x
              fun(j) = x2
            else if (n.eq.3) then
              x2 = x*x
              fun(j) = x*x2
            else if (n.eq.4) then
              x2 = x*x
              fun(j) = x2*x2
            else if (n.eq.5) then
              x2 = x*x
              fun(j) = x2*x2*x
            end if
          end do
          f123 = fun(1)*fun(2)*fun(3)*aexp
          do j= 1,nmo
            cfj = coeff(j,i)
            gun(j) = gun(j) + cfj*f123
          end do
        end do
 2      continue
      end do
    end do
 
    rho = dot_product(occ,gun*gun)

  end subroutine

  subroutine denmos (p,dmos)

    use mod_wfn, only: xyz=>coords, ncent, nlm, ityp, occ, nmo, oexp, coeff, &
                       ngroup, nuexp, nzexp, rcutte
    implicit none
 
    real(kind=rp), intent(in) :: p(3)
    real(kind=rp), intent(out), dimension(nmo) :: dmos

    integer(kind=ip) :: it(3), i, ic, itip, j, jj, k, m, n
    real(kind=rp) :: aexp, cfj, dis2, x2, f123, ori, x
    real(kind=rp) :: xcoor(3), fun(3), dp2
    real(kind=rp), dimension(nmo) :: gun
 
    gun = 0.0_rp
    dmos = 0.0_rp

    do ic = 1,ncent
      xcoor(:) = p(:)-xyz(ic,:)
      dis2 = xcoor(1)*xcoor(1)+xcoor(2)*xcoor(2)+xcoor(3)*xcoor(3)
      do m = 1,ngroup(ic)
        k = nuexp(ic,m,1)
        if (dis2.gt.rcutte(ic,m)) goto 2
        ori = -oexp(k)
        dp2 = ori+ori
        aexp = exp(ori*dis2)
        do jj = 1,nzexp(ic,m)
          i = nuexp(ic,m,jj)
          itip = ityp(i)
          it(:) = nlm(itip,:)
          do j=1,3
            n = it(j)
            x = xcoor(j)
            if (n.eq.0) then
              fun(j) = 1.0_rp
            else if (n.eq.1) then
              fun(j)=x
            else if (n.eq.2) then
              x2 = x*x
              fun(j) = x2
            else if (n.eq.3) then
              x2 = x*x
              fun(j) = x*x2
            else if (n.eq.4) then
              x2 = x*x
              fun(j) = x2*x2
            else if (n.eq.5) then
              x2 = x*x
              fun(j) = x2*x2*x
            end if
          end do
          f123 = fun(1)*fun(2)*fun(3)*aexp
          do j= 1,nmo
            cfj = coeff(j,i)
            gun(j) = gun(j) + cfj*f123
          end do
        end do
 2      continue
      end do
    end do
 
    do i = 1,nmo
      dmos(i) = occ(i)*gun(i)*gun(i)
    end do

  end subroutine

end module mod_fields
