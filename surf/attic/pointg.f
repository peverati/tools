      subroutine pointg (p,rho,gradmod)
c
      real(kind=dp), intent(in) :: p(3)
      integer(kind=ip) :: it(3), k, jj, i, itip, n
      real(kind=dp) :: dis2, ori, dp2, aexp, x, x2, dp2x2, dp2x
      real(kind=dp) :: f12, f13, f23, g21, g23, g32, cfj, tmp(3), vec(3)
      real(kind=dp) :: fun(3), fun1(3), fun2(3), xcoor(3), ggg(3)
      real(kind=dp), dimension(nmo) :: gun
      real(kind=dp), dimension(nmo,3) :: gun1
      real(kind=dp), dimension(nmo,6) :: gun2
c
      rho  = 0.0_dp
      grad = 0.0_dp
      king = 0.0_dp
      gun  = 0.0_dp
      gun1 = 0.0_dp
      gun2 = 0.0_dp
      grad = 0.0_dp
      hess = 0.0_dp
      ggg  = 0.0_dp
      tmp  = 0.0_dp
c
c     Run over centers
c
      do ic=1,ncent
        xcoor(:)=p(:)-xyz(ic,:)
        dis2=dot_product(xcoor,xcoor)
        do m=1,ngroup(ic)
          k=nuexp(ic,m,1)
          if (dis2.gt.rcutte(k)) cycle
          ori=-oexp(k)
          dp2=ori+ori
          aexp=exp(ori*dis2)
          do jj=1,nzexp(ic,m)
            i=nuexp(ic,m,jj)
            itip=ityp(i)
            it(:)=nlm(itip,:)
            do j=1,3
              n=it(j)
              x=xcoor(j)
              x2=x*x
              dp2x2=dp2*x2
              if (n.eq.0) then
                dp2x=dp2*x
                fun2(j)=dp2*(1d0+dp2x2)
                fun1(j)=dp2x
                fun(j)=1d0
              elseif (n.eq.1) then
                fun2(j)=dp2*x*(3d0+dp2x2)
                fun1(j)=1d0+dp2x2
                fun(j)=x
              elseif (n.eq.2) then
                fun2(j)=2d0+dp2x2*(5d0+dp2x2)
                fun1(j)=x*(2d0+dp2x2)
                fun(j)=x2
              elseif (n.eq.3) then
                fun2(j)=x*(6d0+dp2x2*(7d0+dp2x2))
                fun1(j)=x2*(3d0+dp2x2)
                fun(j)=x*x2
              elseif (n.eq.4) then
                fun2(j)=x2*(x2*(dp2*(dp2x2+9d0))+12d0)
                fun1(j)=x2*x*(4d0+dp2x2)
                fun(j)=x2*x2
              elseif (n.eq.5) then
                fun2(j)=x2*x*(x2*(dp2*(dp2x2+11d0))+20d0) 
                fun1(j)=x2*x2*(5d0+dp2x2)
                fun(j)=x2*x2*x
              endif
            enddo
*
            f23=fun(2)*fun(3)*aexp
            f13=fun(1)*fun(3)*aexp
            f12=fun(1)*fun(2)*aexp
            g23=fun1(2)*fun(3)*aexp
            g32=fun1(3)*fun(2)*aexp
            g21=fun1(2)*fun(1)*aexp
c                            
            do j=1,nmo
              cfj=coef(j,i)
              gun(j)=gun(j)+cfj*(fun(1)*f23)
c
              gun1(j,1)=gun1(j,1)+cfj*(fun1(1)*f23)
              gun1(j,2)=gun1(j,2)+cfj*(fun1(2)*f13)
              gun1(j,3)=gun1(j,3)+cfj*(fun1(3)*f12)
c
              gun2(j,1)=gun2(j,1)+cfj*(fun2(1)*f23)
              gun2(j,3)=gun2(j,3)+cfj*(fun2(2)*f13)
              gun2(j,6)=gun2(j,6)+cfj*(fun2(3)*f12)
c
              gun2(j,2)=gun2(j,2)+cfj*(fun1(1)*g23)
              gun2(j,4)=gun2(j,4)+cfj*(fun1(1)*g32)
              gun2(j,5)=gun2(j,5)+cfj*(fun1(3)*g21)
            enddo
          enddo
        enddo
      enddo
c
c     run again over orbitals
c
      do i=1,nmo
        fac=occ(i)
c
        tmp(1)=tmp(1)+fac*gun1(i,1)
        tmp(2)=tmp(2)+fac*gun1(i,2)
        tmp(3)=tmp(3)+fac*gun1(i,3)
c
        hess(1,1)=hess(1,1)+gun2(i,1)
        hess(2,2)=hess(2,2)+gun2(i,3)
        hess(3,3)=hess(3,3)+gun2(i,6)
        hess(1,2)=hess(1,2)+gun2(i,2)
        hess(1,3)=hess(1,3)+gun2(i,4)
        hess(2,3)=hess(2,3)+gun2(i,5)
c
        king=king+fac*(gun1(i,1)*gun1(i,1)+
     &                 gun1(i,2)*gun1(i,2)+
     &                 gun1(i,3)*gun1(i,3))
      enddo
c
      king=king*0.5_dp
c
      rho=king
      vec(1)=hess(1,1)
      vec(2)=hess(1,2)
      vec(3)=hess(1,3)
      grad(1)=dot_product(tmp,vec)
      vec(1)=hess(1,2)
      vec(2)=hess(2,2)
      vec(3)=hess(2,3)
      grad(2)=dot_product(tmp,vec)
      vec(1)=hess(1,3)
      vec(2)=hess(2,3)
      vec(3)=hess(3,3)
      grad(3)=dot_product(tmp,vec)
c
      end subroutine pointg
