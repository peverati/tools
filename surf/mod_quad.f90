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
module mod_quad

  use mod_prec, only: rp, ip
  implicit none
  public
  
contains

  ! Abscissas and weights for the theta quadrature. 
  subroutine weightheta (iq,x,w,n)

    use mod_io, only: ferror, faterr
    implicit none

    ! Arguments
    integer(kind=ip), intent(in) :: iq, n
    real(kind=rp), intent(out) :: x(n), w(n)

    if (iq.eq.1) then ! Gauss-Legendre quadrature 
      call gauleg (-1.0_rp,+1.0_rp,x,w,n)
    else if (iq.eq.2) then ! Clenshaw-Curtis quadrature
      call genclcu (-1.0_rp,+1.0_rp,x,w,n)
    else if (iq.eq.3) then ! Gauss-Chebychev of first kind 
      call gaucheb1 (-1.0_rp,+1.0_rp,x,w,n)
    else if (iq.eq.4) then ! Gauss-Chebychev of second kind 
      call gaucheb2 (x,w,n)
    else if (iq.eq.5) then ! Perez-Jorda (Gauss-Chebychev, 2nd kind) quadrature. 
      call pjt (x,w,n)
    else
      call ferror ('mod_quad', 'bad theta quadrature', faterr)
    end if
 
  end subroutine

  ! Performs three successive rotations R_z, R_y, R_x
  ! about the Z, Y, and X axes, respectively, of NANG points placed on 
  ! an unit radius sphere, and defined by the polar coordinates 
  ! cos(theta)=ct, sin(theta)=st, cos(phi)=cp, and sin(phi)=sp.
  subroutine rotagrid (ct,st,cp,sp,npang,angx,angy,angz)

    implicit none
    integer(kind=ip), intent(in) :: npang
    real(kind=rp), intent(inout) :: ct(npang), st(npang), cp(npang), sp(npang)
    real(kind=rp), intent(in) :: angx, angy, angz
    real(kind=rp), parameter :: eps = 1.0d-07
    real(kind=rp), parameter :: zero = 0.0_rp
    real(kind=rp), parameter :: one = 1.0_rp

    real(kind=rp) :: xmat(3,3), ymat(3,3), zmat(3,3), rot(3,3)
    real(kind=rp) :: prod, x, y, z, xp, yp, zp, r, rxy
    real(kind=rp) :: cangx, cangy, cangz, sangx, sangy, sangz
    integer(kind=ip) :: i, j, k

    ! Rotacion matrix about X axis.
    xmat = zero
    xmat(1,1) = one
    cangx = cos(angx)
    sangx = sin(angx)
    xmat(2,2) = +cangx
    xmat(3,2) = +sangx
    xmat(2,3) = -sangx
    xmat(3,3) = +cangx

    ! Rotacion matrix about Y axis.
    ymat = zero
    ymat(2,2) = one
    cangy = cos(angy)
    sangy = sin(angy)
    ymat(1,1) = +cangy
    ymat(3,1) = +sangy
    ymat(1,3) = -sangy
    ymat(3,3) = +cangy

    ! Rotacion matrix about Z axis.
    zmat = zero
    zmat(3,3) = one
    cangz = cos(angz)
    sangz = sin(angz)
    zmat(1,1) = +cangz
    zmat(2,1) = +sangz
    zmat(1,2) = -sangz
    zmat(2,2) = +cangz

    ! Full rotacion matrix. R = R_X * R_Y * R_Z
    ! TODO: change to matmul
    do i = 1,3
      do j = 1,3
        prod = zero
        do k = 1,3
          prod = prod + ymat(i,k)*zmat(k,j)
        end do
        rot(i,j) = prod
      end do
    end do
    zmat(1:3,1:3) = rot(1:3,1:3)
    do i = 1,3
      do j = 1,3
        prod = zero
        do k = 1,3
          prod = prod + xmat(i,k)*zmat(k,j)
        end do
        rot(i,j) = prod
      end do
    end do

    ! Rotate angles.
    do i = 1,npang
      x = st(i)*cp(i)
      y = st(i)*sp(i)
      z = ct(i)
      xp = rot(1,1)*x + rot(1,2)*y + rot(1,3)*z
      yp = rot(2,1)*x + rot(2,2)*y + rot(2,3)*z
      zp = rot(3,1)*x + rot(3,2)*y + rot(3,3)*z
      ! Rebuild CT,ST,CP, and SP
      rxy = xp*xp + yp*yp
      r = sqrt(rxy+zp*zp)
      if (rxy.lt.eps) then
        if (zp.ge.zero) then
          ct(i) = +one
        else
          ct(i) = -one
        endif
        st(i) = zero
        sp(i) = zero
        cp(i) = one
      else
        rxy = sqrt(rxy)
        ct(i) = zp/r
        st(i) = sqrt((one-ct(i))*(one+ct(i)))
        cp(i) = xp/rxy
        sp(i) = yp/rxy
      end if
    end do
            
  end subroutine

  subroutine gauleg (x1,x2,x,w,n)

    use mod_param, only: pi
    implicit none

    ! Parameters
    real(kind=rp), parameter :: eps = 3.0d-16

    ! Arguments
    integer(kind=ip), intent(in) :: n
    real(kind=rp), intent(in) :: x1, x2
    real(kind=rp), intent(out) :: x(n), w(n)

    ! Local vars
    real(kind=rp) :: p1, p2, p3, pp, xl, xm, z, z1
    integer(kind=ip) :: i, j, m

    ! The roots are symmetric in the interval, so we only have to find
    ! half of them.
    m = (n+1_ip)/2_ip
    xm = 0.5_rp*(x2+x1)
    xl = 0.5_rp*(x2-x1)

    ! Loop over the desired roots.
    do i = 1,m
      z = cos(pi*(i-0.25_rp)/(n+0.5_rp))
1     continue
        p1 = 1.0_rp
        p2 = 0.0_rp
        ! Loop up the recurrence relation to get the Legendre polynomial 
        ! evaluated at z.
        do j = 1,n
          p3 = p2
          p2 = p1
          p1 = ((2.0_rp*j-1.0_rp)*z*p2-(j-1.0_rp)*p3)/j
        end do
        ! p1 is now the desired Legendre polynomial. We next compute pp,
        ! derivative , by a standard relation involving p2, the polyn-
        ! omial of one lower order.
        pp = n*(z*p1-p2)/(z*z-1.0_rp)
        z1 = z
        ! Newton method.
        z = z1 - p1/pp
      if (abs(z-z1).gt.eps) go to 1
      ! Scale the root to the desired interval.
      x(i) = xm - xl*z
      ! and put in its symmetric counterpart.
      x(n+1-i) = xm + xl*z
      ! compute the weight.
      w(i) = 2.0_rp*xl/((1.0_rp-z*z)*pp*pp)
      ! and its symmetric counterpart.
      w(n+1-i)=w(i)
    end do
      
  end subroutine

  subroutine genclcu (a,b,x,w,n)

    use mod_param, only: pi
    implicit none
    real(kind=rp), parameter :: half = 0.5_rp
    real(kind=rp), parameter :: one = 1.0_rp
    real(kind=rp), parameter :: two = 2.0_rp

    ! Arguments
    integer(kind=ip), intent(in) :: n
    real(kind=rp), intent(in) :: a, b
    real(kind=rp), intent(out) :: x(n), w(n) 

    ! Local vars
    integer(kind=ip) :: i, j, nn, nn11, nn12, nnh
    real(kind=rp) :: z, factor, piovernn, term, x1, x2
 
    x1 = (b-a)*half
    x2 = (b+a)*half
    nn = n - 1_ip
    nn11 = nn*nn
    nn12 = nn*(nn-1_ip)
    nnh = (nn-1_ip)/2_ip
    piovernn = pi/real(nn,rp)
    factor = two*(b-a)/real(nn,rp)
    x(1) = b
    x(n) = a
    do i = 1,nn-1
      z = piovernn*i
      x(i+1) = x1*cos(z) + x2
      w(i+1) = 0.0_rp
      do j = 0,nnh
        if (j.eq.0_ip) then
          term = half
        else
          term = one/(1.0_rp-4.0_rp*j*j)*cos(2.0_rp*j*z)
        end if
        w(i+1) = w(i+1) + term
      end do
      w(i+1) = w(i+1)*factor
    end do
    if (mod(n,2).eq.0) then
      w(1) = x1/nn11
      w(n) = w(1)
    else
      w(1) = x1/nn12
      w(n) = w(1)
    end if
 
  end subroutine

  subroutine gaucheb1 (a,b,x,w,n)

    use mod_param, only: pi
    implicit none
    real(kind=rp), parameter :: half = 0.5_rp
 
    integer(kind=ip) :: n, i
    real(kind=rp) :: x(n),w(n)
    real(kind=rp) :: a,b,z,x1,x2
 
    x1 = (b-a)*half
    x2 = (b+a)*half
    do i = 1,n
      z = pi*(real(i,rp)-half)/real(n,rp)
      x(i) = x1*cos(z)+x2
      w(i) = (pi/real(n,rp))*x1*sin(z)
    end do
 
  end subroutine
 
  subroutine gaucheb2 (x,w,n)
 
    use mod_param, only: pi
    implicit none

    integer(kind=ip) :: n, i
    real(kind=rp) :: x(n),w(n)
    real(kind=rp) :: xad, xa

    xad = pi/real(n+1,rp)
    do i = 1,n
      xa = i*xad
      w(i) = xad*sin(xa)
      x(i) = cos(xa)
    end do
 
  end subroutine

  subroutine pjt (x,w,n)

    use mod_param, only: pi
    implicit none
    real(kind=rp), parameter :: twothird = 2.0_rp/3.0_rp
 
    integer(kind=ip) :: n, i
    real(kind=rp) :: x(n),w(n)
    real(kind=rp) :: z,st,ct,dn
 
    do i = 1,n
      dn = real(n+1_ip,rp)
      z = i*pi/dn
      st = sin(z)
      ct = cos(z)
      w(i) = 16.0_rp/3.0_rp/dn*st*st*st*st
      x(i) = 1.0_rp-real(2*i,rp)/dn+2.0_rp/pi*(1.0_rp+twothird*st*st)*ct*st
    end do
 
  end subroutine

end module mod_quad
