! Copyright (c) 2018 
! Jose Luis Casals Sainz <jluiscasalssainz@gmail.com>,
! Ángel Martín Pendás <angel@fluor.quimica.uniovi.es> and 
! Evelio Francisco Miguelez <evelio@uniovi.es>. 
!
! connect is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or (at
! your option) any later version.
! 
! connect is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
module mod_math

  implicit none
  public

  real(kind=8), parameter :: pi = 3.14159265358979323846264338328d0
  real(kind=8), parameter :: rad2deg = 180.0d0/pi

contains

  real(kind=8) function dist (xyza, xyzb)

    implicit none
    real(kind=8), intent(in), dimension(3) :: xyza, xyzb

    real(kind=8) :: a, b, c

    a = (xyzb(1) - xyza(1))*(xyzb(1) - xyza(1)) 
    b = (xyzb(2) - xyza(2))*(xyzb(2) - xyza(2)) 
    c = (xyzb(3) - xyza(3))*(xyzb(3) - xyza(3)) 

    dist = sqrt(a+b+c)
    return

  end function dist

  ! calculate unit vector between to 3-d cartesian coordinates
  subroutine unitv (xyza, xyzb, xyzc)

    implicit none
    real(kind=8), intent(in), dimension(3) :: xyza, xyzb
    real(kind=8), intent(out), dimension(3) :: xyzc

    real(kind=8) :: d
    real(kind=8) :: a, b, c

    d = dist(xyza,xyzb)
    a = (xyzb(1) - xyza(1)) 
    b = (xyzb(2) - xyza(2)) 
    c = (xyzb(3) - xyza(3)) 
    xyzc(1) = a/d
    xyzc(2) = b/d
    xyzc(3) = c/d
  
  end subroutine unitv 

  ! calculate dot product between two unit vectors
  real(kind=8) function unitdp (xyza, xyzb)

    implicit none
    real(kind=8), intent(in), dimension(3) :: xyza, xyzb
    
    real(kind=8) :: a, b, c

    unitdp = 0.0d0
    a = xyzb(1)*xyza(1)
    b = xyzb(2)*xyza(2)
    c = xyzb(3)*xyza(3)
    unitdp = a + b + c
    unitdp = max(min(unitdp, 1.0d0), -1.0d0)

    return
  
  end function unitdp 

  ! calculate unit cross product between two unit vectors
  subroutine unitcp (xyza, xyzb, xyzc)

    implicit none
    real(kind=8), intent(in), dimension(3) :: xyza, xyzb
    real(kind=8), intent(out), dimension(3) :: xyzc

    real(kind=8) :: cos_12, sin_12

    cos_12 = unitdp(xyza, xyzb)
    sin_12 = sqrt(1d0 - cos_12*cos_12)
    xyzc(1) = (xyza(2)*xyzb(3) - xyza(3)*xyzb(2))/sin_12
    xyzc(2) = (xyza(3)*xyzb(1) - xyza(1)*xyzb(3))/sin_12
    xyzc(3) = (xyza(1)*xyzb(2) - xyza(2)*xyzb(1))/sin_12

  end subroutine unitcp

  ! calculate angle between three 3-d cartesian coordinates
  real(kind=8) function a123 (xyza, xyzb, xyzc)

    implicit none
    real(kind=8), intent(in), dimension(3) :: xyza, xyzb, xyzc

    real(kind=8), dimension(3) :: u21, u23
    real(kind=8) :: dp2123

    call unitv(xyzb, xyza, u21)
    call unitv(xyzb, xyzc, u23)
    dp2123 = unitdp(u21, u23)
    a123 = rad2deg*acos(dp2123)
    return 

  end function a123

end module mod_math
