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
module mod_param

  implicit none
  public

  real(kind=8), parameter :: pi = 3.14159265358979324D+00
  integer(kind=4), parameter :: mxfact = 100

  real(kind=8) :: facd(-mxfact:mxfact) 

contains

  subroutine init_param()
  
    implicit none
    integer(kind=4) :: mf21, i

    facd(-1) = 1d0
    facd(0) = 1d0
    do i = 1,mxfact
      facd(i) = real(i,8)*facd(i-2)
    end do
    mf21 = mxfact/2 - 1
    do i = 1,mf21
      facd(-i-i-1) = (-1)**i/facd(i+i-1)
    end do

  end subroutine init_param

end module mod_param
