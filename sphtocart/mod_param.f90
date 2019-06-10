! Copyright (c) 2018 
! Jose Luis Casals Sainz <jluiscasalssainz@gmail.com>,
! Ángel Martín Pendás <angel@fluor.quimica.uniovi.es> and 
! Evelio Francisco Miguelez <evelio@uniovi.es>. 
!
! sph2real is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or (at
! your option) any later version.
! 
! sph2real is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
module mod_param

  implicit none
  public

  integer(kind=4), parameter :: mxfact = 100 
  integer(kind=4), parameter :: mcomb = 12
  
  real(kind=8) :: fact(0:mxfact)
  real(kind=8) :: comb(0:mcomb,0:mcomb) 

contains

  subroutine init_param ()

    implicit none

    integer(kind=4) :: i, j, ij
   
    fact(0) = 1.0d0
    do i = 1,mxfact
      fact(i) = fact(i-1)*real(i,8)
    end do

    do i = 0,mcomb
      do j = 0,i/2
        ij = i-j
        comb(i,j) = fact(i)/fact(j)/fact(ij)
        if (j.ne.ij) comb(i,ij) = comb(i,j)
      end do
    end do

  end subroutine

end module
