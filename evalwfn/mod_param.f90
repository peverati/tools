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
module mod_param

  use mod_prec, only: rp, ip
  implicit none
  public

  ! global parameters
  real(kind=rp), parameter :: pi = 3.14159265358979323846264338328_rp
  real(kind=rp), parameter :: bohrtoa = 0.52917720859d0 !< bohr to angstrom conversion factor

  real(kind=rp), parameter :: vsmall = 1d-80 !< a very small number

  ! some options
  logical :: verbose, debug

  logical :: isdata

contains 

  subroutine init_param ()

    implicit none
    integer(kind=ip) :: i, n, clock
    integer(kind=ip), dimension(:), allocatable :: seed

    isdata = .false.
    verbose = .false.
    debug = .false.

    ! random seed
    call random_seed(size=n)
    allocate(seed(n))
    call system_clock(count=clock)
    seed = clock + 37*(/(i-1, i =1,n)/)
    call random_seed(put=seed)
    deallocate(seed)

  end subroutine init_param

end module mod_param
