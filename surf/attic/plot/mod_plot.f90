! Copyright (c) 2018 
! Jose Luis Casals Sainz <jluiscasalssainz@gmail.com>,
! Ángel Martín Pendás <angel@fluor.quimica.uniovi.es> and 
! Evelio Francisco Miguelez <evelio@uniovi.es>. 
!
! aom is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or (at
! your option) any later version.
! 
! aom is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
module mod_plot

  use mod_prec, only: rp, ip
  implicit none
  private

  integer(kind=ip), parameter :: minter = 10 

  real(kind=rp) :: xnuc(3)
  integer(kind=ip) :: npang
  real(kind=rp) :: rmin, rmax
  real(kind=rp), allocatable, dimension(:) :: ct, cp, st, sp, angw 
  real(kind=rp), allocatable, dimension(:,:) :: rlimsurf
  integer(kind=ip), allocatable, dimension(:) :: nlimsurf

  public :: parse_plot, init_plot

contains

  subroutine parse_plot ()

    use mod_memory, only: free, alloc
    use mod_param, only: pi
    use mod_io, only: ferror, faterr, equal, isword, uin, mline, &
                      string, uout, lgetword, getline, isword
    implicit none

    logical :: ok
    integer(kind=ip) :: lp, lsu, i
    real(kind=rp) :: p(3), xcoor(3), r
    real(kind=rp), allocatable, dimension(:) :: x, y, z
    character(len=:), allocatable :: line, subline, word
    character(len=mline) :: filename

    ! Read the input file
    write (uout,'(1x,a)') string('# *** Reading plot options ***')
    do while (getline(uin,line))
      lp = 1
      word = lgetword(line,lp)
      subline = line(lp:)

      if (equal(word,'#')) then
        continue

      else if (equal(word,'load')) then
        ok = isword(filename, line, lp)
        if (.not.ok) call ferror ('plot', 'wrong load line', faterr) 

      ! End of input
      else if (equal(word,'endplot')) then
        exit

      else
        call ferror ('mod_plot', 'unknown option', faterr)

      end if
    end do
        
    write (uout,'(1x,a,1x,a)') string('# *** Reading data from'), trim(filename)
    call read_surface (filename)
    call alloc ('mod_plot', 'x', x, npang)
    call alloc ('mod_plot', 'y', y, npang)
    call alloc ('mod_plot', 'z', z, npang)
    lsu = 888
    open (lsu,file=trim(filename)//"-ptxt")
    do i = 1,npang
      r = rlimsurf(i,nlimsurf(i)) ! TODO: add multiple intersec
      xcoor(1) = r*st(i)*cp(i)
      xcoor(2) = r*st(i)*sp(i)
      xcoor(3) = r*ct(i)
      p(:) = xnuc(:) + xcoor(:)
      p(:) = p(:)*0.5291772  
      x(i) = p(1)
      y(i) = p(2)
      z(i) = p(3)
      write (lsu,'(3(1x,f12.8))') p
    end do
    close (lsu)
    call free ('mod_plot', 'x', x)
    call free ('mod_plot', 'y', y)
    call free ('mod_plot', 'z', z)
    call deallocate_space_for_surface ()

  end subroutine parse_plot

  subroutine read_surface (filename) 

    use mod_memory, only: alloc, free
    use mod_io, only: fourchar, uout, ferror, faterr, string
    implicit none
    character(len=*), intent(in) :: filename

    integer(kind=ip) :: lsu, ier, idumi, j, nsurf, k

    lsu = 999
    open (lsu,file=trim(filename),form='unformatted',status='old',iostat=ier)
    if (ier == 0) then
      read (lsu) npang, idumi
      write (uout,'(1x,a,1x,i0)') string('# Number of angular points in surface :'), npang
      call allocate_space_for_surface ()
      read (lsu) (nlimsurf(j),j=1,npang)
      do j = 1,npang
        nsurf = nlimsurf(j)
        if (nsurf.gt.minter) then
          call ferror ('mod_plot', 'too many intersections', faterr)
        end if
        read (lsu) ct(j), st(j), cp(j), sp(j), angw(j), (rlimsurf(j,k),k=1,nsurf)
      end do
      read (lsu) rmin, rmax
      write (uout,'(1x,a,2(1x,f8.4))') string('# Limits of ZFS surface :'), rmin, rmax
      read (lsu) xnuc(1), xnuc(2), xnuc(3)
      write (uout,'(1x,a,3(1x,f8.4))') string('# Nuclear coordinates for atom :'), xnuc
      close (lsu)
    else
      call ferror ('mod_plot', 'surface file not available', faterr)
    end if

  end subroutine

  logical function inbasin (r,j)
  
    implicit none
    real(kind=rp), parameter :: eps=1d-6
    integer(kind=ip), intent(in) :: j
    real(kind=rp), intent(in) :: r
 
    integer(kind=ip) :: k
    real(kind=rp) :: rs1, rs2

    inbasin = .false.
    rs1 = 0.0_rp

    do k = 1,nlimsurf(j)
      rs2 = rlimsurf(j,k)
      if (r.ge.rs1-eps .and. r.le.rs2+eps) then
        if (mod(k,2).eq.0) then
          inbasin = .false.
        else
          inbasin = .true.
        end if
        exit
      end if
      rs1 = rs2
    end do
    return

  end function

  subroutine allocate_space_for_surface ()

    use mod_memory, only: alloc
    implicit none

    call alloc ('mod_aom', 'rlimsurf', rlimsurf, npang, minter)
    call alloc ('mod_aom', 'nlimsurf', nlimsurf, npang)
    call alloc ('mod_aom', 'ct', ct, npang)
    call alloc ('mod_aom', 'st', st, npang)
    call alloc ('mod_aom', 'cp', cp, npang)
    call alloc ('mod_aom', 'sp', sp, npang)
    call alloc ('mod_aom', 'angw', angw, npang)

  end subroutine allocate_space_for_surface

  subroutine deallocate_space_for_surface ()

    use mod_memory, only: free
    implicit none

    call free ('mod_aom', 'rlimsurf', rlimsurf)
    call free ('mod_aom', 'nlimsurf', nlimsurf)
    call free ('mod_aom', 'ct', ct)
    call free ('mod_aom', 'st', st)
    call free ('mod_aom', 'cp', cp)
    call free ('mod_aom', 'sp', sp)
    call free ('mod_aom', 'angw', angw)

  end subroutine deallocate_space_for_surface

  subroutine init_plot ()

    implicit none

    xnuc = 0.0_rp
    npang = 0_ip
    rmin = 0.0_rp
    rmax = 0.0_rp

  end subroutine

  subroutine vecprod (a, b, c)
  
    implicit none
    real(kind=rp), intent(in) :: a(3), b(3)
    real(kind=rp), intent(out) :: c(3)

    c(1) = a(2)*b(3) - a(3)*b(2) 
    c(2) = a(3)*b(1) - a(1)*b(3) 
    c(3) = a(1)*b(2) - a(2)*b(1)

  end subroutine
                                                                        
end module mod_plot
