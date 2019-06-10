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
module mod_geom

  use mod_io, only: mline
  implicit none
  public

  real(kind=8) :: cq(3), cm(3)
  real(kind=8) :: moi(3,3), emoi(3), evmoi(3,3)
  real(kind=8) :: covx

  integer(kind=4) :: ncent
  real(kind=8), allocatable, dimension(:,:) :: xyz
  real(kind=8), allocatable, dimension(:) :: charge
  character(len=2), allocatable, dimension(:) :: symbol

  character(len=mline) :: xyzfile 

contains

  subroutine init_geom ()
                
    implicit none
  
    covx = 1.2d0
    
  end subroutine init_geom

  subroutine readxyz (filename)

    use mod_datatm, only: namatm, maxzat
    use mod_io, only: isreal, isword, faterr, equal, uout, &
                      mline, string
    implicit none
    character(len=*), intent(in) :: filename

    integer(kind=4) :: lu, i, j, idx
    character(len=mline) :: aux

    xyzfile = filename

    write (uout,'(1x,a)') string('#')
    write (uout,'(1x,a,1x,a)') string('# Reading xyz file :'), string(xyzfile)

    ! read number of atoms
    lu = 101
    open (lu,file=xyzfile,status='old')
    read (lu,*) ncent

    ! allocate space
    call allocate_space ()

    ! read atomic coordinates and atomic numbers
    write (uout,'(1x,a)') string('# Coordiantes of atoms')
    read (lu,*)
    do i = 1,ncent
      read (lu,*) symbol(i), xyz(i,:)
      do j = 1,maxzat
        if (equal(symbol(i),namatm(j))) charge(i) = real(j,8) 
      end do
      write (uout,'(1x,a,1x,i0,1x,a,1x,f4.1,1x,3f10.6)') &
             '#', i, symbol(i)(1:2), charge(i), xyz(i,:)
    end do
  
    ! close
    close (lu)

    if (index(xyzfile,'.') > 0) then
      idx = index(xyzfile,'.',.true.)
      aux = xyzfile(1:idx-1)
      xyzfile = trim(aux)
    end if
    xyzfile = trim(adjustl(xyzfile))

    call cmcq ()
    call get_connect ()

  end subroutine readxyz

  subroutine cmcq ()
    
    use mod_io, only: ferror, faterr, uout, string, mline
    use mod_linalg, only: jacobi
    use mod_datatm, only: wgatm

    implicit none
    real(kind=8), parameter :: mom_thresh = 1e-3

    real(kind=8) :: wt, zt, twt, tzt, wi, xx, yy, zz
    integer(kind=4) :: i, nrot
    character(len=mline) :: moltype
    logical :: same12, same13, same23, onezero, allzero

    cm = 0d0
    cq = 0d0
    wt = 0d0
    zt = 0d0
    moi = 0d0

    ! compute the center of mass and charge
    do i = 1,ncent
      twt = wgatm(int(charge(i)))
      wt = wt + twt
      cm(1) = cm(1) + twt*xyz(i,1)
      cm(2) = cm(2) + twt*xyz(i,2)
      cm(3) = cm(3) + twt*xyz(i,3)
      tzt = charge(i)
      zt = zt + tzt
      cq(1) = cq(1) + tzt*xyz(i,1)
      cq(2) = cq(2) + tzt*xyz(i,2)
      cq(3) = cq(3) + tzt*xyz(i,3)
    end do
    cm = cm / wt
    cq = cq / zt
    write (uout,'(1x,a,1x,3(1x,f10.6))') '# Center of mass', cm(:)
    write (uout,'(1x,a)') string('# Center of mass translated geometry')
    do i = 1,ncent
      write (uout,'(1x,a,1x,i0,1x,a,1x,f4.1,1x,3f10.6)') &
             '#', i, symbol(i)(1:2), charge(i), xyz(i,:) - cm(:)
    end do
    write (uout,'(1x,a,1x,3(1x,f10.6))') '# Center of charge', cq(:)
    write (uout,'(1x,a)') string('# Center of charge translated geometry')
    do i = 1,ncent
      write (uout,'(1x,a,1x,i0,1x,a,1x,f4.1,1x,3f10.6)') &
             '#', i, symbol(i)(1:2), charge(i), xyz(i,:) - cq(:)
    end do

    ! compute the inertia matrix and diagonalize it
    do i = 1,ncent
      wi = wgatm(int(charge(i)))
      xx = xyz(i,1) - cm(1)
      yy = xyz(i,2) - cm(2)
      zz = xyz(i,3) - cm(3)
      moi(1,1) = moi(1,1) + wi*(yy*yy + zz*zz)
      moi(2,2) = moi(2,2) + wi*(zz*zz + xx*xx)
      moi(3,3) = moi(3,3) + wi*(xx*xx + yy*yy)
      moi(1,2) = moi(1,2) - wi*xx*yy
      moi(2,3) = moi(2,3) - wi*yy*zz
      moi(3,1) = moi(3,1) - wi*zz*xx
    end do
    moi(2,1) = moi(1,2)
    moi(3,2) = moi(2,3)
    moi(1,3) = moi(3,1)
  
    call jacobi (moi,emoi,evmoi,nrot)
    if (nrot .lt. 0) call ferror ('cmcq','fail diag inertia', faterr)
    write (uout,'(1x,a)') '# Moment of inertia tensor (amu * A^2)'
    do i = 1,3
      write (uout,'(1x,a,1x,3(1x,f10.6))') '#', moi(i,:)
    end do
    write (uout,'(1x,a)') '# Principal moments of inertia (amu * A^2)'
    write (uout,'(1x,a,1x,3(1x,f10.6))') '#', emoi(:)

    same12 = are_same(emoi(1), emoi(2), mom_thresh)
    same13 = are_same(emoi(1), emoi(3), mom_thresh)
    same23 = are_same(emoi(2), emoi(3), mom_thresh)
    onezero = are_same(emoi(1), 0.0d0, mom_thresh)
    allzero = are_same(emoi(3), 0.0d0, mom_thresh)
    if (allzero) then
      moltype = 'monatomic'
    else if (onezero) then
      moltype = 'linear'
    else if (same13) then
      moltype = 'a spherical top'
    else if (same12 .or. same23) then
      moltype = 'a symmetric top'
    else
      moltype = 'an asymmetric top'
    end if
    write (uout,'(1x,a,1x,a)') '# The molecule is', string(moltype)

    contains

      logical function are_same(n1, n2, tol)

        implicit none
  
        real(kind=8), intent(in) :: tol
        real(kind=8), intent(in) :: n1, n2

        real(kind=8) :: comp

        are_same = .false.
        comp = abs((n2 - n1))

        if (comp <= tol) then
          are_same = .true.
        end if

        return 

      end function
  
  end subroutine

  subroutine allocate_space ()
  
    use mod_io, only: faterr, ferror
    implicit none

    integer(kind=4) :: ier

    allocate (xyz(ncent,3), stat=ier)
    if (ier.ne.0) then
      call ferror('mod_xyz', 'cannot allocate xyz', faterr)
    end if
    allocate (charge(ncent), stat=ier)
    if (ier.ne.0) then
      call ferror('mod_xyz', 'cannot allocate charge', faterr)
    end if
    allocate (symbol(ncent),stat=ier) 
    if (ier.ne.0) then
      call ferror('mod_xyz', 'cannot allocate symbol', faterr)
    end if
 
  end subroutine allocate_space

  subroutine deallocate_space ()
  
    use mod_io, only: faterr, ferror
    implicit none

    integer(kind=4) :: ier

    deallocate (xyz, stat=ier)
    if (ier.ne.0) then
      call ferror('mod_xyz', 'cannot deallocate xyz', faterr)
    end if
    deallocate (charge, stat=ier)
    if (ier.ne.0) then
      call ferror('mod_xyz', 'cannot deallocate charge', faterr)
    end if
    deallocate (symbol,stat=ier) 
    if (ier.ne.0) then
      call ferror('mod_xyz', 'cannot deallocate symbol', faterr)
    end if
 
  end subroutine deallocate_space

  ! Let x(i), i=1,n, the n eigenvalues of a real symmetric matrix A
  ! of logical size n x n and physical size np x np; i.e. the n solu-
  ! tions of the secular equation of A
  ! det | x I - A | = SUM (k=0,N) c(k) x^k = 0
  ! are x(1), x(2),..., x(n). c(k) are the cofficients of the charac-
  ! teristic polynomial and are the objective of this routine.
  subroutine polich (x, n, np, c)

    use mod_io, only: ferror, faterr
    implicit none
 
    integer(kind=4), intent(in) :: n, np
    real(kind=8) :: x(np), c(0:np)

    integer(kind=4) :: i, k
 
    if (n.gt.0) then
      c(0) = -x(1)
      c(1) = 1.0d0
      do i = 2,n
        c(i) = c(i-1)
        do k = i-1,1,-1
          c(k) = c(k-1)-x(i)*c(k)
        end do
        c(0) = -x(i)*c(0)
      end do
    else
       call ferror('polich','improper dimension', faterr)
    end if

  end subroutine

  subroutine get_connect ()

    use mod_io, only: string, faterr, ferror, uout
    use mod_datatm, only: covr
    use mod_linalg, only: jacobi
    use mod_futils, only: iqcksort
    implicit none
      
    integer(kind=4), parameter :: maxcoord = 20
    integer(kind=4), parameter :: infinity = -1

    logical :: connected
    real(kind=8) :: c(0:ncent), ankl
    character(len=2) :: this, oth(maxcoord)
    integer(kind=4) :: wh(ncent,maxcoord), nu, iclus(ncent), iclaux(ncent)
    real(kind=8) :: cnx(ncent,ncent), catom(ncent), v(ncent,ncent)  
    integer(kind=4) :: coord(ncent), bonded(ncent), iord(ncent)
    real(kind=8) :: xdis(3), rbond, covk, covm, dis2, d(ncent)
    integer(kind=4) :: k, m, ichm, ichk, i, j, nwh, nbonds, nclaux, cdis(ncent)
    integer(kind=4) :: p(0:ncent), nrot, madis(ncent,ncent), nclus
    character(len=1) :: labdis(ncent,ncent), digs(-1:18) 
    data digs / '-','-','1','2','3','4','5','6','7','8','9', &
                'a','b','c','d','e','f','g','h','i'/

    ! Compute connectivity matrix. Two atoms are bonded if the dis-
    ! tance between them is smaller the sum of their covalent radius 
    ! multiplied by a factor covx
    cnx = 0.0d0
    do k = 1,ncent-1
      do m = k+1,ncent
        xdis(:) = xyz(k,:)-xyz(m,:)
        dis2 = xdis(1)**2+xdis(2)**2+xdis(3)**2
        ichk = int(charge(k))
        ichm = int(charge(m))
        covk = covr(ichk)
        covm = covr(ichm)
        rbond = (covk+covm)*covx
        if (dis2 .lt. rbond*rbond) then 
          cnx(k,m) = 1.0d0
          cnx(m,k) = 1.0d0
        end if
      end do
    end do

    ! Compute the coordinations of all the atoms.
    nbonds = 0 
    do i = 1,ncent
      coord(i) = 0
      do j = 1,ncent
        nbonds = nbonds + nint(cnx(i,j))
        coord(i) = coord(i) + nint(cnx(i,j))
      end do
    end do
    nbonds = nbonds/2

    ! Write coordination indices and connectivity matrix.
    write (uout,'(1x,a,f8.4,1x,a)') '# Bonding Criterion', covx,'x Sum of covalent radii'
    write (uout,'(1x,a)') string('# Coordination indices and Connectivity Matrix')
    write (uout,'(1x,a)') string('# --------------------------------------------')
    do k = 1,ncent
      nwh = 0
      do m = 1,ncent
        if (nint(cnx(k,m)).eq.1) then
          nwh = nwh + 1
          if (nwh.gt.MAXCOORD) then
            call ferror ('connect', 'increase the value of maxcoord', faterr)
          end if
          wh(k,nwh) = m
        end if
      end do
      this(1:2) = symbol(k)(1:2)
      do m = 1,nwh
        oth(m)(1:2) = symbol(wh(k,m))(1:2)
      end do
      if (ncent.lt.10) then
        write (uout,221) this,k,coord(k),(oth(m)(1:2),wh(k,m),m=1,nwh)
      else if (ncent.lt.100) then
        write (uout,222) this,k,coord(k),(oth(m)(1:2),wh(k,m),m=1,nwh)
      else if (ncent.lt.1000) then
        write (uout,222) this,k,coord(k),(oth(m)(1:2),wh(k,m),m=1,nwh)
      else
        write (uout,224) this,k,coord(k),(oth(m)(1:2),wh(k,m),m=1,nwh)
      end if
    end do
    bonded = coord
    write (uout,'(1x,a,1x,i0)') '# Number of bonds ', nbonds

    ! Order the coordinations of all the atoms.
    forall (i=1:ncent) iord(i) = i
    call iqcksort (coord,iord,1,ncent)
    forall (i=1:ncent) catom(i) = coord(iord(i))

    ! Diagonalize the connectivity matrix.
    call jacobi (cnx,d,v,nrot)
    if (nrot .lt. 0) call ferror ('connect','fail diag connectivity', faterr)
    ! Detemine the characteristic polynomium.
    call polich (d,ncent,ncent,c)
    forall (i=0:ncent) p(i) = nint(c(i))
 
    ! Compute the distance matrix. Algorithm: See Chemical Graph Theory.
    ! Chapter 2 by O. E. Polansky, Section 2.9, Eq. (51).
    connected = .true.
    do i = 1,ncent-1
      madis(i,i) = 0
      do j = i+1,ncent
        do nu = 1,ncent
          ankl = 0.0d0
          do k = 1,ncent
            ankl = ankl + v(i,k)*v(j,k)*d(k)**nu
          end do
          if (nint(ankl).gt.0) then
            madis(i,j) = nu
            madis(j,i) = nu
            go to 4
          end if
        end do
        ! This cluster is a disconnected graph.
        connected = .false.
        madis(i,j) = infinity
        madis(j,i) = infinity
   4    continue
      end do
    end do
    madis(ncent,ncent) = 0

    ! Compute the sum of distance indices for all the atoms and order them.
    forall (i=1:ncent) coord(i) = sum(madis(i,:))
 
    ! Write the sum of indices of distances and the distance matrix.
    if (connected) then
      write (uout,'(1x,a)') '# Connected graph'
    else
      write (uout,'(1x,a)') '# Non-Connected graph'
    end if
    write (uout,'(1x,a)') '# Distance Matrix: "-" means non connected atoms'
    if (ncent.le.100) then
      labdis(1:ncent,1:ncent) = digs(0)
      do i = 1,ncent
        do j = 1,ncent
          labdis(i,j) = digs(madis(i,j))
        end do
        write (uout,'(1x,a,100a4)') '#', (labdis(i,j),j=1,ncent)
      end do
    end if
 
    ! Order the sums of indices of distances.
    forall (i=1:ncent) iord(i) = i
    call iqcksort (coord,iord,1,ncent)
    forall (i=1:ncent) cdis(i) = coord(iord(i))
 
    ! Determine and identify the different non-connected clusters.
    forall (i=1:ncent) iclus(i) = 0
    nclus = 0
    do i = 1,ncent
      if (iclus(i).eq.0) then
        nclus = nclus + 1
        iclus(i) = nclus
        do j = i+1,ncent
          if (madis(i,j).ne.infinity) iclus(j) = nclus
        end do
      end if
    end do
    write (uout,'(1x,a,1x,i0,1x,a)') '# Molecule made of', nclus, 'non-connected fragmets'

    write (uout,'(1x,a)') string('#')
    do i = 1,nclus
      write (uout,'(1x,a,1x,i0)') string('# *** Internal coordinates of fragment'), i
      nclaux = 0
      do j = 1,ncent
        if (iclus(j).eq.i) then 
          nclaux = nclaux + 1
          iclaux(nclaux) = j
        end if
      end do
      write (uout,'(1x,a,1x,i0,1x,a,1x,i0,1x,a)') '# Fragment', i, 'contains', nclaux, 'atoms'
      write (uout,9) (iclaux(j),j=1,nclaux)
      call internal (i, iclaux, nclaux)
      write (uout,'(1x,a)') string('#')
    end do

    ! Formats
221 format (1x,'# (',a2,i1,') Coor = ',i1,' -->',30(1x,'(',a2,i1,')'))
222 format (1x,'# (',a2,i2,') Coor = ',i2,' -->',30(1x,'(',a2,i2,')'))
224 format (1x,'# (',a2,i4,') Coor = ',i4,' -->',30(1x,'(',a2,i4,')'))
9   format (1x,'#',20(1x,I4))

  end subroutine get_connect

  subroutine internal (id, iclaux, nclaux)

    use mod_io, only: uout, string, mline
    use mod_math, only: dist, a123
    implicit none
    integer(kind=4), intent(in) :: nclaux, id
    integer(kind=4), dimension(nclaux), intent(in) :: iclaux

    !real(kind=8) :: dis, a(3), b(3), c(3), d(3), angle
    !integer(kind=4) :: i, j, k, idx, jdx, kdx, l, ldx, lu
    integer(kind=4) :: i, idx, lu
    character(len=mline) :: filename

    !if (nclaux .eq. 1) then
    !  write (uout,'(1x,a)') string('# Fragment contains only one atom')
    !else
    !  do i = 1,nclaux
    !    idx = iclaux(i)
    !    a = xyz(idx,:)
    !    do j = i+1,nclaux
    !      jdx = iclaux(j)
    !      b = xyz(jdx,:)
    !      dis = dist (a,b)
    !      write (uout,'(1x,a,2(1x,i0),1x,2(1x,a),2x,f10.6)') &
    !             string('# Distance for atoms'), idx, jdx, &
    !             string(symbol(idx)), string(symbol(jdx)), dis
    !    end do
    !  end do
    !  if (nclaux.ge.3) then
    !    do i = 1,nclaux
    !      idx = iclaux(i)
    !      a = xyz(idx,:)
    !      do j = i+1,nclaux
    !        jdx = iclaux(j)
    !        b = xyz(jdx,:)
    !        do k = j+1,nclaux
    !          kdx = iclaux(k)
    !          c = xyz(kdx,:)
    !          angle = a123 (a, c, b)
    !          write (uout,'(1x,a,3(1x,i0),1x,3(1x,a),2x,f10.6)') &
    !                 string('# Angle for atoms'), idx, jdx, kdx, &
    !                 string(symbol(idx)), string(symbol(jdx)), &
    !                 string(symbol(kdx)), angle
    !        end do
    !      end do
    !    end do
    !  end if  
    !  if (nclaux.ge.4) then
    !    do i = 1,nclaux
    !      idx = iclaux(i)
    !      a = xyz(idx,:)
    !      do j = i+1,nclaux
    !        jdx = iclaux(j)
    !        b = xyz(jdx,:)
    !        do k = j+1,nclaux
    !          kdx = iclaux(k)
    !          c = xyz(kdx,:)
    !          do l = k+1,nclaux
    !            ldx = iclaux(l)
    !            d = xyz(ldx,:)
    !          end do
    !        end do
    !      end do
    !    end do
    !  end if
    !end if

    lu = 401
    filename = string(xyzfile)//'_'//string(id,length=2,pad0=.true.)//'.xyz'
    write (uout,'(1x,a,1x,a)') string('# Write xyz file'), string(filename)
    open (unit=lu, file=filename)
    write (lu,'(i0)') nclaux
    write (lu, '(a)') string('0 1')
    do i = 1,nclaux
      idx = iclaux(i)
      write (lu,'(a,1x,3(1x,f12.6))') string(symbol(idx)), xyz(idx,:)
    end do
    close (lu)

  end subroutine internal

end module mod_geom
