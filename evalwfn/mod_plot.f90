module mod_plot

  use mod_prec, only: rp, ip
  implicit none
  private

  real(kind=rp) :: x0(3), x1(3), x2(3), rg(2)
  real(kind=rp) :: sx0, sy0, zx0, zx1, zy0, zy1
  integer(kind=ip) :: nx, ny, xy(2)

  public :: parse_plane

contains

  !> Calculate properties on a plane.
  subroutine parse_plane (fileroot)

    use mod_io, only: ferror, faterr, equal, isinteger, isreal, isword, &
                      uin, string, uout, lgetword, getline
    implicit none

    character(len=*), intent(in) :: fileroot

    integer(kind=ip) :: lp
    character(len=:), allocatable :: line, subline, word
    logical :: ok

    nx = 0
    ny = 0
    xy(1) = 1
    xy(2) = 2
    rg(1) = 0.01
    rg(2) = 0.01
    x0 = 0.0_rp
    x1 = 0.0_rp
    x2 = 0.0_rp
    sx0 = 1.0_rp
    sy0 = 1.0_rp
    zx0 = 0.0_rp
    zx1 = 0.0_rp
    zy0 = 0.0_rp
    zy1 = 0.0_rp

    ! Read the input file
    do while (getline(uin,line))
      lp = 1
      word = lgetword(line,lp)
      subline = line(lp:)

      if (equal(word,'#')) then
        continue
      
      ! Points that define the plane: origin, x1, x2
      else if (equal(word,'points')) then
        ok = isreal(x0(1), line, lp)
        ok = ok .and. isreal(x0(2), line, lp)
        ok = ok .and. isreal(x0(3), line, lp)
        ok = ok .and. isreal(x1(1), line, lp)
        ok = ok .and. isreal(x1(2), line, lp)
        ok = ok .and. isreal(x1(3), line, lp)
        ok = ok .and. isreal(x2(1), line, lp)
        ok = ok .and. isreal(x2(2), line, lp)
        ok = ok .and. isreal(x2(3), line, lp)
        if (ok) then
          write (uout,'(1x,a,1x,9(f5.2,1x))') string('# Plane x0,x1,x2 points'), x0, x1, x2
        else 
          call ferror ('mod_plot', 'wrong poins line', faterr)
        end if

      ! Number of points in each direction
      else if (equal(word,'npoints')) then
        ok = isinteger(nx, line, lp)
        ok = ok .and. isinteger(ny, line, lp)
        if (ok) then
          write (uout,'(1x,a,1x,i0,1x,i0)') string('# Plane wil use (nx,ny) points'), nx, ny
        else 
          call ferror ('mod_plot', 'wrong npoins line', faterr)
        end if
      
      ! Scale plane
      else if (equal(word,'scale')) then
        ok = isreal(sx0, line, lp)
        ok = ok .and. isreal(sy0, line, lp)
        if (ok) then
          write (uout,'(1x,a,1x,2(f4.2,1x))') string('# Scaling factors (x,y)'), sx0, sy0
        else 
          call ferror ('mod_plot', 'wrong scale line', faterr)
        end if

      ! Extend plane
      else if (equal(word,'extendx')) then
        ok = isreal(zx0, line, lp)
        ok = ok .and. isreal(zx1, line, lp)
        if (ok) then
          write (uout,'(1x,a,1x,2(f4.2,1x))') string('# Extend factors (x0,x1)'), zx0, zx1
        else 
          call ferror ('mod_plot', 'wrong extend line', faterr)
        end if

      else if (equal(word,'extendy')) then
        ok = isreal(zy0, line, lp)
        ok = ok .and. isreal(zy1, line, lp)
        if (ok) then
          write (uout,'(1x,a,1x,4(f4.2,1x))') string('# Extend factors (y0,y1)'), zy0, zy1
        else 
          call ferror ('mod_plot', 'wrong extend line', faterr)
        end if

      ! Coords
      else if (equal(word,'coords')) then
        ok = isinteger(xy(1), line, lp)
        ok = ok .and. isinteger(xy(2), line, lp)
        if (ok) then
          write (uout,'(1x,a,1x,4(i0,1x))') string('# Coords direction (x,y)'), xy(1), xy(2)
        else 
          call ferror ('mod_plot', 'wrong coords line', faterr)
        end if

      ! Range values
      else if (equal(word,'range')) then
        ok = isreal(rg(1), line, lp)
        ok = ok .and. isreal(rg(2), line, lp)
        if (ok) then
          write (uout,'(1x,a,1x,4(f4.2,1x))') string('# Coords direction (x,y)'), rg(1), rg(2)
        else 
          call ferror ('mod_plot', 'wrong coords line', faterr)
        end if

      ! End of input
      else if (equal(word,'endplane')) then
        exit

      else
        call ferror ('mod_plot', 'unknown option', faterr)

      end if
    end do

    call plot_plane (fileroot)

  end subroutine parse_plane

  subroutine plot_plane (fileroot)

    use mod_memory, only: alloc, free
    use mod_wfn, only: nmo
    use mod_fields, only: density, denmos
    use mod_io, only: uout, string, fopen_write, fclose
    implicit none

    character(len=*), intent(in) :: fileroot

    integer(kind=ip) :: ix, iy, luout
    real(kind=rp), allocatable :: ff(:,:), ffo(:,:,:)
    real(kind=rp) :: val, du, dv, uu(3), vv(3), xp(3), dmos(nmo)
    character(len=:), allocatable :: outfile

    ! root file name
    outfile = string(trim(fileroot)//".dat")
    write (uout,'(1x,a)') "# Computing data for PLANE "
    write (uout,'(1x,a)') string('# Data file output '//string(outfile))

    ! at least two points 
    nx = max(nx,2)
    ny = max(ny,2)

    ! Extend and scale, and set up the plane vectors.
    call plane_scale_extend (x0,x1,x2,sx0,sy0,zx0,zx1,zy0,zy1)
    uu = (x1-x0)/real(nx-1,rp)
    vv = (x2-x0)/real(ny-1,rp)
    du = norm2(uu)
    dv = norm2(vv)

    ! allocate space for field values on the plane
    call alloc ('plot_plane', 'ff', ff, nx, ny)
    call alloc ('plot_plane', 'ffo', ffo, nmo, nx, ny)
    do ix = 1,nx
      do iy = 1,ny
        xp = x0 + real(ix-1,rp)*uu + real(iy-1,rp)*vv
        call density (xp, val)
        call denmos (xp, dmos)
        ff(ix,iy) = val
        ffo(:,ix,iy) = dmos
      end do
    end do
    ! write info
    luout = fopen_write(outfile)
    write (luout,'(1x,a)') string("# Field x0,x1,x2 values")
    do ix = 1, nx
      do iy = 1,ny
        xp = x0 + real(ix-1,rp)*uu + real(iy-1,rp)*vv
        write (luout,'(1x,4(f15.10,1x),1p,1(e18.10,1x),1p,50(e18.10,1x),0p)') xp, ff(ix,iy), ffo(:,ix,iy)
      end do
      write (luout,*)
    end do
    call free ('plot_plane', 'ff', ff)
    call free ('plot_plane', 'ffo', ffo)
    call fclose (luout)

    ! print templates
    call gnuplot_template (fileroot)
    call latex_template (fileroot)

  end subroutine plot_plane

  !> Scale or extend the plane defined by points x0, x1, x2 
  !> with scale factors sxi, syi (default: 1) and extend factors
  !> zx0i, zx1i in the x direction (default: 0) and zy0i, zy1i in the
  !> y direction (default: 0).
  subroutine plane_scale_extend (x0,x1,x2,sxi,syi,zx0i,zx1i,zy0i,zy1i)

    use mod_io, only: ferror, faterr
    use mod_param, only: vsmall
    implicit none

    real(kind=rp), intent(inout) :: x0(3), x1(3), x2(3)
    real(kind=rp), intent(in), optional :: sxi, syi, zx0i, zx1i, zy0i, zy1i
    
    real(kind=rp) :: sx, sy, zx0, zx1, zy0, zy1
    real(kind=rp) :: ax(3), ay(3), ax0(3), ax1(3), ay0(3), ay1(3), dx, dy

    sx = 1.0_rp
    sy = 1.0_rp
    zx0 = 0.0_rp
    zx1 = 0.0_rp
    zy0 = 0.0_rp
    zy1 = 0.0_rp
    if (present(sxi)) sx = sxi
    if (present(syi)) sy = syi
    if (present(zx0i)) zx0 = zx0i
    if (present(zx1i)) zx1 = zx1i
    if (present(zy0i)) zy0 = zy0i
    if (present(zy1i)) zy1 = zy1i
    
    ax = (x1-x0)
    ay = (x2-x0)
    dx = norm2(ax)
    dy = norm2(ay)
    if (dx < vsmall .or. dy < vsmall) then
      call ferror('plane_scale_extend','zero-area plane', faterr)
    end if

    ax0 = -(ax/dx) * (0.5_rp*(sx-1.0_rp)*dx + zx0)
    ax1 =  (ax/dx) * (0.5_rp*(sx-1.0_rp)*dx + zx1)
    ay0 = -(ay/dy) * (0.5_rp*(sy-1.0_rp)*dy + zy0)
    ay1 =  (ay/dy) * (0.5_rp*(sy-1.0_rp)*dy + zy1)

    x0 = x0 + ax0 + ay0
    x1 = x1 + ax1 + ay0
    x2 = x2 + ax0 + ay1

  end subroutine plane_scale_extend

  subroutine latex_template(fileroot)

    use mod_io, only: uout, fopen_write, fclose, string
    use mod_wfn, only: nmo, occ
    implicit none

    character(len=*), intent(in) :: fileroot

    integer(kind=ip) :: luout, j
    character(len=:), allocatable :: outfile

    ! root file name
    write (uout,'(1x,a)') string('# Latex file output '//string(fileroot)//".tex")

    ! Se puede añadir la distancia
    ! norm2(v1-v2)

    outfile = string(trim(fileroot)//".tex")
    luout = fopen_write(outfile)
    write (luout,'(a)') "\begin{table}[h]"
    write (luout,'(a)') "\small"
    write (luout,'(a)') "\caption{Contribution of each chanel to total $\delta$}"
    write (luout,'(a)') "\label{tbl:nados}"
    write (luout,'(2x,a)',advance="no") "\begin{tabular*}{0.48\textwidth}{@{\extracolsep{\fill}}"
    do j = 1,nmo
      write (luout,'(a)',advance="no") "c"
    end do
    write (luout,'(a)') "}"
    write (luout,'(4x,a)') "\hline"
    write (luout,'(a)',advance="no") "Molecule "
    do j = 1,nmo
      write (luout,'(a,i0)',advance="no") " & ", j
    end do
    write (luout,'(a)') " \\"
    write (luout,'(4x,a)') "\hline"
    write (luout,'(a)',advance="no") "\ce{name} "
    do j = 1,nmo
      write (luout,'(a,f10.6)',advance="no") " & ", occ(j)*2.0
    end do
    write (luout,'(a)') " \\"
    write (luout,'(4x,a)') "\hline"
    write (luout,'(2x,a)') "\end{tabular*}"
    write (luout,'(a)') "\end{table}"
    call fclose (luout)

  end subroutine latex_template

  subroutine gnuplot_template (fileroot)

    use mod_io, only: uout, string, fopen_write, fclose
    use mod_wfn, only: ncent, xyz=>coords, atnam, nmo
    implicit none

    character(len=*), intent(in) :: fileroot

    integer(kind=ip) :: luout, i, j
    character(len=:), allocatable :: outfile
    character(len=2) :: twochar

    ! root file name
    write (uout,'(1x,a)') string('# Gnuplot file output '//string(fileroot)//"_xx.gnu")

    outfile = trim(fileroot)//"_00.gnu" 
    luout = fopen_write(outfile)
    write (luout,'(1x,a)') "#!/usr/bin/env gnuplot"
    write (luout,'(1x,a)') "set terminal pdf enhanced color font 'Helvetica,12'"
    write (luout,'(1x,a)') "set output '"//string(trim(fileroot))//"_00.pdf'"
    write (luout,'(1x,a)') "set pm3d map"
    write (luout,'(1x,a)') "set size ratio -1"
    write (luout,'(1x,a)') "set size square"
    write (luout,'(1x,a)') "set key off"
    write (luout,'(1x,a,f12.8,a,f12.8,a)') "set cbrange [",-1.0*rg(1),":",rg(2),"]"
    write (luout,'(1x,a)') "#set palette rgbformulae 22,13,-31"
    write (luout,'(1x,a)') "load 'paired.pal'"
    do i = 1,ncent
      write (luout,*) "set label '", atnam(i)(1:2), "' at ", &
                      xyz(i,xy(1)), ",", xyz(i,xy(2)), &
                      " front textcolor 'black' offset 0,0.5 point right font ',10'"
    end do
    write (luout,'(1x,a,i0,a,i0,a)') "splot '"//string(trim(fileroot))//".dat' u ($",xy(1),"):($",xy(2),"):($4)"
    call fclose (luout)

    do j = 1,nmo
      twochar = string(j,length=2,pad0=.true.)
      outfile = trim(fileroot)//"_"//string(twochar)//".gnu" 
      luout = fopen_write(outfile)
      write (luout,'(1x,a)') "#!/usr/bin/env gnuplot"
      write (luout,'(1x,a)') "set terminal pdf enhanced color font 'Helvetica,12'"
      write (luout,'(1x,a)') "set output '"//string(trim(fileroot))//"_"//string(twochar)//".pdf'"
      write (luout,'(1x,a)') "set pm3d map"
      write (luout,'(1x,a)') "set size ratio -1"
      write (luout,'(1x,a)') "set size square"
      write (luout,'(1x,a)') "set key off"
      write (luout,'(1x,a,f12.8,a,f12.8,a)') "set cbrange [",-1.0*rg(1),":",rg(2),"]"
      write (luout,'(1x,a)') "#set palette rgbformulae 22,13,-31"
      write (luout,'(1x,a)') "load 'paired.pal'"
      do i = 1,ncent
        write (luout,*) "set label '", atnam(i)(1:2), "' at ", &
                        xyz(i,xy(1)), ",", xyz(i,xy(2)), &
                        " front textcolor 'black' offset 0,0.5 point right font ',10'"
      end do
      write (luout,'(1x,a,i0,a,i0,a,i0,a)') "splot '"//string(trim(fileroot))//".dat' u ($",xy(1),"):($",xy(2),"):($",j+4,")*1.0"
      call fclose (luout)
    end do

  end subroutine gnuplot_template

end module mod_plot
