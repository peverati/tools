! Copyright (c) 2018 
! Jose Luis Casals Sainz <jluiscasalssainz@gmail.com>,
! Ángel Martín Pendás <angel@fluor.quimica.uniovi.es> and 
! Evelio Francisco Miguelez <evelio@uniovi.es>. 
!
! xprom is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or (at
! your option) any later version.
! 
! xprom is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
subroutine pair ()

  use mod_xprom
  use mod_io, only: getline, lgetword, ferror, faterr, &
                    isinteger, uin, equal, string
  implicit none
  
  real(kind=8) :: xnetenergy, xcnetene, conetene, eenetene, totq
  real(kind=8) :: nnnetene, xinterenergy, xclenergy, xcinterenergy 
  real(kind=8) :: coulintenergy, eerepenergy, nnrepenergy, xkinenergy 
  real(kind=8) :: evacuotot, effenergy, addenergy, totel
  real(kind=8) :: corrinterenergy, xcorrenergy, xnetcorrenergy
  integer(kind=4) :: lp, intg, i, j, ig1, ig2, igroup, ii, k, kk, lwrite
  character(len=:), allocatable :: line, subline, word

  do while (getline(uin,line))
    lp = 1
    word = lgetword(line,lp)
    subline = line(lp:)

    if (equal(word,'#')) then
      continue
    
    else if (equal(word,'ga')) then
      i = 0
 3    if (isinteger(intg, line, lp)) then
        if (intg.gt.0 .and. intg.le.ncent) then
          do j = 1,i
            if (intg.eq.ig(j,1)) goto 3
          end do
          i = i + 1
          ig(i,1) = intg
        else 
          call ferror('pair', 'wrong ga format', faterr)
        end if
        goto 3
      end if
      ngroup(1) = i
    
    else if (equal(word,'gb')) then
      i = 0
 4    if (isinteger(intg, line, lp)) then
        if (intg.gt.0 .and. intg.le.ncent) then
          do j = 1,i
            if (intg.eq.ig(j,2)) goto 4
          end do
          i = i + 1
          ig(i,2) = intg
        else 
          call ferror('pair', 'wrong gb format', faterr)
        end if
        goto 4
      end if
      ngroup(2) = i
    
    else if (equal(word,'endgroup')) then
      exit
    
    else
      call ferror('pair', 'wrong group/pair format', faterr)

    end if
  end do

  lwrite = 300
  open (lwrite,file='results.txt')
  write (lwrite,'(1x,64("-"))')
  write (lwrite,'(1x,a)') string('# The following pair will be analyzed')
  write (lwrite,'(1x,a,100(1x,i0))') '# GROUP 1 formed by atoms : ', ig(1:ngroup(1),1)
  write (lwrite,'(1x,a,100(1x,i0))') '# GROUP 2 formed by atoms : ', ig(1:ngroup(2),2)
  write (lwrite,'(1x,64("-"))')

  ! Write interaction 
  do igroup = 1,2
    if (igroup.eq.1) then
      ig1 = 1
      ig2 = 2
    else
      ig1 = 2
      ig2 = 1
    end if
    xnetenergy = 0d0
    xcnetene = 0d0
    conetene = 0d0
    eenetene = 0d0
    nnnetene = 0d0
    xinterenergy = 0d0
    xclenergy = 0d0
    xcinterenergy = 0d0
    coulintenergy = 0d0
    eerepenergy = 0d0
    nnrepenergy = 0d0
    xkinenergy = 0d0
    evacuotot = 0d0
    effenergy = 0d0
    addenergy = 0d0
    corrinterenergy = 0d0
    xcorrenergy = 0d0
    xnetcorrenergy = 0d0
    totel = 0d0
    totq = 0d0
    do i = 1,ngroup(ig1)
      k = ig(i,ig1)
      totel = totel + elec(k)
      totq = totq + qq(k)
      xnetenergy = xnetenergy + fnet(k)
      xkinenergy = xkinenergy + fkin(k)
      xcnetene = xcnetene + feexc(k)
      conetene = conetene + feec(k)
      eenetene = eenetene + fee(k)
      if (t2) xnetcorrenergy = xnetcorrenergy + fcorr(k)
      do ii = 1,ngroup(ig1)
        kk = ig(ii,ig1)
        if (k.ne.kk) then
          xnetenergy = xnetenergy + 0.5d0*fintab(k,kk)
          xcnetene = xcnetene + 0.5d0*fxcab(k,kk)
          conetene = conetene + 0.5d0*fcoulab(k,kk)
          eenetene = eenetene + 0.5d0*feeab(k,kk)
          nnnetene = nnnetene + 0.5d0*fnnab(k,kk)
          if (t2) xnetcorrenergy = xnetcorrenergy + 0.5d0*fcorrab(k,kk)
        end if
      end do
      do ii = 1,ngroup(ig2)
        kk = ig(ii,ig2)
        xinterenergy = xinterenergy + fintab(k,kk)
        xclenergy = xclenergy + fxvclab(k,kk)
        xcinterenergy = xcinterenergy + fxcab(k,kk)
        coulintenergy = coulintenergy + fcoulab(k,kk)
        eerepenergy = eerepenergy + feeab(k,kk)
        nnrepenergy = nnrepenergy + fnnab(k,kk)   
        if (t2) then
          corrinterenergy = corrinterenergy + fcorrab(k,kk)
        end if
      end do
    end do
    if (t2) then
      write (lwrite,803) ig1,totel,                     &
                         ig1,totq,                      &
                         ig1,xnetenergy,                &
                         ig1,xkinenergy,                &
                         ig1,xnetcorrenergy,            &
                         ig1,ig2,xinterenergy,          &
                         ig1,ig2,xclenergy,             &
                         ig1,ig2,xcinterenergy,         &
                         ig1,ig2,coulintenergy,         &
                         ig1,ig2,eerepenergy,           &
                         ig1,ig2,corrinterenergy,       &
                         ig1,ig2,nnrepenergy,           &
                         ig1,xcnetene,                  &
                         ig1,conetene,                  &
                         ig1,eenetene,                  &
                         ig1,nnnetene                    
    else 
      write (lwrite,802) ig1,totel,                     &
                         ig1,totq,                      &
                         ig1,xnetenergy,                &
                         ig1,xkinenergy,                &
                         ig1,ig2,xinterenergy,          &
                         ig1,ig2,xclenergy,             &
                         ig1,ig2,xcinterenergy,         &
                         ig1,ig2,coulintenergy,         &
                         ig1,ig2,eerepenergy,           &
                         ig1,ig2,nnrepenergy,           &
                         ig1,xcnetene,                  &
                         ig1,conetene,                  &
                         ig1,eenetene,                  &
                         ig1,nnnetene                    
    end if
    write (lwrite,'(1x,64("-"))')
  end do

  close (lwrite)

802  format (                                                          &
      ' GROUP ', i2,'   : Number of electrons             = ',F13.6,/, &
      ' GROUP ', i2,'   : Total charge                    = ',F13.6,/, &
      ' GROUP ', i2,'   : Net Energy                      = ',F13.6,/, &
      ' GROUP ', i2,'   : Kinetic Energy                  = ',F13.6,/, &
      ' GROUPS',2i2,' : Total Interaction Energy        = ',F13.6,/,   &
      ' GROUPS',2i2,' : Classic Interaction Energy      = ',F13.6,/,   &
      ' GROUPS',2i2,' : XC Interaction Energy           = ',F13.6,/,   &
      ' GROUPS',2i2,' : Coulomb (ee) Interaction Energy = ',F13.6,/,   &
      ' GROUPS',2i2,' : El-El Interaction Energy        = ',F13.6,/,   &
      ' GROUPS',2i2,' : Nuc-Nuc Interaction Energy      = ',F13.6,/,   &
      ' GROUP ', i2,'   : Intragroup XC Energy            = ',F13.6,/, &
      ' GROUP ', i2,'   : Intragroup Coulomb Energy       = ',F13.6,/, &
      ' GROUP ', i2,'   : Intragroup El-El Repulsion      = ',F13.6,/, &
      ' GROUP ', i2,'   : Intragroup Nuc-Nuc Repulsion    = ',F13.6)

803  format (                                                          &
      ' GROUP ', i2,'   : Number of electrons             = ',F13.6,/, &
      ' GROUP ', i2,'   : Total charge                    = ',F13.6,/, &
      ' GROUP ', i2,'   : Net Energy                      = ',F13.6,/, &
      ' GROUP ', i2,'   : Kinetic Energy                  = ',F13.6,/, &
      ' GROUP ', i2,'   : Correlation Energy              = ',F13.6,/, &
      ' GROUPS',2i2,' : Total Interaction Energy        = ',F13.6,/,   &
      ' GROUPS',2i2,' : Classic Interaction Energy      = ',F13.6,/,   &
      ' GROUPS',2i2,' : XC Interaction Energy           = ',F13.6,/,   &
      ' GROUPS',2i2,' : Coulomb (ee) Interaction Energy = ',F13.6,/,   &
      ' GROUPS',2i2,' : El-El Interaction Energy        = ',F13.6,/,   &
      ' GROUPS',2i2,' : Correlation Interaction Energy  = ',F13.6,/,   &
      ' GROUPS',2i2,' : Nuc-Nuc Interaction Energy      = ',F13.6,/,   &
      ' GROUP ', i2,'   : Intragroup XC Energy            = ',F13.6,/, &
      ' GROUP ', i2,'   : Intragroup Coulomb Energy       = ',F13.6,/, &
      ' GROUP ', i2,'   : Intragroup El-El Repulsion      = ',F13.6,/, &
      ' GROUP ', i2,'   : Intragroup Nuc-Nuc Repulsion    = ',F13.6)

end subroutine pair
