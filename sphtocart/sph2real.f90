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
!
! This program finds the coefficients that gives a real solid 
! harmonic (r^l S_lm (theta,phi) as a linear combination of terms 
! of the form x^i y^j z^k
program sph2real

  use mod_io, only: ferror, faterr
  use mod_param, only: init_param, fact, comb
  implicit none
  
  integer(kind=4), parameter :: maxl = 9
  integer(kind=4), parameter :: maxterm = (maxl+1)*(maxl+2)/2

  real(kind=8) :: coef
  real(kind=8) :: vm, na, nc, xlm, v, c(0:maxl,-maxl:maxl,maxterm)
  integer(kind=4) :: t, u, am, delta, vvm, nterm, i, ia
  integer(kind=4) :: ipow(0:maxl,-maxl:maxl,maxterm,3), ib
  integer(kind=4) :: nterms(0:maxl,-maxl:maxl), ic, j, l, m  
 
  call init_param ()
 
  write (6,*) ' L value ? '
  read (5,*) l
  if (l.lt.0) l = 0
  if (l.gt.maxl) l = maxl
  do m = -l,+l
    am = abs(m)
    na = (l-am)/2d0
    vm = 0d0
    if (m.lt.0) vm = 0.5d0
    nc = int(real(abs(m),4)/2d0-vm)+vm
    delta = 0d0
    if (m.eq.0) delta = 1d0
    xlm = sqrt(2d0*fact(l+am)*fact(l-am)/2**delta)/2**am/fact(l)
    nterm = 0
    do t = 0,int(na)
      do u = 0,t
        do v = vm,nc,1d0
          nterm = nterm + 1
          if (nterm.gt.maxterm) then
            call ferror('sph2real','maxterm parameter', faterr)
          end if
          vvm = int(v-vm)
          coef = (-1)**(t+vvm)/4d0**t*comb(l,t)* &
                 comb(l-t,am+t)*comb(t,u)*comb(am,int(2*v))
          coef = coef*xlm
          ia = 2*t+am-2*u-int(2*v)
          ib = 2*u+int(2*v)
          ic = l-2*t-am
          ipow(l,m,nterm,1) = ia
          ipow(l,m,nterm,2) = ib
          ipow(l,m,nterm,3) = ic
          c(l,m,nterm) = coef
        end do
      end do
    end do
    nterms(l,m) = nterm
    write (6,300) l,m,nterm,(c(l,m,i),&
    (ipow(l,m,i,j),j=1,3),i=1,nterms(l,m))
  end do

 300 format (1x,'(l,m) = (',2I3,' )',15x,'NTERMS = ',I3, &
             100(/1x,F26.20,' * X^',I1,' * Y^',I1,' * Z^',I1))

end program
