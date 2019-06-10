c---------------------------------------------------------------------
      subroutine wrsurf (epsilon,lsu,nr,npang,npth,npph,iqudr,iqudt,
     &       mapr,mapt,mmap,mpr,largwr,alfm,linter,rmesh,
     &       ct,st,cp,sp,files)
c
c...............................................................
c
      USE       mod_prec, only: ip, dp
      USE       space_for_wfnbasis
      USE       space_for_surface
      USE       space_for_sym
      USE       space_for_bicen
      USE       space_for_rinters
      include  'implicit.inc'
      include  'param.inc'
      include  'wfn.inc'
      include  'integ.inc'
      include  'io.inc'
      include  'sym.inc'
      include  'point.inc'  ! For OPENMP
      include  'error.inc'
      include  'parallel.inc'
c
      real(kind=dp) :: rmesh(0:nr,ncent)
      real(kind=dp) :: rsurf(minter,2), ct(npang), st(npang) 
      real(kind=dp) :: cp(npang), sp(npang)
      integer(kind=ip) :: linter(maxtrial,ncent)
      integer(kind=ip) :: iqudr, iqudt, mapr, mapt, mmap, mpr
      real(kind=dp) :: alfm,mid(3) 
      integer(kind=ip) :: leng
      logical :: largwr
      character(len=132) :: fneq
      character(len=132) :: filesh
      character(len=*)   :: files
      character(len=1)   :: digs(0:9)
      character(len=4)   :: d4
      character(len=4)   :: suffix
c
c Di+ variables for surface triangulation
c
      real(kind=dp) :: r
      real(kind=dp),allocatable,dimension (:)    :: x,y,z        
c
      data digs / '0', '1', '2', '3', '4', '5', '6', '7', '8', '9'/
c
c-----------------------------------------------------------------------
c
      files=files(1:leng(files))
      if (largwr) then
        filesh=files(1:leng(files))//".txt"
        filesh=filesh(1:leng(filesh))
        lsh=999
      end if
c
      do i=1,ncent
        if (lstart(i)) then
          do j=1,ncent
            if (idx(j,1).eq.idx(i,1)) then
              lstart(j)=.true.
              do k=1,nrsearch(i)
                rstart(j,k)=rstart(i,k)
              enddo
            endif
          enddo
        endif
      enddo
c
c.....If surfall is .true. the number of atoms for which we want to 
c     obtain the interatomic surface (neqscalc) coincides with the 
c     total number of non equivalent atoms (neq). Then, the indices
c     of the surfaces to be computed are icalcsurf(i)=i, i=1,...,neq.
c
      if (surfall) then
        neqscalc=neq
        do i=1,neqscalc
          icalcsurf(i)=i
        enddo
      else
c
c.....Otherwise, the actual values of icalcsurf(i) (i=1,neqsurf) are
c     obtained from the values of the insurf(i) (i=1,neqsurf) array
c     which has been read in from the input file in 'promolden.f'
c
        neqscalc=0
        do i=1,neqsurf
          ii=insurf(i)
          if (neqscalc.eq.0) then
            neqscalc=1
            icalcsurf(1)=idx(ii,1)
          else
            do j=1,neqscalc
              if (icalcsurf(j).eq.idx(ii,1)) goto 10
            enddo
            neqscalc=neqscalc+1
            icalcsurf(neqscalc)=idx(ii,1)
 10         continue
          endif
        enddo
c
c.......order the indices
c
        do i=1,neqscalc
          do j=i+1,neqscalc
            if (icalcsurf(j).lt.icalcsurf(i)) then
              isave=icalcsurf(i)
              icalcsurf(i)=icalcsurf(j)
              icalcsurf(j)=isave
            endif
          enddo
        enddo
      endif
c
c.....Test whether the file containing the interatomic surfaces of all
c     the ineqs exist or not.
c
      open (lsu,file=files,status='old',form='unformatted',iostat=ierr)
c
c.....if lawgwr open a new file
c
      if (largwr) open (lsh,file=filesh)
c
c.....The file does not exist.
c
      if (ierr.ne.0) then
c
c.......Read surfaces of atoms for which this surface actually exists
c
        do i=1,neq
          i1=i/1000
          i2=mod(i,1000)
          i3=i2/100
          i4=mod(i2,100)
          i5=i4/10
          i6=mod(i4,10)
          d4=digs(i1)//digs(i3)//digs(i5)//digs(i6)
          fneq=files(1:leng(files))//"-"//d4
          open (lsu,file=fneq,status='old',form='unformatted',
     &          iostat=nerr)
          if (nerr.eq.0) then
c
c...........Read surface for ineq i.
c
            read (lsu) ntp1,ntp2,ntp3,ntp4,ntp5,ntp6,ntp7,ntp8,
     &                 ntp9,ntp10,ntp11,temp1
            if (leb) then
              if (ntp3.ne.0.or.ntp4.ne.0) then
                write (uerr,200)
                call error('wrsurf','non compatible surface',faterr)
              endif
            else
              if (ntp3.ne.npth) then
                write (uerr,201) ntp3,npth
                call error('wrsurf','non compatible surface',faterr)
              endif
              if (ntp4.ne.npph) then
                write (uerr,202) ntp4,npph
                call error('wrsurf','non compatible surface',faterr)
              endif
            endif
c
c...........Test similarity between input values and values in files
c
            if (ntp1 .ne.neq  ) write (uerr,203) ntp1,neq
            if (ntp2 .ne.npang) write (uerr,204) ntp2,npang
            if (ntp5 .ne.nr   ) write (uerr,205) ntp5,nr
            if (ntp6 .ne.iqudr) write (uerr,206) ntp6,iqudr
            if (ntp7 .ne.iqudt) write (uerr,207) ntp7,iqudt
            if (ntp8 .ne.mapr ) write (uerr,208) ntp8,mapr
            if (ntp9 .ne.mapt ) write (uerr,209) ntp9,mapt
            if (ntp10.ne.mmap ) write (uerr,210) ntp10,mmap
            if (ntp11.ne.mpr  ) write (uerr,211) ntp11,mpr
            if (temp1.ne.alfm ) write (uerr,212) temp1,alfm
c
            if (ntp1.ne.neq  .or.ntp2 .ne.npang.or.ntp5 .ne.nr  .or.
     &          ntp6.ne.iqudr.or.ntp7 .ne.iqudt.or.ntp8 .ne.mapr.or.
     &          ntp9.ne.mapt .or.ntp10.ne.mmap .or.ntp11.ne.mpr .or.
     &         temp1.ne.alfm) then
              call error('wrsurf','surf not compatible',faterr)
            endif
c
            read (lsu) (nlimsurf(i,j),j=1,npang)
            do j=1,npang
               nsurf=nlimsurf(i,j)
               read (lsu) (rlimsurf(i,j,k),k=1,nsurf)
            enddo
            read (lsu) rmsurftemp
            if (rmsurftemp.ne.rmsurf(i)) then
               write (uerr,213) i,rmsurf(i),rmsurftemp
               call error('wrsurf','surf not compatible',faterr)
            endif
            close (lsu)
            surfineq(i)=.true.
          endif
        enddo
c
c.......Write surface for all the non-equivalent atoms
c
        surfineq = .false.
        do ii=1,neqscalc
          i=icalcsurf(ii)
          i1=i/1000
          i2=mod(i,1000)
          i3=i2/100
          i4=mod(i2,100)
          i5=i4/10
          i6=mod(i4,10)
          d4=digs(i1)//digs(i3)//digs(i5)//digs(i6)
          fneq=files(1:leng(files))//"-"//d4
          open (lsu,file=fneq,status='old',form='unformatted',
     &          iostat=nerr)
c
c.........Write surface for ineq i.
c
          if (nerr.ne.0) then
            rmsurf(i)=RMAXSURF 
            inuc=ineq(i)
            xnuc(1)=xyzrho(inuc,1)
            xnuc(2)=xyzrho(inuc,2)
            xnuc(3)=xyzrho(inuc,3)
            call cpu_time (time1)
!$          time1=omp_get_wtime()
            write (uout,*) '# Finding SURFACE for ineq:', i
            iper=0
            sten=npang/10d0
            xten=sten
! if the loop is over batches be careful, change kk to private
*           do kk=1,npang,ichunk
!$omp parallel  
!$omp& default(none)
!$omp& private(j,k,nsurf,rsurf)
!$omp& shared(npang,xsurface,ct,st,cp,sp,epsilon,i,rmesh,
!$omp& nr,linter,rlimsurf,nlimsurf,sten,xten,iper,uout,
!$omp& ichunk,kk)
!$omp do
!$omp& schedule(dynamic) 
*           do kk=1,npang,ichunk
*
*           block_length=min(ichunk,npang-kk+1)
*           do j=1,block_length
*
*           do j=kk,min(kk+ichunk-1,npang)
            do j=1,npang
              if (.not.xsurface) call surf1 (ct(j),st(j),cp(j),sp(j),
     &             epsilon,rsurf,nsurf)
              if (     xsurface) call surf2 (rmesh(0,i),linter(1,i),
     &             ct(j),st(j),cp(j),sp(j),
     &             epsilon,nr,rsurf,nsurf)
              if (nsurf.gt.minter)
     &          call error('wrsurf',
     &               'increase MINTER parameter in integ.inc',faterr) !explicit stop not allowed in parallel
              do k=1,nsurf
                 rlimsurf(i,j,k)=rsurf(k,2)
              enddo
              nlimsurf(i,j)=nsurf
**            if (j.ge.xten) then
**              xten=xten+sten
**              iper=iper+10
**              write (uout,'(1x,a, i4)') '# surface (%):', iper
**              call flush_unit (uout)
**            endif
            enddo
*           enddo
!$omp end do nowait
!$omp end parallel
*           enddo
            call cpu_time (time2)
!$          time2=omp_get_wtime()
            write (uout,112) time2-time1
            open (lsu,file=fneq,status='new',form='unformatted')
            if (leb) then
               write (lsu) neq,npang,0,0,nr,iqudr,iqudt,mapr,mapt,
     &                     mmap,mpr,alfm
               if (largwr) write (lsh,111) 
     &           i,neq,npang,nr,iqudr,iqudt,mapr,mapt,mmap,mpr,alfm
            else
               write (lsu) neq,npang,npth,npph,nr,iqudr,iqudt,mapr,mapt,
     &                     mmap,mpr,alfm
               if (largwr) write (lsh,113) i,neq,npang,npth,npph,nr,
     &                     iqudr,iqudt,mapr,mapt,mmap,mpr,alfm
            endif
            write (lsu) (nlimsurf(i,j),j=1,npang)
            if (largwr) write (lsh,114) (nlimsurf(i,j),j=1,npang)
            if (largwr) write (lsh,117) 
            do j=1,npang
               nsurf=nlimsurf(i,j)
               write (lsu) (rlimsurf(i,j,k),k=1,nsurf)
               if (largwr) write (lsh,115) ct(j),st(j),cp(j),sp(j),
     &           (rlimsurf(i,j,k),k=1,nsurf)
            enddo
            write (lsu) rmsurf(i)
            if (largwr) write (lsh,116) rmsurf(i)
          endif
          close (lsu)
          surfineq(i)=.true.
        enddo
        if (largwr) close(lsh)
c
c.......If the interatomic surface has been computed for all the ineqs
c       create an unique surface file.
c
        allsurfaces=.true.
        do i=1,neq
          if (.not.surfineq(i)) allsurfaces=.false.
        enddo
        if (allsurfaces) then
          do i=1,neq
            i1=i/1000
            i2=mod(i,1000)
            i3=i2/100
            i4=mod(i2,100)
            i5=i4/10
            i6=mod(i4,10)
            d4=digs(i1)//digs(i3)//digs(i5)//digs(i6)
            fneq=files(1:leng(files))//"-"//d4
            open (lsu,file=fneq,status='old',iostat=nerr)
            if (nerr.ne.0) then
              write (uout,*) '# !!! Warning: Trying to delete a !!!'
              write (uout,'(a)') ' # !!! non existing file : ',fneq
            else
              close (lsu,status='delete')
            endif
          enddo
c
          open (lsu,file=files,status='new',form='unformatted')
          if (     leb) write (lsu) neq,npang,0,0
          if (.not.leb) write (lsu) neq,npang,npth,npph
          do i=1,neq
            write (lsu) (nlimsurf(i,j),j=1,npang)
          enddo
          do i=1,neq
            do j=1,npang
              nsurf=nlimsurf(i,j)
              write (lsu) (rlimsurf(i,j,k),k=1,nsurf)
            enddo
          enddo
          write (lsu) (rmsurf(i),i=1,neq)
          close (lsu)
          write (uout,*) '# Writing Surface for all ineqs done'
C
C  Di+    The following code writes atomic surface files in different
C         formats for visualization purposes. It only works for
C         regular THETA/PHI grids
C
C         Only code for GRASP format is currently operative
C
C         To implement this code in the current version of PROMOLDEN
C         one additional keyword would be needed to select the printing
C         of atomic surface files for visualization
C
C
          if (.not.allocated(x)) then
            allocate (x(npang),stat=ier)
            if (ier.ne.0) stop 'wrsurf.f: Cannot allocate array x()'
          endif
          if (.not.allocated(y)) then
            allocate (y(npang),stat=ier)
            if (ier.ne.0) stop 'wrsurf.f: Cannot allocate array y()'
          endif
          if (.not.allocated(z)) then
            allocate (z(npang),stat=ier)
            if (ier.ne.0) stop 'wrsurf.f: Cannot allocate array z()'
          endif
          idummy=0
          rdummy=0.05
c
          do i=1,neq
            i1=i/1000
            i2=mod(i,1000)
            i3=i2/100
            i4=mod(i2,100)
            i5=i4/10
            i6=mod(i4,10)
            d4=digs(i1)//digs(i3)//digs(i5)//digs(i6)
c           suffix='.sph'
c           suffix='.ms'
            suffix='.srf'
            fneq=files(1:leng(files))//"-"//d4//suffix
c
c   SPH and DSM formats in ASCII files 
c
c           open (lsu,file=fneq,status='new',iostat=nerr)
c
c  SPH format heading
c           write(lsu,'(''PROMOLDEN Surf spheres'')') 
c           write(lsu,'(''cluster'',I5,
c    &      ''   number of spheres in cluster'',I5)') i,npang
c
c   GRASP format (binary unformatted) handled in graspsurf subroutine
c
            do j=1,npang
               r=rlimsurf(i,j,nlimsurf(i,j))
               x(j)=r*st(j)*cp(j) + xyzrho(i,1)
               y(j)=r*st(j)*sp(j) + xyzrho(i,2)
               z(j)=r*ct(j) + xyzrho(i,3)
               x(j)=x(j)*0.5291772
               y(j)=y(j)*0.5291772
               z(j)=z(j)*0.5291772
C
C  DMS format (probably not useful as it can not be edited by Chimera) 
C
C              write (lsu,'(''MOL '',I4,'' Du '',3f10.3,'' SC0  0.0'')') 
C    &                      i,x(j),y(j),z(j)
C
C  SPH format (too heavy for visualization) 
C
c             write (lsu,'(I5,3f10.5,f8.3,I5,I2,I3)') 
c    &               idummy,x(j),y(j),z(j),rdummy,idummy,idummy,idummy
c
            enddo
c
c   GRASP binary format 
c
            do j=1,3
             mid(j)=xyzrho(i,j)
            enddo
            call graspsurf(npang,npth,npph,fneq,lsu,x,y,z,mid)
c
c           close(lsu)
c
          enddo
          if (allocated(x)) then
            deallocate (x,stat=ier)
            if (ier.ne.0) stop 'wrsurf.f: Cannot deallocate array x()'
          endif
          if (allocated(y)) then
            deallocate (y,stat=ier)
            if (ier.ne.0) stop 'wrsurf.f: Cannot deallocate array y()'
          endif
          if (allocated(z)) then
            deallocate (z,stat=ier)
            if (ier.ne.0) stop 'wrsurf.f: Cannot deallocate array z()'
          endif
C
C  Di+
C

        endif
c
c.......The file with the interatomic surfaces of all the ineqs exists.
c
      else
        allsurfaces=.true.
        read (lsu) ntp1,ntp2
        if (ntp1.ne.neq) then
          call error('wrsurf','surface is not compatible',faterr)
        endif
        if (ntp2.ne.npang) then
          write (uerr,214) npang,ntp2
          call error('wrsurf','surface is not compatible',faterr)
        endif
        do i=1,neq
          read (lsu) (nlimsurf(i,j),j=1,npang)
        enddo
        do i=1,neq
          do j=1,npang
             nsurf=nlimsurf(i,j)
             read (lsu) (rlimsurf(i,j,k),k=1,nsurf)
          enddo
          surfineq(i)=.true.
        enddo
        read (lsu) (rmsurf(i),i=1,neq)
        close (lsu)
        write (uout,*) '# Reading Surface done'
      endif
      return
c
c.....Formats
c
 112  format (1x,'# Done using ',F16.6,' seconds')
 111    format (/,' SURFACE INFO FOR INEQ = ',I4,/,
     &  1X,'NEQ   = ',I5,5X,'NPANG = ',I5,5X,'NR   = ',I5,/,
     &  1X,'IQUDR = ',I5,5X,'IQUDT = ',I5,5x,'MAPR = ',I5,5X,/,
     &  1X,'MAPT  = ',I5,5X,'MMAP  = ',I5,5x,'MPR  = ',I5,5x,
     &     'ALFM  =',F12.6)
 113    format (/,' SURFACE INFO INEQ = ',I4,/,
     &  1X,'NEQ   = ',I5,5X,'NPANG = ',I5,5X,'NPTH = ',I5,5x,
     &     'NPPH = ',I5,5X,'NR   = ',I5,/,
     &  1X,'IQUDR = ',I5,5X,'IQUDT = ',I5,5x,'MAPR = ',I5,5X,/,
     &  1X,'MAPT  = ',I5,5X,'MMAP  = ',I5,5x,'MPR  = ',I5,5x,
     &     'ALFM  =',F12.6)
 114  format (1x,'SURFACE INTERSECTIONS FOR ALL ANGLES = ',/,
     &  1000(20I4,/))
 115  format (1x,
     &  'RLIMSURF () FOR CT,ST,CP,SP =',4(1x,F12.6),' = ',20(1x,F12.6))
 116  format (1x,'RMAXSURF = ',F12.6)
 117  format (' CT,ST,CP,SP = cos(theta),sin(theta),cos(phi),sin(phi)') 

 200  format (1x,'# wrsurf.f: non compatible surface')
 201  format (1x,'# wrsurf.f: non compatible surface',/,
     &        1x,'# npth value in file  = ',I6,/,
     &        1x,'# npth value in input = ',I6)
 202  format (1x,'# wrsurf.f: non compatible surface',/,
     &        1x,'# npph value in file  = ',I6,/,
     &        1x,'# npph value in input = ',I6)
 203  format (1x,'# wrsurf.f: non compatible surface',/,
     &        1x,'# neq value in file  = ',I6,/,
     &        1x,'# neq value in input = ',I6)
 204  format (1x,'# wrsurf.f: non compatible surface',/,
     &        1x,'# npang value in file  = ',I6,/,
     &        1x,'# npang value in input = ',I6)
 205  format (1x,'# wrsurf.f: non compatible surface',/,
     &        1x,'# nr value in file  = ',I6,/,
     &        1x,'# nr value in input = ',I6)
 206  format (1x,'# wrsurf.f: non compatible surface',/,
     &        1x,'# iqudr value in file  = ',I6,/,
     &        1x,'# iqudr value in input = ',I6)
 207  format (1x,'# wrsurf.f: non compatible surface',/,
     &        1x,'# iqudt value in file  = ',I6,/,
     &        1x,'# iqudt value in input = ',I6)
 208  format (1x,'# wrsurf.f: non compatible surface',/,
     &        1x,'# mapr value in file  = ',I6,/,
     &        1x,'# mapr value in input = ',I6)
 209  format (1x,'# wrsurf.f: non compatible surface',/,
     &        1x,'# mapt value in file  = ',I6,/,
     &        1x,'# mapt value in input = ',I6)
 210  format (1x,'# wrsurf.f: non compatible surface',/,
     &        1x,'# mmap value in file  = ',I6,/,
     &        1x,'# mmap value in input = ',I6)
 211  format (1x,'# wrsurf.f: non compatible surface',/,
     &        1x,'# mpr value in file  = ',I6,/,
     &        1x,'# mpr value in input = ',I6)
 212  format (1x,'# wrsurf.f: non compatible surface',/,
     &        1x,'# alfm value in file  = ',F14.6,/,
     &        1x,'# alfm value in input = ',F14.6)
 213  format (1x,'# Faterr !! rmsurf value for ineq = ',I4,/,
     &        1x,'# in file different from the input value',/,
     &        1x,'# Input      = ',F14.6,/,
     &        1x,'# File value = ',F14.6)
 214  format (1x,'# Faterr !! Number of angular points',/,
     &        1x,'# in file different from the input value',/,
     &        1x,'# Input is ',I6,' File value is ',I6)
      end
c
c Di+ The following subroutines print out the GRASP atomic surface file 
c    
c
      subroutine  graspsurf(npang,npth,npph,fneq,lsu,x,y,z,mid)
c
      USE       mod_prec, only: ip, dp 
c
      integer(kind=ip) :: npang,npth,npph,lsu  
      real(kind=dp) :: x(npang), y(npang), z(npang) , mid(3)
      character(len=132) :: fneq
c     
      real(kind=dp),allocatable,dimension (:,:) :: vert,vnor
      integer,allocatable,dimension (:) :: vindx , nnorm 
      integer :: vtot,itot,itriang,ivert,nvert,ntriang,it,i,j,k,l
      integer :: inext1,inext2
      real(kind=dp) :: north(3),south(3),center(3),u(3)
      real(kind=dp) :: raver,rdist,unorm,rdummy
c
      character*80 line(5),fname,pname
c
c     vert ---> Cartesian coordinates of the vertex points
c     vnor  --->  average normal vector at the vertex point
c     vindx  ---> array of pointers to the vertex points constituting a triangle
c     
      ntriang=npph*(2*(npth-1))  +   2 * npph
      nvert=npang + 2
      vtot=nvert
      itot=ntriang
c     
      if (.not.allocated(vert)) then
         allocate (vert(3,vtot),stat=ier)
         if (ier.ne.0) stop 'wrsurf.f: Cannot allocate array vert()'
      endif
      if (.not.allocated(vnor)) then
         allocate (vnor(3,vtot),stat=ier)
         if (ier.ne.0) stop 'wrsurf.f: Cannot allocate array vnor()'
      endif
      if (.not.allocated(nnorm)) then
         allocate (nnorm(vtot),stat=ier)
         if (ier.ne.0) stop 'wrsurf.f: Cannot allocate array vnor()'
      endif
      if (.not.allocated(vindx)) then
         allocate (vindx(3*itot),stat=ier)
         if (ier.ne.0) stop 'wrsurf.f: Cannot allocate array vindx()'
      endif
c
      do i=1,vtot
      nnorm(i)=0
      do j=1,3
          vert(j,i)=0.0
          vnor(j,i)=0.0
      enddo
      enddo
c
      itriang=0
      ivert=0
c
      do i=1,npph 

        do it=1,npth-1
         
          ivert=ivert+1
          vert(1,ivert)=x(ivert)
          vert(2,ivert)=y(ivert)
          vert(3,ivert)=z(ivert)

          if ( i .lt. npph) then
             inext1=ivert+npth
             inext2=ivert+npth+1
          else
             inext1=it
             inext2=it+1
          endif

          itriang=itriang+1
          if ( itriang .gt. itot ) then
             print*,'itriang > itot' 
             stop
          endif
          vindx(3*(itriang-1) + 1 ) = ivert
          vindx(3*(itriang-1) + 2 ) = ivert+1
          vindx(3*(itriang-1) + 3 ) = inext2
          call surfnorm ( x , y, z, vnor, nnorm, nvert, 
     &                    ivert, ivert +1, inext2 )
          
          itriang=itriang+1
          if ( itriang .gt. itot ) then
             print*,'itriang > itot' 
             stop
          endif
          vindx(3*(itriang-1) + 1 ) = inext2
          vindx(3*(itriang-1) + 2 ) = inext1
          vindx(3*(itriang-1) + 3 ) = ivert
          call surfnorm ( x , y, z, vnor, nnorm, nvert, 
     &                    inext2 , inext1, ivert )

        enddo
        ivert=ivert+1
        vert(1,ivert)=x(ivert)
        vert(2,ivert)=y(ivert)
        vert(3,ivert)=z(ivert)
         
      enddo
c
c     North and south pole vertex points
c
c
c     The xyz coordinates of poles are estimated
c     by averaging the rim of closer PHI points 
c
c
      do j=1,3
         north(j)=0.0
         south(j)=0.0
      enddo
      do i=1,npph
         north(1)=north(1)+x(npth*(i-1) + 1 ) 
         north(2)=north(2)+y(npth*(i-1) + 1 ) 
         north(3)=north(3)+z(npth*(i-1) + 1 ) 
         south(1)=south(1)+x(npth*(i-1) + npth ) 
         south(2)=south(2)+y(npth*(i-1) + npth ) 
         south(3)=south(3)+z(npth*(i-1) + npth ) 
      enddo
      do j=1,3
         north(j)=north(j)/float(npph)
         south(j)=south(j)/float(npph) 
         center(j)=mid(j) 
      enddo
c
c     The position of the pole is (imperfectly) refined 
c     by adjusting its relative distance to the center
c     of coordinates
c
      raver=0.0
      do i=1,npph
         rdist=0.0
         do j=1,3
            rdist= rdist + ( vert(j,npth*(i-1)+1) - center(j) )**2
         enddo
         rdist=sqrt(rdist)
         raver=raver+rdist
      enddo
      raver=raver/float(npph)
      unorm=0.0
      do j=1,3
        u(j)=north(j)-center(j)
        unorm=unorm+u(j)**2
      enddo  
      unorm=sqrt(unorm)
      do j=1,3
        u(j)=u(j)/unorm
      enddo
      do j=1,3
        north(j)=center(j)+raver*u(j)
      enddo
c
c     The "refined" pole is added to the set of vertex points
c     and a new set of triangles is added
c
      ivert=ivert+1
      vert(1,ivert)=north(1)
      vert(2,ivert)=north(2)
      vert(3,ivert)=north(3)
      do i=1,npph 
c
        if ( i .lt. npph) then
             inext1=npth*(i-1) + 1 
             inext2=npth*i + 1
        else
             inext1=npth*(i-1) + 1 
             inext2=1
        endif
c
        itriang=itriang+1
        if ( itriang .gt. itot ) then
             print*,'itriang > itot' 
             stop
        endif
        vindx(3*(itriang-1) + 1 ) = ivert 
        vindx(3*(itriang-1) + 2 ) = inext1
        vindx(3*(itriang-1) + 3 ) = inext2
        call surfnorm ( x , y, z, vnor, nnorm, nvert, 
     &                ivert, inext1 , inext2 )
      enddo
c
c    Similar things with the south pole .....
c
      raver=0.0
      do i=1,npph
         rdist=0.0
         do j=1,3
            rdist= rdist + ( vert(j,npth*(i-1)+npth) - center(j) )**2
         enddo
         rdist=sqrt(rdist)
         raver=raver+rdist
      enddo
      raver=raver/float(npph)
      unorm=0.0
      do j=1,3
        u(j)=south(j)-center(j)
        unorm=unorm+u(j)**2
      enddo  
      unorm=sqrt(unorm)
      do j=1,3
        u(j)=u(j)/unorm
      enddo
      do j=1,3
        south(j)=center(j)+raver*u(j)
      enddo
      ivert=ivert+1
      vert(1,ivert)=south(1)
      vert(2,ivert)=south(2)
      vert(3,ivert)=south(3)
      do i=1,npph 
c
        if ( i .lt. npph) then
             inext1=npth*(i-1) + npth
             inext2=npth*i + npth
        else
             inext1=npth*(i-1) + npth 
             inext2=npth
        endif
c
        itriang=itriang+1
        if ( itriang .gt. itot ) then
             print*,'itriang > itot' 
             stop
        endif
        vindx(3*(itriang-1) + 1 ) = inext1
        vindx(3*(itriang-1) + 2 ) = ivert  
        vindx(3*(itriang-1) + 3 ) = inext2
        call surfnorm ( x , y, z, vnor, nnorm, nvert, 
     &                inext1, ivert, inext2 )
      enddo
c
      do i=1,nvert
         rdummy=0.0
         do j=1,3
            vnor(j,i)=vnor(j,i)/float(nnorm(i))
            rdummy=rdummy+vnor(j,i)**2
         enddo
         rdummy=sqrt(rdummy)
         do j=1,3
            vnor(j,i)=vnor(j,i)/rdummy 
         enddo
      enddo
c
c     Print out everything to the binary output file
c
      if ( itot .ne. itriang) then
         write(6,*) 'Problem in # of triangle faces'
         stop
      endif
      if ( vtot .ne. ivert) then
         write(6,*) 'Problem in # of vertex points'
         stop
      endif
c
      do i=1,5
          line(i)=" "
      end do
C
      line(1)="format=2"
      line(2)="vertices,normals,triangles"
      line(3)=" "
      igrid=65
      scale=0.33 
      write(line(4),'(3i6,f12.6)') vtot,itot,igrid,scale
      write(line(5),'(3f12.6)') sngl(mid)
c     
      open(lsu,file=fneq,form="unformatted")
c     
      do i=1,5
          write(lsu) line(i)
      enddo
c     
      write(lsu) sngl(vert)
      write(lsu) sngl(vnor)
      write(lsu) vindx
c     
      close(lsu)
c
      if (allocated(vert)) then
         deallocate (vert,stat=ier)
         if (ier.ne.0) stop 'wrsurf.f: Cannot deallocate array vert()'
      endif
      if (allocated(vnor)) then
         deallocate (vnor,stat=ier)
         if (ier.ne.0) stop 'wrsurf.f: Cannot deallocate array vnor()'
      endif
      if (allocated(vindx)) then
         deallocate (vindx,stat=ier)
         if (ier.ne.0) stop 'wrsurf.f: Cannot deallocate array vindx()'
      endif
      if (allocated(nnorm)) then
         deallocate (nnorm,stat=ier)
         if (ier.ne.0) stop 'wrsurf.f: Cannot deallocate array nnorm()'
      endif
      return
      end
C
      subroutine surfnorm ( x , y, z, vnor, nnorm, nvert, i1, i2, i3) 
c
      USE       mod_prec, only: ip, dp 
c
      integer(kind=ip) ::  nvert, i1 i2, i3
      real(kind=dp) :: x(nvert), y(nvert), z(nvert) 

      real(kind=dp) :: vnor(3,nvert) 
      integer(kind=ip) ::  nnorm(nvert) 
c
      real(kind=dp) :: a(3),b(3),c(3) , cnorm
c
      a(1)=x(i2)-x(i1)
      a(2)=y(i2)-y(i1)
      a(3)=z(i2)-z(i1)
      b(1)=x(i3)-x(i1)
      b(2)=y(i3)-y(i1)
      b(3)=z(i3)-z(i1)
      call vecprod(a,b,c)
      cnorm=0.0
      do j=1,3
        cnorm=cnorm+c(j)**2
      enddo
      cnorm=sqrt(cnorm)
      do j=1,3
        c(j)=c(j)/cnorm 
      enddo
      nnorm(i1)=nnorm(i1)+1
      do j=1,3
         vnor(j,i1)=vnor(j,i1)+c(j)
      enddo
c
      b(1)=x(i3)-x(i2)
      b(2)=y(i3)-y(i2)
      b(3)=z(i3)-z(i2)
      call vecprod(a,b,c)
      cnorm=0.0
      do j=1,3
        cnorm=cnorm+c(j)**2
      enddo
      cnorm=sqrt(cnorm)
      do j=1,3
        c(j)=c(j)/cnorm 
      enddo
      nnorm(i2)=nnorm(i2)+1
      do j=1,3
         vnor(j,i2)=vnor(j,i2)+c(j)
      enddo
c
      a(1)=x(i1)-x(i3)
      a(2)=y(i1)-y(i3)
      a(3)=z(i1)-z(i3)
      call vecprod(b,a,c)
      cnorm=0.0
      do j=1,3
        cnorm=cnorm+c(j)**2
      enddo
      cnorm=sqrt(cnorm)
      do j=1,3
        c(j)=c(j)/cnorm 
      enddo
      nnorm(i3)=nnorm(i3)+1
      do j=1,3
         vnor(j,i3)=vnor(j,i3)+c(j)
      enddo
c
      return
      end

