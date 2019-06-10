    else if (equal(word,'surf')) then
      neqsurf = 0_ip
      neqs: do
        ok = isinteger(icsurf,line,lp)
        if (ok) then
          ok = ok .and. icsurf.ne.0_ip
          ok = ok .and. icsurf.le.ncent
          icsurf = abs(icsurf)
          if (ok) then
            do i = 1,neqsurf
              if (icsurf.eq.insurf(i)) cycle neqs
            end do
            neqsurf = neqsurf + 1_ip
            insurf(neqsurf) = icsurf
          else
            call ferror('dosurf', 'wrong surf line', faterr) 
          endif
        else
          exit neqs
        end if
      end do neqs

    else if (equal(word,'ntrial')) then
      ok = isinteger(ntrial, line, lp)
      ok = ok .and. ntrial.ne.0_ip
      if (.not.ok) call ferror('dosurf', 'wrong ntrial line', faterr) 
      ntrial = abs(ntrial)
      if (mod(ntrial,2).eq.0) ntrial = ntrial + 1_ip
      write (uout,'(1x,a,1x,i0)') string('# *** Variable ntrial changed to :'), ntrial

    else if (equal(word,'rprimer')) then
      ok = isreal(rprimer, line, lp)
      if (.not.ok) call ferror('dosurf', 'wrong rprimer line', faterr) 
      rprimer = abs(rprimer)
      write (uout,'(1x,a,1x,e13.6)') string('# *** Variable rprimer changed to :'), rprimer

    else if (equal(word,'rsearch')) then
      ok = isinteger(nsearch,line,lp)
      ok = ok .and. isinteger(nrstart,line,lp)
      ok = ok .and. nsearch.ne.0_ip
      ok = ok .and. nrstart.ne.0_ip
      if (ok) then
        if (nrstart.gt.maxstart) then
          call ferror('dosurf', 'nrstart.gt.maxstart in rsearch order', faterr)
        end if
        nrstart = abs(nrstart)
        if (nsearch.gt.ncent) then
          call ferror('dosurf', 'nsearch.gt.ncent in rsearch order', faterr)
        end if
        if (nsearch.ne.-1_ip) then
          nsearch = abs(nsearch)
          nrsearch(nsearch) = nrstart
          lstart(nsearch) = .true.
          icon = 0_ip
          do while (icon.lt.nrstart)
            if (.not.isreal(rini,line,lp)) then
              call ferror('dosurf', 'bad format in rsearch order', faterr)
            end if
            icon = icon + 1_ip
            rstart(nsearch,icon) = rini
          end do
        else
          forall (ic=1:ncent) nrsearch(ic) = nrstart
          forall (ic=1:ncent) lstart(ic) = .true.
          icon = 0_ip
          do while (icon.lt.nrstart)
            if (.not.isreal(rini,line,lp)) then
              call ferror('dosurf', 'bad format in rsearch order', faterr)
            end if
            icon = icon + 1_ip
            forall (ic=1:ncent) rstart(ic,icon) = rini
          end do
        end if
      end if

