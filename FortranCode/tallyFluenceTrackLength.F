!*******************************************************************
!*                          TALLY                                  *
!*                     FLUENCE SPECTRUM                            *
!*                  (TRACK LENGTH ESTIMATOR)                       *
!*                                                                 *
!* Short description:                                              *
!*   Tallies track lengths in order to estimate the fluence spectra*
!*   of each type of particle in a given material.                 *
!* Rev:                                                            *
!*   2015-01-07                                                    *
!*******************************************************************


!*******************************************************************
!*******************************************************************
!*    Modules section                                              *
!*******************************************************************
!*******************************************************************

      module FTLmod
!*******************************************************************
!*    Vars for the FTL tally.                                      *
!*******************************************************************
      implicit none
      save
      logical active,isLinScale
      integer nbin,matdet
      integer, parameter :: nbinmax=32000
      real*8 flutmp(4,nbinmax),flu(4,nbinmax),flu2(4,nbinmax)
      real*8 egrid(nbinmax),ebingrd(nbinmax)
      real*8 emin,ebin,iebin,eratio,unclimit
      end


!*******************************************************************
!*******************************************************************
!*    Routines section                                             *
!*******************************************************************
!*******************************************************************

      subroutine FTLtally(mode,arg)
!*******************************************************************
!*    Input:                                                       *
!*      mode: identifies the state of the calling procedure        *
!*      arg: different meanings for different modes; see README    *
!*******************************************************************
      use FTLmod
      use PENELOPE_mod, only : e0step,ssoft
      use TRACK_mod, only : e,wght,mat,kpar
      implicit none
      integer mode
      real*8 arg

      logical, save :: inside          ! Yes if inside det mat
      real*8, save :: eini,de          ! Energy at start of flight & soft E loss in substep
      integer*4 i,k

      if (.not.active) return

      if (mode.eq.3) then              ! Called STEP and flew without interface crossings
        if (.not.inside) return        ! Flight outside detection material, nothing to do
        de = ssoft*arg                 ! Eloss for charged particle substep; ignored for photons (DE=0)
        call trackl(arg*wght)          ! Score track length ARG
        return
      endif

      if (mode.eq.4) then              ! A STEP with a crossing has occurred
        if (inside) then               ! Inside det material, do the same as for mode.eq.3
          de = ssoft*arg               ! Eloss for charged particle substep; ignored for photons (DE=0)
          call trackl(arg*wght)        ! Score track length ARG
          inside = .false.             ! If we were inside det material, now we must be outside
        else if (mat.eq.matdet) then   ! If we were not outside det material, find out if we are now
          eini = e                     ! Set energy at start of flight to current value
          inside = .true.
        endif
        return
      endif

      if (mode.lt.0.and.mode.gt.-10) then  ! A KNOCK interaction has occurred
        if (.not.inside) return            ! Outside detection material, nothing to do
        if (kpar.ne.2.and.mode.eq.-1) then ! Soft collision of charged particle
          eini = e0step                    ! E0STEP set by previous call to JUMP
        else
          eini = e                         ! Energy at start of flight for all other cases
        endif
        return
      endif

      if (mode.eq.-99) then         ! Simulation of a new particle begins
        if (mat.eq.matdet) then     ! Inside detector
          eini = e                  ! Energy at start of flight
          inside = .true.           ! Set flag
        else                        ! Outside detector
          inside = .false.
        endif
        return
      endif

      if (mode.eq.6) then           ! End-of-History, dump and clear temp counter
        do i=1,nbin
#ifdef PROTONS
          do k=1,4
#else
          do k=1,3
#endif
            if (flutmp(k,i).eq.0.0) cycle  ! Skip void counters
            flu(k,i)    = flu(k,i) +flutmp(k,i)
            flu2(k,i)   = flu2(k,i)+flutmp(k,i)**2
            flutmp(k,i) = 0.0
          enddo
        enddo
        return
      endif


      contains
        subroutine trackl(ds)
        !*******************************************************************
        !*    Tallies track length distributing its value in E bins.       *
        !*                                                                 *
        !*    Input:                                                       *
        !*      ds: track length (cm) times statistical weight             *
        !*******************************************************************
        use PENELOPE_mod, only : eabs
        use CSIMPHcommonMod, only : eabsh  ! Protons
        implicit none
        real*8 ds

        integer*4 i,binhi,binlo
        real*8 efin,dflu,dsde
        real*8, parameter :: zero=1.0e-30  ! Zero E loss

        if (isLinScale) then                           ! Energy grid in LINEAR SCALE

          binhi = (eini-emin)*iebin+1.0                ! Highest E bin that scores
          if (kpar.eq.2.or.de.lt.+zero) then           ! No need to spread fluence over various bins
            if (binhi.lt.1.or.binhi.gt.nbin) return    ! Not in ROI
            flutmp(kpar,binhi) = flutmp(kpar,binhi)+ds ! Score track length
            return                                     ! Done
          endif

#ifdef PROTONS
          if (kpar.eq.4) then                          ! Its a proton
            efin  = max(eini-de,eabsh(matdet))         ! Final E after flight, cannot be below Eabs(p+)
          else                                         ! Not a proton
#endif
            efin  = max(eini-de,eabs(kpar,matdet))     ! Final E after flight, cannot be below Eabs
#ifdef PROTONS
          endif
#endif
          binlo = (efin-emin)*iebin+1.0                ! Lowest E bin that scores

          dsde  = ds/de                                ! Inverse weighted average stopping power = WGHT/SSOFT
          dflu  = ebin*dsde                            ! Quantity to be scored in each E bin
          do i=max(binlo,1),min(binhi,nbin)            ! Spread over all bins from Efin to Eini, skipping invalid
            flutmp(kpar,i) = flutmp(kpar,i)+dflu       ! Score
          enddo
          if (binhi.gt.0.and.binhi.le.nbin)            ! Subtract excess in Eini (highest) bin
     &      flutmp(kpar,binhi) = flutmp(kpar,binhi)-
     &                           (egrid(binhi)+ebin-eini)*dsde
          if (binlo.gt.0.and.binlo.le.nbin)            ! Subtract excess in Efin (lowest) bin
     &      flutmp(kpar,binlo) = flutmp(kpar,binlo)-
     &                           (efin-egrid(binlo))*dsde

        else      ! Energy grid in LOG SCALE (see detailed comments in code for linear scale)

          binhi = (log(eini)-emin)*iebin+1.0           ! Highest E bin, log scale
          if (kpar.eq.2.or.de.lt.+zero) then           ! Same code as linear grid
            if (binhi.lt.1.or.binhi.gt.nbin) return
            flutmp(kpar,binhi) = flutmp(kpar,binhi)+ds
            return
          endif

#ifdef PROTONS
          if (kpar.eq.4) then                          ! Same code as linear grid
            efin  = max(eini-de,eabsh(matdet))
          else
#endif
            efin  = max(eini-de,eabs(kpar,matdet))
#ifdef PROTONS
          endif
#endif
          binlo = (log(efin)-emin)*iebin+1.0           ! Lowest E bin, log scale

          dsde  = ds/de
          do i=max(binlo,1),min(binhi,nbin)
            flutmp(kpar,i) = flutmp(kpar,i)+ebingrd(i)*dsde  ! Quantity to be scored in bin, log scale
          enddo
          if (binhi.gt.0.and.binhi.le.nbin)
     &      flutmp(kpar,binhi) = flutmp(kpar,binhi)-
     &                           (egrid(binhi)*eratio-eini)*dsde
          if (binlo.gt.0.and.binlo.le.nbin)
     &      flutmp(kpar,binlo) = flutmp(kpar,binlo)-
     &                           (efin-egrid(binlo))*dsde

        endif
        end subroutine
      end subroutine


      subroutine FTLreport(mode,n,cputim,uncdone)
!*******************************************************************
!*    Input:                                                       *
!*      mode: -1 if called at end-of-simulation, 0 for dump only,  *
!*             1 for tally report only, 2 for both.                *
!*      n: no. of histories simulated                              *
!*      cputim: elapsed CPU time                                   *
!*    Output:                                                      *
!*      uncdone: 2 if uncert reached, 1 if not defined, 0 else     *
!*               not computed if mode=0                            *
!*               The reported status is for kpar=1 (electrons)     *
!*******************************************************************
      use FTLmod
      use RSEEDcommonMod
      implicit none
      integer mode,uncdone
      real*8 n,cputim

      character*80 buffer
      integer out,finduf,errorcode,ikpar
      integer*4 i,nchan
      real*8 q,sigma,eff,avesig,maxq,fact,emiddle,uncert,invn

      uncdone = 1
      if (.not.active) return

      if (mode.ne.1) call FTLdump(1)  ! Sim dump file
      if (mode.eq.0) return           ! Dump only, no report

#ifdef PROTONS
      do ikpar=4,1,-1                 ! One report for each particle type
#else
      do ikpar=3,1,-1                 ! One report for each particle type
#endif
        ! Prepare output files:
        select case (ikpar)
        case(4)
          buffer = 'tallyFluenceTrackLength-proton.dat'
        case(3)
          buffer = 'tallyFluenceTrackLength-positron.dat'
        case(2)
          buffer = 'tallyFluenceTrackLength-photon.dat'
        case(1)
          buffer = 'tallyFluenceTrackLength-electron.dat'
        case default
          buffer = 'ERRORkpar'
        end select
        out = finduf()
        open(out,file=buffer,iostat=errorcode)
        if (errorcode.ne.0) then
          write(*,*) ''
          write(*,'(a)')
     &      '*********************************************'
          write(*,'(a)')
     &      'FTLreport:ERROR: cannot open output data file'
          write(*,'(a)')
     &      '*********************************************'
          close(out)  ! Just in case
          return
        endif

        ! Evaluate 1/2 of the max fluence score:
        maxq = 0.0
        do i=1,nbin
          q = flu(ikpar,i)/ebingrd(i)
          maxq = max(maxq,q)
        enddo
        maxq = 0.5*maxq

        write(out,'(a)')
     &'#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
        write(out,'(a)') '# [SECTION REPORT FLUENCE SPECTRUM]'
        if (ikpar.eq.4) buffer = '# PROTONS'
        if (ikpar.eq.3) buffer = '# POSITRONS'
        if (ikpar.eq.2) buffer = '# PHOTONS'
        if (ikpar.eq.1) buffer = '# ELECTRONS'
        write(out,'(a)') buffer
        write(out,'(a)') '#'
        write(out,'(a)')
     &    '# The spectral fluence is integrated over the '//
     &    'detector volume, its units are cm^3/(cm^2*eV)=cm/eV.'
        write(out,'(a)') '#'
        write(out,'(a)') '# Detection material:'
        write(out,'(a,i0)') '#  ',matdet
        write(out,'(a)') '# Number of energy bins:'
        write(out,'(a,i0)') '#  ',nbin
        if (isLinScale) then
          write(out,'(a)') '# Linear energy scale, bin width (eV):'
        else
          write(out,'(a)') '# Log energy scale, E_{i+1}/E_i:'
        endif
        write(out,'(a,es12.5)') '#  ',exp(ebin)
        write(out,'(a)') '#'
        write(out,'(a)')
     &   '# For plotting purposes, two energies per bin are given, '//
     &   'namely, the low end and the middle point of each bin.'
        write(out,'(a)') '#'
        write(out,'(a)')
     &   '# Elow(eV) : Emiddle(eV) : fluence(cm/eV) : +-2sigma : nbin'

        avesig = 0.0
        nchan = 0
        invn = 1.0/n
        do i=1,nbin
          q = flu(ikpar,i)*invn
          sigma = (flu2(ikpar,i)*invn-q**2)*invn
          sigma = sqrt(max(sigma,0.0))
          fact = 1.0/ebingrd(i)
          q = q*fact
          sigma = sigma*fact
          if (isLinScale) then
            emiddle = egrid(i)+ebin*0.5
          else
            emiddle = egrid(i)*(eratio+1.0)*0.5  ! Middle energy of the bin
          endif
          write(out,'(3(1x,es12.5),1x,es8.1,1x,i0)')
     &      egrid(i),emiddle,q,2.0*sigma,i
          ! Evaluate average uncertainty for scores above 1/2 max score:
          if (flu(ikpar,i)*fact.gt.maxq) then
            avesig = avesig+(sigma/q)**2
            nchan = nchan+1
          endif
        enddo

        uncdone = 0
        if (nchan.gt.0) then
          uncert = 200.0*sqrt(avesig/nchan)
          if (uncert.lt.unclimit) uncdone = 2  ! Uncertainty reached
        else
          uncert = 0.0  ! Uncertainty assumed not reached when score is nil
        endif

        ! Generic report:
        write(out,'(a)') ' '
        write(out,'(a)') '# Performance report'
        write(out,'(a)') '#   Random seeds:'
        write(out,'(a,i10)') '#   ',seed1
        write(out,'(a,i10)') '#   ',seed2
        write(out,'(a)') '#   No. of histories simulated [N]:'
        write(out,'(a,f18.0)') '#   ',n
        write(out,'(a)') '#   CPU time [t] (s):'
        write(out,'(a,es12.5)') '#   ',cputim
        if (cputim.gt.0.0) then
          write(out,'(a)') '#   Speed (histories/s):'
          write(out,'(a,es12.5)') '#   ',n/cputim
        endif
        write(out,'(a)')
     & '#   Average uncertainty (above 1/2 max score) in % [uncert]:'
        write(out,'(a,es12.5)') '#   ',uncert
        eff = n*uncert**2
        if (eff.gt.0.0) then
          write(out,'(a)') '#   Intrinsic efficiency [N*uncert^2]^-1:'
          write(out,'(a,es12.5)') '#   ',1.0/eff
        endif
        eff = cputim*uncert**2
        if (eff.gt.0.0) then
          write(out,'(a)') '#   Absolute efficiency [t*uncert^2]^-1:'
          write(out,'(a,es12.5)') '#   ',1.0/eff
        endif
        close(out)
      enddo
      end


      subroutine FTLinitally(mem)
!*******************************************************************
!*    Initializes. To be called before TALLY.                      *
!*                                                                 *
!*    Output:                                                      *
!*      mem: Memory used by tally arrays                           *
!*    Comments:                                                    *
!*      - Must be called *AFTER* PEINIT                            *
!*******************************************************************
      use FTLmod
      use dataTypesMod
      implicit none
      real*8 mem

      character*(*) secid,eos
      parameter (secid=
     &  '[SECTION TALLY FLUENCE TRACK LENGTH v.2012-06-01]')
      parameter (eos='[END OF FTL SECTION]')
      character*80 buffer,upcase
      integer i,errorcode
      real*8 emax

      write(*,*) ' '
      write(*,'(a)')
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      call getline(buffer,0)
      if (index(buffer,secid).eq.0) then
        write(*,'(a)') 'FTLinitally:ERROR: incorrect section header;'
        write(*,'(a,a)') '  expecting to find: ',secid
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif
      write(*,'(a)') secid

      ! Report memory usage:
      mem = sizeOfReal8*(3*4+2)*real(nbinmax)
      write(*,'(a)') 'Memory used by static arrays (MB):'
      write(*,'(1x,es12.5)') mem*1.0e-6

      read(*,'(a80)') buffer
      buffer = adjustl(upcase(buffer))
      buffer = buffer(1:scan(buffer,' ')) ! Clip at 1st blank
      if (buffer.eq.'ON') then
        active = .true.
      else if (buffer.eq.'OFF') then
        active = .false.
        write(*, '(a)')
     &    '>>>> Tally Fluence Track Length is OFF >>>>'
        do
          read(*,'(a80)',iostat=errorcode) buffer
          if (errorcode.ne.0) then
            write(*,'(a,a,a)') 'FTLinitally:ERROR: ',
     &       'Unable to find End-Of-Section mark: ',eos
            stop
          endif
          if (index(buffer,eos).ne.0) return
        enddo
      else
        write(*,'(a)')
     &    'FTLinitally:ERROR: expecting to find ON or OFF'
        write(*,'(a)') 'found instead:'
        write(*,'(a)') buffer
        stop
      endif

      ! Read data:
      write(*,'(a)') 'Detection material set to:'
      read(*,*) matdet
      write(*,'(1x,i0)') matdet

      write(*,'(a)')
     &  'Spectrum limits [Emin,Emax] (eV) and no. bins:'
      read(*,'(a80)') buffer
      read(buffer,*) emin,emax,nbin,buffer
      if (nbin.lt.1.or.nbin.gt.nbinmax-1) then
        write(*,'(a)') 'FTLinitally:ERROR: Too many bins.'
        write(*,'(a,i0)') '  Max no. bins is ',nbinmax-1
        stop
      endif
      buffer = adjustl(upcase(buffer))
      buffer = buffer(1:scan(buffer,' ')) ! Clip at 1st blank
      isLinScale = .true.
      if (buffer.eq.'LOG') isLinScale = .false.
      write(*,'(2(1x,es12.5),1x,i0,$)') emin,emax,nbin

      if (isLinScale) then
        write(*,'(a)') ' linear scale'
        eratio = 0.0                       ! Dummy value, not used, just in case
        ebin = (emax-emin)/nbin
        iebin = 1.0/ebin
        nbin = nbin+1                      ! Add one more bin to tally counts with E=Emax
        do i=1,nbin                        ! Set energy grid
          egrid(i) = emin+(i-1)*ebin
          ebingrd(i) = ebin
        enddo
      else
        write(*,'(a)') ' log scale'
        eratio = (emax/emin)**(1.0d0/nbin) ! Ratio of one grid point to the previous
        ebin = log(eratio)                 ! Bin width in log scale
        iebin = 1.0/ebin
        nbin = nbin+1                      ! Add one more bin to tally counts with E=Emax
        emin = log(emin)                   ! Store values is log scale
        do i=1,nbin                        ! Set energy grid
          egrid(i) = exp(emin+(i-1)*ebin)
          ebingrd(i) = egrid(i)*(eratio-1.0)
        enddo
      endif

      write(*,'(a)') 'Relative uncertainty (%) requested:'
      read(*,*) unclimit
      write(*,'(1x,es12.5)') unclimit

      ! Clear counters:
      flu    = 0.0
      flu2   = 0.0
      flutmp = 0.0

      read(*,'(a80)') buffer
      if (index(buffer,eos).eq.0) then
        write(*,*) 'FTLinitally:ERROR: End-Of-Section mark not found'
        write(*,'(a,a)') '  expecting to find: ',eos
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif

      call FTLdump(0)  ! Sim restart file

      write(*,'(a)') '>>>> FTL tally initialization finished >>>>'
      end


      subroutine FTLdump(mode)
!*******************************************************************
!*    Dumps into or reads data from a dump file.                   *
!*                                                                 *
!*    Input:                                                       *
!*      mode -> 1 to write dump file, else to read from it.        *
!*******************************************************************
      use FTLmod, only : flu,flu2,nbin
      use dumpMod
      implicit none
      integer mode

      if (mode.eq.1) then
        if (dumpfile.lt.0) return  ! No dump file open
        write(dumpfile) flu(1:4,1:nbin),flu2(1:4,1:nbin)
      else
        if (restartfile.lt.0) return  ! No restart file open
        read(restartfile) flu(1:4,1:nbin),flu2(1:4,1:nbin)
      endif
      end


!>>>> End Of File >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

