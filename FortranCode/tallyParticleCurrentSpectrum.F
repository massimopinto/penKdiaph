!*******************************************************************
!*                          TALLY                                  *
!*                    PARTICLE CURRENT SPECTRUM                    *
!*                                                                 *
!* Short description:                                              *
!*   Determines the spectrum and total number of particles (grouped*
!*   according to particle type) that entered or were created by   *
!*   the source in a given material.                               *
!* Rev:                                                            *
!*   2015-01-07                                                    *
!*******************************************************************


!*******************************************************************
!*******************************************************************
!*    Modules section                                              *
!*******************************************************************
!*******************************************************************

      module PCSmod
!*******************************************************************
!*    Vars for the PCS tally.                                      *
!*******************************************************************
      implicit none
      save
      logical active
      integer matdet,nbin
      integer, parameter :: nbinmax=32000
      real*8 countmp(nbinmax,4),counts(nbinmax,4),count2(nbinmax,4)
      real*8 tcount2(4),edeptmp(4),edep(4),edep2(4)
      real*8 tedep2,emin,ebin,iebin,unclimit
      end


!*******************************************************************
!*******************************************************************
!*    Routines section                                             *
!*******************************************************************
!*******************************************************************

      subroutine PCStally(mode)
!*******************************************************************
!*    Input:                                                       *
!*      mode: Identifies the state of the calling procedure        *
!*******************************************************************
      use PCSmod
      use TRACK_mod
      implicit none
      integer mode

      integer bin,k,i
      real*8 totedep,totcount

      if (.not.active) return

      if (mode.eq.4.or.mode.eq.7) then  ! New material OR a primary recovered from stack
        if (mat.ne.matdet) return       ! Not in detector
        edeptmp(kpar) = edeptmp(kpar)+e*wght
        bin = (e-emin)*iebin+1.0
        if (bin.lt.1.or.bin.gt.nbin) return  ! Not in ROI
        countmp(bin,kpar) = countmp(bin,kpar)+wght

      else if (mode.eq.6) then          ! End-of-history bookkeeping

        totedep = 0.0
#ifdef PROTONS
        do k=1,4
#else
        do k=1,3
#endif
          totcount = 0.0
          do i=1,nbin
            if (countmp(i,k).eq.0.0) cycle  ! Skip empty bins
            totcount     = totcount   +countmp(i,k)
            counts(i,k)  = counts(i,k)+countmp(i,k)
            count2(i,k)  = count2(i,k)+countmp(i,k)**2
            countmp(i,k) = 0.0
          enddo
          tcount2(k) = tcount2(k)+totcount**2

          totedep    = totedep +edeptmp(k)
          edep(k)    = edep(k) +edeptmp(k)
          edep2(k)   = edep2(k)+edeptmp(k)**2
          edeptmp(k) = 0.0
        enddo
        tedep2 = tedep2+totedep**2

      endif
      end


      subroutine PCSreport(mode,n,cputim,uncdone)
!*******************************************************************
!*    Input:                                                       *
!*      mode:   -1 if called at end-of-simulation, 0 for dump only,*
!*               1 for tally report only, 2 for both.              *
!*      n:   no. of histories simulated                            *
!*      cputim:   elapsed CPU time                                 *
!*    Output:                                                      *
!*      uncdone:   2 if uncert reached, 1 if not defined, 0 else   *
!*                 not computed if mode=0                          *
!*******************************************************************
      use PCSmod
      use RSEEDcommonMod
      implicit none
      integer mode,uncdone
      real*8 n,cputim

      character*80 buffer
      integer kpar,nchan,out,finduf,errorcode,i
      real*8 q,q2,sigma,eff,avesig,maxq,emiddle,elow,uncert,invn

      uncdone = 1
      if (.not.active) return

      if (mode.ne.1) call PCSdump(1)  ! Sim dump file
      if (mode.eq.0) return           ! Dump only, no report

      !*********************************
      !*** Particle Current Spectrum ***
      !*********************************

#ifdef PROTONS
      do kpar=4,1,-1                 ! One report for each particle type
#else
      do kpar=3,1,-1                 ! One report for each particle type
#endif
        ! Prepare output files:
        buffer = 'ERRORname'
        if (kpar.eq.4)
     &  buffer = 'tallyParticleCurrentSpectrum-proton.dat'
        if (kpar.eq.3)
     &  buffer = 'tallyParticleCurrentSpectrum-positron.dat'
        if (kpar.eq.2)
     &  buffer = 'tallyParticleCurrentSpectrum-photon.dat'
        if (kpar.eq.1)
     &  buffer = 'tallyParticleCurrentSpectrum-electron.dat'
        out = finduf()
        open(out,file=buffer,iostat=errorcode)
        if (errorcode.ne.0) then
          write(*,*)
          write(*,'(a)')
     &      '*********************************************'
          write(*,'(a)')
     &      'PCSreport:ERROR: cannot open output data file'
          write(*,'(a)')
     &      '*********************************************'
          close(out)  ! Just in case
          return
        endif

        ! Evaluate 1/2 of the max fluence score:
        maxq   = 0.0
        do i=1,nbin
          maxq = max(maxq,counts(i,kpar))
        enddo
        maxq = 0.5*maxq

        write(out,'(a)')
     &'#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
        write(out,'(a)')
     &    '# [SECTION REPORT PARTICLE CURRENT SPECTRUM]'
        if (kpar.eq.4) buffer = '# PROTONS'
        if (kpar.eq.3) buffer = '# POSITRONS'
        if (kpar.eq.2) buffer = '# PHOTONS'
        if (kpar.eq.1) buffer = '# ELECTRONS'
        write(out,'(a)') buffer
        write(out,'(a)') '#'
        write(out,'(a)') '# Number of particles, '//
     &    'per unit energy interval and per history, '//
     &    'entering the detector.'
        write(out,'(a)') '# Units are 1/eV per history.'
        write(out,'(a)') '#'
        write(out,'(a)') '# Detection material:'
        write(out,'(a,i0)') '#  ',matdet
        write(out,'(a)') '# No. of energy bins and bin width (eV):'
        write(out,'(a,i0,1x,es12.5)') '#  ',nbin,ebin
        write(out,'(a)') '#'
        write(out,'(a)') '# For plotting purposes, '//
     &    'two values per energy bin are given, '//
     &    'namely, the low end and the middle point of each bin.'
        write(out,'(a)') '#'
        write(out,'(a)')
     & '# Elow(eV) : Emiddle(eV) : counts(1/eV/hist) : +-2sigma : nbin'

        avesig = 0.0
        nchan  = 0
        invn = 1.0/n
        do i=1,nbin
          q = counts(i,kpar)*invn
          q2 = count2(i,kpar)*invn
          sigma = (q2-q**2)*invn
          sigma = sqrt(max(sigma,0.0))*iebin
          q = q*iebin
          elow = emin+ebin*(i-1)
          emiddle = elow+ebin*0.5
          write(out,'(3(1x,es12.5),1x,es8.1,1x,i0)')
     &      elow,emiddle,q,2.0*sigma,i
          ! Evaluate average uncertainty for scores above 1/2 max score:
          if (counts(i,kpar).gt.maxq) then
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

        ! Total counter:
        q = sum(counts(1:nbin,kpar))*invn
        q2 = tcount2(kpar)*invn
        sigma = (q2-q**2)*invn
        sigma = sqrt(max(sigma,0.0))
        write(out,'(a)') '#'
        write(out,'(a)')
     &  '# Total particle current per history : +-2sigma:'
        write(out,'(a,es12.5,1x,es8.1)') '# ',q,2.0*sigma

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

      !**********************
      !*** Energy current ***
      !**********************

      ! Prepare output files:
      out = finduf()
      open(out,file='tallyParticleCurrentSpectrum-totalE.dat',
     &     iostat=errorcode)
      if (errorcode.ne.0) then
        write(*,*)
        write(*,'(a)')
     &    '*********************************************'
        write(*,'(a)')
     &    'PCSreport:ERROR: cannot open output data file'
        write(*,'(a)')
     &    '*********************************************'
        close(out)  ! Just in case
        return
      endif

      write(out,'(a)')
     &'#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      write(out,'(a)') '# [SECTION REPORT ENERGY CURRENT]'
      write(out,'(a)') '# Sum of energies, per primary history, '//
     &  'of particles entering the detector.'
      write(out,'(a)') '# Units are eV per history.'
      write(out,'(a)') '#'
      write(out,'(a,i3)') '# Detection material: ',matdet
      write(out,'(a)') '#'
#ifdef PROTONS
      write(out,'(a)')
     &  '# ParticleType=0,1,2,3,4 indicate TotalSum, '//
     &  'electrons, photons, positrons and protons, respectively.'
#else
      write(out,'(a)')
     &  '# ParticleType=0,1,2,3 indicate TotalSum, '//
     &  'electrons, photons and positrons, respectively.'
#endif
      write(out,'(a)') '#'
      write(out,'(a)')
     &  '# ParticleType : EnergyCurrent(eV/hist) : +-2sigma'

      ! Total counter:
      invn = 1.0/n
      q  = sum(edep)*invn
      q2 = tedep2*invn
      sigma = (q2-q**2)*invn
      sigma = sqrt(max(sigma,0.0))
      write(out,'(i2,1x,es12.5,1x,es8.1)') 0,q,2.0*sigma
      ! Evaluate rel. uncertainty:
      uncert = 200.0
      if (q.gt.0.0) uncert = 200.0*sigma/q

      ! Contribution by particle type:
#ifdef PROTONS
      do kpar=1,4
#else
      do kpar=1,3
#endif
        q = edep(kpar)*invn
        q2 = edep2(kpar)*invn
        sigma = (q2-q**2)*invn
        sigma = sqrt(max(sigma,0.0))
        write(out,'(i2,1x,es12.5,1x,es8.1)') kpar,q,2.0*sigma
      enddo

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
     & '#   Average uncertainty in % [uncert]:'
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
      end


      subroutine PCSinitally(mem)
!*******************************************************************
!*    Initializes. To be called before TALLY.                      *
!*                                                                 *
!*    Output:                                                      *
!*      mem: Memory used by tally arrays                           *
!*******************************************************************
      use PCSmod
      use dataTypesMod
      implicit none
      real*8 mem

      character*(*) secid,eos
      parameter (secid=
     &  '[SECTION TALLY PARTICLE CURRENT SPECTRUM v.2009-06-15]')
      parameter (eos='[END OF PCS SECTION]')
      character*80 buffer,upcase
      integer errorcode
      real*8 emax

      write(*,*) ' '
      write(*,'(a)')
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      call getline(buffer,0)
      if (index(buffer,secid).eq.0) then
        write(*,'(a)') 'PCSinitally:ERROR: incorrect section header;'
        write(*,'(a,a)') '  expecting to find: ',secid
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif
      write(*,'(a)') secid

      ! Report memory usage:
      mem = sizeOfReal8*4*3*real(nbinmax)
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
     &    '>>>> Tally Particle Current Spectrum is OFF >>>>'
        do
          read(*,'(a80)',iostat=errorcode) buffer
          if (errorcode.ne.0) then
            write(*,'(a,a,a)') 'PCSinitally:ERROR: ',
     &       'Unable to find End-Of-Section mark: ',eos
            stop
          endif
          if (index(buffer,eos).ne.0) return
        enddo
      else
        write(*,'(a)')
     &    'PCSinitally:ERROR: expecting to find ON or OFF'
        write(*,'(a)') 'found instead:'
        write(*,'(a)') buffer
        stop
      endif

      write(*,'(a)') 'Detection material set to:'
      read(*,*) matdet
      write(*,'(i3)') matdet

      write(*,'(a)')
     &  'Spectrum limits [Emin,Emax] (eV) and no. bins:'
      read(*,*) emin,emax,nbin
      write(*,'(2(1x,es12.5),1x,i5)') emin,emax,nbin
      ebin = (emax-emin)/nbin
      iebin = 1.0/ebin
      ! Add one bin to store, separately, counts with E=Emax:
      nbin = nbin+1
      if (nbin.gt.nbinmax) then
        write(*,'(a)') 'PCSinitally:ERROR: Too many bins.'
        write(*,'(a,i0)') '  Max no. bins is ',nbinmax-1
        stop
      endif

      write(*,'(a)') 'Relative uncertainty (%) requested:'
      read(*,*) unclimit
      write(*,'(1x,es12.5)') unclimit

      ! Clear counters:
      tedep2  = 0.0
      edeptmp = 0.0
      edep    = 0.0
      edep2   = 0.0
      tcount2 = 0.0
      countmp = 0.0
      counts  = 0.0
      count2  = 0.0

      read(*,'(a80)') buffer
      if (index(buffer,eos).eq.0) then
        write(*,*) 'PCSinitally:ERROR: End-Of-Section mark not found'
        write(*,'(a,a)') '  expecting to find: ',eos
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif

      call PCSdump(0)  ! Sim restart file

      write(*,'(a)') '>>>> PCS tally initialization finished >>>>'
      end


      subroutine PCSdump(mode)
!*******************************************************************
!*    Dumps into or reads data from a dump file.                   *
!*                                                                 *
!*    Input:                                                       *
!*      mode:   1 to write dump file, else to read from it.        *
!*******************************************************************
      use PCSmod, only : counts,count2,tcount2,edep,edep2,tedep2,nbin
      use dumpMod
      implicit none
      integer mode

      if (mode.eq.1) then
        if (dumpfile.lt.0) return  ! No dump file open
        write(dumpfile) counts(1:nbin,1:4),count2(1:nbin,1:4),
     &    tcount2(1:4),edep(1:4),edep2(1:4),tedep2
      else
        if (restartfile.lt.0) return  ! No restart file open
        read(restartfile) counts(1:nbin,1:4),count2(1:nbin,1:4),
     &    tcount2(1:4),edep(1:4),edep2(1:4),tedep2
      endif
      end


!>>>> End Of File >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
