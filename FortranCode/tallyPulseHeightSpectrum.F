!*******************************************************************
!*                   TALLY PULSE HEIGHT SPECTRUM                   *
!*                                                                 *
!* Short description:                                              *
!*   Calculates the spectrum of the energy deposited in a given    *
!*   material.                                                     *
!* Rev:                                                            *
!*   2015-01-07                                                    *
!*******************************************************************


!*******************************************************************
!*******************************************************************
!*    Modules section                                              *
!*******************************************************************
!*******************************************************************

      module PHSmod
!*******************************************************************
!*    Vars for the PHS tally.                                      *
!*******************************************************************
      implicit none
      save
      logical active,isGaussConvo
      integer matdet,nbin
      integer, parameter :: nbinmax=32000
      real*8 counts(nbinmax),count2(nbinmax)
      real*8 edptmp,emin,ebin,iebin,unclimit,a0,a1
      end


!*******************************************************************
!*******************************************************************
!*    Routines section                                             *
!*******************************************************************
!*******************************************************************

      subroutine PHStally(mode,eloss)
!*******************************************************************
!*    Input:                                                       *
!*      mode:   identifies the state of the calling procedure      *
!*      eloss:   energy deposition                                 *
!*******************************************************************
      use PHSmod
      use TRACK_mod, only : mat,wght
      implicit none
      integer mode
      real*8 eloss

      integer bin
      real*8 gausspdf

      if (.not.active) return

      if (mode.le.0) then                ! Energy has been deposited
        if (mat.ne.matdet) return        ! Not in detector material
        edptmp = edptmp+eloss            ! Update temporary E counter

      else if (mode.eq.6) then           ! End-of-history bookkeeping
        if (isGaussConvo)                ! Sample deposited energy from a Gaussian
     &    edptmp = edptmp+sqrt(a0+a1*edptmp)*gausspdf() ! Gaussian spread, variance is linear with Edep

        bin = (edptmp-emin)*iebin+1.0    ! E bin number
        if (bin.gt.0.and.bin.le.nbin) then
          counts(bin) = counts(bin)+wght ! Score one count
          count2(bin) = count2(bin)+wght**2  ! Score variance
        endif
        edptmp = 0.0                     ! Clear counter to start a fresh history

      endif
      end


      subroutine PHSreport(mode,n,cputim,uncdone)
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
      use PHSmod
      use RSEEDcommonMod
      implicit none
      integer mode,uncdone
      real*8 n,cputim

      integer i,nchan,out,finduf,errorcode
      real*8 q,q2,sigma,eff,avesig,maxq,emiddle,elow,uncert,invn

      uncdone = 1
      if (.not.active) return

      if (mode.ne.1) call PHSdump(1)  ! Sim dump file
      if (mode.eq.0) return           ! Dump only, no report

      ! Prepare output files:
      out = finduf()
      open(out,file='tallyPulseHeightSpectrum.dat',iostat=errorcode)
      if (errorcode.ne.0) then
        write(*,*)
        write(*,'(a)')
     &    '*********************************************'
        write(*,'(a)')
     &    'PHSreport:ERROR: cannot open output data file'
        write(*,'(a)')
     &    '*********************************************'
        close(out)  ! Just in case
        return
      endif

      write(out,'(a)')
     &'#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      write(out,'(a)') '# [SECTION REPORT PULSE HEIGHT SPECTRUM]'
      write(out,'(a)') '# Units are 1/eV per history'
      write(out,'(a)') '#'
      write(out,'(a)') '# Detection material:'
      write(out,'(a,i0)') '#  ',matdet
      write(out,'(a)') '# Number of energy bins:'
      write(out,'(a,i0)') '#  ',nbin
      write(out,'(a)') '# Energy bin width (eV):'
      write(out,'(a,es12.5)') '#  ',ebin
      write(out,'(a)') '#'
      write(out,'(a)') '# For plotting purposes, '//
     &  'two values per energy bin are given, '//
     &  'namely, the low end and the middle point of each bin.'
      write(out,'(a)') '#'
      write(out,'(a)')
     &'# Elow(eV) : Emiddle(eV) : counts(1/eV/hist) : +-2sigma : nbin'

      avesig = 0.0
      nchan = 0
      invn = 1.0/n
      maxq = 0.5*maxval(counts)  ! 1/2 of the max score

      do i=1,nbin
        q  = counts(i)*invn
        q2 = count2(i)*invn
        sigma = (q2-q**2)*invn
        sigma = sqrt(max(sigma,0.0))
        q = q*iebin
        sigma = sigma*iebin
        elow = emin+ebin*(i-1)
        emiddle = elow+ebin*0.5
        write(out,'(3(1x,es12.5),1x,es8.1,1x,i0)')
     &    elow,emiddle,q,2.0*sigma,i
        ! Evaluate average uncertainty for scores above 1/2 max score:
        if (counts(i).gt.maxq) then
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
      end


      subroutine PHSinitally(mem)
!*******************************************************************
!*    Initializes. To be called before TALLY.                      *
!*                                                                 *
!*    Output:                                                      *
!*      mem: Memory used by tally arrays                           *
!*******************************************************************
      use PHSmod
      use dataTypesMod
      implicit none
      real*8 mem

      character*(*) secid,eos
      parameter (secid=
     &  '[SECTION TALLY PULSE HEIGHT SPECTRUM v.2012-06-01]')
      parameter (eos='[END OF PHS SECTION]')
      character*80 buffer,upcase
      integer errorcode
      real*8 emax,fwhm2sig2,eps
      parameter (fwhm2sig2=1.0/(8.0*log(2.0d0)),eps=1.0d-3)

      write(*,*) ' '
      write(*,'(a)')
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      call getline(buffer,0)
      if (index(buffer,secid).eq.0) then
        write(*,'(a)') 'PHSinitally:ERROR: incorrect section header;'
        write(*,'(a,a)') '  expecting to find: ',secid
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif
      write(*,'(a)') secid

      ! Report memory usage:
      mem = sizeOfReal8*2*real(nbinmax)
      write(*,'(a)') 'Memory used by static arrays (MB):'
      write(*,'(1x,es12.5)') mem*1.0e-6

      read(*,'(a80)') buffer
      buffer = adjustl(upcase(buffer))
      buffer = buffer(1:scan(buffer,' ')) ! Clip at 1st blank
      if (buffer.eq.'ON') then
        active = .true.
      else if (buffer.eq.'OFF') then
        active = .false.
        write(*, '(a)') '>>>> Tally Pulse Height Spectrum is OFF >>>>'
        do
          read(*,'(a80)',iostat=errorcode) buffer
          if (errorcode.ne.0) then
            write(*,'(a,a,a)') 'PHSinitally:ERROR: ',
     &       'Unable to find End-Of-Section mark: ',eos
            stop
          endif
          if (index(buffer,eos).ne.0) return
        enddo
      else
        write(*,'(a)')
     &    'PHSinitally:ERROR: expecting to find ON or OFF'
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
      write(*,'(2(1x,es12.5),1x,i0)') emin,emax,nbin
      if (emax-emin.lt.eps.or.nbin.lt.1) then
        write(*,'(a)') 'PHSinitally:ERROR: invalid entry.'
        stop
      endif
      ebin = (emax-emin)/nbin
      iebin = 1.0/ebin
      ! Add one bin to store, separately, counts with E=Emax:
      nbin = nbin+1
      if (nbin.gt.nbinmax) then
        write(*,'(a)') 'PHSinitally:ERROR: Too many bins.'
        write(*,'(a,i0)') '  Max no. bins is ',nbinmax-1
        stop
      endif

      write(*,'(a)') 'Parameters for Gaussian convolution:'
      read(*,*) a0,a1
      write(*,'(2(1x,es12.5))') a0,a1
      if (min(a0,a1).lt.0.0) then
        write(*,'(a)') 'PHSinitally:ERROR: Invalid entry.'
        stop
      endif
      isGaussConvo = .false.
      if (max(a0,a1).gt.0.0) isGaussConvo = .true.
      a0 = a0*fwhm2sig2  ! This is in units of sigma^2
      a1 = a1*fwhm2sig2

      write(*,'(a)') 'Relative uncertainty (%) requested:'
      read(*,*) unclimit
      write(*,'(1x,es12.5)') unclimit

      ! Clear counters:
      edptmp = 0.0
      counts = 0.0
      count2 = 0.0

      read(*,'(a80)') buffer
      if (index(buffer,eos).eq.0) then
        write(*,*) 'PHSinitally:ERROR: End-Of-Section mark not found'
        write(*,'(a,a)') '  expecting to find: ',eos
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif

      call PHSdump(0)  ! Sim restart file

      write(*,'(a)') '>>>> PHS tally initialization finished >>>>'
      end


      subroutine PHSdump(mode)
!*******************************************************************
!*    Dumps into or reads data from a dump file.                   *
!*                                                                 *
!*    Input:                                                       *
!*      mode:   1 to write dump file, else to read from it.        *
!*******************************************************************
      use PHSmod, only : counts,count2,nbin
      use dumpMod
      implicit none
      integer mode

      if (mode.eq.1) then
        if (dumpfile.lt.0) return  ! No dump file open
        write(dumpfile) counts(1:nbin),count2(1:nbin)
      else
        if (restartfile.lt.0) return  ! No restart file open
        read(restartfile) counts(1:nbin),count2(1:nbin)
      endif
      end


!>>>> End Of File >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
