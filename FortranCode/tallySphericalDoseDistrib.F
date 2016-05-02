!*******************************************************************
!*                          TALLY                                  *
!*                 SPHERICAL DOSE DISTRIBUTION                     *
!*                                                                 *
!* Short description:                                              *
!*   Dose distribution in a grid of spherical shells. The grid is  *
!*   superimposed on the solid PENGEOM geometry.                   *
!* Rev:                                                            *
!*   2015-01-07                                                    *
!*******************************************************************


!*******************************************************************
!*******************************************************************
!*    Modules section                                              *
!*******************************************************************
!*******************************************************************

      module SPDmod
!*******************************************************************
!*    Vars for the SPD tally.                                      *
!*******************************************************************
      implicit none
      save
      logical active
      integer prtxyz,nbin
      integer, parameter :: nbinmax=32000
      real*8 edptmp(nbinmax),edep(nbinmax),edep2(nbinmax)
      real*8 idens(nbinmax),nlast(nbinmax)
      real*8 rmin,r2min,r2max,dr,idr,unclimit
      end


!*******************************************************************
!*******************************************************************
!*    Routines section                                             *
!*******************************************************************
!*******************************************************************

      subroutine SPDtally(mode,arg)
!*******************************************************************
!*    Input:                                                       *
!*      mode:   identifies the state of the calling routine        *
!*      arg:   energy loss (mode<0) or history no. (mode=1)        *
!*******************************************************************
      use SPDmod
      use TRACK_mod, only : x,y,z,wght
      use ctrsimMod, only : nhist  ! Current history number
      implicit none
      integer mode
      real*8 arg

      integer*4 i
      real*8 r2

      if (.not.active) return

      if (mode.le.0) then
        if (arg.eq.0.0d0) return             ! Nothing to deposit

        r2 = x*x+y*y+z*z
        if (r2.lt.r2min.or.r2.gt.r2max) return ! Not in ROI
        i = (sqrt(r2)-rmin)*idr+1.0          ! Bin index where score will be added

        if (nhist.gt.nlast(i)) then          ! Transfer partial to totals when a new history visits
          edep(i)  = edep(i) +edptmp(i)
          edep2(i) = edep2(i)+edptmp(i)**2
          edptmp(i)= arg*wght
          nlast(i) = nhist+0.5d0             ! Add 1/2 to avoid roundoff errors
        else
          edptmp(i) = edptmp(i)+arg*wght     ! Same history as before, keep adding to temp counter
        endif

      endif
      end


      subroutine SPDreport(mode,n,cputim,uncdone)
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
      use SPDmod
      use RSEEDcommonMod
      implicit none
      integer mode,uncdone
      real*8 n,cputim

      integer out,finduf,errorcode
      integer*4 nchan,i
      real*8 q,sigma,eff,avesig,maxq,fact,uncert,r,rave,rnext
      real*8 threefour,invdr3,pi,inv43pi,invn
      parameter (threefour=3.0d0/4.0d0)
      parameter (pi=3.1415926535897932d0,inv43pi=3.0/(4.0*pi))

      uncdone = 1
      if (.not.active) return

      ! Dump temp counters and obtain max score:
      do i=1,nbin
        if (nlast(i).lt.0.5) cycle      ! Skip empty bins
        edep(i)  = edep(i) +edptmp(i)   ! Transfer temp counter
        edep2(i) = edep2(i)+edptmp(i)**2
        edptmp(i)= 0.0                  ! Reset counter
        nlast(i) = 0.0                  ! Reset last visited to avoid recounting in next report
      enddo

      if (mode.ne.1) call SPDdump(1)  ! Sim dump file, notice that edptmp and nlast have been cleared before
      if (mode.eq.0) return           ! Dump only, no report

      ! Prepare output files:
      out = finduf()
      open(out,file='tallySphericalDoseDistrib.dat',iostat=errorcode)
      if (errorcode.ne.0) then
        write(*,*) ''
        write(*,'(a)')
     &    '*********************************************'
        write(*,'(a)')
     &    'SPDreport:ERROR: cannot open output data file'
        write(*,'(a)')
     &    '*********************************************'
        close(out)  ! Just in case
        return
      endif

      ! Write header:
      write(out,'(a)')
     &'#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      write(out,'(a)') '# [SECTION REPORT SPHERICAL DOSE DISTRIB]'
      write(out,'(a)') '# Dose units are eV/g per history'
      write(out,'(a)') '#'
      write(out,'(a)') '# No. of radial (r) bins:'
      write(out,'(a,1x,i0)') '# ',nbin
      write(out,'(a)') '# Min r and bin width (cm):'
      write(out,'(a,2(1x,es12.5))') '# ',rmin,dr
      if (prtxyz.eq.1) then
        write(out,'(a)') '#'
        write(out,'(a)')
     &    '# For plotting purposes, two radii r for each bin '//
     &    'are given, namely, the low end and an average radius; '//
     &    'the latter is a weighted average with a weight '//
     &    'proportional to r^2.'
      endif
      write(out,'(a)') '#'
      write(out,'(a,$)') '# '
      if (prtxyz.eq.1) then
        write(out,'(a,$)') 'rBinIndex : rLow(cm) : rAve(cm) : '
      endif
      write(out,'(a)') 'dose : +-2sigma'

      ! Write data:
      avesig = 0.0
      nchan  = 0
      invn = 1.0/n
      maxq = 0.5*maxval(edep)         ! 1/2 of the max score

      do i=1,nbin
        r = rmin+dr*(i-1)
        rnext = r+dr
        invdr3 = 1.0/(rnext**3-r**3)
        fact = idens(i)*inv43pi*invdr3             ! This is 1/Delta_mass
        if (prtxyz.eq.1) then
          rave = threefour*(rnext**4-r**4)*invdr3  ! Average with weight(r)~r^2
          write(out,'(1x,i5,2(1x,es12.5),$)') i,r,rave
        endif
        q = edep(i)*invn
        sigma = sqrt(max((edep2(i)*invn-q**2)*invn,0.0))*fact
        q = q*fact
        write(out,'(1x,es12.5,1x,es8.1)') q,2.0*sigma
        if (edep(i).gt.maxq) then                  ! Evaluate average uncertainty for scores above 1/2 max score
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


      subroutine SPDinitally(mem)
!*******************************************************************
!*    Initializes. To be called before TALLY                       *
!*                                                                 *
!*    Output:                                                      *
!*      mem: Memory used by tally arrays                           *
!*******************************************************************
      use SPDmod
      use dataTypesMod
      use TRACK_mod, only : x,y,z,u,v,w
      implicit none
      real*8 mem

      character*(*) secid,eos
      parameter (secid=
     &  '[SECTION TALLY SPHERICAL DOSE DISTRIB v.2009-06-15]')
      parameter (eos='[END OF SPD SECTION]')
      character*80 buffer,upcase
      integer errorcode
      integer*4 i
      real*8 rmax,localdens,flocaldens
      real*8, parameter :: oneplus =1.0d0+1.0d-12  ! A bit more than one
      real*8, parameter :: oneminus=1.0d0-1.0d-12  ! A bit less than one

      write(*,*) ' '
      write(*,'(a)')
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      call getline(buffer,0)
      if (index(buffer,secid).eq.0) then
        write(*,'(a)') 'SPDinitally:ERROR: incorrect section header;'
        write(*,'(a,a)') '  expecting to find: ',secid
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif
      write(*,'(a)') secid

      ! Report memory usage:
      mem = sizeOfReal8*5*real(nbinmax)
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
     &    '>>>> Tally Spherical Dose Distrib is OFF >>>>'
        do
          read(*,'(a80)',iostat=errorcode) buffer
          if (errorcode.ne.0) then
            write(*,'(a,a,a)') 'SPDinitally:ERROR: ',
     &       'Unable to find End-Of-Section mark: ',eos
            stop
          endif
          if (index(buffer,eos).ne.0) return
        enddo
      else
        write(*,'(a)')
     &    'SPDinitally:ERROR: expecting to find ON or OFF'
        write(*,'(a)') 'found instead:'
        write(*,'(a)') buffer
        stop
      endif

      ! Read geometry parameters:
      write(*,'(a)') 'rmin,rmax,nbin:'
      read(*,*) rmin,rmax,nbin
      write(*,'(2(1x,es12.5),1x,i0)') rmin,rmax,nbin
      if (rmin.lt.0.0.or.rmin.ge.rmax.or.nbin.lt.1) then
        write(*,'(a)') 'SPDinitally:ERROR: Invalid entry.'
        stop
      endif
      if (nbin.gt.nbinmax) then
        write(*,'(a)') 'SPDinitally:ERROR: Too many bins.'
        write(*,'(a,i0)') '  Max no. bins is ',nbinmax
        stop
      endif
      dr = (rmax-rmin)/nbin
      idr = 1.0/dr
      r2min = rmin**2*oneplus  ! (1+-eps) tolerances avoid out of range bin indexes in tally routine
      r2max = rmax**2*oneminus

      write(*,'(a)') 'Print coordinates:'
      read(*,*) prtxyz
      if (prtxyz.eq.1) then
        write(*,'(a)') ' yes'
      else
        write(*,'(a)') ' no'
      endif

      write(*,'(a)') 'Relative uncertainty (%) requested:'
      read(*,*) unclimit
      write(*,'(1x,es12.5)') unclimit

      ! Init arrays:
      edptmp = 0.0
      edep   = 0.0
      edep2  = 0.0
      nlast  = 0.0
      idens  = 0.0
      u = 0.0                           ! For locate() to operate properly
      v = 0.0
      w = 1.0
      x = 0.0
      y = 0.0
      do i=1,nbin
        z = rmin+dr*(i-0.5)             ! This is to locate a point and find its material
        localdens = flocaldens()        ! Local mass density--quads or voxs
        if (localdens.gt.0.0) idens(i) = 1.0/localdens
      enddo

      read(*,'(a80)') buffer
      if (index(buffer,eos).eq.0) then
        write(*,'(a)')
     &    'SPDinitally:ERROR: End-Of-Section mark not found'
        write(*,'(a,a)') '  expecting to find: ',eos
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif

      call SPDdump(0)  ! Sim restart file

      write(*,'(a)') '>>>> SPD tally initialization finished >>>>'
      end


      subroutine SPDdump(mode)
!*******************************************************************
!*    Dumps into or reads data from a dump file.                   *
!*                                                                 *
!*    Input:                                                       *
!*      mode:   1 to write dump file, else to read from it.        *
!*******************************************************************
      use SPDmod, only : edep,edep2,nbin
      use dumpMod
      implicit none
      integer mode

      if (mode.eq.1) then
        if (dumpfile.lt.0) return  ! No dump file open
        write(dumpfile) edep(1:nbin),edep2(1:nbin)
      else
        if (restartfile.lt.0) return  ! No restart file open
        read(restartfile) edep(1:nbin),edep2(1:nbin)
      endif
      end


!>>>> End Of File >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
