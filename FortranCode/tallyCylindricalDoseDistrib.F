!*******************************************************************
!*                          TALLY                                  *
!*            CYLINDRICAL (rho-z) DOSE DISTRIBUTION                *
!*                                                                 *
!* Short description:                                              *
!*   Dose distribution in a grid of (rho,z) volumetric bins,       *
!*   where rho^2 is (x^2+y^2). The grid of bins is superimposed on *
!*   the solid PENGEOM geometry.                                   *
!* Rev:                                                            *
!*   2015-01-07                                                    *
!*******************************************************************


!*******************************************************************
!*******************************************************************
!*    Modules section                                              *
!*******************************************************************
!*******************************************************************

      module CDDmod
!*******************************************************************
!*    Vars for the CDD tally.                                      *
!*******************************************************************
      implicit none
      save
      logical active
      integer prtxyz
      integer*4 nr,nz
      integer*4, parameter :: nbinmax=100000
      real*8 edptmp(nbinmax),edep(nbinmax),edep2(nbinmax)
      real*8 idens(nbinmax),nlast(nbinmax)
      real*8 dr,dz,idr,idz,unclimit,rmin,zmin,r2min,r2max
      end


!*******************************************************************
!*******************************************************************
!*    Routines section                                             *
!*******************************************************************
!*******************************************************************

      subroutine CDDtally(mode,arg)
!*******************************************************************
!*    Input:                                                       *
!*      mode:   identifies the state of the calling routine        *
!*      arg:   energy loss (mode<0) or history no. (mode=1)        *
!*******************************************************************
      use CDDmod
      use TRACK_mod
      use ctrsimMod, only : nhist  ! Current history number
      implicit none
      integer mode
      real*8 arg

      integer*4 bin,i,k
      real*8 r2
      real*8, parameter :: oneplus=1.0d0+1.0d-12 ! A bit more than one

      if (.not.active) return

      if (mode.le.0) then
        if (arg.eq.0.0) return              ! Nothing to deposit

        k = (z-zmin)*idz+oneplus            ! ONEPLUS ensures proper truncation when IDZ=0
        if (k.lt.1.or.k.gt.nz) return       ! Not in ROI
        r2 = x*x+y*y
        if (r2.lt.r2min.or.r2.gt.r2max) return ! Not in ROI
        i = (sqrt(r2)-rmin)*idr+1.0
        bin = i+(k-1)*nr                    ! Maps i,k into a single index

        if (nhist.gt.nlast(bin)) then       ! Transfer score to totals only when a new history visits
          edep(bin)  = edep(bin) +edptmp(bin)
          edep2(bin) = edep2(bin)+edptmp(bin)**2
          edptmp(bin)= arg*wght             ! Record new score
          nlast(bin) = nhist+0.5            ! Add 1/2 to avoid roundoff errors
        else
          edptmp(bin) = edptmp(bin)+arg*wght ! Same history as last time, keep adding
        endif

      endif
      end


      subroutine CDDreport(mode,n,cputim,uncdone)
!*******************************************************************
!*    Input:                                                       *
!*      mode:   -1 if called at end-of-simulation, 0 for dump only,*
!*               1 for tally report only, 2 for both.              *
!*      n:   no. of histories simulated                            *
!*      cputim:   elapsed CPU time                                 *
!*    Output:                                                      *
!*      uncdone:   2 if uncert reached, 1 if not defined, 0 else;  *
!*                 not computed if mode=0                          *
!*******************************************************************
      use CDDmod
      use RSEEDcommonMod
      implicit none
      integer mode,uncdone
      real*8 n,cputim

      character*80 buffer
      integer nchan,out,finduf,errorcode
      integer*4 bin,i,k,nzef
      real*8 q,sigma,eff,avesig,maxq,fact,r,z,uncert,zmiddle,rave
      real*8 twothird,invdr2,pi,invpi,invn,factidz
      parameter (twothird=2.0d0/3.0d0)
      parameter (pi=3.1415926535897932d0,invpi=1.0/pi)

      uncdone = 1
      if (.not.active) return

      ! Dump temp counters:
      do k=1,nz
        bin = (k-1)*nr
        do i=1,nr
          bin = bin+1                             ! Maps i,k into a single index
          if (nlast(bin).lt.0.5) cycle            ! Skip empty bins
          edep(bin)  = edep(bin) +edptmp(bin)     ! Transfer temp counter
          edep2(bin) = edep2(bin)+edptmp(bin)**2
          edptmp(bin)= 0.0                        ! Reset counter
          nlast(bin) = 0.0                        ! Reset last visited to avoid recounting in next report
        enddo
      enddo

      if (mode.ne.1) call CDDdump(1)  ! Sim dump file, notice that edptmp and nlast have been cleared before
      if (mode.eq.0) return           ! Dump only, no report

      ! Prepare z factors and filename:
      if (idz.eq.0.0d0) then
        nzef = 0
        factidz = 1.0
        buffer = 'tallyCylindricalDoseDistrib-r.dat'
      else
        nzef = nz
        factidz = idz
        buffer = 'tallyCylindricalDoseDistrib-rz.dat'
      endif

      ! Prepare output file:
      out = finduf()
      open(out,file=buffer,iostat=errorcode)
      if (errorcode.ne.0) then
        write(*,*) ''
        write(*,'(a)')
     &    '*********************************************'
        write(*,'(a)')
     &    'CDDreport:ERROR: cannot open output data file'
        write(*,'(a)')
     &    '*********************************************'
        close(out)  ! Just in case
        return
      endif

      ! Write header:
      write(out,'(a)')
     &'#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      write(out,'(a)') '# [SECTION REPORT CYLINDRICAL DOSE DISTRIB]'
      if (nzef.ne.0) then
        write(out,'(a)') '# Dose units are eV/g per history'
      else
        write(out,'(a)') '# Dose units are eV.cm/g per history'
      endif
      write(out,'(a)') '#'
      write(out,'(a)')
     &  '# No. of radial (r) bins, min value and bin width (cm):'
      write(out,'(a,i0,2(1x,es12.5))') '#  ',nr,rmin,dr
      write(out,'(a)')
     &  '# No. of depth (z) bins, min value and bin width (cm):'
      write(out,'(a,i0,2(1x,es12.5))') '#  ',nz,zmin,dz
      if (prtxyz.eq.1) then
        write(out,'(a)') '#'
        write(out,'(a)')
     &    '# For plotting purposes, two values per bin coordinate '//
     &    'are given, namely, the low end and an average value.'
        write(out,'(a)')
     &    '# For the z coordinate the average is the middle'//
     &    ' point of the bin.'
        write(out,'(a)')
     &    '# For the r coordinate the average is weighted '//
     &    ' with a weight proportional to the radius r.'
      endif
      write(out,'(a)') '#'
      write(out,'(a,$)') '# '
      if (prtxyz.eq.1) then
        write(out,'(a,$)') 'rBinIndex : rLow(cm) : rAve(cm) : '
        if (nzef.ne.0)
     &    write(out,'(a,$)') 'zBinIndex : zLow(cm) : zMiddle(cm) : '
      endif
      write(out,'(a)') 'dose : +-2sigma'

      ! Write data:
      maxq = 0.5*maxval(edep)           ! 1/2 of max score
      nchan  = 0
      avesig = 0.0
      invn = 1.0/n
      do k=1,nz
        z = zmin+dz*(k-1)
        zmiddle = z+dz*0.5
        if (nzef.ne.0.and.prtxyz.ne.1)  ! Since z is not written, give at least a summary
     &    write(out,'(a,i0,a,es12.5)')
     &    '# zBinIndex=',k,' zMiddle(cm)=',zmiddle

        do i=1,nr
          r = rmin+dr*(i-1)
          invdr2 = 1.0/(dr*(dr+2.0*r))  ! Note: dr2 = (r+dr)^2-r^2
          if (prtxyz.eq.1) then
            rave = twothird*((r+dr)**3-r**3)*invdr2  ! Average with weight(r)~r
            write(out,'(1x,i5,2(1x,es12.5),$)') i,r,rave
          endif
          if (nzef.ne.0.and.prtxyz.eq.1)
     &        write(out,'(1x,i5,2(1x,es12.5),$)') k,z,zmiddle
          bin = i+(k-1)*nr
          fact = factidz*idens(bin)*invpi*invdr2
          q = edep(bin)*invn
          sigma = sqrt(max((edep2(bin)*invn-q**2)*invn,0.0))*fact
          q = q*fact
          write(out,'(1x,es12.5,1x,es8.1)') q,2.0d0*sigma

          ! Evaluate average uncertainty for scores above 1/2 max score:
          if (edep(bin).gt.maxq) then
            avesig = avesig+(sigma/q)**2
            nchan = nchan+1
          endif
        enddo

        if (nr.gt.1) write(out,*) ' '  ! Separate data blocks
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
      if (cputim.gt.0.0d0) then
        write(out,'(a)') '#   Speed (histories/s):'
        write(out,'(a,es12.5)') '#   ',n/cputim
      endif
      write(out,'(a)')
     & '#   Average uncertainty (above 1/2 max score) in % [uncert]:'
      write(out,'(a,es12.5)') '#   ',uncert
      eff = n*uncert**2
      if (eff.gt.0.0d0) then
        write(out,'(a)') '#   Intrinsic efficiency [N*uncert^2]^-1:'
        write(out,'(a,es12.5)') '#   ',1.0/eff
      endif
      eff = cputim*uncert**2
      if (eff.gt.0.0d0) then
        write(out,'(a)') '#   Absolute efficiency [t*uncert^2]^-1:'
        write(out,'(a,es12.5)') '#   ',1.0/eff
      endif
      close(out)
      end


      subroutine CDDinitally(mem)
!*******************************************************************
!*    Initializes. To be called before TALLY                       *
!*                                                                 *
!*    Output:                                                      *
!*      mem: Memory used in tally arrays                           *
!*******************************************************************
      use CDDmod
      use dataTypesMod
      use TRACK_mod, only : x,y,z,u,v,w
      implicit none
      real*8 mem

      character*(*) secid,eos
      parameter (secid=
     &  '[SECTION TALLY CYLINDRICAL DOSE DISTRIB v.2009-06-15]')
      parameter (eos='[END OF CDD SECTION]')
      character*80 buffer,upcase
      integer errorcode
      integer*4 i,k,bin
      real*8 rmax,zmax,rnbin,localdens,flocaldens
      real*8, parameter :: oneplus =1.0d0+1.0d-12  ! A bit more than one
      real*8, parameter :: oneminus=1.0d0-1.0d-12  ! A bit less than one

      write(*,*) ' '
      write(*,'(a)')
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      call getline(buffer,0)
      if (index(buffer,secid).eq.0) then
        write(*,'(a)') 'CDDinitally:ERROR: incorrect section header;'
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
     &    '>>>> Tally Cylindrical Dose Distrib is OFF >>>>'
        do
          read(*,'(a80)',iostat=errorcode) buffer
          if (errorcode.ne.0) then
            write(*,'(a,a,a)') 'CDDinitally:ERROR: ',
     &       'Unable to find End-Of-Section mark: ',eos
            stop
          endif
          if (index(buffer,eos).ne.0) return
        enddo
      else
        write(*,'(a)')
     &    'CDDinitally:ERROR: expecting to find ON or OFF'
        write(*,'(a)') 'found instead:'
        write(*,'(a)') buffer
        stop
      endif

      ! Read geometry parameters:
      write(*,'(a)') 'rmin,rmax,nr:'
      read(*,*) rmin,rmax,nr
      write(*,'(2(1x,es12.5),1x,i0)') rmin,rmax,nr
      if (rmin.lt.0.0.or.rmin.ge.rmax.or.nr.lt.1) then
        write(*,'(a)') 'CDDinitally:ERROR: invalid entry.'
        stop
      endif
      dr = (rmax-rmin)/nr
      idr = 1.0/dr
      r2min = rmin**2*oneplus  ! (1+-eps) tolerances avoid out of range bin indexes in tally routine
      r2max = rmax**2*oneminus
      write(*,'(a)') 'zmin,zmax,nz:'
      read(*,*) zmin,zmax,nz
      write(*,'(2(1x,es12.5),1x,i0)') zmin,zmax,nz
      if (zmin.gt.zmax.or.nz.lt.0.or.(nz.gt.0.and.zmin.eq.zmax)) then
        write(*,*) 'CDDinitally:ERROR: zmin >= zmax  or  nz < 0'
        stop
      endif
      if (nz.eq.0) then
        dz = zmax-zmin
        idz = 0.0d0  ! Keep it double, it will be compared later with 0.0d0
        nz = 1
      else
        dz = (zmax-zmin)/nz
        idz = 1.0/dz
      endif
      rnbin = dble(nr)*dble(nz)
      if (rnbin.gt.real(nbinmax)) then
        write(*,'(a)') 'CDDinitally:ERROR: Too many bins.'
        write(*,'(a,i0)') '  Max no. bins is ',nbinmax
        stop
      endif

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
      u = 0.0  ! This is for LOCATE to operate properly
      v = 0.0
      w = 1.0
      y = 0.0
      do k=1,nz
        z = zmin+dz*(k-0.5)
        do i=1,nr
          x = rmin+dr*(i-0.5)
          bin = i+(k-1)*nr
          localdens = flocaldens()  ! Local mass density--quads or voxs
          if (localdens.gt.0.0) idens(bin) = 1.0/localdens
        enddo
      enddo

      read(*,'(a80)') buffer
      if (index(buffer,eos).eq.0) then
        write(*,'(a)')
     &    'CDDinitally:ERROR: End-Of-Section mark not found'
        write(*,'(a,a)') '  expecting to find: ',eos
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif

      call CDDdump(0)  ! Sim restart file

      write(*,'(a)') '>>>> CDD tally initialization finished >>>>'
      end


      subroutine CDDdump(mode)
!*******************************************************************
!*    Dumps into or reads data from a dump file.                   *
!*                                                                 *
!*    Input:                                                       *
!*      mode:   1 to write dump file, else to read from it.        *
!*******************************************************************
      use CDDmod, only : edep,edep2,nr,nz
      use dumpMod
      implicit none
      integer mode

      integer*4 nbin

      nbin = nr*nz
      if (mode.eq.1) then
        if (dumpfile.lt.0) return  ! No dump file open
        write(dumpfile) edep(1:nbin),edep2(1:nbin)
      else
        if (restartfile.lt.0) return  ! No restart file open
        read(restartfile) edep(1:nbin),edep2(1:nbin)
      endif
      end


!>>>> End Of File >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
