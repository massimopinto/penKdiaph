!*******************************************************************
!*                         SOURCE                                  *
!*               BOX-ISOTROPIC-GAUSS-SPECTRUM                      *
!*                                                                 *
!* Short description:                                              *
!*   Generation of primary particle states for PENELOPE.           *
!* Rev:                                                            *
!*   2015-04-30                                                    *
!*******************************************************************


!*******************************************************************
!*******************************************************************
!*    Modules section                                              *
!*******************************************************************
!*******************************************************************

      module BIGSmod
!*******************************************************************
!*    Vars for the BIGS source.                                    *
!*******************************************************************
      implicit none
      logical active,isPosGauss
      integer parsrc,matsrc,nspc,xipol
      integer, parameter :: nemax=32000
      real*8 espc(nemax),pspc(nemax),despc(nemax),esigma
      real*8 rotv(3,3),rotbox(3,3),cossrc0,dcossrc,phi0,dphi
      real*8 boxside(3),ctrans(3),xsigma,ysigma,xsp1,xsp2,xsp3
      real*8, parameter :: kfact=6.0d0   ! Maximum num. of sigmas in Gaussian dispersion
      end


!*******************************************************************
!*******************************************************************
!*    Routines section                                             *
!*******************************************************************
!*******************************************************************

      subroutine BIGSsource
!*******************************************************************
!*    Output:                                                      *
!*      The generated particle is pushed to PENELOPE's stack.      *
!*******************************************************************
      use BIGSmod
      use TRACK_mod
      implicit none

      integer*4 ncross
      real*8 dsef
      real*8, parameter :: infty=1.0e30

      if (.not.active) return

      kpar = parsrc  ! Particle type
      wght   = 1.0   ! Init statistical weight
      ipol = xipol   ! Polarization--effective for photons only
      sp1 = xsp1
      sp2 = xsp2
      sp3 = xsp3
      ilb(1) = 1     ! Tag as source particle (i.e. 1st generation)
      ilb(2) = 0     ! Clear other labels
      ilb(3) = 0
      ilb(4) = 0
      ilb(5) = 0     ! Optional label (transferred to descendants)

      if (matsrc.ge.0) then           ! Local source
        call BIGSsampleDirection
        call BIGSsamplePosition
      else                            ! In-field beam
        call BIGSsampleInfieldBeam
      endif
      if (mat.eq.0) then              ! The particle was born in vacuum
        call stepx(infty,dsef,ncross) ! Advance up to the object or infinity
        call tally(7,dsef)            ! Inform tallies about the maiden flight
      endif

      call BIGSsampleEnergy

      call stores(e,x,y,z,u,v,w,wght,kpar,ilb,ipol)  ! Push particle to stack
      call tally(0,e)                 ! Inform that primary particle was pushed to stack
      end


      subroutine BIGSsampleDirection
!*******************************************************************
!*    Samples a particle's direction of flight.                    *
!*                                                                 *
!*    Output:                                                      *
!*      The generated position is returned via the TRACK module.   *
!*******************************************************************
      use BIGSmod
      use TRACK_mod
      implicit none
      real*8, external :: rand
      real*8 d(3),sintheta,phi

      d(3) = cossrc0+dcossrc*rand(0.0d0)
      sintheta = sqrt(1.0-d(3)**2)
      phi = phi0+dphi*rand(0.0d0)
      d(2) = sintheta*sin(phi)
      d(1) = sintheta*cos(phi)  ! d is the flight direction in Sv frame
      d = matmul(rotv,d)        ! Rotate direction vector from Sv to Slab
      u = d(1)                  ! Transfer to PENELOPE vars
      v = d(2)
      w = d(3)
      end


      subroutine BIGSsamplePosition
!*******************************************************************
!*    Samples a valid particle's position.                         *
!*                                                                 *
!*    Output:                                                      *
!*      The generated position is returned via the TRACK module.   *
!*    Comments:                                                    *
!*      - Requires that direction of flight be set previously for  *
!*        SAMPLEBOX to operate safely.                             *
!*******************************************************************
      use BIGSmod
      use TRACK_mod
      use ctrsimMod, only : nhist
      implicit none

      logical, save :: warned=.false.    ! Warning flag
      integer*8, save :: ntrials=0       ! Total number of trials since program start

      do                                 ! Loop until a valid position is found
        ntrials = ntrials+1              ! Count num. attempts to measure efficiency
        call BIGSsampleBox               ! Sample a particle position inside the box

        ! Test sampling efficiency and issue warnings and errors if needed:
        if (.not.warned.and.dble(ntrials).gt.max(nhist*10.0,1.0d3))
     &    then
          write(*,*) ''
          write(*,'(a)') '***************'
          write(*,'(a)') 'samplePosition:WARNING: '//
     &      'sampling source efficiency is lower than 10%.'
          write(*,'(a)') 'Redefine source parameters appropriately.'
          write(*,'(a)')
     &      'Number of histories and sampling efficiency:'
          write(*,'(2x,f18.0,1x,f5.2,a)')
     &      nhist,nhist/dble(ntrials)*100.0,'%'
          write(*,'(a)') '***************'
          write(*,*) ''
          warned = .true.
        endif
        if (dble(ntrials).gt.max(nhist*1000.0,1.0d5)) then
          write(*,'(a)') 'samplePosition:ERROR: '//
     &      'sampling source efficiency is lower than 0.1%.'
          write(*,'(a)') 'Redefine source parameters appropriately.'
          write(*,'(a)')
     &      'Number of histories and sampling efficiency:'
          write(*,'(2x,f18.0,1x,f5.2,a)')
     &      nhist,nhist/dble(ntrials)*100.0,'%'
          stop
        endif

        if (matsrc.le.0.or.matsrc.eq.mat) return  ! Point accepted, else cycle; note MATSRC<0 for in-beam case
      enddo
      end


      subroutine BIGSsampleBox
!*******************************************************************
!*    Samples a position inside the source box.                    *
!*                                                                 *
!*    Output:                                                      *
!*      The generated position is returned via the TRACK module.   *
!*    Comments:                                                    *
!*      - Requires that direction of flight be set previously for  *
!*        LOCATE to operate safely.                                *
!*******************************************************************
      use BIGSmod
      use TRACK_mod
      implicit none
      integer i
      real*8 r(3),gausspdf
      real*8, external :: rand

      do i=1,3
        r(i) = boxside(i)*(rand(dble(i))-0.5)  ! Coordinates in Sbox
      enddo
      if (isPosGauss) then                     ! Gaussian x,y displacement in Sbox
        r(1) = r(1)+xsigma*gausspdf()          ! Apply displacement in Sbox frame
        r(2) = r(2)+ysigma*gausspdf()
      endif

      r = matmul(rotbox,r)+ctrans              ! Apply rotation and translation
      x = r(1)                                 ! Transfer to PENELOPE vars
      y = r(2)
      z = r(3)

      call locate                              ! Set quadric body and material
      call locatevox                           ! Set voxel indices and material, if inside a voxel geometry
      end


      subroutine BIGSsampleInfieldBeam
!*******************************************************************
!*    Samples a valid particle's position and direction for the    *
!*    in-field beam source mode.                                   *
!*                                                                 *
!*    Output:                                                      *
!*      The generated position is returned via the TRACK module.   *
!*******************************************************************
      use BIGSmod
      use partvoxMod
      use TRACK_mod
      implicit none

      integer*4 ncross,bodytmp,mattmp,xvoxtmp,yvoxtmp,zvoxtmp,voxtmp
      real*8 r(3),dsef
      real*8, parameter :: infty=1.0e30

      direction: do
        call BIGSsampleDirection
        call BIGSsamplePosition

        ! Store original position vars, to be recovered later:
        bodytmp = ibody
        mattmp  = mat
        xvoxtmp = xvox
        yvoxtmp = yvox
        zvoxtmp = zvox
        voxtmp  = absvox
        r = (/x,y,z/)

        step: do
          if (mat.eq.-matsrc) then  ! Inside field, restore position and accept this point
            mat    = mattmp
            xvox   = xvoxtmp
            yvox   = yvoxtmp
            zvox   = zvoxtmp
            absvox = voxtmp
            x = r(1)
            y = r(2)
            z = r(3)
            return
          endif
          call stepx(infty,dsef,ncross)  ! Move forward to check if the field is on the path
          if (mat.eq.0.or.ncross.eq.0) cycle direction  ! Gone to infty, field not in the path; try another
        enddo step
      enddo direction
      end


      subroutine BIGSsampleEnergy
!*******************************************************************
!*    Samples particle's energy.                                   *
!*                                                                 *
!*    Output:                                                      *
!*      The generated energy is returned via the TRACK module.     *
!*******************************************************************
      use BIGSmod
      use TRACK_mod
      implicit none
      integer i,seeki
      real*8 randno,sigf,gausspdf
      real*8, external :: rand

      randno = rand(0.0d0)
      i = seeki(pspc,randno,nspc)
      e = espc(i)+(randno-pspc(i))*despc(i)
      if (esigma.ne.0.0) then             ! Gaussian distrib
 10     sigf = gausspdf()
        if (abs(sigf).gt.kfact) goto 10   ! Truncate Gauss at +-kfact
        e = max(e+esigma*sigf,0.0)        ! Prevent negative values
      endif
      end


      subroutine BIGSinisrc(activated,emax,mem)
!*******************************************************************
!*    Initializes the source.                                      *
!*                                                                 *
!*    Output:                                                      *
!*      activated: TRUE if the source is active.                   *
!*      emax: max source energy (eV)                               *
!*      mem: memory used by source arrays.                         *
!*******************************************************************
      use BIGSmod
      use dataTypesMod
      implicit none
      logical activated
      real*8 emax,mem

      character*80 buffer,upcase,spcfilen
      character*(*) secid,eos
      parameter (secid=
     &  '[SECTION SOURCE BOX ISOTROPIC GAUSS SPECTRUM v.2014-12-21]')
      parameter (eos='[END OF BIGS SECTION]')
      integer j,answer,errorcode,spcfile,finduf
      real*8 omega,theta,phi,prob,norm,theta0,theta1
      real*8 boxctr(3),boxtrans(3),dsrc(3)
      real*8 pi,fwhm2sig,deg2rad,mc2,twomc2,zero,oneeps
      parameter (pi=3.1415926535897932d0,deg2rad=pi/180.0)
      parameter (fwhm2sig=1.0/(8.0*log(2.0d0))**0.5)
      parameter (mc2=5.10998918d5,twomc2=2.0*mc2,zero=1.0d-30)
      parameter (oneeps=1.0d0+1.0d-14)

      write(*,*) ' '
      write(*,'(a)')
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      call getline(buffer,0)
      if (index(buffer,secid).eq.0) then
        write(*,'(a)') 'BIGSinisrc:ERROR: incorrect section header;'
        write(*,'(a,a)') '  expecting to find: ',secid
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif
      write(*,'(a)') secid

      read(*,'(a80)') buffer
      buffer = adjustl(upcase(buffer))    ! Adjust and capitalize
      buffer = buffer(1:scan(buffer,' ')) ! Clip at 1st blank
      if (buffer.eq.'ON') then
        active = .true.
        activated = .true.
      else if (buffer.eq.'OFF') then
        active = .false.
        activated = .false.
        write(*, '(a)')
     &    '>>>> Source Box Isotropic Gauss Spectrum is OFF >>>>'
        do
          read(*,'(a80)',iostat=errorcode) buffer
          if (errorcode.ne.0) then
            write(*,'(a,a,a)') 'BIGSinisrc:ERROR: ',
     &       'Unable to find End-Of-Section mark: ',eos
            stop
          endif
          if (index(buffer,eos).ne.0) return
        enddo
      else
        write(*,'(a)')
     &    'BIGSinisrc:ERROR: expecting to find ON or OFF'
        write(*,'(a)') 'found instead:'
        write(*,'(a)') buffer
        stop
      endif

      ! Report memory usage:
      mem = sizeOfReal8*3*real(nemax)
      write(*,'(a)') 'Memory used by static arrays (MB):'
      write(*,'(1x,es12.5)') mem*1.0e-6

      ! Type of particle:
      write(*,'(a)') 'Particle type:'
      read(*,*) parsrc
      write(*,'(1x,i0)') parsrc
      if (parsrc.lt.1.or.parsrc.gt.4) then
        write(*,*) 'BIGSinisrc:ERROR: invalid particle type'
        stop
      endif
#ifndef PROTONS
      if (parsrc.eq.4) then
        write(*,*) 'BIGSinisrc:ERROR: invalid particle type'
        stop
      endif
#endif

      ! Photon polarization:
      read(*,'(a80)') buffer  ! Read a dummy line
      write(*,'(a)') 'Photon polarization physics:'
      read(*,*) xipol
      if (xipol.ne.1.and.xipol.ne.0) then
        write(*,'(a)')
     &    'BIGSinisrc:ERROR: invalid polarization switch, '//
     &    'should be 0 or 1'
        stop
      endif
      if (xipol.eq.0) then
        write(*,'(a)') '  not activated'
      else
        write(*,'(a)') '  activated'
      endif
      if (xipol.eq.1) write(*,'(a)') 'Stokes parameters:'
      read(*,*) xsp1,xsp2,xsp3
      if (xipol.eq.1) write(*,'(3(1x,es12.5))') xsp1,xsp2,xsp3
      if (xipol.eq.1.and.xsp1**2+xsp2**2+xsp3**2.gt.oneeps) then
        write(*,*)
     &    'BIGSinisrc:ERROR: vector P={P1,P2,P3} must be P^2 <= 1'
        stop
      endif

      ! Position:
      read(*,'(a80)') buffer  ! Read a dummy line
      write(*,'(a)') 'Coordinates of box center (cm):'
      read(*,*) boxctr
      write(*,'(3(1x,es12.5))') boxctr
      write(*,'(a)') 'Box sides (cm):'
      read(*,*) boxside
      write(*,'(3(1x,es12.5))') boxside
      write(*,'(a)') 'FWHMs of x and y distribs. (cm):'
      read(*,*) xsigma,ysigma
      write(*,'(2(1x,es12.5))') xsigma,ysigma
      isPosGauss=.false.
      if (xsigma.ne.0.0d0.or.ysigma.ne.0.0d0) isPosGauss=.true.
      xsigma = xsigma*fwhm2sig
      ysigma = ysigma*fwhm2sig
      write(*,'(a)') 'Euler angles for box rotation (deg):'
      read(*,*) omega,theta,phi
      write(*,'(3(1x,es12.5))') omega,theta,phi
      call rotation(omega*deg2rad,theta*deg2rad,phi*deg2rad,rotbox)  ! Rotation Sbox -> Slab
      write(*,'(a)') 'Translation of box center (cm):'
      read(*,*) boxtrans
      write(*,'(3(1x,es12.5))') boxtrans
      ! Note: the position X of a source point is determined as follows:
      !   X = R.Xo+D
      ! where R is the rotation defined by the Euler angles and D the translation.
      ! Xo is the position before rotation and translation, given by
      !   Xo = C + (RAND-0.5)*dX + dGauss
      ! with C the center of the box, dX its sides and dGauss the [x,y]-Gaussian dispersion.
      ! The factor RAND is different for each component of dX. So,
      !   X = (R.C+D) + R.[(rand-0.5)*dX + dGauss]
      ! The vector R.C+D is computed and stored here for later use:
      ctrans = matmul(rotbox,boxctr)+boxtrans

      ! Material:
      write(*,'(a)') 'Source material:'
      read(*,*) matsrc
      write(*,'(1x,i0)') matsrc
      if (matsrc.lt.0) write(*,'(a,i0)')
     &  'In-field beam mode; '//
     &  'particles will be aiming at material ',-matsrc

      ! Direction:
      read(*,'(a80)') buffer  ! Read a dummy line
      write(*,'(a)') 'Direction vector (u,v,w):'
      read(*,*) dsrc
      write(*,'(3(1x,es12.5))') dsrc
      norm = sqrt(sum(dsrc**2))
      if (norm.lt.zero) then
        write(*,*) 'BIGSinisrc:ERROR: null direction.'
        stop
      endif
      dsrc = dsrc/norm                                                 ! Direction vector normalized to 1
      if (max(abs(dsrc(1)),abs(dsrc(2))).gt.0.0) then                  ! Avoid atan(0,0)
        call rotation(0.0d0,acos(dsrc(3)),atan2(dsrc(2),dsrc(1)),rotv) ! Rotation Sv -> Sbox
      else                                                             ! Case dsrc = +z or -z
        call rotation(0.0d0,acos(dsrc(3)),0.0d0,rotv)                  ! Phi=0 is irrelevant
      endif
      write(*,'(a)') 'Polar angle interval, [theta0,theta1](deg):'
      read(*,*) theta0,theta1
      write(*,'(2(1x,es12.5))') theta0,theta1
      if (theta1.lt.theta0) then
        write(*,*) 'BIGSinisrc:ERROR: theta1 is less than theta0.'
        stop
      endif
      cossrc0 = cos(theta0*deg2rad)          ! Converted to rad
      dcossrc = cos(theta1*deg2rad)-cossrc0  ! Negative value
      write(*,'(a)') 'Azimuthal angle interval, phi0 and Dphi (deg):'
      read(*,*) phi0,dphi
      write(*,'(2(1x,es12.5))') phi0,dphi
      if (phi0.lt.0.0d0.or.phi0.ge.360.0d0.or.
     &    dphi.lt.0.0d0.or.dphi.gt.360.0d0) then
        write(*,'(a)') 'BIGSinisrc:ERROR: invalid interval.'
        stop
      endif
      phi0 = phi0*deg2rad
      dphi = dphi*deg2rad
      write(*,'(a)')
     &  'Apply to direction same rotation use for box position?'
      read(*,*) answer
      select case (answer)
      case (0)
        write(*,'(a)') ' No'
      case (1)
        write(*,'(a)') ' Yes'
        rotv = matmul(rotbox,rotv)  ! Combined rotation Sv -> Sbox -> Slab
      case default
        write(*,'(a)') 'BIGSinisrc:ERROR: '//
     &    'Invalid entry. Must be 0 or 1.'
        stop
      end select

      ! Energy:
      read(*,'(a80)') buffer                 ! Read a dummy line
      read(*,'(a80)') buffer                 ! Spectrum filename
      buffer = adjustl(buffer)
      spcfilen = buffer(1:scan(buffer,' '))  ! Clip at 1st blank
      if (spcfilen.eq.'-') then
        write(*,'(a)') 'Reading energy spectrum from config file.'
        read(*,'(a80)') buffer               ! Dummy header line
        write(*,'(a)') '  Energy(eV)  Relat.Probability   Bin#'
      else
        write(*,'(a)') 'Reading energy spectrum from file:'
        write(*,'(1x,a)') spcfilen
        spcfile = finduf()                   ! Find a valid unit for the file
        open(spcfile,file=spcfilen,status='old',iostat=errorcode)
        if (errorcode.ne.0) then
          write(*,'(a)')
     &      'BIGSinisrc:ERROR: unable to open spectrum file.'
          stop
        endif
      endif

      nspc = 0                               ! Init num. E bins
      pspc(1) = 0.0                          ! Init accumulated prob.
      do
        nspc = nspc+1
        if (spcfilen.eq.'-') then            ! Read from stdin
          read(*,*,iostat=errorcode) espc(nspc),prob
        else
          call getline(buffer,spcfile)       ! Read a non-comment line from file
          read(buffer,*,iostat=errorcode) espc(nspc),prob
        endif
        if (errorcode.ne.0) then
          write(*,'(a)') 'BIGSinisrc:ERROR: '//
     &      'invalid entry at data line:'
          write(*,'(1x,i0)') nspc
          stop
        endif
        if (spcfilen.eq.'-')
     &    write(*,'(2(1x,es12.5),1x,i0)') espc(nspc),prob,nspc  ! Write table to stdout
        if (espc(nspc).lt.0.0) then
          write(*,*) 'BIGSinisrc:ERROR: negative energy'
          stop
        else if (espc(nspc).lt.espc(max(nspc-1,1))) then
          write(*,*) 'BIGSinisrc:ERROR: decreasing energy'
          stop
        endif
        if (prob.lt.0.0) exit                ! End of spectrum
        if (nspc.ge.nemax) then
          write(*,*) 'BIGSinisrc:ERROR: too many bins in spectrum;'
          write(*,*) '              enlarge NEMAX'
          stop
        endif
        pspc(nspc+1) = pspc(nspc)+prob       ! Accumulated probability
      enddo
      if (spcfilen.ne.'-') close(spcfile)

      write(*,'(a)') 'No. of bins read:'
      write(*,'(1x,i0)') nspc-1
      if (nspc-1.lt.1) then
        write(*,*)
     &    'BIGSinisrc:ERROR: at least 1 bin must be defined'
        stop
      endif
      write(*,'(a)')
     &  'Sum of relative probabilities before normalization:'
      write(*,'(1x,es12.5)') pspc(nspc)
      if (pspc(nspc).le.0.0) then
        write(*,'(a)')
     &    'BIGSinisrc:ERROR: all probabilities are zero.'
        stop
      endif
      do j=1,nspc    ! Normalize to unity
        pspc(j) = pspc(j)/pspc(nspc)
      enddo
      do j=1,nspc-1    ! Prepare auxiliar array for sampling
        despc(j) = 0.0
        if (pspc(j+1).gt.pspc(j))
     &    despc(j) = (espc(j+1)-espc(j))/(pspc(j+1)-pspc(j))
      enddo

      write(*,'(a)') 'FWHM of Gaussian energy distrib. (eV):'
      read(*,*) esigma
      write(*,'(1x,es12.5)') esigma
      esigma = esigma*fwhm2sig
      emax = espc(nspc)+kfact*esigma  ! Set max energy
      if (parsrc.eq.3) emax = emax+twomc2  ! Allowance for e+ annihilation

      ! End of section:
      read(*,'(a80)') buffer
      if (index(buffer,eos).eq.0) then
        write(*,*) 'BIGSinisrc:ERROR: End-Of-Section mark not found'
        write(*,'(a,a)') '  expecting to find: ',eos
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif
      write(*,'(a)') '>>>> BIGS source initialization finished >>>>'
      end


!>>>> End Of File >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
