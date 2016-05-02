!*******************************************************************
!*                         SOURCE                                  *
!*                       PHASE SPACE                               *
!*                                                                 *
!* Short description:                                              *
!*   Generation of primary particle states for radiation transport *
!*   calculations with PENELOPE.                                   *
!*                                                                 *
!*   Initial particle states are read from an external phase space *
!*   file (PSF). Notice that, due to the modification of history   *
!*   numbers according to the information in the PSF, this source  *
!*   is *incompatible* with other source models.                   *
!*******************************************************************


      subroutine PSFSource
!*******************************************************************
!*    Output:                                                      *
!*      through secondary stack                                    *
!*      nhist -> top history counter (via module ctrsimMod)        *
!*******************************************************************
      use TRACK_mod
      use ctrsimMod, only : nhist
      implicit none

      integer*4 kpars,ilbs,dns
      real*8 es,xs,ys,zs,us,vs,ws,wghts
      common /srcpsf/ es,xs,ys,zs,us,vs,ws,wghts,kpars,ilbs(5),dns

      logical active,formatiaea
      integer ufile,split
      real*8 rot,xshift,yshift,zshift
      common /srcps1/ rot(3,3),xshift,yshift,zshift,split,ufile,
     &                formatiaea,active

      integer*8 npart,restartnpart,usednpart
      integer*8 nline,restartnline,usednline
      common /srcps2/ npart,restartnpart,usednpart,
     &                nline,restartnline,usednline

      logical isXtraStored(7)
      integer xtraPos(7)
      integer*4 sourceRead
      integer(selected_int_kind(R=18)) :: npiaeaMax,nTopiaea
      common /srciaeaps/ npiaeaMax,nTopiaea,sourceRead,
     &                   xtraPos,isXtraStored

#ifdef IAEAPSF
      logical getparIAEA
#endif

      logical getpar,ispart
      integer i
      integer*4 ncross
      real*8 infty,dsef,norm
      parameter (infty=1.0d30)

      if (.not.active) return

      if (isXtraStored(1)) then         ! IAEA PSF contains DeltaN (incremental history number)
        nhist = nhist+dble(dns-1)       ! Corrects history counter passed from main program
      else
        nhist = anint(nTopiaea*dble(npart)/dble(npiaeaMax)) ! Increment history no. proportionally to no. of particles read
        if (nhist.lt.1.0) nhist = 1.0   ! Avoid N=0 when only a few particles simulated
      endif

      do                                ! For all particles in the same history (while DNS=0)
        ! Load previously stored particle into active common:
        kpar = kpars
        wght = wghts                    ! Could be used by tally(7,.) below; if splitting applied, see further below
        e = es
        ! Rotate and translate position:
        x = xshift+rot(1,1)*xs+rot(1,2)*ys+rot(1,3)*zs
        y = yshift+rot(2,1)*xs+rot(2,2)*ys+rot(2,3)*zs
        z = zshift+rot(3,1)*xs+rot(3,2)*ys+rot(3,3)*zs
        ! Rotate direction and renormalize to double precision:
        u = rot(1,1)*us+rot(1,2)*vs+rot(1,3)*ws
        v = rot(2,1)*us+rot(2,2)*vs+rot(2,3)*ws
        w = rot(3,1)*us+rot(3,2)*vs+rot(3,3)*ws
        norm = 1.0/sqrt(u**2+v**2+w**2)
        u = u*norm
        v = v*norm
        w = w*norm
        ilb = ilbs
        call locate                                  ! Set quadric body and material
        call locatevox                               ! Set voxel indices and material, if in voxels
        if (mat.eq.0) then
          call stepx(infty,dsef,ncross)              ! Vacuum, flight up to object or infinity
          call tally(7,dsef)                         ! Inform tallies about this maiden flight
        endif
        wght = wghts/split                           ! Reduced weight for split particles
        do i=1,split                                 ! Particle splitting
          call stores(e,x,y,z,u,v,w,wght,kpar,ilb,0) ! Store copies in stack
        enddo
        call tally(0,e*split)                        ! Inform that particle(s) has been pushed to stack

        usednpart = npart                            ! Set no. of used PSF particles and lines
        usednline = nline

        ! Read a new particle (stored in /srcpsf/ for later if not used now):
        if (.not.formatiaea) then
          ispart = getpar()
#ifdef IAEAPSF
        else
          ispart = getparIAEA()
          if (.not.ispart) nhist = anint(dble(nTopiaea)) ! PSF exhausted, reset histories to value in header
#endif
        endif

        if (.not.ispart) then
          call simstop                               ! Force simulation to stop
          write(*,*) ''
          write(*,'(a)')
     &      'PSFsource: PSF exhausted, simulation stopped forcefully'
          return
        endif

        if (dns.ne.0) exit                           ! Exit the loop when primary history changes
      enddo
      end


      subroutine PSFSreport(mode)
!*******************************************************************
!*    Reports source data.                                         *
!*                                                                 *
!*    Input:                                                       *
!*      mode -> -1 if called at end-of-simulation, 0 for dump only,*
!*               1 for tally report only, 2 for both.              *
!*******************************************************************
      implicit none
      integer mode

      logical active,formatiaea
      integer ufile,split
      real*8 rot,xshift,yshift,zshift
      common /srcps1/ rot(3,3),xshift,yshift,zshift,split,ufile,
     &                formatiaea,active

      logical isXtraStored(7)
      integer xtraPos(7)
      integer*4 sourceRead
      integer(selected_int_kind(R=18)) :: npiaeaMax,nTopiaea
      common /srciaeaps/ npiaeaMax,nTopiaea,sourceRead,
     &                   xtraPos,isXtraStored

#ifdef IAEAPSF
      integer*4 errorcode
#endif

      if (.not.active) return
      if (mode.ne.1) call PSFSdump(1)  ! Sim dump file

#ifdef IAEAPSF
      if (mode.eq.-1.and.formatiaea) then  ! At EOS
        call iaea_destroy_source(sourceRead,errorcode)
        if (errorcode.lt.0) then
          write(*,'(a)') 'PSFSreport:ERROR: unable to close PSF.'
          ! Already at end-of-sim, so no need to stop
        endif
      endif
#endif
      end


      subroutine PSFSinisrc(activated,emax)
!*******************************************************************
!*    Initializes. To be called before SOURCE.                     *
!*                                                                 *
!*    Output:                                                      *
!*      activated -> TRUE if the source is active.                 *
!*      emax -> max source energy (eV)                             *
!*******************************************************************
      use ctrsimMod, only : nhist
      implicit none
      logical activated
      real*8 emax

      integer*4 kpars,ilbs,dns
      real*8 es,xs,ys,zs,us,vs,ws,wghts
      common /srcpsf/ es,xs,ys,zs,us,vs,ws,wghts,kpars,ilbs(5),dns

      logical active,formatiaea
      integer ufile,split
      real*8 rot,xshift,yshift,zshift
      common /srcps1/ rot(3,3),xshift,yshift,zshift,split,ufile,
     &                formatiaea,active

      integer*8 npart,restartnpart,usednpart
      integer*8 nline,restartnline,usednline
      common /srcps2/ npart,restartnpart,usednpart,
     &                nline,restartnline,usednline

      logical isXtraStored(7)
      integer xtraPos(7)
      integer*4 sourceRead
      integer(selected_int_kind(R=18)) :: npiaeaMax,nTopiaea
      common /srciaeaps/ npiaeaMax,nTopiaea,sourceRead,
     &                   xtraPos,isXtraStored

      logical getpar
      character*80 psfname,buffer,upcase
      character*(*) secid,eos
      parameter (secid=
     &'[SECTION SOURCE PHASE SPACE FILE v.2009-06-15]')
      parameter (eos='[END OF SPSF SECTION]')
      integer finduf,error,validate,formatpsf
      integer*8 npar(4),ntop,nhistpsf
      real*8 omega,theta,phi,pi,deg2rad,mc2,twomc2
      parameter (pi=3.1415926535897932d0,deg2rad=pi/180.0)
      parameter (mc2=5.10998918d5,twomc2=2.0*mc2)
      real*8, parameter:: nmax=1.0d15  ! Max number of allowed histories, to prevent counter overflow

      write(*,*) ''
      write(*,'(a)')
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      call getline(buffer,0)
      if (index(buffer,secid).eq.0) then
        write(*,'(a)') 'PSFinisrc:ERROR: incorrect section header;'
        write(*,'(a,a)') '  expecting to find: ',secid
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif
      write(*,'(a)') secid

      read(*,'(a80)') buffer
      buffer = adjustl(upcase(buffer))
      buffer = buffer(1:scan(buffer,' ')) ! Clip at 1st blank
      if (buffer.eq.'ON') then
        active = .true.
        activated = active
      else if (buffer.eq.'OFF') then
        active = .false.
        activated = active
        write(*, '(a)')
     &    '>>>> Source Phase Space File is OFF >>>>'
        do
          read(*,'(a80)',iostat=error) buffer
          if (error.ne.0) then
            write(*,'(a,a,a)') 'PSFinisrc:ERROR: ',
     &       'Unable to find End-Of-Section mark: ',eos
            stop
          endif
          if (index(buffer,eos).ne.0) return
        enddo
      else
        write(*,'(a)')
     &    'PSFinisrc:ERROR: expecting to find ON or OFF'
        write(*,'(a)') 'found instead:'
        write(*,'(a)') buffer
        stop
      endif

      read(*,*) formatpsf
      if (formatpsf.eq.0) then
        formatiaea = .false.
        write(*,'(a)') 'PSF format: standard penEasy in ASCII.'
        isXtraStored(1) = .true.  ! Informs that DeltaN is included in PSF
      else if (formatpsf.eq.1) then
        formatiaea = .true.
#ifdef IAEAPSF
        write(*,'(a)') 'PSF format: IAEA, in binary.'
        isXtraStored(1) = .false. ! To be set later
#else
        write(*,'(a)')
     &  'PSFinisrc:ERROR: IAEA PSF format requested but not available.'
        write(*,'(a)')
     &    '  Read the instructions provided in the documentation '//
     &    'to activate this option.'
        stop
#endif
      else
        write(*,'(a)') 'PSFinisrc:ERROR: PSF format must be 0 or 1.'
        stop
      endif

      write(*,'(a)') 'PSF filename:'
      read(*,'(a80)') psfname
      psfname = adjustl(psfname)
      psfname = psfname(1:scan(psfname,' ')) ! Clip at 1st blank
      write(*,'(1x,a)') psfname

      write(*,'(a)') 'Splitting factor:'
      read(*,*) split
      write(*,'(1x,i0)') split
      if (split.lt.1) then
        write(*,'(a)') 'PSFinisrc:ERROR: split < 1'
        stop
      endif

      write(*,'(a)')
     & 'Euler angles (deg) to rotate position and direction vectors:'
      read(*,*) omega,theta,phi
      write(*,'(3(1x,es12.5))') omega,theta,phi
      omega = omega*deg2rad
      theta = theta*deg2rad
      phi = phi*deg2rad
      call rotation(omega,theta,phi,rot)  ! Calculate rotation matrix

      write(*,'(a)') 'Cartesian components of position shift (cm):'
      read(*,*) xshift,yshift,zshift
      write(*,'(3(1x,es12.5))') xshift,yshift,zshift

      read(*,*) validate
      if (validate.ne.1.and.validate.ne.0) then
        write(*,'(a)') 'PSFinisrc:ERROR: '//
     &    'invalid entry. VALIDATE field must be 0 or 1.'
        stop
      endif

      if (.not.formatiaea) then    ! Using native penEasy PSF format
        call checkFormat(psfname)  ! Find if it is 2008-compliant
        if (validate.eq.1) then    ! Pre-read PSF, validate and statistics
          write(*,'(a)') 'Starting PSF validation...'
          ufile = finduf()
          open(ufile,file=psfname,status='old',iostat=error)
          if (error.ne.0) then
            write(*,'(a)') 'PSFinisrc:ERROR: cannot open the PSF'
            stop
          endif
          nline = 0
          npart = 0

          emax = 0.0
          ntop = 0
          npar = 0
          do
            if (.not.getpar()) exit  ! EOF reached
            ntop = ntop+dns          ! Count particles in PSF, note that NLINE is already increased by GETPAR
            if (ntop.gt.nmax) then   ! Max num. histories set to prevent counter overflow
              write(*,'(a,es12.5)')
     &          'PSFinisrc:ERROR: No. of histories in PSF exceeds ',
     &          nmax
              stop
            endif
            if (kpars.lt.1.or.kpars.gt.4) then
              write(*,'(a)') 'PSFinisrc:ERROR: invalid KPAR: at line:'
              write(*,'(2(1x,i0))') kpars,nline
              stop
            endif
#ifndef PROTONS
            if (kpars.eq.4) then
              write(*,'(a)') 'PSFinisrc:ERROR: invalid KPAR: at line:'
              write(*,'(2(1x,i0))') kpars,nline
              stop
            endif
#endif
            npar(kpars) = npar(kpars)+1  ! Count particles of each type

            if (es.lt.0.0.or.es.gt.1.0d9) then
              write(*,'(a)')
     &          'PSFinisrc:ERROR: invalid energy(eV): at line:'
              write(*,'(es12.5,1x,i0)') es,nline
              stop
            endif
            emax = max(emax,es)
            if (us**2+vs**2+ws**2.lt.1.0d-30) then
              write(*,'(a)')
     &          'PSFinisrc:ERROR: null vector direction found at line:'
              write(*,'(1x,i0)') nline
              stop
            endif
          enddo
          close(ufile)

          write(*,'(a)') 'PSF validation finished.'
          write(*,'(a)') 'PSF statistics:'
          write(*,'(a)') '  No. electrons:'
          write(*,'(3x,i0)') npar(1)
          write(*,'(a)') '  No. photons:'
          write(*,'(3x,i0)') npar(2)
          write(*,'(a)') '  No. positrons:'
          write(*,'(3x,i0)') npar(3)
#ifdef PROTONS
          write(*,'(a)') '  No. protons:'
          write(*,'(3x,i0)') npar(4)
#endif
          write(*,'(a)') '  No. particles, total:'
          write(*,'(3x,i0)') sum(npar)
          write(*,'(a)') '  No. top primary histories:'
          write(*,'(3x,i0)') ntop
          write(*,'(a)') '  Max energy(eV):'
          write(*,'(3x,es12.5)') emax
          write(*,'(a)') ' (max energy declared in input file ignored)'
          if (npar(3).gt.0) emax = emax+twomc2  ! Allow for e+ annihilation
          read(*,'(a80)') buffer             ! Dummy input line with max E (unused if validate=1)

        else  ! Do not validate the PSF
          write(*,'(a)') 'User opted not to pre-validate the PSF.'
          write(*,'(a)') 'Max energy (eV) taken from input file:'
          read(*,*) emax
          write(*,'(1x,es12.5)') emax
        endif

        ! Prepare for 1st call to SOURCE:
        ufile = finduf()
        open(ufile,file=psfname,status='old',iostat=error)
        if (error.ne.0) then
          write(*,'(a)') 'PSFinisrc:ERROR: cannot open the PSF'
          stop
        endif
        nline = 0
        npart = 0
        usednpart = 0
        usednline = 0

        call PSFSdump(0)           ! Read restart values from restart file, if present; 0 else
        if (restartnpart.ne.0) then
          write(*,'(a)') 'Restarting simulation from dump file; '//
     &                   'no. of particles and file lines already '//
     &                   'read in previous simulation:'
          write(*,'(2(1x,i0))') restartnpart,restartnline
        endif

        nhistpsf = 0               ! Reset the number of histories skipped
        do                         ! Read particles until current particle history exceeds restart particle:
          if (.not.getpar()) then  ! Read a particle
            write(*,'(a)') 'PSFinisrc:ERROR: PSF end-of-file '//
     &        ' reached, not enough particles to initialize.'
            stop
          endif
          nhistpsf = nhistpsf+dns         ! Keep track of histories in PSF for double check
          if (npart.gt.restartnpart) then ! Keep skipping until npart > restart_npart
            if (nhistpsf.lt.nhist+0.5d0) then ! Check for inconsistency when Npsf.le.N (+0.5 for precision)
              write(*,'(a)') 'PSFinisrc:ERROR: Inconsistency '//
     &                       'found in PSF incremental history no.'
              write(*,'(a)') 'Current particle and line nos.:'
              write(*,'(2(1x,i0))') npart,nline
              write(*,'(a)') 'Current history no. N according to '//
     &          'sum(DeltaN) and last N simulated according to '//
     &          'dump file:'
              write(*,'(1x,i0,1x,f18.0)') nhistpsf,nhist+0.1d0
              write(*,'(a)') '  (the former should be larger '//
     &                       'than the latter)'
              stop
            endif
            exit
          endif
        enddo

#ifdef IAEAPSF
      else                              ! IAEA format
        call IAEAiniread(psfname,validate,emax)
        read(*,'(a80)') buffer          ! Dummy input line with max E (unused if IAEA format)
                                        ! The case formatiaea=1 and IAEAPSF not defined has already been blocked above
#endif
      endif

      read(*,'(a80)') buffer
      if (index(buffer,eos).eq.0) then
        write(*,'(a)')
     &    'PSFinisrc:ERROR: End-Of-Section mark not found'
        write(*,'(a,a)') '  expecting to find: ',eos
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif
      write(*,'(a)') '>>>> PSF source initialization finished >>>>'
      end


      subroutine PSFSdump(mode)
!*******************************************************************
!*    Dumps into or reads data from a dump file.                   *
!*                                                                 *
!*    Input:                                                       *
!*      mode -> 1 to write dump file, else to read from it.        *
!*******************************************************************
      use dumpMod
      implicit none
      integer mode

      logical active,formatiaea
      integer ufile,split
      real*8 rot,xshift,yshift,zshift
      common /srcps1/ rot(3,3),xshift,yshift,zshift,split,ufile,
     &                formatiaea,active

      integer*8 npart,restartnpart,usednpart
      integer*8 nline,restartnline,usednline
      common /srcps2/ npart,restartnpart,usednpart,
     &                nline,restartnline,usednline

      if (mode.eq.1) then
        if (dumpfile.lt.0) return     ! No dump file open
        write(dumpfile) usednpart,usednline
      else
        restartnline = 0
        restartnpart = 0
        if (restartfile.lt.0) return  ! No restart file open
        read(restartfile) restartnpart,restartnline
      endif
      end


      logical function getpar()
!*******************************************************************
!*    Reads a new particle from the PSF.                           *
!*                                                                 *
!*    Output:                                                      *
!*      -> returns .false. if EOF has been reached, else .true.    *
!*      -> particle state in /srcpsf/                              *
!*******************************************************************
      implicit none
      integer*4 kpars,ilbs,dns
      real*8 es,xs,ys,zs,us,vs,ws,wghts
      common /srcpsf/ es,xs,ys,zs,us,vs,ws,wghts,kpars,ilbs(5),dns

      logical active,formatiaea
      integer ufile,split
      real*8 rot,xshift,yshift,zshift
      common /srcps1/ rot(3,3),xshift,yshift,zshift,split,ufile,
     &                formatiaea,active

      integer*8 npart,restartnpart,usednpart
      integer*8 nline,restartnline,usednline
      common /srcps2/ npart,restartnpart,usednpart,
     &                nline,restartnline,usednline

      logical format2008
      common /formatVer/ format2008

      character*256 buffer
      integer error

      do
        read(ufile,'(a256)',end=10,iostat=error) buffer
        nline = nline+1               ! Lines read
        if (error.ne.0) then
          write(*,'(a)') 'getpar:ERROR: unable to read PSF line no.:'
          write(*,'(i0)') nline
          stop
        endif
        if (buffer(1:1).ne.'#') exit  ! A non-comment line was found
      enddo
      npart = npart+1                 ! Particles read

      if (format2008) then            ! 2008-compliant format
        read(buffer,*,iostat=error)
     &    kpars,es,xs,ys,zs,us,vs,ws,wghts,dns,
     &    ilbs(1),ilbs(2),ilbs(3),ilbs(4),ilbs(5)
      else                            ! Pre-2008 format
        read(buffer,*,iostat=error)
     &    kpars,es,xs,ys,zs,us,vs,ws,wghts,ilbs(5),dns
        ilbs(1) = 1                   ! Assume all particles are primaries
        ilbs(2) = 0
        ilbs(3) = 0
        ilbs(4) = 0
      endif
      if (error.ne.0) then
        write(*,'(a)')
     &    'getpar:ERROR: invalid or missing datum in PSF line:'
        write(*,'(i0)') nline
        write(*,'(a)') '  line contents:'
        write(*,'(a)') buffer
        stop
      endif
      getpar = .true.
      return

 10   getpar = .false.                ! EOF
      end


      subroutine checkFormat(psfname)
!*******************************************************************
!*    Checks if PSF data format is 2008-compliant.                 *
!*                                                                 *
!*    Input:                                                       *
!*      psfname -> PSF filename.                                   *
!*    Output:                                                      *
!*      -> Sets /formatVer/.                                       *
!*    Comments:                                                    *
!*      The file is opened and closed.                             *
!*******************************************************************
      implicit none
      character*80 psfname

      logical active,formatiaea
      integer ufile,split
      real*8 rot,xshift,yshift,zshift
      common /srcps1/ rot(3,3),xshift,yshift,zshift,split,ufile,
     &                formatiaea,active

      logical format2008
      common /formatVer/ format2008

      character*256 buffer,id2008,id2006
      parameter (id2008=
     & '# [PHASE SPACE FILE FORMAT penEasy v.2008-05-15]')
      parameter (id2006=
     & '# kpar : e : x : y : z : u : v : w : wght : ilb(5) : DeltaN')
      integer error,finduf

      ufile = finduf()
      open(ufile,file=psfname,status='old',iostat=error)
      if (error.ne.0) then
        write(*,'(a)') 'checkFormat:ERROR: cannot open the PSF'
        stop
      endif

      read(ufile,'(a256)',iostat=error) buffer
      if (error.ne.0) then
        write(*,'(a)')
     &    'ckeckFormat:ERROR: unable to read first PSF line.'
        stop
      endif

      if (buffer.eq.id2008) then
        format2008 = .true.
      else if (buffer.eq.id2006) then
        format2008 = .false.
        write(*,*) ''
        write(*,'(a)')
     &    '********************'
        write(*,'(a)')
     &    'checkFormat:WARNING: header for 2008-compliant format '//
     &    'not found;'
        write(*,'(a)')
     &    '  assuming pre-2008 data format, which consists of:'
        write(*,'(a)')
     &    '  KPAR:E:X:Y:Z:U:V:W:WGHT:ILB(5):DeltaN'
        write(*,'(a)')
     &    '********************'
        write(*,*) ''
      else
        format2008 = .false.
        write(*,*) ''
        write(*,'(a)')
     &  'checkFormat:ERROR: unable to identify PSF format;'
        write(*,'(a)')
     &  '  expecting one of these two headers in first line of PSF:'
        write(*,'(a)') id2008
        write(*,'(a)') id2006
        write(*,'(a)')
     &    '  (for 2008 and 2006 formats, respectively) '//
     &    'but neither was found.'
        stop
      endif

      close(ufile)
      end


!*******************************************************************
!*    IAEA specific routines follow.                               *
!*    Not adapted to protons: getparIAEA will issue an error       *
!*    message and abort if called.                                 *
!*******************************************************************
#ifdef IAEAPSF

      logical function getparIAEA()
!*******************************************************************
!*    Reads a new particle from the IAEA PSF.                      *
!*                                                                 *
!*    Output:                                                      *
!*      -> returns .false. if EOF has been reached, else .true.    *
!*      -> particle state in /srcpsf/                              *
!*******************************************************************
      implicit none
      integer*4 kpars,ilbs,dns
      real*8 es,xs,ys,zs,us,vs,ws,wghts
      common /srcpsf/ es,xs,ys,zs,us,vs,ws,wghts,kpars,ilbs(5),dns

      logical active,formatiaea
      integer ufile,split
      real*8 rot,xshift,yshift,zshift
      common /srcps1/ rot(3,3),xshift,yshift,zshift,split,ufile,
     &                formatiaea,active

      integer*8 npart,restartnpart,usednpart
      integer*8 nline,restartnline,usednline
      common /srcps2/ npart,restartnpart,usednpart,
     &                nline,restartnline,usednline

      logical isXtraStored(7)
      integer xtraPos(7)
      integer*4 sourceRead
      integer(selected_int_kind(R=18)) :: npiaeaMax,nTopiaea
      common /srciaeaps/ npiaeaMax,nTopiaea,sourceRead,
     &                   xtraPos,isXtraStored

      integer xtraMax
      parameter (xtraMax=10)  ! Up to 10 xtraInts or xtraReals
      integer*4 errorcode,kpariaea,xtraInt(xtraMax)
      real*4 eiaea,xiaea,yiaea,ziaea,uiaea,viaea,wiaea,wghtiaea
      real*4 xtraReal(xtraMax)

      getparIAEA = .false.
      if (npart.ge.npiaeaMax) return  ! EOF reached

      call iaea_get_particle(sourceRead,errorcode,kpariaea,eiaea,
     &  wghtiaea,xiaea,yiaea,ziaea,uiaea,viaea,wiaea,xtraReal,
     &  xtraInt)
      nline = nline+1  ! Lines read, irrelevant for IAEA format, but just in case
      npart = npart+1  ! Particles read

      select case (errorcode)
      case (-1)
        write(*,'(a,i0)') 'getparIAEA:ERROR: unable to read '//
     &    'particle no.: ',npart
        stop
      case (-2)
        write(*,'(a,i0)')
     &    'getparIAEA:ERROR: EOF reached when attempting to read '//
     &    'particle no.: ',npart
        stop
      case default
        if (errorcode.lt.0) then
          write(*,'(a,i0)')
     &     'getparIAEA:ERROR: undefined error while attempting '//
     &     'to read particle no.: ',npart
          stop
        endif
      end select

      select case (kpariaea)  ! Conversion from IAEA particle index to PENELOPE's
      case (1)                ! Photon
        kpars = 2
      case (2)                ! Electron
        kpars = 1
      case (3)                ! Positron
        kpars = 3
      case default            ! Protons or other particles not implemented, sorry
        write(*,'(a,i0,a,i0)') 'getparIAEA:ERROR: Invalid KPAR: ',
     &    kpariaea,' found at particle no.: ',npart
        stop
      end select

      es = dble(eiaea*1.0e6)  ! Conversion from MeV (IAEA) to eV (PENELOPE)
      xs = dble(xiaea)
      ys = dble(yiaea)
      zs = dble(ziaea)
      us = dble(uiaea)
      vs = dble(viaea)
      ws = dble(wiaea)
      wghts = dble(wghtiaea)

      ilbs(1) = 1             ! By default, assume it's a primary particle
      if (isXtraStored(7)) ilbs(1) = xtraInt(xtraPos(7))
      ilbs(2) = 0             ! By default, assume null value
      if (isXtraStored(6)) ilbs(2) = xtraInt(xtraPos(6))
      ilbs(3) = 0
      if (isXtraStored(5)) ilbs(3) = xtraInt(xtraPos(5))
      ilbs(4) = 0
      if (isXtraStored(4)) ilbs(4) = xtraInt(xtraPos(4))
      ilbs(5) = 0
      if (isXtraStored(3)) ilbs(5) = xtraInt(xtraPos(3))
      ! Note: IAEA xtra type with code=2 reserved for EGS LATCH variable

      if (isXtraStored(1)) then      ! Incremental history was stored in PSF
        dns = xtraInt(xtraPos(1))
      else
        dns = 0                      ! Same history
        if (errorcode.gt.0) dns = 1  ! Arbitrarily set to 1 to break loop in PSFsource()
      endif

      getparIAEA = .true.
      end


      subroutine IAEAiniread(psfnameF,validate,emax)
!*******************************************************************
!*    Initializes.                                                 *
!*                                                                 *
!*    Input:                                                       *
!*      psfnameF ->  PSF filename, Fortran convention              *
!*      validate ->  1 if validation is requested                  *
!*    Output:                                                      *
!*      emax ->  Maximum energy (eV) in the PSF                    *
!*******************************************************************
      use ctrsimMod, only : nhist
      implicit none
      character*80 psfnameF
      integer validate
      real*8 emax

      integer*4 kpars,ilbs,dns
      real*8 es,xs,ys,zs,us,vs,ws,wghts
      common /srcpsf/ es,xs,ys,zs,us,vs,ws,wghts,kpars,ilbs(5),dns

      logical active,formatiaea
      integer ufile,split
      real*8 rot,xshift,yshift,zshift
      common /srcps1/ rot(3,3),xshift,yshift,zshift,split,ufile,
     &                formatiaea,active

      integer*8 npart,restartnpart,usednpart
      integer*8 nline,restartnline,usednline
      common /srcps2/ npart,restartnpart,usednpart,
     &                nline,restartnline,usednline

      logical isXtraStored(7)
      integer xtraPos(7)
      integer*4 sourceRead
      integer(selected_int_kind(R=18)) :: npiaeaMax,nTopiaea
      common /srciaeaps/ npiaeaMax,nTopiaea,sourceRead,
     &                   xtraPos,isXtraStored

      logical getparIAEA,ivarcnt(7)
      character*81 psfname  ! PSF filename, C convention
      integer i,xtraMax
      parameter (xtraMax=10)  ! Up to 10 xtraInt or xtraReals
      integer*4 accessRead,errorcode,nXtraReal,nXtraInt
      parameter (accessRead=1)
      integer*4 xtraRealTypes(xtraMax),xtraIntTypes(xtraMax)
      integer*8 nTopiaeacheck,nhistpsf,nele,npho,npos
      real*4 emaxiaea,constant,varcnt(7)
      real*8 emaxcheck,mc2,twomc2,zero
      parameter (mc2=5.10998918d5,twomc2=2.0*mc2,zero=1.0d-30)
      real*8, parameter:: nmax=1.0d15  ! Max number of allowed histories, to prevent counter overflow

      write(*,'(a)') 'Initializing PSF in IAEA format.'

      psfname = psfnameF
      i = len_trim(psfnameF)+1
      psfname(i:i) = char(0)  ! End string according to C convention

      call iaea_new_source(sourceRead,psfname,accessRead,errorcode)
      if (errorcode.lt.0) then
        write(*,'(a)') 'IAEAiniread:ERROR: unable to open PSF.'
        stop
      endif
      nline = 0  ! Reset PSF pointers
      npart = 0

      call iaea_get_max_particles(sourceRead,-1,npiaeaMax)
      if (npiaeaMax.lt.0) then
        write(*,'(a)') 'IAEAiniread:ERROR: '//
     &    'Unable to get the total number of particles.'
        stop
      endif

      call iaea_get_total_original_particles(sourceRead,nTopiaea)
      if (nTopiaea.lt.0) then
        write(*,'(a)') 'IAEAiniread:ERROR: '//
     &    'Unable to get the total number of histories.'
        stop
      endif
      if (nTopiaea.gt.nmax) then   ! Max num. histories set to prevent counter overflow
        write(*,'(a,es12.5)')
     &    'IAEAiniread:ERROR: No. of histories in PSF exceeds ',
     &    nmax
        stop
      endif

      do i=1,7
        call iaea_get_constant_variable(sourceRead,i-1,constant,
     &                                  errorcode)
        select case (errorcode)
        case (-1)
          write(*,'(a)')
     &      'IAEAiniread:ERROR: Unable to get constant variables.'
          stop
        case(-2)
          write(*,'(a)')
     &      'IAEAiniread:ERROR: Variable constant index out of range.'
          stop
        case (-3)
          ivarcnt(i) = .false.
          varcnt(i) = 0.0
        case default
          ivarcnt(i) = .true.
          varcnt(i) = constant
        end select
      enddo

      if (all(.not.ivarcnt)) then
        write(*,'(a)') '  PSF does not contain constant variables.'
      else
        write(*,'(a)') '  Variables stored as a constant:'
      endif
      if (ivarcnt(1)) then
        write(*,'(2x,a,es12.5,a)') '    x = ',varcnt(1),' cm'
      else if (ivarcnt(2)) then
        write(*,'(2x,a,es12.5,a)') '    y = ',varcnt(2),' cm'
      else if (ivarcnt(3)) then
        write(*,'(2x,a,es12.5,a)') '    z = ',varcnt(3),' cm'
      else if (ivarcnt(4)) then
        write(*,'(2x,a,es12.5)') '    u = ',varcnt(4)
      else if (ivarcnt(5)) then
        write(*,'(2x,a,es12.5)') '    v = ',varcnt(5)
      else if (ivarcnt(6)) then
        write(*,'(2x,a,es12.5)') '    w = ',varcnt(6)
      else if (ivarcnt(7)) then
        write(*,'(2x,a,es12.5)') '    weight = ',varcnt(7)
      endif

      call iaea_get_extra_numbers(sourceRead,nXtraReal,nXtraInt)
      write(*,'(a,i0)') '  No. of extra integers stored: ',nXtraInt
      write(*,'(a,i0)') '  No. of extra reals stored: ',nXtraReal
      if (max(nXtraInt,nXtraReal).gt.xtraMax) then
        write(*,'(a,i0,a)')
     &    'IAEAiniread:ERROR: More than ',xtraMax,
     &    ' extra int or real variables in PSF; '//
     &    'increase array dims or notify authors.'
          stop
      endif

      call iaea_get_type_extra_variables(sourceRead,errorcode,
     &  xtraIntTypes,xtraRealTypes)
      do i=1,nXtraReal
        select case (xtraRealTypes(i))
        case (0)
          write(*,'(a)') '  PSF contains user-defined real variable.'
        case (1)
          write(*,'(a)') '  PSF contains EGSs XLAST variable.'
        case (2)
          write(*,'(a)') '  PSF contains EGSs YLAST variable.'
        case (3)
          write(*,'(a)') '  PSF contains EGSs ZLAST variable.'
        case default
          write(*,'(a)') '  PSF contains unidentified variable type.'
        end select
      enddo

      isXtraStored = .false. ! Init arrays
      xtraPos    = 0
      do i=1,nXtraInt        ! Find what type of var is stored for each xtra int
        select case (xtraIntTypes(i))
        case (0)
          write(*,'(a)') '  PSF contains user-defined int variable.'
        case (1)
          write(*,'(a)') '  PSF contains incremental history number.'
          isXtraStored(1) = .true.
          xtraPos(1) = i
        case (2)
          write(*,'(a)') '  PSF contains EGSs variable LATCH.'
          isXtraStored(2) = .true.  ! Unused by penEasy
          xtraPos(2) = i
        case (3)
          write(*,'(a)') '  PSF contains PENELOPEs variable ILB(5).'
          isXtraStored(3) = .true.
          xtraPos(3) = i
        case (4)
          write(*,'(a)') '  PSF contains PENELOPEs variable ILB(4).'
          isXtraStored(4) = .true.
          xtraPos(4) = i
        case (5)
          write(*,'(a)') '  PSF contains PENELOPEs variable ILB(3).'
          isXtraStored(5) = .true.
          xtraPos(5) = i
        case (6)
          write(*,'(a)') '  PSF contains PENELOPEs variable ILB(2).'
          isXtraStored(6) = .true.
          xtraPos(6) = i
        case (7)
          write(*,'(a)') '  PSF contains PENELOPEs variable ILB(1).'
          isXtraStored(7) = .true.
          xtraPos(7) = i
        case default
          write(*,'(a)') '  PSF contains unidentified int variable.'
        end select
      enddo

      if (validate.eq.1) then  ! Pre-read IAEA PSF, validate and statistics
        write(*,'(a)') 'Starting PSF validation...'

        call iaea_check_file_size_byte_order(sourceRead,errorcode)
        select case (errorcode)
        case (-1)
          write(*,'(a)') 'IAEAiniread:ERROR: '//
     &      'Unable to check file size and byte order of the PSF.'
          stop
        case (-2)
          write(*,'(a)') 'IAEAiniread:ERROR: '//
     &      'The function fseek fails.'
          stop
        case (-3)
          write(*,'(a)') 'IAEAiniread:ERROR: File size'//
     &      'inconsistent with header checksum.'
          stop
        case (-4)
          write(*,'(a)') 'IAEAiniread:ERROR: '//
     &      'There is a byte order mismatch.'
          stop
        case (-5)
          write(*,'(a)') 'IAEAiniread:ERROR: '//
     &      'There is a file size and byte order mismatch.'
          stop
        case default
          if (errorcode.lt.0) then
            write(*,'(a)') 'IAEAiniread:ERROR: '//
     &        'Unidentified error code.'
            stop
          endif
        end select

        ! Reset:
        nele = 0
        npho = 0
        npos = 0
        nTopiaeacheck = 0
        emaxcheck = 0.0
        do                             ! For the whole PSF
          if (.not.getparIAEA()) exit  ! Get a particle and check EOF

          nTopiaeacheck = nTopiaeacheck+dns

          select case (kpars)
          case (1)
            nele = nele+1
          case (2)
            npho = npho+1
          case (3)
            npos = npos+1
          case default
            write(*,'(a)') 'IAEAiniread:internalERROR: invalid KPAR.'
            stop
          end select

          if (kpars.eq.3) es = es+twomc2  ! For e+ add 2*511 keV
          emaxcheck = max(emaxcheck,es)   ! Compute max E
          if (es.lt.0.0.or.es.gt.1.0e9) then
            write(*,'(a,i0)')
     &        'IAEAiniread:ERROR: energy out of range found '//
     &        'in particle no.: ',npart
            stop
          endif

          if (us**2+vs**2+ws**2.lt.zero) then
            write(*,'(a)')
     &        'IAEAiniread:ERROR: null direction vector found '//
     &        'in particle no.: ',npart
            stop
          endif
        enddo

        write(*,'(a)') 'Validation of PSF in IAEA format finished.'
        write(*,'(a)') 'PSF statistics:'
        write(*,'(a)') '  No. electrons:'
        write(*,'(3x,i0)') nele
        write(*,'(a)') '  No. photons:'
        write(*,'(3x,i0)') npho
        write(*,'(a)') '  No. positrons:'
        write(*,'(3x,i0)') npos
        write(*,'(a)') '  No. particles, total:'
        write(*,'(3x,i0)') nele+npho+npos
        write(*,'(a)') '  No. particles declared in header:'
        write(*,'(3x,i0)') npiaeaMax
        write(*,'(a)') '  No. primary histories declared in header:'
        write(*,'(3x,i0)') nTopiaea
        if (npiaeaMax.ne.nele+npho+npos) then
          write(*,'(a)')
     &      'IAEAiniread:ERROR: No. of particles inferred '//
     &      'from PSF differs from that stated in header file.'
            stop
        endif
        if (nTopiaeacheck.gt.nTopiaea) then
          write(*,'(a)')
     &      'IAEAiniread:ERROR: No. of histories inferred '//
     &      'from PSF is greater than that stated in header file.'
            stop
        endif

      else  ! Do not validate the IAEA PSF
        write(*,'(a)') 'User opted not to pre-validate IAEA PSF.'
        emaxcheck = 0.0  ! Dummy value to make header max E prevail
      endif

      call iaea_get_maximum_energy(sourceRead,emaxiaea)
      if (emaxiaea.lt.0.0) then
        write(*,'(a)')
     &    'IAEAiniread:ERROR: Unable to get maximum energy.'
        stop
      endif
      emax = dble(emaxiaea)*1.0d6  ! Conversion from MeV to eV

      if (emaxcheck.gt.emax) then
        write(*,'(a)')
     &   'IAEAiniread:WARNING: Max energy reported by IAEA '//
     &   'routines is less than max value in PSF.'
        write(*,'(a)') '  Emax_IAEA : Emax_PSF (eV):'
        write(*,'(2(1x,es12.5))') emax,emaxcheck
        write(*,'(a)')
     &   '  Taking Emax directly from PSF and continuing.'
        emax = emaxcheck
      endif
      write(*,'(a)') 'Maximum energy in PSF (eV):'
      write(*,'(1x,es12.5)') emax

      ! Prepare for 1st call to SOURCE:
      call iaea_destroy_source(sourceRead,errorcode)                ! Close PSF to destroy allocation
      if (errorcode.lt.0) then
        write(*,'(a)') 'IAEAiniread:ERROR: unable to close PSF.'
        stop
      endif
      call iaea_new_source(sourceRead,psfname,accessRead,errorcode) ! ...and then re-open it.
      if (errorcode.lt.0) then
        write(*,'(a)') 'IAEAiniread:ERROR: unable to open PSF.'
        stop
      endif
      nline = 0  ! Reset pointers
      npart = 0

      call PSFSdump(0)           ! Read restart values from restart file, if present; 0 else
      if (restartnpart.ne.0) then
        write(*,'(a)') 'Restarting simulation from dump file; '//
     &                 'no. of particles already '//
     &                 'read in previous simulation:'
        write(*,'(1x,i0)') restartnpart
      endif

      ! Read particles until current particle exceeds restart particle
      nhistpsf = 0                      ! Reset the number of histories skipped
      do                                ! Read particles until current particle exceeds restart particle
        if (.not.getparIAEA()) then     ! Read a particle
          write(*,'(a)') 'IAEAiniread:ERROR: PSF end-of-file '//
     &      'reached, not enough particles to initialize.'
          stop
        endif
        nhistpsf = nhistpsf+dns         ! Keep track of histories in PSF for double check
        if (npart.gt.restartnpart) then ! Break loop when npart > restart_particle
          if (isXtraStored(1).and.nhistpsf.lt.nhist+0.5d0) then ! Check for inconsistency when Npsf.le.N (+0.5 for precision)
            write(*,'(a)') 'IAEAiniread:ERROR: Inconsistency '//
     &                       'found in PSF incremental history no.'
            write(*,'(a)') '  Current particle no.:'
            write(*,'(1x,i0)') npart
            write(*,'(a)') '  Current history no. N according to '//
     &        'sum(DeltaN) and last N simulated according to '//
     &        'dump file:'
            write(*,'(1x,i0,1x,f18.0)') nhistpsf,nhist+0.1d0  ! 0.1 added for truncation precision
            write(*,'(a)') '  (the former should be larger '//
     &                     'than the latter)'
            stop
          endif
          exit                          ! Break loop
        endif
      enddo

      write(*,'(a)') 'IAEAiniread: Initialization finished.'
      end

#endif
!*******************************************************************
!*    End of IAEA block.                                           *
!*******************************************************************


!>>>> End Of File >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
