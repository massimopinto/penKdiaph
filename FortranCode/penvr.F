!*******************************************************************
!*                           PENVR                                 *
!*                                                                 *
!* Short description:                                              *
!*   Variance reduction techniques for penEasy. Note that the      *
!*   original routines for interaction forcing are included in file*
!*   penpatch.f.                                                   *
!* Rev:                                                            *
!*   2015-05-30                                                    *
!*******************************************************************


!*******************************************************************
!*******************************************************************
!*    Modules section                                              *
!*******************************************************************
!*******************************************************************

      module forcingMod
!*******************************************************************
!*    Interaction forcing vars                                     *
!*******************************************************************
      use PENELOPE_mod, only : maxmat
      implicit none
      save
      logical analog(maxmat,4),isforcing,active
      real*8 minwght
      end


      module splitMod
!*******************************************************************
!*    Particle splitting vars                                      *
!*******************************************************************
      implicit none
      save
      logical active
      integer signsplit,mode,nsplit
      integer*4 matsplit
      real*8 rotS(3,3),rot(3,3),rotz(3),rorx(3,3),rory(3,3),rorxy(3,3)
      real*8 dr(3),phi0,dphi,idphi,minwght
      end


      module russiaMod
!*******************************************************************
!*    Russian roulette vars                                        *
!*******************************************************************
      implicit none
      save
      logical active
      integer*4 matrussia
      real*8 maxwght,psurv
      end


!*******************************************************************
!*******************************************************************
!*    Routines section                                             *
!*******************************************************************
!*******************************************************************

      subroutine iniforce(e0)
!*******************************************************************
!*    Initializes the interaction forcing routines                 *
!*                                                                 *
!*    Input:                                                       *
!*      e0 -> source energy (eV)                                   *
!*******************************************************************
      use forcingMod
      use PENELOPE_mod, only : maxmat
      use PENGEOM_mod, only : mater,nbody
      use PENVARED_mod, only : force
      implicit none
      real*8 e0

      character*(*) secid,eos
      parameter (secid=
     &'[SECTION INTERACTION FORCING v.2009-06-15]')
      parameter (eos='[END OF VRIF SECTION]')
      character*80 buffer
      integer i,j,mat,kpar,icol,error
      real*8 forcing,phmfp,hmfp

      write(*,*) ''
      write(*,'(a)')
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      call getline(buffer,0)
      if (index(buffer,secid).eq.0) then
        write(*,'(a)') 'iniforce:ERROR: incorrect section header;'
        write(*,'(a,a)') '  expecting to find: ',secid
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif
      write(*,'(a)') secid

      read(*,'(a80)') buffer
      buffer = adjustl(buffer)
      buffer = buffer(1:scan(buffer,' ')) ! Clip at 1st blank
      if (buffer.eq.'ON') then
        active = .true.
        write(*, '(a)')
     &    '>>>> Interaction Forcing is ON >>>>'
      else if (buffer.eq.'OFF') then
        active = .false.
        write(*, '(a)')
     &    '>>>> Interaction Forcing is OFF >>>>'
        do
          read(*,'(a80)',iostat=error) buffer
          if (error.ne.0) then
            write(*,'(a,a,a)') 'iniforce:ERROR: ',
     &       'Unable to find End-Of-Section mark: ',eos
            stop
          endif
          if (index(buffer,eos).ne.0) return
        enddo
      else
        write(*,'(a)')
     &    'iniforce:ERROR: expecting to find ON or OFF'
        write(*,'(a)') 'found instead:'
        write(*,'(a)') buffer
        stop
      endif

      write(*,'(a)') 'Min weight to apply:'
      read(*,*) minwght
      write(*,'(1x,es12.5)') minwght

      ! Reset all:
      force = 1.0
      analog = .true.

      write(*,'(a)') 'MAT : KPAR : ICOL : forcing'
      read(*,'(a80)') buffer  ! Table header
      do
        read(*,*,iostat=error) mat,kpar,icol,forcing
        if (error.ne.0) then
          write(*,'(a)')
     &     'iniforce:ERROR: unable to read line containing:'
          write(*,'(a)')
     &     '  MAT,KPAR,ICOL,forcing'
          write(*,'(a)')
     &     '  Recall to use a line with MAT=-1 to end the list.'
          stop
        endif
        if (mat.eq.0.or.mat.eq.-1) exit  ! End of list
        write(*,'(3(1x,i3),1x,es12.5)') mat,kpar,icol,forcing
        if (mat.lt.0.or.mat.gt.maxmat) then
          write(*,'(a)') 'iniforce:ERROR: invalid MAT'
          stop
        endif
        if (kpar.eq.4) then ! Temporary until proton forcing implemented
          write(*,'(a)')
     &      'iniforce:ERROR: proton forcing not implemented yet'
          stop
        endif
        if (kpar.lt.1.or.kpar.gt.3) then
          write(*,'(a)') 'iniforce:ERROR: KPAR must be in [1,3]'
          stop
        endif
        if (icol.lt.0.or.icol.gt.8) then
          write(*,'(a)') 'iniforce:ERROR: ICOL must be in [0,8]'
          stop
        endif
        if (forcing.lt.1.0) then
          write(*,'(a)') 'iniforce:ERROR: FORCING must not be < 1'
          stop
        endif
        do i=1,nbody  ! Sets forcing for all bodies of the selected material
          if (mater(i).eq.mat) then
            if (icol.ne.0) then
              force(i,kpar,icol) = forcing
            else  ! Force all interactions
              do j=1,8
                force(i,kpar,j) = forcing
              enddo
            endif
          endif
        enddo
        analog(mat,kpar) = .false.
      enddo

      write(*,*) ''
      write(*,'(a,es12.5)')
     &  'INFO: Unforced (analog) hard mean free paths at E(eV) =',e0
      write(*,'(a)') 'MAT : KPAR : ICOL : HMFP(cm)'
      do mat=1,maxmat
        do kpar=1,3
          if (analog(mat,kpar)) cycle  ! Skip unforced interactions
          do icol=1,8
            if (kpar.eq.1.and.(icol.lt.2.or.icol.gt.5)) cycle
            if (kpar.eq.2.and.icol.gt.4) cycle
            if (kpar.eq.3.and.(icol.lt.2.or.icol.gt.6)) cycle
            hmfp = phmfp(e0,kpar,mat,icol)
            write(*,'(i3,1x,i1,1x,i1,1x,es12.5)')
     &        mat,kpar,icol,hmfp
          enddo
        enddo
      enddo
      write(*,*) ''
      write(*,'(a)') 'iniforce:WARNING:'//
     &  ' interaction forcing may bias pulse height spectra'

      read(*,'(a80)') buffer
      if (index(buffer,eos).eq.0) then
        write(*,'(a)')
     &    'iniforce:ERROR: End-Of-Section mark not found'
        write(*,'(a,a)') '  expecting to find: ',eos
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif

      write(*,'(a)')
     &  '>>>> Interaction forcing initialization finished >>>>'
      end


      subroutine forcing
!*******************************************************************
!*    Sets interaction forcing status.                             *
!*******************************************************************
      use forcingMod
      use TRACK_mod, only : kpar,mat,wght
      implicit none

      isforcing = .false.
      if (.not.active.or.mat.eq.0.or.wght.lt.minwght) return
      if (analog(mat,kpar)) return
      isforcing = .true.
      end


      subroutine inisplit
!*******************************************************************
!*    Initializes the interaction forcing routines                 *
!*******************************************************************
      use TRACK_mod
      use splitMod
      implicit none
      character*(*) secid,eos
      parameter (secid=
     &'[SECTION SPLITTING v.2015-05-30]')
      parameter (eos='[END OF VRS SECTION]')
      character*80 buffer
      integer error,i
      real*8 omega,theta,phi
      real*8 pi,deg2rad,twopi
      parameter (pi=3.1415926535897932d0,twopi=2.0d0*pi)
      parameter (deg2rad=pi/180.0d0)

      write(*,*) ''
      write(*,'(a)')
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      call getline(buffer,0)
      if (index(buffer,secid).eq.0) then
        write(*,'(a)') 'inisplit:ERROR: incorrect section header;'
        write(*,'(a,a)') '  expecting to find: ',secid
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif
      write(*,'(a)') secid

      read(*,'(a80)') buffer
      buffer = adjustl(buffer)
      buffer = buffer(1:scan(buffer,' ')) ! Clip at 1st blank
      if (buffer.eq.'ON') then
        active = .true.
        write(*, '(a)')
     &    '>>>> Particle splitting is ON >>>>'
      else if (buffer.eq.'OFF') then
        active = .false.
        write(*, '(a)')
     &    '>>>> Particle splitting is OFF >>>>'
        do
          read(*,'(a80)',iostat=error) buffer
          if (error.ne.0) then
            write(*,'(a,a,a)') 'inisplit:ERROR: ',
     &       'Unable to find End-Of-Section mark: ',eos
            stop
          endif
          if (index(buffer,eos).ne.0) return
        enddo
      else
        write(*,'(a)')
     &    'inisplit:ERROR: expecting to find ON or OFF'
        write(*,'(a)') 'found instead:'
        write(*,'(a)') buffer
        stop
      endif

      write(*,'(a)') 'Min weight to apply:'
      read(*,*) minwght
      write(*,'(1x,es12.5)') minwght
      if (minwght.lt.0.0) then
        write(*,'(a,i0)') 'inisplit:ERROR: Min weight must be >0.'
        stop
      endif

      write(*,'(a)') 'Splitting material:'
      read(*,*) matsplit
      write(*,'(1x,i0)') matsplit

      write(*,'(a)') 'Splitting mode:'
      read(*,*) mode
      select case (mode)
      case (1)
        write(*,'(2x,a)') 'simple'
      case (2)
        write(*,'(2x,a)') 'rotational'
      case (3)
        write(*,'(2x,a)') 'xy'
      case default
        write(*,'(a,i0)') 'inisplit:ERROR: invalid mode: ',
     &    mode
        stop
      end select

      if (ipol.eq.1.and.mode.ne.1) then
        write(*,'(a,i0)')
     &   'inisplit:ERROR: Photon polarisation, which'//
     &    'may break some symmetries, is active and,'
        write(*,'(a,i0)')
     &   '  therefore, only simple splitting can be used.'
        stop
      endif

      if (mode.ne.3) then
        write(*,'(a)') 'Splitting factor:'
        read(*,*) nsplit
        write(*,'(1x,i0)') nsplit
        if (nsplit.le.1) then
          write(*,'(a)') 'inisplit:ERROR: invalid factor'
          stop
        endif
      else
        read(*,'(a80)') buffer  ! Dummy line
        write(*,'(a)') '  factor for XY splitting is 4.'
      endif

      if (mode.ne.1) then
        write(*,'(a)') 'Euler angles for splitting body (deg):'
        read(*,*) omega,theta,phi
        write(*,'(3(1x,es12.5))') omega,theta,phi
        omega = omega*deg2rad
        theta = theta*deg2rad
        phi = phi*deg2rad

        write(*,'(a)') 'Shift for splitting body (cm):'
        read(*,*) dr
        write(*,'(3(1x,es12.5))') dr

        write(*,'(a)') 'Sign of z-direction:'
        read(*,'(a80)') buffer
        buffer = adjustl(buffer)  ! Remove leading blanks
        select case (buffer(1:1))
        case ('+')
          write(*,'(1x,a)') 'positive'
          signsplit = +1
        case ('-')
          write(*,'(1x,a)') 'negative'
          signsplit = -1
        case ('0')
          write(*,'(1x,a)') 'any'
          signsplit = 0
        case default
          write(*,'(a,a)') 'inisplit:ERROR: invalid sign: ',
     &      buffer(1:1)
          stop
        end select
      else
        read(*,'(a80)') buffer  ! Dummy line
        read(*,'(a80)') buffer  ! Dummy line
        read(*,'(a80)') buffer  ! Dummy line
      endif

      if (mode.eq.2) then
        write(*,'(a)') 'Azimuthal sector, phi0 and Dphi (deg):'
        read(*,*) phi0,dphi
        write(*,'(2(1x,es12.5))') phi0,dphi
        if (phi0.lt.0.0d0.or.phi0.ge.360.0d0.or.
     &      dphi.le.0.0d0.or.dphi.gt.360.0d0) then
          write(*,'(a)') 'inisplit:ERROR: invalid sector'
          stop
        endif
        if (dphi.eq.360.0d0) then
          write(*,'(a)') ' (2pi azimuthal splitting)'
        else
          write(*,'(a)') ' (restricted azimuthal splitting)'
          mode = 4  ! Rotational with restricted azimuth
        endif
        phi0 = phi0*deg2rad  ! Convert to rad
        dphi = dphi*deg2rad
        idphi = 1.0d0/dphi
      else
        read(*,'(a80)') buffer  ! Dummy line
      endif

      ! Init rotation matrices:
      call rotation(omega,theta,phi,rotS)  ! This yields rot=R(So->S)
      do i=1,3
        rotz(i) = rotS(i,3)  ! 3rd column of R(So->S) = 3rd row of R(S->So)
      enddo

      select case (mode)
      case (2,4)  ! Azimuthal, compose three rotations: R(So->S) . Rz(dphi) . R(S->So)
        nsplit = nint(nsplit*twopi/dphi)  ! Total splitting in 2pi
        call rotation(twopi/nsplit,0.0d0,0.0d0,rot)
        rot = matmul(rotS,matmul(rot,transpose(rotS)))
      case (3)  ! xy, compose rotations and reflection: R(So->S) . Reflect . R(S->So)
        ! Set fixed terms:
        rot = 0.0d0
        rot(3,3) = +1.0d0
        ! Reflection x -> -x:
        rot(2,2) = +1.0d0
        rot(1,1) = -1.0d0
        rorx = matmul(rotS,matmul(rot,transpose(rotS)))
        ! Reflection y -> -y:
        rot(2,2) = -1.0d0
        rot(1,1) = +1.0d0
        rory = matmul(rotS,matmul(rot,transpose(rotS)))
        ! Reflection x -> -x and y -> -y:
        rot(2,2) = -1.0d0
        rot(1,1) = -1.0d0
        rorxy = matmul(rotS,matmul(rot,transpose(rotS)))
      end select

      read(*,'(a80)') buffer
      if (index(buffer,eos).eq.0) then
        write(*,'(a)')
     &    'inisplit:ERROR: End-Of-Section mark not found'
        write(*,'(a,a)') '  expecting to find: ',eos
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif

      write(*,'(a)')
     &  '>>>> Particle splitting initialization finished >>>>'
      end


      subroutine splitting
!*******************************************************************
!*    Applies particle splitting.                                  *
!*******************************************************************
      use TRACK_mod
      use splitMod
      implicit none
      integer i,kmin,kmax,dk,dkout
      real*8 rot1(3,3),r(3),d(3),rtmp(3),dtmp(3),r0(3),d0(3),phi
      real*8, external :: rand
      real*8 pi,twopi,eps,phibin,iphibin
      parameter (pi=3.1415926535897932d0,twopi=2.0*pi,eps=1.0d-10)

      if (.not.active.or.mat.ne.matsplit.or.wght.lt.minwght) return  ! Nothing to do

      select case (mode)
      case (1)                                     ! Simple vanilla splitting
        wght = wght/nsplit
        do i=2,nsplit
          call stores(e,x,y,z,u,v,w,wght,kpar,ilb,ipol) ! Store particle in stack
        enddo
        call tally(-96,e*(nsplit-1))               ! Tally new particles

      case (2)                                     ! Rotational
        d(1) = u
        d(2) = v
        d(3) = w
        if (dot_product(rotz,d)*signSplit.lt.0.0) return ! Don't split if wrong direction
        r(1) = x-dr(1)                             ! Apply inverse shift
        r(2) = y-dr(2)
        r(3) = z-dr(3)
        wght = wght/nsplit
        do i=2,nsplit                              ! Each iteration rotates vectors by DAngle/nsplit
          r = matmul(rot,r)
          d = matmul(rot,d)
          call stores(e,r(1)+dr(1),r(2)+dr(2),r(3)+dr(3),
     &                d(1),d(2),d(3),wght,kpar,ilb,0) ! Store particle in stack
        enddo
        call tally(-96,e*(nsplit-1))               ! Tally new particles

      case (3)                                     ! xy
        d(1) = u
        d(2) = v
        d(3) = w
        if (dot_product(rotz,d)*signSplit.lt.0.0) return ! Don't split if wrong direction
        r(1) = x-dr(1)                             ! Apply inverse shift
        r(2) = y-dr(2)
        r(3) = z-dr(3)
        wght = wght*0.25                           ! Splitting factor is always 4 for XY symmetry

        ! Store x reflection:
        rtmp = matmul(rorx,r)
        dtmp = matmul(rorx,d)
        call stores(e,rtmp(1)+dr(1),rtmp(2)+dr(2),rtmp(3)+dr(3),
     &              dtmp(1),dtmp(2),dtmp(3),wght,kpar,ilb,0)

        ! Store y reflection:
        rtmp = matmul(rory,r)
        dtmp = matmul(rory,d)
        call stores(e,rtmp(1)+dr(1),rtmp(2)+dr(2),rtmp(3)+dr(3),
     &              dtmp(1),dtmp(2),dtmp(3),wght,kpar,ilb,0)

        ! Store xy reflection:
        rtmp = matmul(rorxy,r)
        dtmp = matmul(rorxy,d)
        call stores(e,rtmp(1)+dr(1),rtmp(2)+dr(2),rtmp(3)+dr(3),
     &              dtmp(1),dtmp(2),dtmp(3),wght,kpar,ilb,0)

        call tally(-96,e*3)                        ! Tally new particles

      case (4)                                     ! Rotational with restricted azimuth
        d(1) = u
        d(2) = v
        d(3) = w
        if (dot_product(rotz,d)*signSplit.lt.0.0) return ! Don't split if wrong direction
        r(1) = x-dr(1)                             ! Apply inverse shift
        r(2) = y-dr(2)
        r(3) = z-dr(3)

        rot1 = transpose(rotS)                     ! Temporary buffer for R(S->So)
        r0 = matmul(rot1,r)                        ! Transport position to ref system So
        d0 = matmul(rot1,d)
        if (r0(1).eq.0.0d0.and.r0(2).eq.0.0d0) then ! Split at the origin
          wght = wght/nsplit
          do i=2,nsplit
            call stores(e,x,y,z,u,v,w,wght,kpar,ilb,0)
          enddo
          call tally(-96,e*(nsplit-1))             ! Tally new particles
          return
        endif
        phi = atan2(r0(2),r0(1))
        if (phi.gt.phi0) phi = phi-twopi           ! Back-turn 2pi if not ordered

        phibin = twopi/nsplit
        iphibin = 1.0d0/phibin
        kmin = int((phi0-phi-eps)*iphibin)+1       ! First clone in (phi0,phi0+dphi)
        kmax = int((phi0+dphi-phi+eps)*iphibin)    ! Last clone in (phi0,phi0+dphi)
        call rotation(kmin*phibin,0.0d0,0.0d0,rot1)! Rotates into (phi0,phi0+dphi) in So
        r = matmul(rotS,matmul(rot1,r0))           ! Rotate in So and move back to S
        d = matmul(rotS,matmul(rot1,d0))
        wght = wght/nsplit                         ! Weight takes into account ALL clones
        call stores(e,r(1)+dr(1),r(2)+dr(2),r(3)+dr(3),
     &              d(1),d(2),d(3),wght,kpar,ilb,0)! Store 1st particle in (phi0,phi0+dphi)
        dk = kmax-kmin+1                           ! N# clones inside (phi0,phi0+dphi)
        do i=2,dk                                  ! Store all other clones inside (phi0,phi0+dphi)
          r = matmul(rot,r)                        ! Rotate R(So->S) . Ro(dphi) . R(S->So)
          d = matmul(rot,d)
          call stores(e,r(1)+dr(1),r(2)+dr(2),r(3)+dr(3),
     &                d(1),d(2),d(3),wght,kpar,ilb,0)
        enddo
        call tally(-96,e*dk)                       ! Tally new particles

        ! Now for the russian-rouletted clones, outside (phi0,phi0+dphi):
        dkout = nsplit-dk                          ! N# angles outside (phi0,phi0+dphi)
        if (dkout.gt.0) then                       ! At least 1 part outside
          call rotation((kmax+1+int(rand(0.0d0)*dkout))*phibin,
     &                   0.0d0,0.0d0,rot1)         ! A random angle out of (phi0,phi0+dphi)
          r = matmul(rotS,matmul(rot1,r0))         ! Rotate in So and move back to S
          d = matmul(rotS,matmul(rot1,d0))
          wght = wght*dkout
          call stores(e,r(1)+dr(1),r(2)+dr(2),r(3)+dr(3),
     &                d(1),d(2),d(3),wght,kpar,ilb,0) ! Only one survives
          call tally(-96,e)                        ! Tally new particle
        endif

        wght = -1.0                                ! Will kill current particle

      case default
        write(*,'(a)') 'splitting:ERROR: internal error.'
        stop
      end select
      end


      subroutine inirussia
!*******************************************************************
!*    Initializes the Russian roulette routines                    *
!*                                                                 *
!*    Input:                                                       *
!*      e0 -> source energy (eV)                                   *
!*******************************************************************
      use russiaMod
      implicit none
      character*(*) secid,eos
      parameter (secid=
     &'[SECTION RUSSIAN ROULETTE v.2009-06-15]')
      parameter (eos='[END OF VRRR SECTION]')
      character*80 buffer
      integer error

      write(*,*) ''
      write(*,'(a)')
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      call getline(buffer,0)
      if (index(buffer,secid).eq.0) then
        write(*,'(a)') 'inirussia:ERROR: incorrect section header;'
        write(*,'(a,a)') '  expecting to find: ',secid
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif
      write(*,'(a)') secid

      read(*,'(a80)') buffer
      buffer = adjustl(buffer)
      buffer = buffer(1:scan(buffer,' ')) ! Clip at 1st blank
      if (buffer.eq.'ON') then
        active = .true.
        write(*, '(a)')
     &    '>>>> Russian roulette is ON >>>>'
      else if (buffer.eq.'OFF') then
        active = .false.
        write(*, '(a)')
     &    '>>>> Russian roulette is OFF >>>>'
        do
          read(*,'(a80)',iostat=error) buffer
          if (error.ne.0) then
            write(*,'(a,a,a)') 'inirussia:ERROR: ',
     &       'Unable to find End-Of-Section mark: ',eos
            stop
          endif
          if (index(buffer,eos).ne.0) return
        enddo
      else
        write(*,'(a)')
     &    'inirussia:ERROR: expecting to find ON or OFF'
        write(*,'(a)') 'found instead:'
        write(*,'(a)') buffer
        stop
      endif

      write(*,'(a)') 'Max weight to apply:'
      read(*,*) maxwght
      write(*,'(1x,es12.5)') maxwght

      write(*,'(a)') 'Russian material:'
      read(*,*) matrussia
      write(*,'(1x,i0)') matrussia

      write(*,'(a)') 'Survival probability:'
      read(*,*) psurv
      write(*,'(1x,es12.5)') psurv
      if (psurv.lt.0.0d0.or.psurv.gt.1.0d0) then
        write(*,'(a)') 'inirussia:ERROR: invalid survival prob.'
        stop
      endif

      read(*,'(a80)') buffer
      if (index(buffer,eos).eq.0) then
        write(*,'(a)')
     &    'inirussia:ERROR: End-Of-Section mark not found'
        write(*,'(a,a)') '  expecting to find: ',eos
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif

      write(*,'(a)')
     &  '>>>> Russian roulette initialization finished >>>>'
      end


      subroutine russian
!*******************************************************************
!*    Applies Russian roulette.                                    *
!*******************************************************************
      use TRACK_mod
      use russiaMod
      implicit none
      real*8, external :: rand

      if (.not.active.or.mat.ne.matrussia.or.wght.gt.maxwght) return  ! Nothing to do

      if (rand(0.0d0).gt.psurv) then  ! Kill it
        wght = -1.0                   ! Will kill part at next interaction cycle
      else                            ! It survives
        wght = wght/psurv             ! Increase weight
      endif
      end


!>>>> End Of File >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
