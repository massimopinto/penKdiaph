!*******************************************************************
!*                           PENAUX                                *
!*                                                                 *
!* Short description:                                              *
!*   Miscelaneous routines for penEasy.                            *
!* Rev:                                                            *
!*   2015-02-05                                                    *
!*******************************************************************


!*******************************************************************
!*******************************************************************
!*    Modules section                                              *
!*******************************************************************
!*******************************************************************

      module dataTypesMod
!*******************************************************************
!*    Defines the size in bytes of the data types employed by      *
!*    penEasy. (Superseeded in Fortran2008 by means of the         *
!*    ISO_FORTRAN_ENV module.)                                     *
!*******************************************************************
      implicit none
      save
      private
      character c
      integer i
      integer*4 i4
      integer*8 i8
      real r
      real*8 r8
      !! The following intrinsic STORAGE_SIZE works only on 2008-compliant compilers.
      !! For better compatibility, it is here commented out.
      !! integer, parameter :: bitsPerByte=8
      !! integer, parameter, public ::
      !! & sizeOfChar  = storage_size(c) /bitsPerByte,
      !! & sizeOfInt   = storage_size(i) /bitsPerByte,
      !! & sizeOfInt4  = storage_size(i4)/bitsPerByte,
      !! & sizeOfInt8  = storage_size(i8)/bitsPerByte,
      !! & sizeOfReal  = storage_size(r) /bitsPerByte,
      !! & sizeOfReal8 = storage_size(r8)/bitsPerByte,
      integer, parameter, public ::
     & sizeOfChar  = 1,
     & sizeOfInt   = 4,  ! CAUTION: may not apply to some (older) systems
     & sizeOfInt4  = 4,
     & sizeOfInt8  = 8,
     & sizeOfReal  = 4,
     & sizeOfReal8 = 8,
     & kindChar    = kind(c),
     & kindInt     = kind(i),
     & kindInt4    = kind(i4),
     & kindInt8    = kind(i8),
     & kindReal    = kind(r),
     & kindReal8   = kind(r8)
      end


      module ctrsimMod
!*******************************************************************
!*    Vars for simulation control.                                 *
!*******************************************************************
      implicit none
      save
      real*8 nhist,nhistmax,time0,atime
      real*8 lastComm,lastFresh,refresh,lastDump,refreshDump
      end


      module dsmaxMod
!*******************************************************************
!*    Vars for step length control.                                *
!*******************************************************************
      use PENELOPE_mod, only : maxmat
      implicit none
      save
      real*8 dsmax(maxmat)
      end


      module dumpMod
!*******************************************************************
!*    Vars for dump file.                                          *
!*******************************************************************
      implicit none
      save
      character*80 dumpfilen
      integer restartfile,dumpfile
      end


!*******************************************************************
!*******************************************************************
!*    PENELOPE's 'modularized' commons section.                    *
!*    Makes PENELOPE common blocks available via module blocks for *
!*    better maintainability.                                      *
!*******************************************************************
!*******************************************************************

      module RSEEDcommonMod
!*******************************************************************
      implicit none
      save
      integer*4 seed1,seed2
      common/rseed/seed1,seed2
      end


      module SECSTcommonMod
!*******************************************************************
      use PENELOPE_mod, only : nms
      implicit none
      save
      ! Particle stack:
      integer*4 ks,ibodys,ms,ilbs,ipols,nsec
      real*8 es,xs,ys,zs,us,vs,ws,wghts,sp1s,sp2s,sp3s,pages
      common/secst/es(nms),xs(nms),ys(nms),zs(nms),us(nms),vs(nms),
     &  ws(nms),wghts(nms),sp1s(nms),sp2s(nms),sp3s(nms),pages(nms),
     &  ks(nms),ibodys(nms),ms(nms),ilbs(5,nms),ipols(nms),nsec
      end


      module COMPOScommonMod
!*******************************************************************
      use PENELOPE_mod, only : maxmat
      implicit none
      save
      integer*4 iz,nelem
      real*8 stf,zt,at,rho,vmol
      common/compos/stf(maxmat,30),zt(maxmat),at(maxmat),rho(maxmat),
     &  vmol(maxmat),iz(maxmat,30),nelem(maxmat)
      end


      module PENGEOMcommonsMod
!*******************************************************************
      use pengeom_mod, only : nbody,ns,nb,nxg,mater
      implicit none
      save
      private  ! Prevent messing with PENGEOM vars

      integer*4 nsurf,kplane
      real*8 axx,axy,axz,ayy,ayz,azz,ax,ay,az,a0
      common/qsurf/axx(ns),axy(ns),axz(ns),ayy(ns),ayz(ns),azz(ns),
     1    ax(ns),ay(ns),az(ns),a0(ns),nsurf,kplane(ns)
      public :: ax,ay,az,a0,nsurf,kplane

      integer*4 nbodys,kmoth,kdght,ksurf,kflag,ksp,nwarn
      common/qtree/nbodys,kmoth(nb),kdght(nb,nxg),ksurf(nb,nxg),
     1  kflag(nb,nxg),ksp(ns),nwarn
      public :: kflag,ksurf
      end


      module CSIMPHcommonMod       ! Protons
!*******************************************************************
      use PENELOPE_mod, only : maxmat
      implicit none
      save
      real*8 eabsh,c1h,c2h,wcch,wcrh
      common/csimph/eabsh(maxmat),c1h(maxmat),c2h(maxmat),
     &              wcch(maxmat),wcrh(maxmat)
      end


!*******************************************************************
!*******************************************************************
!*    Routines section                                             *
!*******************************************************************
!*******************************************************************

      subroutine iniconfig(restartf)
!*******************************************************************
!*    Reads simulation configuration.                              *
!*                                                                 *
!*    Output:                                                      *
!*      restartf -> Restart file ID.                               *
!*      Other output is through modules.                           *
!*******************************************************************
      use dumpMod
      use ctrsimMod
      use PENELOPE_mod
      use RSEEDcommonMod
      implicit none
      integer restartf

      character*80 buffer
      character*(*) secid,eos
      parameter (secid='[SECTION CONFIG v.2013-03-18]')
      parameter (eos='[END OF CONFIG SECTION]')
      integer ufile,finduf,errorcode
      real*8 nmax
      parameter (nmax=1.0d15)  ! Max number of allowed histories, to prevent counter overflow

      write(*,*) ''
      write(*,'(a)')
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      call getline(buffer,0)
      if (index(buffer,secid).eq.0) then
        write(*,'(a)') 'iniconfig:ERROR: incorrect section header;'
        write(*,'(a,a)') '  expecting to find: ',secid
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif
      write(*,'(a)') secid

      write(*,'(a)') 'No. of histories requested: '
      read(*,*) nhistmax
      write(*,'(es12.5)') nhistmax
      if (nhistmax.gt.nmax) then
        write(*,'(a,es12.5)')
     &    'iniconfig:ERROR: No. of requested histories exceeds ',nmax
        stop
      endif

      ! Sim time:
      write(*,'(a)') 'Allotted time (s):'
      read(*,*) atime
      write(*,'(es12.5)') atime

      write(*,'(a)') 'Update interval:'
      read(*,*) refresh
      write(*,'(es12.5)') refresh
      if (refresh.lt.0.0) then
        write(*,'(a)') 'iniconfig:ERROR: interval must be positive.'
        stop
      endif

      ! Random seeds:
      write(*,'(a)') 'Random seeds:'
      read(*,*) seed1,seed2
      write(*,'(2(1x,i0))') seed1,seed2
      restartfile = -1       ! Indicates no restart requested
      if (seed1.eq.0.and.seed2.eq.0)  then
        continue             ! External seeds file requested, to be processed later
      else if (seed1.eq.-1.and.seed2.eq.-1) then
        restartfile = 1      ! Restart sim from a previous dump file
      else if (seed1.lt.1.or.seed2.lt.1) then
        write(*,'(a)') 'iniconfig:ERROR: Invalid RNG seeds.'
        stop
      endif

      read(*,'(a80)') buffer                 ! Read seeds filename
      if (seed1.eq.0.and.seed2.eq.0) then
        buffer = adjustl(buffer)             ! Remove leading blanks
        buffer = buffer(1:scan(buffer,' '))  ! Clip at 1st blank
        write(*,'(a)') 'Reading seeds from external file named:'
        write(*,'(1x,a)') buffer
        ufile = finduf()
        open(ufile,file=buffer,status='old',iostat=errorcode)
        if (errorcode.ne.0) then
          write(*,'(a)') 'iniconfig:ERROR: unable to open seeds file.'
          stop
        endif
        read(ufile,*) seed1,seed2
        close(ufile)
        write(*,'(a)') 'Random seeds:'
        write(*,'(2(1x,i0))') seed1,seed2
        if (seed1.lt.1.or.seed2.lt.1) then
          write(*,'(a)') 'iniconfig:ERROR: Invalid RNG seeds.'
          stop
        endif
      endif

      ! Restart file:
      read(*,'(a80)') buffer
      if (restartfile.eq.1) then
        write(*,'(a)')
     &    '*** Restarting simulation from a previous run***'
        buffer = adjustl(buffer)
        buffer = buffer(1:scan(buffer,' '))
        write(*,'(a)') 'Opening restart file:'
        write(*,'(1x,a)') buffer
        restartfile = finduf()
        open(restartfile,file=buffer,status='old',access='stream',
     &       iostat=errorcode)
        if (errorcode.ne.0) then
          write(*,'(a)')
     &      'iniconfig:ERROR: unable to open restart file.'
          stop
        endif
      endif
      restartf = restartfile  ! Output file unit to inform main program

      ! Dump file:
      read(*,'(a80)') buffer
      buffer = adjustl(buffer)
      dumpfilen = buffer(1:scan(buffer,' '))  ! Trim trailing spaces
      dumpfile = -1          ! Indicates no dump file requested
      if (dumpfilen.ne.'-') then
        write(*,'(a)') 'Simulation dump requested. Dump file name:'
        write(*,'(1x,a)') dumpfilen
        dumpfile = 1
      endif

      ! Dump interval:
      write(*,'(a)') 'Interval between dumps:'
      read(*,*) refreshDump
      write(*,'(es12.5)') refreshDump
      if (refreshDump.lt.0.0) then
        write(*,'(a)') 'iniconfig:ERROR: interval must be positive.'
        stop
      endif

      ! Check section integrity:
      read(*,'(a80)') buffer
      if (index(buffer,eos).eq.0) then
        write(*,'(a)')
     &    'iniconfig:ERROR: End-Of-Section mark not found'
        write(*,'(a,a)') '  expecting to find: ',eos
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif
      write(*,'(a)') '>>>> Config initialization finished >>>>'
      end


      subroutine inigeo(nmatgeo,mem)
!*******************************************************************
!*    Reads geometry section and initializes PENGEOM and penVox.   *
!*                                                                 *
!*    Output:                                                      *
!*      nmatgeo: largest material index in GEO & VOX files         *
!*      mem: memory used by voxels arrays                          *
!*******************************************************************
      use PENELOPE_mod
      implicit none
      integer*4 nmatgeo
      real*8 mem

      character*80 buffer
      character*(*) secid,eos
      parameter (secid=
     &'[SECTION PENGEOM+PENVOX v.2009-06-15]')
      parameter (eos='[END OF GEO SECTION]')
      integer ufile,finduf,errorcode
      integer*4 nmatvox,nmatquad,nbody
      real*8 parinp(1)

      write(*,*) ''
      write(*,'(a)')
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      call getline(buffer,0)
      if (index(buffer,secid).eq.0) then
        write(*,'(a)') 'inigeo:ERROR: incorrect section header;'
        write(*,'(a,a)') '  expecting to find: ',secid
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif
      write(*,'(a)') secid

      ! Read quadrics file, if there is one:
      write(*,*) ''
      write(*,'(a)')
     &  '>>>> Initializing quadric geometry >>>>'
      read(*,'(a80)') buffer              ! Quadrics filename
      buffer = adjustl(buffer)
      buffer = buffer(1:scan(buffer,' ')) ! Clip at 1st blank
      if (buffer.eq.'-') then             ! No quadrics
        call setIsQuad(.false.)           ! Flag indicates no quadrics
        buffer = 'default.tmp.geo'
        write(*,'(a)') 'No quadric geometry file defined.'
        write(*,'(a,a,a)')
     &    'An empty quadric geometry will be '//
     &    'automatically created in file ',
     &    trim(buffer),')'
        call autogeo(buffer)
      else                                ! Quadrics present
        call setIsQuad(.true.)            ! Flag indicates there are quadrics
        write(*,'(a)') 'Opening quadric geometry file:'
        write(*,'(1x,a)') buffer
      endif

      ! Init quadric geom (PENGEOM):
      ufile = finduf()
      open(ufile,file=buffer,status='old',iostat=errorcode)
      if (errorcode.ne.0) then
        write(*,'(a)') 'inigeo:ERROR: unable to open quadrics file'
        stop
      endif
      write(*,'(a)') 'Now calling PENGEOM...'
      write(*,*) ''
      call geomin(parinp,1,nmatquad,nbody,ufile,6)
      close(ufile)
      write(*,*) ''
      write(*,'(a)')
     &  '>>>> Quadric geometry initialization finished >>>>'

      ! Init penVox, if there is a voxelized geo:
      call inivox(nmatvox,mem)

      ! Set no. of materials in geometry files:
      nmatgeo = max(nmatvox,nmatquad)
      write(*,*) ''
      write(*,'(a)') 'No. of materials referenced in the geometry:'
      write(*,'(1x,i0)') nmatgeo
      if (nmatgeo.gt.maxmat) then
        write(*,'(a)') 'inigeo:ERROR: too many materials; '//
     &    'enlarge MAXMAT parameter and recompile.'
        stop
      else if (nmatgeo.eq.0) then
        write(*,'(a)') 'inigeo:ERROR: no geometry defined.'
        stop
      endif

      ! Check section integrity:
      read(*,'(a80)') buffer
      if (index(buffer,eos).eq.0) then
        write(*,'(a)') 'inigeo:ERROR: End-Of-Section mark not found'
        write(*,'(a,a)') '  expecting to find: ',eos
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif
      write(*,*) ''
      write(*,'(a)') '>>>> Geometry initialization finished >>>>'
      end


      subroutine inipen(emax,nmatgeo)
!*******************************************************************
!*    Reads simulation parameters and initializes PENELOPE.        *
!*                                                                 *
!*    Input:                                                       *
!*      emax -> max source energy (eV)                             *
!*      nmatgeo -> max number of materials in GEO & VOX files      *
!*******************************************************************
      use dsmaxMod
      use PENELOPE_mod
      use CSIMPHcommonMod  ! Protons
      implicit none
      integer*4 nmatgeo
      real*8 emax

      character*80 buffer
      character*300 lbuffer  ! Longer buffer for transport param table
      character*20 matfile(maxmat),filename
      character*(*) secid,eos
#ifdef PROTONS
      parameter (secid=
     &'[SECTION PENELOPE v.2014-04-01+protons]')
#else
      parameter (secid=
     &'[SECTION PENELOPE v.2009-10-01]')
#endif
      parameter (eos='[END OF PEN SECTION]')
      integer i,errorcode,imat
      integer*4 nmatcfg
      real*8 infty,nowcpu,cputime,zero
      parameter (infty=1.0d30,zero=1.0d-30)

      write(*,*) ''
      write(*,'(a)')
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      call getline(buffer,0)
      if (index(buffer,secid).eq.0) then
        write(*,'(a)') 'inipen:ERROR: incorrect section header;'
        write(*,'(a,a)') '  expecting to find: ',secid
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif
      write(*,'(a)') secid

      write(*,*) ''
      read(*,'(a300)') lbuffer  ! Read header
#ifdef PROTONS
      write(*,'(a)') 'MAT#  FILENAME             EABS(e-)  EABS(ph)'//
     &  '  EABS(e+)  EABS(p)   C1        C2        WCC       WCR'//
     &  '  C1p       C2p       WCCp      WCRp        DSMAX'
#else
      write(*,'(a)') 'MAT#  FILENAME             EABS(e-)  EABS(ph)'//
     &  '  EABS(e+)  C1        C2        WCC       WCR       DSMAX'
#endif

      nmatcfg = 0
      dsmax = +infty   ! Init value to ensure DSMAX is not zero
      do i=1,maxmat+1  ! Allow one extra line to read end-of-list
        read(*,'(a300)') lbuffer

        ! Read mat index:
        read(lbuffer,*) imat
        if (imat.eq.0) exit  ! End-of-list
        if (i.eq.maxmat+1) then
          write(*,'(a,i0)')
     &      'inipen:ERROR: Max number of materials exceeded: ',maxmat
          stop
        endif
        if (imat.ne.i) then
          write(*,'(a,i0,a,i0)')
     &      'inipen:ERROR: Invalid MAT index; found ',imat,
     &      ', should have been ',i
          stop
        endif
        lbuffer = adjustl(lbuffer) ! Remove leading blanks
        lbuffer = lbuffer(scan(lbuffer,' '):) ! Remove 1st datum

        ! Read mat file name:
        lbuffer = adjustl(lbuffer)
        read(lbuffer,'(a20)') filename
        if (scan(filename,' ').ne.0)                ! else there are no spaces in filename
     &    filename = filename(1:scan(filename,' ')) ! Clip at 1st blank
        matfile(i) = filename
        nmatcfg = i
        lbuffer = adjustl(lbuffer) ! Remove leading blanks
        lbuffer = lbuffer(scan(lbuffer,' '):) ! Remove 1st datum

        ! Read transport parameters:
        eabs(1,i) = -1.0
#ifdef PROTONS
        read(lbuffer,*,iostat=errorcode)
     &    eabs(1,i),eabs(2,i),eabs(3,i),eabsh(i),
     &    c1(i),c2(i),wcc(i),wcr(i),c1h(i),c2h(i),wcch(i),wcrh(i),
     &    dsmax(i)
#else
        read(lbuffer,*,iostat=errorcode) eabs(1,i),eabs(2,i),eabs(3,i),
     &    c1(i),c2(i),wcc(i),wcr(i),dsmax(i)
#endif
        if (errorcode.ne.0) then
          if (eabs(1,i).ne.-1.0) then
            write(*,'(a,i0)')
     &       'inipen:ERROR: incomplete list of parameters for MAT ',i
            stop
          endif

          ! Set parameters automatically:
          eabs(1,i) = min(max(emax*1.0d-2,50.0d0),1.0d6)
          eabs(2,i) = min(max(emax*1.0d-3,50.0d0),1.0d6)
          eabs(3,i) = eabs(1,i)
          c1(i)     = 0.1d0
          c2(i)     = 0.1d0
          wcc(i)    = min(emax*1.0d-2,eabs(1,i))
          wcr(i)    = min(emax*1.0d-3,eabs(2,i))
          dsmax(i)  = +infty
          ! Protons (irrelevant if not active):
          eabsh(i)  = min(max(emax*1.0d-2,50.0d0),1.0d7)
          c1h(i)    = c1(i)
          c2h(i)    = c2(i)
          wcch(i)   = wcc(i)
          wcrh(i)   = wcr(i)
        endif

        ! Write data:
#ifdef PROTONS
        write(*,'(i2,4x,a20,13(1x,es9.3))')
     &    i,filename,eabs(1,i),eabs(2,i),eabs(3,i),eabsh(i),
     &    c1(i),c2(i),wcc(i),wcr(i),c1h(i),c2h(i),wcch(i),wcrh(i),
     &    dsmax(i)
#else
        write(*,'(i2,4x,a20,8(1x,es9.3))')
     &    i,filename,eabs(1,i),eabs(2,i),eabs(3,i),
     &    c1(i),c2(i),wcc(i),wcr(i),dsmax(i)
#endif
        if (dsmax(i).lt.zero) then
          write(*,'(a)')
     &      'inipen:ERROR: DSMAX must be larger than zero even if '//
     &      'electrons are not transported.'
          stop
        endif
      enddo

      if (nmatgeo.gt.nmatcfg) then
        write(*,'(a)')
     &    'inipen:ERROR: There are more materials declared in the '//
     &    'geometry file than defined in the config file.'
        stop
      endif

      ! Init PENELOPE kernel:
      write(*,*) ''
      write(*,'(a)') 'Now calling PEINIT...'
      write(*,*) ''
      nowcpu = cputime()                            ! Start stopwatch
      call peinit(emax*1.001d0,nmatcfg,6,1,matfile) ! Increase EMAX to allow absorbents; PEINIT limits EABS=EMAX*(1-1e-4)
#ifdef PROTONS
      call phinit(emax*1.001d0,nmatcfg,6,1,matfile)
#endif

      ! Report time usage:
      write(*,*) ''
      write(*,'(a)') 'CPU time used to init PENELOPE (s):'
      write(*,'(1x,f9.2)') cputime()-nowcpu

      ! Check section integrity:
      read(*,'(a80)') buffer
      if (index(buffer,eos).eq.0) then
        write(*,'(a)') 'inipen:ERROR: End-Of-Section mark not found'
        write(*,'(a,a)') '  expecting to find: ',eos
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif
      write(*,*) ''
      write(*,'(a)') '>>>> PENELOPE initialization finished >>>>'
      end


      logical function endsim()
!*******************************************************************
!*    Determines whether the simulation is done or not and writes  *
!*    progress reports.                                            *
!*                                                                 *
!*    Output:                                                      *
!*      .true. if done                                             *
!*******************************************************************
      use ctrsimMod
      implicit none

      logical notfinal
      parameter (notfinal=.false.)
      integer unc,mode
      real*8 realtime,cputime,nowreal,nowcpu,refreshComm
      parameter (refreshComm=30.0)  ! Command check real time (s) interval

      nowreal = realtime()  ! Ensures one call per day at least
      nowcpu  = cputime()
      endsim  = .false.

      ! Check no. of histories:
      if (nhist-nhistmax.gt.-0.5d0) endsim = .true.  ! 0.5 used for precision

      ! Check time:
      if (atime.gt.0.0) then
        if (nowreal.gt.atime) endsim = .true.    ! Real time
      else
        if (nowcpu.gt.-atime) endsim = .true.    ! User (CPU) time
      endif

      if (endsim) return                         ! Simulation is done

      ! Command check:
      if (nowreal-lastComm.gt.refreshComm) then
        call comand
        lastComm = nowreal
      endif

      mode = -1  ! Initial value, no report

      ! Check if an update is needed:
      if (nowreal-lastFresh.gt.refresh) then
        lastFresh = nowreal
        mode = 1
      endif

      ! Check if a dump is needed:
      if (nowreal-lastDump.gt.refreshDump) then
        lastDump = nowreal
        mode = mode+1                        ! 2 is previously 1, 0 if previously -1
      endif

      if (mode.ne.-1) then
        call tallyreport(mode,nhist,nowcpu,unc)  ! Writes progress report
        if (mode.ne.0) endsim = unc.gt.1     ! Required accuracy attained?
      endif
      end


      subroutine simstop
!*******************************************************************
!*    Forces the simulation to stop after completion of current    *
!*    history by re-setting the no. of histories to simulate to    *
!*    zero.                                                        *
!*******************************************************************
      use ctrsimMod, only : nhistmax
      implicit none
      nhistmax = 0.0d0
      end


      logical function absorb()
!*******************************************************************
!*    Checks whether or not a particle is absorbed given its       *
!*    current energy.                                              *
!*                                                                 *
!*    Input:                                                       *
!*      TRACK_mod                                                  *
!*    Output:                                                      *
!*      .true. if absorbed, .false. else                           *
!*    Comments:                                                    *
!*      -> Tallies deposited energy when necessary; to do so it    *
!*         calls TALLY.                                            *
!*      -> Absorbs particles above 1 GeV; these can be created by  *
!*         SOURCE or be the result of the annihilation of a e+.    *
!*      -> This function is well suited to implement range         *
!*         rejection by making the returned value dependent        *
!*         on the particle's position, material, energy, etc.      *
!*         Notice that in order to preserve the calculation of     *
!*         the fluence unbiased an electron or positron should     *
!*         never be rejected *inside* the detector.                *
!*******************************************************************
      use PENELOPE_mod, only : eabs
      use TRACK_mod
      use CSIMPHcommonMod  ! Protons
      implicit none
      real*8 mc2,twomc2,penemax
      parameter (mc2=5.10998918d5,twomc2=2.0d0*mc2,penemax=1.0d9)

      if (wght.lt.0.0) then     ! Indicates particle must be killed
        call tally(5,0.0d0)
        absorb = .true.
        return
      endif

      if (mat.eq.0) then        ! Particle is gone
        absorb = .true.
        return
      endif

      select case (kpar)
      case (1,2,3)
        if (e.gt.eabs(kpar,mat).and.e.lt.penemax) then
          absorb = .false.
          return
        endif
#ifdef PROTONS
      case (4)
        if (e.gt.eabsh(mat).and.e.lt.penemax) then
          absorb = .false.
          return
        endif
#endif
      case default
        write(*,'(a,i0)')
     &    'absorb:WARNING: unknown particle type: ',kpar
        write(*,'(a)')
     &    '  The particle has been discarded (no E deposition).'
        absorb = .true.
        return
      end select

      absorb = .true.
      call tally(-98,e)         ! Tally remaining kinetic energy
      if (kpar.eq.3.and.e.gt.0.0) then ! e+ absorbed, create 2 gammas; but not if E=0, since KNOCK did it
        call panar(eabs(2,mat)) ! Puts two annihilation photons in the stack if above Eabs
        call tally(-97,twomc2)  ! Tally the extra 2mc^2
      endif

      if (e.ge.penemax) then
        write(*,'(a,es12.5,a)')
     &  'absorb:WARNING: particle with energy out of range: ',e,' eV'
        write(*,'(a)')
     &  '  The particle has been locally absorbed.'
      endif
      end


      subroutine comand
!*******************************************************************
!*    Reads and executes commands from an external file, allowing  *
!*    in-flight steering of the simulation.                        *
!*                                                                 *
!*    Input:                                                       *
!*      -> Command is read from external file command.in           *
!*    Output:                                                      *
!*      -> Command is executed and file command.in reset to default*
!*******************************************************************
      use ctrsimMod
      implicit none
      character*80 buffer
      logical ffound
      integer com,in,finduf,errorcode
      real*8 param,nmax
      parameter (nmax=1.0d15)  ! Max number of allowed histories, to prevent counter overflow

      inquire(file='command.in',exist=ffound)
      if (.not.ffound) return  ! File not found
      in = finduf()
      open(in,file='command.in',iostat=errorcode)
      if (errorcode.ne.0) then
        close(in)  ! Just in case
        return
      endif

      ! Read command file:
      iter: do
        read(in,*,iostat=errorcode) com
        if (errorcode.ne.0.or.com.eq.0) exit iter
        write(*,*) ''
        write(*,'(a,f18.0)')
     &   'comand:WARNING: command received when history number was: ',
     &    nhist
        select case (com)
        case (1)
          read(in,*,iostat=errorcode) param
          if (errorcode.ne.0) exit iter
          nhistmax = param
          if (nhistmax.gt.nmax) then
            write(*,*) ''
            write(*,'(a,f18.0)')
     &      'comand:WARNING: requested max No. of histories exceeds ',
     &        nmax
            nhistmax = nmax
          endif
          write(*,'(a,f18.0)') '  max No. of histories reset to: ',
     &      nhistmax
        case (2)
          read(in,*,iostat=errorcode) param
          if (errorcode.ne.0) exit iter
          atime = param
          write(*,'(a,es12.5)') '  max time reset to: ',atime
        case (3)
          read(in,*,iostat=errorcode) param
          if (errorcode.ne.0) exit iter
          refresh = param
          write(*,'(a,es12.5)') '  update interval reset to: ',refresh
        case (4)
          read(in,*,iostat=errorcode) param
          if (errorcode.ne.0) exit iter
          refreshDump = param
          write(*,'(a,es12.5)') '  dump interval reset to: ',
     &      refreshDump
        case (9)
          read(in,'(a80)',iostat=errorcode) buffer
          if (errorcode.ne.0) exit iter
          write(*,'(a)') '  message from command.in follows:'
          write(*,'(a80)') buffer
        case default
          write(*,'(a,i0,a)')
     &     '  unknown command: ',com,'; remaining commands ignored.'
          exit iter
        end select
      enddo iter
      close(in,iostat=errorcode)

      ! Reset command file:
      open(in,file='command.in',iostat=errorcode)
      if (errorcode.ne.0) then
        close(in)  ! Just in case
        return
      endif
      write(in,'(a)',iostat=errorcode)
     &  ' 0   <- write here the command code'
      write(in,'(a)',iostat=errorcode)
     &  '     <- write here the command parameter, if any'
      write(in,'(a)',iostat=errorcode)
     &  '     <- repeat previous lines as many times as needed'
      write(in,'(a)',iostat=errorcode)
     &  ''
      write(in,'(a)',iostat=errorcode)
     &  '>>>> END OF INPUT >>>>'
      write(in,'(a)',iostat=errorcode)
     &  ''
      write(in,'(a)',iostat=errorcode)
     &  'Code  Parameter  What it does'
      write(in,'(a)',iostat=errorcode)
     &  '-----------------------------'
      write(in,'(a)',iostat=errorcode)
     &  '0     None       Nothing, keep going'
      write(in,'(a)',iostat=errorcode)
     &  '1     N          Reset No. of histories to N'
      write(in,'(a)',iostat=errorcode)
     &  '2     t          Reset simulation time to t seconds'
      write(in,'(a)',iostat=errorcode)
     &  '3     t          Reset update interval to t seconds'
      write(in,'(a)',iostat=errorcode)
     &  '4     t          Reset dump interval to t seconds'
      write(in,'(a)',iostat=errorcode)
     &  '9     <string>   Write <string> (<80 chars) to output'
      write(in,'(a)',iostat=errorcode)
     &  ''
      write(in,'(a)',iostat=errorcode)
     &  '>>>> END OF FILE >>>>'
      write(in,'(a)',iostat=errorcode)
     &  ''
      close(in,iostat=errorcode)
      end


      integer function finduf()
!*******************************************************************
!*    Finds a valid (i.e. unused) file unit.                       *
!*******************************************************************
      implicit none
      logical used
      integer maxuf
      parameter (maxuf=17)

      finduf = 6
      do
        finduf = finduf+1
        if (finduf.gt.maxuf) then
          write(*,*) ''
          write(*,'(a)') 'finduf:ERROR: Unable to find a valid unit'
          stop
        endif
        inquire(finduf,opened=used)
        if (.not.used) return
      enddo
      end


      subroutine getline(buffer,ufile)
!*******************************************************************
!*    Reads a new line from stdin or file. The line is returned    *
!*    only if it is not blank nor a comment line.                  *
!*                                                                 *
!*    Input:                                                       *
!*      ufile: unit file, must be already open; if =0, reads stdin *
!*    Output:                                                      *
!*      buffer: line read                                          *
!*    Comments:                                                    *
!*      - Comment lines start with a '#'                           *
!*******************************************************************
      implicit none
      character*80 buffer
      integer ufile

      integer errorcode

      do                                           ! Read while empty or comment lines
        if (ufile.eq.0) then                       ! Read from std input
          read(*,'(a80)',iostat=errorcode) buffer
        else                                       ! Read from unit file
          read(ufile,'(a80)',iostat=errorcode) buffer
        endif
        if (errorcode.ne.0) then
          write(*,'(a)') 'getline:ERROR: unable to read line.'
          stop
        endif
        if (buffer(1:1).eq.'#') cycle              ! A comment line
        if (verify(buffer,' ').ne.0) return        ! Non-blank line, done
      enddo                                        ! Blank line, read another
      end


      integer function seeki(x,xc,n)
!*******************************************************************
!*    Finds the interval (x(i),x(i+1)] containing the value xc.    *
!*                                                                 *
!*    Input:                                                       *
!*      x(1..n) -> data array                                      *
!*      xc -> point to be located                                  *
!*      n -> no. of data points                                    *
!*    Output:                                                      *
!*      index i of the semiopen interval where xc lies             *
!*    Comments:                                                    *
!*      -> If xc=x(1) then i=1 is returned.                        *
!*      -> If xc is outside the closed interval [x(1),x(n)]  the   *
!*         execution is aborted.                                   *
!*******************************************************************
      implicit none
      integer n
      real*8 xc,x(n)

      integer itop,imid

      if(xc.gt.x(n)) then
        write(*,'(a)') 'seeki:ERROR: value outside range, xc>x(n):'
        write(*,'(2(1x,es12.5))') xc,x(n)
        stop
      endif
      if(xc.lt.x(1)) then
        write(*,'(a)') 'seeki:ERROR: value outside range, xc<x(1):'
        write(*,'(2(1x,es12.5))') xc,x(1)
        stop
      endif

      seeki = 1
      itop = n
      do
        imid = (seeki+itop)/2
        if(xc.gt.x(imid)) then
          seeki = imid
        else
          itop = imid
        endif
        if(itop-seeki.le.1) exit
      enddo
      end


      subroutine rotation(omega,theta,phi,rot)
!*******************************************************************
!*    Computes the rotation matrix determined by the three Euler   *
!*    angles as defined in the PENELOPE manual.                    *
!*                                                                 *
!*    Input:                                                       *
!*      omega -> rotation angle around the z axis (rad)            *
!*      theta -> rotation angle around the y axis (rad)            *
!*      phi   -> rotation angle around the z axis (rad)            *
!*    Output:                                                      *
!*      rot -> 3x3 matrix of the 3D rotation                       *
!*******************************************************************
      implicit none
      real*8 omega,theta,phi,rot(3,3)

      real*8 somega,comega,stheta,ctheta,sphi,cphi

      somega = sin(omega)
      comega = cos(omega)
      stheta = sin(theta)
      ctheta = cos(theta)
      sphi   = sin(phi)
      cphi   = cos(phi)
      rot(1,1) = cphi*ctheta*comega-sphi*somega
      rot(1,2) = -cphi*ctheta*somega-sphi*comega
      rot(1,3) = cphi*stheta
      rot(2,1) = sphi*ctheta*comega+cphi*somega
      rot(2,2) = -sphi*ctheta*somega+cphi*comega
      rot(2,3) = sphi*stheta
      rot(3,1) = -stheta*comega
      rot(3,2) = stheta*somega
      rot(3,3) = ctheta
      end


      function crossProduct(u,w)
!*******************************************************************
!*    Computes the cross product of the two input vectors.         *
!*******************************************************************
      implicit none
      real*8 crossProduct(3),u(3),w(3)

      crossProduct(1) = u(2)*w(3) - u(3)*w(2)
      crossProduct(2) = u(3)*w(1) - u(1)*w(3)
      crossProduct(3) = u(1)*w(2) - u(2)*w(1)
      end


      function upcase(string) result(upper)
!*******************************************************************
!*    Converts a string to uppercase.                              *
!*    Adapted from http://www.star.le.ac.uk/~cgp/fortran.html      *
!*    No explicit interface needed for this version, 80 characters *
!*    assumed.                                                     *
!*******************************************************************
      implicit none
      character*80 string,upper
      integer j

      do j=1,len(string)
        if(string(j:j) >= "a" .and. string(j:j) <= "z") then
             upper(j:j) = achar(iachar(string(j:j))-32)
        else
             upper(j:j) = string(j:j)
        end if
      enddo
      end


      real*8 function gausspdf()
!*******************************************************************
!*    Returns a random value with a Gaussian PDF.                  *
!*    Uses the polar method to avoid expensive trigonometric calls *
!*    implied by the alternative Box-Muller method.                *
!*******************************************************************
      implicit none
      logical, save :: isStored=.false.  ! Flags if the 2nd r.v. of the polar method is ready to be delivered
      real*8, save :: stored             ! Contains the 2nd r.v.
      real*8 x,y,u,s
      real*8, parameter :: zero=1.0d-100
      real*8, external :: rand

      if (isStored) then   ! The 2nd r.v. was ready to be delivered
       gausspdf = stored
       isStored = .false.  ! Set flag for next call
       return
      endif

 10   x = 1.0-2.0*rand(0.0d0)
      y = 1.0-2.0*rand(1.0d0)
      u = x*x+y*y
      if (u.ge.1.0d0.or.u.lt.zero) goto 10  ! Reject point and repeat
      s = sqrt(-2.0*log(u)/u)
      gausspdf = x*s

      stored = y*s        ! Value stored and ready to be delivered in later calls
      isStored = .true.   ! Set flag
      end


!>>>> End Of File >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
