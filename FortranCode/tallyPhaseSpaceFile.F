!*******************************************************************
!*                          TALLY                                  *
!*                    PHASE SPACE FILE                             *
!*                                                                 *
!* Short description:                                              *
!*   Writes to a file the dynamical state of all particles that    *
!*   reach a given detection material.                             *
!* Rev:                                                            *
!*   2015-04-30                                                    *
!*******************************************************************


!*******************************************************************
!*******************************************************************
!*    Modules section                                              *
!*******************************************************************
!*******************************************************************

      module PSFmod
!*******************************************************************
!*    Vars for the PSF tally.                                      *
!*******************************************************************
      implicit none
      save
      logical active
      integer detmat,psfunit,formatiaea
      integer*8 nhisti,nhlast,npar(4)
      end


      module PSFendinoMod
!*******************************************************************
!*    Vars for the endino particle---used for finishing a PSF, see *
!*    tally documentation.                                         *
!*******************************************************************
      implicit none
      save
      integer*4 gkpar,gilb(5)
      real*8 gen,gx,gy,gz,gu,gv,gw
      end


      module PSFiaeaMod
!*******************************************************************
!*    Vars for IAEA-formatted PSFs                                 *
!*******************************************************************
      implicit none
      save
      integer*4 sourceWrite
      end


!*******************************************************************
!*******************************************************************
!*    Routines section                                             *
!*******************************************************************
!*******************************************************************

      subroutine PSFtally(mode,arg)
!*******************************************************************
!*    Input:                                                       *
!*      mode -> Identifies the state of the calling procedure      *
!*      arg -> current history number (when mode=1)                *
!*******************************************************************
      use TRACK_mod
      use PSFmod
      use PSFendinoMod
      implicit none
      integer mode
      real*8 arg

      character*80 fmtstr
      parameter (fmtstr='(i0,8(1x,es14.7),6(1x,i0))')
      integer*4 dn

      if (.not.active) return

      if (mode.eq.4.or.mode.eq.-99) then  ! New material OR particle recovered from stack
        if (mat.ne.detmat) return         ! Not in detector
        dn = int(nhisti-nhlast)           ! Incremental history no.

        if (formatiaea.eq.0) then         ! Use native penEasy format (ASCII)
          write(psfunit,fmtstr) kpar,e,x,y,z,u,v,w,wght,dn,
     &                        ilb(1),ilb(2),ilb(3),ilb(4),ilb(5)

          gkpar = kpar                    ! Store endino, used at end-of-sim to keep sum(dn) consistent with N
          gen = e                         ! endino wght will be set =0, which identifies it
          gx = x
          gy = y
          gz = z
          gu = u
          gv = v
          gw = w
          gilb = ilb

#ifdef IAEAPSF
        else                              ! Use IAEA binary format
          call IAEAwrite(dn)              ! No endino for IAEA PSFs, N is already in a header

#endif
        endif

        nhlast = nhisti                   ! Store last contributing history
        npar(kpar) = npar(kpar)+1         ! Num. particles

      else if (mode.eq.1) then            ! New history begins
        nhisti = nint(arg,kind(nhisti))   ! Set history, use integer*8 NINT()

      endif
      end


      subroutine PSFreport(mode,n,cputim)
!*******************************************************************
!*    Input:                                                       *
!*      mode -> -1 if called at end-of-simulation, 0 for dump only,*
!*               1 for tally report only, 2 for both.              *
!*      n -> no. of histories simulated                            *
!*      cputim -> elapsed CPU time                                 *
!*    Comments:                                                    *
!*      -> 'cputim' should not include initialization procedures;  *
!*         enter 0 or neg. if not available.                       *
!*      -> When the PSF format is penEasy's native (ASCII), then   *
!*         a fictious particle with zero weight will be written to *
!*         the PSF at end-of-simulation to ensure that sum(DeltaN) *
!*         matches the actual history num. N.                      *
!*******************************************************************
      use TRACK_mod
      use RSEEDcommonMod
      use PSFmod
      use PSFendinoMod
      implicit none
      integer mode
      real*8 n,cputim

      integer out,finduf,error

      if (.not.active) return

      if (mode.eq.-1.and.formatiaea.eq.0.and.nhisti.ne.nhlast) then  ! Write a endino to the PSF
        wght = 0.0                ! Identifies the particle as a endino
        kpar = gkpar              ! Set endino state equal to last written part
        e    = gen
        x    = gx
        y    = gy
        z    = gz
        u    = gu
        v    = gv
        w    = gw
        ilb  = gilb
        mat  = detmat
        call PSFtally(-99,0.0d0)  ! Set switch to force PSFtally() to write; 2nd argument is dummy
      endif

      if (mode.ne.1) call PSFdump(1)  ! Write sim dump file
      if (mode.eq.0) return           ! Dump only, no report

      ! Prepare output files:
      out = finduf()
      open(out,file='tallyPhaseSpaceFile.dat',iostat=error)
      if (error.ne.0) then
        write(*,*)
        write(*,'(a)')
     &    '*********************************************'
        write(*,'(a)')
     &    'PSFreport:ERROR: cannot open output data file'
        write(*,'(a)')
     &    '*********************************************'
        close(out)  ! Just in case
        return
      endif

      write(out,'(a)')
     &'#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      write(out,'(a)') '# [SECTION REPORT PSF]'
#ifdef PROTONS
      write(out,'(a)')
     &  '# No. of electrons/photons/positrons/protons/total '//
     &  'written to PSF:'
#else
      write(out,'(a)')
     &  '# No. of electrons/photons/positrons/total '//
     &  'written to PSF:'
#endif
      write(out,'(1x,i0)') npar(1)
      write(out,'(1x,i0)') npar(2)
      write(out,'(1x,i0)') npar(3)
#ifdef PROTONS
      write(out,'(1x,i0)') npar(4)
#endif
      write(out,'(1x,i0)') sum(npar)

      ! Generic report:
      write(out,*) ''
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
      close(out)

#ifdef IAEAPSF
      if (formatiaea.eq.1) call IAEAwriteReport(mode)
#endif

      end


      subroutine PSFinitally(emax)
!*******************************************************************
!*    Initializes.                                                 *
!*                                                                 *
!*    Input:                                                       *
!*      emax -> max source energy (eV).                            *
!*******************************************************************
      use PENELOPE_mod
      use CSIMPHcommonMod, only : eabsh       ! Protons
      use PSFmod
      use PSFendinoMod
      use dumpMod
      implicit none
      real*8 emax

      character*(*) secid,eos
      parameter (secid=
     &'[SECTION TALLY PHASE SPACE FILE v.2009-06-15]')
      parameter (eos='[END OF PSF SECTION]')
      character*80 buffer,upcase
      integer finduf,error

      write(*,*) ' '
      write(*,'(a)')
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      call getline(buffer,0)
      if (index(buffer,secid).eq.0) then
        write(*,'(a)') 'PSFinitally:ERROR: incorrect section header;'
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
      else if (buffer.eq.'OFF') then
        active = .false.
        write(*, '(a)')
     &    '>>>> Tally Phase Space File is OFF >>>>'
        do
          read(*,'(a80)',iostat=error) buffer
          if (error.ne.0) then
            write(*,'(a,a,a)') 'PSFinitally:ERROR: ',
     &       'Unable to find End-Of-Section mark: ',eos
            stop
          endif
          if (index(buffer,eos).ne.0) return
        enddo
      else
        write(*,'(a)')
     &    'PSFinitally:ERROR: expecting to find ON or OFF'
        write(*,'(a)') 'found instead:'
        write(*,'(a)') buffer
        stop
      endif

      read(*,*) formatiaea
      if (formatiaea.eq.0) then
        write(*,'(a)') 'PSF format: standard penEasy in ASCII.'
      else if (formatiaea.eq.1) then
#ifdef IAEAPSF
        write(*,'(a)') 'PSF format: IAEA, in binary.'
#else
        write(*,'(a)')
     &    'PSFinitally:ERROR: IAEA PSF format requested '//
     &    'but not available.'
        write(*,'(a)') '  Read the instructions provided in the '//
     &    'documentation to activate this option.'
        stop
#endif
      else
        write(*,'(a)')
     &    'PSFinitally:ERROR: PSF format must be 0 or 1.'
        stop
      endif

      write(*,'(a)') 'Detection material set to:'
      read(*,*) detmat
      write(*,'(1x,i0)') detmat
      if (detmat.le.0.or.detmat.gt.maxmat) then
        write(*,'(a,i0)')
     &    'PSFinitally:ERROR: detection material out of range: 1,',
     &    maxmat
        stop
      endif
#ifdef PROTONS
      if (eabsh(detmat).lt.emax) then
        write(*,'(a)')
     &    'PSFinitally:ERROR: PSF detection material must be a '//
     &    'perfect absorbent; increase absorption energies above '//
     &    'the maximum source energy.'
        stop
      endif
#endif
      if (minval(eabs(:,detmat)).lt.emax) then
        write(*,'(a)')
     &    'PSFinitally:ERROR: PSF detection material must be a '//
     &    'perfect absorbent;'//
     &    'increase absorption energies above the maximum source '//
     &    'energy.'
        stop
      endif

      write(*,'(a)') 'PSF filename:'
      read(*,'(a80)') buffer
      buffer = adjustl(buffer)
      buffer = buffer(1:scan(buffer,' ')) ! Clip at 1st blank
      write(*,'(1x,a)') buffer

      if (formatiaea.eq.0) then
        psfunit = finduf()
        write(*,'(a)') 'Opening PSF as unit:'
        write(*,'(1x,i0)') psfunit

        if (restartfile.lt.0) then        ! This is not a restart
          open(psfunit,file=buffer,status='replace')
          write(psfunit,'(a)')
     &      '# [PHASE SPACE FILE FORMAT penEasy v.2008-05-15]'
          write(psfunit,'(a)')
     &      '# KPAR : E : X : Y : Z : U : V : W : WGHT : '//
     &      'DeltaN : ILB(1..5)'
        else                              ! This is a restart, append to the PSF
          open(psfunit,file=buffer,status='old',position='append',
     &         iostat=error)
          if (error.ne.0) then
            write(*,'(a)') 'PSFinitally:ERROR: Could not append PSF;'//
     &        ' make sure the file exists and it is readable.'
            stop
          endif
          write(psfunit,'(a)')
     &      '# Restarting simulation, PSF continues...'
        endif

#ifdef IAEAPSF
      else
        call IAEAiniwrite(buffer)
#endif

      endif

      ! Init vars:
      nhisti  = 0
      nhlast  = 0
      npar    = 0
      ! Init endino particle with arbitrary values, used only in case of an empty PSF:
      gkpar = 2
      gen = 1.0e3
      gx = 0.0
      gy = 0.0
      gz = 0.0
      gu = 0.0
      gv = 0.0
      gw = 1.0
      gilb = 0
      gilb(1) = 1

      call PSFdump(0)   ! Read sim dump file

      read(*,'(a80)') buffer
      if (index(buffer,eos).eq.0) then
        write(*,*) 'PSFinitally:ERROR: End-Of-Section mark not found'
        write(*,'(a,a)') '  expecting to find: ',eos
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif
      write(*,'(a)') '>>>> PSF tally initialization finished >>>>'
      end


      subroutine PSFdump(mode)
!*******************************************************************
!*    Dumps into or reads data from a dump file.                   *
!*                                                                 *
!*    Input:                                                       *
!*      mode -> 1 to write dump file, else to read from it.        *
!*******************************************************************
      use dumpMod
      use PSFmod, only : nhlast,npar
      implicit none
      integer mode

      if (mode.eq.1) then
        if (dumpfile.lt.0) return  ! No dump file open
        write(dumpfile) nhlast,npar
      else
        if (restartfile.lt.0) return  ! No restart file open
        read(restartfile) nhlast,npar
      endif
      end


!*******************************************************************
!*    IAEA specific routines follow.                               *
!*    Not adapted to protons: IAEAwrite will issue an error message*
!*    and abort if called.                                         *
!*******************************************************************
#ifdef IAEAPSF

      subroutine IAEAiniwrite(buffer)
!*******************************************************************
!*    Initializes.                                                 *
!*                                                                 *
!*    Input:                                                       *
!*      buffer -> PSF filename                                     *
!*******************************************************************
      use dumpMod
      use PSFiaeaMod
      implicit none
      character*80 buffer

      integer i
      integer*4 accessWrite,accessAppend,nXtraReal,nXtraInt
      integer*4 ind,xtraIntType,errorcode
      parameter (accessWrite=2,accessAppend=3)
      parameter (nXtraReal=0,nXtraInt=6)

      i = len_trim(buffer)+1
      buffer(i:i) = char(0)       ! C convention for end-of-string

      if (restartfile.lt.0) then  ! This is not a restart
        call iaea_new_source(sourceWrite,buffer,accessWrite,
     &                       errorcode)
        if (errorcode.lt.0) then
          write(*,'(a)') 'IAEAiniwrite:ERROR: Unable to open PSF.'
          stop
        endif
      else                        ! This is a restart, append to the PSF
        call iaea_new_source(sourceWrite,buffer,accessAppend,
     &                       errorcode)
        if (errorcode.lt.0) then
          write(*,'(a)') 'IAEAiniwrite:ERROR: '//
     &      'Could not open PSF to append; '//
     &      'make sure the file exists and it is accessible.'
          stop
        endif
      endif

      call iaea_set_extra_numbers(sourceWrite,nXtraReal,nXtraInt)
      errorcode = 0    ! Set to no-error

      ind = 0          ! First extra long-int number
      xtraIntType = 1  ! This IAEA code identifies DeltaN
      call iaea_set_type_extralong_variable(sourceWrite,ind,
     &  xtraIntType)
      errorcode = min(errorcode,xtraIntType)

      xtraIntType = 0  ! IAEA code identifies user-defined extra int, default if IAEAPSFold is active
      ind = 1
#ifndef IAEAPSFold
      ! Using IAEA routines that recognize ILB(1:4) codes
      xtraIntType = 7  ! This IAEA code identifies ilb(1)
#endif
      call iaea_set_type_extralong_variable(sourceWrite,ind,
     &  xtraIntType)
      errorcode = min(errorcode,xtraIntType)

      xtraIntType = 0
      ind = 2
#ifndef IAEAPSFold
      xtraIntType = 6  ! This IAEA code identifies ilb(2)
#endif
      call iaea_set_type_extralong_variable(sourceWrite,ind,
     &  xtraIntType)
      errorcode = min(errorcode,xtraIntType)

      xtraIntType = 0
      ind = 3
#ifndef IAEAPSFold
      xtraIntType = 5  ! This IAEA code identifies ilb(3)
#endif
      call iaea_set_type_extralong_variable(sourceWrite,ind,
     &  xtraIntType)
      errorcode = min(errorcode,xtraIntType)

      xtraIntType = 0
      ind = 4
#ifndef IAEAPSFold
      xtraIntType = 4  ! This IAEA code identifies ilb(4)
#endif
      call iaea_set_type_extralong_variable(sourceWrite,ind,
     &  xtraIntType)
      errorcode = min(errorcode,xtraIntType)

      ind = 5
      xtraIntType = 3  ! This IAEA code identifies ilb(5), already defined in old IAEA codes
      call iaea_set_type_extralong_variable(sourceWrite,ind,
     &  xtraIntType)
      errorcode = min(errorcode,xtraIntType)

      ! Note: xtraIntType=2 identifies LATCH variable in EGS codes

      select case (errorcode)
      case (-1)
        write(*,'(a)') 'IAEAiniwrite:ERROR: '//
     &    'Unable to set an extra variable.'
        stop
      case (-2)
        write(*,'(a)') 'IAEAiniwrite:internalERROR: '//
     &    'Extra variable index is out of range.'
        stop
      case (-3)
        write(*,'(a)') 'IAEAiniwrite:internalERROR: '//
     &    'Extra variable type is out of range.'
        stop
      case default
        if (errorcode.lt.0) then
          write(*,'(a)') 'IAEAiniwrite:internalERROR: '//
     &    'Undefined error while setting an extra variable.'
          stop
        endif
      end select
      end


      subroutine IAEAwrite(dn)
!*******************************************************************
!*    Write a particle to the PSF in IAEA format.                  *
!*                                                                 *
!*    Input:                                                       *
!*      dn -> Incremental history num.                             *
!*******************************************************************
      use TRACK_mod
      use PSFiaeaMod
      implicit none
      integer*4 dn

      integer*4 n_stat,kpariaea,xtraInt(6)
      real*4 xtraFloat(1)
      real*4 eiaea,xiaea,yiaea,ziaea,uiaea,viaea,wiaea,wghtiaea

      ! Convert particle index from PENELOPE to IAEA:
      select case (kpar)
      case (1)
        kpariaea = 2
      case (2)
        kpariaea = 1
      case (3)
        kpariaea = 3
      case (4)
        write(*,'(a)') 'IAEAwrite:ERROR: protons not implemented.'
        stop
      case default
        write(*,'(a)') 'IAEAwrite:ERROR: Invalid KPAR.'
        stop
      end select

      eiaea = real(e)*1.0e-6  ! Energy in MeV
      xiaea = real(x)         ! Vars type should match C++ definition for 'float'
      yiaea = real(y)
      ziaea = real(z)
      uiaea = real(u)
      viaea = real(v)
      wiaea = real(w)
      wghtiaea = real(wght)

      xtraInt(1) = dn
      xtraInt(2) = ilb(1)
      xtraInt(3) = ilb(2)
      xtraInt(4) = ilb(3)
      xtraInt(5) = ilb(4)
      xtraInt(6) = ilb(5)
      xtraFloat(1) = 0.0      ! Dummy value, not written to PSF

      n_stat = dn             ! Should be 0 if same history, >0 else

      call iaea_write_particle(sourceWrite,n_stat,kpariaea,eiaea,
     &  wghtiaea,xiaea,yiaea,ziaea,uiaea,viaea,wiaea,xtraFloat,
     &  xtraInt)

      if (n_stat.eq.-1) then
        write(*,'(a)') 'IAEAwrite:ERROR: Unable to write particle.'
        stop
      endif
      end


      subroutine IAEAwriteReport(mode)
!*******************************************************************
!*    Uptates IAEA header file.                                    *
!*                                                                 *
!*    Input:                                                       *
!*      mode -> -1 if called at end-of-simulation, 0 for dump only,*
!*               1 for tally report only, 2 for both.              *
!*******************************************************************
      use PSFmod
      use PSFiaeaMod
      implicit none
      integer mode

      integer*4 errorcode
      integer(selected_int_kind(R=18)) :: nTopiaea

      call flush(6)  ! Prevents messy output when using Fortran+C

      nTopiaea = nhisti
      call iaea_set_total_original_particles(sourceWrite,nTopiaea)
      if (nTopiaea.lt.0) then
        write(*,'(a)') 'IAEAwriteReport:ERROR: '//
     &    'unable to write no. of histories to header file.'
        stop
      endif

      call iaea_update_header(sourceWrite,errorcode)
      if (errorcode.lt.0) then
        write(*,'(a)') 'IAEAwriteReport:ERROR: '//
     &    'unable to update header file.'
        stop
      endif

      if (mode.eq.-1) then  ! Close the PSF
        call iaea_destroy_source(sourceWrite,errorcode)
        if (errorcode.lt.0) then
          write(*,'(a)') 'IAEAwriteReport:ERROR: '//
     &      'unable to close PSF.'
          ! Since already at end-of-sim, do not stop here
        endif
      endif
      end

#endif
!*******************************************************************
!*    End of IAEA block.                                           *
!*******************************************************************


!>>>> End Of File >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
