!*******************************************************************
!*                    TALLY ENERGY DEPOSITION                      *
!*                                                                 *
!* Short description:                                              *
!*   Calculates the energy deposited in each material.             *
!* Rev:                                                            *
!*   2015-01-07                                                    *
!*******************************************************************


!*******************************************************************
!*******************************************************************
!*    Modules section                                              *
!*******************************************************************
!*******************************************************************

      module EDPmod
!*******************************************************************
!*    Vars for the EDP tally.                                      *
!*******************************************************************
      use PENELOPE_mod, only : maxmat
      implicit none
      save
      logical active
      integer matdet
      real*8 edptmp(maxmat),edep(maxmat),edep2(maxmat)
      real*8 unclimit
      end


!*******************************************************************
!*******************************************************************
!*    Routines section                                             *
!*******************************************************************
!*******************************************************************

      subroutine EDPtally(mode,eloss)
!*******************************************************************
!*    Input:                                                       *
!*      mode:   identifies the state of the calling procedure      *
!*      eloss:   energy deposition                                 *
!*******************************************************************
      use EDPmod
      use TRACK_mod, only : mat,wght
      use PENELOPE_mod, only : nmat
      implicit none
      integer mode
      real*8 eloss

      integer i

      if (.not.active) return

      if (mode.le.0.and.mat.gt.0) then        ! There's energy to be deposited
        edptmp(mat) = edptmp(mat)+eloss*wght  ! Update temporary E counter

      else if (mode.eq.6) then                ! End-of-history bookkeeping
        do i=1,nmat                           ! For all existing materials
          if (edptmp(i).eq.0.0) cycle         ! Skip void counters
          edep(i)   = edep(i) +edptmp(i)      ! Transfer temporary counter to mean and variance
          edep2(i)  = edep2(i)+edptmp(i)**2
          edptmp(i) = 0.0                     ! Clear counter to start a fresh history
        enddo

      endif
      end


      subroutine EDPreport(mode,n,cputim,uncdone)
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
      use EDPmod
      use PENELOPE_mod
      use RSEEDcommonMod
      implicit none
      integer mode,uncdone
      real*8 n,cputim

      integer i,out,finduf,errorcode
      real*8 q,q2,sigma,eff,uncert,invn

      uncdone = 1
      if (.not.active) return

      if (mode.ne.1) call EDPdump(1)  ! Sim dump file
      if (mode.eq.0) return           ! Dump only, no report

      ! Prepare output files:
      out = finduf()
      open(out,file='tallyEnergyDeposition.dat',iostat=errorcode)
      if (errorcode.ne.0) then
        write(*,*)
        write(*,'(a)')
     &    '*********************************************'
        write(*,'(a)')
     &    'EDPreport:ERROR: cannot open output data file'
        write(*,'(a)')
     &    '*********************************************'
        close(out)  ! Just in case
        return
      endif

      write(out,'(a)')
     &'#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      write(out,'(a)') '# [SECTION REPORT ENERGY DEPOSITION]'
      write(out,'(a)') '# Units are eV per history'
      write(out,'(a)') '#'
      write(out,'(a)') '# Material : Energy (eV/hist) : +-2sigma'

      invn = 1.0/n
      do i=1,nmat
        q  = edep(i)*invn
        q2 = edep2(i)*invn
        sigma = (q2-q**2)*invn
        sigma = sqrt(max(sigma,0.0))
        write(out,'(1x,i3,1x,es12.5,1x,es8.1)') i,q,2.0*sigma
      enddo

      ! Evaluate rel. uncertainty:
      q  = edep(matdet)*invn
      q2 = edep2(matdet)*invn
      sigma = (q2-q**2)*invn
      sigma = sqrt(max(sigma,0.0))
      uncert = 200.0
      if (q.gt.0.0) uncert = 200.0*sigma/q

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
     & '#   Uncertainty of energy deposition in detector material, '//
     & 'in % [uncert]:'
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

      subroutine EDPinitally(mem)
!*******************************************************************
!*    Initializes. To be called before TALLY.                      *
!*                                                                 *
!*    Output:                                                      *
!*      mem: Memory used by tally arrays                           *
!*******************************************************************
      use EDPmod
      use dataTypesMod
      use PENELOPE_mod, only : maxmat
      implicit none
      real*8 mem

      character*(*) secid,eos
      parameter (secid=
     &  '[SECTION TALLY ENERGY DEPOSITION v.2012-06-01]')
      parameter (eos='[END OF EDP SECTION]')
      character*80 buffer,upcase
      integer errorcode

      write(*,*) ' '
      write(*,'(a)')
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      call getline(buffer,0)
      if (index(buffer,secid).eq.0) then
        write(*,'(a)') 'EDPinitally:ERROR: incorrect section header;'
        write(*,'(a,a)') '  expecting to find: ',secid
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif
      write(*,'(a)') secid

      ! Report memory usage:
      mem = sizeOfReal8*3*real(maxmat)
      write(*,'(a)') 'Memory used by static arrays (MB):'
      write(*,'(1x,es12.5)') mem*1.0e-6

      read(*,'(a80)') buffer
      buffer = adjustl(upcase(buffer))
      buffer = buffer(1:scan(buffer,' ')) ! Clip at 1st blank
      if (buffer.eq.'ON') then
        active = .true.
      else if (buffer.eq.'OFF') then
        active = .false.
        write(*, '(a)') '>>>> Tally Energy Deposition is OFF >>>>'
        do
          read(*,'(a80)',iostat=errorcode) buffer
          if (errorcode.ne.0) then
            write(*,'(a,a,a)') 'EDPinitally:ERROR: ',
     &       'Unable to find End-Of-Section mark: ',eos
            stop
          endif
          if (index(buffer,eos).ne.0) return
        enddo
      else
        write(*,'(a)')
     &    'EDPinitally:ERROR: expecting to find ON or OFF'
        write(*,'(a)') 'found instead:'
        write(*,'(a)') buffer
        stop
      endif

      write(*,'(a)') 'Detection material set to:'
      read(*,*) matdet
      write(*,'(1x,i0)') matdet

      write(*,'(a)') 'Relative uncertainty (%) requested:'
      read(*,*) unclimit
      write(*,'(1x,es12.5)') unclimit

      read(*,'(a80)') buffer
      if (index(buffer,eos).eq.0) then
        write(*,'(a)')
     &    'EDPinitally:ERROR: End-Of-Section mark not found'
        write(*,'(a,a)') '  expecting to find: ',eos
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif

      ! Clear counters:
      edptmp = 0.0
      edep   = 0.0
      edep2  = 0.0

      call EDPdump(0)  ! Sim restart file

      write(*,'(a)') '>>>> EDP tally initialization finished >>>>'
      end


      subroutine EDPdump(mode)
!*******************************************************************
!*    Dumps into or reads data from a dump file.                   *
!*                                                                 *
!*    Input:                                                       *
!*      mode:   1 to write dump file, else to read from it.        *
!*******************************************************************
      use EDPmod, only : edep,edep2
      use PENELOPE_mod, only : nmat
      use dumpMod
      implicit none
      integer mode

      if (mode.eq.1) then
        if (dumpfile.lt.0) return  ! No dump file open
        write(dumpfile) edep(1:nmat),edep2(1:nmat)
      else
        if (restartfile.lt.0) return  ! No restart file open
        read(restartfile) edep(1:nmat),edep2(1:nmat)
      endif
      end


!>>>> End Of File >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
