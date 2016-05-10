!*******************************************************************
!*                    TALLY ENERGY DEPOSITION                      *
!*                                                                 *
!* Short description:                                              *
!*   Calculates the energy deposited in each material.             *
!* Rev:                                                            *
!     *   2015-01-07                                               *
!     * Revision by Massimo to include Diaphragm effects tallies   *
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

      module DIAPHmod
!*******************************************************************
!     Vars for the calculation of diaphragm correction factors  ****
!     Massimo, from the 2011 version of penEasy, with updates
!     fro the penelope version of 2014                          ****
!*******************************************************************
      use PENELOPE_mod, only : maxmat
      implicit none
      save
      integer*4, parameter ::  maxk=6  ! number of diaphragm correction factors
!      real*8 edptmp,edep,edep2,unclimit,edepdtr,edepdtr2,
!     &     edptmpdiaph
      real*8  edepdtr(maxk),edepdtr2(maxk),edptmpdiaph(maxK)
!     Massimo: energy deposited from transmitted photons
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
      use DIAPHmod
      use TRACK_mod, only : mat,wght,ilb ! ilb needed for kdiaph corr's
      use PENELOPE_mod, only : nmat
      implicit none
      integer mode
      real*8 eloss

      integer i

      if (.not.active) return

      if (mode.le.0.and.mat.gt.0) then ! There's energy to be deposited

!     Writing temporary counters for selected type of events
       if(mat.eq.matdet) then ! if you are scoring in the detection material
                ! ask what is the origin of this event
            if ((ilb(5).eq.0).and.(ilb(1).eq.1)) then
            edptmpdiaph(1) = edptmpdiaph(1) +eloss*wght ! scores "Eprim" from electron's energy
            ! a photo- or Compton-electron was created and the initial energy is stored.
            ! the secondary photon is NOT contemplated in the energy count
            ! in agreement with what the free air chamber is to measure.
            ! endif
            else if (ilb(5).eq.1) then
            edptmpdiaph(1) = edptmpdiaph(1) +eloss*wght
            ! In this case you have that the above-mentioned electron is being transported and got a knockx()
            else if (ilb(5).eq.2) then  ! photon transmitted through diaphragm
            ! careful though because a primary photon might have generated an electron near the scoring volume
            ! which then enters the scoring volume later. This is not contemplated here, and it SHOULD.
            ! You must identify those photo- or Compton-electron that descend directly from a primary photon
                edptmpdiaph(2) = edptmpdiaph(2) +eloss*wght  ! scores transmitted photon's "Edtr"
            else if (ilb(5).eq.3) then
                edptmpdiaph(3) = edptmpdiaph(3) +eloss*wght  ! scores scattered photon's "Edsc"
            ! whereby scattering occurred in the diaphragm and not in air, later on.
            else if (ilb(5).eq.4) then
                edptmpdiaph(4) = edptmpdiaph(4) +eloss*wght  ! scores fluorescence photon's "Edfl"
            else if (ilb(5).eq.5) then
                edptmpdiaph(5) = edptmpdiaph(5) +eloss*wght  ! scores Bremms photon's "Edbr"
                ! This Bremsstrahlung photon was started in the diaphragm
            else if (ilb(5).eq.6) then
                edptmpdiaph(6) = edptmpdiaph(6) +eloss*wght  ! scores electron's Edel
            endif !  and so on...here you can add other ILB(5) conditions
         endif
        edptmp(mat) = edptmp(mat)+eloss*wght  ! Update temporary E counter


      else if (mode.eq.6) then                ! End-of-history bookkeeping
        do i=1,nmat                           ! For all existing materials
          if (edptmp(i).eq.0.0) cycle         ! Skip void counters
          edep(i)   = edep(i) +edptmp(i)      ! Transfer temporary counter to mean and variance
          edep2(i)  = edep2(i)+edptmp(i)**2
          edptmp(i) = 0.0                     ! Clear counter to start a fresh history
       enddo

        do i=1,maxK             ! energy counters specific for each physical process
            if (edptmpdiaph(i).gt.0.0) then  ! Transfer temporary specific counter to mean and variance
                edepdtr(i) = edepdtr(i) + edptmpdiaph(i)
                edepdtr2(i) = edepdtr2(i) + edptmpdiaph(i)**2
                edptmpdiaph(i) = 0.0
                ! Clears counter to start a fresh history
            endif
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
      use DIAPHmod
      use PENELOPE_mod
      use RSEEDcommonMod
      implicit none
      integer mode,uncdone
      real*8 n,cputim

      integer i,out,finduf,errorcode
      real*8 q,q2,sigma,eff,uncert,invn
      real*8 dtr,dtr2,sigmadtr,edepssum,sigmadetmat,sigmaeprim,
     & Kdfactors(6),numer,denom

      ! names of specific energy counters and K_correction factors
      ! Correction factors, see Burns and Kessler 2009
      character*5 Edepnames(maxK),Kdfacnames(maxK)
      Edepnames(1) = 'Eprim'
      Edepnames(2) = 'Edtr'
      Edepnames(3) = 'Edsc'
      Edepnames(4) = 'Edfl'
      Edepnames(5) = 'Edbr'
      Edepnames(6) = 'Edel'

      Kdfacnames(1) = 'Kdtr' ! photons transmitted through diaphragm
      Kdfacnames(2) = 'Kdsc' ! coherent or incoherent scatter from the diaphragm
      Kdfacnames(3) = 'Kdfl' ! fluorescence photons in the diaphragm
      Kdfacnames(4) = 'Kdbr' ! Bremsstrahlung photons in the diaphragm
      Kdfacnames(5) = 'Kdel' ! Electrons exiting from the diaphragm
      Kdfacnames(6) = 'Kdiap'   ! this is the product of all Kfactors = Eprim/Etot

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
       ! Massimo writes output with selective energy deposition counters
       write(out,'(a)')
       write(out,'(a)')
     & '#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
       write(out,'(a)') '# Report on selective energy deposition modes'
       write(out,'(a,1x,i3)') '# Limited to detection material', matdet
       write(out,'(a)') '# All energy deposition values in eV'
       write(out,'(a)')
     & '#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
       write(out,'(a)')
     &   '####   counter   E(eV) :   +-2sig'
        edepssum = 0.0
      do i=1,maxk ! cycle through the energy deposition counters
        dtr  = edepdtr(i)*invn
        edepssum = edepssum + dtr
        dtr2 = edepdtr2(i)*invn
        sigmadtr = (dtr2-dtr**2)*invn
        sigmadtr = sqrt(max(sigmadtr,0.0))
        if (i.eq.1) sigmaeprim = sigmadtr ! stores uncertainty for Eprim separately
        write(out,'(1x,i3,3x,a5,1x,es12.5,1x,es8.1)')
     &   i,Edepnames(i),dtr,2.0*sigmadtr
        enddo
        write(out,'(a)')
     & '#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
        write(out,'(a,1x,es12.5,1x,a,1x,es12.5)')
     &   'Sum of specific depositions',edepssum,', fraction of edep',
     &   edepssum/(edep(matdet)*invn)

      numer = edepdtr(1)
      denom = numer + edepdtr(2)
      Kdfactors(1) = numer / denom
      numer = denom
      denom = denom + edepdtr(3)
      Kdfactors(2) = numer / denom
      numer = denom
      denom = denom + edepdtr(4)
      Kdfactors(3) = numer / denom
      numer = denom
      denom = denom + edepdtr(5)
      Kdfactors(4) = numer / denom
      numer = denom
      denom = denom + edepdtr(6)
      Kdfactors(5) = numer / denom
      Kdfactors(6) = edepdtr(1) / edep(matdet) ! a simplified expression for the 'complete' Kdiaph

      sigma = Kdfactors(6) * sqrt((sigmaeprim/edepdtr(1)*invn)**2+
     &  (sigmadetmat/edep(matdet)*invn)**2)

        ! Massimo writes output with correction factors
       write(out,'(a)')
       write(out,'(a)')
     & '#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
       write(out,'(a)') '# Report on correction factors'
       write(out,'(a,1x,i3)') '# Limited to detection material',matdet
       write(out,'(a)')
     & '# See Burns and Kessler, Phys. Med. Biol. 2009 for def.s'
       write(out,'(a)')
     & '#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
       write(out,'(a)')
     &   '####   counter   Kfactor'
        do i=1,maxk ! cycle through the corr factors
       write(out,'(1x,i3,3x,a5,1x,es12.5)')
     &   i,Kdfacnames(i),Kdfactors(i)
        enddo
       write(out,'(a)')
     & '#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
       write(out,'(a,1x,a,1x,a,1x,es8.1)') '# +/- 2*sigma on',
     &  Kdfacnames(maxK),'=',2.0*sigma

      ! Diaphragm correction factors calculations ended

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
      use DIAPHmod
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

      ! Clear all counters prepared by Massimo
      ! integer i
c$$$      do i=1,maxk
c$$$        edepdtr(i) = 0.0
c$$$        edepdtr2(i) = 0.0
c$$$      enddo
      edepdtr(1:maxk) = 0.0
      edepdtr2(1:maxk) = 0.0
      ! End of diaphragm correction factors counters;

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
      use DIAPHmod
      use PENELOPE_mod, only : nmat
      use dumpMod
      implicit none
      integer mode

      if (mode.eq.1) then
        if (dumpfile.lt.0) return  ! No dump file open
        write(dumpfile) edep(1:nmat),edep2(1:nmat),edepdtr(1:maxk),
     &  edepdtr2(1:maxk),edptmpdiaph(1:maxK)
      else
        if (restartfile.lt.0) return  ! No restart file open
        read(restartfile) edep(1:nmat),edep2(1:nmat),edepdtr(1:maxk),
     &  edepdtr2(1:maxk),edptmpdiaph(1:maxK)
      endif
      end


!>>>> End Of File >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
