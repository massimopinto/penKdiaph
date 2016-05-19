!*******************************************************************
!*                    TALLY DIAPGRAGM CORRECTIONS                  *
!*                                                                 *
!* Short description:                                              *
!*   Calculates the energy deposited in the detector material      *
!*   According to specific particle types and events               *
!* including distinction between primary and secondary photons     *
!* and also including scattering from the chamber diaphragm        *
!* Rev:                                                            *
!*   2016-05                                                       *
!*  by Massimo Pinto of ENEA-INMRI, Italy
!*                                                                 *
!* For diaphragm correction factors see                            *
!*  Burns and Kessler, 2009                                        *
!*******************************************************************


!*******************************************************************
!*******************************************************************
!*    Modules section                                              *
!*******************************************************************
!*******************************************************************

      module DIAPHmod
!*******************************************************************
!     Vars for the calculation of diaphragm correction factors  ****
!     Massimo, from the 2011 version of penEasy, with updates   ****
!     fro the penelope version of 2014                          ****
!*******************************************************************
      ! use PENELOPE_mod, only : maxmat
      ! forse ancha maxmat di PENELOPE_mod non occorre.
      use EDPmod !, only: active ! se davvero serve usarlo qui, non ne son certo.
      implicit none
      save
      integer*4, parameter ::  maxk=6  ! number of diaphragm correction factors
      real*8  edepdtr(maxk),edepdtr2(maxk),edptmpdiaph(maxK)
!     Massimo: energy deposited from transmitted photons
      end

!*******************************************************************
!*******************************************************************
!*    Routines section                                             *
!*******************************************************************
!*******************************************************************

      subroutine DIAPHtally(mode,eloss)
!*******************************************************************
!*    Input:                                                       *
!*      mode:   identifies the state of the calling procedure      *
!*      eloss:   energy deposition                                 *
!*******************************************************************
      use EDPmod, only: matdet
      ! Non penso che questo moduo occorre tutto, qui.
      use DIAPHmod
      use TRACK_mod, only : mat,wght,ilb ! ilb needed for kdiaph corr's
      !use PENELOPE_mod, only : nmat ! lo chiedeva EDPTally().
      implicit none
      integer mode
      real*8 eloss

      integer i

      if (.not.active) return

      if ((mode.le.0).and.(mat.eq.matdet)) then ! There's energy to be deposited in the *detector material*
        ! this is copied from EDPtally but is specific to *matdet* alone, which is air in case of a FAC.
        ! Energy will be deposited in appropriate sub.counters, depending on the type of event.
        ! The same energy value is stored by EDPtally in edptmp(mat).

        ! ho commentato queste due righe che seguono in data 18 Maggio 2016.
        ! pensavo corripsondesse ad una sottrazione di energia nel contatore Eprim.ù
        ! Ma così non è stato!
        !if ((ilb(5).eq.0).and.(ilb(1).eq.1)) then
        !    edptmpdiaph(1) = edptmpdiaph(1) +eloss*wght ! scores "Eprim" from electron's energy

        ! a photo- or Compton-electron was created from a primary photon and the initial energy is stored.
        ! the secondary photon is NOT contemplated in the energy count
        ! in agreement with what the free air chamber is to measure.
        !else
        if (ilb(5).eq.1) then
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
      !endif

      else if (mode.eq.6) then  ! End-of-history bookkeeping
        do i=1,maxK             ! energy counters specific for each physical process
            if (edptmpdiaph(i).gt.0.0) then  ! Skip void counters
                ! Transfer temporary counters historic counters edepdtr and edepdtr2
                edepdtr(i) = edepdtr(i) + edptmpdiaph(i)
                edepdtr2(i) = edepdtr2(i) + edptmpdiaph(i)**2
                edptmpdiaph(i) = 0.0
                ! Clears counter to start a fresh history
            endif
        enddo
      endif
      end

      subroutine DIAPHreport(mode,n,cputim,uncdone)
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

      use EDPmod, only: matdet, edep
      use DIAPHmod  !, only: maxK
      use PENELOPE_mod
      use RSEEDcommonMod
      implicit none
      integer mode,uncdone
      real*8 n,cputim

      integer i,out,finduf,errorcode
      ! real*8 q,q2,sigma,eff,uncert
      real*8 sigma,invn
      real*8 dtr,dtr2,sigmadtr,edepssum,sigmadetmat,sigmaeprim,
     & Kdfactors(6),numer,denom

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
      Kdfacnames(6) = 'Kdiap'! this is the product of all Kfactors = Eprim/Etot

      uncdone = 1
      if (.not.active) return

      if (mode.ne.1) call DIAPHdump(1)  ! Sim dump file
      if (mode.eq.0) return           ! Dump only, no report

      ! Prepare output files:
      out = finduf()
      open(out,file='tallyDiaphragmCorrections.dat',iostat=errorcode)
      if (errorcode.ne.0) then
       write(*,*)
       write(*,'(a)')
     &    '***********************************************'
       write(*,'(a)')
     &    'DIAPHreport ERROR: cannot open output data file'
       write(*,'(a)')
     &    '***********************************************'
       close(out)  ! Just in case
       return
      endif

      ! Massimo writes output with selective energy deposition counters

      write(out,'(a)')
     & '#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      write(out,'(a)') '# [SECTION REPORT DIAPHRAGM CORRECTIONS]'
      write(out,'(a)') '# Report on selective energy deposition modes'
      write(out,'(a,1x,i3)') '# Limited to detection material', matdet
      !write(out,'(a)') '# All energy deposition values in eV/hist'
      write(out,'(a)') '# Units are eV per history'
      write(out,'(a)')
     & '#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      write(out,'(a)')
     &   '####   counter E(eV/hist) +/-2sig'

      edepssum = 0.0
      dtr = 0.0
      dtr2 = 0.0
      invn = 1.0/n ! normalize to number of histories.

      do i=1,maxk ! cycle through the energy deposition counters
        dtr  = edepdtr(i)*invn
        edepssum = edepssum + dtr ! adds up the energy deposition counters
        ! this sum must be equal to the 'general' edep counter for the
        ! detector material, dealt with in tallyEnergyDeposition.
        dtr2 = edepdtr2(i)*invn
        sigmadtr = (dtr2-dtr**2)*invn
        sigmadtr = sqrt(max(sigmadtr,0.0))
        if (i.eq.1) sigmaeprim = sigmadtr !Stores uncertainty for Eprim separately
        ! this is useful if you wish to monitor the uncertainty on the energy
        ! deposited by primary photons alone, as a criteria to halt the sim.
        ! it is not implemented at this time, yet.
        write(out,'(1x,i3,3x,a5,1x,es12.5,1x,es8.1)')
     &   i,Edepnames(i),dtr,2.0*sigmadtr
      enddo
      write(out,'(a)')
     & '#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
!      write(out,'(a,1x,es12.5,1x,a,1x,es12.5)')
!     &   '# Sum of specific depositions',edepssum,', fraction of edep',
!     &   edepssum/(edep(matdet)*invn)
      write(out,'(a,1x,es12.5)')
     &  '# Sum of Exxx counters', edepssum

      write(out,'(a,9x,es12.5)')
     &  '# edep(matdet)', edep(matdet)*invn
      write(out,'(a,5x,es12.5)')
     &  '# Fraction of edep', edepssum/(edep(matdet)*invn)
!      write(*,'(a)')
!     &    '#>>>> with edep being the general energy counter in detector'

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
      Kdfactors(6) = edepdtr(1) / edep(matdet)
      ! a simplified expression for the 'complete' Kdiaph correction factor.

      sigma = Kdfactors(6) * sqrt((sigmaeprim/edepdtr(1)*invn)**2+
     &  (sigmadetmat/edep(matdet)*invn)**2)

      ! Writes output with diaphram correction factors
      write(out,'(a)')
      write(out,'(a)')
     & '#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      write(out,'(a)') '# Report on diaphragm correction factors'
      write(out,'(a,1x,i3)') '# Limited to detection material',matdet
      write(out,'(a)')
     & '# See Burns and Kessler, Phys. Med. Biol. 2009 for def.s'
      write(out,'(a)')
     & '#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      write(out,'(a)')
     &   '####   Kfactor Value'
      do i=1,maxk ! cycle through the corr factors
        write(out,'(1x,i3,3x,a5,2x,es12.5)')
     &   i,Kdfacnames(i),Kdfactors(i)
      enddo
      write(out,'(a)')
     & '#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      write(out,'(a,1x,a,1x,a,1x,es8.1)') '# +/- 2*sigma on',
     &  Kdfacnames(maxK),'=',2.0*sigma

! Generic report, taken from tallyEnergyDeposition.dat
      write(out,'(a)') ' '
      write(out,'(a)') '# Performance report'
      write(out,'(a)') '#   No. of histories simulated [N]:'
      write(out,'(a,f18.0)') '#   ',n
      write(out,'(a)') '#   CPU time [t] (s):'
      write(out,'(a,es12.5)') '#   ',cputim
      if (cputim.gt.0.0) then
        write(out,'(a)') '#   Speed (histories/s):'
        write(out,'(a,es12.5)') '#   ',n/cputim
      endif

    ! Diaphragm correction factors calculations ended
      close(out)
      end

      subroutine DIAPHinitally(mem)
!*******************************************************************
!*    Initializes. To be called before TALLY.                      *
!*                                                                 *
!*    Output:                                                      *
!*      mem: Memory used by tally arrays                           *
!*******************************************************************
      ! use EDPmod
      ! probabilmente EDPmod non occorre. Lì c'è
      ! la variabie unclimit che non stiamo usando in questo modulo
      ! per adesso, almeno.
      ! Potrebbe esser necessario impostare un unclimit specifico per il diaframma
      ! se vorrai arrivare ad una determinata incertezza su un particolare fattore.

      use DIAPHmod
      use dataTypesMod
      use PENELOPE_mod, only : maxmat
      implicit none
      real*8 mem

      character*(*) secid,eos
      parameter (secid=
     &  '[SECTION DIAPHRAGM CORRECTIONS v.2016-05-12]')
      parameter (eos='[END OF DIAPH SECTION]')
      character*80 buffer,upcase
      integer errorcode

      write(*,*) ' '
      write(*,'(a)')
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      call getline(buffer,0)
      if (index(buffer,secid).eq.0) then
        write(*,'(a)') 'DIAPHinitally:ERROR: incorrect section header;'
        write(*,'(a,a)') '  expecting to find: ',secid
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif
      write(*,'(a)') secid

      ! Qui ci sarebbe una parte sul memory usage che forse non
      !  è necessario inserire

      ! Segue una parte sulla lettura delle impostazioni di input
      ! (o del file .in) e la relativa scrittura nello standard output
      ! (o nel file .out).
      read(*,'(a80)') buffer
      buffer = adjustl(upcase(buffer))
      buffer = buffer(1:scan(buffer,' ')) ! Clip at 1st blank
      if (buffer.eq.'ON') then
        active = .true.
      else if (buffer.eq.'OFF') then
        active = .false.
        write(*, '(a)') '>>>> Tally Diaphragm Corrections is OFF >>>>'
        do
          read(*,'(a80)',iostat=errorcode) buffer
          if (errorcode.ne.0) then
            write(*,'(a,a,a)') 'DIAPHinitally:ERROR: ',
     &       'Unable to find End-Of-Section mark: ',eos
            stop
          endif
          if (index(buffer,eos).ne.0) return
        enddo
      else
        write(*,'(a)')
     &    'DIAPHinitally:ERROR: expecting to find ON or OFF'
        write(*,'(a)') 'found instead:'
        write(*,'(a)') buffer
        stop
      endif

      write(*,'(a)') 'Detection material is:'
      read(*,*) matdet
      write(*,'(1x,i0)') matdet

      !write(*,'(a)') 'Relative uncertainty (%) requested:'
      !read(*,*) unclimit
      !write(*,'(1x,es12.5)') unclimit

      read(*,'(a80)') buffer
      if (index(buffer,eos).eq.0) then
        write(*,'(a)')
     &    'DIAPHinitally:ERROR: End-Of-Section mark not found'
        write(*,'(a,a)') '  expecting to find: ',eos
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif


      ! Clear event-sepcific counters:
      edepdtr(1:maxk) = 0.0
      edepdtr2(1:maxk) = 0.0
      ! End of diaphragm correction factors counters;

      call DIAPHdump(0)  ! Sim restart file

      write(*,'(a)') '>>>> DIAPHinitally initialization finished >>>>'
      end

      subroutine DIAPHdump(mode)
!*******************************************************************
!*    Dumps into or reads data from a dump file.                   *
!*                                                                 *
!*    Input:                                                       *
!*      mode:   1 to write dump file, else to read from it.        *
!*******************************************************************
      ! use EDPmod, only : edep,edep2
      use DIAPHmod
      ! use PENELOPE_mod, only : nmat
      use dumpMod
      implicit none
      integer mode

      if (mode.eq.1) then
        if (dumpfile.lt.0) return  ! No dump file open
        write(dumpfile) edepdtr(1:maxk),edepdtr2(1:maxk),
     &   edptmpdiaph(1:maxK)
      else
        if (restartfile.lt.0) return  ! No restart file open
        read(restartfile) edepdtr(1:maxk),edepdtr2(1:maxk),
     &    edptmpdiaph(1:maxK)
      endif
      end

!>>>> End Of File >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
