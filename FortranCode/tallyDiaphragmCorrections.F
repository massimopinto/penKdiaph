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
!*  by Massimo Pinto of ENEA-INMRI, Italy                          *
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
      use PENELOPE_mod, only : maxmat
      ! use EDPmod ! se davvero serve usarlo qui, non ne son certo.
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

      subroutine DIAPHtally(mode,eloss)
!*******************************************************************
!*    Input:                                                       *
!*      mode:   identifies the state of the calling procedure      *
!*      eloss:   energy deposition                                 *
!*******************************************************************
      use EDPmod ! Non penso che questo occorra davvero in questa routine.
      ! forse potrai alleggerire eliminando EDPmod.
      use DIAPHmod
      use TRACK_mod, only : mat,wght,ilb ! ilb needed for kdiaph corr's
      use PENELOPE_mod, only : nmat
      implicit none
      integer mode
      real*8 eloss

      integer i

      if (.not.active) return
      if (mode.le.0.and.mat.eq.matdet) then ! There's energy to be deposited in the *detector material*
        ! this is copied from EDPtally but is specific to *matdet* alone, which is air in case of a FAC.
        ! Energy will be deposited in appropriate sub.counters, depending on the type of event.
        ! The same energy value is stored by EDPtally in edptmp(mat).
        if ((ilb(5).eq.0).and.(ilb(1).eq.1)) then
        edptmpdiaph(1) = edptmpdiaph(1) +eloss*wght ! scores "Eprim" from electron's energy
        ! a photo- or Compton-electron was created from a primary photon and the initial energy is stored.
        ! the secondary photon is NOT contemplated in the energy count
        ! in agreement with what the free air chamber is to measure.
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

      else if (mode.eq.6) then  ! End-of-history bookkeeping
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

      ! use EDPmod
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
    &    '*********************************************'
       write(*,'(a)')
    &    'DIAPHreport ERROR: cannot open output data file'
       write(*,'(a)')
    &    '*********************************************'
       close(out)  ! Just in case
       return
     endif

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
    close(out)
    end

    
