!*******************************************************************
!*                         penEasy  / penKdiaph                    *
!*                                                                 *
!* Short description:                                              *
!*   General-purpose main program for the PENELOPE system.         *
!*   Please refer to the README.txt file for detailed instructions.*
!*                                                                 *
!*   For a list of dependencies from PENELOPE, see                 *
!*   ~/documentation/dependencies.txt                              *
!*                                                                 *
!* Josep Sempau                                                    *
!* email: josep.sempau@upc.es                                      *
!* Universitat Politecnica de Catalunya, Barcelona, Spain          *
!     * SEE COPYRIGHT NOTICE IN FILE README.txt                    *
!                                                                  *
!* Modifications by Massimo Pinto to calculate correction factors  *
!* concerning the effects occurring in the diaphragm of the FACs   *
!* in 2013 based on penEasy for penelope 2011                      *
!* in 2015-2016 for penelope 2014.                                 *
!*******************************************************************

! PENELOPE library:
#include "penelope.f"
#include "pengeom.f"
#include "rita.f"
#include "penvared.f"
! Protons:
#ifdef PROTONS
#include "penh.f"
#endif

! penEasy library:
#include "penaux.F"
#include "penvox.F"
#include "penvr.F"
#include "timing.F"

! Source models:
#include "sourceBoxIsotropicGaussSpectrum.F"
#include "sourcePhaseSpaceFile.F"

! Tallies:
#include "tallyCylindricalDoseDistrib.F"
#include "tallyEnergyDeposition.F"
! Diaphragm Correction Factors from Massimo
#include "tallyDiaphragmCorrections.F"
#include "tallyFluenceTrackLength.F"
! #include "tallyFluencePointDetect.F"  !! Under construction
#include "tallyParticleCurrentSpectrum.F"
#include "tallyParticleTrackStructure.F"
#include "tallyPhaseSpaceFile.F"
#include "tallyPixelImageDetect.F"
#include "tallyPulseHeightSpectrum.F"
#include "tallySpatialDoseDistrib.F"
#include "tallySphericalDoseDistrib.F"
#include "tallyVoxelDoseDistrib.F"


!*******************************************************************
!*    MAIN                                                         *
!*******************************************************************
      program main
      use TRACK_mod, only : mat,e,kpar,ilb
      use ctrsimMod, only : nhist
      use dsmaxMod
      use forcingMod, only : isforcing
      use DIAPHmod
      implicit none
      logical endsim,absorb
      integer*4 ncross,icol,left
      real*8 ds,dsef,de

      write(*,'(a)') '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      write(*,'(a)') '>>> Based on penEasy v.2015-05-30, >>>'
      write(*,'(a)') '>>> This is penKdiaph v.2016-05    >>>'
      write(*,'(a)') '>>> Contributed by Massimo Pinto   >>>'
      write(*,'(a)') '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'

      call init      ! Initialize the PENELOPE/penEasy system and timers

      history: do              ! Each iteration simulates a new history
        nhist = nhist+1.0d0    ! Update history counter
        call cleans            ! Empty the stack, just in case
        call source            ! Put primary particles (from same history) in stack
        call tally(1,nhist)    ! The simulation of this history begins

        particle: do                       ! Each iteration simulates a new particle
          call secpar(left)                ! Retrieve a particle from the stack
          if (left.eq.0) exit particle ! Stack was empty

!****************************************************************
! Now marking particles as they are retrieved from the stack
! Depending on where they are and how they were created
!****************************************************************


!****************************************************************
! Section specific to AIR in the FREE AIR CHAMBER
! either scoring volume or air behind the diaphragm
! note that this is upon retrieval of the particle from the Stack
! it can be either a primary or a secondary particles
!****************************************************************

          if ((mat.eq.1).or.(mat.eq.3)) then  ! Air inside the chamber.
            ! Air surrounding the detection volume (material 1) or is the detection volume itself (material 3)
            ! if ((kpar.eq.1).and.(ilb(1).eq.2).and.(ilb(2).eq.2)) ilb(5)=1
            if ((ilb(1).eq.2).and.(ilb(2).eq.2)) then ! secondary particle descending from a photon
              if (kpar.eq.1) then ! electron
                ilb(5)=1
              else if (kpar.eq.2) then ! photon
                ilb(5)=7 ! Mark that this is not a primary photon anymore.
                ! Note what a free air chamber should really measure.
                ! Qui si potrebbe replicare il lavoro di Pimpinella del 2005, k_fl, k_sc.
                ! con opportuni contatori aggiuntivi
              endif
            endif
          endif

! A secondary particle was generated directly from a *primary photon* in Air
! This is either a photon or an electron.
! the indication ilb(2).eq.2 is perhaps redundant if the primary particles can only be photons.
! subsequent radiative losses with ilb(5)=1 should be substracted from the counters
! of the body in which the electron was first created (not where Bremmsstrahlung occurs)

! Questa routine, così come è, tratta l'aria fuori dalla camera ad aria libera e quella subito dietro
! il diaframma allo stesso modo dell'aria che sta nel volume di misura. E' quello che vuoi?
! Occorre giusto fare attenzione an caso di elettroni generati fuori da CAL.
! Quelli non credo che occorra seguirli con tanto dettaglio, senno' vai troppo piano.
! Potresti magari avere tre tipo di materiale aria: uno esterno e due interni.

!****************************************************************
! Section specific to the events occuring in the DIAPHRAGM
! note that this is upon retrieval of the particle from the Stack
!****************************************************************

          if ((mat.eq.4).and.(ilb(1).gt.1)) then ! Secondary particle in the diaphragm
            ! così però ad ogni generazione si ricontrolla...
            ! puoi definire ilb(5) solo se non è stata definita già.
            ! ovvero controllare if (ilb(5).eq.0) then ilb(5)=....
            if (kpar.eq.2) then ! Photon
                if (ilb(2).eq.2) then ! Fluorescence photon, was photoelectric effect of inner shell
                    ilb(5)=4 ! fluorescence event in the diaphragm
                  else if (ilb(2).eq.1) then
                    ilb(5)=4 ! fluorescence unless...
                    if (ilb(3).eq.4) then
                      ilb(5)=5 !  Bremsstrahlung in the diaphragm
                    endif
                endif
            else if (kpar.eq.1) then ! electron that does not already come from fluorescence or Bremss photons
                if((ilb(5).ne.4).and.(ilb(5).ne.5)) ilb(5)=6
                ! marks yet another type of deposition from electrons who may escape from the diaphragm
                ! and go deposit energy downstream in the collecting region
                ! As at March 20 (2014?) the E_del counter is slightly large, though. Kdel is about 0.994
            endif
          endif

          !!!!!! Codice in vigore fino al 19 Marzo 2013, quando Edsc veniva sempre zero [ilb(5)=3]
          !if ((mat.eq.4).and.(kpar.eq.2)) then ! now starting transport of a photon created in the diaph
          !      if ((ilb(1).eq.2).and.(ilb(2).eq.2)) then ! It's a secondary photon, descending from a photon
          !          ilb(5)=4 ! fluorescence photon unless what happens below
          !          if ((ilb(3).eq.1).or.(ilb(3).eq.2)) ilb(5)=3 ! was coherent of incoherent scatter here in the diaphragm
          !      else if ((ilb(1).gt.2).and.(ilb(2).eq.1)) then ! a photon descending from a secondary electron
          !          ilb(5)=4 ! Fluorescence photon from an electron interaction unless below
          !          if (ilb(3).eq.4) ilb(5)=5 ! Bremsstrahlung in the diaphragm
          !      endif
          !endif
          !!!!! Fine del codice in vigore fino al 19 Marzo. Se rivuoi questo qui, togli il commento


          call locatevox        ! Locate in voxels geometry (info not in PENELOPE's stack)

          call tally(-99,-e)               ! The simulation of this particle begins
          if (absorb()) cycle particle     ! Check particle absorption
          call start                       ! Reset transport mechanics
          call forcing                     ! Set interaction forcing status

          interact: do                     ! Each iteration simulates an interaction
            if (absorb()) cycle particle   ! Check particle absorption

! Following lines will be processed only if compiled with option -DPROTONS:
#ifdef PROTONS
            if (kpar.eq.4) then            ! It's a proton
              call jumph(dsmax(mat),ds)    ! Get distance to next interaction (note: forcing not implemented)
            else                           ! Photon-electron transport
#endif
            if (isforcing) then            ! Interaction forcing is active
              call jumpf(dsmax(mat),ds)    ! Get distance DS until next interaction
            else
              call jump(dsmax(mat),ds)     ! Same without interaction forcing
            endif
#ifdef PROTONS
            endif
#endif

            call stepx(ds,dsef,ncross)     ! Advance up to interaction point or interface
            if (ncross.eq.0) then          ! Moved a distance DS, no interface crossed
              call tally(3,ds)             ! Notify tallies
            else                           ! Moved a distance DSEF, interface found
              call tally(4,dsef)           ! Notify tallies
              if (mat.eq.0) cycle particle ! New material is vacuum => gone
              ! Molta attenzione a questo, perché se hai delle gap vuote, come
              ! nel calorimetro Rx, sarà meglio impostare aria a bassa pressione, piuttosto.
              ! Massimo adds control routine on labeling transmited photon
              ! as discussed via email with Josep Sempau in February 2013
              if ((mat.eq.4).and.(ilb(1).eq.1).and.(kpar.eq.2)) then ! New material is tungsten
                ilb(5)=2 ! label as transmitted photon before any interaction with knockx().
              ! occhio qui a seguire: sperimentazione del 19/05/2016
                ilb(1)=2 ! label as no longer a primary photon. Necessary? over-cautious?
              ! This is now being exploited on line #246 to re-label scattered photons
              ! that were first found inside the diaphragm.
              endif

              call start                   ! New material => reset transport mechanics
              call forcing                 ! Set interaction forcing status
              call russian                 ! Apply Russian roulette
              call splitting               ! Apply particle splitting
              cycle interact
            endif

#ifdef PROTONS
            if (kpar.eq.4) then            ! Proton transport
              call knockh(de,icol)         ! Simulate an interaction
            else                           ! Photon-electron transport
#endif
            if (isforcing) then            ! Interaction forcing (IF) is active
              call knockf(de,icol)         ! Simulate an interaction with IF
            else
              call knock(de,icol)          ! Simulate an interaction but with/o IF
            endif
#ifdef PROTONS
            endif
#endif

! provo ad anticipare qui prima di tally(-int(icol),de), 18 Maggio.
! Quello che vorrei provare a vedere è se lgi eventi di scattering
! non venivano conteggiati perché chiamati solo dopo tally.
!        if((mat.eq.4).and.(kpar.eq.2).and.(ilb(1).eq.1)) then ! *primary photon in tungsten*
!          if ((icol.eq.1).or.(icol.eq.2)) ilb(5)=3
!        endif
            call tally(-int(icol),de)      ! Tally kinetic energy released

            ! AFTER tally (mimicking "cav-INMRI.f", line 507). Why is this after tally(?)
!           if((mat.eq.4).and.(kpar.eq.2).and.(ilb(1).eq.1)) then ! *primary photon in tungsten*
            if((mat.eq.4).and.(kpar.eq.2).and.(ilb(5).eq.2)) then ! *primary photon in tungsten*
              ! this was already labeled as transmitted, earlier up on line #206
              ! now that a call to knock was issued, we re-label the photon, just in case.
              if ((icol.eq.1).or.(icol.eq.2)) ilb(5)=3
            ! ilb(5)=3 per gamma diffusi coerentemente o incoerentemente
            ! This *must* be done here after a call to knock() because Compton
            ! and coherent scatter don't add a new photon to the stack,
            ! but just alter the direction and the energy of the same photon.
            ! if this was a primary photon, such as here, it will stay a primary photon!
            endif

          enddo interact
        enddo particle

        call tally(6,nhist)                ! End-of-history bookkeeping
        if (endsim()) exit history         ! Simulation is finished
      enddo history

      call report(nhist)                   ! Write final report
      end


      subroutine init
!*******************************************************************
!*    Initializes the simulation system.                           *
!*                                                                 *
!*    Comments:                                                    *
!*      -> The history counter is set to >0 in case of a restart   *
!*      -> The order in which init routines are called must not be *
!*         changed, since some of them provide data for others and *
!*         write data to a dump file. This file must be written in *
!*         the same order in which is will be read afterwards.     *
!*******************************************************************
      use ctrsimMod, only : nhist
      use RSEEDcommonMod
      implicit none

      integer restartf
      integer*4 nmatgeo
      real*8 emax,realtime,cputime,rtime,utime,dmem
      real*8 :: mem=0.0         ! Memory usage

      call initime              ! Write date on the screen
      call treset               ! Reset simulation timer to compute init() timing

      call iniconfig(restartf)  ! Simulation config
      nhist = 0.0d0             ! Clear number of simulated histories
      if (restartf.ge.0) read(restartf) nhist,seed1,seed2  ! Simulation restart
      call inisource(emax,dmem) ! Source models
      mem = mem+dmem            ! Used memory
      call inigeo(nmatgeo,dmem) ! Geometry: PENGEOM & penVox
      mem = mem+dmem
      call inipen(emax,nmatgeo) ! PENELOPE init
      call inimassvox(dmem)     ! Init voxels masses (cannot do before PENELOPE inits)
      mem = mem+dmem
      call initally(emax,dmem)  ! Tallies; tallies read from restart file if needed
      mem = mem+dmem
      call iniforce(emax)       ! Interaction forcing
      call inisplit             ! Particle splitting
      call inirussia            ! Russian roulette

      write(*,*) ''
      write(*,*) ''
      write(*,'(a)')

     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      write(*,'(a)') 'init: INITIALIZATION ENDED'
      write(*,'(a,f9.2,a)') 'Elapsed real time:',realtime(),' s'
      write(*,'(a,f9.2,a)') 'Elapsed CPU time :',cputime(),' s'
      write(*,'(a,es7.1,a)') "Memory used by penEasy arrays: ",
     &  mem*1.0e-6,' MB'

      call treset               ! Reset simulation timers

      if (restartf.ge.0) then   ! This is a sim restart
        read(restartf) rtime,utime
        close(restartf)
        call trestart(rtime,utime) ! Re-reset CPU and realtime timers in case of a restart
        write(*,*) ''
        write(*,'(a)') 'Simulation restarted with the following '//
     &                 'data from the previous run:'
        write(*,'(a,f18.0)') '  No. of histories: ',nhist
        write(*,'(a,2(1x,i0))') '  Seed1,Seed2:',seed1,seed2
        write(*,'(a,2(1x,es12.5))') '  RealTime(s),CPUtime(s):',
     &    rtime,utime
      endif

      write(*,'(a)')
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'

      flush(6)                  ! Flush stdout
      end


      subroutine report(n)
!*******************************************************************
!*    Reports final results.                                       *
!*                                                                 *
!*    Input:                                                       *
!*      n -> no. of histories simulated.                           *
!*******************************************************************
      use RSEEDcommonMod
      implicit none
      real*8 n

      integer unc
      real*8 cputime,realtime,nowcpu

      nowcpu = cputime()                 ! Set a reference to measure report timing
      call tallyreport(-1,n,nowcpu,unc)  ! EOS report, each tally will report its data

      write(*,*) ''
      write(*,*) ''
      write(*,'(a)')
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      write(*,'(a)') 'report: SIMULATION ENDED'
      write(*,'(a)')
     & 'Results have been written to the corresponding DAT files.'
      select case(unc)
      case(0)
        write(*,'(a)')
     &   'The requested uncertainty has NOT been reached.'
      case(1)
        continue  ! Uncertainty limit not defined, no message written
      case default
        write(*,'(a)')
     &   'The requested uncertainty has been reached.'
      end select
      write(*,'(a)')
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'

      write(*,*) ''
      write(*,'(a)') 'Last random seeds:'
      write(*,'(2(1x,i0))') seed1,seed2

      write(*,'(a)') 'Elapsed real time (s), excluding init:'
      write(*,'(1x,es12.5)') realtime()

      write(*,'(a)') 'Elapsed CPU time (s), excluding init:'
      write(*,'(1x,es12.5)') nowcpu

      write(*,'(a)') 'Each report update took (in CPU s):'
      write(*,'(1x,es12.5)') cputime()-nowcpu

      write(*,'(a)') 'No. of histories simulated:'
      write(*,'(1x,f18.0)') n

      if (nowcpu.gt.0.0) then
        write(*,'(a)') 'CPU Speed (histories/s):'
        write(*,'(1x,es12.5)') n/nowcpu
      endif

      call endtime  ! Report date and say goodbye
      end


!*******************************************************************
!*******************************************************************
!*    Source routines start here.                                  *
!*    Source models require:                                       *
!*     i) an initialization routine that must be called by         *
!*        INISOURCE                                                *
!*    ii) a particle generation routine that must be called        *
!*        by SOURCE                                                *
!*******************************************************************
!*******************************************************************

      subroutine inisource(emax,mem)
!*******************************************************************
!*    Init routines for source models.                             *
!*                                                                 *
!*    Output:                                                      *
!*      emax: max source energy (eV).                              *
!*      mem: memory used by source arrays.                         *
!*    Comments:                                                    *
!*      - Note that, if more than one source is defined and        *
!*        active, then EMAX must be the maximum of the value       *
!*        reported by all of them.                                 *
!*******************************************************************
      implicit none
      real*8 emax,mem

      logical active,activepsf
      real*8 dmem

      mem = 0.0
      call BIGSinisrc(active,emax,dmem)
      mem = mem+dmem                   ! Used memory
      call PSFSinisrc(activepsf,emax)

      if (active.and.activepsf) then
        write(*,'(a)') 'inisource:ERROR: PSF source ON is '//
     &    'incompatible with other sources ON.'
        stop
      endif
      end


      subroutine source
!*******************************************************************
!*    Source models.                                               *
!*******************************************************************
      implicit none
      call BIGSsource
      call PSFSource
      end


!*******************************************************************
!*******************************************************************
!*    Tally routines start here.                                   *
!*    Tallies require:                                             *
!*    i) an initialization routine that must be called by INITALLY *
!*    ii) a tally routine that must be called by TALLY             *
!*    iii) a reporting routine that must be called by TALLYREPORT  *
!*    iv) a dump routine to allow restarting of the simulation.    *
!*                                                                 *
!*    Notice that the ordering of the tally initialization routines*
!*    must coincide with the ordering of the corresponding sections*
!*    in the input file.                                           *
!*******************************************************************
!*******************************************************************

      subroutine initally(emax,mem)
!*******************************************************************
!*    Init tallying routines.                                      *
!*                                                                 *
!*    Input:                                                       *
!*      emax: max source energy (eV).                              *
!*    Output:                                                      *
!*      mem: memory used by arrays in tallies.                     *
!*******************************************************************
      implicit none
      real*8 mem,emax,dmem

      mem = 0.0
      call VDDinitally(dmem)
      mem = mem+dmem
      call SDDinitally(dmem)
      mem = mem+dmem
      call CDDinitally(dmem)
      mem = mem+dmem
      call SPDinitally(dmem)
      mem = mem+dmem
      call EDPinitally(dmem)  ! The energy deposition tally routines;
      mem = mem+dmem
      call DIAPHinitally(dmem) ! Diaphragm correction factors
      mem = mem+dmem
      call PHSinitally(dmem)
      mem = mem+dmem
      call PIDinitally(dmem)
      mem = mem+dmem
      call FTLinitally(dmem)
      mem = mem+dmem
      ! call FPDinitally(dmem,emax)  !! Under construction
      ! mem = mem+dmem
      call PSFinitally(emax)
      call PCSinitally(dmem)
      mem = mem+dmem
      call PTSinitally
      end


      subroutine tally(mode,arg)
!*******************************************************************
!*    Tallying routines.                                           *
!*******************************************************************
      implicit none
      integer mode
      real*8 arg

      call VDDtally(mode,arg)  ! Must be 1st tally to be called
      call SDDtally(mode,arg)
      call CDDtally(mode,arg)
      call SPDtally(mode,arg)
      call EDPtally(mode,arg)
      call DIAPHtally(mode, arg) ! For mode-specific energy depositons in matdet.
      call PHStally(mode,arg)
      call PIDtally(mode,arg)
      call FTLtally(mode,arg)
      ! call FPDtally(mode,arg)  !! Under construction
      call PSFtally(mode,arg)
      call PCStally(mode)
      call PTStally(mode,arg)
      end

      subroutine tallyreport(mode,n,cputim,unc)
!*******************************************************************
!*    Calls report routines for all tallies and PSF sources.       *
!*                                                                 *
!*    Input:                                                       *
!*      mode -> -1 if called at end-of-simulation, 0 for dump only,*
!*               1 for tally report only, 2 for both.              *
!*      n -> no. of histories simulated                            *
!*      cputim -> elapsed CPU time                                 *
!*    Output:                                                      *
!*      unc -> 0: uncert not reached, 1: undefined, >1: reached    *
!*    Comments:                                                    *
!*      -> The sequence of calls to the report routines MUST be the*
!*         same used for the corresponding initialization routines *
!*         in INISOURCE and INITALLY (see above). This is to ensure*
!*         proper manipulation of the simulation dump file.        *
!*      -> PSF sources need to report the total num. of particles  *
!*         read to the dump file, which may be needed in case of a *
!*         restart.                                                *
!*******************************************************************
      use dumpMod
      use RSEEDcommonMod
      implicit none
      integer mode,unc,uncdone
      real*8 n,cputim

      integer finduf,error
      real*8 realtime,cputime

      ! Write dump file if needed:
      if (mode.ne.1.and.dumpfile.ne.-1) then  ! Create sim dump file
        dumpfile = finduf()              ! Find a valid unit for the file
        open(dumpfile,file=dumpfilen,status='replace',access='stream',
     &       iostat=error)
        if (error.ne.0) then
          write(*,'(a)')
     &      'tallyreport:ERROR: unable to open dump file.'
          dumpfile = -2                  ! Tells others that file coundn't be opened
        else
          write(dumpfile) n,seed1,seed2  ! Write history state to dump file
        endif
      endif

      ! Write source reports:
      call PSFSreport(mode)              ! PSF source dumps the num. of particles read

      ! Write tally reports:
      unc = 1                            ! Unknown uncertainty status at this point
      call VDDreport(mode,n,cputim,uncdone)
      unc = unc*uncdone
      call SDDreport(mode,n,cputim,uncdone)
      unc = unc*uncdone
      call CDDreport(mode,n,cputim,uncdone)
      unc = unc*uncdone
      call SPDreport(mode,n,cputim,uncdone)
      unc = unc*uncdone
      call EDPreport(mode,n,cputim,uncdone)
      unc = unc*uncdone
      call DIAPHreport(mode,n,cputim,uncdone)
      unc = unc*uncdone
      call PHSreport(mode,n,cputim,uncdone)
      unc = unc*uncdone
      call PIDreport(mode,n,cputim,uncdone)
      unc = unc*uncdone
      call FTLreport(mode,n,cputim,uncdone)
      unc = unc*uncdone
      ! call FPDreport(mode,n,cputim,uncdone)  !! Under construction
      ! unc = unc*uncdone
      call PSFreport(mode,n,cputim)      ! No uncertainty for this tally
      call PCSreport(mode,n,cputim,uncdone)
      unc = unc*uncdone

      if (mode.ne.1.and.dumpfile.ge.0) then
        write(dumpfile) realtime(),cputime() ! Write timings to dump file
        close(dumpfile)
      endif
      end


!>>>> End Of File >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
