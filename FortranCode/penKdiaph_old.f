!*******************************************************************
!*                         penEasy                                 *
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
!* SEE COPYRIGHT NOTICE IN FILE README.txt                         *
!*******************************************************************

!*******************************************************************
!*    Includes                                                     *
!*******************************************************************
      ! PENELOPE routines:
      include 'penelope.f'
      include 'pengeom.f'
      include 'rita.f'

      ! penEasy libraries:
      include 'penaux.f'
      include 'penpatch.f'
      include 'penvox.f'
      include 'penvr.f'
      include 'timing.f'

      ! Source models (see documentation for a detailed description):
      include 'sourceBoxIsotropicGaussSpectrum.f'
      include 'sourcePhaseSpaceFile.f'

      ! Tallies (see documentation for a detailed description):
      include 'tallyCylindricalDoseDistrib.f'
      include 'tallyEnergyDeposition.f'
      include 'tallyFluenceTrackLength.f'
      include 'tallyParticleCurrentSpectrum.f'
      include 'tallyParticleTrackStructure.f'
      include 'tallyPhaseSpaceFile.f'
      include 'tallyPulseHeightSpectrum.f'
      include 'tallySpatialDoseDistrib.f'
      include 'tallySphericalDoseDistrib.f'
      include 'tallyVoxelDoseDistrib.f'


!*******************************************************************
!*    MAIN                                                         *
!*******************************************************************
      program main
      implicit none
      integer*4 kpar,ibody,mat,ilb
      real*8 e,x,y,z,u,v,w,wght
      common/track/e,x,y,z,u,v,w,wght,kpar,ibody,mat,ilb(5)
      logical endsim,absorb
      integer*4 ncross,icol,left
      real*8 n,ds,dsef,de,dsmax

      write(*,'(a)')
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      write(*,'(a)')
     & '>>>> This is penEasy v.2012-06-01 >>>>'
      write(*,'(a)')
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'

      call init(n) ! Initialize the PENELOPE/penEasy system and timers;
                   !   history counter n is set to 0 if fresh start, >0 if restart

      history: do              ! Each iteration simulates a new history
        n = n+1.0d0            ! Update history counter
        call noverflow(n)      ! Check that N does not overflow
        call cleans            ! Empty the stack, just in case
        call tally(1,n)        ! The simulation of this history begins
        call source(n)         ! Put primary particles (from the same history) in stack
                               ! in so doing, it calls tally(0,n)
           
        particle: do                       ! Each iteration simulates a new particle
          call secpar(left)
          ! Retrieve a particle from the stack
          if (left.eq.0) exit particle     ! Stack was empty
          if ((mat.eq.1).or.(mat.eq.3)) then ! transport of particles in Air
            if((kpar.eq.1).and.(ilb(1).eq.2).and.(ilb(2).eq.2)) ilb(5)=1
            ! ilb(5)=1 is a label for primary depositions in Air.
            ! subsequent radiative losses with ilb(5)=1 should be substracted from the counters
            ! of the body in which the electron was first created (not where Bremmsstrahlung occurs)
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
          if ((mat.eq.4).and.(ilb(1).gt.1)) then ! Secondary particle in the diaphragm
            if (kpar.eq.2) then ! Photon
                if (ilb(2).eq.2) then ! Fluorescence photon, was photoelectric effect of inner shell
                    ilb(5)=4 ! fluorescence 
                else if (ilb(2).eq.1) then
                    ilb(5)=4 ! fluorescence unless...
                    if (ilb(3).eq.4) ilb(5)=5 !  Bremsstrahlung in the diaphragm
                !else if ((ilb(2).eq.1).and.(ilb(3).eq.4)) then
                !    ilb(5)=5 !  Bremsstrahlung in the diaphragm
                endif
            else if (kpar.eq.1) then ! electron that does not already come from fluorescence or Bremss photons
                if((ilb(5).ne.4).and.(ilb(5).ne.5)) ilb(5)=6
                ! marks yet another type of deposition from electrons who may escape from the diaphragm
                ! and go deposit energy downstream in the collecting region
                ! As at March 20 the E_del counter is slightly large, though. Kdel is about 0.994
            endif
          endif
          call tally(-99,-e)               ! The simulation of this particle begins
                                           ! kinetic energy is subtracted from counters
          if (absorb()) cycle particle     ! Check particle absorption
          call start                       ! Reset transport mechanics
          call forcing                     ! Set interaction forcing status
       
          interact: do                     ! Each iteration simulates an interaction
            if (absorb()) cycle particle   ! Check particle absorption
            call jumpx(dsmax(),ds)         ! Get distance DS until next interaction
            call stepx(ds,dsef,ncross)     ! Advance up to interaction point or interface
            if (ncross.eq.0) then
              call tally(3,ds)             ! Moved a distance DS, no interface crossed
            else
              call tally(4,dsef)           ! Moved a distance DSEF, interface found
              if (mat.eq.0) cycle particle ! New material is vacuum => gone
            ! Massimo adds control routine on labeling transmited photon
            ! as discussed via email with Josep Sempau in February 2013
              if ((mat.eq.4).and.(ilb(1).eq.1).and.(kpar.eq.2)) ilb(5)=2
            ! a primary photon that passes through the diaphragm BEFORE any interaction
            ! is labeled with ilb(5)=2
            ! Massimo stops
              call start                   ! New material => reset transport mechanics
              call forcing                 ! Set interaction forcing status
              call russian                 ! Apply Russian roulette
              call splitting               ! Apply particle splitting
              cycle interact
            endif
            call knockx(de,icol)           ! Simulate an interaction
            call tally(-int(icol),de)      ! Tally kinetic energy released
            ! this is a call to tally specific for knockx() events.
            ! Incoherent or coherent scattering of primary photons here, AFTER tally (mimicking "cav-INMRI.f", line 507)
            if((mat.eq.4).and.(kpar.eq.2).and.(ilb(1).eq.1)) then
                if ((icol.eq.1).or.(icol.eq.2)) ilb(5)=3 ! ilb(5)=3 per gamma diffusi coerentemente o incoerentemente
                ! This *must* be done here after a call to knock because Compton and coherent scatter don't add
                ! a new photon to the stack, but just alter its directions and energy 
            endif
            
          enddo interact
        enddo particle

        call tally(6,n)                    ! End-of-history bookkeeping
        if (endsim(n)) exit history        ! Simulation is finished
      enddo history

      call report(n)                       ! Write final report
      end


      subroutine init(n)
!*******************************************************************
!*    Initializes the simulation system.                           *
!*                                                                 *
!*    Output:                                                      *
!*      n -> history counter, 0 if a fresh start, >0 if restart.   *
!*******************************************************************
      implicit none
      real*8 n

      integer*4 seed1,seed2
      common/rseed/seed1,seed2
      integer restartf
      integer*4 nmat
      real*8 emax,realtime,cputime,rtime,utime

      call initime              ! Write date on the screen
      call treset               ! Reset simulation timer to compute init() timing

      call iniconfig(restartf)  ! Simulation config
      n = 0.0d0                 ! Clear number of simulated histories
      if (restartf.ge.0) read(restartf) n,seed1,seed2  ! Simulation restart
      call inisource(n,emax)    ! Source models
      call inigeo(nmat)         ! Geometry: PENGEOM & penVox
      call inipen(emax,nmat)    ! PENELOPE
      call initally             ! Tallies; tallies read from restart file if needed
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

      call treset               ! Reset simulation timers

      if (restartf.ge.0) then   ! This is a sim restart
        read(restartf) rtime,utime
        close(restartf)
        call trestart(rtime,utime) ! Re-reset CPU and realtime timers in case of a restart
        write(*,*) ''
        write(*,'(a)') 'Simulation restarted with the following '//
     &                 'data from the previous run:'
        write(*,'(a,f18.0)') '  No. of histories: ',n
        write(*,'(a,2(1x,i0))') '  Seed1,Seed2:',seed1,seed2
        write(*,'(a,2(1x,es12.5))') '  RealTime(s),CPUtime(s):',
     &    rtime,utime
      endif

      write(*,'(a)')
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      end


      subroutine report(n)
!*******************************************************************
!*    Reports final results.                                       *
!*                                                                 *
!*    Input:                                                       *
!*      n -> no. of histories simulated.                           *
!*******************************************************************
      implicit none
      real*8 n

      integer*4 seed1,seed2
      common/rseed/seed1,seed2
      integer unc
      real*8 cputime,realtime,nowcpu

      nowcpu = cputime()        ! Set a reference to measure report timing
      call tallyreport(n,nowcpu,unc)     ! Each tally will report its data

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

      subroutine inisource(n,emax)
!*******************************************************************
!*    Init routines for source models.                             *
!*                                                                 *
!*    Input:                                                       *
!*      n -> history counter, >0 if restart.                       *
!*    Output:                                                      *
!*      emax -> max source energy (eV).                            *
!*    Comments:                                                    *
!*      -> Note that, if more than one source is defined and       *
!*         active, then emax must be the maximum of the value      *
!*         reported by all of them.                                *
!*******************************************************************
      implicit none
      real*8 n,emax
      logical active,activepsf

      call BIGSinisrc(active,emax)
      call PSFinisrc(n,activepsf,emax)  ! PSF needs to know N in case of restart

      if (active.and.activepsf) then
        write(*,'(a)')
     &    'inisource:ERROR: the PSF and other source models are '//
     &    'incompatible; turn off PSF or the others.'
        stop
      endif
      end


      subroutine source(n)
!*******************************************************************
!*    Source models.                                               *
!*                                                                 *
!*    Input:                                                       *
!*      n -> top history counter.                                  *
!*******************************************************************
      implicit none
      real*8 n

      call BIGSsource(n)
      call PSFsource(n)
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

      subroutine initally
!*******************************************************************
!*    Init tallying routines.                                      *
!*                                                                 *
!*    Comments:                                                    *
!*      -> VDDinitally sets variables that are needed for          *
!*         proper particle transport in voxelized geometries.      *
!*         Therefore, it should not be removed even if the tally   *
!*         is not used.                                            *
!*******************************************************************
      implicit none
      call VDDinitally
      call SDDinitally
      call CDDinitally
      call SPDinitally
      call EDPinitally ! The energy deposition tally routines;
      call PHSinitally
      call FTLinitally
      call PSFinitally
      call PCSinitally
      call PTSinitally
      end


      subroutine tally(mode,arg)
!*******************************************************************
!*    Tallying routines.                                           *
!*                                                                 *
!*    Comments:                                                    *
!*      -> The order in which tallies are called MUST be the same  *
!*         in which the corresponding initialization routines were *
!*         called in INITALLY (above). This is to ensure proper    *
!*         manipulation of the simulation dump file.               *
!*      -> VDDtally sets state variables that are needed for       *
!*         proper particle transport in voxelized geometries.      *
!*         Therefore, it should not be removed even if the tally   *
!*         is not used.                                            *
!*      -> Furthermore, these variables could be used by other     *
!*         tallies and, in consequence, VDDtally should be the     *
!*         FIRST tally to be called.                               *
!*******************************************************************
      implicit none
      integer mode
      real*8 arg

      call VDDtally(mode,arg)  ! Must be 1st tally to be called
      call SDDtally(mode,arg)
      call CDDtally(mode,arg)
      call SPDtally(mode,arg)
      call EDPtally(mode,arg)  ! This is the Energy deposition tally
      call PHStally(mode,arg)
      call FTLtally(mode,arg)
      call PSFtally(mode,arg)
      call PCStally(mode)
      call PTStally(mode,arg)
      end


      subroutine tallyreport(n,cputim,unc)
!*******************************************************************
!*    Calls report routines for all tallies.                       *
!*                                                                 *
!*    Input:                                                       *
!*      n -> no. of histories simulated                            *
!*      cputim -> elapsed CPU time                                 *
!*    Output:                                                      *
!*      unc -> 0: uncert not reached, 1: undefined, >1: reached    *
!*******************************************************************
      integer unc,uncdone
      real*8 n,cputim

      integer*4 seed1,seed2
      common/rseed/seed1,seed2
      character*80 dumpfilen
      integer restartfile,dumpfile
      common /dump/ restartfile,dumpfile,dumpfilen
      integer finduf,error
      real*8 realtime,cputime

      ! Write dump file if needed:
      if (dumpfile.ne.-1) then           ! Create sim dump file
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

      ! Write partial reports to corresponding data files:
      unc = 1                            ! Unknown uncertainty status at this point
      call VDDreport(n,cputim,uncdone)
      unc = unc*uncdone
      call SDDreport(n,cputim,uncdone)
      unc = unc*uncdone
      call CDDreport(n,cputim,uncdone)
      unc = unc*uncdone
      call SPDreport(n,cputim,uncdone)
      unc = unc*uncdone
      call EDPreport(n,cputim,uncdone)
      unc = unc*uncdone
      call PHSreport(n,cputim,uncdone)
      unc = unc*uncdone
      call FTLreport(n,cputim,uncdone)
      unc = unc*uncdone
      call PSFreport(n,cputim)           ! No uncertainty for this tally
      call PCSreport(n,cputim,uncdone)
      unc = unc*uncdone
      call PTSreport                     ! No arguments for this tally

      if (dumpfile.ge.0) then
        write(dumpfile) realtime(),cputime() ! Write timings to dump file
        close(dumpfile)
      endif
      end


!>>>> End Of File >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
