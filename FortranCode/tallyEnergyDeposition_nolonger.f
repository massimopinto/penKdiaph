!*******************************************************************
!*                    TALLY ENERGY DEPOSITION                      *
!*                                                                 *
!* Short description:                                              *
!*   Calculates the energy deposited in each material.             *
!*   The so-called collision estimator, which is unambiguously     *
!*   defined in the PENELOPE transport algorithm, is used.         *
!*******************************************************************


      subroutine EDPtally(mode,eloss)
!*******************************************************************
!*    Input:                                                       *
!*      mode -> Identifies the state of the calling procedure      *
!*      eloss -> energy deposition                                 *
!*******************************************************************
      implicit none
      integer mode
      real*8 eloss

      integer*4 kpar,ibody,mat,ilb
      real*8 e,x,y,z,u,v,w,wght
      common/track/e,x,y,z,u,v,w,wght,kpar,ibody,mat,ilb(5)
      logical active
      integer maxmat,detmat,maxk
      parameter (maxmat=10)
      parameter (maxk=6) ! number of correction factors to estimate
      real*8 edptmp,edep,edep2,unclimit,edepdtr,edepdtr2,
     &      edptmpdiaph 

      common /scoedp/ edptmp(maxmat),edep(maxmat),edep2(maxmat),
     &    unclimit,detmat,active,
     &    edepdtr(maxk),edepdtr2(maxk),edptmpdiaph(maxK)! Massimo
                ! energy deposited from transmitted photons
      integer i

      if (.not.active) return

      if (mode.le.0.and.mat.gt.0) then        ! There's energy to be deposited
        ! negative mode values correspond to interaction events (-icol())
        ! but also to the retrieval of a particle from the stack (-99)
        ! or absorption (-98); see penEasy/README.txt file for info
        ! if mode=0 a new particle has been created and added to the stack
        ! this can be a primary as well as a secondary

        ! devi verificare se il caso mode.eq.99 necessita di quanto segue
        ! i casi di mode.eq.-icol sono ok
        ! anche il caso mode.eq.98 (particle absorbed) potrebbe necessirare di trattamento diverso
        ! ponendo la condizione if if((mat.eq.detmat).and.(mode.gt.96)) qui sotto non è cambiato il risultato
        
        ! Writing temporary counters for selected type of events
        
        if(mat.eq.detmat) then ! if you are scoring in the detection material
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
          
          edptmp(mat) = edptmp(mat)+eloss*wght  ! Update temporary E counter for general deposition
          ! when mode.le.0 the energy is to be ADDED to the counters
        

      else if (mode.eq.6) then        ! End-of-history bookkeeping
                  
        do i=1,maxmat
          if (edptmp(i).gt.0.0) then  ! Transfer temporary counter to mean and variance
            
            edep(i)   = edep(i) +edptmp(i)
            edep2(i)  = edep2(i)+edptmp(i)**2
            edptmp(i) = 0.0           ! Clear counter to start a fresh history
            ! this stuff used to be here at the end of the history by pure accident.
            ! Do not zero the edptmp(i) counter yet
            !if(i.eq.detmat) then ! if you are scoring in the detection material
            !    ! ask what is the origin of this event
            !    if (ilb(5).eq.6) then ! primary photon all the way until scoring volume
            !    edepdtr(1)  = edepdtr(1)  +edptmp(i) ! scores primary photon's "Eprim"
            !    edepdtr2(1) = edepdtr2(1) +edptmp(i)
            !    else if (ilb(5).eq.2) then      ! photon transmitted through diaphragm
            !    edepdtr(2)  = edepdtr(2)  +edptmp(i) ! scores transmitted photon's "Edtr"
            !    edepdtr2(2) = edepdtr2(2) +edptmp(i)
            !    else if (ilb(5).eq.3) then
            !    edepdtr(3)  = edepdtr(3)  +edptmp(i) ! scores scattered photon's "Edsc"
            !        ! whereby scattering occurred in the diaphragm and not in air, later on.
            !    edepdtr2(3) = edepdtr2(3) +edptmp(i)
            !    endif !  and so on..here you can add other ILB(5) conditions
            !    endif
            endif
        enddo
        
        do i=1,maxK  ! energy counters specific for each physical process
            if (edptmpdiaph(i).gt.0.0) then  ! Transfer temporary specific counter to mean and variance
                edepdtr(i) = edepdtr(i) + edptmpdiaph(i)
                edepdtr2(i) = edepdtr2(i) + edptmpdiaph(i)**2
                edptmpdiaph(i) = 0.0
                ! Clears counter to start a fresh history
            endif
        enddo
      endif
      end


      subroutine EDPreport(n,cputim,uncdone)
!*******************************************************************
!*    Input:                                                       *
!*      n -> no. of histories simulated                            *
!*      cputim -> elapsed CPU time                                 *
!*    Output:                                                      *
!*      uncdone -> 2 if uncert reached, 1 if not defined, 0 else   *
!*******************************************************************
      implicit none
      integer uncdone
      real*8 n,cputim

      integer*4 seed1,seed2
      common/rseed/seed1,seed2
      logical active
      integer maxmat,detmat,maxk
      parameter (maxmat=10)
      parameter (maxk=6) ! number of energy deposition types to be estimated
      real*8 edptmp,edep,edep2,unclimit,edepdtr,edepdtr2,edptmpdiaph  ! Massimo 
      common /scoedp/ edptmp(maxmat),edep(maxmat),edep2(maxmat),
     &         unclimit,detmat,active,
     &         edepdtr(maxk),edepdtr2(maxk),edptmpdiaph(maxK)! Massimo
      integer i,out,finduf,error
      real*8 q,q2,sigma,eff,uncert,invn,dtr,dtr2,sigmadtr,edepssum,
     &    sigmadetmat,sigmaeprim,Kdfactors(6),num,den
      
      ! names of specific energy counters and K_correction factors
      ! Correction factors, see Burns and Kessler 2009
      character*5 Edepnames(maxK),Kdfacnames(maxK)
      Edepnames(1) = 'Eprim'
      Edepnames(2) = 'Edtr'
      Edepnames(3) = 'Edsc'
      Edepnames(4) = 'Edfl'
      Edepnames(5) = 'Edbr'
      Edepnames(6) = 'Edel'
      
      !character*5 Kdnames(5)
      Kdfacnames(1) = 'Kdtr' ! photons transmitted through diaphragm
      Kdfacnames(2) = 'Kdsc' ! coherent or incoherent scatter from the diaphragm
      Kdfacnames(3) = 'Kdfl' ! fluorescence photons in the diaphragm
      Kdfacnames(4) = 'Kdbr' ! Bremsstrahlung photons in the diaphragm
      Kdfacnames(5) = 'Kdel' ! Electrons exiting from the diaphragm
      Kdfacnames(6) = 'Kdiap' ! this is the product of all Kfactors = Eprim/Etot

      uncdone = 1
      if (.not.active) return

      call EDPdump(1)  ! Sim dump file

      ! Prepare output files:
      out = finduf()
      open(out,file='tallyEnergyDeposition.dat',iostat=error)
      if (error.ne.0) then
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
     & '#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      write(out,'(a)') '# [SECTION REPORT ENERGY DEPOSITION]'
      write(out,'(a)') '#  All energy deposition values in eV'
      write(out,'(a)')
     & '#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      write(out,'(a)')
     &   '# Mat :   E(eV)  :   +-2sig'
      write(out,'(a)')
      invn = 1.0/n
      do i=1,maxmat
        q  = edep(i)*invn
        q2 = edep2(i)*invn
        sigma = (q2-q**2)*invn
        sigma = sqrt(max(sigma,0.0))
        if (i.eq.detmat) sigmadetmat = sigma !  detector material
        write(out,'(1x,i3,1x,es12.5,1x,es8.1)') 
     &   i,q,2.0*sigma
       enddo
       
       ! Massimo writes output with selective energy deposition counters
       write(out,'(a)')
       write(out,'(a)')
     & '#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
       write(out,'(a)') '# Report on selective energy deposition modes'
       write(out,'(a,1x,i3)') '# Limited to detection material',detmat
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
     &   edepssum/(edep(detmat)*invn)
      
      num = edepdtr(1)
      den = num + edepdtr(2)
      Kdfactors(1) = num / den
      num = den
      den = den + edepdtr(3)
      Kdfactors(2) = num / den
      num = den
      den = den + edepdtr(4)
      Kdfactors(3) = num / den
      num = den
      den = den + edepdtr(5)
      Kdfactors(4) = num / den
      num = den
      den = den + edepdtr(6)
      Kdfactors(5) = num / den
      Kdfactors(6) = edepdtr(1) / edep(detmat) ! a simplified expression for the 'complete' Kdiaph
      
      sigma = Kdfactors(6) * sqrt((sigmaeprim/edepdtr(1)*invn)**2+
     &  (sigmadetmat/edep(detmat)*invn)**2)

      
        ! Kdfactors(1) = (edepdtr(1))/((edepdtr(1)+edepdtr(2)))
        ! Kdfactors(2) = (edepdtr(1)+edepdtr(2))/
        !&  (edepdtr(1)+edepdtr(2)+edepdtr(3))
        ! Kdfactors(3) = (edepdtr(1)+edepdtr(2)+edepdtr(3))/
        !&  (edepdtr(1)+edepdtr(2)+edepdtr(3)+edepdtr(4))
        ! Kdfactors(4) = (edepdtr(1)+edepdtr(2)+edepdtr(3)+edepdtr(4)) /
        !&  (edepdtr(1)+edepdtr(2)+edepdtr(3)+edepdtr(4)+edepdtr(5))
        !Kdfactors(5) = (edepdtr(1)+edepdtr(2)+edepdtr(3)+edepdtr(4) 
        !&  + edepdtr(5))/(edepdtr(1)+edepdtr(2)+edepdtr(3)+edepdtr(4)
        !&  + edepdtr(5)+edepdtr(6))
        !
        ! Kdfactors(6) = Kdfactors(1)*Kdfactors(2)*Kdfactors(3)*
        !& Kdfactors(4)*Kdfactors(5)
      
        ! Massimo writes output with correction factors
       write(out,'(a)')
       write(out,'(a)')
     & '#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
       write(out,'(a)') '# Report on correction factors'
       write(out,'(a,1x,i3)') '# Limited to detection material',detmat
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
     
      ! Evaluate rel. uncertainty:
      uncert = 200.0
      q  = edep(detmat)*invn
      q2 = edep2(detmat)*invn
      sigma = (q2-q**2)*invn
      sigma = sqrt(max(sigma,0.0))
      if (q.gt.0.0) uncert = 200.0*sigma/q

      ! Generic report:
      write(out,'(a)')
     & '#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
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
      write(out,'(a)') '#'
      write(out,'(a)') '# Have a nice day.'
      close(out)
      end


      subroutine EDPinitally
!*******************************************************************
!*    Initializes. To be called before TALLY.                      *
!*******************************************************************
      implicit none
      logical active
      integer maxmat,detmat,maxk
      parameter (maxmat=10)
      parameter (maxk=6) ! number of correction factors to estimate
      real*8 edptmp,edep,edep2,unclimit,edepdtr,edepdtr2,edptmpdiaph ! Massimo 
      common /scoedp/ edptmp(maxmat),edep(maxmat),edep2(maxmat),
     &   unclimit,detmat,active,
     &   edepdtr(maxk),edepdtr2(maxk),edptmpdiaph(maxK)! Massimo
                ! energy deposited from transmitted photons
      character*(*) secid,eos
      parameter (secid=
     &'[SECTION TALLY ENERGY DEPOSITION v.2012-06-01]')
      parameter (eos='[END OF EDP SECTION]')
      character*80 buffer
      integer error,i

      write(*,*) ' '
      write(*,'(a)')
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      call getline(buffer)
      if (index(buffer,secid).eq.0) then
        write(*,'(a)') 'EDPinitally:ERROR: incorrect section header;'
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
      else if (buffer.eq.'OFF') then
        active = .false.
        write(*, '(a)') '>>>> Tally Energy Deposition is OFF'
        do
          read(*,'(a80)',iostat=error) buffer
          if (error.ne.0) then
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
      read(*,*) detmat
      write(*,'(i3)') detmat

      write(*,'(a)') 'Relative uncertainty (%) requested:'
      read(*,*) unclimit
      write(*,'(1x,es12.5)') unclimit

      ! Clear counters:
      edptmp = 0.0
      edep   = 0.0
      edep2  = 0.0
      ! Clear all counters prepared by Massimo 
      do i=1,maxk
        edepdtr(i) = 0.0
        edepdtr2(i) = 0.0
      enddo 
      ! End of Massimo's counters; 
      read(*,'(a80)') buffer
      if (index(buffer,eos).eq.0) then
        write(*,*) 'EDPinitally:ERROR: End-Of-Section mark not found'
        write(*,'(a,a)') '  expecting to find: ',eos
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif

      call EDPdump(0)  ! Sim restart file

      write(*,'(a)') '>>>> EDP tally initialization finished >>>>'
      end


      subroutine EDPdump(mode)
!*******************************************************************
!*    Dumps into or reads data from a dump file.                   *
!*                                                                 *
!*    Input:                                                       *
!*      mode -> 1 to write dump file, else to read from it.        *
!*******************************************************************
      implicit none
      integer mode

      logical active
      integer maxmat,detmat,maxk
      parameter (maxmat=10)
      parameter (maxk=6) ! number of correction factors to estimate
      real*8 edptmp,edep,edep2,unclimit,edepdtr,edepdtr2,edptmpdiaph ! Massimo 
      common /scoedp/ edptmp(maxmat),edep(maxmat),edep2(maxmat),
     &   unclimit,detmat,active,
     &   edepdtr(maxk),edepdtr2(maxk),edptmpdiaph(maxK)! Massimo
                ! energy deposited from transmitted photons
      character*80 dumpfilen
      integer restartfile,dumpfile
      common /dump/ restartfile,dumpfile,dumpfilen

      if (mode.eq.1) then
        if (dumpfile.lt.0) return  ! No dump file open
        write(dumpfile) edep(1:maxmat),edep2(1:maxmat),edepdtr(1:maxk),
     &  edepdtr2(1:maxk),edptmpdiaph(1:maxK)
      else
        if (restartfile.lt.0) return  ! No restart file open
        read(restartfile) edep(1:maxmat),edep2(1:maxmat),
     &  edepdtr(1:maxk),edepdtr2(1:maxk),edptmpdiaph(maxK)
      endif
      end


!>>>> End Of File >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
