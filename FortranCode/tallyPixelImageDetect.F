!*******************************************************************
!*                          TALLY                                  *
!*               PIXELATED IMAGING DETECTOR                        *
!*                                                                 *
!* Short description:                                              *
!*   Tallies the pixelated image corresponding to the radiation    *
!*   flux arriving at a detector.                                  *
!* Rev:                                                            *
!*   2015-02-06                                                    *
!*******************************************************************


!*******************************************************************
!*******************************************************************
!*    Modules section                                              *
!*******************************************************************
!*******************************************************************

      module PIDmod
!*******************************************************************
!*    Vars for the PID tally.                                      *
!*******************************************************************
      implicit none
      save
      logical active,isLab,isEspread,isEinteg,isPCM
      integer detmode,repformat,collided,selectCollision
      integer*4 bodydet,matdet,nx,ny,nebin,ijmax,ijkmax
      real*8 lab2det(2,3),r0(3),dx,dy,idx,idy,c0,c1
      real*8 unclimit,ethr,emin,emax,ebin,iebin
      real*8, parameter ::
     &       fwhm2sig2=1.0/(8.0*log(2.0d0)),
     &       npixmax=1.024e9    ! Max num pixels=32k squared (do not exceed int*4)
      real*8 :: memPID=0.0      ! Memory used by arrays (bytes)
      real*8, allocatable ::
     &  image(:),image2(:),edeptmp(:),nlast(:),wghtlast(:)
      end


!*******************************************************************
!*******************************************************************
!*    Routines section                                             *
!*******************************************************************
!*******************************************************************

      subroutine PIDtally(mode,arg)
!*******************************************************************
!*    Input:                                                       *
!*      mode: Identifies the state of the calling routine          *
!*      arg: energy loss (mode<0) or history no. (mode=1)          *
!*******************************************************************
      use PIDmod
      use TRACK_mod
      use ctrsimMod, only : nhist  ! Current history number
      use SECSTcommonMod, only : nsec,ilbs
      implicit none
      integer mode
      real*8 arg

      integer*4 i,j,ij,k,ijk
      real*8 rprima(3),xdet,ydet,espread

      if (.not.active) return

      ! Filter photons according to last collision:
      if (kpar.eq.2.and.selectCollision.ne.0) then  ! Photon filtering is active
        select case (mode)
        case (-1,-2)                 ! Rayleigh (-1) or Compton (-2) interaction
          if (collided.eq.-1) then   ! Not interacted before, unscattered
            collided = -mode         ! Tag as collided (Rayleigh or Compton)
          else
            collided = 9             ! Tag as multicollided
          endif
          return                     ! Detector should be a perfect absorbent, so interactions cannot occur in detector
        case (-98,-97)               ! Particle has been absorbed
          if (collided.ne.selectCollision) return  ! Score only selected collided particle events
        case (-99)                   ! Particle retrieved from the stack
          if (ilb(5).eq.ilb(1)) then ! The particle was source-created (could be a secondary if from a PSF)
            collided = -1            ! Tag as unscattered
          else
            collided = 3             ! Tag as secondary particle (not in a PSF)
          endif
        case (1)                     ! A new history is about to begin
          ilbs(5,1:nsec) = ilbs(1,1:nsec) ! Tag all particles in the stack as source-created particles
        end select
      endif

      ! Score:
      if (mat.eq.matdet.and.mode.le.0) then   ! Energy has been deposited in detection material
        if (arg.eq.0.0) return                ! Nothing to deposit (e.g., elastic scattering)

        if (isLab) then                       ! Slab and Sdet basis vectors coincide, no rotation needed
          xdet = x-r0(1)                      ! Position vector w.r.t. pixels origin, only x&y coords are needed
          ydet = y-r0(2)
        else                                  ! Both a translation and a rotatation are needed
          rprima = (/x,y,z/)-r0               ! Translation
          xdet = dot_product(lab2det(1,:),rprima) ! Rotation to obtain coords in Sdet frame
          ydet = dot_product(lab2det(2,:),rprima)
        endif

        ! Check if particle is inside pixels grid:
        i = xdet*idx+1.0                      ! Pixel index
        if (i.lt.1.or.i.gt.nx) return         ! Out of pixels grid
        j = ydet*idy+1.0                      ! Pixel index
        if (j.lt.1.or.j.gt.ny) return         ! Out of pixels grid
        ij = i+(j-1)*nx                       ! Mapping of i,j into a single index

        if (isEinteg) then                    ! E-integrating mode
          if (nhist.gt.nlast(ij)) then               ! Apply "deferred bookkeeping": transfer temporary counter when new history
            if (isEspread) edeptmp(ij)=espread(edeptmp(ij)) ! Gaussian energy spread
            image(ij)   = image(ij) +edeptmp(ij)     ! Transfer temp counter to permanent
            image2(ij)  = image2(ij)+edeptmp(ij)**2  ! ...and the variance too
            edeptmp(ij) = arg*wght                   ! Reset temp counter
            nlast(ij)   = nhist+0.5                  ! Reset last visit--add 1/2 to prevent round-off problems
          else                                       ! Same history as previous event, continue temp scoring
            edeptmp(ij) = edeptmp(ij)+arg*wght       ! Score in temp counter
          endif

        else if (nhist.gt.nlast(ij)) then     ! Photon modes; visit of a new history occurred
          if (isPCM) then                     ! Photon counting mode (PCM)
            if (isEspread) edeptmp(ij)=espread(edeptmp(ij))  ! Gaussian energy spread
            if (edeptmp(ij).gt.ethr) then                    ! Deposited energy exceeds threshold
              ! AutoNote: Gaussian spread could be reduced to prob above the ETHR (Erf) for this case.
              image(ij)  = image(ij) +wghtlast(ij)           ! Add one count to this pixel
              image2(ij) = image2(ij)+wghtlast(ij)**2
            endif

          else                                ! Photon E discriminating mode
            if (isEspread) then                              ! Apply Gaussian energy spread
              k = (espread(edeptmp(ij))-emin)*iebin+1.0      ! E channel where the a count will be added
            else
              k =         (edeptmp(ij) -emin)*iebin+1.0      ! Same without E spread
            endif
            if (k.gt.0.and.k.le.nebin) then                  ! Is a valid channel
              ijk = ij+(k-1)*ijmax                           ! Map i,j,k into a single index
              image(ijk)  = image(ijk) +wghtlast(ij)         ! Add one count
              image2(ijk) = image2(ijk)+wghtlast(ij)**2      ! Variance
            endif

          endif
          edeptmp(ij)  = arg                       ! Reset temp counter
          nlast(ij)    = nhist+0.5                 ! Reset last visit--add 1/2 to prevent round-off problems
          wghtlast(ij) = wght                      ! Reset weight to current value; shouldn't change for same hist

        else                                  ! Photon modes; same history as previous event, continue temp scoring
          edeptmp(ij)  = edeptmp(ij)+arg           ! In photon modes WGHT is applied to the added count, not to deposited energy

        endif

      endif
      end


      subroutine PIDreport(mode,n,cputim,uncdone)
!*******************************************************************
!*    Input:                                                       *
!*      mode: -1 if called at end-of-simulation, 0 for dump only,  *
!*             1 for tally report only, 2 for both                 *
!*      n: no. of histories simulated                              *
!*      cputim: elapsed CPU time                                   *
!*    Output:                                                      *
!*      uncdone: 2 if uncert reached, 1 if not defined, 0 else     *
!*               Not computed if mode=0                            *
!*******************************************************************
      use PIDmod
      use RSEEDcommonMod
      implicit none
      integer mode,uncdone
      real*8 n,cputim

      character*80 buffer
      integer out,finduf,errorcode
      integer*4 nchan,i,j,k,ij,ijk
      real*8 q,sigma,eff,avesig,maxq,uncert,x,y,xmiddle,ymiddle,invn
      real*8 e,emiddle,espread,fact,fact2,factE,fact2E,idarea

      uncdone = 1
      if (.not.active) return

      ! Dump temp counters:
      do ij=1,ijmax
        if (nlast(ij).lt.0.5) cycle                     ! Skip empty pixels
        if (isEspread) edeptmp(ij)=espread(edeptmp(ij)) ! Gaussian energy spread
        select case (detmode)
        case (1)                                      ! E-integrating mode
          image(ij)   = image(ij) +edeptmp(ij)          ! Transfer temp counter to permanent
          image2(ij)  = image2(ij)+edeptmp(ij)**2
        case (2)                                      ! Photon counting mode (PCM)
          if (edeptmp(ij).gt.ethr) then                 ! Deposited energy exceeds threshold
            image(ij)  = image(ij) +wghtlast(ij)        ! Add one count
            image2(ij) = image2(ij)+wghtlast(ij)**2     ! Variance
          endif
          wghtlast(ij) = 0.0                            ! Reset last weight
        case (3)                                      ! Photon energy discrimination mode
          k = (edeptmp(ij)-emin)*iebin+1.0              ! E channel where count will be added
          if (k.gt.0.and.k.le.nebin) then               ! Is a valid channel
            ijk = ij+(k-1)*ijmax                        ! Map i,j,k into a single index
            image(ijk)  = image(ijk) +wghtlast(ij)      ! Add one count
            image2(ijk) = image2(ijk)+wghtlast(ij)**2   ! Variance
          endif
          wghtlast(ij) = 0.0                            ! Reset last weight
        end select
        ! Reset all:
        edeptmp(ij)  = 0.0                              ! Reset counter
        nlast(ij)    = 0.0                              ! Reset last visited to avoid recounting in next report
      enddo

      if (mode.ne.1) call PIDdump(1)  ! Sim dump file, notice that temp counter and nlast have been cleared before
      if (mode.eq.0) return           ! Dump only, no report

      ! Filename:
      select case (detmode)
        case (1)
          buffer = 'tallyPixelImageDetectEI'
        case (2)
          buffer = 'tallyPixelImageDetectPC'
        case (3)
          buffer = 'tallyPixelImageDetectED'
      end select
      select case (repformat)
        case (1)
          buffer = trim(buffer)//'.dat'
        case (2)
          buffer = trim(buffer)//'-matrix.dat'
        case (3)
          buffer = trim(buffer)//'-binary.dat'
      end select
      out = finduf()
      if (repformat.eq.3) then         ! Open data file in binary format
        open(out,file=buffer,status='replace',access='stream',
     &       iostat=errorcode)
      else
        open(out,file=buffer,iostat=errorcode)
      endif
      if (errorcode.ne.0) then
        write(*,*) ''
        write(*,'(a)')
     &    '*********************************************'
        write(*,'(a)')
     &    'PIDreport:ERROR: cannot open output data file'
        write(*,'(a)')
     &    '*********************************************'
        close(out)  ! Just in case
        return
      endif

      if (repformat.ne.3) then ! Write header only if not in binary mode
        write(out,'(a)')
     &  '#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
        write(out,'(a)')
     &  '# [SECTION REPORT PIXELATED IMAGING DETECTOR]'

        select case (selectCollision)
        case (-1)
          write(out,'(a)')
     &     '# Reporting ONLY charged particles and '//
     &     'unscattered photons contributions.'
          write(out,'(a)') ''
        case (1)
          write(out,'(a)')
     &     '# Reporting ONLY charged particles and '//
     &     '(once) Rayleigh-scattered photons contributions.'
          write(out,'(a)') ''
        case (2)
          write(out,'(a)')
     &     '# Reporting ONLY charged particles and '//
     &     '(once) Compton-scattered photons contributions.'
          write(out,'(a)') ''
        case (3)
          write(out,'(a)')
     &     '# Reporting ONLY charged particles and '//
     &     'secondary photons (not in a PSF source) contributions.'
          write(out,'(a)') ''
        case (9)
          write(out,'(a)')
     &     '# Reporting ONLY charged particles and '//
     &     'multi-scattered photons contributions.'
          write(out,'(a)') '#'
        end select

        select case (detmode)
        case (1)
          write(out,'(a)')
     &     '# Energy integrating mode, units are eV/cm^2 per history.'
        case (2)
          write(out,'(a)')
     &     '# Photon counting mode, units are counts/cm^2 per history.'
          write(out,'(a)') '#'
          write(out,'(a)') '# Energy deposition threshold (eV):'
          write(out,'(a,es12.5)') '# ',ethr
        case (3)
          write(out,'(a)')
     &      '# Photon energy discriminating mode, '//
     &      'units are counts/(cm^2.eV) per history.'
          write(out,'(a)') '#'
          write(out,'(a)') '# No. of energy bins and bin width (eV):'
          write(out,'(a,1x,i0,1x,es12.5)') '# ',nebin,ebin
        end select
        write(out,*) '#'
        write(out,'(a)')
     &    '# No. of pixels and pixel size (cm) in x & y directions:'
        write(out,'(a,1x,i0,1x,es12.5)') '# ',nx,dx
        write(out,'(a,1x,i0,1x,es12.5)') '# ',ny,dy
        write(out,'(a)')
     &    '# (x & y with respect to detector reference frame.)'
        write(out,'(a)') '#'
      endif

      select case (repformat)
      case (1)  ! Columnar
        write(out,'(a)')
     &    '# For plotting purposes, two values per '//
     &    'coordinate are given, namely, '//
     &    'the low end and the middle point of each interval.'
        write(out,'(a)') '#'
        write(out,'(a,$)')
     &    '# xPixelIndex : xLow(cm) : xMiddle(cm) : '//
     &      'yPixelIndex : yLow(cm) : yMiddle(cm) : '
        select case (detmode)
        case (1)
          write(out,'(a)') 'image(eV/cm^2/hist) : +-2sigma'
        case (2)
          write(out,'(a)') 'image(counts/cm^2/hist) : +-2sigma'
        case (3)
          write(out,'(a)') 'Eindex : Elow(eV) : Emiddle(eV) : '//
     &                     'image(counts/cm^2/eV/hist) : +-2sigma'
        end select
      case (2)  ! Matrix format
        write(out,'(a)') '# Pixel x index:'
        write(out,'(a,$)') '#'
        do i=1,nx
          write(out,'(i12,a,$)') i,' '  ! Extra padding space to match data columns
        enddo
        write(out,*) ''                 ! Terminates non-advancing output
      end select                        ! Nothing to do for case3 binary

      ! Write data:
      invn   = 1.0/n
      nchan  = 0
      avesig = 0.0
      maxq   = 0.5*maxval(image)                ! 1/2 of the max score
      idarea = 1.0/(dx*dy)                      ! Inverse pixels area
      fact   = invn*idarea                      ! Factors for E integ and PC modes
      fact2  = invn*idarea**2
      factE  = invn*iebin*idarea                ! Factors for E discrim mode
      fact2E = invn*(iebin*idarea)**2

      do j=1,ny
        if (repformat.eq.1) then                ! Compute coords
          y = dy*(j-1)
          ymiddle = y+dy*0.5
        endif
        do i=1,nx
          if (repformat.eq.1) then
            x = dx*(i-1)
            xmiddle = x+dx*0.5
          endif
          ij = i+(j-1)*nx                       ! Mapping of i,j into a single index

          if (detmode.eq.3) then                ! Photon energy discriminating mode
            do k=1,nebin
              if (repformat.eq.1) then
                e = emin+ebin*(k-1)
                emiddle = e+ebin*0.5
              endif
              ijk = ij+(k-1)*ijmax              ! Map i,j,k into a single index
              q = image(ijk)*factE
              sigma = sqrt(max((image2(ijk)*fact2E-q**2)*invn,0.0))
              if (repformat.eq.1) then          ! Columnar
                write(out,'(3(1x,i5,2(1x,es12.5)),1x,es12.5,1x,es7.1)')
     &            i,x,xmiddle,j,y,ymiddle,k,e,emiddle,q,2.0*sigma
              else                              ! Binary; matrix format is disabled in PED mode
                write(out) q,2.0*sigma
              endif
              if (image(ijk).gt.maxq) then      ! Evaluate average uncertainty for scores above 1/2 max score
                avesig = avesig+(sigma/q)**2
                nchan = nchan+1
              endif
            enddo
            if (repformat.ne.3) write(out,*) '' ! Separate data blocks, except in binary mode

          else                                  ! Energy integrating and PC modes
            q = image(ij)*fact
            sigma = sqrt(max((image2(ij)*fact2-q**2)*invn,0.0))
            select case (repformat)
            case (1)                            ! Columnar
              write(out,'(2(1x,i5,2(1x,es12.5)),1x,es12.5,1x,es7.1)')
     &          i,x,xmiddle,j,y,ymiddle,q,2.0*sigma
            case (2)                            ! Matrix
              write(out,'(1x,es12.5,$)') q
            case (3)                            ! Binary
              write(out) q,2.0*sigma
            end select
            if (image(ij).gt.maxq) then         ! Evaluate average uncertainty for scores above 1/2 max score
              avesig = avesig+(sigma/q)**2
              nchan = nchan+1
            endif
          endif
        enddo

        if (repformat.ne.3) write(out,*) ''     ! Separate data blocks, except in binary mode
      enddo

      uncdone = 0
      if (nchan.gt.0) then
        uncert = 200.0*sqrt(avesig/nchan)
        if (uncert.lt.unclimit) uncdone = 2  ! Uncertainty reached
      else
        uncert = 0.0                    ! Uncertainty assumed not reached when score is nil
      endif

      if (repformat.ne.3) then          ! Generic report, except in binary mode
        write(out,'(a)') ' '
        write(out,'(a)') '# Performance report'
        write(out,'(a)') '#   Memory used (MB):'
        write(out,'(a,es12.5)') '#   ',memPID*1.0e-6
        write(out,'(a)') '#   Random seeds:'
        write(out,'(a,i0)') '#   ',seed1
        write(out,'(a,i0)') '#   ',seed2
        write(out,'(a)') '#   No. of histories simulated [N]:'
        write(out,'(a,f18.0)') '#   ',n
        write(out,'(a)') '#   CPU time [t] (s):'
        write(out,'(a,es12.5)') '#   ',cputim
        if (cputim.gt.0.0) then
          write(out,'(a)') '#   Speed (histories/s):'
          write(out,'(a,es12.5)') '#   ',n/cputim
        endif
        write(out,'(a)')
     &  '#   Average uncertainty (above 1/2 max score) in % [uncert]:'
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
      endif

      close(out)
      end


      subroutine PIDinitally(mem)
!*******************************************************************
!*    Initializes. To be called before TALLY                       *
!*                                                                 *
!*    Output:                                                      *
!*      mem: Memory used by tally arrays.                          *
!*******************************************************************
      use PIDmod
      use dataTypesMod
      use PENELOPE_mod, only : eabs,maxmat
      use CSIMPHcommonMod, only : eabsh       ! Protons
      implicit none
      real*8 mem

      character*(*) secid,eos
      parameter (secid=
     &'[SECTION TALLY PIXELATED IMAGING DETECTOR v.2015-02-06]')
      parameter (eos='[END OF PID SECTION]')
      character*80 buffer,upcase
      integer errorcode
      integer*4 ijwght
      real*8 pixlen
      real*8, parameter :: einf=1.1e9  ! Above max E in PENELOPE

      write(*,*) ' '
      write(*,'(a)')
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      call getline(buffer,0)
      if (index(buffer,secid).eq.0) then
        write(*,'(a)') 'SDDinitally:ERROR: incorrect section header;'
        write(*,'(a,a)') '  expecting to find: ',secid
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif
      write(*,'(a)') secid

      read(*,'(a80)') buffer
      buffer = adjustl(upcase(buffer))      ! Adjust and capitalize
      buffer = buffer(1:scan(buffer,' '))   ! Clip at 1st blank
      if (buffer.eq.'ON') then
        active = .true.
      else if (buffer.eq.'OFF') then
        active = .false.
        write(*, '(a)')
     &    '>>>> Tally Pixelated Imaging Detector is OFF >>>>'
        do
          read(*,'(a80)',iostat=errorcode) buffer
          if (errorcode.ne.0) then
            write(*,'(a,a,a)') 'PIDinitally:ERROR: ',
     &       'Unable to find End-Of-Section mark: ',eos
            stop
          endif
          if (index(buffer,eos).ne.0) return
        enddo
      else
        write(*,'(a)')
     &    'PIDinitally:ERROR: expecting to find ON or OFF'
        write(*,'(a)') 'found instead:'
        write(*,'(a)') buffer
        stop
      endif

      write(*,'(a)') 'Detection material set to:'
      read(*,*) matdet
      write(*,'(1x,i0)') matdet
      if (matdet.le.0.or.matdet.gt.maxmat) then
        write(*,'(a)') 'PIDinitally:ERROR: invalid entry.'
        stop
      endif

      write(*,'(a)') 'Filter photon collisions option:'
      read(*,*) selectCollision
      write(*,'(1x,i0)') selectCollision
      select case (selectCollision)
      case (0)             ! Valid entriy, do nothing
      case (-1,1,2,3,9)  ! Valid entries, check if perfect absorbent
#ifdef PROTONS
        if (eabsh(matdet).lt.einf) then
          write(*,'(a,es12.5)')
     &      'PIDinitally:ERROR: detection material must be a '//
     &      'perfect absorbent; increase absorption energies above ',
     &      einf
          stop
        endif
#endif
        if (minval(eabs(:,matdet)).lt.einf) then
          write(*,'(a,es12.5)')
     &      'PIDinitally:ERROR: detection material must be a '//
     &      'perfect absorbent; increase absorption energies above ',
     &      einf
          stop
        endif
      case default
        write(*,'(a)') 'PIDinitally:ERROR: invalid value.'
        stop
      end select

      write(*,'(a)') 'Pixel size (cm) & No. pixels in x direction:'
      read(*,*) dx,nx
      write(*,'(1x,es12.5,1x,i0)') dx,nx
      if (dx.lt.0.0.or.nx.lt.0) then
        write(*,'(a)') 'PIDinitally:ERROR: invalid value.'
        stop
      endif
      if (dx.eq.0.0.and.nx.eq.0) then
        write(*,'(a)') 'PIDinitally:ERROR: '//
     &    'pixel size and no. pixels are both zero.'
        stop
      endif
      pixlen = dx*nx  ! Length of the pixels grid along x
      if (pixlen.ne.0.0) then
        write(*,'(a)') 'Pixels grid x-size is (cm):'
        write(*,'(1x,es12.5)') pixlen
      endif

      write(*,'(a)') 'Pixel size (cm) & No. pixels in y direction:'
      read(*,*) dy,ny
      write(*,'(1x,es12.5,1x,i0)') dy,ny
      if (dy.lt.0.0.or.ny.lt.0) then
        write(*,'(a)') 'PIDinitally:ERROR: invalid value.'
        stop
      endif
      if (dy.eq.0.0.and.ny.eq.0) then
        write(*,'(a)') 'PIDinitally:ERROR: '//
     &    'pixel size and no. pixels are both zero.'
        stop
      endif
      pixlen = dy*ny  ! Length of the pixels grid along y
      if (pixlen.ne.0.0) then
        write(*,'(a)') 'Pixels grid y-size is (cm):'
        write(*,'(1x,es12.5)') pixlen
      endif

      call setDetectFrame              ! Set detector reference frame (Sdet) and pixels size&num
      write(*,'(a)') 'Number of megapixels total:'
      write(*,'(1x,es12.5)') real(nx)*real(ny)*1.0e-6
      if (real(nx)*real(ny).gt.npixmax) then  ! Set a sensible max to ensure there will not be an int*4 overflow
        write(*,'(a,es12.5,a)') 'PIDinitally:ERROR: ',
     &    npixmax*1.0e-6,' max no. megapixels exceeded.'
        stop
      endif

      ! Detector mode:
      write(*,'(a)') 'Detector mode:'
      read(*,*) detmode
      select case (detmode)
      case (1)
        write(*,'(a)') ' Energy integrating'
        isEinteg = .true.       ! Switches help faster execution in tallying routine
        write(*,'(a)') 'E threshold ignored, not used in this mode.'
        read(*,'(a80)') buffer
        write(*,'(a)') 'E bins and no. bins ignored, '//
     &    'not used in this mode.'
        read(*,'(a80)') buffer
        nebin = 1
      case (2)
        write(*,'(a)') ' Photon counting'
        isEinteg = .false.
        isPCM = .true.
        write(*,'(a)') 'Energy deposition threshold (eV):'
        read(*,*) ethr
        write(*,'(1x,es12.5)') ethr
        write(*,'(a)') 'E bins and no. bins ignored, '//
     &    'not used in this mode.'
        read(*,'(a80)') buffer
        nebin = 1
      case (3)
        write(*,'(a)') ' Photon energy discriminating'
        isEinteg = .false.
        isPCM = .false.
        write(*,'(a)') 'E threshold ignored, not used in this mode.'
        read(*,'(a80)') buffer
        write(*,'(a)')
     &    'Spectrum limits [Emin,Emax] (eV) and no. bins:'
        read(*,*) emin,emax,nebin
        write(*,'(2(1x,es12.5),1x,i0)') emin,emax,nebin
        if (emin.lt.0.0.or.emin.ge.emax.or.nebin.lt.1) then
          write(*,'(a)') 'PIDinitally:ERROR: invalid values.'
          stop
        endif
        ebin = (emax-emin)/nebin
        iebin = 1.0/ebin
        nebin = nebin+1  ! Add one bin to store separately counts with E=Emax
      case default
        write(*,'(1x,i0)') detmode
        write(*,'(a)') 'PIDinitally:ERROR: invalid mode.'
        stop
      end select
      if (real(nx*ny)*nebin.gt.npixmax) then  ! Ensure there will not be an int*4 overflow
        write(*,'(a,es12.5,a)') 'PIDinitally:ERROR: ',
     &    npixmax*1.0e-6,' max no. (Mpixels x E bins) exceeded.'
        stop
      endif

      write(*,'(a)') 'Energy resolution parameters A(eV^2),B(eV):'
      read(*,*) c0,c1
      write(*,'(2(1x,es12.5))') c0,c1
      if (min(c0,c1).lt.0.0) then
        write(*,*) 'PIDinitally:ERROR: parameters cannot be negative.'
        stop
      endif
      isEspread = .false.
      if (max(c0,c1).gt.0.0) isEspread = .true.
      c0 = c0*fwhm2sig2  ! This is in units of sigma^2
      c1 = c1*fwhm2sig2

!! >>>> CROSS TALK CURRENTLY NOT AVAILABLE >>>>
!! >>>> Input data in PID section:
!! 0                               SPATIAL SPREAD (aka CROSSTALK) MODEL
!! 1.0                             PARAMETER OF CROSSTALK (MODEL=1,2,3) OR FILE NAME (MODEL=4)
!!   possible PoSF models (r^2=x^2+y^2):
!!     1 UNIFORM: p(x,y) ~ HeavisideStep(R^2-r^2); param: R(cm)
!!     2 GAUSSIAN: p(x,y) ~ exp(-r^2/(2*sigma^2)); param: sigma(cm)
!!     3 RADIAL ATT: p(x,y) ~ exp(-mu*r)/r; param: mu(cm-1)
!!     4 FROM FILE
!! >>>> F code for input:
!!
!!      ! Crosstalk model
!!      write(*,'(a)') 'Crosstalk model:'
!!      read(*,*) crosstalk
!!      select case (crosstalk)
!!        case (0)
!!          write(*,'(a)') ' No spread'
!!        case (1)
!!          write(*,'(a)') ' Uniform up to a radius R: '//
!!     &      'p(x,y) ~ HeavisideStep(R^2-r^2); param: R'
!!        case (2)
!!          write(*,'(a)') ' Gaussian: '//
!!     &      'p(x,y) ~ exp(-r^2/(2*sigma^2)); param: sigma'
!!        case (3)
!!          write(*,'(a)') ' Radial attenuation: '//
!!     &      'p(x,y) ~ exp(-mu*r)/r; param: mu'
!!        case (4)
!!          write(*,'(a)') ' Crosstalk matrix read from file:'
!!        case default
!!          write(*,'(a)') 'PIDinitally:ERROR: invalid model.'
!!          stop
!!      end select
!!
!! >>>> Draft code for crosstalk models:
!!      if (crosstalk.ne.4) then
!!        write(*,'(a)') 'Crosstalk model parameter:'
!!        read(*,*) crosParam
!!        if (crosstalk.eq.0) then
!!          write(*,'(a)') ' Irrelevant, not used in current model.'
!!        else
!!          write(*,'(1x,es12.5)') crosParam
!!          if (crosParam.le.0.0) then
!!            write(*,'(a)') 'PIDinitally:ERROR: invalid parameter.'
!!            stop
!!          endif
!!        endif
!!      else
!!        read(*,'(a80)') buffer               ! Read filename
!!        buffer = adjustl(buffer)             ! Remove leading blanks
!!        buffer = buffer(1:scan(buffer,' '))  ! Clip at 1st blank
!!        write(*,'(1x,a)') buffer
!!        ufile = finduf()
!!        open(ufile,file=buffer,status='old',iostat=errorcode)
!!        if (errorcode.ne.0) then
!!          write(*,'(a)') 'PIDinitally:ERROR: unable to open file.'
!!          stop
!!        endif
!!        read(ufile,*) crosstalkData   To be determined
!!        close(ufile)
!!        write(*,'(a)') 'Crosstalk matrix:'
!!        write(*,*) crosstalkData
!!        if (.not.goodData) then
!!          write(*,'(a)') 'PIDinitally:ERROR: Invalid data.'
!!          stop
!!        endif
!!      endif

      ! Format
      write(*,'(a)') 'Report format:'
      read(*,*) repformat
      select case (repformat)
      case (1)
        write(*,'(a)') ' Columnar'
      case (2)
        write(*,'(a)') ' Matrix'
        if (detmode.eq.3) then
          write(*,'(a)') 'PIDinitally:ERROR: '//
     &      'Matrix format is incompatible with '//
     &      'energy discriminating mode.'
          stop
        endif
      case (3)
        write(*,'(a)') ' Binary'
      case default
        write(*,'(1x,i0)') repformat
        write(*,'(a)') 'PIDinitally:ERROR: invalid format.'
        stop
      end select

      write(*,'(a)') 'Relative uncertainty (%) requested:'
      read(*,*) unclimit
      write(*,'(1x,es12.5)') unclimit

      read(*,'(a80)') buffer
      if (index(buffer,eos).eq.0) then
        write(*,'(a)')
     &    'PIDinitally:ERROR: End-Of-Section mark not found'
        write(*,'(a,a)') '  expecting to find: ',eos
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif

      ! Allocate scoring arrays:
      ijmax = nx*ny                          ! Size of 2D arrays
      select case (detmode)
      case (1)                               ! Energy integrating
        ijkmax = ijmax                       ! Note: DUMP uses this var even for 2D arrays
        ijwght = 1                           ! WGHTLAST defined for simplicity, unused for this mode
      case (2)                               ! Photon counting
        ijkmax = ijmax                       ! All arrays are 2D
        ijwght = ijmax                       ! WGHTLAST used
      case (3)                               ! Photon energy discriminating
        ijkmax = ijmax*nebin                 ! Size of 3D arrays
        ijwght = ijmax
      end select
      ! Memory usage:
      mem = sizeofReal8*real(ijkmax*2.0+ijmax*2.0+ijwght)
      write(*,'(a)') 'Memory used to allocate arrays (MB):'
      write(*,'(1x,es12.5)') mem*1.0e-6
      allocate(image(ijkmax),image2(ijkmax),edeptmp(ijmax),
     &         nlast(ijmax),wghtlast(ijwght),stat=errorcode)
      if (errorcode.ne.0) then
        write(*,'(a)') 'PIDinitally:ERROR: not enough memory.'
        stop
      endif
      memPID = mem  ! Store for later use by report

      ! Init arrays:
      image    = 0.0
      image2   = 0.0
      edeptmp  = 0.0
      nlast    = 0.0
      wghtlast = 0.0

      ! Sim restart file:
      call PIDdump(0)

      write(*,'(a)') '>>>> PID tally initialization finished >>>>'
      end


      subroutine setDetectFrame
!*******************************************************************
!*    Identifies detector reference frame and prepares rot matrix. *
!*******************************************************************
      use PIDmod
      use PENGEOMcommonsMod
      use PENGEOM_mod, only : nxg,mater,nbody
      use dataTypesMod
      implicit none

      integer*4 i,isurf,sidepoint,surf,surfmax
      real*8 normx(3),scalarx,normy(3),scalary,normtop(3),scalartop
      real*8 imodu,detbox,pixgrid
      real*8, parameter :: eps=1.0d-13,oneplus=1.0d0+eps

      ! Search last body containing the detector material
      bodydet = 0                     ! Init void
      do i=nbody,1,-1                 ! Reverse search, start from top
        if (mater(i).eq.matdet) then  ! This body is detector
          write(*,'(a)')
     &      'Detection material found in body/module no. '//
     &      '(PENGEOM internal numbering): '
          write(*,'(1x,i0)') i
          bodydet = i
          exit                        ! Only 1st ocurrence matters
        endif
      enddo
      if (bodydet.eq.0) then
        write(*,'(a)') 'setDetectFrame:ERROR: '//
     &    'detection material not found.'
        stop
      endif
      surfmax = ksurf(bodydet,nxg)    ! Num. surfaces limiting this body

      ! Search for detector box faces:
      write(*,'(a)') 'Searching for detector box faces...'
      isurf = 0                     ! Init
      do i=1,surfmax                ! For all surfaces limiting BODYDET
        sidepoint = 2*kflag(bodydet,i)-3  ! Side pointer -1,+1,+3,etc
        if (sidepoint.gt.+1) cycle  ! Skip surfaces limiting inserted bodies
        isurf = i                   ! Found x-face
        exit
      enddo
      if (isurf.eq.0) then
        write(*,'(a)') 'setDetectFrame:ERROR: '//
     &    'detector box faces not defined.'
        stop
      endif
      surf = ksurf(bodydet,isurf)   ! Surface number
      write(*,'(a)') 'Found x-normal face; it is surface no. '//
     &  '(PENGEOM internal numbering):'
      write(*,'(1x,i0)') surf
      if (kplane(surf).ne.1) then
        write(*,'(a)') 'setDetectFrame:ERROR: '//
     &      'x-normal face is not a plane.'
        stop
      endif
      normx = (/ax(surf),ay(surf),az(surf)/)  ! Surface normal
      imodu = sidepoint/sqrt(sum(normx**2))   ! Signed inverse modulus
      normx = normx*imodu                     ! Inwards-oriented unit normal
      scalarx = a0(surf)*imodu                ! Scalar term, normalized
      write(*,'(a)') 'Searching for the other 3 faces...'

      ! Seach for x-normal top face:
      isurf = isurf+1                         ! Next surface ID should be x-normal top face
      if (isurf.gt.surfmax) then
        write(*,'(a)') 'setDetectFrame:ERROR: '//
     &    'x-normal top face not defined.'
        stop
      endif
      sidepoint = 2*kflag(bodydet,isurf)-3
      if (sidepoint.gt.+1) then               ! This surface is not a valid face, it limits an inserted body
        write(*,'(a)') 'setDetectFrame:ERROR: '//
     &    'x-normal top face not defined.'
        stop
      endif
      surf = ksurf(bodydet,isurf)             ! Surface number
      if (kplane(surf).ne.1) then
        write(*,'(a)') 'setDetectFrame:ERROR: '//
     &      'x-normal top face is not a plane.'
        stop
      endif
      normtop = (/ax(surf),ay(surf),az(surf)/)
      imodu = -sidepoint/sqrt(sum(normtop**2))
      normtop = normtop*imodu                 ! Outwards-oriented unit normal, should coincide with NORMX
      scalartop = a0(surf)*imodu              ! Scalar term, normalized

      ! Check if parallel and correctly oriented:
      if (abs(dot_product(normx,normtop)-1.0d0).gt.eps) then
        write(*,'(a)') 'setDetectFrame:ERROR: '//
     &    'x-normal and x-normal top faces are not parallel or '//
     &    'their side pointers are incorrectly defined.'
        stop
      endif

      ! Set grid size along x:
      detbox = abs(scalarx-scalartop)        ! Detector box size = Distance between planes
      write(*,'(a)') 'Detector box x-size is (cm):'
      write(*,'(1x,es12.5)') detbox
      if (detbox.lt.eps) then
        write(*,'(a)') 'setDetectFrame:ERROR: '//
     &    'detector box size is too small.'
        stop
      endif
      pixgrid = nx*dx                         ! Pixels grid size, as defined by user
      if (pixgrid.ne.0.0) then                ! Pixels grid size explicitly given by user
        if (abs(pixgrid-detbox).gt.eps*detbox) then
          write(*,'(a)') '(pixels grid and detector box sizes '//
     &      'do not coincide; detector box x-size is ignored.)'
        endif
      else if (nx.eq.0) then                  ! Infer num. pixels from geometry
        write(*,'(a)') 'Number of pixels along x set to:'
        if (detbox*oneplus/dx.gt.npixmax) then
          write(*,'(a)') 'setDetectFrame:ERROR: '//
     &    'overflow in no. pixels; redefine detector box parameters.'
          stop
        endif
        nx = detbox*oneplus/dx
        write(*,'(1x,i0)') nx
      else                                    ! user entered dx=0 => infer pixel size
        dx = detbox/nx
        write(*,'(a)') 'Pixel size along x set to (cm):'
        write(*,'(1x,es12.5)') dx
      endif
      idx = 1.0/dx

      ! Seach for y-normal face:
      isurf = isurf+1                         ! Next surface ID should be y-normal face
      if (isurf.gt.surfmax) then
        write(*,'(a)') 'setDetectFrame:ERROR: '//
     &    'y-normal face not defined.'
        stop
      endif
      sidepoint = 2*kflag(bodydet,isurf)-3
      if (sidepoint.gt.+1) then               ! This surface is not a valid face, it limits an inserted body
        write(*,'(a)') 'setDetectFrame:ERROR: '//
     &    'y-normal face not defined.'
        stop
      endif
      surf = ksurf(bodydet,isurf)             ! Surface number
      if (kplane(surf).ne.1) then
        write(*,'(a)') 'setDetectFrame:ERROR: '//
     &      'y-normal face is not a plane.'
        stop
      endif
      normy = (/ax(surf),ay(surf),az(surf)/)
      imodu = sidepoint/sqrt(sum(normy**2))
      normy = normy*imodu                     ! Inwards-oriented unit normal
      scalary = a0(surf)*imodu                ! Scalar term, normalized

      ! Seach for y-normal top face:
      isurf = isurf+1                         ! Next surface ID should be y-normal top face
      if (isurf.gt.surfmax) then
        write(*,'(a)') 'setDetectFrame:ERROR: '//
     &    'y-normal top face not defined.'
        stop
      endif
      sidepoint = 2*kflag(bodydet,isurf)-3
      if (sidepoint.gt.+1) then               ! This surface is not a valid face, it limits an inserted body
        write(*,'(a)') 'setDetectFrame:ERROR: '//
     &    'y-normal top face not defined.'
        stop
      endif
      surf = ksurf(bodydet,isurf)             ! Surface number
      if (kplane(surf).ne.1) then
        write(*,'(a)') 'setDetectFrame:ERROR: '//
     &      'y-normal top face is not a plane.'
        stop
      endif
      normtop = (/ax(surf),ay(surf),az(surf)/)
      imodu = -sidepoint/sqrt(sum(normtop**2))
      normtop = normtop*imodu                 ! Outwards-oriented unit normal, should coincide with NORMY
      scalartop = a0(surf)*imodu              ! Scalar term, normalized

      ! Check if parallel and correctly oriented:
      if (abs(dot_product(normy,normtop)-1.0d0).gt.eps) then
        write(*,'(a)') 'setDetectFrame:ERROR: '//
     &    'y-normal and y-normal top faces are not parallel or '//
     &    'their side pointers are incorrectly defined.'
        stop
      endif

      ! Set grid size along y:
      detbox = abs(scalary-scalartop)         ! Detector box size = Distance between planes
      write(*,'(a)') 'Detector box y-size is (cm):'
      write(*,'(1x,es12.5)') detbox
      if (detbox.lt.eps) then
        write(*,'(a)') 'setDetectFrame:ERROR: '//
     &    'detector box size is too small.'
        stop
      endif
      pixgrid = ny*dy
      if (pixgrid.ne.0.0) then                ! Pixels grid size given by user
        if (abs(pixgrid-detbox).gt.eps*detbox) then
          write(*,'(a)') '(pixels grid and detector box sizes '//
     &      'do not coincide; detector box y-size is ignored.)'
        endif
      else if (ny.eq.0) then                  ! Infer num. pixels
        write(*,'(a)') 'Number of pixels along y set to:'
        if (detbox*oneplus/dy.gt.npixmax) then
          write(*,'(a)') 'setDetectFrame:ERROR: '//
     &    'overflow in no. pixels; redefine detector box parameters.'
          stop
        endif
        ny = detbox*oneplus/dy
        write(*,'(1x,i0)') ny
      else                                    ! Infer pixel size
        dy = detbox/ny
        write(*,'(a)') 'Pixels size along y set to (cm):'
        write(*,'(1x,es12.5)') dy
      endif
      idy = 1.0/dy

      ! Check perpendicularity of x&y faces:
      if (abs(dot_product(normx,normy)).gt.eps) then
        write(*,'(a)') 'setDetectFrame:ERROR: '//
     &    'x & y faces are not perpendicular.'
        stop
      endif

      ! Lab->Pixels change of bases (only x&y coords will be needed):
      lab2det(1,:) = normx
      lab2det(2,:) = normy

      isLab = .false.
      if (abs(normx(1)-1.0d0).lt.eps.and.
     &    abs(normy(2)-1.0d0).lt.eps) then
        isLab = .true.
        write(*,'(a)') 'Detector reference frame coincides with lab '//
     &    'frame; no rotation applied for scoring.'
      endif

      ! Find an arbitrary point r on the line L defined by the two orthogonal planes
      !   defining the detector box. Point r is selected such that its distance to the
      !   origin in minimal, that is, it lies on the plane perpendicular to L
      !   that passes through the origin. Thus, it satisfies the eqs.
      !      normx.r+scalarx  = 0  (x-face plane)
      !      normy.r+scalary  = 0  (y-face plane)
      !     (normx X normy).r = 0  (orthogonal plane through origin)
      !
      !   ('.' represents scalar product, 'X' represents cross product)
      !   The solution of the above system for r is:
      r0 = -scalarx*normx-scalary*normy
      end


      real*8 function espread(x)
!*******************************************************************
!*    Applies a Gaussian dispersion to x.                          *
!*                                                                 *
!*    Input:                                                       *
!*      x: value to be 'dispersed'.                                *
!*         Gaussian parameters taken from global module.           *
!*******************************************************************
      use PIDmod
      implicit none
      real*8 x,gausspdf

      espread = x+sqrt(c0+c1*x)*gausspdf()  ! The squared std deviation depends linearly on x
      end function


      subroutine PIDdump(mode)
!*******************************************************************
!*    Dumps into or reads data from a dump file.                   *
!*                                                                 *
!*    Input:                                                       *
!*      mode: 1 to write dump file, else to read from it.          *
!*******************************************************************
      use PIDmod
      use dumpMod
      implicit none
      integer mode

      if (mode.eq.1) then
        if (dumpfile.lt.0) return  ! No dump file open
        write(dumpfile) image(1:ijkmax),image2(1:ijkmax)
      else
        if (restartfile.lt.0) return  ! No restart file open
        read(restartfile) image(1:ijkmax),image2(1:ijkmax)
      endif
      end


!>>>> End Of File >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
