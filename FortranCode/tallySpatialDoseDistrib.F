!*******************************************************************
!*                          TALLY                                  *
!*             SPATIAL (3D) DOSE DISTRIBUTION                      *
!*                                                                 *
!* Short description:                                              *
!*   Tallies the dose distribution in a grid of 3D bins            *
!*   superimposed on the solid-body PENGEOM geometry.              *
!* Rev:                                                            *
!*   2014-12-30                                                    *
!*******************************************************************


!*******************************************************************
!*******************************************************************
!*    Modules section                                              *
!*******************************************************************
!*******************************************************************

      module SDDmod
!*******************************************************************
!*    Vars for the SDD tally.                                      *
!*******************************************************************
      implicit none
      save
      logical active
      integer prtxyz
      integer*4 :: nx,ny,nz,nxny,nbin
      real*8 xmin,ymin,zmin,dx,dy,dz,idx,idy,idz,unclimit
      real*8, parameter :: nbinmax=1.024e9 ! Max num bins (do not exceed int*4)
      real*8 :: memSDD=0.0                 ! Memory used by arrays (bytes)
      real*8, allocatable ::
     &   edptmp(:),edep(:),edep2(:),idens(:),nlast(:)
      end module


!*******************************************************************
!*******************************************************************
!*    Routines section                                             *
!*******************************************************************
!*******************************************************************

      subroutine SDDtally(mode,arg)
!*******************************************************************
!*    Input:                                                       *
!*      mode: identifies the state of the calling routine          *
!*      arg: energy loss (mode<0) or history no. (mode=1)          *
!*******************************************************************
      use SDDmod
      use TRACK_mod, only : x,y,z,wght
      use ctrsimMod, only : nhist  ! Current history number
      implicit none
      integer mode
      real*8 arg

      integer*4 bin,i,j,k
      real*8, parameter :: oneplus=1.0000000000001d0

      if (.not.active) return

      if (mode.le.0) then
        if (arg.eq.0.0) return         ! Nothing to deposit

        ! Check if particle is inside tally region:
        i = (x-xmin)*idx+oneplus       ! ONEPLUS ensures proper truncation when IDX=0
        if (i.lt.1.or.i.gt.nx) return
        j = (y-ymin)*idy+oneplus
        if (j.lt.1.or.j.gt.ny) return
        k = (z-zmin)*idz+oneplus
        if (k.lt.1.or.k.gt.nz) return
        bin = i+(j-1)*nx+(k-1)*nxny    ! Map i,j,k into a single index

        ! Transfer partial tally to totals only when a new history visits:
        if (nhist.gt.nlast(bin)) then  ! Visit of a new history
          edep(bin)  = edep(bin) +edptmp(bin)
          edep2(bin) = edep2(bin)+edptmp(bin)**2
          edptmp(bin)= arg*wght
          nlast(bin) = nhist+0.5       ! +0.5 to avoid round-off problems
        else
          edptmp(bin) = edptmp(bin)+arg*wght
        endif

      endif
      end


      subroutine SDDreport(mode,n,cputim,uncdone)
!*******************************************************************
!*    Input:                                                       *
!*      mode: -1 if called at end-of-simulation, 0 for dump only,  *
!*             1 for tally report only, 2 for both.                *
!*      n: no. of histories simulated                              *
!*      cputim -> elapsed CPU time                                 *
!*    Output:                                                      *
!*      uncdone -> 2 if uncert reached, 1 if not defined, 0 else   *
!*                 not computed if mode=0                          *
!*******************************************************************
      use SDDmod
      use RSEEDcommonMod
      implicit none
      integer mode,uncdone
      real*8 n,cputim

      logical printx,printy,printz
      character*80 buffer
      integer out,finduf,errorcode,ninteg
      integer*4 bin,i,j,k,kbin,jkbin,nxef,nyef,nzef,nchan
      real*8 q,sigma,eff,avesig,maxq,fact,factint,x,y,z,uncert
      real*8 xmiddle,ymiddle,zmiddle,invn

      uncdone = 1
      if (.not.active) return

      ! Dump temp counters:
      do bin=1,nbin
        if (nlast(bin).lt.0.5) cycle           ! Skip empty bins
        edep(bin)  = edep(bin) +edptmp(bin)    ! Transfer temp counter
        edep2(bin) = edep2(bin)+edptmp(bin)**2
        edptmp(bin)= 0.0                       ! Reset counter
        nlast(bin) = 0.0                       ! Reset last visited to avoid recounting in next report
      enddo

      if (mode.ne.1) call SDDdump(1) ! Sim dump file, notice that edptmp and nlast have been cleared before
      if (mode.eq.0) return          ! Dump only, no report

      ! Prepare x,y,z factors:
      ninteg = 0                     ! No. of dimensions to integrate
      factint = 1.0                  ! Factor of integration
      printx = .false.
      printy = .false.
      printz = .false.
      nxef = nx
      if (idx.eq.0.0d0) then         ! Integrating over x
        nxef = 0
        ninteg = ninteg+1
      else
        factint = factint*idx
        if (prtxyz.eq.1) printx = .true.
      endif
      nyef = ny
      if (idy.eq.0.0d0) then         ! Integrating over y
        nyef = 0
        ninteg = ninteg+1
      else
        factint = factint*idy
        if (prtxyz.eq.1) printy = .true.
      endif
      nzef = nz
      if (idz.eq.0.0d0) then         ! Integrating over z
        nzef = 0
        ninteg = ninteg+1
      else
        factint = factint*idz
        if (prtxyz.eq.1) printz = .true.
      endif

      ! Prepare output file:
      select case (ninteg)
      case (0)
        buffer = 'tallySpatialDoseDistrib-3D.dat'
        if (prtxyz.eq.-1)
     &    buffer = 'tallySpatialDoseDistrib-3D-binary.dat'
      case (1)
        buffer = 'tallySpatialDoseDistrib-2D.dat'
        if (prtxyz.eq.-1)
     &    buffer = 'tallySpatialDoseDistrib-2D-binary.dat'
      case (2)
        buffer = 'tallySpatialDoseDistrib-1D.dat'
        if (prtxyz.eq.-1)
     &    buffer = 'tallySpatialDoseDistrib-1D-binary.dat'
      case default
        write(*,'(a)') 'SDDreport:ERROR: internal error.'
        return
      end select

      out = finduf()
      if (prtxyz.eq.-1) then                    ! Open data file in binary mode
        open(out,file=buffer,status='replace',access='stream',
     &       iostat=errorcode)
      else                                      ! Open data file in text mode
        open(out,file=buffer,iostat=errorcode)
      endif
      if (errorcode.ne.0) then
        write(*,*) ''
        write(*,'(a)')
     &    '*********************************************'
        write(*,'(a)')
     &    'SDDreport:ERROR: cannot open output data file'
        write(*,'(a)')
     &    '*********************************************'
        close(out)                              ! Just in case
        return                                  ! Forget report and continue simulation
      endif

      if (prtxyz.ne.-1) then ! Write header only if not in binary mode
        write(out,'(a)')
     &  '#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
        write(out,'(a)') '# [SECTION REPORT SPATIAL DOSE DISTRIB]'
        select case (ninteg)
        case (0)
          write(out,'(a)') '# Dose units are: eV/g per history'
        case (1)
          write(out,'(a)') '# Dose units are: eV.cm/g per history'
        case (2)
          write(out,'(a)') '# Dose units are: eV.cm^2/g per history'
        end select
        write(out,'(a)') '#'
        write(out,'(a)') '# No. of bins in x,y,z directions and total:'
        write(out,'(a,4(1x,i0))') '# ',nxef,nyef,nzef,nbin
        write(out,'(a)') '# Min values for x,y,z(cm):'
        write(out,'(a,3(1x,es12.5))') '# ',xmin,ymin,zmin
        write(out,'(a)') '# Bin widths for x,y,z(cm):'
        write(out,'(a,3(1x,es12.5))') '# ',dx,dy,dz
        write(out,'(a)') '#'
        write(out,'(a)')
     &    '# For plotting purposes, two values per bin coordinate '//
     &    'are given, namely, the low end and the middle point '//
     &    'of each bin.'
        write(out,'(a)') '#'
        write(out,'(a,$)') '# '
        if (printx)
     &      write(out,'(a,$)') 'xBinIndex : xLow(cm) : xMiddle(cm) : '
        if (printy)
     &      write(out,'(a,$)') 'yBinIndex : yLow(cm) : yMiddle(cm) : '
        if (printz)
     &      write(out,'(a,$)') 'zBinIndex : zLow(cm) : zMiddle(cm) : '
        write(out,'(a)') 'dose : +-2sigma'
      endif

      ! Write data:
      invn   = 1.0/n
      nchan  = 0
      avesig = 0.0
      maxq   = 0.5*maxval(edep)  ! 1/2 of the max score

      do k=1,nz
        z = zmin+dz*(k-1)
        zmiddle = z+dz*0.5
        kbin = (k-1)*nxny                                  ! Used below to map i,j,k into a single index
        if (idz.gt.0.0.and..not.printz.and.prtxyz.ne.-1)   ! Since z is not written, give at least a summary
     &    write(out,'(a,i0,a,es12.5)')
     &    '# zBinIndex=',k,' zMiddle(cm)=',zmiddle

        do j=1,ny
          y = ymin+dy*(j-1)
          ymiddle = y+dy*0.5
          jkbin = (j-1)*nx+kbin                            ! Used below to map i,j,k into a single index
          if (idy.gt.0.0.and..not.printy.and.prtxyz.ne.-1) ! Since y is not written, give at least a summary
     &      write(out,'(a,i0,a,es12.5)')
     &      '# yBinIndex=',j,' yMiddle(cm)=',ymiddle

          do i=1,nx
            x = xmin+dx*(i-1)
            xmiddle = x+dx*0.5
            if (printx)
     &        write(out,'(1x,i5,2(1x,es12.5),$)') i,x,xmiddle
            if (printy)
     &        write(out,'(1x,i5,2(1x,es12.5),$)') j,y,ymiddle
            if (printz)
     &        write(out,'(1x,i5,2(1x,es12.5),$)') k,z,zmiddle
            bin = i+jkbin                               ! Map i,j,k into a single index
            fact = factint*idens(bin)
            q = edep(bin)*invn
            sigma = sqrt(max((edep2(bin)*invn-q**2)*invn,0.0))*fact
            q = q*fact
            if (prtxyz.eq.-1) then                      ! Binary mode
              write(out) q,2.0*sigma
            else                                        ! Text mode
              write(out,'(1x,es12.5,1x,es7.1)') q,2.0*sigma
            endif
            ! Evaluate average uncertainty for scores above 1/2 max score:
            if (edep(bin).gt.maxq) then
              avesig = avesig+(sigma/q)**2
              nchan = nchan+1
            endif
          enddo

          if (prtxyz.eq.-1) cycle        ! No separators in binary mode
          if (nx.gt.1) write(out,*) ' '  ! Separate 2D data blocks
        enddo

        if (prtxyz.eq.-1) cycle          ! No separators in binary mode
        if (ny.gt.1) write(out,*) ' '    ! Separate 3D data blocks
      enddo

      uncdone = 0
      if (nchan.gt.0) then
        uncert = 200.0*sqrt(avesig/nchan)
        if (uncert.lt.unclimit) uncdone = 2  ! Uncertainty reached
      else
        uncert = 0.0                     ! Uncertainty assumed not reached when score is nil
      endif

      if (prtxyz.ne.-1) then ! Generic report if not in binary mode
        ! Generic report:
        write(out,'(a)') ' '
        write(out,'(a)') '# Performance report'
        write(out,'(a)') '#   Memory used (MB):'
        write(out,'(a,es12.5)') '#   ',memSDD*1.0e-6
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
     &   '#   Average uncertainty (above 1/2 max score) in % [uncert]:'
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


      subroutine SDDinitally(mem)
!*******************************************************************
!*    Initializes. To be called before TALLY                       *
!*                                                                 *
!*    Output:                                                      *
!*      mem: Memory used by tally arrays.                          *
!*******************************************************************
      use SDDmod
      use dataTypesMod
      use TRACK_mod, only : x,y,z,u,v,w
      implicit none
      real*8 mem

      character*(*) secid,eos
      parameter (secid=
     &  '[SECTION TALLY SPATIAL DOSE DISTRIB v.2009-06-15]')
      parameter (eos='[END OF SDD SECTION]')
      character*80 buffer,upcase
      integer errorcode
      integer*4 i,j,k,bin,kbin,jkbin
      real*8 xmax,ymax,zmax,rnbin,localdens,flocaldens

      write(*,*) ''
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
      buffer = adjustl(upcase(buffer))
      buffer = buffer(1:scan(buffer,' ')) ! Clip at 1st blank
      if (buffer.eq.'ON') then
        active = .true.
      else if (buffer.eq.'OFF') then
        active = .false.
        write(*, '(a)')
     &    '>>>> Tally Spatial Dose Distrib is OFF >>>>'
        do
          read(*,'(a80)',iostat=errorcode) buffer
          if (errorcode.ne.0) then
            write(*,'(a,a,a)') 'SDDinitally:ERROR: ',
     &       'Unable to find End-Of-Section mark: ',eos
            stop
          endif
          if (index(buffer,eos).ne.0) return
        enddo
      else
        write(*,'(a)')
     &    'SDDinitally:ERROR: expecting to find ON or OFF'
        write(*,'(a)') 'found instead:'
        write(*,'(a)') buffer
        stop
      endif

      ! Read bin parameters sizes:
      write(*,'(a)') 'xmin,xmax,nx:'
      read(*,*) xmin,xmax,nx
      write(*,'(2(1x,es12.5),1x,i0)') xmin,xmax,nx
      if (xmin.gt.xmax.or.nx.lt.0.or.(nx.gt.0.and.xmin.eq.xmax)) then
        write(*,'(a)') 'SDDinitally:ERROR: invalid entry.'
        stop
      endif
      if (nx.eq.0) then  ! Redefine parameters when dose is integrated along x-axis
        dx = xmax-xmin
        idx = 0.0d0  ! Keep it double precision
        nx = 1
      else
        dx = (xmax-xmin)/nx
        idx = 1.0/dx
      endif
      write(*,'(a)') 'ymin,ymax,ny:'
      read(*,*) ymin,ymax,ny
      write(*,'(2(1x,es12.5),1x,i0)') ymin,ymax,ny
      if (ymin.gt.ymax.or.ny.lt.0.or.(ny.gt.0.and.ymin.eq.ymax)) then
        write(*,*) 'SDDinitally:ERROR: invalid entry.'
        stop
      endif
      if (ny.eq.0) then
        dy = ymax-ymin
        idy = 0.0d0  ! Keep it double precision
        ny = 1
      else
        dy = (ymax-ymin)/ny
        idy = 1.0/dy
      endif
      write(*,'(a)') 'zmin,zmax,nz:'
      read(*,*) zmin,zmax,nz
      write(*,'(2(1x,es12.5),1x,i0)') zmin,zmax,nz
      if (zmin.gt.zmax.or.nz.lt.0.or.(nz.gt.0.and.zmin.eq.zmax)) then
        write(*,*) 'SDDinitally:ERROR: invalid entry.'
        stop
      endif
      if (nz.eq.0) then
        dz = zmax-zmin
        idz = 0.0d0  ! Keep it double precision
        nz = 1
      else
        dz = (zmax-zmin)/nz
        idz = 1.0/dz
      endif

      if (max(idx,idy,idz).eq.0.0d0) then
        write(*,*) 'SDDinitally:ERROR: Zero bins defined.'
        stop
      endif
      rnbin = dble(nx)*dble(ny)*dble(nz)
      write(*,'(a)') 'Total number of megabins:'
      write(*,'(1x,es12.5)') rnbin*1.0e-6
      if (rnbin.gt.nbinmax) then    ! Set a sensible max to ensure there will not be an int*4 overflow
        write(*,'(a,es12.5,a)') 'SDDinitally:ERROR: ',
     &    nbinmax*1e-6,' max no. of megabins exceeded.'
        stop
      endif
      nxny = nx*ny
      nbin = nx*ny*nz

      ! Coordinates:
      write(*,'(a)') 'Print coordinates:'
      read(*,*) prtxyz
      select case (prtxyz)
        case (1)
          write(*,'(a)') ' yes'
        case (0)
          write(*,'(a)') ' no'
        case (-1)
          write(*,'(a)') ' no; printing data in binary format'
        case default
          write(*,'(a)') 'SDDinitally:ERROR: expecting 1,0 or -1.'
          stop
      end select

      ! Uncert:
      write(*,'(a)') 'Relative uncertainty (%) requested:'
      read(*,*) unclimit
      write(*,'(1x,es12.5)') unclimit

      ! Allocate arrays:
      mem = sizeOfReal8*5*real(nbin)
      write(*,'(a)') 'Memory used to allocate arrays (MB):'
      write(*,'(1x,es12.5)') mem*1.0e-6
      allocate(edptmp(nbin),edep(nbin),edep2(nbin),idens(nbin),
     &         nlast(nbin),stat=errorcode)
      if (errorcode.ne.0) then
        write(*,'(a)') 'SDDinitally:ERROR: not enough memory.'
        stop
      endif
      memSDD = mem  ! Store for later use by report

      ! Init arrays:
      edptmp = 0.0
      edep   = 0.0
      edep2  = 0.0
      nlast  = 0.0d0
      idens  = 0.0
      u = 0.0  ! For locate() to operate properly
      v = 0.0
      w = 1.0

      do k=1,nz
        z = zmin+dz*(k-0.5)
        kbin = (k-1)*nxny             ! Used below to map i,j,k into a single index
        do j=1,ny
          y = ymin+dy*(j-0.5)
          jkbin = (j-1)*nx+kbin       ! Used below to map i,j,k into a single index
          do i=1,nx
            x = xmin+dx*(i-0.5)       ! Position at the middle of the bin
            bin = i+jkbin             ! Map i,j,k into a single index
            localdens = flocaldens()  ! Local mass density--quads or voxs
            if (localdens.gt.0.0) idens(bin) = 1.0/localdens
          enddo
        enddo
      enddo

      ! Check integrity:
      read(*,'(a80)') buffer
      if (index(buffer,eos).eq.0) then
        write(*,*) 'SDDinitally:ERROR: End-Of-Section mark not found'
        write(*,'(a,a)') '  expecting to find: ',eos
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif

      ! Sim restart file:
      call SDDdump(0)

      write(*,'(a)') '>>>> SDD tally initialization finished >>>>'
      end


      subroutine SDDdump(mode)
!*******************************************************************
!*    Dumps into or reads data from a dump file.                   *
!*                                                                 *
!*    Input:                                                       *
!*      mode: 1 to write dump file, else to read from it.          *
!*******************************************************************
      use SDDmod
      use dumpMod
      implicit none
      integer mode

      if (mode.eq.1) then
        if (dumpfile.lt.0) return  ! No dump file open
        write(dumpfile) edep(1:nbin),edep2(1:nbin)
      else
        if (restartfile.lt.0) return  ! No restart file open
        read(restartfile) edep(1:nbin),edep2(1:nbin)
      endif
      end


!>>>> End Of File >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
