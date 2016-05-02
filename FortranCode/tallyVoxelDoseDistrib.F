!*******************************************************************
!*               TALLY VOXEL DOSE DISTRIBUTION                     *
!*                                                                 *
!* Short description:                                              *
!*   Dose distribution in the volume elements (voxels) defined in  *
!*   a voxels file.                                                *
!*                                                                 *
!* Updates:                                                        *
!*     Modified by                                                 *
!*       J Sempau                                                  *
!*     Created by                                                  *
!*       A Badal & J Sempau, 2008                                  *
!*       Universitat Politecnica de Catalunya, Barcelona, Spain    *
!* Rev:                                                            *
!*   2014-12-30                                                    *
!*******************************************************************


!*******************************************************************
!*******************************************************************
!*    Modules section                                              *
!*******************************************************************
!*******************************************************************

      module VDDmod
!*******************************************************************
!*    Vars for VDD tally                                           *
!*******************************************************************
      implicit none
      save
      logical active
      integer prtxyz,prtdens
      integer*4 xvoxmin,xvoxmax,yvoxmin,yvoxmax,zvoxmin,zvoxmax
      real*8 :: unclimit,memVDD=0.0 ! Memory used by arrays (bytes)
      real*8, allocatable :: nlast(:),edptmp(:),edep(:),edep2(:)
      end


!*******************************************************************
!*******************************************************************
!*    Routines section                                             *
!*******************************************************************
!*******************************************************************

      subroutine VDDtally(mode,arg)
!*******************************************************************
!*    Input:                                                       *
!*      mode -> Identifies the state of the calling routine        *
!*      arg -> energy loss (mode<0) or history no. (mode=1)        *
!*******************************************************************
      use VDDmod
      use partvoxMod
      use TRACK_mod
      use ctrsimMod, only : nhist  ! Current history number
      implicit none
      integer mode
      real*8 arg

      if (.not.active) return                    ! Voxel doses not tallied

      if (mode.le.0) then                        ! Deposit energy
        if (absvox.eq.0.or.arg.eq.0.0d0) return  ! Inside quadrics or nothing to deposit

        if (nhist.gt.nlast(absvox)) then                  ! Visit of a new history
          edep(absvox)  = edep(absvox) +edptmp(absvox)    ! Transfer partial scoring to totals
          edep2(absvox) = edep2(absvox)+edptmp(absvox)**2 ! Score squared values for variance
          edptmp(absvox)= arg*wght                        ! And update temporary counter
          nlast(absvox) = nhist+0.5d0                     ! Update NLAST (+0.5 avoids prec. problems)
        else
          edptmp(absvox) = edptmp(absvox)+arg*wght        ! Same history: update temporary counter only
        endif

      endif
      end


      subroutine VDDreport(mode,n,cputim,uncdone)
!*******************************************************************
!*    Input:                                                       *
!*      mode -> -1 if called at end-of-simulation, 0 for dump only,*
!*               1 for tally report only, 2 for both.              *
!*      n -> no. of histories simulated                            *
!*      cputim -> elapsed CPU time                                 *
!*    Output:                                                      *
!*      uncdone -> 2 if uncert reached, 1 if not defined, 0 else   *
!*                 not computed if mode=0                          *
!*******************************************************************
      use VDDmod
      use geovoxMod
      use RSEEDcommonMod
      implicit none
      integer mode,uncdone
      real*8 n,cputim

      integer nchan,out,finduf,errorcode
      integer*4 vox,i,j,k
      real*8 q,sigma,eff,avesig,maxq,fact,x,y,z,uncert
      real*8 xmiddle,ymiddle,zmiddle,invn,fmassvox,thismassvox

      uncdone = 1
      if (.not.active) return

      ! Dump counters and obtain max score:
      do k=zvoxmin,zvoxmax
        do j=yvoxmin,yvoxmax
          do i=xvoxmin,xvoxmax
            vox = i+(j-1)*nx+(k-1)*nxy           ! Absolute voxel index
            if (nlast(vox).lt.0.5) cycle         ! Skip empty voxels
            edep(vox)  = edep(vox) +edptmp(vox)  ! Transfer temp counter
            edep2(vox) = edep2(vox)+edptmp(vox)**2
            edptmp(vox)= 0.0                     ! Reset counter
            nlast(vox) = 0.0                     ! Reset last visited to avoid recounting in next report
          enddo
        enddo
      enddo

      if (mode.ne.1) call VDDdump(1)  ! Sim dump file, notice that edptmp and nlast have been cleared before
      if (mode.eq.0) return           ! Dump only, no report

      ! Prepare output file:
      out = finduf()
      if (prtxyz.eq.-1) then ! Open data file in binary format
        open(out,file='tallyVoxelDoseDistrib-binary.dat',
     &       status='replace',access='stream',iostat=errorcode)
      else
        open(out,file='tallyVoxelDoseDistrib.dat',iostat=errorcode)
      endif
      if (errorcode.ne.0) then
        write(*,*)
        write(*,'(a)')
     &    '*********************************************'
        write(*,'(a)')
     &    'VDDreport:ERROR: cannot open output data file'
        write(*,'(a)')
     &    '*********************************************'
        close(out)  ! Just in case
        return
      endif

      if (prtxyz.ne.-1) then ! Write header only if not in binary mode
        write(out,'(a)')
     &  '#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
        write(out,'(a)') '# [SECTION REPORT VOXEL DOSE DISTRIB]'
        write(out,'(a)') '# Deposited energy over voxel mass'
        write(out,'(a)')
     &    '# (mass excludes voxel regions overlaid with quadrics)'
        write(out,'(a)') '#'
        write(out,'(a)') '# Dose units are eV/g per history'
        write(out,'(a)') '#'
        write(out,'(a)') '# Number of voxels in x,y,z:'
        write(out,'(a,3(1x,i0))') '# ',nx,ny,nz
        write(out,'(a)') '# Voxels size dx,dy,dz (cm):'
        write(out,'(a,3(1x,es12.5))')'# ',dx,dy,dz
        if (isQuad) then
          write(out,'(a)')
     &      '# The granularity used to evaluate overlaid voxels was:'
          write(out,'(a,1x,i0)') '# ',granul
        endif
        write(out,'(a)') '#'
        if (prtxyz.eq.1) then
          write(out,'(a)')
     &      '# For plotting purposes, two values per voxel '//
     &      'coordinate are given, namely, the low end and the '//
     &      'middle point of each voxel.'
          write(out,'(a)') '#'
        endif
        write(out,'(a,$)') '# '
        if (prtxyz.eq.1) then
          write(out,'(a,$)') 'xLow(cm) : xMiddle(cm) : '
          write(out,'(a,$)') 'yLow(cm) : yMiddle(cm) : '
          write(out,'(a,$)') 'zLow(cm) : zMiddle(cm) : '
        endif
        write(out,'(a,$)') 'dose(eV/g/hist) : +-2sigma'
        if (prtdens.eq.1) write(out,'(a,$)')
     &    ' : voxelMass(g) : Pure(+)/Overlaid(-)'
        write(out,'(a)') ''  ! End of line
      endif

      ! Write data:
      invn = 1.0/n
      avesig = 0.0
      nchan  = 0
      maxq = 0.5*maxval(edep)         ! 1/2 of the max score

      do k=zvoxmin,zvoxmax
        z = dz*(k-1)
        zmiddle = z+dz*0.5
        if (prtxyz.eq.0) then
          write(out,'(a,i0,a,es12.5)')
     &      '# zVoxIndex=',k,' zMiddle(cm)=',zmiddle
        else if (prtxyz.eq.1) then
          write(out,'(a,i0)') '# zVoxIndex=',k
        endif

        do j=yvoxmin,yvoxmax
          y = dy*(j-1)
          ymiddle = y+dy*0.5
          if (prtxyz.eq.0) then
            write(out,'(a,i0,a,es12.5)')
     &        '# yVoxIndex=',j,' yMiddle(cm)=',ymiddle
          else if (prtxyz.eq.1) then
            write(out,'(a,i0)') '# yVoxIndex=',j
          endif

          do i=xvoxmin,xvoxmax
            vox = i+(j-1)*nx+(k-1)*nxy    ! Absolute voxel index
            if (prtxyz.eq.1) then
              x = dx*(i-1)
              xmiddle = x+dx*0.5
              write(out,'(6(1x,es12.5),$)')
     &          x,xmiddle,y,ymiddle,z,zmiddle
            endif
            fact = 0.0
            thismassvox = fmassvox(vox)
            if (thismassvox.gt.0.0) fact = 1.0/thismassvox  ! Voxel mass may be null if partial vol
            q = edep(vox)*invn
            sigma = sqrt(max((edep2(vox)*invn-q**2)*invn,0.0))*fact
            q = q*fact
            if (prtxyz.eq.-1) then    ! Binary mode
              if (prtdens.ne.1) then  ! Do not print voxel mass
                write(out) q,2.0*sigma
              else
                write(out) q,2.0*sigma,thismassvox
              endif
            else                      ! Text mode
              if (prtdens.ne.1) then  ! Do not print voxel mass
                write(out,'(1x,es12.5,1x,es7.1)') q,2.0*sigma
              else if (matvox(vox).ge.0) then
                write(out,'(1x,es12.5,1x,es7.1,1x,es12.5,1x,a1)')
     &            q,2.0*sigma,thismassvox,'+'
              else                    ! Voxel partially overlaid with quadrics
                write(out,'(1x,es12.5,1x,es7.1,1x,es12.5,1x,a1)')
     &            q,2.0*sigma,thismassvox,'-'
              endif
            endif

            ! Evaluate average uncertainty for scores above 1/2 max score:
            if (edep(vox).gt.maxq.and.fact.gt.0.0) then
              avesig = avesig+(sigma/q)**2
              nchan = nchan+1
            endif
          enddo

          if (prtxyz.eq.-1) cycle                 ! No separators in binary mode
          if (xvoxmax.gt.xvoxmin) write(out,*) '' ! Separate 2D data blocks
        enddo

        if (prtxyz.eq.-1) cycle                   ! No separators in binary mode
        if (yvoxmax.gt.yvoxmin) write(out,*) ''   ! Separate 3D data blocks
      enddo

      uncdone = 0
      if (nchan.gt.0) then
        uncert = 200.0*sqrt(avesig/nchan)
        if (uncert.lt.unclimit) uncdone = 2  ! Uncertainty reached
      else
        uncert = 0.0  ! Uncertainty assumed not reached when score is nil
      endif

      if (prtxyz.ne.-1) then ! Generic report if not in binary mode
        write(out,'(a)') ' '
        write(out,'(a)') '# Performance report'
        write(out,'(a)') '#   Memory used (=tally+voxels) (MB):'
        write(out,'(a,3(es12.5,a))') '#   ',(memVDD+memvox)*1.0e-6,
     &    ' (=',memVDD*1.0e-6,' +',memvox*1.0e-6,')'
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


      subroutine VDDinitally(mem)
!*******************************************************************
!*    Initializes the Voxel Dose tally. To be called before TALLY. *
!*    A voxelized geometry must be defined to use this tally.      *
!*                                                                 *
!*    Output:                                                      *
!*      mem: memory used by tally arrays.                          *
!*******************************************************************
      use VDDmod
      use geovoxMod, only : bodymask,nx,ny,nz,nvox
      use dataTypesMod
      implicit none
      real*8 mem

      character*(*) secid,eos
      parameter (secid=
     &  '[SECTION TALLY VOXEL DOSE v.2014-12-27]')
      parameter (eos='[END OF VDD SECTION]')
      character*80 buffer,upcase
      integer errorcode

      write(*,*) ''
      write(*,'(a)')
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      call getline(buffer,0)
      if (index(buffer,secid).eq.0) then
        write(*,'(a)') 'VDDinitally:ERROR: incorrect section header;'
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
        if (bodymask.eq.-1) then
          write(*,'(a)') 'VDDinitally:ERROR: voxel dose tally '//
     &      'is ON but no voxelized geometry has been defined.'
          stop
        endif
      else if (buffer.eq.'OFF') then
        active = .false.
        write(*, '(a)')
     &    '>>>> Tally Voxel Dose Distribution is OFF >>>>'
        do
          read(*,'(a80)',iostat=errorcode) buffer
          if (errorcode.ne.0) then
            write(*,'(a,a,a)') 'VDDinitally:ERROR: ',
     &       'Unable to find End-Of-Section mark: ',eos
            stop
          endif
          if (index(buffer,eos).ne.0) return
        enddo
      else
        write(*,'(a)')
     &    'VDDinitally:ERROR: expecting to find ON or OFF'
        write(*,'(a)') 'found instead:'
        write(*,'(a)') buffer
        stop
      endif

      write(*,'(a)') 'Region Of Interest set to:'
      read(*,*) xvoxmin,xvoxmax
      xvoxmin = max(1,xvoxmin)
      xvoxmax = min(nx,xvoxmax)
      if (xvoxmax.eq.0) xvoxmax = nx
      write(*,'(a,2(2x,i0))') '  x-voxel interval: ',xvoxmin,xvoxmax
      read(*,*) yvoxmin,yvoxmax
      yvoxmin = max(1,yvoxmin)
      yvoxmax = min(ny,yvoxmax)
      if (yvoxmax.eq.0) yvoxmax = ny
      write(*,'(a,2(2x,i0))') '  y-voxel interval: ',yvoxmin,yvoxmax
      read(*,*) zvoxmin, zvoxmax
      zvoxmin = max(1,zvoxmin)
      zvoxmax = min(nz,zvoxmax)
      if (zvoxmax.eq.0) zvoxmax = nz
      write(*,'(a,2(2x,i0))') '  z-voxel interval: ',zvoxmin,zvoxmax
      if (min(xvoxmax-xvoxmin,yvoxmax-yvoxmin,zvoxmax-zvoxmin).lt.0)
     &  then
        write(*,'(a)') 'VDDinitally:ERROR: invalid ROI.'
        stop
      endif

      write(*,'(a)') 'Print voxels mass:'
      read(*,*) prtdens
      if (prtdens.eq.1) then
        write(*,'(a)') ' yes'
      else
        write(*,'(a)') ' no'
      endif

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
          write(*,'(a)') 'VDDinitally:ERROR: expecting 1, 0 or -1.'
          stop
      end select

      write(*,'(a)') 'Relative uncertainty (%) requested:'
      read(*,*) unclimit
      write(*,'(1x,es12.5)') unclimit

      read(*,'(a80)') buffer
      if (index(buffer,eos).eq.0) then
        write(*,'(a)')
     &    'VDDinitally:ERROR: End-Of-Section mark not found'
        write(*,'(a,a)') '  expecting to find: ',eos
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif

      ! Allocate arrays:
      mem = sizeofReal8*4*real(nvox)
      write(*,'(a)') 'Memory used to allocate arrays (MB):'
      write(*,'(1x,es12.5)') mem*1.0e-6
      allocate(nlast(nvox),edptmp(nvox),edep(nvox),edep2(nvox),
     &  stat=errorcode)
      if (errorcode.ne.0) then
        write(*,'(a)') 'VDDinitally:ERROR: not enough memory.'
        stop
      endif
      memVDD = mem  ! Store for later use by report

      ! Init arrays:
      nlast   = 0.0d0
      edptmp  = 0.0
      edep    = 0.0
      edep2   = 0.0

      call VDDdump(0)   ! Sim restart file

      write(*,'(a)') '>>>> VDD tally initialization finished >>>>'
      end


      subroutine VDDdump(mode)
!*******************************************************************
!*    Dumps into or reads data from a dump file.                   *
!*                                                                 *
!*    Input:                                                       *
!*      mode -> 1 to write dump file, else to read from it.        *
!*******************************************************************
      use VDDmod
      use geovoxMod, only : nvox
      use dumpMod
      implicit none
      integer mode

      if (mode.eq.1) then
        if (dumpfile.lt.0) return  ! No dump file open
        write(dumpfile) edep(1:nvox),edep2(1:nvox)
      else
        if (restartfile.lt.0) return  ! No restart file open
        read(restartfile) edep(1:nvox),edep2(1:nvox)
      endif
      end


!>>>> End Of File >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
