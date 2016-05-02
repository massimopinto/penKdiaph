!*******************************************************************
!*                          TALLY                                  *
!*                 PARTICLE TRACK STRUCTURE                        *
!*                                                                 *
!* Short description:                                              *
!*   Writes to a file the position and value of energy loss events *
!*   so that particle tracks can be visualized.                    *
!* Rev:                                                            *
!*   2015-04-30                                                    *
!*******************************************************************


!*******************************************************************
!*******************************************************************
!*    Modules section                                              *
!*******************************************************************
!*******************************************************************

      module PTSmod
!*******************************************************************
!*    Vars for the PTS tally.                                      *
!*******************************************************************
      implicit none
      save
      logical active
      integer ptsunit
      integer*4 ntrack,ntrackmax
      real*8 xlast,ylast,zlast
      end


!*******************************************************************
!*******************************************************************
!*    Routines section                                             *
!*******************************************************************
!*******************************************************************

      subroutine PTStally(mode,arg)
!*******************************************************************
!*    Input:                                                       *
!*      mode -> Identifies the state of the calling routine        *
!*      arg -> energy loss (mode<0) or history no. (mode=1)        *
!*******************************************************************
      use TRACK_mod
      use RSEEDcommonMod
      use dataTypesMod
      use PTSmod
      implicit none
      integer mode
      real*8 arg

      real*8 de

      if (.not.active) return

      if (mode.eq.1) then      ! Start of a new history
        ntrack = ntrack+1
        if (ntrack.gt.ntrackmax) then
          write(ptsunit,*) ''
          write(ptsunit,'(a)') '# Done'
          active = .false.     ! Job done, deactivate tally
          return
        endif
        write(ptsunit,*) ''
        write(ptsunit,'(a)') '# New history started:'
        write(ptsunit,'(a,i0)') '#  ',int(arg,kind=kindInt8)
        write(ptsunit,'(a)') '# Random seeds:'
        write(ptsunit,'(a,i0,1x,i0)') '#  ',seed1,seed2
        return
      endif

      if (mode.lt.0) then
        select case (mode)
        case (-99)             ! A new particle starts
          if (mat.eq.0) return ! Particle not in object, no report
          de = 0.0
          write(ptsunit,*) ''
          write(ptsunit,*) ''  ! Separate data sets with double blank
        case (-96)             ! Particle splitting
          de = 0.0             ! No actual Eloss in this case
        case default
          de = arg             ! Energy scored
        end select

      else if (mode.eq.4) then ! Entered a new material
        if (mat.eq.0) then     ! The particle escaped
          de = 0.0
          x = xlast+arg*u      ! Move the particle up to the exit surface
          y = ylast+arg*v
          z = zlast+arg*w
        else                   ! Still in object, report new material
          de = 0.0
        endif

      else if (mode.eq.5) then ! Particle killed
        de = 0.0               ! No actual Eloss

      else
        return                 ! No need to report
      endif

      xlast = x                ! Store last position for future use
      ylast = y
      zlast = z

      write(ptsunit,'(3(i0,1x),9(es12.5,1x),i3,1x,i0)')
     &  kpar,ibody,mat,x,y,z,u,v,w,e,de,wght,mode,ilb(5)
      end


      subroutine PTSinitally
!*******************************************************************
!*    Initializes. To be called before TALLY                       *
!*******************************************************************
      use PTSmod
      implicit none
      character*(*) secid,eos
      parameter (secid=
     &'[SECTION TALLY PARTICLE TRACK STRUCTURE v.2009-06-15]')
      parameter (eos='[END OF PTS SECTION]')
      character*80 buffer,upcase
      integer finduf,error

      write(*,*) ' '
      write(*,'(a)')
     & '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      call getline(buffer,0)
      if (index(buffer,secid).eq.0) then
        write(*,'(a)') 'PTSinitally:ERROR: incorrect section header;'
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
     &    '>>>> Tally Particle Track Structure is OFF >>>>'
        do
          read(*,'(a80)',iostat=error) buffer
          if (error.ne.0) then
            write(*,'(a,a,a)') 'PTSinitally:ERROR: ',
     &       'Unable to find End-Of-Section mark: ',eos
            stop
          endif
          if (index(buffer,eos).ne.0) return
        enddo
      else
        write(*,'(a)')
     &    'PTSinitally:ERROR: expecting to find ON or OFF'
        write(*,'(a)') 'found instead:'
        write(*,'(a)') buffer
        stop
      endif

      ptsunit = finduf()
      buffer = 'tallyParticleTrackStructure.dat'
      open(ptsunit,file=buffer,iostat=error)
      if (error.ne.0) then
        write(*,'(a)')
     &    'PTSinitally:ERROR: cannot open track data file:'
        write(*,'(a80)') buffer
        stop
      endif
      write(ptsunit,'(a)') '# kpar : body : mat : x : y : z : '//
     &   'u : v : w : e : elost : wght : tallymode : ilb(5)'
      write(*,'(a)') 'Number of history tracks to display:'
      read(*,*) ntrackmax
      write(*,'(1x,i0)') ntrackmax
      if (ntrackmax.gt.10000) then
        write(*,'(a)')
     &  'PTSinitally:WARNING: No. tracks to report set to max=10000'
        ntrackmax = 10000
      endif

      ! Init vars:
      ntrack = 0

      read(*,'(a80)') buffer
      if (index(buffer,eos).eq.0) then
        write(*,*) 'PTSinitally:ERROR: End-Of-Section mark not found'
        write(*,'(a,a)') '  expecting to find: ',eos
        write(*,'(a,a)') '  found instead:     ',buffer
        stop
      endif
      write(*,'(a)') '>>>> PTS tally initialization finished >>>>'
      end


!>>>> End Of File >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

