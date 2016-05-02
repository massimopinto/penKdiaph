!*******************************************************************
!*                         TIMING                                  *
!*                                                                 *
!* Short description:                                              *
!*   Time routines complying with the F95 standard.                *
!*                                                                 *
!*   The following implicit routines are used:                     *
!*    - date_and_time(char*8 sdate,char*10 stime,char*5 zone,      *
!*                  int values(8))                                 *
!*      where                                                      *
!*      sdate = YYYYMMDD                                           *
!*      stime = hhmmss.sss                                         *
!*      zone = hhmm , difference with respect to UTC               *
!*      values(1) = year,                                          *
!*      values(2) = month (1..12)                                  *
!*      values(3) = day (1..31)                                    *
!*      values(4) = difference with respect to UTC in minutes      *
!*      values(5) = hour of the day (0..23)                        *
!*      values(6) = minutes of the hour                            *
!*      values(7) = seconds of the minute                          *
!*      values(8) = milliseconds of the second                     *
!*                                                                 *
!*    - cpu_time(real time) (in s)                                 *
!*******************************************************************


      subroutine initime
!*******************************************************************
!*    Writes the current date to stdout in text format.            *
!*    Intended to be called when the main program starts.          *
!*******************************************************************
      implicit none
      character*100 fdate
      character*8 sdate
      character*10 stime
      character*5 zone
      integer values(8)

      call date_and_time(sdate,stime,zone,values)
      call dateString(values,fdate)
      write(*,*) ''
      write(*,'(a,a)') 'Program started on ',fdate
      end


      subroutine endtime
!*******************************************************************
!*    Writes current date to stdout in text format.                *
!*    Intended to be called when the main program ends.            *
!*******************************************************************
      implicit none
      character*100 fdate
      character*8 sdate
      character*10 stime
      character*5 zone
      integer values(8)

      call date_and_time(sdate,stime,zone,values)
      call dateString(values,fdate)
      write(*,*) ''
      write(*,'(a,a)') 'Program ended on ',fdate
      write(*,*) ''
      write(*,'(a)') 'Have a nice day.'
      end


      subroutine treset
!*******************************************************************
!*    Resets the real and CPU timers                               *
!*******************************************************************
      implicit none
      real utime0
      real*8 rtime0,lasthour,tshift
      common /tim001/ lasthour,tshift,rtime0,utime0
      real*8 drtime,dutime
      common /timrestart/ drtime,dutime
      character*8 sdate
      character*10 stime
      character*5 zone
      integer values(8)

      call date_and_time(sdate,stime,zone,values)
      lasthour = values(5)
      rtime0 = (((values(5)*60.0d0)+values(6))*60.0d0)+values(7)+
     &            1.0d-3*values(8)
      call cpu_time(utime0)
      if (utime0.lt.0.0) then
        write(*,*) 'treset:ERROR: intrinsic time routine failed'
        write(*,*) '        Timings will not be reliable.'
      endif
      tshift = 0.0d0  ! Days passed in s.
      drtime = 0.0d0  ! Elapsed real time in previous runs (restart feature)
      dutime = 0.0d0  ! Elapsed user time in previous runs (restart feature)
      end


      real*8 function cputime()
!*******************************************************************
!*    Returns CPU (user) time elapsed since last call to TRESET    *
!*                                                                 *
!*    Output:                                                      *
!*      CPU time in s                                              *
!*******************************************************************
      implicit none
      real utime0
      real*8 rtime0,lasthour,tshift
      common /tim001/ lasthour,tshift,rtime0,utime0
      real*8 drtime,dutime
      common /timrestart/ drtime,dutime
      real time

      call cpu_time(time)
      if (time.lt.0.0) then
        write(*,*) 'cputime:ERROR: intrinsic time routine failed'
        write(*,*) '         Timings will not be reliable.'
      endif
      cputime = time-utime0+dutime
      end


      real*8 function realtime()
!*******************************************************************
!*    Returns the real time elapsed since last call to TRESET      *
!*                                                                 *
!*    Output:                                                      *
!*      Real time in s                                             *
!*    Comments:                                                    *
!*      -> MUST be called at least once per day or it will not     *
!*         return reliable data                                    *
!*******************************************************************
      implicit none
      real utime0
      real*8 rtime0,lasthour,tshift
      common /tim001/ lasthour,tshift,rtime0,utime0
      real*8 drtime,dutime
      common /timrestart/ drtime,dutime
      character*8 sdate
      character*10 stime
      character*5 zone
      integer values(8)

      call date_and_time(sdate,stime,zone,values)
      if (dble(values(5)).lt.lasthour) tshift = tshift+86400.0d0  ! The day changed
      lasthour = values(5)  ! Last time (hour of day) realtime() was called
      realtime = (((values(5)*60.0d0)+values(6))*60.0d0)+values(7)+
     &            1.0d-3*values(8)+tshift-rtime0+drtime
      end


      subroutine dateString(values,fdate)
!*******************************************************************
!*    Returns a string with the date and time                      *
!*                                                                 *
!*    Input:                                                       *
!*      values(8) -> date array in standard Fortran (see header)   *
!*    Output:                                                      *
!*      fdate -> formatted date (and time) in a char string        *
!*******************************************************************
      implicit none
      integer values(8)
      character*100 fdate
      character*3 month(12)

      month(1)  = 'Jan'
      month(2)  = 'Feb'
      month(3)  = 'Mar'
      month(4)  = 'Apr'
      month(5)  = 'May'
      month(6)  = 'Jun'
      month(7)  = 'Jul'
      month(8)  = 'Aug'
      month(9)  = 'Sep'
      month(10) = 'Oct'
      month(11) = 'Nov'
      month(12) = 'Dec'
      write(fdate,'(i2,1x,a3,1x,i4,2x,i2.2,a1,i2.2,a1,i2.2)')
     &  values(3),month(values(2)),values(1),values(5),':',
     &  values(6),':',values(7)
      end


      subroutine trestart(rtime,utime)
!*******************************************************************
!*    Resets the real and CPU reference times, useful when program *
!*    restarts. TRESET() must be called first.                     *
!*                                                                 *
!*    Input:                                                       *
!*      rtime -> restart initial value of real time (s).           *
!*      utime -> restart initial value of user time (s).           *
!*******************************************************************
      implicit none
      real*8 rtime,utime

      real*8 drtime,dutime
      common /timrestart/ drtime,dutime

      drtime = rtime
      dutime = utime
      end


!>>>> End Of File >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
