C/******************************************************************************
C *
C *  $Id: example_event_generator.F,v 1.2 2006/10/20 23:11:21 capote Exp $
C *
C *  This is an example event generator written in Fortran 77.
C *  It provides a mono-energetic pencil beam with particle type, energy and
C *  position defined in an input file.
C *  This event generator can be loaded only once
C *  (unlike the event generator implemented in example_event_generator.cpp)
C *
C *****************************************************************************/

#ifdef DOUBLE
#define IAEA_Float real*8
#else
#define IAEA_Float real*4
#endif

#define IAEA_I16 integer*2
#define IAEA_I32 integer*4
#define IAEA_I64 integer*8

      subroutine init_source(id,input_file,result)
      implicit none
      IAEA_I32 id,result
      character*(*) input_file
#include "example_event_generator_f.h"
      data initialized/0/
      if( initialized.ne.0 ) then
          write(6,*) 
     &    'Error: this event generator can only be loaded once'
          result = -1
          return
      end if
      open(3,file=input_file,status='old',err=100)
      read(3,*,err=200) type,energy,xo,yo,zo
      result=0
      initialized=1
      id=1
      count=0
      return
100   continue
      write(6,*) 'Error: failed to open ',input_file
      result = -2
      return
200   continue
      write(6,*) 'Error: failed to read ',input_file
      result = -3
      return
      end

      subroutine get_maximum_energy(id,emax)
      implicit none
      IAEA_I32   id
      IAEA_Float emax
#include "example_event_generator_f.h"
      if( id.eq.1 ) emax = energy
      return
      end

      subroutine get_minimum_energy(id,emin)
      implicit none
      IAEA_I32   id
      IAEA_Float emin
#include "example_event_generator_f.h"
      if( id.eq.1 ) emin = energy
      return
      end

      subroutine get_extra_numbers(id,nf,ni)
      implicit none
      IAEA_I32   id,nf,ni
      nf = 0
      ni = 0
      return
      end

      subroutine get_type_extra_long_variable(id,ind,t)
      implicit none
      IAEA_I32   id,ind,t
      t=-1
      return
      end

      subroutine get_type_extra_float_variable(id,ind,t)
      implicit none
      IAEA_I32   id,ind,t
      t=-1
      return
      end
      
      subroutine get_nstat(id,nstat)
      implicit none
      IAEA_I32 id
      IAEA_I64 nstat
#include "example_event_generator_f.h"
      if( id.ne.1.or.initialized.ne.1 ) then
          nstat=-1
      else
          nstat=count
      endif
      return
      end

      subroutine set_parallel_run(id,ipar,ichunk,nchunk,result)
      implicit none
      IAEA_I32 id,ipar,ichunk,nchunk,result
      result = 0
      return
      end

      subroutine get_next_particle(id,nstat,q,E,wt,x,y,z,u,v,w,ef,el)
      implicit none
      IAEA_I32 id,nstat,q,el(*)
      IAEA_Float E,wt,x,y,z,u,v,w,ef(*)
#include "example_event_generator_f.h"
      if( id.ne.1.or.initialized.ne.1 ) return
      nstat=1
      q=type
      E=energy
      wt=1
      x=xo
      y=yo
      z=zo
      u=0
      v=0
      w=1
      return
      end

      subroutine destroy_source(id,result) 
      implicit none
      IAEA_I32 id,result
#include "example_event_generator_f.h"
      if( id.ne.1.or.initialized.ne.1 ) then
          result=-1
          return
      endif
      result=0
      initialized=0
      return
      end
