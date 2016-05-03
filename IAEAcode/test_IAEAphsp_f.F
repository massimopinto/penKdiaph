C Sources files for the interfase (not tested with event generators):
C iaea_header.cpp (iaea_header.h)
C iaea_phsp.cpp (iaea_phsp.h)
C iaea_record.cpp (iaea_record.h)
C utilities.cpp (utilities.h)
C
C Test files: test_IAEAphsp.cpp (C++) or test_IAEAphsp_f.F (FORTRAN)
C
C To compile in Linux without make 
C g77 test_IAEAphsp_f.F iaea_header.cpp iaea_phsp.cpp iaea_record.cpp utilities.cpp -lm -lstdc++ -o test_IAEAphsp
C g77 test_IAEAphsp.c iaea_header.cpp iaea_phsp.cpp iaea_record.cpp utilities.cpp -lm -lstdc++ -o test_IAEAphsp
C Tested in RED HAT LINUX with compilers gcc,cc,g95,icc
C
C If you have GNU make available the corresponding makefiles are provided 
C In that case the executables and DLL can be produced (for a given compiler/OS)  
C make -f Makefile_... (depends on your compiler/OS)
C make -f Makefile_Linux_g++_g77 (an example in Linux using g77/g++ compilers)
C

#ifdef DOUBLE
#define IAEA_Float real*8
#else
#define IAEA_Float real*4
#endif

#define IAEA_I16 integer*2
#define IAEA_I32 integer*4
#define IAEA_I64 integer*8

#define NUM_EXTRA_LONG 1
#define NUM_EXTRA_FLOAT 1

      program test

      character*80  file_name
      integer*4     i
      IAEA_I32      source_read,source_write,access_read,access_write,
     &              result
      parameter    (access_read=1,access_write=2)
      IAEA_I32      typ,extra_ints(NUM_EXTRA_LONG),n_stat,
     &              iextrafloat,iextralong,ind,extralong_type
      IAEA_Float    E, wt, x, y, z, u, v, w, z_constant, 
     &              extra_floats(NUM_EXTRA_FLOAT)
      IAEA_I64      histories,nphot,irecord,read_indep,nrecorded,norig

      do i=1,len(file_name)
          file_name(i:i) = char(0)
      end do

      !
      ! ***** Create a new source for reading
      !       Because file_name is empty, iaea_new_source will ask for 
      !       a name
      !
      call iaea_new_source(source_read,file_name,access_read,result)
      if( result.lt.0 ) then
          write(6,*) 'Error creating new source for reading particles'
          call exit(1)
      end if
      
      !
      ! ***** Get the total number of particles in the file
      !
      typ = -1  ! could be 0
      call iaea_get_max_particles(source_read, typ, histories)
      write(6,*) 'Total number of particles: ',histories
      !
      ! ***** Get the number of photons in the file
      !
      typ = 1  
      call iaea_get_max_particles(source_read, typ, nphot)
      write(6,*) 'Number of photons: ',nphot
      !
      ! ***** Get the number of original statistically independent particles
      !
      call iaea_get_total_original_particles(source_read,norig)
      write(6,*) 'Original particles: ',norig

      do i=1,len(file_name)
          file_name(i:i) = char(0)
      end do
      
      !
      ! ***** Create a new source for writing
      !       Because file_name is empty, iaea_new_source will ask for 
      !       a name
      !
      call iaea_new_source(source_write,file_name,access_write,result)
      if( result.lt.0 ) then
          write(6,*) 'Error creating new source for writing particles'
          call exit(1)
      end if
      !
      ! ***** Set the number of extra variables to be stored:
      !       0 extra floats and 1 extra long
      !
      iextrafloat = 0
      iextralong = 1
      call iaea_set_extra_numbers(source_write,iextrafloat,iextralong)
      !
      ! ***** Set the type of first extra long being stored
      !
      ind = 0
      extralong_type = 2  ! LATCH
      call iaea_set_type_extralong_variable(source_write,ind,
     &                                      extralong_type)
      !
      ! ***** Set z as a constant (i.e. it will not be stored in the file)
      !
      ind = 2
      z_constant = 80
      call iaea_set_constant_variable(source_write,ind,z_constant)
      
      !
      ! ***** Copy the header of the read source to the write source
      !
      call iaea_copy_header(source_read,source_write,result)

      !
      ! ***** Read particles from the read source and write them in the 
      !       write source
      !
      nrecorded = 0
      do irecord=1,histories
          !
          ! ***** read a particle
          !
          call iaea_get_particle(source_read,n_stat,typ,E,wt,x,y,z,
     &              u,v,w,extra_floats,extra_ints)
          if( n_stat.eq.-1 ) exit  ! i.e. end of file
          !
          ! ***** write a particle
          !
          call iaea_write_particle(source_write,n_stat,typ,E,wt,x,y,z,
     &              u,v,w,extra_floats,extra_ints)
          nrecorded = nrecorded + 1
      end do

      call iaea_print_header(source_write,result)

      write(6,*) 'Particles recorded: ',nrecorded

      !
      ! ***** Close sources
      !
      call iaea_destroy_source(source_read,result)
      call iaea_destroy_source(source_write,result)

      end
      
