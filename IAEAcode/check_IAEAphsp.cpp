/*
 * Copyright (C) 2006 International Atomic Energy Agency
 * -----------------------------------------------------------------------------
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is furnished
 * to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 *-----------------------------------------------------------------------------
 *
 *   AUTHORS:
 *
 *   Roberto Capote Noy, PhD
 *   e-mail: R.CapoteNoy@iaea.org (rcapotenoy@yahoo.com)
 *   International Atomic Energy Agency
 *   Nuclear Data Section, P.O.Box 100
 *   Wagramerstrasse 5, Vienna A-1400, AUSTRIA
 *   Phone: +431-260021713; Fax: +431-26007
 *
 *   Iwan Kawrakow, PhD
 *   e-mail iwan@irs.phy.nrc.ca
 *   Ionizing Radiation Standards
 *   Institute for National Measurement Standards 
 *   National Research Council of Canada Ottawa, ON, K1A 0R6 Canada
 *   Phone: +1-613-993 2197, ext.241; Fax: +1-613-952 9865
 *
 **********************************************************************************
 * For documentation
 * see http://www-nds.iaea.org/reports-new/indc-reports/indc-nds/indc-nds-0484.pdf
 **********************************************************************************/
//
// Sources files for the interfase (not tested with event generators):
// iaea_header.cpp (iaea_header.h)
// iaea_phsp.cpp (iaea_phsp.h)
// iaea_record.cpp (iaea_record.h)
// utilities.cpp (utilities.h)
//
// Test files: test_IAEAphsp.cpp (C++) or test_IAEAphsp_f.F (FORTRAN)
//
// Test files: check_IAEAphsp.cpp (C++) 
//
// To compile in Linux without make 
// cc test_IAEAphsp.cpp iaea_header.cpp iaea_phsp.cpp iaea_record.cpp utilities.cpp -lm -lstdc++ -o test_IAEAphsp
// cc check_IAEAphsp.cpp iaea_header.cpp iaea_phsp.cpp iaea_record.cpp utilities.cpp -lm -lstdc++ -o check
//
// Tested in RED HAT LINUX with compilers gcc,cc,g95,icc
//
// Tested in Windows-7,Windows-XP with Microsoft Visual Studio compiler v4.2
//
// If you have GNU make available the corresponding makefiles are provided 
// In that case the executables and DLL can be produced (for a given compiler/OS)  
// make -f Makefile_... (depends on your compiler/OS)
// make -f Makefile_Linux_g++_g77 (an example in Linux using g77/g++ compilers)
//
#ifdef WIN32 
#include <iostream>  // so that namespace std becomes defined
#endif
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cctype> 

#include "iaea_phsp.h"

#if !(defined WIN32) && !(defined WIN64)
using namespace std;
#endif

int main()
{
   // Declaring IAEA phsp and access rights 
   IAEA_I32 source_read = 0, source_write = 1; 
   IAEA_I32 access_read = 1, access_write = 2;

   IAEA_I32 result; 

   // Particle properties
   IAEA_I32 type, index, *extra_ints, *extralong_types, *extrafloat_types;
   IAEA_Float E, wt, x, y, z, u, v, w, *extra_floats;


   // Initializing IAEA source for reading and retrieving logical information
   IAEA_I32 len = 81;
   char fiaearead[81]={'\0'};
   result = -1; 
   iaea_new_source(&source_read , fiaearead , &access_read , &result, len);

   result = -1; // To request for the total number of particles

   iaea_check_file_size_byte_order(&source_read, &result);
   
   if(result != 0) {
      printf("\n iaea_check_file_size_byte_order result %d, ERROR \n",result);
      return(-1);
   }

   IAEA_I64 histories; // Getting total number of histories

   result = -1; // To request for the total number of particles
   iaea_get_max_particles(&source_read, &result, &histories); 
   printf("\n Total number of histories: %d\n",histories);

   // Initializing IAEA source for writing 
   char fiaeawrite[81]={'\0'};
   iaea_new_source(&source_write, fiaeawrite, &access_write, &result, len ); 

   // copying "source_read" header to "source_write" header

   iaea_copy_header(&source_read, &source_write, &result);

   // Setting i/o flags from the source

   IAEA_I32 iextrafloat,iextralong; 

   iaea_get_extra_numbers(&source_read , &iextrafloat, &iextralong);
   iaea_set_extra_numbers(&source_write, &iextrafloat, &iextralong);

   if (iextrafloat>0) extra_floats = (IAEA_Float *) calloc(iextrafloat,sizeof(IAEA_Float));  
   if (iextralong>0 ) extra_ints = (IAEA_I32 *) calloc(iextralong,sizeof(IAEA_I32));  

   if (iextrafloat>0) extrafloat_types = (IAEA_I32 *) calloc(iextrafloat,sizeof(IAEA_I32));  
   if (iextralong>0 ) extralong_types  = (IAEA_I32 *) calloc(iextralong,sizeof(IAEA_I32));  

   IAEA_I32 itmp;
   iaea_get_type_extra_variables(&source_read, &result, extralong_types, extrafloat_types);

   for(itmp=0; itmp<iextralong; itmp++)
	   iaea_set_type_extralong_variable (&source_write, &itmp, &extralong_types[itmp]);
   for(itmp=0; itmp<iextrafloat; itmp++)
   	   iaea_set_type_extrafloat_variable(&source_write, &itmp, &extrafloat_types[itmp]);

   IAEA_Float z_constant;

   for (index=0;index<7;index++) {
     iaea_get_constant_variable(&source_read , &index, &z_constant, &result);
   	 if(result==0) iaea_set_constant_variable(&source_write, &index, &z_constant); 
   }

   IAEA_I64 tot_hist = histories;

   /* The commented block below could be used to retrieve chunks of the phsp
   
   // The starting point set at the second chunk of the phsp (i_chunk =2)
   IAEA_I32 i_parallel = 1, i_chunk = 2, n_chunk = 3;
   iaea_set_parallel(&source_read, 
	      &i_parallel, &i_chunk, &n_chunk, &result);

   histories /= n_chunk; // Only histories/n_chunk number is read

   printf("\n Phase space divided in %d pieces of %d records each \n",
	   n_chunk,histories);   

   printf(" Start reading at record number %i up to record number %i\n",
	   (i_chunk-1)*histories,i_chunk*histories);   

   */

  int iaea_charge[5]={0,-1,+1,0,+1};

   IAEA_I32 n_stat;
   IAEA_I64 read_indep = 0, nrecorded = 0, irecord;
   for(irecord=0; irecord<histories; irecord++)
   {
	   // read IAEA particle
	   iaea_get_particle(&source_read, &n_stat,
              &type, /* particle type */
              &E   , /* kinetic energy in MeV */
              &wt  , /* statistical weight */
              &x   ,
              &y   ,
              &z   , /* position in cartesian coordinates*/
              &u   ,
              &v   ,
              &w   , /* direction in cartesian coordinates*/
              extra_floats, extra_ints);

 	   if(irecord<30) 
		   printf(" %1i %2i %6.2f %6.2f %6.2f %5.3f %5.3f %5.3f %6.3f %5.3f %4i\n",
		   type,iaea_charge[type-1],x,y,z,u,v,w,E,wt,n_stat);

	   if( n_stat > 0 ) read_indep++; 

       if( n_stat == -1 ) break; 
   
	   // write IAEA particle
	   iaea_write_particle(&source_write, &n_stat,
              &type, /* particle type */
              &E   , /* kinetic energy in MeV */
              &wt  , /* statistical weight */
              &x   ,
              &y   ,
              &z   , /* position in cartesian coordinates*/
              &u   ,
              &v   ,
              &w   , /* direction in cartesian coordinates*/
              extra_floats, extra_ints);


           nrecorded++;
           
		   // if(nrecorded>100) break;

   }
   // iaea_print_header(&source_read, &result);  
   printf("\n Number of independent histories read/written %d \n",read_indep);
   printf(" Total number of histories written %d \n",nrecorded);


   iaea_print_header(&source_write, &result);  

   // Closing sources
   iaea_destroy_source(&source_read,&result);

   iaea_destroy_source(&source_write,&result);


   if (iextrafloat>0) free(extrafloat_types);

   if (iextralong>0 ) free(extralong_types);

   if (iextrafloat>0) free(extra_floats);

   if (iextralong>0 ) free(extra_ints);


   // *************************************************************
   printf("\n\n Normal Program Termination\n");

   return(1);
}