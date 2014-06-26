
/*
  Mistral is a constraint satisfaction and optimisation library
  Copyright (C) 2003-2005  Emmanuel Hebrard
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

  The author can be contacted electronically at 
  ehebrard@cse.unsw.edu.au.
*/


#include <mistral_glo.h>

#include <time.h>
#include <sys/times.h>
#include <sys/resource.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
//#include <cstring>
#include <vector>
#include <iostream>

static unsigned mistral_rand_x[56], mistral_rand_y[256], mistral_rand_z;
static int mistral_rand_j, mistral_rand_k;

//using namespace Mistral;

/**********************************************
 * Knuth's Random number generator (code from sp-1.4)
 **********************************************/
/**
   urand(), urand0() return uniformly distributed unsigned random ints
   all available bits are random, e.g. 32 bits on many platforms
   usrand(seed) initializes the generator
   a seed of 0 uses the current time as seed
   
   urand0() is the additive number generator "Program A" on p.27 of Knuth v.2
   urand() is urand0() randomized by shuffling, "Algorithm B" on p.32 of Knuth v.2
   urand0() is one of the fastest known generators that passes randomness tests
   urand() is somewhat slower, and is presumably better
*/

void Mistral::usrand (unsigned seed)
{
  if(seed == 0) {
    std::cerr << "c Warning: 0 is not a valid seed!!" << std::endl;
  }

  int j;
  mistral_rand_x[1] = 1;
  if (seed) mistral_rand_x[2] = seed;
  else mistral_rand_x[2] = time (NULL);
  for (j=3; j<56; ++j) mistral_rand_x[j] = mistral_rand_x[j-1] + mistral_rand_x[j-2];
  
  mistral_rand_j = 24;
  mistral_rand_k = 55;
  for (j=255; j>=0; --j) urand0 ();
  for (j=255; j>=0; --j) mistral_rand_y[j] = urand0 ();
  mistral_rand_z = urand0 ();
}

unsigned Mistral::urand0 (void)
{
  if (--mistral_rand_j == 0) mistral_rand_j = 55;
  if (--mistral_rand_k == 0) mistral_rand_k = 55;
  return mistral_rand_x[mistral_rand_k] += mistral_rand_x[mistral_rand_j];
}

unsigned Mistral::urand (void)
{
  int j;

  j =  mistral_rand_z >> 24;

  mistral_rand_z = mistral_rand_y[j];
  if (--mistral_rand_j == 0) mistral_rand_j = 55;
  if (--mistral_rand_k == 0) mistral_rand_k = 55;
  mistral_rand_y[j] = mistral_rand_x[mistral_rand_k] += mistral_rand_x[mistral_rand_j];
  return mistral_rand_z;
}

int Mistral::randint(const int upto)
//gives a random integer uniformly in [0..upto-1]
{
  return(urand() % upto);
}

double Mistral::randreal()
//gives a random integer uniformly in (0,1)
{
  return (urand()/(double)((double)MAX_URAND+1.0));	
}



/**********************************************
 * Timing Function
 *********************************************/


/// Return the ellapsed cpu time
double Mistral::getRunTime() {
  double df = 0;

#ifdef _UNIX
  struct tms usage;
  static int clock_ticks = sysconf(_SC_CLK_TCK);
  times(&usage);
  df=((double)usage.tms_utime+(double)usage.tms_stime)/clock_ticks;
#endif

  return df;
}

unsigned long int Mistral::getMemory() {
  unsigned long int mem = 0;
//   struct rusage usage;
//   getrusage(RUSAGE_CHILDREN, &usage);
//   mem = usage.ru_maxrss;
//   return mem;

#ifdef _LINUX
  
  char buf[30];
  snprintf(buf, 30, "/proc/%u/statm", (unsigned)getpid());
  FILE* pf = fopen(buf, "r");
  
  if (pf) {
    unsigned int size; //       total program size
    //unsigned resident;//   resident set size
    //unsigned share;//      shared pages
    //unsigned text;//       text (code)
    //unsigned lib;//        library
    //unsigned data;//       data/stack
    //unsigned dt;//         dirty pages (unused in Linux 2.6)
    fscanf(pf, "%u" /* %u %u %u %u %u"*/, &size/*, &resident, &share, &text, &lib, &data*/);
    //DOMSGCAT(MSTATS, std::setprecision(4) << size / (1024.0) << "MB mem used");
    mem = size;
  }  
  fclose(pf);
  
#endif

  return mem;
}



/**********************************************
 * Command line parsing
 *********************************************/

void Mistral::getCommandLine(const char** int_ident, 
		    int*    int_param,
		    int     int_nb,
		    const char** str_ident,
		    const char** str_param,
		    int     str_nb,
		    char**  argv,
		    int     argc)
{
  std::fill(int_param,int_param+int_nb,NOVAL);
  std::fill(str_param,str_param+str_nb,"nil");
  int i,j,k;
  for(i=1; i<argc-1;) {
    k=1;
    j = str_nb;    
    while( j-- ) if( !strcmp(str_ident[j],argv[i]) ) {
      str_param[j] = argv[i+1];
      k=0;
      break;
    }
    if(k) {
      j = int_nb;
      while( j-- ) {
	if( !strcmp(int_ident[j],argv[i]) ) {
	  int_param[j] = atoi(argv[i+1]);
	  break;
	}
      }
    }
    i+=2;
  }
}

