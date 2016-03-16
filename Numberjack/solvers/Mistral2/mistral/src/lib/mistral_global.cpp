
/*
  Mistral 2.0 is a constraint satisfaction and optimisation library
  Copyright (C) 2009  Emmanuel Hebrard
  
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


//#include <mistral_constraint.hpp>
//#include <mistral_structure.hpp>
//include <mistral_backtrack.hpp>
//#include <mistral_variable.hpp>
#include <mistral_global.hpp>


#include <time.h>
#include <sys/times.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdlib.h>
#include <vector>
#include <sstream>
#include <cstring>


/*************** TIMING AND MEMORY USAGE ROUTINES, FROM MINISAT ******/
#define __STDC_FORMAT_MACROS
#include <inttypes.h>


bool probe_macosx() { std::cout << " c Run on mac osx" << std::endl; return true; }
bool probe_linux() { std::cout << " c Run on linux" << std::endl; return true; }
bool probe_win() { std::cout << " c Run on windows" << std::endl; return true; }

bool Mistral::probe() {
#ifdef __APPLE__
  return probe_macosx();
#elif defined __linux__
  return probe_linux();
#elif defined _WIN32 || defined _WIN64
  return probe_win();
#else
#error "unknown platform"
#endif
}

#ifdef _MSC_VER
#include <ctime>

double Mistral::cpu_time(void) {
    return (double)clock() / CLOCKS_PER_SEC; }
#else

#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>

double Mistral::cpu_time(void) {
    struct rusage ru;
    getrusage(RUSAGE_SELF, &ru);
    return (double)ru.ru_utime.tv_sec + (double)ru.ru_utime.tv_usec / 1000000; }
#endif


#if defined(__linux__)
#include <stdio.h>
int Mistral::mem_read_stat(int field)
{
    char    name[256];
    pid_t pid = getpid();
    sprintf(name, "/proc/%d/statm", pid);
    FILE*   in = fopen(name, "rb");
    if (in == NULL) return 0;
    int     value;
    for (; field >= 0; field--)
        fscanf(in, "%d", &value);
    fclose(in);
    return value;
}
uint64_t Mistral::mem_used() { return (uint64_t)mem_read_stat(0) * (uint64_t)getpagesize(); }


#elif defined(__FreeBSD__)
uint64_t Mistral::mem_used(void) {
    struct rusage ru;
    getrusage(RUSAGE_SELF, &ru);
    return ru.ru_maxrss*1024; }


#else
uint64_t Mistral::mem_used() { return 0; }
#endif
/*************** TIMING AND MEMORY USAGE ROUTINES, FROM MINISAT ******/


// std::string Mistral::event2str(Mistral::Event e) {
//   //std::cout
//   std::string str_out = "no event";
//   if(ASSIGNED(e)) {
//     str_out = "value";
//     if(!LB_CHANGED(e)) str_out += " & ub";
//     else if(!UB_CHANGED(e)) str_out += " & lb";
//     else str_out += " & bounds";
//     str_out += " changed";
//   } else if(RANGE_CHANGED(e)) {
//     if(!LB_CHANGED(e)) str_out = "ub";
//     else if(!UB_CHANGED(e)) str_out = "lb";
//     else str_out = "bounds";
//     str_out += " changed";
//   } else if(DOMAIN_CHANGED(e)) {
//     str_out = "domain changed";
//   }
// }

std::string Mistral::event2str(Mistral::Event e) {
  std::string str_evt = "no event";
  if(ASSIGNED(e)) {
    str_evt = "value";
    if(e != VALUE_C) {
      if(!(LB_CHANGED(e))) str_evt += " & ub";
      else if(!(UB_CHANGED(e))) str_evt += " & lb";
      else str_evt += " & bounds";
    } 
  } else if(RANGE_CHANGED(e)) {
    if(!(LB_CHANGED(e))) str_evt = "ub";
    else if(!(UB_CHANGED(e))) str_evt = "lb";
    else str_evt = "bounds";
  } else if(DOMAIN_CHANGED(e)) {
    str_evt = "domain";
  }
  return str_evt;
}

std::string Mistral::event2strc(Mistral::Event e) {
  std::string str_evt = "no";
  if(ASSIGNED(e)) {
    str_evt = "v";
    if(e != VALUE_C) {
      if(!(LB_CHANGED(e))) str_evt += "&u";
      else if(!(UB_CHANGED(e))) str_evt += "&l";
      else str_evt += "&b";
    } 
  } else if(RANGE_CHANGED(e)) {
    if(!(LB_CHANGED(e))) str_evt = "u";
    else if(!(UB_CHANGED(e))) str_evt = "l";
    else str_evt = "b";
  } else if(DOMAIN_CHANGED(e)) {
    str_evt = "d";
  }
  return str_evt;
}


std::string Mistral::outcome2str(Mistral::Outcome e) {
  std::string str_out;
  switch(e) {
  case SAT: str_out = "SAT"; break;
  case OPT: str_out = "OPT"; break;
  case UNSAT: str_out = "UNSAT"; break;
  case UNKNOWN: str_out = "UNKNOWN"; break;
  case LIMITOUT: str_out = "LIMITOUT"; break;
  }
  return str_out;
}

std::string Mistral::domain2str(int d) {
  std::string str_out;
  switch(d) {
  case CONST_VAR: str_out = "constant"; break;
  case BOOL_VAR: str_out = "bool(-)"; break;
  case RANGE_VAR: str_out = "range"; break;
  case BITSET_VAR: str_out = "bitset"; break;
  case LIST_VAR: str_out = "list"; break;
  case VIRTUAL_VAR: str_out = "virtual"; break;
  case EXPRESSION: str_out = "expression"; break;
  case DYN_VAR: str_out = "dynamic"; break;
  default: str_out = "bool";
  }
  return str_out;
}

static unsigned mistral_rand_x[56], mistral_rand_y[256], mistral_rand_z;
static int mistral_rand_j, mistral_rand_k;
static int rand_initialised = false;

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

bool Mistral::random_generator_is_ready() {
  return rand_initialised;
}

void Mistral::usrand (unsigned seed)
{
  int j;

  rand_initialised = true;

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
double Mistral::get_run_time() {
  double df = 0;

#ifdef _UNIX
  struct tms usage;
  static int clock_ticks = sysconf(_SC_CLK_TCK);
  times(&usage);
  df=((double)usage.tms_utime+(double)usage.tms_stime)/clock_ticks;
#elif LINUX
  struct tms usage;
  static int clock_ticks = sysconf(_SC_CLK_TCK);
  times(&usage);
  df=((double)usage.tms_utime+(double)usage.tms_stime)/clock_ticks;
#elif MACOSX
  struct tms usage;
  static int clock_ticks = sysconf(_SC_CLK_TCK);
  times(&usage);
  df=((double)usage.tms_utime+(double)usage.tms_stime)/clock_ticks;
#endif

  return df;
}

unsigned long int Mistral::get_memory() {
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

void Mistral::get_command_line(const char** int_ident, 
			     int*         int_param,
			     int          int_nb,
			     const char** str_ident,
			     const char** str_param,
			     int          str_nb,
			     char**       argv,
			     int          argc)
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


int Mistral::log2_( const unsigned int v ) {
//   union {float f; unsigned int i; } t;
//   unsigned int b = v & -v;
  
//   t.f = (float)b; // cast the least significant bit in v to a float
//   b = t.i >> 23;
//   return b - 0x7f;
// }
  if( !v ) return NOVAL;
  int exponent = -1;
  while( (v >> (++exponent)) > 1 );
  return exponent;
}

std::string Mistral::int2str(const int x) {

 //  unsigned n = x;

  //char buffer[2]; // long enough for 32-bit UINT_MAX + NUL character
 // char *p = buffer + sizeof(buffer);

 // *--p = '\0';
 // do {
 //   *--p = '0' + n % 10;
 //   n /= 10;
 // } while (n);

 // return p;

  std::ostringstream oss;
  
  oss << x;

  //buffer[1] = '\0';
  std::string stb = oss.str();
  
  //const char* buf = stb.c_str();

  return stb;

}


// int __modulo_fct__(const int x, const int m) {
//   int mod = x%m;
//   if(mod && (mod<0) != (m<0))  mod += m;
//   return mod;
// }


// std::string Mistral::toString(const int x)
// {
//   std::ostringstream o;
//   o << x;
//   return o.str();
// }

// // std::string Mistral::toString(const V)
// // {
// //   return (std::string)(*x);
// // }

// std::string Mistral::toString(const Mistral::IntStack& x) {
//   return x.getString();
// }

// std::string Mistral::toString(const Mistral::Queue& x) {
//   return x.getString();
// }

// std::string Mistral::toString(const Mistral::MultiSet& x) {
//   return x.getString();
// }

// std::string Mistral::toString(const Mistral::ConstraintTrigger& x) {
//   return x.getString();
// }

// // std::string Mistral::toString(const Mistral::Goal* x) {
// //   return x->getString();
// // }

// std::string Mistral::toString(const Mistral::Constraint* x) {
//   return x->getString();
// }

// std::string Mistral::toString(const Mistral::SolverStatistics& x) {
//   return x.getString();
// }

// // std::string Mistral::toString(const Mistral::VariableInt* x) {
// //   return x->getString();
// // }


