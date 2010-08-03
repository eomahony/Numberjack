
#ifndef WALKSRC_H
#define WALKSRC_H

/* walksat */
/* Maintained by Henry Kautz <kautz@cs.rochester.edu> */

#define VERSION "walksat-v47" /* 15 Jan 2010 - fixes printing long long ints */

/************************************/
/* Compilation flags                */
/************************************/

/* If the constant DYNAMIC is set, then dynamic arrays are used instead of static arrays.
   This allows very large or very small problems to be handled, without
   having to adjust the constants MAXATOM and/or MAXCLAUSE.  However, on some
   architectures the compiled program is about 15% slower, because not
   as many optimizations can be performed. */

#define DYNAMIC 1

/* Alex Fukunaga watched-literal speedup is enabled if _WATCH is defined.
   Is is not compatible with the DYNAMIC option!  */
/* #define _WATCH */


/********************************************************************/
/* Please select only one of following flags: UNIX, OSX, ANSI or NT */
/* Description:                                                     */
/*    UNIX:  uses some non-POSIX unix time routines, compiles       */
/*           under Linux and OSX                                    */
/*    ANSI:  use POSIX time routines, can compile under all         */
/*           unix and NT, time accuracy is 1sec                     */ 
/*    NT:    use NT/Win32 multimedia routines, compile under        */
/*           NT only, time accuracy is 1ms                          */
/*           Uses NT "timeGetTime" function:                        */
/*             Header: Decl. in Mmsystem.h; include Windows.h.      */
/*             Library: Use Winmm.lib.                              */
/********************************************************************/

#define UNIX 1
#define ANSI 0
#define NT 0

/* Define BIGINT to be the type used for the "cutoff" variable.
   Under gcc "long long int" gives a 64 bit integer.
   Program will still function using a 32-bit value, but it
   limit size of cutoffs that can be specified. */

#if UNIX
#define BIGINT long long int
#define BIGINTSTR "lli"
#endif
#if NT
#define BIGINT __int64
#define BIGINTSTR "I64d"
#endif
#if ANSI
#define BIGINT long int
#define BIGINTSTR "li"
#endif


/************************************/
/* Standard includes                */
/************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <limits.h>
#include <signal.h>

#if UNIX || OSX
  #include <sys/times.h>
  #include <sys/time.h>
  #include <unistd.h>
#endif
#if NT
  #include <time.h>
  #include <windows.h>
  #include <mmsystem.h>
#endif
#if ANSI
  #include <sys/time.h>
#endif

#if ANSI || NT
  #define random() rand()
  #define srandom(seed) srand(seed)
#endif

#ifdef _WATCH 
#include <assert.h>
#endif


/************************************/
/* Constant parameters              */
/************************************/

#define MAXATOM 150000		/* maximum possible number of atoms */

#ifdef DYNAMIC
  #define STOREBLOCK 2000000	/* size of block to malloc each time */
#else
  #define STOREBLOCK 2000000	/* size of block to malloc each time */
  #define MAXCLAUSE 500000	/* maximum possible number of clauses */
#endif

#define TRUE 1
#define FALSE 0

#define MAXLENGTH 500           /* maximum number of literals which can be in any clause */

/************************************/
/* Internal constants               */
/************************************/


enum heuristics { RANDOM, BEST, TABU, NOVELTY, RNOVELTY, NOVELTY_PLUS, RNOVELTY_PLUS};

#define NOVALUE -1
#define INIT_PARTIAL 1
#define HISTMAX 101		/* length of histogram of tail */

#define Var(CLAUSE, POSITION) (ABS(clause[CLAUSE][POSITION]))

static int scratch;
#define ABS(x) ((scratch=(x))>0?(scratch):(-scratch))

#define BIG 100000000


void handle_interrupt(int sig);


/************************************/
/* Main data structures             */
/************************************/

/* Atoms start at 1 */
/* Not a is recorded as -1 * a */
/* One dimensional arrays are statically allocated. */
/* Two dimensional arrays are dynamically allocated in */
/* the second dimension only.  */

class WalksatAlgorithm {

public: 

  //const char* heuristic_names[7] = {"random", "best", "tabu", "novelty", "rnovelty", "novelty+", "rnovelty+"};

  
  //static int (*pickcode[])(void) =
  //{pickrandom, pickbest, picktabu, picknovelty, pickrnovelty,
  //picknoveltyplus, pickrnoveltyplus};


  /* No longer using constant CLK_TCK - it is deprecated! */
  /* Instead: */
  long ticks_per_second;
  
  int numatom;
  int numclause;
  int numliterals;
  
#ifdef DYNAMIC
  int ** clause;			/* clauses to be satisfied */
  /* indexed as clause[clause_num][literal_num] */
  int * size;			/* length of each clause */
  int * falsified;		/* clauses which are false */
  int * lowfalse;
  int * wherefalse;		/* where each clause is listed in falsified */
  int * numtruelit;		/* number of true literals in each clause */
#else
  int * clause[MAXCLAUSE];	/* clauses to be satisfied */
				/* indexed as clause[clause_num][literal_num] */
  int size[MAXCLAUSE];		/* length of each clause */
  int falsified[MAXCLAUSE];	/* clauses which are false */
  int lowfalse[MAXCLAUSE];
  int wherefalse[MAXCLAUSE];	/* where each clause is listed in falsified */
  int numtruelit[MAXCLAUSE];	/* number of true literals in each clause */
#endif
  
  int *occurence[2*MAXATOM+1];	/* where each literal occurs */
				/* indexed as occurence[literal+MAXATOM][occurence_num] */
  
  int numoccurence[2*MAXATOM+1];	/* number of times each literal occurs */
  
  
  int atom[MAXATOM+1];		/* value of each atom */ 
  int lowatom[MAXATOM+1];
  int solution[MAXATOM+1];
  
  int changed[MAXATOM+1];		/* step at which atom was last flipped */
  
  int breakcount[MAXATOM+1];	/* number of clauses that become unsat if var if flipped */
  int makecount[MAXATOM+1];	/* number of clauses that become sat if var if flipped */
  
  int numfalse;			/* number of false clauses */

#ifdef _WATCH 
  int watch1[MAXCLAUSE];
  int watch2[MAXCLAUSE];
#endif

  
  /************************************/
  /* Global flags and parameters      */
  /************************************/
  
  int status_flag;		/* value returned from main procedure */
  int abort_flag;
  
  int heuristic;		/* heuristic to be used */
  int numerator;	/* make random flip with numerator/denominator frequency */
  int denominator;		
  int tabu_length;		/* length of tabu list */
  
  int wp_numerator;	/* walk probability numerator/denominator */
  int wp_denominator;		

  BIGINT numflip;		/* number of changes so far */
  BIGINT numnullflip;		/*  number of times a clause was picked, but no  */
				/*  variable from it was flipped  */
  int numrun;
  BIGINT cutoff;
  BIGINT base_cutoff;
  int target;
  int numtry;			/* total attempts at solutions */
  int numsol;		/* stop after this many tries succeeds */
  int superlinear;
  
  int makeflag;		/* set to true by heuristics that require the make values to be calculated */
  

  /* Printing options */
  int verbosity;
  int printonlysol;
  int printsolcnf;
  int printfalse;
  int printlow;
  int printhist;
  int printtrace;
  int trace_assign;
  
  char outfile[1024];
  
  /* Initialization options */
  
  char initfile[100];
  int initoptions;
  
  /* Randomization */
  
  unsigned int seed;  /* Sometimes defined as an unsigned long int */
  
#if UNIX || OSX
  struct timeval tv;
  struct timezone tzp;
#endif
#if NT
  DWORD win_time;     /* elapsed time in ms, since windows boot up */
#endif

  /* Statistics */
  
  double time_cutoff;
  double expertime;
  BIGINT flips_this_solution;
  long int lowbad;		/* lowest number of bad clauses during try */
  BIGINT totalflip;		/* total number of flips in all tries so far */
  BIGINT totalsuccessflip;	/* total number of flips in all tries which succeeded so far */
  int numsuccesstry;		/* total found solutions */
  
  /* Histogram of tail */
  
  long int tailhist[HISTMAX];	/* histogram of num unsat in tail of run */
  long histtotal;
  int tail;
  int tail_start_flip;
  
  
  BIGINT x;
  BIGINT integer_sum_x;
  double sum_x;
  double sum_x_squared;
  double mean_x;
  double second_moment_x;
  double variance_x;
  double std_dev_x;
  double std_error_mean_x;
  double seconds_per_flip;
  int r;
  int sum_r;
  double sum_r_squared;
  double mean_r;
  double variance_r;
  double std_dev_r;
  double std_error_mean_r;
  
  double avgfalse;
  double sumfalse;
  double sumfalse_squared;
  double second_moment_avgfalse, variance_avgfalse, std_dev_avgfalse, ratio_avgfalse;
  double f;
  double sample_size;

  double sum_avgfalse;
  double sum_std_dev_avgfalse;
  double mean_avgfalse;
  double mean_std_dev_avgfalse;
  int number_sampled_runs;
  double ratio_mean_avgfalse;
  
  double suc_sum_avgfalse;
  double suc_sum_std_dev_avgfalse;
  double suc_mean_avgfalse;
  double suc_mean_std_dev_avgfalse;
  int suc_number_sampled_runs;
  double suc_ratio_mean_avgfalse;
  
  double nonsuc_sum_avgfalse;
  double nonsuc_sum_std_dev_avgfalse;
  double nonsuc_mean_avgfalse;
  double nonsuc_mean_std_dev_avgfalse;
  int nonsuc_number_sampled_runs;
  double nonsuc_ratio_mean_avgfalse;
  
  /* Hamming calcualations */
  char hamming_target_file[512];
  char hamming_data_file[512];
  int hamming_sample_freq;
  int hamming_flag;
  int hamming_distance;
  int hamming_target[MAXATOM+1];
  FILE * hamming_fp;
  
  /* Noise level */
  int samplefreq;


  const char* heuristic_names[7]; 

  WalksatAlgorithm();
  ~WalksatAlgorithm() {}


  /************************************/
  /* Forward declarations             */
  /************************************/

  void read_hamming_file(char initfile[]);
  void open_hamming_data(char initfile[]);
  int calc_hamming_dist(int atom[], int hamming_target[], int numatom);
  
  void parse_parameters(int argc,char *argv[]);
  void print_parameters(int argc, char * argv[]);
  
  int pickrandom(void);
  int pickbest(void);
  int picktabu(void);
  int picknovelty(void);
  int pickrnovelty(void);
  int picknoveltyplus(void);
  int pickrnoveltyplus(void);
  
  double elapsed_seconds(void);
  int countunsat(void);
  void scanone(int argc, char *argv[], int i, int *varptr);
  void scanoneu(int argc, char *argv[], int i, unsigned int *varptr);
  void scanonell(int argc, char *argv[], int i, BIGINT *varptr);
  void init(char initfile[], int initoptions);
  void initprob(void); 
  void flipatom(int toflip);
  
  void print_false_clauses(long int lowbad);
  void save_false_clauses(long int lowbad);
  void print_low_assign(long int lowbad);
  void save_low_assign(void);
  void save_solution(void);
  
  void print_current_assign(void);

  long super(int i);
  
  void print_sol_file(char * filename);
  void print_statistics_header(void);
  void initialize_statistics(void);
  void update_statistics_start_try(void);
  void print_statistics_start_flip(void);
  void update_and_print_statistics_end_try(void);
  void update_statistics_end_flip(void);
  void print_statistics_final(void);
  void print_sol_cnf(void);
  
  int walk_solve();

};

#endif

