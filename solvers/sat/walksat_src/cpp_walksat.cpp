
#include "cpp_walksat.hpp"
#include <iostream>

WalksatAlgorithm *current_solver;
void handle_interrupt(int sig)
{
    if (current_solver->abort_flag) exit(-1);
    current_solver->abort_flag = TRUE;
}



/************************************/
/* Main                             */
/************************************/


WalksatAlgorithm::WalksatAlgorithm() {
  status_flag = 0;		/* value returned from main procedure */
  //abort_flag;
  
  heuristic = BEST;		/* heuristic to be used */
  numerator = NOVALUE;	/* make random flip with numerator/denominator frequency */
  denominator = 100;		
  //tabu_length;		/* length of tabu list */
  
  wp_numerator = NOVALUE;	/* walk probability numerator/denominator */
  wp_denominator = 100;		

  //numflip;		/* number of changes so far */
  //numnullflip;		/*  number of times a clause was picked, but no  */
				/*  variable from it was flipped  */
  numrun = 10000;
  cutoff = 250000;
  base_cutoff = 100000;
  target = 0;
  numtry = 0;			/* total attempts at solutions */
  numsol = 1;		/* stop after this many tries succeeds */
  superlinear = TRUE;
  
  makeflag = FALSE;		/* set to true by heuristics that require the make values to be calculated */
  

  /* Printing options */
  verbosity = 0;
  printonlysol = FALSE;
  printsolcnf = FALSE;
  printfalse = FALSE;
  printlow = FALSE;
  printhist = FALSE;
  printtrace = FALSE;
  trace_assign = FALSE;
  
  outfile[0] = 0;
  
  /* Initialization options */
  
  initfile[0] = 0;
  initoptions = FALSE;
  
  /* Randomization */
  
  //unsigned int seed;  /* Sometimes defined as an unsigned long int */
  
// #if UNIX || OSX
//   struct timeval tv;
//   struct timezone tzp;
// #endif
// #if NT
//   DWORD win_time;     /* elapsed time in ms, since windows boot up */
// #endif

  /* Statistics */
  
  time_cutoff = 0.0;
  expertime = 0.0;
  //flips_this_solution;
  //lowbad;		/* lowest number of bad clauses during try */
  totalflip = 0;		/* total number of flips in all tries so far */
  totalsuccessflip = 0;	/* total number of flips in all tries which succeeded so far */
  numsuccesstry = 0;		/* total found solutions */
  
  /* Histogram of tail */
  
  //tailhist;	/* histogram of num unsat in tail of run */
  //histtotal;
  tail = 3;
  //tail_start_flip;
  
  
  //x;
  integer_sum_x = 0;
  sum_x = 0.0;
  sum_x_squared = 0.0;
  //mean_x;
  //second_moment_x;
  //variance_x;
  //std_dev_x;
  //std_error_mean_x;
  //seconds_per_flip;
  //r;
  sum_r = 0;
  sum_r_squared = 0.0;
  //mean_r;
  //variance_r;
  //std_dev_r;
  //std_error_mean_r;
  
  //avgfalse;
  //sumfalse;
  //sumfalse_squared;
  //second_moment_avgfalse;
  //variance_avgfalse;
  //std_dev_avgfalse;
  //ratio_avgfalse;
  //f;
  //sample_size;

  sum_avgfalse = 0.0;
  sum_std_dev_avgfalse = 0.0;
  //mean_avgfalse;
  //mean_std_dev_avgfalse;
  number_sampled_runs = 0;
  //ratio_mean_avgfalse;
  
  suc_sum_avgfalse = 0.0;
  suc_sum_std_dev_avgfalse = 0.0;
  //suc_mean_avgfalse;
  //suc_mean_std_dev_avgfalse;
  suc_number_sampled_runs = 0;
  //suc_ratio_mean_avgfalse;
  
  nonsuc_sum_avgfalse = 0.0;
  nonsuc_sum_std_dev_avgfalse = 0.0;
  //nonsuc_mean_avgfalse;
  //nonsuc_mean_std_dev_avgfalse;
  nonsuc_number_sampled_runs = 0;
  //nonsuc_ratio_mean_avgfalse;
  
  /* Hamming calcualations */
  
  hamming_target_file[0] = 0;
  hamming_data_file[0] = 0;
  //hamming_sample_freq;
  hamming_flag = FALSE;
  //hamming_distance;
  //hamming_target;
  //hamming_fp;
  
  /* Noise level */
  samplefreq = 1;


#if UNIX || OSX
  ticks_per_second = sysconf(_SC_CLK_TCK);
  gettimeofday(&tv,&tzp);
  seed = (unsigned int)((( tv.tv_sec & 0177 ) * 1000000) + tv.tv_usec);
#endif
#if NT
  seed = (unsigned int)(timeGetTime());
#endif
#if ANSI
  seed = (unsigned int)(time());
#endif
  
  //  heuristic_names = 
  heuristic_names[0] = (char*)"random";
  heuristic_names[1] = (char*)"best";
  heuristic_names[2] = (char*)"tabu";
  heuristic_names[3] = (char*)"novelty";
  heuristic_names[4] = (char*)"rnovelty";
  heuristic_names[5] = (char*)"novelty+";
  heuristic_names[6] = (char*)"rnovelty+";
}

int WalksatAlgorithm::walk_solve() {
  current_solver = this;
  signal(SIGINT, handle_interrupt);
  abort_flag = FALSE;
  (void) elapsed_seconds();
  
  while (! abort_flag && 
	 numsuccesstry < numsol && numtry < numrun &&
	 (time_cutoff == 0.0 ||
	  time_cutoff > expertime)) {

    numtry++;
    init(initfile, initoptions);
    update_statistics_start_try();
    numflip = 0;
    
    if (superlinear) cutoff = base_cutoff * super(numtry);
    
    while((numfalse > target) && (numflip < cutoff)) {
      if(verbosity) print_statistics_start_flip();
      numflip++;
      //flipatom((pickcode[heuristic])());
      flipatom(pickbest());
      update_statistics_end_flip();
    }
    update_and_print_statistics_end_try();

    expertime += elapsed_seconds();
  }
  //expertime = elapsed_seconds();
  if(verbosity) print_statistics_final();
  return status_flag;
}



void WalksatAlgorithm::parse_parameters(int argc,char *argv[])
{
    int i;
    int temp;

    for (i=1;i < argc;i++)
    {
	if (strcmp(argv[i],"-seed") == 0){
	    scanone(argc,argv,++i,&temp);
            seed = (unsigned int)temp;
        }
	else if (strcmp(argv[i],"-out") == 0 && i<argc-1)
	    strcpy(outfile, argv[++i]);
	else if (strcmp(argv[i],"-hist") == 0)
	    printhist = TRUE;
	else if (strcmp(argv[i],"-status") == 0)
	    status_flag = 1;
	else if (strcmp(argv[i],"-cutoff") == 0)
	    scanonell(argc,argv,++i,&cutoff);
	else if (strcmp(argv[i],"-random") == 0)
	    heuristic = RANDOM;
       	else if (strcmp(argv[i],"-novelty") == 0){
	    heuristic = NOVELTY;
	    makeflag = TRUE;
	}
	else if (strcmp(argv[i],"-rnovelty") == 0){
	    heuristic = RNOVELTY;
	    makeflag = TRUE;
	}
      else if (strcmp(argv[i],"-novelty+") == 0){
	heuristic = NOVELTY_PLUS;
	makeflag = TRUE;
      }
      else if (strcmp(argv[i],"-rnovelty+") == 0){
	heuristic = RNOVELTY_PLUS;
	makeflag = TRUE;
      }
	else if (strcmp(argv[i],"-best") == 0)
	    heuristic = BEST;
	else if (strcmp(argv[i],"-noise") == 0){
	    scanone(argc,argv,++i,&numerator);
	    if (i < argc-1 && sscanf(argv[i+1],"%i",&temp)==1){
		denominator = temp;
		i++;
	    }
	}
        else if (strcmp(argv[i],"-wp") == 0){
            scanone(argc,argv,++i,&wp_numerator);
            if (i < argc-1 && sscanf(argv[i+1],"%i",&temp)==1){
                wp_denominator = temp;
                i++;
            }
        }
	else if (strcmp(argv[i],"-init") == 0  && i < argc-1)
	    sscanf(argv[++i], " %s", initfile);
	else if (strcmp(argv[i],"-hamming") == 0  && i < argc-3){
	    sscanf(argv[++i], " %s", hamming_target_file);
	    sscanf(argv[++i], " %s", hamming_data_file);
	    sscanf(argv[++i], " %i", &hamming_sample_freq);
	    hamming_flag = TRUE;
	    numrun = 1;
	}
	else if (strcmp(argv[i],"-partial") == 0)
	    initoptions = INIT_PARTIAL;
	else if (strcmp(argv[i],"-super") == 0)
	    superlinear = TRUE;
	else if (strcmp(argv[i],"-tries") == 0 || strcmp(argv[i],"-restart") == 0)
	    scanone(argc,argv,++i,&numrun);
	else if (strcmp(argv[i],"-target") == 0)
	    scanone(argc,argv,++i,&target);
	else if (strcmp(argv[i],"-tail") == 0)
	    scanone(argc,argv,++i,&tail);
	else if (strcmp(argv[i],"-sample") == 0)
	    scanone(argc,argv,++i,&samplefreq);
	else if (strcmp(argv[i],"-tabu") == 0)
	{
	    scanone(argc,argv,++i,&tabu_length);
	    heuristic = TABU;
	}
	else if (strcmp(argv[i],"-low") == 0)
	    printlow = TRUE;
	else if (strcmp(argv[i],"-sol") == 0)
	{
	    printonlysol = TRUE;
	    printlow = TRUE;
	}
	else if (strcmp(argv[i],"-solcnf") == 0)
	{
	    printsolcnf = TRUE;
	    if (numsol == NOVALUE) numsol = 1;
	}
	else if (strcmp(argv[i],"-bad") == 0)
	    printfalse = TRUE;
	else if (strcmp(argv[i],"-numsol") == 0)
	    scanone(argc,argv,++i,&numsol);
	else if (strcmp(argv[i],"-trace") == 0)
	    scanone(argc,argv,++i,&printtrace);
	else if (strcmp(argv[i],"-assign") == 0){
	    scanone(argc,argv,++i,&printtrace);
	    trace_assign = TRUE;
	}
	else 
	{
	    fprintf(stderr, "General parameters:\n");
	    fprintf(stderr, "  -seed N -cutoff N -restart N\n");
	    fprintf(stderr, "  -numsol N = stop after finding N solutions\n");
	    fprintf(stderr, "  -super = use the Luby series for the cutoff values\n");
	    fprintf(stderr, "  -init FILE = set vars not included in FILE to false\n");
	    fprintf(stderr, "  -partial FILE = set vars not included in FILE randomly\n");
	    fprintf(stderr, "  -status = return fail status if solution not found\n");
	    fprintf(stderr, "  -target N = succeed if N or fewer clauses unsatisfied\n");
	    fprintf(stderr, "Heuristics:\n");
	    fprintf(stderr, "  -random -best -tabu N -novelty -rnovelty\n");
	    fprintf(stderr, "  -noise N or -noise N M (default M = 100)\n");
            fprintf(stderr, "  -novelty+ -rnovelty+\n");
            fprintf(stderr, "  -wp N or -wp N M (default M = 100) = cycle breaking for (r)novelty+\n");
	    fprintf(stderr, "Printing:\n");
	    fprintf(stderr, "  -out FILE = print solution as a list of literals to FILE\n");
	    fprintf(stderr, "  -trace N = print statistics every N flips\n");
	    fprintf(stderr, "  -assign N = print assignments at flip N, 2N, ...\n");
	    fprintf(stderr, "  -sol = print satisfying assignments to stdout\n");
	    fprintf(stderr, "  -solcnf = print sat assign to stdout in DIMACS format, and exit\n");
	    fprintf(stderr, "  -low = print lowest assignment each try\n");
	    fprintf(stderr, "  -bad = print unsat clauses each try\n");
	    fprintf(stderr, "  -hist = print histogram of tail\n");
	    fprintf(stderr, "  -tail N = assume tail begins at nvars*N\n");
	    fprintf(stderr, "  -sample N = sample noise level every N flips\n");
	    fprintf(stderr, "  -hamming TARGET_FILE DATA_FILE SAMPLE_FREQUENCY\n");
	    exit(-1);
	}
    }
    base_cutoff = cutoff;
    if (numsol==NOVALUE || numsol>numrun) numsol = numrun;
    if (numerator==NOVALUE){
	switch(heuristic) {
	case BEST:
	case NOVELTY:
	case RNOVELTY:
        case NOVELTY_PLUS:
        case RNOVELTY_PLUS:
	    numerator = 50;
	    break;
	default:
	    numerator = 0;
	    break;
	}
    }
  if (wp_numerator==NOVALUE){
    switch(heuristic) {
    case NOVELTY_PLUS:
    case RNOVELTY_PLUS:
      wp_numerator = 1;
      break;
    default:
      wp_numerator = 0;
      break;
    }
  }
}


void WalksatAlgorithm::print_parameters(int argc, char * argv[])
{
  int i;

  printf("%s\n", VERSION);
  printf("command line =");
  for (i=0;i < argc;i++){
    printf(" %s", argv[i]);
  }
  printf("\n");
  printf("seed = %u\n",seed);
  printf("cutoff = %Li\n",cutoff);
  printf("tries = %i\n",numrun);

  printf("heuristic = ");
  switch(heuristic)
    {
    case TABU:
      printf("tabu %d", tabu_length);
      break;
    default:
      printf("%s", heuristic_names[heuristic]);
      break;
    }
  if (numerator>0){
    printf(", noise %d / %d", numerator, denominator);
  }
  if (wp_numerator>0){
    printf(", wp %d / %d", wp_numerator, wp_denominator);
  }
  printf("\n");
}

void WalksatAlgorithm::print_statistics_header(void)
{
    printf("numatom = %i, numclause = %i, numliterals = %i\n",numatom,numclause,numliterals);
    printf("wff read in\n\n");
    printf("    lowest     final       avg     noise     noise     total                 avg        mean        mean\n");
    printf("    #unsat    #unsat     noise   std dev     ratio     flips              length       flips       flips\n");
    printf("      this      this      this      this      this      this   success   success       until         std\n");
    printf("       try       try       try       try       try       try      rate     tries      assign         dev\n\n");

    fflush(stdout);

}

void WalksatAlgorithm::initialize_statistics(void)
{
    x = 0; r = 0;
    if (hamming_flag) {
	read_hamming_file(hamming_target_file);
	open_hamming_data(hamming_data_file);
    }
    tail_start_flip = tail * numatom;
    if(verbosity) printf("tail starts after flip = %i\n", tail_start_flip);
    numnullflip = 0;


}

void WalksatAlgorithm::update_statistics_start_try(void)
{
    int i;

    lowbad = numfalse;

    sample_size = 0;
    sumfalse = 0.0;
    sumfalse_squared = 0.0;

    for (i=0; i<HISTMAX; i++)
	tailhist[i] = 0;
    if (tail_start_flip == 0){
	tailhist[numfalse < HISTMAX ? numfalse : HISTMAX - 1] ++;
    }
	      
    if (printfalse) save_false_clauses(lowbad);
    if (printlow) save_low_assign();
}

void WalksatAlgorithm::print_statistics_start_flip(void)
{
    if (printtrace && (numflip % printtrace == 0)){
	printf(" %9li %9i                     %9" BIGINTSTR "\n", lowbad,numfalse,numflip);
	if (trace_assign)
	    print_current_assign();
	fflush(stdout);
    }
}


void WalksatAlgorithm::update_and_print_statistics_end_try(void)
{
    int i;
    int j;

    totalflip += numflip;
    x += numflip;
    r ++;

    if (sample_size > 0){
	avgfalse = sumfalse/sample_size;
	second_moment_avgfalse = sumfalse_squared / sample_size;
	variance_avgfalse = second_moment_avgfalse - (avgfalse * avgfalse);
	if (sample_size > 1) { variance_avgfalse = (variance_avgfalse * sample_size)/(sample_size - 1); }
	std_dev_avgfalse = sqrt(variance_avgfalse);

	ratio_avgfalse = avgfalse / std_dev_avgfalse;

	sum_avgfalse += avgfalse;
	sum_std_dev_avgfalse += std_dev_avgfalse;
	number_sampled_runs += 1;

	if (numfalse <= target){
	    suc_number_sampled_runs += 1;
	    suc_sum_avgfalse += avgfalse;
	    suc_sum_std_dev_avgfalse += std_dev_avgfalse;
	}
	else {
	    nonsuc_number_sampled_runs += 1;
	    nonsuc_sum_avgfalse += avgfalse;
	    nonsuc_sum_std_dev_avgfalse += std_dev_avgfalse;
	}
    }
    else{
	avgfalse = 0;
	variance_avgfalse = 0;
	std_dev_avgfalse = 0;
	ratio_avgfalse = 0;
    }

    if(numfalse <= target){

	status_flag = 0;

	save_solution();
	numsuccesstry++;

	totalsuccessflip += numflip;
	integer_sum_x += x;
	sum_x = (double) integer_sum_x;
	sum_x_squared += ((double)x)*((double)x);
	mean_x = sum_x / numsuccesstry;
	if (numsuccesstry > 1){
	    second_moment_x = sum_x_squared / numsuccesstry;
	    variance_x = second_moment_x - (mean_x * mean_x);
	    /* Adjustment for small small sample size */
	    variance_x = (variance_x * numsuccesstry)/(numsuccesstry - 1);
	    std_dev_x = sqrt(variance_x);
	    std_error_mean_x = std_dev_x / sqrt((double)numsuccesstry);
	}
	sum_r += r;
	mean_r = ((double)sum_r)/numsuccesstry;
	sum_r_squared += ((double)r)*((double)r);

	x = 0;
	r = 0;
    }

    if(verbosity) {  
      printf(" %9li %9i %9.2f %9.2f %9.2f %9" BIGINTSTR " %9i",
	     lowbad,numfalse,avgfalse, std_dev_avgfalse,ratio_avgfalse,numflip, (numsuccesstry*100)/numtry);
      if (numsuccesstry > 0){
	printf(" %9" BIGINTSTR , totalsuccessflip/numsuccesstry);
	printf(" %11.2f", mean_x);
	if (numsuccesstry > 1){
	  printf(" %11.2f", std_dev_x);
	}
      }
      printf("\n");
      
      if (printhist){
	printf("histogram: ");
	for (j=HISTMAX-1; tailhist[j] == 0; j--);
	for (i=0; i<=j; i++){
	  printf(" %li(%i)", tailhist[i], i);
	  if ((i+1) % 10 == 0) printf("\n           ");
	}
	if (j==HISTMAX-1) printf(" +++");
	printf("\n");
      }
      
      if (numfalse>0 && printfalse)
	print_false_clauses(lowbad);
      if (printlow && (!printonlysol || numfalse >= target))
	print_low_assign(lowbad);
      
      if(numfalse == 0 && countunsat() != 0){
	fprintf(stderr, "Program error, verification of solution fails!\n");
	exit(-1);
      }
      
      fflush(stdout);
    }
}

void WalksatAlgorithm::update_statistics_end_flip(void)
{
    if (numfalse < lowbad){
	lowbad = numfalse;
	if (printfalse) save_false_clauses(lowbad);
	if (printlow) save_low_assign();
    }
    if (numflip >= tail_start_flip){
	tailhist[(numfalse < HISTMAX) ? numfalse : (HISTMAX - 1)] ++;
	if ((numflip % samplefreq) == 0){
	  sumfalse += numfalse;
	  sumfalse_squared += numfalse * numfalse;
	  sample_size ++;
	}
    }
}

void WalksatAlgorithm::print_statistics_final(void)
{
    seconds_per_flip = expertime / totalflip;
    printf("\ntotal elapsed seconds = %f\n", expertime);
    printf("average flips per second = %ld\n", (long)(totalflip/expertime));
    if (heuristic == TABU)
	printf("proportion null flips = %f\n", ((double)numnullflip)/totalflip);
    printf("number solutions found = %d\n", numsuccesstry);
    printf("final success rate = %f\n", ((double)numsuccesstry * 100.0)/numtry);
    printf("average length successful tries = %" BIGINTSTR "\n", numsuccesstry ? (totalsuccessflip/numsuccesstry) : 0);
    if (numsuccesstry > 0)
	{
	    printf("average flips per assign (over all runs) = %f\n", ((double)totalflip)/numsuccesstry);
	    printf("average seconds per assign (over all runs) = %f\n", (((double)totalflip)/numsuccesstry)*seconds_per_flip);
	    printf("mean flips until assign = %f\n", mean_x);
	    if (numsuccesstry>1){
		printf("  variance = %f\n", variance_x);
		printf("  standard deviation = %f\n", std_dev_x);
		printf("  standard error of mean = %f\n", std_error_mean_x);
	    }
	    printf("mean seconds until assign = %f\n", mean_x * seconds_per_flip);
	    if (numsuccesstry>1){
		printf("  variance = %f\n", variance_x * seconds_per_flip * seconds_per_flip);
		printf("  standard deviation = %f\n", std_dev_x * seconds_per_flip);
		printf("  standard error of mean = %f\n", std_error_mean_x * seconds_per_flip);
	    }
	    printf("mean restarts until assign = %f\n", mean_r);
	    if (numsuccesstry>1){
		variance_r = (sum_r_squared / numsuccesstry) - (mean_r * mean_r);
		if (numsuccesstry > 1) variance_r = (variance_r * numsuccesstry)/(numsuccesstry - 1);	   
		std_dev_r = sqrt(variance_r);
		std_error_mean_r = std_dev_r / sqrt((double)numsuccesstry);
		printf("  variance = %f\n", variance_r);
		printf("  standard deviation = %f\n", std_dev_r);
		printf("  standard error of mean = %f\n", std_error_mean_r);
	    }
	}

    if (number_sampled_runs){
      mean_avgfalse = sum_avgfalse / number_sampled_runs;
      mean_std_dev_avgfalse = sum_std_dev_avgfalse / number_sampled_runs;
      ratio_mean_avgfalse = mean_avgfalse / mean_std_dev_avgfalse;

      if (suc_number_sampled_runs){
	  suc_mean_avgfalse = suc_sum_avgfalse / suc_number_sampled_runs;
	  suc_mean_std_dev_avgfalse = suc_sum_std_dev_avgfalse / suc_number_sampled_runs;
	  suc_ratio_mean_avgfalse = suc_mean_avgfalse / suc_mean_std_dev_avgfalse;
      }
      else {
	  suc_mean_avgfalse = 0;
	  suc_mean_std_dev_avgfalse = 0;
	  suc_ratio_mean_avgfalse = 0;
      }

      if (nonsuc_number_sampled_runs){
	  nonsuc_mean_avgfalse = nonsuc_sum_avgfalse / nonsuc_number_sampled_runs;
	  nonsuc_mean_std_dev_avgfalse = nonsuc_sum_std_dev_avgfalse / nonsuc_number_sampled_runs;
	  nonsuc_ratio_mean_avgfalse = nonsuc_mean_avgfalse / nonsuc_mean_std_dev_avgfalse;
      }
      else {
	  nonsuc_mean_avgfalse = 0;
	  nonsuc_mean_std_dev_avgfalse = 0;
	  nonsuc_ratio_mean_avgfalse = 0;
      }

      printf("final noise level statistics\n");
      printf("    statistics over all runs:\n");
      printf("      overall mean average noise level = %f\n", mean_avgfalse);
      printf("      overall mean noise std deviation = %f\n", mean_std_dev_avgfalse);
      printf("      overall ratio mean noise to mean std dev = %f\n", ratio_mean_avgfalse);
      printf("    statistics on successful runs:\n");
      printf("      successful mean average noise level = %f\n", suc_mean_avgfalse);
      printf("      successful mean noise std deviation = %f\n", suc_mean_std_dev_avgfalse);
      printf("      successful ratio mean noise to mean std dev = %f\n", suc_ratio_mean_avgfalse);
      printf("    statistics on nonsuccessful runs:\n");
      printf("      nonsuccessful mean average noise level = %f\n", nonsuc_mean_avgfalse);
      printf("      nonsuccessful mean noise std deviation = %f\n", nonsuc_mean_std_dev_avgfalse);
      printf("      nonsuccessful ratio mean noise to mean std dev = %f\n", nonsuc_ratio_mean_avgfalse);
    }

    if (hamming_flag){
	fclose(hamming_fp);
	printf("Final distance to hamming target = %i\n", calc_hamming_dist(atom, hamming_target, numatom));
	printf("Hamming distance data stored in %s\n", hamming_data_file);
    }

    if (numsuccesstry > 0){
	printf("ASSIGNMENT FOUND\n");
	if (printsolcnf == TRUE) print_sol_cnf();
	if (outfile[0]) print_sol_file(outfile);
    }
    else
	printf("ASSIGNMENT NOT FOUND\n");

}


long WalksatAlgorithm::super(int i)
{
    long power;
    int k;

    if (i<=0){
	fprintf(stderr, "bad argument super(%d)\n", i);
	exit(1);
    }
    /* let 2^k be the least power of 2 >= (i+1) */
    k = 1;
    power = 2;
    while (power < (i+1)){
	k += 1;
	power *= 2;
    }
    if (power == (i+1)) return (power/2);
    return (super(i - (power/2) + 1));
}

void WalksatAlgorithm::scanone(int argc, char *argv[], int i, int *varptr)
{
    if (i>=argc || sscanf(argv[i],"%i",varptr)!=1){
	fprintf(stderr, "Bad argument %s\n", i<argc ? argv[i] : argv[argc-1]);
	exit(-1);
    }
}

void WalksatAlgorithm::scanoneu(int argc, char *argv[], int i, unsigned int *varptr)
{
    if (i>=argc || sscanf(argv[i],"%u",varptr)!=1){
	fprintf(stderr, "Bad argument %s\n", i<argc ? argv[i] : argv[argc-1]);
	exit(-1);
    }
}

void WalksatAlgorithm::scanonell(int argc, char *argv[], int i, BIGINT *varptr)
{
    if (i>=argc || sscanf(argv[i],"%Li",varptr)!=1){
	fprintf(stderr, "Bad argument %s\n", i<argc ? argv[i] : argv[argc-1]);
	exit(-1);
    }
}


int WalksatAlgorithm::calc_hamming_dist(int atom[], int hamming_target[], int numatom)
{
    int i;
    int dist = 0;
    
    for (i=1; i<=numatom; i++){
	if (atom[i] != hamming_target[i]) dist++;
    }
    return dist;
}

void WalksatAlgorithm::open_hamming_data(char initfile[])
{
    if ((hamming_fp = fopen(initfile, "w")) == NULL){
	fprintf(stderr, "Cannot open %s for output\n", initfile);
	exit(1);
    }
}


void WalksatAlgorithm::read_hamming_file(char initfile[])
{
    int i;			/* loop counter */
    FILE * infile;
    int lit;    

    printf("loading hamming target file %s ...", initfile);

    if ((infile = fopen(initfile, "r")) == NULL){
	fprintf(stderr, "Cannot open %s\n", initfile);
	exit(1);
    }
    i=0;
    for(i = 1;i < numatom+1;i++)
      hamming_target[i] = 0;

    while (fscanf(infile, " %d", &lit)==1){
	if (ABS(lit)>numatom){
	    fprintf(stderr, "Bad hamming file %s\n", initfile);
	    exit(1);
	}
	if (lit>0) hamming_target[lit]=1;
    }
    printf("done\n");
}


void WalksatAlgorithm::init(char initfile[], int initoptions)
{
    int i;
    int j;
    int thetruelit=0;
#ifdef _WATCH
    int truelit1;
    int truelit2;
#endif
    FILE * infile;
    int lit;

    for(i = 0;i < numclause;i++)
      numtruelit[i] = 0;
    numfalse = 0;

    for(i = 1;i < numatom+1;i++)
      {
	  changed[i] = -BIG;
	  breakcount[i] = 0;
	  makecount[i] = 0;
      }

    if (initfile[0] && initoptions!=INIT_PARTIAL){
	for(i = 1;i < numatom+1;i++)
	  atom[i] = 0;
    }
    else {
	for(i = 1;i < numatom+1;i++)
	  atom[i] = random()%2;
    }

    if (initfile[0]){
	if ((infile = fopen(initfile, "r")) == NULL){
	    fprintf(stderr, "Cannot open %s\n", initfile);
	    exit(1);
	}
	i=0;
	while (fscanf(infile, " %d", &lit)==1){
	    i++;
	    if (ABS(lit)>numatom){
		fprintf(stderr, "Bad init file %s\n", initfile);
		exit(1);
	    }
	    if (lit<0) atom[-lit]=0;
	    else atom[lit]=1;
	}
	if (i==0){
	    fprintf(stderr, "Bad init file %s\n", initfile);
	    exit(1);
	}
	fclose(infile);
    }

    /* Initialize breakcount and makecount in the following: */
    for(i = 0;i < numclause;i++)
      {
#ifdef _WATCH
	truelit1 = 0;
	truelit2 = 0;
#endif
	  for(j = 0;j < size[i];j++)
	    {
		if((clause[i][j] > 0) == atom[ABS(clause[i][j])])
		  {
		      numtruelit[i]++;
		      thetruelit = clause[i][j];
#ifdef _WATCH
		      if (!truelit1)
			truelit1 = clause[i][j];
		      else if (truelit1 && !truelit2)
			truelit2 = clause[i][j];
#endif		      
		  }
	    }
	  if(numtruelit[i] == 0)
	    {
		wherefalse[i] = numfalse;
		falsified[numfalse] = i;
		numfalse++;
		for(j = 0;j < size[i];j++){
		  makecount[ABS(clause[i][j])]++;
		}
	    }
	  else if (numtruelit[i] == 1)
	    {
		breakcount[ABS(thetruelit)]++;
#ifdef _WATCH
		watch1[i] = ABS(thetruelit);
#endif
	    }
#ifdef _WATCH
	  else /*if (numtruelit[i] == 2)*/
	    {
	      watch1[i] = ABS(truelit1);
	      watch2[i] = ABS(truelit2);
	    }
#endif
      }

    if (hamming_flag){
	hamming_distance = calc_hamming_dist(atom, hamming_target, numatom);
	fprintf(hamming_fp, "0 %i\n", hamming_distance);
    }
}

void WalksatAlgorithm::print_false_clauses(long int lowbad)
{
    int i, j;
    int cl;

    printf("Unsatisfied clauses:\n");
    for (i=0; i<lowbad; i++){
	cl = lowfalse[i];
	for (j=0; j<size[cl]; j++){
	    printf("%d ", clause[cl][j]);
	}
	printf("0\n");
    }
    printf("End unsatisfied clauses\n");
}

void WalksatAlgorithm::save_false_clauses(long int lowbad)
{
    int i;

    for (i=0; i<lowbad; i++)
      lowfalse[i] = falsified[i];
}

void WalksatAlgorithm::initprob(void)
{
    int i;
    int j;
    int lastc;
    int nextc;
    int *storeptr=NULL;
    int freestore;
    int lit;

    while ((lastc = getchar()) == 'c')
      {
	  while ((nextc = getchar()) != EOF && nextc != '\n');
      }
    ungetc(lastc,stdin);
    if (scanf("p cnf %i %i",&numatom,&numclause) != 2)
      {
	  fprintf(stderr,"Bad input file\n");
	  exit(-1);
      }
    if(numatom > MAXATOM)
      {
	  fprintf(stderr,"ERROR - too many atoms\n");
	  exit(-1);
      }

#ifdef DYNAMIC
    clause = (int **) malloc(sizeof(int *)*(numclause+1));
    size = (int *) malloc(sizeof(int)*(numclause+1));
    falsified = (int *) malloc(sizeof(int)*(numclause+1));
    lowfalse = (int *) malloc(sizeof(int)*(numclause+1));
    wherefalse = (int *) malloc(sizeof(int)*(numclause+1));
    numtruelit = (int *) malloc(sizeof(int)*(numclause+1));
#else
    if(numclause > MAXCLAUSE)                     
      {                                      
	  fprintf(stderr,"ERROR - too many clauses\n"); 
	  exit(-1);                              
      }                                        
#endif
    freestore = 0;
    numliterals = 0;
    for(i = 0;i < 2*MAXATOM+1;i++)
      numoccurence[i] = 0;
    for(i = 0;i < numclause;i++)
      {
	  size[i] = -1;
	  if (freestore < MAXLENGTH)
	    {
		storeptr = (int *) malloc( sizeof(int) * STOREBLOCK );
		freestore = STOREBLOCK;
		//fprintf(stderr,"allocating memory...\n");
	    }
	  clause[i] = storeptr;
	  do
	    {
		size[i]++;
		if(size[i] > MAXLENGTH)
		  {
		      printf("ERROR - clause too long\n");
		      exit(-1);
		  }
		if (scanf("%i ",&lit) != 1)
		  {
		      fprintf(stderr, "Bad input file\n");
		      exit(-1);
		  }
		if(lit != 0)
		  {
		      *(storeptr++) = lit; /* clause[i][size[i]] = j; */
		      freestore--;
		      numliterals++;
		      numoccurence[lit+MAXATOM]++;
		  }
	    }
	  while(lit != 0);
      }
    if(size[0] == 0)
      {
	  fprintf(stderr,"ERROR - incorrect problem format or extraneous characters\n");
	  exit(-1);
      }

    for(i = 0;i < 2*MAXATOM+1;i++)
      {
	  if (freestore < numoccurence[i])
	    {
		storeptr = (int *) malloc( sizeof(int) * STOREBLOCK );
		freestore = STOREBLOCK;
		//fprintf(stderr,"allocating memory...\n");
	    }
	  occurence[i] = storeptr;
	  freestore -= numoccurence[i];
	  storeptr += numoccurence[i];
	  numoccurence[i] = 0;
      }

    for(i = 0;i < numclause;i++)
      {
	  for(j = 0;j < size[i];j++)
	    {
		occurence[clause[i][j]+MAXATOM]
		  [numoccurence[clause[i][j]+MAXATOM]] = i;
		numoccurence[clause[i][j]+MAXATOM]++;
	    }
      }
}

#ifdef _WATCH
/* new flipping function based on SAT2004 submission work */
void WalksatAlgorithm::flipatom(int toflip)
{
    int i, j;			
    int toenforce;		
    register int cli;
    register int lit;
    int numocc;
    register int sz;
    register int * litptr;
    int * occptr;
    register int v;
    /* printf("flipping %i\n", toflip); */

    if (toflip == NOVALUE){
	numnullflip++;
	return;
    }

    changed[toflip] = numflip;
    if(atom[toflip] > 0)
      toenforce = -toflip;
    else
      toenforce = toflip;
    atom[toflip] = 1-atom[toflip];

    if (hamming_flag){
	if (atom[toflip] == hamming_target[toflip])
	  hamming_distance--;
	else
	  hamming_distance++;
	if ((numflip % hamming_sample_freq) == 0)
	  fprintf(hamming_fp, "%i %i\n", numflip, hamming_distance);
    }
    
    numocc = numoccurence[MAXATOM-toenforce];
    occptr = occurence[MAXATOM-toenforce];
    for(i = 0; i < numocc ;i++)
      {
	  /* cli = occurence[MAXATOM-toenforce][i]; */
	  cli = *(occptr++);

	  if (--numtruelit[cli] == 0){
	      falsified[numfalse] = cli;
	      wherefalse[cli] = numfalse;
	      numfalse++;
	      /* Decrement toflip's breakcount */
	      breakcount[toflip]--;

	      if (makeflag){
		/* Increment the makecount of all vars in the clause */
		sz = size[cli];
		litptr = clause[cli];
		for (j=0; j<sz; j++){
		  /* lit = clause[cli][j]; */
		  lit = *(litptr++);
		  makecount[ABS(lit)]++;
		}
	      }
	  }
	  else if (numtruelit[cli] == 1){
	    if (watch1[cli] == toflip) {
	      assert(watch1[cli] != watch2[cli]);
	      watch1[cli] = watch2[cli];
	    }
	    breakcount[watch1[cli]]++;
	  }
	  else { /* numtruelit[cli] >= 2 */
	    if (watch1[cli] == toflip) {
	      /* find a true literal other than watch1[cli] and watch2[cli] */
	      sz = size[cli];
	      litptr = clause[cli];
	      for (j=0; j<sz; j++) {
		lit = *(litptr++);
		v = ABS(lit);
		if ((lit>0) == atom[v] && v != watch1[cli] && v != watch2[cli]) {
		    watch1[cli] = v;
		    break;
		}
	      }
	    }
	    else if (watch2[cli] == toflip) {
	      /* find a true literal other than watch1[cli] and watch2[cli] */
	      sz = size[cli];
	      litptr = clause[cli];
	      for (j=0; j<sz; j++) {
		lit = *(litptr++);
		v =ABS(lit);
		if ((lit>0) == atom[v] && v != watch1[cli] && v != watch2[cli]) {
		    watch2[cli] = v;
		    break;
		}
	      }
	    }
	  }
      }
    
    numocc = numoccurence[MAXATOM+toenforce];
    occptr = occurence[MAXATOM+toenforce];
    for(i = 0; i < numocc; i++)
      {
	  /* cli = occurence[MAXATOM+toenforce][i]; */
	  cli = *(occptr++);

	  if (++numtruelit[cli] == 1){
	      numfalse--;
	      falsified[wherefalse[cli]] =
		falsified[numfalse];
	      wherefalse[falsified[numfalse]] =
		wherefalse[cli];
	      /* Increment toflip's breakcount */
	      breakcount[toflip]++;

	      if (makeflag){
		/* Decrement the makecount of all vars in the clause */
		sz = size[cli];
		litptr = clause[cli];
		for (j=0; j<sz; j++){
		  /* lit = clause[cli][j]; */
		  lit = *(litptr++);
		  makecount[ABS(lit)]--;
		}
	      }
	      watch1[cli] = toflip;
	  }
	  else if (numtruelit[cli] == 2){
	    watch2[cli] = toflip;
	    breakcount[watch1[cli]]--;
	  }
      }
}

#else
/*
  the "standard" flipatom implementation
*/
void WalksatAlgorithm::flipatom(int toflip)
{
    int i, j;			
    int toenforce;		
    register int cli;
    register int lit;
    int numocc;
    register int sz;
    register int * litptr;
    int * occptr;

    /* printf("flipping %i\n", toflip); */

    if (toflip == NOVALUE){
	numnullflip++;
	return;
    }

    changed[toflip] = numflip;
    if(atom[toflip] > 0)
      toenforce = -toflip;
    else
      toenforce = toflip;
    atom[toflip] = 1-atom[toflip];

    if (hamming_flag){
	if (atom[toflip] == hamming_target[toflip])
	  hamming_distance--;
	else
	  hamming_distance++;
	if ((numflip % hamming_sample_freq) == 0)
	  fprintf(hamming_fp, "%" BIGINTSTR " %i\n", numflip, hamming_distance);
    }
    
    numocc = numoccurence[MAXATOM-toenforce];
    occptr = occurence[MAXATOM-toenforce];
    for(i = 0; i < numocc ;i++)
      {
	  /* cli = occurence[MAXATOM-toenforce][i]; */
	  cli = *(occptr++);

	  if (--numtruelit[cli] == 0){
	      falsified[numfalse] = cli;
	      wherefalse[cli] = numfalse;
	      numfalse++;
	      /* Decrement toflip's breakcount */
	      breakcount[toflip]--;

	      if (makeflag){
		/* Increment the makecount of all vars in the clause */
		sz = size[cli];
		litptr = clause[cli];
		for (j=0; j<sz; j++){
		  /* lit = clause[cli][j]; */
		  lit = *(litptr++);
		  makecount[ABS(lit)]++;
		}
	      }
	  }
	  else if (numtruelit[cli] == 1){
	      /* Find the lit in this clause that makes it true, and inc its breakcount */
	      sz = size[cli];
	      litptr = clause[cli];
	      for (j=0; j<sz; j++){
		  /* lit = clause[cli][j]; */
		  lit = *(litptr++);
		  if((lit > 0) == atom[ABS(lit)]){
		      breakcount[ABS(lit)]++;
		      break;
		  }
	      }
	  }
      }
    
    numocc = numoccurence[MAXATOM+toenforce];
    occptr = occurence[MAXATOM+toenforce];
    for(i = 0; i < numocc; i++)
      {
	  /* cli = occurence[MAXATOM+toenforce][i]; */
	  cli = *(occptr++);

	  if (++numtruelit[cli] == 1){
	      numfalse--;
	      falsified[wherefalse[cli]] =
		falsified[numfalse];
	      wherefalse[falsified[numfalse]] =
		wherefalse[cli];
	      /* Increment toflip's breakcount */
	      breakcount[toflip]++;

	      if (makeflag){
		/* Decrement the makecount of all vars in the clause */
		sz = size[cli];
		litptr = clause[cli];
		for (j=0; j<sz; j++){
		  /* lit = clause[cli][j]; */
		  lit = *(litptr++);
		  makecount[ABS(lit)]--;
		}
	      }
	  }
	  else if (numtruelit[cli] == 2){
	      /* Find the lit in this clause other than toflip that makes it true,
		 and decrement its breakcount */
	      sz = size[cli];
	      litptr = clause[cli];
	      for (j=0; j<sz; j++){
		  /* lit = clause[cli][j]; */
		  lit = *(litptr++);
		  if( ((lit > 0) == atom[ABS(lit)]) &&
		     (toflip != ABS(lit)) ){
		      breakcount[ABS(lit)]--;
		      break;
		  }
	      }
	  }
      }
}
#endif

#if 0  
/* The following is the Watch-1 flip implementation, which watches
   1 true literal; it is slower than the Watch-2 implementation,
   but included for reference*/
void WalksatAlgorithm::flipatom_watch1(int toflip)
{
    int i, j;			
    int toenforce;		
    register int cli;
    register int lit;
    int numocc;
    register int sz;
    register int * litptr;
    int * occptr;

    /* printf("flipping %i\n", toflip); */

    if (toflip == NOVALUE){
	numnullflip++;
	return;
    }

    changed[toflip] = numflip;
    if(atom[toflip] > 0)
      toenforce = -toflip;
    else
      toenforce = toflip;
    atom[toflip] = 1-atom[toflip];

    if (hamming_flag){
	if (atom[toflip] == hamming_target[toflip])
	  hamming_distance--;
	else
	  hamming_distance++;
	if ((numflip % hamming_sample_freq) == 0)
	  fprintf(hamming_fp, "%i %i\n", numflip, hamming_distance);
    }
    
    numocc = numoccurence[MAXATOM-toenforce];
    occptr = occurence[MAXATOM-toenforce];
    for(i = 0; i < numocc ;i++)
      {
	  /* cli = occurence[MAXATOM-toenforce][i]; */
	  cli = *(occptr++);

	  if (--numtruelit[cli] == 0){
	      falsified[numfalse] = cli;
	      wherefalse[cli] = numfalse;
	      numfalse++;
	      /* Decrement toflip's breakcount */
	      breakcount[toflip]--;

	      if (makeflag){
		/* Increment the makecount of all vars in the clause */
		sz = size[cli];
		litptr = clause[cli];
		for (j=0; j<sz; j++){
		  /* lit = clause[cli][j]; */
		  lit = *(litptr++);
		  makecount[ABS(lit)]++;
		}
	      }
	  }
	  else if (numtruelit[cli] == 1){
	      /* Find the lit in this clause that makes it true, and inc its breakcount */
	      sz = size[cli];
	      litptr = clause[cli];
	      for (j=0; j<sz; j++){
		  /* lit = clause[cli][j]; */
		  lit = *(litptr++);
		  if((lit > 0) == atom[ABS(lit)]){
		    watch1[cli] = ABS(lit);
		      breakcount[ABS(lit)]++;
		      break;
		  }
	      }
	  }
      }
    
    numocc = numoccurence[MAXATOM+toenforce];
    occptr = occurence[MAXATOM+toenforce];
    for(i = 0; i < numocc; i++)
      {
	  /* cli = occurence[MAXATOM+toenforce][i]; */
	  cli = *(occptr++);

	  if (++numtruelit[cli] == 1){
	      numfalse--;
	      falsified[wherefalse[cli]] =
		falsified[numfalse];
	      wherefalse[falsified[numfalse]] =
		wherefalse[cli];
	      /* Increment toflip's breakcount */
	      breakcount[toflip]++;
	      watch1[cli] = toflip;

	      if (makeflag){
		/* Decrement the makecount of all vars in the clause */
		sz = size[cli];
		litptr = clause[cli];
		for (j=0; j<sz; j++){
		  /* lit = clause[cli][j]; */
		  lit = *(litptr++);
		  makecount[ABS(lit)]--;
		}
	      }
	  }
	  else if (numtruelit[cli] == 2){
	    breakcount[watch1[cli]]--;
	  }
      }
}
#endif

int WalksatAlgorithm::pickrandom(void)
{
    int tofix;

    tofix = falsified[random()%numfalse];
    return Var(tofix, random()%size[tofix]);
}

int WalksatAlgorithm::pickbest(void)
{
  int numbreak;
    int tofix;
    int clausesize;
    int i;		
    int best[MAXLENGTH];
    register int numbest;
    register int bestvalue;
    register int var;

    tofix = falsified[random()%numfalse];
    clausesize = size[tofix];
    numbest = 0;
    bestvalue = BIG;

    for (i=0; i< clausesize; i++){
      var = ABS(clause[tofix][i]);
      numbreak = breakcount[var];
      if (numbreak<=bestvalue){
	if (numbreak<bestvalue) numbest=0;
	bestvalue = numbreak;
	best[numbest++] = var;
      }
    }

    if (bestvalue>0 && (random()%denominator < numerator))
      return ABS(clause[tofix][random()%clausesize]);

    if (numbest == 1) return best[0];
    return best[random()%numbest];
}

int WalksatAlgorithm::picknovelty(void)
{
  int var, diff, birthdate;
  int youngest=0, youngest_birthdate=0, best=0, second_best=0, best_diff=0, second_best_diff=0;
  int tofix, clausesize, i;

  tofix = falsified[random()%numfalse];
  clausesize = size[tofix];  

  if (clausesize == 1) return ABS(clause[tofix][0]);

  youngest_birthdate = -1;
  best_diff = -BIG;
  second_best_diff = -BIG;

  for(i = 0; i < clausesize; i++){
    var = ABS(clause[tofix][i]);
    diff = makecount[var] - breakcount[var];
    birthdate = changed[var];
    if (birthdate > youngest_birthdate){
      youngest_birthdate = birthdate;
      youngest = var;
    }
    if (diff > best_diff || (diff == best_diff && changed[var] < changed[best])) {
      /* found new best, demote best to 2nd best */
      second_best = best;
      second_best_diff = best_diff;
      best = var;
      best_diff = diff;
    }
    else if (diff > second_best_diff || (diff == second_best_diff && changed[var] < changed[second_best])){
      /* found new second best */
      second_best = var;
      second_best_diff = diff;
    }
  }
  if (best != youngest) return best;
  if ((random()%denominator < numerator)) return second_best;
  return best;
}

int WalksatAlgorithm::pickrnovelty(void)
{
  int var, diff, birthdate;
  int diffdiff;
  int youngest=0, youngest_birthdate=0, best=0, second_best=0, best_diff=0, second_best_diff=0;
  int tofix, clausesize, i;

  tofix = falsified[random()%numfalse];
  clausesize = size[tofix];  

  if (clausesize == 1) return ABS(clause[tofix][0]);
  if ((numflip % 100) == 0) return ABS(clause[tofix][random()%clausesize]);

  youngest_birthdate = -1;
  best_diff = -BIG;
  second_best_diff = -BIG;

  for(i = 0; i < clausesize; i++){
    var = ABS(clause[tofix][i]);
    diff = makecount[var] - breakcount[var];
    birthdate = changed[var];
    if (birthdate > youngest_birthdate){
      youngest_birthdate = birthdate;
      youngest = var;
    }
    if (diff > best_diff || (diff == best_diff && changed[var] < changed[best])) {
      /* found new best, demote best to 2nd best */
      second_best = best;
      second_best_diff = best_diff;
      best = var;
      best_diff = diff;
    }
    else if (diff > second_best_diff || (diff == second_best_diff && changed[var] < changed[second_best])){
      /* found new second best */
      second_best = var;
      second_best_diff = diff;
    }
  }
  if (best != youngest) return best;

  diffdiff = best_diff - second_best_diff;

/*
  if (numerator < 50 && diffdiff > 1) return best;
  if (numerator < 50 && diffdiff == 1){
 */

  if (numerator*2 < denominator && diffdiff > 1) return best;
  if (numerator*2 < denominator && diffdiff == 1){
    if ((random()%denominator) < 2*numerator) return second_best;
    return best;
  }
  if (diffdiff == 1) return second_best;

/*
  if ((random()%denominator) < 2*(numerator-50)) return second_best;
*/

  if ((random()%denominator) < 2*(numerator-(denominator/2))) return second_best;

  return best;
}


/* hh: the following are variants of novelty and rnovelty that are typically
  more robust than the original heuristics.

  For details, see: 

  Holger H. Hoos and Thomas Stuetzle:
  Local Search Algorithms for SAT: An Empirical Evaluation. 
  Journal of Automated Reasoning, Vol.24, pp.421-481, 2000
  A modified version appeared as a book chapter in: 
  I.P.Gent, H.v.Maaren, T.Walsh, editors, SAT 2000, pp.43-86, IOS Press, 2000. 

  Holger H. Hoos:
  On the Run-time Behaviour of Stochastic Local Search Algorithms for SAT. 
  Proc. of AAAI-99, pp.661-666, MIT Press, 1999

  Holger H. Hoos:
  Stochastic Local Search - Methods, Models, Applications. 
  infix-Verlag, Sankt Augustin, Germany, 1999; ISBN 3-89601-215-0.

  Kautz modified slightly, to fix bug when denominator is not 100

*/

int WalksatAlgorithm::picknoveltyplus(void)
{
  int var, diff, birthdate;
  int youngest=0, youngest_birthdate=0, best=0, second_best=0, best_diff=0, second_best_diff=0;
  int tofix, clausesize, i;

  tofix = falsified[random()%numfalse];
  clausesize = size[tofix];  

  if (clausesize == 1) return ABS(clause[tofix][0]);

  /* hh: inserted modified loop breaker: */
  if ((random()%wp_denominator < wp_numerator)) return ABS(clause[tofix][random()%clausesize]);

  youngest_birthdate = -1;
  best_diff = -BIG;
  second_best_diff = -BIG;

  for(i = 0; i < clausesize; i++){
    var = ABS(clause[tofix][i]);
    diff = makecount[var] - breakcount[var];
    birthdate = changed[var];
    if (birthdate > youngest_birthdate){
      youngest_birthdate = birthdate;
      youngest = var;
    }
    if (diff > best_diff || (diff == best_diff && changed[var] < changed[best])) {
      /* found new best, demote best to 2nd best */
      second_best = best;
      second_best_diff = best_diff;
      best = var;
      best_diff = diff;
    }
    else if (diff > second_best_diff || (diff == second_best_diff && changed[var] < changed[second_best])){
      /* found new second best */
      second_best = var;
      second_best_diff = diff;
    }
  }
  if (best != youngest) return best;
  if ((random()%denominator < numerator)) return second_best;
  return best;
}

int WalksatAlgorithm::pickrnoveltyplus(void)
{
  int var, diff, birthdate;
  int diffdiff;
  int youngest=0, youngest_birthdate=0, best=0, second_best=0, best_diff=0, second_best_diff=0;
  int tofix, clausesize, i;

  tofix = falsified[random()%numfalse];
  clausesize = size[tofix];  

  if (clausesize == 1) return ABS(clause[tofix][0]);

/* hh: modified loop breaker: */
  if ((random()%wp_denominator < wp_numerator)) return ABS(clause[tofix][random()%clausesize]);

  /*
  if ((numflip % 100) == 0) return ABS(clause[tofix][random()%clausesize]);
  */

  youngest_birthdate = -1;
  best_diff = -BIG;
  second_best_diff = -BIG;

  for(i = 0; i < clausesize; i++){
    var = ABS(clause[tofix][i]);
    diff = makecount[var] - breakcount[var];
    birthdate = changed[var];
    if (birthdate > youngest_birthdate){
      youngest_birthdate = birthdate;
      youngest = var;
    }
    if (diff > best_diff || (diff == best_diff && changed[var] < changed[best])) {
      /* found new best, demote best to 2nd best */
      second_best = best;
      second_best_diff = best_diff;
      best = var;
      best_diff = diff;
    }
    else if (diff > second_best_diff || (diff == second_best_diff && changed[var] < changed[second_best])){
      /* found new second best */
      second_best = var;
      second_best_diff = diff;
    }
  }
  if (best != youngest) return best;

  diffdiff = best_diff - second_best_diff;

/*
  if (numerator < 50 && diffdiff > 1) return best;
  if (numerator < 50 && diffdiff == 1){
 */

  if (numerator*2 < denominator && diffdiff > 1) return best;
  if (numerator*2 < denominator && diffdiff == 1){
    if ((random()%denominator) < 2*numerator) return second_best;
    return best;
  }
  if (diffdiff == 1) return second_best;

/*
  if ((random()%denominator) < 2*(numerator-50)) return second_best;
*/

  if ((random()%denominator) < 2*(numerator-(denominator/2))) return second_best;

  return best;

}



int WalksatAlgorithm::picktabu(void)
{
    int numbreak[MAXLENGTH];
    int tofix;
    int clausesize;
    int i;			/* a loop counter */
    int best[MAXLENGTH];	/* best possibility so far */
    int numbest;		/* how many are tied for best */
    int bestvalue;		/* best value so far */
    int noisypick;

    tofix = falsified[random()%numfalse];
    clausesize = size[tofix];
    for(i = 0;i < clausesize;i++)
	numbreak[i] = breakcount[ABS(clause[tofix][i])];

    numbest = 0;
    bestvalue = BIG;

    noisypick = (numerator > 0 && random()%denominator < numerator); 
    for (i=0; i < clausesize; i++) {
	if (numbreak[i] == 0) {
	    if (bestvalue > 0) {
		bestvalue = 0;
		numbest = 0;
	    }
	    best[numbest++] = i;
	}
	else if (tabu_length < numflip - changed[ABS(clause[tofix][i])]) {
	    if (noisypick && bestvalue > 0) { 
		best[numbest++] = i; 
	    }
	    else {
		if (numbreak[i] < bestvalue) {
		    bestvalue = numbreak[i];
		    numbest = 0;
		}
		if (numbreak[i] == bestvalue) {
		    best[numbest++] = i; 
		}
	    }
	}
    }
    if (numbest == 0) return NOVALUE;
    if (numbest == 1) return Var(tofix, best[0]);
    return (Var(tofix, best[random()%numbest]));
}

int WalksatAlgorithm::countunsat(void)
	{
	int i, j, unsat, bad, lit, sign;

	unsat = 0;
	for (i=0;i < numclause;i++)
		{
		bad = TRUE;
		for (j=0; j < size[i]; j++)
			{
			lit = clause[i][j];
			sign = lit > 0 ? 1 : 0;
			if ( atom[ABS(lit)] == sign )
				{
				bad = FALSE;
				break;
				}
			}
		if (bad)
			unsat++;
		}
	return unsat;
	}


double WalksatAlgorithm::elapsed_seconds(void) 
{ 
	double answer;

#if ANSI
	static long prev_time = 0;
	time( &long_time );
	/* Note:  time(&t) returns t in seconds, so do not need to /CLK_TCK */
	answer = long_time - prev_time;
	prev_time = long_time;
#endif
#if NT
	static DWORD prev_time = 0;
	win_time = timeGetTime();
	/* Note:  return value of timeGetTime() is ms, so divide by 1000*/
	answer = (double)(win_time - prev_time) / (double)1000; 
	prev_time = win_time;
#endif
#if UNIX || OSX
    static struct tms prog_tms;
    static long prev_times = 0;
    (void) times(&prog_tms);
    answer = ((double)(((long)prog_tms.tms_utime)-prev_times))/((double) ticks_per_second);
    prev_times = (long) prog_tms.tms_utime;
#endif
	return answer; 
}




void WalksatAlgorithm::print_sol_cnf(void)
{
    int i;
    for(i = 1;i < numatom+1;i++)
	printf("v %i\n", solution[i] == 1 ? i : -i);
}


void WalksatAlgorithm::print_sol_file(char * filename)
{
    FILE * fp;
    int i;

    if ((fp = fopen(filename, "w"))==NULL){
	fprintf(stderr, "Cannot open output file\n");
	exit(-1);
    }
    for(i = 1;i < numatom+1;i++){
	fprintf(fp, " %i", solution[i] == 1 ? i : -i);
	if (i % 10 == 0) fprintf(fp, "\n");
    }
    if ((i-1) % 10 != 0) fprintf(fp, "\n");
    fclose(fp);
}


void WalksatAlgorithm::print_low_assign(long int lowbad)
{
    int i;

    printf("Begin assign with lowest # bad = %ld\n", lowbad);
    for (i=1; i<=numatom; i++){
	printf(" %d", lowatom[i]==0 ? -i : i);
	if (i % 10 == 0) printf("\n");
    }
    if ((i-1) % 10 != 0) printf("\n");
    printf("End assign\n");
}

void WalksatAlgorithm::print_current_assign(void)
{
    int i;

    printf("Begin assign at flip = %" BIGINTSTR "\n", numflip);
    for (i=1; i<=numatom; i++){
	printf(" %d", atom[i]==0 ? -i : i);
	if (i % 10 == 0) printf("\n");
    }
    if ((i-1) % 10 != 0) printf("\n");
    printf("End assign\n");
}

void WalksatAlgorithm::save_low_assign(void)
{
    int i;

    for (i=1; i<=numatom; i++)
      lowatom[i] = atom[i];
}

void WalksatAlgorithm::save_solution(void)
{
    int i;

    for (i=1; i<=numatom; i++)
      solution[i] = atom[i];
}


