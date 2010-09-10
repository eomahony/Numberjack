
#include "Walksat.hpp"

#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>



static inline double cpuTime(void) {
    struct rusage ru;
    getrusage(RUSAGE_SELF, &ru);
    return (double)ru.ru_utime.tv_sec + (double)ru.ru_utime.tv_usec / 1000000; }



/**************************************************************
 ********************     Solver        ***********************
 **************************************************************/

WalksatSolver::WalksatSolver() : SatWrapperSolver()
{

#ifdef _DEBUGWRAP
  std::cout << "create a minisat solver" << std::endl;
#endif

  ////////////// Walksat Specific ////////////////
  STARTTIME = cpuTime();
  nbSolutions = 0;
  ////////////// Walksat Specific ////////////////

}

WalksatSolver::~WalksatSolver()
{

#ifdef _DEBUGWRAP
  std::cout << "delete wrapped solver" << std::endl;
#endif

}

void WalksatSolver::store_solution() {
#ifdef _DEBUGWRAP
  std::cout << "store a new solution" << std::endl;
#endif
  
  ++nbSolutions;
  if(!cp_model) cp_model = new int[_expressions.size()];
  for(unsigned int i=0; i<_variables.size(); ++i) {
    if(_variables[i]) {
      cp_model[_variables[i]->_ident] = _variables[i]->get_min();
    } else
      cp_model[_variables[i]->_ident] = 0;
  }
}

void WalksatSolver::initialise(SatWrapperExpArray& arg)
{

  initialise();

}

void WalksatSolver::initialise()
{
  int i, j;
  
  int *storeptr = NULL;
  int freestore;
  int lit;
  //int simplified=0;

  std::vector<Lit> clause;

#ifdef _DEBUGWRAP
  std::cout << "initialise the solver" << std::endl;
#endif

//   for(i=0; (unsigned int)i<clause_base.size(); ++i)
//     {
//       displayClause(clause_base[i]);
//     }
  
  wsat.numatom = _atom_to_domain.size()-1;
  wsat.numclause = clause_base.size();
  
#ifdef DYNAMIC
  wsat.clause = (int **) malloc(sizeof(int *)*(wsat.numclause+1));
  wsat.size = (int *) malloc(sizeof(int)*(wsat.numclause+1));
  wsat.falsified = (int *) malloc(sizeof(int)*(wsat.numclause+1));
  wsat.lowfalse = (int *) malloc(sizeof(int)*(wsat.numclause+1));
  wsat.wherefalse = (int *) malloc(sizeof(int)*(wsat.numclause+1));
  wsat.numtruelit = (int *) malloc(sizeof(int)*(wsat.numclause+1));
#else
  if(wsat.numclause > MAXCLAUSE)                     
    {                                      
      fprintf(stderr,"ERROR - too many clauses\n"); 
      exit(-1);                              
    }                                        
#endif

  freestore = 0;
  wsat.numliterals = 0;
  for(i = 0;i < 2*MAXATOM+1;i++)
    wsat.numoccurence[i] = 0;


  for(i = 0;i < wsat.numclause;i++)
    {
//       //if(!processClause(clause_base[i],clause)) {

//       std::cout << "clause " << (i+1) << " out of " << wsat.numclause << std::endl;
//       displayClause(clause_base[i]);

      wsat.size[i] = -1;
      if (freestore < MAXLENGTH)
	{
	  storeptr = (int *) malloc( sizeof(int) * STOREBLOCK );
	  freestore = STOREBLOCK;
	  //fprintf(stderr,"allocating memory...\n");
	}
      wsat.clause[i] = storeptr;
      if(clause_base[i].size() > MAXLENGTH)
	{
	  printf("ERROR - clause too long\n");
	  exit(-1);
	}
      else 
	{
	  wsat.size[i] = clause_base[i].size();
	  for(j = 0;j<wsat.size[i];j++)
	    {
	      lit = (sign(clause_base[i][j]) ? -1 : 1)*var(clause_base[i][j]);
	      *(storeptr++) = lit;
	      freestore--;
	      wsat.numliterals++;
	      wsat.numoccurence[lit+MAXATOM]++;
	    }
	}
      //} else ++simplified;
    }
  //wsat.numclause -= simplified;

  for(i = 0;i < 2*MAXATOM+1;i++)
    {
      if (freestore < wsat.numoccurence[i])
	{
	  storeptr = (int *) malloc( sizeof(int) * STOREBLOCK );
	  freestore = STOREBLOCK;
	  fprintf(stderr,"allocating memory...\n");
	}
      wsat.occurence[i] = storeptr;
      freestore -= wsat.numoccurence[i];
      storeptr += wsat.numoccurence[i];
      wsat.numoccurence[i] = 0;
    }

  for(i = 0;i < wsat.numclause;i++)
    {
      for(j = 0;j < wsat.size[i];j++)
	{
	  wsat.occurence[wsat.clause[i][j]+MAXATOM]
	    [wsat.numoccurence[wsat.clause[i][j]+MAXATOM]] = i;
	  wsat.numoccurence[wsat.clause[i][j]+MAXATOM]++;
	}
    } 
}

lbool WalksatSolver::truth_value(Lit x)
{
  return (wsat.solution[var(x)] == 1 ? l_True : l_False);
}

int WalksatSolver::solveAndRestart(const int policy, 
				   const unsigned int base, 
				   const double factor,
				   const double decay)
{
  return solve();
}

int WalksatSolver::solve()
{

#ifdef _DEBUGWRAP
  std::cout << "call solve" << std::endl;  
#endif 

#ifdef _DEBUGWRAP
  std::cout << "solve!" << std::endl;  
#endif 

  STARTTIME = cpuTime();

  wsat.initialize_statistics();
  if(wsat.verbosity) wsat.print_statistics_header();
  wsat.walk_solve();

  return is_sat();
}

void WalksatSolver::setFailureLimit(const int cutoff)
{
  wsat.cutoff = cutoff;
}

void WalksatSolver::setNodeLimit(const int cutoff)
{
  wsat.cutoff = cutoff;
}

void WalksatSolver::setTimeLimit(const double cutoff)
{
  wsat.time_cutoff = cutoff;
}

void WalksatSolver::setRestartLimit(const int cutoff)
{
  wsat.numrun = cutoff;
}

void WalksatSolver::setVerbosity(const int degree)
{
  wsat.verbosity = degree;
}

void WalksatSolver::setRandomSeed(const int seed)
{
  srandom(seed);
}

bool WalksatSolver::is_sat()
{
  return wsat.numsuccesstry > 0;
}

void WalksatSolver::printStatistics()
{
  wsat.print_statistics_final();
}

double WalksatSolver::getTime()
{
  return cpuTime() - STARTTIME;
}


