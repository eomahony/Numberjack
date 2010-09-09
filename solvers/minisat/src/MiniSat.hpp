
#ifndef MINISAT_H
#define MINISAT_H


#include <vector>

#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>

#include <cstdio>

#include <Queue.h>
#include <Solver.h>

#include <ctime>
#include <cstring>
#include <stdint.h>
#include <errno.h>

#include <signal.h>
#include <zlib.h>

#include "SimpSolver.hpp"
#include "SatWrapper.hpp"




/**
   The solver itself
*/
class MiniSatSolver : public SatWrapperSolver, public SimpSolver
{

private:
  ////////////// MiniSat Specific ////////////////
  double STARTTIME;
  lbool result;
  int nbSolutions;

  // search stuff
  Clause*     conflict_clause;
  int         backtrack_level;
  int         conflictC;
  //vec<Lit>    forced_decisions;
  std::vector<Lit>    forced_decisions;
  vec<Lit>    learnt_clause;

  int first_decision_level;
  Lit last_decision;
  int saved_level;
  ////////////// MiniSat Specific ////////////////

public:

  MiniSatSolver();
  virtual ~MiniSatSolver();
  //int nbClauses() {return nClauses();}

  // add an expression, in the case of a tree of expressions,
  // each node of the tree is added separately, depth first.
  int create_atom(DomainEncoding* dom, const int type);
  void validate();
  lbool truth_value(Lit x);

  // solving methods
  int solve();
  int solveAndRestart(const int policy = GEOMETRIC, 
		      const unsigned int base = 32, 
		      const double factor = 1.3333333,
		      const double decay = 0.0);

  void reset(bool full);
  bool propagate();
  void save();
  void post(const char* op, SatWrapper_Expression* x, int v);
  bool undo(const int nlevel);
  void deduce();
  bool branch_right();

  void store_solution();
  
  // parameter tuning methods
  void setFailureLimit(const int cutoff);  
  void setNodeLimit(const int cutoff);  
  void setTimeLimit(const double cutoff);
  void setVerbosity(const int degree);
  void setRandomized(const int degree);
  void setRandomSeed(const int seed);

  // statistics methods
  bool is_sat();
  bool is_opt();
  bool is_unsat();
  void printStatistics();
  int getBacktracks();
  int getNodes();
  int getFailures();
  int getPropags();
  double getTime();
};


#endif

