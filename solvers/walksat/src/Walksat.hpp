
#ifndef WALKSAT_H
#define WALKSAT_H


#include <vector>
#include <iostream>

#include "cpp_walksat.hpp"
#include "SatWrapper.hpp"


/**
   The solver itself
*/
class WalksatSolver : public SatWrapperSolver 
{

private:
  ////////////// Walksat Specific ////////////////
  WalksatAlgorithm wsat;

  double STARTTIME;
  int nbSolutions;
  ////////////// Walksat Specific ////////////////

public:

  WalksatSolver();
  virtual ~WalksatSolver();

  lbool truth_value(Lit x);

  // used to initialise search on a given subset of variables
  void initialise(SatWrapperExpArray& arg);
  // initialise the solver before solving (no more calls to add after this)
  void initialise();

  // solving methods
  int solve();
  int solveAndRestart(const int policy = GEOMETRIC, 
		      const unsigned int base = 32, 
		      const double factor = 1.3333333,
		      const double decay = 0.0);
  void store_solution();
  

  void setFailureLimit(const int cutoff);  
  void setNodeLimit(const int cutoff);  
  void setTimeLimit(const double cutoff);
  void setRestartLimit(const int cutoff);
  void setVerbosity(const int degree);
  void setRandomSeed(const int seed);

  // statistics methods
  bool is_sat();
  void printStatistics();
  double getTime();
};


#endif

