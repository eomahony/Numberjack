%module(package="Numberjack.solvers") SCIP 
%import(module="Numberjack.solvers.MipWrapper") "MipWrapper.hpp"

%{
#include <vector>
#include <scip/scip.h>
#include <scip/scipdefplugins.h>
#include <scip/type_stat.h>
#include <scip/struct_stat.h>
#include <scip/clock.h>
#include "MipWrapper.hpp"
%}

%inline %{


/**
   The solver itself
*/
class SCIPSolver : public MipWrapperSolver{
private:
  
  SCIP * _scip;
  int _verbosity;
  SCIP* get_scip();
  void add_in_constraint(LinearConstraint *con, double coef=0);
  void output_orig_problem(const char *filename);
  bool has_been_added;

public:
  int var_counter;

  SCIPSolver();
  virtual ~SCIPSolver();
  
  // initialise the solver before solving (no more calls to add after this)
  void initialise(MipWrapperExpArray& arg);
  void initialise();

  // solving methods
  int solve();
  
  // parameter tuning methods
  void setTimeLimit(const int cutoff);
  void setOptimalityGap(const double gap);
  void setVerbosity(const int degree);

  // statistics methods
  int getNodes();
  bool is_sat();
  bool is_unsat();
  bool is_opt();
  void output_lp(const char *filename);
  void output_mps(const char *filename);
  void printStatistics();
  double getTime();
  double getOptimalityGap();
  
  // Value stuff
  virtual double get_value(void *ptr);

};


%}


%pythoncode %{
import Numberjack

class Solver(Numberjack.NBJ_STD_Solver):
    def __init__(self, model=None, X=None, FD=False, clause_limit=-1, encoding=None):
        Numberjack.NBJ_STD_Solver.__init__(self, "SCIP", "MipWrapper", model, X, FD, clause_limit, encoding)
%}
