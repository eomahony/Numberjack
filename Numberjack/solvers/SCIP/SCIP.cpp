
#include <iostream>

#include "SCIP.hpp"
#include "scip_exception.hpp"
#include <math.h>


/**************************************************************
 ********************     Solver        ***********************
 **************************************************************/

SCIPSolver::SCIPSolver() /*: scipvariables(), scipcons() */{

  DBG("create a scip solver\n%s", "");

  var_counter = 0;
  _verbosity = 0;
  has_been_added = false;

  // Load up SCIP
  SCIP_CALL_EXC( SCIPcreate(& _scip) );
  // load default plugins linke separators, heuristics, etc.
  SCIP_CALL_EXC( SCIPincludeDefaultPlugins(_scip) );
  // create an empty problem
  SCIP_CALL_EXC( SCIPcreateProb(_scip, "Numberjack Model",
				NULL, NULL, NULL, NULL, NULL, NULL, NULL) );
  // set the objective sense to maximize, default is minimize
  SCIP_CALL_EXC( SCIPsetObjsense(_scip, SCIP_OBJSENSE_MAXIMIZE) );
}

SCIPSolver::~SCIPSolver(){
  DBG("delete wrapped solver\n%s", "");

  // Could not figure out the correct way to release SCIP's memory, kept causing seg faults
  //
  // std::cout << "deleting scip varaibles " << scipvariables.size() << std::endl;
  // for(unsigned int i=0; i<scipvariables.size(); ++i){
  //   // std::cout << "releasing " << (&scipvariables->at(i)) << std::endl;
  //   std::cout << "releasing " << (&scipvariables[i]) << std::endl;
  //   // SCIP_CALL_EXC( SCIPreleaseVar(_scip, &(scipvariables->at(i))) );
  //   SCIP_CALL_EXC( SCIPreleaseVar(_scip, &scipvariables[i]) );
  // }
  // std::cout << "done." << std::endl;
  // scipvariables.clear();
  // delete scipvariables;

  //  for (std::vector<SCIP_CONS *>::size_type i = 0; i < scipcons->size(); ++i)
  //     SCIP_CALL_EXC( SCIPreleaseCons(_scip, &(scipcons->at(i))) );
  // scipcons->clear();
  // delete scipcons;

  // SCIP_CALL_EXC( SCIPfree( & _scip) );
}

SCIP* SCIPSolver::get_scip() {return _scip;}

void SCIPSolver::initialise(MipWrapperExpArray& arg){
  DBG("Initialise the Solver with search vars %s\n", "");

  initialise();
}

void SCIPSolver::initialise(){
  DBG("initialise the solver%s\n", "");
  MipWrapperSolver::initialise();
  if(_obj != NULL) add_in_constraint(_obj, _obj_coef);
  for(unsigned int i = 0; i < _constraints.size(); ++i)
    add_in_constraint(_constraints[i]);
  has_been_added = true;
}

void SCIPSolver::add_in_constraint(LinearConstraint *con, double coef){
  DBG("Creating a SCIP representation of a constriant%s\n", "");
  
  double *weights = new double[con->_coefficients.size()];
  SCIP_VAR** vars = new SCIP_VAR*[con->_variables.size()];
  
  for(unsigned int i = 0; i < con->_variables.size(); ++i){
    
    DBG("\tAdding variable to SCIP\n%s", "");
    
    if(con->_variables[i]->_var == NULL){
    
      SCIP_VAR* var_ptr;
      SCIP_VARTYPE type;
    
      if(con->_variables[i]->_continuous) type = SCIP_VARTYPE_CONTINUOUS;
      else type = SCIP_VARTYPE_INTEGER;
    
      SCIP_CALL_EXC(SCIPcreateVar(_scip, &var_ptr,
				  "SCIP_Var",
				  con->_variables[i]->_lower, // LB
				  con->_variables[i]->_upper, // UB
				  coef * con->_coefficients[i], // ective
				  type,
				  TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL_EXC( SCIPaddVar(_scip, var_ptr) );
      // scipvariables.push_back(var_ptr);
    
      con->_variables[i]->_var = (void*) var_ptr;
      vars[i] = var_ptr;
      weights[i] = con->_coefficients[i];
    } else {
      vars[i] = (SCIP_VAR*) con->_variables[i]->_var;
      weights[i] = con->_coefficients[i];
    }
  }
  
  SCIP_CONS *scip_con;
  SCIP_CALL_EXC( SCIPcreateConsLinear(_scip, &scip_con,"constraint",
					con->_variables.size(), // # vars
					vars, // variables
					weights, // values
					con->_lhs, // LHS
					con->_rhs, // RHS
					TRUE, TRUE, TRUE, TRUE, TRUE,
					FALSE, FALSE, FALSE, FALSE, FALSE) );
  SCIP_CALL_EXC( SCIPaddCons(_scip, scip_con) ); 
  // scipcons->push_back(scip_con);

  delete[] weights;
  delete[] vars;
}

int SCIPSolver::solve(){
  DBG("solve!%s\n", "");
  
  if(!has_been_added) initialise();

  if(_verbosity > 0 && _verbosity < 3){
    // Do nothing extra
  } else if(_verbosity >= 3){
    SCIP_CALL_EXC(SCIPprintOrigProblem(_scip, NULL, NULL, FALSE));
    // SCIP_CALL_EXC(SCIPwriteOrigProblem(_scip, "scip.lp", "lp", TRUE));
  } else {
      // disable scip output to stdout
    SCIP_CALL_EXC( SCIPsetMessagehdlr(_scip, NULL) );
  }

  SCIP_CALL_EXC( SCIPsolve(_scip) );
  SCIP_STATUS status = SCIPgetStatus(_scip);
 
  if( status == SCIP_STATUS_OPTIMAL ) return SAT;
  else if( status == SCIP_STATUS_INFEASIBLE ) return UNSAT;
  else return UNKNOWN;
}

void SCIPSolver::setTimeLimit(const int cutoff){
  SCIPsetRealParam(_scip, "limits/time", (double)cutoff);
}

void SCIPSolver::setOptimalityGap(const double gap){
  if(gap < 0.0 || gap > 1.0){
    std::cerr << "Warning: relative optimality gap must be between 0.0 and 1.0, ignoring." << std::endl;
  } else {
    SCIPsetRealParam(_scip, "limits/gap", gap);
  }
}
 
void SCIPSolver::setVerbosity(const int degree){ _verbosity = degree; }

bool SCIPSolver::is_sat(){
  return !( SCIPgetStatus(_scip) == SCIP_STATUS_INFEASIBLE );
}

bool SCIPSolver::is_unsat(){ return ! is_sat(); }

bool SCIPSolver::is_opt(){
  return SCIPgetStatus(_scip) == SCIP_STATUS_OPTIMAL;
}

void SCIPSolver::output_orig_problem(const char *filename){
  // It's only safe to call SCIPwriteOrigProblem in certain cases:
  switch(SCIPgetStage(_scip)){
    case SCIP_STAGE_PROBLEM:
    case SCIP_STAGE_TRANSFORMING:
    case SCIP_STAGE_TRANSFORMED:
    case SCIP_STAGE_INITPRESOLVE:
    case SCIP_STAGE_PRESOLVING:
    case SCIP_STAGE_EXITPRESOLVE:
    case SCIP_STAGE_PRESOLVED:
    case SCIP_STAGE_INITSOLVE:
    case SCIP_STAGE_SOLVING:
    case SCIP_STAGE_SOLVED:
    case SCIP_STAGE_EXITSOLVE:
    case SCIP_STAGE_FREETRANS:
      SCIPwriteOrigProblem(_scip, filename, NULL, true);
      break;
    default:
      std::cerr << "Error: SCIP is not at a stage where it can output the original problem." << std::endl;
  }
}

void SCIPSolver::output_lp(const char *filename){
  output_orig_problem(filename);
}

void SCIPSolver::output_mps(const char *filename){
  output_orig_problem(filename);
}

void SCIPSolver::printStatistics(){
  std::cout << "\td Time: " << getTime() << "\tNodes:" << getNodes()
	    << std::endl;
  if(_verbosity > 1){
    SCIP_CALL_EXC(SCIPprintStatistics(_scip, NULL));
  }
}

int SCIPSolver::getNodes(){ return _scip->stat->nnodes; }

double SCIPSolver::getTime(){
  return SCIPclockGetTime(_scip->stat->solvingtime);
}

double SCIPSolver::getOptimalityGap(){
  return SCIPgetGap(_scip);;
}

double SCIPSolver::get_value(void *ptr){
  return SCIPgetSolVal(_scip, SCIPgetBestSol(_scip), (SCIP_VAR*) ptr);
}
