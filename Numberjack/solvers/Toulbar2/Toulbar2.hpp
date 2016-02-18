
/** \file Toulbar2.hpp
    \brief Header for the PYTHON Wrapper.
*/

#ifndef _PYTHON_H
#define _PYTHON_H

const int MAXCOST = 100000000;

#include <tb2solver.hpp>

/**
   Array of expressions (for nary constraints)
   These are used to pass in a number of variables or expressions to constraints
*/
template<class T>
class Toulbar2Array
{
private:
  std::vector< T > _array;

public:
  Toulbar2Array() {}
  virtual ~Toulbar2Array() {}
  int size() { return _array.size(); }
  void add(T arg) { _array.push_back(arg); }
  void erase(const int i) { _array.erase(_array.begin()+i); }
  T& get_item(const int i) { return _array[i]; }
  void set_item(const int i, T item) { _array[i]=item; }
};

typedef Toulbar2Array< int > Toulbar2IntArray;
typedef Toulbar2Array< double > Toulbar2DoubleArray;
typedef Toulbar2IntArray Toulbar2IntMultiArray; //SdG: multi-arrays are not handled by Numberjack load expr

/**
   Expression (Used to encode variables & constraints)
   
   Everything is an expression. 
   
*/
class Toulbar2Solver;
class Toulbar2_Expression
{
public:
  /**
   * Unique indentifier 
   */
  int nbj_ident;
  Toulbar2Solver *_solver;

  string _name;
  Value _iinf; 
  Value _isup;
  int _size;
  Value *_domain;
  int _wcspIndex;

  /**
   * Returns true if the expression has been added to the underlying solvers
   * used to ensure that things are not added into the solver twice
   */
  bool has_been_added() const;

  /**
   * Creates an expression that has a binary domain
   */
  Toulbar2_Expression();
  
  /**
   * Creates an expression that has a domain [0, nval-1] or a domain with
   * nval values
   */
  Toulbar2_Expression(const int nval);
  
  /**
   * Creates an expression that has a domain lb to ub
   */
  Toulbar2_Expression(const int lb, const int ub);
  
  /**
   * Creates an expression whose domain is specified by a list of values
   */
  Toulbar2_Expression(Toulbar2IntArray& vals);
  
  /**
   * Destructor
   */
  virtual ~Toulbar2_Expression() {}

  /**
   * Returns the identifier of this expression
   */
  int getVariableId() const {return nbj_ident;}
  
  /**
   * Returns the next value after v in the domain of the variable
   */
  int next(int v);

  /**
   * Returns the value of the variable 
   */
  int get_value() const;
  
  /**
   * Returns the current size of the domain
   */
  unsigned int get_size() const;
  
  /**
   * Returns the upper bound of the expression
   */
  int get_min() const;
  
  /**
   * Returs the lower bound of the expression
   */
  int get_max() const;
  
  /**
   * Returns true iff the value v is in the domain of the variable
   */
  bool contain(const int v) const;

  /**
   * VERY IMPORTANT FUNCTION
   *
   * This method adds the object to the underlying solver.
   *
   * @param solver :- Underlying solver wrapper to which the expression will be
   * added
   * @param top_level :- True iff the expresion is the root expression of an
   * expression tree, false otherwise (leaf or within an expression tree)
   *
   * @return The expression thar represents the node in the expression tree.
   * Normally this is the object itself apart from cases when things are encoded
   * or decomposed.
   *
   * Exery class that inherits from expression overloads this method and uses it
   * to ass itself to the underlying solver. Variable does not overload this add
   * method so if Expression::add() is called a variable should be created
   * 
   */
  virtual Toulbar2_Expression* add(Toulbar2Solver *solver, bool top_level);

};

typedef Toulbar2Array< Toulbar2_Expression* > Toulbar2ExpArray;


/**
 * This is the Integer varibale class
 */
#ifdef _DEBUGWRAP
int globalcounter = 0;
#endif
class Toulbar2_IntVar : public Toulbar2_Expression
{
  int counter;
public: 
  /**
   * Creates a binary integer variable
   */
  Toulbar2_IntVar() : Toulbar2_Expression() {
#ifdef _DEBUGWRAP
	globalcounter++;
	counter=globalcounter;
	cout << "INTVAR " << counter << endl; 
#endif
	nbj_ident = -1;
  }
  
  /**
   * Creates an integer variable object
   * whose domain is [lb, ub] with identifier 1
   */
  Toulbar2_IntVar(const int lb, const int ub, const int ident) : Toulbar2_Expression(lb, ub) {
#ifdef _DEBUGWRAP
	globalcounter++;
	counter=globalcounter;
	cout << "INTVAR " << counter << endl; 
#endif
	nbj_ident = ident;
  }
  
  /**
   * Creates an integer variable object whose domain is specified by
   * a list of integers and whose identifier is one
   */
  Toulbar2_IntVar(Toulbar2IntArray& vals, const int ident) : Toulbar2_Expression(vals) {
#ifdef _DEBUGWRAP
	globalcounter++;
	counter=globalcounter;
	cout << "INTVAR " << counter << endl; 
#endif
	nbj_ident = ident;
  }

  virtual ~Toulbar2_IntVar() {
#ifdef _DEBUGWRAP
	cout << "KILL INTVAR " << counter << endl;
#endif
  }
};


/**
 * AllDifferent constraint
 */
class Toulbar2_AllDiff : public Toulbar2_Expression
{
private:
  /**
   * Variables that are in the scope of the constraint
   */
  Toulbar2ExpArray _vars;
  int* _scope;
  string _type;
  string _semantics;
  string _baseCost;
  Cost _costValue;
  
public:
  
  /**
   * All Different constraint on an array of expressions
   */
  Toulbar2_AllDiff(Toulbar2ExpArray& vars);
  
  /**
   * All Different constraint on two expressions, equivalent to a not equal 
   */
  Toulbar2_AllDiff(Toulbar2_Expression *var1, Toulbar2_Expression *var2);

  /**
   * Add the all different constraint to the solver
   *
   * see Expression::add()
   */
  virtual Toulbar2_Expression* add(Toulbar2Solver *solver, bool top_level);
};

/**
 * GCC constraint 
 */
class Toulbar2_Gcc : public Toulbar2_Expression
{
private:
  
  /**
   * Variables in scope of constraint
   */
  Toulbar2ExpArray _vars;
  Toulbar2IntArray _vals;
  Toulbar2IntArray _lb_card;
  Toulbar2IntArray _ub_card;
  int* _scope; 
  int* _values;
  int* _lb;
  int* _ub;
  int _nbValues;

public:
  Toulbar2_Gcc(Toulbar2ExpArray& vars, Toulbar2IntArray& vals, Toulbar2IntArray& lb_card, Toulbar2IntArray& ub_card);

  /**
   * Adds the constraint into the solver
   *
   * see Expression::add()
   */
  virtual Toulbar2_Expression* add(Toulbar2Solver *solver, bool top_level);
};


class Toulbar2_PostNullary : public Toulbar2_Expression
{
private:
  Cost _cost;

public:
  
  /**
   * Create an zero-arity cost function object
   *
   * @param cost : initial cost 
   */
  Toulbar2_PostNullary(int cost);

  /**
   * Adds the constraint into the underlying solver
   *
   * see Expression::add()
   */
  virtual Toulbar2_Expression* add(Toulbar2Solver *solver, bool top_level);
};

class Toulbar2_PostUnary : public Toulbar2_Expression
{
private:
  
  /**
   * Variables in the scope of the constraint
   */
  Toulbar2_Expression* _var;
  std::vector<Cost> _costs;

public:
  
  /**
   * Create an element constraint object
   *
   * @param vars :- 
   */
  Toulbar2_PostUnary(Toulbar2_Expression* var, Toulbar2IntArray& costs);

  /**
   * Adds the constraint into the underlying solver
   *
   * see Expression::add()
   */
  virtual Toulbar2_Expression* add(Toulbar2Solver *solver, bool top_level);
};

class Toulbar2_PostBinary : public Toulbar2_Expression
{
private:
  
  /**
   * Variables in the scope of the constraint
   */
  Toulbar2_Expression* _var1;
  Toulbar2_Expression* _var2;
  std::vector<Cost> _costs;

public:
  
  /**
   * Create an element constraint object
   *
   * @param vars :- 
   */
  Toulbar2_PostBinary(Toulbar2_Expression* var1, Toulbar2_Expression* var2, Toulbar2IntArray& costs);

  /**
   * Adds the constraint into the underlying solver
   *
   * see Expression::add()
   */
  virtual Toulbar2_Expression* add(Toulbar2Solver *solver, bool top_level);
};

class Toulbar2_PostTernary : public Toulbar2_Expression
{
private:
  
  /**
   * Variables in the scope of the constraint
   */
  Toulbar2_Expression* _var1;
  Toulbar2_Expression* _var2;
  Toulbar2_Expression* _var3;
  std::vector<Cost> _costs;

public:
  
  /**
   * Create an element constraint object
   *
   * @param vars :- 
   */
  Toulbar2_PostTernary(Toulbar2ExpArray& vars, Toulbar2IntArray& costs);

  /**
   * Adds the constraint into the underlying solver
   *
   * see Expression::add()
   */
  virtual Toulbar2_Expression* add(Toulbar2Solver *solver, bool top_level);
};

class Toulbar2_PostNary : public Toulbar2_Expression
{
private:
  
  /**
   * Variables in the scope of the constraint
   */
  Toulbar2ExpArray _vars;
  int* _scope;
  int _arity;
  int _defcost; 
  Toulbar2IntMultiArray _values; 
  Toulbar2IntArray _costs;
  
public:
  
  /**
   * Create an element constraint object
   *
   * @param vars :-
   */
  Toulbar2_PostNary(Toulbar2ExpArray& vars, int arity, int _defcost, Toulbar2IntMultiArray& values, Toulbar2IntArray& costs);

  /**
   * Adds the constraint into the underlying solver
   *
   * see Expression::add()
   */
  virtual Toulbar2_Expression* add(Toulbar2Solver *solver, bool top_level);
};

class Toulbar2_PostWSum : public Toulbar2_Expression
{
private:
  
  /**
   * Variables in the scope of the constraint
   */
  Toulbar2ExpArray _vars;
  int* _scope;
  int _arity;
  string _semantics; 
  Cost _baseCost; 
  string _comparator; 
  int _rightRes;
  
public:
  
  /**
   * Create an element constraint object
   *
   * @param vars :-
   */
    //Example of parameter type conversion to Numberjack 
    //virtual void postWSum(int* scopeIndex, int arity, string semantics, Cost baseCost, string comparator, int rightRes) =0;
  Toulbar2_PostWSum(Toulbar2ExpArray& vars, int arity, const char* semantics, const char* baseCost, const char* comparator, int rightRes);

  /**
   * Adds the constraint into the underlying solver
   *
   * see Expression::add()
   */
  virtual Toulbar2_Expression* add(Toulbar2Solver *solver, bool top_level);
};

class Toulbar2_PostWVarSum : public Toulbar2_Expression
{
private:
  
  /**
   * Variables in the scope of the constraint
   */
  Toulbar2ExpArray _vars;
  int* _scope;
  int _arity;
  string _semantics; 
  Cost _baseCost; 
  string _comparator; 
  
public:
  
  /**
   * Create an element constraint object
   *
   * @param vars :-
   */
  Toulbar2_PostWVarSum(Toulbar2ExpArray& vars, int arity, const char* semantics, const char* baseCost, const char* comparator);

  /**
   * Adds the constraint into the underlying solver
   *
   * see Expression::add()
   */
  virtual Toulbar2_Expression* add(Toulbar2Solver *solver, bool top_level);
};

class Toulbar2_PostWAmong : public Toulbar2_Expression
{
private:
  
  /**
   * Variables in the scope of the constraint
   */
  Toulbar2ExpArray _vars;
  int* _scope;
  int _arity;
  string _semantics; 
  Cost _baseCost; 
  int* _values;
  int _nbValues;
  int _lb;
  int _ub;
  
public:
  
  /**
   * Create an element constraint object
   *
   * @param vars :-
   */
    //Example of parameter type conversion to Numberjack 
    //virtual void postWAmong(int* scopeIndex, int arity, string semantics, Cost baseCost, string comparator, int rightRes) =0;
  Toulbar2_PostWAmong(Toulbar2ExpArray& vars, int arity, const char* semantics, const char* baseCost, Toulbar2IntArray& values);
  Toulbar2_PostWAmong(Toulbar2ExpArray& vars, int arity, const char* semantics, const char* baseCost, Toulbar2IntArray& values, int lb, int ub);

  /**
   * Adds the constraint into the underlying solver
   *
   * see Expression::add()
   */
  virtual Toulbar2_Expression* add(Toulbar2Solver *solver, bool top_level);
};

class Toulbar2_Regular : public Toulbar2_Expression
{
private:
  
  /**
   * Variables in the scope of the constraint
   */
  Toulbar2ExpArray _vars;
  int* _scope;
  int _arity;
  int _nbStates; 
  string _type;
  string _measureCost;
  string _baseCost;
  std::vector<std::pair<int, Cost> > _initialStates; 
  Toulbar2IntArray _initialStatesV;
  std::vector<std::pair<int, Cost> > _acceptingStates;
  Toulbar2IntArray _acceptingStatesV;
  int** _transitions;
  Toulbar2IntMultiArray _transitionsV;
  std::vector<Cost> _transitionsCosts;
  
public:
  
  /**
   * Create an element constraint object
   *
   * @param vars :-
   */
    //Example of parameter type conversion to Numberjack 
    //virtual void postWRegular(int* scopeIndex, int arity, int nbStates, int** initial_States, int** accepting_States, int** Wtransitions, int* nbAll) =0; 
  Toulbar2_Regular(Toulbar2ExpArray& vars, int arity, int nbStates, Toulbar2IntArray& initialStates, Toulbar2IntArray& acceptingStates, Toulbar2IntMultiArray& transitions, 
  Toulbar2IntArray& initialCosts, Toulbar2IntArray& acceptingCosts, Toulbar2IntArray& transitionsCosts);
  
  Toulbar2_Regular(Toulbar2ExpArray& vars, int arity, int nbStates, Toulbar2IntArray& initialStates, Toulbar2IntArray& acceptingStates, Toulbar2IntMultiArray& transitions, const char* type,
  const char* measureCost, const char* baseCost);

  /**
   * Adds the constraint into the underlying solver
   *
   * see Expression::add()
   */
  virtual Toulbar2_Expression* add(Toulbar2Solver *solver, bool top_level);
};

class Toulbar2_Same : public Toulbar2_Expression
{
private:
  
  /**
   * Variables in the scope of the constraint
  */
  int* _scope;
  int _arity;
  string _type;
  string _semantics;
  Cost _baseCost;
  string _costValue;
  Toulbar2ExpArray _vars;
  
public:
  
  /**
   * Create an element constraint object
   *
   * @param vars :-
  */
  Toulbar2_Same(Toulbar2ExpArray& vars);
  Toulbar2_Same(Toulbar2ExpArray& vars, const char* type, const char* costValue);
  Toulbar2_Same(Toulbar2ExpArray& vars, const char* type, const char* semantics, const char* baseCost);

  /**
   * Adds the constraint into the underlying solver
   *
   * see Expression::add()
  */
  virtual Toulbar2_Expression* add(Toulbar2Solver *solver, bool top_level);
};

/**
 * WSameGcc constraint 
 */
class Toulbar2_PostWSameGcc : public Toulbar2_Expression
{
private:
  
  /**
   * Variables in scope of constraint
   */
  Toulbar2ExpArray _vars;
  
  /**
   * 
   */
  int* _scope;
  int _arity; 
  string _type;
  string _semantics;
  Cost _baseCost;
  int* _values;
  int* _lb;
  int* _ub;
  int _nbValues;

public:
  Toulbar2_PostWSameGcc(Toulbar2ExpArray& vars, Toulbar2IntArray& vals, Toulbar2IntArray& lb_card, Toulbar2IntArray& ub_card, const char* type, const char* semantics, const char* baseCost);

  /**
   * Adds the constraint into the solver
   *
   * see Expression::add()
   */
  virtual Toulbar2_Expression* add(Toulbar2Solver *solver, bool top_level);
};

class Toulbar2_PostWOverlap : public Toulbar2_Expression
{
private:
  
  /**
   * Variables in the scope of the constraint
   */
  Toulbar2ExpArray _vars;
  int* _scope;
  int _arity;
  string _semantics; 
  Cost _baseCost; 
  string _comparator; 
  int _rightRes;
  
public:
  
  /**
   * Create an element constraint object
   *
   * @param vars :-
   */
    //Example of parameter type conversion to Numberjack 
    //virtual void postWOverlap(int* scopeIndex, int arity, string semantics, Cost baseCost, string comparator, int rightRes) =0;
  Toulbar2_PostWOverlap(Toulbar2ExpArray& vars, int arity, const char* semantics, const char* baseCost, const char* comparator, int rightRes);

  /**
   * Adds the constraint into the underlying solver
   *
   * see Expression::add()
   */
  virtual Toulbar2_Expression* add(Toulbar2Solver *solver, bool top_level);
};


class Toulbar2_Table : public Toulbar2_Expression
{
private:
  Toulbar2ExpArray _vars;
  int* _scope;
  int _arity;
  Toulbar2IntArray _tuples;
  int _spin; 
  
public:
  Toulbar2_Table(Toulbar2ExpArray& vars, Toulbar2IntArray& tuples, const char* type); 
  Toulbar2_Table(Toulbar2_Expression *var1, Toulbar2_Expression *var2, Toulbar2IntArray& tuples, const char* type);
  
  //virtual void add(Toulbar2IntArray& tuple);
  virtual Toulbar2_Expression* add(Toulbar2Solver *solver, bool top_level);
};

// class Toulbar2_Minimise : public Toulbar2_Expression
// {
// private:
//   Toulbar2_Expression* _exp;

// public:
//   Toulbar2_Minimise(Toulbar2_Expression *var);
//   virtual Toulbar2_Expression* add(Toulbar2Solver *solver, bool top_level);
// };


/**
   The solver itself
*/
class Toulbar2Solver
{
public:
  WeightedCSPSolver * solver;
  WeightedCSP* wcsp;
  Cost upperbound;
  Cost optimum;
  Cost costshift;
  bool unsatisfiable;
  bool interrupted;
  vector<Value> solution;

  Toulbar2Solver();
  virtual ~Toulbar2Solver() {}

  // add an expression, in the case of a tree of expressions,
  // each node of the tree is added separately, depth first.
  void add(Toulbar2_Expression* arg);

  // used to initialise search on a given subset of variables
  void initialise(Toulbar2ExpArray& arg) {initialise();}
  // initialise the solver before solving (no more calls to add after this)
  void initialise();

  // solving methods
  bool propagate();
  int solve();
  int solveAndRestart(const int policy = 0, 
		      const unsigned int base = 32, 
		      const double factor = 1.3333333,
		      const double decay = 0.0,
		      const int reinit = -1);

/*
int first_decision_level;
int saved_level;
forceFiniteDomain
startNewSearch
getNextSolution
save
undo(bj)
post
deduce
get_decision_id()
branch_right()
reset
analyze_conflict
get_learnt_clause
nbClauses
get_clause
print_clause()
get_nogood_var
get_nogood_val
get_nogood_type
get_nogood_sign
get_nogood_size
addNogood
sacPreprocess
setRandomSeed
setRandomized
setAntiLex
guide
setLowerBounds
setUpperBounds
setRestartNogood
printStatistics
load_xml
load_gmpl
load_lp
output_cnf
num_vars
extract_graph
numNodes
get_degree
get_neighbor
get_feature
get_feature_name
int get_level() {return 0;}
*/

  void store_solution() {}
  bool is_opt() {return (!interrupted && !unsatisfiable && optimum < upperbound);}
  bool is_sat() {return (!unsatisfiable && (optimum < upperbound));}
  bool is_unsat() {return (unsatisfiable && !interrupted);}
  int getOptimum() {return optimum + costshift;} 
  int getBacktracks() {return solver->getNbBacktracks();}
  int getNodes() {return solver->getNbNodes();}
  int getFailures() {return solver->getNbBacktracks();}
  int getPropags() {return solver->getNbNodes();}
  double getTime() {return cpuTime() - ToulBar2::startCpuTime;}
  int getChecks() {return solver->getNbNodes();}
  void printStatistics() {}
  int getNumVariables() {return wcsp->numberOfVariables();}
  int getNumConstraints() {return wcsp->numberOfConstraints();}
  void printPython();
  
  void setHeuristic(const char* var_heuristic, const char* val_heuristic, const int rand) {}
  void setFailureLimit(const int cutoff) {}
  void setNodeLimit(const int cutoff) {}
  void setTimeLimit(const int cutoff) {timer(cutoff);}
  void setRandomized(const int degree) {}
  void setRandomSeed(const int seed) {}
  void setVerbosity(const int degree) {ToulBar2::verbose = degree - 1;}

  // toulbar2 options
  void debug(const bool debug);
  void writeSolution(const char* write);
  void showSolutions(const bool show);
  void dumpWCSP(const int level, const char* problem); 
  void nopre();
  void updateUb(const char* newUb);  

  void lds(const int maxlds);
  void restart(const long maxrestarts);
  void hbfs(const long hbfsgloballimit);
  void hbfsAlpha(const long hbfsalpha);
  void hbfsBeta(const long hbfsbeta);
  void hbfsOpenNodeLimit(const long openlimit);

  void lcLevel(const int level);
  void QueueComplexity(const bool queue);

  void allSolutions(const bool sol);
  void approximateCountingBTD(const bool aproxim);

  void binaryBranching(const bool boost);
  void staticVariableOrdering(const bool staticOrdering);
  void lastConflict(const bool last);
  void dichotomicBranching(const int dicho);
  void sortDomains(const bool sort);
  void weightedDegree(const int wDegree);
  void weightedTightness(const int wTight);
  void variableEliminationOrdering(const int order);
  void nbDecisionVars(const int nbDecision);
  void partialAssign(const char* certificate);

  void elimDegree(const int degree);
  void elimDegree_preprocessing(const int degree_prepoc);
  void elimSpaceMaxMB(const int size);
  void costfuncSeparate(const bool separate);
  void deadEndElimination(const int level);
  void preprocessTernaryRPC(const int size);
  void preprocessFunctional(const int func);
  void preprocessNary(const int maxnary);

  void btdMode(const int mode);
  void splitClusterMaxSize(const int size);
  void maxSeparatorSize(const int size);
  void minProperVarSize(const int size);
  void boostingBTD(const bool boost);
  void btdRootCluster(const int rCluster);
  void btdSubTree(const int sTree); 

  void vac(const int depth);
  void vacValueHeuristic(const bool vacVal);
  void costThreshold(const char* cost);
  void costThresholdPre(const char* cost);
  void costMultiplier(const char* cost);
  void singletonConsistency(const bool singleCons);
  void minsumDiffusion(const int min);

  void incop(const char* cmd);
};

#endif // _PYTHON_H
