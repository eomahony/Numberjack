
/** \file Mistral2.hpp
    \brief Header for the PYTHON Wrapper.
*/

#include <iostream>
#include <vector>

#include <mistral_solver.hpp>
#include <mistral_search.hpp>


#ifndef _PYTHON_H
#define _PYTHON_H

/**
   Array of expressions (for nary constraints)
   These are used to pass in a number of variables or expressions to constraints
*/
template<class T>
class Mistral2Array
{
private:
  std::vector< T > _array;

public:
  Mistral2Array() {}
  virtual ~Mistral2Array() {}
  int size() { return _array.size(); }
  void add(T arg) { _array.push_back(arg); }
  T& get_item(const int i) { return _array[i]; }
  void set_item(const int i, T item) { _array[i]=item; }
};

typedef Mistral2Array< int > Mistral2IntArray;
typedef Mistral2Array< double > Mistral2DoubleArray;

/**
   Expression (Used to encode variables & constraints)
   
   Everything is an expression. 
   
*/
class Mistral2Solver;
class Mistral2_Expression
{
    
public:

  /**
   * Unique indentifier 
   */
  int nbj_ident;
  Mistral::Variable _self;
  Mistral2Solver *_solver;

  /**
   * Returns true if the expression has been added to the underlying solvers
   * used to ensure that things are not added into the solver twice
   */
  bool has_been_added() const;

  /**
   * Creates an expression that has a binary domain
   */
  Mistral2_Expression();
  
  /**
   * Creates an expression that has a domain [0, nval-1] or a domain with
   * nval values
   */
  Mistral2_Expression(const int nval);
  
  /**
   * Creates an expression that has a domain lb to ub
   */
  Mistral2_Expression(const int lb, const int ub);
  
  /**
   * Creates an expression whose domain is specified by a list of values
   */
  Mistral2_Expression(Mistral2IntArray& vals);
  
  /**
   * Destructor
   */
  virtual ~Mistral2_Expression();

  /**
   * Returns the identifier of this expression
   */
  int getVariableId() const;
  
  /**
   * Returns the value of the variable 
   */
  int get_value() const;
  
  /**
   * Returns the current size of the domain
   */
  int get_size() const;
  
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
  virtual Mistral2_Expression* add(Mistral2Solver *solver, bool top_level);

};

/**
 * This is the Integer varibale class
 */
class Mistral2_IntVar : public Mistral2_Expression
{

public:
  
  /**
   * Creates a binary integer variable
   */
  Mistral2_IntVar() : Mistral2_Expression() {}
  
  /**
   * Creates an integer variable object
   * whose domain is [lb, ub] with identifier 1
   */
  Mistral2_IntVar(const int lb, const int ub, const int ident) : Mistral2_Expression(lb, ub) {nbj_ident = ident;}
  
  /**
   * Creates an integer variable object whose domain is specified by
   * a list of integers and whose identifier is one
   */
  Mistral2_IntVar(Mistral2IntArray& vals, const int ident) : Mistral2_Expression(vals) {nbj_ident = ident;}

};

typedef Mistral2Array< Mistral2_Expression* > Mistral2ExpArray;

/**
 * Expression that represents the minimum of a set of variables
 */
class Mistral2_Min : public Mistral2_Expression
{
private:
  /**
   * The expressions that are the scope of this constraint
   */
  Mistral2ExpArray _vars;

public:
  /**
   * Creates a Min expression over an array of variables
   */
  Mistral2_Min(Mistral2ExpArray& vars);
  
  /**
   * Creates a Min expression over two variabels
   */
  Mistral2_Min(Mistral2_Expression *var1, Mistral2_Expression *var2);
  
  /**
   * Descructor 
   */
  virtual ~Mistral2_Min();

  /**
   * Add the Min expression to the underlying solver
   *
   * See: Expression::add()
   * 
   */
  virtual Mistral2_Expression* add(Mistral2Solver *solver, bool top_level);
};

/**
 * Expression that represents the maximum of a set of variables
 */
class Mistral2_Max : public Mistral2_Expression
{
private:
  
  /**
   * The expressions that are the scope of this constraint
   */
  Mistral2ExpArray _vars;

public:
  
  /**
   * Creates a Max expression over an array of variables
   */
  Mistral2_Max(Mistral2ExpArray& vars);
  
  /**
   * Creates a Max expression over two variabels
   */
  Mistral2_Max(Mistral2_Expression *var1, Mistral2_Expression *var2);
  
  /**
   * Descructor 
   */
  virtual ~Mistral2_Max();

  /**
   * Add the Min expression to the underlying solver
   *
   * See: Expression::add()
   * 
   */
  virtual Mistral2_Expression* add(Mistral2Solver *solver, bool top_level);
};

/**
 * AllDifferent constraint
 */
class Mistral2_AllDiff : public Mistral2_Expression
{
private:
  
  /**
   * Variables that are in the scope of the constraint
   */
  Mistral2ExpArray _vars;

public:
  
  /**
   * All Different constraint on an array of expressions
   */
  Mistral2_AllDiff(Mistral2ExpArray& vars);
  
  /**
   * All Different constraint on two expressions, equivalent to a not equal 
   */
  Mistral2_AllDiff(Mistral2_Expression *var1, Mistral2_Expression *var2);
  
  /**
   * Destructor
   */
  virtual ~Mistral2_AllDiff();

  /**
   * Add the all different constraint to the solver
   *
   * see Expression::add()
   */
  virtual Mistral2_Expression* add(Mistral2Solver *solver, bool top_level);
};

/**
 * GCC constraint 
 */
class Mistral2_Gcc : public Mistral2_Expression
{
private:
  
  /**
   * Variables in scope of constraint
   */
  Mistral2ExpArray _vars;
  
  /**
   * 
   */
  Mistral2IntArray _vals;
  Mistral2IntArray _lb_card;
  Mistral2IntArray _ub_card;

public:
  Mistral2_Gcc(Mistral2ExpArray& vars,
	      Mistral2IntArray& vals,
	      Mistral2IntArray& lb_card,
	      Mistral2IntArray& ub_card);
  
  /**
   * Destructor
   */
  virtual ~Mistral2_Gcc();

  /**
   * Adds the constraint into the solver
   *
   * see Expression::add()
   */
  virtual Mistral2_Expression* add(Mistral2Solver *solver, bool top_level);
};

/**
 * Element constraint
 */
class Mistral2_Element : public Mistral2_Expression
{
private:
  
  /**
   * Variables in the scope of the constraint
   */
  Mistral2ExpArray _vars;

public:
  
  /**
   * Create an element constraint object
   *
   * @param vars :- 
   */
  Mistral2_Element(Mistral2ExpArray& vars);
  
  /**
   * Destructor
   */
  virtual ~Mistral2_Element();

  /**
   * Adds the constraint into the underlying solver
   *
   * see Expression::add()
   */
  virtual Mistral2_Expression* add(Mistral2Solver *solver, bool top_level);
};

/**
 * LeqLex constraint
 */
class Mistral2_LeqLex : public Mistral2_Expression
{
private:
  
  /**
   * Variables in the scope of the constriant
   */
  Mistral2ExpArray _vars;

public:
  
  /**
   * Create an LeqLex constraint object
   */
  Mistral2_LeqLex(Mistral2ExpArray& vars);
  
  /**
   * Destructor
   */
  virtual ~Mistral2_LeqLex();

  /**
   * Add the constraint into the solver
   *
   * see Expression::add()
   */
  virtual Mistral2_Expression* add(Mistral2Solver *solver, bool top_level);
};

/**
 * LessLex constraint
 */
class Mistral2_LessLex : public Mistral2_Expression
{
private:
  
  /**
   * Variables in the scope of the constraint
   */
  Mistral2ExpArray _vars;

public:
  
  /**
   * Create a less lex constraint object
   */
  Mistral2_LessLex(Mistral2ExpArray& vars);
  
  /**
   * Destructor
   */
  virtual ~Mistral2_LessLex();

  /**
   * Add the constraint to the underlying solver
   *
   * see Expression::add()
   */
  virtual Mistral2_Expression* add(Mistral2Solver *solver, bool top_level);
};

/**
 * Weighted sum expression
 */
class Mistral2_Sum: public Mistral2_Expression
{

private:

  /**
   * Variables in the scope of the constraint
   */
  Mistral2ExpArray _vars;
  
  /**
   * Weights of the linear equation
   */
  Mistral2IntArray _weights;
    
  int _offset;


public:
  
  /**
   * Create a weighted sum expression object that is the weighted sum
   * of the expression array vars and the integer array weights with a optional
   * offset
   */
  Mistral2_Sum( Mistral2ExpArray& vars, 
	       Mistral2IntArray& weights, 
	       const int offset=0);
  
  /**
   * Create a weighted sum expression object that is the weighted sum
   * of the two expressions and the integer array weights with an offset
   */
  Mistral2_Sum( Mistral2_Expression* arg1, 
	       Mistral2_Expression* arg2, 
	       Mistral2IntArray& weights,
	       const int offset );
   
  /**
   * Create a weighted sum expression object that is the weighted sum
   * of the one expression and one array of expressions
   * and the integer array weights with an offset
   */ 
  Mistral2_Sum( Mistral2_Expression* arg, 
	       Mistral2IntArray& weights,
	       const int offset );
  
  /**
   * Empty sum object addVar and addWeight methods can be used to add in
   * variables and weights
   */
  Mistral2_Sum();

  /**
   * Add a variable to the weighted sum
   */
  void addVar( Mistral2_Expression* v );
  
  /**
   * Add a weight to the weighted sum
   */
  void addWeight( const int w );

  /**
   * Destructor
   */
  virtual ~Mistral2_Sum();
  
  /**
   * Add the weighted sum to the solver
   *
   * see Expression:add()
   *
   * @return This method needs to return an expression that will represent the
   * value of the weighted sum in the underlying solver
   * be it an intermediate variable or otherwise.
   */
  virtual Mistral2_Expression* add(Mistral2Solver *solver, bool top_level);
};


/**
 * Ordered sum expression
 */
class Mistral2_OrderedSum: public Mistral2_Expression
{

private:

  /**
   * Variables in the scope of the constraint
   */
  Mistral2ExpArray _vars;
  
  int _lb;
	int _ub;


public:
  
  /**
   * Create an ordered sum expression object 
   */
  Mistral2_OrderedSum( Mistral2ExpArray& vars, 
	       const int l, const int u);

  /**
   * Destructor
   */
  virtual ~Mistral2_OrderedSum();
  
  /**
   * Add the ordered sum to the solver
   *
   * see Expression:add()
   *
   * @return This method needs to return an expression that will represent the
   * value of the ordered sum in the underlying solver
   * be it an intermediate variable or otherwise. [todo, currently only useable at top level]
   */
  virtual Mistral2_Expression* add(Mistral2Solver *solver, bool top_level);
};


/**
 * Class that all binary and unary expressions inherit from
 */
class Mistral2_binop: public Mistral2_Expression
{
protected:
  
  /**
   * Two expression array that will represent the scope of the expression
   */
  Mistral2_Expression *_vars[2];
  
  /**
   *
   */
  int _constant;
  
public:
  
  /**
   * @return 1 iff the expression is unary, 2 otherwise
   */
  int arity() { return 1+(_vars[1]==NULL); }
  
  /**
   *
   */
  Mistral2_binop(Mistral2_Expression *var1,
			 Mistral2_Expression *var2);
  
  /**
   *
   */
  Mistral2_binop(Mistral2_Expression *var1, int int_arg);
  
  /**
   * Destructor
   */
  virtual ~Mistral2_binop();

  /**
   *
   */
  virtual Mistral2_Expression* add(Mistral2Solver *solver, bool top_level) = 0;
};

class Mistral2_Abs: public Mistral2_binop
{
public:
  Mistral2_Abs(Mistral2_Expression *var1);
  virtual ~Mistral2_Abs();
  virtual Mistral2_Expression* add(Mistral2Solver *solver, bool top_level);
};

class Mistral2_neg: public Mistral2_binop
{
public:
  Mistral2_neg(Mistral2_Expression *var1);
  virtual ~Mistral2_neg();
  virtual Mistral2_Expression* add(Mistral2Solver *solver, bool top_level);
};


class Mistral2_mul: public Mistral2_binop
{
public:
  Mistral2_mul(Mistral2_Expression *var1, Mistral2_Expression *var2);
  Mistral2_mul(Mistral2_Expression *var1, int int_arg);
  virtual ~Mistral2_mul();
  virtual Mistral2_Expression* add(Mistral2Solver *solver, bool top_level);
};

class Mistral2_div: public Mistral2_binop
{
public:
  Mistral2_div(Mistral2_Expression *var1, Mistral2_Expression *var2);
  Mistral2_div(Mistral2_Expression *var1, int int_arg);
  virtual ~Mistral2_div();
  virtual Mistral2_Expression* add(Mistral2Solver *solver, bool top_level);
};

class Mistral2_mod: public Mistral2_binop
{
public:
  Mistral2_mod(Mistral2_Expression *var1, Mistral2_Expression *var2);
  Mistral2_mod(Mistral2_Expression *var1, int int_arg);
  virtual ~Mistral2_mod();
  virtual Mistral2_Expression* add(Mistral2Solver *solver, bool top_level);
};

class Mistral2_and: public Mistral2_binop
{
public:
  Mistral2_and(Mistral2_Expression *var1, Mistral2_Expression *var2);
  Mistral2_and(Mistral2_Expression *var1, int int_arg);
  virtual ~Mistral2_and();
  virtual Mistral2_Expression* add(Mistral2Solver *solver, bool top_level);
};

class Mistral2_or: public Mistral2_binop
{
public:
  Mistral2_or(Mistral2_Expression *var1, Mistral2_Expression *var2);
  Mistral2_or(Mistral2_Expression *var1, int int_arg);
  virtual ~Mistral2_or();
  virtual Mistral2_Expression* add(Mistral2Solver *solver, bool top_level);
};

class Mistral2_eq: public Mistral2_binop
{
public:
  Mistral2_eq(Mistral2_Expression *var1, Mistral2_Expression *var2);
  Mistral2_eq(Mistral2_Expression *var1, int int_arg);
  virtual ~Mistral2_eq();
  virtual Mistral2_Expression* add(Mistral2Solver *solver, bool top_level);
};

class Mistral2_ne: public Mistral2_binop
{
public:
  Mistral2_ne(Mistral2_Expression *var1, Mistral2_Expression *var2);
  Mistral2_ne(Mistral2_Expression *var1, int int_arg);
  virtual ~Mistral2_ne();
  virtual Mistral2_Expression* add(Mistral2Solver *solver, bool top_level);
};

class Mistral2_NoOverlap: public Mistral2_binop
{
private:
  int _bonstant;
public:
  Mistral2_NoOverlap(Mistral2_Expression *var1, Mistral2_Expression *var2, int d1, int d2);
  virtual ~Mistral2_NoOverlap();
  virtual Mistral2_Expression* add(Mistral2Solver *solver, bool top_level);
};

class Mistral2_le: public Mistral2_binop
{
public:
  Mistral2_le(Mistral2_Expression *var1, Mistral2_Expression *var2);
  Mistral2_le(Mistral2_Expression *var1, int int_arg);
  virtual ~Mistral2_le();
  virtual Mistral2_Expression* add(Mistral2Solver *solver, bool top_level);
};

class Mistral2_ge: public Mistral2_binop
{
public:
  Mistral2_ge(Mistral2_Expression *var1, Mistral2_Expression *var2);
  Mistral2_ge(Mistral2_Expression *var1, int int_arg);
  virtual ~Mistral2_ge();
  virtual Mistral2_Expression* add(Mistral2Solver *solver, bool top_level);
};

class Mistral2_lt: public Mistral2_binop
{
public:
  Mistral2_lt(Mistral2_Expression *var1, Mistral2_Expression *var2);
  Mistral2_lt(Mistral2_Expression *var1, int int_arg);
  virtual ~Mistral2_lt();
  virtual Mistral2_Expression* add(Mistral2Solver *solver, bool top_level);
};

class Mistral2_gt: public Mistral2_binop
{
public:
  Mistral2_gt(Mistral2_Expression *var1, Mistral2_Expression *var2);
  Mistral2_gt(Mistral2_Expression *var1, int int_arg);
  virtual ~Mistral2_gt();
  virtual Mistral2_Expression* add(Mistral2Solver *solver, bool top_level);
};


class Mistral2_Minimise : public Mistral2_Expression
{
private:
  Mistral2_Expression* _exp;

public:
  Mistral2_Minimise(Mistral2_Expression *var);
  virtual ~Mistral2_Minimise();
  virtual Mistral2_Expression* add(Mistral2Solver *solver, bool top_level);
};


class Mistral2_Maximise : public Mistral2_Expression
{
private:
  Mistral2_Expression* _exp;

public:
  Mistral2_Maximise(Mistral2_Expression *var);
  virtual ~Mistral2_Maximise();
  virtual Mistral2_Expression* add(Mistral2Solver *solver, bool top_level);
};


/**
   The solver itself
*/
class Mistral2Solver
{

public:

  Mistral::Solver *solver;

  int first_decision_level;
  int saved_level;


  int _heuristic_randomization;
  std::string _var_heuristic_str;
  std::string _val_heuristic_str;
  std::string _restart_policy_str;

  Mistral::BranchingHeuristic *_branching_heuristic;
  Mistral::RestartPolicy *_restart_policy;
  Mistral::Goal *_search_goal;

  Mistral2Solver();
  virtual ~Mistral2Solver();

  // add an expression, in the case of a tree of expressions,
  // each node of the tree is added separately, depth first.
  void add(Mistral2_Expression* arg);

  // used to initialise search on a given subset of variables
  void initialise(Mistral2ExpArray& arg);
  // initialise the solver before solving (no more calls to add after this)
  void initialise();

  // solving methods
  int solve();
  int solveAndRestart(const int policy = 0, 
		      const unsigned int base = 32, 
		      const double factor = 1.3333333,
		      const double decay = 0.0,
		      const int reinit = -1);
  int startNewSearch();
  int getNextSolution();
  int sacPreprocess(const int type);

  int next(Mistral2_Expression* x, int v);

  int get_level();

  bool propagate();
  void save();
  void post(const char* op, Mistral2_Expression* x, int v);
  bool undo(const int nlevel);
  void deduce(const char* op, Mistral2_Expression* x, int v);
  void deduce();
  bool branch_right();
  void store_solution();

  // parameter tuning methods
  void setHeuristic(const char* var_heuristic, const char* val_heuristic, const int rand);
  void setFailureLimit(const int cutoff);  
  void setNodeLimit(const int cutoff);  
  void setTimeLimit(const int cutoff);
  void setVerbosity(const int degree);
  void setRandomized(const int degree);
  void setRandomSeed(const int seed);
  void forceFiniteDomain(Mistral2ExpArray& vars);
  void addNogood(Mistral2ExpArray& vars, 
		 Mistral2IntArray& vals);
  void guide(Mistral2ExpArray& vars, 
	     Mistral2IntArray& vals,
	     Mistral2DoubleArray& probs);
  void backtrackTo(const int level);
  void upOneLevel();
  void presolve();
  void assign(Mistral2_Expression *X, const int v);
  void increase_init_level(const int i);
  void decrease_init_level(const int i);
  void reset(bool full);
  void setLowerBounds(Mistral2ExpArray& vars, 
		      Mistral2IntArray& vals);
  void setUpperBounds(Mistral2ExpArray& vars, 
		      Mistral2IntArray& vals);
  void setRestartNogood();
	
	void printPython();

  // statistics methods
  bool is_opt();
  bool is_sat();
  bool is_unsat();
  void printStatistics();
  int getBacktracks();
  int getNodes();
  int getFailures();
  int getChecks();
  int getPropags();
  double getTime();

  int getRandomNumber();

  int getNumVariables();
  int getNumConstraints();

  int num_vars();
  int get_degree(int i);
};

#endif // _PYTHON_H
