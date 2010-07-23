
/** \file ExampleInterface.hpp
    \brief Header for the PYTHON Wrapper.
*/

#include <iostream>
#include <vector>

#ifndef _PYTHON_H
#define _PYTHON_H

/**
   Array of expressions (for nary constraints)
   These are used to pass in a number of variables or expressions to constraints
*/
template<class T>
class ExampleInterfaceArray
{
private:
  std::vector< T > _array;

public:
  ExampleInterfaceArray() {}
  virtual ~ExampleInterfaceArray() {}
  int size() { return _array.size(); }
  void add(T arg) { _array.push_back(arg); }
  T& get_item(const int i) { return _array[i]; }
  void set_item(const int i, T item) { _array[i]=item; }
};

typedef ExampleInterfaceArray< int > ExampleInterfaceIntArray;
typedef ExampleInterfaceArray< double > ExampleInterfaceDoubleArray;

/**
   Expression (Used to encode variables & constraints)
   
   Everything is an expression. 
   
*/
class ExampleInterfaceSolver;
class ExampleInterface_Expression
{
    
public:

  /**
   * Unique indentifier 
   */
  int nbj_ident;

  /**
   * Returns true if the expression has been added to the underlying solvers
   * used to ensure that things are not added into the solver twice
   */
  bool has_been_added() const;

  /**
   * Creates an expression that has a binary domain
   */
  ExampleInterface_Expression();
  
  /**
   * Creates an expression that has a domain [0, nval-1] or a domain with
   * nval values
   */
  ExampleInterface_Expression(const int nval);
  
  /**
   * Creates an expression that has a domain lb to ub
   */
  ExampleInterface_Expression(const int lb, const int ub);
  
  /**
   * Creates an expression whose domain is specified by a list of values
   */
  ExampleInterface_Expression(ExampleInterfaceIntArray& vals);
  
  /**
   * Destructor
   */
  virtual ~ExampleInterface_Expression();

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
  virtual ExampleInterface_Expression* add(ExampleInterfaceSolver *solver, bool top_level);

};

/**
 * This is the Integer varibale class
 */
class ExampleInterface_IntVar : public ExampleInterface_Expression
{

public:
  
  /**
   * Creates a binary integer variable
   */
  ExampleInterface_IntVar() : ExampleInterface_Expression() {}
  
  /**
   * Creates an integer variable object
   * whose domain is [lb, ub] with identifier 1
   */
  ExampleInterface_IntVar(const int lb, const int ub, const int ident) : ExampleInterface_Expression(lb, ub) {nbj_ident = ident;}
  
  /**
   * Creates an integer variable object whose domain is specified by
   * a list of integers and whose identifier is one
   */
  ExampleInterface_IntVar(ExampleInterfaceIntArray& vals, const int ident) : ExampleInterface_Expression(vals) {nbj_ident = ident;}

};

typedef ExampleInterfaceArray< ExampleInterface_Expression* > ExampleInterfaceExpArray;

/**
 * Expression that represents the minimum of a set of variables
 */
class ExampleInterface_Min : public ExampleInterface_Expression
{
private:
  /**
   * The expressions that are the scope of this constraint
   */
  ExampleInterfaceExpArray _vars;

public:
  /**
   * Creates a Min expression over an array of variables
   */
  ExampleInterface_Min(ExampleInterfaceExpArray& vars);
  
  /**
   * Creates a Min expression over two variabels
   */
  ExampleInterface_Min(ExampleInterface_Expression *var1, ExampleInterface_Expression *var2);
  
  /**
   * Descructor 
   */
  virtual ~ExampleInterface_Min();

  /**
   * Add the Min expression to the underlying solver
   *
   * See: Expression::add()
   * 
   */
  virtual ExampleInterface_Expression* add(ExampleInterfaceSolver *solver, bool top_level);
};

/**
 * Expression that represents the maximum of a set of variables
 */
class ExampleInterface_Max : public ExampleInterface_Expression
{
private:
  
  /**
   * The expressions that are the scope of this constraint
   */
  ExampleInterfaceExpArray _vars;

public:
  
  /**
   * Creates a Max expression over an array of variables
   */
  ExampleInterface_Max(ExampleInterfaceExpArray& vars);
  
  /**
   * Creates a Max expression over two variabels
   */
  ExampleInterface_Max(ExampleInterface_Expression *var1, ExampleInterface_Expression *var2);
  
  /**
   * Descructor 
   */
  virtual ~ExampleInterface_Max();

  /**
   * Add the Min expression to the underlying solver
   *
   * See: Expression::add()
   * 
   */
  virtual ExampleInterface_Expression* add(ExampleInterfaceSolver *solver, bool top_level);
};

/**
 * AllDifferent constraint
 */
class ExampleInterface_AllDiff : public ExampleInterface_Expression
{
private:
  
  /**
   * Variables that are in the scope of the constraint
   */
  ExampleInterfaceExpArray _vars;

public:
  
  /**
   * All Different constraint on an array of expressions
   */
  ExampleInterface_AllDiff(ExampleInterfaceExpArray& vars);
  
  /**
   * All Different constraint on two expressions, equivalent to a not equal 
   */
  ExampleInterface_AllDiff(ExampleInterface_Expression *var1, ExampleInterface_Expression *var2);
  
  /**
   * Destructor
   */
  virtual ~ExampleInterface_AllDiff();

  /**
   * Add the all different constraint to the solver
   *
   * see Expression::add()
   */
  virtual ExampleInterface_Expression* add(ExampleInterfaceSolver *solver, bool top_level);
};

/**
 * GCC constraint 
 */
class ExampleInterface_Gcc : public ExampleInterface_Expression
{
private:
  
  /**
   * Variables in scope of constraint
   */
  ExampleInterfaceExpArray _vars;
  
  /**
   * 
   */
  ExampleInterfaceIntArray _vals;
  ExampleInterfaceIntArray _lb_card;
  ExampleInterfaceIntArray _ub_card;

public:
  ExampleInterface_Gcc(ExampleInterfaceExpArray& vars,
	      ExampleInterfaceIntArray& vals,
	      ExampleInterfaceIntArray& lb_card,
	      ExampleInterfaceIntArray& ub_card);
  
  /**
   * Destructor
   */
  virtual ~ExampleInterface_Gcc();

  /**
   * Adds the constraint into the solver
   *
   * see Expression::add()
   */
  virtual ExampleInterface_Expression* add(ExampleInterfaceSolver *solver, bool top_level);
};

/**
 * Element constraint
 */
class ExampleInterface_Element : public ExampleInterface_Expression
{
private:
  
  /**
   * Variables in the scope of the constraint
   */
  ExampleInterfaceExpArray _vars;

public:
  
  /**
   * Create an element constraint object
   *
   * @param vars :- 
   */
  ExampleInterface_Element(ExampleInterfaceExpArray& vars);
  
  /**
   * Destructor
   */
  virtual ~ExampleInterface_Element();

  /**
   * Adds the constraint into the underlying solver
   *
   * see Expression::add()
   */
  virtual ExampleInterface_Expression* add(ExampleInterfaceSolver *solver, bool top_level);
};

/**
 * LeqLex constraint
 */
class ExampleInterface_LeqLex : public ExampleInterface_Expression
{
private:
  
  /**
   * Variables in the scope of the constriant
   */
  ExampleInterfaceExpArray _vars;

public:
  
  /**
   * Create an LeqLex constraint object
   */
  ExampleInterface_LeqLex(ExampleInterfaceExpArray& vars);
  
  /**
   * Destructor
   */
  virtual ~ExampleInterface_LeqLex();

  /**
   * Add the constraint into the solver
   *
   * see Expression::add()
   */
  virtual ExampleInterface_Expression* add(ExampleInterfaceSolver *solver, bool top_level);
};

/**
 * LessLex constraint
 */
class ExampleInterface_LessLex : public ExampleInterface_Expression
{
private:
  
  /**
   * Variables in the scope of the constraint
   */
  ExampleInterfaceExpArray _vars;

public:
  
  /**
   * Create a less lex constraint object
   */
  ExampleInterface_LessLex(ExampleInterfaceExpArray& vars);
  
  /**
   * Destructor
   */
  virtual ~ExampleInterface_LessLex();

  /**
   * Add the constraint to the underlying solver
   *
   * see Expression::add()
   */
  virtual ExampleInterface_Expression* add(ExampleInterfaceSolver *solver, bool top_level);
};

/**
 * Weighted sum expression
 */
class ExampleInterface_Sum: public ExampleInterface_Expression
{

private:

  /**
   * Variables in the scope of the constraint
   */
  ExampleInterfaceExpArray _vars;
  
  /**
   * Weights of the linear equation
   */
  ExampleInterfaceIntArray _weights;
    
public:
  
  /**
   * Create a weighted sum expression object that is the weighted sum
   * of the expression array vars and the integer array weights with a optional
   * offset
   */
  ExampleInterface_Sum( ExampleInterfaceExpArray& vars, 
	       ExampleInterfaceIntArray& weights, 
	       const int offset=0);
  
  /**
   * Create a weighted sum expression object that is the weighted sum
   * of the two expressions and the integer array weights with an offset
   */
  ExampleInterface_Sum( ExampleInterface_Expression* arg1, 
	       ExampleInterface_Expression* arg2, 
	       ExampleInterfaceIntArray& weights,
	       const int offset );
   
  /**
   * Create a weighted sum expression object that is the weighted sum
   * of the one expression and one array of expressions
   * and the integer array weights with an offset
   */ 
  ExampleInterface_Sum( ExampleInterface_Expression* arg, 
	       ExampleInterfaceIntArray& weights,
	       const int offset );
  
  /**
   * Empty sum object addVar and addWeight methods can be used to add in
   * variables and weights
   */
  ExampleInterface_Sum();

  /**
   * Add a variable to the weighted sum
   */
  void addVar( ExampleInterface_Expression* v );
  
  /**
   * Add a weight to the weighted sum
   */
  void addWeight( const int w );

  /**
   * Destructor
   */
  virtual ~ExampleInterface_Sum();
  
  /**
   * Add the weighted sum to the solver
   *
   * see Expression:add()
   *
   * @return This method needs to return an expression that will represent the
   * value of the weighted sum in the underlying solver
   * be it an intermediate variable or otherwise.
   */
  virtual ExampleInterface_Expression* add(ExampleInterfaceSolver *solver, bool top_level);
};

/**
 * Class that all binary and unary expressions inherit from
 */
class ExampleInterface_binop: public ExampleInterface_Expression
{
protected:
  
  /**
   * Two expression array that will represent the scope of the expression
   */
  ExampleInterface_Expression *_vars[2];
  
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
  ExampleInterface_binop(ExampleInterface_Expression *var1,
			 ExampleInterface_Expression *var2);
  
  /**
   *
   */
  ExampleInterface_binop(ExampleInterface_Expression *var1, int int_arg);
  
  /**
   * Destructor
   */
  virtual ~ExampleInterface_binop();

  /**
   *
   */
  virtual ExampleInterface_Expression* add(ExampleInterfaceSolver *solver, bool top_level) = 0;
};

class ExampleInterface_mul: public ExampleInterface_binop
{
public:
  ExampleInterface_mul(ExampleInterface_Expression *var1, ExampleInterface_Expression *var2);
  ExampleInterface_mul(ExampleInterface_Expression *var1, int int_arg);
  virtual ~ExampleInterface_mul();
  virtual ExampleInterface_Expression* add(ExampleInterfaceSolver *solver, bool top_level);
};

class ExampleInterface_div: public ExampleInterface_binop
{
public:
  ExampleInterface_div(ExampleInterface_Expression *var1, ExampleInterface_Expression *var2);
  ExampleInterface_div(ExampleInterface_Expression *var1, int int_arg);
  virtual ~ExampleInterface_div();
  virtual ExampleInterface_Expression* add(ExampleInterfaceSolver *solver, bool top_level);
};

class ExampleInterface_mod: public ExampleInterface_binop
{
public:
  ExampleInterface_mod(ExampleInterface_Expression *var1, ExampleInterface_Expression *var2);
  ExampleInterface_mod(ExampleInterface_Expression *var1, int int_arg);
  virtual ~ExampleInterface_mod();
  virtual ExampleInterface_Expression* add(ExampleInterfaceSolver *solver, bool top_level);
};

class ExampleInterface_and: public ExampleInterface_binop
{
public:
  ExampleInterface_and(ExampleInterface_Expression *var1, ExampleInterface_Expression *var2);
  ExampleInterface_and(ExampleInterface_Expression *var1, int int_arg);
  virtual ~ExampleInterface_and();
  virtual ExampleInterface_Expression* add(ExampleInterfaceSolver *solver, bool top_level);
};

class ExampleInterface_or: public ExampleInterface_binop
{
public:
  ExampleInterface_or(ExampleInterface_Expression *var1, ExampleInterface_Expression *var2);
  ExampleInterface_or(ExampleInterface_Expression *var1, int int_arg);
  virtual ~ExampleInterface_or();
  virtual ExampleInterface_Expression* add(ExampleInterfaceSolver *solver, bool top_level);
};

class ExampleInterface_eq: public ExampleInterface_binop
{
public:
  ExampleInterface_eq(ExampleInterface_Expression *var1, ExampleInterface_Expression *var2);
  ExampleInterface_eq(ExampleInterface_Expression *var1, int int_arg);
  virtual ~ExampleInterface_eq();
  virtual ExampleInterface_Expression* add(ExampleInterfaceSolver *solver, bool top_level);
};

class ExampleInterface_ne: public ExampleInterface_binop
{
public:
  ExampleInterface_ne(ExampleInterface_Expression *var1, ExampleInterface_Expression *var2);
  ExampleInterface_ne(ExampleInterface_Expression *var1, int int_arg);
  virtual ~ExampleInterface_ne();
  virtual ExampleInterface_Expression* add(ExampleInterfaceSolver *solver, bool top_level);
};

class ExampleInterface_le: public ExampleInterface_binop
{
public:
  ExampleInterface_le(ExampleInterface_Expression *var1, ExampleInterface_Expression *var2);
  ExampleInterface_le(ExampleInterface_Expression *var1, int int_arg);
  virtual ~ExampleInterface_le();
  virtual ExampleInterface_Expression* add(ExampleInterfaceSolver *solver, bool top_level);
};

class ExampleInterface_ge: public ExampleInterface_binop
{
public:
  ExampleInterface_ge(ExampleInterface_Expression *var1, ExampleInterface_Expression *var2);
  ExampleInterface_ge(ExampleInterface_Expression *var1, int int_arg);
  virtual ~ExampleInterface_ge();
  virtual ExampleInterface_Expression* add(ExampleInterfaceSolver *solver, bool top_level);
};

class ExampleInterface_lt: public ExampleInterface_binop
{
public:
  ExampleInterface_lt(ExampleInterface_Expression *var1, ExampleInterface_Expression *var2);
  ExampleInterface_lt(ExampleInterface_Expression *var1, int int_arg);
  virtual ~ExampleInterface_lt();
  virtual ExampleInterface_Expression* add(ExampleInterfaceSolver *solver, bool top_level);
};

class ExampleInterface_gt: public ExampleInterface_binop
{
public:
  ExampleInterface_gt(ExampleInterface_Expression *var1, ExampleInterface_Expression *var2);
  ExampleInterface_gt(ExampleInterface_Expression *var1, int int_arg);
  virtual ~ExampleInterface_gt();
  virtual ExampleInterface_Expression* add(ExampleInterfaceSolver *solver, bool top_level);
};


class ExampleInterface_Minimise : public ExampleInterface_Expression
{
private:
  ExampleInterface_Expression* _exp;

public:
  ExampleInterface_Minimise(ExampleInterface_Expression *var);
  virtual ~ExampleInterface_Minimise();
  virtual ExampleInterface_Expression* add(ExampleInterfaceSolver *solver, bool top_level);
};


class ExampleInterface_Maximise : public ExampleInterface_Expression
{
private:
  ExampleInterface_Expression* _exp;

public:
  ExampleInterface_Maximise(ExampleInterface_Expression *var);
  virtual ~ExampleInterface_Maximise();
  virtual ExampleInterface_Expression* add(ExampleInterfaceSolver *solver, bool top_level);
};


/**
   The solver itself
*/
class ExampleInterfaceSolver
{

public:

  int first_decision_level;
  int saved_level;

  ExampleInterfaceSolver();
  virtual ~ExampleInterfaceSolver();

  // add an expression, in the case of a tree of expressions,
  // each node of the tree is added separately, depth first.
  void add(ExampleInterface_Expression* arg);

  // used to initialise search on a given subset of variables
  void initialise(ExampleInterfaceExpArray& arg);
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

  int next(ExampleInterface_Expression* x, int v);

  int get_level();

  bool propagate();
  void save();
  void post(const char* op, ExampleInterface_Expression* x, int v);
  bool undo(const int nlevel);
  void deduce(const char* op, ExampleInterface_Expression* x, int v);
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
  void forceFiniteDomain(ExampleInterfaceExpArray& vars);
  void addNogood(ExampleInterfaceExpArray& vars, 
		 ExampleInterfaceIntArray& vals);
  void guide(ExampleInterfaceExpArray& vars, 
	     ExampleInterfaceIntArray& vals,
	     ExampleInterfaceDoubleArray& probs);
  void backtrackTo(const int level);
  void upOneLevel();
  void presolve();
  void assign(ExampleInterface_Expression *X, const int v);
  void increase_init_level(const int i);
  void decrease_init_level(const int i);
  void reset(bool full);
  void setLowerBounds(ExampleInterfaceExpArray& vars, 
		      ExampleInterfaceIntArray& vals);
  void setUpperBounds(ExampleInterfaceExpArray& vars, 
		      ExampleInterfaceIntArray& vals);
  void setRestartNogood();

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
