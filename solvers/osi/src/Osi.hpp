
/** \file Osi.hpp
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
class OsiArray
{
private:
  std::vector< T > _array;

public:
  OsiArray() {}
  virtual ~OsiArray() {}
  int size() { return _array.size(); }
  void add(T arg) { _array.push_back(arg); }
  T& get_item(const int i) { return _array[i]; }
  void set_item(const int i, T item) { _array[i]=item; }
};

typedef OsiArray< int > OsiIntArray;
typedef OsiArray< double > OsiDoubleArray;

/**
   Expression (Used to encode variables & constraints)
   
   Everything is an expression. 
   
*/
class OsiSolver;
class Osi_Expression
{
    
public:

  /**
   * Unique indentifier 
   */
  int nbj_ident;
  std::vector<double> variable_ubs;
  std::vector<double> variable_lbs;

  /**
   * Returns true if the expression has been added to the underlying solvers
   * used to ensure that things are not added into the solver twice
   */
  bool has_been_added() const;
  void initialise();

  /**
   * Creates an expression that has a binary domain
   */
  Osi_Expression();
  
  /**
   * Creates an expression that has a domain [0, nval-1] or a domain with
   * nval values
   */
  Osi_Expression(const int nval);
  
  /**
   * Creates an expression that has a domain lb to ub
   */
  Osi_Expression(const int lb, const int ub);
  
  /**
   * Creates an expression whose domain is specified by a list of values
   */
  Osi_Expression(OsiIntArray& vals);
  
  /**
   * Destructor
   */
  virtual ~Osi_Expression();

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
  virtual Osi_Expression* add(OsiSolver *solver, bool top_level);

};

/**
 * This is the Integer varibale class
 */
class Osi_IntVar : public Osi_Expression
{

public:
  
  /**
   * Creates a binary integer variable
   */
  Osi_IntVar() : Osi_Expression() {}
  
  /**
   * Creates an integer variable object
   * whose domain is [lb, ub] with identifier 1
   */
  Osi_IntVar(const int lb, const int ub, const int ident) : Osi_Expression(lb, ub) {nbj_ident = ident;}
  
  /**
   * Creates an integer variable object whose domain is specified by
   * a list of integers and whose identifier is one
   */
  Osi_IntVar(OsiIntArray& vals, const int ident) : Osi_Expression(vals) {nbj_ident = ident;}

};

typedef OsiArray< Osi_Expression* > OsiExpArray;

/**
 * Expression that represents the minimum of a set of variables
 */
class Osi_Min : public Osi_Expression
{
private:
  /**
   * The expressions that are the scope of this constraint
   */
  OsiExpArray _vars;

public:
  /**
   * Creates a Min expression over an array of variables
   */
  Osi_Min(OsiExpArray& vars);
  
  /**
   * Creates a Min expression over two variabels
   */
  Osi_Min(Osi_Expression *var1, Osi_Expression *var2);
  
  /**
   * Descructor 
   */
  virtual ~Osi_Min();

  /**
   * Add the Min expression to the underlying solver
   *
   * See: Expression::add()
   * 
   */
  virtual Osi_Expression* add(OsiSolver *solver, bool top_level);
};

/**
 * Expression that represents the maximum of a set of variables
 */
class Osi_Max : public Osi_Expression
{
private:
  
  /**
   * The expressions that are the scope of this constraint
   */
  OsiExpArray _vars;

public:
  
  /**
   * Creates a Max expression over an array of variables
   */
  Osi_Max(OsiExpArray& vars);
  
  /**
   * Creates a Max expression over two variabels
   */
  Osi_Max(Osi_Expression *var1, Osi_Expression *var2);
  
  /**
   * Descructor 
   */
  virtual ~Osi_Max();

  /**
   * Add the Min expression to the underlying solver
   *
   * See: Expression::add()
   * 
   */
  virtual Osi_Expression* add(OsiSolver *solver, bool top_level);
};

/**
 * AllDifferent constraint
 */
class Osi_AllDiff : public Osi_Expression
{
private:
  
  /**
   * Variables that are in the scope of the constraint
   */
  OsiExpArray _vars;

public:
  
  /**
   * All Different constraint on an array of expressions
   */
  Osi_AllDiff(OsiExpArray& vars);
  
  /**
   * All Different constraint on two expressions, equivalent to a not equal 
   */
  Osi_AllDiff(Osi_Expression *var1, Osi_Expression *var2);
  
  /**
   * Destructor
   */
  virtual ~Osi_AllDiff();

  /**
   * Add the all different constraint to the solver
   *
   * see Expression::add()
   */
  virtual Osi_Expression* add(OsiSolver *solver, bool top_level);
};

/**
 * GCC constraint 
 */
class Osi_Gcc : public Osi_Expression
{
private:
  
  /**
   * Variables in scope of constraint
   */
  OsiExpArray _vars;
  
  /**
   * 
   */
  OsiIntArray _vals;
  OsiIntArray _lb_card;
  OsiIntArray _ub_card;

public:
  Osi_Gcc(OsiExpArray& vars,
	      OsiIntArray& vals,
	      OsiIntArray& lb_card,
	      OsiIntArray& ub_card);
  
  /**
   * Destructor
   */
  virtual ~Osi_Gcc();

  /**
   * Adds the constraint into the solver
   *
   * see Expression::add()
   */
  virtual Osi_Expression* add(OsiSolver *solver, bool top_level);
};

/**
 * Element constraint
 */
class Osi_Element : public Osi_Expression
{
private:
  
  /**
   * Variables in the scope of the constraint
   */
  OsiExpArray _vars;

public:
  
  /**
   * Create an element constraint object
   *
   * @param vars :- 
   */
  Osi_Element(OsiExpArray& vars);
  
  /**
   * Destructor
   */
  virtual ~Osi_Element();

  /**
   * Adds the constraint into the underlying solver
   *
   * see Expression::add()
   */
  virtual Osi_Expression* add(OsiSolver *solver, bool top_level);
};

/**
 * LeqLex constraint
 */
class Osi_LeqLex : public Osi_Expression
{
private:
  
  /**
   * Variables in the scope of the constriant
   */
  OsiExpArray _vars;

public:
  
  /**
   * Create an LeqLex constraint object
   */
  Osi_LeqLex(OsiExpArray& vars);
  
  /**
   * Destructor
   */
  virtual ~Osi_LeqLex();

  /**
   * Add the constraint into the solver
   *
   * see Expression::add()
   */
  virtual Osi_Expression* add(OsiSolver *solver, bool top_level);
};

/**
 * LessLex constraint
 */
class Osi_LessLex : public Osi_Expression
{
private:
  
  /**
   * Variables in the scope of the constraint
   */
  OsiExpArray _vars;

public:
  
  /**
   * Create a less lex constraint object
   */
  Osi_LessLex(OsiExpArray& vars);
  
  /**
   * Destructor
   */
  virtual ~Osi_LessLex();

  /**
   * Add the constraint to the underlying solver
   *
   * see Expression::add()
   */
  virtual Osi_Expression* add(OsiSolver *solver, bool top_level);
};

/**
 * Weighted sum expression
 */
class Osi_Sum: public Osi_Expression
{

private:

  /**
   * Variables in the scope of the constraint
   */
  OsiExpArray _vars;
  
  /**
   * Weights of the linear equation
   */
  OsiIntArray _weights;
    
public:
  
  /**
   * Create a weighted sum expression object that is the weighted sum
   * of the expression array vars and the integer array weights with a optional
   * offset
   */
  Osi_Sum( OsiExpArray& vars, 
	       OsiIntArray& weights, 
	       const int offset=0);
  
  /**
   * Create a weighted sum expression object that is the weighted sum
   * of the two expressions and the integer array weights with an offset
   */
  Osi_Sum( Osi_Expression* arg1, 
	       Osi_Expression* arg2, 
	       OsiIntArray& weights,
	       const int offset );
   
  /**
   * Create a weighted sum expression object that is the weighted sum
   * of the one expression and one array of expressions
   * and the integer array weights with an offset
   */ 
  Osi_Sum( Osi_Expression* arg, 
	       OsiIntArray& weights,
	       const int offset );
  
  /**
   * Empty sum object addVar and addWeight methods can be used to add in
   * variables and weights
   */
  Osi_Sum();

  /**
   * Add a variable to the weighted sum
   */
  void addVar( Osi_Expression* v );
  
  /**
   * Add a weight to the weighted sum
   */
  void addWeight( const int w );

  /**
   * Destructor
   */
  virtual ~Osi_Sum();
  
  /**
   * Add the weighted sum to the solver
   *
   * see Expression:add()
   *
   * @return This method needs to return an expression that will represent the
   * value of the weighted sum in the underlying solver
   * be it an intermediate variable or otherwise.
   */
  virtual Osi_Expression* add(OsiSolver *solver, bool top_level);
};

/**
 * Class that all binary and unary expressions inherit from
 */
class Osi_binop: public Osi_Expression
{
protected:
  
  /**
   * Two expression array that will represent the scope of the expression
   */
  Osi_Expression *_vars[2];
  
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
  Osi_binop(Osi_Expression *var1,
			 Osi_Expression *var2);
  
  /**
   *
   */
  Osi_binop(Osi_Expression *var1, int int_arg);
  
  /**
   * Destructor
   */
  virtual ~Osi_binop();

  /**
   *
   */
  virtual Osi_Expression* add(OsiSolver *solver, bool top_level) = 0;
};

class Osi_mul: public Osi_binop
{
public:
  Osi_mul(Osi_Expression *var1, Osi_Expression *var2);
  Osi_mul(Osi_Expression *var1, int int_arg);
  virtual ~Osi_mul();
  virtual Osi_Expression* add(OsiSolver *solver, bool top_level);
};

class Osi_div: public Osi_binop
{
public:
  Osi_div(Osi_Expression *var1, Osi_Expression *var2);
  Osi_div(Osi_Expression *var1, int int_arg);
  virtual ~Osi_div();
  virtual Osi_Expression* add(OsiSolver *solver, bool top_level);
};

class Osi_mod: public Osi_binop
{
public:
  Osi_mod(Osi_Expression *var1, Osi_Expression *var2);
  Osi_mod(Osi_Expression *var1, int int_arg);
  virtual ~Osi_mod();
  virtual Osi_Expression* add(OsiSolver *solver, bool top_level);
};

class Osi_and: public Osi_binop
{
public:
  Osi_and(Osi_Expression *var1, Osi_Expression *var2);
  Osi_and(Osi_Expression *var1, int int_arg);
  virtual ~Osi_and();
  virtual Osi_Expression* add(OsiSolver *solver, bool top_level);
};

class Osi_or: public Osi_binop
{
public:
  Osi_or(Osi_Expression *var1, Osi_Expression *var2);
  Osi_or(Osi_Expression *var1, int int_arg);
  virtual ~Osi_or();
  virtual Osi_Expression* add(OsiSolver *solver, bool top_level);
};

class Osi_eq: public Osi_binop
{
public:
  Osi_eq(Osi_Expression *var1, Osi_Expression *var2);
  Osi_eq(Osi_Expression *var1, int int_arg);
  virtual ~Osi_eq();
  virtual Osi_Expression* add(OsiSolver *solver, bool top_level);
};

class Osi_ne: public Osi_binop
{
public:
  Osi_ne(Osi_Expression *var1, Osi_Expression *var2);
  Osi_ne(Osi_Expression *var1, int int_arg);
  virtual ~Osi_ne();
  virtual Osi_Expression* add(OsiSolver *solver, bool top_level);
};

class Osi_le: public Osi_binop
{
public:
  Osi_le(Osi_Expression *var1, Osi_Expression *var2);
  Osi_le(Osi_Expression *var1, int int_arg);
  virtual ~Osi_le();
  virtual Osi_Expression* add(OsiSolver *solver, bool top_level);
};

class Osi_ge: public Osi_binop
{
public:
  Osi_ge(Osi_Expression *var1, Osi_Expression *var2);
  Osi_ge(Osi_Expression *var1, int int_arg);
  virtual ~Osi_ge();
  virtual Osi_Expression* add(OsiSolver *solver, bool top_level);
};

class Osi_lt: public Osi_binop
{
public:
  Osi_lt(Osi_Expression *var1, Osi_Expression *var2);
  Osi_lt(Osi_Expression *var1, int int_arg);
  virtual ~Osi_lt();
  virtual Osi_Expression* add(OsiSolver *solver, bool top_level);
};

class Osi_gt: public Osi_binop
{
public:
  Osi_gt(Osi_Expression *var1, Osi_Expression *var2);
  Osi_gt(Osi_Expression *var1, int int_arg);
  virtual ~Osi_gt();
  virtual Osi_Expression* add(OsiSolver *solver, bool top_level);
};


class Osi_Minimise : public Osi_Expression
{
private:
  Osi_Expression* _exp;

public:
  Osi_Minimise(Osi_Expression *var);
  virtual ~Osi_Minimise();
  virtual Osi_Expression* add(OsiSolver *solver, bool top_level);
};


class Osi_Maximise : public Osi_Expression
{
private:
  Osi_Expression* _exp;

public:
  Osi_Maximise(Osi_Expression *var);
  virtual ~Osi_Maximise();
  virtual Osi_Expression* add(OsiSolver *solver, bool top_level);
};


/**
   The solver itself
*/
class OsiSolver
{

public:

  int first_decision_level;
  int saved_level;

  OsiSolver();
  virtual ~OsiSolver();

  // add an expression, in the case of a tree of expressions,
  // each node of the tree is added separately, depth first.
  void add(Osi_Expression* arg);

  // used to initialise search on a given subset of variables
  void initialise(OsiExpArray& arg);
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

  int next(Osi_Expression* x, int v);

  int get_level();

  bool propagate();
  void save();
  void post(const char* op, Osi_Expression* x, int v);
  bool undo(const int nlevel);
  void deduce(const char* op, Osi_Expression* x, int v);
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
  void forceFiniteDomain(OsiExpArray& vars);
  void addNogood(OsiExpArray& vars, 
		 OsiIntArray& vals);
  void guide(OsiExpArray& vars, 
	     OsiIntArray& vals,
	     OsiDoubleArray& probs);
  void backtrackTo(const int level);
  void upOneLevel();
  void presolve();
  void assign(Osi_Expression *X, const int v);
  void increase_init_level(const int i);
  void decrease_init_level(const int i);
  void reset(bool full);
  void setLowerBounds(OsiExpArray& vars, 
		      OsiIntArray& vals);
  void setUpperBounds(OsiExpArray& vars, 
		      OsiIntArray& vals);
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
