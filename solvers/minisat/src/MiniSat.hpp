
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



/**
   Expressions represent predicate trees.
   
   First, trees are created, bottom-up
     - on creation, the parent indicate the enoding type to the children
     
   Then, each tree is added, the parent adding its children first (dfs)
     - on addition, the domain is first encoded (using the necessary type(s))
     - then, if needed, the relation is encoded using the children's methods:
        - less_than(v)
        - bigger_than(v)
        - equal_to(v)
        - not_equal_to(v)
 */



const int UNSAT     =  0;
const int SAT       =  1;
const int UNKNOWN   =  2;
const int LUBY      =  0;
const int GEOMETRIC =  1;
const int GLUBY     =  2;  



/**
   Array of expressions (for nary constraints)
*/
template<class T>
class MiniSatArray
{
private:
  std::vector< T > _array;

public:
  MiniSatArray() {}
  virtual ~MiniSatArray() {}
  int size() const { return _array.size(); }
  void add(T arg) { _array.push_back(arg); }
  T get_item(const int i) const { return _array[i]; }
  //T& get_item(const int i) { return _array[i]; }
  void set_item(const int i, T item) { _array[i]=item; }
};

typedef MiniSatArray< int > MiniSatIntArray;
typedef MiniSatArray< double > MiniSatDoubleArray;

/**
   Expression (Used to encode variables & constraints)
*/
class MiniSatSolver;
class MiniSat_Expression
{
    
public:

  static const int NO    = 0;
  static const int DIRECT= 1;
  static const int ORDER = 2;

  MiniSatSolver *_solver;

  // unique identifier
  int _ident;
  int nbj_ident;

  // domain's lower bound
  int _lower;

  // domain's upper bound
  int _upper;

  // domain's size
  int _size;

  // domain's values (possibly null)
  int *_values;

  // mapping from values to direc encoding (possibly null)
  Var *_direct_encoding;

  // mapping from values to interval encoding (possibly null)
  Var *_order_encoding;

  // store the required encoding type
  int _encoding_type;


  //virtual int get_lower() const;
  //virtual int get_upper() const;
  virtual int getNextEq(const int v) const;
  virtual int getNext(const int v) const;
  virtual int getPrevEq(const int v) const;
  virtual int getPrev(const int v) const;
  virtual int getval(const int i) const;
  virtual int getind(const int v) const;

  bool has_been_added() const;
  virtual void setDirectEncoding();
  virtual void setOrderEncoding();

  virtual Lit less_or_equal(const int value) const;
  virtual Lit greater_than(const int value) const;
  virtual Lit equal(const int value) const;
  //virtual Lit not_equal(const int value);

  // returns the atom used to represent this value in the direct encoding
  //virtual Var getDirectVar(const int value);
  // returns the atom used to represent this value in the order encoding 
  //virtual Var getInterVar(const int value); 
  virtual void encode(MiniSatSolver* solver);
  virtual void channel(MiniSatSolver* solver);

  void initialise();
  MiniSat_Expression();
  MiniSat_Expression(const int nval);
  MiniSat_Expression(const int lb, const int ub);
  MiniSat_Expression(MiniSatIntArray& vals);
  virtual ~MiniSat_Expression();

  int next(int v);
  int get_value() const;
  virtual int get_min() const;
  virtual int get_max() const;
  virtual int get_size() const;
  bool contain(const int v) const;
  virtual MiniSat_Expression* add(MiniSatSolver *solver, bool top_level);
};

class MiniSat_IntVar : public MiniSat_Expression
{

public:
  
  MiniSat_IntVar() : MiniSat_Expression() {}
  //MiniSat_IntVar(const int nval) : MiniSat_Expression(nval) {}
  MiniSat_IntVar(const int lb, const int ub, const int ident) : MiniSat_Expression(lb, ub) {nbj_ident = ident;}
  MiniSat_IntVar(MiniSatIntArray& vals, const int ident) : MiniSat_Expression(vals) {nbj_ident = ident;}

};

//typedef MiniSat_Expression MiniSat_IntVar;

typedef MiniSatArray< MiniSat_Expression* > MiniSatExpArray;

class MiniSat_AllDiff : public MiniSat_Expression
{
private:
  MiniSatExpArray _vars;
  std::vector< MiniSat_Expression* > _clique;

public:
  MiniSat_AllDiff(MiniSatExpArray& vars);
  MiniSat_AllDiff(MiniSat_Expression* arg1, MiniSat_Expression* arg2);
  void addVar( MiniSat_Expression* v );
  virtual ~MiniSat_AllDiff();

  virtual MiniSat_Expression* add(MiniSatSolver *solver, bool top_level);
};


class MiniSat_Sum: public MiniSat_Expression
{
private:
  MiniSat_Expression *_self;
  int _offset;
  MiniSatExpArray _vars;
  MiniSatIntArray _weights;
  std::vector< MiniSat_Expression* > _subsum;
    
public:
  MiniSat_Sum(MiniSatExpArray& vars, MiniSatIntArray& weights, const int offset=0);
  MiniSat_Sum( MiniSat_Expression* arg1, 
	       MiniSat_Expression* arg2, 
	       MiniSatIntArray& weights,
	       const int offset );
  MiniSat_Sum( MiniSat_Expression* arg, 
	       MiniSatIntArray& weights,
	       const int offset );
  MiniSat_Sum();
  void initialise();


  virtual Lit less_or_equal(const int value) const;
  virtual Lit greater_than(const int value) const;
  virtual Lit equal(const int value) const;
  //virtual int get_min() const;
  //virtual Lit not_equal(const int value);
//   virtual Var getDirectVar(const int value);
//   virtual Var getInterVar(const int value);

  virtual int get_min() const;
  virtual int get_max() const;
  virtual int get_size() const;


  void addVar( MiniSat_Expression* v );
  void addWeight( const int w );
  void set_rhs( const int k );

  virtual ~MiniSat_Sum();
  virtual MiniSat_Expression* add(MiniSatSolver *solver, bool top_level);
};


class MiniSat_binop: public MiniSat_Expression
{
protected:
  MiniSat_Expression *_vars[2];
  int _rhs;
  
public:
  virtual void setDirectEncoding();
  virtual void setOrderEncoding();
  int arity() { return 1+(_vars[1]==NULL); }
  MiniSat_binop(MiniSat_Expression *var1, MiniSat_Expression *var2);
  MiniSat_binop(MiniSat_Expression *var1, int rhs);
  virtual ~MiniSat_binop();

  virtual MiniSat_Expression* add(MiniSatSolver *solver, bool top_level) = 0;
};


class MiniSat_add : public MiniSat_binop
{
public:
  //virtual int get_lower() const;
  //virtual int get_upper() const;
  virtual int getval(const int i) const;
  virtual Lit less_or_equal(const int value) const;
  virtual Lit greater_than(const int value) const;
  virtual Lit equal(const int value) const;
  virtual int get_min() const;
  //virtual Lit not_equal(const int value);
  MiniSat_add( MiniSat_Expression* arg1, MiniSat_Expression* arg2 );
  MiniSat_add( MiniSat_Expression* arg1, const int arg2 );
  virtual ~MiniSat_add();
  virtual MiniSat_Expression* add(MiniSatSolver *solver, bool top_level);
};
// class MiniSat_sub : public MiniSat_binop
// {
// public:
//   virtual Var getDirectVar(const int value);
//   virtual Var getInterVar(const int value); 
//   MiniSat_sub( MiniSat_Expression* arg1, MiniSat_Expression* arg2 );
//   MiniSat_sub( MiniSat_Expression* arg1, const int arg2 );
//   virtual ~MiniSat_sub();
//   virtual MiniSat_Expression* add(MiniSatSolver *solver, bool top_level);
// };
class MiniSat_mul : public MiniSat_binop
{
public:
  //virtual int get_lower() const;
  // virtual int get_upper() const;
  virtual int getval(const int i) const;
  virtual Lit less_or_equal(const int value) const;
  virtual Lit greater_than(const int value) const;
  virtual Lit equal(const int value) const;
  virtual int get_min() const;
  //virtual Lit not_equal(const int value);
  MiniSat_mul( MiniSat_Expression* arg1, MiniSat_Expression* arg2 );
  MiniSat_mul( MiniSat_Expression* arg1, const int arg2 );
  virtual ~MiniSat_mul();
  virtual MiniSat_Expression* add(MiniSatSolver *solver, bool top_level);
};

class MiniSat_or: public MiniSat_binop
{
public:
  MiniSat_or(MiniSat_Expression *var1, MiniSat_Expression *var2);
  MiniSat_or(MiniSat_Expression *var1, int rhs);
  virtual ~MiniSat_or();
  virtual MiniSat_Expression* add(MiniSatSolver *solver, bool top_level);
};

class MiniSat_and: public MiniSat_binop
{
public:
  MiniSat_and(MiniSat_Expression *var1, MiniSat_Expression *var2);
  MiniSat_and(MiniSat_Expression *var1, int rhs);
  virtual ~MiniSat_and();
  virtual MiniSat_Expression* add(MiniSatSolver *solver, bool top_level);
};

class MiniSat_eq: public MiniSat_binop
{
public:
  virtual Lit less_or_equal(const int value) const;
  virtual Lit greater_than(const int value) const;
  virtual Lit equal(const int value) const;
  //virtual Lit not_equal(const int value) const;
  MiniSat_eq(MiniSat_Expression *var1, MiniSat_Expression *var2);
  MiniSat_eq(MiniSat_Expression *var1, int rhs);
  virtual ~MiniSat_eq();
  virtual MiniSat_Expression* add(MiniSatSolver *solver, bool top_level);
};

class MiniSat_ne: public MiniSat_binop
{
public:
  virtual Lit less_or_equal(const int value) const;
  virtual Lit greater_than(const int value) const;
  virtual Lit equal(const int value) const;
  //virtual Lit not_equal(const int value) const;
  MiniSat_ne(MiniSat_Expression *var1, MiniSat_Expression *var2);
  MiniSat_ne(MiniSat_Expression *var1, int rhs);
  virtual ~MiniSat_ne();
  virtual MiniSat_Expression* add(MiniSatSolver *solver, bool top_level);
};

class MiniSat_le: public MiniSat_binop
{
public:
  virtual Lit less_or_equal(const int value) const;
  virtual Lit greater_than(const int value) const;
  virtual Lit equal(const int value) const;
  //virtual Lit not_equal(const int value) const;
  MiniSat_le(MiniSat_Expression *var1, MiniSat_Expression *var2);
  MiniSat_le(MiniSat_Expression *var1, int rhs);
  virtual ~MiniSat_le();
  virtual MiniSat_Expression* add(MiniSatSolver *solver, bool top_level);
};

class MiniSat_ge: public MiniSat_binop
{
public:
  virtual Lit less_or_equal(const int value) const;
  virtual Lit greater_than(const int value) const;
  virtual Lit equal(const int value) const;
  //virtual Lit not_equal(const int value) const;
  MiniSat_ge(MiniSat_Expression *var1, MiniSat_Expression *var2);
  MiniSat_ge(MiniSat_Expression *var1, int rhs);
  virtual ~MiniSat_ge();
  virtual MiniSat_Expression* add(MiniSatSolver *solver, bool top_level);
};

class MiniSat_lt: public MiniSat_binop
{
public:
  virtual Lit less_or_equal(const int value) const;
  virtual Lit greater_than(const int value) const;
  virtual Lit equal(const int value) const;
  //virtual Lit not_equal(const int value) const;
  MiniSat_lt(MiniSat_Expression *var1, MiniSat_Expression *var2);
  MiniSat_lt(MiniSat_Expression *var1, int rhs);
  virtual ~MiniSat_lt();
  virtual MiniSat_Expression* add(MiniSatSolver *solver, bool top_level);
};

class MiniSat_gt: public MiniSat_binop
{
public:
  virtual Lit less_or_equal(const int value) const;
  virtual Lit greater_than(const int value) const;
  virtual Lit equal(const int value) const;
  //virtual Lit not_equal(const int value) const;
  MiniSat_gt(MiniSat_Expression *var1, MiniSat_Expression *var2);
  MiniSat_gt(MiniSat_Expression *var1, int rhs);
  virtual ~MiniSat_gt();
  virtual MiniSat_Expression* add(MiniSatSolver *solver, bool top_level);
};

class MiniSat_Minimise : public MiniSat_Expression
{
protected:
  MiniSat_Expression *_obj;
public:
  MiniSat_Minimise(MiniSat_Expression *var1);
  virtual ~MiniSat_Minimise();
  virtual MiniSat_Expression* add(MiniSatSolver *solver, bool top_level);
};

class MiniSat_Maximise : public MiniSat_Expression
{
protected:
  MiniSat_Expression *_obj;
public:
  MiniSat_Maximise(MiniSat_Expression *var1);
  virtual ~MiniSat_Maximise();
  virtual MiniSat_Expression* add(MiniSatSolver *solver, bool top_level);
};



/**
   The solver itself
*/
class MiniSatSolver : public SimpSolver
{

private:
  //SimpSolver minisolver;
  double STARTTIME;
  lbool result;
  int nbSolutions;

  // search stuff
  Clause*     conflict_clause;
  int         backtrack_level;
  int         conflictC;
  vec<Lit>    forced_decisions;
  vec<Lit>    learnt_clause;

  int first_decision_level;

public:

  MiniSat_Expression *minimise_obj;
  MiniSat_Expression *maximise_obj;
  std::vector< MiniSat_Expression* > _variables;
  std::vector< MiniSat_Expression* > _lit_to_var;
  std::vector< int >                 _lit_to_val;
  int var_counter;
  int *cp_model;
  Lit last_decision;
  
  int saved_level;

  MiniSatSolver();
  virtual ~MiniSatSolver();

  int nbClauses() {return nClauses();}
  void get_clause(int i) { learnt_clause.clear(); for(int j=0; j<clauses[i]->size(); ++j) learnt_clause.push((*(clauses[i]))[j]); }
  int get_nogood_size() {return learnt_clause.size();}
  int get_nogood_var(int i)  {return _lit_to_var[var(learnt_clause[i])]->nbj_ident;}
  int get_nogood_val(int i)  {return ((_lit_to_val[var(learnt_clause[i])]) >> 1);}
  int get_nogood_type(int i) {return (_lit_to_val[var(learnt_clause[i])] % 2);}
  int get_nogood_sign(int i) {return sign(learnt_clause[i]);}

  //int get_value(const int id);
  //void store_solution();
  int declare(MiniSat_Expression *exp);

//   int nVars();
  Var newVar(MiniSat_Expression* var, int val);
  //void addClause(vec<Lit>& cl);
  lbool truth_value(Lit x);

  void reset(bool full);
  bool propagate();
  //void branch_on(const char* op, Mistral_Expression* x, const int v);
  void save();
  void post(const char* op, MiniSat_Expression* x, int v);
  bool undo(const int nlevel);
  //bool backjump();
  void deduce();
  bool branch_right();
  //void analyze_conflict();
  //void learn(vec<Lit>& learnt_cl);

  void store_solution();

  // add an expression, in the case of a tree of expressions,
  // each node of the tree is added separately, depth first.
  void add(MiniSat_Expression* arg);

  // used to initialise search on a given subset of variables
  void initialise(MiniSatExpArray& arg);
  // initialise the solver before solving (no more calls to add after this)
  void initialise();

  // solving methods
  int solve();
  int solveAndRestart(const int policy = GEOMETRIC, 
		      const unsigned int base = 32, 
		      const double factor = 1.3333333,
		      const double decay = 0.0);
  int startNewSearch();
  int getNextSolution();
  int sacPreprocess(const int type);
  
  // parameter tuning methods
  void guide(MiniSatExpArray& vars, 
	     MiniSatIntArray& vals,
	     MiniSatDoubleArray& probs) {}
  void setHeuristic(const char* var_heuristic, const char* val_heuristic, const int rand);
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
  int getChecks();
  int getPropags();
  double getTime();

  int solveDimacs(const char*);
};


#endif

