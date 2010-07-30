
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
const int SELF      =  0;
const int DIRECT    =  1;
const int ORDER     =  2;



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
  void set_item(const int i, T item) { _array[i]=item; }
};

typedef MiniSatArray< int > MiniSatIntArray;
typedef MiniSatArray< double > MiniSatDoubleArray;







class MiniSat_Expression;
class MiniSatSolver;
class AbstractDomain {
public:
  MiniSat_Expression *owner;

  AbstractDomain(MiniSat_Expression *o) {owner = o;}

  virtual int getval(int idx) const = 0;
  virtual int getmin() const = 0;
  virtual int getmax() const = 0;
  virtual int getsize() const = 0;

  virtual void encode(MiniSatSolver *solver) {}

  virtual Lit less_or_equal(const int value, const int index) const = 0;
  virtual Lit equal(const int value, const int index) const = 0;

  virtual void print_lit(Lit p, const int type) const {
    int atom = var(p);
    std::cout << "(";
    if(!sign(p)) std::cout << "~";
    std::cout << atom << ")";
  }
  
  virtual std::ostream& display(std::ostream& o) const = 0;
};

class OffsetDomain : public AbstractDomain {
public:
  AbstractDomain *_dom_ptr;
  int offset;

  OffsetDomain(MiniSat_Expression *os, AbstractDomain *d, const int of);

  virtual int getval(int idx) const;
  virtual int getmin() const;
  virtual int getmax() const;
  virtual int getsize() const {return _dom_ptr->getsize();}

  virtual Lit less_or_equal(const int value, const int index) const;
  virtual Lit equal(const int value, const int index) const;

  virtual std::ostream& display(std::ostream& o) const;
};

class FactorDomain : public AbstractDomain {
public:
  AbstractDomain *_dom_ptr;
  int factor;

  FactorDomain(MiniSat_Expression *o, AbstractDomain *d, const int f);

  virtual int getval(int idx) const;
  virtual int getmin() const;
  virtual int getmax() const;
  virtual int getsize() const {return _dom_ptr->getsize();}

  virtual Lit less_or_equal(const int value, const int index) const;
  virtual Lit equal(const int value, const int index) const;

  virtual std::ostream& display(std::ostream& o) const;
};

class EqDomain : public AbstractDomain {
public:
  AbstractDomain *_dom_ptr;
  int value;
  int spin;

  EqDomain(MiniSat_Expression *o, AbstractDomain *d, const int v, const int s);

  virtual int getval(int idx) const {assert(idx >=0 && idx <= 1); return idx;}
  virtual int getmin() const {return 0;}
  virtual int getmax() const {return 1;}
  virtual int getsize() const {return 2;}

  virtual Lit less_or_equal(const int value, const int index) const;
  virtual Lit equal(const int value, const int index) const;

  virtual std::ostream& display(std::ostream& o) const;
};

class LeqDomain : public AbstractDomain {
public:
  AbstractDomain *_dom_ptr;
  int bound;
  int spin;

  LeqDomain(MiniSat_Expression *o, AbstractDomain *d, const int b, const int s);

  virtual int getval(int idx) const {assert(idx >=0 && idx <= 1); return idx;}
  virtual int getmin() const {return 0;}
  virtual int getmax() const {return 1;}
  virtual int getsize() const {return 2;}

  virtual Lit less_or_equal(const int value, const int index) const;
  virtual Lit equal(const int value, const int index) const;

  virtual std::ostream& display(std::ostream& o) const;
};

class ConstantDomain : public AbstractDomain {
public:
  int value;

  ConstantDomain(MiniSat_Expression *o, const int v) : AbstractDomain(o) {value = v;}

  virtual int getval(int idx) const { assert(idx == 0); return value; }
  virtual int getmin() const {return value;}
  virtual int getmax() const {return value;}
  virtual int getsize() const {return 1;}

  virtual Lit less_or_equal(const int v, const int index) const; 
  virtual Lit equal(const int v, const int index) const; 

  virtual std::ostream& display(std::ostream& o) const;
};

class DomainEncoding : public AbstractDomain{

private:

  // domain's lower bound
  int _lower;

  // domain's upper bound
  int _upper;

  // domain's size
  int _size;

  /** 
      if _values is left NULL, then _direct_encoding keeps the 
      index of the first variable used for the direct encoding 
      (similar for _order_encoding).
      One can get the variable for a value v by adding and
      substracting _lower.
      A particular case is when _size = 2.
      In that case, the variable is boolean and only one
      SAT variable is used.

      If _value is not NULL, then _direct_encoding
      is the index of the first variable in 
   */
  // domain's values (possibly null)
  int *_values;

  // mapping from values to direc encoding (possibly null)
  int _direct_encoding;

  // mapping from values to interval encoding (possibly null)
  int _order_encoding;

  // return the index of 'value', or the next index if value is not in the domain
  int get_index_n(const int v) const {
    int lb = 0, ub = _size-1, x;
    while(lb < ub) {
      x = (lb+ub)/2;
      if(_values[x] == v) return x;
      if(_values[x] < v) lb = x+1;
      else ub = x;
    }
    return lb;
  }

  // return the index of 'value', or the prev index if value is not in the domain
  int get_index_p(const int v) const {
    int lb = 0, ub = _size-1, x;
    while(lb < ub) {
      x = (lb+ub)/2+((lb+ub)%2);
      if(_values[x] == v) return x;
      if(_values[x] > v) ub = x-1;
      else lb = x;
    }
    return lb;
  }


public:

  DomainEncoding(MiniSat_Expression *o);
  DomainEncoding(MiniSat_Expression *o, const int nval);
  DomainEncoding(MiniSat_Expression *o, const int lb, const int ub);
  DomainEncoding(MiniSat_Expression *o, MiniSatIntArray& vals);

  virtual ~DomainEncoding();


  virtual void print_lit(Lit p, const int type) const;


  virtual int getval(int idx) const {
    if(_values) return _values[idx%_size];
    else return _lower+idx;
  }

  virtual int getmin() const {return _lower;}
  virtual int getmax() const {return _upper;}
  virtual int getsize() const {return _size;}

  void encode(MiniSatSolver *solver);

  virtual Lit less_or_equal(const int value, const int index) const;
  virtual Lit equal(const int value, const int index) const;

  virtual std::ostream& display(std::ostream& o) const;

};





/**
   Expression (Used to encode variables & constraints)
*/
class MiniSatSolver;
class MiniSat_Expression
{
    
public:
  
  MiniSatSolver *_solver;

  // unique identifier
  int _ident;
  int nbj_ident;

  AbstractDomain *domain;

  bool has_been_added() const;

  int getval(int idx) const { return domain->getval(idx); }
  int getmin() const { return domain->getmin(); }
  int getmax() const { return domain->getmax(); }
  int getsize() const { return domain->getsize(); }

  int get_value() const ;
  int get_min() const ;
  int get_max() const ;
  int get_size() const ;
  int next(const int value) const ;
  bool contain(const int value) const ;


  Lit greater_than(const int value, const int index=-1) const { return ~(domain->less_or_equal(value,index)); }
  Lit less_or_equal(const int value, const int index=-1) const { return domain->less_or_equal(value,index); }
  Lit equal(const int value, const int index=-1) const { return domain->equal(value,index); }

  void initialise();

  MiniSat_Expression();
  MiniSat_Expression(const int nval);
  MiniSat_Expression(const int lb, const int ub);
  MiniSat_Expression(MiniSatIntArray& vals);
  virtual ~MiniSat_Expression();

  virtual MiniSat_Expression* add(MiniSatSolver *solver, bool top_level);
};

class MiniSat_IntVar : public MiniSat_Expression
{

public:
  
  MiniSat_IntVar() : MiniSat_Expression() {}
  MiniSat_IntVar(const int lb, const int ub, const int ident) : MiniSat_Expression(lb, ub) {nbj_ident = ident;}
  MiniSat_IntVar(MiniSatIntArray& vals, const int ident) : MiniSat_Expression(vals) {nbj_ident = ident;}

};


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
  int arity() { return 1+(_vars[1]==NULL); }
  MiniSat_binop(MiniSat_Expression *var1, MiniSat_Expression *var2);
  MiniSat_binop(MiniSat_Expression *var1, int rhs);
  virtual ~MiniSat_binop();

  virtual MiniSat_Expression* add(MiniSatSolver *solver, bool top_level) = 0;
};


class MiniSat_add : public MiniSat_binop
{
public:
  MiniSat_add( MiniSat_Expression* arg1, MiniSat_Expression* arg2 );
  MiniSat_add( MiniSat_Expression* arg1, const int arg2 );
  virtual ~MiniSat_add();
  virtual MiniSat_Expression* add(MiniSatSolver *solver, bool top_level);
};

class MiniSat_mul : public MiniSat_binop
{
public:
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


bool processClause(std::vector<Lit>& cl_in, std::vector<Lit>& cl_out);


class MiniSat_eq: public MiniSat_binop
{
public:
  MiniSat_eq(MiniSat_Expression *var1, MiniSat_Expression *var2);
  MiniSat_eq(MiniSat_Expression *var1, int rhs);
  virtual ~MiniSat_eq();
  virtual MiniSat_Expression* add(MiniSatSolver *solver, bool top_level);
};

class MiniSat_ne: public MiniSat_binop
{
public:
  MiniSat_ne(MiniSat_Expression *var1, MiniSat_Expression *var2);
  MiniSat_ne(MiniSat_Expression *var1, int rhs);
  virtual ~MiniSat_ne();
  virtual MiniSat_Expression* add(MiniSatSolver *solver, bool top_level);
};

class MiniSat_le: public MiniSat_binop
{
public:
  MiniSat_le(MiniSat_Expression *var1, MiniSat_Expression *var2);
  MiniSat_le(MiniSat_Expression *var1, int rhs);
  virtual ~MiniSat_le();
  virtual MiniSat_Expression* add(MiniSatSolver *solver, bool top_level);
};

class MiniSat_ge: public MiniSat_binop
{
public:
  MiniSat_ge(MiniSat_Expression *var1, MiniSat_Expression *var2);
  MiniSat_ge(MiniSat_Expression *var1, int rhs);
  virtual ~MiniSat_ge();
  virtual MiniSat_Expression* add(MiniSatSolver *solver, bool top_level);
};

class MiniSat_lt: public MiniSat_binop
{
public:
  MiniSat_lt(MiniSat_Expression *var1, MiniSat_Expression *var2);
  MiniSat_lt(MiniSat_Expression *var1, int rhs);
  virtual ~MiniSat_lt();
  virtual MiniSat_Expression* add(MiniSatSolver *solver, bool top_level);
};

class MiniSat_gt: public MiniSat_binop
{
public:
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
  ////////////// MiniSat Specific ////////////////
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
  Lit last_decision;
  int saved_level;
  ////////////// MiniSat Specific ////////////////

public:

  MiniSat_Expression *minimise_obj;
  MiniSat_Expression *maximise_obj;

  // repository for all expressions
  std::vector< MiniSat_Expression* > _expressions;
  std::vector< MiniSat_Expression* > _variables;
  // link each atom to its domain
  std::vector< DomainEncoding* > _atom_to_domain;
  std::vector< int > _atom_to_type;

  std::vector< std::vector<Lit> > clause_base;
  unsigned int current;

  int *cp_model;


  MiniSatSolver();
  virtual ~MiniSatSolver();
  int nbClauses() {return nClauses();}

  // add an expression, in the case of a tree of expressions,
  // each node of the tree is added separately, depth first.
  void add(MiniSat_Expression* arg);
  int declare(MiniSat_Expression* arg, const bool type);
  int create_atom(DomainEncoding* dom, const int type);

  void addClause(std::vector<Lit>& cl);
  void validate();
  void displayClause(std::vector<Lit>& cl); 
  void displayLiteral(Lit p);


  lbool truth_value(Lit x);


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

  void reset(bool full);
  bool propagate();
  void save();
  void post(const char* op, MiniSat_Expression* x, int v);
  bool undo(const int nlevel);
  void deduce();
  bool branch_right();

  void store_solution();
  
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

