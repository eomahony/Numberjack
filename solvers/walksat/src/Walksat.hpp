
#ifndef WALKSAT_H
#define WALKSAT_H


#include <vector>
#include <iostream>

#include "cpp_walksat.hpp"
#include "Literals.h"

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
class WalksatArray
{
private:
  std::vector< T > _array;

public:
  WalksatArray() {}
  virtual ~WalksatArray() {}
  int size() const { return _array.size(); }
  void add(T arg) { _array.push_back(arg); }
  T get_item(const int i) const { return _array[i]; }
  void set_item(const int i, T item) { _array[i]=item; }
};

typedef WalksatArray< int > WalksatIntArray;
typedef WalksatArray< double > WalksatDoubleArray;



class Walksat_Expression;
class WalksatSolver;
class AbstractDomain {
public:
  Walksat_Expression *owner;

  AbstractDomain(Walksat_Expression *o) {owner = o;}

  virtual int getval(int idx) const = 0;
  virtual int getmin() const = 0;
  virtual int getmax() const = 0;
  virtual int getsize() const = 0;

  virtual void encode(WalksatSolver *solver) {}

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

  OffsetDomain(Walksat_Expression *o, AbstractDomain *d, const int o);

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

  FactorDomain(Walksat_Expression *o, AbstractDomain *d, const int f);

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

  EqDomain(Walksat_Expression *o, AbstractDomain *d, const int v, const int s);

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

  LeqDomain(Walksat_Expression *o, AbstractDomain *d, const int b, const int s);

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

  ConstantDomain(Walksat_Expression *o, const int v) : AbstractDomain(o) {value = v;}

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

  DomainEncoding(Walksat_Expression *o);
  DomainEncoding(Walksat_Expression *o, const int nval);
  DomainEncoding(Walksat_Expression *o, const int lb, const int ub);
  DomainEncoding(Walksat_Expression *o, WalksatIntArray& vals);

  virtual ~DomainEncoding();


  virtual void print_lit(Lit p, const int type) const;


  virtual int getval(int idx) const {
    if(_values) return _values[idx%_size];
    else return _lower+idx;
  }

  virtual int getmin() const {return _lower;}
  virtual int getmax() const {return _upper;}
  virtual int getsize() const {return _size;}

  void encode(WalksatSolver *solver);

  virtual Lit less_or_equal(const int value, const int index) const;
  virtual Lit equal(const int value, const int index) const;

  virtual std::ostream& display(std::ostream& o) const;

};





/**
   Expression (Used to encode variables & constraints)
*/
class WalksatSolver;
class Walksat_Expression
{
    
public:
  
  WalksatSolver *_solver;

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

  Walksat_Expression();
  Walksat_Expression(const int nval);
  Walksat_Expression(const int lb, const int ub);
  Walksat_Expression(WalksatIntArray& vals);
  virtual ~Walksat_Expression();

  virtual Walksat_Expression* add(WalksatSolver *solver, bool top_level);
};

class Walksat_IntVar : public Walksat_Expression
{

public:
  
  Walksat_IntVar() : Walksat_Expression() {}
  Walksat_IntVar(const int lb, const int ub, const int ident) : Walksat_Expression(lb, ub) {nbj_ident = ident;}
  Walksat_IntVar(WalksatIntArray& vals, const int ident) : Walksat_Expression(vals) {nbj_ident = ident;}

};


typedef WalksatArray< Walksat_Expression* > WalksatExpArray;

class Walksat_AllDiff : public Walksat_Expression
{
private:
  WalksatExpArray _vars;
  std::vector< Walksat_Expression* > _clique;

public:
  Walksat_AllDiff(WalksatExpArray& vars);
  Walksat_AllDiff(Walksat_Expression* arg1, Walksat_Expression* arg2);
  void addVar( Walksat_Expression* v );
  virtual ~Walksat_AllDiff();

  virtual Walksat_Expression* add(WalksatSolver *solver, bool top_level);
};


class Walksat_Sum: public Walksat_Expression
{
private:
  Walksat_Expression *_self;
  int _offset;
  WalksatExpArray _vars;
  WalksatIntArray _weights;
  std::vector< Walksat_Expression* > _subsum;
    
public:
  Walksat_Sum(WalksatExpArray& vars, WalksatIntArray& weights, const int offset=0);
  Walksat_Sum( Walksat_Expression* arg1, 
	       Walksat_Expression* arg2, 
	       WalksatIntArray& weights,
	       const int offset );
  Walksat_Sum( Walksat_Expression* arg, 
	       WalksatIntArray& weights,
	       const int offset );
  Walksat_Sum();
  void initialise();


  void addVar( Walksat_Expression* v );
  void addWeight( const int w );
  void set_rhs( const int k );

  virtual ~Walksat_Sum();
  virtual Walksat_Expression* add(WalksatSolver *solver, bool top_level);
};


class Walksat_binop: public Walksat_Expression
{
protected:
  Walksat_Expression *_vars[2];
  int _rhs;
  
public:
  int arity() { return 1+(_vars[1]==NULL); }
  Walksat_binop(Walksat_Expression *var1, Walksat_Expression *var2);
  Walksat_binop(Walksat_Expression *var1, int rhs);
  virtual ~Walksat_binop();

  virtual Walksat_Expression* add(WalksatSolver *solver, bool top_level) = 0;
};


class Walksat_add : public Walksat_binop
{
public:
  Walksat_add( Walksat_Expression* arg1, Walksat_Expression* arg2 );
  Walksat_add( Walksat_Expression* arg1, const int arg2 );
  virtual ~Walksat_add();
  virtual Walksat_Expression* add(WalksatSolver *solver, bool top_level);
};

class Walksat_mul : public Walksat_binop
{
public:
  Walksat_mul( Walksat_Expression* arg1, Walksat_Expression* arg2 );
  Walksat_mul( Walksat_Expression* arg1, const int arg2 );
  virtual ~Walksat_mul();
  virtual Walksat_Expression* add(WalksatSolver *solver, bool top_level);
};

class Walksat_or: public Walksat_binop
{
public:
  Walksat_or(Walksat_Expression *var1, Walksat_Expression *var2);
  Walksat_or(Walksat_Expression *var1, int rhs);
  virtual ~Walksat_or();
  virtual Walksat_Expression* add(WalksatSolver *solver, bool top_level);
};

class Walksat_and: public Walksat_binop
{
public:
  Walksat_and(Walksat_Expression *var1, Walksat_Expression *var2);
  Walksat_and(Walksat_Expression *var1, int rhs);
  virtual ~Walksat_and();
  virtual Walksat_Expression* add(WalksatSolver *solver, bool top_level);
};


bool processClause(std::vector<Lit>& cl_in, std::vector<Lit>& cl_out);


class Walksat_eq: public Walksat_binop
{
public:
  Walksat_eq(Walksat_Expression *var1, Walksat_Expression *var2);
  Walksat_eq(Walksat_Expression *var1, int rhs);
  virtual ~Walksat_eq();
  virtual Walksat_Expression* add(WalksatSolver *solver, bool top_level);
};

class Walksat_ne: public Walksat_binop
{
public:
  Walksat_ne(Walksat_Expression *var1, Walksat_Expression *var2);
  Walksat_ne(Walksat_Expression *var1, int rhs);
  virtual ~Walksat_ne();
  virtual Walksat_Expression* add(WalksatSolver *solver, bool top_level);
};

class Walksat_le: public Walksat_binop
{
public:
  Walksat_le(Walksat_Expression *var1, Walksat_Expression *var2);
  Walksat_le(Walksat_Expression *var1, int rhs);
  virtual ~Walksat_le();
  virtual Walksat_Expression* add(WalksatSolver *solver, bool top_level);
};

class Walksat_ge: public Walksat_binop
{
public:
  Walksat_ge(Walksat_Expression *var1, Walksat_Expression *var2);
  Walksat_ge(Walksat_Expression *var1, int rhs);
  virtual ~Walksat_ge();
  virtual Walksat_Expression* add(WalksatSolver *solver, bool top_level);
};

class Walksat_lt: public Walksat_binop
{
public:
  Walksat_lt(Walksat_Expression *var1, Walksat_Expression *var2);
  Walksat_lt(Walksat_Expression *var1, int rhs);
  virtual ~Walksat_lt();
  virtual Walksat_Expression* add(WalksatSolver *solver, bool top_level);
};

class Walksat_gt: public Walksat_binop
{
public:
  Walksat_gt(Walksat_Expression *var1, Walksat_Expression *var2);
  Walksat_gt(Walksat_Expression *var1, int rhs);
  virtual ~Walksat_gt();
  virtual Walksat_Expression* add(WalksatSolver *solver, bool top_level);
};

class Walksat_Minimise : public Walksat_Expression
{
protected:
  Walksat_Expression *_obj;
public:
  Walksat_Minimise(Walksat_Expression *var1);
  virtual ~Walksat_Minimise();
  virtual Walksat_Expression* add(WalksatSolver *solver, bool top_level);
};

class Walksat_Maximise : public Walksat_Expression
{
protected:
  Walksat_Expression *_obj;
public:
  Walksat_Maximise(Walksat_Expression *var1);
  virtual ~Walksat_Maximise();
  virtual Walksat_Expression* add(WalksatSolver *solver, bool top_level);
};


/**
   The solver itself
*/
class WalksatSolver 
{

private:
  ////////////// Walksat Specific ////////////////
  WalksatAlgorithm wsat;

  double STARTTIME;
  int nbSolutions;
  ////////////// Walksat Specific ////////////////

public:

  Walksat_Expression *minimise_obj;
  Walksat_Expression *maximise_obj;

  // repository for all expressions
  std::vector< Walksat_Expression* > _expressions;
  std::vector< Walksat_Expression* > _variables;
  // link each atom to its domain
  std::vector< DomainEncoding* > _atom_to_domain;
  std::vector< int > _atom_to_type;

  std::vector< std::vector<Lit> > clause_base;
  unsigned int current;

  int *cp_model;


  WalksatSolver();
  virtual ~WalksatSolver();
  int nbClauses() {return clause_base.size();}

  // add an expression, in the case of a tree of expressions,
  // each node of the tree is added separately, depth first.
  void add(Walksat_Expression* arg);
  int declare(Walksat_Expression* arg, const bool type);
  int create_atom(DomainEncoding* dom, const int type);

  void addClause(std::vector<Lit>& cl);
  void validate();
  void displayClause(std::vector<Lit>& cl); 
  void displayLiteral(Lit p);


  lbool truth_value(Lit x);


  // used to initialise search on a given subset of variables
  void initialise(WalksatExpArray& arg);
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
  void post(const char* op, Walksat_Expression* x, int v);
  bool undo(const int nlevel);
  void deduce();
  bool branch_right();

  void store_solution();
  
  // parameter tuning methods
  void guide(WalksatExpArray& vars, 
	     WalksatIntArray& vals,
	     WalksatDoubleArray& probs) {}
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

  //int solveDimacs(const char*);
};


#endif

