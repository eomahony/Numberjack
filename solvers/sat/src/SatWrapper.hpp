
#ifndef SATWRAPPER_H
#define SATWRAPPER_H


#include <vector>
#include <iostream>

//#include "Literals.h"
#include "Solver.h"



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
class SatWrapperArray
{
private:
  std::vector< T > _array;

public:
  SatWrapperArray() {}
  virtual ~SatWrapperArray() {}
  int size() const { return _array.size(); }
  void add(T arg) { _array.push_back(arg); }
  T get_item(const int i) const { return _array[i]; }
  void set_item(const int i, T item) { _array[i]=item; }
};

typedef SatWrapperArray< int > SatWrapperIntArray;
typedef SatWrapperArray< double > SatWrapperDoubleArray;







class SatWrapper_Expression;
class SatWrapperSolver;
class AbstractDomain {
public:
  SatWrapper_Expression *owner;

  AbstractDomain(SatWrapper_Expression *o) {owner = o;}

  virtual int getval(int idx) const = 0;
  virtual int getmin() const = 0;
  virtual int getmax() const = 0;
  virtual int getsize() const = 0;
  virtual int contain(const int v) const = 0;
  virtual int next(const int v, const int idx=-1) const = 0;

  virtual void encode(SatWrapperSolver *solver) {}

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

  OffsetDomain(SatWrapper_Expression *os, AbstractDomain *d, const int of);

  virtual int getval(int idx) const;
  virtual int getmin() const;
  virtual int getmax() const;
  virtual int getsize() const {return _dom_ptr->getsize();}
  virtual int contain(const int v) const {return _dom_ptr->contain(v-offset);}
  virtual int next(const int v, const int idx=-1) const {return _dom_ptr->next(v-offset,-1);}

  virtual Lit less_or_equal(const int value, const int index) const;
  virtual Lit equal(const int value, const int index) const;

  virtual std::ostream& display(std::ostream& o) const;
};

class FactorDomain : public AbstractDomain {
public:
  AbstractDomain *_dom_ptr;
  int factor;

  FactorDomain(SatWrapper_Expression *o, AbstractDomain *d, const int f);

  virtual int getval(int idx) const;
  virtual int getmin() const;
  virtual int getmax() const;
  virtual int getsize() const {return _dom_ptr->getsize();}
  virtual int contain(const int v) const {return (!(v%factor) && _dom_ptr->contain(v/factor));}
  virtual int next(const int v, const int idx=-1) const {return factor * (_dom_ptr->next(v/factor, -1));}

  virtual Lit less_or_equal(const int value, const int index) const;
  virtual Lit equal(const int value, const int index) const;

  virtual std::ostream& display(std::ostream& o) const;
};

class EqDomain : public AbstractDomain {
public:
  AbstractDomain *_dom_ptr;
  int value;
  int spin;

  EqDomain(SatWrapper_Expression *o, AbstractDomain *d, const int v, const int s);

  virtual int getval(int idx) const {assert(idx >=0 && idx <= 1); return idx;}
  virtual int getmin() const {return 0;}
  virtual int getmax() const {return 1;}
  virtual int getsize() const {return 2;}
  virtual int contain(const int v) const {return v == 0 || v == 1;}
  virtual int next(const int v, const int idx=-1) const {return (v ? v : 1);}

  virtual Lit less_or_equal(const int value, const int index) const;
  virtual Lit equal(const int value, const int index) const;

  virtual std::ostream& display(std::ostream& o) const;
};

class LeqDomain : public AbstractDomain {
public:
  AbstractDomain *_dom_ptr;
  int bound;
  int spin;

  LeqDomain(SatWrapper_Expression *o, AbstractDomain *d, const int b, const int s);

  virtual int getval(int idx) const {assert(idx >=0 && idx <= 1); return idx;}
  virtual int getmin() const {return 0;}
  virtual int getmax() const {return 1;}
  virtual int getsize() const {return 2;}
  virtual int contain(const int v) const {return v == 0 || v == 1;}
  virtual int next(const int v, const int idx=-1) const {return (v ? v : 1);}

  virtual Lit less_or_equal(const int value, const int index) const;
  virtual Lit equal(const int value, const int index) const;

  virtual std::ostream& display(std::ostream& o) const;
};

class ConstantDomain : public AbstractDomain {
public:
  int value;

  ConstantDomain(SatWrapper_Expression *o, const int v) : AbstractDomain(o) {value = v;}

  virtual int getval(int idx) const { assert(idx == 0); return value; }
  virtual int getmin() const {return value;}
  virtual int getmax() const {return value;}
  virtual int getsize() const {return 1;}
  virtual int contain(const int v) const {return v == value;}
  virtual int next(const int v, const int idx=-1) const {return v;}

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

  DomainEncoding(SatWrapper_Expression *o);
  DomainEncoding(SatWrapper_Expression *o, const int nval);
  DomainEncoding(SatWrapper_Expression *o, const int lb, const int ub);
  DomainEncoding(SatWrapper_Expression *o, SatWrapperIntArray& vals);

  virtual ~DomainEncoding();


  virtual void print_lit(Lit p, const int type) const;


  virtual int getval(int idx) const {
    if(_values) return _values[idx%_size];
    else return _lower+idx;
  }

  virtual int getmin() const {return _lower;}
  virtual int getmax() const {return _upper;}
  virtual int getsize() const {return _size;}
  virtual int contain(const int v) const;
  virtual int next(const int v, const int idx=-1) const;

  void encode(SatWrapperSolver *solver);

  virtual Lit less_or_equal(const int value, const int index) const;
  virtual Lit equal(const int value, const int index) const;

  virtual std::ostream& display(std::ostream& o) const;

};





/**
   Expression (Used to encode variables & constraints)
*/
class SatWrapperSolver;
class SatWrapper_Expression
{
    
public:
  
  SatWrapperSolver *_solver;

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

  SatWrapper_Expression();
  SatWrapper_Expression(const int nval);
  SatWrapper_Expression(const int lb, const int ub);
  SatWrapper_Expression(SatWrapperIntArray& vals);
  virtual ~SatWrapper_Expression();

  virtual SatWrapper_Expression* add(SatWrapperSolver *solver, bool top_level);
};

class SatWrapper_IntVar : public SatWrapper_Expression
{

public:
  
  SatWrapper_IntVar() : SatWrapper_Expression() {}
  SatWrapper_IntVar(const int lb, const int ub, const int ident) : SatWrapper_Expression(lb, ub) {nbj_ident = ident;}
  SatWrapper_IntVar(SatWrapperIntArray& vals, const int ident) : SatWrapper_Expression(vals) {nbj_ident = ident;}

};


typedef SatWrapperArray< SatWrapper_Expression* > SatWrapperExpArray;

class SatWrapper_AllDiff : public SatWrapper_Expression
{
private:
  SatWrapperExpArray _vars;
  std::vector< SatWrapper_Expression* > _clique;

public:
  SatWrapper_AllDiff(SatWrapperExpArray& vars);
  SatWrapper_AllDiff(SatWrapper_Expression* arg1, SatWrapper_Expression* arg2);
  void addVar( SatWrapper_Expression* v );
  virtual ~SatWrapper_AllDiff();

  virtual SatWrapper_Expression* add(SatWrapperSolver *solver, bool top_level);
};


class SatWrapper_Sum: public SatWrapper_Expression
{
private:
  SatWrapper_Expression *_self;
  int _offset;
  SatWrapperExpArray _vars;
  SatWrapperIntArray _weights;
  std::vector< SatWrapper_Expression* > _subsum;
    
public:
  SatWrapper_Sum(SatWrapperExpArray& vars, SatWrapperIntArray& weights, const int offset=0);
  SatWrapper_Sum( SatWrapper_Expression* arg1, 
	       SatWrapper_Expression* arg2, 
	       SatWrapperIntArray& weights,
	       const int offset );
  SatWrapper_Sum( SatWrapper_Expression* arg, 
	       SatWrapperIntArray& weights,
	       const int offset );
  SatWrapper_Sum();
  void initialise();


  void addVar( SatWrapper_Expression* v );
  void addWeight( const int w );
  void set_rhs( const int k );

  virtual ~SatWrapper_Sum();
  virtual SatWrapper_Expression* add(SatWrapperSolver *solver, bool top_level);
};


class SatWrapper_binop: public SatWrapper_Expression
{
protected:
  SatWrapper_Expression *_vars[2];
  int _rhs;
  
public:
  int arity() { return 1+(_vars[1]==NULL); }
  SatWrapper_binop(SatWrapper_Expression *var1, SatWrapper_Expression *var2);
  SatWrapper_binop(SatWrapper_Expression *var1, int rhs);
  virtual ~SatWrapper_binop();

  virtual SatWrapper_Expression* add(SatWrapperSolver *solver, bool top_level) = 0;
};


class SatWrapper_add : public SatWrapper_binop
{
public:
  SatWrapper_add( SatWrapper_Expression* arg1, SatWrapper_Expression* arg2 );
  SatWrapper_add( SatWrapper_Expression* arg1, const int arg2 );
  virtual ~SatWrapper_add();
  virtual SatWrapper_Expression* add(SatWrapperSolver *solver, bool top_level);
};

class SatWrapper_mul : public SatWrapper_binop
{
public:
  SatWrapper_mul( SatWrapper_Expression* arg1, SatWrapper_Expression* arg2 );
  SatWrapper_mul( SatWrapper_Expression* arg1, const int arg2 );
  virtual ~SatWrapper_mul();
  virtual SatWrapper_Expression* add(SatWrapperSolver *solver, bool top_level);
};

class SatWrapper_or: public SatWrapper_binop
{
public:
  SatWrapper_or(SatWrapper_Expression *var1, SatWrapper_Expression *var2);
  SatWrapper_or(SatWrapper_Expression *var1, int rhs);
  virtual ~SatWrapper_or();
  virtual SatWrapper_Expression* add(SatWrapperSolver *solver, bool top_level);
};

class SatWrapper_and: public SatWrapper_binop
{
public:
  SatWrapper_and(SatWrapper_Expression *var1, SatWrapper_Expression *var2);
  SatWrapper_and(SatWrapper_Expression *var1, int rhs);
  virtual ~SatWrapper_and();
  virtual SatWrapper_Expression* add(SatWrapperSolver *solver, bool top_level);
};


bool processClause(std::vector<Lit>& cl_in, std::vector<Lit>& cl_out);


class SatWrapper_eq: public SatWrapper_binop
{
public:
  SatWrapper_eq(SatWrapper_Expression *var1, SatWrapper_Expression *var2);
  SatWrapper_eq(SatWrapper_Expression *var1, int rhs);
  virtual ~SatWrapper_eq();
  virtual SatWrapper_Expression* add(SatWrapperSolver *solver, bool top_level);
};

class SatWrapper_ne: public SatWrapper_binop
{
public:
  SatWrapper_ne(SatWrapper_Expression *var1, SatWrapper_Expression *var2);
  SatWrapper_ne(SatWrapper_Expression *var1, int rhs);
  virtual ~SatWrapper_ne();
  virtual SatWrapper_Expression* add(SatWrapperSolver *solver, bool top_level);
};

class SatWrapper_le: public SatWrapper_binop
{
public:
  SatWrapper_le(SatWrapper_Expression *var1, SatWrapper_Expression *var2);
  SatWrapper_le(SatWrapper_Expression *var1, int rhs);
  virtual ~SatWrapper_le();
  virtual SatWrapper_Expression* add(SatWrapperSolver *solver, bool top_level);
};

class SatWrapper_ge: public SatWrapper_binop
{
public:
  SatWrapper_ge(SatWrapper_Expression *var1, SatWrapper_Expression *var2);
  SatWrapper_ge(SatWrapper_Expression *var1, int rhs);
  virtual ~SatWrapper_ge();
  virtual SatWrapper_Expression* add(SatWrapperSolver *solver, bool top_level);
};

class SatWrapper_lt: public SatWrapper_binop
{
public:
  SatWrapper_lt(SatWrapper_Expression *var1, SatWrapper_Expression *var2);
  SatWrapper_lt(SatWrapper_Expression *var1, int rhs);
  virtual ~SatWrapper_lt();
  virtual SatWrapper_Expression* add(SatWrapperSolver *solver, bool top_level);
};

class SatWrapper_gt: public SatWrapper_binop
{
public:
  SatWrapper_gt(SatWrapper_Expression *var1, SatWrapper_Expression *var2);
  SatWrapper_gt(SatWrapper_Expression *var1, int rhs);
  virtual ~SatWrapper_gt();
  virtual SatWrapper_Expression* add(SatWrapperSolver *solver, bool top_level);
};

class SatWrapper_Minimise : public SatWrapper_Expression
{
protected:
  SatWrapper_Expression *_obj;
public:
  SatWrapper_Minimise(SatWrapper_Expression *var1);
  virtual ~SatWrapper_Minimise();
  virtual SatWrapper_Expression* add(SatWrapperSolver *solver, bool top_level);
};

class SatWrapper_Maximise : public SatWrapper_Expression
{
protected:
  SatWrapper_Expression *_obj;
public:
  SatWrapper_Maximise(SatWrapper_Expression *var1);
  virtual ~SatWrapper_Maximise();
  virtual SatWrapper_Expression* add(SatWrapperSolver *solver, bool top_level);
};


/**
   The solver itself
*/
class SatWrapperSolver 
{

public:

  SatWrapper_Expression *minimise_obj;
  SatWrapper_Expression *maximise_obj;

  // repository for all expressions
  std::vector< SatWrapper_Expression* > _expressions;
  std::vector< SatWrapper_Expression* > _variables;
  // link each atom to its domain
  std::vector< DomainEncoding* > _atom_to_domain;
  std::vector< int > _atom_to_type;

  std::vector< std::vector<Lit> > clause_base;
  unsigned int current;
  
  int clause_limit;

  int *cp_model;

  SatWrapperSolver();
  virtual ~SatWrapperSolver();
  
  int get_cb_size();

  // add an expression, in the case of a tree of expressions,
  // each node of the tree is added separately, depth first.
  virtual void add(SatWrapper_Expression* arg);
  virtual int declare(SatWrapper_Expression* arg, const bool type);
  virtual int create_atom(DomainEncoding* dom, const int type);

  virtual void addClause(std::vector<Lit>& cl);
  virtual void validate();
  virtual void displayClause(std::vector<Lit>& cl); 
  virtual void displayLiteral(Lit p);

  virtual lbool truth_value(Lit x);

  // used to initialise search on a given subset of variables
  virtual void initialise(SatWrapperExpArray& arg);
  // initialise the solver before solving (no more calls to add after this)
  virtual void initialise();

  // solving methods
  virtual int solve();
  virtual int solveAndRestart(const int policy = GEOMETRIC, 
			      const unsigned int base = 32, 
			      const double factor = 1.3333333,
			      const double decay = 0.0);
  virtual int startNewSearch();
  virtual int getNextSolution();
  virtual int sacPreprocess(const int type);

  virtual void reset(bool full);
  virtual bool propagate();
  virtual void save();
  virtual void post(const char* op, SatWrapper_Expression* x, int v);
  virtual bool undo(const int nlevel);
  virtual void deduce();
  virtual bool branch_right();

  virtual void store_solution();
  
  virtual void setClauseLimit(int limit);
  
  // parameter tuning methods
  virtual void guide(SatWrapperExpArray& vars, 
		     SatWrapperIntArray& vals,
		     SatWrapperDoubleArray& probs) {}
  virtual void setHeuristic(const char* var_heuristic, const char* val_heuristic, const int rand);
  virtual void setFailureLimit(const int cutoff);  
  virtual void setNodeLimit(const int cutoff);  
  virtual void setTimeLimit(const double cutoff);
  virtual void setVerbosity(const int degree);
  virtual void setRandomized(const int degree);
  virtual void setRandomSeed(const int seed);

  // statistics methods
  virtual bool is_sat();
  virtual bool is_opt();
  virtual bool is_unsat();
  virtual void printStatistics();
  virtual int getBacktracks();
  virtual int getNodes();
  virtual int getFailures();
  virtual int getChecks();
  virtual int getPropags();
  virtual double getTime();
};


#endif

