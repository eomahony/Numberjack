#ifndef MipWrapper_H
#define MipWrapper_H

//#define _DEBUGWRAP 1

#ifdef _DEBUGWRAP
  #define DBG(fmt, args...) printf("%s:%d "fmt,__FILE__,__LINE__,args)
#else
  #define DBG(fmt, args...) ;
#endif

#include <vector>
#include <stdlib.h>
#include <climits>

const int UNSAT     =  0;
const int SAT       =  1;
const int UNKNOWN   =  2;

const int LUBY      =  0;
const int GEOMETRIC =  1;
const int GLUBY     =  2;  

/**
   Expression (Used to encode variables & constraints)
*/
class MipWrapperSolver;
class MipWrapper_Expression;

typedef struct linear_con_arg{
  MipWrapper_Expression *expr;
  double coef;
  double offset;
} LINEAR_ARG;

class MipWrapper_Expression{    
public:
  int _ident;
  int nbj_ident;
  
  MipWrapperSolver *_solver;
  MipWrapper_Expression **_expr_encoding;
  
  void **_encoding;
  void *_var;

  double _lower;
  double _upper;
  double _coef;
  
  virtual void leq(double value, MipWrapperSolver* solver);
  virtual void lt( double value, MipWrapperSolver* solver);
  virtual void geq(double value, MipWrapperSolver* solver);
  virtual void gt( double value, MipWrapperSolver* solver);
  virtual void eq( double value, MipWrapperSolver* solver);
  virtual void neq(double value, MipWrapperSolver* solver);
  
  bool _continuous;
  
  int get_size();
  int get_max();
  int get_min();
  
  bool has_been_added() const;
  void initialise(bool c);
  void display();

  MipWrapper_Expression();
  virtual ~MipWrapper_Expression();

  virtual void encode(MipWrapperSolver* solver);
  virtual MipWrapper_Expression* add(MipWrapperSolver *solver,
				     bool top_level);
  virtual LINEAR_ARG* for_linear();
  virtual int for_linear_size();
  
  virtual double get_whatever_value();
};

/**
   Array of expressions (for nary constraints)
*/
template<typename T> class MipWrapperArray{
  private:
    std::vector< T > _array;
  public:
    MipWrapperArray() {}
    virtual ~MipWrapperArray() {}
    int size() { return _array.size(); }
    void add(T arg) { _array.push_back(arg); }
    void set_item(const int i, T arg){ _array[i] = arg; } 
    T& get_item(const int i) { return _array[i]; }
};

typedef MipWrapperArray< MipWrapper_Expression* > MipWrapperExpArray;
typedef MipWrapperArray< int > MipWrapperIntArray;
typedef MipWrapperArray< double > MipWrapperDoubleArray;

class MipWrapper_FloatVar : public MipWrapper_Expression {
public:
  MipWrapper_FloatVar();
  MipWrapper_FloatVar(const double lb, const double ub);
  MipWrapper_FloatVar(const double lb, const double ub, const int ident);
  virtual double get_value();
  virtual double get_whatever_value();
};

class MipWrapper_IntVar : public MipWrapper_Expression {
  private:
    MipWrapperIntArray _values;
    bool _has_holes_in_domain;
  public:
    //IntVar(const int nval);
    MipWrapper_IntVar(const int lb, const int ub);
    MipWrapper_IntVar(const int lb, const int ub, const int ident);
    MipWrapper_IntVar(MipWrapperIntArray& values, const int ident);
    virtual void encode(MipWrapperSolver* solver);
    virtual MipWrapper_Expression* add(MipWrapperSolver* solver, bool top_level);
    virtual int get_value();
    virtual double get_whatever_value();
};

class MipWrapper_Flow : public MipWrapper_FloatVar{
protected:
  int max_val;
  int min_val;

  MipWrapperExpArray _vars;

  int *card_lb;
  int *card_ub;

public:
  MipWrapper_Flow(MipWrapperExpArray& vars);
  MipWrapper_Flow();
  virtual ~MipWrapper_Flow();
  virtual MipWrapper_Expression* add(MipWrapperSolver *solver, bool top_level);
  virtual void initbounds();
  void addVar( MipWrapper_Expression* v );
};

class MipWrapper_AllDiff : public MipWrapper_Flow{
public:
  MipWrapper_AllDiff(MipWrapperExpArray& vars);
  MipWrapper_AllDiff(MipWrapper_Expression* arg1, MipWrapper_Expression* arg2);
  virtual ~MipWrapper_AllDiff();
};

class MipWrapper_Gcc : public MipWrapper_Flow{
public:
  MipWrapper_Gcc(MipWrapperExpArray& vars,
		 MipWrapperIntArray& vals,
		 MipWrapperIntArray& lb,
		 MipWrapperIntArray& ub);
  virtual ~MipWrapper_Gcc();
};

class MipWrapper_Sum: public MipWrapper_FloatVar{
public:
  double _offset;
  MipWrapperExpArray _vars;
  MipWrapperIntArray _weights;
  
  MipWrapper_Sum(MipWrapperExpArray& vars, 
	   MipWrapperIntArray& weights, 
	   const int offset=0);
  MipWrapper_Sum(MipWrapper_Expression* arg1, 
	   MipWrapper_Expression* arg2, 
	   MipWrapperIntArray& weights,
	   const int offset=0 );
  MipWrapper_Sum(MipWrapper_Expression* arg, 
	   MipWrapperIntArray& weights,
	   const int offset=0 );
  MipWrapper_Sum();
  
  virtual void leq(double value, MipWrapperSolver* solver);
  virtual void geq(double value, MipWrapperSolver* solver);
  virtual void eq( double value, MipWrapperSolver* solver);
  virtual void neq(double value, MipWrapperSolver* solver);
  
  virtual void encode(MipWrapperSolver* solver);
  
  void initialise();
  void addVar( MipWrapper_Expression* v );
  void addWeight( const int w );
  void set_rhs( const int k );

  virtual ~MipWrapper_Sum();
  virtual MipWrapper_Expression* add(MipWrapperSolver *solver, bool top_level);
  virtual LINEAR_ARG *for_linear();
  virtual int for_linear_size();
  
  virtual double get_value();
};

class MipWrapper_add: public MipWrapper_Sum{
public:
  MipWrapper_add( MipWrapper_Expression* arg1, MipWrapper_Expression* arg2 );
  MipWrapper_add( MipWrapper_Expression* arg1, const int arg2 );
  virtual ~MipWrapper_add();
};

class MipWrapper_sub: public MipWrapper_Sum{
public:
  MipWrapper_sub( MipWrapper_Expression* arg1, MipWrapper_Expression* arg2 );
  MipWrapper_sub( MipWrapper_Expression* arg1, const int arg2 );
  virtual ~MipWrapper_sub();
};


class MipWrapper_binop: public MipWrapper_FloatVar{
protected:
  MipWrapper_Expression *_vars[2];
  MipWrapper_Expression *_repr;
  bool _is_proper_coef;
  double _rhs;
public:
  MipWrapper_binop(MipWrapper_Expression *var1, MipWrapper_Expression *var2);
  MipWrapper_binop(MipWrapper_Expression *var1, double rhs);
  virtual ~MipWrapper_binop();
  virtual MipWrapper_Expression* add(MipWrapperSolver *solver, bool top_level) = 0;
};

class MipWrapper_NoOverlap: public MipWrapper_binop
{
protected:
  MipWrapperIntArray _coefs;
public:
  MipWrapper_NoOverlap(MipWrapper_Expression *var1, MipWrapper_Expression *var2,
		       int d1, int d2);
  virtual ~MipWrapper_NoOverlap();
  virtual MipWrapper_Expression* add(MipWrapperSolver *solver, bool top_level);
};

class MipWrapper_Precedence: public MipWrapper_binop{
protected:
    MipWrapperIntArray _coefs;
public:
    MipWrapper_Precedence(MipWrapper_Expression *var1,
			  MipWrapper_Expression *var2,
			  MipWrapperIntArray &coefs);
    virtual ~MipWrapper_Precedence();
    virtual MipWrapper_Expression* add(MipWrapperSolver *solver,
				       bool top_level);
};

class MipWrapper_eq: public MipWrapper_binop{
public:
  MipWrapper_eq(MipWrapper_Expression *var1, MipWrapper_Expression *var2);
  MipWrapper_eq(MipWrapper_Expression *var1, double rhs);
  virtual ~MipWrapper_eq();
  virtual MipWrapper_Expression* add(MipWrapperSolver *solver, bool top_level);
};

class MipWrapper_ne: public MipWrapper_binop{
public:
  MipWrapper_ne(MipWrapper_Expression *var1, MipWrapper_Expression *var2);
  MipWrapper_ne(MipWrapper_Expression *var1, double rhs);
  virtual ~MipWrapper_ne();
  virtual MipWrapper_Expression* add(MipWrapperSolver *solver, bool top_level);
};

class MipWrapper_le: public MipWrapper_binop{
public:
  MipWrapper_le(MipWrapper_Expression *var1, MipWrapper_Expression *var2);
  MipWrapper_le(MipWrapper_Expression *var1, double rhs);
  virtual ~MipWrapper_le();
  virtual MipWrapper_Expression* add(MipWrapperSolver *solver, bool top_level);
};

class MipWrapper_ge: public MipWrapper_binop{
public:
  MipWrapper_ge(MipWrapper_Expression *var1, MipWrapper_Expression *var2);
  MipWrapper_ge(MipWrapper_Expression *var1, double rhs);
  virtual ~MipWrapper_ge();
  virtual MipWrapper_Expression* add(MipWrapperSolver *solver, bool top_level);
};

class MipWrapper_lt: public MipWrapper_binop{
public:
  MipWrapper_lt(MipWrapper_Expression *var1, MipWrapper_Expression *var2);
  MipWrapper_lt(MipWrapper_Expression *var1, double rhs);
  virtual ~MipWrapper_lt();
  virtual MipWrapper_Expression* add(MipWrapperSolver *solver, bool top_level);
};

class MipWrapper_gt: public MipWrapper_binop{
public:
  MipWrapper_gt(MipWrapper_Expression *var1, MipWrapper_Expression *var2);
  MipWrapper_gt(MipWrapper_Expression *var1, double rhs);
  virtual ~MipWrapper_gt();
  virtual MipWrapper_Expression* add(MipWrapperSolver *solver, bool top_level);
};

class MipWrapper_not: public MipWrapper_binop{
public:
  MipWrapper_not(MipWrapper_Expression *var1);
  virtual ~MipWrapper_not();
  virtual MipWrapper_Expression* add(MipWrapperSolver *solver, bool top_level);
};

class MipWrapper_and: public MipWrapper_binop{
public:
  MipWrapper_and(MipWrapper_Expression *var1, MipWrapper_Expression *var2);
  virtual ~MipWrapper_and();
  virtual MipWrapper_Expression* add(MipWrapperSolver *solver, bool top_level);
};

class MipWrapper_or: public MipWrapper_binop{
public:
  MipWrapper_or(MipWrapper_Expression *var1, MipWrapper_Expression *var2);
  virtual ~MipWrapper_or();
  virtual MipWrapper_Expression* add(MipWrapperSolver *solver, bool top_level);
};

class MipWrapper_Minimise: public MipWrapper_Expression{
private:  
   MipWrapper_Expression* _obj;
public:
    MipWrapper_Minimise(MipWrapper_Expression *var1);
    virtual ~MipWrapper_Minimise();
    virtual MipWrapper_Expression* add(MipWrapperSolver *solver, bool top_level);
};

class MipWrapper_Maximise: public MipWrapper_Expression{
private:  
  MipWrapper_Expression* _obj;
public:
  MipWrapper_Maximise(MipWrapper_Expression *var1);
  virtual ~MipWrapper_Maximise();
  virtual MipWrapper_Expression* add(MipWrapperSolver *solver, bool top_level);

};

class LinearConstraint{
  
  public:
    /**
     * The variables in the scope of the constraint
     */
    std::vector<MipWrapper_Expression*> _variables;
    
    /**
     * The coefficients in the constraint
     */
    std::vector<double> _coefficients;
    
    /**
     * Left hand side of the leq
     */
    double _lhs;
    
    /**
     * Right hand side of the leq
     */
    double _rhs;
    
    /**
     * Creates an empty linear constraint object
     */
    LinearConstraint(double lhs, double rhs);
    
    /**
     * Destructor
     */
    virtual ~LinearConstraint();
    
    /**
     * Add in a normal coefficient
     */
    virtual void add_coef(MipWrapper_Expression* expr,
			  double coef,
			  bool use_encoding=true);
    
    virtual void add_coef(LINEAR_ARG* arg_struct,
			  double coef,
			  bool use_encoding);
    
    /**
     * Prints the linear constraint
     */
    void display();
  
};


/**
   The solver itself
*/
class MipWrapperSolver{
private:
  std::vector< MipWrapper_Expression*> _exprs;
  std::vector< MipWrapperIntArray * > _int_array;
  std::vector< MipWrapperExpArray * > _var_array;
  
  int _verbosity;

public:
  std::vector< LinearConstraint* > _constraints;
  LinearConstraint* _obj;
  double _obj_coef;
  
  int var_counter;
  
  int get_lin_cons_size();

  MipWrapperSolver();
  virtual ~MipWrapperSolver();

  // add an expression, in the case of a tree of expressions,
  // each node of the tree is added separately, depth first.
  virtual void add(MipWrapper_Expression* arg);
  
  void add_expr(MipWrapper_Expression *expr);
  void add_int_array(MipWrapperIntArray *arr);
  void add_var_array(MipWrapperExpArray *arr);

  // used to initialise search on a given subset of variables
  virtual void initialise(MipWrapperExpArray& arg);
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
  
  // parameter tuning methods
  virtual void setHeuristic(const char* var_heuristic,
			    const char* val_heuristic,
			    const int rand);
  virtual void setFailureLimit(const int cutoff);  
  virtual void setTimeLimit(const int cutoff);
  virtual void setNodeLimit(const int cutoff);
  virtual void setVerbosity(const int degree);
  virtual void setRandomized(const int degree);
  virtual void setRandomSeed(const int seed);

  // statistics methods
  virtual bool is_sat();
  virtual bool is_unsat();
  virtual void printStatistics();
  virtual int getBacktracks();
  virtual int getNodes();
  virtual int getFailures();
  virtual int getChecks();
  virtual int getPropags();
  virtual double getTime();
  
  // Value pass back stuff
  virtual double get_value(void *ptr);

};

#endif

