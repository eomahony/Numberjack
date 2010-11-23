
/*
  Mistral is a constraint satisfaction and optimisation library
  Copyright (C) 2003-2005  Emmanuel Hebrard
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

  The author can be contacted electronically at 
  ehebrard@cse.unsw.edu.au.
*/

/** \file mistral_pyt.h
    \brief Header for the PYTHON Wrapper.
*/


#ifndef _PYTHON_H
#define _PYTHON_H


#include <mistral_mod.h>
#include <mistral_glo.h>

//#include "../models/src/xcsp/XMLParser_libxml2.hh"
//#include "../models/src/xcsp/MistralCallback.hh"


  class VarValuation {
    
  public:

    int value;
    char type;
    
    VarValuation(const int v, const int t) {value=v; type=t;}
    virtual ~VarValuation() {}

  };


/**
   Array of expressions (for nary constraints)
*/
template<class T>
class MistralArray
{
private:
  std::vector< T > _array;

public:
  MistralArray() {}
  virtual ~MistralArray() {}
  int size() { return _array.size(); }
  void add(T arg) { _array.push_back(arg); }
  T& get_item(const int i) { return _array[i]; }
  void set_item(const int i, T item) { _array[i]=item; }
};

typedef MistralArray< int > MistralIntArray;
typedef MistralArray< double > MistralDoubleArray;

/**
   Expression (Used to encode variables & constraints)
*/
class MistralSolver;
class Mistral_Expression
{
    
public:

  int nbj_ident;
  MistralSolver *_solver;
  BuildObject *_self;


  bool has_been_added() const;
  void initialise();

  Mistral_Expression();
  Mistral_Expression(BuildObject *x);
  Mistral_Expression(const int nval);
  Mistral_Expression(const int lb, const int ub);
  Mistral_Expression(MistralIntArray& vals);
  virtual ~Mistral_Expression();

  const char* get_type() const;
  int get_arity() const;
  Mistral_Expression *get_child(const int i);
  int get_id() const;

  int next(int v);
  int getVariableId() const;
  int get_value() const;
  int get_size() const;
  int get_min() const;
  int get_max() const;
  bool contain(const int v) const;

  virtual Mistral_Expression* add(MistralSolver *solver, bool top_level);

  void print_python() const ;
};

class Mistral_IntVar : public Mistral_Expression
{

public:
  
  Mistral_IntVar() : Mistral_Expression() {}
  //Mistral_IntVar(const int nval) : Mistral_Expression(nval) {}
  Mistral_IntVar(const int lb, const int ub, const int ident) : Mistral_Expression(lb, ub) {nbj_ident = ident;}
  Mistral_IntVar(MistralIntArray& vals, const int ident) : Mistral_Expression(vals) {nbj_ident = ident;}

};

typedef MistralArray< Mistral_Expression* > MistralExpArray;

class Mistral_Min : public Mistral_Expression
{
private:
  MistralExpArray _vars;

public:
  Mistral_Min(MistralExpArray& vars);
  Mistral_Min(Mistral_Expression *var1, Mistral_Expression *var2);
  virtual ~Mistral_Min();

  virtual Mistral_Expression* add(MistralSolver *solver, bool top_level);
};

class Mistral_Max : public Mistral_Expression
{
private:
  MistralExpArray _vars;

public:
  Mistral_Max(MistralExpArray& vars);
  Mistral_Max(Mistral_Expression *var1, Mistral_Expression *var2);
  virtual ~Mistral_Max();

  virtual Mistral_Expression* add(MistralSolver *solver, bool top_level);
};

class Mistral_AllDiff : public Mistral_Expression
{
private:
  MistralExpArray _vars;

public:
  Mistral_AllDiff(MistralExpArray& vars);
  Mistral_AllDiff(Mistral_Expression *var1, Mistral_Expression *var2);
  virtual ~Mistral_AllDiff();

  virtual Mistral_Expression* add(MistralSolver *solver, bool top_level);
};

class Mistral_Table : public Mistral_Expression
{
private:
  MistralExpArray _vars;
  MistralIntArray _tuples;
  int spin;

public:
  Mistral_Table(MistralExpArray& vars, MistralIntArray& tuples, const char* type); 
  Mistral_Table(Mistral_Expression *var1, Mistral_Expression *var2, MistralIntArray& tuples, const char* type);
  virtual ~Mistral_Table();
  
  virtual void add(MistralIntArray& tuple);
  virtual Mistral_Expression* add(MistralSolver *solver, bool top_level);
};

class Mistral_Gcc : public Mistral_Expression
{
private:
  MistralExpArray _vars;
  MistralIntArray _vals;
  MistralIntArray _lb_card;
  MistralIntArray _ub_card;

public:
  Mistral_Gcc(MistralExpArray& vars,
	      MistralIntArray& vals,
	      MistralIntArray& lb_card,
	      MistralIntArray& ub_card);
  virtual ~Mistral_Gcc();

  virtual Mistral_Expression* add(MistralSolver *solver, bool top_level);
};

class Mistral_Element : public Mistral_Expression
{
private:
  MistralExpArray _vars;

public:
  Mistral_Element(MistralExpArray& vars);
  virtual ~Mistral_Element();

  virtual Mistral_Expression* add(MistralSolver *solver, bool top_level);
};

class Mistral_LeqLex : public Mistral_Expression
{
private:
  MistralExpArray _vars;

public:
  Mistral_LeqLex(MistralExpArray& vars);
  virtual ~Mistral_LeqLex();

  virtual Mistral_Expression* add(MistralSolver *solver, bool top_level);
};

class Mistral_LessLex : public Mistral_Expression
{
private:
  MistralExpArray _vars;

public:
  Mistral_LessLex(MistralExpArray& vars);
  virtual ~Mistral_LessLex();

  virtual Mistral_Expression* add(MistralSolver *solver, bool top_level);
};


class Mistral_Sum: public Mistral_Expression
{

private:

  MistralExpArray _vars;
  MistralIntArray _weights;
    
public:

  Mistral_Sum( MistralExpArray& vars, 
	       MistralIntArray& weights, 
	       const int offset=0);
  Mistral_Sum( Mistral_Expression* arg1, 
	       Mistral_Expression* arg2, 
	       MistralIntArray& weights,
	       const int offset );
  Mistral_Sum( Mistral_Expression* arg, 
	       MistralIntArray& weights,
	       const int offset );
  Mistral_Sum();

  void addVar( Mistral_Expression* v );
  void addWeight( const int w );

  virtual ~Mistral_Sum();
  virtual Mistral_Expression* add(MistralSolver *solver, bool top_level);
};


class Mistral_binop: public Mistral_Expression
{
protected:
  Mistral_Expression *_vars[2];
  int _constant;
  
public:
  int arity() { return 1+(_vars[1]==NULL); }
  Mistral_binop(Mistral_Expression *var1, Mistral_Expression *var2);
  Mistral_binop(Mistral_Expression *var1, int int_arg);
  virtual ~Mistral_binop();

  virtual Mistral_Expression* add(MistralSolver *solver, bool top_level) = 0;
};

class Mistral_mul: public Mistral_binop
{
public:
  Mistral_mul(Mistral_Expression *var1, Mistral_Expression *var2);
  Mistral_mul(Mistral_Expression *var1, int int_arg);
  virtual ~Mistral_mul();
  virtual Mistral_Expression* add(MistralSolver *solver, bool top_level);
};

class Mistral_div: public Mistral_binop
{
public:
  Mistral_div(Mistral_Expression *var1, Mistral_Expression *var2);
  Mistral_div(Mistral_Expression *var1, int int_arg);
  virtual ~Mistral_div();
  virtual Mistral_Expression* add(MistralSolver *solver, bool top_level);
};

class Mistral_mod: public Mistral_binop
{
public:
  Mistral_mod(Mistral_Expression *var1, Mistral_Expression *var2);
  Mistral_mod(Mistral_Expression *var1, int int_arg);
  virtual ~Mistral_mod();
  virtual Mistral_Expression* add(MistralSolver *solver, bool top_level);
};

class Mistral_and: public Mistral_binop
{
public:
  Mistral_and(Mistral_Expression *var1, Mistral_Expression *var2);
  Mistral_and(Mistral_Expression *var1, int int_arg);
  virtual ~Mistral_and();
  virtual Mistral_Expression* add(MistralSolver *solver, bool top_level);
};

class Mistral_or: public Mistral_binop
{
public:
  Mistral_or(Mistral_Expression *var1, Mistral_Expression *var2);
  Mistral_or(Mistral_Expression *var1, int int_arg);
  virtual ~Mistral_or();
  virtual Mistral_Expression* add(MistralSolver *solver, bool top_level);
};

class Mistral_eq: public Mistral_binop
{
public:
  Mistral_eq(Mistral_Expression *var1, Mistral_Expression *var2);
  Mistral_eq(Mistral_Expression *var1, int int_arg);
  virtual ~Mistral_eq();
  virtual Mistral_Expression* add(MistralSolver *solver, bool top_level);
};

class Mistral_ne: public Mistral_binop
{
public:
  Mistral_ne(Mistral_Expression *var1, Mistral_Expression *var2);
  Mistral_ne(Mistral_Expression *var1, int int_arg);
  virtual ~Mistral_ne();
  virtual Mistral_Expression* add(MistralSolver *solver, bool top_level);
};

class Mistral_NoOverlap: public Mistral_binop
{
private:
  int _bonstant;
public:
  Mistral_NoOverlap(Mistral_Expression *var1, Mistral_Expression *var2, int d1, int d2);
  virtual ~Mistral_NoOverlap();
  virtual Mistral_Expression* add(MistralSolver *solver, bool top_level);
};

class Mistral_Precedence: public Mistral_binop
{
public:
  Mistral_Precedence(Mistral_Expression *var1, Mistral_Expression *var2, int constant);
  virtual ~Mistral_Precedence();
  virtual Mistral_Expression* add(MistralSolver *solver, bool top_level);
};

class Mistral_le: public Mistral_binop
{
public:
  Mistral_le(Mistral_Expression *var1, Mistral_Expression *var2);
  Mistral_le(Mistral_Expression *var1, int int_arg);
  virtual ~Mistral_le();
  virtual Mistral_Expression* add(MistralSolver *solver, bool top_level);
};

class Mistral_ge: public Mistral_binop
{
public:
  Mistral_ge(Mistral_Expression *var1, Mistral_Expression *var2);
  Mistral_ge(Mistral_Expression *var1, int int_arg);
  virtual ~Mistral_ge();
  virtual Mistral_Expression* add(MistralSolver *solver, bool top_level);
};

class Mistral_lt: public Mistral_binop
{
public:
  Mistral_lt(Mistral_Expression *var1, Mistral_Expression *var2);
  Mistral_lt(Mistral_Expression *var1, int int_arg);
  virtual ~Mistral_lt();
  virtual Mistral_Expression* add(MistralSolver *solver, bool top_level);
};

class Mistral_gt: public Mistral_binop
{
public:
  Mistral_gt(Mistral_Expression *var1, Mistral_Expression *var2);
  Mistral_gt(Mistral_Expression *var1, int int_arg);
  virtual ~Mistral_gt();
  virtual Mistral_Expression* add(MistralSolver *solver, bool top_level);
};


class Mistral_Minimise : public Mistral_Expression
{
private:
  Mistral_Expression* _exp;
  BuildObjectObjective *_obj;
public:
  Mistral_Minimise(Mistral_Expression *var);
  virtual ~Mistral_Minimise();
  virtual Mistral_Expression* add(MistralSolver *solver, bool top_level);
};


class Mistral_Maximise : public Mistral_Expression
{
private:
  Mistral_Expression* _exp;
  BuildObjectObjective *_obj;
public:
  Mistral_Maximise(Mistral_Expression *var);
  virtual ~Mistral_Maximise();
  virtual Mistral_Expression* add(MistralSolver *solver, bool top_level);
};


/**
   The solver itself
*/
class MistralSolver
{

public:

  Mistral::Solver *solver;
  Mistral::CSP *model;

  //MistralIntArray communicate;
  bool feature_ready;
  bool is_copy;
  //CSPXMLParser::MistralCallback *cb; 
  //CSPXMLParser::XMLParser_libxml2<CSPXMLParser::MistralCallback> *parser;

  Mistral::ConstraintNogoodBase* nogood_base;

  Vector< VarValuation > valuation;
  Vector< int > valuation_size;
  Vector< Mistral_Expression* > decisions;

  std::vector< std::vector<int> > graph;

  int first_decision_level;
  int saved_level;

  MistralSolver();
  virtual ~MistralSolver();

  // add an expression, in the case of a tree of expressions,
  // each node of the tree is added separately, depth first.
  void add(Mistral_Expression* arg);

  Mistral_Expression* get_expression(const int i);
  int num_expression();
  int max_expression_id();
  

  // used to initialise search on a given subset of variables
  void initialise(MistralExpArray& arg);
  // initialise the solver before solving (no more calls to add after this)
  void initialise();

  // solving methods
  int solve();
  int solveAndRestart(const int policy = GEOMETRIC, 
		      const unsigned int base = 32, 
		      const double factor = 1.3333333,
		      const double decay = 0.0,
		      const int reinit = -1);
  int startNewSearch();
  int getNextSolution();
  int sacPreprocess(const int type);

  int get_level(); // { return solver->level; }
  int get_decision_id();

  bool propagate();
  //void branch_on(const char* op, Mistral_Expression* x, const int v);
  void save();
  void post(const char* op, Mistral_Expression* x, int v);
  bool undo(const int nlevel);
  void deduce(const char* op, Mistral_Expression* x, int v);
  void deduce();
  bool branch_right();
  void store_solution();

  // parameter tuning methods
  void setHeuristic(const char* var_heuristic, const char* val_heuristic, const int rand);
  void setAntiLex(MistralExpArray& vars);
  void setFailureLimit(const int cutoff);  
  void setNodeLimit(const int cutoff);  
  void setTimeLimit(const int cutoff);
  void setVerbosity(const int degree);
  void setRandomized(const int degree);
  void setRandomSeed(const int seed);
  void forceFiniteDomain(MistralExpArray& vars);
  void addNogood(MistralExpArray& vars, 
		 MistralIntArray& vals);
  void guide(MistralExpArray& vars, 
	     MistralIntArray& vals,
	     MistralDoubleArray& probs);
  void backtrackTo(const int level);
  void upOneLevel();
  void presolve();
  void assign(Mistral_Expression *X, const int v);
  void increase_init_level(const int i);
  void decrease_init_level(const int i);
  void reset(bool full);
  void setLowerBounds(MistralExpArray& vars, 
		      MistralIntArray& vals);
  void setUpperBounds(MistralExpArray& vars, 
		      MistralIntArray& vals);
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

  void test_x60();
  void printPython();

  void load_xml(const char* filename, const int type=4);
  const char* get_feature_name(int i);
  double get_feature(int i);
  void get_features(MistralDoubleArray& f);
  void get_features();
  int num_vars();
  int get_degree(int i);
  void extract_graph();
  int numNodes();
  int degree(const int x);
  int get_neighbor(const int x, const int y);
};



#endif // _PYTHON_H
