/* -*- mode: C++; c-basic-offset: 2; indent-tabs-mode: nil -*- */
/*
 *  Main authors:
 *     Guido Tack <tack@gecode.org>
 *
 *  Copyright:
 *     Guido Tack, 2007
 *
 *  Last modified:
 *     $Date: 2010-07-21 11:42:47 +0200 (Wed, 21 Jul 2010) $ by $Author: tack $
 *     $Revision: 11243 $
 *
 *  This file is part of Gecode, the generic constraint
 *  development environment:
 *     http://www.gecode.org
 *
 *  Permission is hereby granted, free of charge, to any person obtaining
 *  a copy of this software and associated documentation files (the
 *  "Software"), to deal in the Software without restriction, including
 *  without limitation the rights to use, copy, modify, merge, publish,
 *  distribute, sublicense, and/or sell copies of the Software, and to
 *  permit persons to whom the Software is furnished to do so, subject to
 *  the following conditions:
 *
 *  The above copyright notice and this permission notice shall be
 *  included in all copies or substantial portions of the Software.
 *
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 *  EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 *  MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 *  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
 *  LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 *  OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 *  WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 */

#ifndef __GECODE_FLATZINC_HH__
#define __GECODE_FLATZINC_HH__

#include <iostream>
#include <string>
#include <map>
#include <set>

#include "conexpr.hpp"
#include "ast.hpp"
#include "varspec.hpp"

#include <mistral_solver.hpp>
#include <mistral_search.hpp>

//#define _VERBOSE_PARSER 100
//#define _VERIFICATION
//#define _DEBUG_VERIFICATION
//#define _FLATZINC_OUTPUT
using namespace Mistral;

typedef Vector<Variable> IntVarArray;
typedef Vector<Variable> BoolVarArray;
typedef Vector<Variable> SetVarArray;


/**
 * \namespace Gecode::FlatZinc
 * \brief Interpreter for the %FlatZinc language
 *
 * The Gecode::FlatZinc namespace contains all functionality required
 * to parse and solve constraint models written in the %FlatZinc language.
 *
 */

namespace FlatZinc {


/**
 * \brief Output support class for %FlatZinc interpreter
 *
 */
class Printer {
private:
	AST::Array* _output;
	void printElem(std::ostream& out,
			Solver& solver,
			AST::Node* ai,
			const IntVarArray& iv,
			const BoolVarArray& bv,
			const SetVarArray& sv
	) const;

public:
	Printer(void) : _output(NULL) {}
	void init(AST::Array* output);

	void print(std::ostream& out,
			Solver& solver,
			const IntVarArray& iv,
			const BoolVarArray& bv,
			const SetVarArray& sv
	) const;

	~Printer(void);
private:
	Printer(const Printer&);
	Printer& operator=(const Printer&);
};


#ifdef _VERIFICATION
//we need this only for verification
class SolutionValue
{
public:

	std::string get_string()
	{
		std::string tmp;
		//	std::cout << " \n type " << type << std::endl;
		if (type <= 4)
		{
			//Int, Bool, String, Atom
			tmp = val;
		}
		else
		{
			switch(type)
			{
			case 5:
			{
				if (ai->isSet())
				{
					//tmp = " isSet yess";
					std::ostringstream  oss;
					oss.str("");
					AST::SetLit* s = ai->getSet();
					if (s->interval) {
						oss << s->min << ".." << s->max;
					} else {
						oss << "{";
						for (unsigned int i=0; i<s->s.size(); i++) {
							oss << s->s[i] << (i < s->s.size()-1 ? ", " : "}");
						}
					}
					tmp = oss.str();
				}
			}break;
			case 6:
			{
				if (var != NULL)
				{
					//					tmp = " boolvar yes";
					std::ostringstream  oss;
					//oss.str("");
					int lb = var->get_solution_min();
					int ub = var->get_solution_max();
					if (lb == 1)
					{
						oss << "true";
					}
					else
						if (ub == 0)
						{
							oss << "false";
						} else
						{
							//					oss << "false..true";
							oss << "true";
						}
					tmp = oss.str();
				}
			}break;
			case 7:
			{
				if (var != NULL)
				{
					//	tmp = " intvar yes ";
					std::ostringstream  oss;
					oss.str("");
					//	var->display(std::cout);
					//		var->get_solution_str_value().
					//						std::cout << "\n the domain is \n " << << std::endl;
					//						std::cout << "\n and the value is \n " << var->get_solution_int_value() << std::endl;
					oss << var->get_solution_int_value();
					tmp =oss.str();
				}
			}break;
			case 8:
			{
				if (var != NULL)
				{
					//	tmp = " isSetVar yes";
					std::ostringstream  oss;
					SetExpression *x = (SetExpression*)(var->expression);
					std::set<int> lb;
					std::set<int> ub;
					for(unsigned int i=0; i<x->elts_ub.size; ++i)
					{
						if(x->get_index_var(i).get_solution_min())
						{
							lb.insert(x->elts_ub[i]);
							ub.insert(x->elts_ub[i]);
						}
						else if(x->children[i].get_solution_max())
							ub.insert(x->elts_ub[i]);
					}
					oss << "{";
					for( std::set<int>::const_iterator i = ub.begin(); i != ub.end(); ++i)
					{
						if( i != ub.begin() ) oss << ", ";
						oss << *i;
					}
					oss << "}";
					tmp =oss.str();
				}
			}break;
			case 9:
			{
				//	tmp = "array ";
				int size =a.size();
				std::ostringstream  oss;
				oss << "[";
				for (int i = 0; i < size; i++)
				{
					oss  << a[i].get_string();
					if (i< (size -1))
						oss << ", ";
				}
				oss << "]";
				tmp= oss.str();
			}break;
			}
		}
		return tmp;
	}
	void set_var(Variable * _var){var=_var;}
	void set_type(unsigned int _type){type= _type;}
	void set_val (std::string _val) {val= _val;}
	void set_ai (AST::Node* _ai){ai =_ai;}
	void set_a (std::vector<SolutionValue> _a){a = _a;}

private:
	Variable * var ;
	unsigned int type;
	std::string val;
	AST::Node* ai;
	std::vector<SolutionValue> a;
};
#endif

/**
 * \brief A space that can be initialized with a %FlatZinc model
 *
 */
class FlatZincModel {
public:
	enum Meth {
		SATISFACTION, //< Solve as satisfaction problem
		MINIMIZATION, //< Solve as minimization problem
		MAXIMIZATION  //< Solve as maximization problem
  	};
protected:
	/// Mistral stuff
	Solver &solver;

  // search stuff
  Vector< Vector< Variable > >  fz_search_sequences;
  Vector< BranchingHeuristic *> fz_search_heuristics;
  Vector< RestartPolicy * >     fz_search_policies;
  Vector< Goal * >              fz_search_goals;


	/// Options
	BranchingHeuristic *_option_heuristic;
  std::string _variable_ordering;
  std::string _value_ordering;
  int _randomization;
  //std::string _restart_policy;

	RestartPolicy *_option_policy;
	bool _option_rewriting;
	bool _option_simple_rewriting;
  bool _option_enumerate;
  bool _option_display_mistral_model;
  bool _option_display_solution;
  bool _option_annotations;
  int _option_parity;
	////


	/// Number of integer variables
	int intVarCount;
	/// Number of Boolean variables
	int boolVarCount;
	/// Number of set variables
	int setVarCount;

	/// Index of the integer variable to optimize
	int _optVar;

	/// Whether to solve as satisfaction or optimization problem
	Meth _method;

	/// Annotations on the solve item
	AST::Array* _solveAnnotations;
	struct clause_struct
	{
		//Variables appearing positively in the clause
		Vector<Variable > pos;
		//Variables appearing negatively in the clause
		Vector<Variable > neg;
	};

	//structure used for encoding clauses
	Vector<clause_struct> _clauses;
	Vector< Vector< Literal > > cnf;

public:


	/// The integer variables
	IntVarArray iv;
	/// Indicates whether an integer variable is introduced by mzn2fzn
	std::vector<bool> iv_introduced;
	/// Indicates whether an integer variable aliases a Boolean variable
	std::vector<int> iv_boolalias;
	/// The Boolean variables
	BoolVarArray bv;
	/// Indicates whether a Boolean variable is introduced by mzn2fzn
	std::vector<bool> bv_introduced;
	/// The Set variables
	SetVarArray sv;
	/// Indicates whether a set variable is introduced by mzn2fzn
	std::vector<bool> sv_introduced;

	/// vars fixed to true and false, in case they are encountered often
	//Variable vartrue(1,1);
	//Variable varfalse(0,0);
	//std::map<int, Variable> constants;

	/// Construct empty space
	FlatZincModel(Solver& s);

	/// Destructor
	~FlatZincModel(void);

	/// Initialize space with given number of variables
	void init(int intVars, int boolVars, int setVars);

	/// Create new integer variable from specification
	void newIntVar(IntVarSpec* vs);
	/// Link integer variable \a iv to Boolean variable \a bv
	void aliasBool2Int(int iv, int bv);
	/// Return linked Boolean variable for integer variable \a iv
	int aliasBool2Int(int iv);
	/// Create new Boolean variable from specification
	void newBoolVar(BoolVarSpec* vs);
	/// Create new set variable from specification
	void newSetVar(SetVarSpec* vs);

	/// Post a constraint specified by \a ce
	void postConstraint(const ConExpr& ce, AST::Node* annotation);

	/// Post the solve item
	void solve(AST::Array* annotation);
	/// Post that integer variable \a var should be minimized
	void minimize(int var, AST::Array* annotation);
	/// Post that integer variable \a var should be maximized
	void maximize(int var, AST::Array* annotation);

	/// setup parameters from the command line
	void set_parameters(SolverParameters& p);

	/// setup parameters from the command line
  void set_strategy(std::string var_o, std::string val_o, const int rand, std::string r_pol);

	/// setup the rewriting step
	void set_rewriting(const bool on);
	/// setup a rewriting step (mainly to eliminate bool2int)
	void set_simple_rewriting(const bool on);

	/// setup the rewriting step
	void set_parity_processing(const int lvl);

	/// setup the rewriting step
  void set_enumeration(const bool on);

	/// setup the rewriting step
  void set_display_model(const bool on);

	/// setup the rewriting step
  void set_display_solution(const bool on);

	/// setup annotations
  void set_annotations(const bool on);


  /// get annotations from the flatzinc model
  void get_annotations();

	/// Run the search
	void run(std::ostream& out, const Printer& p);

  /// Produce output on \a out using \a p
  void print_final(std::ostream& out, const Printer& p) const;
  void print_solution(std::ostream& out, const Printer& p) const;

	/**
	 * \brief Remove all variables not needed for output
	 *
	 * After calling this function, no new constraints can be posted through
	 * FlatZinc variable references, and the createBranchers method must
	 * not be called again.
	 *
	 */
	void shrinkArrays(Printer& p);

	/// Return whether to solve a satisfaction or optimization problem
	Meth method(void) const;

	/// Return index of variable used for optimization
	int optVar(void) const;

	/**
	 * \brief Create branchers corresponding to the solve item annotations
	 *
	 * If \a ignoreUnknown is true, unknown solve item annotations will be
	 * ignored, otherwise a warning is written to \a err.
	 */
	void createBranchers(AST::Node* ann, bool ignoreUnknown,
			std::ostream& err = std::cerr);

	/// Return the solve item annotations
	AST::Array* solveAnnotations(void) const;

	bool getAnnotations(AST::Call* , Vector<Variable> &, std::string& , std::string& );

	/// Implement optimization
	//void constrain();

	/// options
	bool findall; // find all solutions
	bool finished()
	{
		Mistral::Outcome outcome = solver.statistics.outcome;
		//std::cout << outcome;
		return ((outcome == SAT) || (outcome == OPT) || ((outcome == LIMITOUT) && (solver.statistics.num_solutions >0)));
	}

	void add_clause(Vector<Variable> pos ,  Vector<Variable> neg);
	void encode_clause(Vector<Variable> pos ,  Vector<Variable> neg);
	void encode_clauses();


#ifdef _VERIFICATION
	//Verification
	std::vector<std::pair<std::string, std::vector<SolutionValue > > > verif_constraints;
	std::ostringstream oss;
	SolutionValue node2SolutionValue(AST::Node * ai );

#endif
};

/// %Exception class for %FlatZinc errors
class Error {
private:
	const std::string msg;
public:
	Error(const std::string& where, const std::string& what)
	: msg(where+": "+what) {}
	const std::string& toString(void) const { return msg; }
};

/**
 * \brief Parse FlatZinc file \a fileName into \a fzs and return it.
 *
 * Creates a new empty FlatZincSpace if \a fzs is NULL.
 */
FlatZincModel* parse(const std::string& fileName,
		Solver& solver,
		Printer& p, std::ostream& err = std::cerr,
		FlatZincModel* fzs=NULL);

/**
 * \brief Parse FlatZinc from \a is into \a fzs and return it.
 *
 * Creates a new empty FlatZincSpace if \a fzs is NULL.
 */
FlatZincModel* parse(std::istream& is,
		Solver& solver,
		Printer& p, std::ostream& err = std::cerr,
		FlatZincModel* fzs=NULL);


  /*! \class SolutionListener
    \brief SolutionListener Class

    * Called whenever the solver solutions *
    
    This is used to implement procedures triggered by solutions
  */
  class SolutionPrinter : public SolutionListener {

  public:

    Printer *p_;
    FlatZincModel *fm_;
    Mistral::Solver *solver_;

    // SolutionPrinter(Printer *p, FlatZincModel *fm, Mistral::Solver *s) : p_(p), fm_(fm), solver_(s) {
    //   solver->add((SolutionPrinter*)this);
    // }

    SolutionPrinter(Printer *p, FlatZincModel *fm, Mistral::Solver *s);
    virtual ~SolutionPrinter();    

    virtual void notify_solution() ;
    //{
    //   fm_->print(std::cout, *p_);
    // };
    
    std::ostream& display(std::ostream& os) { os << "Flatzinc solution-printer"; return os; }    
  };

}

#endif

// STATISTICS: flatzinc-any
