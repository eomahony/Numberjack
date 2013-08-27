/* -*- mode: C++; c-basic-offset: 2; indent-tabs-mode: nil -*- */
/*
 *  Main authors:
 *     Guido Tack <tack@gecode.org>
 *
 *  Copyright:
 *     Guido Tack, 2007
 *
 *  Last modified:
 *     $Date: 2010-05-11 12:33:38 +0200 (Tue, 11 May 2010) $ by $Author: tack $
 *     $Revision: 10940 $
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

#include "flatzinc.hpp"
#include "registry.hpp"
#include <iomanip>

#include <vector>
#include <assert.h>
#include <string>
#include <set>



#include <mistral_variable.hpp>

//#define _DEBUG_FLATZINC true
//#define _VERBOSE_PARSER 1000



using namespace std;

namespace FlatZinc {
inline
set<int> setrange(int min, int max)
{
	set<int> rv;
	for(int i = min; i <= max; ++i)
		rv.insert(i);
	return rv;
}

set<int> vs2is(IntVarSpec* vs) {
	if (vs->assigned) {
		return setrange(vs->i,vs->i);
	}
	if (vs->domain()) {
		AST::SetLit* sl = vs->domain.some();
		if (sl->interval) {
			return setrange(sl->min, sl->max);
		} else {
			set<int> rv;
			for (int i=sl->s.size(); i--;)
				rv.insert(sl->s[i]);
			return rv;
		}
	}
	return setrange(-1000, 1000);
}


Variable vs2var(IntVarSpec* vs) {
	Variable x;
	if (vs->assigned) {
		Variable y(vs->i);
		x = y;
	} else if (vs->domain()) {
		AST::SetLit* sl = vs->domain.some();
		if (sl->interval) {
			Variable y(sl->min, sl->max);
			x = y;
		} else {
			set<int> __tmp_set;
			int size = sl->s.size();
			for (int i=0; i<size;i++)
				__tmp_set.insert(sl->s[i]);
			//	cout << "myset contains:";
			//	for (set<int>::iterator  it=__tmp_set.begin(); it!=__tmp_set.end(); it++)
			//	cout << " " << *it;
			//	cout << endl;
			Vector<int> values;
			for (set<int>::iterator it=__tmp_set.begin(); it!=__tmp_set.end(); it++)
				values.add(*it);

			Variable y(values);
			x = y;
		}
	} else {
		Variable y(-INFTY/1024, +INFTY/1024);
		x = y;
	}

	//std::cout << x << " in " << x.get_domain() << std::endl;

	return x;
}

int vs2bsl(BoolVarSpec* bs) {
	if (bs->assigned) {
		return bs->i;
	}
	if (bs->domain()) {
		AST::SetLit* sl = bs->domain.some();
		assert(sl->interval);
		return std::min(1, std::max(0, sl->min));
	}
	return 0;
}

int vs2bsh(BoolVarSpec* bs) {
	if (bs->assigned) {
		return bs->i;
	}
	if (bs->domain()) {
		AST::SetLit* sl = bs->domain.some();
		assert(sl->interval);
		return std::max(0, std::min(1, sl->max));
	}
	return 1;
}

int ann2ivarsel(AST::Node* ann) {
	if (/*AST::Atom* s =*/ dynamic_cast<AST::Atom*>(ann)) {
		// if (s->id == "input_order")
		//   return TieBreakVarBranch<IntVarBranch>(INT_VAR_NONE);
	}
	std::cerr << "Warning, ignored search annotation: ";
	ann->print(std::cerr);
	std::cerr << std::endl;
	return 0;
}

int ann2ivalsel(AST::Node* ann) {
	if (/*AST::Atom* s =*/ dynamic_cast<AST::Atom*>(ann)) {
		// if (s->id == "indomain_min")
		//   return INT_VAL_MIN;
	}
	std::cerr << "Warning, ignored search annotation: ";
	ann->print(std::cerr);
	std::cerr << std::endl;
	return 0;
}

int ann2asnivalsel(AST::Node* ann) {
	if (/*AST::Atom* s =*/ dynamic_cast<AST::Atom*>(ann)) {
		// if (s->id == "indomain_min")
		//   return INT_ASSIGN_MIN;
	}
	std::cerr << "Warning, ignored search annotation: ";
	ann->print(std::cerr);
	std::cerr << std::endl;
	return 0;
}


FlatZincModel::FlatZincModel(Solver &s)
: solver(s), 
  _option_heuristic(NULL), 
  _option_policy(NULL), 
  _option_rewriting(false),
  _option_simple_rewriting(false),
  _option_enumerate(false),
  _option_display_mistral_model(false),
  _option_annotations(false),
  _option_parity(0),
  intVarCount(-1), boolVarCount(-1), setVarCount(-1), _optVar(-1),
  _solveAnnotations(NULL)
{
}

void
FlatZincModel::init(int intVars, int boolVars, int setVars) {
	intVarCount = 0;
	iv = IntVarArray(intVars);
	iv_introduced = std::vector<bool>(intVars);
	iv_boolalias = std::vector<int>(intVars);
	boolVarCount = 0;
	bv = BoolVarArray(boolVars);
	bv_introduced = std::vector<bool>(boolVars);
	setVarCount = 0;
	sv = SetVarArray(setVars);
	sv_introduced = std::vector<bool>(setVars);
}

void
FlatZincModel::newIntVar(IntVarSpec* vs) {

#ifdef _VERBOSE_PARSER
	if((intVarCount % _VERBOSE_PARSER) == 0) {
		std::cout << "i";
		std::cout.flush();
	}
#endif

	if (vs->alias) {
		iv[intVarCount++] = iv[vs->i];
	} else {
		Variable x = vs2var(vs);
		iv[intVarCount++] = x;
	}
	iv_introduced[intVarCount-1] = vs->introduced;
	iv_boolalias[intVarCount-1] = -1;
}

void
FlatZincModel::newSetVar(SetVarSpec* vs) {

#ifdef _VERBOSE_PARSER
  if((setVarCount % _VERBOSE_PARSER) == 0) {
    std::cout << "s";
    std::cout.flush();
  }
#endif

  if (vs->alias) {
    sv[setVarCount++] = sv[vs->i];
  } else if( vs->assigned) {
    assert(vs->upperBound());
    AST::SetLit* vsv = vs->upperBound.some();
    if (vsv->interval) {

      std::cout << "create a set variable from a ground interval [" << vsv->min << ".." << vsv->max << "]\n";

      Variable x = SetVariable(vsv->min, vsv->max);
      sv[setVarCount++] = x;

      std::cout << x << " in " << x.get_domain() << std::endl;

    } else {
      if( vsv->s.empty() ) {

        // std::cout << "create an empty set variable\n";

        // Variable x = SetVariable();

        // std::cout << x << " in " << x.get_domain() << std::endl;

        // //x.exclude(solver, 0, NO_REASON);
        // //sv[setVarCount++] = x;
      } else {
        int umin = vsv->s[0], umax = vsv->s[0];
        for(size_t i = 1; i != vsv->s.size(); ++i) {
          umin = std::min(umin, vsv->s[i]);
          umax = std::max(umax, vsv->s[i]);
        }
        Variable x = SetVariable(umin, umax);
        sv[setVarCount++] = x;
        // for(size_t i = 0; i != vsv->s.size(); ++i)
        //   x.include(solver, vsv->s[i], NO_REASON);
        // for(int i = x.umin(solver), iend = x.umax(solver); i <= iend; ++i)
        //   if( !x.includes(solver, i) )
        //     x.exclude(solver, i, NO_REASON);
      }
    }
  } else if( vs->upperBound() ) {
    AST::SetLit* vsv = vs->upperBound.some();
    Variable x = SetVariable(vsv->min, vsv->max);
    //setvar x = solver.newSetVar(vsv->min, vsv->max);
    sv[setVarCount++] = x;
    // if( !vsv->interval ) {
    //   int prev = vsv->min;
    //   for(size_t i = 0; i != vsv->s.size(); ++i) {
    //     if( vsv->s[i] > prev+1 ) {
    //       for(int q = prev+1; q != vsv->s[i]; ++q)
    //         x.exclude(solver, q, NO_REASON);
    //     }
    //     prev = vsv->s[i];
    //   }
    // } // otherwise everything is unset and we are done here
  } else {
    // completely free
    //setvar x = solver.newSetVar(-1000, 1000);
    Variable x = SetVariable(-1000, 1000);
    sv[setVarCount++] = x;
  }
  sv_introduced[setVarCount-1] = vs->introduced;
}

void
FlatZincModel::aliasBool2Int(int iv, int bv) {
	iv_boolalias[iv] = bv;
}
int
FlatZincModel::aliasBool2Int(int iv) {
	return iv_boolalias[iv];
}

void
FlatZincModel::newBoolVar(BoolVarSpec* vs) {

#ifdef _VERBOSE_PARSER
	if((boolVarCount % _VERBOSE_PARSER) == 0) {
		std::cout << "b";
		std::cout.flush();
	}
#endif

	if (vs->alias) {
		bv[boolVarCount++] = bv[vs->i];
	} else {
		bv[boolVarCount++] = Variable(vs2bsl(vs), vs2bsh(vs));
	}
	bv_introduced[boolVarCount-1] = vs->introduced;
}


#ifdef _VERIFICATION
SolutionValue
FlatZincModel::node2SolutionValue(AST::Node * ai ) {
	SolutionValue v;
	int k;
	if (ai->isBool())
	{
		oss.str("");
		//		cout<<"isBool ";
		oss << (ai->getBool() ? "true" : "false");
		v.set_type(1);
		v.set_val(oss.str());
	}
	else
		if (ai->isString())
		{
			oss.str("");
			v.set_type(2);
			//	cout<<" Test if node is a isString node";
			std::string s = ai->getString();
			for (unsigned int i=0; i<s.size(); i++)
			{
				if (s[i] == '\\' && i<s.size()-1)
				{
					switch (s[i+1])
					{
					case 'n': oss << "\n"; break;
					case '\\': oss << "\\"; break;
					case 't': oss << "\t"; break;
					default: oss << "\\" << s[i+1];
					}
					i++;
				}
				else
				{
					oss << s[i];
				}
			}

			v.set_val(oss.str());
		}
		else
			if (ai->isInt(k))
			{
				oss.str("");
				oss << k;
				v.set_type(3);
				v.set_val(oss.str());
			}
			else
				if (ai->isAtom())
				{
					oss.str("");
					//cout<<" Test if node is isAtom node";
					v.set_type(4);
					v.set_val("ATOM is not yet supported");
				}
				else
					if (ai->isSet())
					{
						//		cout<<" Test if node is isSet  node";
						AST::SetLit* setLit = new AST::SetLit();
						setLit->interval =ai->getSet()->interval;
						setLit->min =ai->getSet()->min;
						setLit->max =ai->getSet()->max;
						setLit->s =ai->getSet()->s;
						v.set_ai(setLit);
						v.set_type(5);
					}
					else
						if(	ai->isBoolVar() )
						{
							//		cout<<"isBoolVar ";
							v.set_type(6);
							Variable *_var =  & (bv[ai->getBoolVar()]) ;
							v.set_var( _var);
						}
						else
							if (ai->isArray())
							{
								std::vector<SolutionValue> __tmp;
								//	cout<<" Test if node is isArray node" << endl ;
								for (unsigned int i =0; i< ai->getArray()->a.size(); i++)
								{
									AST::Node *n = ai->getArray()->a[i];
									__tmp.push_back(node2SolutionValue(n));
									v.set_a(__tmp);

								}
								v.set_type(9);
							}
							else
								if (ai->isSetVar())
								{
									//		cout<<" Test if node is a isSetVar node";
									Variable *_var =  & (sv[ai->getSetVar()]);
									v.set_var( _var);
									v.set_type(8);
								}
								else
									if (ai->isIntVar())
									{
										//	cout<<"isIntVar ";
										Variable *_var = &( iv[ai->getIntVar()]) ;
										//std::cout << _var.get_domain() << std::endl;
										v.set_var(_var);
										v.set_type(7);
									}
	return v;
}
#endif

void
FlatZincModel::postConstraint(const ConExpr& ce, AST::Node* ann) {
	try {
#ifdef _VERBOSE_PARSER
		if((solver.constraints.size % _VERBOSE_PARSER) == 0) {
			std::cout << "c";
			std::cout.flush();
		}
#endif
		registry().post(solver, *this, ce, ann);

#ifdef _VERIFICATION
		std::vector<SolutionValue > __vars;
		//	std::cout <<" \n posting " << ce.id <<std::endl ;
		for (unsigned int i=0; i<ce.args->a.size(); i++)
		{
			AST::Node * ai = ce.args->a[i];
			__vars.push_back(node2SolutionValue(ai));
		}
		pair<std::string, std::vector<SolutionValue > > pair (ce.id , __vars);
		verif_constraints.push_back(pair);
#endif

	} catch (AST::TypeError& e) {
		throw FlatZinc::Error("Type error", e.what());
	}

}





void flattenAnnotations(AST::Array* ann, std::vector<AST::Node*>& out) {
	for (unsigned int i=0; i<ann->a.size(); i++) {
		if (ann->a[i]->isCall("seq_search")) {
			AST::Call* c = ann->a[i]->getCall();
			if (c->args->isArray())
				flattenAnnotations(c->args->getArray(), out);
			else
				out.push_back(c->args);
		} else {
			out.push_back(ann->a[i]);
		}
	}
}


bool FlatZincModel::getAnnotations( AST::Call*  c , Vector<Variable> &__vars, std::string & _varHeuristic , std::string & _valHeuristic)
{

	if (c->isCall("int_search") || c->isCall("bool_search"))
	{
	if (c->args->isArray())
	{
	//	cout<< "c array size : "<< c->args->getArray()->a.size() ;
		if (c->args->getArray()->a[0]->isArray())
		{
			AST::Array * __varsArray = c->args->getArray()->a[0]->getArray();
			//cout<< " c number of branching variables : "<< __varsArray->a.size() << endl;

			for (unsigned int j=0; j< __varsArray->a.size(); j++ )
			{
                          if (__varsArray->a[j]->isIntVar())
                            __vars.push_back(iv[__varsArray->a[j]->getIntVar()].get_var());
                          else if (__varsArray->a[j]->isSetVar())
                            __vars.push_back( sv[__varsArray->a[j]->getSetVar()]);
                          else if (__varsArray->a[j]->isBoolVar())
                            __vars.push_back( bv[__varsArray->a[j]->getBoolVar()].get_var());
                          else
                            return false;
			}
		}

		else return false;

		//test if atom ?
		std::ostringstream __tmpStream;
		c->args->getArray()->a[1]->print(__tmpStream);
		_varHeuristic = __tmpStream.str();
		__tmpStream.str("");
		c->args->getArray()->a[2]->print(__tmpStream);
		_valHeuristic = __tmpStream.str();

	}
	else
		return false;
	}
	else
		return false;



	return true;
}


void
FlatZincModel::createBranchers(AST::Node* ann, bool ignoreUnknown,
		std::ostream& err) {
	if (ann) {
		err << "Warning, ignored search annotation: ";
		ann->print(err);
		err << std::endl;
	}
}

AST::Array*
FlatZincModel::solveAnnotations(void) const {
	return _solveAnnotations;
}

void
FlatZincModel::solve(AST::Array* ann) {
	_method = SATISFACTION;
	_solveAnnotations = ann;
}

// void
// FlatZincModel::enumerate(AST::Array* ann) {
// 	_method = ENUMERATION;
// 	_solveAnnotations = ann;
// 	// Branch on optimization variable to ensure that it is given a value.
// 	/*
// 	AST::Array* args = new AST::Array(4);
// 	args->a[0] = new AST::Array(new AST::IntVar(_optVar));
// 	args->a[1] = new AST::Atom("input_order");
// 	args->a[2] = new AST::Atom("indomain_min");
// 	args->a[3] = new AST::Atom("complete");
// 	AST::Call* c = new AST::Call("int_search", args);
// 	if (!ann)
// 		ann = new AST::Array(c);
// 	else
// 		ann->a.push_back(c);
// 	*/

// }

void
FlatZincModel::minimize(int var, AST::Array* ann) {
	_method = MINIMIZATION;
	_optVar = var;
	_solveAnnotations = ann;
	// Branch on optimization variable to ensure that it is given a value.
	/*
	AST::Array* args = new AST::Array(4);
	args->a[0] = new AST::Array(new AST::IntVar(_optVar));
	args->a[1] = new AST::Atom("input_order");
	args->a[2] = new AST::Atom("indomain_min");
	args->a[3] = new AST::Atom("complete");
	AST::Call* c = new AST::Call("int_search", args);
	if (!ann)
		ann = new AST::Array(c);
	else
		ann->a.push_back(c);
	*/

}

void
FlatZincModel::maximize(int var, AST::Array* ann) {
	_method = MAXIMIZATION;
	_optVar = var;
	_solveAnnotations = ann;
	// Branch on optimization variable to ensure that it is given a value.
	/*
	AST::Array* args = new AST::Array(4);
	args->a[0] = new AST::Array(new AST::IntVar(_optVar));
	args->a[1] = new AST::Atom("input_order");
	args->a[2] = new AST::Atom("indomain_min");
	args->a[3] = new AST::Atom("complete");
	AST::Call* c = new AST::Call("int_search", args);
	if (!ann)
		ann = new AST::Array(c);
	else
		ann->a.push_back(c);
	 */
}

FlatZincModel::~FlatZincModel(void) {
	delete _solveAnnotations;

        // for(unsigned int i=0; i<iv.size; ++i) {
        
        //   int domain_type = iv[i].domain_type;
        //   if     (domain_type ==  BITSET_VAR) delete iv[i].bitset_domain;
        //   else if(domain_type ==    LIST_VAR) delete iv[i].list_domain;
        //   else if(domain_type ==   RANGE_VAR) delete iv[i].range_domain;
        //   else if(domain_type == VIRTUAL_VAR) delete iv[i].virtual_domain;
        //   else if(domain_type ==  EXPRESSION) delete iv[i].expression;
        //   else if(domain_type !=   CONST_VAR) delete iv[i].variable;
 
        // }
}



void
FlatZincModel::set_parameters(SolverParameters& p) {

}


void
FlatZincModel::set_strategy(string var_o, string val_o, string r_pol) {
  _variable_ordering = var_o;
  _value_ordering = val_o;
	// _option_heuristic = solver.heuristic_factory(var_o, val_o);
	 _option_policy = solver.restart_factory(r_pol);
}

void
FlatZincModel::set_rewriting(const bool on) {
	_option_rewriting = on;
}

void
FlatZincModel::set_simple_rewriting(const bool on) {
	_option_simple_rewriting = on;
}

void
FlatZincModel::set_parity_processing(const int lvl) {
	_option_parity = lvl;
}

void
FlatZincModel::set_display_model(const bool on) {
	_option_display_mistral_model = on;
}

void
FlatZincModel::set_display_solution(const bool on) {
	_option_display_solution = on;
}

void 
FlatZincModel::set_enumeration(const bool on) {
  _option_enumerate = on;
}

void 
FlatZincModel::set_annotations(const bool on) {
  _option_annotations = on;
}


  void 
  FlatZincModel::get_annotations() {
    //cout << _solveAnnotations << endl;

    if (_solveAnnotations!= NULL)
      {
        //int __search_strategies = 1;
        //cout << " c annotations size : " << _solveAnnotations->a.size() << endl;
        //fz_search_goals.add(goal);

        if (_solveAnnotations->a[0]->isCall("seq_search")) {
          //cout << " c SEQ_SEARCH " << endl ;
          
          AST::Call* c = _solveAnnotations->a[0]->getCall();
          if (c->args->isArray())
            {
              int __number_phases =c->args->getArray()->a.size();
              //cout << " c Total number of search strategies:" << __number_phases  << endl;

              for (int j = 0; j< __number_phases; j++)
                {
                  Vector<Variable> __branching_variables;
                  std::string __var_heuristic;
                  std::string __val_heuristic;

                  if (getAnnotations(c->args->getArray()->a[j]->getCall(),__branching_variables, __var_heuristic ,__val_heuristic))
                    {
                      // cout << " c Search strategy number "<< j+1 << endl;
                      // cout << " c var heuristic : "<< __var_heuristic << endl;
                      // cout << " c val heuristic : "<< __val_heuristic << endl;
                      // cout << " c _variables size : "<< __branching_variables.size  << endl;
                      // cout << " c branching variables : \n "<< __branching_variables  << endl;

                      fz_search_sequences.add(__branching_variables);
                      fz_search_heuristics.add(solver.heuristic_factory(__var_heuristic, __val_heuristic));
                      fz_search_policies.add(NULL);
                      if(j) fz_search_goals.add(NULL);

                    }
                  // else
                  //   cout << " c Something wrong with search annotations. The solver will use the default search strategy." << endl;
                }
              //cout << " c seq_search is not yet supported." << endl;
            }
          // else
          //   cout << " c Something wrong with search annotations. The solver will use the default search strategy." << endl;
        }
        else
          {
            Vector<Variable> __branching_variables;
            std::string __var_heuristic;
            std::string __val_heuristic;

            //cout << " c Total number of search strategies:" << __search_strategies  << endl;
            if (getAnnotations(_solveAnnotations->getArray()->a[0]->getCall(), __branching_variables, __var_heuristic ,__val_heuristic ))
              {
                // cout << " c var heuristic : "<< __var_heuristic  << endl;
                // cout << " c val heuristic : "<< __val_heuristic  << endl;
                // cout << " c number of banching variables : "<< __branching_variables.size  << endl;
                // //cout << " c branching variables : \n "<< banching_variables  << endl;

                fz_search_sequences.add(__branching_variables);
                fz_search_heuristics.add(solver.heuristic_factory(__var_heuristic, __val_heuristic));
                fz_search_policies.add(NULL);
                
              }
            // else
            //   cout << " c Something wrong with search annotations. The solver will use the default search strategy." << endl;
          }

        // cout << sequences << std::endl;
        // for(unsigned int i=0; i<heuristics.size; ++i) {
        //   heuristics[i]->display(std::cout);
        //   std::cout << std::endl;
        // }
      }
    // else
    //   cout << " c No specific annotation. The solver will use the default search strategy." << endl;




    //cout << solver.sequence.size << endl;

    //exit(1);
  }


  void
  FlatZincModel::run(std::ostream& out, const Printer& p) {
    using std::setw;
    using std::setfill;
    
#ifdef _DEBUG_FLATZINC
    std::cout << " " << solver.parameters.prefix_comment << " run!" << std::endl;
#endif

    if(_option_rewriting) {
#ifdef _DEBUG_FLATZINC
      std::cout << "before rewriting:\n" << solver << std::endl;
#endif

      solver.rewrite() ;

#ifdef _DEBUG_FLATZINC
      std::cout << "after rewriting:\n" << solver << std::endl;
#endif
    }
    else
        if(_option_simple_rewriting) {
    #ifdef _DEBUG_FLATZINC
          std::cout << "before simple_rewriting:\n" << solver << std::endl;
    #endif

          solver.simple_rewrite() ;

    #ifdef _DEBUG_FLATZINC
          std::cout << "after simple_rewriting:\n" << solver << std::endl;
    #endif
        }


    solver.consolidate();
    solver.sequence.clear();


    if(_option_display_mistral_model)
      std::cout << " " << solver.parameters.prefix_comment << " mistral representation:\n " << solver << std::endl;
    
        
    Outcome result = UNKNOWN;


    Goal *goal = NULL;

    switch (_method) {
    case MINIMIZATION: {

      std::cout << " " << solver.parameters.prefix_comment << " Minimize " << iv[_optVar].get_var() << std::endl;
      
      goal = new Goal(Goal::MINIMIZATION, iv[_optVar].get_var());
      break;
    }
    case MAXIMIZATION: {
      std::cout << " " << solver.parameters.prefix_comment << " Maximize " << iv[_optVar].get_var() << std::endl;
      
      goal = new Goal(Goal::MAXIMIZATION, iv[_optVar].get_var());
      break;
    }
    case SATISFACTION: {
      std::cout << " " << solver.parameters.prefix_comment << " Solve " << std::endl;
      
      if(_option_enumerate) 
        goal = new Goal(Goal::ENUMERATION);
      else 
        goal = new Goal(Goal::SATISFACTION);
      break;
    }
    }


    solver.objective = goal;


    if(solver.is_pseudo_boolean())
      solver.set_learning_on();


    if(_option_annotations) {
      fz_search_goals.add(goal);
      get_annotations();
    }


    if(_option_parity)
      solver.parity_processing(_option_parity);


    _option_heuristic = solver.heuristic_factory(_variable_ordering, _value_ordering);
    

      
    if(fz_search_sequences.size < 2) {

      // there is no annotation, we use the default strategy
      if(fz_search_sequences.empty()) {
        //std::cout << solver.variables << std::endl;
        result = solver.depth_first_search(solver.variables, _option_heuristic, _option_policy, goal);
      } else {
        //std::cout << fz_search_sequences[0] << std::endl;
        //result = solver.depth_first_search(fz_search_sequences[0], _option_heuristic, _option_policy, goal);
        solver.initialise_search(fz_search_sequences[0], _option_heuristic, _option_policy, goal);
        //solver.heuristic->display(std::cout);
        result = solver.depth_first_search();
      }
    } else {
      // follows flatzinc model's annotations

      cout << " " << solver.parameters.prefix_comment 
           << " sequence search on " << fz_search_sequences << std::endl;

      Variable obj = iv[_optVar].get_var();


#ifdef _MONITOR
      for(unsigned int k=0; k<fz_search_sequences.size; ++k) {
        solver.monitor_list << "[ " ;
        for(unsigned int i=0; i<fz_search_sequences[k].size; ++i) {
          solver.monitor_list << fz_search_sequences[k][i] ;
          solver.monitor_list << " " ;
        }
        solver.monitor_list << "]\n" ;
      }
      solver.monitor_list << "objective: " << obj ;
      solver.monitor_list << "\n" ;
#endif

      result = solver.sequence_search(fz_search_sequences, 
                                      fz_search_heuristics, 
                                      fz_search_policies, 
                                      fz_search_goals);


      //cout << outcome2str(result) << " " << outcome2str(solver.statistics.outcome) << endl;


    }
      






#ifdef _DEBUG_VERIFICATION
    if ((result == SAT) || (result == OPT) || ((result == LIMITOUT) && (solver.statistics.num_solutions >0)))
      {
        cout<< " \n \n \n solver variables :  " << solver.variables.size ;
        int __size = solver.variables.size;
        for ( int i = 0; i < __size ; i++)
          {
            cout<< "  " << endl;
            std::cout << solver.variables[i].get_solution_int_value() << " (" << solver.variables[i].id() << ") ";
            std::cout << solver.variables[i] << " in " << solver.variables[i].get_domain() << " ";
          }
      }
#endif
 

  }


FlatZincModel::Meth
FlatZincModel::method(void) const {
	return _method;
}

int
FlatZincModel::optVar(void) const {
	return _optVar;
}

void
FlatZincModel::print_final(std::ostream& out, const Printer& p) const {

  //#ifdef _FLATZINC_OUTPUT
	Mistral::Outcome outcome = solver.statistics.outcome;
        if (outcome == OPT)
          out<<"==========";
        else if (solver.statistics.num_solutions && solver.objective && solver.objective->is_optimization())
          out<<"=====UNBOUNDED=====";
	else if (outcome == UNSAT)
          out<<"=====UNSATISFIABLE=====";
        else if (outcome != SAT)
          out<<"=====UNKNOWN=====";
	out<< std::endl;
        //#else
        //solver.statistics.print_full(out);
        //#endif


}


void
FlatZincModel::print_solution(std::ostream& out, const Printer& p) const {

  //#ifdef _FLATZINC_OUTPUT

  // for(int i=0; i<solver.variables.size; ++i) {
  //   std::cout << solver.variables[i] << ": " ;
  //   if(solver.variables[i].get_solution_min() == solver.variables[i].get_solution_max())
  //     std::cout << solver.variables[i].get_solution_min() << std::endl;
  //   else std::cout << solver.variables[i].get_solution_min() << ".." << solver.variables[i].get_solution_max() << std::endl;
  // }

  if(_option_display_solution) {
    if(solver.statistics.num_solutions) {
      p.print(out, solver, iv, bv, sv);
      
      if(_optVar >= 0)
        out << " " << solver.parameters.prefix_comment << " objective: " << iv[_optVar].get_var() 
            << " in [" << iv[_optVar].get_solution_min() 
            << ".." << iv[_optVar].get_solution_max() << "]" << endl;
    }
    out << "----------" << std::endl;
  }
  //#endif


}



void
Printer::init(AST::Array* output) {
	_output = output;
}

void
Printer::printElem(std::ostream& out,
		Solver& solver,
		AST::Node* ai,
		const IntVarArray& iv,
		const BoolVarArray& bv,
		const SetVarArray& sv
) const {
	int k;
	if (ai->isInt(k)) {
		out << k;
	} else if (ai->isIntVar()) {
          // int lb = iv[ai->getIntVar()].get_solution_min();
          // int ub = iv[ai->getIntVar()].get_solution_max();
          // if( lb == ub )
          //   out << lb;
          // else
          //   out << lb << ".." << ub;
          
          out << iv[ai->getIntVar()].get_solution_int_value();
	} else if (ai->isBoolVar()) {
		/*
		int lb = bv[ai->getBoolVar()].get_solution_min();
		int ub = bv[ai->getBoolVar()].get_solution_max();
		if (lb == 1) {
			out << "true";
		} else if (ub == 0) {
			out << "false";
		} else {
			out << "false..true";
		}
		*/
		int lb = bv[ai->getBoolVar()].get_solution_int_value();
		if (lb == 1)
		{
			out << "true";
		}
		else
			out << "false";
	} else if( ai->isSetVar()) {
		SetExpression *x = (SetExpression*)(sv[ai->getSetVar()].expression);
		set<int> lb;
		set<int> ub;
		for(unsigned int i=0; i<x->elts_lb.size; ++i) {
                  lb.insert(x->elts_lb[i]);
                  ub.insert(x->elts_lb[i]);
                }
		for(unsigned int i=0; i<x->elts_var.size; ++i) {
			if(x->get_index_var(i).get_solution_min()) {
				lb.insert(x->elts_var[i]);
				ub.insert(x->elts_var[i]);
			} else if(x->get_index_var(i).get_solution_max())
				ub.insert(x->elts_var[i]);
		}
		out << "{";
		for( set<int>::const_iterator i = ub.begin(); i != ub.end(); ++i) {
			if( i != ub.begin() ) out << ", ";
			out << *i;
		}
		out << "}";
	} else if (ai->isBool()) {
		out << (ai->getBool() ? "true" : "false");
	} else if (ai->isSet()) {
		AST::SetLit* s = ai->getSet();
		if (s->interval) {
			out << s->min << ".." << s->max;
		} else {
			out << "{";
			for (unsigned int i=0; i<s->s.size(); i++) {
				out << s->s[i] << (i < s->s.size()-1 ? ", " : "}");
			}
		}
	} else if (ai->isString()) {
		std::string s = ai->getString();
		for (unsigned int i=0; i<s.size(); i++) {
			if (s[i] == '\\' && i<s.size()-1) {
				switch (s[i+1]) {
				case 'n': out << "\n"; break;
				case '\\': out << "\\"; break;
				case 't': out << "\t"; break;
				default: out << "\\" << s[i+1];
				}
				i++;
			} else {
				out << s[i];
			}
		}
	}
}

void
Printer::print(std::ostream& out,
		Solver& solver,
		const IntVarArray& iv,
		const BoolVarArray& bv,
		const SetVarArray& sv) const {
	if (_output == NULL)
		return;
	for (unsigned int i=0; i< _output->a.size(); i++) {
		AST::Node* ai = _output->a[i];
		if (ai->isArray()) {
			AST::Array* aia = ai->getArray();
			int size = aia->a.size();
			out << "[";
			for (int j=0; j<size; j++) {
				printElem(out,solver, aia->a[j],iv,bv,sv);
				if (j<size-1)
					out << ", ";
			}
			out << "]";
		} else {
			printElem(out,solver,ai,iv,bv,sv);
		}
	}
}


Printer::~Printer(void) {
	delete _output;
}


void FlatZincModel::add_clause(Vector<Variable> pos ,  Vector<Variable> neg)
{

	clause_struct v;
	v.pos =pos;
	v.neg =neg;
	_clauses.add(v);

	if(! pos.empty())
		for (int i=0; i< pos. size ; i++)
			solver.add(pos[i]);

	if(! neg.empty())
		for (int i=0; i< neg. size ; i++)
			solver.add(neg[i]);
}

//encode the clause having all the positive literals in the vector of variables pos and the negative ones in neg
void FlatZincModel::encode_clause (Vector<Variable> pos ,  Vector<Variable> neg)
{

	/*We use the variable alone only when the size of the clause is 1.
	 *In that case we force it to take the value sign_alone
	 */
	Variable  alone;
	int sign_alone;
	Vector< Literal > clause;
	Literal  lit;
	clause.clear();

	if(! pos.empty())
	{
		for (int i=0; i< pos. size ; i++)
		{
			if (! pos[i].is_ground())
			{
				//	solver.add(pos[i]);
//				lit =  (2* pos[i].id()) +1;
				clause.add(literal(pos[i],1));
				if(clause.size ==1)
				{
					sign_alone=1;
					alone= pos[i];
				}
			}
			else
				if (pos[i].get_value())
					return;
		}
	}
	if(! neg.empty())
	{
		for (int i=0; i< neg.size ; i++)
		{
			if (!neg[i].is_ground())
			{
				//				solver.add(neg[i]);
				//lit =  2* neg[i].id();
				clause.add(literal(neg[i],0));
				if(clause.size ==1)
				{
					sign_alone=0;
					alone=neg[i];
				}
			}
			else
				if (!neg[i].get_value())
					return;
		}
	}
	if (clause.size ==1)
	{
		solver.add(alone==sign_alone);
		return ;
	}
	else
	{
		if (clause.size >0)
		{
			clause.sort();
//			std::cout<< "add clause : " <<clause<< std::endl;
			cnf.add(clause);
		}
		else if (clause.size ==0)
			solver.fail();
	}
}


void FlatZincModel::encode_clauses()
{
	int size =_clauses.size;
	cnf.clear();
	if(size)
		while(--size)
			encode_clause(_clauses[size].pos, _clauses[size].neg);

	for(int i=0; i<cnf.size; ++i)
		solver.add(cnf[i]);
}


}


FlatZinc::SolutionPrinter::SolutionPrinter(Printer *p, FlatZincModel *fm, Mistral::Solver *s) 
  : p_(p), fm_(fm), solver_(s) {
  //solver_->add((SolutionListener*)this);
}

void FlatZinc::SolutionPrinter::notify_solution() {
  fm_->print_solution(std::cout, *p_);
};

// STATISTICS: flatzinc-any
