
/*
  Mistral 2.0 is a constraint satisfaction and optimisation library
  Copyright (C) 2009  Emmanuel Hebrard
  
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

  The author can be contacted electronically at emmanuel.hebrard@gmail.com.
*/


#include <math.h>


#include <assert.h>

#include <mistral_constraint.hpp>
#include <mistral_variable.hpp>
#include <mistral_solver.hpp>






#define PROFILING_HEAD				\
  double __t_tag__ = 0;				\
  if(domain_type != CONST_VAR) {		\
    __t_tag__ = get_run_time();			\
  }						\

#define PROFILING_FOOT(method)						\
  if(domain_type != CONST_VAR) {					\
    __t_tag__ = get_run_time() - __t_tag__;				\
    int __idx__ = VARTYPE[(domain_type < DYN_VAR ? domain_type : BOOL_VAR)]; \
    variable->solver->statistics.prof_time[method][__idx__] += __t_tag__; \
    ++variable->solver->statistics.prof_num[method][__idx__];		\
  }									\



/**
   Constructors for variable
*/

// 
Mistral::Variable::Variable() {
  domain_type = DYN_VAR;
  variable = NULL;
}

// Constant variable
Mistral::Variable::Variable(const int value) {
  domain_type = CONST_VAR;
  variable = NULL;
  constant_value = value;
}

Mistral::Variable::Variable(VariableImplementation* impl, const int type) {
  domain_type = type;
  variable = impl;
}

Mistral::Variable::Variable(Expression* exp) {
  domain_type = EXPRESSION;
  expression = exp;
}

int BOOL_DOM = 3;

Mistral::Variable::Variable(const Vector< int >& values, const int type) {
  initialise_domain(values, type);
}

Mistral::Variable::Variable(const int* values, const int nvalues, const int type) {
	Vector<int> vals(nvalues);
	for(int i=0; i<nvalues; ++i) {
		vals.add(values[i]);
	}
	initialise_domain(vals, type);
}

Mistral::Variable::Variable(const IntStack& values, const int type) {
	Vector<int> vals(values.size);
	for(int i=0; i<values.size; ++i) {
		vals.add(values[i]);
	}
	initialise_domain(vals, type);
}

Mistral::Variable::Variable(const int lo, const int up, const int type) {
  initialise_domain(lo, up, type);
}

Mistral::Variable::Variable(const int lo, const int up, const Vector< int >& values, const int type) {
  initialise_domain(lo, up, values, type);
}

Mistral::Variable::Variable(const Variable& x) {
  initialise(x);
}

void Mistral::Variable::initialise(const Variable& x) {
  bool_domain = x.bool_domain;
  variable = x.variable;
}

Mistral::Variable::Variable(Variable X, bool h) {
  if(X.domain_type == RANGE_VAR) {
    int lb = X.get_initial_min();
    int ub = X.get_initial_max();

    initialise_domain(lb, ub, BITSET_VAR);
    variable->id = X.variable->id;
    variable->solver = X.variable->solver;
    ((VariableRange*)X.variable)->set_history((VariableBitmap*)variable);
  }
}

void Mistral::Variable::free_object() {
  if     (domain_type ==  BITSET_VAR) delete bitset_domain;
  else if(domain_type ==    LIST_VAR) delete list_domain;
  else if(domain_type ==   RANGE_VAR) delete range_domain;
  else if(domain_type == VIRTUAL_VAR) delete virtual_domain;
  else if(domain_type ==  EXPRESSION) delete expression;
  else if(domain_type !=   CONST_VAR) delete variable;
}


void Mistral::Variable::initialise_domain(const int lo, const int up, const int type) {

#ifdef _DEBUG_BUILD
  std::cout << "init domain from bounds [" << lo << "," << up << "] (type = " << domain2str(type) << " -> " ; //<< ")";
#endif 

  if(lo == up) {
    domain_type = CONST_VAR;
    constant_value = lo;
  } else if(type == EXPRESSION) {
    domain_type = EXPRESSION;
    expression = new Expression(lo, up);
  } else if((type & BOOL_VAR) && lo==0 && up==1) {
    bool_domain = &BOOL_DOM;
    variable = new VariableImplementation();
  } else if(type & RANGE_VAR) {
    domain_type = RANGE_VAR;
    range_domain = new VariableRange(lo, up);
  } else if(type & BITSET_VAR) {
    domain_type = BITSET_VAR;

    int nwords = 1+(up >> BitSet::EXP)-(lo >> BitSet::EXP);
#ifdef _BIT64
    if(nwords == 1) bitset_domain = new VariableWord<unsigned long long int, 1>(lo, up);
    else if(nwords == 2) bitset_domain = new VariableWord<unsigned long long int, 2>(lo, up);
    else if(nwords == 3) bitset_domain = new VariableWord<unsigned long long int, 3>(lo, up);
#else
    if(nwords == 1) bitset_domain = new VariableWord<unsigned int, 1>(lo, up);
    else if(nwords == 2) bitset_domain = new VariableWord<unsigned int, 2>(lo, up);
    else if(nwords == 3) bitset_domain = new VariableWord<unsigned int, 3>(lo, up);
    else if(nwords == 4) bitset_domain = new VariableWord<unsigned int, 4>(lo, up);
    else if(nwords == 5) bitset_domain = new VariableWord<unsigned int, 5>(lo, up);
#endif
    else bitset_domain = new VariableBitmap(lo, up);
  } else {
    domain_type = LIST_VAR;
    list_domain = new VariableList(lo, up);
  }

#ifdef _DEBUG_BUILD
  std::cout << domain2str(domain_type) << ": " << get_domain() << ")\n";
#endif

}


void Mistral::Variable::initialise_domain(const Vector< int >& values, const int type) {
  int min = values.front();
  int max = values.front();
  
  for(unsigned int i=1; i<values.size; ++i) {
    if(values[i] < min) min = values[i];
    if(values[i] > max) max = values[i];
  }
  
  initialise_domain(min, max, values, type);
}

void Mistral::Variable::initialise_domain(const int min, const int max, const Vector< int >& values, const int type) {

#ifdef _DEBUG_BUILD
  std::cout << "init domain from set [" << min << "," << max << "] " << values << " (type = " << domain2str(type) << " -> " ; //<< ")";
#endif 
	  //std::cout << type << " " << BITSET_VAR << " " << LIST_VAR << std::endl;

  if(max - min + 1 == (int)(values.size)) {
    initialise_domain(min, max, type);
  } else {
    if(type == EXPRESSION) {
		
		//std::cout << "->build expression" << std::endl;
		
      domain_type = EXPRESSION;
      expression = new Expression(min, max, values);
    } else if(type & BITSET_VAR) {
      domain_type = BITSET_VAR;
      int nwords = 1+(max >> BitSet::EXP)-(min >> BitSet::EXP);
      
#ifdef _BIT64
      if(nwords == 1) bitset_domain = new VariableWord<unsigned long long int, 1>(min, max, values);
      else if(nwords == 2) bitset_domain = new VariableWord<unsigned long long int, 2>(min, max, values);
      else if(nwords == 3) bitset_domain = new VariableWord<unsigned long long int, 3>(min, max, values);
#else
      if(nwords == 1) bitset_domain = new VariableWord<unsigned int, 1>(min, max, values);
      else if(nwords == 2) bitset_domain = new VariableWord<unsigned int, 2>(min, max, values);
      else if(nwords == 3) bitset_domain = new VariableWord<unsigned int, 3>(min, max, values);
      else if(nwords == 4) bitset_domain = new VariableWord<unsigned int, 4>(min, max, values);
      else if(nwords == 5) bitset_domain = new VariableWord<unsigned int, 5>(min, max, values);
#endif
      else {
		  
		  
		  
	bitset_domain = new VariableBitmap(min, max, values);
      }
	  
	  //std::cout << "->build bitsetvar " << bitset_domain->domain << std::endl;
	  
    } else {
      domain_type = LIST_VAR;
      list_domain = new VariableList(min, max, values);
    }
  }

#ifdef _DEBUG_BUILD
  std::cout << domain2str(domain_type) << ": " << get_domain() << ")\n";
#endif

}

bool Mistral::Variable::operator_equal(const Variable x) 
{
  if(x.domain_type == domain_type) {
    if(domain_type ==   CONST_VAR) return (constant_value == x.constant_value);
    else if(domain_type ==   EXPRESSION) return (expression == x.expression);
    else return (variable == x.variable);
  }
  return false;
}

bool Mistral::Variable::is_set_var() { 
  return domain_type == EXPRESSION && expression->is_set(); 
}

Mistral::Variable Mistral::Variable::get_var() {
  
  // SELF CHANGE
  // if(domain_type == EXPRESSION) {
  //   //Variable x = expression->self;
  //   //while(x.is_expression)
  //   return expression->self.get_var();
  // } else 
 
  if(domain_type == CONST_VAR || !variable->solver) {
    return *this;
  } else if(variable->id == -2) return expression->_self.get_var();
  return ((Solver*)(variable->solver))->variables[variable->id];
}

const Mistral::Variable Mistral::Variable::get_var() const {
  // SELF CHANGE
  // if(domain_type == EXPRESSION) {
  //   return expression->self.get_var();
  // } else 
  
  if(domain_type == CONST_VAR || !variable->solver) {
    return *this;
  }
  return ((Solver*)(variable->solver))->variables[variable->id];
}

Mistral::Solver* Mistral::Variable::get_solver() {
  //std::cout << (int*)(variable) << std::endl;
  return (Solver*)(variable->solver);
}



std::ostream& Mistral::Variable::display(std::ostream& os) const {
  if(domain_type == EXPRESSION) {
    expression->display(os);
  } else if(domain_type == CONST_VAR) {
    os << constant_value;
  } else if(variable) {
    int id = variable->id;
    
    if (domain_type == BITSET_VAR) {
      os << "x" ;
    } else if(domain_type == LIST_VAR) {
      os << "y" ;
    } 
    else if(domain_type == RANGE_VAR) {
      os << "r" ;
    }
    else  {
      os << "b" ;
    }

    if(variable->is_initialised()) {
      os << id ;
    } else {
      os << "_";
    }

  } else {

    os << "_";

  }

  return os;
  
  //if(*((int*)domain_type)  ==   CONST_VAR) os << constant_value;
  //else os = implementation->display(os);
  //return os;
}

//#define _DEBUG_BUILD true

void Mistral::Variable::initialise(Solver *s, const int level) {

  //std::cout << "call initialise on " << *this << std::endl;

#ifdef _DEBUG_BUILD
  for(int i=0; i<level; ++i) std::cout << "  " ;
  std::cout << "initialise " << *this << std::endl;
#endif
      

  if(domain_type == EXPRESSION) {

    if(!expression->is_initialised()) {
      
#ifdef _DEBUG_BUILD
      for(int i=0; i<=level; ++i) std::cout << "  " ;
      std::cout << "beg init " << expression->get_name();
      //display(std::cout);
      std::cout << std::endl;
#endif
      
      for(unsigned int i=0; i<expression->children.size; ++i) {
	expression->children[i].initialise(s, level+1);
      }

      if(level == 0 && !expression->children.empty()) {

#ifdef _DEBUG_BUILD
	std::cout << "-> constraint!" << std::endl;
#endif
	expression->extract_constraint(s);
      } else {

#ifdef _DEBUG_BUILD
	for(int i=0; i<=level; ++i) std::cout << "  " ;	
	std::cout << "-> predicate! (extract var)" << std::endl;
#endif

	expression->extract_variable(s);

	// SELF CHANGE
	Variable X = expression->_self;
	expression->id = (X.domain_type == CONST_VAR ? -2 : X.id());
	expression->solver = s;	

#ifdef _DEBUG_BUILD
	for(int i=0; i<=level; ++i) std::cout << "  " ;
	std::cout << " extracted " << X << " -> predicate " << std::endl;
#endif

	expression->extract_predicate(s);//);

      }

#ifdef _DEBUG_BUILD
      for(int i=0; i<=level; ++i) std::cout << "  " ;
      std::cout << "end init " ;
      display(std::cout);
      std::cout << std::endl;
#endif

      s->expression_store.add(expression);
    }
  } else {
    if(domain_type != CONST_VAR && variable->solver != s) {

#ifdef _DEBUG_BUILD
      for(int i=0; i<=level; ++i) std::cout << "  " ;
      std::cout << "add variable " ;
      display(std::cout);
      std::cout << " to the model" << std::endl;
#endif

      s->declare(this->get_var());
      //s->declare(*this);//->get_var());
      s->sequence.declare(*this);
    } 
  }
}



Mistral::Event Mistral::Variable::setValue( const int val ) 
{
  // val should be 1 or 2
  Event evt = VALUE_C;
  int dom = *bool_domain;
  
  if( val == dom ) return NO_EVENT;
  else if( dom<3 // || val<1 || val>2
	   ) return FAIL_EVENT;
  
  *bool_domain = val;

  evt |= (val==1 ? UB_EVENT : LB_EVENT);
  variable->trigger_event_and_save(evt);

  return evt;
}


Mistral::Event Mistral::Variable::setState( const int vals ) 
{
  // vals should be 1, 2 or 3
  Event evt = VALUE_C;
  int dom = *bool_domain;
  int ndom = vals&dom;
  
  if( ndom == dom ) return NO_EVENT;
  else if( !ndom ) return FAIL_EVENT;
  
  *bool_domain = ndom;
  evt |= (ndom==1 ? UB_EVENT : LB_EVENT);

  variable->trigger_event_and_save(evt);

  return evt;
}

int Mistral::Variable::get_solution_int_value() const {
  int value = 0;
  if(is_initialised()) {
    value = variable->get_solution_int_value();
  } else {
    value = get_first();
  }
  return value;
}

std::string Mistral::Variable::get_solution_str_value() const {
  std::ostringstream ret_str;

  if(is_initialised()) {
    ret_str << variable->get_solution_str_value();
  } else {
    ret_str << get_first();
  }

  return ret_str.str();
}

int Mistral::Variable::get_solution_min() const {
  int value = 0;
  if(is_initialised()) {
    value = variable->get_solution_min();
  } else {
    value = get_min();
  }
  return value;
}

int Mistral::Variable::get_solution_max() const {
  int value = 0;
  if(is_initialised()) {
    value = variable->get_solution_max();
  } else {
    value = get_max();
  }
  return value;
}

int Mistral::Variable::get_value() const {
  if     (domain_type ==  BITSET_VAR) return bitset_domain->get_value();
  else if(domain_type ==    LIST_VAR) return list_domain->get_value();
  else if(domain_type ==   RANGE_VAR) return range_domain->get_value();
  //else if(domain_type == VIRTUAL_VAR) return virtual_domain->get_value();
  else if(domain_type ==   CONST_VAR) return constant_value;
  else if(domain_type ==   EXPRESSION) return expression->get_self().get_value();
  else  return (*bool_domain-1);
}

unsigned int Mistral::Variable::get_size() const {

#ifdef _PROFILING_PRIMITIVE
  PROFILING_HEAD
#endif
      
    unsigned int r_size = 0;

  if     (domain_type ==  BITSET_VAR) r_size = bitset_domain->get_size();
  else if(domain_type ==    LIST_VAR) r_size = list_domain->get_size();
  else if(domain_type ==   RANGE_VAR) r_size = range_domain->get_size();
  //else if(domain_type == VIRTUAL_VAR) r_size = virtual_domain->get_size();
  else if(domain_type ==   CONST_VAR) r_size = 1;
  else if(domain_type ==   EXPRESSION) r_size = expression->get_self().get_size();
  else  r_size = ((*bool_domain+1)/2);

#ifdef _PROFILING_PRIMITIVE
  PROFILING_FOOT(_m_get_size_)
#endif

    return r_size;
}

unsigned int Mistral::Variable::get_reduction() const {

#ifdef _PROFILING_PRIMITIVE
  PROFILING_HEAD
#endif
      
    unsigned int r_red = 0;

  if     (domain_type ==  BITSET_VAR) r_red = bitset_domain->get_reduction();
  else if(domain_type ==    LIST_VAR) r_red = list_domain->get_reduction();
  else if(domain_type ==   RANGE_VAR) r_red = range_domain->get_reduction();
  //else if(domain_type == VIRTUAL_VAR) r_size = virtual_domain->get_size();
  else if(domain_type ==   CONST_VAR) r_red = 0;
  else if(domain_type ==   EXPRESSION) r_red = expression->get_self().get_reduction();
  else  r_red = (variable->assigned_at_last_level());

#ifdef _PROFILING_PRIMITIVE
  PROFILING_FOOT(_m_get_red_)
#endif

    return r_red;
}

/// Returns the degree (number of constraints)
unsigned int Mistral::Variable::get_degree() const {

#ifdef _PROFILING_PRIMITIVE
  PROFILING_HEAD
#endif

    unsigned int r_degree = NO_EVENT;

  if(domain_type ==   CONST_VAR) r_degree = 0;
  else if(domain_type ==   EXPRESSION) r_degree = expression->get_self().get_degree();
  else r_degree = ((Solver*)(variable->solver))->constraint_graph[variable->id].size();


#ifdef _PROFILING_PRIMITIVE
  PROFILING_FOOT(_m_get_degree_)
#endif

    return r_degree;
}


std::string Mistral::Variable::get_domain(const bool latex) const {
  std::ostringstream buf;

  if(latex && domain_type ==  BITSET_VAR) buf << "\\";

  if     (domain_type ==  BITSET_VAR) buf << bitset_domain->domain;
  else if(domain_type ==    LIST_VAR) buf << list_domain->domain;
  else if(domain_type ==   RANGE_VAR) {
    //Mistral::VariableRange* r = (Mistral::VariableRange*)implementation;
    if(range_domain->get_min() == range_domain->get_max())
      buf << range_domain->get_min() ;
    else if(range_domain->get_min() == range_domain->get_max()-1)
      buf << "[" << range_domain->get_min() << "," <<  range_domain->get_max() << "]";
    else 
      buf << "[" << range_domain->get_min() << ".." <<  range_domain->get_max() << "]";
  }
  //else if(domain_type == VIRTUAL_VAR) return virtual_domain->get_domain();
  else if(domain_type ==   CONST_VAR)  {
    buf << constant_value ;
  }
  else if(domain_type ==   BOOL_VAR)  {
    buf << "[0,1]";
  } 
  else if(domain_type ==   EXPRESSION) {

    // buf << "e" ;
    // if(expression->is_set()) buf << "s" ;

    if(expression->is_set()) ((SetExpression*)(expression))->display(buf, true);
    else buf << expression->get_self().get_domain();
  } else {
    if(*bool_domain == 3) buf << "[0,1]";
    else if(*bool_domain == 2) buf << "1";
    else buf << "0";
  }

  if(latex && domain_type ==  BITSET_VAR) buf << "\b\\}";


  return buf.str();
}


std::string Mistral::Variable::get_history() const {
  std::ostringstream buf;
  if     (domain_type ==  BITSET_VAR) buf << bitset_domain->get_history();
  else if(domain_type ==   RANGE_VAR) buf << range_domain->get_history();
  else if(domain_type >      DYN_VAR && *bool_domain < 3) buf << "[0,1]";
  return buf.str();
}

bool Mistral::Variable::is_boolean() const {
  bool of_the_living_dead = false;
  if(domain_type > DYN_VAR) of_the_living_dead = true;
  else if(domain_type ==  BITSET_VAR) {
    of_the_living_dead =  (bitset_domain->get_min()==0 && bitset_domain->get_max()==1);
  } else if(domain_type ==    LIST_VAR)  {
    int f = list_domain->get_first(); 
    if(f == 0)
      of_the_living_dead = list_domain->get_last()==1;
    else if(f == 1)
      of_the_living_dead = list_domain->get_last()==0;
  } else if(domain_type ==   RANGE_VAR)  {
    of_the_living_dead =  (range_domain->get_min()==0 && range_domain->get_max()==1);
  } else if(domain_type ==   CONST_VAR)  {
    of_the_living_dead =  constant_value<=1 && constant_value>=0;
  } else if(domain_type ==   EXPRESSION)  {
    of_the_living_dead =  expression->get_self().is_boolean();
  } 

  return of_the_living_dead;
}



int Mistral::Variable::get_first() const {

#ifdef _PROFILING_PRIMITIVE
  PROFILING_HEAD
#endif

    
    int of_the_living_dead = 0;
  if     (domain_type ==  BITSET_VAR) {
    //std::cout << "bitset " << bitset_domain->get_min() << std::endl;
    of_the_living_dead =  bitset_domain->get_first();
  } else if(domain_type ==    LIST_VAR)  {
    //std::cout << "list" << std::endl;
    of_the_living_dead =  list_domain->get_first();
  } else if(domain_type ==   RANGE_VAR)  {
    //std::cout << "range " << range_domain->get_min() << std::endl;
    of_the_living_dead =  range_domain->get_first();
    //else if(domain_type == VIRTUAL_VAR) return virtual_domain->get_min();
  } else if(domain_type ==   CONST_VAR)  {
    //std::cout << "constant" << std::endl;
    of_the_living_dead =  constant_value;
  } else if(domain_type ==   EXPRESSION)  {
    //std::cout << "expression" << expression->get_self().get_min() << std::endl;
    of_the_living_dead =  expression->get_self().get_first();
  } else  of_the_living_dead = !(*bool_domain & 1);


#ifdef _PROFILING_PRIMITIVE
  PROFILING_FOOT(_m_get_first_)
#endif

    return of_the_living_dead;
}

int Mistral::Variable::get_last() const {

#ifdef _PROFILING_PRIMITIVE
  PROFILING_HEAD
#endif

    int of_the_mummy = 0;

  if     (domain_type ==  BITSET_VAR) of_the_mummy = bitset_domain->get_last();
  else if(domain_type ==    LIST_VAR) of_the_mummy = list_domain->get_last();
  else if(domain_type ==   RANGE_VAR) of_the_mummy = range_domain->get_last();
  //else if(domain_type == VIRTUAL_VAR) of_the_mummy = virtual_domain->get_max();
  else if(domain_type ==   CONST_VAR) of_the_mummy = constant_value;
  else if(domain_type ==   EXPRESSION) of_the_mummy = expression->get_self().get_last();
  else  of_the_mummy = (*bool_domain >> 1);

#ifdef _PROFILING_PRIMITIVE
  PROFILING_FOOT(_m_get_last_)
#endif

    return of_the_mummy;
}


int Mistral::Variable::get_min() const {

#ifdef _PROFILING_PRIMITIVE
  PROFILING_HEAD
#endif

    // if     (domain_type ==  BITSET_VAR) {
    //   //std::cout << "bitset " << bitset_domain->get_min() << std::endl;
    //   return bitset_domain->get_min();
    // } else if(domain_type ==    LIST_VAR)  {
    //   //std::cout << "list" << std::endl;
    //   return list_domain->get_min();
    // } else if(domain_type ==   RANGE_VAR)  {
    //   //std::cout << "range " << range_domain->get_min() << std::endl;
    //   return range_domain->get_min();
    // //else if(domain_type == VIRTUAL_VAR) return virtual_domain->get_min();
    // } else if(domain_type ==   CONST_VAR)  {
    //   //std::cout << "constant" << std::endl;
    //   return constant_value;
    // } else if(domain_type ==   EXPRESSION)  {
    //   //std::cout << "expression" << expression->get_self().get_min() << std::endl;
    //   return expression->get_self().get_min();
    // } else  return (!(*bool_domain & 1));

      
    //   std::cout << "get min of "  ;
    // std::cout.flush();
    // display(std::cout);
    // std::cout << std::endl;
    
    int of_the_living_dead = 0;
  if     (domain_type ==  BITSET_VAR) {
    //std::cout << "bitset " << bitset_domain->get_min() << std::endl;
    of_the_living_dead =  bitset_domain->get_min();
  } else if(domain_type ==    LIST_VAR)  {
    //std::cout << "list" << std::endl;
    of_the_living_dead =  list_domain->get_min();
  } else if(domain_type ==   RANGE_VAR)  {
    //std::cout << "range " << range_domain->get_min() << std::endl;
    of_the_living_dead =  range_domain->get_min();
    //else if(domain_type == VIRTUAL_VAR) return virtual_domain->get_min();
  } else if(domain_type ==   CONST_VAR)  {
    //std::cout << "constant" << std::endl;
    of_the_living_dead =  constant_value;
  } else if(domain_type ==   EXPRESSION)  {
    //std::cout << "expression" << expression->get_self().get_min() << std::endl;
    of_the_living_dead =  expression->get_self().get_min();
  } else  of_the_living_dead = !(*bool_domain & 1);


#ifdef _PROFILING_PRIMITIVE
  PROFILING_FOOT(_m_get_min_)
#endif

    return of_the_living_dead;
}

int Mistral::Variable::get_max() const {

#ifdef _PROFILING_PRIMITIVE
  PROFILING_HEAD
#endif

    int of_the_mummy = 0;

  if     (domain_type ==  BITSET_VAR) of_the_mummy = bitset_domain->get_max();
  else if(domain_type ==    LIST_VAR) of_the_mummy = list_domain->get_max();
  else if(domain_type ==   RANGE_VAR) of_the_mummy = range_domain->get_max();
  //else if(domain_type == VIRTUAL_VAR) of_the_mummy = virtual_domain->get_max();
  else if(domain_type ==   CONST_VAR) of_the_mummy = constant_value;
  else if(domain_type ==   EXPRESSION) of_the_mummy = expression->get_self().get_max();
  else  of_the_mummy = (*bool_domain >> 1);

#ifdef _PROFILING_PRIMITIVE
  PROFILING_FOOT(_m_get_max_)
#endif

    return of_the_mummy;
}

int Mistral::Variable::get_initial_min() const {

#ifdef _PROFILING_PRIMITIVE
  PROFILING_HEAD
#endif

    int r_min = 0;

  if     (domain_type ==  BITSET_VAR) r_min = bitset_domain->get_initial_min();
  else if(domain_type ==    LIST_VAR) r_min = list_domain->get_initial_min();
  else if(domain_type ==   RANGE_VAR) r_min = range_domain->get_initial_min();
  //else if(domain_type == VIRTUAL_VAR) r_min = virtual_domain->get_initial_min();
  else if(domain_type ==   CONST_VAR) r_min = constant_value;
  else if(domain_type ==   EXPRESSION) r_min = expression->get_self().get_initial_min();
  //else  r_min = 0;

#ifdef _PROFILING_PRIMITIVE
  PROFILING_FOOT(_m_get_initial_min_)
#endif

    return r_min;
}

int Mistral::Variable::get_initial_max() const {

#ifdef _PROFILING_PRIMITIVE
  PROFILING_HEAD
#endif

    int r_max = 1;

  if     (domain_type ==  BITSET_VAR) r_max = bitset_domain->get_initial_max();
  else if(domain_type ==    LIST_VAR) r_max = list_domain->get_initial_max();
  else if(domain_type ==   RANGE_VAR) r_max = range_domain->get_initial_max();
  //else if(domain_type == VIRTUAL_VAR) r_max = virtual_domain->get_initial_max();
  else if(domain_type ==   CONST_VAR) r_max = constant_value;
  else if(domain_type ==   EXPRESSION) r_max = expression->get_self().get_initial_max();
  //else  r_max = 1;

#ifdef _PROFILING_PRIMITIVE
  PROFILING_FOOT(_m_get_initial_max_)
#endif

    return r_max;
}

int Mistral::Variable::get_min_pos() const {
  
#ifdef _PROFILING_PRIMITIVE
  PROFILING_HEAD
#endif
    
    int r_min = 0;

  if     (domain_type ==  BITSET_VAR) r_min = bitset_domain->get_min_pos();
  else if(domain_type ==    LIST_VAR) r_min = list_domain->get_min_pos();
  else if(domain_type ==   RANGE_VAR) r_min = range_domain->get_min_pos();
  //else if(domain_type == VIRTUAL_VAR) r_min = virtual_domain->get_min_pos();
  else if(domain_type ==   CONST_VAR) r_min = constant_value;
  else if(domain_type ==   EXPRESSION) r_min = expression->get_self().get_min_pos();
  else  r_min = (*bool_domain >> 1); //(!(*bool_domain & 1));

#ifdef _PROFILING_PRIMITIVE
  PROFILING_FOOT(_m_get_min_pos_)
#endif

    return r_min;

}

int Mistral::Variable::get_max_neg() const {

#ifdef _PROFILING_PRIMITIVE
  PROFILING_HEAD
#endif

    int r_max = 0;

  if     (domain_type ==  BITSET_VAR) r_max = bitset_domain->get_max_neg();
  else if(domain_type ==    LIST_VAR) r_max = list_domain->get_max_neg();
  else if(domain_type ==   RANGE_VAR) r_max = range_domain->get_max_neg();
  //else if(domain_type == VIRTUAL_VAR) r_max = virtual_domain->get_max_neg();
  else if(domain_type ==   CONST_VAR) r_max = constant_value;
  else if(domain_type ==   EXPRESSION) r_max = expression->get_self().get_max_neg();
  else  r_max = (!(*bool_domain & 1));

#ifdef _PROFILING_PRIMITIVE
  PROFILING_FOOT(_m_get_max_neg_)
#endif

    return r_max;
}

int Mistral::Variable::next(const int v) const {

#ifdef _PROFILING_PRIMITIVE
  PROFILING_HEAD
#endif

    int r_val = v;

  if     (domain_type ==  BITSET_VAR) r_val = bitset_domain->next(v);
  else if(domain_type ==    LIST_VAR) r_val = list_domain->next(v);
  else if(domain_type ==   RANGE_VAR) r_val = range_domain->next(v);
  //else if(domain_type == VIRTUAL_VAR) r_val = virtual_domain->next(v);
  else if(domain_type ==   CONST_VAR) r_val = constant_value;
  else if(domain_type ==   EXPRESSION) r_val = expression->get_self().next(v);
  else  r_val = (*bool_domain >> 1);

#ifdef _PROFILING_PRIMITIVE
  PROFILING_FOOT(_m_next_)
#endif

    return r_val;
}

int Mistral::Variable::prev(const int v) const {

#ifdef _PROFILING_PRIMITIVE
  PROFILING_HEAD
#endif

    int r_val = v;

  if     (domain_type ==  BITSET_VAR) r_val = bitset_domain->prev(v);
  else if(domain_type ==    LIST_VAR) r_val = list_domain->prev(v);
  else if(domain_type ==   RANGE_VAR) r_val = range_domain->prev(v);
  //else if(domain_type == VIRTUAL_VAR) r_val = virtual_domain->prev(v);
  else if(domain_type ==   CONST_VAR) r_val = constant_value;
  else if(domain_type ==   EXPRESSION) r_val = expression->get_self().prev(v);
  else  r_val = (!(*bool_domain & 1));

#ifdef _PROFILING_PRIMITIVE
  PROFILING_FOOT(_m_prev_)
#endif

    return r_val;
}

bool Mistral::Variable::is_range() const {

#ifdef _PROFILING_PRIMITIVE
  PROFILING_HEAD
#endif

    bool answer = true;

  if     (domain_type ==  BITSET_VAR) answer = bitset_domain->is_range();
  else if(domain_type ==    LIST_VAR) answer = list_domain->is_range();
  //else if(domain_type ==   RANGE_VAR) answer = range_domain->is_range();
  //else if(domain_type == VIRTUAL_VAR) answer = virtual_domain->is_range();
  //else if(domain_type ==   CONST_VAR) answer = true;
  else if(domain_type ==   EXPRESSION) answer = expression->get_self().is_range();
  //else answer = true;

#ifdef _PROFILING_PRIMITIVE
  PROFILING_FOOT(_m_is_range_)
#endif

    return answer;
}

// // bool Mistral::Variable::is_ground(const Expression *x) const {
// //   return x->self.is_ground();
// // }
//   bool Mistral::Variable::is_constant() const {

// #ifdef _PROFILING_PRIMITIVE
//     PROFILING_HEAD
// #endif

//       bool answer = true;

//     // std::cout << (int*)(variable)
//     // 	      << std::endl
//     // 	      << domain_type << " == " 
//     // 	      << BITSET_VAR << "? "
//     // 	      << LIST_VAR << "? "
//     // 	      << RANGE_VAR << "? "
//     // 	      << CONST_VAR << "? " << std::endl;


//     //std::cout << domain2str(domain_type) << std::endl;

//     if     (domain_type ==  BITSET_VAR) answer = bitset_domain->is_ground();
//     else if(domain_type ==    LIST_VAR) answer = list_domain->is_ground();
//     else if(domain_type ==   RANGE_VAR) answer = range_domain->is_ground();
//     //else if(domain_type == VIRTUAL_VAR) answer = virtual_domain->is_ground();
//     else if(domain_type ==   CONST_VAR) answer = true;
//     else if(domain_type ==   EXPRESSION) answer = expression->get_self().is_ground();
//     else  answer = (*bool_domain != 3);

// #ifdef _PROFILING_PRIMITIVE
//     PROFILING_FOOT(_m_is_ground_)
// #endif

//     return answer;
//   }


bool Mistral::Variable::is_ground() const {

#ifdef _PROFILING_PRIMITIVE
  PROFILING_HEAD
#endif

    bool answer = true;

  // std::cout << (int*)(variable)
  // 	      << std::endl
  // 	      << domain_type << " == " 
  // 	      << BITSET_VAR << "? "
  // 	      << LIST_VAR << "? "
  // 	      << RANGE_VAR << "? "
  // 	      << CONST_VAR << "? " << std::endl;


  //std::cout << domain2str(domain_type) << std::endl;

  if     (domain_type ==  BITSET_VAR) answer = bitset_domain->is_ground();
  else if(domain_type ==    LIST_VAR) answer = list_domain->is_ground();
  else if(domain_type ==   RANGE_VAR) answer = range_domain->is_ground();
  //else if(domain_type == VIRTUAL_VAR) answer = virtual_domain->is_ground();
  else if(domain_type ==   CONST_VAR) answer = true;
  else if(domain_type ==   EXPRESSION) answer = (id() == -2 || expression->get_self().is_ground());
  else  answer = (*bool_domain != 3);

#ifdef _PROFILING_PRIMITIVE
  PROFILING_FOOT(_m_is_ground_)
#endif

    return answer;
}

bool Mistral::Variable::equal(const int v) const {

#ifdef _PROFILING_PRIMITIVE
  PROFILING_HEAD
#endif

    bool answer = true;

  if     (domain_type ==  BITSET_VAR) answer = bitset_domain->equal(v);
  else if(domain_type ==    LIST_VAR) answer = list_domain->equal(v);
  else if(domain_type ==   RANGE_VAR) answer = range_domain->equal(v);
  else if(domain_type == VIRTUAL_VAR) answer = virtual_domain->equal(v);
  else if(domain_type ==   CONST_VAR) answer = (constant_value == v);
  else if(domain_type ==   EXPRESSION) answer = expression->get_self().equal(v);
  else  answer = (*bool_domain-1 == v);

#ifdef _PROFILING_PRIMITIVE
  PROFILING_FOOT(_m_equal_)
#endif

    return answer;
}

bool Mistral::Variable::contain(const int v) const {

#ifdef _PROFILING_PRIMITIVE
  PROFILING_HEAD
#endif

    bool answer = true;

  if     (domain_type ==  BITSET_VAR) answer = bitset_domain->contain(v);
  else if(domain_type ==    LIST_VAR) answer = list_domain->contain(v);
  else if(domain_type ==   RANGE_VAR) answer = range_domain->contain(v);
  else if(domain_type == VIRTUAL_VAR) answer = virtual_domain->contain(v);
  else if(domain_type ==   CONST_VAR) answer = (constant_value == v);
  else if(domain_type ==   EXPRESSION) answer = expression->get_self().contain(v);
  else  answer = (!(v >> 1) && (*bool_domain & (v+1)));

#ifdef _PROFILING_PRIMITIVE
  PROFILING_FOOT(_m_contain_)
#endif

    return answer;
}

bool Mistral::Variable::intersect(const int lo, const int up) const {

#ifdef _PROFILING_PRIMITIVE
  PROFILING_HEAD
#endif

    bool answer = true;

  if     (domain_type ==  BITSET_VAR) answer = bitset_domain->intersect(lo, up);
  else if(domain_type ==    LIST_VAR) answer = list_domain->intersect(lo, up);
  else if(domain_type ==   RANGE_VAR) answer = range_domain->intersect(lo, up);
  else if(domain_type == VIRTUAL_VAR) answer = virtual_domain->intersect(lo, up);
  else if(domain_type ==   CONST_VAR) answer = (constant_value >= lo && constant_value <= up);
  else if(domain_type ==   EXPRESSION) answer = expression->get_self().intersect(lo, up);
  else  answer = (((lo<=0) | (2*(up>0))) & *bool_domain);

#ifdef _PROFILING_PRIMITIVE
  PROFILING_FOOT(_m_intersect_range_)
#endif

    return answer;
}

bool Mistral::Variable::intersect(const Interval I) const {
  return intersect(I.min, I.max);
}

bool Mistral::Variable::included(const int lo, const int up) const {

#ifdef _PROFILING_PRIMITIVE
  PROFILING_HEAD
#endif

    bool answer = true;

  if     (domain_type ==  BITSET_VAR) answer = bitset_domain->included(lo, up);
  else if(domain_type ==    LIST_VAR) answer = list_domain->included(lo, up);
  else if(domain_type ==   RANGE_VAR) answer = range_domain->included(lo, up);
  else if(domain_type == VIRTUAL_VAR) answer = virtual_domain->included(lo, up);
  else if(domain_type ==   CONST_VAR) answer = (constant_value >= lo && constant_value <= up);
  else if(domain_type ==   EXPRESSION) answer = expression->get_self().included(lo, up);
  else  {
    int state = *bool_domain;
    answer = ( up >= (state >> 1) && (lo <= !(state & 1)) );
  }

#ifdef _PROFILING_PRIMITIVE
  PROFILING_FOOT(_m_included_range_)
#endif

    return answer;
}

bool Mistral::Variable::included(const Interval I) const {
  return included(I.min, I.max);
}

bool Mistral::Variable::includes(const int lo, const int up) const {
  
#ifdef _PROFILING_PRIMITIVE
  PROFILING_HEAD
#endif

    bool answer = true;

  if     (domain_type ==  BITSET_VAR) answer = bitset_domain->includes(lo, up);
  else if(domain_type ==    LIST_VAR) answer = list_domain->includes(lo, up);
  else if(domain_type ==   RANGE_VAR) answer = range_domain->includes(lo, up);
  else if(domain_type == VIRTUAL_VAR) answer = virtual_domain->includes(lo, up);
  else if(domain_type ==   CONST_VAR) answer = (constant_value == lo && constant_value == up);
  else if(domain_type ==   EXPRESSION) answer = expression->get_self().includes(lo, up);
  else  {
    int state = *bool_domain;
    answer = ( up <= (state >> 1) && (lo >= !(state & 1)) );
  }

#ifdef _PROFILING_PRIMITIVE
  PROFILING_FOOT(_m_includes_range_)
#endif

    return answer;
}

bool Mistral::Variable::includes(const Interval I) const {
  return includes(I.min, I.max);
}

bool Mistral::Variable::intersect(const BitSet& s) const {

#ifdef _PROFILING_PRIMITIVE
  PROFILING_HEAD
#endif

    bool answer = true;

  if     (domain_type ==  BITSET_VAR) answer = bitset_domain->intersect(s);
  else if(domain_type ==    LIST_VAR) answer = list_domain->intersect(s);
  else if(domain_type ==   RANGE_VAR) answer = range_domain->intersect(s);
  else if(domain_type == VIRTUAL_VAR) answer = virtual_domain->intersect(s);
  else if(domain_type ==   CONST_VAR) answer = (s.contain(constant_value));
  else if(domain_type ==   EXPRESSION) answer = expression->get_self().intersect(s);
  else  answer = s.intersect(*bool_domain);

#ifdef _PROFILING_PRIMITIVE
  PROFILING_FOOT(_m_intersect_set_)
#endif

    return answer;
}

bool Mistral::Variable::included(const BitSet& s) const {

#ifdef _PROFILING_PRIMITIVE
  PROFILING_HEAD
#endif
    bool answer = true;

  if     (domain_type ==  BITSET_VAR) answer = bitset_domain->included(s);
  else if(domain_type ==    LIST_VAR) answer = list_domain->included(s);
  else if(domain_type ==   RANGE_VAR) answer = range_domain->included(s);
  else if(domain_type == VIRTUAL_VAR) answer = virtual_domain->included(s);
  else if(domain_type ==   CONST_VAR) answer = (s.contain(constant_value));
  else if(domain_type ==   EXPRESSION) answer = expression->get_self().included(s);
  else  answer = s.includes(*bool_domain);

#ifdef _PROFILING_PRIMITIVE
  PROFILING_FOOT(_m_included_set_)
#endif

    return answer;
}

bool Mistral::Variable::includes(const BitSet& s) const {

#ifdef _PROFILING_PRIMITIVE
  PROFILING_HEAD
#endif
    bool answer = true;

  if     (domain_type ==  BITSET_VAR) answer = bitset_domain->includes(s);
  else if(domain_type ==    LIST_VAR) answer = list_domain->includes(s);
  else if(domain_type ==   RANGE_VAR) answer = range_domain->includes(s);
  else if(domain_type == VIRTUAL_VAR) answer = virtual_domain->includes(s);
  else if(domain_type ==   CONST_VAR) answer = (s.size() == 1 && s.contain(constant_value));
  else if(domain_type ==   EXPRESSION) answer = expression->get_self().includes(s);
  else  answer = s.included(*bool_domain);

#ifdef _PROFILING_PRIMITIVE
  PROFILING_FOOT(_m_includes_set_)
#endif

    return answer;
}

bool Mistral::Variable::intersect(const Mistral::Variable& x) const {

  bool answer = true;

  if(is_ground()) answer = x.contain(get_first());
  else if(x.is_ground()) answer = contain(x.get_first());
  else if(is_range()) answer = x.intersect(get_min(), get_max());
  else if(x.is_range()) answer = intersect(x.get_min(), x.get_max());
  else if(domain_type ==  BITSET_VAR)
    answer = x.intersect(bitset_domain->domain.values);
  // else if(domain_type ==    LIST_VAR) {
  //   std::cout << "TODO! (intersect - list)" << std::endl;
  //   exit(1);
  // }
  // // else if(domain_type == VIRTUAL_VAR) {

  // // }
  else 
    answer = //x.intersect(expression->get_self()); //
      expression->get_self().intersect(x);

  return answer;
}

bool Mistral::Variable::included(const Mistral::Variable& x) const {

  bool answer = true;

  if(is_ground()) answer = x.contain(get_first());
  else if(x.is_ground()) answer = equal(x.get_first());
  else if(is_range()) answer = x.includes(get_min(), get_max());
  else if(x.is_range()) answer = included(x.get_min(), x.get_max()); 
  else if(domain_type ==  BITSET_VAR)
    answer = x.includes(bitset_domain->domain.values);
  //std::cout << "TODO! (included)" << std::endl;
  //answer = true;
  else 
    answer = //x.includes(expression->get_self());//.included(x);
      expression->get_self().included(x);

  return answer;
}

bool Mistral::Variable::includes(const Mistral::Variable& x) const {

  bool answer = true;

  if(is_ground()) answer = x.equal(get_first());
  else if(x.is_ground()) answer = contain(x.get_first());
  else if(is_range()) answer = x.included(get_min(), get_max());
  else if(x.is_range()) answer = includes(x.get_min(), x.get_max()); 
  else if(domain_type ==  BITSET_VAR)
    answer = x.included(bitset_domain->domain.values);
  // std::cout << "TODO!" << std::endl;
  // answer = true;
  else 
    answer = 
      expression->get_self().includes(x);

  return answer;
}

void Mistral::Variable::intersect_to( BitSet& s ) const {

#ifdef _PROFILING_PRIMITIVE
  PROFILING_HEAD
#endif

    if     (domain_type ==  BITSET_VAR) bitset_domain->intersect_to(s);
    else if(domain_type ==    LIST_VAR) list_domain->intersect_to(s);
    else if(domain_type ==   RANGE_VAR) range_domain->intersect_to(s);
    else if(domain_type == VIRTUAL_VAR) virtual_domain->intersect_to(s);
    else if(domain_type ==   CONST_VAR) {
      if(s.contain(constant_value)) {
	s.clear();
	s.add(constant_value);
      } else s.clear();
    }
    else if(domain_type ==   EXPRESSION) expression->get_self().intersect_to(s);
    else  s.intersect_with(*bool_domain);

#ifdef _PROFILING_PRIMITIVE
  PROFILING_FOOT(_m_intersect_to_)
#endif

    }

void Mistral::Variable::union_to( BitSet& s ) const {

#ifdef _PROFILING_PRIMITIVE
  PROFILING_HEAD
#endif
    
    if     (domain_type ==  BITSET_VAR) bitset_domain->union_to(s);
    else if(domain_type ==    LIST_VAR) list_domain->union_to(s);
    else if(domain_type ==   RANGE_VAR) range_domain->union_to(s);
    else if(domain_type == VIRTUAL_VAR) {
      
      // std::cout << domain_type << std::endl;

      // std::cout << "make the union of " << *this // << " in " << get_domain()
      // 		<< " into " << s << std::endl;

      virtual_domain->union_to(s);
    } else if(domain_type ==   CONST_VAR) s.add(constant_value);
    else if(domain_type ==   EXPRESSION) expression->get_self().union_to(s);
    else s.union_with(*bool_domain);

#ifdef _PROFILING_PRIMITIVE
  PROFILING_FOOT(_m_union_to_)
#endif

    }


void Mistral::Variable::put_negation_in( BitSet& s ) const {
  
#ifdef _PROFILING_PRIMITIVE
  PROFILING_HEAD
#endif
    
    if(is_ground()) s.add(-get_value());
    else if(domain_type ==  BITSET_VAR) bitset_domain->put_negation_in(s);
    else if(domain_type ==    LIST_VAR) list_domain->put_negation_in(s);
    else if(domain_type ==   RANGE_VAR) range_domain->put_negation_in(s);
    else if(domain_type == VIRTUAL_VAR) {
      
      // std::cout << domain_type << std::endl;

      // std::cout << "make the union of " << *this // << " in " << get_domain()
      // 		<< " into " << s << std::endl;

      virtual_domain->put_negation_in(s);
    } //else if(domain_type ==   CONST_VAR) s.add(-constant_value);
    else if(domain_type ==   EXPRESSION) expression->get_self().put_negation_in(s);
    else //s.union_with(*bool_domain); 
      {
	s.add(0);
	s.add(-1);
      }

#ifdef _PROFILING_PRIMITIVE
  PROFILING_FOOT(_m_put_negation_in_)
#endif

    }

Mistral::Event Mistral::Variable::remove(const int v) {

#ifdef _PROFILING_PRIMITIVE
  PROFILING_HEAD
#endif
    Event evt = NO_EVENT;

  if     (domain_type ==  BITSET_VAR) evt = bitset_domain->remove(v);
  else if(domain_type ==    LIST_VAR) evt = list_domain->remove(v);
  else if(domain_type ==   RANGE_VAR) evt = range_domain->remove(v);
  else if(domain_type == VIRTUAL_VAR) evt = virtual_domain->remove(v);
  else if(domain_type ==   CONST_VAR) evt = (constant_value == v ? FAIL_EVENT : NO_EVENT);
  else if(domain_type ==   EXPRESSION) evt = expression->get_self().remove(v);
  else {
    evt = (v<0||v>1 ? NO_EVENT : setValue(2-v));
  }

#ifdef _PROFILING_PRIMITIVE
  PROFILING_FOOT(_m_remove_)
#endif

    return evt;
}

Mistral::Event Mistral::Variable::set_domain(const int v) {

#ifdef _PROFILING_PRIMITIVE
  PROFILING_HEAD
#endif
    Event evt = NO_EVENT;

  if     (domain_type ==  BITSET_VAR) evt = bitset_domain->set_domain(v);
  else if(domain_type ==    LIST_VAR) evt = list_domain->set_domain(v);
  else if(domain_type ==   RANGE_VAR) evt = range_domain->set_domain(v);
  else if(domain_type == VIRTUAL_VAR) evt = virtual_domain->set_domain(v);
  else if(domain_type ==   CONST_VAR) evt = (constant_value != v ? FAIL_EVENT : NO_EVENT);
  else if(domain_type ==   EXPRESSION) evt = expression->get_self().set_domain(v);
  else {

    evt = (v<0||v>1 ? FAIL_EVENT : setValue(1+v));
    //int dom = *bool_domain;
    //evt = ((dom==1+v) ? NO_EVENT : setValue(1+v));
  }

#ifdef _PROFILING_PRIMITIVE
  PROFILING_FOOT(_m_set_domain_value_)
#endif

    return evt;
}

Mistral::Event Mistral::Variable::set_min(const int lo) {

#ifdef _PROFILING_PRIMITIVE
  PROFILING_HEAD
#endif
    Event evt = NO_EVENT;

  if     (domain_type ==  BITSET_VAR) evt = bitset_domain->set_min(lo);
  else if(domain_type ==    LIST_VAR) evt = list_domain->set_min(lo);
  else if(domain_type ==   RANGE_VAR) evt = range_domain->set_min(lo);
  else if(domain_type == VIRTUAL_VAR) evt = virtual_domain->set_min(lo);
  else if(domain_type ==   CONST_VAR) evt = (constant_value < lo ? FAIL_EVENT : NO_EVENT);
  else if(domain_type ==   EXPRESSION) evt = expression->get_self().set_min(lo);
  else {

    // 1 [1, inf[
    // 2 [2, inf[
    // 3 [2, inf[
    evt = (lo==1 ? setValue(2) : (lo>1 ? FAIL_EVENT : NO_EVENT));

    //       int dom = *bool_domain;
    //       evt = (lo<=(!(dom&1)) ? NO_EVENT : (lo>(dom>>1) ? FAIL_EVENT : setValue(2)));
    //       //evt = (lo<1 ? NO_EVENT : (lo>1 ? FAIL_EVENT : setValue(2)));
  }

#ifdef _PROFILING_PRIMITIVE
  PROFILING_FOOT(_m_set_min_)
#endif

    return evt;
}

Mistral::Event Mistral::Variable::set_max(const int up) {

#ifdef _PROFILING_PRIMITIVE
  PROFILING_HEAD
#endif
    Event evt = NO_EVENT;

  if     (domain_type ==  BITSET_VAR) evt = bitset_domain->set_max(up);
  else if(domain_type ==    LIST_VAR) evt = list_domain->set_max(up);
  else if(domain_type ==   RANGE_VAR) evt = range_domain->set_max(up);
  else if(domain_type == VIRTUAL_VAR) evt = virtual_domain->set_max(up);
  else if(domain_type ==   CONST_VAR) evt = (constant_value > up ? FAIL_EVENT : NO_EVENT);
  else if(domain_type ==   EXPRESSION) evt = expression->get_self().set_max(up);
  else {

    evt =(up==0 ? setValue(1) : (up<0 ? FAIL_EVENT : NO_EVENT));

    //       int dom = *bool_domain;
    //       evt = (up>=(dom>>1) ? NO_EVENT : (up<(dom&1) ? FAIL_EVENT : setValue(1)));
    //evt = (up>0 ? NO_EVENT : (up<0 ? FAIL_EVENT : setValue(1)));
  }

#ifdef _PROFILING_PRIMITIVE
  PROFILING_FOOT(_m_set_max_)
#endif

    return evt;
}

Mistral::Event Mistral::Variable::set_domain(const BitSet& s) {

#ifdef _PROFILING_PRIMITIVE
  PROFILING_HEAD
#endif
    Event evt = NO_EVENT;

  if     (domain_type ==  BITSET_VAR) evt = bitset_domain->set_domain(s);
  else if(domain_type ==    LIST_VAR) evt = list_domain->set_domain(s);
  else if(domain_type ==   RANGE_VAR) evt = range_domain->set_domain(s);
  else if(domain_type == VIRTUAL_VAR) evt = virtual_domain->set_domain(s);
  else if(domain_type ==   CONST_VAR) evt = (s.contain(constant_value) ? NO_EVENT : FAIL_EVENT);
  else if(domain_type ==   EXPRESSION) evt = expression->get_self().set_domain(s);
  else {
    //evt = ((s.pos_words<1 || s.neg_words>0) ? FAIL_EVENT : ((s.table[0]&3)==3 ? NO_EVENT : setValue(s.table[0])));
    evt = ((s.pos_words<1 || s.neg_words>0) ? FAIL_EVENT : setState(s.table[0]&*bool_domain));
  }

#ifdef _PROFILING_PRIMITIVE
  PROFILING_FOOT(_m_set_domain_set_)
#endif

    return evt;
}

// Mistral::Event Mistral::Variable::set_domain2(Mistral::Variable& x) {

//   Event evt = NO_EVENT;
    
//   if(x.is_ground()) {
//     evt =  set_domain(x.get_min());
//   }
//   else if(x.is_range()) {
//     //evt = (set_min(x.get_min()) | set_max(x.get_max()));
//     Event evt_min = set_min(x.get_min());
//     Event evt_max = set_max(x.get_max());
      
//     evt = (evt_min | evt_max);

//     std::cout << event2str(evt_min) << " | " << event2str(evt_max) << " = " << event2str(evt);
//   }
//   else if(x.domain_type ==  BITSET_VAR) evt = set_domain(x.bitset_domain->domain.values);
//   else if(x.domain_type ==  EXPRESSION) {
//     Variable y = x.expression->get_self();
//     evt = set_domain(y);
//   } else {
//     std::cout << "TODO! (set_domain(var))" << std::endl;
//     exit(1);
//   }

//   return evt;

//   //evt = NO_EVENT;
// }

Mistral::Event Mistral::Variable::set_domain(Mistral::Variable& x) {

  Event evt = NO_EVENT;
    
  if(x.is_ground()) {
    evt =  set_domain(x.get_first());
  }
  else if(x.is_range()) {
    // //evt = (set_min(x.get_min()) | set_max(x.get_max()));
    // Event evt_min = set_min(x.get_min());
    // Event evt_max = set_max(x.get_max());
      
    // evt = (evt_min | evt_max);

    evt = set_domain(x.get_min(), x.get_max());

    //std::cout << event2str(evt_min) << " | " << event2str(evt_max) << " = " << event2str(evt);
  }
  else if(x.domain_type ==  BITSET_VAR) evt = set_domain(x.bitset_domain->domain.values);
  else if(x.domain_type ==  EXPRESSION) {
    Variable y = x.expression->get_self();
    evt = set_domain(y);
  } else {
    std::cout << "TODO! (set_domain(var))" << std::endl;
    exit(1);
  }

  return evt;

  //evt = NO_EVENT;
}


Mistral::Event Mistral::Variable::set_domain(const int lo, const int up) {
  Event evt = set_min(lo);
  if(!(FAILED(evt))) 
    evt |= set_max(up);
  return evt;
}

Mistral::Event Mistral::Variable::set_domain(const Interval& I) {
  return set_domain(I.min, I.max);
}

Mistral::Event Mistral::Variable::set_domain(const BiInterval& I) {
  Event evt = NO_EVENT;

  //std::cout << "set domain to " << I << std::endl;

  if( I.negative.empty() ) {
    if( I.positive.empty() ) {
      if( I.zero ) {
	//std::cout << "set domain to {0}\n";
	evt |= set_domain(0);
      } else evt |= FAIL_EVENT;
    } else {
      evt |= set_max(I.positive.max);
      if( !(FAILED(evt)) ) {
	if( I.zero ) {
	  //std::cout << "set min to 0\n";
	  evt |= set_min(0);
	  if( !(FAILED(evt)) && I.positive.min > 1 ) {
	    //std::cout << "remove interval [1," << (I.positive.min-1) << "]\n";
	    evt |= remove_interval(1,I.positive.min-1);
	  }
	} else {
	  //std::cout << "set min to " << I.positive.min << "\n";
	  evt |= set_min(I.positive.min);
	}
      }
    }
  } else {
    //std::cout << "set min to" << I.negative.min << "\n";
    evt |= set_min(I.negative.min);

    //std::cout << " => "<< get_domain() << std::endl;


    if( !(FAILED(evt)) ) {
      if( I.positive.empty() ) {
	if( I.zero ) {
	  //std::cout << "set max to 0\n";
	  evt |= set_max(0);
	  //std::cout << " => "<< get_domain() << std::endl;

	  if( !(FAILED(evt)) && I.negative.max < -1 ) {
	    //std::cout << "remove interval [" << (I.negative.max+1) << ",-1]\n";
	    evt |= remove_interval(I.negative.max+1,-1);

	    //std::cout << " => "<< get_domain() << std::endl;

	  }
	} else {
	  //std::cout << "set max to " << I.negative.max << "\n";
	  evt |= set_max(I.negative.max);

	  //std::cout << " => "<< get_domain() << std::endl;

	}
      } else {
	//std::cout << "set max to " << I.positive.max << "\n";
	evt |= set_max(I.positive.max);

	//std::cout << " => "<< get_domain() << std::endl;

	if( !(FAILED(evt)) ) {
	  if( I.zero ) {
	    if( I.negative.max < -1 ) {
	      //std::cout << "remove interval [" << (I.negative.max+1) << ",-1]\n";
	      evt |= remove_interval(I.negative.max+1,-1);
	    
	      //std::cout << " => "<< get_domain() << std::endl;

	    }
	    if( !(FAILED(evt)) && I.positive.min > 1 ) {
	      //std::cout << "remove interval [1," << (I.positive.min-1) << "]\n";
	      evt |= remove_interval(1,I.positive.min-1);

	      //std::cout << " => "<< get_domain() << std::endl;

	    }
	  } else {
	    //std::cout << "remove interval [" << (I.negative.max+1) << "," << (I.positive.min-1) << "]\n";
	    evt |= remove_interval(I.negative.max+1,I.positive.min-1);

	    //std::cout << " => "<< get_domain() << std::endl;

	  }
	}
      }
    }
  }

  return evt;
}


Mistral::Event Mistral::Variable::remove_set(const BitSet& s) {

#ifdef _PROFILING_PRIMITIVE
  PROFILING_HEAD
#endif
    Event evt = NO_EVENT;

  if     (domain_type ==  BITSET_VAR) evt = bitset_domain->remove_set(s);
  else if(domain_type ==    LIST_VAR) evt = list_domain->remove_set(s);
  else if(domain_type ==   RANGE_VAR) evt = range_domain->remove_set(s);
  else if(domain_type == VIRTUAL_VAR) evt = virtual_domain->remove_set(s);
  else if(domain_type ==   CONST_VAR) evt = (s.contain(constant_value) ? FAIL_EVENT : NO_EVENT);
  else if(domain_type ==   EXPRESSION) evt = expression->get_self().remove_set(s);
  else {

    evt = ((s.pos_words<1 || s.neg_words>0 || (s.table[0]^3)==3) ? NO_EVENT : setValue(s.table[0]^3));
  }

#ifdef _PROFILING_PRIMITIVE
  PROFILING_FOOT(_m_remove_set_)
#endif

    return evt;
}

Mistral::Event Mistral::Variable::remove_interval(const int lo, const int up) {

#ifdef _PROFILING_PRIMITIVE
  PROFILING_HEAD
#endif
    Event evt = NO_EVENT;

  if     (domain_type ==  BITSET_VAR) evt = bitset_domain->remove_interval(lo, up);
  else if(domain_type ==    LIST_VAR) evt = list_domain->remove_interval(lo, up);
  else if(domain_type ==   RANGE_VAR) evt = range_domain->remove_interval(lo, up);
  else if(domain_type == VIRTUAL_VAR) evt = virtual_domain->remove_interval(lo, up);
  else if(domain_type ==   CONST_VAR) evt = ((constant_value < lo || constant_value > up) ? NO_EVENT : FAIL_EVENT);
  else if(domain_type ==   EXPRESSION) evt = expression->get_self().remove_interval(lo, up);
  else {

    evt = (lo==1 ? setValue(1) : (up==0 ? setValue(2) : ((lo>1 || up<0) ? NO_EVENT : FAIL_EVENT)));
  }

#ifdef _PROFILING_PRIMITIVE
  PROFILING_FOOT(_m_remove_interval_)
#endif

    return evt;
}




Mistral::Event Mistral::Variable::restore() {

#ifdef _PROFILING_PRIMITIVE
  PROFILING_HEAD
#endif
    Event evt = NO_EVENT;

  // if(id() == 17) {
  //   for(int i=0; i<((VariableImplementation*)(variable))->solver->level; ++i)
  // 	std::cout << " ";
  //   std::cout << "restore " ;
  //   display(std::cout);
  //   std::cout << " in " << get_domain() << " " 
  // 		<< ((VariableBitmap*)(variable))->trail_ << std::endl;
  // }

  if     (domain_type ==  BITSET_VAR) evt = bitset_domain->restore();
  //else if(domain_type ==    LIST_VAR) evt = list_domain->restore();
  else if(domain_type ==   RANGE_VAR) evt = range_domain->restore();
  //else if(domain_type == VIRTUAL_VAR) evt = virtual_domain->restore();
  else if(domain_type ==   CONST_VAR) evt = NO_EVENT;
  else if(domain_type ==   EXPRESSION) {
    //std::cout << "RESTORE EXPRESSION" << std::endl;
    //exit(1);
    evt = expression->get_self().restore();
  }
  else {
    *bool_domain = 3;
    //evt = NO_EVENT;
  } 

#ifdef _PROFILING_PRIMITIVE
  PROFILING_FOOT(_m_restore_)
#endif

    return evt;
}


Mistral::Event Mistral::VariableRange::remove_interval(const int lo, const int up) {
  Event removal = DOMAIN_EVENT;

  if(lo <= min) removal = set_min(up+1);
  else if(up >= max) removal = set_max(lo-1);
  else {
    ((Solver*)solver)->make_non_convex(id);
    removal = ((Solver*)solver)->variables[id].remove_interval(lo, up);
  }

  return removal;
}


Mistral::Event Mistral::VariableRange::remove(const int v) {
  Event removal = DOMAIN_EVENT;


  //std::cout << "remove " << v << " from " << this << " in [" << min << ".." << max << "]" << std::endl;

  
  // first check if we can abort early
  if(min>v || max<v) {
    return NO_EVENT;
  }
  if(min!=v && max!=v) {
    ((Solver*)solver)->make_non_convex(id);
    removal = ((Solver*)solver)->variables[id].remove(v);
    return removal;
    //return NO_EVENT;
  }
  if(min==max) return FAIL_EVENT;
  
  save();
  
  if(min==v) {
    ++min;
    removal |= LB_EVENT;
  } else {
    --max;
    removal |= UB_EVENT;
  }
  
  if(min == max) removal |= VALUE_C; 
  solver->trigger_event(id, removal);

  //std::cout << removal << std::endl;



  return removal; 
}


/// Remove all values that do not appear in the set "s"
Mistral::Event Mistral::VariableRange::set_domain(const BitSet& s) {
  Event setdomain = NO_EVENT;

  // std::cout << "include " << s.includes(min, max) << std::endl;
  // std::cout << "intersect " << s.intersect(min, max) << std::endl;
  // std::cout << "[" << s.next(min-1) << ".." << s.prev(max+1) << "]" << std::endl;
  //std::cout << s << std::endl;

  if(s.includes(min, max)) return NO_EVENT;
  if(!s.intersect(min, max)) return FAIL_EVENT;
  int lb = s.next(min-1);
  int ub = s.prev(max+1);
  if(s.includes(lb, ub)) {
    if(lb>min) {
      setdomain |= set_min(lb);
    } 
    if(ub<max) {
      setdomain |= set_max(ub);
    }
  } else {

    // std::cout << "the intersection is not convex" << std::endl;

    ((Solver*)solver)->make_non_convex(id);


    // std::cout << solver->variables[id] << " in " 
    // 	  << solver->variables[id].get_domain() << std::endl;

    return ((Solver*)solver)->variables[id].set_domain(s);
  }

  //return set_domain(s.next(min-1), s.prev(max+1));
  return setdomain;
}


// bool Mistral::VariableImplementation::is_new(Solver *s) {
//   return (solver != s);
// }

// void Mistral::VariableImplementation::initialise(Solver *s) {
//   solver = s;
//   id = s->declare(*this);
// }

int Mistral::VariableImplementation::assigned_at_last_level() const {
  return ((Solver*)solver)->assignment_level[id] == solver->level;
}

void Mistral::VariableList::initialise(Solver *s) {
  VariableImplementation::initialise(s);
  domain.initialise(s);
}

Mistral::Event Mistral::VariableList::remove(const int v) {
  Event removal = NO_EVENT;
  
  // first check if we can abort early
  if(contain(v)) {
    if(domain.size == 1) removal = FAIL_EVENT;
    else {
      removal = DOMAIN_EVENT;
      domain.reversible_remove(v);
      if(domain.size == 1) {
	removal |= VALUE_C; 
	if(domain.head() > v) {
	  removal |= LB_EVENT;
	} else {
	  removal |= UB_EVENT;
	}
      }
      solver->trigger_event(id, removal);
    }
  }  

  //std::cout << "remove " << v << ": " << event2str(removal) << std::endl;

  return removal;
}

/// Remove all values but "v"
Mistral::Event Mistral::VariableList::set_domain(const int v) {
  Event setdomain = VALUE_EVENT;
  
  //std::cout << "set " << domain << " to " << v << " " << contain(v) << " " << is_ground() << std::endl;


  // first check if we can abort early
  if(!contain(v)) setdomain = FAIL_EVENT;
  else if(is_ground()) setdomain = NO_EVENT;
  else {
    domain.reversible_set_to(v);
    solver->trigger_event(id, setdomain);	
  }

  // std::cout << "  ==> " << domain << std::endl;

  // exit(1);
  return setdomain; 
}


void Mistral::VariableImplementation::initialise(Solver *s) {
  id = s->variables.size;
  solver = s;
}

int Mistral::VariableImplementation::get_solution_int_value() const { 
  return ((Solver*)solver)->last_solution_lb[id] ;
}  
std::string Mistral::VariableImplementation::get_solution_str_value() const { 
  std::ostringstream ret_str;
  ret_str <<  ((Solver*)solver)->last_solution_lb[id] ;
  return ret_str.str();
}  
int Mistral::VariableImplementation::get_solution_min() const { 
  return ((Solver*)solver)->last_solution_lb[id] ; 
} 
int Mistral::VariableImplementation::get_solution_max() const { 
  return ((Solver*)solver)->last_solution_ub[id] ; 
}  


// void Mistral::VariableImplementation::trigger_value_event() {
//   solver->trigger_event(id, VALUE_EVENT);
//   solver->save(id);
// }

void Mistral::VariableImplementation::trigger_event_and_save(const Event evt) {
  solver->trigger_event(id, evt);
  solver->save(id);
}

Mistral::BitsetDomain::BitsetDomain(const int lb, const int ub) {
  min = lb;
  max = ub;
  size = ub-lb+1;
}

void Mistral::BitsetDomain::initialise(const int lb, const int ub, const bool vals) {
  min = lb;
  max = ub;
  size = ub-lb+1;
  if(vals) values.initialise(lb, ub, BitSet::full);
}

Mistral::BitsetDomain::BitsetDomain(const Vector< int >& vals) {
  initialise(vals);
}

Mistral::BitsetDomain::BitsetDomain(const int lb, const int ub, const Vector< int >& vals) {
  //initialise(lb, ub, false);
  initialise(lb, ub, vals);
}

void Mistral::BitsetDomain::initialise(const int lb, const int ub, const Vector< int >& vals) {

  //std::cout << "HERE!! " << lb << " " << ub << " " << vals << " " << vals.size << std::endl;
  min = lb;
  max = ub;
  size = vals.size;
  values.initialise(lb, ub, vals);


  //std::cout << values << std::endl;


  //max = vals.back();
  //values.initialise(min, max, vals);
}

void Mistral::BitsetDomain::initialise(const Vector< int >& vals) {

  size = vals.size;

  min = vals.front();
  max = vals.front();
  
  for(int i=1; i<size; ++i) {
    if(vals[i] < min) min = vals[i];
    if(vals[i] > max) max = vals[i];
  }

  values.initialise(min, max, vals);
  
  //max = vals.back();
  //values.initialise(min, max, vals);
}

bool Mistral::Decision::propagateRelation() {
  //return !((Constraint*)(var.implementation))->propagate();
  return !var.constraint->propagate();
}

std::ostream& Mistral::operator<< (std::ostream& os, const Mistral::Goal& x) {
  return x.display(os);
}

std::ostream& Mistral::operator<< (std::ostream& os, const Mistral::Goal* x) {
  return x->display(os);
}

std::ostream& Mistral::operator<< (std::ostream& os, const Mistral::Variable& x) {
  return x.display(os);
}

std::ostream& Mistral::operator<< (std::ostream& os, const Mistral::Variable* x) {
  return x->display(os);
}

std::ostream& Mistral::operator<< (std::ostream& os, const Mistral::Decision& x) {
  return x.display(os);
}

std::ostream& Mistral::operator<< (std::ostream& os, const Mistral::Decision* x) {
  return x->display(os);
}

std::ostream& Mistral::operator<< (std::ostream& os, const Mistral::BitsetDomain& x) {
  return x.display(os);
}

std::ostream& Mistral::operator<< (std::ostream& os, const Mistral::BitsetDomain* x) {
  return x->display(os);
}

std::ostream& Mistral::operator<< (std::ostream& os, const Mistral::Expression& x) {
  return x.display(os);
}

std::ostream& Mistral::operator<< (std::ostream& os, const Mistral::Expression* x) {
  return x->display(os);
}

Mistral::Expression::Expression(const int lo, const int up) 
  : VariableImplementation() {
  id=-1; 
  Variable x(lo, up, DYN_VAR);
  _self = x;
}

Mistral::Expression::Expression(const Vector< int >& values) 
  : VariableImplementation() {
  id=-1; 
  Variable x(values, DYN_VAR);
  _self = x;
}

Mistral::Expression::Expression(const int lo, const int up, const Vector< int >& values) 
  : VariableImplementation() {
  id=-1; 
  Variable x(lo, up, values, DYN_VAR);
  _self = x;
}

Mistral::Expression::Expression(const Vector< Variable >& args) 
  : VariableImplementation() {
  id=-1; 
  for(unsigned int i=0; i<args.size; ++i)
    children.add(args[i]);
}
Mistral::Expression::Expression(const Vector< Variable >& args, const int lo, const int up) 
  : VariableImplementation() {
  id=-1; 
  Variable x(lo, up, DYN_VAR);
  _self = x;
  for(unsigned int i=0; i<args.size; ++i)
    children.add(args[i]);
}
Mistral::Expression::Expression(const std::vector< Variable >& args) 
  : VariableImplementation() {
  id=-1; 
  for(unsigned int i=0; i<args.size(); ++i)
    children.add(args[i]);
}
Mistral::Expression::Expression(const std::vector< Variable >& args, const int lo, const int up) 
  : VariableImplementation() {
  id=-1; 
  Variable x(lo, up, DYN_VAR);
  _self = x;
  for(unsigned int i=0; i<args.size(); ++i)
    children.add(args[i]);
}
Mistral::Expression::Expression(const Variable X, const Variable Y) 
  : VariableImplementation() {
  children.add(X);
  children.add(Y);
}
Mistral::Expression::Expression(Variable X) 
  : VariableImplementation() {
  children.add(X);
}
Mistral::Expression::~Expression() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete expression" << std::endl;
#endif
  // std::cout << "delete exp subvar: " << self << std::endl;

  // int domain_type = self.domain_type;
  // if     (domain_type ==  BITSET_VAR) delete self.bitset_domain;
  // else if(domain_type ==    LIST_VAR) delete self.list_domain;
  // else if(domain_type ==   RANGE_VAR) delete self.range_domain;
  // else if(domain_type == VIRTUAL_VAR) delete self.virtual_domain;
  // else if(domain_type ==  EXPRESSION) delete self.expression;
  // else if(domain_type !=   CONST_VAR) delete self.variable;
}

void Mistral::Expression::extract_variable(Solver *s) {
  // SELF CHANGE

  _self.initialise(s, 1);
  _self = _self.get_var();

  //solver = s;
  //Variable self(lb, ub, );
}


void Mistral::Expression::extract_constraint(Solver *s) {
  std::cerr << "Error: " << get_name() << " predicate can't be used as a constraint" << std::endl;
  exit(0);
}

// void Mistral::Expression::extract_predicate(Solver *s) {
//   std::cerr << "Error: " << get_name() << " constraint can't be used as a predicate" << std::endl;
//   exit(0);
// }


// void Mistral::Expression::extract_variable(Solver *s) {
//   self.initialise(s, 1);
//   self = self.get_var();
// }

Mistral::Variable Mistral::Expression::get_self() { return (id >= 0 ? ((Solver*)solver)->variables[id] : _self); }

std::ostream& Mistral::Expression::display(std::ostream& os) const {
  if(is_initialised())
    os << "e" << id << ":";
  os << get_name() << "(" ;
  os.flush();
  if(children.empty()) {
    // os << id;
    // os.flush();

    os << (id>=0 ? ((Solver*)solver)->variables[id] : _self);
    //os.flush();
  } else {
    os << children[0];
    //os.flush();
    for(unsigned int i=1; i<children.size-is_initialised(); ++i) {
      os << ", " << children[i];
      //os.flush();
    }
  }
  os << ")";
  return os;
}
//   Mistral::BinaryExpression::BinaryExpression(Variable X, Variable Y) 
//     : Expression() {
//     children.add(X);
//     children.add(Y);
//   }
//   Mistral::BinaryExpression::~BinaryExpression() {}


Mistral::AddExpression::AddExpression(Variable X, Variable Y) 
  : Expression(X,Y) {
}
Mistral::AddExpression::~AddExpression() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete add expression" << std::endl;
#endif
}

void Mistral::AddExpression::extract_constraint(Solver *s) {
  std::cerr << "Error: Add predicate can't be used as a constraint" << std::endl;
  exit(0);
}

void Mistral::AddExpression::extract_variable(Solver *s) {
  //void Mistral::AddExpression::reify(Solver *s, Variable X) {
  //int lb = children[0].get_var().get_min()+children[1].get_var().get_min();
  //int ub = children[0].get_var().get_max()+children[1].get_var().get_max();

  int lb = children[0].get_min()+children[1].get_min();
  int ub = children[0].get_max()+children[1].get_max();

  Variable aux(lb, ub, DYN_VAR);
  _self = aux;

  _self.initialise(s, 1);
  _self = _self.get_var();
  children.add(_self);
}

const char* Mistral::AddExpression::get_name() const {
  return "add";
}

void Mistral::AddExpression::extract_predicate(Solver *s) {
  s->add(Constraint(new PredicateAdd(children)));
}


Mistral::IdExpression::IdExpression(Variable X) 
  : Expression(X) {}

Mistral::IdExpression::~IdExpression() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete id expression" << std::endl;
#endif
}
  
// void Mistral::IdExpression::extract_constraint(Solver *s) {
//   std::cerr << "Error: Id predicate can't be used as a constraint" << std::endl;
//   exit(0);
// }

void Mistral::IdExpression::extract_variable(Solver *s) {
  // //void Mistral::IdExpression::reify(Solver *s, Variable X) {
  //   int lb = children[0].get_min()+id;
  //   int ub = children[0].get_max()+id;

  //     Variable aux(lb, ub, DYN_VAR);
  //     _self = aux;

  //   _self.initialise(s, 1);
  _self = children[0]; //_self.get_var();
  //   children.add(_self);
}

const char* Mistral::IdExpression::get_name() const {
  return "id";
}

// void Mistral::IdExpression::extract_predicate(Solver *s) {
//   s->add(Constraint(new PredicateId(children)));
// }


Mistral::OffsetExpression::OffsetExpression(Variable X, const int ofs) 
  : Expression(X) { 
  offset=ofs; 
}
Mistral::OffsetExpression::~OffsetExpression() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete offset expression" << std::endl;
#endif
}
  
void Mistral::OffsetExpression::extract_constraint(Solver *s) {
  std::cerr << "Error: Offset predicate can't be used as a constraint" << std::endl;
  exit(0);
}

void Mistral::OffsetExpression::extract_variable(Solver *s) {
  //void Mistral::OffsetExpression::reify(Solver *s, Variable X) {
  int lb = children[0].get_min()+offset;
  int ub = children[0].get_max()+offset;

  Variable aux(lb, ub, DYN_VAR);
  _self = aux;

  _self.initialise(s, 1);
  _self = _self.get_var();
  children.add(_self);
}

const char* Mistral::OffsetExpression::get_name() const {
  return "offset";
}

void Mistral::OffsetExpression::extract_predicate(Solver *s) {
  s->add(Constraint(new PredicateOffset(children, offset)));
}

Mistral::Variable Mistral::Variable::operator+(Variable x) {
  Variable exp(new AddExpression(*this,x));
  return exp;
}

Mistral::Variable Mistral::Variable::operator+(int k) {
  Variable exp(new OffsetExpression(*this,k));
  return exp;
}

Mistral::Variable Mistral::Variable::operator-(int k) {
  Variable exp(new OffsetExpression(*this,-k));
  return exp;
}

Mistral::MulExpression::MulExpression(Variable X, Variable Y) 
  : Expression(X,Y) {
}
Mistral::MulExpression::~MulExpression() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete mul expression" << std::endl;
#endif
}

void Mistral::MulExpression::extract_constraint(Solver *s) {
  std::cerr << "Error: Mul predicate can't be used as a constraint" << std::endl;
  exit(0);
}

void Mistral::MulExpression::extract_variable(Solver *s) {

  int z[4], lb = INFTY, ub = -INFTY, i = 4, j, k = 0, b[4];


  // std::cout << "extract variable for " << children[0].get_domain() << " x "
  //  	    << children[1].get_domain() << std::endl;



  while(i-->2) { // std::cout << i-2 << std::endl;
    b[i] = children[i-2].get_min(); }
  do { // std::cout << i << std::endl;
    b[i] = children[i].get_max(); } while(i--);

  //0.1
  //0.3
  

  //std::cout << b[0] << " " << b[1] << " " << b[2] << " " << b[3] << std::endl; 

  for(i=0; i<4; i+=2) 
    for(j=1; j<4; j+=2) {
      z[k] = b[i]*b[j];
      if(((b[i]>0 && b[j]>0) || (b[i]<0 && b[j]<0)) && z[k]<0) // int overflow
	z[k] = INFTY;
      else if(((b[i]>0 && b[j]<0) || (b[i]<0 && b[j]>0)) && z[k]>0) // int overflow
	z[k] = -INFTY;


      //std::cout << i << "." << j << " " << b[i] << " * " << b[j] << std::endl;

      ++k;
    }
  
  i = 4;
  while( i-- ) {
    if( z[i] > ub ) ub = z[i];
    if( z[i] < lb ) lb = z[i];
  }

  // std::cout << children[0].get_domain() << " x "
  //  	    << children[1].get_domain() << " = ["
  //  	    << lb << ".." << ub << "]" << std::endl;


  // exit(1);

  Variable aux(lb, ub, DYN_VAR);
  _self = aux;

  _self.initialise(s, 1);
  _self = _self.get_var();
  children.add(_self);
}

const char* Mistral::MulExpression::get_name() const {
  return "mul";
}

void Mistral::MulExpression::extract_predicate(Solver *s) {

  //std::cout << "extract predicate " << children[0] << " * " << children[1] << " = " << children[2] << std::endl;


  s->add(Constraint(new PredicateMul(children)));
}


Mistral::ModExpression::ModExpression(Variable X, Variable Y) 
  : Expression(X,Y) {
}
Mistral::ModExpression::~ModExpression() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete mod expression" << std::endl;
#endif
}

void Mistral::ModExpression::extract_constraint(Solver *s) {
  std::cerr << "Error: Mod predicate can't be used as a constraint" << std::endl;
  exit(0);
}

void Mistral::ModExpression::extract_variable(Solver *s) {
  int maxmod = children[1].get_max();
  int minmod = children[1].get_min();
  int maxabs = (abs(maxmod) > abs(minmod) ? 
		abs(maxmod) : abs(minmod) );


  Variable aux((children[0].get_min() <= 0 ? 1-maxabs : 0), 
	       (children[0].get_max() >= 0 ? maxabs-1 : 0), 
	       DYN_VAR);
  _self = aux;

  _self.initialise(s, 1);
  _self = _self.get_var();
  children.add(_self);
}

const char* Mistral::ModExpression::get_name() const {
  return "mod";
}

void Mistral::ModExpression::extract_predicate(Solver *s) {

  //std::cout << "extract predicate " << children[0] << " * " << children[1] << " = " << children[2] << std::endl;


  s->add(Constraint(new PredicateCMod(children)));
}



Mistral::ModConstantExpression::ModConstantExpression(Variable X, const int mod) 
  : Expression(X) { 
  modulo=mod; 
}
Mistral::ModConstantExpression::~ModConstantExpression() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete modconstant expression" << std::endl;
#endif
}
  
void Mistral::ModConstantExpression::extract_constraint(Solver *s) {
  std::cerr << "Error: ModConstant predicate can't be used as a constraint" << std::endl;
  exit(0);
}

void Mistral::ModConstantExpression::extract_variable(Solver *s) {
  Interval I(children[0].get_min(), children[0].get_max());
  Interval J = I%modulo;

  Variable aux(J.min, J.max, DYN_VAR);
  _self = aux;
  
  _self.initialise(s, 1);
  _self = _self.get_var();
  children.add(_self);
}

const char* Mistral::ModConstantExpression::get_name() const {
  return "modulo";
}

void Mistral::ModConstantExpression::extract_predicate(Solver *s) {
  s->add(Constraint(new PredicateCModConstant(children, modulo)));
}

Mistral::Variable Mistral::Variable::operator%(Variable x) {
  Variable exp(new ModExpression(*this,x));
  return exp;
}

Mistral::Variable Mistral::Variable::operator%(int k) {
  Variable exp(new ModConstantExpression(*this,k));
  return exp;
}


Mistral::AbsExpression::AbsExpression(Variable X) 
  : Expression(X) { }
Mistral::AbsExpression::~AbsExpression() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete abs expression" << std::endl;
#endif
}
  
void Mistral::AbsExpression::extract_constraint(Solver *s) {
  std::cerr << "Error: Abs predicate can't be used as a constraint" << std::endl;
  exit(0);
}

void Mistral::AbsExpression::extract_variable(Solver *s) {
  int lb = 0;
  if(children[0].get_max() < 0) lb = -children[0].get_max();
  if(children[0].get_min() > 0) lb =  children[0].get_min();

  Variable aux(lb, std::max(children[0].get_max(), -children[0].get_min()), DYN_VAR);
  _self = aux;
  
  _self.initialise(s, 1);
  _self = _self.get_var();
  children.add(_self);
}

const char* Mistral::AbsExpression::get_name() const {
  return "abs";
}

void Mistral::AbsExpression::extract_predicate(Solver *s) {
  s->add(Constraint(new PredicateAbs(children)));
}

Mistral::Variable Mistral::Abs(Variable x) {
  Variable exp(new AbsExpression(x));
  return exp;
}



Mistral::DivExpression::DivExpression(Variable X, Variable Y) 
  : Expression(X,Y) {
}
Mistral::DivExpression::~DivExpression() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete div expression" << std::endl;
#endif
}

void Mistral::DivExpression::extract_constraint(Solver *s) {
  std::cerr << "Error: Div predicate can't be used as a constraint" << std::endl;
  exit(0);
}

void Mistral::DivExpression::extract_variable(Solver *s) {
  // z = x/y

  BiInterval X(children[0]);
  BiInterval Y(children[1]);

  BiInterval Z = X/Y;

  Variable aux(Z.get_min(), Z.get_max(), DYN_VAR);
  _self = aux;

  _self.initialise(s, 1);
  _self = _self.get_var();
  children.add(_self);
}

const char* Mistral::DivExpression::get_name() const {
  return "div";
}

void Mistral::DivExpression::extract_predicate(Solver *s) {
  s->add(Constraint(new PredicateDiv(children)));
}


Mistral::QuotientExpression::QuotientExpression(Variable X, const int quo) 
  : Expression(X) { 
  quotient=quo; 
}
Mistral::QuotientExpression::~QuotientExpression() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete quotient expression" << std::endl;
#endif
}
  
void Mistral::QuotientExpression::extract_constraint(Solver *s) {
  std::cerr << "Error: Quotient predicate can't be used as a constraint" << std::endl;
  exit(0);
}

void Mistral::QuotientExpression::extract_variable(Solver *s) {

  BiInterval X(children[0]);
  BiInterval Y = X/quotient;

  Variable aux(Y.get_min(), Y.get_max(), DYN_VAR);
  _self = aux;
  
  _self.initialise(s, 1);
  _self = _self.get_var();
  children.add(_self);
}

const char* Mistral::QuotientExpression::get_name() const {
  return "quotient";
}

void Mistral::QuotientExpression::extract_predicate(Solver *s) {
  s->add(Constraint(new PredicateDivConstant(children, quotient)));
}

Mistral::Variable Mistral::Variable::operator/(Variable x) {
  Variable exp(new DivExpression(*this,x));
  return exp;
}

Mistral::Variable Mistral::Variable::operator/(int k) {
  Variable exp(new QuotientExpression(*this,k));
  return exp;
}




Mistral::FactorExpression::FactorExpression(Variable X, const int fct) 
  : Expression(X) { 
  factor=fct; 
}
Mistral::FactorExpression::~FactorExpression() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete factor expression" << std::endl;
#endif
}
  
void Mistral::FactorExpression::extract_constraint(Solver *s) {
  std::cerr << "Error: Factor predicate can't be used as a constraint" << std::endl;
  exit(0);
}

void Mistral::FactorExpression::extract_variable(Solver *s) {

  int lb = (factor<0 ? children[0].get_max() : children[0].get_min())*factor;
  int ub = (factor<0 ? children[0].get_min() : children[0].get_max())*factor;

  Variable aux(lb, ub, DYN_VAR);
  _self = aux;
  
  _self.initialise(s, 1);
  _self = _self.get_var();
  children.add(_self);
}

const char* Mistral::FactorExpression::get_name() const {
  return "factor";
}

void Mistral::FactorExpression::extract_predicate(Solver *s) {
  s->add(Constraint(new PredicateFactor(children, factor)));
}



Mistral::SquareExpression::SquareExpression(Variable X) 
  : Expression(X) { 
}
Mistral::SquareExpression::~SquareExpression() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete square expression" << std::endl;
#endif
}
  
void Mistral::SquareExpression::extract_constraint(Solver *s) {
  std::cerr << "Error: Square predicate can't be used as a constraint" << std::endl;
  exit(0);
}

void Mistral::SquareExpression::extract_variable(Solver *s) {

	int lb = 0;
	int ub = 0;
	int minabs = 0;
	if(children[0].get_max() > -children[0].get_min()) {
		ub = children[0].get_max()*children[0].get_max();
		minabs = -children[0].get_min();
	} else {
		ub = children[0].get_min()*children[0].get_min();
		minabs = children[0].get_max();
	}
	while(lb<minabs) {
		if(children[0].contain(lb)) {
			break;
		} else ++lb;
	}
	lb *= lb;

  Variable aux(lb, ub, DYN_VAR);
  _self = aux;
  
  _self.initialise(s, 1);
  _self = _self.get_var();
  children.add(_self);
}

const char* Mistral::SquareExpression::get_name() const {
  return "square";
}

void Mistral::SquareExpression::extract_predicate(Solver *s) {
  s->add(Constraint(new PredicateSquare(children)));
}




Mistral::Variable Mistral::Variable::operator*(Variable x) {

  // std::cout << "ADD " ;
  // display(std::cout);
  // std::cout << " * " << x << std::endl;
	Variable exp;
	if(this->same_as(x)) {
		exp = Variable(new SquareExpression(x));
	} else {
		exp = Variable(new MulExpression(*this,x));
	}
  return exp;
}

Mistral::Variable Mistral::Variable::operator*(int k) {
  Variable exp(new FactorExpression(*this,k));
  return exp;
}


// Mistral::Variable Mistral::Variable::operator/(Variable x) {
//   Variable exp(new DivExpression(*this,x));
//   return exp;
// }

// Mistral::Variable Mistral::Variable::operator/(int k) {
//   Variable exp(new QuotientExpression(*this,k));
//   return exp;
// }


Mistral::SubExpression::SubExpression(Variable X, Variable Y) 
  : Expression(X,Y) {
}
Mistral::SubExpression::~SubExpression() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete sub expression" << std::endl;
#endif
}
  
void Mistral::SubExpression::extract_constraint(Solver *s) {
  std::cerr << "Error: Sub predicate can't be used as a constraint" << std::endl;
  exit(0);
}

void Mistral::SubExpression::extract_variable(Solver *s) {
  //void Mistral::SubExpression::reify(Solver *s, Variable X) {
  //int lb = children[0].get_var().get_min()-children[1].get_var().get_max();
  //int ub = children[0].get_var().get_max()-children[1].get_var().get_min();
  
  int lb = children[0].get_min()-children[1].get_max();
  int ub = children[0].get_max()-children[1].get_min();
  
  Variable aux(lb, ub, EXPRESSION);
  _self = aux;

  _self.initialise(s, 1);
  //_self = _self.get_var();
  children.add(_self);
 
}

void Mistral::SubExpression::extract_predicate(Solver *s) {
  VarArray tmp;
  for(int i=3; i;) tmp.add(children[--i]);


  //std::cout << "children: " << tmp << std::endl;

  Constraint sub(new PredicateAdd(tmp));
    
  //sub->initialise();
  //return sub;
  s->add(sub);
}

const char* Mistral::SubExpression::get_name() const {
  return "sub";
}

Mistral::Variable Mistral::Variable::operator-(Variable x) {
  Variable exp(new SubExpression(*this,x));
  return exp;
}


Mistral::NotExpression::NotExpression(Variable X) 
  : Expression(X) { 
}
Mistral::NotExpression::~NotExpression() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete not expression" << std::endl;
#endif
}
  
void Mistral::NotExpression::extract_constraint(Solver *s) {
  children[0].remove(0);
  // std::cerr << "Error: Not predicate can't be used as a constraint" << std::endl;
  // exit(0);
}

void Mistral::NotExpression::extract_variable(Solver *s) {
  Variable aux(0, 1, BOOL_VAR);
  _self = aux;

  _self.initialise(s, 1);
  _self = _self.get_var();
  children.add(_self);
}

const char* Mistral::NotExpression::get_name() const {
  return "not";
}

void Mistral::NotExpression::extract_predicate(Solver *s) {
  s->add(Constraint(new PredicateNot(children)));
}

Mistral::Variable Mistral::Variable::operator!() {
  Variable exp(new NotExpression(*this));
  return exp;
}


Mistral::NegExpression::NegExpression(Variable X) 
  : Expression(X) { 
}
Mistral::NegExpression::~NegExpression() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete neg expression" << std::endl;
#endif
}
  
void Mistral::NegExpression::extract_constraint(Solver *s) {
  //children[0].remove(0);
  std::cerr << "Error: Neg predicate can't be used as a constraint" << std::endl;
  // exit(0);
}

void Mistral::NegExpression::extract_variable(Solver *s) {
  Variable aux(-children[0].get_max(), -children[0].get_min(), DYN_VAR);
  _self = aux;

  _self.initialise(s, 1);
  _self = _self.get_var();
  children.add(_self);
}

const char* Mistral::NegExpression::get_name() const {
  return "neg";
}

void Mistral::NegExpression::extract_predicate(Solver *s) {
  s->add(Constraint(new PredicateNeg(children)));
}

Mistral::Variable Mistral::Variable::operator-() {
  Variable exp(new NegExpression(*this));
  return exp;
}

Mistral::AndExpression::AndExpression(Variable X, Variable Y, const int sp) 
  : Expression(X,Y) { spin=sp; }
Mistral::AndExpression::~AndExpression() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete and expression" << std::endl;
#endif
}
  
void Mistral::AndExpression::extract_constraint(Solver *s) {
  //s->add(new ConstraintAnd(children));
  //std::cout << "not implemented" << std::endl;
  //exit(1);
  if(spin) {

#ifdef _DEBUG_AC
    std::cout << "pre-propagte (" << children[0] << " AND "
	      << children[1] << ") ";
#endif

    if(FAILED(children[0].remove(0))) {
#ifdef _DEBUG_AC
      std::cout << "fail!\n";
#endif
      s->fail();
    }
    else if(FAILED(children[1].remove(0))) {
#ifdef _DEBUG_AC
      std::cout << "fail!\n";
#endif
      s->fail();
    }
#ifdef _DEBUG_AC
    else
      std::cout << "ok\n";
#endif



  } else {
    s->add(Constraint(new ConstraintNotAnd(children)));
  }
}


void Mistral::AndExpression::extract_variable(Solver *s) {
  Variable aux(0, 1, BOOL_VAR);
  _self = aux;

  _self.initialise(s, 1);
  _self = _self.get_var();
  children.add(_self);
}

void Mistral::AndExpression::extract_predicate(Solver *s) {
  if(spin) {
    s->add(Constraint(new PredicateAnd(children)));
  } else {
    //s->add(Constraint(new PredicateNotAnd(children)));
  }
}

const char* Mistral::AndExpression::get_name() const {
  return "and";
}

Mistral::Variable Mistral::NotAnd(Variable x, Variable y) {
  Variable exp(new AndExpression(x,y,false));
  return exp;
}

Mistral::Variable Mistral::Variable::operator&&(Variable x) {
  Variable exp(new AndExpression(*this,x));
  return exp;
}



Mistral::OrExpression::OrExpression(Variable X, Variable Y, const int sp) 
  : Expression(X,Y) { spin=sp; }
Mistral::OrExpression::~OrExpression() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete or expression" << std::endl;
#endif
}
  
void Mistral::OrExpression::extract_constraint(Solver *s) {
  if(spin) {
    s->add(Constraint(new ConstraintOr(children)));
  } else {
    
#ifdef _DEBUG_AC
    std::cout << "pre-propagte NOT(" << children[0] << " OR "
	      << children[1] << ") ";
#endif
    
    if(FAILED(children[0].set_domain(0))) {
#ifdef _DEBUG_AC
      std::cout << "fail!\n";
#endif
      s->fail();
    } else if(FAILED(children[1].set_domain(0))) {
#ifdef _DEBUG_AC
      std::cout << "fail!\n";
#endif
      s->fail();
    }
#ifdef _DEBUG_AC
    else
      std::cout << "ok\n";
#endif
    //s->add(ConstraintOr::ConstraintOr_new(children));
    //std::cout << "not implemented" << std::endl;
    //exit(1);
  }
}
  
void Mistral::OrExpression::extract_variable(Solver *s) {
  Variable aux(0, 1, BOOL_VAR);
  _self = aux;
  
  _self.initialise(s, 1);
  _self = _self.get_var();
  children.add(_self);
}

void Mistral::OrExpression::extract_predicate(Solver *s) {
  if(spin) {
    s->add(Constraint(new PredicateOr(children)));
  } else {
    //s->add(Constraint(new PredicateNotOr(children)));
  }
}

const char* Mistral::OrExpression::get_name() const {
  return "or";
}

Mistral::Variable Mistral::NotOr(Variable x, Variable y) {
  Variable exp(new OrExpression(x,y,false));
  return exp;
}


Mistral::Variable Mistral::Variable::operator||(Variable x) {
  Variable exp(new OrExpression(*this,x));
  return exp;
}




// //   Mistral::NeqExpression::NeqExpression(Variable X, Variable Y) 
// //   : BinaryExpression(X,Y) {}
// //   Mistral::NeqExpression::~NeqExpression() {}
  
// // void Mistral::NeqExpression::extract_constraint(Solver *s) {
// // //     Constraint *neq = new ConstraintNotEqual(children);
// // //     neq->initialise();
// // //     return neq;
// //   s->add(new ConstraintNotEqual(children));
// //   }

// //   void Mistral::NeqExpression::extract_variable(Solver *s) {
// //     Variable aux(0, 1, BOOL_VAR);
// //     _self = aux;

// //     children.add(_self);
// //     _self.initialise(s, 1);
// //   }

// //   void Mistral::NeqExpression::extract_predicate(Solver *s) {
// // //     Constraint *neq = new PredicateEqual(children, false);
// // //     neq->initialise();
// // //     return neq;
// //     s->add(new PredicateEqual(children, false));
// //   }

// // const char* Mistral::NeqExpression::get_name() const {
// //   return "neq";
// // }



Mistral::EqualExpression::EqualExpression(Variable X, Variable Y, const int sp) 
  : Expression(X,Y) { spin=sp; value=NOVAL; }
Mistral::EqualExpression::EqualExpression(Variable X, const int y, const int sp) 
  : Expression() { children.add(X); value=y; spin=sp; }
Mistral::EqualExpression::~EqualExpression() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete equal expression" << std::endl;
#endif
}

void Mistral::EqualExpression::extract_constraint(Solver *s) {

  //std::cout << "extract constraint from expression equal" << std::endl;
  
  if(spin) {
    if(children.size==2)
      s->add(Constraint(new ConstraintEqual(children), (BINARY|IDEMPOTENT)));
    else {

#ifdef _DEBUG_AC
      std::cout << "pre-propagte (" << children[0] << " = "
		<< value << ") ";
#endif

      if(children[0].set_domain(value) == FAIL_EVENT) {
	s->fail();
#ifdef _DEBUG_AC
	std::cout << "fail!\n";
      } else {
	std::cout << "ok\n";
#endif
      }
    }

  } else {
    if(children.size==2) {
      s->add(Constraint(new ConstraintNotEqual(children)));
    } else {

#ifdef _DEBUG_AC
      std::cout << "pre-propagte (" << children[0] << " != "
		<< value << ") ";
#endif
      
      if(children[0].remove(value) == FAIL_EVENT) {

	s->fail();
#ifdef _DEBUG_AC
	std::cout << "fail!\n";
      } else {
	std::cout << "ok\n";
#endif
      }
    }
  }
}


// void Mistral::EqualExpression::extract_constraint(Solver *s) {
//   if(spin) {
    
//     if(children.size==2) s->add(Constraint(new ConstraintEqual(children)// , (BINARY|IDEMPOTENT)
// 					   ));
//     else children[0].set_domain(value);
    
//   } else {
    
//     if(children.size==2) {
//       s->add(Constraint(new ConstraintNotEqual(children)// , (BINARY|IDEMPOTENT)
// 			));
//     } else children[0].remove(value);
//   }
// }

void Mistral::EqualExpression::extract_variable(Solver *s) {
  Variable aux(0, 1, BOOL_VAR);
  _self = aux;
  
  _self.initialise(s, 1);
  _self = _self.get_var();
  children.add(_self);
}

void Mistral::EqualExpression::extract_predicate(Solver *s) {
  
  if(children.size==3) {

    s->add(Constraint(new PredicateEqual(children, spin)));

  } else s->add(Constraint(new PredicateConstantEqual(children, value, spin)));
}

const char* Mistral::EqualExpression::get_name() const {
  return "equal";
}

// Mistral::EqualSetExpression::EqualSetExpression(Variable X, Variable Y, const int sp) 
//   : Expression(X,Y) { spin=sp; value=NOVAL; }
// Mistral::EqualSetExpression::EqualSetExpression(Variable X, const int y, const int sp) 
//   : Expression() { children.add(X); value=y; spin=sp; }
// Mistral::EqualSetExpression::~EqualSetExpression() {
// #ifdef _DEBUG_MEMORY
//   std::cout << "c delete set expression" << std::endl;
// #endif
// }

// void Mistral::EqualSetExpression::extract_constraint(Solver *s) {
//   if(spin) {

//     unsigned int i = 0, j = 0;
//     SetExpression *x = (SetExpression*)(children[0].expression);
//     SetExpression *y = (SetExpression*)(children[1].expression);
//     Variable bx;
//     Variable by;
    
//     s->add(new ConstraintEqual(children));
    
//     Vector<Variable> scp;
    
//     //std::cout << x->elts_ub << " == " << y->elts_ub << std::endl;
    
    
//     while(i<x->elts_ub.size && j<y->elts_ub.size) {
      
//       // std::cout << i << " - " << j << std::endl;
      
//       // std::cout << x->get_index_var(i) << std::endl;
      
//       // std::cout << y->get_index_var(j) << std::endl;
      
//       scp.clear();
//       bx = x->get_index_var(i);
//       by = y->get_index_var(j);
//       if(x->elts_ub[i] == y->elts_ub[j]) {
// 	if(bx.is_ground()) {
// 	  by.set_domain(bx.get_min());
// 	} else if(by.is_ground()) {
// 	  bx.set_domain(by.get_min());
// 	} else {
// 	  scp.add(by);
// 	  scp.add(bx);
// 	  s->add(new ConstraintEqual(scp));
// 	}
// 	++i;
// 	++j;
//       } else if(x->elts_ub[i] <= y->elts_ub[j]) {
// 	// x->elts_ub[i] can't be in Y
// 	bx.set_domain(0);
// 	++i;
//       } else {
// 	// y->elts_ub[j] can't be in X
// 	by.set_domain(0);
// 	++j;
//       }
//     }
    
//     // remaining values of X
//     while(i<x->elts_ub.size) {
//       bx = x->get_index_var(i);
//       bx.set_domain(0);
//       ++i;
//     }
    
//     // remaining values of Y
//     while(j<y->elts_ub.size) {
//       by = y->get_index_var(j);
//       by.set_domain(0);
//       ++j;
//     }
//   } else {

//     if(children[0].domain_type == EXPRESSION && children[0].expression->is_set()
//        && children.size>1 &&
//        children[1].domain_type == EXPRESSION && children[1].expression->is_set()
//        ) {
//       unsigned int i = 0, j = 0;
//       SetExpression *x = (SetExpression*)(children[0].expression);
//       SetExpression *y = (SetExpression*)(children[1].expression);
//       Variable bx;
//       Variable by;

//       Vector<Variable> scp;
//       Vector<Variable> disjunction;

//       while(i<x->elts_ub.size && j<y->elts_ub.size) {
// 	Variable b(0,1);

// 	scp.clear();
// 	bx = x->get_index_var(i);
// 	by = y->get_index_var(j);
// 	if(x->elts_ub[i] == y->elts_ub[j]) {
// 	  if(bx.is_ground()) {
// 	    scp.add(by);
// 	    scp.add(b);
// 	    s->add(Constraint(new PredicateConstantEqual(scp,bx.get_min())));
// 	    disjunction.add(b);
// 	  } else if(by.is_ground()) {
// 	    scp.add(bx);
// 	    scp.add(b);
// 	    s->add(Constraint(new PredicateConstantEqual(scp,by.get_min())));
// 	    disjunction.add(b);
// 	  } else {
// 	    scp.add(by);
// 	    scp.add(bx);
// 	    scp.add(b);
// 	    s->add(Constraint(new PredicateConstantEqual(scp)));
// 	    disjunction.add(b);
// 	  }
// 	  ++i;
// 	  ++j;
// 	} else if(x->elts_ub[i] <= y->elts_ub[j]) {
// 	  // x->elts_ub[i] can't be in Y
// 	  scp.add(bx);
// 	  scp.add(b);
// 	  s->add(Constraint(new PredicateConstantEqual(scp,0)));
// 	  disjunction.add(b);
// 	  ++i;
// 	} else {
// 	  // y->elts_ub[j] can't be in X
// 	  scp.add(by);
// 	  scp.add(b);
// 	  s->add(Constraint(new PredicateConstantEqual(scp,0)));
// 	  disjunction.add(b);
// 	  ++j;
// 	}
//       }
    
//       // remaining values of X
//       while(i<x->elts_ub.size) {
// 	Variable b(0,1);
// 	bx = x->get_index_var(i);
// 	scp.add(bx);
// 	scp.add(b);
// 	s->add(Constraint(new PredicateConstantEqual(scp,0)));
// 	disjunction.add(b);
// 	++i;
//       }
      
//       // remaining values of Y
//       while(j<y->elts_ub.size) {
// 	Variable b(0,1);
// 	by = y->get_index_var(j);
// 	scp.add(by);
// 	scp.add(b);
// 	s->add(Constraint(new PredicateConstantEqual(scp,0)));
// 	disjunction.add(b);
// 	++j;
//       }

//       s->add(Constraint(new ConstraintBoolSumInterval(disjunction,0,disjunction.size-1)));
//     }
//   }
// }

// void Mistral::EqualSetExpression::extract_variable(Solver *s) {
//   Variable aux(0, 1, BOOL_VAR);
//   _self = aux;
  
//   _self.initialise(s, 1);
//   _self = _self.get_var();
//   children.add(_self);
// }

// void Mistral::EqualSetExpression::extract_predicate(Solver *s) {
//   if(children.size==3) {
//     if(spin) {

//       Vector<Variable> scp;



//       unsigned int i = 0, j = 0;
//       int x, y;

//       SetExpression *X = (SetExpression*)(children[0].expression);
//       SetExpression *Y = (SetExpression*)(children[1].expression);
//       Variable bx;
//       Variable by;

//       // std::cout << X->elts_lb << " <= " << children[0] << " <= " << X->elts_ub << std::endl;
//       // std::cout << Y->elts_lb << " <= " << children[1] << " <= " << Y->elts_ub << std::endl;



//       // std::cout << "to be equal, the cardinalities must be equal" << std::endl;

//       Variable equal_card(0,1);
//       for(int i=0; i<2; ++i) scp.add(children[i].get_var());
//       scp.add(equal_card);
//       s->add(Constraint(new PredicateEqual(scp, 1)));

//       scp.clear();
//       scp.add(children[2]);
//       scp.add(equal_card);
//       s->add(Constraint(new ConstraintLess(scp)));

//       //std::cout << s << std::endl;
      
      

      
//       Vector<Variable> aux_pos;
//       Vector<Variable> aux_neg;

      
//       while(i<X->elts_ub.size || j<Y->elts_ub.size) {
// 	// scp.clear();
// 	// Variable b(0,1);
// 	// b.initialise(s, 1);
// 	// aux.add(b);
// 	// scp.add(children[2]);
// 	// scp.add(aux.back());
// 	// s->add(new ConstraintLess(scp));
// 	scp.clear();

// 	x = (i < X->elts_ub.size ? X->elts_ub[i] : INFTY);
// 	y = (j < Y->elts_ub.size ? Y->elts_ub[j] : INFTY);

// 	if(x == y) {

// 	  bx = X->get_index_var(i);
// 	  by = Y->get_index_var(j);

// 	  // std::cout << children[2] << " => (" << x << " in X) iff (" 
// 	  // 	    << y << " in Y)" << std::endl;

// 	  if(bx.is_ground()) {
// 	    scp.add(children[2]);
// 	    scp.add(by);
// 	    s->add(Constraint(new ConstraintLess(scp)));
// 	    //scp.add(aux.back());

// 	    //std::cout << "(" << x << " in X) is true " << std::endl;
	   	    
// 	    //s->add(new PredicateConstantEqual(scp, bx.get_min(), 1));

// 	    aux_pos.add(by);

// 	  } else if(by.is_ground()) {
// 	    scp.add(children[2]);
// 	    scp.add(bx);
// 	    s->add(Constraint(new ConstraintLess(scp)));
// 	    //scp.add(aux.back());

// 	    //std::cout << "(" << y << " in X) is true " << std::endl;

// 	    //s->add(new PredicateConstantEqual(scp, by.get_min(), 1));

// 	    aux_pos.add(bx);

// 	  } else {
// 	    Variable b(0,1);
// 	    b.initialise(s, 1);
// 	    aux_pos.add(b);
// 	    scp.add(children[2]);
// 	    scp.add(aux_pos.back());
// 	    s->add(Constraint(new ConstraintLess(scp)));

// 	    scp.clear();

// 	    scp.add(bx);
// 	    scp.add(by);
// 	    scp.add(aux_pos.back());
// 	    s->add(Constraint(new PredicateEqual(scp, 1)));
// 	  }

// 	  ++i;
// 	  ++j;
// 	} else if(x < y) {
// 	  // x->elts_ub[i] can't be in Y

// 	  //std::cout << "if (" << x << " in X) then X != Y " << std::endl;

// 	  bx = X->get_index_var(i);
// 	  scp.add(bx);
// 	  scp.add(children[2]);
// 	  s->add(Constraint(new ConstraintNotAnd(scp)));

// 	  aux_neg.add(bx);
// 	  //scp.add(aux.back());
// 	  //s->add(new PredicateConstantEqual(scp, 0, 1));
// 	  ++i;
// 	} else {

// 	  //std::cout << "if (" << y << " in Y) then X != Y " << std::endl;

// 	  // y->elts_ub[j] can't be in X
// 	  by = Y->get_index_var(j);
// 	  scp.add(by);
// 	  scp.add(children[2]);
// 	  s->add(Constraint(new ConstraintNotAnd(scp)));

// 	  aux_neg.add(by);
// 	  // scp.add(aux.back());
// 	  // s->add(new PredicateConstantEqual(scp, 0, 1));
// 	  ++j;
// 	}
// 	//std::cout << s << std::endl;
//       }
      

//       if(aux_neg.size) {
// 	if(aux_pos.size) {
// 	  s->add(((BoolSum(aux_pos) == aux_pos.size) && (BoolSum(aux_neg) == 0)) <= children[2]);
// 	} else {
// 	  s->add((BoolSum(aux_neg) == 0) <= children[2]);
// 	}
//       } else if(aux_pos.size) {
// 	s->add((BoolSum(aux_pos) == aux_pos.size) <= children[2]);
//       }


//       // // remaining values of X
//       // while(i<x->elts_ub.size) {
	
//       // 	scp.add(by);
//       // 	scp.add(aux.back());
//       // 	s->add(new PredicateConstantEqual(scp, 0, 1));
//       // 	++i;
//       // }
      
//       // // remaining values of Y
//       // while(j<y->elts_ub.size) {
	
//       // 	scp.add(bx);
//       // 	scp.add(aux.back());
//       // 	s->add(new PredicateConstantEqual(scp, 0, 1));
//       // 	++j;
//       // }
//     } else {
//       std::cerr << "TO DO" << std::endl;
//       exit(1);
//     }
//   } else {
//     std::cerr << "TO DO" << std::endl;
//     exit(1);
//   }
// }

// const char* Mistral::EqualSetExpression::get_name() const {
//   return "set_equal";
// }


Mistral::EqualSetExpression::EqualSetExpression(Variable X, Variable Y, const int sp) 
  : Expression(X,Y) { spin=sp; value=NOVAL; }
Mistral::EqualSetExpression::EqualSetExpression(Variable X, const int y, const int sp) 
  : Expression() { children.add(X); value=y; spin=sp; }
Mistral::EqualSetExpression::~EqualSetExpression() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete set expression" << std::endl;
#endif
}

void Mistral::EqualSetExpression::extract_constraint(Solver *s) {

  //std::cout << "EXTRACT an " << (spin ? "equality" : "inequality") << " constraint" << std::endl;


  SetExpression *X = (SetExpression*)(children[0].variable);
  SetExpression *Y = (SetExpression*)(children[1].variable);


  // X->display(std::cout, true);
  // std::cout <<  std::endl;
  // Y->display(std::cout, true);
  // std::cout <<  std::endl;

  // index in X (var and lb, resp.)
  unsigned int i=0, il=0;
  // index in Y (var and lb, resp.)
  unsigned int j=0, jl=0;

  // x,y: value (in X and Y, resp.)
  int x, y, xl, yl;

  Vector<Variable> disjunction_neg;
  Vector<Variable> disjunction_pos;
  bool satisfied = false;
  int lb_inter = 0;

  // the variable part is the intersection of the upper bounds minus the intersection of the lower bounds
  while(i<X->elts_var.size || il<X->elts_lb.size || 
	j<Y->elts_var.size || jl<Y->elts_lb.size) {

    x = (i<X->elts_var.size ? X->elts_var[i] : INFTY);
    y = (j<Y->elts_var.size ? Y->elts_var[j] : INFTY);
    xl = (il<X->elts_lb.size ? X->elts_lb[il] : INFTY);
    yl = (jl<Y->elts_lb.size ? Y->elts_lb[jl] : INFTY);

    //std::cout << "x=" << x << ", xl=" << xl << ", y=" << y << ", yl=" << yl << std::endl;


    if(xl<x && yl<y) {
      // both in the lower bound, if y is not in lb(x) it must be a fail
      if(yl == xl) {
	// common value in the lower bound
	++lb_inter;
	++il; ++jl;
      } else if(xl<yl) {
	// xl is in the lower bound of x, but cannot be in y
	if(spin) {

	  //std::cout << "cant't be equal since " << xl << " can't be in Y" << std::endl; 

	  s->fail();
	} else {
	  satisfied = true;
	}
	break;
	//++il;
      } else if(xl>yl) {
	// yl is in the lower bound of y, but cannot be in x
	if(spin) {

	  //std::cout << "cant't be equal since " << yl << " can't be in X" << std::endl; 

	  s->fail();
	} else {
	  satisfied = true;
	}
	break;
	//++jl;
      }
    } else {
      if(xl<x) { //(therefore y<yl) 
	if(xl == y) {
	  // common value in the lower bound of x and the var part of y
	  if(spin) {

	    //std::cout << y << " must be in Y" << std::endl;

	    s->add(Y->get_index_var(j) == 1);
	  } else {
	    disjunction_neg.add(Y->get_index_var(j));
	  }
	  ++il; ++j;
	} else if(xl<y) {
	  // xl is in the lower bound of x, but cannot be in y
	  if(spin) {

	    //std::cout << "cant't be equal since " << xl << " can't be in Y" << std::endl; 

	    s->fail();
	  } else {
	    satisfied = true;
	  }
	  break;
	  //++il;
	} else {
	  // y is in the var part of y, but cannot be in x
	  if(spin) {

	    //std::cout << y << " can't be in Y" << std::endl;

	    s->add(Y->get_index_var(j) == 0);
	  } else {
	    disjunction_pos.add(Y->get_index_var(j));
	  }
	  ++j;
	}
      } else if(yl<y) { //(therefore x<xl)
	if(yl == x) {
	  // common value in the lower bound of y and the var part of x
	  if(spin) {

	    //std::cout << x << " must be in X" << std::endl;

	    s->add(X->get_index_var(i) == 1);
	  } else {
	    disjunction_neg.add(X->get_index_var(i));
	  }
	  ++jl; ++i;
	} else if(yl<x) {
	  // yl is in the lower bound of y, but cannot be in x
	  if(spin) {

	    //std::cout << "cant't be equal since " << yl << " can't be in X" << std::endl; 	    

	    s->fail();
	  } else {
	    satisfied = true;
	  }
	  break;
	  //++jl;
	} else {
	  // x is in the var part of x, but cannot be in y
	  if(spin) {

	    //std::cout << x << " must be in X" << std::endl;

	    s->add(X->get_index_var(i) == 0);
	  } else {
	    disjunction_pos.add(X->get_index_var(i));
	  }
	  ++i;
	}
      } else {
	if(x == y) {
	  // common value in the var part
	  if(spin) {

	    //std::cout << x << " is in X iff " << y << " is in Y" << std::endl;

	    s->add(Constraint(new ConstraintEqual(X->get_index_var(i), Y->get_index_var(j))));
	  } else {
	    if(disjunction_pos.size < disjunction_neg.size) {
	      disjunction_neg.add(X->get_index_var(i) == Y->get_index_var(j));
	    } else {
	      disjunction_pos.add(X->get_index_var(i) != Y->get_index_var(j));
	    }
	  }
	  ++i; ++j;
	} else if(x<y) {
	  // x is in the var part of x, but cannot be in y
	  if(spin) {

	    //std::cout << x << " can't be in X" << std::endl;

	    s->add(X->get_index_var(i) == 0);
	  } else {
	    disjunction_pos.add(X->get_index_var(i));
	  }
	  ++i;
	}
	else {
	  // y is in the var part of y, but cannot be in x
	  if(spin) {

	    //std::cout << y << " can't be in Y" << std::endl;

	    s->add(Y->get_index_var(j) == 0);
	  } else {
	    disjunction_pos.add(Y->get_index_var(j));
	  }
	  ++j;
	}
      }
    }
  }
 
  if(spin) {
    int x_offset = X->elts_lb.size-lb_inter;
    int y_offset = Y->elts_lb.size-lb_inter;

    Variable cx = (x_offset ? Variable(new OffsetExpression(children[0],x_offset)) : children[0]);
    Variable cy = (y_offset ? Variable(new OffsetExpression(children[1],y_offset)) : children[1]);

    s->add(Variable(new EqualExpression(cx, cy)));

    //std::cout << "X and Y must have the same cardinality" << std::endl;

    //s->add(Card(children[0]) == Card(children[1]));
  } else {
    if(!satisfied) {
      if(disjunction_neg.size) {
	if(disjunction_pos.size) {
	  s->add((BoolSum(disjunction_neg) < disjunction_neg.size) || 
		 (BoolSum(disjunction_pos) > 0));
	} else {
	  s->add(BoolSum(disjunction_neg) < disjunction_neg.size);
	}
      } else if(disjunction_pos.size) {
	s->add(BoolSum(disjunction_pos) > 0);
      } else {
	s->fail();
      }
    }   
  }
}

void Mistral::EqualSetExpression::extract_variable(Solver *s) {
  Variable aux(0, 1, BOOL_VAR);
  _self = aux;
  
  _self.initialise(s, 1);
  _self = _self.get_var();
  children.add(_self);
}

void Mistral::EqualSetExpression::extract_predicate(Solver *s) {
  SetExpression *X = (SetExpression*)(children[0].variable);
  SetExpression *Y = (SetExpression*)(children[1].variable);
  
  // index in X (var and lb, resp.)
  unsigned int i=0, il=0;
  // index in Y (var and lb, resp.)
  unsigned int j=0, jl=0;
  
  // x,y: value (in X and Y, resp.)
  int x, y, xl, yl;
  
  Vector<Variable> disjunction_neg;
  Vector<Variable> disjunction_pos;
  bool satisfied = false;
  bool violated = false;
  
  // the variable part is the intersection of the upper bounds minus the intersection of the lower bounds
  while(i<X->elts_var.size || il<X->elts_lb.size || 
	j<Y->elts_var.size || jl<Y->elts_lb.size) {
    
    x = (i<X->elts_var.size ? X->elts_var[i] : INFTY);
    y = (j<Y->elts_var.size ? Y->elts_var[j] : INFTY);
    xl = (il<X->elts_lb.size ? X->elts_lb[il] : INFTY);
    yl = (jl<Y->elts_lb.size ? Y->elts_lb[jl] : INFTY);

    //std::cout << "x=" << x << ", xl=" << xl << ", y=" << y << ", yl=" << yl << std::endl;


    if(xl<x && yl<y) {
      // both in the lower bound, if y is not in lb(x) it must be a fail
      if(yl == xl) {
	// common value in the lower bound
	++il; ++jl;
      } else if(xl<yl) {
	// xl is in the lower bound of x, but cannot be in y
	if(spin) {
	  violated = true;
	} else {
	  satisfied = true;
	}
	++il;
      } else if(xl>yl) {
	// yl is in the lower bound of y, but cannot be in x
	if(spin) {
	  violated = true;
	} else {
	  satisfied = true;
	}
	++jl;
      }
    } else {
      if(xl<x) { //(therefore y<yl) 
	if(xl == y) {
	  // common value in the lower bound of x and the var part of y
	  //if(spin) {
	  disjunction_pos.add(Y->get_index_var(j));
	  // } else {
	  //   disjunction_neg.add(Y->get_index_var(j));
	  // }
	  ++il; ++j;
	} else if(xl<y) {
	  // xl is in the lower bound of x, but cannot be in y
	  if(spin) {
	    violated = true;
	  } else {
	    satisfied = true;
	  }
	  ++il;
	} else {
	  // y is in the var part of y, but cannot be in x
	  //if(spin) {
	  disjunction_neg.add(Y->get_index_var(j));
	  // } else {
	  //   disjunction_pos.add(Y->get_index_var(j));
	  // }
	  ++j;
	}
      } else if(yl<y) { //(therefore x<xl)
	if(yl == x) {
	  // common value in the lower bound of y and the var part of x
	  //if(spin) {
	  disjunction_pos.add(X->get_index_var(i));
	  // } else {
	  //   disjunction_neg.add(X->get_index_var(i));
	  // }
	  ++jl; ++i;
	} else if(yl<x) {
	  // yl is in the lower bound of y, but cannot be in x
	  if(spin) {
	    violated = true;
	  } else {
	    satisfied = true;
	  }
	  ++jl;
	} else {
	  // x is in the var part of x, but cannot be in y
	  //if(spin) {
	  disjunction_neg.add(X->get_index_var(i));
	  // } else {
	  //   disjunction_pos.add(X->get_index_var(i));
	  // }
	  ++i;
	}
      } else {
	if(x == y) {
	  // common value in the var part
	  //if(spin) {
	  disjunction_pos.add((X->get_index_var(i) == Y->get_index_var(j)));
	  // } else {
	  //   disjunction_neg.add((X->get_index_var(i) == Y->get_index_var(j)));
	  // }
	  ++i; ++j;
	} else if(x<y) {
	  // x is in the var part of x, but cannot be in y
	  //if(spin) {
	  disjunction_neg.add(X->get_index_var(i));
	  // } else {
	  //   disjunction_pos.add(X->get_index_var(i));
	  // }
	  ++i;
	}
	else {
	  // y is in the var part of y, but cannot be in x
	  //if(spin) {
	  disjunction_neg.add(Y->get_index_var(j));
	  // } else {
	  //   disjunction_pos.add(Y->get_index_var(j));
	  // }
	  ++j;
	}
      }
    }
  }


  if(satisfied) {
    if(violated) {
      s->fail();
    } else {
      if(FAILED(children[2].set_domain(true))) {
	//std::cout << "this fail" << std::endl;
	s->fail();
      }
    }
  } else if(violated) {
    if(FAILED(children[2].set_domain(false))) {
      //std::cout << "that fail" << std::endl;
      s->fail();
    }
  } else {
    if(spin) {
      for(unsigned int i=0; i<disjunction_pos.size; ++i) {
	s->add(Variable(new PrecedenceExpression(children[2],disjunction_pos[i])));
      }
      for(unsigned int i=0; i<disjunction_neg.size; ++i) {
	s->add(Variable(new AndExpression(children[2], disjunction_neg[i], false)));
      }
    } else {
      for(unsigned int i=0; i<disjunction_neg.size; ++i) {
	s->add(Variable(new PrecedenceExpression(disjunction_neg[i], children[2])));
      }
      for(unsigned int i=0; i<disjunction_pos.size; ++i) {
	s->add(Variable(new OrExpression(children[2], disjunction_pos[i])));
      }
    }
    if(disjunction_neg.size) {
      if(disjunction_pos.size) {
	if(spin) {
	  s->add( ((BoolSum(disjunction_pos) == disjunction_pos.size) && 
		   (BoolSum(disjunction_neg) == 0))
		  <=
		  children[2] );
	} else {
	  s->add( children[2] 
		  <= 
		  ((BoolSum(disjunction_neg) > 0) || 
		   (BoolSum(disjunction_pos) < disjunction_pos.size)) );
	}
      } else {
	if(spin) 
	  s->add((BoolSum(disjunction_neg) == 0) <= children[2]);
	else
	  s->add( children[2] <= (BoolSum(disjunction_neg) > 0) );
      }
    } else if(disjunction_pos.size) {
      if(spin)
	s->add((BoolSum(disjunction_pos) == disjunction_pos.size) <= children[2]);
      else
	s->add( children[2] <= (BoolSum(disjunction_pos) < disjunction_pos.size) );
    } else {
      if(FAILED(children[2].set_domain(spin))) s->fail();
    }
  }   
}


const char* Mistral::EqualSetExpression::get_name() const {
  return "set_equal2";
}

Mistral::Variable Mistral::Variable::operator==(Variable x) {
  Variable exp;
  if(is_set_var() && x.is_set_var()) {
    exp = Variable(new EqualSetExpression(*this,x,1));
    // } else if(x.is_set_var()) {
    //   exp = Variable(new EqualSetExpression(x,*this,1));
  } else {
    exp = Variable(new EqualExpression(*this,x,1));
  }
  return exp;
}

Mistral::Variable Mistral::Variable::operator!=(Variable x) {
  Variable exp;
  if(is_set_var() && x.is_set_var()) {
    exp = Variable(new EqualSetExpression(*this,x,0));
    // } else if(x.is_set_var()) {
    //   exp = Variable(new EqualSetExpression(x,*this,0));
  } else {
    exp = Variable(new EqualExpression(*this,x,0));
  }
  return exp;
  // Variable exp(new EqualExpression(*this,x,0));
  // return exp;
}

Mistral::Variable Mistral::Variable::operator==(const int x) {
  Variable exp(new EqualExpression(*this,x,1));
  return exp;
}

Mistral::Variable Mistral::Variable::operator!=(const int x) {
  Variable exp(new EqualExpression(*this,x,0));
  return exp;
}

Mistral::Variable Mistral::Variable::operator!=(const int x[2]) {
  Variable exp(new MemberExpression(*this,x[0],x[1],0));
  return exp;
}

// Mistral::Variable Mistral::Variable::operator==(Variable x) {
//   Variable exp(new EqualExpression(*this,x));
//   return exp;
// }


Mistral::PrecedenceExpression::PrecedenceExpression(Variable X, Variable Y, 
						    const int of, const int sp) 
  : Expression(X,Y) { spin = sp; offset = of; }
Mistral::PrecedenceExpression::PrecedenceExpression(Variable X,  
						    const int of, const int sp) 
  : Expression(X) { spin = sp; offset = of; }
Mistral::PrecedenceExpression::~PrecedenceExpression() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete precedence expression" << std::endl;
#endif
}
  
void Mistral::PrecedenceExpression::extract_constraint(Solver *s) {
  if(children.size==2) {
    s->add(Constraint(new ConstraintLess(children, offset)// , (BINARY|IDEMPOTENT)
		      ) );
  }
  else {
    if(spin) {
      
#ifdef _DEBUG_AC
      std::cout << "pre-propagte (" << children[0] << " <= "
		<< offset << ") ";
#endif
      //std::cout << "HERE" << std::endl;

      if(FAILED(children[0].set_max(offset))) 
	{ 
#ifdef _DEBUG_AC
	  std::cout << "fail!\n";
#endif
	  s->fail(); 
	}
#ifdef _DEBUG_AC
      else
	std::cout << "ok\n";
#endif

    } else {    
#ifdef _DEBUG_AC
      std::cout << "pre-propagte (" << children[0] << " >= "
		<< offset << ") ";
#endif
 
      if(FAILED(children[0].set_min(offset)))
	{ 
#ifdef _DEBUG_AC
	  std::cout << "fail!\n";
#endif
	  s->fail(); 
	}
#ifdef _DEBUG_AC
      else
	std::cout << "ok\n";
#endif
    }
  }
}

void Mistral::PrecedenceExpression::extract_variable(Solver *s) {
  Variable aux(0, 1, BOOL_VAR);
  _self = aux;
  
  _self.initialise(s, 1);
  _self = _self.get_var();
  children.add(_self);
}

void Mistral::PrecedenceExpression::extract_predicate(Solver *s) {
  if(children.size==3) {
    s->add(Constraint(new PredicateLess(children, offset)));
  } else if(spin) {
    s->add(Constraint(new PredicateUpperBound(children, offset)));
  } else {
    s->add(Constraint(new PredicateLowerBound(children, offset)));
  }
}

const char* Mistral::PrecedenceExpression::get_name() const {
  return "prec";
}

Mistral::Variable Mistral::Precedence(Variable X, const int d, Variable Y) 
{
  Variable exp(new PrecedenceExpression(X,Y,d));
  return exp;
}


Mistral::Variable Mistral::Variable::operator<(Variable x) {
  Variable exp(new PrecedenceExpression(*this,x,1));
  return exp;
}

Mistral::Variable Mistral::Variable::operator>(Variable x) {
  Variable exp(new PrecedenceExpression(x,*this,1));
  return exp;
}

Mistral::Variable Mistral::Variable::operator<=(Variable x) {
  Variable exp(new PrecedenceExpression(*this,x));
  return exp;
}

Mistral::Variable Mistral::Variable::operator>=(Variable x) {
  Variable exp(new PrecedenceExpression(x,*this));
  return exp;
}

Mistral::Variable Mistral::Variable::operator<(const int k) {
  Variable exp(new PrecedenceExpression(*this,k-1,1));
  return exp;
}

Mistral::Variable Mistral::Variable::operator>(const int k) {
  Variable exp(new PrecedenceExpression(*this,k+1,0));
  return exp;
}

Mistral::Variable Mistral::Variable::operator<=(const int k) {
  Variable exp(new PrecedenceExpression(*this,k,1));
  return exp;
}

Mistral::Variable Mistral::Variable::operator>=(const int k) {
  Variable exp(new PrecedenceExpression(*this,k,0));
  return exp;
}


Mistral::DisjunctiveExpression::DisjunctiveExpression(const Variable X, const Variable Y, 
						      const int px, const int py) 
  : Expression(X,Y) { processing_time[0] = px; processing_time[1] = py; }
Mistral::DisjunctiveExpression::~DisjunctiveExpression() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete disjunctive expression" << std::endl;
#endif
}
  
void Mistral::DisjunctiveExpression::extract_constraint(Solver *s) {
  s->add(Constraint(new ConstraintDisjunctive(children, processing_time[0], processing_time[1])));
}

void Mistral::DisjunctiveExpression::extract_variable(Solver *s) {
  std::cerr << "Error: Disjunctive constraint can't yet be used as a predicate" << std::endl;
  exit(0);
}

void Mistral::DisjunctiveExpression::extract_predicate(Solver *s) {
  std::cerr << "Error: Disjunctive constraint can't yet be used as a predicate" << std::endl;
  exit(0);
}

const char* Mistral::DisjunctiveExpression::get_name() const {
  return "disjunct";
}

Mistral::Variable Mistral::Disjunctive(Variable X, Variable Y, 
				       const int px, const int py) 
{
  Variable exp(new DisjunctiveExpression(X,Y,px,py));
  return exp;
}



Mistral::ReifiedDisjunctiveExpression::ReifiedDisjunctiveExpression(Variable X, Variable Y, 
								    const int px, const int py) 
  : Expression(X,Y) { processing_time[0] = px; processing_time[1] = py; }

Mistral::ReifiedDisjunctiveExpression::~ReifiedDisjunctiveExpression() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete reified disjunctive expression" << std::endl;
#endif
}
  
void Mistral::ReifiedDisjunctiveExpression::extract_constraint(Solver *s) {
  std::cerr << "Error: ReifiedDisjunctive constraint can't yet be used as a constraint" << std::endl;
  exit(0);
}

void Mistral::ReifiedDisjunctiveExpression::extract_variable(Solver *s) {
  Variable aux(0, 1, BOOL_VAR);
  _self = aux;

  _self.initialise(s, 1);
  _self = _self.get_var();
  children.add(_self);
}

void Mistral::ReifiedDisjunctiveExpression::extract_predicate(Solver *s) {
  //s->add(new ConstraintTernaryDisjunctive(children, processing_time[0], processing_time[1]));
  s->add(Constraint(new ConstraintReifiedDisjunctive(children, processing_time[0], processing_time[1])));
}

const char* Mistral::ReifiedDisjunctiveExpression::get_name() const {
  return "r-disjunct";
}

Mistral::Variable Mistral::ReifiedDisjunctive(Variable X, Variable Y, 
					      const int px, const int py) 
{
  Variable exp(new ReifiedDisjunctiveExpression(X,Y,px,py));
  return exp;
}



Mistral::FreeExpression::FreeExpression(Variable X) 
  : Expression(X) { };

Mistral::FreeExpression::~FreeExpression() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete free expression" << std::endl;
#endif
}

void Mistral::FreeExpression::extract_constraint(Solver *s) {
}

  
const char* Mistral::FreeExpression::get_name() const {
  return "free";
}

Mistral::Variable Mistral::Free(Variable X) 
{
  Variable exp(new FreeExpression(X));
  return exp;
}


Mistral::AllDiffExpression::AllDiffExpression(Vector< Variable >& args, const int ct) 
  : Expression(args) { consistency_level = ct; }

Mistral::AllDiffExpression::~AllDiffExpression() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete alldiff expression" << std::endl;
#endif
}

void Mistral::AllDiffExpression::extract_constraint(Solver *s) { 
  //   Constraint *con = new ConstraintAllDiff(children); 
  //   con->initialise();
  //   return con;
  if(consistency_level == BOUND_CONSISTENCY)
    s->add(Constraint(new ConstraintAllDiff(children))); 
  s->add(Constraint(new ConstraintCliqueNotEqual(children))); 
  //   Vector< Variable > pair;
  //   for(unsigned int i=0; i<children.size-1; ++i)
  //     for(unsigned int j=i+1; j<children.size; ++j) {
  //       pair.clear();
  //       pair.add(children[i]);
  //       pair.add(children[j]);
  //       s->add(new ConstraintNotEqual(pair));
  //     }
}

void Mistral::AllDiffExpression::extract_variable(Solver *s) {
  std::cerr << "Error: AllDiff constraint can't yet be used as a predicate" << std::endl;
  exit(0);
}

void Mistral::AllDiffExpression::extract_predicate(Solver *s) { 
  std::cerr << "Error: AllDiff constraint can't yet be used as a predicate" << std::endl;
  exit(0);
}

const char* Mistral::AllDiffExpression::get_name() const {
  return "alldiff";
}

Mistral::Variable Mistral::AllDiff(Vector< Variable >& args, const int ct) {
  Variable exp(new AllDiffExpression(args,ct));
  return exp;
}



Mistral::OccurrencesExpression::OccurrencesExpression(Vector< Variable >& args, const int first, const int last, const int* lb, const int* ub, const int ct) 
  : Expression(args) { 
  firstval = first;
  lastval = last;
  lower_bounds = lb;
  upper_bounds = ub;
  consistency_level = ct; 
}

Mistral::OccurrencesExpression::~OccurrencesExpression() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete gcc expression" << std::endl;
#endif
}

void Mistral::OccurrencesExpression::extract_constraint(Solver *s) { 
  if(consistency_level == BOUND_CONSISTENCY)
    s->add(Constraint(new ConstraintOccurrences(children, firstval, lastval, lower_bounds, upper_bounds)));

  int lb, ub;
  for(int v=firstval; v<=lastval; ++v) {
    lb = lower_bounds[v-firstval];
    ub = upper_bounds[v-firstval];
    s->add(Occurrence(children, v, lb, ub));
  }
}

void Mistral::OccurrencesExpression::extract_variable(Solver *s) {
  std::cerr << "Error: Occurrences constraint can't yet be used as a predicate" << std::endl;
  exit(0);
}

void Mistral::OccurrencesExpression::extract_predicate(Solver *s) { 
  std::cerr << "Error: Occurrences constraint can't yet be used as a predicate" << std::endl;
  exit(0);
}

const char* Mistral::OccurrencesExpression::get_name() const {
  return "gcc";
}

Mistral::Variable Mistral::Occurrences(Vector< Variable >& args, const int first, const int last, const int* lb, const int* ub, const int ct) {
  Variable exp(new OccurrencesExpression(args, first, last, lb, ub, ct));
  return exp;
}



Mistral::VertexCoverExpression::VertexCoverExpression(Vector< Variable >& args, const Graph& g) 
: _G(g), Expression(args) {
	assert(_G.size() == args.size);
}
  
Mistral::VertexCoverExpression::~VertexCoverExpression() {}

void Mistral::VertexCoverExpression::extract_constraint(Solver*) {
	std::cerr << "Error: Vertex Cover predicate can't be used as a constraint" << std::endl;
	exit(0);
}

void Mistral::VertexCoverExpression::extract_variable(Solver* s) {
	
	Variable aux(0, _G.size(), DYN_VAR);
	_self = aux;
  
	_self.initialise(s, 1);
	_self = _self.get_var();
	children.add(_self);
}

void Mistral::VertexCoverExpression::extract_predicate(Solver* s) {
	s->add(Constraint(new PredicateVertexCover(children, _G)));
}

const char* Mistral::VertexCoverExpression::get_name() const {
	return "vertex cover";
}


Mistral::Variable Mistral::VertexCover(Mistral::Vector< Variable >& args, const Graph& g) {
	Variable exp(new VertexCoverExpression(args, g));
	return exp;
}



Mistral::FootruleExpression::FootruleExpression(Vector< Variable >& args) 
: Expression(args) {
}
  
Mistral::FootruleExpression::~FootruleExpression() {}

void Mistral::FootruleExpression::extract_constraint(Solver*) {
	std::cerr << "Error: Footrule predicate can't be used as a constraint" << std::endl;
	exit(0);
}

void Mistral::FootruleExpression::extract_variable(Solver* s) {
	Variable aux(0, children.size*children.size/2, DYN_VAR);
	_self = aux;
  
	_self.initialise(s, 1);
	_self = _self.get_var();
	children.add(_self);
}

void Mistral::FootruleExpression::extract_predicate(Solver* s) {
	s->add(Constraint(new PredicateFootrule(children)));
}

const char* Mistral::FootruleExpression::get_name() const {
	return "footrule";
}


Mistral::Variable Mistral::Footrule(Vector< Variable >& arg1, Vector< Variable >& arg2) {
	Vector<Variable> args(arg1);
	for(int i=0; i<arg2.size; ++i) args.add(arg2[i]);
	Variable exp(new FootruleExpression(args));
	return exp;
}




Mistral::AtMostSeqCardExpression::AtMostSeqCardExpression(Vector< Variable >& args, const int d, const int p, const int q)
  : Expression(args), _k(1), _d(d) {  
  _p = new int[1];
  _q = new int[1];
  _p[0] = p;
  _q[0] = q;
}

Mistral::AtMostSeqCardExpression::AtMostSeqCardExpression(Vector< Variable >& args, const int d, const Vector< Tuple<2, int> >& c)
  : Expression(args), _k(c.size), _d(d) {  

  _p = new int[_k];
  _q = new int[_k];

  for(int i=0; i<_k; ++i) {
    _p[i] = c[i][0];
    _q[i] = c[i][1];
  }
}

Mistral::AtMostSeqCardExpression::~AtMostSeqCardExpression() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete atmostseqcard expression" << std::endl;
#endif

  delete [] _p;
  delete [] _q;
}

void Mistral::AtMostSeqCardExpression::extract_constraint(Solver *s) { 
  s->add(Constraint(new ConstraintMultiAtMostSeqCard(children, _k, _d, _p, _q))); 
}

void Mistral::AtMostSeqCardExpression::extract_variable(Solver *s) {
  std::cerr << "Error: AtMostSeqCard constraint can't yet be used as a predicate" << std::endl;
  exit(0);
}

const char* Mistral::AtMostSeqCardExpression::get_name() const {
  return "amsc";
}

Mistral::Variable Mistral::AtMostSeqCard(Vector< Variable >& args, const int d, const int p, const int q) {
  Variable exp(new AtMostSeqCardExpression(args,d,p,q));
  return exp;
}

Mistral::Variable Mistral::MultiAtMostSeqCard(Vector< Variable >& args, const int d, const Vector< Tuple<2, int> >& c) {
  Variable exp(new AtMostSeqCardExpression(args,d,c));
  return exp;
}



Mistral::TableExpression::TableExpression(Vector< Variable >& args, Vector<const int*>& rel, const AlgorithmType ct) 
  : Expression(args) { 
  propagator = ct; 
  tuples.copy(rel);
}

Mistral::TableExpression::TableExpression(Vector< Variable >& args, const AlgorithmType ct) 
  : Expression(args) { 
  propagator = ct; 
}

void Mistral::TableExpression::add(int* tuple) 
{
  int n = children.size;
  int *new_tuple = new int[n];
  while(n--) new_tuple[n] = tuple[n];
  tuples.add(new_tuple);
}

Mistral::TableExpression::~TableExpression() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete table expression" << std::endl;
#endif

  tuples.neutralise();
}

void Mistral::TableExpression::extract_constraint(Solver *s) { 
  ConstraintTable *tab;
  switch(propagator) {
  case GAC2001: {tab = new ConstraintGAC2001(children);} break;
  case GAC3: {tab = new ConstraintGAC3(children);} break;
    //case AC3: {tab = new ConstraintAC3(children);} break;
  case GAC4: {tab = new ConstraintGAC4(children);} break;
  default: {tab = new ConstraintGAC2001(children);}
  }
 
  tab->table.copy(tuples);

  s->add(Constraint(tab)); 
}

void Mistral::TableExpression::extract_variable(Solver *s) {
  std::cerr << "Error: Table constraint can't yet be used as a predicate" << std::endl;
  exit(0);
}

void Mistral::TableExpression::extract_predicate(Solver *s) { 
  std::cerr << "Error: Table constraint can't yet be used as a predicate" << std::endl;
  exit(0);
}

const char* Mistral::TableExpression::get_name() const {
  return "alldiff";
}

Mistral::Variable Mistral::Table(Vector< Variable >& args, Vector<const int*>& rel, const Mistral::TableExpression::AlgorithmType ct) {
  Variable exp(new TableExpression(args, rel, ct));
  return exp;
}
// Mistral::Variable Mistral::Table(VarArray& args, const Mistral::TableExpression::AlgorithmType ct) {
//   Variable exp(new TableExpression(args,ct));
//   return exp;
// }


Mistral::LexExpression::LexExpression(Vector< Variable >& r1, Vector< Variable >& r2, const int st_)
  : Expression() { 
  int row_size = r1.size;
  for(int i=0; i<row_size; ++i)
    children.add(r1[i]);
  for(int i=0; i<row_size; ++i)
    children.add(r2[i]);
  for(int i=0; i<=row_size; ++i)
    children.add( Variable(0,1) );
  strict = st_; 
}

Mistral::LexExpression::~LexExpression() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete lex expression" << std::endl;
#endif
}

void Mistral::LexExpression::extract_constraint(Solver *s) { 
  int arity = (children.size-1)/3;
  VarArray scp;
  // scp.add(children[0]);
  // scp.add(children[arity]);
  // scp.add(children[2*arity]);

  // s->add( new ConstraintLexf(scp) );

  for(int i=0; i<arity; ++i) {
    scp.clear();
    
    scp.add(children[i]);
    scp.add(children[i+arity]);
    scp.add(children[i+2*arity]);
    scp.add(children[i+2*arity+1]);

    s->add(Constraint(new ConstraintLex(scp)));
  }

  if(FAILED(children[2*arity].set_domain(0)))
    { s->fail(); }
  if(strict) 
    if(FAILED(children[3*arity].set_domain(1)))
      { s->fail(); }
}

void Mistral::LexExpression::extract_variable(Solver *s) {
  std::cerr << "Error: Lex constraint can't yet be used as a predicate" << std::endl;
  exit(0);
}

void Mistral::LexExpression::extract_predicate(Solver *s) { 
  std::cerr << "Error: Lex constraint can't yet be used as a predicate" << std::endl;
  exit(0);
}

const char* Mistral::LexExpression::get_name() const {
  return "lex";
}

Mistral::Variable Mistral::LexLeq(VarArray& r1, VarArray& r2) {
  Variable exp(new LexExpression(r1, r2, false));
  return exp;
}

Mistral::Variable Mistral::LexLess(VarArray& r1, VarArray& r2) {
  Variable exp(new LexExpression(r1, r2, true));
  return exp;
}

Mistral::Variable Mistral::VarArray::operator<(VarArray& X) {
  return LexLess(*this, X);
}

Mistral::Variable Mistral::VarArray::operator>(VarArray& X) {
  return LexLess(X, *this);
}

Mistral::Variable Mistral::VarArray::operator<=(VarArray& X) {
  return LexLeq(*this, X);
}

Mistral::Variable Mistral::VarArray::operator>=(VarArray& X) {
  return LexLeq(X, *this);
}


Mistral::BoolSumExpression::BoolSumExpression(const int l, const int u) 
  : Expression() {
  lower_bound = l;
  upper_bound = u;
}

Mistral::BoolSumExpression::BoolSumExpression(Vector< Variable >& args, const int l, const int u) 
  : Expression(args) {
  lower_bound = l;
  upper_bound = u;
  //remove_duplicates_and_zeros();
}

Mistral::BoolSumExpression::BoolSumExpression(std::vector< Variable >& args, const int l, const int u) 
  : Expression(args) {
  lower_bound = l;
  upper_bound = u;
  //remove_duplicates_and_zeros();
}

Mistral::BoolSumExpression::BoolSumExpression(Vector< Variable >& args, const Vector< int >& wgts) 
  : Expression(args) {
  lower_bound = -INFTY;
  upper_bound = INFTY;
  for(unsigned int i=0; i<wgts.size; ++i) {
    weight.add(wgts[i]);
  }
  //remove_duplicates_and_zeros();
}

Mistral::BoolSumExpression::BoolSumExpression(Vector< Variable >& args, const Vector< int >& wgts, const int l, const int u) 
  : Expression(args) {
  lower_bound = l;
  upper_bound = u;
  for(unsigned int i=0; i<wgts.size; ++i) {
    weight.add(wgts[i]);
  }
  //remove_duplicates_and_zeros();
}

Mistral::BoolSumExpression::BoolSumExpression(std::vector< Variable >& args, const std::vector< int >& wgts, const int l, const int u) 
  : Expression(args) {
  lower_bound = l;
  upper_bound = u;
  for(unsigned int i=0; i<wgts.size(); ++i) {
    weight.add(wgts[i]);
  }
  //remove_duplicates_and_zeros();
}

void Mistral::BoolSumExpression::remove_duplicates_and_zeros() {
  int i,n=children.size;
  while(n--) {
    if(weight.size && weight[n]==0) {
      children.pop();
      weight.pop();
    } else {
      i = n;
      while(i--) {
	if(children[i].same_as(children[n])) {
	  if(weight.size==0) {
	    for(int j=0; j<children.size; ++j)
	      weight.add(1);
	  }
	  weight[i] += weight[n];
	  children.pop();
	  weight.pop();
	  break;
	}
      }
    }
  }
}


Mistral::BoolSumExpression::~BoolSumExpression() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete boolsum expression" << std::endl;
#endif
}

#define _INCREMENTAL_WBOOLSUM

void Mistral::BoolSumExpression::extract_constraint(Solver *s) { 
  if(weight.empty()) {
    // if(lower_bound == upper_bound) {
    //   s->add(new ConstraintBoolSumEqual(children,lower_bound)); 
    // } else {
    if(lower_bound == (int)(children.size)) {
      
#ifdef _DEBUG_AC
      std::cout << "pre-propagte sum(" << children[0] ;
      for(unsigned int i=0; i<children.size; ++i)
	std::cout << ", " << children[i];
      std::cout << ") = " << children.size;
#endif
      
      for(unsigned int i=0; i<children.size; ++i) {
	if(FAILED(children[i].set_domain(1))) {
#ifdef _DEBUG_AC
	  std::cout << " FAIL!" << std::endl;
#endif
	  s->fail();
	  break;
	}
#ifdef _DEBUG_AC
	else
	  std::cout << "ok\n";
#endif
      }
    } else if(upper_bound == 0) {
      for(unsigned int i=0; i<children.size; ++i) {
	if(FAILED(children[i].set_domain(0))) {
#ifdef _DEBUG_AC
	  std::cout << " FAIL!" << std::endl;
	  //exit(1);
#endif
	  s->fail();
	  break;
	}
#ifdef _DEBUG_AC
	else
	  std::cout << "ok\n";
#endif
      }
    } else {
      s->add(Constraint(new ConstraintBoolSumInterval(children,lower_bound,upper_bound))); 
    }
  } else {

#ifdef _INCREMENTAL_WBOOLSUM
    s->add(Constraint(new ConstraintIncrementalWeightedBoolSumInterval(children,weight,lower_bound,upper_bound)));  
#else
    s->add(Constraint(new ConstraintWeightedBoolSumInterval(children,weight,lower_bound,upper_bound))); 
#endif
  }
}


void Mistral::BoolSumExpression::initialise_bounds() {
  int tlb=0;
  int tub=0;

  if(weight.size) {
    int lb = children[0].get_min()*weight[0];
    int ub = children[0].get_max()*weight[0];
    
    if(lb < ub) {
      tlb += lb;
      tub += ub;
    } else {
      tlb += ub;
      tub += lb;
    }
    
    for(unsigned int i=1; i<children.size; ++i) {
      lb = children[i].get_min()*weight[i];
      ub = children[i].get_max()*weight[i];
      
      if(lb < ub) {
	tlb += lb;
	tub += ub;
      } else {
	tlb += ub;
	tub += lb;
      }
    }
  } else {
    tlb = 0;
    tub = children.size;
  }

  if(tlb > lower_bound) lower_bound = tlb;
  if(tub < upper_bound) upper_bound = tub;
}

void Mistral::BoolSumExpression::extract_variable(Solver *s) {
  remove_duplicates_and_zeros();
  initialise_bounds();

  Variable aux(lower_bound, upper_bound, DYN_VAR);
  _self = aux;
  
  _self.initialise(s, 1);
  _self = _self.get_var();
  if(!weight.empty()) {
    children.add(_self);
  }
  //   weight.add(-1);
  // }
}

void Mistral::BoolSumExpression::extract_predicate(Solver *s) { 
  //s->add(new PredicateBoolSum(children, self)); 

  if(weight.empty()) 
    s->add(Constraint(new PredicateBoolSum(children, s->variables[id]))); 
  else {
    s->add(Constraint(new PredicateWeightedBoolSum(children, weight)));
  }
  //s->add(new PredicateWeightedSum(children, weight, 0, 0));
}

const char* Mistral::BoolSumExpression::get_name() const {
  return "bool_sum";
}

Mistral::Variable Mistral::BoolSum(std::vector< Variable >& args, const int l, const int u) {
  Variable exp(new BoolSumExpression(args,l,(u != -INFTY ? u : l)));
  return exp;
}

Mistral::Variable Mistral::BoolSum(Vector< Variable >& args, const int l, const int u) {
  Variable exp(new BoolSumExpression(args,l,(u != -INFTY ? u : l)));
  return exp;
}

Mistral::Variable Mistral::BoolSum(Vector< Variable >& args) {
  Variable exp(new BoolSumExpression(args,0,args.size));
  return exp;
}

Mistral::Variable Mistral::BoolSum(Vector< Variable >& args, const Vector< int >& w) {
  Variable exp(new BoolSumExpression(args,w));
  return exp;
}

Mistral::Variable Mistral::BoolSum(std::vector< Variable >& args, const std::vector< int >& w, const int l, const int u) {
  Variable exp(new BoolSumExpression(args,w,l,(u != -INFTY ? u : l)));
  return exp;
}

Mistral::Variable Mistral::BoolSum(Vector< Variable >& args, const Vector< int >& w, const int l, const int u) {
  Variable exp(new BoolSumExpression(args,w,l,(u != -INFTY ? u : l)));
  return exp;
}

// Mistral::Variable Mistral::BoolSum(Vector< Variable >& args, Vector< int >& w) {
//   Variable exp(new BoolSumExpression(args,w,0,args.size));
//   return exp;
// }




Mistral::ParityExpression::ParityExpression(Vector< Variable >& args, const int p) 
  : Expression(args) {
  target_parity = p;
}

Mistral::ParityExpression::~ParityExpression() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete boolsum expression" << std::endl;
#endif
}


void Mistral::ParityExpression::extract_constraint(Solver *s) {
  VarArray scope;
  for(unsigned int i=0; i<children.size; ++i)
    if(children[i].is_ground())
      target_parity ^= children[i].get_value();
    else
      scope.add(children[i]);

  if(scope.size)
    s->add(Constraint(new ConstraintParity(scope,target_parity))); 
  else
    if(target_parity) s->fail();
}

void Mistral::ParityExpression::extract_variable(Solver *s) {
  Variable aux(0, 1, DYN_VAR);
  _self = aux;
  
  _self.initialise(s, 1);
  _self = _self.get_var();
  children.add(_self);
}

void Mistral::ParityExpression::extract_predicate(Solver *s) { 
  s->add(Constraint(new ConstraintParity(children,target_parity))); 
}

const char* Mistral::ParityExpression::get_name() const {
  return "parity";
}

Mistral::Variable Mistral::Parity(Vector< Variable >& args, const int p) {
  Variable exp(new ParityExpression(args,p));
  return exp;
}


Mistral::OrderedSumExpression::OrderedSumExpression(Vector< Variable >& args, 
					    const int l, const int u, const int o) 
  : Expression(args) {
  lower_bound = l;
  upper_bound = u;
}

Mistral::OrderedSumExpression::OrderedSumExpression(std::vector< Variable >& args, 
					    const int l, const int u, const int o) 
  : Expression(args) {
  lower_bound = l;
  upper_bound = u;
  offset = o;
}


Mistral::OrderedSumExpression::~OrderedSumExpression() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete osum expression" << std::endl;
#endif
}
  
// -1 -1 -1: - x - y =  z: 
// -1 -1  1: - x - y = -z: z = y + x
// -1  1 -1: - x + y =  z: y = x + z
// -1  1  1: - x + y = -z: x = y + z
//  1 -1 -1: + x - y =  z: x = y + z
//  1 -1  1: + x - y = -z: y = x + z
//  1  1 -1: + x + y =  z: z = x + y
//  1  1  1: + x + y = -z: 

void Mistral::OrderedSumExpression::extract_constraint(Solver *s) {
	s->add(Constraint(new ConstraintOrderedSum(children, lower_bound, upper_bound)));
}


// void Mistral::OrderedSumExpression::initialise_bounds() {
//   int tlb=0;
//   int tub=0;
//
//   for(unsigned int i=0; i<children.size; ++i) {
//     tlb += children[i].get_min();
//     tub += children[i].get_max();
//   }
//
//   if(tlb > lower_bound) lower_bound = tlb;
//   if(tub < upper_bound) upper_bound = tub;
// }

void Mistral::OrderedSumExpression::extract_variable(Solver *s) {
  std::cerr << "Error: OSum constraint can't yet be used as a predicate" << std::endl;
  exit(0);
  // initialise_bounds();
  //
  // Variable aux(lower_bound+offset, upper_bound+offset, DYN_VAR);
  // _self = aux;
  //
  // _self.initialise(s, 1);
  // _self = _self.get_var();
  // children.add(_self);
}

const char* Mistral::OrderedSumExpression::get_name() const {
  return "ordered_sum";
}

void Mistral::OrderedSumExpression::extract_predicate(Solver *s) {
  std::cerr << "Error: OSum constraint can't yet be used as a predicate" << std::endl;
  exit(0);
	//s->add(Constraint(new ConstraintOrderedSum(children, weight, -offset, -offset)));
}


Mistral::Variable Mistral::OSum(Vector< Variable >& args, const int l, const int u, const int offset) {
  Variable exp( new OrderedSumExpression(args, l, u,offset) );
  return exp;
}
Mistral::Variable Mistral::OSum(std::vector< Variable >& args, const int l, const int u, const int offset) {
  Variable exp( new OrderedSumExpression(args, l, u,offset) );
  return exp;
}



Mistral::LinearExpression::LinearExpression(Vector< Variable >& args, 
					    Vector< int >& wgts, 
					    const int l, const int u, const int o) 
  : Expression(args) {
  weighted = 0;
  bool_domains = 0;

  lower_bound = l;
  upper_bound = u;
  offset = o;
  for(unsigned int i=0; i<wgts.size; ++i) {
    weight.add(wgts[i]);
    if(wgts[i]==1) {
      if(weighted<0 || weighted>1) weighted = 2;
      else weighted = 1;
    } else if(wgts[i]==-1) {
      if(weighted>0) weighted = 2;
      else weighted = -1;
    } else weighted = 2;
    if(!children[i].is_boolean()) {
      if(bool_domains==0) bool_domains = i+1;
      else if(bool_domains>0) bool_domains = -1;
    } 
  }
}

Mistral::LinearExpression::LinearExpression(Vector< Variable >& args, 
					    const int l, const int u, const int o) 
  : Expression(args) {
  weighted = 1;
  bool_domains = 0;

  lower_bound = l;
  upper_bound = u;
  offset = o;
  for(unsigned int i=0; i<args.size; ++i) {
    weight.add(1);
    if(!children[i].is_boolean()) {
      if(bool_domains==0) bool_domains = i+1;
      else if(bool_domains>0) bool_domains = -1;
    } 
  }
}

Mistral::LinearExpression::LinearExpression(std::vector< Variable >& args, 
					    std::vector< int >& wgts, 
					    const int l, const int u, const int o) 
  : Expression(args) {
  weighted = 0;
  bool_domains = true;

  lower_bound = l;
  upper_bound = u;
  offset = o;
  for(unsigned int i=0; i<wgts.size(); ++i) {
    weight.add(wgts[i]);
    if(wgts[i]==1) {
      if(weighted<0 || weighted>1) weighted = 2;
      else weighted = 1;
    } else if(wgts[i]==-1) {
      if(weighted>0) weighted = 2;
      else weighted = -1;
    } else weighted = 2;
    if(!children[i].is_boolean()) {
      if(bool_domains==0) bool_domains = i+1;
      else if(bool_domains>0) bool_domains = -1;
    } 
  }
}

Mistral::LinearExpression::LinearExpression(std::vector< Variable >& args, 
					    const int l, const int u, const int o) 
  : Expression(args) {
  weighted = 1;
  bool_domains = 0;

  lower_bound = l;
  upper_bound = u;
  offset = o;
  for(unsigned int i=0; i<args.size(); ++i) {
    weight.add(1);
    if(!children[i].is_boolean()) {
      if(bool_domains==0) bool_domains = i+1;
      else if(bool_domains>0) bool_domains = -1;
    } 
  }
}

// Mistral::LinearExpression::LinearExpression(Variable X, const int coef) 
//   : Expression(X) { 
//   weight.add(coef);
// }
Mistral::LinearExpression::~LinearExpression() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete linear expression" << std::endl;
#endif
}
  
// -1 -1 -1: - x - y =  z: 
// -1 -1  1: - x - y = -z: z = y + x
// -1  1 -1: - x + y =  z: y = x + z
// -1  1  1: - x + y = -z: x = y + z
//  1 -1 -1: + x - y =  z: x = y + z
//  1 -1  1: + x - y = -z: y = x + z
//  1  1 -1: + x + y =  z: z = x + y
//  1  1  1: + x + y = -z: 

void Mistral::LinearExpression::extract_constraint(Solver *s) {
  // check if we can use an 'Add' or 'Sub' predicate
  int post_add = false;
  if(lower_bound == 0 &&
     upper_bound == 0 && 
     children.size == 3 && 
     abs(weight[0]) == 1 &&
     abs(weight[1]) == 1 &&
     abs(weight[2]) == 1
     ) {
    int i=0;
    for(; i<3; ++i)
      if(weight[i] != weight[(i+1)%3] && weight[(i+1)%3] == weight[(i+2)%3]) {
	Variable x = children[i];
	children[i] = children[2];
	children[2] = x;
	weight[0] = weight[1] = 1;
	weight[2] = -1;
	break;
      }
    if(i<3) {
      post_add = true;
      s->add(Constraint(new PredicateAdd(children)));
    }  
  }

  
  if(!post_add) {
    if(!bool_domains) {
      if(weighted) {
	s->add(Constraint(new ConstraintIncrementalWeightedBoolSumInterval(children, weight, lower_bound, upper_bound)));
      } else {
	s->add(Constraint(new ConstraintBoolSumInterval(children, lower_bound, upper_bound)));
      }
    } else if(lower_bound == 0 && upper_bound == 0 && bool_domains>0 && 
	      (weight[bool_domains-1] == 1 || weight[bool_domains-1] == -1)) {
      int n = children.size-1;
      Variable last = children[n];
      int wswap = weight[n];
      
      children[n] = children[bool_domains-1];
      weight[n] = weight[bool_domains-1];
      children[bool_domains-1] = last;
      weight[bool_domains-1] = wswap;
      
      if(weight[n] == 1) {
	for(int i=0; i<n; ++i)
	  weight[i] = -weight[i];
      }

      weight.pop();
      s->add(Constraint(new PredicateWeightedBoolSum(children, weight)));
    } else {
      s->add(Constraint(new PredicateWeightedSum(children, weight, lower_bound, upper_bound)));
    }
  }
}


void Mistral::LinearExpression::initialise_bounds() {

  int tlb=0;
  int tub=0;

  int lb = children[0].get_min()*weight[0];
  int ub = children[0].get_max()*weight[0];

  if(lb < ub) {
    tlb += lb;
    tub += ub;
  } else {
    tlb += ub;
    tub += lb;
  }

  for(unsigned int i=1; i<children.size; ++i) {
    lb = children[i].get_min()*weight[i];
    ub = children[i].get_max()*weight[i];
  
    if(lb < ub) {
      tlb += lb;
      tub += ub;
    } else {
      tlb += ub;
      tub += lb;
    }
  }

  if(tlb > lower_bound) lower_bound = tlb;
  if(tub < upper_bound) upper_bound = tub;

}

void Mistral::LinearExpression::extract_variable(Solver *s) {
  initialise_bounds();

  Variable aux(lower_bound+offset, upper_bound+offset, DYN_VAR);
  _self = aux;

  _self.initialise(s, 1);
  _self = _self.get_var();
  children.add(_self);
  weight.add(-1);
}

const char* Mistral::LinearExpression::get_name() const {
  return "linear_sum";
}

void Mistral::LinearExpression::extract_predicate(Solver *s) {
  s->add(Constraint(new PredicateWeightedSum(children, weight, -offset, -offset)));
}

Mistral::Variable Mistral::Sum(Vector< Variable >& args, Variable T, const int offset) {
  LinearExpression *lexpr = new LinearExpression(args,0,0,offset);
  lexpr->children.add(T);
  lexpr->weight.add(-1);
  if(lexpr->weighted==0)
    lexpr->weighted=-1;
  else if(lexpr->weighted==1)
    lexpr->weighted=2;
  Variable exp(lexpr);
  return exp;
}
Mistral::Variable Mistral::Sum(std::vector< Variable >& args, Variable T, const int offset) {
  LinearExpression *lexpr = new LinearExpression(args,0,0,offset);
  lexpr->children.add(T);
  lexpr->weight.add(-1);
  if(lexpr->weighted==0)
    lexpr->weighted=-1;
  else if(lexpr->weighted==1)
    lexpr->weighted=2;
  Variable exp(lexpr);
  return exp;
}
Mistral::Variable Mistral::Sum(Vector< Variable >& args, const int l, const int u, const int offset) {
  Variable exp( new LinearExpression(args, l, u,offset) );
  return exp;
}
Mistral::Variable Mistral::Sum(std::vector< Variable >& args, const int l, const int u, const int offset) {
  Variable exp( new LinearExpression(args, l, u,offset) );
  return exp;
}


Mistral::Variable Mistral::Sum(Vector< Variable >& args, Vector< int >& wgts, Variable T, const int offset) {
  LinearExpression *lexpr = new LinearExpression(args,wgts,0,0,offset);
  lexpr->children.add(T);
  lexpr->weight.add(-1);
  Variable exp(lexpr);
  return exp;
}
Mistral::Variable Mistral::Sum(std::vector< Variable >& args, std::vector< int >& wgts, Variable T, const int offset) {
  LinearExpression *lexpr = new LinearExpression(args,wgts,0,0,offset);
  lexpr->children.add(T);
  lexpr->weight.add(-1);
  Variable exp(lexpr);
  return exp;
}
Mistral::Variable Mistral::Sum(Vector< Variable >& args, Vector< int >& wgts, const int l, const int u, const int offset) {
  Variable exp( new LinearExpression(args, wgts, l, u,offset) );
  return exp;
}
Mistral::Variable Mistral::Sum(std::vector< Variable >& args, std::vector< int >& wgts, const int l, const int u, const int offset) {
  Variable exp( new LinearExpression(args, wgts, l, u,offset) );
  return exp;
}



Mistral::MinExpression::MinExpression(Vector< Variable >& args) 
  : Expression(args) {}

Mistral::MinExpression::~MinExpression() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete min expression" << std::endl;
#endif
}

void Mistral::MinExpression::extract_constraint(Solver *s) {
  std::cerr << "Error: Min predicate can't be used as a constraint" << std::endl;
  exit(0);
}

void Mistral::MinExpression::extract_variable(Solver *s) {
  int lower_bound = children[0].get_min();
  int upper_bound = children[0].get_max();

  int arity = children.size;
  for(int i=1; i<arity; ++i) {
    if(children[i].get_min() < lower_bound) lower_bound = children[i].get_min();
    if(children[i].get_max() < upper_bound) upper_bound = children[i].get_max();
  }

  Variable aux(lower_bound, upper_bound, DYN_VAR);
  _self = aux;

  _self.initialise(s, 1);
  _self = _self.get_var();
  children.add(_self);
}

const char* Mistral::MinExpression::get_name() const {
  return "min";
}

void Mistral::MinExpression::extract_predicate(Solver *s) {
  s->add(Constraint(new PredicateMin(children)));
}

Mistral::Variable Mistral::Min(Vector<Variable>& X) {
  bool equality = true;
  for(unsigned int i=1; i<X.size && equality; ++i) equality = (X[i-1].operator_equal(X[i]));
  Variable exp;
  if(equality) exp = X[0];
  else exp = Variable( new MinExpression(X) );
  return exp;
}

Mistral::Variable Mistral::Min(VarArray& X) {
  bool equality = true;
  for(unsigned int i=1; i<X.size && equality; ++i) equality = (X[i-1].operator_equal(X[i]));
  Variable exp;
  if(equality) exp = X[0];
  else exp = Variable( new MinExpression(X) );
  return exp;
  // Variable exp( new MinExpression(X) );
  // return exp;
}

Mistral::Variable Mistral::Min(Variable X, Variable Y) {
  Variable exp;
  if(!(X.operator_equal(Y))) {
    MinExpression *mexp = new MinExpression();
    mexp->add(X);
    mexp->add(Y);
    exp = Variable( mexp );
  } else exp = X;
  return exp;
}



Mistral::Variable Mistral::Max(Vector<Variable>& X) {
  Variable exp( new MaxExpression(X) );
  return exp;
}

Mistral::Variable Mistral::Max(VarArray& X) {
  Variable exp( new MaxExpression(X) );
  return exp;
}

Mistral::Variable Mistral::Max(Variable X, Variable Y) {
  MaxExpression *mexp = new MaxExpression();
  mexp->add(X);
  mexp->add(Y);
  Variable exp( mexp );
  return exp;
}


Mistral::MaxExpression::MaxExpression(Vector< Variable >& args) 
  : Expression(args) {}

Mistral::MaxExpression::~MaxExpression() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete max expression" << std::endl;
#endif
}

void Mistral::MaxExpression::extract_constraint(Solver *s) {
  std::cerr << "Error: Max predicate can't be used as a constraint" << std::endl;
  exit(0);
}

void Mistral::MaxExpression::extract_variable(Solver *s) {
  int lower_bound = children[0].get_min();
  int upper_bound = children[0].get_max();

  int arity = children.size;
  for(int i=1; i<arity; ++i) {
    if(children[i].get_min() > lower_bound) lower_bound = children[i].get_min();
    if(children[i].get_max() > upper_bound) upper_bound = children[i].get_max();
  }

  Variable aux(lower_bound, upper_bound, DYN_VAR);
  _self = aux;

  _self.initialise(s, 1);
  _self = _self.get_var();
  children.add(_self);
}

const char* Mistral::MaxExpression::get_name() const {
  return "max";
}

void Mistral::MaxExpression::extract_predicate(Solver *s) {
  s->add(Constraint(new PredicateMax(children)));
}


Mistral::ElementExpression::ElementExpression(const Vector< Variable >& args, 
					      Variable X, int ofs) 
  : Expression(), offset(ofs) {
  for(unsigned int i=0; i<args.size; ++i) children.add(args[i]);
  children.add(X);
}

Mistral::ElementExpression::~ElementExpression() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete element expression" << std::endl;
#endif
}


void Mistral::ElementExpression::extract_constraint(Solver *s) {
  std::cerr << "Error: Element predicate can't be used as a constraint" << std::endl;
  exit(0);
}

void Mistral::ElementExpression::initialise_domain() {

  lower_bound = INFTY;
  upper_bound = -INFTY;

  int arity = children.size-1;
  //Vector< int > values;

  BitSet domain;

  int i, nxt = children[arity].get_min();


  do {
    i = nxt-offset;
    if(i>=0 && i<arity) {
      if(children[i].get_min() < lower_bound) lower_bound = children[i].get_min();
      if(children[i].get_max() > upper_bound) upper_bound = children[i].get_max();
    }

    nxt = children[arity].next(i+offset);
  } while(i+offset<nxt);
  
  domain.initialise(lower_bound, upper_bound, BitSet::empt);
  
  nxt = children[arity].get_min();
  do {
    i = nxt-offset;
    
    if(i>=0 && i<arity) {
      children[i].union_to(domain);
    }

    nxt = children[arity].next(i+offset);
  } while(i+offset<nxt);


  nxt = lower_bound;
  do {
    i = nxt;
    values.add(i);
    nxt = domain.next(i);
  } while(i<nxt);
  
}

void Mistral::ElementExpression::extract_variable(Solver *s) {
  initialise_domain();
  Variable aux(values, DYN_VAR);


  _self = aux;

  _self.initialise(s, 1);
  _self = _self.get_var();
  children.add(_self);
}

const char* Mistral::ElementExpression::get_name() const {
  return "element";
}

//#define _DEBUG_ACELT

void Mistral::ElementExpression::extract_predicate(Solver *s) {
  int arity = children.size-2;
#ifdef _DEBUG_ACELT
  std::cout << "pre-propagte element(" 
	    << children[arity] << " = " << children[arity].get_domain() << " of ["
	    << children[0] ;
  for(unsigned int i=1; i<arity; ++i)
    std::cout << ", " << children[i];
  std::cout << "]) = " << children[arity+1] << " in " << children[arity+1].get_domain() << std::endl;
#endif
  if(FAILED(children[arity].set_min(offset)))
    { 
#ifdef _DEBUG_ACELT
      std::cout << " FAIL!" << std::endl;
      //exit(1);
#endif
      s->fail(); 
    }
  else if(FAILED(children[arity].set_max(arity-1+offset)))
    { 
#ifdef _DEBUG_ACELT
      std::cout << " FAIL!" << std::endl;
      //exit(1);
#endif
      s->fail(); 
    }
#ifdef _DEBUG_ACELT
  else
    std::cout << "ok\n";
#endif
	Constraint con(new PredicateElement(children, offset));
	
  s->add(con);
}


// Mistral::ElementSetExpression::ElementSetExpression(Vector< Variable >& args, 
// 						    Variable X, int ofs) 
//   : SetExpression(), offset(ofs) {
//   for(unsigned int i=0; i<args.size; ++i) children.add(args[i]);
//   children.add(X);
//   num_args = children.size;
//   initialise_domain();
//   initialise_elements();
// }

// Mistral::ElementSetExpression::~ElementSetExpression() {
// #ifdef _DEBUG_MEMORY
//   std::cout << "c delete element set expression" << std::endl;
// #endif
// }


// void Mistral::ElementSetExpression::extract_constraint(Solver *s) {
//   std::cerr << "Error: ElementSet predicate can't be used as a constraint" << std::endl;
//   exit(0);
// }

// // Mistral::Variable Mistral::ElementSetExpression::get_index_var(const int idx) {
// //   ((SetExpression*)(_self.variable))->children[i];
// // }

// void Mistral::ElementSetExpression::initialise_domain() {

//   // lb and ub stand for the cardinality's bounds

//   // val_ub and val_lb stand for the values' bounds
//   int val_lb = INFTY;
//   int val_ub = -INFTY;
//   lb = INFTY;
//   ub = -INFTY;

//   int arity = num_args-1;


//   // first, we go through the children to get the cardinality's and values' bounds

//   SetExpression *S;
//   int i, // j,
//     nxt = children[arity].get_min(), bound;
//   unsigned int j;
//   do {
//     i = nxt-offset;
    
//     if(i>=0 && i<arity) {
//       S = (SetExpression*)(children[i].variable);

//       bound = S->elts_ub.front();
//       if(bound<val_lb) val_lb = bound;

//       bound = S->elts_ub.back();
//       if(bound>val_ub) val_ub = bound;

//       bound = S->lb;
//       if(bound<lb) lb = bound;

//       bound = S->ub;
//       if(bound>ub) ub = bound;
//     }

//     nxt = children[arity].next(i+offset);
//   } while(i+offset<nxt);
  


//   // then we go once more over the children to get the lower and upper bounds of the set var
//   BitSet lb_domain(0,arity-1,BitSet::full);
//   BitSet ub_domain(0,arity-1,BitSet::empt);
//   BitSet aux(0,arity-1,BitSet::empt);
//   nxt = children[arity].get_min();

//   do {
//     i = nxt-offset;
    
//     if(i>=0 && i<arity) {
//       S = (SetExpression*)(children[i].variable);
//       for(j=0; j<S->elts_ub.size; ++j) {
// 	ub_domain.add(S->elts_ub[j]);
//       }

//       for(j=0; j<S->elts_lb.size; ++j) {
// 	aux.add(S->elts_lb[j]);
//       }

//       lb_domain.intersect_with(aux);
//       aux.clear();
//     }

//     nxt = children[arity].next(i+offset);
//   } while(i+offset<nxt);


//   //std::cout << lb_domain << " << S << " << ub_domain << std::endl;

//   i = ub_domain.min();
//   //j = ub_domain.max();

//   do {
//     elts_ub.add(i);
//     nxt = i;
//     i = ub_domain.next(nxt);
//   } while(nxt<i);

  
//  //  for(int elt=0; elt<arity; ++elt) { // for each element 
//  //    nxt = children[arity].get_min();
//  //  do {
//  //    i = nxt-offset;
    
//  //    if(i>=0 && i<arity) {
//  //      if(children[i].get_min()==0) {
//  // 	ub_domain.add()
//  //      }
//  // < lower_bound) lower_bound = children[i].get_min();
//  //      if(children[i].get_max() > upper_bound) upper_bound = children[i].get_max();
//  //    }

//  //    nxt = children[arity].next(i+offset);
//  //  } while(i+offset<nxt);
  

//  //  domain.initialise(lower_bound, upper_bound, BitSet::empt);
  
//  //  nxt = children[arity].get_min();
//  //  do {
//  //    i = nxt-offset;
    
//  //    if(i>=0 && i<arity) {
//  //      children[i].union_to(domain);
//  //    }

//  //    nxt = children[arity].next(i+offset);
//  //  } while(i+offset<nxt);


//  //  nxt = lower_bound;
//  //  do {
//  //    i = nxt;
//  //    values.add(i);
//  //    nxt = domain.next(i);
//  //  } while(i<nxt);
  
// }

// // void Mistral::ElementSetExpression::extract_variable(Solver *s) {
// //   initialise_domain();

// //   Variable aux(new SetExpression(elts_lb, elts_ub, lb, ub));
// //   _self = aux;

// //   _self.initialise(s, 1);
// //   //_self = _self.get_var();
// //   //children.add(_self);
// // }

// const char* Mistral::ElementSetExpression::get_name() const {
//   return "set_element";
// }

// void Mistral::ElementSetExpression::extract_predicate(Solver *s) {
//  int arity = num_args-1;
//  int i;
//  unsigned int j, k;

//  // std::cout << "HERE " << std::endl;
//  // std::cout << children[arity] << " in " << children[arity].get_domain() << std::endl;

//  //int nxt = children[arity].get_min();
//  int lb_index = children[arity].get_min();
//  int ub_index = children[arity].get_max();
//  SetExpression *S;
 
//  //VarArray scp[ub_index-lb_index+1];
//  Vector<Variable> scp[elts_ub.size];
 
//  for(i=lb_index; i<=ub_index; ++i) {
//    if(i-offset>=0 && i-offset<arity) {
//      if(children[arity].contain(i)) {

//        //std::cout << children[i-offset] << std::endl;

//        S = (SetExpression*)(children[i-offset].variable);

// 	   // std::cout << S->children << std::endl;
// 	   // std::cout << S->elts_ub << std::endl;
// 	   // std::cout << S->num_args << std::endl;


//        k=0;
//        for(j=0; j<elts_ub.size; ++j) {

// 	 //   std::cout << k << std::endl;
	 
// 	 // int hh = elts_ub[j];
// 	 // int jj = S->elts_ub[k];

// 	 // if(hh == jj) {

// 	 if(k<S->elts_ub.size && elts_ub[j] == S->elts_ub[k]) {
// 	   scp[j].add(S->get_index_var(k++));
// 	 } else {
// 	   Variable x(0);
// 	   scp[j].add(x);
// 	 }
//        }
//      } else {
//        for(j=0; j<elts_ub.size; ++j) {
// 	 Variable x(0);
// 	 scp[j].add(x);
//        }
//      }
//    }
//  }

 
//  for(j=0; j<elts_ub.size; ++j) {
//    //std::cout << scp[j] << std::endl;

//    scp[j].add(children[arity]);
//    scp[j].add(get_index_var(j));
   
//    s->add(new PredicateElement(scp[j], offset));
   
   
//    //std::cout << std::endl << std::endl << s << std::endl;
   
   
//    // s->add(Element(scp[j], children[arity].get_var(), offset) == 
//    // 	   ((SetExpression*)(_self.variable))->children[j].get_var());
//  }
 
//  scp[0].clear();
//  for(i=0; i<=arity; ++i) {
//    scp[0].add(children[i]);
//  }

//  //scp[0].add(self);
//  scp[0].add(s->variables[id]);
//  s->add(new PredicateElement(scp[0], offset));

//   // do {
//   //   i = nxt-offset;

//   //   if(i>=0 && i<arity) {
//   //     S = (SetExpression*)(children[i].variable);

//   //     k=0;
//   //     for(unsigned int j=0; j<elts_ub.size; ++j) {
//   // 	if(elts_ub[j] == S->elts_ub[k]) {
	  
//   // 	}
//   //     }

//   //   nxt = children[arity].next(i+offset);
//   // } while(i+offset<nxt);

//   //s->add(new PredicateElementSet(children, offset));
// }

Mistral::Variable Mistral::VarArray::operator[](Variable X) const {
  Variable exp( new ElementExpression(*this, X, 0) );
  return exp;
}

Mistral::Variable Mistral::VarArray::operator[](const int X) const {
  return stack_[X];
}

Mistral::Variable& Mistral::VarArray::operator[](const int X) {
  return stack_[X];
}

void Mistral::VarArray::set(const int X, Variable x) {
  stack_[X] = x;
}

Mistral::Variable Mistral::Element(const Vector<Variable>& X, Variable selector, int offset) {
  Variable exp( new ElementExpression(X, selector, offset) );
  return exp;
}

Mistral::Variable Mistral::Element(const VarArray& X, Variable selector, int offset) {
  Variable exp( new ElementExpression(X, selector, offset) );
  return exp;
}


// Mistral::Variable Mistral::ElementSet(Vector<Variable>& X, Variable selector, int offset) {
//   Variable exp( new ElementSetExpression(X, selector, offset) );
//   return exp;
// }

// Mistral::Variable Mistral::ElementSet(VarArray& X, Variable selector, int offset) {
//   Variable exp( new ElementSetExpression(X, selector, offset) );
//   return exp;
// }



// Mistral::IntersectionExpression::IntersectionExpression(Variable X, Variable Y) 
//   : SetExpression() {
//   children.add(X);
//   children.add(Y);
//   num_args = 2;
//   initialise_domain();
//   initialise_elements();
// }

// Mistral::IntersectionExpression::~IntersectionExpression() {
// #ifdef _DEBUG_MEMORY
//   std::cout << "c delete intersection expression" << std::endl;
// #endif
// }


// void Mistral::IntersectionExpression::extract_constraint(Solver *s) {
//   std::cerr << "Error: Intersection predicate can't be used as a constraint" << std::endl;
//   exit(0);
// }

// void Mistral::IntersectionExpression::initialise_domain() { 
//   SetExpression *X = (SetExpression*)(children[0].variable);
//   SetExpression *Y = (SetExpression*)(children[1].variable);

//   unsigned int i=0;
//   unsigned int j=0;

//   int x, y;
  
//   while(i<X->elts_ub.size && j<Y->elts_ub.size) {

//     x = X->elts_ub[i];
//     y = Y->elts_ub[j];
//     if(x == y) {
//       elts_ub.add(x);
//       ++i;++j;
//     } else if(x > y) { 
//       ++j;
//     } else {
//       ++i;
//     }
//   }

//   lb = 0;
//   ub = elts_ub.size;
// }

// const char* Mistral::IntersectionExpression::get_name() const {
//   return "set_intersection";
// }

// void Mistral::IntersectionExpression::extract_predicate(Solver *s) {
//   SetExpression *X = (SetExpression*)(children[0].variable);
//   SetExpression *Y = (SetExpression*)(children[1].variable);

//   unsigned int i=0;
//   unsigned int j=0;
//   unsigned int k=0;

//   int x, y;

//   VarArray scp;
  
//   while(i<X->elts_ub.size && j<Y->elts_ub.size) {

//     x = X->elts_ub[i];
//     y = Y->elts_ub[j];
//     if(x == y) {
//       scp.add(X->get_index_var(i));
//       scp.add(Y->get_index_var(j));
//       scp.add(get_index_var(k++));

//       s->add( Constraint(new PredicateAnd(scp)) );

//       scp.clear();
//       ++i;++j;
//     } else if(x > y) { 
//       ++j;
//     } else {
//       ++i;
//     }
//   }  

//   //scp.add(self);
//   scp.add(s->variables[id]);
//   scp.add(X);
//   s->add( Constraint(new ConstraintLess(scp)) );
 
//   scp.clear();
 
//   //scp.add(self);
//   scp.add(s->variables[id]);
//   scp.add(Y);
//   s->add( Constraint(new ConstraintLess(scp)) );
 
// }


Mistral::IntersectionExpression::IntersectionExpression(Variable X, Variable Y) 
  : SetExpression() {
  children.add(X);
  children.add(Y);
  num_args = 2;
  initialise_domain();
  //initialise_elements();
}

Mistral::IntersectionExpression::~IntersectionExpression() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete intersection expression" << std::endl;
#endif
}


void Mistral::IntersectionExpression::extract_constraint(Solver *s) {
  std::cerr << "Error: Intersection predicate can't be used as a constraint" << std::endl;
  exit(0);
}

void Mistral::IntersectionExpression::initialise_domain() { 
  SetExpression *X = (SetExpression*)(children[0].variable);
  SetExpression *Y = (SetExpression*)(children[1].variable);

  // index in X (var and lb, resp.)
  unsigned int i=0, il=0;
  // index in Y (var and lb, resp.)
  unsigned int j=0, jl=0;

  // x,y: value (in X and Y, resp.)
  int x, y, xl, yl;

  // the variable part is the intersection of the upper bounds minus the intersection of the lower bounds
  while(i<X->elts_var.size || il<X->elts_lb.size || 
	j<Y->elts_var.size || jl<Y->elts_lb.size) {

    x = (i<X->elts_var.size ? X->elts_var[i] : INFTY);
    y = (j<Y->elts_var.size ? Y->elts_var[j] : INFTY);
    xl = (il<X->elts_lb.size ? X->elts_lb[il] : INFTY);
    yl = (jl<Y->elts_lb.size ? Y->elts_lb[jl] : INFTY);

    if(xl<x && yl<y) {
      // both in the lower bound, we put it in the lower bound
      if(xl == yl) {
	elts_lb.add(xl);
	++il; ++jl;
      } else if(xl<yl) { 
	++il;
      } else if(xl>=yl) { 
	++jl;
      }
    } else {
      if(xl<x) { //(therefore y<yl)
	if(xl == y) {
	  elts_var.add(y);
	  children.add(Y->get_index_var(j));
	  map_y.add(j);
	  map_x.add(-1);
	  ++il; ++j;
	} else if(xl<y) ++il;
	else ++j;
      } else if(yl<y) { //(therefore x<xl)
	if(yl == x) {
	  elts_var.add(x);
	  children.add(X->get_index_var(i));
	  map_y.add(-1);
	  map_x.add(i);
	  ++jl; ++i;
	} else if(yl<x) ++jl;
	else ++i;
      } else {
	if(x == y) {
	  elts_var.add(x);
	  Variable x(0, 1);
	  children.add(x);
	  map_y.add(j);
	  map_x.add(i);
	  ++i; ++j;
	} else if(x<y) ++i;
	else ++j;
      }
    }
  }


  lower_bound = 0; //elts_lb.size;
  upper_bound = elts_var.size; //lb+elts_var.size;

  // std::cout << "intersection of " << children[0] << " and " << children[1] << ": " 
  // 	    << elts_lb << " <= this <= " << elts_var << " (car=[" << lb << "," << upper_bound << "])" << std::endl; 
}

const char* Mistral::IntersectionExpression::get_name() const {
  return "set_intersection";
}

void Mistral::IntersectionExpression::extract_predicate(Solver *s) {
  SetExpression *X = (SetExpression*)(children[0].variable);
  SetExpression *Y = (SetExpression*)(children[1].variable);

  Variable x, y, z;

  for(unsigned int i=0; i<elts_var.size; ++i) {
    if(map_x[i]>=0 && map_y[i]>=0) {
      z = get_index_var(i);
      x = X->get_index_var(map_x[i]);
      y = Y->get_index_var(map_y[i]);
      
      s->add(Constraint(new PredicateAnd(x, y, z)));
    }
  } 

  SetExpression::extract_predicate(s);
}


Mistral::Variable Mistral::Intersection(Variable X, Variable Y) {
  Variable exp(new IntersectionExpression(X, Y));
  return exp;
}



Mistral::UnionExpression::UnionExpression(Variable X, Variable Y) 
  : SetExpression() {
  children.add(X);
  children.add(Y);
  num_args = 2;
  initialise_domain();
  //initialise_elements();
}

Mistral::UnionExpression::~UnionExpression() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete union expression" << std::endl;
#endif
}


void Mistral::UnionExpression::extract_constraint(Solver *s) {
  std::cerr << "Error: Union predicate can't be used as a constraint" << std::endl;
  exit(0);
}

void Mistral::UnionExpression::initialise_domain() { 
  SetExpression *X = (SetExpression*)(children[0].variable);
  SetExpression *Y = (SetExpression*)(children[1].variable);

  // index in X (var and lb, resp.)
  unsigned int i=0, il=0;
  // index in Y (var and lb, resp.)
  unsigned int j=0, jl=0;

  // x,y: value (in X and Y, resp.)
  int x, y, xl, yl;

  // the variable part is the union of the upper bounds minus the union of the lower bounds
  while(i<X->elts_var.size || il<X->elts_lb.size || 
	j<Y->elts_var.size || jl<Y->elts_lb.size) {

    x = (i<X->elts_var.size ? X->elts_var[i] : INFTY);
    y = (j<Y->elts_var.size ? Y->elts_var[j] : INFTY);
    xl = (il<X->elts_lb.size ? X->elts_lb[il] : INFTY);
    yl = (jl<Y->elts_lb.size ? Y->elts_lb[jl] : INFTY);

    if(xl<x && yl<y) {
      // both in the lower bound, we put it in the lower bound
      if(xl == yl) {
	// both in the lower bound, if y is not in lb(x) it must be a fail
	elts_lb.add(xl);
	++il; ++jl;
      } else if(xl<yl) {
	// xl is in the lower bound of x, but cannot be in y
	elts_lb.add(xl);
	++il;
      } if(xl>yl) {
	// yl is in the lower bound of y, but cannot be in x
	elts_lb.add(yl);
	++jl;
      }
    } else {
      if(xl<x) { //(therefore y<yl)
	if(xl == y) {
	  // common value in the lower bound of x and the var part of y
	  elts_lb.add(xl);
	  ++il; ++j;
	} else if(xl<y) {
	  // xl is in the lower bound of x, but cannot be in y
	  elts_lb.add(xl);
	  ++il;
	} else {
	  // y is in the var part of y, but cannot be in x
	  elts_var.add(y);
	  children.add(Y->get_index_var(j));
	  map_y.add(j);
	  map_x.add(-1);
	  ++j;
	}
      } else if(yl<y) { //(therefore x<xl)
	if(yl == x) {
	  // common value in the lower bound of y and the var part of x
	  elts_lb.add(x);
	  ++jl; ++i;
	} else if(yl<x) {
	  // yl is in the lower bound of y, but cannot be in x
	  elts_lb.add(yl);
	  ++jl;
	}
	else {
	  // x is in the var part of x, but cannot be in y
	  elts_var.add(x);
	  children.add(X->get_index_var(i));
	  map_y.add(-1);
	  map_x.add(i);
	  ++i;
	}
      } else {
	if(x == y) {
	  // common value in the var part
	  elts_var.add(x);
	  Variable x(0, 1);
	  children.add(x);
	  map_y.add(j);
	  map_x.add(i);
	  ++i; ++j;
	} else if(x<y) {
	  // x is in the var part of x, but cannot be in y
	  elts_var.add(x);
	  children.add(X->get_index_var(i));
	  map_y.add(-1);
	  map_x.add(i);
	  ++i;
	} else {
	  // y is in the var part of y, but cannot be in x
	  elts_var.add(y);
	  children.add(Y->get_index_var(j));
	  map_y.add(j);
	  map_x.add(-1);
	  ++j;
	}
      }
    }
  }


  lower_bound = 0; //elts_lb.size;
  upper_bound = elts_var.size; //lb+elts_var.size;

  // std::cout << "union of " << children[0] << " and " << children[1] << ": " 
  // 	    << elts_lb << " <= this <= " << elts_var << " (car=[" << lb << "," << ub << "])" << std::endl; 
}

const char* Mistral::UnionExpression::get_name() const {
  return "set_union";
}

void Mistral::UnionExpression::extract_predicate(Solver *s) {
  SetExpression *X = (SetExpression*)(children[0].variable);
  SetExpression *Y = (SetExpression*)(children[1].variable);

  Variable x, y, z;

  for(unsigned int i=0; i<elts_var.size; ++i) {
    if(map_x[i]>=0 && map_y[i]>=0) {
      z = get_index_var(i);
      x = X->get_index_var(map_x[i]);
      y = Y->get_index_var(map_y[i]);
      
      s->add(Constraint(new PredicateOr(x, y, z)));
    }
  } 

  SetExpression::extract_predicate(s);
}


Mistral::Variable Mistral::Union(Variable X, Variable Y) {
  Variable exp(new UnionExpression(X, Y));
  return exp;
}




Mistral::SymmetricDifferenceExpression::SymmetricDifferenceExpression(Variable X, Variable Y) 
  : SetExpression() {
  children.add(X);
  children.add(Y);
  num_args = 2;
  initialise_domain();
  //initialise_elements();
}

Mistral::SymmetricDifferenceExpression::~SymmetricDifferenceExpression() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete union expression" << std::endl;
#endif
}


void Mistral::SymmetricDifferenceExpression::extract_constraint(Solver *s) {
  std::cerr << "Error: SymmetricDifference predicate can't be used as a constraint" << std::endl;
  exit(0);
}

void Mistral::SymmetricDifferenceExpression::initialise_domain() { 
  SetExpression *X = (SetExpression*)(children[0].variable);
  SetExpression *Y = (SetExpression*)(children[1].variable);

  // index in X (var and lb, resp.)
  unsigned int i=0, il=0;
  // index in Y (var and lb, resp.)
  unsigned int j=0, jl=0;

  // x,y: value (in X and Y, resp.)
  int x, y, xl, yl;

  // the variable part is the union of the upper bounds minus the union of the lower bounds
  while(i<X->elts_var.size || il<X->elts_lb.size || 
	j<Y->elts_var.size || jl<Y->elts_lb.size) {

    x = (i<X->elts_var.size ? X->elts_var[i] : INFTY);
    y = (j<Y->elts_var.size ? Y->elts_var[j] : INFTY);
    xl = (il<X->elts_lb.size ? X->elts_lb[il] : INFTY);
    yl = (jl<Y->elts_lb.size ? Y->elts_lb[jl] : INFTY);

    if(xl<x && yl<y) {
      // both in the lower bound, we put it in the lower bound
      if(xl == yl) {
	// both in the lower bound, if y is not in lb(x) it must be a fail
	++il; ++jl;
      } else if(xl<yl) {
	// xl is in the lower bound of x, but cannot be in y
	elts_lb.add(xl);
	++il;
      } if(xl>yl) {
	// yl is in the lower bound of y, but cannot be in x
	elts_lb.add(yl);
	++jl;
      }
    } else {
      if(xl<x) { //(therefore y<yl)
	if(xl == y) {
	  // common value in the lower bound of x and the var part of y
	  elts_var.add(xl);
	  children.add(!(Y->get_index_var(j)));
	  ++il; ++j;
	} else if(xl<y) {
	  // xl is in the lower bound of x, but cannot be in y
	  elts_lb.add(xl);
	  ++il;
	} else {
	  // y is in the var part of y, but cannot be in x
	  elts_var.add(y);
	  children.add(Y->get_index_var(j));
	  ++j;
	}
      } else if(yl<y) { //(therefore x<xl)
	if(yl == x) {
	  // common value in the lower bound of y and the var part of x
	  elts_var.add(x);
	  children.add(!(X->get_index_var(i)));
	  ++jl; ++i;
	} else if(yl<x) {
	  // yl is in the lower bound of y, but cannot be in x
	  elts_lb.add(yl);
	  ++jl;
	}
	else {
	  // x is in the var part of x, but cannot be in y
	  elts_var.add(x);
	  children.add(X->get_index_var(i));
	  ++i;
	}
      } else {
	if(x == y) {
	  // common value in the var part
	  map_y.add(j);
	  map_x.add(i);
	  map_z.add(elts_var.size);
	  elts_var.add(x);
	  Variable x(0, 1);
	  children.add(x);
	  ++i; ++j;
	} else if(x<y) {
	  // x is in the var part of x, but cannot be in y
	  elts_var.add(x);
	  children.add(X->get_index_var(i));
	  ++i;
	} else {
	  // y is in the var part of y, but cannot be in x
	  elts_var.add(y);
	  children.add(Y->get_index_var(j));
	  ++j;
	}
      }
    }
  }


  lower_bound = 0; //elts_lb.size;
  upper_bound = elts_var.size; //lb+elts_var.size;

  // std::cout << "union of " << children[0] << " and " << children[1] << ": " 
  // 	    << elts_lb << " <= this <= " << elts_var << " (car=[" << lb << "," << ub << "])" << std::endl; 
}

const char* Mistral::SymmetricDifferenceExpression::get_name() const {
  return "set_symdiff";
}

void Mistral::SymmetricDifferenceExpression::extract_predicate(Solver *s) {
  SetExpression *X = (SetExpression*)(children[0].variable);
  SetExpression *Y = (SetExpression*)(children[1].variable);

  Variable x, y, z;

  for(unsigned int i=0; i<map_x.size; ++i) {
    z = get_index_var(map_z[i]);
    x = X->get_index_var(map_x[i]);
    y = Y->get_index_var(map_y[i]);
    
    s->add(Constraint(new PredicateEqual(x, y, z, false)));
  } 

  SetExpression::extract_predicate(s);
}


Mistral::Variable Mistral::SymmetricDifference(Variable X, Variable Y) {
  Variable exp(new SymmetricDifferenceExpression(X, Y));
  return exp;
}


// Mistral::SetDifferenceExpression::SetDifferenceExpression(Variable X, Variable Y) 
//   : SetExpression() {
//   children.add(X);
//   children.add(Y);
//   num_args = 2;
//   initialise_domain();
//   initialise_elements();
// }

// Mistral::SetDifferenceExpression::~SetDifferenceExpression() {
// #ifdef _DEBUG_MEMORY
//   std::cout << "c delete intersection expression" << std::endl;
// #endif
// }


// void Mistral::SetDifferenceExpression::extract_constraint(Solver *s) {
//   std::cerr << "Error: SetDifference predicate can't be used as a constraint" << std::endl;
//   exit(0);
// }

// void Mistral::SetDifferenceExpression::initialise_domain() { 
//   SetExpression *X = (SetExpression*)(children[0].variable);
//   SetExpression *Y = (SetExpression*)(children[1].variable);

//   unsigned int i=0;
//   unsigned int j=0;

//   int x, y;
  
//   // a value is in the ub of x\y if it is in the ub of x and not in the lb of y
//   while(i<X->elts_ub.size) {
//     x = X->elts_ub[i];
//     y = (j < Y->elts_lb.size ? Y->elts_lb[j] : INFTY);
//     if(x < y) {
//       elts_ub.add(x);
//       ++i;
//     } else if(x > y) { 
//       ++j;
//     } else {
//       ++i; ++j;
//     }
//   }

//   // a value is in the lb of x\y if it is in the lb of x and not in the ub of y
//   while(i<X->elts_lb.size) {
//     x = X->elts_lb[i];
//     y = (j < Y->elts_ub.size ? Y->elts_ub[j] : INFTY);
//     if(x < y) {
//       elts_lb.add(x);
//       ++i;
//     } else if(x > y) { 
//       ++j;
//     } else {
//       ++i; ++j;
//     }
//   }

//   lb = elts_lb.size;
//   ub = elts_ub.size;
// }

// const char* Mistral::SetDifferenceExpression::get_name() const {
//   return "set_difference";
// }

// void Mistral::SetDifferenceExpression::extract_predicate(Solver *s) {
//   SetExpression *X = (SetExpression*)(children[0].variable);
//   SetExpression *Y = (SetExpression*)(children[1].variable);

//   unsigned int i=0;
//   unsigned int j=0;
//   unsigned int k=0;

//   int x, y, z;

//   VarArray scp;

//   // for each element e \in ub(X\Y):
//   // - if e \in ub(Y) then z(e) = y(e)<x(e)
//   // - otherwise z(e) = x(e)
  
//   while(k<elts_ub.size) {
//     z = elts_ub[k];
//     x = X->elts_ub[i];
//     y = (j < Y->elts_ub.size ? Y->elts_ub[j] : INFTY);

//     if(x<z) {
//       ++i;
//     } else {
//       if(z<x) {
// 	std::cerr << "Error when posting set diff" << std::endl;
// 	exit(1);
//       }
//       // we have x = z
//       if(x == y) {
// 	scp.add(Y->get_index_var(j));
// 	scp.add(X->get_index_var(i));
// 	scp.add(get_index_var(k));

// 	s->add( Constraint(new PredicateLess(scp, 1)) );

// 	scp.clear();
// 	++i;++j;++k;
//       } else if(x > y) { 
// 	++j;
//       } else {
// 	scp.add(X->get_index_var(i));
// 	scp.add(get_index_var(k));

// 	s->add( Constraint(new ConstraintEqual(scp)) );

// 	scp.clear();
// 	++i; ++k;
//       }
//     }
//   } 


//   //scp.add(self);
//   scp.add(s->variables[id]);
//   scp.add(children[0]);
//   s->add( Constraint(new ConstraintLess(scp)) );
 
//   // scp.clear();
 
//   // //scp.add(self);
//   // scp.add((X-Y));
//   // scp.add(s->variables[id]);
//   // s->add( Constraint(new ConstraintLess(scp)) );
//   s->add( (children[0]-children[1]) <= s->variables[id] );
 
// }


Mistral::SetDifferenceExpression::SetDifferenceExpression(Variable X, Variable Y) 
  : SetExpression() {
  children.add(X);
  children.add(Y);
  num_args = 2;
  initialise_domain();
  //initialise_elements();
}

Mistral::SetDifferenceExpression::~SetDifferenceExpression() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete intersection expression" << std::endl;
#endif
}


void Mistral::SetDifferenceExpression::extract_constraint(Solver *s) {
  std::cerr << "Error: SetDifference predicate can't be used as a constraint" << std::endl;
  exit(0);
}

void Mistral::SetDifferenceExpression::initialise_domain() { 
  SetExpression *X = (SetExpression*)(children[0].variable);
  SetExpression *Y = (SetExpression*)(children[1].variable);

  // std::cout << "decompose set diff " // << children[0] << " " << children[1] 
  // 	    << std::endl;


  // index in X (var and lb, resp.)
  unsigned int i=0, il=0;
  // index in Y (var and lb, resp.)
  unsigned int j=0, jl=0;

  // x,y: value (in X and Y, resp.)
  int x, y, xl, yl;

  // the variable part is the intersection of the upper bounds minus the intersection of the lower bounds
  while(i<X->elts_var.size || il<X->elts_lb.size || 
	j<Y->elts_var.size || jl<Y->elts_lb.size) {
    
    x = (i<X->elts_var.size ? X->elts_var[i] : INFTY);
    y = (j<Y->elts_var.size ? Y->elts_var[j] : INFTY);
    xl = (il<X->elts_lb.size ? X->elts_lb[il] : INFTY);
    yl = (jl<Y->elts_lb.size ? Y->elts_lb[jl] : INFTY);

    if(xl<x && yl<y) {
      // both in the lower bound, if y is not in lb(x) it must be a fail
      if(yl == xl) {
	// common value in the lower bound
	++il; ++jl;
      } else if(xl<yl) {
	// xl is in the lower bound of x, but cannot be in y
	elts_lb.add(xl);
	++il;
      } else if(xl>yl) {
	// yl is in the lower bound of y, but cannot be in x
	++jl;
      }
    } else {
      if(xl<x) { //(therefore y<yl) 
	if(xl == y) {
	  // common value in the lower bound of x and the var part of y
	  elts_var.add(y);
	  children.add(Variable(new NotExpression(Y->get_index_var(j))));
	  ++il; ++j;
	} else if(xl<y) {
	  // xl is in the lower bound of x, but cannot be in y
	  elts_lb.add(xl);
	  ++il;
	} else {
	  // y is in the var part of y, but cannot be in x
	  ++j;
	}
      } else if(yl<y) { //(therefore x<xl)
	if(yl == x) {
	  // common value in the lower bound of y and the var part of x
	  ++jl; ++i;
	} else if(yl<x) {
	  // yl is in the lower bound of y, but cannot be in x
	  ++jl;
	} else {
	  // x is in the var part of x, but cannot be in y
	  elts_var.add(x);
	  children.add(X->get_index_var(i));
	  ++i;
	}
      } else {
	if(x == y) {
	  // common value in the var part
	  elts_var.add(x);
	  children.add(Variable(new PrecedenceExpression(Y->get_index_var(j), X->get_index_var(i), 1)));	  
	  ++i; ++j;
	} else if(x<y) {
	  // x is in the var part of x, but cannot be in y
	  elts_var.add(x);
	  children.add(X->get_index_var(i));
	  ++i;
	}
	else {
	  // y is in the var part of y, but cannot be in x
	  ++j;
	}
      }
    }
  }    

  //std::cout << " -> "
    
  lower_bound = 0; //elts_lb.size;
  upper_bound = elts_var.size; //lb+elts_var.size;
}

const char* Mistral::SetDifferenceExpression::get_name() const {
  return "set_difference";
}

void Mistral::SetDifferenceExpression::extract_predicate(Solver *s) {
  SetExpression::extract_predicate(s);

  if(elts_var.size) {
    Variable c0 = Card(children[0]);
    Variable c_self = s->variables[id];

    //std::cout << c_self << " <= " << c0 << std::endl;

    s->add( c_self <= c0 );
    s->add( (Card(children[0])-Card(children[1])) <= s->variables[id] );
  }
}

Mistral::Variable Mistral::SetDifference(Variable X, Variable Y) {
  Variable exp(new SetDifferenceExpression(X, Y));
  return exp;
}


// Mistral::SetExpression::SetExpression(const int lelt, const int uelt, 
// 				      const int clb, const int cub) 
//   : BoolSumExpression(clb, cub) {
//   num_args = 0;
//   elts_ub.initialise(0, uelt-lelt+1);
//   for(int elt=lelt; elt<=uelt; ++elt) {
//     elts_ub.add(elt);
//     // Variable x(0, 1, BOOL_VAR);
//     // children.add(x);
//   }
//   initialise_elements();
// }


// // Mistral::SetExpression::SetExpression(const BitSet& lb, const BitSet& ub, const int clb, const int cub);

// Mistral::SetExpression::SetExpression(const Vector<int>& lb, const Vector<int>& ub, const int clb, const int cub)
//   : BoolSumExpression(clb, cub) {
//   // elts_ub.initialise(0, ub.back()-ub.front()+1);
//   // unsigned int i=0, j=0;
//   // while(i<ub.size) {
//   //   elts_ub.add(ub[i]);
//   //   if(lb.size>j && lb[j] == ub[i]) {
//   //     elts_lb.add(lb[j++]);
//   //     Variable c(1, 1, CONST_VAR);
//   //     children.add(c);
//   //   } else {
//   //     Variable x(0, 1, BOOL_VAR);
//   //     children.add(x);
//   //   }
//   //   ++i;
//   // }
//   num_args = 0;
//   elts_ub.initialise(0, ub.back()-ub.front()+1);
//   elts_ub.initialise(0, lb.back()-lb.front()+1);
//   unsigned int i=0;
//   while(i<ub.size) {
//     elts_ub.add(ub[i++]);
//   }

//   i=0;
//   while(i<lb.size) {
//     elts_lb.add(lb[i++]);
//   }
//   initialise_elements();
// }


// void Mistral::SetExpression::initialise_elements() {
//   unsigned int i=0, j=0;
//   while(i<elts_ub.size) {
//     if(elts_lb.size>j && elts_lb[j] == elts_ub[i]) {
//       ++j;
//       Variable c(1);
//       children.add(c);
//     } else {
//       Variable x(0, 1);
//       children.add(x);
//     }
//     ++i;
//   }
// }


// Mistral::SetExpression::~SetExpression() {
// #ifdef _DEBUG_MEMORY
//   std::cout << "c delete set expression" << std::endl;
// #endif
// }

// const char* Mistral::SetExpression::get_name() const {
//   return "set";
// }

// int Mistral::SetExpression::get_element_index(const int vali)  {
//   int target, xlb = 0, xub = elts_ub.size-1, idx = -1;
//   if(elts_ub.back() - elts_ub[0] == elts_ub.size+1) {
//     // if ub is a full interval, go directly to the index
//     idx = vali - elts_ub[0];
//   } else {
//     // else search by dichotomy
//     while(xlb <= xub) {
//       target = (xlb+xub)/2;
//       if(elts_ub[target] == vali) { idx=target; break; }
//       else if(elts_ub[target] > vali) xub = target-1;
//       else xlb = target+1;
//     }
//   }
//   return idx;
// }

// // Mistral::Variable Mistral::SetExpression::get_index_var(const int idx) {
// //   return children[idx];
// // }


// Mistral::Variable Mistral::SetExpression::get_elt_var(const int vali) {
//   // int target = 0, lb = -1, ub = elts_ub.size;
//   // while(lb < ub) {
//   //   target = lb+ub/2;
//   //   if(elts_ub[target] == vali) break;
//   //   else if(elts_ub[target] > vali) ub = target;
//   //   else lb = target;
//   // }
//   // return children[target];
//   int target = get_element_index(vali);
//   if(target>=0) return children[target];
//   Variable c(0);
//   return c;
// }

// std::ostream& Mistral::SetExpression::display(std::ostream& os) const {
//   // os << std::endl << elts_lb << std::endl
//   //    << elts_ub << std::endl;
//   // for(unsigned int i=0; i<elts_ub.size; ++i) {
//   //   os << " " << children[i].get_domain();
//   // }
//   // os << std::endl;

  
//   //std::cout << children << std::endl;


//   // bool first=true;
//   // os << "{" ;
//   // //if(children[0].get_min()) os << elts_ub[0];
//   // for(unsigned int i=0; i<elts_ub.size; ++i) {
//   //   if(children[num_args+i].get_min()) {
//   //     if(!first) os << ", " ;
//   //     else first = false; 
//   //     os << elts_ub[i];
//   //   }
//   // }
//   // first=true;
//   // os << "} <= S" << id << " <= {"; 
//   //if(children[0].get_max()) os << elts_ub[0];
//   // for(unsigned int i=0; i<elts_ub.size; ++i) {
//   //   if(children[num_args+i].get_max()) {
//   //     if(!first) os << ", " ;
//   //     else first = false;
//   //     os << elts_ub[i];
//   //   }
//   // }
//   // os << "}";

//   os << "S" << id;
//   return os;
// }

// int Mistral::SetExpression::get_solution_int_value() const { 
//   int i=elts_ub.size;
//   int t = 0;
//   int m = 1;
//   while(--i>=0) {
//     if((children[i].domain_type == CONST_VAR && children[i].constant_value == 1) ||
//        ((Solver*)solver)->last_solution_lb[children[i].id()]) {
//       t += m*elts_ub[i];
//       m *= 10;
//     }
//   }
//   return t;
// } 

// std::string Mistral::SetExpression::get_solution_str_value() const { 
//   std::ostringstream  ret_str;

//   ret_str << "{";
//   bool first = true;

//   for(unsigned int i=0; i<elts_ub.size; ++i) {
//     //std::cout << i << ": " << elts_ub[i] << " " << children[i].id() << " " << solver->last_solution_lb[children[i].id()] << std::endl;

//     if((children[i].domain_type == CONST_VAR && children[i].constant_value == 1) ||
//        ((Solver*)solver)->last_solution_lb[children[i].id()]) {
//       if(!first) {
// 	ret_str << ", " ;
//       } else {
// 	first = false;
//       }
//       ret_str << elts_ub[i];
//     }
//   }

//   ret_str << "}";

//   return ret_str.str();
// } 


// Mistral::SetExpression::SetExpression() 
//   : BoolSumExpression(0, 0) {
//   num_args = 0;
// }

Mistral::SetExpression::SetExpression(const int lelt, const int uelt, 
				      const int clb, const int cub) 
  : BoolSumExpression(clb, cub) {
  num_args = 0;

  //std::cout << lb << " < " << (uelt-lelt+1) << std::endl;

  if(lower_bound < uelt-lelt+1) {
    elts_var.initialise(0, uelt-lelt+1);
    for(int elt=lelt; elt<=uelt; ++elt) {
      elts_var.add(elt);
    }
  } else {
    elts_lb.initialise(0, uelt-lelt+1);
    for(int elt=lelt; elt<=uelt; ++elt) {
      elts_lb.add(elt);
    }
  }

  //std::cout << elts_lb << " <= this <= " << elts_var << " (car=[" << lb << "," << ub << "])" << std::endl; 

  initialise_elements();
}


Mistral::SetExpression::SetExpression(const Vector<int>& l, const Vector<int>& u, const int clb, const int cub)
  : BoolSumExpression(clb, cub) {
  num_args = 0;
  unsigned int i=0, j;
  //bool need_sorting = false;
  while(i<l.size) {
    //need_sorting = (i && l[i]<l[i-1]);
    elts_lb.add(l[i++]);
  }
  //if(need_sorting)
    

  i=0;
  j=0;
  int x, y;
  while(i<u.size) {
    x = u[i];
    y = (j < l.size ? l[j] : INFTY);
    if(x < y) elts_var.add(x);
    if(x <= y) ++i;
    if(x >= y) ++j;
  }

  lower_bound -= elts_lb.size;
  upper_bound -= elts_lb.size;

  if(lower_bound<0) lower_bound = 0;
  if(upper_bound>(int)(elts_var.size)) upper_bound = elts_var.size;

  // std::cout << "elts_lb: " << elts_lb << std::endl; 
  // std::cout << "elts_var: " << elts_var << std::endl; 

  // //" <= this <= " << elts_lb << "u" << elts_var << " (car=[" << lb+elts_lb.size << "," << ub+elts_lb.size << "])" << std::endl; 

  initialise_elements();
}


void Mistral::SetExpression::initialise_elements() {
  unsigned int i=0;
  while(i++<elts_var.size) {
    Variable x(0, 1);
    children.add(x);
  }
}


Mistral::SetExpression::~SetExpression() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete set expression" << std::endl;
#endif
}


void Mistral::SetExpression::extract_predicate(Solver *s) { 
  //s->add(new PredicateBoolSum(children, self)); 
  Vector<Variable> scp;
  for(unsigned int i=num_args; i<children.size; ++i) scp.add(children[i]);
  
  //std::cout << s->variables.size << " " << id << std::endl;

  //std::cout << lb << " " << ub << std::endl;

  if(elts_var.size) {
    s->add(Constraint(new PredicateBoolSum(scp, s->variables[id]))); 
  } else {
    if(lower_bound) s->fail();
  }

}


const char* Mistral::SetExpression::get_name() const {
  return "set";
}

int Mistral::SetExpression::get_element_index(const int vali)  {
  int target, xlb = 0, xub = elts_var.size-1, idx = -1;
  if(elts_var.size && elts_var[0] <= vali && vali <= elts_var.back()) {
    if(elts_var.back() - elts_var[0] == (int)(elts_var.size+1)) {
      // if var is a full interval, go directly to the index
      idx = vali - elts_var[0];
    } else {
      // else search by dichotomy
      while(xlb <= xub) {
	target = (xlb+xub)/2;
	if(elts_var[target] == vali) { idx=target; break; }
	else if(elts_var[target] > vali) xub = target-1;
	else xlb = target+1;
      }
    }
  }
  return idx;
}

int Mistral::SetExpression::get_element_lb_index(const int vali)  {
  int target, xlb = 0, xub = elts_lb.size-1, idx = -1;
  if(elts_lb.size && elts_lb[0] <= vali && vali <= elts_lb.back()) {
    if(elts_lb.back() - elts_lb[0] == (int)(elts_lb.size+1)) {
      // if lb is a full interval, go directly to the index
      idx = vali - elts_lb[0];
    } else {
      // else search by dichotomy
      while(xlb <= xub) {
	target = (xlb+xub)/2;
	if(elts_lb[target] == vali) { idx=target; break; }
	else if(elts_lb[target] > vali) xub = target-1;
	else xlb = target+1;
      }
    }
  }
  return idx;
}


Mistral::Variable Mistral::SetExpression::get_elt_var(const int vali) {
  int target = get_element_index(vali);
  if(target>=0) return children[target];
  
  Variable x;
  return x;

  /*
    else target = get_element_lb_index(vali);
    if(target>=0) {
    Variable c(1);
    return c;
    }
    Variable c(0);
    return c;
  */
}

std::ostream& Mistral::SetExpression::display(std::ostream& os, const bool domain) const {
  // os << std::endl << elts_lb << std::endl
  //    << elts_ub << std::endl;
  // for(unsigned int i=0; i<elts_ub.size; ++i) {
  //   os << " " << children[i].get_domain();
  // }
  // os << std::endl;

  
  //std::cout << children << std::endl;


  // bool first=true;
  // os << "{" ;
  // //if(children[0].get_min()) os << elts_ub[0];
  // for(unsigned int i=0; i<elts_ub.size; ++i) {
  //   if(children[num_args+i].get_min()) {
  //     if(!first) os << ", " ;
  //     else first = false; 
  //     os << elts_ub[i];
  //   }
  // }
  // first=true;
  // os << "} <= S" << id << " <= {"; 
  //if(children[0].get_max()) os << elts_ub[0];
  // for(unsigned int i=0; i<elts_ub.size; ++i) {
  //   if(children[num_args+i].get_max()) {
  //     if(!first) os << ", " ;
  //     else first = false;
  //     os << elts_ub[i];
  //   }
  // }
  // os << "}";


  if(domain) {
    
    os << "{";
    bool first = true;
    
    unsigned int i=0, j=0;
    int x, y;

    while(i<elts_var.size || j<elts_lb.size) {
      x = (i < elts_var.size ? elts_var[i] : INFTY);
      y = (j < elts_lb.size ? elts_lb[j] : INFTY);
      
      if(x<y) {
	if(children[num_args+i].get_min()) {
	  if(!first) os << ", " ;
	  else first = false;
	  os << elts_var[i];
	} else if(children[num_args+i].get_max()) {
	  if(!first) os << ", " ;
	  else first = false;
	  os << "?" << elts_var[i];
	}
	++i;
      } else {
	if(!first) os << ", " ;
	else first = false;
	os << elts_lb[j];
	++j;
      }
    }
    os << "}";
    
    // os << "[" << lb << " " << ub << "] " << elts_lb << " " << elts_var ;
    // for(unsigned int i=num_args; i<children.size; ++i) {
    //   os << " " << children[i] << " in " << children[i].get_domain() ;
    // }
    // os << " ";
  } else {

    os << "S" << id;

  }

  return os;
}

int Mistral::SetExpression::get_solution_int_value() const { 
  //int i=elts_ub.size;
  int t = 0;
  // int m = 1;
  // while(--i>=0) {
  //   if((children[i].domain_type == CONST_VAR && children[i].constant_value == 1) ||
  //      ((Solver*)solver)->last_solution_lb[children[i].id()]) {
  //     t += m*elts_ub[i];
  //     m *= 10;
  //   }
  // }
  return t;
} 

std::string Mistral::SetExpression::get_solution_str_value() const { 
  std::ostringstream  ret_str;

  ret_str << "{";
  bool first = true;

  unsigned int i=0, j=0;
  int x, y;

  while(i<elts_var.size || j<elts_lb.size) {
    x = (i < elts_var.size ? elts_var[i] : INFTY);
    y = (j < elts_lb.size ? elts_lb[j] : INFTY);

    if(x<y) {
      if(children[num_args+i].get_solution_min()) {
	if(!first) ret_str << ", " ;
	else first = false;
	ret_str << elts_var[i];
      } else if(children[num_args+i].get_solution_max()) {
	if(!first) ret_str << ", " ;
	else first = false;
	ret_str << "?" << elts_var[i];
      }
      ++i;
    } else {
      if(!first) ret_str << ", " ;
      else first = false;
      ret_str << elts_lb[j];
      ++j;
    }
  }
  ret_str << "}";

  return ret_str.str();
} 

Mistral::Variable Mistral::SetVariable() {
  SetExpression *se = new SetExpression();
  Variable exp(se);
  return exp;
}

Mistral::Variable Mistral::SetVariable(const int lelt, const int uelt, const int clb, const int cub) {
  int ub = cub;
  if(cub > uelt-lelt+1) ub = (uelt-lelt+1);
  SetExpression *se = new SetExpression(lelt, uelt, clb, ub);
  //std::cout << "create " << se << std::endl;
  Variable exp(se);
  return exp;
}

Mistral::Variable Mistral::SetVariable(const Vector<int>& slb, const Vector<int>& sub, 
				       const int clb, const int cub) {
  SetExpression *se = new SetExpression(slb, sub, clb, cub);
  Variable exp(se);
  return exp;
}

Mistral::Variable Mistral::Card(Variable S) { 
  int lb_size = ((SetExpression*)(S.expression))->elts_lb.size;
  if(lb_size) {
    // std::cout << "card of ";
    // ((SetExpression*)(S.expression))->display(std::cout, true);
    // std::cout << std::endl;
    Variable C(new OffsetExpression(S, lb_size));
    return C;
  }
  Variable cs = Variable(new IdExpression(S));
  return cs; //(lb_size ? (S+lb_size) : S); 
}

// Mistral::SubsetExpression::SubsetExpression(Variable X, Variable Y) 
//   : Expression(X,Y) { 
// }

// Mistral::SubsetExpression::~SubsetExpression() {
// #ifdef _DEBUG_MEMORY
//   std::cout << "c delete subset expression" << std::endl;
// #endif
// }

// void Mistral::SubsetExpression::extract_constraint(Solver *s) { 
//   unsigned int i = 0, j = 0;
//   SetExpression *x = (SetExpression*)(children[0].expression);
//   SetExpression *y = (SetExpression*)(children[1].expression);


//   // if there are k elements of Y that can't be in X, then card(X) + k <= card(Y)
//   Vector<Variable> extras;
//   int c_offset=0;
  
//   while(i<x->elts_ub.size && j<y->elts_ub.size) {
//     if(x->elts_ub[i] == y->elts_ub[j]) {
//       s->add(x->get_index_var(i) <= y->get_index_var(j));
//       ++i;
//       ++j;
//     } else if(x->elts_ub[i] <= y->elts_ub[j]) {
//       s->add(x->get_index_var(i) == 0);
//       ++i;
//     } else {
//       if(y->get_index_var(j).get_min()>0) {
// 	++c_offset;
//       } else {
// 	extras.add(y->get_index_var(j));
//       }
//       ++j;
//     }
//   }

//   // remaining values of X
//   while(i<x->elts_ub.size) {
//     s->add(x->get_index_var(i) == 0);
//     ++i;
//   }

//   // remaining values of Y
//   while(j<y->elts_ub.size) {
//     if(y->get_index_var(j).get_min()>0) {
//       ++c_offset;
//     } else {
//       extras.add(y->get_index_var(j));
//     }
//     ++j;
//   }


//   if(extras.size>0) {
//     //extras.add(children[0].get_var());
//     s->add(Precedence((BoolSum(extras) + children[0]), c_offset, children[1]));
//     //s->add(children[0].get_var() + Sum() <= children[1].get_var());
//   } else {
//     s->add(Precedence(children[0], c_offset, children[1]));
//     //s->add(children[0].get_var() <= children[1].get_var());
//   }

//   //s->add(children[0].get_var() <= children[1].get_var());
// }

// void Mistral::SubsetExpression::extract_variable(Solver *s) {
//   std::cerr << "Error: Subset constraint can't yet be used as a predicate" << std::endl;
//   exit(0);
// }

// void Mistral::SubsetExpression::extract_predicate(Solver *s) { 
//   std::cerr << "Error: Subset constraint can't yet be used as a predicate" << std::endl;
//   exit(0);
// }

// const char* Mistral::SubsetExpression::get_name() const {
//   return "subset";
// }

Mistral::Variable Mistral::Subset(Variable X, Variable Y) {
  Variable exp(new SubsetExpression(X, Y));
  return exp;
}

Mistral::Variable Mistral::NotSubset(Variable X, Variable Y) {
  Variable exp(new SubsetExpression(X, Y, 0));
  return exp;
}


Mistral::SubsetExpression::SubsetExpression(Variable X, Variable Y, const int sp) 
  : Expression(X,Y) { 
  spin = sp;
}

Mistral::SubsetExpression::~SubsetExpression() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete subset expression" << std::endl;
#endif
}

void Mistral::SubsetExpression::extract_constraint(Solver *s) { 
  SetExpression *X = (SetExpression*)(children[0].variable);
  SetExpression *Y = (SetExpression*)(children[1].variable);
  
  // index in X (var and lb, resp.)
  unsigned int i=0, il=0;
  // index in Y (var and lb, resp.)
  unsigned int j=0, jl=0;

  // x,y: value (in X and Y, resp.)
  int x, y, xl, yl;

  Vector<Variable> extras;
  Vector<Variable> neg;

  // the variable part is the intersection of the upper bounds minus the intersection of the lower bounds
  while(i<X->elts_var.size || il<X->elts_lb.size || 
	j<Y->elts_var.size || jl<Y->elts_lb.size) {

    x = (i<X->elts_var.size ? X->elts_var[i] : INFTY);
    y = (j<Y->elts_var.size ? Y->elts_var[j] : INFTY);
    xl = (il<X->elts_lb.size ? X->elts_lb[il] : INFTY);
    yl = (jl<Y->elts_lb.size ? Y->elts_lb[jl] : INFTY);

    //std::cout << "x=" << x << ", xl=" << xl << ", y=" << y << ", yl=" << yl << std::endl;


    if(xl<x && yl<y) {
      // both in the lower bound, if y is not in lb(x) it must be a fail
      if(yl>xl) {
	s->fail();
	break;
      }
      if(xl<=yl) {
	++il;
      }
      if(xl>=yl) {
	++jl;
      }
    } else {
      if(xl<x) { //(therefore y<yl) 
	if(xl == y) {
	  if(spin) {
	    extras.add(Y->get_index_var(j));
	    s->add(Y->get_index_var(j) == 1);
	  } else {
	    extras.add(Y->get_index_var(j));
	  }
	  ++il; ++j;
	} else if(xl<y) {
	  if(spin) {
	    s->fail();
	  }
	  break;
	  //++il;
	} else {
	  if(spin) 
	    extras.add(Y->get_index_var(j));
	  //s->add(Y->get_index_var(j) == 0);
	  ++j;
	}
      } else if(yl<y) { //(therefore x<xl)
	if(yl == x) {
	  //s->add(X->get_index_var(i) == 1);
	  ++jl; ++i;
	} else if(yl<x) {
	  ++jl;
	} else {
	  if(spin)
	    s->add(X->get_index_var(i) == 0);
	  else
	    neg.add(X->get_index_var(i));
	  ++i;
	}
      } else {
	if(x == y) {
	  if(spin) 
	    s->add(Constraint(new ConstraintLess(X->get_index_var(i), Y->get_index_var(j))));
	  else
	    extras.add(Precedence(X->get_index_var(i), 0, Y->get_index_var(j)));
	  ++i; ++j;
	} else if(x<y) {
	  if(spin)
	    s->add(X->get_index_var(i) == 0);
	  else
	    neg.add(X->get_index_var(i));
	  ++i;
	}
	else {
	  if(spin)
	    extras.add(Y->get_index_var(j));
	  //s->add(Y->get_index_var(j) == 0);
	  ++j;
	}
      }
    }
  }


  if(spin) {
    if(extras.size>0) {
      s->add(Precedence((BoolSum(extras) + children[0]), 0, children[1]));
    } else {
      s->add(Precedence(children[0], 0, children[1]));
    }
  } else {
    int all = extras.size;
    if(extras.size) {
      if(neg.size) {
	s->add((BoolSum(extras) < all) || (BoolSum(neg) > 0));
      } else {
	s->add(BoolSum(extras) < all);
      } 
    } else if(neg.size) {
      s->add(BoolSum(neg) > 0);;
    } else {
      s->fail();
    }
  }
}

void Mistral::SubsetExpression::extract_variable(Solver *s) {
  Variable aux(0, 1, BOOL_VAR);
  _self = aux;
  
  _self.initialise(s, 1);
  _self = _self.get_var();
  children.add(_self);
  // std::cerr << "Error: Subset constraint can't yet be used as a predicate" << std::endl;
  // exit(0);
}

void Mistral::SubsetExpression::extract_predicate(Solver *s) { 
  SetExpression *X = (SetExpression*)(children[0].variable);
  SetExpression *Y = (SetExpression*)(children[1].variable);
  
  // index in X (var and lb, resp.)
  unsigned int i=0, il=0;
  // index in Y (var and lb, resp.)
  unsigned int j=0, jl=0;
  
  // x,y: value (in X and Y, resp.)
  int x, y, xl, yl;
  
  Vector<Variable> implies_true;
  Vector<Variable> implies_false;
  //bool satisfied = false;
  bool violated = false;
  
  // the variable part is the intersection of the upper bounds minus the intersection of the lower bounds
  while(i<X->elts_var.size || il<X->elts_lb.size || 
	j<Y->elts_var.size || jl<Y->elts_lb.size) {
    
    x = (i<X->elts_var.size ? X->elts_var[i] : INFTY);
    y = (j<Y->elts_var.size ? Y->elts_var[j] : INFTY);
    xl = (il<X->elts_lb.size ? X->elts_lb[il] : INFTY);
    yl = (jl<Y->elts_lb.size ? Y->elts_lb[jl] : INFTY);

    //std::cout << "x=" << x << ", xl=" << xl << ", y=" << y << ", yl=" << yl << std::endl;


    if(xl<x && yl<y) {
      // both in the lower bound, if y is not in lb(x) it must be a fail
      if(yl == xl) {
	// common value in the lower bound
	++il; ++jl;
      } else if(xl<yl) {
	// xl is in the lower bound of x, but cannot be in y
	violated = true;
	++il;
      } else if(xl>yl) {
	// yl is in the lower bound of y, but cannot be in x
	++jl;
      }
    } else {
      if(xl<x) { //(therefore y<yl) 
	if(xl == y) {
	  // common value in the lower bound of x and the var part of y
	  implies_true.add(Y->get_index_var(j));
	  ++il; ++j;
	} else if(xl<y) {
	  // xl is in the lower bound of x, but cannot be in y
	  violated = true;
	  ++il;
	} else {
	  // y is in the var part of y, but cannot be in x
	  ++j;
	}
      } else if(yl<y) { //(therefore x<xl)
	if(yl == x) {
	  // common value in the lower bound of y and the var part of x
	  ++jl; ++i;
	} else if(yl<x) {
	  // yl is in the lower bound of y, but cannot be in x
	  ++jl;
	} else {
	  // x is in the var part of x, but cannot be in y
	  implies_false.add(X->get_index_var(i));
	  ++i;
	}
      } else {
	if(x == y) {
	  // common value in the var part
    	  implies_true.add((X->get_index_var(i) <= Y->get_index_var(j)));
	  ++i; ++j;
	} else if(x<y) {
	  // x is in the var part of x, but cannot be in y
	  implies_false.add(X->get_index_var(i));
	  ++i;
	}
	else {
	  // y is in the var part of y, but cannot be in x
	  ++j;
	}
      }
    }
  }

  // std::cout << s << std::endl;
  // std::cout << "disjunction_pos: " << disjunction_pos << std::endl;
  // std::cout << "disjunction_neg: " << disjunction_neg << std::endl;

  // z <-> (t1 & t2 & ... & tn) & (f1 & f2 & ... & fm)
 

  if(violated) {
    if(FAILED(children[2].remove(spin))) s->fail();
  } else {
    int all;
    if(spin) {
      all = implies_true.size;
      for(unsigned int i=0; i<implies_true.size; ++i) {
	s->add(children[2] <= implies_true[i]);
      }
      for(unsigned int i=0; i<implies_false.size; ++i) {
	s->add(Constraint(new ConstraintNotAnd(children[2], implies_false[i])));
      }
      if(implies_true.size) {
	if(implies_false.size) {
	  s->add(((BoolSum(implies_true) == all) && (BoolSum(implies_false) == 0)) <= children[2]);
	} else {
	  s->add((BoolSum(implies_true) == all) <= children[2]);
	}
      } else if(implies_false.size) {
	s->add((BoolSum(implies_false) == 0) <= children[2]);
      } else {
	if(FAILED(children[2].set_domain(spin))) s->fail();
      }
    } else {

      all = implies_true.size;
      for(unsigned int i=0; i<implies_false.size; ++i) {
	s->add(implies_false[i] <= children[2]);
      }
      for(unsigned int i=0; i<implies_true.size; ++i) {
	s->add(Variable(new OrExpression(children[2], implies_true[i])));
      }
      if(implies_true.size) {
	if(implies_false.size) {
	  s->add(children[2] <= ((BoolSum(implies_true) < all) || (BoolSum(implies_false) > 0)));
	} else {
	  s->add(children[2] <= (BoolSum(implies_true) < all));
	}
      } else if(implies_false.size) {
	s->add(children[2] <= (BoolSum(implies_false) > 0));
      } else {
	if(FAILED(children[2].set_domain(spin))) s->fail();
      }
    }
  }
 
  // std::cerr << "Error: Subset constraint can't yet be used as a predicate" << std::endl;
  // exit(0);
}

const char* Mistral::SubsetExpression::get_name() const {
  return "subset";
}




Mistral::MemberExpression::MemberExpression(Variable X, Variable Y, const int sp) 
  : Expression(X,Y) { 
  lb = +INFTY;
  ub = -INFTY;
  size = 0;
  spin = sp;
}

Mistral::MemberExpression::MemberExpression(Variable X, const int lo, const int up, const int sp) 
  : Expression(X) { 
  lb = lo;
  ub = up;
  size = ub-lb+1;
  spin = sp;
}

Mistral::MemberExpression::MemberExpression(Variable X, const BitSet& s, const int sp) 
  : Expression(X) {
  lb = s.min();
  ub = s.max();
  size = s.size();
  spin = sp;
  values.initialise(s);
}

Mistral::MemberExpression::MemberExpression(Variable X, const Vector<int>& s, const int sp) 
  : Expression(X) {
  //lb = s.front();
  //ub = s.back();
  size = s.size; 
  spin = sp;
  lb = +INFTY;
  ub = -INFTY;
  for(unsigned int i=0; i<s.size; ++i) {
    if(s[i] < lb) lb = s[i];
    if(s[i] > ub) ub = s[i];
  }
  values.initialise(lb, ub, BitSet::empt);
  for(unsigned int i=0; i<s.size; ++i) {
    values.add(s[i]);
  }
}

Mistral::MemberExpression::MemberExpression(Variable X, const std::vector<int>& s, const int sp) 
  : Expression(X) {
  //lb = s.front();
  //ub = s.back();
  size = s.size(); 
  spin = sp;
  lb = +INFTY;
  ub = -INFTY;
  for(unsigned int i=0; i<s.size(); ++i) {
    if(s[i] < lb) lb = s[i];
    if(s[i] > ub) ub = s[i];
  }
  values.initialise(lb, ub, BitSet::empt);
  for(unsigned int i=0; i<s.size(); ++i) {
    values.add(s[i]);
  }
}

Mistral::MemberExpression::~MemberExpression() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete member expression" << std::endl;
#endif
}

void Mistral::MemberExpression::extract_constraint(Solver *s) { 

  if(spin) {
    //unsigned int i = 0, j = 0;
    if(children.size == 2) {
      // set variable
      VarArray scp;
      SetExpression *y = (SetExpression*)(children[1].expression);
      s->add(Card(children[1]) >= 1); // at least one element in the set
      int vali, vnxt = children[0].get_min(), idx;
      do {
	scp.clear();
	vali = vnxt;
	idx = y->get_element_index(vali);
	if(idx>=0) {
	  Variable x(0,1);
	  s->add(x);
	  scp.add(children[0]);
	  scp.add(x);
	  s->add(Constraint(new PredicateConstantEqual(scp,vali)));
	  scp.clear();
	  scp.add(x);
	  scp.add(y->get_index_var(idx));
	  s->add(Constraint(new ConstraintLess(scp)));
	  //s->add((children[0] == vali) <= y->children[idx]);
	} else if(y->get_element_lb_index(vali)<0)
	  if(FAILED(children[0].remove(vali))) {
	    s->fail();
	    break;
	  }
	vnxt = children[0].next(vali);
      } while(vali < vnxt);
    } else if(size == (ub-lb+1)) {
      // interval 
      if(FAILED(children[0].set_min(lb))) { 
	s->fail(); }
      if(children[0].set_max(ub) == FAIL_EVENT) { 
	s->fail(); 
      }
    } else {
      // set
      if(FAILED(children[0].set_domain(values))) { 
	s->fail(); 
      }
    }
  } else {
    // not member
    //unsigned int i = 0, j = 0;
    if(children.size == 2) {
      // set variable
      VarArray scp;
      SetExpression *y = (SetExpression*)(children[1].expression);
      //s->add(children[1] >= 1); // at least one element in the set
      int vali, vnxt = children[0].get_min(), idx;
      do {
	scp.clear();
	vali = vnxt;
	idx = y->get_element_index(vali);
	if(idx>=0) {
	  Variable x(0,1);
	  s->add(x);
	  scp.add(children[0]);
	  scp.add(x);
	  s->add(Constraint(new PredicateConstantEqual(scp,vali,0)));
	  scp.clear();
	  scp.add(y->get_index_var(idx));
	  scp.add(x);
	  s->add(Constraint(new ConstraintLess(scp)));
	  //s->add((children[0] == vali) <= y->children[idx]);
	} else if(y->get_element_lb_index(vali)>=0)
	  if(FAILED(children[0].remove(vali))) {
	    s->fail();
	    break;
	  }
	vnxt = children[0].next(vali);
      } while(vali < vnxt);
    } else if(size == (ub-lb+1)) {
      // interval 
      if(FAILED(children[0].remove_interval(lb,ub))) { 
	s->fail(); 
      }
    } else {
      // set
      if(FAILED(children[0].remove_set(values))) {  
	s->fail(); 
      }
    }
  }   

}

void Mistral::MemberExpression::extract_variable(Solver *s) {
  Variable aux(0, 1, BOOL_VAR);
  _self = aux;
  
  _self.initialise(s, 1);
  _self = _self.get_var();
  children.add(_self);
}

void Mistral::MemberExpression::extract_predicate(Solver *s) { 

  //unsigned int i = 0, j = 0;
  if(children.size == 3) {
    // set variable
    
    // go in parallel through the domain of x and the elements of y
    
    int xi, yi, yl;
    unsigned int il=0, i=0;
    Domain dom_xi(children[0]);
    Domain::iterator it = dom_xi.begin();
    Domain::iterator end = dom_xi.end();
    SetExpression *Y = (SetExpression*)(children[1].expression);
    Variable x = children[0];
    Variable b = children[2];
    //VarArray scp;
    VarArray implies_b;
    VarArray implies_notb;

    while(it != end || i<Y->elts_var.size || il<Y->elts_lb.size) {

      xi = (it != end ? dom_xi.get_value(it) : INFTY);
      yi = (i<Y->elts_var.size ? Y->elts_var[i] : INFTY);
      yl = (il<Y->elts_lb.size ? Y->elts_lb[il] : INFTY);
  
      if(yi < yl) {
	// variable element

	if(xi < yi) {
	  // x == xi -> not(b)
	  implies_notb.add( x == xi );

	  ++it;
	} else if(yi == xi) {
	  // x == xi & yi \in Y -> b
	  // x == xi & yi \not\in Y -> not b

	  Variable x_eq_xi = (x == xi);
	  Variable yi_in_Y = Y->get_index_var(i);

	  implies_b.add( x_eq_xi && yi_in_Y );
	  implies_notb.add( x_eq_xi > yi_in_Y );

	  ++it;
	  ++i;
	} else {
	  // not important
	  ++i;
	}

      } else {
	// ground element
	if(xi < yl) {
	  // x == xi -> not(b)
	  implies_notb.add( x == xi );

	  ++it;
	} else if(yl == xi) {
	  // x == xi -> b
	  implies_b.add( x == xi );

	  ++it;
	  ++il;
	} else {
	  // not important
	  ++il;
	}

      }
    }

    if(implies_b.size) {
      for(i=0; i<implies_b.size; ++i) {
	s->add(implies_b[i] <= b);
      }
      s->add(b <= BoolSum(implies_b));
    } else {
      if(FAILED(b.set_domain(false))) {
	s->fail();
      }
    } 

    if(implies_notb.size) {
      for(i=0; i<implies_notb.size; ++i) {
	s->add(NotAnd(implies_notb[i],b));
      }
      implies_notb.add(b);
      s->add(BoolSum(implies_notb, 1, INFTY));
    } else {
      if(FAILED(b.set_domain(true))) {
	s->fail();
      } 
    }


    // VarArray scp;
    // VarArray conjunction;
    // SetExpression *y = (SetExpression*)(children[1].expression);

    // Variable empty(0,1);
    // s->add(empty);
    // scp.add(children[1]);
    // scp.add(empty);
    // s->add(Constraint(new PredicateLowerBound(scp, 1)));

    // scp.clear();
    // scp.add(children[2]);
    // scp.add(empty);
    // s->add(Constraint(new ConstraintLess(scp)));

    // //s->add(children[1] >= 1); // at least one element in the set
    // int vali, vnxt = children[0].get_first(), idx;
    // do {
    //   scp.clear();
    //   vali = vnxt;
    //   idx = y->get_element_index(vali);
    //   if(idx>=0) {
    // 	Variable x(0,1);
    // 	s->add(x);
    // 	scp.add(children[0]);
    // 	scp.add(x);
    // 	s->add(Constraint(new PredicateConstantEqual(scp,vali,1)));
	
    // 	scp.clear();

    // 	Variable b(0,1);
    // 	s->add(b);
    // 	scp.add(x);
    // 	scp.add(y->get_index_var(idx));
    // 	scp.add(b);
    // 	s->add(Constraint(new PredicateLess(scp)));

    // 	conjunction.add(b);

    // 	//s->add((children[0] == vali) <= y->children[idx]);
    //   } else if(y->get_element_lb_index(vali) < 0) {
    //   //} else {
    // 	Variable b(0,1);
    // 	s->add(b);
    // 	scp.add(children[0]);
    // 	scp.add(b);
    // 	s->add(Constraint(new PredicateConstantEqual(scp,vali,0)));

    // 	conjunction.add(b);

    // 	//children[0].remove(vali);

    // 	//s->add(children[0] != vali);
    //   }
    //   vnxt = children[0].next(vali);
    // } while(vali < vnxt);

    // Variable N(0, conjunction.size);
    // s->add(N);
    // conjunction.add(N);
    // s->add(new PredicateBoolSum(conjunction));

    // scp.clear();
    // scp.add(N);
    // scp.add(children[2]);
    // s->add(Constraint(new PredicateConstantEqual(scp,conjunction.size-1,1)));

  } else if(size == (ub-lb+1)) {
    // interval 
    
    s->add(Constraint(new PredicateIntervalMember(children,lb,ub)));
    //children[0]->set_min(lb);
    //children[0]->set_max(ub);
  } else {
    // set
    
    s->add(Constraint(new PredicateSetMember(children,values)));
    //children[0]->set_domain(values);
  }
}

const char* Mistral::MemberExpression::get_name() const {
  return "member";
}

Mistral::Variable Mistral::Member(Variable X, Variable Y) {
  Variable exp(new MemberExpression(X, Y));
  return exp;
}

Mistral::Variable Mistral::Member(Variable X, const int lo, const int up) {
  Variable exp(new MemberExpression(X, lo, up));
  return exp;
}

Mistral::Variable Mistral::Member(Variable X, const Interval I) {
  Variable exp(new MemberExpression(X, I.min, I.max));
  return exp;
}

Mistral::Variable Mistral::Member(Variable X, const BitSet& s) {
  Variable exp(new MemberExpression(X, s));
  return exp;
}

Mistral::Variable Mistral::Member(Variable X, const Vector<int>& s) {
  Variable exp(new MemberExpression(X, s));
  return exp;
}

Mistral::Variable Mistral::Member(Variable X, const std::vector<int>& s) {
  Variable exp(new MemberExpression(X, s));
  return exp;
}


Mistral::Goal::Goal(method t) : type(t) {
  lower_bound = 0;
  upper_bound = 0;
  sub_type = NONE; 
}

Mistral::Goal::Goal(method t, Variable X) : sub_type(t) {
  objective = X;
  type = OPTIMIZATION;

  // std::cout << "OBJECTIVE=" << objective << " in " << objective.get_domain() << std::endl;

  lower_bound = objective.get_min()-1; //(type == MAXIMIZATION);
  upper_bound = objective.get_max()+1; //(type == MINIMIZATION);
}

Mistral::Goal::Goal(method t, method st, Variable X) : type(t), sub_type(st) {
  objective = X;

  // std::cout << "OBJECTIVE=" << objective << " in " << objective.get_domain() << std::endl;

  lower_bound = objective.get_min()-1; //(type == MAXIMIZATION);
  upper_bound = objective.get_max()+1; //(type == MINIMIZATION);
}

Mistral::Goal::~Goal() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete goal" << std::endl;
#endif
}

bool Mistral::Goal::enforce() {
  // if(type == OPTIMIZATION) {
  //   if(sub_type == MINIMIZATION) {
  //     return FAILED(objective.set_max(upper_bound-1));
  //   } else { //if(sub_type == MAXIMIZATION) {
  //     return FAILED(objective.set_min(lower_bound+1));
  //   }
  // } else if(sub_type != NONE) {
  //   if(sub_type == MINIMIZATION) {
  //     return FAILED(objective.set_max(upper_bound));
  //   } else { //if(sub_type == MAXIMIZATION) {
  //     return FAILED(objective.set_min(lower_bound));
  //   }
  // }

  //std::cout << "enforce " << (*this) << std::endl;

  if(sub_type == MINIMIZATION) {

    //std::cout << objective << " in " << objective.get_domain() << " <= " << (upper_bound-1) << std::endl;

    return FAILED(objective.set_max(upper_bound-1));
  } else if(sub_type == MAXIMIZATION) {

    //std::cout << objective << " in " << objective.get_domain() << " >= " << (lower_bound+1) << std::endl;

    return FAILED(objective.set_min(lower_bound+1));
  }

  return false;
}


int Mistral::Goal::value() const {
  return(sub_type == MINIMIZATION ? upper_bound : lower_bound);
  // if(type == MINIMIZATION) {
  //   return upper_bound;
  // } else if(type == MAXIMIZATION) {
  //   return lower_bound;
  // } else return -1;
  // return false;
}

// int Mistral::Goal::best() const {
//   return(type = MINIMIZATION ? upper_bound : lower_bound);
//   // if(type == MINIMIZATION) {
//   //   return upper_bound;
//   // } else if(type == MAXIMIZATION) {
//   //   return lower_bound;
//   // } else return -1;
//   // return false;
// }


bool Mistral::Goal::is_optimization() const {
  return (type == OPTIMIZATION); //(type == MINIMIZATION || type == MAXIMIZATION);
}

bool Mistral::Goal::is_satisfaction() const {
  return (type == SATISFACTION);
}

bool Mistral::Goal::is_enumeration() const {
  return (type == ENUMERATION);
}

bool Mistral::Goal::has_function() const {
  return (sub_type != NONE);
}

bool Mistral::Goal::improving(const int val) const {
  bool value = false;
  if(sub_type == MINIMIZATION) value = val < upper_bound;
  else if(sub_type == MAXIMIZATION) value = val > lower_bound;
  return value;

  // return( type == SATISFACTION ? false :
  // 	  (type == MINIMIZATION ? 
  // 	   (val < upper_bound) : 
  // 	   (val > lower_bound)) );
}
    

Mistral::Outcome Mistral::Goal::notify_exhausted() {
  Outcome out = OPT;
  if(type == SATISFACTION) out = UNSAT;
  //else if(type == ENUMERATION) out = ALL;
  return out;
}

std::ostream& Mistral::Goal::display(std::ostream& os) const {
  if(type == OPTIMIZATION) {
    if(sub_type == MINIMIZATION) {
      os << "minimize " << objective << objective.get_domain();
    } else   if(sub_type == MAXIMIZATION) {
      os << "maximize " << objective << objective.get_domain();
    }
  } else if(type == ENUMERATION) {
    os << "find all solutions" ;
  } else {
    os << "find any solution" ;
    if(sub_type == MAXIMIZATION) 
      os << "(with " << objective << " >= " << lower_bound << ")" ;
    else if(sub_type == MINIMIZATION) 
      os << "(with " << objective << " <= " << upper_bound << ")" ;
  }
  return os;
}

Mistral::Outcome Mistral::Goal::notify_solution(Solver *solver) {

  //std::cout << "notify solution to objective\n";

  if(type == OPTIMIZATION) {
    if(sub_type == MINIMIZATION) {
      upper_bound = objective.get_min();

      //std::cout << "minimization -> lb = " << upper_bound << " " << objective.get_domain() << " (ub = " << lower_bound << ")" << std::endl;
    
      //std::cout << "LEVEL: " << solver->level << std::endl;
      // search the deepest level 
      int level, search_root = solver->search_root;
      do {
	level = solver->level;
	if(level == search_root) {

	  //std::cout << " search tree root reached" << std::endl;

	  lower_bound = upper_bound;
	  return OPT;
	}

	//std::cout << " backtrack one level " ;

	solver->restore(level-1);

	//std::cout << objective.get_domain() << std::endl;

      } while(upper_bound <= objective.get_min());
	
      Decision deduction(objective, Decision::UPPERBOUND, upper_bound-1);


// #ifdef _DEBUG_SEARCH
//       SolverParamaters& parameters(solver->parameters)
//       if(_DEBUG_SEARCH) {
// 	std::cout << "c";
// 	for(unsigned int k=0; k<=solver->decisions.size; ++k) std::cout << " ";
// 	std::cout << "backtrack to lvl " << solver->level << " and deduce " 
// 		  << deduction << " (upper bound)" << std::endl;
//       }
// #endif

      deduction.make();

      return UNKNOWN;

      // upper_bound = objective.get_min();
      // if(!solver->level) lower_bound = upper_bound;
    } else { //if(sub_type == MAXIMIZATION) {

      lower_bound = objective.get_max();

      // search the deepest level 
      int level, search_root = solver->search_root;
      do {
       	level = solver->level;
       	if(level == search_root) {
       	  upper_bound = lower_bound;
       	  return OPT;
       	}
       	solver->restore(level-1);
      } while(lower_bound >= objective.get_max());
	
      Decision deduction(objective, Decision::LOWERBOUND, lower_bound+1);

// #ifdef _DEBUG_SEARCH
//       SolverParamaters& parameters(solver->parameters)
//       if(_DEBUG_SEARCH) {
// 	std::cout << "c";
// 	for(unsigned int k=0; k<=solver->decisions.size; ++k) std::cout << " ";
// 	std::cout << "backtrack to lvl " << solver->level << " and deduce " 
// 		  << deduction << " (lower bound)" << std::endl;
//       }
// #endif
      
      deduction.make();
      return UNKNOWN;

      // lower_bound = objective.get_max();
      // if(!solver->level) upper_bound = lower_bound;
    }
      
    if(upper_bound == lower_bound) return OPT;
    solver->branch_right();
    return UNKNOWN; //(upper_bound == lower_bound ? OPT : UNKNOWN);
  } else if(type == ENUMERATION) {
    //solver->store_solution();

    if(!solver->level) return OPT;
    solver->branch_right();
    return UNKNOWN;
  } else if(sub_type == MAXIMIZATION) {
    lower_bound = objective.get_min()-1;
    if(!solver->level) upper_bound = lower_bound;
  } else if(sub_type == MINIMIZATION){
    upper_bound = objective.get_max()+1;
    if(!solver->level) lower_bound = upper_bound;
  }

  //std::cout << "satisfaction algorithm: new solution" << std::endl;
  return SAT;
}


Mistral::Domain::Domain(const Variable& x, const bool _open) : Variable(x) {
  if(_open) open();
  // _begin_ptr = NULL;
  // _end_ptr = NULL;
  // if(domain_type == LIST_VAR) id = 0;
  // else id = -1;
}

Mistral::Domain::~Domain() {
  close();
}

void Mistral::Domain::open() {

  if(domain_type == LIST_VAR) {
    id = 0;
    _begin_ptr = list_domain->domain.begin();
    _end_ptr = list_domain->domain.end();
  } else if(is_range()) {
    id = -1;
    _begin_ptr = (int*)(uintptr_t)(get_min()*4);
    _end_ptr = (int*)(uintptr_t)((get_max()+1)*4);
  } else {
    Solver *s = get_solver();
    int n = bitset_domain->domain.size;
    // if(id>=0) {
    //   s->iterator_space.release(id);
    // }
    id = s->iterator_space.reserve(n, _begin_ptr, _end_ptr);
    _begin_ptr[0] = bitset_domain->domain.min;
    bitset_domain->domain.values.iterate_into(n, _begin_ptr);
  }
}

void Mistral::Domain::close() {
  if(id>0) {
    Solver *s = get_solver();
    s->iterator_space.release(id);
  }
}
//   if(domain_type != LIST_VAR && !is_range()) {
//     s->release(id);
//     id = -1;
//   }
// }

Mistral::Domain::iterator Mistral::Domain::begin() {
  return _begin_ptr;
}

Mistral::Domain::iterator Mistral::Domain::end() {
  return _end_ptr;
}

int Mistral::DeltaBool::diterator[2] = {0,1};



Mistral::DomainDelta::DomainDelta() {
}

void Mistral::DomainDelta::initialise(const Variable& x) {
  domain_type = x.domain_type;
  if(domain_type == LIST_VAR)  {
    variable = (VariableImplementation*)(new DeltaList(x.list_domain));
  } else if(domain_type == BITSET_VAR) {
    variable = (VariableImplementation*)(new DeltaBitset(x.bitset_domain));
  } else if(domain_type == RANGE_VAR) {
    variable = (VariableImplementation*)(new DeltaRange(x.range_domain));
  } else if(domain_type == CONST_VAR || domain_type == EXPRESSION) {
    std::cerr << "not implemented" << std::endl;
    exit(1);
  } else {
    variable = (VariableImplementation*)(new DeltaBool(x.bool_domain, x.variable->solver));
  }
}

void Mistral::DomainDelta::cleanup() {
  if(domain_type == LIST_VAR)  {
    delete (DeltaList*)variable;
  } else if(domain_type == BITSET_VAR) {
    delete (DeltaBitset*)variable;
  } else if(domain_type == RANGE_VAR) {
    delete (DeltaRange*)variable;
  } else if(domain_type == CONST_VAR || domain_type == EXPRESSION) {
    std::cerr << "not implemented" << std::endl;
    exit(1);
  } else {
    delete (DeltaBool*)variable;
  }
}


Mistral::DomainDelta::~DomainDelta() {
  cleanup();
}

//void open();
void Mistral::DomainDelta::close() {
  if(domain_type == LIST_VAR)  {
    ((DeltaList*)variable)->close();
  } else if(domain_type == BITSET_VAR) {
    ((DeltaBitset*)variable)->close();
  } else if(domain_type == RANGE_VAR) {
    ((DeltaRange*)variable)->close();
  } else {
    ((DeltaBool*)variable)->close();
  }
}

Mistral::DomainDelta::iterator Mistral::DomainDelta::begin() {
  Mistral::DomainDelta::iterator it = NULL;
  if(domain_type == LIST_VAR)  {
    it = ((DeltaList*)variable)->begin();
  } else if(domain_type == BITSET_VAR) {
    it = ((DeltaBitset*)variable)->begin();
  } else if(domain_type == RANGE_VAR) {
    it = ((DeltaRange*)variable)->begin();
  } else {
    it = ((DeltaBool*)variable)->begin();
  }
  return it;
}

Mistral::DomainDelta::iterator Mistral::DomainDelta::end() {
  Mistral::DomainDelta::iterator it = NULL;
  if(domain_type == LIST_VAR)  {
    it = ((DeltaList*)variable)->end();
  } else if(domain_type == BITSET_VAR) {
    it = ((DeltaBitset*)variable)->end();
  } else if(domain_type == RANGE_VAR) {
    it = ((DeltaRange*)variable)->end();
  } else {
    it = ((DeltaBool*)variable)->end();
  }
  return it;
}

Mistral::Literal Mistral::literal(Variable x, const int val) {
  return (x.id()*2+val);
}

Mistral::Literal Mistral::literal(Variable x) {
  return (x.id()*2+x.get_value());
}


// bool Mistral::Decision::make() {

//   //std::cout << _data_ << std::endl;

//   //std::cout << var << " in " << var.get_domain() << " : " << value() << std::endl;
  
//   //if(_data_ == -1) return propagateRelation();
//   switch(type()) {
//   case REMOVAL:    return !FAILED(var.remove(value()));
//   case ASSIGNMENT: return !FAILED(var.set_domain(value()));
//   case LOWERBOUND: return !FAILED(var.set_min(value()+1));
//   case UPPERBOUND: return !FAILED(var.set_max(value()));
//   }
//   return true;
// }
Mistral::Variable Mistral::AtMostSeqCardNaiveReason(Vector< Variable >& args, const int d, const int p, const int q) {
  Variable exp(new AtMostSeqCardExpressionNaiveReason(args,d,p,q));
  return exp;
}


Mistral::Variable Mistral::AtMostSeqCardSimplifiedReason(Vector< Variable >& args, const int d, const int p, const int q) {
  Variable exp(new AtMostSeqCardExpressionSimplifiedReason(args,d,p,q));
  return exp;
}

//left
Mistral::Variable Mistral::AtMostSeqCardLeftExplanationReason(Vector< Variable >& args, const int d, const int p, const int q) {
  Variable exp(new AtMostSeqCardExpressionLeftExplanationReason(args,d,p,q));
  return exp;
}



Mistral::AtMostSeqCardExpressionNaiveReason::AtMostSeqCardExpressionNaiveReason(
										Vector<Variable>& args, const int d, const int p, const int q)
  : AtMostSeqCardExpression(args, d, p, q)
{
}

void Mistral::AtMostSeqCardExpressionNaiveReason::extract_constraint(Solver* s)
{
  s->add(Constraint(new ConstraintNaiveMultiAtMostSeqCard(children, _k, _d, _p, _q)));
}
const char* Mistral::AtMostSeqCardExpressionNaiveReason::get_name() const {
  return "naive amsc";
}


Mistral::AtMostSeqCardExpressionSimplifiedReason::AtMostSeqCardExpressionSimplifiedReason(
											  Vector<Variable>& args, const int d, const int p, const int q)
  : AtMostSeqCardExpression(args, d, p, q)
{
}

void Mistral::AtMostSeqCardExpressionSimplifiedReason::extract_constraint(Solver* s)
{
  s->add(Constraint(new ConstraintSimplifiedExplanationMultiAtMostSeqCard(children, _k, _d, _p, _q)));
}
const char* Mistral::AtMostSeqCardExpressionSimplifiedReason::get_name() const {
  return "simplified reason amsc";
}

//left
Mistral::AtMostSeqCardExpressionLeftExplanationReason::AtMostSeqCardExpressionLeftExplanationReason(
												    Vector<Variable>& args, const int d, const int p, const int q)
  : AtMostSeqCardExpression(args, d, p, q)
{
}

void Mistral::AtMostSeqCardExpressionLeftExplanationReason::extract_constraint(Solver* s)
{
  s->add(Constraint(new ConstraintLeftExplanationMultiAtMostSeqCard(children, _k, _d, _p, _q)));
}
const char* Mistral::AtMostSeqCardExpressionLeftExplanationReason::get_name() const {
  return "left reason amsc";
}







Mistral::OccExpression::OccExpression(Vector< Variable >& args, const int lo, const int up)
  : Expression(args) {
  current_occ = 0;
  lower_bound = lo;
  upper_bound = up;
  //num_args = args.size;
}

Mistral::OccExpression::OccExpression(VarArray& args, const int lo, const int up)
  : Expression(args) {
  current_occ = 0;
  lower_bound = lo;
  upper_bound = up;
  //num_args = args.size;
}

Mistral::OccExpression::~OccExpression() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete min expression" << std::endl;
#endif
}

void Mistral::OccExpression::extract_constraint(Solver *s) {
  encode();
  if(lower_bound > -INFTY || upper_bound < INFTY) {    
    lower_bound -= current_occ;
    upper_bound -= current_occ;
    if(lower_bound > (int)(scope.size) || upper_bound < 0) {
      s->fail();
    } else {
      s->add(Constraint(new ConstraintBoolSumInterval(scope,lower_bound,upper_bound))); 
    }
  } 
}

void Mistral::OccExpression::extract_variable(Solver *s) {
  encode();
  
  int arity = scope.size;
  if(lower_bound < current_occ) lower_bound = current_occ;
  if(upper_bound < arity+current_occ) upper_bound = arity+current_occ;

  Variable aux(lower_bound, upper_bound, DYN_VAR);
  _self = aux;

  _self.initialise(s, 1);
  _self = _self.get_var();
  children.add(_self);
  scope.add(_self);
}

void Mistral::OccExpression::extract_predicate(Solver* s) {
  s->add(Constraint(new PredicateBoolSum(scope,current_occ)));
}

const char* Mistral::OccExpression::get_name() const {
  return "occ";
}


Mistral::ValOccExpression::ValOccExpression(Vector< Variable >& args, const int val, const int lo, const int up)
  : OccExpression(args, lo, up) {
  value = val;
} 

void Mistral::ValOccExpression::encode() {
  for(unsigned int i=0; i<children.size; ++i) {
    if(children[i].equal(value))
      ++current_occ;
    else if(children[i].contain(value)) {
      scope.add( children[i] == value );
    }
  }
}

Mistral::Variable Mistral::Occurrence(Vector<Variable>& X, const int v, const int lo, const int up) {
  Variable exp = Variable( new ValOccExpression(X, v, lo, up) );
  return exp;
}

Mistral::Variable Mistral::Occurrence(VarArray& X, const int v, const int lo, const int up) {
  Variable exp = Variable( new ValOccExpression(X, v, lo, up) );
  return exp;
}



Mistral::VarOccExpression::VarOccExpression(Vector< Variable >& args, Variable x, const int lo, const int up)
  : OccExpression(args, lo, up) {
  X = x;
} 

void Mistral::VarOccExpression::encode() {
  for(unsigned int i=0; i<children.size; ++i) {
    if(X.is_ground() && children[i].equal(X.get_value()))
      ++current_occ;
    else if(children[i].intersect(X)) {
      scope.add( children[i] == X );
      //children.add( )
    }
  }
}

Mistral::Variable Mistral::Occurrence(Vector<Variable>& X, Variable x, const int lo, const int up) {
  Variable exp = Variable( new VarOccExpression(X, x, lo, up) );
  return exp;
}

Mistral::Variable Mistral::Occurrence(VarArray& X, Variable x, const int lo, const int up) {
  Variable exp = Variable( new VarOccExpression(X, x, lo, up) );
  return exp;
}



Mistral::SetOccExpression::SetOccExpression(Vector< Variable >& args, const BitSet& s, const int lo, const int up)
  : OccExpression(args, lo, up) {
  S.initialise(s);
} 

Mistral::SetOccExpression::SetOccExpression(Vector< Variable >& args, const Vector<int>& s, const int lo, const int up)
  : OccExpression(args, lo, up) {
  int lb = +INFTY;
  int ub = -INFTY;
  for(unsigned int i=0; i<s.size; ++i) {
    if(s[i] < lb) lb = s[i];
    if(s[i] > ub) ub = s[i];
  }
  S.initialise(lb, ub, BitSet::empt);
  for(unsigned int i=0; i<s.size; ++i) {
    S.add(s[i]);
  }
}

Mistral::SetOccExpression::SetOccExpression(Vector< Variable >& args, const std::vector<int>& s, const int lo, const int up)
  : OccExpression(args, lo, up) {
  int lb = +INFTY;
  int ub = -INFTY;
  for(unsigned int i=0; i<s.size(); ++i) {
    if(s[i] < lb) lb = s[i];
    if(s[i] > ub) ub = s[i];
  }
  S.initialise(lb, ub, BitSet::empt);
  for(unsigned int i=0; i<s.size(); ++i) {
    S.add(s[i]);
  }
}

void Mistral::SetOccExpression::encode() {
  for(unsigned int i=0; i<children.size; ++i) {
    if(children[i].included(S))
      ++current_occ;
    else if(children[i].intersect(S)) {
      scope.add( Member(children[i], S) );
    }
  }
}

Mistral::Variable Mistral::Occurrence(Vector<Variable>& X, const BitSet& s, const int lo, const int up) {
  Variable exp = Variable( new SetOccExpression(X, s, lo, up) );
  return exp;
}

Mistral::Variable Mistral::Occurrence(VarArray& X, const BitSet& s, const int lo, const int up) {
  Variable exp = Variable( new SetOccExpression(X, s, lo, up) );
  return exp;
}

Mistral::Variable Mistral::Occurrence(Vector<Variable>& X, const Vector<int>& s, const int lo, const int up) {
  Variable exp = Variable( new SetOccExpression(X, s, lo, up) );
  return exp;
}

Mistral::Variable Mistral::Occurrence(VarArray& X, const Vector<int>& s, const int lo, const int up) {
  Variable exp = Variable( new SetOccExpression(X, s, lo, up) );
  return exp;
}

Mistral::Variable Mistral::Occurrence(Vector<Variable>& X, const std::vector<int>& s, const int lo, const int up) {
  Variable exp = Variable( new SetOccExpression(X, s, lo, up) );
  return exp;
}

Mistral::Variable Mistral::Occurrence(VarArray& X, const std::vector<int>& s, const int lo, const int up) {
  Variable exp = Variable( new SetOccExpression(X, s, lo, up) );
  return exp;
}



Mistral::IntOccExpression::IntOccExpression(Vector< Variable >& args, const Interval i, const int lo, const int up)
  : OccExpression(args, lo, up) {
  I = i;
} 

void Mistral::IntOccExpression::encode() {
  for(unsigned int i=0; i<children.size; ++i) {
    if(children[i].included(I))
      ++current_occ;
    else if(children[i].intersect(I)) {
      scope.add( Member(children[i], I) );
    }
  }
}

Mistral::Variable Mistral::Occurrence(Vector<Variable>& X, const Interval i, const int lo, const int up) {
  Variable exp = Variable( new IntOccExpression(X, i, lo, up) );
  return exp;
}

Mistral::Variable Mistral::Occurrence(VarArray& X, const Interval i, const int lo, const int up) {
  Variable exp = Variable( new IntOccExpression(X, i, lo, up) );
  return exp;
}

// Mistral::Variable Mistral::Occurrence(Vector<Variable>& X, const int l, const int u, const int lo=-INFTY, const int up=INFTY) {
//   Variable exp = Variable( new IntOccExpression(X, Interval(l,u), lo, up) );
//   return exp;
// }

// Mistral::Variable Mistral::Occurrence(VarArray& X, const int l, const int u, const int lo=-INFTY, const int up=INFTY) {
//   Variable exp = Variable( new IntOccExpression(X, Interval(l, u), lo, up) );
//   return exp;
// }


/*
  void Mistral::OccExpression::extract_predicate(Solver *s) {
  s->add(new PredicateOcc(children));
  }

  Mistral::Variable Mistral::Occ(Vector<Variable>& X) {
  bool equality = true;
  for(unsigned int i=1; i<X.size && equality; ++i) equality = (X[i-1].operator_equal(X[i]));
  Variable exp;
  if(equality) exp = X[0];
  else exp = Variable( new OccExpression(X) );
  return exp;
  }

  Mistral::Variable Mistral::Occ(VarArray& X) {
  bool equality = true;
  for(unsigned int i=1; i<X.size && equality; ++i) equality = (X[i-1].operator_equal(X[i]));
  Variable exp;
  if(equality) exp = X[0];
  else exp = Variable( new OccExpression(X) );
  return exp;
  // Variable exp( new OccExpression(X) );
  // return exp;
  }
*/
