
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

  The author can be contacted electronically at 
  ehebrard@cse.unsw.edu.au.
*/


#include <mistral_solver.hpp>
#include <mistral_backtrack.hpp>
#include <mistral_constraint.hpp>


//void _save_() { if(first_change) { solver->save(this); first_change=false; trail_.push(value); } }

//void Mistral::Environment::save(Mistral::Reversible *r) { saved_objs.add(r); }

Mistral::Constraint::Constraint() {
  propagator = NULL;
  data = 0;
}

Mistral::Constraint::Constraint(ConstraintImplementation* p) { 
  propagator = p; 
  propagator->ConstraintImplementation::initialise();
  data = propagator->type; 
}

void Mistral::Constraint::re_link_to(Mistral::Trigger* t) {
  propagator->on[index()] = t;
}

void Mistral::Constraint::set_rank(const int idx) {
  int rank = (data&CTYPE);
  propagator->index[rank] = idx;
}

void Mistral::Constraint::set_index(const int idx) {
  data &= ITYPE;
  data |= idx;
}

void Mistral::Constraint::set_id(const int idx) {
  propagator->id = idx;
}

int Mistral::Constraint::id() const {
  return propagator->id;
}

std::string Mistral::Constraint::symbol() const {
  return propagator->name();
}

bool Mistral::Constraint::explained() const { 
  return propagator->explained();
}

void Mistral::Constraint::consolidate() {
  propagator->consolidate();
}

void Mistral::Constraint::consolidate_var() {
  propagator->consolidate_var(index());
}

bool Mistral::Constraint::rewritable() {
  return propagator->rewritable(); 
}

bool Mistral::Constraint::simple_rewritable() {
  return propagator->simple_rewritable();
}

bool Mistral::Constraint::absorb_negation(const int i) {
  return propagator->absorb_negation(i); 
}

Mistral::Constraint Mistral::Constraint::get_negation(const int i, Variable x) {
  return propagator->get_negation(i, x); 
}

void Mistral::Constraint::initialise_activity(double *lvact, double *vact, double norm) {
  propagator->initialise_activity(lvact, vact, norm);
}

void Mistral::Constraint::relax() {
  if(binary()) 
    ((BinaryConstraint*)propagator)->relax(); 
  else if(ternary())
    ((TernaryConstraint*)propagator)->relax(); 
  else
    ((GlobalConstraint*)propagator)->relax(); 
}

void Mistral::Constraint::initialise(Solver *s) {
  if(propagator) {
    //if(propagator->solver != s) {
    propagator->solver = s;
    propagator->initialise_vars(s);
    propagator->initialise();
    //}
  } else {
    std::cerr << "Should raise a contradiction!" << std::endl;
  }
}

int Mistral::Constraint::priority() const {
  return (global() ? ((GlobalConstraint*)propagator)->priority : 2);
}

void Mistral::Constraint::post(Solver* solver) { 
  if(propagator) {
    propagator->initial_post(solver); 
  } else {
    std::cerr << "Should raise a contradiction!" << std::endl;
  }
}


void Mistral::Constraint::notify_assignment() {
  if(binary()) {
    ((BinaryConstraint*)propagator)->notify_assignment(index()); 
  } else if(ternary())
    ((TernaryConstraint*)propagator)->notify_assignment(index()); 
  else
    ((GlobalConstraint*)propagator)->notify_assignment(index()); 
}

void Mistral::Constraint::awaken() { 
  if(binary()) {
    //std::cout << "in awaken" << std::endl;
    ((BinaryConstraint*)propagator)->post(); 
  } else if(ternary())
    ((TernaryConstraint*)propagator)->post(); 
  else
    ((GlobalConstraint*)propagator)->post(); 
}

void Mistral::Constraint::trigger() { 
  if(binary())  {
    //((Solver*)solver)->active_constraints.trigger((BinaryConstraint*)propagator);
    ((BinaryConstraint*)propagator)->trigger(); 
  } else if(ternary()) {
    //((Solver*)solver)->active_constraints.trigger((TernaryConstraint*)propagator);
    ((TernaryConstraint*)propagator)->trigger(); 
  } else {
    //((Solver*)solver)->active_constraints.trigger((GlobalConstraint*)propagator);
    ((GlobalConstraint*)propagator)->trigger(); 
  }
}

int Mistral::Constraint::num_active() const { 
  int n;
  if(binary()) 
    n = size_byte[((BinaryConstraint*)propagator)->active]; 
  else if(ternary())
    n = size_byte[((TernaryConstraint*)propagator)->active]; 
  else
    n = ((GlobalConstraint*)propagator)->active.size; 
  return n;
}


Mistral::ConstraintImplementation* Mistral::Constraint::freeze() {
  if(global()) ((GlobalConstraint*)propagator)->freeze();
  return (idempotent() ? propagator : NULL);
}

Mistral::ConstraintImplementation* Mistral::Constraint::defrost() {
  if(global()) ((GlobalConstraint*)propagator)->defrost();
  return NULL;
}

Mistral::Variable* Mistral::Constraint::get_scope() { 
  return propagator->_scope.stack_;
}

void Mistral::Constraint::set_scope(const int i, Variable x) {
  propagator->set_scope(i, x);
  // propagator->_scope.stack_[i] = x;
  // propagator->on[i] = ;
  if(binary()) 
    ((BinaryConstraint*)propagator)->scope[i] = x; 
  else if(ternary())
    ((TernaryConstraint*)propagator)->scope[i] = x; 
  else
    ((GlobalConstraint*)propagator)->scope[i] = x; 
}

int Mistral::Constraint::arity() const {
  return propagator->_scope.size;
} 

int Mistral::Constraint::get_active(const int i) const { 
  int n=0;
  
  if(binary()) {
    if(i) n = 1; // the second element must be 1
    else n = (((BinaryConstraint*)propagator)->active >> 1); 
    // the first element is the min
  }
  else if(ternary()) {
    if(i==2) n = 2; // the third element must be 2
    else if(i) { // the second element can be either 1 or 2
      n = 1+(((BinaryConstraint*)propagator)->active & 1);
    } else { // the first element 
      int a = ((TernaryConstraint*)propagator)->active; 
      if(!(a&1)) n = (4-(a&2))/2;
    }
  }
  else
    n = ((GlobalConstraint*)propagator)->active[i]; 
  
  return n;
}

int Mistral::Constraint::check(int* sol) {
  return propagator->check(sol);
}

bool Mistral::Constraint::find_support(const int var, const int val) {
  bool exist=false;
  
  if(binary()) {
    exist = ((BinaryConstraint*)propagator)->find_support(var, val);
  } else if(ternary()) {
    exist = ((TernaryConstraint*)propagator)->find_support(var, val);
  } else {
    exist = (((GlobalConstraint*)propagator)->first_support(var, val) ||
	     ((GlobalConstraint*)propagator)->find_support(var, val));
  }

  return exist;
}

bool Mistral::Constraint::find_bound_support(const int var, const int val) {
  bool exist=false;
  
  if(binary()) {
    exist = ((BinaryConstraint*)propagator)->find_bound_support(var, val);
  } else if(ternary()) {
    exist = ((TernaryConstraint*)propagator)->find_bound_support(var, val);
  } else {
    exist = (((GlobalConstraint*)propagator)->first_support(var, val) ||
	     ((GlobalConstraint*)propagator)->find_bound_support(var, val));
  }

  return exist;
}

Mistral::PropagationOutcome Mistral::Constraint::rewrite() {
  return propagator->rewrite();
}

Mistral::PropagationOutcome Mistral::Constraint::propagate() {
  return propagator->propagate();
}

Mistral::PropagationOutcome Mistral::Constraint::propagate(const Event evt) {
  return propagator->propagate(index(), evt);
}


Mistral::PropagationOutcome Mistral::Constraint::checker_propagate() {
  return propagator->checker_propagate();
}

Mistral::PropagationOutcome Mistral::Constraint::checker_propagate(const Event evt) {
  return propagator->checker_propagate(index(), evt);
}

// void Mistral::Constraint::set_idempotent(const bool flag) {
//   std::cout << "SET IDEMPOTENT: to " << flag << " (was " << idempotent() << ")\n";
//   if(idempotent() != flag) {
//     data ^= IDEMPOTENT;
//     if(propagator)
//       for(unsigned int i=0; i<propagator->on.size; ++i) {
//      	if(propagator->self[i].data & IDEMPOTENT)
// 	  propagator->self[i].data ^= IDEMPOTENT;
//       }
//     std::cout << "SET IDEMPOTENT: " << propagator->on.size << std::endl;
//   }
// }

Mistral::PropagationOutcome Mistral::Constraint::bound_checker_propagate() {
  
  //std::cout << "backtrack.cpp: bound checker propagate()" << std::endl; 

  PropagationOutcome wiped_idx = propagator->bound_checker_propagate();

  //std::cout << wiped_idx << std::endl;

  return wiped_idx;
}

Mistral::PropagationOutcome Mistral::Constraint::bound_checker_propagate(const Event evt) {

  //std::cout << "backtrack.cpp: bound checker propagate(evt)" << std::endl; 

  PropagationOutcome wiped_idx = propagator->bound_checker_propagate(index(), evt);

  //std::cout << wiped_idx << std::endl;

  return wiped_idx;
}

void Mistral::Constraint::restore() {
  //unsigned int mytype = index();

  //if(id() == 36)
  //std::cout << "Restore " << (*this) << " (" << (int*)(propagator) << ")" << std::endl;

  if(binary()) {

    // print_bitset(((BinaryConstraint*)propagator)->active, 0, std::cout);
    // std::cout << " -> "; 
    
    //if(mytype&4) ((BinaryConstraint*)propagator)->active = 3;
    //if(mytype&2) ((BinaryConstraint*)propagator)->un_post();
    //if(mytype&1) ((BinaryConstraint*)propagator)->un_relax();
    //((BinaryConstraint*)propagator)->restore(index());
    ((BinaryConstraint*)propagator)->restore(data);

    // print_bitset(((BinaryConstraint*)propagator)->active, 0, std::cout);
    // std::cout << std::endl;

  } else if(ternary()) {
    //if(mytype&4) ((TernaryConstraint*)propagator)->re_activate(); 
    //if(mytype&2) ((TernaryConstraint*)propagator)->un_post(); 
    //if(mytype&1) ((TernaryConstraint*)propagator)->un_relax();
    ((TernaryConstraint*)propagator)->restore(data);
  } else {



#ifdef _DEBUG_BACKTRACK
    GlobalConstraint *c = (GlobalConstraint*)propagator;
    int id = c->id;
    if(_DEBUG_BACKTRACK) {
      std::cout << "c ";
      int lvl=c->solver->level;
      while(--lvl>=0) std::cout << " ";
      std::cout << "[" << std::setw(4) << c->id << "](" << c->name() << "): restore" << std::endl;
      std::cout << "c ";
      lvl=c->solver->level;
      while(--lvl>=0) std::cout << " ";
      std::cout << "[" << std::setw(4) << c->id << "](" << c->name() << "): reset active(?): " ;
      c->print_active();
      std::cout << std::endl;

      std::cout << "posted? " << (data&POSTED) << " / relaxed? " << (data&RELAXED) << std::endl;
  }

#endif




    //if(mytype&4) ((GlobalConstraint*)propagator)->re_activate(); 
    if(data&POSTED) ((GlobalConstraint*)propagator)->un_post(); 
    if(data&RELAXED) ((GlobalConstraint*)propagator)->un_relax();
    //((GlobalConstraint*)propagator)->restore(data);


  }
}

int Mistral::Constraint::get_backtrack_level() {
  if(global()) return ((GlobalConstraint*)propagator)->get_backtrack_level();
  else return propagator->solver->level-1;
}

Mistral::Decision Mistral::Constraint::get_decision() {
 if(global()) return ((GlobalConstraint*)propagator)->get_decision();
 else {
   Decision dec = ((Solver*)(propagator->solver))->decisions.back(0); 
   dec.invert();
   return dec;
 }
 //return propagator->solver->decisions.back(0);
}

std::ostream& Mistral::Constraint::display(std::ostream& os) const {
  propagator->display(os);
  // os 
  //   //<< ":" 
  //   << "[" << id() << "]";
  return os;
}

bool Mistral::Constraint::is_active() {
  return propagator->is_active();
}

void Mistral::Constraint::check_active() {
  if(global()) {
    ((GlobalConstraint*)propagator)->check_active();
  } else if(binary()) {
    ((BinaryConstraint*)propagator)->check_active();
  } else if(ternary()) {
    ((TernaryConstraint*)propagator)->check_active();
  }
}

int Mistral::Constraint::rank() {
  int rank = (data&CTYPE);
  return propagator->index[rank];
}

void Mistral::Constraint::weight_conflict(double unit, Vector<double>& weights) {
  if(!empty()) { // if the backtrack comes from the objective function, there is no culprit
    if(global()) {
      ((GlobalConstraint*)propagator)->weight_conflict(unit, weights);
    } else {
      Variable *scope = get_scope();
      int idx;
      int i = arity();
      while(i--) {
	idx = scope[i].id();
	if(idx>=0) { // this is for constants (which hade id -1)
	  weights[idx] += unit
#ifdef _DIV_ARITY
	    / (binary() ? 2.0 : 3.0)
#endif
	    ;
	}
      }   
    }
  }
}


void Mistral::Environment::_restore_() {
  
  unsigned int previous_level;
  
  previous_level = trail_.pop();
  while( saved_ints.size > previous_level ) 
    saved_ints.pop()->restore();
  
  previous_level = trail_.pop();
  while( saved_lists.size > previous_level ) 
    saved_lists.pop()->restore();
  
  previous_level = trail_.pop();
  while( saved_bools.size > previous_level ) 
    *(saved_bools.pop()) = 3;
  
  previous_level = trail_.pop();
  previous_level = trail_.pop();
  
  --level;
  
}


std::ostream& Mistral::operator<< (std::ostream& os, const Mistral::ReversibleBool& x) {
  return x.display(os);
}
std::ostream& Mistral::operator<< (std::ostream& os, const Mistral::ReversibleBool* x) {
  return x->display(os);
}

std::ostream& Mistral::operator<< (std::ostream& os, const Mistral::VarEvent& x) {
  return x.display(os);
}
std::ostream& Mistral::operator<< (std::ostream& os, const Mistral::VarEvent* x) {
  return x->display(os);
}


    /*!@name Backtrack method*/
    //@{    
// int Mistral::ReversibleSet::get_reduction() const {
//   return (trail_.back() == env->level ? trail_.back(2) - size : 0);
// }
int Mistral::ReversibleSet::get_reduction() const {
  return trail_.back(2) - size;
}

void Mistral::ReversibleSet::restore() { 
  trail_.pop(); size = trail_.pop(); 
} 

void Mistral::ReversibleSet::save() { 
  
  //std::cout << trail_.size << " " << env << std::endl;
  
  if(trail_.back() != env->level) {
    trail_.add(size);
    trail_.add(env->level);
    env->save(this);
  }
}
//@}

    /*!@name Manipulation*/
    //@{  
    // it's either 'add's...
     void Mistral::ReversibleSet::reversible_add(const int elt) {
      save();
      add(elt);
    }

    // ...or remove, but not both!!
     void Mistral::ReversibleSet::reversible_remove(const int elt) {
      save();
      remove(elt);
    }

void Mistral::ReversibleSet::reversible_add(const Vector<int>& elts) {
  save();
  Vector<int>::iterator stop = elts.end();
  for(Vector<int>::iterator elt = elts.begin(); elt!=stop; ++elt) {
    add(*elt);
  }
}


void Mistral::ReversibleSet::reversible_remove(const Vector<int>& elts) {
  save();
  Vector<int>::iterator stop = elts.end();
  for(Vector<int>::iterator elt = elts.begin(); elt!=stop; ++elt) {
    if(contain(*elt))
      remove(*elt);
  }
}


     void Mistral::ReversibleSet::reversible_set_to(const int elt) {
      save();
      set_to(elt);
    }

     int Mistral::ReversibleSet::reversible_pop()
    {
      save();
      return pop();
    }

     int Mistral::ReversibleSet::reversible_pop_head()
    {
      save();
      return pop_head();
    }
    //@}    
