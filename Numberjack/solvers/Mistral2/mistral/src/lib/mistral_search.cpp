 
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

#include <limits>

#include <mistral_search.hpp>

#include <mistral_sat.hpp>
#include <mistral_variable.hpp>



Mistral::ConsolidateListener::ConsolidateListener(Mistral::Solver *s) 
  : VariableListener(), ConstraintListener(), solver(s)
{
	
  sequence = &(solver->sequence);
  constraints.initialise(solver->variables.size);


  // int n = solver->variables.size;
  // for(int i=0; i<n; ++i) {
  //   Vector< Constraint > neighborhood;
  //   neighborhood.initialise(solver->variables[i].get_degree());
  //   for(int k=0; k<3; ++k) {
  //     for(int j=solver->constraint_graph[i].on[k].size-1; j>=0; --j)
  //    	neighborhood.add(solver->constraint_graph[i].on[k][j]);
  //   }
  //   constraints.add(neighborhood);
  // }

  //std::cout << "Variables size: " << solver->variables.size << std::endl;

  Constraint c;
  Variable *scope;
  int i, j, arity, n=solver->variables.size;

  for(i=0; i<n; ++i) {
    constraints[i].initialise(solver->variables[i].get_degree());
  }

  n = solver->constraints.size;
  for(i=0; i<n; ++i) {
    c = solver->constraints[i];
    scope = c.get_scope();
    arity = c.arity();
    for(j=0; j<arity; ++j) {
      if(!scope[j].is_ground()) {
	//std::cout << " " << scope[j] << std::endl;
	c.set_index(j);
	constraints[scope[j].id()].add(c);
      }
    }
  }

  id_obj = -1;
  if(solver->objective) {
    id_obj = solver->objective->objective.id();
  }

  //std::cout << constraints[66] << std::endl;

  solver->add((VariableListener*)this);
  solver->add((ConstraintListener*)this);
}


Mistral::ConsolidateListener::~ConsolidateListener() {
  //std::cout << "in delete consolidate manager" << std::endl;
}

void Mistral::ConsolidateListener::notify_add_var() {
  Variable x = solver->variables.back();
  while((int)(constraints.size)<x.id()) {
    Vector< Constraint > neighborhood;
    constraints.add(neighborhood);
  }
  Vector< Constraint > neighborhood;
  neighborhood.initialise(x.get_degree());
  for(int k=0; k<3; ++k) {
    for(int j=solver->constraint_graph[x.id()].on[k].size-1; j>=0; --j)
      neighborhood.add(solver->constraint_graph[x.id()].on[k][j]);
  }
  constraints.add(neighborhood);
}

void Mistral::ConsolidateListener::notify_post (Constraint c) {};

void Mistral::ConsolidateListener::notify_relax(Constraint c) {};

void Mistral::ConsolidateListener::notify_add_con(Constraint c) {

  Variable *scope = c.get_scope();
  int arity = c.arity();
  for(int i=0; i<arity; ++i) {
		if(!scope[i].is_ground()) {
			constraints[scope[i].id()].add(c);
    	constraints[scope[i].id()].back().set_index(i);
		}
  }
}

void Mistral::ConsolidateListener::notify_change(const int idx) {

  //std::cout << "BEG REACT TO CHANGE ON " << solver->variables[idx] << std::endl;

  Variable X = solver->variables[idx];
  int ids = sequence->index(idx);
  if(ids>=0) sequence->list_[ids] = X;

  //if(idx==66) std::cout << constraints[idx] << std::endl;

  for(int i=constraints[idx].size; --i>=0;) {
    //std::cout << constraints[idx][i] << std::endl;
    constraints[idx][i].consolidate_var();
    //std::cout << constraints[idx][i] << std::endl;
  }

  //std::cout << "IDX=" << idx << " IDO=" << id_obj << std::endl;

  if(idx==id_obj) {
    solver->objective->objective = X;
  }

  //std::cout << "END REACT TO CHANGE ON " << solver->variables[idx] << std::endl;
}


Mistral::HeuristicPoolManager::HeuristicPoolManager(Solver *s) : solver(s) {// }

  //std::cout << " c add restart listener" << std::endl;
  
  heu_index = 1;
  solver->add((RestartListener*)this);
}

Mistral::HeuristicPoolManager::~HeuristicPoolManager() {// }
  for(unsigned int i=0; i<pool.size; ++i) {
    if(solver->heuristic != pool[i]) {
      delete [] pool[i];
    }
  }
  solver->remove((RestartListener*)this);
}


void Mistral::HeuristicPoolManager::notify_restart(const double prog) {
  //std::cout << " c notify restart (3): " << solver->statistics.num_restarts << std::endl;
  if(heu_index<pool.size) {
    if(prog > 0.0) {
      counter = threshold;
    } else if(--counter <= 0) {
      counter = threshold;
      std::cout << " c switch heuristic!\n";
      solver->heuristic = pool[heu_index++];
    }
  }
  //std::cout << " c " << counter << std::endl;
}


// //Mistral::LiteralActivityManager::LiteralActivityManager(Solver *s, void *a) 
// Mistral::LiteralActivityManager::LiteralActivityManager(Solver *s) 
//   : solver(s) {

//   lit_activity = solver->lit_activity.stack_;
//   var_activity = solver->var_activity.stack_;
//   n_vars = solver->variables.size;


//   //   n_vars = solver->base->scope.size;

//   // if(solver->base) {
//   //   lit_activity = solver->base->lit_activity.stack_;
//   //   var_activity = solver->base->var_activity.stack_;
//   //   n_vars = solver->base->scope.size;
//   // } else {
//   //   n_vars = solver->variables.size;
//   //   lit_activity = new double[2*n_vars];
//   //   var_activity = new double[n_vars];
//   //   std::fill(lit_activity, lit_activity+2*n_vars, 0.012);
//   //   std::fill(var_activity, var_activity+n_vars, 0.024);
//   // }


//   //double activity_increment = parameters.activity_increment / (1 << clause.size);
  
//   // if(activity_increment > 0.0) {
//   //   int i=clause.size;
//    //   while(i--) {

//    //     //std::cout << clause << " " << activity_increment << std::endl;

//    //     lit_activity[clause[i]] += activity_increment;
//    //     var_activity[UNSIGNED(clause[i])] += activity_increment;
//    //   }
//    // }

   
//   decay = solver->parameters.activity_decay;
//   solver->add((DecisionListener*)this);
// }

// Mistral::LiteralActivityManager::~LiteralActivityManager() {
//   solver->remove((DecisionListener*)this);
//   // if(!solver->base) {
//   //   delete [] lit_activity;
//   //   delete [] var_activity;
//   // }
// }

// double *Mistral::LiteralActivityManager::get_weight() { return var_activity; }     

// void Mistral::LiteralActivityManager::notify_decision() {
//   int i=n_vars;
//   while(i--) {
//     //std::cout << i << " " << var_activity[i] << " -> ";
//     var_activity[i] *= decay;
//     //std::cout << var_activity[i] << std::endl;
//   }    
//   i=2*n_vars;
//   while(i--) lit_activity[i] *= decay;
// }    

Mistral::RestartPolicy::RestartPolicy(const unsigned int b) {
  base = b;
}

Mistral::RestartPolicy::~RestartPolicy() {
}

Mistral::NoRestart::NoRestart() 
  : RestartPolicy(-1)
{
}

Mistral::NoRestart::~NoRestart() {}

Mistral::Geometric::Geometric(const unsigned int b, const double f) 
  : RestartPolicy(b)
{
  increment = b;
  factor = f;
}

Mistral::Geometric::~Geometric() {}

Mistral::Luby::Luby(const unsigned int b) 
  : RestartPolicy(b)
{
  iteration = 0;
}

Mistral::Luby::~Luby() {}

Mistral::NoOrder::NoOrder(Solver *s) 
  : solver(s) {}

Mistral::NoOrder::~NoOrder() {}

Mistral::Variable Mistral::NoOrder::select() {
  return solver->sequence.back();
}

// Mistral::Lexicographic::Lexicographic(Solver *s) 
//   : solver(s) {
//   index.initialise(s->variables.size);
//   solver->add(this);
//   last.initialise(0,s);
// }

Mistral::Lexicographic::Lexicographic(Solver *s) 
  : solver(s) {
  index.initialise(s->variables.size);

  int n = solver->variables.size;
  std::fill(index.stack_, index.stack_+n, -1);

  //std::cout << n << std::endl;
 
  solver->add(this);
  last.initialise(s,0);
}

void Mistral::Lexicographic::initialise(VarStack< Variable, ReversibleNum<int> >& seq) {
  // int n = solver->variables.size;
  // std::fill(index.stack_, index.stack_+n, -1);
  //for(int i=0; i<seq.size; ++i) {
  if(order.empty())
    for(int i=seq.size; --i>=0;) {
      index[seq[i].id()] = order.size;
      order.add(seq[i]);
    }
}

//void Mistral::Lexicographic::initialise(Solver *s, void *a) {
void Mistral::Lexicographic::initialise(Solver *s) {
  solver = s;
  index.initialise(s->variables.size);

  int n = solver->variables.size;
  std::fill(index.stack_, index.stack_+n, -1);

  //std::cout << n << std::endl;

  solver->add(this);
  last.initialise(s,0);
}

Mistral::Lexicographic::~Lexicographic() {}

Mistral::Variable Mistral::Lexicographic::select() {

  // for(int i=0; i<order.size; ++i) {
  //   std::cout << order[i] << " in " << order[i].get_domain() << " ";
  // }
  // std::cout << std::endl << (int)last << std::endl;
  // for(int i=last; i<order.size; ++i) {
  //   std::cout << order[i] << " in " << order[i].get_domain() << " ";
  // }
  // std::cout << std::endl;

  while(last<(int)(order.size) && order[last].is_ground()) { 
    ++last;
  }
  return order[last];
}

void Mistral::Lexicographic::notify_change(const int idx) {
  int ido = index[idx];
  if(ido>=0) {
    order[ido] = solver->variables[idx];
  }
}

// std::ostream& operator<<(std::ostream& os, Mistral::BranchingHeuristic& x) {
//   return x.display(os);
// }

// std::ostream& operator<<(std::ostream& os, Mistral::BranchingHeuristic* x) {
//   return x->display(os);
// }


std::ostream& operator<<(std::ostream& os, Mistral::DecisionListener& x) {
  return x.display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::RestartListener& x) {
  return x.display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::SuccessListener& x) {
  return x.display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::BacktrackListener& x) {
  return x.display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::ConstraintListener& x) {
  return x.display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::VariableListener& x) {
  return x.display(os);
}

// std::ostream& operator<<(std::ostream& os, Mistral::DecisionListener* x) {
//   return x->display(os);
// }

// std::ostream& operator<<(std::ostream& os, Mistral::RestartListener* x) {
//   return x->display(os);
// }

// std::ostream& operator<<(std::ostream& os, Mistral::SuccessListener* x) {
//   return x->display(os);
// }

// std::ostream& operator<<(std::ostream& os, Mistral::BacktrackListener* x) {
//   return x->display(os);
// }

// std::ostream& operator<<(std::ostream& os, Mistral::ConstraintListener* x) {
//   return x->display(os);
// }

// std::ostream& operator<<(std::ostream& os, Mistral::VariableListener* x) {
//   return x->display(os);
// }

double *weight_sorting_array;

int Mistral::decreasing_weight(const void *x, const void *y) {
    int _x = *(int*)x;
    int _y = *(int*)y;

    int r_value = 0;
    if(weight_sorting_array[_x] > weight_sorting_array[_y]) {
      r_value = -1;
    } else if(weight_sorting_array[_x] < weight_sorting_array[_y]) {
      r_value = 1;
    }

    return r_value;
  }

int log10(const double x){
  long int y = (long int)x;
  int log = 0;
  while(y>0) {
    ++log;
    y/=10;
  }
  return log;
}

std::ostream& Mistral::ImpactManager::display(std::ostream& os, const bool all) const {
  


  return os;
}

std::ostream& Mistral::FailureCountManager::display(std::ostream& os, const bool all) const {
      
      int *all_variables = new int[variable_weight.size];
      int *all_constraints = new int[constraint_weight.size];


      int w, 
	xwidth; //, // = log10(solver->variables[variable_weight.size-1].id()),
	//cwidth; // = log10(solver->constraints[constraint_weight.size-1].id());
    
      for(unsigned int i=0; i<variable_weight.size; ++i) {
	all_variables[i] = i;

	// w = log10(variable_weight[i]);
	// if(w>xwidth) xwidth = w;

      }

      for(unsigned int i=0; i<constraint_weight.size; ++i) {
	all_constraints[i] = i;

	// w = log10(constraint_weight[i]);
	// if(w>cwidth) cwidth = w;

      }


      weight_sorting_array = variable_weight.stack_;
      qsort(all_variables, variable_weight.size, sizeof(int), decreasing_weight);

      weight_sorting_array = constraint_weight.stack_;
      qsort(all_constraints, constraint_weight.size, sizeof(int), decreasing_weight);

      os << " c variable weight: \n c id: ";
      for(unsigned int i=0; i<variable_weight.size; ++i) {
	if(all || solver->sequence.contain(all_variables[i])) {
	  xwidth = log10(variable_weight[all_variables[i]]);
	  w = log10(all_variables[i]);
	  if(w>xwidth) xwidth = w;

	  os << std::setw(xwidth) << all_variables[i] << " ";
	}
      }
      os << "\n c va: ";
      for(unsigned int i=0; i<variable_weight.size; ++i) {
	if(all || solver->sequence.contain(all_variables[i])) {
	  xwidth = log10(variable_weight[all_variables[i]]);
	  w = log10(all_variables[i]);
	  if(w>xwidth) xwidth = w;
	  
	  os << std::setw(xwidth) << variable_weight[all_variables[i]] << " ";
	}
      }
      os << "\n c constraint weight: \n c id: ";
      for(unsigned int i=0; i<constraint_weight.size; ++i) {
	xwidth = log10(constraint_weight[all_constraints[i]]);
	w = log10(all_constraints[i]);
	if(w>xwidth) xwidth = w;

	os << std::setw(xwidth) << all_constraints[i] << " ";
      }
      os << "\n c va: ";
      for(unsigned int i=0; i<constraint_weight.size; ++i) {
	xwidth = log10(constraint_weight[all_constraints[i]]);
	w = log10(all_constraints[i]);
	if(w>xwidth) xwidth = w;

	os << std::setw(xwidth) << constraint_weight[all_constraints[i]] << " ";
      }
      os << std::endl;

      delete [] all_constraints;
      delete [] all_variables;

      return os;
    }    



std::ostream& Mistral::ConflictCountManager::display(std::ostream& os, const bool all) const {
      
      int *all_variables = new int[variable_weight.size];

      int w, 
	xwidth; //, // = log10(solver->variables[variable_weight.size-1].id()),
	//cwidth; // = log10(solver->constraints[constraint_weight.size-1].id());
    
      for(unsigned int i=0; i<variable_weight.size; ++i) {
	all_variables[i] = i;

	// w = log10(variable_weight[i]);
	// if(w>xwidth) xwidth = w;

      }

      weight_sorting_array = variable_weight.stack_;
      qsort(all_variables, variable_weight.size, sizeof(int), decreasing_weight);


      os << " c variable weight: \n c id: ";
      for(unsigned int i=0; i<variable_weight.size; ++i) {
	if(all || solver->sequence.contain(all_variables[i])) {
	  xwidth = log10(variable_weight[all_variables[i]]);
	  w = log10(all_variables[i]);
	  if(w>xwidth) xwidth = w;

	  os << std::setw(xwidth) << all_variables[i] << " ";
	}
      }
      os << "\n c va: ";
      for(unsigned int i=0; i<variable_weight.size; ++i) {
	if(all || solver->sequence.contain(all_variables[i])) {
	  xwidth = log10(variable_weight[all_variables[i]]);
	  w = log10(all_variables[i]);
	  if(w>xwidth) xwidth = w;
	  
	  os << std::setw(xwidth) << variable_weight[all_variables[i]] << " ";
	}
      }
      os << std::endl;

      delete [] all_variables;

      return os;
    }    

std::ostream& Mistral::PruningCountManager::display(std::ostream& os, const bool all) const {
      
      int *all_variables = new int[variable_weight.size];


      int w, 
	xwidth; //, // = log10(solver->variables[variable_weight.size-1].id()),
	//cwidth; // = log10(solver->constraints[constraint_weight.size-1].id());
    
      for(unsigned int i=0; i<variable_weight.size; ++i) {
	all_variables[i] = i;

	// w = log10(variable_weight[i]);
	// if(w>xwidth) xwidth = w;

      }


      weight_sorting_array = variable_weight.stack_;
      qsort(all_variables, variable_weight.size, sizeof(int), decreasing_weight);

      os << " c variable weight: \n c id: ";
      for(unsigned int i=0; i<variable_weight.size; ++i) {
	if(all || solver->sequence.contain(all_variables[i])) {
	  xwidth = log10(variable_weight[all_variables[i]]);
	  w = log10(all_variables[i]);
	  if(w>xwidth) xwidth = w;
	  
	  os << std::setw(xwidth) << all_variables[i] << " ";
	}
      }
      os << "\n c va: ";
      for(unsigned int i=0; i<variable_weight.size; ++i) {
	if(all || solver->sequence.contain(all_variables[i])) {
	  xwidth = log10(variable_weight[all_variables[i]]);
	  w = log10(all_variables[i]);
	  if(w>xwidth) xwidth = w;
	  
	  os << std::setw(xwidth) << variable_weight[all_variables[i]] << " ";
	}
      }
      os << std::endl;

    delete [] all_variables;

      return os;
    }    


Mistral::LearningActivityManager::LearningActivityManager(Solver *s) : solver(s) {
  weight_unit = solver->parameters.activity_increment;
  decay = solver->parameters.activity_decay;
  max_weight = std::numeric_limits<int>::max();
  
  var_activity.initialise(solver->variables.size, solver->variables.size, 0);
  lit_activity.initialise(2*solver->variables.size, 2*solver->variables.size, 0);
  
  int i = solver->constraints.size;
  Constraint *cons = solver->constraints.stack_;
  while(i--) {
    cons[i].initialise_activity(lit_activity.stack_, var_activity.stack_, weight_unit);
  }

  max_activity = 0;
  i = var_activity.size;
  while(i--) {
    if(var_activity[i]>max_activity)
      max_activity = var_activity[i];
  }

#ifdef _DEBUG_ACTIVITY
  for(int a=0; a<solver->variables.size; ++a) {
    std::cout << "init x" << a << " (" << lit_activity[NEG(a)] << "/" << lit_activity[POS(a)] << ")/" << var_activity[a]
	       << std::endl; 
  }
#endif
  
  //solver->lit_activity = lit_activity.stack_;
  //solver->var_activity = var_activity.stack_;
  
  solver->add((BacktrackListener*)this);
}

Mistral::LearningActivityManager::~LearningActivityManager() {
  solver->remove((BacktrackListener*)this);
}


void Mistral::LearningActivityManager::notify_backtrack() {

      // //std::cout << "d " << lit_activity.stack_ << " " << lit_activity[0] << " " << lit_activity[1] << std::endl;


#ifdef _DEBUG_ACTIVITY
  std::cout << "NOTIFY BACKTRACK!" << std::endl;
  std::cout << std::endl;
#endif

  int i;
  Literal q;
  Atom a;

  weight_unit /= decay;
  if(max_weight - weight_unit <= max_activity) {
    // risk of double overflow

#ifdef _DEBUG_ACTIVITY
    std::cout << "\n RISK OF OVERFLOW (" << max_activity << " + " << weight_unit << " >= " << max_weight << ")\n";
#endif

    i=lit_activity.size;
    while(i--) lit_activity[i] /= max_activity;

    i=var_activity.size;
    while(i--) {
#ifdef _DEBUG_ACTIVITY
      std::cout << "x" << i << " (" << lit_activity[NEG(i)] << "/" << lit_activity[POS(i)] << ")/" << var_activity[i] << " -> " 
		<< (lit_activity[POS(i)]+lit_activity[NEG(i)]) << std::endl;
#endif
      var_activity[i] = lit_activity[POS(i)]+lit_activity[NEG(i)];
    }
   


    weight_unit = 1.0/decay;
    max_activity = 1.0;
  }
  i = solver->visited_literals.size;
  while(i--) {

#ifdef _DEBUG_ACTIVITY
    std::cout << "\n UPDATE ACTIVITY BECAUSE OF LITERAL " << q << ":\n";
#endif

    q = solver->visited_literals[i];
    Atom a = UNSIGNED(q);

#ifdef _DEBUG_ACTIVITY
    std::cout << "x" << a << " (" << lit_activity[NEG(a)] << "/" << lit_activity[POS(a)] << ")/" << var_activity[a]
     	      << " -> " ;     
#endif

    lit_activity[q] += weight_unit;
    var_activity[a] += weight_unit;
    if(var_activity[a] > max_activity)
      max_activity = var_activity[a];

    
#ifdef _DEBUG_ACTIVITY
    std::cout << " (" << lit_activity[NEG(a)] << "/" << lit_activity[POS(a)] << ")/" << var_activity[a]
     	      << std::endl ; 
#endif

  }
//   //std::cout << std::endl ; 


// #else

//       if(decay > 0 && decay < 1) {
//       	int i=var_activity.size;
//       	while(i--) {

// // #ifdef _DEBUG_ACTIVITY
// // 	  int a = i;
// // 	  std::cout << "decay x" << a << " (" << lit_activity[NEG(a)] << "/" << lit_activity[POS(a)] << ")/" << var_activity[a]
// // 		 << " -> ";
// // #endif

// 	  var_activity[i] *= decay;

// // #ifdef _DEBUG_ACTIVITY
// // 	  std::cout << " (" << lit_activity[NEG(a)] << "/" << lit_activity[POS(a)] << ")/" << var_activity[a]
// // 		    << std::endl;
// // #endif

// 	}
//       	// i=lit_activity.size;
//       	// while(i--) lit_activity[i] *= decay;
//       }

//   //solver->parameters.activity_increment *= (2.0-decay);

// #endif

}


std::ostream& Mistral::LearningActivityManager::display(std::ostream& os, const bool all) const {

  
  int *all_variables = new int[var_activity.size];


      int w, 
	xwidth; //, // = log10(solver->variables[var_activity.size-1].id()),
	//cwidth; // = log10(solver->constraints[constraint_weight.size-1].id());
    
      for(unsigned int i=0; i<var_activity.size; ++i) {
	all_variables[i] = i;

	// w = log10(var_activity[i]);
	// if(w>xwidth) xwidth = w;

      }


      weight_sorting_array = var_activity.stack_;
      qsort(all_variables, var_activity.size, sizeof(int), decreasing_weight);

      os << " c variable weight: \n c id: ";
      for(unsigned int i=0; i<var_activity.size; ++i) {

	if(!(i%1000)) {

	if(all || solver->sequence.contain(all_variables[i])) {
	  xwidth = log10(var_activity[all_variables[i]])+7;
	  w = log10(all_variables[i]);
	  if(w>xwidth) xwidth = w;
	  
	  os << std::setw(xwidth) << all_variables[i] << " ";
	}

	}

      }
      os << "\n c va: ";
      for(unsigned int i=0; i<var_activity.size; ++i) {

	if(!(i%1000)) {

	if(all || solver->sequence.contain(all_variables[i])) {
	  xwidth = log10(var_activity[all_variables[i]])+7;
	  w = log10(all_variables[i]);
	  if(w>xwidth) xwidth = w;
	  
	  long long int intval = (long long int)(100000 * var_activity[all_variables[i]]);
	  double outputval = ((double)intval)/100000;

	  os << std::setw(xwidth) << outputval << " ";
	}

	}

      }
     os << "\n c  0: ";
      for(unsigned int i=0; i<var_activity.size; ++i) {

	if(!(i%1000)) {

	if(all || solver->sequence.contain(all_variables[i])) {
	  xwidth = log10(var_activity[2*all_variables[i]])+7;
	  w = log10(all_variables[i]);
	  if(w>xwidth) xwidth = w;
	  
	  long long int intval = (long long int)(100000 * lit_activity[2*all_variables[i]]);
	  double outputval = ((double)intval)/100000;

	  os << std::setw(xwidth) << outputval << " ";
	}

	}

      }
     os << "\n c  1: ";
      for(unsigned int i=0; i<var_activity.size; ++i) {

	if(!(i%1000)) {

	if(all || solver->sequence.contain(all_variables[i])) {
	  xwidth = log10(lit_activity[2*all_variables[i]+1])+7;
	  w = log10(all_variables[i]);
	  if(w>xwidth) xwidth = w;
	  
	  long long int intval = (long long int)(100000 * lit_activity[2*all_variables[i]+1]);
	  double outputval = ((double)intval)/100000;

	  os << std::setw(xwidth) << outputval << " ";
	}

	}

      }
      os << std::endl;

    delete [] all_variables;


  return os;
}

std::ostream& operator<<(std::ostream& os, Mistral::MinDomain& x) {
  return x.display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::MinMin& x) {
  return x.display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::MaxMax& x) {
  return x.display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::MaxRegret& x) {
  return x.display(os);
}

// std::ostream& operator<<(std::ostream& os, Mistral::MinDomainMaxDegree& x) {
//   return x.display(os);
// }

std::ostream& operator<<(std::ostream& os, Mistral::MinDomainOverDegree& x) {
  return x.display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::MinDomainOverWeight& x) {
  return x.display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::MinDomainTimesWeight& x) {
  return x.display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::MinNeighborDomainOverWeight& x) {
  return x.display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::MinNeighborDomainOverNeighborWeight& x) {
  return x.display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::MaxWeight& x) {
  return x.display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::MinWeight& x) {
  return x.display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::AnyValue& x) {
  return x.display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::MinValue& x) {
  return x.display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::MaxValue& x) {
  return x.display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::MiddleValue& x) {
  return x.display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::MedianValue& x) {
  return x.display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::RandomValue& x) {
  return x.display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::HalfSplit& x) {
  return x.display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::ReverseSplit& x) {
  return x.display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::RandomSplit& x) {
  return x.display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::RandomMinMax& x) {
  return x.display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::MinWeightValue& x) {
  return x.display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::MinWeightBound& x) {
  return x.display(os);
}

// std::ostream& operator<<(std::ostream& os, Mistral::Guided& x) {
//   return x.display(os);
// }

std::ostream& operator<<(std::ostream& os, Mistral::BoolMinWeightValue& x) {
  return x.display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::BoolMaxWeightValue& x) {
  return x.display(os);
}



std::ostream& operator<<(std::ostream& os, Mistral::MinDomain* x) {
  return x->display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::MinMin* x) {
  return x->display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::MaxMax* x) {
  return x->display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::MaxRegret* x) {
  return x->display(os);
}

// std::ostream& operator<<(std::ostream& os, Mistral::MinDomainMaxDegree* x) {
//   return x->display(os);
// }

std::ostream& operator<<(std::ostream& os, Mistral::MinDomainOverDegree* x) {
  return x->display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::MinDomainOverWeight* x) {
  return x->display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::MinDomainTimesWeight* x) {
  return x->display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::MinNeighborDomainOverWeight* x) {
  return x->display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::MinNeighborDomainOverNeighborWeight* x) {
  return x->display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::MaxWeight* x) {
  return x->display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::MinWeight* x) {
  return x->display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::AnyValue* x) {
  return x->display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::MinValue* x) {
  return x->display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::MaxValue* x) {
  return x->display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::MiddleValue* x) {
  return x->display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::MedianValue* x) {
  return x->display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::RandomValue* x) {
  return x->display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::HalfSplit* x) {
  return x->display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::ReverseSplit* x) {
  return x->display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::RandomSplit* x) {
  return x->display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::RandomMinMax* x) {
  return x->display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::MinWeightValue* x) {
  return x->display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::MinWeightBound* x) {
  return x->display(os);
}

// std::ostream& operator<<(std::ostream& os, Mistral::Guided* x) {
//   return x->display(os);
// }

std::ostream& operator<<(std::ostream& os, Mistral::BoolMinWeightValue* x) {
  return x->display(os);
}

std::ostream& operator<<(std::ostream& os, Mistral::BoolMaxWeightValue* x) {
  return x->display(os);
}
