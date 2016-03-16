
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


#ifndef __SEARCH_HPP
#define __SEARCH_HPP

#include <algorithm>

#include <mistral_global.hpp>
#include <mistral_structure.hpp>
#include <mistral_variable.hpp>
#include <mistral_solver.hpp>
//#include <mistral_sat.hpp>

//#define _DEBUG_IMPACT true

//#define _DEBUG_VARORD

namespace Mistral {


  // class VariableWeightMapping {

  // public:
    
  //   // returns the weight
  //   double operator[](const int x) = 0;

  // };





  
  /**********************************************
   * Listener
   **********************************************/

  /*! \class RestartListener
    \brief RestartListener Class

    * Called whenever the solver restarts *
    
    This is used to implement procedures triggered by restarts
  */
  class RestartListener {
  public:
    int rid;
    virtual void notify_restart(const double prog) = 0;
    
    std::ostream& display(std::ostream& os) { os << "restart-L"; return os; }    
  };

  /*! \class SolutionListener
    \brief SolutionListener Class

    * Called whenever the solver solutions *
    
    This is used to implement procedures triggered by solutions
  */
  class SolutionListener {
  public:
    int mid;
    virtual void notify_solution() = 0;
    
    std::ostream& display(std::ostream& os) { os << "solution-L"; return os; }    
  };

  /*! \class DecisionListener
    \brief DecisionListener Class

    * Called whenever the solver restarts *
    
    This is used to implement procedures triggered by decisions
  */
  class DecisionListener {
  public:
    int did;
    virtual void notify_decision() = 0;

    std::ostream& display(std::ostream& os) { os << "decision-L"; return os; }    
  };

  /*! \class SuccessListener
    \brief SuccessListener Class

    * Called whenever the solver restarts *
    
    This is used to implement procedures triggered on successful propagation steps
  */
  class SuccessListener {
  public:
    int sid;
    virtual void notify_success() = 0;

    std::ostream& display(std::ostream& os) { os << "success-L"; return os; }    
  };

  /*! \class BacktrackListener
    \brief BacktrackListener Class

    * Called whenever the solver restarts *
    
    This is used to implement procedures triggered on failed propagation steps
  */
  class BacktrackListener {
  public:
    int fid;
    virtual void notify_backtrack() = 0;

    std::ostream& display(std::ostream& os) { os << "backtrack-L"; return os; }    
  };

  /*! \class ConstraintListener
    \brief ConstraintListener Class

    * Called whenever the solver restarts *
    
    This is used to implement procedures triggered wheneever a new
    constraint is added to the solver, or an old constraint is posted/relaxed
  */
  class ConstraintListener {
  public:
    int cid;
    virtual void notify_post   (Constraint c) = 0;
    virtual void notify_relax  (Constraint c) = 0;
    virtual void notify_add_con(Constraint c) = 0;

    std::ostream& display(std::ostream& os) { os << "constraint-L"; return os; }    
  };

  /*! \class VariableListener
    \brief VariableListener Class

    * Called whenever the solver restarts *
    
    This is used to implement procedures triggered wheneever a new
    variable is added to the solver, or when its domain representation changes
  */
  class VariableListener {
  public:
    int vid;
    virtual void notify_add_var() = 0;
    virtual void notify_change(const int idx) = 0;

    std::ostream& display(std::ostream& os) { os << "variable-L"; return os; }    
  };

  std::ostream& operator<<(std::ostream& os, DecisionListener& x);
  std::ostream& operator<<(std::ostream& os, RestartListener& x);
  std::ostream& operator<<(std::ostream& os, SuccessListener& x);
  std::ostream& operator<<(std::ostream& os, BacktrackListener& x);
  std::ostream& operator<<(std::ostream& os, ConstraintListener& x);
  std::ostream& operator<<(std::ostream& os, VariableListener& x);

  // std::ostream& operator<<(std::ostream& os, DecisionListener* x);
  // std::ostream& operator<<(std::ostream& os, RestartListener* x);
  // std::ostream& operator<<(std::ostream& os, SuccessListener* x);
  // std::ostream& operator<<(std::ostream& os, BacktrackListener* x);
  // std::ostream& operator<<(std::ostream& os, ConstraintListener* x);
  // std::ostream& operator<<(std::ostream& os, VariableListener* x);


  // /*! \class LiteralActivityManager
  //   \brief LiteralActivityManager Class

  //   * Listener interface for literal activity *

  //   This structure is shared with the clause base constraint.
  //   The constraint updates the activity when doing conflict analysis
  //   while this listener implements the decay.
  // */
  // class LiteralActivityManager : public DecisionListener {

  // public:

  //   /*\ TODO: change double* to vectors* and make it a variable listener \*/
  //   Solver *solver;
  //   double *lit_activity; 
  //   double *var_activity;
  //   int n_vars;
 
  //   double decay;

  //   //LiteralActivityManager(Solver *s, void *a=NULL) ;
  //   LiteralActivityManager(Solver *s) ;
  //   virtual ~LiteralActivityManager() ;

  //   virtual void notify_decision() ;    
  //   double *get_variable_weight() ;
  //   double **get_value_weight() ;
  //   double *get_bound_weight() ;

  //   virtual std::ostream& display(std::ostream& os, const bool all) const {return os;}
    
  // };


  //double *weight_sorting_array;
  int decreasing_weight(const void *x, const void *y);


  /*! \class NoManager
    \brief NoManager Class
    
    * Manager doing nothing, used for genericity *
    */
  class NoManager {
    
  public:

    NoManager(Solver *s) {}
    virtual ~NoManager() {}

    double *get_variable_weight() { return NULL; }  
    double **get_value_weight() { return NULL; }
    double *get_bound_weight() { return NULL; }

    std::ostream& display(std::ostream& os, const bool b=false) { return os; }    
    // void initialise(Solver* s);
    // void initialise_structures(VarComparator &v);
  };



  /*! \class HeuristicPoolManager
    \brief HeuristicPoolManager Class

    * Listener interface for handling a pool of search strategies *
  */
  //class BranchingHeuristic;
  class HeuristicPoolManager : public RestartListener {

  public:

    Solver *solver;

    Vector<BranchingHeuristic*> pool;
    int threshold;
    int counter;
    int heu_index;


    HeuristicPoolManager(Solver *s);//  : solver(s) {// }

    //   //std::cout << " c add restart listener" << std::endl;

    //   heu_index = 1;
    //   solver->add((RestartListener*)this);
    // }

    virtual ~HeuristicPoolManager();//  {// }
    //   // for(unsigned int i=0; i<pool.size; ++i) {
    //   // 	if(solver->heuristic != pool[i]) {
    //   // 	  delete [] pool[i];
    //   // 	}
    //   // }
    //   solver->remove((RestartListener*)this);
    // }


    void add(BranchingHeuristic *h) {
      pool.add(h);
    }

    void set_threshold(const int t) {
      threshold = t;
      counter = t;
    }

    virtual void notify_restart(const double prog);//  {
    //   if(!(solver->statistics.num_restarts % threshold) && ++counter < pool.size) {

    // 	std::cout << " c SWITCH HEURISTIC!!\n";

    // 	//solver->heuristic = pool[counter];
    // 	solver->heuristic = new GenericHeuristic< VSIDS<2>, MaxValue >(solver);
    //   }
    // }

    virtual std::ostream& display(std::ostream& os, const bool all) const {
      for(unsigned int i=0; i<pool.size; ++i) {
	os << pool[i] << std::endl;
      }
      return os;
    }

  };




  /*! \class FailureCountanager
    \brief FailureCountManager Class

    * Listener interface for weighted degree *
    * Counts the number of failures for each constraint *
  */
  //template< float DECAY > 
  class FailureCountManager : public BacktrackListener, public ConstraintListener {

  public:

    Solver *solver;
    double weight_unit;


    /*\ TODO: make it a variable listener \*/
    Vector<double> constraint_weight;
    Vector<double> variable_weight;

    //FailureCountManager(Solver *s, void *a=NULL) : solver(s) {// }
    FailureCountManager(Solver *s) : solver(s) {// }

      weight_unit = solver->parameters.activity_increment;
      
      variable_weight.initialise(solver->variables.size, solver->variables.size);

      for(unsigned int i=0; i<solver->variables.size; ++i) {
	variable_weight[i] = 0;
      // 	variable_weight.add(weight_unit * solver->variables[i].get_degree());
      }
      Variable *scope;
      int arity, j;
      for(unsigned int i=0; i<solver->constraints.size; ++i) {
	constraint_weight.add(weight_unit);
	arity = solver->constraints[i].arity();
	scope = solver->constraints[i].get_scope();
	for(j=0; j<arity; ++j) if(!scope[j].is_ground()) {
	    variable_weight[scope[j].id()] += weight_unit/(double)arity;
	  }
      }

      
      // std::cout << variable_weight.stack_ << ":";
      // for(unsigned int i=0; i<solver->variables.size; ++i) {
      // 	std::cout << " " << variable_weight[i] ;
      // }
      // std::cout << std::endl;


      solver->add((BacktrackListener*)this);
      solver->add((ConstraintListener*)this);
    }

    virtual ~FailureCountManager() {// }

      solver->remove((ConstraintListener*)this);
      solver->remove((BacktrackListener*)this);
    }

    double *get_variable_weight() { return variable_weight.stack_; }   
    double **get_value_weight() { return NULL; }
    double *get_bound_weight() { return NULL; }

    virtual void check_consistency() {
      
      double xweight;
      
      solver->display(std::cout, 1);

      for(unsigned int i=0; i<variable_weight.size; ++i) {
	
	if(!(solver->domain_types[i] & REMOVED_VAR) && solver->sequence.contain(i)) {

	  xweight = 0;
	  for(Event trig = 0; trig<3; ++trig) 
	    for(int cons = solver->constraint_graph[i].on[trig].size; --cons>=0;) {
	      xweight += constraint_weight[solver->constraint_graph[i].on[trig][cons].id()];
	    }

	  if(xweight != variable_weight[i]) {

	    std::cout << "WARNING! inconsistency: on " << solver->variables[i] << ": " 
		      << variable_weight[i] << " should be " << xweight << std::endl;

	  } else {
	    
	    std::cout << "OK!" << std::endl;

	  }
	}
      }

    }
  

    virtual void notify_backtrack() {
      int i;
      Constraint con = solver->culprit;


      //std::cout << "failure on " << con << std::endl;

      if(!con.empty()) {
	Variable *scope = con.get_scope();
	int idx;
	i = con.arity();
	//++constraint_weight[con.id()];
	constraint_weight[con.id()] += weight_unit;
	while(i--) {
	  idx = scope[i].id();
	  if(idx>=0) {
	    //std::cout << " ++x" << idx; 
	    variable_weight[idx] += weight_unit;
	  }
	}
      } 
      // std::cout << std::endl;

      // display(std::cout, false);

    }

    virtual void notify_post(Constraint con) {
      int i = con.num_active(), idx;
      Variable *scope = con.get_scope();
      while(i--) {
	idx = scope[con.get_active(i)].id();
	if(idx>=0) variable_weight[idx] += constraint_weight[con.id()];
      }
    }

    virtual void notify_relax(Constraint con) {
      int i = con.num_active(), idx;
      Variable *scope = con.get_scope();
      while(i--) {
	idx = scope[con.get_active(i)].id();
	if(idx>=0) variable_weight[idx] -= constraint_weight[con.id()];
      }
    }

    virtual void notify_add_con(Constraint con) {
			
      while(constraint_weight.size < solver->constraints.size) {
	constraint_weight.add(weight_unit);
      }

      while(variable_weight.size < solver->variables.size) {
	variable_weight.add(weight_unit*solver->variables[variable_weight.size].get_degree());
      }
    }

    virtual std::ostream& display(std::ostream& os, const bool all) const ;
// {
      
//       int *all_variables = new int[variable_weight.size];
//       int *all_constraints = new int[constraint_weight.size];

//       for(unsigned int i=0; i<variable_weight.size; ++i) {
// 	all_variables[i] = i;
//       }

//       for(unsigned int i=0; i<constraint_weight.size; ++i) {
// 	all_constraints[i] = i;
//       }

//       weight_sorting_array = variable_weight.stack_;
//       qsort(all_variables, variable_weight.size, sizeof(int), decreasing_weight);

//       weight_sorting_array = constraint_weight.stack_;
//       qsort(all_constraints, constraint_weight.size, sizeof(int), decreasing_weight);

//       for(unsigned int i=0; i<variable_weight.size; ++i) {
// 	os << std::setw(5) << solver->variables[all_variables[i]] << " ";
//       }
//       os << std::endl;
//       for(unsigned int i=0; i<variable_weight.size; ++i) {
// 	os << std::setw(5) << variable_weight[i] << " ";
//       }
//       os << std::endl;

//      for(unsigned int i=0; i<constraint_weight.size; ++i) {
// 	os << std::setw(5) << solver->constraints[all_constraints[i]] << " ";
//       }
//       os << std::endl;
//       for(unsigned int i=0; i<constraint_weight.size; ++i) {
// 	os << std::setw(5) << constraint_weight[i] << " ";
//       }
//       os << std::endl;

//     }    

  };



  /*! \class ConflictCountManager
    \brief ConflictCountManager Class

    * Listener interface for a kind of weighted degree *
    * Failed constraints weight variables with an aribtrary policy *
  */
  //template< float DECAY > 
  class ConflictCountManager : public BacktrackListener {

  public:

    Solver *solver;
    double weight_unit;


    /*\ TODO: make it a variable listener \*/
    //Vector<double> constraint_weight;
    Vector<double> variable_weight;

    //ConflictCountManager(Solver *s, void *a=NULL) : solver(s) {// }
    ConflictCountManager(Solver *s) : solver(s) {// }

      //std::cout << "NEW" << std::endl;

      weight_unit = solver->parameters.activity_increment;
      
      variable_weight.initialise(solver->variables.size, solver->variables.size);

      for(unsigned int i=0; i<solver->variables.size; ++i) {
	variable_weight[i] = 0;
      // 	variable_weight.add(weight_unit * solver->variables[i].get_degree());
      }
      Variable *scope;
      int arity, j;
      for(unsigned int i=0; i<solver->constraints.size; ++i) {
	arity = solver->constraints[i].arity();
	scope = solver->constraints[i].get_scope();
	for(j=0; j<arity; ++j) if(!scope[j].is_ground()) {
	    variable_weight[scope[j].id()] += weight_unit/(double)arity;
	  }
      }

      solver->add((BacktrackListener*)this);
      solver->add((ConstraintListener*)this);
    }

    virtual ~ConflictCountManager() {
      solver->remove((ConstraintListener*)this);
      solver->remove((BacktrackListener*)this);
    }

    double *get_variable_weight() { return variable_weight.stack_; }   
    double **get_value_weight() { return NULL; }
    double *get_bound_weight() { return NULL; }
  

    virtual void notify_backtrack() {
      int i;
      Constraint con = solver->culprit;

      con.weight_conflict(weight_unit, variable_weight);
    }

    virtual std::ostream& display(std::ostream& os, const bool all) const ;
// {
      
//       int *all_variables = new int[variable_weight.size];
//       int *all_constraints = new int[constraint_weight.size];

//       for(unsigned int i=0; i<variable_weight.size; ++i) {
// 	all_variables[i] = i;
//       }

//       for(unsigned int i=0; i<constraint_weight.size; ++i) {
// 	all_constraints[i] = i;
//       }

//       weight_sorting_array = variable_weight.stack_;
//       qsort(all_variables, variable_weight.size, sizeof(int), decreasing_weight);

//       weight_sorting_array = constraint_weight.stack_;
//       qsort(all_constraints, constraint_weight.size, sizeof(int), decreasing_weight);

//       for(unsigned int i=0; i<variable_weight.size; ++i) {
// 	os << std::setw(5) << solver->variables[all_variables[i]] << " ";
//       }
//       os << std::endl;
//       for(unsigned int i=0; i<variable_weight.size; ++i) {
// 	os << std::setw(5) << variable_weight[i] << " ";
//       }
//       os << std::endl;

//      for(unsigned int i=0; i<constraint_weight.size; ++i) {
// 	os << std::setw(5) << solver->constraints[all_constraints[i]] << " ";
//       }
//       os << std::endl;
//       for(unsigned int i=0; i<constraint_weight.size; ++i) {
// 	os << std::setw(5) << constraint_weight[i] << " ";
//       }
//       os << std::endl;

//     }    

  };


  /*! \class PruningCountManager
    \brief PruningCountManager Class

    * Listener interface for ABS *
    * Counts the number of nodes in which at least one pruning event occurred for each variable *
  */
  class PruningCountManager : public SuccessListener// , public VariableListener
  {

  public:

    Solver *solver;
    double weight_unit;

    Vector<double> variable_weight;

    //PruningCountManager(Solver *s, void *a=NULL) : solver(s) {
    PruningCountManager(Solver *s) : solver(s) {

      weight_unit = solver->parameters.activity_increment;
      variable_weight.initialise(solver->variables.size, solver->variables.size);

      for(unsigned int i=0; i<solver->variables.size; ++i) {
	//variable_weight.add(weight_unit * (double)(solver->variables[i].get_degree()));
	variable_weight[i] = (weight_unit * (double)(solver->variables[i].get_degree()));
      }

      // std::cout << variable_weight.stack_ << ":";
      // for(unsigned int i=0; i<solver->variables.size; ++i) {
      // 	std::cout << " " << solver->variables[i] << ": " << variable_weight[i] ;
      // }
      // std::cout << std::endl;


      solver->add((SuccessListener*)this);
      //solver->add((VariableListener*)this);
    }

    virtual ~PruningCountManager() {
      solver->remove((SuccessListener*)this);
      //solver->remove((VariableListener*)this);
    }

    double *get_variable_weight() { return variable_weight.stack_; }     
    double **get_value_weight() { return NULL; }
    double *get_bound_weight() { return NULL; }

    virtual void notify_success() {
      int id;
      int i = solver->trail_.back(5), n=solver->saved_vars.size;

      //std::cout << "increment weight of ";
      while(++i<n) {	
	id = solver->saved_vars[i]; 
	variable_weight[id] += weight_unit;
      }
    }

   // virtual void notify_add_variable() {
   //   while(variable_weight.size < solver->variables.size) {
   //     variable_weight.add(weight_unit*solver->variables[variable_weight.size].get_degree());
   //   }
   // }

    virtual std::ostream& display(std::ostream& os, const bool all) const ;
  };



  /*! \class PruningActivityManager
    \brief PruningActivityManager Class
    
    * Listener interface for ABS *
    * Activitys the number of times each variable was visited when computing a nogood *
    */
  class LearningActivityManager : public BacktrackListener {

  public:

    Solver *solver;
    double weight_unit;
    double max_activity;
    double max_weight;

    Vector<double> var_activity;
    Vector<double> lit_activity;

    double decay;

    //LearningActivityManager(Solver *s, void *a=NULL) : solver(s) {
    LearningActivityManager(Solver *s);//  : solver(s) {
    //   weight_unit = solver->parameters.activity_increment;
    //   decay = solver->parameters.activity_decay;

    //   var_activity.initialise(solver->variables.size, solver->variables.size, 0);
    //   lit_activity.initialise(2*solver->variables.size, 2*solver->variables.size, 0);

    //   int i = solver->constraints.size;
    //   Constraint *cons = solver->constraints.stack_;
    //   while(i--) {
    //    	cons[i].initialise_activity(lit_activity.stack_, var_activity.stack_, weight_unit);
    //   }

    //   solver->lit_activity = lit_activity.stack_;
    //   solver->var_activity = var_activity.stack_;

    //   solver->add((DecisionListener*)this);
    // }

    virtual ~LearningActivityManager();//  {
    //   solver->remove((DecisionListener*)this);
    //   //solver->remove((VariableListener*)this);
    // }

    double *get_variable_weight() { return var_activity.stack_; }     
    double *get_bound_weight() { return lit_activity.stack_; }
    double **get_value_weight() { return NULL; }

    virtual void notify_backtrack() ;//{

    //   //std::cout << "d " << lit_activity.stack_ << " " << lit_activity[0] << " " << lit_activity[1] << std::endl;

    //   if(decay > 0 && decay < 1) {
    // 	int i=var_activity.size;
    // 	while(i--) var_activity[i] *= decay;
    // 	i=lit_activity.size;
    // 	while(i--) lit_activity[i] *= decay;
    //   }

    // }

    virtual std::ostream& display(std::ostream& os, const bool all) const ;
  };




//   /*! \class ImpactManager
//     \brief ImpactManager Class

//     * Listener interface for Impact *
//     * NB: this is a simplified version of Impact, where only the impact of left vs right branches are distinguished
//   */
//   class ImpactManager : public BacktrackListener, public SuccessListener, public DecisionListener, public VariableListener {

//   public:

//     Solver *solver;
//     double weight_unit;
//     int left;

//     /*\ TODO: make it a variable listener \*/
//     Vector<double> variable_weight;
//     Vector<double> left_weight;
//     Vector<double> right_weight;
//     Vector<double> avg_right_branches;
//     Vector<int> num_right_branches;
//     Vector<int> num_probes;
//     Vector<int> num_left_probes;
//     Vector<int> num_right_probes;

//     //ImpactManager(Solver *s, void *a=NULL) : solver(s) {// }
//     ImpactManager(Solver *s) : solver(s) {// }

//       left = -1;
//       weight_unit = solver->parameters.activity_increment;
      
//       variable_weight.initialise(solver->variables.size, solver->variables.size);
//       left_weight.initialise(solver->variables.size, solver->variables.size);
//       right_weight.initialise(solver->variables.size, solver->variables.size);
//       num_probes.initialise(solver->variables.size, solver->variables.size);
//       num_left_probes.initialise(solver->variables.size, solver->variables.size);
//       num_right_probes.initialise(solver->variables.size, solver->variables.size);
//       num_right_branches.initialise(solver->variables.size, solver->variables.size);
//       avg_right_branches.initialise(solver->variables.size, solver->variables.size);


//       std::cout << solver->variables << std::endl;

//       for(unsigned int i=0; i<solver->variables.size; ++i) {
// 	variable_weight[i] = 0;
// 	left_weight[i] = 0;
// 	right_weight[i] = 0;
// 	num_probes[i] = 0;
// 	num_left_probes[i] = 0;
// 	num_right_probes[i] = 0;
// 	num_right_branches[i] = 0;
// 	avg_right_branches[i] = solver->variables[i].get_size();
//       }

//       solver->add((BacktrackListener*)this);
//       solver->add((SuccessListener*)this);
//       solver->add((DecisionListener*)this);
//       solver->add((VariableListener*)this);
//     }

//     virtual ~ImpactManager() {// }

//       solver->remove((VariableListener*)this);
//       solver->remove((SuccessListener*)this);
//       solver->remove((DecisionListener*)this);
//       solver->remove((BacktrackListener*)this);
//     }

//     double *get_variable_weight() { return variable_weight.stack_; }   
//     double **get_value_weight() { return NULL; }
//     double *get_bound_weight() { return NULL; }

//     virtual void check_consistency() {
//     }

//     virtual void notify_add_var() {
// #ifdef _DEBUG_IMPACT
//       std::cout << "NEW VAR!!\n";
// #endif
//     }

//     virtual void notify_change(const int id) {
// #ifdef _DEBUG_IMPACT
//       std::cout << solver->variables[id] << "'s domain implementation has changed!!\n";
// #endif
//     }

//     virtual void notify_decision() {
//       left = 1;
//     }
  
//     virtual void notify_success() {
//       // propagation went without wipe-out
//       // - check if it was after a left or a right branch
//       // - find out what was the decision/refutation
//       int dec=-1, nxt=-1;
//       int id;
//       int i, n;
//       Variable x;
//       double residual_space;
//       int size;
//       if(left>=0) {
// 	i = solver->trail_.back(5), n=solver->saved_vars.size;
// 	if(left) {
// 	  // left branch
// 	  dec = solver->decisions.back().var.id();

// 	  // if(num_right_branches[dec]>=0) {
// 	  //   avg_right_branches[dec] = avg_right_branches[dec]
// 	  // }
// 	  //num_right_branches[dec] = 0;

// #ifdef _DEBUG_IMPACT
// 	  std::cout << " notified of success after a left branch on " << solver->variables[dec] << "\n";
// #endif

// 	} else {
// 	  // right branch
// 	  if(!solver->decisions.empty()) dec = solver->decisions.back().var.id();
// 	  nxt = solver->decisions.back(0).var.id();
// 	  ++num_right_branches[nxt];

// #ifdef _DEBUG_IMPACT
// 	  std::cout << " notified of success after the " << num_right_branches[nxt] 
// 		    << "th right branch on " << solver->variables[nxt] << "\n";
// #endif

// 	}

// 	residual_space = 1.0;
// 	//std::cout << "increment weight of ";
// 	while(i<n) {	
// 	  id = solver->saved_vars[i];
// 	  x = solver->variables[id];

// #ifdef _DEBUG_IMPACT
// 	  std::cout << " -> " << x << x.get_domain() << " lost " << x.get_reduction() << " values (";
// #endif

// 	  size = x.get_size();
// 	  residual_space *= (((double)size))/((double)(size+x.get_reduction()));

// #ifdef _DEBUG_IMPACT
// 	  std::cout << residual_space << ")\n";
// #endif

// 	  ++i;
// 	} 

// 	if(left) {

// #ifdef _DEBUG_IMPACT
// 	  std::cout << " ==> left-weight[" << solver->variables[dec] << "] was " << left_weight[dec] ;
// #endif

// 	  left_weight[dec] = (((double)(num_left_probes[dec]) * left_weight[dec]) + residual_space)/(double)(++num_left_probes[dec]);

// #ifdef _DEBUG_IMPACT
// 	  std::cout << " now " << left_weight[dec] << std::endl;
// #endif

// 	  variable_weight[dec] = avg_right_branches[dec] * (left_weight[dec] + right_weight[dec]);

// #ifdef _DEBUG_IMPACT
// 	  std::cout << " ===> total weight of " << solver->variables[dec] << " = " << avg_right_branches[dec] << " * (" << left_weight[dec] << " + " << right_weight[dec] << ") = " << variable_weight[dec] << " \n";
// #endif

// 	} else {

// #ifdef _DEBUG_IMPACT
// 	  std::cout << " ==> right-weight[" << solver->variables[nxt] << "] was " << right_weight[nxt] ;
// 	  //std::cout << " [" << dec << "]";
// 	  std::cout << " (tot=" << residual_space << ")";
// #endif

// 	  if(dec>=0) 
// 	    residual_space /= left_weight[dec];

// #ifdef _DEBUG_IMPACT
// 	  std::cout << " (/" << "lw[x" << dec << "/" << left_weight[dec] << "]=" << residual_space << ")";
// #endif

// 	  for(int i=0; i<num_right_branches[nxt]-1; ++i)
// 	    residual_space /= right_weight[nxt];

// #ifdef _DEBUG_IMPACT
// 	  std::cout << " (slf=" << residual_space << ")";
// #endif

// 	  right_weight[nxt] = (((double)(num_right_probes[nxt]) * right_weight[nxt]) + residual_space)/(double)(++num_right_probes[nxt]);

// #ifdef _DEBUG_IMPACT
// 	  std::cout << " now " << right_weight[nxt] << std::endl;
// #endif

// 	  variable_weight[nxt] = avg_right_branches[nxt] * (left_weight[nxt] + right_weight[nxt]);
	  
// #ifdef _DEBUG_IMPACT
// 	  std::cout << " ===> total weight of " << solver->variables[nxt] << " = " << avg_right_branches[nxt] << " * (" << left_weight[nxt] << " + " << right_weight[nxt] << ") = " << variable_weight[nxt] << " \n";
// #endif

// 	}
//       }

// #ifdef _DEBUG_IMPACT
//       else {
// 	std::cout << "initial propagate!!\n";
//       }
// #endif

//       left = 0;
//     }

//     virtual void notify_backtrack() {
//       // propagation produced a wipe-out
//       // - check if it was after a left or a right branch
//       // - find out what was the decision/refutation

//       int dec;
//       //if(!solver->decisions.empty()) {
//       if(left==1) {
// 	// left branch
// 	dec = solver->decisions.back().var.id();

// #ifdef _DEBUG_IMPACT
// 	std::cout << " notified of backtrack after a left branch on " << solver->variables[dec] << "\n";
// 	std::cout << " ==> left-weight[" << solver->variables[dec] << "] was " << left_weight[dec] ;
// #endif


// 	left_weight[dec] = (((double)(num_left_probes[dec]) * left_weight[dec]))/(double)(++num_left_probes[dec]);

// #ifdef _DEBUG_IMPACT
// 	std::cout << " now " << left_weight[dec] << std::endl;
// #endif

// 	variable_weight[dec] = avg_right_branches[dec] * (left_weight[dec] + right_weight[dec]);
	  
// #ifdef _DEBUG_IMPACT
// 	  std::cout << " ===> total weight of " << solver->variables[dec] << " = " << avg_right_branches[dec] << " * (" << left_weight[dec] << " + " << right_weight[dec] << ") = " << variable_weight[dec] << " \n";
// #endif

//       } else if(left==0) {
// 	// right branch
// 	dec = solver->decisions.back(0).var.id();

// #ifdef _DEBUG_IMPACT
// 	double rwb = right_weight[dec];
// 	double nbr = avg_right_branches[dec];
// #endif

// 	right_weight[dec] = (((double)(num_right_probes[dec]) * right_weight[dec]))/(double)(++num_right_probes[dec]);

// 	++num_right_branches[dec];
// 	avg_right_branches[dec] = (avg_right_branches[dec]*(double)(num_probes[dec]) + (double)(num_right_branches[dec]))/((double)(++num_probes[dec]));

// #ifdef _DEBUG_IMPACT
// 	std::cout << " notified of backtrack after the " << (num_right_branches[dec]) 
// 		  << "th right branch on " << solver->variables[dec] << "\n";
// 	std::cout << " ==> right-weight[" << solver->variables[dec] << "] was " << rwb
// 		  << " now " << right_weight[dec] << std::endl;
// 	std::cout << " ==> num branches on " << solver->variables[dec] << " was " << nbr
// 		  << " now " << avg_right_branches[dec] << std::endl;
// #endif

// 	num_right_branches[dec] = 0;

// 	variable_weight[dec] = avg_right_branches[dec] * (left_weight[dec] + right_weight[dec]);

// #ifdef _DEBUG_IMPACT
// 	  std::cout << " ===> total weight of " << solver->variables[dec] << " = " << avg_right_branches[dec] << " * (" << left_weight[dec] << " + " << right_weight[dec] << ") = " << variable_weight[dec] << " \n";
// #endif

// 	//}
//       } 

// #ifdef _DEBUG_IMPACT
//       else {
// 	std::cout << "initial propagate!!\n";
//       }
// #endif

//       left = 0;
//     }

 




//     virtual std::ostream& display(std::ostream& os, const bool all) const ;
// // {
      
// //       int *all_variables = new int[variable_weight.size];
// //       int *all_constraints = new int[constraint_weight.size];

// //       for(unsigned int i=0; i<variable_weight.size; ++i) {
// // 	all_variables[i] = i;
// //       }

// //       for(unsigned int i=0; i<constraint_weight.size; ++i) {
// // 	all_constraints[i] = i;
// //       }

// //       weight_sorting_array = variable_weight.stack_;
// //       qsort(all_variables, variable_weight.size, sizeof(int), decreasing_weight);

// //       weight_sorting_array = constraint_weight.stack_;
// //       qsort(all_constraints, constraint_weight.size, sizeof(int), decreasing_weight);

// //       for(unsigned int i=0; i<variable_weight.size; ++i) {
// // 	os << std::setw(5) << solver->variables[all_variables[i]] << " ";
// //       }
// //       os << std::endl;
// //       for(unsigned int i=0; i<variable_weight.size; ++i) {
// // 	os << std::setw(5) << variable_weight[i] << " ";
// //       }
// //       os << std::endl;

// //      for(unsigned int i=0; i<constraint_weight.size; ++i) {
// // 	os << std::setw(5) << solver->constraints[all_constraints[i]] << " ";
// //       }
// //       os << std::endl;
// //       for(unsigned int i=0; i<constraint_weight.size; ++i) {
// // 	os << std::setw(5) << constraint_weight[i] << " ";
// //       }
// //       os << std::endl;

// //     }    

//   };





  /*! \class ImpactManager
    \brief ImpactManager Class

    * Listener interface for Impact *
    * NB: this is a simplified version of Impact, where only the impact of left vs right branches are distinguished
  */
#define INIT_IMPACT .001
  class ImpactManager : public BacktrackListener, public SuccessListener, public DecisionListener, public VariableListener {

  public:

    Solver *solver;
    double weight_unit;
    int left;

    Vector<double> variable_weight;
    Vector<double> impact;
    Vector<double> avg_branches;
    Vector<int> num_probes;
    Vector<int> tot_probes;
    Vector<int> tot_fails;


    ImpactManager(Solver *s) : solver(s) {// }

      left = -1;
      weight_unit = solver->parameters.activity_increment;
      
      // actual weight put on the variable (for xi)
      variable_weight.initialise(solver->variables.size, solver->variables.size);

      // average impact (stored as the ratio size_after/size_before) (for xi)
      impact.initialise(solver->variables.size, solver->variables.size);

      // total number of probes (on xi)
      tot_probes.initialise(solver->variables.size, solver->variables.size); 

      // number of probes since the last failure (on xi)
      num_probes.initialise(solver->variables.size, solver->variables.size);

      // total number of failures (on xi)
      tot_fails.initialise(solver->variables.size, solver->variables.size); 
      
      // average number of branches explored (for xi)
      avg_branches.initialise(solver->variables.size, solver->variables.size);


      for(unsigned int i=0; i<solver->variables.size; ++i) {
	avg_branches[i] = solver->variables[i].get_size();
	impact[i] = INIT_IMPACT/(double)(solver->variables[i].get_degree());
	num_probes[i] = 1;
	tot_probes[i] = 0;
	tot_fails[i] = 0;
	variable_weight[i] = avg_branches[i] * impact[i] ;
      }

      solver->add((BacktrackListener*)this);
      solver->add((SuccessListener*)this);
      solver->add((DecisionListener*)this);
      solver->add((VariableListener*)this);
    }

    virtual ~ImpactManager() {// }

      solver->remove((VariableListener*)this);
      solver->remove((SuccessListener*)this);
      solver->remove((DecisionListener*)this);
      solver->remove((BacktrackListener*)this);
    }

    double *get_variable_weight() { return variable_weight.stack_; }   
    double **get_value_weight() { return NULL; }
    double *get_bound_weight() { return NULL; }

    virtual void check_consistency() {
    }

    virtual void notify_add_var() {
#ifdef _DEBUG_IMPACT
      std::cout << "NEW VAR!!\n";
#endif
    }

    virtual void notify_change(const int id) {
#ifdef _DEBUG_IMPACT
      std::cout << solver->variables[id] << "'s domain implementation has changed!!\n";
#endif
    }

    virtual void notify_decision() {
      left = 1;
    }
  
    virtual void notify_success() {
      // propagation went without wipe-out
      // - check if it was after a left or a right branch
      // - find out what was the decision/refutation
      int dec;
      int id;
      int i, n;
      Variable x;
      double residual_space;
      int size;
      if(left==1) {
	i = solver->trail_.back(5), n=solver->saved_vars.size;
	if(left) {
	  // left branch
	  dec = solver->decisions.back().var.id();

#ifdef _DEBUG_IMPACT
	  std::cout << " notified of success after the " << num_probes[dec] << "th left branch on " << solver->variables[dec] << "\n";
#endif

	  residual_space = 1.0;
	  while(i<n) {	
	    id = solver->saved_vars[i];
	    x = solver->variables[id];
	    
#ifdef _DEBUG_IMPACT
	    std::cout << " -> " << x << x.get_domain() << " lost " << x.get_reduction() << " values (";
#endif
	    
	    size = x.get_size();
	    residual_space *= (((double)size))/((double)(size+x.get_reduction()));
	    
#ifdef _DEBUG_IMPACT
	    std::cout << residual_space << ")\n";
#endif
	    
	    ++i;
	  } 

#ifdef _DEBUG_IMPACT
	  std::cout << " ==> impact[" << solver->variables[dec] << "] was " << impact[dec] ;
#endif
	  
	  impact[dec] = (((double)(tot_probes[dec]) * impact[dec]) + residual_space)/(double)(++tot_probes[dec]);
	  
#ifdef _DEBUG_IMPACT
	  std::cout << " now " << impact[dec] << std::endl;
#endif
	  
	  ++num_probes[dec];
	  variable_weight[dec] = avg_branches[dec] * impact[dec];
	  
#ifdef _DEBUG_IMPACT
	  std::cout << " ===> total weight of " << solver->variables[dec] << " = " << avg_branches[dec] << " * " << impact[dec] << " = " << variable_weight[dec] << " \n";
#endif
	  
	}
      }
      
      left = 0;
    }
      
    virtual void notify_backtrack() {
      // propagation produced a wipe-out
      // - check if it was after a left or a right branch
      // - find out what was the decision/refutation

      int dec;
      //if(!solver->decisions.empty()) {
      if(left==1) {
	// left branch
	dec = solver->decisions.back().var.id();

#ifdef _DEBUG_IMPACT
	std::cout << " notified of backtrack after the " << num_probes[dec] << " left branch on " << solver->variables[dec] << "\n";
	std::cout << " ==> left-weight[" << solver->variables[dec] << "] was " << impact[dec] ;
#endif

	impact[dec] = (((double)(tot_probes[dec]) * impact[dec]))/(double)(++tot_probes[dec]);

#ifdef _DEBUG_IMPACT
	std::cout << " now " << impact[dec] << std::endl;
#endif

	variable_weight[dec] = avg_branches[dec] * impact[dec] ;
	++num_probes[dec];

#ifdef _DEBUG_IMPACT
	  std::cout << " ===> total weight of " << solver->variables[dec] << " = " << avg_branches[dec] << " * " << impact[dec] << " = " << variable_weight[dec] << " \n";
#endif

      } else if(left==0) {
	// right branch
	dec = solver->decisions.back(0).var.id();

#ifdef _DEBUG_IMPACT
	double nbr = avg_branches[dec];
#endif

	avg_branches[dec] = (avg_branches[dec]*(double)(tot_fails[dec]) + (double)(num_probes[dec]))/((double)(++tot_fails[dec]));

#ifdef _DEBUG_IMPACT
	std::cout << " notified of backtrack after a right branch on " << solver->variables[dec] << "\n";
	std::cout << " ==> num branches on " << solver->variables[dec] << " was " << nbr
		  << " now " << avg_branches[dec] << std::endl;
#endif

	variable_weight[dec] = avg_branches[dec] * impact[dec];

#ifdef _DEBUG_IMPACT
	  std::cout << " ===> total weight of " << solver->variables[dec] << " = " << avg_branches[dec] << " * " << impact[dec] << " = " << variable_weight[dec] << " \n";
#endif

	num_probes[dec] = 0;
      } 

      left = 0;
    }
    
    virtual std::ostream& display(std::ostream& os, const bool all) const ;
    
  };




  // /*! \class ProgressSavingManager
  //   \brief ProgressSavingManager Class

  //   * Listener interface for progress saving *
  // */
  // class ProgressSavingManager : public BacktrackListener {

  // public:

  //   Solver *solver;

  //   int max_solution_length;
  //   int best_objective;

  //   Vector<int> progress;
    

  //   //ProgressSavingManager(Solver *s, void *a=NULL) : solver(s) {
  //   ProgressSavingManager(Solver *s) : solver(s) {
  //     best_objective = objective->value();
  //     max_solution_length = 0;
  //     for(unsigned int i=0; i<solver->variables.size; ++i) {
  // 	progress.add(solver->variables[i].get_min());
  //     }
  //     solver->add((BacktrackListener*)this);
  //   }

  //   virtual ~ProgressSavingManager() {
  //     solver->remove((BacktrackListener*)this);
  //   }

  //    virtual void notify_backtrack() {

  //   }
  // };


  // /*! \class GuidedSearchManager
  //   \brief GuidedSearchManager Class

  //   * Listener interface for progress saving *
  // */
  // class GuidedSearchManager : public BacktrackListener {

  // public:

  //   Solver *solver;

  //   int best_objective;

  //   Vector<int> best_solution;
    

  //   //GuidedSearchManager(Solver *s, void *a=NULL) : solver(s) {
  //   GuidedSearchManager(Solver *s) : solver(s) {
  //     best_objective = objective->value();
  //     max_solution_length = 0;
  //     for(unsigned int i=0; i<solver->variables.size; ++i) {
  // 	progress.add(solver->variables[i].get_min());
  //     }
  //     solver->add((BacktrackListener*)this);
  //   }

  //   virtual ~GuidedSearchManager() {
  //     solver->remove((BacktrackListener*)this);
  //   }

  //    virtual void notify_backtrack() {

  //   }
  // };


  /*! \class RestartPolicy
    \brief  Interface RestartPolicy

    super class for restart-cutoff sequence generators
  */
  class RestartPolicy {
    
  public:

    unsigned int base;
    
    RestartPolicy(const unsigned int b=256);
    virtual ~RestartPolicy();

    virtual void reset(unsigned int& limit) = 0;
    virtual void initialise(unsigned int& limit) = 0;
    
  };


  class NoRestart : public RestartPolicy {
    
  public:
    
    NoRestart();
    virtual ~NoRestart();
    
    void reset(unsigned int& limit) {
      limit = base;
    }

    void initialise(unsigned int& limit) {
      limit = base;
    }
    
  };


  class Geometric : public RestartPolicy {
    
  public:
    
    unsigned int increment;
    double factor;

    Geometric(const unsigned int b=256, const double f=1.333);
    virtual ~Geometric();
    
    void reset(unsigned int& limit) {
      limit += increment;
      increment = (unsigned int)((double)increment * factor);
    }

    void initialise(unsigned int& limit) {
      limit = 0;
      increment = base;
      reset(limit);
    }
    
  };

  class Luby : public RestartPolicy {

  private:
    
    unsigned int luby_seq(const int iter) {
      unsigned int thelog = log2_(iter);
      if( iter == (1 << (thelog + 1))-1 )
	return (1 << thelog);
      return luby_seq(iter - (1 << thelog) + 1);
    }
    
    
  public:
    
    unsigned int iteration;

    Luby(const unsigned int b=100);
    virtual ~Luby();
    
    void reset(unsigned int& limit) {
      // unsigned int increment = (base * luby_seq(++iteration));
      // std::cout << "restart for " << increment << std::endl;
      // limit += increment;

      limit += (base * luby_seq(++iteration));
    }

    void initialise(unsigned int& limit) {
      iteration = 0;
      reset(limit);
    }
    
    
  };


  /**********************************************
   * Search Strategies
   **********************************************/

  /*! \class BranchingHeuristic
    \brief  Interface BranchingHeuristic

    The branching heuristic is queried by the
    solver on each new node in order to take
    a branching decision.

    It implements one method: branch(), which return a 
    Decision object representing the branching decision
    to make.
  */
  class BranchingHeuristic {

  public:
    
    Solver *solver;

    BranchingHeuristic() {}
    BranchingHeuristic(Solver *s) {solver = s;}
    virtual ~BranchingHeuristic() {}

    //virtual void initialise() {}
    virtual void initialise(VarStack< Variable, ReversibleNum<int> >& seq) {}
    virtual void close() {}

    virtual Decision branch() = 0;

    virtual std::ostream& display(std::ostream& os) = 0;
  };


// std::ostream& operator<<(std::ostream& os, BranchingHeuristic& x) {
//   return x.display(os);
// }

// std::ostream& operator<<(std::ostream& os, BranchingHeuristic* x) {
//   return x->display(os);
// }

  
  //std::ostream& operator<<(std::ostream& os, BranchingHeuristic& x);
  //std::ostream& operator<<(std::ostream& os, BranchingHeuristic* x);


  /**********************************************
   * Generic Variable/Value Ordering heuristics
   **********************************************/
  /*! \class GenericHeuristic
    \brief  Class GenericHeuristic

    Generic branching heuristic: 
    1/ a method to select the next variable to branch on (var)
    2/ a method to reduce the domain of this variable (choice)
  */
  template < class VarSelector, class ValSelector >
  class GenericHeuristic : public BranchingHeuristic {
  public:

    VarSelector var;
    ValSelector choice;

    GenericHeuristic(Solver *s) 
      : BranchingHeuristic(s) {

      //std::cout << (int*)(s) << std::endl;
      //std::cout << s << std::endl;

      var.initialise(s);
      choice = ValSelector(s,var.get_value_weight(),var.get_bound_weight());
    }

    virtual ~GenericHeuristic() {}

    // GenericHeuristic(Solver *s, void *a) 
    //   : BranchingHeuristic(s) {
    //   var.initialise(s,a);
    //   choice = ValSelector(s,a);
    // }

    virtual void initialise(VarStack< Variable, ReversibleNum<int> >& seq) {var.initialise(seq);}

    virtual Decision branch() {
      return choice.make(var.select());
    }

    virtual std::ostream& display(std::ostream& os) {
      //os << "Branch on the variable " << std::endl;
      var.display(os);
      os << " c ";
      choice.display(os);
      os << " " << std::endl ;
      return os;
    }
  };


  template < class VarSelector, class ValSelector >
  std::ostream& operator<<(std::ostream& os, BranchingHeuristic* x) {
    return x->display(os);
  }



  template< class VarComparator >
  class Identifiable {
    
  public :
    
    VarComparator criterion;
    int id;
    
    /**@name Utils*/
    //@{
    inline double value() { return criterion.value(); } 
    inline bool operator<( const Identifiable<VarComparator>& x ) const { 
      if(criterion < x.criterion) return true; 
      else if(x.criterion < criterion) return false;
      return (id > x.id);
    }

    inline bool operator>( const Identifiable<VarComparator>& x ) const { 
      if(criterion < x.criterion) return false; 
      else if(x.criterion < criterion) return true;
      return (id < x.id);
    }
    inline bool operator==( const Identifiable<VarComparator>& x ) const { 
      return(criterion == x.criterion && id == x.id);
    }
    inline void operator=( const Identifiable<VarComparator>& x ) { criterion = x.criterion; id = x.id; }
    inline void operator=( const Variable x ) { criterion = x; id = x.id(); }
    //@}  

    std::ostream& display_criterion(std::ostream& os) const {
      return criterion.display_criterion(os);
    }

    std::ostream& display(std::ostream& os) const {
      return criterion.display(os);
    }

  };




  /*! \class GenericRandomDVO
    \brief  Class GenericRandomDVO

    Randomized Generic (dynamic) variable ordering heuristic 
    - A parameterized comparison method is used to select 
      the k 'best' variables, then one of them is randomly selected
  */
  template < class VarComparator, int RAND = 1, class WeightManager = NoManager >
  class GenericDVO 
  {
  public: 

    /**@name Parameters*/
    //@{ 
    Solver         *solver;
    WeightManager *manager;
    VarComparator  current;
    VarComparator    bests[RAND+1];
    Variable      bestvars[RAND+1];
    //@}

    /**@name Constructors*/
    //@{
    GenericDVO() { solver = NULL; manager = NULL; }
    GenericDVO(Solver* s) : solver(s) { manager = NULL; }
   
    virtual void initialise(Solver *s) { 
      solver = s; 
      initialise_manager();
    }
    virtual void initialise(VarStack< Variable, ReversibleNum<int> >& seq) {}
    virtual void initialise_manager() {

      //std::cout << "initialise manager (" << manager << ") \n";

      if(!manager) {
	manager = new WeightManager(solver);
	current.initialise(manager->get_variable_weight());
	for(int i=0; i<=RAND; ++i)
	  bests[i].initialise(manager->get_variable_weight());
      }
    }

    virtual ~GenericDVO() { delete manager; }
    //@}

    /**@name Utils*/
    //@{ 
    double *get_variable_weight() { return manager->get_variable_weight(); }  
    double **get_value_weight() { return manager->get_value_weight(); }
    double *get_bound_weight() { return manager->get_bound_weight(); }


    Variable select()
    {

#ifdef _DEBUG_VARORD
      std::cout << std::endl;
#endif

      Variable *variables = solver->sequence.list_;
      unsigned int length = solver->sequence.size-1;
      unsigned int realsize=1, i, j;
      bests[0] = bestvars[0] = variables[length];
      for(j=length; j--;)
	{  
	  current = variables[j];
	  i = realsize;
	  while( i && current < bests[i-1] ) {
	    bests[i] = bests[i-1];
	    bestvars[i] = bestvars[i-1];
	    --i;
	  }

#ifdef _DEBUG_VARORD
	  if(i<RAND) {
	    std::cout << "*";
	  }
	  std::cout << std::endl;
#endif

	  bests[i] = current;
	  bestvars[i] = variables[j];
	  
	  if(realsize<RAND) ++realsize;
	}
      return bestvars[(realsize>1 ? randint(realsize) : 0)];
    }
    //@}



    virtual std::ostream& display(std::ostream& os) const { //,  const int n, double* weights) {

      manager->display(os, false);

      //double* weights = (manager ? manager->get_variable_weight() : NULL);
      //int n = RAND;

      os << " c Select the " ;
      if(RAND>1) os << RAND << " ";
      os << "variable" << (RAND > 1 ? "s " : " ") ;//<< "with minimal value of ";
      VarComparator v;
      v.display_criterion(os);
      if(RAND>1) os << " and pick one uniformly at random" ;
      os << std::endl;

      if(solver->sequence.size > 1) {
	Variable *variables = solver->sequence.list_;
	unsigned int length = solver->sequence.size-1;
	Variable var = variables[length];
	
	
	os << "--> branch in [";
	
	std::vector< Identifiable< VarComparator > > all_vars;
	for(unsigned int i=0; i<=length; ++i) {
	  Identifiable<VarComparator> vc;
	  //if(weights) vc.criterion.weight = weights;
	  vc = variables[i];
	  vc.id = i;
	  all_vars.push_back(vc);
	  
	  os << variables[i] << " in " << variables[i].get_domain() << " ";
	  
	}
	
	os << std::endl;


	sort(all_vars.begin(), all_vars.end());
	
	os << " c [" << variables[all_vars[0].id].id() << ":";
	all_vars[0].display(os);
	
	for(int i=1; i<RAND; ++i) {
	  os << " " << variables[all_vars[i].id].id() << ":";
	  all_vars[i].display(os);
	}
	
	os << "]";

	for(unsigned int i=RAND; i<all_vars.size(); ++i) {
	  os << " " << variables[all_vars[i].id].id() << ":";
	  all_vars[i].display(os);
	}
	os << std::endl;
      }
      return os;
    }

  };





  /*! \class GenericDVO
    \brief  Class GenericDVO

    Generic (dynamic) variable ordering heuristic with a weight function. 
    A 'WeightManager' is used to compute each variable's weight.
    This is the randomized version. i.e., the k best variables are computed
    then one of them is randomly selected
  */
  template <template< class T > class Aggregator, class VarComparator, int RAND = 1, class WeightManager = NoManager >
  class GenericNeighborDVO : public GenericDVO< Aggregator< VarComparator >, RAND, WeightManager >
  {
  public: 
    
    /**@name Parameters*/
    //@{ 
    Solver *solver;
    Vector< Variable > *neighborhood;
    //@}

    /**@name Constructors*/
    //@{
    GenericNeighborDVO() : GenericDVO< Aggregator< VarComparator >, RAND, WeightManager >() { neighborhood=NULL; }
    GenericNeighborDVO(Solver* s) : GenericDVO< Aggregator< VarComparator >, RAND, WeightManager >(s) {
      neighborhood = NULL;
    }
    virtual ~GenericNeighborDVO() {
      delete [] neighborhood;
    }
    //virtual void initialise(Solver *s, void *a=NULL) { 
    virtual void initialise(Solver *s) { 
      GenericDVO< Aggregator< VarComparator >, RAND, WeightManager >::initialise(s);
    }
    //virtual void initialise(Solver *s) { manager->initialise(s); }
    virtual void initialise(VarStack< Variable, ReversibleNum<int> >& seq) 
    {
      int i, j, k, cons, self_idx;
      Constraint constraint;
      Variable *scope;
      Event trig;
      bool is_in;

      //std::cout << "NEIGHBORHOOD " << neighborhood << std::endl;

      if(!neighborhood) {

	//std::cout << "NEIGHBORHOOD" << std::endl;

	int n = GenericDVO< Aggregator< VarComparator >, RAND, WeightManager >::solver->variables.size;
	neighborhood = new Vector< Variable >[n];
	std::fill(neighborhood, neighborhood+n, NULL);
	for(i=seq.size; --i>=0;) {

	  self_idx = seq[i].id();
	  for(trig = 0; trig<3; ++trig) {
	    
	    for(cons = GenericDVO< Aggregator< VarComparator >, RAND, WeightManager >::solver->constraint_graph[self_idx].on[trig].size; --cons>=0;) {
	      
	      constraint = GenericDVO< Aggregator< VarComparator >, RAND, WeightManager >::solver->constraint_graph[self_idx].on[trig][cons];

	      for(scope=constraint.get_scope(), j=constraint.arity(); --j>=0;) if(scope[j].id() != seq[i].id()) {
		  is_in = false;
		  for(k = neighborhood[self_idx].size; --k>=0 && !is_in;)
		    if(neighborhood[self_idx][k].id() == scope[j].id()) is_in = true;
		  if(!is_in) neighborhood[self_idx].add(scope[j]);
		}
	    }
	  }	  //std::cout << seq[i] << ": " << neighborhood[self_idx] << std::endl;
	}
      }

      GenericDVO< Aggregator< VarComparator >, RAND, WeightManager >::current.map = neighborhood;
      for(int i=0; i<=RAND; ++i)
	GenericDVO< Aggregator< VarComparator >, RAND, WeightManager >::bests[i].map = neighborhood;
      
    }

    // virtual void initialise_manager() {
    //   manager = new WeightManager(GenericRandomDVO< Aggregator< VarComparator >, RAND >::solver);
    //   //GenericRandomDVO< VarComparator >::solver->add(manager);
    //   if(manager->get_variable_weight()) {
    // 	GenericRandomDVO< Aggregator< VarComparator >, RAND >::current.weight = manager->get_variable_weight();
    // 	for(int i=0; i<=RAND; ++i)
    // 	  GenericRandomDVO< Aggregator< VarComparator >, RAND >::bests[i].weight = manager->get_variable_weight();
    //   }
    // }

    virtual std::ostream& display(std::ostream& os) const {
      //GenericDVO< Aggregator< VarComparator >, RAND, WeightManager >::manager->display(os, false);
      return GenericDVO< Aggregator< VarComparator >, RAND, WeightManager >::display(os); //, RAND, GenericDVO< Aggregator< VarComparator >, RAND, WeightManager >::manager->get_variable_weight());
    }
    
    //@}
    
  };



  /*! \class NoOrder
    \brief  Class NoOrder

    This heuristic selects the first variable in the sequence, 
    there is no way of knowing which it will be (in particular
    this is NOT a lexicographic heuristic, and this is not
    even a static heuristic)
  */
  class NoOrder {

  public: 
    Solver *solver;

    NoOrder() { solver = NULL; }
    NoOrder(Solver *s);
    virtual ~NoOrder();
    void initialise(Solver *s) { solver = s; }
    //void initialise(Solver *s, void *a) { solver = s; }
    void initialise(VarStack< Variable, ReversibleNum<int> >& seq) {}

    double *get_variable_weight() { return NULL; }
    double **get_value_weight() { return NULL; }
    double *get_bound_weight() { return NULL; }

    
    Variable select();

    virtual std::ostream& display(std::ostream& os) const {
      os << "Go by the default sequence: " << solver->sequence.back(); 
      return os;
    }
  };


  /*! \class Lexicographic
    \brief  Class Lexicographic

    This heuristic selects the variable with lowest rank
    in the initial sequence of search variables.
  */
  class Lexicographic : public VariableListener {

  public: 
    
    Solver             *solver;
    Vector< Variable >   order;
    Vector< int >        index;
    ReversibleNum< int >  last;

    Lexicographic() { solver = NULL; }
    Lexicographic(Solver *s);
    void initialise(Solver *s);
    //void initialise(Solver *s, void *a=NULL);
    void initialise(VarStack< Variable, ReversibleNum<int> >& seq);
    virtual ~Lexicographic();

    virtual void notify_add_var() {};
    virtual void notify_change(const int idx);

    double *get_variable_weight() { return NULL; }
    double **get_value_weight() { return NULL; }
    double *get_bound_weight() { return NULL; }
    
    Variable select();

    virtual std::ostream& display(std::ostream& os) const {
      os << "Go by lexicographic order: " ;

      int i = last;
      while(i<(int)(order.size) && order[i].is_ground()) { 
	++i;
      }
      os << order[i];
      return os;
    }

  };


  /*! \class Lexicographic
    \brief  Class Lexicographic

    This heuristic selects the variable with lowest rank
    in the initial sequence of search variables.
  */
  class ConsolidateListener : public VariableListener, public ConstraintListener {

  public: 
    
    Solver                                      *solver;
    VarStack < Variable, ReversibleNum<int> > *sequence;
    Vector< Vector< Constraint > >          constraints;
    int                                          id_obj;

    ConsolidateListener(Solver *s);
    virtual ~ConsolidateListener();

    virtual void notify_post   (Constraint c);
    virtual void notify_relax  (Constraint c);
    virtual void notify_add_con(Constraint c);
    virtual void notify_add_var();
    virtual void notify_change (const int idx);

  };


  /**********************************************
   * Variable Comparators
   **********************************************/



  /*! \class MinDomain
    \brief  Class MinDomain

    Order two variables by their domain sizes
  */
  class MinDomain 
  {
  public: 

    /**@name Constructors*/
    //@{
    MinDomain() {dom_ = LARGE_VALUE;}
    void initialise(const double* _w) {dom_ = LARGE_VALUE;}
    //@}

    /**@name Parameters*/
    //@{ 
    int dom_;
    //@}  

    /**@name Utils*/
    //@{
    inline double value() { return (double)dom_; } 
    inline bool operator<( const MinDomain& x ) const { return dom_ < x.dom_; }
    inline bool operator>( const MinDomain& x ) const { return dom_ > x.dom_; }
    inline bool operator==( const MinDomain& x ) const { return dom_ == x.dom_; }

    // inline MinDomain& operator*( const int x ) const { MinDomain d; d.dom_=(dom_ * x); return d; }
    // inline MinDomain& operator+( const int x ) const { MinDomain d; d.dom_=(dom_ + x); return d; }
    // inline MinDomain& operator-( const int x ) const { MinDomain d; d.dom_=(dom_ - x); return d; }
    // inline MinDomain& operator/( const int x ) const { MinDomain d; d.dom_=(dom_ / x); return d; }

    inline MinDomain& operator*=( const int x ) { dom_ *= x; return *this; }
    inline MinDomain& operator+=( const int x ) { dom_ += x; return *this; }
    inline MinDomain& operator-=( const int x ) { dom_ -= x; return *this; }
    inline MinDomain& operator/=( const int x ) { dom_ /= x; return *this; }

    // inline MinDomain& operator*( const MinDomain& x ) const { MinDomain d; d.dom_=(dom_ * x.dom_); return d; }
    // inline MinDomain& operator+( const MinDomain& x ) const { MinDomain d; d.dom_=(dom_ + x.dom_); return d; }
    // inline MinDomain& operator-( const MinDomain& x ) const { MinDomain d; d.dom_=(dom_ - x.dom_); return d; }
    // inline MinDomain& operator/( const MinDomain& x ) const { MinDomain d; d.dom_=(dom_ / x.dom_); return d; }

    inline MinDomain& operator*=( const MinDomain& x ) { dom_ *= x.dom_; return *this; }
    inline MinDomain& operator+=( const MinDomain& x ) { dom_ += x.dom_; return *this; }
    inline MinDomain& operator-=( const MinDomain& x ) { dom_ -= x.dom_; return *this; }
    inline MinDomain& operator/=( const MinDomain& x ) { dom_ /= x.dom_; return *this; }

    inline void operator=( const MinDomain& x ) { dom_ = x.dom_; }
    inline void operator=( const Variable x ) { dom_ = x.get_size(); }

    //@}  

    std::ostream& display_criterion(std::ostream& os) const {
      os << "with minimum domain size";
      return os;
    }

    std::ostream& display(std::ostream& os) const {
      os << dom_;
      return os;
    }
  };

  std::ostream& operator<<(std::ostream& os, MinDomain& x);
  std::ostream& operator<<(std::ostream& os, MinDomain* x);


  /*! \class MinMin
    \brief  Class MinMin

    Order two variables by their minimum value
  */
  class MinMin 
  {
  public: 

    /**@name Constructors*/
    //@{
    MinMin() {min_ = LARGE_VALUE;}
    void initialise(const double* _w) {min_ = LARGE_VALUE;}
    //@}

    /**@name Parameters*/
    //@{ 
    int min_;
    //@}  

    /**@name Utils*/
    //@{
    inline double value() { return (double)min_; } 
    inline bool operator<( const MinMin& x ) const { return min_ < x.min_; }
    inline bool operator>( const MinMin& x ) const { return min_ > x.min_; }
    inline bool operator==( const MinMin& x ) const { return min_ == x.min_; }

    inline MinMin& operator*=( const int x ) { min_ *= x; return *this; }
    inline MinMin& operator+=( const int x ) { min_ += x; return *this; }
    inline MinMin& operator-=( const int x ) { min_ -= x; return *this; }
    inline MinMin& operator/=( const int x ) { min_ /= x; return *this; }

    inline MinMin& operator*=( const MinMin& x ) { min_ *= x.min_; return *this; }
    inline MinMin& operator+=( const MinMin& x ) { min_ += x.min_; return *this; }
    inline MinMin& operator-=( const MinMin& x ) { min_ -= x.min_; return *this; }
    inline MinMin& operator/=( const MinMin& x ) { min_ /= x.min_; return *this; }

    inline void operator=( const MinMin& x ) { min_ = x.min_; }
    inline void operator=( const Variable x ) { min_ = x.get_min(); }
    //@}  

    std::ostream& display_criterion(std::ostream& os) const {
      os << "with minimum mininimum";
      return os;
    }

    std::ostream& display(std::ostream& os) const {
      os << min_;
      return os;
    }
  };

  std::ostream& operator<<(std::ostream& os, MinMin& x);
  std::ostream& operator<<(std::ostream& os, MinMin* x);

  /*! \class MaxMax
      \brief  Class MaxMax

      Order two variables by their maximum value
    */
    class MaxMax
    {
    public:

      /**@name Constructors*/
      //@{
      MaxMax() {max_ = SMALL_VALUE;}
      void initialise(const double* _w) {max_ = SMALL_VALUE;}
      //@}

      /**@name Parameters*/
      //@{
      int max_;
      //@}

      /**@name Utils*/
      //@{
      inline double value() { return (double)max_; }
      inline bool operator<( const MaxMax& x ) const { return max_ > x.max_; }
      inline bool operator>( const MaxMax& x ) const { return max_ < x.max_; }
      inline bool operator==( const MaxMax& x ) const { return max_ == x.max_; }

      inline MaxMax& operator*=( const int x ) { max_ *= x; return *this; }
      inline MaxMax& operator+=( const int x ) { max_ += x; return *this; }
      inline MaxMax& operator-=( const int x ) { max_ -= x; return *this; }
      inline MaxMax& operator/=( const int x ) { max_ /= x; return *this; }

      inline MaxMax& operator*=( const MaxMax& x ) { max_ *= x.max_; return *this; }
      inline MaxMax& operator+=( const MaxMax& x ) { max_ += x.max_; return *this; }
      inline MaxMax& operator-=( const MaxMax& x ) { max_ -= x.max_; return *this; }
      inline MaxMax& operator/=( const MaxMax& x ) { max_ /= x.max_; return *this; }

      inline void operator=( const MaxMax& x ) { max_ = x.max_; }
      inline void operator=( const Variable x ) { max_ = x.get_max(); }
      //@}

      std::ostream& display_criterion(std::ostream& os) const {
        os << "with maximum maximum";
        return os;
      }

      std::ostream& display(std::ostream& os) const {
        os << max_;
        return os;
      }
    };

    std::ostream& operator<<(std::ostream& os, MaxMax& x);
    std::ostream& operator<<(std::ostream& os, MaxMax* x);


  /*! \class MaxRegret
    \brief  Class MaxRegret

    Order two variables by their minimum value
  */
  class MaxRegret 
  {
  public: 

    /**@name Constructors*/
    //@{
    MaxRegret() {reg_ = SMALL_VALUE;}
    void initialise(const double* _w) {reg_ = SMALL_VALUE;}
    //@}

    /**@name Parameters*/
    //@{ 
    int reg_;
    //@}  

    /**@name Utils*/
    //@{
    inline double value() { return (double)reg_; } 
    inline bool operator<( const MaxRegret& x ) const { return reg_ > x.reg_; }
    inline bool operator>( const MaxRegret& x ) const { return reg_ < x.reg_; }
    inline bool operator==( const MaxRegret& x ) const { return reg_ == x.reg_; }

    inline MaxRegret& operator*=( const int x ) { reg_ *= x; return *this; }
    inline MaxRegret& operator+=( const int x ) { reg_ += x; return *this; }
    inline MaxRegret& operator-=( const int x ) { reg_ -= x; return *this; }
    inline MaxRegret& operator/=( const int x ) { reg_ /= x; return *this; }

    inline MaxRegret& operator*=( const MaxRegret& x ) { reg_ *= x.reg_; return *this; }
    inline MaxRegret& operator+=( const MaxRegret& x ) { reg_ += x.reg_; return *this; }
    inline MaxRegret& operator-=( const MaxRegret& x ) { reg_ -= x.reg_; return *this; }
    inline MaxRegret& operator/=( const MaxRegret& x ) { reg_ /= x.reg_; return *this; }

    inline void operator=( const MaxRegret& x ) { reg_ = x.reg_; }
    inline void operator=( const Variable x ) { 
      int min_ = x.get_min(); 
      int nxt_ = x.next(min_);
      reg_ = (nxt_ - min_);
    }
    //@}  

    std::ostream& display_criterion(std::ostream& os) const {
      os << "with maximum regret";
      return os;
    }

    std::ostream& display(std::ostream& os) const {
      os << reg_;
      return os;
    }
  };

  std::ostream& operator<<(std::ostream& os, MaxRegret& x);
  std::ostream& operator<<(std::ostream& os, MaxRegret* x);


  /*! \class Anti
    \brief  Class Anti

    Make the opposite decision of the parameter comparator
  */
  template<class VarComparator>
  class Anti 
  {
  public: 

    /**@name Constructors*/
    //@{
    Anti() : crit() {}
    void initialise(const double* _w) {crit.initialise(_w);}
    //@}

    /**@name Parameters*/
    //@{ 
    VarComparator crit;
    //@}  

    /**@name Utils*/
    //@{
    inline double value() { return 1.0/(double)(crit.value()); } 
    inline bool operator<( const Anti<VarComparator>& x ) const { return crit > x.crit; }
    inline bool operator>( const Anti<VarComparator>& x ) const { return crit < x.crit; }
    inline bool operator==( const Anti<VarComparator>& x ) const { return crit == x.crit; }

    inline Anti<VarComparator>& operator*=( const int x ) { crit *= x; return *this; }
    inline Anti<VarComparator>& operator+=( const int x ) { crit += x; return *this; }
    inline Anti<VarComparator>& operator-=( const int x ) { crit -= x; return *this; }
    inline Anti<VarComparator>& operator/=( const int x ) { crit /= x; return *this; }

    inline Anti<VarComparator>& operator*=( const Anti<VarComparator>& x ) { crit *= x.crit; return *this; }
    inline Anti<VarComparator>& operator+=( const Anti<VarComparator>& x ) { crit += x.crit; return *this; }
    inline Anti<VarComparator>& operator-=( const Anti<VarComparator>& x ) { crit -= x.crit; return *this; }
    inline Anti<VarComparator>& operator/=( const Anti<VarComparator>& x ) { crit /= x.crit; return *this; }

    inline void operator=( const Anti<VarComparator>& x ) { crit = x.crit; }
    inline void operator=( const Variable x ) { crit = x; }
    //@}  

    std::ostream& display_criterion(std::ostream& os) const {
      os << "NOT<";
      crit.display_criterion(os);
      os << ">";
      return os;
    }

    std::ostream& display(std::ostream& os) const {
      os << "-" ; //<< crit ;
      crit.display(os);
      return os;
    }
  };

  template<class VarComparator>
  std::ostream& operator<<(std::ostream& os, Anti<VarComparator>& x) {
    return x.display(os);
  }

  template<class VarComparator>
  std::ostream& operator<<(std::ostream& os, Anti<VarComparator>* x) {
    return x->display(os);
  }



  /*! \class LexCombination
    \brief  Class LexCombination

    Combines two comparators lexicographically
  */
  template<class ComparatorA, class ComparatorB>
  class LexCombination 
  {
  public: 

    /**@name Constructors*/
    //@{
    LexCombination() : critA(), critB() {}
    void initialise(const double* _w) {critA.initialise(_w); critB.initialise(_w);}
    //@}

    /**@name Parameters*/
    //@{ 
    ComparatorA critA;
    ComparatorB critB;
    //@}  

    /**@name Utils*/
    //@{
    inline double value() { return critA.value() + critB.value()*0.00001; } 
    inline bool operator<( const LexCombination<ComparatorA, ComparatorB>& x ) const { return (critA < x.critA || (critA == x.critA && critB < x.critB)); }
    inline bool operator>( const LexCombination<ComparatorA, ComparatorB>& x ) const { return (critA > x.critA || (critA == x.critA && critB > x.critB)); }
    inline bool operator==( const LexCombination<ComparatorA, ComparatorB>& x ) const { return (critA == x.critA && critB == x.critB); }

    inline LexCombination<ComparatorA, ComparatorB>& operator*=( const int x ) { critA *= x; return *this; }
    inline LexCombination<ComparatorA, ComparatorB>& operator+=( const int x ) { critA += x; return *this; }
    inline LexCombination<ComparatorA, ComparatorB>& operator-=( const int x ) { critA -= x; return *this; }
    inline LexCombination<ComparatorA, ComparatorB>& operator/=( const int x ) { critA /= x; return *this; }

    inline LexCombination<ComparatorA, ComparatorB>& operator*=( const LexCombination<ComparatorA, ComparatorB>& x ) { critA *= x.critA; critB *= x.critB; return *this; }
    inline LexCombination<ComparatorA, ComparatorB>& operator+=( const LexCombination<ComparatorA, ComparatorB>& x ) { critA += x.critA; critB += x.critB; return *this; }
    inline LexCombination<ComparatorA, ComparatorB>& operator-=( const LexCombination<ComparatorA, ComparatorB>& x ) { critA -= x.critA; critB -= x.critB; return *this; }
    inline LexCombination<ComparatorA, ComparatorB>& operator/=( const LexCombination<ComparatorA, ComparatorB>& x ) { critA /= x.critA; critB /= x.critB; return *this; }

    inline void operator=( const LexCombination<ComparatorA, ComparatorB>& x ) { critA = x.critA; critB = x.critB; }
    inline void operator=( const Variable x ) { critA = x; critB = x; }
    //@}  

    std::ostream& display_criterion(std::ostream& os) const {
      os << "<";
      critA.display_criterion(os);
      os << "->";
      critB.display_criterion(os);
      os << ">";
      return os;
    }

    std::ostream& display(std::ostream& os) const {
      //os <<  ; //<< crit ;
      critA.display(os);
      os << "->";
      critB.display(os);
      return os;
    }
  };

  template<class ComparatorA, class ComparatorB>
  std::ostream& operator<<(std::ostream& os, LexCombination<ComparatorA, ComparatorB>& x) {
    return x.display(os);
  }

  template<class ComparatorA, class ComparatorB>
  std::ostream& operator<<(std::ostream& os, LexCombination<ComparatorA, ComparatorB>* x) {
    return x->display(os);
  }



  // /*! \class DivCombination
  //   \brief  Class DivCombination

  //   Combines two comparators by dividing the former by the later
  // */
  // template<class ComparatorA, class ComparatorB>
  // class DivCombination 
  // {
  // public: 

  //   /**@name Constructors*/
  //   //@{
  //   DivCombination() : crit() {}
  //   void initialise(const double* _w) {critA.initialise(_w); critB.initialise(_w);}
  //   //@}

  //   /**@name Parameters*/
  //   //@{ 
  //   ComparatorA critA;
  //   ComparatorB critB;
  //   //@}  

  //   /**@name Utils*/
  //   //@{
  //   inline double value() { return critA.value() / critB.value(); } 
  //   inline bool operator<( const LexCombination<ComparatorA, ComparatorB>& x ) const { return critA * x
  //   inline bool operator>( const LexCombination<ComparatorA, ComparatorB>& x ) const { return (critA > x.critA || (critA == x.critA && critB > x.critB)); }
  //   inline bool operator==( const LexCombination<ComparatorA, ComparatorB>& x ) const { return (critA == x.critA && critB == x.critB); }

  //   inline LexCombination<ComparatorA, ComparatorB>& operator*=( const int x ) { crit *= x; return *this; }
  //   inline LexCombination<ComparatorA, ComparatorB>& operator+=( const int x ) { crit += x; return *this; }
  //   inline LexCombination<ComparatorA, ComparatorB>& operator-=( const int x ) { crit -= x; return *this; }
  //   inline LexCombination<ComparatorA, ComparatorB>& operator/=( const int x ) { crit /= x; return *this; }

  //   inline LexCombination<ComparatorA, ComparatorB>& operator*=( const LexCombination<ComparatorA, ComparatorB>& x ) { crit *= x.crit; return *this; }
  //   inline LexCombination<ComparatorA, ComparatorB>& operator+=( const LexCombination<ComparatorA, ComparatorB>& x ) { crit += x.crit; return *this; }
  //   inline LexCombination<ComparatorA, ComparatorB>& operator-=( const LexCombination<ComparatorA, ComparatorB>& x ) { crit -= x.crit; return *this; }
  //   inline LexCombination<ComparatorA, ComparatorB>& operator/=( const LexCombination<ComparatorA, ComparatorB>& x ) { crit /= x.crit; return *this; }

  //   inline void operator=( const LexCombination<ComparatorA, ComparatorB>& x ) { crit = x.crit; }
  //   inline void operator=( const Variable x ) { crit = x; }
  //   //@}  

  //   std::ostream& display_criterion(std::ostream& os) const {
  //     os << "<";
  //     critA.display_criterion(os);
  //     os << "->"
  //     critB.display_criterion(os);
  //     os << ">";
  //     return os;
  //   }

  //   std::ostream& display(std::ostream& os) const {
  //     //os <<  ; //<< crit ;
  //     critA.display(os);
  //     os << "->"
  //     critB.display(os);
  //     return os;
  //   }
  // };

  // template<class ComparatorA, ComparatorB>
  // std::ostream& operator<<(std::ostream& os, LexCombination<ComparatorA, ComparatorB>& x) {
  //   return x.display(os);
  // }

  // template<class ComparatorA, ComparatorB>
  // std::ostream& operator<<(std::ostream& os, LexCombination<ComparatorA, ComparatorB>* x) {
  //   return x.display(os);
  // }



  /*! \class MaxDegree
    \brief  Class MaxDegree

    Order two variables by their domain sizes
  */
  class MaxDegree 
  {
  public: 

    /**@name Constructors*/
    //@{
    MaxDegree() {deg_ = LARGE_VALUE;}
    void initialise(const double* _w) {deg_ = LARGE_VALUE;}
    //@}

    /**@name Parameters*/
    //@{ 
    int deg_;
    //@}  

    /**@name Utils*/
    //@{
    inline double value() { return (double)deg_; } 
    inline bool operator<( const MaxDegree& x ) const { return deg_ > x.deg_; }
    inline bool operator>( const MaxDegree& x ) const { return deg_ < x.deg_; }
    inline bool operator==( const MaxDegree& x ) const { return deg_ == x.deg_; }

    // inline MaxDegree& operator*( const int x ) const { MaxDegree d; d.deg_=(deg_ * x); return d; }
    // inline MaxDegree& operator+( const int x ) const { MaxDegree d; d.deg_=(deg_ + x); return d; }
    // inline MaxDegree& operator-( const int x ) const { MaxDegree d; d.deg_=(deg_ - x); return d; }
    // inline MaxDegree& operator/( const int x ) const { MaxDegree d; d.deg_=(deg_ / x); return d; }

    inline MaxDegree& operator*=( const int x ) { deg_ *= x; return *this; }
    inline MaxDegree& operator+=( const int x ) { deg_ += x; return *this; }
    inline MaxDegree& operator-=( const int x ) { deg_ -= x; return *this; }
    inline MaxDegree& operator/=( const int x ) { deg_ /= x; return *this; }

    // inline MaxDegree& operator*( const MaxDegree& x ) const { MaxDegree d; d.deg_=(deg_ * x.deg_); return d; }
    // inline MaxDegree& operator+( const MaxDegree& x ) const { MaxDegree d; d.deg_=(deg_ + x.deg_); return d; }
    // inline MaxDegree& operator-( const MaxDegree& x ) const { MaxDegree d; d.deg_=(deg_ - x.deg_); return d; }
    // inline MaxDegree& operator/( const MaxDegree& x ) const { MaxDegree d; d.deg_=(deg_ / x.deg_); return d; }

    inline MaxDegree& operator*=( const MaxDegree& x ) { deg_ *= x.deg_; return *this; }
    inline MaxDegree& operator+=( const MaxDegree& x ) { deg_ += x.deg_; return *this; }
    inline MaxDegree& operator-=( const MaxDegree& x ) { deg_ -= x.deg_; return *this; }
    inline MaxDegree& operator/=( const MaxDegree& x ) { deg_ /= x.deg_; return *this; }

    inline void operator=( const MaxDegree& x ) { deg_ = x.deg_; }
    inline void operator=( const Variable x ) { deg_ = x.get_degree(); }

    //@}  

    std::ostream& display_criterion(std::ostream& os) const {
      os << "with maximum degree";
      return os;
    }

    std::ostream& display(std::ostream& os) const {
      os << deg_;
      return os;
    }
  };

  std::ostream& operator<<(std::ostream& os, MaxDegree& x);

  std::ostream& operator<<(std::ostream& os, MaxDegree* x);


  // /*! \class MinDomainMaxDegree
  //   \brief  Class MinDomainMaxDegree

  //   Order two variables by their domain sizes, with ties broken by maximum degree
  // */
  // class MinDomainMaxDegree 
  // {
  // public: 

  //   /**@name Constructors*/
  //   //@{
  //   MinDomainMaxDegree() {dom_ = LARGE_VALUE; deg_ = 0;}
  //   void initialise(const double* _w) {dom_ = LARGE_VALUE; deg_ = 0;}
  //   //@}

  //   /**@name Parameters*/
  //   //@{ 
  //   int dom_;
  //   int deg_;
  //   //@}  

  //   /**@name Utils*/
  //   //@{
  //   inline double value() { return (deg_ ? (double)dom_ : (double)dom_ + 1.0/(double)deg_); } 
  //   inline bool operator<( const MinDomainMaxDegree& x ) const { return dom_ < x.dom_ || (dom_ == x.dom_ && x.deg_ < deg_); }
  //   inline bool operator<( const MinDomainMaxDegree& x ) const { return dom_ < x.dom_ || (dom_ == x.dom_ && x.deg_ < deg_); }
  //   inline bool operator<( const MinDomainMaxDegree& x ) const { return dom_ < x.dom_ || (dom_ == x.dom_ && x.deg_ < deg_); }


  //   inline MinDomainMaxDegree& operator*=( const int x ) { dom_ *= x; return *this; }
  //   inline MinDomainMaxDegree& operator+=( const int x ) { dom_ += x * deg_; return *this; }
  //   inline MinDomainMaxDegree& operator-=( const int x ) { dom_ -= x * deg_; return *this; }
  //   inline MinDomainMaxDegree& operator/=( const int x ) { deg_ *= x; return *this; }


  //   inline MinDomainMaxDegree& operator*=( const MinDomainMaxDegree& x ) { dom_ *= x.dom_; deg_ /= x.deg_; return *this; }
  //   inline MinDomainMaxDegree& operator+=( const MinDomainMaxDegree& x ) { dom_ += x.dom_; deg_ -= x.deg_; return *this; }
  //   inline MinDomainMaxDegree& operator-=( const MinDomainMaxDegree& x ) { dom_ -= x.dom_; deg_ += x.deg_; return *this; }
  //   inline MinDomainMaxDegree& operator/=( const MinDomainMaxDegree& x ) { dom_ *= x.dom_; deg_ /= x.deg_; return *this; }

  //   inline void operator=( const MinDomainMaxDegree& x ) { dom_ = x.dom_; deg_ = x.deg_; }
  //   inline void operator=( const Variable x ) { dom_ = x.get_size(); deg_ = x.get_degree(); }
  //   //@}  


  //   std::ostream& display_criterion(std::ostream& os) const {
  //     os << "with minimum domain size (ties broken with max degree)";
  //     return os;
  //   }

  //   std::ostream& display(std::ostream& os) const {
  //     os << dom_ << "+1/" << deg_;
  //     return os;
  //   }
  // };

  // std::ostream& operator<<(std::ostream& os, MinDomainMaxDegree& x);

  // std::ostream& operator<<(std::ostream& os, MinDomainMaxDegree* x);


  /*! \class MinDomainOverDegree
    \brief  Class MinDomainOverDegree

    Order two variables by the ratio of their domain sizes and degree
  */
  class MinDomainOverDegree 
  {
  public: 

    /**@name Constructors*/
    //@{
    MinDomainOverDegree() {dom_ = LARGE_VALUE; deg_ = 0;}
    void initialise(const double* _w) {dom_ = LARGE_VALUE; deg_ = 0;}
    //@}

    /**@name Parameters*/
    //@{ 
    int dom_;
    int deg_;
    //@}  

    /**@name Utils*/
    //@{
    inline double value() { return (deg_ ? (double)dom_/(double)deg_ : (double)INFTY); } 
    inline bool operator<( const MinDomainOverDegree& x ) const { return dom_*x.deg_ < x.dom_*deg_; }
    inline bool operator>( const MinDomainOverDegree& x ) const { return dom_*x.deg_ > x.dom_*deg_; }
    inline bool operator==( const MinDomainOverDegree& x ) const { return dom_*x.deg_ == x.dom_*deg_; }

    // inline MinDomainOverDegree& operator*( const MinDomainOverDegree& x ) const { MinDomainOverDegree d; d.dom_=(dom_ * x.dom_); d.deg_=(deg_ * x.deg_); return d; }
    // inline MinDomainOverDegree& operator+( const MinDomainOverDegree& x ) const { MinDomainOverDegree d; d.dom_=(dom_ + x.dom_); d.deg_=(deg_ + x.deg_); return d; }
    // inline MinDomainOverDegree& operator-( const MinDomainOverDegree& x ) const { MinDomainOverDegree d; d.dom_=(dom_ - x.dom_); d.deg_=(deg_ - x.deg_); return d; }
    // inline MinDomainOverDegree& operator/( const MinDomainOverDegree& x ) const { MinDomainOverDegree d; d.dom_=(dom_ / x.dom_); d.deg_=(deg_ / x.deg_); return d; }

  

    inline MinDomainOverDegree& operator*=( const int x ) { dom_ *= x; return *this; }
    inline MinDomainOverDegree& operator+=( const int x ) { dom_ += x * deg_; return *this; }
    inline MinDomainOverDegree& operator-=( const int x ) { dom_ -= x * deg_; return *this; }
    inline MinDomainOverDegree& operator/=( const int x ) { deg_ *= x; return *this; }


    inline MinDomainOverDegree& operator*=( const MinDomainOverDegree& x ) { dom_ *= x.dom_; deg_ *= x.deg_; return *this; }
    inline MinDomainOverDegree& operator+=( const MinDomainOverDegree& x ) { 
      dom_ = (dom_ * x.deg_ + x.dom_ * deg_); 
      deg_ *= x.deg_; return *this; }
    inline MinDomainOverDegree& operator-=( const MinDomainOverDegree& x ) { 
      dom_ = (dom_ * x.deg_ - x.dom_ * deg_); 
      deg_ *= x.deg_; return *this; }
    inline MinDomainOverDegree& operator/=( const MinDomainOverDegree& x ) { dom_ /= x.dom_; deg_ /= x.deg_; return *this; }

    // inline MinDomainOverDegree& operator*=( const int x ) { dom_ *= x; deg_ *= x; return *this; }
    // inline MinDomainOverDegree& operator+=( const int x ) { dom_ += x; deg_ += x; return *this; }
    // inline MinDomainOverDegree& operator-=( const int x ) { dom_ -= x; deg_ -= x; return *this; }
    // inline MinDomainOverDegree& operator/=( const int x ) { dom_ /= x; deg_ /= x; return *this; }


    // inline MinDomainOverDegree& operator*=( const MinDomainOverDegree& x ) { dom_ *= x.dom_; deg_ *= x.deg_; return *this; }
    // inline MinDomainOverDegree& operator+=( const MinDomainOverDegree& x ) { dom_ += x.dom_; deg_ += x.deg_; return *this; }
    // inline MinDomainOverDegree& operator-=( const MinDomainOverDegree& x ) { dom_ -= x.dom_; deg_ -= x.deg_; return *this; }
    // inline MinDomainOverDegree& operator/=( const MinDomainOverDegree& x ) { dom_ /= x.dom_; deg_ /= x.deg_; return *this; }

    inline void operator=( const MinDomainOverDegree& x ) { dom_ = x.dom_; deg_ = x.deg_; }
    inline void operator=( const Variable x ) { dom_ = x.get_size(); deg_ = x.get_degree(); }
    //@}  


    std::ostream& display_criterion(std::ostream& os) const {
      os << "with minimum (domain size / degree)";
      return os;
    }

    std::ostream& display(std::ostream& os) const {
      os << dom_ << "/" << deg_;
      return os;
    }
  };

  std::ostream& operator<<(std::ostream& os, MinDomainOverDegree& x);

  std::ostream& operator<<(std::ostream& os, MinDomainOverDegree* x);


  /*! \class MinDomainOverWeight
    \brief  Class MinDomainOverWeight

    Order two variables by the ratio of their domain sizes and weight
  */
  class MinDomainOverWeight
  {
  public: 

    /**@name Constructors*/
    //@{
    MinDomainOverWeight() {dom_ = LARGE_VALUE; wei_ = 0;}
    void initialise(double* _w=NULL) {dom_ = LARGE_VALUE; wei_ = 0; weight = _w;}
    //MinDomainOverWeight(void *w) : weight((double*)w) {dom_ = LARGE_VALUE; wei_ = 0;}
    //@}

    /**@name Parameters*/
    //@{ 
    double *weight;
    int dom_;
    double wei_;
    //@}  

    /**@name Utils*/
    //@{
    inline double value() { return (wei_ ? (double)dom_/wei_ : (double)INFTY); } 
    inline bool operator<( const MinDomainOverWeight& x ) const { return dom_*x.wei_ < x.dom_*wei_; }
    inline bool operator>( const MinDomainOverWeight& x ) const { return dom_*x.wei_ > x.dom_*wei_; }
    inline bool operator==( const MinDomainOverWeight& x ) const { return dom_*x.wei_ == x.dom_*wei_; }


    // inline MinDomainOverWeight& operator*( const MinDomainOverWeight& x ) const { MinDomainOverWeight d; d.dom_=(dom_ * x.dom_); d.wei_=(wei_ * x.wei_); return d; }
    // inline MinDomainOverWeight& operator+( const MinDomainOverWeight& x ) const { MinDomainOverWeight d; d.dom_=(dom_ + x.dom_); d.wei_=(wei_ + x.wei_); return d; }
    // inline MinDomainOverWeight& operator-( const MinDomainOverWeight& x ) const { MinDomainOverWeight d; d.dom_=(dom_ - x.dom_); d.wei_=(wei_ - x.wei_); return d; }
    // inline MinDomainOverWeight& operator/( const MinDomainOverWeight& x ) const { MinDomainOverWeight d; d.dom_=(dom_ / x.dom_); d.wei_=(wei_ / x.wei_); return d; }


    inline void reduce() {
      while(dom_>100000 || !(dom_&1)) {
	dom_ /=2;
	wei_ /=2;
      }
    }


    inline MinDomainOverWeight& operator*=( const int x ) { dom_ *= x; return *this; }
    inline MinDomainOverWeight& operator+=( const int x ) { dom_ += x * wei_; return *this; }
    inline MinDomainOverWeight& operator-=( const int x ) { dom_ -= x * wei_; return *this; }
    inline MinDomainOverWeight& operator/=( const int x ) { wei_ *= x; return *this; }


    inline MinDomainOverWeight& operator*=( const MinDomainOverWeight& x ) { dom_ *= x.dom_; wei_ *= x.wei_; return *this; }
    inline MinDomainOverWeight& operator+=( const MinDomainOverWeight& x ) { 
      dom_ = (dom_ * x.wei_ + x.dom_ * wei_); 
      wei_ *= x.wei_; 
      reduce();
      return *this; }
    inline MinDomainOverWeight& operator-=( const MinDomainOverWeight& x ) { 
      dom_ = (dom_ * x.wei_ - x.dom_ * wei_); 
      wei_ *= x.wei_; return *this; }
    inline MinDomainOverWeight& operator/=( const MinDomainOverWeight& x ) { dom_ /= x.dom_; wei_ /= x.wei_; return *this; }


    inline void operator=( const MinDomainOverWeight& x ) { dom_ = x.dom_; wei_ = x.wei_; }
    inline void operator=( const Variable x ) { 
      dom_ = x.get_size(); wei_ = weight[x.id()]; 
      //std::cout << x << ": " << dom_ << "/" << wei_ << std::endl;
    }
    //@}  

    std::ostream& display_criterion(std::ostream& os) const {
      os << "with minimum (domain size / weight)";
      return os;
    }

    std::ostream& display(std::ostream& os) const {
      os << dom_ << "/" << wei_;
      return os;
    }
  };

  std::ostream& operator<<(std::ostream& os, MinDomainOverWeight& x);

  std::ostream& operator<<(std::ostream& os, MinDomainOverWeight* x);


  /*! \class MinDomainTimesWeight
    \brief  Class MinDomainTimesWeight

    Order two variables by the ratio of their domain sizes and weight
  */
  class MinDomainTimesWeight
  {
  public: 

    /**@name Constructors*/
    //@{
    MinDomainTimesWeight() {sco_ = LARGE_VALUE; }
    void initialise(double* _w=NULL) {sco_ = LARGE_VALUE; weight = _w;}
    //MinDomainTimesWeight(void *w) : weight((double*)w) {dom_ = LARGE_VALUE; wei_ = 0;}
    //@}

    /**@name Parameters*/
    //@{ 
    double *weight;
    double sco_;
    //@}  

    /**@name Utils*/
    //@{
    inline double value() { return sco_; } 
    inline bool operator<( const MinDomainTimesWeight& x ) const { return sco_ < x.sco_; }
    inline bool operator>( const MinDomainTimesWeight& x ) const { return sco_ > x.sco_; }
    inline bool operator==( const MinDomainTimesWeight& x ) const { return sco_ == x.sco_; }


    inline MinDomainTimesWeight& operator*=( const int x ) { sco_ *= x; return *this; }
    inline MinDomainTimesWeight& operator+=( const int x ) { sco_ += x; return *this; }
    inline MinDomainTimesWeight& operator-=( const int x ) { sco_ -= x; return *this; }
    inline MinDomainTimesWeight& operator/=( const int x ) { sco_ /= x; return *this; }


    inline MinDomainTimesWeight& operator*=( const MinDomainTimesWeight& x ) { sco_ *= x.sco_; return *this; }
    inline MinDomainTimesWeight& operator+=( const MinDomainTimesWeight& x ) { sco_ += x.sco_; return *this; }
    inline MinDomainTimesWeight& operator-=( const MinDomainTimesWeight& x ) { sco_ -= x.sco_; return *this; }
    inline MinDomainTimesWeight& operator/=( const MinDomainTimesWeight& x ) { sco_ /= x.sco_; return *this; }


    inline void operator=( const MinDomainTimesWeight& x ) { sco_ = x.sco_; }
    inline void operator=( const Variable x ) { sco_ = (double)(x.get_size()) * weight[x.id()]; }
    //@}  

    std::ostream& display_criterion(std::ostream& os) const {
      os << "with minimum (domain size * weight)";
      return os;
    }

    std::ostream& display(std::ostream& os) const {
      os << sco_;
      return os;
    }
  };

  std::ostream& operator<<(std::ostream& os, MinDomainTimesWeight& x);

  std::ostream& operator<<(std::ostream& os, MinDomainTimesWeight* x);


  /*! \class MinNeighborDomainOverNeighborWeight
    \brief  Class MinNeighborDomainOverNeighborWeight

    Order two variables by the ratio of the domain sizes and weight
    of a vector of other variables
  */
  class MinNeighborDomainOverNeighborWeight
  {
  public: 

    /**@name Constructors*/
    //@{
    MinNeighborDomainOverNeighborWeight() {dom_ = LARGE_VALUE; wei_ = 0;}
    void initialise(double* _w=NULL) {dom_ = LARGE_VALUE; wei_ = 0; weight = _w;}
    //@}

    /**@name Parameters*/
    //@{ 
    double *weight;
    Vector< Variable > *map;
    int dom_;
    double wei_;
    //@}  

    /**@name Utils*/
    //@{
    inline double value() { return (wei_ ? (double)dom_/wei_ : (double)INFTY); } 
    inline bool operator<( const MinNeighborDomainOverNeighborWeight& x ) const { return dom_*x.wei_ < x.dom_*wei_; }
    inline bool operator>( const MinNeighborDomainOverNeighborWeight& x ) const { return dom_*x.wei_ > x.dom_*wei_; }
    inline bool operator==( const MinNeighborDomainOverNeighborWeight& x ) const { return dom_*x.wei_ == x.dom_*wei_; }
    inline void operator=( const MinNeighborDomainOverNeighborWeight& x ) { dom_ = x.dom_; wei_ = x.wei_; }
    inline void operator=( const Variable x ) { 
      int idx = x.id();
      int i = map[idx].size;
      Variable y;
      dom_ = 0;
      wei_ = 0;
      while(i--) {
  	y = map[idx][i];
  	dom_ += y.get_size(); 
  	wei_ += weight[y.id()];
      } 

      //std::cout << x << ": " << dom_ << "/" << wei_ << std::endl;

    }
    //@}  

    std::ostream& display_criterion(std::ostream& os) const {
      os << "with minimum (Sum of neighrbors' domain sizes / Sum of neighrbors' weights)";
      return os;
    }

    std::ostream& display(std::ostream& os) const {
      os << dom_ << "/" << wei_;
      return os;
    }
  };

  //std::ostream& operator<<(std::ostream& os, MinNeighborDomainOverNeighborWeight& x);

  //std::ostream& operator<<(std::ostream& os, MinNeighborDomainOverNeighborWeight& x);


  /*! \class MinNeighborDomainOverNeighborWeight
    \brief  Class MinNeighborDomainOverNeighborWeight

    Order two variables by the ratio of the domain sizes and weight
    of a vector of other variables
  */
  template< class VarComparator >
  class SelfPlusAverage
  {
  public: 

    /**@name Constructors*/
    //@{
    SelfPlusAverage() { map = NULL; }
    void initialise(double* _w=NULL) { crit.initialise(_w); map = NULL; weight = _w; }
    //@}

    /**@name Parameters*/
    //@{ 
    double *weight;
    VarComparator crit;
    Vector< Variable > *map;
    //@}  

    /**@name Utils*/
    //@{
    inline double value() { return 1; } 
    inline bool operator<( const SelfPlusAverage<VarComparator> & x ) const { return crit < x.crit; }
    inline bool operator>( const SelfPlusAverage<VarComparator> & x ) const { return crit > x.crit; }
    inline bool operator==( const SelfPlusAverage<VarComparator> & x ) const { return crit == x.crit; }
    inline void operator=( const SelfPlusAverage<VarComparator>& x ) { crit = x.crit; }
    inline void operator=( const Variable x ) { 



      crit.weight = weight;
      crit = x;


      // std::cout << "\n => ";
      // crit.display(std::cout);
      // std::cout << std::endl;
      // std::cout << "scan neighborhood" << std::endl;



      VarComparator crit_neighbor, aux;
      crit_neighbor.weight = weight;
      aux.weight = weight;

      int idx = x.id();

      int n = map[idx].size;

      int i = n-1;
      Variable y = map[idx][i];
      crit_neighbor = y;

      // std::cout << "   ---> ";
      // crit_neighbor.display(std::cout);
      // std::cout << std::endl;

      while(--i>=0) {
  	y = map[idx][i];
	aux = y;
	crit_neighbor += aux;

	// std::cout << "      + ";
	// aux.display(std::cout);
	// std::cout << " = " ;
	// crit_neighbor.display(std::cout);
	// std::cout << std::endl;
      } 
      crit_neighbor /= n;

      // std::cout << "   ---> ";
      //  crit_neighbor.display(std::cout);
      // std::cout << std::endl;
      
      crit += crit_neighbor;


      // std::cout << "   ---> " ;
      //  crit.display(std::cout);
      // std::cout << std::endl;


      //std::cout << x << ": " << dom_ << "/" << wei_ << std::endl;

    }
    //@}  

    std::ostream& display_criterion(std::ostream& os) const {
      os << "with sum of neighbors' " ;
      crit.display_criterion(os);
      os << "plus self" ;
      return os;
    }

    std::ostream& display(std::ostream& os) const {
      //os << crit;
      crit.display(os);
      return os;
    }
  };

  //std::ostream& operator<<(std::ostream& os, MinNeighborDomainOverNeighborWeight& x);

  //std::ostream& operator<<(std::ostream& os, MinNeighborDomainOverNeighborWeight* x);


  /*! \class MinNeighborDomainOverWeight
    \brief  Class MinNeighborDomainOverWeight

    Order two variables by the ratio of the domain sizes and weight
    of a vector of other variables
  */
  class MinNeighborDomainOverWeight
  {
  public: 

    /**@name Constructors*/
    //@{
    MinNeighborDomainOverWeight() {dom_ = LARGE_VALUE; wei_ = 0;}
    void initialise(double* _w=NULL) {dom_ = LARGE_VALUE; wei_ = 0; weight = _w;}
    //@}

    /**@name Parameters*/
    //@{ 
    double *weight;
    Vector< Variable > *map;
    double dom_;
    double wei_;
    //@}  

    /**@name Utils*/
    //@{
    inline double value() { return (wei_ ? (double)dom_/wei_ : (double)INFTY); } 
    inline bool operator<( const MinNeighborDomainOverWeight& x ) const { return dom_*x.wei_ < x.dom_*wei_; }
    inline bool operator>( const MinNeighborDomainOverWeight& x ) const { return dom_*x.wei_ > x.dom_*wei_; }
    inline bool operator==( const MinNeighborDomainOverWeight& x ) const { return dom_*x.wei_ == x.dom_*wei_; }
    inline void operator=( const MinNeighborDomainOverWeight& x ) { dom_ = x.dom_; wei_ = x.wei_; }
    inline void operator=( const Variable x ) { 

      int idx = x.id();
      int i = map[idx].size;
      Variable y;

#ifdef _DEBUG_VARORD
      std::cout << "check " << x << " (" << weight[idx] 
		<< "): " << map[idx][0] << " in " <<  map[idx][0].get_domain() 
		<< " <> " << map[idx][1] << " in " <<  map[idx][1].get_domain() ; //<< std::endl;
#endif

      wei_ = weight[idx];
      dom_ = 0;
      while(i--) {
  	y = map[idx][i];
  	dom_ += y.get_size(); 
      } 
    }
    //@}  

    std::ostream& display_criterion(std::ostream& os) const {
      os << " with minimum (Sum of neighrbors' domain sizes / weight)";
      return os;
    }

    std::ostream& display(std::ostream& os) const {
      os << dom_ << "/" << wei_;
      return os;
    }
  };

  //std::ostream& operator<<(std::ostream& os, MinNeighborDomainOverWeight& x);


  /*! \class MaxWeight
    \brief  Class MaxWeight

    Order two variables by their weights
  */
  class MaxWeight
  {
  public: 

    /**@name Constructors*/
    //@{
    MaxWeight() {wei_ = 0;}
    void initialise(double* _w=NULL) {wei_ = 0; weight = _w;}
    //MaxWeight(int *w) : weight(w) {wei_ = 0;}
    //@}

    /**@name Parameters*/
    //@{ 
    double *weight;
    double wei_;
    //@}  

    /**@name Utils*/
    //@{
    inline double value() { return wei_ ; } 
    inline bool operator<( const MaxWeight& x ) const { return x.wei_ < wei_; }
    inline bool operator>( const MaxWeight& x ) const { return x.wei_ > wei_; }
    inline bool operator==( const MaxWeight& x ) const { return x.wei_ == wei_; }
    inline void operator=( const MaxWeight& x ) { wei_ = x.wei_; }
    inline void operator=( const Variable x ) { wei_ = weight[x.id()]; }
    //@}  


    std::ostream& display_criterion(std::ostream& os) const {
      os << " with maximum weight";
      return os;
    }

    std::ostream& display(std::ostream& os) const {
      os << wei_;
      return os;
    }
  };

  std::ostream& operator<<(std::ostream& os, MaxWeight& x);

  std::ostream& operator<<(std::ostream& os, MaxWeight* x);



  /*! \class MinWeight
    \brief  Class MinWeight

    Order two variables by their weights
  */
  class MinWeight
  {
  public: 

    /**@name Constructors*/
    //@{
    MinWeight() {wei_ = 0;}
    void initialise(double* _w=NULL) {wei_ = 0; weight = _w;}
    //MinWeight(int *w) : weight(w) {wei_ = 0;}
    //@}

    /**@name Parameters*/
    //@{ 
    double *weight;
    double wei_;
    //@}  

    /**@name Utils*/
    //@{
    inline double value() { return wei_ ; } 
    inline bool operator<( const MinWeight& x ) const { return x.wei_ > wei_; }
    inline bool operator>( const MinWeight& x ) const { return x.wei_ < wei_; }
    inline bool operator==( const MinWeight& x ) const { return x.wei_ == wei_; }
    inline void operator=( const MinWeight& x ) { wei_ = x.wei_; }
    inline void operator=( const Variable x ) { wei_ = weight[x.id()]; }
    //@}  


    std::ostream& display_criterion(std::ostream& os) const {
      os << " with maximum weight";
      return os;
    }

    std::ostream& display(std::ostream& os) const {
      os << wei_;
      return os;
    }
  };

  std::ostream& operator<<(std::ostream& os, MinWeight& x);

  std::ostream& operator<<(std::ostream& os, MinWeight* x);


  /**********************************************
   * Branching Decisions
   **********************************************/

  /*! \class AnyValue
    \brief  Class AnyValue

    Assigns the variable to a value in its domain.
  */
  class AnyValue {
    
  public: 
    
    AnyValue() {}
    AnyValue(Solver *s, double **vw, double *bw) {}
    void initialise(Solver *s, double **vw, double *bw) {}
    //AnyValue(Solver *s, void *a) {}
    virtual ~AnyValue() {};
    
    inline Decision make(Variable x) {
      Decision d(x, Decision::ASSIGNMENT, x.get_first());
      return d;
    }

     std::ostream& display(std::ostream& os) const {
       os << "assign it to a value in its domain";
       return os;
     }

  };

  std::ostream& operator<<(std::ostream& os, AnyValue& x);

  std::ostream& operator<<(std::ostream& os, AnyValue* x);


  /*! \class MinValue
    \brief  Class MinValue

    Assigns the variable to its minimum value.
  */
  class MinValue {
    
  public: 
    
    MinValue() {}
    MinValue(Solver *s, double **vw, double *bw) {}
    void initialise(Solver *s, double **vw, double *bw) {}
    //MinValue(Solver *s, void *a) {}
    virtual ~MinValue() {};
    
    inline Decision make(Variable x) {

      //std::cout << "(MV) make a decision on " << x << " in " << x.get_domain() << std::endl;

      Decision d(x, Decision::ASSIGNMENT, x.get_min());

      //std::cout << d << std::endl;

      return d;
    }

     std::ostream& display(std::ostream& os) const {
       os << "assign it to the minimum value in its domain";
       return os;
     }

  };

  std::ostream& operator<<(std::ostream& os, MinValue& x);

  std::ostream& operator<<(std::ostream& os, MinValue* x);


  /*! \class MaxValue
    \brief  Class MaxValue

    Assigns the variable to its maximum value.
  */
  class MaxValue {

  public: 
    
    MaxValue() {}
    MaxValue(Solver *s, double **vw, double *bw) {}
    void initialise(Solver *s, double **vw, double *bw) {}
    //MaxValue(Solver *s, void *a) {}
    virtual ~MaxValue() {};
    
    inline Decision make(Variable x) {
      Decision d(x, Decision::ASSIGNMENT, x.get_max());
      return d;
    }

     std::ostream& display(std::ostream& os) const {
       os << "assign it to the maximum value in its domain";
       return os;
     }

  };

  std::ostream& operator<<(std::ostream& os, MaxValue& x);

  std::ostream& operator<<(std::ostream& os, MaxValue* x);


  /*! \class MiddleValue
    \brief  Class MiddleValue

    Assigns the variable to its maximum value.
  */
  class MiddleValue {

  public: 
    
    MiddleValue() {}
    MiddleValue(Solver *s, double **vw, double *bw) {}
    void initialise(Solver *s, double **vw, double *bw) {}
    //MiddleValue(Solver *s, void *a) {}
    virtual ~MiddleValue() {};
    
    inline Decision make(Variable x) {
      int val = (x.get_min()+x.get_max())/2;
      if(!x.contain(val)) val = x.next(val);
      Decision d(x, Decision::ASSIGNMENT, val);
      return d;
    }

     std::ostream& display(std::ostream& os) const {
       os << "assign it to the mean of its bounds";
       return os;
     }

  };

  std::ostream& operator<<(std::ostream& os, MiddleValue& x);

  std::ostream& operator<<(std::ostream& os, MiddleValue* x);


 /*! \class MedianValue
    \brief  Class MedianValue

    Assigns the variable randomly a value in its domain.
  */
  class MedianValue {

  public: 
    
    MedianValue() {}
    MedianValue(Solver *s, double **vw, double *bw) {}
    void initialise(Solver *s, double **vw, double *bw) {}
    virtual ~MedianValue() {};
    
    inline Decision make(Variable x) {
      int val = x.get_min(); //, rank = randint(x.get_size());
      Decision d(x, Decision::ASSIGNMENT, INFTY);

      if(x.is_range())
	d.set_value((x.get_min()+x.get_max())/2);
      else {
	int prev_val, rank = x.get_size()/2; // = nxt_val-1;
	do {
	  prev_val = val;
	  val = x.next(prev_val);
	} while (rank-- && val > prev_val);
	d.set_value(prev_val);
      }

      return d;
    }

     std::ostream& display(std::ostream& os) const {
       os << "assign it to its median value";
       return os;
     }

  };

  std::ostream& operator<<(std::ostream& os, MedianValue& x);

  std::ostream& operator<<(std::ostream& os, MedianValue* x);


  /*! \class HalfSplit
    \brief  Class HalfSplit

    Set the upper bound of the variable to (ub+lb)/2.
  */
  class HalfSplit {

  public: 
    
    HalfSplit() {}
    HalfSplit(Solver *s, double **vw, double *bw) {}
    void initialise(Solver *s, double **vw, double *bw) {}
    //HalfSplit(Solver *s, void *a) {}
    virtual ~HalfSplit() {};
    
    inline Decision make(Variable x) {

      //std::cout << "(HS) make a decision on " << x << " in " << x.get_domain() << std::endl;

      Decision d(x, Decision::UPPERBOUND, (x.get_min()+x.get_max())/2);
      return d;
    }

     std::ostream& display(std::ostream& os) const {
       os << "set its upper bound to (min+max)/2";
       return os;
     }

  };

  std::ostream& operator<<(std::ostream& os, HalfSplit& x);

  std::ostream& operator<<(std::ostream& os, HalfSplit* x);



  /*! \class ReverseSplit
    \brief  Class ReverseSplit

    Set the upper bound of the variable to (ub+lb)/2.
  */
  class ReverseSplit {

  public: 
    
    ReverseSplit() {}
    ReverseSplit(Solver *s, double **vw, double *bw) {}
    void initialise(Solver *s, double **vw, double *bw) {}
    //ReverseSplit(Solver *s, void *a) {}
    virtual ~ReverseSplit() {};
    
    inline Decision make(Variable x) {
      Decision d(x, Decision::LOWERBOUND, (x.get_min()+x.get_max()+1)/2);
      return d;
    }

     std::ostream& display(std::ostream& os) const {
       os << "set its lower bound to (min+max+1)/2";
       return os;
     }

  };

  std::ostream& operator<<(std::ostream& os, ReverseSplit& x);

  std::ostream& operator<<(std::ostream& os, ReverseSplit* x);


  /*! \class RandomSplit
    \brief  Class RandomSplit

    Set the upper bound of the variable to (ub+lb)/2.
  */
  class RandomSplit {

  public: 
    
    RandomSplit() {}
    RandomSplit(Solver *s, double **vw, double *bw) {}
    void initialise(Solver *s, double **vw, double *bw) {}
    //RandomSplit(Solver *s, void *a) {}
    virtual ~RandomSplit() {};
    
    inline Decision make(Variable x) {
      int the_min = x.get_min();
      Decision d(x, Decision::UPPERBOUND, (the_min+randint(x.get_max()-the_min)));
      return d;
    }

     std::ostream& display(std::ostream& os) const {
       os << "set its upper bound to a random value in [min,max]";
       return os;
     }

  };

  std::ostream& operator<<(std::ostream& os, RandomSplit& x);

  std::ostream& operator<<(std::ostream& os, RandomSplit* x);


  /*! \class RandomMinMax
    \brief  Class RandomMinMax

    Assigns the variable randomly its minimum or maximum value.
  */
  class RandomMinMax {

  public: 
    
    RandomMinMax() {}
    RandomMinMax(Solver *s, double **vw, double *bw) {}
    void initialise(Solver *s, double **vw, double *bw) {}
    //RandomMinMax(Solver *s, void *a) {}
    virtual ~RandomMinMax() {};
    
    inline Decision make(Variable x) {
      Decision d(x, Decision::ASSIGNMENT, 
		 (randint(2) ? x.get_min() : x.get_max()));
      return d;
    }

     std::ostream& display(std::ostream& os) const {
       os << "assign it to its minimum or maximum value at random with equal probability";
       return os;
     }

  };

  std::ostream& operator<<(std::ostream& os, RandomMinMax& x);

  std::ostream& operator<<(std::ostream& os, RandomMinMax* x);


  /*! \class RandomValue
    \brief  Class RandomValue

    Assigns the variable randomly a value in its domain.
  */
  class RandomValue {

  public: 
    
    RandomValue() {}
    RandomValue(Solver *s, double **vw, double *bw) {}
    void initialise(Solver *s, double **vw, double *bw) {}
    virtual ~RandomValue() {};
    
    inline Decision make(Variable x) {
      int val = x.get_min(), rank = randint(x.get_size());
      Decision d(x, Decision::ASSIGNMENT, INFTY);

      if(x.is_range())
	d.set_value(val+rank);
      else {
	int prev_val; // = nxt_val-1;
	do {
	  prev_val = val;
	  val = x.next(prev_val);
	} while (rank-- && val > prev_val);
	d.set_value(prev_val);
      }

      return d;
    }

     std::ostream& display(std::ostream& os) const {
       os << "assign it to a value at random";
       return os;
     }

  };

  std::ostream& operator<<(std::ostream& os, RandomValue& x);

  std::ostream& operator<<(std::ostream& os, RandomValue* x);


  /*! \class MinWeightValue
    \brief  Class MinWeightValue

    Assigns the variable to its value with minimum weight for some weight matrix.
  */
  class MinWeightValue {

  public: 

    double **weight;
    
    MinWeightValue() {}
    MinWeightValue(Solver *s, double **vw, double *bw) {
      initialise(s, vw, bw);
    }
    void initialise(Solver *s, double **vw, double *bw) {
      weight = vw;
    }
    //MinWeightValue(Solver *s, void *a) { weight = (double**)a; }
    virtual ~MinWeightValue() {};
    
    inline Decision make(Variable x) {
      int // id_x = x.id(),
	best_val = x.get_min();
      double *wgt = weight[x.id()];
      //double min_weight = weight[id_x][best_val]// [best_val][id_x]
      double min_weight = wgt[best_val]// [best_val][id_x]
	, aux_weight;
      int vali, vnxt=x.next(best_val);
      do {
	vali = vnxt;
	vnxt = x.next(vali);
	aux_weight = wgt[vali]; //weight[id_x][vali]; //weight[vali][id_x];
	if(aux_weight < min_weight) {
	  min_weight = aux_weight;
	  best_val = vali;
	}
      } while(vali<vnxt);
      
      Decision d(x, Decision::ASSIGNMENT, best_val);
      return d;
    }

     std::ostream& display(std::ostream& os) const {
       os << "assign it to the value with minimum weight in its domain";
       return os;
     }

  };

  std::ostream& operator<<(std::ostream& os, MinWeightValue& x);

  std::ostream& operator<<(std::ostream& os, MinWeightValue* x);


  /*! \class MaxWeightValue
    \brief  Class MaxWeightValue

    Assigns the variable to its value with minimum weight for some weight matrix.
  */
  class MaxWeightValue {

  public: 

    double **weight;
    
    MaxWeightValue() {}
    MaxWeightValue(Solver *s, double **vw, double *bw) {
      initialise(s, vw, bw);
    }
    void initialise(Solver *s, double **vw, double *bw) {
      weight = vw;
    }
    //MaxWeightValue(Solver *s, void *a) { weight = (double**)a; }
    virtual ~MaxWeightValue() {};
    
    inline Decision make(Variable x) {
      int // id_x = x.id(),
	best_val = x.get_min();
      double *wgt = weight[x.id()];
      //double min_weight = weight[id_x][best_val]// [best_val][id_x]
      double max_weight = wgt[best_val]// [best_val][id_x]
	, aux_weight;
      int vali, vnxt=x.next(best_val);
      do {
	vali = vnxt;
	vnxt = x.next(vali);
	aux_weight = wgt[vali]; //weight[id_x][vali]; //weight[vali][id_x];
	if(aux_weight > max_weight) {
	  max_weight = aux_weight;
	  best_val = vali;
	}
      } while(vali<vnxt);
      
      Decision d(x, Decision::ASSIGNMENT, best_val);
      return d;
    }

     std::ostream& display(std::ostream& os) const {
       os << "assign it to the value with minimum weight in its domain";
       return os;
     }

  };

  std::ostream& operator<<(std::ostream& os, MaxWeightValue& x);

  std::ostream& operator<<(std::ostream& os, MaxWeightValue* x);



  /*! \class MinWeightBound
    \brief  Class MinWeightBound

    Assigns the variable to its bound with minimum weight for some weight matrix.
  */
  class MinWeightBound {

  public: 

    double *weight;
    
    MinWeightBound() {}
    MinWeightBound(Solver *s, double **vw, double *bw) {
      initialise(s, vw, bw);
    }
    void initialise(Solver *s, double **vw, double *bw) {
      weight = bw;
    }
    //MinWeightBound(Solver *s, void *a) { weight = (double**)a; }
    virtual ~MinWeightBound() {};
    
    inline Decision make(Variable x) {
      int id_x = 2*x.id();
      int best_val = (weight[id_x+1] < weight[id_x] ? x.get_max() : x.get_min());
      Decision d(x, Decision::ASSIGNMENT, best_val);
      return d;
    }

     std::ostream& display(std::ostream& os) const {
       os << "assign it to the bound with minimum weight";
       return os;
     }

  };

  std::ostream& operator<<(std::ostream& os, MinWeightBound& x);

  std::ostream& operator<<(std::ostream& os, MinWeightBound* x);


  // /*! \class Guided
  //   \brief  Class Guided

  //   Assigns the variable to its minimum value.
  // */
  // class Guided {
    
  // public: 
    
  //   Solver *solver;
    
  //   Guided() {}
  //   Guided(Solver *s) {solver=s;}
  //   virtual ~Guided() {};
    
  //   inline Decision make(Variable x) {
  //     int val = solver->last_solution_lb[x.id()];
  //     Decision d(x, Decision::ASSIGNMENT, val);
  //     if(val == -INFTY || !x.contain(val)) 
  // 	//d.set_value((randint(2) ? x.get_min() : x.get_max()));
  // 	d.set_value(x.get_min());
  //     return d;
  //   }


  //    std::ostream& display(std::ostream& os) const {
  //      os << "assign it to the value of this variable in the last solution";
  //      return os;
  //    }

  // };

  // std::ostream& operator<<(std::ostream& os, Guided& x);

  // std::ostream& operator<<(std::ostream& os, Guided* x);


  /*! \class Guided
    \brief  Class Guided

    Assigns the variable to its minimum value.
  */
  template< class Default >
  class Guided {
    
  public: 
    
    Solver *solver;
    Default init_choice;
    
    int cool;

    Guided() {cool=150;}
    Guided(Solver *s, double **vw, double *bw) { initialise(s, vw, bw); cool=150;}
    void initialise(Solver *s, double **vw, double *bw) { solver=s, init_choice.initialise(s, vw, bw); }
    virtual ~Guided() {};
    
    inline Decision make(Variable x) {


      //std::cout << "(G) make a decision on " << x << " in " << x.get_domain() << std::endl;
 
      Decision d;
      
      int val = solver->last_solution_lb[x.id()];
      if(val != -INFTY && x.contain(val)) {
	d = Decision(x, Decision::ASSIGNMENT, val);
      } else
	d = init_choice.make(x);
 
      //std::cout << d << std::endl;

      return d;

    }
    
    std::ostream& display(std::ostream& os) const {
      os << "assign it to the value of this variable in the last solution";
      return os;
    }
    
  };

  template< class Default >
  std::ostream& operator<<(std::ostream& os, Guided<Default>& x) {
    return x.display(os);
  }


  template< class Default >
  std::ostream& operator<<(std::ostream& os, Guided<Default>* x) {
    return x->display(os);
  }




  /*! \class GuidedSplit
    \brief  Class GuidedSplit

    Restricts the variable to the half that contains the solution value
  */
  template< class Default >
  class GuidedSplit {
    
  public: 
    
    Solver *solver;
    Default init_choice;
    
    GuidedSplit() {}
    GuidedSplit(Solver *s, double **vw, double *bw) { initialise(s, vw, bw); }
    void initialise(Solver *s, double **vw, double *bw) { solver=s; init_choice.initialise(s, vw, bw); }
    virtual ~GuidedSplit() {};
    
    inline Decision make(Variable x) {

      //std::cout << "(GS) make a decision on " << x << " in " << x.get_domain() << std::endl;

      int val = solver->last_solution_lb[x.id()];
      Decision d;
      if(val == -INFTY || !x.contain(val)) 
	d = init_choice.make(x);
      else {
	int half = (x.get_min()+x.get_max())/2;
	if(half < val)
	  d = Decision(x, Decision::LOWERBOUND, half);
	else
	  d = Decision(x, Decision::UPPERBOUND, half);
      }
      return d;
    }
    
    std::ostream& display(std::ostream& os) const {
      os << "halves the domain so that it keeps the value of this variable in the last solution";
      return os;
    }
    
  };


  template< class Default >
  std::ostream& operator<<(std::ostream& os, GuidedSplit<Default>& x) {
    return x.display(os);
  }

  template< class Default >
  std::ostream& operator<<(std::ostream& os, GuidedSplit<Default>* x) {
    return x->display(os);
  }




  /*! \class ConditionalOnSize
    \brief  Class ConditionalOnSize

    - uses range_branching if the domain size is continuous and larger than threshold 
    - uses fd_branching otherwise
  */
  template< class RangeBranching, class FDBranching >
  class ConditionalOnSize {
    
  public: 
    
    Solver *solver;
    RangeBranching range_branching;
    FDBranching fd_branching;

    int threshold;
    
    ConditionalOnSize() {}
    ConditionalOnSize(Solver *s, double **vw, double *bw) { initialise(s, vw, bw); }
    void initialise(Solver *s, double **vw, double *bw) {

      //std::cout << "initialise COS" << std::endl;
      //std::cout << (int*)s << std::endl;

      solver=s; 

      //std::cout << (int*)solver << std::endl;
      //std::cout << solver << std::endl;


      threshold=10; 
      range_branching.initialise(solver, vw, bw); 
      fd_branching.initialise(solver, vw, bw);  
    }
    virtual ~ConditionalOnSize() {};
    
    inline Decision make(Variable x) {


      //std::cout << "(COS) make a decision on " << x << " in " << x.get_domain() << std::endl;

      Decision d;
      if(x.is_range() && x.get_size()>=threshold) {

	//std::cout << "  -> RANGE!\n";

	d = range_branching.make(x);
      } else {

	//std::cout << "  -> FD!\n";
	

	d = fd_branching.make(x);
      }


      //std::cout << "  ==> " << d << std::endl;

      return d;
    }
    
    std::ostream& display(std::ostream& os) const {
      os << "uses ";
      range_branching.display(os);
      os << " on large intervals and ";
      fd_branching.display(os); 
      os << " otherwise";
      return os;
    }
    
  };


  template< class RangeBranching, class FDBranching >
  std::ostream& operator<<(std::ostream& os, ConditionalOnSize< RangeBranching, FDBranching >& x) {
    return x.display(os);
  }

  template< class RangeBranching, class FDBranching >
  std::ostream& operator<<(std::ostream& os, ConditionalOnSize< RangeBranching, FDBranching >* x) {
    return x->display(os);
  }




#define NEG(a) ((2*a))
#define POS(a) ((2*a+1))
  /*! \class BoolMinWeightValue
    \brief  Class BoolMinWeightValue

    Assigns the variable to its value with minimum weight for some weight matrix.
    Assumes that the variable is Boolean
  */
  class BoolMinWeightValue {

  public: 

    //double **weight;
    double *weight;
    
    BoolMinWeightValue() {}
    BoolMinWeightValue(Solver *s, double **vw, double *bw) {
      initialise(s, vw, bw);
    }
    void initialise(Solver *s, double **vw, double *bw) {
      weight = bw;
    }
    // BoolMinWeightValue(Solver *s, void *a) {
    //   weight = (double*)a;
    // }
    virtual ~BoolMinWeightValue() {};
    
    inline Decision make(Variable x) {

      // this thing is tricky: the value 'lit_activity' of a literal gets incremented when we post a new clause 
      // involving this literal. Now, this means that the opposite literal is found to be more constrained, 
      // hence when weighting clause one should use the opposite literal, and when branching one should use the literal
      // with highest activity

      Atom a = x.id();
      Decision d(x, Decision::ASSIGNMENT, (weight[NEG(a)]>weight[POS(a)]));
      return d;
    }

     std::ostream& display(std::ostream& os) const {
       os << "assign it to the value with minimum weight in its (Boolean) domain";
       return os;
     }
        
  };
  
  std::ostream& operator<<(std::ostream& os, BoolMinWeightValue& x);

  std::ostream& operator<<(std::ostream& os, BoolMinWeightValue* x);



  /*! \class BoolMaxWeightValue
    \brief  Class BoolMaxWeightValue

    Assigns the variable to its value with minimum weight for some weight matrix.
    Assumes that the variable is Boolean
  */
  class BoolMaxWeightValue {

  public: 

    //double **weight;
    double *weight;
    
    BoolMaxWeightValue() {}
    BoolMaxWeightValue(Solver *s, double **vw, double *bw) {
      initialise(s, vw, bw);
    }
    void initialise(Solver *s, double **vw, double *bw) {
      weight = bw;
    }
    // BoolMaxWeightValue(Solver *s, void *a) {
    //   weight = (double*)a;
    // }
    virtual ~BoolMaxWeightValue() {};
    
    inline Decision make(Variable x) {

      // this thing is tricky: the value 'lit_activity' of a literal gets incremented when we post a new clause 
      // involving this literal. Now, this means that the opposite literal is found to be more constrained, 
      // hence when weighting clause one should use the opposite literal, and when branching one should use the literal
      // with highest activity

      Atom a = x.id();
      Decision d(x, Decision::ASSIGNMENT, (weight[NEG(a)]<weight[POS(a)]));
      return d;
    }

     std::ostream& display(std::ostream& os) const {
       os << "assign it to the value with maximum weight in its (Boolean) domain";
       return os;
     }
        
  };
  
  std::ostream& operator<<(std::ostream& os, BoolMaxWeightValue& x);

  std::ostream& operator<<(std::ostream& os, BoolMaxWeightValue* x);




  template< class VarComparator >
  std::ostream& operator<<(std::ostream& os, Identifiable<VarComparator>& x) {
    x.criterion.display(os);
    return os;
  }

  
  template < int R = 1 >
  class VSIDS : public GenericDVO< MaxWeight, R, LearningActivityManager > {
  public:
    VSIDS() : GenericDVO< MaxWeight, R, LearningActivityManager >() {}
    VSIDS(Solver *s) : GenericDVO< MaxWeight, R, LearningActivityManager >(s) {}
  };

  template < int R = 1 >
  class WDEG : public GenericDVO< MaxWeight, R, FailureCountManager > {
  public:
    WDEG() : GenericDVO< MaxWeight, R, FailureCountManager >() {}
    WDEG(Solver *s) : GenericDVO< MaxWeight, R, FailureCountManager >(s) {}
  };

  template < int R = 1 >
  class DWDEG : public GenericDVO< MinDomainOverWeight, R, FailureCountManager > {
  public:
    DWDEG() : GenericDVO< MinDomainOverWeight, R, FailureCountManager >() {}
    DWDEG(Solver *s) : GenericDVO< MinDomainOverWeight, R, FailureCountManager >(s) {}
  };

  template < int R = 1 >
  class ABS : public GenericDVO< MinDomainOverWeight, R, PruningCountManager > {
  public:
    ABS() : GenericDVO< MinDomainOverWeight, R, PruningCountManager >() {}
    ABS(Solver *s) : GenericDVO< MinDomainOverWeight, R, PruningCountManager >(s) {}
  };

  template < int R = 1 >
  class Impact : public GenericDVO< MinWeight, R, ImpactManager > {
  public:
    Impact() : GenericDVO< MinWeight, R, ImpactManager >() {}
    Impact(Solver *s) : GenericDVO< MinWeight, R, ImpactManager >(s) {}
  };

  template < int R = 1 >
  class IBS : public GenericDVO< MinDomainTimesWeight, R, ImpactManager > {
  public:
    IBS() : GenericDVO< MinDomainTimesWeight, R, ImpactManager >() {}
    IBS(Solver *s) : GenericDVO< MinDomainTimesWeight, R, ImpactManager >(s) {}
  };


  // typedef GenericDVO< MinDomainOverWeight, 1, FailureCountManager > WDEG;
  
  // typedef GenericDVO< MinDomainOverWeight, 1, PruningCountManager > ABS;
  
}

#endif // __SEARCH_HPP
