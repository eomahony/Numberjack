
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


/*! \file mistral_solver.hpp
    \brief Header file for the solver.
*/


#ifndef __SOLVER_HPP
#define __SOLVER_HPP


#include <mistral_constraint.hpp>

#include <tclap/CmdLine.h>



namespace Mistral {


  class Solution {
    
  public: 
    
    Vector< Variable > variables;
    int min_id;
    int max_id;
    int *values;

    Solution( Vector< Variable >& vars );
    virtual ~Solution();

    virtual std::ostream& display(std::ostream&) const ;
    int& operator[](Variable x) ;

  };


  class SolverParameters {

  public:

    SolverParameters();
    SolverParameters(const SolverParameters&);
    virtual ~SolverParameters();
    void initialise();
    void copy(const SolverParameters&);

    /// random seed
    int seed;

    /// degree of verbosity
    int verbosity;
    /// Number of solutions to find, if equal to -1, then all solutions are listed
    int find_all;

    /// Limit on the number of nodes
    unsigned int node_limit;
    /// Limit on the number of backtracks
    unsigned int backtrack_limit;
    /// Limit on the number of failures
    unsigned int fail_limit;
    /// Limit on the number of failures
    unsigned int propagation_limit;
    /// fail limit used for restarts
    unsigned int restart_limit;
    /// flag to know if there is a limit
    unsigned int limit;

    /// Limit on cpu time
    double time_limit;  

    /// restart policy settings 
    int restart_policy;
    /// restart policy settings (base limit)
    unsigned int restart_base;
    /// restart policy settings (geometric increment)
    double restart_factor;

    /// type of randomization
    unsigned int randomization;
    /// variables sequence shuffle between restarts
    bool shuffle;


    int backjump;

    /// whether solutions are checked
    // 0 -> not checked
    // 1 -> check constraints which variables are all assigned
    // 2 -> 1+check that there exist a support for other constraints
    // 3 -> 2+check that all values are consistent for other constraints
    // NOT IMPLEMENTED! 4 -> check that the solution can be extended to all variables
    int checked; 
     


    /////// PARAMETERS FOR SAT SOLVING ///////
    double            activity_increment;
    double            normalize_activity;
    int               init_activity;
    double            forgetfulness;
    double            activity_decay;

    int               value_selection;
    int               dynamic_value;


    /// MISC
    std::string prefix_comment;
    std::string prefix_statistics;
    std::string prefix_objective;
    std::string prefix_solution;
    std::string prefix_outcome;


  };

  class SolverStatistics {
  
  public:

    SolverStatistics(Solver *s=NULL);
    SolverStatistics(const SolverStatistics&);
    virtual ~SolverStatistics();
    void initialise(Solver *s);
    void copy(const SolverStatistics&);
    void update(const SolverStatistics&);


    Solver *solver;

    /// Number of nodes, that is recursive calls to  the dfs algo
    unsigned long int num_nodes; 
    /// Number of backtracks, that is unsuccesful recursive calls 
    unsigned long int num_backtracks;
    /// Number of constraint failures
    unsigned long int num_failures; 
    /// Number of constraint failures
    unsigned long int num_restarts; 
    /// Number of calls to a constraint propagator
    unsigned long int num_propagations;
    /// Number of solutions found so far
    unsigned long int num_solutions;
    /// Number of inference steps
    unsigned long int num_filterings;
    /// Search outcome
    Outcome outcome;
    /// Objective value (ub for minimization, lb for maximization, -1 otherwise)
    long int objective_value;

    /// maximum length of a partial solution
    unsigned long int max_depth;

    /// Number of nodes, that is recursive calls to  the dfs algo
    unsigned long int num_variables; 
    /// Number of nodes, that is recursive calls to  the dfs algo
    unsigned long int num_values; 
    /// Number of nodes, that is recursive calls to  the dfs algo
    unsigned long int num_constraints; 
    ///
    unsigned long int num_clauses; 
    ///
    unsigned long int max_arity;
    ///
    unsigned long int num_learned;
    ///
    unsigned long int size_learned;
    ///
    double avg_learned_size;
    ///
    unsigned long int avg_amsc_expl_size;
    /// 
    unsigned long int num_amsc_explanations;

    /// 
    unsigned long int num_branch_on_large_domains;
    
    /// timestamp
    double creation_time;
    /// timestamp
    double start_time;
    /// timestamp
    double end_time;

    //unsigned int memory;

    // /////// STATISTICS FOR SAT SOLVING ///////
    // unsigned int      literals;
    // unsigned int      small;
    // double            base_avg_size;
    // double            learnt_avg_size;

    // MISC
    bool pseudo_boolean;


#ifdef _PROFILING

    void init_prof();

    double total_propag_time;
    //std::vector< double > constraint_propag_time;
    
    double total_branching_time;

    double total_restore_time;

    //int vartype_index[17];

#ifdef _PROFILING_PRIMITIVE

    double prof_time[NUM_METHODS][NUM_VARTYPES];
    unsigned long long int prof_num[NUM_METHODS][NUM_VARTYPES];
    virtual std::ostream& print_profile(std::ostream&) const ;

#endif

#endif



    //virtual std::string getString() const;
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::ostream& print_full(std::ostream&) const ;
    virtual std::ostream& print_short(std::ostream&) const ;
  };
  

  class Solver;
  class SearchMonitor {
    
  public:
    
    SearchMonitor(Solver* s) { solver = s; }
    virtual ~SearchMonitor() {}

    Solver *solver;

    Vector<int> sequence;
    Vector<Variable> sequence_var;
    Vector<Constraint> sequence_con;
    Vector<int> sequence_type;

    std::vector<const char*> strs;

    void add(Variable);
    void add(Constraint);
    void add(const char*);

    virtual std::ostream& display( std::ostream& os ) const;

  };


  /**********************************************
   * MultiQueue
   **********************************************/
  
  /*! \class MultiQueue
    \brief MultiQueue Class
  */
  class ConstraintQueue // : public Vector<Constraint*>
  {

  private:
    Solver *solver;
    int min_priority;
    /// The constraint being processed
    //ConstraintImplementation* taboo_constraint;

  public:

    int      cardinality;
    Queue      *triggers;
    int  higher_priority;

    BitSet _set_;
    
    ConstraintQueue();
    void initialise(Solver *s);
    void initialise(const int min_p, const int max_p, const int size);
    void declare(Constraint c, Solver *s);
    virtual ~ConstraintQueue();
    inline bool empty() { return higher_priority<min_priority; }
    
    void add( Constraint cons );
    void add( ConstraintImplementation *cons);
    void trigger( BinaryConstraint *cons);
    void trigger(TernaryConstraint *cons);
    void trigger( GlobalConstraint *cons);
    void trigger( GlobalConstraint *cons, const int var, const Event evt);
    /*
    void trigger(Constraint cons);
    void trigger(Constraint cons, const Event evt);
    */

    void reset_higher_priority();
    // inline void reset_higher_priority() {
    //   while(--higher_priority>=min_priority && triggers[higher_priority].empty());
    // }

    Constraint select(Vector<Constraint>& constraints);
    // inline Constraint select(Vector<Constraint>& constraints) {
    //   int cons_id = triggers[higher_priority].pop();
    //   Constraint cons = constraints[cons_id];
    //   _set_.fast_remove(cons_id);
    //   if(triggers[higher_priority].empty()) reset_higher_priority();
    //   //taboo_constraint = cons.freeze();
    //   return cons;
    // }
    // // inline void select(Constraint cons) {
    // //   //_set_.remove(cons->id);
    // //   taboo_constraint = cons.freeze();
    // // }

    void clear();
    // inline void clear() {
    //   while(higher_priority>=min_priority) triggers[higher_priority--].clear();
    //   _set_.clear();
    //   //taboo_constraint = NULL;
    // }
    virtual std::ostream& display(std::ostream&) ;
  };



   /**********************************************
   * Solver
   **********************************************/

  /*! \class Solver
    \brief Solver Class

    The solver class do the following:
    * it keeps a list of the constraints and variables 
    * it controls 
    - the propagation queue
    - the bracktracking structures
    - the search stack
  */
  //typedef Varstack< ConstraintElement > ConstraintStack;

  class ConsolidateListener;
  class SolutionListener;
  class RestartListener;
  class SuccessListener;
  class BacktrackListener;
  class DecisionListener;
  class VariableListener;
  class ConstraintListener;
  class BranchingHeuristic;
  class RestartPolicy;
  class Reversible;
  class Expression;
  class Decision;
  class Variable;
  class VarArray;
  class Goal;
  class VariableImplementation;
  class ConstraintClauseBase;
  class Solver : public Environment {

  public:
    /*!@name Parameters*/
    //@{
    /// The current level in the search tree (number of left branches from the root)
    /*! The current level in the search tree (number of left branches from the root)
      level -1: state of the domains as stated when modelling
      level 0: level at wich constraints are added (+ initial propagation)
      level k: after k decisions
    */
    //int level;
    bool     search_started;
    int      backtrack_level;
    int      search_root;
    //Decision deduction;

    /// The set of variables, in the initial order, that is as loaded from the model
    Vector< Variable >   variables;
    /// The set of variables, as originally declared (i.e., as expressions)
    //Vector< Variable >   declared_variables;
    Vector< Variable >   removed_variables;
    Vector< int >     domain_types;
    Vector< int > assignment_level;
    Vector< int > assignment_order;
    ReversibleNum<int> assignment_rank;
    Vector< int > assigned;

#ifdef _MONITOR
    SearchMonitor monitor_list;

    Vector< unsigned int > monitored;
    Vector< unsigned int > monitored_index;
    void monitor(Vector<Variable>& X);
    void monitor(Variable X);
#endif

#ifdef _CHECK_NOGOOD
    Vector< Explanation* >       nogood_origin;
    Vector < Vector< Literal > >  nogood_clause;
    Vector< int >  solution;
    Vector< int >  node_num;
    Vector< Atom >  atom;

    void store_reason(Explanation *expl, Atom a);
    void store_nogood(Vector< Literal >& lc);
    void read_solution(const char* fname);

    void check_nogoods();
#endif



    /// The set of constraints, with special accessors for triggers
    Vector< Constraint >         constraints;
    IntStack              posted_constraints;
    ConstraintQueue       active_constraints;
    //VariableQueue           active_variables;
    
    
    Vector< ConstraintTriggerArray > constraint_graph;


    /// For each level, the list of reversible objects that changed at this level, 
    // /// and will need to be restored
    // //Vector< Variable > saved_vars;
    // Vector< int > saved_vars;
    // /// For each level, the list of reversible objects that changed at this level, 
    // /// and will need to be restored
    // Vector< Constraint > saved_post;
    // Vector< Constraint > saved_relax;

    /// Stores the last solution
    Vector< int > last_solution_lb;
    Vector< int > last_solution_ub;

    /// Set of parameters for the solver
    SolverParameters parameters;
    /// Set of statistics collected by the solver
    SolverStatistics statistics;
    //@}

    /// The search part
    /// These are the search variables.
    VarStack < Variable, ReversibleNum<int> >   sequence;
    //ReversibleNum<int> sequence_size;

    Vector< Decision > decisions;
    //Vector< Clause* > reason;
    Vector< Explanation* > reason_for;
    // // stores the index of the variable in the guilty constraint
    // Vector< int > reason_index;
    //Vector< DomainExplanation* > reason_for;
    Vector< Literal > learnt_clause;
    Vector< Literal > visited_literals;
    BitSet visited;
    int num_search_variables;

    ConstraintClauseBase *base;
    // Vector< double > lit_activity;
    // Vector< double > var_activity;
    //double * lit_activity;
    //double * var_activity;

    /// Variable selection and branching
    BranchingHeuristic *heuristic;
    /// Restart policy
    RestartPolicy *policy; 
    /// Goal 
    Goal *objective;

    Vector<SolutionListener*>     solution_triggers;
    Vector<RestartListener*>       restart_triggers;
    Vector<SuccessListener*>       success_triggers;
    Vector<BacktrackListener*>   backtrack_triggers;
    Vector<DecisionListener*>     decision_triggers;
    Vector<VariableListener*>     variable_triggers;
    Vector<ConstraintListener*> constraint_triggers;

    /// the constraint responsible for the last fail (NULL otherwise) 
    Constraint culprit;
    /// the constraint-index of the last wiped_out variable in the culprit constraint (-1 otherwise)
    int wiped_idx;
    /// keep the index of the variable that produced the last conflict
    int prev_wiped_idx;
    /// the overall-index of the variable responsible for the last trigger before a failure
    int wiper_idx;

    Vector< Expression* > expression_store;

    ConsolidateListener *consolidate_manager;

    class MemoryManager : public Vector<int> {

    private:
      
      Vector<int> allocation;

    public:

      MemoryManager() {
	initialise(4096);
	allocation.initialise(8);
      }

      int reserve(const int n, int*& beg, int*& end) {
	if(size+n > capacity) {
	  extendStack(size+n-capacity);
	}
	allocation.add(size);
	beg = stack_ + size;
	end = stack_ + size + n;

	// return the id of the allocated space. 
	return allocation.size;
      }

      void release(const int id) {
	if((unsigned int)id == allocation.size) {
	  //do {
	  size = allocation.pop();
	  //} while(size < 0);
	} else {
	  allocation[id-1] = -1;
	}
      }

    };

    MemoryManager iterator_space;

    Outcome satisfied();
    Outcome exhausted();
    Outcome interrupted();

    /*!@name Constructor*/
    //@{
    Solver();

    void parse_dimacs(const char* filename);
    void parse_pbo(const char* filename);

    void set_parameters(SolverParameters& p);

    class BooleanMemoryManager {
    public:
      Vector< unsigned int > size;
      Vector<int*> slots;

      BooleanMemoryManager() {
	size.add(0);
	int *nslot = new int[1024];
	std::fill(nslot, nslot+1024, 3);
	slots.add(nslot);
      }
      virtual ~BooleanMemoryManager() {
	while(!slots.empty()) 
	  delete [] slots.pop();
      }

      void add(Variable *x);
      void add(Vector< Variable >& bool_vars);
    };

    unsigned int initialised_vars;
    unsigned int initialised_cons;
    BooleanMemoryManager booleans;
    //@}


    // Vector< Variable > non_boolean_variables;
    // Vector< Constraint > non_explained_constraints;
    // //statistics.pseudo_boolean &= (x.is_boolean() || (objective && objective->objective == x));

    /*!@name Destructor*/
    //@{
    virtual ~Solver();
    bool is_initialised() {return variables.size <= initialised_vars;}
    //@}
  
    /*!@name Declarations*/
    //@{
    int declare(Variable x);
    /// add a variable (prior to search!!!)
    void add(Variable x);
    void remove(Variable x);
    void add(VarArray& x);
    void add(Constraint x); 
    void add(Vector< Literal >& clause); 
    //void add(ConstraintW x); 
    //void add(BranchingHeuristic* h);


    void minimize(Variable X);
    void maximize(Variable X);

    void add(SolutionListener* l);
    void add(RestartListener* l);
    void add(SuccessListener* l);
    void add(BacktrackListener* l);
    void add(DecisionListener* l);
    void add(VariableListener* l);
    void add(ConstraintListener* l);
    
    void remove(SolutionListener* l);
    void remove(RestartListener* l);
    void remove(SuccessListener* l);
    void remove(BacktrackListener* l);
    void remove(DecisionListener* l);
    void remove(VariableListener* l);
    void remove(ConstraintListener* l);

    void forbid(const int var, const int m) {
      domain_types[var] &= (~m);
    }
    void mark_non_convex(const int var) { 
      domain_types[var] &= (~RANGE_VAR); 
    }
    //void add(Vector<Literal>& clause);
    //@}

    /*!@name Trail accessors*/
    //@{
    // inline void save() {
    //   trail_.add(sequence.size);
    //   trail_.add(saved_objs.size);
    //   trail_.add(saved_cons.size);
    //   trail_.add(saved_vars.size);
    //   trail_.add(saved_post.size);
    //   trail_.add(saved_relax.size);

    //   ++statistics.num_nodes;
    //   ++level;
    // }
    /// called by reversible objects when they want to be restored when backtrack to this level
    // inline void save(Reversible* object) { saved_objs.add(object); }
    // inline void save(Constraint c) { saved_cons.add(c); }
    // //void save(Variable x); 
    // void save(const int idx); 
    // //void save(VariableImplementation *x, int dtype); 
    //@}

    /*!@name Propagation accessors*/
    //@{
    void notify_backtrack();
    void notify_success();
    void notify_decision();
    void notify_restart(const double prog);
    void notify_relax(Constraint c);
    void notify_post(Constraint c);
    void notify_add_constraint(Constraint c);
    void notify_change_variable(const int idx);
    void notify_add_variable();
    /// called when the var'th variable of constraint cons changes (with event type evt)
    //void trigger_event(const int var, const Event evt);

    /// achieve propagation closure
    PropagationOutcome propagate(Constraint c, const bool ft=true, const bool ts=false);
    PropagationOutcome checker_propagate(Constraint c, const bool ft=true, const bool ts=false);
    PropagationOutcome bound_checker_propagate(Constraint c, const bool ft=true, const bool ts=false);
    bool propagate(); 
    
    void fail();

    void parity_processing(const int k=1);
    bool simple_rewrite();
    bool rewrite(); 
    bool is_pseudo_boolean() const; 
    void consolidate(); 
    void make_non_convex(const int idx);
    //void specialise(); 
    //void initialise_domains(); 
    //@}

    /*!@name Search accesors*/
    //@{
    /// make a branching decision
    void branch_left();
    /// undo the previous branching decision and enforce the complementary constraint
    Outcome branch_right();
    void backjump();

    /// backtrack one level
    void restore();
    /// backtrack to a given level
    void restore(const int);
    

    //Solution *get_solution(Vector< Variable >& seq);
    void initialise_search(Vector< Variable >& seq, 
			   BranchingHeuristic *heu=NULL, 
			   RestartPolicy *pol=NULL,
			   Goal *goal=NULL);
    
    void initialise_search(VarStack < Variable, ReversibleNum<int> >& seq, 
			   BranchingHeuristic *heu=NULL, 
			   RestartPolicy *pol=NULL,
			   Goal *goal=NULL);

    Outcome sequence_search(Vector< Vector< Variable > >& sequences,
			    Vector< BranchingHeuristic * >& heuristics,
			    Vector< RestartPolicy * >& policies,
			    Vector< Goal * >& goals
			    );

    Outcome get_next_solution();

    /*!
      Launches a depth first search on the sequence 'seq' of variables
      with variable ordering 'heu' and restart policy 'pol'.
      Returns an outcome (SAT/UNSAT/OPT/UNKNOWN) and undo all decisions before returning.
    */
    Outcome depth_first_search(Vector< Variable >& seq, 
			       BranchingHeuristic *heu=NULL, 
			       RestartPolicy *pol=NULL,
			       Goal *goal=NULL,
			       bool _restore_=true);
    /*!
      Launches a depth first search on all variables
      with variable ordering 'heu' and restart policy 'pol'.
      Returns an outcome (SAT/UNSAT/OPT/UNKNOWN) and undo all decisions before returning.
    */
    Outcome depth_first_search(BranchingHeuristic *heu=NULL, 
			       RestartPolicy *pol=NULL,
			       Goal *goal=NULL,
			       bool _restore_=true);

    Outcome restart_search(const int root=0, const bool _restore_=true);

    /*!
      Black box search.
    */
    Outcome solve();
    Outcome search_minimize(Variable X);
    Outcome search_maximize(Variable X);

    ///
    bool limits_expired();

    /// depth first search algorithm
    Outcome chronological_dfs(const int root=0);

    // /// sat search algorithm
    Outcome conflict_directed_backjump();
    void learn_nogood();
    void forget();
    double get_current_target();
    // //@}


    BranchingHeuristic *heuristic_factory(std::string var_ordering, std::string branching, const int randomness=1);
    // {
    //   BranchingHeuristic *heu = NULL;
    //   if(var_ordering == "dom/wdeg") {
    // 	if(branching == "minval") {
    // 	  heu = new GenericHeuristic < GenericWeightedDVO < FailureCountManager, MinDomainOverWeight >, MinValue > (this); 
    // 	} else if(branching == "maxval") {
    // 	  heu = new GenericHeuristic < GenericWeightedDVO < FailureCountManager, MinDomainOverWeight >, MaxValue > (this); 
    // 	} else if(branching == "halfsplit") {
    // 	  heu = new GenericHeuristic < GenericWeightedDVO < FailureCountManager, MinDomainOverWeight >, HalfSplit > (this); 
    // 	} else if(branching == "randminmax") {
    // 	  heu = new GenericHeuristic < GenericWeightedDVO < FailureCountManager, MinDomainOverWeight >, RandMinMax > (this); 
    // 	} else if(branching == "minweight") {
    // 	  heu = new GenericHeuristic < GenericWeightedDVO < FailureCountManager, MinDomainOverWeight >, MinWeightValue > (this); 
    // 	} 
    //   } else if(var_ordering == "dom/activity") {
    // 	if(branching == "minval") {
    // 	  heu = new GenericHeuristic < GenericWeightedDVO < PruningCountManager, MinDomainOverWeight >, MinValue > (this); 
    // 	} else if(branching == "maxval") {
    // 	  heu = new GenericHeuristic < GenericWeightedDVO < PruningCountManager, MinDomainOverWeight >, MaxValue > (this); 
    // 	} else if(branching == "halfsplit") {
    // 	  heu = new GenericHeuristic < GenericWeightedDVO < PruningCountManager, MinDomainOverWeight >, HalfSplit > (this); 
    // 	} else if(branching == "randminmax") {
    // 	  heu = new GenericHeuristic < GenericWeightedDVO < PruningCountManager, MinDomainOverWeight >, RandMinMax > (this); 
    // 	} else if(branching == "minweight") {
    // 	  heu = new GenericHeuristic < GenericWeightedDVO < PruningCountManager, MinDomainOverWeight >, MinWeightValue > (this); 
    // 	} 
    //   }						
    //   return heu;
    // }

    RestartPolicy *restart_factory(std::string rpolicy);
    // {
    //   RestartPolicy pol = new NoRestart();
    //   if(policy == "luby") pol = new Luby(); 
    //   else if(policy == "geom") pol = new Geometric(); 
    //   return pol;
    // }

    double *get_literal_activity();

    //void check_constraints();

    void check_constraint_graph_integrity();

    void extract_instance_statistics();

    void initialise_random_seed(const int seed);
    void set_time_limit(const double limit);
    void set_learning_on();
    void close_propagation();


    /*!@name Printing*/
    //@{
    void debug_print();
    void full_print();
    void print_clist(int k) ;
    virtual std::ostream& display(std::ostream&, const int current=0)  ;
    //@}
  };






  class SolverCmdLine : public TCLAP::CmdLine {

  private:
    
    TCLAP::ValueArg<std::string> *fileArg;
    TCLAP::ValueArg<int>         *seedArg;
    TCLAP::ValueArg<double>      *timeArg;
    //TCLAP::MultiArg<std::string> *printArg;
    //TCLAP::MultiArg<std::string> *nameArg;
    TCLAP::ValueArg<int>         *verbosityArg;
    TCLAP::ValueArg<int>         *randomizationArg;
    TCLAP::SwitchArg             *rewriteArg;
    TCLAP::ValueArg<std::string> *restartArg;
    TCLAP::ValueArg<double>      *factorArg;
    TCLAP::ValueArg<int>         *baseArg;
    TCLAP::ValueArg<double>      *decayArg;
    TCLAP::ValueArg<double>      *forgetArg;
    TCLAP::ValueArg<double>      *incrementArg;
    TCLAP::SwitchArg             *learningArg;
    TCLAP::ValueArg<std::string> *branchingArg;
    TCLAP::ValueArg<std::string> *orderingArg;
    TCLAP::SwitchArg             *printsolArg;
    TCLAP::SwitchArg             *printparArg;
    TCLAP::SwitchArg             *printmodArg;
    TCLAP::SwitchArg             *printinsArg;
    TCLAP::SwitchArg             *printstaArg;
    //TCLAP::ValueArg<std::string> *commentArg;
    TCLAP::ValueArg<std::string> *pcommentArg;
    TCLAP::ValueArg<std::string> *pstatArg;
    TCLAP::ValueArg<std::string> *pobjectiveArg;
    TCLAP::ValueArg<std::string> *psolutionArg;
    TCLAP::ValueArg<std::string> *poutcomeArg;
    TCLAP::SwitchArg             *allsolArg;
    

    TCLAP::ValuesConstraint<std::string> * r_allowed;
    TCLAP::ValuesConstraint<std::string> *vo_allowed;
    TCLAP::ValuesConstraint<std::string> *bo_allowed;
    // int print_sol;
    // int print_sta;
    // int print_mod;
    // int print_par;
    // int print_ins;

    // void init_print() {
    //   if(print_sol<0) {

    // 	print_sol = print_sta = print_mod = print_ins = print_par = 0;
	
    // 	// std::vector<std::string> to_print = printArg->getValue();
    // 	// for(std::vector<std::string>::iterator it=to_print.begin(); it!=to_print.end(); ++it) {
	  
    // 	//   if(*it == "sol") {
    // 	//     print_sol = 1;
    // 	//   } else if(*it == "model") {
    // 	//     print_mod = 1;
    // 	//   } else if(*it == "instance") {
    // 	//     print_ins = 1;
    // 	//   } else if(*it == "stat") {
    // 	//     print_sta = 1;
    // 	//   } else if(*it == "params") {
    // 	//     print_par = 1;
    // 	//   }
    // 	// }
    //   }
    //   std::cout << print_sol << print_mod << print_sta << print_ins << print_par << std::endl;

    // }
    
public:
    
    SolverCmdLine(const std::string& message,
		  const char delimiter = ' ',
		  const std::string& version = "none",
		  bool helpAndVersion = true) 
    // : CmdLine(message, delimiter, version, helpAndVersion)//  {
    //   initialise();
    // }
	;
    
    ~SolverCmdLine()//  {
    //   delete fileArg;
    //   delete seedArg;
    //   delete timeArg;
    //   //delete printArg;
    //   // delete printsolArg;
    //   // delete printstaArg;
    //   // delete printmodArg;
    //   // delete printparArg;
    //   // delete printinsArg;
    //   delete verbosityArg;
    //   delete randomizationArg;
    //   delete rewriteArg;
    //   delete restartArg;
    //   delete factorArg;
    //   delete baseArg;
    //   delete decayArg;
    //   delete forgetArg;
    //   delete incrementArg;
    //   delete learningArg;
    //   delete branchingArg;
    // }
      ;

  void initialise()//  {

  //   // INPUT FILE
  //   fileArg = new TCLAP::ValueArg<std::string>("f","file","instance file",true,"data/example.opb","string");
  //   add( *fileArg );

  //   // DECAY FACTOR
  //   timeArg = new TCLAP::ValueArg<double>("t","time_limit","time limit",false,0,"double");
  //   add( *timeArg );



  //   // // PRINTING OPTIONS
  //   // std::vector<std::string> pallowed;
  //   // pallowed.push_back("params");
  //   // pallowed.push_back("instance");
  //   // pallowed.push_back("model");
  //   // pallowed.push_back("sol");
  //   // pallowed.push_back("stat");
  //   // TCLAP::ValuesConstraint<std::string> p_allowed( pallowed );
  //   // printArg = new TCLAP::MultiArg<std::string>("p","print","objects to print",false,&p_allowed);
  //   // add( *printArg );

  //   // VERBOSITY LEVEL
  //   verbosityArg = new TCLAP::ValueArg<int>("v","verbosity","verbosity level",false,1,"int");
  //   add( *verbosityArg );
    
  //   // RANDOM SEED
  //   seedArg = new TCLAP::ValueArg<int>("s","seed","random seed",false,12345,"int");
  //   add( *seedArg );
    
  //   // HEURISTIC RANDOMIZATION
  //   randomizationArg = new TCLAP::ValueArg<int>("z","randomization","randomization level",false,0,"int");
  //   add( *randomizationArg );
    
  //   // WHETHER WE USE REWRITING OPTIONS
  //   rewriteArg = new TCLAP::SwitchArg("w","rewrite","use rewriting during preprocessing", false);
  //   add( *rewriteArg );
    
  //   // RESTART POLICY
  //   std::vector<std::string> rallowed;
  //   rallowed.push_back("no");
  //   rallowed.push_back("geom");
  //   rallowed.push_back("luby");
  //   TCLAP::ValuesConstraint<std::string> r_allowed( rallowed );
  //   restartArg = new TCLAP::ValueArg<std::string>("r","restart","restart policy",false,"no",&r_allowed);
  //   add( *restartArg );
    
  //   // RESTART FACTOR
  //   factorArg = new TCLAP::ValueArg<double>("m","factor","restart factor",false,1.05,"double");
  //   add( *factorArg );
    
  //   // RESTART BASE
  //   baseArg = new TCLAP::ValueArg<int>("e","base","restart base",false,200,"int");
  //   add( *baseArg );
    
  //   // DECAY FACTOR
  //   decayArg = new TCLAP::ValueArg<double>("d","decay","decay factor",false,.96,"double");
  //   add( *decayArg );
    
  //   // FORGETFULNESS
  //   forgetArg = new TCLAP::ValueArg<double>("g","forget","clause forgetfulness",false,.75,"double");
  //   add( *forgetArg );
    
  //   // ACTIVITY INCREMENT
  //   incrementArg = new TCLAP::ValueArg<double>("i","increment","activity increment",false,.012,"double");
  //   add( *incrementArg );
    
  //   // USE CLAUSE LEARNING
  //   learningArg = new TCLAP::SwitchArg("l","learning","Switch on clause learning (CDCL)", false);
  //   add( *learningArg );

  // // // PRINT SOLUTION
  // //   printsolArg = new TCLAP::SwitchArg("p", "print_sol","Print the solution, if found", false);
  // //   add( *printsolArg );

  // // // PRINT MODEL
  // //  printmodArg = new TCLAP::SwitchArg("x", "print_mod","Print the model", false);
  // //   add( *printmodArg );

  // // // PRINT STATISTICS
  // //  printstaArg = new TCLAP::SwitchArg("y", "print_sta","Print the statistics", false);
  // //   add( *printstaArg );

  // // // PRINT INSTANCE
  // //  printinsArg = new TCLAP::SwitchArg("w", "print_ins","Print the instance", false);
  // //   add( *printinsArg );

  // // // PRINT PARAMETERS
  // //  printparArg = new TCLAP::SwitchArg("q", "print_par","Print the parameters", false);
  // //   add( *printparArg );

  //   // // VARIABLE ORDERING
  //   // std::vector<std::string> voallowed;
  //   // rallowed.push_back("dom/deg");
  //   // rallowed.push_back("dom/wdeg");
  //   // rallowed.push_back("dom/pruning");
  //   // rallowed.push_back("dom/activity");
  //   // rallowed.push_back("activity");
  //   // TCLAP::ValuesConstraint<std::string> vo_allowed( voallowed );
  //   // TCLAP::ValueArg<std::string> choiceArg("c","choice","variable ordering",false,"dom/activity",&vo_allowed);
  //   // add( *choiceArg ); 

  //   // VALUE ORDERING
  //   std::vector<std::string> boallowed;
  //   boallowed.push_back("min");
  //   boallowed.push_back("max");
  //   boallowed.push_back("boand");
  //   boallowed.push_back("halfsplit");
  //   boallowed.push_back("boandsplit");
  //   boallowed.push_back("guided");
  //   TCLAP::ValuesConstraint<std::string> bo_allowed( boallowed );
  //   branchingArg = new TCLAP::ValueArg<std::string>("b","branching","value ordering",false,"guided",&bo_allowed);
  //   add( *branchingArg );    

  // }

    ;


  void set_parameters(Solver& s)//  {

  //   s.parameters.verbosity = verbosityArg->getValue();
  //   s.parameters.restart_factor = factorArg->getValue();
  //   s.parameters.restart_base = baseArg->getValue();
  //   s.parameters.restart_limit = s.parameters.restart_base;
  //   s.parameters.time_limit = timeArg->getValue();
  //   s.parameters.activity_decay = decayArg->getValue();
  //   s.parameters.backjump = learningArg->getValue();

  // }
    ;


    std::string get_value_ordering();

    std::string get_variable_ordering();

    std::string get_restart_policy();

    std::string get_filename()//  {
    //   return fileArg->getValue().c_str();
    // }
      ;
    int get_seed() // { 
    //   return seedArg->getValue();
    // }
      ;

    int get_randomization(); // { 

    bool print_model() // {
    //   //init_print();
    //   //return printmodArg->getValue();
    //   return true;
    // }
      ;

    bool print_solution() // {
    //   //init_print();
    //   //return printsolArg->getValue();
    //   return true;
    // }
      ;

    bool print_parameters() // {
    //   //init_print();
    //   //return printparArg->getValue();
    //   return true;
    // }
      ;

    bool print_instance() // {
    //   //init_print();
    //   //return printinsArg->getValue();
    //   return true;
    // }
      ;

    bool print_statistics()//  {
    //   //init_print();
    //   //return printstaArg->getValue();
    //   return true;
    // }
      ;


    bool use_rewrite();

    bool enumerate_solutions();


};



  std::ostream& operator<< (std::ostream& os, Solution& x);
  std::ostream& operator<< (std::ostream& os, Solution* x);

  std::ostream& operator<< (std::ostream& os, Solver& x);
  std::ostream& operator<< (std::ostream& os, Solver* x);

  std::ostream& operator<< (std::ostream& os, ConstraintQueue& x);
  std::ostream& operator<< (std::ostream& os, ConstraintQueue* x);

  std::ostream& operator<< (std::ostream& os, const SolverStatistics& x);
  std::ostream& operator<< (std::ostream& os, const SolverStatistics* x);

  std::ostream& operator<< (std::ostream& os, ConstraintTriggerArray& x);
  std::ostream& operator<< (std::ostream& os, ConstraintTriggerArray* x);



  SearchMonitor& operator<< (SearchMonitor& os, Variable& x);
  SearchMonitor& operator<< (SearchMonitor& os, VarArray& x);
  SearchMonitor& operator<< (SearchMonitor& os, Constraint& x);
  //SearchMonitor& operator<< (SearchMonitor& os, std::string& x);
  SearchMonitor& operator<< (SearchMonitor& os, const char* x);
  SearchMonitor& operator<< (SearchMonitor& os, const int x);


}

#endif // __SOLVER_HPP
