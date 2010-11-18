
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

/** \file sol.h
    \brief Header for the generic Solver.
*/


#ifndef _SOLVER_H
#define _SOLVER_H

#include <mistral_mod.h>
#include <mistral_glo.h>
#include <string>
#include <stdexcept>
#include <stack>
#include <iostream>
#include <fstream>
#include <string.h>

//using namespace std;

/**********************************************************************/

namespace Mistral {

  class SolutionMethod {

  public:

    Solver *solver;

    SolutionMethod(Solver *s) { solver = s; }
    virtual ~SolutionMethod() {}

    virtual void execute() = 0;
    virtual void initialise() = 0;
  };


  class abort_search : public std::exception {
  public:
    abort_search() : exception() {}
  };


  /**********************************************
   * MistralGacVarStack
   **********************************************/

  /// Stack of Variable for computing the GAC closure.
  class MistralGacVarStack {
  public:
    /**@name Parameters*/
    //@{
    VariableInt **MistralGacVarStack_;
    int *trigger;
    int     size;
    //@}

    /**@name Constructors*/
    //@{
    MistralGacVarStack()
    {
      size=0;
    }
    MistralGacVarStack(int n)
    {
      initStack(n);
    }
    ~MistralGacVarStack()
    {
      delete [] trigger;
      delete [] MistralGacVarStack_;
    }
    void initStack(int n)
    {
      trigger = new int[n];
      std::fill( trigger, trigger+n, 0 );
      MistralGacVarStack_ = new VariableInt*[n];
      size = 0;
    }
    //@}

    /**@name Accessors*/
    //@{
    inline void push(VariableInt* v, int event)
    {
      if( !trigger[v->id] )
	MistralGacVarStack_[size++] = v;
      trigger[v->id] |= event;
    }

    inline void push(VariableInt* v, VariableInt* x)
    {
      if( !trigger[v->id] )
	MistralGacVarStack_[size++] = v;
      trigger[v->id] = x->id;
    }

    inline VariableInt* pop( int& event )
    {
      int j, i = size, mindom = NOVAL;
      VariableInt *best, *last = MistralGacVarStack_[--size];
      while( i-- )
	if( mindom > MistralGacVarStack_[i]->domsize() )
	  {
	    j = i;
	    mindom = MistralGacVarStack_[i]->domsize();
	  }
      event = trigger[MistralGacVarStack_[j]->id];
      trigger[MistralGacVarStack_[j]->id] = 0;
      best = MistralGacVarStack_[j];
      MistralGacVarStack_[j] = last;
      return best;
    }

    inline VariableInt* back( int& event )
    {
      VariableInt* last = MistralGacVarStack_[size-1];
      event = trigger[last->id];
      return last;
    }

    inline bool empty() const
    {
      return (!size);
    }

    inline void clear()
    {
      while( size )
	trigger[MistralGacVarStack_[--size]->id] = 0;
    }
    //@}
  };




//   class VarValuation {
    
//   public:

//     int value;
//     char type;
    
//     VarValuation(const int v, const int t) {value=v; type=t;}
//     virtual ~VarValuation() {}

//   };


  /**********************************************
   * Solver
   **********************************************/

  /// A representation of the search process.
  /**
     Maintains data structures
     and Implements the generic backtracking method
  */
  class Solver {

  public:
    
#ifdef _DRAW_EXPLORATION
    std::ofstream* expFile;
#endif

    /*!@name Static Parameters*/
    //@{
    /// Whether the constraint are optimized for bound or ac reasoning
    bool boundOptimised;
    /// Whether we use domain splitting as branching strategy
    bool domainSplitting;
    /// Whether the variables are randomly shuffled before a restart
    bool randomizedRestart;
    /// Whether we are in branch&bound mode or not
    bool optimisation;

    /// A method to execute whenever a solution is found (optional)
    SolutionMethod *function;
    /// The objective function, for optimisation problem (optional)
    ObjectiveFunction *goal;
    /// The dynamic variable ordering heuristic
    DVO  *heuristic;
    /// The different weight and quantities "learned" during search
    Vector<Weighter*> learners;

    /// The current status (UNSAT/SAT/LIMITOUT)
    int         status;
    /// The current level in the search tree (number of left branches from the root)
    int          level;
    int     init_level; // ususally:
    //                       level 0: kept for the original problem
    //                       level 1: the 'global' level, no decision taken
    //                       level 1+k: k decisions taken
    /// The target backtrack level
    int backtrackLevel;
    /// The number of search variables
    int         length;
    /// The number of variables
    int        numvars;

    /// The instantiation sequence, an array containing all variables, ordered until the index level-1
    VariableInt   **sequence;
    /// A pointer to the first unassigned variable in sequence
    VariableInt     **future;
    /// A pointer to the end of the array 'sequence'
    VariableInt      **empty;
    /// A pointer to the variables that do not participate in the search
    VariableInt **auxilliary;


    /// The array of variables for which a decision has been taken
    Vector< VariableInt* > decision;
    Vector< SimpleUnaryConstraint > branching_decision;

    //Vector< VarValuation > valuation;
    /// An array to store the mins for last solution found
    int     *solution;
    int     *max_solution;

    /// The list of Variables
    Vector < VariableInt* > variables;
    /// The list of Constraints
    Vector < Constraint* > constraints;
    /// The list of unary Constraints
    Vector < UnaryConstraint* > unaryCons;
    /// The list of unary Constraints (new implementation)
    Vector < SimpleUnaryConstraint > sUnaryCons;
    int unaryCons_size;
    /// The list of constants used in the model
    std::map<int, VariableInt*> constants;


    /******************************************
     * Backtracking structures
     ******************************************/

  private:

    ///// Pool of Boolean domains
    //ReversibleSparseSet domainPool;
    /// Trail of reversible objects
    Vector<ReversibleObj*> store;
    /// Trail of levels
    Vector<int> lvl_;
    /// Trail of unlinked search variables
    Vector<VariableInt **> past;
    /// Trail of unlinked non-search variables
    Vector<VariableInt **> auxv;

  public:


    /// GAC Stack
    MistralGacList<VariableInt*> gacvarstack;
    /// GAC Stack
    MistralGacStack<Constraint*> gacconstack;
    /// List of unary constraints

    /// Number of solutions to find, if equal to -1, then all solutions are listed
    int FIND_ALL;
    /// Used for lds search
    int DISCREPANCY;
    /// Number of nodes, that is recursive calls to  the dfs algo
    unsigned long int NODES;
    /// Limit on the number of nodes
    unsigned long int NDSLIMIT;
    /// Number of backtracks, that is unsuccesful recursive calls
    unsigned long int BACKTRACKS;
    /// Limit on the number of backtracks, used for restarting
    unsigned long int BTSLIMIT;
    /// Number of constraint failures
    unsigned long int FAILURES;
    /// Limit on the number of failures, used for probing
    unsigned long int FAILLIMIT;
    /// Number of calls to a constraint propagator
    unsigned long int PROPAGS;
    /// Number of constraint checks (unused)
    unsigned long int CHECKS;
    /// Number of solutions found so far
    unsigned long int SOLUTIONS;
    /// Number of solutions found so far
    unsigned long int FILTERINGS;
    /// Number of solutions found so far
    unsigned long int RESTORES;
    /// flag to know if there is a limit
    unsigned int LIMIT;
    /// Number of solutions found so far
    unsigned long int MISC;

    /// Number of backtracks for each restart
    Vector< unsigned long int > BTSLIST;

#ifdef _WEIGHT_STATS
    Vector< double > rts_time;
    Vector< double > rts_node;
    Vector< double > gini_clist;
    Vector< double > gini_vlist;
    Vector< double > hoover_list;
#endif

    /// timestamp
    double INITTIME;
    /// timestamp
    double STARTTIME;
    /// timestamp
    double ENDTIME;
    /// timestamp
    double TOTTIME;
    /// timestamp
    double SOLTIME;
    /// timestamp
    double TIMELIMIT;

    /// level of verbosity
    int verbosity;

    ///
    unsigned long int fail_increment;
    //@}


    ConstraintClauseBase *sat;


    /**@name Constructors*/
    //@{
    Solver(CSP& c, VarArray& x, VariableOrdering& v);
    void init(CSP& c, VarArray& x, VariableOrdering& v);
    Solver(CSP& c, VariableOrdering& v);
    void init(CSP& c, VariableOrdering& v);
    Solver(CSP& c, VarArray& x);
    void init(CSP& c, VarArray& x);
    Solver(CSP& c, BuildObject **x, const int l);
    void init(CSP& c, BuildObject **x, const int l);
    Solver(CSP& c);
    void init(CSP& c);
    Solver();
    void init();

    bool is_built();
    void prebuild(CSP&);
    void build(CSP&, const int n=NOVAL);
    void initSearch( VarArray& X, const int s=0 );
    void initSearch( BuildObject **x, const int n, const int s );
    void initSearch( const int n=0,
		     const int s=NOVAL,
		     const bool r=true );

    void add(VariableOrdering& g);
    void add(BranchingStrategy& b);
    virtual ~Solver();
    //@}


    /**@name Solving methods*/
    //@{
    /// Returns UNSAT/SAT/LIMITOUT
    /*
      inline int iterative_dfs()
      {
      VariableInt *lastDecision;

      while( status == UNKNOWN ) {
      if( filtering() ) {
      if( future == empty ) {
      solutionFound();
      } else {
      newNode();
      }
      } else {
      if( !level ) {
      status = UNSAT;
      } else if( limitsExpired() ) {
      status = LIMITOUT;
      } else {
      lastDecision = decision[level];
      backtrackTo( backtrackLevel );
      lastDecision->branch->right();
      }
      }
      }
      return status;
      }
    */


    
    inline void newNode( )
    {
      if(level++ > init_level) learnSuccess();

      past.push( future );
      auxv.push( auxilliary  );
      lvl_.push( store.size );


      VariableInt *x = heuristic->select();
      decision.push( x );

      SimpleUnaryConstraint dec( x );
      branching_decision.push( dec );

      //decision[level]->branch->left();

      branching_decision[level].make();
      branching_decision[level].left();

#ifdef _DEBUGSEARCH
      if(verbosity > 2) {
	std::cout << "c";
	for(int k=0; k<=level; ++k) std::cout << " ";
	branching_decision[level].print(std::cout);
	std::cout << std::endl;
      }
#endif

#ifdef _DRAW_EXPLORATION
      if((*expFile)) {
	(*expFile) << "+ newNode\n";
	branching_decision[level].print(*expFile);
      }
#endif

      ++NODES;
      int i=learners.size;
      while( i-- )
	learners[i]->notifyChoice( );
    }


    inline void newNode(VariableInt *currentDecision)
    {
      if(level++ > init_level) learnSuccess();

      past.push( future );
      auxv.push( auxilliary  );
      lvl_.push( store.size );

      SimpleUnaryConstraint dec(currentDecision);
      branching_decision.push( dec );

      decision.push( currentDecision );

      //decision[level]->branch->left();

      branching_decision[level].make();
      branching_decision[level].left();

#ifdef _DEBUGSEARCH
      if(verbosity > 2) {
	std::cout << "c";
	for(int k=0; k<=level; ++k) std::cout << " ";
	branching_decision[level].print(std::cout);
	std::cout << std::endl;
      }
#endif

#ifdef _DRAW_EXPLORATION
      if((*expFile)) {
	(*expFile) << "+ newNode\n";
	branching_decision[level].print(*expFile);
      }
#endif

      ++NODES;
      int i=learners.size;
      while( i-- )
	learners[i]->notifyChoice( );
    }


    inline void save()
    {
      if(level++ > init_level) learnSuccess();
      
      past.push( future );
      auxv.push( auxilliary  );
      lvl_.push( store.size );

      ++NODES;
      decision.push(NULL);
      SimpleUnaryConstraint no_decision;
      branching_decision.push( no_decision );
    }


    inline void reverseNewNode( )
    {
      if(level++ > init_level) learnSuccess();

      past.push( future );
      auxv.push( auxilliary  );
      lvl_.push( store.size );

      decision.push( heuristic->select() );


#ifdef _DEBUGSEARCH
      std::cout << "c";
      for(int k=0; k<=level; ++k) std::cout << " ";
      decision[level]->print(std::cout);
#endif

#ifdef _DRAW_EXPLORATION
      if((*expFile)) {
	(*expFile) << "+ reverseNewNode\n";
	decision[level]->print(*expFile);
      }
#endif

      decision[level]->branch->reverse_left();

#ifdef _DEBUGSEARCH
      decision[level]->branch->printRight( std::cout );
      std::cout << std::endl;
#endif

#ifdef _DRAW_EXPLORATION
      if((*expFile)) {
	decision[level]->branch->printRight(*expFile);
	(*expFile) << std::endl;
      }
#endif

      ++NODES;
      int i=learners.size;
      while( i-- )
	learners[i]->notifyChoice( );
    }

    inline void reverseNewNode(VariableInt *currentDecision)
    {
      if(level++ > init_level) learnSuccess();

      past.push( future );
      auxv.push( auxilliary  );
      lvl_.push( store.size );

      decision.push( currentDecision );


#ifdef _DEBUGSEARCH
      std::cout << "c";
      for(int k=0; k<=level; ++k) std::cout << " ";
      decision[level]->print(std::cout);
#endif

#ifdef _DRAW_EXPLORATION
      if((*expFile)) {
	(*expFile) << "+ reverseNewNode\n";
	decision[level]->print(*expFile);
      }
#endif

      decision[level]->branch->reverse_left();

#ifdef _DEBUGSEARCH
      decision[level]->branch->printRight( std::cout );
      std::cout << std::endl;
#endif

#ifdef _DRAW_EXPLORATION
      if((*expFile)) {
	decision[level]->branch->printRight(*expFile);
	(*expFile) << std::endl;
      }
#endif

      ++NODES;
      int i=learners.size;
      while( i-- )
	learners[i]->notifyChoice( );
    }


    inline int iterative_dfs()
    {
      //VariableInt *lastDecision;
      SimpleUnaryConstraint last_decision;
      
      while( status == UNKNOWN ) {

	//checkDecisions();

	if( filtering() ) {

// 	  for(int i=0; i<length; ++i) {
// 	    variables[i]->print(std::cout);
// 	    std::cout << std::endl;
// 	  }
	  
	  if( future == empty ) {
	    solutionFound(init_level);
	  } else {

// 	    variables[1]->print(std::cout);
// 	    std::cout << " " << variables[1]->weight << std::endl;
// 	    variables[18]->print(std::cout);
// 	    std::cout << " " << variables[18]->weight << std::endl;

// 	    if(verbosity > 2) {
// 	      for(int i=0; i<(empty-future); ++i) {
// 		std::cout << "future: ";
// 		future[i]->print(std::cout);
// 		std::cout << heuristic->get_value(future[i]) << std::endl;
// 	      }
// 	    }


	    newNode();
	  }
	} else {
	  
	  //	  std::cout << "FAIL" << std::endl;

	  if( level <= init_level ) {
	    
#ifdef _DEBUGSEARCH
	    if(verbosity > 2) {
	      std::cout << "c UNSAT!" << std::endl;
	    }
#endif
	    
	    status = UNSAT;
	  } else if( limitsExpired() ) {
	    
#ifdef _DEBUGSEARCH
	    if(verbosity > 2) {
	      std::cout << "c";
	      for(int k=0; k<=level; ++k) std::cout << " ";
	      SimpleUnaryConstraint d = branching_decision[level];
	      d.revert();
	      d.print(std::cout);
	      std::cout << " (limit expired at level " << level << ")" << std::endl;
	    }
#endif
	    
	    status = LIMITOUT;
	  } else {
	    //lastDecision = decision[level];
	    last_decision = branching_decision[level];

#ifdef _DEBUGSEARCH
	    if(verbosity > 2) {
	      if( level > backtrackLevel+1 ) {
		std::cout << "c";
		for(int k=0; k<=level; ++k) std::cout << " ";
		std::cout << " backjump to level " << backtrackLevel << std::endl;
	      }
	    }
#endif
	    
	    backtrackTo( backtrackLevel );
	    
	    last_decision.right();
	    //lastDecision->branch->right();
	    //lastDecision->branch->left();
	    //lastDecision->branch->reverse_right();

#ifdef _DEBUGSEARCH
	    if(verbosity > 2) {
	      std::cout << "c";
	      for(int k=0; k<=level; ++k) std::cout << " ";
	      last_decision.print( std::cout );
	      std::cout << std::endl;
	    }
#endif

	  }
	}
      }
      return status;
    }


    int pseudoRldsProbe(int K, ReversibleNum<int>& numDiscrepancies)
    {
      std::cout << std::endl << "pseudoRldsProbe -=- K=" << K << std::endl;
      //		std::cout << "1 - ";

#ifdef _DRAW_EXPLORATION
      std::stringstream ossK;
      ossK << K;
      std::string stringK = ossK.str();
      expFile = new std::ofstream(("./pseudoRldsProbe_" + stringK + "_exploration.txt").c_str(), std::ios::out | std::ios::trunc);
      if(!(*expFile)) {
	std::cout << "c fail to create exploration file\n";
      }
#endif

      VariableInt *lastDecision;
      bool enoughDiscrepanciesAllowed = true;

      while( status == UNKNOWN ) {
	// Fail either if
	//  - the maximum number of discrepancies is reached
	//  - the problem is inconsistent
	if( numDiscrepancies <= K && filtering() ) {
	  if( future == empty ) {
	    solutionFound(init_level);
	  }
	  else {
	    newNode();
	  }
	} else {
	  if( level <= init_level ) {
	    // distinguish if we exit because of lds or not
	    if(enoughDiscrepanciesAllowed) {
	      status = UNSAT;
	    }
	    break;
	  } else if( limitsExpired() ) {
	    status = LIMITOUT;
	    //					std::cout << "limitout" << std::endl;
	  } else {
	    lastDecision = decision[level];
	    backtrackTo( level-1 );
	    ++numDiscrepancies;

	    // do the right branch only if the threshold is not met
	    if(numDiscrepancies <= K) {
#ifdef _DEBUGSEARCH
	      std::cout << "c";
	      for(int k=0; k<=level; ++k) std::cout << " ";
	      lastDecision->print(std::cout);
#endif
	      lastDecision->branch->right();
#ifdef _DEBUGSEARCH
	      lastDecision->branch->printRight( std::cout );
	      std::cout << std::endl;
#endif
	    } else {
	      enoughDiscrepanciesAllowed = false;
	    }
	  }
	}
      }
      return status;
    }

    int pseudoRldsProbe_level(unsigned int K, ReversibleIntList& discrepancy_level)
    {
      //		std::cout << "\n\npseudoRldsProbe_level -=- K=" << K << std::endl;

#ifdef _DRAW_EXPLORATION
      std::stringstream ossK;
      ossK << K;
      std::string stringK = ossK.str();
      expFile = new std::ofstream(("./pseudoRldsProbe_level_" + stringK + "_exploration.txt").c_str(), std::ios::out | std::ios::trunc);
      if(!(*expFile)) {
	std::cout << "c fail to create exploration file\n";
      }
#endif

      VariableInt *lastDecision;
      bool enoughDiscrepanciesAllowed = true;

      while( status == UNKNOWN ) {
	// Fail either if
	//  - the maximum number of discrepancies is reached
	//  - the problem is inconsistent
	if( discrepancy_level.size <= K && filtering() ) {
	  if( future == empty ) solutionFound(init_level);
	  else newNode();
	} else {
	  if( discrepancy_level.size > K) {
	    enoughDiscrepanciesAllowed = false;
	  }
	  if( level <= init_level ) {
	    // distinguish if we exit because of lds or not
	    if(enoughDiscrepanciesAllowed) {
	      status = UNSAT;
	    }
	    break;
	  } else if( limitsExpired() ) {
	    status = LIMITOUT;
	  } else {
	    lastDecision = decision[level];
	    backtrackTo( level-1 );
	    if (!discrepancy_level.member(level)) {
	      discrepancy_level.insert(level);
	    }

	    // do the right branch only if the threshold is not met
	    if(discrepancy_level.size <= K) {
#ifdef _DEBUGSEARCH
	      std::cout << "c";
	      for(int k=0; k<=level; ++k) std::cout << " ";
	      lastDecision->print(std::cout);
#endif
	      lastDecision->branch->right();
#ifdef _DEBUGSEARCH
	      lastDecision->branch->printRight( std::cout );
	      std::cout << std::endl;
#endif
	    }
	  }
	}
      }
      return status;
    }

    int pseudoRldsProbe_variable(unsigned int K, ReversibleIntList& discrepancy_variable)
    {
      //		std::cout << "\n\npseudoRldsProbe_variable -=- K=" << K << std::endl;

#ifdef _DRAW_EXPLORATION
      std::stringstream ossK;
      ossK << K;
      std::string stringK = ossK.str();
      expFile = new std::ofstream(("./pseudoRldsProbe_variable_" + stringK + "_exploration.txt").c_str(), std::ios::out | std::ios::trunc);
      if(!(*expFile)) {
	std::cout << "c fail to create exploration file\n";
      }
#endif

      VariableInt *lastDecision;
      bool enoughDiscrepanciesAllowed = true;

      while( status == UNKNOWN ) {
	// Fail either if
	//  - the maximum number of discrepancies is reached
	//  - the problem is inconsistent
	if( discrepancy_variable.size <= K && filtering() ) {
	  if( future == empty ) solutionFound(init_level);
	  else newNode();
	} else {
	  if( discrepancy_variable.size > K) {
	    enoughDiscrepanciesAllowed = false;
	  }
	  if( level <= init_level ) {
	    // distinguish if we exit because of lds or not
	    if(enoughDiscrepanciesAllowed) {
	      status = UNSAT;
	    }
	    break;
	  } else if( limitsExpired() ) {
	    status = LIMITOUT;
	  } else {
	    lastDecision = decision[level];
	    backtrackTo( level-1 );
	    if (!discrepancy_variable.member(lastDecision->id)) {
	      discrepancy_variable.insert(lastDecision->id);
	    }

	    // do the right branch only if the threshold is not met
	    if(discrepancy_variable.size <= K) {
#ifdef _DEBUGSEARCH
	      std::cout << "c";
	      for(int k=0; k<=level; ++k) std::cout << " ";
	      lastDecision->print(std::cout);
#endif
	      lastDecision->branch->right();
#ifdef _DEBUGSEARCH
	      lastDecision->branch->printRight( std::cout );
	      std::cout << std::endl;
#endif
	    }
	  }
	}
      }
      return status;
    }

    int ldsProbe(int K, ReversibleNum<int>& numDiscrepancies)
    {
      //		std::cout << "\n\nldsProbe -=- K=" << K << std::endl;

#ifdef _DRAW_EXPLORATION
      std::stringstream ossK;
      ossK << K;
      std::string stringK = ossK.str();
      expFile = new std::ofstream(("./ldsProbe_" + stringK + "_exploration.txt").c_str(), std::ios::out | std::ios::trunc);
      if(!(*expFile)) {
	std::cout << "c fail to create exploration file\n";
      }
#endif

      VariableInt *lastDecision;
      bool enoughDiscrepanciesAllowed = true;

      while( status == UNKNOWN ) {
	// Fail either if
	//  - the maximum number of discrepancies is reached
	//  - the problem is inconsistent
	if( numDiscrepancies <= K && filtering() ) {
	  if( future == empty ) solutionFound(init_level);
	  else {
	    reverseNewNode();
	    ++numDiscrepancies;
	  }
	} else {
	  if (numDiscrepancies > K) {
	    enoughDiscrepanciesAllowed = false;
	  }
	  if( level <= init_level ) {
	    // distinguish if we exit because of lds or not
	    if(enoughDiscrepanciesAllowed) {
	      status = UNSAT;
	    }
	    break;
	  } else if( limitsExpired() ) {
	    status = LIMITOUT;
	  } else {
	    lastDecision = decision[level];
	    backtrackTo( level-1 );
#ifdef _DEBUGSEARCH
	    std::cout << "c";
	    for(int k=0; k<=level; ++k) std::cout << " ";
	    lastDecision->print(std::cout);
#endif
	    lastDecision->branch->reverse_right();
#ifdef _DEBUGSEARCH
	    lastDecision->branch->printLeft( std::cout );
	    std::cout << std::endl;
#endif
	  }
	}
      }
      return status;
    }

    int ldsStackProbe(int K, ReversibleNum<int>& numDiscrepancies)
    {
      //		std::cout << "\n\nldsStackProbe -=- K=" << K << std::endl;

#ifdef _DRAW_EXPLORATION
      std::stringstream ossK;
      ossK << K;
      std::string stringK = ossK.str();
      expFile = new std::ofstream(("./ldsStackProbe_" + stringK + "_exploration.txt").c_str(), std::ios::out | std::ios::trunc);
      if(!(*expFile)) {
	std::cout << "c fail to create exploration file\n";
      }
#endif

      std::stack<int> backtrackLevels;
      VariableInt *lastDecision;
      bool enoughDiscrepanciesAllowed = true;

      while( status == UNKNOWN ) {
	// Fail either if
	//  - the maximum number of discrepancies is reached
	//  - the problem is inconsistent
	if( filtering() ) {
	  if( future == empty ) solutionFound(init_level);
	  else {
	    if( numDiscrepancies < K ){
	      backtrackLevels.push(level);
	      reverseNewNode();
	      ++numDiscrepancies;
	    } else {
	      enoughDiscrepanciesAllowed = false;
	      newNode();
	    }
	  }
	} else {
	  if(backtrackLevels.empty()) {
	    // distinguish if we exit because of lds or not
	    if(enoughDiscrepanciesAllowed) {
	      status = UNSAT;
	    }
	    break;
	  } else if( limitsExpired() ) {
	    status = LIMITOUT;
	  } else {

	    int levelForBacktrack = backtrackLevels.top();
	    backtrackLevels.pop();
	    lastDecision = decision[levelForBacktrack+1];
	    backtrackTo( levelForBacktrack );
#ifdef _DEBUGSEARCH
	    std::cout << "c";
	    for(int k=0; k<=level; ++k) std::cout << " ";
	    lastDecision->print(std::cout);
#endif
	    lastDecision->branch->reverse_right();
#ifdef _DEBUGSEARCH
	    lastDecision->branch->printLeft( std::cout );
	    std::cout << std::endl;
#endif
	  }
	}
      }
      return status;
    }

    int ldsDeltaProbe(int K, ReversibleNum<int>& numDiscrepancies)
    {
      //		std::cout << "\n\nldsDeltaProbe -=- K=" << K << std::endl;

#ifdef _DRAW_EXPLORATION
      std::stringstream ossK;
      ossK << K;
      std::string stringK = ossK.str();
      expFile = new std::ofstream(("./ldsDeltaProbe_" + stringK + "_exploration.txt").c_str(), std::ios::out | std::ios::trunc);
      if(!(*expFile)) {
	std::cout << "c fail to create exploration file\n";
      } else {
	(*expFile) << init_level << std::endl;
      }
#endif

      int deltaBacktrack = 1;
      VariableInt *lastDecision;
      bool enoughDiscrepanciesAllowed = true;

      while( status == UNKNOWN ) {
	// Fail either if
	//  - the maximum number of discrepancies is reached
	//  - the problem is inconsistent
	if( filtering() ) {
	  if( future == empty ) solutionFound(init_level);
	  else {
	    if( numDiscrepancies < K ){
	      reverseNewNode();
	      deltaBacktrack = 1;
	      ++numDiscrepancies;
	    } else {
	      enoughDiscrepanciesAllowed = false;
	      newNode();
	      deltaBacktrack++;
	    }
	  }
	} else {
	  if( limitsExpired() ) {
	    status = LIMITOUT;
	  } else {
	    int levelForBacktrack = level - deltaBacktrack;
	    if( levelForBacktrack < init_level ) {
	      // distinguish if we exit because of lds or not
	      if(enoughDiscrepanciesAllowed) {
		status = UNSAT;
	      }
	      break;
	    }
	    lastDecision = decision[levelForBacktrack+1];

#ifdef _DRAW_EXPLORATION
	    if((*expFile)) {
	      (*expFile) << "- backtrackTo\n";
	      (*expFile) << levelForBacktrack;
	      (*expFile) << std::endl;
	    }
#endif
	    backtrackTo( levelForBacktrack );

#ifdef _DEBUGSEARCH
	    std::cout << "c";
	    for(int k=0; k<=level; ++k) std::cout << " ";
	    lastDecision->print(std::cout);
#endif

#ifdef _DRAW_EXPLORATION
	    if((*expFile)) {
	      (*expFile) << "+ reverse_right\n";
	      lastDecision->print(*expFile);
	    }
#endif
	    lastDecision->branch->reverse_right();
#ifdef _DEBUGSEARCH
	    lastDecision->branch->printLeft( std::cout );
	    std::cout << std::endl;
#endif

#ifdef _DRAW_EXPLORATION
	    if((*expFile)) {
	      lastDecision->branch->printLeft(*expFile);
	      (*expFile) << std::endl;
	    }
#endif
	  }
	}
      }
#ifdef _DRAW_EXPLORATION
      (*expFile).flush();
      (*expFile).close();
      delete expFile;
#endif
      return status;
    }

    int ldsProbe_level(unsigned int K,ReversibleIntList& discrepancy_levels, ReversibleNum<int>& myLevel)
    {
      //		std::cout << "\n\nldsProbe_level -=- K=" << K << std::endl;

#ifdef _DRAW_EXPLORATION
      std::stringstream ossK;
      ossK << K;
      std::string stringK = ossK.str();
      expFile = new std::ofstream(("./ldsProbe_level_" + stringK + "_exploration.txt").c_str(), std::ios::out | std::ios::trunc);
      if(!(*expFile)) {
	std::cout << "c fail to create exploration file\n";
      }
#endif

      VariableInt *lastDecision;
      bool enoughDiscrepanciesAllowed = true;
      //		int myLevel = 1;
      int deltaBacktrack = 1;

      while( status == UNKNOWN ) {
	if( filtering() ) {
	  if( future == empty ) {
	    solutionFound(init_level);
	  }
	  else {
	    if (discrepancy_levels.member(myLevel)) {
	      reverseNewNode();
	      deltaBacktrack = 1;
	    }else if (discrepancy_levels.size < K) {
	      reverseNewNode();
	      discrepancy_levels.insert(myLevel);
	      deltaBacktrack = 1;
	    }else{
	      enoughDiscrepanciesAllowed = false;
	      newNode();
	      deltaBacktrack++;
	    }
	  }
	} else {
	  if( limitsExpired() ) {
	    status = LIMITOUT;
	  } else {
	    int levelForBacktrack = level - deltaBacktrack;
	    if( levelForBacktrack < init_level ) {
	      // distinguish if we exit because of lds or not
	      if(enoughDiscrepanciesAllowed) {
		status = UNSAT;
	      }
	      break;
	    }
	    lastDecision = decision[levelForBacktrack+1];
	    backtrackTo( levelForBacktrack );
	    deltaBacktrack = 1;
	    ++myLevel;
#ifdef _DEBUGSEARCH
	    std::cout << "c";
	    for(int k=0; k<=level; ++k) std::cout << " ";
	    lastDecision->print(std::cout);
#endif
	    lastDecision->branch->reverse_right();
#ifdef _DEBUGSEARCH
	    lastDecision->branch->printLeft( std::cout );
	    std::cout << std::endl;
#endif
	  }
	}
      }
      return status;
    }

    int ldsProbe_variable(unsigned int K, ReversibleIntList& discrepancy_variable)
    {
      //		std::cout << "\n\nldsProbe_variable -=- K=" << K << std::endl;

#ifdef _DRAW_EXPLORATION
      std::stringstream ossK;
      ossK << K;
      std::string stringK = ossK.str();
      expFile = new std::ofstream(("./ldsProbe_variable" + stringK + "_exploration.txt").c_str(), std::ios::out | std::ios::trunc);
      if(!(*expFile)) {
	std::cout << "c fail to create exploration file\n";
      }
#endif

      std::stack<int> backtrackLevels;
      VariableInt *currentDecision;
      bool enoughDiscrepanciesAllowed = true;

      while( status == UNKNOWN ) {
	if( filtering() ) {
	  if( future == empty ) solutionFound(init_level);
	  else {
	    currentDecision = heuristic->select();
	    if (discrepancy_variable.member(currentDecision->id)) {
	      backtrackLevels.push(level);
	      reverseNewNode(currentDecision);
	    } else if( discrepancy_variable.size < K ){
	      backtrackLevels.push(level);
	      reverseNewNode(currentDecision);
	      discrepancy_variable.insert(currentDecision->id);
	    } else {
	      enoughDiscrepanciesAllowed = false;
	      newNode(currentDecision);
	    }
	  }
	} else {
	  if(backtrackLevels.empty()) {
	    // distinguish if we exit because of lds or not
	    if(enoughDiscrepanciesAllowed) {
	      status = UNSAT;
	    }
	    else break;
	  } else if( limitsExpired() ) {
	    status = LIMITOUT;
	  } else {
	    int levelForBacktrack = backtrackLevels.top();
	    backtrackLevels.pop();
	    currentDecision = decision[levelForBacktrack+1];
	    backtrackTo( levelForBacktrack );

#ifdef _DEBUGSEARCH
	    std::cout << "c";
	    for(int k=0; k<=level; ++k) std::cout << " ";
	    currentDecision->print(std::cout);
#endif
	    currentDecision->branch->reverse_right();
#ifdef _DEBUGSEARCH
	    currentDecision->branch->printLeft( std::cout );
	    std::cout << std::endl;
#endif
	  }
	}
      }
      return status;
    }

    int ddsProbe(signed int K, ReversibleNum<int>& numDiscrepancies)
    {
      //		std::cout << "\n\nddsProbe -=- K=" << K << std::endl;

#ifdef _DRAW_EXPLORATION
      std::stringstream ossK;
      ossK << K;
      std::string stringK = ossK.str();
      expFile = new std::ofstream(("./ddsProbe_" + stringK + "_exploration.txt").c_str(), std::ios::out | std::ios::trunc);
      if(!(*expFile)) {
	std::cout << "c fail to create exploration file\n";
      }
#endif

      int deltaBacktrack = 1;
      bool enoughDiscrepanciesAllowed = true;
      VariableInt *lastDecision;

      while( status == UNKNOWN ) {
	if( filtering() ) {
	  if( future == empty ) {
	    solutionFound(init_level);
	  } else {
	    if (level-init_level+numDiscrepancies > K-1) {
	      enoughDiscrepanciesAllowed = false;
	      newNode();
	      deltaBacktrack++;
	    } else {
	      newNode();
	      deltaBacktrack = 1;
	    }
	  }
	} else {
	  if( limitsExpired() ) {
	    status = LIMITOUT;
	  } else {
	    int levelForBacktrack = level - deltaBacktrack;

	    if( levelForBacktrack < init_level ) {
	      // distinguish if we exit because of lds or not
	      if(enoughDiscrepanciesAllowed) {
		status = UNSAT;
	      }
	      break;
	    }
	    lastDecision = decision[levelForBacktrack+1];
	    backtrackTo( levelForBacktrack );
#ifdef _DEBUGSEARCH
	    std::cout << "c";
	    for(int k=0; k<=level; ++k) std::cout << " ";
	    lastDecision->print(std::cout);
#endif
	    lastDecision->branch->right();
#ifdef _DEBUGSEARCH
	    lastDecision->branch->printRight( std::cout );
	    std::cout << std::endl;
#endif
	    ++numDiscrepancies;
	  }
	}
      }
      return status;
    }

    int ddsProbe_level(signed int K)
    {
      //		std::cout << "\n\nddsProbe_level -=- K=" << K << std::endl;

#ifdef _DRAW_EXPLORATION
      std::stringstream ossK;
      ossK << K;
      std::string stringK = ossK.str();
      expFile = new std::ofstream(("./ddsProbe_level_" + stringK + "_exploration.txt").c_str(), std::ios::out | std::ios::trunc);
      if(!(*expFile)) {
	std::cout << "c fail to create exploration file\n";
      }
#endif

      std::stack<int> backtrackLevels;
      bool enoughDiscrepanciesAllowed = true;
      VariableInt *lastDecision;

      while( status == UNKNOWN ) {
	if( filtering() ) {
	  if( future == empty ) {
	    solutionFound(init_level);
	  } else {
	    if (level > init_level+K-1) {
	      enoughDiscrepanciesAllowed = false;
	      newNode();
	    } else {
	      backtrackLevels.push(level);
	      newNode();
	    }
	  }
	} else {
	  if( limitsExpired() ) {
	    status = LIMITOUT;
	  } else {
	    if (backtrackLevels.empty()) {
	      if (enoughDiscrepanciesAllowed) {
		status = UNSAT;
	      }
	      break;
	    } else {
	      int levelForBacktrack = backtrackLevels.top();
	      backtrackLevels.pop();
	      lastDecision = decision[levelForBacktrack+1];
	      backtrackTo( levelForBacktrack );
#ifdef _DEBUGSEARCH
	      std::cout << "c";
	      for(int k=0; k<=level; ++k) std::cout << " ";
	      lastDecision->print(std::cout);
#endif
	      lastDecision->branch->right();
#ifdef _DEBUGSEARCH
	      lastDecision->branch->printRight( std::cout );
	      std::cout << std::endl;
#endif
	    }
	  }
	}
      }
      return status;
    }

    int ddsProbe_variable(signed int K, ReversibleIntList& discrepancy_variable)
    {
      //		std::cout << "\n\nddsProbe_variable -=- K=" << K << std::endl;

#ifdef _DRAW_EXPLORATION
      std::stringstream ossK;
      ossK << K;
      std::string stringK = ossK.str();
      expFile = new std::ofstream(("./ddsProbe_variable_" + stringK + "_exploration.txt").c_str(), std::ios::out | std::ios::trunc);
      if(!(*expFile)) {
	std::cout << "c fail to create exploration file\n";
      }
#endif

      std::stack<int> backtrackLevels;
      bool enoughDiscrepanciesAllowed = true;
      VariableInt *currentDecision;

      while( status == UNKNOWN ) {
	if( filtering() ) {
	  if( future == empty ) {
	    solutionFound(init_level);
	  } else {
	    currentDecision = heuristic->select();
	    if (level-init_level< K
		&& (discrepancy_variable.member(currentDecision->id)
		    || discrepancy_variable.size < (unsigned int)K)) {
	      backtrackLevels.push(level);
	    } else {
	      enoughDiscrepanciesAllowed = false;
	    }
	    newNode(currentDecision);
	  }
	} else {
	  if(backtrackLevels.empty()) {
	    // distinguish if we exit because of lds or not
	    if(enoughDiscrepanciesAllowed) {
	      status = UNSAT;
	    }
	    break;
	  } else if( limitsExpired() ) {
	    status = LIMITOUT;
	  } else {
	    int levelForBacktrack = backtrackLevels.top();
	    backtrackLevels.pop();
	    if( levelForBacktrack < init_level ) {
	      // distinguish if we exit because of lds or not
	      if(enoughDiscrepanciesAllowed) {
		status = UNSAT;
	      } else {
		break;
	      }
	    }
	    currentDecision = decision[levelForBacktrack+1];
	    backtrackTo( levelForBacktrack );
#ifdef _DEBUGSEARCH
	    std::cout << "c";
	    for(int k=0; k<=level; ++k) std::cout << " ";
	    currentDecision->print(std::cout);
#endif
	    currentDecision->branch->right();
#ifdef _DEBUGSEARCH
	    currentDecision->branch->printRight( std::cout );
	    std::cout << std::endl;
#endif
	    if (!discrepancy_variable.member(currentDecision->id)) {
	      discrepancy_variable.insert(currentDecision->id);
	    }
	  }
	}
      }
      return status;
    }


    void checkDecisions() {
      
      if(decision.size != branching_decision.size) {
	std::cout << level << " inconsistency (size)" << std::endl;
	std::cout << decision.size << " / "
		  << branching_decision.size << std::endl;
      }

      for(int i=0; i<decision.size; ++i) {
	if(decision[i] != branching_decision[i].var) {
	  std::cout << i << " " << level << " inconsistency (var)" << std::endl;
	  if(decision[i]) decision[i]->print(std::cout);
	  else std::cout << "null";
	  std::cout << " / " ;
	  if(branching_decision[i].var) branching_decision[i].var->print(std::cout);
	  else std::cout << "null" << std::endl;
	}

	if(decision[i]) {
	  if(decision[i]->value() != branching_decision[i].value()) {
	    std::cout << i << " " << level << " inconsistency (val)" << std::endl;
	  std::cout << decision[i]->value() << " / "
		    << branching_decision[i].value() << std::endl;
	}
	}
      }

    }

    /*!
      Decomposed dfs
    */
    int startNewSearch();
    /*!
      Decomposed dfs
    */
    int getNextSolution();
    /*!
      Preprocessing step (closure computed with filtering())
    */
    bool preprocess();

    /// Returns UNSAT/SAT/LIMITOUT
    /*!
      presolving method
    */
    int presolve();
    /*!
      Simple solving method
    */
    int solve();
    int pseudoRldsSolve(int deltaK);
    int pseudoRldsSolve_level(int deltaK);
    int pseudoRldsSolve_variable(int deltaK);
    int ldsSolve(int deltaK);
    int ldsStackSolve(int deltaK);
    int ldsDeltaSolve(int deltaK);
    int ldsSolve_level(int deltaK);
    int ldsSolve_variable(int deltaK);
    int ddsSolve(int deltaK);
    int ddsSolve_level(int deltaK);
    int ddsSolve_variable(int deltaK);
    /*!
      LDS solving method
    */
    int ldSolve( const int step=1, const int limit=NOVAL );
    int ldSolve( const int*, const int step=1, const int limit=NOVAL );
    int ldSolve( const int, BuildObject**, int*, const int step=1,
		 const int limit=NOVAL );
    int ldSolve( VarArray&, int*, const int step=1, const int limit=NOVAL );
    int ldSolve( ConstraintLDS*, const int step=1, const int limit=NOVAL );
    bool restart(const double decay=0.0, const int reinit=-1);
    /// Returns UNSAT/SAT/LIMITOUT
    /*!
      Solve using restarts. There is no randomization besides
      that induced by the variables and values ordering.
    */
    int solve_and_restart( const int policy = GEOMETRIC,
			   const unsigned int base = 32,
			   const double factor = 1.3333333,
			   const double decay = 0.0,
			   const int reinit=-1);
    /*!
      Diarmuid's random probes
    */
    int random_probe( const unsigned int iterations = 100,
		      const unsigned int limit = 30 );
    /*!
      singleton arc consistency
    */
    int sacPreprocess( const bool complete = false,
		       const int wtype = (Weighter::IPT | Weighter::WLD) );
    /// GAC closure
    /*!
      Implementation of the AC3 algorithm
      return the level to which we should backtrack
    */
    bool filtering();
    /*!
      Called when a solution is found.
    */
    int solutionFound(int init_level);
    void store_solution();
    /*!
      Check the solution for consistency
    */
    int checkSolution();
    /*!
      Scale down the weights
    */
    inline int weightDecay( double decay)
    {
      int old_weight, new_weight, i=constraints.size, j;
      int threshold = (int)ceil(1.0 / (1.0 - decay));

      while( i-- )
	{
	  old_weight = constraints[i]->weight;
	  if( old_weight >= threshold )
	    {
	      new_weight = (int)ceil(decay * (double)old_weight);
	      j = constraints[i]->arity;
	      while( j-- )
		constraints[i]->scope[j]->weight += (new_weight - old_weight);
	    }
	}

      return 1;
    }
    /*!
      Call the individual restore() methods of Variable
      and Constraint
    */
    void restore()
    {

      decision.pop();
      branching_decision.pop();
      VariableInt **vp, **va;
      int lvl;

//       if(level == 175) {
// 	std::cout << 11 << std::endl;
// 	test_x60();
//       }

      past.pop( vp );
      auxv.pop( va );
      lvl_.pop( lvl );

//       if(level == 175) {
// 	std::cout << 22 << std::endl;
// 	test_x60();
//       }

      while( store.size > lvl ) {
	store.back()->restore();
	store.pop();
      }

//       if(level == 175) {
// 	std::cout << 33 << std::endl;
// 	test_x60();
//       }


      while(vp != future) {

// 	if(level == 175) {
// 	  std::cout  << "link ";
// 	  future[-1]->print( std::cout);
// 	  std::cout << std::endl;
// 	}
	
	future[-1]->link( );

 	//test_x60();


      }



      while(va != auxilliary)
	auxilliary[-1]->link( );

//       if(level == 175) {
// 	std::cout << 55 << std::endl;
// 	test_x60();
//       }


    }

    void upOneLevel();
    void reset_stats();
    void reset_trail(const bool full=false);
    void reset(const bool full=false);
    /// Undo search
    inline void backtrackTo( const int lvl )
    {
      while( level > lvl ) {
	++BACKTRACKS;
	restore();
	--level;
      }
    }

    /// Triggers an event e for the variable x
    inline void triggerEvent(const int x, const int e)
    {
      gacvarstack.push( x, e );
    }

    /// Check if one limit has expired
    inline bool limitsExpired()
    {

#ifdef _DEBUGSEARCH
      if(verbosity > 5) {
      if(TIMELIMIT > .0 && getRunTime() - STARTTIME >= TIMELIMIT)
	std::cout << "TIME LIMIT " << TIMELIMIT << " " << (getRunTime() - STARTTIME) << std::endl;

      if(BTSLIMIT  >  0 && BACKTRACKS   >=  BTSLIMIT)
	std::cout << "BTS LIMIT " << BACKTRACKS << " " << BTSLIMIT << std::endl;

      if(NDSLIMIT  >  0 && NODES   >=  NDSLIMIT)
	std::cout << "NDS LIMIT " << NODES << " " << NDSLIMIT << std::endl;

      if(FAILLIMIT == 1 && FAILURES   >=  FAILLIMIT)
	std::cout << "FAIL LIMIT " << FAILURES << " " << FAILLIMIT << std::endl;
      }
#endif

      //		std::cout << "TIMELIMIT=" << TIMELIMIT << "; getRunTime()-STARTTIME=" << getRunTime() - STARTTIME << std::endl;

      return (
	      (TIMELIMIT > .0 && getRunTime() - STARTTIME >= TIMELIMIT) ||
	      (BTSLIMIT  >  0 && BACKTRACKS   >=  BTSLIMIT)      ||
	      (NDSLIMIT  >  0 && NODES   >=  NDSLIMIT)      ||
	      (FAILLIMIT  >  0 && FAILURES   >=  FAILLIMIT)      //DG Change
	      );
    }
    //@}

    /**@name Helpers*/
    //@{
    void binds( ReversibleObj& x )
    {
      x.level = &level;
      x.store = &store;
    }


    inline void closeSearch()
    {
      ENDTIME = getRunTime() - STARTTIME;
      TOTTIME = (ENDTIME + STARTTIME - INITTIME);
      if((goal || FIND_ALL > 0) && SOLUTIONS) {
	if(status == UNSAT) { status = OPT; }
	else { status = SAT; }
      }
    }

    inline void learnFailure( Constraint *con )
    {
      ++FAILURES;
      int i=learners.size;
      while( i-- )
	learners[i]->notifyFailure( con );
    }

    inline void learnSuccess( )
    {
      int i=learners.size;
      while( i-- ) {
	learners[i]->notifySuccess( );
      }
    }


    //void setHeuristic(string hname, const int rdz);
    void setHeuristic(const char* var_name, const char* val_name, const int rdz);
    void setBranching(BuildObject **x, const int l, const char* val_name);

    void setRandomValueOrdering();

    void setLex();
    void setAntiLex();
    void setSplit();
    void setRandSplit();
    void setRandMinMax();

    void setLex(VarArray&);
    void setAntiLex(VarArray&);
    void setSplit(VarArray&);
    void setRandSplit(VarArray&);
    void setRandMinMax(VarArray&);

    void setLex(BuildObject **, const int);
    void setAntiLex(BuildObject **, const int);
    void setSplit(BuildObject **, const int);
    void setRandSplit(BuildObject **, const int);
    void setRandMinMax(BuildObject **, const int);

    /* flags:
     * std: x == ideal / x != ideal (min if ideal is not in the domain)
     */
    void setGuidedOrdering(VarArray&, int*, const char* flag="std");
    void setRandGuidedOrdering(VarArray&, int*, int*, int*);
    void setGuidedOrdering(BuildObject **x, const int l, int*, int pb=0);
    void setRandGuidedOrdering(BuildObject **x, const int l, int*, int*, int*);
    void setGuidedSplitOrdering(BuildObject **bvar, const int l, int* ideal, int* proba);

    void setLowerBounds(VarArray& x, int*);
    void setUpperBounds(VarArray& x, int*);
    void setLowerBounds(BuildObject **x, const int l, int*);
    void setUpperBounds(BuildObject **x, const int l, int*);


    /// Cutoff in number of backtracks
    void setNodeLimit     (unsigned long nl){ NDSLIMIT  = nl; LIMIT = 1; }
    /// Cutoff in number of backtracks
    void setBacktrackLimit(unsigned long nl){ BTSLIMIT  = nl; LIMIT = 1; }
    /// Cutoff in number of failures
    void setFailureLimit  (unsigned long nl){ FAILLIMIT = nl; LIMIT = 1; }
    /// Cutoff in cpu time (s)
    void setTimeLimit     (double tl)       { TIMELIMIT = tl; LIMIT = 1; }
    /// Print satisfiability, number of backtracks and cpu time
    void setVerbosity     (const int v)     { verbosity =  v; }
    /// Add a learner
    Weighter* setLearner( int wtype );
    /// Use domain splitting
    void setDomainSplitting() { domainSplitting = true; }
    /// Shuffle the array of variables
    void randomizeSequence();
    /// Reinitialise the array of variables
    void reorderSequence();
    /// Randomize the array of variables before restarts
    void setRandomized( const bool r=true ) { randomizedRestart = r; }
    /// Set the random seed
    void setRandomSeed( const unsigned int seed ) {

      usrand(seed);
    }
    /// Setup the nogoods on restart
    void setRestartNogood();
    //ConstraintGenNogoodBase* setRestartGenNogood();
    WeighterRestartGenNogood* setRestartGenNogood();
    //void addNogood();
    /// Setup the forget utility for the clause base
    void setForgetfulness( const double f );
    /// Print statistics corresponding to the bitset stats
    void printStatistics(std::ostream &log=std::cout, const int stats=ALL);
    unsigned long int getBacktracks() const ;
    unsigned long int getNodes     () const ;
    unsigned long int getFailures  () const ;
    unsigned long int getChecks    () const ;
    unsigned long int getPropags   () const ;
    double            getTime      () const ;
 
    void test_x60() {
      if(variables[59]->weight > 100) {
	variables[59]->print(std::cout);
	std::cout << std::endl;
      }
    }

    void printSequence()
    {
      int i=0;
      std::cout << "[ ";
      for(; i<length; ++i) {
	if(sequence+i == future)
	  std::cout << " | ";
	else if(sequence+i == empty)
	  std::cout << " ] ";
	sequence[i]->print(std::cout);
	std::cout << " ";
      }
      std::cout << std::endl;
    }
    void print(std::ostream& o) const;
    void printPython() const;
    void printXML(std::ostream& o) const;
    void printWeightProfile(std::ostream& o, int limit, int t) const;
    double getGiniVarCoef() const;
    double getGiniConCoef() const;
    double getHooverCoef() const;
    //@}

#ifdef _WEIGHT_STATS

    void print_search_stats(std::string& key);

#endif

  };

};


#endif // _SOLVER_H
