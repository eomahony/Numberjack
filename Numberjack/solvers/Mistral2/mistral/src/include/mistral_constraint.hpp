
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


/** \file mistral_constraint.hpp
    \brief Header for the Constraint library.
*/


#include <string>
#include <vector>

#include <mistral_variable.hpp>


#ifndef __CONSTRAINT_HPP
#define __CONSTRAINT_HPP



#define _ALLDIFF_WC
#define _ELT_WC
#define _CNE_WC

//#define _DIV_ARITY

//#define _DIV_WEIGHT

#define _PWBS_WC
#define _PWS_WC
#define _CIWBSI_WC
#define _CWBSI_WC
#define _PBS_WC
#define _CBSI_WC



/*
#define _PWS_WC_ALT
#define _CWBSI_WC_ALT
#define _CIWBSI_WC_ALT
#define _PWBS_WC_ALT
*/

#define FILTER1( var, method )				    \
  event_type[(var)] = scope[(var)].method ;		     \
  if(event_type[(var)] == FAIL_EVENT) wiped = FAILURE(var);		\
  else if(event_type[(var)] != NO_EVENT && !changes.contain(var)) changes.add(var); 


#define FILTER2( evt, var, method )		\
  evt = scope[(var)].method ; \
  if(evt == FAIL_EVENT) wiped = FAILURE(var); \
  else if(evt != NO_EVENT && !changes.contain(var)) { \
  changes.add(var); \
  event_type[(var)] |= evt; \
  }

#define FILTER4( var, method )				    \
  event_type[(var)] |= scope[(var)].method ;		     \
  if(FAILED(event_type[(var)])) wiped = FAILURE(var);			\
  else if(event_type[(var)] != NO_EVENT && !changes.contain(var)) changes.add(var); 

//returning the corresponding index ogf the variable from right to left
#define INVERSE(size, index) ((size - index - 1))


#define FILTER3( var, method ) \
  Event evt = scope[(var)].method ; \
  if(FAILED(evt)) wiped = FAILURE(var);	\
  else if(evt != NO_EVENT) { \
    if(changes.contain(var)) { \
      event_type[(var)] |= evt; \
    } else { \
      event_type[(var)] = evt; \
      changes.add(var); \
    } \
  } \


//#define  _DEBUG_TABLE (id == 6)
//((scope[0].id() == 9 || scope[1].id() == 9 || scope[2].id() == 9))
//#define _OS_PRUNING true


int integer_div_up(const int x, const int y);
int integer_div_lo(const int x, const int y);



/**
 *Constraint:
 -propagator         (pointer to the constraint)
 -postponed          (should it be propagated on changed and/or postponed)


 *ConstraintImplementation: 
 -scope              (the variables)
 -id                 (an id)
 -trigger / index    (post/relax operations)

 +virtual propagate()
 +virtual propagate(const int changedIdx, const Event evt)
 +virtual check()

 *BinaryConstraintImplementation:

 *TernaryConstraintImplementation:
 
 *GlobalConstraintImplementation:
 */

/**
   Constraint posting/relaxing
   ===========================

   * For each variable in their scope, a constraint has a pointer to the list they belong to for that variable [on]
   - If the constraint should not be called on any modification of this variable, the pointer is null
   - If the constraint should be called on value events, it points to the value list of the given variable
   - If the constraint should be called on range events, it points to the range list of the given variable
   - If the constraint should be called on domain events, it points to the domain list of the given variable

   * Again for each variable in the scope, it keeps its current index in the given list [index]
   - If the constraint is not currently in the list of var $i, the value of index[$i] should be -1
   - If the constraint is currently posted in the list of var $i, the value of index[$i] should be its rank in the list

   * For each variable, a version of a MistralPointer to itself is kept. [self]
   Moreover, pointer info is given the rank of that variable in the scope. 
   This is this pointer which is added to the list of triggers of a variables (via [on])
   This is important when propagating: when scanning the list of triggers for a given variable, 
   it indicates the rank of the variable in the constraint's scope, allowing to use the version of
   propagate() informed about the caller.

   * A constraint has flag to indicate whether it achieves at least a level of propagation equal to NFC1 [enforce_nfc1]
   - If it is set to true, then the constraint will be automatically relaxed from its trigger list in the last
   active variable when there is only one of them
   - If it is set to false, it is never relaxed

   * The set of "active" variables is kept in a specific, however non-virtual attribute "active"


   Automatic relax:
   ================
   
   When the solver process a value event on a variable, it notifies it to all its constraints.

   [notify_assignment]: The variable is removed from the set [active]. If [active] is reduced to a signleton,
   the constraint is relaxed from the last variable in [active].

   [relax_from($i)]: if the constraint is currently posted on the $ith variable, then:
   the constraint is relaxed from that variable.


   Backtracks:
   ===========

   The method [save(Constraint c)] indicates that the constraint should be saved.
   That is, when backtracking upon this level, the restore() method of $c will be called.
   Some information about what should be undone is passed trough the [data] field of $c.

   1/ First, it contains the index of the variable
   2/ Second, it can be complemented with the flags [ACTIVITY], [RELAXED] or [POSTED].

   FixedArityConstraints:

   * notify_assignment() changes the value of the set [active], hence it turn on the flag [ACTIVITY]
   if activity becomes a singleton, the constraint is relaxed from the last remaining variable and the flag [RELAXED] is set before saving.

   * When relaxing the flag [RELAXED] is set
   
   * When posting the flag [POSTED] is set

   * Restore: the constrint is posted/relaxed from the corresponding variable according to the flag
       If the flag [ACTIVITY] was set, then
              a/ For Binary constraints, activity is set back to {0,1}
	      b/ For Ternary constraints, activity is reversed using a specific mechanism

   [TODO:] we could have a different save operation for acticity and post/relax
   The advantage is that the activity set could be reversed by simply adding 
   the 



 */


namespace Mistral {

  /********************************************
   * ConstraintWrapper 
   ********************************************/
  /*! \class ConstraintWrapper
    \brief Wrapper for constraints. 
    
    Holds a pointer to a constraint and the index 
    (in the scope) of the variable that owns the list 
  */  

  
  class Variable;
  class ConstraintImplementation : public Explanation {
    
  public:
    
    /*!@name Parameters*/
    //@{
    /// Pointer to its solver
    Environment *solver;
    /// An unique id
    int id;
  
    // used to post/relax the constraints on the right triggers
    // the list on wich it should be added [among the 3 triggers of the variable]
    Vector< Trigger* > on;

    // the position of that constraint in the list
    int            *index;
    
    // what to put on the list: a pointer to the constraint and the rank of the variable 
    // (so that it knows who triggered it)
    Constraint      *self;

    // the variable to which the trigger applies to 
    Vector< Variable > _scope;

    
    //int arity;
    //int num_triggers;
    unsigned int type;
    bool enforce_nfc1;
    //@}
    

    /*!@name Constructors*/
    //@{
    /// The _scope is build by copying an array "scp" containing "l" variables
    ConstraintImplementation();
    //ConstraintImplementation(const int a);
    virtual Constraint clone() = 0;
    virtual ~ConstraintImplementation();


    // virtual Explanation::iterator begin(Atom a) { return NULL; }
    // virtual Explanation::iterator end  (Atom a) { return NULL; } 
    // TODO: implement is for all constraints, and leave that one void
    virtual iterator get_reason_for(const Atom a, const int lvl, iterator& end) {
      return (end = NULL);
    }
    
    virtual void initialise() { type = get_type(); }
    virtual void initialise_vars(Solver*) = 0;

    virtual int get_type() = 0;
    virtual int postponed() = 0;
    virtual int idempotent() = 0;
    //virtual int enforce_nfc1() = 0;

    virtual void mark_domain() {}    

    virtual void desactivate(const int var) = 0;

    virtual void print_active() = 0;

    //virtual void remove_duplicates();

    bool is_active() {
      bool active = false;
      for(unsigned int i=0; !active && i<on.size; ++i) {
	active = index[i]>=0;
      }
      return active;
    }

   void un_post_from(const int var) {
     //std::cout << "relax [" << id << "] from " << scope[var] << ": " << on[var] << " -> " ;      
      
      on[var]->relax(index[var]);
      index[var] = -1;    
      
      //std::cout << on[var]  << std::endl;
    }
    
    void un_relax_from(const int var) {
      index[var] = on[var]->post(self[var]);
    }

    // does not post a constraint that was posted on the same variable recently
    void check_and_un_relax_from(const int var) {
      index[var] = on[var]->check_and_post(self[var]);
    }


    void trigger_on(const int t, Variable x) ;//{
    //   trigger.add(solver->constraint_trigger_graph[x.id()][t]);
    // }

    bool is_triggered_on(const int i, const int t);

    void initial_post(Solver *s);
   
    int get_trigger_type(const int i) ;
    void set_scope(const int i, Variable x);

    /*!@name Propagators*/
    //@{
    /*!
     * This methods returns 0 if the constraint is satisfied
     * by the combination and any positive value otherwise.  
     */
    virtual int check(const int*) const = 0;
    /*!
     *  This method is called when the domain of at least one   
     *  variable in _scope has been modified. The list of 
     *  changed variables is in 'modified', and the type
     *  of event in 'event_type'.
     *  Some domains in _scope may be reduced, and NULL is
     *  returned on success (there is still at least one consistent
     *  assignment) or a ptr to the wiped-out variable otherwise. 
     */
    //virtual PropagationOutcome propagate_and_explain() { return CONSISTENT; }
    virtual PropagationOutcome propagate() = 0; // { return NULL; }
    virtual PropagationOutcome checker_propagate() = 0;
    virtual PropagationOutcome bound_checker_propagate() = 0;
    virtual PropagationOutcome propagate(const int changed_idx, 
					 const Event evt) = 0;
    virtual PropagationOutcome checker_propagate(const int changed_idx, 
					 const Event evt) = 0;
    virtual PropagationOutcome bound_checker_propagate(const int changed_idx, 
					 const Event evt) = 0;
    //  { return CONSISTENT; }
    virtual bool absorb_negation(const int var) { return false; }
    virtual Constraint get_negation(const int var, Variable x) { return Constraint(); }
    virtual bool rewritable() { return false; }
    virtual bool simple_rewritable() { return false; }
    virtual bool explained() { return false; }
    virtual RewritingOutcome rewrite() { return NO_EVENT; }
    virtual void consolidate() = 0; 
    virtual void consolidate_var(const int idx) = 0; 

    virtual void initialise_activity(double *lact, double *vact, double norm) = 0;
    //virtual int get_backtrack_level();// {return solver->level-1;}
    //virtual Decision get_decision();// { return solver->decisions.back(0); }
    //@}


    inline Solver* get_solver() { return (Solver*)solver; }

    /*!@name Miscellanous methods*/
    //@{
    /// Print the constraint
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "c"; }
    virtual bool is_clause() {return false;}
    //virtual bool verify_state() = 0;
    //@}
  };





  class ConstraintTriggerArray {
    
  public:
    
    Trigger on[3];
    //Vector< ConstraintWrapper > range_trigger;
    //Vector< ConstraintWrapper > domain_trigger;

    ConstraintTriggerArray() ;
    ConstraintTriggerArray(const int size) ;
    void initialise(const int size) ;

    virtual ~ConstraintTriggerArray() ;

    inline int size() { return (on[0].size+on[1].size+on[2].size); }

    // each constraint keeps its index in the array it appears in
    // to remove: trigger.remove(index)
    // to add: index[i] = trigger[i].size; trigger.add(self[i]);

    std::ostream& display(std::ostream& os) const ;


    // class friend Iterator {
    
    // public :

    //   int T;
    //   int i;

    // };

  };


  /**
     
     POST/RELAX

     relax -> remove the constraint from all the active variables' lists
     relax_from(var) -> remove the constraint from var
     -> flag == active

     un_relax -> add the constraint to all the active vars' lists
     


   */

  template< int ARITY >
  class ValTuple {
  public:

    ValTuple() { std::fill(data, data+ARITY, NOVAL); }

    int data[ARITY];
    
    inline int& operator[](const int i) { return data[i]; }
    inline int operator[](const int i) const { return data[i]; }
  };



  template< int ARITY >
  class FixedArityConstraint : public ConstraintImplementation {

  public:

    int solution[ARITY];
    ValTuple<ARITY> *support[ARITY];
    Variable scope[ARITY];
    int active;


    FixedArityConstraint() : ConstraintImplementation() { 
      for(int i=0; i<ARITY; ++i) support[i] = NULL;
      //std::cout << "there" << std::endl;
      active = (1 << ARITY)-1; 
      //init_type(); //type = get_type();
    }

    FixedArityConstraint(Variable x, Variable y) : ConstraintImplementation() { 
      for(int i=0; i<ARITY; ++i) support[i] = NULL;
      //std::cout << "there" << std::endl;
      active = (1 << ARITY)-1; 
      //init_type(); //type = get_type();
      scope[0] = x;
      scope[1] = y;
    }

    FixedArityConstraint(Variable x, Variable y, Variable z) : ConstraintImplementation() { 
      for(int i=0; i<ARITY; ++i) support[i] = NULL;
      //std::cout << "there" << std::endl;
      active = (1 << ARITY)-1; 
      //init_type(); //type = get_type();
      scope[0] = x;
      scope[1] = y;
      scope[2] = z;
    }

    FixedArityConstraint(Vector< Variable >& scp) : ConstraintImplementation() { 
      for(int i=0; i<ARITY; ++i) support[i] = NULL;
      //std::cout << "there" << std::endl;
      active = (1 << ARITY)-1; 

      //std::cout << "there " << active << std::endl;

      //init_type(); //type = get_type();
      for(int i=0; i<ARITY; ++i) {
	scope[i] = scp[i];
      }
    }

    FixedArityConstraint(std::vector< Variable >& scp) : ConstraintImplementation() { 
      for(int i=0; i<ARITY; ++i) support[i] = NULL;
      //std::cout << "there" << std::endl;
      active = (1 << ARITY)-1; 
      //init_type(); //type = get_type();
      for(int i=0; i<ARITY; ++i) {
	scope[i] = scp[i];
      }
    }


    virtual void initialise_vars(Solver *s) {
      for(unsigned int i=0; i<ARITY; ++i) {
	scope[i].initialise(s, false);
      }
    }

    virtual void initialise_activity(double *lact, double *vact, double norm) {
      double w = norm/ARITY;
      int idx, i = ARITY;
      while(i--) {
	idx = scope[i].id();
	vact[idx] += w;
	lact[2*idx] += w/2;
	lact[2*idx+1] += w/2;
      }
    }


    //int init_type() { std::cout << "here" << std::endl; type=get_type(); }
    virtual int get_type() {return 0;}
    virtual int postponed() {return 0;}
    //virtual void initialise() {type = get_type();}
    //virtual int idempotent() {return 0;}

    virtual void desactivate(const int var) { 


// #ifdef _DEBUG_RELAX
//       std::cout // << "[" << std::setw(4) << id << "](" << name() << "): desactivate "
// 	<< "-"
// 	<< scope[var] << " " ; //<< std::endl;
// #endif

      int var_elt = (1 << var);
      if(active&var_elt) active ^= var_elt;
    }

    void initialise_supports() {
      int vari, vali, vnext, min_vari;
      for(vari=0; vari<ARITY; ++vari) {
	min_vari = scope[vari].get_min();
	
	support[vari] = new ValTuple<ARITY>[scope[vari].get_max() - min_vari + 1];
	support[vari] -= min_vari;
	
	vnext = scope[vari].get_min();
	do {
	  vali = vnext;
	  //support[vari][vali] = new int[ARITY];
	  support[vari][vali][vari] = vali;
	  vnext = scope[vari].next(vali);

	} while( vali < vnext );
      }
    }

    virtual void consolidate() {
      
      // std::cout << "consolidate " ;
      // display(std::cout);
      // std::cout << std::endl;

      for(unsigned int i=0; i<_scope.size; ++i) {
	_scope[i] = _scope[i].get_var();
      }
      
      for(unsigned int i=0; i<ARITY; ++i) {
	scope[i] = scope[i].get_var();
      }

      //std::cout << _scope << std::endl << scope[0] << " " << scope[1] << std::endl;

    }

   virtual void consolidate_var( const int idx ) {
      
     // std::cout << "consolidate " ;
     // display(std::cout);
     // std::cout << std::endl;
     
     _scope[idx] = _scope[idx].get_var();
     scope[idx] = scope[idx].get_var();
     
     //std::cout << _scope << std::endl << scope[0] << " " << scope[1] << std::endl;
   }


    void check_active() {
      for(int i=on.size; i--;) {
	if(index[i]>=0) {
	  if(active & (1 << i)) {
	    if(scope[i].is_ground()) {
	      std::cout << "[" << std::setw(4) << id << "] " ;
	      display(std::cout);
	      std::cout << std::endl;
	      std::cout << "Warning: " << scope[i] << " = " << scope[i].get_domain()
			<< " is ground and active!!" << std::endl;
	      print_active();
	      std::cout << " :: ";
	      for(unsigned int j=0; j<on.size; ++j)
		std::cout << scope[j] << " in " << scope[j].get_domain() << " ";
	      std::cout << " (exit on check_active())" << std::endl;
	      //exit(1);
	    }
	  } else {
	    if(!(scope[i].is_ground())) {
	      std::cout << "[" << std::setw(4) << id << "] " ;
	      display(std::cout);
	      std::cout << std::endl;
	      std::cout << "Warning: " << scope[i] << " in " << scope[i].get_domain()
			<< " is not ground and not active!!" << std::endl;
	      print_active();
	      std::cout << " :: ";
	      for(unsigned int j=0; j<on.size; ++j)
		std::cout << scope[j] << " in " << scope[j].get_domain() << " ";
	      std::cout << " (exit on check_active())" << std::endl;
	      //exit(1);
	    }
	  }
	}
      }
    }
      
    virtual void print_active() {
      

      int k=0, x = (active & (7<<k)) >> k;
      
      std::cout << active << " -> " << x << " ";
      
      while(x) {
	
	if(x==7) std::cout << "{" << scope[0] << "," << scope[1] << "," << scope[2] << "}" ;
	else if(x==6) std::cout << "{"  << scope[1] << "," << scope[2] << "}" ;
	else if(x==5) std::cout << "{" << scope[0] << "," << scope[2] << "}" ;
	else if(x==3) std::cout << "{" << scope[0] << "," << scope[1] << "}" ;
	else if(x==4) std::cout << "{" << scope[2] << "}" ;
	else if(x==2) std::cout << "{" << scope[1] << "}" ;
	else if(x==1) std::cout << "{" << scope[0] << "}" ;
	else std::cout << "{}" ;
	
	k+=3;
	x = (active & (7<<k)) >> k;

      }
    }
    

    void post_on(const int var) {

      // std::cout << "po post " ;
      // display(std::cout);
      // std::cout << " on " << scope[var] << std::endl;

//       if(scope[var].is_ground()) {
// 	active^=(1<<var);
// 	index[var] = -1;
//       } else {

// #ifdef _DEBUG_RELAX
//  if(_DEBUG_RELAX) {
// 	std::cout << "[" << std::setw(4) << id << "](" << name() << "): post on " << scope[var] << std::endl;
// #endif
// 	//std::cout << "(yep)" << std::endl;
	
// 	Constraint c = self[var];
// 	c.data |= POSTED;
// 	solver->save( c );
// 	un_relax_from(var);
//       }



      if(scope[var].is_ground()) {
	active^=(1<<var);
	index[var] = -1;
      } else {

#ifdef _DEBUG_RELAX
  if(_DEBUG_RELAX) {
	std::cout << "[" << std::setw(4) << id << "](" << name() << "): post on " << scope[var] << std::endl;
  }
#endif
	//std::cout << "(yep)" << std::endl;
	
	Constraint c = self[var];
	c.data |= POSTED;
	solver->save( c );
	un_relax_from(var);
      }
    }
  
    void post() {

      // std::cout << "p post " ;
      //  display(std::cout);
      //  std::cout << std::endl;

      active = (1 << ARITY)-1;
      int nb_actives = on.size;
      int i = nb_actives;

      while(i--) {
	index[i] = -1;
	if(scope[i].is_ground()) {
	  active^=(1<<i);
	  --nb_actives;
	}
      }

      if(!enforce_nfc1 || nb_actives>1) {
	Constraint c;
	
	for(i=on.size; i--;) {
	  //if(!scope[i].is_ground()) {
	  if(on[i]) {

	    //std::cout << "c *** " << scope[i] << " in " << scope[i].get_domain() << std::endl;

	    c = self[i];
	    c.data |= POSTED;
	    solver->save( c );
	    un_relax_from(i);	  
	  }
	}
      }
    }


    void relax() {
      
      // std::cout << "active variables: ";
      // print_bitset(active, 0, std::cout);
      // std::cout << "( ";
      // for(int i=on.size; i--;) {
      // 	if(active & (1 << i)) {
      // 	  std::cout << scope[i] << " ";
      // 	}
      // }
      // std::cout << ")" << std::endl;

      for(int i=on.size; i--;) {
	if(active & (1 << i)) relax_from(i);
      }    
    }

    void relax_from(const int var) {

      //if(active & (1 << var)) {
      if(index[var] >= 0) {

#ifdef _DEBUG_RELAX
  if(_DEBUG_RELAX) {
      std::cout << "[" << std::setw(4) << id << "](" << name() << "): f relax from " << scope[var] << std::endl;
  }
#endif

     Constraint c = self[var];
      c.data |= RELAXED;
      solver->save( c ); //self[var]|RELAXED );
	un_post_from(var);
      }
    }    

    std::ostream& display(std::ostream& os) const {
      os << name() << "(" << scope[0] << ", " << scope[1];
      if(ARITY > 2)
	os << ", " << scope[2];
      os << ")";
      return os;
    }
  
  };


  class BinaryConstraint : public FixedArityConstraint<2> {
    
  public:

    BinaryConstraint() : FixedArityConstraint<2>() {}
    BinaryConstraint(Variable x, Variable y)  : FixedArityConstraint<2>(x,y) {}
    BinaryConstraint(Vector< Variable >& scp) : FixedArityConstraint<2>(scp) {}
    BinaryConstraint(std::vector< Variable >& scp) : FixedArityConstraint<2>(scp) {}
    virtual int get_type() { return BINARY|(IDEMPOTENT*idempotent()); }

    virtual PropagationOutcome propagate(); // { return NULL; }
    virtual PropagationOutcome bound_propagate(); // { return NULL; }
    virtual PropagationOutcome checker_propagate() { return BinaryConstraint::propagate(); }
    virtual PropagationOutcome bound_checker_propagate() { return BinaryConstraint::bound_propagate(); }
    virtual PropagationOutcome propagate(const int changed_idx, 
					 const Event evt);
    virtual PropagationOutcome bound_propagate(const int changed_idx, 
					 const Event evt);
    virtual PropagationOutcome checker_propagate(const int changed_idx, 
						 const Event evt) { return BinaryConstraint::propagate(changed_idx, evt); }
    virtual PropagationOutcome bound_checker_propagate(const int changed_idx, 
						       const Event evt) { return BinaryConstraint::bound_propagate(changed_idx, evt); }
    bool find_support(const int revise_idx, const int vli);
    bool find_bound_support(const int revise_idx, const int vli);


   inline void notify_assignment(const int var) {

#ifdef _DEBUG_RELAX
  if(_DEBUG_RELAX) {
     std::cout << "[" << std::setw(4) << id << "](" << name() << "): notify assignment of " << scope[var] << std::endl;
  }
#endif
     
     
     // if(!(active & (1 << var))) {
       
     //   std::cout << "SHOULD NOT NOTIFY TWICE (" << scope[var] << " in " << scope[var].get_domain() << ")" << std::endl;
       
     // } else {

     if(active & (1 << var)) {
       
       active ^= (1 << var);
     
       int last = (active&7)/2;
       Constraint c = self[last];
       c.data |= ACTIVITY;

       // print_active();
       // std::cout << " " << last << std::endl;

       if(enforce_nfc1 && size_byte[active] == 1) {
	 
	 if(index[last] >= 0) {
	   
#ifdef _DEBUG_RELAX
  if(_DEBUG_RELAX) {
	   std::cout << "[" << std::setw(4) << id << "](" << name() << "): f relax from " << scope[last] << std::endl;
  }
#endif
	   c.data |= RELAXED;
	   un_post_from(last);
	 }
       }
       
       solver->save( c ); //self[var]|RELAXED );
     }
   }

    void restore(const int rtype) {

#ifdef _DEBUG_BACKTRACK
      if(_DEBUG_BACKTRACK) {
      std::cout << "c ";
      int lvl=solver->level;
      while(--lvl>=0) std::cout << " ";
      std::cout << "[" << std::setw(4) << id << "](" << name() << "): restore" << std::endl;
  }
#endif
     
      int var = rtype&CTYPE;

      // if(var<0 || var>1) {
      //   std::cout << var << std::endl;
      //   exit(1);
      // }

      // #ifdef _DEBUG_RELAX 
      //      std::cout << "RESTORE [" << id << "] " ;
      //      display(std::cout);
      //      std::cout << " (" << (rtype&2) << " / " << (rtype&CTYPE) << ")" << std::endl;
      // #endif

      //     if(index[var]<0) {

      // #ifdef _DEBUG_RELAX
      //        std::cout << "[" << std::setw(4) << id << "](" << name() << "): repost on " << scope[var] << std::endl;
      // #endif

      //        un_relax_from(var);
      //      } else {

      // #ifdef _DEBUG_RELAX
      //        std::cout << "[" << std::setw(4) << id << "](" << name() << "): relax from " << scope[var] << std::endl;
      // #endif

      //        un_post_from(var);
      //      }
      //      active = 3;

      // #ifdef _DEBUG_RELAX
      //      std::cout << "[" << std::setw(4) << id << "](" << name() << "): reset active: " ;
      //      print_active();
      //      std::cout << std::endl;
      // #endif


      if(rtype&RELAXED) {

      	// if(index[var]>=0) {

	//   std::cout << "SHOULD NOT REPOST TWICE!" << std::endl;
	//   exit(1);

	// }

	if(index[var] < 0) {

#ifdef _DEBUG_BACKTRACK
	  if(_DEBUG_BACKTRACK) {
	    std::cout << "c ";
	    int lvl=solver->level;
	    while(--lvl>=0) std::cout << " ";
	    std::cout << "[" << std::setw(4) << id << "](" << name() << "): repost on " << scope[var] << std::endl;
	  }
#endif
	  
	  un_relax_from(var);
	} 
      }


      if(rtype&POSTED) {

	// if(index[var]<0) {
       
	//   std::cout << "SHOULD NOT RELAX TWICE!" << std::endl;
	//   exit(1);

	// }

	if(index[var] >= 0) {

#ifdef _DEBUG_BACKTRACK
      if(_DEBUG_BACKTRACK) {
	  std::cout << "c ";
	  int lvl=solver->level;
	  while(--lvl>=0) std::cout << " ";
	  std::cout << "[" << std::setw(4) << id << "](" << name() << "): relax from " << scope[var] << std::endl;
  }
#endif

	  un_post_from(var);
	}
      }
    

      if(rtype&ACTIVITY) active = 3;

#ifdef _DEBUG_BACKTRACK
      if(_DEBUG_BACKTRACK) {
      std::cout << "c ";
      int lvl=solver->level;
      while(--lvl>=0) std::cout << " ";
      std::cout << "[" << std::setw(4) << id << "](" << name() << "): reset active: " ;
      print_active();
      std::cout << std::endl;
  }
#endif



      // if(scope[0].is_ground() || scope[1].is_ground()) {

      //   std::cout << "ERROR RESTORE BINARY CONSTRAINT" << std::endl;
       
      //   std::cout << scope[0].get_domain() << " " << scope[1].get_domain() << std::endl;

      // }

    }  
    
    void trigger();
  };


  class TernaryConstraint : public FixedArityConstraint<3> {

  public:

    int lvl;

    TernaryConstraint() : FixedArityConstraint<3>() { lvl = 3; }
    TernaryConstraint(Variable x, Variable y, Variable z) : FixedArityConstraint<3>(x,y,z) { lvl = 3; }
    TernaryConstraint(Vector< Variable >& scp) : FixedArityConstraint<3>(scp) { lvl = 3; }
    TernaryConstraint(std::vector< Variable >& scp) : FixedArityConstraint<3>(scp) { lvl = 3; }
    virtual int get_type() {return TERNARY|(IDEMPOTENT*idempotent());}

    virtual PropagationOutcome checker_propagate() 
    { return TernaryConstraint::propagate(); }
    virtual PropagationOutcome bound_checker_propagate() { 
      return TernaryConstraint::bound_propagate(); 
    }
    virtual PropagationOutcome propagate(); // { return NULL; }
    virtual PropagationOutcome bound_propagate(); // { return NULL; }
    virtual PropagationOutcome propagate(const int changed_idx, 
					 const Event evt);
    virtual PropagationOutcome bound_propagate(const int changed_idx, 
					 const Event evt);
    virtual PropagationOutcome checker_propagate(const int changed_idx, 
						 const Event evt) 
    { return TernaryConstraint::propagate(changed_idx, evt); }
    virtual PropagationOutcome bound_checker_propagate(const int changed_idx, 
						       const Event evt) { 

      //std::cout << "constraint.cpp: bound checker propagate(evt)" << std::endl; 

      return TernaryConstraint::bound_propagate(changed_idx, evt); }
    bool find_support(const int revise_idx, const int vli);
    bool find_bound_support(const int revise_idx, const int vli);


    /*
    // the 3 first bits stand for the current set of active vars, 
    // the 3 next bits stand for the intermediate state 
    // lvl stands for the level at wich the intermediate


    -------- FIRST POSSIBILITY -----------
    // when assigning a new variable, 
    -- if lvl = -1 :  do the remove and change lvl (save)
    -- if lvl = s->level : do the remove
    -- if -1 < lvl < s->level : store the current state, do the remove (save)

    // when backtracking
    -- if s->level > lvl : copy the stored state
    -- otherwise : set to {0,1,2}

              -1|a    b
    {0,1,2} | {1,2} | {2}

             -1|a|a
    {0,1,2} | {1}  

    ------- SECOND POSSIBILITY -----------
    lvl store the current position where values should eb stored
    
    Whenever a variable is assigned, store the current state into the next slot and save.

    {1,2}{0,1,2}

    */

   inline void notify_assignment(const int var) {

#ifdef _DEBUG_RELAX
  if(_DEBUG_RELAX) {
     std::cout << "[" << std::setw(4) << id << "](" << name() << "): notify assignment of " << scope[var] << std::endl;
  }
#endif
     
     
     // if(!(active & (1 << var))) {
       
     //   std::cout << "SHOULD NOT NOTIFY TWICE (" << scope[var] << " in " << scope[var].get_domain() << ")" << std::endl;
       
     // } else {

     if(active & (1 << var)) {
       
       int elt = (1 << var);
       int tmp = active&7;
       tmp ^= elt;
       

       //std::cerr << active << " "<< tmp << " "<< last << " "<< elt << std::endl; 

       active <<= 3;
       active |= tmp;


       Constraint c;
       
       //solver->save(Constraint(this, type|ACTIVITY));
       
       // print_active();
       // std::cout << " " << last << std::endl;

       if(enforce_nfc1 && size_byte[active] == 1) {
	 
	 int last = tmp/2;
	 c = self[last];

	 if(index[last] >= 0) {
	   
#ifdef _DEBUG_RELAX
  if(_DEBUG_RELAX) {
	   std::cout << "[" << std::setw(4) << id << "](" << name() << "): f relax from " << scope[last] << std::endl;
  }
#endif
	   c.data |= (RELAXED|ACTIVITY);
	   un_post_from(last);
	 }
       } else {
	 c = self[0];
	 c.data |= ACTIVITY;
       }
       
       solver->save( c ); //self[var]|RELAXED );
     }
   }

    bool assign(const int var) {

      // //if(id==36) {
      // std::cout << std::endl;
      // print_active();
      // std::cout << " Assign " << scope[var] << " " ;
      // //}

#ifdef _DEBUG_RELAX
  if(_DEBUG_RELAX) {
     std::cout << "[" << std::setw(4) << id << "](" << name() << "): notify assignment of " << scope[var] << std::endl;
  }
#endif
	
      bool ret_val = false;
      int elt = (1 << var);
      if(active & elt){


#ifdef _DEBUG_RELAX
  if(_DEBUG_RELAX) {
     std::cout << "[" << std::setw(4) << id << "](" << name() << "): -> deseactivate" << std::endl;
  }
#endif

	int tmp = active&7;
	active <<= 3;
	active |= tmp;
	active ^= elt;

	if(active&448) ret_val = true;
	solver->save(Constraint(this, type|ACTIVITY));
	//solver->save(self[0]|ACTIVITY));
      
      // if(active&64) ret_val = true;
	// else solver->save(Constraint(this, type|ACTIVITY));
      }

      // //if(id==36) {
      // print_active();
      // std::cout << std::endl;
      // //}

      return ret_val;
    }


    void restore(const int rtype) {
      
#ifdef _DEBUG_BACKTRACK
      if(_DEBUG_BACKTRACK) {
            std::cout << "c ";
      int lvl=solver->level;
      while(--lvl>=0) std::cout << " ";
std::cout << "[" << std::setw(4) << id << "](" << name() << "): restore" << std::endl;
  }
#endif
      
      // //if(id==36) {
      // //std::cout << std::endl;
      
      // // std::cout << scope[0] << " " << on[0] << std::endl;
      // // std::cout << scope[1] << " " << on[1] << std::endl;
      // // std::cout << scope[2] << " " << on[2] << std::endl;
      // print_active();
      // std::cout << " restore " ;
      // //}

      if(rtype&ACTIVITY) active >>= 3;	
      
      if(rtype&RELAXED) {
	int var = rtype&CTYPE;
	
	// if(index[var]>=0) {
	  
	//   std::cout << "SHOULD NOT REPOST TWICE!" << std::endl;
	//   exit(1);

	// }
	// std::cout << on[var] << " repost " ;
	// display(std::cout);
	// //std::cout << std::endl; 
	
	if(index[var] < 0) {

#ifdef _DEBUG_BACKTRACK
      if(_DEBUG_BACKTRACK) {
	  std::cout << "c ";
	  lvl=solver->level;
	  while(--lvl>=0) std::cout << " ";
	  std::cout << "[" << std::setw(4) << id << "](" << name() << "): repost on " << scope[var] << std::endl;
  }
#endif
	
	  un_relax_from(var);
	}
      }
      
      
      if(rtype&POSTED) {
	
	int var = rtype&CTYPE;
	
	// if(index[var]<0) {
	  
	//   std::cout << "SHOULD NOT RELAX TWICE!" << std::endl;
	//   exit(1);

	// }
	
	if(index[var] >= 0) {

#ifdef _DEBUG_BACKTRACK
      if(_DEBUG_BACKTRACK) {
	  std::cout << "c ";
	  int lvl=solver->level;
	  while(--lvl>=0) std::cout << " ";
	  std::cout << "[" << std::setw(4) << id << "](" << name() << "): relax from " << scope[var] << std::endl;
  }
#endif
	  
	  un_post_from(var);
	}
      }
      
      //       else {
      //       int var = rtype&CTYPE;
      
  //       //std::cout << var << " " << index[var] << " ";
  
  
  //       if(index[var]<0) {
  
  // 	// std::cout << on[var] << " repost " ;
  // 	// display(std::cout);
  // 	// //std::cout << std::endl; 
  
  // #ifdef _DEBUG_RELAX
  //        std::cout << "[" << std::setw(4) << id << "](" << name() << "): repost on " << scope[var] << std::endl;
  // #endif
  
  // 	un_relax_from(var);
  
  // 	// std::cout << " " << on[var] << " ";

  //       } else {
  
  
  // #ifdef _DEBUG_RELAX
  //        std::cout << "[" << std::setw(4) << id << "](" << name() << "): relax from " << scope[var] << std::endl;
  // #endif
  
  // 	un_post_from(var);
  //       }
  //       }
  // //if(id==36) {
  // print_active();
  // std::cout << std::endl ;
  
  // // std::cout << scope[0] << " " << on[0] << std::endl;
  // // std::cout << scope[1] << " " << on[1] << std::endl;
  // // std::cout << scope[2] << " " << on[2] << std::endl;
  // std::cout << std::endl;
  
  // //}
  
#ifdef _DEBUG_BACKTRACK
      if(_DEBUG_BACKTRACK) {
      std::cout << "c ";
      lvl=solver->level;
      while(--lvl>=0) std::cout << " ";
  std::cout << "[" << std::setw(4) << id << "](" << name() << "): reset active: " ;
  print_active();
  std::cout << std::endl;
  }
#endif
   
}

    void update(const int changed_idx, const Event evt) {

      if(ASSIGNED(evt) && assign(changed_idx)) {

	
#ifdef _DEBUG_RELAX
  if(_DEBUG_RELAX) {
     std::cout << "[" << std::setw(4) << id << "](" << name() << "): ff notify assignment of " << scope[changed_idx] << std::endl;
  }
#endif

	// std::cout << " (" << ((active&7)/2) << ") ";
	// std::cout.flush();
	// std::cout << on[((active&7)/2)] << " relax " ;
	// display(std::cout) ;
	// std::cout.flush();


	relax_from((active&7)/2);

	// std::cout << " from " << scope[((active&7)/2)] 
	// 	  << on[((active&7)/2)] << std::endl;

      }
    }
    
    void trigger();
  };


  // class Contradiction : public ConstraintImplementation {

  // public:
  
  //   Contradiction() : ConstraintImplementation() {}
  //   //void initialise();
  //   virtual int get_type() {return 0;}
  //   virtual int postponed() {return 0;}
  //   virtual int idempotent() { return 1;}
  //   virtual ~Contradiction() {}

  //   virtual Constraint clone() { return Constraint(new Contradiction()); }
  //   virtual void initialise_vars(Solver *s) {}

  //   virtual void desactivate(const int var) {}
  //   virtual void consolidate() {}
  //   virtual void consolidate_var(const int idx) {}

  //   virtual int check(const int* sol) const { return 1; }

  //   virtual PropagationOutcome checker_propagate() { return FAILURE(0); }
  //   virtual PropagationOutcome bound_checker_propagate() { return FAILURE(0); }
  //   virtual PropagationOutcome propagate() { return FAILURE(0); }
  //   virtual PropagationOutcome propagate(const int changed_idx, const Event evt) { return FAILURE(0); }
  //   virtual PropagationOutcome bound_propagate() { return FAILURE(0); }
  //   virtual PropagationOutcome bound_propagate(const int changed_idx, const Event evt) { return FAILURE(0); }
  //   virtual PropagationOutcome checker_propagate(const int changed_idx, const Event evt)  { return FAILURE(0); }
  //   virtual PropagationOutcome bound_checker_propagate(const int changed_idx, const Event evt) { return FAILURE(0); }



  //   std::ostream& display(std::ostream& os) const { os << "False"; return os; }  
  // };




  class GlobalConstraint : public ConstraintImplementation {

  public:

    Vector< Variable > scope;
    /// We use two lists so that the active constraint can add events to its own 
    /// list (events) without changing the one used during propagation (changes)
    IntStack changes; // this is the list that is accessible from a propagator
    IntStack events; // this is the list that collects the events
    /// The type of event for each modified variable
    Event *event_type;
    /// Set of non-ground variables
    ReversibleSet active;

    int priority;

    ////
    int   *solution;
    int ***supports;
    ////

  
    GlobalConstraint() : ConstraintImplementation() {}
    GlobalConstraint(Vector< Variable > scp);
    GlobalConstraint(std::vector< Variable > scp);
    GlobalConstraint(Variable* scp, const int n);
    void initialise();
    virtual int get_type() {return (PUSHED*pushed())|(POSTPONED*postponed())|(IDEMPOTENT*idempotent());};
    virtual int pushed() = 0;
    virtual int postponed() = 0;
    virtual int idempotent() = 0;
    virtual ~GlobalConstraint();


    virtual void initialise_activity(double *lact, double *vact, double norm);

    virtual void initialise_vars(Solver*);

    virtual void desactivate(const int var) { 

// #ifdef _DEBUG_RELAX
//       std::cout // << "[" << std::setw(4) << id << "](" << name() << "): desactivate "
// 	<< "-"
// 	<< scope[var] << " " ; //<< std::endl;
// #endif
      
      //active.remove(var);
      active.reversible_remove(var);
   }

    void trigger();

    virtual void consolidate() {
      for(unsigned int i=0; i<scope.size; ++i) {
	scope[i] = scope[i].get_var();
      }
     for(unsigned int i=0; i<_scope.size; ++i) {
	_scope[i] = _scope[i].get_var();
      }
    }

    virtual void consolidate_var(const int idx) {

      //std::cout << "consolidate " << scope[idx] << "/" << _scope[idx] << " -> ";

      scope[idx] = scope[idx].get_var();
      _scope[idx] = _scope[idx].get_var();

      //      std::cout << scope[idx] << "/" << _scope[idx] << std::endl;

    }
  
    /// An idempotent constraint should not be called on events triggered by itself.
    /// To forbid that, the lists 'events' and 'changes' are merged
    /// during its propagation, events will be added to the events list of the constraint 
    /// after the propagation, the lists events and changes are swapped back
    /// and the change list is cleared. When idempotent, since the two lists
    /// point to the same object, they are both cleared.
    inline void set_idempotent() {
      if(idempotent()) {
	events.size = 0;
	events.index_capacity = changes.index_capacity;
	events.list_capacity = changes.list_capacity;
	events.list_ = changes.list_;
	events.index_ = changes.index_;
	events.start_ = NULL;
      } else {
	events.initialise(0, changes.index_capacity-1, changes.index_capacity, false);
      }
    }


    virtual void print_active() {
      std::cout << active << std::endl;
    }


    // called before propagation, the events strored in the list 'events'
    // are copied onto the list 'changes'.
    inline void freeze() {
      // if the constraint is idempotent, the two lists point to the same
      // elements, we just set the size up to what it should be.
      changes.size = events.size;

      // otherwise,
      // before each propagation, the lists events and changes are swapped 
      if(changes.list_ != events.list_) {

	// make the changes list points to the events list, and clear the other
	events.size = 0;
	
	int *iaux = events.list_;
	events.list_ = changes.list_;
	changes.list_ = iaux;
	
	unsigned int *uaux = events.index_;
	events.index_ = changes.index_;
	changes.index_ = uaux;
      }
    }

    inline void defrost() {
      //      if(is_posted && active.size <= stress) {

      // #ifdef _DEBUG_CGRAPH
      // 	std::cout << " ---> relax " ;
      // 	display(std::cout);
      // 	std::cout << std::endl;
      // #endif

      // 	relax();
      //       }

      if(enforce_nfc1 && active.size == 1) {
	relax_from(active.back());
      }

      if(changes.list_ == events.list_)
	// if idempotent, clear the events list
	events.size = 0;	
    }

    inline void un_relax() {

#ifdef _DEBUG_RELAX
  if(_DEBUG_RELAX) {
	  std::cout << "[" << std::setw(4) << id << "](" << name() << "): UNRELAX " << std::endl;
  }
#endif

      int i, n;
      for(n=active.size; --n>=0;) {
	i = active[n];

	if(index[i]<0) {
	  
#ifdef _DEBUG_RELAX
  if(_DEBUG_RELAX) {
	  std::cout << "[" << std::setw(4) << id << "](" << name() << "): post on " << scope[i] << std::endl;
  }
#endif
	  
	  index[i] = on[i]->post(self[i]);
	}
	
// #ifdef _DEBUG_RELAX
// 	else {
// 	  std::cout << "SHOULD NOT TRY TO post on " << scope[i] << std::endl;
// 	}
// #endif
	
      }
    }

    inline void post() {

      // save 
      solver->save( Constraint(this, type|POSTED) );
      //solver->save( Constraint(this, type|2) );
    
      active.fill();
      for(int i=on.size; --i>=0;) {
	if(scope[i].is_ground()) {
	  active.reversible_remove(i);
	  index[i] = -1;
	} else {

#ifdef _DEBUG_RELAX
  if(_DEBUG_RELAX) {
      std::cout << "[" << std::setw(4) << id << "](" << name() << "): post on " << scope[i] << std::endl;
  }
#endif

	  index[i] = on[i]->post(self[i]);
	}
      }
    }

    inline void un_post() {

// #ifdef _DEBUG_RELAX
// 	  std::cout << "[" << std::setw(4) << id << "](" << name() << "): UNPOST " << std::endl;
// #endif

//       for(int i=on.size; --i;) {
// 	if(index[i]>=0) {

// #ifdef _DEBUG_RELAX
// 	  std::cout << "[" << std::setw(4) << id << "](" << name() << "): g1 relax from " << scope[i] << std::endl;
// #endif

// 	  un_post_from(i);
// 	  // on[i]->relax(index[i]);
// 	  // index[i] = -1;
// 	}
// #ifdef _DEBUG_RELAX
// 	else {
// 	  std::cout << "SHOULD NOT TRY TO relax from " << scope[i] << std::endl;
// 	}
// #endif
//       }


#ifdef _DEBUG_RELAX
  if(_DEBUG_RELAX) {
      std::cout << "[" << std::setw(4) << id << "](" << name() << "): UNPOST " << std::endl;
  }
#endif
      
      int var;
      for(int i=active.size; --i>=0;) {
	//if(active.contain(i)) {
	var = active[i];
	
	if(index[var]>=0) {
	  
#ifdef _DEBUG_RELAX
  if(_DEBUG_RELAX) {
	  std::cout << "[" << std::setw(4) << id << "](" << name() << "): g1 relax from " << scope[var] << std::endl;
  }
#endif
	  
	  un_post_from(var);
	  // on[i]->relax(index[i]);
	  // index[i] = -1;
	}
	
// #ifdef _DEBUG_RELAX
// 	else {
// 	  std::cout << "SHOULD NOT TRY TO relax from " << scope[var] << std::endl;
// 	}
// #endif
      }
      
    }

    inline void relax() {

#ifdef _DEBUG_RELAX
  if(_DEBUG_RELAX) {
	  std::cout << "[" << std::setw(4) << id << "](" << name() << "): RELAX " << std::endl;
  }
#endif

      //solver->save( Constraint(this, type|1) );
      solver->save( Constraint(this, type|RELAXED) );

//       for(int i=on.size; --i;) {
// 	if(active.contain(i)) {
	  
// #ifdef _DEBUG_RELAX
// 	  std::cout << "[" << std::setw(4) << id << "](" << name() << "): g1 relax from " << scope[i] << std::endl;
// #endif
// 	  if(index[i]>=0) {
// 	    un_post_from(i);
// 	    // on[i]->relax(index[i]);
// 	    // index[i] = -1;
// 	  }
// #ifdef _DEBUG_RELAX
// 	else {
// 	  std::cout << "SHOULD NOT TRY TO relax from " << scope[i] << std::endl;
// 	}
// #endif
// 	}    
//       }

      int var;
      for(int i=active.size; --i>=0;) {
	//if(active.contain(i)) {
	var = active[i];
	  
	if(index[var]>=0) {
	  
#ifdef _DEBUG_RELAX
  if(_DEBUG_RELAX) {
	std::cout << "[" << std::setw(4) << id << "](" << name() << "): g1 relax from " << scope[var] << std::endl;
  }
#endif

	  un_post_from(var);
	  // on[i]->relax(index[i]);
	  // index[i] = -1;
	}
// #ifdef _DEBUG_RELAX
// 	else {
// 	  std::cout << "SHOULD NOT TRY TO relax from " << scope[i] << std::endl;
// 	}
// #endif
      } 
    }


    inline void relax_from(const int var) {
      //if(index[var] > -1) {

#ifdef _DEBUG_RELAX
  if(_DEBUG_RELAX) {
	  std::cout << "[" << std::setw(4) << id << "](" << name() << "): RELAX FROM "  << std::endl;
  }
#endif

      if(index[var]>=0) {
	
	//solver->save( Constraint(this, type|1) );
	solver->save( Constraint(this, type|RELAXED) );
	
	
#ifdef _DEBUG_RELAX
  if(_DEBUG_RELAX) {
	std::cout << "[" << std::setw(4) << id << "](" << name() << "): g2 relax from " << scope[var] << std::endl;
  }
#endif
	
	un_post_from(var);
      
      }    
// #ifdef _DEBUG_RELAX
// 	else {
// 	  std::cout << "SHOULD NOT TRY TO relax from " << scope[var] << std::endl;
// 	}
// #endif
    }

    virtual PropagationOutcome checker_propagate() { return GlobalConstraint::propagate(); }
    virtual PropagationOutcome bound_checker_propagate() { return GlobalConstraint::bound_propagate(); }
    virtual PropagationOutcome propagate();
    virtual PropagationOutcome propagate(const int changed_idx, const Event evt);
    virtual PropagationOutcome bound_propagate();
    virtual PropagationOutcome bound_propagate(const int changed_idx, const Event evt);
    virtual PropagationOutcome checker_propagate(const int changed_idx, 
						 const Event evt) 
    { return GlobalConstraint::propagate(changed_idx, evt); } 
    virtual PropagationOutcome bound_checker_propagate(const int changed_idx, 
						       const Event evt) 
    { return GlobalConstraint::bound_propagate(changed_idx, evt); } 
    virtual int get_backtrack_level();
    virtual Decision get_decision();// { return solver->decisions.back(0); }

    bool first_support(const int vri, const int vli);
    bool find_support(const int vri, const int vli);
    bool find_bound_support(const int vri, const int vli);

    inline void notify_assignment(const int var) {
      
#ifdef _DEBUG_RELAX
  if(_DEBUG_RELAX) {
      std::cout << "[" << std::setw(4) << id << "](" << name() << "): notify assignment of " << scope[var] << " => ";//<< std::endl;
  }
#endif
      
      // if(!active.contain(var)) {
      
      //   std::cout << "SHOULD NOT NOTIFY TWICE (" << scope[var] << " in " << scope[var].get_domain() << ")" << std::endl;
      
      //   //exit(1);
      // } else {
      
      if(active.contain(var)) {
	
	active.reversible_remove(var);
       
	// if(enforce_nfc1 && active.size == 1) {
	//   relax_from(active.back());
	// }
	
      }
      
#ifdef _DEBUG_RELAX
  if(_DEBUG_RELAX) {
      print_active();
      std::cout << std::endl;
  }
#endif
      
    }

    inline void notify_other_event(const int var, 
				   const Mistral::Event evt) {

      

      // if(ASSIGNED(evt) && active.contain(var)) {

      // 	// std::cout << active << " -> " << scope[var] << " in " << scope[var].get_domain() 
      // 	// 	  << " is assigned -> " ;
	
      // 	active.reversible_remove(var);
	
      // 	// std::cout << active << std::endl;

      // }

      if(events.contain(var)) {
	event_type[var] |= evt;
      } else {
	events.add(var);
	event_type[var] = evt;
      }
    }

    inline void notify_first_event(const int var, 
				   const Mistral::Event evt) {

      // if(ASSIGNED(evt) && active.contain(var)) {

      // 	// std::cout << active << " -> " << scope[var] << " in " << scope[var].get_domain() 
      // 	// 	  << " is assigned -> " ;

      // 	active.reversible_remove(var);

      // 	// std::cout << active << std::endl;

      // }

      events.set_to(var);
      event_type[var] = evt;
    }


    void check_active() {
      for(int i=on.size; i--;) {
	if(index[i]>=0) {
	  if(active.contain(i)) {
	    if(scope[i].is_ground()) {
	      std::cout << "[" << std::setw(4) << id << "] " ;
	      display(std::cout);
	      std::cout << std::endl;
	      std::cout << "Warning: " << scope[i] << " = " << scope[i].get_domain()
			<< " is ground and active!!" << std::endl;
	      print_active();
	      std::cout << " :: ";
	      for(unsigned int j=0; j<on.size; ++j)
		std::cout << scope[j] << " in " << scope[j].get_domain() << " ";
	      std::cout << " (exit on check_active())" << std::endl;
	      //exit(1);
	    }
	  } else {
	    if(!(scope[i].is_ground())) {
	      std::cout << "[" << std::setw(4) << id << "] " ;
	      display(std::cout);
	      std::cout << std::endl;
	      std::cout << "Warning: " << scope[i] << " in " << scope[i].get_domain()
			<< " is not ground and not active!!" << std::endl;
	      print_active();
	      std::cout << " :: ";
	      for(unsigned int j=0; j<on.size; ++j)
		std::cout << scope[j] << " in " << scope[j].get_domain() << " ";
	      std::cout << " (exit on check_active())" << std::endl;
	      //exit(1);
	    }
	  }
	}
      }
    }



    virtual void weight_conflict(double unit, Vector<double>& weights) {
      int idx;
      double w = unit
#ifdef _DIV_ARITY
	    / (double)(scope.size)
#endif
      ;

      for(int i=0; i<scope.size; ++i) {
	idx = scope[i].id();

	// std::cout << i << " " << idx << "/" << weights.size << std::endl;
	// if(idx<0 || idx>weights.size) {
	//   std::cout << this << std::endl;
	//   std::cout << scope << std::endl;
	//   std::cout << scope[i] << std::endl;
	//   exit(1);
	// }



	if(idx>=0) // constants have negative ids
	  weights[idx] += w;
      }
    }


    std::ostream& display(std::ostream& os) const;  
  };




  // class ConstraintImplementation;
  // class ConstraintWrapper {




  /**********************************************
   * == Constraint
   **********************************************/ 
  /*! \class ConstraintEqual
    \brief  Binary equal Constraint (x0 == x1).
  */
  class ConstraintEqual : public BinaryConstraint {

  public:  
    /**@name Constructors*/
    //@{
    ConstraintEqual() : BinaryConstraint() {}
    ConstraintEqual(Variable x, Variable y) 
      : BinaryConstraint(x, y) {}
    ConstraintEqual(Vector< Variable >& scp) 
      : BinaryConstraint(scp[0], scp[1]) {}
    ConstraintEqual(std::vector< Variable >& scp) 
      : BinaryConstraint(scp[0], scp[1]) {}
    virtual Constraint clone() { return Constraint(new ConstraintEqual(scope[0], scope[1])// , type
						   ); }
    virtual void initialise();
    virtual void mark_domain();
    virtual int idempotent() { return 1;}
    virtual ~ConstraintEqual() {}
    //@}

    /**@name Solving*/
    //@{
    virtual int check(const int* sol) const { return (sol[0] != sol[1]); }
    virtual PropagationOutcome propagate(const int changed_idx, const Event evt);
    virtual PropagationOutcome propagate();
    virtual bool rewritable() { return true; }
    virtual bool simple_rewritable() { return true; }
    virtual RewritingOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "="; }
    //@}
  };

  /**********************************************
   * =/= Constraint
   **********************************************/ 
  /*! \class ConstraintNotEqual
    \brief  Binary not-equal Constraint (x0 =/= x1).
  */
  class ConstraintNotEqual : public BinaryConstraint {

  public:  
    /**@name Constructors*/
    //@{
    ConstraintNotEqual() : BinaryConstraint() {}
    ConstraintNotEqual(Variable x, Variable y)
      : BinaryConstraint(x, y) {}
    ConstraintNotEqual(Vector< Variable >& scp) 
      : BinaryConstraint(scp// [0], scp[1]
			 ) {}
    ConstraintNotEqual(std::vector< Variable >& scp) 
      : BinaryConstraint(scp[0], scp[1]) {}
    virtual Constraint clone() { return Constraint(new ConstraintNotEqual(scope[0], scope[1])// , type
						   ); }
    virtual void initialise();
    virtual void mark_domain();
    virtual int idempotent() { return 1;}
    virtual ~ConstraintNotEqual() {}
    //@}


    /**@name Solving*/
    //@{
    virtual int check(const int* sol) const { return (sol[0] == sol[1]); }
    virtual PropagationOutcome propagate();
    virtual PropagationOutcome propagate(const int changed_idx, const Event evt);  
    virtual Constraint get_negation(const int var, Variable x) { 
      return Constraint( new ConstraintEqual( (var?scope[0]:x), (var?x:scope[1]) ) );
    }
    virtual bool absorb_negation(const int var) { 
      return (scope[0].get_min()==0 &&
	      scope[1].get_min()==0 &&
	      scope[0].get_max()==1 &&
	      scope[1].get_max()==1);
    }
    virtual bool rewritable() { return true; }
    virtual RewritingOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "=/="; }
    //@}
  };




  /**********************************************
   * Equality Predicate
   **********************************************/
  /*! \class PredicateEqual
    \brief  Truth value of a binary equality ((x0 = x1) <-> y)
  */
  class PredicateEqual : public TernaryConstraint
  {

  public:
    /**@name Parameters*/
    //@{ 
    int spin;
    //@}

    /**@name Constructors*/
    //@{
    PredicateEqual() : TernaryConstraint() {}
    PredicateEqual(Variable x, Variable y, Variable z, const int sp=1) 
      : TernaryConstraint(x, y, z) { spin = sp; }
    PredicateEqual(Vector< Variable >& scp, const int sp=1) 
      : TernaryConstraint(scp) { spin = sp; }
    PredicateEqual(std::vector< Variable >& scp, const int sp=1) 
      : TernaryConstraint(scp) { spin = sp; }
    virtual Constraint clone() { return Constraint(new PredicateEqual(scope[0], scope[1], scope[2], spin)// , type
						   ); }
    virtual void initialise();
    virtual bool rewritable() { return true; }
    virtual void mark_domain();
    virtual int idempotent() { return 1;}
    virtual ~PredicateEqual() {}
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const { return((sol[0] == sol[1]) == (sol[2] ^ spin)); }
    virtual PropagationOutcome propagate();
    virtual PropagationOutcome propagate(const int changed_idx, const Event evt);
    //virtual bool rewritable() { return true; }
    virtual RewritingOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "=?="; }
    //@}
  };


  /**********************************************
   * ConstantEquality Predicate
   **********************************************/
  /*! \class PredicateConstantEqual
    \brief  Truth value of a binary equality ((x0 = x1) <-> y)
  */
  class PredicateConstantEqual : public BinaryConstraint
  {

  public:
    /**@name Parameters*/
    //@{ 
    int spin;
    int value;
    //@}

    /**@name Constructors*/
    //@{
    PredicateConstantEqual() : BinaryConstraint() {}
    PredicateConstantEqual(Variable x, Variable y, const int val=0, const int sp=1) 
      : BinaryConstraint(x,y) { value = val; spin = sp; }
    PredicateConstantEqual(Vector< Variable >& scp, const int val=0, const int sp=1) 
      : BinaryConstraint(scp) { value = val; spin = sp; }
    PredicateConstantEqual(std::vector< Variable >& scp, const int val=0, const int sp=1) 
      : BinaryConstraint(scp) { value = val; spin = sp; }
    virtual Constraint clone() { return Constraint(new PredicateConstantEqual(scope[0], scope[1], value, spin)// , type
						   ); }
    virtual Constraint get_negation(const int var, Variable x) { 
      return Constraint( new PredicateConstantEqual( scope[0], x, value, !spin ) );
    }
    virtual void initialise();
    virtual bool rewritable() { return true; }
    virtual void mark_domain();
    virtual int idempotent() { return 1;}
    virtual bool absorb_negation(const int var) { return var==1; }
    virtual ~PredicateConstantEqual() {}
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const { return((sol[0] == value) == (sol[1] ^ spin)); }
    virtual PropagationOutcome propagate();
    virtual PropagationOutcome propagate(const int changed_idx, const Event evt);
    virtual RewritingOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "=k?"; }
    //@}
  };


  /**********************************************
   * Not Predicate
   **********************************************/
  /*! \class PredicateNot
    \brief  Truth value of a conjunction (!x0 <-> y)
  */
  class PredicateNot : public BinaryConstraint
  {

  public:

    Literal explanation[2];


    /**@name Constructors*/
    //@{
    PredicateNot() : BinaryConstraint() {}
    PredicateNot(Variable x, Variable y) 
      : BinaryConstraint(x,y) {}
    PredicateNot(Vector< Variable >& scp) 
      : BinaryConstraint(scp) {}
    PredicateNot(std::vector< Variable >& scp) 
      : BinaryConstraint(scp) {}
    virtual Constraint clone() { return Constraint(new PredicateNot(scope[0], scope[1])); }
    virtual Constraint get_negation(const int var, Variable x) { 
      return Constraint( new ConstraintEqual( (var?scope[0]:x), (var?x:scope[1]) ) );
    }
    virtual bool explained() { return true; }
    virtual void initialise();
    virtual int idempotent() { return 1; }
    virtual bool absorb_negation(const int var) { return true; }
    virtual ~PredicateNot() {}
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const { 
      return((sol[0]>0) == sol[1]);
    }
    virtual PropagationOutcome propagate();
    virtual PropagationOutcome propagate(const int changed_idx, const Event evt);
    virtual bool rewritable() { return true; }
    virtual RewritingOutcome rewrite();
    //virtual RewritingOutcome rewrite();

    virtual iterator get_reason_for(const Atom a, const int lvl, iterator& end);
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "not?"; }
    //@}
  };


  /**********************************************
   * Neg Predicate
   **********************************************/
  /*! \class PredicateNeg
    \brief  Truth value of a conjunction (!x0 <-> y)
  */
  class PredicateNeg : public BinaryConstraint
  {

  public:
    /**@name Constructors*/
    //@{
    PredicateNeg() : BinaryConstraint() {}
    PredicateNeg(Variable x, Variable y) 
      : BinaryConstraint(x,y) {}
    PredicateNeg(Vector< Variable >& scp) 
      : BinaryConstraint(scp) {}
    PredicateNeg(std::vector< Variable >& scp) 
      : BinaryConstraint(scp) {}
    virtual Constraint clone() { return Constraint(new PredicateNeg(scope[0], scope[1])); }
    virtual void initialise();
    virtual void mark_domain();
    virtual int idempotent() { return 1; }
    virtual ~PredicateNeg() {}
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const { 
      return(sol[0] != -sol[1]);
    }
    virtual PropagationOutcome propagate();
    virtual PropagationOutcome propagate(const int changed_idx, const Event evt);
    //virtual RewritingOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "neg?"; }
    //@}
  };


  /**********************************************
   * Interval Membership Predicate
   **********************************************/
  /*! \class PredicateIntervalMember
    \brief  Truth value of a binary equality ((x0 = x1) <-> y)
  */
  class PredicateIntervalMember : public BinaryConstraint
  {

  public:
    /**@name Parameters*/
    //@{ 
    int spin;
    int lower_bound;
    int upper_bound;
    //@}

    /**@name Constructors*/
    //@{
    PredicateIntervalMember() : BinaryConstraint() {}
    PredicateIntervalMember(Variable x, Variable y, const int lb=-INFTY, const int ub=+INFTY, const int sp=1) 
      : BinaryConstraint(x,y) { spin = sp; lower_bound=lb; upper_bound=ub; }
    PredicateIntervalMember(Vector< Variable >& scp, const int lb=-INFTY, const int ub=+INFTY, const int sp=1) 
      : BinaryConstraint(scp) { spin = sp; lower_bound=lb; upper_bound=ub; }
    PredicateIntervalMember(std::vector< Variable >& scp, const int lb, const int ub, const int sp=1) 
      : BinaryConstraint(scp) { spin = sp; lower_bound=lb; upper_bound=ub; }
    virtual Constraint clone() { return Constraint(new PredicateIntervalMember(scope[0], scope[1], lower_bound, upper_bound, spin)); }
    virtual Constraint get_negation(const int var, Variable x) { 
      return Constraint( new PredicateIntervalMember( scope[0], x, lower_bound, upper_bound, !spin ) );
    }
    virtual void initialise();
    virtual void mark_domain();
    virtual int idempotent() { return 1; }
    virtual bool absorb_negation(const int var) { return var==1; }
    virtual ~PredicateIntervalMember() {}
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const { return(((sol[0]>=lower_bound) && (sol[0]<=upper_bound)) 
						       == (sol[1] ^ spin)); }
    virtual PropagationOutcome propagate();
    virtual PropagationOutcome propagate(const int changed_idx, const Event evt);
    //virtual RewritingOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "in[]?"; }
    //@}
  };


  /**********************************************
   * Set Membership Predicate
   **********************************************/
  /*! \class PredicateSetMember
    \brief  Truth value of a binary equality ((x0 = x1) <-> y)
  */
  class PredicateSetMember : public BinaryConstraint
  {

  public:
    /**@name Parameters*/
    //@{ 
    int spin;
    BitSet values;
    BitSet non_values;
    //@}

    /**@name Constructors*/
    //@{
    PredicateSetMember() : BinaryConstraint() {}
    PredicateSetMember(Vector< Variable >& scp, const int sp=1) 
      : BinaryConstraint(scp) { spin = sp; }
    PredicateSetMember(Variable x, Variable y, const BitSet& vals, const int sp=1) 
      : BinaryConstraint(x,y) { spin = sp; values=vals; }
    PredicateSetMember(Vector< Variable >& scp, const BitSet& vals, const int sp=1) 
      : BinaryConstraint(scp) { spin = sp; values=vals; }
    PredicateSetMember(std::vector< Variable >& scp, const BitSet& vals, const int sp=1) 
      : BinaryConstraint(scp) { spin = sp; values=vals; }
    virtual Constraint clone() { return Constraint(new PredicateSetMember(scope[0], scope[1], values, spin)); }
    virtual Constraint get_negation(const int var, Variable x) { 
      return Constraint( new PredicateSetMember( scope[0], x, values, !spin ) );
    }
    virtual void initialise();
    virtual void mark_domain();
    virtual int idempotent() { return 1; }
    virtual bool absorb_negation(const int var) { return var==1; }
    virtual ~PredicateSetMember() {}
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const { return((values.contain(sol[0]) == (sol[1] ^ spin))); }
    virtual PropagationOutcome propagate(const int changed_idx, const Event evt);
    virtual PropagationOutcome propagate();
    //virtual RewritingOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "in{}?"; }
    //@}
  };


  /**********************************************
   * <= Constraint
   **********************************************/ 
  /*! \class ConstraintLess
    \brief  Binary Less Than Constraint (x0 + k <= x1).
  */
  class ConstraintLess : public BinaryConstraint {

  public: 
    /**@name Parameters*/
    //@{  
    int offset;
    //@}

    /**@name Constructors*/
    //@{
    ConstraintLess() : BinaryConstraint() {}
    ConstraintLess(const int ofs=0) 
      : BinaryConstraint() { offset = ofs; }
    ConstraintLess(Variable x, Variable y, const int ofs=0) 
      : BinaryConstraint(x, y) { offset = ofs; }
    ConstraintLess(Vector< Variable >& scp, const int ofs=0) 
      : BinaryConstraint(scp) { offset = ofs; }
    ConstraintLess(std::vector< Variable >& scp, const int ofs=0) 
      : BinaryConstraint(scp) { offset = ofs; }
    virtual Constraint clone() { return Constraint(new ConstraintLess(scope[0], scope[1], offset)// , type
						   ); }
    virtual void initialise();
    virtual void mark_domain();
    virtual int idempotent() { return 1;}
    // virtual bool absorb_negation(const int var) { 
    //   return (offset = 0 &&
    // 	      scope[0].get_min()==0 &&
    // 	      scope[1].get_min()==0 &&
    // 	      scope[0].get_max()==1 &&
    // 	      scope[1].get_max()==1);
    // }
    virtual ~ConstraintLess() {}
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const { return (sol[0]+offset > sol[1]); }
    virtual PropagationOutcome propagate();
    virtual PropagationOutcome propagate(const int changed_idx, const Event evt);
    //virtual RewritingOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "<"; }
    //@}
  };

  /**********************************************
   * <= Predicate
   **********************************************/ 
  /*! \class PredicateLess

    \brief  Truth value of a precedence ((x0 + k <= x1) <-> y)
  */
  class PredicateLess : public TernaryConstraint {

  public: 
    /**@name Parameters*/
    //@{  
    int offset;
    //@}

    /**@name Constructors*/
    //@{
    PredicateLess() : TernaryConstraint(), offset(0) {}
    PredicateLess(Variable x, Variable y, Variable z, const int ofs=0)
      : TernaryConstraint(x, y, z), offset(ofs) {}
    PredicateLess(Vector< Variable >& scp, const int ofs=0) 
      : TernaryConstraint(scp), offset(ofs) {}
    PredicateLess(std::vector< Variable >& scp, const int ofs=0) 
      : TernaryConstraint(scp), offset(ofs) {}
    virtual Constraint clone() { return Constraint(new PredicateLess(scope[0], scope[1], scope[2], offset)// , type
						   ); }
    virtual Constraint get_negation(const int var, Variable x) { 
      return Constraint( new PredicateLess( scope[1], scope[0], x, 1-offset ) );
    }
    virtual void initialise();
    virtual void mark_domain();
    virtual int idempotent() { return 1;}
    virtual bool absorb_negation(const int var) { return var==2; }
    virtual ~PredicateLess() {}
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const { 
      return ((sol[0]+offset <= sol[1]) != sol[2]); 
    }
    virtual PropagationOutcome propagate();
    virtual PropagationOutcome propagate(const int changed_idx, const Event evt);
    //@}
    
    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "<?"; }
    //@}
  };



  /**********************************************
   * <=K Predicate
   **********************************************/ 
  /*! \class PredicateUpperBound

    \brief  Truth value of a precedence ((x0 <= k) <-> y)
  */
  class PredicateUpperBound : public BinaryConstraint {

  public: 

    int bound;

    /**@name Constructors*/
    //@{
    PredicateUpperBound() : BinaryConstraint() {}
    PredicateUpperBound(Variable x, Variable y, const int b=0) 
      : BinaryConstraint(x,y) { bound = b; }
    PredicateUpperBound(Vector< Variable >& scp, const int b=0) 
      : BinaryConstraint(scp) { bound = b; }
    PredicateUpperBound(std::vector< Variable >& scp, const int b=0) 
      : BinaryConstraint(scp) { bound = b; }
    virtual Constraint clone() { return Constraint(new PredicateUpperBound(scope[0], scope[1], bound)// , type
						   ); }
    virtual Constraint get_negation(const int var, Variable x);
    // { 
    //   return Constraint( new PredicateLowerBound( scope[0], x, bound-1 ) );
    // }
    virtual void initialise();
    virtual void mark_domain();
    virtual int idempotent() { return 1;}
    virtual bool absorb_negation(const int var) { return var==1; }
    virtual ~PredicateUpperBound() {}
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const { 
      return ((sol[0] <= bound) != sol[1]); 
    }
    virtual PropagationOutcome propagate();
    virtual PropagationOutcome propagate(const int changed_idx, const Event evt);
    //virtual RewritingOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "<u?"; }
    //@}
  };


  /**********************************************
   * >=K Predicate
   **********************************************/ 
  /*! \class PredicateLowerBound

    \brief  Truth value of a precedence ((x0 >= k) <-> y)
  */
  class PredicateLowerBound : public BinaryConstraint {

  public: 

    int bound;

    /**@name Constructors*/
    //@{
    PredicateLowerBound() : BinaryConstraint() {}
    PredicateLowerBound(Variable x, Variable y, const int b=0) 
      : BinaryConstraint(x,y) { bound = b; }
    PredicateLowerBound(Vector< Variable >& scp, const int b=0) 
      : BinaryConstraint(scp) { bound = b; }
    PredicateLowerBound(std::vector< Variable >& scp, const int b=0) 
      : BinaryConstraint(scp) { bound = b; }
    virtual Constraint clone() { return Constraint(new PredicateLowerBound(scope[0], scope[1], bound)// , type
						   ); }
    virtual Constraint get_negation(const int var, Variable x) { 
      return Constraint( new PredicateUpperBound( scope[0], x, bound+1 ) );
    }
    virtual void initialise();
    virtual void mark_domain();
    virtual int idempotent() { return 1;}
    virtual bool absorb_negation(const int var) { return var==1; }
    virtual ~PredicateLowerBound() {}
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const { 
      return ((sol[0] >= bound) != sol[1]); 
    }
    virtual PropagationOutcome propagate();
    virtual PropagationOutcome propagate(const int changed_idx, const Event evt);
    //virtual RewritingOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return ">l?"; }
    //@}
  };



  /**********************************************
   * And Predicate
   **********************************************/
  /*! \class PredicateAnd
    \brief  Truth value of a conjunction ((x0 and x1) <-> y)
  */
  class PredicateAnd : public TernaryConstraint
  {

  public:
    /**@name Constructors*/
    //@{
    PredicateAnd() : TernaryConstraint() {}
    PredicateAnd(Vector< Variable >& scp) 
      : TernaryConstraint(scp) {}
    PredicateAnd(Variable x, Variable y, Variable z) 
      : TernaryConstraint(x,y,z) {}
    PredicateAnd(std::vector< Variable >& scp) 
      : TernaryConstraint(scp) {}
    virtual Constraint clone() { return Constraint(new PredicateAnd(scope[0], scope[1], scope[2])); }
    // virtual Constraint get_negation(const int var, Variable x) { 
    //   return Constraint( new PredicateUpperBound( scope[0], x, bound+1 ) );
    // }
    virtual void initialise();
    //virtual void mark_domain();
    virtual int idempotent() { return 1;}
    //virtual bool absorb_negation(const int var) { return true; }
    virtual ~PredicateAnd() {}
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const { 
      return((sol[0] && sol[1]) != (sol[2])); 
    }
    virtual PropagationOutcome propagate();
    virtual PropagationOutcome propagate(const int changed_idx, const Event evt);
    //virtual RewritingOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "and?"; }
    //@}
  };

  /**********************************************
   * And Constraint
   **********************************************/
  /*! \class ConstraintAnd
    \brief  Conjunction (x0 and x1)
  */
  class ConstraintAnd : public BinaryConstraint
  {

  public:
    /**@name Constructors*/
    //@{
    ConstraintAnd() : BinaryConstraint() {}
    ConstraintAnd(Vector< Variable >& scp) 
      : BinaryConstraint(scp) {}
    ConstraintAnd(Variable x, Variable y) 
      : BinaryConstraint(x,y) {}
    ConstraintAnd(std::vector< Variable >& scp) 
      : BinaryConstraint(scp) {}
    virtual Constraint clone() { return Constraint(new ConstraintAnd(scope[0], scope[1])); }
    virtual void initialise();
    //virtual void mark_domain();
    virtual int idempotent() { return 1;}
    //virtual bool absorb_negation(const int var) { return true; }
    virtual ~ConstraintAnd() {}
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const { 
      return(!(sol[0] && sol[1])); 
    }
    virtual PropagationOutcome propagate();
    virtual PropagationOutcome propagate(const int changed_idx, const Event evt);
    //virtual RewritingOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "and"; }
    //@}
  };


  /**********************************************
   * NotAnd Constraint
   **********************************************/
  /*! \class ConstraintNotAnd
    \brief  Conjunction (x0 and x1)
  */
  class ConstraintNotAnd : public BinaryConstraint
  {

  public:
    /**@name Constructors*/
    //@{
    ConstraintNotAnd() : BinaryConstraint() {}
    ConstraintNotAnd(Vector< Variable >& scp) 
      : BinaryConstraint(scp) {}
    ConstraintNotAnd(Variable x, Variable y) 
      : BinaryConstraint(x,y) {}
    ConstraintNotAnd(std::vector< Variable >& scp) 
      : BinaryConstraint(scp) {}
    virtual Constraint clone() { return Constraint(new ConstraintNotAnd(scope[0], scope[1])); }
    virtual void initialise();
    //virtual void mark_domain();
    virtual int idempotent() { return 1;}
    //virtual bool absorb_negation(const int var) { return true; }
    virtual ~ConstraintNotAnd() {}
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const { 
      return(sol[0] && sol[1]); 
    }
    virtual PropagationOutcome propagate();
    virtual PropagationOutcome propagate(const int changed_idx, const Event evt);
    //virtual RewritingOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "!and"; }
    //@}
  };


  /**********************************************
   * Or Predicate
   **********************************************/
  /*! \class PredicateOr
    \brief  Truth value of a disjunction ((x0 or x1) <-> y)
  */
  class PredicateOr : public TernaryConstraint
  {

  public:
    /**@name Constructors*/
    //@{
    PredicateOr() : TernaryConstraint() {}
    PredicateOr(Vector< Variable >& scp) 
      : TernaryConstraint(scp) {}
    PredicateOr(Variable x, Variable y, Variable z) 
      : TernaryConstraint(x,y,z) {}
    PredicateOr(std::vector< Variable >& scp) 
      : TernaryConstraint(scp) {}
    virtual Constraint clone() { return Constraint(new PredicateOr(scope[0], scope[1], scope[2])); }
    virtual Constraint get_negation(const int var, Variable x) { 
      return Constraint( (var<2 ? 
			  new PredicateLess( x, (var?scope[0]:scope[1]), scope[2], 0 ) : 
			  NULL
			  //new PredicateNotOr(scope[0], scope[1], scope[2]) 
			  ) );
    }
    virtual void initialise();
    //virtual void mark_domain();
    virtual int idempotent() { return 1;}
    virtual bool absorb_negation(const int var) { return true; }
    virtual ~PredicateOr() {}
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const { 
      return((sol[0] || sol[1]) != (sol[2])); 
    }
    virtual PropagationOutcome propagate();
    virtual PropagationOutcome propagate(const int changed_idx, const Event evt);
    //virtual RewritingOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "or?"; }
    //@}
  };

  /**********************************************
   * Or Constraint
   **********************************************/
  /*! \class ConstraintOr
    \brief  Disjunction (x0 or x1)
  */
  class ConstraintOr : public BinaryConstraint
  {

  public:
    /**@name Constructors*/
    //@{
    ConstraintOr(Vector< Variable >& scp) 
      : BinaryConstraint(scp) { }
    ConstraintOr(Variable x, Variable y) 
      : BinaryConstraint(x,y) {}
    ConstraintOr(std::vector< Variable >& scp) 
      : BinaryConstraint(scp) { }
    virtual Constraint clone() { return Constraint(new ConstraintOr(scope[0], scope[1])); }
    virtual Constraint get_negation(const int var, Variable x) { 
      return Constraint( new ConstraintLess( x, (var?scope[0]:scope[1]), 0 ) );
    }
    virtual void initialise();
    //virtual void mark_domain();
    virtual int idempotent() { return 1;}
    virtual bool absorb_negation(const int var) { return true; }
    virtual ~ConstraintOr() {}
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const { 
      return(!(sol[0] || sol[1])); 
    }
    virtual PropagationOutcome propagate();
    virtual PropagationOutcome propagate(const int changed_idx, const Event evt);
    //virtual RewritingOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "or"; }
    //@}
  };


  /**********************************************
   * Lex Constraint
   **********************************************/
  /*! \class ConstraintLex
    \brief  Basic bloc of a lex lt/leq constraint

    let x0 and x1 being two cells with same rank on two rows/columns 
    and b0, b1 be two Boolean variables.
    This constraint ensures that 
       - x0 =/= x1 => b1
       - b0 < b1 => x0 < x1
       - b0 <= b1 
  */
  class ConstraintLex : public GlobalConstraint
  {

  public:
    /**@name Constructors*/
    //@{
    ConstraintLex() : GlobalConstraint() { priority=1; }
    ConstraintLex(Vector< Variable >& scp) 
      : GlobalConstraint(scp) { priority=1; }
    ConstraintLex(std::vector< Variable >& scp) 
      : GlobalConstraint(scp) { priority=1; }
    virtual Constraint clone() { return Constraint(new ConstraintLex(scope)); }
    virtual void initialise();
    virtual void mark_domain();
    virtual int idempotent() { return 1;}
    virtual int postponed() { return 1;}
    virtual int pushed() { return 1;}
    virtual ~ConstraintLex() {}
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const { 
      return( ((!sol[2] && !sol[3]) > (sol[0] == sol[1])
	       || (sol[2] < sol[3]) > (sol[0] <  sol[1])
	       || sol[2] > sol[3])
	      ); 
    }
    virtual PropagationOutcome propagate();
    //virtual RewritingOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "lex"; }
    //@}
  };


  /**********************************************
   * Offset Predicate
   **********************************************/
  /*! \class PredicateOffset
    \brief  Offset (x+k = y)
  */
  class PredicateOffset : public BinaryConstraint
  {

  public:

    /**@name Parameters*/
    //@{  
    int offset;
    //@}

    /**@name Constructors*/
    //@{
    PredicateOffset() : BinaryConstraint() {}
    PredicateOffset(Vector< Variable >& scp, const int ofs) 
      : BinaryConstraint(scp) { offset=ofs; }
    PredicateOffset(Variable x, Variable y, const int ofs) 
      : BinaryConstraint(x,y) { offset=ofs; }
    PredicateOffset(std::vector< Variable >& scp, const int ofs) 
      : BinaryConstraint(scp) { offset=ofs; }
    virtual Constraint clone() { return Constraint(new PredicateOffset(scope[0], scope[1], offset)); }
    virtual void initialise();
    virtual void mark_domain();
    virtual int idempotent() { return 1;}
    virtual ~PredicateOffset() {}
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const { 
      return((sol[0]+offset) != sol[1]);
    }
    virtual PropagationOutcome propagate(const int changed_idx, const Event evt);
    virtual PropagationOutcome propagate();
    //virtual RewritingOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "+k"; }
    //@}
  };


  /**********************************************
   * Addition Predicate
   **********************************************/
  /*! \class PredicateAdd
    \brief  Binary addition predicate (x0 + x1 = y)
  */
  class PredicateAdd : public TernaryConstraint
  {
    
  public:
    /**@name Constructors*/
    //@{
    PredicateAdd() : TernaryConstraint() {}
    PredicateAdd(Variable x, Variable y, Variable z)
      : TernaryConstraint(x, y, z) {}
    PredicateAdd(Vector< Variable >& scp)
      : TernaryConstraint(scp) {}
    PredicateAdd(std::vector< Variable >& scp)
      : TernaryConstraint(scp) {}
    virtual Constraint clone() { return Constraint(new PredicateAdd(scope[0], scope[1], scope[2])// , type
						   ); }
    virtual void initialise();
    virtual void mark_domain();
    virtual int idempotent() { return 0;}
    virtual ~PredicateAdd() {}
    //@}
    
    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const { return (sol[2] != (sol[0]+sol[1])); }
    virtual PropagationOutcome propagate(const int changed_idx, const Event evt);
    virtual PropagationOutcome propagate();
    virtual RewritingOutcome rewrite();
    //@}
    
    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "+="; }
    //@}
  };


  /**********************************************
   * Factor Predicate
   **********************************************/
  /*! \class PredicateFactor
    \brief  Truth value of a conjunction (!x0 <-> y)
  */
  class PredicateFactor : public BinaryConstraint
  {

  public:

    /**@name Parameters*/
    //@{  
    int factor;
    //@}

    /**@name Constructors*/
    //@{
    PredicateFactor() : BinaryConstraint() {}
    PredicateFactor(Vector< Variable >& scp, const int fct) 
      : BinaryConstraint(scp) { factor=fct; }
    PredicateFactor(Variable x, Variable y, const int fct) 
      : BinaryConstraint(x,y) { factor=fct; }
    PredicateFactor(std::vector< Variable >& scp, const int fct) 
      : BinaryConstraint(scp) { factor=fct; }
    virtual Constraint clone() { return Constraint(new PredicateFactor(scope[0], scope[1], factor)); }
    virtual void initialise();
    virtual void mark_domain();
    virtual int idempotent() { return 1;}
    virtual ~PredicateFactor() {}
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const { 
      return((sol[0]*factor) != sol[1]);
    }
    virtual PropagationOutcome propagate(const int changed_idx, const Event evt);
    virtual PropagationOutcome propagate();
    //virtual RewritingOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "*k"; }
    //@}
  };


 /**********************************************
  * Square Predicate
  **********************************************/
 /*! \class PredicateSquare
   \brief  Truth value of a conjunction (!x0 <-> y)
 */
 class PredicateSquare : public BinaryConstraint
 {

 public:

   /**@name Parameters*/
   //@{  
   //@}

   /**@name Constructors*/
   //@{
   PredicateSquare() : BinaryConstraint() {}
   PredicateSquare(Vector< Variable >& scp) 
     : BinaryConstraint(scp) {  }
   PredicateSquare(Variable x, Variable y) 
     : BinaryConstraint(x,y) {  }
   PredicateSquare(std::vector< Variable >& scp) 
     : BinaryConstraint(scp) {  }
   virtual Constraint clone() { return Constraint(new PredicateSquare(scope[0], scope[1])); }
   virtual void initialise();
   virtual void mark_domain();
   virtual int idempotent() { return 1;}
   virtual ~PredicateSquare() {}
   //@}

   /**@name Solving*/
   //@{
   virtual int check( const int* sol ) const { 
     return((sol[0]*sol[0]) != sol[1]);
   }
   virtual PropagationOutcome propagate(const int changed_idx, const Event evt);
   virtual PropagationOutcome propagate();
   //virtual RewritingOutcome rewrite();
   //@}

   /**@name Miscellaneous*/
   //@{  
   virtual std::ostream& display(std::ostream&) const ;
   virtual std::string name() const { return "^2"; }
   //@}
 };



  /**********************************************
   * Abs Predicate
   **********************************************/
  /*! \class PredicateAbs
    \brief  Absolute value of a variable
  */
  class PredicateAbs : public BinaryConstraint
  {

  public:

    /**@name Parameters*/
    //@{  
    BitSet buffer;
    //@}


    /**@name Constructors*/
    //@{
    PredicateAbs() : BinaryConstraint() {}
    PredicateAbs(Vector< Variable >& scp) 
      : BinaryConstraint(scp) { }
    PredicateAbs(Variable x, Variable y) 
      : BinaryConstraint(x,y) { }
    PredicateAbs(std::vector< Variable >& scp) 
      : BinaryConstraint(scp) { }
    virtual Constraint clone() { return Constraint(new PredicateAbs(scope[0], scope[1])); }
    virtual void initialise();
    virtual void mark_domain();
    virtual int idempotent() { return 1;}
    virtual ~PredicateAbs() {}
    //@}

    /**@name Solving*/
    //@{
    PropagationOutcome propagate_change_on_X( const Event evt );
    PropagationOutcome propagate_change_on_absX( const Event evt );
    virtual int check( const int* sol ) const { 
      return(abs(sol[0]) != sol[1]);
    }
    virtual PropagationOutcome propagate(const int changed_idx, const Event evt);
    virtual PropagationOutcome propagate();
    //virtual RewritingOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "abs"; }
    //@}
  };


  /**********************************************
   * Modulo Predicate
   **********************************************/
  /*! \class PredicateMod
    \brief  Binary Modulo predicate (x0 % x1 = y)
  */
  class PredicateMod : public TernaryConstraint
  {

 public:
    /**@name Constructors*/
    //@{
    PredicateMod() : TernaryConstraint() {}
    PredicateMod(Variable x, Variable y, Variable z)
      : TernaryConstraint(x, y, z) {}
    PredicateMod(Vector< Variable >& scp)
      : TernaryConstraint(scp) {}
    PredicateMod(std::vector< Variable >& scp)
      : TernaryConstraint(scp) {}
    virtual Constraint clone() { return Constraint(new PredicateMod(scope[0], scope[1], scope[2])); }
    virtual void initialise();
    virtual void mark_domain();
    virtual int idempotent() { return 0; }
    virtual ~PredicateMod() {}
    //@}
    
    /**@name Solving*/
    //@{
    PropagationOutcome filter();
    virtual int check( const int* sol ) const { return (!sol[1] || sol[2] != __modulo_fct__(sol[0],sol[1])); }
    virtual PropagationOutcome propagate(const int changed_idx, const Event evt);
    virtual PropagationOutcome propagate();
    virtual RewritingOutcome rewrite();
    //@}
    
    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "%="; }
    //@}

  };



  /**********************************************
   * CModulo Predicate
   **********************************************/
  /*! \class PredicateCMod
    \brief  Binary CModulo predicate (x0 % x1 = y)
  */
  class PredicateCMod : public TernaryConstraint
  {

 public:
    /**@name Constructors*/
    //@{
    PredicateCMod() : TernaryConstraint() {}
    PredicateCMod(Variable x, Variable y, Variable z)
      : TernaryConstraint(x, y, z) {}
    PredicateCMod(Vector< Variable >& scp)
      : TernaryConstraint(scp) {}
    PredicateCMod(std::vector< Variable >& scp)
      : TernaryConstraint(scp) {}
    virtual Constraint clone() { return Constraint(new PredicateCMod(scope[0], scope[1], scope[2])); }
    virtual void initialise();
    virtual void mark_domain();
    virtual int idempotent() { return 0; }
    virtual ~PredicateCMod() {}
    //@}
    
    /**@name Solving*/
    //@{
    PropagationOutcome filter();
    virtual int check( const int* sol ) const { return (!sol[1] || sol[2] != sol[0] %sol[1]); }
    virtual PropagationOutcome propagate(const int changed_idx, const Event evt);
    virtual PropagationOutcome propagate();
    virtual RewritingOutcome rewrite();
    //@}
    
    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "%="; }
    //@}

  };


  /**********************************************
   * Modulo Predicate
   **********************************************/
  /*! \class PredicateMod
    \brief  Binary Modulo predicate (x0 % x1 = y)
  */
  class PredicateModConstant : public BinaryConstraint
  {

  public:
   /**@name Parameters*/
    //@{  
    int modulo;
    //@}

    /**@name Constructors*/
    //@{
    PredicateModConstant() : BinaryConstraint() {}
    PredicateModConstant(Vector< Variable >& scp, const int mod=2) 
      : BinaryConstraint(scp) { modulo=mod; }
    PredicateModConstant(Variable x, Variable y, const int mod=2) 
      : BinaryConstraint(x,y) { modulo=mod; }
    PredicateModConstant(std::vector< Variable >& scp, const int mod=2) 
      : BinaryConstraint(scp) { modulo=mod; }
    virtual Constraint clone() { return Constraint(new PredicateModConstant(scope[0], scope[1], modulo)); }
    virtual void initialise();
    virtual void mark_domain();
    virtual int idempotent() { return 0;}
    virtual ~PredicateModConstant() {}
    //@}

    /**@name Solving*/
    //@{
    PropagationOutcome filter();
    virtual int check( const int* sol ) const { 
      //return((sol[0] % modulo) != sol[1]);
      return(__modulo_fct__(sol[0],modulo) != sol[1]);
    }
    virtual PropagationOutcome propagate(const int changed_idx, const Event evt);
    virtual PropagationOutcome propagate();
    //virtual RewritingOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "%k"; }
    //@}
  };



  /**********************************************
   * CModulo Predicate
   **********************************************/
  /*! \class PredicateCMod
    \brief  Binary CModulo predicate (x0 % x1 = y)
  */
  class PredicateCModConstant : public BinaryConstraint
  {

  public:
   /**@name Parameters*/
    //@{  
    int modulo;
    //@}

    /**@name Constructors*/
    //@{
    PredicateCModConstant() : BinaryConstraint() {}
    PredicateCModConstant(Vector< Variable >& scp, const int mod=2) 
      : BinaryConstraint(scp) { modulo=mod; }
    PredicateCModConstant(Variable x, Variable y, const int mod=2) 
      : BinaryConstraint(x,y) { modulo=mod; }
    PredicateCModConstant(std::vector< Variable >& scp, const int mod=2) 
      : BinaryConstraint(scp) { modulo=mod; }
    virtual Constraint clone() { return Constraint(new PredicateCModConstant(scope[0], scope[1], modulo)); }
    virtual void initialise();
    virtual void mark_domain();
    virtual int idempotent() { return 0;}
    virtual ~PredicateCModConstant() {}
    //@}

    /**@name Solving*/
    //@{
    PropagationOutcome filter();
    virtual int check( const int* sol ) const { 
      //return((sol[0] % modulo) != sol[1]);
      return(sol[0] % modulo != sol[1]);
    }
    virtual PropagationOutcome propagate(const int changed_idx, const Event evt);
    virtual PropagationOutcome propagate();
    //virtual RewritingOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "%k"; }
    //@}
  };



  /**********************************************
   * Quotient Predicate
   **********************************************/
  /*! \class PredicateMod
    \brief  Binary Quotient predicate (x0 % x1 = y)
  */
  class PredicateDivConstant : public BinaryConstraint
  {

  public:
   /**@name Parameters*/
    //@{  
    int quotient;
    //@}

    /**@name Constructors*/
    //@{
    PredicateDivConstant() : BinaryConstraint() {}
    PredicateDivConstant(Vector< Variable >& scp, const int q=2) 
      : BinaryConstraint(scp) { quotient=q; }
    PredicateDivConstant(Variable x, Variable y, const int q=2) 
      : BinaryConstraint(x,y) { quotient=q; }
    PredicateDivConstant(std::vector< Variable >& scp, const int q=2) 
      : BinaryConstraint(scp) { quotient=q; }
    virtual Constraint clone() { return Constraint(new PredicateDivConstant(scope[0], scope[1], quotient)); }
    virtual void initialise();
    virtual void mark_domain();
    virtual int idempotent() { return 0;}
    virtual ~PredicateDivConstant() {}
    //@}

    /**@name Solving*/
    //@{
    PropagationOutcome filter();
    virtual int check( const int* sol ) const { 
      //return((sol[0] % quotient) != sol[1]);
      return((sol[0] / quotient) != sol[1]);
    }
    virtual PropagationOutcome propagate(const int changed_idx, const Event evt);
    virtual PropagationOutcome propagate();
    //virtual RewritingOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "%k"; }
    //@}
  };



  /**********************************************
   * Division Predicate
   **********************************************/
  /*! \class PredicateDiv
    \brief  Binary Division predicate (x0 / x1 = x2)
  */
  class PredicateDiv : public TernaryConstraint
  {

 public:
    /**@name Constructors*/
    //@{
    PredicateDiv() : TernaryConstraint() {}
    PredicateDiv(Variable x, Variable y, Variable z)
      : TernaryConstraint(x, y, z) {}
    PredicateDiv(Vector< Variable >& scp)
      : TernaryConstraint(scp) {}
    PredicateDiv(std::vector< Variable >& scp)
      : TernaryConstraint(scp) {}
    virtual Constraint clone() { return Constraint(new PredicateDiv(scope[0], scope[1], scope[2])); }
    virtual void initialise();
    virtual void mark_domain();
    virtual int idempotent() { return 0; }
    virtual ~PredicateDiv() {}
    //@}
    
    /**@name Solving*/
    //@{
    PropagationOutcome filter();
    virtual int check( const int* sol ) const { return (!sol[1] || sol[2] != (sol[0] / sol[1])); }
    virtual PropagationOutcome propagate(const int changed_idx, const Event evt);
    virtual PropagationOutcome propagate();
    //virtual RewritingOutcome rewrite();
    //@}
    
    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "%="; }
    //@}

  };



  /**********************************************
   * Multiplication Predicate
   **********************************************/
  /*! \class PredicateMul
    \brief  Binary mulition predicate (x0 + x1 = y)
  */
  class PredicateMul : public GlobalConstraint
  {
    
  public:

    int min_pos[3];
    int max_pos[3];
    int min_neg[3];
    int max_neg[3];
    int zero[3];

    /**@name Constructors*/
    //@{
    PredicateMul() : GlobalConstraint() { priority = 1; }
    PredicateMul(Vector< Variable >& scp) 
      : GlobalConstraint(scp) { priority = 1; }
    // PredicateMul(Variable x, Variable y, Variable z) 
    //   : GlobalConstraint(x,y,z) {}
    PredicateMul(std::vector< Variable >& scp) 
      : GlobalConstraint(scp) { priority = 1; }
    //virtual Constraint clone() { return Constraint(new PredicateMul(scope[0], scope[1], scope[2])); }
    virtual Constraint clone() { return Constraint(new PredicateMul(scope)); }
    virtual void initialise();
    virtual void mark_domain();
    virtual int idempotent() { return 1;}
    virtual int postponed() { return 1;}
    virtual int pushed() { return 1;}
    virtual ~PredicateMul() {}
    //@}
    
    /**@name Solving*/
    //@{
    //
    PropagationOutcome revise_division(const int X, const int Y, const int Z);
    //
    PropagationOutcome revise_multiplication(const int X, const int Y, const int Z);
    //
    PropagationOutcome prune(const int lb_neg, 
			     const int ub_neg, 
			     const int lb_pos, 
			     const int ub_pos,
			     const bool pzero,
			     const int Z);

    virtual int check( const int* sol ) const { return (sol[2] != (sol[0]*sol[1])); }
    virtual PropagationOutcome propagate();
    virtual RewritingOutcome rewrite();
    //@}
    
    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "*="; }
    //@}
  };




  /**********************************************
   * Disjunctive Constraint (implemented with two precedence constraints)
   **********************************************/ 
  /*! \class ConstraintDisjunctive
    \brief  Binary Disjunctive Constraint (x0 + p0 <= x1 || x1 + p1 <= x0).
  */
  class ConstraintDisjunctive : public BinaryConstraint {
    
  public: 
    /**@name Parameters*/
    //@{
    int processing_time[2];
    Constraint precedence[2];
    //@}
    
    /**@name Constructors*/
    //@{
    ConstraintDisjunctive() : BinaryConstraint() {}
    ConstraintDisjunctive(Variable x, Variable y, const int p0, const int p1); 
    ConstraintDisjunctive(Vector< Variable >& scp, const int p0, const int p1); 
    ConstraintDisjunctive(std::vector< Variable >& scp, const int p0, const int p1); 
    virtual Constraint clone() { 
      return Constraint(new ConstraintDisjunctive(scope[0], scope[1], 
						  processing_time[0], processing_time[1])); }
    virtual void initialise();
    virtual void mark_domain();
    virtual int idempotent() { return 1; }
    virtual ~ConstraintDisjunctive() {}
    //@}

    /**@name Solving*/
    //@{
    void decide(const int choice);
    virtual int check( const int* sol ) const { 
      return ((sol[0]+processing_time[0] > sol[1])
	      &&
	      (sol[1]+processing_time[1] > sol[0])); 
    }
    virtual PropagationOutcome propagate(const int changed_idx, const Event evt);
    virtual PropagationOutcome propagate();
    //virtual RewritingOutcome rewrite();
    virtual void consolidate();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "<>"; }
    //@}
  };


  /**********************************************
   * ReifiedDisjunctive Constraint (implemented with two precedence constraints)
   **********************************************/ 
  /*! \class ConstraintReifiedDisjunctive
    \brief  Binary ReifiedDisjunctive Constraint (x0 + p0 <= x1 || x1 + p1 <= x0).
  */
  class ConstraintReifiedDisjunctive : public TernaryConstraint {

  public: 
    /**@name Parameters*/
    //@{
    int processing_time[2];
    int *min_t0_ptr;
    int *max_t0_ptr;
    int *min_t1_ptr;
    int *max_t1_ptr;
    int *state;
    //@}

    /**@name Constructors*/
    //@{
    ConstraintReifiedDisjunctive() : TernaryConstraint() {}
    ConstraintReifiedDisjunctive(Variable x, Variable y, Variable z, const int p0, const int p1); 
    ConstraintReifiedDisjunctive(Vector< Variable >& scp, const int p0, const int p1); 
    ConstraintReifiedDisjunctive(std::vector< Variable >& scp, const int p0, const int p1); 

    virtual Constraint clone() { return Constraint(new ConstraintReifiedDisjunctive(scope[0], scope[1], scope[2], processing_time[0], processing_time[1])); }
    virtual Constraint get_negation(const int var, Variable x) { 
      return Constraint( new ConstraintReifiedDisjunctive( scope[0], scope[1], x, processing_time[0], processing_time[1] ) );
    }
    virtual void initialise();
    virtual void mark_domain();
    virtual ~ConstraintReifiedDisjunctive() {}
    virtual int idempotent() { return 1; }
    virtual bool absorb_negation(const int var) { return var==2; }
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const { 
      int ret_value = 0;
      if(sol[2]) {
	ret_value = (sol[0] + processing_time[0] > sol[1]);
      } else {
	ret_value = (sol[1] + processing_time[1] > sol[0]);
      }
      return ret_value;

      // return ( (!sol[2] && (sol[1] + processing_time[1] > sol[0])) 
      // 	       ||
      // 	       (sol[2] && (sol[0] + processing_time[0] > sol[1])) );
    }
    virtual PropagationOutcome propagate(const int changed_idx, const Event evt);
    virtual PropagationOutcome propagate();
    //virtual RewritingOutcome rewrite();
    //virtual void consolidate();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "<>="; }
    //@}
  };


//   template< int ARITY >
//   /**********************************************
//    *Tuple Constraint
//    **********************************************/
//   //  <x1 , ... , xn> == y
//   /// associate a value in y with each tuple in a given list
//   /// and ensures that the tuple corresponding to the value of y
//   /// is equal to <x1 , ... , xn>
//   class PredicateTuple : public GlobalConstraint {
 
//   public:

//     typedef Tuple< ARITY, int >  tuple;
//     typedef Tuple< 2, int > assignment;

//     /**@name Parameters*/
//     //@{
//     // // given the rank $i$ of a tuple in $table$, and a variable $j$, the id of this tuple in the set $support_of[j]$
//     // tuple<ARITY> *table2support;
    
//     // // given a variable $j$ and a value $k$ in the set $support_of[j]$, the id of the corresponding tuple in $table$
//     // int **support2table;

//     // for each variable and each value, the set of support. In fact, the rank of the support in $tabl$ is stored in the set.
//     // since support are partitioned between values, all set for a given variable share the same $index_$
//     ReversibleSet          *support_of[ARITY];
//     // values initialy in the domains (used for memory deallocation)
//     Vector<int>                 values[ARITY];

//     // the list of tuples
//     Vector< tuple >             table;
    
//     Vector< assignment >       pruned;

//     DomainDelta            value_delta[ARITY];

//     //@}

//     /**@name Constructors*/
//     //@{
//     PredicateTuple() : GlobalConstraint() { }
//     PredicateTuple(Vector< Variable >& scp);
//     PredicateTuple(std::vector< Variable >& scp);
//     virtual Constraint clone() { return Constraint(new PredicateTuple(scope)); }
//     virtual void initialise();
//     virtual int idempotent() { return 1;}
//     virtual int postponed() { return 1;}
//     virtual int pushed() { return 1;}
//     //virtual bool absorb_negation(const int var) { return true; }
//     virtual ~PredicateTuple();
//     //@}

//     /**@name Solving*/
//     //@{
//     virtual int check( const int* sol ) const ;
//     virtual PropagationOutcome propagate();
//     //virtual RewritingOutcome rewrite();
//     //@}
  
//     /**@name Miscellaneous*/
//     //@{  
//     virtual std::ostream& display(std::ostream&) const ;
//     virtual std::string name() const { return "bsum=k"; }
//     //@}
//   };


// template< int ARITY >
// PredicateTuple<ARITY>::PredicateTuple(Vector< Variable >& scp)
//   : GlobalConstraint(scp) { 
// }

// template< int ARITY >
// void PredicateTuple<ARITY>::initialise() {
//   ConstraintImplementation::initialise();

//   for(unsigned int i=0; i<scope.size; ++i)
//     trigger_on(_DOMAIN_, scope[i]);

//   //Vector<int> supports(0, max(256, table.size/ARITY/10));
//   //support_of = new ReversibleSet*[ARITY];
//   for(int i=0; i<ARITY; ++i) {
//     support_of[i] = new ReversibleSet[scope[i].get_initial_max() - scope[i].get_initial_min() + 1];
//     support_of[i] -= scope[i].get_initial_min();
//     values[i].initialise(0,scope[i].get_size());

//     int vnxt = scope[i].get_first(), val;
//     do {
//       val = vnxt;

//       if(values[i].empty()) {
//   	support_of[i][val].initialise(0, table.size, 8, false);
//       } else {
// 	support_of[i][val].initialise(support_of[i][values[i][0]], 8);
//       }

//       values[i].add(val);
      
//       vnxt = scope[i].next(val);
//     } while( val != vnxt );
//   }

//   for(unsigned int k=0; k<table.size; ++k) {
//     //std::cout << "add " << table[k] << std::endl;
//     for(int i=0; i<ARITY; ++i) {
//       support_of[i][table[k][i]].init_add(k);
//     }
//   }

//   GlobalConstraint::initialise();


//   // initialise the domain_delta
//   for(int i=0; i<ARITY; ++i) {
//     value_delta[i].initialise(scope[i]);
//   }


//   std::cout << changes << std::endl;
//   std::cout << events << std::endl;

//   // then we "propagate"
//   for(int xi=0; xi<ARITY; ++xi) {
//     Domain dom_xi(scope[xi]);
    
//     std::cout << "iterate over " << scope[xi].get_domain() << std::endl;
     
//     Domain::iterator xstop = dom_xi.begin();

//     int valj;
//     for(Domain::iterator xit = dom_xi.end(); --xit>=xstop; ) {
//       valj = dom_xi.get_value(xit); 

//       if(support_of[xi][valj].empty())  {
// 	scope[xi].remove(valj);
// 	pruned.add( assignment(xi, valj) );
//       }
//       std::cout << scope[xi] << " = " << valj << ":" ;
//       for(int k=0; k<support_of[xi][valj].size; ++k) {
// 	std::cout << " " << table[support_of[xi][valj][k]] ;
//       }
//       std::cout << std::endl;

//     }
//   }

//   std::cout << changes << std::endl;
//   std::cout << events << std::endl;





//   // initialise the domain_delta
//   for(int i=0; i<ARITY; ++i) {
//     DomainDelta::iterator xit = value_delta[i].begin();
//     DomainDelta::iterator xend = value_delta[i].end();

//     std::cout << scope[i] << " in " << scope[i].get_domain() << ": " << std::endl;
//     while(xit != xend) {
//       std::cout << "removed " << scope[i] << " = " << *xit << std::endl;
//       ++xit;
//     }
//     std::cout << std::endl;

//     value_delta[i].close();
//   }


//   for(int i=0; i<ARITY; ++i) {
//     DomainDelta::iterator xit = value_delta[i].begin();
//     DomainDelta::iterator xend = value_delta[i].end();

//     std::cout << scope[i] << " in " << scope[i].get_domain() << ": " << std::endl;
//     while(xit != xend) {
//       std::cout << "removed " << scope[i] << " = " << *xit << std::endl;
//       ++xit;
//     }
//     std::cout << std::endl;

//     value_delta[i].close();
//   }


//   scope[0].remove(13);
//   scope[0].remove(9);


//   scope[2].remove(3);
//   scope[2].remove(17);
//   scope[2].remove(8);
//   scope[2].remove(9);

//   for(int i=0; i<ARITY; ++i) {
//     DomainDelta::iterator xit = value_delta[i].begin();
//     DomainDelta::iterator xend = value_delta[i].end();

//     std::cout << scope[i] << " in " << scope[i].get_domain() << ": " << std::endl;
//     while(xit != xend) {
//       std::cout << "removed " << scope[i] << " = " << *xit << std::endl;
//       ++xit;
//     }
//     std::cout << std::endl;

//     value_delta[i].close();
//   }



//   std::cout << pruned << std::endl;

// }

// template< int ARITY >
// PredicateTuple<ARITY>::~PredicateTuple() 
// {
// #ifdef _DEBUG_MEMORY
//   std::cout << "c delete boolsumequal constraint" << std::endl;
// #endif
// }

// template< int ARITY >
// PropagationOutcome PredicateTuple<ARITY>::propagate() 
// {
//   PropagationOutcome wiped = CONSISTENT;

//   int xi, xj, vali, valj, sk, tk;
//   Domain::iterator xstop;
//   Domain::iterator xit;

//   while(!changes.empty()) {
//     xi = changes.pop();
//     xstop = value_delta[xi].begin();
//     for(xit = value_delta[xi].end(); --xit>=xstop; ) {
//       vali = *xit;
//       for(sk = support_of[xi][vali].size; sk--;) {
// 	tk = support_of[xi][vali][sk];
// 	for(xj = ARITY; xj--;) {
// 	  valj = table[tk][xj];
// 	  support_of[xj][valj].remove(tk);
// 	  if(support_of[xj][valj].empty()) {
// 	    if(FAILED(scope[xj].remove(valj))) { wiped = FAILURE(xj); goto FAIL; }
// 	    if(!changes.contain(xj)) changes.add(xj);
// 	  }
// 	}
//       }
//     }
//   }
  
//  FAIL: 
//   return wiped;
// }

// template< int ARITY >
// int PredicateTuple<ARITY>::check( const int* s ) const 
// {
//   // int i=scope.size, t=0;
//   // while(--i) t+=s[i];
//   // return total != t; 
//   return true;
// }

// template< int ARITY >
// std::ostream& PredicateTuple<ARITY>::display(std::ostream& os) const {
//   // os << "(" << scope[0]/*.get_var()*/ ;
//   // for(unsigned int i=1; i<scope.size; ++i) 
//   //   os << " + " << scope[i]/*.get_var()*/;
//   // os << ") == " << total ;
//   return os;
// }




  /**********************************************
   *Table Constraint
   **********************************************/
  //  <x1 , ... , xn> == y
  /// associate a value in y with each tuple in a given list
  /// and ensures that the tuple corresponding to the value of y
  /// is equal to <x1 , ... , xn>
  class ConstraintTable : public GlobalConstraint {
 
  public:

    /**@name Parameters*/
    //@{
    // the list of tuples
    Vector< const int* >          table;
    //@}

    /**@name Constructors*/
    //@{
    ConstraintTable() : GlobalConstraint() { }
    ConstraintTable(Vector< Variable >& scp);
    ConstraintTable(std::vector< Variable >& scp);
    virtual Constraint clone() { return Constraint(new ConstraintTable(scope)); }
    virtual void initialise();
    virtual void mark_domain();
    virtual int idempotent() { return 1;}
    virtual int postponed() { return 1;}
    virtual int pushed() { return 1;}
    virtual ~ConstraintTable();
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const ;
    //@}

    void add(const int* tuple);
  
    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "table"; }
    //@}
  };




  /**********************************************
   *GAC4 Constraint
   **********************************************/
  //  <x1 , ... , xn> == y
  /// associate a value in y with each tuple in a given list
  /// and ensures that the tuple corresponding to the value of y
  /// is equal to <x1 , ... , xn>
  class ConstraintGAC4 : public ConstraintTable {
 
  public:

    //typedef Tuple< 2, int > assignment;

    /**@name Parameters*/
    //@{
    // for each variable and each value, the set of support. In fact, the rank of the support in $tabl$ is stored in the set.
    // since support are partitioned between values, all set for a given variable share the same $index_$
#ifdef _OS_PRUNING
    ReversibleSet         ***support_of;
#else
    ReversibleSet          **support_of;
#endif
    // values initialy in the domains (used for memory deallocation)
    Vector<int>                 *values;

    // the list of tuples
    //Vector< const int* >          table;
    
    // Handlers for domain-delta
    DomainDelta            *value_delta;


    // 
    // 
    
#ifdef _OS_PRUNING
    ReversibleSet        *support_lists;
    int                      **list_map;
    int                        *var_map;
    int                        *val_map;
    IntStack               pruned_lists;
    Vector<int>           *list_pruning;
#endif

    //@}

#ifdef _DEBUG_TABLE
    int* init_support_size;
#endif

    /**@name Constructors*/
    //@{
    ConstraintGAC4() : ConstraintTable() { }
    ConstraintGAC4(Vector< Variable >& scp);
    ConstraintGAC4(std::vector< Variable >& scp);
    virtual Constraint clone() { return Constraint(new ConstraintGAC4(scope)); }
    virtual void initialise();
    //virtual int idempotent() { return 1;}
    //virtual int postponed() { return 1;}
    //virtual int pushed() { return 1;}
    //virtual bool absorb_negation(const int var) { return true; }
    virtual ~ConstraintGAC4();
    //@}



    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const ;
    virtual PropagationOutcome propagate();
    //virtual RewritingOutcome rewrite();
    //@}

#ifdef _DEBUG_TABLE
    bool check_integrity( ) const ;
#endif
  
    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display_supports(std::ostream&) const ;
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "bsum=k"; }
    //@}
  };


  /**********************************************
   *GAC2001 Constraint
   **********************************************/
  //  <x1 , ... , xn> == y
  /// associate a value in y with each tuple in a given list
  /// and ensures that the tuple corresponding to the value of y
  /// is equal to <x1 , ... , xn>
  class ConstraintGAC2001 : public ConstraintTable {
 
  public:

    //typedef Tuple< 2, int > assignment;

    /**@name Parameters*/
    //@{
    ReversibleNum<int> **firstSupport;
    Vector< const int* > **supportList;
    int *order;
    //BitSet *D_X;
    const int** supports_X;
    int* themins;
    //Vector< const int* > tuples;
    //bool isClone;
    //@}

    /**@name Constructors*/
    //@{
    ConstraintGAC2001() : ConstraintTable() { }
    ConstraintGAC2001(Vector< Variable >& scp);
    ConstraintGAC2001(std::vector< Variable >& scp);
    virtual Constraint clone() { return Constraint(new ConstraintGAC2001(scope)); }
    virtual void initialise();
    //virtual int idempotent() { return 1;}
    //virtual int postponed() { return 1;}
    //virtual int pushed() { return 1;}
    //virtual bool absorb_negation(const int var) { return true; }
    virtual ~ConstraintGAC2001();
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const ;
    virtual PropagationOutcome propagate();
    //virtual RewritingOutcome rewrite();
    //@}

    //void add(const int* tuple);
  
    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display_supports(std::ostream&) const ;
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "gac2001"; }
    //@}
  };



  /**********************************************
   *GAC3 Constraint
   **********************************************/
   /*! \class ConstraintGAC3Valid
    \brief  implementation of a n-ary extensional constraint

    The conflicts are encoded using
    a N-dimensional flatened matrix:
    an entry of the matrix equal to 0
    denotes an allowed combination (support), whilst
    1 denotes a forbidden combination (conflict).
  */
  class ConstraintGAC3 : public ConstraintTable {
 
  private:
    /**@name Parameters*/
    //@{ 
    int *var_sizes; // NOTE: this has arity-1 elements
    BitSet matrix;
    //bool isClone;
    //@}

 public:
 
    /**@name Constructors*/
    //@{
    ConstraintGAC3() : ConstraintTable() {}
    ConstraintGAC3(Vector< Variable >& scp);
    ConstraintGAC3(std::vector< Variable >& scp);
    virtual Constraint clone() { return Constraint(new ConstraintGAC3(scope)); }
    virtual void initialise();
    //virtual int idempotent() { return 1;}
    //virtual int postponed() { return 1;}
    //virtual int pushed() { return 1;}
    //virtual bool absorb_negation(const int var) { return true; }
    virtual ~ConstraintGAC3();
    //@}

    int getpos(const int *vals) const;
    int getpos() const;
    bool isValid(const int *tuple) const;

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const ;
    //virtual PropagationOutcome propagate();
    //virtual RewritingOutcome rewrite();
    //@}

    //void add(const int* tuple);
  
    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display_supports(std::ostream&) const ;
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "gac3"; }
    //@}
  };




  // /**********************************************
  //  * BoolSum Equal Constraint
  //  **********************************************/
  // //  b1 + ... + bn = k
  // /// Constraint on the value of the sum of a set of variables.
  // class ConstraintBoolSumEqual : public GlobalConstraint {
 
  // public:
  //   /**@name Parameters*/
  //   //@{
  //   // last "gap" from the tightest bound to the total, used to know what cause pruning/failure
  //   //int gap;
  //   int total;
  //   int lower_bound;
  //   int upper_bound;
  //   ReversibleNum<int> min_;
  //   ReversibleNum<int> max_;
  //   // used to store the explanation when "get_reason_for()" is called
  //   Vector<Literal> explanation;
  //   //@}

  //   /**@name Constructors*/
  //   //@{
  //   ConstraintBoolSumEqual() : GlobalConstraint() { priority = 1; total = 0; }
  //   ConstraintBoolSumEqual(Vector< Variable >& scp, const int t);
  //   ConstraintBoolSumEqual(std::vector< Variable >& scp, const int t);
  //   virtual Constraint clone() { return Constraint(new ConstraintBoolSumEqual(scope, total)); }
  //   virtual void initialise();
  //   virtual void mark_domain();
  //   virtual int idempotent() { return 1;}
  //   virtual int postponed() { return 1;}
  //   virtual int pushed() { return 1;}
  //   //virtual bool absorb_negation(const int var) { return true; }
  //   virtual ~ConstraintBoolSumEqual();
  //   //@}


  //   // virtual Explanation::iterator begin(Atom a);// { return NULL; }
  //   // virtual Explanation::iterator end  (Atom a);// { return NULL; } 
  //   virtual iterator get_reason_for(const Atom a, const int lvl, iterator& end);
  //   virtual void initialise_activity(double *lact, double *vact, double norm);

  //   /**@name Solving*/
  //   //@{
  //   virtual int check( const int* sol ) const ;
  //   virtual PropagationOutcome propagate();
  //   //virtual PropagationOutcome propagate_and_explain(Vector<Explanation*>);
  //   //virtual RewritingOutcome rewrite();
  //   //@}
  
  //   /**@name Miscellaneous*/
  //   //@{  
  //   virtual std::ostream& display(std::ostream&) const ;
  //   virtual std::string name() const { return "bsum=k"; }
  //   //@}
  // };


  /**********************************************
   * BoolSum Interval Constraint
   **********************************************/
  //  lb <= b1 + ... + bn <= ub
  /// Constraint on the value of the sum of a set of variables.
  class ConstraintBoolSumInterval : public GlobalConstraint {
 
  public:
    /**@name Parameters*/
    //@{
    int lower_bound;
    int upper_bound;
    ReversibleNum<int> min_;
    ReversibleNum<int> max_;
    // used to store the explanation when "get_reason_for()" is called
    Vector<Literal> explanation;
    bool init_prop;
    //@}

    /**@name Constructors*/
    //@{
    ConstraintBoolSumInterval() : GlobalConstraint() { priority = 1; }
    ConstraintBoolSumInterval(Vector< Variable >& scp, const int l, const int u);
    ConstraintBoolSumInterval(std::vector< Variable >& scp, const int l, const int u);
    virtual Constraint clone() { return Constraint(new ConstraintBoolSumInterval(scope, lower_bound, upper_bound)); }
    virtual void initialise();
    virtual void mark_domain();
    virtual bool explained() { return true; }
    virtual int idempotent() { return 1;}
    virtual int postponed() { return 1;}
    virtual int pushed() { return 1;}
    //virtual bool absorb_negation(const int var) { return true; }
    virtual ~ConstraintBoolSumInterval();
#ifdef _CBSI_WC
    void weight_conflict(double unit, Vector<double>& weights);//  {
#endif
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const ;
    virtual PropagationOutcome propagate();
    virtual iterator get_reason_for(const Atom a, const int lvl, iterator& end);
    virtual void initialise_activity(double *lact, double *vact, double norm);
    //@}
  
    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "bsum[]"; }
    //@}
  };


  /**********************************************
   * WeightedBoolSum Interval Constraint
   **********************************************/
  //  lb <= a1 * b1 + ... + an * bn <= ub
  /// Constraint on the value of the sum of a set of variables.
  class ConstraintWeightedBoolSumInterval : public GlobalConstraint {
 
  public:
    /**@name Parameters*/
    //@{
    // Lower bound of the linear expression
    int lower_bound;

    // Upper bound of the linear expression
    int upper_bound;

    // coefficients should be in decreasing order
    Vector< int > weight;
    // from index 0 to wpos (not included), the coefficients are all 1s
    int wpos;
    // from index wpos to wneg (not included), the coefficients are all >0
    int wneg;
    // from index wneg to size (not included), the coefficients are all <0

    // utils for the propagation
    BoolDomain *domains;
    ReversibleNum<int> parity;
    //ReversibleIntStack unknown_parity;
    ReversibleSet unknown_parity;


    // used to store the explanation when "get_reason_for()" is called
    // variable assignments that increase the value of the sum
    Vector<Literal> positive_contributors;
    // variable assignments that decrease the value of the sum
    Vector<Literal> negative_contributors;
    //@}

    /**@name Constructors*/
    //@{
    ConstraintWeightedBoolSumInterval() : GlobalConstraint() { priority = 1; }

    ConstraintWeightedBoolSumInterval(Vector< Variable >& scp,
				      const int L=0, const int U=0);
    ConstraintWeightedBoolSumInterval(Vector< Variable >& scp,
				      Vector< int >& coefs,
				      const int L=0, const int U=0);
    ConstraintWeightedBoolSumInterval(std::vector< Variable >& scp,
				      std::vector< int >& coefs,
				      const int L=0, const int U=0);
    //ConstraintWeightedBoolSumInterval(Vector< Variable >& scp, const int l, const int u);
    //ConstraintWeightedBoolSumInterval(std::vector< Variable >& scp, const int l, const int u);
    virtual Constraint clone() { return Constraint(new ConstraintWeightedBoolSumInterval(scope, lower_bound, upper_bound)); }
    virtual void initialise();
    virtual void mark_domain();
    virtual bool explained() { return true; }
    virtual int idempotent() { return 1;}
    virtual int postponed() { return 1;}
    virtual int pushed() { return 1;}
    //virtual bool absorb_negation(const int var) { return true; }
    virtual ~ConstraintWeightedBoolSumInterval();
#ifdef _CWBSI_WC
    void weight_conflict(double unit, Vector<double>& weights);//  {
#endif
#ifdef _CWBSI_WC_ALT
    void weight_conflict(double unit, Vector<double>& weights);//  {
#endif

    //@}


    virtual iterator get_reason_for(const Atom a, const int lvl, iterator& end);
    //virtual void initialise_activity(double *lact, double *vact, double norm);

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const ;
    virtual PropagationOutcome propagate();
    //@}
  
    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "wbsum[]"; }
    //@}
  };



  /**********************************************
   * WeightedBoolSum Interval Constraint
   **********************************************/
  //  lb <= a1 * b1 + ... + an * bn <= ub
  /// Constraint on the value of the sum of a set of variables.
  class ConstraintIncrementalWeightedBoolSumInterval : public GlobalConstraint {
 
  public:
    /**@name Parameters*/
    //@{
    // Lower bound of the linear expression
    int lower_bound;

    // Upper bound of the linear expression
    int upper_bound;

    // coefficients should be in decreasing order
    Vector< int > weight;
 
    bool init_prop;

    // utils for the propagation
    BoolDomain *domains;
    ReversibleNum<int> bound_[2];
    //ReversibleNum<int> max_;
    // points to the unassigned variable with maximum coefficient (in absolute value)
    ReversibleNum<int> index_;
    //ReversibleNum<int> max_index;

    // used to store the explanation when "get_reason_for()" is called
    Vector<Literal> explanation;
    //@}

    /**@name Constructors*/
    //@{
    ConstraintIncrementalWeightedBoolSumInterval() : GlobalConstraint() { priority = 1; }

    ConstraintIncrementalWeightedBoolSumInterval(Vector< Variable >& scp,
				      const int L=0, const int U=0);
    ConstraintIncrementalWeightedBoolSumInterval(Vector< Variable >& scp,
				      Vector< int >& coefs,
				      const int L=0, const int U=0);
    ConstraintIncrementalWeightedBoolSumInterval(std::vector< Variable >& scp,
				      std::vector< int >& coefs,
				      const int L=0, const int U=0);
    //ConstraintIncrementalWeightedBoolSumInterval(Vector< Variable >& scp, const int l, const int u);
    //ConstraintIncrementalWeightedBoolSumInterval(std::vector< Variable >& scp, const int l, const int u);
    virtual Constraint clone() { return Constraint(new ConstraintIncrementalWeightedBoolSumInterval(scope, lower_bound, upper_bound)); }
    virtual void initialise();
    virtual void mark_domain();
    virtual bool explained() { return true; }
    virtual int idempotent() { return 1;}
    virtual int postponed() { return 1;}
    virtual int pushed() { return 1;}
    //virtual bool absorb_negation(const int var) { return true; }
    virtual ~ConstraintIncrementalWeightedBoolSumInterval();
    //@}

#ifdef _CIWBSI_WC
    void weight_conflict(double unit, Vector<double>& weights);//  {
#endif
#ifdef _CIWBSI_WC_ALT
    void weight_conflict(double unit, Vector<double>& weights);//  {
#endif
    virtual iterator get_reason_for(const Atom a, const int lvl, iterator& end);
    virtual void initialise_activity(double *lact, double *vact, double norm);

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const ;
    virtual PropagationOutcome propagate();
    //@}
  
    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "wbsum[]"; }
    //@}
  };



  /**********************************************
   * WeightedBoolSum Predicate
   **********************************************/
  //  lb <= a1 * b1 + ... + an * bn <= ub
  /// Constraint on the value of the sum of a set of variables.
  class PredicateWeightedBoolSum : public GlobalConstraint {
 
  public:
    /**@name Parameters*/
    //@{

    // coefficients should be in decreasing order
    Vector< int > weight;

    // offset (constant term)
    int offset;
 

    // utils for the propagation
    BoolDomain *domains;
    ReversibleNum<int> bound_[2];
    //ReversibleNum<int> max_;
    // points to the unassigned variable with maximum coefficient (in absolute value)
    ReversibleNum<int> index_;
    //ReversibleNum<int> max_index;

    bool init_prop;

    // used to store the explanation when "get_reason_for()" is called
    Vector<Literal> explanation;
    //@}

    /**@name Constructors*/
    //@{
    PredicateWeightedBoolSum() : GlobalConstraint() { priority = 1; }

    PredicateWeightedBoolSum(Vector< Variable >& scp, const int o=0);
    PredicateWeightedBoolSum(Vector< Variable >& scp,
			     Vector< int >& coefs, const int o=0);
    PredicateWeightedBoolSum(std::vector< Variable >& scp,
			     std::vector< int >& coefs, const int o=0);
    //PredicateWeightedBoolSum(Vector< Variable >& scp, const int l, const int u);
    //PredicateWeightedBoolSum(std::vector< Variable >& scp, const int l, const int u);
    virtual Constraint clone() { return Constraint(new PredicateWeightedBoolSum(scope, weight, offset)); }
    virtual void initialise();
    virtual void mark_domain();
    virtual bool explained() { return true; }
    virtual int idempotent() { return 1;}
    virtual int postponed() { return 1;}
    virtual int pushed() { return 1;}
    //virtual bool absorb_negation(const int var) { return true; }
    virtual ~PredicateWeightedBoolSum();
    //@}

    virtual void initialise_activity(double *lact, double *vact, double norm);
    virtual iterator get_reason_for(const Atom a, const int lvl, iterator& end);
#ifdef _PWBS_WC
    void weight_conflict(double unit, Vector<double>& weights);//  {
#endif
#ifdef _PWBS_WC_ALT
    void weight_conflict(double unit, Vector<double>& weights);//  {
#endif

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const ;
    virtual PropagationOutcome propagate();
    //@}
  
    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "wbsum[]"; }
    //@}
  };


  /**********************************************
   * BoolSum  Predicate
   **********************************************/
  //  b1 + ... + bn-1 = xn
  /// predicate on the value of the sum of a set of variables.
  class PredicateBoolSum : public GlobalConstraint {
 
  public:
    /**@name Parameters*/
    //@{

    // offset (constant term)
    int offset;
 
    int lb;
    int ub;
    ReversibleNum<int> min_;
    ReversibleNum<int> max_;

    Vector< Literal > explanation;
    bool init_prop;
    //@}

    /**@name Constructors*/
    //@{
    PredicateBoolSum() : GlobalConstraint() { priority = 1; }
    PredicateBoolSum(Vector< Variable >& scp, const int o=0);
    PredicateBoolSum(std::vector< Variable >& scp, const int o=0);
    PredicateBoolSum(Vector< Variable >& scp, Variable tot, const int o=0);
    PredicateBoolSum(std::vector< Variable >& scp, Variable tot, const int o=0);
    virtual Constraint clone() { return Constraint(new PredicateBoolSum(scope,offset)); }
    virtual void initialise();
    virtual void mark_domain();
    virtual bool explained() { return true; }
    virtual int idempotent() { return 1;}
    virtual int postponed() { return 1;}
    virtual int pushed() { return 1;}
    //virtual bool absorb_negation(const int var) { return true; }
    virtual ~PredicateBoolSum();
    //@}

#ifdef _PBS_WC
    void weight_conflict(double unit, Vector<double>& weights);//  {
#endif
    virtual iterator get_reason_for(const Atom a, const int lvl, iterator& end);

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const ;
    virtual PropagationOutcome propagate();
    //@}
  
    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "bsum="; }
    //@}
  };


  /**********************************************
   * WeightedSum Constraint
   **********************************************/
  /*! \class ConstraintWeightedSum
    \brief  Constraint on a sum of variables (a1 * x1 + ... + an * xn = total)
  */
  class PredicateWeightedSum : public GlobalConstraint {

  public:
    /**@name Parameters*/
    //@{ 
    // Lower bound of the linear expression
    int lower_bound;

    // Upper bound of the linear expression
    int upper_bound;

    Vector< int > weight;
    // from index 0 to wpos (not included), the coefficients are all 1s
    int wpos;
    // from index wpos to wneg (not included), the coefficients are all >0
    int wneg;
    // from index wneg to size (not included), the coefficients are all <0

    // utils for the propagation
    int *lo_bound;
    int *up_bound;
    int *span;
    ReversibleNum<int> parity;
    //ReversibleIntStack unknown_parity;
    ReversibleSet unknown_parity;
    
    //@}

    /**@name Constructors*/
    //@{
    PredicateWeightedSum(Vector< Variable >& scp,
			 const int L=0, const int U=0);
    PredicateWeightedSum(Vector< Variable >& scp,
			 Vector< int >& coefs,
			 const int L=0, const int U=0);
    PredicateWeightedSum(std::vector< Variable >& scp,
			 std::vector< int >& coefs,
			 const int L=0, const int U=0);
    virtual ~PredicateWeightedSum();
    virtual Constraint clone() { return Constraint(new PredicateWeightedSum(scope, weight)); }
    virtual bool rewritable() { return true; }
    virtual int idempotent() { return 1;}
    virtual int postponed() { return 1;}
    virtual int pushed() { return 1;}
    virtual void initialise();
    virtual void mark_domain();
#ifdef _PWS_WC
    void weight_conflict(double unit, Vector<double>& weights);//  {
#endif
#ifdef _PWS_WC_ALT
    void weight_conflict(double unit, Vector<double>& weights);//  {
#endif
//@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const ;
    virtual PropagationOutcome propagate();
    virtual RewritingOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "sum="; }
    //@}
  };



  /**********************************************
   * WeightedSum Constraint
   **********************************************/
  /*! \class ConstraintWeightedSum
    \brief  Constraint on a sum of variables (a1 * x1 + ... + an * xn = total)
  */
  class ConstraintOrderedSum : public GlobalConstraint {

  public:
    /**@name Parameters*/
    //@{ 
    // Lower bound of the linear expression
    int lower_bound;

    // Upper bound of the linear expression
    int upper_bound;

    // utils for the propagation
    int *lo_bound;
    int *up_bound;
		
    int *offset_min;
    int *offset_max;
    //@}

    /**@name Constructors*/
    //@{
    ConstraintOrderedSum(Vector< Variable >& scp,
			 const int L=0, const int U=0);
    ConstraintOrderedSum(std::vector< Variable >& scp,
			 const int L=0, const int U=0);
    virtual ~ConstraintOrderedSum();
    virtual Constraint clone() { return Constraint(new ConstraintOrderedSum(scope, lower_bound, upper_bound)); }
    virtual void initialise();
		virtual bool rewritable() { return false; }
    virtual int idempotent() { return 1;}
    virtual int postponed() { return 1;}
    virtual int pushed() { return 1;}
    //virtual void mark_domain();
		//@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const ;
    virtual PropagationOutcome propagate();
    //virtual RewritingOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "osum="; }
    //@}
  };





  /**********************************************
   * Parity Constraint
   **********************************************/
  /*! \class ConstraintParity
    \brief  Constraint on a sum modulo 2 of variables (x1 + ... + xn = p)
    with p in {0,1}
  */
  class ConstraintParity : public GlobalConstraint {

  public:
    /**@name Parameters*/
    //@{ 
    // the wanted parity 
    int target_parity;

    //
    Vector< Literal > explanation;
    //@}

    /**@name Constructors*/
    //@{
    ConstraintParity(Vector< Variable >& scp, const int p=0);
    virtual ~ConstraintParity();
    virtual Constraint clone() { return Constraint(new ConstraintParity(scope, target_parity)); }
    //virtual bool rewritable() { return true; }
    virtual bool explained() { return true; }
    virtual int idempotent() { return 1;}
    virtual int postponed() { return 1;}
    virtual int pushed() { return 1;}
    virtual void initialise();
    //virtual void mark_domain();
    //@}

    virtual iterator get_reason_for(const Atom a, const int lvl, iterator& end);

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const ;
    virtual PropagationOutcome propagate();
    //virtual RewritingOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "parity"; }
    //@}
  };



  /**********************************************
   * AtMostSeqCard Constraint
   **********************************************/
  //  
  /// 
  class ConstraintMultiAtMostSeqCard : public GlobalConstraint {
 
  public:
    /**@name Parameters*/
    //@{
    int _k;
    int _d;
    int *_p;
    int *_q;

    int *wl; //[arity+2*_q];
    int *wr; //[arity+2*_q];
    int **occurrences; //[k*arity+_p+1];
    int **cardinality; //[k*2*_q];
    
    int *lcumulated; //[arity+1];
    int *rcumulated; //[arity+1];
    VarArray reverse;
    // used to store the explanation when "get_reason_for()" is called
    Vector<Literal> explanation;
    //we need this for the explanation to check if the maximum cardinality of all subsequences at position i is equal to p.
    Vector< bool> max_equal_to_p ;
	//Vector<int> sequence_image;
	Vector<int> left_right_intersection;
    //@}

    /**@name Constructors*/
    //@{
    ConstraintMultiAtMostSeqCard();
    ConstraintMultiAtMostSeqCard(Vector< Variable >& scp, const int k, const int d, const int* p, const int* q);
    ConstraintMultiAtMostSeqCard(std::vector< Variable >& scp, const int k, const int d, const int* p, const int* q);
    void initialise_struct(const int k=0, const int d=0, const int* p=NULL, const int* q=NULL);
    virtual Constraint clone() { return Constraint(new ConstraintMultiAtMostSeqCard(scope, _k, _d, _p, _q)); }
    virtual void initialise();
    virtual void mark_domain();
    virtual bool explained() { return true; }
    virtual int idempotent() { return 1;}
    virtual int postponed() { return 1;}
    virtual int pushed() { return 1;}
    virtual ~ConstraintMultiAtMostSeqCard();
    //@}
    
    virtual iterator get_reason_for(const Atom a, const int lvl, iterator& end);
    void greedy_assign_for_explanation(Vector<Variable>& X, int __size, int __rank);
    void set_max_equal_to_p_at_rank(int __rank, int __size,  Vector<Variable>& X);
    /**@name Solving*/
    //@{
    bool greedy_assign(int *w, int *cumulated, Vector<Variable>& X) ;
    virtual int check( const int* sol ) const ;
    virtual PropagationOutcome propagate();
    //@}
  
    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "bsum=k"; }
    //@}
  };

  /**********************************************
   * AtMostSeqCard Constraint
   * The same constraint as ConstraintMultiAtMostSeqCard with the exeption of generating default nogoods without reductions
   **********************************************/
  //
  ///
  class ConstraintNaiveMultiAtMostSeqCard : public ConstraintMultiAtMostSeqCard {

  public:
	  /**@name Constructors*/
	  //@{
	  //ConstraintNaiveMultiAtMostSeqCard();
	  ConstraintNaiveMultiAtMostSeqCard(Vector< Variable >& scp, const int k, const int d, const int* p, const int* q);
	  //ConstraintNaiveMultiAtMostSeqCard(std::vector< Variable >& scp, const int k, const int d, const int* p, const int* q);
	  virtual iterator get_reason_for(const Atom a, const int lvl, iterator& end);
  };



//Simplified Explanation
  class ConstraintSimplifiedExplanationMultiAtMostSeqCard : public ConstraintMultiAtMostSeqCard {

  public:
	  /**@name Constructors*/
	  //@{
	  //ConstraintNaiveMultiAtMostSeqCard();
	  ConstraintSimplifiedExplanationMultiAtMostSeqCard(Vector< Variable >& scp, const int k, const int d, const int* p, const int* q);
	  //ConstraintNaiveMultiAtMostSeqCard(std::vector< Variable >& scp, const int k, const int d, const int* p, const int* q);
	  virtual iterator get_reason_for(const Atom a, const int lvl, iterator& end);
  };


  //lazyleft Explanation
    class ConstraintLeftExplanationMultiAtMostSeqCard : public ConstraintMultiAtMostSeqCard {

    public:
  	  /**@name Constructors*/
  	  //@{
  	  //ConstraintNaiveMultiAtMostSeqCard();
	  ConstraintLeftExplanationMultiAtMostSeqCard(Vector< Variable >& scp, const int k, const int d, const int* p, const int* q);
  	  //ConstraintNaiveMultiAtMostSeqCard(std::vector< Variable >& scp, const int k, const int d, const int* p, const int* q);
  	  virtual iterator get_reason_for(const Atom a, const int lvl, iterator& end);
  	  void simple_greedy_assign_for_explanation(Vector<Variable>& X, int __rank, int index_a);
    };


   /**********************************************
    * VertexCover Predicate
    **********************************************/
   /*! \class PredicateVertexCover
     \brief  Predicate on the size of the vertex cover {X} of a graph G
   */
   class PredicateVertexCover : public GlobalConstraint {

   public:
     /**@name Parameters*/
     //@{ 
	 Graph _G;
     //@}

     /**@name Constructors*/
     //@{
     PredicateVertexCover(Vector< Variable >& scp, Graph& g);
     virtual ~PredicateVertexCover();
     virtual Constraint clone() { return Constraint(new PredicateVertexCover(scope,_G)); }
     virtual int idempotent() { return 1;}
     virtual int postponed() { return 1;}
     virtual int pushed() { return 1;}
     virtual void initialise();
     //virtual void mark_domain();
     //@}

     /**@name Solving*/
     //@{
     virtual int check( const int* sol ) const ;
     virtual PropagationOutcome propagate();
     //@}

     /**@name Miscellaneous*/
     //@{  
     virtual std::ostream& display(std::ostream&) const ;
     virtual std::string name() const { return "|vertex cover|="; }
     //@}
 };
 
 
  /**********************************************
   * Footrule Predicate
   **********************************************/
  /*! \class PredicateFootrule
    \brief  Predicate on the size of the vertex cover {X} of a graph G
  */
  class PredicateFootrule : public GlobalConstraint {

  public:
    /**@name Parameters*/
    //@{ 
  	int  N;
  	int  uncorrelated_distance;
  	int* values;
	
	
	bool init_prop;
	
	// structure for the DP
	
	// state <=> value of the sum [N+1]
	IntStack* state;
	
	// transition <=> manhattan distance value [N]
	IntStack* distance;
	
	
	// util bitset/intstacks
	BitSet   util_bitset;
	IntStack util_stack;
	
	
	int **twitness[2];
	int **switness;
	
	
#ifdef _CHECKED_MODE
	IntStack *domains;/////.....
#endif
	
	
	// for each transition
	
    //@}

    /**@name Constructors*/
    //@{
    PredicateFootrule(Vector< Variable >& scp);
    virtual ~PredicateFootrule();
    virtual Constraint clone() { return Constraint(new PredicateFootrule(scope)); }
    virtual int idempotent() { return 1;}
    virtual int postponed() { return 1;}
    virtual int pushed() { return 1;}
    virtual void initialise();
    //virtual void mark_domain();
    //@}

    /**@name Solving*/
    //@{
	int max_md(const int n, const int k) const;
	PropagationOutcome prune_from_transitions(const int i);
	PropagationOutcome compute_DP_from_scratch();
	PropagationOutcome initial_propagate();
    virtual int check( const int* sol ) const ;
    virtual PropagationOutcome propagate();
    //@}

    /**@name Miscellaneous*/
    //@{  
	void print_automaton() const;
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "footrule="; }
    //@}
};
 

  /**********************************************
   * Min Predicate
   **********************************************/
  /*! \class PredicateMin
    \brief  Predicate on the mininum of a set of variables z = (x1 , ... , xn)
  */
  class PredicateMin : public GlobalConstraint {

  public:
    /**@name Parameters*/
    //@{ 
    int last_min;
    ReversibleSet candidates;
    //@}

    /**@name Constructors*/
    //@{
    PredicateMin(Vector< Variable >& scp);
    virtual ~PredicateMin();
    virtual Constraint clone() { return Constraint(new PredicateMin(scope)); }
    virtual int idempotent() { return 1;}
    virtual int postponed() { return 1;}
    virtual int pushed() { return 1;}
    virtual void initialise();
    virtual void mark_domain();
    //@}

    /**@name Solving*/
    //@{
    void react_to(PropagationOutcome& wiped, const int changed_idx, const Event evt);
    virtual int check( const int* sol ) const ;
    virtual PropagationOutcome propagate();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "min="; }
    //@}
  };


  /**********************************************
   * Max Predicate
   *********************************************/
  /*! \class PredicateMax
    \brief  Predicate on the maxinum of a set of variables z = (x1 , ... , xn)
  */
  class PredicateMax : public GlobalConstraint {

  public:
    /**@name Parameters*/
    //@{ 
    int last_max;
    ReversibleSet candidates;
    //@}

    /**@name Constructors*/
    //@{
    PredicateMax(Vector< Variable >& scp);
    virtual ~PredicateMax();
    virtual Constraint clone() { return Constraint(new PredicateMax(scope)); }
    virtual int idempotent() { return 1;}
    virtual int postponed() { return 1;}
    virtual int pushed() { return 1;}
    virtual void initialise();
    virtual void mark_domain();
    //@}

    /**@name Solving*/
    //@{
    void react_to(PropagationOutcome& wiped, const int changed_idx, const Event evt);
    virtual int check( const int* sol ) const ;
    virtual PropagationOutcome propagate();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "max="; }
    //@}
  };


  /**********************************************
   * Element Predicate
   **********************************************/
  /*! \class PredicateElement
    \brief  
  */
  class PredicateElement : public GlobalConstraint {

  public:
    /**@name Parameters*/
    //@{ 
    BitSet aux_dom;
    int offset;
    //@}

    /**@name Constructors*/
    //@{
    PredicateElement(Vector< Variable >& scp, const int o=0);
    PredicateElement(std::vector< Variable >& scp, const int o=0);
    virtual ~PredicateElement();
    virtual Constraint clone() { return Constraint(new PredicateElement(scope, offset)); }
    virtual int idempotent() { return 1;}
    virtual int postponed() { return 1;}
    virtual int pushed() { return 1;}
    virtual void initialise();
    virtual void mark_domain();
#ifdef _ELT_WC
    void weight_conflict(double unit, Vector<double>& weights);//  {
#endif
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const ;
    virtual PropagationOutcome propagate();
    //virtual RewritingOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "[x]="; }
    //@}
  };


  /**********************************************
   * BoolElement Predicate
   **********************************************/
  /*! \class PredicateBoolElement
    \brief  
  */
  class PredicateBoolElement : public GlobalConstraint {

  public:
    /**@name Parameters*/
    //@{ 
    int aux_dom;
    int offset;
    //@}

    /**@name Constructors*/
    //@{
    PredicateBoolElement(Vector< Variable >& scp, const int o=0);
    PredicateBoolElement(std::vector< Variable >& scp, const int o=0);
    virtual ~PredicateBoolElement();
    virtual Constraint clone() { return Constraint(new PredicateBoolElement(scope, offset)); }
    virtual int idempotent() { return 1;}
    virtual int postponed() { return 1;}
    virtual int pushed() { return 1;}
    virtual void initialise();
    virtual void mark_domain();
    //void weight_conflict(double unit, Vector<double>& weights);//  {
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const ;
    virtual PropagationOutcome propagate();
    //virtual RewritingOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "[x]="; }
    //@}
  };



  /***********************************************
   * AllDifferent Constraint (forward checking).
   ***********************************************/
  /*! \class ConstraintCliqueNotEqual
    \brief  Clique of NotEqual Constraint.
  */
  class ConstraintCliqueNotEqual : public GlobalConstraint {

  public:

    //Vector<int> assigned;
    int culprit;
    
    /**@name Constructors*/
    //@{
    ConstraintCliqueNotEqual() : GlobalConstraint() { priority = 2; }
    ConstraintCliqueNotEqual(Vector< Variable >& scp);
    ConstraintCliqueNotEqual(std::vector< Variable >& scp);
    ConstraintCliqueNotEqual(Variable* scp, const int n);
    virtual Constraint clone() { return Constraint(new ConstraintCliqueNotEqual(scope)// , type
						   ); }
    virtual void initialise();
    virtual void mark_domain();
    virtual int idempotent() { return 1;}
    virtual int postponed() { return 1;}
    virtual int pushed() { return 1;}
    virtual ~ConstraintCliqueNotEqual();
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const ;
    virtual PropagationOutcome propagate();
    //virtual PropagationOutcome propagate(const int changed_idx, const Event evt);
    //virtual RewritingOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "{=/=}"; }
    //@}
#ifdef _CNE_WC
    void weight_conflict(double unit, Vector<double>& weights);//  {
#endif
    //   std::cout << std::endl;
    //   for(int i=0; i<scope.size; ++i) {
    // 	std::cout << scope[i].get_domain() << std::endl;
    // 	weights[scope[i].id()] += unit;
    //   }
    // }
    
  };



  /***********************************************
   * All Different Constraint (bounds consistency).
   ***********************************************/
  typedef struct {
    int min, max;		// start, end of interval
    int minrank, maxrank; // rank of min & max in bounds[] of an adcsp
  } AD_Interval;
  
  /*! \class ConstraintAllDiff
    \brief  AllDifferent Constraint.

    Constraint of difference on a set of variables.
    Ensures that a set of variables are assigned distinct 
    values (only Bounds Consistency is implemented)
    The code is from Lopez-Ortiz, Quimper, Tromp and van Beek
  */
  class ConstraintAllDiff : public GlobalConstraint {
    
  private:
    /**@name Parameters*/
    //@{  
    //int *level;
    int lastLevel;
    int *t;		// tree links
    int *d;		// diffs between critical capacities
    int *h;		// hall interval links
    AD_Interval *iv;
    AD_Interval **minsorted;
    AD_Interval **maxsorted;
    int *bounds;  // bounds[1..nb] hold set of min & max in the niv intervals
                  // while bounds[0] and bounds[nb+1] allow sentinels
    int nb;

    
    int expl_note;
    // used to store 
    //  -- which interval order was used to find the Hall interval (least significant bit)
    //  -- at which rank was is discovered (rest)
    // When computing an explanation, suppose that it was when exploring the maxsorted intervals an at rank r
    //  -> then maxsorted[r]->max is the max of the Hall interval AND it is a HI of the subset of variables whose rank is <= r
    //     we then use the following array:
    //Vector<int> count_bound;
    std::vector<int> other_bounds;
    // let m = maxsorted[r]->max
    // for each pos s in the interval [0, m] we store how many times an interval finishing before or at m starts at s
    // then, going from m toward 0, we compute how many of these intervals start at a pos s or before (call this b[s])
    // IF s+b[s]>m+1 THEN [s, m] is a Hall interval
    // if the Hall interval was found on minsorted intervals, then we do the same thing in reverse order
    // i.e., in the interval [m, n-1]


    void sortit();
    int filterlower();
    int filterupper();
    void propagateValue();
    //@}

  public:
    /**@name Constructors*/
    //@{
    ConstraintAllDiff() : GlobalConstraint() { priority = 0; }
    ConstraintAllDiff(Vector< Variable >& scp);
    ConstraintAllDiff(std::vector< Variable >& scp);
    ConstraintAllDiff(Variable* scp, const int n);
    virtual void mark_domain();
    virtual Constraint clone() { return Constraint(new ConstraintAllDiff(scope)// , type
						   ); }
    virtual void initialise();
    virtual ~ConstraintAllDiff();
    //void init();
    virtual int idempotent() { return 1;}
    virtual int postponed() { return 1;}
    virtual int pushed() { return 1;}
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const ;
    virtual PropagationOutcome propagate();
    //virtual PropagationOutcome propagate(const int changed_idx, const Event evt);
    //virtual RewritingOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    //virtual std::string getString() const ;
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "alldiff"; }
    //@}

    void print_structs();

#ifdef _ALLDIFF_WC
    void weight_conflict(double unit, Vector<double>& weights);//  {
#endif
  };



  /***********************************************
   * Global Cardinality Constraint (bounds consistency).
   ***********************************************/
  typedef struct {
    int firstValue;
    int lastValue;
    int* sum;
    int* ds;
  } partialSum;
   

  /*! \class ConstraintOccurrences
    \brief  Global Cardinality Constraint.

    User defined propagator for enforcing bounds consistency
    on the restricted gcc constraint when bounds on
    occurrences are [a_i,b_i].
    A value "v" must be assigned to at least
    minOccurrences[v - firstDomainValue] variables and at most
    maxOccurrences[v - firstDomainValue] variables
    The code is from Lopez-Ortiz, Quimper, Tromp and van Beek
  */   
  class ConstraintOccurrences : public GlobalConstraint {
    
  private:
    /**@name Parameters*/
    //@{  
    int lastLevel;
    int *t;			// tree links
    int *d;			// diffs between critical capacities
    int *h;			// hall interval links
    int *stableInterval;	// stable sets
    int *potentialStableSets;	// links elements that potentialy belong to same stable set
    int *newMin;
    AD_Interval *iv;
    AD_Interval **minsorted;
    AD_Interval **maxsorted;
    int *bounds;  // bounds[1..nb] hold set of min & max of the n intervals
    // while bounds[0] and bounds[nb+1] allow sentinels
    int nb;
  
    partialSum* l; 
    partialSum* u;
    partialSum* initializePartialSum(const int firstValue, 
				     int count, const int* elements);
    void destroyPartialSum(partialSum *p);
    int  sum(partialSum *p, int from, int to);
    int  searchValue(partialSum *p, int value);
    int  minValue(partialSum *p);
    int  maxValue(partialSum *p);
    int  skipNonNullElementsRight(partialSum *p, int value);
    int  skipNonNullElementsLeft(partialSum *p, int value);
  
    void sortit();
    int  filterLowerMax();
    int  filterUpperMax();
    int  filterLowerMin(int *tl, int *c,
			int* stableAndUnstableSets,
			int* stableInterval,
			int* potentialStableSets,
			int* newMin);
    int  filterUpperMin(int *tl, int *c,
			int* stableAndUnstableSets,
			int* stableInterval,
			int* newMax);
    int getlb(const int val) const;
    int getub(const int val) const;
    //@}  

  public:
    /**@name Constructors*/
    //@{
    ConstraintOccurrences() : GlobalConstraint() { priority = 0; }
    ConstraintOccurrences(Vector< Variable >& scp,
			  const int firstDomainValue,
			  const int lastDomainValue,
			  const int* minOccurrences,
			  const int* maxOccurrences);
    ConstraintOccurrences(std::vector< Variable >& scp,
			  const int firstDomainValue,
			  const int lastDomainValue,
			  const int* minOccurrences,
			  const int* maxOccurrences);
    // ConstraintOccurrences(Variable* scp, const int n,
    // 			  const int firstDomainValue,
    // 			  const int lastDomainValue,
    // 			  const int* minOccurrences,
    // 			  const int* maxOccurrences);
    virtual void mark_domain();
    virtual Constraint clone() { return Constraint(new ConstraintOccurrences(scope, l->firstValue+3, l->lastValue-2, l->sum+3, u->sum+3)// , type
						   ); }
    virtual void initialise();
    virtual ~ConstraintOccurrences();
    //void init();
    virtual int idempotent() { return 1;}
    virtual int postponed() { return 1;}
    virtual int pushed() { return 1;}
    //@}

    /**@name Solving*/
    //@{
    virtual int check( const int* sol ) const ;
    virtual PropagationOutcome propagate();
    //virtual PropagationOutcome propagate(const int changed_idx, const Event evt);
    //virtual RewritingOutcome rewrite();
    //@}

    /**@name Miscellaneous*/
    //@{  
    //virtual std::string getString() const ;
    virtual std::ostream& display(std::ostream&) const ;
    virtual std::string name() const { return "occurrences"; }
    //@}
  };




  std::ostream& operator<< (std::ostream& os, const Constraint& x);
  std::ostream& operator<< (std::ostream& os, const Constraint* x);

  std::ostream& operator<< (std::ostream& os, const ConstraintImplementation& x);
  std::ostream& operator<< (std::ostream& os, const ConstraintImplementation* x);


  std::ostream& operator<< (std::ostream& os,  ConstraintTriggerArray& x);
  std::ostream& operator<< (std::ostream& os,  ConstraintTriggerArray* x);


}

#endif //__CONSTRAINT_HPP




