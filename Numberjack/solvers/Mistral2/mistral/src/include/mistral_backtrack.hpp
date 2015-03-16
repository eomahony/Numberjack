
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


/*! \file mistral_backtrack.hpp
  \brief Header for the reversible structures
*/


#ifndef _MISTRAL_BACKTRACK_HPP
#define _MISTRAL_BACKTRACK_HPP


#include <mistral_global.hpp>


namespace Mistral {


  /*!
    A reversible structure keeps a pointer to an
    "Environment" (i.e., Solver) that manages the
    backtracking process.
    Second, it implements a virtual function restore
    that "undo" the last change. 
    Most reversible structures work this way:
    A Vector<int> trail_ is used to encode the
    delta-information used to undo. The last
    integer stored on the trail_ Vector is the
    value of solver->level when the last change
    occured. It is used, when changing the object,
    to decide if one should "save" the current state,
    or simply replace the current value.

    It is still undecided if this process can be 
    made generic enough so that we can make 
    restore() a static method.
  */


  /*! \class Trigger
    \brief List of constraints
  */
  /********************************************
   * Trigger Objects
   ********************************************/
  /**
     A trigger is list of constraints that should be triggerred on a given event
     (there will be three lists per variable, one for domain events, one for range events, and one for value events)
  
     The difference with vectors is that when a constraint is added to such a list, the index in the vector is returned, 
     then stored into the constraint object :index: (and linked to the variable :index[var]:) for further reference (i.e., relax itself)
  
     Similarly, when removing a constraint, since it is swapped with the last element in the list, the index of this last element 
     must also be updated (set to the value of the relaxed constraint). A constraint from a trigger of variable will give index -1 to this variable
  */
  class Constraint;
  class Trigger : public Vector< Constraint > {

  public:

    virtual ~Trigger() {
#ifdef _DEBUG_MEMORY
  std::cout << "c delete trigger" << std::endl;
#endif
    }
    
    int check_and_post(Constraint ct) ;

    int post(Constraint ct) ;

    void relax(const int idx) ;

  };




  /*! \class Constraint
    \brief A wrapper class for constraints
  */
  /********************************************
   * Constraint Objects
   ********************************************/
  /**
     This class is mainly used to
     1/ resolve subtypes at runtime without virtual overloading 
     (there are only three main types of constraints: binary, ternary and global)

     2/ manipulate constraints as objects instead of pointers

     The "real" constraints are objects of the class ConstraintImplementation, pointed by "propagator"
     The "data" field stores some data about the constraint (mainly its type), but also some
     context dependent information. For instance, as element of a trigger, the field "data" 
     will store the index of that variable in its scope 
  */
  template<class T> class ReversibleNum;
  class ReversibleSet;
  class Reversible;
  class Variable;
  class Decision;
  class Solver;
  class ConstraintImplementation;
  class Constraint  {
    
  public:
    
    ConstraintImplementation* propagator;
    unsigned int data;
    
    Constraint(); // {propagator = NULL; data=0;}
    Constraint(ConstraintImplementation* p) ;
    Constraint(ConstraintImplementation* p, const int t) { propagator = p; data = t; }
    void initialise(Solver*);
    virtual ~Constraint() {}
    
    // return the context-dependent information [TODO: list all types]
    inline int  index()      const { return data&CTYPE; }
    // return the context-dependent information [TODO: list all types]
    inline int  info()      const { return data&ITYPE; }

    // [context-independent:] whether the constraint should be pushed on the constraint stack on  var event
    inline bool pushed()     const { return (data&PUSHED); }    

    // [context-independent:] whether the constraint should not be propagated on a variable event 
    inline bool postponed()  const { return (data&POSTPONED); }
    // [context-independent:] whether the constraint is of the "global" subtype (warning: it might be binary or ternary)
    inline bool global()     const { return !(data&0xc0000000); }
    // [context-independent:] whether the constraint is of the "binary" subtype 
    inline bool binary()     const { return (data&BINARY); }
    // [context-independent:] whether the constraint is of the "binary" subtype 
    inline bool ternary()    const { return (data&TERNARY); }
    // [context-independent:] whether the constraint is idempotent (it should not be triggerred on its own changes)
    inline bool idempotent() const { return (data&IDEMPOTENT); }
    // [context-independent:] whether the constraint is able to explain its pruning
    bool explained() const;//  { 
    //   return propagator->explained();
    //   //return (data&EXPLAINED); }
    // }
    //void set_idempotent(const bool flag);

    // stop pointing to a constraint
    inline void clear() { propagator = NULL; data = 0; }
    // whether the wrapper actually points to a constraint
    inline bool empty() const { return !propagator; }

    // returns the unique id of the pointed constraint (equal to its rank in the "constraints" array of the Solver class)
    int id() const ;
    // returns the symbol of the pointed constraint 
    std::string symbol() const ;
    // return the priority of the constraint during AC closure (the constraint list will have as many sub-lists as there are priorities)
    int priority() const ;
    // set the context-dependent part of the data field to the value idx
    void set_index(const int idx);
    // set the value of the "index[index()]" array of the pointed constraint to the value idx
    void set_rank(const int idx);
    // set the unique id of the pointed constraint to the value idx
    void set_id(const int idx);
    // post the constraint to the solver (initial post)
    void post(Solver*);
    // re-post the constraint
    void awaken();
    // activate the constraint (it will be propagated during the next call to solver.propagate())
    void trigger();
    // check that the assignment sol satisfies the constraint
    int check(int* sol);
    // change the variables of the pointed constraint to their "final" objects
    void consolidate();
    // change the variable corresponding to "index()" of the pointed constraint to its "final" object
    void consolidate_var();
    // set the constraint to be triggerred with t
    void re_link_to(Trigger* t);
    // whether the pointed constraint is succeptible to be rewritten
    bool rewritable(); 
    // a temporal solution for rewriting in flatzinc
    bool simple_rewritable();
    // whether the constraint can be rewritten to accommodate the not(x[i]) instead of x[i]
    bool absorb_negation(const int i);
    // return the constraint C such that C(x1,..,xi-1,not(xi),i+1,..,xk) <=> this(x1,..,xk)
    // and where the occurrence of not(xi) is taken by x
    Constraint get_negation(const int i, Variable x);
    // get weights for vars/literals
    // the weight represents (sometimes roughly) the probability that 
    // a random assignment of the constraint containing this literal will not be allowed
    // The variable weight is the sum of all the weights of its values/literal
    void initialise_activity(double *lvact, double *vact, double norm);

    // to be called before propagation, returns the constraint itself if it is idempotent and the NULL pointer otherwise
    ConstraintImplementation *freeze();
    // to be called after propagation
    ConstraintImplementation *defrost();


    // returns the arity of the pointed constraint
    int arity() const ;
    // returns the number of active (i.e. non-ground) variables in the constraint's scope
    int num_active() const ;
    // returns the index of the ith active variable of the constraint
    int get_active(const int i) const ;
    // returns an array containing the constrained variables 
    Variable* get_scope() ;
    // set the ith variable of the scope to x
    void set_scope(const int i, Variable x) ;

    // relax the constraint from all its active triggers
    void relax();

    // void relax();
     void notify_assignment();
    bool is_active();

// int num_active();
    // Variable* get_scope();
    // void set_scope(int i, Variable);
    // int arity();

    void check_active();
    int rank();


    // weight the variables according to how much they are responsible for the last conflict
    void weight_conflict(double unit, Vector<double>& weights);

    // call the rewriting procedure of the constraint
    RewritingOutcome rewrite();

    bool find_support(const int var, const int val);
    bool find_bound_support(const int var, const int val);
    PropagationOutcome propagate();
    PropagationOutcome propagate(const Event evt);

    PropagationOutcome checker_propagate();
    PropagationOutcome checker_propagate(const Event evt);

    PropagationOutcome bound_checker_propagate();
    PropagationOutcome bound_checker_propagate(const Event evt);

    void restore();

    int get_backtrack_level();
    Decision get_decision();
    
    bool operator==(Constraint c) const {return propagator == c.propagator;}
    bool operator!=(Constraint c) const {return propagator != c.propagator;}

    std::ostream& display(std::ostream& os) const;

  };


  /*! \class VarEvent
    \brief Stores the variable id, the type of change, and the constraint that triggered this event
  */
  /********************************************
   * Environment Objects
   ********************************************/
  class VarEvent : public Triplet < int, Event, ConstraintImplementation*> {
    
  public:

    VarEvent(int t1, Event t2, ConstraintImplementation* t3) : Triplet < int, Event, ConstraintImplementation*>(t1, t2, t3) {}

    operator int() {
      return first;
    }

    void update(const Event evt, const ConstraintImplementation* con) {
      second |= evt;
      if(third != con) third = NULL;
    }

    std::ostream& display(std::ostream& os) const {
      os << event2strc(second) << "(" << first << ")";
      return os;
    }
  };

  // /*! \class VarEvent
  //   \brief Stores the variable id, the type of change, and the constraint that triggered this event
  // */
  // /********************************************
  //  * Environment Objects
  //  ********************************************/
  // class VarEvent : public Triplet < int, Event, Constraint > {
    
  // public:

  //   VarEvent(int t1, Event t2, Constraint t3) : Triplet < int, Event, Constraint >(t1, t2, t3) {}

  //   operator int() {
  //     return first;
  //   }

  //   void update(const Event evt, const Constraint con) {
  //     second |= evt;
  //     if(third != con) third.clear(); // = Constraint();
  //   }

  //   std::ostream& display(std::ostream& os) const {
  //     os << event2strc(second) << "(" << first << ")";
  //     return os;
  //   }
  // };


  typedef TwoWayStack< VarEvent > VariableQueue;
  //typedef TwoWayStack< Triplet < int, Event, ConstraintImplementation*> > VariableQueue;


  /*! \class Environment
    \brief The minimal structures used to control the backtracking process
  */
  /********************************************
   * Environement Objects
   ********************************************/
  class Environment {
  public:
    
    /*!@name Parameters*/
    //@{
    int level;

    Vector< int >                 saved_vars;
    Vector< Constraint >          saved_cons;
    Vector< int* >                saved_bools;
    Vector< ReversibleSet* >      saved_lists;
    Vector< ReversibleNum<int>* > saved_ints;

    /// The delimitation between different levels is kept by this vector of integers
    Vector< int > trail_;

    VariableQueue active_variables;

    ConstraintImplementation *taboo_constraint;
    //@}

    /*!@name Constructors*/
    //@{
    Environment() { 
      level = 0;
      taboo_constraint = NULL;
    }
    virtual ~Environment() {}
    //@}


    /*!@name Backtrack method*/
    //@{
    inline void save() {

      trail_.add(saved_vars.size);
      trail_.add(saved_bools.size);
      trail_.add(saved_lists.size);
      trail_.add(saved_ints.size);
      trail_.add(saved_cons.size);

      ++level;

    }


    void trigger_event(const int var, const Event evt) {
      if(active_variables.contain(var)) {
	active_variables[var].update( evt, taboo_constraint );
      } else {
	//Triplet< int, int, ConstraintImplementation* > 
	VarEvent t(var, evt, taboo_constraint);
	active_variables.push_back(t);
      }
    }

    void _restore_();

    inline void save(ReversibleNum<int> *r) {saved_ints.add(r);}
    inline void save(ReversibleSet *r) {saved_lists.add(r);}
    inline void save(int *r) {saved_bools.add(r);}
    
    inline void save(int r) {saved_vars.add(r);}
    inline void save(Constraint r) {

#ifdef _DEBUG_BACKTRACK
      int id = r.id();
      if(_DEBUG_BACKTRACK) {
	for(int i=0; i<level; ++i) std::cout << " ";
	std::cout << "c save [" << r.id() << "]" ;
	r.display(std::cout);
	
	int info = r.info();
	if(info & ACTIVITY) {
	  std::cout << " (change on activity)";
	}
	if(info & RELAXED) {
	  std::cout << " (was relaxed)";
	}
	if(info & POSTED) {
	  std::cout << " (was posted)";
	}
	
	std::cout << std::endl;
      }
#endif

      saved_cons.add(r);


    }
    //@}

  };

 


  /********************************************
   * Reversible Objects
   ********************************************/
  /*! \class Reversible
    \brief Backtrackable data structures.

    All structures that need to be restored to previous
    state during search implement the methods 
    save() and restore(). 
  */
  //class Environment;
  class Reversible {
  public:
    /*!@name Parameters*/
    //@{
    Environment *env;
    //@}

    /*!@name Constructors*/
    //@{
    Reversible() { env=NULL; }
    Reversible(Environment *s) : env(s) {;}
    void initialise(Environment *s) {env=s;}
    virtual ~Reversible() {}
    //@}

  };



  /********************************************
   * Reversible Boolean Domain
   ********************************************/
  /*! \class ReversibleBool
    \brief Backtrackable Primitive
  */
  class ReversibleBool : public Reversible
  {
  
  public:
    /*!@name Parameters*/
    //@{  
    int value;
    //@}

    /*!@name Constructors*/
    //@{ 
    ReversibleBool() : Reversible() {
      value = 3;
    }
    ReversibleBool(Environment *s) 
      : Reversible(s)
    {
      value = 3;
    }

    void initialise(Environment *s) {
      Reversible::initialise(s);
    }
    
    virtual ~ReversibleBool() {}
    //@}


    /*!@name Backtrack method*/
    //@{
    inline void save() { 
      if(value == 3) env->save(&value); 
    }
    inline void restore() { 
      value = 3;
    }
    //@}

    /*!@name Manipulation*/
    //@{  
    inline void remove(const int v) { 
      save();
      value = 2-v;
    }

    inline void operator=(const int v) { 
      save();
      value = 1+v;
    }

    inline bool contain(const int v) {
      return value&(v+1);
    }

    inline bool equal(const int v) {
      return (value&(v+1)) == value;
    }
    
    inline int size() { return (1+(value==3)); }
    //@}
    
    /*!@name Printing*/
    //@{
    std::ostream& display(std::ostream& os) const { os << (value == 3 ? "{0,1}" : (value == 1 ? "0" : "1")) ; return os; }
    //@}
  };

  /********************************************
   * Reversible Primitive Type
   ********************************************/
  /*! \class ReversibleNum
    \brief Backtrackable Primitive
  */
  template < class PRIMITIVE_TYPE >
  class ReversibleNum : public Reversible
  {
  
  public:
    /*!@name Parameters*/
    //@{  
    /// value trail
    Vector< PRIMITIVE_TYPE > trail_;
    /// current value
    PRIMITIVE_TYPE value;
    //@}

    /*!@name Constructors*/
    //@{ 
    ReversibleNum() : Reversible() {
    }
    ReversibleNum(const PRIMITIVE_TYPE v) 
    {
      initialise(v);
    }
    ReversibleNum(Environment *s) 
      : Reversible(s)
    {
      Reversible::initialise(s);
    }
    ReversibleNum(Environment *s, const PRIMITIVE_TYPE v) 
      : Reversible(s)
    {
      Reversible::initialise(s);
      initialise(v);
    }

    void initialise(Environment *s, const PRIMITIVE_TYPE v) 
    {
      Reversible::initialise(s);
      initialise(v);
    }

    void initialise(Environment *s) 
    {
      Reversible::initialise(s);
    }

    void initialise(const PRIMITIVE_TYPE v) 
    {
      value = v;
      trail_.initialise(0, 16);
      trail_.add(value);
      trail_.add(-1);
    }
    virtual ~ReversibleNum() {}
    //@}

    /*!@name Accessors*/
    //@{  
    inline operator const PRIMITIVE_TYPE() const { return value; }
    //@}  

    /*!@name Backtrack method*/
    //@{
    inline void save() { 
      if((int)trail_.back() != env->level) { 
	env->save(this); 
	trail_.add((int)value); 
	trail_.add(env->level); 
      } 
    }
    inline void restore() { 
      trail_.pop();//current_level); 
      value = (PRIMITIVE_TYPE)(trail_.pop()); 
    }
    //@}

    /*!@name Manipulation*/
    //@{  
    inline bool operator!() { return !value; }
    inline PRIMITIVE_TYPE operator-() { return -value; }
    inline void operator=  ( const PRIMITIVE_TYPE x ) { save(); value  = x; }
    inline void operator+= ( const PRIMITIVE_TYPE x ) { save(); value += x; }
    inline void operator-= ( const PRIMITIVE_TYPE x ) { save(); value -= x; }
    inline void operator*= ( const PRIMITIVE_TYPE x ) { save(); value *= x; }
    inline void operator/= ( const PRIMITIVE_TYPE x ) { save(); value /= x; }
    inline void operator|= ( const PRIMITIVE_TYPE x ) { save(); value |= x; }
    inline void operator&= ( const PRIMITIVE_TYPE x ) { save(); value &= x; }
    inline void operator^= ( const PRIMITIVE_TYPE x ) { save(); value ^= x; }
    inline ReversibleNum< PRIMITIVE_TYPE >& operator++ () { save(); ++value; return *this; }
    inline ReversibleNum< PRIMITIVE_TYPE >& operator-- () { save(); --value; return *this; } 
    //@}

    std::ostream& display(std::ostream& os) const {
      os << value;
      return os;
    }
  };


  /********************************************
   * Reversible Set
   ********************************************/
  /*! \class ReversibleSet
    \brief Backtrackable IntStack
  */
  class ReversibleSet : public Reversible, public IntStack {
    
  public:

    typedef int* iterator;

    /*!@name Parameters*/
    //@{  
    /// value trail
    Vector< int > trail_;
    //@}

    /*!@name Constructors*/
    //@{ 
    ReversibleSet() : Reversible(), IntStack() {}
    ReversibleSet(Environment *s, const int lb=0, const int ub=0, const int sz=-1, const bool full=true)
      : Reversible(s)
    {
      if(sz < 0)
	initialise(lb, ub, ub-lb+1, full);
      else
	initialise(lb, ub, sz, full);
    }

    virtual ~ReversibleSet()
    {
    }

    using IntStack::initialise;
    virtual void initialise(Environment *s)
    {
      Reversible::initialise(s);
    }

    // virtual void initialise(Environment *s, ReversibleSet *e, const int)
    // {
    //   Reversible::initialise(s);
    //   initialise(e);
    //   // trail_.add(size);
    //   // trail_.add(-1);
    // }

    virtual void initialise(Environment *s, const int lb, const int ub, const int sz, const bool full)
    {
      Reversible::initialise(s);
      initialise(lb, ub, sz, full);
      // trail_.add(size);
      // trail_.add(-1);
    }

    virtual void initialise(const int lb, const int ub, const int sz, const bool full)
    {
      int l = lb;
      int u = ub;
      if(l>u) {
	u = l-1;
	l = 0;
      }
      IntStack::initialise(l, u, sz, full);
      trail_.initialise(0, 2*sz);
      trail_.add(size);
      trail_.add(-1);
    }

    virtual void initialise(Environment *s, const int lb, const int ub, const Vector<int>& vals)
    {
      Reversible::initialise(s);
      initialise(lb, ub, vals);
      // trail_.add(size);
      // trail_.add(-1);
    }

    virtual void initialise(ReversibleSet& shared, const int sz)
    {
      Reversible::initialise(shared.env);

      index_capacity = shared.index_capacity;
  
      list_capacity = sz;
      list_ = new int[list_capacity];
  
      start_ = shared.start_;
      index_ = shared.index_;
  
      size = 0;

      trail_.initialise(0, 2*sz);
      trail_.add(size);
      trail_.add(-1);
      //IntStack::initialise((IntStack)shared, sz);
      // trail_.add(size);
      // trail_.add(-1);
    }

    virtual void initialise(const int lb, const int ub, const Vector<int>& vals)
    {
      IntStack::initialise(lb, ub, vals.size, false);
      for(unsigned int i=0; i<vals.size; ++i)
	init_add(vals[i]);
      trail_.initialise(0, 2*size);
      trail_.add(size);
      trail_.add(-1);
    }
    //@}

    /*!@name Backtrack method*/
    //@{    
    int get_reduction() const;
    void restore();
    void save();
    //@}

    /*!@name Manipulation*/
    //@{  
    // it's either 'add's...
    void reversible_add(const int elt);
    void reversible_add(const Vector<int>& elt);

    // ...or remove, but not both!!
    void reversible_remove(const int elt);
    void reversible_remove(const Vector<int>& elt);

    void reversible_set_to(const int elt);
    int reversible_pop();

    int reversible_pop_head();
    //@}    

    // /*!@name Backtrack method*/
    // //@{    
    //  void restore() { 
    //   trail_.pop(); size = trail_.pop(); 
    // } 
    //  void save() { 

    //   //std::cout << trail_.size << " " << env << std::endl;

    //   if(trail_.back() != env->level) {
    // 	trail_.add(size);
    // 	trail_.add(env->level);
    // 	env->save(this);
    //   }
    // }
    // //@}

    // /*!@name Manipulation*/
    // //@{  
    // // it's either 'add's...
    //  void reversible_add(const int elt) {
    //   save();
    //   add(elt);
    // }

    // // ...or remove, but not both!!
    //  void reversible_remove(const int elt) {
    //   save();
    //   remove(elt);
    // }


    //  void reversible_set_to(const int elt) {
    //   save();
    //   set_to(elt);
    // }

    //  int reversible_pop()
    // {
    //   save();
    //   return pop();
    // }

    //  int reversible_pop_head()
    // {
    //   save();
    //   return pop_head();
    // }
    // //@}    

  };

  // /********************************************
  //  * Reversible Bitset
  //  ********************************************/
  // /*! \class ReversibleBitset
  //   \brief Backtrackable IntStack
  // */
  // template< class WORD_TYPE, class FLOAT_TYPE >
  // class ReversibleBitset : public Reversible, public Bitset< WORD_TYPE, FLOAT_TYPE > {
    
  // public:

  //   /*!@name Parameters*/
  //   //@{  
  //   Vector< int > trail_;

  //   /// trail for the bitset representation
  //   WORD_TYPE **delta_;
  //   int **level_;
  //   WORD_TYPE **delta_abs;
  //   int **level_abs;
  //   //@}

  //   /*!@name Constructors*/
  //   //@{ 
  //   ReversibleBitset() : Reversible(), IntStack() {}
  //   ReversibleBitset(Environment *s, const int lb=0, const int ub=0, const bool full=true)
  //     : Reversible(s)
  //   {
  //     initialise(lb, ub, full);
  //   }

  //   virtual ~ReversibleBitset()
  //   {
  //   }

  //   virtual void initialise(Environment *s)
  //   {
  //     Reversible::initialise(s);
  //   }

  //   // virtual void initialise(Environment *s, ReversibleBitset *e, const int)
  //   // {
  //   //   Reversible::initialise(s);
  //   //   initialise(e);
  //   //   // trail_.add(size);
  //   //   // trail_.add(-1);
  //   // }

  //   virtual void initialise(Environment *s, const int lb, const int ub, const bool full)
  //   {
  //     Reversible::initialise(s);
  //     initialise(lb, ub, full);
  //     // trail_.add(size);
  //     // trail_.add(-1);
  //   }

  //   virtual void initialise(const int lb, const int ub, const bool full)
  //   {
  //     int l = lb;
  //     int u = ub;
  //     if(l>u) {
  // 	u = l-1;
  // 	l = 0;
  //     }
  //     IntStack::initialise(l, u, u-l+1, full);
  //     trail_.initialise(0, 2*(u-l+1));
  //     trail_.add(size);
  //     trail_.add(-1);
  //   }

  //   virtual void initialise(Environment *s, const int lb, const int ub, const Vector<int>& vals)
  //   {
  //     Reversible::initialise(s);
  //     initialise(lb, ub, vals);
  //     // trail_.add(size);
  //     // trail_.add(-1);
  //   }

  //   virtual void initialise(ReversibleBitset& shared, const int sz)
  //   {
  //     index_capacity = shared.index_capacity;
  
  //     list_capacity = sz;
  //     list_ = new int[list_capacity];
  
  //     start_ = shared.start_;
  //     index_ = shared.index_;
  
  //     size = 0;

  //     trail_.initialise(0, 2*sz);
  //     trail_.add(size);
  //     trail_.add(-1);
  //     //IntStack::initialise((IntStack)shared, sz);
  //     // trail_.add(size);
  //     // trail_.add(-1);
  //   }

  //   virtual void initialise(const int lb, const int ub, const Vector<int>& vals)
  //   {
  //     IntStack::initialise(lb, ub, vals.size, false);
  //     for(unsigned int i=0; i<vals.size; ++i)
  // 	init_add(vals[i]);
  //     trail_.initialise(0, 2*size);
  //     trail_.add(size);
  //     trail_.add(-1);
  //   }
  //   //@}

  //   /*!@name Backtrack method*/
  //   //@{    
  //   inline void restore() { 
  //     trail_.pop(); size = trail_.pop(); 
  //   } 
  //   inline void save() { 
  //     if(trail_.back() != env->level) {
  // 	trail_.add(size);
  // 	trail_.add(env->level);
  // 	env->save(this);
  //     }
  //   }
  //   //@}

  //   /*!@name Manipulation*/
  //   //@{  
  //   // it's either 'add's...
  //   inline void reversible_add(const int elt) {
  //     save();
  //     add(elt);
  //   }

  //   // ...or remove, but not both!!
  //   inline void reversible_remove(const int elt) {
  //     save();
  //     remove(elt);
  //   }


  //   inline void reversible_set_to(const int elt) {
  //     save();
  //     set_to(elt);
  //   }

  //   inline int reversible_pop()
  //   {
  //     save();
  //     return pop();
  //   }

  //   inline int reversible_pop_head()
  //   {
  //     save();
  //     return pop_head();
  //   }
  //   //@}    

  // };


  template < class PRIMITIVE_TYPE > 
  std::ostream& operator<< (std::ostream& os, const ReversibleNum< PRIMITIVE_TYPE >& x) {
    return x.display(os);
  }

  template < class PRIMITIVE_TYPE > 
  std::ostream& operator<< (std::ostream& os, const ReversibleNum< PRIMITIVE_TYPE >* x) {
    return x->display(os);
  }

  std::ostream& operator<< (std::ostream& os, const ReversibleBool& x);
  std::ostream& operator<< (std::ostream& os, const ReversibleBool* x);


  std::ostream& operator<< (std::ostream& os, const VarEvent& x);
  std::ostream& operator<< (std::ostream& os, const VarEvent* x);

}

#endif //__BACKTRACK_HPP
