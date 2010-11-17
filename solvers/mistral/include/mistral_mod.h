
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



/** \file mod.h
    \brief Header for the variables and constraints wrappers and build objects.
*/

/**********************************************
 * Wrappers & Build object
 **********************************************/

#ifndef _Wrapper_H
#define _Wrapper_H

#include <mistral_dvo.h>
#include <mistral_sat.h>
//#include <Mistral.hpp>

namespace Mistral {
  
  //class Mistral_Expression;
  class CSP;
  class Variable;
  class VarArray;
  class SatSolver;
  class ConstraintClauseBase;

  /**********************************************
   * MistralGacStack
   **********************************************/

  /// Stack of Variable for computing the GAC closure.
  template <class T>
  class MistralGacStack {
  public:
    T *MistralGacStack_;
    int        *trigger;
    int            size;
  
    /**@name Constructors*/
    //@{
    ///
    MistralGacStack()
    {
      size=0;
      MistralGacStack_=NULL;
      trigger=NULL;
    }
    ///
    MistralGacStack(int n)
    {
      initStack(n);
    }
    ///
    ~MistralGacStack()
    {
      delete [] trigger;
      delete [] MistralGacStack_;
    }
    ///
    void initStack(int n)
    { 
      trigger = new int[n];
      std::fill( trigger, trigger+n, 0 );
      MistralGacStack_ = new T[n];
      size = 0;
    }
    //@}
  
    void push(T v, int event)
    {
      if( !trigger[v->id] )
	MistralGacStack_[size++] = v;    
      trigger[v->id] |= event;
    }

    void push(T v, VariableInt* x)
    {
      if( !trigger[v->id] )
	MistralGacStack_[size++] = v;    
      trigger[v->id] = 1;
    }


    T pop( int& event )
    {
      T last = MistralGacStack_[--size];
      event = trigger[last->id];
      trigger[last->id] = 0;
      return last;
    }

    T back( int& event )
    {
      T last = MistralGacStack_[size-1];
      event = trigger[last->id];
      return last;
    }

    bool empty() const
    {
      return (!size);
    }

    void clear()
    {
      while( size )
	trigger[MistralGacStack_[--size]->id] = 0; 
    }
  };


  ////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////
  //////                                                                            //////
  //////  There are 3 layers:                                                       //////
  //////  1/ Wrappers are simply syntactic sugar that allow to                      //////
  //////  a/ hide 'new' operations in models                                        //////
  //////  b/ allow to overload operator (+, [], ...)                                //////
  //////  2/ Build objects are 'temporary' variables and/or constraints.             //////
  //////  They contain only the necessary information (domain, scope, ...)          //////
  //////  and the method 'build' actually create the real object                    //////
  //////  used during search. They allow to do some re-modelling                    //////
  //////  before the solving phase.                                                 //////
  //////  For instance, if a variable is an interval and is                         //////
  //////  constrained only by BC constraints, the RangeVar                          //////
  //////  implementation will be used instead of BitVar.                            //////
  //////   + BuildObject stand for variables. They keep a representation            //////
  //////     of their domain, as well as a number of informations in the            //////
  //////     struct state.                                                          //////
  //////   + BuildObjectConstraints stand for constraints' *relation*,              //////
  //////     that is, at most one of each is created, no matter how many            //////
  //////     such constraints are added to the model.                               //////
  //////   + BuildObjectPredicate stand for constraint scope. They ihnerit          //////
  //////     from BuildObject and therefore can also represent variables            //////
  //////     themselves.                                                            //////
  //////     - When a constraint, or variable wrapper is declared,                  //////
  //////       the corresponding BuildObject is created                             //////
  //////     - When a (tree of) objects are added to the model                      //////
  //////       the top-level one is marked as constraint only, whilst               //////
  //////       all descendants are marked as variables.                             //////


  //////     - During the prebuilding phase, two things happen                      //////
  //////       * The method propagateInfo() is called on all                        //////
  //////         top-level objects, which recursively call their                    //////
  //////         children.                                                          //////
  //////         . Constraints that become incorrect when reasoning on bounds       //////
  //////         . Constraints that can profit from holes                           //////
  //////       * The method prebuild() is called on all objects from                //////
  //////         top to bottom to perform rewriting when needed.                    //////

  /*
    + transformations: (limited propagation from top to bottom + rewriting thereof)
    > general:
    - if all arguments are constants, delete the declaration (check consistency?)
    unsetRel(): the relation part should not be build, because it is 
    a totology
    deReference() all arguments. Each argument loses a reference 
    > Not
    - if set to true or at toplevel (!isVar()):
    unsetRel(), set argument to false (0).
    - if set to false:
    unsetRel(), set argument to false (0).
    > And 
    - if set to true or at toplevel (!isVar()):
    unsetRel(), set arguments to true (1).
    > Or
    - if set to false:
    unsetRel(), set arguments to false (0).
    > Add 
    - if an argument is  a constant:
    turn into an offset predicate
    > Sub
    - if an argument is  a constant:
    turn into an offset predicate
    > EqualConstant
    - if set to true or at toplevel (!isVar()):
    unsetRel(), set argument to val.
    - if set to false:
    unsetRel(), remove val in arg's domain.
    > NotEqualConstant
    - if set to true or at toplevel (!isVar()):
    unsetRel(), remove val in arg's domain.
    - if set to false:
    unsetRel(), set argument to val.
    > Equal
    - if set to true or at toplevel (!isVar()):
    unsetRel(), merge arguments.
    - if an argument is  a constant:
    turn into an EqualConstant predicate
    > NotEqual
    - if set to true or at toplevel (!isVar()), and if one of the agrs is a constant
    unsetRel(), and remove the constant from the other agr's domain.
    - if an argument is  a constant:
    turn into an NotEqualConstant predicate
    > SumEq
    - absorb constants.
    > Sum
    - if the last arg is a constant and is not referenced otherwise:
    turns into a SumEq constraint
    - absorb constants then according to the arity,
    change to offset/factor/add/sub/sumeq
    > if N is constant:
    - turns into an equal predicate

  */

  //////                                                                            //////
  //////  3/ Actual variables and constraints (csp.h / var.h / con.h)               //////
  //////                                                                            //////
  ////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////
  /* 

     propagate info

     * needed because when rewriting, the state of an object can change,
     the change must be propagated downstream.


     mechanism to decide if the representation should be predicate 
     or constraint.
     * if not at top level -> always variable
     * some object are variables even at the top level 
     (actual variables, element, sum, etc)
     * top level known when adding, but it can change 
     (for instance with a AND constraint)
  */




  class BuildObjectConstraint;
  /**********************************************
   * ConstraintStore
   **********************************************/
  /*! \class ConstraintStore
    \brief A pool of constraint build objects

    We need only one instance of each constraint
    buld object.
  */
  class ConstraintStore
  {
  public: 
    // LOGIC
    static const int AND         =  1 ;
    static const int OR          =  2 ;
    static const int NOT         =  3 ;

    // ARITHMETIC
    static const int ABS         =  4 ;
    static const int DIV         =  5 ;
    static const int MOD         =  6 ;
    static const int MUL         =  7 ;
    static const int NEG         =  8 ;

    // RELATION
    static const int DISJUNCTIVE =  9 ;
    static const int EQUAL       = 10 ;
    static const int MEMBER      = 11 ;
    static const int PRECEDENCE  = 12 ;
    static const int OVERLAP     = 13 ;

    // GLOBAL
    static const int ALLDIFF     = 14 ;
    static const int ELEMENT     = 15 ;
    static const int GCC         = 16 ;
    static const int LEX         = 17 ;
    static const int DLEX        = 18 ;
    static const int MAX         = 19 ;
    static const int MIN         = 20 ;
    static const int CARD        = 21 ;
    static const int SUM         = 22 ;
    static const int IFTHENELSE  = 23 ;
    static const int TREE        = 24 ;
    static const int TDAG        = 25 ;
    static const int CLAUSE      = 26 ;

    // COUNT
    static const int NUMCONS     = 27 ;


    // private:
    Vector<BuildObjectConstraint *> store;
    BuildObjectConstraint *getConstraint(const int c);

  public:
    ConstraintStore() ;
    virtual ~ConstraintStore() ;
    BuildObjectConstraint *operator[]( const int c ) ;
    int addConstraint( BuildObjectConstraint *c ) ;
    void remConstraint( int idx ) ;
    void flush();
  };


  //static ConstraintStore ENV;


  /********************************************************
   *
   *  BuildObject for Variables, Constraints 
   *  and Predicates.
   *
   ********************************************************/

  class BuildObjectPredicate;
  /**********************************************
   * BuildObject
   **********************************************/
  class BuildObject
  {

  protected:
    //public:
  
    // Domain:
    int size_;
    int min_;
    int max_;
    BitSet *domain_;  

  public: 

    // pointer to the object if reference (NULL otherwise)
    BuildObject *veqptr_;

    // info
    int state;

    // number of parents
    int refs;


    // identifiant
    int id;

    // constraints on this variable
    Vector<BuildObjectPredicate*> parent;

    // pointer to the model
    CSP *model;

    // pointer to the actual variable (NULL as long as not build)
    VariableInt *varptr_;

    BuildObject *getBuildObject();
    const BuildObject *getConstBuildObject() const ;

    virtual const char* get_type() const {return (veqptr_ ? veqptr_->get_type() : "var");}
    virtual int get_arity() {return 0;}

    int getState() const;
    void setState( const int s ) ;
    
    int getVariableId() const ;
    int getVarId() const ;
    //   int getState() ;
    //   void setState( int nstate ) ;

    static const int RANGE  =    1;
    static const int BLIST =     2;
    static const int ILIST =     4;
    static const int TOP    =    8;
    static const int NEQ    =   16;
    static const int VIRTUAL=   32;

    void releaseEvent( const int );

    /*!
     * set to 1 iff it can be implemented as an interval
     */
    virtual bool isRange() const ;
    virtual void setRange() ;
    virtual void unsetRange() ;

    /*!
     * set to 1 iff it can be a bit list
     */
    virtual bool isBList() const ;
    virtual void setBList() ;
    virtual void unsetBList() ;

    /*!
     * set to 1 iff it can be an int list
     */
    virtual bool isIList() const ;
    virtual void setIList() ;
    virtual void unsetIList() ;

    /*!
     * set to 1 iff it is a constraint
     */
    virtual bool isRel() const ;
    virtual void unsetRel() ;

    /*!
     * set to 1 iff it is strictly tighter than a neq constraint 
     */
    virtual bool isNeq() const ;
    virtual void setNeq() ;
    virtual void unsetNeq() ;

    /*!
     * set to 1 iff it is strictly tighter than a neq constraint 
     */
    virtual bool isVirtual() const ;
    virtual void setVirtual() ;
    virtual void unsetVirtual() ;

    /*!
     * set to 1 iff it is a top level constraint 
     */
    virtual bool isTop() const ;
    //virtual void setTop() ;
    //virtual void unsetTop() ;


  public:

    void reference() ;
    void _reference() ;
    void deReference() ;
    void _deReference() ;
    int isReferenced() const;
    int _isReferenced() const;

    virtual void add( CSP *model, const int l );
    virtual void close();

    // prebuild generally does the following things:
    /*
     * - if all variables are constants, 
     * - absorb constants as much as possible
     * -
     */
    /// return false if the predicate is "ground", true otherwise.
    virtual int propagateDownward() { return true; }
    virtual int propagateUpward() { return true; }
    virtual void build( Solver *s );
    virtual void build( SatSolver *s, const int*idx );

    //inline bool notIn() ;
    virtual bool isVar() ;
    virtual bool isBoolean() ;
    virtual bool isFDomain() ;
    virtual bool isInterval() ;
    virtual bool isConstant() ;
    virtual bool isSearchable() ;
    virtual bool isConstraint() ;
    virtual VariableInt *_getVariable() ;
    virtual VariableInt *getVariable() ;
  
    void initVar(const int l, const int u);
    BuildObject() ;
    BuildObject( const int lo, const int up ) ;
    BuildObject( const int* array, const int length ) ;
    BuildObject( const int* array, const int length, const int lo, const int up ) ;
    BuildObject( const BitSet& d, const int length, const int lo, const int up ) ;
    virtual ~BuildObject() ;

    int getValue() const ;
    int getMin() const ;
    int getMax() const ;
    bool isIn(const int v) const ;

    virtual int getSize() const ;
    virtual int size() const ;
    virtual int domsize() const ;
    virtual int value() const ;
    virtual int first() const ;
    virtual int last() const ;
    virtual int min() const ;
    virtual int max() const ;
    virtual int minCapacity() const ;
    virtual int maxCapacity() const ;
    virtual int minPosAbs() const ;
    virtual int minNegAbs() const ;
    virtual bool remove(const int v);
    virtual bool setDomain(const int v);
    virtual bool setMin(const int l);
    virtual bool setMax(const int u);
    virtual bool setDomain(const BitSet& s);
    virtual bool setEqual(BuildObject *v); 
    virtual bool removeSet(const BitSet& s);
    virtual bool removeSet(const int* vals, const int n);
    virtual bool removeRange(const int l, const int u);

    virtual int getNext(const int v) const ;
    virtual bool good(const int v) const ;
    virtual bool setNext(int& v) const ;

    virtual bool isInterval() const ;
    virtual bool isGround() const ;

    virtual int _size() const ;
    virtual int _domsize() const ;
    virtual int _value() const ;
    virtual int _first() const ;
    virtual int _last() const ;
    virtual int _minCapacity() const ;
    virtual int _maxCapacity() const ;
    virtual int _min() const ;
    virtual int _max() const ;
    virtual int _minPosAbs() const ; 
    virtual int _minNegAbs() const ; 
    virtual bool _remove(const int v);
    virtual bool _setDomain(const int v);
    virtual bool _setMin(const int l);
    virtual bool _setMax(const int u);
    virtual bool _setDomain(const BitSet& s);
    virtual bool _setEqual(BuildObject *v); 
    virtual bool _removeSet(const BitSet& s);
    virtual bool _removeSet(const int* vals, const int n);
    virtual bool _removeRange(const int l, const int u);

    virtual int _getNext(const int v) const ;
    virtual bool _good(const int v) const ;
    virtual bool _setNext(int& v) const ;

    virtual bool _isInterval() const ;
    virtual bool _isGround() const ;


    virtual bool setDomain(const int* a, const int l) { return true; }
    virtual bool setDomain(VariableInt* x) const { return true; }
    virtual void makeDecision(int&) {}
    virtual bool makeComplementary(int&) { return true; }


    virtual bool contain(const int v) const ;
    virtual bool equal(const int v) const ;
    virtual bool intersect(const int lo, const int up) const { return true; } // TODO
    virtual bool included(const int lo, const int up) const { return true; } // TODO
    virtual bool intersect(const BuildObject* x) const ;
    virtual bool included(const BuildObject* x) const ;
    virtual bool intersect(const BitSet* s) const { return true; } // TODO
    virtual bool included(const BitSet* s) const { return true; } // TODO
    virtual bool intersect(const BitSet& s) const { return true; } // TODO
    virtual bool intersect(const BitSet& s, int& r) const { return true; } // TODO
    virtual bool included(const BitSet& s) const { return true; } // TODO
    virtual void intersectTo( BitSet* s ) const {} // TODO
    virtual void unionTo( BitSet* s ) const {} // TODO
    virtual void intersectTo( BitSet& s ) const {} // TODO
    virtual void unionTo( BitSet& s ) const {} // TODO

    virtual bool _contain(const int v) const ;
    virtual bool _equal(const int v) const ;
    /*   virtual bool _intersect(const int lo, const int up) const ; */
    /*   virtual bool _included(const int lo, const int up) const ; */
    virtual bool _intersect(const BuildObject* x) const ; 
    virtual bool _included(const BuildObject* x) const ; 
    /*   virtual bool _intersect(const BitSet* s) const ; */
    /*   virtual bool _included(const BitSet* s) const ; */
    /*   virtual bool _intersect(const BitSet& s) const ; */
    /*   virtual bool _intersect(const BitSet& s, int& r) const ; */
    /*   virtual bool _included(const BitSet& s) const ; */
    /*   virtual void _intersectTo( BitSet* s ) const ; */
    /*   virtual void _unionTo( BitSet* s ) const ; */
    /*   virtual void _intersectTo( BitSet& s ) const ; */
    /*   virtual void _unionTo( BitSet& s ) const ; */

    virtual void printshort(std::ostream& o) const ;
    virtual void print_python()  ;
    virtual void print(std::ostream& o) const ;
    virtual std::string toString() const ;

    virtual std::string xmlPred( int& idx, int level, int *sub ) const ;
    virtual std::string xmlCon( BitSet& scope, int *index ) const ;
    virtual std::string xmlDom( ) const ;
    virtual void printstate( std::ostream& o ) const ;
    bool operator==( BuildObject& x ) const ;


//     virtual Mistral_Expression* ExpressionWrap() ;
// //     {
// //       Mistral_Expression *x = new Mistral_Expression(getMin(), getMax());
// //       return x;
// //     }

  };


  class BuildObjectPredicate;
  /**********************************************
   * Constraint BuildObject
   **********************************************/
  class BuildObjectConstraint 
  {

  public:

    BuildObjectConstraint() ;
    virtual ~BuildObjectConstraint() ;

    virtual const char* get_type(const BuildObjectPredicate *p) const {return "not implemented";}
  
    virtual void add( CSP *model, const int l, BuildObjectPredicate *pred ) ;
    virtual void build( Solver *s, VariableInt **tmp, BuildObjectPredicate *pred ) = 0;  
    virtual void self_build( Solver *s, VariableInt **tmp, BuildObjectPredicate *pred );  
    virtual void build( SatSolver *s, const int*idx, BuildObjectPredicate *pred );  
    virtual int propagateDownward( BuildObjectPredicate *pred ) const { return true; }
    virtual int propagateUpward( BuildObjectPredicate *pred ) const { return true; }
    virtual void close( BuildObjectPredicate *pred ) {}

    virtual std::string toString(const BuildObjectPredicate*) const ;
    virtual void print(std::ostream&, const BuildObjectPredicate*) const ;
    virtual std::string xmlPred(int&, int, const BuildObjectPredicate*, int*) const ;
  };



  /**********************************************
   * Predicate BuildObject
   **********************************************/
  class BuildObjectPredicate : public BuildObject
  {

  public:
 
    BuildObjectConstraint *relation;
    BuildObject **scope;
    int *params;

    // includes itself? NO
    int arity;
  
    BuildObjectPredicate(BuildObject **x, 
			 const int, const int, const int, 
			 BuildObjectConstraint*, int*); 
    virtual ~BuildObjectPredicate() ;

    virtual const char* get_type() const {return (relation ? relation->get_type(this) : "none");}
    virtual int get_arity() {return arity;}


    virtual bool isRel() const;
    virtual void unsetRel();

    virtual void build( Solver *s ) ;
    virtual void build( SatSolver *s, const int*idx ) ;
    virtual void add( CSP *model, const int l ) ;  
    virtual void close() ;
    virtual int propagateDownward( ) ; 
    virtual int propagateUpward( ) ; 

    virtual void printshort(std::ostream& o) const ;
    virtual std::string toString() const ;
    virtual void print(std::ostream& o) const ;
    virtual std::string xmlPred( int& idx, int level, int *sub ) const ;
    virtual std::string xmlCon( BitSet& scp, int *index ) const ;

  };




  /********************************************************
   *
   *  Unary Predicates BuildObjects & Wrappers
   *
   ********************************************************/

  /**********************************************
   *  Negation Predicate BuildObject
   **********************************************/ 
  class BuildObjectNegation : public BuildObjectConstraint {
 
  public:

    //virtual int state();
    //virtual int getMin( BuildObjectPredicate *pred ) const ;
    //virtual int getMax( BuildObjectPredicate *pred ) const ;

    virtual const char* get_type(const BuildObjectPredicate *p) const {return "neg";}
    virtual int propagateUpward( BuildObjectPredicate *pred ) const ;
    virtual int propagateDownward( BuildObjectPredicate *pred ) const ;
    virtual void close( BuildObjectPredicate *pred ) ;
    virtual void   build  (Solver *, VariableInt **, BuildObjectPredicate *) ;
    virtual void   print  (std::ostream&, const BuildObjectPredicate*) const ;
    virtual std::string xmlPred(int&, int, const BuildObjectPredicate*, int*) const ;
    virtual std::string toString(const BuildObjectPredicate*) const ;
    //virtual int computeValue( BuildObjectPredicate *pred ) const;
  };

  /**********************************************
   *  Not Predicate BuildObject
   **********************************************/ 
  class BuildObjectNot : public BuildObjectConstraint {
 
  public:

    //virtual int state();
    //virtual int getMin( BuildObjectPredicate *pred ) const ;
    //virtual int getMax( BuildObjectPredicate *pred ) const ;
    virtual const char* get_type(const BuildObjectPredicate *p) const {return "not";}
    virtual int propagateUpward( BuildObjectPredicate *pred ) const ;
    virtual int propagateDownward( BuildObjectPredicate *pred ) const ;
    virtual void close( BuildObjectPredicate *pred ) ;
    virtual void   build  (Solver *, VariableInt **, BuildObjectPredicate *) ;
    virtual void   build(SatSolver *, const int *, BuildObjectPredicate *);  
    virtual void   print  (std::ostream&, const BuildObjectPredicate*) const ;
    virtual std::string xmlPred(int&, int, const BuildObjectPredicate*, int*) const ;
    virtual std::string toString(const BuildObjectPredicate*) const ;
    //virtual int computeValue( BuildObjectPredicate *pred ) const;
  };

  /**********************************************
   *  Abs Predicate BuildObject
   **********************************************/ 
  class BuildObjectAbs : public BuildObjectConstraint {
 
  public:

    //virtual int state();
    //virtual int getMin( BuildObjectPredicate *pred ) const ;
    //virtual int getMax( BuildObjectPredicate *pred ) const ;
    virtual const char* get_type(const BuildObjectPredicate *p) const {return "abs";}
    virtual int propagateDownward( BuildObjectPredicate *pred ) const ;
    virtual int propagateUpward( BuildObjectPredicate *pred ) const ;
    virtual void close( BuildObjectPredicate *pred ) ;
    virtual void   build  (Solver *, VariableInt **, BuildObjectPredicate *) ;
    virtual void   print  (std::ostream&, const BuildObjectPredicate*) const ;
    virtual std::string xmlPred(int&, int, const BuildObjectPredicate*, int*) const ;
    virtual std::string toString(const BuildObjectPredicate*) const ;
    //virtual int computeValue( BuildObjectPredicate *pred ) const;
  };


  /********************************************************
   *
   *  Binary Predicates BuildObjects 
   *
   ********************************************************/

  //////////// LOGIC PREDICATES /////////////////
  /**********************************************
   * And Constraint BuildObject
   **********************************************/
  class BuildObjectAnd : public BuildObjectConstraint {

  public:

    //virtual int state();
    //virtual int getMin( BuildObjectPredicate *pred ) const ;
    //virtual int getMax( BuildObjectPredicate *pred ) const ;
    virtual const char* get_type(const BuildObjectPredicate *p) const {return "abs";}
    virtual int propagateUpward( BuildObjectPredicate *pred ) const ;
    virtual int propagateDownward( BuildObjectPredicate *pred ) const ;
    virtual void close( BuildObjectPredicate *pred ) ;
    virtual void   build  (Solver *, VariableInt **, BuildObjectPredicate *) ;
    virtual void   build(SatSolver *, const int *, BuildObjectPredicate *);  
    virtual void   print  (std::ostream&, const BuildObjectPredicate*) const ;
    virtual std::string xmlPred(int&, int, const BuildObjectPredicate*, int*) const ;
    virtual std::string toString(const BuildObjectPredicate*) const ;
    //virtual int computeValue( BuildObjectPredicate *pred ) const;
  };

  /**********************************************
   * Or Constraint BuildObject
   **********************************************/
  class BuildObjectOr : public BuildObjectConstraint {

  public:

    //virtual int state();
    //virtual int getMin( BuildObjectPredicate *pred ) const ;
    //virtual int getMax( BuildObjectPredicate *pred ) const ;
    virtual const char* get_type(const BuildObjectPredicate *p) const {return "or";}
    virtual int propagateUpward( BuildObjectPredicate *pred ) const ;
    virtual int propagateDownward( BuildObjectPredicate *pred ) const ;
    virtual void close( BuildObjectPredicate *pred ) ;
    virtual void   build  (Solver *, VariableInt **, BuildObjectPredicate *) ;
    virtual void   build(SatSolver *, const int *, BuildObjectPredicate *);  
    virtual void   print  (std::ostream&, const BuildObjectPredicate*) const ;
    virtual std::string xmlPred(int&, int, const BuildObjectPredicate*, int*) const ;
    virtual std::string toString(const BuildObjectPredicate*) const ;
    //virtual int computeValue( BuildObjectPredicate *pred ) const;
  };

  /**********************************************
   * IfThenElse Constraint BuildObject
   **********************************************/
  class BuildObjectIfThenElse : public BuildObjectConstraint {

  public:

    //virtual int state();
    //virtual int getMin( BuildObjectPredicate *pred ) const ;
    //virtual int getMax( BuildObjectPredicate *pred ) const ;
    virtual const char* get_type(const BuildObjectPredicate *p) const {return "ite";}
    virtual int propagateUpward( BuildObjectPredicate *pred ) const ;
    virtual int propagateDownward( BuildObjectPredicate *pred ) const ;
    virtual void close( BuildObjectPredicate *pred ) ;
    virtual void   build  (Solver *, VariableInt **, BuildObjectPredicate *) ;
    virtual void   print  (std::ostream&, const BuildObjectPredicate*) const ;
    virtual std::string xmlPred(int&, int, const BuildObjectPredicate*, int*) const ;
    virtual std::string toString(const BuildObjectPredicate*) const ;
    //virtual int computeValue( BuildObjectPredicate *pred ) const;
  };


  //////////// COMPARISON PREDICATES ////////////
  /**********************************************
   *  Disjunctive Predicate BuildObject
   **********************************************/ 
  class BuildObjectDisjunctive : public BuildObjectConstraint {
 
  public:

    virtual const char* get_type(const BuildObjectPredicate *p) const {return "dis";}
    virtual int propagateDownward( BuildObjectPredicate *pred ) const ;
    virtual int propagateUpward( BuildObjectPredicate *pred ) const ;
    virtual void close( BuildObjectPredicate *pred ) ;
    virtual void   build  (Solver *, VariableInt **, BuildObjectPredicate *) ;
    virtual void   print  (std::ostream&, const BuildObjectPredicate*) const ;
    virtual std::string xmlPred(int&, int, const BuildObjectPredicate*, int*) const ;
    virtual std::string toString(const BuildObjectPredicate*) const ;
  };

  //////////// COMPARISON PREDICATES ////////////
  /**********************************************
   *  Overlap Predicate BuildObject
   **********************************************/ 
  class BuildObjectOverlap : public BuildObjectConstraint {
 
  public:

    virtual const char* get_type(const BuildObjectPredicate *p) const {return "ovl";}
    virtual int propagateDownward( BuildObjectPredicate *pred ) const ;
    virtual int propagateUpward( BuildObjectPredicate *pred ) const ;
    virtual void close( BuildObjectPredicate *pred ) ;
    virtual void   build  (Solver *, VariableInt **, BuildObjectPredicate *) ;
    virtual void   print  (std::ostream&, const BuildObjectPredicate*) const ;
    virtual std::string xmlPred(int&, int, const BuildObjectPredicate*, int*) const ;
    virtual std::string toString(const BuildObjectPredicate*) const ;
  };

  /**********************************************
   *  Precedence Predicate BuildObject
   **********************************************/ 
  class BuildObjectPrecedence : public BuildObjectConstraint {
 
  public:

    //virtual int state();
    //virtual int getMin( BuildObjectPredicate *pred ) const ;
    //virtual int getMax( BuildObjectPredicate *pred ) const ;
    virtual const char* get_type(const BuildObjectPredicate *p) const {return "prec";}
    virtual int propagateUpward( BuildObjectPredicate *pred ) const ;
    virtual int propagateDownward( BuildObjectPredicate *pred ) const ;
    virtual void close( BuildObjectPredicate *pred ) ;
    virtual void   build  (Solver *, VariableInt **, BuildObjectPredicate *) ;
    //virtual void   build(SatSolver *, const int *, BuildObjectPredicate *);  
    virtual void   print  (std::ostream&, const BuildObjectPredicate*) const ;
    virtual std::string xmlPred(int&, int, const BuildObjectPredicate*, int*) const ;
    virtual std::string toString(const BuildObjectPredicate*) const ;
    //virtual int computeValue( BuildObjectPredicate *pred ) const;
  };

  /**********************************************
   *  Mul Predicate BuildObject
   **********************************************/
  class BuildObjectMul : public BuildObjectConstraint {
 
  public:

    //virtual int state();
    //virtual int getMin( BuildObjectPredicate *pred ) const ;
    //virtual int getMax( BuildObjectPredicate *pred ) const ;
    virtual const char* get_type(const BuildObjectPredicate *p) const {return "mul";}
    virtual void self_build( Solver *s, VariableInt **tmp, BuildObjectPredicate *pred );  
    virtual int propagateDownward( BuildObjectPredicate *pred ) const ;
    virtual int propagateUpward( BuildObjectPredicate *pred ) const ;
    virtual void close( BuildObjectPredicate *pred ) ;
    virtual void   build  (Solver *, VariableInt **, BuildObjectPredicate *) ;
    //virtual void   build(SatSolver *, const int *, BuildObjectPredicate *);  
    virtual void   print  (std::ostream&, const BuildObjectPredicate*) const ;
    virtual std::string xmlPred(int&, int, const BuildObjectPredicate*, int*) const ;
    virtual std::string toString(const BuildObjectPredicate*) const ;
    //virtual int computeValue( BuildObjectPredicate *pred ) const;
  };

  /**********************************************
   *  Div Predicate BuildObject
   **********************************************/
  class BuildObjectDiv : public BuildObjectConstraint {
 
  public:

    //virtual int state();
    //virtual int getMin( BuildObjectPredicate *pred ) const ;
    //virtual int getMax( BuildObjectPredicate *pred ) const ;
    virtual const char* get_type(const BuildObjectPredicate *p) const {return "div";}
    virtual void self_build( Solver *s, VariableInt **tmp, BuildObjectPredicate *pred );  
    virtual int propagateDownward( BuildObjectPredicate *pred ) const ;
    virtual int propagateUpward( BuildObjectPredicate *pred ) const ;
    virtual void close( BuildObjectPredicate *pred ) ;
    virtual void   build  (Solver *, VariableInt **, BuildObjectPredicate *) ;
    virtual void   print  (std::ostream&, const BuildObjectPredicate*) const ;
    virtual std::string xmlPred(int&, int, const BuildObjectPredicate*, int*) const ;
    virtual std::string toString(const BuildObjectPredicate*) const ;
    //virtual int computeValue( BuildObjectPredicate *pred ) const;
  };


  /**********************************************
   *  Mod Predicate BuildObject
   **********************************************/
  class BuildObjectMod : public BuildObjectConstraint {
 
  public:

    //virtual int state();
    //virtual int getMin( BuildObjectPredicate *pred ) const ;
    //virtual int getMax( BuildObjectPredicate *pred ) const ;
    virtual const char* get_type(const BuildObjectPredicate *p) const {return "mod";}
    virtual int propagateDownward( BuildObjectPredicate *pred ) const ;
    virtual int propagateUpward( BuildObjectPredicate *pred ) const ;
    virtual void close( BuildObjectPredicate *pred ) ;
    virtual void   build  (Solver *, VariableInt **, BuildObjectPredicate *) ;
    virtual void   print  (std::ostream&, const BuildObjectPredicate*) const ;
    virtual std::string xmlPred(int&, int, const BuildObjectPredicate*, int*) const ;
    virtual std::string toString(const BuildObjectPredicate*) const ;
    //virtual int computeValue( BuildObjectPredicate *pred ) const;
  };

  /**********************************************
   *  Equal Predicate BuildObject
   **********************************************/ 
  class BuildObjectEqual : public BuildObjectConstraint {
 
  public:

    //virtual int state();
    //virtual int getMin( BuildObjectPredicate *pred ) const ;
    //virtual int getMax( BuildObjectPredicate *pred ) const ;
    virtual const char* get_type(const BuildObjectPredicate *p) const {return (p->params[0] ? "eq" : "ne");}
    virtual int propagateUpward( BuildObjectPredicate *pred ) const ;
    virtual int propagateDownward( BuildObjectPredicate *pred ) const ;
    virtual void close( BuildObjectPredicate *pred ) ;
    virtual void   build  (Solver *, VariableInt **, BuildObjectPredicate *) ;
    virtual void   build(SatSolver *, const int *, BuildObjectPredicate *);  
    virtual void   print  (std::ostream&, const BuildObjectPredicate*) const ;
    virtual std::string xmlPred(int&, int, const BuildObjectPredicate*, int*) const ;
    virtual std::string toString(const BuildObjectPredicate*) const ;
    //virtual int computeValue( BuildObjectPredicate *pred ) const;
  };

  /**********************************************
   *  Member Predicate BuildObject
   **********************************************/ 
  class BuildObjectMember : public BuildObjectConstraint {
 
  public:

    //virtual int state();
    //virtual int getMin( BuildObjectPredicate *pred ) const ;
    //virtual int getMax( BuildObjectPredicate *pred ) const ;
    virtual const char* get_type(const BuildObjectPredicate *p) const {return "mem";}
    virtual int propagateUpward( BuildObjectPredicate *pred ) const ;
    virtual int propagateDownward( BuildObjectPredicate *pred ) const ;
    virtual void close( BuildObjectPredicate *pred ) ;
    virtual void   build  (Solver *, VariableInt **, BuildObjectPredicate *) ;
    virtual void   print  (std::ostream&, const BuildObjectPredicate*) const ;
    virtual std::string xmlPred(int&, int, const BuildObjectPredicate*, int*) const ;
    virtual std::string toString(const BuildObjectPredicate*) const ;
    //virtual int computeValue( BuildObjectPredicate *pred ) const;
  };

  /********************************************************
   *
   *  N-ary Predicates BuildObjects & Wrappers
   *
   ********************************************************/

  /**********************************************
   * Sum Constraint BuildObject
   **********************************************/
  class BuildObjectSum : public BuildObjectConstraint {

  public:

    //virtual int state();
    //virtual int getMin( BuildObjectPredicate *pred ) const ;
    //virtual int getMax( BuildObjectPredicate *pred ) const ;
    virtual const char* get_type(const BuildObjectPredicate *p) const {return "sum";}
    virtual int propagateUpward( BuildObjectPredicate *pred ) const ;
    virtual int propagateDownward(BuildObjectPredicate *pred) const ;
    virtual void close( BuildObjectPredicate *pred ) ;
    virtual void self_build( Solver *s, VariableInt **tmp, BuildObjectPredicate *pred );  
    virtual void   build  (Solver *, VariableInt **, BuildObjectPredicate *) ; 
    virtual void   print  (std::ostream&, const BuildObjectPredicate*) const ;
    virtual std::string xmlPred(int&, int, const BuildObjectPredicate*, int*) const ;
    virtual std::string toString(const BuildObjectPredicate*) const ;
    //virtual int computeValue( BuildObjectPredicate *pred ) const;
  };

  /**********************************************
   * Min Constraint BuildObject
   **********************************************/
  class BuildObjectMin : public BuildObjectConstraint {

  public:

    //virtual int state();
    //virtual int getMin( BuildObjectPredicate *pred ) const ;
    //virtual int getMax( BuildObjectPredicate *pred ) const ;
    virtual const char* get_type(const BuildObjectPredicate *p) const {return "min";}
    virtual int propagateUpward( BuildObjectPredicate *pred ) const ;
    virtual int propagateDownward(BuildObjectPredicate *pred ) const ;
    virtual void close( BuildObjectPredicate *pred ) ;
    virtual void   build  (Solver *, VariableInt **, BuildObjectPredicate *) ; 
    //virtual void   build(SatSolver *, const int *, BuildObjectPredicate *);  
    virtual void   print  (std::ostream&, const BuildObjectPredicate*) const ;
    virtual std::string xmlPred(int&, int, const BuildObjectPredicate*, int*) const ;
    virtual std::string toString(const BuildObjectPredicate*) const ;
    //virtual int computeValue( BuildObjectPredicate *pred ) const;
  };

  /**********************************************
   * Max Constraint BuildObject
   **********************************************/
  class BuildObjectMax : public BuildObjectConstraint {

  public:

    //virtual int state();
    //virtual int getMin( BuildObjectPredicate *pred ) const ;
    //virtual int getMax( BuildObjectPredicate *pred ) const ;
    virtual const char* get_type(const BuildObjectPredicate *p) const {return "max";}
    virtual int propagateUpward( BuildObjectPredicate *pred ) const ;
    virtual int propagateDownward( BuildObjectPredicate *pred ) const ;
    virtual void close( BuildObjectPredicate *pred ) ;
    virtual void   build  (Solver *, VariableInt **, BuildObjectPredicate *) ; 
    //virtual void   build(SatSolver *, const int *, BuildObjectPredicate *);  
    virtual void   print  (std::ostream&, const BuildObjectPredicate*) const ;
    virtual std::string xmlPred(int&, int, const BuildObjectPredicate*, int*) const ;
    virtual std::string toString(const BuildObjectPredicate*) const ;
    //virtual int computeValue( BuildObjectPredicate *pred ) const;
  };

  /**********************************************
   * Element Constraint BuildObject
   **********************************************/
  class BuildObjectElement : public BuildObjectConstraint {

  public:

    //virtual int state();
    //virtual int getMin( BuildObjectPredicate *pred ) const ;
    //virtual int getMax( BuildObjectPredicate *pred ) const ;
    virtual const char* get_type(const BuildObjectPredicate *p) const {return "Element";}
    virtual int propagateUpward( BuildObjectPredicate *pred ) const ;
    virtual int propagateDownward( BuildObjectPredicate *pred ) const ;
    virtual void close( BuildObjectPredicate *pred ) ;
    virtual void   build  (Solver *, VariableInt **, BuildObjectPredicate *) ; 
    virtual void   print  (std::ostream&, const BuildObjectPredicate*) const ;
    virtual std::string xmlPred(int&, int, const BuildObjectPredicate*, int*) const ;
    virtual std::string toString(const BuildObjectPredicate*) const ;
    //virtual int computeValue( BuildObjectPredicate *pred ) const;
  };

  /**********************************************
   * Cardinality Constraint BuildObject
   **********************************************/
  class BuildObjectCardinality : public BuildObjectConstraint {
 
  public:
 
    virtual const char* get_type(const BuildObjectPredicate *p) const {return "Cardinality";}
    virtual void add( CSP *model, const int l, BuildObjectPredicate *pred ) ;
    //virtual int state();
    //virtual int propagateUpward( BuildObjectPredicate *pred ) const ;
    //virtual int propagateDownward( BuildObjectPredicate *pred ) const ;
    //virtual void close( BuildObjectPredicate *pred ) ;
    virtual void   build  (Solver *, VariableInt **, BuildObjectPredicate *) ;
    virtual void   print  (std::ostream&, const BuildObjectPredicate*) const ; 
    virtual std::string xmlPred(int&, int, const BuildObjectPredicate*, int*) const ; 
    virtual std::string toString(const BuildObjectPredicate*) const ;
    //virtual int computeValue( BuildObjectPredicate *pred ) const;
 };


  /********************************************************
   *
   *  Pure Constraints BuildObjects
   *
   ********************************************************/

  /**********************************************
   * Clause BuildObject
   **********************************************/ 
  class BuildObjectClause : public BuildObjectConstraint {
 
  public:

    virtual const char* get_type(const BuildObjectPredicate *p) const {return "Clause";}
    virtual void close( BuildObjectPredicate *pred ) ;
    virtual void   build  (Solver *, VariableInt **, BuildObjectPredicate *) ;
    virtual void   print  (std::ostream&, const BuildObjectPredicate*) const ; 
    virtual std::string xmlPred(int&, int, const BuildObjectPredicate*, int*) const ; 
    virtual std::string toString(const BuildObjectPredicate*) const ;
  };


  /**********************************************
   *  Extensional Constraint BuildObject
   **********************************************/ 
  class BuildObjectTable : public BuildObjectConstraint {
 
  public:

    //Constraint *table;
    int spin;
    int arity;
    Vector<int*> tuples;

    std::vector< Constraint* > tables;
    std::vector< BuildObjectPredicate* > instance;
    std::vector< int > leader;
    int built_count;

    void addInstance( BuildObjectPredicate* pred );

    BuildObjectTable(const int a) ;
    BuildObjectTable(const int a, const int nbt) ;
    virtual ~BuildObjectTable() ;

    void printInstances();
    void selectType       ( BuildObjectPredicate *pred );
    void setNotEqual      ( BuildObjectPredicate *pred );

    virtual const char* get_type(const BuildObjectPredicate *p) const {return "Table";}
    virtual void close( BuildObjectPredicate *pred ) ;
    virtual int propagateUpward( BuildObjectPredicate *pred ) const ;
    virtual void   build  (Solver *, VariableInt **, BuildObjectPredicate *) ;
    virtual void   build(SatSolver *, const int *, BuildObjectPredicate *);  
    virtual void   print  (std::ostream&, const BuildObjectPredicate*) const ; 
    virtual std::string xmlPred(int&, int, const BuildObjectPredicate*, int*) const ; 
    virtual std::string toString(const BuildObjectPredicate*) const ;
    void add( const int* sol ) ;
  };

  /**********************************************
   * AllDiff Constraint BuildObject
   **********************************************/
  class BuildObjectAllDiff : public BuildObjectConstraint {

  public:

    virtual const char* get_type(const BuildObjectPredicate *p) const {return "AllDiff";}
    virtual void add( CSP *model, const int l, BuildObjectPredicate *pred ) ;
    virtual void close( BuildObjectPredicate *pred ) ;
    virtual void   build  (Solver *, VariableInt **, BuildObjectPredicate *) ; 
    virtual void   print  (std::ostream&, const BuildObjectPredicate*) const ; 
    virtual std::string xmlPred(int&, int, const BuildObjectPredicate*, int*) const ; 
    virtual std::string toString(const BuildObjectPredicate*) const ;
  };
 
  /**********************************************
   * Tree Constraint BuildObject
   **********************************************/
  class BuildObjectTree : public BuildObjectConstraint {

  public:

    //virtual void add( CSP *model, const int l, BuildObjectPredicate *pred ) ;
    //virtual void close( BuildObjectPredicate *pred ) ;
    virtual const char* get_type(const BuildObjectPredicate *p) const {return "Tree";}
    virtual void   build  (Solver *, VariableInt **, BuildObjectPredicate *) ;
    virtual void   print  (std::ostream&, const BuildObjectPredicate*) const ; 
    virtual std::string xmlPred(int&, int, const BuildObjectPredicate*, int*) const ; 
    virtual std::string toString(const BuildObjectPredicate*) const ;
  };

  /**********************************************
   * TDAG Constraint BuildObject
   **********************************************/
  class BuildObjectTDAG : public BuildObjectConstraint {

  public:

    //virtual void add( CSP *model, const int l, BuildObjectPredicate *pred ) ;
    virtual const char* get_type(const BuildObjectPredicate *p) const {return "TDAG";}
    virtual void close( BuildObjectPredicate *pred ) ;
    virtual void   build  (Solver *, VariableInt **, BuildObjectPredicate *) ;
    virtual void   print  (std::ostream&, const BuildObjectPredicate*) const ; 
    virtual std::string xmlPred(int&, int, const BuildObjectPredicate*, int*) const ; 
    virtual std::string toString(const BuildObjectPredicate*) const ;
  };

  /**********************************************
   * Gcc Constraint BuildObject
   **********************************************/
  class BuildObjectGcc : public BuildObjectConstraint {

  public:

    virtual const char* get_type(const BuildObjectPredicate *p) const {return "Gcc";}
    virtual void add( CSP *model, const int l, BuildObjectPredicate *pred ) ;
    virtual void close( BuildObjectPredicate *pred ) ;
    virtual void   build  (Solver *, VariableInt **, BuildObjectPredicate *) ;
    virtual void   print  (std::ostream&, const BuildObjectPredicate*) const ; 
    virtual std::string xmlPred(int&, int, const BuildObjectPredicate*, int*) const ; 
    virtual std::string toString(const BuildObjectPredicate*) const ;
  };

  /**********************************************
   * Lex Constraint BuildObject
   **********************************************/
  class BuildObjectLex : public BuildObjectConstraint {
 
  public:

    virtual const char* get_type(const BuildObjectPredicate *p) const {return "Lex";}
    virtual void add( CSP *model, const int l, BuildObjectPredicate *pred ) ;
    //virtual int state();
    //virtual void close( BuildObjectPredicate *pred ) ;
    virtual void   build  (Solver *, VariableInt **, BuildObjectPredicate *) ;
    virtual void   print  (std::ostream&, const BuildObjectPredicate*) const ; 
    virtual std::string xmlPred(int&, int, const BuildObjectPredicate*, int*) const ; 
    virtual std::string toString(const BuildObjectPredicate*) const ;
  };

  /**********************************************
   * Lex Constraint BuildObject
   **********************************************/
  class BuildObjectDLex : public BuildObjectConstraint {
 
  public:

    //virtual int state();  
    //virtual int getMin( BuildObjectPredicate *pred ) const ;
    //virtual int getMax( BuildObjectPredicate *pred ) const ;
    virtual const char* get_type(const BuildObjectPredicate *p) const {return "DLex";}
    virtual int propagateUpward(BuildObjectPredicate *pred) const ;
    virtual int propagateDownward(BuildObjectPredicate *pred) const ;
    virtual void close( BuildObjectPredicate *pred )  ;
    virtual void   build  (Solver *, VariableInt **, BuildObjectPredicate *) ;
    virtual void   print  (std::ostream&, const BuildObjectPredicate*) const ; 
    virtual std::string xmlPred(int&, int, const BuildObjectPredicate*, int*) const ; 
    virtual std::string toString(const BuildObjectPredicate*) const ;
  };


  /**********************************************
   *  Objectives declaration
   **********************************************/
  class BuildObjectObjective {
  public:
    BuildObject* X;
    virtual void build( Solver *s ) {}

    BuildObjectObjective(BuildObject *x) ;
    virtual ~BuildObjectObjective() {}

    virtual void print( std::ostream& o ) const ;
  };

  /// Maximisation declaration
  class BuildObjectMaximiseVar : public BuildObjectObjective {
  public:
    BuildObjectMaximiseVar(BuildObject *x) ;
    virtual ~BuildObjectMaximiseVar() {}
  
    virtual void build( Solver *s ) ;
    virtual void print( std::ostream& o ) const ;
  };

  /// Minimisation declaration
  class BuildObjectMinimiseVar : public BuildObjectObjective {
  public:
    BuildObjectMinimiseVar(BuildObject *x) ;
    virtual ~BuildObjectMinimiseVar() {}
  
    virtual void build( Solver *s ) ;
    virtual void print( std::ostream& o ) const ;
  };


  /**********************************************
   * Wrappers
   **********************************************/

  /**********************************************
   * Variable
   **********************************************/
  /*! \class Variable
    \brief Syntactic sugar used to hide the implementation of variables

    Arithmetic and logic C++ operators are overloaded
    so that an Ilog-like syntax is possible. 
  */
  class Variable
  {

  public:

    BuildObject *var_ptr_;
  
    Variable() ;
    Variable( const Variable& v ) ;
    Variable( BuildObject *v ) ;
    Variable( const int nvals ) ;
    Variable( const int lo, const int up ) ;
    Variable( const int* array, const int length ) ;
    Variable( const int* array, const int length, const int lo, const int up ) ;
    Variable( const BitSet& d ) ;
    Variable( const BitSet& d, const int length, const int lo, const int up ) ;
    virtual ~Variable() ;

    //VariableInt *getVarArray( VarArray& X );
    //VariableInt *getVariable();

    inline operator const BuildObject*() const
    {
      return var_ptr_;
    }


    int value() ;
    int size() ;
    int min() ;
    int max() ;

    /*!@name Operators*/
    Variable& operator= ( const Variable& x ) ;
    Variable operator+ ( Variable x ) ;
    Variable operator- ( Variable x ) ;
    Variable operator* ( Variable x ) ; 
    Variable operator/ ( Variable x ) ; 
    Variable operator% ( Variable x ) ; 
    Variable operator^ ( Variable x ) ;
    Variable operator+ ( const int k ) ;
    Variable operator- ( const int k ) ;
    Variable operator* ( const int k ) ; 
    Variable operator/ ( const int k ) ; 
    Variable operator% ( const int k ) ; 
    Variable operator^ ( const int k ) ; 
    Variable operator- ( ) ;
    Variable operator&&( Variable x ) ;
    Variable operator||( Variable x ) ; 
    Variable operator! ( ) ; 
    Variable operator==( Variable x ) ;
    Variable operator!=( Variable x ) ;
    Variable operator<=( Variable x ) ;
    Variable operator< ( Variable x ) ;
    Variable operator>=( Variable x ) ;
    Variable operator> ( Variable x ) ;
    Variable operator==( const BitSet& s ) ;
    Variable operator!=( const BitSet& s ) ;
    Variable operator==( const int k ) ;
    Variable operator!=( const int k ) ;
    Variable operator<=( const int k ) ;
    Variable operator< ( const int k ) ;
    Variable operator>=( const int k ) ;
    Variable operator> ( const int k ) ;
    VarArray operator,( Variable x ) ;

    VariableInt *getVariable() ;
    virtual void printshort(std::ostream& o) const ;
    virtual void print(std::ostream& o) const ;

  };


  /**********************************************
   * VarArray
   **********************************************/
  /*! \class VarArray
    \brief An array of variables (no wonder there)

  */
  class VarArray 
  {

  public:
  
#ifdef _NARRAY
    Vector<BuildObject*> array_;
    Variable buffer;
#else  
    std::vector<BuildObject*> tmp_;
    std::vector<Variable> array_;
#endif

    VarArray( ) ;
    VarArray( const int nvars ) ;
    VarArray( const int nvars, const int nvals ) ;
    VarArray( const int nvars, const int lo, const int up ) ;
    VarArray( const int nvars, const int* array, const int length ) ;
    VarArray( const int nvars, const int* array, const int length, const int lo, const int up ) ;

    virtual ~VarArray() ;

    BuildObject **getArgs( ) ;
    VariableInt **getVarArray( );

    virtual void build( Solver *s ) ;

    void resize( const int n ) ;
    void clear() ;
    int size() const ;
    void add( const int k );
    void add( Variable x );
    void add( VarArray& a );
  
    /*!@name Operators*/
    Variable& operator[]( const int i ) ;
    Variable operator[]( Variable x ) ;  
    VarArray operator,( Variable x ) ;
    VarArray operator,( VarArray& x ) ;

    Variable operator< ( VarArray& x ) ; 
    Variable operator> ( VarArray& x ) ; 
    Variable operator<=( VarArray& x ) ; 
    Variable operator>=( VarArray& x ) ;

    VarArray& operator=( VarArray& x ) ;

    virtual void printshort(std::ostream& o) const ;
    virtual void print(std::ostream& o) const ;


  };


  /**********************************************
   * Matrix
   **********************************************/
  /*! \class Matrix
    \brief A matrix of variables (no wonder there)

  */
  class Matrix 
  {

  public:

    int nRows;
    int nCols;
    Vector<BuildObject*> array_;
    VarArray tmp_;

    Matrix( const int nr, const int nc ) ;
    Matrix( const int nr, const int nc, const int nv ) ;
    Matrix( const int nr, const int nc, const int lo, const int up ) ;
    Matrix( const int nr, const int nc, const int* matrix, const int length ) ;
    Matrix( const int nr, const int nc, const int* matrix, const int length, const int lo, const int up ) ;

    virtual ~Matrix() ;

    //void addRow( VarArray& x );
    //void addColumn( VarArray& x );
  
    /*!@name Operators*/
    Variable cell( const int i, const int j ) ;
    VarArray& row( const int i ) ;
    VarArray& column( const int j ) ;


    virtual void printshort(std::ostream& o) const ;
    virtual void print(std::ostream& o) const ;
  };

  /**********************************************
   *  Negation Predicate Wrapper
   **********************************************/ 
  class Negation : public Variable {
  public:
    Negation( Variable& x_ );
  };
  /**********************************************
   *  Not Predicate Wrapper
   **********************************************/ 
  class Not : public Variable {
  public:
    Not( Variable& x_ );
  };
  /**********************************************
   *  Abs Predicate Wrapper
   **********************************************/ 
  class Abs : public Variable {
  public:
    Abs( Variable x_ );
  };
  /**********************************************
   *  And Predicate Wrapper
   **********************************************/ 
  class And : public Variable {
  public:
    And( Variable& x_, Variable& y_ );
  };
  /**********************************************
   *  Or Predicate Wrapper
   **********************************************/ 
  class Or : public Variable {
  public:
    Or( Variable& x_, Variable& y_ );
  };
  /**********************************************
   *  IfThenElse Predicate Wrapper
   **********************************************/ 
  class IfThenElse : public Variable {
  public:
    IfThenElse( Variable& x_, Variable& y_, Variable& z_ );
  };
  /**********************************************
   *  Disjunctive Predicate Wrapper
   **********************************************/ 
  class Disjunctive : public Variable {
  public:
    Disjunctive( Variable& x, int dx, Variable& y, int dy, const int tp=0 );
    //Disjunct( Variable& x, int dx, Variable& y, int dy );
  };
  /**********************************************
   *  Overlap Predicate Wrapper
   **********************************************/ 
  class Overlap : public Variable {
  public:
    Overlap( Variable& x, int dx, Variable& y, int dy );
    //Disjunct( Variable& x, int dx, Variable& y, int dy );
  };
  /**********************************************
   *  Precedence Predicate Wrapper
   **********************************************/ 
  class Precedence : public Variable {
  public:
    Precedence( Variable& x, const int d, Variable& y, const int tp=0 );
    Precedence( Variable& x, const int v );
    Precedence( const int v, Variable& x );
  };
  /**********************************************
   *  Mul Predicate Wrapper
   **********************************************/
  class Mul : public Variable {
  public:
    Mul( Variable& x, Variable& y );
    Mul( Variable& x, const int v );
  };
  /**********************************************
   *  Div Predicate Wrapper
   **********************************************/
  class Div : public Variable {
  public:
    Div( Variable& x, Variable& y );
  };
  /**********************************************
   *  Mod Predicate Wrapper
   **********************************************/
  class Mod : public Variable {
  public:
    Mod( Variable& x_, Variable& y_ );
    Mod( Variable& x_, const int m );
  };
  /**********************************************
   *  Equal Predicate Wrapper
   **********************************************/ 
  class Equal : public Variable {
  public:
    Equal( Variable& x_, Variable& y_, const int spin=1 );
    Equal( Variable& x_, const int v, const int spin=1 );
  };
  /**********************************************
   *  Member Predicate Wrapper
   **********************************************/ 
  class Member : public Variable {
  public:
    Member( Variable& x_, const BitSet& k, const int spin=1 );
  };
  /**********************************************
   * Sum Predicate Wrapper
   **********************************************/
  class Sum : public Variable {
  public:
    Sum( Variable& x, const int k );
    Sum( VarArray& x, const int *os=NULL );
    Sum( Variable& x_, Variable& y_, const int t=1 );
  };
  /**********************************************
   * Min Predicate Wrapper
   **********************************************/
  class Min : public Variable {
  public:
    Min( VarArray& x );
  };
  /**********************************************
   * Max Predicate Wrapper
   **********************************************/
  class Max : public Variable {
  public:
    Max( VarArray& x );
  };
  /**********************************************
   * Element Predicate Wrapper
   **********************************************/
  class Element : public Variable {
  public:
    Element( VarArray& x );
    Element( VarArray& x, Variable& y );
  };
  // /**********************************************
  //  *  Clause Wrapper
  //  **********************************************/ 
  // class BClause : public Variable {
  // public:
  //   BClause( VarArray& x, const int* c );
  // };
  /**********************************************
   *  Extensional Constraint Wrapper
   **********************************************/ 
  class Table : public Variable {
  private:
    BuildObjectTable *tabptr_;

  public:
    Table();
    Table( VarArray& x, const bool sp=true );
    Table( VarArray& x, const int tidx, const bool sp );
  
    void add( const int* sol );
  };
  /**********************************************
   * AllDifferent Constraint Wrapper
   **********************************************/
  class AllDifferent : public Variable {
  public:
    AllDifferent( VarArray& x, const int c=3, const int k=NOVAL/2 );
  };

  /**********************************************
   * Gcc Constraint Wrapper
   **********************************************/
  class Gcc : public Variable {
  public:
    Gcc( VarArray& x,
	 const int firstDomainValue, 
	 const int lastDomainValue,
	 const int* minOccurrences,
	 const int* maxOccurrences );
  };
  /**********************************************
   * TDAG Constraint Wrapper
   **********************************************/
  class TransitiveDAG : public Variable {
  public:
    TransitiveDAG( VarArray& x );
  };
  /**********************************************
   * Tree Constraint Wrapper
   **********************************************/
  class Tree : public Variable {
  public:
    Tree( VarArray& x );
    Tree( VarArray& x, VarArray& y );
  };
  /**********************************************
   * Lex < Constraint Wrapper
   **********************************************/
  class LexOrder : public Variable {
  public:
    LexOrder( VarArray& x, VarArray& y, int eq=1 );
  };
  /**********************************************
   * Cardinality < Constraint Wrapper
   **********************************************/
  class Cardinality : public Variable {
  public:
    Cardinality( VarArray& x, const int k );
  };
  // /**********************************************
  //  * MinOccurrence < Constraint Wrapper
  //  **********************************************/
  // class MinOccurrence : public Variable {
  // public:
  //   MinOccurrence( VarArray& x, const int k );
  // };



  class Objective {
  public:
    BuildObject *X;
    BuildObjectObjective *obj;
  
    Objective(BuildObject *x);// : X(x) { obj = NULL; }
  };

  /// Maximisation wrapper
  class Maximise : public Objective {
  public:
    Maximise(Variable x) ;
  };
  /// Minimisation wrapper
  class Minimise : public Objective {
  public:
    Minimise(Variable x) ;
  };




  /**********************************************
   * CSP
   **********************************************/
  /*! \class CSP
    \brief Implementation of a Constraint network.

    This implementation is composed
    of an array of pointers to Variables, and
    a list of pointers to Constraints
  */
  class CSP {

  public:
  
    //static ConstraintStore ENVIRONMENT;

    /*!@name Parameters*/
    //@{

    ConstraintClauseBase *clauseBase;
    bool satCompatible;

    // stats
    unsigned int nList;
    unsigned int eList;
    unsigned int nBit;
    unsigned int eBit;
    unsigned int nRange;
    unsigned int eRange;
    unsigned int nBoolean;
    unsigned int eBoolean;
    unsigned int nConstant;
    unsigned int eConstant;
    unsigned int nValues;
    unsigned int eValues;
    unsigned int nConstraint;
    unsigned int eConstraint;
    unsigned int nNeq;

    // objective function
    BuildObjectObjective *goald; 

    // Model declaration
    Vector<BuildObject*> declarations; 
    MistralGacList<BuildObject*> buildqueue;

    Vector<BuildObject*> toplevel_expressions; 
    List variables;

    // Infered not-equals
    Vector<BuildObject*> notequals;


    bool closed;
    // 
    //Vector< BuildObject** > cl_scope;
    //Vector< int* > cl_lit;
    //Vector< int > cl_size;
    //@}

    /*!@name Constructors*/
    //@{
    CSP() ;
    virtual ~CSP() ; 
   
    void cleanUp() ;
    //@}


    int inferAllDiffs( const int all, const int limit );
    inline void triggerEvent(int id, int e) ;
    bool isSatCompatible() const ;
    void unsetSat() ;
    void setGAC( const int g );

  
    /*!@name Accessors*/
    //@{
    int size() ;
    int numVars() const ;
    int numEVars() const ;
    /// Add a single variable x 
    void add( BuildObject *x ) ;
    void add( BuildObjectObjective *x ) ;
    void add( BuildObject **x, const int n ) ;

    void add( Variable x ) ;
    /// Add an array of variables 
    void add( VarArray& x ) ;
    /// Add an objective function g 
    void add( Objective g ) ;
    /// Build the data structure used during search;
    void build( Solver *s ) ;
    void build( SatSolver *s ) ;
    bool preBuild( ) ;
    void close( ) ;
    //@}


    /*!@name Utils*/
    //@{
    static BuildObject* _Variable( const int nv );
    static BuildObject* _Variable( const int lb, const int ub );
    static BuildObject* _Variable( const int *array, const int l );
    static BuildObject* _Variable( const int *array, const int l,
				   const int lb, const int ub );
    static BuildObject* _Variable( const BitSet& d, const int l,
				   const int lb, const int ub );


    static BuildObject** _VarArray( const int n, const int nv );
    static BuildObject** _VarArray( const int n, const int lb, const int ub );
    static BuildObject** _VarArray( const int n, const int *array, const int length );
    static BuildObject** _VarArray( const int n, const int *array, const int length, 
				    const int lb, const int ub );
  
    static BuildObject* _Constraint( const int t, BuildObject **x, const int n, int*p );

    static BuildObject* _Negation( BuildObject *x );
    static BuildObject* _Not( BuildObject *x );
    static BuildObject* _Abs( BuildObject *x );
    static BuildObject* _And( BuildObject *x, BuildObject *y );
    static BuildObject* _Or( BuildObject *x, BuildObject *y );
    static BuildObject* _IfThenElse( BuildObject *x, BuildObject *y, BuildObject *z );
    static BuildObject* _Disjunctive( BuildObject *x, const int dx, BuildObject *y, int const dy, const int t );
    static BuildObject* _Overlap( BuildObject *x, const int dx, BuildObject *y, int const dy );
    static BuildObject* _Precedence( BuildObject *x, const int d, BuildObject *y );
    static BuildObject* _Precedence( const int d, BuildObject *y );
    static BuildObject* _Precedence( BuildObject *x, const int d );
    static BuildObject* _Add( BuildObject *x, BuildObject *y );
    static BuildObject* _Sub( BuildObject *x, BuildObject *y );
    static BuildObject* _Mul( BuildObject *x, BuildObject *y );
    static BuildObject* _Mul( BuildObject *x, const int v );
    static BuildObject* _Div( BuildObject *x, BuildObject *y );
    static BuildObject* _Mod( BuildObject *x, BuildObject *y );
    static BuildObject* _Mod( BuildObject *x, const int m );
    static BuildObject* _Equal( BuildObject *x, BuildObject *y, const int spin );
    static BuildObject* _Equal( BuildObject *x, const int v, const int spin );
    static BuildObject* _Member( BuildObject *x, const BitSet& k, const int spin );
    static BuildObject* _Clause( BuildObject **x, const int nx, const int *pol );
    static BuildObject* _Sum( BuildObject **x, const int nx, const int *os=NULL );
    static BuildObject* _Sum( BuildObject *x, const int o );
    static BuildObject* _Sum( BuildObject *x, BuildObject *y, const int t );
    static BuildObject* _Min( BuildObject *x, BuildObject *y );
    static BuildObject* _Max( BuildObject *x, BuildObject *y );
    static BuildObject* _Min( BuildObject **x, const int n );
    static BuildObject* _Max( BuildObject **x, const int n );
    static BuildObject* _Element( BuildObject **x, const int nx, const int k );
    static BuildObject* _Tree( BuildObject **x, const int nx );
    static BuildObject* _TDAG( BuildObject **x, const int nx );
    static BuildObject* _AllDifferent( BuildObject **x, const int nx );
    static BuildObject* _AllDifferent( BuildObject **x, const int nx, const int c, const int k );
    static BuildObject* _Gcc( BuildObject **x, const int nx,
			      const int firstDomainValue, 
			      const int lastDomainValue,
			      const int* minOccurrences,
			      const int* maxOccurrences );
    static BuildObject* _LexOrder( BuildObject **x, const int nx, const int eq=1 );
    static BuildObject* _Cardinality( BuildObject **x, const int nx, const int k );

    static BuildObject* _Table( BuildObject **x, const int nx, const bool sp=true );
    static BuildObject* _Table( BuildObject **x, const int nx, const int tidx, const bool sp );
    static void _addTuple( BuildObject *tab, const int* sol );
			      
    static int addTable( const int a, const int n );
    static void addTuple( const int idx, int *tuple );

    static BuildObjectObjective* _Maximise( BuildObject *x );
    static BuildObjectObjective* _Minimise( BuildObject *x );
    //@}
  
    /*!@name Miscellaneous*/
    //@{
    /// Printshort the CSP
    void print(std::ostream&) const ;
    void printPython() const ;
    void printXML(std::ostream&) const ;
    //@}
  };



  int registerConstraint(BuildObjectConstraint*);
  BuildObjectConstraint* getConstraint(const int);

};



std::ostream& operator<< (std::ostream& os, const Mistral::MistralGraph& g) ;
std::ostream& operator<< (std::ostream& os, const Mistral::Variable& x) ;
std::ostream& operator<< (std::ostream& os, const Mistral::VarArray& x) ;
std::ostream& operator<< (std::ostream& os, const Mistral::Matrix& x) ;
std::ostream& operator<< (std::ostream& os, const Mistral::CSP& p) ;



#endif // __Wrapper_H



