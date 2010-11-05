
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

/** \file con.h
    \brief Header for the Constraint library.
*/

/**********************************************
 * Constraint and subclasses
 **********************************************/


#ifndef __CONSTRAINT_H
#define __CONSTRAINT_H

#include <mistral_csp.h>
#include <mistral_var.h>
#include <mistral_sat.h>

#ifdef _COINLP

#include "ClpSimplex.hpp"
#include "ClpInterior.hpp"
#include "CoinHelperFunctions.hpp"
#include "CoinTime.hpp"
#include "CoinBuild.hpp"
#include "CoinModel.hpp"
#include "ClpCholeskyBase.hpp"

#endif


#ifdef _CPLEX

#include <ilcplex/cplex.h>

#endif


#include <assert.h>



namespace Mistral {

  typedef int Literal;
  typedef int Atom;
  typedef Array<Literal> Clause;
  
  /**********************************************
   * Top level Constraints
   **********************************************/


  /**********************************************
   * Unary Constraint
   **********************************************/
  /*! \class UnaryConstraint
    \brief  Unary Constraints.
  */
  class UnaryConstraint {

  public:
    /**@name Parameters*/
    //@{ 
    const char* name;
    Solver *solver;
    VariableInt *X;
    //@}

    /**@name Constructors*/
    //@{
    UnaryConstraint(Solver *s, VariableInt *v, const char* n);
    virtual ~UnaryConstraint();
    //@}

    /**@name Solving*/
    //@{
    virtual bool propagate() = 0;
    void activate();
    //@}

    /**@name Miscellaneous*/
    //@{
    virtual void print(std::ostream& o) const = 0;
    //@}
  };


  /**********************************************
   * Unary Lower Bound Constraint
   **********************************************/
  /*! \class UnaryConstraintMore
    \brief  Unary Constraints (lower bound).
  */
  class UnaryConstraintMore : public UnaryConstraint {

  public:
    /**@name Parameters*/
    //@{ 
    int bound;
    //@}

    /**@name Constructors*/
    //@{
    UnaryConstraintMore(Solver*, VariableInt*, int);
    virtual ~UnaryConstraintMore() {};
    //@}

    /**@name Solving*/
    //@{
    virtual bool propagate();
    //@}

    /**@name Miscellaneous*/
    //@{
    virtual void print(std::ostream& o) const;
    //@}
  };


  /**********************************************
   * Unary Lower Bound Constraint
   **********************************************/
  /*! \class UnaryConstraintLess
    \brief  Unary Constraints (upper bound).
  */
  class UnaryConstraintLess : public UnaryConstraint {

  public:
    /**@name Parameters*/
    //@{ 
    int bound;
    //@}

    /**@name Constructors*/
    //@{
    UnaryConstraintLess(Solver*, VariableInt*, int);
    virtual ~UnaryConstraintLess() {};
    //@}

    /**@name Solving*/
    //@{
    virtual bool propagate();
    //@}

    /**@name Miscellaneous*/
    //@{
    virtual void print(std::ostream& o) const;
    //@}
  };



  /**********************************************
   * Clause Constraint
   **********************************************/ 
  /*! \class ConstraintClause
    \brief  implementation of a disjunction of literals "x =/= v"
  */
  class ConstraintClause : public Constraint {

  public:
    /**@name Parameters*/
    //@{ 
    int* conflict;
    //@}

    /**@name Constructors*/
    //@{
    ConstraintClause(Solver*, VariableInt** v, const int n, const int* c);
    virtual ~ConstraintClause();
    //@}

    /**@name Solving*/
    //@{
    inline int check( const int* ) const ;  
    inline bool propagate();
    //  inline bool propagate(VariableInt *v, const int e);
    inline bool propagate(const int changedIdx, const int e);
    //@}

    /**@name Miscellaneous*/
    //@{
    virtual void print(std::ostream& o) const;
    //@}
  };


  /**********************************************
   * Extensional Constraint
   **********************************************/ 
  /*! \class ConstraintGAC3Valid
    \brief  implementation of a n-ary extensional constraint

    The conflicts are encoded using
    a N-dimensional flatened matrix:
    an entry of the matrix equal to 0
    denotes an allowed combination (support), whilst
    1 denotes a forbidden combination (conflict).
  */
  class ConstraintGAC3Valid : public Constraint {

  private:
    /**@name Parameters*/
    //@{ 
    int *var_sizes; // NOTE: this has arity-1 elements
    BitSet matrix;
    bool isClone;
    //@}

  public:
    /**@name Constructors*/
    //@{
    //   ConstraintGAC3Valid(Solver*, VariableInt** v, 
    // 		      const int l, BitSet& m, 
    // 		      int*);
    ConstraintGAC3Valid(Solver* s, VariableInt** v, const int l);
    void init( Vector<int*>& tuples, const int spin, ConstraintGAC3Valid *con );
    virtual ~ConstraintGAC3Valid();
    //@}

    /**@name Accessors*/
    //@{
    void fillNogood () ;
    void fillSupport() ;

    void addNogood(const int* vals) ;
    void addSupport(const int* vals) ;

    bool isValid(const int *tuple) const;
    int getpos(const int *vals) const;
    int getpos() const;
    //@}

    /**@name Solving*/
    //@{
    /*!
      Check the corresponding cell in the matrix
    */
    inline int check( const int* ) const ;
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual void print(std::ostream& o) const;
    //@}
  };




  /**********************************************
   * GAC2001 Extensional Constraint
   **********************************************/ 
  /*! \class ConstraintGAC2001Allowed
    \brief  implementation of a n-ary extensional constraint

    The conflicts are encoded using
    a list of supporting tuples for each
    value of each variable.
    This encoding is less time efficient
    but more space efficient than the matrix encoding
  */
  class ConstraintGAC2001Allowed : public Constraint {

  private:
    /**@name Parameters*/
    //@{ 
    ReversibleNum<int> **firstSupport;
    Vector< int* > **supportList;
    int *order;
    BitSet *D_X;
    int** supports_X;
    int* themins;
    bool isClone;
    //@}

  public:
    /**@name Constructors*/
    //@{
    //   ConstraintGAC2001Allowed(Solver*, VariableInt** v, 
    // 			   const int n, Vector<int*>** sl);
    ConstraintGAC2001Allowed(Solver* s, VariableInt** v, const int l);
    void init( Vector<int*>& tuples, const int spin, ConstraintGAC2001Allowed *con );
    virtual ~ConstraintGAC2001Allowed();
    //@}

    /**@name Solving*/
    //@{
    /*!
      Check the corresponding cell in the matrix
    */
    inline int check( const int* ) const ;
    inline bool propagate(const int changedIdx, const int e);
    //@}

    /**@name Miscellaneous*/
    //@{ 
    int getNumSupports( const int var, const int val ) const;
    virtual void print(std::ostream& o) const;
    virtual void printlong(std::ostream& o) const;
    virtual std::string name() { return "GAC2001"; }
    //@}
  };


  /**********************************************
   * Extensional Binary Constraint
   **********************************************/ 
  /*! \class ConstraintAC3Bitset
    \brief  implementation of a binary extensional constraint

    Uses bit-wise operation to speed-up support checking.
    Lists of supports are stored in bitsets
  */
  class ConstraintAC3Bitset : public Constraint {

  private:

    /**@name Parameters*/
    //@{ 
    BitSet *supportList[2];
    int size[2];
    int minSup[2];
    int *residualSupport[2];
    bool isClone;
    //@}

  public:
    /**@name Constructors*/
    //@{
    //ConstraintAC3Bitset(Solver*, VariableInt**, BitSet**);
    ConstraintAC3Bitset(Solver* s, VariableInt** v, const int l);
    void init( Vector<int*>& tuples, const int spin, ConstraintAC3Bitset *con );
    virtual ~ConstraintAC3Bitset();
    void init();  
    //@}

    /**@name Accessors*/
    //@{
    void fillNogood () ;
    void fillSupport() ;

    void addNogood(const int* vals) ;
    void addSupport(const int* vals) ;
    //@}

    /**@name Solving*/
    //@{
    /*!
      Check the corresponding cell in the matrix
    */
    inline int check( const int* ) const ;
    inline bool propagate(const int changedIdx, const int e);
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual void print(std::ostream& o) const;
    //@}
  };


  /**********************************************
   * =/= Constraint
   **********************************************/ 
  /*! \class ConstraintNotEqual
    \brief  Binary not-equal Constraint (x0 =/= x1).
  */
  class ConstraintNotEqual : public Constraint {

  public:  
    /**@name Constructors*/
    //@{
    ConstraintNotEqual(Solver*, VariableInt** v);
    virtual ~ConstraintNotEqual();
    //@}

    /**@name Solving*/
    //@{
    inline int check( const int* ) const ;
    inline bool propagate();
    inline bool propagate(const int changedIdx, const int e);
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual void print(std::ostream& o) const;
    //@}
  };


  /**********************************************
   * =/= Constraint with a wildcard
   **********************************************/ 
  /*! \class ConstraintNotEqualWildcard
    \brief  Binary not-equal Constraint with a wildcard
    (x0 =/= x1 or x0 == w).
  */
  class ConstraintNotEqualWildcard : public Constraint {

  public:  
    /**@name Parameters*/
    //@{ 
    int wildcard;
    //@}

    /**@name Constructors*/
    //@{
    ConstraintNotEqualWildcard(Solver*, VariableInt**, const int);
    virtual ~ConstraintNotEqualWildcard();
    //@}

    /**@name Solving*/
    //@{
    inline int check( const int* ) const ;
    inline bool propagate();
    inline bool propagate(const int changedIdx, const int e);
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual void print(std::ostream& o) const;
    //@}
  };


  /**********************************************
   * Nogood Constraint 
   **********************************************/ 
  /*! \class ConstraintNogood
    \brief  Nogood base

  */
  class GenLiteral {
  public:
    VariableInt* var;
    int val;
    
    GenLiteral(VariableInt*x) {
      var = x;
      val = x->min();
    }
    GenLiteral(VariableInt*x, int v) {
      var = x;
      val = v;
    }
    GenLiteral() {
      var = NULL;
      val = NOVAL;
    }
    virtual ~GenLiteral() {}

    bool satisfied() { return !var->contain(val); }
    bool violated() { return var->equal(val); }

    void devirtualize() {
      // when the variable is virtual, override the current literal by a non-virtual one
      if(var->getType() == VariableInt::VIRTUAL) {

// 	std::cout << "DE-VIRTUALIZE ";
// 	print(std::cout);
// 	std::cout << " ==> ";

	VariableInt *real_var = var->getVar();
	int real_val = ((VariableVirtual*)var)->conversion->backward(val);
	var = real_var;
	val = real_val;

// 	print(std::cout);
// 	std::cout << std::endl;
      }
    }


    void print(std::ostream& o) const {
      var->printshort(o);
      o << " =/= " << val; 
    }
  };

  class ConstraintNogoodBase : public Constraint {

  public:  
    /**@name Parameters*/
    //@{ 
    //VariableInt **scp;
    // minimum values, used as an offset when accessing the base
    int *minimums;
    // list of nogoods
    Vector< Array < GenLiteral >* > nogood;
    // the watched literals data structure
    Vector< Array < GenLiteral >* > **watched;
    //@}

    /**@name Constructors*/
    //@{
    ConstraintNogoodBase(Solver*);
    virtual ~ConstraintNogoodBase();
    void add( Vector<GenLiteral>& );
    //@}

    /**@name Solving*/
    //@{
    inline int check( const int* ) const ;
    inline bool propagate();
    inline bool propagate(const int changedIdx, const int e);
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual void print(std::ostream& o) const;
    //@}
  };


  class SimpleUnaryConstraint;
  // can be upper/lower bound or equality/disequality
  typedef SimpleUnaryConstraint GeneralLiteral; 
  class ConstraintGenNogoodBase : public Constraint {

  public:  
    /**@name Parameters*/
    //@{ 

    /// Attributes:
    /// - A list of (pointer to) clauses per variable/value (the "watched literal")
    /// - The (flat) list of clauses

    // minimum values, used as an offset when accessing the base
    int *minimums;
    // list of nogoods
    Vector< Array < GeneralLiteral >* > nogood;
    // the watched literals data structure
    Vector< Array < GeneralLiteral >* > **watched_values;
    Vector< Array < GeneralLiteral >* > *watched_domains;
    Vector< Array < GeneralLiteral >* > *watched_bounds;
    //@}

    /**@name Constructors*/
    //@{
    ConstraintGenNogoodBase(Solver*);
    virtual ~ConstraintGenNogoodBase();
    void add( Vector<GeneralLiteral>& );
    //@}

    /**@name Solving*/
    //@{
    inline int check( const int* ) const ;
    inline bool propagate();
    inline bool propagate(const int changedIdx, const int e);
    //@}

    /**@name Miscellaneous*/
    //@{  
    void forget(const int);
    void removeClause(const int);
    virtual void print(std::ostream& o) const;
    //@}
  };

  /**********************************************
   * Inverse Constraint Decomposition
   **********************************************/ 
  /*! \class ConstraintInverse
    \brief  Binary inverse Constraint (xi = j <-> xj = i)

    Use to encode primal/dual channeling
  */
  class ConstraintInverse : public Constraint {

  public:  
    /**@name Parameters*/
    //@{ 
    int i_[2];
    //@}

    /**@name Constructors*/
    //@{
    ConstraintInverse(Solver*, VariableInt** v, int, int);
    virtual ~ConstraintInverse();
    //@}

    /**@name Solving*/
    //@{
    inline int check( const int* ) const ;
    inline bool propagate(const int changedIdx, const int e);
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual void print(std::ostream& o) const;
    //@}
  };


  /**********************************************
   * Lex Constraint Decomposition
   **********************************************/ 
  /*! \class ConstraintLex
    \brief  Quaternary constraint used to encode LexOrder
  */
  class ConstraintLex : public Constraint {

  public:  

    //ReversibleBool& domain_b1;
    //ReversibleBool& domain_b2;
    int& domain_b1;
    int& domain_b2;

    /**@name Constructors*/
    //@{
    ConstraintLex(Solver*, VariableInt** v);
    virtual ~ConstraintLex();
    //@}

    /**@name Solving*/
    //@{
    inline bool rule1forward();
    inline bool rule2forward();
    //  inline bool rule3forward();
    inline bool rule1backward();
    inline bool rule2backward();
    //  inline bool rule3backward();
    virtual int check( const int* ) const ;
    virtual bool propagate(const int idx, const int e);
    //@}

    /**@name Miscellaneous*/
    //@{    
    virtual void print(std::ostream& o) const ;
    //@}
  };


  /**********************************************
   * <= Constraint
   **********************************************/ 
  /*! \class ConstraintLess
    \brief  Binary Less Than Constraint (x0 + k <= x1).
  */
  class ConstraintLess : public Constraint {

  public: 
    /**@name Parameters*/
    //@{  
    int offset;
    //@}

    /**@name Constructors*/
    //@{
    ConstraintLess(Solver*, VariableInt** v);
    ConstraintLess(Solver*, VariableInt** v, int ofs);
    virtual ~ConstraintLess();
    //@}

    /**@name Solving*/
    //@{
    inline int check( const int* ) const ;
    inline bool propagate(const int changedIdx, const int e);
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual void print(std::ostream& o) const;
    //@}
  };


  // /***********************************************
  //  * All Different Constraint (arc consistency)
  //  ***********************************************/

  // /*! \class ConstraintAllDiff
  //  \brief  AllDifferent Constraint.

  //  Constraint of difference on a set of variables.
  //  Ensures that a set of variables are assigned distinct 
  //  values 
  // */
  // class ConstraintAllDiff : public Constraint {
  
  //  private:
  //   /**@name Parameters*/
  //   //@{  
  //   // store back edges for the scc and also the matching
  //   int *matching; // value to which variable i is matched
  //   int *back; // back[i] = -1 => value i is not matched
  //   BitSet freeValues; // values that do not participate in the matching
  //   //@}

  //   // need to implement:
  //   // 


  //  public:
  //   /**@name Constructors*/
  //   //@{
  //   ConstraintAllDiff(Solver*, VariableInt** x, const int n);
  //   virtual ~ConstraintAllDiff();
  //   //@}

  //   /**@name Solving*/
  //   //@{
  //   inline int check( const int* ) const ;
  //   inline bool propagate();
  //   inline bool propagate(const int changedIdx, const int e); 
  //   //@}

  //   /**@name Miscellaneous*/
  //   //@{  
  //   virtual void print(std::ostream& o) const ;
  //   //@}
  // };


  /***********************************************
   * All Different Constraint (bounds consistency).
   ***********************************************/
  typedef struct {
    int min, max;		// start, end of interval
    int minrank, maxrank; // rank of min & max in bounds[] of an adcsp
  } Interval;

  /*! \class ConstraintAllDiff
    \brief  AllDifferent Constraint.

    Constraint of difference on a set of variables.
    Ensures that a set of variables are assigned distinct 
    values (only Bounds Consistency is implemented)
    The code is from Lopez-Ortiz, Quimper, Tromp and van Beek
  */
  class ConstraintAllDiff : public Constraint {
  
  private:
    /**@name Parameters*/
    //@{  
    int lastLevel;
    int *t;		// tree links
    int *d;		// diffs between critical capacities
    int *h;		// hall interval links
    Interval *iv;
    Interval **minsorted;
    Interval **maxsorted;
    int *bounds;  // bounds[1..nb] hold set of min & max in the niv intervals
    // while bounds[0] and bounds[nb+1] allow sentinels
    int nb;
    void sortit();
    int filterlower();
    int filterupper();
    void propagateValue();
    //@}

  public:
    /**@name Constructors*/
    //@{
    ConstraintAllDiff(Solver*, VariableInt** x, const int n);
    virtual ~ConstraintAllDiff();
    void init();
    //@}

    /**@name Solving*/
    //@{
    inline int check( const int* ) const ;
    inline bool propagate();
    inline bool propagate(const int changedIdx, const int e); 
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual void print(std::ostream& o) const ;
    //@}
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
   

  /*! \class ConstraintGlobalCardinality
    \brief  Global Cardinality Constraint.

    User defined propagator for enforcing bounds consistency
    on the restricted gcc constraint when bounds on
    occurrences are [a_i,b_i].
    A value "v" must be assigned to at least
    minOccurrences[v - firstDomainValue] variables and at most
    maxOccurrences[v - firstDomainValue] variables
    The code is from Lopez-Ortiz, Quimper, Tromp and van Beek
  */   
  class ConstraintGlobalCardinality : public Constraint
  {

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
    Interval *iv;
    Interval **minsorted;
    Interval **maxsorted;
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
    //@}  

  public:
    /**@name Constructors*/
    //@{
    ConstraintGlobalCardinality( Solver *s,
				 VariableInt **v,
				 const int n,
				 const int firstDomainValue,
				 const int lastDomainValue,
				 const int* minOccurrences,
				 const int* maxOccurrences);
    ~ConstraintGlobalCardinality();
    //@}

    /**@name Solving*/
    //@{  
    inline int check( const int* ) const ;
    inline bool propagate();
    inline bool propagate(const int changedIdx, const int e); 
    //@}

    /**@name Miscellaneous*/
    //@{    
    virtual void print(std::ostream& o) const ;
    //@}
  };


  /**********************************************
   * Transitive DAG Constraint
   **********************************************/ 
  /*! \class ConstraintTDAG
    \brief  TDAG Constraint.

  */
  class ConstraintTDAG : public Constraint {
  
  private:
    /**@name Parameters*/
    //@{  
    // 
    int num_nodes;
    int *node_index;

    VariableInt **intervals;
    int *first_index;
    int *second_index;
    int **pair_index;
    // for each node, successor is the list of all descendants
    ReversibleIntList *successor;
    // for each node, predecessor is the list of all ancestors
    ReversibleIntList *predecessor;
    //@}

  public:
    /**@name Constructors*/
    //@{
    ConstraintTDAG(Solver*, VariableInt** x, const int n, 
		   VariableInt** interval1, 
		   VariableInt** interval2);
    virtual ~ConstraintTDAG();
    //@}

    /**@name Solving*/
    //@{
    inline int check( const int* ) const ;
    inline bool propagate();
    inline bool propagate(const int changedIdx, const int e); 
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual void print(std::ostream& o) const ;
    //@}
  };


  /**********************************************
   * Tree Constraint
   **********************************************/ 
  /*! \class ConstraintTree
    \brief  Tree Constraint.

  */
  class ConstraintTree : public Constraint {
  
  private:
    /**@name Parameters*/
    //@{  
    //ReversibleInt *root;
    ReversibleNum<int> *root;
    ReversibleSet *descendant;
    //ReversibleInt *leftChild;
    ReversibleNum<int> *leftChild;
    //ReversibleInt *rightChild;
    ReversibleNum<int> *rightChild;
    //@}

  public:
    /**@name Constructors*/
    //@{
    ConstraintTree(Solver*, VariableInt** x, const int n);
    virtual ~ConstraintTree();
    //@}

    /**@name Solving*/
    //@{
    inline int check( const int* ) const ;
    inline bool propagate();
    inline bool propagate(const int changedIdx, const int e); 
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual void print(std::ostream& o) const ;
    //@}
  };


  /**********************************************
   * Tree Constraint
   **********************************************/ 
  /*! \class ConstraintTreeExplicit
    \brief  Tree Constraint with reified descendant information.

  */
  class ConstraintTreeExplicit : public Constraint {
  
  private:
    /**@name Parameters*/
    //@{  
    int N;
    ReversibleNum<int> *root;
    VariableInt **descendant;
    //@}

  public:
    /**@name Constructors*/
    //@{
    ConstraintTreeExplicit(Solver*, VariableInt** x, const int n);
    virtual ~ConstraintTreeExplicit();
    //@}

    /**@name Solving*/
    //@{
    inline int check( const int* ) const ;
    inline bool propagate();
    inline bool propagate(const int changedIdx, const int e); 
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual void print(std::ostream& o) const ;
    //@}
  };


  /**********************************************
   * EdgeFinder Constraint
   **********************************************/ 
  /*! \class ConstraintEdgeFinder
    \brief  Constraint on a set of tasks sharing a unary resource.
  */
  class ConstraintEdgeFinder : public Constraint {

  public:
    /**@name Parameters*/
    //@{ 
    /*!
     * The durations of the tasks
     */
    int *durations;
    /*!
     * The tasks' indices sorted by decreasing max due date
     */
    int *sorted;
    /*!
     * The residual durations at each release date computed 
     * during the JPS
     */
    int **Pk;
    //@}

    /**@name Constructors*/
    //@{
    ConstraintEdgeFinder(Solver*, VariableInt** v, const int n, const int *d);
    virtual ~ConstraintEdgeFinder();
    //@}

    /**@name Utils*/
    //@{
    /*!
     * Compute the Jackson Preemptive Schedule necessary 
     * to the propagation method
     */
    bool jacksonPSchedule();
    //@}

    /**@name Solving*/
    //@{
    inline int check( const int* ) const ;
    inline bool propagate();
    inline bool propagate(const int changedIdx, const int e);
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual void print(std::ostream& o) const;
    //@}
  };


  /**********************************************
   * Or Constraint
   **********************************************/ 
  /*! \class ConstraintOr
    \brief  Binary Or Constraint (x0 or x1).
  */
  class ConstraintOr : public Constraint
  {

  public:
    /**@name Constructors*/
    //@{
    ConstraintOr(Solver*, VariableInt**);
    virtual ~ConstraintOr();
    //@}

    /**@name Solving*/
    //@{
    inline int check( const int* ) const ;
    inline bool propagate(const int changedIdx, const int e);
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual void print(std::ostream& o) const ;
    //@}
  };


  /* /\********************************************** */
  /*  * Nested Predicates */
  /*  **********************************************\/ */

  /**********************************************
   * BoolSum Equal Constraint
   **********************************************/
  //  lb <= x1 + ... + xn <= ub
  /// Constraint on the value of the sum of a set of variables.
  class ConstraintBoolSumEqual : public Constraint {
 
  public:
    /**@name Parameters*/
    //@{
    int total;
    //ReversibleInt min_;
    //ReversibleInt max_;
    ReversibleNum<int> min_;
    ReversibleNum<int> max_;
    ReversibleIntList unassigned;
    //@}

    /**@name Constructors*/
    //@{
    ConstraintBoolSumEqual(Solver*, VariableInt** v, 
			   const int n, const int t);
    virtual ~ConstraintBoolSumEqual();
    //@}

    /**@name Solving*/
    //@{
    inline int check( const int* ) const ;
    inline bool propagate();
    inline bool propagate(const int changedIdx, const int e);
    //@}
  
    /**@name Miscellaneous*/
    //@{  
    virtual void print(std::ostream& o) const ;
    //@}
  };


  /**********************************************
   * BoolSum Less Constraint
   **********************************************/
  //  lb <= x1 + ... + xn <= ub
  /// Constraint on the value of the sum of a set of variables.
  class ConstraintBoolSumLess : public Constraint {
 
  public:
    /**@name Parameters*/
    //@{
    int total;
    //ReversibleInt min_;
    ReversibleNum<int> min_;
    //@}

    /**@name Constructors*/
    //@{
    ConstraintBoolSumLess(Solver*, VariableInt** v, 
			  const int n, const int t);
    virtual ~ConstraintBoolSumLess();
    //@}

    /**@name Solving*/
    //@{
    inline int check( const int* ) const ;
    inline bool propagate();
    inline bool propagate(const int changedIdx, const int e);
    //@}
  
    /**@name Miscellaneous*/
    //@{  
    virtual void print(std::ostream& o) const ;
    //@}
  };


  /**********************************************
   * BoolSum More Constraint
   **********************************************/
  //  lb <= x1 + ... + xn <= ub
  /// Constraint on the value of the sum of a set of variables.
  class ConstraintBoolSumMore : public Constraint {
 
  public:
    /**@name Parameters*/
    //@{
    int total;
    //ReversibleInt max_;
    ReversibleNum<int> max_;
    //@}

    /**@name Constructors*/
    //@{
    ConstraintBoolSumMore(Solver*, VariableInt** v, 
			  const int n, const int t);
    virtual ~ConstraintBoolSumMore();
    //@}

    /**@name Solving*/
    //@{
    inline int check( const int* ) const ;
    inline bool propagate();
    inline bool propagate(const int changedIdx, const int e);
    //@}
  
    /**@name Miscellaneous*/
    //@{  
    virtual void print(std::ostream& o) const ;
    //@}
  };


  // /**********************************************
  //  * SimpleSum Predicate
  //  **********************************************/
  // /*! \class PredicateSimpleSum
  //  \brief  Sum of variables (x1 + ... + xn = y)
  // */
  // class PredicateSimpleSum : public Constraint {

  //  public:
  //   /**@name Constructors*/
  //   //@{
  //   PredicateSimpleSum(Solver*, VariableInt** v, const int n);
  //   virtual ~PredicateSimpleSum();
  //   //@}

  //   /**@name Solving*/
  //   //@{
  //   inline int check( const int* ) const ;
  //   inline bool propagate();
  //   inline bool propagate(const int changedIdx, const int e);
  //   //@}

  //   /**@name Miscellaneous*/
  //   //@{  
  //   virtual void print(std::ostream& o) const ;
  //   //@}
  // };


  // /**********************************************
  //  * SumEqual Constraint
  //  **********************************************/
  // /*! \class ConstraintSumEqual
  //  \brief  Constraint on a sum of variables (x1 + ... + xn = total)
  // */
  // class ConstraintSumEqual : public Constraint {

  //  public:
  //   /**@name Parameters*/
  //   //@{ 
  //   int total;
  //   //@}

  //   /**@name Constructors*/
  //   //@{
  //   ConstraintSumEqual(Solver*, 
  // 		     VariableInt** v,
  // 		     const int n,
  // 		     const int s);
  //   virtual ~ConstraintSumEqual();
  //   //@}

  //   /**@name Solving*/
  //   //@{
  //   inline int check( const int* ) const ;
  //   inline bool propagate();
  //   inline bool propagate(const int changedIdx, const int e);
  //   //@}

  //   /**@name Miscellaneous*/
  //   //@{  
  //   virtual void print(std::ostream& o) const ;
  //   //@}
  // };


  // /**********************************************
  //  * WeightedSumEqual Constraint
  //  **********************************************/
  // /*! \class ConstraintWeightedSumEqual
  //  \brief  Constraint on a sum of variables (a1 * x1 + ... + an * xn = total)
  // */
  // class ConstraintWeightedSumEqual : public Constraint {

  //  public:
  //   /**@name Parameters*/
  //   //@{ 
  //   int *weights;
  //   //@}

  //   /**@name Constructors*/
  //   //@{
  //   ConstraintWeightedSumEqual(Solver*, 
  // 			     VariableInt** v,
  // 			     const int n,
  // 			     const int *s);
  //   virtual ~ConstraintWeightedSumEqual();
  //   //@}

  //   /**@name Solving*/
  //   //@{
  //   inline int check( const int* ) const ;
  //   inline bool propagate();
  //   inline bool propagate(const int changedIdx, const int e);
  //   //@}

  //   /**@name Miscellaneous*/
  //   //@{  
  //   virtual void print(std::ostream& o) const ;
  //   //@}
  // };


  /**********************************************
   * WeightedSum Constraint
   **********************************************/
  /*! \class ConstraintWeightedSum
    \brief  Constraint on a sum of variables (a1 * x1 + ... + an * xn = total)
  */
  class ConstraintWeightedSum : public Constraint {

  public:
    /**@name Parameters*/
    //@{ 
    int *weights;
    int UB;
    int LB;
    int *up_bound;
    int *lo_bound;
    int wpos;
    int wneg;
    //@}

    /**@name Constructors*/
    //@{
    ConstraintWeightedSum(Solver*, 
			  VariableInt** v,
			  const int n,
			  const int u,
			  const int l,
			  const int *s,
			  const int wp,
			  const int wn);
    virtual ~ConstraintWeightedSum();
    //@}

    /**@name Solving*/
    //@{
    void checkBounds();
    inline int check( const int* ) const ;
    inline bool propagate();
    inline bool propagate(const int changedIdx, const int e);
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual void print(std::ostream& o) const ;
    void printFull(std::ostream& o) const ;
    //@}
  };



  /**********************************************
   * LinearCoef Constraint
   **********************************************/
  /*! \class ConstraintLinearCoef
    \brief  Constraint on a sum of variables (a1 * x1 + ... + an * xn = total)
  */
  class ConstraintLinearCoef : public Constraint {

  public:
    /**@name Parameters*/
    //@{ 
    double *weights;
    double UB;
    double LB;
    double *up_bound;
    double *lo_bound;
    int wpos;
    int wneg;
    //@}

    /**@name Constructors*/
    //@{
    ConstraintLinearCoef(Solver*, 
			  VariableInt** v,
			  const int n,
			  const double u,
			  const double l,
			  const double *s,
			  const int wp,
			  const int wn);
    virtual ~ConstraintLinearCoef();
    //@}

    /**@name Solving*/
    //@{
    void checkBounds();
    inline int check( const int* ) const ;
    inline bool propagate();
    inline bool propagate(const int changedIdx, const int e);
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual void print(std::ostream& o) const ;
    void printFull(std::ostream& o) const ;
    //@}
  };



  // /**********************************************
  //  * SumLess Constraint
  //  **********************************************/
  // /*! \class ConstraintSumLess
  //  \brief  Constraint on a sum of variables (x1 + ... + xn = total)
  // */
  // class ConstraintSumLess : public Constraint {

  //  public:
  //   /**@name Parameters*/
  //   //@{ 
  //   int total;
  //   //@}

  //   /**@name Constructors*/
  //   //@{
  //   ConstraintSumLess(Solver*, 
  // 		     VariableInt** v,
  // 		     const int n,
  // 		     const int s);
  //   virtual ~ConstraintSumLess();
  //   //@}

  //   /**@name Solving*/
  //   //@{
  //   inline int check( const int* ) const ;
  //   inline bool propagate();
  //   inline bool propagate(const int changedIdx, const int e);
  //   //@}

  //   /**@name Miscellaneous*/
  //   //@{  
  //   virtual void print(std::ostream& o) const ;
  //   //@}
  // };


  // /**********************************************
  //  * WeightedSum Predicate
  //  **********************************************/
  // /*! \class PredicateWeightedSum
  //  \brief  WeightedSum of variables with coefficients (a1 * x1 + ... + an * xn == y)
  // */
  // class PredicateWeightedSum : public Constraint 
  // {

  //  public:
  //   /**@name Parameters*/
  //   //@{ 
  //   int *weights;
  //   //@}

  //   /**@name Constructors*/
  //   //@{
  //   PredicateWeightedSum( Solver*, 
  // 			VariableInt** nvars,
  // 			const int n,
  // 			const int* os );
  //   virtual ~PredicateWeightedSum();
  //   //@}

  //   /**@name Solving*/
  //   //@{
  //   inline int check( const int* ) const ;
  //   inline bool propagate();
  //   inline bool propagate(const int changedIdx, const int e);
  //   //@}

  //   /**@name Miscellaneous*/
  //   //@{  
  //   virtual void print(std::ostream& o) const ;
  //   //@}
  // };


  /**********************************************
   * Equality with offset Predicate
   **********************************************/
  /*! \class PredicateOffset
    \brief  Equality with an offset (x + k = y)
  */
  class PredicateOffset : public Constraint//, public BitVar
  {

  private:
    BitSet xdom;

  public:
    /**@name Parameters*/
    //@{ 
    int offset;
    //@}

    /**@name Constructors*/
    //@{
    PredicateOffset(Solver*, VariableInt**, const int);
    virtual ~PredicateOffset();
    //@}

    /**@name Solving*/
    //@{
    inline int check( const int* ) const ;
    inline bool propagate(const int changedIdx, const int e);
    inline bool awakeOnRange ( VariableInt * );
    inline bool awakeOnDomain( VariableInt * );  
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual void print(std::ostream& o) const ;
    //@}
  };


  /**********************************************
   * Equality by a factor Predicate
   **********************************************/
  /*! \class PredicateFactor
    \brief  Equality by a factor (x * k = y)
  */
  class PredicateFactor : public Constraint//, public BitVar
  {

  public:
    /**@name Parameters*/
    //@{ 
    int factor;
    //@}

    /**@name Constructors*/
    //@{
    PredicateFactor(Solver*, VariableInt**, const int);
    virtual ~PredicateFactor();
    //@}

    /**@name Solving*/
    //@{
    inline int check( const int* ) const ;
    inline bool propagate(const int changedIdx, const int e);
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual void print(std::ostream& o) const ;
    //@}
  };


  /**********************************************
   * Negation Predicate
   **********************************************/
  /*! \class PredicateNegation
    \brief  Arithmetic negation (-x = y)
  */
  class PredicateNegation : public Constraint
  {

  private:
    /**@name Parameters*/
    //@{ 
    BitSet dom[2];
    //@}

  public:
    /**@name Constructors*/
    //@{
    PredicateNegation(Solver*, VariableInt**);
    virtual ~PredicateNegation();
    //@}

    /**@name Solving*/
    //@{
    inline int check( const int* ) const ;
    inline bool propagate(const int changedIdx, const int e);
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual void print(std::ostream& o) const ;
    //@}
  };


  /**********************************************
   * Absolute value Predicate
   **********************************************/
  /*! \class PredicateAbs
    \brief  Absolute value (|x| = y)
  */
  class PredicateAbs : public Constraint
  {

  public:
    /**@name Constructors*/
    //@{
    PredicateAbs(Solver*, VariableInt**);
    virtual ~PredicateAbs();
    //@}

    /**@name Solving*/
    //@{
    inline int check( const int* ) const ;
    inline bool propagate(const int changedIdx, const int e);
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual void print(std::ostream& o) const ;
    //@}
  };


  /**********************************************
   * Addition Predicate
   **********************************************/
  /*! \class PredicateAdd
    \brief  Binary addition predicate (x0 + x1 = y)
  */
  class PredicateAdd : public Constraint
  {

  private:
    /**@name Parameters*/
    //@{ 
    BitSet xdom;
    BitSet tdom;
    //@}

  public:
    /**@name Constructors*/
    //@{
    PredicateAdd(Solver*, VariableInt**);
    virtual ~PredicateAdd();
    //@}

    /**@name Solving*/
    //@{
    inline int check( const int* ) const ;
    inline bool propagate(const int changedIdx, const int e);
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual void print(std::ostream& o) const ;
    //@}
  };


  /**********************************************
   * Substraction Predicate
   **********************************************/
  /*! \class PredicateSub
    \brief  Binary Substraction predicate (x0 - x1 = y)
  */
  class PredicateSub : public Constraint
  {

  private:
    /**@name Parameters*/
    //@{ 
    BitSet xdom;
    BitSet tdom;
    //@}

  public:
    /**@name Constructors*/
    //@{
    PredicateSub(Solver*, VariableInt**);
    virtual ~PredicateSub();
    //@}

    /**@name Solving*/
    //@{
    inline int check( const int* ) const ;
    inline bool propagate(const int changedIdx, const int e);
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual void print(std::ostream& o) const ;
    //@}
  };


  /**********************************************
   * Multiplication Predicate
   **********************************************/
  /*! \class PredicateMul
    \brief  Binary Multiplication predicate (x0 * x1 = y)
  */
  class PredicateMul : public Constraint
  {

  private:
    /**@name Parameters*/
    //@{ 
    BitSet xdom;
    //@}

  public:
    /**@name Constructors*/
    //@{
    PredicateMul(Solver*, VariableInt**);
    virtual ~PredicateMul();
    //@}

    /**@name Solving*/
    //@{
    inline int check( const int* ) const ;
    inline bool propagate(const int changedIdx, const int e);

    inline int pruneZeros(const int changedIdx);
    inline bool pruneUnary(const int j);
    inline bool pruneBinary(const int i, const int j, const int v);
    inline bool pruneTernary(const int changedIdx);
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual void print(std::ostream& o) const ;
    //@}
  };


  /**********************************************
   * Division Predicate
   **********************************************/
  /*! \class PredicateDiv
    \brief  Binary Division predicate (x0 / x1 = y)
  */
  class PredicateDiv : public Constraint
  {

  private:
    /**@name Parameters*/
    //@{ 
    BitSet xdom;
    //@}

  public:
    /**@name Constructors*/
    //@{
    PredicateDiv(Solver*, VariableInt**);
    virtual ~PredicateDiv();
    //@}

    /**@name Solving*/
    //@{
    inline int check( const int* ) const ;
    inline bool propagate(const int changedIdx, const int e);
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual void print(std::ostream& o) const ;
    //@}
  };


  /**********************************************
   * Modulo Predicate
   **********************************************/
  /*! \class PredicateMod
    \brief  Binary Modulo predicate (x0 % x1 = y)
  */
  class PredicateMod : public Constraint
  {

  public:
    /**@name Constructors*/
    //@{
    PredicateMod(Solver*, VariableInt**);
    virtual ~PredicateMod();
    //@}

    /**@name Solving*/
    //@{
    inline int check( const int* ) const ;
    inline bool propagate(const int changedIdx, const int e);
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual void print(std::ostream& o) const ;
    //@}
  };


  /**********************************************
   * Modulo Predicate
   **********************************************/
  /*! \class PredicateMod
    \brief  Binary Modulo predicate (x0 % x1 = y)
  */
  class PredicateModConstant : public Constraint
  {

  public:

    int modulo;

    /**@name Constructors*/
    //@{
    PredicateModConstant(Solver*, VariableInt**, const int mod);
    virtual ~PredicateModConstant();
    //@}

    /**@name Solving*/
    //@{
    inline int check( const int* ) const ;
    inline bool propagate(const int changedIdx, const int e);
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual void print(std::ostream& o) const ;
    //@}
  };


  /**********************************************
   * Power Predicate
   **********************************************/
  /*! \class PredicatePow
    \brief  Binary Exponant predicate (x0 ^ x1 = y)
  */
  class PredicatePow : public Constraint
  {

  public:
    /**@name Constructors*/
    //@{
    PredicatePow(Solver*, VariableInt**);
    virtual ~PredicatePow();
    //@}

    /**@name Solving*/
    //@{
    inline int check( const int* ) const ;
    inline bool propagate(const int changedIdx, const int e);
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual void print(std::ostream& o) const ;
    //@}
  };


  /**********************************************
   * Minimum Predicate
   **********************************************/
  /*! \class PredicateMin
    \brief  Minimum of a set of variables (min(x0,...x1) = y)
  */
  class PredicateMin : public Constraint
  {

  public:
    /**@name Constructors*/
    //@{
    PredicateMin(Solver*, VariableInt**, const int);
    virtual ~PredicateMin();
    //@}

    /**@name Solving*/
    //@{
    inline int check( const int* ) const ;
    inline bool propagate(const int changedIdx, const int e);
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual void print(std::ostream& o) const ;
    //@}
  };


  /**********************************************
   * Maximum Predicate
   **********************************************/
  /*! \class PredicateMax
    \brief  Maximum of a set of variables (max(x0,...x1) = y)
  */
  class PredicateMax : public Constraint
  {

  public:
    /**@name Constructors*/
    //@{
    PredicateMax(Solver*, VariableInt**, int);
    virtual ~PredicateMax();
    //@}

    /**@name Solving*/
    //@{
    inline int check( const int* ) const ;
    inline bool propagate(const int changedIdx, const int e);
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual void print(std::ostream& o) const ;
    //@}
  };


  /**********************************************
   * Equality Predicate
   **********************************************/
  /*! \class PredicateEqual
    \brief  Truth value of a binary equality ((x0 = x1) <-> y)
  */
  class PredicateEqual : public Constraint
  {

  public:
    /**@name Parameters*/
    //@{ 
    int spin;
    //@}

    /**@name Constructors*/
    //@{
    PredicateEqual(Solver*, VariableInt**, int);
    virtual ~PredicateEqual();
    //@}

    /**@name Solving*/
    //@{
    inline int check( const int* ) const ;
    inline bool propagate(const int changedIdx, const int e);
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual void print(std::ostream& o) const ;
    //@}
  };


  /**********************************************
   * Equality with a constant Predicate
   **********************************************/
  /*! \class PredicateEqualConstant
    \brief  Truth value of a unary equality ((x0 = k) <-> y)
  */
  class PredicateEqualConstant : public Constraint
  {

  public:
    /**@name Parameters*/
    //@{ 
    int K;
    int spin;
    //@}

    /**@name Constructors*/
    //@{
    PredicateEqualConstant(Solver*, 
			   VariableInt**, 
			   const int,
			   const int);
    virtual ~PredicateEqualConstant();
    //@}

    /**@name Solving*/
    //@{
    inline int check( const int* ) const ;
    inline bool propagate(const int changedIdx, const int e);
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual void print(std::ostream& o) const ;
    //@}
  };



  /**********************************************
   * Membership to a constant set Predicate
   **********************************************/
  /*! \class PredicateMember
    \brief  Truth value of membership ((x0 \in k) <-> y)
  */
  class PredicateMember : public Constraint
  {

  public:
    /**@name Parameters*/
    //@{ 
    BitSet K;
    int spin;
    //@}

    /**@name Constructors*/
    //@{
    PredicateMember(Solver*, 
		    VariableInt**, 
		    const int*,
		    const int,
		    const int);
    virtual ~PredicateMember();
    //@}

    /**@name Solving*/
    //@{
    inline int check( const int* ) const ;
    inline bool propagate(const int changedIdx, const int e);
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual void print(std::ostream& o) const ;
    //@}
  };


  /**********************************************
   * Precedence Predicate
   **********************************************/
  /*! \class PredicatePrecedence
    \brief  Truth value of a binary precedence ((x0 + k <= x1) <-> y)
  */
  class PredicateLess : public Constraint
  {

  public:
    /**@name Parameters*/
    //@{ 
    int offset;
    //@}

    /**@name Constructors*/
    //@{
    PredicateLess(Solver*, VariableInt**, const int);
    virtual ~PredicateLess();
    //@}

    /**@name Solving*/
    //@{
    inline int check( const int* ) const ;
    inline bool propagate(const int changedIdx, const int e);
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual void print(std::ostream& o) const ;
    //@}
  };


  /**********************************************
   * Lower bound Predicate
   **********************************************/
  /*! \class PredicatelowerBound
    \brief  Truth value of a binary precedence ((x0 >= k) <-> y)
  */
  class PredicateLowerBound : public Constraint
  {

  public:
    /**@name Parameters*/
    //@{ 
    int bound;
    //@}

    /**@name Constructors*/
    //@{
    PredicateLowerBound(Solver*, VariableInt**, const int);
    virtual ~PredicateLowerBound();
    //@}

    /**@name Solving*/
    //@{
    inline int check( const int* ) const ;
    inline bool propagate(const int changedIdx, const int e);
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual void print(std::ostream& o) const ;
    //@}
  };


  /**********************************************
   * Upper bound Predicate
   **********************************************/
  /*! \class PredicateUpperBound
    \brief  Truth value of a binary precedence ((x0 <= k) <-> y)
  */
  class PredicateUpperBound : public Constraint
  {

  public:
    /**@name Parameters*/
    //@{ 
    int bound;
    //@}

    /**@name Constructors*/
    //@{
    PredicateUpperBound(Solver*, VariableInt**, const int);
    virtual ~PredicateUpperBound();
    //@}

    /**@name Solving*/
    //@{
    inline int check( const int* ) const ;
    inline bool propagate(const int changedIdx, const int e);
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual void print(std::ostream& o) const ;
    //@}
  };


  /**********************************************
   * Disjunctive Predicate
   **********************************************/
  /*! \class PredicateDisjunctive
    \brief  Ensures that two intervals do not overlap

    The last variable takes 0 if the first interval precedes
    the second, and 1 otherwise : (y -> (x0 + k0 <= x1) and !y -> (x1 + k1 <= x0))
  */
  class PredicateDisjunctive : public Constraint//, public BoolVar
  {

  public:
    /**@name Parameters*/
    //@{ 
    int duration[2];
    int& state;
    int& min_0;
    int& min_1;
    int& max_0;
    int& max_1;
    //@}

    /**@name Constructors*/
    //@{
    PredicateDisjunctive(Solver*, VariableInt**, const int*);
    virtual ~PredicateDisjunctive();
    //@}

    /**@name Solving*/
    //@{
    inline int check( const int* ) const ;
    inline bool propagate(const int changedIdx, const int e);
    //@}

    /**@name Miscellaneous*/
    //@{  

//     inline int intersection() {
//       int a = (min_0 < min_1 ? min_1 : min_0);
//       int b = (max_0+duration[0] > max_1+duration[1] ? max_1+duration[1] : max_0+duration[0]);
//       return( b > a ? b-a : 0 );
//     }

    inline int XltY() {
      int a = (max_0 - max_1 + duration[0]);
      int b = (min_0 + duration[0] - min_1);
      int c = (a > 0 ? a : 0) + (b > 0 ? b : 0);
      return c;
    }
    inline int YltX() {
      int a = (max_1 - max_0 + duration[1]);
      int b = (min_1 + duration[1] - min_0);
      int c = (a > 0 ? a : 0) + (b > 0 ? b : 0);
      return c;
    }
    inline int domsize() {
      return (max_1 + max_0 - min_1 - min_0 + 2);
    }
    virtual void print(std::ostream& o) const ;
    //@}
  };


  /**********************************************
   * Overlap Predicate | 0 = overlap, 1 = (t0 < t1), 2 = (t1 < t0)
   **********************************************/
  /*! \class PredicateOverlap
    \brief  Reifies the overlapping relation between two intervals

  */
  class PredicateOverlap : public Constraint//, public BoolVar
  {

  public:
    /**@name Parameters*/
    //@{ 
    int duration[2];
    unsigned int& state;
    int& min_0;
    int& min_1;
    int& max_0;
    int& max_1;
    //@}

    /**@name Constructors*/
    //@{
    PredicateOverlap(Solver*, VariableInt**, const int*);
    virtual ~PredicateOverlap();
    //@}

    /**@name Solving*/
    //@{
    inline int check( const int* ) const ;
    inline bool propagate(const int changedIdx, const int e);
    //@}

    /**@name Miscellaneous*/
    //@{  
    inline int XltY() {
      int a = (max_0 - max_1 + duration[0]);
      int b = (min_0 + duration[0] - min_1);
      int c = (a > 0 ? a : 0) + (b > 0 ? b : 0);
      return c;
    }
    inline int YltX() {
      int a = (max_1 - max_0 + duration[1]);
      int b = (min_1 + duration[1] - min_0);
      int c = (a > 0 ? a : 0) + (b > 0 ? b : 0);
      return c;
    }
    inline int domsize() {
      return (max_1 + max_0 - min_1 - min_0 + 2);
    }
    virtual void print(std::ostream& o) const ;
    //@}
  };


  /**********************************************
   * SafeDisjunct Predicate
   **********************************************/
  /*! \class PredicateSafeDisjunct
    \brief  Ensures that two intervals do not overlap

    The last variable takes 0 if the first interval precedes
    the second, and 1 otherwise : (y -> (x0 + k0 <= x1) and !y -> (x1 + k1 <= x0))
  */
  class PredicateSafeDisjunct : public Constraint//, public BoolVar
  {

  public:
    /**@name Parameters*/
    //@{ 
    int duration[2];
    //@}

    /**@name Constructors*/
    //@{
    PredicateSafeDisjunct(Solver*, VariableInt**, const int*);
    virtual ~PredicateSafeDisjunct();
    //@}

    /**@name Solving*/
    //@{
    inline int check( const int* ) const ;
    inline bool propagate(const int changedIdx, const int e);
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual void print(std::ostream& o) const ;
    //@}
  };


  /**********************************************
   * Hole Predicate
   **********************************************/
  /*! \class PredicateHole
    \brief  Ensures that two intervals do not overlap

    The last variable takes 0 if the first interval precedes
    the second, and 1 otherwise : (y -> (x0 + k0 <= x1) and !y -> (x1 + k1 <= x0))
  */
  class PredicateHole : public Constraint//, public BoolVar
  {

  public:
    /**@name Parameters*/
    //@{ 
    int bounds[2];
    int& state;
    int& min_x;
    int& max_x;
    //@}

    /**@name Constructors*/
    //@{
    PredicateHole(Solver*, VariableInt**, const int*);
    virtual ~PredicateHole();
    //@}

    /**@name Solving*/
    //@{
    inline int check( const int* ) const ;
    inline bool propagate(const int changedIdx, const int e);
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual void print(std::ostream& o) const ;
    //@}
  };


  /**********************************************
   * Not Predicate
   **********************************************/
  /*! \class PredicateNot
    \brief  Logic negation (!x0 = x1)
  */
  class PredicateNot : public Constraint
  {

  public:
    /**@name Constructors*/
    //@{
    PredicateNot(Solver*, VariableInt**);
    virtual ~PredicateNot();
    //@}

    /**@name Solving*/
    //@{
    inline int check( const int* ) const ;
    inline bool propagate();
    inline bool propagate(const int changedIdx, const int e);
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual void print(std::ostream& o) const ;
    //@}
  };


  /**********************************************
   * And Predicate
   **********************************************/
  /*! \class PredicateAnd
    \brief  Truth value of a binary conjunction (x0 & x1 <-> y)
  */
  class PredicateAnd : public Constraint
  {

  public:

    //static const int is_true = 0;
    //static const int is_false = 0;

    //ReversibleBool& domain_x;
    //ReversibleBool& domain_y;
    //ReversibleBool& domain_z;
    int& domain_x;
    int& domain_y;
    int& domain_z;

    /**@name Constructors*/
    //@{
    PredicateAnd(Solver*, VariableInt**);
    virtual ~PredicateAnd();
    //@}

    /**@name Solving*/
    //@{
    inline int check( const int* ) const ;
    inline bool propagate(const int changedIdx, const int e);
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual void print(std::ostream& o) const ;
    //@}
  };


  /**********************************************
   * Or Predicate
   **********************************************/
  /*! \class PredicateOr
    \brief  Truth value of a binary conjunction (x0 or x1 <-> y)
  */
  class PredicateOr : public Constraint
  {

  public:
    /**@name Constructors*/
    //@{
    PredicateOr(Solver*, VariableInt**);
    virtual ~PredicateOr();
    //@}

    /**@name Solving*/
    //@{
    inline int check( const int* ) const ;
    inline bool propagate(const int changedIdx, const int e);
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual void print(std::ostream& o) const ;
    //@}
  };


  /**********************************************
   * IfThenElse Predicate
   **********************************************/
  /*! \class PredicateIfThenElse
    \brief  Truth value of a if then else statement (x0 ? x1 : x2)
  */
  class PredicateIfThenElse : public Constraint
  {

  public:
    /**@name Constructors*/
    //@{
    PredicateIfThenElse(Solver*, VariableInt**);
    virtual ~PredicateIfThenElse();
    //@}

    /**@name Solving*/
    //@{
    inline int check( const int* ) const;
    inline bool propagate(const int changedIdx, const int e);
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual void print(std::ostream& o) const ;
    //@}
  };


  /**********************************************
   * Element Predicate
   **********************************************/
  /*! \class PredicateElement
    \brief  Element predicate

    Predicate on a vector x of variables and two variables y and z, ensuring that x[y] = z.
  */
  class PredicateElement : public Constraint
  {

  private:
    /**@name Parameters*/
    //@{ 
    int offset;
    BitSet aux_dom;
    //@}

  public:
    /**@name Constructors*/
    //@{
    PredicateElement(Solver*, VariableInt** v, const int n, const int o=0);
    virtual ~PredicateElement();
    //@}

    /**@name Solving*/
    //@{
    inline int check( const int* ) const ;
    inline bool propagate(const int changedIdx, const int e);
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual void print(std::ostream& o) const;
    //@}
  };


  /**********************************************
   * LDS Constraint
   **********************************************/
  /*! \class ConstraintLDS
    \brief  Constraint on a sum of variables (x1 + ... + xn = total)
  */
  class ConstraintLDS : public Constraint {

  public:
    int lb_threshold;
    int ub_threshold;
    ConstraintLDS(Solver*, VariableInt**, const int);
    virtual int getMaxThreshold() = 0;

  };


  /**********************************************
   * Hamming Constraint
   **********************************************/
  /*! \class ConstraintHamming
    \brief  Constraint on a sum of variables (x1 + ... + xn = total)
  */
  class BuildObject;
  class ConstraintHamming : public ConstraintLDS {

  public:
    /**@name Parameters*/
    //@{
    int                *ideal;
    //   ReversibleInt     mindiff;
    //   ReversibleInt     maxdiff;
    ReversibleNum<int>     mindiff;
    ReversibleNum<int>     maxdiff;
    ReversibleSet  maybeequal;
    ReversibleSet cannotequal;
    //@}

    /**@name Constructors*/
    //@{
    ConstraintHamming(Solver*);
    ConstraintHamming(Solver*, const int*);
    ConstraintHamming(Solver*, const int, BuildObject**, int*);
    void initialise(Solver*, const int*);
    virtual ~ConstraintHamming();
    //@}

    /**@name Solving*/
    //@{
    virtual int getMaxThreshold() { return arity; }
    inline int check( const int* ) const ;
    inline bool propagate(const int changedIdx, const int e);
    //@}

    /**@name Miscellaneous*/
    //@{
    virtual void print(std::ostream& o) const ;
    //@}
  };


  /**********************************************
   * Distance Constraint
   **********************************************/
  /*! \class ConstraintDistance
    \brief  Constraint on a sum of variables (x1 + ... + xn = total)
  */
  class ConstraintDistance : public ConstraintLDS {

  public:
    /**@name Parameters*/
    //@{
    int        **distance_to_val;
    int        **val_to_distance;
    int             maxthreshold;
    ReversibleNum<int>  *mindiff;
    ReversibleNum<int>  *maxdiff;
    //@}

    /**@name Constructors*/
    //@{
    ConstraintDistance(Solver*);
    virtual ~ConstraintDistance();
    //@}

    /**@name Solving*/
    //@{
    virtual int getMaxThreshold() { return maxthreshold; }
    inline int check( const int* ) const ;
    inline bool propagate(const int changedIdx, const int e);
    //@}

    /**@name Miscellaneous*/
    //@{
    virtual void print(std::ostream& o) const ;
    //@}
  };



#ifdef _CPLEX

  /**********************************************
   * CPlex Constraint
   **********************************************/
  class CPlexConstraint : public Constraint {

  protected:
    /**@name Parameters*/
    //@{ 
    int nvals;
    int nRows;
    int nColumns;
    int lb;
    int ub;
    int type;

    int *values2index;
    int *index2values;

    int **asgn2index;
    int *index2valasgn;
    int *index2varasgn;

    int *lbs;

    double epsilon;

    void updateModel();


    /*************************
     * CPlex structs
     *************************/
    int status;

    char     *probname;  
    int      objsen;
    double   *obj;
    double   *rhs;
    char     *sense;
    int      *matbeg;
    int      *matcnt;
    int      *matind;
    double   *matval;
    double   *zlb;
    double   *zub;
    int      *qmatbeg;
    int      *qmatcnt;
    int      *qmatind;
    double   *qmatval;
  

    int      solstat;
    double   objval;
    double   *x;//[NUMCOLS];
    double   *pi;//[NUMROWS];
    double   *slack;//[NUMROWS];
    double   *dj;//[NUMCOLS];
  
    CPXENVptr     env;
    CPXLPptr      lp;
    /*************************
     *                       *
     *************************/
    //@}

  public:
    /**@name Constructors*/
    //@{
    CPlexConstraint(Solver*, VariableInt** v, const int n);
    virtual ~CPlexConstraint();
    //@}

    /**@name Solving*/
    //@{
    inline int check( const int* ) const ;
    virtual bool propagate(const int changedIdx, const int e);
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual void print(std::ostream& o) const ;
    //@}
  };


  /**********************************************
   * CPlexMinDiff Predicate
   **********************************************/
  /// Counts the maximum number of differences in a set of vars
  class PredicateCPlexMinDiff : public CPlexConstraint
  {  
  public:

    /**@name Parameters*/
    //@{ 
    int K;
    int *occurrence;
    int *orderedvals;
    VariableInt **vars;
    int  filtering;
    //@}

    /**@name Constructors*/
    //@{
    PredicateCPlexMinDiff(Solver*, const int f, VariableInt** x, const int n);
    virtual ~PredicateCPlexMinDiff();
    //@}

    /**@name Solving*/
    //@{
    int approximation(int&);

    inline bool propagate(const int changedIdx, const int e);
    //@}
  
    /**@name Miscellaneous*/
    //@{  
    virtual void print(std::ostream& o) const ;
    //@}
  };

#endif


#ifdef _COINLP

  /**********************************************
   * LP Constraint
   **********************************************/
  class LpConstraint : public Constraint {

  public:

    /**@name Parameters*/
    //@{ 
    static const int SIMPLEX_PRIMAL     = 0;
    static const int SIMPLEX_DUAL       = 1;
    static const int SIMPLEX_BARRIER    = 2;
    static const int SIMPLEX_PRIMAL_IT  = 3;
    static const int SIMPLEX_DUAL_IT    = 4;
    static const int SIMPLEX_BARRIER_IT = 5;
    static const int INTERIOR           = 6;

    static const int LPO = 0;
    static const int LPA = 1;
    static const int LPS = 2;
    static const int APX = 3;
    static const int UPB = 4;


  protected:
  
    ClpSimplex  smodel;
    ClpInterior imodel;

    int nvals;
    int nRows;
    int nColumns;
    int lb;
    int ub;
    int type;

    int *values2index;
    int *index2values;

    int **asgn2index;
    int *index2valasgn;
    int *index2varasgn;

    int *lbs;

    double epsilon;

    void updateModel();
    //@}

  public:
    /**@name Constructors*/
    //@{
    LpConstraint(Solver*);
    LpConstraint(Solver*, VariableInt** v, const int n, const int t=6);
    void init(const int t);
    virtual ~LpConstraint();
    //@}

    /**@name Solving*/
    //@{
    inline int check( const int* ) const ;
    virtual bool propagate(const int changedIdx, const int e);
    //@}

    /**@name Miscellaneous*/
    //@{    
    virtual void print(std::ostream& o) const ;
    //@}
  };


  /**********************************************
   * MaxDiff Predicate
   **********************************************/
  /// Counts the maximum number of differences in a set of vars
  class PredicateMaxDiff : public LpConstraint
  {  
  
  public:

    /**@name Parameters*/
    //@{ 
    BitSet freevals;
    int *groundvals;
    int *orderedvals;
    int  filtering;
    //@}

    /**@name Constructors*/
    //@{
    PredicateMaxDiff(Solver*, VariableInt** v, const int n, const int t=6);
    PredicateMaxDiff(Solver*, const int f, VariableInt** v, const int n, const int t=6);
    void initStruct(const int t);
    virtual ~PredicateMaxDiff();
    //@}

    /**@name Solving*/
    //@{
    int upperbound();
    int approximation(int&);

    inline int check(int *sol) const;
    inline bool propagate(const int changedIdx, const int e);
    //@}
  
    /**@name Miscellaneous*/
    //@{  
    virtual void print(std::ostream& o) const ;
    //@}
  };


  /**********************************************
   * MinDiff Predicate
   **********************************************/
  /// Counts the maximum number of differences in a set of vars
  class PredicateMinDiff : public LpConstraint
  {  
  public:

    /**@name Parameters*/
    //@{ 
    int K;
    int *occurrence;
    int *orderedvals;
    VariableInt **vars;
    int  filtering;

    ClpCholeskyBase * cholesky;
    //@}

    /**@name Constructors*/
    //@{
    PredicateMinDiff(Solver*, VariableInt** x, const int n, const int t=6);
    PredicateMinDiff(Solver*, const int f, VariableInt** x, const int n, const int t=6);
    void initStruct(const int t);
    virtual ~PredicateMinDiff();
    //@}

    /**@name Solving*/
    //@{
    int approximation(int&);

    inline bool propagate(const int changedIdx, const int e);
    //@}

    /**@name Miscellaneous*/
    //@{  
    virtual void print(std::ostream& o) const ;
    //@}
  };

#endif // _COINLP


  /**********************************************
   * Decomposition of Slide Constraint
   **********************************************/ 
  /*! \class ConstraintSlideDecomp
    \brief  implementation of slide
  */
  class ConstraintSlideDecomp : public Constraint {

  public:
    /**@name Parameters*/
    //@{ 
    int nbVals; 
    int sequenceWidth;
    int* conflicts;
    int* int2ind;
    int* ind2int;
    int* exposant;
    //@}

    /**@name Utilities*/
    //@{ 
    inline void tuple2int(int& x, int* tuple, int n) {
      x = 0;
      while( n-- )
	x += (tuple[n] * exposant[n]);  
    }  
    inline int ktuple( const int left, int const right ) const {
      return ind2int[left]%nbVals + ind2int[right]*nbVals;
    }
    //@{ 

    /**@name Constructors*/
    //@{
    ConstraintSlideDecomp(Solver*, VariableInt**, 
			  const int, const int, const int, 
			  const int, const int,  
			  int*, int*, int*, int*);
    ~ConstraintSlideDecomp();
    //@{ 

    /**@name Solving*/
    //@{
    inline int check( const int* ) const {return 0;}  
    inline bool propagate(const int changedIdx, const int e) {return true;}
    //@}
  };


  /**********************************************
   * Slide Constraint
   **********************************************/ 
  /*! \class ConstraintSlide
    \brief  implementation of tuple sliding
  */
  class ConstraintSlide : public Constraint {

  protected:

    /**@name Parameters*/
    //@{ 
    int nbVals; 
    int sequenceWidth;
    int* conflicts;
    int* int2ind;
    int* ind2int;
    int* exposant;
    //@}

    /**@name Utilities*/
    //@{ 
    inline int appendLeft(const int t, const int e) const {
      return ind2int[t]*nbVals+e;
    }  
    inline int appendRight(const int t, const int e) const {
      return ind2int[t]+e*exposant[sequenceWidth-1];
    }
    //@}

  public:

    /**@name Constructors*/
    //@{
    ConstraintSlide(Solver*, VariableInt**, 
		    ConstraintSlideDecomp*);
    ~ConstraintSlide();
    //@}

    /**@name Solving*/
    //@{
    inline int check( const int* ) const ;  
    inline bool propagate(const int changedIdx, const int e);
    //@}

    /**@name Miscellaneous*/
    //@{
    virtual void print(std::ostream& o) const;
    //@}
  };


  /**********************************************
   * Tuple Constraint
   **********************************************/ 
  /*! \class PredicateTuple
    \brief  implementation of tuple
  */
  class PredicateTuple : public Constraint {

  protected:

    /**@name Parameters*/
    //@{ 
    int           nbVals;
    BitSet         *vals;
    int         *int2ind;
    int         *ind2int;
    int        *exposant;
    //@}

  public:

    /**@name Constructors*/
    //@{
    PredicateTuple(Solver*, VariableInt**,
		   const int, const int,
		   int*, int*, int*);
    ~PredicateTuple();
    //@}

    /**@name Solving*/
    //@{
    inline int check( const int* ) const ;  
    inline bool propagate(const int changedIdx, const int e);
    //@}

    /**@name Miscellaneous*/
    //@{
    virtual void print(std::ostream& o) const;
    //@}
  };






  /**********************************************
   * Disjunctive Predicate bis
   **********************************************/
  /*! \class PredicateDisjunctive
    \brief  Ensures that two intervals do not overlap

    The last variable takes 0 if the first interval precedes
    the second, and 1 otherwise : (y -> (x0 + k0 <= x1) and !y -> (x1 + k1 <= x0))
  */
  class ConstraintClauseBase;

  class ConstraintSatLess : public Constraint//, public BoolVar
  {

  public:
    /**@name Parameters*/
    //@{   
    int duration;
    int& min_0;
    int& min_1;
    int& max_0;
    int& max_1;

    int atom_id;
    int task_id[2];

    // the reason of the last pruning on max(X0)
    ReversibleNum<int> *reason_for_ub_decrease;
    // the reason of the last pruning on min(X1)
    ReversibleNum<int> *reason_for_lb_increase;

    // 
    int *responsible_ub;
    int *responsible_lb;
    int *responsible_disjunct;
    ConstraintClauseBase *sat;
    //@}

    /**@name Constructors*/
    //@{
    ConstraintSatLess(Solver*, VariableInt**, const int, ConstraintClauseBase*);
    virtual ~ConstraintSatLess();
    //@}

    /**@name Solving*/
    //@{
    inline int check( const int* ) const ;
    inline bool propagate(const int changedIdx, const int e);
    //@}
    virtual void print(std::ostream& o) const ;
    //@}
  };

  class PredicateDisjunct : public Constraint//, public BoolVar
  {

  public:
    /**@name Parameters*/
    //@{   
    int duration[2];
    int& state;
    int& min_0;
    int& min_1;
    int& max_0;
    int& max_1;

    int atom_id;
    int task_id[2];
    // the two reasons (lb and ub) why this decision was taken
    int *reason_for_choice[2];

    // the reason of the last pruning on max(X0) and max(X1)
    //int *reason_for_ub_decrease[2];
    ReversibleNum<int> *reason_for_ub_decrease[2];
    // the reason of the last pruning on min(X0) and min(X1)
    //int *reason_for_lb_increase[2];
    ReversibleNum<int> *reason_for_lb_increase[2];

    // 
    //int *failed_task[2];
    Clause **reason;
    int *responsible_ub;
    int *responsible_lb;
    int *responsible_disjunct;
    ConstraintClauseBase *sat;
    //@}

    /**@name Constructors*/
    //@{
    PredicateDisjunct(Solver*, VariableInt**, const int*, ConstraintClauseBase*);
    virtual ~PredicateDisjunct();
    //@}

    /**@name Solving*/
    //@{
    inline int check( const int* ) const ;
    inline bool propagate(const int changedIdx, const int e);
    //@}

    /**@name Miscellaneous*/
    //@{  
    inline int XltY() {
      int a = (max_0 - max_1 + duration[0]);
      int b = (min_0 + duration[0] - min_1);
      int c = (a > 0 ? a : 0) + (b > 0 ? b : 0);
      return c;
    }
    inline int YltX() {
      int a = (max_1 - max_0 + duration[1]);
      int b = (min_1 + duration[1] - min_0);
      int c = (a > 0 ? a : 0) + (b > 0 ? b : 0);
      return c;
    }
    inline int domsize() {
      return (max_1 + max_0 - min_1 - min_0 + 2);
    }
    virtual void print(std::ostream& o) const ;
    //@}
  };





  /**********************************************
   * Clause Base
   **********************************************/ 
  /*! \class ConstraintClauseBase
    \brief  implementation of a cnf
  */
  class SatSolver;
  class ConstraintClauseBase : public Constraint, public SatSolver {

  public:

    static const int CSP=0;
    static const int SAT=1;
    static const int DTP=2;
  
    /**@name Parameters*/
    //@{ 
    int mode;
    bool is_choice;

    BitSet visitedUbs;
    BitSet visitedLbs;
    BitSet visitedDisjuncts;
    List   lbToExplain;
    List   ubToExplain;
    List   disjunctToExplain;


    Solver *solver;
    VariableInt **X;
    Clause *conflict;
    bool use_backjump;
    bool use_learning;



    int failure[3];
    Vector<PredicateDisjunct*> disjuncts;
    Vector<ConstraintSatLess*> precedences;
    int nTasks;
    ReversibleNum<int> *reason_lb;
    ReversibleNum<int> *reason_ub;
    int *reason_lit[2];
    int *disjunct[2];
  
    VariableInt **the_task;


    Vector< VariableInt** > cl_scope;
    Vector< int* > cl_lits;
    Vector< int > cl_size;

    void addDelayedClause(VariableInt **x, const int n, int *l);
//  {
//       VariableInt** scp = new VariableInt*[n];
//       int* lits = new int[n];
//       for(int i=0; i<n; ++i) {
// 	scp[i] = x[i];
// 	lits[i] = l[i]; 
//       }
//       cl_scope.push(scp);
//       cl_lits.push(lits);
//       cl_size.push(n);
//     }
    //@}

    /**@name Constructors*/
    //@{
    ConstraintClauseBase(Solver*);
    void initialise();
    virtual ~ConstraintClauseBase();
    void addDisjunct( PredicateDisjunct *d )
    {    
      disjuncts.push( d );
    }
    void addPrecedence( ConstraintSatLess *p )
    {    
      precedences.push( p );
    }
    //@}

    /**@name Solving*/
    //@{
    inline int check( const int* ) const ;  
    inline bool propagate();
    inline bool propagate(const int changedIdx, const int e);
    void resolve( Clause *c, Atom a );
    void analyzeConflict();
    void exploreImplicationGraph();
    void extractNogood();
    void explainBounds();
    void update();
    inline void save() { store->push( this ); } 
    void restore();

    // explore the explanations of the pruning on task's upper bound
    void explain_upperbound( const int task );
    // explore the explanations of the pruning on task's lower bound
    void explain_lowerbound( const int task );
    // explore the explanations of the polarity of this disjunct
    void explain_disjunct( const int d );
    int explain_disjunct2( const int d, int& count );

    inline bool isDecision(const int a)
    {
      return( !reason[a] && !reason_lit[0][a] && !reason_lit[1][a] );
    }
    //@}

    /**@name Miscellaneous*/
    //@{
    void checkConsistency();
    virtual void print(std::ostream& o) const;
    //@}
  };









  // /**********************************************
  //  * Disjunctive Predicate bis
  //  **********************************************/
  // /*! \class PredicateDisjunctive
  //  \brief  Ensures that two intervals do not overlap

  //  The last variable takes 0 if the first interval precedes
  //  the second, and 1 otherwise : (y -> (x0 + k0 <= x1) and !y -> (x1 + k1 <= x0))
  // */
  // class ConstraintClauseBase;
  // class PredicateDisjunct : public Constraint//, public BoolVar
  // {

  //  public:
  //   /**@name Parameters*/
  //   //@{   
  //   int duration[2];
  //   int& state;
  //   int& min_0;
  //   int& min_1;
  //   int& max_0;
  //   int& max_1;

  //   int atom_id;
  //   int task_id[2];
  //   // the two reasons (lb and ub) why this decision was taken
  //   int *reason_for_choice[4];

  //   // the reason of the last pruning on max(X0) and max(X1)
  //   //int *reason_for_ub_decrease[2];
  //   ReversibleNum<int> *reason_for_ub_decrease[2];
  //   // the reason of the last pruning on min(X0) and min(X1)
  //   //int *reason_for_lb_increase[2];
  //   ReversibleNum<int> *reason_for_lb_increase[2];

  //   // 
  //   //int *failed_task[2];
  //   Clause **reason;
  //   int *responsible_ub;
  //   int *responsible_lb;
  //   int *responsible_disjunct;
  //   ConstraintClauseBase *sat;
  //   //@}

  //   /**@name Constructors*/
  //   //@{
  //   PredicateDisjunct(Solver*, VariableInt**, const int*, ConstraintClauseBase*);
  //   virtual ~PredicateDisjunct();
  //   //@}

  //   /**@name Solving*/
  //   //@{
  //   inline int check( const int* ) const ;
  //   inline bool propagate(const int changedIdx, const int e);
  //   //@}

  //   /**@name Miscellaneous*/
  //   //@{  
  //   inline int XltY() {
  //     int a = (max_0 - max_1 + duration[0]);
  //     int b = (min_0 + duration[0] - min_1);
  //     int c = (a > 0 ? a : 0) + (b > 0 ? b : 0);
  //     return c;
  //   }
  //   inline int YltX() {
  //     int a = (max_1 - max_0 + duration[1]);
  //     int b = (min_1 + duration[1] - min_0);
  //     int c = (a > 0 ? a : 0) + (b > 0 ? b : 0);
  //     return c;
  //   }
  //   inline int domsize() {
  //     return (max_1 + max_0 - min_1 - min_0 + 2);
  //   }
  //   virtual void print(std::ostream& o) const ;
  //   //@}
  // };




  // /**********************************************
  //  * Clause Base
  //  **********************************************/ 
  // /*! \class ConstraintClauseBase
  //  \brief  implementation of a cnf
  // */
  // class SatSolver;
  // class ConstraintClauseBase : public Constraint, public SatSolver {

  // public:

  //   static const int CSP=0;
  //   static const int SAT=1;
  //   static const int DTP=2;
  
  //   /**@name Parameters*/
  //   //@{ 
  //   int mode;
  //   bool is_choice;

  //   BitSet visitedUbs;
  //   BitSet visitedLbs;
  //   BitSet visitedDisjuncts;
  //   List   lbToExplain;
  //   List   ubToExplain;
  //   List   disjunctToExplain;


  //   Solver *solver;
  //   VariableInt **X;
  //   Clause *conflict;
  //   bool use_backjump;
  //   bool use_learning;



  //   int failure[3];
  //   Vector<PredicateDisjunct*> disjuncts;
  //   int nTasks;
  //   ReversibleNum<int> *reason_lb;
  //   ReversibleNum<int> *reason_ub;
  //   int *reason_lit[4];
  //   int *disjunct[2];
  
  //   VariableInt **the_task;
  //   //@}

  //   /**@name Constructors*/
  //   //@{
  //   ConstraintClauseBase(Solver*);
  //   void initialise();
  //   virtual ~ConstraintClauseBase();
  //   void addDisjunct( PredicateDisjunct *d )
  //   {    
  //     disjuncts.push( d );
  //   }
  //   //@}

  //   /**@name Solving*/
  //   //@{
  //   inline int check( const int* ) const ;  
  //   inline bool propagate();
  //   inline bool propagate(const int changedIdx, const int e);
  //   void resolve( Clause *c, Atom a );
  //   void analyzeConflict();
  //   void exploreImplicationGraph();
  //   //void update();
  //   inline void save() { store->push( this ); } 
  //   void restore();

  //   // explore the explanations of the pruning on task's upper bound
  //   void explain_upperbound( const int task );
  //   // explore the explanations of the pruning on task's lower bound
  //   void explain_lowerbound( const int task );
  //   // explore the explanations of the polarity of this disjunct
  //   void explain_disjunct( const int d );

  //   inline bool isDecision(const int a)
  //   {
  //     return( lvl[a] && !reason[a] && !reason_lit[0][a] && !reason_lit[1][a] );
  //   }
  //   //@}

  //   /**@name Miscellaneous*/
  //   //@{
  //   void checkConsistency();
  //   virtual void print(std::ostream& o) const;
  //   //@}
  // };

};

#endif // _CONSTRAINT_H
