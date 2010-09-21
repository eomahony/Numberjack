
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



/*! \file csp.h
    \brief Header for the data structures representing a CSP.
*/

/**********************************************
 * The core: CSP, Constraint, Variable,
 **********************************************/

#ifndef _CSP_H
#define _CSP_H

#include <iostream>
#include <iomanip>
#include <vector>
#include <mistral_set.h>
#include <mistral_glo.h>

#include <map>


namespace Mistral {

  class ValSelector;
  class VariableInt;
  class Constraint;
  class Solver;
  class CSP;

  class UnaryConstraint;
  class ObjectiveFunction;



  /********************************************
   * Reversible Objects
   ********************************************/
  /*! \class ReversibleObj
    \brief Backtrackable data structures.

    All structures that need to be restored to previous
    state during search implement the method
    restore(). 
  */
  class ReversibleObj {
  public:
    /*!@name Parameters*/
    //@{
    Vector<ReversibleObj*> *store;
    int *level;
    //@}

    /*!@name Constructors*/
    //@{
    ReversibleObj() {}
    virtual ~ReversibleObj() {}
    //@}

    /*!@name Backtrack method*/
    //@{
    virtual void restore() = 0; 
    //virtual void save() {} 
    //@}

    /*!@name Miscellaneous*/
    //@{
    /// Print on out stream
    virtual void print(std::ostream& o) const {}
    //@}
  };



  /********************************************
   * Domain Iterator
   ********************************************/
  /*! \class DomainIterator
    \brief Values Iteration Accessors.
  */
  class DomainIterator {

  public:

    /*!@name Constructors*/
    //@{
    DomainIterator() {}
    virtual ~DomainIterator() {}
    //@}

    /*!@name Accessors*/
    //@{
    virtual bool next()=0;
    virtual operator const int() const=0;
    //@}
  };

  template< class T >
  class MistralGacList;

  /**********************************************
   * VariableInt
   **********************************************/
  /*! \class VariableInt
    \brief Abstract Representation of Integer Variables.

  */
  class VariableInt {
  
  public:

    /********************************************
     * Constraint graph attributes.
     ********************************************/
    /*!@name Constraint Graph Parameters*/
    //@{
    /// The list of constraints that should be triggered when this variable is assigned
    MistralList<Constraint*> valueTrigger;
    /// The list of constraints that should be triggered when the bounds of this variable change
    MistralList<Constraint*> rangeTrigger;
    /// The list of constraints that should be triggered when the domain of this variable changes
    MistralList<Constraint*> domainTrigger;
    /// The list of constraints whose watchers include this variable
    MistralList<Constraint*> watchedConstraints;
    //@}

    /*!@name Constraint Graph Accessors*/
    //@{
    /// The list of constraints triggered on value assignment
    virtual MistralList<Constraint*>* triggerOnValue() {return &valueTrigger;}
    MistralNode<Constraint*>* constraintsOnValue() {return &(valueTrigger.head);}

    /// The list of constraints triggered on bound(s) reduction
    virtual MistralList<Constraint*>* triggerOnRange() {return &rangeTrigger;}
    MistralNode<Constraint*>* constraintsOnRange() {return &(rangeTrigger.head);}
  
    /// The list of constraints triggered on value(s) removal
    virtual MistralList<Constraint*>* triggerOnDomain() {return &domainTrigger;}
    MistralNode<Constraint*>* constraintsOnDomain() {return &(domainTrigger.head);}
    //@}



    /********************************************
     * Abstract variable attributes.
     ********************************************/
    /*!@name Static variable types*/
    //@{
    /// Constant type
    static const int CONST   = 0;
    /// Bitset type
    static const int BIT     = 1;
    /// Boolean type
    static const int BOOL    = 2;
    /// Range type
    static const int RANGE   = 3;
    /// Integer List type
    static const int LIST    = 4;
    /// Integer List type
    static const int INTLIST = 5;
    /// View
    static const int VIRTUAL = 6;
    virtual int getType() const = 0;
    virtual int& getIntDomain() = 0;
    //@}

    /*!@name Static value for non-value*/
    //@{
    //static const int NOVAL = INT_MAX;
    //@}

    /*!@name Abstract variable Parameters and Methods*/
    //@{
    /// Direct pointer to the solver
    Solver *solver;
    /// Pointer to the variable Queue
    MistralGacList<VariableInt*> *ACQueue;
    /// Value (branching) selection;
    ValSelector *branch;
    /// The index of this Variable.
    int id;
    /// The agreggated constraint weight
    long unsigned int weight;
    /// Number active Constraints on this Variable
    int degree;


    /// Pointer to its cell in the array 'sequence' of the solver
    /// Used when linking/unlinking the variable.
    VariableInt** seqIdx;
    /// Pointer to either the array 'sequence', or 'auxilliary'
    /// Used when linking/unlinking the variable.
    VariableInt*** seqBeg;
    /// Whether or not the variable is in the current constraint network
    inline bool isLinked() const { return (seqIdx >= *seqBeg); } 
    /// unLink when it becomes assigned and link when backtracking
    /// Acts as if the Variable was removed from the Constraint network
    void unLink();
    /// Restore the Variable in the Constraint network
    void link();
    /// Return itself
    virtual VariableInt* getVar() { return this; }
    //@}


    /*!@name Constructors*/
    //@{
    VariableInt();
    VariableInt(Solver*);
    virtual ~VariableInt();
    //@}
 

    /*!@name Accessors and Iterators*/
    //@{
    /// Returns the assigned value if it exists NOVAL otherwise
    virtual int value() const = 0;
    /// the same method, in testing
    int getValue() const;
    int getMin() const;
    int getMax() const;

    virtual int getNext(const int) const=0;

    /// DomainIterator
    virtual DomainIterator *begin()=0;

    /// Returns a random value in the domain
    virtual int random() const = 0;

    /// Returns the first (not necessarily min) value in the domain
    virtual int first() const = 0;

    /// Returns the minimum value in the domain
    virtual int min() const = 0;

    /// Returns the maximum value in the domain
    virtual int max() const = 0;

    /// Returns the minimum absolute value in [0..max] \\inter domain
    virtual int minPosAbs() const = 0;

    /// Returns the minimum absolute value in [min..0] \\inter domain
    virtual int minNegAbs() const = 0;

    /// Returns the size of the domain
    virtual int domsize() const = 0;

    /// Returns the minimum value that could belong to the domain
    virtual int minCapacity() const = 0;

    /// Returns 1 + the maximum value that could belong to the domain
    virtual int maxCapacity() const = 0;

    /// Value increment 
    virtual bool setNext(int& v) const = 0;
    //@}


    /*!@name Query methods*/
    //@{
    /// Whether or not the Variable is currently an interval
    virtual bool isRange() const = 0;

    /// Whether or not the Variable is bound to a ground value
    virtual bool isGround() const = 0;

    /// Whether or not the Variable is bound to a given ground value
    virtual bool equal(const int v) const = 0;

    /*!
      Whether "v" is currently contained in the domain
    */
    virtual bool contain(const int v) const = 0;

    // set operation functions
    /*!
      Whether the domain has a nonempty intersection
      with the interval [l..u]
    */
    virtual bool intersect(const int lo, const int up) const = 0;

    /*!
      Whether the domain is included in
      the interval [l..u]
    */
    virtual bool included(const int lo, const int up) const = 0;

    /*!
      Whether the domain has a nonempty intersection
      with the Variable x
    */
    virtual bool intersect(const VariableInt* x) const = 0;

    /*!
      Whether the domain is included
      in the Variable x
    */
    virtual bool included(const VariableInt* x) const = 0;

    /*!
      Whether the Variable x is included
      in the domain
    */
    virtual bool include(const VariableInt* x) const = 0;

    /*!
      Whether the domain has a nonempty intersection
      with the set s
    */
    virtual bool intersect(const BitSet& s) const = 0;

    /*!
      Whether the domain is included
      in the set s
    */
    virtual bool included(const BitSet& s) const = 0;

    /*!
      Whether the set s is included
      in the domain
    */
    virtual bool include(const BitSet& s) const = 0;

    /*!
      Intersect its domain with a set s
    */
    virtual void intersectTo( BitSet& s ) const = 0;

    /*!
      Do the union of its domain with a set s
    */
    virtual void unionTo( BitSet& s ) const = 0;

    /*!
      Copy its domain into a set s
    */
    virtual void copyTo( BitSet& s ) const = 0;
    //@}


    /*!@name Domain handling methods*/
    //@{
    /// Remove the value "v"
    virtual bool remove(const int v) = 0;
    /// Remove all values but "v"
    virtual bool setDomain(const int v) = 0;
    /// Remove all values strictly lower than l
    virtual bool setMin(const int l) = 0;
    /// Remove all values strictly greater than u
    virtual bool setMax(const int u) = 0;
    /// Remove all values that do not appear in the array "a" of length "l"
    virtual bool setDomain(const int* a, const int l) = 0;
    /// Remove all values that do not appear in the set "s"
    virtual bool setDomain(const BitSet& s) = 0;
    /// Remove all values that do not appear in the current domain of the Variable "x"
    virtual bool setDomain( VariableInt* x) const = 0;
    /// Remove all values that belong to the set "s"
    virtual bool removeSet(const BitSet& s) = 0;
    /// Remove all values in the interval [l..u]
    virtual bool removeRange(const int l, const int u) = 0;
    //@}

    /*!@name Miscellanous*/
    //@{ 
    virtual void print(std::ostream&) const = 0;
    virtual void printshort(std::ostream&) const = 0;
    virtual void printDomain(std::ostream&) const = 0;

    /// Dirty hack for efficiency issue on binary extensional constraints
    bool isWord;
    /// Dirty hack for efficiency issue on binary extensional constraints
    virtual bool revise( BitSet*, int*, VariableInt* ) { return true; }
    /// Dirty hack for efficiency issue on binary extensional constraints
    virtual bool revise( Constraint*, const int ) { return true; }
    /// Dirty hack for efficiency issue on binary extensional constraints
    virtual bool fastContain(const int v) const { return contain(v); }
    /// Dirty hack for efficiency issue on binary extensional constraints
    virtual bool intersect(const BitSet& s, int& r) const { return true; };
    /// Dirty hack for efficiency issue on binary extensional constraints
    virtual bool wordIntersect(const BitSet& s) const { return true; };
    /// Change the value explored in the next right branch
    void postCut( const int p );
    //@}
  };



  /**********************************************
   * MistralGacList
   **********************************************/
  /// List of Variable for computing the GAC closure (fifo).
  template< class T >
  class MistralGacList {
  public:

    /**@name Constructors*/
    //@{
    T   *X;
    int *next;
    int *trigger;
    int head;
    int tail;
    //@}

    /**@name Constructors*/
    //@{
    MistralGacList() {
      trigger = NULL;
      next = NULL;
    }
    ~MistralGacList()
    {
      delete [] trigger;
      delete [] next;
    }
    void initList(int n, T *x)
    {
      X = x;
      trigger = new int[n];
      std::fill( trigger, trigger+n, 0 );
      next = new int[n+2];
      std::fill( next, next+n, NOVAL );
      head = tail = n;
      next[n] = n;
    }
    //@}
  
    /*!@name Accessors*/
    //@{
    /// Add an element to the queue
    inline void push(const int i, const int event)
    {
      //int i=v->id;
      if( !trigger[i] ) {
	next[tail] = i;
	tail = i;
	next[tail] = head;
      }
      trigger[i] |= event;
    }

    /// Pop the first event of the queue
    inline T pop( int& event )
    {
      int i=next[head];
      T x = X[i];
      next[head] = next[i];
      event = trigger[i];
      trigger[i] = 0;

      if( next[head] == head )
	tail = head;
      return x;
    }

    /// Is the queue empty?
    inline bool empty() const
    {
      return (next[head] == head);
    }

    /// Clear all elements of the queue
    inline void clear()
    {
      while( next[head] != head ) {
	trigger[next[head]] = 0;
	next[head] = next[next[head]];
      }
      tail = head;
    }
    //@}

    /*!@name Miscellanous*/
    //@{ 
    void print(std::ostream& o) const
    {
      int i;
      std::cout << "{";
      i = next[head];
      while( i != head ) {
	X[i]->printshort( std::cout );
	std::cout << " ";
	std::cout.flush();
	i = next[i];
      }
      std::cout << "}";
    }
    //@}

  };



  /**********************************************
   * Constraint
   **********************************************/
  /*! \class Constraint
    \brief Representation of Constraints.

    The constrained variables are stored in the
    array _scope[]. 
    The method propagate() is called suring preprocessing.
    The method propagate(const int, const int) is
    called when computing the GAC closure. 
    Generic propagate methods, using AC3
    with residual supports are implemented in
    the class Constraint.
  */
  class Constraint : public ReversibleObj {

  public:

    /*!@name Static members*/
    //@{
    static const int WEIGHTONE      = 1;
    static const int WEIGHTTWO      = 1;
    static const int WEIGHTTHREE    = 1;
    static const int DOMAINTRIGGER  = 1;
    static const int RANGETRIGGER   = 2;
    static const int VALUETRIGGER   = 4;
    static const int DELAYTRIGGER   = 8;
    //@}


    /*!@name Parameters*/
    //@{
    /// Whether the propagation is delayed
    bool delayed;
    /// An unique id
    int id;
    /// The number of constrained variables
    int arity;
    /// Weight, used in some heuristics
    long unsigned int weight;
    /// Used to keep the weight consistent with current constraints
    long unsigned int oldweight;
    /// Hack for the linkage
    int lastlinked;

    /// The constrained 'real' variables.
    VariableInt** _scope;
    /// The variables whose domains are filtered (may be references)
    VariableInt** scope;
    /// The indices of the two variables that are watching us
    VariableInt* watch[2];

    /// The type of trigger for each constrained variable
    /*!
      The possible triggers, for the ith variable in the _scope are:
      - NULL
      - &(_scope[i]->valueTrigger)
      - &(_scope[i]->rangeTrigger)
      - &(_scope[i]->domainTrigger)
    */
    MistralList<Constraint*> **triggers;
    /// The element of the list
    /*!
      This node points to this constraint and
      can be added or removed from the corresponding
      list in constant time.
    */
    MistralNode<Constraint*> *elements;
  
    /// used to manipulate supports in the generic AC algorithm
    int *supp; 
    /// structure to keep residual supports
    int ***supports; 
    //@}


    /*!@name Constructors*/
    //@{
    /// The _scope is build by copying an array "scp" containing "l" variables
    //void relax();
    void initScope(VariableInt **scp, const int l, const int t, const int w=1);
    Constraint() {}
    Constraint(Solver *s);
    Constraint(Solver *s, VariableInt **scp, const int l, const int t=DOMAINTRIGGER, const int w=1);
    virtual ~Constraint();
    //@}

    /*!@name Propagators*/
    //@{
    /*!
     * This methods returns 0 if the constraint is satisfied
     * by the combination and any positive value otherwise.  
     */
    virtual int check(const int*) const = 0;
    /*!
     * Called during preprocessing
     */
    virtual bool propagate();
    /*! 
     *  This method is called when the domain of the 
     *  variable _scope[idx] has been modified. Some domains in 
     *  _scope may be reduced, and true 
     *  is returned on success (there are still at least one 
     *  consistent assignment) or false otherwise. 
     */
    virtual bool propagate(const int idx, const int e);
    /*!
     *  Check if the cached support is valid.
     *  Otherwise initialize an iterator on supports
     */
    bool firstSupport(const int, const int);
    /*!
     *  Find the first valid support. 
     *  Return true if a support has been found, 
     *  or false if no such support exists.
     */
    bool findSupport(const int, const int);
    //@}


    /*!@name Manipulation methods*/
    //@{
    /// entailment
    void restore();
    /// forced disentailment
    void relax();

    bool isLinked();
    //   {
    //     return ( watch[0]->isLinked() && watch[1]->isLinked() );
    //   }

    inline void stopWatching( VariableInt* );
    inline void resumeWatching( VariableInt* );
    /// Whether the residual support should be used
    void useResidualSupports() ;
    /// Print the constraint
    virtual void print(std::ostream & os) const ;
    virtual void printshort(std::ostream & os) const ;
    bool constrains( VariableInt *x ) const ;
    virtual std::string name() { return "noname"; }
    //@}
  };


  /**********************************************
   * MistralGraph
   **********************************************/
  /*! \class MistralGraph
    \brief  A directed graph implementation.
  */
  class MistralGraph
  {
  public:

    /*!@name Parameters*/
    //@{
    int n;
    int **neighbor;
    int ***node;
    int **first;
    int **last;
    //@}

    /*!@name Constructors*/
    //@{
    MistralGraph( VariableInt **vars, const int l );
    MistralGraph( CSP& model ); //TODO
    MistralGraph( const int l );
    MistralGraph( const int l, std::vector<int> E );
    MistralGraph( ) { n = 0; }
    void init(const int);
    void cleanUp();
    virtual ~MistralGraph();
    //@}

  private:

    int getAndRemoveMaxCand( BitSet& cand );
    void removeNeighbours  ( const int x, BitSet& neighbours );
    void updateGammaK      ( const int x, BitSet& res );
    void storeClique       ( BitSet& K, std::vector< std::vector< int > >& cliques );
    void BronKerbosch      ( BitSet& cand,
			     BitSet& K,
			     BitSet& gammaK,
			     std::vector< std::vector< int > >& cliques,
			     const int limit );

  public:
    /*!@name Algorithms*/
    //@{
    void getAllCliques( std::vector< std::vector< int > >& cliques, const int limit=NOVAL );
    void getSomeCliques( std::vector< std::vector< int > >& cliques );
    void dfs(const int) const;
    int AltBlumMehlorn( int nv );
    int greedyMaxClique() const;
    //@}

    /*!@name Accessors*/
    //@{
    inline bool isIn(const int x, const int y) const
    {
      return (bool)(node[x][y]);
    }

    int next(const int i, const int j)
    {
      return *(node[i][j]+1);
    }
    int degree(const int) const;
    void add(const int, const int);
    void remove(const int, const int);
    //@}

    /*!@name Manipulation methods*/
    //@{  
    void print(std::ostream& o) const;
    //@}
  };

};

#endif // __CSP_H




/*! \mainpage Mistral
 *
 * \section intro_sec Introduction
 *
 * Welcome to <i><b>Mistral</b></i>'s web page.
 * <i><b>Mistral</b></i> is a C++ library for solving
 * Constraint satisfaction and optmisation problems.
 * If you are unsure about Constraint satisfaction,
 * you may want to check the <a href="http://kti.ms.mff.cuni.cz/~bartak/Constraints/index.html">online guide to Constraint Programming</a> by Roman Bart&aacute;k first.
 * This documentation has been
 * generated using <a href="http://www.doxygen.org/index.html">doxygen</a>.
 *
 * 
 *
 *
 * \section install_sec Installation
 * 
 * In principle, the makefile should
 * work ok as is.
 * <br>
 * To compile the library and all examples:
 * <br>
 * $ <tt>make</tt>
 * <br>
 * To compile the model <tt>model-name</tt> (if
 * the file <tt>models/src/model-name.cpp</tt> exists):
 * <br>
 * $ <tt>make model-name</tt>
 * <br>
 * To compile the library only:
 * <br>
 * $ <tt>make lib</tt>
 * <br><br>
 * The directory <b>include/</b> contains all headers 
 * and <b>lib/src/</b> all sources for the library.
 * The object files go into <b>lib/obj/</b> and the binaries
 * into <b>bin/</b>. Last, the directory <b>models/src/</b>
 * contains the models (examples).
 *
 * \section features_sec Features
 * 
 * Some features implemented in Mistral:
 * <ul>
 * <li> Neighborhood based heuristic (Bessiere, Chmeiss and Saïs)
 * <li> Weighted Degree heuristic (Lestd::coutre et al.)
 * <li> Bound Consistency for the AllDifferent Constraint 
 * (Code from Lopez-Ortiz, Quimper, Tromp and van Beek)
 * <li> NValue Constraint (Beldiceanu) <!--<li> Singleton Arc Consistency (Bessiere and Debruyne) and Shaving (Torres and Lopez)-->
 * <li> EdgeFinder Algorithm implemented as a global Constraint (Carlier) <!--<li> Operation Resource Reliance Heuristic for job-shop scheduling (Sadeh)-->
 * <li> PCP branching algorithm for job-shop scheduling (Smith and Cheng)
 * <li> ...
 *
 *</ul>
 * \section tutorial_sec Tutorial
 *
 * This \ref tutorial shows step by step how to write a 
 * model for the N queens problem using Mistral.
 * It also shows how to declare new Constraints
 * and propagators.
 *
 * \section copyright_sec Copyright
 *
 * This software is under GNU 
 *  <a href="http://www.gnu.org/copyleft/gpl.html">General Public License</a>.
 *  <br><br>
 *  Mistral, a Constraint satisfaction and optimisation library
 *  Copyright (C) 2003-2005  <br> Emmanuel Hebrard
 *  <br><br>
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *  <br><br>
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  <br><br>
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 * <br><br>
 * The author can be contacted electronically at 
 * ehebrard@cse.unsw.edu.au.
 */

/*! \page tutorial Tutorial

  <A NAME="top"></A>
  \section golomb_sec Golomb Ruler Problem

  The following is a simple model for the Golomb Ruler problem.


  \code

#include <sol.h>
 
using namespace std;

// Known optimal rulers. Rulers of lower order are used as lower bound for distances.
static const int ruler[] = {
  0,0,1,3,6,11,17,25,34,
  44,55,72,85,106,127};


int main(int argc, char *argv[])
{  
  int nbMarks = (argc > 1 ? atoi(argv[1]) : 10);
  int rulerSize = ( 2 << (nbMarks-1) ); // naive upper bound for the ruler size

  // variables
  CSP model;
  VarArray mark(nbMarks, 0, rulerSize-1);  
  VarArray distance(nbMarks*(nbMarks-1)/2, 1, rulerSize-1);

  // constraints
  int i, j, k=0;
  for(i=1; i<nbMarks; ++i) {
    
    // strict ordering of the marks
    model.add( mark[i-1] < mark[i] );

    // lower bounds from earlier golomb rulers
    model.add( mark[i-1] >= ruler[i] );

    for(j=0; j<i; ++j) {

      // lower bounds from earlier golomb rulers
      model.add( distance[k] >= ruler[i-j] );

      // set up the distances
      model.add( mark[i] == (mark[j] + distance[k++]) );
    }
  }

  // distances are all different
  model.add( BoundAllDifferent(distance) );

  // symmetry breaking
  model.add( mark[0] == 0 );
  if( nbMarks > 2 )
    model.add( distance[0] < distance[k-1] );

  // set up objective function
  model.add( Minimise( mark[nbMarks-1] ) );


  // solver declaration
  Solver s( model, mark );
  s.setDomainSplitting();
  s.solve();


  s.printStatistics( std::cout, (Solver::NDS | Solver::RUNTIME | Solver::SPEED) );
  std::cout << std::endl;
}

  \endcode

  Here is an example of output 

  <pre>

o         82 NDS         0 s    solution(80): {1 3 7 12 20 30 44 65 80}
o         89 NDS         0 s    solution(75): {1 3 7 12 20 34 49 59 75}
o        101 NDS         0 s    solution(73): {1 3 7 12 22 35 49 65 73}
o        123 NDS      0.01 s    solution(72): {1 3 7 12 26 41 54 62 72}
o        248 NDS      0.02 s    solution(70): {1 3 7 15 24 34 54 59 70}
o        419 NDS      0.03 s    solution(68): {1 3 7 15 31 36 49 58 68}
o        605 NDS      0.04 s    solution(66): {1 3 7 17 22 35 46 58 66}
o        715 NDS      0.05 s    solution(62): {1 3 7 18 30 38 43 52 62}
o       2480 NDS      0.18 s    solution(60): {1 3 11 17 29 36 51 56 60}
o      12392 NDS      0.87 s    solution(55): {1 6 10 23 26 34 41 53 55}
    38826 NDS    14326 BTS/s    2.71 s

  </pre>


*/
