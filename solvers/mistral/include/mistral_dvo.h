 
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

/** \file dvo.h
    \brief Header for DVOs (Dynamic Variable/Value Ordering) and Goals (Objective functions).
*/


#ifndef _DVO_H
#define _DVO_H

#include <mistral_con.h>
#include <mistral_gen.h>
#include <mistral_var.h>

void free_disjuncts();

namespace Mistral {



  class SimpleUnaryConstraint {
  public:
    // 0 -> removal 'n'
    static const int    REMOVAL = 0;
    // 1 -> assignment 'e'
    static const int ASSIGNMENT = 1;
    // 2 -> lower bound 'g'
    static const int LOWERBOUND = 2;
    // 3 -> upper bound 'l'
    static const int UPPERBOUND = 3;
    
    /**
       2 bits for the type 
       30 bits for the value
     */

    unsigned int _data_;
    VariableInt* var;

    inline int type() const {return _data_&3; }
    inline int value() const {return _data_>>2; }

    inline bool satisfied() {
      switch(type()) {
      case REMOVAL: return !var->contain(value());
      case ASSIGNMENT: return var->equal(value());
      case LOWERBOUND: return var->min() >  value(); // > v
      case UPPERBOUND: return var->max() <= value();   // <= v
      }
      return true;
    }

    inline bool violated() {
      switch(type()) {
      case REMOVAL: return var->equal(value());
      case ASSIGNMENT: return !var->contain(value());
      case LOWERBOUND: return var->max() <= value(); // > v
      case UPPERBOUND: return var->min() >  value();   // <= v
      }
      return true;
    }
    

    SimpleUnaryConstraint(const char t, const int v, VariableInt *x) {
      init_data(t,v);
      var = x;
    }
    
    void init_data(const char t, const int v) {
      _data_ = ((v - (t == 'g')) << 2);
      switch(t) {
      case 'n': _data_ |= REMOVAL;
      case 'e': _data_ |= ASSIGNMENT;
      case 'g': _data_ |= LOWERBOUND;
      case 'l': _data_ |= UPPERBOUND;
      }
      //type = t;
      //val = v;
    }

    SimpleUnaryConstraint(VariableInt *x) {
      _data_ = 0xffffffff;
      //type = 'x';
      //val = NOVAL;
      var = x;
    }

    SimpleUnaryConstraint() {
      _data_ = 0xffffffff;
      //type = 'x';
      //val = NOVAL;
      var = NULL;
    }

    virtual ~SimpleUnaryConstraint() {
    }

    inline void revert() { _data_^=1; }
    inline void make();
//  { 
//       int v;
//       int t;
//       var->branch->make(t,v); 
//       _data_ = ((v<<2)+t);
//     }
    inline bool left() {
      switch(type()) {
      case REMOVAL: return var->remove(value());
      case ASSIGNMENT: return var->setDomain(value());
      case LOWERBOUND: return var->setMin(value()+1); // > v
      case UPPERBOUND: return var->setMax(value());   // <= v
      }
      return true;
    }
    inline bool right() { 
      revert();
      return left();
    }

    inline bool propagate() {
      return left();
    }

    std::ostream& print(std::ostream& os) const {
      var->printshort(os);
      switch(type()) {
      case REMOVAL: {
	os << " =/= " << value();
      } break;
      case ASSIGNMENT: {
	os << " == " << value();
      } break;
      case LOWERBOUND: {
	os << " > " << value();
      } break;
      case UPPERBOUND: {
	os << " <= " << value();
      } break;
      }

      return os;
    }

  };
  typedef SimpleUnaryConstraint Decision;

  /**********************************************
   * Branching Strategy
   **********************************************/

  /// Wrapper used to add a value heuristic to a bunch of variables
  /**
     It implements a method "extract()" that 
     assign a value heuristic to the variables
     in the array scope
  */
  class BranchingStrategy {

  public:
    BranchingStrategy();
    virtual ~BranchingStrategy();
    
    virtual void extract() = 0;
  };

  class VarArray;
  class Contention : public BranchingStrategy {
    
  public:
//     double ***proba;
//     double ***tmp_proba;
//     VariableInt ***scope;
    int width;
    int length;

    int *min_dx;
    double **proba_; // vars x vals : probability
    double **tmp_proba_; // vars x vals : probability
    VariableInt **scope_; // vars
    Vector<int> active; // non-singleton vars
    Vector<int>*neighborhood; // var indices x neighbor index : neighbors
    int num_vars;


    Contention(VarArray* scp, int w);
    //Contention(BuildObject*** scp, int w, int l);
    virtual ~Contention();
    
    void compute_proba();
    void reset();
    int get_best(const int id);
    int update(const int id);
    std::ostream& print(std::ostream&) const;
    std::ostream& print(std::ostream&, const int) const;

    virtual void extract();
  };




  class Weighter
  {
  public:

    /**@name Parameters*/
    //@{ 
    static const int NO=0;
    static const int WDG=1;
    static const int WLD=3;
    static const int IPT=8;
    static const int SAC=16;
    static const int PRU=32;
    static const int ISAC=64;
    static const int RNGD=128;
    static const int RGNGD=256;

    int& level;

    int init_level;
    //@}

    /**@name Constructors*/
    //@{
    Weighter( Solver* );
    virtual ~Weighter() {}
    //@}  

    /**@name Utils*/
    //@{ 
    inline void notify_init_level( const int il ) { init_level = il; }
    virtual void notifyFailure( Constraint *con ) {}
    virtual void notifySuccess() {}
    virtual void notifyChoice() {}  
    virtual void notifyRestart() {}  
    virtual int getType() { return NO; }
    virtual void print( std::ostream& o ) const {}
    //@}  
  };

  class WeighterDegree : public Weighter
  {
  public:

#ifdef _WEIGHT_STATS

    int **choice;
    int numvars;
    int numchoices;
    Solver *solver;

#endif

    int threshold;

    /**@name Constructors*/
    //@{  
    WeighterDegree( Solver* s );  
    virtual ~WeighterDegree( ) {


#ifdef _WEIGHT_STATS
      //       delete [] order;

      for(int i=0; i<=numvars; ++i)
	delete [] choice[i];
      delete [] choice;

#endif 
      //

    }
    //@}  

    /**@name Utils*/
    //@{
#ifdef _WEIGHT_STATS

    virtual void notifyChoice();  
    void print_choices(std::ostream& o);
    void init_choices();

#endif 

    virtual void notifyFailure( Constraint *con ) ;
    virtual int getType() { return WDG; }
    //@}  
  };


  class WeighterLevelDegree : public Weighter
  {
  public:

    /**@name Parameters*/
    //@{ 
    /// Depth of the search tree
    int lmax;
    //@}

    /**@name Constructors*/
    //@{
    WeighterLevelDegree( Solver* s );  
    virtual ~WeighterLevelDegree( ) {}
    //@}  

    /**@name Utils*/
    //@{ 
    virtual void notifyFailure( Constraint *con ) ;
    virtual int getType() { return WLD; }
    //@}  
  };


  // class WeighterPruning : public Weighter
  // {
  //  public:

  //   /**@name Constructors*/
  //   //@{
  //   WeighterPruning( Solver* s );  
  //   virtual ~WeighterPruning( ) {}
  //   //@}  
 
  //   /**@name Utils*/
  //   //@{ 
  //   virtual void notifyFailure( Constraint *con ) ;
  //   virtual void notifySuccess() ;
  //   virtual void notifyChoice() ; 
  //   virtual int getType() { return PRU; }
  //   //@}  
  // };


  class WeighterImpact : public Weighter
  {
  public:

    /**@name Parameters*/
    //@{ 
    // repetition of some of the dvo's params
    /// The list of decisions
    VariableInt**& decision;
    /// A reference to the beginning of the array of future variables
    VariableInt**& first;
    /// A reference to the end of the array of future variables
    VariableInt**&  last;

    // original length of the variable array
    unsigned int ilength;

    // let Sb be the size of the CSP before decision d,
    // and Sa be the size of the CSP after this decision.
    // I(d) = Sb/Sa, we store 1/I(d), that is Sa/Sb in [0..1].
    double **decision_impact;
    bool *isRange;

    // number of time this decision was taken (to compute the average)
    unsigned int **decision_count;

    // use to store the variables and domain size 
    // when exploring the future vars.
    // when we are notified that the filtering occurred, 
    // we can compare the size in order to update impacts
    VariableInt  **X;
    unsigned int *domain_size;
    unsigned int nbvars;

    int *varmin;
    bool needUpdate;

    VariableInt** variables;
    //@}

    /**@name Constructors*/
    //@{
    WeighterImpact( Solver* s );  
    virtual ~WeighterImpact( );
    //@}

    /**@name Utils*/
    //@{ 
    virtual void notifyFailure( Constraint *con ) ;
    virtual void notifySuccess() ;
    virtual void notifyChoice() ; 
    virtual int getType() { return IPT; }
    virtual void print( std::ostream& o ) const;
    //@}
  };


  class WeighterSAC : public Weighter
  {
  public:

    /**@name Parameters*/
    //@{ 
    bool         complete;
    bool              ins;
    bool         *isRange;
    bool      needPruning;
    int           ilength;
    int           domsize;
    int           sacvals;
    Solver        *solver;
    BitSet     *SACdomain;
    //@}

    /**@name Constructors*/
    //@{
    WeighterSAC( Solver *s );
    virtual ~WeighterSAC( );
    //@}  

    /**@name Utils*/
    //@{ 
    bool isNotSac();

    virtual void notifyFailure( Constraint *con ) ;
    virtual void notifyChoice() ; 
    virtual int getType() { return SAC; }
    //@}  
  };


  class WeighterISAC : public WeighterSAC
  {
  public:

    /**@name Constructors*/
    //@{
    WeighterISAC( Solver *s );
    virtual ~WeighterISAC( );
    //@}
 
    /**@name Utils*/
    //@{
    virtual void notifyFailure( Constraint *con ) ;
    //@} 

  };


  class WeighterSAT : public Weighter
  {
  public:

    /**@name Parameters*/
    //@{
    /// Pointer to the sat solver
    ConstraintClauseBase *sat;
    /// A reference to the beginning of the array of future variables
    VariableInt**& first;
    //@}

    /**@name Constructors*/
    //@{
    WeighterSAT( Solver* s, ConstraintClauseBase *c );  
    virtual ~WeighterSAT( );
    //@}  

    /**@name Utils*/
    //@{ 
    virtual void notifyFailure( Constraint *con ) ;
    virtual void notifyChoice() ;
    //@}  
  };



  class WeighterRestartNogood : public Weighter
  {
  public:

    /**@name Parameters*/
    //@{
    ConstraintClauseBase *sat;
    VariableInt **decision;

    int       *choices;
    Vector< int > path;
    int lvl;
    //int         *nBranches;
    //@}

    /**@name Constructors*/
    //@{
    WeighterRestartNogood( Solver* s );  
    virtual ~WeighterRestartNogood( );
    //@}  

    /**@name Utils*/
    //@{ 
    virtual void notifyFailure( Constraint *con ) ;
    virtual void notifyChoice() ;
    virtual void notifyRestart() ; 
    //@}  
  };


  class WeighterRestartGenNogood : public Weighter
  {
  public:

    /**@name Parameters*/
    //@{
    ConstraintGenNogoodBase *base;
    Decision *decision;

    // record the bad decisions at a given level
    Vector< Decision >* bad_choices;

    // depth of the tree to parse
    int depth;
    //@}

    /**@name Constructors*/
    //@{
    WeighterRestartGenNogood( Solver* s );  
    virtual ~WeighterRestartGenNogood( );
    //@}  

    /**@name Utils*/
    //@{ 
    void forget(const int l);
    void reinit();

    virtual void notifyChoice() ; 
    virtual void notifyFailure( Constraint *con ) ;
    virtual void notifyRestart() ; 
    //@}  
  };


  class WeighterSwitchPromise : public Weighter
  {
  public:

    /**@name Parameters*/
    //@{
    int length;
    int threshold;
    int n_restart;
    VariableInt **sequence;
    //@}

    /**@name Constructors*/
    //@{
    WeighterSwitchPromise( Solver* s, const int thresh );  
    virtual ~WeighterSwitchPromise( );
    //@}  

    /**@name Utils*/
    //@{ 
    virtual void notifyRestart() ; 
    //@}  
  };


  class VarSelectorDomain 
  {
  public: 

    /**@name Constructors*/
    //@{
    VarSelectorDomain() {d_ = NOVAL;}
    //@}

    double value() { return 1.0/((double)d_); }

    /**@name Parameters*/
    //@{ 
    int d_;
    //@}  

    /**@name Utils*/
    //@{ 
    inline bool operator<( VarSelectorDomain& x ) const { return d_ < x.d_; }
    inline void operator=( VarSelectorDomain& x ) { d_ = x.d_; }
    inline void operator=( VariableInt    *x ) { d_ = x->domsize(); }
    //@}  

    std::ostream& print(std::ostream& os) const {
      os << "-" ;
      return os;
    }
  };

  class VarSelectorDomainMinVal 
  {
  public: 

    /**@name Constructors*/
    //@{
    VarSelectorDomainMinVal() {d_ = NOVAL; v_ = NOVAL;}
    //@}

    double value() { return 1.0/((double)d_); }

    /**@name Parameters*/
    //@{ 
    int d_;
    int v_;
    //@}  

    /**@name Utils*/
    //@{ 
    inline bool operator<( VarSelectorDomainMinVal& x ) const { return d_ < x.d_  || ( d_ <= x.d_ && v_ < x.v_) ;}
    inline void operator=( VarSelectorDomainMinVal& x ) { d_ = x.d_; v_ = x.v_; }
    inline void operator=( VariableInt    *x ) { d_ = x->domsize(); v_ = x->min(); }
    //@}  

    std::ostream& print(std::ostream& os) const {
      os << "-" ;
      return os;
    }
  };

  class VarSelectorDegree 
  {
  public: 

    /**@name Constructors*/
    //@{
    VarSelectorDegree() {d_ = 0;}
    //@}

    double value() { return 1.0/((double)d_); }

    /**@name Parameters*/
    //@{ 
    int d_;
    //@}    

    /**@name Utils*/
    //@{ 
    inline bool operator<( VarSelectorDegree& x ) const { return d_ > x.d_; }
    inline void operator=( VarSelectorDegree& x ) { d_ = x.d_; }
    inline void operator=( VariableInt    *x ) { d_ = x->degree; }
    //@}  

    std::ostream& print(std::ostream& os) const {
      os << "-" ;
      return os;
    }
  };

  class VarSelectorDomainDegree 
  {
  public: 

    /**@name Constructors*/
    //@{
    VarSelectorDomainDegree() {dom_ = NOVAL; deg_ = 0;}
    //@}

    /**@name Parameters*/
    //@{ 
    int deg_;
    int dom_;
    //@}  

    double value() { return 1.0/((double)dom_); }
  
    /**@name Utils*/
    //@{ 
    inline bool operator<( VarSelectorDomainDegree& x ) const { 
      return dom_ < x.dom_  || ( dom_ <= x.dom_ && deg_ > x.deg_) ; }   
    inline void operator=( VarSelectorDomainDegree& x ) { deg_ = x.deg_; dom_ = x.dom_; }
    inline void operator=( VariableInt    *x ) { deg_ = x->degree; dom_ = x->domsize(); }
    //@}  

    std::ostream& print(std::ostream& os) const {
      os << "-" ;
      return os;
    }
  };

  class VarSelectorDomainOverDegree 
  {
  public: 

    /**@name Constructors*/
    //@{
    VarSelectorDomainOverDegree() {dom_ = NOVAL; deg_ = 0;}
    //@}

    double value() { return 1.0/((double)dom_ / (double)deg_); }

    /**@name Parameters*/
    //@{ 
    int deg_;
    int dom_;
    //@}  
  
    /**@name Utils*/
    //@{ 
    inline bool operator<( VarSelectorDomainOverDegree& x ) const { return dom_ * x.deg_ < x.dom_ * deg_ ; }
    inline void operator=( VarSelectorDomainOverDegree& x ) { deg_ = x.deg_; dom_ = x.dom_; }
    inline void operator=( VariableInt    *x ) { deg_ = x->degree; dom_ = x->domsize(); }
    //@}  

    std::ostream& print(std::ostream& os) const {
      os << "-" ;
      return os;
    }
  };

  class VarSelectorWeight 
  {
  public: 

    /**@name Constructors*/
    //@{
    VarSelectorWeight() {w_ = 0;}
    //@}

    double value() { return ((double)w_); }

    /**@name Parameters*/
    //@{ 
    int w_;
    //@}  

    /**@name Utils*/
    //@{ 
    inline bool operator<( VarSelectorWeight& x ) const { return w_ > x.w_; }
    inline void operator=( VarSelectorWeight& x ) { w_ = x.w_; }
    inline void operator=( VariableInt    *x ) { w_ = x->weight; }
    //@}  

    std::ostream& print(std::ostream& os) const {
      os << "-" ;
      return os;
    }
  };

  class VarSelectorDomainOverWeight 
  {
  public: 

    /**@name Constructors*/
    //@{
    VarSelectorDomainOverWeight() { dom_ = NOVAL; wgt_ = 0; }
    //@}

    double value() { return ((double)wgt_); }

    /**@name Parameters*/
    //@{ 
    int wgt_;
    int dom_;
    //@}  
  
    /**@name Utils*/
    //@{ 
    inline bool operator<( VarSelectorDomainOverWeight& x ) const { return dom_ * x.wgt_ < x.dom_ * wgt_ ; }
    inline void operator=( VarSelectorDomainOverWeight& x ) { wgt_ = x.wgt_; dom_ = x.dom_; }
    inline void operator=( VariableInt    *x ) { 
      wgt_ = x->weight; 
      dom_ = x->domsize(); 
//       std::cout << " =" ;
//       x->print(std::cout);
//       std::cout << " " ;
//       std::cout.flush();
    }
    //@}  

    std::ostream& print(std::ostream& os) const {
      os << dom_ << "/" ;
      os << wgt_ ;
      return os;
    }
  };


  class VarSelectorImpact 
  {
  public: 

    /**@name Constructors*/
    //@{
    VarSelectorImpact() {impact_ = NOVAL;}
    //@}

    double value() { return 1/((double)impact_); }

    /**@name Parameters*/
    //@{ 
    double impact_;
    double **decision_impact;
    //@}  
  
    /**@name Utils*/
    //@{ 
    inline bool operator<( VarSelectorImpact& x ) const { return impact_ < x.impact_; }
    inline void operator=( VarSelectorImpact& x ) { impact_ = x.impact_; }
    inline void operator=( VariableInt    *x ) 
    { 
      impact_ = 0.0;
      int idx = x->id;

      if( x->getType() != VariableInt::RANGE ) {
	DomainIterator *valit = x->begin();
	do {
	  impact_ += decision_impact[idx][*valit]; 
	} while( valit->next() );    
      } else {
	impact_ = ((decision_impact[idx][0] + decision_impact[idx][1]) * x->domsize() / 2);
      }
    }
    //@}  

    std::ostream& print(std::ostream& os) const {
      os << "-" ;
      return os;
    }
  };

  class VarSelectorImpactOverDegree 
  {
  public: 

    /**@name Constructors*/
    //@{
    VarSelectorImpactOverDegree() {impact_ = NOVAL; deg_ = 0;}
    //@}

    double value() { return 1/((double)impact_ / (double)deg_); }

    /**@name Parameters*/
    //@{ 
    int deg_;
    double impact_;
    double **decision_impact;
    //@}  
  
    /**@name Utils*/
    //@{ 
    inline bool operator<( VarSelectorImpactOverDegree& x ) const { 
      return impact_ * x.deg_ < x.impact_ * deg_; }
    inline void operator=( VarSelectorImpactOverDegree& x ) { 
      impact_ = x.impact_; deg_ = x.deg_; }
    inline void operator=( VariableInt    *x ) 
    { 
      impact_ = 0.0;
      int idx = x->id;

      DomainIterator *valit = x->begin();
      do impact_ += decision_impact[idx][*valit]; 
      while( valit->next() ); 

      deg_ = x->degree;
    }
    //@}  

    std::ostream& print(std::ostream& os) const {
      os << "-" ;
      return os;
    }
  };

  class VarSelectorImpactOverWeight 
  {
  public: 

    /**@name Constructors*/
    //@{
    VarSelectorImpactOverWeight() {impact_ = NOVAL; wgt_ = 0;}
    //@}

    double value() { return 1/((double)impact_ / (double)wgt_); }

    /**@name Parameters*/
    //@{ 
    int wgt_;
    double impact_;
    double **decision_impact;
    //@}  
  
    /**@name Utils*/
    //@{ 
    inline bool operator<( VarSelectorImpactOverWeight& x ) const { 
      return impact_ * x.wgt_ < x.impact_ * wgt_; }
    inline void operator=( VarSelectorImpactOverWeight& x ) { 
      impact_ = x.impact_; wgt_ = x.wgt_; }
    inline void operator=( VariableInt    *x ) 
    { 
      impact_ = 0.0;
      int idx = x->id;

      DomainIterator *valit = x->begin();
      do impact_ += decision_impact[idx][*valit]; 
      while( valit->next() ); 

      wgt_ = x->weight;
    }
    //@}

    std::ostream& print(std::ostream& os) const {
      os << "-" ;
      return os;
    }
  };

  class VarSelectorNeighborDD 
  {
  public: 

    /**@name Constructors*/
    //@{ 
    VarSelectorNeighborDD() {val_ = NOVAL;}
    //@}

    double value() { return 1/((double)val_); }

    /**@name Parameters*/
    //@{ 
    double val_;
    VariableInt **scp;
    int j;
    //@}  

    /**@name Utils*/
    //@{ 
    inline bool operator<( VarSelectorNeighborDD& x ) const { 
      return val_ < x.val_ ; }
    inline void operator=( VarSelectorNeighborDD& x ) { 
      val_ = x.val_; }
    inline void operator=( VariableInt    *x ) 
    { 
      val_ = 0.0;
      //MistralNode<Constraint*> *nd = &(x->valueTrigger.head);
      MistralNode<Constraint*> *nd = x->constraintsOnValue();
      while( nextNode(nd) ) {
	scp = nd->elt->scope;
	j=nd->elt->arity;
	while( j-- )
	  if( scp[j]->isLinked() && 
	      scp[j]->degree && 
	      scp[j] != x )
	    val_ += ( (double)(scp[j]->domsize()) / 
		      (double)(scp[j]->degree) );
      }
      val_ = ((double)(x->domsize() + val_) /
	      (double)(x->degree*x->degree));
    }
    //@}

    std::ostream& print(std::ostream& os) const {
      os << "-" ;
      return os;
    }
  };


  class VarSelectorOSP 
  {
  public: 

    /**@name Constructors*/
    //@{
    VarSelectorOSP() {impact_ = 0; weight_ = 0; weight_task_ = 0;}
    //@}

    double value() { return ((double)weight_task_); }

    /**@name Parameters*/
    //@{ 
    int domsize_;
    int impact_;
    int weight_;
    int weight_task_;
    PredicateDisjunctive **disjuncts;
    int strategy;
    //@}  
  
    /**@name Utils*/
    //@{ 
    inline bool operator<( VarSelectorOSP& x ) const { 
      //return ( (weight_ > x.weight_) || ( weight_ == x.weight_ && (domsize_ * x.impact_ < x.domsize_ * impact_)) ) ; 
      //return ( (domsize_ - impact_) * x.weight_ < (x.domsize_ - x.impact_) * weight_ );
      //return ( (weight_ > x.weight_) || (weight_ == x.weight_ && (domsize_ - impact_) < (x.domsize_ - x.impact_)) ) ;
      //return ( (weight_ > x.weight_) || ( weight_ == x.weight_ && (domsize_ < x.domsize_)) ) ;  
    
      //return (domsize_ * x.weight_ < x.domsize_ * weight_) ;

      //return (domsize_ * (x.weight_ + x.weight_task_) < x.domsize_ * (weight_ + weight_task_)) ;

      //std::cout << "osp d/tw" << std::endl;

      return (domsize_ * (x.weight_task_ // - x.weight_
			  ) < x.domsize_ * (weight_task_ // - weight_
					    )) ;


    }
    inline void operator=( VarSelectorOSP& x ) { // impact_ = x.impact_;
      weight_task_ = x.weight_task_;
      weight_ = x.weight_; 
      domsize_ = x.domsize_; 
    }
    inline void operator=( VariableInt    *x ) 
    { 
      domsize_ = disjuncts[x->id]->domsize();
      weight_task_ = (disjuncts[x->id]->scope[0]->weight
		      +disjuncts[x->id]->scope[1]->weight);
      weight_ = x->weight;
    }
    //@}  

    std::ostream& print(std::ostream& os) const {
      os << "-" ;
      return os;
    }
  };


  class VarSelectorPFSP
  {
  public: 
    
    /**@name Constructors*/
    //@{
    VarSelectorPFSP() {domsize_ = 0; weight_ = 0;}
    //@}

    double value() { return (double)(weight_); }

    /**@name Parameters*/
    //@{ 
    int domsize_;
    int weight_;
    PredicateDisjunctive ***disjuncts;
    int *sdegree;
    //@}  
  
    /**@name Utils*/
    //@{ 
    inline bool operator<( VarSelectorPFSP& x ) const { 
      return (domsize_ * x.weight_ < x.domsize_ * weight_) ;
    }
    inline void operator=( VarSelectorPFSP& x ) { 
      weight_ = x.weight_; domsize_ = x.domsize_; 
    }
    inline void operator=( VariableInt    *x ) 
    { 
      int i=x->id, j;
      domsize_ = 0;
      
      // sum of the domain size:
      for(j=0; j<sdegree[i]; ++j) 
	domsize_ += disjuncts[i][j]->domsize();

      weight_ = x->weight;
    }
    //@}  

    std::ostream& print(std::ostream& os) const {
      os << "-" ;
      return os;
    }
  };

  class VarSelectorOSP_Domain
  {
  public: 

    /**@name Constructors*/
    //@{
    VarSelectorOSP_Domain() {domsize_ = 0;}
    //@}

    double value() { return (double)domsize_; }

    /**@name Parameters*/
    //@{ 
    int domsize_;
    PredicateDisjunctive **disjuncts;
    //@}  
  
    /**@name Utils*/
    //@{ 
    inline bool operator<( VarSelectorOSP_Domain& x ) const { 
      return (domsize_ < x.domsize_) ;
    }
    inline void operator=( VarSelectorOSP_Domain& x ) { 
      domsize_ = x.domsize_; 
    }
    inline void operator=( VariableInt    *x ) 
    { 
      int idx = x->id;
      if(disjuncts[idx]) {
	domsize_ = disjuncts[idx]->domsize();
      } else {
	domsize_ = x->domsize();
      }

      //domsize_ = disjuncts[x->id]->domsize();
    }
    //@}  

    std::ostream& print(std::ostream& os) const {
      os << "-" ;
      return os;
    }
  };


  class VarSelectorOSP_DomainType
  {
  public: 

    /**@name Constructors*/
    //@{
    VarSelectorOSP_DomainType() {domsize_ = 0;}
    //@}

    double value() { return (double)domsize_; }

    /**@name Parameters*/
    //@{ 
    int type_;
    int domsize_;
    PredicateDisjunctive **disjuncts;
    //@}  
  
    /**@name Utils*/
    //@{ 
    inline bool operator<( VarSelectorOSP_DomainType& x ) const { 
      return (type_ < x.type_ || (type_ == x.type_ && (domsize_ < x.domsize_)));
    }
    inline void operator=( VarSelectorOSP_DomainType& x ) { 
      domsize_ = x.domsize_; 
    }
    inline void operator=( VariableInt    *x ) 
    { 
      
      int idx = x->id;
      if(disjuncts[idx]) {
	type_ = 1;
	domsize_ = disjuncts[idx]->domsize();
      } else {
	domsize_ = x->domsize();
	type_ = 2*(domsize_ > 2);
      }

      //domsize_ = disjuncts[x->id]->domsize();
    }
    //@}  

    std::ostream& print(std::ostream& os) const {
      os << "-" ;
      return os;
    }
  };



  class VarSelectorOSP_DoBoolWeight 
  {
  public: 

    /**@name Constructors*/
    //@{
    VarSelectorOSP_DoBoolWeight() {domsize_ = 0; weight_ = 0;}
    //@}

    double value() { return ((double)domsize_ / (double)weight_); }

    /**@name Parameters*/
    //@{ 
    int domsize_;
    int weight_;
    PredicateDisjunctive **disjuncts;
    //@}  
  
    /**@name Utils*/
    //@{ 
    inline bool operator<( VarSelectorOSP_DoBoolWeight& x ) const { 
      return (domsize_ * x.weight_ < x.domsize_ * weight_) ;
    }
    inline void operator=( VarSelectorOSP_DoBoolWeight& x ) { 
      weight_ = x.weight_; domsize_ = x.domsize_; 
    }
    inline void operator=( VariableInt    *x ) 
    { 

      int idx = x->id;
      if(disjuncts[idx]) {
	domsize_ = disjuncts[idx]->domsize();
      } else {
	domsize_ = x->domsize();
      }
      weight_ = x->weight;

    }
    //@}  

    std::ostream& print(std::ostream& os) const {
      os << "-" ;
      return os;
    }
  };


  class VarSelectorOSP_DoBoolWeightType 
  {
  public: 

    /**@name Constructors*/
    //@{
    VarSelectorOSP_DoBoolWeightType() {domsize_ = 0; weight_ = 0;}
    //@}

    double value() { return ((double)weight_); }

    /**@name Parameters*/
    //@{ 
    int type_;
    int domsize_;
    int weight_;
    PredicateDisjunctive **disjuncts;
    //@}  
  
    /**@name Utils*/
    //@{ 
    inline bool operator<( VarSelectorOSP_DoBoolWeightType& x ) const { 
      return (type_ < x.type_ || (type_ == x.type_ && (domsize_ * x.weight_ < x.domsize_ * weight_))) ;
    }
    inline void operator=( VarSelectorOSP_DoBoolWeightType& x ) { 
      weight_ = x.weight_; domsize_ = x.domsize_; 
    }
    inline void operator=( VariableInt    *x ) 
    { 

      //x->print(std::cout);

      int idx = x->id;
      if(disjuncts[idx]) {
	type_ = 1;
	domsize_ = disjuncts[idx]->domsize();
      } else {
	domsize_ = x->domsize();
	type_ = 2*(domsize_ > 2);
      }
      weight_ = x->weight;

      //std::cout << " dom = " << domsize_ << " weight = " << weight_ << std::endl;

    }
    //@}  

    std::ostream& print(std::ostream& os) const {
      os << "-" ;
      return os;
    }
  };


  class VarSelectorOSP_DoTaskWeight 
  {
  public: 
    /**@name Constructors*/
    //@{
    VarSelectorOSP_DoTaskWeight() {domsize_ = 0; task_weight_ = 0;}
    //@}

    double value() { return ((double)task_weight_); }

    /**@name Parameters*/
    //@{ 
    int domsize_;
    int task_weight_;
    PredicateDisjunctive **disjuncts;
    //@}  
  
    /**@name Utils*/
    //@{ 
    inline bool operator<( VarSelectorOSP_DoTaskWeight& x ) const { 
      return (domsize_ * x.task_weight_ < x.domsize_ * task_weight_) ;
    }
    inline void operator=( VarSelectorOSP_DoTaskWeight& x ) { 
      task_weight_ = x.task_weight_; domsize_ = x.domsize_; 
    }
    inline void operator=( VariableInt    *x ) 
    { 
      PredicateDisjunctive* d = disjuncts[x->id];
      if(d) {
	domsize_ = d->domsize();
	task_weight_ = (d->scope[0]->weight + d->scope[1]->weight);
      } else {
	domsize_ = x->domsize();
	task_weight_ = x->weight;
      }
    }
      
    //@}  

    std::ostream& print(std::ostream& os) const {
      os << "-" ;
      return os;
    }
  };

  class VarSelectorOSP_DoTaskWeightType
  {
  public: 
    /**@name Constructors*/
    //@{
    VarSelectorOSP_DoTaskWeightType() {domsize_ = 0; task_weight_ = 0;}
    //@}

    double value() { return ((double)task_weight_); }

    /**@name Parameters*/
    //@{ 
    int type_;
    int domsize_;
    int task_weight_;
    PredicateDisjunctive **disjuncts;
    //@}  
  
    /**@name Utils*/
    //@{ 
    inline bool operator<( VarSelectorOSP_DoTaskWeightType& x ) const { 
      return (type_ < x.type_ || (type_ == x.type_ && (domsize_ * x.task_weight_ < x.domsize_ * task_weight_)));
    }
    inline void operator=( VarSelectorOSP_DoTaskWeightType& x ) { 
      task_weight_ = x.task_weight_; domsize_ = x.domsize_; 
    }
    inline void operator=( VariableInt    *x ) 
    { 
      PredicateDisjunctive* d = disjuncts[x->id];
      if(d) {
	type_ = 1;
	domsize_ = d->domsize();
	task_weight_ = (d->scope[0]->weight + d->scope[1]->weight);
	//std::cout << domsize_ << " / " << task_weight_ << std::endl;
      } else {
	domsize_ = x->domsize();
	task_weight_ = x->weight;
	type_ = 2*(domsize_ > 2);
      }
    }
      
    //@}  

    std::ostream& print(std::ostream& os) const {
      os << "-" ;
      return os;
    }
  };

  class VarSelectorOSP_DomainWeight 
  {
  public: 
    /**@name Constructors*/
    //@{
    VarSelectorOSP_DomainWeight() {domsize_ = 0; weight_ = 0;}
    //@}

    double value() { return ((double)weight_); }

    /**@name Parameters*/
    //@{ 
    int domsize_;
    //int intersize_;
    //int prunesize_;
    int weight_;
    PredicateDisjunctive **disjuncts;
    //@}  
  
    /**@name Utils*/
    //@{ 
    inline bool operator<( VarSelectorOSP_DomainWeight& x ) const { 
      
      //std::cout << "osp d/bw" << std::endl;

      return (domsize_ < x.domsize_ || ( domsize_ == x.domsize_ && x.weight_ < weight_ )) ;

//       int self = ((domsize_ * domsize_) / intersize_);
//       int other = ((x.domsize_ * x.domsize_) / x.intersize_);
//       return( self < other || (self == other && x.weight_ < weight_ ) ); 
    }
    inline void operator=( VarSelectorOSP_DomainWeight& x ) { 
      weight_ = x.weight_; domsize_ = x.domsize_; 
    }
    inline void operator=( VariableInt    *x ) 
    { 
      PredicateDisjunctive* d = disjuncts[x->id];
      //domsize_ = (d->domsize() - (std::min(d->XltY(), d->YltX())));

      
      //intersize_ = 1+(d->intersection());
      
      domsize_ = (d->domsize());
      weight_ = x->weight;
      //weight_ = (d->scope[0]->weight + d->scope[1]->weight);
    }
    //@}  

    std::ostream& print(std::ostream& os) const {
      os << "-" ;
      return os;
    }
  };


  class VarSelectorOSP_DoBoolTaskWeight 
  {
  public: 

    /**@name Constructors*/
    //@{
    VarSelectorOSP_DoBoolTaskWeight() {domsize_ = 0; weight_ = 0; task_weight_ = 0;}
    //@}

    double value() { return ((double)(weight_ + task_weight_)); }


    /**@name Parameters*/
    //@{ 
    int domsize_;
    int weight_;
    int task_weight_;
    PredicateDisjunctive **disjuncts;
    //@}  
  
    /**@name Utils*/
    //@{ 
    inline bool operator<( VarSelectorOSP_DoBoolTaskWeight& x ) const { 
      return (domsize_ * (x.weight_ + x.task_weight_) < x.domsize_ * (weight_ + task_weight_)) ;
    }
    inline void operator=( VarSelectorOSP_DoBoolTaskWeight& x ) { 
      weight_ = x.weight_; task_weight_ = x.task_weight_; domsize_ = x.domsize_; 
    }
    inline void operator=( VariableInt    *x ) 
    { 
      PredicateDisjunctive* d = disjuncts[x->id];
      if(d) {
	domsize_ = d->domsize();
	task_weight_ = (d->scope[0]->weight + d->scope[1]->weight);
	//std::cout << domsize_ << " / " << task_weight_ << std::endl;
      } else {
	domsize_ = x->domsize();
	task_weight_ = x->weight;
      }
      weight_ = x->weight;
    }
    //@}  public: 

    std::ostream& print(std::ostream& os) const {
      os << "-" ;
      return os;
    }
  };


  class VarSelectorOSP_DoBoolTaskWeightType
  {
  public: 

    /**@name Constructors*/
    //@{
    VarSelectorOSP_DoBoolTaskWeightType() {domsize_ = 0; weight_ = 0; task_weight_ = 0;}
    //@}

    double value() { return ((double)(weight_ + task_weight_)); }


    /**@name Parameters*/
    //@{ 
    int type_;
    int domsize_;
    int weight_;
    int task_weight_;
    PredicateDisjunctive **disjuncts;
    //@}  
  
    /**@name Utils*/
    //@{ 
    inline bool operator<( VarSelectorOSP_DoBoolTaskWeightType& x ) const { 
      return (type_ < x.type_ || (type_ == x.type_ && (domsize_ * (x.weight_ + x.task_weight_) < x.domsize_ * (weight_ + task_weight_))));
    }
    inline void operator=( VarSelectorOSP_DoBoolTaskWeightType& x ) { 
      weight_ = x.weight_; task_weight_ = x.task_weight_; domsize_ = x.domsize_; 
    }
    inline void operator=( VariableInt    *x ) 
    { 
      PredicateDisjunctive* d = disjuncts[x->id];
      if(d) {
	type_ = 1;
	domsize_ = d->domsize();
	task_weight_ = (d->scope[0]->weight + d->scope[1]->weight);
	//std::cout << domsize_ << " / " << task_weight_ << std::endl;
      } else {
	domsize_ = x->domsize();
	task_weight_ = x->weight;
	type_ = 2*(domsize_ > 2);
      }
      weight_ = x->weight;
    }
    //@}  public: 

    std::ostream& print(std::ostream& os) const {
      os << "-" ;
      return os;
    }
  };



  class VarSelectorOSP_DoWeakWeight 
  {
  public: 

    /**@name Constructors*/
    //@{
    VarSelectorOSP_DoWeakWeight() {domsize_ = 0; weight_ = 0; task_weight_ = 0; weak_ = 0;}
    //@}

    double value() { return ((double)(weight_ + weak_)); }


    /**@name Parameters*/
    //@{ 
    int weak_;
    int domsize_;
    int weight_;
    int task_weight_;
    PredicateDisjunctive **disjuncts;
    //@}  
  
    /**@name Utils*/
    //@{ 
    inline bool operator<( VarSelectorOSP_DoWeakWeight& x ) const { 
      return (domsize_ * (x.weight_ + weak_) < x.domsize_ * (weight_ + weak_)) ;
      //return (domsize_ < x.domsize_) ;
    }
    inline void operator=( VarSelectorOSP_DoWeakWeight& x ) { 
      weight_ = x.weight_; task_weight_ = x.task_weight_; domsize_ = x.domsize_; 
    }
    inline void operator=( VariableInt    *x ) 
    { 
      PredicateDisjunctive* d = disjuncts[x->id];
      domsize_ = d->domsize();
      task_weight_ = (d->scope[0]->weight + d->scope[1]->weight);
      weight_ = x->weight;
    }
    //@}  public: 

    std::ostream& print(std::ostream& os) const {
      os << "-" ;
      return os;
    }
  };


  class VarSelectorOSPSAT
  {
  public: 

    /**@name Constructors*/
    //@{
    VarSelectorOSPSAT() {activity_ = 0.0; domsize_ = 0, weight_ = 0; }
    //@}

    double value() { return (double)((activity_*20+1) + weight_); }

    /**@name Parameters*/
    //@{ 
    int domsize_;
    int weight_;
    double activity_;
    double *activity;
    PredicateDisjunct **disjuncts;
    //@}  
  
    /**@name Utils*/
    //@{ 
    inline bool operator<( VarSelectorOSPSAT& x ) const { 
      //return (domsize_ * x.weight_ < x.domsize_ * weight_) ;
      //return (domsize_ * x.activity_ < x.domsize_ * activity_) ;

      //     if( x.weight_ < x.activity_ )
      //       std::cout << 11 << std::endl;

      //     if( weight_ < activity_ )
      //       std::cout << 22 << std::endl;

      //     //std::cout << weight_ << " " << activity_ << std::endl;


      //     if(x.activity_ > x.weight_)
      //       std::cout << x.activity_ << std::endl;

      int domact_a = domsize_ * ((int)(x.activity_*20+1) + x.weight_);
      int domact_b = x.domsize_ * ((int)(activity_*20+1) + weight_);

      //assert( domact_a >= 0 );
      //assert( domact_b >= 0 );

      return domact_a < domact_b;
    

      //return (domsize_ * (x.weight_ + x.activity_) < x.domsize_ * (weight_ + activity_)) ;
    }
    inline void operator=( VarSelectorOSPSAT& x ) {
      weight_ = x.weight_; 
      activity_ = x.activity_; 
      domsize_ = x.domsize_; 
    }
    inline void operator=( VariableInt    *x ) 
    { 
      int idx = x->id;
      domsize_ = disjuncts[idx]->domsize();
      weight_ = x->weight;
      ++idx;
      activity_ = (activity[idx] + activity[-idx]);

      //std::cout << activity_ << " " << weight_ << std::endl;
    }
    //@}  

    std::ostream& print(std::ostream& os) const {
      os << "-" ;
      return os;
    }
  };


  class VarSelectorFPP 
  {
  public: 

    /**@name Constructors*/
    //@{
    VarSelectorFPP() {domsize_ = 0; vweight_ = 0; tweight_ = 0; cweight_ = 0; }
    //@}

    double value() { return (double)(vweight_); }

    /**@name Parameters*/
    //@{ 
    int domsize_;
    int vweight_;
    int tweight_;
    int cweight_;
    int *isInterval;
    PredicateLess **precs;
    PredicateUpperBound **inters;
    //@}  
  
    /**@name Utils*/
    //@{ 
    inline bool operator<( VarSelectorFPP& x ) const { 
      //     //return (domsize_ * x.weight_ < x.domsize_ * weight_) ;

      //     //int r1 = domsize_ * x.vweight_;
      //     //int r2 = x.domsize_ * vweight_;

      //     int r1 = domsize_ * x.vweight_ + x.tweight_;
      //     int r2 = x.domsize_ * vweight_ + tweight_;

      //     return (r1 < r2);
      //     //bool ret = (r1 < r2 || (r1 == r2 && vweight_ > x.vweight_));

      // //     if( ret ) {
      // //       std::cout << x.domsize_ << " " << x.vweight_ << " " << x.tweight_ << " " << x.cweight_ ;
      // //       std::cout << " -> ";
      // //       std::cout << domsize_ << " " << vweight_ << " " << tweight_ << " " << cweight_ << std::endl;
      // //     }

      //     //return ret;
    
      //return ((domsize_ < x.domsize_) || (domsize_ == x.domsize_ && vweight_ > x.vweight_)) ;
      return ((vweight_ > x.vweight_) || 
	      (vweight_ == x.vweight_ && domsize_ < x.domsize_));
      //  ||
      // 	    (vweight_ == x.vweight_ && domsize_ == x.domsize_ && tweight_ > x.tweight_)) ;

    }
    inline void operator=( VarSelectorFPP& x ) { 
      tweight_ = x.tweight_; cweight_ = x.cweight_; 
      vweight_ = x.vweight_; domsize_ = x.domsize_; 
    }
    inline void operator=( VariableInt    *x ) 
    { 
      //     int i, idx = x->id;
      //     idx -= (idx%4);
      //     for(i=0; i<4; ++i)
      //       domsize_ = prec[idx+i]->sprecscope[0]->domsize();

      //     precs[x->id]->print( std::cout );
      //     std::cout << std::endl;

      if( isInterval[x->id] )
	exit(0);

      VariableInt **squares = precs[x->id]->scope;
      domsize_ = (squares[0]->domsize() + squares[1]->domsize());
      vweight_ = x->weight;
      tweight_ = (squares[0]->weight + squares[1]->weight);
      cweight_ = precs[x->id]->weight;
      //    std::cout << domsize_ << "/" << weight_ << std::endl;
    }
    //@}  

    std::ostream& print(std::ostream& os) const {
      os << "-" ;
      return os;
    }
  };


  class ValSelector
  {
  public:

    /**@name Parameters*/
    //@{   
    VariableInt *_X;
    //@}

    /**@name Constructors*/
    //@{
    virtual ~ValSelector() {}
    ValSelector( VariableInt *x ) : _X(x) {}
    ValSelector( ) : _X(NULL) {}
    //@}

    /**@name Utils*/
    //@{ 
    //virtual double value() { return 1.0; }
    virtual void make(int& t, int& v) = 0;
    virtual int getBest() = 0;
    virtual void getOrder( int* dtv, int* vtd ) {

      DomainIterator *valit = _X->begin();
      int i=0;
      do {
	dtv[i] = *valit;
	vtd[*valit] = i;
	++i;
      } while( valit->next() ); 

    }
    virtual void left() = 0;
    virtual void right() = 0;
    virtual void reverse_left() = 0;
    virtual void reverse_right() = 0;
    virtual void postCut( const int p ) = 0;
    //@}  

    /**@name Miscellanous*/
    //@{
    virtual void printLeft(std::ostream& o) const {}
    virtual void printRight(std::ostream& o) const {}
    //@}
  };


  class ValSelectorMin : public ValSelector
  {
  public:

    /**@name Parameters*/
    //@{
    int val;
    Vector< int > decision;
    //@}

    /**@name Constructors*/
    //@{
    ValSelectorMin( VariableInt *x ) : ValSelector(x) { }
    //@}

    /**@name Utils*/
    //@{ 
    void make(int& t, int& v) { v=_X->min(); t=Decision::ASSIGNMENT;}
    inline int getBest() { return (val = _X->min()); }
    inline void left() {
      val = _X->min();
      _X->setDomain( val ); 
    }
    inline void right() { 
      //val = _X->min();
      _X->remove( val ); 
    }
    inline void reverse_left() {
      val = _X->min();
      _X->remove( val ); 
    }
    inline void reverse_right() { 
      val = _X->min();
      _X->setDomain( val ); 
    }
    void postCut( const int p ) {
      val = p;
    }
    //@}  

    /**@name Miscellanous*/
    //@{
    virtual void printLeft(std::ostream& o) const { 
      std::cout << " == " << val;
    }
    virtual void printRight(std::ostream& o) const {
      std::cout << " =/= " << val;
    }
    //@}
  };


  class ValSelectorFirst : public ValSelector
  {
  public:

    /**@name Parameters*/
    //@{
    int val; 
    Vector< int > decision;
    //@}   
  
    /**@name Constructors*/
    //@{
    ValSelectorFirst( VariableInt *x ) : ValSelector(x) { decision.init(0,8); }
    //@}

    /**@name Utils*/
    //@{ 
    void make(int& t, int& v) { v=_X->first(); t=Decision::ASSIGNMENT;}
    inline int getBest() { return _X->first(); }
    inline void left() { 
      val = _X->first(); 
      _X->setDomain( val ); 
    }
    inline void right() { 
      _X->remove( val ); 
    }
    inline void reverse_left() {
      val = _X->first();
      decision.push( val );
      _X->remove( val ); 
    }
    inline void reverse_right() { 
      val = decision.pop();
      _X->setDomain( val );
    }
    void postCut( const int p ) {
      val = p;
    }
    //@}  

    /**@name Miscellanous*/
    //@{
    virtual void printLeft(std::ostream& o) const { 
      std::cout << " == " << val;
    }
    virtual void printRight(std::ostream& o) const {
      std::cout << " =/= " << val;
    }
    //@}
  };

  class ValSelectorMax : public ValSelector
  {
  public:

    /**@name Parameters*/
    //@{
    int val; 
    //@}   

    /**@name Constructors*/
    //@{
    ValSelectorMax( VariableInt *x ) : ValSelector(x) {}
    //@}

    /**@name Utils*/
    //@{ 
    void make(int& t, int& v) { v=_X->max(); t=Decision::ASSIGNMENT;}
    inline int getBest() { return (val = _X->max()); }
    inline void left() { 
      val = _X->max(); 
      _X->setDomain( val ); 
    }
    inline void right() { 
      _X->remove( val ); 
    }
    inline void reverse_left() {
      val = _X->max(); 
      _X->remove( val ); 
    }
    inline void reverse_right() { 
      val = _X->max(); 
      _X->setDomain( val );
    }
    void postCut( const int p ) { val = p; }
    //@}  

    /**@name Miscellanous*/
    //@{
    virtual void printLeft(std::ostream& o) const { 
      std::cout << " == " << val;
    }
    virtual void printRight(std::ostream& o) const {
      std::cout << " =/= " << val;
    }
    //@}
  };

  class ValSelectorRand : public ValSelector 
  { 
  public: 

    /**@name Constructors*/
    //@{
    ValSelectorRand( VariableInt *x ) : ValSelector(x) { decision.init(0,8); } 
    //@}   

    /**@name Parameters*/
    //@{
    int val; 
    Vector< int > decision;
    //@}   

    /**@name Utils*/ 
    //@{  
    void make(int& t, int& v) { v=_X->random(); t=Decision::ASSIGNMENT;}
    inline int getBest() { return _X->max(); }
    inline void left() { 
      val = _X->random(); 
      _X->setDomain( val ); 
    } 
    inline void right() { 
      _X->remove( val ); 
    } 
    inline void reverse_left() {
      val = _X->random();
      decision.push( val );
      _X->remove( val ); 
    }
    inline void reverse_right() { 
      val = decision.pop();
      _X->setDomain( val );
    }

    /**@name Miscellanous*/
    //@{
    virtual void printLeft(std::ostream& o) const { 
      std::cout << " == " << val;
    }
    virtual void printRight(std::ostream& o) const {
      std::cout << " =/= " << val;
    }
    //@}

    void postCut( const int p ) { val = p; }
    //@}   
  }; 

  class ValSelectorGuided : public ValSelector 
  { 
  public: 

    /**@name Constructors*/
    //@{
    ValSelectorGuided( VariableInt *x, int i, const int pb ) : ValSelector(x) { 
      ideal = i; 
      planB = pb;
      decision.init(0,8);
    } 
    //@}   

    /**@name Parameters*/
    //@{
    Vector< int > decision;
    int ideal;
    int val;
    int planB;
    /// 0: min
    /// 1: nearest
    /// 2: nearest split
    /// 3: nearest bound
    //@}   

    /**@name Utils*/ 
    //@{  
    void make(int& t, int& v) {
     
      if(_X->getType() == VariableInt::RANGE) {
	if(_X->max()-ideal < ideal-_X->min())
	  v = _X->max();
	else
	  v = _X->min();
	t=Decision::ASSIGNMENT;
      } else {
	if(_X->contain(ideal)) {
	  v = ideal;
	  t=Decision::ASSIGNMENT;
	} else {
	  if(planB == 0) {
	    v = _X->min();
	    t=Decision::ASSIGNMENT;
	  } else if(planB == 1) {
	    if(_X->min() > ideal) {
	      v = _X->min();
	    } else if(_X->max() < ideal) {
	      v = _X->max();
	    } else {
	      int i=1, dir=1;
	      while(true) {
		v = ideal + (dir*i);
		if(_X->contain(v)) break;
		if(dir < 0) ++i;
		dir*=-1;
	      }
	    }
	    t=Decision::ASSIGNMENT;
	  } else if(planB == 2) {
	    v = ((_X->max() + _X->min()) >> 1); 
	    if(ideal <= val) t=Decision::UPPERBOUND;
	    else t=Decision::LOWERBOUND;
	  } else if(planB == 3) {
	    int lb = _X->min();
	    int ub = _X->max();
	    if(ub < ideal) v = ub;
	    else if(lb > ideal) v = lb;
	    else if(ideal-lb < ub-ideal) v = lb;
	    else v = ub;
	    t=Decision::ASSIGNMENT;
	  }
	}
      }
    }
    inline int getBest() { return ideal; }
    inline void setSecondBest() {
      int i=1, dir=1;
      while(true) {
	val = ideal + (dir*i);
	if(_X->contain(val)) break;
	if(dir < 0) ++i;
	dir*=-1;
      }
      decision.push(val);

      //_X->print(std::cout);
      //std::cout << " <- " << val << " (" << ideal << ")" << std::endl;

      _X->setDomain(val);
    }
    inline void bestSplit() {
      val = ((_X->max() + _X->min()) >> 1); 
      decision.push(val);
      if(ideal < val) _X->setMax( val ); 
      else _X->setMin( val+1 ); 
    }
    inline void worstSplit() {
      val = ((_X->max() + _X->min()) >> 1); 
      if(ideal < val) _X->setMin( val+1 ); 
      else _X->setMax( val ); 
    }
    inline void left() { 
      if(_X->contain(ideal)) {
	decision.push(ideal);
	_X->setDomain(ideal);
      } else {
	if(planB == 0) {
	  val = _X->min();
	  decision.push(val);
	  _X->setDomain(val);
	} else if(planB == 1) {
	  setSecondBest();
	} else if(planB == 2) {
	  bestSplit();
	}
      }
    } 
    inline void right() {
      val = decision.pop();
      if(val != ideal && planB == 2) worstSplit();
      else _X->remove( val ); 
    } 
    inline void reverse_left() {
      val = (_X->contain(ideal) ? ideal : _X->min());
      decision.push( val );
      _X->remove( val ); 
    }
    inline void reverse_right() { 
      val = decision.pop();
      _X->setDomain( val );
    }
    void postCut( const int p ) { val = p; }
    //@}   
  };

  class ValSelectorBoolean : public ValSelector 
  { 
  public: 

    /**@name Constructors*/
    //@{
    ValSelectorBoolean( VariableInt *x, VariableInt *b ) : ValSelector(x) { 
      earlybool=b; 
      decision.init(0,8);
    } 
    //@}   

    /**@name Parameters*/
    //@{
    Vector< int > decision;
    VariableInt* earlybool;
    int val;
    //@}   

    /**@name Utils*/ 
    //@{
    void make(int& t, int& v) { 
      if(earlybool->min()) v = _X->max();
      else v = _X->max();
      t=Decision::ASSIGNMENT;
    }
    inline int getBest() { if(earlybool->min()) return _X->max(); else return _X->min(); }
    inline void left() { 
      val = getBest();
      _X->setDomain(val);
    } 
    inline void right() { 
      _X->remove( val ); 
    } 
    inline void reverse_left() {
      val = getBest();
      decision.push( val );
      _X->remove( val ); 
    }
    inline void reverse_right() { 
      val = decision.pop();
      _X->setDomain( val );
    }
    void postCut( const int p ) { val = p; }
    //@}   
  };


  class ValSelectorContention : public ValSelector 
  { 
  public: 

    /**@name Constructors*/
    //@{
    ValSelectorContention( VariableInt *x, Contention *c, const int i ) 
      : ValSelector(x) { 
      decision.init(0,8);
      probability_base = c;
      id = i;
    } 
    //@}   

    /**@name Parameters*/
    //@{
    Vector< int > decision;
    Contention *probability_base;
    int id;
    int val;
    //@}   

    /**@name Utils*/ 
    //@{  
    void make(int& t, int& v) { 
      v = probability_base->update(id);
      t=Decision::ASSIGNMENT;
    }
    inline int getBest() { return val; }
    inline void left() { 

//       DomainIterator *valit = _X->begin();
//       do {
	
// 	if(proba[*valit] > proba[val]) val = *valit;
//       } while(valit->next());
      //val = probability_base->update(id);
//       probability_base->reset();
//       for(int i=0; i<1; ++i) {
// 	probability_base->print(std::cout, id);
// 	probability_base->compute_proba();
//       }
//       probability_base->print(std::cout, id);

      //val = probability_base->get_best(id);
      //std::cout << "\t-> " << val << std::endl;
      val = probability_base->update(id);
      //std::cout << "\t-> " << val << std::endl << std::endl;
      
      _X->setDomain( val ); 
    } 
    inline void right() { 
      _X->remove( val ); 
    } 
    inline void reverse_left() {
      val = probability_base->update(id);
      decision.push( val );
      _X->remove( val ); 
    }
    inline void reverse_right() { 
      val = decision.pop();
      _X->setDomain( val );
    }
    void postCut( const int p ) { val = p; }
    //@}   
  };

  class ValSelectorRandGuided : public ValSelector 
  { 
  public: 

    /**@name Constructors*/
    //@{
    ValSelectorRandGuided( VariableInt *x, int i, int p, int r ) : ValSelector(x) 
    {ideal=i; proba=p; range=r; decision.init(0,8);} 
    //@}   

    /**@name Parameters*/
    //@{
    int ideal;
    int proba;
    int range;
    int val;
    Vector< int > decision;
    //@}   

    /**@name Utils*/ 
    //@{  
    void make(int& t, int& v) { 
      if(_X->contain(ideal) && (randint(range) < proba)) v=ideal;
      else v = (randint(2) ? _X->min() : _X->max());
      t=Decision::ASSIGNMENT;
    }
    inline int getBest() { return ideal; }
    inline void left() { 
      if(_X->contain(ideal) && (randint(range) < proba)) val=ideal;
      else val = (randint(2) ? _X->min() : _X->max());
      _X->setDomain( val );
    } 
    inline void right() { 
      _X->setDomain( val ); 
    } 
    inline void reverse_left() {
      if(_X->contain(ideal) && (randint(range) < proba)) val=ideal;
      else val = (randint(2) ? _X->min() : _X->max());
      decision.push( val );
      _X->remove( val ); 
    }
    inline void reverse_right() { 
      val = decision.pop();
      _X->setDomain( val );
    }
    void postCut( const int p ) { val = p; }
    //@}   
  };  


  class ValSelectorRandMinMax : public ValSelector 
  { 
  public: 

    /**@name Constructors*/
    //@{
    ValSelectorRandMinMax( VariableInt *x ) : ValSelector(x) { decision.init(0,8); } 
    //@}   

    /**@name Parameters*/
    //@{
    int val; 
    Vector< int > decision;
    //@}   

    /**@name Utils*/ 
    //@{  
    void make(int& t, int& v) { 
      v = (randint(2) ? _X->min() : _X->max());
      t=Decision::ASSIGNMENT;
    }
    inline int getBest() { return _X->first(); }
    inline void left() { 
      val = (randint(2) ? _X->min() : _X->max()); 
      _X->setDomain( val ); 
    }
    inline void right() { 
      _X->remove( val ); 
    }
    inline void reverse_left() {
      val = (randint(2) ? _X->min() : _X->max()); 
      decision.push( val );
      _X->remove( val ); 
    }
    inline void reverse_right() {
      val = decision.pop(); 
      _X->setDomain( val );
    }
    void postCut( const int p ) { val = p; }
    //@}  
   
    /**@name Miscellanous*/
    //@{
    virtual void printLeft(std::ostream& o) const { 
      std::cout << " == " << val;
    }
    virtual void printRight(std::ostream& o) const {
      std::cout << " =/= " << val;
    }
    //@}
   
  }; 

  class ValSelectorSplit : public ValSelector
  {
  public:

    int val; 
    

    /**@name Constructors*/
    //@{
    ValSelectorSplit( VariableInt *x ) : ValSelector(x) {}
    //@}   
  
    /**@name Utils*/
    //@{ 
    void make(int& t, int& v) { 
      v = ((_X->max() + _X->min()) >> 1); 
      t=Decision::UPPERBOUND;
    }
    inline int getBest() { return (val = _X->min()); }
    inline void left() { 
      val = ((_X->max() + _X->min()) >> 1); 
      _X->setMax( val ); 
    }
    inline void right() { 
      val = ((_X->max() + _X->min()) >> 1); 
      _X->setMin( val+1 ); 
    }
    inline void reverse_left() {
      val = ((_X->max() + _X->min()) >> 1); 
      _X->setMin( val+1 ); 
    }
    inline void reverse_right() {
      val = ((_X->max() + _X->min()) >> 1); 
      _X->setMax( val ); 
    }
    void postCut( const int p ) { }
    //@}  

    /**@name Miscellanous*/
    //@{
    virtual void printLeft(std::ostream& o) const { 
      std::cout << " <= " << val;
    }
    virtual void printRight(std::ostream& o) const {
      std::cout << " > " << val;
    }
    //@}
  };


  class ValSelectorRandomSplit : public ValSelector
  {
  public:

    Vector< int > decision;

    /**@name Constructors*/
    //@{
    ValSelectorRandomSplit( VariableInt *x ) : ValSelector(x) { decision.init(0,8); }
    //@}   
  
    /**@name Utils*/
    //@{ 
    void make(int& t, int& v) { 
      v = _X->min();
      v += randint(_X->max() - v);
      if(randint(2)) t=Decision::UPPERBOUND;
      else t=Decision::LOWERBOUND;
    }
    inline int getBest() { return (_X->min()); }
    inline void left() { 
      int val = _X->min();
      int dir = randint(2);
      val += randint(_X->max() - val);
      if(dir) _X->setMax( val ); 
      else _X->setMin( val+1 );
      decision.push((val << 1) + dir);
    }
    inline void right() { 
      int val = decision.pop();
      int dir = (val & 1);
      val >>= 1;
      if(dir) _X->setMin( val+1 ); else _X->setMax( val ); 
    }
    inline void reverse_left() {
      int val = _X->min();
      int dir = randint(2);
      val += randint(_X->max() - val);
      if(dir) _X->setMin( val+1 ); 
      else _X->setMax( val );
      decision.push((val << 1) + dir);
    }
    inline void reverse_right() { 
      int val = decision.pop();
      int dir = (val & 1);
      val >>= 1;
      if(dir) _X->setMax( val ); else _X->setMin( val ); 
    }
    void postCut( const int p ) { }
    //@}  

   /**@name Miscellanous*/
    //@{
    virtual void printLeft(std::ostream& o) const { 

      int val = decision.back();
      int dir = (val & 1);
      val >>= 1;

      if(dir) {
	std::cout << " <= " << val ;
      } else {
	std::cout << " >= " << (val+1) ;
      }
    }
    virtual void printRight(std::ostream& o) const {

      int val = decision[decision.size];
      int dir = (val & 1);
      val >>= 1;

      if(dir) {
	o << " >= " << (val+1) ;
      } else {
	o << " <= " << (val) ;
      }
    }
    //@}
  };


  class ValSelectorRandGuidedSplit : public ValSelector
  {
  public:
    
    Vector< int > decision;
    int value;
    int proba;
    
    /**@name Constructors*/
    //@{
    ValSelectorRandGuidedSplit( ) : ValSelector() { value=0; proba=1; }
    ValSelectorRandGuidedSplit( VariableInt *x, int v, int p ) : ValSelector(x) { 

//       std::cout << "================== ";
//       x->print(std::cout);
//       std::cout << " <> " << v << " (" << p << ")" << std::endl;

      decision.init(0,8); value=v; proba=p; 
    }
    //@}   
  
    /**@name Utils*/
    //@{ 
    void make(int& t, int& v) { 
      v = _X->min();
      v += randint(_X->max() - v);
      if((randint(1000) < proba)^(value>v)) t=Decision::UPPERBOUND;
      else t=Decision::LOWERBOUND;
    }
    inline int getBest() { return (value); }
    inline void left() { 
      int val = _X->min();
      int dir = (randint(1000) < proba);

      val += randint(_X->max() - val);
      dir ^= (value > val);

      if(dir) _X->setMax( val ); 
      else _X->setMin( val+1 );
      
      decision.push((val << 1) + dir);
    }
    inline void right() { 
      int val = decision.pop();
      int dir = (val & 1);
      val >>= 1;
      if(dir) _X->setMin( val+1 ); else _X->setMax( val ); 
    }
    inline void reverse_left() {
      int val = _X->min();
      int dir = (randint(1000) < proba);

      val += randint(_X->max() - val);
      dir ^= (value > val);

      if(dir) _X->setMin( val+1 ); 
      else _X->setMax( val );
      
      decision.push((val << 1) + dir);
    }
    inline void reverse_right() { 
      int val = decision.pop();
      int dir = (val & 1);
      val >>= 1;
      if(dir) _X->setMax( val ); else _X->setMin( val+1 ); 
    }
    void postCut( const int p ) { }
    //@}  



    /**@name Miscellanous*/
    //@{
    virtual void printLeft(std::ostream& o) const { 

      int val = decision.back();
      int dir = (val & 1);
      val >>= 1;

      if(dir) {
	std::cout << " <= " << val << " (" << value << ")";
      } else {
	std::cout << " >= " << (val+1) << " (" << value << ")";
      }
    }
    virtual void printRight(std::ostream& o) const {
      //_X->printshort(o);

      int val = decision[decision.size];
      int dir = (val & 1);
      val >>= 1;

      if(dir) {
	o << " >= " << (val+1) << " (" << value << ")";
      } else {
	o << " <= " << (val) << " (" << value << ")";
      }
    }
    //@}
  };



//   class ValSelectorRandMultiGuidedSplit : public ValSelector
//   {
//   public:
    
//     Vector< int > decision;
//     int value;
//     int proba;
    
//     /**@name Constructors*/
//     //@{
//     ValSelectorRandGuidedSplit( VariableInt *x, int v, int p ) : ValSelector(x) { 
//       decision.init(0,8); value=v; proba=p; 
//     }
//     //@}   
  
//     /**@name Utils*/
//     //@{ 
//     inline int getBest() { return (value); }
//     inline void left() { 
//       int val = _X->min();
//       int dir = (randint(1000) < proba);

//       val += randint(_X->max() - val);
//       dir ^= (value > val);

//       if(dir) _X->setMax( val ); 
//       else _X->setMin( val+1 );
      
//       decision.push((val << 1) + dir);
//     }
//     inline void right() { 
//       int val = decision.pop();
//       int dir = (val & 1);
//       val >>= 1;
//       if(dir) _X->setMin( val+1 ); else _X->setMax( val ); 
//     }
//     inline void reverse_left() {
//       int val = _X->min();
//       int dir = (randint(1000) < proba);

//       val += randint(_X->max() - val);
//       dir ^= (value > val);

//       if(dir) _X->setMin( val+1 ); 
//       else _X->setMax( val );
      
//       decision.push((val << 1) + dir);
//     }
//     inline void reverse_right() { 
//       int val = decision.pop();
//       int dir = (val & 1);
//       val >>= 1;
//       if(dir) _X->setMax( val ); else _X->setMin( val+1 ); 
//     }
//     void postCut( const int p ) { }
//     //@}  



//     /**@name Miscellanous*/
//     //@{
//     virtual void printLeft(std::ostream& o) const { 

//       int val = decision.back();
//       int dir = (val & 1);
//       val >>= 1;

//       if(dir) {
// 	std::cout << " <= " << val << " (" << value << ")";
//       } else {
// 	std::cout << " >= " << (val+1) << " (" << value << ")";
//       }
//     }
//     virtual void printRight(std::ostream& o) const {
//       //_X->printshort(o);

//       int val = decision[decision.size];
//       int dir = (val & 1);
//       val >>= 1;

//       if(dir) {
// 	o << " >= " << (val+1) << " (" << value << ")";
//       } else {
// 	o << " <= " << (val) << " (" << value << ")";
//       }
//     }
//     //@}
//   };



  // class ValSelectorWeight : public ValSelector
  // {
  //  public:

  //   /**@name Parameters*/
  //   //@{   
  //   double *weight;
  //   int val;
  //   //@}

  //   /**@name Constructors*/
  //   //@{
  //   ValSelectorWeight( VariableInt *x, double *w ) : ValSelector(x) { weight = w; }
  //   //@}

  //   /**@name Utils*/
  //   //@{ 
  //   inline int getBest() 
  //   { 
  //     DomainIterator *valit = _X->begin();
  //     val = *valit;
  //     double w = weight[val];
  //     while( valit->next() ) {
  //       if( weight[*valit] > w )
  // 	{
  // 	  w = weight[*valit];
  // 	  val = *valit;
  // 	}
  //     }
  //     return val;
  //   }
  //   virtual void getOrder( int* dtv, int* vtd ) {
  //     ValSelector::getOrder( dtv, vtd );

  //     int i, j, aux, n=_X->domsize();
  //     for(i=1; i<n; ++i) {
  //       j = i;
  //       while( j && weight[dtv[j]] > weight[dtv[j-1]] ) {
  // 	aux = dtv[j] ;
  // 	dtv[j] = dtv[j-1] ;
  // 	dtv[j-1] = aux;
  // 	--j;
  //       }
  //     }
  //     for(i=0; i<n; ++i) {
  //       vtd[dtv[i]] = i;
  //     }
  //   }
  //   inline int getBestConst() const 
  //   { 
  //     DomainIterator *valit = _X->begin();
  //     int v_ = *valit;
  //     double w = weight[v_];
  //     while( valit->next() ) {
  //       if( weight[*valit] > w )
  // 	{
  // 	  w = weight[*valit];
  // 	  v_ = *valit;
  // 	}
  //     }
  //     return v_;
  //   }
  //   inline void left() 
  //   { 
  //     _X->setDomain( getBest() );    
  //   }
  //   inline bool right() { return _X->remove( val ); }
  //   //@}

  //   /**@name Miscellanous*/
  //   //@{
  //   virtual void printLeft(std::ostream& o) const { 
  //     _X->printshort(o);
  //     std::cout << " == " << getBestConst();
  //   }
  //   virtual void printRight(std::ostream& o) const {
  //     _X->printshort(o);
  //     std::cout << " =/= " << getBestConst();
  //   }
  //   //@}  

  // };


  class ValSelectorWeight : public ValSelector
  {
  public:

    /**@name Parameters*/
    //@{   
    double *weight;
    int val;
    Vector< int > decision;
    //@}

    /**@name Constructors*/
    //@{
    ValSelectorWeight( VariableInt *x, double *w ) : ValSelector(x) { 
      weight = w; 
      decision.init(0,8);
    }
    //@}

    /**@name Utils*/
    //@{ 
    void make(int& t, int& v) { 
      DomainIterator *valit = _X->begin();
      v = *valit;
      double w = weight[v];
      while( valit->next() ) {
	if( weight[*valit] > w )
	  {
	    w = weight[*valit];
	    v = *valit;
	  }
      }
      t=Decision::ASSIGNMENT;
    }
    inline int getBest() 
    { 
      DomainIterator *valit = _X->begin();
      val = *valit;
      double w = weight[val];
      while( valit->next() ) {
	if( weight[*valit] > w )
	  {
	    w = weight[*valit];
	    val = *valit;
	  }
      }
      return val;
    }
    virtual void getOrder( int* dtv, int* vtd ) {
      ValSelector::getOrder( dtv, vtd );

      int i, j, aux, n=_X->domsize();
      for(i=1; i<n; ++i) {
	j = i;
	while( j && weight[dtv[j]] > weight[dtv[j-1]] ) {
	  aux = dtv[j] ;
	  dtv[j] = dtv[j-1] ;
	  dtv[j-1] = aux;
	  --j;
	}
      }
      for(i=0; i<n; ++i) {
	vtd[dtv[i]] = i;
      }
    }
    inline int getBestConst() const 
    { 
      DomainIterator *valit = _X->begin();
      int v_ = *valit;
      double w = weight[v_];
      while( valit->next() ) {
	if( weight[*valit] > w )
	  {
	    w = weight[*valit];
	    v_ = *valit;
	  }
      }
      return v_;
    }
    inline void left() { 
      _X->setDomain( getBest() ); 
    }
    inline void right() { 
      _X->remove( val ); 
    }
    inline void reverse_left() {
      val = getBest();
      decision.push( val );
      _X->remove( val ); 
    }
    inline void reverse_right() { 
      val = decision.pop();
      _X->setDomain( val );
    }
    void postCut( const int p ) { val = p; }
    //@}

    /**@name Miscellanous*/
    //@{
    virtual void printLeft(std::ostream& o) const { 
      std::cout << " == " << val;
    }
    virtual void printRight(std::ostream& o) const {
      std::cout << " =/= " << val;
    }
    //@}  

  };

  class ValSelectorSac : public ValSelector
  {
  public:

    /**@name Parameters*/
    //@{   
    BitSet *visited;
    int val;
    //@}

    /**@name Constructors*/
    //@{
    ValSelectorSac( VariableInt *x, BitSet *s ) : ValSelector(x) { visited = s; }
    //@}

    /**@name Utils*/
    //@{ 
    void make(int& t, int& v) { 
      DomainIterator *valit = _X->begin();
      v = *valit;
      bool isin = true;
      while( visited->member( *valit ) && ( isin = valit->next() ) );
      if( isin ) v = *valit;
      t=Decision::ASSIGNMENT;
    }
    inline int getBest() 
    {
      DomainIterator *valit = _X->begin();
      val = *valit;
      bool isin = true;
      while( visited->member( *valit ) && ( isin = valit->next() ) );
      if( isin )
	val = *valit;
      return val;
    }
    inline void left() 
    { 
      _X->setDomain( getBest() );
    }
    inline void right() { 
      _X->remove( val ); 
    }
    inline void reverse_left() {
//       val = getBest();
//       //decision.push( val );
//       _X->remove( val ); 
    }
    inline void reverse_right() { 
//       //val = decision.pop();
//       _X->setDomain( val );
    }
    void postCut( const int p ) { val = p; }
    //@}

    /**@name Miscellanous*/
    //@{
    virtual void printLeft(std::ostream& o) const { 
      DomainIterator *valit = _X->begin();
      bool isin = true;
      while( visited->member( *valit ) && ( isin = valit->next() ) );

      std::cout << " == " << (*valit);
      std::cout << " : ";
      visited->print( std::cout );
    }
    virtual void printRight(std::ostream& o) const {
      std::cout << " =/= " << val;
    }
    //@}  
  };

  class ValSelectorSacRange : public ValSelector
  {
  public:

    /**@name Parameters*/
    //@{   
    BitSet *visited;
    int val;
    //@}

    /**@name Constructors*/
    //@{
    ValSelectorSacRange( VariableInt *x, BitSet *s ) : ValSelector(x) { visited = s; }
    //@}

    /**@name Utils*/
    //@{
    void make(int& t, int& v) { 
      if( visited->member(1) ) v=_X->min();
      else v=_X->max();
      t=Decision::ASSIGNMENT;
    } 
    inline int getBest() { return _X->min(); }
    inline void left() 
    { 
      val = ( visited->member(1) );
      if( val ) 
	_X->setMax( _X->min() );
      else
	_X->setMin( _X->max() );
    }
    //   inline bool right() 
    //   {
    //     bool c = true;
    //     if( val )
    //       c = _X->setMin( _X->min()+1 );
    //     else
    //       c = _X->setMax( _X->max()-1 );
    //     return c;
    //   }
    inline void right() 
    {
      if( val )
	_X->setMin( _X->min()+1 );
      else
	_X->setMax( _X->max()-1 );
    }
    inline void reverse_left() {

    }
    inline void reverse_right() { 
 
    }
    void postCut( const int p ) { val = p; }
    //@}

    /**@name Miscellanous*/
    //@{
    virtual void printLeft(std::ostream& o) const { 
      _X->printshort(o);
      int v = ( visited->member(1) );
      if( v )
	std::cout << " >= " << _X->max();
      else 
	std::cout << " <= " << _X->min();
      std::cout << " : ";
      visited->print( std::cout );
    }
    virtual void printRight(std::ostream& o) const {
      _X->printshort(o);
      int v = ( visited->member(1) );
      if( v )
	std::cout << " < " << _X->max();
      else 
	std::cout << " > " << _X->min();
    }
    //@}  
  };


  class ValSelectorMOSP : public ValSelector
  {
  public:

    /**@name Parameters*/
    //@{
    int val;
    PredicateDisjunctive *disjunct;
    //@}

    /**@name Constructors*/
    //@{
    ValSelectorMOSP( VariableInt *x, PredicateDisjunctive *d ) : ValSelector(x), disjunct(d) {}
    //@}

    /**@name Utils*/
    //@{
    void make(int& t, int& v) { 
      v = ( disjunct->XltY() < disjunct->YltX() );
      t=Decision::ASSIGNMENT;
    } 
    inline int getBest() { 
      val = ( disjunct->XltY() < disjunct->YltX() );
      return val;
    }
    inline void left() { 
      val = ( disjunct->XltY() < disjunct->YltX() );
      _X->setDomain( val );
    }
    inline void right() { 
      _X->remove( val ); 
    }
    inline void reverse_left() {
      val = ( disjunct->XltY() < disjunct->YltX() );
      _X->remove( val ); 
    }
    inline void reverse_right() { 
      _X->setDomain( val );
    }
    void postCut( const int p ) {
      val = p;
    }
    //@}  

    /**@name Miscellanous*/
    //@{
    virtual void printLeft(std::ostream& o) const { 
      std::cout << " == " << val;
    }
    virtual void printRight(std::ostream& o) const {
      std::cout << " =/= " << val;
    }
    //@}
  };


  class ValSelectorLOSP : public ValSelector
  {
  public:

    /**@name Parameters*/
    //@{
    int val;
    PredicateDisjunctive *disjunct;
    //@}

    /**@name Constructors*/
    //@{
    ValSelectorLOSP( VariableInt *x, PredicateDisjunctive *d ) : ValSelector(x), disjunct(d) {}
    //@}

    /**@name Utils*/
    //@{
 
    void make(int& t, int& v) { 
      v = ( disjunct->XltY() > disjunct->YltX() );
      t=Decision::ASSIGNMENT;
    } 
    inline int getBest() { 
      val = ( disjunct->XltY() > disjunct->YltX() );
      return val;
    }
    inline void left() { 
//       disjunct->print(std::cout);
//       std::cout << " ";
//       disjunct->scope[0]->print(std::cout);
//       std::cout << " ";
//       disjunct->scope[1]->print(std::cout);
//       std::cout << " ";


      val = ( disjunct->XltY() > disjunct->YltX() );
      _X->setDomain( val );

//       disjunct->print(std::cout);
//       std::cout << std::endl;
    }
    inline void right() { 
      _X->remove( val ); 
    }
    inline void reverse_left() {
      val = ( disjunct->XltY() > disjunct->YltX() );
      _X->remove( val ); 
    }
    inline void reverse_right() { 
      _X->setDomain( val );
    }
    void postCut( const int p ) {
      val = p;
    }
    //@}  

    /**@name Miscellanous*/
    //@{
    virtual void printLeft(std::ostream& o) const { 
      std::cout << " == " << val;
    }
    virtual void printRight(std::ostream& o) const {
      std::cout << " =/= " << val;
    }
    //@}
  };



  class ValSelectorPFSP : public ValSelector
  {
  public:

    /**@name Parameters*/
    //@{
    int val;
    PredicateDisjunctive **disjunct;
    int ndsj;
    //@}

    /**@name Constructors*/
    //@{
    ValSelectorPFSP( VariableInt *x, PredicateDisjunctive **d, const int n ) : ValSelector(x), disjunct(d), ndsj(n) {}
    //@}

    /**@name Utils*/
    //@{

    void make(int& t, int& v) { 
      int i, XltY = 0, YltX = 0;
      for(i=0; i<ndsj; ++i) {
	XltY += disjunct[i]->XltY();
	YltX += disjunct[i]->YltX();
      }
      v = ( XltY > YltX );
      t=Decision::ASSIGNMENT;
    }  
    inline int getBest() { 
      int i, XltY = 0, YltX = 0;
      for(i=0; i<ndsj; ++i) {
	XltY += disjunct[i]->XltY();
	YltX += disjunct[i]->YltX();
      }

      val = ( XltY > YltX );
      return val;
    }
    inline void left() { 
      int i, XltY = 0, YltX = 0;
      for(i=0; i<ndsj; ++i) {
	XltY += disjunct[i]->XltY();
	YltX += disjunct[i]->YltX();
      }

      val = ( XltY > YltX );
      _X->setDomain( val );
    }
    inline void right() { 
      _X->remove( val ); 
    }
    inline void reverse_left() {
      int i, XltY = 0, YltX = 0;
      for(i=0; i<ndsj; ++i) {
	XltY += disjunct[i]->XltY();
	YltX += disjunct[i]->YltX();
      }

      val = ( XltY > YltX );
      _X->remove( val ); 
    }
    inline void reverse_right() { 
      _X->setDomain( val );
    }
    void postCut( const int p ) {
      val = p;
    }
    //@}  

    /**@name Miscellanous*/
    //@{
    virtual void printLeft(std::ostream& o) const { 
      _X->printshort(o);
      std::cout << " == " << val;
    }
    virtual void printRight(std::ostream& o) const {
      _X->printshort(o);
      std::cout << " =/= " << val;
    }
    //@}
  };



  class ValSelectorSOSP : public ValSelector
  {
  public:

    /**@name Parameters*/
    //@{
    int mode;
    int val;
    PredicateDisjunctive *disjunct;
    //@}

    /**@name Constructors*/
    //@{
    ValSelectorSOSP( VariableInt *x, PredicateDisjunctive *d, const int m ) : ValSelector(x), disjunct(d) { mode = m; }
    //@}

    /**@name Utils*/
    //@{

    void make(int& t, int& v) { 
      switch(mode) {
      case 1: {
	v = ( disjunct->XltY() > disjunct->YltX() );
	break;
      }
      case -1: {
	v = ( disjunct->XltY() < disjunct->YltX() );
	break;
      }
      case 0: {
	v = 0;
	break;
      }
      default: {
	v = randint(2);
      }
      }
      t=Decision::ASSIGNMENT;
    } 
    inline int getBest() { 
      switch(mode) {
      case 1: {
	//std::cout << "promise" << std::endl;
	val = ( disjunct->XltY() > disjunct->YltX() );
	break;
      }
      case -1: {
	//std::cout << "anti" << std::endl;
	val = ( disjunct->XltY() < disjunct->YltX() );
	break;
      }
      case 0: {
	//std::cout << "lex" << std::endl;
	val = 0;
	break;
      }
      default: {
	//std::cout << "rand" << std::endl;
	val = randint(2);
      }
      }
      return val;
    }
    inline void left() { 
//       switch(mode) {
//       case 1: {
// 	//std::cout << "promise" << std::endl;
// 	val = ( disjunct->XltY() > disjunct->YltX() );
// 	break;
//       }
//       case -1: {
// 	//std::cout << "anti" << std::endl;
// 	val = ( disjunct->XltY() < disjunct->YltX() );
// 	break;
//       }
//       case 0: {
// 	//std::cout << "lex" << std::endl;
// 	val = 0;
// 	break;
//       }
//       default: {
// 	//std::cout << "rand" << std::endl; 
// 	val = randint(2);
//       }
//       }
      val = getBest();
      _X->setDomain( val );
    }
    inline void right() { 
      _X->remove( val ); 
    }
    inline void reverse_left() {
      val = getBest();
      _X->remove( val ); 
    }
    inline void reverse_right() { 
      _X->setDomain( val );
    }
    void postCut( const int p ) {
      val = p;
    }
    //@}  

    /**@name Miscellanous*/
    //@{
    virtual void printLeft(std::ostream& o) const { 
      _X->printshort(o);
      std::cout << " == " << val;
    }
    virtual void printRight(std::ostream& o) const {
      _X->printshort(o);
      std::cout << " =/= " << val;
    }
    //@}
  };



  /**********************************************
   * DVOs Dynamic Variable Ordering
   **********************************************/
  /*! \class DVO
    \brief  Representation of Dynamic Variable Ordering.

    "Future variables", that is, those that are not
    yet assigned are stored between the pointers 'first'
    and 'last' The function select() return the  
    next variable to branch on.
  */
  class DVO {
  public:

    /**@name Parameters*/
    //@{ 
    /// A reference to the variables in lex static order
    VariableInt**& variables;
    /// A reference to the beginning of the array of future variables
    VariableInt**& first;
    /// A reference to the end of the array of future variables
    VariableInt**&  last;
    /// The list of decisions
    VariableInt**& decision;
    /// The number of decisions;
    int& level;
    //int verbosity;
    int limit;
    /// Structure used to learn weights
    //Weighter *learner;
    //@}
    
    PredicateDisjunctive **get_disjuncts();
    virtual double get_value(VariableInt *X) { return 1.0; }

    /**@name Constructors*/
    //@{
    // virtual void print_weights() {}
    
    DVO(Solver*, const int l=NOVAL);
    virtual ~DVO() {}
    //@}  

    /**@name Utils*/
    //@{ 
    /*!
      The variable selection method.
      Returns the next variable to pick,
      a value may also be selected by assigning the parameter, if no
      value is chosen then the default value ordering will be used
      (depending on the Variable's implementation).
    */
    virtual VariableInt* select() = 0;
    //virtual void print( std::ostream& o ) const { o << "DVO" << std::endl; }  
    //@} 
  };



  /**********************************************
   * Generic heuristic
   **********************************************/

  template <class T>
  class GenericDVO : public DVO
  {
  public: 
    
    /**@name Parameters*/
    //@{ 
    T best;
    T current;
    //@}

    virtual double get_value(VariableInt *X) {
      current = X;
      return current.value();
    }

    /**@name Constructors*/
    //@{
    GenericDVO(Solver* s, const int l=NOVAL) : DVO(s,l) { }
    virtual ~GenericDVO()
    {
      //delete [] _garbage_disjuncts;
      free_disjuncts();
    }
    //@}

    /**@name Utils*/
    //@{
    
//     virtual void print_weights() {
//       VariableInt **var_iterator;
//       for( var_iterator = first; var_iterator != last; ++var_iterator ) 
// 	{
// 	  current = (*var_iterator);
// 	  std::cout << "branch: " << (current.value()) << std::endl;
// 	}
//     }
 
    inline VariableInt* select()
    {    
      VariableInt *var=*first, **var_iterator;
      VariableInt **last_var = (last-first > limit ? first+limit : last);
      best = var; 

      //std::cout << "number of unassigned vars: " << (last-first) << std::endl;

      //for( var_iterator = first+1; var_iterator != last; ++var_iterator ) 
      for( var_iterator = first+1; var_iterator != last_var; ++var_iterator ) 
	{  
	
	  // 	  	VariableInt *x = (*var_iterator);
	  // 	  	int var_weight = x->weight;
	  // 	  	int con_weight = 0;
	  // 	  	MistralNode<Constraint*> *nd = x->constraintsOnValue();
	  // 	  	while( nextNode(nd) ) {
	  // 	  	  //nd->elt->print( std::cout );
	  // 	  	  //std::cout << std::endl;
	  // 	  	  con_weight += nd->elt->weight;
	  // 	  	}

	  // 		if(var_weight != 1)
	  // 		  std::cout << var_weight << " " << con_weight << std::endl;

	  // 	  	assert( var_weight == con_weight );

	  current = (*var_iterator);
	  if( current < best ) 
	    {
	      best = current;
	      var = (*var_iterator);


// 	      //std::cout << "jjjjjjj " << verbosity << std::endl;
// 	      if(verbosity > 3) {
// 		var->print(std::cout);
// 		std::cout << " " << best.value() << std::endl;
// 	      } 
	    }
	}

      return var;
    }
    //@}
  };


  /**********************************************
   * Generic Scheduling heuristic
   **********************************************/

  template <class T>
  class GenericSchedulingDVO : public GenericDVO< T >
  {
  public: 

    /**@name Parameters*/
    //@{ 
    PredicateDisjunctive **the_disjuncts;
    //@}

    /**@name Constructors*/
    //@{
    GenericSchedulingDVO(Solver* s, const int l=NOVAL) : GenericDVO< T >(s,l) {}
    virtual ~GenericSchedulingDVO()
    {
      delete [] the_disjuncts;
    }
    //@}
  };


  /**********************************************
   * GenericRandom heuristic
   **********************************************/

  template <class T>
  class GenericRandomDVO : public DVO
  {
  public: 

    /**@name Parameters*/
    //@{ 
    T *bests;
    T current;
    VariableInt **bestvars;
    int size;
    //@}

    virtual double get_value(VariableInt *X) {
      current = X;
      return current.value();
    }

    /**@name Constructors*/
    //@{
    GenericRandomDVO(Solver* s, const int sz, const int l=NOVAL)  : DVO(s,l) 
    {
      size = sz;
      bests = new T[size+1];
      bestvars = new VariableInt*[size+1];
      for(int i=0; i<=size; ++i) bestvars[i] = NULL;
    }

    virtual ~GenericRandomDVO()
    {
      delete [] bests;
      delete [] bestvars;
      free_disjuncts();
      //delete [] _garbage_disjuncts;
    }
    /*
      {
      delete [] bests;
      delete [] bestvars;
      delete [] disjuncts;
      }
    */
    //@}

    /**@name Utils*/
    //@{ 
    inline VariableInt* select()
    {
      int realsize=1, i;

      //std::cout << " ===> first " << (*first) << std::endl;
      bests[0] = bestvars[0] = *first;
      VariableInt **var_iterator;
      VariableInt **last_var = (last-first > limit ? first+limit : last);


//       for( var_iterator = first; var_iterator != last; ++var_iterator ) {
// 	if((*var_iterator)->weight > 1000) {
// 	  (*var_iterator)->print(std::cout);
// 	  std::cout << " > 1000" <<std::endl;
// 	} 
//       }

      //bests[0].print(std::cout);
      //std::cout << " ";
      //(*first)->print(std::cout);
      //std::cout << std::endl;
      //Solver *solver = (*first)->solver;
      //std::cout << solver->getNodes() << std::endl;

      //
      //std::cout << " ===> last " << (last) << std::endl;
      //for( var_iterator = first+1; var_iterator != last; ++var_iterator ) 
      for( var_iterator = first+1; var_iterator != last_var; ++var_iterator ) 
	{  
	  //std::cout << " ===> iter " << (var_iterator) << std::endl;

//   	  std::cout << "|";
//  	  std::cout.flush();
// //  	  (*var_iterator)->print(std::cout);
// // 	  std::cout.flush();
	  current = (*var_iterator);
// // 	  std::cout.flush();
//  	  std::cout << " ";
//  	  std::cout.flush();
//  	  current.print(std::cout);
//  	  std::cout.flush();
//  	  std::cout << "| ";
//  	  std::cout.flush();

// 	  if( (*var_iterator)->getType() != VariableInt::BIT ) {
// 	    (*var_iterator)->print(std::cout);
// 	    std::cout << "is weird" << std::endl;
// 	  }


	  // 	VariableInt *x = (*var_iterator);
	  // 	int var_weight = x->weight;
	  // 	int con_weight = 0;
	  // 	MistralNode<Constraint*> *nd = x->constraintsOnDomain();
	  // 	while( nextNode(nd) ) {
	  // 	  con_weight += nd->elt->weight;
	  // 	}
	  // 	assert( var_weight == con_weight );

	  //std::cout << " ===> rs " << (realsize) << std::endl;

	  i = realsize;
	  //std::cout << " ===> i " << (i-1) << std::endl;
	  while( i && current < bests[i-1] ) {
	    bests[i] = bests[i-1];
	    bestvars[i] = bestvars[i-1];
	    --i;
	    //std::cout << " ===> i " << (i-1) << std::endl;
	  }

	  //std::cout << " ===> iter " << (*var_iterator) << std::endl;
	  bests[i] = current;
	  bestvars[i] = (*var_iterator);

// 	  if(i<realsize || realsize<size) {
// 	    std::cout << "*";
// 	  }
// 	  std::cout << std::endl;
	

	  //std::cout << " ===> siz " << size << std::endl;
	  if(realsize<size) ++realsize;
	}

      //std::cout << std::endl;    

      //std::cout << std::endl;
      //return bestvars[rand() % realsize];
      return bestvars[randint(realsize)];
    }
    //@}
  };


  /**********************************************
   * Generic Scheduling Random heuristic
   **********************************************/

  template <class T>
  class GenericSchedulingRandomDVO : public GenericRandomDVO< T >
  {
  public: 

    /**@name Parameters*/
    //@{ 
    PredicateDisjunctive **the_disjuncts;
    //@}

    /**@name Constructors*/
    //@{
    GenericSchedulingRandomDVO(Solver* s, const int sz) : GenericRandomDVO< T >(s, sz) {}
    virtual ~GenericSchedulingRandomDVO()
    {
      delete [] the_disjuncts;
    }
    //@}
  };

  /*! \class VariableOrdering 
    \brief  Wrapper for DVO's
  */
  class VariableOrdering {
  public:

    virtual ~VariableOrdering() {}

    /**@name Utils*/
    //@{
    virtual DVO* extract( Solver* ) = 0 ;
    //@}
  };


  /**********************************************
   * No Ordering
   **********************************************/

  /*! \class DVONoOrder
    \brief  No variable ordering

    Variables are explored in no order 
  */
  class DVONoOrder : public DVO {
  public:
    /**@name Constructors*/
    //@{
    DVONoOrder(Solver* s) : DVO(s) {}
    //@} 

    /**@name Utils*/
    //@{    
    /**!
       The variable selection method.
       Return the first unassigned variable
    */
    VariableInt* select();
    //@} 
  };

  /*! \class No
    \brief  Wrapper for DVONo
  */
  class NoOrder : public VariableOrdering {
  public:

    virtual ~NoOrder() {}

    /**@name Utils*/
    //@{  
    DVO* extract( Solver* s )
    {
      return new DVONoOrder(s);
    }
    //@}
  };


  /**********************************************
   * Lexicographic Ordering
   **********************************************/

  /*! \class DVOLexicographic
    \brief  Lexico variable ordering

    Variables are explored in lexicographic order 
  */
  class DVOLexicographic : public DVO {
  public:

    //ReversibleInt last;
    VariableInt     **sequence;
    ReversibleNum<int> lastIdx;

    /**@name Constructors*/
    //@{
    DVOLexicographic(Solver* s) ;
    virtual ~DVOLexicographic() ;
    //@} 

    /**@name Utils*/
    //@{    
    /**!
       The variable selection method.
       Return the first unassigned variable
       in lexico order
    */
    VariableInt* select();
    //@} 
  };


  /**********************************************
   * AntiLex Ordering
   **********************************************/

  /*! \class DVOAntiLex
    \brief  Lexico variable ordering

    Variables are explored in lexicographic order 
  */
  class DVOAntiLex : public DVO {
  public:

    //ReversibleInt last;
    VariableInt     **sequence;
    ReversibleNum<int> lastIdx;

    /**@name Constructors*/
    //@{
    DVOAntiLex(Solver* s) ;
    virtual ~DVOAntiLex() ;
    //@} 

    /**@name Utils*/
    //@{    
    /**!
       The variable selection method.
       Return the first unassigned variable
       in lexico order
    */
    VariableInt* select();
    //@} 
  };

  /*! \class Lexicographic
    \brief  Wrapper for DVOLexicographic
  */
  class Lexicographic : public VariableOrdering {
  public:

    /**@name Utils*/
    //@{  
    DVO* extract( Solver* s )
    {
      return new DVOLexicographic(s);
    }
    //@}
  };

  /*! \class Lexicographic
    \brief  Wrapper for DVOLexicographic
  */
  class AntiLex : public VariableOrdering {
  public:

    /**@name Utils*/
    //@{  
    DVO* extract( Solver* s )
    {
      return new DVOAntiLex(s);
    }
    //@}
  };


  /**********************************************
   * Random Ordering
   **********************************************/

  /*! \class DVORandom
    \brief  Random variable ordering

    Variables are explored in lexicographic order
  */
  class DVORandom : public DVO {
  public:
    /**@name Constructors*/
    //@{
    DVORandom(Solver* s);//, unsigned int seed) ;
    //@}

    /**@name Utils*/
    //@{
    /**!
       The variable selection method.
       Return the first unassigned variable
    */
    virtual VariableInt* select();
    //@}
  };

  /*! \class Random
    \brief  Wrapper for DVORandom
  */
  class Random : public VariableOrdering {
  public:

    /**@name Parameters*/
    //@{ 
    //unsigned int seed;
    //@}

    /**@name Utils*/
    //@{
    DVO* extract( Solver* s )
    {
      return new DVORandom(s);//, seed);
    }
    //@}
  };


  /**********************************************
   * Min Domain
   **********************************************/

  /*! \class MinDomain
    \brief  Wrapper for DVOMinDomain
  */
  class MinDomain : public VariableOrdering {
  public:

    int size;
    MinDomain( const int sz=0 ) { size=sz; }
  
    /**@name Utils*/
    //@{
    DVO* extract( Solver* s )
    {
      if( size > 1 )
	return new GenericRandomDVO<VarSelectorDomain>(s, size);
      else 
	return new GenericDVO<VarSelectorDomain>(s);
    }
    //@}
  };


  /**********************************************
   * Min Domain, Min Value
   **********************************************/

  /*! \class MinDomMin
    \brief  Wrapper for DVOMinDomMin
  */
  class MinDomMin : public VariableOrdering {
  public:
  
    int size;
    MinDomMin( const int sz=0 ) { size=sz; }

    /**@name Utils*/
    //@{
    DVO* extract( Solver* s )
    {
      if( size > 1 )
	return new GenericRandomDVO<VarSelectorDomainMinVal>(s, size);
      else
	return new GenericDVO<VarSelectorDomainMinVal>(s);
    }
    //@}
  };


  /**********************************************
   * Max Degree
   **********************************************/

  /*! \class MaxDegree
    \brief  Wrapper for DVOMaxDegree
  */
  class MaxDegree : public VariableOrdering {
  public:
  
    int size;
    MaxDegree( const int sz=0 ) { size=sz; }

    /**@name Utils*/
    //@{
    DVO* extract( Solver* s )
    {
      if( size > 1 )
	return new GenericRandomDVO<VarSelectorDegree>(s, size);
      else 
	return new GenericDVO<VarSelectorDegree>(s);
    }
    //@}
  };


  /**********************************************
   * Min Domain , Max Degree
   **********************************************/

  /*! \class MinDomMaxDeg
    \brief  Wrapper for DVOMinDomMaxDeg
  */
  class MinDomMaxDeg : public VariableOrdering {
  public:
  
    int size;
    MinDomMaxDeg( const int sz=0 ) { size=sz; }
  
    /**@name Utils*/
    //@{
    DVO* extract( Solver* s )
    {
      if( size > 1 )
	return new GenericRandomDVO<VarSelectorDomainDegree>(s, size);
      else
	return new GenericDVO<VarSelectorDomainDegree>(s);
    }
    //@}
  };


  /**********************************************
   * Min domain/degree
   **********************************************/

  /*! \class DomOverDeg
    \brief  Wrapper for DVODomOverDeg
  */
  class DomOverDeg : public VariableOrdering {
  public:
  
    int size;
    DomOverDeg( const int sz=0 ) { size=sz; }

    /**@name Utils*/
    //@{
    DVO* extract( Solver* s )
    {
      if( size > 1 )
	return new GenericRandomDVO<VarSelectorDomainOverDegree>(s, size);
      else
	return new GenericDVO<VarSelectorDomainOverDegree>(s);
    }
    //@}
  };


  /**********************************************
   * Min Domain/degree on neighbors
   **********************************************/

  /*! \class Neighbor
    \brief  Wrapper for DVONeighbor
  */
  class Neighbor : public VariableOrdering {
  public:
  
    int size;
    Neighbor( const int sz=0 ) { size=sz; }

    /**@name Utils*/
    //@{
    DVO* extract( Solver* s )
    {
      if( size > 1 )
	return new GenericRandomDVO<VarSelectorNeighborDD>(s, size);
      else
	return new GenericDVO<VarSelectorNeighborDD>(s);
    }
    //@}
  };


  /**********************************************
   * Min domain/weighted degree (level)
   **********************************************/

  /*! \class DomOverWLDeg
    \brief  Wrapper for DVODomOverWLDeg
  */
  class DomOverWLDeg : public VariableOrdering {
  public:

    int size;
    DomOverWLDeg( const int sz=0 ) { size=sz; }

    /**@name Utils*/
    //@{
    DVO* extract( Solver* s );
    //@}
  };


  /**********************************************
   * Min domain/weighted degree
   **********************************************/

  /*! \class DomOverWDeg
    \brief  Wrapper for DVODomOverWDeg
  */
  class DomOverWDeg : public VariableOrdering {
  public:

    int size;
    DomOverWDeg( const int sz=0 ) { size=sz; }

    /**@name Utils*/
    //@{
    DVO* extract( Solver* s );
    //@}
  };


  /**********************************************
   * Max impact
   **********************************************/

  /*! \class Impact
    \brief  Wrapper for DVOImpact
  */
  class Impact : public VariableOrdering {
  public:
  
    int size;
    Impact( const int sz=0 ) { size=sz; }

    /**@name Utils*/
    //@{
    DVO* extract( Solver* s );
    //@}
  };


  /**********************************************
   * Max impact * degree
   **********************************************/

  /*! \class ImpactOverDeg
    \brief  Wrapper for DVOImpactOverDeg
  */
  class ImpactOverDeg : public VariableOrdering {
  public:

    int size;
    ImpactOverDeg( const int sz=0 ) { size=sz; }

    /**@name Utils*/
    //@{
    DVO* extract( Solver* s );
    //@}
  };


  /**********************************************
   * Max impact * weighted degree
   **********************************************/

  /*! \class ImpactOverWDeg
    \brief  Wrapper for DVOImpactOverWDeg
  */
  class ImpactOverWDeg : public VariableOrdering {
  public:

    int size;
    ImpactOverWDeg( const int sz=0 ) { size=sz; }

    /**@name Utils*/
    //@{
    DVO* extract( Solver* s );
    //@}
  };


  /**********************************************
   * Max impact * weighted degree (level)
   **********************************************/

  /*! \class ImpactOverWLDeg
    \brief  Wrapper for DVOImpactOverWLDeg
  */
  class ImpactOverWLDeg : public VariableOrdering {
  public:

    int size;
    ImpactOverWLDeg( const int sz=0 ) { size=sz; }

    /**@name Utils*/
    //@{  
    DVO* extract( Solver* s );
    //@}
  };


  /**********************************************
   * Probing heuristic
   **********************************************/

  /*! \class Probedvo
    \brief  Wrapper for Probing heuristic
  */
  class Probedvo : public VariableOrdering {
  public:

    /**@name Parameters*/
    //@{ 
    //unsigned int seed;  
    int ltype;
    int htype;
    //@}

    /**@name Constructors*/
    //@{
    Probedvo( /* unsigned int graine=11041979, */ int lt=Weighter::WDG, int ht=0 )
    {
      //seed = graine;
      ltype = lt;
      htype = ht;
    }

    Probedvo( /* unsigned int graine=11041979, */ bool ulw, bool uip, int ht )
    {
      //seed = graine;
      htype = ht;

      ltype = Weighter::WDG;
      if( ulw ) 
	ltype |= Weighter::WLD;
      if( uip )
	ltype |= Weighter::IPT;
    }
    //@}

    /**@name Utils*/
    //@{
    virtual DVO* extract( Solver* s );
    //@}
  };


  /**********************************************
   * Singleton Arc Consistency
   **********************************************/

  /*! \class DVOSingletonAC
    \brief SAC variable ordering

    Choose first the variables and values not
    yet SAC.
  */
  class DVOSingletonAC : public DVO {
  public:

    /**@name Parameters*/
    //@{ 
    BitSet *SACdomain;
    bool *isRange;
    //  unsigned int *sacsize;
    //@}

    /**@name Constructors*/
    //@{
    DVOSingletonAC(Solver*, WeighterSAC*) ;//: DVO(s) {}
    //@}

    /**@name Utils*/
    //@{
    VariableInt* select();
    //@}
  };


  /**********************************************
   * OSP Heuristic
   **********************************************/

  /*! \class OSP
    \brief  Wrapper for DVOOSP
  */
  class OSP : public VariableOrdering {
  public:
    static const int DOMAIN_P_TWEIGHT = 4;
    static const int DOMAIN_O_NOT = 3;
    static const int DOM_O_BOOLWEIGHT = 0;
    static const int DOM_O_TASKWEIGHT = 1;
    static const int DOM_O_BOOLTASKWEIGHT = 2;
    static const int DOMAIN_O_NOTTYPE = 8;
    static const int DOM_O_BOOLWEIGHTTYPE = 5;
    static const int DOM_O_TASKWEIGHTTYPE = 6;
    static const int DOM_O_BOOLTASKWEIGHTTYPE = 7;
    //static const int DOM_O_WEAKWEIGHT = 3;

    int size;
    int promise;
    int strategy;
    OSP( const int sz=0, const int pr=0, const int st=DOM_O_BOOLWEIGHT ) { size=sz; promise=pr; strategy=st;}
    virtual ~OSP();
  
    /**@name Utils*/
    //@{
    PredicateDisjunctive** get_disjuncts(Solver *s);
    DVO* extract( Solver* s );
    //@}
  };


  /**********************************************
   * PFSP Heuristic
   **********************************************/

  /*! \class PFSP
    \brief  Wrapper for DVOPFSP
  */
  class PFSP : public VariableOrdering {
  public:

    int size;
    int promise;

    PFSP( const int sz=0, const int pr=0 ) { size=sz; promise=pr; }
  
    /**@name Utils*/
    //@{
    void get_disjuncts(Solver* s,
		       PredicateDisjunctive***& disjunct,
		       int*& sdegree);
    DVO* extract( Solver* s );
    //@}
  };


  /**********************************************
   * OSPSAT Heuristic
   **********************************************/

  /*! \class OSPSAT
    \brief  Wrapper for DVOOSPSAT
  */
  class OSPSAT : public VariableOrdering {
  public:

    int size;
    OSPSAT( const int sz=0 ) { size=sz; }
  
    /**@name Utils*/
    //@{
    DVO* extract( Solver* s );
    //@}
  };


  /**********************************************
   * FPP Heuristic
   **********************************************/

  /*! \class FPP
    \brief  Wrapper for DVOFPP
  */
  class FPP : public VariableOrdering {
  public:

    int size;
    int promise;
    FPP( const int sz=0, const int pr=0 ) { size=sz; promise=pr; }
  
    /**@name Utils*/
    //@{
    DVO* extract( Solver* s );
    //@}
  };


  /**********************************************
   * Objective Functions 
   **********************************************/

  /*! \class ObjectiveFunction 
    \brief  Abstract objective function class

    Always minimise a criterion. The upper bound is 
    updated when a solution is found. 
  */
  class ObjectiveFunction {
  public:
    /**@name Parameters*/
    //@{ 
    /// The current upper bound on the objective to minimise
    int upper_bound;
    //@}

    virtual ~ObjectiveFunction() {}

    /**@name Utils*/
    //@{   
    /// Return an optimistic evaluation of the best objective value currently achievable
    virtual int score() = 0;
    /// Return the current value of the objective
    virtual int solution_score() = 0;
    /// Called by the solver whenever a solution is found 
    virtual int update();
    //@}

    virtual void print( std::ostream& o ) const {}  
  };


  /*! \class MaximiseVar
    \brief  Maximise a variable.
  */
  class MaximiseVar : public ObjectiveFunction {

  public:
    /**@name Parameters*/
    //@{ 
    /// The unary constraint used to bound search
    UnaryConstraintMore *umore;
    /// The best achievable value, used to terminate early
    int maxX;
    //@}  

    /**@name Constructors*/
    //@{ 
    MaximiseVar(Solver*, VariableInt*);
    virtual ~MaximiseVar();
    //@}  

    /**@name Utils*/
    //@{ 
    /*!
      returns an optimistic estimate of the bound
    */
    virtual int score();
    virtual int solution_score();
    /*!
      The score of the solution is computed, then
      the unary constraint is updated consequently.
    */
    virtual int update();
    //@}  

    virtual void print( std::ostream& o ) const { 
      o << "maximise " ;
      umore->X->print(o);
    }  
  };


  /*! \class MinimiseVar
    \brief  Minimise a variable.
  */
  class MinimiseVar : public ObjectiveFunction {

  public:
    /**@name Parameters*/
    //@{ 
    /// The unary constraint used to bound search
    UnaryConstraintLess *uless;
    /// The best achievable value, used to terminate early
    int minX;
    //@}  

    /**@name Constructors*/
    //@{ 
    MinimiseVar(Solver*, VariableInt*);
    virtual ~MinimiseVar();
    //@}  

    /**@name Utils*/
    //@{ 
    /*!
      returns an optimistic estimate of the bound
    */
    virtual int score();
    virtual int solution_score();
    /*!
      The score of the solution is computed, then
      the unary constraint is updated consequently.
    */
    virtual int update();
    //@}  

    virtual void print( std::ostream& o ) const { 
      o << "minimise " ;
      uless->X->print(o);
    }  
  };

};


inline void Mistral::SimpleUnaryConstraint::make() { 
  int v;
  int t;
  var->branch->make(t,v); 
  _data_ = ((v<<2)+t);
}



#endif // _DVO_H

