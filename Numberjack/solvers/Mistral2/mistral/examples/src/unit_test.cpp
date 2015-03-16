
#include <vector>
#include <iomanip>
#include <fstream>

#include <math.h>

#include <mistral_search.hpp>
#include <mistral_variable.hpp>
#include <mistral_constraint.hpp>


using namespace std;
using namespace Mistral;


#define LOW     1
#define MEDIUM  2
#define HIGH    3
#define EXTREME 4

int cur_iteration = 0;
//

class UnitTest {

protected:

  void checkDomainIntegrity(Variable X);

public:

  int Verbosity;
  int Quality;
  int Quantity;
  
  UnitTest();
  UnitTest(const int vb, const int ql, const int qt);
  virtual ~UnitTest();

  virtual void run() = 0;

};


class RandomDomainRandomRemove : public UnitTest {

public:
  
  RandomDomainRandomRemove(const int ql=MEDIUM, const int qt=MEDIUM);
  ~RandomDomainRandomRemove();

  virtual void run();
};


class RandomIntervalTest : public UnitTest {

public:
  
  RandomIntervalTest(const int ql=MEDIUM, const int qt=MEDIUM);
  ~RandomIntervalTest();

  virtual void run();
};


class RandomDomainRandomSetDomainAndRestore : public UnitTest {

public:
  
  RandomDomainRandomSetDomainAndRestore(const int ql=MEDIUM, 
					const int qt=MEDIUM);
  ~RandomDomainRandomSetDomainAndRestore();

  virtual void run();

};


class RandomDomainRandomSetMinAndRestore : public UnitTest {

public:
  
  RandomDomainRandomSetMinAndRestore(const int ql=MEDIUM, 
				     const int qt=MEDIUM);
  ~RandomDomainRandomSetMinAndRestore();

  virtual void run();
};


class RandomDomainRandomSetMaxAndRestore : public UnitTest {

public:
  
  RandomDomainRandomSetMaxAndRestore(const int ql=MEDIUM, 
				     const int qt=MEDIUM);
  ~RandomDomainRandomSetMaxAndRestore();

  virtual void run();
};



class RandomDomainRandomSetDomainBitsetAndRestore : public UnitTest {

public:
  
  RandomDomainRandomSetDomainBitsetAndRestore(const int ql=MEDIUM, 
					      const int qt=MEDIUM);
  ~RandomDomainRandomSetDomainBitsetAndRestore();

  virtual void run();
};



class RandomDomainRandomRemoveRangeAndRestore : public UnitTest {

public:
  
  RandomDomainRandomRemoveRangeAndRestore(const int ql=MEDIUM, 
					  const int qt=HIGH);
  ~RandomDomainRandomRemoveRangeAndRestore();

  virtual void run();
};

// template< int NUM_HEAD >
// class RandomCListRandomRemoveAndRestore : public UnitTest {

// public:

//   int NUM_ELT;
  
//   RandomCListRandomRemoveAndRestore(const int ql=HIGH, 
// 				   const int qt=MEDIUM); 
//   ~RandomCListRandomRemoveAndRestore();

//   class NaiveImplementation : public Reversible {

//   public:
//     Vector<int> elements;
//     Vector<int> whichlist;
//     Vector<int> whichleveladded;
//     Vector<int> whichlevelremoved;

//     NaiveImplementation(Solver *s) ;

//     void create(int elt, int wl=0);
//     void reversible_remove(int elt, int wl);
//     void reversible_add(int elt, int wl=0);
//     void save();
//     void restore();

//     void print(int num_lists);
//   };


//   void checkEquality(ReversibleMultiList<int,NUM_HEAD>& alist, 
// 		     NaiveImplementation& blist);
//   virtual void run();
// };



template< class T >
class RandomRevNumAffectations : public UnitTest {

public:
  
  RandomRevNumAffectations(const int ql=LOW, const int qt=LOW); 
  ~RandomRevNumAffectations();

  virtual void run();  
};


class CostasAllDiffAllSolutions : public UnitTest {

public:
  
  static const int num_sol[13];

  int size;
  int consistency;
  int domain;
  
  CostasAllDiffAllSolutions(const int sz, const int ct, const int dt=DYN_VAR);
  ~CostasAllDiffAllSolutions();

  virtual void run();
};

const int CostasAllDiffAllSolutions::num_sol[] = {0,0,0,4,12,40,116,200,444,760,2160,0,0};


class CostasNotEqualAllSolutions : public UnitTest {

public:

  static const int num_sol[13];
    
  int size;
  
  CostasNotEqualAllSolutions(const int sz);
  ~CostasNotEqualAllSolutions();

  virtual void run();
};

const int CostasNotEqualAllSolutions::num_sol[] = {0,0,0,2,12,40,116,200,444,760,2160,0,0};

class Reset : public UnitTest {

public:
  
  Reset();
  ~Reset();

  virtual void run();
};

class Pigeons : public UnitTest {

public:
  
  static const int num_bts[13];

  int size;
  
  Pigeons(const int sz);
  ~Pigeons();

  virtual void run();
};

const int Pigeons::num_bts[] = {0,0,0,0,0,0,119,719,5039,40319,362879,3628799,0};


class BoolPigeons : public UnitTest {

public:

  static const int num_bts[13];
  
  int var_type;
  int size;
  
  BoolPigeons(const int sz, const int vt);
  ~BoolPigeons();

  virtual void run();
};

const int BoolPigeons::num_bts[] = {0,0,0,0,0,0,119,719,5039,40319,362879,3628799,0};

class VarStackDynamicTest : public UnitTest {

public:
  
  VarStackDynamicTest();
  ~VarStackDynamicTest();

  virtual void run();
};

class ModelTest : public UnitTest {

public:
  
  ModelTest();
  ~ModelTest();

  virtual void run();
};

class RewriteTest1 : public UnitTest {

public:
  
  RewriteTest1();
  ~RewriteTest1();

  virtual void run();
};

class SubsetTest : public UnitTest {

public:
  
  SubsetTest();
  ~SubsetTest();
  
  void run1();
  void run2();
  void run3();
  void run4();
  void run5();
  void run6();
  void run7();
  void run8();
  virtual void run();
};

class EqualSetTest : public UnitTest {

public:
  
  EqualSetTest();
  ~EqualSetTest();
  
  void run1();
  void run2();
  void run3();
  void run4();
  void run5();
  void run6();
  virtual void run();
};

class MemberTest : public UnitTest {

public:
  
  MemberTest();
  ~MemberTest();

  void run1();
  void run2();
  virtual void run();
};


class CardTest : public UnitTest {

public:
  
  CardTest();
  ~CardTest();

  //void run1();
  //void run2();
  virtual void run();
};


class DivTest : public UnitTest {

public:
  
  DivTest();
  ~DivTest();

  virtual void run();
};

class SatTest : public UnitTest {

public:
  
  SatTest();
  ~SatTest();

  virtual void run();
};

class IntersectionTest : public UnitTest {

public:
  
  IntersectionTest();
  ~IntersectionTest();

  void run1();
  void run2();
  void run3();
  virtual void run();
};

class UnionTest : public UnitTest {

public:
  
  UnionTest();
  ~UnionTest();

  void run1();
  void run2();
  void run3();
  virtual void run();
};

class SymmetricDifferenceTest : public UnitTest {

public:
  
  SymmetricDifferenceTest();
  ~SymmetricDifferenceTest();

  void run1();
  void run2();
  void run3();
  virtual void run();
};

class SetDifferenceTest : public UnitTest {

public:
  
  SetDifferenceTest();
  ~SetDifferenceTest();

  void run1();
  void run2();
  void run3();
  virtual void run();
};


class LexTest : public UnitTest {

public:
  
  LexTest();
  ~LexTest();

  virtual void run();
};


class WeightedSumTest : public UnitTest {

public:
  
  WeightedSumTest();
  ~WeightedSumTest();

void run1();
void run2();
void run3();
void run4();

  virtual void run();
};

class ElementTest : public UnitTest {

public:
  
  ElementTest();
  ~ElementTest();

  void run1();
  void run2();
//void run3();
//void run4();

  virtual void run();
};

class MinMaxTest : public UnitTest {

public:
  
  MinMaxTest();
  ~MinMaxTest();

  void run1();
  void run2();

  virtual void run();
};


class OpshopTest : public UnitTest {

public:
  
  OpshopTest();
  ~OpshopTest();

  void run1();
  void run2();
  void run3();

  virtual void run();
};

class RewriteTest : public UnitTest {

public:
  
  RewriteTest();
  ~RewriteTest();

  virtual void run();
};




template<class CON_TYPE>
class ConChecker {

public:

  std::string name;
  
  bool AC;
  bool BC;
  // AC = true, BC = false: = AC
  // AC = false, BC = true: = BC
  // AC = true, BC = true: BC <= P <= AC

  int seed;
  int arity;
  int domsize;
  int nbool;
  
  // whether we force the checker propagator to trigger on its own changes
  // (should be set to true if the propagator is idempotent, but the checker propagator is not)
  bool trigger_self;

  Solver s1;
  Solver s2;
  Solver s3;

  VarArray scope1;
  VarArray scope2;
  VarArray scope3;

  Constraint con1;// = new CON_TYPE(scope1);
  Constraint con2;// = new CON_TYPE(scope2);
  Constraint con3;// = new CON_TYPE(scope3);

  ConChecker(//CON_TYPE *cons, 
	     std::string nm,
	     const bool ac,
	     const bool bc,
	     const int s,
	     const int a, 
	     const int d, 
	     const int n=0,
	     const bool t=false) : name(nm), AC(ac), BC(bc), seed(s), arity(a), domsize(d), nbool(n), trigger_self(t) {
    
    VarArray X1(arity-nbool, -domsize/2, +domsize/2);
    VarArray b1(nbool, 0, 1);
    
    VarArray X2(arity-nbool, -domsize/2, +domsize/2);
    VarArray b2(nbool, 0, 1);

    VarArray X3(arity-nbool, -domsize/2, +domsize/2);
    VarArray b3(nbool, 0, 1);
    
    for(int i=0; i<arity-nbool; ++i) {
      scope1.add(X1[i]);
      scope2.add(X2[i]);
      scope3.add(X3[i]);
    }
    
    for(int i=0; i<nbool; ++i) {
      scope1.add(b1[i]);
      scope2.add(b2[i]);
      scope3.add(b3[i]);
    }
    
    
    con1 = Constraint(new CON_TYPE(scope1));
    con2 = Constraint(new CON_TYPE(scope2));
    //con2.set_idempotent(false);
    con3 = Constraint(new CON_TYPE(scope3));
    //con3.set_idempotent(false);


    usrand(seed);
  }
  


  void init() {

    s1.add( (con1) );
    //s1.rewrite();
    s1.consolidate();
    
    s2.add( (con2) );
    //s2.rewrite();
    s2.consolidate();
    
    s3.add( (con3) );
    //s3.rewrite();
    s3.consolidate();
    
    //con2.set_idempotent(false);
    //con3.set_idempotent(false);


    // std::cout << s1.constraints[0] << (s1.constraints[0].idempotent() ? " is " : " isn't " ) << "idempotent\n";
    // std::cout << s2.constraints[0] << (s2.constraints[0].idempotent() ? " is " : " isn't " ) << "idempotent\n";
    // std::cout << s3.constraints[0] << (s3.constraints[0].idempotent() ? " is " : " isn't " ) << "idempotent\n";
    

    // exit(1);
  }


  void check(PropagationOutcome wiped1,
	     PropagationOutcome wiped2,
	     PropagationOutcome wiped3,
	     const int iteration) {


#ifdef _DEBUG_CHECKER
    if(_DEBUG_CHECKER) {
    Variable tmp;
	      std::cout << "=> ";
	      for(int j=0; j<arity; ++j) {
		std::cout << scope1[j].get_var().get_domain() << " ";
	      }
	      if(AC) {
		std::cout << std::endl << "=> ";
		for(int j=0; j<arity; ++j) {
		  std::cout << scope2[j].get_var().get_domain() << " " ;
		}	
	      }
	      if(BC) {
		std::cout << std::endl << "=> ";
		for(int j=0; j<arity; ++j) {
		  tmp = scope3[j].get_var();
		  std::cout << tmp.get_domain() << " " ;
		}	
	      }
	      std::cout << std::endl;
    }
#endif
	      
	      //test if there is a fail
	      if( AC && BC ) {
		if(IS_OK(wiped2) && !IS_OK(wiped1)) {
		  cout << "at iteration " << iteration 
		       << " Error - inconsistency in propag of \"" << con1 << "\"! (result propag </= AC)" 
		       << endl << " propagator: fail " << endl ;
		  std::cout << std::endl << " generic AC: "; 
		  for(int j=0; j<arity; ++j) {
		    std::cout << scope2[j].get_var().get_domain() << " " ;
		  }
		  std::cout << std::endl;
		  exit(1);
		}
		if(IS_OK(wiped1) && !IS_OK(wiped3)) {
		  cout << "at iteration " << iteration 
		       << " Error - inconsistency in propag of \"" << con1 << "\"! (result BC </= propag)" 
		       << endl << " propagator: " ;
		  for(int j=0; j<arity; ++j) {
		    std::cout << scope1[j].get_domain() << " " ;
		  }
		  std::cout << std::endl << " generic BC: fail" << endl;
		  exit(1);
		}
	      } else if(AC) {
		if(IS_OK(wiped1) != IS_OK(wiped2)) {
		  cout << "at iteration " << iteration 
		       << " Error - inconsistency in propag of \"" << con1 << "\" (result propag =/= AC)!" 
		       << endl << " propagator: " ;
		  if(IS_OK(wiped1)) {
		    for(int j=0; j<arity; ++j) {
		      std::cout << scope1[j].get_var().get_domain() << " " ;
		    }
		  } else std::cout << "fail" ;
		
		  std::cout << std::endl << " generic AC: "; 
		  if(IS_OK(wiped2)) {
		    for(int j=0; j<arity; ++j) {
		      std::cout << scope2[j].get_var().get_domain() << " " ;
		    }
		  } else std::cout << "fail" ;
		  std::cout << endl;
		  exit(1);
		}
	      } else if(BC) {
		if(IS_OK(wiped1) != IS_OK(wiped3)) {
		  cout << "at iteration " << iteration 
		       << " Error - inconsistency in propag of \"" << con1 << "\" (result propag =/= BC)!" 
		       << endl << " propagator: " ;
		  if(IS_OK(wiped1)) {
		    for(int j=0; j<arity; ++j) {
		      std::cout << scope1[j].get_var().get_domain() << " " ;
		    }
		  } else std::cout << "fail" ;
		
		  std::cout << std::endl << " generic BC: "; 
		  if(IS_OK(wiped3)) {
		    for(int j=0; j<arity; ++j) {
		      std::cout << scope3[j].get_var().get_domain() << " " ;
		    }
		  } else std::cout << "fail"  ;
		  std::cout << endl;
		  exit(1);
		}
	      }
	    
	      //std::cout << std::endl;
	      if(IS_OK(wiped1) && (!AC || IS_OK(wiped2)) && (!BC || IS_OK(wiped3))) {
		for(int i=0; i<arity; ++i) {
		  int vnxt1 = scope1[i].get_var().get_min();
		  int vnxt2 = scope2[i].get_var().get_min();
		  int vnxt3 = scope3[i].get_var().get_min();
		  int val1 = vnxt1-1;
		  int val2 = vnxt2-1;
		  int val3 = vnxt3-1;
		
		  BitSet d1(vnxt1, scope1[i].get_var().get_max(), BitSet::empt);
		  BitSet d2(vnxt2, scope2[i].get_var().get_max(), BitSet::empt);
		  BitSet d3(vnxt3, scope3[i].get_var().get_max(), BitSet::empt);
		

		  //std::cout << std::endl << scope1[i].get_var().get_domain()  << " "  << scope1[i].get_var().get_min()  << std::endl;
		  while(val1<vnxt1) {
		    val1 = vnxt1;
		    vnxt1 = scope1[i].get_var().next(val1);
		    d1.add(val1);
		  }
		
		  if(AC) {
		    while(val2<vnxt2) {
		      val2 = vnxt2;
		      vnxt2 = scope2[i].get_var().next(val2);
		      d2.add(val2);
		    }
		  }
		
		  if(BC) {
		    //std::cout << scope3[i].get_var().get_domain() << " " << scope3[i].get_var().get_min() << std::endl;
		    while(val3<vnxt3) {
		      val3 = vnxt3;
		      vnxt3 = scope3[i].get_var().next(val3);
		      //std::cout << " " << val3;
		      d3.add(val3);
		    }
		    //std::cout << std::endl << d3 << std::endl;
		  
		  }
		
		  if( AC && BC ) {
		    if(!d1.includes(d2)) {
		      cout << "at iteration " << iteration 
			   << " Error - inconsistency in propag of \"" << con1 << "\"! (pruning propag </= AC)" 
			   << endl << " propagator: " ;
		      for(int j=0; j<arity; ++j)
			std::cout << scope1[j].get_var() << " " << scope1[j].get_var().get_domain() << ( i==j ? "* " : " ") ;
		    
		      std::cout << std::endl << " generic AC: ";		  
		      for(int j=0; j<arity; ++j)
			std::cout << scope2[j].get_var() << " " << scope2[j].get_var().get_domain() << ( i==j ? "* " : " ") ;
		      std::cout << std::endl << std::endl;
		      exit(1);
		    }
		    if(!d3.includes(d1)) {
		      cout << "at iteration " << iteration 
			   << " Error - inconsistency in propag of \"" << con1 << "\"! (pruning BC </= propag)" 
			   << endl << " propagator: " ;
		      for(int j=0; j<arity; ++j)
			std::cout << scope1[j].get_var() << " " << scope1[j].get_var().get_domain() << ( i==j ? "* " : " ") ;
		    
		      std::cout << std::endl << " generic BC: ";		  
		      for(int j=0; j<arity; ++j)
			std::cout << scope3[j].get_var() << " " << scope3[j].get_var().get_domain() << ( i==j ? "* " : " ") ;
		      std::cout << std::endl << std::endl;
		      exit(1);
		    }
		  } else if(AC) {
		    if(d1 != d2) {
		      cout << "at iteration " << iteration 
			   << " Error - inconsistency in propag of \"" << con1 << "\"! (pruning AC =/= propag)" 
			   << endl << " propagator: " ;
		      for(int j=0; j<arity; ++j)
			std::cout << scope1[j].get_var() << " " << scope1[j].get_var().get_domain() << ( i==j ? "* " : " ") ;
		  
		      std::cout << std::endl << " generic AC: ";		  
		      for(int j=0; j<arity; ++j)
			std::cout << scope2[j].get_var() << " " << scope2[j].get_var().get_domain() << ( i==j ? "* " : " ") ;
		      std::cout << std::endl << std::endl;
		      exit(1);		  
		    }
		  } else if(BC) {
		    if(d1 != d3) {
		      cout << d1 << " - " << d3 << std::endl;
		      cout << "at iteration " << iteration 
			   << " Error - inconsistency in propag of \"" << con1 << "\"! (pruning BC =/= propag)" 
			   << endl << " propagator: " ;
		      for(int j=0; j<arity; ++j)
			std::cout << scope1[j].get_var() << " " << scope1[j].get_var().get_domain() << ( i==j ? "* " : " ") ;
		  
		      std::cout << std::endl << " generic BC: ";		  
		      for(int j=0; j<arity; ++j)
			std::cout << scope3[j].get_var() << " " << scope3[j].get_var().get_domain() << ( i==j ? "* " : " ") ;
		      std::cout << std::endl << std::endl;
		      exit(1);
		    }
		  } 
		}
	      }
  }


  void run(int n_iterations=1) {

    std::cout << ".";
    std::cout.flush();
    cur_iteration = 0;
    while(cur_iteration++ < n_iterations) {      
      int n_reductions = 0;
      
      PropagationOutcome wiped1 = CONSISTENT;
      PropagationOutcome wiped2 = CONSISTENT;
      PropagationOutcome wiped3 = CONSISTENT;
      
      s1.save();
      if(AC) s2.save();
      if(BC) s3.save();
      
#ifdef _DEBUG_CHECKER
    if(_DEBUG_CHECKER) {
      //std::cout << "s1.statistics.num_propagations: " << s1.statistics.num_propagations << std::endl;
      std::cout << std::endl << cur_iteration << " check " << name << " " << con1 << std::endl;
      for(int j=0; j<arity; ++j) {
	std::cout << scope1[j].get_var().get_domain() << " " ;
      }
      if(AC) {
	std::cout << std::endl;
	for(int j=0; j<arity; ++j) {
	  std::cout << scope2[j].get_var().get_domain() << " " ;
	}	
      }
      if(BC) {
	std::cout << std::endl;
	for(int j=0; j<arity; ++j) {
	  std::cout << scope3[j].get_var().get_domain() << " " ;
	}	
      }
      std::cout << std::endl << "initial propagate: " << std::endl;
    }
#endif
      
      
      bool finished = false;


      //std::cout << "\nspecific propagate" << std::endl;
      wiped1 = s1.propagate(con1, true, false);
      if(AC && IS_OK(wiped2)) {
	//std::cout << "\nchecker propagate" << std::endl;
	wiped2 = s2.checker_propagate(con2, true, trigger_self);
      }
      if(BC && IS_OK(wiped3)) {
	//std::cout << "\nbound checker propagate" << std::endl;
	wiped3 = s3.bound_checker_propagate(con3, true, trigger_self);
      }

      check(wiped1, wiped2, wiped3, cur_iteration);

      while(!finished) {
	
	for(int j=0; j<arity; ++j) {
	  
	  bool is_ground = true;
	  for(int i=0; is_ground && i<arity; ++i) {
	    is_ground = scope1[i].get_var().is_ground();
	  }
	  if(is_ground) {
	    finished = true;
	    break;
	  } else if(!scope1[j].get_var().is_ground()) {
	    
	    int type = 0;
	    if(j<(arity-nbool) && randint(3)<1) type = 1;
	    
	    int k = randint(scope1[j].get_var().get_max() - scope1[j].get_var().get_min() + 1) + scope1[j].get_var().get_min();
	    
	    ++n_reductions;
	    
#ifdef _DEBUG_CHECKER
    if(_DEBUG_CHECKER) {
	    std::cout << std::endl << std::endl;;
    }
#endif
	    
	    if(type) {
	      
	      if(randint(2)) {
		
		scope1[j] = scope1[j].get_var();
		scope2[j] = scope2[j].get_var();
		scope3[j] = scope3[j].get_var();

#ifdef _DEBUG_CHECKER
    if(_DEBUG_CHECKER) {
		std::cout << scope1[j].get_var() << " in " << scope1[j].get_var().get_domain() << ">=" << k << std::endl;;
    }
#endif

		if( FAILED(scope1[j].get_var().set_min(k)) ) wiped1 = FAILURE(j);
		if( AC && FAILED(scope2[j].get_var().set_min(k)) ) wiped2 = FAILURE(j);
		if( BC && FAILED(scope3[j].get_var().set_min(k)) ) wiped3 = FAILURE(j);

		scope1[j] = scope1[j].get_var();
		scope2[j] = scope2[j].get_var();
		scope3[j] = scope3[j].get_var();

	      } else {
	    
		scope1[j] = scope1[j].get_var();
		scope2[j] = scope2[j].get_var();
		scope3[j] = scope3[j].get_var();

#ifdef _DEBUG_CHECKER
    if(_DEBUG_CHECKER) {
		std::cout << scope1[j].get_var() << " in " << scope1[j].get_var().get_domain() << "<=" << k << std::endl;;
    }
#endif

		if( FAILED(scope1[j].get_var().set_max(k)) ) wiped1 = FAILURE(j);
		if( AC && FAILED(scope2[j].get_var().set_max(k)) ) wiped2 = FAILURE(j);
		if( BC && FAILED(scope3[j].get_var().set_max(k)) ) wiped3 = FAILURE(j);

		scope1[j] = scope1[j].get_var();
		scope2[j] = scope2[j].get_var();
		scope3[j] = scope3[j].get_var();

	      }
	    } else {

	      if(randint(5)<1) {

		scope1[j] = scope1[j].get_var();
		scope2[j] = scope2[j].get_var();
		scope3[j] = scope3[j].get_var();
	    
#ifdef _DEBUG_CHECKER
    if(_DEBUG_CHECKER) {
		std::cout << scope1[j].get_var() << " in " << scope1[j].get_var().get_domain() << "==" << k << std::endl;;
    }
#endif
	    
		if( FAILED(scope1[j].get_var().set_domain(k)) ) wiped1 = FAILURE(j);
		if( AC && FAILED(scope2[j].get_var().set_domain(k)) ) wiped2 = FAILURE(j);
		if( BC && FAILED(scope3[j].get_var().set_domain(k)) ) wiped3 = FAILURE(j);

		scope1[j] = scope1[j].get_var();
		scope2[j] = scope2[j].get_var();
		scope3[j] = scope3[j].get_var();

	      } else {

		scope1[j] = scope1[j].get_var();
		scope2[j] = scope2[j].get_var();
		scope3[j] = scope3[j].get_var();
	    
#ifdef _DEBUG_CHECKER
    if(_DEBUG_CHECKER) {
		std::cout << scope1[j].get_var() << " in " << scope1[j].get_var().get_domain() << "!=" << k << std::endl;;
    }
#endif
	    
		if( FAILED(scope1[j].get_var().remove(k)) ) wiped1 = FAILURE(j);
		if( AC && FAILED(scope2[j].get_var().remove(k)) ) wiped2 = FAILURE(j);
		if( BC && FAILED(scope3[j].get_var().remove(k)) ) wiped3 = FAILURE(j);

		scope1[j] = scope1[j].get_var();
		scope2[j] = scope2[j].get_var();
		scope3[j] = scope3[j].get_var();

	      }
	    }

#ifdef _DEBUG_CHECKER
    if(_DEBUG_CHECKER) {
	    Variable tmp;
	    for(int j=0; j<arity; ++j) {
	      tmp =  scope1[j].get_var();
	      std::cout // <<tmp << " in " 
			<<tmp.get_domain() << " " ;//<< " |" <<tmp.get_size() << "| + ";
	      //if(tmp.domain_type == BITSET_VAR) std::cout <<((VariableBitmap*)(tmp.variable))->trail_ << " " ;
	    }
	    std::cout << std::endl;
    }
#endif

	    if( IS_OK(wiped1) ) {

	      //std::cout << "\nspecific propagate" << std::endl;
	      wiped1 = s1.propagate(con1, false, false);

	      if(AC && IS_OK(wiped2)) {

		//std::cout << "\nchecker propagate" << std::endl;
		wiped2 = s2.checker_propagate(con2, true, trigger_self);
		
	      }

	      if(BC && IS_OK(wiped3)) {
		
		//std::cout << "\nbound checker propagate" << std::endl;
		wiped3 = s3.bound_checker_propagate(con3, true, trigger_self);

	      }
	    
	      check(wiped1, wiped2, wiped3, cur_iteration);
	      
	    } else {
	      s1.active_variables.clear();
	      s1.active_constraints.clear();
	      s2.active_variables.clear();
	      s2.active_constraints.clear();
	      s3.active_variables.clear();
	      s3.active_constraints.clear();
	      finished = true;
	      break;
	    }
	  }
	}
      }
    

      s1.restore(0);
      s1.wiped_idx = CONSISTENT;
      if(AC) {
	s2.restore(0);
	s2.wiped_idx = CONSISTENT;
      }
      if(BC) {
	s3.restore(0);
	s3.wiped_idx = CONSISTENT;
      }

    }

  }
  
  virtual ~ConChecker() {};
};

class CheckerTest : public UnitTest {
  
public:
  
  CheckerTest() : UnitTest() {}
  virtual ~CheckerTest() {}

 
  void run() {
    
    if(Verbosity) cout << "Run Checker test: "; 


    ConChecker< PredicateLess > cc1a("<=", false,true,12345,3,20,1);
    cc1a.init();
    cc1a.run(1000);


    ConChecker< PredicateLess > cc1b("<=", false,true,12345,3,20,1);
    ((PredicateLess*)(cc1b.con1.propagator))->offset = 3;
    ((PredicateLess*)(cc1b.con2.propagator))->offset = 3;
    ((PredicateLess*)(cc1b.con3.propagator))->offset = 3;
    cc1b.init();
    cc1b.run(1000);


    ConChecker< PredicateLess > cc1c("<=", false,true,12345,3,20,1);
    ((PredicateLess*)(cc1c.con1.propagator))->offset = -3;
    ((PredicateLess*)(cc1c.con2.propagator))->offset = -3;
    ((PredicateLess*)(cc1c.con3.propagator))->offset = -3;   
    cc1c.init();
    cc1c.run(1000);

    ConChecker< PredicateEqual > cc2a("==", true,false,12345,3,20,1);
    cc2a.init();
    cc2a.run(1000);

    ConChecker< PredicateEqual > cc2b("!=", true,false,12345,3,20,1);
    ((PredicateEqual*)(cc2b.con1.propagator))->spin = 0;
    ((PredicateEqual*)(cc2b.con2.propagator))->spin = 0;
    ((PredicateEqual*)(cc2b.con3.propagator))->spin = 0;
    cc2b.init();
    cc2b.run(1000);


    ConChecker< PredicateLowerBound > cc3a(">=-3", false,true,12345,2,20,1);
    ((PredicateLowerBound*)(cc3a.con1.propagator))->bound = -3;
    ((PredicateLowerBound*)(cc3a.con2.propagator))->bound = -3;
    ((PredicateLowerBound*)(cc3a.con3.propagator))->bound = -3;
    cc3a.init();
    cc3a.run(1000);

    ConChecker< PredicateLowerBound > cc3b(">=3", false,true,12345,2,20,1);
    ((PredicateLowerBound*)(cc3b.con1.propagator))->bound = 3;
    ((PredicateLowerBound*)(cc3b.con2.propagator))->bound = 3;
    ((PredicateLowerBound*)(cc3b.con3.propagator))->bound = 3;
    cc3b.init();
    cc3b.run(1000);

    ConChecker< PredicateUpperBound > cc4a("<=-3", false,true,12345,2,20,1);
    ((PredicateUpperBound*)(cc4a.con1.propagator))->bound = -3;
    ((PredicateUpperBound*)(cc4a.con2.propagator))->bound = -3;
    ((PredicateUpperBound*)(cc4a.con3.propagator))->bound = -3;
    cc4a.init();
    cc4a.run(1000);

    ConChecker< PredicateUpperBound > cc4b("<=3", false,true,12345,2,20,1);
    ((PredicateUpperBound*)(cc4b.con1.propagator))->bound = 3;
    ((PredicateUpperBound*)(cc4b.con2.propagator))->bound = 3;
    ((PredicateUpperBound*)(cc4b.con3.propagator))->bound = 3;
    cc4b.init();
    cc4b.run(1000);

    ConChecker< PredicateConstantEqual > cc5a("==0", true,false,12345,2,20,1);
    cc5a.init();
    cc5a.run(1000);

    ConChecker< PredicateConstantEqual > cc5b("!=0", true,false,12345,2,20,1);
    ((PredicateConstantEqual*)(cc5b.con1.propagator))->spin = 0;
    ((PredicateConstantEqual*)(cc5b.con2.propagator))->spin = 0;
    ((PredicateConstantEqual*)(cc5b.con3.propagator))->spin = 0;
    cc5b.init();
    cc5b.run(1000);

    ConChecker< PredicateAnd > cc6a("and", true,false,12345,3,1,3);
    cc6a.init();
    cc6a.run(1000);

    ConChecker< PredicateOr > cc7a("or", true,false,12345,3,1,3);
    cc7a.init();
    cc7a.run(1000);

    ConChecker< PredicateAdd > cc8a("+", false,true,12345,3,20,0);
    cc8a.init();
    cc8a.run(1000);



    // ConChecker< PredicateSub > cc8a("-", false,true,12345,3,20,0);
    // cc8a.run(100);
   

    ConChecker< PredicateWeightedSum > cc9a("sum(1)", false,true,12345,5,10,0);
    ((PredicateWeightedSum*)(cc9a.con1.propagator))->lower_bound = 10;
    ((PredicateWeightedSum*)(cc9a.con1.propagator))->upper_bound = 10;
    ((PredicateWeightedSum*)(cc9a.con2.propagator))->lower_bound = 10;
    ((PredicateWeightedSum*)(cc9a.con2.propagator))->upper_bound = 10;
    ((PredicateWeightedSum*)(cc9a.con3.propagator))->lower_bound = 10;
    ((PredicateWeightedSum*)(cc9a.con3.propagator))->upper_bound = 10;
    cc9a.init();
    cc9a.run(100);



    ConChecker< PredicateWeightedSum > cc9b("sum(+)", false,true,12345,5,10,0);
    for(int i=2; i<5; ++i) {
      ((PredicateWeightedSum*)(cc9b.con1.propagator))->weight[i] = randint(4)+2;
      ((PredicateWeightedSum*)(cc9b.con2.propagator))->weight[i] = ((PredicateWeightedSum*)(cc9b.con1.propagator))->weight[i];
      ((PredicateWeightedSum*)(cc9b.con3.propagator))->weight[i] = ((PredicateWeightedSum*)(cc9b.con1.propagator))->weight[i];
      // if(!(cc9b.con1->weight[i]%2)) {
      // 	cc9b.con1->unknown_parity.reversible_remove(i);
      // 	cc9b.con2->unknown_parity.reversible_remove(i);
      // 	cc9b.con3->unknown_parity.reversible_remove(i);
      // }
    }
    ((PredicateWeightedSum*)(cc9b.con1.propagator))->lower_bound = 15;
    ((PredicateWeightedSum*)(cc9b.con1.propagator))->upper_bound = 15;
    ((PredicateWeightedSum*)(cc9b.con2.propagator))->lower_bound = 15;
    ((PredicateWeightedSum*)(cc9b.con2.propagator))->upper_bound = 15;
    ((PredicateWeightedSum*)(cc9b.con3.propagator))->lower_bound = 15;
    ((PredicateWeightedSum*)(cc9b.con3.propagator))->upper_bound = 15;
    ((PredicateWeightedSum*)(cc9b.con1.propagator))->wpos = 2;
    ((PredicateWeightedSum*)(cc9b.con2.propagator))->wpos = 2;
    ((PredicateWeightedSum*)(cc9b.con3.propagator))->wpos = 2;
    ///std::cout << std::endl<< std::endl<< cc9b.con1 << std::endl;
    //std::cout << cc9b.con2 << std::endl;
    cc9b.init();
    cc9b.run(100);

    ConChecker< PredicateWeightedSum > cc9c("sum(-)", false,true,12345,5,10,0);
    for(int i=3; i<5; ++i) {
      ((PredicateWeightedSum*)(cc9c.con1.propagator))->weight[i] = -(randint(4)+1);
      ((PredicateWeightedSum*)(cc9c.con2.propagator))->weight[i] = ((PredicateWeightedSum*)(cc9c.con1.propagator))->weight[i];
      ((PredicateWeightedSum*)(cc9c.con3.propagator))->weight[i] = ((PredicateWeightedSum*)(cc9c.con1.propagator))->weight[i];
      // if(!(cc9c.con1->weight[i]%2)) {
      // 	cc9c.con1->unknown_parity.reversible_remove(i);
      // 	cc9c.con2->unknown_parity.reversible_remove(i);
      // 	cc9c.con3->unknown_parity.reversible_remove(i);
      // }
    }
    ((PredicateWeightedSum*)(cc9c.con1.propagator))->wpos = 3;
    ((PredicateWeightedSum*)(cc9c.con2.propagator))->wpos = 3;
    ((PredicateWeightedSum*)(cc9c.con3.propagator))->wpos = 3;
    ((PredicateWeightedSum*)(cc9c.con1.propagator))->wneg = 3;
    ((PredicateWeightedSum*)(cc9c.con2.propagator))->wneg = 3;
    ((PredicateWeightedSum*)(cc9c.con3.propagator))->wneg = 3;
    ((PredicateWeightedSum*)(cc9c.con1.propagator))->lower_bound = -10;
    ((PredicateWeightedSum*)(cc9c.con1.propagator))->upper_bound = -10;
    ((PredicateWeightedSum*)(cc9c.con2.propagator))->lower_bound = -10;
    ((PredicateWeightedSum*)(cc9c.con2.propagator))->upper_bound = -10;
    ((PredicateWeightedSum*)(cc9c.con3.propagator))->lower_bound = -10;
    ((PredicateWeightedSum*)(cc9c.con3.propagator))->upper_bound = -10;
    cc9c.init();
    cc9c.run(100);


    ConChecker< PredicateWeightedSum > cc9d("sum(+/-)", false,true,12345,5,10,0);
    for(int i=1; i<3; ++i) {
      ((PredicateWeightedSum*)(cc9d.con1.propagator))->weight[i] = randint(4)+2;
      ((PredicateWeightedSum*)(cc9d.con2.propagator))->weight[i] = ((PredicateWeightedSum*)(cc9d.con1.propagator))->weight[i];
      ((PredicateWeightedSum*)(cc9d.con3.propagator))->weight[i] = ((PredicateWeightedSum*)(cc9d.con1.propagator))->weight[i];
      // if(!(cc9d.con1->weight[i]%2)) {
      // 	cc9d.con1->unknown_parity.reversible_remove(i);
      // 	cc9d.con2->unknown_parity.reversible_remove(i);
      // 	cc9d.con3->unknown_parity.reversible_remove(i);
      // }
    }
    for(int i=3; i<5; ++i) {
      ((PredicateWeightedSum*)(cc9d.con1.propagator))->weight[i] = -(randint(4)+1);
      ((PredicateWeightedSum*)(cc9d.con2.propagator))->weight[i] = ((PredicateWeightedSum*)(cc9d.con1.propagator))->weight[i];
      ((PredicateWeightedSum*)(cc9d.con3.propagator))->weight[i] = ((PredicateWeightedSum*)(cc9d.con1.propagator))->weight[i];
      // if(!(cc9d.con1.propagator))->weight[i]%2)) {
      // 	cc9d.con1.propagator))->unknown_parity.reversible_remove(i);
      // 	cc9d.con2.propagator))->unknown_parity.reversible_remove(i);
      // 	cc9d.con3.propagator))->unknown_parity.reversible_remove(i);
      // }
    }
    ((PredicateWeightedSum*)(cc9d.con1.propagator))->wpos = 1;
    ((PredicateWeightedSum*)(cc9d.con2.propagator))->wpos = 1;
    ((PredicateWeightedSum*)(cc9d.con3.propagator))->wpos = 1;
    ((PredicateWeightedSum*)(cc9d.con1.propagator))->wneg = 3;
    ((PredicateWeightedSum*)(cc9d.con2.propagator))->wneg = 3;
    ((PredicateWeightedSum*)(cc9d.con3.propagator))->wneg = 3;
    cc9d.init();
    cc9d.run(100);

    ConChecker< PredicateWeightedSum > cc9e("sum(1)", false,true,12345,5,10,0);
    ((PredicateWeightedSum*)(cc9e.con1.propagator))->lower_bound = 5;
    ((PredicateWeightedSum*)(cc9e.con1.propagator))->upper_bound = 10;
    ((PredicateWeightedSum*)(cc9e.con2.propagator))->lower_bound = 5;
    ((PredicateWeightedSum*)(cc9e.con2.propagator))->upper_bound = 10;
    ((PredicateWeightedSum*)(cc9e.con3.propagator))->lower_bound = 5;
    ((PredicateWeightedSum*)(cc9e.con3.propagator))->upper_bound = 10;
    cc9e.init();
    cc9e.run(100);

    ConChecker< PredicateWeightedSum > cc9f("sum(+)", false,true,12345,5,10,0);
    for(int i=2; i<5; ++i) {
      ((PredicateWeightedSum*)(cc9f.con1.propagator))->weight[i] = randint(4)+2;
      ((PredicateWeightedSum*)(cc9f.con2.propagator))->weight[i] = ((PredicateWeightedSum*)(cc9f.con1.propagator))->weight[i];
      ((PredicateWeightedSum*)(cc9f.con3.propagator))->weight[i] = ((PredicateWeightedSum*)(cc9f.con1.propagator))->weight[i];
      // if(!(cc9f.con1.propagator))->weight[i]%2)) {
      // 	cc9f.con1.propagator))->unknown_parity.reversible_remove(i);
      // 	cc9f.con2.propagator))->unknown_parity.reversible_remove(i);
      // 	cc9f.con3.propagator))->unknown_parity.reversible_remove(i);
      // }
    }
    ((PredicateWeightedSum*)(cc9f.con1.propagator))->lower_bound = 10;
    ((PredicateWeightedSum*)(cc9f.con1.propagator))->upper_bound = 15;
    ((PredicateWeightedSum*)(cc9f.con2.propagator))->lower_bound = 10;
    ((PredicateWeightedSum*)(cc9f.con2.propagator))->upper_bound = 15;
    ((PredicateWeightedSum*)(cc9f.con3.propagator))->lower_bound = 10;
    ((PredicateWeightedSum*)(cc9f.con3.propagator))->upper_bound = 15;
    ((PredicateWeightedSum*)(cc9f.con1.propagator))->wpos = 2;
    ((PredicateWeightedSum*)(cc9f.con2.propagator))->wpos = 2;
    ((PredicateWeightedSum*)(cc9f.con3.propagator))->wpos = 2;
    ///std::cout << std::endl<< std::endl<< cc9b.con1 << std::endl;
    //std::cout << cc9b.con2 << std::endl;
    cc9f.init();
    cc9f.run(100);

    ConChecker< PredicateWeightedSum > cc9g("sum(-)", false,true,12345,5,10,0);
    for(int i=3; i<5; ++i) {
      ((PredicateWeightedSum*)(cc9g.con1.propagator))->weight[i] = -(randint(4)+1);
      ((PredicateWeightedSum*)(cc9g.con2.propagator))->weight[i] = ((PredicateWeightedSum*)(cc9g.con1.propagator))->weight[i];
      ((PredicateWeightedSum*)(cc9g.con3.propagator))->weight[i] = ((PredicateWeightedSum*)(cc9g.con1.propagator))->weight[i];
      // if(!(cc9g.con1.propagator))->weight[i]%2)) {
      // 	cc9g.con1.propagator))->unknown_parity.reversible_remove(i);
      // 	cc9g.con2.propagator))->unknown_parity.reversible_remove(i);
      // 	cc9g.con3.propagator))->unknown_parity.reversible_remove(i);
      // }
    }
    ((PredicateWeightedSum*)(cc9g.con1.propagator))->wpos = 3;
    ((PredicateWeightedSum*)(cc9g.con2.propagator))->wpos = 3;
    ((PredicateWeightedSum*)(cc9g.con3.propagator))->wpos = 3;
    ((PredicateWeightedSum*)(cc9g.con1.propagator))->wneg = 3;
    ((PredicateWeightedSum*)(cc9g.con2.propagator))->wneg = 3;
    ((PredicateWeightedSum*)(cc9g.con3.propagator))->wneg = 3;
    ((PredicateWeightedSum*)(cc9g.con1.propagator))->lower_bound = -9;
    ((PredicateWeightedSum*)(cc9g.con1.propagator))->upper_bound = -9;
    ((PredicateWeightedSum*)(cc9g.con2.propagator))->lower_bound = -9;
    ((PredicateWeightedSum*)(cc9g.con2.propagator))->upper_bound = -9;
    ((PredicateWeightedSum*)(cc9g.con3.propagator))->lower_bound = -9;
    ((PredicateWeightedSum*)(cc9g.con3.propagator))->upper_bound = -9;
    cc9g.init();
    cc9g.run(100);


    ConChecker< PredicateWeightedSum > cc9h("sum(-)", false,true,12345,5,10,0);
    for(int i=3; i<5; ++i) {
      ((PredicateWeightedSum*)(cc9h.con1.propagator))->weight[i] = -(randint(4)+1);
      ((PredicateWeightedSum*)(cc9h.con2.propagator))->weight[i] = ((PredicateWeightedSum*)(cc9h.con1.propagator))->weight[i];
      ((PredicateWeightedSum*)(cc9h.con3.propagator))->weight[i] = ((PredicateWeightedSum*)(cc9h.con1.propagator))->weight[i];
      // if(!(cc9h.con1.propagator))->weight[i]%2)) {
      // 	cc9h.con1.propagator))->unknown_parity.reversible_remove(i);
      // 	cc9h.con2.propagator))->unknown_parity.reversible_remove(i);
      // 	cc9h.con3.propagator))->unknown_parity.reversible_remove(i);
      // }
    }
    ((PredicateWeightedSum*)(cc9h.con1.propagator))->wpos = 3;
    ((PredicateWeightedSum*)(cc9h.con2.propagator))->wpos = 3;
    ((PredicateWeightedSum*)(cc9h.con3.propagator))->wpos = 3;
    ((PredicateWeightedSum*)(cc9h.con1.propagator))->wneg = 3;
    ((PredicateWeightedSum*)(cc9h.con2.propagator))->wneg = 3;
    ((PredicateWeightedSum*)(cc9h.con3.propagator))->wneg = 3;
    ((PredicateWeightedSum*)(cc9h.con1.propagator))->lower_bound = -9;
    ((PredicateWeightedSum*)(cc9h.con1.propagator))->upper_bound = -3;
    ((PredicateWeightedSum*)(cc9h.con2.propagator))->lower_bound = -9;
    ((PredicateWeightedSum*)(cc9h.con2.propagator))->upper_bound = -3;
    ((PredicateWeightedSum*)(cc9h.con3.propagator))->lower_bound = -9;
    ((PredicateWeightedSum*)(cc9h.con3.propagator))->upper_bound = -3;
    cc9h.init();
    cc9h.run(100);

    ConChecker< PredicateWeightedSum > cc9j("sum(+/-)", false,true,12345,5,10,0);
    for(int i=1; i<3; ++i) {
      ((PredicateWeightedSum*)(cc9j.con1.propagator))->weight[i] = randint(4)+2;
      ((PredicateWeightedSum*)(cc9j.con2.propagator))->weight[i] = ((PredicateWeightedSum*)(cc9j.con1.propagator))->weight[i];
      ((PredicateWeightedSum*)(cc9j.con3.propagator))->weight[i] = ((PredicateWeightedSum*)(cc9j.con1.propagator))->weight[i];
      // if(!(cc9j.con1.propagator))->weight[i]%2)) {
      // 	cc9j.con1.propagator))->unknown_parity.reversible_remove(i);
      // 	cc9j.con2.propagator))->unknown_parity.reversible_remove(i);
      // 	cc9j.con3.propagator))->unknown_parity.reversible_remove(i);
      // }
    }
    for(int i=3; i<5; ++i) {
      ((PredicateWeightedSum*)(cc9j.con1.propagator))->weight[i] = -(randint(4)+1);
      ((PredicateWeightedSum*)(cc9j.con2.propagator))->weight[i] = ((PredicateWeightedSum*)(cc9j.con1.propagator))->weight[i];
      ((PredicateWeightedSum*)(cc9j.con3.propagator))->weight[i] = ((PredicateWeightedSum*)(cc9j.con1.propagator))->weight[i];
      // if(!(cc9j.con1.propagator))->weight[i]%2)) {
      // 	cc9j.con1.propagator))->unknown_parity.reversible_remove(i);
      // 	cc9j.con2.propagator))->unknown_parity.reversible_remove(i);
      // 	cc9j.con3.propagator))->unknown_parity.reversible_remove(i);
      // }
    }
    ((PredicateWeightedSum*)(cc9j.con1.propagator))->wpos = 1;
    ((PredicateWeightedSum*)(cc9j.con2.propagator))->wpos = 1;
    ((PredicateWeightedSum*)(cc9j.con3.propagator))->wpos = 1;
    ((PredicateWeightedSum*)(cc9j.con1.propagator))->wneg = 3;
    ((PredicateWeightedSum*)(cc9j.con2.propagator))->wneg = 3;
    ((PredicateWeightedSum*)(cc9j.con3.propagator))->wneg = 3;
    ((PredicateWeightedSum*)(cc9j.con1.propagator))->lower_bound = -1;
    ((PredicateWeightedSum*)(cc9j.con1.propagator))->upper_bound = -1;
    ((PredicateWeightedSum*)(cc9j.con2.propagator))->lower_bound = -1;
    ((PredicateWeightedSum*)(cc9j.con2.propagator))->upper_bound = -1;
    ((PredicateWeightedSum*)(cc9j.con3.propagator))->lower_bound = -1;
    ((PredicateWeightedSum*)(cc9j.con3.propagator))->upper_bound = -1;
    cc9j.init();
    cc9j.run(100);


    ConChecker< PredicateWeightedSum > cc9k("sum(+/-)", false,true,12345,5,10,0);
    for(int i=1; i<3; ++i) {
      ((PredicateWeightedSum*)(cc9k.con1.propagator))->weight[i] = randint(4)+2;
      ((PredicateWeightedSum*)(cc9k.con2.propagator))->weight[i] = ((PredicateWeightedSum*)(cc9k.con1.propagator))->weight[i];
      ((PredicateWeightedSum*)(cc9k.con3.propagator))->weight[i] = ((PredicateWeightedSum*)(cc9k.con1.propagator))->weight[i];
      // if(!(cc9k.con1.propagator))->weight[i]%2)) {
      // 	cc9k.con1.propagator))->unknown_parity.reversible_remove(i);
      // 	cc9k.con2.propagator))->unknown_parity.reversible_remove(i);
      // 	cc9k.con3.propagator))->unknown_parity.reversible_remove(i);
      // }
    }
    for(int i=3; i<5; ++i) {
      ((PredicateWeightedSum*)(cc9k.con1.propagator))->weight[i] = -(randint(4)+1);
      ((PredicateWeightedSum*)(cc9k.con2.propagator))->weight[i] = ((PredicateWeightedSum*)(cc9k.con1.propagator))->weight[i];
      ((PredicateWeightedSum*)(cc9k.con3.propagator))->weight[i] = ((PredicateWeightedSum*)(cc9k.con1.propagator))->weight[i];
      // if(!(cc9k.con1.propagator))->weight[i]%2)) {
      // 	cc9k.con1.propagator))->unknown_parity.reversible_remove(i);
      // 	cc9k.con2.propagator))->unknown_parity.reversible_remove(i);
      // 	cc9k.con3.propagator))->unknown_parity.reversible_remove(i);
      // }
    }
    ((PredicateWeightedSum*)(cc9k.con1.propagator))->wpos = 1;
    ((PredicateWeightedSum*)(cc9k.con2.propagator))->wpos = 1;
    ((PredicateWeightedSum*)(cc9k.con3.propagator))->wpos = 1;
    ((PredicateWeightedSum*)(cc9k.con1.propagator))->wneg = 3;
    ((PredicateWeightedSum*)(cc9k.con2.propagator))->wneg = 3;
    ((PredicateWeightedSum*)(cc9k.con3.propagator))->wneg = 3;
    ((PredicateWeightedSum*)(cc9k.con1.propagator))->lower_bound = -1;
    ((PredicateWeightedSum*)(cc9k.con1.propagator))->upper_bound = 3;
    ((PredicateWeightedSum*)(cc9k.con2.propagator))->lower_bound = -1;
    ((PredicateWeightedSum*)(cc9k.con2.propagator))->upper_bound = 3;
    ((PredicateWeightedSum*)(cc9k.con3.propagator))->lower_bound = -1;
    ((PredicateWeightedSum*)(cc9k.con3.propagator))->upper_bound = 3;
    cc9k.init();
    cc9k.run(100);


    ConChecker< PredicateElement > cc10a("[]", true,false,12345,6,10,0);
    cc10a.init();
    cc10a.run(100);


    /// propagator for MUL does NOT enforce BC
    ConChecker< PredicateMul > cc11a("*", true,true,12345,3,20,0);
    cc11a.init();
    cc11a.run(100);


    ConChecker< PredicateIntervalMember > cc12a("member[]", true,false,12345,2,20,1);
    ((PredicateIntervalMember*)(cc12a.con1.propagator))->spin = 1;
    ((PredicateIntervalMember*)(cc12a.con1.propagator))->lower_bound = -5;
    ((PredicateIntervalMember*)(cc12a.con1.propagator))->upper_bound =  5;
    ((PredicateIntervalMember*)(cc12a.con2.propagator))->spin = 1;
    ((PredicateIntervalMember*)(cc12a.con2.propagator))->lower_bound = -5;
    ((PredicateIntervalMember*)(cc12a.con2.propagator))->upper_bound =  5;
    ((PredicateIntervalMember*)(cc12a.con3.propagator))->spin = 1;
    ((PredicateIntervalMember*)(cc12a.con3.propagator))->lower_bound = -5;
    ((PredicateIntervalMember*)(cc12a.con3.propagator))->upper_bound =  5;
    cc12a.init();
    cc12a.run(1000);

    ConChecker< PredicateIntervalMember > cc12b("notmember[]", true,false,12345,2,20,1);
    ((PredicateIntervalMember*)(cc12b.con1.propagator))->spin = 0;
    ((PredicateIntervalMember*)(cc12b.con1.propagator))->lower_bound = -5;
    ((PredicateIntervalMember*)(cc12b.con1.propagator))->upper_bound =  5;
    ((PredicateIntervalMember*)(cc12b.con2.propagator))->spin = 0;
    ((PredicateIntervalMember*)(cc12b.con2.propagator))->lower_bound = -5;
    ((PredicateIntervalMember*)(cc12b.con2.propagator))->upper_bound =  5;
    ((PredicateIntervalMember*)(cc12b.con3.propagator))->spin = 0;
    ((PredicateIntervalMember*)(cc12b.con3.propagator))->lower_bound = -5;
    ((PredicateIntervalMember*)(cc12b.con3.propagator))->upper_bound =  5;
    cc12b.init();
    cc12b.run(1000);


    BitSet vals(-10,10,BitSet::empt);
    vals.add(-10);
    vals.add(-5);
    vals.add(-1);
    vals.add(0);
    vals.add(3);
    vals.add(7);
    vals.add(8);
    vals.add(9);
    ConChecker< PredicateSetMember > cc13a("member[]", true,false,12345,2,20,1);
    ((PredicateSetMember*)(cc13a.con1.propagator))->spin = 1;
    ((PredicateSetMember*)(cc13a.con1.propagator))->values = vals;
    ((PredicateSetMember*)(cc13a.con2.propagator))->spin = 1;
    ((PredicateSetMember*)(cc13a.con2.propagator))->values = vals;
    ((PredicateSetMember*)(cc13a.con3.propagator))->spin = 1;
    ((PredicateSetMember*)(cc13a.con3.propagator))->values = vals;
    cc13a.init();
    cc13a.run(1000);


    ConChecker< PredicateSetMember > cc13b("notmember[]", true,false,12345,2,20,1);
    ((PredicateSetMember*)(cc13b.con1.propagator))->spin = 0;
    ((PredicateSetMember*)(cc13b.con1.propagator))->values = vals;
    ((PredicateSetMember*)(cc13b.con2.propagator))->spin = 0;
    ((PredicateSetMember*)(cc13b.con2.propagator))->values = vals;
    ((PredicateSetMember*)(cc13b.con3.propagator))->spin = 0;
    ((PredicateSetMember*)(cc13b.con3.propagator))->values = vals;
    cc13b.init();
    cc13b.run(1000);



   ConChecker< PredicateMin > cc14("min", true,true,12345,5,10,0);
   cc14.init();
   cc14.run(100);


   ConChecker< PredicateMax > cc15("max", true,true,12345,5,10,0);
   cc15.init();
   cc15.run(100);
   



    ConChecker< PredicateModConstant > cc16a("%",true,true,12345,2,30,0,true);
    ((PredicateModConstant*)(cc16a.con1.propagator))->modulo = 3;
    ((PredicateModConstant*)(cc16a.con2.propagator))->modulo = 3;
    ((PredicateModConstant*)(cc16a.con3.propagator))->modulo = 3;
    cc16a.init();
    cc16a.run(1000);
    


    ConChecker< PredicateModConstant > cc16b("%",true,true,12345,2,30,0,true);
    ((PredicateModConstant*)(cc16b.con1.propagator))->modulo = 5;
    ((PredicateModConstant*)(cc16b.con2.propagator))->modulo = 5;
    ((PredicateModConstant*)(cc16b.con3.propagator))->modulo = 5;
    cc16b.init();
    cc16b.run(1000);



    ConChecker< PredicateModConstant > cc16c("%",true,true,12345,2,30,0,true);
    ((PredicateModConstant*)(cc16c.con1.propagator))->modulo = 8;
    ((PredicateModConstant*)(cc16c.con2.propagator))->modulo = 8;
    ((PredicateModConstant*)(cc16c.con3.propagator))->modulo = 8;
    cc16c.init();
    cc16c.run(1000);


   ConChecker< PredicateModConstant > cc16na("%",true,true,12345,2,30,0,true);
    ((PredicateModConstant*)(cc16na.con1.propagator))->modulo = -3;
    ((PredicateModConstant*)(cc16na.con2.propagator))->modulo = -3;
    ((PredicateModConstant*)(cc16na.con3.propagator))->modulo = -3;
    cc16na.init();
    cc16na.run(1000);
    

    ConChecker< PredicateModConstant > cc16nb("%",true,true,12345,2,30,0,true);
    ((PredicateModConstant*)(cc16nb.con1.propagator))->modulo = -5;
    ((PredicateModConstant*)(cc16nb.con2.propagator))->modulo = -5;
    ((PredicateModConstant*)(cc16nb.con3.propagator))->modulo = -5;
    cc16nb.init();
    cc16nb.run(1000);


    ConChecker< PredicateModConstant > cc16nc("%",true,true,12345,2,30,0,true);
    ((PredicateModConstant*)(cc16nc.con1.propagator))->modulo = -8;
    ((PredicateModConstant*)(cc16nc.con2.propagator))->modulo = -8;
    ((PredicateModConstant*)(cc16nc.con3.propagator))->modulo = -8;
    cc16nc.init();
    cc16nc.run(1000);
    

    ConChecker< PredicateMod > cc17("%",true,true,12345,3,30,0,true);
    cc17.init();
    cc17.run(100);



   ConChecker< PredicateCModConstant > cc18a("%",true,true,12345,2,30,0,true);
    ((PredicateCModConstant*)(cc18a.con1.propagator))->modulo = 3;
    ((PredicateCModConstant*)(cc18a.con2.propagator))->modulo = 3;
    ((PredicateCModConstant*)(cc18a.con3.propagator))->modulo = 3;
    cc18a.init();
    cc18a.run(10000);
    


    ConChecker< PredicateCModConstant > cc18b("%",true,true,12345,2,30,0,true);
    ((PredicateCModConstant*)(cc18b.con1.propagator))->modulo = 5;
    ((PredicateCModConstant*)(cc18b.con2.propagator))->modulo = 5;
    ((PredicateCModConstant*)(cc18b.con3.propagator))->modulo = 5;
    cc18b.init();
    cc18b.run(1000);



    ConChecker< PredicateCModConstant > cc18c("%",true,true,12345,2,30,0,true);
    ((PredicateCModConstant*)(cc18c.con1.propagator))->modulo = 8;
    ((PredicateCModConstant*)(cc18c.con2.propagator))->modulo = 8;
    ((PredicateCModConstant*)(cc18c.con3.propagator))->modulo = 8;
    cc18c.init();
    cc18c.run(1000);


   ConChecker< PredicateCModConstant > cc18na("%",true,true,12345,2,30,0,true);
    ((PredicateCModConstant*)(cc18na.con1.propagator))->modulo = -3;
    ((PredicateCModConstant*)(cc18na.con2.propagator))->modulo = -3;
    ((PredicateCModConstant*)(cc18na.con3.propagator))->modulo = -3;
    cc18na.init();
    cc18na.run(1000);
    

    ConChecker< PredicateCModConstant > cc18nb("%",true,true,12345,2,30,0,true);
    ((PredicateCModConstant*)(cc18nb.con1.propagator))->modulo = -5;
    ((PredicateCModConstant*)(cc18nb.con2.propagator))->modulo = -5;
    ((PredicateCModConstant*)(cc18nb.con3.propagator))->modulo = -5;
    cc18nb.init();
    cc18nb.run(1000);


    ConChecker< PredicateCModConstant > cc18nc("%",true,true,12345,2,30,0,true);
    ((PredicateCModConstant*)(cc18nc.con1.propagator))->modulo = -8;
    ((PredicateCModConstant*)(cc18nc.con2.propagator))->modulo = -8;
    ((PredicateCModConstant*)(cc18nc.con3.propagator))->modulo = -8;
    cc18nc.init();
    cc18nc.run(1000);
    

    ConChecker< PredicateCMod > cc19("%",true,true,12345,3,30,0,true);
    cc19.init();
    cc19.run(100);
    
    
    ConChecker< PredicateAbs > cc20("abs",true,false,12345,2,100,0,true);
    cc20.init();
    cc20.run(1000);
  

    ConChecker< PredicateDivConstant > cc21pa("/", true,true,12345,2,50,0);
    cc21pa.init();
    cc21pa.run(1000);
    ((PredicateDivConstant*)(cc21pa.con1.propagator))->quotient = 2;
    ((PredicateDivConstant*)(cc21pa.con2.propagator))->quotient = 2;
    ((PredicateDivConstant*)(cc21pa.con3.propagator))->quotient = 2;

    ConChecker< PredicateDivConstant > cc21pb("/", true,true,12345,2,50,0);
    cc21pb.init();
    cc21pb.run(1000);
    ((PredicateDivConstant*)(cc21pb.con1.propagator))->quotient = 3;
    ((PredicateDivConstant*)(cc21pb.con2.propagator))->quotient = 3;
    ((PredicateDivConstant*)(cc21pb.con3.propagator))->quotient = 3;

    ConChecker< PredicateDivConstant > cc21pc("/", true,true,12345,2,50,0);
    cc21pc.init();
    cc21pc.run(1000);
    ((PredicateDivConstant*)(cc21pc.con1.propagator))->quotient = 5;
    ((PredicateDivConstant*)(cc21pc.con2.propagator))->quotient = 5;
    ((PredicateDivConstant*)(cc21pc.con3.propagator))->quotient = 5;

    ConChecker< PredicateDivConstant > cc21pd("/", true,true,12345,2,50,0);
    cc21pd.init();
    cc21pd.run(1000);
    ((PredicateDivConstant*)(cc21pd.con1.propagator))->quotient = 12;
    ((PredicateDivConstant*)(cc21pd.con2.propagator))->quotient = 12;
    ((PredicateDivConstant*)(cc21pd.con3.propagator))->quotient = 12;


    ConChecker< PredicateDivConstant > cc21na("/", true,true,12345,2,50,0);
    cc21na.init();
    cc21na.run(1000);
    ((PredicateDivConstant*)(cc21na.con1.propagator))->quotient = -2;
    ((PredicateDivConstant*)(cc21na.con2.propagator))->quotient = -2;
    ((PredicateDivConstant*)(cc21na.con3.propagator))->quotient = -2;

    ConChecker< PredicateDivConstant > cc21nb("/", true,true,12345,2,50,0);
    cc21nb.init();
    cc21nb.run(1000);
    ((PredicateDivConstant*)(cc21nb.con1.propagator))->quotient = -3;
    ((PredicateDivConstant*)(cc21nb.con2.propagator))->quotient = -3;
    ((PredicateDivConstant*)(cc21nb.con3.propagator))->quotient = -3;

    ConChecker< PredicateDivConstant > cc21nc("/", true,true,12345,2,50,0);
    cc21nc.init();
    cc21nc.run(1000);
    ((PredicateDivConstant*)(cc21nc.con1.propagator))->quotient = -5;
    ((PredicateDivConstant*)(cc21nc.con2.propagator))->quotient = -5;
    ((PredicateDivConstant*)(cc21nc.con3.propagator))->quotient = -5;

    ConChecker< PredicateDivConstant > cc21nd("/", true,true,12345,2,50,0);
    cc21nd.init();
    cc21nd.run(1000);
    ((PredicateDivConstant*)(cc21nd.con1.propagator))->quotient = -12;
    ((PredicateDivConstant*)(cc21nd.con2.propagator))->quotient = -12;
    ((PredicateDivConstant*)(cc21nd.con3.propagator))->quotient = -12;


    ConChecker< PredicateDiv > cc22("%",true,true,12345,3,50,0,true);
    cc22.init();
    cc22.run(100);


   /*

    ConChecker< PredicateDiv > cc16a("/", true,true,12345,3,20,0);
    cc16a.init();
    cc16a.run(100);

   */

    std::cout << " ";

  }

};


class URCSPGenerator{
private:
  int comb(int b, int d);  
  int i2comb(int& code, int ind, int base, int dim, int stop);
  
  //Random_ generator;
  
  int Var;
  int Dom;
  int Con; 
  int Ari; 
  int Ngd;
  
  
  int c;
  int t;
  int selectedNG;
  int selectedCT; 
  int PossibleCTs; 
  int PossibleNGs;
  
  int *CTarray;
  int *NGarray;
  
  int uv;
  int *unconnected_var; 
  
public: 
  long Seed;
  
  /**@name Constructors*/
  //@{
  /// 
  URCSPGenerator(int S, int v, int d, int c, int a, int n);
  void initURCSP( int v, int d, int c, int a, int n );
  ///
  ~URCSPGenerator();
  //@}
  
  /// An array of indices correponding to the constrained variables
  int *vars;
  /// An array corresponding to a nogood, such that vals[vars[i]] is the ith value of the nogood
  int *vals;
  /**
   *  Initialisation function, to be called after every CSP generation
   */
  void reInit();
  void reInit(int v, int d, int c, int a, int n);
  /**
   *  Generate a constraint scope and store it in vars
   */
  bool erateConstraint();
  /**
   *  Generate a nogood tuple and store it in vals
   */
  bool erateNogood();
  
};


//#define


int main(int argc, char *argv[])
{  

  // // BitSet s(-100, 200, BitSet::empt);

  // // s.add(-43);
  // // s.add(-9);
  // // s.add(-7);
  // // s.add(12);
  // // s.add(13);
  // // s.add(14);
  // // s.add(16);
  // // s.add(32);
  // // s.add(63);
  // // s.add(127);
  // // s.add(128);
  // // s.add(187);

  
  // // std::cout << s << std::endl;

  // // int *b = new int[12];
  // // b[0] = -43;
  // // s.iterate_into(12, b);

  // // for(int i=0; i<12; ++i) {
  // //   std::cout << b[i] << std::endl;
  // // }




  // Vector<int> vals;

  // vals.add(-43);
  // vals.add(-9);
  // vals.add(-7);
  // vals.add(12);
  // vals.add(13);
  // vals.add(14);
  // vals.add(16);
  // vals.add(32);
  // vals.add(63);
  // vals.add(127);
  // vals.add(128);
  // vals.add(187);

  // Variable X(vals, BITSET_VAR);

  // Variable Y(vals, LIST_VAR);

  // Variable Z(-13, 7, RANGE_VAR);

  // Solver s;
  // s.add(X);
  // s.add(Y);
  // s.add(Z);

  // Domain dom_x(X);
  // Domain dom_y(Y);
  // Domain dom_z(Z);


  // Domain::iterator xit = dom_x.begin();
  // Domain::iterator xend = dom_x.end();

  // Domain::iterator yit = dom_y.begin();
  // Domain::iterator yend = dom_y.end();

  // Domain::iterator zit = dom_z.begin();
  // Domain::iterator zend = dom_z.end();

  // cout << "bitset domain: " << X.get_domain() << "\n";
  // while(xit != xend) {
  //   cout << dom_x.get_value(xit) << endl;
  //   ++xit;
  // }

  // cout << "\nlist domain: " << Y.get_domain() << "\n";
  // while(yit != yend) {
  //   cout << dom_y.get_value(yit) << endl;
  //   ++yit;
  // }

  // cout << "\nrange domain: " << Z.get_domain() << "\n";
  // int k = 0;
  // while(zit != zend && k < 30) {
  //   cout << dom_z.get_value(zit) << endl;
  //   ++zit;
  //   ++k;
  // }


  // exit(1);



  usrand(12345);

  std::vector<UnitTest*> tests;

  int N = 8; //atoi(argv[1]);
  if(argc>1) N=atoi(argv[1]);

  /*
  tests.push_back(new CheckerTest());
  tests.push_back(new SymmetricDifferenceTest());
  tests.push_back(new LexTest());
  tests.push_back(new UnionTest());
  tests.push_back(new IntersectionTest());
  tests.push_back(new SubsetTest());
  tests.push_back(new EqualSetTest());
  tests.push_back(new SetDifferenceTest());
  tests.push_back(new CardTest());
  tests.push_back(new DivTest());
  tests.push_back(new MinMaxTest());
  tests.push_back(new WeightedSumTest());
  tests.push_back(new ElementTest());
  tests.push_back(new MemberTest());
  tests.push_back(new RewriteTest1());
  tests.push_back(new OpshopTest());
  tests.push_back(new BoolPigeons(N+1, EXPRESSION));
  tests.push_back(new BoolPigeons(N+1, BITSET_VAR));
  */
  tests.push_back(new SatTest());
  /*
  tests.push_back(new Pigeons(N+2)); 
  tests.push_back(new CostasAllDiffAllSolutions(N+1, FORWARD_CHECKING));
  tests.push_back(new CostasAllDiffAllSolutions(N+1, BOUND_CONSISTENCY, RANGE_VAR));
  tests.push_back(new CostasAllDiffAllSolutions(N+1, BOUND_CONSISTENCY));
  tests.push_back(new CostasNotEqualAllSolutions(N+1));
  //tests.push_back(new RandomCListRandomRemoveAndRestore<4>());
  tests.push_back(new RandomDomainRandomRemoveRangeAndRestore());
  tests.push_back(new RandomDomainRandomSetDomainBitsetAndRestore());
  tests.push_back(new RandomDomainRandomSetDomainAndRestore());
  tests.push_back(new RandomDomainRandomSetMaxAndRestore());
  tests.push_back(new RandomDomainRandomSetMinAndRestore());
  tests.push_back(new RandomDomainRandomRemove());
  // tests.push_back(new RandomRevNumAffectations<int>());
  // //tests.push_back(new ConstraintArrayTest());
  // tests.push_back(new RandomIntervalTest());
  */
 
  //tests[0]->Verbosity = HIGH;
  //tests[0]->Quality = HIGH;
  //tests[0]->Quantity = EXTREME;
  //tests[0]->Verbosity = EXTREME;


  while(tests.size() > 0) {
    UnitTest *t = tests.back();

    //t->Verbosity = EXTREME;
    
    double TIME = get_run_time();
    t->run();
    cout << (get_run_time() - TIME) << endl ;
    
    delete t;
    tests.pop_back();
  }

}


UnitTest::UnitTest() {
  Verbosity=LOW;
  Quality=MEDIUM;
  Quantity=MEDIUM;
}
  
UnitTest::UnitTest(const int vb, const int ql, const int qt) 
{
  Verbosity = vb;
  Quality = ql;
  Quantity = qt;
}

UnitTest::~UnitTest() {
}

void UnitTest::checkDomainIntegrity(Variable X) {
  // get all values, min, max and size through iteration
  int nxt = X.get_min();
  int v=nxt;
  
  int xmin = INFTY;
  int xmax = -INFTY;
  Vector<int> values;
  int max_iteration = 10000;

  do {
    if(--max_iteration < 0) {
      cout << "Error while iterating (infinite loop?)!" << endl
	   << X << " in " << X.get_domain() << " / " << values << endl;
      exit(1);
    }
    
    for(int k=v+1; k<nxt; ++k) {
      if(X.contain(k)) {
	cout << "Error while iterating (some values were missed)!" << endl
	     << X << " in " << X.get_domain() << " / " << values << endl;
	exit(1);
      }
    }

    v=nxt;
    if(xmin > v) xmin = v;
    if(xmax < v) xmax = v;
    values.add(v);

    if(!X.contain(v)) {
      cout << "Error while iterating (got non-member values)!" << endl
	   << X << " in " << X.get_domain() << " / " << values << endl;
      exit(1);
    }

    nxt = X.next(v);    
  } while(nxt != v);

  if(xmin != X.get_min()) {
    cout << "lower bound not properly maintained!" << endl
	 << X << " in " << X.get_domain() << " / " << values 
	 << ": " << xmin << " vs " << X.get_min() << endl;
    exit(1);
  }

  if(xmax != X.get_max()) {
    cout << "upper bound not properly maintained!" << endl
	 << X << " in " << X.get_domain() << " / " << values 
	 << ": " << xmax << " vs " << X.get_max() << endl;
    exit(1);
  }

  if(values.size != X.get_size()) {
    cout << "domain size not properly maintained!" << endl
	 << X << " in " << X.get_domain() << " / " << values 
	 << ": " << values.size << " vs " << X.get_size() << endl;
    exit(1);
  }
}



// int _modulo_fct_(const int x, const int m) {
//   int mod = x%m;
//   if(mod && (mod<0) != (m<0))  mod += m;
//   return mod;
// }


RandomIntervalTest::RandomIntervalTest(const int ql, 
				       const int qt) 
  : UnitTest(LOW, ql, qt) {}
RandomIntervalTest::~RandomIntervalTest() {}

void RandomIntervalTest::run() {

  int n=100000;
  int m=51;

  if(Verbosity) cout << "Run " << n << " random interval checks ";// << endl;

  int lb, ub, val;
  bool all_in, min_reached, max_reached;
  Interval I, J, K;    
  double dval;
  for(int count=0; count<n; ++count) {

    if(!(count % (n/10))) {
      std::cout << ".";
      std::cout.flush();
    }

    lb = randint(m)-m/2;
    ub = lb+randint(m);
    I = Interval(lb, ub);

    lb = randint(m)-m/2;
    ub = lb+randint(m);
    J = Interval(lb, ub);

    K = (I*J);

    // std::cout << "[" << I.min << "," << I.max 
    // 	      << "] * [" << J.min << "," << J.max << "] = ["
    // 	      << K.min << "," << K.max << "]\n";

    all_in = true;
    min_reached = false;
    max_reached = false;
    for(int i=I.min; i<=I.max; ++i)
      for(int j=J.min; j<=J.max; ++j)
	{
	  val = i*j;
	  if(val == K.max) max_reached = true;
	  if(val == K.min) min_reached = true;
	  if(val > K.max || val < K.min) all_in =false;
	  
	  if(!all_in) {
	    std::cout << "Error: " << i << "*" << j << " is not in [" << K.min << "," << K.max 
		      << "] != ([" << I.min << "," << I.max 
		      << "] * [" << J.min << "," << J.max << "])\n";
	    exit(1);
	  }
	}
    if(!min_reached) {
      std::cout << "Error: " << K.min << " is not min([" << I.min << "," << I.max 
		<< "] * [" << J.min << "," << J.max << "])\n";
      exit(1);
    }
    if(!max_reached) {
      std::cout << "Error: " << K.max << " is not max([" << I.min << "," << I.max 
		<< "] * [" << J.min << "," << J.max << "])\n";
      exit(1);
    }


    K = (I/J);

    // std::cout << "[" << I.min << "," << I.max 
    // 	      << "] / [" << J.min << "," << J.max << "] = ["
    // 	      << K.min << "," << K.max << "]\n";

    all_in = true;
    min_reached = false;
    max_reached = false;
    for(int i=I.min; i<=I.max; ++i)
      for(int j=J.min; j<=J.max; ++j)
	if(j) {
	  val = i/j;
	  if(val == K.max) max_reached = true;
	  if(val == K.min) min_reached = true;
	  if(val > K.max || val < K.min) all_in =false;
	  
	  if(!all_in) {
	    std::cout << "Error: " << i << "/" << j << " is not in [" << K.min << "," << K.max 
		      << "] != ([" << I.min << "," << I.max 
		      << "] / [" << J.min << "," << J.max << "])\n";
	    exit(1);
	  }
	}

    if(!K.empty()) {
      if(!min_reached) {
	std::cout << "Error: " << K.min << " is not min([" << I.min << "," << I.max 
		  << "] / [" << J.min << "," << J.max << "])\n";
	exit(1);
      }
      if(!max_reached) {
	std::cout << "Error: " << K.max << " is not max([" << I.min << "," << I.max 
		  << "] / [" << J.min << "," << J.max << "])\n";
	exit(1);
      }
    }



    K = (I.anti_mul(J)); // K*J = I

    // std::cout << "[" << I.min << "," << I.max 
    //  	      << "] // [" << J.min << "," << J.max << "] = ["
    //  	      << K.min << "," << K.max << "]\n";

    if(!K.empty()) { //[TODO: some knd of check when it is empty]
      all_in = true;
      min_reached = false;
      max_reached = false;
      for(int i=I.min; i<=I.max; ++i)
	for(int j=J.min; j<=J.max; ++j) 
	  if(j) {
	    //if(i/j || (J.min <= 0 && 0 <= J.max)) {
	      dval = (double)i/(double)j;
	      //if(val == I.max) max_reached = true;
	      //if(val == I.min) min_reached = true;
	      
	      //std::cout << i << "//" << j << " <-> " << (double)(K.min-1) << " <=? " << dval << " <=? " << (double)(K.max+1) << std::endl;
	      
	      if(dval >= (double)(K.max+1) || dval <= (double)(K.min-1)) all_in = false;
	      
	      if(!all_in) {
		std::cout << "Error: " << i << "//" << j << " is not in [" << K.min << "," << K.max 
			  << "] != ([" << I.min << "," << I.max 
			  << "] // [" << J.min << "," << J.max << "])\n";
		exit(1);
	      }
	      //}
	  }
      
      if(K.min > -INFTY) {
	// we check that K.min-1 is not supported
	for(int j=J.min; j<=J.max; ++j) {
	  val = ((K.min-1) * j);
	  if(I.min <= val && val <= I.max) min_reached = true;
	}
	if(min_reached) {
	  std::cout << "Error: " << K.min-1 << " is a valid min([" << I.min << "," << I.max 
		    << "] // [" << J.min << "," << J.max << "])\n";
	  exit(1);
	}
      }
      
      if(K.max < INFTY) {
	// we check that K.max+1 is not supported
	for(int j=J.min; j<=J.max; ++j) {
	  val = ((K.max+1) * j);
	  if(I.min <= val && val <= I.max) max_reached = true;
	}
	if(max_reached) {
	  std::cout << "Error: " << K.max+1 << " is a valid max([" << I.min << "," << I.max 
		    << "] // [" << J.min << "," << J.max << "])\n";
	  exit(1);
	}
      }
    }



    for(int modulo=1; modulo<m/2; ++modulo) {
      K = (I.operator_modulo(modulo));

      // std::cout << "[" << I.min << "," << I.max 
      // 		<< "] % " << modulo << " = ["
      // 		<< K.min << "," << K.max << "]\n";
      
      all_in = true;
      min_reached = false;
      max_reached = false;
      for(int i=I.min; i<=I.max; ++i) {
	val = __modulo_fct__(i,modulo);

	//std::cout << "[" << K.min << "][" << K.max << "] " << i << "%" << modulo << "=" << val << std::endl;

	if(val == K.max) max_reached = true;
	if(val == K.min) min_reached = true;
	if(val > K.max || val < K.min) all_in =false;
	  
	if(!all_in) {
	  std::cout << "Error: " << i << "%" << modulo << " is not in [" << K.min << "," << K.max 
		    << "] != ([" << I.min << "," << I.max 
		    << "] % " << modulo << ")\n";
	  exit(1);
	}
      }

      if(!min_reached) {
	std::cout << "Error: " << K.min << " is not min([" << I.min << "," << I.max 
		  << "] % " << modulo << ")\n";
	exit(1);
      }
      if(!max_reached) {
	std::cout << "Error: " << K.max << " is not max([" << I.min << "," << I.max 
		  << "] % " << modulo << ")\n";
	exit(1);
      }
    }

    for(int modulo=1; modulo<m/2; ++modulo) {
      K = (I%modulo);

      // std::cout << "[" << I.min << "," << I.max 
      // 		<< "] % " << modulo << " = ["
      // 		<< K.min << "," << K.max << "]\n";
      
      all_in = true;
      min_reached = false;
      max_reached = false;
      for(int i=I.min; i<=I.max; ++i) {
	val = i%modulo;

	//std::cout << "[" << K.min << "][" << K.max << "] " << i << "%" << modulo << "=" << val << std::endl;

	if(val == K.max) max_reached = true;
	if(val == K.min) min_reached = true;
	if(val > K.max || val < K.min) all_in =false;
	  
	if(!all_in) {
	  std::cout << "Error: " << i << "%" << modulo << " is not in [" << K.min << "," << K.max 
		    << "] != ([" << I.min << "," << I.max 
		    << "] % " << modulo << ")\n";
	  exit(1);
	}
      }

      if(!min_reached) {
	std::cout << "Error: " << K.min << " is not min([" << I.min << "," << I.max 
		  << "] % " << modulo << ")\n";
	exit(1);
      }
      if(!max_reached) {
	std::cout << "Error: " << K.max << " is not max([" << I.min << "," << I.max 
		  << "] % " << modulo << ")\n";
	exit(1);
      }
    }
  }

  //std::cout << std::endl;
}


RandomDomainRandomRemove::RandomDomainRandomRemove(const int ql, 
						   const int qt) 
  : UnitTest(LOW, ql, qt) {}
RandomDomainRandomRemove::~RandomDomainRandomRemove() {}

void RandomDomainRandomRemove::run() {

  //Verbosity = EXTREME;

  /// create vars with random domains and remove a random set of values.
  Vector<int> values;
  Vector<int> rand_values;

  int N = (1<<(1<<Quantity));

  int max_size = (1<<(1<<Quality)); // 2, 4, 16, 256, 65536
  int max_min = (1<<(1<<Quality));
  
  int step = (Quantity ? (1<<Quality)/(Quantity) : INFTY);

  if(Verbosity) cout << "Run random removal checks [2.." 
		     << max_size << "] x [2.." << max_min 
		     << "] (" << step << ") ";// << endl;
  if(Verbosity>LOW) cout << endl;

  for(int dom_size = 2; dom_size <= max_size ; dom_size += step) {
    if(Verbosity > LOW) {
      cout << " do " << (N*(max_min-2)/step) 
	   << " checks with domain size = " << dom_size << endl;
    }
      
    for(int dom_min = 2; dom_min <= max_min ; dom_min += step) {
      
      if(Verbosity > LOW) {
	cout << "  do " << N << " checks with domains in [-" << (dom_min/2) 
	     << ".." << (dom_size-(dom_min/2)) << "] to [" << (dom_min/2) 
	     << ".." << (dom_size+(dom_min/2)) << "]" << endl;
      }

      for(int i=0; i<N; ++i) {

	int lb = randint(dom_min) - (dom_min/2);
	int ub = randint(dom_size) + lb + 1;
	
	int nvalues = dom_size;
	for(int j=0; j<nvalues; ++j) 
	  rand_values.add(randint(2*dom_size)-(ub-lb+1));
	
	Solver s;
	Variable X(lb,ub);

	//std::cout << "X.domain_type: " << X.domain_type << " / " << X.variable << std::endl;

	//X.add_to(&s);
	//s.add(X);
	//X = X.get_var();
	
	//std::cout << "X.domain_type: " << X.domain_type << " / " << X.variable << std::endl;

	//s.initialise();
	
	s.add(X);
	X = X.get_var();

	//cout << s << endl;

	if(Verbosity > MEDIUM) {
	  cout << "    " << X << " in " << X.get_domain() << ": " << endl;
	  cout << "    remove " << rand_values << endl; 
	}

	for(int j=0; j<nvalues; ++j) { 
	  
	  if(Verbosity > HIGH) {
	    cout << "      remove " << rand_values[j] << " from " 
		 << X << " in " << X.get_domain() << " => "; 
	    cout.flush();
	  }
	  Event evt = X.remove(rand_values[j]);
	  X = X.get_var();
	  
	  if( evt != FAIL_EVENT) {
	    if(Verbosity > HIGH) {
	      cout << "    " << X << " in " << X.get_domain() << endl; 
	    }
	  
	    checkDomainIntegrity(X);
	  
	    if(X.contain(rand_values[j])) {
	      cout << "Error on remove " << rand_values[j]
		   << " from " << X << " in " << X.get_domain() << endl;
	      exit(1);
	    }
	  } else {
	    if(!X.equal(rand_values[j]))  {
	      cout << "Error on remove " << rand_values[j]
		   << " from " << X << " in " << X.get_domain() << endl;
	      exit(1);
	    }
	      
	    if(Verbosity > HIGH) {
	      cout << "    Wipe out!" << endl;
	    }
	  }
	}
	
	if(Verbosity > MEDIUM) {
	  cout << " OK (" << X << " in " << X.get_domain() << ")" << endl;
	}

	rand_values.clear();
      }      
    }
  }
}


  
RandomDomainRandomSetDomainAndRestore::RandomDomainRandomSetDomainAndRestore(const int ql, 
									     const int qt) 
  : UnitTest(LOW, ql, qt) {}
RandomDomainRandomSetDomainAndRestore::~RandomDomainRandomSetDomainAndRestore() {}

void RandomDomainRandomSetDomainAndRestore::run() {
  int N = (1<<(1<<Quantity));

  if(Verbosity) cout << "Run " << N << " random set_domain checks ";// << endl;
  if(Verbosity>LOW) cout << endl;

  for(int i=0; i<N; ++i) {
    int dom_size = randint(1<<(1<<Quality))+1;

    if(Verbosity > LOW) cout << "  Domain size: " << dom_size << endl;

    for(int j=0; j<dom_size; ++j) {

      Solver s;
      int dom_min = randint(dom_size) - (dom_size/2);
      int lb = dom_min;
      int ub = dom_min+dom_size;
      Variable X(lb,ub);
      s.add(X);
      //s.initialise();

      if(Verbosity > MEDIUM) cout << "    " << X << " in " << X.get_domain() << endl;

      for(int k=0; k<50; ++k) {

	if(Verbosity > HIGH) cout << "      <- " << lb-k-1 ;

	if(X.set_domain(lb-k-1) != FAIL_EVENT) {
	  cout << "Error on [set_domain] " << lb-k-1
	       << " on " << X << " in " << X.get_domain() << endl;
	  exit(1);
	}
	if(X.get_min() != lb || X.get_max() != ub || X.get_size() != (unsigned int)(ub-lb+1)) {
	  cout << "Error on [FAIL] " << X << " in " << X.get_domain() 
	       << " / [" << lb << ".." << ub << "] sz=" << (ub-lb+1) << endl;
	  exit(1);
	}   

	if(Verbosity > HIGH) cout << " (" << X << " in " 
				  << X.get_domain() << ")" << endl;  

	checkDomainIntegrity(X);
      }

      for(int k=0; k<50; ++k) {

	if(Verbosity > HIGH) cout << "      <- " << lb+k+1 ;

	if(X.set_domain(ub+k+1) != FAIL_EVENT) {
	  cout << "Error on [set_domain] " << ub+k+1
	       << " on " << X << " in " << X.get_domain() << endl;
	  exit(1);
	}
	if(X.get_min() != lb || X.get_max() != ub || X.get_size() != (unsigned int)(ub-lb+1)) {
	  cout << "Error on [FAIL] " << X << " in " << X.get_domain() 
	       << " / [" << lb << ".." << ub << "] sz=" << (ub-lb+1) << endl;
	  exit(1);
	}

	if(Verbosity > HIGH) cout << " (" << X << " in " 
				  << X.get_domain() << ")" << endl;       

	checkDomainIntegrity(X);
      }

      for(int k=0; k<dom_size; ++k) {

	if(Verbosity > HIGH) cout << "      <- " << lb+k ;

	X.set_domain(lb+k);
	if(!X.equal(lb+k)) {
	  cout << "Error on [set_domain] " << lb+k
	       << " on " << X << " in " << X.get_domain() << endl;
	  exit(1);
	}

	if(Verbosity > HIGH) cout << " (" << X << " in " 
				  << X.get_domain() << ")" << endl;
	  
	checkDomainIntegrity(X);

	X.restore();
	  
	if(X.get_min() != lb || X.get_max() != ub || X.get_size() != (unsigned int)(ub-lb+1)) {
	  cout << "Error on [restore] " << X << " in " << X.get_domain() << " / [" << lb << ".." << ub << "]" << endl;
	  exit(1);
	}
	  
	checkDomainIntegrity(X);
      }
    }
  }
}


  
RandomDomainRandomSetMinAndRestore::RandomDomainRandomSetMinAndRestore(const int ql, 
								       const int qt) 
  : UnitTest(LOW, ql, qt) {}
RandomDomainRandomSetMinAndRestore::~RandomDomainRandomSetMinAndRestore() {}

void RandomDomainRandomSetMinAndRestore::run() {
  int N = (1<<(1<<Quantity));

  //Verbosity = EXTREME;

  if(Verbosity) cout << "Run " << N << " random set_min checks " ;// << endl;
  if(Verbosity>LOW) cout << endl;

  for(int i=0; i<N; ++i) {
    int dom_size = randint(1<<(1<<Quality))+1;

    if(Verbosity > LOW) cout << "  Domain size: " << dom_size << endl;

    for(int j=0; j<dom_size; ++j) {

      Solver s;
      int dom_min = randint(dom_size) - (dom_size/2);
      int lb = dom_min;
      int ub = dom_min+dom_size;
      Variable X(lb,ub);
      s.add(X);
      //s.initialise();

      if(Verbosity > MEDIUM) cout << "    " << X << " in " << X.get_domain() << endl;

      for(int k=0; k<50; ++k) {

	if(Verbosity > HIGH) cout << "      >= " << lb-k ;

	if(X.set_min(lb-k) != NO_EVENT) {
	  cout << "Error on [set_min] " << lb-k
	       << " on " << X << " in " << X.get_domain() << endl;
	  exit(1);
	}
	if(X.get_min() != lb || X.get_max() != ub || X.get_size() != (unsigned int)(ub-lb+1)) {
	  cout << "Error on [FAIL] " << X << " in " << X.get_domain() 
	       << " / [" << lb << ".." << ub << "] sz=" << (ub-lb+1) << endl;
	  exit(1);
	}
	  
	if(Verbosity > HIGH) cout << " (" << X << " in " 
				  << X.get_domain() << ")" << endl;  
	  
	checkDomainIntegrity(X);
      }
	
      for(int k=0; k<50; ++k) {
	  
	if(Verbosity > HIGH) cout << "      >= " << ub+k+1 ;
	  
	if(X.set_min(ub+k+1) != FAIL_EVENT) {
	  cout << "Error on [set_min] " << ub+k+1
	       << " on " << X << " in " << X.get_domain() << endl;
	  exit(1);
	}
	if(X.get_min() != lb || X.get_max() != ub || X.get_size() != (unsigned int)(ub-lb+1)) {
	  cout << "Error on [FAIL] " << X << " in " << X.get_domain() 
	       << " / [" << lb << ".." << ub << "] sz=" << (ub-lb+1) << endl;
	  exit(1);
	}

	if(Verbosity > HIGH) cout << " (" << X << " in " 
				  << X.get_domain() << ")" << endl;       

	checkDomainIntegrity(X);
      }

      for(int k=0; k<dom_size; ++k) {

	if(Verbosity > HIGH) cout << "      >= " << lb+k+1 ;

	X.set_min(lb+k+1);
	if(X.get_min() < lb+k+1) {
	  cout << "Error on [set_min] " << lb+k+1
	       << " on " << X << " in " << X.get_domain() << endl;
	  exit(1);
	}

	if(Verbosity > HIGH) cout << " (" << X << " in " 
				  << X.get_domain() << ")" << endl;
	  
	checkDomainIntegrity(X);

	X.restore();
	  
	if(X.get_min() != lb || X.get_max() != ub || X.get_size() != (unsigned int)(ub-lb+1)) {
	  cout << "Error on [restore] " << X << " in " << X.get_domain() << " / [" << lb << ".." << ub << "]" << endl;
	  exit(1);
	}
	  
	checkDomainIntegrity(X);
      }
    }
  }
}

  
RandomDomainRandomSetMaxAndRestore::RandomDomainRandomSetMaxAndRestore(const int ql, 
								       const int qt) 
  : UnitTest(LOW, ql, qt) {}
RandomDomainRandomSetMaxAndRestore::~RandomDomainRandomSetMaxAndRestore() {}

void RandomDomainRandomSetMaxAndRestore::run() {
  int N = (1<<(1<<Quantity));

  //Verbosity = EXTREME;

  if(Verbosity) cout << "Run " << N << " random set_max checks " ; // << endl;
  if(Verbosity>LOW) cout << endl;

  for(int i=0; i<N; ++i) {
    int dom_size = randint(1<<(1<<Quality))+1;

    if(Verbosity > LOW) cout << "  Domain size: " << dom_size << endl;

    for(int j=0; j<dom_size; ++j) {

      Solver s;
      int dom_min = randint(dom_size) - (dom_size/2);
      int lb = dom_min;
      int ub = dom_min+dom_size;
      Variable X(lb,ub);
      s.add(X);
      //s.initialise();

      if(Verbosity > MEDIUM) cout << "    " << X << " in " << X.get_domain() << endl;

      for(int k=0; k<50; ++k) {

	if(Verbosity > HIGH) cout << "      <= " << lb-k-1 ;

	if(X.set_max(lb-k-1) != FAIL_EVENT) {
	  cout << "Error on [set_max] " << lb-k-1
	       << " on " << X << " in " << X.get_domain() << endl;
	  exit(1);
	}
	if(X.get_min() != lb || X.get_max() != ub || X.get_size() != (unsigned int)(ub-lb+1)) {
	  cout << "Error on [FAIL] " << X << " in " << X.get_domain() 
	       << " / [" << lb << ".." << ub << "] sz=" << (ub-lb+1) << endl;
	  exit(1);
	}   

	if(Verbosity > HIGH) cout << " (" << X << " in " 
				  << X.get_domain() << ")" << endl;

	checkDomainIntegrity(X);  
      }

      for(int k=0; k<50; ++k) {

	if(Verbosity > HIGH) cout << "      <= " << ub+k ;

	if(X.set_max(ub+k) != NO_EVENT) {
	  cout << "Error on [set_max] " << ub+k
	       << " on " << X << " in " << X.get_domain() << endl;
	  exit(1);
	}
	if(X.get_min() != lb || X.get_max() != ub || X.get_size() != (unsigned int)(ub-lb+1)) {
	  cout << "Error on [NO EVENT] " << X << " in " << X.get_domain() 
	       << " / [" << lb << ".." << ub << "] sz=" << (ub-lb+1) << endl;
	  exit(1);
	}

	if(Verbosity > HIGH) cout << " (" << X << " in " 
				  << X.get_domain() << ")" << endl;

	checkDomainIntegrity(X);       
      }

      for(int k=0; k<dom_size; ++k) {

	if(Verbosity > HIGH) cout << "      <= " << lb+k ;

	X.set_max(lb+k);
	if(X.get_max() != lb+k) {
	  cout << "Error on [set_max] " << lb+k
	       << " on " << X << " in " << X.get_domain() << endl;
	  exit(1);
	}

	if(Verbosity > HIGH) cout << " (" << X << " in " 
				  << X.get_domain() << ")" << endl;
	  
	checkDomainIntegrity(X);

	X.restore();
	  
	if(X.get_min() != lb || X.get_max() != ub || X.get_size() != (unsigned int)(ub-lb+1)) {
	  cout << "Error on [restore] " << X << " in " << X.get_domain() << " / [" << lb << ".." << ub << "]" << endl;
	  exit(1);
	}
	  
	checkDomainIntegrity(X);
      }
    }
  }
}

  
RandomDomainRandomSetDomainBitsetAndRestore::RandomDomainRandomSetDomainBitsetAndRestore(const int ql, 
											 const int qt) 
  : UnitTest(LOW, ql, qt) {}
RandomDomainRandomSetDomainBitsetAndRestore::~RandomDomainRandomSetDomainBitsetAndRestore() {}

void RandomDomainRandomSetDomainBitsetAndRestore::run() {

  //Verbosity = EXTREME;
  

  int N = (1<<(1<<Quantity));

  if(Verbosity) cout << "Run " << N << " random set_domain (bitset) checks " ; // << endl;
  if(Verbosity>LOW) cout << endl;

  for(int i=0; i<N; ++i) {
    int dom_size = randint(1<<(1<<Quality))+1;

    if(Verbosity > LOW) cout << "  Domain size: " << dom_size << endl;

    for(int j=0; j<dom_size; ++j) {

      int dom_min = randint(dom_size) - (dom_size/2);
      int lb = dom_min;
      int ub = dom_min+dom_size;

      Solver s;
      Variable X(lb,ub);
      s.add(X);
      //s.initialise();

      if(Verbosity > LOW) cout << "    " << X << " in " << X.get_domain() << endl;

      BitSet vals(lb-50, ub+50, BitSet::empt);

	
      for(int k=0; k<dom_size; ++k) {

	int set_size = randint(dom_size+10);
	  
	for(int l=0; l<set_size; ++l) {
	  vals.add(randint(ub-lb+101) + (lb-50));
	}
	  
	if(Verbosity > MEDIUM) cout << "    Set: " << vals << ": ";

	Event evt = X.set_domain(vals);
	X = X.get_var();

	if(Verbosity > MEDIUM) cout << X << " in " << X.get_domain() << endl;

	if(evt != FAIL_EVENT) {
	  int z = 0;
	  int nxt = X.get_min();
	  do {
	    z = nxt;
	    if(!vals.contain(z)) {
	      cout << "Error on [set_domain] " 
		   << " on " << X << " in " << X.get_domain() << endl;
	      exit(1);
	    }
	    nxt = X.next(z);
	  } while( nxt != z );
	    
	  
	  vals.flip();
	  z = 0;
	  nxt = vals.min();
	  do {
	    z = nxt;
	    if(X.contain(z)) {
	      cout << "Error on [set_domain] " 
		   << " on " << X << " in " << X.get_domain() << endl;
	      exit(1);
	    }
	    nxt = vals.next(z);
	  } while( nxt != z );
	    
	  checkDomainIntegrity(X);
	  
	  if(evt != NO_EVENT) X.restore();

	  checkDomainIntegrity(X);

	}
	  
	if(X.get_min() != lb || X.get_max() != ub || X.get_size() != (unsigned int)(ub-lb+1)) {
	  cout << X.get_min() << " " << lb << endl;
	  cout << X.get_max() << " " << ub << endl;
	  cout << X.get_size() << " " << (ub-lb+1) << endl;
	  cout << "Error on [restore] " << X << " in " << X.get_domain() << " / [" << lb << ".." << ub << "]" << endl;
	  exit(1);
	}
	  
	checkDomainIntegrity(X);
      }
    }
  }
}

  
RandomDomainRandomRemoveRangeAndRestore::
RandomDomainRandomRemoveRangeAndRestore(const int ql, 
					const int qt) 
  : UnitTest(LOW, ql, qt) {}
RandomDomainRandomRemoveRangeAndRestore::
~RandomDomainRandomRemoveRangeAndRestore() {}

void RandomDomainRandomRemoveRangeAndRestore::run() {
  int N = (1<<(1<<Quantity));

  if(Verbosity) cout << "Run " << N << " random set_domain + remove_interval checks " ;// << endl;
  if(Verbosity>LOW) cout << endl;

  //Verbosity = EXTREME;


  for(int i=0; i<N; ++i) {
    int dom_size = randint(1<<(1<<Quality))+1;

    if(Verbosity > LOW) cout << "  Domain size: " << dom_size << endl;

    for(int j=0; j<dom_size; ++j) {

      int dom_min = randint(dom_size) - (dom_size/2);
      int lb = dom_min;
      int ub = dom_min+dom_size;

      Solver s;
      Variable X(lb,ub);
      //s.add(X);
      X.add_to(&s);
      //s.initialise();

      if(Verbosity > MEDIUM) cout << "    " << X << " in " << X.get_domain() << endl;

      BitSet vals(lb-50, ub+50, BitSet::empt);
	
      int set_size = randint(dom_size+10);
	  
      for(int l=0; l<set_size; ++l) {
	vals.add(randint(ub-lb+101) + (lb-50));
      }
	  
      if(Verbosity > MEDIUM) cout << "    Set: " << vals << ": " ;
      s.save();
	
      Event evt1 = X.set_domain(vals);

      X = X.get_var();

      if(Verbosity > MEDIUM) cout << X << " in " << X.get_domain() << endl;

      if(evt1 != FAIL_EVENT) {	
	int z = 0;
	int nxt = X.get_min();
	do {
	  z = nxt;
	  if(!vals.contain(z)) {
	    cout << "Error on [set_domain] " 
		 << " on " << X << " in " << X.get_domain() << endl;
	    exit(1);
	  }
	  nxt = X.next(z);
	} while( nxt != z );
	  
	  
	vals.flip();
	z = 0;
	nxt = vals.min();
	do {
	  z = nxt;
	  if(X.contain(z)) {
	    cout << "Error on [set_domain] " 
		 << " on " << X << " in " << X.get_domain() << endl;
	    exit(1);
	  }
	  nxt = vals.next(z);
	} while( nxt != z );
	  
	checkDomainIntegrity(X);
	  
	for(int l=lb; l<ub; ++l) {
	  for(int m=l+1; m<=ub; ++m) {
	    s.save();
	    if(Verbosity > HIGH) cout << "        remove [" << l << ".." 
				      << m << "] from " << X << " in " 
				      << X.get_domain() << ": ";
	    Event evt2 = X.remove_interval(l,m);

	    X = X.get_var();

	    if(Verbosity > HIGH) cout << X << " in " << X.get_domain() << endl;
	    if(evt2 != NO_EVENT && evt2 != FAIL_EVENT) {
	      for(int v=l; v<=m; ++v) {
		if(X.contain(v)) {
		  cout << "Error on [remove_interval] " 
		       << " on " << X << " in " << X.get_domain() << endl;
		  exit(1);
		}
	      }
	      checkDomainIntegrity(X);
	    }

	    s.restore();
	  }
	}
      }
      
      s.restore();
      
      if(Verbosity > MEDIUM) {
	std::cout << "restore" << std::endl << s << std::endl;
      }

	
      if(X.get_min() != lb || X.get_max() != ub || X.get_size() != (unsigned int)(ub-lb+1)) {
	cout << X.get_min() << " " << lb << endl;
	cout << X.get_max() << " " << ub << endl;
	cout << X.get_size() << " " << (ub-lb+1) << endl;
	cout << "Error on [restore] " << X << " in " << X.get_domain() << " / [" << lb << ".." << ub << "]" << endl;
	exit(1);
      }
	
      checkDomainIntegrity(X);
    }
  }
}



// template< int NUM_HEAD >
// RandomCListRandomRemoveAndRestore< NUM_HEAD >::
// RandomCListRandomRemoveAndRestore(const int ql, 
// 				 const int qt) 
//   : UnitTest(LOW, ql, qt) { NUM_ELT = 100*(1 << (1 << Quality)); }

// template< int NUM_HEAD >
// RandomCListRandomRemoveAndRestore< NUM_HEAD >::
// ~RandomCListRandomRemoveAndRestore() {}

// template< int NUM_HEAD >
// void RandomCListRandomRemoveAndRestore< NUM_HEAD >::run() {

//   Solver s;
//   //s.initialise();

//   int N = (1<<(1<<Quantity));
//   if(Verbosity) cout << "Run " << N << " random ConstraintList checks " ;
//   if(Verbosity>LOW) cout << endl;
    
//   for(int iteration=0; iteration<N; ++iteration) {
      
//     ReversibleMultiList<int,NUM_HEAD> alist(&s);
//     //alist.env = &s;
      
//     NaiveImplementation blist(&s);
      
//     Vector<unsigned int> wl;
      
//     Vector<unsigned int> elts;
      
//     Vector<unsigned int> idx1;
//     Vector<unsigned int> idx2;
//     for(int i=0; i<NUM_ELT; ++i)
//       elts.add(i+1);
//     for(int i=0; i<NUM_ELT; ++i) {
//       int j=randint(NUM_ELT-i)+i;
//       elts[i] = elts[j];
//       elts[j] = i+1;
//     }
      
      
//     int i_elts = (1 << (1 << Quality));
//     int n_elts = i_elts;
//     int head;
//     for(int i=0; i<n_elts; ++i) {
//       head = randint(NUM_HEAD);
//       wl.add(head);

//       if( i >= (int)(elts.size) ) {
// 	cout << i << " " << NUM_ELT << endl;
// 	exit(1);
//       }

//       unsigned int the_elt = alist.create( elts[i], head );
//       idx1.add( the_elt );
//       idx2.add( i );
//       blist.create( elts[i], head );
//     }
      
//     checkEquality(alist, blist);
      
      
//     int branch_length = randint(2*i_elts);
      
//     BitSet isIn(0, NUM_ELT, BitSet::full);
//     isIn.set_max(i_elts-1);
      

//     if(Verbosity>LOW) cout << "Do a run on a branch of length " << branch_length << endl;

//     for(int i=0; i<branch_length; ++i) {
//       s.save();
	
//       if(!randint(2))
// 	for(int k=0; k<2; ++k) {
// 	  if(isIn.empty()) break;
	    
// 	  int j = randint(n_elts);
// 	  while(!isIn.contain(j))
// 	    j = randint(n_elts);
	    
// 	  alist.reversible_remove(idx1[j], wl[j]);
// 	  blist.reversible_remove(idx2[j], wl[j]);
	    
// 	  if(Verbosity>MEDIUM) {
// 	    for(int l=0; l<s.level; ++l) cout << " ";
// 	    cout << s.level << " remove " << (elts[idx2[j]]) << endl;
// 	    for(int l=0; l<s.level; ++l) cout << " ";
// 	    //alistcout << alist << endl;
// 	  }

// 	  isIn.remove(j);
	  
// 	}

//       if(!randint(3) && n_elts<NUM_ELT-1) {
// 	// add
// 	head = randint(NUM_HEAD);
// 	wl.add(head);
// 	idx1.add( alist.reversible_add( elts[n_elts], head ) );
// 	idx2.add( n_elts );
// 	blist.create( elts[n_elts], head );

// 	if(Verbosity>MEDIUM) {
// 	  for(int l=0; l<s.level; ++l) cout << " ";
// 	  cout << s.level << " add " << (elts[idx2[n_elts]]) << " to " 
// 	       << wl[n_elts] << "th list" << endl;
// 	  for(int l=0; l<s.level; ++l) cout << " ";
// 	  //cout << alist << endl;
// 	}

// 	isIn.add(n_elts);

// 	++n_elts;
//       }

//       checkEquality(alist, blist);
//     }
    
//     isIn.fill();
//     isIn.set_max(i_elts-1);

//     for(int i=0; i<branch_length; i++) {
//       s.restore();
//       //int j1 = (10+(50-i-1))%n_elts;
//       //int j2 = (10+(50-i-2))%n_elts;

//       //cout << "add " << (elts[idx2[j1]]) << " and " << (elts[idx2[j2]]) << endl;
//       if(Verbosity>MEDIUM) {
// 	for(int l=0; l<s.level; ++l) cout << " ";
// 	cout << s.level << " restore" << endl;
// 	for(int l=0; l<s.level; ++l) cout << " ";
// 	//cout << alist << endl;      
//       }
      
//       checkEquality(alist, blist);
//     }
//   }
// }

// template< int NUM_HEAD >
// void RandomCListRandomRemoveAndRestore< NUM_HEAD >::
// checkEquality(ReversibleMultiList<int,NUM_HEAD>& alist, 
// 	      RandomCListRandomRemoveAndRestore< NUM_HEAD >::
// 	      NaiveImplementation& blist) {
//   BitSet inA(0, NUM_ELT, BitSet::empt);
//   BitSet inB(0, NUM_ELT, BitSet::empt);

//   for(int i=0; i<NUM_HEAD; ++i) {
//     Node<int> nd = alist.first(i);
//     while(alist.next(nd)) {
//       inA.add((int)nd);
//     }
//     for(unsigned int j=0; j<blist.elements.size; ++j) {
//       if(blist.whichlist[j] >= i && 
// 	 blist.whichlevelremoved[j] > blist.env->level &&
// 	 blist.whichleveladded[j] <= blist.env->level) 
// 	inB.add(blist.elements[j]);
//     }

//     if(inA != inB) {
//       cout << "Discrepancy between the lists!" << endl
// 	   << inA << endl
// 	   << inB << endl;
//       exit(1);
//     }      
//     inA.clear();
//     inB.clear();
//   }
// }

// template< int NUM_HEAD >
// RandomCListRandomRemoveAndRestore< NUM_HEAD >::NaiveImplementation::NaiveImplementation(Solver *s) 
//   : Reversible(s) {}

// template< int NUM_HEAD >
// void RandomCListRandomRemoveAndRestore< NUM_HEAD >::NaiveImplementation::create(int elt, int wl) {
//   elements.add(elt);
//   whichlist.add(wl);
//   whichleveladded.add(env->level);
//   whichlevelremoved.add(INFTY);
// }

// template< int NUM_HEAD >
// void RandomCListRandomRemoveAndRestore< NUM_HEAD >::NaiveImplementation::reversible_remove(int elt, int wl) {
//   whichlevelremoved[elt] = env->level;
// }

// template< int NUM_HEAD >
// void RandomCListRandomRemoveAndRestore< NUM_HEAD >::NaiveImplementation::reversible_add(int elt, int wl) {
//   create(elt, wl);
// }

// template< int NUM_HEAD >
// void RandomCListRandomRemoveAndRestore< NUM_HEAD >::NaiveImplementation::save() {}
// template< int NUM_HEAD >
// void RandomCListRandomRemoveAndRestore< NUM_HEAD >::NaiveImplementation::restore() {}

// template< int NUM_HEAD >
// void RandomCListRandomRemoveAndRestore< NUM_HEAD >::NaiveImplementation::print(int num_lists) {
//   for(int i=0; i<num_lists; ++i) {
//     cout << "[";
//     for(unsigned int j=0; j<elements.size; ++j)
//       if(whichlist[j] >= i)
// 	cout << " " << elements[j];
//     cout << " ] ";
//   }
// }



template< class T >
RandomRevNumAffectations< T >::
RandomRevNumAffectations(const int ql, const int qt) 
  : UnitTest(LOW, ql, qt) { }
template< class T >
RandomRevNumAffectations< T >::
~RandomRevNumAffectations() {}

template< class T >
void RandomRevNumAffectations< T >::run() {
  if(Verbosity) cout << "Run " << 10 << " random Reversible<int> checks " ;

  Solver s;
  T x = 1000;
  int val = 0;

  ReversibleNum<T> rx(&s, x);

  int x1 = x;
  int x2;

  //s.initialise();

  for(int k=0; k<5; ++k) {
    s.save();
	
    for(int i=0; i<10; ++i) {
      rx += ((i+k)%20)-10;
      for(int j=0; j<5; ++j)
	val += rx;
    }

    x2 = rx;
      
    for(int i=0; i<999990; ++i) {
      if(!(i%10))s.save();
      rx += ((i+k)%20)-10;
      for(int j=0; j<5; ++j)
	val += rx;
    }
      
    for(int i=0; i<999990; i+=10) {
      s.restore();
    }
      
    if(x2 != rx) {
      cout << "Error: discrepancy between the values! " 
	   << x << " " << rx << endl;
      exit(1);
    }

    s.restore();

    if(x1 != rx) {
      cout << "Error: discrepancy between the values! " 
	   << x << " " << rx << endl;
      exit(1);
    }
  }
}
  
  
CostasAllDiffAllSolutions::CostasAllDiffAllSolutions(const int sz, const int ct, const int dt) 
  : UnitTest(LOW, LOW, LOW) { size=sz; consistency=ct; domain=dt; }
CostasAllDiffAllSolutions::~CostasAllDiffAllSolutions() {}

void CostasAllDiffAllSolutions::run() {

  if(Verbosity) cout << "Run costas array of order " << size << " (alldiff - "
		     << (consistency==FORWARD_CHECKING ? "FC" : "BC" )
		     << "), look for all solutions "; 

  int i, j;
  VarArray X(size, 1, size, domain);
  Solver s;

  s.add( AllDiff(X) );


  Vector< Variable > *distance = new Vector< Variable >[size-2];
  for(i=1; i<size-1; ++i) {
    for(j=0; j<size-i; ++j) {
      distance[i-1].add(X[j] - X[j+i]);
    }
    s.add( AllDiff(distance[i-1], consistency) );
  }

  // Vector< Variable > scope;
  // for(i=1; i<size-1; ++i) {
  //   scope.clear();
  //   for(j=0; j<size-i; ++j) {
  //     scope.add(X[j]-X[j+i]);
  //   }
  //   s.add( AllDiff(scope, consistency) );
  // }


  s.consolidate();
  //std::cout << s << std::endl;
  //cout << s << endl;

#ifdef _MONITOR
  s.monitor(X);
  for(i=0; i<size-2; ++i)
    s.monitor(distance[i]);
#endif



  //s.initialise();
  s.initialise_search(X, 
		      new GenericHeuristic< Lexicographic, MinValue >(&s), 
		      new NoRestart());
  int num_solutions = 0;
  while(s.get_next_solution() == SAT) {
    ++num_solutions;
     // for(i=0; i<size; ++i) {
     //   cout << " " << X[i].get_solution_int_value() ;
     // }
     // cout << endl;
  }

  if(Verbosity) cout << "(" << num_solutions << ") " ;
  if(num_solutions != num_sol[size] /*444*/) {
    cout << "Error: wrong number of solutions! (should be " << num_sol[size] << ")" << endl;
    //exit(1);
  }

  delete [] distance;
}

  
CostasNotEqualAllSolutions::CostasNotEqualAllSolutions(const int sz) 
  : UnitTest(LOW, LOW, LOW) { size=sz; }
CostasNotEqualAllSolutions::~CostasNotEqualAllSolutions() {}

void CostasNotEqualAllSolutions::run() {
  if(Verbosity) cout << "Run costas array of order " << size << " (not equals), look for all solutions "; 
		       

  VarArray X(size, 1, size);

  Solver s;
  int i, j, k;
  for(i=0; i<size; ++i)
    for(j=i+1; j<size; ++j)
      s.add( X[i] != X[j] );

  Vector< Variable > *distance =  new Vector< Variable >[size-2];
  for(i=1; i<size-1; ++i) {
    for(j=0; j<size-i; ++j) {
      distance[i-1].add(X[j] - X[j+i]);
    }
    for(j=1; j<size-i; ++j)
      for(k=0; k<j; ++k) {
	s.add( distance[i-1][j] != distance[i-1][k] );
      }
  }

  //cout << s << endl;

  //s.initialise();
  s.consolidate();
  
  //std::cout << s << std::endl;

#ifdef _DEBUG_PRUNING
  s.monitor(X);
  for(i=0; i<size-2; ++i)
    s.monitor(distance[i]);
#endif


  s.initialise_search(X, 
		      new GenericHeuristic< Lexicographic, MinValue >(&s), 
		      new NoRestart());

  int num_solutions = 0;
  while(s.get_next_solution() == SAT) {
//     for(i=0; i<size; ++i) {
//       cout << " " << setw(2) << X[i].get_solution_str_value() ;
//     }
//     cout << endl;
    ++num_solutions;

    //    if(num_solutions == 2) exit(1);
  }

  if(Verbosity) cout << "(" << num_solutions << ") " ;
  if(num_solutions != num_sol[size] /*444*/) {
    cout << "Error: wrong number of solutions! (should be " << num_sol[size] << ")" << endl;
    //cout << "Error: wrong number of solutions!" << endl;
    //exit(1);
  }

  delete [] distance;
}


Reset::Reset() 
  : UnitTest(LOW, LOW, LOW) { }
Reset::~Reset() {}

void Reset::run() {
    
  VarArray X(3,0,4);

  Solver s;
    
  for(int i=0; i<2; ++i) {
    s.add(X[i] < X[i+1]);
  }

  //s.initialise();

  for(int i=0; i<3; ++i) 
    std::cout << X[i] << " in " << (X[i].get_domain()) << std::endl;

  std::cout << s << std::endl;
    
  s.propagate();

  for(int i=0; i<3; ++i) 
    std::cout << X[i] << " in " << (X[i].get_domain()) << std::endl;
    
  std::cout << s << std::endl;

  s.restore();

  for(int i=0; i<2; ++i) {
    s.add(X[i] < X[i+1]);
  }

  for(int i=0; i<3; ++i) 
    std::cout << X[i] << " in " << (X[i].get_domain()) << std::endl;

    
  std::cout << s << std::endl;

  s.initialise_search(X);

  //     std::cout << s.get_next_solution() << std::endl;

  while(s.get_next_solution() == SAT) {

    //s.depth_first_search(X);

    std::cout << X[0].get_solution_str_value() << "<"
	      << X[1].get_solution_str_value() << "<"
	      << X[2].get_solution_str_value() << std::endl;

  }
    
    

}
  
Pigeons::Pigeons(const int sz) 
  : UnitTest(LOW, LOW, LOW) { size=sz; }
Pigeons::~Pigeons() {}

void Pigeons::run() {
  if(Verbosity) cout << "Run pigeon-hole of order " << size << ": "; 

  VarArray X(size, 1, size-1);
  Solver s;
  int i, j;
  for(i=0; i<size; ++i)
    for(j=i+1; j<size; ++j)
      s.add( X[i] != X[j] );

  s.consolidate();
  //std::cout << s << std::endl;
  //s.initialise();

  s.depth_first_search(X, 
		       new GenericHeuristic< NoOrder, MinValue >(&s), 
		       new NoRestart);
  
  if((int)s.statistics.num_backtracks != num_bts[size] /*362879*/) {
    cout << "Error: wrong number of backtracks! (" 
	 << (s.statistics.num_backtracks) << ")" << endl;
    exit(1);
  }

}


BoolPigeons::BoolPigeons(const int sz, const int vt) 
  : UnitTest(LOW, LOW, LOW) { size=sz; var_type=vt; }
BoolPigeons::~BoolPigeons() {}

void BoolPigeons::run() {
  if(Verbosity) cout << "Run pigeon-hole (" 
		     << (var_type == BITSET_VAR ? "bitset" : "boolean")
		     << ") of order " << size << ": "; 

  cout.flush();

  VarArray X(size*(size-1), 0, 1, var_type);
  VarArray domain;
  Solver s;
  int i, j, k;
  for(i=0; i<size; ++i) {
    domain.clear();
    for(j=i*(size-1); j<(i+1)*(size-1); ++j) {
      domain.add( X[j] );
    }

    s.add( BoolSum(domain, size-2, size-2) );
    //s.add( BoolSum(X[i*(size-1), (i+1)*(size-1)], size-2) );
    for(j=i+1; j<size; ++j) {
      for(k=0; k<size-1; ++k) {
	s.add( X[i*(size-1)+k] || X[j*(size-1)+k] );
      }
    }
  }
  //s.initialise();

  //cout << s << endl;
  
    s.consolidate();
  
    //std::cout << s << std::endl;

#ifdef _MONITOR
   for(i=0; i<size; ++i) {
     for(j=i*(size-1); j<(i+1)*(size-1); ++j) {
       s.monitor_list << X[j];
     }
     s.monitor_list << "\n";
   }
#endif

  //exit(1);

  s.depth_first_search(X, 
		       new GenericHeuristic< Lexicographic, MinValue >(&s), 
		       new NoRestart);
  
  if((int)s.statistics.num_backtracks != num_bts[size] /*40319*/ /*362879*/) {
    cout << "Error: wrong number of backtracks! (" 
	 << (s.statistics.num_backtracks) << ")" << endl;
    //exit(1);
  }

}


// VariableTest::VariableTest() : UnitTest() {}
// VariableTest::~VariableTest() {}

// void VariableTest::run() {

//   Variable X(0,10);
//   Variable Y(-5,5);

//   Solver s;

//   s.add(X);
//   s.add(Y);

//   cout << s << endl;

// }


ModelTest::ModelTest() : UnitTest() {}
ModelTest::~ModelTest() {}

void ModelTest::run() {

  VarArray X(5,1,7);
  Solver s;

  s.add(X);

  cout << s << endl;
  
  s.add( X[0]<=5 );
  s.add( X[1]<=4 );
  s.add( X[2]<=4 );
  s.add( X[3]<=3 );
  s.add( X[4]>=6 );
  s.add( X[0]!=4 );


  cout << s << endl;

  s.add( (X[0]+2)==X[4] || (X[1]+3)==X[4] );

  s.propagate();
  cout << s << endl;

  s.add( (!(X[0]==4) || !(X[1]==4)) == (X[2] >= X[3]) );

  s.propagate();
  cout << s << endl;

  s.add( (X[2] + X[3]) <= X[4] );

//   s.propagate();
//   cout << s << endl;

//   s.add( (X[2] + X[3]) >= X[4] );

  s.propagate();
  cout << s << endl;

  s.save();
  s.rewrite();
  
  cout << s << endl;

  s.restore();

  cout << s << endl;


  s.initialise_search(X, 
		      new GenericHeuristic< Lexicographic, MinValue >(&s), 
		      new NoRestart());
  int num_solutions = 0;
  while(s.get_next_solution() == SAT) {    
    ++num_solutions;
    for(unsigned int i=0; i<5; ++i) {
      cout << " " // << setw(2) << X[i] << ":"
	   << X[i].get_solution_str_value() ;
    }
    cout << endl;
  }

// 

//   cout << s << endl;


}



RewriteTest::RewriteTest() : UnitTest() {}
RewriteTest::~RewriteTest() {}

void RewriteTest::run() {

  VarArray X(4,1,3);
  Solver s;

  s.add(X);

  Variable Y(1,2);

  s.add(Y);


  cout << "Init\n" << s << endl;

  s.add( (Y==1) <= (X[0]==X[2] || X[1]==X[3]) );
  s.add( (Y==1) <= (X[0]!=X[2] || X[1]!=X[3]) );
  s.add( Y!=2 );

  //s.propagate();

  cout << "Add constraints\n" << s << endl;

  s.consolidate();

  cout << "Consolidate\n" << s << endl;  

  s.save();
  s.rewrite();
  
  cout << "Rewrite\n" << s << endl;

  s.save();
  s.add( X[0]+1 < X[1] );

  cout << "Add a new constraint\n" << s << endl;

  s.consolidate();
  s.rewrite();

  cout << "Rewrite\n" << s << endl;

  s.restore();

  cout << "Restore\n" << s << endl;

  s.restore();

  cout << "Restore\n" << s << endl;



//   s.initialise_search(X, 
// 		      new GenericHeuristic< Lexicographic, MinValue >(&s), 
// 		      new NoRestart());
//   int num_solutions = 0;
//   while(s.get_next_solution() == SAT) {    
//     ++num_solutions;
//     for(unsigned int i=0; i<5; ++i) {
//       cout << " " // << setw(2) << X[i] << ":"
// 	   << X[i].get_solution_str_value() ;
//     }
//     cout << endl;
//   }

// 

//   cout << s << endl;


}




SubsetTest::SubsetTest() : UnitTest() {}
SubsetTest::~SubsetTest() {}

void SubsetTest::run() {

  if(Verbosity) cout << "Run Subset test: "; 
  cout.flush();


  run1();
  cout << "1 ";
  cout.flush();

  run2();
  cout << "2 ";
  cout.flush();
  
  run3();
  cout << "3 ";
  cout.flush();
  
  run4();
  cout << "4 ";
  cout.flush();
  
  run5();
  cout << "5 ";
  cout.flush();

  run6();
  cout << "6 ";
  cout.flush();

  run7();
  cout << "7 ";
  cout.flush();

  run8();
  cout << "8 ";
}

// 714 + 81942 = 82656 / 82656
//  50 + 21510 = 21560 / 21560
 

void SubsetTest::run1() {

  Solver s;

  Variable X = SetVariable(1,9,3,5);
  Variable Y = SetVariable(1,9,2,4);

  //std::cout << X << " " << Y << std::endl;

  s.add( Subset(X,Y) );

  // std::cout << s << std::endl;

  s.rewrite();
  
  //cout << "Rewrite\n" << s << endl;

  s.consolidate();

  
  //cout << "Consolidate\n" << s << endl;  


  s.initialise_search(s.variables,
		      new GenericHeuristic< Lexicographic, MinValue >(&s), 
		      new NoRestart());

  int num_solutions = 0;
  while(s.get_next_solution() == SAT) {
    ++num_solutions;
    //cout << X.get_solution_str_value() << " subset of " << Y.get_solution_str_value() << endl;
  }

  if(s.statistics.num_backtracks != 713) {
    cout << "Error: wrong number of backtracks! (" 
	 << (s.statistics.num_backtracks) << ")" << endl;
    //exit(1);
  }
  if(num_solutions != 714) {
    cout << "Error: wrong number of solutions! (" 
	 << (num_solutions) << ")" << endl;
    //exit(1);
  }
}


void SubsetTest::run2() {

  //if(Verbosity) cout << "Run Subset test: "; 

  Solver s;

  Vector<int> lbx;
  Vector<int> ubx;
  Vector<int> lby;
  Vector<int> uby;
  
  ubx.add(-3);
  ubx.add(2);
  ubx.add(5);
  ubx.add(10);
  ubx.add(100);
  ubx.add(123);
  ubx.add(1000);
  ubx.add(1001);

  lbx.add(2);
  lbx.add(100);


  uby.add(-3);
  uby.add(2);
  uby.add(7);
  uby.add(10);
  uby.add(12);
  uby.add(15);
  uby.add(100);
  uby.add(123);
  uby.add(300);
  uby.add(1000);
  uby.add(1001);

  lby.add(300);



  Variable X = SetVariable(lbx,ubx,3,6);
  Variable Y = SetVariable(lby,uby,2,5);

  //std::cout << X << " " << Y << std::endl;

  s.add( Subset(X,Y) );

  //std::cout << s << std::endl;

  s.rewrite();
  
  //cout << "Rewrite\n" << s << endl;

  s.consolidate();

  
  //cout << "Consolidate\n" << s << endl;  


  s.initialise_search(s.variables,
		      new GenericHeuristic< Lexicographic, MinValue >(&s), 
		      new NoRestart());

  int num_solutions = 0;
  while(s.get_next_solution() == SAT) {
    ++num_solutions;
    //cout << X.get_solution_str_value() << " subset of " << Y.get_solution_str_value() << endl;
  }

  //if(s.statistics.num_backtracks != 52) {
  if(s.statistics.num_backtracks != 49) {
    cout << "Error: wrong number of backtracks! (" 
	 << (s.statistics.num_backtracks) << ")" << endl;
    //exit(1);
  }
  if(num_solutions != 50) {
    cout << "Error: wrong number of solutions! (" 
	 << (num_solutions) << ")" << endl;
    //exit(1);
  }
}


void SubsetTest::run3() {

  Solver s;

  Variable X = SetVariable(1,9,3,5);
  Variable Y = SetVariable(1,9,2,4);
  Variable B(0,1);

  //std::cout << X << " " << Y << std::endl;

  s.add( Subset(X,Y) == B );

  // std::cout << s << std::endl;

  s.rewrite();
  
  //cout << "Rewrite\n" << s << endl;

  s.consolidate();

  
  //cout << "Consolidate\n" << s << endl;  


  s.initialise_search(s.variables,
		      new GenericHeuristic< Lexicographic, MinValue >(&s), 
		      new NoRestart());

  int num_solutions = 0;
  while(s.get_next_solution() == SAT) {
    ++num_solutions;
    //cout << X.get_solution_str_value() << " subset of " << Y.get_solution_str_value() << ": " << B.get_solution_str_value() << endl;
  }

  if(s.statistics.num_backtracks != 82655) {
    cout << "Error: wrong number of backtracks! (" 
	 << (s.statistics.num_backtracks) << ")" << endl;
    //exit(1);
  }
  if(num_solutions != 82656) {
    cout << "Error: wrong number of solutions! (" 
	 << (num_solutions) << ")" << endl;
    //exit(1);
  }
}


void SubsetTest::run4() {

  //if(Verbosity) cout << "Run Subset test: "; 

  Solver s;

  Vector<int> lbx;
  Vector<int> ubx;
  Vector<int> lby;
  Vector<int> uby;
  
  ubx.add(-3);
  ubx.add(2);
  ubx.add(5);
  ubx.add(10);
  ubx.add(100);
  ubx.add(123);
  ubx.add(1000);
  ubx.add(1001);

  lbx.add(2);
  lbx.add(100);


  uby.add(-3);
  uby.add(2);
  uby.add(7);
  uby.add(10);
  uby.add(12);
  uby.add(15);
  uby.add(100);
  uby.add(123);
  uby.add(300);
  uby.add(1000);
  uby.add(1001);

  lby.add(300);



  Variable X = SetVariable(lbx,ubx,3,6);
  Variable Y = SetVariable(lby,uby,2,5);
  Variable B(0,1);

  //std::cout << X << " " << Y << std::endl;

  s.add( Subset(X,Y) == B );

  //std::cout << s << std::endl;

  s.rewrite();
  
  //cout << "Rewrite\n" << s << endl;

  s.consolidate();

  
  //cout << "Consolidate\n" << s << endl;  


  s.initialise_search(s.variables,
		      new GenericHeuristic< Lexicographic, MinValue >(&s), 
		      new NoRestart());

  int num_solutions = 0;
  while(s.get_next_solution() == SAT) {
    ++num_solutions;
    //cout << X.get_solution_str_value() << " subset of " << Y.get_solution_str_value() << ": " << B.get_solution_str_value() << endl;
  }

  //if(s.statistics.num_backtracks != 52) {
  if(s.statistics.num_backtracks != 21559) {
    cout << "Error: wrong number of backtracks! (" 
	 << (s.statistics.num_backtracks) << ")" << endl;
    //exit(1);
  }
  if(num_solutions != 21560) {
    cout << "Error: wrong number of solutions! (" 
	 << (num_solutions) << ")" << endl;
    //exit(1);
  }
}


void SubsetTest::run5() {

  Solver s;

  Variable X = SetVariable(1,9,3,5);
  Variable Y = SetVariable(1,9,2,4);

  //std::cout << X << " " << Y << std::endl;

  s.add( NotSubset(X,Y) );

  // std::cout << s << std::endl;

  s.rewrite();
  
  //cout << "Rewrite\n" << s << endl;

  s.consolidate();

  
  //cout << "Consolidate\n" << s << endl;  


  s.initialise_search(s.variables,
		      new GenericHeuristic< Lexicographic, MinValue >(&s), 
		      new NoRestart());

  int num_solutions = 0;
  while(s.get_next_solution() == SAT) {
    ++num_solutions;
    //cout << X.get_solution_str_value() << " not subset of " << Y.get_solution_str_value() << endl;
  }

  if(s.statistics.num_backtracks != 81941) {
    cout << "Error: wrong number of backtracks! (" 
	 << (s.statistics.num_backtracks) << ")" << endl;
    //exit(1);
  }
  if(num_solutions != 81942) {
    cout << "Error: wrong number of solutions! (" 
	 << (num_solutions) << ")" << endl;
    //exit(1);
  }
}


void SubsetTest::run6() {

  //if(Verbosity) cout << "Run Subset test: "; 

  Solver s;

  Vector<int> lbx;
  Vector<int> ubx;
  Vector<int> lby;
  Vector<int> uby;
  
  ubx.add(-3);
  ubx.add(2);
  ubx.add(5);
  ubx.add(10);
  ubx.add(100);
  ubx.add(123);
  ubx.add(1000);
  ubx.add(1001);

  lbx.add(2);
  lbx.add(100);


  uby.add(-3);
  uby.add(2);
  uby.add(7);
  uby.add(10);
  uby.add(12);
  uby.add(15);
  uby.add(100);
  uby.add(123);
  uby.add(300);
  uby.add(1000);
  uby.add(1001);

  lby.add(300);



  Variable X = SetVariable(lbx,ubx,3,6);
  Variable Y = SetVariable(lby,uby,2,5);

  //std::cout << X << " " << Y << std::endl;

  s.add( NotSubset(X,Y) );

  //std::cout << s << std::endl;

  s.rewrite();
  
  //cout << "Rewrite\n" << s << endl;

  s.consolidate();

  
  //cout << "Consolidate\n" << s << endl;  


  s.initialise_search(s.variables,
		      new GenericHeuristic< Lexicographic, MinValue >(&s), 
		      new NoRestart());

  int num_solutions = 0;
  while(s.get_next_solution() == SAT) {
    ++num_solutions;
    //cout << X.get_solution_str_value() << " not subset of " << Y.get_solution_str_value() << endl;
  }

  //if(s.statistics.num_backtracks != 52) {
  if(s.statistics.num_backtracks != 21509) {
    cout << "Error: wrong number of backtracks! (" 
	 << (s.statistics.num_backtracks) << ")" << endl;
    //exit(1);
  }
  if(num_solutions != 21510) {
    cout << "Error: wrong number of solutions! (" 
	 << (num_solutions) << ")" << endl;
    //exit(1);
  }
}


void SubsetTest::run7() {

  Solver s;

  Variable X = SetVariable(1,9,3,5);
  Variable Y = SetVariable(1,9,2,4);
  Variable B(0,1);

  //std::cout << X << " " << Y << std::endl;

  s.add( NotSubset(X,Y) == B );

  // std::cout << s << std::endl;

  s.rewrite();
  
  //cout << "Rewrite\n" << s << endl;

  s.consolidate();

  //cout << "Consolidate\n" << s << endl;  

  //exit(1);

  s.initialise_search(s.variables,
		      new GenericHeuristic< Lexicographic, MinValue >(&s), 
		      new NoRestart());

  int num_solutions = 0;
  while(s.get_next_solution() == SAT) {
    ++num_solutions;
    //cout << X.get_solution_str_value() << " not subset of " << Y.get_solution_str_value() << ": " << B.get_solution_str_value() << endl;
  }

  if(s.statistics.num_backtracks != 82655) {
    cout << "Error: wrong number of backtracks! (" 
	 << (s.statistics.num_backtracks) << ")" << endl;
    //exit(1);
  }
  if(num_solutions != 82656) {
    cout << "Error: wrong number of solutions! (" 
	 << (num_solutions) << ")" << endl;
    //exit(1);
  }
}


void SubsetTest::run8() {

  //if(Verbosity) cout << "Run Subset test: "; 

  Solver s;

  Vector<int> lbx;
  Vector<int> ubx;
  Vector<int> lby;
  Vector<int> uby;
  
  ubx.add(-3);
  ubx.add(2);
  ubx.add(5);
  ubx.add(10);
  ubx.add(100);
  ubx.add(123);
  ubx.add(1000);
  ubx.add(1001);

  lbx.add(2);
  lbx.add(100);


  uby.add(-3);
  uby.add(2);
  uby.add(7);
  uby.add(10);
  uby.add(12);
  uby.add(15);
  uby.add(100);
  uby.add(123);
  uby.add(300);
  uby.add(1000);
  uby.add(1001);

  lby.add(300);



  Variable X = SetVariable(lbx,ubx,3,6);
  Variable Y = SetVariable(lby,uby,2,5);
  Variable B(0,1);

  //std::cout << X << " " << Y << std::endl;

  s.add( NotSubset(X,Y) == B );

  //std::cout << s << std::endl;

  s.rewrite();
  
  //cout << "Rewrite\n" << s << endl;

  s.consolidate();

  
  //cout << "Consolidate\n" << s << endl;  


  s.initialise_search(s.variables,
		      new GenericHeuristic< Lexicographic, MinValue >(&s), 
		      new NoRestart());

  int num_solutions = 0;
  while(s.get_next_solution() == SAT) {
    ++num_solutions;
    //cout << X.get_solution_str_value() << " not subset of " << Y.get_solution_str_value() << ": " << B.get_solution_str_value() << endl;
  }

  //if(s.statistics.num_backtracks != 52) {
  if(s.statistics.num_backtracks != 21559) {
    cout << "Error: wrong number of backtracks! (" 
	 << (s.statistics.num_backtracks) << ")" << endl;
    //exit(1);
  }
  if(num_solutions != 21560) {
    cout << "Error: wrong number of solutions! (" 
	 << (num_solutions) << ")" << endl;
    //exit(1);
  }
}



EqualSetTest::EqualSetTest() : UnitTest() {}
EqualSetTest::~EqualSetTest() {}

void EqualSetTest::run() {

  if(Verbosity) cout << "Run EqualSet test: "; 
  cout.flush();

  run1();
  cout << "1 ";
  cout.flush();

  run2();
  cout << "2 ";
  cout.flush();

  run3();
  cout << "3 ";
  cout.flush();

  run4();
  cout << "4 ";
  cout.flush();

  run5();
  cout << "5 ";
  cout.flush();

  run6();
  cout << "6 ";
  //cout.flush();

}

void EqualSetTest::run1() {
  
  Solver s;

  Variable X = SetVariable(1,9,3,5);
  Variable Y = SetVariable(1,9,2,4);

  //std::cout << X << " " << Y << std::endl;

  s.add( X == Y );

  // std::cout << s << std::endl;

  s.rewrite();
  
  //cout << "Rewrite\n" << s << endl;

  s.consolidate();

  
  //cout << "Consolidate\n" << s << endl;  


  s.initialise_search(s.variables,
		      new GenericHeuristic< Lexicographic, MinValue >(&s), 
		      new NoRestart());

  int num_solutions = 0;
  while(s.get_next_solution() == SAT) {
    ++num_solutions;
    //cout << X.get_solution_str_value() << " == " << Y.get_solution_str_value() << endl;
  }

  if(s.statistics.num_backtracks != 209) {
    cout << "Error: wrong number of backtracks! (" 
	 << (s.statistics.num_backtracks) << ")" << endl;
    //exit(1);
  }
  if(num_solutions != 210) {
    cout << "Error: wrong number of solutions! (" 
	 << (num_solutions) << ")" << endl;
    //exit(1);
  }
}


void EqualSetTest::run2() {

  //if(Verbosity) cout << "Run EqualSet test: "; 

  Solver s;

  Vector<int> lbx;
  Vector<int> ubx;
  Vector<int> lby;
  Vector<int> uby;
  
  ubx.add(-3);
  ubx.add(2);
  ubx.add(5);
  ubx.add(10);
  ubx.add(100);
  ubx.add(123);
  ubx.add(1000);
  ubx.add(1001);

  lbx.add(2);
  lbx.add(100);


  uby.add(-3);
  uby.add(2);
  uby.add(7);
  uby.add(10);
  uby.add(12);
  uby.add(15);
  uby.add(100);
  uby.add(123);
  uby.add(300);
  uby.add(1000);
  uby.add(1001);

  lby.add(10);
  lby.add(100);



  Variable X = SetVariable(lbx,ubx,3,6);
  Variable Y = SetVariable(lby,uby,2,5);

  //std::cout << X << " " << Y << std::endl;

  s.add( (X == Y) );

  //s.add( X != Y );

  //std::cout << s << std::endl;

  s.rewrite();
  
  //cout << "Rewrite\n" << s << endl;

  s.consolidate();

  
  //cout << "Consolidate\n" << s << endl;  


  s.initialise_search(s.variables,
		      new GenericHeuristic< Lexicographic, MinValue >(&s), 
		      new NoRestart());

  int num_solutions = 0;
  while(s.get_next_solution() == SAT) {
    ++num_solutions;
    // cout << X.get_solution_str_value() << " != {" << Y.get_solution_str_value() << "}: " 
    // 	 << (B.get_solution_int_value() ? "true" : "false") << endl;
  }

  //if(s.statistics.num_backtracks != 52) {
  if(s.statistics.num_backtracks != 10) {
    cout << "Error: wrong number of backtracks! (" 
	 << (s.statistics.num_backtracks) << ")" << endl;
    //exit(1);
  }
  if(num_solutions != 11) {
    cout << "Error: wrong number of solutions! (" 
	 << (num_solutions) << ")" << endl;
    //exit(1);
  }
}




void EqualSetTest::run3() {

  //if(Verbosity) cout << "Run EqualSet test: "; 

  Solver s;

  Vector<int> lbx;
  Vector<int> ubx;
  Vector<int> lby;
  Vector<int> uby;
  
  ubx.add(-3);
  ubx.add(2);
  ubx.add(5);
  ubx.add(10);
  ubx.add(100);
  ubx.add(123);
  ubx.add(1000);
  ubx.add(1001);

  lbx.add(2);
  lbx.add(100);


  uby.add(-3);
  uby.add(2);
  uby.add(7);
  uby.add(10);
  uby.add(12);
  uby.add(15);
  uby.add(100);
  uby.add(123);
  uby.add(300);
  uby.add(1000);
  uby.add(1001);

  lby.add(10);
  lby.add(100);



  Variable X = SetVariable(lbx,ubx,3,6);
  Variable Y = SetVariable(lby,uby,2,5);
  Variable B(0, 1);

  //std::cout << X << " " << Y << std::endl;

  s.add( (X == Y) == B );

  //s.add( X != Y );

  //std::cout << s << std::endl;

  s.rewrite();
  
  //cout << "Rewrite\n" << s << endl;

  s.consolidate();

  
  //cout << "Consolidate\n" << s << endl;  


  s.initialise_search(s.variables,
		      new GenericHeuristic< Lexicographic, MinValue >(&s), 
		      new NoRestart());

  int num_solutions = 0;
  while(s.get_next_solution() == SAT) {
    ++num_solutions;
    // cout << X.get_solution_str_value() << " != {" << Y.get_solution_str_value() << "}: " 
    // 	 << (B.get_solution_int_value() ? "true" : "false") << endl;
  }

  //if(s.statistics.num_backtracks != 52) {
  if(s.statistics.num_backtracks != 7279) {
    cout << "Error: wrong number of backtracks! (" 
	 << (s.statistics.num_backtracks) << ")" << endl;
    //exit(1);
  }
  if(num_solutions != 7280) {
    cout << "Error: wrong number of solutions! (" 
	 << (num_solutions) << ")" << endl;
    //exit(1);
  }
}


void EqualSetTest::run4() {
  
  Solver s;

  Variable X = SetVariable(1,9,3,5);
  Variable Y = SetVariable(1,9,2,4);

  //std::cout << X << " " << Y << std::endl;

  s.add( X != Y );

  // std::cout << s << std::endl;

  s.rewrite();
  
  //cout << "Rewrite\n" << s << endl;

  s.consolidate();

  
  //cout << "Consolidate\n" << s << endl;  


  s.initialise_search(s.variables,
		      new GenericHeuristic< Lexicographic, MinValue >(&s), 
		      new NoRestart());

  int num_solutions = 0;
  while(s.get_next_solution() == SAT) {
    ++num_solutions;
    //cout << X.get_solution_str_value() << " == " << Y.get_solution_str_value() << endl;
  }

  if(s.statistics.num_backtracks != 82515) {
    cout << "Error: wrong number of backtracks! (" 
	 << (s.statistics.num_backtracks) << ")" << endl;
    //exit(1);
  }
  if(num_solutions != 82446) {
    cout << "Error: wrong number of solutions! (" 
	 << (num_solutions) << ")" << endl;
    //exit(1);
  }
}


void EqualSetTest::run5() {

  //if(Verbosity) cout << "Run EqualSet test: "; 

  Solver s;

  Vector<int> lbx;
  Vector<int> ubx;
  Vector<int> lby;
  Vector<int> uby;
  
  ubx.add(-3);
  ubx.add(2);
  ubx.add(5);
  ubx.add(10);
  ubx.add(100);
  ubx.add(123);
  ubx.add(1000);
  ubx.add(1001);

  lbx.add(2);
  lbx.add(100);


  uby.add(-3);
  uby.add(2);
  uby.add(7);
  uby.add(10);
  uby.add(12);
  uby.add(15);
  uby.add(100);
  uby.add(123);
  uby.add(300);
  uby.add(1000);
  uby.add(1001);

  lby.add(10);
  lby.add(100);



  Variable X = SetVariable(lbx,ubx,3,6);
  Variable Y = SetVariable(lby,uby,2,5);

  //std::cout << X << " " << Y << std::endl;

  s.add( (X != Y) );

  //s.add( X != Y );

  //std::cout << s << std::endl;

  s.rewrite();
  
  //cout << "Rewrite\n" << s << endl;

  s.consolidate();

  
  //cout << "Consolidate\n" << s << endl;  


  s.initialise_search(s.variables,
		      new GenericHeuristic< Lexicographic, MinValue >(&s), 
		      new NoRestart());

  int num_solutions = 0;
  while(s.get_next_solution() == SAT) {
    ++num_solutions;
    // cout << X.get_solution_str_value() << " != {" << Y.get_solution_str_value() << "}: " 
    // 	 << (B.get_solution_int_value() ? "true" : "false") << endl;
  }

  //if(s.statistics.num_backtracks != 52) {
  if(s.statistics.num_backtracks != 7271) {
    cout << "Error: wrong number of backtracks! (" 
	 << (s.statistics.num_backtracks) << ")" << endl;
    //exit(1);
  }
  if(num_solutions != 7269) {
    cout << "Error: wrong number of solutions! (" 
	 << (num_solutions) << ")" << endl;
    //exit(1);
  }
}




void EqualSetTest::run6() {

  //if(Verbosity) cout << "Run EqualSet test: "; 

  Solver s;

  Vector<int> lbx;
  Vector<int> ubx;
  Vector<int> lby;
  Vector<int> uby;
  
  ubx.add(-3);
  ubx.add(2);
  ubx.add(5);
  ubx.add(10);
  ubx.add(100);
  ubx.add(123);
  ubx.add(1000);
  ubx.add(1001);

  lbx.add(2);
  lbx.add(100);


  uby.add(-3);
  uby.add(2);
  uby.add(7);
  uby.add(10);
  uby.add(12);
  uby.add(15);
  uby.add(100);
  uby.add(123);
  uby.add(300);
  uby.add(1000);
  uby.add(1001);

  lby.add(10);
  lby.add(100);



  Variable X = SetVariable(lbx,ubx,3,6);
  Variable Y = SetVariable(lby,uby,2,5);
  Variable B(0, 1);

  //std::cout << X << " " << Y << std::endl;

  s.add( (X != Y) == B );

  //s.add( X != Y );

  //std::cout << s << std::endl;

  s.rewrite();
  
  //cout << "Rewrite\n" << s << endl;

  s.consolidate();

  
  //cout << "Consolidate\n" << s << endl;  


  s.initialise_search(s.variables,
		      new GenericHeuristic< Lexicographic, MinValue >(&s), 
		      new NoRestart());

  int num_solutions = 0;
  while(s.get_next_solution() == SAT) {
    ++num_solutions;
    // cout << X.get_solution_str_value() << " != {" << Y.get_solution_str_value() << "}: " 
    // 	 << (B.get_solution_int_value() ? "true" : "false") << endl;
  }

  //if(s.statistics.num_backtracks != 52) {
  if(s.statistics.num_backtracks != 7279) {
    cout << "Error: wrong number of backtracks! (" 
	 << (s.statistics.num_backtracks) << ")" << endl;
    //exit(1);
  }
  if(num_solutions != 7280) {
    cout << "Error: wrong number of solutions! (" 
	 << (num_solutions) << ")" << endl;
    //exit(1);
  }
}




MemberTest::MemberTest() : UnitTest() {}
MemberTest::~MemberTest() {}

void MemberTest::run() {
  if(Verbosity) cout << "Run Member test: "; 
  run1();
  cout << "1 ";
  run2();
  cout << "2 ";
}

void MemberTest::run1() {

  //if(Verbosity) cout << "Run Member test: "; 

  Solver s;

  Variable X = Variable(0,10);
  Variable Y = SetVariable(1,9,0,4);

  //std::cout << X << " " << Y << std::endl;

  s.add( Member(X,Y) );

  //std::cout << s << std::endl;

  s.rewrite();
  
  //cout << "Rewrite\n" << s << endl;

  s.consolidate();
  
  //cout << "Consolidate\n" << s << endl;  

#ifdef _MONITOR
  
  s.monitor_list << X ;
  
  s.monitor_list << " in " ;
  
  s.monitor_list << Y << "\n";

  /*
  s.monitor( X );
  s.monitor( " in ");
  s.monitor( Y );
  s.monitor( "\n");
  */
#endif

  s.initialise_search(s.variables,
		      new GenericHeuristic< Lexicographic, MinValue >(&s), 
		      new NoRestart());

  int num_solutions = 0;
  while(s.get_next_solution() == SAT) {
    ++num_solutions;
    //cout << X.get_solution_str_value() << " member of " << Y.get_solution_str_value() << endl;
  }

  if(s.statistics.num_backtracks != 836) {
    cout << "Error: wrong number of backtracks! (" 
	 << (s.statistics.num_backtracks) << ")" << endl;
    //exit(1);
  }
  if(num_solutions != 837) {
    cout << "Error: wrong number of solutions! (" 
	 << (num_solutions) << ")" << endl;
    //exit(1);
  }
}


void MemberTest::run2() {

  //if(Verbosity) cout << "Run Member test: "; 

  Solver s;

  Variable X = Variable(0,10);
  Variable Y = SetVariable(1,9,0,4);
  Variable b(0,1);

  //std::cout << X << " " << Y << std::endl;

  s.add( Member(X,Y) != b );

  //std::cout << s << std::endl;

  s.rewrite();
  
  // cout << "Rewrite\n" << s << endl;

  s.consolidate();
  
  //cout << "Consolidate\n" << s << endl;  

  //cout << s.variables << endl;

#ifdef _MONITOR
  
  s.monitor_list << X ;
  
  s.monitor_list << " in " ;
  
  s.monitor_list << Y ;

  s.monitor_list << " <-/-> " ;
  
  s.monitor_list << b ;
  
  s.monitor_list << "\n";

  /*
  s.monitor( X );
  s.monitor( " in ");
  s.monitor( Y );
  s.monitor( "\n");
  */
#endif


  s.initialise_search(s.variables,
		      new GenericHeuristic< Lexicographic, MinValue >(&s), 
		      new NoRestart());


  int num_solutions = 0;
  while(s.get_next_solution() == SAT) {
    ++num_solutions;
    // cout << X.get_solution_str_value() << " member of " << Y.get_solution_str_value() 
    // 	 << (b.get_solution_int_value() ? " (true)" : " (false)") << endl;
  }

  //if(s.statistics.num_backtracks != 2906) {
  if(s.statistics.num_backtracks != 2815) {
    cout << "Error: wrong number of backtracks! (" 
	 << (s.statistics.num_backtracks) << ")" << endl;
    //exit(1);
  }
  if(num_solutions != 2816) {
    cout << "Error: wrong number of solutions! (" 
	 << (num_solutions) << ")" << endl;
    //exit(1);
  }
}



SatTest::SatTest() : UnitTest() {}
SatTest::~SatTest() {}

void SatTest::run() {

  if(Verbosity) cout << "Run Sat test: "; 

  Solver solver;

  solver.parse_dimacs("cnf/gen-1.3/unif-c1000-v250-s542677735.cnf");
  //solver.parameters.verbosity = 2;
  solver.parameters.backjump = 1;
  
  std::cout << solver << std::endl;
  

  solver.depth_first_search(solver.variables,
			    new GenericHeuristic< VSIDS<4>, Guided<RandomMinMax> >(&solver),
			    new Geometric()
			    );
  
  if(solver.statistics.num_backtracks != 54349) {
    cout << "Error: wrong number of backtracks! (" 
	 << (solver.statistics.num_backtracks) << ")" << endl;
    //exit(1);
  }
}


DivTest::DivTest() : UnitTest() {}
DivTest::~DivTest() {}

void DivTest::run() {

  if(Verbosity) cout << "Run Div test: "; 

  Solver s;

  Variable X(-450, 1834);
  Variable Y(-12, 7);

  Vector<int> vals;
  vals.add(-45);
  vals.add(-23);
  vals.add(-5);
  vals.add(-4);
  vals.add(-2);
  vals.add(0);
  vals.add(3);
  vals.add(4);
  vals.add(6);
  vals.add(12);
  vals.add(20);


  Variable Z(vals);

  s.add(X/Y == Z);

  s.rewrite();
  s.consolidate();

  //std::cout << s << std::endl;

  s.initialise_search(s.variables,
		      new GenericHeuristic< Lexicographic, MinValue >(&s), 
		      new NoRestart());

  int num_solutions = 0;
  while(s.get_next_solution() == SAT) {
    Solution sol(s.variables);
    ++num_solutions;
    
    int x = X.get_solution_int_value();
    int y = Y.get_solution_int_value();
    int z = Z.get_solution_int_value();

    /*
    std::cout <<  x
	      << "/" << y 
	      << "=" << z
	      << std::endl;
    */

    if(x/y != z) {
      std::cout << "error!" << std::endl;
      exit(1);
    }
  }

  if(s.statistics.num_backtracks != 1473) {
    cout << "Error: wrong number of backtracks! (" 
	 << (s.statistics.num_backtracks) << ")" << endl;
    //exit(1);
  }
  if(num_solutions != 1253) {
    cout << "Error: wrong number of solutions! (" 
	 << (num_solutions) << ")" << endl;
    //exit(1);
  }
}


CardTest::CardTest() : UnitTest() {}
CardTest::~CardTest() {}

void CardTest::run() {

  if(Verbosity) cout << "Run Card test: "; 

  Solver s;

  Variable X = Variable(1,8);
  Variable Y = SetVariable(1,9,0,9);


  s.add( Card(Y) == X );

  //std::cout << s << std::endl;

  s.rewrite();
  
  //cout << "Rewrite\n" << s << endl;

  s.consolidate();
  
  //cout << "Consolidate\n" << s << endl;  


  s.initialise_search(s.variables,
		      new GenericHeuristic< Lexicographic, MinValue >(&s), 
		      new NoRestart());

  int num_solutions = 0;
  while(s.get_next_solution() == SAT) {
    ++num_solutions;
    //cout << X.get_solution_str_value() << " == |" << Y.get_solution_str_value() << "|" << endl;
  }

  if(s.statistics.num_backtracks != 509) {
    cout << "Error: wrong number of backtracks! (" 
	 << (s.statistics.num_backtracks) << ")" << endl;
    //exit(1);
  }
  if(num_solutions != 510) {
    cout << "Error: wrong number of solutions! (" 
	 << (num_solutions) << ")" << endl;
    //exit(1);
  }
}



IntersectionTest::IntersectionTest() : UnitTest() {}
IntersectionTest::~IntersectionTest() {}

void IntersectionTest::run() {

  if(Verbosity) cout << "Run Intersection test: "; 
  cout.flush();

  run1();
  cout << "1 ";
  cout.flush();

  run2();
  cout << "2 ";
  cout.flush();

  run3();
  cout << "3 ";

}


void IntersectionTest::run1() {

  Solver s;

  Variable Z = SetVariable(1,9,2,3);
  Variable Y = SetVariable(1,9,1,5);
  Variable X = SetVariable(1,9,1,5);


  s.add( Intersection(X,Y) == Z );

  //std::cout << s << std::endl;

  s.rewrite();
  
  //cout << "Rewrite\n" << s << endl;

  s.consolidate();
  
  //cout << "Consolidate\n" << s << endl;  

  //cout << s << endl;

  s.initialise_search(s.variables,
		      new GenericHeuristic< Lexicographic, MinValue >(&s), 
		      new NoRestart());


  //cout << s.sequence << std::endl;

  //exit(1);

  int num_solutions = 0;
  while(s.get_next_solution() == SAT) {
    ++num_solutions;

    // // //cout << s.last_solution_lb << endl;
    // cout << Z.get_solution_str_value() << " == " 
    //  	 << X.get_solution_str_value() << " inter " 
    //  	 << Y.get_solution_str_value() << endl;

  }

  //if(s.statistics.num_backtracks != 75215) {
  if(s.statistics.num_backtracks != 75881) {
    cout << "Error: wrong number of backtracks! (" 
	 << (s.statistics.num_backtracks) << ")" << endl;
    //exit(1);
  }
  if(num_solutions != 75216) {
    cout << "Error: wrong number of solutions! (" 
	 << (num_solutions) << ")" << endl;
    //exit(1);
  }
}


void IntersectionTest::run2() {

  //if(Verbosity) cout << "Run Subset test: "; 

  Solver s;

  Vector<int> lbx;
  Vector<int> ubx;
  Vector<int> lby;
  Vector<int> uby;
  
  ubx.add(-3);
  ubx.add(2);
  ubx.add(5);
  ubx.add(10);
  ubx.add(100);
  ubx.add(123);
  ubx.add(1000);
  ubx.add(1001);

  lbx.add(2);
  lbx.add(100);


  uby.add(-3);
  uby.add(2);
  uby.add(7);
  uby.add(10);
  uby.add(12);
  uby.add(15);
  uby.add(100);
  uby.add(123);
  uby.add(300);
  uby.add(1000);
  uby.add(1001);

  lby.add(100);
  lby.add(300);



  Variable X = SetVariable(lbx,ubx,3,6);
  Variable Y = SetVariable(lby,uby,2,5);
 

  //std::cout << X << " " << Y << std::endl;

  s.add( Card(Intersection(X,Y)) == 3 );

  //std::cout << s << std::endl;

  s.rewrite();
  
  //cout << "Rewrite\n" << s << endl;

  s.consolidate();

  
  //cout << "Consolidate\n" << s << endl;  


  s.initialise_search(s.variables,
		      new GenericHeuristic< Lexicographic, MinValue >(&s), 
		      new NoRestart());

  int num_solutions = 0;
  while(s.get_next_solution() == SAT) {
    ++num_solutions;
    // cout << "Card of " << X.get_solution_str_value() << " inter " << Y.get_solution_str_value() 
    // 	 << " = 3" << endl;
  }

  //if(s.statistics.num_backtracks != 52) {
  if(s.statistics.num_backtracks != 1899) {
    cout << "Error: wrong number of backtracks! (" 
	 << (s.statistics.num_backtracks) << ")" << endl;
    //exit(1);
  }
  if(num_solutions != 1470) {
    cout << "Error: wrong number of solutions! (" 
	 << (num_solutions) << ")" << endl;
    //exit(1);
  }
}


void IntersectionTest::run3() {

  //if(Verbosity) cout << "Run Subset test: "; 

  Solver s;

  Vector<int> lbx;
  Vector<int> ubx;
  Vector<int> lby;
  Vector<int> uby;
  Vector<int> lbz;
  Vector<int> ubz;
  
  ubx.add(-3);
  ubx.add(2);
  ubx.add(5);
  ubx.add(10);
  ubx.add(100);
  ubx.add(123);
  ubx.add(1000);
  ubx.add(1001);

  lbx.add(2);
  lbx.add(100);


  uby.add(-3);
  uby.add(2);
  uby.add(7);
  uby.add(10);
  uby.add(12);
  uby.add(15);
  uby.add(100);
  uby.add(123);
  uby.add(300);
  uby.add(1000);
  uby.add(1001);

  lby.add(100);
  lby.add(300);

  
  lbz.add(2);
  lbz.add(10);

  ubz.add(-3);
  ubz.add(2);
  ubz.add(5);
  ubz.add(7);
  ubz.add(10);
  ubz.add(12);
  ubz.add(15);
  ubz.add(100);
  ubz.add(123);
  ubz.add(300);
  ubz.add(1000);




  Variable X = SetVariable(lbx,ubx,3,6);
  Variable Y = SetVariable(lby,uby,2,5);
  Variable Z = SetVariable(lbz,ubz,2,5);
 

  //std::cout << X << " " << Y << std::endl;

  s.add( Intersection(X,Y) == Z );

  //std::cout << s << std::endl;

  s.rewrite();
  
  //cout << "Rewrite\n" << s << endl;

  s.consolidate();

  
  //cout << "Consolidate\n" << s << endl;  


  s.initialise_search(s.variables,
		      new GenericHeuristic< Lexicographic, MinValue >(&s), 
		      new NoRestart());

  int num_solutions = 0;
  while(s.get_next_solution() == SAT) {
    ++num_solutions;
    // cout << X.get_solution_str_value() << " inter " << Y.get_solution_str_value() << " = " 
    // 	 << Z.get_solution_str_value() << endl;
  }

  //if(s.statistics.num_backtracks != 52) {
  if(s.statistics.num_backtracks != 196) {
    cout << "Error: wrong number of backtracks! (" 
	 << (s.statistics.num_backtracks) << ")" << endl;
    //exit(1);
  }
  if(num_solutions != 197) {
    cout << "Error: wrong number of solutions! (" 
	 << (num_solutions) << ")" << endl;
    //exit(1);
  }
}



UnionTest::UnionTest() : UnitTest() {}
UnionTest::~UnionTest() {}

void UnionTest::run() {

  if(Verbosity) cout << "Run Union test: "; 
  cout.flush();

  run1();
  cout << "1 ";
  cout.flush();

  run2();
  cout << "2 ";
  cout.flush();

  run3();
  cout << "3 ";

}


void UnionTest::run1() {

  Solver s;

  Variable Z = SetVariable(1,9,2,3);
  Variable Y = SetVariable(1,9,1,5);
  Variable X = SetVariable(1,9,1,5);


  s.add( Union(X,Y) == Z );

  //std::cout << s << std::endl;

  s.rewrite();
  
  //cout << "Rewrite\n" << s << endl;

  s.consolidate();
  
  //cout << "Consolidate\n" << s << endl;  

  //cout << s << endl;

  s.initialise_search(s.variables,
		      new GenericHeuristic< Lexicographic, MinValue >(&s), 
		      new NoRestart());


  //cout << s.sequence << std::endl;

  //exit(1);

  int num_solutions = 0;
  while(s.get_next_solution() == SAT) {
    ++num_solutions;

    // //cout << s.last_solution_lb << endl;
    // cout << Z.get_solution_str_value() << " == " 
    //  	 << X.get_solution_str_value() << " union " 
    //  	 << Y.get_solution_str_value() << endl;

  }

  //if(s.statistics.num_backtracks != 75215) {
  if(s.statistics.num_backtracks != 2351) {
    cout << "Error: wrong number of backtracks! (" 
	 << (s.statistics.num_backtracks) << ")" << endl;
    //exit(1);
  }
  if(num_solutions != 2352) {
    cout << "Error: wrong number of solutions! (" 
	 << (num_solutions) << ")" << endl;
    //exit(1);
  }
}


void UnionTest::run2() {

  //if(Verbosity) cout << "Run Subset test: "; 

  Solver s;

  Vector<int> lbx;
  Vector<int> ubx;
  Vector<int> lby;
  Vector<int> uby;
  
  ubx.add(-3);
  ubx.add(2);
  ubx.add(5);
  ubx.add(10);
  ubx.add(100);
  ubx.add(123);
  ubx.add(1000);
  ubx.add(1001);

  lbx.add(2);
  lbx.add(100);


  uby.add(-3);
  uby.add(2);
  uby.add(7);
  uby.add(10);
  uby.add(12);
  uby.add(15);
  uby.add(100);
  uby.add(123);
  uby.add(300);
  uby.add(1000);
  uby.add(1001);

  lby.add(100);
  lby.add(300);



  Variable X = SetVariable(lbx,ubx,3,6);
  Variable Y = SetVariable(lby,uby,2,5);
 

  //std::cout << X << " " << Y << std::endl;

  s.add( Card(Union(X,Y)) == 5 );

  //std::cout << s << std::endl;

  s.rewrite();
  
  //cout << "Rewrite\n" << s << endl;

  s.consolidate();

  
  //cout << "Consolidate\n" << s << endl;  


  s.initialise_search(s.variables,
		      new GenericHeuristic< Lexicographic, MinValue >(&s), 
		      new NoRestart());

  int num_solutions = 0;
  while(s.get_next_solution() == SAT) {
    ++num_solutions;
    // cout << "Card of " << X.get_solution_str_value() << " union " << Y.get_solution_str_value() 
    //  	 << " = 5" << endl;
  }

  //if(s.statistics.num_backtracks != 52) {
  if(s.statistics.num_backtracks != 255) {
    cout << "Error: wrong number of backtracks! (" 
	 << (s.statistics.num_backtracks) << ")" << endl;
    //exit(1);
  }
  if(num_solutions != 256) {
    cout << "Error: wrong number of solutions! (" 
	 << (num_solutions) << ")" << endl;
    //exit(1);
  }
}


void UnionTest::run3() {

  //if(Verbosity) cout << "Run Subset test: "; 

  Solver s;

  Vector<int> lbx;
  Vector<int> ubx;
  Vector<int> lby;
  Vector<int> uby;
  Vector<int> lbz;
  Vector<int> ubz;
  
  ubx.add(-3);
  ubx.add(2);
  ubx.add(5);
  ubx.add(10);
  ubx.add(100);
  ubx.add(123);
  ubx.add(1000);
  ubx.add(1001);

  lbx.add(2);
  lbx.add(100);


  uby.add(-3);
  uby.add(2);
  uby.add(7);
  uby.add(10);
  uby.add(12);
  uby.add(15);
  uby.add(100);
  uby.add(123);
  uby.add(300);
  uby.add(1000);
  uby.add(1001);

  lby.add(100);
  lby.add(300);

  
  lbz.add(2);
  lbz.add(10);

  ubz.add(-3);
  ubz.add(2);
  ubz.add(5);
  ubz.add(7);
  ubz.add(10);
  ubz.add(12);
  ubz.add(15);
  ubz.add(100);
  ubz.add(123);
  ubz.add(300);
  ubz.add(1000);




  Variable X = SetVariable(lbx,ubx,3,6);
  Variable Y = SetVariable(lby,uby,2,5);
  Variable Z = SetVariable(lbz,ubz,2,5);
 

  //std::cout << X << " " << Y << std::endl;

  s.add( Union(X,Y) == Z );

  //std::cout << s << std::endl;

  s.rewrite();
  
  //cout << "Rewrite\n" << s << endl;

  s.consolidate();

  
  //cout << "Consolidate\n" << s << endl;  


  s.initialise_search(s.variables,
		      new GenericHeuristic< Lexicographic, MinValue >(&s), 
		      new NoRestart());

  int num_solutions = 0;
  while(s.get_next_solution() == SAT) {
    ++num_solutions;
    // cout << X.get_solution_str_value() << " union " << Y.get_solution_str_value() << " = " 
    //  	 << Z.get_solution_str_value() << endl;
  }

  //if(s.statistics.num_backtracks != 52) {
  if(s.statistics.num_backtracks != 69) {
    cout << "Error: wrong number of backtracks! (" 
	 << (s.statistics.num_backtracks) << ")" << endl;
    //exit(1);
  }
  if(num_solutions != 70) {
    cout << "Error: wrong number of solutions! (" 
	 << (num_solutions) << ")" << endl;
    //exit(1);
  }
}




SymmetricDifferenceTest::SymmetricDifferenceTest() : UnitTest() {}
SymmetricDifferenceTest::~SymmetricDifferenceTest() {}

void SymmetricDifferenceTest::run() {

  if(Verbosity) cout << "Run SymmetricDifference test: "; 
  cout.flush();

  run1();
  cout << "1 ";
  cout.flush();

  run2();
  cout << "2 ";
  cout.flush();

  run3();
  cout << "3 ";

}


void SymmetricDifferenceTest::run1() {

  Solver s;

  Variable Z = SetVariable(1,9,2,3);
  Variable Y = SetVariable(1,9,1,5);
  Variable X = SetVariable(1,9,1,5);


  s.add( SymmetricDifference(X,Y) == Z );

  //std::cout << s << std::endl;

  s.rewrite();
  
  //cout << "Rewrite\n" << s << endl;

  s.consolidate();
  
  //cout << "Consolidate\n" << s << endl;  

  //cout << s << endl;

  s.initialise_search(s.variables,
		      new GenericHeuristic< Lexicographic, MinValue >(&s), 
		      new NoRestart());


  //cout << s.sequence << std::endl;

  //exit(1);

  int num_solutions = 0;
  while(s.get_next_solution() == SAT) {
    ++num_solutions;

    // //cout << s.last_solution_lb << endl;
    // cout << Z.get_solution_str_value() << " == " 
    //  	 << X.get_solution_str_value() << " symdiff " 
    //  	 << Y.get_solution_str_value() << endl;

  }

  //if(s.statistics.num_backtracks != 75215) {
  if(s.statistics.num_backtracks != 40174) {
    cout << "Error: wrong number of backtracks! (" 
	 << (s.statistics.num_backtracks) << ")" << endl;
    //exit(1);
  }
  if(num_solutions != 36360) {
    cout << "Error: wrong number of solutions! (" 
	 << (num_solutions) << ")" << endl;
    //exit(1);
  }
}


void SymmetricDifferenceTest::run2() {

  //if(Verbosity) cout << "Run Subset test: "; 

  Solver s;

  Vector<int> lbx;
  Vector<int> ubx;
  Vector<int> lby;
  Vector<int> uby;
  
  ubx.add(-3);
  ubx.add(2);
  ubx.add(5);
  ubx.add(10);
  ubx.add(100);
  ubx.add(123);
  ubx.add(1000);
  ubx.add(1001);

  lbx.add(2);
  lbx.add(100);


  uby.add(-3);
  uby.add(2);
  uby.add(7);
  uby.add(10);
  uby.add(12);
  uby.add(15);
  uby.add(100);
  uby.add(123);
  uby.add(300);
  uby.add(1000);
  uby.add(1001);

  lby.add(100);
  lby.add(300);



  Variable X = SetVariable(lbx,ubx,3,6);
  Variable Y = SetVariable(lby,uby,2,5);
 

  //std::cout << X << " " << Y << std::endl;

  s.add( Card(SymmetricDifference(X,Y)) == 5 );

  //std::cout << s << std::endl;

  s.rewrite();
  
  //cout << "Rewrite\n" << s << endl;

  s.consolidate();

  
  //cout << "Consolidate\n" << s << endl;  


  s.initialise_search(s.variables,
		      new GenericHeuristic< Lexicographic, MinValue >(&s), 
		      new NoRestart());

  int num_solutions = 0;
  while(s.get_next_solution() == SAT) {
    ++num_solutions;
    // cout << "Card of " << X.get_solution_str_value() << " symdiff " << Y.get_solution_str_value() 
    //  	 << " = 5" << endl;
  }

  //if(s.statistics.num_backtracks != 52) {
  if(s.statistics.num_backtracks != 3084) {
    cout << "Error: wrong number of backtracks! (" 
	 << (s.statistics.num_backtracks) << ")" << endl;
    //exit(1);
  }
  if(num_solutions != 1853) {
    cout << "Error: wrong number of solutions! (" 
	 << (num_solutions) << ")" << endl;
    //exit(1);
  }
}


void SymmetricDifferenceTest::run3() {

  //if(Verbosity) cout << "Run Subset test: "; 

  Solver s;

  Vector<int> lbx;
  Vector<int> ubx;
  Vector<int> lby;
  Vector<int> uby;
  Vector<int> lbz;
  Vector<int> ubz;
  
  ubx.add(-3);
  ubx.add(2);
  ubx.add(5);
  ubx.add(10);
  ubx.add(100);
  ubx.add(123);
  ubx.add(1000);
  ubx.add(1001);

  lbx.add(2);
  lbx.add(100);


  uby.add(-3);
  uby.add(2);
  uby.add(7);
  uby.add(10);
  uby.add(12);
  uby.add(15);
  uby.add(100);
  uby.add(123);
  uby.add(300);
  uby.add(1000);
  uby.add(1001);

  lby.add(100);
  lby.add(300);

  
  lbz.add(2);
  lbz.add(10);

  ubz.add(-3);
  ubz.add(2);
  ubz.add(5);
  ubz.add(7);
  ubz.add(10);
  ubz.add(12);
  ubz.add(15);
  ubz.add(100);
  ubz.add(123);
  ubz.add(300);
  ubz.add(1000);




  Variable X = SetVariable(lbx,ubx,3,6);
  Variable Y = SetVariable(lby,uby,2,5);
  Variable Z = SetVariable(lbz,ubz,2,5);
 

  //std::cout << X << " " << Y << std::endl;

  s.add( SymmetricDifference(X,Y) == Z );

  //std::cout << s << std::endl;

  s.rewrite();
  
  //cout << "Rewrite\n" << s << endl;

  s.consolidate();

  
  //cout << "Consolidate\n" << s << endl;  


  s.initialise_search(s.variables,
		      new GenericHeuristic< Lexicographic, MinValue >(&s), 
		      new NoRestart());

  int num_solutions = 0;
  while(s.get_next_solution() == SAT) {
    ++num_solutions;
    // cout << X.get_solution_str_value() << " symdiff " << Y.get_solution_str_value() << " = " 
    //   	 << Z.get_solution_str_value() << endl;
  }

  //if(s.statistics.num_backtracks != 52) {
  if(s.statistics.num_backtracks != 608) {
    cout << "Error: wrong number of backtracks! (" 
	 << (s.statistics.num_backtracks) << ")" << endl;
    //exit(1);
  }
  if(num_solutions != 500) {
    cout << "Error: wrong number of solutions! (" 
	 << (num_solutions) << ")" << endl;
    //exit(1);
  }
}




SetDifferenceTest::SetDifferenceTest() : UnitTest() {}
SetDifferenceTest::~SetDifferenceTest() {}

void SetDifferenceTest::run() {

  if(Verbosity) cout << "Run SetDifference test: "; 

  cout.flush();

  run1();
  cout << "1 ";
  cout.flush();

  run2();
  cout << "2 ";
  cout.flush();

  run3();
  cout << "3 ";

}


void SetDifferenceTest::run1() {
  Solver s;

  Variable Z = SetVariable(1,9,2,3);
  Variable Y = SetVariable(1,9,1,5);
  Variable X = SetVariable(1,9,1,5);


  // Variable Z = SetVariable(1,3,2,3);
  // Variable Y = SetVariable(1,3,1,3);
  // Variable X = SetVariable(1,3,1,3);
  
  // //std::cout << X << " - " << Y << " = " << Z << std::endl;


  s.add( SetDifference(X,Y) == Z );

  //std::cout << s << std::endl;



  s.rewrite();
  
  //cout << "Rewrite\n" << s << endl;

  s.consolidate();
  
  //cout << "Consolidate\n" << s << endl;  

  // cout << s << endl;

  // exit(1);

  s.initialise_search(s.variables,
		      new GenericHeuristic< Lexicographic, MinValue >(&s), 
		      new NoRestart());


  //cout << s.sequence << std::endl;

  //exit(1);

  int num_solutions = 0;
  while(s.get_next_solution() == SAT) {
    ++num_solutions;

    // // //cout << s.last_solution_lb << endl;
    // cout << Z.get_solution_str_value() << " == " 
    //  	 << X.get_solution_str_value() << " \\ " 
    //  	 << Y.get_solution_str_value() << endl;

  }

  //if(s.statistics.num_backtracks != 75215) {
  if(s.statistics.num_backtracks != 92036) {
    cout << "Error: wrong number of backtracks! (" 
	 << (s.statistics.num_backtracks) << ")" << endl;
    //exit(1);
  }
  if(num_solutions != 91896) {
    cout << "Error: wrong number of solutions! (" 
	 << (num_solutions) << ")" << endl;
    //exit(1);
  }
}

void SetDifferenceTest::run2() {

  //if(Verbosity) cout << "Run Subset test: "; 

  Solver s;

  Vector<int> lbx;
  Vector<int> ubx;
  Vector<int> lby;
  Vector<int> uby;
  
  ubx.add(-3);
  ubx.add(2);
  ubx.add(5);
  ubx.add(10);
  ubx.add(100);
  ubx.add(123);
  ubx.add(1000);
  ubx.add(1001);

  lbx.add(2);
  lbx.add(100);


  uby.add(-3);
  uby.add(2);
  uby.add(7);
  uby.add(10);
  uby.add(12);
  uby.add(15);
  uby.add(100);
  uby.add(123);
  uby.add(300);
  uby.add(1000);
  uby.add(1001);

  lby.add(100);
  lby.add(300);



  Variable X = SetVariable(lbx,ubx,3,6);
  Variable Y = SetVariable(lby,uby,2,5);
 

  //std::cout << X << " " << Y << std::endl;

  s.add( Card(SetDifference(X,Y)) == 3 );

  //std::cout << s << std::endl;

  s.rewrite();
  
  //cout << "Rewrite\n" << s << endl;

  s.consolidate();


  //cout << "|" << X.get_domain() << " \\ " << Y.get_domain() << "| = 3" << endl;

  
  //cout << "Consolidate\n" << s << endl;  


  s.initialise_search(s.variables,
		      new GenericHeuristic< Lexicographic, MinValue >(&s), 
		      new NoRestart());

  int num_solutions = 0;
  while(s.get_next_solution() == SAT) {
    ++num_solutions;
    // cout << "Card of " << X.get_solution_str_value() << " \\ " << Y.get_solution_str_value() 
    //  	 << " = 3" << endl;
  }

  //if(s.statistics.num_backtracks != 52) {
  if(s.statistics.num_backtracks != 2958) {
    cout << "Error: wrong number of backtracks! (" 
	 << (s.statistics.num_backtracks) << ")" << endl;
    //exit(1);
  }
  if(num_solutions != 2650) {
    cout << "Error: wrong number of solutions! (" 
	 << (num_solutions) << ")" << endl;
    //exit(1);
  }
}


void SetDifferenceTest::run3() {

  //if(Verbosity) cout << "Run Subset test: "; 

  Solver s;

  Vector<int> lbx;
  Vector<int> ubx;
  Vector<int> lby;
  Vector<int> uby;
  Vector<int> lbz;
  Vector<int> ubz;
  
  ubx.add(-3);
  ubx.add(2);
  ubx.add(5);
  ubx.add(10);
  ubx.add(100);
  ubx.add(123);
  ubx.add(1000);
  ubx.add(1001);

  lbx.add(2);
  lbx.add(100);


  uby.add(-3);
  uby.add(2);
  uby.add(7);
  uby.add(10);
  uby.add(12);
  uby.add(15);
  uby.add(100);
  uby.add(123);
  uby.add(300);
  uby.add(1000);
  uby.add(1001);

  lby.add(100);
  lby.add(300);

  
  lbz.add(2);
  lbz.add(10);

  ubz.add(-3);
  ubz.add(2);
  ubz.add(5);
  ubz.add(7);
  ubz.add(10);
  ubz.add(12);
  ubz.add(15);
  ubz.add(100);
  ubz.add(123);
  ubz.add(300);
  ubz.add(1000);




  Variable X = SetVariable(lbx,ubx,3,6);
  Variable Y = SetVariable(lby,uby,2,5);
  Variable Z = SetVariable(lbz,ubz,2,5);
 

  //std::cout << X << " " << Y << std::endl;

  s.add( SetDifference(X,Y) == Z );

  //std::cout << s << std::endl;

  s.rewrite();
  
  //cout << "Rewrite\n" << s << endl;

  s.consolidate();

  //cout << X.get_domain() << " \\ " << Y.get_domain() << " = " << Z.get_domain() << endl;  
  //cout << "Consolidate\n" << s << endl;  


  s.initialise_search(s.variables,
		      new GenericHeuristic< Lexicographic, MinValue >(&s), 
		      new NoRestart());

  int num_solutions = 0;
  while(s.get_next_solution() == SAT) {
    ++num_solutions;
    // cout << X.get_solution_str_value() << " inter " << Y.get_solution_str_value() << " = " 
    // 	 << Z.get_solution_str_value() << endl;

    // cout 
    // 	 << X.get_solution_str_value() << " \\ " 
    // 	 << Y.get_solution_str_value() << " == "
    // 	 << Z.get_solution_str_value() << endl;

  }

  //if(s.statistics.num_backtracks != 52) {
  if(s.statistics.num_backtracks != 1201) {
    cout << "Error: wrong number of backtracks! (" 
	 << (s.statistics.num_backtracks) << ")" << endl;
    //exit(1);
  }
  if(num_solutions != 1202) {
    cout << "Error: wrong number of solutions! (" 
	 << (num_solutions) << ")" << endl;
    //exit(1);
  }
}




LexTest::LexTest() : UnitTest() {}
LexTest::~LexTest() {}

void LexTest::run() {

  if(Verbosity) cout << "Run Lex test: "; 
  cout.flush();


  Solver s;

  VarArray X(5,1,3);
  VarArray Y(5,1,3);

  s.add( X < Y );

  //std::cout << s << std::endl;

  s.rewrite();
  
  //  cout << "Rewrite\n" << s << endl;

  s.consolidate();
  
  //  cout << "Consolidate\n" << s << endl;  

  s.initialise_search(s.variables,
		      new GenericHeuristic< Lexicographic, MinValue >(&s), 
		      new NoRestart());

  int num_solutions = 0;
  while(s.get_next_solution() == SAT) {
    ++num_solutions;

    // for(unsigned int i=0; i<X.size; ++i)
    //   cout << " " << X[i].get_solution_str_value();
    // cout << " < " << endl;
    // for(unsigned int i=0; i<X.size; ++i)
    //   cout << " " << Y[i].get_solution_str_value();
    // cout << endl << endl;

    // // //cout << s.last_solution_lb << endl;
    // // cout << Z.get_solution_str_value() << " == " 
    // // 	 << X.get_solution_str_value() << " inter " 
    // // 	 << Y.get_solution_str_value() << endl;

  }

  //if(s.statistics.num_backtracks != 29402) {
  if(s.statistics.num_backtracks != 29686) {
    cout << "Error: wrong number of backtracks! (" 
	 << (s.statistics.num_backtracks) << ")" << endl;
    //exit(1);
  }
  if(num_solutions != 29403) {
    cout << "Error: wrong number of solutions! (" 
	 << (num_solutions) << ")" << endl;
    //exit(1);
  }
}



WeightedSumTest::WeightedSumTest() : UnitTest() {}
WeightedSumTest::~WeightedSumTest() {}

void WeightedSumTest::run1() {

  Solver s;

  VarArray X(5,-3,3);
  Vector<int> coefs;
  coefs.add(-1);
  coefs.add(3);
  coefs.add(-2);
  coefs.add(1);
  coefs.add(1);

  s.add( Sum(X, coefs, -2, 2) );

  //std::cout << s << std::endl;

  s.consolidate();

  //cout << "Consolidate\n" << s << endl;  

  s.rewrite();

  //cout << "Rewrite\n" << s << endl;

  s.initialise_search(X,
		      new GenericHeuristic< Lexicographic, MinValue >(&s), 
		      new NoRestart());

  int num_solutions = 0;
  while(s.get_next_solution() == SAT) {
    ++num_solutions;
    int total = coefs[0] * X[0].get_solution_int_value();
    //cout << coefs[0] << "*" << X[0].get_solution_int_value() ;
    for(int i=1; i<5; ++i) {
      total += coefs[i] * X[i].get_solution_int_value();
      //cout << " + " << coefs[i] << "*" << X[i].get_solution_int_value() ;
    }
    //cout << " = " << total << " (in [-1,1])" << endl; 
    //cout << X.get_solution_value() << " subset of " << Y.get_solution_value() << endl;

    if(total < -2 || total > 2) 
      cout << "Error: wrong solution! (" 
	   << (total) << ")" << endl;
  }

  //std::cout << num_solutions << " " << s.statistics.num_backtracks << std::endl;

  if(s.statistics.num_backtracks != 3820) {
    cout << "Error: wrong number of backtracks! (" 
	 << (s.statistics.num_backtracks) << ")" << endl;
    //exit(1);
  }
  if(num_solutions != 3821) {
    cout << "Error: wrong number of solutions! (" 
	 << (num_solutions) << ")" << endl;
    //exit(1);
  }
}


void WeightedSumTest::run2() {



  Solver s;

  VarArray X(5,-2,2);
  Vector<int> coefs;
  coefs.add(-1);
  coefs.add(3);
  coefs.add(-2);
  coefs.add(1);
  coefs.add(1);

  Variable the_sum(-10, 10);
  
  s.add( (the_sum < -6) || (the_sum > 6) );

  s.add( Sum(X, coefs) == the_sum );

  //std::cout << s << std::endl;

  s.rewrite();
  
  //cout << "Rewrite\n" << s << endl;

  s.consolidate();

  //cout << "Consolidate\n" << s << endl;  


  //exit(1);

  s.initialise_search(X,
		      new GenericHeuristic< Lexicographic, MinValue >(&s), 
		      new NoRestart());

  int num_solutions = 0;
  while(s.get_next_solution() == SAT) {
    ++num_solutions;
    int total = coefs[0] * X[0].get_solution_int_value();
    //cout << coefs[0] << "*" << X[0].get_solution_value() ;
    for(int i=1; i<5; ++i) {
      total += coefs[i] * X[i].get_solution_int_value();
      //cout << " + " << coefs[i] << "*" << X[i].get_solution_value() ;
    }
    //cout << " = " << total << " (in [-1,1])" << endl; 
    //cout << X.get_solution_value() << " subset of " << Y.get_solution_value() << endl;

    if(total >= -6 && total <= 6) 
      cout << "Error: wrong solution! (" 
	   << (total) << ")" << endl;
  }

  //std::cout << num_solutions << " " << s.statistics.num_backtracks << std::endl;

  if(s.statistics.num_backtracks != 673) {
    cout << "Error: wrong number of backtracks! (" 
	 << (s.statistics.num_backtracks) << ")" << endl;
    //exit(1);
  }
  if(num_solutions != 674) {
    cout << "Error: wrong number of solutions! (" 
	 << (num_solutions) << ")" << endl;
    //exit(1);
  }
}


void WeightedSumTest::run3() {

  Solver s;

  VarArray X(5,-3,3);
  Vector<int> coefs;
  coefs.add(-1);
  coefs.add(3);
  coefs.add(-2);
  coefs.add(1);
  coefs.add(1);

  s.add( Sum(X, coefs) <= -10 );

  //std::cout << s << std::endl;

  s.rewrite();
  
  //cout << "Rewrite\n" << s << endl;

  s.consolidate();

  //cout << "Consolidate\n" << s << endl;  


  //exit(1);

  s.initialise_search(X,
		      new GenericHeuristic< Lexicographic, MinValue >(&s), 
		      new NoRestart());

  int num_solutions = 0;
  while(s.get_next_solution() == SAT) {
    ++num_solutions;
    int total = coefs[0] * X[0].get_solution_int_value();
    //cout << coefs[0] << "*" << X[0].get_solution_value() ;
    for(int i=1; i<5; ++i) {
      total += coefs[i] * X[i].get_solution_int_value();
      //cout << " + " << coefs[i] << "*" << X[i].get_solution_str_value() ;
    }
    //cout << " = " << total << " (in [-1,1])" << endl; 
    //cout << X.get_solution_str_value() << " subset of " << Y.get_solution_str_value() << endl;

    if(total > -10) 
      cout << "Error: wrong solution! (" 
	   << (total) << ")" << endl;
  }

  //std::cout << num_solutions << " " << s.statistics.num_backtracks << std::endl;

  if(s.statistics.num_backtracks != 2137) {
    cout << "Error: wrong number of backtracks! (" 
	 << (s.statistics.num_backtracks) << ")" << endl;
    //exit(1);
  }
  if(num_solutions != 2138) {
    cout << "Error: wrong number of solutions! (" 
	 << (num_solutions) << ")" << endl;
    //exit(1);
  }
}


void WeightedSumTest::run4() {

  Solver s;

  VarArray X(5,-3,3);
  Vector<int> coefs;
  coefs.add(-1);
  coefs.add(3);
  coefs.add(-2);
  coefs.add(1);
  coefs.add(1);

  s.add( Sum(X, coefs, -INFTY, -10) );

  //std::cout << s << std::endl;

  s.rewrite();
  
  //cout << "Rewrite\n" << s << endl;

  s.consolidate();

  //cout << "Consolidate\n" << s << endl;  


  //exit(1);

  s.initialise_search(X,
		      new GenericHeuristic< Lexicographic, MinValue >(&s), 
		      new NoRestart());

  int num_solutions = 0;
  while(s.get_next_solution() == SAT) {
    ++num_solutions;
    int total = coefs[0] * X[0].get_solution_int_value();
    //cout << coefs[0] << "*" << X[0].get_solution_value() ;
    for(int i=1; i<5; ++i) {
      total += coefs[i] * X[i].get_solution_int_value();
      //cout << " + " << coefs[i] << "*" << X[i].get_solution_str_value() ;
    }
    //cout << " = " << total << " (in [-1,1])" << endl; 
    //cout << X.get_solution_str_value() << " subset of " << Y.get_solution_str_value() << endl;

    if(total > -10) 
      cout << "Error: wrong solution! (" 
	   << (total) << ")" << endl;
  }

  //std::cout << num_solutions << " " << s.statistics.num_backtracks << std::endl;

  if(s.statistics.num_backtracks != 2137) {
    cout << "Error: wrong number of backtracks! (" 
	 << (s.statistics.num_backtracks) << ")" << endl;
    //exit(1);
  }
  if(num_solutions != 2138) {
    cout << "Error: wrong number of solutions! (" 
	 << (num_solutions) << ")" << endl;
    //exit(1);
  }
}

void WeightedSumTest::run() {
  if(Verbosity) cout << "Run WeightedSum test: "; 
  run1();
  cout << "1 ";
  run2();
  cout << "2 ";
  run3();
  cout << "3 ";
  run4();
  cout << "4 ";

  //  cout << " ";
}


ElementTest::ElementTest() : UnitTest() {}
ElementTest::~ElementTest() {}

void ElementTest::run1() {

  Solver s;

  VarArray X;

  Variable x1(-10,10);
  X.add(x1);

  Variable x2(0,100);
  X.add(x2);

  Variable x3(-100,0);
  X.add(x3);

  Vector< int > vals;
  vals.add(-250);
  vals.add(-12);
  vals.add(250);

  Variable x4(vals);
  X.add(x4);

  Variable Y(0,3);
  //X.add(Y);

  vals.add(12);
  Variable Z(vals);

  s.add( X[Y] == Z );

  //std::cout << s << std::endl;

  s.rewrite();
  
  //cout << "Rewrite\n" << s << endl;


  s.consolidate();

  //cout << "Consolidate\n" << s << endl;  

  //std::cout << s << std::endl;

#ifdef _MONITOR
  s.monitor_list << X ;
  s.monitor_list << "[";
  s.monitor_list << Y ;
  s.monitor_list << "] = ";
  s.monitor_list << Z ;
  s.monitor_list << "\n";
#endif


  s.initialise_search(s.variables,
		      new GenericHeuristic< Lexicographic, MinValue >(&s), 
		      new NoRestart());

  int num_solutions = 0;
  while(s.get_next_solution() == SAT) {
    ++num_solutions;

    int index = Y.get_solution_int_value();
    int equal_value = Z.get_solution_int_value();
    int array_value = X[index].get_solution_int_value();
    
    if(array_value != equal_value) {

      cout << "[ ";
      for(int i=0; i<4; ++i) {
	std::cout << X[i].get_var().get_solution_str_value() << " " ;
      }
      cout  << "]: X[" << Y.get_var().get_solution_str_value() << "] = " << Z.get_var().get_solution_str_value() 
	    << " (" << X[Y.get_var().get_solution_int_value()].get_var().get_solution_str_value() << ")" << endl;
      
      exit(1);

    }

//     cout << Y.get_solution_str_value() << endl;
//     cout << Z.get_solution_str_value() << endl;
//     //cout << X[0].get_solution_str_value() << endl;
// cout << X[1].get_solution_str_value() << endl;
// cout << X[2].get_solution_str_value() << endl;
// cout << X[3].get_solution_str_value() << endl;

    //std::cout << num_solutions << std::endl;

    // cout << "[ ";
    // for(int i=0; i<4; ++i) {
    //   std::cout << X[i].get_var().get_solution_str_value() << " " ;
    // }
    // cout  << ": X[" << Y.get_var().get_solution_str_value() << "] = " << Z.get_var().get_solution_str_value() 
    // 	  << " (" << X[Y.get_var().get_solution_int_value()].get_var().get_solution_str_value() << ")" << endl;
  }

  //std::cout << num_solutions << " " << s.statistics.num_backtracks << std::endl;


  //cout << "After\n" << s << endl;  

  if(s.statistics.num_backtracks != 655388) {
    //if(s.statistics.num_backtracks != 428442) {
    cout << "Error: wrong number of backtracks! (" 
	 << (s.statistics.num_backtracks) << ")" << endl;
    //exit(1);
  }
  if(num_solutions != 655389) {
  //if(num_solutions != 428442) {
    cout << "Error: wrong number of solutions! (" 
	 << (num_solutions) << ")" << endl;
    //exit(1);
  }
}


// void ElementTest::run2() {

//   Solver s;

//   VarArray X;

//   // Variable s1 = SetVariable(1,9,3,5);
//   // X.add(s1);

//   // std::cout << s1 << std::endl;


//   // Variable s2 = SetVariable(2,8,2,4);
//   // X.add(s2);

//   // std::cout << s2 << std::endl;


//   // Variable s3 = SetVariable(0,7,2,5);
//   // X.add(s3);

//   // std::cout << s3 << std::endl;


//   // Variable s4 = SetVariable(5,9,1,3);
//   // X.add(s4);

//   // std::cout << s4 << std::endl;


//   // Variable Y(0,3);

//   // std::cout << Y << std::endl;


//   // Variable Z = SetVariable(0,9,1,5);

//   //std::cout << Z << std::endl;
  



//   Variable s1 = SetVariable(1,5,2,4);
//   X.add(s1);

//   //  std::cout << s1 << std::endl;


//   Variable s2 = SetVariable(5,9,1,3);
//   X.add(s2);

//   //  std::cout << s2 << std::endl;


//   Variable s3 = SetVariable(0,3,2,3);
//   X.add(s3);

//   //  std::cout << s3 << std::endl;


//   Variable Y(0,2);

//   //  std::cout << Y << std::endl;


//   Variable Z = SetVariable(0,9,1,4);

//   //std::cout << Z << std::endl;



//   //s.add( X[Y] == Z );

//   s.add( ElementSet(X,Y) == Z );

  
//   //std::cout << s << std::endl;

//   s.rewrite();
  
//   //cout << "Rewrite\n" << s << endl;


//   s.consolidate();

//   //cout << "Consolidate\n" << s << endl;  


//   //exit(1);

//   s.initialise_search(s.variables,
// 		      new GenericHeuristic< Lexicographic, MaxValue >(&s), 
// 		      new NoRestart());

//   int num_solutions = 0;
//   while(s.get_next_solution() == SAT) {
//     ++num_solutions;

//  //    //cout << Y << endl;

// //     cout << "[ ";
// //     for(int i=0; i<3; ++i) {
// //       std::cout << X[i].get_solution_str_value() << " " ;
// //     }

// //     cout << "] : X[" << Y.get_solution_str_value() << "] = "
// // 	 << Z.get_solution_str_value() // << " / "
// // 	 // << X[Y.get_solution_int_value()].get_solution_str_value()
// // 	 << endl;
// // //     //cout << X[0].get_solution_str_value() << endl;
// // // cout << X[1].get_solution_str_value() << endl;
// // // cout << X[2].get_solution_str_value() << endl;
// // // cout << X[3].get_solution_str_value() << endl;

// //     // cout << "X[" << Y.get_solution_str_value() << "] = " << Z.get_solution_str_value() 
// //     //  	 << " (" << X[Y.get_solution_str_value()].get_solution_str_value() << ")" << endl;
//   }

//   //std::cout << num_solutions << " " << s.statistics.num_backtracks << std::endl;


//   //cout << "After\n" << s << endl;  

//   //if(s.statistics.num_backtracks != 18750) {
//   if(s.statistics.num_backtracks != 18749) {
//     cout << "Error: wrong number of backtracks! (" 
// 	 << (s.statistics.num_backtracks) << ")" << endl;
//     //exit(1);
//   }
//   if(num_solutions != 18750) {
//     cout << "Error: wrong number of solutions! (" 
// 	 << (num_solutions) << ")" << endl;
//     //exit(1);
//   }
// }


void ElementTest::run() {
  if(Verbosity) cout << "Run Element test: "; 
  run1();
  cout << "1 ";
  // run2();
  // cout << "2 ";
  // run3();
  // cout << " 3";
  // run4();
  // cout << " 4";

  //cout << " ";
}


MinMaxTest::MinMaxTest() : UnitTest() {}
MinMaxTest::~MinMaxTest() {}

void MinMaxTest::run1() {

  Solver s;

  VarArray X;

  Variable x1(-10,10);
  X.add(x1);

  Variable x2(0,20);
  X.add(x2);

  Variable x3(-10,0);
  X.add(x3);


  Vector< int > vals;
  vals.add(-250);
  vals.add(-12);
  vals.add(-2);
  vals.add(2);
  vals.add(250);
  Variable x4(vals);

  X.add(x4);

  Variable x5(0,3);
  X.add(x5);

  Variable Z(-250,250);


  s.add( Min(X) == Z );

  s.rewrite();
  
  //cout << "Rewrite\n" << s << endl;


  s.consolidate();

  //cout << "Consolidate\n" << s << endl;  


  s.initialise_search(X, //s.variables,
		      new GenericHeuristic< Lexicographic, MinValue >(&s), 
		      new NoRestart());

  int num_solutions = 0;
  while(s.get_next_solution() == SAT) {
    ++num_solutions;

    // cout << "min(" << X[0].get_solution_str_value();
    // for(unsigned int i=1; i<X.size; ++i)
    //   cout << ", " << X[i].get_solution_str_value();
    // cout << ") = " << Z.get_solution_str_value() << endl;

  }

  if(s.statistics.num_backtracks != 97019) {
    cout << "Error: wrong number of backtracks! (" 
	 << (s.statistics.num_backtracks) << ")" << endl;
    //exit(1);
  }
  if(num_solutions != 97020) {
    cout << "Error: wrong number of solutions! (" 
	 << (num_solutions) << ")" << endl;
    //exit(1);
  }
}


void MinMaxTest::run2() {

  Solver s;

  VarArray X;

  Variable x1(-10,10);
  X.add(x1);

  Variable x2(0,20);
  X.add(x2);

  Variable x3(-10,0);
  X.add(x3);

  Vector< int > vals;
  vals.add(-250);
  vals.add(-12);
  vals.add(-2);
  vals.add(2);
  vals.add(250);
  Variable x4(vals);
  X.add(x4);

  Variable x5(0,3);
  X.add(x5);

  Variable Z(-250,250);


  s.add( Max(X) == Z );

  //std::cout << s << std::endl;

  s.rewrite();
  
  //cout << "Rewrite\n" << s << endl;


  s.consolidate();

  //cout << "Consolidate\n" << s << endl;  


  s.initialise_search(X, //s.variables,
		      new GenericHeuristic< Lexicographic, MinValue >(&s), 
		      new NoRestart());

  int num_solutions = 0;
  while(s.get_next_solution() == SAT) {
    ++num_solutions;

    // cout << "max(" << X[0].get_solution_str_value();
    // for(unsigned int i=1; i<X.size; ++i)
    //   cout << ", " << X[i].get_solution_str_value();
    // cout << ") = " << Z.get_solution_str_value() << endl;

  }

  if(s.statistics.num_backtracks != 97019) {
    cout << "Error: wrong number of backtracks! (" 
	 << (s.statistics.num_backtracks) << ")" << endl;
    //exit(1);
  }
  if(num_solutions != 97020) {
    cout << "Error: wrong number of solutions! (" 
	 << (num_solutions) << ")" << endl;
    //exit(1);
  }
}

void MinMaxTest::run() {
  if(Verbosity) cout << "Run Min/Max test: "; 
  run1();
  cout << "1 ";
  run2();
  cout << "2 ";

  //  cout << " ";
}


VarStackDynamicTest::VarStackDynamicTest() : UnitTest() {}

VarStackDynamicTest::~VarStackDynamicTest() {}

void VarStackDynamicTest::run() {

  Queue q;
  q.initialise(7,13);

  q.add(8);

  q.add(7);

  q.add(13);

  q.add(11);

  cout << q << endl;

  cout << q.pop() << endl;

  cout << q << endl;

  cout << q.pop() << endl;

  cout << q << endl;


  q.extend(-3);

  cout << q << endl;
  
  q.add(-3);

  cout << q << endl;


  q.declare(20);

  cout << q << endl;

  cout << q.pop() << endl;

  cout << q << endl;

  q.declare(19);

  cout << q << endl;

  cout << q.pop() << endl;

  cout << q << endl;

  cout << q.pop() << endl;

  cout << q << endl;

  cout << q.pop() << endl;

  cout << q << endl;

 cout << q.pop() << endl;

  cout << q << endl;


  q.add(-1);
  cout << q << endl;

  q.add(5);
  cout << q << endl;

  q.add(15);
  cout << q << endl;


//   IntStack witness;
//   witness.initialise(0,10,true);

//   cout << witness << endl;

//   for(int i=0; i<10; ++i) {
//     int elt = 10+(i+1)*5;
//     cout << "add " << elt << ": ";
//     witness.declare(elt);
//     cout << witness << endl;
//   }

//   for(int i=0; i<10; ++i) {
//     int elt = (i+1)*-5;
//     cout << "add " << elt << ": ";
//     witness.declare(elt);
//     cout << witness << endl;
//   }
  
}


OpshopTest::OpshopTest() : UnitTest() {}
OpshopTest::~OpshopTest() {}


void OpshopTest::run1() {

  // read file
  int nJobs, nMachines;
  int i=0, j, k, dur=0, lb, opt, bufsize=1000;
  char buf[bufsize];
  std::ifstream infile( "./examples/data/GP03-01.txt", std::ios_base::in );
  //std::ifstream infile( "./examples/data/tiny.txt", std::ios_base::in );
	
  do {
    infile.getline( buf, bufsize, '\n' );
  } while( buf[0] == '#' );
	
  while( buf[i] != ' ' ) ++i;
  buf[i] = '\0';
  lb = atoi( buf );
	
  while( buf[i] == '\0' || buf[i] == ' ' || buf[i] == '*' ) ++i;
  j = i;
  while( buf[i] != ' ' && buf[i] != '\n' && buf[i] != '\0' ) ++i;
  buf[i] = '\0';
  opt = atoi( &(buf[j]) );
	
  do {
    infile.get( buf[0] );
    if( buf[0] != '#' ) break;
    infile.getline( buf, bufsize, '\n' );
  } while( true );
  infile.unget();
	
  infile >> nJobs;
  infile >> nMachines;

  int *jobs[nJobs];
  for(i=0; i<nJobs; ++i) {
    jobs[i] = new int[nMachines];
    std::fill(jobs[i], jobs[i]+nMachines, 0);
  }

  infile.getline( buf, bufsize, '\n' );
	
  do {
    infile.get( buf[0] );
    if( buf[0] != '#' ) break;
    infile.getline( buf, bufsize, '\n' );
  } while( true );
  infile.unget();
	
  k = 0;
  for(i=0; i<nJobs; ++i) {
    for(j=0; j<nMachines; ++j) {
      infile >> k;
      jobs[i][j] = k;
      dur += jobs[i][j];
    }
  }

  VarArray task;

  int makespan = (int)((double)opt*1.2);

  for(i=0; i<nJobs; ++i) {
    for(j=0; j<nMachines; ++j) {
      Variable t(0, makespan-jobs[i][j]);
      task.add(t);
    }
  }

  Solver s;

  s.add(task);
  
  for(int i=0; i<nJobs; ++i) 
    {
      for(int j=0; j<nMachines-1; ++j)
	for(int k=j+1; k<nMachines; ++k) {

	  //s.add(Precedence(task[i*nMachines+j], jobs[i][j], task[i*nMachines+k]));

	  s.add(Disjunctive(task[i*nMachines+j], 
	  		    task[i*nMachines+k], 
	  		    jobs[i][j],
	  		    jobs[i][k]
	  		    ));
	}
    }

  for(int j=0; j<nMachines; ++j)
    for(int i=0; i<nJobs-1; ++i) 
      {
	for(int k=i+1; k<nJobs; ++k) {
	  s.add(Disjunctive(task[i*nMachines+j], 
			    task[k*nMachines+j], 
			    jobs[i][j],
			    jobs[k][j]
			    ));
	}
    }

  s.consolidate();
  makespan = (int)((double)opt*1.0001);
  for(int i=0; i<nJobs; ++i)
    for(int j=0; j<nMachines; ++j) 
      s.add( task[i*nMachines+j] <= makespan-jobs[i][j] );


  s.consolidate();

  //std::cout << s << std::endl;
  //s.specialise();
  //std::cout << s << std::endl;

#ifdef _MONITOR
  //s.monitor_list << task;

  Variable x;

  s.monitor_list << " t0 ";
  x = task[0];
  s.monitor_list << x;

  s.monitor_list << " t1 ";
  x = task[1];
  s.monitor_list << x;

  s.monitor_list << " t2 ";
  x = task[2];
  s.monitor_list << x;

  s.monitor_list << " t3 ";
  x = task[3];
  s.monitor_list << x;

  s.monitor_list << " t4 ";
  x = task[4];
  s.monitor_list << x;

  s.monitor_list << " t5 ";
  x = task[5];
  s.monitor_list << x;


#endif

  s.initialise_search(task, 
		      new GenericHeuristic< GenericDVO < MinDomain >, HalfSplit >(&s), 
		      new NoRestart());

  int num_solutions = 0;
  while(s.get_next_solution() == SAT) {
    ++num_solutions;
    // for(int i=0; i<task.size; ++i) {
    //   std::cout << task[i].get_solution_int_value() << " ";
    // }
    // std::cout << std::endl;
  }


  if(s.statistics.num_backtracks != 774518) {
    cout << "Error: wrong number of backtracks! (" 
	 << (s.statistics.num_backtracks) << ")" << endl;
    //exit(1);
  }
  if(num_solutions != 771316) {
    cout << "Error: wrong number of solutions! (" 
	 << (num_solutions) << ")" << endl;
    //exit(1);
  }

}


void OpshopTest::run2() {

  // read file
  int nJobs, nMachines;
  int i=0, j, k, dur=0, lb, opt, bufsize=1000;
  char buf[bufsize];
  std::ifstream infile( "./examples/data/GP03-01.txt", std::ios_base::in );
	
  do {
    infile.getline( buf, bufsize, '\n' );
  } while( buf[0] == '#' );
	
  while( buf[i] != ' ' ) ++i;
  buf[i] = '\0';
  lb = atoi( buf );
	
  while( buf[i] == '\0' || buf[i] == ' ' || buf[i] == '*' ) ++i;
  j = i;
  while( buf[i] != ' ' && buf[i] != '\n' && buf[i] != '\0' ) ++i;
  buf[i] = '\0';
  opt = atoi( &(buf[j]) );
	
  do {
    infile.get( buf[0] );
    if( buf[0] != '#' ) break;
    infile.getline( buf, bufsize, '\n' );
  } while( true );
  infile.unget();
	
  infile >> nJobs;
  infile >> nMachines;

  int *jobs[nJobs];
  for(i=0; i<nMachines; ++i)
    jobs[i] = new int[nMachines];

  infile.getline( buf, bufsize, '\n' );
	
  do {
    infile.get( buf[0] );
    if( buf[0] != '#' ) break;
    infile.getline( buf, bufsize, '\n' );
  } while( true );
  infile.unget();
	
  k = 0;
  for(i=0; i<nJobs; ++i) {
    for(j=0; j<nMachines; ++j) {
      infile >> jobs[i][j];
      dur += jobs[i][j];
    }
  }

  VarArray task;
  VarArray disjuncts;

  int makespan = (int)((double)opt*1.2);

  for(i=0; i<nJobs; ++i) {
    for(j=0; j<nMachines; ++j) {
      Variable t(0, makespan-jobs[i][j]);
      task.add(t);
    }
  }

  Solver s;

  s.add(task);
  
  for(int i=0; i<nJobs; ++i) 
    {
      for(int j=0; j<nMachines-1; ++j)
	for(int k=j+1; k<nMachines; ++k) {
	  disjuncts.add(ReifiedDisjunctive(task[i*nMachines+j], 
					   task[i*nMachines+k], 
					   jobs[i][j],
					   jobs[i][k]
					   ));
	}
    }

  for(int j=0; j<nMachines; ++j)
    for(int i=0; i<nJobs-1; ++i) 
      {
	for(int k=i+1; k<nJobs; ++k) {
	  disjuncts.add(ReifiedDisjunctive(task[i*nMachines+j], 
					   task[k*nMachines+j], 
					   jobs[i][j],
					   jobs[k][j]
					   ));
	}
      }

  for(unsigned int i=0; i<disjuncts.size; ++i)
    s.add(Free(disjuncts[i]));


  s.consolidate();
  makespan = (int)((double)opt*1.0001);
  for(int i=0; i<nJobs; ++i)
    for(int j=0; j<nMachines; ++j) 
      s.add( task[i*nMachines+j] <= makespan-jobs[i][j] );

  //std::cout << s << std::endl;

  //exit(1);

#ifdef _DEBUG_PRUNING
  s.monitor(task);
#endif

  s.initialise_search(disjuncts, 
		      new GenericHeuristic< GenericDVO < MinDomain >, HalfSplit >(&s), 
		      new NoRestart());

  int num_solutions = 0;
  while(s.get_next_solution() == SAT) {
    ++num_solutions;
  }


  if(s.statistics.num_backtracks != 9) {
    cout << "Error: wrong number of backtracks! (" 
	 << (s.statistics.num_backtracks) << ")" << endl;
    //exit(1);
  }
  if(num_solutions != 4) {
    cout << "Error: wrong number of solutions! (" 
	 << (num_solutions) << ")" << endl;
    //exit(1);
  }

}


void OpshopTest::run3() {

  // read file
  int nJobs, nMachines;
  int i=0, j, k, dur=0, lb, opt, bufsize=1000;
  char buf[bufsize];
  std::ifstream infile( "./examples/data/GP03-01.txt", std::ios_base::in );
	
  do {
    infile.getline( buf, bufsize, '\n' );
  } while( buf[0] == '#' );
	
  while( buf[i] != ' ' ) ++i;
  buf[i] = '\0';
  lb = atoi( buf );
	
  while( buf[i] == '\0' || buf[i] == ' ' || buf[i] == '*' ) ++i;
  j = i;
  while( buf[i] != ' ' && buf[i] != '\n' && buf[i] != '\0' ) ++i;
  buf[i] = '\0';
  opt = atoi( &(buf[j]) );
	
  do {
    infile.get( buf[0] );
    if( buf[0] != '#' ) break;
    infile.getline( buf, bufsize, '\n' );
  } while( true );
  infile.unget();
	
  infile >> nJobs;
  infile >> nMachines;

  int *jobs[nJobs];
  for(i=0; i<nMachines; ++i)
    jobs[i] = new int[nMachines];

  infile.getline( buf, bufsize, '\n' );
	
  do {
    infile.get( buf[0] );
    if( buf[0] != '#' ) break;
    infile.getline( buf, bufsize, '\n' );
  } while( true );
  infile.unget();
	
  k = 0;
  for(i=0; i<nJobs; ++i) {
    for(j=0; j<nMachines; ++j) {
      infile >> jobs[i][j];
      dur += jobs[i][j];
    }
  }

  VarArray task;
  VarArray disjuncts;

  int makespan = (int)((double)opt*1.2);

  for(i=0; i<nJobs; ++i) {
    for(j=0; j<nMachines; ++j) {
      Variable t(0, makespan-jobs[i][j]);
      task.add(t);
    }
  }

  Solver s;

  s.add(task);
  
  for(int i=0; i<nJobs; ++i) 
    {
      for(int j=0; j<nMachines-1; ++j)
	for(int k=j+1; k<nMachines; ++k) {
	  disjuncts.add(ReifiedDisjunctive(task[i*nMachines+j], 
					   task[i*nMachines+k], 
					   jobs[i][j],
					   jobs[i][k]
					   ));
	}
    }

  for(int j=0; j<nMachines; ++j)
    for(int i=0; i<nJobs-1; ++i) 
      {
	for(int k=i+1; k<nJobs; ++k) {
	  disjuncts.add(ReifiedDisjunctive(task[i*nMachines+j], 
					   task[k*nMachines+j], 
					   jobs[i][j],
					   jobs[k][j]
					   ));
	}
      }

  for(unsigned int i=0; i<disjuncts.size; ++i)
    s.add(Free(disjuncts[i]));


  s.consolidate();
  makespan = (int)((double)opt*1.0001);
  for(int i=0; i<nJobs; ++i)
    for(int j=0; j<nMachines; ++j) 
      s.add( task[i*nMachines+j] <= makespan-jobs[i][j] );

  //std::cout << s << std::endl;

  //exit(1);

#ifdef _DEBUG_PRUNING
  s.monitor(task);
#endif

  s.initialise_search(task, 
		      new GenericHeuristic< GenericDVO < MinDomain >, HalfSplit >(&s), 
		      new NoRestart());

  //int last_sol[disjuncts.size], new_sol[disjuncts.size];

  int num_solutions = 0;
  while(s.get_next_solution() == SAT) {
    

    // bool diff = false;
    // for(unsigned int i=0; i<disjuncts.size; ++i) {
    //   new_sol[i] = disjuncts[i].get_solution_int_value();
    //   if(!num_solutions || new_sol[i] != last_sol[i]) diff = true;
    // }


    // if(diff) {
    //   for(unsigned int i=0; i<disjuncts.size; ++i) {
    // 	last_sol[i] = new_sol[i] ;
    // 	std::cout << last_sol[i] ;
    //   }
    //   std::cout << std::endl;
    // }

    ++num_solutions;
  }


  if(s.statistics.num_backtracks != 774518) { //775654) {
    cout << "Error: wrong number of backtracks! (" 
	 << (s.statistics.num_backtracks) << ")" << endl;
    //exit(1);
  }
  if(num_solutions != 771316) {
    cout << "Error: wrong number of solutions! (" 
	 << (num_solutions) << ")" << endl;
    //exit(1);
  }

}


void OpshopTest::run() {
  if(Verbosity) cout << "Run open-shop (GP03-01): "; 

  //double TIME = get_run_time();
  cout << "1 ";
  cout.flush();
  run1();
  //cout << (get_run_time() - TIME) << endl ;
  
  // //TIME = get_run_time();
  // cout << "2 ";
  // cout.flush();
  // run2();
  // //cout << (get_run_time() - TIME) << endl ;
    
  // //TIME = get_run_time();
  // cout << "3 ";
  // cout.flush();
  // run3();
  // //cout << (get_run_time() - TIME) << endl ;

}


RewriteTest1::RewriteTest1() : UnitTest() {}
RewriteTest1::~RewriteTest1() {}

void RewriteTest1::run() {

  if(Verbosity) cout << "Run rewritting test ";

  Variable X(0,10);
  Variable Y(-5,5);
  Variable Z(-100,10);
  Variable A(0, 10);

  int N=10;
  VarArray E(N, -50, 50);

  Vector< Variable > all_vars;
  for(int i=0; i<N; ++i)
    all_vars.add(E[i]);
  all_vars.add(X);
  all_vars.add(Y);
  all_vars.add(Z);

  Solver s;
  s.add(X);
  s.add(Y);
  s.add(Z);

  s.add(X != Z);
  s.add(Y <= Z);

  s.add(X == Y);
  for(int i=0; i<N; ++i) s.add(X == E[i]);


  s.add(E[3] <= A);

  s.add(E[7] != A);

  s.add((E[8] > A) != (A == 3));

  //std::cout << s << std::endl;

  s.rewrite();

  // int sum_degree = X.get_degree() + Y.get_degree() + Z.get_degree();
  // for(int i=0; i<N; ++i)
  //   sum_degree += E[i].get_degree();

  // std::cout << sum_degree << std::endl;
  

  // if(sum_degree != 26) {
  //   cout << "Error while rewritting (wrong number of constraints)!" << endl
  // 	 << "degree(X) = " << X.get_degree() << endl
  // 	 << "degree(Y) = " << Y.get_degree() << endl
  // 	 << "degree(Z) = " << Z.get_degree() << endl
  // 	 << "degree(E[5]) = " << E[5].get_degree() << endl;
    
  //   exit(1);
  // }
  
  // //std::cout << s << std::endl;

  // //s.consolidate();

  //std::cout << s << std::endl;

 s.add(E[4] > 1);
 s.propagate();
 s.consolidate();

 //std::cout << "add a constraint on x[0]\n" << s << std::endl;

 s.initialise_search(all_vars, 
		     new GenericHeuristic< Lexicographic, MinValue >(&s), 
		     new NoRestart());
 
 int num_solutions = 0;
 while(s.get_next_solution() == SAT) {
   Solution sol(all_vars);

   //std::cout << sol << std::endl;

   ++num_solutions;
   if(num_solutions > 50) break;
 }
 
 cout << "(" << num_solutions << ") " ;

 if(num_solutions != 8) {
   cout << "Error wrong number of solution! (" << num_solutions << ")" << endl;
   //exit(1);
 }

}



/******************************************
 * URCSPGenerator
 ******************************************/

URCSPGenerator::URCSPGenerator(int S, int v, int d, int c, int a, int n) 
{
  Seed = S;
  usrand(Seed);
  initURCSP( v, d, c, a, n );  
}

void URCSPGenerator::initURCSP( int v, int d, int c, int a, int n ) {
  Var  = v;
  Dom  = d;
  Con  = c;
  Ari  = a;
  Ngd  = n;
  
  if (Var < 2) 
    std::cerr<<" ***Illegal number of variables: "<<Var<<std::endl;
  if (Dom < 2) 
    std::cerr<<" ***Illegal domain size:"<<Dom<<std::endl;
  if (Con < 0 || Con > comb(Var,Ari)) 
    std::cerr<<" ***Illegal number of constraints: "<<Con << std::endl;
  if (Ari < 2)
    std::cerr<<" ***Illegal arity: "<<Ari<<std::endl;
  if (Ngd < 0 || Ngd > ((int)(pow((double)Dom,Ari) - 1)))
    std::cerr<<" ***Illegal tightness: "<<Ngd<<std::endl;

  int i;
  PossibleCTs = comb(Var,Ari);
  CTarray = new int[PossibleCTs];
  PossibleNGs = 1;
  if( Ngd ) {
    for( i=0; i<Ari; ++i )
      PossibleNGs = PossibleNGs*Dom;
  }
  NGarray = new int[PossibleNGs];
  vars = new int[Ari+1];
  //vals = new int[Var];
  vals = new int[Ari];
  vars[0]=-1;
  uv=Var;
  unconnected_var = new int[Var];
  reInit();
}

URCSPGenerator::~URCSPGenerator()
{
  delete [] CTarray;
  delete [] NGarray;
  delete [] vars;
  delete [] vals;
  delete [] unconnected_var;
}

int URCSPGenerator::comb(int b, int d){
  int res = 1;
  int x=1;
  for(int i=0;i<d;i++){
    res = res*b--;
    x = x*(i+1);
  }
  return res/x;
}

int URCSPGenerator::i2comb(int& code, int ind, int base, int dim, int stop){ 
  int c=comb(base,dim); 
  if(ind == stop || code < c)       
    return ind;        
  else{     
    code=(code-c);        
    return i2comb(code, ind+1, base-1, dim,stop);
  }
}

void URCSPGenerator::reInit() 
{
  int i;
  for( i=0; i<PossibleCTs; ++i )
    CTarray[i]=i;
  for( i=0; i<Var; ++i )
    unconnected_var[i]=1;
  c = 0;
}

void URCSPGenerator::reInit( int v, int d, int c, int a, int n ) 
{
  delete [] CTarray;
  delete [] NGarray;
  delete [] vars;
  delete [] vals;
  delete [] unconnected_var;
  initURCSP( v, d, c, a, n );
}

bool URCSPGenerator::erateConstraint() 
{
  t = 0;
  if( c == Con ) return false;
  /* Choose a random_ number between c and PossibleCTs - 1, inclusive. */
  int qq = randint(PossibleCTs - c);
  int i, r =  c + qq;
  /* Swap elements [c] and [r]. */
  selectedCT = CTarray[r];
  CTarray[r] = CTarray[c];
  CTarray[c] = selectedCT;
  /* Create the constraint. */
  for(i=0;i<Ari;i++) {
    vars[i+1] = i2comb(CTarray[c],vars[i]+1,Var-(vars[i]+1)-1,Ari-1-i,Var-Ari+i);
    if(unconnected_var[vars[i+1]]) {
      unconnected_var[vars[i+1]]=0;	
      --uv;
    }
  }
  /* Initialize the NGarray. */
  for (i=0; i<PossibleNGs; ++i)
    NGarray[i] = i;
  ++c;
  return true;
}


bool URCSPGenerator::erateNogood() {
  if( t == Ngd ) return false;
  int r =  t + randint(PossibleNGs - t);
  selectedNG = NGarray[r];
  NGarray[r] = NGarray[t];
  NGarray[t] = selectedNG;
  /* Add the nogood to the constraint. */
  for(int i=0;i<Ari;i++) {
    int x = (int)(pow((double)Dom,(Ari-i-1)));	    
    vals[i] = NGarray[t] / x;
    NGarray[t] -= (vals[i]*x);
  }
  ++t;
  return true;
}  
