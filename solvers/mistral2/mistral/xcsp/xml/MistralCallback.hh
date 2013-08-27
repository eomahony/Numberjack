#ifndef _MistralCallback_h_
#define _MistralCallback_h_

#include "CSPParserCallback.hh"
#include <mistral_solver.hpp>
#include <mistral_search.hpp>
//#include <mistral_sat.h>

#define ShowNumericalIDS


using namespace Mistral;
using namespace std;

namespace CSPXMLParser {

  // definition of the functions which are called by the parser when it
  // reads new information in the file. These functions are in charge of
  // creating the data structure that the solver needs to do his job
  class MistralCallback : public CSPParserCallback
  {
  public:
    
    // parser parameters
    int verbose;
    
    // timing
    double time_start;
    double time_now;
    double time_last;
    double time_solving;

    // solver
    Solver solver;
    // variables, as stated in the xml file
    VarArray X;

    // constraints

    // network params
    int numVariables;
    int numExtraVars;
    int numDomains;
    int numRelations;
    int numConstraints;
    int numPredicates;
    int numGlobal;
    int numWeightedSum;
    int numAllDiff;
    int numElement;
    int numCumulative;
    int numGlobalCardinality;
    int numEXTCons;
    int numPREDCons;

    int numBinExt;
    int numNaryExt;
    int numLargeExt;

    // domains
    Vector<int> *domains;
    //int *domainSize;
    int *domainMax;
    int *domainMin;
    int *D;    

    // constraints
    string constraintReference;
    VarArray constraintScope;
    VarArray predicateScope;
    //int *params;
    AST** predicates;
    int consArity;
    //int *predicateCount;
    //int *predicateScore;
    class Relation {
   
    public:
      
      RelType rel_type;
      int arity;
      Vector< int > tuples;

      Relation() {}
      // Relation(const int a, const int n) {
      // 	initialise(a, n, r);
      // }
      void initialise(const int a, const int n, const RelType r) {
	rel_type = r;
	arity = a;
	tuples.initialise(a*n);
      }

      int get_table_size() {
	return tuples.size/arity;
      }

      int* get_tuple(const int i) {
	return &(tuples.stack_[i*arity]);
      }

      void add(int* T) {
	for(int i=0; i<arity; ++i) {
	  tuples.add(T[i]);
	}
      }
    };
    
    Relation *relations;

    // utilities
    int countI, countJ, countK, maxCount, curCount, idcur;
    CSPDefinitionType curtype;

  public: 

    MistralCallback() {
      numExtraVars = numVariables = numDomains = numGlobal
	= numRelations = numConstraints = numPredicates 
	= numWeightedSum = numAllDiff = numElement 
	= numCumulative = numGlobalCardinality
	= numBinExt = numNaryExt = numLargeExt 
	= numEXTCons = numPREDCons = 0;
      verbose = 1;
    }

    virtual ~MistralCallback() 
      {
	cleanup();
      }

    virtual void beginInstance(const string & name)
    {
      if( verbose )
	cout << "c \nc start parsing xml file (" << name << ")" << endl;
      time_start = get_run_time();
    }

    virtual void beginDomainsSection(int nbDomains)
    {
      if( verbose ) {
	time_last = get_run_time();
	cout << "c" << setw(10) << " " << "domains....." ;
	cout.flush();
      }

      numDomains = nbDomains;
      //domainSize = new int[numDomains];
      //domains = new int*[numDomains];
      domains = new Vector<int>[numDomains];
      domainMin = new int[numDomains];
      domainMax = new int[numDomains];

      countI = 0;
      curCount = 0;
      maxCount = nbDomains;
    }
  
    virtual void beginDomain(const string & name,
			     int idDomain,
			     int nbValue)
    {
      //domainSize[countI] = nbValue;
      //domains[countI] = new int[nbValue];
      domainMax[countI] = -INFTY;
      domainMin[countI] =  INFTY;
      domains[countI].initialise(nbValue);
      countJ = 0;
    }

    void addDomainValue(int v)
    {
      if( v > domainMax[countI] )
	domainMax[countI] = v;
      if( v < domainMin[countI] )
	domainMin[countI] = v;
      domains[countI].add(v);
    }

    virtual void addDomainValue(int first,int last)
    {
      if( last  > domainMax[countI] )
	domainMax[countI] = last;
      if( first < domainMin[countI] )
	domainMin[countI] = first;
      for(int k=first; k<=last; ++k)
	domains[countI].add(k);
    }

    virtual void endDomain()
    {
      ++countI;

      while(curCount < ((int)(10 * countI / maxCount))) {
	++curCount;
	if( verbose ) {
	  cout << ".";
	  cout.flush();
	}	
      }
    }

    /*!
     * end the definition of all domains
     */
    virtual void endDomainsSection()
    {
      assert( countI == numDomains );

      if( verbose ) {
	time_now = get_run_time();
	cout << setw(6) << setprecision(5)
	     << (time_now - time_last) << "\t(" << numDomains 
	     << " domains)" << endl;
	if( verbose > 2 )
	  for(int i=0; i<numDomains; ++i)
	    {
	      cout << "D" << i << ":" ;
	      for(unsigned int j=0; j<domains[i].size; ++j)
		cout << " " << domains[i][j];
	      cout << endl;
	    }
	time_last = time_now;
      }
    }


    virtual void beginVariablesSection(int nbVariables)
    {
      if( verbose ) {
	cout << "c" << setw(10) << " " << "variables..." ;
	cout.flush();
      }

      numVariables = nbVariables;
      X.initialise( numVariables );
      //D = new int[numVariables];

      countI = 0;
      curCount = 0;
      maxCount = nbVariables;
    }
  
    virtual void addVariable(const string & name, int idVar,
			     const string & domain, int idDomain)
    {

      int dtype = 
	LIST_VAR;
	//BITSET_VAR;

      ++countI;
      //D[idVar] = idDomain;
      //solution[idVar] = domainMin[idDomain];
      if((int)(domains[idDomain].size) == (domainMax[idDomain] - domainMin[idDomain] + 1)) {
	X.add( Variable(domainMin[idDomain], domainMax[idDomain], dtype) );
      } else {
	X.add( Variable(domainMin[idDomain], domainMax[idDomain], domains[idDomain], dtype) );
      }

      //std::cout << X.back() << " in " << X.back().get_domain() << std::endl;


      while(curCount < ((int)(10 * countI / maxCount))) {
	++curCount;
	if( verbose ) {
	  cout << ".";
	  cout.flush();
	}
      }
    }
      

    virtual void endVariablesSection()
    {
      assert( countI == numVariables );
      
      solver.add( X );

      if( verbose ) {
	time_now = get_run_time();
	cout << setw(6) << setprecision(5)
	     << (time_now - time_last) << "\t(" << numVariables 
	     << " variables)" << endl;
	time_last = time_now;
      }
    }


    virtual void beginRelationsSection(int nbRelations)
    {
      if( verbose ) {
	cout << "c" << setw(10) << " " << "relations..." ;
	cout.flush();
      }

      numRelations = nbRelations;
      relations = new Relation[numRelations];
      // relationsIndex = new int[numRelations];
      // relationSpin = new int[numRelations];

      countI = 0;
      curCount = 0;
      maxCount = nbRelations;
    }
      
  
    virtual void beginRelation(const string & name, int idRel,
			       int arity, int nbTuples, RelType relType) 
    {
      //relationsIndex[countI] = CSP::addTable( arity, nbTuples );
      //relationSpin[countI] = (relType == REL_SUPPORT);
      relations[countI].initialise(arity, nbTuples, relType);
      countJ=0;
    }
    
    virtual void addRelationTuple(int arity, int tuple[]) 
    {
      //CSP::addTuple( relationsIndex[countI], tuple );
      relations[countI].add( tuple );
      ++countJ;
    }

    virtual void addRelationTuple(int arity, int tuple[], int cost) 
    {      
    }

    virtual void endRelation()
    {
      ++countI;

      while(curCount < ((int)(10 * countI / maxCount))) {
	++curCount;
	if( verbose ) {
	  cout << ".";
	  cout.flush();
	}
      }
    }

    virtual void endRelationsSection()
    {
      assert( countI == numRelations );
      
      if( verbose ) {
	time_now = get_run_time();
	cout << setw(6) << setprecision(5)
	     << (time_now - time_last) << "\t(" << numRelations 
	     << " relations)" << endl;
	time_last = time_now;
      }
    }

    virtual void beginPredicatesSection(int nbPredicates)
    {
      if( verbose ) {
	cout << "c" << setw(10) << " " << "predicates.." ;
	cout.flush();
      }
      
      numPredicates = nbPredicates;
      predicates = new AST*[numPredicates];
      countI = 0;
      curCount = 0;
      maxCount = nbPredicates;
    }
  
    virtual void beginPredicate(const string & name, int idPred)
    {
      //cout << name << " " << idPred << endl;
    }

    virtual void addFormalParameter(int pos, const string & name,
				    const string & type)
    {
    }

    virtual void predicateExpression(AST *tree)
    {
      predicates[countI] = tree;
    }

    virtual void predicateExpression(const string &expr)
    {
      cerr << "cannot handle expressions in this form" << endl;
      exit(0);
    }

    virtual void endPredicate()
    {
      ++countI;

      while(curCount < ((int)(10 * countI / maxCount))) {
	++curCount;
	if( verbose ) {
	  cout << ".";
	  cout.flush();
	}
      }
    }

    virtual void endPredicatesSection()
    {
      assert( countI == numPredicates );

      if( verbose ) {
	time_now = get_run_time();
	cout << setw(6) << setprecision(5)
	     << (time_now - time_last) << "\t(" << numPredicates 
	     << " predicates)" << endl;
	time_last = time_now;
      }
    }

    virtual void beginConstraintsSection(int nbConstraints)
    {
      if( verbose ) {
	cout << "c" << setw(10) << " " << "constraints." ;
	cout.flush();
      }

      numConstraints = nbConstraints;
      countI = 0;
      curCount = 0;
      maxCount = nbConstraints;
    }
  
    virtual void beginConstraint(const string & name, int idConstr,
				 int arity, 
				 const string & reference, 
				 CSPDefinitionType type, int id,
				 const ASTList &scope)
    {
      //if( arity > maxArity ) maxArity = arity;

      idcur = id;
      curtype = type;
      consArity = arity;
      constraintScope.clear();
      predicateScope.clear();
      for(int i=0; i<arity; ++i) 
       	constraintScope.add( X[scope[i].getVarId()] );
      constraintReference=reference;

    }

    virtual void constraintParameters(const ASTList &args) 
    {
      if (constraintReference=="global:cumulative")
	{
	}
      else if (constraintReference=="global:element")
	{
	}
      else if (constraintReference=="global:weightedsum")
	{
	}
      else if (constraintReference=="global:global_cardinality")
	{
	}
      else 
	{
	}
    }


    virtual void endConstraint()
    {
      
      if( curtype == RelationType ) {
	++numEXTCons;
	if( consArity == 2 )
	  ++numBinExt;
	else if( consArity < 5 )
	  ++numNaryExt;
	else
	  ++numLargeExt;
	
	//ConstraintTable *tab = new ConstraintTable(constraintScope);
	//ConstraintTable *tab = new ConstraintGAC3(constraintScope);
	ConstraintTable *tab = new ConstraintGAC2001(constraintScope);

	for(int i=0; i<relations[idcur].get_table_size(); ++i) {
	  tab->add( relations[idcur].get_tuple(i) );
	}

	solver.add( Constraint(tab) );

	// Table con( constraintScope, relationsIndex[idcur], relationSpin[idcur] ;
	// model.add( con );  
      } else if( curtype == PredicateType ) { 
	++numPREDCons;
      } else if( curtype == GlobalConstraintType ) {
	++numGlobal;
	if (constraintReference=="global:alldifferent") {
	  ++numAllDiff;
	  //model.add( AllDifferent(constraintScope) ); 
	} else if (constraintReference=="global:element") {
	  ++numElement;
	  //model.add( Element(predicateScope) );
	} else if (constraintReference=="global:weightedsum") {
	  ++numWeightedSum;
	} else if (constraintReference=="global:cumulative" ) {
	  ++numCumulative;
	} else if (constraintReference=="global:global_cardinality" ) {
	  ++numGlobalCardinality;
	} else {
	  cout << endl << "s NOT SUPPORTED" << endl;
	  exit(0);
	}
      }

      ++countI;
      
      while(curCount < ((int)(10 * countI / maxCount))) {
	++curCount;
	if( verbose ) {
	  cout << ".";
	  cout.flush();
	}
      }
    }

    /**
     * end the definition of all constraints
     */
    virtual void endConstraintsSection()
    {
      // if( numRelations ) {
      // 	delete [] relationSpin;
      // 	delete [] relationsIndex;
      // }

      // for(int k=0; k<numDomains; ++k)
      // 	delete [] domains[k];
      delete [] domains;
      // if( numPredicates ) {
      // 	delete [] predicateCount;
      // 	delete [] predicateScore;
      // }
      if( verbose ) {
	time_now = get_run_time();
	cout << setw(6) << setprecision(5)
	     << (time_now - time_last) << "\t(" << numConstraints 
	     << " constraints)" << endl;
	time_last = time_now;
      }
    }

    /********************************************************************/


    /*!
     * signal the end of parsing
     */
    virtual void endInstance()
    {
      if( verbose > 1 ) {
	cout << "c\n"
	     << "c" << setw(10) << numEXTCons << " extensional constraints" << endl
	     << "c" << setw(10) << numPREDCons << " predicates " << endl;
      }
      if( verbose ) {
	time_last = get_run_time();
	time_now = (time_last - time_start);
	cout << "c Parsing time: " << setw(6) 
	     << setprecision(5) << time_now 
	     << "\nc\nc Building model" << endl;
      }
    }

    Outcome solve() {
      
      std::cout << solver << std::endl;
      
      //BranchingHeuristic *heuristic = new GenericHeuristic< GenericDVO< MinDomainOverWeight, 1, FailureCountManager >, MinValue >(&solver);
      BranchingHeuristic *heuristic = new GenericHeuristic< GenericDVO< MinDomainOverDegree >, MinValue >(&solver);

      RestartPolicy *restart = new Geometric();
      
      Outcome result = solver.depth_first_search(X, heuristic, restart);

      Solution solution(X);

      cout << solution << std::endl;
      

      return result;
    }


    void print_outcome(int result) {
    }

    void cleanup() {
      if( numPredicates ) {
	for(int i=0; i<numPredicates; ++i)
	  delete predicates[i];
	delete [] predicates;
      }
      delete [] relations;
      // delete [] solution;
      // delete [] params;
      //delete [] D;
      //delete [] domainSize;
      delete [] domainMax;
      delete [] domainMin;
      // delete [] features_vec;
      // int i=trash.size;
      // while( i-- )
      // 	delete trash[i];
      // delete cp_solver;
      // delete sat_solver;
    }


  };

}; // namespace

#endif

