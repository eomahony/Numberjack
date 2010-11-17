#ifndef _MistralCallback_h_
#define _MistralCallback_h_

#include "CSPParserCallback.hh"
#include <mistral_sol.h>
#include <mistral_sat.h>

#define ShowNumericalIDS


Vector<BuildObjectConstraint*> trash;

using namespace Mistral;
using namespace std;

namespace CSPXMLParser {

  /**********************************************
   * Intentional Constraint
   **********************************************/ 

  class IntentionalConstraint : public Constraint {

  public:
    AST *predicate;
    int *assignment;
    Vector<int>* parameters;
    int nbparam;

    /**@name Constructors*/
    //@{
    IntentionalConstraint(Solver *s, VariableInt **v, VariableInt **vp, const int l, const int lp, AST *p)
      : Constraint(s, v, l)
    {
      nbparam = lp;
      predicate = p;
      assignment = new int[lp];
      std::fill( assignment, assignment+lp, NOVAL );
      parameters = new Vector<int>[arity];
      bool intparam;
      for(int i=0; i<lp; ++i) {
	intparam = true;
	for(int j=0; j<arity && intparam; ++j)
	  if( vp[i] == v[j] ) {
	    parameters[j].push( i );
	    intparam = false;
	  }
	if( intparam ) 
	  assignment[i] = vp[i]->first();
      }
      
      unsigned long int memusage = 1;
      for(int i=0; memusage < 1024 && i<l; ++i)
	{
	  memusage *= (v[i]->max() - v[i]->min() + 1);
	}
      if( memusage < 1024 )
	useResidualSupports(); 
    }
    virtual ~IntentionalConstraint()
    {
      delete [] assignment;
      delete [] parameters;
    }
    //@}
    

    /*!
      Check the corresponding cell in the matrix
    */
    inline int check( const int* s ) const 
    { 
      int i, j;
      for(i=0; i<arity; ++i) 
	for(j=0; j<parameters[i].size; ++j) 
	  assignment[parameters[i][j]] = s[i];
      int res = !(predicate->value( assignment ));
      return res;
    }

    virtual void print(std::ostream& o) const
    {
//       o << id << " ";
//       for(int i=0; i<arity; ++i) {
// 	scope[i]->print( o );
// 	o << " ";
//       }
      predicate->infixExpression( o );
      o << "(" ;
      int i,j,k;
      bool isvar;
      for(k=0; k<nbparam; ++k) {
	if( k ) o << ", ";
	isvar = false;
	for(i=0; i<arity; ++i)
	  for(j=0; j<parameters[i].size; ++j) {
	    if(k == parameters[i][j]) {
	      o << " " ;
	      scope[i]->printshort( o );
	      isvar = true;
	    }
	  }
	if( !isvar )
	  o << assignment[k] ;
      }
      o << " )";
    }
  
  };



  /**********************************************
   * Predicate Constraint Declaration
   **********************************************/
  class BuildObjectPredTree : public BuildObjectConstraint {
    
  private:
    int arity;
    AST *tree;
    
  public:

    BuildObjectPredTree( AST *t, int n ) { arity = n; tree = t; }
    
    virtual void build( Solver *s, VariableInt **tmp, BuildObjectPredicate *pred )
    {    
      new IntentionalConstraint( s, tmp, &tmp[arity], arity, pred->arity-arity, tree ); 
    }
 
    virtual void close( BuildObjectPredicate *pred ) 
    {
      //bool res = BuildObjectConstraint::preBuild( pred );
      bool continuous = ( ((pred->scope[0]->max() - pred->scope[0]->min() + 1) / 
			   pred->scope[0]->size()) <= 2 && 
			  ((pred->scope[1]->max() - pred->scope[1]->min() + 1) / 
			   pred->scope[1]->size()) <= 2 );
      if( continuous // && res
	  && arity == 2 ) {
	int i, j, k, length = pred->arity-arity;

	//std::cout << pred->arity << " " << arity << std::endl;

	//pred->print(std::cout);
	
	//std::cout << std::endl;

	int assignment[length];
	int vars[length];
	BuildObject **scope = pred->scope;
	k=0;
	for( i=0; i<length; ++i )
	  if( scope[arity+i] == scope[0] ||
	      scope[arity+i] == scope[1] )
	    vars[k++] = i;
	  else assignment[i] = scope[arity+i]->min();
	
	int borneinf = std::max( pred->scope[0]->min(), pred->scope[1]->min() );
	int bornesup = std::min( pred->scope[0]->max(), pred->scope[1]->max() );
	for( i=borneinf; i<=bornesup; ++i )
	  {
	    for( j=0; j<k; ++j )
	      assignment[vars[j]] = i;
	    if( tree->value( assignment ) ) break;
	  }
	if( i > bornesup ) 
	  pred->setNeq();
      }
      BuildObject **x = pred->scope;
      int n=pred->arity, i;  
      pred->deReference();
      pred->model->unsetSat();
      for(i=0; i<n; ++i) {
	x[i]->unsetRange();
      }
    }

    virtual void print(std::ostream& o, const BuildObjectPredicate *pred) const
    {
      tree->infixExpression( o );
    }

  };
  /**********************************************
   * Predicate Constraint Wrapper
   **********************************************/
  class Predicate : public Variable {

  public:
 
    Predicate( VarArray& x, VarArray& y, AST *t)
      {
	int n = x.size(), m = y.size();
	BuildObject **scope = new BuildObject*[n+m+1];
	for(int i=0; i<n; ++i) 
	  scope[i] = x[i].var_ptr_;
	for(int i=0; i<m; ++i) 
	  scope[i+n] = y[i].var_ptr_;
	trash.push( new BuildObjectPredTree(t,n) ); 
	//var_ptr_ = new BuildObjectPredicate( scope, n+m, 0, 1, trash.back(), NULL)
	// the following might be buggy!!!;
	var_ptr_ = new BuildObjectPredicate( scope, n+m, 0, 1, trash.back(), NULL);
      }
  
  };


// /**********************************************
//  * Predicate Tree Constraint BuildObject
//  **********************************************/

//   BuildObjectPredTree::BuildObjectPredTree( AST *t, int n ) 
//   { 
//     arity = n; 
//     tree = t; 
//   }
  
//   // virtual void add( CSP *model, const int l, 
//   // 		  BuildObject **scope, 
//   // 		  const int n ) 
//   // {
//   //   for(int i=arity; i<n; ++i)
//   //     scope[i]->add( model, l+1, (BuildObject::CRANGE | BuildObject::CLIST) );
//   // }
  
//   void BuildObjectPredTree::build( Solver *s, VariableInt **tmp, 
// 				   BuildObjectPredicate *pred )
//   {    
//     Constraint *con = new IntentionalConstraint( s, tmp, &tmp[arity], arity, 
// 						 pred->arity-arity, tree ); 
//   }
  
//   bool BuildObjectPredTree::preBuild( BuildObjectPredicate *pred ) 
//   {
//     bool res = BuildObjectConstraint::preBuild( pred );
//     bool continuous = ( ((pred->scope[0]->max() - pred->scope[0]->min() + 1) / 
// 			 pred->scope[0]->size()) <= 2 && 
// 			((pred->scope[1]->max() - pred->scope[1]->min() + 1) / 
// 			 pred->scope[1]->size()) <= 2 );
//     if( continuous && res && arity == 2 ) {
//       int i, j, k, length = pred->arity-arity;
//       int assignment[length];
//       int vars[length];
//       BuildObject **scope = pred->scope;
//       k=0;
//       for( i=0; i<length; ++i )
// 	if( scope[arity+i] == scope[0] ||
// 	    scope[arity+i] == scope[1] )
// 	  vars[k++] = i;
// 	else assignment[i] = scope[arity+i]->min();
      
//       int borneinf = std::max( pred->scope[0]->min(), pred->scope[1]->min() );
//       int bornesup = std::min( pred->scope[0]->max(), pred->scope[1]->max() );
//       for( i=borneinf; i<=bornesup; ++i )
// 	{
// 	  for( j=0; j<k; ++j )
// 	    assignment[vars[j]] = i;
// 	  if( tree->value( assignment ) ) break;
// 	}
//       if( i > bornesup ) 
// 	pred->setNeq();
//     }
//     return res;
//   }
  
//   void BuildObjectPredTree::print(std::ostream& o, 
// 				  const BuildObjectPredicate *pred) const
//   {
//     tree->infixExpression( o );
//   }
  

  // definition of the functions which are called by the parser when it
  // reads new information in the file. These functions are in charge of
  // creating the data structure that the solver needs to do his job
  class MistralCallback : public CSPParserCallback
  {
  public:
    static const int DYNAMIC = 13;
    
    // timing
    double time_start;
    double time_now;
    double time_last;
    double time_solving;
    double time_limit;
    int node_limit;
    int probe_sol;
    int probe_nodes;
    double probe_time;
    double probe_start;
    double probe_end;	  
	  
    // parameters/strategies
    unsigned int randomSeed;
    int    randomized;
    string heuristic;
    int    restart_policy;
    int    restart_base;
    double restart_factor;
    int    model_type;
    int    verbose;
    int    evar_limit;
    int    evar_size;
    int    infer_cliques;
    int    cliques_limit;
    int    revert_sat;

    int    use_sac;
    int    use_lds;
    int    probing;
    int    p_heuristic;
    int    p_valheur;
    int    p_iteration;
    int    p_limit;
    bool   updateLW;
    bool   updateIP;
    int    ext_feature;
    int    dom_split;
    int    all_solution;
    int    valheur;

    string model_name;
    string constraintReference;
    bool isEdgeFinder;

    // model
    CSP model;
    VarArray X;

    // constraints
    VarArray constraintScope;
    VarArray predicateScope;
    int *params;
    AST** predicates;
    int consArity;
    int *predicateCount;
    int *predicateScore;

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
    int numGACPred;
    int numDECPred;
    int numEXTCons;

    int numBinExt;
    int numNaryExt;
    int numLargeExt;

    int avgPredShape;
    int avgPredSize;
    int avgPredArity;

    // domains
    int **domains;
    //BitSet *domains;
    int *domainSize;
    int *domainMax;
    int *domainMin;

    // variables;
    int *solution;
    int *D;

    /// old relation stuff
    int *relationSpin;

    /// New relation stuff
    int *relationsIndex;



    // models stats
    int nlistvars;
    int nbitvars;
    int nboolvars;
    int nrangevars;
    int nconstants;

    int elistvars;
    int ebitvars;
    int eboolvars;
    int erangevars;
    int econstants;

    int maxArity;

    // utilities
    int countI, countJ, countK, maxCount, curCount, idcur;
    CSPDefinitionType curtype;

    Solver *cp_solver;
    SatSolver * sat_solver;

  public: 

    // features
    double *features_vec;

    MistralCallback() {
      numExtraVars = numVariables = numDomains = numGlobal
	= numRelations = numConstraints = numPredicates = 
	numWeightedSum = numAllDiff = numElement = numCumulative = numGlobalCardinality
	= avgPredShape = avgPredSize = avgPredArity = 
	numBinExt = numNaryExt = numLargeExt = 0;

      nlistvars = nbitvars = nboolvars = nrangevars = nconstants =
	elistvars = ebitvars = eboolvars = erangevars = econstants = 
	numGACPred = numDECPred = numEXTCons = maxArity = 0;

      features_vec = NULL;
      cp_solver = NULL;
      sat_solver = NULL;
      domains = NULL;
      domainSize = NULL;
      domainMax = NULL;
      domainMin = NULL;
      solution = NULL;
      D = NULL;
      relationSpin = NULL;
      relationsIndex = NULL;
      params = NULL;
      predicates = NULL;
      predicateCount = NULL;
      predicateScore = NULL;
    }
    virtual ~MistralCallback() 
      {
	//std::cout << "c (mistral) delete parser" << std::endl;

	cleanup();
      }

    void setRandomize       ( const int rd )   { randomized   = rd; }
    void setRandomSeed      ( const unsigned int rs ) { randomSeed  = rs; }
    void setFeatureExt      ( const int   fe ) { ext_feature  = fe; }
    void setModel           ( const int   mt ) { model_type   = mt; }
    void setTimeLimit       ( const double t ) { time_limit   =  t; }
    void setNodeLimit       ( const int    t ) { node_limit   =  t; }
    void setHeuristic       ( const string h ) { heuristic    =  h; }
    void setProbeHeuris     ( const int    h ) { p_heuristic  =  h; }
    void setPValHeur        ( const int  pvh ) { p_valheur  =  pvh; }
    void setEVarLimit       ( const int    t ) { evar_limit   =  t; }
    void setEVarSize        ( const int    s ) { evar_size    =  s; }
    void setDomainSplit     ( const int    t ) { dom_split    =  t; }
    void setAllSolution     ( const int    t ) { all_solution =  t; }
    void setValHeur         ( const int   vh ) { valheur  =  vh; }
    void setName            ( const string n ) { model_name   =  n; }
    void setRestartPolicy   ( string r ) {
      if( r == "no" ) restart_policy = -1;
      else if( r == "geom" ) restart_policy = GEOMETRIC;
      else if( r == "luby" ) restart_policy = LUBY;
      else if( r == "dyn" ) restart_policy = DYNAMIC;
      else restart_policy = NOVAL;
    }
    void setRestartBase  ( int b ) { restart_base = b; }
    void setRestartFactor( double f ) { restart_factor = f;	
      //1+((double)1 / (double)f);
    }
    void setVerbosity    ( int v ) { verbose       = v; }
    void setProbing      ( int p ) { probing       = p; }
    void setPIteration   ( int i ) { p_iteration   = i; }
    void setPLimit       ( int l ) { p_limit       = l; }
    void setUpdateLW     ( int l ) { updateLW      = l; }
    void setUpdateIP     ( int l ) { updateIP      = l; }
    void setSAC          ( int l ) { use_sac       = l; }
    void setLDS          ( int l ) { use_lds       = l; }
    void setInferCliques ( int l ) { infer_cliques = l; }
    void setCliquesLimit ( int l ) { cliques_limit = l; }
    void setRevertToSAT  ( int l ) { revert_sat    = l; }


    virtual void beginInstance(const string & name)
    {

//		 probesol = 0;
//		 probenodes = 0;
//		 probetime = 0.0;
		
		
      if( verbose )
	cout << "c \nc Parsing xml file" << endl;
//       numExtraVars = numVariables = numDomains = numGlobal
// 	= numRelations = numConstraints = numPredicates = 
// 	numWeightedSum = numAllDiff = numElement = numCumulative = numGlobalCardinality
// 	= avgPredShape = avgPredSize = avgPredArity = 
// 	numBinExt = numNaryExt = numLargeExt = 0;

//       nlistvars = nbitvars = nboolvars = nrangevars = nconstants =
// 	elistvars = ebitvars = eboolvars = erangevars = econstants = 
// 	numGACPred = numDECPred = numEXTCons = maxArity = 0;

      time_start = getRunTime();
      
    }

    virtual void beginDomainsSection(int nbDomains)
    {
      if( verbose ) {
	time_last = getRunTime();
	cout << "c" << setw(10) << " " << "domains....." ;
	cout.flush();
      }

      numDomains = nbDomains;
      domains = new int*[numDomains];
      //domains = new BitSet[numDomains];
      domainSize = new int[numDomains];
      domainMax  = new int[numDomains];
      domainMin  = new int[numDomains];
   
      countI = 0;
      curCount = 0;
      maxCount = nbDomains;
    }
  
    virtual void beginDomain(const string & name,
			     int idDomain,
			     int nbValue)
    {
      domainSize[countI] = nbValue;
      domains[countI] = new int[nbValue];
      countJ = 0;
    }

    void addDomainValue(int v)
    {
      domains[countI][countJ++] = v;
    }

    virtual void addDomainValue(int first,int last)
    {
      //if(last-first+1 != domainSize[countI])
      for(int k=first; k<=last; ++k)
	domains[countI][countJ++] = k;
      //else	
    }

    virtual void endDomain()
    {
      domainMax[countI] = -1*NOVAL;
      domainMin[countI] =    NOVAL;
      int k=domainSize[countI];
      while( k-- ) {
	if( domains[countI][k] > domainMax[countI] )
	  domainMax[countI] = domains[countI][k];
	if( domains[countI][k] < domainMin[countI] )
	  domainMin[countI] = domains[countI][k];
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

    /*!
     * end the definition of all domains
     */
    virtual void endDomainsSection()
    {
      assert( countI == numDomains );

      if( verbose ) {
	time_now = getRunTime();
	cout << setw(6) << setprecision(5)
	     << (time_now - time_last) << "\t(" << numDomains 
	     << " domains)" << endl;
	if( verbose > 2 )
	  for(int i=0; i<numDomains; ++i)
	    {
	      cout << "D" << i << ":" ;
	      for(int j=0; j<domainSize[i]; ++j)
		cout << " " << domains[i][j];
	      cout << endl;
	    }
	time_last = time_now;
      }
    }


    virtual void beginVariablesSection(int nbVariables)
    {
      if( verbose ) {
	//time_now = getRunTime();
	cout << "c" << setw(10) << " " << "variables..." ;
	cout.flush();
      }

      numVariables = nbVariables;
      X.resize( numVariables );
      solution = new int[numVariables];
      params = new int[numVariables+1];
      D = new int[numVariables];

      countI = 0;
      curCount = 0;
      maxCount = nbVariables;
    }
  
    virtual void addVariable(const string & name, int idVar,
			     const string & domain, int idDomain)
    {

      ++countI;
      D[idVar] = idDomain;
      solution[idVar] = domainMin[idDomain];



      X[idVar] = Variable( domains[idDomain], domainSize[idDomain],
			   domainMin[idDomain], domainMax[idDomain] );

      //if( idVar == 10 )
      //	cout << idVar << " " << solution[idVar] << endl;

      //cout << X[idVar] << endl;
      //X[idVar].print( cout );

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

      model.add( X );

      if( verbose ) {
	time_now = getRunTime();
	//if(time_now < 0.01)
	//time_now = 0;
	cout << setw(6) << setprecision(5)
	     << (time_now - time_last) << "\t(" << numVariables 
	     << " variables)" << endl;

// 	if( verbose > 2 )
// 	  {
// 	    cout << X << endl;
// 	  }
	time_last = time_now;
      }
    }


    virtual void beginRelationsSection(int nbRelations)
    {
      if( verbose ) {
	//time_now = getRunTime();
	cout << "c" << setw(10) << " " << "relations..." ;
	cout.flush();
      }

      numRelations = nbRelations;
      relationsIndex = new int[numRelations];
      relationSpin = new int[numRelations];

      countI = 0;
      curCount = 0;
      maxCount = nbRelations;
    }
      
  
    virtual void beginRelation(const string & name, int idRel,
			       int arity, int nbTuples, RelType relType) 
    {
      relationsIndex[countI] = //CSP::Table( arity, nbTuples );
	CSP::addTable( arity, nbTuples );
      //ENV.addConstraint( new BuildObjectTable(arity, nbTuples) );
      relationSpin[countI] = (relType == REL_SUPPORT);
      countJ=0;
    }
    
    virtual void addRelationTuple(int arity, int tuple[]) 
    {
      CSP::addTuple( relationsIndex[countI], tuple );
      //BuildObjectTable *tab = (BuildObjectTable*)(ENV[relationsIndex[countI]]);
      //tab->add( tuple );
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
	time_now = getRunTime();
	//if(time_now < 0.01)
	//time_now = 0;
	cout << setw(6) << setprecision(5)
	     << (time_now - time_last) << "\t(" << numRelations 
	     << " relations)" << endl;
	time_last = time_now;
      }
    }

    virtual void beginPredicatesSection(int nbPredicates)
    {
      if( verbose ) {
	//time_now = getRunTime();
	cout << "c" << setw(10) << " " << "predicates.." ;
	cout.flush();
      }
      
      numPredicates = nbPredicates;
      predicates = new AST*[numPredicates];
      predicateCount = new int[numPredicates];
      predicateScore = new int[numPredicates];
      std::fill(predicateCount, predicateCount+numPredicates, 0);
      std::fill(predicateScore, predicateScore+numPredicates, 0);
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

//        cout << "P" << countI << " = " << tree << " ";
//        tree->infixExpression( cout );
//        cout << endl;

      predicates[countI] = tree;
      if( tree )
	tree->eval(predicateCount[countI], 
		   predicateScore[countI]);
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
	time_now = getRunTime();
	//if(time_now < 0.01)
	//time_now = 0;
	cout << setw(6) << setprecision(5)
	     << (time_now - time_last) << "\t(" << numPredicates 
	     << " predicates)" << endl;
	time_last = time_now;
      }
    }

    virtual void beginConstraintsSection(int nbConstraints)
    {
      if( verbose ) {
	//time_now = getRunTime();
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

      //cout << "======= " << reference << " ========" << endl;

      if( arity > maxArity ) maxArity = arity;

      idcur = id;
      curtype = type;
      consArity = arity;
      constraintScope.clear();
      predicateScope.clear();
      for(int i=0; i<arity; ++i) 
	constraintScope.add( X[scope[i].getVarId()] );
      constraintReference=reference;

      //cout << constraintScope << endl;

    }

    virtual void constraintParameters(const ASTList &args) 
    {
      if (constraintReference=="global:cumulative")
	{
	  const AST &tasks=args[0];
	  const AST &limit=args[1];


	  int i=0, j, k, ntasks = tasks.size();
	  int minHorizon=NOVAL, maxHorizon=-NOVAL;
	  VarArray booltask;
	  int height[ntasks];
	  int duration[ntasks];
	  int **boolindex = new int*[ntasks];
	 
	  for(i=0; i<ntasks; ++i)
	    {
	      const AST &desc=tasks[i];
	      
	      if( desc.hasKey("origin") && desc["origin"].isVar() ) {
		params[i] = desc["origin"].getVarId();
// 		cout << params[i] << "[" << domainMin[D[params[i]]] << ".." 
// 		     << domainMax[D[params[i]]] << "] + ";
	      } else break;
	      
	      if( desc.hasKey("duration") && desc["duration"].isInteger() ) {
		duration[i] = desc["duration"].getInteger();
		if( minHorizon > domainMin[D[params[i]]] )
		  minHorizon = domainMin[D[params[i]]];
		if( maxHorizon < domainMax[D[params[i]]]+duration[i] )
		  maxHorizon = domainMax[D[params[i]]]+duration[i];
// 		cout << duration[i] << " / " ; 
	      } else break;
	      
	      if( desc.hasKey("height") && desc["height"].isInteger() ) {
		height[i] = desc["height"].getInteger();
// 		cout << height[i] << " " ;
	      } else break;

// 	      cout << endl;
	    }

	  if( i < ntasks ) 
	    {	    
	      cerr << "not implemented" << endl;
	      exit( 0 );	      
	    }


// 	  cout << "<" << minHorizon << ".." << maxHorizon << ">" << endl;

	  for(i=0; i<ntasks; ++i)
	    {
	      boolindex[i] = new int[maxHorizon - minHorizon + 1];
	      boolindex[i] -= minHorizon;
	      std::fill( boolindex[i]+minHorizon, boolindex[i]+maxHorizon+1, -1);
	    }

	  for(i=0; i<ntasks; ++i)
	    {
	      k=domainMax[D[params[i]]]+duration[i];
// 	      cout << "[" << domainMin[D[params[i]]] << ".." << k-1 << "]" ;
	      for(j=domainMin[D[params[i]]]; j<k; ++j)
		{ 
// 		  cout << ((j-duration[i])) << " ";
		  boolindex[i][j] = booltask.size();
		  if( (j-duration[i]) >= domainMin[D[params[i]]] && j < domainMax[D[params[i]]] ) {
		    if( duration[i] == 1 )
		      booltask.add( X[params[i]] == j );
		    else 
		      booltask.add( ((X[params[i]] <= j) && (X[params[i]] > (j-duration[i])))  );
		  }
		  else if( j < domainMax[D[params[i]]] )
		    booltask.add( X[params[i]] <= j );
		  else 
		    booltask.add( X[params[i]] > (j-duration[i]) );
		    //booltask.add( (j-duration[i]) < X[params[i]] );
		}
// 	      cout << endl;
	    }

	  //MistralSet printtasks(0, booltask.size()-1, MistralSet::empt);


	  int weights[maxHorizon-minHorizon+1];
	  for(j=minHorizon; j<maxHorizon; ++j)
	    {

 	      //cout << " ---- " << j << " ---- " << endl;

	      VarArray sumscope;
	      for(i=0; i<ntasks; ++i)
		{
		  if( boolindex[i][j] >= 0 ) {
		    weights[sumscope.size()] = height[i];
		    sumscope.add( booltask[boolindex[i][j]] );
		    
		    //printtasks.insert( boolindex[i][j] );
		  }
		}
	      weights[sumscope.size()]=0;

	      //cout << sumscope << endl;

	      if( limit.isVar() )
		model.add( Sum( sumscope, weights ) <= X[limit.getVarId()] );
	      else 
		model.add( Sum( sumscope, weights ) <= limit.getInteger() );

	    }

	  for(i=0; i<ntasks; ++i)
	    {
	      boolindex[i] += minHorizon;
	      delete [] boolindex[i];
	    }
	  delete [] boolindex;

	}
      else if (constraintReference=="global:element")
	{
	  const AST &index=args[0];
	  const AST &table=args[1];
	  const AST &value=args[2];

	  if( index.isVar() )
	    predicateScope.add( X[index.getVarId()] );
	  else
	    predicateScope.add( index.getInteger() );

	  for(int i=0;i<table.size();++i)
	    if (table[i].isVar())
	      predicateScope.add( X[table[i].getVarId()] );
	    else
	      predicateScope.add( table[i].getInteger() );

	  if( value.isVar() )
	    predicateScope.add( X[value.getVarId()] );
	  else
	    predicateScope.add( value.getInteger() );
	
	}
      else if (constraintReference=="global:weightedsum")
	{
	  const AST &sum=args[0];
	  const AST &op=args[1];
	  const AST &rhs=args[2];


	  for(int i=0; i<sum.size(); ++i) {
	    params[i] = sum[i]["coef"].getInteger() ;
	    predicateScope.add( X[sum[i]["var"].getVarId()] );
	  }
	  params[sum.size()] = 0;


	  if( op.getType() == SYMB_EQ )
	    {
	      if( rhs.isVar() )
		{
		  model.add( Sum( predicateScope, params ) == X[rhs.getVarId()] );
		}
	      else
		{
		  model.add( Sum( predicateScope, params ) == rhs.getInteger() );
		}
	    } 
	  else if( op.getType() == SYMB_GE )
	    {
	      if( rhs.isVar() )
		{
		  model.add( Sum( predicateScope, params ) >= X[rhs.getVarId()] );
		}
	      else
		{
		  model.add( Sum( predicateScope, params ) >= rhs.getInteger() );

		  //cout << "HERE" << endl;
		  //exit( 0 );

		}
	    }
	  else if( op.getType() == SYMB_GT )
	    {
	      if( rhs.isVar() )
		{
		  model.add( Sum( predicateScope, params ) > X[rhs.getVarId()] );
		}
	      else
		{
		  model.add( Sum( predicateScope, params ) > rhs.getInteger() );
		}
	    }
	  else if( op.getType() == SYMB_LE )
	    {
	      if( rhs.isVar() )
		{
		  model.add( Sum( predicateScope, params ) <= X[rhs.getVarId()] );
		}
	      else
		{
		  model.add( Sum( predicateScope, params ) <= rhs.getInteger() );
		}
	    }
	  else if( op.getType() == SYMB_LT )
	    {
	      if( rhs.isVar() )
		{
		  model.add( Sum( predicateScope, params ) < X[rhs.getVarId()] );
		}
	      else
		{
		  model.add( Sum( predicateScope, params ) < rhs.getInteger() );
		}
	    }
	}
      else if (constraintReference=="global:global_cardinality")
	{
// 	  const AST &vars=args[0];
// 	  const AST &cards=args[1];
// 	  int i, j, nvars = vars.size(), ncards = cards.size();
// 	  const AST &desc=cards[0];	      
// 	  if( desc.hasKey("value") ) //&& desc["value"].isInteger() ) 
// 	    cout << "c ok for values" << endl;
// 	  if( desc.hasKey("card") && desc["card"].isVar() ) 
// 	    cout << "c ok for cards" << endl;
	}
      else 
	{
// 	  // default
//  	  cout << "constraint parameters=";
//  	  args.postfixExpression(cout);
//  	  cout << endl;
	  for(int i=0;i<args.size();++i) 
	    if( args[i].isVar() )
	      predicateScope.add( X[args[i].getVarId()] );
	    else if( args[i].isInteger() )
	      predicateScope.add( args[i].getInteger() );
	    else {
	      
	      for(int j=0;j<args[i].size();++j) 
		if( args[i][j].isVar() )
		  predicateScope.add( X[args[i][j].getVarId()] );
		else if( args[i][j].isInteger() )
		  predicateScope.add( args[i][j].getInteger() );
	      
	    }
	}


      //cout << predicateScope << endl;
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
	// THOSE ARE THE EXTENSIONAL CONSTRAINTS
	
	Table con( constraintScope, relationsIndex[idcur], relationSpin[idcur] );
	model.add( con );
  
      } else if( curtype == PredicateType ) { 


	//cout << " PREDICATE " << endl;
	
	/******************************************************
	 * We decide what representation to use amongst:
	 *  + tree of predicate constraints
	 *  + generic AC algorithm using a predicate
	 *  + both methods simultaneously
	 * We use the following criteria:
	 *  - arity
	 *  - size of the cartesian product of the domains
	 *  - size of the predicate tree
	 *******************************************************/

	switch(model_type) {
	case 0 : {
	  ++numDECPred;
	  model.add( predicates[idcur]->
		     getPredicate( predicateScope ) );
// 	  model.add( Variable( predicates[idcur]->
// 			       getPredicate( predicateScope ) ) );
	  break;
	}
	case 1 : {
	  ++numGACPred;
	  model.add( Predicate( constraintScope,
				predicateScope,
				predicates[idcur] ) );
	  break;
	}
	case 2 : {
	  ++numGACPred;
	  ++numDECPred;
	  model.add( Variable( predicates[idcur]->
			       getPredicate( predicateScope ) ) );
	  if( consArity < 7 )
	    model.add( Predicate( constraintScope,
				  predicateScope,
				  predicates[idcur] ) );
	  break;
	}
	default : {

	  
	  // helpers
	  int aux, i = consArity;
	  long unsigned cart_cont = 1;

	  /// Attributes:
	  int predArity = predicateScope.size();
	  int rPredArity = predArity;
	  int maxDomSize = 0;
	  long unsigned cartesian_p = 1;
	  int continuity;

	  /// compute cartesian_p, maxDomSize and  continuity
	  while( i-- && cartesian_p <= 10000000 ) {
	    aux = constraintScope[i].var_ptr_->size();
	    if(maxDomSize < aux) maxDomSize = aux;
	    cartesian_p *= aux;
	    cart_cont *= (constraintScope[i].var_ptr_->max() -
			  constraintScope[i].var_ptr_->min() + 1);
	  }
	  continuity = (int)((double)(cart_cont)/(double)(cartesian_p));

	  // compute rPredArity
	  i = predArity;
	  while( i-- )
	    if( predicateScope[i].var_ptr_->size() < 2 )
	      --rPredArity;
	  

	  int predscore = 0;
	  int igacscore = 0;
	  int lastps = 0;
	  int lastgs = 0;


	  if( verbose > 2 ) {
	    predicates[idcur]->infixExpression( cout );
	    cout << endl;
	  }
// 	  if( verbose > 2 ) {
// 	    cout << predicateScope << endl;
// 	    cout << constraintScope << " " << cartesian_p << endl;
// 	  }

	  // case where we use predicate trees unconditionally
	  if( consArity == 1 || cartesian_p > 10000000 )
	    predscore = NOVAL;
	  else {


	    /********************************************
	     * Arity
	     *
	     * GAC is better for low arity constraint
	     ********************************************/

	    predscore += (consArity - 2);

	    if( verbose > 2 ) {
	      cout << "Arity:             (" << (predscore - lastps)
		   << " / " << (igacscore - lastgs) << ")" << endl;
	      lastps = predscore;
	      lastgs = igacscore;
	    }

	    avgPredArity += consArity;

	    /********************************************
	     * "Bad" predicates
	     *
	     * Mul, Div, Mod, Or etc are not very
	     * efficient nor very well implemented,
	     * we should rather use GAC
	     ********************************************/

	    igacscore += ( 3 * (int)((double)(predicateScore[idcur]) / 
				     (double)(predicateCount[idcur])) );

	    if( verbose > 2 ) {
	      cout << "Bad predicates:    (" << (predscore - lastps)
		   << " / " << (igacscore - lastgs) << ")" << endl;
	      lastps = predscore;
	      lastgs = igacscore;
	    }

	    avgPredShape += (igacscore - lastgs);

	    /********************************************
	     * Shape of tree
	     *
	     * Flat trees are better popagated by the decomposition
	     * whilst large trees with a few leaves are better 
	     * propagated using GAC
	     ********************************************/
	    
	    if( predicateCount[idcur] < consArity ) predscore += 1000;
	    igacscore += ( (int)((double)(predicateCount[idcur]) / 
				 (double)(consArity)) );

	    if( verbose > 2 ) {
	      cout << "Shape of the tree: (" << (predscore - lastps)
		   << " / " << (igacscore - lastgs) << ")" << endl;
	      lastps = predscore;
	      lastgs = igacscore;
	    }

	    avgPredSize += predicateCount[idcur];

	    /********************************************
	     * Continuity
	     *
	     * Domains with holes are often not very well
	     * handled by the decomposition as several
	     * predicates enforces BC instead of GAC.
	     * GAC is much better in this case
	     ********************************************/

	    igacscore += (continuity - 1);

	    if( verbose > 2 ) {
	      cout << "Continuity:        (" << (predscore - lastps)
		   << " / " << (igacscore - lastgs) << ")" << endl;
	      lastps = predscore;
	      lastgs = igacscore;
	    }


	    /********************************************
	     * Search Size
	     *
	     * Large search space, i.e., numerous and/or
	     * big domains may be very costly for a 
	     * generic GAC algorithm. 
	     * The decomposition is better in this case
	     ********************************************/
	    
	    predscore += (int)(log2( cartesian_p ) - 3);
	    
	    if( verbose > 2 ) {
	      cout << "Search Size:       (" << (predscore - lastps)
		   << " / " << (igacscore - lastgs) << ")" << endl;
	      lastps = predscore;
	      lastgs = igacscore;
	    }


	    /********************************************
	     * Domain size
	     *
	     * On Boolean variables the predicates
	     * that enforce BC are equivalent to GAC
	     ********************************************/

	    predscore += 5 * ( maxDomSize <= 2 );
	    
	    if( verbose > 2 ) {
	      cout << "Domain size:       (" << (predscore - lastps)
		   << " / " << (igacscore - lastgs) << ")" << endl;
	      lastps = predscore;
	      lastgs = igacscore;
	    }


	    /********************************************
	     * Number of constraints
	     *
	     * GAC is more memory efficient
	     ********************************************/

	    if( model_type > 3 ) {
	      igacscore += (numConstraints > 32000)*12 + (numConstraints > 64000)*3;
	      
	      if( verbose > 2 ) {
		cout << "Problem size:       (" << (predscore - lastps)
		     << " / " << (igacscore - lastgs) << ")" << endl;
		lastps = predscore;
		lastgs = igacscore;
	      }
	    }


	  }

	    
	  if( predscore > igacscore ) {
	    
// 	    cout << " ==> DECOMPOSITION " << endl << endl; 
// 	    predicates[idcur]->infixExpression( cout );
// 	    cout << endl;

	    if( verbose > 2 ) {
	      cout << " ==> DECOMPOSITION " << endl << endl; 
	    }
	    
	    ++numDECPred;
	    if( predicates[idcur] )
	      model.add( Variable( predicates[idcur]->
				   getPredicate( predicateScope ) ) );
	  } else {
	    
	    if( verbose > 2 ) {
	      cout << " ==> GAC/PREDICATE" << endl << endl; 
	    }

// 	    cout << " ==> GAC/PREDICATE" << endl << endl; 
	    
	    ++numGACPred;	    
	    model.add( Predicate( constraintScope,
				  predicateScope,
				  predicates[idcur] ) );
	  }
 	}
	}

      } else if( curtype == GlobalConstraintType ) {

	++numGlobal;

	if (constraintReference=="global:alldifferent") {
	  ++numAllDiff;
	  model.add( AllDifferent(constraintScope) ); 
	} else if (constraintReference=="global:element") {
	  ++numElement;
	  model.add( Element(predicateScope) );
	} else if (constraintReference=="global:weightedsum") {
	  ++numWeightedSum;
	} else if (constraintReference=="global:cumulative" ) {
	  ++numCumulative;
	} else if (constraintReference=="global:global_cardinality" ) {
	  ++numGlobalCardinality;
	  cout << endl << "s NOT SUPPORTED" << endl;
	  exit(0);
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
      if( numRelations ) {
	delete [] relationSpin;
	delete [] relationsIndex;
      }

      for(int k=0; k<numDomains; ++k)
	delete [] domains[k];
      delete [] domains;
      if( numPredicates ) {
	delete [] predicateCount;
	delete [] predicateScore;
      }
      if( verbose ) {
	time_now = getRunTime();
	//if(time_now < 0.01)
	//time_now = 0;
	cout << setw(6) << setprecision(5)
	     << (time_now - time_last) << "\t(" << numConstraints 
	     << " constraints)" << endl;
	time_last = time_now;
      }
    }

    /********************************************************************/


    void addHeuristic( Solver& s ) {
      if( heuristic == "lex" ) {
	Lexicographic h;
	s.add( h );
      } else if( heuristic == "deg") {
	MaxDegree h(abs(randomized));
	s.add( h );
      } else if( heuristic == "rand") {
	Random h;
	s.add( h );
      } else if( heuristic == "dom+deg") {
	MinDomMaxDeg h(abs(randomized));
	s.add( h );
      } else if( heuristic == "dom/deg") {
	DomOverDeg h(abs(randomized));
	s.add( h );
      } else if( heuristic == "dom/wldeg") {
	DomOverWLDeg h(abs(randomized));
	s.add( h );
      } else if( heuristic == "dom/wdeg") {
	DomOverWDeg h(abs(randomized));
	s.add( h );
      } else if( heuristic == "neighbor") {
	Neighbor h(abs(randomized));
	s.add( h );
      } else if( heuristic == "impact") {
	Impact h(abs(randomized));
	s.add( h );
      } else if( heuristic == "impact/wdeg") {
	ImpactOverWDeg h(abs(randomized));
	s.add( h );
      } else if( heuristic == "impact/wldeg") {
	ImpactOverWLDeg h(abs(randomized));
	s.add( h );
      }
      else {
	MinDomain h(abs(randomized));
	s.add( h );
      }
    }
    
    
    void solveWithSAT( int& result ) {
      sat_solver = new SatSolver( model );
      sat_solver->params.verbosity = verbose+1;
      sat_solver->setPolicy( GEOMETRIC );
      result = sat_solver->solve();
    }

//       if( result == SAT ) {
// 	cout << "s SATISFIABLE" << endl << "v";
// 	for(int k=0; k<X.size(); ++k) {
// 	  if( domainMin[D[k]] == domainMax[D[k]] ) {
// 	    cout << " " << domainMin[D[k]];
// 	  } else cout << " " << ( sat_solver.polarity[k+1] == k+1 ) ;
// 	}
// 	cout << endl;
//       } else if ( result == UNSAT ) {
// 	cout << "s UNSATISFIABLE" << endl;
//       } else {
// 	cout << "s UNKNOWN" << endl;
//       }
//       double time = (getRunTime() - sat_solver.stats.start_time);
//       cout << "d CHECKS " << sat_solver.stats.unit_props << endl
// 	   << "d ASSIGNMENTS " << sat_solver.stats.nodes << endl
// 	   << "d NODES " << sat_solver.stats.nodes 
// 	   << " BACKTRACKS " << sat_solver.stats.conflicts
// 	   << " FAILURES " << sat_solver.stats.conflicts
// 	   << " RUNTIME " << time
// 	   << " MISC " << 0
// 	   << " TOTALTIME " << time
// 	   << " NODES/s " << (long unsigned int)(time > 0.001 ? (double)(sat_solver.stats.nodes)/time : 100*sat_solver.stats.nodes) 
// 	   << " CHECKS/s " << (long unsigned int)(time > 0.001 ? (double)(sat_solver.stats.unit_props)/time : 100*sat_solver.stats.unit_props) << endl ;
      
//    }

	  
    void solveWithCP( int& result, Solver *s=NULL ) {
      
      cp_solver = (s ? s : new Solver());

      if( dom_split ) 
	cp_solver->setDomainSplitting();
	      
      if( infer_cliques && model.nNeq > 3 && numGlobal == 0 // && model.nNeq < 8192
	  ) {
	int nalldiff = model.inferAllDiffs( (infer_cliques > 1), cliques_limit );
	numGlobal += nalldiff;
	if( verbose ) 
	  cout << "c" << setw(10) << nalldiff << " AllDifferent constraints infered" << endl;
      }
	      
      if( model.numEVars() > evar_limit ) {
	cp_solver->build( model, X.size() );
	cp_solver->initSearch( X, 0 );
      } else {
	cp_solver->build( model );
	cp_solver->initSearch( X, evar_size );
      }
	  
      cp_solver->INITTIME = time_start;
      cp_solver->setTimeLimit( time_limit ); 
      if( node_limit > 0 )
	cp_solver->setNodeLimit( node_limit );
      cp_solver->setVerbosity( verbose );
	  
      if( restart_base <= 0 ) 
	restart_base = std::min(1000, X.size());

      unsigned int nvalues = 0;


	  
      /********************************************
       *
       * Dynamic heuristics, branching, 
       *
       ********************************************/
	  
      if( heuristic == "dyn" ) 
	{
	  nvalues = model.nValues + model.eValues;
	      
	  if( numEXTCons == numConstraints || nvalues > 16000 )
	    heuristic = "dom/wdeg" ;
	  else
	    heuristic = "impact/wdeg" ;
	      
	  if( use_sac == 3 && nvalues < 12000 )
	    use_sac = 1;
	      
	  if( dom_split == 2 && true )
	    dom_split = 1;
	}
	  
	  
      if( verbose ) {
	cout << "c" << setw(10) << (model.nConstant) << " constants\n";
	cout << "c" << setw(10) << (model.nBoolean) << " Boolean variables    (+" 
	     << (model.eBoolean) << ")\n";
	cout << "c" << setw(10) << (model.nRange)   << " range variables      (+" 
	     << (model.eRange) << ")\n";
	cout << "c" << setw(10) << (model.nBit)     << " domain (b) variables (+" 
	     << (model.eBit) << ")\n";
	cout << "c" << setw(10) << (model.nList)       << " domain (l) variables (+" 
	     << (model.eList) << ")\n";
	cout << "c" << setw(10) << (model.nConstraint) << " constraints          (+" 
	     << (cp_solver->constraints.size - model.nConstraint) << ")\n";
	cout << "c" << setw(10) << (model.nValues)  << " values \n";


#ifdef _LINUX
	memory_usage = getMemory();
#endif
	    
	time_now = getRunTime();
	cout << "c Building time: " << setw(6)
	     << setprecision(5) << (time_now - time_last)
#ifdef _LINUX
	     << "  (" << (memory_usage >> 8) << " Mb)"  
#endif
	     << "\nc\nc Solving instance" << endl;
	time_solving = time_last = time_now; 
      }


      if( restart_policy == DYNAMIC )
	{
	  //restart_policy = ( (numGACPred + numDECPred) ? GEOMETRIC : -1 );
	  restart_policy = ( (numPredicates || numGlobal) ? GEOMETRIC : -1 );
	}
	  

      /********************************************
       *
       * Solving
       *
       ********************************************/
	  
      if( verbose ) {
	cout << "c" << setw(10) << " " << "heuristic = " << heuristic 
	     << (dom_split ? " (domain splitting)" : "" ) << endl;
	if( probing ) {
	  cout << "c" << setw(10) << " " << (randomized ? "Randomized " : "" ) 
	       << "probing iteration = " << p_iteration << endl;
	  cout << "c" << setw(10) << " " << "probing limit = " << p_limit << endl;
	  cout << "c" << setw(10) << " " << "level based weights = " << (updateLW ? "yes" : "no") << endl;
	  cout << "c" << setw(10) << " " << "impact probing = " << (updateIP ? "yes" : "no") << endl;
	}
	if( use_sac )
	  cout << "c" << setw(10) << " " << "singleton arc consistency pass " << endl;
	cout << "c" << setw(10) << " " << "restart policy = " ;
	switch( restart_policy ) {
	case -1 : cout << "No restart" << endl; break;
	case GEOMETRIC : cout << (randomized ? "Randomized / " : "" ) << "Geometric increment (seed = " << randomSeed << ")" << endl; break;
	case LUBY : cout << (randomized ? "Randomized / " : "" ) << "Luby policy (seed = " << randomSeed << ")" << endl; break;
	default : cout << "Unknown policy" << endl; break;
	}
	if( restart_policy >= 0 )
	  cout << "c" << setw(10) << " " << "restart base = " << restart_base << endl;
	if( restart_policy == GEOMETRIC )
	  cout << "c" << setw(10) << " " << "restart factor = " << restart_factor << endl;	
      }
      
      cp_solver->setRandomSeed( randomSeed );
      cp_solver->setRandomized( (randomized > 0) );
      
//       cp_solver->print( cout );
//       cout << endl;
	      
	  
      /********************************************
       * Probing for cpzilla
       ********************************************/
	  
      if( all_solution )
	{
	      
	  cp_solver->FIND_ALL = -1;
	  cp_solver->solve();
	  cout << "d "
	       << "SOLUTIONS " << cp_solver->SOLUTIONS
	       << " CHECKS " << cp_solver->CHECKS 
	       << " ASSIGNMENTS " << cp_solver->NODES << endl;
	      
	}
      else if( ext_feature ) 
	{
	  features_vec = new double[36];
	      
	  features_vec[0] = (nconstants ? (log2(nconstants)) : -1);
	  features_vec[1] = ((model.nBoolean) ? (log2(model.nBoolean)) : -1);
	  features_vec[2] = ((model.nRange) ? (log2(model.nRange)) : -1);
	  features_vec[3] = ((model.nBit) ? (log2(model.nBit)) : -1);
	  features_vec[4] = ((model.nList) ? (log2(model.nList)) : -1);
	  features_vec[5] = ((model.nValues) ? (log2(model.nValues)) : -1);
	  features_vec[6] = ((model.eBoolean) ? (log2(model.eBoolean)) : -1);
	  features_vec[7] = ((model.eRange) ? (log2(model.eRange)) : -1);
	  features_vec[8] = ((model.eBit+model.eList) ? (log2(model.eBit+model.eList)) : -1);
	  features_vec[9] = ((model.eValues) ? (log2(model.eValues)) : -1);
	  features_vec[10] = ((cp_solver->constraints.size) ? (log2(cp_solver->constraints.size)) : -1);
	  features_vec[11] = ((cp_solver->length) ? (log2(cp_solver->length)) : -1);

	  /*
	  std::cout << (ext_feature > 2 ? "\nc (mistral) 0 constants:        " : "") << features_vec[0] << " ";
	  std::cout << (ext_feature > 2 ? "\nc (mistral) 1 booleans:         " : "") << features_vec[1] << " ";
	  std::cout << (ext_feature > 2 ? "\nc (mistral) 2 ranges:          " : "") << features_vec[2] << " ";
	  std::cout << (ext_feature > 2 ? "\nc (mistral) 3 bits:             " : "") << features_vec[3] << " ";
	  std::cout << (ext_feature > 2 ? "\nc (mistral) 4 lists:            " : "") << features_vec[4] << " ";
	  std::cout << (ext_feature > 2 ? "\nc (mistral) 5 values:           " : "") << features_vec[5] << " ";
	  std::cout << (ext_feature > 2 ? "\nc (mistral) 6 extra bools:      " : "") << features_vec[6] << " ";
	  std::cout << (ext_feature > 2 ? "\nc (mistral) 7 extra ranges:     " : "") << features_vec[7] << " ";
	  std::cout << (ext_feature > 2 ? "\nc (mistral) 8 extra domains:    " : "") << features_vec[8] << " ";
	  std::cout << (ext_feature > 2 ? "\nc (mistral) 9 extra values:     " : "") << features_vec[9] << " ";
	  std::cout << (ext_feature > 2 ? "\nc (mistral)10 constraints:      " : "") << features_vec[10] << " ";
	  std::cout << (ext_feature > 2 ? "\nc (mistral)11 search variables: " : "") << features_vec[11] << " ";
	  */
	  addHeuristic( *cp_solver );
	      
	  cp_solver->setTimeLimit( 1.457 );
	  result = cp_solver->solve();
	      
	  int i;
	      
	  double avgweight=0, stddev=0;
	  int n=cp_solver->length; 
	  for(i=0; i<n; ++i)
	    avgweight += cp_solver->variables[i]->weight ;	   
	  avgweight /= (double)n;
	  for(i=0; i<n; ++i) {
	    //cout << "\t\t\t + " << ((avgweight - (double)(cp_solver->variables[i]->weight)) * (avgweight - (double)(cp_solver->variables[i]->weight)))
	    //<< " = " ;
	    stddev += ((avgweight - (double)(cp_solver->variables[i]->weight)) * (avgweight - (double)(cp_solver->variables[i]->weight)));
	    //cout << stddev << endl;
	  }
	      
	  //cout << " ==> " << ((stddev)/(double)n) << " :=: " ;
	  stddev = (sqrt((stddev)/(double)n));
	  //cout << stddev << endl;

	  features_vec[12] = ((avgweight) ? (log2(avgweight)) : -1);
	  features_vec[13] = ((stddev) ? (log2(stddev)) : -1);
	  features_vec[14] = ((cp_solver->NODES) ? (log2(cp_solver->NODES)) : -1);
	  features_vec[15] = ((cp_solver->PROPAGS) ? (log2(cp_solver->PROPAGS)) : -1);
	  features_vec[16] = maxArity;
	  features_vec[17] = ((double)numEXTCons*100.0/(double)numConstraints);
	  features_vec[18] = ((double)numGlobal*100.0/(double)numConstraints);
	  features_vec[19] = ((double)(numGACPred)*100.0/(double)numConstraints);
	  features_vec[20] = ((double)(numDECPred)*100.0/(double)numConstraints);

	  /*
	  cerr << (ext_feature > 2 ? "\n12 avg weight:       " : "") << features_vec[12] << " ";
	  cerr << (ext_feature > 2 ? "\n13 std deviation:    " : "") << features_vec[13] << " ";
	  cerr << (ext_feature > 2 ? "\n14 nodes:            " : "") << features_vec[14] << " ";
	  cerr << (ext_feature > 2 ? "\n15 checks:           " : "") << features_vec[15] << " ";
	  cerr << (ext_feature > 2 ? "\n16 max arity:        " : "") << features_vec[16] << " ";
	  cerr << (ext_feature > 2 ? "\n17 % extensional:    " : "") << features_vec[17] << " ";
	  cerr << (ext_feature > 2 ? "\n18 % global:         " : "") << features_vec[18] << " ";	  
	  cerr << (ext_feature > 2 ? "\n19 % bad preds:      " : "") << features_vec[19] << " ";
	  cerr << (ext_feature > 2 ? "\n20 % good preds:     " : "") << features_vec[20] << " ";
	  */

	  double avgDomSize = 0;
	  double maxDomSize = 0;
	  double avgDomCont = 0;
	  double minDomCont = NOVAL;
	  double auxCont;
	  for(i=0; i<numVariables; ++i) {
	    avgDomSize += domainSize[D[i]];
	    if(domainSize[D[i]] > maxDomSize) maxDomSize = domainSize[D[i]]; 
	    auxCont = ((double)domainSize[D[i]]/(double)(domainMax[D[i]]-domainMin[D[i]]+1));
	    avgDomCont += auxCont;
	    if( auxCont < minDomCont ) minDomCont = auxCont;
	  }	 
	  avgDomSize /= (double)(numVariables);
	  avgDomCont /= (double)(numVariables);
	      

	  features_vec[21] = (sqrt(avgDomSize));
	  features_vec[22] = (sqrt(maxDomSize));
	  features_vec[23] = (100.0 * avgDomCont);
	  features_vec[24] = (100.0 *minDomCont);
	  features_vec[25] = numAllDiff;
	  features_vec[26] = (numEXTCons ? 
			      ((double)numBinExt*10.0/(double)numEXTCons)
			      : -1);
	  features_vec[27] = (numEXTCons ? 
			      ((double)numNaryExt*10.0/(double)numEXTCons)
			      : -1);
	  features_vec[28] = (numEXTCons ? 
			      ((double)numLargeExt*10.0/(double)numEXTCons)
			      : -1);
	  features_vec[29] = (numGlobal ? 
			      ((double)numAllDiff*100.0/(double)numGlobal)
			      : -1);
	  features_vec[30] = (numGlobal ? 
			      ((double)numWeightedSum*100.0/(double)numGlobal)
			      : -1);
	  features_vec[31] = (numGlobal ? 
			      ((double)numElement*100.0/(double)numGlobal)
			      : -1);
	  features_vec[32] = (numGlobal ? 
			      ((double)numCumulative*100.0/(double)numGlobal)
			      : -1);
	  features_vec[33] = ((numDECPred+numGACPred) ?
			      ((double)avgPredShape*10.0 / (double)(numDECPred+numGACPred)) 
			      : -1);
	  features_vec[34] = ((numDECPred+numGACPred) ?
			      ((double)avgPredSize*10.0 / (double)(numDECPred+numGACPred)) 
			      : -1);
	  features_vec[35] = ((numDECPred+numGACPred) ?
			      ((double)avgPredArity*10.0 / (double)(numDECPred+numGACPred)) : -1);


	  //cerr << "log_constants log_booleans log_ranges log_bits log_lists log_values log_extra_booleans log_extra_ranges log_extra_bits log_extra_values log_constraints log_search_vars log_avg_weight log_stdev_weight log_nodes log_propags max_arity percent_ext percent_global percent_gac_predicate percent_dec_predicate sqrt_avg_domsize sqrt_max_domsize percent_avg_continuity percent_min_continuity num_alldiff perten_naryext perten_largeext percent_alldiff percent_wsum percent_element percent_cumulative perten_avg_predshape perten_avg_predsize perten_avg_predarity"

	  /*
	  cerr << (ext_feature > 2 ? "\n21 avg dom size:     " : "") << features_vec[21] << " " ;
	  cerr << (ext_feature > 2 ? "\n22 max dom size:     " : "") << features_vec[22] << " ";
	  cerr << (ext_feature > 2 ? "\n23 avg continuity:   " : "") << features_vec[23] << " " ;
	  cerr << (ext_feature > 2 ? "\n24 min continuity:   " : "") << features_vec[24] << " ";
	  cerr << (ext_feature > 2 ? "\n25 num alldiff:      " : "") << features_vec[25] << " ";
	  cerr << (ext_feature > 2 ? "\n26 ratio binary ext  " : "") << features_vec[26] << " ";
	  cerr << (ext_feature > 2 ? "\n27 ratio nary ext:   " : "") << features_vec[27] << " ";
	  cerr << (ext_feature > 2 ? "\n28 ratio large ext:  " : "") << features_vec[28] << " ";
	  cerr << (ext_feature > 2 ? "\n29 ratio alldiff:    " : "") << features_vec[29] << " ";
	  cerr << (ext_feature > 2 ? "\n30 ratio wsum:       " : "") << features_vec[30] << " ";
	  cerr << (ext_feature > 2 ? "\n31 ratio element:    " : "") << features_vec[31] << " ";
	  cerr << (ext_feature > 2 ? "\n32 ratio cumulative: " : "") << features_vec[32] << " ";
	  cerr << (ext_feature > 2 ? "\n33 predicate score:  " : "") << features_vec[33] << " " ;
	  cerr << (ext_feature > 2 ? "\n34 predicate size:   " : "") << features_vec[34] << " " ;
	  cerr << (ext_feature > 2 ? "\n35 predicate arity:  " : "") << features_vec[35] << " " ;
	  cerr << "mistral 1801 abscon 1801 choco 1801" << endl;	  
	  */
	}
      else  
	{

	  /********************************************
	   * Singleton Arc Consistency
	   ********************************************/
	      
	  if( verbose > 1 && use_lds ) {
	    cout << "c Initial impact: " << endl;
	    for(int i=0; i<cp_solver->learners.size; ++i)
	      {
		cp_solver->learners[i]->print( cout );
		cout << "c " << endl;
	      }
	  }
	  
	  if( use_sac )
	    {      
	      if( verbose ) {
		cout << "c" << setw(10) << " " << "singleton arc consistency " ;
		cout.flush();
	      }
	      
	      result = cp_solver->sacPreprocess( (use_sac > 1) );

	      if( verbose ) {
		time_now = getRunTime(); 
		cout << setw(6)
		     << setprecision(5) << (time_now - time_last) << endl;
		time_last = time_now;

		if( use_lds && verbose > 1 ) {
		  cout << "c Impact learned during sac: " << endl;
		  for(int i=0; i<cp_solver->learners.size; ++i)
		    {
		      cp_solver->learners[i]->print( cout );
		      cout << "c " << endl;
		    }
		}	      
	      }
	    }
	  
	  if( result == UNKNOWN ) 
	    {
	      /********************************************
	       * Probing
	       ********************************************/
	      if( probing ) 
		{
		  /*********************************************/
		  if( verbose ) {
		    cout << "c" << setw(10) << " " << "probing " ;
		    cout.flush();
		  }
		  /*********************************************/

		  Probedvo probe_hrst(updateLW, updateIP, p_heuristic);
		  cp_solver->add( probe_hrst );
                  if(p_valheur == 2)
                    cp_solver->setRandomValueOrdering();
                  else if (p_valheur == 1)
                    cp_solver->setAntiLex();
                  else if (p_valheur == 0)
                    cp_solver->setLex();

                  probe_start = getRunTime();
		  result = cp_solver->random_probe(p_iteration, p_limit);
                  probe_end = getRunTime();

                   if((result == LIMITOUT) || (result == UNKNOWN)){ // DG Added for counting num times prob solved in probing
                     probe_sol = 0;
                     }
                   else{
                     probe_sol = 1;
                     }

                  probe_nodes = cp_solver->NODES; // DG Added, ask Emmanuel about time changes
                  probe_time = (probe_end - probe_start);

		  
		  /*********************************************/
		  if( verbose ) {
		    time_now = getRunTime(); 
		    cout << setw(6)
			 << setprecision(5) << (time_now - time_last) << endl;
		    time_last = time_now;

		    if( use_lds && verbose > 1 ) {
		      cout << "c Impact learned during probing: " << endl;
		      for(int i=0; i<cp_solver->learners.size; ++i)
			{
			  cp_solver->learners[i]->print( cout );
			  cout << "c " << endl;
			}
		    }

		  }
		  /*********************************************/
		}
	      
	      /********************************************
	       * Solving
	       ********************************************/

	      if( use_lds ) 
		{
		  
		  if( verbose ) {
		    cout << "c" << setw(10) << " " << "impact reinforcement search " ;
		    cout.flush();
		  }
		  
		  Impact h;
		  cp_solver->add( h );
		  cp_solver->setTimeLimit( use_lds+getRunTime() );
		  // 	  result = cp_solver->solve_and_restart(
		  // 				       restart_policy,
		  // 				       restart_base,
		  // 				       restart_factor
		  // 				       );
		  result = cp_solver->solve();
		  
		  if( verbose ) {
		    time_now = getRunTime(); 
		    cout << setw(6)
			 << setprecision(5) << (time_now - time_last) << endl;
		    time_last = time_now;

		    if( use_lds && verbose > 1 ) {
		      cout << "c Impact learned during reinforcement phase: " << endl;
		      for(int i=0; i<cp_solver->learners.size; ++i)
			{
			  cp_solver->learners[i]->print( cout );
			  cout << "c " << endl;
			}
		    }
		  }
		  
		  if( verbose ) {
		    cout << "c" << setw(10) << " " << "starting lds " << endl;
		  }	  
		  
		  if( result == LIMITOUT ) {
		    cp_solver->setTimeLimit( 0 );
		    cp_solver->status = UNSAT;
		    ImpactOverWDeg h;
		    cp_solver->add( h );
		    result = cp_solver->ldSolve();
		  }
		  
		  if( verbose ) {
		    cout << "c" << setw(10) << " " << "lds time " ;
		    time_now = getRunTime(); 
		    cout << setw(6)
			 << setprecision(5) << (time_now - time_last) << endl;
		    cout.flush();
		    time_last = time_now;
		  }
		}
	      else if( result == UNKNOWN ) 
		{
		  addHeuristic( *cp_solver );

		  /*********************************************/
		  if( verbose ) {
		    if( restart_policy < 0 ) 
		      cout << "c" << setw(10) << " " << "final run ";
		    else
		      cout << "c" << setw(10) << " " << "starting restarts " << endl;
		    cout.flush();
		  }
		  /*********************************************/

                  if(valheur == 2)
                    cp_solver->setRandomValueOrdering();
                  else if (valheur == 1)
                    cp_solver->setAntiLex();
                  else if (valheur == 0)
                    cp_solver->setLex();

		  if( restart_policy < 0 ) 
		    {
		      result = cp_solver->solve();
		    } 
		  else 
		    {
		      result = cp_solver->solve_and_restart(
							   restart_policy,
							   restart_base,
							   restart_factor
							   );
		    }

		  /*********************************************/
		  if( verbose ) {
		    if( restart_policy >= 0 )		      
		      {
			cout << "c" << setw(10) << " " << "restart time " ;
		      }
		    time_now = getRunTime(); 
		    cout << setw(6)
			 << setprecision(5) << (time_now - time_last) << endl;
		    cout.flush();
		    time_last = time_now;		  
		  }
		  /*********************************************/
		}
	    }
	  /*********************************************/
	  if( verbose ) {
	    
#ifdef _LINUX
	    memory_usage = (getMemory() - memory_usage);
#endif
	    
	    time_now = getRunTime(); 
	    cout << "c Solving time: " << setw(6)
		 << setprecision(5) << (time_now - time_solving) 

#ifdef _LINUX
		 << "  (+" << (memory_usage >> 8) << " Mb)"  
#endif
	      
		 << "\nc\n";
	
	  }
	  /*********************************************/
	}
       
    }

    /*!
     * signal the end of parsing
     */
    virtual void endInstance()
    {


      if( verbose > 1 ) {
	cout << "c\n"
	     << "c" << setw(10) << numEXTCons << " extensional constraints" << endl
	     << "c" << setw(10) << numGACPred << " predicates (GAC)" << endl
	     << "c" << setw(10) << numDECPred << " predicates (decomposition)" << endl;
      }
      if( verbose ) {
	time_last = getRunTime();
	time_now = (time_last - time_start);
	cout << "c Parsing time: " << setw(6) 
	     << setprecision(5) << time_now 
	     << "\nc\nc Building model" << endl;
      }

      //solve();
    }

    int solve() {
      /********************************************
       *
       * Build the model
       *
       ********************************************/

      //Solver s;

      //s.prebuild( model );
      int result = UNKNOWN;
      if( model.preBuild() != UNSAT )
	{
	  if( revert_sat && model.isSatCompatible() ) {

	    std::cout << "SOLVE WITH SAT" << std::endl;

	    solveWithSAT( result );
	  } else {
	    solveWithCP( result );
	  }
	}
      else
	{
	  result = UNSAT;
	  if( ext_feature ) 
	    {
	      features_vec = new double[36];
	      std::fill(features_vec, features_vec+36, -1);
	    }
	}

      return result;
    }


    void print_outcome(int result) {

      /********************************************
       *
       * Print result
       *
       ********************************************/
      
      bool outputinfile=( model_name != "default" );
      ofstream solfile;
      if( outputinfile) {
	string solname( "sol_" + model_name );
	solfile.open( solname.c_str(), ios_base::out );
      }
      
      int k, finalvalue;
      

      if(cp_solver) {
	if( result == SAT ) {
	  if(outputinfile) solfile << "SAT" << endl;
	  cout << "s SATISFIABLE" << endl << "v";
	  for(k=0; k<X.size(); ++k) {
	    finalvalue = X[k].value();
	    if( domainMin[D[k]] != domainMax[D[k]] &&
		finalvalue != NOVAL ) {
	      solution[k] = finalvalue;
	    }
	    cout << " " << solution[k] ;
	    if(outputinfile) solfile << solution[k] << " ";
	  }
	  //if(outputinfile) solfile << endl;
	  cout << endl;
	  // 	cout << endl;
	  // 	for(k=0; k<X.size(); ++k) 
	  // 	  cout << " " << X[k].var_ptr_->getVariable()->id;
	  // 	cout << endl;
	  
	} else if ( result == UNSAT ) {
	  cout << "s UNSATISFIABLE" << endl;
	  if(outputinfile) solfile << "UNSAT" << endl;
	} else {
	  cout << "s UNKNOWN" << endl;
	  if(outputinfile) solfile << "UNKNOWN" << endl;
	}
	
	cout << "d CHECKS " << cp_solver->PROPAGS << endl
	     << "d ASSIGNMENTS " << cp_solver->NODES << endl
	     << "d NODES " << cp_solver->NODES 
	     << " BACKTRACKS " << cp_solver->BACKTRACKS
	     << " FAILURES " << cp_solver->FAILURES
	     << " MISC " << cp_solver->MISC
	     << " RUNTIME " << cp_solver->ENDTIME
	     << " TOTALTIME " << cp_solver->TOTTIME
	     << " NODES/s " << (cp_solver->ENDTIME ? (double)(cp_solver->NODES)/cp_solver->ENDTIME : 100*cp_solver->NODES) 
	     << " CHECKS/s " << (cp_solver->ENDTIME ? (double)(cp_solver->PROPAGS)/cp_solver->ENDTIME : 100*cp_solver->CHECKS); 
         if(probing){
           cout << " ProbingSol " << probe_sol // DG Added next 3 for probing stats
           << " ProbeNodes " << probe_nodes
           << " ProbeTime " << probe_time;
           }
       cout << endl;
	
	if(outputinfile) {  
	  cp_solver->printStatistics( solfile, (OUTCOME |
						TOTALTIME |
						BTS |
						CKS |
						NDS) );
	  solfile << endl;
	}
      }
      else if(sat_solver) {

	if( result == SAT ) {
	  cout << "s SATISFIABLE" << endl << "v";
	  for(k=0; k<X.size(); ++k) {
	    if( domainMin[D[k]] == domainMax[D[k]] ) {
	      cout << " " << domainMin[D[k]];
	    } else cout << " " << ( sat_solver->polarity[k+1] == k+1 ) ;
	  }
	  cout << endl;
	} else if ( result == UNSAT ) {
	  cout << "s UNSATISFIABLE" << endl;
	} else {
	  cout << "s UNKNOWN" << endl;
	}
	double time = (getRunTime() - sat_solver->stats.start_time);
	cout << "d CHECKS " << sat_solver->stats.unit_props << endl
	     << "d ASSIGNMENTS " << sat_solver->stats.nodes << endl
	     << "d NODES " << sat_solver->stats.nodes 
	     << " BACKTRACKS " << sat_solver->stats.conflicts
	     << " FAILURES " << sat_solver->stats.conflicts
	     << " RUNTIME " << time
	     << " MISC " << 0
	     << " TOTALTIME " << time
	     << " NODES/s " << (long unsigned int)(time > 0.001 ? (double)(sat_solver->stats.nodes)/time : 100*sat_solver->stats.nodes) 
	     << " CHECKS/s " << (long unsigned int)(time > 0.001 ? (double)(sat_solver->stats.unit_props)/time : 100*sat_solver->stats.unit_props) << endl ;
      }
      
      //       if(result == UNSAT) {
	
// 	if( verbose ) {
// 	  cout << "c Found inconsistent while building" << endl;
// 	}
	
// 	cout << "s UNSATISFIABLE" << endl;
// 	cout << "d CHECKS " << 0 << endl
// 	     << "d ASSIGNMENTS " << 0 << endl
// 	     << "d NODES " << 0 
// 	     << " BACKTRACKS " << 0
// 	     << " FAILURES " << 0
// 	     << " RUNTIME " << 0
// 	     << " MISC " << 0
// 	     << " TOTALTIME " << 0
// 	     << " NODES/s " << 0
// 	     << " CHECKS/s " << 0 << endl;

//     } else {

//       }
      
    }

    void cleanup() {
      if( numPredicates ) {
	for(int i=0; i<numPredicates; ++i)
	  delete predicates[i];
	delete [] predicates;
      }
      delete [] solution;
      delete [] params;
      delete [] D;
      delete [] domainSize;
      delete [] domainMax;
      delete [] domainMin;
      delete [] features_vec;
      int i=trash.size;
      while( i-- )
	delete trash[i];
      delete cp_solver;
      delete sat_solver;
    }


  };

}; // namespace

#endif

