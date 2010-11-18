

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

#include <cstring>
#include <signal.h>

#include <mistral_sol.h>
#include <mistral_dvo.h>

using namespace Mistral;
using namespace std;

/**********************************************
I/   compute length
II/  build the model
III/ create the array 'sequence'
**********************************************/


// prebuild :-> basic initialisation, invoque prebuild on the model

// build :->  n is an upperbound on level max, invoque build
// on the model and initialise sequence

// initSearch( X ) :-> put the searchable variables in X in front
// of the array variables and returns the number of searched variables.

// initsearch( n, type ) :-> put all variables of the given type
// in front of the array 'variables+n' and the rest at the back.
// return the number of variables whose size is under a given limit


double RPTIMELIMIT=60.0;

Solver* mistral_solver;
static void Mistral_SIGINT_handler(int signum) {
  cout <<"\nc ****************** INTERRUPTED ******************\n";
  //mistral_solver->stopWatch();
  mistral_solver->closeSearch();
  mistral_solver->printStatistics();
  cout <<"c ****************** INTERRUPTED ******************\n";
  exit(1);
  //throw abort_search();
}



/*********************************************
 * Solver
 **********************************************/

void Solver::randomizeSequence()
{

  //cout << "c shuffle the variables" << endl ;

  int i, j, n=(empty-future);
  VariableInt *aux;
  for(i=0; i<n-1; ++i)
    {
      j = i + (randint(n-i));
      aux = future[j];
      future[j] = future[i];
      future[i] = aux;
      future[i]->seqIdx = &future[i];
      future[j]->seqIdx = &future[j];
    }

  //   for(i=0; i<n; ++i) {
  //     future[i]->printshort( cout );
  //     cout << " ";
  //   }
  //  cout << endl << endl;
}

void Solver::reorderSequence()
{
  for(int i=0; i<length; ++i)
    {
      sequence[i] = variables[i];
      sequence[i]->seqIdx = &sequence[i];
    }
  for(int i=length; i<numvars; ++i)
    {
      sequence[i] = variables[i];
      sequence[i]->seqIdx = &sequence[i];
    }
}

Weighter* Solver::setLearner( int wtype )
{

  int i=learners.size;
  if( wtype & Weighter::WLD )
    while( i-- )
      if( learners[i]->getType() == Weighter::WDG ) {
	delete learners[i];
	learners.erase( i );
      }

  i=learners.size;
  if( wtype & Weighter::WDG )
    while( i-- )
      if( learners[i]->getType() == Weighter::WLD ) {
	delete learners[i];
	learners.erase( i );
      }

  i=learners.size;
  while( i-- )
    if( learners[i]->getType() & wtype )
      return learners[i];
  if( (wtype & Weighter::WLD) == Weighter::WLD ) {
    learners.push( new WeighterLevelDegree( this ) );
  }
  else if( (wtype & Weighter::WDG) == Weighter::WDG ) {
    learners.push( new WeighterDegree( this ) );
  }
  if( wtype & Weighter::IPT) {
    learners.push( new WeighterImpact( this ) );
  }
  if( wtype & Weighter::SAC) {
    learners.push( new WeighterSAC( this ) );
  }
  if( wtype & Weighter::ISAC) {
    learners.push( new WeighterISAC( this ) );
  }
  if( wtype & Weighter::RNGD) {
    learners.push( new WeighterRestartNogood( this ) );
  }
  if( wtype & Weighter::RGNGD) {
    learners.push( new WeighterRestartGenNogood( this ) );
  }

  return learners.back();
}


void Solver::init()
{
  usrand( 11041979 );
  randomizedRestart=false;
  boundOptimised=true;
  domainSplitting=false;
  optimisation = true;
  verbosity = 0;
  fail_increment = 0;
  status = UNKNOWN;
  FIND_ALL = 0;
  DISCREPANCY = -1;
  SOLUTIONS = NODES = BACKTRACKS = CHECKS = FAILURES = PROPAGS = FILTERINGS = RESTORES = MISC = 0;
  INITTIME = getRunTime();
  ENDTIME = TOTTIME = SOLTIME = 0;


  BTSLIMIT  =  0;
  NDSLIMIT  =  0;
  TIMELIMIT = .0;
  FAILLIMIT =  0;
  LIMIT = 0;

  level = 0;
  init_level = 1;

  unaryCons_size = 0;
  solution = NULL;
  max_solution = NULL;
  sequence = NULL;
  heuristic = NULL;
  goal = NULL;
  function = NULL;
  sat = NULL;
}

bool Solver::is_built()
{
  return(sequence != NULL);
}

void Solver::prebuild(CSP& c)
{
  if( !c.preBuild() ) {
    status = UNSAT;
  }
}

void Solver::build(CSP& c, const int n)
{
  length = std::min(n, c.numVars());

  c.build( this );
  numvars    =  variables.size;
  sequence = new VariableInt*[numvars];
  decision.init(0,numvars+1);
  branching_decision.init(0,numvars+1);
}

void Solver::initSearch( BuildObject **x, const int l, const int maxsize )
{
  int i=0, j=0, n=l, idx;

//   std::cout << "INIT SEARCH" << std::endl;

//   for(i=0; i<n; ++i) {
//     x[i]->print_python();
//   }
  //i=0;
  //std::cout << std::endl << i << " " << n << std::endl;


  // remove from x the variables that are "non-searchable"
  VariableInt *temp;
  BuildObject *bvar;
  while( i<n ) {
    bvar = x[i];

    //     bvar->print( cout );
    //     //     cout << " / ";
    //     //     if( bvar->getVariable() )
    //     //       bvar->getVariable()->print( cout );
    //     //     cout << " ";
    //     //     cout << bvar;

    //     cout << "(" << bvar << " / " << (bvar->getBuildObject()) << ") ";

    //     cout << " " << (bvar->isSearchable()) << " "
    // 	 << (bvar->getVariable()) << endl;



    if(bvar != bvar->getBuildObject())
      bvar = bvar->getBuildObject();

  //   bvar->print(std::cout);
//     std::cout <<  " " << (bvar->isSearchable()) << " " << bvar->getVariable() << std::endl;
    
    if( bvar->isSearchable() &&
	((temp = bvar->getVariable()) != NULL) ) {

      if(temp->getType() == VariableInt::VIRTUAL) {
// 	std::cout << "     swap virtual (v" << temp->id << " = x" ;
 	temp = ((VariableVirtual*)temp)->reference;
// 	std::cout << temp->id << ") with " << j << std::endl;
//       } else {
// 	std::cout << "not virtual " << (temp->getType()) << " =/= " << (VariableInt::VIRTUAL) << std::endl;
       }




      idx = temp->id;
      if( idx != j ) {
	variables[idx] = variables[j];
	variables[j] = temp;

	variables[j]->id = j;
	variables[idx]->id = idx;
      }


//        bvar->print(std::cout);
//        std::cout << std::endl << std::endl;

      ++j;
    }
    ++i;
  }
  n=j;


 //  std::cout << j << std::endl;
//      for(i=0; i<j; ++i)
//        {
//          variables[i]->print( std::cout );
//   //       cout << endl;
//        }
//      std::cout << std::endl;

//   //   //  exit(0);

//   //  cout << "INIT SEARCH " << n << endl;

  initSearch( n, maxsize, false );
}



void Solver::initSearch( VarArray& X, const int maxsize )
{
  int i=0, j=0, n=X.size(), idx;

  // remove from x the variables that are "non-searchable"
  VariableInt *temp;
  BuildObject *bvar;

  //cout << variables.size << endl;

  while( i<n ) {
#ifdef _NARRAY
    bvar = X.array_[i];
#else
    bvar = X.array_[i].var_ptr_;
#endif


    if(bvar != bvar->getBuildObject())
      bvar = bvar->getBuildObject();


    if( bvar->isSearchable() &&
	((temp = bvar->getVariable()) != NULL) ) {
      
      bool is_not_in = true;
      for(int x=0; x<j && is_not_in; ++x)
	is_not_in = (variables[x] != temp);

      if(is_not_in) {
// 	temp->print( cout );
// 	cout << " / ";
// 	cout.flush();
	
	idx = temp->id;
	if( idx != j ) {
	  
// 	  cout << " swap ";
// 	  variables[idx]->print(cout);
	  
// 	  cout << " and ";
// 	  variables[j]->print(cout);
	  
// 	  cout << endl;
	  
	  
	  variables[idx] = variables[j];
	  variables[j] = temp;
	  
	  variables[j]->id = j;
	  variables[idx]->id = idx;
	}

// 	temp->print( cout );
// 	cout << endl;
      

	++j;
      }
    }
    ++i;
  }
  n=j;

  //   for(i=0; i<j; ++i)
  //     {
  //       variables[i]->print( cout );
  //       cout << endl;
  //     }

  //   //  exit(0);

  //cout << "HERE" << endl;

  initSearch( n, maxsize, false );
}

void Solver::initSearch( const int n,
			 const int maxsize,
			 const bool range )
{
  int i, searchable = n, //= ( maxsize < NOVAL ? n : numvars ),
    end = numvars;


  //cout << maxsize << " " << searchable << endl;
  //   exit(0);


  if( maxsize ) {
    VariableInt* aux;

    //cout << searchable << " " << end << endl;
    //exit(0);

    while( searchable<end )
      {
	aux = variables[searchable];
	if( aux->isGround() || aux->domsize() > maxsize || aux->getType() == VariableInt::VIRTUAL ||
	    (!range && aux->getType() == VariableInt::RANGE) )
	  {
	    variables[searchable] = variables[--end];
	    variables[searchable]->id = searchable;
	    variables[end] = aux;
	    aux->id = end;
	  }
	else {
	  ++searchable;
	}
      }
  }

  // then copy the rest of variables into sequence
  std::memcpy(sequence, variables.stack_, numvars*sizeof(VariableInt*));


  empty = (sequence+searchable);
  future = sequence;
  auxilliary = empty;


  //   cout << searchable << endl;
  //   for(VariableInt **iter=future; iter<empty; ++iter) {
  //     (*iter)->print( cout );
  //     cout << endl;
  //   }
  //   cout << future << " " << empty << endl;
  //   exit(0);



  length   = (empty-sequence);


  gacvarstack.initList(numvars, variables.stack_);
  gacconstack.initStack(constraints.size);


  solution     = new int[numvars];
  max_solution = new int[numvars];
  std::fill( solution,     solution    +numvars, NOVAL );
  std::fill( max_solution, max_solution+numvars, NOVAL );


  decision.push( NULL );
  SimpleUnaryConstraint no_decision;
  branching_decision.push( no_decision );
  lvl_.push(0);
  past.push( future );
  auxv.push( auxilliary );


  for(i=0; i<length; ++i) {
    sequence[i]->seqBeg = &future;
    sequence[i]->seqIdx = sequence+i;
  }

  for(i=0; i<(numvars-length); ++i) {
    auxilliary[i]->seqBeg = &auxilliary;
    auxilliary[i]->seqIdx = auxilliary+i;
  }
}

Solver::Solver(CSP& c, VarArray& x, VariableOrdering& v)
{
  init();
  init(c, x, v);
}

void Solver::init(CSP& c, VarArray& x, VariableOrdering& v)
{
  prebuild( c );
  if( status == UNKNOWN ) {
    build( c, x.size() );
    initSearch( x );
    heuristic = v.extract( this );
  }
}

Solver::Solver(CSP& c, VarArray& x)
{
  init();
  init(c, x);
}

void Solver::init(CSP& c, VarArray& x)
{
  prebuild( c );
  if( status == UNKNOWN ) {
    build( c, x.size() );
    initSearch( x );
  }
}

Solver::Solver(CSP& c, BuildObject **x, const int l)
{
  init();
  init(c, x, l);
}

void Solver::init(CSP& c, BuildObject **x, const int l)
{
  prebuild( c );
  if( status == UNKNOWN ) {
    build( c, l );
    initSearch( x, l, 0 );
  }
}

Solver::Solver(CSP& c, VariableOrdering& v)
{
  init();
  init(c, v);
}

void Solver::init(CSP& c, VariableOrdering& v)
{
  prebuild( c );
  if( status == UNKNOWN ) {
    build( c );
    initSearch();
    heuristic = v.extract( this );
  }
}

Solver::Solver(CSP& c)
{
  init();
  init(c);
}

void Solver::init(CSP& c)
{
  prebuild( c );
  if( status == UNKNOWN ) {
    build( c );
    initSearch();
  }

  if(sat) sat->initialise();
}

Solver::Solver()
{
  init();
}

void Solver::add(VariableOrdering& v)
{
  if(status == UNKNOWN) {
    DVO *old_heuristic = heuristic;
    heuristic = v.extract( this );
    if( old_heuristic ) delete old_heuristic;
  }
}

void Solver::add(BranchingStrategy& b)
{
  b.extract();
}


Solver::~Solver()
{

  //std::cout  << "*** DELETE SOLVER ***" << endl;

  delete heuristic;
  delete goal;
  delete function;

  delete [] solution;
  delete [] max_solution;
  delete [] sequence;

  int i;
  for(i=0; i<constraints.size; ++i)
    delete constraints[i];

  for(i=0; i<variables.size; ++i)
    delete variables[i];


  i=learners.size;
  while( i-- )
    {
      delete learners[i];
    }
}

void Solver::upOneLevel() {
  past.push( future );
  auxv.push( auxilliary  );
  lvl_.push( store.size );
  decision.push( NULL );
  SimpleUnaryConstraint no_decision;
  branching_decision.push( no_decision );
  ++level;
}

int Solver::presolve()
{

#ifdef _DEBUGAC
      if(verbosity > 2) {
  cout << "\tPRESOLVE" << endl;
   for(int k=0; k<numvars; ++k) {
     variables[k]->print( cout );
     cout << endl;
   }
   cout << endl;
      }
#endif

#ifdef _DEBUGPROPAG
  //std::ostringstream o_propag;
  std::string prop_str;
  int size_before = 0;
#endif

  //cout << "\tPRESOLVE init:" << init_level << " level " << level << endl;

  if( level < init_level )
    {

      // create nodes to buffer the first round of pruning
      // (in case we want to resume to the original formulation)
      while(level < init_level) {
	past.push( future );
	auxv.push( auxilliary  );
	lvl_.push( store.size );
	decision.push( NULL );
	SimpleUnaryConstraint no_decision;
	branching_decision.push( no_decision );

	++level;
      }

      int i = learners.size;
      while( i-- )
	learners[i]->notify_init_level(init_level);

      if(sat)
	sat->initialise();

      if( heuristic == NULL )
	//heuristic = new DVONoOrder( this );
	heuristic = new GenericDVO<VarSelectorDomain>( this );

      i = constraints.size;
      while( i-- ) {

#ifdef _DEBUGPROPAG
	std::ostringstream o_propag;
	size_before = 0;
      if(verbosity > 2) {
	o_propag << i << " PROPAGATE \n" ;
	constraints[i]->print( o_propag );
	o_propag << endl;
	for(int a=0; a<constraints[i]->arity; ++a) {
	  size_before += constraints[i]->scope[a]->domsize();
	  constraints[i]->scope[a]->print( o_propag );
	  o_propag << " ";
	}
	o_propag << endl;
	prop_str = o_propag.str();
      }
#endif

	// 	bool is42 = false;
	// 	for(int x=0; !is42 && x<constraints[i]->arity; ++x)
	// 	  is42 |= (constraints[i]->scope[x]->id == 41);
	// 	if(is42) {
	// 	  std::ostringstream o_propag;
	// 	  size_before = 0;

	// 	  o_propag << i << " PROPAGATE \n" ;
	// 	  constraints[i]->print( o_propag );
	// 	  o_propag << endl;
	// 	  for(int a=0; a<constraints[i]->arity; ++a) {
	// 	    size_before += constraints[i]->scope[a]->domsize();
	// 	    constraints[i]->scope[a]->print( o_propag );
	// 	    o_propag << " ";
	// 	  }
	// 	  o_propag << endl;
	// 	  prop_str = o_propag.str();
	// 	}

	++PROPAGS;
	if( !(constraints[i]->propagate( )) ) {

	  // 	  if(is42) {
	  // 	    cout << prop_str;
	  // 	    for(int a=0; a<constraints[i]->arity; ++a) {
	  // 	      constraints[i]->scope[a]->print( cout );
	  // 	      cout << " ";
	  // 	    }
	  // 	    cout << endl;
	  // 	    constraints[i]->print( cout );
	  // 	    cout << " FAIL" << endl << endl;
	  // 	  }

#ifdef _DEBUGPROPAG
      if(verbosity > 2) {
	  cout << prop_str;
	  for(int a=0; a<constraints[i]->arity; ++a) {
	    constraints[i]->scope[a]->print( cout );
	    cout << " ";
	  }
	  cout << endl;
	  constraints[i]->print( cout );
	  cout << " FAIL" << endl << endl;
      }
#endif


	  status = UNSAT;
	  return status;
	}

	// 	if(is42) {
	// 	  for(int a=0; a<constraints[i]->arity; ++a) {
	// 	    size_before -= constraints[i]->scope[a]->domsize();
	// 	  }

	// 	  if(size_before) {
	// 	    cout << prop_str;
	// 	    for(int a=0; a<constraints[i]->arity; ++a) {
	// 	      constraints[i]->scope[a]->print( cout );
	// 	      cout << " ";
	// 	    }
	// 	    cout << endl;
	// 	    constraints[i]->print( cout );
	// 	    cout << " OK" << endl << endl;
	// 	  }
	// 	}

#ifdef _DEBUGPROPAG
      if(verbosity > 2) {
	for(int a=0; a<constraints[i]->arity; ++a) {
	  size_before -= constraints[i]->scope[a]->domsize();
	}

	if(size_before) {
	  cout << prop_str;
	  for(int a=0; a<constraints[i]->arity; ++a) {
	    constraints[i]->scope[a]->print( cout );
	    cout << " ";
	  }
	  cout << endl;
	  constraints[i]->print( cout );
	  cout << " OK" << endl << endl;
	}
      }
#endif

      }

      i = unaryCons.size;
      while( i-- ) {
	++PROPAGS;


#ifdef _DEBUGPROPAG
      if(verbosity > 2) {
	cout << i << " PROPAGATE \n" ;
	unaryCons[i]->print( cout );
	cout << endl;
      }
#endif

	if( !(unaryCons[i]->propagate( )) ) {


#ifdef _DEBUGPROPAG
      if(verbosity > 2) {
	  unaryCons[i]->print( cout );
	  cout << " FAIL" << endl;
      }
#endif

	  status = UNSAT;
	  return status;
	} else {

#ifdef _DEBUGPROPAG
      if(verbosity > 2) {
	  unaryCons[i]->print( cout );
	  cout << " OK" << endl;
      }
#endif

	}
      }
      unaryCons_size = unaryCons.size;
      unaryCons.size = 0;

      i = numvars;
      while( i-- )
	if( variables[i]->isGround() )
	  //triggerEvent( variables[i], Constraint::VALUETRIGGER );
	  triggerEvent( variables[i]->id, Constraint::VALUETRIGGER );

      //std::cout << "START FILTER IN PRESOLVE" << std::endl;
      if( !filtering() ) {
	status = UNSAT;
	return status;
      }
      //std::cout << "END FILTER IN PRESOLVE" << std::endl;


      i = length;
      while( i-- )
	if( !variables[i]->isGround() && !variables[i]->branch ) {
	  if( domainSplitting && variables[i]->domsize() > 3 && variables[i]->getType() != VariableInt::LIST ) {
	    variables[i]->branch = new ValSelectorSplit( variables[i] );
	  } else if( variables[i]->getType() == VariableInt::BOOL ) {
	    //variables[i]->branch = new ValSelectorRandMinMax( variables[i] );
	    variables[i]->branch = new ValSelectorMin( variables[i] );
	  } else if( variables[i]->getType() == VariableInt::RANGE ) {
	    //variables[i]->branch = new ValSelectorRandMinMax( variables[i] );
	    variables[i]->branch = new ValSelectorMin( variables[i] );
	  } else if( variables[i]->getType() == VariableInt::BIT ) {
	    //variables[i]->branch = new ValSelectorRandMinMax( variables[i] );
	    variables[i]->branch = new ValSelectorMin( variables[i] );
	  } else if( variables[i]->getType() == VariableInt::LIST ) {
	    //variables[i]->branch = new ValSelectorFirst( variables[i] );
	    variables[i]->branch = new ValSelectorMin( variables[i] );
	    //variables[i]->branch = new ValSelectorMax( variables[i] );
	    //variables[i]->branch = new ValSelectorRand( variables[i] );
	    //variables[i]->branch = new ValSelectorRandMinMax( variables[i] );
	  } else {
	    variables[i]->branch = new ValSelectorMin( variables[i] );
	  }
	}
    }

#ifdef _DEBUGAC
      if(verbosity > 2) {
  for(int k=0; k<numvars; ++k) {
    variables[k]->print( cout );
    cout << endl;
  }
  cout << endl;
      }
#endif

  return status;
}

int Solver::solve()
{
  mistral_solver = this;
  signal(SIGINT,Mistral_SIGINT_handler);

  try {
    if( verbosity > 0 )
      cout << "c ===============[ Mistral - solve ]===============" << endl;
    if( randomizedRestart ) randomizeSequence();
    if( init_level < level ) init_level = level;
    if( status == UNKNOWN ) {
      STARTTIME = getRunTime();
      if( presolve() == UNKNOWN ) {
	if(function)
	  function->initialise();
	iterative_dfs( );
      }
      closeSearch();
    }
  } catch (const abort_search& ex) {
    std::cout << "c *        Abort search exception caught!         *" << std::endl
	      << "c *************************************************" << std::endl;
    status = UNKNOWN;
    //exit(1);
  } catch (...) {
    std::cout << "Fatal error, unkown exception" << std::endl;
    exit(1);
  }
  return status;
}


int Solver::pseudoRldsSolve(int deltaK)
{
  mistral_solver = this;
  signal(SIGINT,Mistral_SIGINT_handler);

  try {
    if( verbosity > 0 )
      cout << "c ===============[ Mistral - pseudoRldsSolve ]===============" << endl;
    if( status == UNKNOWN ) {
      int threshold = deltaK / 2;

      // need to set up reversible structures before search
      ReversibleNum<int> numDiscrepancies;
      this->binds(numDiscrepancies);
      numDiscrepancies.setValue(0);

      STARTTIME = getRunTime();
      while (status == UNKNOWN) {

	// the presolve has to be done before each iteration
	if( presolve() == UNKNOWN ) {
	  pseudoRldsProbe(threshold, numDiscrepancies);
	}
	if(status == UNKNOWN) { // the probe ended because of the discrepancies
	  threshold += deltaK;
	  // need to reset the problem to its initial state
	  reset_trail(true);
	}
      }
    }
    closeSearch();
  } catch (const abort_search& ex) {
    std::cout << "c *        Abort search exception caught!         *" << std::endl
	      << "c *************************************************" << std::endl;
    status = UNKNOWN;
  } catch (...) {
    std::cout << "Fatal error, unkown exception" << std::endl;
    exit(1);
  }
  return status;
}

int Solver::pseudoRldsSolve_level(int deltaK)
{
  mistral_solver = this;
  signal(SIGINT,Mistral_SIGINT_handler);

  try {
    if( verbosity > 0 )
      cout << "c ===============[ Mistral - pseudoRldsSolve_level ]===============" << endl;
    if( status == UNKNOWN ) {
      unsigned int threshold = deltaK / 2;

      // need to set up reversible structures before search
      ReversibleIntList discrepancyLevels;
      this->binds(discrepancyLevels);
      // size of discrepency_level is the sum of domains from all search variables
      int size = 0;
      for (int variableId=0; variableId<length; variableId++) {
	size += sequence[variableId]->domsize();
      }
      discrepancyLevels.setValue(0, size, 0);

      STARTTIME = getRunTime();
      while (status == UNKNOWN) {

	// the presolve has to be done before each iteration
	if( presolve() == UNKNOWN ) {
	  pseudoRldsProbe_level(threshold, discrepancyLevels);
	}
	if(status == UNKNOWN) { // the probe ended because of the discrepancies
	  threshold += deltaK;
	  // need to reset the problem to its initial state
	  reset_trail(true);
	}
      }
    }
    closeSearch();
  } catch (const abort_search& ex) {
    std::cout << "c *        Abort search exception caught!         *" << std::endl
	      << "c *************************************************" << std::endl;
    status = UNKNOWN;
  } catch (...) {
    std::cout << "Fatal error, unkown exception" << std::endl;
    exit(1);
  }
  return status;
}

int Solver::pseudoRldsSolve_variable(int deltaK)
{
  mistral_solver = this;
  signal(SIGINT,Mistral_SIGINT_handler);

  try {
    if( verbosity > 0 )
      cout << "c ===============[ Mistral - pseudoRldsSolve_variable ]===============" << endl;
    if( status == UNKNOWN ) {
      unsigned int threshold = deltaK / 2;

      // need to set up reversible structures before search
      ReversibleIntList discrepency_variables;
      this->binds(discrepency_variables);
      discrepency_variables.setValue(0, variables.size-1, 0);

      STARTTIME = getRunTime();
      while (status == UNKNOWN) {
	// the presolve has to be done before each iteration
	if( presolve() == UNKNOWN ) {
	  pseudoRldsProbe_variable(threshold, discrepency_variables);
	}
	if(status == UNKNOWN) { // the probe ended because of the discrepancies
	  threshold += deltaK;
	  // need to reset the problem to its initial state
	  reset_trail(true);
	}
      }
    }
    closeSearch();
  } catch (const abort_search& ex) {
    std::cout << "c *        Abort search exception caught!         *" << std::endl
	      << "c *************************************************" << std::endl;
    status = UNKNOWN;
  } catch (...) {
    std::cout << "Fatal error, unkown exception" << std::endl;
    exit(1);
  }
  return status;
}

int Solver::ldsSolve(int deltaK)
{
  mistral_solver = this;
  signal(SIGINT,Mistral_SIGINT_handler);

  try {
    if( verbosity > 0 )
      cout << "c ===============[ Mistral - ldsSolve ]===============" << endl;
    if( status == UNKNOWN ) {
      int threshold = deltaK / 2;

      // need to set up reversible structures before search
      ReversibleNum<int> numDiscrepancies;
      this->binds(numDiscrepancies);
      numDiscrepancies.setValue(0);

      STARTTIME = getRunTime();
      while (status == UNKNOWN) {
	// the presolve has to be done before each iteration
	if( presolve() == UNKNOWN ) {
	  ldsProbe(threshold, numDiscrepancies);
	}
	if(status == UNKNOWN) { // the probe ended because of the discrepancies
	  threshold += deltaK;
	  // need to reset the problem to its initial state
	  reset_trail(true);
	}
      }
    }
    closeSearch();
  } catch (const abort_search& ex) {
    std::cout << "c *        Abort search exception caught!         *" << std::endl
	      << "c *************************************************" << std::endl;
    status = UNKNOWN;
  } catch (...) {
    std::cout << "Fatal error, unknown exception" << std::endl;
    exit(1);
  }
  return status;
}

int Solver::ldsStackSolve(int deltaK)
{
  mistral_solver = this;
  signal(SIGINT,Mistral_SIGINT_handler);

  try {
    if( verbosity > 0 )
      cout << "c ===============[ Mistral - ldsStackSolve ]===============" << endl;
    if( status == UNKNOWN ) {
      int threshold = deltaK / 2;

      // need to set up reversible structures before search
      ReversibleNum<int> numDiscrepancies;
      this->binds(numDiscrepancies);
      numDiscrepancies.setValue(0);

      STARTTIME = getRunTime();
      while (status == UNKNOWN) {
	// the presolve has to be done before each iteration
	if( presolve() == UNKNOWN ) {
	  ldsStackProbe(threshold, numDiscrepancies);
	}
	if(status == UNKNOWN) { // the probe ended because of the discrepancies
	  threshold += deltaK;
	  // need to reset the problem to its initial state
	  reset_trail(true);
	}
      }
    }
    closeSearch();
  } catch (const abort_search& ex) {
    std::cout << "c *        Abort search exception caught!         *" << std::endl
	      << "c *************************************************" << std::endl;
    status = UNKNOWN;
  } catch (...) {
    std::cout << "Fatal error, unknown exception" << std::endl;
    exit(1);
  }
  return status;
}

int Solver::ldsDeltaSolve(int deltaK)
{
  mistral_solver = this;
  signal(SIGINT,Mistral_SIGINT_handler);

  try {
    if( verbosity > 0 )
      cout << "c ===============[ Mistral - ldsDeltaSolve ]===============" << endl;
    if( status == UNKNOWN ) {
      int threshold = deltaK / 2;

      // need to set up reversible structures before search
      ReversibleNum<int> numDiscrepancies;
      this->binds(numDiscrepancies);
      numDiscrepancies.setValue(0);

      STARTTIME = getRunTime();
      while (status == UNKNOWN) {
	// the presolve has to be done before each iteration
	if( presolve() == UNKNOWN ) {
	  ldsDeltaProbe(threshold, numDiscrepancies);
	}
	if(status == UNKNOWN) { // the probe ended because of the discrepancies
	  threshold += deltaK;
	  // need to reset the problem to its initial state
	  reset_trail(true);
	}
      }
    }
    closeSearch();
  } catch (const abort_search& ex) {
    std::cout << "c *        Abort search exception caught!         *" << std::endl
	      << "c *************************************************" << std::endl;
    status = UNKNOWN;
  } catch (...) {
    std::cout << "Fatal error, unknown exception" << std::endl;
    exit(1);
  }
  return status;
}

int Solver::ldsSolve_level(int deltaK)
{
  mistral_solver = this;
  signal(SIGINT,Mistral_SIGINT_handler);

  try {
    if( verbosity > 0 )
      cout << "c ===============[ Mistral - ldsSolve_level ]===============" << endl;
    if( status == UNKNOWN ) {
      int threshold = deltaK / 2;

      // need to set up reversible structures before search
      ReversibleIntList discrepancyLevels;
      this->binds(discrepancyLevels);
      // size of discrepency_level is the sum of domains from all search variables
      int size = 1;
      for (int variableId=0; variableId<length; variableId++) {
	//FIXME: ici on fait le produit et non plus la somme -> ca fait bcp!
	//       changer l'algo pour qu'il y ait moins de niveaux??
	//       on pourrait rendre myLevel reversible, et la somme suffirait
	size += sequence[variableId]->domsize();
      }
      discrepancyLevels.setValue(0, size, 0);

      ReversibleNum<int> myLevel;
      this->binds(myLevel);
      myLevel.setValue(1);

      STARTTIME = getRunTime();
      while (status == UNKNOWN) {
	// the presolve has to be done before each iteration
	if( presolve() == UNKNOWN ) {
	  ldsProbe_level(threshold, discrepancyLevels, myLevel);
	}
	if(status == UNKNOWN) { // the probe ended because of the discrepancies
	  threshold += deltaK;
	  // need to reset the problem to its initial state
	  reset_trail(true);
	}
      }
    }
    closeSearch();
  } catch (const abort_search& ex) {
    std::cout << "c *        Abort search exception caught!         *" << std::endl
	      << "c *************************************************" << std::endl;
    status = UNKNOWN;
  } catch (...) {
    std::cout << "Fatal error, unknown exception" << std::endl;
    exit(1);
  }
  return status;
}

int Solver::ldsSolve_variable(int deltaK)
{
  mistral_solver = this;
  signal(SIGINT,Mistral_SIGINT_handler);

  try {
    if( verbosity > 0 )
      cout << "c ===============[ Mistral - ldsSolve_variable ]===============" << endl;
    if( status == UNKNOWN ) {
      int threshold = deltaK / 2;

      // need to set up reversible structures before search
      ReversibleIntList discrepency_variables;
      this->binds(discrepency_variables);
      discrepency_variables.setValue(0, variables.size-1, 0);

      STARTTIME = getRunTime();
      while (status == UNKNOWN) {
	// the presolve has to be done before each iteration
	if( presolve() == UNKNOWN ) {
	  ldsProbe_variable(threshold, discrepency_variables);
	}
	if(status == UNKNOWN) { // the probe ended because of the discrepancies
	  threshold += deltaK;
	  // need to reset the problem to its initial state
	  reset_trail(true);
	}
      }
    }
    closeSearch();
  } catch (const abort_search& ex) {
    std::cout << "c *        Abort search exception caught!         *" << std::endl
	      << "c *************************************************" << std::endl;
    status = UNKNOWN;
  } catch (...) {
    std::cout << "Fatal error, unknown exception" << std::endl;
    exit(1);
  }
  return status;
}

int Solver::ddsSolve(int deltaK)
{
  mistral_solver = this;
  signal(SIGINT,Mistral_SIGINT_handler);

  try {
    if( verbosity > 0 )
      cout << "c ===============[ Mistral - ddsSolve ]===============" << endl;
    if( status == UNKNOWN ) {
      int threshold = deltaK / 2;

      // need to set up reversible structures before search
      ReversibleNum<int> numDiscrepancies;
      this->binds(numDiscrepancies);
      numDiscrepancies.setValue(0);

      STARTTIME = getRunTime();
      while (status == UNKNOWN) {
	// the presolve has to be done before each iteration
	if( presolve() == UNKNOWN ) {
	  ddsProbe(threshold, numDiscrepancies);
	}
	if(status == UNKNOWN) { // the probe ended because of the discrepancies
	  threshold += deltaK;
	  // need to reset the problem to its initial state
	  reset_trail(true);
	}
      }
    }
    closeSearch();
  } catch (const abort_search& ex) {
    std::cout << "c *        Abort search exception caught!         *" << std::endl
	      << "c *************************************************" << std::endl;
    status = UNKNOWN;
  } catch (...) {
    std::cout << "Fatal error, unknown exception" << std::endl;
    exit(1);
  }
  return status;
}

int Solver::ddsSolve_level(int deltaK)
{
  mistral_solver = this;
  signal(SIGINT,Mistral_SIGINT_handler);

  try {
    if( verbosity > 0 )
      cout << "c ===============[ Mistral - ddsSolve_level ]===============" << endl;
    if( status == UNKNOWN ) {
      int threshold = deltaK / 2;

      STARTTIME = getRunTime();
      while (status == UNKNOWN) {
	// the presolve has to be done before each iteration
	if( presolve() == UNKNOWN ) {
	  ddsProbe_level(threshold);
	}
	if(status == UNKNOWN) { // the probe ended because of the discrepancies
	  threshold += deltaK;
	  // need to reset the problem to its initial state
	  reset_trail(true);
	}
      }
    }
    closeSearch();
  } catch (const abort_search& ex) {
    std::cout << "c *        Abort search exception caught!         *" << std::endl
	      << "c *************************************************" << std::endl;
    status = UNKNOWN;
  } catch (...) {
    std::cout << "Fatal error, unknown exception" << std::endl;
    exit(1);
  }
  return status;
}

int Solver::ddsSolve_variable(int deltaK)
{
  mistral_solver = this;
  signal(SIGINT,Mistral_SIGINT_handler);

  try {
    if( verbosity > 0 )
      cout << "c ===============[ Mistral - ddsSolve_variable ]===============" << endl;
    if( status == UNKNOWN ) {
      int threshold = deltaK / 2;

      // need to set up reversible structures before search
      ReversibleIntList discrepency_variables;
      this->binds(discrepency_variables);
      discrepency_variables.setValue(0, variables.size-1, 0);

      STARTTIME = getRunTime();
      while (status == UNKNOWN) {
	// the presolve has to be done before each iteration
	if( presolve() == UNKNOWN ) {
	  ddsProbe_variable(threshold, discrepency_variables);
	}
	if(status == UNKNOWN) { // the probe ended because of the discrepancies
	  threshold += deltaK;
	  // need to reset the problem to its initial state
	  reset_trail(true);
	}
      }
    }
    closeSearch();
  } catch (const abort_search& ex) {
    std::cout << "c *        Abort search exception caught!         *" << std::endl
	      << "c *************************************************" << std::endl;
    status = UNKNOWN;
  } catch (...) {
    std::cout << "Fatal error, unknown exception" << std::endl;
    exit(1);
  }
  return status;
}

int Solver::ldSolve(VarArray& scp, int *ideal, const int step, const int limit)
{
  BuildObject* scope[scp.size()];
  for(int i=0; i<scp.size(); ++i)
    scope[i] = scp[i].var_ptr_;
  return ldSolve(scp.size(), scope, ideal, step, limit);
}

int Solver::ldSolve(const int n, BuildObject** scope, int *ideal, const int step, const int limit)
{
  if( status == UNKNOWN ) {
    STARTTIME = getRunTime();
    if( presolve() == UNKNOWN ) {
      ConstraintHamming *lds = new ConstraintHamming( this, n, scope, ideal );
      //constraints.push( lds );
      status = ldSolve( lds, step, limit );
    }
    closeSearch();
    //stopWatch();
  }
  return status;
}

int Solver::ldSolve(const int *ideal, const int step, const int limit)
{
  if( status == UNKNOWN ) {
    STARTTIME = getRunTime();
    if( presolve() == UNKNOWN ) {
      for(int i=0; i<length; ++i) {
	std::cout << " " << ideal[i] ;
      }
      std::cout << std::endl;

      ConstraintHamming *lds = new ConstraintHamming( this, ideal );
      //constraints.push( lds );
      status = ldSolve( lds, step, limit );
    }
    closeSearch();
    //stopWatch();
  }
  return status;
}
int Solver::ldSolve(const int step, const int limit)
{
  if( status == UNKNOWN ) {
    STARTTIME = getRunTime();
    if( presolve() == UNKNOWN ) {
      ConstraintDistance *lds = new ConstraintDistance( this );
      //constraints.push( lds );
      //int ub = lds->getMaxThreshold() ;
      status = ldSolve( lds, step, limit );
    }
    closeSearch();
    //stopWatch();
  }
  return status;
}

int Solver::ldSolve(ConstraintLDS *lds, const int step, const int limit)
{
  mistral_solver = this;
  signal(SIGINT,Mistral_SIGINT_handler);

  try {
    if( status == UNKNOWN ) {
      int ub = std::min(lds->getMaxThreshold(), limit);

      status = UNSAT;
      init_level = 1;

      ++DISCREPANCY;
      STARTTIME = getRunTime();
      while( status == UNSAT && DISCREPANCY < ub ) {

	++level;
	past.push( future );
	auxv.push( auxilliary  );
	lvl_.push( store.size );
	decision.push( NULL );
	SimpleUnaryConstraint no_decision;
	branching_decision.push( no_decision );


	lds->lb_threshold = DISCREPANCY;
	lds->ub_threshold = DISCREPANCY+step-1;

	if( verbosity > 0 )
	  {
	    ENDTIME = getRunTime() - STARTTIME;
	    cout << "c" << setw(10) << " " << "discrepancy: "
		 << setw(5) << DISCREPANCY << " to "
		 << (DISCREPANCY+step-1) << "  ( " << setw(4) << ENDTIME
		 << " s)" << endl;
	  }


	//decision.push( NULL );
	//lvl_.push(0);
	//past.push( future );
	//auxv.push( auxilliary );


	//recursive_dfs();


#ifdef _DEBUGPROPAG
      if(verbosity > 2) {
	cout << " PROPAGATE \n" ;
	lds->print( cout );
	cout << endl;
	for(int a=0; a<lds->arity; ++a) {
	  lds->scope[a]->print( cout );
	  cout << " ";
	}
	cout << endl;
      }
#endif

	++PROPAGS;
	bool lds_propagation = lds->propagate();

#ifdef _DEBUGPROPAG
      if(verbosity > 2) {
	for(int a=0; a<lds->arity; ++a) {
	  lds->scope[a]->print( cout );
	  cout << " ";
	}
	cout << endl;
	lds->print( cout );
	cout << endl;
	if(lds_propagation)
	  cout << " OK" << endl << endl;
	else
	  cout << " FAIL" << endl << endl;
      }
#endif

	if( !lds_propagation || !filtering() ) {
	  status = UNSAT;
	} else {

	  if(function)
	    function->initialise();

	  status = UNKNOWN;
	  iterative_dfs();

	}

	if( status == UNSAT ) {
	  reset(true);
	}

	DISCREPANCY += step;
      }
    }
  } catch (const abort_search& ex) {
    std::cout << "c *        Abort search exception caught!         *" << std::endl
	      << "c *************************************************" << std::endl;
    status = UNKNOWN;
    //exit(1);
  } catch (...) {
    std::cout << "Fatal error, unkown exception" << std::endl;
    exit(1);
  }
  return status;
}


int binlog( unsigned long int x ) {
  if( !x ) return NOVAL;
  if( x == 1 ) return 0;
  int exponent = 0;
  while( (x >> (++exponent)) > 1 );
  return exponent;
}
unsigned long int luby_seq( unsigned long int iteration ) {
  int thelog = binlog( iteration );
  if( (int)iteration == (1 << (thelog + 1))-1 )
    return (1 << thelog);
  return luby_seq( iteration - (1 << thelog) + 1 );
}

bool Solver::restart(const double decay, const int reinit) {

//   std::cout << "X[0,0] just before a restart: ";
//   variables[0]->print(std::cout);
//   std::cout << std::endl;

#ifdef _DEBUGSEARCH
      if(verbosity > 2) {
	std::cout << "c";
	for(int k=0; k<=level; ++k) std::cout << " ";
	std::cout << "restart!" << std::endl;
      }
#endif



  if( randomizedRestart ) randomizeSequence();
  int i=learners.size;

  if(function)
    function->initialise();

  while( i-- )
    learners[i]->notifyRestart( );

  if( iterative_dfs() == LIMITOUT ) {
    status = UNKNOWN;

//     // BUGGY!! it should be init_level, but the roadef challenge model doesn't like it.
//      std::cout << "c " << reinit << " " << BTSLIST.size << std::endl
//  	      << "c reset structs (" 
// 	       << ((reinit >= 0 && BTSLIST.size >= reinit) ? "full" : "standard") 
// 	       << ")" << std::endl; 

    
//     std::cout 
//       << ((reinit >= 0 && BTSLIST.size >= reinit) ? "f" : "s") ;
//     std::cout.flush();
    
//    backtrackTo(init_level-(reinit >= 0 && BTSLIST.size >= reinit));
    backtrackTo(init_level);

//     std::cout << "X[0,0] after a reset: ";
//     variables[0]->print(std::cout);
//     std::cout << std::endl;

    if( decay > 0.0 )
      weightDecay( decay );
    if( sat ) {
      sat->forget();
    }
  }

  if(BTSLIST.empty())
    BTSLIST.push(BACKTRACKS);
  else
    BTSLIST.push(BACKTRACKS - BTSLIST.back());

#ifdef _WEIGHT_STATS

  ENDTIME = (getRunTime() - STARTTIME);

  // 	  outfile << iteration << "\t"
  // 		  << ENDTIME << "\t"
  // 		  << NODES << "\t"
  // 		  << setw(10) << (getGiniCoef()) << "\t"
  // 		  << setw(10) << (getHooverCoef()) << std::endl;

  rts_time.push( ENDTIME );
  rts_node.push( NODES );
  gini_vlist.push( getGiniVarCoef() );
  gini_clist.push( getGiniConCoef() );
  hoover_list.push( getHooverCoef() );

  //	  std::cout << hoover_list.back() <<std::endl;

#endif


  FAILLIMIT = (FAILURES + fail_increment);
  return( (NDSLIMIT <= 0 || NODES < NDSLIMIT) &&
	  (TIMELIMIT <= .0 || (getRunTime() - STARTTIME < TIMELIMIT)) );
}

int Solver::solve_and_restart( const int policy,
			       const unsigned int base,
			       const double factor,
			       const double decay,
			       const int reinit)
{

  //std::cout << "unary constraints " << sUnaryCons.size << std::endl;

  mistral_solver = this;
  signal(SIGINT,Mistral_SIGINT_handler);
  unsigned int iteration = 0;

  try {
    if( init_level < level ) init_level = level;

    if( status == UNKNOWN ) {
      ////////////////////////////////////////////////////////////
      //double PRESOLVETIME = getRunTime();
      ////////////////////////////////////////////////////////////

      int numvar, numval, numcon;

      if( verbosity > 1 )
	cout << endl
	     << "c  ==================================[ Mistral ]=====================================" << endl
	     << "c  |      SEARCH  STATISTICS           |           PROBLEM STATISTICS               |" << endl
	     << "c  |     Fails  |       Nodes  CPUTime |  #Variables       #values     #constraints |" << endl
	     << "c  ==================================================================================" << endl;
      else if( verbosity > 0 )
	cout << "c ==============[ Mistral - restarts ]=============" << endl
	     << "c | restart |    limit |   # of nodes |  cpu time |" << endl
	     << "c =================================================" << endl;


      //       std::cout << "BEFORE PRESOLVE AT LEVEL " << level << " ";
      //       variables[41]->print(std::cout);
      //       std::cout << std::endl;

      STARTTIME = getRunTime();
      if( presolve() == UNKNOWN ) {

	// 	variables[41]->print(std::cout);
	// 	std::cout << std::endl;

	

	//  	print(cout);
	// 	cout << endl;

#ifdef _WEIGHT_STATS

	ENDTIME = (getRunTime() - STARTTIME);

	if(ENDTIME < 0.0001) ENDTIME = 0;

	rts_time.push( ENDTIME );
	rts_node.push( NODES );
	gini_vlist.push( getGiniVarCoef() );
	gini_clist.push( getGiniConCoef() );
	hoover_list.push( getHooverCoef() );

	std::cout << "\n\nSTORE STATS AFTER RESTART: " << (rts_time.size) << "\n" ;

#endif




	FAILLIMIT = base+FAILURES;
	fail_increment = base;
	while( status == UNKNOWN ) {
	  ++iteration;
	  switch( policy ) {
	  case GEOMETRIC : fail_increment =
	      lround((double)(fail_increment * factor));
	    break;
	  case LUBY : fail_increment =
	      lround((double)(base * luby_seq(iteration+1)));
	    break;
	  default : fail_increment =
	      lround((double)(base * luby_seq(iteration) * factor));
	  }


	  if( verbosity > 1 ) {

	    numvar = (empty - future);
	    numcon = constraints.size;
	    numval = 0;
	    for(int i=0; i<length; ++i)
	      numval += sequence[i]->domsize();

	    ENDTIME = getRunTime() - STARTTIME;
	    if( ENDTIME < 0.001 ) ENDTIME = 0.001;

	    cout << "c  | " << setw(10) << (FAILLIMIT-FAILURES)
		 << " | " << setw(10) << NODES
		 << " " << setw(9) << ENDTIME
		 << " | " << setw(11) << numvar
		 << "  " << setw(13) << numval
		 << "  " << setw(13) << numcon
		 << "  |" << endl;
	  } else if( verbosity > 0 ) {

	    ENDTIME = getRunTime() - STARTTIME;
	    if( ENDTIME < 0.001 ) ENDTIME = 0.001;


	    cout << "c | " << left  << setw(7)  << (iteration)
		 << " | "  << right << setw(8)  << (FAILLIMIT-FAILURES)
		 << " | "  << right << setw(12) << NODES
		 << " | "  << right << setw(9) << ENDTIME
		 << " |"   << endl;

	  }

	  // 	  if( randomizedRestart ) randomizeSequence();
	  // 	  int i=learners.size;

	  // 	  if(function)
	  // 	    function->initialise();

	  // 	  while( i-- )
	  // 	    learners[i]->notifyRestart( );

	  // 	  if( iterative_dfs() == LIMITOUT ) {
	  // 	    status = UNKNOWN;
	  // 	    backtrackTo(0);
	  // 	    if( decay > 0.0 )
	  // 	      weightDecay( decay );
	  // 	    if( sat ) {
	  // 	      sat->forget();
	  // 	    }
	  // 	  }

	  // 	  //if( sat )
	  // 	  //BTSLIST.push( sat->stats.learnt_avg_size );
	  // 	  if(BTSLIST.empty())
	  // 	    BTSLIST.push(BACKTRACKS);
	  // 	  else
	  // 	    BTSLIST.push(BACKTRACKS - BTSLIST.back());


	  // 	  FAILLIMIT = (FAILURES + fail_increment);
	  // 	  if( (NDSLIMIT > 0 && NODES >= NDSLIMIT) ||
	  // 	      (TIMELIMIT > .0 && getRunTime() - STARTTIME >= TIMELIMIT) ) break;

	  // #ifdef _PLOT_GINI


	  // #endif

#ifdef _WEIGHT_STATS

	  int i=0;
	  while(i < learners.size && learners[i]->getType() != Weighter::WDG) ++i;
	  if(i < learners.size) {
	    ((WeighterDegree*)(learners[i]))->init_choices();
	  }

#endif

	  if(!restart(decay, reinit)) {
	    break;
	  }

	}
      } 

      closeSearch();
      //stopWatch();
    }


    if( verbosity > 1 ) {
      cout << "c  ==================================================================================" << endl;
      if( status == SAT )
	cout << "c  |                               SATISFIABLE!                                     |" << endl;
      else if( status == UNSAT )
	cout << "c  |                              UNSATISFIABLE!                                    |" << endl;
      else if( status == UNKNOWN )
	cout << "c  |                                 UNKNOWN!                                       |" << endl;
      if( ENDTIME < 0.001 ) ENDTIME = 0.001;
      cout << "c  ==================================================================================" << endl ;
      cout << "c  |      Fails  |      Nodes  |       Propags  |  CPU Time |   Fails/s |   Nodes/s |" << endl
	   << "c  | " << setw(10) << FAILURES
	   << "  | " << setw(10) << NODES
	   << "  | " << setw(13) << PROPAGS
	   << "  | " << setw(9) << ENDTIME
	   << " | " << setw(9) << int(double(FAILURES)/ENDTIME)
	   << " | " << setw(9) << int(double(NODES)/ENDTIME)
	   << " |" << endl
	   << "c  ==================================================================================" << endl << endl;
    } else if( verbosity > 0 )
      cout << "c =================================================" << endl;

  } catch (const abort_search& ex) {
    std::cout << "c *        Abort search exception caught!         *" << std::endl
	      << "c *************************************************" << std::endl;
    status = UNKNOWN;
    //exit(1);
  } catch (...) {
    std::cout << "Fatal error, unkown exception" << std::endl;
    exit(1);
  }


  // #ifdef _PLOT_GINI

  //   outfile.close();

  // #endif

  FAILLIMIT = 0;

  return status;
}





#ifdef _WEIGHT_STATS

void Solver::print_search_stats(std::string& key)
{
  std::cout << "\n\nPRINT TO " << key << std::endl;

  int i=0;
  while(i < learners.size && learners[i]->getType() != Weighter::WDG) ++i;
  if(i < learners.size) {

    //int bts = BACKTRACKS;
    //if(BTSLIST.size)
    //bts = (BACKTRACKS - BTSLIST.back());

    std::string chname = ("tmp_choices_"+key);
    std::ofstream outch(chname.c_str(), ios_base::out);
    ((WeighterDegree*)(learners[i]))->print_choices(outch);
    outch.close();

    std::string coname = ("tmp_coefs_"+key);
    std::ofstream outco(coname.c_str(), ios_base::out);


    assert( rts_time.size == rts_node.size );
    assert( rts_time.size == BTSLIST.size+1 );

    for(i=0; i<rts_time.size; ++i) {
      outco << i << "\t"
	    << rts_time[i] << "\t"
	    << rts_node[i] << "\t"
	    << (((double)(rts_node[i])) / ((double)NODES)) << "\t"
	    << setw(10) << gini_vlist[i] << "\t"
	    << setw(10) << gini_clist[i] << "\t"
	    << setw(10) << hoover_list[i] << std::endl;
    }
    outco.close();
  }
}
#endif


int Solver::random_probe( const unsigned int iterations,
			  const unsigned int limit )
{
  mistral_solver = this;
  signal(SIGINT,Mistral_SIGINT_handler);

  if( status == UNKNOWN ) {
    if( presolve() == UNKNOWN ) {

      int i = iterations;
      status = LIMITOUT;

      while( i-- && status == LIMITOUT ) {
	status = UNKNOWN;

	if( verbosity > 0 && !((i+1) % (iterations/10)))
	  cout << ".";

	FAILLIMIT = FAILURES+limit;

	if( randomizedRestart ) randomizeSequence();
	iterative_dfs();
	if( status == LIMITOUT )
	  backtrackTo(0);

	if( (NDSLIMIT > 0 && NODES >= NDSLIMIT)
	    || ((getRunTime() - STARTTIME) >= RPTIMELIMIT)
	    )
	  break;
      }

      FAILLIMIT = 0;

      if( status == LIMITOUT )
	status = UNKNOWN;
    }

    if(status != UNKNOWN) {
      closeSearch();
      //  stopWatch();
    }

  }

  return status;
}


int Solver::sacPreprocess( const bool comp, const int wtype )
{
  mistral_solver = this;
  signal(SIGINT,Mistral_SIGINT_handler);

  WeighterSAC *ws = NULL;
  WeighterImpact *wi = NULL;

  if( presolve() == UNKNOWN ) {

    if( (wtype & Weighter::WLD) == Weighter::WLD )
      setLearner( Weighter::WLD );
    else if( wtype & Weighter::WDG )
      setLearner( Weighter::WDG );
    if( wtype & Weighter::IPT )
      wi = ((WeighterImpact*)(setLearner( Weighter::IPT )));

    if( wtype & Weighter::IPT ){
      ws = (WeighterISAC*)(setLearner( Weighter::ISAC )); // DG Added for doing sac initialization of impacts
      ws->complete = comp;


      for( int i=0; i<length; ++i ){
	if( variables[i]->getType() != VariableInt::RANGE )
	  std::fill( wi->decision_impact[i]+variables[i]->minCapacity(), wi->decision_impact[i]+variables[i]->maxCapacity(), 0 );
	else
	  wi->decision_impact[i][0] = wi->decision_impact[i][1] = 0;
      }
    }
    else
      {  ws = (WeighterSAC*)(setLearner( Weighter::SAC ));
	ws->complete = comp;
	cout << "SAC if statement" << endl;
      }

    DVOSingletonAC *h = new DVOSingletonAC( this, ws );

    DVO *oldh = heuristic;
    heuristic = h;


    ValSelector* oldvo[length];
    for( int i=0; i<length; ++i )
      {
	oldvo[i] = variables[i]->branch;
	if( h->isRange[i] )
	  variables[i]->branch = new ValSelectorSacRange( variables[i], &(h->SACdomain[i]) );
	else
	  variables[i]->branch = new ValSelectorSac( variables[i], &(h->SACdomain[i]) );
      }


    status = LIMITOUT;
    while( status == LIMITOUT && ws->isNotSac() )
      {
	status = UNKNOWN;

	FAILLIMIT = FAILURES+1;

	iterative_dfs();

	//cout << endl << (status == SAT ? "SAT" : (status == UNSAT ? "UNSAT" : (status == LIMITOUT ? "LIMITOUT" : "UNKNOWN"))) << endl;

	if( status == LIMITOUT )
	  backtrackTo(0);
      }


    FAILLIMIT = 0;


    delete h;
    heuristic = NULL;

    if( oldh )
      heuristic = oldh;

    for( int i=0; i<length; ++i )
      {
	delete variables[i]->branch;
	variables[i]->branch = oldvo[i];
      }

    delete learners.back();
    learners.pop();
  }

  if( status == LIMITOUT ) status = UNKNOWN;

  return status;
}


int Solver::checkSolution()
{
  int max_arity = 0;
  int i, j ;//, *aux;
  Vector<int> inconsistencies;
  bool assigned;

  for(i =0; i<constraints.size; ++i) {
    if(max_arity < constraints[i]->arity)
      max_arity = constraints[i]->arity;
  }

  int sol[max_arity];

  for(i =0; i<constraints.size; ++i) {
    j = constraints[i]->arity;

    //    while(j--)
    //       sol[j] = constraints[i]->scope[j]->min();

    //     if( constraints[i]->check(sol) ) {
    //       inconsistencies.push(i);

    //       exit(0);
    //       }

    assigned = true;
    while(j--) {
      if(!(constraints[i]->scope[j]->isGround())) assigned = false;
      sol[j] = constraints[i]->scope[j]->first();
    }

    if(assigned) {
      if( constraints[i]->check(sol) )
	inconsistencies.push( i );
    } else {
      //       aux = constraints[i]->supp;
      //       constraints[i]->supp = sol;

      //       cout << "find support" << endl;
      //       if( !constraints[i]->findSupport(0, sol[0]) )
      //       inconsistencies.push( i );
      //       constraints[i]->supp = aux;
    }
  }

  while( inconsistencies.size ) {
    inconsistencies.pop( j );
    cout << "c Error: Solution does not verify\nc\t";
    constraints[j]->print( cout );
    for(i=0; i<constraints[j]->arity; ++i) {
      cout << endl << "c\t" << (constraints[j]->scope[i]->getType()) << " ";
      constraints[j]->scope[i]->print( cout );
    }
    cout << endl << "s NOT SUPPORTED" << endl << endl;
    exit( 0 );
  }

  return 1;
}

unsigned long int Solver::getBacktracks() const {return BACKTRACKS;}
unsigned long int Solver::getNodes     () const {return NODES     ;}
unsigned long int Solver::getFailures  () const {return FAILURES  ;}
unsigned long int Solver::getChecks    () const {return CHECKS    ;}
unsigned long int Solver::getPropags   () const {return PROPAGS   ;}
double            Solver::getTime      () const {return TOTTIME   ;}


void Solver::printStatistics(std::ostream &log, const int stats)
{

  if(ENDTIME < .001)
    ENDTIME = .0;

  if(TOTTIME < .001)
    TOTTIME = .0;

  if( stats & OUTCOME ) {
    if(goal) {
      log << std::left << std::setw(30) << "c Best objective" << ":"
	  << std::right << std::setw(20) << (goal->upper_bound) ;
      if(status == OPT)
	log << " (optimal)" ;
      else
	log << " (may not be optimal)" ;
      log << std::endl;
    } else if(FIND_ALL || SOLUTIONS > 1) {
      log << std::left << std::setw(30) << "c #Solutions" << ":"
	  << std::right << std::setw(20) << SOLUTIONS ;
      if(status == OPT || status == UNSAT)
	log << " (complete)" ;
      else
	log << " (may not be complete)" ;
      log << std::endl;
    } else
      log << std::left << std::setw(30) << "c Result" << ":"
	  << std::right << std::setw(20) << ((status == SAT) ? "SATISFIABLE" : (status == UNSAT ? "UNSATISFIABLE" : "UNKNOWN")) << std::endl;
  }

  if( BTSLIST.empty() )
    BTSLIST.push( BACKTRACKS );

  if( stats & BTS )
    log << std::left << std::setw(30) << "c #Backtracks" << ":"
	<< std::right << std::setw(20) << BACKTRACKS << std:: endl;
  if( stats & NDS )
    log << std::left << std::setw(30) << "c #Nodes" << ":"
	<< std::right << std::setw(20) << NODES << std::endl;
  if( stats & FLS )
    log << std::left << std::setw(30) << "c #Failures" << ":"
	<< std::right << std::setw(20) << FAILURES << std::endl;
  if( stats & PPGS )
    log << std::left << std::setw(30) << "c #Propagations" << ":"
	<< std::right << std::setw(20) << PROPAGS << std::endl;
  if( stats & CKS )
    log << std::left << std::setw(30) << "c #Constraint checks" << ":"
	<< std::right << std::setw(20) << CHECKS << std::endl;
  if( stats & REST )
    log << std::left << std::setw(30) << "c #Restarts" << ":"
	<< std::right << std::setw(20) << BTSLIST.size << std::endl;
  if( stats & LASTBTS )
    log << std::left << std::setw(30) << "c #Backtracks on final restart" << ":"
	<< std::right << std::setw(20) << BTSLIST.back() << std::endl;
  if( stats & RUNTIME )
    log << std::left << std::setw(30) << "c Solving time (in sec.)" << ":"
	<< std::right << std::setw(20) << std::setprecision(5) << ENDTIME << std::endl ;
  if( stats & OPTTIME )
    log << std::left << std::setw(30) << "c Optimisation time (in sec.)" << ":"
	<< std::right << std::setw(20) << std::setprecision(5) << SOLTIME << std::endl ;
  if( stats & PROOFTIME )
    log << std::left << std::setw(30) << "c Proof time (in sec.)" << ":"
	<< std::right << std::setw(20) << std::setprecision(5) << (ENDTIME-SOLTIME) << std::endl ;
  if( stats & TOTALTIME )
    log << std::left << std::setw(30) << "c Total time (incl. build)" << ":"
	<< std::right << std::setw(20) << std::setprecision(5) << TOTTIME << std::endl ;
  if( stats & SPEEDB ) {
    log << std::left << std::setw(30) << "c Nodes/second" << ":";
    if(ENDTIME)
      log << std::right << std::setw(20) << (unsigned int)(((double)NODES)/ENDTIME) << std::endl;
    else
      log << std::right << "N/A" << std::endl;
  }
  if( stats & SPEEDP ) {
    log << std::left << std::setw(30) << "c Propagations/second" << ":";
    if(ENDTIME)
      log << std::right << std::setw(20) << (unsigned int)(((double)PROPAGS)/ENDTIME) << std::endl;
    else
      log << std::right << "N/A" << std::endl;
  }
  log.flush();
}


double Solver::getGiniVarCoef() const
{

  //std::cout << (heuristic->get_disjuncts()) << std::endl;

  int i, j, aux, orderedvar[length];
  double avg = variables[0]->weight, Gini=0;
  for(i=0; i<length; ++i)
    orderedvar[i] = i;
  for(i=1; i<length; ++i)
    {
      //       std::cout << (variables[i]->weight) << " =?= "
      // 		<< heuristic->get_value(variables[i]) << std::endl;

      //avg += variables[i]->weight;
      avg += heuristic->get_value(variables[i]);
      j=i;
      //while( j && variables[orderedvar[j]]->weight < variables[orderedvar[--j]]->weight )
      while( j && heuristic->get_value(variables[orderedvar[j]]) < heuristic->get_value(variables[orderedvar[j-1]]) )
	{
	  --j;
	  aux = orderedvar[j];
	  orderedvar[j] = orderedvar[j+1];
	  orderedvar[j+1] = aux;
	}
    }
  avg /= length;


  for(i=1; i<length; ++i) {
    assert(heuristic->get_value(variables[orderedvar[i]]) >= heuristic->get_value(variables[orderedvar[i-1]]));
  }

  for(i=0; i<length; ++i)
    Gini += ((i+1) * (((double)(heuristic->get_value(variables[orderedvar[i]]))) - avg));

  Gini /= (length * length * avg);
  Gini *= 2;

  return Gini;
}

double Solver::getGiniConCoef() const
{

  //std::cout << (heuristic->get_disjuncts()) << std::endl;

  int i, j, aux, orderedcon[constraints.size];
  double avg = constraints[0]->weight, Gini=0;
  for(i=0; i<constraints.size; ++i)
    orderedcon[i] = i;
  for(i=1; i<constraints.size; ++i)
    {
      //       std::cout << (variables[i]->weight) << " =?= "
      // 		<< heuristic->get_value(variables[i]) << std::endl;

      avg += constraints[i]->weight;
      //avg += heuristic->get_value(variables[i]);
      j=i;
      while( j && constraints[orderedcon[j]]->weight < constraints[orderedcon[j-1]]->weight )
	//while( j && heuristic->get_value(variables[orderedvar[j]]) < heuristic->get_value(variables[orderedvar[--j]]) )
	{
	  --j;
	  aux = orderedcon[j];
	  orderedcon[j] = orderedcon[j+1];
	  orderedcon[j+1] = aux;
	}
    }
  avg /= constraints.size;


  //   for(i=1; i<length; ++i) {
  //     assert(heuristic->get_value(variables[orderedvar[i]]) >= heuristic->get_value(variables[orderedvar[i-1]]));
  //   }

  for(i=0; i<constraints.size; ++i)
    //Gini += ((i+1) * (((double)(heuristic->get_value(variables[orderedvar[i]]))) - avg));
    Gini += ((i+1) * (((double)(constraints[orderedcon[i]]->weight)) - avg));

  Gini /= (constraints.size * constraints.size * avg);
  Gini *= 2;

  return Gini;
}


double Solver::getHooverCoef() const
{
  int i;
  double avg = 0.0, aux, Hoover = 0.0;

  for(i=0; i<length; ++i)
    //avg += sequence[i]->weight;
    avg += heuristic->get_value(variables[i]);

  //std::cout << std::endl << avg << "/" << length << std::endl;

  avg /= length;

  //std::cout << avg << std::endl;


  for(i=0; i<length; ++i) {
    //aux = sequence[i]->weight;
    aux = heuristic->get_value(variables[i]);
    if(aux > avg)
      Hoover += (aux - avg);
  }

  //std::cout << Hoover << std::endl;

  Hoover /= (avg * length);

  //std::cout << Hoover << std::endl;

  return Hoover;
}


void Solver::printWeightProfile(std::ostream& o, int limit, int threshold) const
{
  int i, j, aux, orderedvar[length], total=0;
  for(i=0; i<length; ++i)
    orderedvar[i] = i;
  for(i=1; i<length; ++i)
    {
      j=i;
      while( j && variables[orderedvar[j]]->weight > variables[orderedvar[j-1]]->weight )
	{
	  --j;
	  aux = orderedvar[j];
	  orderedvar[j] = orderedvar[j+1];
	  orderedvar[j+1] = aux;
	}
    }
  if( length < limit )
    limit = length;
  double n=0;
  for(i=0; i<limit; ++i) {
    if( variables[orderedvar[i]]->weight > (unsigned int)threshold )
      ++n;
    total += (variables[orderedvar[i]]->weight-1);
  }

  cout << "d TOTAL " << total << " WEIGHTEDVARS " << (n/((double)limit)) << "\nc ";
}

bool Solver::preprocess()
{
  for(int i=0; i<numvars; ++i)
    gacvarstack.push(variables[i]->id, Constraint::RANGETRIGGER);
  return filtering();
}

void Solver::store_solution()
{
  ++SOLUTIONS;

  for(int i=0; i<length; ++i) {
    solution[i] = max_solution[i] = variables[i]->value() ;
  }
  for(int i=length; i<numvars; ++i) {
    solution[i] = variables[i]->value() ;
    max_solution[i] = variables[i]->max() ;
  }

  checkSolution();

  SOLTIME = getRunTime() - STARTTIME;
  if(SOLTIME < .0001)
    SOLTIME = .0;

  status = SAT;
}

int Solver::solutionFound(int init_level)
{

  store_solution();

//   //   for(int i=0; i<numvars; ++i) {
//   //     cout << setw(4) << (i+1) << " " << variables[i]->value() << endl;
//   //   }

//   for(int i=0; i<length; ++i) {
//     solution[i] = max_solution[i] = variables[i]->value() ;
//   }
//   for(int i=length; i<numvars; ++i) {
//     solution[i] = variables[i]->value() ;
//     max_solution[i] = variables[i]->max() ;
//   }

//   checkSolution();
//   ++SOLUTIONS;

//   SOLTIME = getRunTime() - STARTTIME;
//   if(SOLTIME < .0001)
//     SOLTIME = .0;

//   status = SAT;
  if(function) function->execute();

  if(verbosity > 2  || (goal && verbosity > 0)) {

    if(verbosity > 1) {
      if(FAILLIMIT)
	cout << "c  ----------------------------------------------------------------------------------" << endl;
      cout << "c  | " << std::setw(10) << NODES << " NDS"
	   << std::setw(10) << std::setprecision(5)
	   << SOLTIME << " s" ;
      cout << "            solution: " << std::setw(8)
	   << (int)(goal ? goal->solution_score() : SOLUTIONS) ;

      if( verbosity > 2 ) {
	cout << ": {" ;
	for(int i=0; i<length-1; ++i)
	  cout << solution[i] << " ";
	cout << solution[length-1] << "}" ;
      }
      cout << endl;
      if(FAILLIMIT)
	cout << "c  ----------------------------------------------------------------------------------" << endl;

    } else {
      cout << left  << setw(13) << "c | solution: "
	   << left  << setw(8)  << (goal ? goal->solution_score() : SOLUTIONS)
	   << " | " << right << setw(12) << NODES
	   << " | "  << right << setw(9) << SOLTIME
	   << " |"   << endl;
    }
  }

  if(goal) {
    if( (status = goal->update()) != OPT ) {
      if(!optimisation) {
	status = SAT;
      } else if( level > init_level ) {
	backtrackLevel = level-1;
	//VariableInt *lastDecision = decision[level];
	SimpleUnaryConstraint last_decision = branching_decision[level];
	backtrackTo( backtrackLevel );
#ifdef _DEBUGSEARCH
      if(verbosity > 2) {
	std::cout << "c";
	for(int k=0; k<=level; ++k) std::cout << " ";
	last_decision.print(std::cout);
      }
#endif
      
      //lastDecision->branch->right();
      last_decision.right();
#ifdef _DEBUGSEARCH
      if(verbosity > 2) {
	//lastDecision->branch->printRight( std::cout );
	last_decision.print( std::cout );
	std::cout << std::endl;
      }
#endif

      } else {
	status = UNSAT;
      }
    }

    return status;
  }
  if(FIND_ALL && (int)SOLUTIONS != FIND_ALL) {
    if( level > init_level ) {
      status = UNKNOWN;
      backtrackLevel = level-1;
      //VariableInt *lastDecision = decision[level];
      SimpleUnaryConstraint last_decision = branching_decision[level];
      backtrackTo( backtrackLevel );
      //lastDecision->branch->right();
      last_decision.right();
    } else status = UNSAT;
  }

  return status;
}


/*
VariableInt *x, *y;

Constraint *con = NULL;
MistralNode<Constraint*> *nd;
nd = x->constraintsOnValue();
while( nextNode(nd) ) {
  con = nd->elt;
  for(int i=0; i<con->arity; ++i) {
    y = con->_scope[i];
    
    y->id;
    //
  }
 }


0 
cp_solver->variables.size;
*/

bool Solver::filtering()
{
  backtrackLevel = NOVAL;

#ifdef _DEBUGAC
      if(verbosity > 2) {
  cout << "\tGAC" << endl;
  for(int k=0; k<numvars; ++k) {
    variables[k]->print( cout );
    cout << endl;
  }
  cout << endl;
      }
#endif

#ifdef _DEBUGPROPAG
  //std::ostringstream o_propag;
  std::string prop_str;
  int size_before;
#endif

  bool consistent = ( !goal || goal->score() <= goal->upper_bound );
  Constraint *con = NULL;

  if( consistent ) {
    int sz = unaryCons.size, domain_event;
    while( consistent && sz-- ) {
      ++PROPAGS;
      if( !(unaryCons[sz]->propagate( )) ) consistent = false;
    }


    //int nuaryprop = 0;
    while( consistent && !sUnaryCons.empty() ) {
      ++PROPAGS;
      SimpleUnaryConstraint cons = sUnaryCons.pop();

      if(verbosity > 2) {
        std::cout << "PROPAGATE UNARY CONSTRAINT AT LEVEL " << level << ": ";
        cons.print(std::cout);
        std::cout << std::endl;
      }

//       if(cons.val > cons.var->min()) {
// 	++nuaryprop;
//       }
      consistent &= cons.propagate();
    }
//     if(nuaryprop)
//       //std::cout << "# unary propags: " 
//       std::cout << nuaryprop ;
// 	//<< std::endl; 

    if( consistent ) {
      VariableInt *pvar;
      MistralNode<Constraint*> *nd;

      bool fixedPoint = (gacvarstack.empty() && gacconstack.empty());

      while( consistent && !fixedPoint ) {

	if( !gacvarstack.empty() ) {

#ifdef _DEBUGPROPAGSTACK
      if(verbosity > 2) {
	  cout << "stack:";
	  int q = gacvarstack.head;
	  do {
	    q = gacvarstack.next[q];
	    cout << " " << q;
	  } while( q != gacvarstack.tail );
	  cout << endl;
      }
#endif

	  pvar = gacvarstack.pop( domain_event );

#ifdef _DEBUGPROPAGSTACK
      if(verbosity > 2) {
	  if( domain_event & Constraint::VALUETRIGGER ) {
	    cout << "\tValue event on " ;
	    pvar->print( cout );
	    cout << endl;
	    nd = pvar->constraintsOnValue();
	  } else if( domain_event & Constraint::RANGETRIGGER ) {
	    cout << "\tRange event on " ;
	    pvar->print( cout );
	    cout << endl;
	    nd = pvar->constraintsOnRange();
	  } else {
	    cout << "\tDomain event on " ;
	    pvar->print( cout );
	    cout << endl;
	    nd = pvar->constraintsOnDomain();
	  }
	  cout << endl;
      }
#endif

	  if( domain_event & Constraint::VALUETRIGGER ) {
	    nd = pvar->constraintsOnValue();
	    pvar->unLink();
	  } else if( domain_event & Constraint::RANGETRIGGER )
	    nd = pvar->constraintsOnRange();
	  else
	    nd = pvar->constraintsOnDomain();

	  while( nextNode(nd) ) {
	    con = nd->elt;
	    if( con->delayed ) {
	      gacconstack.push( con, pvar );
	      continue;
	    }

#ifdef _DEBUGPROPAG
	    std::ostringstream o_propag;
	    size_before = 0;
      if(verbosity > 2) {
	    o_propag << "PROPAGATE \n" ;
	    con->print( o_propag );
	    o_propag << endl;
	    for(int a=0; a<con->arity; ++a) {
	      size_before += con->scope[a]->domsize();
	      con->scope[a]->print( o_propag );
	      o_propag << " ";
	    }
	    o_propag << endl;
	    prop_str = o_propag.str();
      }
#endif

	    // 	bool is42 = false;
	    // 	for(int x=0; !is42 && x<con->arity; ++x)
	    // 	  is42 |= (con->scope[x]->id == -41);
	    // 	if(is42) {
	    // 	  std::ostringstream o_propag;
	    // 	  size_before = 0;

	    // 	  o_propag << " PROPAGATE \n" ;
	    // 	  con->print( o_propag );
	    // 	  o_propag << endl;
	    // 	  for(int a=0; a<con->arity; ++a) {
	    // 	    size_before += con->scope[a]->domsize();
	    // 	    con->scope[a]->print( o_propag );
	    // 	    o_propag << " ";
	    // 	  }
	    // 	  o_propag << endl;
	    // 	  prop_str = o_propag.str();
	    // 	}

	    ++PROPAGS;
	    if( !(consistent = con->propagate(nd->index, domain_event)) ) {
	      learnFailure( con );
	      break;
	    }

	    // 	if(is42) {
	    //  	  for(int a=0; a<con->arity; ++a) {
	    //  	    size_before -= con->scope[a]->domsize();
	    //  	  }

	    //  	  if(size_before) {
	    // 	    cout << prop_str;
	    // 	    for(int a=0; a<con->arity; ++a) {
	    // 	      con->scope[a]->print( cout );
	    // 	      cout << " ";
	    // 	    }
	    // 	    cout << endl;
	    // 	    con->print( cout );
	    // 	    cout << " OK" << endl << endl;
	    // 	  }
	    // 	}

#ifdef _DEBUGPROPAG
      if(verbosity > 2) {
	    for(int a=0; a<con->arity; ++a) {
	      size_before -= con->scope[a]->domsize();
	    }
	    if(size_before) {
	      cout << prop_str;
	      for(int a=0; a<con->arity; ++a) {
		con->scope[a]->print( cout );
		cout << " ";
	      }
	      cout << endl;
	      con->print( cout );
	      cout << " OK" << endl << endl;
	    }
      }
#endif
	  }

	} else if(consistent) {
	  if( gacconstack.empty() ) fixedPoint = true;
	  else {
	    con = gacconstack.pop( sz );

#ifdef _DEBUGPROPAG
	    std::ostringstream o_propag;
	    size_before = 0;
      if(verbosity > 2) {
	    o_propag << "DELAYED: " ;
	    variables[sz]->print( o_propag );
	    o_propag << endl;
	    for(int a=0; a<con->arity; ++a) {
	      con->scope[a]->print( o_propag );
	      o_propag << " ";
	    }
	    o_propag << endl;
	    con->print( o_propag );
	    o_propag << endl;
      }
#endif
	    ++PROPAGS;
	    if( !(consistent = con->propagate()) ) {
	      learnFailure( con );
	      break;
	    }

#ifdef _DEBUGPROPAG
      if(verbosity > 2) {
	    for(int a=0; a<con->arity; ++a) {
	      size_before -= con->scope[a]->domsize();
	    }
	    if(size_before) {
	      cout << prop_str;
	      for(int a=0; a<con->arity; ++a) {
		con->scope[a]->print( cout );
		cout << " ";
	      }
	      cout << endl;
	      con->print( cout );
	      cout << " OK" << endl << endl;
	    }
      }
#endif

	  }
	}
      }
    }
  }

#ifdef _DEBUGPROPAG
      if(verbosity > 2) {
  if( !consistent ) {
    cout << prop_str;
    if( con )
      con->print( cout );
    cout << " FAILS" << endl << endl;
  }
      }
#endif

  gacvarstack.clear();
  gacconstack.clear();

#ifdef _DEBUGAC
      if(verbosity > 2) {
  if( consistent ) {
    for(int k=0; k<numvars; ++k) {
      variables[k]->print( cout );
      cout << endl;
    }
    cout << endl;
  }
  else cout << "FAILURE" << endl;
      }
#endif


  if( !consistent && backtrackLevel > level )
    backtrackLevel = level-1;

  return consistent;
}



int Solver::startNewSearch()
{
  optimisation = false;
  mistral_solver = this;
  signal(SIGINT,Mistral_SIGINT_handler);
  if( status == UNKNOWN) presolve();
  else {
    STARTTIME = getRunTime();
    closeSearch();
    //stopWatch();
  }
  return status;
}


int Solver::getNextSolution()
{
  if( SOLUTIONS ) {
    if( level <= init_level ) {
      status = UNSAT;
      return UNSAT;
    }
    status = UNKNOWN;
    backtrackLevel = level-1;
    //VariableInt *lastDecision = decision[level];
    SimpleUnaryConstraint last_decision = branching_decision[level];

#ifdef _DEBUGSEARCH
    if( level > backtrackLevel+1 ) {
      cout << "c";
      for(int k=0; k<=level; ++k) cout << " ";
      cout << " backjump to level " << backtrackLevel << endl;
    }
#endif

    backtrackTo( backtrackLevel );

#ifdef _DEBUGSEARCH
    cout << "c";
    //for(int k=0; k<=level; ++k) cout << " ";
    cout << " ";
    //lastDecision->branch->printRight( cout );
    last_decision.print( cout );
    cout << endl;
#endif

    //lastDecision->branch->right();
    last_decision.right();
  }

  int res = iterative_dfs();
  closeSearch();
  //stopWatch();
  return res;
}

void Solver::setForgetfulness( const double f )
{
  if( !sat ) {
    sat = new ConstraintClauseBase( this );
    sat->initialise();
  }
  sat->params.forgetfulness = f;
}

void Solver::setRestartNogood()
{
  if(status == UNKNOWN) {
    if( !sat ) {
      sat = new ConstraintClauseBase( this );
      sat->initialise();
    }
    setLearner( Weighter::RNGD );
    ((WeighterRestartNogood*)(learners.back()))->sat = sat;
  }
}

WeighterRestartGenNogood* Solver::setRestartGenNogood()
{

  ConstraintGenNogoodBase* base = new ConstraintGenNogoodBase(this);
  setLearner( Weighter::RGNGD );
  ((WeighterRestartGenNogood*)(learners.back()))->base = base;

  return ((WeighterRestartGenNogood*)(learners.back()));
}

void Solver::printPython() const
{
  print( cout );

//   for(int i = 0; i < variables.size; ++i)
//     {
//       std::cout << "\t" ;
//       variables[i]->print(std::cout);
//       std::cout << std::endl;
//     }
//   //o << endl;

}

void Solver::print(std::ostream& o) const
{
  if(goal) {
    goal->print(o);
    o << endl << "subject to" << endl;
  }

  int i;
  o << "Variables: " << endl;
  for(i = 0; i < variables.size; ++i)
    {
      o << "\t" ;
      variables[i]->print(o);
      o << std::endl;
    }
  o << endl;
  o << "Constraints:" << endl;
  for(i = 0; i < constraints.size; ++i)
    {
      o << "\t" ;
      constraints[i]->print(o);
      o << std::endl;
    }
  o << endl;
}

void Solver::printXML(std::ostream& o) const
{
  int i,j;

  // header
  o << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>"
    << endl << endl << "<instance>" << endl
    << "<presentation name=\"?\" format=\"XCSP 2.0\"/>"
    << endl << endl;


  // domains
  int nDomains=0;
  bool known;
  BitSet domain[length];
  int domainSize[length];
  int domainIndex[length];
  for(i=0; i<length; ++i)
    {
      known = false;
      for(j=0; !known && j<nDomains; ++j)
	known = ( sequence[i]->included( domain[j] ) &&
		  sequence[i]->include( domain[j] ) );
      if(!known) {
	domainIndex[i] = nDomains;
	domain[nDomains].init( sequence[i]->min(),
			       sequence[i]->max(),
			       BitSet::empt);
	sequence[i]->unionTo( domain[nDomains] );
	domainSize[i] = sequence[i]->domsize();
	++nDomains;
      } else {
	domainIndex[i] = j;
      }
    }
  o << "<domains nbDomains=\"" << nDomains << "\">" << endl;
  for(i=0; i<nDomains; ++i)
    {
      cout << "<domain name=\"D" << i
	   << "\" nbValues=\"" << domainSize[i] << "\">" ;
      domain[i].print( cout, "   " );
      cout << "</domain>" << endl;
    }
  cout << "</domains>" << endl << endl;

}



void Solver::setBranching(BuildObject **x, const int l, const char* val_name)
{
  string Val = val_name;
  VariableInt *temp;
  int i=l;
  if(Val == "Lex") {
    while( i-- ) {
      if( x[i]->isSearchable() && ((temp = x[i]->getVariable()) != NULL) )
	if( !temp->isGround() && !temp->branch ) {
	  if( temp->getType() != VariableInt::INTLIST ) {
	    temp->branch = new ValSelectorMin( temp );
	  } else {
	    temp->branch = new ValSelectorFirst( temp );
	  }
	}
    }
  } else if(Val == "AntiLex") {
    while( i-- )
      if( x[i]->isSearchable() && ((temp = x[i]->getVariable()) != NULL) )
	if( !temp->isGround() && !temp->branch ) {
	  if( temp->getType() != VariableInt::INTLIST ) {
	    temp->branch = new ValSelectorMax( temp );
	  } else {
	    temp->branch = new ValSelectorFirst( temp );
	  }
	}
  } else if(Val == "Random") {
    while( i-- )
      if( x[i]->isSearchable() && ((temp = x[i]->getVariable()) != NULL) )
	if( !temp->isGround() && !temp->branch ) {
	  if( temp->getType() == VariableInt::RANGE ) {
	    temp->branch = new ValSelectorRandMinMax( temp );
	  } else if( temp->getType() != VariableInt::INTLIST ) {
	    temp->branch = new ValSelectorRand( temp );
	  } else {
	    temp->branch = new ValSelectorFirst( temp );
	  }
	}
  } else if(Val == "RandomMinMax") {
    while( i-- )
      if( x[i]->isSearchable() && ((temp = x[i]->getVariable()) != NULL) )
	if( !temp->isGround() && !temp->branch ) {
	  if( temp->getType() != VariableInt::INTLIST ) {
	    temp->branch = new ValSelectorRandMinMax( temp );
	  } else {
	    temp->branch = new ValSelectorFirst( temp );
	  }
	}
  } else if(Val == "DomainSplit") {
    while( i-- ) {
      if( x[i]->isSearchable() && ((temp = x[i]->getVariable()) != NULL) ) {
	if( !temp->isGround() && !temp->branch ) {
	  if( temp->getType() == VariableInt::INTLIST ) {
	    temp->branch = new ValSelectorFirst( temp );
	  } else if( temp->domsize() > 3 && temp->getType() != VariableInt::LIST )
	    temp->branch = new ValSelectorSplit( temp );
	} else {
	  temp->branch = new ValSelectorMin( temp );
	}
      }
    }
  }
}

//void Solver::setHeuristic(string Heu, const int rdz)
void Solver::setHeuristic(const char* var_name, const char* val_name, const int rdz)
{
  if(status == UNKNOWN) {

    string Heu = var_name;
    string Val = val_name;

    if( rdz > 1 ) {
      if( Heu == "MinDomain" )
	heuristic = new GenericRandomDVO<VarSelectorDomain>(this, rdz);
      else if( Heu == "Lex" )
	heuristic = new DVOLexicographic(this);
      else if( Heu == "AntiLex" )
	heuristic = new DVOAntiLex(this);
      else if( Heu == "MaxDegree")
	heuristic = new GenericRandomDVO<VarSelectorDegree>(this, rdz);
      else if( Heu == "MinDomainMinVal")
	heuristic = new GenericRandomDVO<VarSelectorDomainMinVal>(this, rdz);
      else if( Heu == "Random")
	heuristic = new DVORandom(this);
      else if( Heu == "MinDomainMaxDegree")
	heuristic = new GenericRandomDVO<VarSelectorDomainDegree>(this, rdz);
      else if( Heu == "DomainOverDegree")
	heuristic = new GenericRandomDVO<VarSelectorDomainOverDegree>(this, rdz);
      else if( Heu == "DomainOverWLDegree") {
	setLearner( Weighter::WLD );
	heuristic =  new GenericRandomDVO<VarSelectorDomainOverWeight>(this, rdz);
      }
      else if( Heu == "DomainOverWDegree") {
	setLearner( Weighter::WDG );
	heuristic =  new GenericRandomDVO<VarSelectorDomainOverWeight>(this, rdz);
      }
      else if( Heu == "Neighbour")
	heuristic =  new GenericRandomDVO<VarSelectorNeighborDD>(this, rdz);
      else if( Heu == "Impact") {
	Impact H(rdz);
	heuristic = H.extract(this);
      }
      else if( Heu == "ImpactOverDegree") {
	ImpactOverDeg H(rdz);
	heuristic = H.extract(this);
      }
      else if( Heu == "ImpactOverWDegree") {
	ImpactOverWDeg H(rdz);
	heuristic = H.extract(this);
      }
      else if( Heu == "ImpactOverWLDegree") {
	ImpactOverWLDeg H(rdz);
	heuristic = H.extract(this);
      }
      else if( Heu == "Scheduling") {
	OSP H(rdz, (Val == "Promise" ? 1 : 0), OSP::DOM_O_TASKWEIGHT);
	heuristic = H.extract(this);
      }
    } else {
      if( Heu == "No" )
	heuristic = new DVONoOrder(this);
      else if( Heu == "MinDomain" )
	heuristic = new GenericDVO<VarSelectorDomain>(this);
      else if( Heu == "Lex" )
	heuristic = new DVOLexicographic(this);
      else if( Heu == "AntiLex" )
	heuristic = new DVOAntiLex(this);
      else if( Heu == "MaxDegree")
	heuristic = new GenericDVO<VarSelectorDegree>(this);
      else if( Heu == "MinDomainMinVal")
	heuristic = new GenericDVO<VarSelectorDomainMinVal>(this);
      else if( Heu == "Random")
	heuristic = new DVORandom(this);
      else if( Heu == "MinDomainMaxDegree")
	heuristic = new GenericDVO<VarSelectorDomainDegree>(this);
      else if( Heu == "DomainOverDegree")
	heuristic = new GenericDVO<VarSelectorDomainOverDegree>(this);
      else if( Heu == "DomainOverWLDegree") {
	setLearner( Weighter::WLD );
	heuristic =  new GenericDVO<VarSelectorDomainOverWeight>(this);
      }
      else if( Heu == "DomainOverWDegree") {
	setLearner( Weighter::WDG );
	heuristic =  new GenericDVO<VarSelectorDomainOverWeight>(this);
      }
      else if( Heu == "Neighbour")
	heuristic =  new GenericDVO<VarSelectorNeighborDD>(this);
      else if( Heu == "Impact") {
	Impact H;
	heuristic = H.extract(this);
      }
      else if( Heu == "ImpactOverDeg") {
	ImpactOverDeg H;
	heuristic = H.extract(this);
      }
      else if( Heu == "ImpactOverWDeg") {
	ImpactOverWDeg H;
	heuristic = H.extract(this);
      }
      else if( Heu == "ImpactOverWLDeg") {
	ImpactOverWLDeg H;
	heuristic = H.extract(this);
      }
      else if( Heu == "Scheduling") {
	OSP H(1, (Val == "Promise" ? 1 : 0), OSP::DOM_O_TASKWEIGHT);
	heuristic = H.extract(this);
      }
    }


    //   cout << endl << Val << endl;

    int i;
    if(Val == "Lex") {
      i = length;
      while( i-- ) {
	if( !variables[i]->isGround() && !variables[i]->branch ) {
	  if( variables[i]->getType() != VariableInt::INTLIST ) {
	    variables[i]->branch = new ValSelectorMin( variables[i] );
	  } else {
	    variables[i]->branch = new ValSelectorFirst( variables[i] );
	  }
	}
      }
    } else if(Val == "AntiLex") {
      i = length;
      while( i-- )
	if( !variables[i]->isGround() && !variables[i]->branch ) {
	  if( variables[i]->getType() != VariableInt::INTLIST ) {
	    variables[i]->branch = new ValSelectorMax( variables[i] );
	  } else {
	    variables[i]->branch = new ValSelectorFirst( variables[i] );
	  }
	}
    } else if(Val == "Random") {
      i = length;
      while( i-- ) {

	// 	variables[i]->print( std::cout );
	//	std::cout << " ";

	if( !variables[i]->isGround() && !variables[i]->branch ) {
	  if( variables[i]->getType() == VariableInt::RANGE ) {

	    //	    std::cout << "R -> randminmax" <<std::endl;

	    variables[i]->branch = new ValSelectorRandMinMax( variables[i] );
	  } else if( variables[i]->getType() != VariableInt::INTLIST ) {

	    //	    std::cout << "* -> rand" <<std::endl;

	    variables[i]->branch = new ValSelectorRand( variables[i] );
	  } else {

	    //	    std::cout << "* -> first" <<std::endl;

	    variables[i]->branch = new ValSelectorFirst( variables[i] );
	  }
	} //else std::cout << "ground or already done" << std::endl;
      }
    } else if(Val == "RandomMinMax") {
      i = length;
      while( i-- ) {

	//	variables[i]->print( std::cout );
	//	std::cout << " ";

	if( !variables[i]->isGround() && !variables[i]->branch ) {
	  if( variables[i]->getType() != VariableInt::INTLIST ) {

	    //	    std::cout << "* -> randminmax" <<std::endl;

	    variables[i]->branch = new ValSelectorRandMinMax( variables[i] );
	  } else {

	    //	    std::cout << "* -> first" <<std::endl;

	    variables[i]->branch = new ValSelectorFirst( variables[i] );
	  }
	} //else std::cout << "ground or already done" << std::endl;
      }
    } else if(Val == "RandomSplit") {
      i = length;
      while( i-- )
	if( !variables[i]->isGround() && !variables[i]->branch ) {
	  if( variables[i]->getType() == VariableInt::INTLIST ) {
	    variables[i]->branch = new ValSelectorFirst( variables[i] );
	  } else {
	    if( variables[i]->getType() != VariableInt::LIST ) {
	      if( variables[i]->domsize() <= 2 )
		variables[i]->branch = new ValSelectorRand( variables[i] );
	      else
		variables[i]->branch = new ValSelectorRandomSplit( variables[i] );
	    } else
	      variables[i]->branch = new ValSelectorFirst( variables[i] );
	  }
	}
    } else if(Val == "DomainSplit") {
      i = length;
      while( i-- )
	if( !variables[i]->isGround() && !variables[i]->branch ) {
	  if( variables[i]->getType() == VariableInt::INTLIST ) {
	    variables[i]->branch = new ValSelectorFirst( variables[i] );
	  } else if( variables[i]->domsize() > 3 && variables[i]->getType() != VariableInt::LIST )
	    variables[i]->branch = new ValSelectorSplit( variables[i] );
	} else {
	  variables[i]->branch = new ValSelectorMin( variables[i] );
	}
    }
  }
}



void Solver::setAntiLex() {
  for(VariableInt **x = sequence; x!=empty; ++x) {
    delete (*x)->branch;
    (*x)->branch = new ValSelectorMax(*x);
  }
}

void Solver::setLex() {
  for(VariableInt **x = sequence; x!=empty; ++x) {
    delete (*x)->branch;
    (*x)->branch = new ValSelectorMin(*x);
  }
}

void Solver::setSplit() {
  for(VariableInt **x = sequence; x!=empty; ++x) {
    delete (*x)->branch;
    (*x)->branch = new ValSelectorSplit(*x);
  }
}

void Solver::setRandSplit() {
  for(VariableInt **x = sequence; x!=empty; ++x) {
    delete (*x)->branch;
    (*x)->branch = new ValSelectorRandomSplit(*x);
  }
}

void Solver::setRandMinMax() {
  for(VariableInt **x = sequence; x!=empty; ++x) {
    delete (*x)->branch;
    (*x)->branch = new ValSelectorRandMinMax(*x);
  }
}


void Solver::setLex(VarArray& scope) {
  BuildObject *bvar[scope.size()];
  int i, l = scope.size();
  for(i=0; i<l; ++i)
    bvar[i] = scope[i].var_ptr_;
  setLex(bvar, l);
}

void Solver::setAntiLex(VarArray& scope) {
  BuildObject *bvar[scope.size()];
  int i, l = scope.size();
  for(i=0; i<l; ++i)
    bvar[i] = scope[i].var_ptr_;
  setAntiLex(bvar, l);
}

void Solver::setSplit(VarArray& scope) {
  BuildObject *bvar[scope.size()];
  int i, l = scope.size();
  for(i=0; i<l; ++i)
    bvar[i] = scope[i].var_ptr_;
  setSplit(bvar, l);
}

void Solver::setRandSplit(VarArray& scope) {
  BuildObject *bvar[scope.size()];
  int i, l = scope.size();
  for(i=0; i<l; ++i)
    bvar[i] = scope[i].var_ptr_;
  setRandSplit(bvar, l);
}

void Solver::setRandMinMax(VarArray& scope) {
  BuildObject *bvar[scope.size()];
  int i, l = scope.size();
  for(i=0; i<l; ++i)
    bvar[i] = scope[i].var_ptr_;
  setRandMinMax(bvar, l);
}

void Solver::setLex(BuildObject **bvar, const int l) {
  VariableInt *x;
  BuildObject *bv;
  for(int i=0; i<l; ++i) {
    bv = bvar[i];
    if(bv != bv->getBuildObject())
      bv = bv->getBuildObject();
    if( bv->isSearchable() &&
	((x = bv->getVariable()) != NULL)) {
      delete x->branch;
      x->branch = new ValSelectorMin( x );
    }
  }
}

void Solver::setAntiLex(BuildObject **bvar, const int l) {
  VariableInt *x;
  BuildObject *bv;
  for(int i=0; i<l; ++i) {
    bv = bvar[i];
    if(bv != bv->getBuildObject())
      bv = bv->getBuildObject();
    if( bv->isSearchable() &&
	((x = bv->getVariable()) != NULL)) {
      delete x->branch;
      x->branch = new ValSelectorMax( x );
    }
  }
}

void Solver::setSplit(BuildObject **bvar, const int l) {
  VariableInt *x;
  BuildObject *bv;
  for(int i=0; i<l; ++i) {
    bv = bvar[i];
    if(bv != bv->getBuildObject())
      bv = bv->getBuildObject();
    if( bv->isSearchable() &&
	((x = bv->getVariable()) != NULL)) {
      delete x->branch;
      x->branch = new ValSelectorSplit( x );
    }
  }
}

void Solver::setRandSplit(BuildObject **bvar, const int l) {
  VariableInt *x;
  BuildObject *bv;
  for(int i=0; i<l; ++i) {
    bv = bvar[i];
    if(bv != bv->getBuildObject())
      bv = bv->getBuildObject();
    if( bv->isSearchable() &&
	((x = bv->getVariable()) != NULL)) {
      delete x->branch;
      x->branch = new ValSelectorRandomSplit( x );
    }
  }
}

void Solver::setRandMinMax(BuildObject **bvar, const int l) {
  VariableInt *x;
  BuildObject *bv;
  for(int i=0; i<l; ++i) {
    bv = bvar[i];
    if(bv != bv->getBuildObject())
      bv = bv->getBuildObject();
    if( bv->isSearchable() &&
	((x = bv->getVariable()) != NULL)) {
      delete x->branch;
      x->branch = new ValSelectorRandMinMax( x );
    }
  }
}

void Solver::setRandomValueOrdering() {
  int i = length;
  while( i-- )
    if( !variables[i]->isGround() && !variables[i]->branch ) {
      if( variables[i]->getType() == VariableInt::RANGE ) {
	variables[i]->branch = new ValSelectorRandMinMax( variables[i] );
      } else if( variables[i]->getType() != VariableInt::INTLIST ) {
	variables[i]->branch = new ValSelectorRand( variables[i] );
      } else {
	variables[i]->branch = new ValSelectorFirst( variables[i] );
      }
    }
}

void Solver::reset_stats() {
  SOLUTIONS = NODES = BACKTRACKS = CHECKS
    = FAILURES = PROPAGS = FILTERINGS = RESTORES = MISC = 0;
  status = UNKNOWN;
  SOLTIME = 0.0;
  STARTTIME = getRunTime();
  BTSLIST.clear();
}

void Solver::reset_trail(const bool full) {
  backtrackTo(init_level-full);
  if(full) unaryCons.size = unaryCons_size;
}

void Solver::reset(const bool full) {
  reset_trail(full);
  //status = UNKNOWN;
  //if(full) 
  reset_stats();
}

void Solver::setLowerBounds(VarArray& vars, int* lb) {
  int i, n=vars.size();
  BuildObject **bvar = new BuildObject*[n];
  for(i=0; i<n; ++i) {
    bvar[i] = vars[i].var_ptr_;
  }
  setLowerBounds(bvar, n, lb);
}

void Solver::setUpperBounds(VarArray& vars, int* ub) {
  int i, n=vars.size();
  BuildObject **bvar = new BuildObject*[n];
  for(i=0; i<n; ++i) {
    bvar[i] = vars[i].var_ptr_;
  }
  setUpperBounds(bvar, n, ub);
}

//void Solver::clearUnaryConstraints()

void Solver::setLowerBounds(BuildObject **bvar, const int l, int* lb) {

  //   /// map
  //   int mapping[numvars];
  //   std::fill(mapping, mapping+numvars, -1);
  //   for(int j=0; j<unaryCons.size; ++j)
  //     if(!strcmp(unaryCons[j]->name, "lower_bound"))
  //       mapping[unaryCons[j]->X->id] = j;

  //   VariableInt *x;
  //   BuildObject *bv;
  //   for(int i=0; i<l; ++i) {
  //     bv = bvar[i];
  //     if(bv != bv->getBuildObject())
  //       bv = bv->getBuildObject();
  //     if( bv->isSearchable() && ((x = bv->getVariable()) != NULL) ) {
  //       if(mapping[x->id] >= 0) {
  // 	((UnaryConstraintMore*)(unaryCons[mapping[x->id]]))->bound = lb[i];
  //       } else {
  // 	unaryCons.push( new UnaryConstraintMore(this, x, lb[i]) );
  //       }
  //     }
  //   }

  VariableInt *x;
  BuildObject *bv;
  for(int i=0; i<l; ++i) {
    bv = bvar[i]->getBuildObject();
    if( bv->isSearchable() && ((x = bv->getVariable()) != NULL) ) {

      //      if(x->id == -41 ) {
      // 	x->print(std::cout);
      // 	std::cout << " >= " << lb[i] << std::endl;
      //       }

      SimpleUnaryConstraint cons('l', lb[i], x);
      sUnaryCons.push(cons);
    }
  }

  //   std::cout << "SET LOWER BOUND AT LEVEL " << level << " ";
  //   variables[41]->print(std::cout);
  //   std::cout << std::endl;

}

void Solver::setUpperBounds(BuildObject **bvar, const int l, int* ub) {

  /// map
  int mapping[numvars];
  std::fill(mapping, mapping+numvars, -1);
  for(int j=0; j<unaryCons.size; ++j)
    if(!strcmp(unaryCons[j]->name, "upper_bound"))
      mapping[unaryCons[j]->X->id] = j;

  VariableInt *x;
  BuildObject *bv;
  for(int i=0; i<l; ++i) {
    bv = bvar[i];
    if(bv != bv->getBuildObject())
      bv = bv->getBuildObject();
    if( bv->isSearchable() && ((x = bv->getVariable()) != NULL) ) {
      //x->setMax(ub[i]);
      if(mapping[x->id] >= 0) {
	((UnaryConstraintLess*)(unaryCons[mapping[x->id]]))->bound = ub[i];
      } else {
	unaryCons.push( new UnaryConstraintLess(this, x, ub[i]) );
      }
    }
  }
}


void Solver::setGuidedOrdering(VarArray& scope, int* ideal, const char*planb) {
  BuildObject *bvar[scope.size()];
  int i, l = scope.size(),pb=0;
  for(i=0; i<l; ++i)
    bvar[i] = scope[i].var_ptr_;
  if(!strcmp(planb,"std")) pb = 0;
  else if(!strcmp(planb,"2nd")) pb = 1;
  else if(!strcmp(planb,"spl")) pb = 2;
  else if(!strcmp(planb,"nbd")) pb = 3;
  setGuidedOrdering(bvar, l, ideal, pb);
}

void Solver::setGuidedOrdering(BuildObject **bvar, const int l, int* ideal, const int pb) {
  VariableInt *x;
  BuildObject *bv;

  //  std::cout << "GUIDE2" <<std::endl;

  for(int i=0; i<l; ++i) {

    bv = bvar[i];
    if(bv != bv->getBuildObject())
      bv = bv->getBuildObject();

    //     bv->print(std::cout);
    //     std::cout << std::endl;

    if( bv->isSearchable() &&
	((x = bv->getVariable()) != NULL) // &&
	// 	(x->getType() == VariableInt::BOOL)
	) {

      //       std::cout << '        ' << std::endl;
      //       bv->print(std::cout);
      //       std::cout << std::endl;

      delete x->branch;
      x->branch = new ValSelectorGuided( x, ideal[i], pb );
    }
  }
}

void Solver::setGuidedSplitOrdering(BuildObject **bvar, const int l, int* ideal, int* proba) {
  VariableInt *x;
  BuildObject *bv;

  //std::cout << "GUIDE" << std::endl;

  for(int i=0; i<l; ++i) {
    bv = bvar[i];


//     if(!bv) {
//       std::cout << "\n\n **** build object is null **** \n\n" ;
//     }


    if(bv && bv != bv->getBuildObject())
      bv = bv->getBuildObject();

    //std::cout << (bvar[i]->isSearchable()) << " " << (bvar[i]->getVariable()) << std::endl;


    if( bv->isSearchable() &&
	((x = bv->getVariable()) != NULL) ) {

 //      //std::cout << "guide" << std::endl;

//       if(x->getType() == VariableInt::VIRTUAL) {
// 	std::cout << "\n\n **** build object is virtual **** \n\n" ;
//       }


//      if(x->getType() == VariableInt::CONST) {
//        std::cout << "\n\n **** build object is const **** \n\n" ;
//      }

      if(x->getType() != VariableInt::VIRTUAL && x->getType() != VariableInt::CONST) {
	
	delete x->branch;
	x->branch = new ValSelectorRandGuidedSplit( x, ideal[i], proba[i] );
      }
    }
  }
}

void Solver::setRandGuidedOrdering(VarArray& scope, int* ideal, int* proba, int* range) {
  BuildObject *bvar[scope.size()];
  int i, l = scope.size();
  for(i=0; i<l; ++i)
    bvar[i] = scope[i].var_ptr_;
  setRandGuidedOrdering(bvar, l, ideal, proba, range);
}

void Solver::setRandGuidedOrdering(BuildObject **bvar, const int l, int* ideal, int* proba, int* range) {
  //BuildObject *bvar;
  VariableInt *x;
  BuildObject *bv;

  for(int i=0; i<l; ++i) {

    bv = bvar[i];
    if(bv != bv->getBuildObject())
      bv = bv->getBuildObject();

    //bvar = scope[i].var_ptr_;
    if( bv->isSearchable() &&
	((x = bv->getVariable()) != NULL) // &&
	// 	(x->getType() == VariableInt::BOOL)
	)
      x->branch = new ValSelectorRandGuided( x, ideal[i], proba[i], range[i] );
  }
}
