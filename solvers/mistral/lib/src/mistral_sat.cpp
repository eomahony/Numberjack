
#include <mistral_sat.h>
#include <mistral_mod.h>

#include <fstream>
#include <stdio.h>
#include <stdlib.h>
 

double *sorting_array;
int compar(const void *a, const void *b)
{
  if (sorting_array[*(int*)a]==sorting_array[*(int*)b])
    return 0;
  else
    if (sorting_array[*(int*)a] < sorting_array[*(int*)b])
      return 1;
    else
      return -1;
}
void initSort(double *sa)
{
  sorting_array = sa;
}

int atom( const Mistral::Literal l )
{
  return( l < 0 ? -l : l );
}

using namespace Mistral;
using namespace std;

// int value( const Atom a )
// {
//   return( polarity[a] == a );
// }

// int literal( const Atom a )
// {
//   return polarity[a];
// }

// int ground( const Atom a )
// {
//   return assumptions.member(a);
// }

void SatSolver::initParams() {

  numAtoms = 0;

  isWatchedBy = NULL;
  polarity = NULL;
  activity = NULL;
  reason = NULL;
  lvl = NULL;
  restart_policy = NULL;

  stats.unit_props = 0;
  stats.nodes = 0;
  stats.conflicts = 0;
  stats.base_avg_size = 0;
  stats.learnt_avg_size = 0;
  stats.literals = 0;

  params.verbosity = 4;
  params.time_limit = -1;
  params.seed = 11041979;
  params.policy = GEOMETRIC;
  params.restart_base = 200;
  params.restart_limit = 200;
  params.restart_factor = 1.05;
  params.activity_increment = 1e-2;
  params.forgetfulness = .75;
  params.randomization = 2;
  params.shuffle = true;
  params.decay = 0.96;

  restart_policy = NULL;
}

void SatSolver::setPolicy( const int policy )
{
  if( policy == LUBY )
    restart_policy = new Luby(this);
  else
    restart_policy = new Geometric(this);
}

SatSolver::SatSolver() 
{
  initParams();
}

SatSolver::SatSolver(const char* filename) 
{
  initParams();
  parseDimacs(filename);
}

SatSolver::SatSolver(CSP& model) 
{
  initParams();
  model.build( this );
  numClauses = base.size;
}

SatSolver::~SatSolver() 
{
  isWatchedBy -= (numAtoms+1);
  delete [] isWatchedBy;
  delete [] polarity;
  activity -= (numAtoms+1);
  delete [] activity;
  delete [] reason;
  delete [] lvl;
  for(int i=0; i<learnt.size; ++i)
    free(learnt[i]);
  for(int i=0; i<base.size; ++i)
    free(base[i]);
  for(int i=0; i<original.size; ++i)
    free(original[i]);
  delete restart_policy;
}


RestartPolicy::RestartPolicy(SatSolver* s)
  : restart_base(s->params.restart_base), 
    restart_limit(s->params.restart_limit), 
    restart_factor(s->params.restart_factor)
{
}

RestartPolicy::~RestartPolicy() {}


Geometric::Geometric(SatSolver* s)
  : RestartPolicy(s)
{
}

Geometric::~Geometric() {}

void Geometric::update()
{
  restart_base = (int)((double)restart_base * (double)restart_factor);
  restart_limit += restart_base;
}


Luby::Luby(SatSolver* s)
  : RestartPolicy(s)
{
  iteration = 0;
}

Luby::~Luby() {}

unsigned int Luby::binlog( unsigned int x ) {
  if( !x ) return NOVAL;
  if( x == 1 ) return 0;
  unsigned int exponent = 0;
  while( (x >> (++exponent)) > 1 );
  return exponent;
}

unsigned int Luby::luby_seq( unsigned int iteration ) {
  unsigned int thelog = binlog( iteration );
  if( iteration == (unsigned int)((1 << (thelog + 1))-1) )
    return (1 << thelog);
  return luby_seq( iteration - (1 << thelog) + 1 );
}

void Luby::update()
{
  restart_limit += (restart_base * luby_seq(++iteration));  
}


int SatSolver::solve()
{
  //srand(params.seed);
  usrand(params.seed);
  stats.start_time = getRunTime();
  if(!restart_policy)
    setPolicy( params.policy );
  if(params.verbosity>1)
    cout << endl 
	 << "c  ==================================[ Mistral (Sat module) ]===================================" << endl
	 << "c  |      SEARCH  STATISTICS           |                  PROBLEM STATISTICS                   |" << endl
	 << "c  |  Conflicts |       Nodes  CPUTime |  Atoms  Clauses (Size)   Learnt (Size)   Added (Size) |" << endl
	 << "c  =============================================================================================" << endl;
  int result = (unsigned int)(stats.base_avg_size)+1;
  while(1) {
    if(params.shuffle) shuffle();
    result = iterative_search();
    if(restart_policy) restart_policy->update();
    if(params.verbosity>1)
      cout << "c  | " << setw(10) << stats.conflicts
	   << " | " << setw(10) << stats.nodes
	   << " " << setw(9) << (getRunTime() - stats.start_time)
	   << " | " << setw(6) << (numAtoms - stats.literals)
	   << " " << setw(8) << base.size
	   << " " << setw(5) << setprecision(3) << stats.base_avg_size     
	   << "  " << setw(8) << learnt.size
	   << " " << setw(5) << setprecision(3) << stats.learnt_avg_size ;
    forget();
    if(params.verbosity>1)
      cout << "  " << setw(7) << learnt.size
	   << " " << setw(5) << setprecision(3) << stats.learnt_avg_size 
	   << "  |" << endl;
    if(result == LIMITOUT)
      status = UNKNOWN;
    else break;
  } 
  
  cout << "c  =============================================================================================" << endl;
  if( result == SAT ) 
    cout << "c  |                                    SATISFIABLE!                                           |" << endl;
  else if( result == UNSAT ) 
    cout << "c  |                                     UNSATISFIABLE!                                        |" << endl;
  else if( result == UNKNOWN ) 
    cout << "c  |                                       UNKNOWN!                                            |" << endl;
  double cputime = (getRunTime() - stats.start_time);
  if( cputime < 0.001 ) cputime = 0.001;
  cout << "c  =============================================================================================" << endl ;
  if(params.verbosity>0)    
    cout << "c  |      Conflicts  |          Nodes  |    CPU Time  | Conflicts/s |     Nodes/s | node/confl |" << endl
	 << "c  | " << setw(14) << stats.conflicts 
	 << "  | " << setw(14) << stats.nodes 
	 << "  | " << setw(11) << cputime
	 << "  | " << setw(11) << int(double(stats.conflicts)/cputime) 
	 << " | " << setw(11) << int(double(stats.nodes)/cputime) 
	 << " | " << setw(10) << (double(stats.nodes)/double(stats.conflicts))
	 << " |" << endl
	 << "c  =============================================================================================" << endl << endl;
  if(params.verbosity > 3) {
    if(result == SAT)
      {
	cout << "s SATISFIABLE\nv";
	for(int i=1; i<=numAtoms; ++i) 
	  cout << " " << (polarity[i] == i ? 1 : 0) ;
	cout << endl;
      } else cout << "s UNSATISFIABLE" << endl;
    cout << "d ASSIGNMENTS " << stats.nodes << endl
	 << "d FAILS " << stats.conflicts << endl
	 << "d BACKTRACKS " << stats.conflicts << endl
	 << "d NODES/s " << int(double(stats.nodes)/cputime) << endl;
  }


//   if( result == SAT )
//     printDecisions(cout);

  return status;
}

void SatSolver::printAll(ostream& o) const
{
  printDecisions( o );
  printClauses( o );
  //printWatchers( o );
}

void SatSolver::printWatchers(ostream& o, int beg, int end) const
{
  if( beg == NOVAL ) beg = -numAtoms;
  if( end == NOVAL ) end = numAtoms;
  for(int i=beg; i<=end; ++i)
    {
      o << i << " is watched by";
      for(int j=0; j<isWatchedBy[i].size; ++j) {
	assert(isWatchedBy[i][j]);
	Clause& clause = *(isWatchedBy[i][j]);
	printClause(cout, isWatchedBy[i][j]);
	assert( clause[0] == i || 
		clause[1] == i );
	assert( clause[0] != clause[1] );
      }
      o << endl;
    }
}

void SatSolver::printDecisions(ostream& o, bool mode) const
{
  if( mode ) {
    cout << "c literals: ";
    for(unsigned int i=0; i<assumptions.size; ++i)
      o << " " <<setw(3)<< polarity[assumptions[i]];
    o << endl;
    cout << "c   level: ";
    for(unsigned int i=0; i<assumptions.size; ++i)
      o << " " << setw(3) << lvl[assumptions[i]];
    o << endl;
  } else {
    for(unsigned int i=0; i<assumptions.size; ++i) {
      o << lvl[assumptions[i]] << "\t" << polarity[assumptions[i]] << "\t";
      if( reason[assumptions[i]] )
	printClause(o, reason[assumptions[i]]);
      else 
	cout << "decision";
      cout << endl;
    }
  }
}

void SatSolver::printClause(ostream& o, Clause *cl) const
{
  Clause& clause = *cl;
  o //<< " " << cl 
    << "(";
  for(int i=0; i<clause.size-1; ++i)
    o << clause[i] << " ";
  o << clause[clause.size-1] << ")";
}

void SatSolver::printClauses(ostream& o) const
{
  //   int cl[numAtoms], cs, flag, a;
  //   for(int i=0; i<clauses.size; ++i)
  //     {
  //       Clause& clause = *(clauses[i]);
  //       cs = 0;
  //       flag = 1;
  //       for(int j=0; j<clause.size; ++j)
  // 	{
  // 	  a = atom( clause[j] );
  // 	  if( !assumptions.member(a) )
  // 	    cl[cs++] = clause[j];
  // 	  else if( polarity[a] == clause[j] ) {
  // 	    flag = 0;
  // 	    break;
  // 	  }
  // 	}
  //       if( flag ) {
  // 	cout << "c" << i << " (";
  // 	for(int j=0; j<cs; ++j)
  // 	  o << " " << cl[j] ;
  // 	o << ")" << endl;
  //       }
  //     }
  //   o << endl;

  o << "base (" << base.size << "):" << endl;
  for(int i=0; i<base.size; ++i)
    {
      o << "c" << i ;
      printClause( o, base[i] );
      o << endl;
    }
  o << endl;
  o << "learnt (" << learnt.size << "):" << endl;
  for(int i=0; i<learnt.size; ++i)
    {
      o << "c" << i ;
      printClause( o, learnt[i] );
      o << endl;
    }
  o << endl;
}

void SatSolver::init(const int n, const int m)
{
  status = UNKNOWN;
  numAtoms = n;
  numClauses = m;
  nextDeduction = 0;

  visited.init(1, numAtoms, BitSet::empt);

  decisions.init(0, numAtoms);
  assumptions.init(1, numAtoms);
  learnt.init(0, m);

  isWatchedBy = new Vector<Clause*>[2*numAtoms+3];
  isWatchedBy += (numAtoms+1);

  polarity = new int[numAtoms+2];
  activity = new double[2*numAtoms+3];
  activity += (numAtoms+1);
  reason = new Clause*[numAtoms+2];
  lvl = new int[numAtoms+2];
  for(int i=0; i<=numAtoms+1; ++i) {
    polarity[i] = i;
    activity[i] = 0;
    activity[-i] = 0;
    reason[i] = NULL;
    lvl[i] = numAtoms+2;
    if(i) {
      isWatchedBy[i].init(0,512);
      isWatchedBy[-i].init(0,512);
    }
  }
  if( base.empty() )
    base.init(0, m);
  else
    assert( base.size == m );
  int i, j=base.size;
  for(i=0; i<j; ++i) {
    Clause& clause = *(base[i]);

//     std::cout << clause.size << " " << numAtoms << std::endl;

//     std::cout << clause[0] << std::endl;
//     std::cout << clause[1] << std::endl;
    
//     //std::cout << (isWatchedBy+clause[0]) << std::endl;
//     //std::cout << (isWatchedBy+clause[1]) << std::endl;

    isWatchedBy[clause[0]].push(base[i]);
    isWatchedBy[clause[1]].push(base[i]);

//     std::cout << "here" << std::endl;
  }
} 

void SatSolver::parseDimacs(const char* filename) 
{
  unsigned int LARGENUMBER = 131072;
  ifstream infile( filename );
  char c=' ';
  string word;
  int N, M, l=0;

  // skip comments
  infile >> c;
  while( c != 'p' ) {
    infile.ignore( LARGENUMBER, '\n' );
    infile >> c;
  }

  infile >> word;
  assert( word == "cnf" );
  
  // get number of atoms and clauses
  infile >> N;
  infile >> M;

  init(N, M);
  learnt_clause.init(0, N);

  for(int i=0; i<M; ++i)
    {
      learnt_clause.clear();
      do {
	infile >> l;
	learnt_clause.push(l);
      } while(l && infile.good());
      learnt_clause.pop();
      addClause( base, learnt_clause, stats.base_avg_size );
      addOriginalClause( learnt_clause );
    }
}

int SatSolver::checkSolution()
{
  Atom x;
  bool correct = true;
  for(int i=0; correct && i<numClauses; ++i) {
    bool satisfied = false;
    Clause& clause = *(original[i]);

    for(int j=0; !satisfied && j<clause.size; ++j) {
      x = atom(clause[j]);
      satisfied = (polarity[x] == clause[j] && assumptions.member( x ) );
    }
    if(!satisfied) {
      printClause(cerr, original[i]);
      cerr << endl;
    }

    correct = satisfied;
  }
  if(!correct) {
    cerr << "/!\\ The solution is not correct /!\\" << endl; 
    return UNKNOWN;
  }
  return SAT;
}

