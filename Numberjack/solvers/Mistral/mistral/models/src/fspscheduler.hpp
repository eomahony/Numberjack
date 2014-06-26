 
#include <mistral_sol.h>
#include <iostream>
#include <fstream>

using namespace std;


int bufsize = 10000;

int **duration=NULL, **machine=NULL, **family=NULL, nMachines, nJobs, nFamilies,
  lb, opt, *jsp_duration=NULL, *jsp_machine=NULL, 
  ***setup_time=NULL, *max_setup=NULL, **family_matrix=NULL,
  **release_date=NULL, **due_date=NULL, **time_lag[2] = {NULL, NULL};
int *best_solution=NULL, *zeros=NULL, *ones=NULL, *probability=NULL, *range=NULL;
vector< int* > solutions;
vector< int* > disjunct_weights;
vector< int* > normalised_weights;
vector< int > min_weights;
vector< int > max_weights;
int *disjunct_index=NULL;
int Tprob = (1 << 24);



static const int LEX     =  0;
static const int PROMISE =  1;
static const int ANTI    = -1;
static const int GUIDED  =  2;
static const int RGUIDED =  3;
static const int RAND    =  4;

int init_ub;
int init_lb;
int maxfsble;
int minfsble;
int makespan;
int max_infeasible;
int ndisjuncts=0;

int nogood_size=0;
int nb_nogood=0;


long unsigned int total_solutions     = 0;
double            avg_distance        = 0;
int               min_distance        = NOVAL;
int               max_distance        = 0;

int        num_weight_updates  = 0;
double avg_disjunct_avg_weight = 0;
double avg_disjunct_min_weight = 0;
double avg_disjunct_max_weight = 0;
double min_disjunct_avg_weight = NOVAL;
int min_disjunct_min_weight = NOVAL;
int min_disjunct_max_weight = NOVAL;
double max_disjunct_avg_weight = 0;
int max_disjunct_min_weight = 0;
int max_disjunct_max_weight = 0;

double avg_task_avg_weight = 0;
double avg_task_min_weight = 0;
double avg_task_max_weight = 0;
double  min_task_avg_weight = NOVAL;
int min_task_min_weight = NOVAL;
int min_task_max_weight = NOVAL;
double max_task_avg_weight = 0;
int max_task_min_weight = 0;
int max_task_max_weight = 0;

int total_restarts = 0;
int first_makespan = 0;

long unsigned int total_nodes     = 0;
long unsigned int total_bts       = 0;
long unsigned int total_fails     = 0;
long unsigned int total_propags   = 0;
double opt_time = 0.0;
double dsrun_time = 0.0;
double proof_time = 0.0;
double total_time = 0.0;


int Seed, Cutoff, Dichotomy, Base, PIterations, PLimit, Probe, UpdateW, UpdateI, 
		  SAC, Randomized, Verbose, Step, Weight, Optimise, Rngd, Gap, 
		  Reset, VO, D_VO, I_VO, JobRank, Wprof, Multiplier, PolicyRestart, 
		  Solved, Result, Pool, Rprob, UBinit;
double Factor, Decay, Skew;
string Policy, Heuristic, Type, Value, Algo, DValue, IValue;

const int nia = 19;
const char* int_ident[nia] = {"-seed", "-cutoff", "-dichotomy", "-base", "-sac", 
			      "-randomized", "-verbose", "-step", "-weight", "-optimise", 
			      "-restart_ngd", "-gap", "-reset", "-wprof", "-mult", 
			      "-jobrank", "-pool", "-rprob", "-ub"};
int int_param[nia];

const int nsa = 10;
const char* str_ident[nsa] = {"-heuristic", "-restart", "-factor", "-decay", "-type", "-value", "-algo", "-dvalue", "-ivalue", "-skew"};
const char* str_param[nsa];


using namespace Mistral;

void addHeuristic( Solver& s, string Heu, const int rdz, const int val_ord ) {

  //   if(ValueO == "lex") {
  //     for(int i=0; i<s.length; ++i)
  //       s.sequence[i]->branch = new ValSelectorMin( s.sequence[i] );
  //   } else if(ValueO == "rand") {
  //     for(int i=0; i<s.length; ++i)
  //       s.sequence[i]->branch = new ValSelectorRand( s.sequence[i] );
  //   } else if(ValueO == "antilex") {
  //     for(int i=0; i<s.length; ++i)
  //       s.sequence[i]->branch = new ValSelectorMax( s.sequence[i] );
  //   } 

  
  if( Heu == "dom" ) {
    MinDomain h(abs(rdz));
    s.add( h );
  }
  if( Heu == "lex" ) {
    Lexicographic h;
    s.add( h );
  } 
  else if( Heu == "deg") {
    MaxDegree h(abs(rdz));
    s.add( h );
  } 
  else if( Heu == "rand") {
    Random h;
    s.add( h );
  } 
  else if( Heu == "dom+deg") {
    MinDomMaxDeg h(abs(rdz));
    s.add( h );
  } 
  else if( Heu == "domodeg") {
    DomOverDeg h(abs(rdz));
    s.add( h );
  } 
  else if( Heu == "domowldeg") {
    DomOverWLDeg h(abs(rdz));
    s.add( h );
  }
  else if( Heu == "domowdeg") {
    DomOverWDeg h(abs(rdz));
    s.add( h );
  }
  else if( Heu == "neighbor") {
    Neighbor h(abs(rdz));
    s.add( h );
  } 
  else if( Heu == "impact") {
    Impact h(abs(rdz));
    s.add( h );
  }
  else if( Heu == "impactodeg") {
    ImpactOverDeg h(abs(rdz));
    s.add( h );
  }
  else if( Heu == "impactowdeg") {
    ImpactOverWDeg h(abs(rdz));
    s.add( h );
  }
  else if( Heu == "impactowldeg") {
    ImpactOverWLDeg h(abs(rdz));
    s.add( h );
  }
  else if( Heu == "pfsp") {
    PFSP h(abs(rdz), val_ord);
    s.add( h );
  }
  else if( Heu == "osp") {
    OSP h(abs(rdz), val_ord);
    s.add( h );
  }
  else if( Heu == "osp-b") {
    OSP h(abs(rdz), val_ord, OSP::DOM_O_BOOLWEIGHT);
    s.add( h );
  }
  else if( Heu == "osp-t") {
    OSP h(abs(rdz), val_ord, OSP::DOM_O_TASKWEIGHT);
    s.add( h );
  }
  else if( Heu == "osp-bt") {
    OSP h(abs(rdz), val_ord, OSP::DOM_O_BOOLTASKWEIGHT);
    s.add( h );
  }
  else if( Heu == "dom-t") {
	  OSP h(abs(rdz), val_ord, OSP::TASKDOM);
	  s.add( h );
  }
  else {
    NoOrder h;
    s.add( h );
  }
}


// int greedy_jsp() {
//   int i;

//   // those are the tasks that we can currently execute
//   int nopen=nJobs, *open_task[nJobs];
//   for(i=0; i<nJobs; ++i) 
//     open_task[i] = i*nMachines;
  
//   while(nopen) // as long as at least one task remains
//     {
//       // select the task with least due date;
      
//     }

// }



void tsp_readData( char* filename )
{
  int i, j;
  string tag;
  ifstream infile( filename, ios_base::in );
  
  infile >> nJobs;
  nMachines = 1;

  jsp_duration = new int[nJobs*nMachines];
  jsp_machine = new int[nJobs*nMachines];

  duration = new int*[nJobs];
  machine = new int*[nJobs];
  for(i=0; i<nJobs; ++i) {
    duration[i] = new int[nMachines];
    machine[i] = new int[nMachines];
  }

  max_setup = new int[nMachines];
  setup_time = new int**[nMachines];
  for(i=0; i<nMachines; ++i) {
    max_setup[i] = 0;
    setup_time[i] = new int*[nJobs];
    for(j=0; j<nJobs; ++j)
      setup_time[i][j] = new int[nJobs];
  }

  double X[nJobs];
  double Y[nJobs];
  double demand[nJobs];
  double D;

  due_date = new int*[nJobs];
  release_date = new int*[nJobs];

  for(i=0; i<nJobs; ++i) {
    due_date[i] = new int[1];
    release_date[i] = new int[1];

    infile >> j;

    //cout << j << endl;

    //assert( j == i+1 );
    
    infile >> X[i];
    X[i] *= Multiplier;
    infile >> Y[i];
    Y[i] *= Multiplier;
    infile >> demand[i];
    infile >> D;
    release_date[i][0] = (int)(D*Multiplier);
    infile >> D;
    due_date[i][0] = (int)(D*Multiplier);
    infile >> D;
    duration[i][0] = (int)(D*Multiplier);

  }


  for(i=1; i<nJobs; ++i) {
    jsp_duration[i] = duration[i][0];
    machine[i][0] = 0;
    jsp_machine[i] = 0;
    for(j=0; j<i; ++j) {
      setup_time[0][i][j] = setup_time[0][j][i] = 
	(int)(sqrt((X[i] - X[j])*(X[i] - X[j]) + (Y[i] - Y[j])*(Y[i] - Y[j])));
      if(max_setup[0] < setup_time[0][i][j]) max_setup[0] = setup_time[0][i][j];
      cout << j << "." << i << ": " << setup_time[0][i][j] << " ";
    }
    cout << endl;
  }

  ndisjuncts = nMachines * (nJobs * (nJobs-1) / 2);

}



void jsp_readData( char* filename )
{
  int i, j;
  long int dump;
  string tag;
  ifstream infile( filename, ios_base::in );
  
  infile >> nJobs;
  infile >> nMachines;

  jsp_duration = new int[nJobs*nMachines];
  jsp_machine = new int[nJobs*nMachines];

  duration = new int*[nJobs];
  machine = new int*[nJobs];
  for(i=0; i<nJobs; ++i) {
    duration[i] = new int[nMachines];
    machine[i] = new int[nMachines];
  }

  infile >> dump;
  infile >> dump;

  infile >> opt ;
  infile >> lb;

  infile >> tag;
  assert( tag == "Times" );

  for(i=0; i<nJobs; ++i) 
    for(j=0; j<nMachines; ++j) {
      infile >> duration[i][j];
      jsp_duration[i*nMachines+j] = duration[i][j];
    }

  infile >> tag;
  assert( tag == "Machines" );

  for(i=0; i<nJobs; ++i) 
    for(j=0; j<nMachines; ++j) {
      infile >> machine[i][j];
      --machine[i][j];
      jsp_machine[i*nMachines+j] = machine[i][j];
    }


  //   setup_time = new int**[nMachines];
  //   for(int k=0; k<nMachines; ++k) {
  //     setup_time[k] = new int*[nJobs];
  //     for(i=0; i<nJobs; ++i) {
  //       setup_time[k][i] = new int[nJobs];
  //       for(j=0; j<nJobs; ++j)
  // 	setup_time[k][i][j] = randint(50);
  //     }
  //   }

  ndisjuncts = nMachines * (nJobs * (nJobs-1) / 2);

}

void fsp_readData( char* filename )
{
  int i, j;
  long int dump;
  string tag;
  ifstream infile( filename, ios_base::in );
  
  infile.ignore(1000, '\n');

  infile >> nJobs;
  infile >> nMachines;

  infile >> dump;
  infile >> opt ;
  infile >> lb;



  infile >> tag;
  assert( tag == "processing" );

  infile >> tag;
  assert( tag == "times" );

  infile >> tag;
  assert( tag == ":" );

  //infile.ignore(1000, '\n');


  jsp_duration = new int[nJobs*nMachines];
  jsp_machine = new int[nJobs*nMachines];

  duration = new int*[nJobs];
  machine = new int*[nJobs];
  for(i=0; i<nJobs; ++i) {
    duration[i] = new int[nMachines];
    machine[i] = new int[nMachines];
  }

  for(j=0; j<nMachines; ++j) 
    for(i=0; i<nJobs; ++i) {
      infile >> duration[i][j];
      jsp_duration[i*nMachines+j] = duration[i][j];
      machine[i][j] = j;
      jsp_machine[i*nMachines+j] = machine[i][j];
    }

  ndisjuncts = (Type != "pfsp" ? nMachines : 1) * (nJobs * (nJobs-1) / 2);

}


void jtl_readData( char* filename )
{
  int i, j;
  string tag;
  ifstream infile( filename, ios_base::in );
  
  infile >> nJobs;
  infile >> nMachines;

  jsp_duration = new int[nJobs*nMachines];
  jsp_machine = new int[nJobs*nMachines];

  duration = new int*[nJobs];
  machine = new int*[nJobs];
  for(i=0; i<nJobs; ++i) {
    duration[i] = new int[nMachines];
    machine[i] = new int[nMachines];
  }

  for(i=0; i<nJobs; ++i) {
    for(j=0; j<nMachines; ++j) {

      infile >> machine[i][j];
      jsp_machine[i*nMachines+j] = machine[i][j];

      infile >> duration[i][j];
      jsp_duration[i*nMachines+j] = duration[i][j];
    }
  }


  infile >> tag;
  infile >> opt;
  infile.ignore( 100, '\n' );

  char c;
  infile.get(c);
  assert( c == 'T' );
  infile.get(c);
  assert( c == 'L' );
  infile.get(c);
  assert( c == '=' );

  time_lag[0] = new int*[nJobs];
  time_lag[1] = new int*[nJobs];

  for(i=0; i<nJobs; ++i) {
    time_lag[0][i] = new int[nMachines];
    time_lag[1][i] = new int[nMachines];
    for(j=0; j<nMachines; ++j) {
      infile >> time_lag[0][i][j];
      infile >> time_lag[1][i][j];
    }
  }

  //   setup_time = new int**[nMachines];
  //   for(int k=0; k<nMachines; ++k) {
  //     setup_time[k] = new int*[nJobs];
  //     for(i=0; i<nJobs; ++i) {
  //       setup_time[k][i] = new int[nJobs];
  //       for(j=0; j<nJobs; ++j)
  // 	setup_time[k][i][j] = randint(50);
  //     }
  //   }

  ndisjuncts = nMachines * (nJobs * (nJobs-1) / 2);
}


void jla_readData( char* filename )
{
  int i, j;
  string tag;
  ifstream infile( filename, ios_base::in );
  
  infile >> nJobs;
  infile >> nMachines;

  jsp_duration = new int[nJobs*nMachines];
  jsp_machine = new int[nJobs*nMachines];

  duration = new int*[nJobs];
  machine = new int*[nJobs];
  for(i=0; i<nJobs; ++i) {
    duration[i] = new int[nMachines];
    machine[i] = new int[nMachines];
  }

  for(i=0; i<nJobs; ++i) {
    for(j=0; j<nMachines; ++j) {

      infile >> machine[i][j];
      jsp_machine[i*nMachines+j] = machine[i][j];

      infile >> duration[i][j];
      jsp_duration[i*nMachines+j] = duration[i][j];
    }
  }

  ndisjuncts = nMachines * (nJobs * (nJobs-1) / 2);
}


void sds_readData( char* filename )
{
  int i, j;
  string tag;
  ifstream infile( filename, ios_base::in );
  
  infile >> nMachines;
  infile >> nJobs;
  infile >> nFamilies;

  jsp_duration = new int[nJobs*nMachines];
  jsp_machine = new int[nJobs*nMachines];
  
  family_matrix = new int*[nFamilies+1];
  for(i=0; i<=nFamilies; ++i)
    family_matrix[i] = new int[nFamilies];

  max_setup = new int[nMachines];
  setup_time = new int**[nMachines];
  for(i=0; i<nMachines; ++i) {
    max_setup[i] = 0;
    setup_time[i] = new int*[nJobs];
    for(j=0; j<nJobs; ++j)
      setup_time[i][j] = new int[nJobs];
  }

  family = new int*[nJobs];
  duration = new int*[nJobs];
  machine = new int*[nJobs];
  for(i=0; i<nJobs; ++i) {
    duration[i] = new int[nMachines];
    machine[i] = new int[nMachines];
    family[i] = new int[nMachines];
  }

  //infile >> dump;
  //infile >> dump;

  //infile >> opt ;
  //infile >> lb;

  //infile >> tag;
  //assert( tag == "Times" );

  for(i=0; i<nJobs; ++i) {
    
    infile >> j;
    assert(j==nMachines);

    for(j=0; j<nMachines; ++j) {
      infile >> duration[i][j];
      jsp_duration[i*nMachines+j] = duration[i][j];
      
      infile >> machine[i][j];
      --machine[i][j];
      jsp_machine[i*nMachines+j] = machine[i][j];

      infile >> family[i][j];
      --family[i][j];

      //       cout << "(" << duration[i][j] << ","
      // 	   << (machine[i][j]+1) << ","
      // 	   << (family[i][j]+1) << ")" << endl;

    }
  }

  for(i=0; i<=nFamilies; ++i)
    for(j=0; j<nFamilies; ++j)
      infile >> family_matrix[i][j];

  for(int k=0; k<nMachines; ++k) {
    //    cout << "Machine " << k << endl;
    for(i=0; i<nJobs; ++i) 
      for(j=0; j<nJobs; ++j) {
	//	cout << "(" << i << "," << j << ") ";
	setup_time[k][i][j] = family_matrix[1+family[i][k]][family[j][k]] ;
	if(max_setup[k] < setup_time[k][i][j]) max_setup[k] = setup_time[k][i][j];
	//	cout << setup_time[k][i][j] << "  ";
      }
    //    cout << endl;
  }


  release_date = new int*[nJobs];
  for(i=0; i<nJobs; ++i) {
    release_date[i] = new int[nMachines];
    for(j=0; j<nMachines; ++j)
      release_date[i][j] = family_matrix[0][family[i][j]];
  }

  //   //   setup_time = new int**[nMachines];
  //   //   for(int k=0; k<nMachines; ++k) {
  //   //     setup_time[k] = new int*[nJobs];
  //   //     for(i=0; i<nJobs; ++i) {
  //   //       setup_time[k][i] = new int[nJobs];
  //   //       for(j=0; j<nJobs; ++j)
  //   // 	setup_time[k][i][j] = randint(50);
  //   //     }
  //   //   }

  ndisjuncts = nMachines * (nJobs * (nJobs-1) / 2);
}


void osp_readData( char* filename )
{
  char buf[bufsize];
  ifstream infile( filename, ios_base::in );

  do {
    infile.getline( buf, bufsize, '\n' );
  } while( buf[0] == '#' );

  int i=0; 
  while( buf[i] != ' ' ) ++i;
  buf[i] = '\0';
  lb = atoi( buf );

  while( buf[i] == '\0' || buf[i] == ' ' || buf[i] == '*' ) ++i;
  int j = i;
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
  infile.getline( buf, bufsize, '\n' );

  duration = new int*[nJobs];
  for(i=0; i<nJobs; ++i)
    duration[i] = new int[nMachines];

  do {
    infile.get( buf[0] );
    if( buf[0] != '#' ) break;
    infile.getline( buf, bufsize, '\n' );
  } while( true );
  infile.unget();

  for(i=0; i<nJobs; ++i) {
    for(j=0; j<nMachines; ++j) {
      infile >> duration[i][j];
    }
  }

  //ndisjuncts = nMachines * (nJobs * (nJobs-1) / 2);

  ndisjuncts = (nMachines * (nJobs * (nJobs-1) / 2) + 
		nJobs * (nMachines * (nMachines-1) / 2));
}


int osp_upperbound() {
  int i, j, mkp=0, task, boundJob[nJobs], boundMachine[nJobs]; 
  std::fill(boundJob, boundJob+nJobs, 0);
  std::fill(boundMachine, boundMachine+nJobs, 0);

  assert( nJobs == nMachines );

  for(i=0; i<nJobs; ++i)
    {
      for(j=0; j<nJobs; ++j)
	{
	  task = std::max( boundJob[j], boundMachine[(i+j) % nJobs] );
	  task += duration[j][(i+j) % nJobs];
	  boundJob[j] = task;
	  boundMachine[(i+j) % nJobs] = task;
	  if( task > mkp ) mkp = task;
	}
    }
  return mkp;
}


int osp_lowerbound()
{
  int i, j, mkp=0, boundJob[nJobs], boundMachine[nMachines]; 
  std::fill(boundJob, boundJob+nJobs, 0);
  std::fill(boundMachine, boundMachine+nMachines, 0);

  for(i=0; i<nJobs; ++i) {
    for(j=0; j<nMachines; ++j)
      {
	boundJob[i] += duration[i][j];
	boundMachine[j] += duration[i][j];
      }
  }

  for(i=0; i<nJobs; ++i) {
    if(Verbose > 1) 
      cout << "c job" << i << " " << boundJob[i] << endl;
    if( mkp < boundJob[i] ) mkp = boundJob[i];
  }
  for(i=0; i<nMachines; ++i) {
    if(Verbose > 1) 
      cout << "c machine" << i << " " << boundMachine[i] << endl;
    if( mkp < boundMachine[i] ) mkp = boundMachine[i];
  }

  return mkp;
}

int jsp_upperbound() {
  int i, j, mkp=0, task;


  if(!time_lag[0] && !time_lag[1]) {
    int boundJob[nJobs], boundMachine[nJobs]; 
    std::fill(boundJob, boundJob+nJobs, 0);
    std::fill(boundMachine, boundMachine+nJobs, 0);
    
    
    for(i=0; i<nMachines; ++i)
      {
	for(j=0; j<nJobs; ++j)
	  {
	    task = std::max( boundJob[j], boundMachine[machine[j][i]] );
	    task += duration[j][i];
	    boundJob[j] = task;
	    boundMachine[machine[j][i]] = task;
	    if( task > mkp ) mkp = task;
	  }
      }
  } else {

    for(i=0; i<nMachines; ++i)
      for(j=0; j<nJobs; ++j)
	mkp += duration[j][i];
    
  }

  if(setup_time)
    for(i=0; i<nMachines; ++i)
      mkp += (nJobs * max_setup[i]);


  //return ((setup_time || time_lag[0]) ? mkp + 1000*Multiplier : mkp);
  return mkp;
}


int jsp_lowerbound()
{
  int i, j, mkp=0, boundJob[nJobs], boundMachine[nMachines]; 
  std::fill(boundJob, boundJob+nJobs, 0);
  std::fill(boundMachine, boundMachine+nMachines, 0);

  for(i=0; i<nJobs; ++i) {
    for(j=0; j<nMachines; ++j)
      {
	boundJob[i] += duration[i][j];
	boundMachine[machine[i][j]] += duration[i][j];
      }
  }

  for(i=0; i<nJobs; ++i)
    if( mkp < boundJob[i] ) mkp = boundJob[i];
  for(i=0; i<nMachines; ++i)
    if( mkp < boundMachine[i] ) mkp = boundMachine[i];

  return mkp;
}


int jsp_setup(CSP& model, VarArray& disjuncts, VarArray& tasks) {
  int i,j,k, lb, ub, per_machine[nJobs*nMachines];

  // tasks' start-time
  for(i=0; i<nJobs; ++i) 
    for(j=0; j<nMachines; ++j) {
      //tasks.add( Variable((family_matrix ? family_matrix[0][family[i][j]] : 0), 
      //makespan-duration[i][j]) );
      lb = (release_date ? release_date[i][j] : 0);
      ub = makespan;
      if(due_date && due_date[i][j] < ub) ub = due_date[i][j];
      ub -= duration[i][j];

      if(lb <= ub)
	tasks.add( Variable(lb, ub) );
      else
	return UNSAT;
    }
 
  // order within a job
  int lag_min = -1;
  int lag_max = -1;

  for(i=0; i<nJobs; ++i)
    for(j=1; j<nMachines; ++j) {

      lag_min = duration[i][j-1];
      if(time_lag[0])
	lag_min += time_lag[0][i][j-1];

      if(time_lag[1])
	lag_max = (duration[i][j-1] + time_lag[1][i][j-1]);
      
      model.add( Precedence(tasks[i*nMachines+j-1], lag_min, tasks[i*nMachines+j]) );

      if( lag_max != -1 )
	model.add( Precedence(tasks[i*nMachines+j], -lag_max, tasks[i*nMachines+j-1]) );
      
      //model.add( Precedence(tasks[i*nMachines+j-1], duration[i][j-1], tasks[i*nMachines+j]) );
      //cout << duration[i][j-1] << endl;
    }
  // non-overlapping reified constraints
  for(i=0; i<nJobs; ++i)
    for(j=0; j<nMachines; ++j) 
      per_machine[i*nMachines+machine[i][j]] = (i*nMachines+j);
  for(k=0; k<nMachines; ++k) 
    for(i=0; i<nJobs; ++i)
      for(j=i+1; j<nJobs; ++j) {
	int x = per_machine[i*nMachines+k]/nMachines;
	int y = per_machine[i*nMachines+k]%nMachines;
	int u = per_machine[j*nMachines+k]/nMachines;
	int v = per_machine[j*nMachines+k]%nMachines;

	int release_xy=0;
	for(int z=0; z<y; ++z)
	  release_xy += duration[x][z];
	
	int release_uv=0;
	for(int z=0; z<v; ++z)
	  release_uv += duration[u][z];

	//cout << release_xy << " <? " << release_uv << endl;

	//if (ValueO != "jobrank" || y < v){
	if(!JobRank || 
	   ((JobRank == 1) && y < v) || 
	   ((JobRank == 2) && release_xy < release_uv)) {
          disjuncts.add( Disjunctive( tasks[per_machine[i*nMachines+k]],
                                      duration[x][y]+(setup_time ? setup_time[k][i][j] : 0),
                                      tasks[per_machine[j*nMachines+k]],
                                      duration[u][v]+(setup_time ? setup_time[k][j][i] : 0) ) );
        } // else if (y == v){
	//           if (randint(2)){
	//             disjuncts.add( Disjunctive( tasks[per_machine[i*nMachines+k]],
	//                                         duration[x][y]+(setup_time ? setup_time[k][i][j] : 0),
	//                                         tasks[per_machine[j*nMachines+k]],
	//                                         duration[u][v]+(setup_time ? setup_time[k][j][i] : 0) ) );
	//           }  else {
	//             disjuncts.add( Disjunctive( tasks[per_machine[j*nMachines+k]],
	//                                         duration[u][v]+(setup_time ? setup_time[k][j][i] : 0),
	//                                         tasks[per_machine[i*nMachines+k]],
	//                                         duration[x][y]+(setup_time ? setup_time[k][i][j] : 0) ) );
	//           }
	//         }
	else {
          disjuncts.add( Disjunctive( tasks[per_machine[j*nMachines+k]],
                                      duration[u][v]+(setup_time ? setup_time[k][j][i] : 0),
                                      tasks[per_machine[i*nMachines+k]],
                                      duration[x][y]+(setup_time ? setup_time[k][i][j] : 0) ) );
        }
	
	// 	//assert( duration[x][y] == jsp_duration[per_machine[i*nMachines+k]] );
	// 	disjuncts.add( Disjunctive( tasks[per_machine[i*nMachines+k]], 
	// 				    duration[x][y],
	// 				    tasks[per_machine[j*nMachines+k]], 
	// 				    duration[u][v] ) );
      }

  model.add( disjuncts );
  return UNKNOWN;
  //exit(0);
}


int permfsp_setup(CSP& model, VarArray& disjuncts, VarArray& tasks) {
  int i,j,k, lb, ub, per_machine[nJobs*nMachines];
	
  // tasks' start-time
  for(i=0; i<nJobs; ++i) 
    for(j=0; j<nMachines; ++j) {
      //tasks.add( Variable((family_matrix ? family_matrix[0][family[i][j]] : 0), 
      //makespan-duration[i][j]) );
      lb = (release_date ? release_date[i][j] : 0);
      ub = makespan;
      if(due_date && due_date[i][j] < ub) ub = due_date[i][j];
      ub -= duration[i][j];
			
      if(lb <= ub)
	tasks.add( Variable(lb, ub) );
      else
	return UNSAT;
    }
	
  // order within a job
  int lag_min = -1;
  int lag_max = -1;
	
  for(i=0; i<nJobs; ++i)
    for(j=1; j<nMachines; ++j) {
			
      lag_min = duration[i][j-1];
      if(time_lag[0])
	lag_min += time_lag[0][i][j-1];
			
      if(time_lag[1])
	lag_max = (duration[i][j-1] + time_lag[1][i][j-1]);
			
      model.add( Precedence(tasks[i*nMachines+j-1], lag_min, tasks[i*nMachines+j]) );
			
      if( lag_max != -1 )
	model.add( Precedence(tasks[i*nMachines+j], -lag_max, tasks[i*nMachines+j-1]) );
			
      //model.add( Precedence(tasks[i*nMachines+j-1], duration[i][j-1], tasks[i*nMachines+j]) );
      //cout << duration[i][j-1] << endl;
    }
  // non-overlapping reified constraints
  for(i=0; i<nJobs; ++i)
    for(j=0; j<nMachines; ++j) 
      per_machine[i*nMachines+machine[i][j]] = (i*nMachines+j);

  
  //EH We need only one Boolean var for all the disjuncts of the ith tasl of a job:

  for(i=0; i<ndisjuncts; ++i) {
    Variable B(0,1);
    disjuncts.add(B);
  }

  model.add(disjuncts);

  int dsj;
  for(k=0; k<nMachines; ++k) {
    dsj = 0;
    for(i=0; i<nJobs; ++i)
      for(j=i+1; j<nJobs; ++j) {
	int x = per_machine[i*nMachines+k]/nMachines;
	int y = per_machine[i*nMachines+k]%nMachines;
	int u = per_machine[j*nMachines+k]/nMachines;
	int v = per_machine[j*nMachines+k]%nMachines;
				
	int release_xy=0;
	for(int z=0; z<y; ++z)
	  release_xy += duration[x][z];
				
	int release_uv=0;
	for(int z=0; z<v; ++z)
	  release_uv += duration[u][z];

	model.add( disjuncts[dsj++] ==
		   Disjunctive( tasks[per_machine[i*nMachines+k]],
				duration[x][y]+(setup_time ? setup_time[k][i][j] : 0),
				tasks[per_machine[j*nMachines+k]],
				duration[u][v]+(setup_time ? setup_time[k][j][i] : 0) )
		   );
      }
  }

//   // DG Create a chain of disjuncts where each set of disjuncts i from (0 to nMachines-1) modulo nMachines is set to the next disjunct
//   for(i=1; i<ndisjuncts; ++i) {
//     j = i-1;
//     if (i%nMachines != 0) // 
//       disjuncts[i]==disjuncts[j]; //Probably wrong way to state var disjuncts.
//   }	
  
				
	
  return UNKNOWN;
}
			  

int permfspjb_setup(CSP& model, VarArray& disjuncts, VarArray& tasks, VarArray& jobdisjuncts) {
  int i,j,k, lb, ub, per_machine[nJobs*nMachines];
	
  // tasks' start-time
  for(i=0; i<nJobs; ++i) 
    for(j=0; j<nMachines; ++j) {
      //tasks.add( Variable((family_matrix ? family_matrix[0][family[i][j]] : 0), 
      //makespan-duration[i][j]) );
      lb = (release_date ? release_date[i][j] : 0);
      ub = makespan;
      if(due_date && due_date[i][j] < ub) ub = due_date[i][j];
      ub -= duration[i][j];
			
      if(lb <= ub)
	tasks.add( Variable(lb, ub) );
      else
	return UNSAT;
    }
	
  // order within a job
  int lag_min = -1;
  int lag_max = -1;
	
  for(i=0; i<nJobs; ++i)
    for(j=1; j<nMachines; ++j) {
			
      lag_min = duration[i][j-1];
      if(time_lag[0])
	lag_min += time_lag[0][i][j-1];
			
      if(time_lag[1])
	lag_max = (duration[i][j-1] + time_lag[1][i][j-1]);
			
      model.add( Precedence(tasks[i*nMachines+j-1], lag_min, tasks[i*nMachines+j]) );
			
      if( lag_max != -1 )
	model.add( Precedence(tasks[i*nMachines+j], -lag_max, tasks[i*nMachines+j-1]) );
			
      //model.add( Precedence(tasks[i*nMachines+j-1], duration[i][j-1], tasks[i*nMachines+j]) );
      //cout << duration[i][j-1] << endl;
    }
  // non-overlapping reified constraints
  for(i=0; i<nJobs; ++i)
    for(j=0; j<nMachines; ++j) 
      per_machine[i*nMachines+machine[i][j]] = (i*nMachines+j);
  for(k=0; k<nMachines; ++k) 
    for(i=0; i<nJobs; ++i)
      for(j=i+1; j<nJobs; ++j) {
	int x = per_machine[i*nMachines+k]/nMachines;
	int y = per_machine[i*nMachines+k]%nMachines;
	int u = per_machine[j*nMachines+k]/nMachines;
	int v = per_machine[j*nMachines+k]%nMachines;
				
	int release_xy=0;
	for(int z=0; z<y; ++z)
	  release_xy += duration[x][z];
				
	int release_uv=0;
	for(int z=0; z<v; ++z)
	  release_uv += duration[u][z];
				
	disjuncts.add( Disjunctive( tasks[per_machine[i*nMachines+k]],
				    duration[x][y]+(setup_time ? setup_time[k][i][j] : 0),
				    tasks[per_machine[j*nMachines+k]],
				    duration[u][v]+(setup_time ? setup_time[k][j][i] : 0) ) );				
      }
	
  model.add( disjuncts );
	
  // DG Create a chain of disjuncts where each set of disjuncts i from (0 to nMachines-1) modulo nMachines is set to the next disjunct
  for(i=1; i<ndisjuncts; ++i) {
    j = i-1;
    if (i%nMachines != 0) // 
      disjuncts[i]==disjuncts[j]; //Probably wrong way to state var disjuncts.
  }		

  // DG Create Job disjuncts (basically repeats the disjuncts on the first machine between tasks)
  // i.e. says that if disjunct true then task 1 of job i comes before task 1 of job j and vice versa if false.
  for(i=0; i<nJobs; ++i)
    for(j=0; j<nMachines; ++j){ 
      int x = per_machine[i*nMachines+1]/nMachines;
      int y = per_machine[i*nMachines+1]%nMachines;
      int u = per_machine[j*nMachines+1]/nMachines;
      int v = per_machine[j*nMachines+1]%nMachines;
			
      int release_xy=0;
      for(int z=0; z<y; ++z)
	release_xy += duration[x][z];
			
      int release_uv=0;
      for(int z=0; z<v; ++z)
	release_uv += duration[u][z];
			
      jobdisjuncts.add( Disjunctive( tasks[per_machine[i*nMachines+1]],
				     duration[x][y]+(setup_time ? setup_time[1][i][j] : 0),
				     tasks[per_machine[j*nMachines+1]],
				     duration[u][v]+(setup_time ? setup_time[1][j][i] : 0) ) );	
    }
	
  model.add( jobdisjuncts ); // This adds job disjuncts, but how do we force it to just branch on these?
	
  return UNKNOWN;
}

		
int jtl0_setup(CSP& model, VarArray& disjuncts, VarArray& tasks) {
  int i,j,k, lb, ub, per_machine[nJobs*nMachines];

  // tasks' start-time
  for(i=0; i<nJobs; ++i) {
    lb = 0;
    ub = makespan;
    for(j=0; j<nMachines; ++j) 
      ub -= duration[i][j];

    if(lb <= ub) {
      //std::cout << "job_" << i << " [" << lb << ".." << ub << "]" << std::endl;
      tasks.add( Variable(lb, ub) );
    } else
      return UNSAT;
  }
 
  // non-overlapping reified constraints
  for(i=0; i<nJobs; ++i)
    for(j=0; j<nMachines; ++j) 
      per_machine[i*nMachines+machine[i][j]] = (i*nMachines+j);
  for(k=0; k<nMachines; ++k) 
    for(i=0; i<nJobs; ++i)
      for(j=i+1; j<nJobs; ++j) {
	int x = per_machine[i*nMachines+k]/nMachines;
	int y = per_machine[i*nMachines+k]%nMachines;

	int u = per_machine[j*nMachines+k]/nMachines;
	int v = per_machine[j*nMachines+k]%nMachines;

	int release_xy=0;
	for(int z=0; z<y; ++z)
	  release_xy += duration[x][z];
	
	int release_uv=0;
	for(int z=0; z<v; ++z)
	  release_uv += duration[u][z];
	
	//task[i]+release_xy+duration[x][y]-release_uv <= task[j] \/ task[j]+release_uv+duration[u][v]-release_xy <= task[i]


// 	std::cout << "job_" << i << " + " << (release_xy + duration[x][y] - release_uv) 
// 		  << " <= job_" << j << " \\/ job_" 
		  
// 		  << j << " + " << (release_uv + duration[u][v] - release_xy)
// 		  << " <= job_" << i << std::endl;


	if(!JobRank || 
	   ((JobRank == 1) && y < v) || 
	   ((JobRank == 2) && release_xy < release_uv)) {
          disjuncts.add( Disjunctive( tasks[i],
                                      duration[x][y]+release_xy-release_uv,
                                      tasks[j],
                                      duration[u][v]+release_uv-release_xy ) );
        } else {
          disjuncts.add( Disjunctive( tasks[j],
                                      duration[u][v]+release_uv-release_xy,
				      tasks[i],
                                      duration[x][y]+release_xy-release_uv ) );
        }
      }

  model.add( disjuncts );
  return UNKNOWN;
}


void osp_setup(CSP& model, VarArray& disjuncts, VarArray& tasks) {
  int i,j,k;
  for(i=0; i<nJobs; ++i) 
    for(j=0; j<nMachines; ++j) 
      tasks.add( Variable(0, makespan-duration[i][j]) );

  for(i=0; i<nMachines; ++i) {
    for(j=0; j<nJobs; ++j) 
      for(k=j+1; k<nJobs; ++k) 
	disjuncts.add( Disjunctive( tasks[j*nMachines+i], duration[j][i],
				    tasks[k*nMachines+i], duration[k][i] ) );
  }

  for(i=0; i<nJobs; ++i) {
    for(j=0; j<nMachines; ++j) 
      for(k=j+1; k<nMachines; ++k) 
	disjuncts.add( Disjunctive( tasks[i*nMachines+j], duration[i][j],
				    tasks[i*nMachines+k], duration[i][k] ) );
  }

  model.add( disjuncts );

}


int get_jsp_makespan(VarArray& tasks)
{
  int d, mkp = 0;
  for(int i=0; i<nJobs; ++i) {
    d = (tasks[i*nMachines+nMachines-1].value() + duration[i][nMachines-1]);
    if(d > mkp) mkp = d;
  }
  return mkp;
}

int get_jtl0_makespan(VarArray& tasks)
{
  int d, mkp = 0;
  for(int i=0; i<nJobs; ++i) {
    d = tasks[i].value();
    for(int j=0; j<nMachines; ++j)
      d += duration[i][j];
    if(d > mkp) mkp = d;
  }
  return mkp;
}

int get_osp_makespan(VarArray& tasks)
{
  int d, mkp = 0;
  for(int i=0; i<nJobs; ++i) 
    for(int j=0; j<nMachines; ++j) {
      d = (tasks[i*nMachines+j].value() + duration[i][j]);
      if(d > mkp) mkp = d;
    }
  return mkp;
}

int get_makespan(VarArray& tasks)
{
  int mkp=0;
  if(Type == "osp")
    mkp = get_osp_makespan(tasks);
  else if(Type == "jtl0")
    mkp = get_jtl0_makespan(tasks);
  else
    mkp = get_jsp_makespan(tasks);
  return mkp;
}


class ResetBase : public SolutionMethod {

public:
  ResetBase(Solver *s) : SolutionMethod(s) {}
  virtual ~ResetBase() {}
  
  virtual void execute() 
  { 
    solver->fail_increment = Base;

    solver->FAILLIMIT = solver->FAILURES+Base;

    //cout << endl << "Solution Found!!!" << endl; 
  }

  virtual void initialise() 
  {
  }
};


void setWeights(Solver& s, VarArray& tasks) {
  if(disjunct_weights.size() > 0) {

    int i, n=s.length, m, nw, min_weight=NOVAL, max_weight=0;
    double avg_weight = 0;
    MistralNode<Constraint*>* nd;
    for(i=0; i<n; ++i) {
      //s.variables[i]->weight = disjunct_weights.back()[disjunct_index[i]];
      if(Weight > 1)
	nw = normalised_weights.back()[disjunct_index[i]];
      else
	nw = disjunct_weights.back()[disjunct_index[i]];
      nd = s.variables[i]->constraintsOnValue();
      if(nextNode(nd) && Rngd && nd->elt->arity !=3)
	nextNode(nd);
      //       nd->elt->print(std::cout);
      //       std::cout << std::endl;
      
      nd->elt->weight = nw;
      //nw = s.variables[i]->weight;
      avg_weight += (double)nw;
      if(min_weight > nw) min_weight = nw;
      if(max_weight < nw) max_weight = nw;
    }

    //     cout << endl << "c  min weight: "<< min_weight << "  total weight: " << avg_weight << "  max weight: " << max_weight ;
    //     avg_weight /= (double)n;
    //     cout << "  avg weight: " << avg_weight << endl;

    //     cout << endl << min_weight << " " << avg_weight 
    // 	 << " " << max_weight ;
    //     avg_weight /= (double)n;
    //     cout << " " << avg_weight << endl;


    double factor = (double)num_weight_updates / (double)(num_weight_updates+1);
    ++num_weight_updates;
    avg_disjunct_avg_weight = (avg_disjunct_avg_weight * factor + 
			       avg_weight / (double)num_weight_updates);
    avg_disjunct_min_weight = (avg_disjunct_min_weight * factor + 
			       (double)min_weight / (double)num_weight_updates);
    avg_disjunct_max_weight = (avg_disjunct_max_weight * factor + 
			       (double)max_weight / (double)num_weight_updates);

    if(min_disjunct_max_weight > max_weight) 
      min_disjunct_max_weight = max_weight;
    if(min_disjunct_avg_weight > avg_weight) 
      min_disjunct_avg_weight = avg_weight;
    if(min_disjunct_min_weight > min_weight) 
      min_disjunct_min_weight = min_weight;

    if(max_disjunct_max_weight < max_weight) 
      max_disjunct_max_weight = max_weight;
    if(max_disjunct_avg_weight < avg_weight) 
      max_disjunct_avg_weight = avg_weight;
    if(max_disjunct_min_weight < min_weight) 
      max_disjunct_min_weight = min_weight;


    min_weight=NOVAL;
    max_weight=0;
    avg_weight = 0;
    n = s.numvars;
    m = s.length;
    for(i=0; i<n; ++i) {
      //s.variables[i]->print( cout );
      //cout << ": ";
      
      nw = 0;
      nd = s.variables[i]->constraintsOnValue();
      while( nextNode(nd) ) {
	//nd->elt->print(cout);
	//cout << " ";
	nw += nd->elt->weight;
      }
      s.variables[i]->weight = nw;

      //cout << " -> " << nw << endl;

      if(i>=m) {
	avg_weight += (double)nw;
	if(min_weight > nw) min_weight = nw;
	if(max_weight < nw) max_weight = nw;
      }
    }

    //     cout << "c  min weight: "<< min_weight << "  total weight: " << avg_weight << "  max weight: " << max_weight ;
    //     avg_weight /= (double)(n-m);
    //     cout << "  avg weight: " << avg_weight << endl;


    //     min_weight=NOVAL;
    //     max_weight=0;
    //     avg_weight = 0;
    //     BuildObject *bvar;
    //     VariableInt *temp;
    //     for(i=0; i<m; ++i) {
    //       bvar = tasks[i].var_ptr_;
    //       if(bvar->isSearchable() && (temp = bvar->getVariable())) {
    // 	nw = task_weights.back()[i];
    // 	//temp->weight = nw; 

    // 	avg_weight += (double)nw;
    // 	if(min_weight > nw) min_weight = nw;
    // 	if(max_weight < nw) max_weight = nw;
    //       }
    //     }

    //     cout << min_weight << " " << avg_weight << " " << max_weight ;
 
    //     avg_weight /= (double)n;

    //     cout << " " << avg_weight << endl;

    //     factor = (double)(num_weight_updates-1) / (double)num_weight_updates;
    
    avg_task_avg_weight = (avg_task_avg_weight * factor + 
			   avg_weight / (double)num_weight_updates);
    avg_task_min_weight = (avg_task_min_weight * factor + 
 			   min_weight / (double)num_weight_updates);
    avg_task_max_weight = (avg_task_max_weight * factor + 
			   max_weight / (double)num_weight_updates);
    if(min_task_max_weight > max_weight) 
      min_task_max_weight = max_weight;
    if(min_task_avg_weight > avg_weight) 
      min_task_avg_weight = avg_weight;
    if(min_task_min_weight > min_weight) 
      min_task_min_weight = min_weight;
    if(max_task_max_weight < max_weight) 
      max_task_max_weight = max_weight;
    if(max_task_avg_weight < avg_weight) 
      max_task_avg_weight = avg_weight;
    if(max_task_min_weight < min_weight) 
      max_task_min_weight = min_weight;

  }
}


void updateBestSolution(VarArray& disjuncts, VarArray& tasks) {
  int i,
    *new_solution = new int[ndisjuncts],
    n_discrepancies=0;
  //*new_dweight  = new int[ndisjuncts],
  //*new_tweight  = new int[ntasks];
  //BuildObject *bvar;
  //VariableInt *temp;
  for(i=0; i<ndisjuncts; ++i) {
    //     bvar = disjuncts[i].var_ptr_;
    //     if(bvar->isSearchable() && (temp = bvar->getVariable())) {
    //       disjunct_index[j++] = i;
    //       new_dweight[i] = temp->weight;
    //     } else new_dweight[i] = 1;
    new_solution[i] = disjuncts[i].value();
    if(new_solution[i]) {
      zeros[i] /= 2;
      if(zeros[i] < 1) zeros[i] = 1;
      ones[i] = Tprob-zeros[i];
      probability[i] = ones[i];
      n_discrepancies += (1-best_solution[i]);
    } else {
      ones[i] /= 2;
      if(ones[i] < 1) ones[i] = 1;
      zeros[i] = Tprob-ones[i];
      probability[i] = zeros[i];
      n_discrepancies += best_solution[i];
    }
  }

  cout << left << setw(30) << "c Distance " << ":" << right << setw(20) << n_discrepancies << endl;

  //   for(i=0; i<ntasks; ++i) {
  //     bvar = tasks[i].var_ptr_;
  //     if(bvar->isSearchable() && (temp = bvar->getVariable())) {
  //       new_tweight[i] = temp->weight;
  //     } else new_tweight[i] = temp->degree;
  //   }

  //   task_weights.push_back(new_tweight);
  //   disjunct_weights.push_back(new_dweight);

  

  solutions.push_back(new_solution);
  best_solution = new_solution;

  if(solutions.size() > 1) {
    double factor = ((double)total_solutions / (double)(total_solutions+1));
    ++total_solutions;
    //    cout << "avg distance " << factor << " * " << avg_distance << " + " 
    // 	 << ((double)n_discrepancies) << " / " << ((double)total_solutions) << endl;

    avg_distance = factor * avg_distance + 
      (double)n_discrepancies/(double)total_solutions;
    if(min_distance > n_discrepancies) min_distance = n_discrepancies;
    if(max_distance < n_discrepancies) max_distance = n_discrepancies;
  }
}



void updateWeights(VarArray& disjuncts, VarArray& tasks) {
  int i, j=0, //, ntasks=tasks.size(), 
    *new_dweight  = new int[ndisjuncts],
    *new_norm_weight  = new int[ndisjuncts],
    min_weight=NOVAL, max_weight=0;
  //*new_tweight  = new int[ntasks];
  BuildObject *bvar;
  VariableInt *temp;
  MistralNode<Constraint*>* nd;
  bool is_weighted;

  

  for(i=0; i<ndisjuncts; ++i) {
    bvar = disjuncts[i].var_ptr_;
    if(bvar->isSearchable() && (temp = bvar->getVariable())) {
      
      is_weighted= false;
      nd = temp->constraintsOnValue();
      if(nd) while( nextNode(nd) )
	       if(nd->elt->arity == 3) {
		 is_weighted = true;
		 break;	    
	       }
     
      if(is_weighted) {
	new_dweight[i] = nd->elt->weight;
      } else
	new_dweight[i] = 0;

      if(new_dweight[i] > max_weight) max_weight = new_dweight[i];
      if(new_dweight[i] < min_weight) min_weight = new_dweight[i];

      disjunct_index[j++] = i;
    } else {
      new_dweight[i] = 1;
    }
  }

  //   for(i=0; i<ntasks; ++i) {
  //     bvar = tasks[i].var_ptr_;
  //     if(bvar->isSearchable() && (temp = bvar->getVariable())) {
  //       new_tweight[i] = temp->weight;
  //     } else new_tweight[i] = temp->degree;
  //   }

  //   task_weights.push_back(new_tweight);
  disjunct_weights.push_back(new_dweight);

  if(Weight > 1) {
    max_weights.push_back(max_weight);
    min_weights.push_back(min_weight);
    j = (max_weight - min_weight);
    if(j<1) j=1;
    for(i=0; i<ndisjuncts; ++i) {
      new_norm_weight[i] = 1 + ((Weight*(new_dweight[i]-min_weight)) / j);
    }
    normalised_weights.push_back(new_norm_weight);
  }
}



class SolutionGuidedSearch : public SolutionMethod {

protected:
  int *_index;
  vector<int*> _solutions;
  vector<int*> _pool;

  int guide_prob; // over 1000
  int max_pool_size;

  ValSelector **not_used;  
  int guide;
  int same_search;

public:

  SolutionGuidedSearch(Solver *s, VarArray& X, int ps=4, int gp=750) : SolutionMethod(s) 
  {
    int i, j=0, n=solver->length;

    max_pool_size = ps;
    guide_prob = gp;
    same_search = false;

    _index = new int[n];
 
    BuildObject *bvar;
    for(i=0; i<ndisjuncts; ++i) {
      bvar = X[i].var_ptr_;
      if( bvar->isSearchable() && (bvar->getVariable() != NULL) ) 
	_index[j++] = i;
    }

    for(i=0; (unsigned int)i<solutions.size(); ++i) {
      int* copy = new int[n];
      for(j=0; j<n; ++j)
	copy[j] = solutions[i][_index[j]];
      _solutions.push_back(copy);
    }

    for(i=(int)(_solutions.size())-1; i>=0 && i>=(int)(_solutions.size())-max_pool_size; --i) {
      _pool.push_back(_solutions[i]);
    }

    //assert(_pool.size() > 0);
    assert(_pool.size() <= (unsigned int)max_pool_size);
    guide = _pool.size()-1;
    solver->setGuidedOrdering(X, best_solution);

    if(guide_prob < 1000) {
      not_used = new ValSelector*[n];
      for(i=0; i<n; ++i)
	not_used[i] = new ValSelectorRand(solver->variables[i]);
    } else not_used = NULL;
  }

  virtual ~SolutionGuidedSearch() 
  {
    delete [] _index;     
    for(unsigned int i=0; i<_solutions.size(); ++i)
      delete [] _solutions[i];
    if(not_used) {
      for(int i=0; i<solver->length; ++i)
	delete not_used[i];
    }
    delete [] not_used;
  }
  
  virtual void execute() 
  { 
    if(Verbose>1)
      cout << "c found an improving solution" << endl;

    unsigned int i, n = solver->length;
    int *new_solution = new int[n];
    
    ValSelectorGuided *vo;
    for(i=0; i<n; ++i) {
      //vo = (ValSelectorGuided *)(solver->variables[i]->branch);
      /* vo->ideal = */ 

      new_solution[i] = solver->solution[i];
      if(guide_prob >= 1000) {
	vo = (ValSelectorGuided *)(solver->variables[i]->branch);
	vo->ideal = new_solution[i];
      }
    }

    _solutions.push_back(new_solution);

    if(_pool.size() < (unsigned int)max_pool_size) { // not enough solution -> fill the pool!
      if(Verbose>1)
	cout << "c add new solution to the pool" << endl;
      _pool.push_back(new_solution);
    } else {
      if(!same_search) { // within the same search, we replace the same solution
	int start = std::max(guide, 0);
	if(Verbose>1)
	  cout << "c replace solution " << start << endl;
	for(i=start; i<(unsigned int)max_pool_size-1; ++i)
	  _pool[i] = _pool[i+1];
      }
      _pool[max_pool_size-1] = new_solution;
    }

    if( Reset ) {
      
      if(Verbose>1)
	cout << "c reset the restarts" << endl;

      solver->fail_increment = Base;
      solver->FAILLIMIT = solver->FAILURES+Base;
    }

    if(Verbose>1) {
      cout << "c ";
      for(i=1; i<_pool.size(); ++i) 
	for(unsigned int j=0; j<i; ++j) {
	  int disc = 0;
	  for(int k=0; k<solver->length; ++k)
	    if(_pool[i][k] != _pool[j][k]) ++disc;
	  cout << " " << disc;
	}
      cout << endl;
    }

    same_search = true;
  }


  virtual void initialise() 
  {
    same_search = false;

    if(guide_prob < 1000) {

      int i, n=solver->length;
      ValSelector *aux;

      if(Verbose>1)
	cout << "\nc initialise solution guided search method, restarts for " 
	     << (solver->FAILLIMIT - solver->FAILURES) << " fails"  << endl;

      int p=randint(1000);
      
      if(Verbose>1)
	cout << "c (" << ((double)p/1000.0) << " / " << ((double)guide_prob/1000.0) << ") ";
      
      if(_pool.size() > 0 && p < guide_prob) {

	if(Verbose>1)
	  cout << "guide with solution ";

	// solution guided
	if(guide < 0)
	  for(i=0; i<n; ++i) {
	    aux = solver->variables[i]->branch;
	    solver->variables[i]->branch = not_used[i];
	    not_used[i] = aux;
	  }

	guide = randint((_pool.size()));
	for(i=0; i<n; ++i) 
	  ((ValSelectorGuided*)(solver->variables[i]->branch))->ideal = _pool[guide][i];

	if(Verbose>1)
	  cout << guide << " over " << (_pool.size()) << endl;

      } else {
	// not guided
	if(guide >= 0)
	  for(i=0; i<n; ++i) {
	    aux = solver->variables[i]->branch;
	    solver->variables[i]->branch = not_used[i];
	    not_used[i] = aux;
	  }
	
	guide = -1;

	if(Verbose>1)
	  cout << "start random search" << endl;
	
      }
    }
  }
};




class SolutionRandGuidedSearch : public SolutionMethod {

protected:
  int *_index;
  vector<int*> _solutions;
  int *_zeros;
  int *_ones;

 
public:
  SolutionRandGuidedSearch(Solver *s, VarArray& X) : SolutionMethod(s) 
  {
    unsigned int i, j=0, n = solver->length;
    _index = new int[n];
 
    BuildObject *bvar;
    for(i=0; i<(unsigned int)ndisjuncts; ++i) {
      bvar = X[i].var_ptr_;
      if( bvar->isSearchable() && (bvar->getVariable() != NULL) ) 
	_index[j++] = i;
    }

    for(i=0; i<solutions.size(); ++i) {
      int* copy = new int[n];
      for(j=0; j<n; ++j)
	copy[j] = solutions[i][_index[j]];
      _solutions.push_back(copy);
    }

    _zeros = new int[n];
    _ones = new int[n];
    for(i=0; i<n; ++i) {
      _zeros[i] = zeros[_index[i]];
      _ones[i] = ones[_index[i]];
    }

    solver->setRandGuidedOrdering(X, best_solution, probability, range);
  }

  virtual ~SolutionRandGuidedSearch() 
  {
    delete [] _index;     
    for(unsigned int i=0; i<_solutions.size(); ++i)
      delete [] _solutions[i];
    
    delete [] _zeros;
    delete [] _ones;
  }
  
  virtual void execute() 
  { 

    if(Verbose>1)
      cout << "c  found an improving solution, update probabilities:" << endl;

    int n = solver->length;
    int i, val, *new_solution = new int[n];

    ValSelectorRandGuided *vo;
    for(i=0; i<n; ++i) {
      vo = (ValSelectorRandGuided *)(solver->variables[i]->branch);
      val = solver->solution[i];

      vo->ideal = val;
      
      if(val) {
	_zeros[i] /= 2;
	if(_zeros[i] < 1) _zeros[i] = 1;
	_ones[i] = Tprob-_zeros[i];
	vo->proba = _ones[i];
      } else {
	_ones[i] /= 2;
	if(_ones[i] < 1) _ones[i] = 1;
	_zeros[i] = Tprob-_ones[i];
	vo->proba = _zeros[i];
      }
      
      new_solution[i] = val;
    }

    if( Reset ) {
      solver->fail_increment = Base;
      solver->FAILLIMIT = solver->FAILURES+Base;
    }
  }
  

  virtual void initialise() 
  {
  }
  
};


void dichotomic_search()
{
  int new_makespan = 0;
  double dsopttime = 0.0,  maxidstime=0.0;
  bool first_step = true;

  ////////// dichotomic search ///////////////
  while( Gap < (maxfsble-minfsble) && Dichotomy-- && minfsble < maxfsble ) {
    
    if(first_step)
      makespan = (int)(floor((Skew * (double)minfsble + (double)maxfsble)/(1.0+Skew)));
    else
      makespan = (int)(floor(((double)minfsble + (double)maxfsble)/2));
      

    cout << "c ============[ start dichotomic step ]============" << endl;
    cout << left << setw(30) << "c current dichotomic range" << ":" 
	 << right << setw(6) << " " << setw(5) << minfsble << " to " << setw(5) << maxfsble << std::endl;

    cout << left << setw(30) << "c target makespan" << ":"  << right << setw(20) << makespan << std::endl;
    cout << "c ===========[ start dichotomic search ]===========" << endl;
      
    VarArray disjuncts;
    VarArray tasks;
    CSP model;
    Result = UNKNOWN;

    if(Type == "osp")
      osp_setup(model, disjuncts, tasks);
    else if(Type == "jtl0")
      Result = jtl0_setup(model, disjuncts, tasks);     
<<<<<<< .mine
    else if(Type == "fsp")
      Result = permfsp_setup(model, disjuncts, tasks);
    else
=======
    else if(Type == "pfsp")
      Result = permfsp_setup(model, disjuncts, tasks);     
    else
>>>>>>> .r396
      Result = jsp_setup(model, disjuncts, tasks);     


    std::cout << UNKNOWN << " " << Result << " " << (disjuncts.size()) << " " << ndisjuncts << std::endl;

    assert(disjuncts.size() == ndisjuncts);    


    if(Result == UNKNOWN) {
      Solver s(model, disjuncts);
 
//       s.print(cout);
      
//       exit(1);


      if(s.status == UNKNOWN) {
	
	if( Rngd ) {
	  s.setRestartNogood();
	  s.setForgetfulness( 0.0 );
	}
	if( Randomized > 0 ) s.setRandomized();
	s.setTimeLimit( Cutoff );
	if(first_step) {
	  addHeuristic( s, Heuristic, Randomized, I_VO );
	} else
	  addHeuristic( s, Heuristic, Randomized, D_VO );
	s.setVerbosity( Verbose-1 );
	s.setRandomSeed( Seed );
	
	Result = UNKNOWN;
	if( SAC ) Result = s.sacPreprocess( SAC > 1 );
      }
      if( Result != UNSAT ) {
	if(first_step) {
	  if(I_VO == RAND)
	    s.setRandomValueOrdering();
	} else
	  switch(D_VO) {
	  case RAND: {
	    s.setRandomValueOrdering();
	  } break;
	  case GUIDED: {
	    s.setGuidedOrdering(disjuncts, best_solution);
	  } break;
	  case RGUIDED: {
	    s.setRandGuidedOrdering(disjuncts, best_solution, probability, range);
	  } break;
	  }
       
	if(Weight>0)
	  setWeights(s,tasks);

	Result = s.solve_and_restart(PolicyRestart, Base, Factor, Decay);
      }

      if( Result == SAT ) {

	updateBestSolution(disjuncts, tasks);

	new_makespan = (get_makespan(tasks));
	cout << left << setw(30) << "c solutions's makespan" << ":" << right << setw(20) << new_makespan << endl;
	dsopttime = (getRunTime() - total_time);
      } // else {
      // 	cout << left << setw(30) << "c " << ":" << right << setw(20) << "no solution";
      //       }

      if((Result != UNKNOWN) && (s.ENDTIME > maxidstime))
	maxidstime = s.ENDTIME;
      
      if(Weight>0)
	updateWeights(disjuncts, tasks);


      s.printStatistics(std::cout, (RUNTIME + BTS + PPGS + OUTCOME) );
      
      total_nodes += s.NODES;
      total_fails += s.FAILURES;
      total_bts += s.BACKTRACKS;
      total_propags += s.PROPAGS;
      total_restarts += s.BTSLIST.size;
      
      cout << "c =============[ end dichotomic step ]=============" << endl;
      
    } 

    cout << endl;

    if( Result == SAT ) {
      maxfsble = new_makespan;
	  first_step = false;
    } else {
      minfsble = makespan+1;
      if( Result != SAT && Result != UNSAT )
	Solved = false;
      else
	max_infeasible = makespan;
    }
    
  }

  Solved = (maxfsble == max_infeasible+1);

  opt_time = dsopttime;
  dsrun_time = (getRunTime() - total_time);
  proof_time = (dsrun_time - opt_time);

  cout << "c =================[ statistics ]==================" << endl;
  cout << left << setw(30) << "d DSLOWERBOUND "    << right << setw(21) << (max_infeasible+1) << endl
       << left << setw(30) << "d DSOBJECTIVE "     << right << setw(21) << maxfsble << endl 
       << left << setw(30) << "d DSRUNTIME "       << right << setw(21) << dsrun_time << endl
       << left << setw(30) << "d DSOPTTIME "       << right << setw(21) << dsopttime << endl
       << left << setw(30) << "d DSIDSTIME "       << right << setw(21) << maxidstime << endl
       << left << setw(30) << "d DSNODES "         << right << setw(21) << total_nodes << endl
       << left << setw(30) << "d DSBACKTRACKS "    << right << setw(21) << total_bts << endl
       << left << setw(30) << "d DSFAILS "         << right << setw(21) << total_fails << endl
       << left << setw(30) << "d DSPROPAGS "       << right << setw(21) << total_propags << endl
       << left << setw(30) << "d DSRESTARTS "      << right << setw(21) << total_restarts << endl

       << left << setw(30) << "d DSSOLUTIONS "     << right << setw(21) << total_solutions << endl
       << left << setw(30) << "d AVGDISTANCE "     << right << setw(21) << avg_distance << endl
       << left << setw(30) << "d MINDISTANCE "     << right << setw(21) << min_distance << endl
       << left << setw(30) << "d MAXDISTANCE "     << right << setw(21) << max_distance << endl;

  if(Weight) 
    cout << left << setw(30) << "d AVGAVGDWEIGHT " << right << setw(21) << avg_disjunct_avg_weight << endl
	 << left << setw(30) << "d AVGMINDWEIGHT " << right << setw(21) << avg_disjunct_min_weight << endl
	 << left << setw(30) << "d AVGMAXDWEIGHT " << right << setw(21) << avg_disjunct_max_weight << endl
	 << left << setw(30) << "d MINAVGDWEIGHT " << right << setw(21) << min_disjunct_avg_weight << endl
	 << left << setw(30) << "d MINMINDWEIGHT " << right << setw(21) << min_disjunct_min_weight << endl
	 << left << setw(30) << "d MINMAXDWEIGHT " << right << setw(21) << min_disjunct_max_weight << endl
	 << left << setw(30) << "d MAXAVGDWEIGHT " << right << setw(21) << max_disjunct_avg_weight << endl
	 << left << setw(30) << "d MAXMINDWEIGHT " << right << setw(21) << max_disjunct_min_weight << endl
	 << left << setw(30) << "d MAXMAXDWEIGHT " << right << setw(21) << max_disjunct_max_weight << endl
      
	 << left << setw(30) << "d AVGAVGTWEIGHT " << right << setw(21) << avg_task_avg_weight << endl
	 << left << setw(30) << "d AVGMINTWEIGHT " << right << setw(21) << avg_task_min_weight << endl
	 << left << setw(30) << "d AVGMAXTWEIGHT " << right << setw(21) << avg_task_max_weight << endl
	 << left << setw(30) << "d MINAVGTWEIGHT " << right << setw(21) << min_task_avg_weight << endl
	 << left << setw(30) << "d MINMINTWEIGHT " << right << setw(21) << min_task_min_weight << endl
	 << left << setw(30) << "d MINMAXTWEIGHT " << right << setw(21) << min_task_max_weight << endl
	 << left << setw(30) << "d MAXAVGTWEIGHT " << right << setw(21) << max_task_avg_weight << endl
	 << left << setw(30) << "d MAXMINTWEIGHT " << right << setw(21) << max_task_min_weight << endl
	 << left << setw(30) << "d MAXMAXTWEIGHT " << right << setw(21) << max_task_max_weight << endl;      

  cout << left << setw(30) << "d DSOPTIMAL " << right << setw(21) << Solved << endl;
  cout << "c =================[ statistics ]==================" << endl
       << endl;
  
  ////// END DIHOTOMIC SEARCH /////////////

}



void mks_dec_search()
{
  int new_makespan = 0, mks_dec = abs(Dichotomy);
  double dsopttime = 0.0,  maxidstime=0.0;
  bool first_step = true;
	
  ////////// dichotomic search ///////////////
  while( minfsble < maxfsble ) {
		
    if((maxfsble+1) > (minfsble+mks_dec))
      makespan = maxfsble-mks_dec;
    else
      makespan = (int)(floor(((double)minfsble + (double)maxfsble)/2));
		
		
    cout << "c ============[ start dichotomic step ]============" << endl;
    cout << left << setw(30) << "c current dichotomic range" << ":" 
	 << right << setw(6) << " " << setw(5) << minfsble << " to " << setw(5) << maxfsble << std::endl;
		
    cout << left << setw(30) << "c target makespan" << ":"  << right << setw(20) << makespan << std::endl;
    cout << "c ===========[ start dichotomic search ]===========" << endl;
		
    VarArray disjuncts;
    VarArray tasks;
    CSP model;
    Result = UNKNOWN;
		
    if(Type == "osp")
      osp_setup(model, disjuncts, tasks);
    else if(Type == "jtl0")
      Result = jtl0_setup(model, disjuncts, tasks);
    else //if(Type == "jsp")
      Result = jsp_setup(model, disjuncts, tasks);     
		
    assert(disjuncts.size() == ndisjuncts);
		
		
    if(Result == UNKNOWN) {
      Solver s(model, disjuncts);
			
      if(s.status == UNKNOWN) {
				
	if( Rngd ) {
	  s.setRestartNogood();
	  s.setForgetfulness( 0.0 );
	}
	if( Randomized > 0 ) s.setRandomized();
	s.setTimeLimit( Cutoff );
	if(first_step) {
	  addHeuristic( s, Heuristic, Randomized, I_VO );
	} else
	  addHeuristic( s, Heuristic, Randomized, D_VO );
	s.setVerbosity( Verbose-1 );
	s.setRandomSeed( Seed );
				
	Result = UNKNOWN;
	if( SAC ) Result = s.sacPreprocess( SAC > 1 );
      }
      if( Result != UNSAT ) {
	if(first_step) {
	  if(I_VO == RAND)
	    s.setRandomValueOrdering();
	} else
	  switch(D_VO) {
	  case RAND: {
	    s.setRandomValueOrdering();
	  } break;
	  case GUIDED: {
	    s.setGuidedOrdering(disjuncts, best_solution);
	  } break;
	  case RGUIDED: {
	    s.setRandGuidedOrdering(disjuncts, best_solution, probability, range);
	  } break;
	  }
				
	if(Weight>0)
	  setWeights(s,tasks);
				
	Result = s.solve_and_restart(PolicyRestart, Base, Factor, Decay);
      }
			
      if( Result == SAT ) {
				
	updateBestSolution(disjuncts, tasks);
				
	new_makespan = (get_makespan(tasks));
	cout << left << setw(30) << "c solutions's makespan" << ":" << right << setw(20) << new_makespan << endl;
	dsopttime = (getRunTime() - total_time);
      } // else {
      // 	cout << left << setw(30) << "c " << ":" << right << setw(20) << "no solution";
      //       }
			
      if((Result != UNKNOWN) && (s.ENDTIME > maxidstime))
	maxidstime = s.ENDTIME;
			
      if(Weight>0)
	updateWeights(disjuncts, tasks);
			
			
      s.printStatistics(std::cout, (RUNTIME + BTS + PPGS + OUTCOME) );
			
      total_nodes += s.NODES;
      total_fails += s.FAILURES;
      total_bts += s.BACKTRACKS;
      total_propags += s.PROPAGS;
      total_restarts += s.BTSLIST.size;
			
      cout << "c =============[ end dichotomic step ]=============" << endl;
			
    } 
		
    cout << endl;
		
    if( Result == SAT ) {
      maxfsble = new_makespan;
    } else {
      minfsble = makespan+1;
      if( Result != SAT && Result != UNSAT )
	Solved = false;
      else
	max_infeasible = makespan;
    }
		
    first_step = false;
  }
	
  Solved = (maxfsble == max_infeasible+1);

  opt_time = dsopttime;
  dsrun_time = (getRunTime() - total_time);
  proof_time = (dsrun_time - opt_time);

	
  cout << "c =================[ statistics ]==================" << endl;
  cout << left << setw(30) << "d DSLOWERBOUND "    << right << setw(21) << (max_infeasible+1) << endl
       << left << setw(30) << "d DSOBJECTIVE "     << right << setw(21) << maxfsble << endl 
       << left << setw(30) << "d DSRUNTIME "       << right << setw(21) << dsrun_time << endl
       << left << setw(30) << "d DSOPTTIME "       << right << setw(21) << dsopttime << endl
       << left << setw(30) << "d DSMAXTIME "       << right << setw(21) << maxidstime << endl
       << left << setw(30) << "d DSNODES "         << right << setw(21) << total_nodes << endl
       << left << setw(30) << "d DSBACKTRACKS "    << right << setw(21) << total_bts << endl
       << left << setw(30) << "d DSFAILS "         << right << setw(21) << total_fails << endl
       << left << setw(30) << "d DSPROPAGS "       << right << setw(21) << total_propags << endl
       << left << setw(30) << "d DSRESTARTS "      << right << setw(21) << total_restarts << endl
	
       << left << setw(30) << "d DSSOLUTIONS "     << right << setw(21) << total_solutions << endl
       << left << setw(30) << "d AVGDISTANCE "     << right << setw(21) << avg_distance << endl
       << left << setw(30) << "d MINDISTANCE "     << right << setw(21) << min_distance << endl
       << left << setw(30) << "d MAXDISTANCE "     << right << setw(21) << max_distance << endl;
	
  if(Weight) 
    cout << left << setw(30) << "d AVGAVGDWEIGHT " << right << setw(21) << avg_disjunct_avg_weight << endl
	 << left << setw(30) << "d AVGMINDWEIGHT " << right << setw(21) << avg_disjunct_min_weight << endl
	 << left << setw(30) << "d AVGMAXDWEIGHT " << right << setw(21) << avg_disjunct_max_weight << endl
	 << left << setw(30) << "d MINAVGDWEIGHT " << right << setw(21) << min_disjunct_avg_weight << endl
	 << left << setw(30) << "d MINMINDWEIGHT " << right << setw(21) << min_disjunct_min_weight << endl
	 << left << setw(30) << "d MINMAXDWEIGHT " << right << setw(21) << min_disjunct_max_weight << endl
	 << left << setw(30) << "d MAXAVGDWEIGHT " << right << setw(21) << max_disjunct_avg_weight << endl
	 << left << setw(30) << "d MAXMINDWEIGHT " << right << setw(21) << max_disjunct_min_weight << endl
	 << left << setw(30) << "d MAXMAXDWEIGHT " << right << setw(21) << max_disjunct_max_weight << endl
		
	 << left << setw(30) << "d AVGAVGTWEIGHT " << right << setw(21) << avg_task_avg_weight << endl
	 << left << setw(30) << "d AVGMINTWEIGHT " << right << setw(21) << avg_task_min_weight << endl
	 << left << setw(30) << "d AVGMAXTWEIGHT " << right << setw(21) << avg_task_max_weight << endl
	 << left << setw(30) << "d MINAVGTWEIGHT " << right << setw(21) << min_task_avg_weight << endl
	 << left << setw(30) << "d MINMINTWEIGHT " << right << setw(21) << min_task_min_weight << endl
	 << left << setw(30) << "d MINMAXTWEIGHT " << right << setw(21) << min_task_max_weight << endl
	 << left << setw(30) << "d MAXAVGTWEIGHT " << right << setw(21) << max_task_avg_weight << endl
	 << left << setw(30) << "d MAXMINTWEIGHT " << right << setw(21) << max_task_min_weight << endl
	 << left << setw(30) << "d MAXMAXTWEIGHT " << right << setw(21) << max_task_max_weight << endl;      
	
  cout << left << setw(30) << "d DSOPTIMAL " << right << setw(21) << Solved << endl;
  cout << "c =================[ statistics ]==================" << endl
       << endl;
	
  ////// END DIHOTOMIC SEARCH /////////////
	
}

void branch_and_bound()
{
  int i, j;

  VarArray disjuncts;
  VarArray tasks;
  //VarArray endlasttask( (Type == "osp" ? nJobs*nMachines : nJobs), max_infeasible+1, maxfsble-1 );
  Variable endlasttask(max_infeasible+1, maxfsble-1);
  //new_makespan = maxfsble-1;

  CSP model;
  makespan = maxfsble-1;

  if(Type == "osp")
    osp_setup(model, disjuncts, tasks);
  else if(Type == "jtl0")
    jtl0_setup(model, disjuncts, tasks);     
  else if(Type == "pfsp")
    permfsp_setup(model, disjuncts, tasks);     
  else 
    jsp_setup(model, disjuncts, tasks);

  if( Type == "osp" )
    for(i=0; i<nJobs; ++i)
      for(j=0; j<nMachines; ++j) 
	model.add( tasks[i*nMachines+j] + duration[i][j] <= endlasttask );
  else if( Type == "jtl0" )
    for(i=0; i<nJobs; ++i) {
      int duration_ = 0;
      for(j=0; j<nMachines; ++j)
	duration_ += duration[i][j];
      model.add( tasks[i] + duration_ <= endlasttask );
    }
  else
    for(i=0; i<nJobs; ++i)
      model.add( tasks[(i+1)*nMachines-1] + duration[i][nMachines-1] <= endlasttask );
    
  model.add( Minimise(endlasttask) );

  Solver s(model, disjuncts);



  //if( Reset ) s.function = new ResetBase( &s );
  if(VO == GUIDED)
    s.function = new SolutionGuidedSearch( &s, disjuncts, Pool, Rprob );
  else if(VO == RGUIDED)
    s.function = new SolutionRandGuidedSearch( &s, disjuncts );

  if( Rngd ) {
    s.setRestartNogood();      
    s.setForgetfulness( 0.0 );
  }
  if( Randomized > 0 ) s.setRandomized();
  addHeuristic( s, Heuristic, Randomized, VO );
  s.setVerbosity( Verbose );
  //s.setVerbosity( 0 );
  s.setRandomSeed( Seed );

  double TIME = total_time + (Optimise - getRunTime());

  if( TIME > 2 ) {

    switch(VO) {
    case RAND: {
      s.setRandomValueOrdering();
    } break;
      //       case GUIDED: {
      // 	s.setGuidedOrdering(disjuncts, best_solution);
      //       } break;
      //       case RGUIDED: {
      // 	s.setRandGuidedOrdering(disjuncts, best_solution, probability, range);
      //       } break;
    }

    if(Weight>0)
      setWeights(s,tasks);

    cout << "c =============[ start branch&bound ]==============" << endl;
    cout 
      << left << setw(30) << "c initial bounds" << ":" 
      << right << setw(6) << " " << setw(5) << (max_infeasible+1) 
      << " to " << setw(5) << (maxfsble-1) << endl
      << left << setw(30) << "c search for (in sec.)" << ":" << right << setw(20) << TIME << endl;
    s.setTimeLimit( TIME );
    s.setVerbosity(Verbose);
    Result = s.solve_and_restart(PolicyRestart, Base, Factor, Decay);

    Solved = ( Result != LIMITOUT && Result != UNKNOWN );

    if( s.SOLUTIONS > 0 )
      maxfsble = s.goal->upper_bound;
    if( Solved )
      max_infeasible = maxfsble-1;

    if(s.SOLTIME > 0)
      opt_time = dsrun_time + s.SOLTIME;
    //if( Solved )
    //proof_time = s.ENDTIME - s.SOLTIME; 


    //cout << "c ";
    s.printStatistics(std::cout);
    //cout << endl;

    //s.printWeightProfile( cout, NOVAL, Wprof );
    //cout << endl;

    cout << "c ==============[ end branch&bound ]===============" << endl << endl;
      
    total_nodes += s.NODES;
    total_fails += s.FAILURES;
    total_bts += s.BACKTRACKS;
    total_propags += s.PROPAGS;
    total_restarts += s.BTSLIST.size;

    //// END BRANCH & BOUND ///////////////
  } else {
    //if(Verbose)
    cout << "d USEBNB 0" << endl;
  }    
}


 
void climbing_discrepancy_search()
{
  double TIME = (Optimise - getRunTime());

  cout << "c start climbing discrepancy search in [ " 
       << (max_infeasible+1) 
       << ".." << maxfsble-1 << " ] for " 
       << TIME << " seconds" << endl;

  while(max_infeasible+1 < maxfsble) {

    VarArray disjuncts;
    VarArray tasks;

    CSP model;
    makespan = maxfsble-1;

    if(Type == "osp")
      osp_setup(model, disjuncts, tasks);
    else 
      jsp_setup(model, disjuncts, tasks);

    Solver s(model, disjuncts);

    if( Randomized > 0 ) s.setRandomized();
    addHeuristic( s, Heuristic, Randomized, VO );
    s.setVerbosity( Verbose );
    s.setRandomSeed( Seed );

    switch(VO) {
    case RAND: {
      s.setRandomValueOrdering();
    } break;
    case GUIDED: {
      s.setGuidedOrdering(disjuncts, best_solution);
    } break;
    case RGUIDED: {
      s.setRandGuidedOrdering(disjuncts, best_solution, probability, range);
    } break;
    }
    
    if(Weight>0)
      setWeights(s,tasks);

    TIME = (Optimise - getRunTime());
    s.setTimeLimit( TIME );
    Result = s.ldSolve(disjuncts, best_solution, Step);

    if(Result == SAT) {
      updateBestSolution(disjuncts, tasks);
      
      //for(i=0; i<ndisjuncts; ++i)
      //best_solution[i] = disjuncts[i].value();
      maxfsble = (get_makespan(tasks));
      cout << "mkp: " << maxfsble << "  " << endl;
    } else Solved = (Result == UNSAT);

    cout << "c ";
    s.printStatistics(std::cout);
    cout << endl;      
      
    total_nodes += s.NODES;
    total_fails += s.FAILURES;
    total_bts += s.BACKTRACKS;
    total_propags += s.PROPAGS;
    total_restarts += s.BTSLIST.size;

  }
}







