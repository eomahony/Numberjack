
#include <mistral_sol.h>
#include <iostream>
#include <fstream>

using namespace std;

int bufsize = 10000;
int **duration, nMachines, nJobs, lb, opt;

int maxfsble;
int minfsble;
int makespan;
int max_infeasible;


long unsigned int total_nodes = 0;
long unsigned int total_bts = 0;
long unsigned int total_fails = 0;
long unsigned int total_propags = 0;
double total_time = 0.0;

int Seed, Cutoff, Dichotomy, Base, PIterations, PLimit, Probe, UpdateW, UpdateI, SAC, Randomized, ReuseWeight, Verbose, Step, Weight, Optimise, Rngd;
double Factor, Decay;
string Policy, Heuristic;

const int nia = 12;
const char* int_ident[nia] = {"-seed", "-cutoff", "-dichotomy", "-base", "-sac", "-randomized", "-reuse_weight", "-verbose", "-step", "-weight", "-optimise", "-restart_ngd"};
int int_param[nia];

const int nsa = 4;
const char* str_ident[nsa] = {"-heuristic", "-restart", "-factor", "-decay"};
const char* str_param[nsa];




void addHeuristic( Solver& s, string Heu, const int rdz ) {
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
  else if( Heu == "osp") {
    OSP h(abs(rdz));
    s.add( h );
  }
  else {
    NoOrder h;
    s.add( h );
  }
}



void readData( char* filename )
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
}


int upperbound() {
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


int lowerbound()
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

void setup(CSP& model, VarArray& disjuncts, VarArray& tasks) {
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




int main( int argc, char** argv )
{
  if( argc < 2 )
    {
      cerr << "need a data file" << endl;
      exit( 0 );
    }

  getCommandLine(int_ident,int_param,nia,str_ident,str_param,nsa,&argv[1],argc-1);

  Seed        = ( int_param[0]  != NOVAL ? int_param[0] : 11041979 );
  Cutoff      = ( int_param[1]  != NOVAL ? int_param[1] : 300);
  Dichotomy   = ( int_param[2]  != NOVAL ? int_param[2] : 64 );
  Base        = ( int_param[3]  != NOVAL ? int_param[3] : 64 );
  SAC         = ( int_param[4]  != NOVAL ? int_param[4] : 0  );
  Randomized  = ( int_param[5]  != NOVAL ? int_param[5] : 2  );
  ReuseWeight = ( int_param[6]  != NOVAL ? int_param[6] : 0  );
  Verbose     = ( int_param[7]  != NOVAL ? int_param[7] : 0  );
  Step        = ( int_param[8]  != NOVAL ? int_param[8] : 10 );
  Weight      = ( int_param[9]  != NOVAL ? int_param[9] : 0  );
  Optimise    = ( int_param[10] != NOVAL ? int_param[10]: 3000 );
  Rngd        = ( int_param[11] != NOVAL ? int_param[11]: 0 );

  Heuristic = ( strcmp(str_param[0],"nil") ? str_param[0] : "domowdeg" );
  Policy    = ( strcmp(str_param[1],"nil") ? str_param[1] : "geom"      );		
  Factor    = ( strcmp(str_param[2],"nil") ? atof(str_param[2]) : 1.2 );
  Decay     = ( strcmp(str_param[3],"nil") ? atof(str_param[3]) : 0.0 );


  //if(Verbose > 1) {
    cout << "c seed: " << Seed << endl;
    cout << "c cutoff: " << Cutoff << endl;
    cout << "c dichotomy: " << (Dichotomy ? "yes" : "no") << endl;
    cout << "c restart policy: " << Policy << endl;
    cout << "c base: " << Base << endl;
    cout << "c factor: " << Factor << endl;
    cout << "c sac: ";
    if(SAC == 0) cout << "no" ;
    else if(SAC == 1) cout << "approx" ;
    else if(SAC == 2) cout << "full" ;
    cout << endl;
    cout << "c reuse weights: " << (ReuseWeight ? "yes" : "no") << endl;
    cout << "c heuristic: " << Heuristic << " (" << abs(Randomized) << ")" << endl;
    //}

  int policy = NO;
  if( Policy == "geom" )
    policy = GEOMETRIC;
  else if( Policy == "luby" )
    policy = LUBY;

  int i, j;
  readData( argv[1] );


  total_time = getRunTime();


  maxfsble = upperbound();
  minfsble = lowerbound();
  //maxfsble = 311;
  //minfsble = 310;
  max_infeasible = minfsble-1;

  cout << "c " << minfsble << " <= " << lb << " <= " << opt << " <= " << maxfsble << endl;
  
  // Model

  bool solved = true;
  int N=0, M=0, *var_weights=NULL, *con_weights=NULL, k=0;
  double ratio=2, last_ratio=1;
  int result;

  while( Dichotomy-- && minfsble < maxfsble ) {
    //if( Dichotomy ) 
    makespan = (int)(floor(((double)minfsble + (double)maxfsble)/2));
    //else makespan = maxfsble-1;
      
    cout << "c \t" << setw(5) << minfsble << " " << setw(5) << maxfsble 
	 << " (" << setw(5) << makespan << "):\t";
    cout.flush();
      
    VarArray disjuncts;
    VarArray tasks;
    CSP model;
    setup(model, disjuncts, tasks);
      
    Solver s(model, disjuncts);
      
    //s.setLearner( Weighter::RNGD );
      
    if( Rngd ) {
      s.setRestartNogood();      
      s.setForgetfulness( 0.0 );
    }
    if( Randomized > 0 ) s.setRandomized();
    s.setTimeLimit( Cutoff );
    addHeuristic( s, Heuristic, Randomized );
    s.setVerbosity( Verbose );
    s.setRandomSeed( Seed );

    if(ReuseWeight) { 
      if(N) {
	if( ReuseWeight > 1 ) {
	  for(i=0; i<k; ++i) {
	    s.variables[var_weights[i]]->weight += (Weight+(ReuseWeight-i)*Step);
	  }
	} else {
	  for(i=0; i<N; ++i)
	    s.variables[i]->weight = var_weights[i];
	  for(i=0; i<M; ++i)
	    s.constraints[i]->weight = con_weights[i];
	}
      }
    }

    //     cout << endl;
    //     s.printWeightProfile( cout, 25 ); //s.length );

    result = UNKNOWN;
    if( SAC ) result = s.sacPreprocess( SAC );
    if( result != UNSAT )
      result = s.solve_and_restart(policy, Base, Factor, Decay);

    if( result == SAT ) {

      //       for(int i=0; i<disjuncts.size(); ++i)
      // 	cout << " " << ( disjuncts[i].value() ? i+1 : -i-1 );
      //       cout << endl;

      maxfsble = makespan;
    } else {
      minfsble = makespan+1;
      if( result != SAT && result != UNSAT )
	solved = false;
      else
	max_infeasible = makespan;
    }
    
    last_ratio = ((double)(s.BTSLIST.back())/(double)(s.FAILURES * log2(s.FAILURES)));
    
    //cout << "c ";
    s.printStatistics(std::cout, (RUNTIME + BTS + PPGS + OUTCOME + REST + LASTBTS) );
    cout << "  " << last_ratio << endl;

    total_nodes += s.NODES;
    total_fails += s.FAILURES;
    total_bts += s.BACKTRACKS;
    total_propags += s.PROPAGS;
    

    if(ReuseWeight) {      
      if(!var_weights) {
	N = s.length;
	M = s.constraints.size;

	if(ReuseWeight <= 1) {
	  var_weights = new int[N];
	  con_weights = new int[M];
	} else {
	  var_weights = new int[1+ReuseWeight];
	}
      }
      if( last_ratio < ratio && (result == SAT || result == UNSAT) ) {
	ratio = last_ratio;
	if(ReuseWeight > 1) {
	  k=0;
	  for(i=0; i<N; ++i) {
	    j = k;
	    var_weights[j] = i;
	    while( j && s.variables[var_weights[j]]->weight > s.variables[var_weights[j-1]]->weight ) {
	      var_weights[j] = var_weights[j-1];
	      var_weights[--j] = i;
	    }
	    if(k < ReuseWeight) ++k;
	  }
	  cout << "c ";
	  for(i=0; i<k; ++i) 
	    cout << setw(4) << var_weights[i] << " ";	   
	  cout << endl;
	} else {
	  for(i=0; i<N; ++i)
	    var_weights[i] = s.variables[i]->weight;
	  for(i=0; i<M; ++i)
	    con_weights[i] = s.constraints[i]->weight;
	}
      }
    }
  }

 

  solved = (maxfsble == max_infeasible+1);

  //if(Verbose) {
  cout << "d DSLOWERBOUND " << max_infeasible << endl
       << "d DSOBJECTIVE " << maxfsble << endl 
       << "d DSRUNTIME " << (getRunTime() - total_time) << endl
       << "d DSNODES " << total_nodes << endl
       << "d DSBACKTRACKS " << total_bts << endl
       << "d DSFAILS " << total_fails << endl
       << "d DSPROPAGS " << total_propags << endl
       << "d DSOPTIMAL " << solved << endl
       << endl;
  //}

  if( Optimise && !solved ) {
    VarArray disjuncts;
    VarArray tasks;
    VarArray endlasttask( nJobs*nMachines, max_infeasible+1, maxfsble-1 );
    
    CSP model;
    makespan = maxfsble-1;
    setup(model, disjuncts, tasks);
    for(i=0; i<nJobs; ++i)
      for(j=0; j<nMachines; ++j) 
 	model.add( tasks[i*nMachines+j] + duration[i][j] <= endlasttask[i*nMachines+j] );
    
    model.add( Minimise(Max(endlasttask)) );

    Solver s(model, disjuncts);


    s.function = new ResetBase( &s );

    if( Rngd ) {
      s.setRestartNogood();      
      s.setForgetfulness( 0.0 );
    }
    if( Randomized > 0 ) s.setRandomized();
    addHeuristic( s, Heuristic, Randomized );
    s.setVerbosity( Verbose );
    //s.setVerbosity( 0 );
    s.setRandomSeed( Seed );

    int TIME = total_time + (Optimise - getRunTime());

    if( TIME > 2 ) {
      
      cout << "c start branch&bound in [ " 
	   << (max_infeasible+1) 
	   << ".." << maxfsble-1 << " ] for " 
	   << TIME << " seconds" << endl << "d USEBNB 1" << endl;
      
      s.setTimeLimit( TIME );
      result = s.solve_and_restart(policy, Base, Factor, Decay);
      //result = s.solve();

 //      cout << "c " ;
//       switch( result ) {
//       case SAT: cout << "SAT"; break;
//       case UNSAT: cout << "UNSAT"; break;
//       case LIMITOUT: cout << "LIMITOUT"; break;
//       case UNKNOWN: cout << "UNKNOWN"; break;
//       }

      solved = ( result != LIMITOUT && result != UNKNOWN );

      if( s.SOLUTIONS > 0 )
	maxfsble = s.goal->upper_bound;
      //       else if( solved ) {
      // 	max_infeasible = s.goal->upper_bound;
      //       }

      if( solved )
	max_infeasible = maxfsble-1;

      
      cout << "c ";
      s.printStatistics(std::cout);
      cout << endl;
      
      total_nodes += s.NODES;
      total_fails += s.FAILURES;
      total_bts += s.BACKTRACKS;
      total_propags += s.PROPAGS;
    } else {
      //if(Verbose)
      cout << "d USEBNB 0" << endl;
    }
  } else {
    //if(Verbose)
    cout << "d USEBNB 0" << endl;
  }
  
  total_time = (getRunTime() - total_time);
  
  //if(Verbose) {
  cout << "s" << (solved ? " " : " UN" ) << "SATISFIABLE" << endl 
       << "d LOWERBOUND " << max_infeasible << endl
       << "d OBJECTIVE " << maxfsble << endl 
       << "d BESTKNOWN " << opt << endl
       << "d RUNTIME " << total_time << endl
       << "d NODES " << total_nodes << endl
       << "d BACKTRACKS " << total_bts << endl
       << "d FAILS " << total_fails << endl
       << "d PROPAGS " << total_propags << endl
       << "d OPTIMAL " << solved << endl
       << endl;
  //}
}



  


  




