
#include <mistral_sol.h>
#include <iostream>
#include <fstream>

using namespace std;

int bufsize = 10000;
int *duration, *machine, nMachines, nJobs, lb, opt;

int maxfsble;
int minfsble;
int max_infeasible;
int makespan;



long unsigned int total_nodes = 0;
long unsigned int total_bts = 0;
long unsigned int total_fails = 0;
long unsigned int total_propags = 0;
double total_time = 0.0;

int Seed, Cutoff, Dichotomy, Base, PIterations, PLimit, Probe, UpdateW, UpdateI, SAC, Randomized, ReuseWeight, Verbose, Step, Weight, Optimise;
double Factor, Decay;
string Policy, Heuristic;

const int nia = 11;
const char* int_ident[nia] = {"-seed", "-cutoff", "-dichotomy", "-base", "-sac", "-randomized", "-reuse_weight", "-verbose", "-step", "-weight", "-optimise"};
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
  int i, j;
  long int dump;
  string tag;
  ifstream infile( filename, ios_base::in );
  
  infile >> nJobs;
  infile >> nMachines;

  duration = new int[nJobs*nMachines];
  machine = new int[nJobs*nMachines];

  infile >> dump;
  infile >> dump;

  infile >> opt ;
  infile >> lb;

  infile >> tag;
  assert( tag == "Times" );

  for(i=0; i<nJobs; ++i) 
    for(j=0; j<nMachines; ++j) 
      infile >> duration[i*nMachines+j];

  infile >> tag;
  assert( tag == "Machines" );

  for(i=0; i<nJobs; ++i) 
    for(j=0; j<nMachines; ++j) {
      infile >> machine[i*nMachines+j];
      --machine[i*nMachines+j];
    }
  
}

 
int upperbound() {
  int i, j, mkp=0, task, boundJob[nJobs], boundMachine[nJobs]; 
  std::fill(boundJob, boundJob+nJobs, 0);
  std::fill(boundMachine, boundMachine+nJobs, 0);


  for(i=0; i<nMachines; ++i)
    {
      for(j=0; j<nJobs; ++j)
	{
	  task = std::max( boundJob[j], boundMachine[machine[j*nMachines+i]] );
	  task += duration[j*nMachines+i];
	  boundJob[j] = task;
	  boundMachine[machine[j*nMachines+i]] = task;
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
	boundJob[i] += duration[i*nMachines+j];
	boundMachine[machine[i*nMachines+j]] += duration[i*nMachines+j];
      }
  }

  for(i=0; i<nJobs; ++i)
    if( mkp < boundJob[i] ) mkp = boundJob[i];
  for(i=0; i<nMachines; ++i)
    if( mkp < boundMachine[i] ) mkp = boundMachine[i];

  return mkp;
}

void setup(CSP& model, VarArray& disjuncts, VarArray& tasks) {
  int i,j,k, per_machine[nJobs*nMachines];

  // tasks' start-time
  for(i=0; i<nJobs; ++i) 
    for(j=0; j<nMachines; ++j) {
      tasks.add( Variable(0, makespan-duration[i*nMachines+j]) );
    }
 
  // order within a job
  for(i=0; i<nJobs; ++i)
    for(j=1; j<nMachines; ++j) {
      model.add( Precedence(tasks[i*nMachines+j-1], duration[i*nMachines+j-1], tasks[i*nMachines+j]) );
    }

  // non-overlapping reified constraints
  for(i=0; i<nJobs; ++i)
    for(j=0; j<nMachines; ++j) 
      per_machine[i*nMachines+machine[i*nMachines+j]] = (i*nMachines+j);
  for(k=0; k<nMachines; ++k) 
    for(i=0; i<nJobs; ++i)
      for(j=i+1; j<nJobs; ++j) {       
	disjuncts.add( Disjunctive( tasks[per_machine[i*nMachines+k]], 
				    duration[per_machine[i*nMachines+k]],
				    tasks[per_machine[j*nMachines+k]], 
				    duration[per_machine[j*nMachines+k]] ) );
      }

  model.add( disjuncts );
}


int main( int argc, char** argv )
{
  if( argc < 2 )
    {
      cerr << "need a data file" << endl;
      exit( 0 );
    }



  getCommandLine(int_ident,int_param,nia,str_ident,str_param,nsa,&argv[1],argc-1);

  Seed        = ( int_param[0]  != NOVAL ? int_param[0] : 11041979 );
  Cutoff      = ( int_param[1]  != NOVAL ? int_param[1] : 300 );
  Dichotomy   = ( int_param[2]  != NOVAL ? int_param[2] : 1   );
  Base        = ( int_param[3]  != NOVAL ? int_param[3] : 64  );
  SAC         = ( int_param[4]  != NOVAL ? int_param[4] : 0   );
  Randomized  = ( int_param[5]  != NOVAL ? int_param[5] : 2   );
  ReuseWeight = ( int_param[6]  != NOVAL ? int_param[6] : 0   );
  Verbose     = ( int_param[7]  != NOVAL ? int_param[7] : 0   );
  Step        = ( int_param[8]  != NOVAL ? int_param[8] : 10  );
  Weight      = ( int_param[9]  != NOVAL ? int_param[9] : 0   );
  Optimise    = ( int_param[10] != NOVAL ? int_param[10]: 3000 );

  Heuristic = ( strcmp(str_param[0],"nil") ? str_param[0] : "domowdeg" );
  Policy    = ( strcmp(str_param[1],"nil") ? str_param[1] : "geom"      );		
  Factor    = ( strcmp(str_param[2],"nil") ? atof(str_param[2]) : 1.2 );
  Decay     = ( strcmp(str_param[3],"nil") ? atof(str_param[3]) : 0.0 );


    //  if(Verbose > 1) {
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
    //  }

  int policy = NO;
  if( Policy == "geom" )
    policy = GEOMETRIC;
  else if( Policy == "luby" )
    policy = LUBY;


//   int randc = ( argc >= 3 ? atoi(argv[2]) : 1 );
//   int dichotomy_cutoff = ( argc >= 4 ? atoi(argv[3]) : 30 );
//   int total_cutoff = ( argc >= 5 ? atoi(argv[4]) : 1800 );

  total_time = getRunTime();



  int i, j;
  readData( argv[1] );

  cout << "c "<< nJobs << " x " << nMachines << endl;

  maxfsble = upperbound();
  minfsble = lowerbound();
  max_infeasible = minfsble-1;

  cout << "c " << minfsble << " <= " << lb << " <= " << opt << " <= " << maxfsble << endl;


  bool solved = false;
  int result, N=0, M=0, *var_weights=NULL, *con_weights=NULL, k=0;
  double ratio=2, last_ratio=1;



  while( minfsble < maxfsble ) {
    if( Dichotomy ) makespan = (int)(floor(((double)minfsble + (double)maxfsble)/2));
    else makespan = maxfsble-1;

    cout << "c\t" << setw(5) << minfsble << " " << setw(5) << maxfsble 
	 << " (" << setw(5) << makespan << "):\t";
    cout.flush();
                      

    // model and solve
    VarArray disjuncts;
    VarArray tasks;
    CSP model;
    setup(model, disjuncts, tasks);
    Solver s(model, disjuncts);

    if( Randomized ) s.setRandomized();
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

    //int result = s.solve_and_restart(policy, Base, Factor, Decay);

    // update minfsble and maxfsble
    if( result == SAT ) {
      maxfsble = makespan;
    } else {
      minfsble = makespan+1;
      if( result != SAT && result != UNSAT )
	solved = false;
      else 
	max_infeasible = makespan;
    }

    last_ratio = ((double)(s.BTSLIST.back())/(double)(s.BACKTRACKS * log2(s.BACKTRACKS)));

    //cout << "c ";
    s.printStatistics(std::cout, (RUNTIME + BTS + PPGS + OUTCOME + REST + LASTBTS ) );
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


  cout << "d DSLOWERBOUND " << max_infeasible << endl
       << "d DSOBJECTIVE " << maxfsble << endl 
       << "d DSRUNTIME " << (getRunTime() - total_time) << endl
       << "d DSNODES " << total_nodes << endl
       << "d DSBACKTRACKS " << total_bts << endl
       << "d DSFAILS " << total_fails << endl
       << "d DSPROPAGS " << total_propags << endl
       << "d DSOPTIMAL " << solved << endl
       << endl;


  // if not solved, start a complete branch & bound search
  if( Optimise && !solved ) {

    VarArray disjuncts;
    VarArray tasks;
    VarArray endlasttask( nJobs, 0, maxfsble-1 );

    CSP model;
    makespan = maxfsble-1;
    setup(model, disjuncts, tasks);
    for(i=0; i<nJobs; ++i)
      model.add( tasks[(i+1)*nMachines-1] + duration[(i+1)*nMachines-1] <= endlasttask[i] );
    model.add( Minimise(Max(endlasttask)) );


    //DomOverWDeg heuristic(abs(Randomized));
    Solver s(model, disjuncts);
    if( Randomized ) s.setRandomized();
    addHeuristic( s, Heuristic, Randomized );
    s.setVerbosity( 0 );
    s.setRandomSeed( Seed );
    s.setRandomized();

    int TIME = total_time + (Optimise - getRunTime());

    if( TIME > 2 ) {

      cout << "c start branch&bound in [ " 
	   << (max_infeasible+1) 
	   << ".." << maxfsble-1 << " ] for " 
	   << TIME << " seconds" << endl << "d USEBNB 1" << endl;
      
      s.setTimeLimit( TIME );
      result = s.solve();

  //     cout << "c " ;
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
     cout << "d USEBNB 0" << endl;
    }
  } else {
    cout << "d USEBNB 0" << endl;
  }
    
  total_time = (getRunTime() - total_time);

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


}



  


  




