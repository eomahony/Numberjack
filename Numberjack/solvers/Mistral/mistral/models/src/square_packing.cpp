
#include <mistral_sol.h>
#include <iostream>
#include <fstream>

using namespace std;


int maxfsble;
int minfsble;
int makespan;
int max_infeasible;
int M;


long unsigned int total_nodes = 0;
long unsigned int total_bts = 0;
long unsigned int total_fails = 0;
long unsigned int total_propags = 0;
double total_time = 0.0;

int Seed, Cutoff, Dichotomy, Base, PIterations, PLimit, Probe, UpdateW, UpdateI, 
  SAC, Randomized, ReuseWeight, Verbose, Step, Weight, Optimise, Rngd, Gap, Reset;
double Factor, Decay;
string Policy, Heuristic, Type, ValueO;

const int nia = 14;
const char* int_ident[nia] = {"-seed", "-cutoff", "-dichotomy", "-base", "-sac", "-randomized", "-reuse_weight", "-verbose", "-step", "-weight", "-optimise", "-restart_ngd", "-gap", "-reset"};
int int_param[nia];

const int nsa = 6;
const char* str_ident[nsa] = {"-heuristic", "-restart", "-factor", "-decay", "-type", "-vo"};
const char* str_param[nsa];




void addHeuristic( Solver& s, string Heu, const int rdz ) {

  if(ValueO == "lex") {
    for(int i=0; i<s.length; ++i)
      s.sequence[i]->branch = new ValSelectorMin( s.sequence[i] );
  } else if(ValueO == "rand") {
    for(int i=0; i<s.length; ++i)
      s.sequence[i]->branch = new ValSelectorRand( s.sequence[i] );
  } else if(ValueO == "antilex") {
    for(int i=0; i<s.length; ++i)
      s.sequence[i]->branch = new ValSelectorMax( s.sequence[i] );
  } 
  
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
  else if( Heu == "fpp") {
    FPP h(abs(rdz));
    s.add( h );
  }
  else {
    NoOrder h;
    s.add( h );
  }
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
};




void setup(CSP& model, VarArray& precedence, VarArray& X, VarArray& Y) {
  int i,j,k;

  model.add( precedence );

  for(i=1; i<=M; ++i) {
    X.add( Variable(0, makespan-i) );
    Y.add( Variable(0, makespan-i) );
  }

  k=0;
  for(i=0; i<M; ++i) {
    for(j=i+1; j<M; ++j) {
      model.add( precedence[k++] == Precedence( X[i], (i+1), X[j] ) );
      model.add( precedence[k++] == Precedence( X[j], (j+1), X[i] ) );
      model.add( precedence[k++] == Precedence( Y[i], (i+1), Y[j] ) );
      model.add( precedence[k++] == Precedence( Y[j], (j+1), Y[i] ) );
    }
  }

  for(i=0; i<precedence.size(); i+=4) {
    VarArray overlap;
    for(j=0; j<4; ++j)
      overlap.add(precedence[i+j]);
    model.add( Sum(overlap) > 0 );
  }
}




int main( int argc, char** argv )
{
  if( argc < 2 )
    {
      cerr << "need a number" << endl;
      exit( 0 );
    }

  getCommandLine(int_ident,int_param,nia,str_ident,str_param,nsa,&argv[1],argc-1);

  Seed        = ( int_param[0]  != NOVAL ? int_param[0] : 11041979 );
  Cutoff      = ( int_param[1]  != NOVAL ? int_param[1] : 60 );
  Dichotomy   = ( int_param[2]  != NOVAL ? int_param[2] : 64  );
  Base        = ( int_param[3]  != NOVAL ? int_param[3] : 256 );
  SAC         = ( int_param[4]  != NOVAL ? int_param[4] : 0   );
  Randomized  = ( int_param[5]  != NOVAL ? int_param[5] : 1   );
  ReuseWeight = ( int_param[6]  != NOVAL ? int_param[6] : 0   );
  Verbose     = ( int_param[7]  != NOVAL ? int_param[7] : 0   );
  Step        = ( int_param[8]  != NOVAL ? int_param[8] : 10  );
  Weight      = ( int_param[9]  != NOVAL ? int_param[9] : 0   );
  Optimise    = ( int_param[10] != NOVAL ? int_param[10]: 3600 );
  Rngd        = ( int_param[11] != NOVAL ? int_param[11]: 0   );
  Gap         = ( int_param[12] != NOVAL ? int_param[12]: 0   );
  Reset       = ( int_param[13] != NOVAL ? int_param[13]: 0   );

  Heuristic = ( strcmp(str_param[0],"nil") ? str_param[0] : "domowdeg" );
  Policy    = ( strcmp(str_param[1],"nil") ? str_param[1] : "geom"      );		
  Factor    = ( strcmp(str_param[2],"nil") ? atof(str_param[2]) : 1.3 );
  Decay     = ( strcmp(str_param[3],"nil") ? atof(str_param[3]) : 0.0 );
  Type      = ( strcmp(str_param[4],"nil") ? str_param[4] : "jsp" );
  ValueO    = ( strcmp(str_param[5],"nil") ? str_param[5] : "rand" );

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

  int i;
  M = atoi( argv[1] );

  total_time = getRunTime();

  int total_surface = 0;
  for(i=1; i<=M; ++i)
    total_surface += (i*i);

  minfsble = sqrt(total_surface);
  maxfsble = 3*minfsble/2;
  
  // Model
  bool solved = true;
  int result;

  while( Gap < (maxfsble-minfsble) && Dichotomy-- && minfsble < maxfsble ) {
    makespan = (int)(floor(((double)minfsble + (double)maxfsble)/2));
      
    cout << "c \t" << setw(5) << minfsble << " " << setw(5) << maxfsble 
	 << " (" << setw(5) << makespan << "):\t";
    cout.flush();
      
    VarArray X;
    VarArray Y;
    VarArray precedence(M*(M-1)*2, 0, 1);
    CSP model;
    setup(model, precedence, X, Y);
      
    Solver s(model, precedence);

//     cout << endl << s.status << " " << UNSAT << endl;
//      s.print( cout );
//      cout << endl;

// //     //exit(0);

    result = s.status;
    if( s.status == UNKNOWN ) {

      if( Rngd ) {
	s.setRestartNogood();      
	s.setForgetfulness( 0.0 );
      }
      if( Randomized > 0 ) s.setRandomized();
      s.setTimeLimit( Cutoff );
      addHeuristic( s, Heuristic, Randomized );
      s.setVerbosity( Verbose );
      s.setRandomSeed( Seed );
     
      if( SAC ) result = s.sacPreprocess( SAC );
      if( result != UNSAT )
	result = s.solve_and_restart(policy, Base, Factor, Decay);
    }      


    if( result == SAT ) {
      maxfsble = makespan;
    } else {
      minfsble = makespan+1;
      if( result != SAT && result != UNSAT )
	solved = false;
      else
	max_infeasible = makespan;
    }

    s.printStatistics(std::cout, (RUNTIME + BTS + PPGS + OUTCOME + REST + LASTBTS) );
    cout << endl;

    total_nodes += s.NODES;
    total_fails += s.FAILURES;
    total_bts += s.BACKTRACKS;
    total_propags += s.PROPAGS;
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
  
  total_time = (getRunTime() - total_time);


  if( Optimise && !solved ) {
    makespan = maxfsble-1;

    VarArray X;
    VarArray Y;
    VarArray precedence(M*(M-1)*2, 0, 1);
    Variable EdgeSize(max_infeasible+1, maxfsble-1);
    CSP model;
    setup(model, precedence, X, Y);
    for(i=0; i<M; ++i) {
      model.add( Precedence(X[i], (i+1), EdgeSize) );
      model.add( Precedence(Y[i], (i+1), EdgeSize) );
    }
    model.add( Minimise(EdgeSize) );

    Solver s(model, precedence);

//     s.print( cout );
//     cout << endl;
//     exit(0);

    result = s.status;
    if( s.status == UNKNOWN ) {
      if( Rngd ) {
	s.setRestartNogood();      
	s.setForgetfulness( 0.0 );
      }
      if( Randomized > 0 ) s.setRandomized();
      addHeuristic( s, Heuristic, Randomized );
      s.setVerbosity( Verbose );
      s.setRandomSeed( Seed );
     
      if( SAC ) result = s.sacPreprocess( SAC );
      if( result != UNSAT )
	result = s.solve_and_restart(policy, Base, Factor, Decay);

      double TIME = total_time + (Optimise - getRunTime());

      if( TIME > 2 ) {
      
	cout << "c start branch&bound in [ " 
	     << (max_infeasible+1) 
	     << ".." << maxfsble-1 << " ] for " 
	     << TIME << " seconds" << endl << "d USEBNB 1" << endl;
      
	s.setTimeLimit( TIME );

	result = s.solve_and_restart(policy, Base, Factor, Decay);
	solved = ( result != LIMITOUT && result != UNKNOWN );

	if( s.SOLUTIONS > 0 )
	  maxfsble = s.goal->upper_bound;
	if( solved )
	  max_infeasible = maxfsble-1;
	
	total_nodes += s.NODES;
	total_fails += s.FAILURES;
	total_bts += s.BACKTRACKS;
	total_propags += s.PROPAGS;
      } 
    }      
  }
  
  //if(Verbose) {
  cout << "s" << (solved ? " " : " UN" ) << "SATISFIABLE" << endl 
       << "d LOWERBOUND " << max_infeasible << endl
       << "d OBJECTIVE " << maxfsble << endl 
       << "d RUNTIME " << total_time << endl
       << "d NODES " << total_nodes << endl
       << "d BACKTRACKS " << total_bts << endl
       << "d FAILS " << total_fails << endl
       << "d PROPAGS " << total_propags << endl
       << "d OPTIMAL " << solved << endl
       << endl;
  //}
}



  


  




