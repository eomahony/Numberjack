
#include <mistral_sol.h>

/* 
 * a template model
 */

 
using namespace std;

int Seed, Cutoff, Base, SAC, Randomized, Verbose;
double Factor, Decay;
string Policy, VarOrdering, ValOrdering;

const int nia = 6;
const char* int_ident[nia] = {"-seed", "-cutoff", "-base", "-sac", 
			      "-randomized", "-verbose"};
int int_param[nia];

const int nsa = 5;
const char* str_ident[nsa] = {"-var_order", "-algo", "-factor", 
			      "-decay",  "-val_order"};
const char* str_param[nsa];

Contention *proba;


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
  else {
    NoOrder h;
    s.add( h );
  }
}

void addBranching( Solver& s, string Heu) {
  if( Heu == "lex" ) {
    s.setLex();
  } else if( Heu == "antilex" ) {
    s.setAntiLex();
  } else if( Heu == "random" ) {
    s.setSplit();
  } else if( Heu == "randminmax" ) {
    s.setRandMinMax();
  } else if( Heu == "split" ) {
    s.setSplit();
  } else if( Heu == "randsplit" ) {
    s.setRandSplit();
  } else if( Heu == "contention" ) {
    s.add(*proba);
  } 
}

int get_data(char *filename, int**& squares) {
  ifstream infile(filename, ios_base::in);
  string buf;
  int i, j, order;
  
  infile >> buf;
  assert(buf == "order");

  infile >> order;
  squares = new int*[order];
  for(i=0; i<order; ++i) {
    squares[i] = new int[order];
    for(j=0; j<order; ++j) 
      infile >> squares[i][j];
  }

  return order;
}


int main(int argc, char *argv[])
{
  if(argc > 1) {

    getCommandLine(int_ident,int_param,nia,str_ident,str_param,nsa,&argv[1],argc-1);

    Seed        = ( int_param[0]  != NOVAL ? int_param[0] : 11041979 );
    Cutoff      = ( int_param[1]  != NOVAL ? int_param[1] : 60 );
    Base        = ( int_param[2]  != NOVAL ? int_param[2] : 256 );
    SAC         = ( int_param[3]  != NOVAL ? int_param[3] : 0   );
    Randomized  = ( int_param[4]  != NOVAL ? int_param[4] : 2   );
    Verbose     = ( int_param[5]  != NOVAL ? int_param[5] : 1   );

    VarOrdering = ( strcmp(str_param[0],"nil") ? str_param[0] : "dom" );
    Policy      = ( strcmp(str_param[1],"nil") ? str_param[1] : "geom"      );	
    Factor      = ( strcmp(str_param[2],"nil") ? atof(str_param[2]) : 1.3 );
    Decay       = ( strcmp(str_param[3],"nil") ? atof(str_param[3]) : 0.0 );
    ValOrdering = ( strcmp(str_param[4],"nil") ? str_param[4] : "randminmax" );

    int **squares;
    int i, j, N=get_data(argv[1], squares);
  
    // Create the csp
    CSP model;
    VarArray X(N*N, 0, N-1);
    VarArray rows[N];
    VarArray cols[N];
    
    model.add( X );
    
    for(i=0; i<N; ++i) {
      VarArray row( N );
      VarArray col( N );
      for(j=0; j<N; ++j) {
	row[j] = X[i*N+j];
	col[j] = X[j*N+i];
      }
      model.add( AllDifferent(row) );
      model.add( AllDifferent(col) );
      rows[i] = row;
      cols[i] = col;
    }
 

  for(i=0; i<N; ++i) 
    for(j=0; j<N; ++j) 
      if(squares[i][j] >= 0)
	model.add( rows[i][j] == squares[i][j] );

  int Result = UNKNOWN;
  Solver s( model );

  proba = new Contention(rows, N);

  s.setRandomSeed(Seed);
  s.setVerbosity(Verbose);
  if( Randomized > 0 ) s.setRandomized();
  
  if( SAC ) Result = s.sacPreprocess( SAC > 1 );

  if( Result == UNKNOWN ) {

    addHeuristic(s, VarOrdering, Randomized);
    addBranching(s, ValOrdering);

    if( Policy == "geom" )
      Result = s.solve_and_restart(GEOMETRIC, Base, Factor, Decay);
    else if( Policy == "luby" )
      Result = s.solve_and_restart(LUBY, Base, Factor, Decay);
    else if( Policy == "lds" ) {
      cout << "to do" << endl;
      exit(1);
    } else if( Policy == "lds" ) {
      cout << "to do" << endl;
      exit(1);
    } else if( Policy == "lds-lvl" ) {
      cout << "to do" << endl;
      exit(1);
    } else if( Policy == "lds-var" ) {
      cout << "to do" << endl;
      exit(1);
    } else if( Policy == "rds" ) {
      cout << "to do" << endl;
      exit(1);
    } else if( Policy == "rds-lvl" ) {
      cout << "to do" << endl;
      exit(1);
    } else if( Policy == "lds-var" ) {
      cout << "to do" << endl;
      exit(1);
    } else if( Policy == "dds" ) {
      cout << "to do" << endl;
      exit(1);
    } else if( Policy == "dds-lvl" ) {
      cout << "to do" << endl;
      exit(1);
    } else if( Policy == "dds-var" ) {
      cout << "to do" << endl;
      exit(1);
    } else
      Result = s.solve();
  }

  if(Result == SAT )
    {
      for(i=0; i<N; ++i) {
	for(j=0; j<N; ++j) 
	  cout << setw(3) << (rows[i][j].value());
	cout << endl;
      }
    }

  s.printStatistics( cout );
  cout << endl;
  }
}


