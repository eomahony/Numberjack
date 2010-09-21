
#include <mistral_sat.h>
#include <assert.h>
#include <stdlib.h>
#include <fstream>
 
using namespace std;



const int nia = 6;
const char* int_ident[nia] = {
  "-h",
  "-verbose",
  "-time_limit",
  "-seed",
  "-randomize", 
  "-restart_base"
};

const char* int_man[nia] = {
  "help message",
  "verbosity     {0,1}", 
  "time limit",
  "random seed",
  "0 -> none, >0 -> shuffle vars on restart, k -> randomized heuristic", 
  "restart base"
   };

			      
int int_param[nia];

const int nsa = 4;
const char* str_ident[nsa] = {
  "-restart_factor",
  "-decay",
  "-forget",
  "-policy"
};
const char* str_man[nsa] = {
  "restart factor",
  "1 - activity decay rate",
  "forgetfulness ratio",
  "restart policy {luby, geom}"
};
const char* str_param[nsa];


void outputHelpMessage()
{
  cerr << "\nUsage: \t bin/sat [path to the cnf file] <args>" << endl << endl;
  for(int i=0; i<nia; ++i)
    cerr  << setw(20) << int_ident[i] << ": \t" << int_man[i] << endl;
  for(int i=0; i<nsa; ++i)
    cerr  << setw(20) << str_ident[i] << ": \t" << str_man[i] << endl;
  cerr << endl;
}


int main(int argc, char **argv)
{

  str_param[0] = "";
  getCommandLine(int_ident, int_param, nia,
		 str_ident, str_param, nsa,
		 &(argv[1]), argc-1);
  
  if( argc < 2 || int_param[0] != NOVAL || !strcmp( argv[argc-1], "-h" ) )
    {
      outputHelpMessage();
    }
  else
    {
      SatSolver solver;

      solver.params.verbosity      = ( int_param[1] != NOVAL ? int_param[1] : 4 );
      solver.params.time_limit     = ( (double)(int_param[2] != NOVAL ? int_param[2] : -1 ) );
      solver.params.seed           = ( int_param[3] != NOVAL ? int_param[3] : 11041979 );	  
      solver.params.randomization  = ( int_param[4] != NOVAL ? abs(int_param[4]) : 2 );
      if(solver.params.randomization == 0)
	solver.params.randomization = 1;
      solver.params.shuffle        = ( int_param[4] > 0 );
      solver.params.restart_base   = ( int_param[5] != NOVAL ? int_param[5] : 200 );
      solver.params.restart_limit  = ( int_param[5] != NOVAL ? int_param[5] : 200 );
      solver.params.restart_factor = ( strcmp(str_param[0],"nil") ? atof(str_param[0]) : 1.05 );
      solver.params.decay          = ( strcmp(str_param[1],"nil") ? atof(str_param[1]) : .96 );
      solver.params.forgetfulness  = ( strcmp(str_param[2],"nil") ? atof(str_param[2]) : .75 );
      solver.setPolicy             ( ( strcmp(str_param[3],"luby") ? GEOMETRIC : LUBY ) );


      solver.parseDimacs(argv[1]);
      solver.solve();
    }
  return 0;
}



