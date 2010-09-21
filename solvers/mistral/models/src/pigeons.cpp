

#include <mistral_sol.h>

/* 
 * a template model
 */

 
using namespace std;

int main(int argc, char *argv[])
{  
  int N = (argc > 1 ? atoi(argv[1]) : 8);

  // create the csp
  CSP model;
  VarArray pigeons(N, 1, N-1);
  
  for(int i=0; i<N; ++i)
    for(int j=i+1; j<N; ++j) 
      model.add( pigeons[i] != pigeons[j] );

  Solver s( model );

  s.print(cout);
  cout << endl;

  s.solve();

  s.printStatistics( cout );
  cout << endl;
}


