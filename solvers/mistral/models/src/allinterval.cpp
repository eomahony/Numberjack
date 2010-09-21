
#include <mistral_sol.h>

/* 
 * a template model (All Interval serie)
 */

using namespace std;

int main(int argc, char *argv[])
{  

  //cout << (-3 / 32) << endl
  //<< (-3 & 31) << endl;

  //exit( 0 );

  int N = (argc > 1 ? atoi(argv[1]) : 100);

  // declare a CSP
  CSP model;
  VarArray X(N, 0, N-1);
  VarArray distance(N-1, 1, N-1);

  // post the constraints
  for(int i=1; i<N; ++i)
    model.add( distance[i-1] == Abs(X[i-1] - X[i]) );
  model.add( AllDifferent(X) );
  model.add( AllDifferent(distance) );

  // initialise the solver
  Lexicographic heuristic;
  //MinDomain heuristic;
  //DomOverWDeg heuristic;
  Solver s( model, X, heuristic );
  s.solve();

  // print the result
  cout << X << endl;
  s.printStatistics( cout );
  cout << endl;
}


