
#include <mistral_sol.h>

/* 
 * a model of the costas array problem
 */

using namespace Mistral;
using namespace std;

int main(int argc, char *argv[])
{  
  /// input
  int i, j, N = ( argc > 1 ? atoi(argv[1]) : 11 );

  /// vars
  CSP model;
  VarArray X( N, 1, N );
  VarArray differences[N-2];

  /// constraints
  model.add( AllDifferent( X ) );  
  for(i=1; i<N-1; ++i) {
    for(j=i; j<N; ++j) {
      differences[i-1].add( (X[j-i] - X[j]) );
      model.add( differences[i-1][j-i] != 0 );
    }
    model.add( AllDifferent( differences[i-1] ) );
  }

  //model.add( X[(N+1)/2] >= N/2 );

  /// search
  ImpactOverWDeg heuristic(2);
  //MinDomain heuristic(3);
  Solver s(model, heuristic);
  s.setVerbosity(1);
  s.solve_and_restart(GEOMETRIC, 64, 1.33);

//   /// search
//   Lexicographic heuristic;
//   Solver s(model, X, heuristic);
//   s.FIND_ALL = -1;
//   s.setVerbosity(4);
//   s.solve();

  /// output
  cout << endl << endl;    
  for(i=0; i<N; ++i)
    cout << setw(3) << X[i].value() ;
  cout << endl ;
  for(i=0; i<N-2; ++i) {
    for(j=0; j<(differences[i].size()); ++j)
      cout << setw(3) << differences[i][j].value() ;
    cout << endl;
  }
  cout << setw(3) << (X[0].value() - X[N-1].value()) << endl << endl;
  s.printStatistics(std::cout);
  cout<<endl<<endl;
}



//  16 13  4 18 12 15  8  7  3  5 17  6  1 10 11  9 14  2
//   3  9-14  6 -3  7  1  4 -2-12 11  5 -9 -1  2 -5 12
//  12 -5 -8  3  4  8  5  2-14 -1 16 -4-10  1 -3  7
//  -2  1-11 10  5 12  3-10 -3  4  7 -5 -8 -4  9
//   4 -2 -4 11  9 10 -9  1  2 -5  6 -3-13  8
//   1  5 -3 15  7 -2  2  6 -7 -6  8 -8 -1
//   8  6  1 13 -5  9  7 -3 -8 -4  3  4
//   9 10 -1  1  6 14 -2 -4 -6 -9 15
//  13  8-13 12 11  5 -3 -2-11  3
//  11 -4 -2 17  2  4 -1 -7  1
//  -1  7  3  8  1  6 -6  5
//  10 12 -6  7  3  1  6
//  15  3 -7  9 -2 13
//   6  2 -5  4 10
//   5  4-10 16
//   7 -1  2
//   2 11
//  14

//       SAT  1270166 BTS  1270441 NDS  1270218 FAILS     4298 BTS/s     2187427 PPGS/s     646319143 PPGS             0 CKS  295.47 s  295.47 s



//  15  7 18  8  2  6 11 14 16 13 12  3  4 17 10  5  1  9 19
//   8-11 10  6 -4 -5 -3 -2  3  1  9 -1-13  7  5  4 -8-10
//  -3 -1 16  2 -9 -8 -5  1  4 10  8-14 -6 12  9 -4-18
//   7  5 12 -3-12-10 -2  2 13  9 -5 -7 -1 16  1-14
//  13  1  7 -6-14 -7 -1 11 12 -4  2 -2  3  8 -9
//   9 -4  4 -8-11 -6  8 10 -1  3  7  2 -5 -2
//   4 -7  2 -5-10  3  7 -3  6  8 11 -6-15
//   1 -9  5 -4 -1  2 -6  4 11 12  3-16
//  -1 -6  6  5 -2-11  1  9 15  4 -7
//   2 -5 15  4-15 -4  6 13  7 -6
//   3  4 14 -9 -8  1 10  5 -3
//  12  3  1 -2 -3  5  2 -5
//  11-10  8  3  1 -3 -8
//  -2 -3 13  7 -7-13
//   5  2 17 -1-17
//  10  6  9-11
//  14 -2 -1
//   6-12
//  -4

//       SAT  7034549 BTS  7034912 NDS  7034610 FAILS     4532 BTS/s     2369542 PPGS/s    3677506804 PPGS             0 CKS    1552 s    1552 s



//  11 10  6 16 19  4 20 12 14  8  3  1 15  5  2  7 13 17 18  9
//   1  4-10 -3 15-16  8 -2  6  5  2-14 10  3 -5 -6 -4 -1  9
//   5 -6-13 12 -1 -8  6  4 11  7-12 -4 13 -2-11-10 -5  8
//  -5 -9  2 -4  7-10 12  9 13 -7 -2 -1  8 -8-15-11  4
//  -8  6-14  4  5 -4 17 11 -1  3  1 -6  2-12-16 -2
//   7-10 -6  2 11  1 19 -3  9  6 -4-12 -2-13 -7
//  -9 -2 -8  8 16  3  5  7 12  1-10-16 -3 -4
//  -1 -4 -2 13 18-11 15 10  7 -5-14-17  6
//  -3  2  3 15  4 -1 18  5  1 -9-15 -8
//   3  7  5  1 14  2 13 -1 -3-10 -6
//   8  9 -9 11 17 -3  7 -5 -4 -1
//  10 -5  1 14 12 -9  3 -6  5
//  -4  5  4  9  6-13  2  3
//   6  8 -1  3  2-14 11
//   9  3 -7 -1  1 -5
//   4 -3-11 -2 10
//  -2 -7-12  7
//  -6 -8 -3
//  -7  1
//   2

//       SAT 21432843 BTS 21433281 NDS 21432920 FAILS     3680 BTS/s     2079584 PPGS/s   12108505213 PPGS             0 CKS  5822.6 s  5822.6 s

