
#include <mistral_sol.h>

using namespace Mistral;
using namespace std;

/*****************************************
 * Main
 *****************************************/
int main(int argc, char *argv[]) { 
 
  // command line
  int i, j, N = ( argc > 1 ? 3*atoi(argv[1]) : 9 );

  CSP model;
  VarArray QueenX(4*N/3, 0, N-1);
  VarArray QueenY(4*N/3, 0, N-1);

  for(i=0; i<N/3; ++i) {
    // the pairs on columns
    model.add( QueenX[i] == QueenX[i+N/3] );

    // the pairs on rows
    model.add( QueenY[i+2*N/3] == QueenY[i+N] );
  }

  VarArray column;
  VarArray row;
  for(i=0; i<N; ++i) {
    column.add( QueenX[N/3+i] );
    row.add( QueenY[i] );
  }

  model.add( AllDifferent(column) );
  model.add( AllDifferent(row) );


  for(i=0; i<4*N/3; ++i)
    for(j=i+1; j<4*N/3; ++j) {
      if( ( i < N/3 && j == i+N/3 ) ||
	  ( j >= N && i == j-N/3 ) ) continue;
      // make sure that queen[i] is not on the same diagonal as queen[j]
      model.add( (Abs( QueenX[i] - QueenX[j] ) != Abs( QueenY[i] - QueenY[j] )) );
    }

  //symmetry breaking
  for(i=1; i<N/3; ++i) {
    model.add( QueenX[i-1] < QueenX[i] );    
    model.add( QueenY[2*N/3+i-1] < QueenY[2*N/3+i] );
  }
  for(i=0; i<N/3; ++i) {
    model.add( QueenY[i] < QueenY[i+N/3] );
    model.add( QueenX[i+2*N/3] < QueenX[i+N] );
  }
  //model.add( QueenY[0] <= QueenX[2*N/3] );
  //model.add( QueenY[0] < N/2 );
  //model.add( QueenX[2*N/3] < N/2 );

  // < < 16735
  // < > 22417
  // > < 10580
  // > >  9203
  // no:  1163

  DomOverWDeg heuristic;
  //MinDomMin heuristic;
  //Impact heuristic;
  //MinDomain heuristic;
  //Lexicographic heuristic;

  Solver s( model, heuristic );

  //s.setVerbosity(1);
  //s.solve_and_restart(Solver::GEOMETRIC, 32, 1.2);
  //s.solve();


  s.print( cout );
  cout << endl;

  int chessboard[N*N];
  

  s.startNewSearch();

  while( s.getNextSolution() ) {
    std::fill( chessboard, chessboard+N*N, 0 );

    for(i=0; i<4*N/3; ++i) { 
      cout << (QueenX[i].value()) << " " << (QueenY[i].value()) << endl;  
      chessboard[ (QueenX[i].value())*N + QueenY[i].value() ] = 1;
    }
    
    for(i=0; i<N; ++i) {
      for(j=0; j<N; ++j)
	cout << chessboard[i*N+j] << " ";
      cout << endl;
    }
    cout << endl;

    //break;
  }

  cout << endl;
  s.printStatistics( cout );
  cout << endl;
}







// 0 1
// 0 2

// 1 4
// 1 5

// 2 7
// 2 8

// 6 6
// 8 3
// 4 0
// 7 6
// 5 3
// 3 0
// 0 1 1 0 0 0 0 0 0 
// 0 0 0 0 1 1 0 0 0 
// 0 0 0 0 0 0 0 1 1 
// 1 0 0 0 0 0 0 0 0 
// 1 0 0 0 0 0 0 0 0 
// 0 0 0 1 0 0 0 0 0 
// 0 0 0 0 0 0 1 0 0 
// 0 0 0 0 0 0 1 0 0 
// 0 0 0 1 0 0 0 0 0 

// 0 0 1 0 1 0 0 0 0 0 0 0 
// 0 0 0 0 0 0 1 0 1 0 0 0 
// 0 0 0 1 0 0 0 0 0 0 0 1 
// 1 0 0 0 0 0 0 0 0 0 0 0 
// 0 0 0 0 0 0 0 0 0 0 1 0 
// 0 1 0 0 0 0 0 0 0 0 0 0 
// 0 0 0 0 0 1 0 0 0 0 0 0 
// 0 0 0 0 0 0 0 0 0 0 1 0 
// 1 0 0 0 0 0 0 0 0 0 0 0 
// 0 0 0 0 0 0 0 1 0 1 0 0 
// 0 0 0 0 0 1 0 0 0 0 0 0 
// 0 1 0 0 0 0 0 0 0 0 0 0 

// 0 0 0 0 0 0 0 0 1 0 0 0 1 0 0 
// 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 
// 0 0 0 0 0 0 0 1 0 0 0 1 0 0 0 
// 0 1 0 0 0 0 0 0 0 0 0 0 0 1 0 
// 0 0 0 1 0 0 1 0 0 0 0 0 0 0 0 
// 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 
// 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
// 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 
// 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 
// 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 
// 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 
// 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
// 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 
// 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 
// 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 

// 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 
// 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 
// 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 
// 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
// 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
// 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 
// 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 
// 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
// 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 1 0 0 
// 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 
// 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
// 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
// 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
// 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 
// 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
// 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
// 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 
// 0 0 0 0 0 0 0 1 0 0 0 1 0 0 0 0 0 0 
