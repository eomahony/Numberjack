
#include <mistral_sol.h>

using namespace Mistral;
using namespace std;

/*****************************************
 * A model of the NQueens problem
 *****************************************/


int binlog2( unsigned long int x ) {
  if( !x ) return NOVAL;
  if( x == 1 ) return 0;
  int exponent = 0;
  while( (x >> (++exponent)) > 1 );
  return exponent;
}
unsigned long int luby_seq2( unsigned long int iteration ) {
  int thelog = binlog2( iteration );
  if( (int)iteration == (1 << (thelog + 1))-1 )
    return (1 << thelog);
  return luby_seq2( iteration - (1 << thelog) + 1 );
}




/*****************************************
 * Print a solution
 *****************************************/
void printChessboard( VarArray x )
{
  for(int q = 0; q < x.size(); ++q) 
    {
      cout<<"\n\t";
      for(int i = 0; i < x.size(); ++i) 
	{
	  if(x[q].value() == i)
	    cout<<"Q";
	  else
	    cout<<".";
	}   
    }
  cout<<"\n";
}

/*****************************************
 * Model
 *****************************************/
void initChessBoard1( CSP& model, VarArray& queens ) {  
  
  // extra variables  
  VarArray first_diagonal( queens.size() );
  VarArray second_diagonal( queens.size() );
  for(int i=0; i<queens.size(); ++i) {
    first_diagonal[i]  = (queens[i] + i);
    second_diagonal[i] = (queens[i] - i);
  }
  
  // the constraints
  model.add( AllDifferent(first_diagonal) );
  model.add( AllDifferent(second_diagonal) );
  model.add( AllDifferent(queens) );
}

/*****************************************
 * Model
 *****************************************/
void initChessBoard2( CSP& model, VarArray& queens ) {  
 
  // the constraints
  for(int i=0; i<queens.size(); ++i) 
    for(int j=i+1; j<queens.size(); ++j) {
      model.add( queens[i] != queens[j] );
      model.add( (queens[i] + i) != (queens[j] + j) );
      model.add( (queens[i] - i) != (queens[j] - j) );
    }
}


// /*****************************************
//  * Model
//  *****************************************/
// void initChessBoard3( CSP& model, VarArray& queens, int t ) {  
  
//   // extra variables  
//   VarArray first_diagonal( queens.size() );
//   VarArray second_diagonal( queens.size() );
//   for(int i=0; i<queens.size(); ++i) {
//     first_diagonal[i]  = (queens[i] + i);
//     second_diagonal[i] = (queens[i] - i);
//   }
  
//   // the constraints
//   model.add( LPAllDifferent(queens, t) );
//   model.add( LPAllDifferent(first_diagonal, t) );
//   model.add( LPAllDifferent(second_diagonal, t) );
// }

/*****************************************
 * Main
 *****************************************/
int main(int argc, char *argv[]) { 
 
  for(int i=0; i<100; ++i)
    std::cout << luby_seq2(i) << " ";
  std::cout << std::endl;
  
  exit(1);


  // command line
  int N = ( argc > 1 ? atoi(argv[1]) : 20 );
  int E = ( argc > 2 ? atoi(argv[2]) : 20 );
  int M = ( argc > 2 ? atoi(argv[2]) : 1 );

  while( N <= E ) {

    // the model
    VarArray queens(N,N);
    CSP model;

    switch( M ) {
    case 1 : { initChessBoard1( model, queens ); break; }
    case 2 : { initChessBoard2( model, queens ); break; }
    default : { initChessBoard1( model, queens ); break; }
    }

    // the solver
    MinDomMin heuristic;
    //Solver s( model, queens, heuristic );
    Solver s( model );
    s.add( heuristic );
    

//      s.print( cout );
//      cout<< endl;

//     s.startNewSearch();
//     while( s.getNextSolution() != Solver::UNSAT )
//       {
// 	printChessboard( queens );
// 	cout << endl;
//       }

    //s.print( cout );
    //cout << endl;


     //s.setVerbosity(3);

    //s.solve_and_restart(GEOMETRIC, 32, 2);
    s.solve();

    // the solution
    if( N <= 50 )
      printChessboard( queens );
    
    // the statistics
    //printChessboard( queens );
    
    cout << N << " queens: \t" ;
    s.printStatistics(cout);
    cout << endl ;


    ++N;
    }
}


