
#include <mistral_sol.h>

/*
 * a model of the magic square problem
 */


using namespace std;

void print_square(VarArray& square, const int N) {
  for(int i=0; i<N; ++i) {
    std::cout<<"\t+";
    for(int j=0; j<N; ++j)
      std::cout<<"----+";
    std::cout<<std::endl<<"\t";
    for(int j=0; j<N; ++j)
      std::cout<<"|"<<std::setw(3)<<square[i*N+j].value()<<" ";
    std::cout<<"|"<<std::endl;
  }
  cout<<"\t+";
  for(int j=0; j<N; ++j)
    cout<<"----+";
  cout<<endl<<endl;
}

int main(int argc, char *argv[])
{
  int N = ( argc > 1 ? atoi(argv[1]) : 6 );
  int sum = N*(N*N-1)/2;

  std::cout << endl << "\tMagic Square (order "
       << N << ")" << endl << endl;

  // create the csp
  CSP model;
  VarArray square( N*N, 0, N*N-1 );
  VarArray alignment(N);

  // Constraints
  model.add( AllDifferent(square) );

  for(int i=0; i<N; ++i) {

    // Sums on rows
    for(int j=0; j<N; ++j)
      alignment[j] = square[j+N*i];
    model.add( Sum( alignment ) == sum );

    // Sums on columns
    for(int j=0; j<N; ++j)
      alignment[j] = square[j*N + i];
    model.add( Sum( alignment ) == sum );
  }

  // Sum on the first diagonal
  for(int i=0; i<N; ++i)
    alignment[i] = square[i*(N+1)];
  model.add( Sum( alignment ) == sum );

  // Sum on the second diagonal
  for(int i=0; i<N; ++i)
    alignment[i] = square[(i+1)*(N-1)];
  model.add( Sum( alignment ) == sum );


  model.add( square[0] > square[N-1] );
  model.add( square[0] > square[N*(N-1)] );


//  MinDomain heuristic(2);
  Lexicographic heuristic;
  Solver s(model, heuristic);

//  s.print(cout);
//  cout << endl;

  s.setVerbosity( 2 );

//  	s.solve();
//  	s.pseudoRldsSolve();
//  	s.pseudoRldsSolve_level();
//  	s.pseudoRldsSolve_variable();
//  	s.ldsSolve();
//  	s.ldsStackSolve();
//  	s.ldsDeltaSolve();
//  	s.ldsSolve_level();
//  	s.ldsSolve_variable();
//  	s.ddsSolve();
//  	s.ddsSolve_level();
  	s.ddsSolve_variable();

  print_square(square, N);

}
