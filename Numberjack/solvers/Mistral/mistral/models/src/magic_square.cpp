
#include <mistral_sol.h>

/* 
 * a model of the magic square problem
 */


using namespace std;

int main(int argc, char *argv[])
{  
  int N = ( argc > 1 ? atoi(argv[1]) : 6 );
  int sum = N*(N*N-1)/2;

  cout << endl << "\tMagic Square (order " 
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


  //MinDomain heuristic(2);
  //Solver s(model, heuristic);
  Solver s(model);
  s.setHeuristic("MinDomain", "RandomSplit", 2);


  s.setVerbosity( 1 );
  //s.setDomainSplitting();
  //s.setRandomized();
  s.solve_and_restart(GEOMETRIC, 32, 1.1);
  //s.solve_and_restart(LUBY, 256);
  //s.solve();
  
  for(int i=0; i<N; ++i) {
    std::cout<<"\t+";
    for(int j=0; j<N; ++j)
      std::cout<<"----+";
    std::cout<<std::endl<<"\t";
    for(int j=0; j<N; ++j)
      std::cout<<"|"<<std::setw(3)<<square[i*N+j].value()+1<<" ";
    std::cout<<"|"<<std::endl;
  }
  cout<<"\t+";
  for(int j=0; j<N; ++j)
    cout<<"----+";
  cout<<endl<<endl;
  
  s.printStatistics(std::cout);
  cout<<endl<<endl;
}




