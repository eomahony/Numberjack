
#include <mistral_solver.hpp>
#include <mistral_variable.hpp>
#include <mistral_search.hpp>


using namespace std;
using namespace Mistral;


int main(int argc, char *argv[])
{  

  // double time, start = get_run_time();
  // for(int i=0; i<10000000; ++i) {
  //   time = get_run_time();
  // }
  // std::cout << (get_run_time() - start) << std::endl;


  // // double time, start = cpu_time();
  // // for(int i=0; i<10000000; ++i) {
  // //   time = cpu_time();
  // // }
  // // std::cout << (cpu_time() - start) << std::endl;

  // exit(1);


  int i, j, N=8;
  if(argc>1) N = atoi(argv[1]);

  VarArray X(N, 1, N);
  Vector< VarArray > differences;
  Solver s;

  s.add( AllDiff(X) );



  VarArray scope;
  for(i=1; i<N-1; ++i) {
    scope.clear();
    for(j=0; j<N-i; ++j) {
      scope.add(X[j]-X[j+i]);
    }
    s.add( AllDiff(scope, FORWARD_CHECKING) );
    differences.add(scope);
  }


  //cout << s << endl;

  s.rewrite();

  //cout << s << endl;

  s.consolidate();
  
  cout << s << endl;

  s.parameters.verbosity = 0;

#ifdef _MONITOR
  s.monitor_list << X << "\n";
  for(i=0; i<N-2; ++i) {
    s.monitor_list << differences[i] << "\n";
  }
#endif

  ImpactManager *im = new ImpactManager(&s);


  
  /*
    if(s.depth_first_search(X, 
			  new GenericHeuristic< 
			  // GenericNeighborDVO< 
			  //   FailureCountManager, 
			  //   //PruningCountManager, 
			  
			  //   SelfPlusAverage,
			  
			  //   MinDomainOverWeight 
			  //   //MinNeighborDomainOverNeighborWeight
			  //   , 1
			  
			  GenericDVO< 
			  MinDomainOverWeight, 1,
			  FailureCountManager
			  >,
			  MinValue >(&s), 
			  new Geometric()) == SAT) {
			  for(i=0; i<N; ++i)
			  cout << setw(3) << X[i].get_solution_int_value() << " " ;
    cout << endl;
    for(i=0; i<N-2; ++i) {
    for(j=0; j<N-i-1; ++j)
	cout << setw(3) << differences[i][j].get_solution_int_value() << " " ;
      cout << endl;
    }
    cout << setw(3) << (X[0].get_solution_int_value() - X[N-1].get_solution_int_value()) << endl << endl;
  }
  //std::cout << s.statistics << std::endl;
  */



  s.initialise_search(X, 
		      new GenericHeuristic< 
			GenericDVO< 
			  MinDomainOverWeight, 1,
			  FailureCountManager
			  >,
			MinValue >(&s));

  while(s.get_next_solution() == SAT) {
    for(i=0; i<N; ++i)
      cout << setw(3) << X[i].get_solution_int_value() << " " ;
    cout << endl;
    for(i=0; i<N-2; ++i) {
      for(j=0; j<N-i-1; ++j)
	cout << setw(3) << differences[i][j].get_solution_int_value() << " " ;
      cout << endl;
    }
    cout << setw(3) << (X[0].get_solution_int_value() - X[N-1].get_solution_int_value()) << endl << endl;
  }
  //std::cout << s.statistics << std::endl;


}


