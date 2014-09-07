

#include <mistral_sol.h>

/*
 * a template model
 */


using namespace std;

int main(int argc, char *argv[])
{
	int N = (argc > 1 ? atoi(argv[1]) : 8);

	cout << "pigeons_lds " << N << endl;

	// create the csp
	CSP model;
	VarArray pigeons(N, 1, N-1);

	for(int i=0; i<N; ++i)
		for(int j=i+1; j<N; ++j)
			model.add( pigeons[i] != pigeons[j] );

	Solver s( model );
	Lexicographic heuristic;
	s.add(heuristic);

	//  ValSelector **valSelectors;
	//  valSelectors = new ValSelector*[N];
	//  for(int i=0; i<N; ++i)
	//	  valSelectors[i] = new ValSelectorMax(((VariableInt *)&pigeons[i]));

//	s.solve();
//	s.pseudoRldsSolve();
//	s.pseudoRldsSolve_level();
//	s.pseudoRldsSolve_variable();
//	s.ldsSolve();
//	s.ldsStackSolve();
//	s.ldsDeltaSolve();
//	s.ldsSolve_level();
//	s.ldsSolve_variable();
//	s.ddsSolve();
//	s.ddsSolve_level();
	s.ddsSolve_variable();
}


