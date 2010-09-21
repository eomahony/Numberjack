#include <mistral_sol.h>
#include <mistral_glo.h>

/*  Golomb ruler */

using namespace std;
using namespace Mistral;

// Known optimal rulers. Rulers of lower order are used as lower bound for distances.
static const int ruler[] = {0,1,3,6,11,17,25,34,44,55,72,85,106,127};


int main(int argc, char *argv[])
{  
  int nbMarks = (argc > 1 ? atoi(argv[1]) : 8);
  int rulerSize = ( 2 << (nbMarks-1) ); // naive upper bound for the ruler size
  int filtering = (argc > 2 ? atoi(argv[2]) : 2);
  int K = (argc > 3 ? atoi(argv[3]) : 4);

  // variables
  CSP model;
  VarArray mark(nbMarks, 0, rulerSize-1);  
  VarArray distance(nbMarks*(nbMarks-1)/2, 1, rulerSize-1);

  // constraints
  int i, j, k=0;
  for(i=1; i<nbMarks; ++i) {
    
    // strict ordering of the marks
    model.add( mark[i-1] < mark[i] );

    // lower bounds from earlier golomb rulers
    //model.add( mark[i] >= ruler[i] );

    for(j=0; j<i; ++j) {

      // lower bounds from earlier golomb rulers
      //model.add( distance[k] >= ruler[i-j] );

      // set up the distances
      model.add( mark[i] == (mark[j] + distance[k++]) );
    }
  }

  // distances are all different
  model.add( AllDifferent(distance, filtering, K) );
 
   // symmetry breaking
   model.add( mark[0] == 0 );
   if( nbMarks > 2 )
     model.add( distance[0] < distance[k-1] );

  // set up objective function
  model.add( Minimise( mark[nbMarks-1] ) );

  NoOrder heuristic;
  Solver s( model, mark, heuristic );
  s.setVerbosity(1);

  s.solve();

  s.printStatistics( cout, (NDS | RUNTIME | OPTTIME | PROOFTIME | SPEEDB) );
  cout << endl;
}


