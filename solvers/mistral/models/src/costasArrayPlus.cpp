
#include <mistral_sol.h>

/* 
 * a model of the costas array problem
 */

using namespace std;


void print(VarArray& X)
{
  /// output
  int i, j, N=X.size();
  cout << endl << endl;    
  for(i=0; i<N; ++i)
    cout << setw(3) << X[i].value() ;
  cout << endl ;
  for(i=1; i<N; ++i) {
    for(j=0; j<N-i; ++j)
      cout << setw(3) << (X[j].value() - X[j+i].value());
    cout << endl;
  }
}

int main(int argc, char *argv[])
{  
  /// input
  int i, j, N = ( argc > 1 ? atoi(argv[1]) : 11 );
  int Vertical = (argc > 2 ? atoi(argv[2]) : 1);
  int Triangle = (argc > 3 ? atoi(argv[3]) : 1);
  int Symmetry = (argc > 4 ? atoi(argv[4]) : 1);
  int Card     = (argc > 5 ? atoi(argv[5]) : 1);
  int AC       = (argc > 6 ? atoi(argv[6]) : 1);

  /// vars
  CSP model;
  VarArray horizontal[N];
  VarArray opposite[N];
  VarArray vertical[N];
  VarArray absValues;


  for(j=0; j<N; ++j) {
    Variable x(1,N);
    horizontal[0].add(x);
  }
  for(i=1; i<N; ++i) {
    for(j=0; j<N-i; ++j) {
      horizontal[i].add(horizontal[0][j] - horizontal[0][j+i]);
      model.add( horizontal[i][j] != 0 );
    }
    model.add( AllDifferent(horizontal[i-1]) );    
  }

  if( Vertical || Triangle) 
    for(i=1; i<N; ++i) 
      for(j=0; j<N-i; ++j) {
	opposite[i].add( -(horizontal[i][j]) );
	model.add( opposite[i][j] != 0 );
      }

  if(Vertical) {
    // vertical alldiffs
    for(j=0; j<N; ++j) {
      for(i=0; i<N-j; ++i) 
	vertical[j].add(horizontal[i][j]);
      for(i=1; i<=j; ++i) 
	vertical[j].add(opposite[i][j-i]);
      model.add( AllDifferent(vertical[j]) );
    }
  } 
  if(Triangle) {
    //triangles
    int x, y, z;
    for(x=1; x<N-1; ++x) 
      for(y=x+1; y<N; ++y) 
	for(z=y+1; z<=N; ++z) 
	  model.add( (horizontal[y-x][x-1] + horizontal[z-y][y-1] + opposite[z-x][x-1]) == 0 );
  } 
  if(Symmetry) {
    //symmetry breaking
    if(N%2) {
      if(Symmetry < 0)
	model.add( horizontal[0][N/2] >= (N/2+1) );
      else
	model.add( horizontal[0][N/2] <= (N/2+1) );
    } else {
      if(Symmetry < 0)
	model.add( horizontal[0][N/2] > N/2 );
      else
	model.add( horizontal[0][N/2] <= N/2 );
    }
    if(Symmetry < 0)
      model.add( horizontal[0][N/2] > horizontal[0][N/2-1] );
    else
      model.add( horizontal[0][N/2] < horizontal[0][N/2-1] );
  }


  if(Card) {
    for(i=1; i<N; ++i) 
      for(j=0; j<N-i; ++j) 
	absValues.add( Abs(horizontal[i][j]) );

    int occurrence[N];
    for(i=1; i<N; ++i)
      occurrence[i] = N-i;
    occurrence[0] = 0;
    
    model.add( Gcc(absValues, 0, N-1, occurrence, occurrence) );
  }

  /// search
  if(argc < 8) {
    ImpactOverWDeg heuristic(2);
    //DomOverDeg heuristic(2);
    Solver s;
    if(AC)
      s.boundOptimised = false;
    s.init(model, heuristic);
    s.setVerbosity(1);

//     s.print(cout);
//     cout << endl;

    s.solve_and_restart(GEOMETRIC, 64, 1.33);
    print(horizontal[0]);
    s.printStatistics(std::cout);
  } else {
    Lexicographic heuristic;
    Solver s;
    if(AC) {
      s.boundOptimised = false;
    } else {
      s.boundOptimised = true;
    }


    VarArray searchVars;
    
    if(N%2) {// odd
      searchVars.add( horizontal[0][N/2] );
      for(i=1; i<=N/2; ++i) {
	searchVars.add( horizontal[0][N/2-i] );
	searchVars.add( horizontal[0][N/2+i] );
      }
    } else {
      for(i=0; i<N/2; ++i) {
	searchVars.add( horizontal[0][N/2+i] );
	searchVars.add( horizontal[0][N/2-i-1] );
      }
    }
    

    s.init(model, searchVars, heuristic);
    if( argc == 9 )
      s.FIND_ALL = -1;
    s.setVerbosity(4);


//     s.print( cout );
//     cout << endl;

    s.solve();

    if( argc < 9 )
      print(horizontal[0]);

    s.printStatistics(std::cout);
   }

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


//   7  5  8 12 13  1  9  3 18 11 16  6 15 10  2  4 17 14
//   2 -3 -4 -1 12 -8  6-15  7 -5 10 -9  5  8 -2-13  3
//  -1 -7 -5 11  4 -2 -9 -8  2  5  1 -4 13  6-15-10
//  -5 -8  7  3 10-17 -2-13 12 -4  6  4 11 -7-12
//  -6  4 -1  9 -5-10 -7 -3  3  1 14  2 -2 -4
//   6 -4  5 -6  2-15  3-12  8  9 12-11  1
//  -2  2-10  1 -3 -5 -6 -7 16  7 -1 -8
//   4-13 -3 -4  7-14 -1  1 14 -6  2
// -11 -6 -8  6 -2 -9  7 -1  1 -3
//  -4-11  2 -3  3 -1  5-14  4
//  -9 -1 -7  2 11 -3 -8-11
//   1-10 -2 10  9-16 -5
//  -8 -5  6  8 -4-13
//  -3  3  4 -5 -1
//   5  1 -9 -2
//   3-12 -6
// -10 -9
//  -7
//       SAT  1302091 BTS  1302387 NDS  1302142 FAILS     3546 BTS/s     3072661 PPGS/s    1128035615 PPGS             0 CKS  367.12 s  367.12 s



 
