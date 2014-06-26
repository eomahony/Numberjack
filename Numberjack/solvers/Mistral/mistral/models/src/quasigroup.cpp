
#include <mistral_sol.h>

/* 
 * a template model
 */

 
using namespace std;

//static const int NO=0;
static const int QG1=1;
static const int QG2=2;
static const int QG3=3;
static const int QG4=4;
static const int QG5=5;
static const int QG6=6;
static const int QG7=7;


void postQG3( CSP& model, VarArray& X, const int N )
{
  int i, j;
  for(i=0; i<N; ++i)
    for(j=0; j<N; ++j)
      model.add( X[X[i*N+j]*N+X[j*N+i]]==i );
}

void postQG4( CSP& model, VarArray& X, const int N )
{
  int i, j;
  for(i=0; i<N; ++i)
    for(j=0; j<N; ++j)
      model.add( X[X[j*N+i]*N+X[i*N+j]]==i );
}

void postQG5( CSP& model, VarArray *rows, const int N )
{
  int i, j;
  for(i=0; i<N; ++i)
    for(j=0; j<N; ++j)
      model.add( rows[i][rows[i][rows[j][i]]] == j );  
}

void postQG5b( CSP& model, VarArray& X, const int N )
{
  int i, j;
  for(i=0; i<N; ++i)
    for(j=0; j<N; ++j)
      model.add(X[X[X[j*N+i]*N+j]*N+j] == i);
}

void postQG6( CSP& model, VarArray& X, const int N )
{
  int i, j;
  for(i=0; i<N; ++i)
    for(j=0; j<N; ++j)
      model.add( X[X[i*N+j]*N+j] == X[X[i*N+j]+i*N] );
}


void postQG7( CSP& model, VarArray& X, const int N )
{
  int i, j;
  for(i=0; i<N; ++i)
    for(j=0; j<N; ++j)
      model.add( X[X[j*N+i]*N+j]==X[X[j*N+i]+i*N] );
}


int main(int argc, char *argv[])
{
  int i, j, N = (argc > 1 ? atoi(argv[1]) : 9);
  int type = (argc > 2 ? atoi(argv[2]) : 5);
  
  // Create the csp
  CSP model;
  VarArray X(N*N, 0, N-1);
  VarArray rows[N];
  VarArray cols[N];

  model.add( X );
  
  for(i=0; i<N; ++i) {
    VarArray row( N );
    VarArray col( N );
    for(j=0; j<N; ++j) {
      row[j] = X[i*N+j];
      col[j] = X[j*N+i];
    }
    model.add( AllDifferent(row) );
    model.add( AllDifferent(col) );
    rows[i] = row;
    cols[i] = col;
  }
 
  //  for (i=0; i<N; ++i) model.add( X[i*N+i] == i );   // Idempotency

  switch( type ) {
  case QG1 : cerr << "QG1 not handled" << endl; break;
  case QG2 : cerr << "QG2 not handled" << endl; break;
  case QG3 : postQG3(model, X, N); break;
  case QG4 : postQG4(model, X, N); break;
  case QG5 : postQG5(model, rows, N); break;
  case QG6 : postQG6(model, X, N); break;
  case QG7 : postQG7(model, X, N); break;
  default : ; 
  }

  Impact h;
  Solver s( model, h );

  if( s.solve() == SAT )
    {
      for(i=0; i<N; ++i) {
	for(j=0; j<N; ++j) 
	  cout << setw(3) << (rows[i][j].value());
	cout << endl;
      }
    }

  s.printStatistics( cout );
  cout << endl;
}


