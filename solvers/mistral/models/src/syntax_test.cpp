
#include <mistral_sol.h>

using namespace std;

void setupMix( CSP& model, VarArray& X )
{
  model.add( X[0] >= -10 && X[0] <= 10 );
  model.add( X[1] >= -10 && X[1] <= 10 );
  model.add( X[2] >= -10 && X[2] <= 10 );
  model.add( X[3] >= -10 && X[3] <= 10 );
  model.add( X[4] >= -10 && X[4] <= 10 );
  model.add( X[5] >= -10 && X[5] <= 10 );

  VarArray S;
  S.add( X[0] );
  S.add( X[1] );
  S.add( X[2] );


  model.add( -X[3] + X[1] + Sum(S) == X[4] - X[5] );
  model.add( X[0] != X[1] );
  model.add( X[4] != X[5] );
  model.add( X[0] > X[4] + X[5] );
  model.add( X[3] <= X[0] + X[1] + X[2] );


  model.add( (X[0] == 3) == (X[1]+2 < 4) );

  model.add( X[X[X[4]]] > X[0] );

  //  model.add( X[X[X[X[4]]]] > X[0] );
  //  model.add( Maximise( X[X[5]] + (X[0] - X[1]) ) );
}


// void setupMixLo( CSP& model, BuildObject **X )
// {
//   model.add( CSP::_And( CSP::_Precedence(-10, X[0]),
// 			CSP::_Precedence(X[0], 10) ) );
//   model.add( CSP::_And( CSP::_Precedence(-10, X[1]),
// 			CSP::_Precedence(X[1], 10) ) );
//   model.add( CSP::_And( CSP::_Precedence(-10, X[2]),
// 			CSP::_Precedence(X[2], 10) ) );
//   model.add( CSP::_And( CSP::_Precedence(-10, X[3]),
// 			CSP::_Precedence(X[3], 10) ) );
//   model.add( CSP::_And( CSP::_Precedence(-10, X[4]),
// 			CSP::_Precedence(X[4], 10) ) );
//   model.add( CSP::_And( CSP::_Precedence(-10, X[5]),
// 			CSP::_Precedence(X[5], 10) ) );

//   BuildObject **scope = new BuildObject*[4];
//   scope[0] = X[0];
//   scope[1] = X[1];
//   scope[2] = X[2];

//   model.add( CSP::_Equal( CSP::_Add( CSP::_Add( CSP::_Negation( X[3] ), X[1] ), CSP::_Sum(scope, 3) ),
// 			  CSP::_Sub( X[4], X[5] )
// 			  ) ) ;
//   model.add( CSP::_NotEqual( X[0], X[1] ) );
//   model.add( CSP::_NotEqual( X[4], X[5] ) );
//   model.add( CSP::_Precedence( CSP::_Add( X[4], X[5] ), 1, X[0] ) );
//   model.add( CSP::_Precedence( X[3], 0, CSP::_Add( X[0], CSP::_Add( X[1], X[2] ) ) ) );

// //   BuildObject **Y = new BuildObject*[7];
// //   Y[0] = X[0];
// //   Y[1] = X[1];
// //   Y[2] = X[2];
// //   Y[3] = X[3];
// //   Y[4] = X[4];
// //   Y[5] = X[5];

//   model.add( CSP::_Precedence( X[0], 1, CSP::_Element(X, 6, CSP::_Element(X, 6, X[4]) ) ) );


//   //  model.add( CSP::_Precedence( X[0], 1,
//   //			       CSP::_Element(X, 6, CSP::_Element(X, 6, CSP::_Element(X, 6, X[4]) ) ) ) );

//   //  model.add( CSP::_Maximise( CSP::_Add( CSP::_Element(X, 6, X[5]), 
//   //					CSP::_Sub(X[0], X[1]) ) ) );
// }


// sum equal to some constant
void setupSum1( CSP& model, VarArray& X )
{
  model.add( X[0] >= -10 && X[0] <= 10 );
  model.add( X[1] >= -10 && X[1] <= 10 );
  model.add( X[2] >= -10 && X[2] <= 10 );
  model.add( X[3] >= -10 && X[3] <= 10 );
  model.add( X[4] >= -10 && X[4] <= 10 );
  model.add( X[5] >= -10 && X[5] <= 10 );
  
  model.add( X[1] == 1 );

  model.add( X[2] != -4 );

  model.add( X[3] == 3 );

  model.add( X[5] == 5 );

  model.add( Sum( X ) == 25 );

}

// sum as variable
void setupSum2( CSP& model, VarArray& X )
{
  model.add( X[0] >= -10 && X[0] <= 10 );
  model.add( X[1] >= -10 && X[1] <= 10 );
  model.add( X[2] >= -10 && X[2] <= 10 );
  model.add( X[3] >= -10 && X[3] <= 10 );
  model.add( X[4] >= -10 && X[4] <= 10 );
  model.add( X[5] >= -10 && X[5] <= 10 );

  VarArray Y;
  Y.add(X[0]);
  Y.add(X[1]);
  VarArray Z;
  Z.add(X[2]);
  Z.add(X[3]);
  Z.add(X[4]);
  Z.add(X[5]);
  
  model.add( Sum( X ) == 30 );

  model.add( Sum( Y ) == Sum( Z ) );
}

// offset, factor, negation
void setupOFN( CSP& model, VarArray& X )
{
  model.add( X[0] >= -100 && X[0] <= 100 );
  model.add( X[1] >= -100 && X[1] <= 100 );
  model.add( X[2] >= -100 && X[2] <= 100 );
  model.add( X[3] >= -100 && X[3] <= 100 );
  model.add( X[4] >= -100 && X[4] <= 100 );
  model.add( X[5] >= -100 && X[5] <= 100 );

  model.add( X[0] == (X[1] + 3) );

  model.add( X[2] == (X[3] * 3) );

  model.add( X[3] == -X[4] );

  model.add( X[4] == -X[5] );

}

// abs value
void setupAbs( CSP& model, VarArray& X )
{
  model.add( X[0] >= -5 && X[0] <= 5 );
  model.add( X[1] >= -5 && X[1] <= 5 );
  model.add( X[2] >= -5 && X[2] <= 5 );
  model.add( X[3] >= -5 && X[3] <= 5 );
  model.add( X[4] >= -5 && X[4] <= 5 );
  model.add( X[5] >= -5 && X[5] <= 5 );

  model.add( Sum( X ) == 15 );

  model.add( Abs( (X[1] - X[0]) ) == Abs( (X[2] - X[3]) ) );

  model.add( Abs( (X[1] - X[0]) ) > Abs( (X[4] - X[5]) ) );
}

// and or not
void setupLog( CSP& model, VarArray& X )
{
  model.add( X[0] >= -10 && X[0] <= 10 );
  model.add( X[1] >= -10 && X[1] <= 10 );
  model.add( X[2] >= -10 && X[2] <= 10 );
  model.add( X[3] >= -10 && X[3] <= 10 );
  model.add( X[4] >= -10 && X[4] <= 10 );
  model.add( X[5] >= -10 && X[5] <= 10 );

  model.add( ( X[0] < X[1] ) || ( X[1] < X[0] ) );

  model.add( !(( X[2] == X[3] ) && ( X[4] == X[5] )) );

}

// lower/upper bounds
void setupBound( CSP& model, VarArray& X )
{
  model.add( X[0] >= -10 && X[0] <= 10 );
  model.add( X[1] >= -10 && X[1] <= 10 );
  model.add( X[2] >= -10 && X[2] <= 10 );
  model.add( X[3] >= -10 && X[3] <= 10 );
  model.add( X[4] >= -10 && X[4] <= 10 );
  model.add( X[5] >= -10 && X[5] <= 10 );

  model.add( ( X[0] > 3 ) || ( X[1] >= 5 ) || ( X[2] > 7 ) );

  model.add( ( X[3] > 3 ) || ( X[4] >= 5 ) || ( X[5] > 7 ) );

  model.add( X[2] < 3 || X[5] < 3 );

}


// Ranges
void setupRanges( CSP& model, VarArray& X )
{
  model.add( Sum(X) == 1000 );

  model.add( X[0] == X[1] || X[0] == X[5] );

  model.add( X[2] == X[3] );
}


// Ranges
void setupOffset( CSP& model, VarArray& X )
{

  model.add( X[0] >= -3 && X[0] <= 3 );
  model.add( X[1] >= -3 && X[1] <= 3 );
  model.add( X[2] >= -3 && X[2] <= 3 );
  model.add( X[3] >= -3 && X[3] <= 3 );
  model.add( X[4] >= -3 && X[4] <= 3 );
  model.add( X[5] >= -3 && X[5] <= 3 );

  model.add( ( X[0] + 1 ) == Abs( X[1] - X[2] ) );

  model.add( AllDifferent( X ) );
}

void setupQueens( CSP& model, VarArray& X )
{

  model.add( X[0] >= -10 && X[0] <= 10 );
  model.add( X[1] >= -10 && X[1] <= 10 );
  model.add( X[2] >= -10 && X[2] <= 10 );
  model.add( X[3] >= -10 && X[3] <= 10 );
  model.add( X[4] >= -10 && X[4] <= 10 );
  model.add( X[5] >= -10 && X[5] <= 10 );

  for(int i=1; i<X.size(); ++i)
    model.add( ( X[i-1] != X[i] ) && ( Abs( X[i-1] - X[i] ) != i ) );

  model.add( (X[5] + 0) == 1 );
}


void setupPreced( CSP& model, VarArray& X )
{

  model.add( X[0] >= -10 && X[0] <= 10 );
  model.add( X[1] >= -10 && X[1] <= 10 );
  model.add( X[2] >= -10 && X[2] <= 10 );
  model.add( X[3] >= -10 && X[3] <= 10 );
  model.add( X[4] >= -10 && X[4] <= 10 );
  model.add( X[5] >= -10 && X[5] <= 10 );

  model.add( X[0] < X[1] );
  model.add( X[1] < X[2] );
  model.add( X[2] == 3 );
  model.add( !(X[4] >= X[5]) );
  model.add( !(X[2] >= X[4]) );
}

// offset, factor, negation
void setupOFN2( CSP& model, VarArray& X )
{
  model.add( X[0] >= -100 && X[0] <= 100 );
  model.add( X[1] >= -100 && X[1] <= 100 );
  model.add( X[2] >= -100 && X[2] <= 100 );
  model.add( X[3] >= -100 && X[3] <= 100 );
  model.add( X[4] >= -100 && X[4] <= 100 );
  model.add( X[5] >= -100 && X[5] <= 100 );

  VarArray Y;
  Y.add(X[0]);
  Y.add(X[1]);
  VarArray Z;
  Z.add(X[2]);
  Z.add(X[3]);
  Z.add(X[4]);
  Z.add(X[5]);
  
  model.add( Sum( X ) == 30 );

  model.add( Sum( Y ) == Sum( Z ) );
  
  model.add( X[0] == (X[1] + 3) );

}


// offset, factor, negation
void setupSum3( CSP& model, VarArray& X )
{
  model.add( X[0] >= -100 && X[0] <= 100 );
  model.add( X[1] >= -100 && X[1] <= 100 );
  model.add( X[2] >= -100 && X[2] <= 100 );
  model.add( X[3] >= -100 && X[3] <= 100 );
  model.add( X[4] >= -100 && X[4] <= 100 );
  model.add( X[5] >= -100 && X[5] <= 100 );

  VarArray Y;
  Y.add(X[0]);
  Y.add(X[1]);
  VarArray Z;
  Z.add(X[2]);
  Z.add(X[3]);
  Z.add(X[4]);
  Z.add(X[5]);
  Z.add(X[5]);

  int wx[7] = {1, -2, 3, -4, 1, 3, 0};

  model.add( Sum( X, wx ) == 30 );

  int wz[6] = {3, -2, 0, -2, 2, 0};

  model.add( Sum( Y ) == Sum( Z, wz ) );

}



void model_and_solve( const int M ) 
{
  int n=6;
  CSP model;

  VarArray X(n, -1000, 1000);

  if( M < 0 ) 
    {
  
    }
  else
    {    

      switch(M) {
      case 0 :  setupMix( model, X );    break;
      case 1 :  setupSum1( model, X );   break;
      case 2 :  setupSum2( model, X );   break;
      case 3 :  setupOFN( model, X );    break;
      case 4 :  setupAbs( model, X );    break;
      case 5 :  setupLog( model, X );    break;
      case 6 :  setupBound( model, X );  break;
      case 7 :  setupRanges( model, X ); break;
      case 8 :  setupOffset( model, X ); break;
      case 9 :  setupQueens( model, X ); break;
      case 10:  setupPreced( model, X ); break;
      case 11:  setupOFN2( model, X );    break;
      case 12:  setupSum3( model, X );    break;
      }
    }

  MinDomain heuristic;
  Solver s(model, heuristic);  
 
  cout << endl << "MODEL_" << M << endl;
  s.print( cout );
  cout << endl;
  
  s.solve();

  for(int i=0; i<6; ++i)
    cout << " " << X[i].value() ;
  cout << endl;
  
  s.printStatistics(cout);
  cout << endl;
}


int main(int argc, char *argv[])
{  
  int i, M = (argc > 1 ? atoi(argv[1]) : 11);

  if( argc > 1 )
    model_and_solve( M );
  else {
    for( i=0; i<M; ++i )
      model_and_solve( i );
  }

}






