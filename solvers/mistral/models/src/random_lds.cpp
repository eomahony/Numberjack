
#include <mistral_sol.h>
#include <mistral_gen.h>

using namespace std;
using namespace Mistral;


int timelimit, ins, Var, Dom, Con, Ari, Ngd, Ext, verbose;
long Seed;




const int nia = 10;
const char* int_ident[nia] = {"-s", "-i", "-t", "-v", "-d", "-c", "-a", "-n", "-e", "-verbose"};
int int_param[nia];

const int nsa = 1;
const char* str_ident[nsa] = {"-heu"};
const char* str_param[nsa];





int main(int argc, char *argv[])
{  
  getCommandLine(int_ident,int_param,nia,str_ident,str_param,nsa,argv,argc);


  Seed      = (long int)(-1* ( int_param[0]  != NOVAL ? int_param[0] : 11041979 ));
  ins       = ( int_param[1 ]  != NOVAL ? int_param[1 ] : 10  );
  timelimit = ( int_param[2 ]  != NOVAL ? int_param[2 ] : -1  );
  Var       = ( int_param[3 ]  != NOVAL ? int_param[3 ] : 20 );
  Dom       = ( int_param[4 ]  != NOVAL ? int_param[4 ] : 10  );
  Con       = ( int_param[5 ]  != NOVAL ? int_param[5 ] : 50 );
  Ari       = ( int_param[6 ]  != NOVAL ? int_param[6 ] : 2   );
  Ngd       = ( int_param[7 ]  != NOVAL ? int_param[7 ] : 30  );
  Ext       = ( int_param[8 ]  != NOVAL ? int_param[8 ] : 2   );
  verbose   = ( int_param[9 ]  != NOVAL ? int_param[9 ] : 0   );


  int nSAT = 1;
  vector<int> distance;
  int allConstraints[(Con+Ext)*Ari];
  int mngd = 1;
  for(int i=0; i<Ari; ++i)
    mngd *= Dom;

  VarArray scope(Ari);
  int allNogoods[(Con+Ext)*mngd*Ari];

  while( nSAT ) {
    nSAT = 0;
    distance.clear();

    URCSPGenerator gen(Seed, Var, Dom, (Con+Ext), Ari, Ngd);

    if( verbose > 1 )  
      cout<<endl<<"\t"<<ins<<" Random Problem(s)"
	  <<"\n\t\t["<<Var<<" variables, "<<Dom<<" values, "
	  <<Con<<" "<<Ari<<"-ary consraints, "<<Ngd
	  <<" nogoods] +- " << Ext << " constraints"<<endl
	  <<endl;

    for(int inst=0; inst<ins; ++inst) {
            
      int i,j,k=0;
      while( gen.erateConstraint() ) {
	for(i=0; i<Ari; ++i) 
	  allConstraints[k*Ari+i] = gen.vars[i+1];
	j=0;
	while( gen.erateNogood() ) {
	  for(i=0; i<Ari; ++i) 
	    allNogoods[k*Ngd*Ari+j*Ari+i] = gen.vals[i];
	  ++j;
	}
	++k;
      }
          
      CSP model1;
      VarArray X1(Var, 0, Dom-1);
      model1.add( X1 );
      for(k=0; k<Con; ++k) {
	for(i=0; i<Ari; ++i) 
	  scope[i] = X1[allConstraints[k*Ari+i]];
	Table con(scope, false);
	for(j=0; j<Ngd; ++j)
	  con.add( &(allNogoods[k*Ngd*Ari+j*Ari]) );
	model1.add( con );
      }    
      Solver solver1(model1);   
      DomOverWDeg heuristic1;
      solver1.add( heuristic1 );
      int status = solver1.solve();
      
      if( verbose )
	cout << setw(4) << inst << " " << setw(3) << Var << " " << setw(2) << Dom << " " << setw(3) << Con 
	     << " " << setw(1) << Ari << " " << setw(3) << Ngd << " " 
	     << setw(1) << Ext << " " ;
      
      if( verbose ) {
	solver1.printStatistics(cout, (OUTCOME | 
				       TOTALTIME | 
				       PPGS |
				       NDS) );
	cout.flush();
      }
      
      if( status == SAT ) {
      
	if( verbose > 1 ) {
	  for(i=0; i<Var; ++i) 
	    cout << (solver1.solution[i]) << " ";
	  cout << endl;
	  
	  solver1.printStatistics(cout, (OUTCOME | 
					 TOTALTIME | 
					 SPEEDB |
					 PPGS |
					 NDS) );
	  cout << endl;
	}
	
	
	CSP model2;
	VarArray scope(Ari);
	VarArray X2(Var, 0, Dom-1);
	model2.add( X2 );
	for(k=Ext; k<Con+Ext; ++k) {
	  for(i=0; i<Ari; ++i) 
	    scope[i] = X2[allConstraints[k*Ari+i]];
	  Table con(scope, false);
	  for(j=0; j<Ngd; ++j)
	    con.add( &(allNogoods[k*Ngd*Ari+j*Ari]) );
	  model2.add( con );
	}
	Solver solver2(model2);   
	DomOverWDeg heuristic2;
	solver2.add( heuristic2 );
	solver2.setVerbosity( verbose ? verbose-1 : 0 );
	int status = solver2.ldSolve( solver1.solution );
      
	nSAT += ( status == SAT );

	if( verbose > 1 ) {
	  for(i=0; i<Var; ++i) 
	    cout << (solver2.solution[i]) << " ";
	  cout << endl;
	  
	  solver2.printStatistics(cout, (OUTCOME | 
					 TOTALTIME | 
					 SPEEDB |
					 PPGS |
					 NDS) );
	  cout << endl;
	}

      
	if( verbose ) {
	  cout << " ";
	  solver2.printStatistics(cout, (OUTCOME | 
					 TOTALTIME | 
					 PPGS |
					 NDS) );
      
	  cout << " " << setw(5) << solver2.DISCREPANCY ;
	}

	distance.push_back( solver2.DISCREPANCY );
	
      } 
      if( verbose )
	cout << endl;


      //cout << "reinit" << endl;
      gen.reInit(Var, Dom, (Con+Ext), Ari, Ngd);
    }

    int i;
    int j, k, n=distance.size();
    for(i=1; i<n; ++i)
      {
	j = i;
	while( j > 0 && distance[j] < distance[j-1] )
	  {
	    k = distance[j];
	    distance[j] = distance[--j];
	    distance[j] = k;
	  }
      }
//     for(i=0; i<n; ++i)
//       cout << " " << distance[i];
//     cout << endl;

    double avgDist = 0;
    int medDist = distance[((distance.size()-1)/2)], maxDist = distance.back(), minDist = distance[0];
    for(i=0; i<n; ++i)
      avgDist += distance[i];
    avgDist /= distance.size();

    double tightness = (double)Ngd;
    for(i=0; i<Ari; ++i)
      tightness /= (double)Dom;

    cerr << Ngd << " " << tightness << " " 
	 << minDist << " " << medDist << " "
	 << avgDist << " " << maxDist << endl;

    ++Ngd;
  }
    
 
}
