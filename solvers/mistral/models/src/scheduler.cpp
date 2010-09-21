 
#include <mistral_sol.h>
#include <iostream>
#include <fstream>

#include "scheduler.hpp"





int main( int argc, char** argv )
{
  if( argc < 2 )
    {
      cerr << "need a data file" << endl;
      exit( 0 );
    }

  int i=0;
  while(argv[1][i] != '\0') ++i;
  while(argv[1][i] != '/') --i;

  data_file = &(argv[1][i+1]);

  getCommandLine(int_ident,int_param,nia,str_ident,str_param,nsa,&argv[1],argc-1);

  Seed        = ( int_param[0]  != NOVAL ? int_param[0] : 11041979 );
  Dichotomy   = ( int_param[2]  != NOVAL ? int_param[2] : 64  );
  Base        = ( int_param[3]  != NOVAL ? int_param[3] : 256 );
  SAC         = ( int_param[4]  != NOVAL ? int_param[4] : 0   );
  Randomized  = ( int_param[5]  != NOVAL ? int_param[5] : 1   );

  Verbose     = ( int_param[6]  != NOVAL ? int_param[6] : 1   );
  Step        = ( int_param[7]  != NOVAL ? int_param[7] : 1   );
  Weight      = ( int_param[8]  != NOVAL ? int_param[8] : 0   );
  Optimise    = ( int_param[9]  != NOVAL ? int_param[9] : 3600);
  Rngd        = ( int_param[10] != NOVAL ? int_param[10]: 1   );
  Gap         = ( int_param[11] != NOVAL ? int_param[11]: 0   );
  Reset       = ( int_param[12] != NOVAL ? int_param[12]: 0   );
  Wprof       = ( int_param[13] != NOVAL ? int_param[13]: 1   );
  Multiplier  = ( int_param[14] != NOVAL ? int_param[14]: 1   );
  JobRank     = ( int_param[15] != NOVAL ? int_param[15]: 1   ); 
  Pool        = ( int_param[16] != NOVAL ? int_param[16]: 1   ); 
  Rprob       = ( int_param[17] != NOVAL ? int_param[17]: 1000); 
  UBinit      = ( int_param[18] != NOVAL ? int_param[18]: -1  ); 
  //OSP_Type    = ( int_param[19] != NOVAL ? int_param[19]: 1  ); 
  NBJ         = ( int_param[20] != NOVAL ? int_param[20]: 0  ); 


  Policy    = ( strcmp(str_param[1],"nil") ? str_param[1] : "geom"      );		
  Factor    = ( strcmp(str_param[2],"nil") ? atof(str_param[2]) : 1.3 );
  Decay     = ( strcmp(str_param[3],"nil") ? atof(str_param[3]) : 0.0 );
  Type      = ( strcmp(str_param[4],"nil") ? str_param[4] : "jsp" );
  Value     = ( strcmp(str_param[5],"nil") ? str_param[5] : "guided" );
  Algo      = ( strcmp(str_param[6],"nil") ? str_param[6] : "b&b" );
  DValue    = ( strcmp(str_param[7],"nil") ? str_param[7] : "guided" );
  IValue    = ( strcmp(str_param[8],"nil") ? str_param[8] : "promise" );
  Skew      = ( strcmp(str_param[9],"nil") ? atof(str_param[9]) : -1.0 );
//  Key       = ( strcmp(str_param[10],"nil") ? str_param[10] : "nokey" );

  Cutoff    = ( int_param[1]  != NOVAL ? int_param[1] : (Type == "osp" || Type == "sds" ? 3 : 300) );
  Heuristic = ( strcmp(str_param[0],"nil") ? str_param[0] : 
		(Type == "pfsp" ?  "pfsp" : 
		 (Type == "jtl0" ? "osp-dw" : 
		  ((Type == "jsp" || Type == "fsp" || Type == "jtl0") ? "osp-t" : "osp-b"))) );
  
  //if(  Value == "guided" )   
  VO = GUIDED;
  if(  Value == "rguided")   VO = RGUIDED;
  if(  Value == "lex"    )   VO = LEX;
  if(  Value == "rand"   )   VO = RAND;
  if(  Value == "promise")   VO = PROMISE;
  if(  Value == "anti"   )   VO = ANTI;
  
  //if( DValue == "guided" ) 
  D_VO = GUIDED;
  if( DValue == "rguided") D_VO = RGUIDED;
  if( DValue == "lex"    ) D_VO = LEX;
  if( DValue == "rand"   ) D_VO = RAND;
  if( DValue == "promise") D_VO = PROMISE;
  if( DValue == "anti"   ) D_VO = ANTI;

  //if( IValue == "promise") 
  I_VO = PROMISE;
  if( IValue == "lex"    ) I_VO = LEX;
  if( IValue == "rand"   ) I_VO = RAND;
  if( IValue == "anti"   ) I_VO = ANTI;
    

  usrand( Seed );

  //if(Verbose > 1) {
  cout << endl << "c =================[ parameters ]==================" << endl;
  cout << left << setw(30) << "c seed " << ":" << right << setw(20) << Seed << endl;
  cout << left << setw(30) << "c cutoff " << ":" << right << setw(20) << Cutoff << endl;
  cout << left << setw(30) << "c dichotomy " << ":" << right << setw(20) << (Dichotomy ? "yes" : "no") << endl;
  cout << left << setw(30) << "c reuse learned weights " << ":" << right << setw(20) << (Weight ? "yes" : "no") << endl;
  cout << left << setw(30) << "c restart policy " << ":" << right << setw(20) << Policy << endl;
  cout << left << setw(30) << "c base " << ":" << right << setw(20) << Base << endl;
  cout << left << setw(30) << "c factor " << ":" << right << setw(20) << Factor << endl;
  cout << left << setw(30) << "c sac " << ":" << right << setw(20) << ( SAC ? (SAC==1 ? "approx" : "full") : "no") << endl ;
  cout << left << setw(30) << "c alogrithm " << ":" << right << setw(20) << Algo << endl;
  cout << left << setw(30) << "c heuristic " << ":" << right << setw(20) << Heuristic << " (" << abs(Randomized) << ")" << endl;
  cout << left << setw(30) << "c jobrank " << ":" << right << setw(20) << JobRank << endl;
  cout << left << setw(30) << "c value ordering (init step) " << ":" << right << setw(20) << IValue << endl;
  cout << left << setw(30) << "c value ordering (dichotomy) " << ":" << right << setw(20) << DValue << endl;
  cout << left << setw(30) << "c value ordering (optim) " << ":" << right << setw(20) << Value << endl;
  cout << "c =================[ parameters ]==================" << endl;
  //}

  PolicyRestart = NO;
  if( Policy == "geom" )
    PolicyRestart = GEOMETRIC;
  else if( Policy == "luby" )
    PolicyRestart = LUBY;

  if(Type == "osp") {
    osp_readData( argv[1] );
  } else if(Type == "sds") {
    sds_readData( argv[1] );
  } else if(Type == "jtl") {
    jtl_readData( argv[1] );
  } else if(Type == "jtl0") {
    jtl_readData( argv[1] );
  } else if(Type == "jla") {
    jla_readData( argv[1] );
  } else if(Type == "tsp") {
    tsp_readData( argv[1] );
  } else if(Type == "fsp") {
    fsp_readData( argv[1] );
  } else if(Type == "pfsp") {
    fsp_readData( argv[1] );
  } else if(Type == "jsp") {
    jsp_readData( argv[1] );
  }

  total_time = getRunTime();

  maxfsble = UBinit;
  if(Type == "osp") {
    if(UBinit < 0) maxfsble = osp_upperbound();
    minfsble = osp_lowerbound();
  } else if(Type == "jtl0") {
    if(UBinit < 0) maxfsble = jsp_upperbound();
    minfsble = jsp_lowerbound();
  } else { //if(Type == "jsp") {
    if(UBinit < 0) maxfsble = jsp_upperbound();
    minfsble = jsp_lowerbound();
  }
  max_infeasible = minfsble-1;

  init_lb = minfsble;
  init_ub = maxfsble;

  best_solution  = new int[ndisjuncts];
  probability    = new int[ndisjuncts];
  range          = new int[ndisjuncts];
  zeros          = new int[ndisjuncts];
  ones           = new int[ndisjuncts];
  disjunct_index = new int[ndisjuncts];
  
  std::fill(best_solution, best_solution+ndisjuncts, 0);
  std::fill(probability, probability+ndisjuncts, Tprob/2);
  std::fill(range, range+ndisjuncts, Tprob);
  std::fill(zeros, zeros+ndisjuncts, Tprob/2);
  std::fill(ones , ones +ndisjuncts, Tprob/2);

  solutions.push_back(best_solution);

  cout << endl << "c ===================[ bounds ]====================" << endl;
  cout 
       << left << setw(30) << "c initial bounds" << ":" 
       << right << setw(6) << " " << setw(5) << minfsble << " to " << setw(5) << maxfsble << endl
       << left << setw(30) << "c values in file" << ":" 
       << right << setw(6) << " " << setw(5) << lb << " to " << setw(5) << opt 
       << endl ;
  if(Dichotomy) {
    if(Skew < 0) Skew = ((double)(maxfsble-minfsble)/(double)minfsble);
    if(Skew != 1.0) {
      if(Skew < 1.0) Skew = 1.0; 
      double rs = (randreal() - .5)/2.0;
      cout << left << setw(30) << "c Dichotomic skew " << ":" ;
      if(rs<0) cout << right << setw(6) << setprecision(2) << (-rs) << " - " ;
      else cout << right << setw(6) << setprecision(1) << rs << " + ";
      cout << right << setw(4) << setprecision(3) << Skew;
      Skew += rs;
      cout << " = " << right << setw(4) << Skew << endl;
    } else {
      cout << left << setw(30) << "c Dichotomic skew " << ":" << right << setw(20) <<  "no" << endl;
    }
  }
  cout << "c ===================[ bounds ]====================" << endl << endl;
  
  if( Dichotomy < 0)
    mks_dec_search();
  else if(NBJ)
    nbj_dichotomic_search();
  else
    dichotomic_search();
	
  if( Optimise && !Solved ) {
    if(Algo == "b&b") {
      if(NBJ) nbj_branch_and_bound();
      else branch_and_bound();
    } else if(Algo == "cds")
      climbing_discrepancy_search();
    else 
      cout << "c stop after dichotomic search (unknown algo?)" << endl;
  } 
  
  total_time = (getRunTime() - total_time);
  
  cout << "c =================[ statistics ]==================" << endl;
  cout << left << setw(30) << "s" << right << setw(21) << (Solved ? "SATISFIABLE" : "UNSATISFIABLE" ) << endl 
       << left << setw(30) << "d INITUB "      << right << setw(21) << init_ub << endl
       << left << setw(30) << "d INITLB "      << right << setw(21) << init_lb << endl
       << left << setw(30) << "d LOWERBOUND "  << right << setw(21) << (max_infeasible+1) << endl
       << left << setw(30) << "d OBJECTIVE "   << right << setw(21) << maxfsble << endl 
       << left << setw(30) << "d BESTKNOWN "   << right << setw(21) << opt << endl
       << left << setw(30) << "d RUNTIME "     << right << setw(21) << total_time << endl
       << left << setw(30) << "d OPTTIME "     << right << setw(21) << opt_time << endl	
       << left << setw(30) << "d PROOFTIME "   << right << setw(21) << proof_time << endl	
       << left << setw(30) << "d NODES "       << right << setw(21) << total_nodes << endl
       << left << setw(30) << "d BACKTRACKS "  << right << setw(21) << total_bts << endl
       << left << setw(30) << "d FAILS "       << right << setw(21) << total_fails << endl
       << left << setw(30) << "d PROPAGS "     << right << setw(21) << total_propags << endl
       << left << setw(30) << "d RESTARTS "    << right << setw(21) << total_restarts << endl
       << left << setw(30) << "d OPTIMAL "     << right << setw(21) << Solved << endl;
    //<< "d NOGOODSIZE " << (double)(nogood_size)/(double)(nb_nogood) << endl
  cout << "c =================[ statistics ]==================" << endl
       << endl;

  for(int i=0; i<nJobs; ++i) {
    delete [] duration[i];
    if(machine)
      delete [] machine[i];
    if(due_date) {
      delete [] due_date[i];
      delete [] release_date[i];
    }
  }

  if(setup_time)
    for(int i=0; i<nMachines; ++i) {
      for(int j=0; j<nJobs; ++j) 
	delete [] setup_time[i][j];
      delete [] setup_time[i];
    }

  delete [] probability;
  delete [] range;
  delete [] zeros;
  delete [] ones;

  delete [] duration;
  delete [] machine;
  delete [] jsp_duration;
  delete [] jsp_machine;

  for(unsigned int i=0; i<solutions.size(); ++i)
    delete [] solutions[i];

  for(unsigned int i=0; i<disjunct_weights.size(); ++i)
    delete [] disjunct_weights[i];

  for(unsigned int i=0; i<normalised_weights.size(); ++i)
    delete [] normalised_weights[i];

  delete [] disjunct_index;


}

/*
 1 ->  798*
 2 ->  784*
 3 ->  749*
 4 ->  730*
 5 ->  691*
 6 -> 1009*
 7 ->  970*
 8 ->  963*
 9 -> 1060
10 -> 1018*
11 -> 1626
12 -> 1442
13 -> 1544
14 -> 1552
15 -> 1607
*/

  




