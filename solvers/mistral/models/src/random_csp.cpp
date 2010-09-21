
#include <mistral_sol.h>
#include <mistral_gen.h>
#include <mistral_glo.h>

using namespace std;
using namespace Mistral;


int timelimit, ins, Var, Dom, Con, Ari, Ngd, Base, pi, pl, 
  probing, sac, rdz, lds, si;
bool ulw, uip;
string Heu, Pol, ValHeu;
long Seed;
double Factor, Decay;

void addHeuristic( Solver& s ) {
  if( Heu == "dom" ) {
    MinDomain h(abs(rdz));
    s.add( h );
  }
  if( Heu == "lex" ) {
    Lexicographic h;
    s.add( h );
  } 
  else if( Heu == "deg") {
    MaxDegree h(abs(rdz));
    s.add( h );
  } 
  else if( Heu == "rand") {
    Random h;
    s.add( h );
  } 
  else if( Heu == "dom+deg") {
    MinDomMaxDeg h(abs(rdz));
    s.add( h );
  } 
  else if( Heu == "dom/deg") {
    DomOverDeg h(abs(rdz));
    s.add( h );
  } 
  else if( Heu == "dom/wldeg") {
    DomOverWLDeg h(abs(rdz));
    s.add( h );
  }
  else if( Heu == "dom/wdeg") {
    DomOverWDeg h(abs(rdz));
    s.add( h );
  }
  else if( Heu == "neighbor") {
    Neighbor h(abs(rdz));
    s.add( h );
  } 
  else if( Heu == "impact") {
    Impact h(abs(rdz));
    s.add( h );
  }
  else if( Heu == "impact/deg") {
    ImpactOverDeg h(abs(rdz));
    s.add( h );
  }
  else if( Heu == "impact/wdeg") {
    ImpactOverWDeg h(abs(rdz));
    s.add( h );
  }
  else if( Heu == "impact/wldeg") {
    ImpactOverWLDeg h(abs(rdz));
    s.add( h );
  }
  else {
    NoOrder h;
    s.add( h );
  }
  
  if(ValHeu == "anti") {
    s.setAntiLex();
  } else if(ValHeu == "rand") {
    s.setRandMinMax();
  } else if(ValHeu == "split") {
    s.setSplit();
  } else if(ValHeu == "rsplit") {
    s.setRandSplit();
  } else if(ValHeu != "lex") {
    std::cerr << "Warning: unknown value selection heuristic" << std::endl;
  } 

}


const int nia = 18;
const char* int_ident[nia] = {"-seed", "-iterations", "-timelimit", "-variable", "-domain", "-constraint", "-arity", "-nogood", "-base", "-probe_iterations", "-probe_limit", "-probe", "-update_lweight", "-update_impact", "-sac", "-randomized", "-lds", "-single_instance"};
int int_param[nia];

const int nsa = 5;
const char* str_ident[nsa] = {"-heuristic", "-restart", "-factor", "-decay", "-value"};
const char* str_param[nsa];





int main(int argc, char *argv[])
{  

  getCommandLine(int_ident,int_param,nia,str_ident,str_param,nsa,argv,argc);


  Seed      = (long int)(-1* ( int_param[0]  != NOVAL ? int_param[0] : 11041979 ));//2 ));
  ins       = ( int_param[1 ]  != NOVAL ? int_param[1 ] : 3   );//10  );
  timelimit = ( int_param[2 ]  != NOVAL ? int_param[2 ] : -1  );
  Var       = ( int_param[3 ]  != NOVAL ? int_param[3 ] : 120 );//100 );
  Dom       = ( int_param[4 ]  != NOVAL ? int_param[4 ] : 10  );//20  );
  Con       = ( int_param[5 ]  != NOVAL ? int_param[5 ] : 400 );//250 );
  Ari       = ( int_param[6 ]  != NOVAL ? int_param[6 ] : 2   );
  Ngd       = ( int_param[7 ]  != NOVAL ? int_param[7 ] : 46  );//263  );
  Base      = ( int_param[8 ]  != NOVAL ? int_param[8 ] : 32  );
  pi        = ( int_param[9 ]  != NOVAL ? int_param[9 ] : 100 );
  pl        = ( int_param[10]  != NOVAL ? int_param[10] : 30  );
  probing   = ( int_param[11]  != NOVAL ? int_param[11] : 0   );
  ulw       = ( int_param[12]  != NOVAL ? int_param[12] : 0   );
  uip       = ( int_param[13]  != NOVAL ? int_param[13] : 0   );
  sac       = ( int_param[14]  != NOVAL ? int_param[14] : 0   );
  rdz       = ( int_param[15]  != NOVAL ? int_param[15] : 0   );
  lds       = ( int_param[16]  != NOVAL ? int_param[16] : 0   );
  si        = ( int_param[17]  != NOVAL ? int_param[17] : -1  );
  
  Heu       = ( strcmp(str_param[0],"nil") ? str_param[0] : "dom/wdeg" );
  Pol       = ( strcmp(str_param[1],"nil") ? str_param[1] : "no"      );		

  Factor    = ( strcmp(str_param[2],"nil") ? atof(str_param[2]) : 1.3 );
  Decay     = ( strcmp(str_param[3],"nil") ? atof(str_param[3]) : 0.0 );
  ValHeu    = ( strcmp(str_param[4],"nil") ? str_param[4] : "lex" );


  
  URCSPGenerator gen(Seed, Var, Dom, Con, Ari, Ngd);
  
  cout<<endl<<"\t"<<ins<<" Random Problem(s)"
      <<"\n\t\t["<<Var<<" variables, "<<Dom<<" values, "
      <<Con<<" "<<Ari<<"-ary consraints, "<<Ngd
      <<" nogoods]"<<endl
      << "\t heuristic = " << Heu << endl
      << "\t Restart policy = " << Pol << endl
      <<endl;


  unsigned int avg_bts=0, avg_nds=0, avg_fls=0, avg_pgs=0, solved=0, sat=0;
  double avg_time=0;

  for(int inst=0; inst<ins; inst++) {

    CSP model;
    VarArray X(Var, 0, Dom-1);
    VarArray scope(Ari);

    model.add( X );

    while( gen.erateConstraint() ) {
      for(int i=0; i<Ari; ++i) 
	scope[i] = X[gen.vars[i+1]];
      Table con(scope, false);
      while( gen.erateNogood() ) 
	con.add(gen.vals);

      model.add( con );
    }
    
    if( si < 0 || si == inst ) {
      Solver s(model);   
      
      //s.printXML( cout );
      //cout << endl;


 //      ReversibleNum<int> testint;
//       testint.init( &s, 0 );
//       cout << (testint) << endl;

//       ++testint;
//       cout << "++testint " << (testint) << endl;


//       testint*=3;
//       cout << "testint*=3 " << (testint) << endl;


//       testint+=5;
//       cout << "testint+=5 " << (testint) << endl;

//       testint/=4;
//       cout << "testint/=4 " << (testint) << endl;



//       exit( 0 );
      
      
      s.setRandomized( (rdz != 0) );
      
      int Policy = NOVAL;
      if( Pol == "geom" ) {
	Policy = GEOMETRIC;
      } else if( Pol == "luby" ) {
	Policy = LUBY;
      }
      
      if( sac )
	s.sacPreprocess( (sac > 1) );
      
      int status = UNKNOWN;
      
      
      
      if( probing ) {
// 	int wtype = Weighter::WDG;
// 	if( ulw ) 
// 	  wtype |= Weighter::WLD;
// 	if( uip )
// 	  wtype |= Weighter::IPT;
	
	Probedvo heuristic(ulw, uip, 0);
	s.add( heuristic );
	
	status = s.random_probe(pi, pl);
      } 
      
      s.setVerbosity(0);
      

//       s.print( cout );
//       cout << endl;

      if( status == UNKNOWN )
	{
	  s.setTimeLimit( timelimit );
	  addHeuristic(s);
	  if( Policy != NOVAL ) 
	    status = s.solve_and_restart(Policy, Base, Factor, Decay);
	  else if( lds )
	    status = s.ldSolve();
	  else
	    status = s.solve();
	}
      
      if( status != LIMITOUT ) {
	++solved;
	if( status == SAT )
	  ++sat;
      }
      
      s.printStatistics(cout, (OUTCOME | 
			       TOTALTIME | 
			       SPEEDB |
			       SPEEDP |
			       PPGS |
			       NDS | BTS) );
      cout << endl;
      
      avg_bts += s.getBacktracks();
      avg_nds += s.getNodes();
      avg_fls += s.getFailures();
      avg_pgs += s.getPropags();
      avg_time += s.getTime();
    }

    gen.reInit(Var, Dom, Con, Ari, Ngd);
  }

  cout << endl << setw(4) << ((double)solved/(double)ins) << "/"
       << setw(4) << ((double)sat/(double)solved) 
    //<< setw(9) << (avg_bts/ins) << " BTS"
       << setw(9) << (avg_nds/ins) << " NDS"
    //<< setw(9) << (avg_fls/ins) << " FAILS"
       << setw(9) << ( avg_time ? (unsigned long int)((double)avg_bts/avg_time) : 0 ) << " BTS/s"
       << setw(14) << (avg_pgs/ins) << " PPGS"
       << setw(8) << (avg_time/ins) << " s" << endl << endl;

}

// 100 Random Problem(s)      [100 variables, 10 values, 250 2-ary consraints, 55 nogoods]
// neighbor     1/0.33      367 BTS      380 NDS      371 FAILS     4731 BTS/s       1139426 CKS   0.078 s
// dom/wdeg     1/0.33      447 BTS      460 NDS      452 FAILS     4952 BTS/s       1468505 CKS    0.09 s
// dom/wldeg    1/0.33      479 BTS      492 NDS      484 FAILS     4905 BTS/s       1597022 CKS   0.098 s
// impact/wldeg 1/0.33      516 BTS      530 NDS      521 FAILS     4643 BTS/s       1635613 CKS    0.11 s
// impact/wdeg  1/0.33      539 BTS      553 NDS      545 FAILS     4730 BTS/s       1667841 CKS    0.11 s
// dom/deg      1/0.33      607 BTS      620 NDS      611 FAILS     5294 BTS/s       1875343 CKS    0.11 s
// impact       1/0.33     1829 BTS     1843 NDS     1834 FAILS     5258 BTS/s       5055444 CKS    0.35 s
// min domain   1/0.33    37666 BTS    37680 NDS    37670 FAILS     6124 BTS/s      13880388 CKS     6.2 s


////////// OLD IMPLEMENTATION
// bit vars / filo
// 100 Random Problem(s)      [120 variables, 10 values, 400 2-ary consraints, 46 nogoods] (seed: 44444)
// probes       1/0.53    30708 BTS    31792 NDS    30716 FAILS     4484 BTS/s      22192410 CKS     6.8 s  
// neighbor     1/0.53    31081 BTS    31107 NDS    31092 FAILS     4088 BTS/s      19938774 CKS     7.6 s
// dom/wdeg     1/0.53    38237 BTS    38261 NDS    38245 FAILS     4309 BTS/s      12424051 CKS     8.9 s


// list vars / filo
// 100 Random Problem(s)      [120 variables, 10 values, 400 2-ary consraints, 46 nogoods] (seed: 44444)
// probes       1/0.53    30708 BTS    31792 NDS    30716 FAILS     7303 BTS/s      20409833 CKS  4.2047 s
// neighbor     1/0.53    31081 BTS    31107 NDS    31092 FAILS     6170 BTS/s      20020918 CKS  5.0368 s
// dom/wdeg     1/0.53    38237 BTS    38261 NDS    38245 FAILS     7000 BTS/s      26584207 CKS  5.4619 s


// list vars / fifo
// 100 Random Problem(s)      [120 variables, 10 values, 400 2-ary consraints, 46 nogoods] (seed: 44444)
// probes       1/0.53    31839 BTS    32927 NDS    31847 FAILS     8574 BTS/s      16240080 CKS  3.7134 s
// neighbor     1/0.53    31082 BTS    31108 NDS    31093 FAILS     7082 BTS/s      15098776 CKS  4.3883 s
// dom/wdeg     1/0.53    42630 BTS    42655 NDS    42639 FAILS     8232 BTS/s      22816663 CKS  5.1785 s
// dom/wldeg    1/0.53    41122 BTS    41147 NDS    41131 FAILS     8248 BTS/s      21957240 CKS  4.9853 s
// sac+dom/wdeg 1/0.53    34064 BTS    42775 NDS    35194 FAILS     7397 BTS/s      20128777 CKS  4.6047 s
// sac+dom/wdeg 1/0.53    34692 BTS    40068 NDS    35381 FAILS     7659 BTS/s      19788762 CKS  4.5296 s



//intersection: 27.2
// 	19.9%	19.9%	random_csp	MistralSet::wordIntersect(MistralSet const&) const	
// 	7.3%	7.3%	random_csp	VariableDomain::wordIntersect(MistralSet const&) const	

//domain iteration: 16.7
// 	16.7%	16.7%	random_csp	VariableList::revise(MistralSet*, int*, VariableInt*)	

//search/trailing/backtrack 16.6
// 	4.5%	4.5%	random_csp	VariableInt::unLink()	
// 	3.2%	3.2%	random_csp	VariableInt::link()	
// 	2.5%	2.5%	random_csp	ReversibleDList::save()	
// 	2.1%	2.1%	random_csp	ReversibleDList::restore()	
// 	0.0%	0.8%	random_csp	 (invert) ReversibleDList::restore()	
// 	1.3%	1.3%	random_csp	ReversibleInt::restore()	
// 	1.2%	1.2%	random_csp	Solver::depthFirstSearch()	
// 	1.0%	1.0%	random_csp	Solver::save(ReversibleObj*)	

//domain modification: 14.7
// 	8.6%	8.6%	random_csp	VariableList::remove(int)	
// 	0.0%	1.2%	random_csp	 (invert) VariableList::remove(int)	
// 	0.0%	1.3%	random_csp	 (member) VariableDomain::contain(int) const	
// 	0.0%	0.1%	random_csp	 (member) VariableList::remove(int)	
// 	0.0%	1.3%	random_csp	 (RevInt operator=) VariableList::remove(int)	
// 	0.0%	0.8%	random_csp	 (contain) VariableList::remove(int)	
// 	0.0%	0.7%	random_csp	 (next) VariableList::remove(int)	
// 	0.0%	0.7%	random_csp	 (prev) VariableList::remove(int)	

//AC3 11.6
// 	11.5%	11.5%	random_csp	Solver::filtering()	
// 	0.0%	0.1%	random_csp	 Solver::filtering()	

//revision condition + virtual call: 10.6
// 	7.8%	7.8%	random_csp	ConstraintBinaryExtensional::propagate(int, int)
// 	0.0%	2.8%	random_csp	 (domsize) ConstraintBinaryExtensional::propagate(int, int)	

//variable heuristic (dom/wldeg) 0.9
// 	0.0%	0.3%	random_csp	 (domsize) GenericDVO<VarSelectorDomainOverWeight>::select()	
// 	0.6%	0.6%	random_csp	GenericDVO<VarSelectorDomainOverWeight>::select()	




// 	19.9%	19.9%	random_csp	MistralSet::wordIntersect(MistralSet const&) const	
// 	16.7%	16.7%	random_csp	VariableList::revise(MistralSet*, int*, VariableInt*)	
// 	11.5%	11.5%	random_csp	Solver::filtering()	
// 	8.6%	8.6%	random_csp	VariableList::remove(int)	
// 	7.8%	7.8%	random_csp	ConstraintBinaryExtensional::propagate(int, int)
// 	7.3%	7.3%	random_csp	VariableDomain::wordIntersect(MistralSet const&) const	
// 	4.5%	4.5%	random_csp	VariableInt::unLink()	
// 	3.2%	3.2%	random_csp	VariableInt::link()	
//// 	3.1%	3.1%	random_csp	VariableList::domsize() const	
// 	0.0%	2.8%	random_csp	 ConstraintBinaryExtensional::propagate(int, int)	
// 	0.0%	0.3%	random_csp	 GenericDVO<VarSelectorDomainOverWeight>::select()	
// 	0.0%	0.1%	random_csp	 Solver::filtering()	
// 	2.5%	2.5%	random_csp	ReversibleDList::save()	
// 	2.1%	2.1%	random_csp	ReversibleDList::restore()	
//// 	2.1%	2.1%	random_csp	MistralSet::fastInvert(int)	
// 	0.0%	1.2%	random_csp	 VariableList::remove(int)	
// 	0.0%	0.8%	random_csp	 ReversibleDList::restore()	
// 	0.0%	0.0%	random_csp	 Solver::depthFirstSearch()	
// 	0.0%	0.0%	random_csp	 VariableList::revise(MistralSet*, int*, VariableInt*)	
//// 	1.4%	1.4%	random_csp	MistralSet::member(int) const	
// 	0.0%	1.3%	random_csp	 VariableDomain::contain(int) const	
// 	0.0%	0.1%	random_csp	 VariableList::remove(int)	
//// 	1.3%	1.3%	random_csp	ReversibleInt::operator=(int)	
// 	0.0%	1.2%	random_csp	 VariableList::remove(int)	
// 	0.0%	0.0%	random_csp	 VariableList::setDomain(int)	
// 	0.0%	0.0%	random_csp	 VariableList::revise(MistralSet*, int*, VariableInt*)	
// 	1.3%	1.3%	random_csp	ReversibleInt::restore()	
// 	1.2%	1.2%	random_csp	Solver::depthFirstSearch()	
// 	1.0%	1.0%	random_csp	Solver::save(ReversibleObj*)	
//// 	0.8%	0.8%	random_csp	VariableDomain::contain(int) const	
// 	0.0%	0.8%	random_csp	 VariableList::remove(int)	
//// 	0.7%	0.7%	random_csp	MistralSet::next(int) const	
// 	0.0%	0.7%	random_csp	 VariableList::remove(int)	
// 	0.0%	0.0%	random_csp	 VariableList::revise(MistralSet*, int*, VariableInt*)	
//// 	0.7%	0.7%	random_csp	MistralSet::prev(int) const	
// 	0.0%	0.7%	random_csp	 VariableList::remove(int)	
// 	0.6%	0.6%	random_csp	GenericDVO<VarSelectorDomainOverWeight>::select()	



