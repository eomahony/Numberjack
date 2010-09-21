
#include <mistral_sol.h>
#include <mistral_gen.h>
#include <mistral_glo.h>

#include <sstream>

using namespace std;
using namespace Mistral;




// hack to go from lower to upper case
double STARTTIME, TIME;
int shift;
int *dictionary;

// list of word-slots  
vector< vector < int > > slots;
// white cells are 0s white cells are 1s
vector< vector < int > > grid;

// parameters, accessible from the command line
int Timelimit, Randomized, Base, Height, Width, Sac, All, Policy, Verbose;
string Heu, Pol, Dico, Grid, Squares;
long Seed;
double Factor, Decay;


int toInt( const char c ) {
  //return (int)c - (65+shift);
  int res = (int)c - 65;
  if( res > 25 )
    res -= 32;
  if( res < 0 )
    cout << c << endl;

  assert( res >= 0 );
  assert( res < 26 );
  return res;
}

char toChar( const int i ) {
  return (char)(i+65);
}


string toString( const int i ) {
  std::stringstream out;
  out << i;
  return out.str();
}

void snail( int xlb, int xub, int x,
	    int ylb, int yub, int y,
	    int direction,
	    vector< int >& coords,
	    int count
	    );


/// Custom value selector
class ValSelectorCrossword : public ValSelector
{
public:

  /**@name Parameters*/
  //@{
  double total_proba;
  int val;
  int deg;
  ConstraintGAC2001Allowed **word;
  int *idxWord;
  int *nWords;
  //@}

  /**@name Constructors*/
  //@{
  ValSelectorCrossword( VariableInt *x );
  virtual ~ValSelectorCrossword();
  //@}

  /**@name Utils*/
  //@{
  double getProba() ;
  int getBest() ;
  void left() ;
  void right() ; 
  void postCut( const int p ) ;
  //@}  

  /**@name Miscellanous*/
  //@{
  virtual void printLeft(std::ostream& o) const ;
  virtual void printRight(std::ostream& o) const ;
  //@}
};



void addHeuristic( Solver& s ) {
  if( Heu == "dom" ) {
    if(Verbose > 0) cout << "c using minimum domain variable ordering" << endl;    
    MinDomain h(abs(Randomized));
    s.add( h );
  }
  if( Heu == "lex" ) {
    if(Verbose > 0) cout << "c using lexicographic variable ordering" << endl;
    Lexicographic h;
    s.add( h );
  } 
  else if( Heu == "deg") {
    if(Verbose > 0) cout << "c using maximum degree variable ordering" << endl;
    MaxDegree h(abs(Randomized));
    s.add( h );
  } 
  else if( Heu == "rand") {
    if(Verbose > 0) cout << "c using random variable ordering" << endl;
    Random h;
    s.add( h );
  } 
  else if( Heu == "dom+deg") {
    if(Verbose > 0) cout << "c using min domain -> max degree variable ordering" << endl;
    MinDomMaxDeg h(abs(Randomized));
    s.add( h );
  } 
  else if( Heu == "dom/deg") {
    if(Verbose > 0) cout << "c using minimum domain/degree variable ordering" << endl;
    DomOverDeg h(abs(Randomized));
    s.add( h );
  } 
  else if( Heu == "dom/wldeg") {
    if(Verbose > 0) cout << "c using minimum domain/weighted degree (l) variable ordering" << endl;
    DomOverWLDeg h(abs(Randomized));
    s.add( h );
  }
  else if( Heu == "dom/wdeg") {
    if(Verbose > 0) cout << "c using minimum domain/weighted degree variable ordering" << endl;
    DomOverWDeg h(abs(Randomized));
    s.add( h );
  }
  else if( Heu == "neighbor") {
    if(Verbose > 0) cout << "c using neighbor variable ordering" << endl;
    Neighbor h(abs(Randomized));
    s.add( h );
  } 
  else if( Heu == "impact") {
    if(Verbose > 0) cout << "c using impact variable ordering" << endl;
    Impact h(abs(Randomized));
    s.add( h );
  }
  else if( Heu == "impact/deg") {
    if(Verbose > 0) cout << "c using impact/degree variable ordering" << endl;
    ImpactOverDeg h(abs(Randomized));
    s.add( h );
  }
  else if( Heu == "impact/wdeg") {
    if(Verbose > 0) cout << "c using impact/weighted degree variable ordering" << endl;
    ImpactOverWDeg h(abs(Randomized));
    s.add( h );
  }
  else if( Heu == "impact/wldeg") {
    if(Verbose > 0) cout << "c using impact/weighted degree (l) variable ordering" << endl;
    ImpactOverWLDeg h(abs(Randomized));
    s.add( h );
  }
//   else if( Heu == "proba/wdeg") {
//     if(Verbose > 0) cout << "c using minimum proba/weighted degree variable ordering" << endl;
//     ProbaOverWDeg h(abs(Randomized));
//     s.add( h );
//   }
  else {
    NoOrder h;
    s.add( h );
  }
}


const int nia = 9;
const char* int_ident[nia] = {"-seed", "-timelimit", "-base",  "-randomized", "-height", "-width", "-sac", "-all", "-verbose"};
int int_param[nia];

const int nsa = 6;
const char* str_ident[nsa] = {"-heuristic", "-restart", "-factor", "-dico", "-grid", "-squares"};
const char* str_param[nsa];


/// Initialise black cells from a file
void initGrid()
{
  ifstream grid_data(Grid.c_str(), ios_base::in);
  if(!grid_data.good()) {
    cerr << "No grid found in " << Grid << endl;
    exit(0);
  }
  char line[1000];
  int i, j; 
  
  grid_data.getline( line, 1000 );

  while(true) {
    grid_data.getline( line, 1000 );
    if( !grid_data.good() ) break;
    vector< int > x;
    for(i=0; line[i] != ')'; ++i) {
      if( line[i] == '_' ) x.push_back(1);
      else if( line[i] == '*' ) x.push_back(0);
    }
    if( x.size() > 0 )
      grid.push_back( x );
  }  

  Width = grid[0].size();
  Height = grid.size();

  // horizontal slots:
  for(i=0; i<Height; ++i) {
    for(j=0; j<Width; ++j) {
      if( grid[i][j] ) {
	if(!j || !grid[i][j-1]) {
	  if(slots.size() && slots.back().size() <= 1) {
	    slots.pop_back();
	  }
	  vector< int > x;
	  slots.push_back( x );
	}
	slots.back().push_back( i*Width+j );
      }
    }
  }

  // vertical slots:
  for(j=0; j<Width; ++j) {
    for(i=0; i<Height; ++i) {
      if( grid[i][j] ) {
	if(!i || !grid[i-1][j]) {
	  if(slots.size() && slots.back().size() <= 1) {
	    slots.pop_back();
	  }
	  vector< int > x;
	  slots.push_back( x );
	}
	slots.back().push_back( i*Width+j );
      }
    }
  }

  for(i=0; i<Height; ++i) {
    cout << "c  (";
    for(j=0; j<Width; ++j) {
      cout << (grid[i][j] ? " _" : " *");
    }
    cout << " )"<< endl;
  }

}

/// initialise a grid without black cells
void initBlank()
{
  int i, j;
  for(i=0; i<Height; ++i) {
    vector< int > x;
    for(j=0; j<Width; ++j) {
      x.push_back(1);
    }
    grid.push_back(x);
  }

  for(i=0; i<Height; ++i) {
    vector< int > x;
    for(j=0; j<Width; ++j) {
      x.push_back( i*Width+j );
    }
    slots.push_back(x);
  }
  for(j=0; j<Width; ++j) {
    vector< int > x;
    for(i=0; i<Height; ++i) {
      x.push_back( i*Width+j );
    }
    slots.push_back(x);
  }
}

/// read dictionary from a file, format: one word per line wihtout, either all in lower case or all in upper case
void readDictionary()
{

  TIME = getRunTime();  

  unsigned int i, j, k, edge=std::max(Width, Height)+1;
  int *size_index = new int[edge+1];
  vector<int> word_sizes;
  vector<int> *w_list = new vector<int>[edge+1];


  ifstream allowed(Dico.c_str(), ios_base::in);
  if(!allowed.good()) {
    cerr << "No dictionary found in " << Dico << endl;
    exit(0);
  } else {
    cout << "c " << setw(30) << "read dictionary from" << setw(30) << Dico << "   ";
  }
  char word[100];
  dictionary = new int[edge+1];
  std::fill( size_index, size_index+edge+1, 0);

  for(i=0; (unsigned int)i<slots.size(); ++i) 
    if( ++size_index[slots[i].size()] == 1) word_sizes.push_back(slots[i].size());

  char c = allowed.get();
  if( (int)c >= 97 )
    shift = 32;
  else 
    shift = 0;
  allowed.unget();

  while( true ) {
    allowed.getline( word, 100 );
    if( !allowed.good() ) break;    
    for(i=0; word[i]!='\0' && i<=edge; ++i) {
      j = toInt(word[i]);
      if(j<0 || j>25) break;
    }
    if(i<=edge && size_index[i]) {
      for(k=0; k<i; ++k) {
	j = toInt(word[k]);
	w_list[i].push_back( j );
      }
    }
  }

  for(k=0; k<word_sizes.size(); ++k) {
    dictionary[word_sizes[k]] = CSP::addTable(word_sizes[k], w_list[word_sizes[k]].size()/word_sizes[k]);
    for(i=0; i<w_list[word_sizes[k]].size(); i+=word_sizes[k]) 
      CSP::addTuple( dictionary[word_sizes[k]], &w_list[word_sizes[k]][i] );
  }

  delete [] w_list;
  delete [] size_index;

  cout << setw(6) << setprecision(2) << (getRunTime() - TIME) << "s" << endl;
}


/// add the 'square' implied constraints
void addSquares(VarArray& X, CSP& model) 
{
  TIME = getRunTime();
  cout << "c " << setw(30) << "adding implied squares from" << setw(30) << Squares << "   ";

  int i, j, k, tuple[4];
  char word[100];
  for(i=0; i<Height-1; ++i)
    for(j=0; j<Width-1; ++j)
      {
	string s = (Squares+"/support"+toString(Width)+"x"+toString(Height)+"_"+toString(i)+toString(i+1)+toString(j)+toString(j+1)+".txt");
	//cout << s << endl;
	ifstream squarefile( s.c_str(), ios_base::in );
	VarArray scope;
	scope.add( X[i*Width+j] );
	scope.add( X[i*Width+j+1] );
	scope.add( X[(i+1)*Width+j] );
	scope.add( X[(i+1)*Width+j+1] );
	Table implied(scope, true);
	while( true ) {
	  squarefile.getline( word, 100 );
	  if( !squarefile.good() ) break;    
	  
	  for(k=0; k<4; ++k)
	    tuple[k] = toInt(word[k]);
	  
	  implied.add( tuple );
	}

	//cout << implied << endl;

	model.add( implied );
	squarefile.close();
      }

  cout << setw(6) << setprecision(2) << (getRunTime() - TIME) << "s" << endl;
}


void printGrid(VarArray& X)
{
  int i, j;
  for(i=0; i<Height; ++i) {
    cout << "c  ";
    for(j=0; j<Width; ++j) { 
      if(grid[i][j])
	cout << (toChar(X[i*Width+j].value())) << " ";
      else
	cout << "  ";
    }
    cout << endl;
  }
}


int main(int argc, char *argv[])
{  

  getCommandLine(int_ident,int_param,nia,str_ident,str_param,nsa,argv,argc);


  Seed       = (long int)(-1* ( int_param[0]  != NOVAL ? int_param[0] : 11041979 ));
  Timelimit  = ( int_param[1 ]  != NOVAL ? int_param[1 ] : -1  );
  Base       = ( int_param[2 ]  != NOVAL ? int_param[2 ] : 32  );
  Randomized = ( int_param[3 ]  != NOVAL ? int_param[3 ] : 0   );
  Height     = ( int_param[4 ]  != NOVAL ? int_param[4 ] : 5   );
  Width      = ( int_param[5 ]  != NOVAL ? int_param[5 ] : 5   );
  Sac        = ( int_param[6 ]  != NOVAL ? int_param[6 ] : 0   );
  All        = ( int_param[7 ]  != NOVAL ? int_param[7 ] : 0   );  
  Verbose    = ( int_param[8 ]  != NOVAL ? int_param[8 ] : 1   );  
  
  Heu        = ( strcmp(str_param[0],"nil") ? str_param[0] : "dom/wdeg" );
  Pol        = ( strcmp(str_param[1],"nil") ? str_param[1] : "geom"      );
  Factor     = ( strcmp(str_param[2],"nil") ? atof(str_param[2]) : 1.3 );
  Dico       = ( strcmp(str_param[3],"nil") ? str_param[3] : "./data/dictionary/ogd.txt" );
  Grid       = ( strcmp(str_param[4],"nil") ? str_param[4] : "blank" );
  Squares    = ( strcmp(str_param[5],"nil") ? str_param[5] : "no" );

  unsigned int i, j, k;




  //exit( 0 );


  STARTTIME = getRunTime();


  if( Grid == "blank" )
    initBlank();
  else
    initGrid();

  cout << "c " << Height << " x " << Width << endl;

  readDictionary();

  /// The actual CP model!
  CSP model;
  VarArray X(Height*Width, 0, 25);

  for(i=0; i<slots.size(); ++i) {
    // set an array of variables called 'word' containing all X's of slot i
    VarArray word;
    for(j=0; j<slots[i].size(); ++j)
      word.add( X[slots[i][j]] );
    
    // post the corresponding table constraint (word has to be in dictionary)
    Table mustBeInDico(word, dictionary[slots[i].size()], true);
    model.add( mustBeInDico );
  }
  
  // make sur that all word are distincts
  for(i=1; i<slots.size(); ++i) {
    for(j=0; j<i; ++j) {
      if( slots[i].size() == slots[j].size() ) {
	// an array of Boolean variable to store the reified not-equal constraints
	VarArray differences;
	for(k=0; k<slots[i].size(); ++k) {
	  if( slots[i][k] != slots[j][k] )
	    differences.add( X[slots[i][k]] != X[slots[j][k]] );
	}

	// there must be at least one difference
	model.add( Sum(differences) > 0 );
      }
    }
  }


  if(Squares != "no")
    addSquares(X, model);

 //  VarArray searchVars;
//   vector< int > shell;
//   snail( 0, Width-1, Width-1, 0, Height-1, Height-1, 0, shell, Width*Height );

//   for(i=0; i<shell.size(); ++i)
//     searchVars.add( X[shell[i]] );


//   // search part
//   Solver solver(model, searchVars);



  Solver solver(model, X);


  //solver.setRandomized( Randomized );
  solver.setRandomSeed( Seed );
  solver.randomizeSequence();

  // variable heuristic
  addHeuristic(solver);

  // value heuristic
//   for(i=0; i<solver.length; ++i)
//     solver.sequence[i]->branch = new ValSelectorCrossword( solver.sequence[i] );


//   solver.print( cout );
//   cout << endl;
//   exit(0);

  if( Sac ) {
    cout << "c " << "SAC preprocessing (" << Sac << ")" << endl;
    solver.sacPreprocess(Sac > 1);
  }



  solver.setVerbosity(1);
   
  if( All ) {
    cout << "c " << "solving (print " ;
    if(All < 0) cout << "all";
    else cout << All;
    cout << " solutions)" << endl;    
    solver.startNewSearch();
    while(All-- && solver.getNextSolution() != UNSAT)
      {
	printGrid(X);
	cout << endl;
	solver.printStatistics( cout, RUNTIME | NDS );
	cout << endl;
      }
  } else {

    if( Pol == "no" ) {
      if( Verbose > 0 ) cout << "c no restarts" << endl;
      Policy = NO;
    } else if( Pol == "geom" ) {
      if( Verbose > 0 ) cout << "c geometric restarts (" << Base << " x " << Factor << ")" << endl;
      Policy = GEOMETRIC;
    } else if( Pol == "luby" ) {
      if( Verbose > 0 ) cout << "c Luby restarts (" << Base << ")" << endl;
      Policy = LUBY;
    }
    
    if( Policy != NO ) {
      cout << "c solving (w/ restart)" << endl;
      if(solver.solve_and_restart(Policy, Base, Factor)) {
	printGrid(X);	
	cout << endl << "s SATISFIABLE" << endl;
      } else
	cout << endl << "s UNSATISFIABLE" << endl;
    } else {
      cout << "c solving" << endl;
      if( solver.solve() == SAT) {
	printGrid(X);
	cout << endl << "s SATISFIABLE" << endl;
      } else {
	cout << endl << "s UNSATISFIABLE" << endl;
	printGrid(X);
      }
    }
  }

  cout << "d CHECKS " << solver.PROPAGS << endl
       << "d ASSIGNMENTS " << solver.NODES << endl
       << "d NODES " << solver.NODES 
       << " BACKTRACKS " << solver.BACKTRACKS
       << " FAILURES " << solver.FAILURES
       << " RUNTIME " << solver.ENDTIME
       << " TOTALTIME " << solver.TOTTIME
       << " NODES/s " << (solver.ENDTIME ? (double)(solver.NODES)/solver.ENDTIME : 100*solver.NODES) 
       << " CHECKS/s " << (solver.ENDTIME ? (double)(solver.PROPAGS)/solver.ENDTIME : 100*solver.CHECKS) << endl 
       << "d TOTALTIME " << setw(6) << setprecision(2) << (getRunTime() - STARTTIME) << "s" << endl;
}




ValSelectorCrossword::ValSelectorCrossword( VariableInt *x ) : ValSelector(x)
{
  MistralNode<Constraint*> *nd = _X->constraintsOnDomain();
  int i=0, max_deg=_X->degree;
  word = new ConstraintGAC2001Allowed*[max_deg];
  idxWord = new int[max_deg];
  nWords = new int[max_deg];

  while( nextNode(nd) ) {
    if( nd->elt->name() == "GAC2001" /*&& nd->elt->arity > 4*/ ) {
      word[i] = (ConstraintGAC2001Allowed*)(nd->elt);
      idxWord[i] = nd->index;
      ++i;
    }
    deg = i;
  }

  std::fill(nWords, nWords+max_deg, 0);
  //nWords[0] = nWords[1] = 0;
  DomainIterator *valit = _X->begin();
  do {   
    val = *valit;
    for(i=0; i<deg; ++i) {
      nWords[i] += word[i]->getNumSupports( idxWord[i], val );
    }
  } while( valit->next() );
}

ValSelectorCrossword::~ValSelectorCrossword() 
{
  delete [] nWords;
  delete [] word;
  delete [] idxWord;
}

double ValSelectorCrossword::getProba() {
  getBest();
  return total_proba;
}
 
int ValSelectorCrossword::getBest() {
  total_proba = 0;
  int v;
  double proba=0, best=0;
  DomainIterator *valit = _X->begin();
  do {
    v = *valit;
    proba = 1.0;
    for(int i=0; i<deg; ++i)
      proba *= ((double)(word[i]->getNumSupports(idxWord[i],v)) / (double)(nWords[i]));
    //if(deg == 2) proba *= ((double)(word[1]->getNumSupports(idxWord[1],v)) / (double)(nWords[1]));
    
    total_proba += proba;
    if(proba > best) {
      best = proba;
      val = v;
    }
  } while( valit->next() );
  return val;
}

void ValSelectorCrossword::left() { 
  getBest();
  _X->setDomain( val );
}

void ValSelectorCrossword::right() { 
  _X->remove( val ); 
}

void ValSelectorCrossword::postCut( const int p ) {
  val = p;
}

void ValSelectorCrossword::printLeft(std::ostream& o) const { 
  _X->printshort(o);
  cout << " == " << val;
}

void ValSelectorCrossword::printRight(std::ostream& o) const {
  _X->printshort(o);
  cout << " =/= " << val;
}



void snail( int xlb, int xub, int x,
	    int ylb, int yub, int y,
	    int direction,
	    vector< int >& coords,
	    int count
	    )
{
  if(count) {
    if( direction == 0 ) // up
      {      
	if( y == ylb ) {
	  --xub;
	  direction = 1; // left
	} else {
	  //coords.push_back(x);
	  //coords.push_back(y);
	  //cout << x << " " << y << " " << (y*Width+x+1) << endl;
	  coords.push_back( y*Width+x );
	  --count;
	  --y;
	}
      }
    else if( direction == 1 ) // left
      {
	if( x == xlb ) {
	  ++ylb;
	  direction = 2; // down
	} else {
// 	  coords.push_back(x);
// 	  coords.push_back(y);
// 	  cout << x << " " << y << " " << (y*Width+x+1) << endl;
	  coords.push_back( y*Width+x );
	  --count;
	  --x;
	}
      }
    else if( direction == 2 ) // down
      {
	if( y == yub ) {
	  ++xlb;
	  direction = 3; // right
	} else {
// 	  coords.push_back(x);
// 	  coords.push_back(y);
// 	  cout << x << " " << y << " " << (y*Width+x+1) << endl;
	  coords.push_back( y*Width+x );
	  --count;
	  ++y;
	}
      }
    else if( direction == 3 ) // right
      {
	if( x == xub ) {
	  --yub;
	  direction = 0; // up
	} else {
// 	  coords.push_back(x);
// 	  coords.push_back(y);
// 	  cout << x << " " << y << " " << (y*Width+x+1) << endl;
	  coords.push_back( y*Width+x );
	  --count;
	  ++x;
	}
      }
    snail( xlb, xub, x, ylb, yub, y, direction, coords, count );
  }
}





// class VarSelectorProbaOverWeight
// {
// public: 

//   /**@name Constructors*/
//   //@{
//   VarSelectorProbaOverWeight() {dom_ = NOVAL; pro_ = NOVAL; wgt_ = 0;}
//   //@}

//   /**@name Parameters*/
//   //@{ 
//   int wgt_;
//   int dom_;
//   double pro_;
//   //@}  
  
//   /**@name Utils*/
//   //@{ 
//   inline bool operator<( VarSelectorProbaOverWeight& x ) const 
//   { 
//     return (((double)dom_ + 100*pro_) * (double)(x.wgt_)) < (((double)(x.dom_) + 100*x.pro_) * (double)wgt_) ; 
//     //return pro_ * (double)(x.wgt_) < x.pro_ * (double)wgt_ ; 
//   }
//   inline void operator=( VarSelectorProbaOverWeight& x ) { wgt_ = x.wgt_; pro_ = x.pro_; dom_ = x.dom_; }
//   inline void operator=( VariableInt    *x ) { 
//     wgt_ = x->weight; 
//     dom_ = x->domsize(); 
//     pro_ = ((ValSelectorCrossword*)x->branch)->getProba();
//     //cout << pro_ << "/" << wgt_ << endl;
//   }
//   //@}  
// };

// class ProbaOverWDeg : public VariableOrdering {
//  public:

//   int size;
//   ProbaOverWDeg( const int sz=0 ) { size=sz; }

//   /**@name Utils*/
//   //@{
//   DVO* extract( Solver* s );
//   //@}
// };


// DVO* ProbaOverWDeg::extract( Solver* s )
// {
//   s->setLearner( Weighter::WDG );
//   if( size > 1 )
//     return new GenericRandomDVO<VarSelectorProbaOverWeight>(s, size);
//   else 
//     return new GenericDVO<VarSelectorProbaOverWeight>(s);
// }


