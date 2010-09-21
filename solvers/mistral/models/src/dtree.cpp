
#include <mistral_sol.h>
#include <mistral_sat.h>

#ifdef _CPLEX

#include "check.h"
#include <ilcplex/cplex.h>

#endif

using namespace std;

int nFeat;
int mFeat;
int NOTHING = 0;
int LP_BARRIER = 1;
int LP_SIMPLEX = 2;
int MIP = 3;
double eps = 1e-5;

int ENDDEPTH =0;

long int ENDCPLEXCALLS =0;
long int ENDPROPAGS =0;
long int ENDNODES =0;
long int ENDBACKTRACKS =0;
long int ENDFAILURES =0;
long int FIRSTCPLEXCALLS =-1;
long int FIRSTPROPAGS =-1;
long int FIRSTNODES =-1;
long int FIRSTBACKTRACKS =-1;
long int FIRSTFAILURES =-1;
long int OPTCPLEXCALLS =0;
long int OPTPROPAGS =-1;
long int OPTNODES =-1;
long int OPTBACKTRACKS =-1;
long int OPTFAILURES =-1;
double ENDTIME = -1;
double OPTTIME = -1;
double FIRSTTIME = -1;
double STARTTIME ;

string OutTree;
int bestTree;
int firstTree=-1;

int stopValue=0;


string input_file;

// int ntExamples;
// int npExamples;
// int nnExamples;


int lexico( const void* _x, const void* _y )
{

  //cout << (*((string*)_x)) << " <? " << (*((string*)_y)) << endl;

  const char* cx = ((string*)_x)->c_str();
  const char* cy = ((string*)_y)->c_str();
  
  float fx = atof(cx);
  float fy = atof(cy);
  
  //int res =  ((string*)_x)->compare( *((string*)_y) );

  //cout << fx << " " << fy << " ";


  int res = 0;
  if( fx < fy ) {
    res = -1;
  } else if( fx > fy ) {
    res = 1;
  } else res = ((string*)_x)->compare( *((string*)_y) );
  return res;
}



class DataSet
{
public:

  DataSet() { discrepancies = NULL; }
  virtual ~DataSet() { 
    delete [] binExamples; 
    delete [] discrepancies;
    delete [] bin2feature;
    delete [] bin2value;
    for(int i=0; i<numExamples; ++i)
      delete [] examples[i];
    delete [] examples;
  }
  
  // index of the class, given its name
  map< string,int > class2index;
  // index of the given valuation of the given feature
  vector< map<string,int> > choice2index;
  // name of the classs
  vector< string > index2class; 
  // valuations of the features
  vector< vector < string > > index2choice;
  // final examples used (all numeric data)
  //vector< vector < int > > examples;
  int** examples;
  int numExamples;
  int nFeatures;

  /////////// binary encoding
  BitSet *binExamples;
  vector< int > positive_example;
  vector< int > negative_example;
  unsigned int nBinFeatures;

  int *bin2feature;
  int *bin2value;


  int *exampleIndex;
  BitSet *discrepancies;
  vector<int> *columns;
  int nDiscrepancies;


  int getFeature( const int binFeature, string& value )
  {
    int feature = bin2feature[binFeature];
    value = index2choice[feature][bin2value[binFeature]];
    return feature;
  }  


  void printNumeric( ostream& o )
  {
    int i, j, space, unknown;

    // , maxlog = 2;
    //     for(i=0; i<nFeatures; ++i) {
    //       int nlog = (int)(log10(index2choice[i].size()));
    //       if( nlog > maxlog ) maxlog = nlog;
    //     }

    
    for(i=0; i<numExamples; ++i)
      {
	unknown = false;
	for(j=0; j<nFeatures; ++j)
	  if(index2choice[j][examples[i][j]] == "?")
	    unknown = true;

	if( unknown ) continue;
	
	space = (int)(log10(index2choice[j].size()))+1;
	if(space < 1) space = 1;
	for(j=0; j<nFeatures; ++j)
	  o << setw(space) << examples[i][j] << ", "; 
	o << setw(space) << (examples[i][j] != 0) << endl; 
	//o << setw(space) << (examples[i][j] == 2) << endl; 
      }
  }


  void binaryEncode( int classe )
  {

    int i,j,k,l;
    nBinFeatures = 0; 
    for(i=0; i<nFeatures; ++i)
      if(index2choice[i].size() > 2)
	nBinFeatures += index2choice[i].size();
      else if(index2choice[i].size() == 2)
	++nBinFeatures;

    //cout << mFeat << " " << nFeat << endl;

    nFeat *= nBinFeatures;
    mFeat *= nBinFeatures;

    //cout << mFeat << " " << nFeat << endl;

    bin2feature = new int[nBinFeatures];
    bin2value = new int[nBinFeatures];
    
    k=0;
    for(i=0; i<nFeatures; ++i) {
      if(index2choice[i].size() > 2)
	for(j=0; j<index2choice[i].size(); ++j) {
	  bin2feature[k] = i;
	  bin2value[k] = j;
	  ++k;
	}
      else if(index2choice[i].size() == 2) {
	bin2feature[k] = i;
	bin2value[k] = 0;
	++k;
      }
    }

    assert( k == nBinFeatures );


    //    cout << nBinFeatures << endl;

    binExamples = new BitSet[numExamples];
    for(i=0; i<numExamples; ++i)
      {

	//	cout << i << ": " << endl;

	//cout << mFeat << " " << nFeat << endl;


	binExamples[i].init(mFeat, nFeat-1, BitSet::empt);

	k=0;
	for(j=0; j<nFeatures; ++j) {
	  if(index2choice[j].size() > 2) {
	    
	    //	    cout << "\t" << (k+examples[i][j]) << endl;

	    if(!mFeat)
	      binExamples[i].insert( k+examples[i][j] );
	    if(nFeat > nBinFeatures)
	      binExamples[i].addInterval( nBinFeatures+k+examples[i][j], nBinFeatures+k+index2choice[j].size()-1 );

	    k += index2choice[j].size();
	  } else if(index2choice[j].size() == 2) {
	    //	    cout << examples[i][j] << " ";
	    if(!examples[i][j]) {
	      if(nFeat > nBinFeatures)
		binExamples[i].insert(nBinFeatures+k);
	      if(!mFeat)
		binExamples[i].insert(k);
	    }
	    ++k;
	  }
	}
	if(examples[i][nFeatures] == classe) 
	  //positive_example.push_back( binExamples[i] );
	  positive_example.push_back( i );
	else 
	  //negative_example.push_back( binExamples[i] );
	  negative_example.push_back( i );
      }


    for(i=0; i<numExamples; ++i)
      for(j=0; j<nBinFeatures; ++j) {
	if(!mFeat)
	  assert( binExamples[i][j] == (examples[i][bin2feature[j]] == bin2value[j]) );
	if(nFeat > nBinFeatures)
	  assert( binExamples[i][nBinFeatures+j] == (examples[i][bin2feature[j]] <= bin2value[j]) );
      }

//     npExamples = positive_example.size();
//     nnExamples = negative_example.size();

//     cout << "positive_example: " << positive_example.size() << endl;
//     cout << "negative_example: " << negative_example.size() << endl;
//     //exit(0);
    
  }


  void computeXor( )
  {

    //    MinHittingSet mhs(nBinFeatures, negative_example, positive_example);
    int incl, k, i, j, l, n=positive_example.size(), m=negative_example.size();
    discrepancies = new BitSet[n*m];
    exampleIndex = new int[n*m];
    columns = new vector<int>[nFeat-mFeat];
    for(k=n*m; k;)       
      discrepancies[--k].init(mFeat, nFeat-1, BitSet::empt);
    for(i=0; i<n; ++i)
      for(j=0; j<m; ++j) {
	discrepancies[k].clear();
	discrepancies[k].unionWith(binExamples[positive_example[i]]);
	discrepancies[k].xorWith(binExamples[negative_example[j]]);
	// 	discrepancies[k].unionWith(positive_example[i]);
	// 	discrepancies[k].xorWith(negative_example[j]);
	exampleIndex[k] = i;

	// 	incl = false;
	// 	for(l=k; !incl && l;) {
	// 	  if( discrepancies[--l].included( discrepancies[k] ) ) incl = true;
	// 	} if( incl ) {
	// 	  --k;
	// 	  continue;
	// 	} else {
	// 	  incl = 0;
	// 	  for(l=k; l;) {
	// 	    if( discrepancies[k].included( discrepancies[--l] ) ) {
	// 	      discrepancies[l].swap(discrepancies[k-(++incl)]);
	// 	int aux = exampleIndex[l];
	// 	exampleIndex[l] = exampleIndex[k-incl];
	// 	exampleIndex[k-incl] = aux;
	//      }
	// 	  }
	// 	  if( incl ) {    
	// 	    discrepancies[k-incl].swap( discrepancies[k] );
	// 	    k -= incl;
	// 	  }
	// 	}

	++k;
      }
    

    nDiscrepancies = k;
    for(i=0; i<k; ++i) {
      //       cout << setw(6) << i << " ";
      //       discrepancies[i].print( cout );
      //       cout << endl;

      BitsetIterator bit( discrepancies[i] );
      do {
	j = bit-mFeat;
	columns[j].push_back(i);
      } while( ++bit );


    }
    //    cout << endl;


  }


  void getData( const char *filename )
  {
    int i,j,k,val;//n=0;
    numExamples = 0;
    char buf[1024];

    ifstream datafile( filename );  
    assert( datafile.good() );

    nFeatures = 0;
    vector<string> data;
  
    while( true ) {
      datafile.getline( buf, 1024 );
      if( datafile.eof() || !datafile.good() ) break;

//       // new example
//       vector< int > ex;
//       examples.push_back( ex );

      j=0;
      k=0;
      for( i=0; buf[i] != '\0'; ++i ) {
	if( buf[i] == ',' ) { 
	  // new feature?
	  if( nFeatures <= k ) {
	    ++nFeatures;
	    map<string, int>  valuations_map;
	    vector<string> valuations_vector;
	    choice2index.push_back( valuations_map );
	    index2choice.push_back( valuations_vector );
	  }
	 
	  buf[i] = '\0';
	  string valuation(&buf[j]);

	  data.push_back( valuation );

	  // new valuation?
	  if( choice2index[k].find( valuation ) == choice2index[k].end() ) {
	    choice2index[k][valuation] = index2choice[k].size();
	    index2choice[k].push_back( valuation );
	  }
	  // build example
	  //examples.back().push_back( choice2index[k][valuation] );

	  j=i+1;	
	  ++k;
	}
      }	
      
      if(j) {
	++numExamples;
	string theclass(&buf[j]);

	data.push_back( theclass );
	
	// new class?
	if( class2index.find( theclass ) == class2index.end() ) {
	  class2index[theclass] = index2class.size();
	  index2class.push_back( theclass );
	}
	
	//examples.back().push_back( class2index[theclass] );
	
      } 
//       else {
// 	examples.pop_back();
//       }
    }
    
//     for(i=0; i<nFeatures; ++i) {
//       for(j=0; j<index2choice[i].size(); ++j)
// 	cout << setw(3) << index2choice[i][j] ;
//       cout << endl;
//     }

    nFeatures = index2choice.size();

    for(i=0; i<nFeatures; ++i) {
      qsort( &(index2choice[i][0]), index2choice[i].size(), sizeof(const char*), lexico );
      for(j=0; j<index2choice[i].size(); ++j) 
	choice2index[i][index2choice[i][j]] = j;
    }

//     for(i=0; i<nFeatures; ++i) {
//       for(j=0; j<index2choice[i].size(); ++j)
// 	cout << setw(3) << index2choice[i][j] ;
//       cout << endl;
//     }
    
    //exit(0);

    //ntExamples = numExamples;
    //cout << "numExamples: " << numExamples << endl;

    examples = new int*[numExamples];
       
    k=0;
    for(i=0; i<numExamples; ++i) {
      // new example
      examples[i] = new int[nFeatures+1];
      for(j=0; j<nFeatures; ++j)
	examples[i][j] = choice2index[j][data[k++]];
      examples[i][nFeatures] = class2index[data[k++]];
    }
    
//     for(i=0; i<numExamples; ++i)
//       {
// 	for(j=0; j<nFeatures; ++j) 
// 	  cout << setw(3) << index2choice[j][examples[i][j]];
// 	cout << setw(3) << index2class[examples[i][j]] << endl;
//       }
    
//     exit(0);
    
  }



  void printData( )
  {
    cout << endl << "print" << endl;

    int i, j;
    cout << "classes: " << index2class.size() ;
    for(i=0; i<index2class.size(); ++i)
      cout << " " << index2class[i];
    cout << endl;

    cout << "features: " << nFeatures << endl;
  
    for(i=0; i<nFeatures; ++i)
      {
	cout << " feature_" << i << ": ("
	     << index2choice[i].size() << ") ";
	for(j=0; j<index2choice[i].size(); ++j)
	  {
	    cout << " " << index2choice[i][j];
	  }
	cout << endl;
      }
    cout << endl;

    cout << "examples: " << numExamples << endl;
    for(i=0; i<numExamples; ++i)
      {
	for(j=0; j<nFeatures; ++j)
	  {
	    cout << index2choice[j][examples[i][j]] << " ";
	  }
	cout << " -> " << index2class[examples[i][nFeatures]] << endl;
      }
    cout << endl;

    
    cout << endl << setw(6) << negative_example.size() << " "
	 << "negative examples:" << endl;
    for( int i=0; i<positive_example.size(); ++i ) {
      cout << setw(6) << i << " ";
      binExamples[negative_example[i]].print( cout );
      cout << endl;
    }
    
    cout << endl << setw(6) << positive_example.size() << " "
	 << "positive examples:" << endl;
    for( int i=0; i<positive_example.size(); ++i ) {
      cout << setw(6) << i << " ";
      binExamples[positive_example[i]].print( cout );
      cout << endl;
    }
    cout << endl;
  }
 

  void satEncode( const int depth, SatSolver *dll=NULL )
  {
    

    string output_file = input_file+".cnf";
    //cout << output_file << endl;

    ofstream outfile( output_file.c_str(), ios_base::out );


    int i, j, k, l, node, incr, nLeaves = (1 << depth);
    int nNodes = nLeaves-1;
    int **path = new int*[nLeaves];
    for(i=0; i<nLeaves; ++i)
      path[i] = new int[depth+1];

    /// pre-compute all paths
    incr = nLeaves;
    node = 0;
    for(i=0; i<=depth; ++i)
      {
	for(k=0; k<nLeaves; k+=incr) {
	  for(j=k; j<k+incr; ++j)
	    path[j][i] = node;
	  ++node;
	}
	incr /= 2;
      }

 //    for(i=0; i<nLeaves; ++i) {
//       for(j=0; j<=depth; ++j)
// 	cout << path[i][j] << " ";
//       cout << endl;
//     }


    ///header
    int nClauses = //nNodes + 
      (nBinFeatures * (nBinFeatures-1) / 2) * nNodes + 
      nLeaves * nDiscrepancies 
      + nLeaves * (depth-1) * nBinFeatures;
      
    if( dll )
      dll->init((nNodes*nBinFeatures), nClauses);
    else
      outfile<< "p cnf " << (nNodes*nBinFeatures) << " " << nClauses << endl; 

    Vector<Literal> clause;

    /// domain clauses

//     long long int size = ((long long int)nNodes * (long long int)nBinFeatures * (long long int)(nBinFeatures-1));
//     size += ((long long int)nBinFeatures * (long long int)nLeaves * (long long int)nDiscrepancies * (long long int)depth);
//     size += ((long long int)nBinFeatures * (long long int)nBinFeatures * (long long int)nLeaves * (long long int)depth);


//     cerr << "c " << depth << " "  << nNodes << " " << nBinFeatures << " " << size << endl;

//     exit( 0 );

    for(i=0; i<nNodes; ++i)
      {
// 	for(j=0; j<nBinFeatures; ++j)
// 	  {

// 	    if(dll) clause.push( i*nBinFeatures+j+1 );
// 	    else {		      
// 	      cout << i*nBinFeatures+j+1 << " ";
// 	      string value;
// 	      int feature = getFeature(j, value);
// 	      cerr << i*nBinFeatures+j+1 << " " << i 
// 		//<< " example[" 
// 		   << " " 
// 		   << feature 
// 		//<< "] == " 
// 		   << " " 
// 		   << value << endl;
// 	    }

// 	  }	       

// 	if(dll) {
// 	  dll->addClause( dll->base, clause, dll->stats.base_avg_size );
// 	  dll->addOriginalClause( clause );
// 	  clause.clear();
// 	} else cout << "0\n";

// 	--nClauses;
	for(j=1; j<nBinFeatures; ++j)
	  for(k=0; k<j; ++k) {

	    if(dll) {
	      clause.push(((int)(-(i*nBinFeatures+j+1))));
	      clause.push(((int)(-(i*nBinFeatures+k+1))));
	      dll->addClause( dll->base, clause, dll->stats.base_avg_size );
	      dll->addOriginalClause( clause );
	      clause.clear();
	    } else {
	      outfile<< ((int)(-(i*nBinFeatures+j+1))) << " " 
		   << ((int)(-(i*nBinFeatures+k+1))) << " 0\n";
	    }

	    --nClauses;
	  }
      }

    for(i=0; i<nDiscrepancies; ++i) {

//       cout << " ---- " ;
//       discrepancies[i].print(cout);
//       cout << endl;

      for(j=0; j<nLeaves; ++j) {

	//	cout << "leaf " << j << endl;

	for(k=0; k<nBinFeatures; ++k) {

	  //	  cout << i << " " << j << " " << k << " ";

	  bool goToLeft = !(binExamples[exampleIndex[i]].member(k));

// 	  cout << endl;
// 	  binExamples[exampleIndex[i]].print( cout );
// 	  cout << endl;
// 	  discrepancies[i].print(cout);
// 	  cout << endl;
// 	  cout << k << " goes to " << (goToLeft ? "left" : "right") << endl;

	  if(discrepancies[i].member(k))
	    for(l=0; l<depth; ++l) {
	      //cout << "\t";
	      if(dll) clause.push((path[j][l]*nBinFeatures+k+1));
	      else outfile<< (path[j][l]*nBinFeatures+k+1) << " ";
	      //cout << " diff" << endl;
	    }
	  else
	    for(l=0; l<depth; ++l) {
// 	      cout << "\t";
// 	      cout << path[j][l] << " " << (depth-l-1) << " ";

	      int largestLeftLeaf = path[j][l]*(1 << (depth-l))+3*(1 << (depth-l-1))-2;

// 	      cout << largestLeftLeaf << ": " << j << " is on the " 
// 		   << ((nNodes+j) <= largestLeftLeaf ? "left" : "right") << endl;
	      
	      if( ((nNodes+j) <= largestLeftLeaf) != goToLeft )
		if(dll) clause.push((path[j][l]*nBinFeatures+k+1));
		else outfile<< (path[j][l]*nBinFeatures+k+1) << " ";
	    }
	}

	if(dll) {
	  dll->addClause( dll->base, clause, dll->stats.base_avg_size );
	  dll->addOriginalClause( clause );
	  clause.clear();
	} else outfile<< "0\n";

	--nClauses;
      }
    }


    //cout << nClauses << endl
    //<< ()

    // avoid internal nodes with no tests
    for(i=0; i<nLeaves; ++i) {
      for(j=0; j<depth-1; ++j) {
	//cout << path[i][j] << " ";
	
 	for(k=0; k<nBinFeatures; ++k) {
	 
	  //if( !dll )
	  // cout << "== ";
	
	  for(l=0; l<nBinFeatures; ++l) 
	    if( k != l ) {

	      if(dll) {
		clause.push(((int)(path[i][j]*nBinFeatures+l+1)));
	      } else {
		outfile<< ((int)(path[i][j]*nBinFeatures+l+1)) << " " ;
	      }
	      
	    }

	  if(dll) {
	    clause.push(((int)(-(path[i][j+1]*nBinFeatures+k+1))));
	    dll->addClause( dll->base, clause, dll->stats.base_avg_size );
	    dll->addOriginalClause( clause );
	    clause.clear();
	  } else {
	    outfile<< ((int)(-(path[i][j+1]*nBinFeatures+k+1))) << " 0\n";
	  }

	  --nClauses;

	}
       
      }
      //cout << endl;
    }
    
    assert ( !nClauses );

  }


};





class ExampleStruct;
// M subsets of U, with |U| = N
class MinHittingSet {

 public:
  /**@name Parameters*/
  //@{ 
  int nColumns;
  int nRows;

  double epsilon;

  int status;
  char     *probname;  
  int      objsen;
  double   *obj;
  double   *rhs;
  char     *sense;
  int      *matbeg;
  int      *matcnt;
  int      *matind;
  double   *matval;
  double   *zlb;
  double   *zub;
  char     *ctype;
  
  int      solstat;
  int      lpstat;
  double   objval;
  double   *x;//[NUMCOLS:N];
  double   *pi;//[NUMROWS:M];
  double   *slack;//[NUMROWS:M];
  double   *dj;//[NUMCOLS:N];
  
#ifdef _CPLEX

   CPXENVptr     the_env;
   CPXLPptr      the_lp;

#endif

  //@}

  MinHittingSet(DataSet& ds);
  MinHittingSet(ExampleStruct *es);
  virtual ~MinHittingSet();

  void init(const int nFeatures, 
	    vector<int> *columns,
	    const int nDiscrepancies);

  void print();
  double getSimplexBound();
  double getBarrierBound();
  double solve();

};




class ExampleStruct {

public:

  static const int left=0;
  static const int right=1;
  static const int positive=1;
  static const int negative=0;
  
  ReversibleIntList  *Example[2];
  
  VariableInt**     Feature;
  unsigned int       nNodes;
  unsigned int nBinFeatures;
  unsigned int    nFeatures;
  unsigned int nExamples[2];
  BitSet*      bin_examples;
  int**        int_examples;
  int*      idx_examples[2];
  int*          bin2feature;
  int*            bin2value;
  string**     index2choice;


  /********************************************************************
   * data structure that keeps, for every node, the number of positive 
   * and negative examples that would be on the left or right handside
   ********************************************************************/

  // for each node 'n', feature 'f', polarity 'p' and direction 'd'
  // how many examples of polarity 'p' would be routed to the direction 'd' 
  // at node 'n' if it was assigned with feature 'f' 
  int **remainder[2][2];

  int xorThreshold;
  int Mode;
  BitSet *discrepancies;
  int *exampleIndex;
  vector<int> *columns;
  int nDiscrepancies;


public:

  int lowerBound1( const int node )
  {
    int n=Example[positive][node].size, m=Example[negative][node].size;

    //assert( n && m ) ;

    if( n && m ) {
      
      if( !(Feature[node]->isGround()) && (n*m <= xorThreshold) ) {
	double result;
	computeXor( node, n, m );
	MinHittingSet *mhs = new MinHittingSet(this);
	
	//mhs->print();
	
	if(Mode == MIP)
	  result = mhs->solve();
	else if(Mode == LP_BARRIER)
	  result = mhs->getBarrierBound();
	else if(Mode == LP_SIMPLEX)
	  result = mhs->getSimplexBound();
	else {
	  cerr << "unknown mode!" << endl;
	  exit(0);
	}
	
	delete mhs;
	
	/// HACK: it seems that result could be slightly larger than the maximum value, 
	///       so we remove a small bit before ceiling it
	return (int)(ceil(result-eps));
	
      } else {

	return 1;
	
      }
    } else {

      return 0;

    }
  }


  void computeXor( const int node, const int n, const int m )
  {
    //    MinHittingSet mhs(nBinFeatures, negative_example, positive_example);
    int i, j, k=0, i_idx, j_idx;//, n=Example[positive][node].size, m=Example[negative][node].size;
    
    for(i=0; i<(nFeat-mFeat); ++i)
      columns[i].clear();
    
    for(i=0; i<n; ++i) {
      i_idx = Example[positive][node][i];
      for(j=0; j<m; ++j) {
	j_idx = Example[negative][node][j];
	discrepancies[k].clear();
	discrepancies[k].unionWith(bin_examples[idx_examples[positive][i_idx]]);
	discrepancies[k].xorWith(bin_examples[idx_examples[negative][j_idx]]);
	exampleIndex[k] = i_idx;
	++k;
      } 
    }
    nDiscrepancies = k;
    for(i=0; i<k; ++i) {
      for(j=mFeat; j<nFeat; ++j)
	if( discrepancies[i].member(j) )
	  columns[j-mFeat].push_back(i);

//       discrepancies[i].print( cout );
//       cout << endl;
//       BitsetIterator bit( discrepancies[i] );
//       int lastbit = -1;
//       do {
// 	j = bit-mFeat;

// 	cout << (j+mFeat) << " ";
// 	cout.flush();

// 	assert( j != lastbit );
// 	lastbit = j;

// 	//cout << mFeat << " " << j << " " << (nFeat-mFeat) << " " << i << endl;
// 	columns[j].push_back(i);
//       } while( ++bit );
//       cout << endl;
    }

//     for(i=0; i<nFeatures; ++i)
//       {
// 	for(j=0; j<columns[i].size(); ++j)
// 	  cout << " " << columns[i][j];
// 	cout << endl;
//       }

    //exit(0);
  }


  void checkExamples(const int node, const int id)
  {
    int numEx = (nExamples[positive] + nExamples[negative]);
    for(int polarity=negative; polarity<=positive; ++polarity)
      {
	int size = Example[polarity][node].size;
	for(int i=0; i<size; ++i) {
	  if( Example[polarity][node][i] < 0 )
	    cout << id << endl;
	  if( Example[polarity][node][i] > numEx )
	    cout << id << endl;
	  assert( Example[polarity][node][i] >= 0 );
	  assert( Example[polarity][node][i] <= numEx );
	}
      }
  }


  void init( const int t )
  {
    xorThreshold = t;
    discrepancies = new BitSet[xorThreshold];

    for(int k=0; k<xorThreshold; ++k)       
      discrepancies[k].init(mFeat, nFeat, BitSet::empt);

    exampleIndex = new int[xorThreshold];
    columns = new vector<int>[nFeat-mFeat];
  }


  /**@name Constructors*/
  //@{
  //ExampleStruct(const int ln, int **pe, int np, int **ne, int nn, int nf)
  //ExampleStruct(const int ln, const int nf, vector< BitSet >& pe, vector< BitSet >& ne)
  ExampleStruct(DataSet& ds, const int n)
  {
    xorThreshold = -1;

    nExamples[positive] = ds.positive_example.size(); 
    nExamples[negative] = ds.negative_example.size();
    idx_examples[positive] = &(ds.positive_example[0]); 
    idx_examples[negative] = &(ds.negative_example[0]);

    bin2value = ds.bin2value;
    bin2feature = ds.bin2feature;
    nBinFeatures = ds.nBinFeatures;
    nFeatures = ds.nFeatures;
    nNodes = n;
    Example[positive] = new ReversibleIntList[nNodes];
    Example[negative] = new ReversibleIntList[nNodes];

   
    bin_examples = ds.binExamples;
    int_examples = ds.examples;

// new int*[ds.numExamples];
//     int i=ds.numExamples;
//     while( i-- )
//       int_examples[i] = &(ds.examples[i][0]);
    index2choice = new string*[ds.nFeatures];
    int i=ds.nFeatures;
    while( i-- )
      index2choice[i] = &(ds.index2choice[i][0]);
//     for( int polarity = negative; polarity <= positive; ++polarity ) {
//       cout << endl << setw(6) << examples[polarity].size() << " "
// 	   << (polarity ? "positive examples:" : "negative examples:") << endl;
//       for( int i=0; i<examples[polarity].size(); ++i ) {
// 	cout << setw(6) << i << " ";
// 	examples[polarity][i].print( cout );
// 	cout << endl;
//       }
//       cout << endl;
//     }

    for( int polarity = negative; polarity <= positive; ++polarity )
      for( int direction = left; direction <= right; ++direction ) {
	remainder[polarity][direction] = new int*[nNodes];
	for( int node = 0; node < nNodes; ++node ) {
	  remainder[polarity][direction][node] = new int[nFeat-mFeat];
	  remainder[polarity][direction][node] -= mFeat;
	  std::fill( remainder[polarity][direction][node]+mFeat, remainder[polarity][direction][node]+nFeat, 0 );
	}
      }

    discrepancies = NULL;
    exampleIndex = NULL;
    columns = NULL;
  }

  void init(Solver *s, VariableInt **feat) {
    Feature = feat;
    s->binds( Example[positive][0] ); 
    s->binds( Example[negative][0] ); 
    Example[positive][0].setValue( 0, nExamples[positive]-1, nExamples[positive] );
    Example[negative][0].setValue( 0, nExamples[negative]-1, nExamples[negative] );
    checkExamples(0,0);
    for(int i=1; i<nNodes; ++i) {
      s->binds( Example[positive][i] ); 
      s->binds( Example[negative][i] ); 
      Example[positive][i].setValue( 0, nExamples[positive]-1, 0 );
      Example[negative][i].setValue( 0, nExamples[negative]-1, 0 );
      checkExamples(i,0);
    }

    
  }

  virtual ~ExampleStruct()
  {
    delete [] discrepancies;
    delete [] exampleIndex;
    delete [] columns;
    delete [] Feature;
    delete [] index2choice;
    delete [] Example[negative];
    delete [] Example[positive];
    for( int polarity = negative; polarity <= positive; ++polarity ) {
      for( int direction = left; direction <= right; ++direction ) {
	for( int node = 0; node < nNodes; ++node ) {
	  remainder[polarity][direction][node] += mFeat;
	  delete [] remainder[polarity][direction][node];
	}
	delete [] remainder[polarity][direction];
      }
    }
  }
  //@}

  void printInfo( int node, bool fullq=false ) 
  {
    cout << "n" << node << ": "
	 << (Example[positive][node].size) << "/"
	 << (Example[negative][node].size) << endl;
    if( fullq ) {
      for( int feature = 0; feature < nFeat; ++feature ) 
	if( Feature[node]->contain( feature ) ) {
	  cout << "\tleft: " << (feature) << " (" << (bin2feature[feature%nBinFeatures]+1)
	       << ( feature < nBinFeatures ? " = " : " <= " ) 
	       << index2choice[bin2feature[feature%nBinFeatures]][bin2value[feature%nBinFeatures]] << ") : "
	       << (remainder[positive][left][node][feature]) << "/"
	       << (remainder[negative][left][node][feature]) 
	       << "  .   right: "
	       << (remainder[positive][right][node][feature]) << "/"
	       << (remainder[negative][right][node][feature])  << " => "
	       << informationGain(node, feature) << " "
	       << (isRelevant(node, feature) ? (decides(node, feature) ? "*" : " " ) : "irrelevant" ) 
	       << endl;	
	} 
      cout << endl;
    }
  }

  void printExamples( int node ) 
  {
    cout << endl;
    for( int polarity=negative; polarity<=positive; ++polarity ) {
      cout << (polarity ? "positive:" : "negative:" ) << endl;
      for( int i=0; i<Example[polarity][node].size; ++i ) {
	for( int j=0; j<nFeatures; ++j ) {
	  cout << " ";
	  int value = (int_examples[idx_examples[polarity][Example[polarity][node][i]]][j]);
	  cout << index2choice[j][value];
	}
	cout << " " << (int_examples[idx_examples[polarity][Example[polarity][node][i]]][nFeatures]) << endl;
      }
      cout << endl;
    }
  }





  inline int needChild( const int direction, const int node, const int feature )
  {
    
    assert( feature < nFeat );

    return ( remainder[negative][direction][node][feature] && remainder[positive][direction][node][feature] );
  }
  
  inline int isRelevant( const int node, const int feature )
  {

    //cout << "isRelevant(" << node << ", " << feature << ") - [" << nNodes << ", " << nFeatures << "]" << endl;

    return ( (remainder[negative][left][node][feature] +remainder[positive][left][node][feature]) && 
	     (remainder[negative][right][node][feature]+remainder[positive][right][node][feature]) );
  }

  inline int decides( const int node, const int feature )
  {
    return ( !((remainder[negative][left][node][feature]*remainder[positive][left][node][feature]) + 
	       (remainder[negative][right][node][feature]*remainder[positive][right][node][feature])) );
  }
  
  inline double entropy( const int direction, const int node, const int feature ) // compute entropy(S_direction) for the given feature
  {
    //    assert(feature < nFeat);
    double npos = remainder[positive][direction][node][feature];
    double nneg = remainder[negative][direction][node][feature];

//     if(remainder[positive][direction][node][feature] > 1000000)
//       cout << direction << " " << node << " " << feature << " " 
// 	   << remainder[positive][direction][node][feature] << " "
// 	   << remainder[negative][direction][node][feature] << " " << endl;

    double tot = npos+nneg;
    double etp = 0;
    if( npos > 0 && nneg > 0 ) 
      etp = -(((npos/tot) * log2(npos/tot)) + ((nneg/tot) * log2(nneg/tot))); 
    return etp;
  }

  inline double entropy( const int node ) // compute the entropy entropy(S)
  {
    double npos = Example[positive][node].size, nneg = Example[negative][node].size, tot=npos+nneg;
    double etp = 0;
    if( npos > 0 && nneg > 0 ) 
      etp = -(((npos/tot) * log2(npos/tot)) + ((nneg/tot) * log2(nneg/tot)));
    return etp;
  }

  inline double informationGain( const int node, const int feature ) // Gain(S, feature)
  {
    double S = Example[positive][node].size + Example[negative][node].size;
    double Sl = remainder[positive][left][node][feature] + remainder[negative][left][node][feature];
    double Sr = remainder[positive][right][node][feature] + remainder[negative][right][node][feature];

    Sl*=entropy(left,node,feature);
    Sr*=entropy(right,node,feature);
    Sl+=Sr;
    Sl/=S;

    return (entropy(node) - Sl);
  }

  void updateEntropy( const int node )
  {
    int i, j, k, n;//, nchild;
    int direction, polarity;
    
    for( polarity = negative; polarity <= positive; ++polarity )
      for( direction = left; direction <= right; ++direction )
	std::fill( remainder[polarity][direction][node]+mFeat, remainder[polarity][direction][node]+nFeat, 0 );

    int nc, val, feat;
    DomainIterator *feature = Feature[node]->begin();
    do {
      k = *feature;
      if(k<nBinFeatures) {
	for( polarity = negative; polarity <= positive; ++polarity ) {
	  n = Example[polarity][node].size;
	  for(i=0; i<n; ++i) {	  
	    j = Example[polarity][node][i];
	    ++remainder[polarity][bin_examples[idx_examples[polarity][j]][k]][node][k];	  
	  }
	}
      } else if(k<nFeat) {

	assert( (k-nBinFeatures) < nBinFeatures );
	assert( (k-nBinFeatures) >= 0 );

	feat = bin2feature[k-nBinFeatures];
	val = bin2value[k-nBinFeatures];

// 	if( feat )
// 	  {
// 	    k++;
// 	    --k;
// 	  }


	assert( feat >= 0 );
	assert( feat < nFeatures );
	//cout << feat << endl;

	for( polarity = negative; polarity <= positive; ++polarity ) {
	  n = Example[polarity][node].size;
	  for(i=0; i<n; ++i) {
	    j = Example[polarity][node][i];

// 	    if( polarity == positive )
// 	      assert( j < npExamples );
// 	    else
// 	      assert( j < nnExamples );

// 	    assert( idx_examples[polarity][j] < ntExamples );

// 	    //cout << k << " " << feat << " " << j << " " << idx_examples[polarity][j]

	    int x1 = idx_examples[polarity][j];

// 	    if( int_examples[x1] > int_examples[5] )
// 	      {
// 		x1 += 1;
// 		--x1;
// 	      }	    


// 	    if( x1 == 1 )
// 	      cerr << x1 << " " << ntExamples << " " << feat << " " << nFeatures << endl;
	    if( int_examples[x1][feat] > val ) {

// 	    //if( int_examples[idx_examples[polarity][j]][feat] > val ) {
// 	      //cerr << node << " < " << nNodes << "  " << mFeat << " <= " << k << " < " << nFeat << endl;
// 	      //cerr << (remainder[polarity][left][node] + k) << endl;
// 	      assert(k >= mFeat);
// 	      assert(k < nFeat);
	      ++remainder[polarity][left][node][k];
	    } else
	      ++remainder[polarity][right][node][k];
	  }
	}
      }
    } while( feature->next() );
  }
  
};


// This constraint does the following:
// 1 - it keeps a list of examples for each node.
// 2 - it prunes the left and right childs of a node according to the projected entropy information
class ConstraintLearning : public Constraint {

public:

static const int left=0;
static const int right=1;
static const int positive=1;
static const int negative=0;

private:
  
  unsigned int         nFeatures; 
  unsigned int            nNodes; 
  VariableInt       **   Feature;
  VariableInt       **    Parent;
  VariableInt       **  Child[2];
  VariableInt       **  NumChild;

  ExampleStruct *ES;

public:
  

  // returns the left (direction = left) or right (direction = right) child of a node
  // if it is known, and -1 otherwise.
  // notice that if the method returns 'node', it means that node
  // cannot have a child for this direction
  inline int child( const int node, const int direction )
  {
    int i = (1+direction)*arity+node;
    if( scope[i]->isGround() ) return scope[i]->min();
    return -1;
  }


  /**@name Constructors*/
  //@{
  ConstraintLearning(Solver *s, VariableInt **x,// const int ln, 
		     ExampleStruct *es)
    : Constraint(s, x, ((es->nNodes)*5), Constraint::VALUETRIGGER)
  {
    nFeatures = es->nBinFeatures;

    //assert( nf == es->nBinFeatures );

    nNodes = es->nNodes;
    Parent = scope;
    Child[left]  = scope+nNodes;
    Child[right] = scope+2*nNodes;
    NumChild     = scope+3*nNodes;
    Feature      = scope+4*nNodes;
    ES = es;

    for(int i=0; i<nNodes; ++i) {
      int j = 3*nNodes+i;
      elements[j].erase();
      triggers[j] = scope[j]->triggerOnDomain();
      triggers[j]->insert( elements[j] );
    }
  }

  virtual ~ConstraintLearning()
  {
  }
  //@}

  inline bool propagate()
  {

#ifdef _DEBUGSEARCH_
    cout << "beg propagate" << endl;
#endif

    bool consistent = true;
    for(int i=0; consistent && i<arity; ++i) {
      if( scope[i]->isGround() ) {
	if( i<nNodes && scope[i]->value() == i ) {
	  ES->updateEntropy(i);
	  consistent = ( pruneFeature( i ) );
	} else if( i>=4*nNodes )
	  consistent = propagate(i, VALUETRIGGER);
      }
    }

#ifdef _DEBUGSEARCH_
    cout << "end propagate" << endl;
#endif

    return consistent;
  } 

  /// precond: examples are known
  inline bool pruneFeature(const int node)
  {

#ifdef _DEBUGSEARCH_    
    cout << "\tprune features for " << node << endl;
#endif

    bool consistent = Feature[node]->remove(nFeat);
    int feature, direction, nc=NumChild[node]->max(), has[2], no[2], need[3] = {1,1,1};
    no[left] = Child[left][node]->equal(node);
    no[right] = Child[right][node]->equal(node);
    DomainIterator *valit = Feature[node]->begin();
    do {
      feature = *valit;

      if( !ES->isRelevant(node, feature) ) {
	consistent = Feature[node]->remove( feature );
      } else if( ES->decides(node, feature) ) {
	consistent = ( Feature[node]->setDomain( feature ) &&
		       NumChild[node]->setDomain( 0 ) &&
		       Child[left][node]->setDomain( node ) &&
		       Child[right][node]->setDomain( node )
		       );
	need[left] = need[right] = need[2] = 0;
	break;
      } else {
	for(direction=left; direction<=right; ++direction) {
	  has[direction] = ES->needChild(direction, node, feature);
	  if(no[direction] && has[direction]) {
	    consistent = Feature[node]->remove(feature);
	  }
	}

	if( (nc < 2) && has[left] && has[right] ) {
	  consistent = Feature[node]->remove(feature);
	}
	if( consistent )
	  {
	    need[left] &= has[left];
	    need[right] &= has[right];
	    need[2] &= (has[left] | has[right]);
	  }
      }
    } while( consistent && valit->next() );

#ifdef _DEBUGSEARCH_
    cout << "\t\t["
	 << Feature[node]->min()
	 << ".." << Feature[node]->max() << "]";
      
    cout << endl;
#endif

    if(consistent && need[2]) {
#ifdef _DEBUGSEARCH_
      cout << "\tpre-open new nodes for " << node << endl;
#endif
      
#ifdef _DEBUGSEARCH_
      if( need[left] && need[right] )
	cout << "\t\t pre-open " << 2 << " nodes" << endl;
      else
	cout << "\t\t pre-open " << need[2] << " nodes" << endl;
#endif
      
      if( need[left] && need[right] )
	consistent = NumChild[node]->setMin( 2 );
      else
	consistent = NumChild[node]->setMin( need[2] );
      
#ifdef _DEBUGSEARCH_
      cout << "\t\t NumChild[" << node << "] = ";
      NumChild[node]->print( cout );
      cout << endl;      
      if( !consistent )
	cout << "\t\t not enough nodes!" << endl;
#endif
      
    }

#ifdef _DEBUGSEARCH_
    cout << "\t" << (consistent ? "ok" : "inconsistent!") << endl;
#endif

    return consistent;
  }

  /// precond: the feature is known
  inline bool openNewNodes(const int node, const int feature)
  {
#ifdef _DEBUGSEARCH_
    cout << "open new nodes for " << node << " = " << feature << endl;
#endif

    int i, consistent = true, direction, child[3];
    child[2] = node;
    child[left] = -1;
    child[right] = -1;

    for(direction=right; consistent && direction>=left; --direction) {

#ifdef _DEBUGSEARCH_
      cout << endl << "\t";
      Child[direction][node]->print(cout);
      cout << endl;
#endif

      if(ES->needChild(direction, node, feature)) {

#ifdef _DEBUGSEARCH_
	cout << "\tneed a " << (direction == left ? "left" : "right") << " child" << endl;
#endif

	for(i=2; consistent && i>direction && child[i]>=0; --i)	 	
	  consistent = Child[direction][node]->remove(child[i]);

#ifdef _DEBUGSEARCH_
	cout << "\tremove non-children: (" 
	     << child[0] << " "
      	     << child[1] << " "
	     << child[2] << ") ";
	Child[direction][node]->print(cout);
	cout << endl;
#endif

	if(consistent) {
	  child[direction] = Child[direction][node]->min();

#ifdef _DEBUGSEARCH_
	  cout << "\tselect " << child[direction] << endl;
#endif

	  consistent = Child[direction][node]->setDomain(child[direction]);
	}
      } 
      
#ifdef _DEBUGSEARCH_
      cout << "\t" ;
      Child[direction][node]->print(cout);
      cout << endl;
#endif
      
    }

    //exit(0);
#ifdef _DEBUGSEARCH_
    cout << (consistent ? "ok" : "inconsistent!") << endl;
#endif

    return consistent;
  }

  /// precond: the feature is known, the childs are known
  inline bool dispatchExamples(const int node, const int feature)
  {
#ifdef _DEBUGSEARCH_
    cout << "dispatch examples for " << node << " = " << feature << endl;
#endif


    int consistent = true, direction, polarity, i, j, k, n, child[2] = {-1,-1};
    for(direction=left; consistent && direction<=right; ++direction)
      if(Child[direction][node]->isGround()) {
	child[direction] = Child[direction][node]->value();		
	if(child[direction] == node)
	  child[direction] = -1;
#ifdef _DEBUGSEARCH_
	cout << "\tchild[" << (direction == left ? "left" : "right") << "] = " << child[direction] << endl;
#endif
      }

    if(child[left] >= 0 || child[right] >= 0) {
      // dispatch the examples on newly opened nodes

      ES->checkExamples(node, 1);
      if(child[left] >= 0)
	ES->checkExamples(child[left], 1);
      if(child[right] >= 0)
	ES->checkExamples(child[right], 1);
      
      for(polarity = negative; polarity <= positive; ++polarity) {
	n = ES->Example[polarity][node].size;
	for(i=0; i<n; ++i) {
	  j = ES->Example[polarity][node][i];
	  k = ES->idx_examples[polarity][j];

	  if( feature < nFeatures ) {
	    direction = ES->bin_examples[k][feature];
	  }
	  else
	    direction = (ES->int_examples[k][ES->bin2feature[feature-nFeatures]] <= ES->bin2value[feature-nFeatures]);


	  if( child[direction] >= 0 ) {
	    ES->Example[polarity][child[direction]].insert( j );
	  }
	}
      }
      
#ifdef _DEBUGSEARCH_
      //ES->checkExamples(node, 2);
      if(child[left] >= 0) {
	cout << "Node_" << child[left] << ": " << endl;
	ES->printExamples(child[left]);
      }
      if(child[right] >= 0) {
	cout << "Node_" << child[right] << ": " << endl;
	ES->printExamples(child[right]);
      }
#endif
      
      // update the entropy of the newly opened nodes
      for(direction = left; consistent && direction <= right; ++direction)
	if( child[direction] >= 0 ) {
	  
	  ES->updateEntropy( child[direction] );
	  consistent = ( pruneFeature( child[direction] ) );
	  
#ifdef _DEBUGSEARCH_
	  if(consistent) {
	    cout << "\tupdate entropy";	
	    ES->printInfo( child[direction] );
	  }
#endif
	}
    }

#ifdef _DEBUGSEARCH_
      cout << (consistent ? "ok" : "inconsistent!") << endl;
#endif

    return consistent;
  }

  inline bool propagateFeature(const int node, const int feature)
  {
    return ( openNewNodes(node, feature) &&
	     dispatchExamples(node, feature) );
  }

  inline bool propagate(const int changedIdx, const int e) 
  { 
    int consistent = true;
    int node = changedIdx%nNodes;
    int value = scope[changedIdx]->value();

#ifdef _DEBUGSEARCH_
    cout << "propagate " ;
    if(changedIdx < nNodes)
      cout << "Parent";
    else if(changedIdx < 2*nNodes)
      cout << "Left Child";
    else if(changedIdx < 3*nNodes)
      cout << "Right Child";
    else if(changedIdx < 4*nNodes)
      cout << "Num Child";
    else 
      cout << "Feature";
    cout << "[" << node << "] = ";
    cout <<     scope[changedIdx]->value();
    cout << endl;
#endif

    if(changedIdx < 3*nNodes) {

    } else if(changedIdx < 4*nNodes) {
      int i, nchilds = 0;
      for(i=0; i<nNodes; ++i)
	nchilds += NumChild[i]->min();
      for(i=1; consistent && i<=nchilds; ++i)
	consistent = Feature[i]->remove( nFeatures );
    } else {
      consistent = (value == nFeat || propagateFeature(node, value));
    }

    //exit(0);
    return consistent;
  }
  
  inline int check( const int* s ) const 
  { 
    return 0;
  }

  virtual void print(std::ostream& o) const
  {
    o << "Global Learning Constraint";
  }
  
};




MinHittingSet::MinHittingSet(DataSet& ds) {
  init( ds.nBinFeatures, ds.columns, ds.nDiscrepancies );
}

MinHittingSet::MinHittingSet(ExampleStruct *es) {
  init( (nFeat-mFeat), es->columns, es->nDiscrepancies );
}

MinHittingSet::~MinHittingSet()
{

#ifdef _CPLEX

  free(obj);
  free(rhs);
  free(sense);
//   free(matbeg);
//   free(matcnt);
//   free(matind);
//   free(matval);
  free(zlb);
  free(zub);
  free(ctype);
  delete [] x;
  delete [] dj;
  delete [] pi;
  delete [] slack;

  CPXcloseCPLEX( &the_env );

#endif
}

void MinHittingSet::init( const int nFeatures,
			  vector<int> *columns,
			  const int nDiscrepancies) {
  {
      int i, j, k;
      
      epsilon = 0.000001;
      nColumns=(nFeat-mFeat);
      nRows=nDiscrepancies;

      x     = new double[nColumns];
      dj    = new double[nColumns];
      pi    = new double[nRows];
      slack = new double[nRows];


      /*************************
       * CPlex structs
       *************************/  


      obj     = (double *) malloc (nColumns * sizeof(double));
      rhs     = (double *) malloc (nRows    * sizeof(double));
      sense   = (char *)   malloc (nRows    * sizeof(char)); 
      matbeg  = (int *)    malloc (nColumns * sizeof(int));
      matcnt  = (int *)    malloc (nColumns * sizeof(int));   
      zlb     = (double *) malloc (nColumns * sizeof(double));
      zub     = (double *) malloc (nColumns * sizeof(double));
      ctype   = (char *)   malloc (nColumns * sizeof(char));

#ifdef _CPLEX

      the_env = NULL;
      the_env = CPXopenCPLEX (&status);

      if ( the_env == NULL ) {
      	char  errmsg[1024];
      	fprintf (stderr, "Could not open CPLEX environment.\n");
      	CPXgeterrorstring (the_env, status, errmsg);
      	fprintf (stderr, "%s", errmsg);
	exit(0);
      } 
//       else {
//       	fprintf (stderr, 
//       		 "CPlex environment open.\n");
//       }
      
      status = CPXsetintparam (the_env, CPX_PARAM_SCRIND, CPX_ON);
      if ( status ) {
      	fprintf (stderr, 
      		 "Failure to turn on screen indicator, error %d.\n", status);
	exit(0);
      } 
//       else {
//       	fprintf (stderr, 
//       		 "Param initialised.\n");
//       }
      
      the_lp = CPXcreateprob (the_env, &status, "bidule");
      if ( the_lp == NULL ) {
      	fprintf (stderr, "Failed to create problem.\n");
	exit(0);
      } 
//       else {
//       	fprintf (stderr, "Problem created.\n");
//       }

      
      //////////////// CREATING THE LP /////////////////

      CPXsetintparam( the_env, CPX_PARAM_BARDISPLAY, CPX_OFF );
      CPXsetintparam( the_env, CPX_PARAM_CONFLICTDISPLAY, CPX_OFF );
      CPXsetintparam( the_env, CPX_PARAM_MIPDISPLAY, CPX_OFF );
      CPXsetintparam( the_env, CPX_PARAM_NETDISPLAY, CPX_OFF );
      CPXsetintparam( the_env, CPX_PARAM_SIFTDISPLAY, CPX_OFF );
      CPXsetintparam( the_env, CPX_PARAM_SIMDISPLAY, CPX_OFF );

#endif      

      // minimise
      objsen = 1;
      k=0;
      for(i=0; i<nColumns; ++i) { // for each columns	

	// all variables are Boolean
	ctype[i] = 'B';

	// [0..1] variables
	zlb[i] = 0.0;
	zub[i] = 1.0;
	
	// objective: minimise the size of the hitting set
	obj[i] = 1.0;

	//
	matbeg[i] = k;

	// number of elements in each row
	matcnt[i] = columns[i].size();
	k += matcnt[i];
      }

      for(i=0; i<nRows; ++i) {
	// exactly one value per domain
	rhs[i] = 1.0;
	
	// all equalities
	sense[i] = 'G';
      }
      

      matind = (int*)   malloc(k * sizeof(int));   
      matval = (double*)malloc(k * sizeof(double));
      std::fill(matval, matval+k, 1.0);
      
      for(i=0; i<nColumns; ++i) { // for each columns	
	memcpy(matind+matbeg[i], &(columns[i][0]), matcnt[i]*sizeof(int));
      }


      //assert( j == k );

#ifdef _CPLEX

//       status = CPXcheckcopylp(the_env, the_lp, nColumns, nRows, objsen, 
// 			      obj, rhs, sense, matbeg, matcnt, matind, 
// 			      matval, zlb, zub, NULL);

      status = CPXcopylp(the_env, the_lp, nColumns, nRows, objsen, 
			 obj, rhs, sense, matbeg, matcnt, matind, 
			 matval, zlb, zub, NULL);

      if ( status ) {
      	fprintf (stderr, 
      		 "Failure to copy data, error %d.\n", status);
	exit(0);
      } 
//       else {
//       	fprintf (stderr, 
//       		 "Data copied.\n");
//       }


      free(matbeg);
      free(matcnt);
      free(matind);
      free(matval);


#endif

    }
  }


  void MinHittingSet::print()
  {

#ifdef _CPLEX

    status = CPXwriteprob(the_env, the_lp, "myprob.lp", NULL);

      if ( status ) {
      	fprintf (stderr, 
      		 "Failure to write problem to file.\n", status);
	exit(0);
      } 
//       else {
//       	fprintf (stderr, 
//       		 "Problem written.\n");
//       }

#endif

  }

  double MinHittingSet::getBarrierBound()
  {

#ifdef _CPLEX

    //print();

    ++ENDCPLEXCALLS;
    status = CPXbaropt (the_env, the_lp);

    if ( status ) {
      fprintf (stderr, 
	       "Failure to solve linear relaxation.\n", status);
      exit(0);
    } 
//     else {
//       fprintf (stderr, 
// 	       "Linear relaxation solved.\n");
//     }

    status = CPXsolution (the_env, the_lp, &lpstat, &objval, x, pi, slack, dj);

    //cout << 33 << endl;

    if ( status ) {
      fprintf (stderr, 
	       "Failure to access solution.\n", status);
      exit(0);
    } 
//     else {
//       fprintf (stderr, 
// 	       "Solution loaded.\n");
//     }

    CPXfreeprob(the_env, &the_lp);

#endif

    return objval;

  }



  double MinHittingSet::getSimplexBound()
  {

#ifdef _CPLEX

    //print();
    ++ENDCPLEXCALLS;
    status = CPXlpopt (the_env, the_lp);
    
    if ( status ) {
      fprintf (stderr, 
	       "Failure to solve linear relaxation.\n", status);
      exit(0);
    } 
//     else {
//       fprintf (stderr, 
// 	       "Linear relaxation solved.\n");
//     }

    status = CPXsolution (the_env, the_lp, &lpstat, &objval, x, pi, slack, dj);

    //cout << 33 << endl;

    if ( status ) {
      fprintf (stderr, 
	       "Failure to access solution.\n", status);
      exit(0);
    } 
//     else {
//       fprintf (stderr, 
// 	       "Solution loaded.\n");
//     }

    CPXfreeprob(the_env, &the_lp);

#endif

    return objval;

  }

  double MinHittingSet::solve()
  {

#ifdef _CPLEX

    status = CPXcopyctype (the_env, the_lp, ctype);
    
    if ( status ) {
      fprintf (stderr, 
	       "Failure to set integer vars.\n", status);
      exit(0);
    } 
//     else {
//       fprintf (stderr, 
// 	       "Integrality constraint.\n");
//     }

    ++ENDCPLEXCALLS;
    status = CPXmipopt (the_env, the_lp);

    if ( status ) {
      fprintf (stderr, 
	       "Failure to solve the problem.\n", status);
      exit(0);
    }
//     else {
//       fprintf (stderr, 
// 	       "Problem solved.\n");
//     }

    status = CPXgetbestobjval (the_env, the_lp, &objval);

    if ( status ) {
      fprintf (stderr, 
	       "Failure to access objective value.\n", status);
      exit(0);
    } 
//     else {
//       fprintf (stderr, 
// 	       "Objective value loaded.\n");
//     }

    status = CPXgetx (the_env, the_lp, x, 0, CPXgetnumcols(the_env, the_lp)-1);

    if ( status ) {
      fprintf (stderr, 
	       "Failure to access solution.\n", status);
      exit(0);
    } 
//     else {
//       fprintf (stderr, 
// 	       "Solution loaded.\n");
//     }

    CPXfreeprob(the_env, &the_lp);

#endif

    return objval;

  }




/**********************************************
 * Learning Constraint Declaration
 **********************************************/
class BuildObjectLearning : public BuildObjectConstraint {

 public:
  
  ConstraintLearning *conptr_;
  ExampleStruct *exptr_;
  
  BuildObjectLearning( ExampleStruct *es )
  {
    exptr_ = es;
  }

  virtual void build( Solver *s, VariableInt **tmp, 
		      BuildObjectPredicate *pred )
  {
    conptr_ = new ConstraintLearning( s, tmp, exptr_ );
  }

  virtual void print(std::ostream& o,
		     BuildObject **scope, 
		     const int *params,
		     const int n, const int info) const 
  {
    o << "Global Learning Constraint";
  }

};
/**********************************************
 * Learning Constraint Wrapper
 **********************************************/
class Learning : public Variable {

 public:
 
  Learning( VarArray& x, VarArray& y, VarArray& z, 
	    VarArray& t, VarArray& r, 
	    ExampleStruct *es )
      {

	int i, n = x.size();
	assert( n == y.size() );
	BuildObject **scope = new BuildObject*[5*n+1];
	for(i=0; i<n; ++i) 
	  scope[i] = x[i].var_ptr_;
	for(i=0; i<n; ++i) 
	  scope[n+i] = y[i].var_ptr_;
	for(i=0; i<n; ++i) 
	  scope[2*n+i] = z[i].var_ptr_;
	for(i=0; i<n; ++i) 
	  scope[3*n+i] = t[i].var_ptr_;
	for(i=0; i<n; ++i) 
	  scope[4*n+i] = r[i].var_ptr_;

	int idx = registerConstraint(new BuildObjectLearning(es));
	var_ptr_ = new BuildObjectPredicate( scope, 5*n, 0, -1, 
					     getConstraint(idx), NULL );
      }  
};





// This constraint does the following:
// 1 - it keeps a list of examples for each node.
// 2 - it prunes the left and right childs of a node according to the projected entropy information
class ConstraintMhs : public Constraint {

public:

           int         nNodes; 
  unsigned int         nFeatures; 
  VariableInt       ** Parent;
  VariableInt        * Objective; 
  ExampleStruct      * ES;
  ReversibleNum<int> * nDescendants;
  BitSet               internal;


public:
  

  /**@name Constructors*/
  //@{
  ConstraintMhs(Solver *s, VariableInt **x, 
		const int ln, 
		ExampleStruct *es,
		const int mode, const int thresh)
    : Constraint(s, x, ln, Constraint::VALUETRIGGER)
  {
    nNodes = ln-1;
    Parent = scope;
    ES = es;
    ES->xorThreshold = thresh;
    ES->Mode = mode;
    Objective = scope[nNodes];
    nDescendants = new ReversibleNum<int>[nNodes];
    for(int i=0; i<nNodes; ++i)
      {
	s->binds( nDescendants[i] );
	nDescendants[i].setValue( 0 );
      }
    internal.init(0, nNodes-1, BitSet::empt);
  }

  virtual ~ConstraintMhs()
  {
    delete [] nDescendants;
  }
  //@}

  inline bool propagate()
  {
    bool consistent = true;
    return consistent;
  } 

  inline bool propagate(const int changedIdx, const int e) 
  { 
    int consistent = true;

#ifdef _DEBUGSEARCH_
    cout << endl << "== propagate " ;
    if(changedIdx < nNodes) 
      cout << "Parent";
    else 
      cout << "Objective";
    cout << "[" << changedIdx << "] = ";
    cout <<     scope[changedIdx]->value();
    cout << endl;

#endif

    if(changedIdx < nNodes) 
      {
	int reqNodes = 0;//changedIdx+1;
	int node = changedIdx;
	int parent;// = scope[node]->value();
	//assert( parent != node )
	//cout << 11 << endl;

	//int lb2 = ES->lowerBound2( node );
	int lb = ES->lowerBound1( node );

// 	cout << lb1 << " " << lb2 << endl;

// 	assert( lb1 == lb2 );
	 

// 	int lb = lb1;

// 	//cout << 22 << endl;
	
	nDescendants[node] = lb;

	node = nNodes;

	while( node-- && reqNodes <= nNodes )
	  {
	    if( !scope[node]->isGround() ) continue;

	    parent = scope[node]->value(); 
	    if(parent != node)
	      internal.insert( parent );

	    if( internal.member(node) ) continue;

	    
	    //cout << "open node (" << node << ") < " << parent << endl;

	    lb = nDescendants[node];
	    reqNodes += lb;
	    //int lb = ES->lowerBound( node );

	    //cout << node << "  ===> " << lb << endl << endl;
	    

	    //--node;
	    //parent = scope[node]->value(); 
	  }
	reqNodes += internal.size();
	internal.clear();

//  	cout << "  =======> " << reqNodes << " / " << (nNodes-reqNodes) << endl ;
// 	Objective->print(cout);
// 	cout << endl << endl;


	consistent = Objective->setMax(nNodes - reqNodes);




      }

    //if(!consistent) exit(0);
    return consistent;
  }
  
  inline int check( const int* s ) const 
  { 
    return 0;
  }

  virtual void print(std::ostream& o) const
  {
    o << "Global Mhs Constraint";
  }
  
};


/**********************************************
 * Mhs Constraint Declaration
 **********************************************/
class BuildObjectMhs : public BuildObjectConstraint {

 public:
  
  ConstraintMhs *conptr_;
  ExampleStruct *exptr_;
  
  BuildObjectMhs( ExampleStruct *es )
  {
    exptr_ = es;
  }

  virtual void build( Solver *s, VariableInt **tmp, 
		      BuildObjectPredicate *pred )
  {
    conptr_ = new ConstraintMhs( s, tmp, pred->arity, exptr_, pred->params[0], pred->params[1] );
  }

  virtual void print(std::ostream& o,
		     BuildObject **scope, 
		     const int *params,
		     const int n, const int info) const 
  {
    o << "Global Mhs Constraint";
  }

};
/**********************************************
 * Mhs Constraint Wrapper
 **********************************************/
class Mhs : public Variable {

 public:

  Mhs( VarArray& x, Variable& z, ExampleStruct *es, const int mode, const int thresh )
  {
    int i, n = x.size();
    BuildObject **scope = new BuildObject*[n+2];
    for(i=0; i<n; ++i) 
      scope[i] = x[i].var_ptr_;
    scope[n] = z.var_ptr_;
    int *params = new int[2];
    params[0] = mode;
    params[1] = thresh;
    var_ptr_ = new BuildObjectPredicate( scope, n+1, 0, -1, 
					 getConstraint( registerConstraint( new BuildObjectMhs(es) ) ),
					 params );
  }
};






class ValSelectorDtree : public ValSelector
{
public:

  ExampleStruct *ES;
  int node;
  int val;
  int randomized;
  int *bestVal;
  double *bestEnt;

  virtual ~ValSelectorDtree()
  {
    delete [] bestVal;
    delete [] bestEnt;
  }

  ValSelectorDtree( VariableInt *x, ExampleStruct *es, const int n, const int r )
    : ValSelector(x) 
  { 
    ES = es;
    node=n; 
    randomized = r;
    if( randomized < 1 )
      randomized = 1;

    bestVal = new int[randomized+1];
    bestEnt = new double[randomized+1];
  }

  /**@name Utils*/
  //@{ 
//   inline int getBest() 
//   { 
//     int second_val; 
//     if( _X->contain(nFeat) ) {
//       val = nFeat;
//       return val;
//     } else {
//       double mig1 = -1, mig2 = -1, ig;
//       int feature = _X->min();
//       second_val = feature; 
//       val = feature;
//       do {
// 	ig = ES->informationGain(node, feature);
// 	if( ig > mig1 ) { 
// 	  mig2 = mig1;
// 	  second_val = val;
// 	  mig1 = ig;
// 	  val = feature;	  
// 	} else if( ig > mig2 ) {
// 	  mig2 = ig;
// 	  second_val = feature;
// 	}
//       } while( _X->setNext( feature ) );
//     }

//     cout << val << " " << second_val << " " << 2 << " " ;

//     if(randomized & 2) 
//       if(rand()%2) {
// 	cout << second_val << endl;
// 	return second_val;
//       }
//     cout << val << endl;
//     return val;
//   }


  inline int getBest() 
  { 
    int i, n_val = 0;   
    if( _X->contain(nFeat) ) {
      val = nFeat;
      return val;
    } else {

//       _X->print( cout );
//       cout << endl;

      std::fill(bestEnt, bestEnt+randomized+1, -1);      
      int feature = _X->min();
      double aux;
      do {
	i = n_val;
	
	bestEnt[i] = ES->informationGain(node, feature);

	//cout << "Value: " << feature << " - " << bestEnt[i]  << endl;

	while( i && bestEnt[i-1] <= bestEnt[i] )
	  {
	    aux = bestEnt[i-1];
	    bestEnt[i-1] = bestEnt[i];
	    bestEnt[i] = aux;
	    bestVal[i] = bestVal[i-1];

// 	    for(int j=0; j<n_val; ++j)
// 	      if( i == j )
// 		cout << " " << feature;
// 	      else
// 		cout << " " << bestVal[j];
// 	    cout << endl;
// 	    for(int j=0; j<n_val; ++j)
// 	      cout << " " << setprecision(4) << bestEnt[j] ;
// 	    cout << endl;

	    --i;
	  }
	bestVal[i] = feature;
       
	if( n_val < randomized )
	  ++n_val;
	

// 	cout << " => ";
// 	for(int j=0; j<n_val; ++j)
// 	  cout << " " << bestVal[j];
// 	cout << endl;
// 	for(int j=0; j<n_val; ++j)
// 	  cout << " " << setprecision(4) << bestEnt[j] ;
// 	cout << endl << endl;


      } while( _X->setNext( feature ) );
    }

    //cout << bestVal[0] << " " << bestVal[1] << " " << n_val << " ";
    
    val = bestVal[rand()%n_val];

    //cout << val << endl;

    return( val );
  }


  inline int getBest() const 
  { 
    int second_val, val_=0;
    if( _X->contain(nFeat) ) {
      val_ = nFeat;
      return val_;
    } else {
      double mig1 = -1, mig2 = -1, ig;
      int feature = _X->min();
      second_val = feature; 
      val_ = feature;
      do {
	ig = ES->informationGain(node, feature);
	if( ig > mig1 ) { 
	  mig2 = mig1;
	  second_val = val_;
	  mig1 = ig;
	  val_ = feature;	  
	} else if( ig > mig2 ) {
	  mig2 = ig;
	  second_val = feature;
	}
      } while( _X->setNext( feature ) );
    }
    return val_; 
  }
  inline void left() 
  {
    _X->setDomain( getBest() ); 
  }
  inline void right() 
  { 
    _X->remove( val ); 
  }
  void postCut( const int p ) {
    val = p;
  }
  //@}  

  virtual void printLeft(std::ostream& o) const { 
    _X->printshort(o);
    cout << " == " << getBest();
  }
  virtual void printRight(std::ostream& o) const {
    _X->printshort(o);
    cout << " =/= " << val;
  }
};





/**********************************************
 * Search strategy:
 * 1 - Pick an unused node
 *      . 
 **********************************************/

class DTreeVO : public DVO
{
 public: 

  /**@name Parameters*/
  //@{ 
  int nNodes;
  VariableInt **Parent;
  VariableInt **LeftChild;
  VariableInt **RightChild;
  VariableInt **NumChild;
  VariableInt **Feature;
  VariableInt **NotUsed;
  ExampleStruct *ES;
  int randomized;

  int *openNode;
  int nOpen;
  int nClose;
  //@}

  /**@name Constructors*/
  //@{
  DTreeVO(Solver* s,
	  const int n, 
	  VariableInt **p,
	  VariableInt **l,
	  VariableInt **r,
	  VariableInt **k,
	  VariableInt **f,
	  VariableInt **u,
	  ExampleStruct *es,
	  const int rd) : DVO(s) 
  {
    nNodes = n;
    Parent = p;
    LeftChild = l;
    RightChild = r;
    NumChild = k;
    Feature = f;
    NotUsed = u;
    openNode = new int[nNodes];
    nOpen = 0;
    nClose = 0;
    ES = es;
    randomized = rd;

    for(int i=0; i<nNodes; ++i)
      Feature[i]->branch = new ValSelectorDtree( Feature[i], ES, i, randomized );

  }

  virtual ~DTreeVO() {
    delete [] Parent;
    delete [] LeftChild;
    delete [] RightChild;
    delete [] NumChild;
    delete [] Feature;
    delete [] NotUsed;
    delete [] openNode;
  }
  //@}

  /**@name Utils*/
  //@{ 
  inline VariableInt* select()
  {
#ifdef _DEBUGSEARCH_
    cout << "\nselect:\nnode | parent |  left  |  right | NumChild | feature | NotUsed?" << endl;
    int encore = true;
    for(int i=0; encore && i<nNodes; ++i)  {
      cout << setw(4) << i << " |";
      cout << " ";
      Parent[i]->printDomain( cout );
      cout << " |";
      cout << " ";
      LeftChild[i]->printDomain( cout );
      cout << " |";
      cout << " ";
      RightChild[i]->printDomain( cout );
      cout << " |";
      cout << " ";
      NumChild[i]->printDomain( cout );
      cout << " |";
      cout << " ";
      if( Feature[i]->isGround())
	cout << "{" << Feature[i]->value() << "}";
      else 
	cout << "[ " << Feature[i]->min() << ".." << Feature[i]->max() << "]";
      //Feature[i]->printDomain( cout );
      cout << " |";
      cout << " ";
      NotUsed[i]->printDomain( cout );
      cout << "  ";
      ES->printInfo( i, false );
    }

#endif
 
    // find an open node, and assign a feature.
    nOpen = 0;
    nClose = 0;
    int i, j;
    for(i=0; i<nNodes; ++i) 
      if( !Feature[i]->isGround() && Parent[i]->isGround() ) {
	j = Parent[i]->value();
	if( !j || Feature[j]->isGround() ) {
	  openNode[nOpen++] = i;
	  //assert(!i || Parent[i]->value() != i);
	  //cout << i << " ";
	} else {
	  ++nClose;
	  openNode[nNodes-nClose] = i;

	  //cout << " <<== " << i << endl; 

	}
      }
    

    if( !nOpen ) {
      if( !nClose ) {
	//cout << "FIRST!!!" << endl;
	return *first;
      } else return Feature[openNode[nNodes-1]];      
    }


    int ridx;
    //     if(randomized & 1) {

    //       // doing this (assigning nodes in non-lexico order) does not work
    //       // well because of the symmetry breaking constraints
    //       // if we do not assign nodes in order, then we can't be sure what 
    //       // node we should open next

    //       cerr << "BUGGY!!" << endl;
    //       //exit(0)

    //       ridx = (rand()%nOpen);
    //     } else
    ridx = 0;
    //cout << " => " << openNode[ridx] << endl;

    //if(openNode[ridx] == 39) exit(0);

#ifdef _DEBUGSEARCH_

    cout << setw(4) << openNode[ridx] << " |";
    cout << " ";
    Parent[openNode[ridx]]->printDomain( cout );
    cout << " |";
    cout << " ";
    LeftChild[openNode[ridx]]->printDomain( cout );
    cout << " |";
    cout << " ";
    RightChild[openNode[ridx]]->printDomain( cout );
    cout << " |";
    cout << " ";
    NumChild[openNode[ridx]]->printDomain( cout );
    cout << " |";
    cout << " ";
    if( Feature[openNode[ridx]]->isGround())
      cout << "{" << Feature[openNode[ridx]]->value() << "}";
    else 
      cout << "[ " << Feature[openNode[ridx]]->min() << ".." << Feature[openNode[ridx]]->max() << "]";
    //Feature[openNode[ridx]]->printDomain( cout );
    cout << " |";
    cout << " ";
    NotUsed[openNode[ridx]]->printDomain( cout );
    cout << "  ";
    ES->printInfo( openNode[ridx], true );
    cout << endl;      
    ES->printExamples(openNode[ridx]);

#endif

    return Feature[openNode[ridx]];
  }
  //@}
};


class DecisionTreeStrategy : public VariableOrdering {
 public:

  int nNodes;
  VariableInt **tmpPrnt;
  VariableInt **tmpLeft;
  VariableInt **tmpRigh;
  VariableInt **tmpNumC;
  VariableInt **tmpFeat;
  VariableInt **tmpUsed;
  ExampleStruct *ES;
  int randomized;

  DecisionTreeStrategy( const int n, 
			VariableInt **p,
			VariableInt **l,
			VariableInt **r,
			VariableInt **k,
			VariableInt **f,
			VariableInt **u,
			ExampleStruct *es,
			const int rd
			) { 
    nNodes=n; 
    tmpPrnt=p;
    tmpLeft=l;
    tmpRigh=r;
    tmpNumC=k;
    tmpFeat=f;
    tmpUsed=u;
    ES = es;
    randomized = rd;
  }

  /**@name Utils*/
  //@{
  DVO* extract( Solver* s ) {    
    return new DTreeVO( s, nNodes, tmpPrnt, tmpLeft, tmpRigh, tmpNumC, tmpFeat, tmpUsed, ES, randomized );
  }
  //@}
};


void removeColumn( const char *filename )
{
  int i;
  char buf[1024]; 
  ifstream datafile( filename );  
  assert( datafile.good() );
  
  while( true ) {
    datafile.getline( buf, 1024 );
    if( datafile.eof() || !datafile.good() ) break;

    for(i=0; buf[i]!=','; ++i);
    cout << &(buf[i+1]) << endl;

  }
}

void firstToLastColumn( const char *filename )
{
  int i, j, k;
  char buffer[2048]; 
  ifstream datafile( filename );  
  assert( datafile.good() );
  
  bool first;

  while( true ) {
    datafile.getline( buffer, 2048 );
    if( datafile.eof() || !datafile.good() ) break;

    for(i=0; buffer[i]!=','; ++i) {}
//     //cout << i << endl;
//     //cout << "+ " << buffer[i] << endl;
//     //assert( buffer[i] == ',' );
//     //cout << &(buffer[i+1]) << endl;
//     //cout << i << endl;
    buffer[i] = '\0';
    cout << &(buffer[i+1]) << "," << buffer << endl;

//     //    cout << "&buffer[i+1]: " << &(buffer[i+1]) << endl;
// //     cout.flush();
// //     cout << "&buffer[0]: " << &(buffer[0]) << endl << endl;
    
//     //cout << buffer << endl;
//     //0,0,0,1,1,0,0,1,0,0,0,0,1,1,1,0,0,0,0,0,0,1,1
//     for(j=0; buffer[j]=='0' || buffer[j]=='1' || buffer[j]==','; ++j) {}
//     //cout << j << endl;

//     //cout << buffer[45] << endl;

//     buffer[j] = ',';
//     ++j;

//     for(k=0; k<i; ++k) {
//       buffer[j+k] = buffer[k];
//       //cout << "buffer[" << (j+k) << "] = buffer[" << k << "];" << endl;
//     }
//     buffer[j+k] = '\0';
    
//     cout << &(buffer[2]) << endl;
  }
}


int  printTree( int lvl,
		int x,
		DataSet& ds,
		int *solution,
		int n,
		ExampleStruct *es,
		ofstream& outfile);


int solveWithLP( DataSet& ds )
{

  ds.computeXor();

#ifdef _CPLEX
   MinHittingSet mhs(ds);
   mhs.print();
   double result = mhs.solve();
   cout << result << endl;

   //for(i=0; i<ds.nBinFeatures; ++i)
#endif

}


int solveWithMiniSAT( DataSet& ds, int Depth)
{
  ds.computeXor();
  ds.satEncode(Depth);
}


int solveWithSAT( DataSet& ds, int Depth)
{

  ds.computeXor();

  int result;
  do {
    SatSolver dll;
    dll.params.verbosity = 4;
    dll.setPolicy( GEOMETRIC );
    //ds.satEncode(Depth);
    ds.satEncode(Depth, &dll);
    result = dll.solve();
  } while(Depth-- && result == SAT);
  
}


int solveWithCP( DataSet& ds, int nNodes, 
		 const int mode, 
		 const int thresh,
		 const int randomized,
		 const double timelimit,
		 const int verbose )
{

  int i, j;

  // create the csp
  CSP model;
  //cout << split << " " << nFeat << endl;

  VarArray Parent( nNodes, 0, nNodes-1 );
  VarArray LeftChild( nNodes, 0, nNodes-1 );
  VarArray RightChild( nNodes, 0, nNodes-1 );
  VarArray NumChild( nNodes, 0, 2 );
  VarArray Feature( nNodes, mFeat, nFeat );
  VarArray Descendant( nNodes*nNodes, 0, 1 );
  VarArray Children[nNodes];
  VarArray NotUsed( nNodes, 0, 1 );


  //////////////////////////////////////////////////////////////
  //////////// THE TREE PART ///////////////////////////////////
  //////////////////////////////////////////////////////////////
  // The array 'Parent' should form a tree.
  model.add( Tree( Parent, Descendant ) );

  // Moreover, it should be binary.
  for( i=0; i<nNodes; ++i ) {

    // symmetry breaking?
    if(i) model.add( Parent[i-1] <= Parent[i] );
    model.add( Parent[i] <= i );
    model.add( LeftChild[i] <= (2*i+2) );
    model.add( RightChild[i] <= (2*i+1) );
    model.add( LeftChild[i] >= i );
    model.add( RightChild[i] >= i );
    model.add( (LeftChild[i] == i) || (LeftChild[i] >= RightChild[i]) );

    //VarArray children_i;
    for( j=0; j<nNodes; ++j ) 
      if( i != j ) {
	Variable x = (Parent[j] == i);
	model.add( x == ((LeftChild[i] == j) + (RightChild[i] == j)) );
	Children[i].add( x );
      } else {
	Variable x(0,0);
	Children[i].add( x );
      }    
    model.add( Sum( Children[i] ) == NumChild[i] );
    model.add( (NumChild[i] + (LeftChild[i] == i) + (RightChild[i] == i)) <= 2 );


    for( j=0; j<i; ++j )
      model.add( Descendant[i*nNodes+j] == 0 );
  }
  //////////////////////////////////////////////////////////////
  //////////// END OF THE TREE PART ////////////////////////////
  //////////////////////////////////////////////////////////////


  //////////////////////////////////////////////////////////////
  //////////// THE FEATURES PART ///////////////////////////////
  //////////////////////////////////////////////////////////////
  for( i=0; i<nNodes; ++i )
    for( j=0; j<nNodes; ++j ) 
      if( i != j ) {
	//model.add( Descendant[i*nNodes+j] <= (Feature[i] != Feature[j]) );
	model.add( (NotUsed[i] + Children[j][i]) < 2 );
      }

  // The array 'Used' indicates which nodes will indeed represent a feature check.
  for( i=0; i<nNodes; ++i ) {
    model.add( (NotUsed[i] + Descendant[i]) < 2 );
    model.add( NotUsed[i] == (Feature[i] == nFeat) );
    model.add( NotUsed[i] <= (NumChild[i] == 0) );
  }

  // We cannot use the child of an unused node.
  for( i=0; i<nNodes; ++i ) 
    model.add( NotUsed[Parent[i]] <= NotUsed[i] );
  //////////////////////////////////////////////////////////////
  //////////// END OF THE FEATURES PART ////////////////////////
  //////////////////////////////////////////////////////////////
  
  
  //////////////////////////////////////////////////////////////
  //////////// THE CLASSIFICATION PART /////////////////////////
  //////////////////////////////////////////////////////////////
  ExampleStruct *ES = new ExampleStruct(ds, nNodes);
  model.add( Learning(Parent, LeftChild, RightChild, NumChild, Feature, ES) );
  //////////////////////////////////////////////////////////////
  //////////// END OF THE CLASSIFICATION PART //////////////////
  //////////////////////////////////////////////////////////////
  
  
  //////////////////////////////////////////////////////////////
  //////////// THE OBJECTIVE PART //////////////////////////////
  //////////////////////////////////////////////////////////////
  //The objective is to minimise the number of used nodes.

  Variable Obj = Sum( NotUsed );
  model.add( Obj + Sum( NumChild ) == (nNodes-1) );

  if( mode ) {
    ES->init(thresh);
    model.add( Mhs(Parent, Obj, ES, mode, thresh) );
  }
  //  model.add( Maximise( Obj ) );
  //////////////////////////////////////////////////////////////
  //////////// END OF THE OBJECTIVE PART ///////////////////////
  //////////////////////////////////////////////////////////////

  VarArray searchVars;
  searchVars.add( Feature );
  searchVars.add(LeftChild);
  searchVars.add(RightChild);
  
  Solver s( model, searchVars );
  //////////////////////////////////////////////////////////////
  //////////// THE HEURISTIC PART //////////////////////////////
  //////////////////////////////////////////////////////////////
  ES->init(&s, Feature.getVarArray());
  DecisionTreeStrategy strategy( nNodes, 
				 Parent.getVarArray(),
				 LeftChild.getVarArray(),
				 RightChild.getVarArray(),
				 NumChild.getVarArray(),
				 Feature.getVarArray(),
				 NotUsed.getVarArray(),
				 ES, randomized );
  s.add( strategy );
  //////////////////////////////////////////////////////////////
  //////////// END OF THE HEURISTIC PART ///////////////////////
  //////////////////////////////////////////////////////////////
  

  s.setVerbosity( verbose );
  cout << "c  ATTEMPT TO SOLVE " << 2*nNodes+1 << endl; 
  bool satisfiable;

  s.setTimeLimit(timelimit);

  int returnVal;

//   if( randomized )
//     satisfiable = (s.solve_and_restart(GEOMETRIC, 100, 1.5) == SAT);
//   else
//     satisfiable = (s.solve() == SAT);


  if( randomized )
    returnVal = s.solve_and_restart(GEOMETRIC, 100, 1.5);
  else
    returnVal = s.solve();


  satisfiable = (returnVal == SAT);
  
  ENDPROPAGS += s.PROPAGS ;
  ENDNODES += s.NODES;
  ENDBACKTRACKS += s.BACKTRACKS;
  ENDFAILURES += s.FAILURES;
  if( FIRSTPROPAGS < 0 ) {
    FIRSTPROPAGS = ENDPROPAGS;
    FIRSTNODES = ENDNODES;
    FIRSTBACKTRACKS = ENDBACKTRACKS;
    FIRSTFAILURES = ENDFAILURES;
    FIRSTTIME = getRunTime();
    FIRSTTIME -= STARTTIME;
    FIRSTCPLEXCALLS = ENDCPLEXCALLS;
    firstTree = 2*( nNodes-Obj.value() )+1;
  }
  if( satisfiable ) {
    OPTPROPAGS = ENDPROPAGS;
    OPTNODES = ENDNODES;
    OPTBACKTRACKS = ENDBACKTRACKS;
    OPTFAILURES = ENDFAILURES;
    OPTTIME = getRunTime();
    OPTTIME -= STARTTIME;
    OPTCPLEXCALLS = ENDCPLEXCALLS;
  }


  if( satisfiable )
    { 
      ofstream outfile( OutTree.c_str(), ios_base::trunc );
      //outfile << Obj.value() << endl;
      outfile << "#include \"d_tree.h\"" << endl << endl;
      outfile << "int classify( int *example ) {" << endl;

      int level = printTree( 1, 0, ds, s.solution, nNodes, ES, outfile );
      
      ///cout << "c DEPTH: " << level << endl;
      if( level > ENDDEPTH )
	ENDDEPTH = level;

      outfile << "}" << endl << endl
	   << "int main(int argc, const char *argv[]) {" << endl
	   << "if(argc>1) checkTree( argv[1] );" << endl
	   << "}" << endl << endl;
      returnVal = ( nNodes-Obj.value() );
    } else {
    cout << "c  FAILS" << endl;
    cout << "d OPTIMAL " << ( returnVal == UNSAT ) << endl;      
    returnVal = nNodes+1;    
  }

  delete ES;

  return returnVal;
}

int printTree(int lvl, int node, DataSet& ds, int *solution, int n, ExampleStruct *es, ofstream& outfile)
{
  int aux, rval = lvl;

  int i, leftChild=node, rightChild=node, binFeature = solution[node];
  if( binFeature < nFeat ) {

    if( node < n-1 ) {
      leftChild = solution[n+node];
      rightChild = solution[2*n-1+node];
    }

//     int value, feature=0, threshold=ds.index2choice[0].size();
//     while( binFeature >= threshold ) {
//       threshold += ds.index2choice[++feature].size();	
//     }
//     threshold -= ds.index2choice[feature].size();
//     value = (binFeature - threshold);

    int value = ds.bin2value[(binFeature%(ds.nBinFeatures))];
    int feature = ds.bin2feature[(binFeature%(ds.nBinFeatures))];    
    string val = (ds.index2choice[feature][value]);

    for(i=0; i<lvl; ++i) outfile << " ";

    if( binFeature < ds.nBinFeatures )
      outfile << "if( example[" << feature 
	      << "] == " << val << " ) {" << endl;
    else
      outfile << "if( example[" << feature 
	      << "] <= " << val << " ) {" << endl;
    
    if( rightChild > node ) {     
      aux = printTree( lvl+1, rightChild, ds, solution, n, es, outfile );
      if( aux > rval ) rval = aux;
    } else {
      for(i=0; i<=lvl; ++i ) outfile << " ";
      outfile << "return " << (ds.index2class[(es->remainder[1][1][node][binFeature] == 0)]) << ";" << endl; 
    }

    for(i=0; i<lvl; ++i ) outfile << " ";
    outfile << "} else {" << endl;    
    
    if( leftChild > node ) {
      aux = printTree( lvl+1, leftChild, ds, solution, n, es, outfile );   
      if( aux > rval ) rval = aux;
    } else {
      for(i=0; i<=lvl; ++i ) outfile << " ";
      outfile << "return " << (ds.index2class[(es->remainder[1][0][node][binFeature] == 0)]) << ";" << endl;
    }
    
    for(i=0; i<lvl; ++i ) outfile << " ";
    outfile << "}" << endl;    
  }

  return rval;
}







const int nia = 10;
const char* int_ident[nia] = {
  "-h",
  "-verbose",
  "-seed",
  "-nodes",
  "-depth",
  "-threshold",
  "-rand",
  "-split",
  "-verbose",
  "-stop"
};
const char* int_man[nia] = {
  "help message",
  "verbosity {0,1}", 
  "random seed",
  "maximum number of nodes",
  "maximum depth",
  "Size threshold for the min Hitting Set",
  "Type of randomization {0: none, 1: vars, 2: vals, 3: vars and vals}",
  "Uses domain splitting",
  "Degree of verbosity",
  "stop when finding a tree with that number of nodes"
};			      
int int_param[nia];


const int nsa = 3;
const char* str_ident[nsa] = {
  "-algo",
  "-timelimit",
  "-out"
};
const char* str_man[nsa] = {
  "algorithm to use",
  "time limit",
  "output file"
};
const char* str_param[nsa];


void outputHelpMessage()
{
  cerr << "\nUsage: \t bin/sat [path to the cnf file] <args>" << endl << endl;
  for(int i=0; i<nia; ++i)
    cerr  << setw(20) << int_ident[i] << ": \t" << int_man[i] << endl;
  for(int i=0; i<nsa; ++i)
    cerr  << setw(20) << str_ident[i] << ": \t" << str_man[i] << endl;
  cerr << endl;
}


int main(int argc, char **argv)
{

  //firstToLastColumn( argv[1] );
  //exit(0);

  int Verbosity;
  long int Seed;
  int Nodes;
  int Depth;
  int Threshold;
  int Rand;
  int Split;
  int Verbose;
  string Algo;  
  double TimeLimit;


  ENDPROPAGS = 0;
  ENDNODES = 0;
  ENDBACKTRACKS = 0;
  ENDFAILURES = 0;
  STARTTIME = getRunTime();
  

  if( argc <= 1 ) {
    cerr << "Need an argument" << endl;
    exit( 0 );
  }

  int i, j, k, l;
  int lastNodes = -1, curNodes;
  DataSet ds;

  str_param[0] = "";
  getCommandLine(int_ident, int_param, nia,
		 str_ident, str_param, nsa,
		 &(argv[1]), argc-1);
  
  if( argc < 2 || int_param[0] != NOVAL || !strcmp( argv[argc-1], "-h" ) )
    {
      outputHelpMessage();
    }
  else
    {
      Verbosity = ( int_param[1] != NOVAL ? int_param[1] : -1 );
      Seed = ( int_param[2] != NOVAL ? int_param[2] : 11041979 );	  
      Nodes = ( int_param[3] != NOVAL ? int_param[3] : 15 );
      Depth = ( int_param[4] != NOVAL ? int_param[4] : 4 );
      Threshold = ( int_param[5] != NOVAL ? int_param[5] : 10000 );
      Rand = ( int_param[6] != NOVAL ? int_param[6] : 0 );
      Split = ( int_param[7] != NOVAL ? int_param[7] : 0 );
      Verbose = ( int_param[8] != NOVAL ? int_param[8] : 0 );
      Algo = ( strcmp(str_param[0],"nil") ? str_param[0] : "cp" );
      TimeLimit = ( strcmp(str_param[1],"nil") ? (double)(atof(str_param[1])) : -1 );
      OutTree = ( strcmp(str_param[2],"nil") ? str_param[2] : "d_tree.cpp" );


      input_file = string(argv[1]);

      if(Split)
	nFeat = 2;
      else 
	nFeat = 1;

      if(Split>1)
	mFeat = 1;
      else
	mFeat = 0;

      // get the data
      ds.getData( argv[1] );
      ds.binaryEncode( 0 );


      srand(Seed);


      if( Algo == "cp" ) {
	
	if(Verbosity > 1)
	  cerr << "cp: " << Depth << " " << Nodes << endl;

	lastNodes = Nodes+1;
	do {
	  curNodes = lastNodes-1;
	  lastNodes = solveWithCP( ds, curNodes, NOTHING, 0, Rand, TimeLimit, Verbose );
	  bestTree = (2*lastNodes+1);
	  cout << "c  SOLUTION FOUND: " << (2*lastNodes+1) << endl;
	} while( bestTree > stopValue && curNodes >= lastNodes && lastNodes > 3 );

	if( bestTree <= 2*Nodes+1 )
	  cout << "s SATISFIABLE" << endl;
	else 
	  cout << "s UNKNOWN" << endl;

      } else if( Algo == "first" ) {
	
	if(Verbosity > 1)
	  cerr << "first: " << Depth << " " << Nodes << endl;

	lastNodes = Nodes+1;
	curNodes = lastNodes-1;
	lastNodes = solveWithCP( ds, curNodes, NOTHING, 0, Rand, TimeLimit, Verbose );
	bestTree = (2*lastNodes+1);
	cout << "c  SOLUTION FOUND: " << (2*lastNodes+1) << endl;
	if( bestTree <= 2*Nodes+1 )
	  cout << "s SATISFIABLE" << endl;
	else 
	  cout << "s UNKNOWN" << endl;

      } else if( Algo == "minisat" ) {

	bestTree = 111111;

	if(Verbosity > 1)
	  cerr << "sat: " << Depth << " " << Nodes << endl;

	solveWithMiniSAT( ds, Depth );

      } else if( Algo == "sat" ) {

	//bestTree = INT_MAX;

	if(Verbosity > 1)
	  cerr << "sat: " << Depth << " " << Nodes << endl;

	solveWithSAT( ds, Depth );

      } else if( Algo == "lp" ) {

	if(Verbosity > 1)
	  cerr << "lp: " << Depth << " " << Nodes << endl;

	solveWithLP( ds );

      } else if( Algo == "cp+simplex" ) {

	if(Verbosity > 1)
	  cerr << "cp+simplex: " << Depth << " " << Nodes << endl;	

	lastNodes = Nodes+1;	
	do {
	  curNodes = lastNodes-1;
	  lastNodes = solveWithCP( ds, curNodes, LP_SIMPLEX, Threshold, Rand, TimeLimit, Verbose );
	  bestTree = (2*lastNodes+1);
	  cout << "c  SOLUTION FOUND: " << (2*lastNodes+1) << endl;
	} while( bestTree > stopValue && curNodes >= lastNodes && lastNodes > 3 );


	if( bestTree <= 2*Nodes+1 )
	  cout << "s SATISFIABLE" << endl;
	else 
	  cout << "s UNKNOWN" << endl;

      } else if( Algo == "cp+barrier" ) {

	if(Verbosity > 1)
	  cerr << "cp+barrier: " << Depth << " " << Nodes << endl;	

	lastNodes = Nodes+1;	
	do {
	  curNodes = lastNodes-1;
	  lastNodes = solveWithCP( ds, curNodes, LP_BARRIER, Threshold, Rand, TimeLimit, Verbose );
	  bestTree = (2*lastNodes+1);
	  cout << "c  SOLUTION FOUND: " << (2*lastNodes+1) << endl;
	} while( curNodes >= lastNodes && lastNodes > 3 );


	if( bestTree <= 2*Nodes+1 )
	  cout << "s SATISFIABLE" << endl;
	else 
	  cout << "s UNKNOWN" << endl;

      } else if( Algo == "cp+mip" ) {

	if(Verbosity > 1)
	  cerr << "cp+mip: " << Depth << " " << Nodes << endl;	

	lastNodes = Nodes+1;       
	do {
	  curNodes = lastNodes-1;
	  lastNodes = solveWithCP( ds, curNodes, MIP, Threshold, Rand, TimeLimit, Verbose );
	  cout << "c  SOLUTION FOUND: " << (2*lastNodes+1) << endl;
	  bestTree = (2*lastNodes+1);
	} while( bestTree > stopValue && curNodes >= lastNodes && lastNodes > 3 );

	if( bestTree <= 2*Nodes+1 )
	  cout << "s SATISFIABLE" << endl;
	else 
	  cout << "s UNKNOWN" << endl;

      } else if( Algo == "print" ) {

	ds.printNumeric( cout );
	exit(0);

      } else {
	cerr << "Algorithm not defined!" << endl;
      }
    }


  ENDTIME = (getRunTime() - STARTTIME);

  if( bestTree <= 2*Nodes+1 ) 
    cout << "d TREESIZE " << (2*lastNodes+1)
	 << " FIRSTTREE " << firstTree
	 << " FIRSTPROPAGS " << FIRSTPROPAGS
	 << " FIRSTNODES " << FIRSTNODES
	 << " FIRSTBACKTRACKS " << FIRSTBACKTRACKS
	 << " FIRSTFAILURES " << FIRSTFAILURES
	 << " FIRSTTIME " << FIRSTTIME
	 << " FIRSTCPLEXCALLS " << FIRSTCPLEXCALLS
	 << " OPTPROPAGS " << OPTPROPAGS
	 << " OPTNODES " << OPTNODES
	 << " OPTBACKTRACKS " << OPTBACKTRACKS
	 << " OPTFAILURES " << OPTFAILURES
	 << " OPTTIME " << OPTTIME
	 << " OPTCPLEXCALLS " << OPTCPLEXCALLS
	 << " ENDPROPAGS " << ENDPROPAGS
	 << " ENDNODES " << ENDNODES
	 << " ENDBACKTRACKS " << ENDBACKTRACKS
	 << " ENDFAILURES " << ENDFAILURES
	 << " TOTALTIME " << ENDTIME 
	 << " ENDCPLEXCALLS " << ENDCPLEXCALLS
	 << " DEPTH " << ENDDEPTH
	 << endl;

  return 0;
}

// income [0.01 .. 0.015]
// breast [0.2 .. 1]
// car [0.05 .. 0.9]
// forest-fire : probably not sound

