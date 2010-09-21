
#include <mistral_sol.h>



#ifdef _CPLEX

#include "check.h"
#include <ilcplex/cplex.h>

#endif

using namespace std;
using namespace Mistral;

class DataSet
{
public:

  DataSet() { discrepancies = NULL; }
  virtual ~DataSet() { 
    delete [] binExamples; 
    delete [] discrepancies;
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
  vector< vector < int > > examples;
  

  /////////// binary encoding
  Bitset < unsigned int > *binExamples;
  vector< Bitset < unsigned int > > positive_example;
  vector< Bitset < unsigned int > > negative_example;
  unsigned int nBinFeatures;


  int *exampleIndex;
  Bitset< unsigned int > *discrepancies;
  vector<int>            *columns;
  int nDiscrepancies;

  int getFeature( const int binFeature, string& value )
  {
    int feature=0, threshold=index2choice[0].size();
    while( binFeature >= threshold ) 
      threshold += index2choice[++feature].size();	
    threshold -= index2choice[feature].size();
    value = index2choice[feature][(binFeature - threshold)];

    return feature;
  }
  

  void binaryEncode( int classe )
  {

    unsigned int i,j,k;
    nBinFeatures = 0; 
    for(i=0; i<index2choice.size(); ++i)
      if(index2choice[i].size() > 2)
	nBinFeatures += index2choice[i].size();
      else if(index2choice[i].size() == 2)
	++nBinFeatures;
    
    //    cout << nBinFeatures << endl;

    binExamples = new Bitset < unsigned int >[examples.size()];
    for(i=0; i<examples.size(); ++i)
      {

	//	cout << i << ": " << endl;

	binExamples[i].init(0, nBinFeatures-1, Bitset< unsigned int >::empt);

	k=0;
	for(j=0; j<index2choice.size(); ++j) {
	  if(index2choice[j].size() > 2) {
	    
	    //	    cout << "\t" << (k+examples[i][j]) << endl;

	    binExamples[i].insert( k+examples[i][j] );
	    k += index2choice[j].size();
	  } else if(index2choice[j].size() == 2) {
	    //	    cout << examples[i][j] << " ";
	    if(examples[i][j]) {
 	      binExamples[i].insert(k);
	      //	      cout << "\t*" << (k) << endl;
	    }
	    ++k;
	  }
	}
	if(examples[i][index2choice.size()] == classe) 
	  positive_example.push_back( binExamples[i] );
	else 
	  negative_example.push_back( binExamples[i] );
      }

    //exit(0);
    
  }


  void computeXor( )
  {

    //    MinHittingSet mhs(nBinFeatures, negative_example, positive_example);
    int k, i, j, n=positive_example.size(), m=negative_example.size();
    discrepancies = new Bitset< unsigned int >[n*m];
    exampleIndex = new int[n*m];
    columns = new vector<int>[nBinFeatures];
    for(k=n*m; k;)       
      discrepancies[--k].init(0, nBinFeatures-1, Bitset< unsigned int >::empt);
    for(i=0; i<n; ++i)
      for(j=0; j<m; ++j) {
	discrepancies[k].clear();
	discrepancies[k].unionWith(positive_example[i]);
	discrepancies[k].xorWith(negative_example[j]);
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
	j = bit;
	columns[j].push_back(i);
      } while( ++bit );


    }
    //    cout << endl;


  }


  void getData( const char *filename )
  {
    int i,j,k;
    char buf[1024];

    ifstream datafile( filename );  
    assert( datafile.good() );
  
    while( true ) {
      datafile.getline( buf, 1024 );

      // cout << endl << buf << endl;

      if( datafile.eof() || !datafile.good() ) break;

      //  cout << "process" << endl;

      // new example
      vector< int > ex;
      examples.push_back( ex );

      j=0;
      k=0;
      for( i=0; buf[i] != '\0'; ++i ) {
	if( buf[i] == ',' ) { 
	  // new feature?
	  if( index2choice.size() <= (unsigned int)k ) {
	    map<string, int>  valuations_map;
	    vector<string> valuations_vector;
	    choice2index.push_back( valuations_map );
	    index2choice.push_back( valuations_vector );
	  }
	 
	  buf[i] = '\0';
	  string valuation(&buf[j]);

	  // new valuation?
	  if( choice2index[k].find( valuation ) == choice2index[k].end() ) {
	    choice2index[k][valuation] = index2choice[k].size();
	    index2choice[k].push_back( valuation );
	  }
	  // build example
	  examples.back().push_back( choice2index[k][valuation] );

	  j=i+1;	
	  ++k;
	}	
      }

//       cout << endl;
//       for(int feat=0; feat<index2choice.size(); ++feat) {
// 	cout << "feature_" << feat << ": ";
// 	for(int vidx=0; vidx<index2choice[feat].size(); ++vidx) {
// 	  cout << index2choice[feat][vidx] << " ";
// 	  assert(choice2index[feat][index2choice[feat][vidx]] == vidx);
// 	}
// 	cout << endl;
//       }


      if(j) {
	string theclass(&buf[j]);

	// new class?
	if( class2index.find( theclass ) == class2index.end() ) {
	  class2index[theclass] = index2class.size();
	  index2class.push_back( theclass );
	}

	examples.back().push_back( class2index[theclass] );
      
// 	for(int classe=0; classe<index2class.size(); ++classe)
// 	    cout << index2class[classe] << " ";
// 	cout << endl << endl;

      } else {
	examples.pop_back();
      }
    }
  }



  void printData( )
  {
    cout << endl << "print" << endl;

    unsigned int i, j;
    cout << "classes: " << index2class.size() ;
    for(i=0; i<index2class.size(); ++i)
      cout << " " << index2class[i];
    cout << endl;

    cout << "features: " << index2choice.size() << endl;
  
    for(i=0; i<index2choice.size(); ++i)
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

    cout << "examples: " << examples.size() << endl;
    for(i=0; i<examples.size(); ++i)
      {
	for(j=0; j<index2choice.size(); ++j)
	  {
	    cout << index2choice[j][examples[i][j]] << " ";
	  }
	cout << " -> " << index2class[examples[i].back()] << endl;
      }
    cout << endl;

    
    cout << endl << setw(6) << negative_example.size() << " "
	 << "negative examples:" << endl;
    for( unsigned int i=0; i<positive_example.size(); ++i ) {
      cout << setw(6) << i << " ";
      negative_example[i].print( cout );
      cout << endl;
    }
    
    cout << endl << setw(6) << positive_example.size() << " "
	 << "positive examples:" << endl;
    for( unsigned int i=0; i<positive_example.size(); ++i ) {
      cout << setw(6) << i << " ";
      positive_example[i].print( cout );
      cout << endl;
    }
    cout << endl;
  }
 

  void satEncode( const int depth )
  {
    unsigned int i, j, k, l, node, incr, nLeaves = (1 << depth);
    unsigned int nNodes = nLeaves-1;
    int **path = new int*[nLeaves];
    for(i=0; i<nLeaves; ++i)
      path[i] = new int[depth+1];

    /// pre-compute all paths
    incr = nLeaves;
    node = 0;
    for(i=0; i<=(unsigned int)depth; ++i)
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
    int nClauses = nNodes + 
      (nBinFeatures * (nBinFeatures-1) / 2) * nNodes + 
      nLeaves * nDiscrepancies;
    cout << "p cnf " << (nNodes*nBinFeatures) << " " << nClauses << endl; 

    /// domain clauses
    cerr << depth << endl;
    for(i=0; i<nNodes; ++i)
      {
	for(j=0; j<nBinFeatures; ++j)
	  {
	    cout << i*nBinFeatures+j+1 << " ";

	    string value;
	    int feature = getFeature(j, value);
	    cerr << i*nBinFeatures+j+1 << " " << i 
	      //<< " example[" 
		 << " " 
		 << feature 
	      //<< "] == " 
		 << " " 
		 << value << endl;

	  }
	cout << "0\n";
	--nClauses;
	for(j=1; j<nBinFeatures; ++j)
	  for(k=0; k<j; ++k) {
	    cout << ((int)(-(i*nBinFeatures+j+1))) << " " 
		 << ((int)(-(i*nBinFeatures+k+1))) << " 0\n";
	    --nClauses;
	  }
      }

    for(i=0; i<(unsigned int)nDiscrepancies; ++i) {

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
	    for(l=0; l<(unsigned int)depth; ++l) {
	      //cout << "\t";
	      cout << (path[j][l]*nBinFeatures+k+1) << " ";
	      //cout << " diff" << endl;
	    }
	  else
	    for(l=0; l<(unsigned int)depth; ++l) {
// 	      cout << "\t";
// 	      cout << path[j][l] << " " << (depth-l-1) << " ";

	      int largestLeftLeaf = path[j][l]*(1 << (depth-l))+3*(1 << (depth-l-1))-2;

// 	      cout << largestLeftLeaf << ": " << j << " is on the " 
// 		   << ((nNodes+j) <= largestLeftLeaf ? "left" : "right") << endl;
	      
	      if( ((nNodes+j) <= (unsigned int)largestLeftLeaf) != goToLeft )
		cout << (path[j][l]*nBinFeatures+k+1) << " ";
	    }
	}
	cout << "0\n";
	--nClauses;
      }
    }
    

    assert ( !nClauses );

  }


};




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

  MinHittingSet(DataSet& ds) {
    //init( ds.nBinFeatures, ds.discrepancies, ds.nDiscrepancies );
    init( ds.nBinFeatures, ds.columns, ds.nDiscrepancies );
  }

  void init(const int nFeatures, 
	    vector<int> *columns,
	    const int nDiscrepancies) {
    {
      int i, j, k;
      
      epsilon = 0.000001;
      nColumns=nFeatures;
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
      matbeg  = (int *)    malloc (nRows    * sizeof(int));
      matcnt  = (int *)    malloc (nRows    * sizeof(int));   
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
      } else {
      	fprintf (stderr, 
      		 "CPlex environment open.\n");
      }
      
      status = CPXsetintparam (the_env, CPX_PARAM_SCRIND, CPX_ON);
      if ( status ) {
      	fprintf (stderr, 
      		 "Failure to turn on screen indicator, error %d.\n", status);
      } else {
      	fprintf (stderr, 
      		 "Param initialised.\n");
      }
      
      the_lp = CPXcreateprob (the_env, &status, "bidule");
      if ( the_lp == NULL ) {
      	fprintf (stderr, "Failed to create problem.\n");
      } else {
      	fprintf (stderr, 
      		 "Problem created.\n");
      }

#endif
      
      //////////////// CREATING THE LP /////////////////

      
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


      assert( j == k );

#ifdef _CPLEX

      status = CPXcheckcopylp(the_env, the_lp, nColumns, nRows, objsen, 
			      obj, rhs, sense, matbeg, matcnt, matind, 
			      matval, zlb, zub, NULL);

      status = CPXcopylp(the_env, the_lp, nColumns, nRows, objsen, 
			 obj, rhs, sense, matbeg, matcnt, matind, 
			 matval, zlb, zub, NULL);

      if ( status ) {
      	fprintf (stderr, 
      		 "Failure to copy data, error %d.\n", status);
      } else {
      	fprintf (stderr, 
      		 "Data copied.\n");
      }

#endif

    }
  }


  void print()
  {

#ifdef _CPLEX

    status = CPXwriteprob(the_env, the_lp, "myprob.lp", NULL);

      if ( status ) {
      	fprintf (stderr, 
      		 "Failure to write problem to file.\n", status);
      } else {
      	fprintf (stderr, 
      		 "Problem written.\n");
      }

#endif

  }

  double getLpBound()
  {

#ifdef _CPLEX

    status = CPXlpopt (the_env, the_lp);
    
    if ( status ) {
      fprintf (stderr, 
	       "Failure to solve linear relaxation.\n", status);
    } else {
      fprintf (stderr, 
	       "Linear relaxation solved.\n");
    }

    status = CPXsolution (the_env, the_lp, &lpstat, &objval, x, pi, slack, dj);

    if ( status ) {
      fprintf (stderr, 
	       "Failure to access solution.\n", status);
    } else {
      fprintf (stderr, 
	       "Solution loaded.\n");
    }

#endif

    return objval;

  }

  double solve()
  {

#ifdef _CPLEX

    status = CPXcopyctype (the_env, the_lp, ctype);
    
    if ( status ) {
      fprintf (stderr, 
	       "Failure to set integer vars.\n", status);
    } else {
      fprintf (stderr, 
	       "Integrality constraint.\n");
    }

    status = CPXmipopt (the_env, the_lp);

    if ( status ) {
      fprintf (stderr, 
	       "Failure to solve the problem.\n", status);
    } else {
      fprintf (stderr, 
	       "Problem solved.\n");
    }

    status = CPXgetbestobjval (the_env, the_lp, &objval);

    if ( status ) {
      fprintf (stderr, 
	       "Failure to access objective value.\n", status);
    } else {
      fprintf (stderr, 
	       "Objective value loaded.\n");
    }

    status = CPXgetx (the_env, the_lp, x, 0, CPXgetnumcols(the_env, the_lp)-1);

    if ( status ) {
      fprintf (stderr, 
	       "Failure to access solution.\n", status);
    } else {
      fprintf (stderr, 
	       "Solution loaded.\n");
    }

#endif

    return objval;

  }

};




class ExampleStruct {

public:

  static const int left=0;
  static const int right=1;
  static const int positive=1;
  static const int negative=0;
  
  ReversibleIntList  *Example[2];
  
  VariableInt      **    Feature;
  unsigned int            nNodes;
  unsigned int         nFeatures;
  unsigned int      nExamples[2]; 
  vector< Bitset< unsigned int> > examples[2];


  /********************************************************************
   * data structure that keeps, for every node, the number of positive 
   * and negative examples that would be on the left or right handside
   ********************************************************************/

  // for each node 'n', feature 'f', polarity 'p' and direction 'd'
  // how many examples of polarity 'p' would be routed to the direction 'd' 
  // at node 'n' if it was assigned with feature 'f' 
  int **remainder[2][2];

public:


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


  /**@name Constructors*/
  //@{
  //ExampleStruct(const int ln, int **pe, int np, int **ne, int nn, int nf)
  ExampleStruct(const int ln, const int nf, vector< Bitset< unsigned int > >& pe, vector< Bitset< unsigned int > >& ne)
  {
    nExamples[positive] = pe.size(); 
    nExamples[negative] = ne.size();
    examples[positive] = pe; 
    examples[negative] = ne;
    nFeatures = nf;
    nNodes = ln;
    Example[positive] = new ReversibleIntList[nNodes];
    Example[negative] = new ReversibleIntList[nNodes];


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
	for( unsigned int node = 0; node < nNodes; ++node )
	  remainder[polarity][direction][node] = new int[nFeatures];
      }
  }

  void init(Solver *s, VariableInt **feat) {
    Feature = feat;
    s->binds( Example[positive][0] ); 
    s->binds( Example[negative][0] ); 
    Example[positive][0].setValue( 0, nExamples[positive]-1, nExamples[positive] );
    Example[negative][0].setValue( 0, nExamples[negative]-1, nExamples[negative] );
    checkExamples(0,0);
    for(unsigned int i=1; i<nNodes; ++i) {
      s->binds( Example[positive][i] ); 
      s->binds( Example[negative][i] ); 
      Example[positive][i].setValue( 0, nExamples[positive]-1, 0 );
      Example[negative][i].setValue( 0, nExamples[negative]-1, 0 );
      checkExamples(i,0);
    }

    
  }

  virtual ~ExampleStruct()
  {
    delete [] Example[negative];
    delete [] Example[positive];
    for( int polarity = negative; polarity <= positive; ++polarity ) {
      for( int direction = left; direction <= right; ++direction ) {
	for( unsigned int node = 0; node < nNodes; ++node )
	  delete [] remainder[polarity][direction][node];
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
      for( unsigned int feature = 0; feature < nFeatures; ++feature ) 
	if( Feature[node]->contain( feature ) ) {
	  cout << "\tleft: " << (feature) << ": "
	       << (remainder[positive][left][node][feature]) << "/"
	       << (remainder[negative][left][node][feature]) 
	       << "  .   right: "
	       << (remainder[positive][right][node][feature]) << "/"
	       << (remainder[negative][right][node][feature]) 
	       << (isRelevant(node, feature) ? (decides(node, feature) ? "*" : " " ) : "irrelevant" ) 
	       << endl;	
	} 
      cout << endl;
    }
  }

  void printExamples( int node ) 
  {
    for( int polarity=negative; polarity<=positive; ++polarity ) {
      cout << (polarity ? "positive:" : "negative:" ) ;
      for( unsigned int i=0; i<Example[polarity][node].size; ++i ) {
	cout << " ";
	cout << (Example[polarity][node][i]);
      }
      cout << endl;
    }
  }





  inline int needChild( const int direction, const int node, const int feature )
  {
    return ( remainder[negative][direction][node][feature] && remainder[positive][direction][node][feature] );
  }
  
  inline int isRelevant( const int node, const int feature )
  {
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
    double npos = remainder[positive][direction][node][feature], nneg = remainder[negative][direction][node][feature], tot=npos+nneg;
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
	std::fill( remainder[polarity][direction][node], remainder[polarity][direction][node]+nFeatures, 0 );

    DomainIterator *feature = Feature[node]->begin();
    do {
      k = *feature;
      if(k<(int)nFeatures) {
	for( polarity = negative; polarity <= positive; ++polarity ) {
	  n = Example[polarity][node].size;
	  for(i=0; i<n; ++i) {	  
	    j = Example[polarity][node][i];
	    ++remainder[polarity][examples[polarity][j][k]][node][k];	  
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
  ConstraintLearning(Solver *s, VariableInt **x, const int ln, 
		     ExampleStruct *es, int nf)
    : Constraint(s, x, ln, Constraint::VALUETRIGGER)
  {
    nFeatures = nf;
    nNodes = (ln / 5);
    Parent = scope;
    Child[left]  = scope+nNodes;
    Child[right] = scope+2*nNodes;
    NumChild     = scope+3*nNodes;
    Feature      = scope+4*nNodes;
    ES = es;

    for(unsigned int i=0; i<nNodes; ++i) {
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
	if( (unsigned int)i<nNodes && scope[i]->value() == i ) {
	  ES->updateEntropy(i);
	  consistent = ( pruneFeature( i ) );
	} else if( (unsigned int)i>=4*nNodes )
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

    bool consistent = true;
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

    int consistent = true, direction, polarity, i, j, n, child[2] = {-1,-1};
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
	  direction = ES->examples[polarity][j][feature];
	  if( child[direction] >= 0 ) {
	    ES->Example[polarity][child[direction]].insert( j );
	  }
	}
      }
      
      ES->checkExamples(node, 2);
      if(child[left] >= 0)
	ES->checkExamples(child[left], 2);
      if(child[right] >= 0)
	ES->checkExamples(child[right], 2);
      
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

    if((unsigned int)changedIdx < 3*nNodes) {

    } else if((unsigned int)changedIdx < 4*nNodes) {
      int i, nchilds = 0;
      for(i=0; (unsigned int)i<nNodes; ++i)
	nchilds += NumChild[i]->min();
      for(i=1; consistent && i<=nchilds; ++i)
	consistent = Feature[i]->remove( nFeatures );
    } else {
      consistent = propagateFeature(node, value);
    }

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



/**********************************************
 * Learning Constraint Declaration
 **********************************************/
class BuildObjectLearning : public BuildObjectConstraint {

 public:
  
  int   n_feat;
  ConstraintLearning *conptr_;
  ExampleStruct *exptr_;
  
  BuildObjectLearning( 
		      ExampleStruct *es,
		      int nf )
  {
    n_feat = nf;
    exptr_ = es;
  }

  virtual void build( Solver *s, VariableInt **tmp, 
		      BuildObjectPredicate *pred )
  {
    conptr_ = new ConstraintLearning( s, tmp, pred->arity, 
				      exptr_, n_feat );
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
  int   n_feat;
 
  Learning( VarArray& x, VarArray& y, VarArray& z, 
	    VarArray& t, VarArray& r, 
	    ExampleStruct *es,
	    const int nFeatures )
      {
	n_feat = nFeatures;

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
	var_ptr_ = new BuildObjectPredicate( scope, 5*n, 0, -1, 
					     new BuildObjectLearning(es, n_feat ),
					     NULL );					     
      }  
};






class ValSelectorDtree : public ValSelector
{
public:

  ExampleStruct *ES;
  int node;
  int val;

  ValSelectorDtree( VariableInt *x, ExampleStruct *es, const int n )
    : ValSelector(x) 
  { 
    ES = es;
    node=n; 
  }

  /**@name Utils*/
  //@{ 
  inline int getBest() 
  { 
    if( _X->contain(ES->nFeatures) )
      val = ES->nFeatures;
    else {
      double mig = -1, ig;
      int feature = _X->min();
      val = feature;
      do {
	ig = ES->informationGain(node, feature);
	if( ig > mig ) { 
	  mig =  ig;
	  val = feature;
      }
      } while( _X->setNext( feature ) );
    }
    return val; 
  }
  inline int getBest() const 
  { 
    int val_=0;
    if( _X->contain(ES->nFeatures) )
      val_ = ES->nFeatures;
    else {
      double mig = -1, ig;
      int feature = _X->min();
      val_ = feature;
      do {
	ig = ES->informationGain(node, feature);
	if( ig > mig ) { 
	  mig =  ig;
	  val_ = feature;
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
  
  int *openNode;
  int nOpen;
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
	  ExampleStruct *es) : DVO(s) 
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
    ES = es;

    for(int i=0; i<nNodes; ++i)
      Feature[i]->branch = new ValSelectorDtree( Feature[i], ES, i );

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
      //ES->printInfo( i, true );
      //cout << endl;      
      //ES->printExamples(i);
      
      encore = !(NotUsed[i]->max());

    }
#endif
 
    // find an open node, and assign a feature.
    nOpen = 0;
    for(int i=0; i<nNodes; ++i) 
      if( !Feature[i]->isGround() && Parent[i]->isGround() )
	openNode[nOpen++] = i;

    if( !nOpen ) {
      return *first;
    }

    return Feature[openNode[0]];
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

  DecisionTreeStrategy( const int n, 
			VariableInt **p,
			VariableInt **l,
			VariableInt **r,
			VariableInt **k,
			VariableInt **f,
			VariableInt **u,
			ExampleStruct *es
			) { 
    nNodes=n; 
    tmpPrnt=p;
    tmpLeft=l;
    tmpRigh=r;
    tmpNumC=k;
    tmpFeat=f;
    tmpUsed=u;
    ES = es;
  }

  /**@name Utils*/
  //@{
  DVO* extract( Solver* s ) {    
    return new DTreeVO( s, nNodes, tmpPrnt, tmpLeft, tmpRigh, tmpNumC, tmpFeat, tmpUsed, ES );
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


void printTree( int lvl,
		int x,
		DataSet& ds,
		int *solution,
		int n,
		ExampleStruct *es );

int main(int argc, const char *argv[])
{ 

  //removeColumn( argv[1] );
  //exit(0);

  // get the data
  int i, j;
  DataSet ds;
  if( argc <= 1 ) {
    cerr << "Need an argument" << endl;
    exit( 0 );
  } else
    ds.getData( argv[1] );
  ds.binaryEncode( 0 );

  //ds.printData();   
  //cout << endl;

  if( argc > 3 ) {
    ds.computeXor();
    ds.satEncode(atoi(argv[2]));
  } else {

    int nNodes = (argc > 2 ? atoi( argv[2] ) : 10 );

#ifdef _CPLEX
   MinHittingSet mhs(ds);
   mhs.print();
//    //cout << "solve:" << endl;
//    double result = mhs.getLpBound();
//    cout << result << endl;

   double result = mhs.solve();
   cout << result << endl;

   for(i=0; i<ds.nBinFeatures; ++i)
 //     //cout << (mhs.x[i] > 0.001 ? (mhs.x[i] < 0.99 ? 0.5 : 1) : 0 ) << " ";
     //     cout << (mhs.x[i] > 0.001 ? (mhs.x[i] < 0.99 ? mhs.x[i] : 1) : 0 ) << " ";
 //   cout << endl;

   exit(0);
#endif

  // create the csp
  CSP model;
  
  VarArray Parent( nNodes, 0, nNodes-1 );
  VarArray LeftChild( nNodes, 0, nNodes-1 );
  VarArray RightChild( nNodes, 0, nNodes-1 );
  VarArray NumChild( nNodes, 0, 2 );
  VarArray Feature( nNodes, 0, ds.nBinFeatures );
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
	model.add( Descendant[i*nNodes+j] <= (Feature[i] != Feature[j]) );
	model.add( (NotUsed[i] + Children[j][i]) < 2 );
      }

  // The array 'Used' indicates which nodes will indeed represent a feature check.
  for( i=0; i<nNodes; ++i ) {
    model.add( (NotUsed[i] + Descendant[i]) < 2 );
    model.add( NotUsed[i] == (Feature[i] == ds.nBinFeatures) );
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
  ExampleStruct *ES = new ExampleStruct(nNodes,
 					ds.nBinFeatures,
					ds.positive_example, 
 					ds.negative_example);
  model.add( Learning(Parent, LeftChild, RightChild, NumChild, Feature, ES, ds.nBinFeatures) );
  //////////////////////////////////////////////////////////////
  //////////// END OF THE CLASSIFICATION PART //////////////////
  //////////////////////////////////////////////////////////////
  
  
//////////////////////////////////////////////////////////////
//////////// THE OBJECTIVE PART //////////////////////////////
//////////////////////////////////////////////////////////////
  //The objective is to minimise the number of used nodes.

  Variable Obj = Sum( NotUsed );
  model.add( Obj + Sum( NumChild ) == (nNodes-1) );

  //  model.add( Maximise( Obj ) );
//////////////////////////////////////////////////////////////
//////////// END OF THE OBJECTIVE PART ///////////////////////
//////////////////////////////////////////////////////////////




   VarArray searchVars;
   searchVars.add( Feature );
//   searchVars.add(Parent);
   searchVars.add(LeftChild);
   searchVars.add(RightChild);
//   Solver s( model, searchVars );
//   Lexicographic heuristic;
//   s.add(heuristic);
//   s.FIND_ALL = -1;

  
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
				 ES );
  s.add( strategy );
  //////////////////////////////////////////////////////////////
  //////////// END OF THE HEURISTIC PART ///////////////////////
  //////////////////////////////////////////////////////////////
  

  //  s.setVerbosity(4);
  //s.setTimeLimit(1);
  s.solve();

  cout << "#include \"tree.h\"" << endl << endl;

  cout << "string classify( string *example ) {" << endl;

  printTree( 1, 0, ds, s.solution, nNodes, ES );

  cout << "}" << endl << endl
       << "int main(int argc, const char *argv[]) {" << endl
       << "if(argc>1) checkTree( argv[1] );" << endl
       << "}" << endl << endl;
  }

}




void printTree( //int dir, 
	       int lvl,
	       int node,
	       DataSet& ds,
	       int *solution,
	       int n,
	       ExampleStruct *es
		)
  
{

  int i, leftChild=node, rightChild=node, binFeature = solution[node];
  if( (unsigned int)binFeature < ds.nBinFeatures ) {

    if( node < n-1 ) {
      leftChild = solution[n+node];
      rightChild = solution[2*n-1+node];
    }

//     cout << binFeature << endl;
//     int oldt = 0;

    int value, feature=0, threshold=ds.index2choice[0].size();
    //if( leftChild || rightChild ) {
      // work out which feature we are talking about    
      while( binFeature >= threshold ) {
	//cout << "[" << oldt << ".." << (threshold-1) << "]" << endl;
	//oldt = threshold;
	threshold += ds.index2choice[++feature].size();	
      }
      threshold -= ds.index2choice[feature].size();
      value = (binFeature - threshold);
      //}


//     cout << (ds.index2choice.size()) << " " << feature << endl
// 	 << (ds.index2choice[feature].size()) << " " << value << endl;
    string val = (ds.index2choice[feature][value]);
    for(i=0; i<lvl; ++i) cout << " ";
    cout << "if( example[" << feature 
	 << "] == \"" << val << "\") {" << endl;
    //<< "] == " << value << ") {" << endl;
    
    if( rightChild > node ) {     
      printTree( lvl+1, rightChild, ds, solution, n, es );      
    } else {
      for(i=0; i<=lvl; ++i ) cout << " ";
      cout << "return \"" << (ds.index2class[(es->remainder[1][1][node][binFeature] == 0)]) << "\";" << endl; 
      //cout << "return " << (es->remainder[1][0][node][binFeature] ? "true;" : "false;") << endl;   
    }

    for(i=0; i<lvl; ++i ) cout << " ";
    cout << "} else {" << endl;    
    
    if( leftChild > node ) {
      printTree( lvl+1, leftChild, ds, solution, n, es );          
    } else {
      for(i=0; i<=lvl; ++i ) cout << " ";
      cout << "return \"" << (ds.index2class[(es->remainder[1][0][node][binFeature] == 0)]) << "\";" << endl;
      //cout << "return " << (es->remainder[1][0][node][binFeature] ? "true;" : "false;") << endl;
    }
    
    for(i=0; i<lvl; ++i ) cout << " ";
    cout << "}" << endl;
    
  }

}




// void printTree( DataSet& ds,
// 		VarArray& Feature, 
// 		VarArray& LeftChild,
// 		VarArray& RightChild )
// {
//   recPrintTree( 0, Feature, LeftChild, RightChild );
// }

// void recPrintTree( int x, 
// 		   DataSet& ds,
// 		   VarArray& Feature, 
// 		   VarArray& LeftChild,
// 		   VarArray& RightChild )
// {

//   binFeature = Example[x]->value();

//   // work out which feature we are talking about
//   int i, value, feature=0, threshold=index2choice[0].size();
//   while( binFeature > threshold ) threshold += index2choice[++feature].size();
//   threshold -= index2choice[feature].size();
//   value = (binFeature - threshold);

//   for(i=0; i<x; ++i )
//   cout << "if( example[" << feature 
//        << "] == " << value << ")"

// }



// void printTree( int dir, 
// 		int x,
// 		DataSet& ds,
// 		int *solution,
// 		int n )

// {
//   int i;
//   if(x<n) {

//     int binFeature = solution[x];

//     if( binFeature < ds.nBinFeatures ) {

//       // work out which feature we are talking about
//       int value, feature=0, threshold=ds.index2choice[0].size();
//       while( binFeature > threshold ) threshold += ds.index2choice[++feature].size();
//       threshold -= ds.index2choice[feature].size();
//       value = (binFeature - threshold);
      
//       for(i=0; i<x; ++i ) cout << " ";
//       cout << "if( example[" << feature 
// 	   << "] == " << value << ") {" << endl;
//       printTree( 0, (solution[2*n+x] == x ? n : solution[2*n+x]), ds, solution, n );
//       for(i=0; i<x; ++i ) cout << " ";
//       cout << "} else {" << endl;
//       printTree( 1, (solution[n+x] == x ? n : solution[n+x]), ds, solution, n );  
//       for(i=0; i<x; ++i ) cout << " ";
//       cout << "}" << endl;
      
//     } else {
//       for(i=0; i<x; ++i ) cout << " ";
//       cout << (dir ? "True" : "False") << endl;
//     }
//   } else  {
//     for(i=0; i<x; ++i ) cout << " ";
//     cout << (dir ? "True" : "False") << endl;
//   }
// }

