
/*
  Mistral is a constraint satisfaction and optimisation library
  Copyright (C) 2003-2005  Emmanuel Hebrard
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

  The author can be contacted electronically at 
  ehebrard@cse.unsw.edu.au.
*/

#include <mistral_gen.h>
#include <mistral_glo.h>
//#include <mod.h>


// Return the binomial Coefficient C(N,K).
unsigned int Mistral::binomialCoeff(int N, int K) 
{
  unsigned int res = 1;
  unsigned int x=1;
  for(int i=0;i<K;i++){
    res = res*N--;
    x = x*(i+1);
  }
  return res/x;
}

// return the int corresponding to a given combination
unsigned int Mistral::comb2int(int *cm, int K) 
{
  int i=K;
  unsigned int res=0;
  while(i--)
    res+=binomialCoeff(cm[i],i+1);
  return res;
}

int hamming(int * s1, int * s2, int n) 
{
  int dist=0;
  for(int i=0; i<n; i++)
    dist+=(s1[i] != s2[i]);
  return dist;
}

using namespace Mistral;


/**********************************************
 * Combination Iterator
 *********************************************/

void Combination_iterator::create(int n, int k) 
{
  N=n;
  comb = new int[k+2];
}

Combination_iterator::Combination_iterator(int n, int k) 
{
  create(n, k);
}

void Combination_iterator::init(int n, int k) 
{
  N=n;
  init(k);
}

void Combination_iterator::init(int k) 
{
  K=k;
  for(j=0; j<K; j++)
    comb[j] = j;
  comb[K]   = N;
  comb[K+1] = 0;
  j=0;    
}

Combination_iterator::~Combination_iterator() 
{
  delete [] comb;
}

bool Combination_iterator::next() 
{    
  j=0;
  while( comb[j]+1 == comb[j+1] )
    comb[j] = j++;     
  comb[j]++;
    return (j < K);
}

bool Combination_iterator::operator++() 
{    
  j=0;
  while( comb[j]+1 == comb[j+1] )
    comb[j] = j++;     
  comb[j]++;
    return (j < K);
}

void Combination_iterator::print(std::ostream& o) {
  for(int i=0; i<K; i++)
    o << comb[i] << " ";
}


// /**********************************************
//  * Timing Function
//  *********************************************/

// #include <time.h>
// #include <sys/times.h>
// #include <sys/resource.h>
// #include <unistd.h>


// /// Return the ellapsed cpu time
// double getRunTime() {
//   double df = 0;

// #ifdef _UNIX
//   struct tms usage;
//   static int clock_ticks = sysconf(_SC_CLK_TCK);
//   times(&usage);
//   df=((double)usage.tms_utime+(double)usage.tms_stime)/clock_ticks;
// #endif

//   return df;
// }

// unsigned long int getMemory() {
//   unsigned long int mem = 0;
// //   struct rusage usage;
// //   getrusage(RUSAGE_CHILDREN, &usage);
// //   mem = usage.ru_maxrss;
// //   return mem;

// #ifdef _LINUX
  
//   char buf[30];
//   snprintf(buf, 30, "/proc/%u/statm", (unsigned)getpid());
//   FILE* pf = fopen(buf, "r");
  
//   if (pf) {
//     unsigned int size; //       total program size
//     //unsigned resident;//   resident set size
//     //unsigned share;//      shared pages
//     //unsigned text;//       text (code)
//     //unsigned lib;//        library
//     //unsigned data;//       data/stack
//     //unsigned dt;//         dirty pages (unused in Linux 2.6)
//     fscanf(pf, "%u" /* %u %u %u %u %u"*/, &size/*, &resident, &share, &text, &lib, &data*/);
//     //DOMSGCAT(MSTATS, std::setprecision(4) << size / (1024.0) << "MB mem used");
//     mem = size;
//   }  
//   fclose(pf);
  
// #endif

//   return mem;
// }



// /**********************************************
//  * Command line parsing
//  *********************************************/

// void getCommandLine(const char** int_ident, 
// 		    int*    int_param,
// 		    int     int_nb,
// 		    const char** str_ident,
// 		    const char** str_param,
// 		    int     str_nb,
// 		    char**  argv,
// 		    int     argc)
// {
//   std::fill(int_param,int_param+int_nb,NOVAL);
//   std::fill(str_param,str_param+str_nb,"nil");
//   int i,j,k;
//   for(i=1; i<argc-1;) {
//     k=1;

//     j = str_nb;    
//     while( j-- ) if( !strcmp(str_ident[j],argv[i]) ) {
//       str_param[j] = argv[i+1];
//       k=0;
//       break;
//     }
//     if(k) {
//       j = int_nb;
//       while( j-- ) if( !strcmp(int_ident[j],argv[i]) ) {
// 	int_param[j] = atoi(argv[i+1]);
// 	break;
//       }
//     }
//     i+=2;
//   }
// }

/**********************************************
 * Random_ number generator
 **********************************************/

float Random_::ran2(long *idum) const
{
  int j;
  long k;
  static long idum2 = 123456789;
  static long iy = 0;
  static long iv[NTAB];
  float temp;
  
  if (*idum <= 0) {                             /* initialize */
    if (-(*idum) < 1)                           /* prevent idum == 0 */
      *idum = 1;
    else
      *idum = -(*idum);                         /* make idum positive */
    idum2 = (*idum);
    for (j = NTAB + 7; j >= 0; j--) {           /* load the shuffle table */
      k = (*idum) / IQ1;
      *idum = IA1 * (*idum - k*IQ1) - k*IR1;
      if (*idum < 0)
	*idum += IM1;
      if (j < NTAB)
	iv[j] = *idum;
    }
    iy = iv[0];
  }
  
  k = (*idum) / IQ1;
  *idum = IA1 * (*idum - k*IQ1) - k*IR1;
  if (*idum < 0)
    *idum += IM1;
  k = idum2/IQ2;
  idum2 = IA2 * (idum2 - k*IQ2) - k*IR2;
  if (idum2 < 0)
    idum2 += IM2;
  j = iy / NDIV;
  iy = iv[j] - idum2;
  iv[j] = *idum;
  if (iy < 1)
    iy += IMM1;
  if ((temp = AM * iy) > RNMX)
    return RNMX;                                /* avoid endpoint */
  else
    return temp;
}

/*************************************************
 * Random_ warehouse allocation problems Generator 
 *************************************************/

WarehouseGenerator::WarehouseGenerator(long S, int ns, int nw)
{
  Seed = S;
  nbStores = ns;
  nbWarehouses = nw;
  capacity = new int[nbWarehouses];
  supplyCost = new int*[nbStores];
  for(int i=0; i<nbStores; ++i)
    supplyCost[i] = new int[nbWarehouses];
}
/*
WarehouseGenerator::WarehouseGenerator(const std::string filename)
{
  infile.open(filename, std::ifstream::in);
  
  infile >> nbStores;      
  infile >> nbWarehouses;
  capacity = new int[nbWarehouses];
  supplyCost = new int*[nbStores];
  for(int i=0; i<nbStores; ++i)
    supplyCost[i] = new int[nbWarehouses];
  //infile.close();
}
*/
WarehouseGenerator::WarehouseGenerator(const char *filename)
{
  infile.open(filename, std::ifstream::in);
  
  infile >> nbStores;      
  infile >> nbWarehouses;
  capacity = new int[nbWarehouses];
  supplyCost = new int*[nbStores];
  for(int i=0; i<nbStores; ++i)
    supplyCost[i] = new int[nbWarehouses];
  //infile.close();
}

WarehouseGenerator::~WarehouseGenerator()
{
  delete [] capacity;
  for(int i=0; i<nbStores; ++i)
    delete [] supplyCost[i];
  delete [] supplyCost;
}

void WarehouseGenerator::gen(int mind, int maxd, int minc, int maxc)
{
  // capacity:
  for(int i=0; i<nbWarehouses; ++i) 
    capacity[i] = (int)((float)(maxc-minc) * generator.ran2(&Seed)) + minc;

  // supply costs:
  for(int i=0; i<nbStores; ++i)
    for(int j=0; j<nbWarehouses; ++j)
      supplyCost[i][j] = (int)((float)(maxd-mind) * generator.ran2(&Seed)) + mind;
}

void WarehouseGenerator::getFromFile()
{
  int i,j;
  //std::ifstream infile(filename);
  
  for(j=0; j<nbWarehouses; ++j)
    infile >> capacity[j];

  for(i=0; i<nbStores; ++i)
    for(j=0; j<nbWarehouses; ++j)
      infile >> supplyCost[i][j];

  infile.close();
}


/******************************************
 * URCSPGenerator
 ******************************************/

URCSPGenerator::URCSPGenerator(int S, int v, int d, int c, int a, int n) 
{
  Seed = S;
  initURBCSP( v, d, c, a, n );  
}

void URCSPGenerator::initURBCSP( int v, int d, int c, int a, int n ) {
  Var  = v;
  Dom  = d;
  Con  = c;
  Ari  = a;
  Ngd  = n;
  
  if (Var < 2) 
    std::cerr<<" ***Illegal number of variables: "<<Var<<std::endl;
  if (Dom < 2) 
    std::cerr<<" ***Illegal domain size:"<<Dom<<std::endl;
  if (Con < 0 || Con > comb(Var,Ari)) 
    std::cerr<<" ***Illegal number of constraints: "<<Con << std::endl;
  if (Ari < 2)
    std::cerr<<" ***Illegal arity: "<<Ari<<std::endl;
  if (Ngd < 0 || Ngd > ((int)(pow((double)Dom,Ari) - 1)))
    std::cerr<<" ***Illegal tightness: "<<Ngd<<std::endl;

  int i;
  PossibleCTs = comb(Var,Ari);
  CTarray = new int[PossibleCTs];
  PossibleNGs = 1;
  if( Ngd ) {
    for( i=0; i<Ari; ++i )
      PossibleNGs = PossibleNGs*Dom;
  }
  NGarray = new int[PossibleNGs];
  vars = new int[Ari+1];
  //vals = new int[Var];
  vals = new int[Ari];
  vars[0]=-1;
  uv=Var;
  unconnected_var = new int[Var];
  reInit();
}

URCSPGenerator::~URCSPGenerator()
{
  delete [] CTarray;
  delete [] NGarray;
  delete [] vars;
  delete [] vals;
  delete [] unconnected_var;
}

int URCSPGenerator::comb(int b, int d){
  int res = 1;
  int x=1;
  for(int i=0;i<d;i++){
    res = res*b--;
    x = x*(i+1);
  }
  return res/x;
}

int URCSPGenerator::i2comb(int& code, int ind, int base, int dim, int stop){ 
  int c=comb(base,dim); 
  if(ind == stop || code < c)       
    return ind;        
  else{     
    code=(code-c);        
    return i2comb(code, ind+1, base-1, dim,stop);
  }
}

void URCSPGenerator::reInit() 
{

  int i;
  for( i=0; i<PossibleCTs; ++i )
    CTarray[i]=i;
  for( i=0; i<Var; ++i )
    unconnected_var[i]=1;
  c = 0;
}

void URCSPGenerator::reInit( int v, int d, int c, int a, int n ) 
{

//   std::cout << (ENV.store.size()) << std::endl;
//   ENV.flush();
//   std::cout << (ENV.store.size()) << std::endl << std::endl;

  delete [] CTarray;
  delete [] NGarray;
  delete [] vars;
  delete [] vals;
  delete [] unconnected_var;
  initURBCSP( v, d, c, a, n );
}

bool URCSPGenerator::erateConstraint() 
{
  t = 0;
  if( c == Con ) return false;
  /* Choose a random_ number between c and PossibleCTs - 1, inclusive. */
  int qq = (int) (generator.ran2(&Seed) * (PossibleCTs - c));
  int i, r =  c + qq;
  /* Swap elements [c] and [r]. */
  selectedCT = CTarray[r];
  CTarray[r] = CTarray[c];
  CTarray[c] = selectedCT;
  /* Create the constraint. */
  for(i=0;i<Ari;i++) {
    vars[i+1] = i2comb(CTarray[c],vars[i]+1,Var-(vars[i]+1)-1,Ari-1-i,Var-Ari+i);
    if(unconnected_var[vars[i+1]]) {
      unconnected_var[vars[i+1]]=0;	
      --uv;
    }
  }
  /* Initialize the NGarray. */
  for (i=0; i<PossibleNGs; ++i)
    NGarray[i] = i;
  ++c;
  return true;
}


bool URCSPGenerator::erateNogood() {
  if( t == Ngd ) return false;
  int r =  t + (int) (generator.ran2(&Seed) * (PossibleNGs - t));
  selectedNG = NGarray[r];
  NGarray[r] = NGarray[t];
  NGarray[t] = selectedNG;
  /* Add the nogood to the constraint. */
  for(int i=0;i<Ari;i++) {
    int x = (int)(pow((double)Dom,(Ari-i-1)));	    
    vals[i] = NGarray[t] / x;
    NGarray[t] -= (vals[i]*x);
  }
  ++t;
  return true;
}  



/******************************************
 * SATGenerator
 ******************************************/

SATGenerator::SATGenerator(int S, int v, int d, int c, int a) 
{
  Seed = S;
  initSAT( v, d, c, a );  
}

void SATGenerator::initSAT( int v, int d, int c, int a ) {
  Var  = v;
  Dom  = d;
  Con  = c;
  Ari  = a;

  _c = 0;
  
  if (Var < 2) 
    std::cerr<<" ***Illegal number of variables: "<<Var<<std::endl;
  if (Dom < 2) 
    std::cerr<<" ***Illegal domain size:"<<Dom<<std::endl;
  if (Con < 0) 
    std::cerr<<" ***Illegal number of constraints: "<<Con << std::endl;
  if (Ari < 2)
    std::cerr<<" ***Illegal arity: "<<Ari<<std::endl;

  vars = new int[Var];
  int i = Var;
  while( i-- )
    vars[i] = i;
  vals = new int[Ari];
}

SATGenerator::~SATGenerator()
{
  delete [] vars;
  delete [] vals;
}

bool SATGenerator::erateClause() 
{
  int i, qq, r, aux;
  if( _c >= Con ) return false;
  for(i=0;i<Ari;i++) {
    qq = (int) (generator.ran2(&Seed) * (Var - i));    
    r =  i + qq;
    aux = vars[i];
    vars[i] = vars[r];
    vars[r] = aux;
    vals[i] = (int) (generator.ran2(&Seed) * Dom); 
  }
  ++_c;
  return true;
}


/******************************************
 * JSPGenerator
 ******************************************/


JSPGenerator::JSPGenerator(int Nj, int Nm) 
{

  optAStdDevLowerBound=1;
  optAStdDevUpperBound=10;

  optAMeanLowerBound=35;
  optAMeanUpperBound=0;
  optAAlpha=0.75;     
  
  optBStdDevLowerBound=1;
  optBStdDevUpperBound=10;
  
  optBMeanLowerBound=35;
  optBMeanUpperBound=0; // computed below
  optBAlpha=0.75;      

  
  JSPGenerator::GBL_RandNumSeed=0;

 optTaillardLowerBound=20;
 optTaillardUpperBound=10000;
 
 optTaillardLowerBoundKeyword="lowerBound";
 optTaillardUpperBoundKeyword="upperBound";

 optTaillardGaussianMean=3000;
 optTaillardGaussianSigma=1000;

 optDepSetupGenerateMode=true;
 
 
 optAStdDevLowerBoundKeyword="stdDevLowerBound";
 optAStdDevUpperBoundKeyword="stdDevUpperBound";
 optAMeanLowerBoundKeyword="meanLowerBound";
 optAAlphaKeyword="alpha";

 
 optBStdDevLowerBoundKeyword="stdDevLowerBound";
 optBStdDevUpperBoundKeyword="stdDevUpperBound";
 optBMeanLowerBoundKeyword="meanLowerBound";
 optBAlphaKeyword="alpha";
 
 optDepSetupGenerateModeKeyword="depSetupComputeMode";
 optDepSetupGenerateModeFixed="fixed";
 optDepSetupGenerateModeOpDependent="opdependent";



  Njobs = Nj;
  Nmach = Nm;
  Durations = new int*[Nmach];
  Machines  = new int*[Nmach];
  for(int k=0; k<Nmach; k++) {
    Durations[k] = new int[Njobs];
    Machines [k] = new int[Njobs];
  }
  JSPDurations = Durations;
  JSPMachines  = Machines;
}  

JSPGenerator::~JSPGenerator()
{
  for(int k=0; k<Nmach; k++) {
    delete[] JSPDurations[k];
    delete[] JSPMachines [k];
  }  
  delete[] JSPDurations;
  delete[] JSPMachines;
}

//static std::map<std::string,std::string> GBL_OptionalParameterValueMap;

void JSPGenerator::storeOptionalParameters(char **parms,int numParms)
{
  for(int i=0;i<numParms;i++)
    {
      char *option=parms[i];
      char *equalPos=strchr(option,'=');
      if(equalPos==0)
	{
	  std::cout << "Invalid option specified encountered: " << option << std::endl;
	  std::cout << "Option was ignored" << std::endl;
	}
      else
	{
	  // form the keyword=value components from the option std::string.
	  (*equalPos)='\0';
	  std::string key(option+1);     // skip the dash
	  std::string value(equalPos+1); // skip the equal sign itself
	  JSPGenerator::GBL_OptionalParameterValueMap[key]=value;
	}
    }
}

bool JSPGenerator::optionPresent(const std::string &option)
{
  std::map<std::string,std::string>::iterator iter;
  iter=JSPGenerator::GBL_OptionalParameterValueMap.find(option);
  return iter!=JSPGenerator::GBL_OptionalParameterValueMap.end();
}

int JSPGenerator::optionAsInteger(const std::string &option)
{
  assert(optionPresent(option));

  const std::string &value=JSPGenerator::GBL_OptionalParameterValueMap[option];
  const char *stringData=value.c_str();

  return atoi(stringData);
}

double JSPGenerator::optionAsDouble(const std::string &option)
{
  assert(optionPresent(option));

  const std::string &value=JSPGenerator::GBL_OptionalParameterValueMap[option];
  const char *stringData=value.c_str();

  return atof(stringData);
}

std::string JSPGenerator::optionAsString(const std::string &option)
{
  assert(optionPresent(option));

  return JSPGenerator::GBL_OptionalParameterValueMap[option];
}

////////////////////////////////////////////////////////////////////////////////////
// RANDOM_ NUMBER GENERATORS DERIVED FROM TAILLARD'S ORIGINAL PAPER
// AND NUMERICAL RECIPIES IN 'C'
////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// set the current random_ number seed
////////////////////////////////////////////////////////////

//static long GBL_RandNumSeed=0;

void JSPGenerator::setSeed(long newSeed)
{
  JSPGenerator::GBL_RandNumSeed=newSeed;
}

// generate a number from 0.0 (inclusive) to 1.0 (exclusive)
double JSPGenerator::unifZeroOne(void)
{
  static long m = 2147483647, a = 16807, b = 127773, c = 2836;
  double  value_0_1;              

  long k = JSPGenerator::GBL_RandNumSeed / b;
  JSPGenerator::GBL_RandNumSeed = a * (JSPGenerator::GBL_RandNumSeed % b) - k * c;
  if(JSPGenerator::GBL_RandNumSeed < 0) JSPGenerator::GBL_RandNumSeed = JSPGenerator::GBL_RandNumSeed + m;
  value_0_1 =  JSPGenerator::GBL_RandNumSeed / (double) m;

  assert(value_0_1>=0.0);
  assert(value_0_1<1.0);

  return value_0_1;
}

// generate a number from low to high, inclusive
int JSPGenerator::unifSupplied(int low, int high)
{
  assert(low<=high);

  double value_0_1=unifZeroOne();

  int theNum=low + int(floor(value_0_1 * (high - low + 1)));

  assert(theNum>=low);
  assert(theNum<=high);

  return theNum;
}

// generate a random_ number from a normal distribution
// with the supplied parameters.
double JSPGenerator::normal(double mean, double stdDev)
{
  // straight from Numerical Recipies in C.
  // (the Box-Muller method)
  static int iset=0;
  double fac,r,v1,v2;
  static double gset;
  
  if(iset==0)
    {
      do
	{
	  v1=2.0*unifZeroOne()-1.0;
	  v2=2.0*unifZeroOne()-1.0;
	  r=v1*v1+v2*v2;
	}
      while(r>=1.0);
      fac=sqrt(-2.0*log(r)/r);
      gset=v1*fac;
      iset=1;
      return (v2*fac)*stdDev+mean;
    }
  else
    {
      iset=0;
      return gset*stdDev+mean;
    }

  // UNREACHABLE
}

                      
//////////////////////////////////////////////////////////////////////////////////////
// ROUTINES FOR GENERATION MODE 'TAILLARD' (the standard problem taillard generator):
//////////////////////////////////////////////////////////////////////////////////////

// static int optTaillardLowerBound=8;
// static int optTaillardUpperBound=50;

// const char *optTaillardLowerBoundKeyword="lowerBound";
// const char *optTaillardUpperBoundKeyword="upperBound";

void JSPGenerator::generateTaillardJSP_real(int numJobs, int numMachines)
{
  for(int j=0;j<numJobs;j++)
    {
      for(int i=0;i<numMachines;i++)
	{
	  JSPDurations[i][j]=unifSupplied(optTaillardLowerBound,
					  optTaillardUpperBound);
	}
    }
}

void JSPGenerator::generateTaillardJSP(int numJobs,int numMachines)
{
  // process any command-line over-rides of the default
  // problem generation parameters.
  if(optionPresent(optTaillardLowerBoundKeyword))
    {
      optTaillardLowerBound=optionAsInteger(optTaillardLowerBoundKeyword);
    }
  if(optionPresent(optTaillardUpperBoundKeyword))
    {
      optTaillardUpperBound=optionAsInteger(optTaillardUpperBoundKeyword);
    }
  assert(optTaillardLowerBound<optTaillardUpperBound);

  // generate the problem
  generateTaillardJSP_real(numJobs,numMachines);
}

////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////
// ROUTINES FOR GENERATION MODE 'TAILLARD-GAUSSIAN':
// instead of taking operation durations from a uniform distribution, take them
// from a Gaussian with mean=50, std. dev. =16. this should still cover
// the same range.
//////////////////////////////////////////////////////////////////////////////////////

// static int optTaillardGaussianMean=50;
// static int optTaillardGaussianSigma=16;

void JSPGenerator::generateTaillardJSPGaussian_real(int numJobs, int numMachines)
{
  for(int i=0;i<numMachines;i++)
    {
      for(int j=0;j<numJobs;j++)
	{
	  double num=normal(optTaillardGaussianMean,
			    optTaillardGaussianSigma);
	  if(num<1.0)
	    {
	      num=1.0;
	    }

	  JSPDurations[i][j]=int(num);
	}
    }
}

void JSPGenerator::generateTaillardJSPGaussian(int numJobs,int numMachines)
{
  // currently not handling any command-line over-rides.

  // generate the problem
  generateTaillardJSPGaussian_real(numJobs,numMachines);
}

////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////
// ROUTINES FOR GENERATION MODE 'JOB-CORRELATED':
// all operations for a given job are taken from a normal
// distribution with a given mean and variance
////////////////////////////////////////////////////////////////////////////////////

// std::vector<int> optAJobOpMeans;     // mean of the normal distribution
// std::vector<int> optAJobOpStdDevs;   // std. deviation of a normal distribution

// int optAStdDevLowerBound=1;
// int optAStdDevUpperBound=10;

// int    optAMeanLowerBound=35;
// int    optAMeanUpperBound=0;
// double optAAlpha=0.75;        

// const char *optAStdDevLowerBoundKeyword="stdDevLowerBound";
// const char *optAStdDevUpperBoundKeyword="stdDevUpperBound";
// const char *optAMeanLowerBoundKeyword="meanLowerBound";
// const char *optAAlphaKeyword="alpha";

void JSPGenerator::genOptionAJobOpDurations(int numJobs,int numMachines)
{
  optAJobOpMeans   = std::vector<int>(numJobs,0);
  optAJobOpStdDevs = std::vector<int>(numJobs,0);

  // select the machine distribution standard deviations.
  int sumStdDevs=0;
  for(int i=0;i<numJobs;i++)
    {
      optAJobOpStdDevs[i]=unifSupplied(optAStdDevLowerBound,
				       optAStdDevUpperBound);
      sumStdDevs+=optAJobOpStdDevs[i];
    }

  // the majority (~95%) of the distribution is contained
  // within 2 standard deviations of the mean; if you
  // take this as the width of a normal distribution,
  // the sum of such widths is, to some level of approximation,
  // the minimum interval size which can be expected to produce
  // a minimum number of overlaps in the distributions.
  int intervalWidth=4*sumStdDevs;
  
  // the optAAlpha parameter dictates the degree of 
  // overlap likely in the final distributions; 
  optAMeanUpperBound=optAMeanLowerBound+(int)(optAAlpha*intervalWidth);
  
  // select the machine distribution means.
  for(int j=0;j<numJobs;j++)
    {
      int num=unifSupplied(optAMeanLowerBound,
			   optAMeanUpperBound);
      optAJobOpMeans[j]=num;
    }
}

void JSPGenerator::generateJobCorrelatedJSP_real(int numJobs, int numMachines)
{
  for(int i=0;i<numMachines;i++)
    {
      for(int j=0;j<numJobs;j++)
	{
	  double num=normal(optAJobOpMeans[j],
			    optAJobOpStdDevs[j]);
	  if(num<1.0)
	    {
	      num=1.0;
	    }
	  JSPDurations[i][j]=int(num);
	}
    }
}

void JSPGenerator::generateJobCorrelatedJSP(int numJobs,int numMachines)
{
  // process any command-line over-rides of the default
  // problem generation parameters.
  if(optionPresent(optAStdDevLowerBoundKeyword))
    {
      optAStdDevLowerBound=optionAsInteger(optAStdDevLowerBoundKeyword);
    }
  if(optionPresent(optAStdDevUpperBoundKeyword))
    {
      optAStdDevUpperBound=optionAsInteger(optAStdDevUpperBoundKeyword);
    }
  assert(optAStdDevLowerBound<optAStdDevUpperBound);
  if(optionPresent(optAMeanLowerBoundKeyword))
    {
      optAMeanLowerBound=optionAsInteger(optAMeanLowerBoundKeyword);
      assert(optAMeanLowerBound>0);
    }
  if(optionPresent(optAAlphaKeyword))
    {
      optAAlpha=optionAsDouble(optAAlphaKeyword);
      assert((optAAlpha<=1.0)&&(optAAlpha>0.0));
    }

  // generate the problem
  genOptionAJobOpDurations(numJobs,numMachines);
  generateJobCorrelatedJSP_real(numJobs,numMachines);
}

////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////
// ROUTINES FOR GENERATION MODE 'MACHINE-CORRELATED':
// all operations on given machine are taken from a normal
// distribution with a given mean and variance
////////////////////////////////////////////////////////////////////////////////////

// std::vector<int> optBMachineMeans;     // mean of the normal distribution
// std::vector<int> optBMachineStdDevs;   // std. deviation of a normal distribution

// int optBStdDevLowerBound=1;
// int optBStdDevUpperBound=10;

// int    optBMeanLowerBound=35;
// int    optBMeanUpperBound=0; // computed below
// double optBAlpha=0.75;        

// const char *optBStdDevLowerBoundKeyword="stdDevLowerBound";
// const char *optBStdDevUpperBoundKeyword="stdDevUpperBound";
// const char *optBMeanLowerBoundKeyword="meanLowerBound";
// const char *optBAlphaKeyword="alpha";

void JSPGenerator::genOptionBMachineDurations(int numJobs,int numMachines)
{
  optBMachineMeans   = std::vector<int>(numMachines,0);
  optBMachineStdDevs = std::vector<int>(numMachines,0);

  // select the machine distribution standard deviations.
  int sumStdDevs=0;
  for(int i=0;i<numMachines;i++)
    {
      optBMachineStdDevs[i]=unifSupplied(optBStdDevLowerBound,
					 optBStdDevUpperBound);
      sumStdDevs+=optBMachineStdDevs[i];
    }

  // the majority (~95%) of the distribution is contained
  // within 2 standard deviations of the mean; if you
  // take this as the width of a normal distribution,
  // the sum of such widths is, to some level of approximation,
  // the minimum interval size which can be expected to produce
  // a minimum number of overlaps in the distributions.
  int intervalWidth=4*sumStdDevs;
  
  // the optBAlpha parameter dictates the degree of 
  // overlap likely in the final distributions; 
  optBMeanUpperBound=optBMeanLowerBound+(int)(optBAlpha*intervalWidth);
  
  // select the machine distribution means.
  for(int j=0;j<numMachines;j++)
    {
      int num=unifSupplied(optBMeanLowerBound,
			   optBMeanUpperBound);
      optBMachineMeans[j]=num;
    }
}

void JSPGenerator::generateMachineCorrelatedJSP_real(int numJobs, int numMachines)
{
  for(int i=0;i<numMachines;i++)
    {
      for(int j=0;j<numJobs;j++)
	{
	  double num=normal(optBMachineMeans[i],
			    optBMachineStdDevs[i]);
	  if(num<1.0)
	    {
	      // a statistical anamoly can happen; in this
	      // case, make sure the number is above one.
	      num=1.0;
	      //	      std::cout << "function generateMachineCorrelatedJSP_real: num was : " << num << std::endl;
	      //	      std::cout << "mean  : " << optBMachineMeans[i] << std::endl;
	      //	      std::cout << "stddev: " << optBMachineStdDevs[i] << std::endl;
	      //	      assert(num>=1.0);
	    }

	  JSPDurations[i][j]=(int)num;
	}
    }
}

void JSPGenerator::generateMachineCorrelatedJSP(int numJobs,int numMachines)
{
  // process any command-line over-rides of the default
  // problem generation parameters.
  if(optionPresent(optBStdDevLowerBoundKeyword))
    {
      optBStdDevLowerBound=optionAsInteger(optBStdDevLowerBoundKeyword);
    }
  if(optionPresent(optBStdDevUpperBoundKeyword))
    {
      optBStdDevUpperBound=optionAsInteger(optBStdDevUpperBoundKeyword);
    }
  assert(optBStdDevLowerBound<optBStdDevUpperBound);
  if(optionPresent(optBMeanLowerBoundKeyword))
    {
      optBMeanLowerBound=optionAsInteger(optBMeanLowerBoundKeyword);
      assert(optBMeanLowerBound>0);
    }
  if(optionPresent(optBAlphaKeyword))
    {
      optBAlpha=optionAsDouble(optBAlphaKeyword);
      assert((optBAlpha<=1.0)&&(optBAlpha>0.0));
    }

  // generate the problem
  genOptionBMachineDurations(numJobs,numMachines);
  generateMachineCorrelatedJSP_real(numJobs,numMachines);
}

////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////
// ROUTINES FOR SETUP TIME MATRICES
// ONLY ONE OPTION CURRENTLY EXISTS:
// A - SETUP TIMES ON MACHINES DEPEND ON THE JOB AND PREVIOUS JOB ON THAT MACHINE
////////////////////////////////////////////////////////////////////////////////////

// true  - sample dependent setup times uniformly from an interval of fixed
//         width, independent of operation.
// false - compute dependent setup times by first sampling an integer x between
//         1 and 10, and then taking the setup time as the operation duration
//         multiplied by x/100. 
// static bool optDepSetupGenerateMode=true;

// const char *optDepSetupGenerateModeKeyword="depSetupComputeMode";
// const char *optDepSetupGenerateModeFixed="fixed";
// const char *optDepSetupGenerateModeOpDependent="opdependent";

void JSPGenerator::loadDependentSetupTimesGenerationOptions(void)
{
  if(optionPresent(optDepSetupGenerateModeKeyword))
    {
      if(optionAsString(optDepSetupGenerateModeKeyword)==optDepSetupGenerateModeFixed)
	{
	  optDepSetupGenerateMode=true;
	}
      else if(optionAsString(optDepSetupGenerateModeKeyword)==optDepSetupGenerateModeOpDependent)
	{
	  optDepSetupGenerateMode=false;
	}
      else
	{
	  std::cerr << "***Invalid value specified for 'depSetupComputeMode' keyword: " << optionAsString(optDepSetupGenerateModeKeyword) << std::endl;
	  std::cerr << "***Defaulting to 'fixed' mode" << std::endl;
	}
    }
}

// the 3-D matrix of operation-dependent setup times,
// indexed by :
// 0th dim.  - job of interest
// 1rst dim. - machine of iterest
// 2nd dim.  - the job that was there previously
//static std::vector<std::vector<std::vector<int> > > optDependentSetupTimes;

// job and machine are assumed to be 0-based.
int JSPGenerator::generateDependentSetupTime(int job, int machine)
{
  if(optDepSetupGenerateMode==true)
    {
      const int depSetupLowerBound=1;
      const int depSetupUpperBound=10;
      return unifSupplied(depSetupLowerBound,depSetupUpperBound);
    }
  else
    {
      const int pctLowerBound=1;
      const int pctUpperBound=10;
      int percentage=unifSupplied(pctLowerBound,pctUpperBound);
      return std::max(1,int(ceil(double(percentage)/100.0*double(JSPDurations[machine][job]))));
    }
  // UNREACHABLE
}

void JSPGenerator::generateDependentSetupTimes(int numJobs,int numMachines)
{
  // 'null' the existing setup times matrix.
  optDependentSetupTimes.erase(optDependentSetupTimes.begin(),
			       optDependentSetupTimes.end());

  for(int i=0;i<numJobs;i++)
    {
      std::vector<std::vector<int> > jobSubmatrix(numMachines);
      for(int j=0;j<numMachines;j++)
	{
	  std::vector<int> machineVector(numJobs,0);
	  for(int k=0;k<numJobs;k++)
	    {
	      if(i!=k)
		{
		  machineVector[k]=generateDependentSetupTime(i,j);
		}
	    }
	  jobSubmatrix[j]=machineVector;
	}
      optDependentSetupTimes.push_back(jobSubmatrix);
    }
}

void JSPGenerator::writeDependentSetupTimes(int numJobs,int numMachines)
{
  std::cout << std::endl;
  for(int i=0;i<numJobs;i++)
    {
      for(int j=0;j<numMachines;j++)
	{
	  for(int k=0;k<numJobs;k++)
	    {
	      if(i==k)
		{
		  std::cout << std::setw(3) << (int)0 << " ";
		}
	      else
		{
		  std::cout << std::setw(3) << optDependentSetupTimes[i][j][k] << " ";
		}
	    }
	  std::cout << std::endl;
	}
      std::cout << std::endl;
    }
  std::cout << std::endl;
}

// computes the minimum setup time for the input job
// on the input machine, over all previous jobs possibly
// placed on that machine. theJob and theMachine are assumed
// to be 0-based.
int JSPGenerator::minJobSetupTime(const int theJob,const int theMachine)
{
  // if the dependent setup times have not been computed, just return 0 -
  // otherwise, compute it from the matrix.
  if(optDependentSetupTimes.size()==0)
    {
      return 0;
    }
  else
    {
      assert(theJob>=0);
      assert(theMachine>=0);

      const std::vector<int> &subVector=optDependentSetupTimes[theJob][theMachine];
      int minSetup=NOVAL;
      for(int i=0;i<int(subVector.size());i++)
	{
	  if(i!=theJob)
	    {
	      minSetup=std::min(minSetup,subVector[i]);
	    }
	}
      return minSetup;
    }
  // UNREACHABLE
}

////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////
// ROUTINES COMMON TO ALL GENERATION MODES
////////////////////////////////////////////////////////////////////////////////////

// assign machines to each of the job operations - implements Taillard's method.
void JSPGenerator::computeOperationMachineAssignments(int numJobs,int numMachines)
{
  // default the operation machine to the operation #.
  for(int i=0;i<numJobs;i++)
    {
      for(int j=0;j<numMachines;j++)
	{
	  JSPMachines[j][i]=j;
	}
    }

  // randomly swap orderings to create the actual assignments.
  for(int i=0;i<numJobs;i++)
    {
      for(int j=0;j<numMachines;j++)
	{
	  int swapIndex=unifSupplied(j,numMachines-1);
	  int temp=JSPMachines[swapIndex][i];
	  JSPMachines[swapIndex][i]=JSPMachines[j][i];
	  JSPMachines[j][i]=temp;
	}
    }
}

// computes the index of the operation of the input job
// that is to be processed on the input machine - just
// computes a 1-d scan of a row of JSPMachines.
int JSPGenerator::operationIndex(int job,int machine,int numMachines)
{
  for(int i=0;i<numMachines;i++)
    {
      if(JSPMachines[i][job]==machine)
	{
	  return i;
	}
    }

  // to force failure - you should never get this far.
  assert(true);
  return 0;
}

// straight from Taillard's paper 'Benchmarks for basic scheduling problems'.
// currently does not include consideration of setup times.
long JSPGenerator::taillardLowerBound(int numJobs,int numMachines)
{
  std::vector<long> b(numMachines,MAXLNG); // Bi's
  std::vector<long> a(numMachines,MAXLNG); // Ai's
  std::vector<long> t(numMachines,0);        // Ti's

  for(int i=0;i<numMachines;i++)
    {
      // compute Bi's.
      long minSoFar=MAXLNG;
      for(int j=0;j<numJobs;j++)
	{
	  long thisMin=0;
	  for(int k=0;k<operationIndex(j,i,numMachines);k++)
	    {
	      thisMin+=JSPDurations[k][j];
	    }
	  minSoFar=std::min(minSoFar,thisMin);
	}
      b[i]=minSoFar;

      // compute Ai's.
      minSoFar=MAXLNG;
      for(int j=0;j<numJobs;j++)
	{
	  long thisMin=0;
	  for(int k=operationIndex(j,i,numMachines)+1;k<numMachines;k++)
	    {
	      thisMin+=JSPDurations[k][j];
	    }
	  minSoFar=std::min(minSoFar,thisMin);
	}
      a[i]=minSoFar;

      // compute Ti's.
      for(int j=0;j<numJobs;j++)
	{
	  t[i]+=JSPDurations[operationIndex(j,i,numMachines)][j];
	}
    }

  // compute the longest-duration job.
  std::vector<long> d(numJobs,0); // d_{j}'s
  long maxDur=0;
  for(int j=0;j<numJobs;j++)
    {
      long sum=0;
      for(int k=0;k<numMachines;k++)
	{
	  sum+=JSPDurations[k][j];
	}
      maxDur=std::max(maxDur,sum);
    }

  // compute max 'compressed' machine duration
  long maxCompTime=0;
  for(int k=0;k<numMachines;k++)
    {
      maxCompTime=std::max(maxCompTime,a[k]+b[k]+t[k]);
    }

  return std::max(maxCompTime,maxDur);
}

void JSPGenerator::writeProblem(int numJobs,int numMachines)
{
  //  std::cout << std::endl;
  //std::cout << std::setw(3) << numJobs << " " << std::setw(3) << numMachines << std::endl << std::endl;
  //  int act;
  for(int j=0;j<numJobs;j++)
    {
      for(int i=0;i<numMachines;i++) 
	{
	  //  std::cout << std::setw(3) << JSPMachines[i][j]+FIRMACIND << " " << std::setw(3) << JSPDurations[i][j] << " " ;
	  //j_for_res[i].push_back(JSPMachines[i][j]);
	  //std::cout<<11<<std::endl;
	  //act = numJobs*JSPMachines[i][j]+j;
	  //std::cout<<22<<std::endl;
	  //a_for_res[i].push_back(act);
	  //std::cout<<33<<" "<<JSPMachines[i][j]<<" "<<j<<" "<<JSPDurations[i][j]<<std::endl;
	  //activity[j][i]=JSPDurations[i][j];
	  //std::cout<<44<<std::endl;

	}
      //std::cout << std::endl;
    }
  //std::cout << std::endl;
}

void JSPGenerator::writeLowerBounds(int numJobs,int numMachines, int* mkp)
{
  long taillardLB=taillardLowerBound(numJobs,numMachines);
  *mkp = (int)taillardLB;
  //std::cout << "Lower bound: " << taillardLB << std::endl;
}

////////////////////////////////////////////////////////////////////////////////////

int JSPGenerator::jspgen( std::string optionStr, int numJobs, int numMachines,
			  int timeSeed, int machineSeed, int *mkp )
{
  //std::std::cout<<numJobs<<" "<<numMachines<<" "<<*mkp<<" "<<std::std::endl;

  if(numJobs<=0)
    {
      std::cout << "Illegal number of jobs specified: " << numJobs << std::endl;
      return 0;
    }
  else if(numMachines<=0)
    {
      std::cout << "Illegal number of machines specified: " << numMachines << std::endl;
      return 0;
    }
  else if(timeSeed<=0)
    {
      std::cerr << "Illegal time seed specified: " << timeSeed << std::endl;
      return 0;      
    }
  else if(machineSeed<=0)
    {
      std::cerr << "Illegal machine seed specified: " << machineSeed << std::endl;
      return 0;      
    }

  // PHASE 1 - GENERATE THE OPERATION DURATIONS ACCORDING TO 
  //           THE SPECIFIED GENERATION OPTION - THE 'TIME' SEED
  //           IS USED IN THIS PHASE. ALSO, SETUP TIMES, IF SPECIFIED,
  //           WILL BE GENERATED.
  setSeed(timeSeed);

  const char *optionTaillardKeyword="taillard";
  const char *optionTaillardGaussianKeyword="taillard-gaussian";
  const char *optionJobCorrelatedKeyword="job-correlated";
  const char *optionMachineCorrelatedKeyword="machine-correlated";

  if(optionStr==optionTaillardKeyword)
    {
      generateTaillardJSP(numJobs,numMachines);
    }
  else if(optionStr==optionTaillardGaussianKeyword)
    {
      generateTaillardJSPGaussian(numJobs,numMachines);
    }
  else if(optionStr==optionJobCorrelatedKeyword)
    {
      generateJobCorrelatedJSP(numJobs,numMachines);
    }
  else if(optionStr==optionMachineCorrelatedKeyword)
    {
      generateMachineCorrelatedJSP(numJobs,numMachines);
    }
  else
    {
      std::cout << "Unknown generation option " << optionStr << " specified" << std::endl;
      return 0;
    }

  const char *setupOptionKeyword="setup";
  const char *dependentSetupTimeKeyword="dependent";      

  if(optionPresent(setupOptionKeyword))
    {
      std::string setupOption(optionAsString(setupOptionKeyword));

      if(setupOption==dependentSetupTimeKeyword)
	{
	  loadDependentSetupTimesGenerationOptions();
	  generateDependentSetupTimes(numJobs,numMachines);
	}
      else
	{
	  std::cerr << "***Illegal setup option specified: " << setupOption << std::endl;
	  std::cerr << "***No setup times will be computed" << std::endl;
	}
    }

  // PHASE 2 - GENERATE THE OPERATION ORDERING FOR EACH JOB -
  //           THE 'MACHINE' SEED IS USED IN THIS PHASE. 
  setSeed(machineSeed);

  computeOperationMachineAssignments(numJobs,numMachines);

  // PHASE 3 - OUTPUT THE OPERATION DURATIONS/MACHINE ORDER, 
  //           SETUP TIMES, AND LOWER BOUND.
  writeProblem(numJobs,numMachines/*,activity,a_for_res,j_for_res*/);

  if(optionPresent(setupOptionKeyword))
    {
      writeDependentSetupTimes(numJobs,numMachines);
    }

  writeLowerBounds(numJobs,numMachines,mkp);

  return 1;
}


int JSPGenerator::jspgen( const char* filename, 
			  int *mkp )
{
  std::ifstream infile( filename, std::ifstream::in );

  int i,j;
  infile >> Njobs;
  infile >> Nmach;
  //Durations = new int*[Nmach];
  //Machines  = new int*[Nmach];
  for(i=0; i<Nmach; ++i) {
    //Machines [i] = new int[Njobs];
    //Durations[i] = new int[Njobs];
    for(j=0; j<Njobs; ++j) {
      infile >> JSPMachines[j][i];      
      infile >> JSPDurations[j][i];
      //std::cout << Machines[i][j] << " " << Durations[i][j] << " ";
    }
    //std::cout << std::endl;
  }

  infile.close();

  // PHASE 3 - OUTPUT THE OPERATION DURATIONS/MACHINE ORDER, 
  //           SETUP TIMES, AND LOWER BOUND.
  writeProblem(Njobs,Nmach/*,activity,a_for_res,j_for_res*/);

  writeLowerBounds(Njobs,Nmach,mkp);

  //std::cout << "****" << *mkp << "****" << std::endl; 
  
  

  return 1;
}

