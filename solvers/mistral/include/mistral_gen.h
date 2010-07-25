
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

/** \file gen.h
    \brief Header for Problem Generators.
*/

#ifndef __GENERATOR_H
#define __GENERATOR_H

#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <vector>
#include <utility>
#include <map>
#include <string>
#include <algorithm>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <mistral_dvo.h>


// /**********************************************
//  * Knuth's Random number generator (code from sp-1.4)
//  **********************************************/

// void usrand (unsigned seed);
// unsigned urand0 (void);
// unsigned urand (void);
// int randint(int upto);
// double randreal();

// #define MAX_URAND 0xFFFFFFFFL

/**********************************************
 * Other Random number generator
 **********************************************/

#define IM1   2147483563
#define IM2   2147483399
#define AM    (1.0/IM1)
#define IMM1  (IM1-1)
#define IA1   40014
#define IA2   40692
#define IQ1   53668
#define IQ2   52774
#define IR1   12211
#define IR2   3791
#define NTAB  32
#define NDIV  (1+IMM1/NTAB)
#define EPS   1.2e-7
#define RNMX  (1.0 - EPS)
#define FIRMACIND 0  /* 0,1: first machine index */ 


namespace Mistral {

  /// Random number generator.
  class Random_ {
  public:
    /** 
	Return a random floating point value between 0.0 and
	1.0 exclusive.  If idum is negative, a new series starts (and
	idum is made positive so that subsequent calls using an unchanged
	idum will continue in the same sequence). 
    */
    float ran2( long * ) const;
  };

  /// Generate random warehouse allocation problems. 
  /**
   *  The code is mine :)
   */
  class WarehouseGenerator{
  private:
    std::ifstream infile;
    Random_ generator;

  public: 
    long Seed;

    /**@name Constructors*/
    //@{
    /// 
    WarehouseGenerator(long, int, int);
    ///
    WarehouseGenerator(const char*);
    ///
    //WarehouseGenerator(const std::string);
    ///
    ~WarehouseGenerator();
    //@}

    /// The number os stores
    int nbStores;
    /// The number of warehouses
    int nbWarehouses;
    /// the matrix of costs
    int **supplyCost;
    /// the capacity of the warehouses
    int *capacity;  

    /// the generator
    void gen(int, int, int, int);

    /// the generator
    void getFromFile();

  };


  /// Generate random uniform CSPs. 
  /**
   *  The code is from Frost and Bessiere's generator
   */
  class URCSPGenerator{
  private:
    int comb(int b, int d);  
    int i2comb(int& code, int ind, int base, int dim, int stop);
  
    Random_ generator;

    int Var;
    int Dom;
    int Con; 
    int Ari; 
    int Ngd;


    int c;
    int t;
    int selectedNG;
    int selectedCT; 
    int PossibleCTs; 
    int PossibleNGs;

    int *CTarray;
    int *NGarray;

    int uv;
    int *unconnected_var; 

  public: 
    long Seed;

    /**@name Constructors*/
    //@{
    /// 
    URCSPGenerator(int S, int v, int d, int c, int a, int n);
    void initURBCSP( int v, int d, int c, int a, int n );
    ///
    ~URCSPGenerator();
    //@}

    /// An array of indices correponding to the constrained variables
    int *vars;
    /// An array corresponding to a nogood, such that vals[vars[i]] is the ith value of the nogood
    int *vals;
    /**
     *  Initialisation function, to be called after every CSP generation
     */
    void reInit();
    void reInit(int v, int d, int c, int a, int n);
    /**
     *  Generate a constraint scope and store it in vars
     */
    bool erateConstraint();
    /**
     *  Generate a nogood tuple and store it in vals
     */
    bool erateNogood();

  };


  /// Generate random uniform SAT instances. 
  class SATGenerator{
  private:
    Random_ generator;
    long Seed;

    int Var;
    int Dom;
    int Con; 
    int Ari; 

  public: 
    /**@name Constructors*/
    //@{
    /// 
    SATGenerator(int S, int v, int d, int c, int a);
    void initSAT( int v, int d, int c, int a );
    ///
    ~SATGenerator();
    //@}

    int _c;
    /// An array of indices correponding to the constrained variables
    int *vars;
    /// An array corresponding to a nogood, such that vals[vars[i]] is the ith value of the nogood
    int *vals;
    /**
     *  Generate a clause
     */
    bool erateClause();

  };

  /// Generate random job-shop scheduling problems. 
  /**
   *  The code is from Beck's generator
   */
  class JSPGenerator{
  
  public:
  
    std::map<std::string,std::string> GBL_OptionalParameterValueMap;

    long GBL_RandNumSeed;

    int optTaillardLowerBound;
    int optTaillardUpperBound;
  
    const char *optTaillardLowerBoundKeyword;
    const char *optTaillardUpperBoundKeyword;
  
    int optTaillardGaussianMean;
    int optTaillardGaussianSigma;
  
    std::vector<int> optAJobOpMeans;     // mean of the normal distribution
    std::vector<int> optAJobOpStdDevs;   // std. deviation of a normal distribution
  
    int optAStdDevLowerBound;
    int optAStdDevUpperBound;
  
    int    optAMeanLowerBound;
    int    optAMeanUpperBound;
    double optAAlpha;        
  
    const char *optAStdDevLowerBoundKeyword;
    const char *optAStdDevUpperBoundKeyword;
    const char *optAMeanLowerBoundKeyword;
    const char *optAAlphaKeyword;

  
    std::vector<int> optBMachineMeans;     // mean of the normal distribution
    std::vector<int> optBMachineStdDevs;   // std. deviation of a normal distribution

    int optBStdDevLowerBound;
    int optBStdDevUpperBound;

    int    optBMeanLowerBound;
    int    optBMeanUpperBound; // computed below
    double optBAlpha;        

    bool optDepSetupGenerateMode;

    const char *optBStdDevLowerBoundKeyword;
    const char *optBStdDevUpperBoundKeyword;
    const char *optBMeanLowerBoundKeyword;
    const char *optBAlphaKeyword;

    const char *optDepSetupGenerateModeKeyword;
    const char *optDepSetupGenerateModeFixed;
    const char *optDepSetupGenerateModeOpDependent;

    std::vector<std::vector<std::vector<int> > > optDependentSetupTimes;



    // INDEXED BY MACHINE/OPERATION FIRST, JOB SECOND
    int **JSPDurations; // JSPDurations[x][y] = D: 
    // the xth operation of the yth job has duration D
    int **JSPMachines;  // JSPMachines[x][y]  = R: 
    // the xth operation of the yth job need resource R


  public:
    /// Number of jobs
    int Njobs;
    /// Number of machines
    int Nmach;

    /// JSPDurations[x][y] = D: the xth operation of the yth job has duration D
    int **Durations; 
    /// JSPMachines[x][y]  = R: the xth operation of the yth job need resource R
    int **Machines; 
  
    /**
     *  Generate a job-shop scheduling problem and store it in 
     *  JSPDurations and JSPMachines
     */
    int jspgen( std::string optionStr, int numJobs, int numMachines,
		int timeSeed, int machineSeed, int *mkp );

    int jspgen( const char* fname, int *mkp );

    /**@name Constructors*/
    //@{
    ///
    JSPGenerator(int Nj, int Nm);
    ///
    ~JSPGenerator();
    //@}

    void storeOptionalParameters(char **parms,int numParms);
    bool optionPresent(const std::string &option);
    int optionAsInteger(const std::string &option);
    double optionAsDouble(const std::string &option);
    std::string optionAsString(const std::string &option);

    void setSeed(long newSeed);
    double unifZeroOne(void);
    int unifSupplied(int low, int high);
    double normal(double mean, double stdDev);

    void generateTaillardJSP_real(int numJobs, int numMachines);
    void generateTaillardJSP(int numJobs,int numMachines);

    void generateTaillardJSPGaussian_real(int numJobs, int numMachines);
    void generateTaillardJSPGaussian(int numJobs,int numMachines);

    void genOptionAJobOpDurations(int numJobs,int numMachines);
    void generateJobCorrelatedJSP_real(int numJobs, int numMachines);
    void generateJobCorrelatedJSP(int numJobs,int numMachines);

    void genOptionBMachineDurations(int numJobs,int numMachines);
    void generateMachineCorrelatedJSP_real(int numJobs, int numMachines);
    void generateMachineCorrelatedJSP(int numJobs,int numMachines);

    void loadDependentSetupTimesGenerationOptions(void);
    int generateDependentSetupTime(int job, int machine);
    void generateDependentSetupTimes(int numJobs,int numMachines);
    void writeDependentSetupTimes(int numJobs,int numMachines);
    int minJobSetupTime(const int theJob,const int theMachine);

    void computeOperationMachineAssignments(int numJobs,int numMachines);
    int operationIndex(int job,int machine,int numMachines);
    long taillardLowerBound(int numJobs,int numMachines);
    void writeProblem(int numJobs,int numMachines);
    void writeLowerBounds(int numJobs,int numMachines, int* mkp);

  };


  /**********************************************
   * Combination Iterator
   *********************************************/

  ///   An iterator on all combinations of K elements among N.
  /**
     The algo implemented is Knuth's.
     The combination is stored in the array "comb"
     the following combination is obtained 
     by calling "next".
   
     Note that the combination comb among C(N,K)
     is obtained after comb2int(comb, K) calls to next.
     Therefore a combination can be referrenced 
     with a integer.
  */
  class Combination_iterator {

  private:
    int N;
    int K;
    int j;
  
  public:
    /// The combination 
    int *comb;

    /**@name Constructors*/
    //@{
    ///
    Combination_iterator() {}
    ///
    Combination_iterator(int, int);
    ///
    ~Combination_iterator();
    //@}

    /// Allocate the array comb
    void create(int, int);
  
    /// Initialise N and K
    void init(int, int);
    /// Initialise a new K and keep N as it is
    void init(int);
    /**
       Compute the combination following comb and put the result in comb
    */
    bool operator++();
    bool next();
    void print(std::ostream& o);
  };
  /**
     Compute the binomiall coeff
  */
  unsigned int binomialCoeff(int N, int K);
  /**
     Compute the rank of a combination in the iteration
  */
  unsigned int comb2int(int *cm, int K);


  // /**********************************************
  //  * Timing Memory and Command line utilities 
  //  *********************************************/
  
  // double getRunTime();
  // unsigned long int getMemory();
  // void getCommandLine(const char**,int*,int,const char**,const char**,int,char**,int);
  
};

#endif // _GENERATOR_H





