/** \file tb2pedigree.hpp
 *  \brief Pedigree data structure
 * 
 */

#ifndef TB2PEDIGREE_HPP_
#define TB2PEDIGREE_HPP_

#include "tb2wcsp.hpp"

class Genotype {
public:
    int allele1;
    int allele2;
    bool fixed;
};

typedef enum {MALE=1, FEMALE=2} Sex;

class Individual {
public:
    int individual;
    int varindex;
    int father;
    int mother;
    int sex;
    Genotype genotype;
    vector<Genotype> genotypes;
    bool typed; // true if one of its descendant children (or itself) is typed
    int generation;
    int nbtyped; // number of genotyped children

    Individual(int ind);

    void print(ostream& os);
};

class Pedigree {
    int locus;    				  // same locus for all the genotypes
    vector<Individual> pedigree;    // list of individuals
    vector<int> genotypes;          // list of genotyped individuals id.
    vector<Genotype> genoconvert;   // convert domain value to genotype
    map<int, int> individuals;      // sorted list of pair<individual id, pedigree id>
    map<int, int> alleles;          // sorted list of pair<allele number, encoding consecutive number>
    int nbtyped;  				  // number of individuals with a genotyped descendant
    int generations;
    bool bayesian;
    vector<TProb> foundersprob;
    map<int, int> freqalleles;      // frequencies of original alleles: freqalleles[allele number] = frequency
    map<int, int> gencorrects;

    void typeAscendants(int individual);
    int fixGenerationNumber(int index);

public:
    Pedigree() : locus(-1), nbtyped(0), generations(0), bayesian(false) {alleles[0] = 0;}

    void read(const char *fileName, WCSP *wcsp);
    void read_bayesian(const char *fileName, WCSP *wcsp);
    void save(const char *fileName, WCSP *wcsp, bool corrected, bool reduced);

    void readPedigree(const char *fileName, WCSP *wcsp);
    void buildWCSP(const char *fileName, WCSP *wcsp);
    void buildWCSP_bayesian(const char *fileName, WCSP *wcsp );
    void iniProb( WCSP* wcsp );

    int convertgen( int allele1, int allele2 );

    void printSol(WCSP *wcsp);
    void printCorrectSol(WCSP *wcsp);
    void printCorrection(WCSP *wcsp);

    void printGenotype(ostream& os, Value value);
};

#endif /*TB2PEDIGREE_HPP_*/

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

