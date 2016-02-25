/** \file tb2haplotype.hpp
 *  \brief Haplotype data structure
 *
 */

#ifndef TB2HAPLOTYPE_HPP_
#define TB2HAPLOTYPE_HPP_

#include "tb2wcsp.hpp"
#include "tb2pedigree.hpp"

struct classcomp {
    bool operator() (const pair<int,int>& lp, const pair<int,int>& rp) const
    {
        if(lp.first>rp.first) return false;
        else
            if(lp.first<rp.first) return true;
            else
                if(lp.second>rp.second) return false;
                else
                    if(lp.second<rp.second) return true;
                    else return false;
    }
};

class Haplotype {
    int family;    				  // same locus for all the genotypes
    int nbloci;
    vector<Individual> pedigree;    // list of individuals
    vector<int> genotypes;          // list of genotyped individuals id.
    map<int,vector<Genotype> > genoconvert;   // convert domain value to genotype
    map<int, int> individuals;      // sorted list of pair<individual id, pedigree id>
    map<int, map<int,int> > alleles;          // sorted list of pair<locus number, allele number, encoding consecutive number>
    int nbtyped;  				  // number of individuals with a genotyped descendant
    int generations;
    bool bayesian;
    vector<TProb> foundersprob;
    map< int, map<int,int> > freqalleles;      // frequencies of original alleles: freqalleles[allele number] = frequency
    map<int, int> gencorrects;

    vector<Double> maplocus; // marker map
    map<int, vector<int> > transmission; //pair <individual id, transmission vector>
    map< pair<int,int>,Double,classcomp > W;
    int sire;
    Double K;
    //Double multiplier; // for conversion in integer cost


    void typeAscendants(int individual);
    int fixGenerationNumber(int index);

public:
    Haplotype() : family(-1), nbtyped(0), generations(0), bayesian(false),sire(-1), K(0.0) {/*alleles[0] = 0;*/}

    void read(const char *fileName, WCSP *wcsp);
    void read_bayesian(const char *fileName, WCSP *wcsp);
    void save(const char *fileName, WCSP *wcsp, bool corrected, bool reduced);

    void readPedigree(const char *fileName, WCSP *wcsp);
    void readMap(const char *fileName);
    void buildWCSP(const char *fileName, WCSP *wcsp);
    void buildWCSP_haplotype(const char *fileName, WCSP *wcsp);
    void buildWCSP_bayesian(const char *fileName, WCSP *wcsp );
    void iniProb( WCSP* wcsp );

    int convertgen( int locus, int allele1, int allele2 );

    void printSol(WCSP *wcsp);
    void printCorrectSol(WCSP *wcsp);
    void printCorrection(WCSP *wcsp);

    void printGenotype(ostream& os, Value value, int locus);
    void printHaplotype(ostream& os, Value value, int locus);

    void initTransmission();
    void sparse_matrix();
    Double haldane(Double x){return 0.5*(1-exp(-2.0*abs(x)));}
    Double getK(){return K;}
    Double Cost2LogProb(Cost c)const{return K - 4*to_double(c)/ ToulBar2::NormFactor;}
};

#endif /*TB2HAPLOTYPE_HPP_*/

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

