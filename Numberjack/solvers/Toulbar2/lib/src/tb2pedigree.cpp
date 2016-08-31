/*
 * **************** Read pedigree in pre format files **************************
 *
 */

#include "toulbar2lib.hpp"
#include "tb2enumvar.hpp"
#include "tb2pedigree.hpp"

void Pedigree::iniProb( WCSP* wcsp ) {
    TProb TopProb = 0.;

    ToulBar2::NormFactor = (-1/Log1p(-Exp10(-(TProb)ToulBar2::resolution)));
    if (ToulBar2::NormFactor > (Pow( (TProb)2., (TProb)INTEGERBITS)-1)/(TProb)ToulBar2::resolution) {
        cerr << "This resolution cannot be ensured on the data type used to represent costs." << endl;
        exit(EXIT_FAILURE);
    }

    int nballeles = alleles.size()-1;

    TopProb = 0.;

    if(nballeles > 1) {
        int ngenotyped = genotypes.size();
        while(ngenotyped) {
            TopProb += -Log(ToulBar2::errorg / (TProb)(nballeles-1)) * ToulBar2::NormFactor;
            ngenotyped--;
        }
    }

    for (vector<Individual>::iterator iter = pedigree.begin(); iter != pedigree.end(); ++iter) {
        Individual& individual = *iter;
        if(individual.typed) {
            if(individual.mother && individual.father) {
                TopProb += -Log(0.25) * ToulBar2::NormFactor;
            }
            else if(individual.mother || individual.father) {
                TopProb += -Log(0.50) * ToulBar2::NormFactor;
            }
            else
            {
                TProb minp = 1.;
                switch(ToulBar2::foundersprob_class)
                {
                case 0:  TopProb += -Log(1./(nballeles * nballeles)) * ToulBar2::NormFactor;
                break;

                case 1:  for (map<int,int>::iterator iter = freqalleles.begin(); iter != freqalleles.end(); ++iter) {
                    TProb p = (TProb) (iter->second * iter->second)/ (TProb)(genotypes.size() * genotypes.size() * 4);
                    if(p < minp) minp = p;
                }

                TopProb += -Log(minp) * ToulBar2::NormFactor;
                break;

                default: foundersprob.clear();
                assert((int) ToulBar2::allelefreqdistrib.size() == nballeles);
                for (int i=1; i<=nballeles; i++) { /* i = first allele of child  */
                    for (int j=i; j<=nballeles; j++) { /* j = second allele of child */
                        // case i!=j must be counted twice (i,j and j,i)
                        foundersprob.push_back((TProb)((i!=j)?2.:1.) * (TProb)ToulBar2::allelefreqdistrib[i-1] * (TProb)ToulBar2::allelefreqdistrib[j-1]);
                    }
                }
                for (vector<TProb>::iterator iter = foundersprob.begin(); iter != foundersprob.end(); ++iter) {
                    if(*iter < minp) minp = *iter;
                }
                TopProb += -Log(minp) * ToulBar2::NormFactor;
                }
            }
        }
    }
    if(TopProb > to_double(MAX_COST)) {
        cerr << "Overflow: product of min probabilities < size of used datatype." << endl;
        exit(EXIT_FAILURE);
    }
    wcsp->updateUb((Cost) ((Long) TopProb));
}

typedef struct {
    EnumeratedVariable *var;
    vector<Cost> costs;
} TemporaryUnaryConstraint;

Individual::Individual(int ind)
{
    individual = ind;
    varindex = -1;
    father = 0;
    mother = 0;
    sex = MALE;
    genotype.allele1 = 0;
    genotype.allele2 = 0;
    genotype.fixed = false;
    typed = false;
    generation = -1;
    nbtyped = 0;
}

void Individual::print(ostream& os)
{
    os << individual << " " << father << " " << mother << " " << sex;
    if (genotypes.size()>0) {
        for (unsigned int i = 0; i < genotypes.size(); ++i)
            cout << " " << genotypes[i].allele1 << " " << genotypes[i].allele2;
    } else {
        os << " " << genotype.allele1 << " " << genotype.allele2;
    }
    os << endl;
}

void Pedigree::typeAscendants(int individual)
{
    if (individual > 0) {
        assert(individuals.count(individual)!=0);
        int index = individuals[individual];
        if (!pedigree[index].typed) {
            pedigree[index].typed = true;
            nbtyped++;
            typeAscendants(pedigree[index].father);
            typeAscendants(pedigree[index].mother);
        }
    }
}

void Pedigree::read(const char *fileName, WCSP *wcsp)
{
    bayesian = false;
    readPedigree(fileName, wcsp);
    buildWCSP(fileName, wcsp);
}

void Pedigree::read_bayesian(const char *fileName, WCSP *wcsp )
{
    bayesian = true;
    readPedigree(fileName, wcsp);
    buildWCSP_bayesian(fileName, wcsp);
}

inline bool cmp_generation(Individual i1, Individual i2) { return i1.generation < i2.generation || (i1.generation == i2.generation && i1.individual < i2.individual); }

int Pedigree::fixGenerationNumber(int index)
{
    if (pedigree[index].generation!=-1) return pedigree[index].generation;
    else {
        if (pedigree[index].father>0 && pedigree[index].mother>0) {
            pedigree[index].generation = max(fixGenerationNumber(individuals[pedigree[index].father]), fixGenerationNumber(individuals[pedigree[index].mother]))+1;
        } else if (pedigree[index].father>0 && pedigree[index].mother==0) {
            pedigree[index].generation = fixGenerationNumber(individuals[pedigree[index].father])+1;
        } else if (pedigree[index].mother>0 && pedigree[index].father==0) {
            pedigree[index].generation = fixGenerationNumber(individuals[pedigree[index].mother])+1;
        } else {
            assert(pedigree[index].mother==0 && pedigree[index].father==0);
            pedigree[index].generation = 1;
        }
        return pedigree[index].generation;
    }
}

// warning! locus information is not used: assume that only one locus is defined in the pedigree file
void Pedigree::readPedigree(const char *fileName, WCSP *wcsp)
{
    int individual;
    int nbindividuals = 0;
    int nballeles = 0;
    int nbtypings = 0;
    map<int, int> allelesInv;
    int maxallele = 0;

    string strfile(fileName);
    int pos = strfile.find_last_of(".");
    string errorfilename = strfile.substr(0,pos) + ".errors";

    // open the file
    ifstream file(fileName);
    if (!file) {
        cerr << "Could not open file " << fileName << endl;
        exit(EXIT_FAILURE);
    }

    ifstream fileErrors( errorfilename.c_str() );
    if (fileErrors) ToulBar2::consecutiveAllele = true;


    while (file) {
        int cur_locus = -1;
        individual = 0;

        file >> cur_locus;
        if (!file) break;
        if (locus == -1) locus = cur_locus;
        if (locus != cur_locus) {
            cerr << "Pedigree datafile contains more than one locus!" << endl;
            exit(EXIT_FAILURE);
        }

        file >> individual;
        if (!file) {
            cerr << "Wrong data after individual " << individual << endl;
            exit(EXIT_FAILURE);
        }
        assert(individual != 0);
        if (individuals.count(individual)==0) {
            individuals[individual] = nbindividuals;
            Individual geno(individual);
            pedigree.push_back(geno);
            nbindividuals++;
        }

        file >> pedigree[individuals[individual]].father;
        if (!file) {
            cerr << "Wrong data after individual " << individual << endl;
            exit(EXIT_FAILURE);
        }
        if (pedigree[individuals[individual]].father > 0 && individuals.count(pedigree[individuals[individual]].father)==0) {
            individuals[pedigree[individuals[individual]].father] = nbindividuals;
            Individual geno(pedigree[individuals[individual]].father);
            geno.sex = MALE;
            pedigree.push_back(geno);
            nbindividuals++;
        }

        file >> pedigree[individuals[individual]].mother;
        if (!file) {
            cerr << "Wrong data after individual " << individual << endl;
            exit(EXIT_FAILURE);
        }
        if (pedigree[individuals[individual]].mother > 0 && individuals.count(pedigree[individuals[individual]].mother)==0) {
            individuals[pedigree[individuals[individual]].mother] = nbindividuals;
            Individual geno(pedigree[individuals[individual]].mother);
            geno.sex = FEMALE;
            pedigree.push_back(geno);
            nbindividuals++;
        }

        if (pedigree[individuals[individual]].father==0 && pedigree[individuals[individual]].mother==0) {
            pedigree[individuals[individual]].generation = 1;
        } else if (pedigree[individuals[individual]].father>0 && pedigree[individuals[individual]].mother>0 &&
                pedigree[individuals[pedigree[individuals[individual]].father]].generation!=-1 &&
                pedigree[individuals[pedigree[individuals[individual]].mother]].generation!=-1) {
            pedigree[individuals[individual]].generation = max(pedigree[individuals[pedigree[individuals[individual]].father]].generation, pedigree[individuals[pedigree[individuals[individual]].mother]].generation) + 1;
        } else if (pedigree[individuals[individual]].father>0 && pedigree[individuals[individual]].mother==0 &&
                pedigree[individuals[pedigree[individuals[individual]].father]].generation!=-1) {
            pedigree[individuals[individual]].generation = pedigree[individuals[pedigree[individuals[individual]].father]].generation + 1;
        } else if (pedigree[individuals[individual]].father==0 && pedigree[individuals[individual]].mother>0 &&
                pedigree[individuals[pedigree[individuals[individual]].mother]].generation!=-1) {
            pedigree[individuals[individual]].generation = pedigree[individuals[pedigree[individuals[individual]].mother]].generation + 1;
        }

        file >> pedigree[individuals[individual]].sex;
        if (!file) {
            cerr << "Wrong data after individual " << individual << endl;
            exit(EXIT_FAILURE);
        }

        int allele = 0;
        file >> allele;
        if (!file) { cerr << "Wrong data after individual " << individual << endl; exit(EXIT_FAILURE); }
        if (allele < 0) {
            pedigree[individuals[individual]].genotype.fixed = true;
            allele = -allele;
        }
        pedigree[individuals[individual]].genotype.allele1 = allele;

        allele = 0;
        file >> allele;
        if (!file) { cerr << "Wrong data after individual " << individual << endl; exit(EXIT_FAILURE); }
        if (allele < 0) {
            pedigree[individuals[individual]].genotype.fixed = true;
            allele = -allele;
        }
        pedigree[individuals[individual]].genotype.allele2 = allele;

        int allele1 = pedigree[individuals[individual]].genotype.allele1;
        int allele2 = pedigree[individuals[individual]].genotype.allele2;

        if (alleles.count(allele1)==0) {
            nballeles++;
            alleles[allele1] = nballeles;
            freqalleles[ allele1 ] = 1;
            if(allele1 > maxallele) maxallele =  allele1;
        }
        else { freqalleles[ allele1 ]++; }


        if (alleles.count(allele2)==0) {
            nballeles++;
            alleles[allele2] = nballeles;
            freqalleles[ allele2 ] = 1;
            if(allele2 > maxallele) maxallele =  allele2;
        }
        else { freqalleles[ allele2 ]++; }


        if (allele1>0 || allele2>0) {
            nbtypings++;
            genotypes.push_back(individual);
        }
    }

    if(ToulBar2::consecutiveAllele) {
        cout << "Make alleles consecutive and supose 4 allele." << endl;
        maxallele = 4;
        for(int i=1;i<=maxallele;i++) {
            map<int,int>::iterator it = alleles.find(i);
            if(it == alleles.end()) {
                nballeles++;
                alleles[i] = nballeles;
                freqalleles[i] = 0;
            }
        }
    }

    /* re-encoding of alleles */
    int nb = 0;
    if (ToulBar2::verbose >= 2) cout << "Alleles encoding:" << endl;
    for (map<int,int>::iterator iter = alleles.begin(); iter != alleles.end(); ++iter) {
        if ((*iter).first == 0) continue;
        nb++;
        if (ToulBar2::verbose >= 2) cout << (*iter).first << ": " << nb << endl;
        iter->second = nb;
        allelesInv[nb] = (*iter).first;
    }

    assert(nballeles == nb);


    if (ToulBar2::verbose >= 1) cout << "Genotype encoding:" << endl;
    for (int i=1; i<=nballeles; i++) { /* i = first allele of child */
        for (int j=i; j<=nballeles; j++) { /* j = second allele of child */
            Genotype geno;
            geno.allele1 = allelesInv[i];
            geno.allele2 = allelesInv[j];
            genoconvert.push_back(geno);
            if (ToulBar2::verbose >= 1) {
                cout << genoconvert.size()-1 << ": ";
                printGenotype(cout, genoconvert.size()-1);
                cout << endl;
            }
        }
    }

    /* individuals re-ordering by generation levels */
    for (unsigned int i=0; i<pedigree.size(); i++) {
        if (pedigree[i].generation==-1) fixGenerationNumber(i);
        if (pedigree[i].generation > generations) generations = pedigree[i].generation;
    }
    if (ToulBar2::generation) {
        stable_sort(pedigree.begin(), pedigree.end(), cmp_generation);
        for (unsigned int i=0; i<pedigree.size(); i++) {
            individuals[pedigree[i].individual] = i;
        }
    }

    for (unsigned int i=0; i<genotypes.size(); i++) {
        typeAscendants(genotypes[i]);
        if (genotypes[i]>=0 && pedigree[individuals[genotypes[i]]].father > 0 &&
                pedigree[individuals[pedigree[individuals[genotypes[i]]].father]].genotype.allele1>0 &&
                pedigree[individuals[pedigree[individuals[genotypes[i]]].father]].genotype.allele2>0)
            pedigree[individuals[pedigree[individuals[genotypes[i]]].father]].nbtyped++;
        if (genotypes[i]>=0 && pedigree[individuals[genotypes[i]]].mother > 0 &&
                pedigree[individuals[pedigree[individuals[genotypes[i]]].mother]].genotype.allele1>0 &&
                pedigree[individuals[pedigree[individuals[genotypes[i]]].mother]].genotype.allele2>0)
            pedigree[individuals[pedigree[individuals[genotypes[i]]].mother]].nbtyped++;
    }
    cout << nbtyped << " informative individuals found (either genotyped or having a genotyped descendant)." << endl;


    gencorrects.clear();
    if (fileErrors) {
        while(fileErrors) {
            int bonallele1;
            fileErrors >> individual;
            fileErrors >> bonallele1;
            int bonallele2 = pedigree[individuals[individual]].genotype.allele2;
            if(bonallele1 == pedigree[individuals[individual]].genotype.allele1) cout << "error file with a mistake" << endl;
            gencorrects[individual] = convertgen(bonallele1, bonallele2);
        }
    }
    assert(wcsp->numberOfVariables() == 0);
    assert(wcsp->numberOfConstraints() == 0);
}

int Pedigree::convertgen( int allele1, int allele2 ) {
    int nballeles = alleles.size() - 1;
    int bongen = 0;

    if(allele1 > allele2) {
        int alleleaux = allele1;
        allele1 = allele2;
        allele2 = alleleaux;
    }

    for (int i=1; i<=nballeles; i++) {
        for (int j=i; j<=nballeles; j++) {
            if((i==allele1) && (j==allele2)) return bongen;
            bongen++;
        }}
    return -1;
}



void Pedigree::buildWCSP(const char *fileName, WCSP *wcsp)
{
    ifstream file(fileName);

    vector<TemporaryUnaryConstraint> unaryconstrs;

    int nbindividuals = individuals.size();
    int nballeles = alleles.size() - 1;
    int nbfounders = 0;
    int nbtypings = genotypes.size();

    wcsp->updateUb(nbtypings+1);

    /* create variables */
    int nbvar = 0;
    for (int i=0; i<nbindividuals; i++) {
        if (pedigree[i].father == 0 && pedigree[i].mother == 0) nbfounders++;
        if (pedigree[i].typed) {
            string varname;
            varname = to_string(pedigree[i].individual);
            wcsp->makeEnumeratedVariable(varname, 0, nballeles*(nballeles+1)/2 - 1);
            pedigree[i].varindex = nbvar;
            nbvar++;
        }
    }

    /* create ternary Mendelian hard constraint table */
    vector<Cost> costs3;
    for (int k=1; k<=nballeles; k++) { /* k = first allele of father */
        for (int l=k; l<=nballeles; l++) { /* l = second allele of father */
            for (int m=1; m<=nballeles; m++) { /* m = first allele of mother */
                for (int n=m; n<=nballeles; n++) { /* n = second allele of mother */
                    for (int i=1; i<=nballeles; i++) { /* i = first allele of child */
                        for (int j=i; j<=nballeles; j++) { /* j = second allele of child */
                            costs3.push_back(((i==k && j==m) || (i==k && j==n) || (i==l && j==m) || (i==l && j==n) || (i==m && j==k) || (i==m && j==l) || (i==n && j==k) || (i==n && j==l))?0:wcsp->getUb()*MEDIUM_COST);
                        }
                    }
                }
            }
        }
    }

    /* create binary Mendelian hard constraint table */
    vector<Cost> costs2;
    for (int k=1; k<=nballeles; k++) { /* k = first allele of father or mother */
        for (int l=k; l<=nballeles; l++) { /* l = second allele of father or mother */
            for (int i=1; i<=nballeles; i++) { /* i = first allele of child */
                for (int j=i; j<=nballeles; j++) { /* j = second allele of child */
                    costs2.push_back((i==k || i==l || j==k || j==l)?0:wcsp->getUb()*MEDIUM_COST);
                }
            }
        }
    }

    /* create constraint network */
    while (file) {
        int cur_locus = -1;
        int individual = 0;
        int father = 0;
        int mother = 0;
        int sex = -1;
        int allele1 = 0;
        int allele2 = 0;

        file >> cur_locus;
        if (!file) break;
        file >> individual;
        file >> father;
        file >> mother;
        file >> sex;
        file >> allele1;
        if (allele1 < 0) allele1 = -allele1;
        allele1 = alleles[allele1];
        file >> allele2;
        if (allele2 < 0) allele2 = -allele2;
        allele2 = alleles[allele2];
        if (!pedigree[individuals[individual]].typed) continue;

        /* add unary costs (soft constraint) if genotyping is given */
        if (allele1 > 0 || allele2 > 0) {
            EnumeratedVariable *var = (EnumeratedVariable *) wcsp->getVar(pedigree[individuals[individual]].varindex);
            TemporaryUnaryConstraint unaryconstr;
            unaryconstr.var = var;
            for (int i=1; i<=nballeles; i++) { /* i = first allele of child */
                for (int j=i; j<=nballeles; j++) { /* j = second allele of child */
                    if ((allele1>0 && allele2>0 && ((i==allele1 && j==allele2) || (i==allele2 && j==allele1)))
                            || ((allele1==0 || allele2==0) && (i==allele1 || i==allele2 || j==allele1 || j==allele2))) {
                        unaryconstr.costs.push_back(0);
                    } else {
                        unaryconstr.costs.push_back((pedigree[individuals[individual]].genotype.fixed)?wcsp->getUb():1);
                    }
                }
            }
            unaryconstrs.push_back(unaryconstr);
            var->queueNC();
        }

        /* add ternary or binary Mendelian hard constraint */
        if (father > 0 || mother > 0) {
            if (father > 0 && mother > 0) {
                assert(pedigree[individuals[pedigree[individuals[individual]].father]].typed);
                assert(pedigree[individuals[pedigree[individuals[individual]].mother]].typed);
                wcsp->postTernaryConstraint(pedigree[individuals[pedigree[individuals[individual]].father]].varindex,pedigree[individuals[pedigree[individuals[individual]].mother]].varindex,pedigree[individuals[individual]].varindex,costs3);
            } else if (father > 0) {
                wcsp->postBinaryConstraint(pedigree[individuals[pedigree[individuals[individual]].father]].varindex,pedigree[individuals[individual]].varindex,costs2);
            } else {
                wcsp->postBinaryConstraint(pedigree[individuals[pedigree[individuals[individual]].mother]].varindex,pedigree[individuals[individual]].varindex,costs2);
            }
        }
    }
    wcsp->sortConstraints();

    // apply basic initial propagation AFTER complete network loading
    for (unsigned int u=0; u<unaryconstrs.size(); u++) {
        for (unsigned int a = 0; a < unaryconstrs[u].var->getDomainInitSize(); a++) {
            if (unaryconstrs[u].costs[a] > 0) unaryconstrs[u].var->project(unaryconstrs[u].var->toValue(a), unaryconstrs[u].costs[a]);
        }
        unaryconstrs[u].var->findSupport();
    }

    if (ToulBar2::verbose >= 0) {
        cout << "Read pedigree with " << nbindividuals << " individuals, " << nbfounders << " founders, " << nballeles << " alleles, " << nbtypings << " genotypings and " << generations << " generations." << endl;
    }
}



void Pedigree::buildWCSP_bayesian( const char *fileName, WCSP *wcsp )
{
    ifstream file(fileName);

    vector<TemporaryUnaryConstraint> unaryconstrs;
    int nbindividuals = individuals.size();
    int nballeles = alleles.size() - 1;
    int domsize = nballeles * ( nballeles + 1 ) / 2;
    int nbfounders = 0;
    int nbtypings = genotypes.size();

    iniProb( wcsp );

    /* create variables */
    int nbvar = 0;
    for (int i=0; i<nbindividuals; i++) {
        if (pedigree[i].father == 0 && pedigree[i].mother == 0) nbfounders++;
        if (pedigree[i].typed) {
            string varname;
            varname = to_string(pedigree[i].individual);
            wcsp->makeEnumeratedVariable(varname, 0, nballeles*(nballeles+1)/2 - 1);
            pedigree[i].varindex = nbvar;
            nbvar++;
        }
    }

    /* create ternary Mendelian hard constraint table */
    vector<Cost> costs3;
    for (int k=1; k<=nballeles; k++) { /* k = first allele of father */
        for (int l=k; l<=nballeles; l++) { /* l = second allele of father */
            for (int m=1; m<=nballeles; m++) { /* m = first allele of mother */
                for (int n=m; n<=nballeles; n++) { /* n = second allele of mother */
                    for (int i=1; i<=nballeles; i++) { /* i = first allele of child */
                        for (int j=i; j<=nballeles; j++) { /* j = second allele of child */
                            TProb p = 0.;
                            if((i==k && j==m) || (i==m && j==k)) p += 0.25;
                            if((i==k && j==n) || (i==n && j==k)) p += 0.25;
                            if((i==l && j==m) || (i==m && j==l)) p += 0.25;
                            if((i==l && j==n) || (i==n && j==l)) p += 0.25;
                            costs3.push_back( wcsp->Prob2Cost(p) );
                        }
                    }
                }
            }
        }
    }

    map<int, int> allelesInv;

    int i,j;
    switch(ToulBar2::foundersprob_class) {
    case 0: foundersprob.clear();
    for (i=1; i<=nballeles; i++) { /* i = first allele of child  */
        for (j=i; j<=nballeles; j++) { /* j = second allele of child */
            foundersprob.push_back( ((i!=j)?2.:1.)/(TProb)(nballeles * nballeles) );
        }
    }
    break;

    case 1: foundersprob.clear();
    for (map<int,int>::iterator iter = alleles.begin(); iter != alleles.end(); ++iter) {
        allelesInv[iter->second] = iter->first;
    }
    for (i=1; i<=nballeles; i++) { /* i = first allele of child  */
        for (j=i; j<=nballeles; j++) { /* j = second allele of child */
            // case i!=j must be counted twice (i,j and j,i)
            foundersprob.push_back((TProb)((i!=j)?2.:1.) * (TProb)freqalleles[ allelesInv[i] ] * (TProb)freqalleles[ allelesInv[j] ] / (TProb)(nbtypings * nbtypings * 4) );
        }
    }
    break;

    default: foundersprob.clear();
    assert((int) ToulBar2::allelefreqdistrib.size() == nballeles);
    for (i=1; i<=nballeles; i++) { /* i = first allele of child  */
        for (j=i; j<=nballeles; j++) { /* j = second allele of child */
            // case i!=j must be counted twice (i,j and j,i)
            foundersprob.push_back((TProb)((i!=j)?2.:1.) * (TProb)ToulBar2::allelefreqdistrib[i-1] * (TProb)ToulBar2::allelefreqdistrib[j-1]);
        }
    }
    break;
    }
    if (ToulBar2::verbose >= 1) {
        cout << "Genotype prior:" << endl;
        for (unsigned int n = 0; n < genoconvert.size(); n++) {
            printGenotype(cout, n);
            cout << " " << foundersprob[n] << endl;
        }
    }

    /* create binary Mendelian hard constraint table */
    vector<Cost> costs2;
    for (int k=1; k<=nballeles; k++) { /* k = first allele of father or mother */
        for (int l=k; l<=nballeles; l++) { /* l = second allele of father or mother */
            for (int i=1; i<=nballeles; i++) { /* i = first allele of child */
                for (int j=i; j<=nballeles; j++) { /* j = second allele of child */
                    TProb p = 0;
                    if(i==k || i==l || j==k || j==l) {
                        if(k == l)   p += 1. / (TProb)nballeles;
                        else 		p += 1. / ((TProb)nballeles + (TProb)nballeles - 1);
                    }
                    costs2.push_back( wcsp->Prob2Cost(p) );
                }
            }
        }
    }


    /* create constraint network */
    while (file) {
        int cur_locus = -1;
        int individual = 0;
        int father = 0;
        int mother = 0;
        int sex = -1;
        int allele1 = 0;
        int allele2 = 0;

        file >> cur_locus;
        if (!file) break;
        file >> individual;
        file >> father;
        file >> mother;
        file >> sex;
        file >> allele1;
        if (allele1 < 0) allele1 = -allele1;
        allele1 = alleles[allele1];
        file >> allele2;
        if (allele2 < 0) allele2 = -allele2;
        allele2 = alleles[allele2];
        if (!pedigree[individuals[individual]].typed) continue;

        EnumeratedVariable *var = (EnumeratedVariable *) wcsp->getVar(pedigree[individuals[individual]].varindex);

        /* add unary costs (soft constraint) if genotyping is given */
        if (allele1 > 0 || allele2 > 0) {
            TemporaryUnaryConstraint unaryconstr;
            unaryconstr.var = var;
            for (int i=1; i<=nballeles; i++) { /* i = first allele of child  */
                for (int j=i; j<=nballeles; j++) { /* j = second allele of child */
                    bool typed = allele1>0 && allele2>0;
                    bool halftyped = allele1>0 || allele2>0;

                    bool theone = (i==allele1 && j==allele2) || (i==allele2 && j==allele1);
                    bool posible = (i==allele1 || i==allele2 || j==allele1 || j==allele2);

                    bool fixed = pedigree[individuals[individual]].genotype.fixed;

                    TProb p = 0;
                    Cost penalty = 0;
                    if (typed) {
                        if(theone) p = 1. - ToulBar2::errorg;
                        else	   {p = ToulBar2::errorg / (TProb)(domsize-1);
                        penalty = pedigree[individuals[individual]].nbtyped; }
                    }
                    else if(halftyped) {
                        if(posible) p = (1. - ToulBar2::errorg) / (TProb)nballeles;
                        else {p = ToulBar2::errorg / (TProb)(domsize - nballeles);
                        penalty = pedigree[individuals[individual]].nbtyped; }
                    }
                    if (ToulBar2::pedigreePenalty>0 && ToulBar2::verbose >= 1) cout << individual << ": "  << penalty << " nbtyped " << ((penalty>ToulBar2::pedigreePenalty)?wcsp->Cost2LogProb(-((penalty>0)?wcsp->Prob2Cost(to_double(penalty)):MIN_COST))/Log(10.):0.) << " log10like " << -((penalty>ToulBar2::pedigreePenalty)?wcsp->Prob2Cost(to_double(penalty)):MIN_COST) << " cost" << endl;
                    unaryconstr.costs.push_back((typed && fixed && !theone)?wcsp->getUb():(wcsp->Prob2Cost(p) - ((ToulBar2::pedigreePenalty>0 && penalty>ToulBar2::pedigreePenalty)?wcsp->Prob2Cost(to_double(penalty)):MIN_COST)) );
                }
            }
            unaryconstrs.push_back(unaryconstr);
            var->queueNC();
        }

        int id_father = individuals[pedigree[individuals[individual]].father];
        int id_mother = individuals[pedigree[individuals[individual]].mother];

        /* add ternary or binary Mendelian hard constraint */
        if (father > 0 || mother > 0) {
            if (father > 0 && mother > 0) {
                assert(pedigree[id_father].typed);
                assert(pedigree[id_mother].typed);
                wcsp->postTernaryConstraint(pedigree[id_father].varindex,pedigree[id_mother].varindex,pedigree[individuals[individual]].varindex,costs3);
            } else if (father > 0) {
                wcsp->postBinaryConstraint(pedigree[id_father].varindex,pedigree[individuals[individual]].varindex,costs2);
            } else {
                wcsp->postBinaryConstraint(pedigree[id_mother].varindex,pedigree[individuals[individual]].varindex,costs2);
            }
        }
        else {
            TemporaryUnaryConstraint unaryconstr;
            unaryconstr.var = var;
            for (vector<TProb>::iterator iter = foundersprob.begin(); iter != foundersprob.end(); ++iter) {
                unaryconstr.costs.push_back( wcsp->Prob2Cost( *iter) );
            }
            unaryconstrs.push_back(unaryconstr);
            var->queueNC();
        }
    }
    wcsp->sortConstraints();

    // apply basic initial propagation AFTER complete network loading
    for (unsigned int u=0; u<unaryconstrs.size(); u++) {
        for (unsigned int a = 0; a < unaryconstrs[u].var->getDomainInitSize(); a++) {
            if (unaryconstrs[u].costs[a] > 0) unaryconstrs[u].var->project(unaryconstrs[u].var->toValue(a), unaryconstrs[u].costs[a]);
        }
        unaryconstrs[u].var->findSupport();
    }

    if (ToulBar2::verbose >= 0) {
        int nbtypings = genotypes.size();
        cout << "Read pedigree with " << nbindividuals << " individuals, " << nbfounders << " founders, " << nballeles << " alleles, " << nbtypings << " genotypings and " << generations << " generations." << endl;
        cout << "Bayesian MPE (genotyping error rate: " <<  ToulBar2::errorg << ", genotype prior: " << ToulBar2::foundersprob_class << ", precision(1-10^-p): " << ToulBar2::resolution << ", normalization: " << ToulBar2::NormFactor << ", ub: " << wcsp->getUb() << ")" << endl;

    }
}


void Pedigree::printCorrectSol(WCSP *wcsp)
{
    if(!gencorrects.size()) return;

    string problemname = ToulBar2::problemsaved_filename;
    if (problemname.rfind( ".wcsp" ) != string::npos) problemname.replace( problemname.rfind( ".wcsp" ), 5, "_correct.sol" );
    if (problemname.rfind( ".pre" ) != string::npos) problemname.replace( problemname.rfind( ".pre" ), 4, "_correct.sol" );
    if (problemname.rfind( "_correct.sol" ) == string::npos) problemname = problemname + to_string("_correct.sol");
    ofstream file(problemname.c_str());
    if (!file) {
        cerr << "Could not write file " << "solution" << endl;
        exit(EXIT_FAILURE);
    }

    for(vector<Individual>::iterator it = pedigree.begin(); it != pedigree.end(); ++it )
    {
        Individual& ind = *it;
        int allele1 = ind.genotype.allele1;
        int allele2 = ind.genotype.allele2;
        if((allele1 > 0) || (allele2 > 0)) {
            map<int,int>::iterator it = gencorrects.find(ind.individual);
            if(it != gencorrects.end()) file << " " << it->second;
            else file << " " << convertgen(allele1, allele2);
        }
        else file << " " << -1;
    }
    file << endl;
}


void Pedigree::printSol(WCSP *wcsp)
{
    string problemname = ToulBar2::problemsaved_filename;
    if (problemname.rfind( ".wcsp" ) != string::npos) problemname.replace( problemname.rfind( ".wcsp" ), 5, ".sol" );
    if (problemname.rfind( ".pre" ) != string::npos) problemname.replace( problemname.rfind( ".pre" ), 4, ".sol" );
    if (problemname.rfind( ".sol" ) == string::npos) problemname = problemname + to_string(".sol");
    ofstream file(problemname.c_str());
    if (!file) {
        cerr << "Could not write file " << "solution" << endl;
        exit(EXIT_FAILURE);
    }

    for(vector<Individual>::iterator it = pedigree.begin(); it != pedigree.end(); ++it )
    {
        Individual& ind = *it;
        if(ind.typed) file << " " << wcsp->getValue(ind.varindex);
        else file << " " << -1;
    }
    file << endl;
}

void Pedigree::printCorrection(WCSP *wcsp)
{
    bool errorinfo = gencorrects.size() > 0;

    int ncorrect = 0;
    int ncorrections = 0;
    int ncorrectok = 0;

    int nbcorrection = 0;

    TProb penalty = 0.;

    cout << "Correction:";
    for (unsigned int i=0; i<genotypes.size(); i++) {
        int sol = wcsp->getValue(pedigree[individuals[genotypes[i]]].varindex);
        int a1 = genoconvert[sol].allele1;
        int a2 = genoconvert[sol].allele2;
        int allele1 = pedigree[individuals[genotypes[i]]].genotype.allele1;
        int allele2 = pedigree[individuals[genotypes[i]]].genotype.allele2;
        if (!((allele1>0 && allele2>0 && ((a1==allele1 && a2==allele2) || (a1==allele2 && a2==allele1)))
                || ((allele1==0 || allele2==0) && (a1==allele1 || a1==allele2 || a2==allele1 || a2==allele2)))) {
            cout << " " << genotypes[i];
            nbcorrection++;

            if (ToulBar2::pedigreePenalty>0 && pedigree[individuals[genotypes[i]]].nbtyped>ToulBar2::pedigreePenalty) {
                cout << "*";
                penalty -= Log(pedigree[individuals[genotypes[i]]].nbtyped);
            }

            if(errorinfo) {
                map<int,int>::iterator it = gencorrects.find(genotypes[i]);
                if(it != gencorrects.end()) {
                    ncorrect++;
                    int allellereal = it->second;
                    if((allellereal == a1) && (a1 != allele1)) ncorrectok++;
                }
            }
            ncorrections++;
        }
    }
    cout << " (" << nbcorrection;
    if (ToulBar2::pedigreePenalty>0) {
        cout << " " << wcsp->Cost2LogProb(wcsp->getLb()) - penalty << " " << penalty;
    }
    cout << ")";
    cout << endl;
    if(errorinfo) {
        cout << "Info errorfile: ";
        cout << ncorrect << "/"  << ncorrections << ", " <<  ncorrectok << "/" << ncorrect << ", " << gencorrects.size() << " original" << endl;
        cout << endl;
    }
}


void Pedigree::printGenotype(ostream& os, Value value)
{
    os << genoconvert[value].allele1 << "/" << genoconvert[value].allele2;
}

void Pedigree::save(const char *fileName, WCSP *wcsp, bool corrected, bool reduced)
{
    assert(!(corrected && reduced));

    // open the file
    ofstream file(fileName);
    if (!file) {
        cerr << "Could not open file " << fileName << endl;
        exit(EXIT_FAILURE);
    }

    for (map<int,int>::iterator iter = individuals.begin(); iter != individuals.end(); ++iter) {
        if ((*iter).first == 0) continue;
        if (reduced && (pedigree[(*iter).second].varindex < 0 || pedigree[(*iter).second].varindex >= (int) wcsp->numberOfVariables() || wcsp->assigned(pedigree[(*iter).second].varindex))) {
            continue;
        }
        file << locus << " ";
        if (corrected && pedigree[(*iter).second].varindex >= 0 && pedigree[(*iter).second].varindex < (int) wcsp->numberOfVariables() && wcsp->assigned(pedigree[(*iter).second].varindex)) {
            int sol = wcsp->getValue(pedigree[(*iter).second].varindex);
            int a1 = genoconvert[sol].allele1;
            int a2 = genoconvert[sol].allele2;
            int allele1 = pedigree[(*iter).second].genotype.allele1;
            int allele2 = pedigree[(*iter).second].genotype.allele2;
            if ((allele1>0 && allele2>0 && !((a1==allele1 && a2==allele2) || (a1==allele2 && a2==allele1)))
                    || ((allele1==0 || allele2==0) && !(allele1==0 && allele2==0) && !(a1==allele1 || a1==allele2 || a2==allele1 || a2==allele2))
                    || ToulBar2::pedigreeCorrectionMode == 2) {
                if (ToulBar2::pedigreeCorrectionMode > 0) {
                    pedigree[(*iter).second].genotype.allele1 = a1;
                    pedigree[(*iter).second].genotype.allele2 = a2;
                } else {
                    pedigree[(*iter).second].genotype.allele1 = 0;
                    pedigree[(*iter).second].genotype.allele2 = 0;
                }
            }
            pedigree[(*iter).second].print(file);
            pedigree[(*iter).second].genotype.allele1 = allele1;
            pedigree[(*iter).second].genotype.allele2 = allele2;
        } else {
            pedigree[(*iter).second].print(file);
        }
    }
}

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

