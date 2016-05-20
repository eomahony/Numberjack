/*
 * **************** Read pedigree in pre format files **************************
 *
 */

#include "toulbar2lib.hpp"
#include "tb2enumvar.hpp"
#include "tb2haplotype.hpp"

// A MODIFIER POUR PLUSIEURS LOCUS (freqalleles.find(locus))
void Haplotype::iniProb( WCSP* wcsp ) {
    TProb TopProb = 0.;

    int locus = 0;

    ToulBar2::NormFactor = (-1.0/Log1p(-Exp10(-(TProb)ToulBar2::resolution)));
    if (ToulBar2::NormFactor > (Pow( (TProb)2., (TProb)INTEGERBITS)-1)/((TProb)ToulBar2::resolution)) {
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

                case 1:  for (map<int,int>::iterator iter = freqalleles.find(locus)->second.begin(); iter != freqalleles.find(locus)->second.end(); ++iter) {
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

void Haplotype::typeAscendants(int individual)
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

void Haplotype::read(const char *fileName, WCSP *wcsp)
{
    double time=cpuTime();
    bayesian = false;
    readPedigree(fileName, wcsp);
    if(ToulBar2::haplotype){
        readMap(fileName);
        initTransmission();
        buildWCSP_haplotype(fileName, wcsp);
    }
    else{
        ToulBar2::map_file="";
        buildWCSP(fileName, wcsp);
    }
    cout << "Reading and creating problem time :\t " << cpuTime()-time << endl;
}

void Haplotype::read_bayesian(const char *fileName, WCSP *wcsp )
{
    bayesian = true;
    readPedigree(fileName, wcsp);
    buildWCSP_bayesian(fileName, wcsp);
}

void Haplotype::readMap(const char *fileName){

    Double position;
    bool ok = true;
    ifstream fmap (ToulBar2::map_file.c_str());
    if (!fmap){
        int pos = string(fileName).find_last_of(".");
        string strmap(string(fileName).substr(0,pos) + string(".map"));
        fmap.open(strmap.c_str());
        cerr << "No markers map file specified. Trying " << strmap << endl;
        if(!fmap){
            cerr << "No markers map file found." << endl;
            exit(EXIT_FAILURE);
        }
    }
    while(fmap && ok){

        fmap >> position;
        if(fmap)
            maplocus.push_back(position);
    }
    fmap.close();
    /*	int i=1;
	for(vector<Double>::iterator it=maplocus.begin(); it!=maplocus.end();++it){
		cout << i << " " << *it << " " << maplocus[i-1] << endl;// *it << endl;
		++i;
	}
     */

}

inline bool cmp_generation(Individual i1, Individual i2) { return i1.generation < i2.generation || (i1.generation == i2.generation && i1.individual < i2.individual); }

int Haplotype::fixGenerationNumber(int index)
{

    if (pedigree[index].generation!=-1)
        return pedigree[index].generation;
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
void Haplotype::readPedigree(const char *fileName, WCSP *wcsp)
{
    int individual;
    int nbindividuals = 0;
    int nbtypings = 0;
    //  map<int, int> allelesInv;
    map<int, map<int,int> > allelesInv; // allelesInv[locus,val]=allele
    map<int,int> nballeles;	// nballeles[i] number of allele for the marker i
    map<int,int> maxallele;
    int numal = 0;
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
        int cur_family = -1;
        individual = 0;

        file >> cur_family;
        if (!file) break;
        if (family == -1) family = cur_family;
        if (family != cur_family) {
            cerr << "Pedigree datafile contains more than one locus!" << endl;
            exit(EXIT_FAILURE);
        }

        file >> individual;
        if (!file) {
            cerr << "(1) Wrong data after individual " << individual << endl;
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
            cerr << "(2) Wrong data after individual " << individual << endl;
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
            cerr << "(3) Wrong data after individual " << individual << endl;
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
            cerr << "(4) Wrong data after individual " << individual << endl;
            exit(EXIT_FAILURE);
        }
        numal=0;
        do{
            int allele = 0;
            Genotype gen;
            file >> allele;
            if (!file) { cerr << "(5) Wrong data after individual " << individual << endl; exit(EXIT_FAILURE); }
            if (allele < 0) {
                gen.fixed=true;
                //pedigree[individuals[individual]].genotypes[numal].fixed = true;
                allele = -allele;
            }
            gen.allele1 = allele;
            //pedigree[individuals[individual]].genotypes[numal].allele1 = allele;

            allele = 0;
            file >> allele;
            if (!file) { cerr << "(6) Wrong data after individual " << individual << endl; exit(EXIT_FAILURE); }
            if (allele < 0) {
                gen.fixed=true;
                allele = -allele;
            }
            gen.allele2 = allele;
            pedigree[individuals[individual]].genotypes.push_back(gen);
            numal++;
        }while(file.peek()!='\n');
        nbloci=numal;


        for(int a=0; a<nbloci; ++a){
            int allele1 = pedigree[individuals[individual]].genotypes[a].allele1;
            int allele2 = pedigree[individuals[individual]].genotypes[a].allele2;

            if (alleles.find(a)->second.count(allele1) == 0) {
                nballeles[a]++;
                alleles[a].insert(pair<int,int> (allele1,nballeles[a]));
                //alleles.insert(pair< pair<int,int>,int >(pair<int,int>(a,allele1),nballeles[a]));
                //alleles[allele1] = nballeles;
                freqalleles[a].insert(pair<int,int>(pair<int,int>(allele1,1)));

                if(allele1 > maxallele[a]) maxallele[a] =  allele1;
            }
            else {
                freqalleles[a].find(allele1)->second++;
                //freqalleles.find(pair<int,int> (a,allele1))->second ++;
            }



            if (alleles.find(a)->second.count(allele2) == 0) {
                nballeles[a]++;

                alleles[a].insert(pair<int,int> (allele2,nballeles[a]));
                //alleles.insert(pair<pair<int,int>,int >(pair<int,int>(a,allele2),nballeles[a]));
                //alleles[allele1] = nballeles;
                freqalleles[a].insert(pair<int,int>(pair<int,int>(allele2,1)));
                //freqalleles.insert(pair<pair<int,int>,int >(pair<int,int>(a,allele2),1));
                if(allele2 > maxallele[a]) maxallele[a] =  allele2;
            }
            else {
                freqalleles[a].find(allele2)->second++;
                //freqalleles.find(pair<int,int> (a,allele2))->second ++;
            }
            // A MODIFIER : met dans genotypes les individu qui le sont pour au moins 1 marqueur => les ajouter si tout les marqueurs ou faire une variable genotypes pour chaque marqueur (proposition la plus pertinante?)
            if (allele1>0 || allele2>0) {
                nbtypings++;
                genotypes.push_back(individual);
            }
        }


    }
#ifdef MENDELSOFT
    //  if(ToulBar2::consecutiveAllele) {
    //  	cout << "Make alleles consecutive and supose 4 allele." << endl;
    //    maxallele = 4;
    //  	for(int i=1;i<=maxallele;i++) {
    //  		map<int,int>::iterator it = alleles.find(i);
    //  		if(it == alleles.end()) {
    //	        nballeles++;
    //  			alleles[i] = nballeles;
    // 	        freqalleles[i] = 0;
    //  		}
    //  	}
    //  }
#endif

    /* re-encoding of alleles */
    map<int,int> nb;
    if (ToulBar2::verbose >= 2) cout << "Alleles encoding:" << endl;
    for (map<int,map<int,int> >::iterator iter = alleles.begin(); iter != alleles.end(); ++iter) {
        //if ((*iter).first.second == 0) continue;
        //nb[(*iter).first]++;
        for(map< int, int>::iterator iter2 = iter->second.begin();iter2 != iter->second.end(); ++iter2){
            nb[(*iter).first]++;
            int n= nb.find((*iter).first)->second;
            if (ToulBar2::verbose >= 2) cout << "locus " << (*iter).first << ", " << (*iter2).first /*<< ": " << n*/ << endl;
            iter2->second = n;
            //allelesInv[nb] = (*iter).first;
            //allelesInv.insert(pair<pair<int,int>,int >(pair<int,int>((*iter).first,n),(*iter2).first));
            allelesInv[(*iter).first].insert(pair<int,int>(n,(*iter2).first));
        }
    }
    for(int i=0; i<nbloci; ++i)
        assert(nballeles[i] == nb.find(i)->second);

    if(!ToulBar2::haplotype){
        for(int l=0; l<nbloci;++l){
            if (ToulBar2::verbose >= 1) cout << "Genotype encoding for locus " << l << " :" << endl;
            for (int i=1; i<=nballeles[l]; i++) { /* i = first allele of child */
                for (int j=i; j<=nballeles[l]; j++) { /* j = second allele of child */
                    Genotype geno;
                    geno.allele1 = allelesInv.find(l)->second.find(i)->second;//allelesInv.find(pair<int,int>(l,i))->second;
                    geno.allele2 = allelesInv.find(l)->second.find(j)->second;//allelesInv.find(pair<int,int>(l,j))->second; //[j];
                    genoconvert[l].push_back(geno);
                    if (ToulBar2::verbose >= 1) {
                        cout << genoconvert[l].size()-1 << ": ";
                        printGenotype(cout, genoconvert[l].size()-1,l);
                        cout << endl;
                    }
                }
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

    int l = 0; //  A MODIFIER POUR LA PRISE EN COMPTE DE PLUSIEURS LOCUS
    for (unsigned int i=0; i<genotypes.size(); i++) {
        typeAscendants(genotypes[i]);
        if (genotypes[i]>=0 && pedigree[individuals[genotypes[i]]].father > 0 &&
                pedigree[individuals[pedigree[individuals[genotypes[i]]].father]].genotypes[l].allele1>0 &&
                pedigree[individuals[pedigree[individuals[genotypes[i]]].father]].genotypes[l].allele2>0)
            pedigree[individuals[pedigree[individuals[genotypes[i]]].father]].nbtyped++;
        if (genotypes[i]>=0 && pedigree[individuals[genotypes[i]]].mother > 0 &&
                pedigree[individuals[pedigree[individuals[genotypes[i]]].mother]].genotypes[l].allele1>0 &&
                pedigree[individuals[pedigree[individuals[genotypes[i]]].mother]].genotypes[l].allele2>0)
            pedigree[individuals[pedigree[individuals[genotypes[i]]].mother]].nbtyped++;
    }
    cout << nbtyped << " informative individuals found (either genotyped or having a genotyped descendant)." << endl;


    gencorrects.clear();
    if (fileErrors) {
        while(fileErrors) {
            int bonallele1;
            fileErrors >> individual;
            fileErrors >> bonallele1;
            int bonallele2 = pedigree[individuals[individual]].genotypes[numal].allele2;
            if(bonallele1 == pedigree[individuals[individual]].genotypes[numal].allele1) cout << "error file with a mistake" << endl;
            gencorrects[individual] = convertgen(l,bonallele1, bonallele2);
        }
    }
    assert(wcsp->numberOfVariables() == 0);
    assert(wcsp->numberOfConstraints() == 0);
}

int Haplotype::convertgen( int locus, int allele1, int allele2 ) {
    int nballeles = alleles.find(locus)->second.size() - 1;
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


// A MODIFIER POUR PLUSIEURS LOCUS
void Haplotype::buildWCSP(const char *fileName, WCSP *wcsp)
{
    int locus=0;
    ifstream file(fileName);

    vector<TemporaryUnaryConstraint> unaryconstrs;

    int nbindividuals = individuals.size();
    int nballeles = alleles.find(locus)->second.size() - 1;
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
        allele1 = alleles.find(locus)->second.find(allele1)->second;
        file >> allele2;
        if (allele2 < 0) allele2 = -allele2;
        allele2 = alleles.find(locus)->second.find(allele2)->second;
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
                        unaryconstr.costs.push_back((pedigree[individuals[individual]].genotypes[locus].fixed)?wcsp->getUb():1);
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

    wcsp->histogram();

    if (ToulBar2::verbose >= 0) {
        cout << "Read pedigree with " << nbindividuals << " individuals, " << nbfounders << " founders, " << nballeles << " alleles, " << nbtypings << " genotypings and " << generations << " generations." << endl;
    }

}


// A MODIFIER POUR PLUSIEURS LOCUS
void Haplotype::buildWCSP_bayesian( const char *fileName, WCSP *wcsp )
{

    int locus = 0;
    ifstream file(fileName);

    vector<TemporaryUnaryConstraint> unaryconstrs;
    int nbindividuals = individuals.size();
    int nballeles = alleles.find(locus)->second.size() - 1;
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

    map<int, int> allelesInv; // A MODIFIER POUR PLUSIEURS LOCUS voir readPedigree;

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
    for (map<int,int>::iterator iter = alleles[locus].begin(); iter != alleles[locus].end(); ++iter) {
        allelesInv[iter->second] = iter->first;
    }
    for (i=1; i<=nballeles; i++) { /* i = first allele of child  */
        for (j=i; j<=nballeles; j++) { /* j = second allele of child */
            // case i!=j must be counted twice (i,j and j,i)
            foundersprob.push_back((TProb)((i!=j)?2.:1.) * (TProb)freqalleles[locus].find( allelesInv[i] )->second * freqalleles[locus].find( allelesInv[j] )->second / (TProb)(nbtypings * nbtypings * 4) );
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
            printGenotype(cout, n,locus);
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
        allele1 = alleles.find(locus)->second.find(allele1)->second;
        file >> allele2;
        if (allele2 < 0) allele2 = -allele2;
        allele2 = alleles.find(locus)->second.find(allele2)->second;
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

                    bool fixed = pedigree[individuals[individual]].genotypes[locus].fixed;

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

    wcsp->histogram();

    if (ToulBar2::verbose >= 0) {
        int nbtypings = genotypes.size();
        cout << "Read pedigree with " << nbindividuals << " individuals, " << nbfounders << " founders, " << nballeles << " alleles, " << nbtypings << " genotypings and " << generations << " generations." << endl;
        cout << "Bayesian MPE (genotyping error rate: " <<  ToulBar2::errorg << ", genotype prior: " << ToulBar2::foundersprob_class << ", precision(1-10^-p): " << ToulBar2::resolution << ", normalization: " << ToulBar2::NormFactor << ", ub: " << wcsp->getUb() << ")" << endl;

    }


}

void Haplotype::buildWCSP_haplotype( const char *fileName, WCSP *wcsp )
{
    //create the sparse matrix
    sparse_matrix();


    //create Boolean variables
    for (int i=0; i<nbloci; i++) {
        wcsp->makeEnumeratedVariable(to_string(i), 0, 1);
    }

    vector<Cost> unaryCosts0(nbloci, 0);
    vector<Cost> unaryCosts1(nbloci, 0);

    // find total cost
    Double sumcost = 0.;

    for(map<pair<int,int>,Double,classcomp >::iterator w = W.begin(); w != W.end(); ++w)
        sumcost += abs(w->second);
    ToulBar2::NormFactor = (-1.0/Log1p(-Exp10(-(TProb)ToulBar2::resolution)));
    wcsp->updateUb((Cost) (ToulBar2::NormFactor*sumcost));
    //	  Double constante = 0.;
    for(map<pair<int,int>,Double,classcomp >::iterator w = W.begin(); w != W.end(); ++w){
        if(w->first.first != w->first.second){
            vector<Cost> costs(4,0);
            if(w->second > 0){
                costs[1] = (Cost) ( ToulBar2::NormFactor * w->second);
                costs[2] = costs[1];
                K += 2. * w->second;
                //				  constante += 2. * w->second;
            }else{
                costs[0] = (Cost) ( -ToulBar2::NormFactor * w->second);
                costs[3] = costs[0];
                K += -2. * w->second;
                //				  constante += -2. * w->second;
            }
            if(w->second != 0)
                wcsp->postBinaryConstraint(w->first.first, w->first.second, costs);
        }else{
            //			  if(w->second > 0)
            //				  unaryCosts1[w->first.first] += (Cost) (multiplier * w->second);
            //			  else
            //				  unaryCosts0[w->first.first] += (Cost) (multiplier * -w->second);
        }
    }
    // create weighted unary clauses
    for (int i=0; i<nbloci; i++) {
        if (unaryCosts0[i] > 0 || unaryCosts1[i] > 0) {
            vector<Cost> costs(2, 0);
            costs[0] = unaryCosts0[i];
            costs[1] = unaryCosts1[i];
            wcsp->postUnary(i, costs);
        }
    }
    cout << "Read " << nbloci << " variables, with " << 2 << " values at most, and " << W.size() << " constraints." << endl;
    // special data structure to be initialized for variable ordering heuristics
    if(ToulBar2::verbose == 1) cout << "pedigree ub: " << wcsp->getUb() << endl;

}

// A MODIFIER POUR PLUSIEURS LOCUS
void Haplotype::printCorrectSol(WCSP *wcsp)
{
    int locus=0;
    if(!gencorrects.size()) return;

    ofstream file("sol_correct");
    if (!file) {
        cerr << "Could not write file " << "solution" << endl;
        exit(EXIT_FAILURE);
    }

    for(vector<Individual>::iterator it = pedigree.begin(); it != pedigree.end(); ++it )
    {
        Individual& ind = *it;
        int allele1 = ind.genotypes[locus].allele1;
        int allele2 = ind.genotypes[locus].allele2;
        if((allele1 > 0) || (allele2 > 0)) {
            map<int,int>::iterator it = gencorrects.find(ind.individual);
            if(it != gencorrects.end()) file << " " << it->second;
            else file << " " << convertgen(locus,allele1, allele2);
        }
        else file << " " << -1;
    }
    file << endl;
}


void Haplotype::printSol(WCSP *wcsp)
{

    if(!ToulBar2::haplotype){
        ofstream file("sol");
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
        file.close();
    }
    else{
        ofstream file("haplotypes");
        cout << "haplotypes\n";
        file << "sire " << sire << endl;
        for(int locus = 0; locus < nbloci; ++locus){
            file << " ";
            if(wcsp->getValue(locus)==1) file << pedigree[individuals.find(sire)->second].genotypes[locus].allele1 << " ";
            else file << pedigree[individuals.find(sire)->second].genotypes[locus].allele2 << " ";
        }
        file << endl;
        for(int locus = 0; locus < nbloci; ++locus){
            file << " ";
            if(wcsp->getValue(locus)==0) file << pedigree[individuals.find(sire)->second].genotypes[locus].allele1 << " ";
            else file << pedigree[individuals.find(sire)->second].genotypes[locus].allele2 << " ";
        }
        file << endl;
        file.close();
    }
}


// A MODIFIER POUR PLUSIEURS LOCUS
void Haplotype::printCorrection(WCSP *wcsp)
{
    int locus = 0;
    bool errorinfo = gencorrects.size() > 0;

    int ncorrect = 0;
    int ncorrections = 0;
    int ncorrectok = 0;

    int nbcorrection = 0;

    TProb penalty = 0.;

    cout << "Correction:";
    for (unsigned int i=0; i<genotypes.size(); i++) {
        int sol = wcsp->getValue(pedigree[individuals[genotypes[i]]].varindex);
        int a1 = (genoconvert[locus])[sol].allele1;
        int a2 = (genoconvert[locus])[sol].allele2;
        int allele1 = pedigree[individuals[genotypes[i]]].genotypes[locus].allele1;
        int allele2 = pedigree[individuals[genotypes[i]]].genotypes[locus].allele2;
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


void Haplotype::printGenotype(ostream& os, Value value, int locus)
{
    os << (genoconvert[locus])[value].allele1 << "/" << (genoconvert[locus])[value].allele2;
}

void Haplotype::printHaplotype(ostream& os, Value value, int locus)
{
    if(value == 1) os << pedigree[individuals.find(sire)->second].genotypes[locus].allele1 << "|" << pedigree[individuals.find(sire)->second].genotypes[locus].allele2 << " ";
    else os << pedigree[individuals.find(sire)->second].genotypes[locus].allele2 << "|" << pedigree[individuals.find(sire)->second].genotypes[locus].allele1 << " ";
}

// A MODIFIER POUR PLUSIEURS LOCUS
void Haplotype::save(const char *fileName, WCSP *wcsp, bool corrected, bool reduced)
{
    int locus = 0;
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
        file << family << " ";
        if (corrected && pedigree[(*iter).second].varindex >= 0 && pedigree[(*iter).second].varindex < (int) wcsp->numberOfVariables() && wcsp->assigned(pedigree[(*iter).second].varindex)) {
            int sol = wcsp->getValue(pedigree[(*iter).second].varindex);
            int a1 = (genoconvert[locus])[sol].allele1;
            int a2 = (genoconvert[locus])[sol].allele2;
            int allele1 = pedigree[(*iter).second].genotypes[locus].allele1;
            int allele2 = pedigree[(*iter).second].genotypes[locus].allele2;
            if ((allele1>0 && allele2>0 && !((a1==allele1 && a2==allele2) || (a1==allele2 && a2==allele1)))
                    || ((allele1==0 || allele2==0) && !(allele1==0 && allele2==0) && !(a1==allele1 || a1==allele2 || a2==allele1 || a2==allele2))
                    || ToulBar2::pedigreeCorrectionMode == 2) {
                if (ToulBar2::pedigreeCorrectionMode > 0) {
                    pedigree[(*iter).second].genotypes[locus].allele1 = a1;
                    pedigree[(*iter).second].genotypes[locus].allele2 = a2;
                } else {
                    pedigree[(*iter).second].genotypes[locus].allele1 = 0;
                    pedigree[(*iter).second].genotypes[locus].allele2 = 0;
                }
            }
            pedigree[(*iter).second].print(file);
            pedigree[(*iter).second].genotypes[locus].allele1 = allele1;
            pedigree[(*iter).second].genotypes[locus].allele2 = allele2;
        } else {
            pedigree[(*iter).second].print(file);
        }
    }
}

void Haplotype::initTransmission()
{
    assert(generations == 2);
    sire = -1;
    vector<int> sons;
    // search the sire
    for(vector<Individual>::iterator ind=pedigree.begin(); ind != pedigree.end(); ++ind)
    {
        if((*ind).father == 0 && (*ind).mother == 0 && (*ind).sex == 1)
            sire = (*ind).individual;
    }
    assert(sire != -1);

    if(ToulBar2::verbose >= 1) cout << "sire is individual " << sire << endl;
    for(vector<Individual>::iterator fils=pedigree.begin(); fils != pedigree.end(); ++fils)
    {
        vector<int> T;
        if((*fils).father == sire){
            sons.push_back((*fils).individual);
            int ind_mother = (*fils).mother;
            Individual father = pedigree[individuals[(*fils).father]];
            for(int locus = 0; locus < nbloci; ++locus){
                int trans=5;
                if(father.genotypes[locus].allele1 == father.genotypes[locus].allele2) trans=0;//T.push_back(0);
                else
                {
                    if((*fils).genotypes[locus].allele1 == (*fils).genotypes[locus].allele2){
                        if((*fils).genotypes[locus].allele1 == father.genotypes[locus].allele1)
                            trans = -1;//T.push_back(-1);
                        else trans = 1;//T.push_back(-1);
                    }
                    else
                        if(ind_mother != 0){
                            Individual mother = pedigree[individuals[(*fils).mother]];
                            if(mother.genotypes[locus].allele1 != 0 &&  mother.genotypes[locus].allele2 != 0){ // genotyped mother at this locus
                                if( mother.genotypes[locus].allele1 != mother.genotypes[locus].allele2) trans=0;//T.push_back(0);
                                else{
                                    if( (*fils).genotypes[locus].allele1 != mother.genotypes[locus].allele1){
                                        if((*fils).genotypes[locus].allele1 == father.genotypes[locus].allele1)
                                            trans = -1; //T.push_back(-1);
                                        else trans = 1;//T.push_back(-1);
                                    }
                                    else if( (*fils).genotypes[locus].allele2 != mother.genotypes[locus].allele1){
                                        if((*fils).genotypes[locus].allele2 == father.genotypes[locus].allele1)
                                            trans = -1;//T.push_back(-1);
                                        else trans = 1
                                                ;//T.push_back(-1);
                                    }

                                }
                            }
                            else trans = 0;
                        }
                        else trans = 0;

                }
                T.push_back(trans);
            }
            transmission.insert(pair<int,vector<int> >(fils->individual,T));
        }
    }
    if(ToulBar2::verbose > 1){
        cout << "Transmission vectors : \n";
        for(map<int, vector<int> >::iterator it=transmission.begin(); it!=transmission.end(); ++it){
            cout << "offspring " << it->first << ":\t ";
            for(vector<int>::iterator it2 = it->second.begin();it2 != it->second.end();++it2){
                if((*it2) == 0) cout << " * ";
                if((*it2) == 1) cout << " " << *it2 << " " ;
                if((*it2) == -1) cout << *it2 << " " ;
            }
            cout << endl;
        }
    }
    /*
	ofstream file("transmission");
	for(map<int, vector<int> >::iterator it=transmission.begin(); it!=transmission.end(); ++it){
		file << "offspring " << it->first << ":\t ";
		for(vector<int>::iterator it2 = it->second.begin();it2 != it->second.end();++it2){
			if((*it2) == 0 && pedigree[individuals[sons[it->first].father]].genotypes[*it2].allele1 == pedigree[individuals[sons[it->first].father]].genotypes[*it2].allele2) file << "+ ";
			else if((*it2) == 0) file << " * ";
			if((*it2) == 1) file << " " << *it2 << " " ;
			if((*it2) == -1) file << *it2 << " " ;
		}
		file << endl;
	}
     */

}

void Haplotype::sparse_matrix()
{
    int nbDesc = 0;
    for(vector<Individual>::iterator fils=pedigree.begin(); fils != pedigree.end(); ++fils)
    {
        if((*fils).father == sire){
            int locus_prec=-1;		// first left informative locus with respect to current locus
            bool first=true; // if first informative locus
            nbDesc++;
            for(int locus=0; locus<nbloci; ++locus){
                if(first && transmission.find(fils->individual)->second[locus] != 0){
                    first=false;
                    locus_prec=locus;
                }
                else
                    if(!first && transmission.find(fils->individual)->second[locus] != 0){
                        Double recombination_frac = haldane(maplocus[locus]-maplocus[locus_prec]);
                        Double coef = 0.25*log((1-recombination_frac)/recombination_frac);
                        K += log((1-recombination_frac)*recombination_frac);
                        if(transmission.find(fils->individual)->second[locus] == transmission.find(fils->individual)->second[locus_prec])
                        {
                            if(W.count(pair<int,int>(locus_prec,locus)) == 0)
                                W.insert(pair<pair<int,int>,Double >(pair<int,int>(locus_prec,locus),coef ));
                            else
                                W.find(pair<int,int>(locus_prec,locus))->second+=coef;
                        }
                        else
                        {
                            if(W.count(pair<int,int>(locus_prec,locus)) == 0)
                                W.insert(pair<pair<int,int>,Double >(pair<int,int>(locus_prec,locus),-coef ));
                            else
                                W.find(pair<int,int>(locus_prec,locus))->second+=-coef;
                            locus_prec=locus;
                        }
                        locus_prec=locus;
                    }
            }
        }
    }
    K=0.5*K+nbDesc*log(0.5);
    //affichage
    if(ToulBar2::verbose >= 1){
        cout << "sparse matrix : \n";
        for(map< pair<int,int>,Double,classcomp >::iterator it = W.begin(); it != W.end(); ++it)
            if((*it).second>=0)
                cout << "W" <<(*it).first.first << "," << (*it).first.second << "\t =\t  " << (*it).second  << endl;
            else
                cout << "W" <<(*it).first.first << "," << (*it).first.second << "\t =\t " << (*it).second  << endl;
        cout << "constant K =\t " << K << endl;
    }
}

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

