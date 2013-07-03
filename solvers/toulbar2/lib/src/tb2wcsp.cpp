/*
 * ****** Weighted constraint satisfaction problem modeling and local reasoning ********
 *
 * Contains also ToulBar2 options expressed by global variable definitions
 */

#include "tb2wcsp.hpp"
#include "tb2enumvar.hpp"
#include "tb2intervar.hpp"
#include "tb2binconstr.hpp"
#include "tb2ternaryconstr.hpp"
#include "tb2naryconstr.hpp"
#include "tb2arithmetic.hpp"
#include "tb2pedigree.hpp"
#include "tb2haplotype.hpp"
#include "tb2vac.hpp"
#include "tb2clusters.hpp"
#include "ToulbarVersion.hpp"

#include "tb2globaldecomposable.hpp"
#include "tb2globalconstr.hpp"
#include "tb2flowbasedconstr.hpp"
#include "tb2alldiffconstr.hpp"
#include "tb2globalcardinalityconstr.hpp"
#include "tb2sameconstr.hpp"
#include "tb2regularconstr.hpp"

/*
 * Global variables with their default value
 *
 */

string ToulBar2::version = Toulbar_VERSION;
int ToulBar2::verbose = 0;
int ToulBar2::debug = 0;
bool ToulBar2::showSolutions = false;
bool ToulBar2::writeSolution = false;
bool ToulBar2::allSolutions = false;
int ToulBar2::dumpWCSP = 0;
bool ToulBar2::approximateCountingBTD = false;
int ToulBar2::elimDegree = 3;
int ToulBar2::elimDegree_preprocessing = -1;
int ToulBar2::elimDegree_ = -1;
int ToulBar2::elimDegree_preprocessing_ = -1;
int ToulBar2::elimSpaceMaxMB = 0;
int ToulBar2::preprocessTernaryRPC = 0;
int ToulBar2::preprocessFunctional = 1;
bool ToulBar2::costfuncSeparate = true;
int ToulBar2::preprocessNary = 10;
LcLevelType ToulBar2::LcLevel = LC_EDAC;
bool ToulBar2::QueueComplexity = false;
bool ToulBar2::binaryBranching = true;
bool ToulBar2::lastConflict = true;
int ToulBar2::dichotomicBranching = 1;
unsigned int ToulBar2::dichotomicBranchingSize = 10;
bool ToulBar2::sortDomains = false;
map<int, ValueCost *> ToulBar2::sortedDomains;
int ToulBar2::lds = 0;
bool ToulBar2::limited = false;
Long ToulBar2::restart = -1;
bool ToulBar2::generation = false;
int ToulBar2::minsumDiffusion = 0;
bool ToulBar2::Static_variable_ordering=false;
int ToulBar2::weightedDegree = 10000;
int ToulBar2::weightedTightness = 0;
bool ToulBar2::MSTDAC = false;
int ToulBar2::DEE = 1;
int ToulBar2::nbDecisionVars = 0;
bool ToulBar2::singletonConsistency = false;
bool ToulBar2::vacValueHeuristic = false;

externalevent ToulBar2::setvalue = NULL;
externalevent ToulBar2::setmin = NULL;
externalevent ToulBar2::setmax = NULL;
externalevent ToulBar2::removevalue = NULL;
externalcostevent ToulBar2::setminobj = NULL;
externalsolution ToulBar2::newsolution = NULL;
Pedigree *ToulBar2::pedigree = NULL;
Haplotype *ToulBar2::haplotype = NULL;

bool ToulBar2::bayesian = false;
bool ToulBar2::uai = false;
string ToulBar2::evidence_file;
ofstream ToulBar2::solution_file;
bool ToulBar2::uai_firstoutput = true;
TProb ToulBar2::markov_log = 0;
bool ToulBar2::xmlflag = false;
string ToulBar2::map_file;
bool ToulBar2::maxsateval = false;

int ToulBar2::resolution = 7;
TProb ToulBar2::errorg = 0.05;
TProb ToulBar2::NormFactor = 1;
/// Allele frequencies of founders
/// - 0: 			equal frequencies
/// - 1: 			probs depending on the frequencies found in the problem
/// - otherwise:  read probability distribution from command line
int ToulBar2::foundersprob_class = 0;
vector<TProb> ToulBar2::allelefreqdistrib;
bool ToulBar2::consecutiveAllele = false;
int ToulBar2::pedigreeCorrectionMode = 0;
int ToulBar2::pedigreePenalty = 0;

int ToulBar2::vac = 0;
Cost ToulBar2::costThreshold = UNIT_COST;
Cost ToulBar2::costThresholdPre = UNIT_COST;
double ToulBar2::costMultiplier = UNIT_COST;
Cost ToulBar2::relaxThreshold = MIN_COST;

ElimOrderType ToulBar2::elimOrderType = ELIM_NONE;

BEP *ToulBar2::bep = NULL;
bool ToulBar2::wcnf = false;
bool ToulBar2::qpbo = false;

int WCSP::wcspCounter = 0;

char* ToulBar2::varOrder = NULL;
int ToulBar2::btdMode = 0;
int ToulBar2::btdSubTree = -1;
int ToulBar2::btdRootCluster = -1;

double ToulBar2::startCpuTime = 0;

int ToulBar2::splitClusterMaxSize = 0;
bool ToulBar2::boostingBTD = false;
int ToulBar2::maxSeparatorSize = -1;
int ToulBar2::minProperVarSize = 0;

int ToulBar2::smallSeparatorSize = 4;

bool ToulBar2::isZ = false;
TProb ToulBar2::logZ = -numeric_limits<TProb>::infinity();
TProb ToulBar2::logU = -numeric_limits<TProb>::infinity();
TProb ToulBar2::logepsilon = -3;
int ToulBar2::Berge_Dec=0; // berge decomposition flag  > 0 if wregular found in the problem
int ToulBar2::nbvar=0; // berge decomposition flag  > 0 if wregular found in the problem

externalfunc ToulBar2::timeOut = NULL;
bool ToulBar2::interrupted = false;

/*
 * WCSP constructors
 *
 */

WCSP::WCSP(Store *s, Cost upperBound, void *_solver_) :
	solver(_solver_), storeData(s), lb(MIN_COST, &s->storeCost), ub(upperBound), negCost(MIN_COST, &s->storeCost), NCBucketSize(cost2log2gub(upperBound) + 1),
			NCBuckets(NCBucketSize, VariableList(&s->storeVariable)), PendingSeparator(&s->storeSeparator),
			objectiveChanged(false), nbNodes(0), nbDEE(0), lastConflictConstr(NULL), maxdomainsize(0), isDelayedNaryCtr(true),
			elimOrder(0, &s->storeInt), elimBinOrder(0, &s->storeInt), elimTernOrder(0, &s->storeInt),
	        maxDegree(-1), elimSpace(0) {
	instance = wcspCounter++;
	if (ToulBar2::vac) vac = new VACExtension(this);
	else vac = NULL;

	td = NULL;
}

WCSP::~WCSP() {
	if (vars.size()) for (unsigned int i = 0; i < vars.size(); i++)
		delete vars[i];
	if (constrs.size()) for (unsigned int i = 0; i < constrs.size() - 1; i++)
		delete constrs[i]; // Warning! The last constraint may be badly allocated due to an exception occuring in its constructor (because of propagate) // If there is no constraint then (constrs.size()-1) overflow!
	if (elimBinConstrs.size()) for (unsigned int i = 0; i < elimBinConstrs.size(); i++)
		delete elimBinConstrs[i];
	if (elimTernConstrs.size()) for (unsigned int i = 0; i < elimTernConstrs.size(); i++)
		delete elimTernConstrs[i];
}

WeightedCSP *WeightedCSP::makeWeightedCSP(Store *s, Cost upperBound, void *solver) {
	WeightedCSP * W = new WCSP(s, upperBound, solver);
	return W;
}

/// \brief create an enumerated variable with its domain bounds
int WCSP::makeEnumeratedVariable(string n, Value iinf, Value isup) {
	EnumeratedVariable *x;
	if (!ToulBar2::vac) {
		x = new EnumeratedVariable(this, n, iinf, isup);
	} else {
		x = new VACVariable(this, n, iinf, isup);
	}
	if (maxdomainsize < isup - iinf + 1) maxdomainsize = isup - iinf + 1;
	listofsuccessors.push_back(vector<int>()); // add new variable in the topological order list;
	return x->wcspIndex;
}

/// \brief create an enumerated variable with its domain values
int WCSP::makeEnumeratedVariable(string n, Value *d, int dsize) {
	EnumeratedVariable *x;
	if (!ToulBar2::vac) {
		x = new EnumeratedVariable(this, n, d, dsize);
	} else {
		x = new VACVariable(this, n, d, dsize);
	}
	if (maxdomainsize < dsize) maxdomainsize = dsize;
	listofsuccessors.push_back(vector<int>()); // add new variable in the topological order list;
	return x->wcspIndex;
}

/// \brief create an interval variable with its domain bounds
int WCSP::makeIntervalVariable(string n, Value iinf, Value isup) {
	IntervalVariable *x = new IntervalVariable(this, n, iinf, isup);
	listofsuccessors.push_back(vector<int>()); // add new variable in the topological order list;
	return x->wcspIndex;
}

/// \brief create a binary cost function from a vector of costs
/// \param xIndex index of enumerated variable x as returned by makeEnumeratedVariable
/// \param yIndex index of enumerated variable y
/// \param costs a flat vector of costs (y indexes moving first)
///
/// \note It looks for an existing constraint
/// (even not connected constraints). It also allocates
/// memory for a new constraint. For this two reasons it should
/// ONLY be called before search.
///
/// \warning Vector costs must have the same size as Cartesian product of original domains.
int WCSP::postBinaryConstraint(int xIndex, int yIndex, vector<Cost> &costs) {
	EnumeratedVariable* x = (EnumeratedVariable *) vars[xIndex];
	EnumeratedVariable* y = (EnumeratedVariable *) vars[yIndex];

    if(ToulBar2::vac) {
        for (unsigned int a = 0; a < x->getDomainInitSize(); a++) {
            for (unsigned int b = 0; b < y->getDomainInitSize(); b++) {
     			Cost c = costs[a * y->getDomainInitSize() + b];
                histogram(c);
            }
        }
    }

	BinaryConstraint* ctr = x->getConstr(y);
	if (ctr) {
		ctr->reconnect();
		ctr->addCosts(x, y, costs);
		ctr->propagate();
	} else {
		if (!ToulBar2::vac) {
			ctr
					= new BinaryConstraint(this, (EnumeratedVariable *) vars[xIndex], (EnumeratedVariable *) vars[yIndex], costs, &storeData->storeCost);
		} else {
			ctr
					= new VACBinaryConstraint(this, (EnumeratedVariable *) vars[xIndex], (EnumeratedVariable *) vars[yIndex], costs, &storeData->storeCost);
		}
	}

	return ctr->wcspIndex;
}

/// \brief create a ternary cost function from a flat vector of costs (z indexes moving first)
int WCSP::postTernaryConstraint(int xIndex, int yIndex, int zIndex, vector<Cost> &costs) {
	EnumeratedVariable* x = (EnumeratedVariable *) vars[xIndex];
	EnumeratedVariable* y = (EnumeratedVariable *) vars[yIndex];
	EnumeratedVariable* z = (EnumeratedVariable *) vars[zIndex];

	if(ToulBar2::vac) {
        for (unsigned int a = 0; a < x->getDomainInitSize(); a++) {
            for (unsigned int b = 0; b < y->getDomainInitSize(); b++) {
                for (unsigned int c = 0; c < z->getDomainInitSize(); c++) {
         			Cost co = costs[a * y->getDomainInitSize() * z->getDomainInitSize() + b * z->getDomainInitSize() + c];
                    histogram(co);
                }
            }
        }
    }

	TernaryConstraint* ctr = x->getConstr(y, z);

	if (!ctr) {
		unsigned int a, b;
		vector<Cost> zerocostsxy;
		vector<Cost> zerocostsxz;
		vector<Cost> zerocostsyz;

		for (a = 0; a < x->getDomainInitSize(); a++) {
			for (b = 0; b < y->getDomainInitSize(); b++) {
				zerocostsxy.push_back(MIN_COST);
			}
		}
		for (a = 0; a < x->getDomainInitSize(); a++) {
			for (b = 0; b < z->getDomainInitSize(); b++) {
				zerocostsxz.push_back(MIN_COST);
			}
		}
		for (a = 0; a < y->getDomainInitSize(); a++) {
			for (b = 0; b < z->getDomainInitSize(); b++) {
				zerocostsyz.push_back(MIN_COST);
			}
		}

		BinaryConstraint* xy = x->getConstr(y);
		BinaryConstraint* xz = x->getConstr(z);
		BinaryConstraint* yz = y->getConstr(z);

		if (!ToulBar2::vac) {
			if (!xy) {
				xy = new BinaryConstraint(this, x, y, zerocostsxy, &storeData->storeCost);
				xy->deconnect(true);
			}
			if (!xz) {
				xz = new BinaryConstraint(this, x, z, zerocostsxz, &storeData->storeCost);
				xz->deconnect(true);
			}
			if (!yz) {
				yz = new BinaryConstraint(this, y, z, zerocostsyz, &storeData->storeCost);
				yz->deconnect(true);
			}
		} else {
			if (!xy) {
				xy = new VACBinaryConstraint(this, x, y, zerocostsxy, &storeData->storeCost);
				xy->deconnect(true);
			}
			if (!xz) {
				xz = new VACBinaryConstraint(this, x, z, zerocostsxz, &storeData->storeCost);
				xz->deconnect(true);
			}
			if (!yz) {
				yz = new VACBinaryConstraint(this, y, z, zerocostsyz, &storeData->storeCost);
				yz->deconnect(true);
			}
		}

		ctr = new TernaryConstraint(this, x, y, z, xy, xz, yz, costs, &storeData->storeCost);
	} else ctr->addCosts(x, y, z, costs);

	return ctr->wcspIndex;
}

/// \brief create a global cost function in extension using a default cost (tuples with a different cost will be enter later using WCSP::postNaryConstraintTuple)
/// \param scopeIndex array of enumerated variable indexes (as returned by makeEnumeratedVariable)
/// \param arity size of scopeIndex
/// \param defval default cost for any tuple
/// \note should not be used for binary or ternary cost functions
/// \warning do not forget to initially propagate the global cost function using WCSP::postNaryConstraintEnd
int WCSP::postNaryConstraintBegin(int* scopeIndex, int arity, Cost defval) {
	EnumeratedVariable** scopeVars = new EnumeratedVariable*[arity];
	for (int i = 0; i < arity; i++)
		scopeVars[i] = (EnumeratedVariable *) vars[scopeIndex[i]];
	NaryConstraint *ctr = new NaryConstraintMap(this, scopeVars, arity, defval);
	delete[] scopeVars;

	if (isDelayedNaryCtr) delayedNaryCtr.push_back(ctr->wcspIndex);
	else {
		BinaryConstraint* bctr;
		TernaryConstraint* tctr = new TernaryConstraint(this, &storeData->storeCost);
		elimTernConstrs.push_back(tctr);
		for (int j = 0; j < 3; j++) {
			if (!ToulBar2::vac) bctr = new BinaryConstraint(this, &storeData->storeCost);
			else bctr = new VACBinaryConstraint(this, &storeData->storeCost);
			elimBinConstrs.push_back(bctr);
		}
	}
	return ctr->wcspIndex;
}

/// \brief set one tuple with a specific cost for global cost function in extension
/// \param ctrindex index of cost function as returned by WCSP::postNaryConstraintBegin
/// \param tuple array of values assigned to variables ordered by following the original scope order
/// \param arity size of the array
/// \param cost new cost for this tuple
/// \warning valid only for global cost function in extension
void WCSP::postNaryConstraintTuple(int ctrindex, Value* tuple, int arity, Cost cost) {
	if(ToulBar2::vac) histogram(cost);
	Constraint *ctr = getCtr(ctrindex);
	assert(ctr->extension()); // must be an NaryConstraint
	int indextuple[arity];
	for (int i=0; i<arity; i++) indextuple[i] = ((EnumeratedVariable *) ctr->getVar(i))->toIndex(tuple[i]);
	ctr->setTuple(indextuple, cost, NULL);
}

/// \brief set one tuple with a specific cost for global cost function in extension
/// \param ctrindex index of cost function as returned by WCSP::postNaryConstraintBegin
/// \param tuple String encoding of values (indexes plus CHAR_FIRST) assigned to variables ordered by following the original scope order
/// \param cost new cost for this tuple
/// \warning valid only for global cost function in extension
/// \warning string encoding of tuples is for advanced users only!
void WCSP::postNaryConstraintTuple(int ctrindex, String& tuple, Cost cost) {
	if(ToulBar2::vac) histogram(cost);
	Constraint *ctr = getCtr(ctrindex);
	assert(ctr->extension()); // must be an NaryConstraint
	ctr->setTuple(tuple, cost, NULL);
}

/// \brief create a global cost function in intension with a particular semantic
/// \param scopeIndex array of enumerated variable indexes (as returned by makeEnumeratedVariable)
/// \param arity size of scopeIndex
/// \param gcname specific \e keyword name of the global cost function (\e eg salldiff, sgcc, sregular, ssame)
/// \param file problem file (\see \ref wcspformat)
int WCSP::postGlobalConstraint(int* scopeIndex, int arity, string &gcname, istream &file) {
	assert(arity >= 4); // does not work for binary or ternary cost functions!!!
	if (ToulBar2::verbose >= 2) cout << "Number of global constraints = " << globalconstrs.size() << endl;
	GlobalConstraint* gc = NULL;
	EnumeratedVariable** scopeVars = new EnumeratedVariable*[arity];
	for (int i = 0; i < arity; i++)
		scopeVars[i] = (EnumeratedVariable *) vars[scopeIndex[i]];

	if (gcname == "salldiff") {
		gc = new AllDiffConstraint(this, scopeVars, arity);
	} else if (gcname == "sgcc") {
		gc = new GlobalCardinalityConstraint(this, scopeVars, arity);
	} else if (gcname == "ssame") {
		gc = new SameConstraint(this, scopeVars, arity);
	} else if (gcname == "sregular") {
		gc = new RegularConstraint(this, scopeVars, arity);
	} else {
		cout << gcname << " undefined" << endl;
		exit(1);
	}
	if (gc != NULL) globalconstrs.push_back(gc);
	if (file != NULL) gc->read(file);

	return gc->wcspIndex;
}

void WCSP::postWSum(int* scopeIndex, int arity, string semantics, Cost baseCost, string comparator, int rightRes){
    string gcname = "wsum";
    WeightedSum* decomposableGCF = new WeightedSum(arity, scopeIndex);
    decomposableGCF->setSemantics(semantics); 				
    decomposableGCF->setBaseCost(baseCost); 				
    decomposableGCF->setComparator(comparator); 
	decomposableGCF->setRightRes(rightRes);
    decomposableGCF->addToCostFunctionNetwork(this);
}

void WCSP::postWVarSum(int* scopeIndex, int arity, string semantics, Cost baseCost, string comparator, int varIndex)
{
    string gcname = "wvarsum";
    WeightedVarSum* decomposableGCF = new WeightedVarSum(arity, scopeIndex);
    decomposableGCF->setSemantics(semantics);
    decomposableGCF->setBaseCost(baseCost);
    decomposableGCF->setComparator(comparator);
	decomposableGCF->setIndex(varIndex);
    decomposableGCF->addToCostFunctionNetwork(this);
}

void WCSP::postWAmong(int* scopeIndex, int arity, string semantics, Cost baseCost, Value* values, int nbValues, int lb, int ub)
{
    WeightedAmong* decomposableGCF = new WeightedAmong(arity, scopeIndex);
    decomposableGCF->setSemantics(semantics); 				
    decomposableGCF->setBaseCost(baseCost); 
    for(int i = 0; i < nbValues; i++)
    {
        decomposableGCF->addValue(values[i]);
    }
    decomposableGCF->setBounds(lb,ub);
    decomposableGCF->addToCostFunctionNetwork(this); 
    delete [] values;   
    //delete [] decomposableGCF;
}

void WCSP::postWVarAmong(int* scopeIndex, int arity, string semantics, Cost baseCost, Value* values, int nbValues, int varIndex)
{
    WeightedVarAmong* decomposableGCF = new WeightedVarAmong(arity, scopeIndex);
    decomposableGCF->setSemantics(semantics);
    decomposableGCF->setBaseCost(baseCost);
    for(int i = 0; i < nbValues; i++)
    {
        decomposableGCF->addValue(values[i]);
    }
    decomposableGCF->setIndex(varIndex);
    decomposableGCF->addToCostFunctionNetwork(this);
    delete [] values;
}

void WCSP::postWRegular(int* scopeIndex, int arity, int nbStates, vector<pair<int, Cost> > initial_States, vector<pair<int, Cost> > accepting_States, int** Wtransitions,
  vector<Cost> transitionsCosts)
{
    WFA* automaton=new WFA(nbStates);
    for(unsigned int i=0; i<initial_States.size(); i++)
    {
      automaton->getInitialStates().push_back(initial_States[i]);
    }
    for(unsigned int i=0; i<accepting_States.size(); i++)
    {
      automaton->getAcceptingStates().push_back(accepting_States[i]);
    }
    for(unsigned int i=0; i<transitionsCosts.size(); i++)
    {
      automaton->getTransitions().push_back(new WTransition(Wtransitions[i][0],Wtransitions[i][1],Wtransitions[i][2],transitionsCosts[i]));
    }
    WeightedRegular regular(arity,scopeIndex);
    regular.setWFA(automaton);
    regular.addToCostFunctionNetwork(this);   
    for(unsigned int i=0; i < transitionsCosts.size(); i++)
    {
      delete [] Wtransitions[i];    
    }
    delete [] Wtransitions;
}

void WCSP::postWGcc(int* scopeIndex, int arity, string semantics, Cost baseCost, Value* values, int nbValues, int* lb, int* ub)
{
    WeightedGcc* decomposableGCF = new WeightedGcc(arity, scopeIndex);
    decomposableGCF->setSemantics(semantics);
    decomposableGCF->setBaseCost(baseCost);
    decomposableGCF->setNbValue(nbValues);
    for(int i = 0; i < nbValues; i++)
    {
        decomposableGCF->setBounds(values[i],lb[i],ub[i]);
    }
    decomposableGCF->addToCostFunctionNetwork(this);
    delete [] values;
    delete [] lb;
    delete [] ub;
}

void WCSP::postWSame(int* scopeIndex, int arity, string semantics, Cost baseCost)
{
    WeightedSame* decomposableGCF = new WeightedSame(arity, scopeIndex);
    decomposableGCF->setSemantics(semantics);
    decomposableGCF->setBaseCost(baseCost);
    decomposableGCF->addToCostFunctionNetwork(this);
}

void WCSP::postWAllDiff(int* scopeIndex, int arity, string semantics, Cost baseCost)
{
    WeightedAllDifferent* decomposableGCF = new WeightedAllDifferent(arity, scopeIndex);
    decomposableGCF->setSemantics(semantics);
    decomposableGCF->setBaseCost(baseCost);
    decomposableGCF->addToCostFunctionNetwork(this);
}

void WCSP::postWSameGcc(int* scopeIndex, int arity, string semantics, Cost baseCost, Value* values, int nbValues, int* lb, int* ub)
{
    WeightedSameGcc* decomposableGCF = new WeightedSameGcc(arity, scopeIndex);
    decomposableGCF->setSemantics(semantics);
    decomposableGCF->setBaseCost(baseCost);
    decomposableGCF->setNbValue(nbValues);
    for(int i = 0; i < nbValues; i++)
    {
        decomposableGCF->setBounds(values[i],lb[i],ub[i]);
    }
    decomposableGCF->addToCostFunctionNetwork(this);
    delete [] values;
    delete [] lb;
    delete [] ub;
}

void WCSP::postWOverlap(int* scopeIndex, int arity, string semantics, Cost baseCost, string comparator, int rightRes)
{
    WeightedOverlap* decomposableGCF = new WeightedOverlap(arity, scopeIndex);
    decomposableGCF->setSemantics(semantics);
    decomposableGCF->setBaseCost(baseCost);
    decomposableGCF->setComparator(comparator);
	decomposableGCF->setRightRes(rightRes);
    decomposableGCF->addToCostFunctionNetwork(this);
}


/// \brief add unary costs to enumerated variable \e xIndex
/// \note a unary cost function associated to an enumerated variable is not a Constraint object, it is directly managed inside the EnumeratedVariable class, this is why this function does not return any Constraint index. By doing so, unary costs are better shared inside the cost function network.
void WCSP::postUnary(int xIndex, vector<Cost> &costs) {
	assert(vars[xIndex]->enumerated());
	EnumeratedVariable* x = (EnumeratedVariable *) vars[xIndex];

	if(ToulBar2::vac) {
		for (unsigned int a = 0; a < x->getDomainInitSize(); a++) {
			Cost c = costs[a];
			histogram(c);
		}
	}
	for (unsigned int a = 0; a < x->getDomainInitSize(); a++) {
		if (costs[a] > MIN_COST) x->project(x->toValue(a), costs[a]);
	}
	x->findSupport();
    x->queueNC();
}

/// \brief add unary costs to interval variable \e xIndex
int WCSP::postUnary(int xIndex, Value *d, int dsize, Cost penalty) {
	assert(!vars[xIndex]->enumerated());
	Unary *ctr = new Unary(this, (IntervalVariable *) vars[xIndex], d, dsize, penalty, &storeData->storeValue);
	return ctr->wcspIndex;
}

/// \brief create a soft constraint \f$x \geq y + cst\f$ with the associated cost function \f$max( (y + cst - x \leq delta)?(y + cst - x):top , 0 )\f$
int WCSP::postSupxyc(int xIndex, int yIndex, Value cst, Value delta) {
  if (!vars[xIndex]->enumerated() && !vars[yIndex]->enumerated()) {
	Supxyc
			*ctr =
					new Supxyc(this, (IntervalVariable *) vars[xIndex], (IntervalVariable *) vars[yIndex], cst, delta, &storeData->storeCost, &storeData->storeValue);
	return ctr->wcspIndex;
  } else if (vars[xIndex]->enumerated() && vars[yIndex]->enumerated()) {
	EnumeratedVariable *x = (EnumeratedVariable *) vars[xIndex];
	EnumeratedVariable *y = (EnumeratedVariable *) vars[yIndex];
	vector<Cost> costs;
	for (unsigned int a = 0; a <  x->getDomainInitSize(); a++) {
	  for (unsigned int b = 0; b < y->getDomainInitSize(); b++) {
		costs.push_back(max( (y->toValue(b) + cst - x->toValue(a) <= delta)?((Cost)(y->toValue(b) + cst - x->toValue(a))):getUb(), MIN_COST ));
	  }
	}
	return postBinaryConstraint(xIndex, yIndex, costs);
  } else {
	cerr << "Cannot mix variables with interval and enumerated domains!!!" << endl;
	exit(EXIT_FAILURE);
  }
}

/// \brief soft binary disjunctive constraint \f$x \geq y + csty \vee y \geq x + cstx\f$ with associated cost function \f$(x \geq y + csty \vee y \geq x + cstx)?0:penalty\f$
int WCSP::postDisjunction(int xIndex, int yIndex, Value cstx, Value csty, Cost penalty) {
  if (!vars[xIndex]->enumerated() && !vars[yIndex]->enumerated()) {
	Disjunction
			*ctr =
					new Disjunction(this, (IntervalVariable *) vars[xIndex], (IntervalVariable *) vars[yIndex], cstx, csty, penalty, &storeData->storeValue);
	return ctr->wcspIndex;
  } else if (vars[xIndex]->enumerated() && vars[yIndex]->enumerated()) {
	EnumeratedVariable *x = (EnumeratedVariable *) vars[xIndex];
	EnumeratedVariable *y = (EnumeratedVariable *) vars[yIndex];
	vector<Cost> costs;
	for (unsigned int a = 0; a < x->getDomainInitSize(); a++) {
	  for (unsigned int b = 0; b < y->getDomainInitSize(); b++) {
		costs.push_back(((x->toValue(a) >= y->toValue(b) + csty) || (y->toValue(b) >= x->toValue(a) + cstx))?MIN_COST:penalty);
	  }
	}
	return postBinaryConstraint(xIndex, yIndex, costs);
  } else {
	cerr << "Cannot mix variables with interval and enumerated domains!!!" << endl;
	exit(EXIT_FAILURE);
  }
}

/// \brief special disjunctive constraint with three implicit hard constraints \f$x \leq xinfty\f$ and \f$y \leq yinfty\f$ and \f$x < xinfty \wedge y < yinfty \Rightarrow (x \geq y + csty \vee y \geq x + cstx)\f$ and an additional cost function \f$((x = xinfty)?costx:0) + ((y= yinfty)?costy:0)\f$
int WCSP::postSpecialDisjunction(int xIndex, int yIndex, Value cstx, Value csty, Value xinfty, Value yinfty, Cost costx, Cost costy) {
  if (!vars[xIndex]->enumerated() && !vars[yIndex]->enumerated()) {
	SpecialDisjunction
			*ctr =
					new SpecialDisjunction(this, (IntervalVariable *) vars[xIndex], (IntervalVariable *) vars[yIndex], cstx, csty, xinfty, yinfty, costx, costy, &storeData->storeCost, &storeData->storeValue);
	return ctr->wcspIndex;
  } else if (vars[xIndex]->enumerated() && vars[yIndex]->enumerated()) {
	EnumeratedVariable *x = (EnumeratedVariable *) vars[xIndex];
	EnumeratedVariable *y = (EnumeratedVariable *) vars[yIndex];
	vector<Cost> costs;
	for (unsigned int a = 0; a < x->getDomainInitSize(); a++) {
	  for (unsigned int b = 0; b < y->getDomainInitSize(); b++) {
		costs.push_back((x->toValue(a) <= xinfty && y->toValue(b) <= yinfty && (x->toValue(a) == xinfty || y->toValue(b) == yinfty || (x->toValue(a) >= y->toValue(b) + csty || y->toValue(b) >= x->toValue(a) + cstx)))?(((x->toValue(a) == xinfty)?costx:MIN_COST) + ((y->toValue(b) == yinfty)?costy:MIN_COST)):getUb());
	  }
	}
	return postBinaryConstraint(xIndex, yIndex, costs);
  } else {
	cerr << "Cannot mix variables with interval and enumerated domains!!!" << endl;
	exit(EXIT_FAILURE);
  }
}

void WCSP::sortConstraints()
{
  for (vector<int>::iterator idctr = delayedNaryCtr.begin(); idctr != delayedNaryCtr.end(); ++idctr) {
	BinaryConstraint* bctr;
	TernaryConstraint* tctr = new TernaryConstraint(this, &storeData->storeCost);
	elimTernConstrs.push_back(tctr);
	for (int j = 0; j < 3; j++) {
		if (!ToulBar2::vac) bctr = new BinaryConstraint(this, &storeData->storeCost);
		else bctr = new VACBinaryConstraint(this, &storeData->storeCost);
		elimBinConstrs.push_back(bctr);
	}
	getCtr(*idctr)->propagate();
  }
  delayedNaryCtr.clear();
  isDelayedNaryCtr = false;
  if(ToulBar2::Berge_Dec > 0 ) {
	// flag pour indiquer si une variable a deja ete visitee initialement a faux
	vector<bool> marked(numberOfVariables(), false);
	vector<int> revdac;
	// nouvel ordre DAC inverse
	//	for (int i = numberOfVariables()-1 ; i >= 0; i--) { if (!marked[i]){ visit(i,revdac,marked,listofsuccessors); }}
	//	for (unsigned int i = 0; i <  numberOfVariables(); i++) { if (!marked[i]){ visit(i,revdac,marked,listofsuccessors); }}

	//Mark native variable
	for (int i = ToulBar2::nbvar-1; i >= 0; i--) { if (!marked[i]){ visit(i,revdac,marked,listofsuccessors); }}
	//Mark q variable only 
	for (int i = numberOfVariables()-1 ; i > ToulBar2::nbvar-1 ; i--) { if (!marked[i]){ visit(i,revdac,marked,listofsuccessors); }}

	// listofsuccessors.clear(); // appel a la methode clear de l'objet vector

	if( ToulBar2::verbose >= 1 ) { 
	  cout << "BERGE DAC reverse order:";
	  for (unsigned int i = 0; i < numberOfVariables(); i++) {
		cout << " " << revdac[i];
	  }
	  cout << endl;
	}

	assert( revdac.size() == numberOfVariables() );

	setDACOrder(revdac);
  }
  if (ToulBar2::varOrder) {
	vector<int> order;
	elimOrderFile2Vector(ToulBar2::varOrder, order);
	setDACOrder(order);
  }
  for (unsigned int i = 0; i < vars.size(); i++) {
	vars[i]->sortConstraints();
  }
}

void WCSP::updateCurrentVarsId() {
	int pos = 0;
	for (unsigned int i = 0; i < vars.size(); i++) {
		if (vars[i]->unassigned()) {
			vars[i]->setCurrentVarId(pos);
			pos++;
		}
	}
}

bool cmpTernaryConstraint(TernaryConstraint *c1, TernaryConstraint *c2)
{
    int v1 = c1->getVar(c1->getDACScopeIndex())->getDACOrder();
    int v2 = c2->getVar(c2->getDACScopeIndex())->getDACOrder();
    if (v1 < v2) return true;
    else if (v1==v2) {
    	v1 = min(c1->getVar((c1->getDACScopeIndex()+1)%3)->getDACOrder(),c1->getVar((c1->getDACScopeIndex()+2)%3)->getDACOrder());
        v2 = min(c2->getVar((c2->getDACScopeIndex()+1)%3)->getDACOrder(),c2->getVar((c2->getDACScopeIndex()+2)%3)->getDACOrder());
        if (v1 < v2) return true;
        else if (v1==v2) {
        	v1 = max(c1->getVar((c1->getDACScopeIndex()+1)%3)->getDACOrder(),c1->getVar((c1->getDACScopeIndex()+2)%3)->getDACOrder());
            v2 = max(c2->getVar((c2->getDACScopeIndex()+1)%3)->getDACOrder(),c2->getVar((c2->getDACScopeIndex()+2)%3)->getDACOrder());
        	return (v1 < v2);
        }
    }
    return false;
}

void WCSP::processTernary() {
	// double maxtight = 0;
	// Variable* var;
	// TernaryConstraint *tctr1max = NULL, *tctr2max = NULL;

	// if (ToulBar2::preprocessTernaryHeuristic) {

	//     for (unsigned int i=0; i<vars.size(); i++) {
	// 		TernaryConstraint *tctr1, *tctr2;
	//         double tight = vars[i]->strongLinkedby(var,tctr1,tctr2);
	//         if(tight > maxtight) { maxtight = tight; tctr1max = tctr1; tctr2max = tctr2; }
	// 		if(tctr1 && tctr2 && (tctr1 != tctr2)) {
	// 			tctr1->extendTernary();
	// 			tctr2->extendTernary();
	// 			BinaryConstraint* b = tctr1->commonBinary(tctr2);
	// 			if(b->connected())
	// 			{
	//     			tctr1->projectTernaryBinary(b);
	//     			tctr2->projectTernaryBinary(b);
	//     			b->propagate();
	// 			}
	// 		}
	//     }

	//     if(ToulBar2::verbose > 0) {
	//         cout << "Strongest part has mean cost: " << maxtight;
	//         if(var) cout << "  Variable: " << var->wcspIndex;
	//         if(tctr1max) cout << ", 1. ternary with tight: " << tctr1max->getTightness();
	//         if(tctr2max) cout << ", 2. ternary with tight: " << tctr2max->getTightness();
	//         cout << endl;
	//     }
	// }

	vector<TernaryConstraint*> ternaries;
	for (unsigned int i = 0; i < constrs.size(); i++)
		if (constrs[i]->connected() && !constrs[i]->isSep() && constrs[i]->extension() && constrs[i]->arity() == 3) {
			TernaryConstraint* t = (TernaryConstraint*) constrs[i];
			ternaries.push_back(t);
		}
	for (int i = 0; i < elimTernOrder; i++)
		if (elimTernConstrs[i]->connected()) {
			TernaryConstraint* t = (TernaryConstraint*) elimTernConstrs[i];
			ternaries.push_back(t);
		}
    sort(ternaries.begin(), ternaries.end(), cmpTernaryConstraint);
	for (int i = ternaries.size()-1; i >= 0; i--) {
		TernaryConstraint* t = ternaries[i];
//		cout << "PROJECT&SUBTRACT tern(" << t->getVar(0)->getName() << "," << t->getVar(1)->getName() << "," << t->getVar(2)->getName() << ")" << endl;
		t->extendTernary();
		if (ToulBar2::costfuncSeparate) t->decompose();
		if (t->connected()) t->projectTernary();
	}
}

/// \defgroup preprocessing Preprocessing techniques
/// Depending on toulbar2 options, the sequence of preprocessing techniques applied before the search is:
/// -# \e i-bounded variable elimination with user-defined \e i bound
/// -# pairwise decomposition of cost functions (binary cost functions are implicitly decomposed by soft AC and empty cost function removals)
/// -# MinSumDiffusion propagation (see VAC)
/// -# projects\&substracts n-ary cost functions in extension on all the binary cost functions inside their scope (3 < n < max, see toulbar2 options)
/// -# functional variable elimination (see \ref varelim)
/// -# projects\&substracts ternary cost functions in extension on their three binary cost functions inside their scope (before that, extends the existing binary cost functions to the ternary cost function and applies pairwise decomposition)
/// -# creates new ternary cost functions for all triangles (\e ie occurences of three binary cost functions \e xy, \e yz, \e zx)
/// -# removes empty cost functions while repeating #1 and #2 until no new cost functions can be removed
///
/// \note the propagation loop is called after each preprocessing technique (see \ref WCSP::propagate)

void WCSP::preprocessing() {
	Eliminate.clear();
	if (ToulBar2::elimDegree_preprocessing <= -3) {
	  int deg = medianDegree();
	  Double domsize = medianDomainSize();
	  Double size = (Double) numberOfUnassignedVariables() * (sizeof(Char)*deg+sizeof(Cost)) * Pow(domsize,deg+1);
	  if (ToulBar2::debug>=2) cout << "MAX ESTIMATED ELIM SIZE: " << size << endl;
	  if (deg >= 3 && deg <= -ToulBar2::elimDegree_preprocessing && size < (Double) ToulBar2::elimSpaceMaxMB * 1024. * 1024.) {
		ToulBar2::elimDegree_preprocessing = deg;
		cout << "Generic variable elimination of degree " << deg << endl;
	  } else {
		ToulBar2::elimDegree_preprocessing = -1;
		cout << "Generic variable elimination disabled." << endl;
	  }
	}
	if (ToulBar2::elimDegree >= 0 || ToulBar2::elimDegree_preprocessing >= 0 || ToulBar2::preprocessFunctional > 0) {
		initElimConstrs();
		if (ToulBar2::elimDegree_preprocessing >= 0) {
			if (ToulBar2::verbose >= 1) cout << "Variable elimination in preprocessing of true degree <= "
					<< ToulBar2::elimDegree_preprocessing << endl;
			ToulBar2::elimDegree_preprocessing_ = ToulBar2::elimDegree_preprocessing;
			maxDegree = -1;
			propagate();
			cout << "Maximum degree of generic variable elimination: " << maxDegree << endl;
		} else if (ToulBar2::elimDegree >= 0) {
			ToulBar2::elimDegree_ = ToulBar2::elimDegree;
		}
	}

	int posConstrs = 0;
	int posElimTernConstrs = 0;
	if (ToulBar2::costfuncSeparate) {
		double time = cpuTime();
		for (unsigned int i = 0; i < constrs.size(); i++) {
			if (constrs[i]->connected() && !constrs[i]->isSep()) {
				constrs[i]->decompose();
			}
		}
		for (int i = 0; i < elimTernOrder; i++) {
			if (elimTernConstrs[i]->connected() && !elimTernConstrs[i]->isSep()) {
				elimTernConstrs[i]->decompose();
			}
		}
		posConstrs = constrs.size();
		posElimTernConstrs = elimTernOrder;
		if (ToulBar2::verbose >= 0) cout << "Cost function decomposition time : " << cpuTime() - time << " seconds.\n";
	}

	propagate();

	// recompute current DAC order and its reverse
	vector<int> elimorder(numberOfVariables(), -1);
	vector<int> revelimorder(numberOfVariables(), -1);
	for (unsigned int i = 0; i < numberOfVariables(); i++) {
	  revelimorder[getVar(i)->getDACOrder()] = i;
	  elimorder[numberOfVariables() - getVar(i)->getDACOrder() -1] = i;
	}	
//	cout << "DAC:";
//	for (int i = 0; i < numberOfVariables(); i++) {
//		cout << " " << revelimorder[i];
//	}
//	cout << endl;
//	cout << "REVDAC:";
//	for (int i = 0; i < numberOfVariables(); i++) {
//		cout << " " << elimorder[i];
//	}
//	cout << endl;
	setDACOrder(revelimorder);
	setDACOrder(elimorder);

	if (ToulBar2::preprocessNary > 0) {
		for (unsigned int i = 0; i < constrs.size(); i++) {
			if (constrs[i]->connected() && !constrs[i]->isSep() && (constrs[i]->arity() > 3) && (constrs[i]->arity()
					<= ToulBar2::preprocessNary) && constrs[i]->extension()) {
				NaryConstraintMap* nary = (NaryConstraintMap*) constrs[i];
				if (nary->size() >= 2 || nary->getDefCost() > MIN_COST) {
					nary->keepAllowedTuples(getUb());
					nary->preprojectall2();
					//			if (nary->connected()) nary->preproject3();
				}
			}
		}
//		processTernary();
		propagate();
	}

	// Merge functional (binary bijection only) variables in decision variables
	bool merged = (ToulBar2::preprocessFunctional > 0);
	while (merged) {
		merged = false;
		for (unsigned int i = 0; i < constrs.size(); i++)
			if (constrs[i]->connected() && !constrs[i]->isSep() && constrs[i]->arity() == 2 && constrs[i]->extension()
					&& (ToulBar2::preprocessFunctional == 1 || constrs[i]->getVar(0)->getDomainSize()
							== constrs[i]->getVar(1)->getDomainSize())) {
				BinaryConstraint *xy = (BinaryConstraint *) constrs[i];
				EnumeratedVariable *x = (EnumeratedVariable*) xy->getVar(0);
				EnumeratedVariable *y = (EnumeratedVariable*) xy->getVar(1);
				map<Value, Value> functional;
				if (xy->isFunctional(x, y, functional) && y->canbeMerged()) {
					y->mergeTo(xy, functional);
					merged = true;
					propagate();
				} else if (xy->isFunctional(y, x, functional) && x->canbeMerged()) {
					x->mergeTo(xy, functional);
					merged = true;
					propagate();
				}
			}
		for (int i = 0; i < elimBinOrder; i++)
			if (elimBinConstrs[i]->connected() && !elimBinConstrs[i]->isSep()
					&& elimBinConstrs[i]->getVar(0)->getDomainSize() == elimBinConstrs[i]->getVar(1)->getDomainSize()) {
				assert(elimBinConstrs[i]->arity()==2);
				assert(elimBinConstrs[i]->extension());
				BinaryConstraint *xy = (BinaryConstraint *) elimBinConstrs[i];
				EnumeratedVariable *x = (EnumeratedVariable*) xy->getVar(0);
				EnumeratedVariable *y = (EnumeratedVariable*) xy->getVar(1);
				map<Value, Value> functional;
				if (xy->isFunctional(x, y, functional) && y->canbeMerged()) {
					y->mergeTo(xy, functional);
					merged = true;
					propagate();
				} else if (xy->isFunctional(y, x, functional) && x->canbeMerged()) {
					x->mergeTo(xy, functional);
					merged = true;
					propagate();
				}
			}
	}


	if (ToulBar2::preprocessTernaryRPC) {
		ternaryCompletion();
		setDACOrder(revelimorder);
		processTernary();
		propagate();
		setDACOrder(elimorder);
		processTernary();
		propagate();
	} else if (ToulBar2::preprocessNary > 0) {
		processTernary();
		propagate();
	}

	if (ToulBar2::minsumDiffusion && ToulBar2::vac) vac->minsumDiffusion();
	if (ToulBar2::vac) {
		if (ToulBar2::verbose >= 0) cout << "Preprocessing ";
		vac->printStat(true);
		//    	vac->afterPreprocessing();
		for (unsigned int i = 0; i < vars.size(); i++)
			vars[i]->queueEliminate();
		propagate();
	}

	// Deconnect empty cost functions
	unsigned int nbunvar;
	do {
		nbunvar = numberOfUnassignedVariables();
		for (unsigned int i = 0; i < constrs.size(); i++)
			if (constrs[i]->connected() && !constrs[i]->isSep() && constrs[i]->universal()) {
				if (ToulBar2::verbose >= 3) cout << "deconnect empty cost function: " << *constrs[i];
				constrs[i]->deconnect(true);
			}
		for (int i = 0; i < elimBinOrder; i++)
			if (elimBinConstrs[i]->connected() && !elimBinConstrs[i]->isSep() && elimBinConstrs[i]->universal()) {
				if (ToulBar2::verbose >= 3) cout << "deconnect empty cost function: " << *elimBinConstrs[i];
				elimBinConstrs[i]->deconnect(true);
			}
		for (int i = 0; i < elimTernOrder; i++)
			if (elimTernConstrs[i]->connected() && !elimTernConstrs[i]->isSep() && elimTernConstrs[i]->universal()) {
				if (ToulBar2::verbose >= 3) cout << "deconnect empty cost function: " << *elimTernConstrs[i];
				elimTernConstrs[i]->deconnect(true);
			}
		if (ToulBar2::costfuncSeparate) {
			for (unsigned int i = posConstrs; i < constrs.size(); i++) {
				if (constrs[i]->connected() && !constrs[i]->isSep()) {
					constrs[i]->decompose();
				}
			}
			for (int i = posElimTernConstrs; i < elimTernOrder; i++) {
				if (elimTernConstrs[i]->connected() && !elimTernConstrs[i]->isSep()) {
					elimTernConstrs[i]->decompose();
				}
			}
			posConstrs = constrs.size();
			posElimTernConstrs = elimTernOrder;
		}
		propagate();
	} while (numberOfUnassignedVariables() < nbunvar);

	if (ToulBar2::elimDegree >= 0 || ToulBar2::elimDegree_preprocessing >= 0 || ToulBar2::preprocessFunctional > 0) {
		ToulBar2::elimDegree_preprocessing_ = -1;
		if (ToulBar2::elimDegree >= 0) {
			ToulBar2::elimDegree_ = ToulBar2::elimDegree;
			if (ToulBar2::elimDegree_preprocessing < min(2,ToulBar2::elimDegree)) {
				for(int i=numberOfVariables()-1; i>=0; i--) vars[i]->queueEliminate();
				propagate();
			}
			if (ToulBar2::verbose >= 1) cout << "Variable elimination during search of degree <= "
					<< ToulBar2::elimDegree << endl;
		}
	}

#ifdef BOOST
	if (ToulBar2::MSTDAC) spanningTreeOrdering();
#endif
}

Value WCSP::getDomainSizeSum() {
	//    cout << " " << connectedComponents() << endl;
	Value sum = 0;
	for (unsigned int i = 0; i < vars.size(); i++) {
		if (vars[i]->unassigned()) sum += vars[i]->getDomainSize();
	}
	return sum;
}

bool WCSP::getEnumDomain(int varIndex, Value *array) {
	if (EnumeratedVariable *var = dynamic_cast<EnumeratedVariable*>(vars[varIndex])) {
		var->getDomain(array);
		return true;
	} else return false;
}

bool WCSP::getEnumDomainAndCost(int varIndex, ValueCost *array) {
	if (EnumeratedVariable *var = dynamic_cast<EnumeratedVariable*>(vars[varIndex])) {
		var->getDomainAndCost(array);
		return true;
	} else return false;
}

unsigned int WCSP::numberOfConnectedConstraints() const {
	int res = 0;
	for (unsigned int i = 0; i < constrs.size(); i++)
		if (constrs[i]->connected() && !constrs[i]->isSep()) res++;
	for (int i = 0; i < elimBinOrder; i++)
		if (elimBinConstrs[i]->connected() && !elimBinConstrs[i]->isSep()) res++;
	for (int i = 0; i < elimTernOrder; i++)
		if (elimTernConstrs[i]->connected() && !elimTernConstrs[i]->isSep()) res++;
	return res;
}

unsigned int WCSP::numberOfConnectedBinaryConstraints() const {
	int res = 0;
	for (unsigned int i = 0; i < constrs.size(); i++)
		if (constrs[i]->connected() && constrs[i]->arity()==2 && !constrs[i]->isSep()) res++;
	for (int i = 0; i < elimBinOrder; i++)
		if (elimBinConstrs[i]->connected() && !elimBinConstrs[i]->isSep()) res++;
	return res;  
}

unsigned int WCSP::medianDomainSize() const {
  unsigned int nbunvars = numberOfUnassignedVariables();
  if (nbunvars==0) return 0;
  int domain[nbunvars];
  unsigned int pos = 0;
  for (unsigned int i=0; i<vars.size(); i++) if (unassigned(i)) {domain[pos]=getDomainSize(i); pos++;}
  assert(pos == numberOfUnassignedVariables() && pos == nbunvars);
  return stochastic_selection<int>(domain, 0, nbunvars-1, nbunvars / 2);
}

unsigned int WCSP::medianDegree() const {
  unsigned int nbunvars = numberOfUnassignedVariables();
  if (nbunvars==0) return 0;
  int degree[nbunvars];
  unsigned int pos = 0;
  for (unsigned int i=0; i<vars.size(); i++) if (unassigned(i)) {degree[pos]=getTrueDegree(i); pos++;}
  assert(pos == numberOfUnassignedVariables() && pos == nbunvars);
  return stochastic_selection<int>(degree, 0, nbunvars-1, nbunvars / 2);
}

void WCSP::printNCBuckets() {
	int lastbucket = -1;
	for (int bucket = 0; bucket < NCBucketSize; bucket++) {
		if (NCBuckets[bucket].begin() != NCBuckets[bucket].end()) lastbucket = bucket;
	}

	for (int bucket = 0; bucket <= lastbucket; bucket++) {
		cout << "NC " << bucket << ":";
		for (VariableList::iterator iter = NCBuckets[bucket].begin(); iter != NCBuckets[bucket].end(); ++iter) {
			cout << " " << (*iter)->getName() << "," << (*iter)->getMaxCostValue() << "," << (*iter)->getMaxCost();

			assert((*iter)->canbe((*iter)->getMaxCostValue()));
			assert((*iter)->getCost((*iter)->getMaxCostValue()) == (*iter)->getMaxCost() || !LUBTEST((*iter)->getMaxCost(),(*iter)->getCost((*iter)->getMaxCostValue())));
			assert((bucket && !PARTIALORDER) ? (to_double((*iter)->getMaxCost()) >= (Long) pow(2.,bucket)) : ((*iter)->getMaxCost() > MIN_COST));
			assert(PARTIALORDER || to_double((*iter)->getMaxCost()) < (Long) pow(2.,bucket+1));
		}
		cout << endl;
	}
}

/** \defgroup verbosity Output messages and verbosity options
 *
 * Depending on verbosity level given as option "-v=level", \p toulbar2 will output:
 * - (level=0, no verbosity) default output mode: shows version number, number of variables and cost functions read in the problem file,
 *  number of unassigned variables and cost functions after preprocessing,
 *  problem upper and lower bounds after preprocessing.
 *  Outputs current best solution cost found,
 *  ends by giving the optimum or "No solution".
 *  Last output line should always be: "end."
 * - (level=-1, no verbosity) restricted output mode: do not print current best solution cost found
 * -# (level=1) shows also search choices ("["\e search_depth \e problem_lower_bound \e problem_upper_bound \e sum_of_current_domain_sizes"] Try" \e variable_index \e operator \e value)
 *  with \e operator being assignment ("=="), value removal ("!="), domain splitting ("<=" or ">=", also showing EAC value in parenthesis)
 * -# (level=2) shows also current domains (\e variable_index \e list_of_current_domain_values "/" \e number_of_cost_functions (see approximate degree in \ref varelim) "/" \e weighted_degree \e list_of_unary_costs "s:" \e support_value) before each search choice
 *  and reports problem lower bound increases, NC bucket sort data (see \ref ncbucket), and basic operations on domains of variables
 * -# (level=3) reports also basic arc EPT operations on cost functions (see \ref softac)
 * -# (level=4) shows also current list of cost functions for each variable and reports more details on arc EPT operations (showing all changes in cost functions)
 * -# (level=5) reports more details on cost functions defined in extension giving their content (cost table by first increasing values in the current domain of the last variable in the scope)
 *
 * For debugging purposes, another option "-Z=level" allows to monitor the search:
 * -# (level 1) shows current search depth (number of search choices from the root of the search tree) and reports statistics on nogoods for BTD-like methods
 * -# (level 2) idem
 * -# (level 3) also saves current problem into a file before each search choice
 *
 * \note \p toulbar2, compiled in debug mode, can be more verbose and it checks a lot of assertions (pre/post conditions in the code)
 *
 * \note \p toulbar2 will output an help message giving available options if run without any parameters
 *
 */

void WCSP::print(ostream& os) {
	os << "Objective: [" << getLb() << "," << getUb() << "]" << endl;
	os << "Variables:" << endl;
	for (unsigned int i = 0; i < vars.size(); i++)
		os << *vars[i] << endl;
	if (ToulBar2::verbose >= 4) {
		os << "Constraints:" << endl;
		for (unsigned int i = 0; i < constrs.size(); i++)
			if (constrs[i]->connected()) os << *constrs[i];
		for (int i = 0; i < elimBinOrder; i++)
			if (elimBinConstrs[i]->connected()) os << *elimBinConstrs[i];
		for (int i = 0; i < elimTernOrder; i++)
			if (elimTernConstrs[i]->connected()) os << *elimTernConstrs[i];
	}
}

void printClique(ostream& os, int arity, Constraint *ctr) {
	assert(arity >= 2);
	for (int i = 0; i < arity - 1; i++) {
		for (int j = i + 1; j < arity; j++) {
			os << "v" << ctr->getVar(i)->wcspIndex + 1 << " -- v" << ctr->getVar(j)->wcspIndex + 1 << " [len= " << ctr->getTightness() << "];" << endl;
		}
	}
}

// Warning! make the assumption that all initial domains start at zero!!!
void WCSP::dump(ostream& os, bool original) {
	Value maxdomsize = 0;
	unsigned int maxdomsizeUI = 0;
	Value xcosts = 0;
	// dump filename
	char Pb_basename[80];
	char Pb_graph[80];
	char Pb_degree[80];

	strcpy(Pb_basename, "problem");
	strcpy(Pb_graph, Pb_basename);
	strcpy(Pb_degree, Pb_basename);

	if (getLb() > MIN_COST) xcosts++;
	for (unsigned int i = 0; i < vars.size(); i++) {
		if (vars[i]->getInf() < 0) {
			cerr << "Cannot save domain of variable " << vars[i]->getName() << " with negative values!!!" << endl;
			exit(EXIT_FAILURE);
		}
		if (original) {
			if (vars[i]->getSup() + 1 > maxdomsize) maxdomsize = vars[i]->getSup() + 1;
		} else {
			if (vars[i]->unassigned() && vars[i]->getDomainSize() > maxdomsizeUI) maxdomsizeUI
					= vars[i]->getDomainSize();
		}
		if (vars[i]->enumerated() && (original || vars[i]->unassigned())) xcosts++;
		//          else if (vars[i]->getInfCost() > MIN_COST || vars[i]->getSupCost() > MIN_COST) {
		//              cerr << "Cannot save interval variable " << vars[i]->getName() << " with bound unary costs!!!" << endl;
		//              exit(EXIT_FAILURE);
		//          }
	}
	os << "wcsp " << ((original) ? numberOfVariables() : numberOfUnassignedVariables()) << " "
			<< ((original) ? maxdomsize : maxdomsizeUI) << " " << numberOfConnectedConstraints() + xcosts << " "
			<< getUb() << endl;
	unsigned int nbvar = 0;
	for (unsigned int i = 0; i < vars.size(); i++) {
		if (original) {
			if (!vars[i]->enumerated()) os << "-";
			os << vars[i]->getSup() + 1;
			if (i < vars.size() - 1) os << " ";
		} else if (vars[i]->unassigned()) {
			nbvar++;
			if (!vars[i]->enumerated()) os << "-";
			os << vars[i]->getDomainSize();
			if (nbvar < numberOfUnassignedVariables()) os << " ";
		}
	}
	if (((original) ? numberOfVariables() : numberOfUnassignedVariables()) > 0) os << endl;
	for (unsigned int i = 0; i < constrs.size(); i++)
		if (constrs[i]->connected() && !constrs[i]->isSep()) constrs[i]->dump(os, original);
	for (int i = 0; i < elimBinOrder; i++)
		if (elimBinConstrs[i]->connected() && !elimBinConstrs[i]->isSep()) elimBinConstrs[i]->dump(os, original);
	for (int i = 0; i < elimTernOrder; i++)
		if (elimTernConstrs[i]->connected() && !elimTernConstrs[i]->isSep()) elimTernConstrs[i]->dump(os, original);
	for (unsigned int i = 0; i < vars.size(); i++) {
		if (vars[i]->enumerated() && (original || vars[i]->unassigned())) {
			int size = vars[i]->getDomainSize();
			ValueCost domcost[size]; // replace size by MAX_DOMAIN_SIZE in case of compilation problem
			getEnumDomainAndCost(i, domcost);
			os << "1 " << ((original) ? i : vars[i]->getCurrentVarId()) << " " << getUb() << " " << size << endl;
			for (int v = 0; v < size; v++) {
				os << ((original) ? (domcost[v].value) : v) << " "
						<< ((original) ? domcost[v].cost : min(getUb(), domcost[v].cost)) << endl;
			}
		}
	}
	if (getLb() > MIN_COST) os << "0 " << getLb() << " 0" << endl;

	//####################" dump dot file ###############################""
	cout << " Graph structure saved in problem.dot " << endl;
	if (original) {
		strcat(Pb_graph, "_original.dot");
	} else {
		strcat(Pb_graph, ".dot");
	}
	//ofstream pb("current_problem.graph");
	ofstream pb(Pb_graph);
	pb << " graph graphname {\n " << endl;
	int res = 0;
	for (unsigned int i = 0; i < constrs.size(); i++)
		if (constrs[i]->connected()) res += (constrs[i]->arity() * (constrs[i]->arity() - 1) / 2);
	for (int i = 0; i < elimBinOrder; i++)
		if (elimBinConstrs[i]->connected()) res += (elimBinConstrs[i]->arity() * (elimBinConstrs[i]->arity() - 1) / 2);
	for (int i = 0; i < elimTernOrder; i++)
		if (elimTernConstrs[i]->connected()) res += (elimTernConstrs[i]->arity() * (elimTernConstrs[i]->arity() - 1)
				/ 2);
	pb << "// number of constraint = " << res << " number of variable=  " << numberOfVariables() << endl;
	for (unsigned int i = 0; i < constrs.size(); i++)
		if (constrs[i]->connected()) {
			//            pb << constrs[i]->getVar(0)->wcspIndex + 1;
			//            for (int j=1; j<constrs[i]->arity(); j++) {
			//                pb << " " << constrs[i]->getVar(j)->wcspIndex + 1;
			//            }
			//            pb << endl;
			printClique(pb, constrs[i]->arity(), constrs[i]);
		}
	for (int i = 0; i < elimBinOrder; i++)
		if (elimBinConstrs[i]->connected()) {
			//            pb << elimBinConstrs[i]->getVar(0)->wcspIndex + 1;
			//            for (int j=1; j<elimBinConstrs[i]->arity(); j++) {
			//                pb << " " << elimBinConstrs[i]->getVar(j)->wcspIndex + 1;
			//            }
			//            pb << endl;
			printClique(pb, elimBinConstrs[i]->arity(), elimBinConstrs[i]);
		}
	for (int i = 0; i < elimTernOrder; i++)
		if (elimTernConstrs[i]->connected()) {
			//            pb << elimTernConstrs[i]->getVar(0)->wcspIndex + 1;
			//            for (int j=1; j<elimTernConstrs[i]->arity(); j++) {
			//                pb << " " << elimTernConstrs[i]->getVar(j)->wcspIndex + 1;
			//            }
			//            pb << endl;
			printClique(pb, elimTernConstrs[i]->arity(), elimTernConstrs[i]);
		}
	pb << "}" << endl;

	//####################" end dump dot file ###############################""
	//#######################dump degree distribution ###################
#ifdef BOOST
	cout << "Connected components: " << connectedComponents() << endl;
    cout << "Biconnected components: " << biConnectedComponents() << endl;
	cout << "Diameter : " << diameter() << endl;
#endif

	int* degDistrib = new int[vars.size()];

	for (unsigned int i = 0; i < vars.size(); i++)
		degDistrib[i] = 0;
	for (unsigned int i = 0; i < vars.size(); i++)
		if (unassigned(i)) degDistrib[getTrueDegree(i)]++;

	unsigned int lastnonzero = 0;
	for (unsigned int i = 0; i < vars.size(); i++)
		if (degDistrib[i]) lastnonzero = i;

	if (original) {
		strcat(Pb_degree, "_original.degree"); // original distribution
	} else {
		strcat(Pb_degree, ".degree"); // after preprocessing
	}
	ofstream file(Pb_degree);
	for (unsigned int i = 0; i <= lastnonzero; i++)
		if (degDistrib[i]) file << i << " " << degDistrib[i] << endl;
	delete[] degDistrib;

	//#######################dump degree distribution ###################


	if (ToulBar2::pedigree) ToulBar2::pedigree->save("problem.pre", this, false, true);
}

ostream& operator<<(ostream& os, WCSP &wcsp) {
	wcsp.print(os);
	return os;
}

ostream& operator<<(ostream& os, WeightedCSP &wcsp) {
	wcsp.print(os);
	return os;
}

/*
 * WCSP propagation methods
 *
 */

bool WCSP::verify() {
	for (unsigned int i = 0; i < vars.size(); i++) {
		if (vars[i]->unassigned()) {
			if (td) {
				if (td->isActiveAndInCurrentClusterSubTree(vars[i]->getCluster())) {
					if (!vars[i]->verifyNC()) return false;
#ifdef DEECOMPLETE
					if (ToulBar2::DEE && !vars[i]->verifyDEE()) return false;
#endif
				}
			} else {
				if (!vars[i]->verifyNC()) return false;
#ifdef DEECOMPLETE
				if (ToulBar2::DEE && !vars[i]->verifyDEE()) return false;
#endif
			}
		}
		// Warning! in the CSP case, EDAC is no equivalent to GAC on ternary constraints due to the combination with binary constraints
		//isEAC may change current support for variables and constraints
		if (ToulBar2::LcLevel == LC_EDAC && vars[i]->unassigned() && !CSP(getLb(), getUb()) && !vars[i]->isEAC()) {
			cout << "support of variable " << vars[i]->getName() << " not EAC!" << endl;
			return false;
		}
	}
	if (ToulBar2::LcLevel >= LC_AC) {
		for (unsigned int i = 0; i < constrs.size(); i++) {
			if (constrs[i]->connected() && !constrs[i]->verify()) return false;
		}
		for (int i = 0; i < elimBinOrder; i++) {
			if (elimBinConstrs[i]->connected() && !elimBinConstrs[i]->verify()) return false;
		}
		for (int i = 0; i < elimTernOrder; i++) {
			if (elimTernConstrs[i]->connected() && !elimTernConstrs[i]->verify()) return false;
		}
	}
	return true;
}

void WCSP::whenContradiction() {
	NC.clear();
	IncDec.clear();
	AC.clear();
	DAC.clear();
	EAC1.clear();
	EAC2.clear();
	Eliminate.clear();
	DEE.clear();
	objectiveChanged = false;
	nbNodes++;
}

///\defgroup ncbucket NC bucket sort
/// maintains a sorted list of variables having non-zero unary costs in order to make NC propagation incremental.\n
/// - variables are sorted into buckets
/// - each bucket is associated to a single interval of non-zero costs (using a power-of-two scaling, first bucket interval is [1,2[, second interval is [2,4[, etc.)
/// - each variable is inserted into the bucket corresponding to its largest unary cost in its domain
/// - variables having all unary costs equal to zero do not belong to any bucket
///
/// NC propagation will revise only variables in the buckets associated to costs sufficiently large wrt current objective bounds.

void WCSP::propagateNC() {
	if (ToulBar2::verbose >= 2) cout << "NCQueue size: " << NC.getSize() << endl;
	while (!NC.empty()) {
		Variable *x = NC.pop();
		if (x->unassigned()) x->propagateNC();
	}
	if (ToulBar2::verbose >= 3) {
		for (unsigned int i = 0; i < vars.size(); i++)
			cout << *vars[i] << endl;
	}
	if (ToulBar2::verbose >= 2) printNCBuckets();

	if (objectiveChanged) {
		objectiveChanged = false;
		int bucket = cost2log2glb(getUb() - getLb());
		if (bucket < 0) bucket = 0;
		for (; bucket < NCBucketSize; bucket++) {
			for (VariableList::iterator iter = NCBuckets[bucket].begin(); iter != NCBuckets[bucket].end();) {
				Variable *x = *iter;
				++iter; // Warning! the iterator could be moved to another place by propagateNC
				if (x->unassigned() && CUT(x->getMaxCost() + getLb(), getUb())) {
					if (td) {
						if (td->isActiveAndInCurrentClusterSubTree(x->getCluster())) x->propagateNC();
					} else x->propagateNC();
				}
			}
		}
	}
	if (objectiveChanged || !NC.empty()) propagateNC();
}

void WCSP::propagateIncDec() {
	if (ToulBar2::verbose >= 2) cout << "IncDecQueue size: " << IncDec.getSize() << endl;
	while (!IncDec.empty()) {
		int incdec;
		Variable *x = IncDec.pop(&incdec);
		if (x->unassigned()) x->propagateIncDec(incdec);
	}
}

void WCSP::propagateAC() {
	if (ToulBar2::verbose >= 2) cout << "ACQueue size: " << AC.getSize() << endl;
	while (!AC.empty()) {
		EnumeratedVariable *x = (EnumeratedVariable *) ((ToulBar2::QueueComplexity) ? AC.pop_min() : AC.pop());
		if (x->unassigned()) x->propagateAC();
		// Warning! propagateIncDec() necessary to transform inc/dec event into remove event
		propagateIncDec(); // always examine inc/dec events before remove events
	}
}

void WCSP::propagateDAC() {
	if (ToulBar2::verbose >= 2) cout << "DACQueue size: " << DAC.getSize() << endl;
	while (!DAC.empty()) {
		EnumeratedVariable *x = (EnumeratedVariable *) ((ToulBar2::QueueComplexity) ? DAC.pop_max() : DAC.pop());
		if (x->unassigned()) x->propagateDAC();
		propagateIncDec(); // always examine inc/dec events before projectFromZero events
	}
}

void WCSP::fillEAC2() {
	assert(EAC2.empty());
	if (ToulBar2::verbose >= 2) cout << "EAC1Queue size: " << EAC1.getSize() << endl;
	while (!EAC1.empty()) {
		EnumeratedVariable *x = (EnumeratedVariable *) ((ToulBar2::QueueComplexity) ? EAC1.pop_min() : EAC1.pop());
		if (x->unassigned()) x->fillEAC2(true);
	}
}

void WCSP::propagateEAC() {
	fillEAC2();
	if (ToulBar2::verbose >= 2) cout << "EAC2Queue size: " << EAC2.getSize() << endl;
	while (!EAC2.empty()) {
		EnumeratedVariable *x = (EnumeratedVariable *) ((ToulBar2::QueueComplexity) ? EAC2.pop_min() : EAC2.pop());
		if (x->unassigned()) x->propagateEAC();
		propagateIncDec(); // always examine inc/dec events before projectFromZero events
	}
}

void WCSP::propagateSeparator() {
	if (ToulBar2::verbose >= 2) cout << "PendingSeparator size: " << PendingSeparator.getSize() << endl;
	for (SeparatorList::iterator iter = PendingSeparator.begin(); iter != PendingSeparator.end(); ++iter) {
		(*iter)->propagate();
	}
}

void WCSP::propagateDEE() {
	if (ToulBar2::verbose >= 2) cout << "DEEQueue size: " << DEE.getSize() << endl;
	assert(NC.empty());
	while (!DEE.empty()) {
		EnumeratedVariable *x = (EnumeratedVariable *) DEE.pop();
		if (x->unassigned()) {
			if (ToulBar2::DEE>=3 || (ToulBar2::DEE==2 && getStore()->getDepth()==0)) {
				for (EnumeratedVariable::iterator itera = x->begin(); itera != x->end(); ++itera ) {
					for (EnumeratedVariable::iterator iterb = x->lower_bound(*itera +1); iterb != x->end(); ++iterb ) {
						assert(x->canbe(*itera));
					    assert(x->canbe(*iterb));
						assert(*itera != *iterb);
						x->propagateDEE(*itera, *iterb, false);
						if (x->cannotbe(*itera)) break;
					}
				}
			} else {
				Value a = x->getSupport();
				Value b = x->getMaxCostValue();
				assert(x->canbe(a) && x->getCost(a)==MIN_COST);
			    assert(x->canbe(b) && x->getCost(b)==x->getMaxCost());
				if (a == b) {
					assert(x->getMaxCost() == MIN_COST);
					if (a != x->getSup()) b = x->getSup();
					else b = x->getInf();
				}
				assert(a != b);
				x->propagateDEE(a,b);
			}
			propagateNC(); // DEE assumes NC already done
		}
	}
}

/// \defgroup varelim Variable elimination
/// - \e i-bounded variable elimination eliminates all variables with a degree less than or equal to \e i.
///		It can be done with arbitrary i-bound in preprocessing only and iff all their cost functions are in extension.
/// - \e i-bounded variable elimination with i-bound less than or equal to two can be done during the search.
/// - functional variable elimination eliminates all variables which have a bijective or functional binary hard constraint (\e ie ensuring a one-to-one or several-to-one value mapping) and iff all their cost functions are in extension.
///		It can be done without limit on their degree, in preprocessing only.
/// \note Variable elimination order used in preprocessing is either lexicographic or given by an external file *.order (see toulbar2 options)
/// \note 2-bounded variable elimination during search is optimal in the sense that any elimination order should result in the same final graph
/// \warning It is not possible to display/save solutions when bounded variable elimination is applied in preprocessing
/// \warning toulbar2 maintains a list of current cost functions for each variable.
///		It uses the size of these lists as an approximation of variable degrees.
///		During the search, if variable \e x has three cost functions \e xy, \e xz, \e xyz, its true degree is two but its approximate degree is three.
///		In toulbar2 options, it is the approximate degree which is given by the user for variable elimination during the search (thus, a value at most three).
///		But it is the true degree which is given by the user for variable elimination in preprocessing.

void WCSP::eliminate() {
	while (!Eliminate.empty()) {
		EnumeratedVariable *x = (EnumeratedVariable *) Eliminate.pop();
		if (x->unassigned()) {
			if (td) {
				if (td->isInCurrentClusterSubTree(x->getCluster())) x->eliminate();
			} else x->eliminate();
		}
	}
}

/// \defgroup softac Soft arc consistency and problem reformulation
/// Soft arc consistency is an incremental lower bound technique for optimization problems.
/// Its goal is to move costs from high-order (typically arity two or three) cost functions towards the problem lower bound and unary cost functions.
/// This is achieved by applying iteratively local equivalence-preserving problem transformations (EPTs) until some terminating conditions are met.
/// \note \e eg an EPT can move costs between a binary cost function and a unary cost function such that the sum of the two functions remains the same for any complete assignment.
/// \see <em> Arc consistency for Soft Constraints. </em> T. Schiex. Proc. of CP'2000. Singapour, 2000.
/// \note Soft Arc Consistency in toulbar2 is limited to binary and ternary and some global cost functions (\e eg alldifferent, gcc, regular, same).
///		Other n-ary cost functions are delayed for propagation until their number of unassigned variables is three or less.
/// \see <em> Towards Efficient Consistency Enforcement for Global Constraints in Weighted Constraint Satisfaction. </em> Jimmy Ho-Man Lee, Ka Lun Leung. Proc. of IJCAI 2009, pages 559-565. Pasadena, USA, 2009.

/// \defgroup propagation Propagation loop
/// Propagates soft local consistencies and bounded variable elimination until all the propagation queues are empty or a contradiction occurs.\n
/// While (queues are not empty or current objective bounds have changed):
/// -# queue for bounded variable elimination of degree at most two (except at preprocessing)
/// -# BAC queue
/// -# EAC queue
/// -# DAC queue
/// -# AC queue
/// -# global cost function propagation (not incremental)
/// -# NC queue
/// -# returns to #1 until all the previous queues are empty
/// -# VAC propagation (not incremental)
/// -# returns to #1 until all the previous queues are empty
/// -# exploits goods in pending separators for BTD-like methods
///
/// Queues are first-in / first-out lists of variables (avoiding multiple insertions).
/// In case of a contradiction, queues are explicitly emptied by WCSP::whenContradiction

void WCSP::propagate() {
	revise(NULL);
	if (ToulBar2::vac) vac->iniThreshold();

	for (vector<GlobalConstraint*>::iterator it = globalconstrs.begin(); it != globalconstrs.end(); it++) {
		(*(it))->init();
	}
	if (isGlobal()) {
		for (unsigned int i = 0; i < vars.size(); i++) {
			EnumeratedVariable* x = (EnumeratedVariable*) vars[i];
			if (x->unassigned()) {
				x->setCostProvidingPartition(); // For EAC propagation
			}
		}
	}

	do {
		do {
			do {
				do {
					eliminate();
					while (objectiveChanged || !NC.empty() || !IncDec.empty() || ((ToulBar2::LcLevel == LC_AC
							|| ToulBar2::LcLevel >= LC_FDAC) && !AC.empty()) || (ToulBar2::LcLevel >= LC_DAC
									&& !DAC.empty()) || (ToulBar2::LcLevel == LC_EDAC && !CSP(getLb(), getUb()) && !EAC1.empty())) {
						propagateIncDec();
						if (ToulBar2::LcLevel == LC_EDAC && !CSP(getLb(), getUb())) propagateEAC();
						assert(IncDec.empty());
						if (ToulBar2::LcLevel >= LC_DAC) propagateDAC();
						assert(IncDec.empty());
						if (ToulBar2::LcLevel == LC_AC || ToulBar2::LcLevel >= LC_FDAC) propagateAC();
						assert(IncDec.empty());

						Cost oldLb = getLb();
						bool cont = true;
						while (cont) {
							oldLb = getLb();
							cont = false;
							for (vector<GlobalConstraint*>::iterator it = globalconstrs.begin(); it != globalconstrs.end(); it++) {
								(*(it))->propagate();
								if (ToulBar2::LcLevel == LC_SNIC) if (!IncDec.empty()) cont = true; //For detecting value removal during SNIC enforcement
								propagateIncDec();
							}
							if (ToulBar2::LcLevel == LC_SNIC) if (!NC.empty() || objectiveChanged) cont = true; //For detecting value removal and upper bound change
							propagateNC();
							if (ToulBar2::LcLevel == LC_SNIC) if (oldLb != getLb() || !AC.empty()) {
								cont = true;
								AC.clear();//For detecting value removal and lower bound change
							}
						}
						propagateNC();
					}
				} while (!Eliminate.empty());

				if (ToulBar2::LcLevel < LC_EDAC || CSP(getLb(), getUb())) EAC1.clear();
				if (ToulBar2::vac) {
					//				assert(verify());
					if (vac->firstTime()) {
						vac->init();
						cout << "Lb before VAC: " << getLb() << endl;
					}
					vac->propagate();
				}
			} while (ToulBar2::vac && !vac->isVAC());
			if (ToulBar2::DEE) propagateDEE();
		} while (objectiveChanged || !NC.empty() || !IncDec.empty()
				 || ((ToulBar2::LcLevel == LC_AC || ToulBar2::LcLevel >= LC_FDAC) && !AC.empty())
				 || (ToulBar2::LcLevel >= LC_DAC && !DAC.empty())
			     || (ToulBar2::LcLevel == LC_EDAC && !CSP(getLb(), getUb()) && !EAC1.empty())
			     || !Eliminate.empty()
			     || (ToulBar2::vac && !vac->isVAC()));
		// TO BE DONE AFTER NORMAL PROPAGATION
		if (td) propagateSeparator();
	} while (objectiveChanged);
	revise(NULL);

	for (vector<GlobalConstraint*>::iterator it = globalconstrs.begin(); it != globalconstrs.end(); it++) {
		(*(it))->end();
	}
	assert(verify());
	assert(!objectiveChanged);
	assert(NC.empty());
	assert(IncDec.empty());
	if (ToulBar2::LcLevel == LC_AC || ToulBar2::LcLevel >= LC_FDAC) assert(AC.empty());
	else AC.clear();
	if (ToulBar2::LcLevel >= LC_DAC) assert(DAC.empty());
	else DAC.clear();
	assert(EAC1.empty());
	assert(EAC2.empty());
	assert(Eliminate.empty());
	if (ToulBar2::DEE) assert(DEE.empty());
	else DEE.clear();
	nbNodes++;
}

void WCSP::restoreSolution(Cluster* c) {
	int elimo = getElimOrder();
	for (int i = elimo - 1; i >= 0; i--) {
		elimInfo ei = elimInfos[i];
		EnumeratedVariable* x = (EnumeratedVariable*) ei.x;
		EnumeratedVariable* y = (EnumeratedVariable*) ei.y;
		EnumeratedVariable* z = (EnumeratedVariable*) ei.z;
		assert(x);
		assert(x->assigned());
		if (c && !c->isVar(x->wcspIndex)) continue;
		if (y && y->unassigned()) continue;
		if (z && z->unassigned()) continue;
		BinaryConstraint* xy = ei.xy;
		BinaryConstraint* xz = ei.xz;
		TernaryConstraint* xyz = ei.xyz;
		Constraint* ctr = ei.ctr;
		Value vy = -1;
		Value vz = -1;
		String tctr;
		Cost cctr;
		int xctrindex = -1;
		if (y) vy = getValue(y->wcspIndex);
		if (z) vz = getValue(z->wcspIndex);
		if (ctr) {ctr->firstlex(); ctr->nextlex(tctr, cctr); xctrindex=ctr->getIndex(x); assert(xctrindex>=0);}

		int minv = -1;
		Cost mincost = MAX_COST;
		for (unsigned int vxi = 0; vxi < x->getDomainInitSize(); vxi++) {
			Value vx = x->toValue(vxi);
			if (!x->canbeAfterElim(vx)) continue;
			Cost cxy = MIN_COST;
			Cost cxz = MIN_COST;
			Cost cxyz = MIN_COST;
			cctr = MIN_COST;
			if (xy) {
				if (xy->getIndex(y) == 0) cxy = xy->getCost(vy, vx);
				else cxy = xy->getCost(vx, vy);
			}
			if (xz) {
				if (xz->getIndex(z) == 0) cxz = xz->getCost(vz, vx);
				else cxz = xz->getCost(vx, vz);
			}
			if (xyz) cxyz = xyz->getCost(x, y, z, vx, vy, vz);
			if (ctr) {tctr[xctrindex]=vxi+CHAR_FIRST; cctr = ctr->evalsubstr(tctr, ctr);}
			Cost c = x->getCost(vx) + cxy + cxz + cxyz + cctr;
			//cout << "test " << vx << "," << x->getCost(vx) << "," << cxy << "," << cxz << "," << cxyz << endl;
			if (c < mincost) {
				mincost = c;
				minv = vx;
			}
		}
		//cout << i << ": elim " << x->getName() << "_" << minv << ", y= " << ((y)?y->getName():"-") << "_" << vy << ", z= " << ((z)?z->getName():"-") << "_" << vz << endl;
		x->assignWhenEliminated(minv);
	}
}

// -----------------------------------------------------------
// Methods for Variable Elimination

// Creates n fake empty constraints and puts them in the pool 'elimBinConstrs'

void WCSP::initElimConstr() {
	BinaryConstraint* xy = NULL;
	if (!ToulBar2::vac) xy = new BinaryConstraint(this, &storeData->storeCost);
	else xy = new VACBinaryConstraint(this, &storeData->storeCost);
	elimBinConstrs.push_back(xy);
	elimInfo ei = { NULL, NULL, NULL, NULL, NULL, NULL, NULL };
	elimInfos.push_back(ei);
}

void WCSP::initElimConstrs()
{
	for (unsigned int i = 0; i < vars.size(); i++) initElimConstr();

	vector<int> order;
	elimOrderFile2Vector(ToulBar2::varOrder, order);
	for (int i = vars.size()-1; i >= 0; --i) vars[order[i]]->queueEliminate();
	elimSpace = 0;
}

// Function that adds a new binary constraint from the pool of fake constraints
BinaryConstraint* WCSP::newBinaryConstr(EnumeratedVariable* x, EnumeratedVariable* y, Constraint *from1, Constraint *from2) {
	int newIndex = (int) elimBinOrder;
	BinaryConstraint* ctr = (BinaryConstraint*) elimBinConstrs[newIndex];
	ctr->fillElimConstr(x, y, from1, from2);
	if (ToulBar2::vac) ((VACBinaryConstraint *) ctr)->VACfillElimConstr();
	ctr->isDuplicate_ = false;
	return ctr;
}

// warning! you must create beforehand three binary constraints in fake pool (elimBinConstrs)
// if they do not exist in the main pool (constrs)
TernaryConstraint* WCSP::newTernaryConstr(EnumeratedVariable* x, EnumeratedVariable* y, EnumeratedVariable* z, Constraint *from1) {
	int newIndex = (int) elimTernOrder;
	TernaryConstraint* ctr = (TernaryConstraint*) elimTernConstrs[newIndex];
	ctr->fillElimConstr(x, y, z, from1);
	ctr->isDuplicate_ = false;
	return ctr;
}

TernaryConstraint* WCSP::newTernaryConstr(EnumeratedVariable* x, EnumeratedVariable* y, EnumeratedVariable* z, vector<
		Cost> costs) {
	unsigned int a, b;
	vector<Cost> zerocostsxy;
	vector<Cost> zerocostsxz;
	vector<Cost> zerocostsyz;

	for (a = 0; a < x->getDomainInitSize(); a++) {
		for (b = 0; b < y->getDomainInitSize(); b++) {
			zerocostsxy.push_back(MIN_COST);
		}
	}
	for (a = 0; a < x->getDomainInitSize(); a++) {
		for (b = 0; b < z->getDomainInitSize(); b++) {
			zerocostsxz.push_back(MIN_COST);
		}
	}
	for (a = 0; a < y->getDomainInitSize(); a++) {
		for (b = 0; b < z->getDomainInitSize(); b++) {
			zerocostsyz.push_back(MIN_COST);
		}
	}

	BinaryConstraint* xy = x->getConstr(y);
	BinaryConstraint* xz = x->getConstr(z);
	BinaryConstraint* yz = y->getConstr(z);

	if (!ToulBar2::vac) {
		if (!xy) {
			xy = new BinaryConstraint(this, x, y, zerocostsxy, &storeData->storeCost);
			xy->deconnect(true);
		}
		if (!xz) {
			xz = new BinaryConstraint(this, x, z, zerocostsxz, &storeData->storeCost);
			xz->deconnect(true);
		}
		if (!yz) {
			yz = new BinaryConstraint(this, y, z, zerocostsyz, &storeData->storeCost);
			yz->deconnect(true);
		}
	} else {
		if (!xy) {
			xy = new VACBinaryConstraint(this, x, y, zerocostsxy, &storeData->storeCost);
			xy->deconnect(true);
		}
		if (!xz) {
			xz = new VACBinaryConstraint(this, x, z, zerocostsxz, &storeData->storeCost);
			xz->deconnect(true);
		}
		if (!yz) {
			yz = new VACBinaryConstraint(this, y, z, zerocostsyz, &storeData->storeCost);
			yz->deconnect(true);
		}
	}

	TernaryConstraint *ctr = new TernaryConstraint(this, x, y, z, xy, xz, yz, costs, &storeData->storeCost);
	return ctr;
}

Constraint* WCSP::sum(Constraint* ctr1, Constraint* ctr2) {
	assert( ctr1 != ctr2 );
	if (ctr1->order(ctr2) < 0) {
		Constraint* ctraux = ctr1;
		ctr1 = ctr2;
		ctr2 = ctraux;
	}
	if (ToulBar2::verbose >= 1) cout << endl << "Sum of constraints: " << *ctr1 << " " << *ctr2 << endl;
	assert(ctr1->connected());
	assert(ctr2->connected());
	ctr1->deconnect();
	ctr2->deconnect(true);

	TSCOPE scopeUinv;
	TSCOPE scopeIinv;
	ctr1->scopeUnion(scopeUinv, ctr2);
	ctr1->scopeCommon(scopeIinv, ctr2);
	int arityU = scopeUinv.size();
	int arityI = scopeIinv.size();

	if (arityU == ctr2->arity() && ctr2->extension()) {
		ctr2->sumScopeIncluded(ctr1);
		ctr2->reconnect();
		ctr2->propagate();
		if (ToulBar2::verbose >= 1) cout << endl << "Scopes Included.  Has result: " << *ctr2 << endl;
		return ctr2;
	}

	EnumeratedVariable** scopeU = new EnumeratedVariable*[arityU];
	EnumeratedVariable** scopeI = new EnumeratedVariable*[arityI];
	int* scopeUi = new int[arityU];

	int i = 0;
	TSCOPE::iterator it = scopeIinv.begin();
	while (it != scopeIinv.end()) {
		int xi = it->first;
		assert(xi == vars[xi]->wcspIndex);
		scopeU[i] = (EnumeratedVariable*) vars[xi];
		scopeI[i] = scopeU[i];
		scopeUi[i] = vars[xi]->wcspIndex;
		it++;
		i++;
		scopeUinv.erase(xi);
	}
	it = scopeUinv.begin();
	while (it != scopeUinv.end()) {
		int xi = it->first;
		scopeU[i] = (EnumeratedVariable*) vars[xi];
		scopeUi[i] = vars[xi]->wcspIndex;
		it++;
		i++;
	}

	EnumeratedVariable* x = scopeU[0];
	EnumeratedVariable* y = scopeU[1];

	Cost Top = getUb();
	unsigned int vxi, vyi, vzi;
	String tuple, tuple1, tuple2;
	Cost cost, cost1, cost2;
	int ctrIndex = -INT_MAX;
	Constraint* ctr = NULL;
	vector<Cost> costs;

	if (arityU > 3) { // || isGlobal()) {
		ctrIndex = postNaryConstraintBegin(scopeUi, arityU, Top);
		ctr = getCtr(ctrIndex);
		assert(ctr->extension());
		NaryConstraintMap* nary = (NaryConstraintMap*) ctr;

		nary->fillFilters();

		bool tupleXtuple = (ctr1->getDefCost() >= Top) && (ctr2->getDefCost() >= Top);

		if (tupleXtuple) {
			ctr1->first();
			while (ctr1->next(tuple1, cost1)) {
				ctr2->first();
				while (ctr2->next(tuple2, cost2)) {
					nary->insertSum(tuple1, cost1, ctr1, tuple2, cost2, ctr2, true);
				}
			}
		} else {
			nary->firstlex();
			while (nary->nextlex(tuple, cost)) {
				cost1 = ctr1->evalsubstr(tuple, nary);
				cost2 = ctr2->evalsubstr(tuple, nary);
				if (cost1 + cost2 < Top) nary->setTuple(tuple, cost1 + cost2);
			}
		}
	} else if (arityU == 3) {
		EnumeratedVariable* z = scopeU[2];
		for (vxi = 0; vxi < x->getDomainInitSize(); vxi++)
			for (vyi = 0; vyi < y->getDomainInitSize(); vyi++)
				for (vzi = 0; vzi < z->getDomainInitSize(); vzi++) {
					Value vx = x->toValue(vxi);
					Value vy = y->toValue(vyi);
					Value vz = z->toValue(vzi);
					Cost costsum = Top;
					if (x->canbe(vx) && y->canbe(vy) && z->canbe(vz)) {
						costsum = MIN_COST;
						if (arityI == 1) costsum += ((BinaryConstraint*) ctr1)->getCost(x, y, vx, vy)
								+ ((BinaryConstraint*) ctr2)->getCost(x, z, vx, vz);
						else if (arityI == 2) costsum += ((BinaryConstraint*) ctr1)->getCost(x, y, vx, vy)
								+ ((TernaryConstraint*) ctr2)->getCost(x, y, z, vx, vy, vz);
						else if (arityI == 3) costsum += ((TernaryConstraint*) ctr1)->getCost(x, y, z, vx, vy, vz)
								+ ((TernaryConstraint*) ctr2)->getCost(x, y, z, vx, vy, vz);
else						assert(false);
						if(costsum > Top) costsum = Top;
					}
					costs.push_back(costsum);
				}
				ctrIndex = postTernaryConstraint( x->wcspIndex, y->wcspIndex, z->wcspIndex, costs );
			}
			else if(arityU == 2) {
				BinaryConstraint* bctr1 = (BinaryConstraint*) ctr1;
				BinaryConstraint* bctr2 = (BinaryConstraint*) ctr2;
				for (vxi = 0; vxi < x->getDomainInitSize(); vxi++)
				for (vyi = 0; vyi < y->getDomainInitSize(); vyi++) {
					Value vx = x->toValue(vxi);
					Value vy = y->toValue(vyi);
					Cost costsum = Top;
					if(x->canbe(vx) && y->canbe(vy)) {
						Cost costsum = bctr1->getCost(x,y,vx,vy) + bctr2->getCost(x,y,vx,vy);
						if(costsum > Top) costsum = Top;
					}
					costs.push_back(costsum);
				}
				ctrIndex= postBinaryConstraint( x->wcspIndex, y->wcspIndex, costs );
			}
			assert(ctrIndex > -INT_MAX);
			delete [] scopeU;
			delete [] scopeUi;
			delete [] scopeI;
			ctr = getCtr(ctrIndex);
			ctr->propagate();
			if (ToulBar2::verbose >= 1) cout << endl << "Has result: " << *ctr << endl;
			return ctr;
		}

void WCSP::project(Constraint* &ctr_inout, EnumeratedVariable* var) {
	if (ctr_inout->getIndex(var) < 0) return;
	unsigned int vxi, vyi, vzi;
	int arity = ctr_inout->arity();
	int truearity = arity;
	if (arity >= 5) { // n-ary (n>=5) with only 4 unassigned variables is projected on a ternary!
		truearity = 0;
		for (int i = 0; i < arity; i++)
			if (ctr_inout->getVar(i)->unassigned()) truearity++;
		assert(truearity>=4);
	}

	if (ToulBar2::verbose >= 1) cout << endl << "Projection of var " << var->wcspIndex << " in ctr: " << *ctr_inout
			<< endl;

	if (truearity - 1 > 3) {
		((NaryConstraint*) ctr_inout)->project(var);
		ctr_inout->propagate();
		if (ToulBar2::verbose >= 1) cout << endl << "   has result*: " << *ctr_inout << endl;
		return;
	}
	ctr_inout->deconnect();

	int i, j;
	int ivars[3];
	EnumeratedVariable* evars[3];

	j = 0;
	for (i = 0; i < arity; i++) {
		EnumeratedVariable* v = (EnumeratedVariable*) ctr_inout->getVar(i);
		if (v != var && (arity <= 4 || v->unassigned())) {
			ivars[j] = v->wcspIndex;
			evars[j] = v;
			j++;
		}
	}

	Constraint* ctr = NULL;
	Cost Top = MAX_COST; // getUb();
	int ctrIndex;
	Char t[arity + 1];
	vector<Cost> costs;
	Cost negcost = 0;

	switch (truearity - 1) {
	case 3: {
		NaryConstraint *nctr = (NaryConstraint*) ctr_inout;
		if (truearity == 4 && arity >= 5) {
			for (i = 0; i < arity; i++) {
				if (ctr_inout->getVar(i)->assigned()) t[ctr_inout->getIndex(ctr_inout->getVar(i))] = ((EnumeratedVariable*) ctr_inout->getVar(i))->toIndex(ctr_inout->getVar(i)->getValue()) + CHAR_FIRST;
			}
		}
		for (vxi = 0; vxi < evars[0]->getDomainInitSize(); vxi++) {
		  for (vyi = 0; vyi < evars[1]->getDomainInitSize(); vyi++) {
			for (vzi = 0; vzi < evars[2]->getDomainInitSize(); vzi++) {
			  Value v0 = evars[0]->toValue(vxi);
			  Value v1 = evars[1]->toValue(vyi);
			  Value v2 = evars[2]->toValue(vzi);
			  Cost mincost = Top;
			  if (evars[0]->canbe(v0) && evars[1]->canbe(v1) && evars[2]->canbe(v2)) {
				t[ctr_inout->getIndex(evars[0])] = vxi + CHAR_FIRST;
				t[ctr_inout->getIndex(evars[1])] = vyi + CHAR_FIRST;
				t[ctr_inout->getIndex(evars[2])] = vzi + CHAR_FIRST;

				for (EnumeratedVariable::iterator itv = var->begin(); itv != var->end(); ++itv) {
				  t[ctr_inout->getIndex(var)] = var->toIndex(*itv) + CHAR_FIRST;
				  t[arity] = '\0';
				  String strt(t);
				  Cost c = nctr->eval(strt) + var->getCost(*itv);
				  if (ToulBar2::isZ) mincost = SumLogLikeCost(mincost,c);
				  else if (c < mincost) mincost = c;
				}
			  }
			  if (ToulBar2::isZ && mincost < negcost) negcost = mincost;
			  costs.push_back(mincost);
			}
		  }
		}
		assert(negcost <= 0);
		if (negcost < 0) {
		  for (vxi = 0; vxi < evars[0]->getDomainInitSize(); vxi++) {
			for (vyi = 0; vyi < evars[1]->getDomainInitSize(); vyi++) {
			  for (vzi = 0; vzi < evars[2]->getDomainInitSize(); vzi++) {
				costs[vxi * evars[1]->getDomainInitSize() * evars[2]->getDomainInitSize() + vyi * evars[2]->getDomainInitSize() + vzi] -= negcost;
			  }
			}
		  }
		  decreaseLb(negcost);
		}
		ctrIndex = postTernaryConstraint(ivars[0], ivars[1], ivars[2], costs);
		ctr = getCtr(ctrIndex);
	    }
		break;
	case 2: {
		TernaryConstraint *tctr = (TernaryConstraint*) ctr_inout;
		for (vxi = 0; vxi < evars[0]->getDomainInitSize(); vxi++)
			for (vyi = 0; vyi < evars[1]->getDomainInitSize(); vyi++) {
				Value v0 = evars[0]->toValue(vxi);
				Value v1 = evars[1]->toValue(vyi);
				Cost mincost = Top;
				if (evars[0]->canbe(v0) && evars[1]->canbe(v1)) {
					for (EnumeratedVariable::iterator itv = var->begin(); itv != var->end(); ++itv) {
						Cost c = tctr->getCost(evars[0], evars[1], var, v0, v1, *itv) + var->getCost(*itv);
						if (ToulBar2::isZ) mincost = SumLogLikeCost(mincost,c);
						else if (c < mincost) mincost = c;
					}
				}
				if (ToulBar2::isZ && mincost < negcost) negcost = mincost;
				costs.push_back(mincost);
			}
		assert(negcost <= 0);
		if (negcost < 0) {
		  for (vxi = 0; vxi < evars[0]->getDomainInitSize(); vxi++) {
			for (vyi = 0; vyi < evars[1]->getDomainInitSize(); vyi++) { 
			  costs[vxi * evars[1]->getDomainInitSize() + vyi] -= negcost;
			}
		  }
		  decreaseLb(negcost);
		}
		ctrIndex = postBinaryConstraint(ivars[0], ivars[1], costs);
		ctr = getCtr(ctrIndex);
	    }
		break;
	case 1: {
	    BinaryConstraint *bctr =  ((BinaryConstraint*) ctr_inout);
		for (vxi = 0; vxi < evars[0]->getDomainInitSize(); vxi++) {
			Value v0 = evars[0]->toValue(vxi);
			Cost mincost = Top;
			if (evars[0]->canbe(v0)) {
				for (EnumeratedVariable::iterator itv = var->begin(); itv != var->end(); ++itv) {
					Cost c = bctr->getCost(evars[0], var, v0, *itv) + var->getCost(*itv);
					if (ToulBar2::isZ) mincost = SumLogLikeCost(mincost,c);
					else if (c < mincost) mincost = c;
				}
			}
			if (ToulBar2::isZ && mincost < negcost) negcost = mincost;
			costs.push_back(mincost);
		}
		assert(negcost <= 0);
		for (EnumeratedVariable::iterator itv0 = evars[0]->begin(); itv0 != evars[0]->end(); ++itv0) {
			vxi = evars[0]->toIndex(*itv0);
			if (costs[vxi] - negcost > MIN_COST) evars[0]->project(*itv0, costs[vxi] - negcost);
		}
		evars[0]->findSupport();
		if (negcost < 0) decreaseLb(negcost);
	    }
		break;
	default:
		cerr << "Bad resulting cost function arity during generic variable elimination!";
		exit(EXIT_FAILURE);
	}
	ctr_inout = ctr;
	if (ctr) {
		ctr->propagate();
		if (ToulBar2::verbose >= 1) cout << endl << "   has result: " << *ctr_inout << endl;
	}
}

void WCSP::variableElimination(EnumeratedVariable* var) {
	int degree = var->getTrueDegree();
	if (ToulBar2::verbose >= 1) cout << endl << "Generic variable elimination of " << var->getName() << "    degree: "
			<< var->getDegree() << " true degree: " << degree << endl;
	if (degree > maxDegree) maxDegree = degree;

	if (var->getDegree() > 0) {

		ConstraintList::iterator it1 = var->getConstrs()->begin();
		ConstraintList::iterator it2;
		Constraint* c1 = (*it1).constr;
		Constraint* c2 = NULL;
		Constraint* csum = c1;
		Constraint* csumcopy = NULL;

		while (var->getDegree() > 1) {
			it1 = var->getConstrs()->begin();
			it2 = var->getConstrs()->rbegin();
			c1 = (*it1).constr;
			c2 = (*it2).constr;
			csum = sum(c1, c2);

			if (getTreeDec()) {
				csum->setCluster(var->getCluster());
			}
		}

		assert(csum->getIndex(var) >= 0);
		csumcopy = csum->copy();
		assert(csumcopy != NULL);
		elimInfo ei = {var,NULL,NULL,NULL,NULL,NULL,csumcopy};
		elimInfos[getElimOrder()] = ei;
		elimOrderInc();
		elimSpace += csumcopy->space();
		project(csum, var);
	} else {
	  if (ToulBar2::isZ) { // add all unary loglike into lowerbound or negCost
		Cost clogz = getUb();
		for (EnumeratedVariable::iterator itv = var->begin(); itv != var->end(); ++itv) {
		  clogz = SumLogLikeCost(clogz, var->getCost(*itv));
		}
		if (clogz < 0) decreaseLb(clogz);
		else increaseLb(clogz);
	  }
	}
	assert(var->getDegree() == 0);

	var->assign(var->getSupport()); // warning! dummy assigned value
}

bool WCSP::kconsistency(int xIndex, int yIndex, int zIndex, BinaryConstraint* xy, BinaryConstraint* yz, BinaryConstraint* xz) {
	if ((xIndex == yIndex) || (xIndex == zIndex) || (yIndex == zIndex)) return false;
	EnumeratedVariable* x = (EnumeratedVariable *) vars[xIndex];
	EnumeratedVariable* y = (EnumeratedVariable *) vars[yIndex];
	EnumeratedVariable* z = (EnumeratedVariable *) vars[zIndex];
	TernaryConstraint* tctr = x->getConstr(y, z);
	if (tctr) return false;

	bool added = false;
	vector<Cost> costs;
	Cost ub = getUb();
	Cost minc = ub;
	for (unsigned int a = 0; a < x->getDomainInitSize(); a++) {
		Value va = x->toValue(a);
		Cost costa = x->getCost(va);
		for (unsigned int b = 0; b < y->getDomainInitSize(); b++) {
			Value vb = y->toValue(b);
			Cost costb = y->getCost(vb);
			for (unsigned int c = 0; c < z->getDomainInitSize(); c++) {
				Value vc = z->toValue(c);
				Cost costc = z->getCost(vc);
				Cost ctuple = ub;
				if (x->canbe(va) && y->canbe(vb) && z->canbe(vc)) {
					ctuple = costa + costb + costc + xy->getCost(x, y, va, vb) + xz->getCost(x, z, va, vc)
							+ yz->getCost(y, z, vb, vc);
				}
				if (ctuple < minc) minc = ctuple;
				costs.push_back(ctuple);
			}
		}
	}

	if (minc > MIN_COST) {
		tctr = new TernaryConstraint(this, &storeData->storeCost);
		elimTernConstrs.push_back(tctr);
		tctr = newTernaryConstr(x, y, z);
		tctr->fillElimConstrBinaries();
		tctr->reconnect();
		elimTernOrderInc();

		vector<Cost>::iterator it = costs.begin();
		for (unsigned int a = 0; a < x->getDomainInitSize(); a++) {
			Value va = x->toValue(a);
			Cost costa = x->getCost(va);
			if (x->canbe(va)) x->extend(va, costa);
			for (unsigned int b = 0; b < y->getDomainInitSize(); b++) {
				Value vb = y->toValue(b);
				Cost costb = y->getCost(vb);
				if (y->canbe(vb)) y->extend(vb, costb);
				for (unsigned int c = 0; c < z->getDomainInitSize(); c++) {
					Value vc = z->toValue(c);
					Cost costc = z->getCost(vc);
					if (z->canbe(vc)) z->extend(vc, costc);
					if (x->canbe(va) && y->canbe(vb)) {
						Cost costab = xy->getCost(x, y, va, vb);
						xy->addcost(x, y, va, vb, -costab);
					}
					if (y->canbe(vb) && z->canbe(vc)) {
						Cost costbc = yz->getCost(y, z, vb, vc);
						yz->addcost(y, z, vb, vc, -costbc);
					}
					if (x->canbe(va) && z->canbe(vc)) {
						Cost costac = xz->getCost(x, z, va, vc);
						xz->addcost(x, z, va, vc, -costac);
					}
					tctr->setcost(x, y, z, va, vb, vc, *it - minc);
					++it;
				}
			}
		}
		tctr->projectTernary();
		increaseLb(minc);
		if (ToulBar2::verbose >= 1) cout << "new ternary(" << x->wcspIndex << "," << y->wcspIndex << ","
				<< z->wcspIndex << ")  newLb: " << getLb() << endl;
		added = true;
	}
	return added;
}

void WCSP::ternaryCompletion() {
    if (numberOfUnassignedVariables()<3) return;
	Double nbunvars  = numberOfUnassignedVariables();
	Double connectivity = 2. * numberOfConnectedBinaryConstraints() / (nbunvars * (nbunvars-1));
	Double domsize =  medianDomainSize();
	Double size =  domsize;
	size = sizeof(StoreCost) * size * size * size * nbunvars * (nbunvars-1) * (nbunvars-2) * connectivity * connectivity * connectivity;
	if (ToulBar2::debug>=2) cout << "MAX ESTIMATED RPC SIZE: " << size << " (" << nbunvars << "," << connectivity <<")" << endl;
	if (size > 1024. * 1024. * ToulBar2::preprocessTernaryRPC) {
	  cout << "Restricted path consistency disabled (" << size/1024./1024. << " >= " << ToulBar2::preprocessTernaryRPC << " MB)" << endl;
	  ToulBar2::preprocessTernaryRPC = 0;
	  return;
	}
	Double totalsize = 0.;
	Double maxsize = 1024. * 1024. * ToulBar2::preprocessTernaryRPC;
	int ntern = 0;
	for (unsigned int i = 0; i < vars.size(); i++) {
		EnumeratedVariable* x = (EnumeratedVariable*) vars[i];
		for (ConstraintList::iterator it = x->getConstrs()->begin(); it != x->getConstrs()->end(); ++it) {
			Constraint* ctr = (*it).constr;
			if (ctr->arity() == 2) {
				BinaryConstraint* bctr = (BinaryConstraint*) ctr;
				EnumeratedVariable* y = (EnumeratedVariable*) bctr->getVarDiffFrom(x);
				if (y->wcspIndex > x->wcspIndex) for (ConstraintList::iterator it2 = y->getConstrs()->begin(); it2
						!= y->getConstrs()->end(); ++it2) {
					Constraint* ctr2 = (*it2).constr;
					if (ctr != ctr2 && ctr2->arity() == 2) {
						BinaryConstraint* bctr2 = (BinaryConstraint*) ctr2;
						EnumeratedVariable* z = (EnumeratedVariable*) bctr2->getVarDiffFrom(y);
						if (z->wcspIndex > y->wcspIndex) for (ConstraintList::iterator it3 = z->getConstrs()->begin(); it3
								!= z->getConstrs()->end(); ++it3) {
							Constraint* ctr3 = (*it3).constr;
							if (ctr2 != ctr3 && ctr3->arity() == 2) {
								BinaryConstraint* bctr3 = (BinaryConstraint*) ctr3;
								EnumeratedVariable * xx = (EnumeratedVariable*) bctr3->getVarDiffFrom(z);
								if (x == xx) {
									// bool added = kconsistency(x->wcspIndex,
									// 		y->wcspIndex, z->wcspIndex, bctr,
									// 		bctr2, bctr3);
									// if (added)
								    int xsize = x->getDomainInitSize();
									int ysize = y->getDomainInitSize();
									int zsize = z->getDomainInitSize();
									totalsize += (Double) sizeof(StoreCost) * xsize * ysize * zsize;
									if (totalsize > maxsize) {
									  cout << "Restricted path consistency stopped (" << totalsize/1024./1024. << " >= " << ToulBar2::preprocessTernaryRPC << " MB)" << endl;
									  cout << "add " << ntern << " ternary ctrs." << endl;
									  ToulBar2::preprocessTernaryRPC = 0;
									  return;
									}
									vector<Cost> costs(xsize * ysize * zsize, MIN_COST);
									//									cout << x->wcspIndex << y->wcspIndex << z->wcspIndex << endl;
									postTernaryConstraint(x->wcspIndex, y->wcspIndex, z->wcspIndex, costs);
									ntern++;
								}
							}
						}
					}
				}
			}

		}
	}
	//	cout << "total ctrs: " << numberOfConnectedConstraints() << endl;
	cout << "add " << ntern << " ternary ctrs." << endl;
}

// -----------------------------------------------------------
// Methods for Virtual Arc Consistency

void WCSP::histogram(Cost c) {
	if (vac) vac->histogram(c);
}
void WCSP::histogram() {
	if (vac) vac->histogram();
}
void WCSP::iniSingleton() {
	if (vac) vac->iniSingleton();
}
void WCSP::updateSingleton() {
	if (vac) vac->updateSingleton();
}
void WCSP::removeSingleton() {
	if (vac) vac->removeSingleton();
}
void WCSP::printVACStat() {
	if (vac) vac->printStat();
}
int WCSP::getVACHeuristic() {
	if (vac) return vac->getHeuristic();
	else return -1;
}

// -----------------------------------------------------------
// Methods for Cluster Tree Decomposition

void WCSP::buildTreeDecomposition() {
	td = new TreeDecomposition(this);
	double time = cpuTime();
	if (ToulBar2::approximateCountingBTD) td->buildFromOrderForApprox();
	else td->buildFromOrder();
	if(ToulBar2::verbose >= 0) cout << "Tree decomposition time: " << cpuTime() - time << " seconds." << endl;
	if (!ToulBar2::approximateCountingBTD) {
		vector<int> order;
		td->getElimVarOrder(order);
		// allows propagation to operate on the whole problem without modifying tree decomposition local lower bounds and delta costs
		// it is important for RDS-BTD which assumes zero cluster lower bounds and no delta cost moves
		TreeDecomposition *tmptd = td;
		td = NULL;
		setDACOrder(order);
		td = tmptd;
		// new constraints may be produced by variable elimination that must be correctly assigned to a cluster
		for (unsigned int i=0; i<numberOfConstraints(); i++) if (constrs[i]->getCluster()==-1) constrs[i]->assignCluster();
		for (int i=0; i<elimBinOrder; i++) if (elimBinConstrs[i]->connected() && elimBinConstrs[i]->getCluster()==-1) elimBinConstrs[i]->assignCluster();
		for (int i=0; i<elimTernOrder; i++) if (elimTernConstrs[i]->connected() && elimTernConstrs[i]->getCluster()==-1) elimTernConstrs[i]->assignCluster();
		// check if ternary constraint cluster assignments are valid
#ifndef NDEBUG
		for (unsigned int i=0; i<numberOfConstraints(); i++) {
		  Constraint* ctr = getCtr(i);
		  if (ctr->connected() && !ctr->isSep()) {
			if(ctr->arity() == 3) {
			  TernaryConstraint* tctr = (TernaryConstraint*) ctr;
			  //			  tctr->setDuplicates();
			  assert(tctr->xy->getCluster() == tctr->getCluster() &&
					 tctr->xz->getCluster() == tctr->getCluster() &&
					 tctr->yz->getCluster() == tctr->getCluster() );
			}
		  }
		}
		for (int i=0; i<elimTernOrder; i++) if (elimTernConstrs[i]->connected()) {
			Constraint* ctr = elimTernConstrs[i];
			if (ctr->connected() && !ctr->isSep()) {
			  if(ctr->arity() == 3) {
				TernaryConstraint* tctr = (TernaryConstraint*) ctr;
				//				tctr->setDuplicates();
				assert(tctr->xy->getCluster() == tctr->getCluster() &&
					   tctr->xz->getCluster() == tctr->getCluster() &&
					   tctr->yz->getCluster() == tctr->getCluster() );
			  }
			}
		  }
#endif
	}
}

void WCSP::elimOrderFile2Vector(char *elimVarOrder, vector<int> &order)
{
	ifstream file;
	if (elimVarOrder) file.open(elimVarOrder);
	if (!elimVarOrder || !file) {
		if (ToulBar2::verbose >= 1) {
			cout << "Variable elimination order file missing or unreadable... takes reverse problem file variable index order." << endl;
		}
		//		for(unsigned int i=0;i<numberOfVariables();i++) order.push_back(i);
		for(int i=numberOfVariables()-1; i>=0; i--) order.push_back(i);
	} else {
		while (file) {
			int ix;
			file >> ix;
			if (file) order.push_back(ix);
		}
	}
	if (order.size() != numberOfVariables()) {
		cerr << "Variable elimination order file has incorrect number of variables." << endl;
		exit(EXIT_FAILURE);
	}
}

/// \param order variable elimination order (reverse of DAC order)
void WCSP::setDACOrder(vector<int> &order) {
	if (order.size() != numberOfVariables()) {
		cerr << "DAC order has incorrect number of variables." << endl;
		exit(EXIT_FAILURE);
	}

	// set DAC order to the inverse of the elimination variable ordering
	if (ToulBar2::verbose >= 1) cout << "DAC order:";
	for (int i = order.size() - 1; i >= 0; i--) {
		if (ToulBar2::verbose >= 1) cout << " " << getVar(order[i])->getName();
		getVar(order[i])->setDACOrder(order.size() - 1 - i);
		if (ToulBar2::DEE>=2) getVar(order[i])->queueDEE();
	}
	if (ToulBar2::verbose >= 1) cout << endl;

	for (unsigned int i = 0; i < numberOfConstraints(); i++) {
		Constraint* ctr = getCtr(i);
		ctr->setDACScopeIndex();
		 // Postpone global constraint propagation at the end (call to WCSP::propagate())
		if (ctr->connected() && !ctr->isGlobal()) ctr->propagate();
	}
	for (int i = 0; i < elimBinOrder; i++) {
		Constraint* ctr = elimBinConstrs[i];
		ctr->setDACScopeIndex();
		if (ctr->connected()) ctr->propagate();
	}
	for (int i = 0; i < elimTernOrder; i++) {
		Constraint* ctr = elimTernConstrs[i];
		ctr->setDACScopeIndex();
		if (ctr->connected()) ctr->propagate();
	}
	propagate();
}

// -----------------------------------------------------------
// Functions for dealing with probabilities
// Warning: ToulBar2::NormFactor has to be initialized

Cost WCSP::Prob2Cost(TProb p) const {
	if (p == 0.0) return getUb();
	TProb res = -Log10(p) * ToulBar2::NormFactor;
	if (res > to_double(MAX_COST)) {
		cerr << "Overflow when converting probability to cost." << endl;
		exit(EXIT_FAILURE);
	}
	Cost c = (Cost) res;
	if (c > getUb()) return getUb();
	return c;
}

Cost WCSP::LogLike2Cost(TProb p) const {
	TProb res = -p * ToulBar2::NormFactor;
	Cost c = (Cost) res;
	return c;	
}
TProb WCSP::Cost2LogLike(Cost c) const {
	return -to_double(c) / ToulBar2::NormFactor;
}
TProb WCSP::Cost2Prob(Cost c) const {
  return Exp10(-to_double(c) / ToulBar2::NormFactor);
}
Cost  WCSP::SumLogLikeCost(Cost c1, Cost c2) const {
  if (c1 >= getUb()) return c2;
  else if (c2 >= getUb()) return c1;
  else if (c1 == c2) return c1+LogLike2Cost(Log10(2.));
  else {
	if (c1 < c2) return c1 + LogLike2Cost(Log1p(Exp10(Cost2LogLike(c2 - c1)))/Log(10.));
	else return c2 + LogLike2Cost(Log1p(Exp10(Cost2LogLike(c1 - c2)))/Log(10.));
  }
}
TProb WCSP::SumLogLikeCost(TProb logc1, Cost c2) const {
  TProb logc2 = Cost2LogLike(c2);
  if (logc1 == -numeric_limits<TProb>::infinity()) return logc2;
  else if (c2 >= getUb()) return logc1;
  else {
	if (logc1 >= logc2) return logc1 + (Log1p(Exp10(logc2 - logc1))/Log(10.));
	else return logc2 + (Log1p(Exp10(logc1 - logc2))/Log(10.));
  }
}
TProb WCSP::SumLogLikeCost(TProb logc1, TProb logc2) const {
  if (logc1 == -numeric_limits<TProb>::infinity()) return logc2;
  else if (logc2 == -numeric_limits<TProb>::infinity()) return logc1;
  else {
	if (logc1 >= logc2) return logc1 + (Log1p(Exp10(logc2 - logc1))/Log(10.));
	else return logc2 + (Log1p(Exp10(logc1 - logc2))/Log(10.));
  }
}

//----------------------------------------
//procedure when berge acycl constant are present in the problem
// toulbar2::Berge_Dec has to be initialized > 0;

void WCSP::visit(int i, vector <int>&revdac, vector <bool>& marked, const vector< vector<int> >&listofsuccessors ) {
   marked[i] = true;
  for (unsigned int  j = 0 ; j < listofsuccessors[i].size(); j++) {
//   for (int  j = listofsuccessors[i].size()-1 ; j >= 0 ; j--) {
      if (!marked[listofsuccessors[i][j]]) visit(listofsuccessors[i][j],revdac,marked,listofsuccessors) ;
    }
   revdac.push_back(i);
}



