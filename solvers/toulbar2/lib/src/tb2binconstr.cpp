/*
 * ****** Binary constraint applied on variables with enumerated domains ******
 */

#include "tb2binconstr.hpp"
#include "tb2wcsp.hpp"
#include "tb2clusters.hpp"

/*
 * Constructors and misc.
 *
 */

BinaryConstraint::BinaryConstraint(WCSP *wcsp, EnumeratedVariable *xx, EnumeratedVariable *yy, vector<Cost> &tab, StoreStack<Cost, Cost> *storeCost) :
        AbstractBinaryConstraint<EnumeratedVariable,EnumeratedVariable>(wcsp, xx, yy), sizeX(xx->getDomainInitSize()), sizeY(yy->getDomainInitSize())
{
    deltaCostsX = vector<StoreCost>(sizeX,StoreCost(MIN_COST,storeCost));
    deltaCostsY = vector<StoreCost>(sizeY,StoreCost(MIN_COST,storeCost));
    assert(tab.size() == sizeX * sizeY);
    supportX = vector<Value>(sizeX,y->getInf());
    supportY = vector<Value>(sizeY,x->getInf());

	costs = vector<StoreCost>(sizeX*sizeY,StoreCost(MIN_COST,storeCost));

    for (unsigned int a = 0; a < x->getDomainInitSize(); a++)
         for (unsigned int b = 0; b < y->getDomainInitSize(); b++)
                costs[a * sizeY + b] = tab[a * sizeY + b];

    propagate();
}

BinaryConstraint::BinaryConstraint(WCSP *wcsp, StoreStack<Cost, Cost> *storeCost)
 	: AbstractBinaryConstraint<EnumeratedVariable,EnumeratedVariable>(wcsp), sizeX(0), sizeY(0)
{
//	unsigned int maxdomainsize = wcsp->getMaxDomainSize();
//    deltaCostsX = vector<StoreCost>(maxdomainsize,StoreCost(MIN_COST,storeCost));
//    deltaCostsY = vector<StoreCost>(maxdomainsize,StoreCost(MIN_COST,storeCost));
//    supportX = vector<Value>(maxdomainsize,0);
//    supportY = vector<Value>(maxdomainsize,0);
    linkX = new DLink<ConstraintLink>;
    linkY = new DLink<ConstraintLink>;

//    costs = vector<StoreCost>(maxdomainsize*maxdomainsize,StoreCost(MIN_COST,storeCost));
//    for (unsigned int a = 0; a < maxdomainsize; a++)
//         for (unsigned int b = 0; b < maxdomainsize; b++)
//                costs[a * maxdomainsize + b] = MIN_COST;
}

void BinaryConstraint::print(ostream& os)
{
    os << this << " BinaryConstraint(" << x->getName() << "," << y->getName() << ")";
	if (ToulBar2::weightedDegree) os << "/" << getConflictWeight();
    if(wcsp->getTreeDec()) cout << "   cluster: " << getCluster() << endl;
    else cout << endl;
    if (ToulBar2::verbose >= 5) {
        for (EnumeratedVariable::iterator iterX = x->begin(); iterX != x->end(); ++iterX) {
            for (EnumeratedVariable::iterator iterY = y->begin(); iterY != y->end(); ++iterY) {
                os << " " << getCost(*iterX, *iterY);
            }
            os << endl;
        }
    }
}

void BinaryConstraint::dump(ostream& os, bool original)
{
  os << "2 " << ((original)?(x->wcspIndex):x->getCurrentVarId()) << " " << ((original)?(y->wcspIndex):y->getCurrentVarId()) << " " << MIN_COST << " " << x->getDomainSize() * y->getDomainSize() << endl;
  int i=0;
  for (EnumeratedVariable::iterator iterX = x->begin(); iterX != x->end(); ++iterX, i++) {
	int j=0;
	for (EnumeratedVariable::iterator iterY = y->begin(); iterY != y->end(); ++iterY, j++) {
	  os << ((original)?(*iterX):i) << " " << ((original)?(*iterY):j) << " " << ((original)?getCost(*iterX, *iterY):min(wcsp->getUb(),getCost(*iterX, *iterY))) << endl;
	}
  }
}

/*
 * Propagation methods
 *
 */
bool BinaryConstraint::project(EnumeratedVariable *x, Value value, Cost cost, vector<StoreCost> &deltaCostsX)
{
	assert(ToulBar2::verbose < 4 || ((cout << "project(C" << getVar(0)->getName() << "," << getVar(1)->getName() << ", (" << x->getName() << "," << value << "), " << cost << ")" << endl), true));
    
    // hard binary constraint costs are not changed
    if (!CUT(cost + wcsp->getLb(), wcsp->getUb())) {
	    TreeDecomposition* td = wcsp->getTreeDec();
	    if(td) td->addDelta(cluster,x,value,cost);
    	deltaCostsX[x->toIndex(value)] += cost;  // Warning! Possible overflow???
    }
    	
    Cost oldcost = x->getCost(value);
    x->project(value, cost);
#ifdef DEECOMPLETE
    getVarDiffFrom(x)->queueDEE();
#endif
    return (x->getSupport() == value || SUPPORTTEST(oldcost, cost));
}


void BinaryConstraint::extend(EnumeratedVariable *x, Value value, Cost cost, vector<StoreCost> &deltaCostsX)
{
	assert(ToulBar2::verbose < 4 || ((cout << "extend(C" << getVar(0)->getName() << "," << getVar(1)->getName() << ", (" << x->getName() << "," << value << "), " << cost << ")" << endl), true));

    TreeDecomposition* td = wcsp->getTreeDec();
    if(td) td->addDelta(cluster,x,value,-cost);

    deltaCostsX[x->toIndex(value)] -= cost;  // Warning! Possible overflow???
    x->extend(value, cost);
}

template <typename T>
void BinaryConstraint::findSupport(T getCost, EnumeratedVariable *x, EnumeratedVariable *y,
        vector<Value> &supportX, vector<StoreCost> &deltaCostsX)
{
	assert(connected());
    wcsp->revise(this);
    if (ToulBar2::verbose >= 3) cout << "findSupport C" << x->getName() << "," << y->getName() << endl;
    bool supportBroken = false;
    for (EnumeratedVariable::iterator iterX = x->begin(); iterX != x->end(); ++iterX) {
        unsigned int xindex = x->toIndex(*iterX);
        Value support = supportX[xindex];
        if (y->cannotbe(support) || getCost(x,y, *iterX, support) > MIN_COST) {
            Value minCostValue = y->getInf();
            Cost minCost = getCost(x,y, *iterX, minCostValue);
            EnumeratedVariable::iterator iterY = y->begin();
            for (++iterY; minCost > MIN_COST && iterY != y->end(); ++iterY) {
			  Cost cost = getCost(x,y, *iterX, *iterY);
                if (GLB(&minCost, cost)) {
                    minCostValue = *iterY;
                }
            }
            if (minCost > MIN_COST) {
                supportBroken |= project(x, *iterX, minCost, deltaCostsX);
                if (deconnected()) return;
            }
            supportX[xindex] = minCostValue;
        }
    }
    if (supportBroken) {
        x->findSupport();
    }
}

template <typename T>
void BinaryConstraint::findFullSupport(T getCost, EnumeratedVariable *x, EnumeratedVariable *y,
        vector<Value> &supportX, vector<StoreCost> &deltaCostsX,
        vector<Value> &supportY, vector<StoreCost> &deltaCostsY)
{
	assert(connected());
    wcsp->revise(this);
    if (ToulBar2::verbose >= 3) cout << "findFullSupport C" << x->getName() << "," << y->getName() << endl;
    bool supportBroken = false;
    for (EnumeratedVariable::iterator iterX = x->begin(); iterX != x->end(); ++iterX) {
    	unsigned int xindex = x->toIndex(*iterX);
        Value support = supportX[xindex];
        if (y->cannotbe(support) || getCost(x,y, *iterX, support) + y->getCost(support) > MIN_COST) {
            Value minCostValue = y->getInf();
            Cost minCost = getCost(x,y, *iterX, minCostValue) + y->getCost(minCostValue);
            EnumeratedVariable::iterator iterY = y->begin();
            for (++iterY; minCost > MIN_COST && iterY != y->end(); ++iterY) {
                Cost cost = getCost(x,y, *iterX, *iterY) + y->getCost(*iterY);
                if (GLB(&minCost, cost)) {
                    minCostValue = *iterY;
                }
            }
            if (minCost > MIN_COST) {
                // extend unary to binary
                for (EnumeratedVariable::iterator iterY = y->begin(); iterY != y->end(); ++iterY) {
                    Cost cost = getCost(x,y, *iterX, *iterY);
                    if (GLBTEST(minCost, cost)) {
						extend(y, *iterY, minCost - cost, deltaCostsY);
                        supportY[y->toIndex(*iterY)] = *iterX;
//                         if (ToulBar2::vac) {
//                             x->queueVAC2();
//                             y->queueVAC2();
//                         }
                     }
                }
                supportBroken |= project(x, *iterX, minCost, deltaCostsX);
                if (deconnected()) return;
            }
            supportX[xindex] = minCostValue;
        }
    }
    if (supportBroken) {
        x->findSupport();
    }
}

template <typename T>
void BinaryConstraint::projection(T getCost, EnumeratedVariable *x, EnumeratedVariable *y, Value valueY, vector<StoreCost> &deltaCostsX)
{
	x->queueDEE();
    bool supportBroken = false;
    wcsp->revise(this);
    for (EnumeratedVariable::iterator iterX = x->begin(); iterX != x->end(); ++iterX) {
        Cost cost = getCost(x,y, *iterX, valueY);
        if (cost > MIN_COST) {
            supportBroken |= project(x, *iterX, cost, deltaCostsX);
        }
    }
    if (supportBroken) {
        x->findSupport();
    }
}

template <typename T>
bool BinaryConstraint::verify(T getCost, EnumeratedVariable *x, EnumeratedVariable *y)
{
    for (EnumeratedVariable::iterator iterX = x->begin(); iterX != x->end(); ++iterX) {
        Cost minCost = getCost(x,y, *iterX, y->getInf());
		if (ToulBar2::LcLevel>=LC_DAC && getDACScopeIndex() == getIndex(x)) minCost += y->getCost(y->getInf());
        EnumeratedVariable::iterator iterY = y->begin();
        for (++iterY; minCost > MIN_COST && iterY != y->end(); ++iterY) {
            Cost cost = getCost(x,y, *iterX, *iterY);
			if (ToulBar2::LcLevel>=LC_DAC && getDACScopeIndex() == getIndex(x)) cost += y->getCost(*iterY);
            GLB(&minCost, cost);
        }
        if (minCost > MIN_COST) {
            cout << *this;
            return false;
        }
    }
    return true;
}


void BinaryConstraint::permute(EnumeratedVariable *xin, Value a, Value b)
{
	EnumeratedVariable *yin = y;
	if(xin != x) yin = x;

	vector<Cost> aux;
    for (EnumeratedVariable::iterator ity = yin->begin(); ity != yin->end(); ++ity)  aux.push_back( getCost(xin, yin, a, *ity) );
    for (EnumeratedVariable::iterator ity = yin->begin(); ity != yin->end(); ++ity)  setcost(xin, yin, a, *ity, getCost(xin, yin, b, *ity));

	vector<Cost>::iterator itc = aux.begin();
    for (EnumeratedVariable::iterator ity = yin->begin(); ity != yin->end(); ++ity) {
    	setcost(xin, yin, b, *ity, *itc);
    	++itc;
    }
}

bool BinaryConstraint::isFunctional(EnumeratedVariable* xin, EnumeratedVariable* yin, map<Value, Value> &functional)
{
  assert(xin != yin);
  assert(getIndex(xin) >= 0);
  assert(getIndex(yin) >= 0);
  bool isfunctional = true;
  functional.clear();
  for (EnumeratedVariable::iterator itx = xin->begin(); isfunctional && itx != xin->end(); ++itx) {
	bool first = true;
	for (EnumeratedVariable::iterator ity = yin->begin(); isfunctional && ity != yin->end(); ++ity) {
	  if (!CUT(getCost(xin,yin,*itx,*ity) + wcsp->getLb(), wcsp->getUb())) {
		if (first) {
		  functional[*itx] = *ity;
		  first = false;
		} else {
		  isfunctional = false;
		}
	  }
	}
	assert(!first); // assumes it is SAC already
  }
  if (isfunctional) return true;
  functional.clear();
  return false;
}
