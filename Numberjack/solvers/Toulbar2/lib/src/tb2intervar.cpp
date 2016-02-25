/*
 * ****** Variable with domain represented by an interval *******
 */

#include "tb2intervar.hpp"
#include "tb2wcsp.hpp"
#include "tb2clusters.hpp"


/*
 * Constructors and misc.
 * 
 */


IntervalVariable::IntervalVariable(WCSP *w, string n, Value iinf, Value isup) : Variable(w, n, iinf, isup),
        infCost(MIN_COST, &w->getStore()->storeCost), supCost(MIN_COST, &w->getStore()->storeCost)
{
}

void IntervalVariable::print(ostream& os)
{
    os << " [" << inf << "," << sup << "]";
    os << "/" << getDegree();
    if (ToulBar2::weightedDegree) os << "/" << getWeightedDegree();
    if (unassigned()) {
        os << " < " << getInfCost() << "," << getSupCost() << " >";
    }
}


/*
 * Propagation methods
 * 
 */

void IntervalVariable::projectInfCost(Cost cost)
{
    infCost += cost;
    assert(infCost >= MIN_COST);
    if (getInf() == maxCostValue || infCost > maxCost) queueNC();
    if (CUT(infCost + wcsp->getLb(),wcsp->getUb())) increaseFast(getInf() + 1);
}

void IntervalVariable::projectSupCost(Cost cost)
{
    supCost += cost;
    assert(supCost >= MIN_COST);
    if (getSup() == maxCostValue || supCost > maxCost) queueNC();
    if (CUT(supCost + wcsp->getLb(), wcsp->getUb())) decreaseFast(getSup() - 1);
}

void IntervalVariable::propagateNC()
{
    if (ToulBar2::verbose >= 3) cout << "propagateNC for " << getName() << endl;
    if (CUT(getInfCost() + wcsp->getLb(), wcsp->getUb())) increaseFast(getInf() + 1);
    if (CUT(getSupCost() + wcsp->getLb(), wcsp->getUb())) decreaseFast(getSup() - 1);
    if (getInfCost() > getSupCost()) {
        setMaxUnaryCost(getInf(), getInfCost());
    } else {
        setMaxUnaryCost(getSup(), getSupCost());
    }
}

bool IntervalVariable::verifyNC()
{
    if (CUT(getInfCost() + wcsp->getLb(),wcsp->getUb())) {
        cout << *this << " has inf cost not NC!" << endl;
        return false;
    }
    if (CUT(getSupCost() + wcsp->getLb(),wcsp->getUb())) {
        cout << *this << " has sup cost not NC!" << endl;
        return false;
    }
    return true;
}

void IntervalVariable::increaseFast(Value newInf)
{
    if (ToulBar2::verbose >= 2) cout << "increase " << getName() << " " << inf << " -> " << newInf << endl;
    assert(!wcsp->getIsPartOfOptimalSolution() || ((wcsp->getTreeDec())?wcsp->getTreeDec()->getRoot()->getUb():wcsp->getUb()) <= ToulBar2::verifiedOptimum || wcsp->getBestValue(wcspIndex) >= newInf);
    if (newInf > inf) {
        if (newInf > sup) {THROWCONTRADICTION;
        } else {
            if (newInf == sup) {assign(newInf);
            } else {
                inf = newInf;
                infCost = MIN_COST;
                queueInc();
                if (ToulBar2::setmin) (*ToulBar2::setmin)(wcsp->getIndex(), wcspIndex, newInf, wcsp->getSolver());
            }
        }
    }
}

void IntervalVariable::increase(Value newInf, bool isDecision)
{
    if (ToulBar2::verbose >= 2) cout << "increase " << getName() << " " << inf << " -> " << newInf << endl;
#ifndef NDEBUG
    if (isDecision && wcsp->getIsPartOfOptimalSolution() && wcsp->getBestValue(wcspIndex) < newInf) wcsp->setIsPartOfOptimalSolution(false);
    assert(isDecision || !wcsp->getIsPartOfOptimalSolution() || ((wcsp->getTreeDec())?wcsp->getTreeDec()->getRoot()->getUb():wcsp->getUb()) <= ToulBar2::verifiedOptimum || wcsp->getBestValue(wcspIndex) >= newInf);
#endif
    if (newInf > inf) {
        if (newInf > sup) {THROWCONTRADICTION;
        } else {
            if (newInf == sup) {assign(newInf);
            } else {
                inf = newInf;
                infCost = MIN_COST;
                if (newInf > maxCostValue) queueNC();           // single diff with increaseFast
                queueInc();
                if (ToulBar2::setmin) (*ToulBar2::setmin)(wcsp->getIndex(), wcspIndex, newInf, wcsp->getSolver());
            }
        }
    }
}

void IntervalVariable::decreaseFast(Value newSup)
{
    if (ToulBar2::verbose >= 2) cout << "decrease " << getName() << " " << sup << " -> " << newSup << endl;
    assert(!wcsp->getIsPartOfOptimalSolution() || ((wcsp->getTreeDec())?wcsp->getTreeDec()->getRoot()->getUb():wcsp->getUb()) <= ToulBar2::verifiedOptimum || wcsp->getBestValue(wcspIndex) <= newSup);
    if (newSup < sup) {
        if (newSup < inf) {THROWCONTRADICTION;
        } else {
            if (inf == newSup) {assign(newSup);
            } else {
                sup = newSup;
                supCost = MIN_COST;
                queueDec();
                if (ToulBar2::setmax) (*ToulBar2::setmax)(wcsp->getIndex(), wcspIndex, newSup, wcsp->getSolver());
            }
        }
    }
}

void IntervalVariable::decrease(Value newSup, bool isDecision)
{
    if (ToulBar2::verbose >= 2) cout << "decrease " << getName() << " " << sup << " -> " << newSup << endl;
#ifndef NDEBUG
    if (isDecision && wcsp->getIsPartOfOptimalSolution() && wcsp->getBestValue(wcspIndex) > newSup) wcsp->setIsPartOfOptimalSolution(false);
    assert(isDecision || !wcsp->getIsPartOfOptimalSolution() || ((wcsp->getTreeDec())?wcsp->getTreeDec()->getRoot()->getUb():wcsp->getUb()) <= ToulBar2::verifiedOptimum || wcsp->getBestValue(wcspIndex) <= newSup);
#endif
    if (newSup < sup) {
        if (newSup < inf) {THROWCONTRADICTION;
        } else {
            if (inf == newSup) {assign(newSup);
            } else {
                sup = newSup;
                supCost = MIN_COST;
                if (newSup < maxCostValue) queueNC();           // single diff with decreaseFast
                queueDec();
                if (ToulBar2::setmax) (*ToulBar2::setmax)(wcsp->getIndex(), wcspIndex, newSup, wcsp->getSolver());
            }
        }
    }
}

void IntervalVariable::assign(Value newValue, bool isDecision)
{
    if (ToulBar2::verbose >= 2) cout << "assign " << *this << " -> " << newValue << endl;
#ifndef NDEBUG
    if (isDecision && wcsp->getIsPartOfOptimalSolution() && wcsp->getBestValue(wcspIndex) != newValue) wcsp->setIsPartOfOptimalSolution(false);
    assert(isDecision || !wcsp->getIsPartOfOptimalSolution() || ((wcsp->getTreeDec())?wcsp->getTreeDec()->getRoot()->getUb():wcsp->getUb()) <= ToulBar2::verifiedOptimum || wcsp->getBestValue(wcspIndex) == newValue);
#endif
    if (unassigned() || getValue() != newValue) {
        if (cannotbe(newValue)) THROWCONTRADICTION;
        changeNCBucket(-1);
        maxCostValue = newValue;
        maxCost = MIN_COST;
        inf = newValue;
        sup = newValue;
        infCost = MIN_COST;
        supCost = MIN_COST;
        if (ToulBar2::setvalue) (*ToulBar2::setvalue)(wcsp->getIndex(), wcspIndex, newValue, wcsp->getSolver());
        for (ConstraintList::iterator iter=constrs.begin(); iter != constrs.end(); ++iter) {
            (*iter).constr->assign((*iter).scopeIndex);
        }
    }
}

/// assign a variable with delayed constraint propagation
void IntervalVariable::assignLS(Value newValue, set<Constraint *>& delayedCtrs)
{
    if (ToulBar2::verbose >= 2) cout << "assignLS " << *this << " -> " << newValue << endl;
    if (unassigned() || getValue() != newValue) {
        if (cannotbe(newValue)) THROWCONTRADICTION;
        changeNCBucket(-1);
        maxCostValue = newValue;
        maxCost = MIN_COST;
        inf = newValue;
        sup = newValue;
        infCost = MIN_COST;
        supCost = MIN_COST;
        if (ToulBar2::setvalue) (*ToulBar2::setvalue)(wcsp->getIndex(), wcspIndex, newValue, wcsp->getSolver());
        for (ConstraintList::iterator iter=constrs.begin(); iter != constrs.end(); ++iter) {
            delayedCtrs.insert((*iter).constr);
        }
    }
}

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

