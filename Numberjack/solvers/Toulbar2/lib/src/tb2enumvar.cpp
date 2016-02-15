
/*
 * ****** Variable with domain represented by an enumerated domain *******
 */

#include "tb2enumvar.hpp"
#include "tb2wcsp.hpp"
#include "tb2binconstr.hpp"
#include "tb2ternaryconstr.hpp"
#include "tb2clusters.hpp"

/*
 * Constructors and misc.
 *
 */


EnumeratedVariable::EnumeratedVariable(WCSP *w, string n, Value iinf, Value isup) : Variable(w, n, iinf, isup),
        domain(iinf, isup, &w->getStore()->storeDomain), deltaCost(MIN_COST, &w->getStore()->storeCost),
        support(iinf, &w->getStore()->storeValue), watchForIncrease(false), watchForDecrease(false)
{
    init();
}

EnumeratedVariable::EnumeratedVariable(WCSP *w, string n, Value *d, int dsize) : Variable(w, n, min(d,dsize), max(d, dsize)),
        domain(d, dsize, &w->getStore()->storeDomain), deltaCost(MIN_COST, &w->getStore()->storeCost),
        support(min(d,dsize), &w->getStore()->storeValue), watchForIncrease(false), watchForDecrease(false)
{
    init();
}

void EnumeratedVariable::init()
{
    if (wcsp->getStore()->getDepth() > 0) {
        cerr << "You cannot create a variable during the search!" << endl;
        exit(EXIT_FAILURE);
    }

    costs = vector<StoreCost>(getDomainInitSize(), StoreCost(MIN_COST, &wcsp->getStore()->storeCost));
    linkACQueue.content.var = this;
    linkACQueue.content.timeStamp = -1;
    linkDACQueue.content.var = this;
    linkDACQueue.content.timeStamp = -1;
    linkEAC1Queue.content.var = this;
    linkEAC1Queue.content.timeStamp = -1;
    linkEAC2Queue.content.var = this;
    linkEAC2Queue.content.timeStamp = -1;
    linkDEEQueue.content.var = this;
    linkDEEQueue.content.timeStamp = -1;
    DEE.constr = NULL;
    DEE.scopeIndex = -1;
    if (ToulBar2::DEE >= 2) DEE2 = vector<ConstraintLink>(getDomainInitSize()*getDomainInitSize(), DEE);
    queueDEE();
}

void EnumeratedVariable::getDomain(Value *array)
{
    for (iterator iter = begin(); iter != end(); ++iter) {
        *array = *iter;
        ++array;
    }
}

void EnumeratedVariable::getDomainAndCost(ValueCost *array)
{
    for (iterator iter = begin(); iter != end(); ++iter) {
        array->value = *iter;
        array->cost = getCost(*iter);
        ++array;
    }
}

Cost EnumeratedVariable::getBinaryCost(ConstraintLink c, Value myvalue, Value itsvalue) {
    return (c.scopeIndex == 0)?((BinaryConstraint *) c.constr)->getCost(myvalue, itsvalue):((BinaryConstraint *) c.constr)->getCost(itsvalue, myvalue);
}

Cost EnumeratedVariable::getBinaryCost(BinaryConstraint* c, Value myvalue, Value itsvalue) {
    return (c->getIndex(this) == 0)?c->getCost(myvalue, itsvalue):c->getCost(itsvalue, myvalue);
}


void EnumeratedVariable::print(ostream& os)
{
    if (unassigned()) {
        os << " " << domain;
    } else {
        os << " [" << inf << "," << sup << "]";
    }
    os << "/" << getDegree();
    if (ToulBar2::weightedDegree) os << "/" << getWeightedDegree();
    if (unassigned()) {
        os << " <";
        for (iterator iter=begin(); iter != end(); ++iter) {
            os << " " << getCost(*iter);
        }
        os << " > s:" << support;
    }
}

/*
 * Propagation methods
 *
 */

void EnumeratedVariable::queueAC()
{
    wcsp->queueAC(&linkACQueue);
}

void EnumeratedVariable::queueDAC()
{
    wcsp->queueDAC(&linkDACQueue);
}

void EnumeratedVariable::queueEAC1()
{
    wcsp->queueEAC1(&linkEAC1Queue);
}

void EnumeratedVariable::queueEAC2()
{
    wcsp->queueEAC2(&linkEAC2Queue);
}

void EnumeratedVariable::queueDEE()
{
    wcsp->queueDEE(&linkDEEQueue);
}

void EnumeratedVariable::project(Value value, Cost cost, bool delayed)
{
    assert(cost >= MIN_COST);
    Cost oldcost = getCost(value);
    costs[toIndex(value)] += cost;
    Cost newcost = oldcost + cost;
    if (value == maxCostValue || LUBTEST(maxCost, newcost)) queueNC();
    if (DACTEST(oldcost, cost)) {
        queueDAC();
        queueEAC1();
    }
    if (CUT(newcost + wcsp->getLb(), wcsp->getUb())) {
        if (delayed) queueNC();
        else removeFast(value);     // Avoid any unary cost overflow
    }
}

void EnumeratedVariable::projectInfCost(Cost cost)
{
    assert(cost >= MIN_COST);
    Value value = getInf();
    Cost oldcost = getCost(value);
    project(value, cost);
    if (support == value || SUPPORTTEST(oldcost, cost)) findSupport();
}

void EnumeratedVariable::projectSupCost(Cost cost)
{
    assert(cost >= MIN_COST);
    Value value = getSup();
    Cost oldcost = getCost(value);
    project(value, cost);
    if (support == value || SUPPORTTEST(oldcost, cost)) findSupport();
}

void EnumeratedVariable::extend(Value value, Cost cost)
{
    assert( ToulBar2::verbose < 4 || ((cout << "extend " << getName() << " (" << value << ") -= " << cost << endl), true) );
    assert(cost >= MIN_COST);
    assert(CUT(costs[toIndex(value)], cost));
    costs[toIndex(value)] -= cost;
    if (value == maxCostValue || PARTIALORDER) queueNC();
}

void EnumeratedVariable::extendAll(Cost cost)
{
    assert(cost > MIN_COST);
    deltaCost += cost;          // Warning! Possible overflow???
    queueNC();
}

void EnumeratedVariable::findSupport()
{
    if (cannotbe(support) || getCost(support) > MIN_COST) {
        Value newSupport = getInf();
        Cost minCost = getCost(newSupport);
        iterator iter = begin();
        for (++iter; minCost > MIN_COST && iter != end(); ++iter) {
            Cost cost = getCost(*iter);
            if (GLB(&minCost, cost)) {
                newSupport = *iter;
            }
        }
        if (minCost > MIN_COST) {
            extendAll(minCost);
            projectLB(minCost);
        }
        assert(canbe(newSupport) && (getCost(newSupport) == MIN_COST || SUPPORTTEST(getCost(newSupport))));
        if (support != newSupport) queueDEE();
        support = newSupport;
    }
}

void EnumeratedVariable::propagateNC()
{
    wcsp->revise(NULL);
    if (ToulBar2::verbose >= 3) cout << "propagateNC for " << getName() << endl;
    Value maxcostvalue = getSup()+1;
    Cost maxcost = MIN_COST;
    bool supportBroken = false;
    // Warning! the first value must be visited because it may be removed
    for (iterator iter = begin(); iter != end(); ++iter) {
        Cost cost = getCost(*iter);
        if (CUT(cost + wcsp->getLb(), wcsp->getUb())) {
            if (SUPPORTTEST(cost)) supportBroken = true;
            removeFast(*iter);
        } else if (LUB(&maxcost, cost) || cannotbe(maxcostvalue)) {
            maxcostvalue = *iter;
        }
    }
    assert(getCost(maxcostvalue) == maxcost || !LUBTEST(maxcost, getCost(maxcostvalue)));
    setMaxUnaryCost(maxcostvalue, maxcost);
    if (supportBroken) findSupport();
}

bool EnumeratedVariable::verifyNC()
{
    bool supported = true;
    Cost minCost = MAX_COST;
    //    Value minCostValue = getSup()+1;
    for (iterator iter = begin(); iter != end(); ++iter) {
        Cost cost = getCost(*iter);
        if (CUT(cost + wcsp->getLb(), wcsp->getUb())) {
            cout << *this << " not NC!" << endl;
            return false;
        }
        //		if (cost < minCost) minCostValue = *iter;
        GLB(&minCost, cost);
    }
    if (minCost > MIN_COST) {
        cout << *this << " not NC*!" << endl;
        supported = false;
    }
    if (cannotbe(support) || (getCost(support)>MIN_COST && !SUPPORTTEST(getCost(support)))) { // || minCostValue != support) {
        cout << *this << " has an unvalid NC support!" << endl;
        supported = false;
    }
    return supported;
}

void EnumeratedVariable::propagateAC()
{
    for (ConstraintList::iterator iter=constrs.begin(); iter != constrs.end(); ++iter) {
        (*iter).constr->remove((*iter).scopeIndex);
    }
}

void EnumeratedVariable::propagateDAC()
{
    for (ConstraintList::iterator iter=constrs.rbegin(); iter != constrs.rend(); --iter) {
        (*iter).constr->projectFromZero((*iter).scopeIndex);
    }
}

void EnumeratedVariable::fillEAC2(bool self)
{
    if (self) queueEAC2();
    for (ConstraintList::iterator iter=constrs.begin(); iter != constrs.end(); ++iter) {
        if ((*iter).constr->isDuplicate()) continue;
        (*iter).constr->fillEAC2((*iter).scopeIndex);
    }
}

void EnumeratedVariable::setCostProvidingPartition() 
{	
    vector<bool> used(wcsp->numberOfVariables(), false);
    //  set<int> used;
    int maxArity = 4;
    // binary and ternary constraints are not under consideration
    for (ConstraintList::iterator iter = constrs.begin(); iter != constrs.end(); ++iter) {
        int arity = (*iter).constr->arity();
        maxArity = max(maxArity, arity);
        if (arity <= 3) {
            for (int i=0;i<arity;i++) {
                Variable *var = (*iter).constr->getVar(i);
                if (var != this) {
                    used[var->wcspIndex] = true;
                    //		  used.insert(var->wcspIndex);
                }
            }
        }
    }

    for (int scopeSize = maxArity;scopeSize>3;scopeSize--) {
        for (ConstraintList::iterator iter = constrs.begin(); iter != constrs.end(); ++iter) {
            int arity = (*iter).constr->arity();
            if (arity == scopeSize) {
                for (int i=0;i<arity;i++) {
                    Variable *var = (*iter).constr->getVar(i);
                    if (var != this && !used[var->wcspIndex]) {
                        //		  if (var != this && (used.find(var->wcspIndex) == used.end())) {
                        used[var->wcspIndex] = true;
                        //			used.insert(var->wcspIndex);
                        (*iter).constr->linkCostProvidingPartition((*iter).scopeIndex, var);
                    }
                }
            }
        }
    }
}

bool EnumeratedVariable::isEAC(Value a)
{
    if (getCost(a)==MIN_COST) {
        for (ConstraintList::iterator iter=constrs.begin(); iter != constrs.end(); ++iter) {
            if ((*iter).constr->isDuplicate()) continue;

            if (!(*iter).constr->isEAC((*iter).scopeIndex, a)) {
#ifndef NDEBUG
                if (ToulBar2::verbose >=4) {
                    cout << getName() << "(" << a << ") is not EAC due to constraint " << *(*iter).constr << endl;
                    if ((*iter).constr->arity()==3) {
                        TernaryConstraint *c = (TernaryConstraint *) (*iter).constr;
                        if (c->xy->connected()) cout << *c->xy;
                        if (c->xz->connected()) cout << *c->xz;
                        if (c->yz->connected()) cout << *c->yz;
                    }
                }
#endif
                return false;
            }
        }
        if (support != a) queueDEE();
        support = a;
#ifndef NDEBUG
        if (ToulBar2::verbose >=4) cout << getName() << "(" << a << ") is EAC!" << endl;
#endif
        return true;
    }
#ifndef NDEBUG
    if (ToulBar2::verbose >=4) cout << getName() << "(" << a << ") is not EAC due to unary cost " << getCost(a) << endl;
#endif
    return false;
}

bool EnumeratedVariable::isEAC()
{
    assert(canbe(support));
    if (isEAC(support)) return true;
    else {
        for (iterator iter = begin(); iter != end(); ++iter) {
            if (*iter != support && isEAC(*iter)) return true;
        }
    }
    return false;
}

void EnumeratedVariable::propagateEAC()
{
    if (!isEAC()) {
#ifndef NDEBUG
        Cost beforeLb = wcsp->getLb();
        assert(verify());
#endif
        for (ConstraintList::iterator iter = constrs.begin(); iter != constrs.end(); ++iter) {
            if ((*iter).constr->isDuplicate()) continue;
            (*iter).constr->findFullSupportEAC((*iter).scopeIndex);
        }
        fillEAC2(false);
        if (unassigned()) {
            queueDEE();
            // findFullSupportEAC may have inserted current variable in EAC1
            if (!linkEAC1Queue.removed) {
                if (linkEAC1Queue.content.timeStamp == wcsp->getNbNodes()) {
                    assert(((BTList<VariableWithTimeStamp> *) wcsp->getQueueEAC1())->inBTList(&linkEAC1Queue));
                    wcsp->getQueueEAC1()->remove(&linkEAC1Queue);
                } else {
                    linkEAC1Queue.removed = true;
                    linkEAC1Queue.content.timeStamp = -1;
                    linkEAC1Queue.content.incdec = NOTHING_EVENT;
                }
            }
#ifndef NDEBUG
            // check if lb has effectively been increased
            if (wcsp->getLb() == beforeLb)
                if (ToulBar2::verbose >= 1) cout << "EAC failed on " << getName() << endl;
#endif
        }
    }
}

void EnumeratedVariable::propagateDEE(Value a, Value b, bool dee)
{   //cout << "check DEE for " << *this << " " << a << " -> " << b << " " << dee << endl;
    if (a == b) return;
    Cost costa = getCost(a);
    Cost costb = getCost(b);
    assert(!dee || costa <= costb);
    Cost totalmaxcosta = costa;
    Cost totaldiffcosta = costa;
    Cost totalmaxcostb = costb;
    Cost totaldiffcostb = costb;
    ConstraintLink residue = ((dee)?DEE:DEE2[a * getDomainInitSize() + b]);
    ConstraintLink residue2 = ((dee)?DEE:DEE2[b * getDomainInitSize() + a]);
    if (costa <= costb && residue.constr && residue.constr->connected() && residue.scopeIndex < residue.constr->arity() && residue.constr->getVar(residue.scopeIndex) == this) {
        pair< pair<Cost,Cost>, pair<Cost,Cost> > costs = residue.constr->getMaxCost(residue.scopeIndex, a, b);
        if (totalmaxcosta <= getMaxCost()) totalmaxcosta += costs.first.first;
        if (totalmaxcostb <= getMaxCost()) totalmaxcostb += costs.second.first;
        if (totaldiffcosta <= getMaxCost()) totaldiffcosta += costs.first.second;
        if (totaldiffcostb <= getMaxCost()) totaldiffcostb += costs.second.second;
        if (totaldiffcosta > costb && totaldiffcostb > costa) return;
    }
    if (costb <= costa && residue2.constr && (residue2.constr != residue.constr || costa > costb) && residue2.constr->connected() && residue2.scopeIndex < residue2.constr->arity() && residue2.constr->getVar(residue2.scopeIndex) == this) {
        pair< pair<Cost,Cost>, pair<Cost,Cost> > costs = residue2.constr->getMaxCost(residue2.scopeIndex, a, b);
        if (totalmaxcosta <= getMaxCost()) totalmaxcosta += costs.first.first;
        if (totalmaxcostb <= getMaxCost()) totalmaxcostb += costs.second.first;
        if (totaldiffcosta <= getMaxCost()) totaldiffcosta += costs.first.second;
        if (totaldiffcostb <= getMaxCost()) totaldiffcostb += costs.second.second;
        if (totaldiffcosta > costb && totaldiffcostb > costa) return;
    }
    for (ConstraintList::iterator iter=constrs.begin(); iter != constrs.end(); ++iter) {
        if (costa <= costb && residue.constr == (*iter).constr && residue.scopeIndex == (*iter).scopeIndex) continue;
        if (costb <= costa && residue2.constr == (*iter).constr && residue2.scopeIndex == (*iter).scopeIndex) continue;
        pair< pair<Cost,Cost>, pair<Cost,Cost> > costs = (*iter).constr->getMaxCost((*iter).scopeIndex, a, b);
        if (costs.second.second > costa) {
            if (dee) {
                if (costa == costb) DEE = (*iter);
            } else DEE2[b * getDomainInitSize() + a] = (*iter);
        }
        if (costs.first.second > costb) {
            if (dee) DEE = (*iter);
            else DEE2[a * getDomainInitSize() + b] = (*iter);
        }
        if (totalmaxcosta <= getMaxCost()) totalmaxcosta += costs.first.first;
        if (totalmaxcostb <= getMaxCost()) totalmaxcostb += costs.second.first;
        if (totaldiffcosta <= getMaxCost()) totaldiffcosta += costs.first.second;
        if (totaldiffcostb <= getMaxCost()) totaldiffcostb += costs.second.second;
        if (totaldiffcosta > costb && totaldiffcostb > costa) return;
    }
    assert((totalmaxcosta >= totaldiffcosta || (totalmaxcosta > getMaxCost() && totaldiffcosta > getMaxCost())));
    assert((totalmaxcostb >= totaldiffcostb || (totalmaxcostb > getMaxCost() && totaldiffcostb > getMaxCost())));
    assert(totaldiffcosta <= costb || totaldiffcostb <= costa);
    if (totalmaxcosta == MIN_COST) {
        if (ToulBar2::verbose >= 2) cout << "DEE " << *this << " (" << a << "," << totalmaxcosta << ") -> (*,*)" << endl;
        wcsp->incNbDEE(getDomainSize()-1);
        assign(a);
        return;
    }
    if (totalmaxcostb == MIN_COST) {
        if (ToulBar2::verbose >= 2) cout << "DEE " << *this << " (" << b << "," << totalmaxcostb << ") -> (*,*)" << endl;
        wcsp->incNbDEE(getDomainSize()-1);
        assign(b);
        return;
    }
    if (totaldiffcosta <= costb) {
        if (ToulBar2::verbose >= 2) cout << "DEE " << *this << " (" << a << "," << totaldiffcosta << ") -> (" << b << "," << costb << ")" << endl;
        wcsp->incNbDEE();
        remove(b);
        if (assigned()) return;
    } else {
        if (ToulBar2::verbose >= 2) cout << "DEE " << *this << " (" << b << "," << totaldiffcostb << ") -> (" << a << "," << costa << ")" << endl;
        wcsp->incNbDEE();
        remove(a);
        if (assigned()) return;
    }
    assert(unassigned());
    Cost totalmaxcost = min(totalmaxcosta, totalmaxcostb);
    if (totalmaxcost > getMaxCost()) return;
    for (iterator iter = begin(); iter != end(); ++iter) {
        if (*iter == a) continue;
        if (*iter == b) continue;
        assert(unassigned());
        if (getCost(*iter) >= totalmaxcost) {
            if (ToulBar2::verbose >= 2) cout << "DEE " << *this << " (" << ((totalmaxcosta<totalmaxcostb)?a:b) << "," << totalmaxcost << ") -> (" << *iter << "," << getCost(*iter) << ")" << endl;
            wcsp->incNbDEE();
            remove(*iter);
        }
    }
}

bool EnumeratedVariable::verifyDEE(Value a, Value b)
{
    if (a == b) return true;
    Cost totalmaxcost = getCost(a);
    Cost totaldiffcost = getCost(a);
    for (ConstraintList::iterator iter=constrs.begin(); iter != constrs.end(); ++iter) {
        pair< pair<Cost,Cost>, pair<Cost,Cost> > costs = (*iter).constr->getMaxCost((*iter).scopeIndex, a, b);
        if (totalmaxcost + wcsp->getLb() < wcsp->getUb()) totalmaxcost += costs.first.first;
        if (totaldiffcost + wcsp->getLb() < wcsp->getUb()) totaldiffcost += costs.first.second;
    }
    if (getCost(b) >= ((ToulBar2::DEE_>=3 || (ToulBar2::DEE_==2 && wcsp->getStore()->getDepth()==0))?totaldiffcost:totalmaxcost)) {
        cout << *this << " has missed dominated value (" << a << "," << ((ToulBar2::DEE_>=3 || (ToulBar2::DEE_==2 && wcsp->getStore()->getDepth()==0))?totaldiffcost:totalmaxcost) << ") -> (" << b << "," << getCost(b) << ")"  << endl;
        return true; // should be false but we need to queue all variables each time LB or UB change
    }
    for (iterator iter = begin(); iter != end(); ++iter) {
        if (*iter == a) continue;
        if (*iter == b) continue;
        if (getCost(*iter) >= totalmaxcost) {
            cout << *this << " has missed dominated value (" << a << "," << totalmaxcost << ") -> (" << (*iter) << "," << getCost(*iter) << ")"  << endl;
            return true; // should be false but we need to queue all variables each time LB or UB change
        }
    }
    return true;
}

bool EnumeratedVariable::verifyDEE()
{
    if (ToulBar2::DEE_>=3 || (ToulBar2::DEE_==2 && wcsp->getStore()->getDepth()==0)) {
        for (iterator itera = begin(); itera != end(); ++itera) {
            for (iterator iterb = begin(); iterb != end(); ++iterb) {
                if (!verifyDEE(*itera, *iterb)) return false;
            }
        }
    } else {
        if (getSupport() == getMaxCostValue()) {
            if (getSupport() != getSup()) {
                if (!verifyDEE(getSupport(), getSup())) return false;
            } else {
                if (!verifyDEE(getSupport(), getInf())) return false;
            }
        } else if (!verifyDEE(getSupport(), getMaxCostValue())) return false;
    }
    return true;
}

void EnumeratedVariable::increaseFast(Value newInf)
{
    if (ToulBar2::verbose >= 2) cout << "increase " << getName() << " " << inf << " -> " << newInf << endl;
    assert(!wcsp->getIsPartOfOptimalSolution() || ((wcsp->getTreeDec())?wcsp->getTreeDec()->getRoot()->getUb():wcsp->getUb()) <= ToulBar2::verifiedOptimum || wcsp->getBestValue(wcspIndex) >= newInf);
    if (newInf > inf) {
        if (newInf > sup) {THROWCONTRADICTION;
        } else {
            newInf = domain.increase(newInf);
            if (newInf == sup) {assign(newInf);
            } else {
                inf = newInf;
                if (watchForIncrease) queueInc(); else queueAC();
                if (PARTIALORDER) queueDAC();
                if (ToulBar2::setmin) (*ToulBar2::setmin)(wcsp->getIndex(), wcspIndex, newInf, wcsp->getSolver());
            }
        }
    }
}

void EnumeratedVariable::increase(Value newInf, bool isDecision)
{
    if (ToulBar2::verbose >= 2) cout << "increase " << getName() << " " << inf << " -> " << newInf << endl;
#ifndef NDEBUG
    if (isDecision && wcsp->getIsPartOfOptimalSolution() && wcsp->getBestValue(wcspIndex) < newInf) wcsp->setIsPartOfOptimalSolution(false);
    assert(isDecision || !wcsp->getIsPartOfOptimalSolution() || ((wcsp->getTreeDec())?wcsp->getTreeDec()->getRoot()->getUb():wcsp->getUb()) <= ToulBar2::verifiedOptimum || wcsp->getBestValue(wcspIndex) >= newInf);
#endif
    if (newInf > inf) {
        if (newInf > sup) {THROWCONTRADICTION;
        } else {
            newInf = domain.increase(newInf);
            if (newInf == sup) {assign(newInf);
            } else {
                inf = newInf;
                if (newInf > maxCostValue || PARTIALORDER) queueNC();           // diff with increaseFast
                if (newInf > support || PARTIALORDER) findSupport();            // diff with increaseFast
                queueDAC();                                     // diff with increaseFast
                queueEAC1();                                     // diff with increaseFast
                if (watchForIncrease) queueInc(); else queueAC();
                if (ToulBar2::setmin) (*ToulBar2::setmin)(wcsp->getIndex(), wcspIndex, newInf, wcsp->getSolver());
            }
        }
    }
}

void EnumeratedVariable::decreaseFast(Value newSup)
{
    if (ToulBar2::verbose >= 2) cout << "decrease " << getName() << " " << sup << " -> " << newSup << endl;
    assert(!wcsp->getIsPartOfOptimalSolution() || ((wcsp->getTreeDec())?wcsp->getTreeDec()->getRoot()->getUb():wcsp->getUb()) <= ToulBar2::verifiedOptimum || wcsp->getBestValue(wcspIndex) <= newSup);
    if (newSup < sup) {
        if (newSup < inf) {THROWCONTRADICTION;
        } else {
            newSup = domain.decrease(newSup);
            if (inf == newSup) {assign(newSup);
            } else {
                sup = newSup;
                if (watchForDecrease) queueDec(); else queueAC();
                if (PARTIALORDER) queueDAC();
                if (ToulBar2::setmax) (*ToulBar2::setmax)(wcsp->getIndex(), wcspIndex, newSup, wcsp->getSolver());
            }
        }
    }
}

void EnumeratedVariable::decrease(Value newSup, bool isDecision)
{
    if (ToulBar2::verbose >= 2) cout << "decrease " << getName() << " " << sup << " -> " << newSup << endl;
#ifndef NDEBUG
    if (isDecision && wcsp->getIsPartOfOptimalSolution() && wcsp->getBestValue(wcspIndex) > newSup) wcsp->setIsPartOfOptimalSolution(false);
    assert(isDecision || !wcsp->getIsPartOfOptimalSolution() || ((wcsp->getTreeDec())?wcsp->getTreeDec()->getRoot()->getUb():wcsp->getUb()) <= ToulBar2::verifiedOptimum || wcsp->getBestValue(wcspIndex) <= newSup);
#endif
    if (newSup < sup) {
        if (newSup < inf) {THROWCONTRADICTION;
        } else {
            newSup = domain.decrease(newSup);
            if (inf == newSup) {assign(newSup);
            } else {
                sup = newSup;
                if (newSup < maxCostValue || PARTIALORDER) queueNC();           // diff with decreaseFast
                if (newSup < support || PARTIALORDER) findSupport();            // diff with decreaseFast
                queueDAC();                                     // diff with decreaseFast
                queueEAC1();                                     // diff with decreaseFast
                if (watchForDecrease) queueDec(); else queueAC();
                if (ToulBar2::setmax) (*ToulBar2::setmax)(wcsp->getIndex(), wcspIndex, newSup, wcsp->getSolver());
            }
        }
    }
}

void EnumeratedVariable::removeFast(Value value)
{
    if (ToulBar2::verbose >= 2) cout << "remove " << *this << " <> " << value << endl;
    assert(!wcsp->getIsPartOfOptimalSolution() || ((wcsp->getTreeDec())?wcsp->getTreeDec()->getRoot()->getUb():wcsp->getUb()) <= ToulBar2::verifiedOptimum || wcsp->getBestValue(wcspIndex) != value);
    if (value == inf) increaseFast(value + 1);
    else if (value == sup) decreaseFast(value - 1);
    else if (canbe(value)) {
        domain.erase(value);
        queueAC();
        if (PARTIALORDER) queueDAC();
        if (ToulBar2::removevalue) (*ToulBar2::removevalue)(wcsp->getIndex(), wcspIndex, value, wcsp->getSolver());
    }
}

void EnumeratedVariable::remove(Value value, bool isDecision)
{
    if (ToulBar2::verbose >= 2) cout << "remove " << *this << " <> " << value << endl;
#ifndef NDEBUG
    if (isDecision && wcsp->getIsPartOfOptimalSolution() && wcsp->getBestValue(wcspIndex) == value) wcsp->setIsPartOfOptimalSolution(false);
    assert(isDecision || !wcsp->getIsPartOfOptimalSolution() || ((wcsp->getTreeDec())?wcsp->getTreeDec()->getRoot()->getUb():wcsp->getUb()) <= ToulBar2::verifiedOptimum || wcsp->getBestValue(wcspIndex) != value);
#endif
    if (value == inf) increase(value + 1);
    else if (value == sup) decrease(value - 1);
    else if (canbe(value)) {
        domain.erase(value);
        if (value == maxCostValue || PARTIALORDER) queueNC();
        if (value == support || PARTIALORDER) findSupport();
        queueDAC();
        queueEAC1();
        queueAC();
        if (ToulBar2::removevalue) (*ToulBar2::removevalue)(wcsp->getIndex(), wcspIndex, value, wcsp->getSolver());
    }
}

// this function is used ONLY for restoring the solution when
// variable elimination is tuned on
void EnumeratedVariable::assignWhenEliminated(Value newValue)
{
    //          assert(NCBucket == -1); // may be not already assigned in BTD newsolution
    inf = newValue;
    sup = newValue;
    support = newValue;
    maxCostValue = newValue;
    maxCost = MIN_COST;
}

void EnumeratedVariable::assign(Value newValue, bool isDecision)
{
    if (ToulBar2::verbose >= 2) cout << "assign " << *this << " -> " << newValue << endl;
#ifndef NDEBUG
    if (isDecision && wcsp->getIsPartOfOptimalSolution() && wcsp->getBestValue(wcspIndex) != newValue) wcsp->setIsPartOfOptimalSolution(false);
    assert(isDecision || !wcsp->getIsPartOfOptimalSolution() || ((wcsp->getTreeDec())?wcsp->getTreeDec()->getRoot()->getUb():wcsp->getUb()) <= ToulBar2::verifiedOptimum || wcsp->getBestValue(wcspIndex) == newValue);
#endif
    if (unassigned() || getValue() != newValue) {
        if (cannotbe(newValue)) THROWCONTRADICTION;
        changeNCBucket(-1);
        inf = newValue;
        sup = newValue;
        support = newValue;
        maxCostValue = newValue;
        maxCost = MIN_COST;

        Cost cost = getCost(newValue);
        if (cost > MIN_COST) {
            deltaCost += cost;
            projectLB(cost);
        }

        if (ToulBar2::setvalue) (*ToulBar2::setvalue)(wcsp->getIndex(), wcspIndex, newValue, wcsp->getSolver());
        for (ConstraintList::iterator iter=constrs.begin(); iter != constrs.end(); ++iter) {
            (*iter).constr->assign((*iter).scopeIndex);
        }
        //        for (ConstraintList::iterator iter=triangles.begin(); iter != triangles.end(); ++iter) {
        //            (*iter).constr->deconnect();
        //        }
    }
}

void EnumeratedVariable::assignLS(Value newValue, set<Constraint *>& delayedCtrs)
{
    if (ToulBar2::verbose >= 2) cout << "assignLS " << *this << " -> " << newValue << endl;
#ifndef NDEBUG
    if (wcsp->getIsPartOfOptimalSolution() && wcsp->getBestValue(wcspIndex) != newValue) wcsp->setIsPartOfOptimalSolution(false);
#endif
    if (unassigned() || getValue() != newValue) {
        if (cannotbe(newValue)) THROWCONTRADICTION;
        changeNCBucket(-1);
        inf = newValue;
        sup = newValue;
        support = newValue;
        maxCostValue = newValue;
        maxCost = MIN_COST;

        Cost cost = getCost(newValue);
        if (cost > MIN_COST) {
            deltaCost += cost;
            projectLB(cost);
        }

        if (ToulBar2::setvalue) (*ToulBar2::setvalue)(wcsp->getIndex(), wcspIndex, newValue, wcsp->getSolver());
        for (ConstraintList::iterator iter=constrs.begin(); iter != constrs.end(); ++iter) {
            delayedCtrs.insert((*iter).constr);
        }
    }
}

// eliminates the current (this) variable that participates
// in a single binary constraint ctr
bool EnumeratedVariable::elimVar( BinaryConstraint* ctr )
{
    assert(getDegree() == 1);

    EnumeratedVariable *x = (EnumeratedVariable *) ctr->getVarDiffFrom(this);

    TreeDecomposition* td = wcsp->getTreeDec();
    if(td && cluster != ctr->getCluster()) return false;

    if (ToulBar2::verbose >= 2) cout << "   elim linked to one binary " << ctr << endl;

    // deconnect first to be sure the current var is not involved in future propagation
    ctr->deconnect();
    // to be done before propagation
    WCSP::elimInfo ei = {this, x, NULL, ctr, NULL, NULL, NULL};
    wcsp->elimInfos[wcsp->getElimOrder()] = ei;
    wcsp->elimOrderInc();

    Cost mycosts[x->getDomainInitSize()];
#ifndef NDEBUG
    for (unsigned int i=0; i<x->getDomainInitSize(); i++) mycosts[i]=-MAX_COST;
#endif
    Cost negcost = MIN_COST;
    bool supportBroken = false;
    for (iterator iter1 = x->begin(); iter1 != x->end(); ++iter1) {
        Cost mincost = MAX_COST;
        for (iterator iter = begin(); iter != end(); ++iter) {
            Cost curcost = getCost(*iter) + getBinaryCost(ctr, *iter, *iter1);
            if (ToulBar2::isZ) mincost = wcsp->LogSumExp(mincost,curcost);
            else if (curcost < mincost) mincost = curcost;
        }
        if (ToulBar2::isZ) {
            if (mincost < negcost) negcost = mincost;
            mycosts[x->toIndex(*iter1)] = mincost;
        } else if (mincost > MIN_COST) {
            if(td) td->addDelta(cluster,x,*iter1,mincost);
            if (x->getSupport() == *iter1) supportBroken = true;
            x->project(*iter1, mincost);
        }
    }
    assert(negcost <= MIN_COST);
    if (ToulBar2::isZ) {
        for (iterator iter1 = x->begin(); iter1 != x->end(); ++iter1) {
            assert(mycosts[x->toIndex(*iter1)] != -MAX_COST);
            Cost mincost = mycosts[x->toIndex(*iter1)] - negcost;
            assert(mincost >= MIN_COST);
            if (mincost > MIN_COST) {
                if(td) td->addDelta(cluster,x,*iter1,mincost);
                if (x->getSupport() == *iter1) supportBroken = true;
                x->project(*iter1, mincost);
            }
        }
        if (negcost < 0) wcsp->decreaseLb(negcost);
    }
    if (supportBroken) {  x->findSupport();  }

    return true;
}


// eliminates the current (this) variable that participates
// in two binary constraints (its links are xylink and xzlink)
bool EnumeratedVariable::elimVar( ConstraintLink  xylink,  ConstraintLink xzlink )
{
    EnumeratedVariable *y = (EnumeratedVariable *) wcsp->getVar(xylink.constr->getSmallestVarIndexInScope(xylink.scopeIndex));
    EnumeratedVariable *z = (EnumeratedVariable *) wcsp->getVar(xzlink.constr->getSmallestVarIndexInScope(xzlink.scopeIndex));

    TreeDecomposition* td = wcsp->getTreeDec();
    if(td) {
        if(y->isSep() &&  z->isSep()) return false;
        if((cluster != xylink.constr->getCluster()) ||
                (cluster != xzlink.constr->getCluster()) ||
                (xylink.constr->getCluster() != xzlink.constr->getCluster()))
            return false;
    }

    assert(getDegree() == 2);
    xylink.constr->deconnect();
    xzlink.constr->deconnect();

    BinaryConstraint* yz = y->getConstr(z);

    if(td && yz && (getCluster() != yz->getCluster())) {
        BinaryConstraint* yz_ =  y->getConstr(z, getCluster());
        if(yz_) {
            yz = yz_;
            if (ToulBar2::verbose > 1) cout << "(" << y->wcspIndex << "," << z->wcspIndex << ") var elim binary is already duplicate and has same cluster" << endl;
        }
    }

    BinaryConstraint* yznew = wcsp->newBinaryConstr(y,z,xylink.constr,xzlink.constr);
    wcsp->elimBinOrderInc();

    Cost negcost = MIN_COST;
    for (iterator itery = y->begin(); itery != y->end(); ++itery) {
        for (iterator iterz = z->begin(); iterz != z->end(); ++iterz) {
            Cost mincost = MAX_COST;

            for (iterator iter = begin(); iter != end(); ++iter) {
                Cost curcost = getCost(*iter) +
                        getBinaryCost(xylink, *iter, *itery) +
                        getBinaryCost(xzlink, *iter, *iterz);

                if (ToulBar2::isZ) mincost = wcsp->LogSumExp(mincost,curcost);
                else if (curcost < mincost) mincost = curcost;
            }
            if (mincost < negcost) negcost = mincost;
            yznew->setcost(*itery,*iterz,mincost); // Warning! it can set a negative cost temporally
        }}
    assert(negcost <= MIN_COST);
    if (negcost < MIN_COST) {
        for (iterator itery = y->begin(); itery != y->end(); ++itery) {
            for (iterator iterz = z->begin(); iterz != z->end(); ++iterz) {
                yznew->addcost(*itery,*iterz, -negcost);
            }
        }
        wcsp->decreaseLb(negcost);
    }

    if(yz) {
        if(td && yz->getCluster() != cluster) {
            yz = yznew;
            yz->reconnect();
            yz->setDuplicate();
        } else {
            yz->addCosts( yznew );
            if (y->unassigned() && z->unassigned()) yz->reconnect();
        }
    } else {
        yz = yznew;
        yz->reconnect();
    }

    if(td) yz->setCluster( getCluster() );
    // to be done before propagation
    WCSP::elimInfo ei = {this, y,z, (BinaryConstraint*) xylink.constr, (BinaryConstraint*) xzlink.constr, NULL, NULL};
    wcsp->elimInfos[wcsp->getElimOrder()] = ei;
    wcsp->elimOrderInc();
    yz->propagate();
    return true;
}

// eliminates the current (this) variable that participates
// in the ternary constraint 'xyz'
// the function can fail to eliminate the current variable
// if it is linked to more than (in total) two variables.
// It returns true if the current variable was eliminated
bool EnumeratedVariable::elimVar( TernaryConstraint* xyz )
{
    BinaryConstraint* yz = NULL;
    if(xyz->xy->getIndex(this) < 0) yz = xyz->xy;
    else if(xyz->xz->getIndex(this) < 0) yz = xyz->xz;
    else if(xyz->yz->getIndex(this) < 0) yz = xyz->yz;
    assert(yz != NULL);

    int n2links = 0;
    int n3links = 0;

    ConstraintLink links[2] = {{NULL, 0},{NULL, 0}};
    for(ConstraintList::iterator iter=constrs.begin(); iter != constrs.end(); ++iter) {
        if((*iter).constr->arity() == 2) links[n2links++] =  (*iter);
        else if((*iter).constr->arity() == 3) n3links++;
        else return false;
    }

    if(n3links > 1) return false;

    if(wcsp->getTreeDec()) {
        if(cluster != xyz->getCluster()) return false;
    }

    for(int i=0; i<n2links; i++) {
        int idvar = links[i].constr->getSmallestVarIndexInScope(links[i].scopeIndex);
        if(xyz->getIndex( wcsp->getVar(idvar) ) < 0) return false;
    }

    xyz->deconnect();
    if(n2links > 0) links[0].constr->deconnect();
    if(n2links > 1) links[1].constr->deconnect();

    EnumeratedVariable* y = (EnumeratedVariable*) yz->getVar(0);
    EnumeratedVariable* z = (EnumeratedVariable*) yz->getVar(1);

    bool flag_rev = false;
    if(n2links > 0) {
        if (links[0].constr->getSmallestVarIndexInScope(links[0].scopeIndex) != y->wcspIndex) {
            assert(links[0].constr->getSmallestVarIndexInScope(links[0].scopeIndex) == z->wcspIndex);
            flag_rev = true;
        }
    }

    Cost negcost = MIN_COST;
    for (iterator itery = y->begin(); itery != y->end(); ++itery) {
        for (iterator iterz = z->begin(); iterz != z->end(); ++iterz) {
            Cost mincost = MAX_COST;
            for (iterator iter = begin(); iter != end(); ++iter) {
                Cost curcost = getCost(*iter) + xyz->getCost(this,y,z,*iter,*itery,*iterz);

                if(!flag_rev) {
                    if(n2links > 0) { assert(links[0].constr->getIndex(y) >= 0);
                    curcost += getBinaryCost(links[0], *iter, *itery); }
                    if(n2links > 1) { assert(links[1].constr->getIndex(z) >= 0);
                    curcost += getBinaryCost(links[1], *iter, *iterz); }
                } else {
                    if(n2links > 0) { assert(links[0].constr->getIndex(z) >= 0);
                    curcost += getBinaryCost(links[0], *iter, *iterz); }
                    if(n2links > 1) { assert(links[1].constr->getIndex(y) >= 0);
                    curcost += getBinaryCost(links[1], *iter, *itery); }
                }
                if (ToulBar2::isZ) mincost = wcsp->LogSumExp(mincost,curcost);
                else if (curcost < mincost) mincost = curcost;
            }
            if (mincost < negcost) negcost = mincost;
            yz->addcost(*itery,*iterz,mincost); // Warning! it can add a negative cost temporally
        }
    }
    assert(negcost <= MIN_COST);
    if (negcost < MIN_COST) {
        for (iterator itery = y->begin(); itery != y->end(); ++itery) {
            for (iterator iterz = z->begin(); iterz != z->end(); ++iterz) {
                yz->addcost(*itery,*iterz, -negcost);
            }
        }
        wcsp->decreaseLb(negcost);
    }

    if (y->unassigned() && z->unassigned()) yz->reconnect();

    // to be done before propagation
    WCSP::elimInfo ei = {this,y,z,(BinaryConstraint*) links[(flag_rev)?1:0].constr, (BinaryConstraint*) links[(flag_rev)?0:1].constr, xyz, NULL};
    wcsp->elimInfos[wcsp->getElimOrder()] = ei;
    wcsp->elimOrderInc();
    yz->propagate();
    return true;
}

void EnumeratedVariable::eliminate()
{
    if(isSep_) return;
    if (ToulBar2::nbDecisionVars>0 && wcspIndex < ToulBar2::nbDecisionVars) return;
    if (ToulBar2::allSolutions && ToulBar2::btdMode!=1 && wcspIndex < ToulBar2::nbvar) return;

    assert(!wcsp->getTreeDec() || wcsp->getTreeDec()->getCluster( cluster )->isActive() );


    if (ToulBar2::elimDegree_preprocessing_ >= 0 &&
            (getDegree() <= min(1,ToulBar2::elimDegree_preprocessing_) ||
                    getTrueDegree() <= ToulBar2::elimDegree_preprocessing_)) {
        if (ToulBar2::elimSpaceMaxMB && wcsp->elimSpace > (Long) ToulBar2::elimSpaceMaxMB * 1024 * 1024) {
            cout << "Generic variable elimination stopped (" <<  wcsp->elimSpace/1024./1024. << " >= " << ToulBar2::elimSpaceMaxMB << " MB)" << endl;
            return;
        }
        wcsp->variableElimination(this);
        return;
    }

    if(ToulBar2::allSolutions && ToulBar2::btdMode==1){
        if (wcsp->getStore()->getDepth()==0) return; // tree decomposition not available yet!
        if (getDegree() > ToulBar2::elimDegree_) return;
        assert(getDegree() == 0);
        assert(wcsp->getTreeDec());
        assert(getDomainSize()>1);
        wcsp->getTreeDec()->getCluster( cluster )->multCount(getDomainSize());
        //		cout << "** "<< wcsp->getTreeDec()->getCluster( cluster )->getCount() << " **" << endl;
    } else {
        if (getDegree() > ToulBar2::elimDegree_) return;
        // if (getDegree()==1 && CSP(wcsp->getLb(),wcsp->getUb()) && constrs->begin().constr->arity() >= 4) {
        // 	// special case: there is only one n-ary constraint on this variable
        // 	// and the current domain size is larger than the number of forbidden tuples
        // 	NaryConstraintMap *ctr = (NaryConstraintMap *) constrs->begin().constr;
        // 	if  (ctr->getDefCost()==0 && ctr->getpf()->size() < getDomainSize()) {
        // 	  // TO BE DONE: find zero-cost unary support for this variable if it exists
        //    // if not, cannot deconnect the constraint!
        // 	  ctr->deconnect();
        // }
        if(getDegree() > 0) {
            TernaryConstraint* ternCtr = existTernary();

            if(ternCtr) { if(!elimVar(ternCtr)) return; }
            else {
                if(getDegree() > 2) return;

                ConstraintLink xylink = *constrs.begin();
                ConstraintLink xzlink = {NULL,0};

                if(xylink.constr->arity() > 2 || !(xylink.constr->extension())) return;

                if(getDegree() == 2) {
                    xzlink = *constrs.rbegin();
                    if(xzlink.constr->arity() > 2 || !(xzlink.constr->extension())) return;

                    if(!elimVar(xylink,xzlink)) return;
                } else {
                    BinaryConstraint* xy = (BinaryConstraint*) xylink.constr;
                    if(!elimVar( xy )) return;
                }
            }
        } else {
            if (ToulBar2::isZ) { // add all unary loglike into lowerbound or negCost
                Cost clogz = wcsp->getUb();
                for (EnumeratedVariable::iterator itv = begin(); itv != end(); ++itv) {
                    clogz = wcsp->LogSumExp(clogz, getCost(*itv));
                }
                if (clogz < MIN_COST) wcsp->decreaseLb(clogz);
                else wcsp->increaseLb(clogz);
            }
        }
    }
    assert(getDegree() == 0);
    if (ToulBar2::verbose >= 2) cout << "Eliminate End of var " << getName() << endl;
    assert(getCost(support) == MIN_COST); // it is ensured by previous calls to findSupport
    assign(support); // warning! dummy assigned value
}

void EnumeratedVariable::permuteDomain(int nperm)
{
    while(nperm) {
        Value a = myrand() % getDomainInitSize();
        Value b = myrand() % getDomainInitSize();
        if(canbe(a) && canbe(b)) permuteDomain(a,b);
        nperm--;
    }
}

void EnumeratedVariable::permuteDomain(Value a, Value b)
{
    assert(canbe(a) && canbe(b));
    for(ConstraintList::iterator iter=constrs.begin(); iter != constrs.end(); ++iter) {
        Constraint* ctr = (*iter).constr;
        if(ctr->arity() == 2 &&  !ctr->isSep() && ctr->extension()) {
            BinaryConstraint* bctr = (BinaryConstraint*) ctr;
            bctr->permute(this,a,b);
        } else {
            cerr <<  "Sorry! cannot sort domains with non-binary cost functions!" << endl;
            exit(EXIT_FAILURE);
        }
    }
}

static int cmpValueCost(const void *p1, const void *p2)
{
    Cost c1 = ((ValueCost *) p1)->cost;
    Cost c2 = ((ValueCost *) p2)->cost;
    Value v1 = ((ValueCost *) p1)->value;
    Value v2 = ((ValueCost *) p2)->value;
    if (c1 < c2) return -1;
    else if (c1 > c2) return 1;
    else if (v1 < v2) return -1;
    else if (v1 > v2) return 1;
    else return 0;
}

/// \warning Must be done just after loading the problem and before any propagate (even projection of unary costs)
/// \warning Does not work with negative domain values!
ValueCost *EnumeratedVariable::sortDomain(vector<Cost> &costs)
{
    cout << "sort variable " << getName() << " (size="<< getDomainSize() << ")" << endl;

    //  // Buble-sort (too slow!!!)
    //  bool swap = false;
    //  do {
    //	swap = false;
    //	iterator iter = begin();
    //	Value a = *iter;
    //	Cost costa = getCost(a);
    //	++iter;
    //	while (iter != end()) {
    //	  Value b = *iter;
    //	  Cost costb = getCost(b);
    //	  if (costa > costb) {
    //		swap = true;
    //		permuteDomain(a,b);
    //	  } else {
    //		a = b;
    //		costa=costb;
    //	  }
    //	  ++iter;
    //	}
    //  } while (swap);

    int size = getDomainInitSize();
    int position[size];
    int value[size];
    for (int i=0; i<size; i++) {
        position[i]=i;
        value[i]=i;
    }
    //  ValueCost sorted[size];
    ValueCost* sorted = new ValueCost[size];
    for (int i=0; i<size; i++) {
        sorted[i].value = i;
        sorted[i].cost = costs[i];
    }
    qsort(sorted, size, sizeof(ValueCost), cmpValueCost);
    for (int i=0; i<size; i++) {
        if (position[toIndex(sorted[i].value)] != i) {
            //		  cout << "swap " << value[i] << "(" << i << ") with " << sorted[i].value << "(" << position[toIndex(sorted[i].value)] << ")" << endl;
            assert(i < position[toIndex(sorted[i].value)]);
            permuteDomain(i,position[toIndex(sorted[i].value)]);
            position[value[i]] = position[toIndex(sorted[i].value)];
            value[position[toIndex(sorted[i].value)]] = value[i];
            value[i] = sorted[i].value;
            position[toIndex(sorted[i].value)] = i;
        }
    }
    for (int i=0; i<size; i++) {
        costs[i] = sorted[i].cost;
    }
    return sorted;
}

bool EnumeratedVariable::canbeMerged(EnumeratedVariable *x)
{
    if (ToulBar2::nbDecisionVars>0 && wcspIndex < ToulBar2::nbDecisionVars) return false;
    if (ToulBar2::allSolutions && wcspIndex < ToulBar2::nbvar && (ToulBar2::elimDegree >= 0 || ToulBar2::elimDegree_preprocessing >= 0)) return false;
    double mult = (1.0 * x->getDomainSize()) / getDomainSize();
    for(ConstraintList::iterator iter=constrs.begin(); iter != constrs.end(); ++iter) {
        Constraint* ctr = (*iter).constr;
        if (!ctr->extension() || ctr->isSep()) return false;
        if (mult>1.1 && ctr->extension() && (ctr->arity() >= 4) && (mult * ((NaryConstraintMap *) ctr)->size() > MAX_NB_TUPLES)) return false;
    }
    return true;
}

// only in preprocessing
void EnumeratedVariable::mergeTo(BinaryConstraint *xy, map<Value, Value> &functional)
{
    assert(wcsp->getStore()->getDepth()==0);
    assert(unassigned());
    EnumeratedVariable *x = (EnumeratedVariable *) xy->getVarDiffFrom(this);
    assert(x->unassigned());
    assert(x != this);
    if (ToulBar2::verbose >= 1) {
        cout << "Merge variable " << *this << " into variable " << *x << endl;
        for (map<Value, Value>::iterator itermap = functional.begin(); itermap != functional.end(); ++itermap) {
            cout << (*itermap).first << "," << (*itermap).second << endl;
        }
    }
    for(ConstraintList::iterator iter=constrs.begin(); iter != constrs.end(); ++iter) {
        Constraint* ctr = (*iter).constr;
        assert(ctr->extension());
        assert(!ctr->isSep());
        if (ctr != xy) ctr->deconnect();
        bool noduplicate = (ctr->getIndex(x) < 0);
        int scopeIndex[ctr->arity()];
        int scopeSize = 0;
        for (int i=0;i<ctr->arity();i++) {
            if (ctr->getVar(i) == this) {
                if (noduplicate) {
                    scopeIndex[scopeSize++] = x->wcspIndex;
                }
            } else {
                scopeIndex[scopeSize++] = ctr->getVar(i)->wcspIndex;
            }
        }
        switch (scopeSize) {
        case 0: {
            cerr << "Error: empty scope from " << *ctr << " when merging functional variable " << *this << " to variable " << *x << endl;
            exit(EXIT_FAILURE);
            break; }
        case 1: {
            assert(scopeIndex[0] == x->wcspIndex);
            assert(ctr->arity() == 2);
            assert(ctr->getIndex(x) >= 0);
            if (ctr != xy) {
                bool modified = false;
                for (iterator iterx = x->begin(); iterx != x->end(); ++iterx) {
                    assert(canbe(functional[*iterx]));
                    Cost cost = ((BinaryConstraint *)ctr)->getCost(x,this,*iterx,functional[*iterx]);
                    if (cost > MIN_COST) {
                        x->project(*iterx,cost);
                        modified = true;
                    }
                }
                if (modified) x->findSupport();
            } // else xy will be eliminated at the end of this function
            break; }
        case 2: {
            assert(wcsp->unassigned(scopeIndex[0]));
            assert(wcsp->unassigned(scopeIndex[1]));
            EnumeratedVariable* u = (EnumeratedVariable *) wcsp->getVar(scopeIndex[0]);
            EnumeratedVariable* v = (EnumeratedVariable *) wcsp->getVar(scopeIndex[1]);
            assert(x==u || x==v);
            vector<Cost> costs(u->getDomainInitSize() * v->getDomainInitSize(), MIN_COST);
            bool empty = true;
            String oldtuple(ctr->arity(),CHAR_FIRST);
            for (EnumeratedVariable::iterator iterU = u->begin(); iterU != u->end(); ++iterU) {
                for (EnumeratedVariable::iterator iterV = v->begin(); iterV != v->end(); ++iterV) {
                    for(int i=0;i<ctr->arity();i++) {
                        if (ctr->getVar(i) == this) {
                            oldtuple[i] = toIndex(functional[(x==u)?(*iterU):(*iterV)]) + CHAR_FIRST;
                        } else {
                            assert(ctr->getVar(i)==u || ctr->getVar(i)==v);
                            oldtuple[i] = ((ctr->getVar(i)==u)?(u->toIndex(*iterU)):(v->toIndex(*iterV))) + CHAR_FIRST;
                        }
                    }
                    Cost cost = ctr->evalsubstr(oldtuple, ctr);
                    if (cost > MIN_COST) {
                        empty = false;
                        costs[u->toIndex(*iterU) * v->getDomainInitSize() + v->toIndex(*iterV)] = cost;
                    }
                }
            }
            if (!empty) wcsp->postBinaryConstraint(scopeIndex[0], scopeIndex[1], costs);
            break; }
        case 3: {
            assert(wcsp->unassigned(scopeIndex[0]));
            assert(wcsp->unassigned(scopeIndex[1]));
            assert(wcsp->unassigned(scopeIndex[2]));
            EnumeratedVariable* u = (EnumeratedVariable *) wcsp->getVar(scopeIndex[0]);
            EnumeratedVariable* v = (EnumeratedVariable *) wcsp->getVar(scopeIndex[1]);
            EnumeratedVariable* w = (EnumeratedVariable *) wcsp->getVar(scopeIndex[2]);
            assert(x==u || x==v || x==w);
            vector<Cost> costs(u->getDomainInitSize() * v->getDomainInitSize() * w->getDomainInitSize(), MIN_COST);
            bool empty = true;
            String oldtuple(ctr->arity(),CHAR_FIRST);
            for (EnumeratedVariable::iterator iterU = u->begin(); iterU != u->end(); ++iterU) {
                for (EnumeratedVariable::iterator iterV = v->begin(); iterV != v->end(); ++iterV) {
                    for (EnumeratedVariable::iterator iterW = w->begin(); iterW != w->end(); ++iterW) {
                        for(int i=0;i<ctr->arity();i++) {
                            if (ctr->getVar(i) == this) {
                                oldtuple[i] = toIndex(functional[(x==u)?(*iterU):((x==v)?(*iterV):(*iterW))]) + CHAR_FIRST;
                            } else {
                                assert(ctr->getVar(i)==u || ctr->getVar(i)==v || ctr->getVar(i)==w);
                                oldtuple[i] = ((ctr->getVar(i)==u)?(u->toIndex(*iterU)):((ctr->getVar(i)==v)?(v->toIndex(*iterV)):(w->toIndex(*iterW)))) + CHAR_FIRST;
                            }
                        }
                        Cost cost = ctr->evalsubstr(oldtuple, ctr);
                        if (cost > MIN_COST) {
                            empty = false;
                            costs[u->toIndex(*iterU) * v->getDomainInitSize() * w->getDomainInitSize() + v->toIndex(*iterV) * w->getDomainInitSize() + w->toIndex(*iterW)] = cost;
                        }
                    }
                }
            }
            if (!empty) wcsp->postTernaryConstraint(scopeIndex[0], scopeIndex[1], scopeIndex[2], costs);
            break; }
        default: {
            //	  int res = wcsp->postNaryConstraintBegin(scopeIndex, scopeSize, MIN_COST);
            //	  NaryConstraintMap *newctrok = (NaryConstraintMap*) wcsp->getCtr(res);
            //	  assert(newctrok->arity() == scopeSize);
            //	  bool empty = true;
            //	  String oldtuple(ctr->arity(),CHAR_FIRST);
            //	  EnumeratedVariable* scopeNewCtr[newctrok->arity()];
            //	  for(int i=0;i<newctrok->arity();i++) {
            //		scopeNewCtr[i] = (EnumeratedVariable*) newctrok->getVar(i);
            //	  }
            //	  String tuple;
            //	  Cost cost;
            //	  newctrok->firstlex();
            //	  while (newctrok->nextlex(tuple,cost)) {
            //		for(int i=0;i<ctr->arity();i++) {
            //		  if (ctr->getVar(i) == this) {
            //			assert(newctrok->getIndex(x)>=0);
            //			oldtuple[i] = toIndex(functional[x->toValue(tuple[newctrok->getIndex(x)] - CHAR_FIRST)]) + CHAR_FIRST;
            //		  } else {
            //			assert(newctrok->getIndex(ctr->getVar(i))>=0);
            //			oldtuple[i] = tuple[newctrok->getIndex(ctr->getVar(i))];
            //		  }
            //		}
            //		cost = ctr->evalsubstr(oldtuple, ctr);
            //		if (cost > MIN_COST) {
            //		  empty = false;
            //		  newctrok->addtoTuple(tuple,cost,scopeNewCtr);
            //		}
            //	  }

            assert(ctr->arity() >= 4);
            NaryConstraintMap *oldctr = (NaryConstraintMap*) ctr;
            Cost defcost = oldctr->getDefCost();
            int res = wcsp->postNaryConstraintBegin(scopeIndex, scopeSize, defcost);
            NaryConstraintMap *newctr = (NaryConstraintMap*) wcsp->getCtr(res);
            assert(newctr->arity() == scopeSize);
            String newtuple(scopeSize,CHAR_FIRST);
            EnumeratedVariable* scopeNewCtr[scopeSize];
            for(int i=0;i<scopeSize;i++) {
                scopeNewCtr[i] = (EnumeratedVariable*) newctr->getVar(i);
            }
            int posx = newctr->getIndex(x);
            assert(posx >= 0);
            int posold = oldctr->getIndex(this);
            assert(posold >= 0);
            int posxold = oldctr->getIndex(x); // check if x is already in the scope of the old cost function
            String tuple;
            Cost cost;
            oldctr->first();
            while (oldctr->next(tuple,cost)) {
                if (cost != defcost) {
                    for(int i=0;i<scopeSize;i++) {
                        if (i != posx || posxold>=0) {
                            assert(oldctr->getIndex(scopeNewCtr[i])>=0);
                            newtuple[i] = tuple[oldctr->getIndex(scopeNewCtr[i])];
                        }
                    }
                    if (posxold>=0) {
                        if (functional[x->toValue(tuple[posxold] - CHAR_FIRST)]==toValue(tuple[posold] - CHAR_FIRST)) {
                            newctr->addtoTuple(newtuple,cost,scopeNewCtr);
                        }
                    } else {
                        for(EnumeratedVariable::iterator iterX = x->begin(); iterX != x->end(); ++iterX) {
                            if (functional[*iterX]==toValue(tuple[posold] - CHAR_FIRST)) {
                                newtuple[posx] = x->toIndex(*iterX) + CHAR_FIRST;
                                newctr->addtoTuple(newtuple,cost,scopeNewCtr);
                            }
                        }
                    }
                }
            }

            //	  newctrok->firstlex();
            //	  while (newctrok->nextlex(tuple,cost)) {
            //		  Cost newcost = newctr->evalsubstr(tuple, newctr);
            //		  if (cost != newcost) {
            //				  int tmpv = ToulBar2::verbose;
            //				  ToulBar2::verbose = 7;
            //				  cout << "BUG: " << cost << " " << newcost << " (" << defcost << ")" << endl;
            //				  for(unsigned int i=0;i<tuple.size();i++) {
            //					  cout << tuple[i] - CHAR_FIRST;
            //				  }
            //				  cout << endl;
            //				  cout << *oldctr << endl;
            //				  cout << *newctrok << endl;
            //				  cout << *newctr << endl;
            //				  ToulBar2::verbose = tmpv;
            //		  }
            //	  }
            //	  newctrok->deconnect(true);

            assert(newctr->connected());
            if (newctr->universal()) newctr->deconnect(true);
            else newctr->propagate();
            break; }
        }
    }
    assert(xy->connected());
    assert(getDegree() == 1);
    elimVar(xy);
    assert(getDegree() == 0);
    assert(getCost(support) == MIN_COST);
    assign(support); // warning! dummy assigned value
}

bool EnumeratedVariable::verify()
{
    TreeDecomposition* td = wcsp->getTreeDec();
    if(!td) return true;
    for(ConstraintList::iterator iter=constrs.begin(); iter != constrs.end(); ++iter) {
        Constraint* ctr1 = (*iter).constr;
        if(ctr1->isSep()) continue;
        if(ctr1->arity() > 3) continue;
        for(ConstraintList::iterator iter2=iter; iter2 != constrs.end(); ++iter2) {
            Constraint* ctr2 = (*iter2).constr;
            if(ctr1 == ctr2) continue;
            if(ctr2->isSep()) continue;
            if(ctr2->arity() > 3) continue;
            if(ctr1->arity() == 3 && ctr2->arity() == 2) {
                TernaryConstraint* tctr1 = (TernaryConstraint*) ctr1;
                BinaryConstraint* bctr2 = (BinaryConstraint*) ctr2;
                if( (tctr1->getIndex( bctr2->getVar(0) ) >= 0 ) &&
                        (tctr1->getIndex( bctr2->getVar(1) ) >= 0 ) ) {
                    if(tctr1->getCluster() != bctr2->getCluster()) {
                        bool t1 = tctr1->xy != bctr2;
                        bool t2 = tctr1->xz != bctr2;
                        bool t3 = tctr1->yz != bctr2;
                        bool t4 = tctr1->isDuplicate() || (!tctr1->xy->isDuplicate() && !tctr1->xz->isDuplicate() && !tctr1->yz->isDuplicate());
                        if (!t1 || !t2 || !t3 || !t4) {
                            cout << "isDuplicate problem: " << tctr1->isDuplicate() << "," << tctr1->xy->isDuplicate() << "," << tctr1->xz->isDuplicate() << "," << tctr1->yz->isDuplicate() << " , " << bctr2->isDuplicate() << endl;
                            cout << *tctr1;
                            cout << *tctr1->xy;
                            cout << *tctr1->xz;
                            cout << *tctr1->yz;
                            cout << *bctr2;
                            assert(false);
                        }
                    } else {
                        bool t1 = tctr1->xy == bctr2 || tctr1->xz == bctr2 || tctr1->yz == bctr2;
                        if (!t1) {
                            cout << "isDuplicate problem: " << tctr1->isDuplicate() << "," << bctr2->isDuplicate() << ")" << endl;
                            cout << *tctr1;
                            cout << *tctr1->xy;
                            cout << *tctr1->xz;
                            cout << *tctr1->yz;
                            cout << *bctr2;
                            assert(false);
                        }
                    }
                }
            }
            if(ctr1->arity() == 2 && ctr2->arity() == 3) {
                BinaryConstraint* bctr1 =  (BinaryConstraint*) ctr1;
                TernaryConstraint* tctr2 = (TernaryConstraint*) ctr2;
                if( (tctr2->getIndex( bctr1->getVar(0) ) >= 0 ) &&
                        (tctr2->getIndex( bctr1->getVar(1) ) >= 0 ) ) {
                    if(bctr1->getCluster() != tctr2->getCluster()) {
                        bool t1 = tctr2->xy != bctr1;
                        bool t2 = tctr2->xz != bctr1;
                        bool t3 = tctr2->yz != bctr1;
                        bool t4 = tctr2->isDuplicate() || (!tctr2->xy->isDuplicate() && !tctr2->xz->isDuplicate() && !tctr2->yz->isDuplicate());
                        if (!t1 || !t2 || !t3 || !t4) {
                            cout << "isDuplicate problem tern/bin: " << tctr2->isDuplicate() << "," << tctr2->xy->isDuplicate() << "," << tctr2->xz->isDuplicate() << "," << tctr2->yz->isDuplicate() << " , " << bctr1->isDuplicate() << endl;
                            cout << *tctr2;
                            cout << *tctr2->xy;
                            cout << *tctr2->xz;
                            cout << *tctr2->yz;
                            cout << *bctr1;
                            assert(false);
                        }
                    } else {
                        bool t1 = tctr2->xy == bctr1 || tctr2->xz == bctr1 || tctr2->yz == bctr1;
                        if (!t1) {
                            cout << "isDuplicate problem: " << tctr2->isDuplicate() << " , " << bctr1->isDuplicate() << endl;
                            cout << *tctr2;
                            cout << *bctr1;
                            assert(false);
                        }
                    }
                }
            }
            if(ctr1->arity() == 3 && ctr2->arity() == 3) {
                TernaryConstraint* tctr1 = (TernaryConstraint*) ctr1;
                TernaryConstraint* tctr2 = (TernaryConstraint*) ctr2;
                BinaryConstraint* bctr1 = tctr1->commonBinary(tctr2);
                BinaryConstraint* bctr2 = tctr2->commonBinary(tctr1);
                if(bctr1) {
                    if(bctr1 != bctr2) {
                        bool t1 = tctr1->getCluster() != tctr2->getCluster();
                        bool t2 = bctr1->getCluster() != bctr2->getCluster();
                        bool t3 = (tctr1->isDuplicate() && bctr1->isDuplicate()) || (tctr2->isDuplicate() && bctr2->isDuplicate());
                        if(!t1 || !t2 || !t3) {
                            cout << "isDuplicate problem tern/tern/bin/bin: " << tctr1->isDuplicate() << "," << tctr2->isDuplicate() << "," << bctr1->isDuplicate() << "," << bctr2->isDuplicate() << endl;
                            cout << *tctr1;
                            cout << *tctr2;
                            cout << *bctr1;
                            cout << *bctr2;
                            assert(false);
                        }

                    } else {
                        assert( tctr1->getCluster() == tctr2->getCluster() );
                        assert( bctr1->getCluster() == bctr2->getCluster() );
                    }
                }
            }
            if(ctr1->arity() == 2 && ctr2->arity() == 2) {
                BinaryConstraint* bctr1 = (BinaryConstraint*) ctr1;
                BinaryConstraint* bctr2 = (BinaryConstraint*) ctr2;
                if( (bctr1->getIndex( bctr2->getVar(0) ) >= 0 ) &&
                        (bctr1->getIndex( bctr2->getVar(1) ) >= 0 ) ) {
                    assert(bctr1->isDuplicate() || bctr2->isDuplicate());
                }
            }
        }
    }
    return true;
}





/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

