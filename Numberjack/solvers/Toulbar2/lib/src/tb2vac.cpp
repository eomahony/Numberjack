/** \file tb2vac.cpp
 *  \brief VAC implementation for binary cost functions.
 *
 *      \defgroup VAC Virtual Arc Consistency enforcing
 *  The three phases of VAC are enforced in three different "Pass".
 *  Bool(P) is never built. Instead specific functions (getVACCost) booleanize the WCSP on the fly.
 *  The domain variables of Bool(P) are the original variable domains (saved and restored using trailing at each iteration) 
 *  All the counter data-structures (k) are timestamped to avoid clearing them at each iteration.
 *  \note Simultaneously AC (and potentially DAC, EAC) are maintained by proper queuing.
 *  \see <em> Soft Arc Consistency Revisited. </em> Cooper et al. Artificial Intelligence. 2010.
 */

#include "tb2vac.hpp"
#include "tb2clusters.hpp"
#include <list>
#include <algorithm>

class tVACStat {
public:

    int var;
    Cost sumlb;
    Long nlb;

    tVACStat(int varin) {
        var = varin;
        sumlb = MIN_COST;
        nlb = 0;
    }};

bool cmp_function(tVACStat * v1, tVACStat * v2)
{
    return v1->sumlb > v2->sumlb;
}

VACExtension::VACExtension(WCSP * w):wcsp(w), VAC2(&w->getStore()->storeVariable), nbIterations(0),
        inconsistentVariable(-1)
{
    queueP = new stack < pair < int, int > >;
    queueR = new stack < pair < int, int > >;
    minlambda = MAX_COST;
    sumlb = MIN_COST;
    sumvars = 0;
    sumk = 0;
    theMaxK = 0;
    nlb = 0;
    // varAssign = -1;
}

VACExtension::~VACExtension()
{
    delete queueP;
    delete queueR;
}

void VACExtension::init()
{
    VACVariable *xi;
    for (unsigned int i = 0; i < wcsp->numberOfVariables(); i++) {
        xi = (VACVariable *) wcsp->getVar(i);
        xi->setThreshold(MIN_COST);
    }
    iniThreshold();
    nearIncVar = NULL;

    // for (unsigned int i = 0; i < wcsp->numberOfVariables(); i++) {
    //    tVACStat* vacinfo = new tVACStat(i);
    //    heapAccess[i] = vacinfo;
    //    heap.insert(heap.end(), vacinfo);
    // }
}

void VACExtension::histogram(Cost c)
{
    if (c != MIN_COST) {
        tScale::iterator it = scaleCost.find(c);
        if (it == scaleCost.end())
            scaleCost[c] = 0;
        else
            it->second++;
    }
}

void VACExtension::histogram()
{
    int cumulus = 0;
    int packetsize = 50;
    bool toomany = true;
    while (toomany) {
        scaleVAC.clear();
        tScale::iterator it = scaleCost.begin();
        while (it != scaleCost.end()) {
            cumulus += it->second;
            if (cumulus > packetsize) {
                scaleVAC.push_front(it->first);
                cumulus = 0;
            }
            ++it;
        }
        toomany = scaleVAC.size() > 20;
        if (toomany)
            packetsize *= 2;
    }

    if (ToulBar2::verbose >= 1) {
        cout << "Reverse Costs Range and Scale: " << scaleCost.rbegin()->first;
        list < Cost >::iterator itl = scaleVAC.begin();
        while (itl != scaleVAC.end()) {
            cout << " " << *itl;
            ++itl;
        }
        cout <<  " " << scaleCost.begin()->first << endl;
    }
}

void VACExtension::iniThreshold()
{
    list < Cost >::iterator it = (scaleVAC.begin());
    Cost c = ((it == scaleVAC.end())?UNIT_COST:(*it));
    if (wcsp->getUb() < c)
        c = wcsp->getUb();
    itThreshold = c;
    nearIncVar = NULL;
}

void VACExtension::nextScaleCost()
{
    Cost c = MAX_COST;
    bool done = false;
    list < Cost >::iterator it = scaleVAC.begin();
    while ((it != scaleVAC.end()) && !done) {
        c = *it;
        done = c < itThreshold;
        ++it;
    }
    if (!done)
        c = itThreshold / (UNIT_COST + UNIT_COST);

    if (wcsp->getStore()->getDepth() == 0) {
        if (c < ToulBar2::costThresholdPre)
            c = MIN_COST;
    } else if (c < ToulBar2::costThreshold)
        c = MIN_COST;

    itThreshold = c;
}

// do not need to revise all variables because we assume soft AC already done
void VACExtension::reset()
{
    VACVariable *x;
    TreeDecomposition *td = wcsp->getTreeDec();
    clear();
    while (!queueP->empty())
        queueP->pop();
    while (!queueR->empty())
        queueR->pop();
    //  int bucket = cost2log2glb(ToulBar2::costThreshold);
    int bucket = cost2log2glb(itThreshold);
    if (bucket < 0)
        bucket = 0;
    for (; bucket < wcsp->getNCBucketSize(); bucket++) {
        VariableList *varlist = wcsp->getNCBucket(bucket);
        for (VariableList::iterator iter = varlist->begin();
                iter != varlist->end();) {
            x = (VACVariable *) * iter;
            if (x->unassigned() && (x->getMaxCost() >= itThreshold)) {
                if (td) {
                    if (td->
                            isActiveAndInCurrentClusterSubTree
                            (x->getCluster()))
                        x->queueVAC();
                } else
                    x->queueVAC();
            }
            ++iter;
        }
    }
    //   for (BTQueue::iterator it = VAC2.begin(); it != VAC2.end(); ++it) {
    //     x = (VACVariable*) (*it);
    //     if (td) { if(td->isActiveAndInCurrentClusterSubTree(x->getCluster())) x->queueVAC(); }
    //  else x->queueVAC();
    //   }
}

bool VACExtension::propagate()
{

    if (wcsp->getStore()->getDepth() >= ToulBar2::vac) {
        return false;
    }
    // if(getVarTimesStat(varAssign) > 100) {
    //    long double m = to_double(getVarCostStat(varAssign))/getVarTimesStat(varAssign);
    //    if(m < to_double(ToulBar2::costMultiplier)/10.) {
    //        inconsistentVariable = -1;
    //        return false;
    //    }
    // }

    Cost ub = wcsp->getUb();
    Cost lb = wcsp->getLb();

    bool isvac = true;
    bool util = true;

    breakCycles = 0;

    vector < pair < VACVariable *, Value > >acSupport;
    bool acSupportOK = false;

    while ((!util || isvac) && itThreshold != MIN_COST) {
        minlambda = ub - lb;
        nbIterations++;
        reset();
        //		if (ToulBar2::verbose>=8) cout << *wcsp;
        wcsp->getStore()->store();
        enforcePass1();
        isvac = isVAC();

        if (ToulBar2::vacValueHeuristic && isvac) {
            acSupportOK = true;
            acSupport.clear();
            // fill SeekSupport with ALL variables if in preprocessing (i.e. before the search)
            if (wcsp->getStore()->getDepth() <= 1 || ToulBar2::debug) {
                for (unsigned int i = 0;
                        i < wcsp->numberOfVariables(); i++) {
                    ((VACVariable *) wcsp->getVar(i))->
                            queueSeekSupport();
                }
            }
            // remember first arc consistent domain values in Bool(P) before restoring domains
            int nbassigned = 0;
            int nbassignedzero = 0;
            while (!SeekSupport.empty()) {
                VACVariable *x =
                        (VACVariable *) SeekSupport.pop();
                if (x->assigned()) nbassigned++;
                pair < VACVariable *, Value > p;
                p.first = x;
                p.second = x->getSup() + 1;
                for (EnumeratedVariable::iterator iterX =
                        x->begin(); iterX != x->end(); ++iterX) {
                    if (x->getCost(*iterX) == MIN_COST) {
                        if (x->assigned()) nbassignedzero++;
                        p.second = *iterX;
                        break;
                    }
                }
                if (x->canbe(p.second))
                    acSupport.push_back(p);
            }
            if (ToulBar2::debug && nbassignedzero>0) cout << "[" << wcsp->getStore()->getDepth() << "] " << nbassignedzero << "/" << nbassigned-nbassignedzero << "/" << wcsp->numberOfUnassignedVariables() << " fixed/singletonnonzerocost/unassigned" << endl;
        }
        wcsp->getStore()->restore();

        if (!isvac) {
            enforcePass2();
            if (ToulBar2::verbose > 0)
                cout << "VAC Lb: " << wcsp->
                getLb() << "    incvar: " <<
                inconsistentVariable << "    minlambda: " <<
                minlambda << "      itThreshold: " <<
                itThreshold << endl;
            util = enforcePass3();
        } else {
            nextScaleCost();
            //if(nearIncVar) cout << "var: " << nearIncVar->wcspIndex << "  at Cost: " << atThreshold << endl;
        }
    }

    if (ToulBar2::vacValueHeuristic && acSupportOK) {
        // update current unary support if possible && needed
        for (vector < pair < VACVariable *, Value > >::iterator iter =
                acSupport.begin(); iter != acSupport.end(); ++iter) {
            VACVariable *x = iter->first;
            if (x->canbe(iter->second)) {
                assert(x->getCost(iter->second) == MIN_COST);
                if (ToulBar2::verbose > 0
                        && x->getSupport() != iter->second)
                    cout << "CHANGE SUPPORT " << x->
                    getName() << " from " << x->
                    getSupport() << " to " << iter->
                    second << endl;
                x->setSupport(iter->second);
            } else {
                if (ToulBar2::verbose > 0)
                    cout << "WARNING: BAD AC SUPPORT " <<
                    iter->
                    second << " FOR VARIABLE " << *x <<
                    endl;
            }
        }
    }
    //  updateStat(wcsp->getLb() - lb);
    //if(isvac) assert(checkPass1());
    return util;
}

// DO NOTHING!!! See in tb2solver.cpp variable heuristics
int VACExtension::getHeuristic()
{
    //if(nearIncVar) return nearIncVar->wcspIndex; else
    return -1;
}

bool VACExtension::enforcePass1(VACVariable * xj, VACBinaryConstraint * cij)
{
    bool wipeout = false;
    VACVariable *xi;
    xi = (VACVariable *) cij->getVarDiffFrom(xj);
    for (EnumeratedVariable::iterator it = xi->begin(); it != xi->end();
            ++it) {
        Value v = *it;
        if (xi->getVACCost(v) != MIN_COST) {
            xi->removeVAC(v);
        }		// xi->queueVAC(); }
        else if (cij->revise(xi, v)) {
            wipeout = xi->removeVAC(v);
            xi->setKiller(v, xj->wcspIndex);
            xj->killedOne();	// HEUR
            queueP->push(pair < int, int >(xi->wcspIndex, v));
            xi->queueVAC();
            if (ToulBar2::vacValueHeuristic)
                xi->queueSeekSupport();
            if (wipeout) {
                inconsistentVariable = xi->wcspIndex;
                return true;
            }
        }
    }
    //  if((xi->getDomainSize() == 1) && (!nearIncVar)) {
    //    nearIncVar = xi;
    //    atThreshold = itThreshold;
    //  }
    //  if((xj->getDomainSize() == 1) && (!nearIncVar)) {
    //    nearIncVar = xj;
    //    atThreshold = itThreshold;
    //  }

    return false;
}

void VACExtension::enforcePass1()
{
    //  VACVariable* xi;
    VACVariable *xj;
    VACBinaryConstraint *cij;
    //if (ToulBar2::verbose > 1) cout << "VAC Enforce Pass 1" << endl;

    while (!VAC.empty()) {
        xj = (VACVariable *) VAC.pop_first();
        //list<Constraint*> l;
        for (ConstraintList::iterator itc = xj->getConstrs()->begin();
                itc != xj->getConstrs()->end(); ++itc) {
            cij = (VACBinaryConstraint *) (*itc).constr;
            if (cij->arity() == 2 && !cij->isSep()) {
                //        xi = (VACVariable *)cij->getVarDiffFrom(xj);
                //if(xj->getMaxK(nbIterations) > 2) l.push_back(cij); else
                if (enforcePass1(xj, cij))
                    return;
            }
        }

        /*for (list<Constraint*>::iterator itl = l.begin(); itl != l.end(); ++itl) {
		   cij = (VACConstraint *) *itl;
		   if(enforcePass1(xj,cij)) return;
		   } */
    }
    inconsistentVariable = -1;
}

bool VACExtension::checkPass1() const
{
    VACBinaryConstraint *cij;
    VACVariable *xi, *xj;
    bool supportFound;
    TreeDecomposition *td = wcsp->getTreeDec();

    for (unsigned int i = 0; i < wcsp->numberOfVariables(); i++) {
        xi = (VACVariable *) wcsp->getVar(i);
        if (td
                && !td->isActiveAndInCurrentClusterSubTree(xi->
                        getCluster()))
            continue;
        for (ConstraintList::iterator iter = xi->getConstrs()->begin();
                iter != xj->getConstrs()->end(); ++iter) {
            cij = (VACBinaryConstraint *) (*iter).constr;
            if (cij->arity() == 2 && !cij->isSep()) {
                xj = (VACVariable *) cij->getVarDiffFrom(xi);
                for (EnumeratedVariable::iterator iti =
                        xi->begin(); iti != xi->end(); ++iti) {
                    Value v = *iti;
                    supportFound = false;
                    for (EnumeratedVariable::iterator itj =
                            xj->begin(); itj != xj->end();
                            ++itj) {
                        Value w = *itj;
                        if ((xi->getVACCost(v) ==
                                MIN_COST)
                                && (xj->getVACCost(w) ==
                                        MIN_COST)
                                        && (cij->
                                                getVACCost(xi, xj, v,
                                                        w) ==
                                                                MIN_COST)) {
                            supportFound = true;
                            break;
                        }
                    }
                    if (!supportFound) {
                        return false;
                    }
                }
            }
        }
    }
    return true;
}

bool VACExtension::isVAC() const
{
    return (inconsistentVariable == -1);
}

void VACExtension::enforcePass2()
{
    int i0 = inconsistentVariable;
    int i, j;
    VACVariable *xi0, *xi, *xj;
    VACBinaryConstraint *cij;
    Cost tmplambda;
    Value v;

    //if (ToulBar2::verbose > 0)  cout << "VAC Enforce Pass 2" << endl;

    assert(i0 >= 0);
    xi0 = (VACVariable *) wcsp->getVar(i0);

    for (EnumeratedVariable::iterator iti0 = xi0->begin();
            iti0 != xi0->end(); ++iti0) {
        v = *iti0;
        xi0->addToK(v, 1, nbIterations);
        Cost cost = xi0->getVACCost(v);
        if (cost > MIN_COST) {
            if (cost < minlambda) {
                minlambda = cost;	// NB: we don't need to check for bottleneck here as k=1 neesssarily
            }
        } else
            xi0->setMark(v, nbIterations);
    }

    while (!queueP->empty()) {
        i = queueP->top().first;
        v = queueP->top().second;
        queueP->pop();
        xi = (VACVariable *) wcsp->getVar(i);
        if (xi->isMarked(v, nbIterations)) {
            j = xi->getKiller(v);
            xj = (VACVariable *) wcsp->getVar(j);
            queueR->push(pair < int, int >(i, v));
            cij = (VACBinaryConstraint *) xi->getConstr(xj);
            assert(cij);
            //if (ToulBar2::verbose > 6) cout << "x" << xi->wcspIndex << "," << v << "   killer: " << xj->wcspIndex << endl;

            for (EnumeratedVariable::iterator itj = xj->begin();
                    itj != xj->end(); ++itj) {
                Value w = *itj;
                Cost costij = cij->getVACCost(xi, xj, v, w);
                if (costij > MIN_COST) {
                    int tmpK = xi->getK(v, nbIterations);
                    if (xj->getKiller(w) == i
                            && xj->isMarked(w, nbIterations))
                        tmpK +=
                                xj->getK(w, nbIterations);
                    if (!CUT
                            (wcsp->getLb() + costij,
                                    wcsp->getUb())) {
                        if ((costij / tmpK) < minlambda) {
                            minlambda =
                                    costij / tmpK;
                            if (minlambda < UNIT_COST) {	//A cost function bottleneck here !
                                bneckCost =
                                        costij;
                                bneckCF = cij;
                                bneckVar = -1;
                            }
                        }
                    }
                } else {
                    int tmpK =
                            xi->getK(v,
                                    nbIterations) -
                                    cij->getK(xj, w, nbIterations);
                    if (tmpK > 0) {
                        xj->addToK(w, tmpK,
                                nbIterations);
                        cij->setK(xj, w,
                                xi->getK(v,
                                        nbIterations),
                                        nbIterations);
                        Cost cost = xj->getVACCost(w);
                        if (cost == MIN_COST)
                            xj->setMark(w,
                                    nbIterations);
                        else if (!CUT
                                (wcsp->getLb() + cost,
                                        wcsp->getUb())) {
                            tmplambda =
                                    cost / xj->getK(w,
                                            nbIterations);
                            if (tmplambda <
                                    minlambda) {
                                minlambda =
                                        tmplambda;
                                if (minlambda < UNIT_COST) {	// A unary cost bottleneck here
                                    bneckVar
                                    = j;
                                    bneckCF
                                    =
                                            NULL;
                                    bneckCost
                                    =
                                            cost;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    //if (maxK == 0) {
    //  maxK = wcsp->getUb() - wcsp->getLb();
    //}
    if (ToulBar2::verbose > 1)
        cout << "minLambda: " << minlambda << "\t\t (lb = " << wcsp->
        getLb() << ", ub = " << wcsp->getUb() << ")" << endl;
}

bool VACExtension::enforcePass3()
{
    bool util = (minlambda >= UNIT_COST);
    /*if(util) {
	   Cost ub = wcsp->getUb();
	   Cost lb = wcsp->getLb();
	   util = ( (ub - lb)/(10000*ToulBar2::costMultiplier) ) < minlambda;
	   } */
    Cost lambda = minlambda;

    //if (ToulBar2::verbose > 2) cout << "VAC Enforce Pass 3.   minlambda " << minlambda << " , var: " << inconsistentVariable << endl;
    int i, j;
    VACVariable *xi, *xj;
    VACBinaryConstraint *cij;
    int i0 = inconsistentVariable;
    VACVariable *xi0 = (VACVariable *) wcsp->getVar(i0);
    Value w;

    int maxk = 0;

    if (!util) {		// Empty R ?
        assert(bneckVar != -1 || bneckCF != NULL);
        if (bneckVar != -1) {
            xi = (VACVariable *) wcsp->getVar(bneckVar);
            xi->setThreshold(bneckCost + UNIT_COST);
        } else {
            bneckCF->setThreshold(bneckCost + UNIT_COST);
        }
        breakCycles++;
        //if (ToulBar2::verbose > 1) cout << "BreakCycle   (var: " << minvar << ", val= " << minval << ")   thr: " <<  xi->getThreshold() << endl;
        if (breakCycles > 5) {
            inconsistentVariable = -1;
            itThreshold = MIN_COST;
        }
        //inconsistentVariable = -1; itThreshold = MIN_COST;
        return false;
    }
    // update general stats
    nlb++;
    sumlb += lambda;
    sumvars += queueR->size();
    maxk = 0;

    while (!queueR->empty()) {
        j = queueR->top().first;
        w = queueR->top().second;
        queueR->pop();
        xj = (VACVariable *) wcsp->getVar(j);
        i = xj->getKiller(w);
        xi = (VACVariable *) wcsp->getVar(i);
        cij = (VACBinaryConstraint *) xi->getConstr(xj);
        assert(cij);

        int xjk = xj->getK(w, nbIterations);
        if (maxk < xjk)
            maxk = xjk;

        for (EnumeratedVariable::iterator iti = xi->begin();
                iti != xi->end(); ++iti) {
            Value v = *iti;
            if (cij->getK(xi, v, nbIterations) != 0) {
                Cost ecost =
                        lambda * cij->getK(xi, v, nbIterations);
                cij->setK(xi, v, 0, nbIterations);
                cij->VACextend(xi, v, ecost);
                // extention from unary to binary cost function may break soft AC/DAC in both directions due to isNull/itThreshold
                if (ToulBar2::LcLevel == LC_AC) {
                    xi->queueAC();
                    xj->queueAC();
                } else {
                    if (cij->getDACScopeIndex() ==
                            cij->getIndex(xi)) {
                        xi->queueAC();
                        xi->queueEAC1();
                        xj->queueDAC();
                    } else {
                        xi->queueDAC();
                        xj->queueAC();
                        xj->queueEAC1();
                    }
                }
            }
        }
        cij->VACproject(xj, w, lambda * xj->getK(w, nbIterations));
    }
    sumk += maxk;
    if (maxk > theMaxK)
        theMaxK = maxk;

    xi0->extendAll(lambda);
    xi0->projectLB(lambda);
    return true;
}

// void VACExtension::updateStat(Cost lambda)
// {
//   //tVACStat* v = heapAccess[inconsistentVariable];  
//   if(varAssign >= 0) {
//    tVACStat* v = heapAccess[varAssign];  
//    v->sumlb += lambda;
//    v->nlb++;
//   }
// }

// Cost VACExtension::getVarCostStat( int i )
// {
//   tVACStat* v = heap[i];
//   return v->sumlb;   
// }

// Long VACExtension::getVarTimesStat( int i )
// {
//   if(i < 0) return 0; 
//   tVACStat* v = heap[i];
//   return v->nlb; 
// }

void VACExtension::afterPreprocessing()
{
    int discarded = 0;
    for (unsigned int i = 0; i < wcsp->numberOfConstraints(); i++) {
        Constraint *c = wcsp->getCtr(i);
        if (c->connected() && (c->arity() <= 3) && !c->isSep()) {
            if (c->getTightness() <
                    to_double(ToulBar2::relaxThreshold)) {
                c->deconnect();
                discarded++;
            }
        }
    }
    if (discarded)
        cout << "WARNING num of discarded ctrs: " << discarded << endl;
}

void VACExtension::iniSingleton()
{
    singletonI.clear();
    for (unsigned int i = 0; i < wcsp->numberOfVariables(); i++) {
        int size = wcsp->getDomainSize(i);
        for (int a = 0; a < size; a++)
            singletonI.insert(MAX_DOMAIN_SIZE * i + a);
    }
}

void VACExtension::updateSingleton()
{
    set < int >&s1 = singleton;
    set < int >s2(singletonI);
    singletonI.clear();
    set_intersection(s1.begin(), s1.end(),
            s2.begin(), s2.end(),
            inserter(singletonI, singletonI.begin()));
    singleton.clear();
}

void VACExtension::removeSingleton()
{
    set < int >&s = singletonI;
    set < int >::iterator it = s.begin();
    while (it != s.end()) {
        int ivar = *it / MAX_DOMAIN_SIZE;
        Value a = *it % MAX_DOMAIN_SIZE;
        Variable *var = wcsp->getVar(ivar);
        var->remove(a);
        var->queueNC();
        ++it;
    }
    wcsp->propagate();
}

void VACExtension::clear()
{
    while (!VAC.empty())
        VAC.pop();
    if (ToulBar2::vacValueHeuristic)
        while (!SeekSupport.empty())
            SeekSupport.pop();
}

void VACExtension::queueVAC(DLink < VariableWithTimeStamp > *link)
{
    assert(ToulBar2::vac);
    VAC.push(link, wcsp->getNbNodes());
}

void VACExtension::queueSeekSupport(DLink < VariableWithTimeStamp > *link)
{
    assert(ToulBar2::vac);
    SeekSupport.push(link, wcsp->getNbNodes());
}

void VACExtension::queueVAC2(DLink < Variable * >*link)
{
    assert(ToulBar2::vac);
    VAC2.push(link);
}

void VACExtension::dequeueVAC2(DLink < Variable * >*link)
{
    assert(ToulBar2::vac);
    VAC2.remove(link);
}

void VACExtension::printStat(bool ini)
{
    long double mean = to_double(sumlb) / (long double)nlb;
    cout << "VAC mean lb/incr: " << mean << "     total increments: " << nlb
            << "     cyclesize: " << (double)sumvars /
            (double)nlb << "     k: " << (double)sumk /
            (double)nlb << " (mean), " << theMaxK << " (max)" << endl;
    if (ini)
        cout << "Lb after VAC: " << wcsp->getLb() << endl;

    //sort(heap.begin(), heap.end(), cmp_function);
    /*cout << "Vars: ";
	   vector<tVACStat*>::iterator it = heap.begin();
	   while(it != heap.end()) {
	   tVACStat* v = *it;
	   if(v->sumlb != MIN_COST) cout << "(" << v->var << "," << v->sumlb << ") "; 
	   ++it;
	   }
	   cout << endl; */

    sumk = 0;
    theMaxK = 0;
    sumvars = 0;
    sumlb = MIN_COST;
    nlb = 0;
}

void VACExtension::printTightMatrix()
{
    ofstream ofs("problem.dat");

    Cost Top = wcsp->getUb();
    for (unsigned int i = 0; i < wcsp->numberOfVariables(); i++) {
        for (unsigned int j = 0; j < wcsp->numberOfVariables(); j++) {
            if (i != j) {
                EnumeratedVariable *x =
                        (EnumeratedVariable *) wcsp->getVar(i);
                EnumeratedVariable *y =
                        (EnumeratedVariable *) wcsp->getVar(j);
                Constraint *bctr = x->getConstr(y);
                double t = 0;
                if (bctr)
                    t = bctr->getTightness();
                if (t > to_double(Top))
                    t = to_double(Top);
                t = t * 256.0 / to_double(Top);
                ofs << t << " ";
            } else
                ofs << 0 << " ";
        }
        ofs << endl;
    }
}

/* Min-Sum diffusion algorithm */
void VACExtension::minsumDiffusion()
{
    for (int times = 0; times < 2; times++) {
        bool change = true;
        int maxit = ToulBar2::minsumDiffusion;
        cout << "MinSumDiffusion: " << endl;
        cout << "   max iterations " << maxit << endl;
        cout << "   C0 = " << wcsp->getLb() << endl;
        int ntimes = 0;
        while (change && (ntimes < maxit)) {
            change = false;
            int nchanged = 0;
            for (unsigned int i = 0; i < wcsp->numberOfVariables();
                    i++)
                if (wcsp->unassigned(i)) {
                    VACVariable *evar =
                            (VACVariable *) wcsp->getVar(i);
                    if (evar->averaging()) {
                        change = true;
                        nchanged++;
                        evar->findSupport();
                    }
                }
            ntimes++;
            //cout << "it " << ntimes << "   changed: " << nchanged << endl;
        }
        cout << "   done iterations: " << ntimes << endl;
        for (unsigned int i = 0; i < wcsp->numberOfVariables(); i++)
            if (wcsp->unassigned(i)) {
                EnumeratedVariable *evar =
                        (EnumeratedVariable *) wcsp->getVar(i);
                evar->findSupport();
            }
        for (unsigned int i = 0; i < wcsp->numberOfConstraints(); i++)
            if (wcsp->getCtr(i)->connected())
                wcsp->getCtr(i)->propagate();
        for (int i = 0; i < wcsp->getElimBinOrder(); i++)
            if (wcsp->getElimBinCtr(i)->connected()
                    && !wcsp->getElimBinCtr(i)->isSep())
                wcsp->getElimBinCtr(i)->propagate();
        for (int i = 0; i < wcsp->getElimTernOrder(); i++)
            if (wcsp->getElimTernCtr(i)->connected()
                    && !wcsp->getElimTernCtr(i)->isSep())
                wcsp->getElimTernCtr(i)->propagate();
        wcsp->propagate();
        cout << "   C0 = " << wcsp->getLb() << endl;
        //    printTightMatrix();
    }
}

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

