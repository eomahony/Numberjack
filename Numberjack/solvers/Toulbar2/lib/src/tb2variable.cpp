/*
 * ****** Variable with domain represented by an interval or an enumerated domain *******
 */

#include "tb2variable.hpp"
#include "tb2wcsp.hpp"
#include "tb2binconstr.hpp"
#include "tb2ternaryconstr.hpp"
#include "tb2clusters.hpp"

/*
 * Constructors and misc.
 * 
 */

Variable::Variable(WCSP *w, string n, Value iinf, Value isup) : WCSPLink(w,w->numberOfVariables()), name(n), dac(w->numberOfVariables()),
        timestamp(-1), pos(-1),
        inf(iinf, &w->getStore()->storeValue), sup(isup, &w->getStore()->storeValue), 
        constrs(&w->getStore()->storeConstraint),
        //triangles(&w->getStore()->storeConstraint),
        maxCost(MIN_COST, &w->getStore()->storeCost), maxCostValue(iinf, &w->getStore()->storeValue), 
        NCBucket(-1, &w->getStore()->storeInt),
        cluster(-1, &w->getStore()->storeInt)
{
    if (w->getStore()->getDepth() > 0) {
        cerr << "You cannot create a variable during the search!" << endl;
        exit(EXIT_FAILURE);
    }
    w->link(this);

    linkNCBucket.content = this;
    linkNCQueue.content.var = this;
    linkNCQueue.content.timeStamp = -1;
    linkIncDecQueue.content.var = this;
    linkIncDecQueue.content.timeStamp = -1;
    linkIncDecQueue.content.incdec = NOTHING_EVENT;
    linkEliminateQueue.content.var = this;
    linkEliminateQueue.content.timeStamp = -1;
    isSep_ = false;
}

DLink<ConstraintLink> *Variable::link(Constraint *c, int index)
{
    ConstraintLink e;
    e.constr = c;
    e.scopeIndex = index;
    DLink<ConstraintLink> *elt = new DLink<ConstraintLink>;
    elt->content = e;
    //    if (c->isTriangle()) triangles.push_back(elt,true);
    //    else
    constrs.push_back(elt,true);
    return elt;
}

int Variable::getCurrentVarId()
{
    if (assigned()) return -1;
    if (wcsp->getNbNodes() > timestamp) wcsp->updateCurrentVarsId();
    assert(pos>=0);
    assert(wcsp->getNbNodes() == timestamp);
    return pos;
}

void Variable::setCurrentVarId(int idx)
{
    pos=idx;
    timestamp=wcsp->getNbNodes();
}

int cmpConstraint(const void *p1, const void *p2)
{
    DLink<ConstraintLink> *c1 = *((DLink<ConstraintLink> **) p1);
    DLink<ConstraintLink> *c2 = *((DLink<ConstraintLink> **) p2);
    int v1 = c1->content.constr->getSmallestVarIndexInScope(c1->content.scopeIndex);
    int v2 = c2->content.constr->getSmallestVarIndexInScope(c2->content.scopeIndex);
    if (v1 < v2) return -1;
    else if (v1 > v2) return 1;
    else return 0;
}

int cmpConstraintDAC(const void *p1, const void *p2)
{
    DLink<ConstraintLink> *c1 = *((DLink<ConstraintLink> **) p1);
    DLink<ConstraintLink> *c2 = *((DLink<ConstraintLink> **) p2);
    int v1 = c1->content.constr->getSmallestDACIndexInScope(c1->content.scopeIndex);
    int v2 = c2->content.constr->getSmallestDACIndexInScope(c2->content.scopeIndex);
    if (v1 < v2) return 1;
    else if (v1 > v2) return -1;
    else return 0;
}

int cmpConstraintTightness(const void *p1, const void *p2)
{
    DLink<ConstraintLink> *c1 = *((DLink<ConstraintLink> **) p1);
    DLink<ConstraintLink> *c2 = *((DLink<ConstraintLink> **) p2);
    double v1 = c1->content.constr->getTightness();
    double v2 = c2->content.constr->getTightness();
    if (v1 < v2) return 1;
    else if (v1 > v2) return -1;
    else return 0;
}

void Variable::sortConstraints()
{
    int size = constrs.getSize();
    DLink<ConstraintLink> **sorted = new DLink<ConstraintLink> * [size]; // replace size by MAX_DOMAIN_SIZE in case of compilation problem
    int i=0;
    for (ConstraintList::iterator iter = constrs.begin(); iter != constrs.end(); ++iter) {
        sorted[i++] = iter.getElt();
    }
    qsort(sorted, size, sizeof(DLink<ConstraintLink> *), cmpConstraintDAC);
    for (int i = 0; i < size; i++) {
        constrs.erase(sorted[i],true);
        constrs.push_back(sorted[i],true);
    }
    delete [] sorted;
}

void Variable::deconnect(DLink<ConstraintLink> *link, bool reuse)
{
    if (!link->removed) {
        //        if (link->content.constr->isTriangle()) getTriangles()->erase(link, true);
        //        else
        getConstrs()->erase(link, true);

        if (getDegree() <= ToulBar2::elimDegree_ ||
                (ToulBar2::elimDegree_preprocessing_ >= 0 &&
                        (getDegree() <= min(1,ToulBar2::elimDegree_preprocessing_) ||
                                getTrueDegree() <= ToulBar2::elimDegree_preprocessing_))) queueEliminate();
    }
    if (reuse) {
        assert(wcsp->getStore()->getDepth()==0);
        link->prev = NULL;
        link->next = NULL;
    }
}

int Variable::getTrueDegree()
{
    //	if (constrs.getSize() >= ToulBar2::weightedDegree) return getDegree(); ///\warning returns an approximate degree if the constraint list is too large!
    //    TSCOPE scope1,scope2,scope3;
    set<int> scope1;
    for (ConstraintList::iterator iter=constrs.begin(); iter != constrs.end(); ++iter) {
        if((*iter).constr->isSep()) continue;
        for(int k=0; k < (*iter).constr->arity(); k++) {
            scope1.insert((*iter).constr->getVar(k)->wcspIndex);
        }
        //		(*iter).constr->getScope(scope2);
        //		for (TSCOPE::iterator iter2=scope2.begin(); iter2 != scope2.end(); ++iter2) {
        //		   scope1.insert( iter2->first );
        //		}
        //     	set_union( scope1.begin(), scope1.end(),
        //	  		   	   scope2.begin(), scope2.end(),
        //			  	   inserter(scope3, scope3.begin()) );
        //		scope1 = scope3;
        //		scope3.clear();
    }
    if (scope1.size() >= 1) return scope1.size() - 1;
    else return 0;
}

Long Variable::getWeightedDegree()
{
    Long res = 0;
    for (ConstraintList::iterator iter=constrs.begin(); iter != constrs.end(); ++iter) {
        //    	if((*iter).constr->isSep()) continue;
        res += (*iter).constr->getConflictWeight((*iter).scopeIndex);
        if((*iter).constr->isSep()) res--;  // do not count unused separators
    }
    return res;
}

void Variable::resetWeightedDegree()
{
    for (ConstraintList::iterator iter=constrs.begin(); iter != constrs.end(); ++iter) {
        if((*iter).constr->isSep()) continue;
        (*iter).constr->resetConflictWeight();
    }
}

void Variable::conflict()
{
    wcsp->conflict();
}

/*
 * Propagation methods
 * 
 */

void Variable::queueNC()
{
    wcsp->queueNC(&linkNCQueue);
}

void Variable::queueInc()
{
    wcsp->queueInc(&linkIncDecQueue);
}

void Variable::queueDec()
{
    wcsp->queueDec(&linkIncDecQueue);
}

void Variable::queueEliminate()
{
    wcsp->queueEliminate(&linkEliminateQueue);
}

void Variable::changeNCBucket(int newBucket)
{
    if (NCBucket != newBucket) {
        if (ToulBar2::verbose >= 3) cout << "changeNCbucket " << getName() << ": " << NCBucket << " -> " << newBucket << endl;
        wcsp->changeNCBucket(NCBucket, newBucket, &linkNCBucket);
        NCBucket = newBucket;
    }
}

void Variable::setMaxUnaryCost(Value a, Cost cost)
{
    assert(canbe(a));
    maxCostValue = a;
    assert(cost >= MIN_COST);
    if (maxCost != cost) {
        if (cost > maxCost) queueDEE();
        maxCost = cost;
        int newbucket = min(cost2log2gub(cost), wcsp->getNCBucketSize() - 1);
        changeNCBucket(newbucket);
    }
}

void Variable::projectLB(Cost cost)
{
    if (cost == 0) return;
    if (ToulBar2::verbose >= 2) cout << "lower bound increased " << wcsp->getLb() << " -> " << wcsp->getLb()+cost << endl;
    wcsp->increaseLb(cost); // done before cluster LB because of #CSP (assuming a contradiction will occur here)
    if (wcsp->td) {
        if (ToulBar2::verbose >= 2) cout << " in cluster C" << getCluster() << " (from " << wcsp->td->getCluster(getCluster())->getLb() << " to " << wcsp->td->getCluster(getCluster())->getLb() + cost << ")" << endl;
        wcsp->td->getCluster(getCluster())->increaseLb(cost);
    }
}

void Variable::propagateIncDec(int incdec)
{
    for (ConstraintList::iterator iter=constrs.begin(); iter != constrs.end(); ++iter) {
        if (incdec & INCREASE_EVENT) {
            (*iter).constr->increase((*iter).scopeIndex);
        }
        if ((*iter).constr->connected() && (incdec & DECREASE_EVENT)) {
            (*iter).constr->decrease((*iter).scopeIndex);
        }
    }
}

// Looks for the constraint that links this variable with x
BinaryConstraint* Variable::getConstr( Variable* x )
{
    BinaryConstraint* ctr2;
    TernaryConstraint* ctr3;
    for (ConstraintList::iterator iter=constrs.begin(); iter != constrs.end(); ++iter) {
        if((*iter).constr->isSep()) continue;

        if ((*iter).constr->arity() == 2) {
            ctr2 = (BinaryConstraint*) (*iter).constr;
            if(ctr2->getIndex(x) >= 0) return ctr2;
        }
        else if ((*iter).constr->arity() == 3) {
            ctr3 = (TernaryConstraint*) (*iter).constr;
            int idx = ctr3->getIndex(x);
            if(idx >= 0) {
                int idt = (*iter).scopeIndex;
                if((0 != idx) && (0 != idt)) return ctr3->yz;
                else if((1 != idx) && (1 != idt)) return ctr3->xz;
                else return ctr3->xy;
            }
        }
    }
    return NULL;
}     


BinaryConstraint* Variable::getConstr( Variable* x, int cid )
{
    BinaryConstraint* res;
    BinaryConstraint* ctr2;
    TernaryConstraint* ctr3;
    for (ConstraintList::iterator iter=constrs.begin(); iter != constrs.end(); ++iter) {
        if((*iter).constr->isSep()) continue;

        if ((*iter).constr->arity() == 2) {
            ctr2 = (BinaryConstraint*) (*iter).constr;
            if(ctr2->getIndex(x) >= 0) {
                res = ctr2;
                if(res->getCluster() == cid) return res;
            }
        }
        else if ((*iter).constr->arity() == 3) {
            ctr3 = (TernaryConstraint*) (*iter).constr;
            int idx = ctr3->getIndex(x);
            if(idx >= 0) {
                int idt = (*iter).scopeIndex;
                if((0 != idx) && (0 != idt)) res = ctr3->yz;
                else if((1 != idx) && (1 != idt)) res = ctr3->xz;
                else res = ctr3->xy;

                if(res && res->getCluster() == cid) return res;
            }
        }
    }
    return NULL;
}     



// Looks for the ternary constraint that links this variable with x and y
TernaryConstraint* Variable::getConstr( Variable* x, Variable* y )
{
    TernaryConstraint* ctr;
    for (ConstraintList::iterator iter=constrs.begin(); iter != constrs.end(); ++iter) {
        if((*iter).constr->isSep()) continue;

        if ((*iter).constr->arity() == 3) {
            ctr = (TernaryConstraint*) (*iter).constr;
            if((ctr->getIndex(x)  >= 0) && (ctr->getIndex(y)  >= 0)) return ctr;
        }
    }
    return NULL;
}

TernaryConstraint* Variable::getConstr( Variable* x, Variable* y, int cid )
{
    TernaryConstraint* ctr = NULL;
    for (ConstraintList::iterator iter=constrs.begin(); iter != constrs.end(); ++iter) {
        if((*iter).constr->isSep()) continue;

        if ((*iter).constr->arity() == 3) {
            ctr = (TernaryConstraint*) (*iter).constr;
            if((ctr->getIndex(x)  >= 0) && (ctr->getIndex(y)  >= 0)) {
                if(ctr->getCluster() == cid) return ctr;
            }
        }
    }
    return NULL;
}



// returns a ternary constraint if the current variable is linked to one
TernaryConstraint* Variable::existTernary()
{
    TernaryConstraint* ctr;
    for (ConstraintList::iterator iter=constrs.begin(); iter != constrs.end(); ++iter) {
        if((*iter).constr->isSep()) continue;
        if ((*iter).constr->extension() && (*iter).constr->arity() == 3) {
            ctr = (TernaryConstraint*) (*iter).constr;
            return ctr;
        }
    }
    return NULL;
}







double Variable::strongLinkedby( Variable* &strvar,  TernaryConstraint* &tctr1max, TernaryConstraint* &tctr2max  ) {
    double maxtight = -1;
    strvar = NULL; tctr1max = NULL; tctr2max = NULL;

    TernaryConstraint *tctr1 = NULL;

    for(ConstraintList::iterator iter=constrs.begin(); iter != constrs.end(); ++iter) {
        if((*iter).constr->isSep()) continue;
        if((*iter).constr->arity() == 2) {
            BinaryConstraint* bctr = (BinaryConstraint*) (*iter).constr;
            double bintight = bctr->getTightness();
            if(bintight > maxtight) { maxtight = bintight; strvar = wcsp->getVar(bctr->getSmallestVarIndexInScope((*iter).scopeIndex)); tctr1max = NULL; tctr2max = NULL; }
        }
        else if((*iter).constr->arity() == 3) {
            double terntight;
            tctr1 = (TernaryConstraint*) (*iter).constr;
            terntight = tctr1->getTightness() +
                    tctr1->xy->getTightness() +
                    tctr1->xz->getTightness() +
                    tctr1->yz->getTightness();

            Variable *x1 = NULL, *x2 = NULL;
            switch((*iter).scopeIndex) {
            case 0: x1 = tctr1->getVar(1); x2 = tctr1->getVar(2); break;
            case 1: x1 = tctr1->getVar(0); x2 = tctr1->getVar(2); break;
            case 2: x1 = tctr1->getVar(0); x2 = tctr1->getVar(1); break;
            default:;
            }

            if(terntight > maxtight) { maxtight = terntight; strvar = x1; tctr1max = tctr1; tctr1max = NULL; }

            for(ConstraintList::iterator iter2=iter; iter2 != constrs.end(); ++iter2) {
                if((*iter2).constr->arity() == 3) {
                    TernaryConstraint* tctr2 = (TernaryConstraint*) (*iter2).constr;
                    Variable* commonvar = NULL;
                    if(tctr2->getIndex(x1) >= 0) commonvar = x1;
                    else if(tctr2->getIndex(x2) >= 0) commonvar = x2;

                    if(commonvar) {
                        terntight += tctr2->getTightness() +
                                tctr2->xy->getTightness() +
                                tctr2->xz->getTightness() +
                                tctr2->yz->getTightness();

                        if(tctr1->xy->getIndex(commonvar) >= 0) terntight -= tctr1->xy->getTightness();
                        else if(tctr1->xz->getIndex(commonvar) >= 0) terntight -= tctr1->xz->getTightness();
                        else if(tctr1->yz->getIndex(commonvar) >= 0) terntight -= tctr1->yz->getTightness();

                        if(terntight > maxtight) { maxtight = terntight; strvar = commonvar; tctr1max = tctr1; tctr2max = tctr2; }
                    }
                }
            }
        }
    }

    return maxtight;
}






ostream& operator<<(ostream& os, Variable &var) {
    os << var.name; // << " #" << var.dac;
    var.print(os);
    if (ToulBar2::verbose >= 3) {
        for (ConstraintList::iterator iter=var.constrs.begin(); iter != var.constrs.end(); ++iter) {
            os << " (" << (*iter).constr << "," << (*iter).scopeIndex << ")";
        }
    }
    return os;
}

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

