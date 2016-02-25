/*
 * **************** Abstract constraint **************
 */

#include "tb2constraint.hpp"
#include "tb2wcsp.hpp"
#include "tb2clusters.hpp"

/*
 * Constructor
 *
 */

Constraint::Constraint(WCSP *w) : WCSPLink(w,w->numberOfConstraints()), conflictWeight(1), fromElim1(NULL), fromElim2(NULL)
{
    w->link(this);
    tight = -1;
    isSep_ = false;
    isDuplicate_ = false;
    cluster = -1;
}

Constraint::Constraint(WCSP *w, int elimCtrIndex) : WCSPLink(w,elimCtrIndex), conflictWeight(1), fromElim1(NULL), fromElim2(NULL)
{
    tight = -1;
    isSep_ = false;
    isDuplicate_ = false;
    cluster = -1;
}

/// \return size of the cartesian product of all domains in the constraint scope.
/// \warning use deprecated MAX_DOMAIN_SIZE for performance.
Long Constraint::getDomainSizeProduct()
{
    if (arity()==0) return 0;
    Long cartesianProduct = 1;
    for (int i=0; i<arity(); i++) {
        // trap overflow numbers
        if (cartesianProduct > LONGLONG_MAX / MAX_DOMAIN_SIZE) return LONGLONG_MAX;
        cartesianProduct *= getVar(i)->getDomainSize();
    }
    return cartesianProduct;
}

void Constraint::conflict()
{
    wcsp->conflict();
}

void Constraint::projectLB(Cost cost)
{
    if (cost == 0) return;
    if (ToulBar2::verbose >= 2) cout << "lower bound increased " << wcsp->getLb() << " -> " << wcsp->getLb()+cost << endl;
    wcsp->increaseLb(cost); // done before cluster LB because of #CSP (assuming a contradiction will occur here)
    if (wcsp->td) {
        if (ToulBar2::verbose >= 2) cout << " in cluster C" << getCluster() << " (from " << wcsp->td->getCluster(getCluster())->getLb() << " to " << wcsp->td->getCluster(getCluster())->getLb() + cost << ")" << endl;
        wcsp->td->getCluster(getCluster())->increaseLb(cost);
    }
}

void Constraint::assigns()
{
    for(int i=0;connected() && i<arity();i++) if (getVar(i)->assigned()) assign(i);
}

void Constraint::sumScopeIncluded( Constraint* ctr )
{
    int ar = arity();
    EnumeratedVariable** scopethis = new EnumeratedVariable * [arity()];
    for(int i=0;i<ar;i++) scopethis[i] = (EnumeratedVariable*) getVar(i);

    Cost Top = wcsp->getUb();
    Cost c;
    String t;

    if(getDefCost() < Top) {      // enumeration case
        firstlex();
        while( nextlex(t,c) ) {
            Cost cplus = ctr->evalsubstr(t, this);
            if(c + cplus < Top) {
                if (arity() > 3 && getDefCost() > MIN_COST) setTuple( t, c + cplus, scopethis);
                else addtoTuple( t, cplus, scopethis);
            } else {
                setTuple( t, Top, scopethis);
            }
        }
    } else {
        first();
        while( next(t,c) ) {
            Cost cplus = ctr->evalsubstr(t, this);
            if(c + cplus < Top) addtoTuple( t, cplus, scopethis);
            else setTuple( t, Top, scopethis);
        }
    }

    delete [] scopethis;
}


void Constraint::assignCluster() {
    TreeDecomposition* td = wcsp->getTreeDec();
    if(!td) return;
    Cluster* lowest = td->getRoot();
    for(int i=0;i<arity();i++) if (getVar(i)->unassigned()) {
        Variable* x = getVar(i);
        Cluster* c = td->getCluster( x->getCluster() );
        if(lowest->isDescendant(c)) lowest = c;
    }
    cluster = lowest->getId();
}

/// \warning always returns 0 for cost functions in intention
Cost Constraint::getMinCost()
{
    if (!extension()) return MIN_COST;

    // 	    Cost minc = MAX_COST;
    //         String tuple;
    //         Cost cost;
    //         firstlex();
    //         while (nextlex(tuple,cost)) {
    //             if (cost < minc) minc = cost;
    //         }
    //         return minc;

    Cost minc = MAX_COST;
    String tuple;
    Cost cost;
    Long nbtuples = 0;
    first();
    while (next(tuple,cost)) {
        nbtuples++;
        if (cost < minc) minc = cost;
    }
    if (getDefCost() < minc && nbtuples < getDomainSizeProduct()) minc = getDefCost();
    return minc;
}

/// \warning always returns false for cost functions in intention
bool Constraint::universal()
{
    if (!extension()) return false;

    //   String tuple;
    //   Cost cost;
    //   firstlex();
    //   while (nextlex(tuple,cost)) {
    // 	if (cost > MIN_COST) return false;
    //   }
    //   return true;

    String tuple;
    Cost cost;
    Long nbtuples = 0;
    first();
    while (next(tuple,cost)) {
        nbtuples++;
        if (cost > MIN_COST) return false;
    }
    if (getDefCost() > MIN_COST && nbtuples < getDomainSizeProduct()) return false;
    return true;
}

/// \warning always returns MAX_COST for cost functions in intention
Cost Constraint::getMaxFiniteCost()
{
    if (!extension()) return MAX_COST;

    Cost maxcost = MIN_COST;
    String tuple;
    Cost cost;
    Long nbtuples = 0;
    first();
    while (next(tuple,cost)) {
        nbtuples++;
        if (cost < wcsp->getUb() && cost > maxcost) maxcost = cost;
    }
    if (getDefCost() < wcsp->getUb() && getDefCost() > maxcost && nbtuples < getDomainSizeProduct()) maxcost = getDefCost();
    return maxcost;
}

/// \warning always returns false for cost functions in intention
bool Constraint::ishard()
{
    if (!extension()) return false;

    String tuple;
    Cost cost;
    firstlex();
    while (nextlex(tuple,cost)) {
        if (cost > MIN_COST && !CUT(cost,wcsp->getUb())) return false;
    }
    return true;
}

bool Constraint::verifySeparate(Constraint * ctr1, Constraint * ctr2)
{
    assert(scopeIncluded(ctr1));
    assert(scopeIncluded(ctr2));
    String tuple;
    Cost cost, c1,c2;
    firstlex();
    if(ToulBar2::verbose >= 3) { cout << "[ ";
    for(int i = 0; i< arity(); ++i)
        cout << getVar(i)->getName() << " ";
    cout << " ]\n" ;
    }
    while(nextlex(tuple,cost)){
        c1 = ctr1->evalsubstr(tuple,this);
        c2 = ctr2->evalsubstr(tuple,this);
        if(ToulBar2::verbose >= 3){
            for(int i = 0; i< arity(); ++i)
                cout << tuple[i] - CHAR_FIRST << " ";
            //cout << endl;
            cout << " : " << cost << " =? " << c1 << " + " << c2 << " : " << c1 + c2 << endl;
        }
        if(cost < wcsp->getUb() && c1+c2 != cost) return false;
        if(cost >= wcsp->getUb() && c1+c2 < wcsp->getUb()) return false;
    }
    return true;
}

bool Constraint::decompose()
{
    bool sep = false;
    if(extension() && !universal() && (arity() == 3 || (arity()>=4 && (getDefCost()>MIN_COST || ((NaryConstraint *) this)->size()>1)))) {
        TSCOPE scopeinv;
        getScope(scopeinv);
        EnumeratedVariable * vx = NULL;
        EnumeratedVariable * vz = NULL;
        for(TSCOPE::reverse_iterator it1 = scopeinv.rbegin(); it1 != scopeinv.rend() && !sep; ++it1) {
            TSCOPE::reverse_iterator it2 = it1;
            for(++it2; it2 != scopeinv.rend() && !sep; ++it2) {
                vx = (EnumeratedVariable *) wcsp->getVar((*it2).first);
                vz = (EnumeratedVariable *) wcsp->getVar((*it1).first);
                if(ToulBar2::verbose >= 1) cout << /*"\n" <<*/ vx->getName() << " and " << vz->getName() << " are separable in ";
                sep = separability(vx,vz);
                if(sep && ToulBar2::verbose >= 1) {
                    cout << " YES";
#ifndef NDEBUG
                    if (!ishard()) cout << " with finite costs";
#endif
                    cout << endl;
                }
                if(!sep && ToulBar2::verbose >= 1) cout << " NO" << endl;
            }
        }
        if(sep) separate(vx,vz);
        if(ToulBar2::verbose >= 3) cout << "=====================================================" << endl;
    }
    return sep;
}

Constraint *Constraint::copy() 
{
    int scope[arity()];
    for (int i=0; i<arity(); i++) scope[i] = getVar(i)->wcspIndex;
    int ctrIndex = wcsp->postNaryConstraintBegin(scope, arity(), getDefCost());
    Cost c;
    String t;
    first();
    while(next(t,c)) {
        wcsp->postNaryConstraintTuple(ctrIndex, t, c);
    }
    wcsp->getCtr(ctrIndex)->deconnect();
    return wcsp->getCtr(ctrIndex);
}

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

