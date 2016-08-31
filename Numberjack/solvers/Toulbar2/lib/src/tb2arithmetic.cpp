/*
 * ****** Binary arithmetic soft constraints ******
 */

#include "tb2arithmetic.hpp"
#include "tb2wcsp.hpp"


/*
 * Constructors and misc.
 * 
 */

Unary::Unary(WCSP *wcsp, IntervalVariable *xx, Value *d, int dsize, Cost c) :
        AbstractUnaryConstraint<IntervalVariable>(wcsp, xx), permitted(d, d+dsize), penalty(c),
        deltaValueXinf(xx->getSup()+1), deltaValueXsup(xx->getSup()+1)
{
    xx->queueInc();
    xx->queueDec();
}

void Unary::print(ostream& os)
{
    os << this << " " << x->getName() << " var in {";
    for (set<Value>::iterator it=permitted.begin(); it != permitted.end(); ++it) {
        os << " " << *it;
    }
    os  << " } (" << penalty << ")" << endl;
}

void Unary::dump(ostream& os, bool original)
{
    os << "1 " << ((original)?(x->wcspIndex):x->getCurrentVarId()) << " " << penalty << " " << permitted.size() << endl;
    for (set<Value>::iterator it=permitted.begin(); it != permitted.end(); ++it) {
        os << *it << " 0" << endl;
    }
}

Supxyc::Supxyc(WCSP *wcsp, IntervalVariable *xx, IntervalVariable *yy, Value c, Value delta) :
        AbstractBinaryConstraint<IntervalVariable,IntervalVariable>(wcsp, xx, yy),
        cst(c), deltamax(delta), deltaCost(MIN_COST),
        deltaValueXinf(xx->getSup()+1), deltaValueYsup(yy->getSup()+1),
        deltaCostXinf(MIN_COST), deltaCostYsup(MIN_COST)
{
    xx->queueInc();
    xx->queueDec();
    yy->queueInc();
    yy->queueDec();
}

void Supxyc::print(ostream& os)
{
    os << this << " " << x->getName() << " >= " << y->getName() << " + " << cst << " (" << deltamax << ")" << endl;
}

void Supxyc::dump(ostream& os, bool original)
{
    os << "2 " << ((original)?(x->wcspIndex):x->getCurrentVarId()) << " " << ((original)?(y->wcspIndex):y->getCurrentVarId()) << " -1 >= " << cst << " " << deltamax << endl;
}

Disjunction::Disjunction(WCSP *wcsp, IntervalVariable *xx, IntervalVariable *yy, Value cxx, Value cyy,
        Cost cost) :
        AbstractBinaryConstraint<IntervalVariable,IntervalVariable>(wcsp, xx, yy),
        cstx(cxx), csty(cyy), penalty(cost),
        deltaValueXinf(xx->getSup()+1),deltaValueYinf(yy->getSup()+1),
        deltaValueXsup(xx->getSup()+1),deltaValueYsup(yy->getSup()+1)
{
    xx->queueInc();
    xx->queueDec();
    yy->queueInc();
    yy->queueDec();
}

void Disjunction::print(ostream& os)
{
    os << this << " " << x->getName() << " >= " << y->getName() << " + " << csty << " or " << y->getName() << " >= " << x->getName() << " + " << cstx << " ("  << penalty<< ")" << endl;
}

void Disjunction::dump(ostream& os, bool original)
{
    os << "2 " << ((original)?(x->wcspIndex):x->getCurrentVarId()) << " " << ((original)?(y->wcspIndex):y->getCurrentVarId()) << " -1 disj " << cstx << " " << csty << " " << penalty << endl;
}

SpecialDisjunction::SpecialDisjunction(WCSP *wcsp, IntervalVariable *xx, IntervalVariable *yy, 
        Value cxx, Value cyy, Value xmax, Value ymax,
        Cost xcost, Cost ycost) :
        AbstractBinaryConstraint<IntervalVariable,IntervalVariable>(wcsp, xx, yy),
        cstx(cxx), csty(cyy), xinfty(xmax), yinfty(ymax), costx(xcost), costy(ycost),
        deltaCost(MIN_COST),
        deltaValueXinf(xx->getSup()+1),deltaValueYinf(yy->getSup()+1),
        deltaValueXsup(xx->getSup()+1),deltaValueYsup(yy->getSup()+1),
        deltaCostXinf(MIN_COST), deltaCostYinf(MIN_COST),
        deltaCostXsup(MIN_COST), deltaCostYsup(MIN_COST)
{
    xx->queueInc();
    xx->queueDec();
    yy->queueInc();
    yy->queueDec();
}

void SpecialDisjunction::print(ostream& os)
{
    os << this << " " << x->getName() << " < " << xinfty << " and " << y->getName() << " < " << yinfty << " and (" << x->getName() << " >= " << y->getName() << " + " << csty << " or " << y->getName() << " >= " << x->getName() << " + " << cstx << ") ("  << costx << "," << costy << ")" << endl;
}

void SpecialDisjunction::dump(ostream& os, bool original)
{
    os << "2 " << ((original)?(x->wcspIndex):x->getCurrentVarId()) << " " << ((original)?(y->wcspIndex):y->getCurrentVarId()) << " -1 sdisj " << cstx << " " << csty << " " << xinfty << " " << yinfty << " " << costx << " " << costy << endl;
}

/*
 * Propagation methods
 * 
 */

void Unary::propagate()
{
    if (ToulBar2::verbose >= 3) {
        print(cout);
        cout << " dxinf=" << deltaValueXinf << " dxsup=" << deltaValueXsup << endl;
    }
    wcsp->revise(this);
    set<Value>::iterator itinf = permitted.lower_bound(x->getInf());
    set<Value>::iterator itsup = permitted.upper_bound(x->getSup());
    --itsup;
    if (itinf == permitted.end() || itsup == permitted.end()) {
        // IC0 propagatation (increase global lower bound)
        deconnect();
        projectLB(penalty);
    } else {
        // propagate hard constraint
        if (CUT(wcsp->getLb()+penalty, wcsp->getUb())) {
            if (x->getInf() < *itinf) x->increase(*itinf);
            if (x->getSup() > *itsup) x->decrease(*itsup);
        }
        // BAC* propagation (increase unary costs of domain bounds)
        Value xinf = x->getInf();
        if (xinf != deltaValueXinf && xinf != deltaValueXsup && permitted.find(xinf) == permitted.end()) {
            deltaValueXinf = xinf;
            x->projectInfCost(penalty);
        }
        Value xsup = x->getSup();
        if (xsup != deltaValueXinf && xsup != deltaValueXsup && permitted.find(xsup) == permitted.end()) {
            deltaValueXsup = xsup;
            x->projectSupCost(penalty);
        }
    }
}

bool Unary::verify()
{
    bool support = (permitted.lower_bound(x->getInf()) != permitted.end() || (--permitted.upper_bound(x->getSup())) != permitted.end());
    support = support && (permitted.find(x->getInf()) != permitted.end() || deltaValueXinf == x->getInf());
    support = support && (permitted.find(x->getSup()) != permitted.end() || deltaValueXsup == x->getSup());
    if (!support) {
        print(cout);
        x->print(cout);
        cout << endl;
        cout << " dxinf=" << deltaValueXinf << " dxsup=" << deltaValueXsup << endl;
    }
    return support;
}

void Supxyc::propagate()
{
    if (ToulBar2::verbose >= 3) {
        print(cout);
        cout << " delta=" << deltaCost << " dxinf=" << deltaValueXinf << " dxcost=" << deltaCostXinf << " dysup=" << deltaValueYsup << " dycost=" << deltaCostYsup << endl;
    }
    wcsp->revise(this);
    // deconnect the constraint if always satisfied
    if (x->getInf() >= y->getSup() + cst) {
        deconnect();
    } else {
        // propagate hard constraint
        Cost gap = wcsp->getUb() - wcsp->getLb() - 1 + deltaCost;
        Value newInf = ceil(y->getInf() + cst - ((gap < deltamax)?gap:deltamax));
        if (x->getInf() < newInf) x->increase(newInf);

        Value newSup = floor(x->getSup() - cst + ((gap < deltamax)?gap:deltamax));
        if (y->getSup() > newSup) y->decrease(newSup);

        // IC0 propagatation (increase global lower bound)
        Cost cost = y->getInf() + cst - x->getSup() - deltaCost;
        if (cost > MIN_COST) {
            deltaCost += cost;
            projectLB(cost);
        }

        // BAC* propagation (increase unary costs of domain bounds for unassigned variables)
        if (x->unassigned()) {
            Value xinf = x->getInf();
            Value yinf = y->getInf();
            cost = max(yinf + cst - xinf - deltaCost, MIN_COST);
            if (xinf == deltaValueXinf) {
                Cost delta = cost - deltaCostXinf;
                if (delta != MIN_COST) {
                    deltaCostXinf = cost;
                    x->projectInfCost(delta);
                }
            } else {
                deltaValueXinf = xinf;
                deltaCostXinf = cost;
                x->projectInfCost(cost);
            }
        }

        if (y->unassigned()) {
            Value xsup = x->getSup();
            Value ysup = y->getSup();
            cost = max(ysup + cst - xsup - deltaCost, MIN_COST);
            if (ysup == deltaValueYsup) {
                Cost delta = cost - deltaCostYsup;
                if (delta != MIN_COST) {
                    deltaCostYsup = cost;
                    y->projectSupCost(delta);
                }
            } else {
                deltaValueYsup = ysup;
                deltaCostYsup = cost;
                y->projectSupCost(cost);
            }
        }
    }
}

bool Supxyc::verify()
{
    Cost cmin=MIN_COST;
    Cost cxinf=MIN_COST;
    Cost cysup=MIN_COST;

    cmin = y->getInf() + cst - x->getSup() - deltaCost;
    if (cmin > MIN_COST) cout << "cmin=" << cmin << endl;
    cxinf = y->getInf() + cst - x->getInf() - deltaCost
        - ((y->getInf() == deltaValueYsup)?(Cost)deltaCostYsup:MIN_COST)
        - ((x->getInf() == deltaValueXinf)?(Cost)deltaCostXinf:MIN_COST);
    if (cxinf > MIN_COST) cout << "cxinf=" << cxinf << endl;
    cysup = y->getSup() + cst - x->getSup() - deltaCost
        - ((y->getSup() == deltaValueYsup)?(Cost)deltaCostYsup:MIN_COST)
        - ((x->getSup() == deltaValueXinf)?(Cost)deltaCostXinf:MIN_COST);
    if (cysup > MIN_COST) cout << "cysup=" << cysup << endl;
    bool icbac = (cmin <= MIN_COST) && (cysup <= MIN_COST) && (cxinf <= MIN_COST);
    if (!icbac) {
        print(cout);
        x->print(cout);
        cout << endl;
        y->print(cout);
        cout << endl;
        cout << " delta=" << deltaCost << " dxinf=" << deltaValueXinf << " dxcost=" << deltaCostXinf << " dysup=" << deltaValueYsup << " dycost=" << deltaCostYsup << endl;
    }
    return icbac;
}

void Disjunction::propagate()
{
    if (ToulBar2::verbose >= 3) {
        print(cout);
        cout << " dxinf=" << deltaValueXinf << " dxsup=" << deltaValueXsup << " dyinf=" << deltaValueYinf << " dysup=" << deltaValueYsup << endl;
    }
    wcsp->revise(this);
    // deconnect the constraint if always satisfied
    if (x->getInf() >= y->getSup() + csty || y->getInf() >= x->getSup() + cstx) {
        deconnect();
    } else if (x->getSup() < y->getInf() + csty && y->getSup() < x->getInf() + cstx) {
        // IC0 propagatation (increase global lower bound if always unsatisfied)
        deconnect();
        if (x->unassigned()) {
            if (x->getInf() == deltaValueXinf) x->projectInfCost(-penalty);
            if (x->getSup() == deltaValueXsup) x->projectSupCost(-penalty);
        }
        if (y->unassigned()) {
            if (y->getInf() == deltaValueYinf) y->projectInfCost(-penalty);
            if (y->getSup() == deltaValueYsup) y->projectSupCost(-penalty);
        }
        projectLB(penalty);
    } else {
        // propagate hard constraint
        if (CUT(wcsp->getLb()+penalty, wcsp->getUb())) {
            if (x->getSup() < y->getInf() + csty) {
                Value newInf = x->getInf() + cstx;
                if (y->getInf() < newInf) y->increase(newInf);
                Value newSup = y->getSup() - cstx;
                if (x->getSup() > newSup) x->decrease(newSup);
            } else if (y->getSup() < x->getInf() + cstx) {
                Value newInf = y->getInf() + csty;
                if (x->getInf() < newInf) x->increase(newInf);
                Value newSup = x->getSup() - csty;
                if (y->getSup() > newSup) y->decrease(newSup);
            }
        }

        // BAC* propagation (increase unary costs of domain bounds)
        if (x->unassigned()) {
            Value xinf = x->getInf();
            if (xinf != deltaValueXinf && xinf != deltaValueXsup && xinf < y->getInf() + csty && xinf > y->getSup() - cstx) {
                deltaValueXinf = xinf;
                x->projectInfCost(penalty);
            }
            Value xsup = x->getSup();
            if (xsup != deltaValueXsup && xsup != deltaValueXinf && xsup < y->getInf() + csty && xsup > y->getSup() - cstx) {
                deltaValueXsup = xsup;
                x->projectSupCost(penalty);
            }
        }
        if (y->unassigned()) {
            Value yinf = y->getInf();
            if (yinf != deltaValueYinf && yinf != deltaValueYsup && yinf < x->getInf() + cstx && yinf > x->getSup() - csty) {
                deltaValueYinf = yinf;
                y->projectInfCost(penalty);
            }
            Value ysup = y->getSup();
            if (ysup != deltaValueYsup && ysup != deltaValueYinf && ysup < x->getInf() + cstx && ysup > x->getSup() - csty) {
                deltaValueYsup = ysup;
                y->projectSupCost(penalty);
            }
        }
    }
}

bool Disjunction::verify()
{
    bool support = (x->getSup() >= y->getInf() + csty || y->getSup() >= x->getInf() + cstx);
    Value xinf = x->getInf();
    Value yinf = y->getInf();
    Value xsup = x->getSup();
    Value ysup = y->getSup();
    support = support && (xinf >= yinf + csty || xinf <= ysup - cstx || deltaValueXinf == xinf);
    support = support && (yinf >= xinf + cstx || yinf <= xsup - csty || deltaValueYinf == yinf);
    support = support && (xsup >= yinf + csty || xsup <= ysup - cstx || deltaValueXsup == xsup);
    support = support && (ysup >= xinf + cstx || ysup <= xsup - csty || deltaValueYsup == ysup);
    if (!support) {
        print(cout);
        x->print(cout);
        cout << endl;
        y->print(cout);
        cout << endl;
        cout << " dxinf=" << deltaValueXinf << " dxsup=" << deltaValueXsup << " dyinf=" << deltaValueYinf << " dysup=" << deltaValueYsup << endl;
    }
    return support;
}

void SpecialDisjunction::propagate()
{
    if (ToulBar2::verbose >= 3) {
        print(cout);
        cout << "deltaCost= " << deltaCost << " dxinf=" << deltaValueXinf << " dxsup=" << deltaValueXsup << " dyinf=" << deltaValueYinf << " dysup=" << deltaValueYsup << endl;
    }
    wcsp->revise(this);
    if (x->getSup()>xinfty) x->decrease(xinfty);
    if (y->getSup()>yinfty) y->decrease(yinfty);
    if ((x->getSup() < xinfty && y->getSup() < yinfty &&
            (x->getInf() >= y->getSup() + csty || y->getInf() >= x->getSup() + cstx)) ||
            (x->getInf() >= xinfty && y->getSup() < yinfty) ||
            (y->getInf() >= yinfty && x->getSup() < xinfty)) {
        // deconnect the constraint if always satisfied
        assert(x->getInf() < xinfty || y->getSup() >= yinfty || deltaCost == costx);
        assert(y->getInf() < yinfty || x->getSup() >= xinfty || deltaCost == costy);
        deconnect();
    } else if (x->getSup() < xinfty && y->getSup() < yinfty &&
            x->getSup() < y->getInf() + csty && y->getSup() < x->getInf() + cstx) {
        // backtrack if the constraint if always unsatisfied
        THROWCONTRADICTION;
    } else {
        if (x->getInf() < xinfty && y->getInf() < yinfty) {
            // IC0 propagatation (increase global lower bound)
            if (min(x->getSup(),xinfty-1) < y->getInf() + csty &&
                    min(y->getSup(),yinfty-1) < x->getInf() + cstx) {
                Cost cost = min(costx,costy) - deltaCost;
                if (cost > MIN_COST) {
                    deltaCost += cost;
                    projectLB(cost);
                }
            }
            // propagate hard constraint
            if ((x->getSup() < xinfty || CUT(wcsp->getLb()+costx-deltaCost, wcsp->getUb())) &&
                    min(x->getSup(),xinfty-1) < y->getInf() + csty) {
                Value newInf = min(x->getInf() + cstx, yinfty);
                if (y->getInf() < newInf) y->increase(newInf);
                if (y->getSup() < yinfty) {
                    Value newSup = y->getSup() - cstx;
                    if (x->getSup() > newSup) x->decrease(newSup);
                }
            }
            if ((y->getSup() < yinfty || CUT(wcsp->getLb()+costy-deltaCost, wcsp->getUb())) &&
                    min(y->getSup(),yinfty-1) < x->getInf() + cstx) {
                Value newInf = min(y->getInf() + csty, xinfty);
                if (x->getInf() < newInf) x->increase(newInf);
                if (x->getSup() < xinfty) {
                    Value newSup = x->getSup() - csty;
                    if (y->getSup() > newSup) y->decrease(newSup);
                }
            }
            // BAC* propagation (increase unary costs of domain bounds)
            if (x->unassigned() && y->getInf() < yinfty) {
                Value xinf = x->getInf();
                Cost cost = -((Cost) deltaCost);
                if (xinf < y->getInf() + csty && xinf > min(y->getSup(),yinfty-1) - cstx) cost += costy;
                if (xinf == deltaValueXinf) {
                    Cost delta = cost - deltaCostXinf;
                    if (delta != MIN_COST) {
                        deltaCostXinf = cost;
                        x->projectInfCost(delta);
                    }
                } else {
                    deltaValueXinf = xinf;
                    deltaCostXinf = cost;
                    x->projectInfCost(cost);
                }
            }
            if (x->unassigned() && y->getInf() < yinfty) {
                Value xsup = x->getSup();
                Cost cost = -((Cost) deltaCost);
                if (xsup==xinfty) cost += costx;
                else if (xsup < y->getInf() + csty && xsup > min(y->getSup(),yinfty-1) - cstx) cost += costy;
                if (xsup == deltaValueXsup) {
                    Cost delta = cost - deltaCostXsup;
                    if (delta != MIN_COST) {
                        deltaCostXsup = cost;
                        x->projectSupCost(delta);
                    }
                } else {
                    deltaValueXsup = xsup;
                    deltaCostXsup = cost;
                    x->projectSupCost(cost);
                }
            }
            if (y->unassigned() && x->getInf() < xinfty) {
                Value yinf = y->getInf();
                Cost cost = -((Cost) deltaCost);
                if (yinf < x->getInf() + cstx && yinf > min(x->getSup(),xinfty-1) - csty) cost += costx;
                if (yinf == deltaValueYinf) {
                    Cost delta = cost - deltaCostYinf;
                    if (delta != MIN_COST) {
                        deltaCostYinf = cost;
                        y->projectInfCost(delta);
                    }
                } else {
                    deltaValueYinf = yinf;
                    deltaCostYinf = cost;
                    y->projectInfCost(cost);
                }
            }
            if (y->unassigned() && x->getInf() < xinfty) {
                Value ysup = y->getSup();
                Cost cost = -((Cost) deltaCost);
                if (ysup==yinfty) cost += costy;
                else if (ysup < x->getInf() + cstx && ysup > min(x->getSup(),xinfty-1) - csty) cost += costx;
                if (ysup == deltaValueYsup) {
                    Cost delta = cost - deltaCostYsup;
                    if (delta != MIN_COST) {
                        deltaCostYsup = cost;
                        y->projectSupCost(delta);
                    }
                } else {
                    deltaValueYsup = ysup;
                    deltaCostYsup = cost;
                    y->projectSupCost(cost);
                }
            }
        }
    }
}

bool SpecialDisjunction::verify()
{
    bool support = (x->getSup()>=xinfty || y->getSup()>=yinfty || x->getSup() >= y->getInf() + csty || y->getSup() >= x->getInf() + cstx);
    Value xinf = x->getInf();
    Value yinf = y->getInf();
    Value xsup = x->getSup();
    Value ysup = y->getSup();
    support = support && (xinf >= xinfty || yinf >= yinfty || xinf >= yinf + csty || xinf <= ysup - cstx || deltaValueXinf == xinf);
    support = support && (xinf >= xinfty || yinf >= yinfty || yinf >= xinf + cstx || yinf <= xsup - csty || deltaValueYinf == yinf);
    support = support && (xsup >= xinfty || xsup >= yinf + csty || xsup <= ysup - cstx || deltaValueXsup == xsup);
    support = support && (ysup >= yinfty || ysup >= xinf + cstx || ysup <= xsup - csty || deltaValueYsup == ysup);
    if (!support) {
        print(cout);
        x->print(cout);
        cout << endl;
        y->print(cout);
        cout << endl;
        cout << "deltaCost= " << deltaCost << " dxinf=" << deltaValueXinf << " dxsup=" << deltaValueXsup << " dyinf=" << deltaValueYinf << " dysup=" << deltaValueYsup << endl;
    }
    return support;
}

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

