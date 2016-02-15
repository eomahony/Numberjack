/** \file tb2arithmetic.hpp
 *  \brief Binary arithmetic soft constraints.
 *  \warning Value and Cost must be of the same type.
 *
 */

#ifndef TB2ARITHMETIC_HPP_
#define TB2ARITHMETIC_HPP_

#include "tb2abstractconstr.hpp"
#include "tb2intervar.hpp"

#include <set>

/** semantic: \f$soft(penalty, x \in SetOfValues) : (x \in SetOfValues)?0:penalty\f$
 */
class Unary : public AbstractUnaryConstraint<IntervalVariable>
{
    set<Value> permitted;
    Cost penalty;
    StoreValue deltaValueXinf;
    StoreValue deltaValueXsup;

public:
    Unary(WCSP *wcsp, IntervalVariable *xx, Value *d, int dsize, Cost penalty, StoreStack<Value, Value> *storeValue);

    ~Unary() {}

    void propagate();

    void remove(int varIndex) {}

    void assign(int varIndex) {
        assert(connected());
        deconnect();  // Warning! deconnection has to be done before the projection
        if (permitted.find(x->getValue()) == permitted.end()) {
            projectLB(penalty);
        }
    }

    bool verify();

    double  computeTightness() {
        tight = (double) penalty * abs((int) permitted.size() - (int) x->getDomainSize()) / x->getDomainSize();
        return tight;
    }

    void print(ostream& os);
    void dump(ostream& os, bool original = true);
};

/** semantic: \f$soft(x \geq y + cst) : max( (y + cst - x \leq deltamax)?(y + cst - x):top , 0)\f$
 */
class Supxyc : public AbstractBinaryConstraint<IntervalVariable,IntervalVariable>
{
    Value cst;
    Value deltamax;
    StoreCost deltaCost;
    StoreValue deltaValueXinf;
    StoreValue deltaValueYsup;
    StoreCost deltaCostXinf;
    StoreCost deltaCostYsup;

public:
    Supxyc(WCSP *wcsp, IntervalVariable *xx, IntervalVariable *yy, Value c, Value delta,
            StoreStack<Cost, Cost> *storeCost, StoreStack<Value, Value> *storeValue);

    ~Supxyc() {}

    void propagate();

    void remove(int varIndex) {}

    void assign(int varIndex) {
        if (x->assigned() && y->assigned()) deconnect();
        propagate();
    }

    bool verify();

    double  computeTightness() { return 0; }

    void print(ostream& os);
    void dump(ostream& os, bool original = true);
};

/** semantic: \f$soft(penalty, x \geq y + csty \vee y \geq x + cstx) : (x \geq y + csty \vee y \geq x + cstx)?0:penalty\f$
 */
class Disjunction : public AbstractBinaryConstraint<IntervalVariable, IntervalVariable>
{
    Value cstx;
    Value csty;
    Cost penalty;
    StoreValue deltaValueXinf;
    StoreValue deltaValueYinf;
    StoreValue deltaValueXsup;
    StoreValue deltaValueYsup;

public:
    Disjunction(WCSP *wcsp, IntervalVariable *xx, IntervalVariable *yy, Value cxx, Value cyy,
            Cost penalty, StoreStack<Value, Value> *storeValue);

    ~Disjunction() {}

    void propagate();

    void remove(int varIndex) {}

    bool verify();

    double  computeTightness() { return 0; }

    void print(ostream& os);
    void dump(ostream& os, bool original = true);
};

/** semantic:
 * implicit hard constraint: \f$x \leq xinfty\f$
 *
 * implicit hard constraint: \f$y \leq yinfty\f$
 *
 * cost function definition:
 * - \f$(x<xinfty \wedge y<yinfty \wedge (x \geq y + csty \vee y \geq x + cstx)): cost = 0\f$ (mutual exclusion between tasks x and y)
 * - \f$x=xinfty \wedge y=yinfty: cost = costx + costy\f$ (both tasks x and y unselected)
 * - \f$x=xinfty: cost = costx\f$ (task x unselected)
 * - \f$y=yinfty: cost = costy\f$ (task y unselected)
 * - else \f$cost = UB\f$ (task x and y must be selected but cannot be scheduled)
 */
class SpecialDisjunction : public AbstractBinaryConstraint<IntervalVariable, IntervalVariable>
{
    Value cstx;
    Value csty;
    Value xinfty;
    Value yinfty;
    Cost costx;
    Cost costy;
    StoreCost deltaCost;
    StoreValue deltaValueXinf;
    StoreValue deltaValueYinf;
    StoreValue deltaValueXsup;
    StoreValue deltaValueYsup;
    StoreCost deltaCostXinf;
    StoreCost deltaCostYinf;
    StoreCost deltaCostXsup;
    StoreCost deltaCostYsup;

public:
    SpecialDisjunction(WCSP *wcsp, IntervalVariable *xx, IntervalVariable *yy, Value cxx, Value cyy, 
            Value xmax, Value ymax, Cost xcost, Cost ycost,
            StoreStack<Cost, Cost> *storeCost, StoreStack<Value, Value> *storeValue);

    ~SpecialDisjunction() {}

    void propagate();

    void remove(int varIndex) {}

    void assign(int varIndex) {
        assert(connected());
        wcsp->revise(this);
        if (x->assigned() && y->assigned()) {
            deconnect();
            if (x->getValue() >= xinfty && y->getValue() >= yinfty && costx + costy > deltaCost) {
                projectLB(costx + costy - deltaCost);
            } else if (x->getValue() >= xinfty && y->getValue() < yinfty && costx > deltaCost) {
                projectLB(costx - deltaCost);
            } else if (y->getValue() >= yinfty && x->getValue() < xinfty && costy > deltaCost) {
                projectLB(costy - deltaCost);
            } else if (x->getValue() < xinfty && y->getValue() < yinfty && x->getValue() < y->getValue() + csty && y->getValue() < x->getValue() + cstx) {
                THROWCONTRADICTION;
            }
        } else {
            if (varIndex==0 && x->getValue()>=xinfty) {
                if (y->getInf()==deltaValueYinf) {
                    Cost cost = deltaCostYinf;
                    deltaCostYinf = MIN_COST;
                    y->projectInfCost(-cost);
                }
                if (y->getSup()==deltaValueYsup) {
                    Cost cost = deltaCostYsup;
                    deltaCostYsup = MIN_COST;
                    y->projectSupCost(-cost);
                }
                Cost cost = costx - deltaCost;
                if (cost > MIN_COST) {
                    deltaCost += cost;
                    projectLB(cost);
                }
                if (y->getSup() < yinfty) {
                    deconnect();
                } else {
                    assert(y->getSup() == yinfty);
                    deltaValueYsup = yinfty;
                    deltaCostYsup = costy;
                    y->projectSupCost(costy);
                }
            } else if (varIndex==1 && y->getValue()>=yinfty) {
                if (x->getInf()==deltaValueXinf) {
                    Cost cost = deltaCostXinf;
                    deltaCostXinf = MIN_COST;
                    x->projectInfCost(-cost);
                }
                if (x->getSup()==deltaValueXsup) {
                    Cost cost = deltaCostXsup;
                    deltaCostXsup = MIN_COST;
                    x->projectSupCost(-cost);
                }
                Cost cost = costy - deltaCost;
                if (cost > MIN_COST) {
                    deltaCost += cost;
                    projectLB(cost);
                }
                if (x->getSup() < xinfty) {
                    deconnect();
                } else {
                    assert(x->getSup() == xinfty);
                    deltaValueXsup = xinfty;
                    deltaCostXsup = costx;
                    x->projectSupCost(costx);
                }
            } else {
                propagate();
            }
        }
    }

    bool verify();

    double  computeTightness() { return 0; }

    void print(ostream& os);
    void dump(ostream& os, bool original = true);
};

#endif /*TB2ARITHMETIC_HPP_*/

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

