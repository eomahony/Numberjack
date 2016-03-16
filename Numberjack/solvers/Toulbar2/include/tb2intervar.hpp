/** \file tb2variable.hpp
 *  \brief Variable with domain represented by an interval.
 * 
 */

#ifndef TB2INTERVAR_HPP_
#define TB2INTERVAR_HPP_

#include "tb2variable.hpp"

class IntervalVariable : public Variable
{
    StoreCost infCost;
    StoreCost supCost;

    void increaseFast(Value newInf);        // Do not insert in NC queue
    void decreaseFast(Value newSup);        // Do not insert in NC queue

public:    
    IntervalVariable(WCSP *wcsp, string n, Value iinf, Value isup);

    bool enumerated() const {return false;}

    unsigned int getDomainSize() const {
        return sup - inf + 1;
    }

    void increase(Value newInf, bool isDecision = false);
    void decrease(Value newSup, bool isDecision = false);
    void remove(Value newValue, bool isDecision = false) {if (newValue==inf) increase(newValue+1, isDecision); else if (newValue==sup) decrease(newValue-1, isDecision);}
    void assign(Value newValue, bool isDecision = false);
    void assignLS(Value newValue, set<Constraint *>& delayedCtrs);

    Cost getInfCost() const {return infCost;}
    Cost getSupCost() const {return supCost;}
    void projectInfCost(Cost cost);
    void projectSupCost(Cost cost);

    // this method can be applied to interval or enumerated domain
    Cost getCost(const Value value) const {
        if (value == inf) return getInfCost();
        else if (value == sup) return getSupCost();
        else return MIN_COST;
    }

    void propagateNC();    
    bool verifyNC();

    class iterator;
    friend class iterator;
    class iterator {    // : public Variable::iterator {
        IntervalVariable *var;
        Value value;
    public:
        iterator(IntervalVariable *v, Value vv) : var(v), value(vv) {}

        Value operator*() const {return value;}

        inline iterator &operator++() {    // Prefix form
            if (value < var->sup) ++value;
            else value = var->sup + 1;
            return *this;
        }

        iterator &operator--() {    // Prefix form
            if (value > var->inf) --value;
            else value = var->sup + 1;
            return *this;
        }

        // To see if you're at the end:
        bool operator==(const iterator &iter) const {return value == iter.value;}
        bool operator!=(const iterator &iter) const {return value != iter.value;}
    };
    iterator begin() {
        return iterator(this, inf);
    }
    iterator end() {
        return iterator(this, sup + 1);
    }
    iterator rbegin() {
        return iterator(this, sup);
    }
    iterator rend() {return end();}

    //Finds the first available element whose value is greater or equal to v
    iterator lower_bound(Value v) {
        if (v <= sup) return iterator(this, max(getInf(), v));
        else return end();
    }

    //Finds the first available element whose value is lower or equal to v
    iterator upper_bound(Value v) {
        if (v >= inf) return iterator(this, min(getSup(), v));
        else return end();
    }

    void print(ostream& os);
};

#endif /*TB2INTERVAR_HPP_*/

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

