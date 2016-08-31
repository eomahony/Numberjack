/** \file tb2enumvar.hpp
 *  \brief Variable with domain represented by an enumerated domain.
 * 
 */

#ifndef TB2ENUMVAR_HPP_
#define TB2ENUMVAR_HPP_

#include "tb2variable.hpp"
#include "tb2domain.hpp"

class EnumeratedVariable : public Variable
{
protected:
    Domain domain;
    vector<StoreCost> costs;
    StoreCost deltaCost;
    StoreValue support;     // Warning! the unary support has to be backtrackable 

    DLink<VariableWithTimeStamp> linkACQueue;
    DLink<VariableWithTimeStamp> linkDACQueue;
    DLink<VariableWithTimeStamp> linkEAC1Queue;
    DLink<VariableWithTimeStamp> linkEAC2Queue;
    DLink<VariableWithTimeStamp> linkDEEQueue;

    bool watchForIncrease;	///< \warning should be true if there exists a cost function on this variable watching for increase events
    bool watchForDecrease;	///< \warning should be true if there exists a cost function on this variable watching for decrease events
    ConstraintLink DEE;		///< \brief residue for dead-end elimination
    vector<ConstraintLink> DEE2;		///< \brief residue for generalized dead-end elimination

    void init();

    virtual void increaseFast(Value newInf);        // Do not check for a support nor insert in NC and DAC queue
    virtual void decreaseFast(Value newSup);        // Do not check for a support nor insert in NC and DAC queue
    virtual void removeFast(Value val);             // Do not check for a support nor insert in NC and DAC queue

public:    
    EnumeratedVariable(WCSP *wcsp, string n, Value iinf, Value isup);
    EnumeratedVariable(WCSP *wcsp, string n, Value *d, int dsize);

    bool enumerated() const FINAL {return true;}

    unsigned int getDomainInitSize() const {return domain.getInitSize();}
#if defined(WCSPFORMATONLY) && !defined(NUMBERJACK)
    unsigned int toIndex(Value v) const {return (unsigned int) v;}
    Value toValue(unsigned int idx) const {return idx;}
#else
    unsigned int toIndex(Value v) const {return domain.toIndex(v);}
    Value toValue(unsigned int idx) const {return domain.toValue(idx);}
#endif
    unsigned int toCurrentIndex(Value v) {return domain.toCurrentIndex(v);} // return value position in current domain
    unsigned int getDomainSize() const FINAL {
        if (assigned()) return 1; 
        else return domain.getSize(); ///< \warning can return a negative size in the case of a wrong list utilization
    }
    void getDomain(Value *array);
    void getDomainAndCost(ValueCost *array);

    bool canbe(Value v) const FINAL {return v >= inf && v <= sup && domain.canbe(v);}
    bool canbeAfterElim(Value v) const {return domain.canbe(v);}
    bool cannotbe(Value v) const FINAL {return v < inf || v > sup || domain.cannotbe(v);}

    virtual void increase(Value newInf, bool isDecision = false);
    virtual void decrease(Value newSup, bool isDecision = false);
    virtual void remove(Value value, bool isDecision = false);
    virtual void assign(Value newValue, bool isDecision = false);
    void assignWhenEliminated(Value newValue);
    void assignLS(Value newValue, ConstraintSet& delayedCtrs);

    virtual void project(Value value, Cost cost, bool delayed = false); ///< \param delayed if true, it does not check for forbidden cost/value and let node consistency do the job later
    virtual void extend(Value value, Cost cost);
    virtual void extendAll(Cost cost);
    Value getSupport() const FINAL {return support;}
    void setSupport(Value val) {support = val;}    
    inline Cost getCost(const Value value) const FINAL {
        return costs[toIndex(value)] - deltaCost;
    }
    Cost getBinaryCost(ConstraintLink c,    Value myvalue, Value itsvalue);
    Cost getBinaryCost(BinaryConstraint* c, Value myvalue, Value itsvalue);

    Cost getInfCost() const FINAL {return costs[toIndex(getInf())] - deltaCost;}
    Cost getSupCost() const FINAL {return costs[toIndex(getSup())] - deltaCost;}
    void projectInfCost(Cost cost);
    void projectSupCost(Cost cost);

    void propagateNC();    
    bool verifyNC();
    void queueAC();                     // public method used also by tb2binconstr.hpp
    void queueDAC();
    void propagateAC();
    void propagateDAC();
    void findSupport();
    bool verify();

    void queueEAC1();
    void queueEAC2();
    void fillEAC2(bool self);
    bool isEAC(Value a);
    bool isEAC();
    void propagateEAC();
    void setCostProvidingPartition();
    virtual void queueVAC2() {}

    void eliminate();
    bool elimVar( BinaryConstraint* xy );
    bool elimVar( ConstraintLink xylink,  ConstraintLink xzlink );
    bool elimVar( TernaryConstraint* xyz );

    void queueDEE();
    void propagateDEE(Value a, Value b, bool dee = true);
    bool verifyDEE(Value a, Value b);
    bool verifyDEE();

    // merge current cost functions to x's list by replacing current variable y by x thanks to functional constraint xy (i.e., y := functional[x])
    void mergeTo( BinaryConstraint *xy, map<Value, Value> &functional);
    bool canbeMerged(EnumeratedVariable *x);

    class iterator;
    friend class iterator;
    class iterator {
        EnumeratedVariable *var;
        Domain::iterator diter;
    public:
        iterator() { var = NULL; }
        iterator(EnumeratedVariable *v, Domain::iterator iter) : var(v), diter(iter) {}

        Value operator*() const {return *diter;}

        iterator &operator++() {    // Prefix form //TODO: add a const_iterator to speed-up iterations (should be inlined?)
            if (var->unassigned()) ++diter;
            else {
                if (*diter < var->getValue()) diter = var->domain.lower_bound(var->getValue());
                else diter = var->domain.end();
            }
            return *this;
        }

        iterator &operator--() {    // Prefix form
            if (var->unassigned()) --diter;
            else {
                if (*diter > var->getValue()) diter = var->domain.lower_bound(var->getValue());
                else diter = var->domain.end();
            }
            return *this;
        }

        // To see if you're at the end:
        bool operator==(const iterator &iter) const {return diter == iter.diter;}
        bool operator!=(const iterator &iter) const {return diter != iter.diter;}
    };
    iterator begin() {
        if (assigned()) return iterator(this, domain.lower_bound(getValue()));
        else return iterator(this, domain.begin());
    }
    iterator end() {return iterator(this, domain.end());}
    iterator rbegin() {
        if (assigned()) return iterator(this, domain.upper_bound(getValue()));
        else return iterator(this, domain.rbegin());
    }
    iterator rend() {return end();}

    //Finds the first available element whose value is greater or equal to v
    iterator lower_bound(Value v) {
        if (assigned()) {
            if (v <= getValue()) return iterator(this, domain.lower_bound(getValue()));
            else return end();
        } else if (v > sup) {
            return end();
        } else return iterator(this, domain.lower_bound(max(getInf(), v)));
    }

    //Finds the first available element whose value is lower or equal to v
    iterator upper_bound(Value v) {
        if (assigned()) {
            if (v >= getValue()) return iterator(this, domain.upper_bound(getValue()));
            else return end();
        } else if (v < inf) {
            return end();
        } else return iterator(this, domain.upper_bound(min(getSup(), v)));
    }


    void permuteDomain(int numberOfPermutations);
    void permuteDomain(Value a, Value b);
    ValueCost *sortDomain(vector<Cost> &costs);

    void print(ostream& os);
};

#endif /*TB2ENUMVAR_HPP_*/

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

