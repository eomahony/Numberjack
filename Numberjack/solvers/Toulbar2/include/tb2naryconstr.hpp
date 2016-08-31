#ifndef TB2NARYCONSTR_HPP_
#define TB2NARYCONSTR_HPP_

#include "tb2abstractconstr.hpp"
#include "tb2ternaryconstr.hpp"
#include "tb2enumvar.hpp"
#include "tb2wcsp.hpp"

#include <set>


class NaryConstraint : public AbstractNaryConstraint
{
    typedef map<String,Cost> TUPLES;
    TUPLES* pf;
    Cost *costs;
    ptrdiff_t costSize;
    Cost default_cost;          // default cost returned when tuple t is not found in TUPLES (used by function eval(t))
    StoreInt nonassigned;       // nonassigned variables during search, must be backtrackable (StoreInt)!
    ConstraintSet* filters;
    TUPLES::iterator  tuple_it;
    String evalTuple;
    vector<Long> conflictWeights; // used by weighted degree heuristics

public:
    NaryConstraint(WCSP *wcsp, EnumeratedVariable** scope_in, int arity_in, Cost defval, Long nbtuples = 0);
    NaryConstraint(WCSP *wcsp);
    virtual ~NaryConstraint();

    bool extension() const FINAL {return true;}

    void reconnect() {
        if (deconnected()) {
            nonassigned = arity_;
            AbstractNaryConstraint::reconnect();
        }
    }
    int getNonAssigned() const {return nonassigned;}

    Long getConflictWeight() const {return Constraint::getConflictWeight();} 
    Long getConflictWeight(int varIndex) const {assert(varIndex>=0);assert(varIndex<arity_);return conflictWeights[varIndex]+Constraint::getConflictWeight();}
    void incConflictWeight(Constraint *from) {
        //assert(fromElim1==NULL);
        //assert(fromElim2==NULL);
        if (from==this) {
            Constraint::incConflictWeight(1);
        } else if (deconnected()) {
            for (int i=0; i<from->arity(); i++) {
                int index = getIndex(from->getVar(i));
                if (index>=0) { // the last conflict constraint may be derived from two binary constraints (boosting search), each one derived from an n-ary constraint with a scope which does not include parameter constraint from
                    assert(index < arity_);
                    conflictWeights[index]++;
                }
            }
        }
    }

    ptrdiff_t getCostsIndex( const String& s) const {
        ptrdiff_t index = 0;
        ptrdiff_t base = 1;
        for (int i=arity_-1; i>=0; --i) {
            index += (s[i] - CHAR_FIRST)*base;
            base *= ((EnumeratedVariable *) getVar(i))->getDomainInitSize();
        }
        assert(base == costSize);
        assert(index < costSize);
        assert(index >= 0);
        return index;
    }
    Long size() const FINAL {return (Long) (pf)?pf->size():((costs)?costSize:0);}
    Long space() const FINAL {return ((pf)?((Long) pf->size() * (sizeof(Cost) + arity_*sizeof(Char))):((costs)?((Long) costSize * sizeof(Cost)):0));} // actual memory space (not taking into account map space overhead)
    Long space(Long nbtuples) const {return (nbtuples < LONGLONG_MAX / ((Long) (sizeof(Cost) + arity_*sizeof(Char))))?(nbtuples * (sizeof(Cost) + arity_*sizeof(Char))):LONGLONG_MAX;} // putative memory space
    bool expandtodo() {return space() > getDomainInitSizeProduct();} // should be getDomainInitSizeProduct() * sizeof(Cost) ?
    bool expandtodo(Long nbtuples) {return space(nbtuples) > getDomainInitSizeProduct();} // getDomainInitSizeProduct() * sizeof(Cost) ?
    void expand();

    bool consistent( const String& t );
    Cost eval( const String& s );
    Cost eval( const String& s, EnumeratedVariable** scope_in );
    Cost evalsubstr( const String& s, Constraint* ctr ) FINAL {return evalsubstrAny(s,ctr);}
    Cost evalsubstr( const String& s, NaryConstraint* ctr ) FINAL {return evalsubstrAny(s,ctr);}
    template<class T> Cost evalsubstrAny( const String& s, T* ctr ) {
        int count = 0;

        for(int i=0;i<arity_;i++) {
            int ind = ctr->getIndex( getVar(i) );
            if(ind >= 0) { evalTuple[i] = s[ind]; count++; }
        }
        assert(count <= arity_);

        Cost cost;
        if(count == arity_) cost = eval( evalTuple );
        else cost = MIN_COST;

        return cost;
    }


    Cost getDefCost() { return default_cost; }
    void keepAllowedTuples( Cost df );

    void resetFilters();
    void fillFilters();

    void project( EnumeratedVariable* x );
//    void sum( NaryConstraint* nary );
    double computeTightness();

    void first();
    bool next( String& t, Cost& c);

    void first(EnumeratedVariable* a, EnumeratedVariable* b);
    bool separability(EnumeratedVariable* a, EnumeratedVariable* b);
    void separate(EnumeratedVariable *a, EnumeratedVariable *c);

    void setTuple( const String& tin, Cost c ) FINAL {
        if (pf) (*pf)[tin] = c;
        else costs[getCostsIndex(tin)] = c;
    }
    void addtoTuple( const String& tin, Cost c ) FINAL {
        if (pf) (*pf)[tin] += c;
        else costs[getCostsIndex(tin)] += c;
    }
//    void setTuple( const String& tin, Cost c, EnumeratedVariable** scope_in );
//    void addtoTuple( const String& tin, Cost c, EnumeratedVariable** scope_in );

    void addtoTuples( Cost c ); // c can be positive or negative (if greater than the minimum cost)
    void addtoTuples( EnumeratedVariable* x, Value v, Cost c); // the same operation but restricted to tuples with x assigned to v

    void setInfiniteCost(Cost ub);
    void insertSum( const String& t1, Cost c1, Constraint* ctr1, const String& t2, Cost c2, Constraint* ctr2, bool bFilters = false );
//    void permute( EnumeratedVariable** scope_in );

    void projectxy( EnumeratedVariable* x, EnumeratedVariable* y, TUPLES& fproj);
//    void projectxyz( EnumeratedVariable* x, EnumeratedVariable* y, EnumeratedVariable* z, TUPLES& fproj);
//    void preproject3();
    void preprojectall2();

    void assign(int varIndex);

    void projectNary();
    void projectNaryTernary(TernaryConstraint* xyz);
    void projectNaryBinary(BinaryConstraint* xy);

    void propagate() {
        for(int i=0;connected() && i<arity_;i++) {
            if (getVar(i)->assigned()) assign(i);
        }
    };

    bool verify() {return true;}
    void increase(int index) {}
    void decrease(int index) {}
    void remove(int index) {}

//    void starrule(const String& t, Cost minc);
    void projectFromZero(int index) {}

    void fillRandom();
    void print(ostream& os);
    void dump(ostream& os, bool original = true);
};
#endif /*TB2NARYCONSTR_HPP_*/

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

