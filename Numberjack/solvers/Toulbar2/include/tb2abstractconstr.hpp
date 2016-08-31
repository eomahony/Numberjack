/** \file tb2abstractconstr.hpp
 *  \brief Abstract constraints of predefined arities
 *
 */

#ifndef TB2ABSTRACTCONSTR_HPP_
#define TB2ABSTRACTCONSTR_HPP_

#include "tb2constraint.hpp"
#include "tb2variable.hpp"
#include "tb2enumvar.hpp"
#include "tb2wcsp.hpp"

#include <set>

// Warning! Constraint reconnect() may interfer with variable constraints list iterator used by EnumerateVariable assign method. Do not reconnect a constraint with assigned variables in project functions.

template<class T1>
class AbstractUnaryConstraint : public Constraint
{
protected:
    T1 *x;
    DLink<ConstraintLink> *linkX;

public:
    AbstractUnaryConstraint(WCSP *wcspin, T1 *xx) : Constraint(wcspin), x(xx), linkX(NULL) {
        linkX = xx->link(this,0);
    }

    virtual ~AbstractUnaryConstraint() {delete linkX;}

    bool connected() const FINAL {return !linkX->removed;}
    bool deconnected() const FINAL {return linkX->removed;}
    void deconnect(bool reuse = false) FINAL {
        if (connected()) {
            if (ToulBar2::verbose >= 3) cout << "deconnect " << this << endl;
            x->deconnect(linkX, reuse);
        }
    }
    void reconnect() FINAL {
        if (deconnected()) {
            if (ToulBar2::verbose >= 3) cout << "reconnect " << this << endl;
            assert(linkX->prev == NULL && linkX->next == NULL);
            x->getConstrs()->push_back(linkX, true);
        }
    }

    int arity() const FINAL {return 1;}

    Variable *getVar(int varCtrIndex) const FINAL {return x;}
    Variable *getDACVar(int varCtrIndex) const FINAL {return getVar(varCtrIndex);}

    Variable *getVarDiffFrom( Variable* v ) const {
        if(v != x) return x;
        else exit(EXIT_FAILURE);
    }

    int getIndex(Variable* var) const FINAL
    {
        if(var == x) return 0;
        return -1;
    }

    int getSmallestVarIndexInScope(int forbiddenScopeIndex) FINAL {assert(forbiddenScopeIndex >= 0); assert(forbiddenScopeIndex < 1); return x->wcspIndex;}
    int getSmallestVarIndexInScope() FINAL {return x->wcspIndex;}
    int getSmallestDACIndexInScope(int forbiddenScopeIndex) FINAL {assert(forbiddenScopeIndex >= 0); assert(forbiddenScopeIndex < 1); return x->getDACOrder();}

    void getScope( TSCOPE& scope_inv ) FINAL {
        scope_inv.clear();
        scope_inv[ x->wcspIndex ] = 0;
    }

    ConstraintSet subConstraint() FINAL {
        ConstraintSet subcstr;
        return subcstr;
    }
};


template<class T1, class T2>
class AbstractBinaryConstraint : public Constraint
{
protected:
    T1 *x;
    T2 *y;
    DLink<ConstraintLink> *linkX;
    DLink<ConstraintLink> *linkY;
    int dacvar;

public:
    AbstractBinaryConstraint(WCSP *wcspin, T1 *xx, T2 *yy) : Constraint(wcspin), x(xx), y(yy), linkX(NULL), linkY(NULL) {
        assert(xx != yy);
        linkX = xx->link(this,0);
        linkY = yy->link(this,1);
        setDACScopeIndex();
    }

    AbstractBinaryConstraint(WCSP *wcspin) : Constraint(wcspin,-wcspin->elimBinConstrs.size()-1), x(NULL), y(NULL), linkX(NULL), linkY(NULL)
    { }

    virtual ~AbstractBinaryConstraint() {delete linkX; delete linkY;}

    bool connected() const FINAL {return !linkX->removed && !linkY->removed;}
    bool deconnected() const FINAL {return linkX->removed || linkY->removed;}
    void deconnect(bool reuse = false) FINAL {
        if (connected()) {
            if (ToulBar2::verbose >= 3) cout << "deconnect " << this << endl;
            x->deconnect(linkX, reuse);
            y->deconnect(linkY, reuse);
        }
    }
    void reconnect() FINAL {
        if (deconnected()) {
            if (ToulBar2::verbose >= 3) cout << "reconnect " << this << endl;
            assert(linkX->prev == NULL && linkX->next == NULL);
            x->getConstrs()->push_back(linkX, true);
            assert(linkY->prev == NULL && linkY->next == NULL);
            y->getConstrs()->push_back(linkY, true);
        }
    }

    int arity() const FINAL {return 2;}

    Variable *getVar(int varCtrIndex) const FINAL {return (varCtrIndex == 0)?x:y;}

    Variable *getVarDiffFrom( Variable* v ) const {
        if(v == x) return y;
        else if(v == y) return x;
        else exit(EXIT_FAILURE);
    }

    int getIndex(Variable* var) const FINAL
    {
        if(var == x) return 0;
        else if(var == y) return 1;
        return -1;
    }

    int getSmallestVarIndexInScope(int forbiddenScopeIndex) FINAL {assert(forbiddenScopeIndex >= 0); assert(forbiddenScopeIndex < 2); return (forbiddenScopeIndex)?x->wcspIndex:y->wcspIndex;}
    int getSmallestVarIndexInScope() FINAL {return min(x->wcspIndex, y->wcspIndex);}
    int getDACScopeIndex() const FINAL {return dacvar;}
    void setDACScopeIndex() FINAL {if (x->getDACOrder() < y->getDACOrder()) dacvar = 0; else dacvar = 1;}
    int getSmallestDACIndexInScope(int forbiddenScopeIndex) FINAL {assert(forbiddenScopeIndex >= 0); assert(forbiddenScopeIndex < 2); return (forbiddenScopeIndex)?x->getDACOrder():y->getDACOrder();}
    Variable *getDACVar(int varCtrIndex) const FINAL {return ((varCtrIndex == 0) == (dacvar == 0))?x:y;}

    void getScope( TSCOPE& scope_inv ) FINAL {
        scope_inv.clear();
        scope_inv[ x->wcspIndex ] = 0;
        scope_inv[ y->wcspIndex ] = 1;
    }

    ConstraintSet subConstraint() FINAL {
        ConstraintSet subcstr;
        return subcstr;
    }
};


template<class T1, class T2, class T3>
class AbstractTernaryConstraint : public Constraint
{
protected:
    T1 *x;
    T2 *y;
    T3 *z;
    DLink<ConstraintLink> *linkX;
    DLink<ConstraintLink> *linkY;
    DLink<ConstraintLink> *linkZ;
    int dacvar;

public:
    AbstractTernaryConstraint(WCSP *wcsp, T1 *xx, T2 *yy, T2 *zz) : Constraint(wcsp), x(xx), y(yy), z(zz), linkX(NULL), linkY(NULL), linkZ(NULL) {
        assert(xx != yy);
        assert(xx != zz);
        assert(yy != zz);
        linkX = xx->link(this,0);
        linkY = yy->link(this,1);
        linkZ = zz->link(this,2);
        setDACScopeIndex();
    }

    AbstractTernaryConstraint(WCSP *wcspin) : Constraint(wcspin,-wcspin->elimTernConstrs.size()-MAX_ELIM_BIN-1), x(NULL), y(NULL), z(NULL), linkX(NULL), linkY(NULL), linkZ(NULL)
    { }


    virtual ~AbstractTernaryConstraint() {delete linkX; delete linkY; delete linkZ;}

    bool connected() const FINAL {return !linkX->removed && !linkY->removed && !linkZ->removed;}
    bool deconnected() const FINAL {return linkX->removed || linkY->removed || linkZ->removed;}
    void deconnect(bool reuse = false) FINAL {
        if (connected()) {
            if (ToulBar2::verbose >= 3) cout << "deconnect " << this << endl;
            x->deconnect(linkX, reuse);
            y->deconnect(linkY, reuse);
            z->deconnect(linkZ, reuse);
        }
    }
    void reconnect() FINAL {
        if (deconnected()) {
            if (ToulBar2::verbose >= 3) cout << "reconnect " << this << endl;
            assert(linkX->prev == NULL && linkX->next == NULL);
            //			if (linkX->content.constr->isTriangle()) x->getTriangles()->push_back(linkX, true);
            //			else
            x->getConstrs()->push_back(linkX, true);
            assert(linkY->prev == NULL && linkY->next == NULL);
            //			if (linkY->content.constr->isTriangle()) y->getTriangles()->push_back(linkY, true);
            //			else
            y->getConstrs()->push_back(linkY, true);
            assert(linkZ->prev == NULL && linkZ->next == NULL);
            //			if (linkZ->content.constr->isTriangle()) z->getTriangles()->push_back(linkZ, true);
            //			else
            z->getConstrs()->push_back(linkZ, true);
        }
    }

    int arity() const FINAL {return 3;}

    Variable *getVar(int varCtrIndex) const FINAL
    {
        switch(varCtrIndex) { case 0: return x; break;
        case 1: return y; break;
        case 2: return z; break;
        default: exit(EXIT_FAILURE); }
    }

    Variable *getVarDiffFrom( Variable* v1, Variable* v2 ) const
    {
        if      ((x == v1) && (y == v2)) return z;
        else if ((x == v2) && (y == v1)) return z;
        else if ((x == v1) && (z == v2)) return y;
        else if ((x == v2) && (z == v1)) return y;
        else if ((y == v1) && (z == v2)) return x;
        else if ((y == v2) && (z == v1)) return x;
        else exit(EXIT_FAILURE);
    }

    int getIndex(Variable* var) const FINAL
    {
        if(var == x) return 0;
        else if(var == y) return 1;
        else if(var == z) return 2;
        return -1;
    }


    int getSmallestVarIndexInScope(int forbiddenScopeIndex) FINAL
    {
        assert(forbiddenScopeIndex >= 0); assert(forbiddenScopeIndex < 3);
        switch (forbiddenScopeIndex) {
        case 0: return min(y->wcspIndex,z->wcspIndex); break;
        case 1: return min(x->wcspIndex,z->wcspIndex); break;
        case 2: return min(x->wcspIndex,y->wcspIndex); break;
        default: exit(EXIT_FAILURE);
        }
    }
    int getSmallestVarIndexInScope() FINAL
    {
        int res = min(x->wcspIndex,y->wcspIndex);
        return min(res, z->wcspIndex);
    }
    int getDACScopeIndex() const FINAL {return dacvar;}
    void setDACScopeIndex() FINAL {
        if (x->getDACOrder() < y->getDACOrder() && x->getDACOrder() < z->getDACOrder()) dacvar = 0;
        else if (y->getDACOrder() < x->getDACOrder() && y->getDACOrder() < z->getDACOrder()) dacvar = 1;
        else dacvar = 2;
    }
    int getSmallestDACIndexInScope(int forbiddenScopeIndex) FINAL
    {
        assert(forbiddenScopeIndex >= 0); assert(forbiddenScopeIndex < 3);
        switch (forbiddenScopeIndex) {
        case 0: return min(y->getDACOrder(),z->getDACOrder()); break;
        case 1: return min(x->getDACOrder(),z->getDACOrder()); break;
        case 2: return min(x->getDACOrder(),y->getDACOrder()); break;
        default: exit(EXIT_FAILURE);
        }
    }
    Variable *getDACVar(int varCtrIndex) const FINAL
    {
        if (varCtrIndex == 0) {
            switch (dacvar) {
            case 0: return x;
            case 1: return y;
            case 2: return z;
            default: exit(EXIT_FAILURE);
            }
        } else { // varCtrIndex >= 1
            switch (dacvar) {
            case 0: return ((varCtrIndex == 1) == (y->getDACOrder() < z->getDACOrder()))?y:z;
            case 1: return ((varCtrIndex == 1) == (x->getDACOrder() < z->getDACOrder()))?x:z;
            case 2: return ((varCtrIndex == 1) == (x->getDACOrder() < y->getDACOrder()))?x:y;
            default: exit(EXIT_FAILURE);
            }
        }
    }

    void getScope( TSCOPE& scope_inv ) FINAL {
        scope_inv.clear();
        scope_inv[ x->wcspIndex ] = 0;
        scope_inv[ y->wcspIndex ] = 1;
        scope_inv[ z->wcspIndex ] = 2;
    }

    ConstraintSet subConstraint() FINAL {
        ConstraintSet subcstr;
        set<int> scope;
        for(int k=0; k < arity(); k++) {
            scope.insert(getVar(k)->wcspIndex);
        }
        for(set<int>::iterator itx = scope.begin(); itx != scope.end(); ++itx){
            ConstraintList* xctrs = (wcsp->getVar(*itx))->getConstrs();
            for (ConstraintList::iterator it=xctrs->begin(); it != xctrs->end(); ++it){
                Constraint * ctr = (*it).constr;
                if(ctr->arity() == 2 && scopeIncluded(ctr))subcstr.insert(ctr);
            }
        }

        return subcstr;
    }

};

extern int cmpDAC(const void *var1, const void *var2);

class AbstractNaryConstraint : public Constraint
{
protected:

    int arity_;

    EnumeratedVariable** scope;
    EnumeratedVariable** scope_dac; // scope sorted by increasing DAC order
    TSCOPE scope_inv;

    DLink<ConstraintLink>** links;

    vector<EnumeratedVariable::iterator> it_values; // used by firstlex/nextlex and for separator decomposition tests
    String iterTuple;

public:
    AbstractNaryConstraint(WCSP *wcsp, EnumeratedVariable** scope_in, int arity_in) : Constraint(wcsp), arity_(arity_in)
    {
        scope = new EnumeratedVariable* [arity_];
        scope_dac = new EnumeratedVariable* [arity_];
        links = new DLink<ConstraintLink>* [arity_];

        for(int i=0; i < arity_; i++) {
            EnumeratedVariable* var = scope_in[i];
            scope_inv[ var->wcspIndex ] = i;
            scope[i] = var;
            scope_dac[i] = var;
            links[i] = var->link(this,i);
            iterTuple.append(1, CHAR_FIRST);
        }
        setDACScopeIndex();
    }

    AbstractNaryConstraint(WCSP *wcsp) : Constraint(wcsp)
    {
    }

    virtual ~AbstractNaryConstraint() {}

    int arity() const FINAL {return arity_;}
    Long getDomainInitSizeProduct(); // warning! return LONGLONG_MAX if overflow occurs

    Variable *getVar(int varCtrIndex) const FINAL {
        assert(varCtrIndex < arity_);
        return scope[varCtrIndex];
    }

    int getIndex(Variable* var) const FINAL {
        int index = var->wcspIndex;
        map<int,int>::const_iterator it = scope_inv.find(index);
        if(it == scope_inv.end()) return -1;
        else return it->second;
    }

    bool connected(int varIndex) const {return !links[varIndex]->removed;}
    bool deconnected(int varIndex) const {return links[varIndex]->removed;}

    bool connected() const FINAL {
        for(int i=0;i<arity_;i++) if(!links[i]->removed) return true;
        return false;
    }

    bool deconnected() const FINAL {
        for(int i=0;i<arity_;i++) if(!links[i]->removed) return false;
        return true;
    }

    void deconnect(int varIndex, bool reuse = false) {
        scope[varIndex]->deconnect( links[varIndex], reuse );
    }

    void deconnect(bool reuse = false) FINAL {
        if (connected()) {
            if (ToulBar2::verbose >= 3) cout << "deconnect " << this << endl;
            for(int i=0;i<arity_;i++) deconnect(i, reuse);
        }
    }

    virtual void reconnect() {
        if (deconnected()) {
            if (ToulBar2::verbose >= 3) cout << "reconnect " << this << endl;
            for(int i=0;i<arity_;i++) {
                assert(links[i]->prev == NULL && links[i]->next == NULL);
                scope[i]->getConstrs()->push_back(links[i], true);
            }
        }
    }

    virtual Cost eval( const String& t ) {cout << "dummy eval on: " << *this << endl; return MIN_COST;}
//    virtual void insertTuple( String t, Cost c, EnumeratedVariable** scope_in ) { }

    int getSmallestVarIndexInScope(int forbiddenScopeIndex) FINAL
    {
        assert(forbiddenScopeIndex >= 0); assert(forbiddenScopeIndex < arity_);
        int indexmin = INT_MAX;
        for(int i=0; i < arity_; i++) if (i != forbiddenScopeIndex) {
            if (scope[i]->wcspIndex < indexmin) {
                indexmin = scope[i]->wcspIndex;
            }
        }
        return indexmin;
    }
    int getSmallestVarIndexInScope() FINAL
    {
        int indexmin = INT_MAX;
        for(int i=0; i < arity_; i++) {
            if (scope[i]->wcspIndex < indexmin) {
                indexmin = scope[i]->wcspIndex;
            }
        }
        return indexmin;
    }
    void getScope( TSCOPE& scope_inv_in ) FINAL {
        scope_inv_in = scope_inv;
    }

    /// \warning scope_dac uses pointers not indexes!
    /// \warning additional side-effect: updates scope_inv with current variable wcspIndexes
    void setDACScopeIndex() FINAL {
        scope_inv.clear();
        for(int i=0; i < arity_; i++) {
            scope_inv[ scope[i]->wcspIndex ] = i;
        }
        qsort(scope_dac, arity_, sizeof(EnumeratedVariable *), cmpDAC);
    }
    int getSmallestDACIndexInScope(int forbiddenScopeIndex) FINAL
    {
        assert(forbiddenScopeIndex >= 0); assert(forbiddenScopeIndex < arity_);
        int indexmin = INT_MAX;
        for(int i=0; i < arity_; i++) if (i != forbiddenScopeIndex) {
            if (scope[i]->getDACOrder() < indexmin) {
                indexmin = scope[i]->getDACOrder();
            }
        }
        return indexmin;
    }
    Variable *getDACVar(int varCtrIndex) const FINAL {
        assert(varCtrIndex < arity_);
        return scope_dac[varCtrIndex];
    }

    ConstraintSet subConstraint() FINAL {
        ConstraintSet subcstr;
        set<int> scope;
        for(int k=0; k < arity_; k++) {
            scope.insert(getVar(k)->wcspIndex);
        }
        for(set<int>::iterator itx = scope.begin(); itx != scope.end(); ++itx){
            ConstraintList* xctrs = (wcsp->getVar(*itx))->getConstrs();
            for (ConstraintList::iterator it=xctrs->begin(); it != xctrs->end(); ++it){
                Constraint * ctr = (*it).constr;
                if(ctr->arity() < arity() && scopeIncluded(ctr)) subcstr.insert(ctr);
            }
        }

        return subcstr;
    }

    void firstlex();
    bool nextlex( String& t, Cost& c);

    void projectNaryBeforeSearch();
};

#define AbstractGlobalConstraint AbstractNaryConstraint


#endif /*TB2ABSTRACTCONSTR_HPP_*/

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

