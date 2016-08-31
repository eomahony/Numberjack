/** \file tb2variable.hpp
 *  \brief Abstract Variable class extended with unary costs.
 * 
 */

#ifndef TB2VARIABLE_HPP_
#define TB2VARIABLE_HPP_

#include "tb2btlist.hpp"
#include "tb2queue.hpp"
#include "tb2domain.hpp"

#include <set>
/*
 * Main class
 * 
 */

class Variable : public WCSPLink
{
protected:
    string name;
    int dac;

    Long timestamp;
    int pos; // current position in the list of unassigned variables

    StoreValue inf;
    StoreValue sup;

    ConstraintList constrs;
    //ConstraintList triangles;

    // incremental NC data
    StoreCost maxCost;
    StoreValue maxCostValue;
    StoreInt NCBucket;
    DLink< Variable * > linkNCBucket;

    DLink<VariableWithTimeStamp> linkNCQueue;
    DLink<VariableWithTimeStamp> linkIncDecQueue;
    DLink<VariableWithTimeStamp> linkEliminateQueue;

    void setMaxUnaryCost(Value a, Cost cost);
    void changeNCBucket(int newBucket);
    void conflict();

    // make it private because we don't want copy nor assignment
    Variable(const Variable &x);
    Variable& operator=(const Variable &x);


public:    
    Variable(WCSP *w, string n, Value iinf, Value isup);

    virtual ~Variable() {}

    virtual bool enumerated() const =0;

    string getName() const {return name;}
    int getDACOrder() const {return dac;}
    void setDACOrder(int order) {dac = order;}
    Value getInf() const {return inf;}
    Value getSup() const {return sup;}
    Value getValue() const {assert(assigned()); return inf;}
    virtual unsigned int getDomainSize() const =0;
    int getCurrentVarId();
    void setCurrentVarId(int idx);

    bool assigned() const {return inf == sup;}
    bool unassigned() const {return inf != sup;}
    virtual bool canbe(Value v) const {return v >= inf && v <= sup;}
    virtual bool cannotbe(Value v) const {return v < inf || v > sup;}

    virtual void increase(Value newInf, bool isDecision = false) =0;
    virtual void decrease(Value newSup, bool isDecision = false) =0;
    virtual void remove(Value remValue, bool isDecision = false) =0;
    virtual void assign(Value newValue, bool isDecision = false) =0;
    virtual void assignLS(Value newValue, ConstraintSet& delayedCtrs) =0;

    //    ConstraintList *getTriangles() {return &triangles;}
    ConstraintList *getConstrs() {return &constrs;}
    int getDegree() {return constrs.getSize();}
    int getTrueDegree();
    Double getMaxElimSize(); /// \brief returns estimated size of the resulting cost function (including this variable) to eliminate itself
    Long getWeightedDegree();
    void resetWeightedDegree();
    DLink<ConstraintLink> *link(Constraint *c, int index);
    void sortConstraints();
    virtual void eliminate() {cout << "variable elimination not implemented!" << endl;};

    BinaryConstraint* getConstr( Variable* x ); 
    TernaryConstraint* getConstr( Variable* x, Variable* y ); 
    TernaryConstraint* existTernary(); 
    double strongLinkedby( Variable* &strvar,  TernaryConstraint* &tctr1, TernaryConstraint* &tctr2  );
    void deconnect(DLink<ConstraintLink> *link, bool reuse = false);

    void projectLB(Cost cost);

    virtual Cost getInfCost() const =0;
    virtual Cost getSupCost() const =0;
    virtual void projectInfCost(Cost cost) =0;
    virtual void projectSupCost(Cost cost) =0;
    virtual Cost getCost(const Value value) const =0;

    virtual Value getSupport() const {return inf;}      // If there is no defined support then return inf

    Cost getMaxCost() const {return maxCost;}
    Value getMaxCostValue() const {return maxCostValue;}

    virtual void propagateNC() =0;    
    virtual bool verifyNC() =0;
    virtual bool isEAC() {return true;}
    virtual bool verifyDEE() {return true;}

    void queueNC();
    void queueInc();
    void queueDec();
    void queueEliminate();
    virtual void queueDEE() {}

    void propagateIncDec(int incdec);

    /**********************************************************************/
    //   added for tree decomposition stuff	
    StoreInt cluster;
    void setCluster( int c ) { cluster = c; }
    int  getCluster()        { return cluster; }

    BinaryConstraint* getConstr( Variable* x, int cid );
    TernaryConstraint* getConstr( Variable* x, Variable* y, int cid );


    bool isSep_;
    void setSep() 		   {isSep_ = true;}
    bool isSep() 		   {return isSep_;}

    typedef set< pair<int,int> >   TSepLink;   // set of pairs <cluster in wihch the variable appears,
    //  			    position of the variable in the delta structure>
    TSepLink clusters;

    void addCluster( int c, int pos ) {
        clusters.insert( pair<int,int>(c,pos) );
    }

    TSepLink::iterator itclusters;

    int nbSeparators() { return clusters.size(); }

    void beginCluster() { itclusters = clusters.begin(); }

    bool nextCluster(int& c, int& pos) {
        if(itclusters != clusters.end()) {
            c = (*itclusters).first;
            pos = (*itclusters).second;
            ++itclusters;
            return true;
        }
        else return false;
    }

    /*
    class iterator;
    friend class iterator;
    class iterator {
    public:
        virtual Value operator*() const =0;

        virtual iterator &operator++() =0;     // Prefix form
        virtual iterator &operator--() =0;     // Prefix form

        // To see if you're at the end:
        virtual bool operator==(const iterator &iter) const =0;
        virtual bool operator!=(const iterator &iter) const =0;
    };
    virtual iterator begin() =0;
    virtual iterator end() =0;
    virtual iterator rbegin() =0;
    virtual iterator rend() =0;

    //Finds the first available element whose value is greater or equal to v
    virtual iterator lower_bound(Value v) =0;

    //Finds the first available element whose value is lower or equal to v
    virtual iterator upper_bound(Value v) =0;
     */

    virtual void print(ostream& os) =0;

    friend ostream& operator<<(ostream& os, Variable &var);
};

#endif /*TB2VARIABLE_HPP_*/

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

