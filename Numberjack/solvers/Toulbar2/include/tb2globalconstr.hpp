/** \file tb2globalconstr.hpp
 *  \brief Global Constraint using enumerated variables with parameters read from file
 * 
 */

#ifndef TB2GLOBALCONSTR_HPP_
#define TB2GLOBALCONSTR_HPP_

#include "tb2abstractconstr.hpp"
#include "tb2naryconstr.hpp"
#include "tb2enumvar.hpp"
#include "tb2wcsp.hpp"

#include <map>
#include <set>
using namespace std;

class GlobalConstraint : public AbstractGlobalConstraint
{

protected:	

    vector<StoreCost> *deltaCost; // the cost transferred from/to nary-constraint, must be backtractable
    vector<StoreCost> *extendedCost; // the cost extended nary-constraint, must be backtractable
    StoreCost projectedCost; // the cost projected to the C_null constraint, must be backtractable
    StoreInt nonassigned;       // nonassigned variables during search, must be backtrackable (storeint) !

    set<int> *fullySupportedSet;

    map<Value, Cost> EACCost;
    vector<vector<Cost> > preUnaryCosts;
    int currentVar;
    bool needPropagateAC, needPropagateDAC, needPropagateEAC;

    // mode : the cost measure
    // def : the cost of the violation edge
    int currentDepth;
    Cost def;
    int mode;
    map<string, int> modeEnum;

    int count_nic, count_gac, count_fdac, count_edac, error;

    // find the minimum cost of the tuple when varindex = v for each v in
    // D(varindex)
    virtual void findProjection(int varindex, map<Value, Cost> &delta) {}
    // check and remove from the constraint structure the value already removed
    // by others
    virtual void checkRemoved(vector<int> &rmv) {}
    // extend the cost stored in deltas[i] from the unary constraint of
    // supports[i] to the constraint struture
    virtual void changeAfterExtend(vector<int> &supports, vector<map<Value, Cost> > &deltas){}
    virtual void changeAfterExtend(int support, map<Value, Cost> &delta){
        vector<int> supports; supports.push_back(support);
        vector<map<Value, Cost> > deltas; deltas.push_back(delta);
        changeAfterExtend(supports, deltas);
    }
    // project the cost stored in deltas[i] to the unary constraint of
    // supports[i] from the constraint struture
    virtual void changeAfterProject(vector<int> &supports, vector<map<Value, Cost> > &deltas){}
    virtual void changeAfterProject(int support, map<Value, Cost> &delta){
        vector<int> supports; supports.push_back(support);
        vector<map<Value, Cost> > deltas; deltas.push_back(delta);
        changeAfterProject(supports, deltas);
    }
    void project(int index, Value value, Cost cost, bool delayed = false);
    void extend(int index, Value value, Cost cost);
    // undo the previous extension
    virtual void undoExtend() {}
    // compute the original cost of the tuple s (i.e. cost without projection)
    virtual Cost evalOriginal( String s ) {return 0;}        	

    // compute the minimum cost of the tuple from all feasible tuples



public:
    // construtor
    GlobalConstraint(WCSP *wcsp, EnumeratedVariable** scope_in, int arity_in, Cost defval);
    // destrutor
    virtual ~GlobalConstraint();

    virtual void setBaseCost(Cost cost)  {def = cost;}
    virtual void setSemantics(const string &semantic) {mode = modeEnum[semantic];}

    bool isGlobal() const {return true;}

    // evaluate the cost of the tuple
    virtual Cost eval( String s );

    double computeTightness() { return 0; }
    virtual string getName() {return "global constraint";}
    virtual void print(ostream& os);

    // initialize the constraint structure for enforing consistency
    void init();
    virtual void initStructure() {}
    // clear up the structure
    virtual void end() {}

    // used for enforcing "EDGAC", still have some bugs
    virtual bool isEAC(int index, Value a);
    virtual void fillEAC2(int index);
    //virtual void getCostsWithUnary(int index, map<Value, Cost> &costs);
    virtual void propagateEAC();
    virtual void findFullSupportEAC(int index);
    virtual void linkCostProvidingPartition(int index, Variable* support) {
        int sindex = -1;
        for (int i=0;i<arity_ && sindex == -1;i++) {
            if (getVar(i) == support) sindex = i;
        }
        if ((sindex != index) && (index != -1)) {
            fullySupportedSet[index].insert(sindex);
        }
    }
    virtual void showCostProvidingPartition(int i) {
        cout << getVar(i)->getName() << ": ";
        for (set<int>::iterator j = fullySupportedSet[i].begin(); j !=
                fullySupportedSet[i].end();j++) {
            if (getVar(*j)->unassigned()) {
                cout << getVar(*j)->getName() << " ";
            }
        } cout << endl;
        EnumeratedVariable* x = (EnumeratedVariable*)getVar(i);
        for (EnumeratedVariable::iterator v = x->begin();v != x->end();++v) {
            cout << EACCost[*v] << " ";
        } cout << endl;
    }
    virtual void showCostProvidingPartition() {
        for (int i=0;i<arity_;i++) {
            showCostProvidingPartition(i);
        }
    }

    //Trick! instead of doing remove(index) now, let AC queue do the job.
    //So several incdec events on the same constraint can be merged into one AC event
    virtual void increase(int index) {((EnumeratedVariable*)getVar(index))->queueAC();}
    virtual void decrease(int index) {((EnumeratedVariable*)getVar(index))->queueAC();}

    // read the parameter of the constraint parameter from the file
    virtual void read(istream &file) {}

    // return the minimum cost of the tuples
    virtual Cost getMinCost() {return 0;}

    virtual bool universal() {return false;}

    /*virtual void valueRemoved(int index, Value value) {
	  if (ToulBar2::consistencyLevel == FINE_IC) {
	  propagateStrongNIC();
	  } 
	  }*/
    // Still consider whether we should reduce to binary, as done in nary
    // constraints
    virtual void assign(int varIndex);

    // function used for propagation
    virtual void remove(int index);
    virtual void projectFromZero(int index);
    void pushAll() {
        for (int i=0;i<arity_;i++) {
            EnumeratedVariable* x = (EnumeratedVariable*)getVar(i);
            if (x->unassigned()) {
                x->queueEAC1();
                x->queueDAC();
                x->queueAC();
            }
        }
    }

    // propagation for the whole constraint
    virtual void propagate();
    virtual void propagateDAC();
    virtual void propagateAC();
    virtual void propagateStrongNIC();
    virtual void propagateNIC();

    // used for FDGAC*
    virtual void findFullSupport(int index) {
        vector<int> provide;
        for (int i=index+1;i < arity_;i++)
            if (getVar(i)->unassigned()) provide.push_back(i);
        findFullSupport(index, provide, false);
    }
    // used for FDGAC* and EDGAC*
    // find the full support w.r.t. the set of variable support
    virtual void findFullSupport(int index, vector<int> &support, bool isEAC);
    // used for GAC*
    virtual void findSupport(int varindex);
    // used for Strong NIC*
    virtual void checkMinCost(int varindex);

    bool verify() {return true;}
};



#endif /*TB2GLOBALCONSTR_HPP_*/

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

