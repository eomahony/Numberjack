/** \file tb2dpglobalconstr.hpp
 *  \brief base class of polynomially decomposable cost functions
 */

#ifndef TB2DPGLOBALCONSTR_HPP_
#define TB2DPGLOBALCONSTR_HPP_

#include <vector>
#include <algorithm>
#include "tb2globalconstr.hpp"
using namespace std;

class DPGlobalConstraint : public GlobalConstraint {
private:
    vector<bool> * zero;
    vector<Cost> * preUnaryCosts;

    bool initialized;

    void clear();
    void record(Value *tuple);
    void findSupport(int var, bool &changed);

protected:
    DPGlobalConstraint(WCSP * wcsp, EnumeratedVariable ** scope, int arity);
    virtual ~DPGlobalConstraint();

    virtual void initMemoization() {}

    virtual void initStructure() {
        if (!initialized) {
            initMemoization();
            initialized = true;
        }
    }

    typedef pair<Cost, Value*> Result;
    virtual Cost minCostOriginal() = 0;
    virtual Cost minCostOriginal(int var, Value val, bool changed) = 0;
    virtual Result minCost(int var, Value val, bool changed) = 0;

    virtual void propagateNIC();
    virtual void propagateStrongNIC();
    virtual void propagateAC();
    virtual void propagateDAC();

    //EAC
    virtual bool isEAC(int var, Value val);
    virtual void findFullSupportEAC(int var);

};

#endif //TB2GLOBALCONSTR3_HPP_

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

