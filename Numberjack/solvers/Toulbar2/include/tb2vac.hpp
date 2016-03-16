/** \file tb2vac.hpp
 *  \brief Enforce VAC in a WCSP.
 */

#ifndef TB2VAC_HPP_
#define TB2VAC_HPP_

#include <stack>
#include <set>
#include <list>
#include "tb2btqueue.hpp"
#include "tb2vacutils.hpp"


class tVACStat;

typedef map<Cost,int> tScale;


/**
 * The class that enforces VAC
 */
class VACExtension {

private:

    WCSP *wcsp;               /**< Reference to the WCSP that will be processed */
    Queue VAC;                /**< Non backtrackable list; AC2001 queue used for Pass1 inside VAC */
    Queue SeekSupport;        /**< Non backtrackable list; collect all variables with a deletion caused by binary constraints during Pass1 */
    BTQueue VAC2;             /**< Backtrackable list; updated during AC and EDAC */
    Long nbIterations;        /**< Incremented at each pass, used as a TimeStamp */
    int inconsistentVariable; /**< WipeOut variable, Used also to check after enforcePass1() if the network is VAC */

    Cost itThreshold;         /**< The cost threshold (theta) for the iterative threshold descent */
    int breakCycles;          /**< Number of iterations with no c0 increase */
    tScale scaleCost;         /**w The list of all costs used in the WCSP ? */
    list<Cost> scaleVAC;      /**< The scale of costs used for the thresold descent */

    Cost minlambda;           /**< The amount of cost that will go to c0 (lambda) */

    stack< pair<int, Value> > *queueP;     /**< Values removed by hard AC (created in Pass1, used in Pass2) */
    stack< pair<int, Value> > *queueR;     /**< Minimal set of deletions needed to increase c0 (created in Pass2, used in Pass3) */

    void enforcePass1 ();     /**< Enforces instrumented hard AC (Phase 1) */
    bool enforcePass1( VACVariable *xj, VACBinaryConstraint* cij);  /**< Revises /a xj wrt /a cij and updates /a k */
    bool checkPass1 () const; /**< Checks if Bool(P) is AC */
    void enforcePass2 ();     /**< Finds a minimal set of deletions needed for wipeout and computes k and lambda */
    bool enforcePass3 ();     /**< Project and extends costs to increase c0 according to the plan */
    void enforcePass3VACDecomposition ();    /**< Enforces VAC decomposition pass 3 (substract cost and decrease top) */

    void reset();             /**< Cleanup for next iteration: clean Q, selects variables for AC2001 queue */

    map<int,tVACStat*> heapAccess;
    vector<tVACStat*>   heap;
    Cost sumlb;
    Long nlb;
    Long sumvars;
    int sumk;
    int theMaxK;

    int bneckVar;
    VACBinaryConstraint *bneckCF;
    Cost bneckCost;

    EnumeratedVariable* nearIncVar; /**< A variable whose domain is reduced to 1 value during Pass1 -- Unused apparently */
    Cost atThreshold;               /**< The value of itThreshold when this single value variable appeared -- Unused apparently */

public:

    VACExtension (WCSP *w);
    ~VACExtension ();

    bool firstTime() { return nbIterations == 0; }  /**< Is it the first iteration ? */

    bool isVAC () const;  /**< Is the WCSP VAC-epsilon ? Pass1 must be enforced */
    bool propagate ();    /**< Starts one VAC iteration */
    bool isNull (Cost c) const { return (c < itThreshold); }  /**< is the Cost significant (above itThreshold aka theta) */

    void clear();         /**< empty VAC queue */
    void queueVAC(DLink<VariableWithTimeStamp> *link);
    void queueSeekSupport(DLink<VariableWithTimeStamp> *link);
    void queueVAC2(DLink<Variable *> *link);
    void dequeueVAC2(DLink<Variable *> *link);

    void init();
    void iniThreshold();  /**< Initialize itThreshold to the strongest cost in the cost scale */
    Cost getThreshold() { return itThreshold; }
    void nextScaleCost(); /**< Sets ItThreshold to the next scale */
    void histogram( Cost c );
    void histogram();     /**< Computes the ScaleVAC splitting the cost scale in 20 buckets or less */

    set<int> singletonI;
    set<int> singleton;

    // int varAssign;
    // Value valAssign;
    // void assign(int varIndex, Value newValue) { varAssign = varIndex; valAssign = newValue; }

    void afterPreprocessing();
    void iniSingleton();
    void updateSingleton();
    void removeSingleton();

    int  getHeuristic();
    // int  getVarACDom( int i );
    // Cost getVarCostStat( int i );
    // Long getVarTimesStat( int i );
    // void updateStat(Cost lambda);
    void printStat(bool ini = false);
    void printTightMatrix();

    void minsumDiffusion(); /**< MinSumDiffusion implementation */
};


#endif /*TB2VAC_HPP_*/

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

