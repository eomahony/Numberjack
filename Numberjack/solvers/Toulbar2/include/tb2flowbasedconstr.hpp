/** \file tb2flowbasedconstr.hpp
 *  \brief Global constraint using flow networks structure for propagation
 *
 */

#ifndef TB2FLOWBASEDCONSTR_HPP_
#define TB2FLOWBASEDCONSTR_HPP_

#include "tb2globalconstr.hpp"
#include "tb2graph.hpp"

class FlowBasedGlobalConstraint : public GlobalConstraint 
{
protected:

    // graph : the flow network corresponding to current cost measure
    Graph* graph;
    StoreCost cost;

    // zeroEdges : store the edges containing in a zero cycle.
    // zeroEdges[i][j] returns true if (i,j) lies in a zero cycle.
    bool **zeroEdges;

    bool hasConfigOrganized;

    // mapval : map the value to the node in the network
    map<Value, int> mapval;

    // mapto : map the assignment to the corresponding edge in the network
    virtual pair<int,int> mapto(int varindex, Value val) {
        return make_pair(varindex+1, mapval[val]);
    }

    // compute the projection from the network. store the projected cost in
    // the map delta
    virtual void findProjection(Graph &graph, StoreCost &cost, int varindex, map<Value, Cost> &delta);
    void findProjection(int varindex, map<Value, Cost> &delta) {
        findProjection(*graph, cost, varindex, delta);
    }

    // check whether the network corresponding to the current domains
    // remove any edge which is corresponded to an infeasible assignment
    virtual void checkRemoved(Graph &graph, StoreCost &cost, vector<int> &rmv);
    void checkRemoved(vector<int> &rmv) {
        checkRemoved(*graph, cost, rmv);
    }

    virtual void changeAfterExtend(vector<int> &supports, vector<map<Value, Cost> > &deltas);
    virtual void changeAfterProject(vector<int> &supports, vector<map<Value, Cost> > &deltas);
    virtual void undoExtend() {}

    // construct the flow network
    virtual size_t GetGraphAllocatedSize() = 0;
    virtual void buildGraph(Graph &g) = 0;
    //inline void buildGraph() {buildGraph(graph);}
    inline void augmentGraph(int varindex, map<Value, Cost> &delta) {
        augmentStructure(*graph, cost, varindex, delta);
        graph->removeNegativeCycles(cost);
    }

    // construct the flow in the network
    // the network changed to the residual network
    virtual Cost constructFlow(Graph &g);
    //inline Cost constructFlow() {return constructFlow(graph);}

    // compute the domains of a variable from the network
    virtual void getDomainFromGraph(Graph &graph, int varindex, vector<int> &domain);

    // augment the network by increasing the weight of the edegs
    // corresponding to the assignment of varindex according to delta
    virtual void augmentStructure(Graph &graph, StoreCost &cost, int varindex, map<Value, Cost> &delta);

    // compute the cost according to the original cost structure
    virtual Cost evalOriginal(const String& s) {return MIN_COST;}
    virtual Cost getMinCost () {
        return cost;
    }

public:
    FlowBasedGlobalConstraint(WCSP *wcsp, EnumeratedVariable** scope_in, int
            arity_in);

    ~FlowBasedGlobalConstraint() { }

    virtual void read(istream &file) {}
    virtual void organizeConfig() {}
    virtual void initStructure();
    virtual void end() {}

    //void propagate();

    // check whether the consistency is achieved
    bool verify(){return true;}

#ifndef NDEBUG
    bool isStrongNIC();
    bool isGAC();
    bool isFDGAC();
    bool isEDGAC();
#endif

};


#endif /*TB2FLOWBASEDCONSTR_HPP_*/

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

