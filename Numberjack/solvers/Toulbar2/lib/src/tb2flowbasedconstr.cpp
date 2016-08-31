#include "tb2flowbasedconstr.hpp"
#include "tb2wcsp.hpp"

FlowBasedGlobalConstraint::FlowBasedGlobalConstraint(WCSP *wcsp, EnumeratedVariable** scope_in, int arity_in) : GlobalConstraint(wcsp, scope_in, arity_in, 0),
        graph(NULL), cost(MIN_COST), zeroEdges(NULL), hasConfigOrganized(false)
{
}

Cost FlowBasedGlobalConstraint::constructFlow(Graph &g) {
    pair<int, Cost> result = g.minCostFlow(MIN_COST, g.size()-1);
    return result.second;
}

void FlowBasedGlobalConstraint::initStructure() {

    if (!hasConfigOrganized) {                
        organizeConfig();
        hasConfigOrganized = true;
    }

    if (graph == NULL) {
        size_t graphSize = GetGraphAllocatedSize();
        graph = new Graph(graphSize, arity_);

        if (zeroEdges == NULL) {
            zeroEdges = new bool*[graph->size()];
            for (int i=0;i<graph->size();i++) zeroEdges[i] = new bool[graph->size()];
        }
        for (int i=0;i<graph->size();i++) {
            for (int j=0;j<graph->size();j++) {
                zeroEdges[i][j] = false;
            }
        }

        buildGraph(*graph);
        cost = constructFlow(*graph);
    }


}

void FlowBasedGlobalConstraint::checkRemoved(Graph &graph, StoreCost &cost, vector<int> &rmv) {

    //if (ToulBar2::GCLevel == LC_NC) return;

    pair<Cost, bool> result;
    vector<int> cDomain, cDomain2;
    bool deleted = false;
    //for (int i=0;i<arity_;i++) {
    for (vector<int>::iterator it = rmv.begin();it != rmv.end();it++)
    {
        int i = *it;
        cDomain.clear();
        getDomainFromGraph(graph, i, cDomain);
        sort(cDomain.begin(), cDomain.end());
        EnumeratedVariable* y = (EnumeratedVariable*)getVar(i);
        for (EnumeratedVariable::iterator v = y->begin(); v != y->end();++v) {
            vector<int>::iterator it = find(cDomain.begin(), cDomain.end(), *v);
            if (it == cDomain.end()) {
                cout << "non exist a value ?" << endl;
                for (vector<int>::iterator v=cDomain.begin();v != cDomain.end();v++) {
                    cout << *v << " ";
                } cout << endl;
                for (EnumeratedVariable::iterator v = y->begin(); v != y->end();++v) {
                    cout << *v << " ";
                } cout << endl;
                graph.print();
                exit(0);
            }
            cDomain.erase(it);
            deleted = true;
        }
        if (!cDomain.empty()) {
            //bool flag = false;
            cDomain2.clear();
            //rmv.push_back(i);
            for (vector<int>::iterator v=cDomain.begin();v != cDomain.end();v++) {
                pair<int, int> edge = mapto(i, *v);
                if (!graph.removeEdge(edge.first, edge.second)) {
                    cDomain2.push_back(*v);
                }
            }
            for (vector<int>::iterator v=cDomain2.begin();v != cDomain2.end();v++) {
                pair<int, int> edge = mapto(i, *v);
                vector<Cost> weight = graph.getWeight(edge.second, edge.first);
                if (weight.size() == 0) {
                    cout << "error for non-existence of edge (" << edge.second << "," << edge.first << ")\n";
                    graph.print();
                    exit(0);
                }
                result = graph.augment(edge.first, edge.second, true);
                if (result.second) {
                    //flag = true;
                    cost += weight[0]+result.first;
                    result.second = graph.removeEdge(edge.first, edge.second);
                }
                if (!result.second) {
                    cout << "ERROR cannot delete edge (" << edge.second << "," << edge.first << ")\n";
                    graph.print();
                    exit(0);
                }
            }
            if (cost > 0) graph.removeNegativeCycles(cost);
            deleted = true;
        }
    }
    if (deleted) {
        for (int i=0;i<graph.size() && (zeroEdges != NULL);i++) {
            for (int j=0;j<graph.size();j++) {
                zeroEdges[i][j] = false;
            }
        }
    }

}

void FlowBasedGlobalConstraint::findProjection(Graph &graph, StoreCost &cost, int varindex, map<Value, Cost> &delta) {

    //if (ToulBar2::GCLevel == LC_NC) return;

    pair<Cost, bool> result;
    delta.clear();
    EnumeratedVariable* x = (EnumeratedVariable*)getVar(varindex);
    for (EnumeratedVariable::iterator j = x->begin(); j != x->end(); ++j) {
        pair<int,int> edge = mapto(varindex, *j);
        Cost tmp = cost;
        //vector<Cost> weight = graph.getWeight(edge.first, edge.second);
        //if (!weight.empty()) {
        if (graph.edgeExist(edge.first, edge.second)) {
            if (zeroEdges[edge.first][edge.second]) {
                //cout << "good\n";
                tmp = cost;
            } else {
                vector<pair<int, int> > edges;
                result = graph.augment(edge.second, edge.first, false, edges);
                /*if (!result.second) {
				  printf("error! no shortest path\n");
				  exit(0);
				  }*/
                //tmp = cost+result.first+weight[0];
                tmp = cost+result.first + graph.getMinWeight(edge.first, edge.second);
                zeroEdges[edge.first][edge.second] = true;
                for (vector<pair<int,int> >::iterator i = edges.begin();i !=
                        edges.end();i++) {
                    zeroEdges[i->first][i->second] = true;
                }
            }
        }
        assert(tmp >= 0);
        delta[*j] = tmp;
    }

}

void FlowBasedGlobalConstraint::augmentStructure(Graph &graph, StoreCost &cost, int varindex, map<Value, Cost> &delta) {

    for (map<Value, Cost>::iterator i = delta.begin(); i != delta.end();i++) {
        pair<int,int> edge = mapto(varindex, i->first);
        if (!graph.increaseCost(edge.first, edge.second, -i->second)) {
            graph.increaseCost(edge.second, edge.first, i->second);
            cost -= i->second;
        }
    }

}

void FlowBasedGlobalConstraint::changeAfterExtend(vector<int> &supports, vector<map<Value, Cost> > &deltas){	

    for (unsigned int i=0;i<supports.size();i++) {
        for (map<Value, Cost>::iterator v = deltas[i].begin();v != deltas[i].end();v++)
            v->second *= -1;
        augmentStructure(*graph, cost, supports[i], deltas[i]);
        for (map<Value, Cost>::iterator v = deltas[i].begin();v != deltas[i].end();v++)
            v->second *= -1;
    }
    graph->removeNegativeCycles(cost);
    for (int i=0;i<graph->size();i++) {
        for (int j=0;j<graph->size();j++) {
            zeroEdges[i][j] = false;
        }
    }
}

void FlowBasedGlobalConstraint::changeAfterProject(vector<int> &supports, vector<map<Value, Cost> > &deltas){

    for (unsigned int i=0;i<supports.size();i++) {
        augmentStructure(*graph, cost, supports[i], deltas[i]);
    }
    graph->removeNegativeCycles(cost);

}

void FlowBasedGlobalConstraint::getDomainFromGraph(Graph &graph, int varindex, vector<int> &domain) {

    domain.clear();
    for (map<Value, int>::iterator i = mapval.begin();i != mapval.end();i++) {
        if ((graph.edgeExist(i->second, varindex+1)) ||  (graph.edgeExist(varindex+1, i->second)))  {
            domain.push_back(i->first);
        }
    }

}


/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

