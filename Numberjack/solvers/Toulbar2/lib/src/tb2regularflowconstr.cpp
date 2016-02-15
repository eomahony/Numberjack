#include "tb2regularflowconstr.hpp"
#include "tb2wcsp.hpp"

#include <queue>
#include <functional>

RegularFlowConstraint::RegularFlowConstraint(WCSP *wcsp, EnumeratedVariable** scope_in,
        int arity_in) : FlowBasedGlobalConstraint(wcsp, scope_in, arity_in), subdef(0), insdef(0), deldef(0), epsilonChar(0) {

    tempdomain.resize(arity_);
    predomain.resize(arity_);
    curdomain.resize(arity_);

    mapedge.resize(arity_);
    for (int i = 0; i < arity_; i++) mapedge[i].resize(((EnumeratedVariable*) getVar(i))->getDomainInitSize());

    modeEnum["var"] = RegularFlowConstraint::VAR;
    modeEnum["edit"] = RegularFlowConstraint::EDIT;    
}

void RegularFlowConstraint::read(istream &file) {

    string str;
    file >> str >> def;
    /*if (str == "var") {
        subdef = def;
        insdef = deldef = 0;
    } else if (str == "edit") {
        subdef = insdef = deldef = def;
    }*/
    setSemantics(str);
    // input the automaton
    int nstate;
    file >> nstate;
    dfa.setNumStates(nstate);
    int nstart;
    file >> nstart;
    for (int i = 0; i < nstart; i++) {
        int t;
        file >> t;
        dfa.init.push_back(t);
    }
    int nfinal;
    file >> nfinal;
    for (int i = 0; i < nfinal; i++) {
        int t;
        file >> t;
        dfa.final.push_back(t);
    }
    int ntransition;
    file >> ntransition;
    for (int i = 0; i < ntransition; i++) {
        int start, end, symbol;
        file >> start;
        file >> symbol;
        file >> end;
        dfa.symbol.insert(symbol);
        dfa.addTransition(start, symbol, end, 0);
    }

}

void RegularFlowConstraint::organizeConfig() {

    if (mode == RegularFlowConstraint::VAR) {   
        subdef = def;
        insdef = deldef = 0;
    } else if (mode == RegularFlowConstraint::EDIT) {
        subdef = insdef = deldef = def;
    }

    //int nstate = dfa.size();
    //graph.setSize(nstate * (arity_ + 1) + 2);

    if (insdef > 0) {
        /*Graph g;
        g.setSize(dfa.size());
        for (int start = 0; start < dfa.size(); start++) {
            for (int end = 0; end < dfa.size(); end++) {
                set<int> sym = dfa.getSymbolNeed(start, end);
                if (sym.size() != 0 && (start != end)) g.addEdge(start, end, insdef);
            }
        }
        table.resize(dfa.size());
        for (int s = 0; s < dfa.size(); s++) {
            g.shortest_path(s, table[s]);
        } */
        vector<vector<Cost> > weightTable;
        weightTable.resize(dfa.size());
        for (int start = 0; start < dfa.size(); start++) {
            weightTable[start].resize(dfa.size());
            fill_n(weightTable[start].begin(), dfa.size(), MAX_COST);
            for (int end = 0; end < dfa.size(); end++) {
                if (start != end) {
                    set<int> sym = dfa.getSymbolNeed(start, end);
                    if (!sym.empty()) weightTable[start][end] = insdef;
                }
            }
        }

        min_priority_queue<pair<int,Cost> > q;

        table.resize(dfa.size());
        //g.shortest_path(s, table[s]);
        for (int s = 0; s < dfa.size(); s++) {
            vector<Cost> &d = table[s];
            d.resize(dfa.size());
            fill_n(d.begin(), dfa.size(), MAX_COST);
            q.push(make_pair(0,s));
            while (!q.empty()) {
                pair<int,Cost> ele = q.top();
                q.pop();
                for (int i = 0; i < dfa.size(); i++) {
                    Cost tmp = ele.first + d[s];
                    if (d[i] >= tmp) {
                        d[i] = tmp;
                        q.push(make_pair(d[i], i));
                    }
                }
            }
        }		        

    }    

    buildWeightedDFATable();

}

void RegularFlowConstraint::buildWeightedDFATable() {

    int nstate = dfa.size();

    set<Value> sigma;
    for (int i = 0; i < arity_; i++) {
        EnumeratedVariable* x = (EnumeratedVariable*) getVar(i);     
        for (EnumeratedVariable::iterator v = x->begin(); v != x->end(); ++v) {
            sigma.insert(*v);
        }
    }

    epsilonChar = (*(dfa.symbol.begin()))-1;
    costTb.resize(nstate);

    // [start][char][end]
    for (int i=0;i<nstate;i++) {    	
        for (vector<pair<int, int> >::iterator it = dfa.transition[i].begin();
                it != dfa.transition[i].end();it++) {
            for (set<int>::iterator jt = sigma.begin();
                    jt != sigma.end();jt++) {
                costTb[i][*jt][it->second] = subdef;
            }
            costTb[i][it->first][it->second] = 0;
        }
    }

    if (insdef > 0) {
        for (int i=0;i<nstate;i++) {
            for (vector<pair<int, int> >::iterator it = dfa.transition[i].begin();
                    it != dfa.transition[i].end();it++) {
                costTb[i][epsilonChar][it->second] = insdef;
            }
            for (set<int>::iterator jt = sigma.begin();
                    jt != sigma.end();jt++) {
                map<int, Cost>::iterator pos = costTb[i][*jt].find(i);
                if (pos == costTb[i][*jt].end()) {
                    costTb[i][*jt][i] = deldef;
                } else if (pos->second > deldef) {
                    pos->second = deldef;
                }
            }
        }
    }      

}

Cost RegularFlowConstraint::evalOriginal(String s) {

    typedef pair<Cost,pair<int,int> > Element;
    //priority_queue<Element, vector<Element>, greater<Element> > minqueue;
    min_priority_queue<Element> minqueue;
    for (vector<int>::iterator i = dfa.init.begin();i != dfa.init.end();i++) {
        minqueue.push(make_pair(0, make_pair(0, *i)));
    }

    Cost myResult = wcsp->getUb();
    while (!minqueue.empty()) {
        Element ele = minqueue.top();
        minqueue.pop();
        Cost weight = ele.first;
        int curIndex = ele.second.first;
        int curState = ele.second.second;

        if (curIndex == arity_) {
            if (find(dfa.final.begin(), dfa.final.end(), curState) != dfa.final.end()) {
                myResult = 	weight;
                break;
            }
        } else {
            int curValue = s[curIndex] - CHAR_FIRST;
            for (map<int, Cost>::iterator i = costTb[curState][curValue].begin();
                    i != costTb[curState][curValue].end();i++) {
                int nextWeight = weight + i->second;
                int nextState = i->first;
                int nextIndex = curIndex + 1;
                minqueue.push(make_pair(nextWeight, make_pair(nextIndex, nextState)));
            }

            if (insdef > 0) {
                for (map<int, Cost>::iterator i = costTb[curState][epsilonChar].begin();
                        i != costTb[curState][epsilonChar].end();i++) {
                    int nextWeight = weight + i->second;
                    int nextState = i->first;
                    int nextIndex = curIndex;
                    minqueue.push(make_pair(nextWeight, make_pair(nextIndex, nextState)));
                }
            }
        }

    }

    return myResult;

}

size_t RegularFlowConstraint::GetGraphAllocatedSize() {
    return dfa.size()*(arity_ + 1) + 2;
}

void RegularFlowConstraint::buildGraph(Graph &g) {
    for (int i = 0; i < arity_; i++) {
        EnumeratedVariable* x = (EnumeratedVariable*) getVar(i);
        tempdomain[i].clear();
        predomain[i].clear();
        for (EnumeratedVariable::iterator v = x->begin(); v != x->end(); ++v) {
            tempdomain[i].insert(*v);
            predomain[i].insert(*v);
        }
    }
    buildGraphBasic(g, true);  

    /* g.print();

    exit(0);*/

}

void RegularFlowConstraint::buildGraphBasic(Graph &g, bool needRebuildIndex) {

    //if (g.size() == 0) g.setSize(dfa.size()*(arity_ + 1) + 2);

    //g.clearEdge();
    for (vector<int>::iterator i = dfa.init.begin(); i != dfa.init.end(); i++) {
        g.addEdge(0, (*i) + 1, 0, 1, NO_TAG, false);
    }
    if (needRebuildIndex) {
        for (int i = 0; i < arity_; i++) {
            for (unsigned int j = 0; j < mapedge[i].size(); j++) mapedge[i][j].clear();
        }
    }

    for (int i = 0; i < arity_; i++) {
        EnumeratedVariable* x = (EnumeratedVariable*) getVar(i);
        for (int start = 0; start < dfa.size(); start++) {
            for (int end = 0; end < dfa.size(); end++) {
                set<int> sym = dfa.getSymbolNeed(start, end);
                if (sym.size() != 0) {
                    for (set<int>::iterator v = tempdomain[i].begin(); v !=
                            tempdomain[i].end(); v++) {                        
                        Cost w = -deltaCost[i][x->toIndex(*v)];
                        if (sym.find(*v) == sym.end()) w += subdef;
                        g.addEdge(i * dfa.size() + start + 1, (i + 1) * dfa.size() + end + 1, w, 1, *v, false);
                        if (needRebuildIndex) mapedge[i][x->toIndex(*v)].push_back(make_pair(i * dfa.size() + start + 1, (i + 1) * dfa.size() + end + 1));
                    }
                    if (insdef > 0) 
                        g.addEdge(i * dfa.size() + start + 1, i * dfa.size() + end + 1, insdef, 1, INS_TAG, false);
                }
            }
        }
    }

    if (deldef > 0) {
        for (int i = 0; i < arity_; i++) {
            EnumeratedVariable* x = (EnumeratedVariable*) getVar(i);
            for (int start = 0; start < dfa.size(); start++) {
                vector<Cost> weight = g.getWeight(i * dfa.size() + start + 1, (i + 1) * dfa.size() + start + 1);
                if (weight.empty()) {
                    for (set<int>::iterator v = tempdomain[i].begin(); v !=
                            tempdomain[i].end(); v++) {                        
                        Cost w = -deltaCost[i][x->toIndex(*v)];
                        g.addEdge(i * dfa.size() + start + 1, (i + 1) * dfa.size() + start + 1, deldef + w, 1, *v, false);
                        if (needRebuildIndex) mapedge[i][x->toIndex(*v)].push_back(make_pair(i * dfa.size() + start + 1, (i + 1) * dfa.size() + start + 1));
                    }
                }
            }
        }       
    }

    for (vector<int>::iterator i = dfa.final.begin(); i != dfa.final.end(); i++) {
        g.addEdge((arity_) * dfa.size() + (*i) + 1, g.size() - 1, 0, 1, NO_TAG, false);
    }

}

void RegularFlowConstraint::findProjection(Graph &graph, StoreCost &cost, int varindex, map<Value, Cost> &delta) {

    //pair<int, bool> result;
    EnumeratedVariable* x = (EnumeratedVariable*) getVar(varindex);

    computeShortestPath(graph, cost);

    for (EnumeratedVariable::iterator v = x->begin(); v != x->end(); ++v) {
        Cost mincost = INF;
        vector<pair<int, int> > &edges = mapedge[varindex][x->toIndex(*v)];
        for (vector<pair<int, int> >::iterator i = edges.begin(); i !=
                edges.end(); i++) {
            pair<int, int> edge = *i;
            vector<Cost> weight = graph.getWeight(edge.first, edge.second, *v);
            if (weight.size() > 1) cout << "multiple edges?\n";
            if (!weight.empty()) {
                Cost t = weight[0] + fromSource[edge.first] + toSink[edge.second];
                if (mincost > t) mincost = t;
            }
        }
        if (mincost == INF) {
            delta[*v] = wcsp->getUb();
        } else {
            delta[*v] = mincost;
        }     
    }    

}

void RegularFlowConstraint::checkRemoved(Graph &graph, StoreCost &cost, vector<int> &rmv) {		

    //for (int varindex = 0; varindex < arity_; varindex++) {
    for (vector<int>::iterator i = rmv.begin();i != rmv.end();i++)
    {
        int varindex = *i;
        EnumeratedVariable* x = (EnumeratedVariable*) getVar(varindex);
        for (unsigned int valIndex = 0;valIndex < mapedge[varindex].size();valIndex++) {
            if (x->cannotbe(x->toValue(valIndex))) {
                vector<pair<int, int> > &edges = mapedge[varindex][valIndex];
                for (vector<pair<int, int> >::iterator i = edges.begin(); i !=
                        edges.end(); i++) {
                    pair<int, int> edge = *i;
                    graph.removeEdge(edge.first, edge.second, x->toValue(valIndex));
                }
            }
        }
    }

}

void RegularFlowConstraint::augmentStructure(Graph &g, StoreCost &cost, int varindex, map<Value, Cost> &delta) {

    EnumeratedVariable* x = (EnumeratedVariable*) getVar(varindex);
    for (EnumeratedVariable::iterator v = x->begin(); v != x->end(); ++v) {
        vector<pair<int, int> > &edges = mapedge[varindex][x->toIndex(*v)];
        for (vector<pair<int, int> >::iterator e = edges.begin(); e !=
                edges.end(); e++) {
            vector<Cost> weight = g.getWeight(e->first, e->second, *v);
            if (weight.size() != 0) {
                g.increaseCost(e->first, e->second, -delta[*v], *v);
            }
        }
    }

}

void RegularFlowConstraint::computeShortestPath(Graph &g, StoreCost &cost) {

    fromSource.resize(g.size());
    toSink.resize(g.size());
    for (unsigned int i = 0; i < fromSource.size(); i++) fromSource[i] = toSink[i] = INF;
    for (vector<int>::iterator i = dfa.init.begin(); i != dfa.init.end(); i++) {
        fromSource[(*i) + 1] = MIN_COST;
    }
    fromSource[0] = MIN_COST;
    for (int i = 0; i < arity_ + 1; i++) {        
        bool change = true;
        while (change) {
            change = false;
            for (int j = 0; j < dfa.size(); j++) {
                int s = i * dfa.size() + j + 1;                
                for (Graph::node_iterator node = g.node_begin(s); node != g.node_end(s); ++node) {
                    for (Graph::edge_iterator v = g.begin(s, *node); v != g.end(s, *node); ++v) {
                        if (fromSource[v.adjNode()] > fromSource[s] + v.weight()) {
                            fromSource[v.adjNode()] = fromSource[s] + v.weight();
                            change = true;
                        }
                    }
                }
            }
        }        
    }
    toSink[g.size() - 1] = MIN_COST;
    for (vector<int>::iterator i = dfa.final.begin(); i != dfa.final.end(); i++) {
        toSink[dfa.size()*(arity_)+(*i) + 1] = MIN_COST;
    }
    if (insdef > 0) {
        for (int j = 0; j < dfa.size(); j++) {
            int s = (arity_) * dfa.size() + j + 1;
            for (int k = 0; k < dfa.size(); k++) {
                int e = (arity_) * dfa.size() + k + 1;
                if (toSink[e] > toSink[s] + table[j][k]) {
                    toSink[e] = toSink[s] + table[j][k];
                }
            }
        }
    }
    for (int i = arity_ - 1; i >= 0; i--) {
        bool change = true;
        while (change) {
            change = false;
            for (int j = 0; j < dfa.size(); j++) {
                int s = i * dfa.size() + j + 1;              
                for (Graph::node_iterator node = g.node_begin(s); node != g.node_end(s); ++node) {
                    for (Graph::edge_iterator v = g.begin(s, *node); v != g.end(s, *node); ++v) {
                        if (toSink[s] > toSink[v.adjNode()] + v.weight()) {
                            toSink[s] = toSink[v.adjNode()] + v.weight();
                            change = true;
                        }
                    }
                }
            }
        }      
    }
    //cost = graph.augment(0, graph.size() - 1, false).first;

}


// void RegularFlowConstraint::dump(ostream& os, bool original)
// {
//   if (original) {
//     os << arity_;
//     for(int i = 0; i < arity_;i++) os << " " << scope[i]->wcspIndex;
//   } else {
// 	os << nonassigned;
//     for(int i = 0; i < arity_; i++) if (scope[i]->unassigned()) os << " " << scope[i]->getCurrentVarId();
//   }
//   os << " -1 sregular" << endl << ((insdef==0 && deldef==0)?"var":"edit") << " " << def << endl;
//   os << endl;
// }

void RegularFlowConstraint::print(ostream& os) {
    os << "sregular(";
    for (int i = 0; i < arity_; i++) {
        os << scope[i]->wcspIndex;
        if (i < arity_ - 1) os << ",";
    }
    os << ")[" << subdef << "," << insdef << "," << deldef << "]";
}

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

