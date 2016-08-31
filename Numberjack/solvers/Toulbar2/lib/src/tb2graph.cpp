#include "tb2graph.hpp"

Graph::Graph(int n, int depth_)  : adjlist(n)
        , vertexList(n)
        //, potential(n, StoreCost(0))
        , p(n)
        , counter(n)
        , d(n)
        , gsize(n)
        , depth(depth_)
        , intDLinkStore(n*n)
{
    for (int i=0;i<n;i++) vertexList[i] = new Vertex(n, depth_, &intDLinkStore);
}

Graph::~Graph() {
    for (int i=0;i<gsize;i++) {
        delete vertexList[i];
        for (vector<List_Node*>::iterator it = adjlist[i].begin();
                it != adjlist[i].end(); it++) {
            delete *it;
        }
    }
}

int Graph::addEdgeInternal(int u, int v, Cost w, Cost capacity, int tag, 
        bool addReverse, int index) {

    if ((u < 0) || (u >= size())) return -1;
    if ((v < 0) || (v >= size())) return -1;

    int eIndex = -1;
    if (tag == NO_TAG) {
        adjlist[u].push_back(new List_Node(depth, v, w, capacity, tag, index));
        eIndex = adjlist[u].size() - 1;
        if (capacity > 0) {
            vertexList[u]->edgeList[v]->push_back(adjlist[u].size()-1);
            if (vertexList[u]->edgeList[v]->size() == 1) {
                vertexList[u]->neighbor.push_back(v);
            }
        }
    } else {
        bool exist = false;
        for (BTListWrapper<int>::iterator i = vertexList[u]->edgeList[v]->begin(); i !=
                vertexList[u]->edgeList[v]->end() && !exist;++i) {
            List_Node &edge = *(adjlist[u][*i]);
            if (edge.tag == tag) exist = true;
        }

        if (!exist) {
            adjlist[u].push_back(new List_Node(depth, v, w, capacity, tag, index));
            eIndex = adjlist[u].size() - 1;
            if (capacity > 0) {
                vertexList[u]->edgeList[v]->push_back(adjlist[u].size()-1);
                if (vertexList[u]->edgeList[v]->size() == 1) {
                    vertexList[u]->neighbor.push_back(v);
                }
            }
        }
    }

    if (addReverse) {
        int rEdgeIndex = addEdgeInternal(v, u, -w, 0, tag, false, eIndex);
        adjlist[u][eIndex]->rEdgeIndex = rEdgeIndex;
    }

    return eIndex;
}

bool Graph::removeEdge(int u, int v, int tag) {

    bool exist = false;

    if ((u < 0) || (u >= size())) return exist;
    if ((v < 0) || (v >= size())) return exist;

    for (BTListWrapper<int>::iterator i = vertexList[u]->edgeList[v]->begin(); i !=
            vertexList[u]->edgeList[v]->end() && !exist;++i) {
        List_Node &edge = *(adjlist[u][*i]);
        if ((tag == NO_TAG) || (tag == edge.tag)) {
            edge.cap = 0;
            vertexList[u]->edgeList[v]->erase(i);
            if (vertexList[u]->edgeList[v]->size() == 0) {
                vertexList[u]->neighbor.remove(v);
            }
            exist = true;
            break;
        }
    }

    return exist;

}

bool Graph::modifyCost(int u, int v, Cost cost, int tag) {		

    bool exist = false;

    if ((u < 0) || (u >= size())) return exist;
    if ((v < 0) || (v >= size())) return exist;

    //int originalWeight = -1;
    int rEdgeIndex = -1;

    for (BTListWrapper<int>::iterator i = vertexList[u]->edgeList[v]->begin(); i !=
            vertexList[u]->edgeList[v]->end();++i) {
        List_Node &edge = *(adjlist[u][*i]);
        if ((tag == NO_TAG) || (tag == edge.tag)) {
            //originalWeight = edge.weight;
            rEdgeIndex = edge.rEdgeIndex;
            edge.weight = cost;
            exist = true;
            break;
        }
    }

    if (exist) {
        if (rEdgeIndex >= 0) {
            List_Node &edge = *(adjlist[v][rEdgeIndex]);
            if (edge.cap == 0) {
                edge.weight = -cost;
            }
        }
    }

    return exist;

}

bool Graph::increaseCost(int u, int v, Cost cost, int tag) {

    bool exist = false;

    if ((u < 0) || (u >= size())) return exist;
    if ((v < 0) || (v >= size())) return exist;

    //int originalWeight = -1;
    int rEdgeIndex = -1;
    for (BTListWrapper<int>::iterator i = vertexList[u]->edgeList[v]->begin(); i !=
            vertexList[u]->edgeList[v]->end();++i) {
        List_Node &edge = *(adjlist[u][*i]);
        if ((tag == NO_TAG) || (tag == edge.tag)) {
            //originalWeight = edge.weight;
            rEdgeIndex = edge.rEdgeIndex;
            edge.weight += cost;
            exist = true;
            break;
        }
    }

    if (exist) {
        if (rEdgeIndex >= 0) {
            List_Node &edge = *(adjlist[v][rEdgeIndex]);
            if (edge.cap == 0) {
                edge.weight -= cost;
            }
        }
    }

    return exist;

}

bool Graph::edgeExist(int u, int v) {
    return !(vertexList[u]->edgeList[v]->empty());
}

vector<Cost> Graph::getWeight(int u, int v, int tag) {

    vector<Cost> weight;
    for (BTListWrapper<int>::iterator i = vertexList[u]->edgeList[v]->begin(); i !=
            vertexList[u]->edgeList[v]->end();++i) {
        List_Node &edge = *(adjlist[u][*i]);
        if ((tag == NO_TAG) || (tag == edge.tag)) {
            weight.push_back(edge.weight);
        }
    }

    return weight;

}

Cost Graph::getMinWeight(int u, int v, int tag) {

    Cost minWeight = MAX_COST+2;
    for (BTListWrapper<int>::iterator i = vertexList[u]->edgeList[v]->begin(); i !=
            vertexList[u]->edgeList[v]->end();++i) {
        List_Node &edge = *(adjlist[u][*i]);
        if ((tag == NO_TAG) || (tag == edge.tag)) {
            minWeight = min(minWeight, (Cost)edge.weight);
        }
    }

    return minWeight;
}

void Graph::addFlow(int u, int v, Cost flowval) {	

    Cost target = getMinWeight(u,v);
    for (BTListWrapper<int>::iterator i = vertexList[u]->edgeList[v]->begin(); i !=
            vertexList[u]->edgeList[v]->end();++i) {
        List_Node &edge = *(adjlist[u][*i]);
        if (edge.weight == target) {
            edge.cap -= flowval;
            if (edge.cap <= 0) {
                vertexList[u]->edgeList[v]->erase(i);
                if (vertexList[u]->edgeList[v]->size() == 0) {
                    vertexList[u]->neighbor.remove(v);
                }
            }
            assert(edge.rEdgeIndex >= 0);
            List_Node &redge = *(adjlist[v][edge.rEdgeIndex]);
            if ((redge.cap == 0) && (flowval > 0)) {
                vertexList[v]->edgeList[u]->push_back(edge.rEdgeIndex);
                if (vertexList[v]->edgeList[u]->size() == 1) {
                    vertexList[v]->neighbor.push_back(u);
                }
            }
            redge.cap += flowval;
            break;
        }
    }

}

pair<int, Cost> Graph::minCostFlow(int s, int t) {

    //int n = size();
    Cost flow = 0;
    Cost cost = 0;
    pair<Cost, bool> result;
    bool stopped = false;

    //shortest_path(s);
    //for (int i=0;i<n;i++) potential[i] = d[i];
    int iterationCount = 0;
    while (!stopped) {
        iterationCount++;
        stopped = false;
        //shortest_path_with_potential(s);
        shortest_path(s);
        int u = t;
        Cost minc = MAX_COST+1;
        while (p[u] != u && !stopped) {
            int v = p[u];
            if (v < 0) {
                stopped = true;
                break;
            } else {
                Cost minw = MAX_COST+2;
                for (BTListWrapper<int>::iterator j = vertexList[v]->edgeList[u]->begin(); j !=
                        vertexList[v]->edgeList[u]->end();++j) {
                    List_Node &edge = *(adjlist[v][*j]);
                    if ((minw > edge.weight) && (minc > edge.cap)) minc = edge.cap;
                }
                u = v;
            }
        }

        if (!stopped) {
            flow += minc;
            cost += d[t]*minc;
            int u = t;
            while (p[u] != u) {
                addFlow(p[u], u, minc);
                u = p[u];
            }
            //for (int i=0;i<n;i++) potential[i] += d[i];
        }
    }

    return make_pair(flow, cost);
}

pair<Cost, bool> Graph::augment(int s, int t, bool can_change) {

    shortest_path(s);
    pair<Cost, bool> result;
    result.first = d[t];
    result.second = p[t] != -1;

    int u = t;
    Cost minc = INF+1;
    bool exist = false;

    Cost minw = MAX_COST+3;
    for (BTListWrapper<int>::iterator i = vertexList[t]->edgeList[s]->begin(); i !=
            vertexList[t]->edgeList[s]->end();++i) {
        List_Node &edge = *(adjlist[t][*i]);
        assert(edge.cap != 0);
        if (minw >= edge.weight) {
            minc = edge.cap;
            exist = true;
        }
    }
    int count = 0;

    assert(minc > 0);

    std::ostringstream ss;
    while (p[u] != u && result.second) {
        int v = p[u];
        count++;
        if (count >= size()+4) break;
        if (v < 0) {
            result.second = false;
            break;
        } else {
            Cost minw = MAX_COST+2;
            for (BTListWrapper<int>::iterator i = vertexList[v]->edgeList[u]->begin(); i !=
                    vertexList[v]->edgeList[u]->end();++i) {
                List_Node &edge = *(adjlist[v][*i]);
                if (minw >= edge.weight) {
                    if (minc > edge.cap) minc = edge.cap;
                }
            }
            u = v;
        }
    }

    result.first = d[t]*minc;
    result.second = p[t] != -1;

    assert(count < size() + 4);
    assert(minc > 0);

    if (can_change && result.second) {
        u = t;
        while (p[u] != u) {
            addFlow(p[u], u, minc);
            u = p[u];
        }
        if (exist) addFlow(t, s, minc);
        //for (int i=0;i<size();i++) potential[i] = d[i];
    }

    return result;

}

void Graph::removeNegativeCycles(StoreCost &cost) {

    int n = size();
    int pass[n];

    while (1) {

        bool nevloop = false;
        list<int> Q;
        for (int i=0;i<n;i++) Q.push_back(i);
        shortest_path(Q, nevloop);

        if (!nevloop) break;

        int t = -1;
        for (int i=0;i<n;i++) {
            if (counter[i] > n) {
                t = i;
            }
        }
        for (int i=0;i<n;i++) pass[i] = 0;
        int s = p[t], u = s, v = s;
        vector<int> path;
        while (pass[u] == 0){
            pass[u] = 1;
            path.push_back(u);
            u = p[u];
        } path.push_back(u);

        t = path.back();
        for (vector<int>::iterator i = path.begin(); i != path.end();i++) {
            if (*i == t) {
                path.erase(path.begin(), i+1);
                break;
            }
        }

        reverse(path.begin(), path.end());
        path.push_back(*(path.begin()));

        Cost weight = 0, minc = INF;

        for (vector<int>::iterator i = path.begin(); i != path.end()-1;i++) {
            u = *i;  v = *(i+1);
            Cost w = INF;
            Cost c = INF;
            for (BTListWrapper<int>::iterator j = vertexList[u]->edgeList[v]->begin(); j !=
                    vertexList[u]->edgeList[v]->end();++j) {
                List_Node &edge = *(adjlist[u][*j]);
                if (edge.weight < w) {
                    w = edge.weight;
                    c = edge.cap;
                }
            }
            weight += w;
            if (minc > c) minc = c;
        }

        if (weight < 0) {
            for (vector<int>::iterator i = path.begin(); i != path.end()-1;i++) {
                u = *i;  v = *(i+1);
                Cost w = INF;
                for (BTListWrapper<int>::iterator j = vertexList[u]->edgeList[v]->begin(); j !=
                        vertexList[u]->edgeList[v]->end();++j) {
                    List_Node &edge = *(adjlist[u][*j]);
                    if (edge.weight < w) {
                        w = edge.weight;
                    }
                }
                addFlow(u, v, minc);
            }
            cost += minc*weight;
        }

    }

}

void Graph::print(ostream &os) {

    for (int u=0;u<gsize;u++) {
        os << u << ": ";
        for (BTListWrapper<int>::iterator j = vertexList[u]->neighbor.begin(); j !=
                vertexList[u]->neighbor.end();++j) {
            os << *j << "(" << vertexList[u]->edgeList[*j]->size() << ") ";
        }
        os << "\n";
    }

    /*os << "==potential==\n";

	for (int u=0;u<gsize;u++) {
		os << u << ": " << potential[u] << " ";
	}
	os << "\n";
     */
    os << "==graph===\n";

    for (int i = 0; i < gsize; i++) {
        os << i << ":";
        for (vector<List_Node*>::iterator ptr = adjlist[i].begin(); ptr != adjlist[i].end(); ++ptr) {
            List_Node* j = *ptr;
            if (j->cap > 0) {
                if (j->tag != NO_TAG) {
                    os << "(" << j->adj << "," << j->weight << "," << j->cap << "," << j->tag << ") ";
                } else {
                    os << "(" << j->adj << "," << j->weight << "," << j->cap << ",-) ";
                }
            } else {
                if (j->tag != NO_TAG) {
                    os << "[[" << j->adj << "," << j->weight << "," << j->cap << "," << j->tag << "]] ";
                } else {
                    os << "[[" << j->adj << "," << j->weight << "," << j->cap << ",-]] ";
                }
            }
        }
        os << endl;
    }

}

void Graph::shortest_path(list<int> &sources, bool &nevloop) {

    int n = size();
    for (int i=0;i<n;i++) {
        p[i] = -1;  d[i] = INF; counter[i] = 0;
    }

    nevloop = false;
    list<int> &Q = sources;
    for (list<int>::iterator i = sources.begin();i!=sources.end();i++) {
        counter[*i]++;
        d[*i] = 0;   p[*i] = *i;
    }

    while (!Q.empty()) {
        int u = Q.front();
        if (counter[u] > n+2) {
            nevloop = true;
            break;
        }
        Q.pop_front();
        for (BTListWrapper<int>::iterator j = vertexList[u]->neighbor.begin(); j !=
                vertexList[u]->neighbor.end();++j) {
            BTListWrapper<int> &edgeList = *(vertexList[u]->edgeList[*j]);
            for (BTListWrapper<int>::iterator k = edgeList.begin(); k !=
                    edgeList.end();++k) {
                List_Node &edge = *(adjlist[u][*k]);
                if ((d[u] + edge.weight < d[edge.adj])) {
                    d[edge.adj] = d[u] + edge.weight;
                    p[edge.adj] = u;
                    Q.push_back(edge.adj);
                    counter[edge.adj]++;
                }
            }
        }
    }

}

void Graph::printPath(int s, int t) {

    int u = t;
    cout << u << " <- ";
    while (p[u] != u) {
        cout << p[u] << " <- ";
        u = p[u];
    }
    cout << u << " <- " << endl;

}

// Not used due to error in updating potentials after argumentation
/*
void Graph::shortest_path_with_potential(int s) {

	int n = size();
	for (int i=0;i<n;i++) {
		p[i] = -1;  
		d[i] = INF; 
		counter[i] = 0; 
	}
	d[s] = 0;   p[s] = s;

	priority_queue<pair<Cost, int> > Q;
	Q.push(make_pair(0, s));	
	for (int i=0;i<n;i++) {
		int u = -1;
		while (!Q.empty()) {
			pair<Cost, int> top = Q.top();
			Q.pop();	
			if (counter[top.second] == 0) {
				u = top.second;
				break;
			}
		}
		if (u == -1) break;
		counter[u]++;
		for (BTListWrapper<int>::iterator j = vertexList[u]->neighbor.begin(); j !=
						vertexList[u]->neighbor.end();++j) {
			BTListWrapper<int> &edgeList = *(vertexList[u]->edgeList[*j]);
			for (BTListWrapper<int>::iterator k = edgeList.begin(); k !=
				edgeList.end();++k) {		
				List_Node &edge = *(adjlist[u][*k]);
				Cost weight = edge.weight + potential[u] - potential[edge.adj];
				if (weight < 0) {
					cout << "(u,v): " << u << " " << edge.adj << " weight: " <<
					edge.weight << "\n";
				}
				assert(weight >= 0);
				if ((d[u] + weight < d[edge.adj])) {
					d[edge.adj] = d[u] + weight;
					p[edge.adj] = u;
					Q.push(make_pair(-d[edge.adj], edge.adj));
				}
			}
		}
	}	

	for (int i=0;i<n;i++) {
		d[i] = d[i] - potential[s] + potential[i];
	}

}
 */

/*
   set<set<int> >& Graph::compute_scc() {
   if (!p) p = new int[gsize];	
   if (!d) d = new Cost[gsize];	
   if (!low) low = new Cost[gsize];	
   if (!out) out = new bool[gsize];	

   for (int i=0;i<gsize;i++) {
   d[i] = low[i] = 0;
   p[i] = -1;
   out[i] = false;
   }	

   vector<int> stk;
   scc.clear();
   int ts = 0;

   scc_dfs(0, ts, stk);

   return scc;	

   }

   int Graph::scc_dfs(int i, int &ts, vector<int> &stk) {
   p[i] = 0;
   ts++;
   low[i] = d[i] = ts;
   stk.push_back(i);
   for (vector<List_Node >::iterator j = adjlist[i].begin();j !=
   adjlist[i].end();j++) {
   if (p[j->adj] == -1) {
   low[i] = min(low[i], scc_dfs(j->adj, ts, stk));
   } else if (!out[j->adj]) {
   low[i] = min(low[i], d[i]);
   }
   }
   if (low[i] == d[i]) {
   set<int> a;
   while (true) {
   int t = stk.back();
   stk.pop_back();
   out[t] = true;
   a.insert(t);
   if (t == i) break;
   }
   scc.insert(a);
   }
   return low[i];

   }
 */

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

