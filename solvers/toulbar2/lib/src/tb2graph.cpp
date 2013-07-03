#include "tb2graph.hpp"

Graph::Graph(const Graph &g) {
	setSize(g.size());
	for (int i=0;i<gsize;i++) adjlist[i] = g.adjlist[i];
}

void Graph::setSize(int n) {
	gsize = n;
	counter = new int[n];
	p = new int[n];
	d = new Cost[n];
	adjlist = new vector<List_Node>[n];
	for (int i=0;i<n;i++) adjlist[i].resize(n);
	potential = new Cost[n];
	clearEdge();
}	

void Graph::clearEdge() {
	for (int i=0;i<gsize;i++) {
		adjlist[i].clear();
	}
} 

void Graph::addEdge(int u, int v, Cost w, Cost capacity, int tag) {

	if ((u < 0) || (u >= size())) return; 
	if ((v < 0) || (v >= size())) return;

	if (tag == -1) {
		adjlist[u].push_back(List_Node(v, w, capacity, tag));
	} else {
		bool exist = false;
		for (vector<List_Node>::iterator i = adjlist[u].begin(); i !=
				adjlist[u].end() && !exist;i++) {
			if ((i->adj == v) && (i->tag == tag)) exist = true;
		}
		if (!exist) adjlist[u].push_back(List_Node(v, w, capacity, tag));
	}

}

bool Graph::removeEdge(int u, int v, int tag) {

	bool exist = false;

	if ((u < 0) || (u >= size())) return exist; 
	if ((v < 0) || (v >= size())) return exist;

	for (vector<List_Node>::iterator i = adjlist[u].begin(); i !=
			adjlist[u].end();i++) {
		if (i->adj == v) {
			if ((tag == -1) || (tag == i->tag)) {
				adjlist[u].erase(i);
				exist = true;
				break;
			}
		}
	}

	return exist;

}

bool Graph::modifyCost(int u, int v, Cost cost, int tag) {

	bool exist = false;

	if ((u < 0) || (u >= size())) return exist; 
	if ((v < 0) || (v >= size())) return exist;

	for (vector<List_Node>::iterator i = adjlist[u].begin(); i !=
			adjlist[u].end();i++) {
		if ((i->adj == v) && ((tag == -1) || (tag == i->tag))) {
			i->weight = cost;
			exist = true;
			break;
		}
	}

	return exist;

}

bool Graph::increaseCost(int u, int v, Cost cost, int tag) {

	bool exist = false;

	if ((u < 0) || (u >= size())) return exist; 
	if ((v < 0) || (v >= size())) return exist;

	for (vector<List_Node>::iterator i = adjlist[u].begin(); i !=
			adjlist[u].end();i++) {
		if ((i->adj == v) && ((tag == -1) || (tag == i->tag))) {
			i->weight += cost;
			exist = true;
			break;
		}
	}

	return exist;

}

vector<Cost> Graph::getWeight(int u, int v, int tag) {

	vector<Cost> weight;

	for (vector<List_Node>::iterator i = adjlist[u].begin(); i !=
			adjlist[u].end();i++) {
		if (i->adj == v) {
			if ((tag == -1) || (tag == i->tag)) weight.push_back(i->weight);
		}
	}

	return weight;

}

void Graph::addFlow(int u, int v, Cost flowval) {
	vector<Cost> vec = getWeight(u,v);
	Cost target = *min_element(vec.begin(), vec.end());
	for (vector<List_Node>::iterator i = adjlist[u].begin(); i !=
			adjlist[u].end();i++) {
		if ((i->adj == v) && (i->weight == target)) {
			i->cap -= flowval;
			bool found = false;
			for (vector<List_Node>::iterator j = adjlist[v].begin(); j
					!= adjlist[v].end() && !found;j++) {
				if ((j->adj == u) && (j->weight == -(i->weight))  
						&& (i->tag == j->tag) ) {
					j->cap += flowval;
					found = true;
				}
			}
			if (!found) {
				addEdge(v, u, -(i->weight), flowval, i->tag);
			}
			if (i->cap == 0) adjlist[u].erase(i);
			break;
		}
	}
}

pair<int, Cost> Graph::minCostFlow(int s, int t, Cost maxValue) {

	int n = size();
	Cost flow = 0;
	Cost cost = 0;
	pair<Cost, bool> result;
	bool stopped = false;

	for (int i=0;i<gsize;i++) sort(adjlist[i].begin(), adjlist[i].end());
	shortest_path(s);
	for (int i=0;i<n;i++) potential[i] = d[i]; 
	for (Cost i=0;i<maxValue && !stopped;i++) {
		/*stopped = false;
		shortest_path(s);
		int u = t;
		Cost minc = MAX_COST;
		while (p[u] != u && !stopped) {
			int v = p[u];	
			if (v < 0) {
				stopped = true;
				break;
			} else {
				for (vector<List_Node>::iterator j = adjlist[v].begin(); j !=
						adjlist[v].end(); j++) {
					if ((j->adj == u) && (minc > j->cap)) {
						minc = j->cap;
						break;
					}
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
		}*/
		stopped = false;
		shortest_path_with_potential(s);
		int u = t;
		Cost minc = MAX_COST;
		while (p[u] != u && !stopped) {
			int v = p[u];	
			if (v < 0) {
				stopped = true;
				break;
			} else {
				for (vector<List_Node>::iterator j = adjlist[v].begin(); j !=
						adjlist[v].end(); j++) {
					if ((j->adj == u) && (minc > j->cap)) {
						minc = j->cap;
						break;
					}
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
			for (int i=0;i<n;i++) potential[i] = d[i]; 
		}
	}

	return make_pair(flow, cost);

}

pair<Cost, bool> Graph::augment(int s, int t, bool can_change) {

	if (!can_change) {
		for (int i=0;i<gsize;i++) {
			sort(adjlist[i].begin(), adjlist[i].end());
		}
	}

	//shortest_path_with_potential(s);
	shortest_path(s);
	pair<Cost, bool> result;
	result.first = d[t];
	result.second = p[t] != -1;

	//	assert(!nevloop);
	/*if (nevloop) {
	  cout << "negative loop exists from " << s << " to " << t <<endl;
	  print();
	  exit(1);
	  }*/

	int u = t;
	Cost minc = INF+1;
	bool exist = false;
	for (vector<List_Node>::iterator i = adjlist[t].begin(); i !=
			adjlist[t].end() && !exist;i++) {
		if (i->adj == s) {
			minc = i->cap;
			exist = true;
		}
	}
	int count = 0;
	while (p[u] != u && result.second) {
		int v = p[u];	
		count++; 
		if (count >= size()+4) break;	
		if (v < 0) {
			result.second = false;
			break;
		} else {
			for (vector<List_Node>::iterator j = adjlist[v].begin(); j !=
					adjlist[v].end(); j++) {
				if ((j->adj == u) && (minc > j->cap)) {
					minc = j->cap;
					break;
				}
			}
			u = v;
		}
	}

	//pair<Cost, bool> result;
	result.first = d[t]*minc;
	//result.second = p[t] != -1;

	assert(count < size() + 4);
	assert(minc > 0);

	if (can_change && result.second) {
		u = t;
		while (p[u] != u) {
			addFlow(p[u], u, minc);
			u = p[u];
		}
		if (exist) addFlow(t, s, minc); 
		for (int i=0;i<size();i++) potential[i] = d[i]; 
	}

	return result;

}

void Graph::removeNegativeCycles(Cost &cost) {

	int n = size();
	//	int *pass = new int[n];
	int pass[n];

	while (1) { 
		for (int i=0;i<n;i++) {
			p[i] = -1;  d[i] = 0; 
		}
		bool nevloop = false;
		bool changed = true;
		int count = 0;

		while (changed && count < n+2) {
			changed = false;
			count++;
			for (int i=0;i<n;i++) pass[i] = 0;
			for (int i=0;i<n;i++) {
				for (vector<List_Node>::iterator j = adjlist[i].begin(); j !=
						adjlist[i].end(); j++) {
					if (d[i] + j->weight < d[j->adj]) {
						d[j->adj] = d[i] + j->weight;
						p[j->adj] = i;
						pass[j->adj] = 1;
						changed = true;
					}
				}
			}
		}

		nevloop = count >= n+1;

		if (!nevloop) break;

		int t = -1;
		for (int i=0;i<n;i++) {
			if (pass[i] == 1){
				if (t == -1) t = i;
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
			for (vector<List_Node>::iterator j = adjlist[u].begin(); j !=
					adjlist[u].end();j++) {
				if ((j->adj == v) && (j->weight < w)) {
					w = j->weight;
					c = j->cap;
				}
			}
			weight += w;
			if (minc > c) minc = c;
		}

		if (weight < 0) {
			for (vector<int>::iterator i = path.begin(); i != path.end()-1;i++) {
				u = *i;  v = *(i+1);
				Cost w = INF;
				for (vector<List_Node>::iterator j = adjlist[u].begin(); j !=
						adjlist[u].end();j++) {
					if ((j->adj == v) && (j->weight < w)) {
						w = j->weight;
					}
				}
				addFlow(u, v, minc);
			}
			cost += minc*weight;
		}

	}

	//	delete [] pass;
}

void Graph::print() {
	int u;
	u = 0;
	for (vector<List_Node> *i = adjlist; i != adjlist+gsize; i++, u++) {
		cout << u << ":";
		for (vector<List_Node>::iterator j = i->begin(); j != i->end(); j++) {
			cout << "(" << j->adj << "," << j->weight << "," << j->tag << ") ";
		}
		cout << endl;
	}
}

void Graph::shortest_path(int s) {

	int n = size();
	for (int i=0;i<n;i++) {
		p[i] = -1;  d[i] = INF; counter[i] = 0;
	}
	d[s] = 0;   p[s] = s;

	bool nevloop = false;
	vector<int> Q;
	Q.push_back(s);
	counter[s]++;
	while (!Q.empty()) {
		int u = *(Q.begin());
		if (counter[u] > n+2) {
			nevloop = true;
			break;
		}
		Q.erase(Q.begin());
		for (vector<List_Node>::iterator j = adjlist[u].begin(); j !=
				adjlist[u].end(); j++) {
			if ((d[u] + j->weight < d[j->adj])) {
				d[j->adj] = d[u] + j->weight;
				p[j->adj] = u;
				Q.push_back(j->adj);
				counter[j->adj]++;
			}
		}
	}

	if (nevloop) {
		cout << "negative loop exists from " <<endl;
		print();
		exit(1);
	}


	/*pathCost.resize(n);
	  for (int i=0;i<n;i++) pathCost[i] = d[i];*/

}

void Graph::shortest_path_with_potential(int s) {

	int n = size();
	for (int i=0;i<n;i++) {
		p[i] = -1;  d[i] = INF; 
	}
	d[s] = 0;   p[s] = s;

	priority_queue<pair<Cost, int> > Q;
	//	vector<bool> Q(n);
	//	for (int i=0;i<n;i++) Q[i]= false;
	Q.push(make_pair(0, s));
	//Q[s] = true;
	while (!Q.empty()) {
		pair<Cost, int> top = Q.top();
		Q.pop();
		int u = top.second;
		/*int u = -1;
		  int minCost = INF;
		  for (int i=0;i<n;i++) {
		  if (Q[i] && (d[i] < INF)) {
		  u = i;
		  minCost = d[i];
		  }
		  }
		  if (u == -1) break;
		  Q[u] = false;*/
		for (vector<List_Node>::iterator j = adjlist[u].begin(); j !=
				adjlist[u].end(); j++) {
			Cost weight = j->weight + potential[u] - potential[j->adj];
			if ((d[u] + weight < d[j->adj])) {
				d[j->adj] = d[u] + weight;
				p[j->adj] = u;
				//Q[j->adj] = true;
				Q.push(make_pair(d[j->adj], j->adj));
			}
		}
	}

	//pathCost.resize(n);
	for (int i=0;i<n;i++) {
		d[i] = d[i] - potential[s] + potential[i];
		//pathCost[i] = d[i];
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
