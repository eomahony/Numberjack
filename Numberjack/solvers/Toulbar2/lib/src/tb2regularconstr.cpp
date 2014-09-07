#include "tb2regularconstr.hpp"
#include "tb2wcsp.hpp"

RegularConstraint::RegularConstraint(WCSP *wcsp, EnumeratedVariable** scope_in, 
		int arity_in) : FlowBasedGlobalConstraint(wcsp, scope_in, arity_in) {
	tempdomain.resize(arity_);
	mapedge.resize(arity_);
	for (int i=0;i<arity_;i++) mapedge[i].resize(((EnumeratedVariable*)getVar(i))->getDomainInitSize());
}

void RegularConstraint::read(istream &file) {
	/*file >> subdef;
	file >> insdef;
	file >> deldef;
	def = subdef;
	if (insdef > def) def = insdef;
	if (deldef > def) def = deldef;*/
	string str;
	file >> str >> def;
	if (str == "var") {
		subdef = def;
		insdef = deldef = 0;
	} else if (str == "edit") {
		subdef = insdef = deldef = def;
	}
	// input the automaton
	int nstate;
	file >> nstate;
	dfa.setSize(nstate);
	int nstart;
	file >> nstart;
	for (int i=0;i<nstart;i++) {
		int t; file >> t;
		dfa.init.push_back(t);
		//cout << t << " ";
	} //cout << endl;
	int nfinal;
	file >> nfinal;
	for (int i=0;i<nfinal;i++) {
		int t; file >> t;
		dfa.final.push_back(t);
		//cout << t << " ";
	} //cout << endl;
	int ntransition;
	file >> ntransition;
	for (int i=0;i<ntransition;i++) {
		int start, end, symbol;
		file >> start;
		file >> symbol;
		file >> end;
		dfa.symbol.insert(symbol);
		dfa.addTransition(start, end, symbol);
	}
	//dfa.print();
	graph.setSize(nstate*(arity_+1)+2);
	tgraph.setSize(nstate*(arity_+1)+2);

	if (insdef > 0) {
		Graph g;
		g.setSize(dfa.size());
		for (int start=0;start<dfa.size();start++) {
			for (int end=0;end<dfa.size();end++) {
				set<int> sym = dfa.getSymbolNeed(start, end);
				if (sym.size() != 0 && (start != end)) g.addEdge(start, end, insdef);
			}
		}
		table.resize(dfa.size());
		for (int s=0;s<dfa.size();s++) {
			g.shortest_path(s, table[s]);
		}
	} 

}

Cost RegularConstraint::eval( String s ) {	
	for (int i=0;i<arity_;i++) {
		tempdomain[i].clear();
		tempdomain[i].insert(s[i]-CHAR_FIRST);
	}
	buildGraphBasic(tgraph, false);
	pair<Cost, bool> result = tgraph.augment(0, tgraph.size()-1, false);
	if (!result.second) result.first = wcsp->getUb();
	return result.first - projectedCost;
}

void RegularConstraint::buildGraph(Graph &g) {
	for (int i=0;i<arity_;i++) {
		EnumeratedVariable* x = (EnumeratedVariable*)getVar(i);
		tempdomain[i].clear();
		for (EnumeratedVariable::iterator v = x->begin();v != x->end();++v) {
			tempdomain[i].insert(*v);
		}
	}
	buildGraphBasic(g, true);
	/*for (int i=0;i<arity_;i++) {
	  EnumeratedVariable* x = (EnumeratedVariable*)getVar(i);
	  cout << i << " : " << endl;
	  for (EnumeratedVariable::iterator v = x->begin();v != x->end();++v) {
	  cout << "   " << *v << " : " << endl;
	  vector<pair<int,int> > edges = mapedge[i][*v];
	  for (vector<pair<int,int> >::iterator j = edges.begin(); j !=
	  edges.end();j++) {
	  cout << "(" << j->first << "," << j->second << ") ";
	  } cout << endl;
	  }
	  }*/
	//g.print();
	//exit(0);
}

void RegularConstraint::buildGraphBasic(Graph &g, bool needRebuildIndex) {

	if (g.size() == 0) g.setSize(dfa.size()*(arity_+1)+2);

	g.clearEdge();
	for (vector<int>::iterator i = dfa.init.begin(); i != dfa.init.end();i++) {
		g.addEdge(0, (*i)+1, 0);
	}
	if (needRebuildIndex) {
		for (int i=0;i<arity_;i++) {
			for (unsigned int j=0;j<mapedge[i].size();j++) mapedge[i][j].clear();
		}
	}

	for (int i=0;i<arity_;i++) {
		for (int start=0;start<dfa.size();start++) {
			for (int end=0;end<dfa.size();end++) {
				set<int> sym = dfa.getSymbolNeed(start, end);
				if (sym.size() != 0) {
					for (set<int>::iterator v = tempdomain[i].begin();v !=
							tempdomain[i].end();v++) {
						Cost w = -deltaCost[i][*v];
						if (sym.find(*v) == sym.end()) w += subdef;
						g.addEdge(i*dfa.size()+start+1, (i+1)*dfa.size()+end+1, w, 1, *v);
						if (needRebuildIndex) mapedge[i][*v].push_back(make_pair(i*dfa.size()+start+1, (i+1)*dfa.size()+end+1));
					}
					if (insdef > 0) g.addEdge(i*dfa.size()+start+1, i*dfa.size()+end+1, insdef, 1, -2);
				}
			}
		}
	}

	if (deldef > 0) {
		for (int i=0;i<arity_;i++) {
			for (int start=0;start<dfa.size();start++) {
				vector<Cost> weight = g.getWeight(i*dfa.size()+start+1, (i+1)*dfa.size()+start+1);
				if (weight.empty()) {
					for (set<int>::iterator v = tempdomain[i].begin();v !=
							tempdomain[i].end();v++) {
						Cost w = -deltaCost[i][*v];
						g.addEdge(i*dfa.size()+start+1, (i+1)*dfa.size()+start+1, deldef + w, 1, *v);
						if (needRebuildIndex) mapedge[i][*v].push_back(make_pair(i*dfa.size()+start+1, (i+1)*dfa.size()+start+1));
					}
				}
			}
		}
		/*for (int i=0;i<arity_;i++) {
		  for (int start=0;start<dfa.size();start++) {
		  set<int> sym = dfa.getSymbolNeed(start, start);
		  if (sym.size() == 0) {
		  for (set<int>::iterator v = tempdomain[i].begin();v !=
		  tempdomain[i].end();v++) {
		  int w = -deltaCost[i][*v];
		  if (sym.find(*v) == sym.end()) w += deldef;
		  g.addEdge(i*dfa.size()+start+1, (i+1)*dfa.size()+start+1, w, 1, *v);
		  if (needRebuildIndex) mapedge[i][*v].push_back(make_pair(i*dfa.size()+start+1, (i+1)*dfa.size()+start+1));
		  }
		  }
		  }
		  }*/
	}

	for (vector<int>::iterator i = dfa.final.begin(); i != dfa.final.end();i++) {
		g.addEdge((arity_)*dfa.size() + (*i) + 1, g.size()-1, 0);
	}

	//g.print();
	//exit(0);

}

void RegularConstraint::findProjection(Graph &graph, Cost cost, int varindex, map<Value, Cost> &delta) {

	pair<int, bool> result;
	EnumeratedVariable* x = (EnumeratedVariable*)getVar(varindex);	

	computeShortestPath(graph, cost);	

	for (EnumeratedVariable::iterator v = x->begin(); v != x->end(); ++v) {
		Cost mincost = INF;
		vector<pair<int, int> > &edges = mapedge[varindex][*v];
		for (vector<pair<int,int> >::iterator i = edges.begin();i !=
				edges.end();i++) {
			pair<int,int> edge = *i;
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
		//cout << *v << delta[*v] << endl;
	} 
	//cout << "====\n";

}

void RegularConstraint::checkRemoved(Graph &graph, Cost &cost, vector<int> &rmv) {

	//bool rebuild = false;
	//bool rebuild = true;
	/*for (int i=0;i<arity_ && !rebuild;i++) {
	  EnumeratedVariable* x = (EnumeratedVariable*)getVar(i);
	  tempdomain[i].clear();
	  for (EnumeratedVariable::iterator v = x->begin();v != x->end();++v) {
	  tempdomain[i].insert(*v);
	  }
	  for (int v=0;v<mapedge[i].size();v++) {
	  if ((tempdomain[i].find(v) == tempdomain[i].end()) &&
	  (mapedge[i][v].size() != 0)) {
	  vector<pair<int,int> > &edges = mapedge[i][v];
	  for (vector<pair<int,int> >::iterator e = edges.begin();e !=
	  edges.end();e++) {
	  vector<int> weight = graph.getWeight(e->first, e->second, v);
	  if (weight.size() != 0) {
	  graph.removeEdge(e->first, e->second, v);
	//cout << "edeg " << e->first << "," << e->second << " of label " << v << " removed" << endl; 
	} else {
	rebuild = true;
	}
	}
	mapedge[i][v].clear();
	}
	}
	}*/
	//if (rebuild) {
	for (int i=0;i<arity_;i++) {
		EnumeratedVariable* x = (EnumeratedVariable*)getVar(i);
		tempdomain[i].clear();
		for (EnumeratedVariable::iterator v = x->begin();v != x->end();++v) {
			tempdomain[i].insert(*v);
		}
		for (unsigned int v=0;v<mapedge[i].size();v++) {
			if ((tempdomain[i].find(v) == tempdomain[i].end()) &&
					(mapedge[i][v].size() != 0)) {
				rmv.push_back(i);
				break;
			}
		}
	}
	buildGraphBasic(graph, true);
	//computeShortestPath(graph, cost);	
	//cost = constructFlow(graph);
	//}

	//cout << "No. of SCC = " << graph.compute_scc().size() << endl;

}

void RegularConstraint::augmentStructure(Graph &g, Cost &cost, int varindex, map<Value, Cost> &delta) {

	EnumeratedVariable* x = (EnumeratedVariable*)getVar(varindex);
	for (EnumeratedVariable::iterator v = x->begin();v != x->end();++v) {
		vector<pair<int,int> > &edges = mapedge[varindex][*v];
		for (vector<pair<int,int> >::iterator e = edges.begin();e !=
				edges.end();e++) {
			vector<Cost> weight = graph.getWeight(e->first, e->second, *v);
			if (weight.size() != 0) {
				g.increaseCost(e->first, e->second, -delta[*v], *v);
			} 
			/*else {
			  g.increaseCost(e->second, e->first, delta[*v], *v);
			  cost -= delta[*v];
			  }*/
		}
	}

	//graph.removeNegativeCycles(cost);
	//computeShortestPath(g, cost);

}

void RegularConstraint::computeShortestPath(Graph &g, Cost &cost) {

	fromSource.resize(graph.size());
	toSink.resize(graph.size());
	for (unsigned int i=0;i<fromSource.size();i++) fromSource[i] = toSink[i] = INF;
	for (vector<int>::iterator i = dfa.init.begin(); i != dfa.init.end();i++) {
		fromSource[(*i)+1] = MIN_COST;
	}
	fromSource[0] = MIN_COST;
	for (int i=0;i<arity_+1;i++) {
		/*if (insdef > 0) {
		  for (int j=0;j<dfa.size();j++) {
		  int s = i*dfa.size()+j+1;
		  for (int k=0;k<dfa.size();k++) {
		  int e = i*dfa.size()+k+1;
		  if (fromSource[e] > fromSource[s] + table[j][k]) {
		  fromSource[e] = fromSource[s] + table[j][k];
		  }
		  }
		  }
		  }
		  for (int j=0;j<dfa.size();j++) {
		  int s = i*dfa.size()+j+1;
		  for (vector<List_Node >::iterator v = graph[s].begin();v !=
		  graph[s].end();v++){
		  if ((v->tag != -2) && (fromSource[v->adj] > fromSource[s] + v->weight)) {
		  fromSource[v->adj] = fromSource[s] + v->weight;
		  }
		  }
		  }*/
		bool change = true;
		while (change) {
			change = false;
			for (int j=0;j<dfa.size();j++) {
				int s = i*dfa.size()+j+1;
				for (Graph::iterator v = graph.begin(s);v !=
						graph.end(s);++v){
					if (fromSource[v.adjNode()] > fromSource[s] + v.weight()) {
						fromSource[v.adjNode()] = fromSource[s] + v.weight();
						change = true;
					}
				}
			}
		}
	}
	toSink[graph.size()-1] = MIN_COST;
	for (vector<int>::iterator i = dfa.final.begin(); i != dfa.final.end();i++) {
		toSink[dfa.size()*(arity_)+(*i)+1] = MIN_COST;
	}
	if (insdef > 0) {
		for (int j=0;j<dfa.size();j++) {
			int s = (arity_)*dfa.size()+j+1;
			for (int k=0;k<dfa.size();k++) {
				int e = (arity_)*dfa.size()+k+1;
				if (toSink[e] > toSink[s] + table[j][k]) {
					toSink[e] = toSink[s] + table[j][k];
				}
			}
		}
	}
	for (int i=arity_-1;i>=0;i--) {
		bool change = true;
		while (change) {
			change = false;
			for (int j=0;j<dfa.size();j++) {
				int s = i*dfa.size()+j+1;
				for (Graph::iterator v = graph.begin(s);v !=
						graph.end(s);++v){
					if (toSink[s] > toSink[v.adjNode()] + v.weight()) {
						toSink[s] = toSink[v.adjNode()] + v.weight();
						change = true;
					}
				}
			}
		}
		/*for (int j=0;j<dfa.size();j++) {
		  int s = i*dfa.size()+j+1;
		  for (vector<List_Node >::iterator v = graph[s].begin();v !=
		  graph[s].end();v++){
		  if ((v->tag != -2) && (toSink[s] > toSink[v->adj] + v->weight)) {
		  toSink[s] = toSink[v->adj] + v->weight;
		  }
		  }
		  }
		  if (insdef > 0) {
		  for (int j=0;j<dfa.size();j++) {
		  int s = i*dfa.size()+j+1;
		  for (int k=0;k<dfa.size();k++) {
		  int e = i*dfa.size()+k+1;
		  if (toSink[e] > toSink[s] + table[j][k]) {
		  toSink[e] = toSink[s] + table[j][k];
		  }
		  }
		  }
		  }*/
	}
	cost = graph.augment(0, graph.size()-1, false).first;

}


// void RegularConstraint::dump(ostream& os, bool original)
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

void RegularConstraint::print(ostream& os)
{
  os << "sregular(";
  for(int i = 0; i < arity_;i++) {
	os << scope[i]->wcspIndex;
	if(i < arity_-1) os << ",";
  }
  os << ")[" << subdef << "," << insdef << "," << deldef << "]";
}
