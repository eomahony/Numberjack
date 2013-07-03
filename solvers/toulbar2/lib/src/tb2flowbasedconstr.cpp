#include "tb2flowbasedconstr.hpp"
#include "tb2wcsp.hpp"

#define verify

FlowBasedGlobalConstraint::FlowBasedGlobalConstraint(WCSP *wcsp, EnumeratedVariable** scope_in, 
		int arity_in) : GlobalConstraint(wcsp, scope_in, arity_in, 0), zeroEdges(NULL) {
}

Cost FlowBasedGlobalConstraint::constructFlow(Graph &g) {
	pair<int, Cost> result = g.minCostFlow(0, g.size()-1);
	return result.second;
}

void FlowBasedGlobalConstraint::initStructure() {

	if (bugraph.size() == 0) bugraph.setSize(graph.size());	
	
	if (zeroEdges == NULL) {
		zeroEdges = new bool*[graph.size()];
		for (int i=0;i<graph.size();i++) zeroEdges[i] = new bool[graph.size()];
	}
	for (int i=0;i<graph.size();i++) {
		for (int j=0;j<graph.size();j++) {	
			zeroEdges[i][j] = false;
		}
	}

	buildGraph();
	cost = constructFlow();
	propagate();

}

void FlowBasedGlobalConstraint::end() {

	if (deconnected()) return;

#ifdef verify
	switch (ToulBar2::LcLevel) {
		case LC_EDAC :
			/*if (!isEDGAC()) {
				cout << "The constraint is not DGAC !\n";
				cout << "cost = " << cost << endl;
				graph.print();
				//exit(0);
			}*/
		case LC_DAC :
			if (!isFDGAC()) {
				cout << "The constraint is not DGAC !\n";
				cout << "cost = " << cost << endl;
				graph.print();
				//GlobalConstraint::propagate();
				exit(0);
			}
			break;
		case LC_FDAC :
			if (!isFDGAC()) {
				cout << "The constraint is not DGAC !\n";
				cout << "cost = " << cost << endl;
				graph.print();
				//GlobalConstraint::propagate();
				exit(0);
			}
		case LC_AC :
			if (!isGAC()) {
				cout << "The constraint is not GAC !\n";
				cout << "cost = " << cost << endl;
				graph.print();
				cout << "++++" << endl;
				for (int i=0;i<arity_;i++) {			
					EnumeratedVariable* y = (EnumeratedVariable*)getVar(i);
					for (EnumeratedVariable::iterator v = y->begin(); v != y->end();++v) {
						cout << *v << " ";
					} cout << endl;
				}
				//GlobalConstraint::propagate();
				exit(0);
			}
			break;
		case LC_SNIC :
			if (!isStrongNIC()) {
				cout << "The constraint is not strong NIC !\n";
				cout << "cost = " << cost << endl;
				graph.print();
				exit(0);
			}
		default:
			break;
	}
#endif

}

bool FlowBasedGlobalConstraint::isStrongNIC() {

	Cost ccost;
	map<Value, Cost> delta;

	buildGraph(bugraph);
	ccost = constructFlow(bugraph);

	if (ccost - projectedCost > 0) {
		cout << ccost << ", " << projectedCost << endl;
		cout << "reduced cost > 0!\n";
		return false;
	}

	bool isNIC = true;
	for (int index=0;index<arity_ && isNIC;index++) {
		EnumeratedVariable* x = (EnumeratedVariable*)getVar(index);
		findProjection(bugraph, ccost, index, delta);
		for (EnumeratedVariable::iterator j = x->begin(); j != x->end(); ++j) {
			if (wcsp->getLb()+ delta[*j] + x->getCost(*j) - projectedCost >= wcsp->getUb()) {
				cout << "value " << *j << " in variable " << x->getName() << " should be removed" << endl;

				cout << "delta = " << delta[*j] << endl;
				cout << "unary cost = " << x->getCost(*j) << endl;
				cout << "lb = " << wcsp->getLb() << endl;
				cout << projectedCost << endl;

				//cout << "sum = " << /*wcsp->getLb()+*/ x->getCost(*j) + delta[*j] << endl;
				cout << "ub = " << wcsp->getUb() << endl;
				cout << "correspond to edge (" << mapto(index, *j).first << "," << mapto(index, *j).second << ")" << endl;
				isNIC = false;
			}
		}
	}
	return isNIC;
}

bool FlowBasedGlobalConstraint::isGAC() {

	Cost ccost;
	map<Value, Cost> delta;

	buildGraph(bugraph);
	ccost = constructFlow(bugraph);

	bool isGAC = true;
	for (int index=0;index<arity_ && isGAC;index++) {
		EnumeratedVariable* x = (EnumeratedVariable*)getVar(index);
		if (x->unassigned()) {
			delta.clear();
			findProjection(bugraph, ccost, index, delta);
			for (EnumeratedVariable::iterator j = x->begin(); j != x->end(); ++j) {
				if (delta[*j] > 0) {
					cout << "Variable " << index << " is not GAC*" << endl;
					cout << "A cost of " << delta[*j] << " should be projected to the value " << *j << " in variable " << index << endl;
					pair<int, int> edge = mapto(index, *j);
					//bugraph.augment(*j, index, false);
					bugraph.printPath(edge.second, edge.first);
					isGAC = false;
				}
			}
		}
	}
	return isGAC;

}

bool FlowBasedGlobalConstraint::isFDGAC() {

	Cost ccost;
	map<Value, Cost> delta;

	bool isDGAC = true;
	for (int index=0;index<arity_ && isDGAC;index++) {
		EnumeratedVariable* x = (EnumeratedVariable*)getVar(index);
		if (x->unassigned()) {
			buildGraph(bugraph);
			for (int j=index+1;j<arity_;j++) {
				EnumeratedVariable* y = (EnumeratedVariable*)getVar(j);
				delta.clear();
				for (EnumeratedVariable::iterator v = x->begin(); v != x->end(); ++v) {
					delta[*v] = -y->getCost(*v);
				}
				augmentStructure(bugraph, ccost, j, delta);	
			}
			ccost = constructFlow(bugraph);
			delta.clear();
			findProjection(bugraph, ccost, index, delta);
			for (EnumeratedVariable::iterator j = x->begin(); j != x->end(); ++j) {
				if (delta[*j] > 0) {
					isDGAC = false;
				}
			}
			if (!isDGAC) cout << "Variable " << index << " is not FDGAC*\n";
		}
	}
	if (!isDGAC) {
		bugraph.print();
		cout << "++++++++++++++++++++++\n";
	}
	return isDGAC;

}

bool FlowBasedGlobalConstraint::isEDGAC() {
	return true;
}

void FlowBasedGlobalConstraint::checkRemoved(Graph &graph, Cost &cost, vector<int> &rmv) {

	//if (ToulBar2::GCLevel == LC_NC) return;

	pair<Cost, bool> result;
	vector<int> cDomain, cDomain2;
	bool deleted = false;
	for (int i=0;i<arity_;i++) {
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
			bool flag = false;
			cDomain2.clear();
			rmv.push_back(i);
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
					cout << "error for non-existence of egde (" << edge.second << "," << edge.first << ")\n";
					graph.print();
					exit(0);
				}
				result = graph.augment(edge.first, edge.second, true);
				if (result.second) {
					flag = true;
					cost += weight[0]+result.first;
					result.second = graph.removeEdge(edge.first, edge.second);
				}
				if (!result.second) {
					cout << "ERROR cannot delete egde (" << edge.second << "," << edge.first << ")\n";
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

void FlowBasedGlobalConstraint::findProjection(Graph &graph, Cost cost, int varindex, map<Value, Cost> &delta) {

	//if (ToulBar2::GCLevel == LC_NC) return;

	pair<Cost, bool> result;
	delta.clear();
	EnumeratedVariable* x = (EnumeratedVariable*)getVar(varindex);
	for (EnumeratedVariable::iterator j = x->begin(); j != x->end(); ++j) {
		pair<int,int> edge = mapto(varindex, *j);
		Cost tmp = cost;
		vector<Cost> weight = graph.getWeight(edge.first, edge.second);
		if (!weight.empty()) {
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
				tmp = cost+result.first+weight[0];
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

void FlowBasedGlobalConstraint::augmentStructure(Graph &graph, Cost &cost, int varindex, map<Value, Cost> &delta) {

	for (map<Value, Cost>::iterator i = delta.begin(); i != delta.end();i++) {
		pair<int,int> edge = mapto(varindex, i->first);
		if (!graph.increaseCost(edge.first, edge.second, -i->second)) {
			graph.increaseCost(edge.second, edge.first, i->second);
			cost -= i->second;
		}
	}

}

void FlowBasedGlobalConstraint::changeAfterExtend(vector<int> &supports, vector<map<Value, Cost> > &deltas){	

	bucost = cost;
	bugraph.clearEdge();
	for (int i=0;i<bugraph.size();i++) bugraph[i] = graph[i];
	for (unsigned int i=0;i<supports.size();i++) {
		for (map<Value, Cost>::iterator v = deltas[i].begin();v != deltas[i].end();v++)
			v->second *= -1;
		augmentStructure(graph, cost, supports[i], deltas[i]);
		for (map<Value, Cost>::iterator v = deltas[i].begin();v != deltas[i].end();v++)
			v->second *= -1;
	}
	graph.removeNegativeCycles(cost);
	for (int i=0;i<graph.size();i++) {
		for (int j=0;j<graph.size();j++) {	
			zeroEdges[i][j] = false;
		}
	}
}

void FlowBasedGlobalConstraint::changeAfterProject(vector<int> &supports, vector<map<Value, Cost> > &deltas){

	for (unsigned int i=0;i<supports.size();i++) {
		augmentStructure(graph, cost, supports[i], deltas[i]);
	}
	graph.removeNegativeCycles(cost);

}

void FlowBasedGlobalConstraint::getDomainFromGraph(Graph &graph, int varindex, vector<int> &domain) {

	domain.clear();
	for (Graph::iterator k = graph.begin(varindex+1); 
			k != graph.end(varindex+1); ++k) {
		if (k.adjNode() > 0) {
			for (map<Value, int>::iterator i = mapval.begin();i !=
					mapval.end();i++) {
				if (i->second == k.adjNode()) domain.push_back(i->first);
			}
		}
	}
	for (map<Value, int>::iterator i = mapval.begin();i !=
			mapval.end();i++) {
		for (Graph::iterator k = graph.begin(i->second); 
				k != graph.end(i->second); ++k) {
			if (k.adjNode() == varindex+1) {
				domain.push_back(i->first);
			}
		}
	}

}

