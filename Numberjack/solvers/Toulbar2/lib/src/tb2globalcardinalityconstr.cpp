#include "tb2globalcardinalityconstr.hpp"
#include "tb2wcsp.hpp"

#define upper_bound first
#define lower_bound second

GlobalCardinalityConstraint::GlobalCardinalityConstraint(WCSP *wcsp, EnumeratedVariable** scope_in, int arity_in) : FlowBasedGlobalConstraint(wcsp, scope_in, arity_in) {
    buildIndex();

    modeEnum["var"] = GlobalCardinalityConstraint::VAR;
    modeEnum["dec"] = GlobalCardinalityConstraint::VALUE;
    modeEnum["wdec"] = GlobalCardinalityConstraint::WVALUE;
}

void GlobalCardinalityConstraint::buildIndex() {
    vector<Value> D;
    mapval.clear();
    for (int i=0;i<arity_;i++) {
        EnumeratedVariable* x = (EnumeratedVariable*)getVar(i);
        for (EnumeratedVariable::iterator iterx = x->begin(); iterx != x->end(); ++iterx) {
            D.push_back(*iterx);
        }
    }
    sort(D.begin(), D.end());
    D.erase(unique(D.begin(), D.end()), D.end());

    for (vector<Value>::iterator i = D.begin(); i != D.end();i++) {
        mapval[*i] = arity_+(int)(i-D.begin())+1;
    }
    nDistinctDomainValue = D.size();
    //graph.setSize(arity_+D.size()+4);
}

void GlobalCardinalityConstraint::read(istream &file) {
    // "var" => softvar
    // "dec" => softdec
    // "wdec" => sigmadec
    string str;
    int nvalues;
    //int sumlow = 0, sumhigh = 0;
    file >> str;
    //JP Start// alteration
    /*mode = EMPTY;
	if (strcmp(str.c_str(), "var") 	== 0) mode = VAR;
	if (strcmp(str.c_str(), "dec") 	== 0) mode = VALUE;
	if (strcmp(str.c_str(), "wdec") == 0) mode = WVALUE;
	if (mode == EMPTY) {
		cerr << "Error occur in reading gcc() : No violation measure" << endl;
		exit(1);
	}*/
    setSemantics(str);
    //JP End//
    file >> def;
    file >> nvalues;
    //JP End//
    for (int i=0;i<nvalues;i++) {
        int d, high, low;
        file >> d >> low >> high;
        if (high < low) {
            cerr << "Error occur in reading gcc()" << endl;
            exit(1);
        }
        //JP Start//
        int wshortage = def;
        int wexcess = def;
        if (mode == WVALUE) {
            file >> wshortage >> wexcess;
        }
        //JP End//
        bound[d] = make_pair(high, low);
        weights[d] = make_pair(wshortage, wexcess);
        //sumlow += low;
        //sumhigh += high;
    }

}

void GlobalCardinalityConstraint::organizeConfig() {       

    int sumlow = 0, sumhigh = 0;

    for (map<Value, pair<int, int> >::iterator i = bound.begin(); i != bound.end();i++) {
        sumlow += i->second.lower_bound;
        sumhigh += i->second.upper_bound;
    }

    for (map<Value,int>::iterator i = mapval.begin();i != mapval.end();i++) {
        if (bound.find(i->first) == bound.end()) {
            bound[i->first] = make_pair(arity_+4, 0);
            sumhigh += arity_+4;
        }
    }

    for (map<Value,int>::iterator i = mapval.begin();i != mapval.end();i++) {
        if (weights.find(i->first) == weights.end()) {
            weights[i->first] = make_pair(def, def);
        }
    }

    if ((mode == VAR) && ((arity_ < sumlow) || (arity_ > sumhigh))) {
        cerr << "Error occur in gcc() model using variable-based measure : " << endl;
        cerr << "sum of lower bound is too high / sum of upper bound is too low\n." << endl;
        cerr << "sum high = " << sumhigh << endl;
        cerr << "arity_ = " << arity_ << endl;
        cerr << "sum low = " << sumlow << endl;
        exit(1);
    }

}

Cost GlobalCardinalityConstraint::evalOriginal( String s ) {

    Cost excess = 0, shortage = 0, cost = 0;
    map<Value ,int> appear;
    for (unsigned int i=0;i<s.length();i++) {
        appear[s[i]-CHAR_FIRST]++;
    }
    for (map<Value, pair<int,int> >::iterator i = bound.begin(); i != bound.end();i++) {
        if (appear[i->first] < i->second.lower_bound) {
            //JP Start// Alteration
            Cost lshortage = i->second.lower_bound - appear[i->first];
            shortage += lshortage;
            cost += lshortage*weights[i->first].first;
            //JP End//
        }
        if (appear[i->first] > i->second.upper_bound) {
            //JP Start// Alteration
            Cost lexcess = appear[i->first] - i->second.upper_bound;
            excess += lexcess;
            cost += lexcess*weights[i->first].second;
            //JP End//
        }
    }
    //JP Start// Alteration
    if (mode == VAR) {
        cost = (excess>shortage)?excess*def:shortage*def;
    }
    //JP End//
    return cost;
}

size_t GlobalCardinalityConstraint::GetGraphAllocatedSize() {
    return arity_+nDistinctDomainValue+4;
}

void GlobalCardinalityConstraint::buildGraph(Graph &g) {

    int n = g.size();
    int t = n-3;

    int ss = n-1;
    int st = n-2;

    //g.clearEdge();
    g.addEdge(t, 0, 0, INF);
    g.addEdge(0, st, 0, arity_);

    for (int i=0;i<arity_;i++) {
        g.addEdge(ss, i+1, 0, 1);
        EnumeratedVariable* x = (EnumeratedVariable*)getVar(i);
        for (EnumeratedVariable::iterator v = x->begin(); v != x->end(); ++v) {
            g.addEdge(i+1, mapval[*v], -deltaCost[i][x->toIndex(*v)]);
        }
    }

    int sumlow = 0;
    for (map<Value, int>::iterator i = mapval.begin(); i != mapval.end();i++) {
        if (bound[i->first].lower_bound != 0) g.addEdge(i->second, st, 0, bound[i->first].lower_bound);
        if (bound[i->first].upper_bound != bound[i->first].lower_bound)
            g.addEdge(i->second, t, 0, bound[i->first].upper_bound-bound[i->first].lower_bound);
        sumlow += bound[i->first].lower_bound;
    }
    if (sumlow > 0) g.addEdge(ss, t, 0, sumlow);

    if (mode == VAR) {
        for (map<Value, int>::iterator i = mapval.begin(); i != mapval.end();i++) {
            for (map<Value, int>::iterator j = mapval.begin(); j != mapval.end();j++) {
                if (i->first != j->first) g.addEdge(i->second, j->second, def, arity_);
            }
        }
    } else {
        for (map<Value, int>::iterator i = mapval.begin(); i != mapval.end();i++) {
            //JP Start// Alteration
            if( bound[i->first].lower_bound > 0) g.addEdge(0, i->second, weights[i->first].first, bound[i->first].lower_bound);
            g.addEdge(i->second, t, weights[i->first].second, arity_);
            //JP End//
        }
    }

}

Cost GlobalCardinalityConstraint::constructFlow(Graph &g) {

    //cout << "use the one\n";
    /*pair<int, bool> result;
	int n = g.size();
	int ss = n-1;
	int st = n-2;
	int fcost = 0;
	//int fcost = -projectedCost;

	do {
		int minc = 0;
		result = g.augment(ss, st, true, minc);
		if (result.second) fcost += minc*result.first;
	} while (result.second);*/
    //checker(g, fcost);

    pair<int, Cost> result = g.minCostFlow(g.size()-1, g.size()-2);
    return result.second;
    //return fcost;
}

/*void GlobalCardinalityConstraint::getDomainFromGraph(Graph &graph, int varindex, vector<int> &domain) {

	domain.clear();
	for (vector<List_Node >::iterator k = graph[varindex+1].begin(); 
			k != graph[varindex+1].end(); k++) {
		if (k->adj > 0) {
			for (map<Value, int>::iterator i = mapval.begin();i !=
					mapval.end();i++) {
				if (i->second == k->adj) domain.push_back(i->first);
			}
		}
	}
	for (map<Value, int>::iterator i = mapval.begin();i !=
			mapval.end();i++) {
		for (vector<List_Node >::iterator k = graph[i->second].begin(); 
				k != graph[i->second].end(); k++) {
			if (k->adj == varindex+1) {
				domain.push_back(i->first);
			}
		}
	}

}*/

void GlobalCardinalityConstraint::dump(ostream& os, bool original) {
    int nvalues = 0;
    if (original) {
        os << arity_;
        for(int i = 0; i < arity_;i++) os << " " << scope[i]->wcspIndex;
    } else {
        os << nonassigned;
        for(int i = 0; i < arity_; i++) if (scope[i]->unassigned()) os << " " << scope[i]->getCurrentVarId();
    }
    for (map<Value, pair<int,int> >::iterator i = bound.begin(); i !=	bound.end();i++) nvalues++;
    os << " -1 sgcc" << " ";
    if (mode == VAR   ) os << "var";
    if (mode == VALUE ) os << "dec";
    if (mode == WVALUE) os << "wdec";
    os << " " << def << " " <<  nvalues << endl;
    for (map<Value, pair<int,int> >::iterator i = bound.begin(); i !=	bound.end();i++) {
        os << i->first << " " << i->second.lower_bound << " " << i->second.upper_bound;
        if (mode == WVALUE) os << " " << weights[i->first].first << " " << weights[i->first].second;
        os << endl;
    }
}

void GlobalCardinalityConstraint::print(ostream& os) {
    int nvalues = 0;

    os << "sgcc(";
    for(int i = 0; i < arity_;i++) {
        os << scope[i]->wcspIndex;
        if(i < arity_-1) os << ",";
    }
    for (map<Value, pair<int,int> >::iterator i = bound.begin(); i !=	bound.end();i++) nvalues++;
    os << ")[" ;
    if (mode == VAR   ) os << "var";
    if (mode == VALUE ) os << "dec";
    if (mode == WVALUE) os << "wdec";
    os << "," << def << "," << nvalues;
    for (map<Value, pair<int,int> >::iterator i = bound.begin(); i !=	bound.end();i++) {
        os << "," << i->first << "," << i->second.lower_bound << "," << i->second.upper_bound;
        if (mode == WVALUE) os << "," << weights[i->first].first << "," << weights[i->first].second;
    }
    os << "]";
}

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

