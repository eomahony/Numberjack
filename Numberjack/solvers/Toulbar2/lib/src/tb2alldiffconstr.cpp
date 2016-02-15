#include "tb2alldiffconstr.hpp"
#include "tb2wcsp.hpp"

AllDiffConstraint::AllDiffConstraint(WCSP *wcsp, EnumeratedVariable** scope_in, 
        int arity_in) : FlowBasedGlobalConstraint(wcsp, scope_in, arity_in) {
    buildIndex();

    modeEnum["var"] = AllDiffConstraint::VAR;
    modeEnum["dec"] = AllDiffConstraint::DEC;
    modeEnum["decbi"] = AllDiffConstraint::DECBI;
}

void AllDiffConstraint::buildIndex() {
    vector<Value> D;
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
    //graph.setSize(arity_+D.size()+2);
}

void AllDiffConstraint::read(istream &file) {

    // Only two pararmeters for control :

    // 1) the cost measure : "var" = mu_var ; "dec" = mu_dec
    string str;
    file >> str;
    // 2) the cost of the violation edge : assume to be a constant function
    // mapping to def
    file >> def;
    //cout << "str = " << str << endl;
    /*if (str == "var") {
		mode = VAR;
	} else if (str == "dec") {
		mode = DEC;
	} else if (str == "decbi") {
		mode = DECBI;
		decompose();
	} else {
		cout << "unknown mode?\n";
		exit(0);
	}*/
    setSemantics(str);

}

void AllDiffConstraint::organizeConfig() {
    if (mode == DECBI) decompose();
}

Cost AllDiffConstraint::evalOriginal( String s ) {
    Cost tuple_cost = 0;
    if (mode == DEC) {
        for (unsigned int i=0;i<s.length();i++) {
            for (unsigned int j=i+1;j<s.length();j++) {
                if (s[i] == s[j]) tuple_cost += def;
            }
        }
    } else {
        set<char> count;
        for (unsigned int i=0;i<s.length();i++) {
            count.insert(s[i]);
        }
        tuple_cost = (s.length() - count.size())*def;
    }
    return tuple_cost;
}

size_t AllDiffConstraint::GetGraphAllocatedSize() {
    return mapval.size() + arity_  + 2;
}

void AllDiffConstraint::buildGraph(Graph &g) {

    // if (g.size() == 0) g.setSize(mapval.size() + arity_  + 2);
    // g.clearEdge();
    for (int i=0;i<arity_;i++) {
        g.addEdge(0, i+1, 0);
        EnumeratedVariable* x = (EnumeratedVariable*)getVar(i);
        for (EnumeratedVariable::iterator j = x->begin(); j != x->end(); ++j) {
            int index = mapval[*j];
            if (index != 0) {
                g.addEdge(i+1, index, -deltaCost[i][x->toIndex(*j)]);
                vector<Cost> weight = g.getWeight(index, g.size()-1);
                Cost count = 0;
                if (weight.size() != 0) {
                    if (mode == DEC) {
                        count = *max_element(weight.begin(), weight.end())+def;
                    } else {
                        count = def;
                    }
                }
                g.addEdge(index, g.size()-1, count);
                //g.print();
            }
        }
    }

}

/*void AllDiffConstraint::getDomainFromGraph(Graph &graph, int varindex, vector<int> &domain) {

	domain.clear();
	for (vector<List_Node >::iterator k = graph[varindex+1].begin(); 
			k != graph[varindex+1].end(); k++) {
		if (k->adj > 0) {
			for (map<Value, Cost>::iterator i = mapval.begin();i !=
					mapval.end();i++) {
				if (i->second == k->adj) domain.push_back(i->first);
			}
		}
	}
	for (map<Value, Cost>::iterator i = mapval.begin();i !=
			mapval.end();i++) {
		for (vector<List_Node >::iterator k = graph[i->second].begin(); 
				k != graph[i->second].end(); k++) {
			if (k->adj == varindex+1) {
				domain.push_back(i->first);
			}
		}
	}

}*/

void AllDiffConstraint::decompose() {
    deconnect();
    for (int i=0;i<arity_;i++) {
        for (int j=i+1;j<arity_;j++) {
            EnumeratedVariable* x = (EnumeratedVariable*)getVar(i);
            EnumeratedVariable* y = (EnumeratedVariable*)getVar(j);
            vector<Cost> costs;
            for (unsigned int a = 0; a < x->getDomainInitSize(); a++) {
                for (unsigned int b = 0; b < y->getDomainInitSize(); b++) {
                    if (a == b) {
                        costs.push_back(def);
                    } else {
                        costs.push_back(0);
                    }
                }
            }
            if(ToulBar2::vac) {
                for (unsigned int a = 0; a < x->getDomainInitSize(); a++) {
                    for (unsigned int b = 0; b < y->getDomainInitSize(); b++) {
                        Cost c = costs[a * y->getDomainInitSize() + b];
                        wcsp->histogram(c);
                    }
                }
            }
            BinaryConstraint* ctr = x->getConstr(y);
            if(ctr)	{
                ctr->reconnect();
                ctr->addCosts(x,y,costs);
                ctr->propagate();
            }
            else {
                if (!ToulBar2::vac) {
                    ctr = new BinaryConstraint(wcsp, x, y, costs, &wcsp->getStore()->storeCost);
                } else {
                    ctr = new VACBinaryConstraint(wcsp, x, y, costs, &wcsp->getStore()->storeCost);
                }
            }
        }
    }
}

void AllDiffConstraint::dump(ostream& os, bool original)
{
    if (original) {
        os << arity_;
        for(int i = 0; i < arity_;i++) os << " " << scope[i]->wcspIndex;
    } else {
        os << nonassigned;
        for(int i = 0; i < arity_; i++) if (scope[i]->unassigned()) os << " " << scope[i]->getCurrentVarId();
    }
    os << " -1 salldiff" << endl << ((mode==VAR)?"var":"dec") << " " << def << endl;
}

void AllDiffConstraint::print(ostream& os) {
    os << "salldiff(";
    for(int i = 0; i < arity_;i++) {
        os << scope[i]->wcspIndex;
        if(i < arity_-1) os << ",";
    }
    os << ")[" << ((mode==VAR)?"var":"dec") << "," << def << "]";
}

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

