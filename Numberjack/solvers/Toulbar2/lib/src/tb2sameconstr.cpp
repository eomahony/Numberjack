#include "tb2sameconstr.hpp"
#include "tb2wcsp.hpp"

SameConstraint::SameConstraint(WCSP *wcsp, EnumeratedVariable** scope_in, 
        int arity_in) : FlowBasedGlobalConstraint(wcsp, scope_in, arity_in) {
    buildIndex();
}

void SameConstraint::buildIndex() {
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
    nDistinctDomainValues = D.size();
    //graph.setSize(arity_+D.size()+2);
}

pair<int,int> SameConstraint::mapto(int varindex, Value val) {
    if (binary_search(group[0].begin(), group[0].end(), varindex)) {
        return make_pair(varindex+1, mapval[val]);
        //return make_pair(0,0);
    } else if (binary_search(group[1].begin(), group[1].end(), varindex)) {
        return make_pair(mapval[val], varindex+1);
        //return make_pair(0,0);
    } else {
        cout << "no group belong ?" << endl;
        exit(0);
    }
}

void SameConstraint::read(istream &file) {
    file >> def;
    int size[2];
    //cout << "def. value = " << def << endl;
    //cout << "consistency level = " << ToulBar2::consistencyLevel << endl;
    file >> size[0];
    file >> size[1];
    for (int g=0;g<2;g++) {
        for (int i=0;i<size[g];i++) {
            int var;
            file >> var;
            for (int j=0;j<arity_;j++) {
                if (wcsp->getVar(var) == getVar(j)) {
                    group[g].push_back(j);
                    break;
                }
            }
        }
        sort(group[g].begin(), group[g].end());
    }
}

Cost SameConstraint::evalOriginal(const String& s) {
    Cost tuple_cost = 0;
    map<char, int> appear;
    for (vector<int>::iterator i = group[0].begin();i != group[0].end();i++) {
        appear[s[*i]] += def;
    }
    for (vector<int>::iterator i = group[1].begin();i != group[1].end();i++) {
        appear[s[*i]] -= def;
    }
    int sum = 0;
    for (map<char, int>::iterator i = appear.begin();i != appear.end();i++) {
        sum += (i->second<0)?(-(i->second)):i->second;
    }
    tuple_cost += sum/2;
    /*for (int i=0;i<s.length();i++) {
		if (tuple_cost < wcsp->getUb()) {
			tuple_cost -= deltaCost[i][s[i]-CHAR_FIRST];
		}
	}
	tuple_cost -= projectedCost;
	if (tuple_cost < 0) {
		cout << "Error ! " << s << " has -ve cost of " << tuple_cost << endl;
		exit(0);
	}*/
    return tuple_cost;
}

size_t SameConstraint::GetGraphAllocatedSize() {
    return arity_+nDistinctDomainValues+2;
}

void SameConstraint::buildGraph(Graph &g) {
    //g.clearEdge();
    for (vector<int>::iterator i = group[0].begin();i != group[0].end();i++) {
        EnumeratedVariable *x = (EnumeratedVariable*)getVar(*i);
        for (EnumeratedVariable::iterator v = x->begin();v != x->end();++v) {
            g.addEdge((*i)+1, mapval[*v], -deltaCost[*i][x->toIndex(*v)], 1, *v);
        }
        g.addEdge(0, (*i)+1, 0);
    }
    for (vector<int>::iterator i = group[1].begin();i != group[1].end();i++) {
        EnumeratedVariable *x = (EnumeratedVariable*)getVar(*i);
        for (EnumeratedVariable::iterator v = x->begin();v != x->end();++v) {
            g.addEdge(mapval[*v], (*i)+1, -deltaCost[*i][x->toIndex(*v)], 1, *v);
        }
        g.addEdge((*i)+1, g.size()-1, 0);
    }
    for (map<Value, int>::iterator i = mapval.begin(); i != mapval.end();i++) {
        map<Value, int>::iterator j = i;
        for (j++; j != mapval.end();j++) {
            g.addEdge(i->second, j->second, def, arity_);
            g.addEdge(j->second, i->second, def, arity_);
        }
    }
}

/*void SameConstraint::getDomainFromGraph(Graph &graph, int varindex, vector<int> &domain) {

	domain.clear();
	for (vector<List_Node >::iterator k = graph[varindex+1].begin(); 
			k != graph[varindex+1].end(); k++) {
		for (map<Value, int>::iterator i = mapval.begin();i !=
				mapval.end();i++) {
			if (i->second == k->adj) domain.push_back(i->first);
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


void SameConstraint::dump(ostream& os, bool original)
{
    if (original) {
        os << arity_;
        for(int i = 0; i < arity_;i++) os << " " << scope[i]->wcspIndex;
    } else {
        os << nonassigned;
        for(int i = 0; i < arity_; i++) if (scope[i]->unassigned()) os << " " << scope[i]->getCurrentVarId();
    }
    os << " -1 ssame " << def << " " <<  group[0].size() << " " << group[1].size() << endl;
    for (int g=0; g<2; g++) {
        for (unsigned int i = 0; i < group[g].size(); i++) {
            os << " " << getVar(group[g][i])->wcspIndex;
        }
    }
    os << endl;
}

//void SameConstraint::print(ostream& os)
//{
//    os << "ssame(";
//    for(int i = 0; i < arity_;i++) {
//        os << scope[i]->wcspIndex;
//        if(i < arity_-1) os << ",";
//    }
//    os << ")[" << def << "," << group[0].size() << "," << group[1].size();
//    for (int g=0; g<2; g++) {
//        for (unsigned int i = 0; i < group[g].size(); i++) {
//            os << "," << getVar(group[g][i])->wcspIndex;
//        }
//    }
//    os << "]";
//}

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

