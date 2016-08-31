#include "tb2regulardpconstr.hpp"
#include <vector>
#include <fstream>
#include <string>
using namespace std;

RegularDPConstraint::RegularDPConstraint(WCSP * wcsp, EnumeratedVariable ** scope, int arity)
: DPGlobalConstraint(wcsp, scope, arity)
, top(0) {
}


RegularDPConstraint::~RegularDPConstraint(void) {
	deleteTable(f);
	deleteTable(u);
	deleteTable(invf);
	deleteTable(curf);	
}

void RegularDPConstraint::read(istream & file) {
	string str;
	file >> str >> def;

	// input the automaton
	int nstate;
	file >> nstate;
	dfa.setNumStates(nstate);
	int nstart;
	file >> nstart;
	for (int i=0;i<nstart;i++) {
		int t; file >> t;
		dfa.init.push_back(t);		
	} 
	int nfinal;
	file >> nfinal;
	for (int i=0;i<nfinal;i++) {
		int t; file >> t;
		dfa.final.push_back(t);	
	} 
	int ntransition;
	file >> ntransition;
	for (int i=0;i<ntransition;i++) {
		int start, end, symbol;
		file >> start;
		file >> symbol;
		file >> end;		
		dfa.addTransition(start, symbol, end, 0);
	}
	
}

void RegularDPConstraint::initMemoization() {
    dfa.finalize();

    resizeTable(f, arity()+1, dfa.size());
	resizeTable(curf, arity()+1, dfa.size());
	resizeTable(invf, arity()+1, dfa.size());
	resizeTable(u, arity()+1, dfa.symbol.size());

	top = max(wcsp->getUb(), MAX_COST);
}

Cost RegularDPConstraint::minCostOriginal() {
	int n = arity();
	for (int i=1;i<=n;i++) {
		for (unsigned int j=0;j<dfa.symbol.size();j++) {									
			u[i][j].val = top;
			EnumeratedVariable *x = scope[i-1];				
			for(EnumeratedVariable::iterator it = x->begin(); it != x->end(); ++it) {
				if (u[i][j].val > unary(dfa.symbol[j], i-1, *it))
					u[i][j].val = unary(dfa.symbol[j], i-1, *it);				
			}						  
		}
	}

	recomputeTable(curf);

	int minCost = top;
	for (vector<int>::iterator s = dfa.final.begin();s != dfa.final.end();s++)
		if (minCost > curf[n][*s].val) minCost = curf[n][*s].val;

	return minCost;
}

Cost RegularDPConstraint::minCostOriginal(int var, Value val, bool changed) {	
	return minCost(var, val, changed).first;
}

Cost RegularDPConstraint::eval(const String& s) {
	int n = arity();
	for (int i=1;i<=n;i++) {
		for (unsigned int j=0;j<dfa.symbol.size();j++) {												
			u[i][j].val = unary(dfa.symbol[j], i-1, s[i-1]-CHAR_FIRST);			
		}
	}

	recomputeTable(curf);

	int minCost = top;
	for (vector<int>::iterator s = dfa.final.begin();s != dfa.final.end();s++)
		if (minCost > curf[n][*s].val) minCost = curf[n][*s].val;

	return minCost - projectedCost;
}

void RegularDPConstraint::recompute() {
	int n = arity();
	for (int i=1;i<=n;i++) {
		for (unsigned int j=0;j<dfa.symbol.size();j++) {									
			u[i][j].val = top;			
			u[i][j].source	= -1;
			EnumeratedVariable *x = scope[i-1];				
			for(EnumeratedVariable::iterator it = x->begin(); it != x->end(); ++it) {
				if (u[i][j].val > unary(dfa.symbol[j], i-1, *it)) {
					u[i][j].val = unary(dfa.symbol[j], i-1, *it);					
					u[i][j].source	= *it;
				}
			}						  
		}
	}
	recomputeTable(f, invf);
}

DPGlobalConstraint::Result RegularDPConstraint::minCost(int var, Value val, bool changed) {

  	if (changed) recompute();

	int minCost = wcsp->getUb();        			
	for (int qk=0;qk<dfa.size();qk++) {
		for (vector<pair<int, int> >::iterator qj = dfa.transition[qk].begin();qj != dfa.transition[qk].end();qj++) {             
			int curCost = f[var][qk].val + unary(qj->first, var, val)  + invf[var+1][qj->second].val;
			if (minCost > curCost) minCost = curCost;						
		}	
	}	                             

	return DPGlobalConstraint::Result(minCost, NULL);
}

void RegularDPConstraint::recomputeTable(DPTableCell** table, DPTableCell** invTable, int startRow) {	
	int n = arity();

	if (startRow == 0)
	{
		for (int j=0;j<dfa.size();j++) {
			table[0][j].val = top;				
			table[0][j].source = make_pair(-1, -1);		
		}	
		for (vector<int>::iterator it = dfa.init.begin();it != dfa.init.end();it++) {
			table[0][*it].val = 0;	
			table[0][*it].source = make_pair(-1, *it);
		}	
		startRow++;
	}

	for (int i=startRow;i<=n;i++) {
		for (int j=0;j<dfa.size();j++) {
			table[i][j].val = top;
			table[i][j].source = make_pair(-1,-1);				
			for (vector<pair<int, int> >::iterator qk = dfa.invTransition[j].begin();qk != dfa.invTransition[j].end();qk++) {                        
				int curCost = table[i-1][qk->second].val + u[i][dfa.symbolIndex[qk->first]].val;
				if (table[i][j].val > curCost) {
					table[i][j].val = curCost;
					table[i][j].source =  make_pair(u[i][dfa.symbolIndex[qk->first]].source, qk->second);
				}                                
			}

		}	
	}

	if (invTable != NULL) {
		for (int j=0;j<dfa.size();j++) invTable[n][j].val = top; 	        
		for (vector<int>::iterator it = dfa.final.begin();it != dfa.final.end();it++) invTable[n][*it].val = 0;							

		for (int i=n-1;i>=0;i--) {
			for (int j=0;j<dfa.size();j++) {
				invTable[i][j].val = top;				
				for (vector<pair<int, int> >::iterator qj = dfa.transition[j].begin();qj != dfa.transition[j].end();qj++) {									
					int curCost = invTable[i+1][qj->second].val + u[i+1][dfa.symbolIndex[qj->first]].val;
					if (invTable[i][j].val > curCost) {
						invTable[i][j].val = curCost;  
						invTable[i][j].source =  make_pair(u[i+1][dfa.symbolIndex[qj->first]].source, qj->second);
					}
				}						
			}	
		}
	}
}

Cost RegularDPConstraint::unary(int ch, int var, Value v) {
	Cost ucost = (v==ch)?0:def;	
        EnumeratedVariable *x = scope[var];	
	return ucost - deltaCost[var][x->toIndex(v)];
}

void RegularDPConstraint::dump(ostream& os, bool original)
{
    if (original) {
        os << arity_;
        for(int i = 0; i < arity_;i++) os << " " << scope[i]->wcspIndex;
    } else {
        os << nonassigned;
        for(int i = 0; i < arity_; i++) if (scope[i]->unassigned()) os << " " << scope[i]->getCurrentVarId();
    }
    os << " -1 sregulardp var " << def << endl;
    dfa.dump(os, original);
}

void RegularDPConstraint::print(ostream& os) {
    os << "sregulardp(";
    for (int i = 0; i < arity_; i++) {
        os << scope[i]->wcspIndex;
        if (i < arity_ - 1) os << ",";
    }
    os << ")[" << def << "]";
    dfa.dump(os, true);
}

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

