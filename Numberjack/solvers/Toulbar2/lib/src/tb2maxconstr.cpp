#include "tb2maxconstr.hpp"
#include "tb2enumvar.hpp"
#include <vector>
#include <fstream>
#include <algorithm>
#include <set>

using namespace std;

MaxConstraint::MaxConstraint(WCSP * wcsp, EnumeratedVariable ** scope, int arity):
            DPGlobalConstraint(wcsp, scope, arity), top(MIN_COST), largest(MIN_COST) {
    weightMap.resize(arity);
}

MaxConstraint::~MaxConstraint(){
}

void MaxConstraint::read(istream &file){
    //    int n = arity();
    // weightMap.resize(n);

    file >> def;
    /*for(int i = 0; i < n; i++){
		EnumeratedVariable * x = scope[i];        
        for(EnumeratedVariable::iterator it = x->begin(); it != x->end(); ++it)
            weightMap[i][*it] = def; 
    }*/

    int nTuple;
    file >> nTuple;
    top = def;
    for (int it=0;it<nTuple;it++)
    {
        int varID;
        unsigned int v;
        Cost w;
        file >> varID >> v >> w;
        setAssignmentWeight((EnumeratedVariable*)(wcsp->getVar(varID)), v, w);
    }

}

void MaxConstraint::initMemoization() {
    int n = arity();
    for(int i = 0; i < n; i++){
        EnumeratedVariable * x = scope[i];
        for(EnumeratedVariable::iterator it = x->begin(); it != x->end(); ++it) {
            map<Value, Cost>::iterator pos = weightMap[i].find(*it);
            if (pos == weightMap[i].end()) weightMap[i][*it] = def; 
        }
    }    
    top = max(MAX_COST, wcsp->getUb());              
    mincosts.resize(n);       
    stack.resize(n);
    last.resize(n);
    best.resize(n);
}

Cost MaxConstraint::evalOriginal( String s ) {
    int largeComp = 0;
    int n = arity();
    for(int i = 0; i < n; i++){
        if (largeComp < weightMap[i][s[i] - CHAR_FIRST]) largeComp = weightMap[i][s[i] - CHAR_FIRST];    		
    }
    return largeComp;
}

Cost MaxConstraint::minCostOriginal(){		
    findLargest();
    return largest;	
}

Cost MaxConstraint::minCostOriginal(int var, Value val, bool changed){	
    if(changed) findLargest();
    return max(weightMap[var][val], largest);	
}

void MaxConstraint::findLargest(){
    largest = 0;
    for(int i = 0; i <  arity(); i++){
        EnumeratedVariable * x = scope[i];
        Cost tmp = top;
        for(EnumeratedVariable::iterator it = x->begin(); it != x->end(); ++it)           			
            if (tmp > weightMap[i][*it]) tmp = weightMap[i][*it];
        if (largest < tmp) largest = tmp;
    }
}

void MaxConstraint::recompute() {

    int n = arity();
    sorted.clear();
    for(int i = 0; i < n; i++){
        EnumeratedVariable *x = scope[i];
        for(EnumeratedVariable::iterator it = x->begin(); it != x->end(); ++it)
            sorted.push_back(Entry(i, *it, weightMap[i][*it]));
    }

    sort(sorted.begin(), sorted.end());
    int m = sorted.size();

    int cnt = 0;
    Cost sum = 0;
    query.resize(m);
    cost.resize(m);    
    for(int i = 0; i < n; i++){
        last[i] = -1;
        stack[i].clear();
    }

    for(int i = 0; i < m; i++){
        int var = sorted[i].var;
        if(last[var] < 0){
            last[var] = i;
            stack[var].push_back(i);
            best[var] = unary(sorted[i].var, sorted[i].val);
            sum += best[var];
            cnt++;
        }else{
            query[last[var]] = i-1;
            last[var] = i;
            if(best[var] > unary(sorted[i].var, sorted[i].val)){
                sum -= best[var];
                best[var] = unary(sorted[i].var, sorted[i].val);
                sum += best[var];
                stack[var].push_back(i);
            }
        }
        if(cnt < n)
            cost[i] = top;
        else
            cost[i] = sum + sorted[i].weight - best[var] + unary(sorted[i].var, sorted[i].val);
    }

    link.resize(m);
    for(int i = 0; i < m ; i++)
        link[i] = i;
    for(int i = 0; i < n; i++){
        query[last[i]] = m - 1;
        best[i] = top;
    }

    tree.resize(0);
    for(int i = m-1; i >= 0; i--){
        int var = sorted[i].var;
        int val = sorted[i].val;
        if(query[i] > i){
            Cost cur = cost[ancestor(query[i])];
            int curVar = sorted[stack[var].back()].var;
            int curVal = sorted[stack[var].back()].val;
            cur -= unary(curVar, curVal);
            best[var] = min(best[var], cur);
        }
        if(stack[var].back() == i) stack[var].pop_back();
        mincosts[var][val] = min(cost[i], best[var] + unary(sorted[i].var, sorted[i].val));
        while(tree.size() != 0 && cost[*tree.rbegin()] >= cost[i]){
            link[*tree.rbegin()] = i;
            tree.pop_back();
        }
        tree.push_back(i);
    }
}

DPGlobalConstraint::Result MaxConstraint::minCost(int var, Value val, bool changed){	
    if(changed) recompute();
    return DPGlobalConstraint::Result(mincosts[var][val], NULL);	
}


Cost MaxConstraint::unary(int var, int val){
    return -deltaCost[var][scope[var]->toIndex(val)];
}

int MaxConstraint::ancestor(int i){
    int ret = i;
    while(link[ret] != ret)
        ret = link[ret];
    for(int j = i; j != ret; j = i){
        i = link[i];
        link[j] = ret;
    }
    return ret;
}


/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

