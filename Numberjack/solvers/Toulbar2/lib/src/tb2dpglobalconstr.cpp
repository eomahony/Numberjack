#include "tb2dpglobalconstr.hpp"
#include <vector>
#include <algorithm>
using namespace std;

DPGlobalConstraint::DPGlobalConstraint(WCSP * wcsp, EnumeratedVariable ** scope, int arity) : GlobalConstraint(wcsp, scope, arity, 0), initialized(false){
    zero = new vector<bool>[arity];
    for(int i = 0; i < arity; i++)
        zero[i] = vector<bool>(scope[i]->getDomainInitSize(), false);
    preUnaryCosts = new vector<Cost>[arity];
    for(int i = 0; i < arity; i++)
        preUnaryCosts[i] = vector<Cost>(scope[i]->getDomainInitSize(), 0);
}

DPGlobalConstraint::~DPGlobalConstraint(){
    delete[] zero;
    delete[] preUnaryCosts;
}

void DPGlobalConstraint::clear(){
    for(int i = 0; i < arity(); i++) {
        fill(zero[i].begin(), zero[i].end(), false);            
        fill(preUnaryCosts[i].begin(), preUnaryCosts[i].end(), 0);            
    }
}

void DPGlobalConstraint::record(Value *tuple){
    if(tuple == NULL) return;
    for(int i = 0; i < arity(); i++)
        zero[i][scope[i]->toIndex(tuple[i])] = true;
    if(ToulBar2::verbose >= 3){
        cout << "tuple(";
        for(int i = 0 ; i < arity(); i++) cout << tuple[i] << ",";
        cout << ")" << endl;
    }
    delete[] tuple;
}

void DPGlobalConstraint::propagateNIC(){

    Cost least = minCostOriginal();
    if(least > projectedCost){
        wcsp->increaseLb(least - projectedCost);
        projectedCost = least;        	
    }
}

void DPGlobalConstraint::propagateStrongNIC(){

    propagateNIC();
    Cost ub = wcsp->getUb();
    Cost lb = wcsp->getLb();
    bool changed = true;
    for(int i = 0; i < arity_; i++){
        EnumeratedVariable * x = scope[i];
        if(x->assigned()) continue;
        bool first = true;
        for(EnumeratedVariable::iterator it = x->begin(); it != x->end(); ++it){
            Cost cost = minCostOriginal(i, *it, changed && first) - projectedCost;
            changed = first = false;
            if(cost + x->getCost(*it) + lb >= ub) {           
                x->remove(*it);
                changed = true;
            }
        }
        x->findSupport();
    }
}

void DPGlobalConstraint::propagateAC(){	

    bool changed = true;
    clear();  

    //Cost thisUb;
    //thisUb = wcsp->getUb();
    for(int i = 0; i < arity(); i++){
        EnumeratedVariable * x = scope[i];
        bool first = true;
        for(EnumeratedVariable::iterator it = x->begin(); it != x->end(); ++it){            
            if(zero[i][x->toIndex(*it)]) continue;
            Result r = minCost(i, *it, changed && first);            
            if(changed && first) changed = false;
            first = false;

            Cost cost = r.first;
            /*if(cost >= thisUb || cost + wcsp->getLb() + scope[i]->getCost(*it) >= thisUb){            
			  x->remove(*it);
			  changed = true;
            }else */
            if(cost > 0){
                project(i, *it, cost);
                changed = true;
            } else if(cost < 0) {
                /* Should not happen*/
                printf("Warning: AC propagation get negative cost\n");
                extend(i, *it, -cost);
                changed = true;
            }
            if(x->canbe(*it)) record(r.second);
        }
        x->findSupport();
    }
}

void DPGlobalConstraint::findSupport(int var, bool &changed){
    EnumeratedVariable * x = scope[var];
    bool first = true;
    vector<Value> remove;

    for(EnumeratedVariable::iterator it = x->begin(); it != x->end(); ++it){
        Cost cost;
        Result r = pair<Cost, Value*>(0, NULL);       
        if(zero[var][x->toIndex(*it)]){
            cost = 0;
        }else{
            r = minCost(var, *it, changed && first);
            if(changed && first) changed = false;
            first = false;
            cost = r.first;
        }      
        //deltaCost[var][x->toIndex(*it)] += x->getCost(*it);
        deltaCost[var][x->toIndex(*it)] += 	preUnaryCosts[var][x->toIndex(*it)];
        //Cost delta = cost - x->getCost(*it);
        Cost delta = cost - preUnaryCosts[var][x->toIndex(*it)];
        if(delta > 0){
            project(var, *it, delta, true); // NC will be delayed (avoid forward checking on binary/ternay cost functions)
            assert(x->canbe(*it));
        } else if (delta < 0) {
            extend(var, *it, -delta);
        }
        if (!zero[var][x->toIndex(*it)] && x->getCost(*it) + wcsp->getLb() < wcsp->getUb()) record(r.second);
    }
    x->findSupport();
    changed = true;  //Detect any change in variable domain or unary costs
}

void DPGlobalConstraint::propagateDAC(){
    if (ToulBar2::verbose >= 3) cout << "propagateDAC for " << *this << endl;

    clear();

    for(int ii = 0; ii < arity_; ii++){
        EnumeratedVariable * x = scope_dac[ii];
        int i = scope_inv[x->wcspIndex];
        for(EnumeratedVariable::iterator it = x->begin(); it != x->end(); ++it){
            if (x->unassigned()) {
                deltaCost[i][x->toIndex(*it)] -= x->getCost(*it);
                preUnaryCosts[i][x->toIndex(*it)] = x->getCost(*it);
            }
        }
    }

    bool changed = true;    
    for(int ii = 0; ii < arity_; ii++){
        EnumeratedVariable * x = scope_dac[ii];
        int i = scope_inv[x->wcspIndex];
        if (x->unassigned()) {
            findSupport(i, changed);
        }
    }

}

bool DPGlobalConstraint::isEAC(int var, Value val){		

    for(set<int>::iterator it = fullySupportedSet[var].begin(); it !=
            fullySupportedSet[var].end(); ++it){
        EnumeratedVariable *x = scope[*it];
        if (x->unassigned() && (*it != var)) {
            for(EnumeratedVariable::iterator jt = x->begin(); jt != x->end(); ++jt){
                deltaCost[*it][x->toIndex(*jt)] -= x->getCost(*jt);
            }
        }
    }
    bool ret = (minCost(var, val, true).first == 0);
    for(set<int>::iterator it = fullySupportedSet[var].begin(); it !=
            fullySupportedSet[var].end(); ++it){
        EnumeratedVariable *x = scope[*it];
        if (x->unassigned() && (*it != var)) {
            for(EnumeratedVariable::iterator jt = x->begin(); jt != x->end(); ++jt){
                deltaCost[*it][x->toIndex(*jt)] += x->getCost(*jt);
            }
        }
    }
    return ret;
}

//TODO: applies DAC order when enumerating variables (fullySupportedSet does not preserve DAC order)
void DPGlobalConstraint::findFullSupportEAC(int var){
    assert(fullySupportedSet[var].find(var) == fullySupportedSet[var].end());

    clear();
    //fullySupportedSet[var].insert(var);
    for(set<int>::iterator it = fullySupportedSet[var].begin(); it !=
            fullySupportedSet[var].end(); ++it){
        EnumeratedVariable *x = scope[*it];
        if (x->unassigned() && (*it != var)) {
            for(EnumeratedVariable::iterator jt = x->begin(); jt != x->end(); ++jt){
                /* fix the problem in EAC */
                preUnaryCosts[*it][x->toIndex(*jt)] = x->getCost(*jt);
                deltaCost[*it][x->toIndex(*jt)] -= x->getCost(*jt);
            }
        }
    }
    //fullySupportedSet[var].erase(var);
    EnumeratedVariable *cur = scope[var];
    for(EnumeratedVariable::iterator jt = cur->begin(); jt != cur->end(); ++jt){
        /* fix the problem in EAC */
        preUnaryCosts[var][cur->toIndex(*jt)] = cur->getCost(*jt);
        deltaCost[var][cur->toIndex(*jt)] -= cur->getCost(*jt);
    }


    bool changed = true;
    findSupport(var, changed);
    for(set<int>::iterator it = fullySupportedSet[var].begin(); it !=
            fullySupportedSet[var].end(); ++it)
    {
        EnumeratedVariable *x = scope[*it];
        if (x->unassigned() && (*it != var)) {
            findSupport(*it, changed);
        }
    }
}

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

