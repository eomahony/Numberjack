#include "tb2globalconstr.hpp"
#include "tb2wcsp.hpp"
#include "tb2vac.hpp"
#include "tb2clusters.hpp"

GlobalConstraint::GlobalConstraint(WCSP *wcsp, EnumeratedVariable** scope_in, int arity_in, Cost defval) : AbstractGlobalConstraint(wcsp, scope_in, arity_in),
        extendedCost(NULL), projectedCost(0),
        nonassigned(arity_in),
        currentVar(-1), needPropagateAC(false), needPropagateDAC(false), needPropagateEAC(false),
        currentDepth(-1), def(1), mode(0),
        count_nic(0), count_gac(0), count_fdac(0), count_edac(0), error(0)
{
    deltaCost = new vector<StoreCost>[arity_];
    for (int i=0;i<arity_in;i++) {
        deltaCost[i] = vector<StoreCost>(scope_in[i]->getDomainInitSize(), StoreCost(0));
        fill(deltaCost[i].begin(), deltaCost[i].end(), 0);
    }

    fullySupportedSet = new set<int>[arity_in];

    preUnaryCosts.resize(arity_in);
    for(int i = 0; i < arity_in; i++)
        preUnaryCosts[i].resize(scope[i]->getDomainInitSize(), 0);
}

GlobalConstraint::~GlobalConstraint()
{
    if (deltaCost != NULL) delete[] deltaCost;
    if (fullySupportedSet != NULL) delete[] fullySupportedSet;
}

void GlobalConstraint::init() {
    if (deconnected()) return;
    needPropagateAC = true;
    needPropagateDAC = true;
    needPropagateEAC = false;
    EACCost.clear();
    for (int i=0;i<arity_;i++) fullySupportedSet[i].clear();
    /*if ((currentDepth == -1) || (currentDepth >= Store::getDepth())) {
		initStructure();
	} else {
		vector<int> rmv;
		checkRemoved(rmv);
	}
	currentDepth = Store::getDepth();*/
    initStructure();
    propagate();
}

void GlobalConstraint::print(ostream& os) {
    os << this << " " << getName() << "(";
    int unassigned_ = 0;
    for(int i = 0; i < arity();i++) {
        if(scope[i]->unassigned()) unassigned_++;
        os << scope[i]->wcspIndex;
        if(i < arity_-1) os << ",";
    }
    // string smode = "";
    // for (map<string,int>::iterator iter = modeEnum.begin(); iter != modeEnum.begin(); ++iter) {
    // 	if ((*iter).second == mode) {
    // 	  smode = (*iter).first;
    // 	  break;
    // 	}
    // }
    os << ") mode: " << mode << ",defcost: " << def << ",";
    os << " unassigned: " << unassigned_ << "/" << nonassigned << "/" << arity();
    os << endl;
    if (ToulBar2::verbose >= 4) showCostProvidingPartition();
    if (ToulBar2::verbose >= 5) {
        os << "projectedCost: " << projectedCost << endl;
        for(int i = 0; i < arity();i++) {
            os << "deltaCost[" << i << "]:";
            for(unsigned int j = 0; j < deltaCost[i].size();j++) {os << " " << deltaCost[i][j];}
            os << endl;
            // os << "extendedCost[" << i << "]:";
            // for(unsigned int j = 0; j < extendedCost[i].size();j++) {os << " " << extendedCost[i][j];}
            // os << endl;
        }
    }
    if (ToulBar2::verbose >= 6) for(int i = 0; i < arity();i++) os << *scope[i] << endl;
}

Cost GlobalConstraint::eval(const String& s) {

    Cost tcost = evalOriginal(s);
    for (unsigned int i=0;i<s.length();i++) {
        EnumeratedVariable* x = (EnumeratedVariable*)getVar(i);
        if (tcost < wcsp->getUb()) {
            tcost -= deltaCost[i][x->toIndex(s[i]-CHAR_FIRST)];
        }
    }
    if (tcost < wcsp->getUb()) {
        tcost -= projectedCost;
    }
    assert(tcost >= 0);
    return tcost;

}

void GlobalConstraint::assign(int varIndex) {

    if (connected(varIndex)) {
        deconnect(varIndex);
        nonassigned = nonassigned - 1;
        if(nonassigned == 0) {
            deconnect();
            String t(arity_, CHAR_FIRST);
            for(int i=0;i<arity_;i++) {
                t[i] = getVar(i)->getValue() + CHAR_FIRST;
            }
            wcsp->revise(this);
            projectLB(eval(t));
        } else {
            pushAll();
        }
    }
}

void GlobalConstraint::project(int index, Value value, Cost cost, bool delayed) {
    if (deconnected()) return;
    assert(ToulBar2::verbose < 4 || ((cout << "[" << Store::getDepth() << "] project(" << getName() << ", " << getVar(index)->getName() << ", " << value << ", " << cost << ")" << endl), true));
    EnumeratedVariable* x = (EnumeratedVariable*)getVar(index);
    // hard binary constraint costs are not changed
    if (!CUT(cost + wcsp->getLb(), wcsp->getUb())) {
        TreeDecomposition* td = wcsp->getTreeDec();
        if(td) td->addDelta(cluster,x,value,cost);
        deltaCost[index][x->toIndex(value)] += cost;  // Warning! Possible overflow???
    }
    x->project(value, cost, delayed);
}

void GlobalConstraint::extend(int index, Value value, Cost cost) {
    assert(ToulBar2::verbose < 4 || ((cout << "[" << Store::getDepth() << "] extend(" << getName() << "," << getVar(index)->getName() << "," << value << ", " << cost << ")" << endl), true));
    EnumeratedVariable* x = (EnumeratedVariable*)getVar(index);
    TreeDecomposition* td = wcsp->getTreeDec();
    if(td) td->addDelta(cluster,x,value,-cost);
    deltaCost[index][x->toIndex(value)] -= cost;  // Warning! Possible overflow???
    x->extend(value, cost);
}

// function used for propagation
void GlobalConstraint::remove(int index) {
    //vector<int> rmv;
    currentVar = -1;
    needPropagateDAC = true;
    needPropagateAC = true;
}

void GlobalConstraint::projectFromZero(int index) {
    //vector<int> rmv;
    currentVar = -1;
    needPropagateDAC = true;
}

void GlobalConstraint::propagate() {

    for(int i=0;connected() && i<arity_;i++) {
        if (getVar(i)->assigned()) assign(i);
    }

    if (deconnected()) return;

    needPropagateEAC = false;
    currentVar = -1;

    wcsp->revise(this);

    vector<int> rmv;
    for (int i=0;i<arity_;i++) rmv.push_back(i);
    checkRemoved(rmv);

    switch (ToulBar2::LcLevel) {
    case LC_DAC:
        if (needPropagateDAC || needPropagateAC)
        {
            propagateDAC();
            needPropagateDAC = false;
            needPropagateAC = false;
        }
        break;
    case LC_EDAC:
        propagateEAC(); // not implemented yet (uses isEAC and findFullSupportEAC instead)
    case LC_FDAC:
        if (needPropagateDAC) {
            propagateDAC();
            needPropagateDAC = false;
        }
    case LC_AC:
        if (needPropagateAC) {
            propagateAC();
            needPropagateAC = false;
        }
        break;
    case LC_SNIC:
        propagateStrongNIC();
        break;
    default:
        break;
    }
}

void GlobalConstraint::propagateEAC() {
    //not done here but wait for variable's propagateEAC

    //	vector<int> provide;
    //	provide.resize(arity_);
    //	vector<map<Value, Cost> > deltas(arity_);
    //	for (int i=0;i<arity_;i++) {
    //		EnumeratedVariable *x = (EnumeratedVariable*)getVar(i);
    //		for (EnumeratedVariable::iterator v = x->begin();v != x->end();++v) {
    //			provide[i] = i;
    //			deltas[i][*v] = x->getCost(*v);
    //		}
    //	}
    //	changeAfterExtend(provide, deltas);
    //
    //	for (int i=0;i<arity_;i++) {
    //		if (getVar(i)->unassigned()) {
    //			map<Value, Cost> delta;
    //			EnumeratedVariable *x = (EnumeratedVariable*)getVar(i);
    //			findProjection(i, delta);
    //			for (map<Value, Cost>::iterator p = delta.begin(); p != delta.end();p++) {
    //				if (wcsp->getLb() + p->second >= wcsp->getUb()) {
    //					x->remove(p->first);
    //				}
    //			}
    //		}
    //	}
    //	//undoExtend();
    //    changeAfterProject(provide, deltas);
}

void GlobalConstraint::propagateDAC() {
    if (ToulBar2::verbose >= 3) cout << "propagateDAC for " << *this << endl;

    vector<map<Value, Cost> > deltas;
    vector<int> vars;
    vector<int> rmv;
    //checkRemoved(rmv);
    for(int ii = 0; ii < arity_; ii++){
        EnumeratedVariable * x = scope_dac[ii];
        int i = scope_inv[x->wcspIndex];
        if (x->unassigned()) {
            map<Value, Cost> delta;
            for(EnumeratedVariable::iterator it = x->begin(); it != x->end(); ++it){
                delta[*it] = x->getCost(*it);
                //deltaCost[i][x->toIndex(*it)] -= x->getCost(*it);
                preUnaryCosts[i][x->toIndex(*it)] = x->getCost(*it);
            }
            vars.push_back(i);
            deltas.push_back(delta);
        }
    }
    changeAfterExtend(vars, deltas);

    for(int ii = 0; ii < arity_; ii++){
        EnumeratedVariable * x = scope_dac[ii];
        int varindex = scope_inv[x->wcspIndex];
        if (x->unassigned()) {
            //delta.clear();
            map<Value, Cost> delta;
            findProjection(varindex, delta);
            for (map<Value, Cost>::iterator i = delta.begin(); i != delta.end();i++) {
                //deltaCost[varindex][x->toIndex(i->first)] += preUnaryCosts[varindex][x->toIndex(i->first)];
                Cost costDelta = i->second - preUnaryCosts[varindex][x->toIndex(i->first)];
                if(costDelta > 0){
                    if (CUT(costDelta + wcsp->getLb(), wcsp->getUb()))
                    {
                        i->second = preUnaryCosts[varindex][x->toIndex(i->first)];
                    }
                    project(varindex, i->first, costDelta, true);
                } else if (costDelta < 0) {
                    extend(varindex, i->first, -costDelta);
                }
            }
            changeAfterProject(varindex, delta);
            x->findSupport();

            rmv.clear();
            rmv.push_back(varindex);
            checkRemoved(rmv);
        }
    }
    //vector<int> rmv; //
    /*for (int i=0;i<arity_;i++) {
		if (getVar(i)->unassigned()) {
			checkRemoved(rmv);
			findFullSupport(i);
		}
	}*/
}

void GlobalConstraint::propagateAC() {
    vector<int> rmv;
    for (int i=0;i<arity_;i++) {
        if (getVar(i)->unassigned()) {
            findSupport(i);
            rmv.clear();
            rmv.push_back(i);
            checkRemoved(rmv);
        }
    }
}

void GlobalConstraint::propagateStrongNIC() {
    vector<int> rmv;
    checkRemoved(rmv);
    bool cont = false;
    do {
        cont = false;
        for (int varindex=0;varindex < arity_;varindex++) {
            EnumeratedVariable* x = (EnumeratedVariable*)getVar(varindex);
            if (x->unassigned()) {
                unsigned int oldD = x->getDomainSize();
                checkMinCost(varindex);
                if (oldD != x->getDomainSize()) cont = true;
            }
        }
        if (cont) checkRemoved(rmv);
    } while (cont);

    propagateNIC();
}

void GlobalConstraint::propagateNIC() {

    if (deconnected()) return;
    //wcsp->revise(this);
    vector<int> rmv;
    checkRemoved(rmv);
    Cost mincost = getMinCost();
    if (mincost-projectedCost > 0) {
        Cost diff = mincost - projectedCost;
        projectedCost += diff;
        wcsp->increaseLb(diff);
    }

}

/*void GlobalConstraint::getCostsWithUnary(int index, map<Value, Cost> &costs) {

	vector<int> rmv, support;
	checkRemoved(rmv);
	vector<map<Value, Cost> > deltas;
	for (set<int>::iterator i = fullySupportedSet[index].begin();i !=
			fullySupportedSet[index].end();i++) {
		EnumeratedVariable *x = (EnumeratedVariable*)getVar(*i);
		if (x->unassigned() && (*i != index)) {
			support.push_back(*i);
			map<Value, Cost> delta;
			for (EnumeratedVariable::iterator v = x->begin();v != x->end();++v) delta[*v] = x->getCost(*v);
			deltas.push_back(delta);
		}
	}
	changeAfterExtend(support, deltas);		

	map<Value, Cost> delta;
	EnumeratedVariable *x = (EnumeratedVariable*)getVar(index);
	findProjection(index, delta);
	for (map<Value, Cost>::iterator p = delta.begin(); p != delta.end();p++) {
		costs[p->first] += p->second;
	}
	undoExtend();

}*/

bool GlobalConstraint::isEAC(int index, Value a) {
    if (currentVar != index)
    {
        currentVar = index;

        vector<int> rmv;
        for (int i=0;i<arity_;i++) rmv.push_back(i);
        checkRemoved(rmv);

        vector<int> support;
        vector<map<Value, Cost> > deltas;
        for (set<int>::iterator i = fullySupportedSet[index].begin();i !=
                fullySupportedSet[index].end();i++) {
            EnumeratedVariable *x = (EnumeratedVariable*)getVar(*i);
            if (x->unassigned() && (*i != index)) {
                support.push_back(*i);
                map<Value, Cost> delta;
                for (EnumeratedVariable::iterator v = x->begin();v != x->end();++v) delta[*v] = x->getCost(*v);
                deltas.push_back(delta);
            }
        }
        changeAfterExtend(support, deltas);
        EACCost.clear();
        findProjection(index, EACCost);
        //undoExtend();
        changeAfterProject(support, deltas);
    }
    if (EACCost[a] < 0) {
        cout << "EAC error\n";
        EACCost[a] = 0;
    }
    //if (EACCost[a] > 0) needPropagateEAC = true;
    return (EACCost[a] == 0);

}

void GlobalConstraint::fillEAC2(int varindex) {
    currentVar = -1;
    for (int index=0;index<arity_;index++) {
        if (index != varindex) {
            EnumeratedVariable* var = (EnumeratedVariable*)getVar(index);
            if (var->unassigned()) var->queueEAC2();
        }
    }
}

//TODO: applies DAC order when enumerating variables (fullySupportedSet does not preserve DAC order)
void GlobalConstraint::findFullSupportEAC(int index) {
    if (ToulBar2::verbose >= 3) cout << "findFullSupportEAC for variable " << index << endl;
    currentVar = -1;
    //if (needPropagateEAC)
    {
        //needPropagateEAC = false;

        vector<map<Value, Cost> > deltas;
        vector<int> vars;
        vector<int> rmv;
        //vars.push_back(index);
        vars.push_back(index);
        for (set<int>::iterator i = fullySupportedSet[index].begin();i !=
                fullySupportedSet[index].end();i++) {
            EnumeratedVariable * x = scope[*i];
            if (x->unassigned() && (*i != index)) {
                vars.push_back(*i);
            }
        }

        for(vector<int>::iterator i = vars.begin(); i != vars.end(); ++i){
            EnumeratedVariable* x = scope[*i];
            map<Value, Cost> delta;
            for(EnumeratedVariable::iterator it = x->begin(); it != x->end(); ++it){
                delta[*it] = x->getCost(*it);
                //deltaCost[i][x->toIndex(*it)] -= x->getCost(*it);
                preUnaryCosts[*i][x->toIndex(*it)] = x->getCost(*it);
            }
            //vars.push_back(*i);
            deltas.push_back(delta);
        }
        changeAfterExtend(vars, deltas);

        for(vector<int>::iterator it = vars.begin(); it != vars.end(); ++it){
            int varindex = *it;
            EnumeratedVariable * x = scope[varindex];
            if (x->unassigned()) {
                //delta.clear();
                map<Value, Cost> delta;
                findProjection(varindex, delta);
                for (map<Value, Cost>::iterator i = delta.begin(); i != delta.end();i++) {
                    //deltaCost[varindex][x->toIndex(i->first)] += preUnaryCosts[varindex][x->toIndex(i->first)];
                    Cost costDelta = i->second - preUnaryCosts[varindex][x->toIndex(i->first)];
                    if(costDelta > 0){
                        if (CUT(costDelta + wcsp->getLb(), wcsp->getUb()))
                        {
                            i->second = preUnaryCosts[varindex][x->toIndex(i->first)];
                        }
                        project(varindex, i->first, costDelta, true);
                    } else if (costDelta < 0) {
                        extend(varindex, i->first, -costDelta);
                    }
                }
                changeAfterProject(varindex, delta);
                x->findSupport();

                rmv.clear();
                rmv.push_back(varindex);
                checkRemoved(rmv);
            }
        }

    }
}

// THIS CODE IS NEVER USED!!!
//TODO: applies DAC order when enumerating variables
void GlobalConstraint::findFullSupport(int varindex, vector<int> &support, bool isEAC)
{
    if (ToulBar2::verbose >= 3) cout << "findFullSupport for variable " << varindex << endl;
    //wcsp->revise(this);
    EnumeratedVariable* var = (EnumeratedVariable*)getVar(varindex);
    vector<map<Value, Cost> > deltas(support.size());
    int count = 0;
    for (vector<int>::iterator i = support.begin();i !=
            support.end();i++, count++) {
        EnumeratedVariable *x = (EnumeratedVariable*)getVar(*i);
        for (EnumeratedVariable::iterator v = x->begin();v != x->end();++v) {
            deltas[count][*v] = x->getCost(*v);
        }
    }
    changeAfterExtend(support, deltas);

    map<Value, Cost> delta;
    findProjection(varindex, delta);

    map<Value, Cost> P;
    bool allzero = true;
    for (map<Value, Cost>::iterator i = delta.begin(); i != delta.end();i++) {
        if (i->second > 0) {
            P[i->first] = i->second;
            allzero = false;
        } else if (i->second < 0) {
            P[i->first] = 0;
            cout << "projection error\n";
        }
    }

    if (allzero) {
        //undoExtend();
        changeAfterProject(support, deltas);
        return;
    } else {
        count_fdac++;
    }

    changeAfterProject(varindex, delta);

    //	map<Value,Cost> *E = new map<Value, Cost>[support.size()];
    vector< map<Value,Cost> > E;
    E.resize(support.size());

    int index = 0;
    for (vector<int>::iterator i = support.begin();i !=
            support.end();i++, index++) {
        EnumeratedVariable *x = (EnumeratedVariable*)getVar(*i);
        if (x->unassigned()) {
            delta.clear();
            findProjection(*i, delta);
            for (map<Value, Cost>::iterator j = delta.begin(); j != delta.end();j++) {
                E[index][j->first] = min(x->getCost(j->first) - j->second, x->getCost(j->first));
                /*if (j->second < 0) {
				  j->second = 0;
				  }*/
                //if (E[index][j->first] > 0) extendedCost[*i][j->first] += j->second;
            }
            changeAfterProject(*i, delta);
        }
    }

    index = 0;
    for (vector<int>::iterator i = support.begin();i !=
            support.end();i++, index++) {
        EnumeratedVariable *x = (EnumeratedVariable*)getVar(*i);
        if (x->unassigned()) {
            for (map<Value, Cost>::iterator j = E[index].begin(); j != E[index].end();j++) {
                if (j->second > 0) extend(*i, j->first, j->second);
            }
        }
    }
    index = 0;
    for (vector<int>::iterator i = support.begin();i !=
            support.end();i++, index++) {
        EnumeratedVariable *x = (EnumeratedVariable*)getVar(*i);
        if (x->unassigned()) {
            bool supportbroken = false;
            for (map<Value, Cost>::iterator j = E[index].begin(); j != E[index].end();j++) {
                if (j->second < 0) {
                    if (j->first == x->getSupport()) supportbroken = true;
                    project(*i, j->first, -j->second);
                }
            }
            if (supportbroken) x->findSupport();
        }
    }

    bool supportbroken = false;
    for (map<Value, Cost>::iterator i = P.begin(); i != P.end();i++) {
        if (i->second > 0) {
            if (i->first == var->getSupport()) supportbroken = true;
            project(varindex, i->first, i->second);
        }
    }
    if (supportbroken) var->findSupport();
    //	delete[] E;

}

void GlobalConstraint::findSupport(int varindex)
{	
    //wcsp->revise(this);
    if (ToulBar2::verbose >= 3) cout << "findSupport for variable " << varindex << endl;
    map<Value, Cost> delta;
    findProjection(varindex, delta);
    bool allzero = true;
    for (map<Value, Cost>::iterator i = delta.begin(); i != delta.end();i++) {
        if (i->second > 0) allzero = false;
    }
    if (!allzero) {
        count_gac++;
        bool supportbroken = false;
        EnumeratedVariable* x = (EnumeratedVariable*)getVar(varindex);
        x->queueAC();
        for (map<Value, Cost>::iterator i = delta.begin(); i != delta.end();i++) {
            if (i->second > 0) {
                if (i->first == x->getSupport()) supportbroken = true;
                project(varindex, i->first, i->second);
            }
        }
        if (supportbroken) {
            x->findSupport();
        }
        changeAfterProject(varindex, delta);
    }

}

void GlobalConstraint::checkMinCost(int varindex) {

    count_nic++;

    vector<Value> removed;
    EnumeratedVariable* x = (EnumeratedVariable*)getVar(varindex);

    map<Value, Cost> delta;
    findProjection(varindex, delta);

    for (EnumeratedVariable::iterator j = x->begin(); j != x->end(); ++j) {
        if (wcsp->getLb() + delta[*j] + x->getCost(*j) - projectedCost >= wcsp->getUb()) {
            removed.push_back(*j);
        }
    }
    for (vector<Value>::iterator i = removed.begin(); i != removed.end();i++) {
        x->remove(*i);
    }

    if (!removed.empty()) x->queueAC();

}

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

