#include "tb2globalconstr.hpp"
#include "tb2wcsp.hpp"
#include "tb2vac.hpp"
#include "tb2clusters.hpp"


GlobalConstraint::GlobalConstraint(WCSP *wcsp, EnumeratedVariable** scope_in, int arity_in, Cost defval)
			: AbstractGlobalConstraint(wcsp, scope_in, arity_in),
				projectedCost(0, &wcsp->getStore()->storeCost), 
				nonassigned(arity_in, &wcsp->getStore()->storeInt),
				currentDepth(-1) 
{
	deltaCost = new vector<StoreCost>[arity_];
	for (int i=0;i<arity_in;i++) {
		deltaCost[i] = vector<StoreCost>(scope_in[i]->getDomainInitSize(), StoreCost(0,&(wcsp->getStore())->storeCost));
		fill(deltaCost[i].begin(), deltaCost[i].end(), 0);
	}

	count_nic = 0;
	count_gac = 0;
	count_fdac = 0;
	count_edac = 0;
	error = 0;

	fullySupportedSet = new set<int>[arity_in];

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
	/*if ((currentDepth == -1) || (currentDepth >= wcsp->getStore()->getDepth())) {
		initStructure();
	} else {
		vector<int> rmv;
		checkRemoved(rmv);
	}
	currentDepth = wcsp->getStore()->getDepth();*/
	initStructure();
}

Cost GlobalConstraint::eval(String s) {

	Cost tcost = evalOriginal(s);
	if (tcost < wcsp->getUb()) {
		for (unsigned int i=0;i<s.length();i++) tcost -= deltaCost[i][s[i]-CHAR_FIRST];
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
			Char tbuf[arity_ + 1]; 
			for(int i=0;i<arity_;i++) { 		
				tbuf[i] = getVar(i)->getValue() + CHAR_FIRST;
			}
			tbuf[arity_] =  '\0';
			String t = String(tbuf);
			projectLB(eval(t));
		} else { 
			pushAll();
		}
	}
}

void GlobalConstraint::project(int index, Value value, Cost cost) {
	EnumeratedVariable* x = (EnumeratedVariable*)getVar(index);
	if (deconnected()) return;
    // hard binary constraint costs are not changed
    if (!CUT(cost + wcsp->getLb(), wcsp->getUb())) {
	    TreeDecomposition* td = wcsp->getTreeDec();
	    if(td) td->addDelta(cluster,x,value,cost);
    	deltaCost[index][value] += cost;  // Warning! Possible overflow???
    }
	x->project(value, cost);
}

void GlobalConstraint::extend(int index, Value value, Cost cost) {
	EnumeratedVariable* x = (EnumeratedVariable*)getVar(index);
    TreeDecomposition* td = wcsp->getTreeDec();
    if(td) td->addDelta(cluster,x,value,-cost);
	deltaCost[index][value] -= cost;  // Warning! Possible overflow???
	x->extend(value, cost);
}

// function used for propagation
void GlobalConstraint::remove(int index) {
	vector<int> rmv;
	currentVar = -1;
	needPropagateDAC = true;
	needPropagateAC = true;
}

void GlobalConstraint::projectFromZero(int index) {
	vector<int> rmv;
	currentVar = -1;
	needPropagateDAC = true;
}

void GlobalConstraint::propagate() {

	if (deconnected()) return;
	for(int i=0;connected() && i<arity_;i++) {
	  if (getVar(i)->assigned()) assign(i);
	}

	needPropagateEAC = false;
	currentVar = -1;

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
		    // propagateEAC(); // Not implemented yet!
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

/// \warning Not implemented yet!
void GlobalConstraint::propagateEAC() {
	/*vector<int> provide;
	provide.resize(arity_);
	vector<map<Value, Cost> > deltas(arity_);
	for (int i=0;i<arity_;i++) {
		EnumeratedVariable *x = (EnumeratedVariable*)getVar(i);
		for (EnumeratedVariable::iterator v = x->begin();v != x->end();++v) {
			provide[i] = i;
			deltas[i][*v] = x->getCost(*v);
		}
	}
	changeAfterExtend(provide, deltas);

	for (int i=0;i<arity_;i++) {
		if (getVar(i)->unassigned()) {
			map<Value, Cost> delta;
			EnumeratedVariable *x = (EnumeratedVariable*)getVar(i);
			findProjection(i, delta);
			for (map<Value, Cost>::iterator p = delta.begin(); p != delta.end();p++) {
				if (wcsp->getLb() + p->second >= wcsp->getUb()) {
					x->remove(p->first);
				}
			}
		}
	}
	undoExtend();*/
}

void GlobalConstraint::propagateDAC() {
	vector<int> rmv;
	for (int i=0;i<arity_;i++) {
		if (getVar(i)->unassigned()) {
			checkRemoved(rmv);
			findFullSupport(i);
		}
	}
}

void GlobalConstraint::propagateAC() {
	vector<int> rmv;
	for (int i=0;i<arity_;i++) {
		if (getVar(i)->unassigned()) {
			checkRemoved(rmv);
			findSupport(i); 
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
	wcsp->revise(this);
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
	if (currentVar != index) {
		currentVar = index;
		vector<int> rmv;
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
		undoExtend();
	}
	if (EACCost[a] < 0) {
		cout << "EAC error\n";
		EACCost[a] = 0;
	}
	if (EACCost[a] > 0) needPropagateEAC = true;
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

void GlobalConstraint::findFullSupportEAC(int index) {
	currentVar = -1;
	if (needPropagateEAC) {
		needPropagateEAC = false;
		vector<int> provide;
		for (set<int>::iterator i = fullySupportedSet[index].begin();i !=
				fullySupportedSet[index].end();i++) {
			if (getVar(*i)->unassigned() && (*i != index)) provide.push_back(*i);
		}
		findFullSupport(index, provide, true);
	} 
}

void GlobalConstraint::findFullSupport(int varindex, vector<int> &support, bool isEAC)
{
	wcsp->revise(this);
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
		undoExtend();
		return;
	} else {
		count_fdac++;
	}

	changeAfterProject(varindex, delta);

	map<Value,Cost> *E = new map<Value, Cost>[support.size()];

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
	delete[] E;

}

void GlobalConstraint::findSupport(int varindex)
{	
	wcsp->revise(this);
	if (ToulBar2::verbose >= 3) cout << "findSupport for variable " << varindex << endl;
	map<Value, Cost> delta;	
	findProjection(varindex, delta);
	bool allzero = true;
	for (map<Value, Cost>::iterator i = delta.begin(); i != delta.end();i++) {
		if (i->second > 0) allzero = false;
	} 
	if (!allzero) {
		count_gac++;
		changeAfterProject(varindex, delta);
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
