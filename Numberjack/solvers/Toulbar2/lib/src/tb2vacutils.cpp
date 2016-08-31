/*
 * ****** Set of useful classes to enforce VAC
 */

#include "tb2vacutils.hpp"
#include "tb2vac.hpp"
#include "tb2clusters.hpp"



VACVariable::VACVariable (WCSP *wcsp, string n, Value iinf, Value isup) : EnumeratedVariable(wcsp, n, iinf, isup), vac(wcsp->vac), myThreshold(MIN_COST)
{
    init();
}

VACVariable::VACVariable (WCSP *wcsp, string n, Value *d, int dsize) : EnumeratedVariable(wcsp, n, d, dsize), vac(wcsp->vac), myThreshold(MIN_COST)
{
    init();
}

VACVariable::~VACVariable () {
}

void VACVariable::init () {
    killed =0; // HEUR
    maxk_timeStamp = 0;
    maxk = 0;
    for (unsigned int a = 0; a < getDomainInitSize(); a++) {
        mark.push_back(0);
        k_timeStamp.push_back(0);
        k.push_back(0);
        killer.push_back(0);
    }
    linkVACQueue.content.var = this;
    linkVACQueue.content.timeStamp = -1;
    linkSeekSupport.content.var = this;
    linkSeekSupport.content.timeStamp = -1;
    linkVAC2Queue.content = this;
}


bool VACVariable::increaseVAC(Value newInf) {
    if (newInf > inf) {
        if(newInf > sup) return true;
        else {
            newInf = domain.increase(newInf);
            inf = newInf;
        }
    }
    return false;
}

bool VACVariable::decreaseVAC(Value newSup) {
    if (newSup < sup) {
        if(newSup < inf) return true;
        else {
            newSup = domain.decrease(newSup);
            sup = newSup;
        }
    }
    return false;
}

bool VACVariable::removeVAC ( Value v )
{
    if (v == inf) return increaseVAC(v + 1);
    else if (v == sup) return decreaseVAC(v - 1);
    else if (canbe(v)) domain.erase(v);
    return false;
}

void VACVariable::decreaseCost(Value v, Cost c) {	
    assert(c > MIN_COST);
    Cost cini = getCost(v);
    if (wcsp->getLb() + cini < wcsp->getUb()) {
        costs[toIndex(v)] -= c;
    }
}


void VACVariable::VACproject (Value v, const Cost c) {
    //   Cost oldCost = getVACCost(v);
    costs[toIndex(v)] += c;
    //   Cost newCost = getVACCost(v);

    //   if ((v == maxCostValue) || (newCost > maxCost) || CUT(wcsp->getLb() + newCost,wcsp->getUb())) {
    if (CUT(wcsp->getLb() + getCost(v), wcsp->getUb())) {
        queueNC();
    }
    //   if (oldCost == MIN_COST) {
    //     queueNC();
    //     queueDAC();
    //     queueEAC1();
    //   }
    //   if ((isNull(oldCost)) && (!isNull(newCost))) {
    //     queueVAC2();
    //   }

    //   if(v == getSupport()) {
    //     Value newSupport = getInf();
    //     Cost minCost = getCost(newSupport);
    //     EnumeratedVariable::iterator iter = begin();
    //     for (++iter; minCost > MIN_COST && iter != end(); ++iter) {
    //         Cost cost = getCost(*iter);
    //         if (cost < minCost) {
    //             minCost = cost;
    //             newSupport = *iter;
    //         }
    //     }
    // 	assert(canbe(newSupport));
    // 	//	cout << "setsupport " << wcspIndex << " " << newSupport << endl;
    // 	setSupport(newSupport);
    //   }
}

void VACVariable::VACextend (Value v, const Cost c) {
    decreaseCost(v,c);
    if (v == maxCostValue) queueNC();
    assert(canbe(getSupport()));
    //   if(cannotbe(getSupport()) || getCost(getSupport())>MIN_COST) { // TO BE REMOVED ???
    //     Value newSupport = getInf();
    //     Cost minCost = getCost(newSupport);
    //     EnumeratedVariable::iterator iter = begin();
    //     for (++iter; minCost > MIN_COST && iter != end(); ++iter) {
    //         Cost cost = getCost(*iter);
    //         if (cost < minCost) {
    //             minCost = cost;
    //             newSupport = *iter;
    //         }
    //     }
    // 	assert(canbe(newSupport));
    // 	//	cout << "setsupport " << wcspIndex << " " << newSupport << endl;
    // 	setSupport(newSupport);
    //   }
}

bool VACVariable::isSimplyNull(Cost c) 
{
    return (vac->isNull(c));
}

bool VACVariable::isNull (Cost c) 
{
    return (vac->isNull(c) || (c < myThreshold));
}

void VACVariable::queueVAC() {
    wcsp->vac->queueVAC(&linkVACQueue);
}

void VACVariable::queueSeekSupport() {
    wcsp->vac->queueSeekSupport(&linkSeekSupport);
}

void VACVariable::queueVAC2() {
    wcsp->vac->queueVAC2(&linkVAC2Queue);
}

void VACVariable::dequeueVAC2() {
    wcsp->vac->dequeueVAC2(&linkVAC2Queue);
}

// void VACVariable::extendAll(Cost cost) {
//   VACVariable *xj;
//   if (ToulBar2::vac) {
//     for (ConstraintList::iterator iter = getConstrs()->begin(); iter != getConstrs()->end(); ++iter) {
// 	   Constraint *c = (*iter).constr;
//        if(c->arity() == 2 && !c->isSep()) {
// 		   xj = (VACVariable *) c->getVar(1 - (*iter).scopeIndex);
// 	       xj->queueVAC2();
//        }
//     }
//   }
//   EnumeratedVariable::extendAll(cost);
// }

// void VACVariable::assign(Value newValue) {
//   vac->assign(wcspIndex, newValue);

//   if (ToulBar2::vac) {
//     for (ConstraintList::iterator iter = getConstrs()->begin(); iter != getConstrs()->end(); ++iter) {
// 	   Constraint *c = (*iter).constr;
//        if(c->arity() == 2 && !c->isSep()) {
// 	       VACVariable *xj = (VACVariable *) c->getVar(1 - (*iter).scopeIndex);
// 	       xj->queueVAC2();
//        }
//     }
//   }
//   EnumeratedVariable::assign(newValue);
// }


void VACVariable::remove (Value value) {
    // if (canbe(value)) {
    //   queueVAC2();
    // }
    if (ToulBar2::singletonConsistency) vac->singleton.insert(MAX_DOMAIN_SIZE*wcspIndex + value);
    EnumeratedVariable::remove(value);
}


void VACVariable::removeFast (Value value) {
    if (ToulBar2::singletonConsistency) vac->singleton.insert(MAX_DOMAIN_SIZE*wcspIndex + value);
    EnumeratedVariable::removeFast(value);
}


// void VACVariable::project (Value value, Cost cost) {
//   assert(cost >= MIN_COST);
//   Cost oldcost = getCost(value);
//   Cost newcost = oldcost + cost;
//   if ((isNull(oldcost)) && (!isNull(newcost))) {
//     queueVAC2();
//   }
//   EnumeratedVariable::project(value, cost);
// }

// void VACVariable::extend (Value value, Cost cost) {
//   VACVariable *xj;
//   queueVAC2();
//   for (ConstraintList::iterator iter = getConstrs()->begin(); iter != getConstrs()->end(); ++iter) {
// 	Constraint *c = (*iter).constr;
// 	if(c->arity() == 2 && !c->isSep()) {
// 	  xj = (VACVariable *) c->getVar(1 - (*iter).scopeIndex);
// 	  xj->queueVAC2();
// 	}
//   }
//   EnumeratedVariable::extend(value, cost);
// }

void VACVariable::increase(Value newInf) {
    // if ((newInf > inf) && (newInf < sup)) {
    //   queueVAC2();
    // }
    if (ToulBar2::singletonConsistency) for(int i=inf;i<=newInf;i++) vac->singleton.insert(MAX_DOMAIN_SIZE*wcspIndex + i);
    EnumeratedVariable::increase(newInf);
}

void VACVariable::decrease(Value newSup) {
    // if ((newSup < sup) && (newSup > inf)) {
    //   queueVAC2();
    // }
    if (ToulBar2::singletonConsistency) for(int i=sup;i>=newSup;i--) vac->singleton.insert(MAX_DOMAIN_SIZE*wcspIndex + i);
    EnumeratedVariable::decrease(newSup);
}


/******************************
 * min-sum diffusion algorithm
 */

bool VACVariable::averaging()
{
    Cost Top = wcsp->getUb();
    bool change = false;
    EnumeratedVariable* x;
    EnumeratedVariable* y;
    Constraint* ctr = NULL;
    ConstraintList::iterator itc = getConstrs()->begin();
    if(itc != getConstrs()->end())	ctr = (*itc).constr;
    while(ctr) {
        if(ctr->arity() == 2 && !ctr->isSep()) {
            BinaryConstraint* bctr = (BinaryConstraint*) ctr;
            x = (EnumeratedVariable*) bctr->getVarDiffFrom( (Variable*) this );
            for (iterator it = begin(); it != end(); ++it) {
                Cost cu = getCost(*it);
                Cost cmin = Top;
                for (iterator itx = x->begin(); itx != x->end(); ++itx) {
                    Cost cbin = bctr->getCost(this,x,*it,*itx);
                    if(cbin < cmin) cmin = cbin;
                }
                assert(cmin < Top);
                Double mean = to_double(cmin + cu) / 2.;
                Double extc = to_double(cu) - mean;
                if(abs(extc) >= 1) {
                    Cost costi = (Long) extc;
                    for (iterator itx = x->begin(); itx != x->end(); ++itx) {
                        bctr->addcost(this,x,*it,*itx,costi);
                    }
                    if(mean > to_double(cu)) project(*it, -costi);
                    else extend(*it, costi);
                    change = true;
                }
            }
        } else if(ctr->arity() == 3 && !ctr->isSep()) {
            TernaryConstraint* tctr = (TernaryConstraint*) ctr;
            x = (EnumeratedVariable*) tctr->getVar( 0 );
            if(x == this) x = (EnumeratedVariable*) tctr->getVar( 1 );
            y = (EnumeratedVariable*) tctr->getVarDiffFrom((Variable*) this, (Variable*)x);
            for (iterator it = begin(); it != end(); ++it) {
                Cost cu = getCost(*it);
                Cost cmin = Top;
                for (iterator itx = x->begin(); itx != x->end(); ++itx) {
                    for (iterator ity = y->begin(); ity != y->end(); ++ity) {
                        Cost ctern = tctr->getCost(this,x,y,*it,*itx,*ity);
                        if(ctern < cmin) cmin = ctern;
                    }}
                assert(cmin < Top);
                Double mean = to_double(cmin + cu) / 2.;
                Double extc = to_double(cu) - mean;
                if(abs(extc) >= 1) {
                    Cost costi = (Long) extc;
                    for (iterator itx = x->begin(); itx != x->end(); ++itx) {
                        for (iterator ity = y->begin(); ity != y->end(); ++ity) {
                            tctr->addCost(this,x,y,*it,*itx,*ity,costi);
                        }}
                    if(mean > to_double(cu)) project(*it, -costi);
                    else extend(*it, costi);
                    change = true;
                }
            }
        } else if(ctr->arity() >= 4 && ctr->extension() && !ctr->isSep()) {
            NaryConstraint* nctr = (NaryConstraint*) ctr;
            for (iterator it = begin(); it != end(); ++it) {
                Cost cu = getCost(*it);
                Cost cmin = Top;
                int tindex = nctr->getIndex(this);
                String tuple;
                Cost cost;
                Long nbtuples = 0;
                nctr->first();
                while (nctr->next(tuple,cost)) {
                    nbtuples++;
                    if (toValue(tuple[tindex] - CHAR_FIRST)==(*it) && cost < cmin) cmin = cost;
                }
                if (nctr->getDefCost() < cmin && nbtuples < nctr->getDomainSizeProduct()/getDomainSize()) cmin = nctr->getDefCost();
                //				assert(cmin < Top);
                Double mean = to_double(cmin + cu) / 2.;
                Double extc = to_double(cu) - mean;
                if(abs(extc) >= 1) {
                    Cost costi = (Cost) extc;
                    nctr->addtoTuples(this, *it, costi);
                    if(mean > to_double(cu)) project(*it, -costi);
                    else extend(*it, costi);
                    change = true;
                }
            }
        }
        ++itc;
        if(itc != getConstrs()->end()) ctr = (*itc).constr;
        else ctr = NULL;
    }
    return change;
}


/************************************************************
 * VACBinaryConstraint:
 *   A class that stores information about a binary cost function
 */

VACBinaryConstraint::VACBinaryConstraint (WCSP *wcsp, EnumeratedVariable *xx, EnumeratedVariable *yy, vector<Cost> &tab) :  BinaryConstraint(wcsp, xx, yy, tab), myThreshold(MIN_COST)
{
    for (unsigned int a = 0; a < xx->getDomainInitSize(); a++) {
        kX.push_back(0);
        kX_timeStamp.push_back(0);
    }
    for (unsigned int b = 0; b < yy->getDomainInitSize(); b++) {
        kY.push_back(0);
        kY_timeStamp.push_back(0);
    }
}

VACBinaryConstraint::VACBinaryConstraint (WCSP *wcsp) : BinaryConstraint(wcsp) , myThreshold(MIN_COST)
{}

void VACBinaryConstraint::VACfillElimConstr ()
{
    for (unsigned int a = kX.size(); a < sizeX; a++) {
        kX.push_back(0);
        kX_timeStamp.push_back(0);
    }
    for (unsigned int b = kY.size(); b < sizeY; b++) {
        kY.push_back(0);
        kY_timeStamp.push_back(0);
    }
}

VACBinaryConstraint::~VACBinaryConstraint ()
{
}

void VACBinaryConstraint::VACproject (VACVariable* x, Value v, Cost c) {
    assert(ToulBar2::verbose < 4 || ((cout << "project(C" << getVar(0)->getName() << "," << getVar(1)->getName() << ", (" << x->getName() << "," << v << "), " << c << ")" << endl), true));

    TreeDecomposition* td = wcsp->getTreeDec();
    if(td) td->addDelta(cluster,x,v,c);

    int index = x->toIndex(v);
    // TO BE REPLACED BY A LOOP ON THE DOMAIN IN ORDER TO AVOID SUBTRACTING TOP???
    if(!getIndex(x)) deltaCostsX[index] += c;
    else             deltaCostsY[index] += c;
    assert(getCost((EnumeratedVariable *) x,(EnumeratedVariable *) getVarDiffFrom(x),v,getVarDiffFrom(x)->getInf()) >= MIN_COST);
    assert(getCost((EnumeratedVariable *) x,(EnumeratedVariable *) getVarDiffFrom(x),v,getVarDiffFrom(x)->getSup()) >= MIN_COST);
    x->VACproject(v, c);
}

void VACBinaryConstraint::VACextend(VACVariable* x, Value v, Cost c) {
    assert(ToulBar2::verbose < 4 || ((cout << "extend(C" << getVar(0)->getName() << "," << getVar(1)->getName() << ", (" << x->getName() << "," << v << "), " << c << ")" << endl), true));

    TreeDecomposition* td = wcsp->getTreeDec();
    if(td) td->addDelta(cluster,x,v,-c);

    int index = x->toIndex(v);
    // TO BE REPLACED BY A LOOP ON THE DOMAIN IN ORDER TO AVOID SUBTRACTING TOP???
    if(!getIndex(x)) deltaCostsX[index] -= c;
    else             deltaCostsY[index] -= c;
    x->VACextend(v, c);
}

int VACBinaryConstraint::getK (VACVariable* var, Value v, Long timeStamp) {
    if(var == (VACVariable*) getVar(0)) {
        if(kX_timeStamp[var->toIndex(v)] < timeStamp) return 0;
        else return kX[var->toIndex(v)];
    }  else  {
        if(kY_timeStamp[var->toIndex(v)] < timeStamp) return 0;
        else return kY[var->toIndex(v)];
    }
}

void VACBinaryConstraint::setK (VACVariable* var, Value v, int c, Long timeStamp) {
    if(var == getVar(0)) {
        kX[var->toIndex(v)] = c;
        kX_timeStamp[var->toIndex(v)] = timeStamp;
    } else {
        kY[var->toIndex(v)] = c;
        kY_timeStamp[var->toIndex(v)] = timeStamp;
    }
}

bool VACBinaryConstraint::isNull (Cost c)
{
    VACVariable* xi = (VACVariable*) getVar(0);
    return (xi->isSimplyNull(c) || (c < myThreshold));
}

bool VACBinaryConstraint::revise (VACVariable* var, Value v) {
    //  bool wipeout = false;
    VACVariable* xi = (VACVariable*) getVar(0);
    VACVariable* xj = (VACVariable*) getVar(1);
    Value sup = getSupport(var,v);
    Value minsup = sup;
    if(var != xi) {  xi = (VACVariable*)getVar(1); xj = (VACVariable*)getVar(0); }
    Cost cost, minCost = wcsp->getUb();

    if(xj->canbe(sup)) {
        if(xj->getVACCost(sup) != MIN_COST) { xj->removeVAC(sup);  } // wipeout = xj->removeVAC(sup);
        else {
            if (getVACCost(xi,xj,v,sup) == MIN_COST) {
                return false;
            }
        }
    }

    for (EnumeratedVariable::iterator it = xj->lower_bound(sup); it != xj->end(); ++it) {
        Value w = *it;
        if(xj->getVACCost(w) != MIN_COST) { xj->removeVAC(w); xj->queueVAC(); } // wipeout = xj->removeVAC(w); xj->queueVAC();
        else {
            cost = getVACCost(xi,xj,v, w);
            if (cost == MIN_COST) {
                setSupport(xi,v,w);
                return false;
            } else if (cost < minCost) {
                minCost = cost;
                minsup = w;
            }
        }
    }
    for (EnumeratedVariable::iterator it = xj->begin(); it != xj->lower_bound(sup); ++it) {
        Value w = *it;
        if(xj->getVACCost(w) != MIN_COST) { xj->removeVAC(w); xj->queueVAC(); } // wipeout = xj->removeVAC(w); xj->queueVAC();
        else {
            cost = getVACCost(xi,xj,v, w);
            if (cost == MIN_COST) {
                setSupport(xi,v,w);
                return false;
            } else if (cost < minCost) {
                minCost = cost;
                minsup = w;
            }
        }
    }

    setSupport(xi,v,minsup);
    return true;
}


/************************************************************
 * VACTernaryConstraint:
 *   A class that stores information about a ternary cost function
 */
/*
VACTernaryConstraint::VACTernaryConstraint (WCSP *wcsp, EnumeratedVariable *xx, EnumeratedVariable *yy, EnumeratedVariable *zz, BinaryConstraint *xy, BinaryConstraint *xz, BinaryConstraint *yz, vector<Cost> &tab) :  TernaryConstraint(wcsp, xx, yy, zz, xy, xz, yz, tab, storeCost), myThreshold(MIN_COST)
{
   for (int a = 0; a < xx->getDomainInitSize(); a++) {
	   	kX.push_back(0);
	   	kX_timeStamp.push_back(0);
   }
   for (int b = 0; b < yy->getDomainInitSize(); b++) {
	   	kY.push_back(0);
	   	kY_timeStamp.push_back(0);
   }
   for (int c = 0; c < zz->getDomainInitSize(); c++) {
	   	kZ.push_back(0);
	   	kZ_timeStamp.push_back(0);
   }
}

VACTernaryConstraint::VACTernaryConstraint (WCSP *wcsp) : TernaryConstraint(wcsp) , myThreshold(MIN_COST)
{
   for (int a = 0; a < wcsp->getMaxDomainSize(); a++) {
	   	kX.push_back(0);
	   	kX_timeStamp.push_back(0);
   }
   for (int b = 0; b < wcsp->getMaxDomainSize(); b++) {
	   	kY.push_back(0);
	   	kY_timeStamp.push_back(0);
   }
   for (int c = 0; c < wcsp->getDomainInitSize(); c++) {
	   	kZ.push_back(0);
	   	kZ_timeStamp.push_back(0);
   }
}

VACTernaryConstraint::~VACTernaryConstraint ()
{
}

void VACTernaryConstraint::VACproject (VACVariable* x, Value v, Cost c) {
  assert(ToulBar2::verbose < 4 || ((cout << "project(C" << getVar(0)->getName() << "," << getVar(1)->getName() << "," << getVar(2)->getName() << ", (" << x->getName() << "," << v << "), " << c << ")" << endl), true));

  TreeDecomposition* td = wcsp->getTreeDec();
  if(td) td->addDelta(cluster,x,v,c);

  int index = x->toIndex(v);
  // TO BE REPLACED BY A LOOP ON THE DOMAIN IN ORDER TO AVOID SUBTRACTING TOP???
  if(getIndex(x)==0) deltaCostsX[index] += c;
  else if(getIndex(x)==1) deltaCostsY[index] += c;
  else deltaCostsZ[index] += c;
  x->VACproject(v, c);
}

void VACTernaryConstraint::VACextend(VACVariable* x, Value v, Cost c) {
  assert(ToulBar2::verbose < 4 || ((cout << "extend(C" << getVar(0)->getName() << "," << getVar(1)->getName() << "," << getVar(2)->getName() << ", (" << x->getName() << "," << v << "), " << c << ")" << endl), true));

  TreeDecomposition* td = wcsp->getTreeDec();
  if(td) td->addDelta(cluster,x,v,-c);

  int index = x->toIndex(v);
  // TO BE REPLACED BY A LOOP ON THE DOMAIN IN ORDER TO AVOID SUBTRACTING TOP???
  if(getIndex(x)==0) deltaCostsX[index] -= c;
  else if(getIndex(x)==1) deltaCostsY[index] -= c;
  else deltaCostsZ[index] -= c;
  x->VACextend(v, c);
}

int VACTernaryConstraint::getK (VACVariable* var, Value v, Long timeStamp) {
  if(var == (VACVariable*) getVar(0)) {
  	if(kX_timeStamp[var->toIndex(v)] < timeStamp) return 0;
  	else return kX[var->toIndex(v)];
  }  else if(var == (VACVariable*) getVar(1)) {
  	if(kY_timeStamp[var->toIndex(v)] < timeStamp) return 0;
  	else return kY[var->toIndex(v)];
  }  else {
	if(kZ_timeStamp[var->toIndex(v)] < timeStamp) return 0;
	else return kZ[var->toIndex(v)];
  }
}

void VACTernaryConstraint::setK (VACVariable* var, Value v, int c, Long timeStamp) {
  if(var == getVar(0)) {
  	kX[var->toIndex(v)] = c;
  	kX_timeStamp[var->toIndex(v)] = timeStamp;
  } else if(var == getVar(1)) {
    kY[var->toIndex(v)] = c;
   	kY_timeStamp[var->toIndex(v)] = timeStamp; 					   
  }	 else {
	kZ[var->toIndex(v)] = c;
	kZ_timeStamp[var->toIndex(v)] = timeStamp;
  }
}

bool VACTernaryConstraint::isNull (Cost c)
{
  VACVariable* xi = (VACVariable*) getVar(0);
  return (xi->isSimplyNull(c) || (c < myThreshold));
}

bool VACTernaryConstraint::revise (VACVariable* var, Value v) {
  bool wipeout = false;
  VACVariable* xi = (VACVariable*) getVar(0);
  VACVariable* xj = (VACVariable*) getVar(1);
  VACVariable* xl = (VACVariable*) getVar(2);
  pair<Value,Value> sup = getSupport(var,v);
  pair<Value,Value> minsup = sup;
  if(var != xi) {
	  if (var != xj) {
		 xi = (VACVariable*)getVar(2); xj = (VACVariable*)getVar(0); xl = (VACVariable*)getVar(1);
	  } else {
		 xi = (VACVariable*)getVar(1); xj = (VACVariable*)getVar(0);
	  }
  }
  Cost cost, minCost = wcsp->getUb();

  assert(getindex(xj) < getindex(xl)); // check support is correctly oriented w.r.t. xj/first and xl/second
  if(xj->canbe(sup.first) && xl->canbe(sup.second)) {
	  bool unarytest = true;
	  if(xj->getVACCost(sup.first) != MIN_COST) { wipeout = xj->removeVAC(sup.first);  unarytest= false;}
	  if(xl->getVACCost(sup.second) != MIN_COST) { wipeout = xl->removeVAC(sup.second); unarytest= false;}
	  if (unarytest) {
		  if (getVACCost(xi,xj,xl,v,sup.first,sup.second) == MIN_COST) {
		    return false;
		  }
	  }
  }

  for (EnumeratedVariable::iterator it = xj->lower_bound(sup); it != xj->end(); ++it) {
	  Value w = *it;	
	  if(xj->getVACCost(w) != MIN_COST) { wipeout = xj->removeVAC(w); xj->queueVAC(); }
	  else {
	      cost = getVACCost(xi,xj,v, w);
	      if (cost == MIN_COST) {		
	      	setSupport(xi,v,w);
	        return false;
	      } else if (cost < minCost) {
	      	  minCost = cost;
	          minsup = w;
	      }
	  }
  }
  for (EnumeratedVariable::iterator it = xj->begin(); it != xj->lower_bound(sup); ++it) {
	  Value w = *it;	
	  if(xj->getVACCost(w) != MIN_COST) { wipeout = xj->removeVAC(w); xj->queueVAC(); }
	  else {
	      cost = getVACCost(xi,xj,v, w);
	      if (cost == MIN_COST) {		
	      	setSupport(xi,v,w);
	        return false;
	      } else if (cost < minCost) {
	      	  minCost = cost;
	          minsup = w;
	      }
	  }
  }

  setSupport(xi,v,minsup);
  return true;
}
 */

/* Local Variables: */
/* c-basic-offset: 4 */
/* tab-width: 4 */
/* indent-tabs-mode: nil */
/* c-default-style: "k&r" */
/* End: */

