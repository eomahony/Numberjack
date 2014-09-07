/*
 * ****** Binary constraint applied on variables with enumerated domains ******
 */

#include "tb2ternaryconstr.hpp"
#include "tb2wcsp.hpp"
#include "tb2vac.hpp"
#include "tb2clusters.hpp"

/*
 * Constructors and misc.
 * 
 */

TernaryConstraint::TernaryConstraint(WCSP *wcsp, 
									 EnumeratedVariable *xx, 
									 EnumeratedVariable *yy,
									 EnumeratedVariable *zz,
									 BinaryConstraint *xy_,
									 BinaryConstraint *xz_,
									 BinaryConstraint *yz_,
									 vector<Cost> &tab, 
									 StoreStack<Cost, Cost> *storeCost) 
									 
					: AbstractTernaryConstraint<EnumeratedVariable,EnumeratedVariable,EnumeratedVariable>(wcsp, xx, yy, zz), 
					  sizeX(xx->getDomainInitSize()), sizeY(yy->getDomainInitSize()), sizeZ(zz->getDomainInitSize()),
					  functionalX(true), functionalY(true), functionalZ(true)
{
    assert(tab.size() == sizeX * sizeY * sizeZ);
    deltaCostsX = vector<StoreCost>(sizeX,StoreCost(MIN_COST,storeCost));
    deltaCostsY = vector<StoreCost>(sizeY,StoreCost(MIN_COST,storeCost));
    deltaCostsZ = vector<StoreCost>(sizeZ,StoreCost(MIN_COST,storeCost));
	assert(getIndex(x) < getIndex(y) && getIndex(y) < getIndex(z));
	functionX = vector<Value>(sizeY * sizeZ, WRONG_VAL);
	functionY = vector<Value>(sizeX * sizeZ, WRONG_VAL);
	functionZ = vector<Value>(sizeX * sizeY, WRONG_VAL);
    supportX = vector< pair<Value,Value> >(sizeX,make_pair(y->getInf(),z->getInf()));
    supportY = vector< pair<Value,Value> >(sizeY,make_pair(x->getInf(),z->getInf()));
    supportZ = vector< pair<Value,Value> >(sizeZ,make_pair(x->getInf(),y->getInf()));

//    costs = vector<StoreCost>(sizeX*sizeY*sizeZ,StoreCost(MIN_COST,storeCost));
    
    for (unsigned int a = 0; a < x->getDomainInitSize(); a++) {
	  for (unsigned int b = 0; b < y->getDomainInitSize(); b++) {
		for (unsigned int c = 0; c < z->getDomainInitSize(); c++) {
		  Cost cost = tab[a * sizeY * sizeZ + b * sizeZ + c];
//		  costs[a * sizeY * sizeZ + b * sizeZ + c] = cost;
		  if (!CUT(cost, wcsp->getUb()) && x->canbe(x->toValue(a)) && y->canbe(y->toValue(b)) && z->canbe(z->toValue(c))) {
			if (functionalX) {
			  if (functionX[b * sizeZ + c] == WRONG_VAL) functionX[b * sizeZ + c] = x->toValue(a);
			  else functionalX = false;
			}
			if (functionalY) {
			  if (functionY[a * sizeZ + c] == WRONG_VAL) functionY[a * sizeZ + c] = y->toValue(b);
			  else functionalY = false;
			}
			if (functionalZ) {
			  if (functionZ[a * sizeY + b] == WRONG_VAL) functionZ[a * sizeY + b] = z->toValue(c);
			  else functionalZ = false;
			}
		  }
		}
	  }
	}
 
    xy = xy_;
    xz = xz_;
    yz = yz_;
                
    if (functionalX) {
    	costsYZ = vector<StoreCost>(sizeY*sizeZ,StoreCost(MIN_COST,storeCost));
		for (unsigned int b = 0; b < y->getDomainInitSize(); b++) {
			for (unsigned int c = 0; c < z->getDomainInitSize(); c++) {
				if (functionX[b * sizeZ + c] != WRONG_VAL) costsYZ[b * sizeZ + c] = tab[x->toIndex(functionX[b * sizeZ + c]) * sizeY * sizeZ + b * sizeZ + c];
			}
		}
//    	costs.free_all();
    } else {
    	costs = vector<StoreCost>(sizeX*sizeY*sizeZ,StoreCost(MIN_COST,storeCost));
        for (unsigned int a = 0; a < x->getDomainInitSize(); a++) {
    	  for (unsigned int b = 0; b < y->getDomainInitSize(); b++) {
    		for (unsigned int c = 0; c < z->getDomainInitSize(); c++) {
    		  costs[a * sizeY * sizeZ + b * sizeZ + c] = tab[a * sizeY * sizeZ + b * sizeZ + c];
    		}
    	  }
        }
    }

    vector<int> &vecX = wcsp->getListSuccessors()->at(xx->wcspIndex);
    vector<int> &vecY = wcsp->getListSuccessors()->at(yy->wcspIndex);
    vector<int> &vecZ = wcsp->getListSuccessors()->at(zz->wcspIndex);
    if ((std::find(vecX.begin(), vecX.end(), yy->wcspIndex)==vecX.end()) && (std::find(vecX.begin(), vecX.end(), zz->wcspIndex)==vecX.end()) &&
    	(std::find(vecY.begin(), vecY.end(), xx->wcspIndex)==vecY.end()) && (std::find(vecY.begin(), vecY.end(), zz->wcspIndex)==vecY.end()) &&
    	(std::find(vecZ.begin(), vecZ.end(), xx->wcspIndex)==vecZ.end()) && (std::find(vecZ.begin(), vecZ.end(), yy->wcspIndex)==vecZ.end())) {
		switch (functionalX + functionalY + functionalZ) {
		case 1:
			if (functionalX) {
				vecX.push_back(yy->wcspIndex);
				vecX.push_back(zz->wcspIndex);
			} else if (functionalY) {
				vecY.push_back(xx->wcspIndex);
				vecY.push_back(zz->wcspIndex);
			} else if (functionalZ) {
				vecZ.push_back(xx->wcspIndex);
				vecZ.push_back(yy->wcspIndex);
			}
			ToulBar2::Berge_Dec = 1;
			break;
		case 2:
			if (functionalX && functionalY) {
				if (xx->wcspIndex < yy->wcspIndex) {
					vecX.push_back(zz->wcspIndex);
					vecZ.push_back(yy->wcspIndex);
				} else {
					vecY.push_back(zz->wcspIndex);
					vecZ.push_back(xx->wcspIndex);
				}
			} else if (functionalX && functionalZ) {
				if (xx->wcspIndex < zz->wcspIndex) {
					vecX.push_back(yy->wcspIndex);
					vecY.push_back(zz->wcspIndex);
				} else {
					vecZ.push_back(yy->wcspIndex);
					vecY.push_back(xx->wcspIndex);
				}
			} else if (functionalY && functionalZ) {
				if (yy->wcspIndex < zz->wcspIndex) {
					vecY.push_back(xx->wcspIndex);
					vecX.push_back(zz->wcspIndex);
				} else {
					vecZ.push_back(xx->wcspIndex);
					vecX.push_back(yy->wcspIndex);
				}
			}
			ToulBar2::Berge_Dec = 1;
			break;
		case 3:
			if (xx->wcspIndex < yy->wcspIndex && xx->wcspIndex < zz->wcspIndex) {
				if (yy->wcspIndex < zz->wcspIndex) {
					vecY.push_back(xx->wcspIndex);
					vecX.push_back(zz->wcspIndex);
				} else {
					vecZ.push_back(xx->wcspIndex);
					vecX.push_back(yy->wcspIndex);
				}
			} else if (yy->wcspIndex < xx->wcspIndex && yy->wcspIndex < zz->wcspIndex) {
				if (xx->wcspIndex < zz->wcspIndex) {
					vecX.push_back(yy->wcspIndex);
					vecY.push_back(zz->wcspIndex);
				} else {
					vecZ.push_back(yy->wcspIndex);
					vecY.push_back(xx->wcspIndex);
				}
			} else if (zz->wcspIndex < xx->wcspIndex && zz->wcspIndex < yy->wcspIndex) {
				if (xx->wcspIndex < yy->wcspIndex) {
					vecX.push_back(zz->wcspIndex);
					vecZ.push_back(yy->wcspIndex);
				} else {
					vecY.push_back(zz->wcspIndex);
					vecZ.push_back(xx->wcspIndex);
				}
			}
			ToulBar2::Berge_Dec = 1;
			break;
		default:
			break;
		}
    }

    propagate();
}


TernaryConstraint::TernaryConstraint(WCSP *wcsp, StoreStack<Cost, Cost> *storeCost) 
					: AbstractTernaryConstraint<EnumeratedVariable,EnumeratedVariable,EnumeratedVariable>(wcsp), 
					  sizeX(0), sizeY(0), sizeZ(0),
					  functionalX(false), functionalY(false), functionalZ(false)
{
//	unsigned int maxdom = wcsp->getMaxDomainSize();
//    deltaCostsX = vector<StoreCost>(maxdom,StoreCost(MIN_COST,storeCost));
//    deltaCostsY = vector<StoreCost>(maxdom,StoreCost(MIN_COST,storeCost));
//    deltaCostsZ = vector<StoreCost>(maxdom,StoreCost(MIN_COST,storeCost));
//    supportX = vector< pair<Value,Value> >(maxdom);
//    supportY = vector< pair<Value,Value> >(maxdom);
//    supportZ = vector< pair<Value,Value> >(maxdom);
	linkX = new DLink<ConstraintLink>;
    linkY = new DLink<ConstraintLink>;    
    linkZ = new DLink<ConstraintLink>;    
    
//    costs = vector<StoreCost>(maxdom*maxdom*maxdom,StoreCost(MIN_COST,storeCost));
//    for (unsigned int a = 0; a < maxdom; a++)
//       for (unsigned int b = 0; b < maxdom; b++)
//           for (unsigned int c = 0; c < maxdom; c++)
//               costs[a * maxdom * maxdom + b * maxdom + c] = MIN_COST;
 	xy = NULL;
 	xz = NULL;
 	yz = NULL;
 	
}


double TernaryConstraint::computeTightness()
{
   int count = 0;
   double sum = 0;
   Cost costs[x->getDomainSize()*y->getDomainSize()*z->getDomainSize()];
   for (EnumeratedVariable::iterator iterX = x->begin(); iterX != x->end(); ++iterX) {
      for (EnumeratedVariable::iterator iterY = y->begin(); iterY != y->end(); ++iterY) {
	      for (EnumeratedVariable::iterator iterZ = z->begin(); iterZ != z->end(); ++iterZ) {
	    	Cost c = getCost(*iterX, *iterY, *iterZ);
			sum += to_double(min(wcsp->getUb(), c));
    	    costs[count] = min(wcsp->getUb(), c);
			count++;
       }
     }
   }
   
   if (ToulBar2::weightedTightness == 2) {
   	tight = to_double(stochastic_selection<Cost>(costs, 0, count-1, count / 2));
   } else {
   	tight = sum / (double) count;
   }
	return tight;
}


void TernaryConstraint::print(ostream& os)
{
  os << this << " TernaryConstraint(" << x->getName() << ((functionalX)?"!":"") << "," << y->getName() << ((functionalY)?"!":"") << "," << z->getName() << ((functionalZ)?"!":"") << ")";
	if (ToulBar2::weightedDegree) os << "/" << getConflictWeight();
    if(wcsp->getTreeDec()) cout << "   cluster: " << getCluster() << endl;
    else cout << endl;
    
    if (ToulBar2::verbose >= 5) {
        for (EnumeratedVariable::iterator iterX = x->begin(); iterX != x->end(); ++iterX) {
            for (EnumeratedVariable::iterator iterY = y->begin(); iterY != y->end(); ++iterY) {
			   for (EnumeratedVariable::iterator iterZ = z->begin(); iterZ != z->end(); ++iterZ) {
		         os << " " << getCost(*iterX, *iterY, *iterZ);
			   }
                os << " ; ";
            }
            os << endl;
        }
    }
}

void TernaryConstraint::dump(ostream& os, bool original)
{
  os << "3 " << ((original)?(x->wcspIndex):x->getCurrentVarId()) << " " << ((original)?(y->wcspIndex):y->getCurrentVarId()) << " " << ((original)?(z->wcspIndex):z->getCurrentVarId()) << " " << MIN_COST << " " << x->getDomainSize() * y->getDomainSize() * z->getDomainSize() << endl;
  int i=0;
  for (EnumeratedVariable::iterator iterX = x->begin(); iterX != x->end(); ++iterX, i++) {
	int j=0;
	for (EnumeratedVariable::iterator iterY = y->begin(); iterY != y->end(); ++iterY, j++) {
	  int k=0;
	  for (EnumeratedVariable::iterator iterZ = z->begin(); iterZ != z->end(); ++iterZ, k++) {
		os << ((original)?(*iterX):i) << " " << ((original)?(*iterY):j) << " " << ((original)?(*iterZ):k) << " " << ((original)?getCost(*iterX, *iterY, *iterZ):min(wcsp->getUb(),getCost(*iterX, *iterY, *iterZ))) << endl;
	  }
	}
  }
}

/*
 * Propagation methods
 * 
 */
bool  TernaryConstraint::project(EnumeratedVariable *x, Value value, Cost cost, vector<StoreCost> &deltaCostsX)
{
	assert(ToulBar2::verbose < 4 || ((cout << "project(C" << getVar(0)->getName() << "," << getVar(1)->getName() << "," << getVar(2)->getName() << ", (" << x->getName() << "," << value << "), " << cost << ")" << endl), true));

    // hard ternary constraint costs are not changed
    if (!CUT(cost + wcsp->getLb(), wcsp->getUb())) {
	    TreeDecomposition* td = wcsp->getTreeDec();
	    if(td) td->addDelta(cluster,x,value,cost);
    	deltaCostsX[x->toIndex(value)] += cost;  // Warning! Possible overflow???
    }
    	
    Cost oldcost = x->getCost(value);
    x->project(value, cost);
#ifdef DEECOMPLETE
    int xindex = getIndex(x);
    getVar((xindex+1) % 3)->queueDEE();
    getVar((xindex+2) % 3)->queueDEE();
#endif
    return (x->getSupport() == value || SUPPORTTEST(oldcost, cost));
}

void  TernaryConstraint::extend(EnumeratedVariable *x, Value value, Cost cost, vector<StoreCost> &deltaCostsX)
{
	assert(ToulBar2::verbose < 4 || ((cout << "extend(C" << getVar(0)->getName() << "," << getVar(1)->getName() << "," << getVar(2)->getName() << ", (" << x->getName() << "," << value << "), " << cost << ")" << endl), true));
    TreeDecomposition* td = wcsp->getTreeDec();
    if(td) td->addDelta(cluster,x,value,-cost);
    deltaCostsX[x->toIndex(value)] -= cost;  // Warning! Possible overflow???
    x->extend(value, cost);
}

template <typename T1, typename T2, typename T3>
void TernaryConstraint::project(T1 getCost, T2 addCost, bool functionalZ, T3 getFunctionZ, BinaryConstraint* xy, EnumeratedVariable *x, EnumeratedVariable *y, EnumeratedVariable *z, Value valx, Value valy, Cost cost) {
	assert(ToulBar2::verbose < 4 || ((cout << "project(C" << getVar(0)->getName() << "," << getVar(1)->getName() << "," << getVar(2)->getName() << ", ((" << x->getName() << "," << valx << "),(" << y->getName() << "," << valy << ")), " << cost << ")" << endl), true));
	assert(cost >= MIN_COST);
	// BUG!
	// if (functionalZ) {
	//   Value valz = getFunctionZ(x,y, valx, valy);
	//   if (valz != WRONG_VAL && z->canbe(valz)) {
	// 	if (!CUT(getCost(x,y,z, valx,valy,valz), wcsp->getUb())) { // keeps forbidden costs into ternaries to get strong GAC3
	// 	  addCost(x,y,z,valx,valy,valz,-cost);
	// 	}
	//   }
	// } else {
	  for (EnumeratedVariable::iterator iterZ = z->begin(); iterZ != z->end(); ++iterZ) {
		if (!CUT(getCost(x,y,z, valx,valy,*iterZ), wcsp->getUb())) { // keeps forbidden costs into ternaries to get strong GAC3
		  addCost(x,y,z,valx,valy,*iterZ,-cost);
		}
	  }
	// }
    xy->addcost(x,y,valx,valy,cost);
#ifdef DEECOMPLETE
    getVar(0)->queueDEE();
    getVar(1)->queueDEE();
    getVar(2)->queueDEE();
#endif
}

template <typename T1, typename T2, typename T3>
void TernaryConstraint::extend(T1 getCost, T2 addCost, bool functionalZ, T3 getFunctionZ, BinaryConstraint* xy, EnumeratedVariable *x, EnumeratedVariable *y, EnumeratedVariable *z, Value valx, Value valy, Cost cost) {
	assert(ToulBar2::verbose < 4 || ((cout << "extend(C" << getVar(0)->getName() << "," << getVar(1)->getName() << "," << getVar(2)->getName() << ", ((" << x->getName() << "," << valx << "),(" << y->getName() << "," << valy << ")), " << cost << ")" << endl), true));
	assert(cost >= MIN_COST);
	// BUG!
	// if (functionalZ) {
	//   Value valz = getFunctionZ(x,y, valx, valy);
	//   if (valz != WRONG_VAL && z->canbe(valz)) {
    //     addCost(x,y,z,valx,valy,valz,cost);
	//   }
 	// } else {
	  for (EnumeratedVariable::iterator iterZ = z->begin(); iterZ != z->end(); ++iterZ) {
        addCost(x,y,z,valx,valy,*iterZ,cost);
	  }
	// }
    assert(xy->connected());
    xy->addcost(x,y,valx,valy,-cost);
}

template <typename T1, typename T2, typename T3>
void TernaryConstraint::findSupport(T1 getCost, bool functionalY, T2 getFunctionY, bool functionalZ, T3 getFunctionZ,
            EnumeratedVariable *x, EnumeratedVariable *y, EnumeratedVariable *z,
            int getIndexX, int getIndexY, int getIndexZ,
            vector< pair<Value,Value> > &supportX, vector<StoreCost> &deltaCostsX, 
			vector< pair<Value,Value> > &supportY, vector< pair<Value,Value> > &supportZ)
{
    assert(getIndex(y) < getIndex(z));  // check that support.first/.second is consistent with y/z parameters
    assert(connected());
    wcsp->revise(this);
    if (ToulBar2::verbose >= 3) cout << "findSupport C" << x->getName() << "?" << "," << y->getName() << ((functionalY)?"!":"") << "," << z->getName() << ((functionalZ)?"!":"") << endl;
    if (ToulBar2::verbose >= 7) cout << *x << endl << *y << endl << *z << endl << *this;
    bool supportBroken = false;
    for (EnumeratedVariable::iterator iterX = x->begin(); iterX != x->end(); ++iterX) {
    	unsigned int xindex = x->toIndex(*iterX);
        pair<Value,Value> support = supportX[xindex];
        if (y->cannotbe(support.first) || z->cannotbe(support.second) || 
            getCost(x,y,z,*iterX,support.first,support.second) > MIN_COST) {
        	support = make_pair(y->getInf(),z->getInf());
            Cost minCost = MAX_COST;
			if (functionalZ) {
			  for (EnumeratedVariable::iterator iterY = y->begin(); minCost > MIN_COST && iterY != y->end(); ++iterY) {
				Value valZ = getFunctionZ(x,y, *iterX, *iterY);
				if (valZ != WRONG_VAL && z->canbe(valZ)) {
					Cost cost = getCost(x,y,z, *iterX,*iterY,valZ);
                    if (GLB(&minCost, cost)) {
                        support = make_pair(*iterY,valZ);
                    }
				}
			  }
			} else if (functionalY) {
			  for (EnumeratedVariable::iterator iterZ = z->begin(); minCost > MIN_COST && iterZ != z->end(); ++iterZ) {
				Value valY = getFunctionY(x,z, *iterX, *iterZ);
				if (valY != WRONG_VAL && y->canbe(valY)) {
					Cost cost = getCost(x,y,z, *iterX,valY,*iterZ);
                    if (GLB(&minCost, cost)) {
					  support = make_pair(valY,*iterZ);
                    }
				}
			  }
			} else {
			  for (EnumeratedVariable::iterator iterY = y->begin(); minCost > MIN_COST && iterY != y->end(); ++iterY) {
                for (EnumeratedVariable::iterator iterZ = z->begin(); minCost > MIN_COST && iterZ != z->end(); ++iterZ) {
                	Cost cost = getCost(x,y,z, *iterX,*iterY,*iterZ);
                    if (GLB(&minCost, cost)) {
                        support = make_pair(*iterY,*iterZ);
                    }
                }
			  }
            }
            if (minCost > MIN_COST) {
                supportBroken |= project(x, *iterX, minCost, deltaCostsX);
                if (deconnected()) return;
            }
            supportX[xindex] = support;
            assert(getIndexY < getIndexZ);
            // warning! do not break DAC support for the variable in the constraint scope having the smallest wcspIndex
            if (getIndexY!=getDACScopeIndex()) {
              unsigned int yindex = y->toIndex(support.first);
			  if (getIndexX < getIndexZ) supportY[yindex] = make_pair(*iterX, support.second);
			  else supportY[yindex] = make_pair(support.second, *iterX);
            }
            if (getIndexZ!=getDACScopeIndex()) {
              unsigned int zindex = z->toIndex(support.second);
			  if (getIndexX < getIndexY) supportZ[zindex] = make_pair(*iterX, support.first);
			  else supportZ[zindex] = make_pair(support.first, *iterX);
            }
        }
    }
    if (supportBroken) {
        x->findSupport();
    }
}

// take into account associated binary constraints and perform unary extension to binary instead of ternary
template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10> void TernaryConstraint::findFullSupport(T1 getCost, T2 getCostXZY, T3 getCostYZX, T4 getCostWithBinaries, T5 addCost, T6 addCostXZY, T7 addCostYZX, bool functionalX, T8 getFunctionX, bool functionalY, T9 getFunctionY, bool functionalZ, T10 getFunctionZ,
            EnumeratedVariable *x, EnumeratedVariable *y, EnumeratedVariable *z,
            int getIndexX, int getIndexY, int getIndexZ,
            vector< pair<Value,Value> > &supportX, vector<StoreCost> &deltaCostsX, 
            vector< pair<Value,Value> > &supportY, vector<StoreCost> &deltaCostsY,
            vector< pair<Value,Value> > &supportZ, vector<StoreCost> &deltaCostsZ,
			BinaryConstraint* xy, BinaryConstraint* xz, BinaryConstraint* yz)
{
    assert(connected());
    wcsp->revise(this);   
    if (ToulBar2::verbose >= 3) cout << "findFullSupport C"<< x->getName() << ((functionalX)?"!":"") << "," << y->getName() << ((functionalY)?"!":"") << "," << z->getName() << ((functionalZ)?"!":"") << endl;
    if (ToulBar2::verbose >= 7) cout << *x << endl << *y << endl << *z << endl << *this << *xy << *xz << *yz;
    bool supportBroken = false;
    bool supportReversed = (getIndexY > getIndexZ);
    for (EnumeratedVariable::iterator iterX = x->begin(); iterX != x->end(); ++iterX) {
        unsigned int xindex = x->toIndex(*iterX);
        pair<Value,Value> support = (supportReversed)?make_pair(supportX[xindex].second,supportX[xindex].first):make_pair(supportX[xindex].first,supportX[xindex].second);
        if (y->cannotbe(support.first) || z->cannotbe(support.second) || 
            getCostWithBinaries(x,y,z, *iterX,support.first,support.second) + y->getCost(support.first) + z->getCost(support.second) > MIN_COST) {
        	support = make_pair(y->getInf(),z->getInf());
            Cost minCost = MAX_COST;
			if (functionalZ) {
			  for (EnumeratedVariable::iterator iterY = y->begin(); minCost > MIN_COST && iterY != y->end(); ++iterY) {
				Value valZ = getFunctionZ(x,y, *iterX, *iterY);
				if (valZ != WRONG_VAL && z->canbe(valZ)) {
				  Cost cost = getCostWithBinaries(x,y,z, *iterX,*iterY,valZ) + y->getCost(*iterY) + z->getCost(valZ);
				  if (GLB(&minCost, cost)) {
					support = make_pair(*iterY,valZ);
				  }
				}
			  }
			} else if (functionalY) {
			  for (EnumeratedVariable::iterator iterZ = z->begin(); minCost > MIN_COST && iterZ != z->end(); ++iterZ) {
				Value valY = getFunctionY(x,z, *iterX, *iterZ);
				if (valY != WRONG_VAL && y->canbe(valY)) {
				  Cost cost = getCostWithBinaries(x,y,z, *iterX,valY,*iterZ) + y->getCost(valY) + z->getCost(*iterZ);
				  if (GLB(&minCost, cost)) {
					support = make_pair(valY,*iterZ);
				  }
				}
			  }
			} else {
			  for (EnumeratedVariable::iterator iterY = y->begin(); minCost > MIN_COST && iterY != y->end(); ++iterY) {
                for (EnumeratedVariable::iterator iterZ = z->begin(); minCost > MIN_COST && iterZ != z->end(); ++iterZ) {
				  Cost cost = getCostWithBinaries(x,y,z, *iterX,*iterY,*iterZ) + y->getCost(*iterY) + z->getCost(*iterZ);
				  if (GLB(&minCost, cost)) {
					support = make_pair(*iterY,*iterZ);
				  }
                }
			  }
			}
			if (CUT(minCost +  wcsp->getLb(), wcsp->getUb())) {
			  supportBroken |= project(x, *iterX, minCost, deltaCostsX);
			  if (deconnected()) return;
			  continue;
			}
			assert(minCost < MAX_COST);

            if (minCost > MIN_COST) {
                // extend unary to binary
                for (EnumeratedVariable::iterator iterY = y->begin(); iterY != y->end(); ++iterY) {
                    if (y->getCost(*iterY) > MIN_COST) {
                        Cost costfromy = MIN_COST;
                        for (EnumeratedVariable::iterator iterZ = z->begin(); iterZ != z->end(); ++iterZ) {
						  Cost cost = getCost(x,y,z, *iterX,*iterY,*iterZ);
                            if (LUBTEST(cost, minCost)) {
                                Cost zcost = z->getCost(*iterZ);
                                Cost xycost = (xy->connected())?xy->getCost(x,y,*iterX,*iterY):MIN_COST;
                                Cost xzcost = (xz->connected())?xz->getCost(x,z,*iterX,*iterZ):MIN_COST;
                                Cost yzcost = (yz->connected())?yz->getCost(y,z,*iterY,*iterZ):MIN_COST;
                                Cost remain = minCost - (cost + xycost + xzcost + yzcost + zcost);
                                LUB(&costfromy, remain);
                            }
                        }
                        assert(costfromy <= y->getCost(*iterY));
                        if (costfromy > MIN_COST) {
						    assert(x->unassigned() && y->unassigned());
                            xy->reconnect();    // must be done before using the constraint
                            xy->extend(xy->getIndex(y), *iterY, costfromy);
                        }
                    }
                }
                for (EnumeratedVariable::iterator iterZ = z->begin(); iterZ != z->end(); ++iterZ) {
                    if (z->getCost(*iterZ) > MIN_COST) {
                        Cost costfromz = MIN_COST;
                        for (EnumeratedVariable::iterator iterY = y->begin(); iterY != y->end(); ++iterY) {
						  Cost cost = getCost(x,y,z, *iterX,*iterY,*iterZ);
                            if (LUBTEST(cost, minCost)) {
                                Cost xycost = (xy->connected())?xy->getCost(x,y,*iterX,*iterY):MIN_COST;
                                Cost xzcost = (xz->connected())?xz->getCost(x,z,*iterX,*iterZ):MIN_COST;
                                Cost yzcost = (yz->connected())?yz->getCost(y,z,*iterY,*iterZ):MIN_COST;
                                Cost remain = minCost - (cost + xycost + xzcost + yzcost);
                                LUB(&costfromz, remain);
                            }
                        }
                        assert(costfromz <= z->getCost(*iterZ));
                        if (costfromz > MIN_COST) {
                            assert(x->unassigned() && z->unassigned());
						    xz->reconnect();    // must be done before using the constraint
                            xz->extend(xz->getIndex(z), *iterZ, costfromz);
                        }
                    }
                }
                // extend binary to ternary
                if (yz->connected()) {
                    for (EnumeratedVariable::iterator iterY = y->begin(); iterY != y->end(); ++iterY) {
                        Cost xycost = (xy->connected())?xy->getCost(x,y,*iterX,*iterY):MIN_COST;
                        for (EnumeratedVariable::iterator iterZ = z->begin(); iterZ != z->end(); ++iterZ) {
                            Cost yzcost = yz->getCost(y,z,*iterY,*iterZ);
                            if (yzcost > MIN_COST) {
							  Cost cost = getCost(x,y,z, *iterX,*iterY,*iterZ);
                                if (LUBTEST(cost, minCost)) {
                                    Cost xzcost = (xz->connected())?xz->getCost(x,z,*iterX,*iterZ):MIN_COST;
                                    Cost remain = minCost - (cost + xycost + xzcost);
                                    if (remain > MIN_COST) extend(getCostYZX,addCostYZX,functionalX,getFunctionX,yz,y,z,x,*iterY,*iterZ, remain);
                                }
                            }
                        }
                    }
                }
                if (xy->connected()) {
                    for (EnumeratedVariable::iterator iterY = y->begin(); iterY != y->end(); ++iterY) {
                        Cost xycost = xy->getCost(x,y,*iterX,*iterY);
                        if (xycost > MIN_COST) {
                            Cost costfromxy = MIN_COST;
                            for (EnumeratedVariable::iterator iterZ = z->begin(); iterZ != z->end(); ++iterZ) {
							  Cost cost = getCost(x,y,z, *iterX,*iterY,*iterZ);
                                if (LUBTEST(cost, minCost)) {
                                    Cost xzcost = (xz->connected())?xz->getCost(x,z,*iterX,*iterZ):MIN_COST;
//                                    Cost yzcost = (yz->connected())?yz->getCost(y,z,*iterY,*iterZ):MIN_COST;
//                                    Cost remain = minCost - (cost + xzcost + yzcost);
                                    Cost remain = minCost - (cost + xzcost);
                                    LUB(&costfromxy, remain);
                                }
                            }
                            assert(costfromxy <= xycost);
                            if (costfromxy > MIN_COST) {
							  extend(getCost,addCost,functionalZ,getFunctionZ,xy,x,y,z,*iterX,*iterY,costfromxy);
                            }
                        }
                    }
                }
                if (xz->connected()) {
                    for (EnumeratedVariable::iterator iterZ = z->begin(); iterZ != z->end(); ++iterZ) {
                        Cost xzcost = xz->getCost(x,z,*iterX,*iterZ);
                        if (xzcost > MIN_COST) {
                            Cost costfromxz = MIN_COST;
                            for (EnumeratedVariable::iterator iterY = y->begin(); iterY != y->end(); ++iterY) {
							  Cost cost = getCost(x,y,z, *iterX,*iterY,*iterZ);
                                if (LUBTEST(cost, minCost)) {
//                                    Cost yzcost = (yz->connected())?yz->getCost(y,z,*iterY,*iterZ):MIN_COST;
//                                    Cost remain = minCost - (cost + yzcost);
                                    Cost remain = minCost - cost;
                                    LUB(&costfromxz, remain);
                                }
                            }
                            assert(costfromxz <= xzcost);
                            if (costfromxz > MIN_COST) {
							  extend(getCostXZY,addCostXZY,functionalY,getFunctionY,xz,x,z,y,*iterX,*iterZ,costfromxz);
                            }
                        }
                    }
                }
//                if (yz->connected()) {
//                    for (EnumeratedVariable::iterator iterY = y->begin(); iterY != y->end(); ++iterY) {
//                        for (EnumeratedVariable::iterator iterZ = z->begin(); iterZ != z->end(); ++iterZ) {
//                            Cost yzcost = yz->getCost(y,z,*iterY,*iterZ);
//                            if (yzcost > MIN_COST) {
//                                Cost cost = getCost(x,y,z,*iterX,*iterY,*iterZ);
//                                if (LUBTEST(cost, minCost)) {
//                                    assert(yzcost >= minCost - cost);
//                                    extend(getCostYZX,addCostYZX,functionalX,getFunctionX,yz,y,z,x,*iterY,*iterZ,minCost - cost);
//                                }
//                            }
//                        }
//                    }
//                }
                supportBroken |= project(x, *iterX, minCost, deltaCostsX);
                if (deconnected()) return;
            }
            
            supportX[xindex] = (supportReversed)?make_pair(support.second, support.first):make_pair(support.first, support.second);
            
            unsigned int yindex = y->toIndex(support.first);
            unsigned int zindex = z->toIndex(support.second);
            if (getIndexX < getIndexZ) supportY[yindex] = make_pair(*iterX, support.second);
            else supportY[yindex] = make_pair(support.second, *iterX);
            if (getIndexX < getIndexY) supportZ[zindex] = make_pair(*iterX, support.first);
            else supportZ[zindex] = make_pair(support.first, *iterX);
        }
    }
    if (supportBroken) {
        x->findSupport();
    }
}

pair< pair<Cost,Cost>, pair<Cost,Cost> > TernaryConstraint::getMaxCost(int varIndex, Value a, Value b)
{
	Cost maxcosta = MIN_COST;
	Cost diffcosta = MIN_COST;
	Cost maxcostb = MIN_COST;
	Cost diffcostb = MIN_COST;
	if (varIndex == 0) {
		Cost ucosta = x->getCost(a);
		Cost ucostb = x->getCost(b);
		for (EnumeratedVariable::iterator iterY = y->begin(); iterY != y->end(); ++iterY) {
			Cost ucosty = y->getCost(*iterY);
			for (EnumeratedVariable::iterator iterZ = z->begin(); iterZ != z->end(); ++iterZ) {
				Cost costa = getCost(a, *iterY, *iterZ);
				Cost costb = getCost(b, *iterY, *iterZ);
				if (costa > maxcosta) maxcosta = costa;
				if (costb > maxcostb) maxcostb = costb;
    			Cost ucostz = z->getCost(*iterZ);
				if (!CUT(ucostb + getCostWithBinaries(b, *iterY, *iterZ) + ucosty + ucostz + wcsp->getLb(), wcsp->getUb())) {
					if (costa-costb > diffcosta) diffcosta = costa-costb;
				}
				if (!CUT(ucosta + getCostWithBinaries(a, *iterY, *iterZ) + ucosty + ucostz + wcsp->getLb(), wcsp->getUb())) {
					if (costb-costa > diffcostb) diffcostb = costb-costa;
				}
			}
		}
	} else if (varIndex == 1) {
		Cost ucosta = y->getCost(a);
		Cost ucostb = y->getCost(b);
		for (EnumeratedVariable::iterator iterX = x->begin(); iterX != x->end(); ++iterX) {
			Cost ucostx = x->getCost(*iterX);
			for (EnumeratedVariable::iterator iterZ = z->begin(); iterZ != z->end(); ++iterZ) {
				Cost costa = getCost(*iterX, a, *iterZ);
				Cost costb = getCost(*iterX, b, *iterZ);
				if (costa > maxcosta) maxcosta = costa;
				if (costb > maxcostb) maxcostb = costb;
    			Cost ucostz = z->getCost(*iterZ);
				if (!CUT(ucostb + getCostWithBinaries(*iterX, b, *iterZ) + ucostx + ucostz + wcsp->getLb(), wcsp->getUb())) {
					if (costa-costb > diffcosta) diffcosta = costa-costb;
				}
				if (!CUT(ucosta + getCostWithBinaries(*iterX, a, *iterZ) + ucostx + ucostz + wcsp->getLb(), wcsp->getUb())) {
					if (costb-costa > diffcostb) diffcostb = costb-costa;
				}
			}
		}
	} else {
		assert(varIndex == 2);
		Cost ucosta = z->getCost(a);
		Cost ucostb = z->getCost(b);
		for (EnumeratedVariable::iterator iterX = x->begin(); iterX != x->end(); ++iterX) {
			Cost ucostx = x->getCost(*iterX);
    		for (EnumeratedVariable::iterator iterY = y->begin(); iterY != y->end(); ++iterY) {
				Cost costa = getCost(*iterX, *iterY, a);
				Cost costb = getCost(*iterX, *iterY, b);
				if (costa > maxcosta) maxcosta = costa;
				if (costb > maxcostb) maxcostb = costb;
    			Cost ucosty = y->getCost(*iterY);
				if (!CUT(ucostb + getCostWithBinaries(*iterX, *iterY, b) + ucostx + ucosty + wcsp->getLb(), wcsp->getUb())) {
					if (costa-costb > diffcosta) diffcosta = costa-costb;
				}
				if (!CUT(ucosta + getCostWithBinaries(*iterX, *iterY, a) + ucostx + ucosty + wcsp->getLb(), wcsp->getUb())) {
					if (costb-costa > diffcostb) diffcostb = costb-costa;
				}
			}
		}
	}
	assert(maxcosta >= diffcosta);
	assert(maxcostb >= diffcostb);
	return make_pair(make_pair(maxcosta,diffcosta), make_pair(maxcostb,diffcostb));
}

bool TernaryConstraint::separability(EnumeratedVariable* vy, EnumeratedVariable* vz)
{
    	Cost c1,c;
    	Char tch[4];
    	bool neweq = true; // true if we have  not a difference value
    	bool sep = true; // false if vy and vz are not separable
    	Cost diff = 0;
    	first(vy,vz);
    	EnumeratedVariable::iterator itvyfirst = yvar->begin();
    	if(ToulBar2::verbose >= 1) cout << " [ " << zvar->getName()  << " " << xvar->getName() << " " << yvar->getName() << " ] ?"; // << endl;
    	while (sep && itvyfirst != yvar->end()){
    		itvx = xvar->begin();
    		itvy = yvar->begin();
    		itvz = zvar->begin();
    		while(sep && itvx != xvar->end()) {
    			unsigned int ix = xvar->toIndex(*itvx);
    			tch[0] = ix + CHAR_FIRST;
    			while(sep && itvy != yvar->end()) {
					unsigned int iy = yvar->toIndex(*itvy);
    				tch[1] = iy + CHAR_FIRST;
    				while(sep && itvy != itvyfirst && itvz != zvar->end()) {
    					unsigned int iz = zvar->toIndex(*itvz);
    					tch[2] = iz + CHAR_FIRST;
    					tch[3] = '\0';

    					c1 = getCost(xvar,yvar,zvar,*itvx, *(itvyfirst), *itvz);
    					c = getCost(xvar,yvar,zvar,*itvx, *itvy, *itvz);

    					if(!universe(c, c1, wcsp->getUb())){
    						if(ToulBar2::verbose >= 3) {if(!neweq)  cout << "= \n"; else cout << endl;}
    						if(ToulBar2::verbose >= 3) cout << " C" << tch[2]-CHAR_FIRST << "."  << tch[0]-CHAR_FIRST << "." << tch[1]-CHAR_FIRST << " -  C" << tch[2]-CHAR_FIRST << "." << tch[0]-48  << "."  << yvar->toIndex(*itvyfirst) << " = " << c << " - " << c1;

    						if(neweq) {diff = squareminus(c,c1,wcsp->getUb()); neweq = false; }
    						else sep = (diff == squareminus(c,c1,wcsp->getUb()));
    						if(ToulBar2::verbose >= 3) cout << " = " << squareminus(c,c1,wcsp->getUb()) <<  endl;
    					}
    					else{
    						if(ToulBar2::verbose >= 3) cout << "universe\n";
    					}
    					++itvz;
    				}
    				++itvy;
    				itvz = zvar->begin();
    				neweq = true;


    			}
    			++itvx;
    			itvy = yvar->begin();
    			neweq = true;
    		}
    		++itvyfirst;
    		itvx = xvar->begin();
    		neweq = true;
    		if(ToulBar2::verbose >= 3) cout << "---\n";
    	}
    	return sep;
}

void TernaryConstraint::separate(EnumeratedVariable *vy, EnumeratedVariable *vz)
{
	Cost cost,minCost = MAX_COST;
	//assert(separability(vy,vz));
	first(vy,vz);
    vector<Cost> costsZX(zvar->getDomainInitSize() * xvar->getDomainInitSize(), MIN_COST);
    vector<Cost> costsXY(xvar->getDomainInitSize() * yvar->getDomainInitSize(), MIN_COST);
	string xv(xvar->getName()), yv(yvar->getName()), zv(zvar->getName());
	if(ToulBar2::verbose == 1) cout << "\n";

	if(ToulBar2::verbose >= 3) cout << "[ " << zvar->getName() << " " << xvar->getName() << " ]" << endl;
	while(itvz != zvar->end()) {
		while(itvx != xvar->end()) {
			minCost = MAX_COST;
			while(itvy != yvar->end()) {
				cost = getCost(xvar,yvar,zvar,*itvx, *itvy, *itvz);
				if(cost < minCost) minCost = cost;
				if(minCost>= wcsp->getUb()) minCost = wcsp->getUb();
				++itvy;
			}
			if(ToulBar2::verbose >= 3) cout << *itvx << " " << *itvz << " : " << minCost << endl;
			costsZX[zvar->toIndex(*itvz)*xvar->getDomainInitSize()+xvar->toIndex(*itvx)] = minCost;
			++itvx;
			itvy = yvar->begin();
		}
		++itvz;
		itvx = xvar->begin();
	}
	BinaryConstraint* existZX = xvar->getConstr(zvar);
	assert(existZX);
	BinaryConstraint  *zx = new BinaryConstraint(wcsp,zvar, xvar,costsZX, &wcsp->getStore()->storeCost);
	if(ToulBar2::verbose >= 3)cout << "-------------\n";
	if(ToulBar2::verbose >= 3) cout << "[ " << xvar->getName() << " " << yvar->getName() << " ]" <<  endl;

	first(vy,vz);
	Cost costzx;
	while(itvx != xvar->end()) {
		while(itvy != yvar->end()) {
			itvz = zvar->begin();
			do {
				cost = getCost(xvar,yvar,zvar,*itvx, *itvy, *itvz);
				costzx = zx->getCost(*itvz,*itvx);
				++itvz;
			}while(itvz != zvar->end() && cost>= wcsp->getUb() && costzx >= wcsp->getUb() );
			costsXY[xvar->toIndex(*itvx)*yvar->getDomainInitSize()+yvar->toIndex(*itvy)] = squareminus(cost,costzx,wcsp->getUb());
			if(ToulBar2::verbose >= 3) cout << *itvx << " " << *itvy << " : " << squareminus(cost,costzx,wcsp->getUb()) << endl;
			++itvy;
		}
		++itvx;
		itvy = yvar->begin();
	}
	BinaryConstraint* existXY = xvar->getConstr(yvar);
	assert(existXY);
	BinaryConstraint * xy = new BinaryConstraint(wcsp,xvar, yvar,costsXY, &wcsp->getStore()->storeCost);

	assert(verifySeparate(zx,xy));

	// fusion with the existing constraint (xz)
	if(!zx->universal()){
		if(ToulBar2::verbose >= 1) cout << "[ " << zv << " " << xv << " ]" << endl;
		existZX->addCosts(zx);
		existZX->reconnect();
		existZX->propagate();
	}
	zx->deconnect();	//  unsafe to delete zx due to x and z lists of cost functions

	// fusion with the existing constraint (xy)
	if(!xy->universal()){
		if(ToulBar2::verbose >= 1) cout << "[ " << xv << " " << yv << " ]" <<  endl;
		existXY->addCosts(xy);
		existXY->reconnect();
		existXY->propagate();
	}
	xy->deconnect();	//  unsafe to delete xy due to x and y lists of cost functions
	deconnect();
}

template <typename T1, typename T2, typename T3, typename T4>
void TernaryConstraint::projectTernaryBinary(T1 getCost, T2 getCostYZX, T3 addCostYZX, bool functionalX, T4 getFunctionX,
		 EnumeratedVariable* x, EnumeratedVariable* y, EnumeratedVariable* z, BinaryConstraint* yzin)
{
	//	cout << "PROJECT " << *this <<endl;
	//	cout << "on " << y->getName() << "," << z->getName() << endl;

    bool flag = false;
    for (EnumeratedVariable::iterator itery = y->begin(); itery != y->end(); ++itery) {
    for (EnumeratedVariable::iterator iterz = z->begin(); iterz != z->end(); ++iterz) {
        Cost mincost = MAX_COST;
		for (EnumeratedVariable::iterator iterx = x->begin(); iterx != x->end(); ++iterx) {
		  Cost curcost = getCost(x,y,z, *iterx, *itery, *iterz);                            
		  GLB(&mincost, curcost);
		}
        if (mincost > MIN_COST) {
            flag = true;
            project(getCostYZX,addCostYZX,functionalX,getFunctionX,yzin,y,z,x,*itery,*iterz,mincost);
        }
    }}
 
	if(wcsp->getTreeDec()) {
		yzin->setCluster( getCluster() );
	}
    
    if (flag) {
        if (y->unassigned() && z->unassigned()) yzin->reconnect();
        yzin->propagate();
    }
}

template <typename T1, typename T2, typename T3>
bool TernaryConstraint::isEAC(T1 getCostWithBinaries, bool functionalY, T2 getFunctionY, bool functionalZ, T3 getFunctionZ,
            EnumeratedVariable *x, Value a, EnumeratedVariable *y, EnumeratedVariable *z,
            vector< pair<Value,Value> > &supportX)
{
    unsigned int xindex = x->toIndex(a);
    assert(getIndex(y) < getIndex(z));
    pair<Value,Value> support = supportX[xindex];
    if (y->cannotbe(support.first) || z->cannotbe(support.second) || 
        getCostWithBinaries(x,y,z,a,support.first,support.second) + y->getCost(support.first) + z->getCost(support.second) > MIN_COST) {
	  if (functionalZ) {
        for (EnumeratedVariable::iterator iterY = y->begin(); iterY != y->end(); ++iterY) {
		  Value valZ = getFunctionZ(x,y, a, *iterY);
		  if (valZ != WRONG_VAL && z->canbe(valZ)) {
			if (getCostWithBinaries(x,y,z,a,*iterY,valZ) + y->getCost(*iterY) + z->getCost(valZ) == MIN_COST) {
			  supportX[xindex] = make_pair(*iterY,valZ);
			  return true;
			}
		  }
		}
	  } else if (functionalY) {
		for (EnumeratedVariable::iterator iterZ = z->begin(); iterZ != z->end(); ++iterZ) {
		  Value valY = getFunctionY(x,z, a, *iterZ);
		  if (valY != WRONG_VAL && y->canbe(valY)) {
			if (getCostWithBinaries(x,y,z,a,valY,*iterZ) + y->getCost(valY) + z->getCost(*iterZ) == MIN_COST) {
			  supportX[xindex] = make_pair(valY,*iterZ);
			  return true;
			}
		  }
		}
	  } else {
        for (EnumeratedVariable::iterator iterY = y->begin(); iterY != y->end(); ++iterY) {
            for (EnumeratedVariable::iterator iterZ = z->begin(); iterZ != z->end(); ++iterZ) {
                if (getCostWithBinaries(x,y,z,a,*iterY,*iterZ) + y->getCost(*iterY) + z->getCost(*iterZ) == MIN_COST) {
                    supportX[xindex] = make_pair(*iterY,*iterZ);
                    return true;
                }
            }
        }
	  }
//        cout << x->getName() << " = " << a << " not EAC due to constraint " << *this << endl; 
	  return false;
    }
    return true;
}


void TernaryConstraint::fillxy() {
	TreeDecomposition* td = wcsp->getTreeDec();
    BinaryConstraint* xy_ = NULL;
    xy_ = x->getConstr(y); 
    if(td && xy_ && (getCluster() != xy_->getCluster())) {  
    	BinaryConstraint* xy__ =  x->getConstr(y, getCluster());
    	if(xy__) xy_ = xy__; // we have found another constraint of the same cluster
    }
    if(!xy_ || (xy_ && td && (getCluster() != xy_->getCluster())) ) {
	    xy = wcsp->newBinaryConstr(x,y,this); 
		xy->setCluster( getCluster() );
		if(td && xy_ && (getCluster() != xy_->getCluster())) xy->setDuplicate();
		wcsp->elimBinOrderInc(); 
    } else xy = xy_; 
	if(xy->isDuplicate()) setDuplicate();
}

void TernaryConstraint::fillxz() {
	TreeDecomposition* td = wcsp->getTreeDec();
    BinaryConstraint* xz_ = NULL;
    xz_ = x->getConstr(z); 
    if(td && xz_ && (getCluster() != xz_->getCluster())) {
    	BinaryConstraint* xz__ =  x->getConstr(z, getCluster());
    	if(xz__) xz_ = xz__; // we have found another constraint of the same cluster
    }
    if(!xz_|| (xz_ && td && getCluster() != xz_->getCluster()) ) {
	    xz = wcsp->newBinaryConstr(x,z,this); 
		xy->setCluster( getCluster() );
		if(td && xz_ && (getCluster() != xz_->getCluster())) xz->setDuplicate();
		wcsp->elimBinOrderInc(); 
		if(td) xz->setCluster( getCluster() );
    } else xz = xz_; 
	if(xz->isDuplicate()) setDuplicate();
}

void TernaryConstraint::fillyz() {
	TreeDecomposition* td = wcsp->getTreeDec();
    BinaryConstraint* yz_ = NULL;
    yz_ = y->getConstr(z); 
    if(td && yz_ && (getCluster() != yz_->getCluster())) {
    	BinaryConstraint* yz__ =  y->getConstr(z, getCluster());
    	if(yz__) yz_ = yz__; 
    }
    if(!yz_ || (yz_ && td && getCluster() != yz_->getCluster()) ) {
	    yz = wcsp->newBinaryConstr(y,z,this); 
		yz->setCluster( getCluster() );
		if(td && yz_ && (getCluster() != yz_->getCluster())) yz->setDuplicate();
		wcsp->elimBinOrderInc(); 
    } else yz = yz_; 
	if(yz->isDuplicate()) setDuplicate();
}

void TernaryConstraint::fillElimConstrBinaries()
{
	fillxy();
	fillxz();
	fillyz();
	
    if (ToulBar2::verbose > 1) cout << " fillElimConstrBinaries (" << x->wcspIndex << "," << y->wcspIndex << "," << z->wcspIndex << ")  ";
	
	
}



void TernaryConstraint::setDuplicates() {
	assert(wcsp->getTreeDec());
	if(xy->getCluster() != cluster) {
		BinaryConstraint* xy_ =  x->getConstr(y, getCluster());
		if(xy_) {
			if(xy_->isDuplicate()) setDuplicate();
			xy = xy_;
		} else { 
			wcsp->initElimConstr();
			xy = wcsp->newBinaryConstr(x,y); 
			xy->setCluster( cluster );
			xy->setDuplicate();
			wcsp->elimBinOrderInc(); 
			setDuplicate();
		}
	}
	if(xz->getCluster() != cluster) {
		BinaryConstraint* xz_ =  x->getConstr(z, getCluster());
		if(xz_) {
			xz = xz_;
			if(xz_->isDuplicate()) setDuplicate();
		} else { 
			wcsp->initElimConstr();
			xz = wcsp->newBinaryConstr(x,z); 
			xz->setCluster( getCluster() );
			xz->setDuplicate();
			wcsp->elimBinOrderInc(); 
			setDuplicate();
		}
	}
	if(yz->getCluster() != cluster) {
		BinaryConstraint* yz_ =  y->getConstr(z, getCluster());
		if(yz_) {
			yz = yz_;
			if(yz_->isDuplicate()) setDuplicate();
		} else { 
			wcsp->initElimConstr();
			yz = wcsp->newBinaryConstr(y,z); 
			yz->setCluster( getCluster() );
			yz->setDuplicate();
			wcsp->elimBinOrderInc(); 
			setDuplicate();
		}
	}
}


bool TernaryConstraint::verify(EnumeratedVariable *x, EnumeratedVariable *y, EnumeratedVariable *z)
{
    for (EnumeratedVariable::iterator iterX = x->begin(); iterX != x->end(); ++iterX) {
        Cost minCost = MAX_COST;
        for (EnumeratedVariable::iterator iterY = y->begin(); minCost > MIN_COST && iterY != y->end(); ++iterY) {
            for (EnumeratedVariable::iterator iterZ = z->begin(); minCost > MIN_COST && iterZ != z->end(); ++iterZ) {
                Cost cost = getCost(x,y,z,*iterX,*iterY,*iterZ);
                if (ToulBar2::LcLevel>=LC_DAC && getIndex(x)==getDACScopeIndex()) cost += y->getCost(*iterY) + z->getCost(*iterZ);
                GLB(&minCost, cost);
            }
        }
        if (minCost > MIN_COST) {
            cout << "not FDAC: variable " << x->getName() << " value " << *iterX << " of " << *this;
            return false;
        }
    }
    return true;
}





bool TernaryConstraint::verify() {
  TreeDecomposition* td = wcsp->getTreeDec();
  
  if(td) {
	  if ( cluster != xy->getCluster() || cluster != xz->getCluster() ||  cluster != yz->getCluster())  return false;
  }

  if (ToulBar2::LcLevel==LC_DAC) {
    switch(getDACScopeIndex()) {
        case 0: return verifyX(); break;
        case 1: return verifyY(); break;
        case 2: return verifyZ(); break;
	    default: return false;
    }
  } else {
	return verifyX() && verifyY() && verifyZ();
  }
}
