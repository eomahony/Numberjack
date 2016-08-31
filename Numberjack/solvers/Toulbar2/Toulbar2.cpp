
/** \file Toulbar2.cpp
    \brief Solver interface for PYTHON Wrapper.
*/

//#define _DEBUGWRAP 1

#include "Toulbar2.hpp"

static Toulbar2Solver* MyWCSPSolver = NULL;

void timeout()
{
  assert(MyWCSPSolver);
  MyWCSPSolver->interrupted = true;
  double time = cpuTime() - ToulBar2::startCpuTime;
  if (ToulBar2::verbose>=0) {
	if (MyWCSPSolver->wcsp->getUb() < MyWCSPSolver->upperbound) {
	  cout << "Best upper bound: " << MyWCSPSolver->costshift + MyWCSPSolver->wcsp->getUb() << " in " << MyWCSPSolver->solver->getNbBacktracks() << " backtracks and " << MyWCSPSolver->solver->getNbNodes() << " nodes and " << time << " seconds." << endl;
	} else {
	  cout << "No solution found after " << MyWCSPSolver->solver->getNbBacktracks() << " backtracks and " << MyWCSPSolver->solver->getNbNodes() << " nodes and " << time << " seconds." << endl;
	}
  }
  ToulBar2::interrupted = true; // throw TimeOut();
}

/**************************************************************
 ********************     EXPRESSION        *******************
 **************************************************************/

Toulbar2_Expression::Toulbar2_Expression()
{
    _solver = NULL;
    _iinf = 0;
    _isup = 1; 
	_size = 2;
	_domain = NULL;
	_wcspIndex = -1;
}

Toulbar2_Expression::Toulbar2_Expression(const int nval)
{
    _solver = NULL;
    _iinf = 0;
    _isup = (nval-1); 
	_size = nval;
	_domain = NULL;
	_wcspIndex = -1;
}

Toulbar2_Expression::Toulbar2_Expression(const int lb, const int ub)
{
    _solver = NULL;
    _iinf = lb;
    _isup = ub;
	_size = ub - lb + 1;
	_domain = NULL;
	_wcspIndex = -1;
}

Toulbar2_Expression::Toulbar2_Expression(Toulbar2IntArray& vals)
{
    _solver = NULL;
	Value min = vals.get_item(0);
	Value max = vals.get_item(0);
	_size = vals.size();
	_domain = new Value[vals.size()];
	for (int i=0; i<vals.size(); i++) {
	  Value v = vals.get_item(i);
	  _domain[i] = v;
	  if (v < min) min = v;
	  if (v > max) max = v;
	}
    _iinf = min;
    _isup = max;
	_wcspIndex = -1;
}

int Toulbar2_Expression::get_value() const
{
  Value var = -1;  
  if(_solver->is_sat()) {
    var = _solver->solution[_wcspIndex];
  } else if (_solver->wcsp->assigned(_wcspIndex)) {
    var = _solver->wcsp->getValue(_wcspIndex);    
  }
  return var;
}

unsigned int Toulbar2_Expression::get_size() const
{
  return _solver->wcsp->getDomainInitSize(_wcspIndex);
}

int Toulbar2_Expression::get_max() const
{
  return _isup;
}

int Toulbar2_Expression::get_min() const
{
  return _iinf;
}

int Toulbar2_Expression::next(int v)
{
  int nxt = v;

  if(_wcspIndex >= 0 && _wcspIndex < _solver->wcsp->numberOfVariables())  {
	nxt = _solver->wcsp->nextValue(_wcspIndex, v);
  }

  return nxt;
}

bool Toulbar2_Expression::contain(const int v) const
{
  return _solver->wcsp->canbe(_wcspIndex, v);
}

bool Toulbar2_Expression::has_been_added() const
{
    return (_solver != NULL);
}

Toulbar2_Expression* Toulbar2_Expression::add(Toulbar2Solver *solver, bool top_level)
{
  if(!has_been_added()) {
	_solver = solver;
	string index = to_string(_solver->wcsp->numberOfVariables());
	if (_domain) {
	  _wcspIndex = _solver->wcsp->makeEnumeratedVariable(index,_domain,_size);
	} else {
	  _wcspIndex = _solver->wcsp->makeEnumeratedVariable(index,_iinf,_isup);
	}
#ifdef _DEBUGWRAP    
	cout << "CREATE " << _wcspIndex << endl;
#endif
  } else {
#ifdef _DEBUGWRAP    
	cout << "REUSE "  << _wcspIndex << endl;
#endif
  }
  if(top_level) {
#ifdef _DEBUGWRAP    
	cout << "Adding at top level an expression variable which must be True, i.e., different from zero: " << _wcspIndex << endl; 
#endif
	_solver->wcsp->remove(_wcspIndex, 0);
  }
  
  return this;
}


/**
 * Constraints 
 */

Toulbar2_AllDiff::Toulbar2_AllDiff( Toulbar2ExpArray& vars ) 
  : Toulbar2_Expression() 
{
  _vars = vars;
  _scope = new int[_vars.size()];
}

Toulbar2_AllDiff::Toulbar2_AllDiff( Toulbar2_Expression *var1, Toulbar2_Expression *var2 ) 
  : Toulbar2_Expression() 
{
  _vars.add(var1);
  _vars.add(var2);
  _scope = new int[2];
}

Toulbar2_Expression* Toulbar2_AllDiff::add(Toulbar2Solver *solver, bool top_level)
{
  if (!has_been_added()) {
    _solver = solver;
    if(top_level) {
	  int maxdom = 0;
	  for(int i = 0; i < _vars.size(); i++) {  
        _vars.get_item(i)->add(solver, false);  
#ifdef _DEBUGWRAP    
		cout << "add AllDiff var: " << _vars.get_item(i)->_wcspIndex << " " << ((_vars.get_item(i)->_wcspIndex>=0)?_solver->wcsp->getName(_vars.get_item(i)->_wcspIndex):".") << " " << ((_vars.get_item(i)->_wcspIndex>=0)?_solver->wcsp->getDomainInitSize(_vars.get_item(i)->_wcspIndex):0) << endl;
#endif
        _scope[i] = _vars.get_item(i)->_wcspIndex;
		if (_vars.get_item(i)->get_size() > maxdom) maxdom = _vars.get_item(i)->get_size();
	  }
       //  string alldiff = "salldiff";
       //  string parameters_("hard " + cost2string(_solver->wcsp->getUb()));
       //  istringstream parameters(parameters_);
       // _solver->wcsp->postGlobalConstraint(_scope, _vars.size(), alldiff, parameters);
      if (_vars.size() == 2 || maxdom > ((_vars.size()-1)/4)) {
		// decompose AllDiff into a clique of binary disequality constraints
        for(int i = 0; i < _vars.size(); i++) {
		  for(int j = i+1; j < _vars.size(); j++) {
			_solver->wcsp->postDisjunction(_vars.get_item(i)->_wcspIndex, _vars.get_item(j)->_wcspIndex, 1, 1, _solver->wcsp->getUb());
		  }
		}
	  } else {
		// decompose AllDiff into a chain of ternary constraints for each domain value
        _solver->wcsp->postWAllDiff(_scope, _vars.size(), "hard", _solver->wcsp->getUb());
      }
    }
  }
  return this;
}


Toulbar2_Gcc::Toulbar2_Gcc(Toulbar2ExpArray& vars, Toulbar2IntArray& vals, Toulbar2IntArray& lb_card, Toulbar2IntArray& ub_card) : Toulbar2_Expression() 
{    
  _vars = vars;
  _vals = vals;
  _lb_card = lb_card;
  _ub_card = ub_card;
  _scope = new int[_vars.size()];
  _nbValues = vals.size();
  _values = new int[_nbValues];
  for(int i = 0; i < _nbValues; i++) {
	_values[i] = vals.get_item(i);  
  }
  _lb = new int[_nbValues];
  for(int i = 0; i < _nbValues; i++) {
	_lb[i] = lb_card.get_item(i);  
  }
  _ub = new int[_nbValues];
  for(int i = 0; i < _nbValues; i++) {
	_ub[i] = ub_card.get_item(i);  
  }  
}

Toulbar2_Expression* Toulbar2_Gcc::add(Toulbar2Solver *solver, bool top_level)
{
  if (!has_been_added()) {
    _solver = solver;  
    if(top_level) {
	  for(int i = 0; i < _vars.size(); i++) {  
		_vars.get_item(i)->add(solver, false);  
#ifdef _DEBUGWRAP    
		cout << "add Gcc var: " << _vars.get_item(i)->_wcspIndex << " " << ((_vars.get_item(i)->_wcspIndex>=0)?_solver->wcsp->getName(_vars.get_item(i)->_wcspIndex):".") << " " << ((_vars.get_item(i)->_wcspIndex>=0)?_solver->wcsp->getDomainInitSize(_vars.get_item(i)->_wcspIndex):0) << endl;
#endif
		_scope[i] = _vars.get_item(i)->_wcspIndex;
	  }
	  _solver->wcsp->postWGcc(_scope, _vars.size(), "hard", _solver->wcsp->getUb(), _values, _nbValues, _lb, _ub);
	  // string gcc = "sgcc";
	  // string parameters_("hard " + cost2string(_solver->wcsp->getUb()) + " ");
	  // parameters_ = parameters_ + to_string(_vals.size()) + " ";
	  // for(int i = 0; i < _vals.size(); i++) {
	  // 	parameters_ = parameters_ + to_string(_vals.get_item(i)) + " " + to_string(_lb_card.get_item(i)) + " " + to_string(_ub_card.get_item(i)) + " ";    
	  // }            
	  // istringstream parameters(parameters_);
	  // _solver->wcsp->postGlobalConstraint(_scope, arity, gcc, parameters);
	}
  }
  return this;
}

Toulbar2_PostNullary::Toulbar2_PostNullary(int cost)
  : Toulbar2_Expression() 
{
  _cost = cost;
}

Toulbar2_Expression* Toulbar2_PostNullary::add(Toulbar2Solver *solver, bool top_level)
{
  if(!has_been_added()) 
  {
    _solver = solver;
    if(top_level) {
#ifdef _DEBUGWRAP    
	  cout << "increasLb by " << _cost << endl;
#endif
	  if (_cost > 0) {
		_solver->wcsp->increaseLb(_cost);
	  } else {
		solver->costshift += _cost;
	  }
    } 
  }
  return this;
}

Toulbar2_PostUnary::Toulbar2_PostUnary(Toulbar2_Expression* var, Toulbar2IntArray& costs)
  : Toulbar2_Expression() 
{
  _var = var;
  for(int i=0;i<costs.size();i++)
  {
    _costs.push_back(costs.get_item(i));
  }
}

Toulbar2_Expression* Toulbar2_PostUnary::add(Toulbar2Solver *solver, bool top_level)
{
  if(!has_been_added()) 
  {
    _solver = solver;
    _var->add(_solver,false);
    if(top_level) {
#ifdef _DEBUGWRAP    
	  cout << "postUnary on varIndex " << _var->_wcspIndex << " with costs: " << _costs[0] << " " << _costs[1] << ".." << endl;
#endif
	  Cost mincost = 0;
	  for (unsigned int i=0; i<_costs.size(); i++) {
		if (_costs[i] < mincost) {
		  mincost = _costs[i];
		}
	  }
	  if (mincost < 0) {
		for (unsigned int i=0; i<_costs.size(); i++) {
		  _costs[i] -= mincost;
		}
		_solver->costshift += mincost;
	  }
      _solver->wcsp->postUnary(_var->_wcspIndex,_costs);
    } 
  }
  return this;
}

Toulbar2_PostBinary::Toulbar2_PostBinary(Toulbar2_Expression* var1, Toulbar2_Expression* var2, Toulbar2IntArray& costs)
  : Toulbar2_Expression() 
{
  _var1 = var1;
  _var2 = var2;
  for(int i=0;i<costs.size();i++)
  {
    _costs.push_back(costs.get_item(i));
  }
}

Toulbar2_Expression* Toulbar2_PostBinary::add(Toulbar2Solver *solver, bool top_level)
{
  if(!has_been_added()) {
    _solver = solver;
    _var1->add(_solver,false);
    _var2->add(_solver,false);
    if(top_level) {
	  Cost mincost = 0;
	  for (unsigned int i=0; i<_costs.size(); i++) {
		if (_costs[i] < mincost) {
		  mincost = _costs[i];
		}
	  }
	  if (mincost < 0) {
		for (unsigned int i=0; i<_costs.size(); i++) {
		  _costs[i] -= mincost;
		}
		_solver->costshift += mincost;
	  }
	  if (_solver->wcsp->getDomainInitSize(_var1->_wcspIndex) * _solver->wcsp->getDomainInitSize(_var2->_wcspIndex) != _costs.size()) {
		vector<Cost> newcosts(_solver->wcsp->getDomainInitSize(_var1->_wcspIndex) * _solver->wcsp->getDomainInitSize(_var2->_wcspIndex), 0);
		int id1 = 0;
		for (Value v1 = _solver->wcsp->getInf(_var1->_wcspIndex); v1 <= _solver->wcsp->getSup(_var1->_wcspIndex) && id1 < _solver->wcsp->getDomainSize(_var1->_wcspIndex); id1++, v1 = _solver->wcsp->nextValue(_var1->_wcspIndex, v1)) {
		  int id2 = 0;
		  for (Value v2 = _solver->wcsp->getInf(_var2->_wcspIndex); v2 <= _solver->wcsp->getSup(_var2->_wcspIndex) && id2 < _solver->wcsp->getDomainSize(_var2->_wcspIndex); id2++, v2 = _solver->wcsp->nextValue(_var2->_wcspIndex, v2)) {
			newcosts[_solver->wcsp->toIndex(_var1->_wcspIndex, v1) * _solver->wcsp->getDomainInitSize(_var2->_wcspIndex) + _solver->wcsp->toIndex(_var2->_wcspIndex, v2)] = _costs[id1 * _solver->wcsp->getDomainSize(_var2->_wcspIndex) + id2];
		  }
		}
		_solver->wcsp->postBinaryConstraint(_var1->_wcspIndex, _var2->_wcspIndex, newcosts);
	  } else {
		_solver->wcsp->postBinaryConstraint(_var1->_wcspIndex, _var2->_wcspIndex, _costs);
	  }
	}
  }
  return this;
}

Toulbar2_PostTernary::Toulbar2_PostTernary(Toulbar2ExpArray& vars, Toulbar2IntArray& costs)
  : Toulbar2_Expression() 
{
  _var1 = vars.get_item(0);
  _var2 = vars.get_item(1);
  _var3 = vars.get_item(2);
  for(int i=0;i<costs.size();i++)
  {
    _costs.push_back(costs.get_item(i));
  }
}

Toulbar2_Expression* Toulbar2_PostTernary::add(Toulbar2Solver *solver, bool top_level)
{
  if(!has_been_added()) {
    _solver = solver;
    _var1->add(_solver,false);
    _var2->add(_solver,false);
    _var3->add(_solver,false);
    if(top_level) {
	  Cost mincost = 0;
	  for (unsigned int i=0; i<_costs.size(); i++) {
		if (_costs[i] < mincost) {
		  mincost = _costs[i];
		}
	  }
	  if (mincost < 0) {
		for (unsigned int i=0; i<_costs.size(); i++) {
		  _costs[i] -= mincost;
		}
		_solver->costshift += mincost;
	  }
	  if (_solver->wcsp->getDomainInitSize(_var1->_wcspIndex) * _solver->wcsp->getDomainInitSize(_var2->_wcspIndex) * _solver->wcsp->getDomainInitSize(_var3->_wcspIndex) != _costs.size()) {
		vector<Cost> newcosts(_solver->wcsp->getDomainInitSize(_var1->_wcspIndex) * _solver->wcsp->getDomainInitSize(_var2->_wcspIndex) * _solver->wcsp->getDomainInitSize(_var3->_wcspIndex), 0);
		int id1 = 0;
		for (Value v1 = _solver->wcsp->getInf(_var1->_wcspIndex); v1 <= _solver->wcsp->getSup(_var1->_wcspIndex) && id1 < _solver->wcsp->getDomainSize(_var1->_wcspIndex); id1++, v1 = _solver->wcsp->nextValue(_var1->_wcspIndex, v1)) {
		  int id2 = 0;
		  for (Value v2 = _solver->wcsp->getInf(_var2->_wcspIndex); v2 <= _solver->wcsp->getSup(_var2->_wcspIndex) && id2 < _solver->wcsp->getDomainSize(_var2->_wcspIndex); id2++, v2 = _solver->wcsp->nextValue(_var2->_wcspIndex, v2)) {
			int id3 = 0;
			for (Value v3 = _solver->wcsp->getInf(_var3->_wcspIndex); v3 <= _solver->wcsp->getSup(_var3->_wcspIndex) && id3 < _solver->wcsp->getDomainSize(_var3->_wcspIndex); id3++, v3 = _solver->wcsp->nextValue(_var3->_wcspIndex, v3)) {
			  newcosts[_solver->wcsp->toIndex(_var1->_wcspIndex, v1) * _solver->wcsp->getDomainInitSize(_var2->_wcspIndex) * _solver->wcsp->getDomainInitSize(_var3->_wcspIndex) + _solver->wcsp->toIndex(_var2->_wcspIndex, v2) * _solver->wcsp->getDomainInitSize(_var3->_wcspIndex) + _solver->wcsp->toIndex(_var3->_wcspIndex, v3)] = _costs[id1 * _solver->wcsp->getDomainSize(_var2->_wcspIndex) * _solver->wcsp->getDomainSize(_var3->_wcspIndex) + id2 * _solver->wcsp->getDomainSize(_var3->_wcspIndex) + id3];
			}
		  }
		}
		_solver->wcsp->postTernaryConstraint(_var1->_wcspIndex, _var2->_wcspIndex, _var3->_wcspIndex, newcosts);
	  } else {
		_solver->wcsp->postTernaryConstraint(_var1->_wcspIndex, _var2->_wcspIndex, _var3->_wcspIndex, _costs);
	  }
	} 
  }
  return this;
}

Toulbar2_PostNary::Toulbar2_PostNary(Toulbar2ExpArray& vars, int arity, int default_cost, Toulbar2IntMultiArray& values, Toulbar2IntArray& costs)
  : Toulbar2_Expression() 
{
  _vars = vars;
  _arity = arity; 
  _defcost = default_cost; 
  _values = values;
  _costs = costs;
  _scope = new int[arity];
}

Toulbar2_Expression* Toulbar2_PostNary::add(Toulbar2Solver *solver, bool top_level)
{
  if(!has_been_added()) {
    _solver = solver;
    for(int i = 0; i < _vars.size(); i++) {
      _vars.get_item(i)->add(solver, false);
      _scope[i] = _vars.get_item(i)->_wcspIndex;  
    }
    if(top_level) {
	  Cost mincost = _defcost;
	  for (unsigned int i=0; i<_costs.size(); i++) {
		if (_costs.get_item(i) < mincost) {
		  mincost = _costs.get_item(i);
		}
	  }
	  if (mincost < 0) {
		_defcost -= mincost;
		for (unsigned int i=0; i<_costs.size(); i++) {
		  _costs.set_item(i, _costs.get_item(i) - mincost);
		}
		_solver->costshift += mincost;
	  }	  
	  int ctrIndex = _solver->wcsp->postNaryConstraintBegin(_scope,_arity,_defcost);
	  if(_costs.size() != 0) { 
		for(int i = 0; i < _costs.size(); i++) {
		  Value val[_arity];  
		  for(int j = 0; j < _arity; j++) {
			val[j] = _values.get_item(i*_arity + j);
		  }
		  _solver->wcsp->postNaryConstraintTuple(ctrIndex,val,_arity,_costs.get_item(i));
		}
	  }
	  _solver->wcsp->postNaryConstraintEnd(ctrIndex);
    } 
  }
  return this;
}

Toulbar2_PostWSum::Toulbar2_PostWSum(Toulbar2ExpArray& vars, int arity, const char* semantics, const char* baseCost, const char* comparator, int rightRes): Toulbar2_Expression() 
{
  _vars = vars;
  _arity = arity; 
  string seman(semantics);
  _semantics = seman; 
  _baseCost = string2Cost(baseCost);
  string comp(comparator); 
  _comparator = comp;
  _rightRes = rightRes;
  _scope = new int[_vars.size()];
}

Toulbar2_Expression* Toulbar2_PostWSum::add(Toulbar2Solver *solver, bool top_level)
{
  if(!has_been_added()) {
    _solver = solver;
    for(int i = 0; i < _vars.size(); i++) {
      _vars.get_item(i)->add(solver, false);
      _scope[i] = _vars.get_item(i)->_wcspIndex;
    }
    if(top_level){ 
      _solver->wcsp->postWSum(_scope, _arity, _semantics, _baseCost, _comparator, _rightRes);
    } 
  }
  return this;
}

Toulbar2_PostWVarSum::Toulbar2_PostWVarSum(Toulbar2ExpArray& vars, int arity, const char* semantics, const char* baseCost, const char* comparator): Toulbar2_Expression() 
{
  _vars = vars;
  _arity = arity; 
  string seman(semantics);
  _semantics = seman; 
  _baseCost = string2Cost(baseCost);
  string comp(comparator); 
  _comparator = comp;
  _scope = new int[_vars.size()];
}

Toulbar2_Expression* Toulbar2_PostWVarSum::add(Toulbar2Solver *solver, bool top_level)
{
  if(!has_been_added()) {
    _solver = solver;
    for(int i = 0; i < _vars.size(); i++) {
      _vars.get_item(i)->add(solver, false);
      _scope[i] = _vars.get_item(i)->_wcspIndex;
    }
    if(top_level){ 
      _solver->wcsp->postWVarSum(_scope, _arity, _semantics, _baseCost, _comparator, _scope[_vars.size()-1]);
    } 
  }
  return this;
}

Toulbar2_PostWAmong::Toulbar2_PostWAmong(Toulbar2ExpArray& vars, int arity, const char* semantics, const char* baseCost, Toulbar2IntArray& values, int lb, int ub): Toulbar2_Expression() 
{
  _vars = vars;
  _arity = arity; 
  string seman(semantics);
  _semantics = seman; 
  _baseCost = string2Cost(baseCost);
  _values = new int[values.size()];
  _nbValues = values.size();
  for(int i = 0; i < _nbValues; i++) {
    _values[i] = values.get_item(i);
  }
  _lb = lb;
  _ub = ub;
  _scope = new int[_vars.size()];
}

Toulbar2_PostWAmong::Toulbar2_PostWAmong(Toulbar2ExpArray& vars, int arity, const char* semantics, const char* baseCost, Toulbar2IntArray& values): Toulbar2_Expression() 
{
  _vars = vars;
  _arity = arity; 
  string seman(semantics);
  _semantics = seman; 
  _baseCost = string2Cost(baseCost);
  _values = new int[values.size()];
  _nbValues = values.size();
  for(int i = 0; i < _nbValues; i++)
  {
    _values[i] = values.get_item(i);
  }
  _lb = -1;
  _ub = -1;
  _scope = new int[_vars.size()];
}

Toulbar2_Expression* Toulbar2_PostWAmong::add(Toulbar2Solver *solver, bool top_level)
{
  if(!has_been_added()) {
    _solver = solver;
    for(int i = 0; i < _vars.size(); i++) {
      _vars.get_item(i)->add(solver, false);
      _scope[i] = _vars.get_item(i)->_wcspIndex;
    }
    if (top_level) { 
      if(_lb != -1 && _ub != -1) {  
        _solver->wcsp->postWAmong(_scope, _arity, _semantics, _baseCost, _values, _nbValues, _lb, _ub);
      } else {
        _solver->wcsp->postWVarAmong(_scope, _arity, _semantics, _baseCost, _values, _nbValues, _scope[_vars.size()-1]);
      }
    } 
  }
  return this;
}

Toulbar2_Regular::Toulbar2_Regular(Toulbar2ExpArray& vars, int arity, int nbStates, Toulbar2IntArray& initialStates, Toulbar2IntArray& acceptingStates, Toulbar2IntMultiArray& transitions, 
Toulbar2IntArray& initialCosts, Toulbar2IntArray& acceptingCosts, Toulbar2IntArray& transitionsCosts): Toulbar2_Expression() 
{
  _vars = vars;
  _arity = arity; 
  _nbStates = nbStates;
  _transitions = new int* [transitions.size()/3];
  for(int i = 0; i < transitions.size()/3; i++)
  {
    _transitions[i] = new int[3];
  }
  for(int i = 0; i < initialStates.size(); i++)
  {
    _initialStates.push_back(make_pair(initialStates.get_item(i), initialCosts.get_item(i)));
  }
  for(int i = 0; i < acceptingStates.size(); i++)
  {
    _acceptingStates.push_back(make_pair(acceptingStates.get_item(i), acceptingCosts.get_item(i)));
  }
  for(int i = 0; i < transitions.size()/3; i++)
  {
    for(int j =0; j < 3; j++)
    {  
      _transitions[i][j] = transitions.get_item(i*3+j);
    }
    _transitionsCosts.push_back(transitionsCosts.get_item(i));
  }
  _scope = new int[_vars.size()];
  _type = "w";
}

Toulbar2_Regular::Toulbar2_Regular(Toulbar2ExpArray& vars, int arity, int nbStates, Toulbar2IntArray& initialStates, Toulbar2IntArray& acceptingStates, Toulbar2IntMultiArray& transitions, 
const char* type, const char* measureCost,  const char* baseCost): Toulbar2_Expression() 
{
  _vars = vars;
  _arity = arity; 
  _nbStates = nbStates;
  _initialStatesV = initialStates;
  _acceptingStatesV = acceptingStates;
  _transitionsV = transitions;
  _transitions = NULL;
  _scope = new int[_vars.size()];      
  string temp(type);
  _type = temp; 
  string tempCost(measureCost);
  _measureCost = tempCost;   
  string tempVal(baseCost);
  _baseCost = tempVal;     
}

Toulbar2_Expression* Toulbar2_Regular::add(Toulbar2Solver *solver, bool top_level)
{
  if(!has_been_added()) {
    _solver = solver;
    for(int i = 0; i < _vars.size(); i++) {
      _vars.get_item(i)->add(solver, false);
      _scope[i] = _vars.get_item(i)->_wcspIndex;
    }
    if (top_level) { 
	  _solver->wcsp->postWRegular(_scope, _arity, _nbStates, _initialStates, _acceptingStates, _transitions, _transitionsCosts);
	  // string regular = "sregular";
	  // string parameters_(_measureCost + " " + _baseCost + " " + to_string(_nbStates) + " " + to_string(_initialStatesV.size()) + " ");
	  // for(int i = 0; i < _initialStatesV.size(); i++)
      //   {
      //     parameters_ = parameters_ + to_string(_initialStatesV.get_item(i)) + " ";    
      //   }
	  // parameters_ = parameters_ + to_string(_acceptingStatesV.size()) + " ";
	  // for(int i = 0; i < _acceptingStatesV.size(); i++)
      //   {
      //     parameters_ = parameters_ + to_string(_acceptingStatesV.get_item(i)) + " ";    
      //   }
	  // parameters_ = parameters_ + to_string(_transitionsV.size()) + " ";
	  // for(int i = 0; i < _transitionsV.size(); i++)
      //   {
      //     parameters_ = parameters_ + to_string(_transitionsV.get_item(i).get_item(0)) + " " + to_string(_transitionsV.get_item(i).get_item(1)) 
      //     + " " + to_string(_transitionsV.get_item(i).get_item(2)) + " ";    
      //   }
	  // istringstream parameters(parameters_);
	  // _solver->wcsp->postGlobalConstraint(_scope, _arity, regular, parameters);           
	}
  } 
  return this;
}

Toulbar2_Same::Toulbar2_Same(Toulbar2ExpArray& vars)
{
  _vars = vars;
  _arity = _vars.size();
  _scope = new int[_arity];
  _type = "n";
}

Toulbar2_Same::Toulbar2_Same(Toulbar2ExpArray& vars, const char* type, const char* costValue)
{
  _vars = vars;
  _arity = _vars.size();
  _scope = new int[_arity];
  string tempType(type);
  _type = tempType;
  string tempCost(costValue);
  _costValue = tempCost;    
}

Toulbar2_Same::Toulbar2_Same(Toulbar2ExpArray& vars, const char* type, const char* semantics, const char* baseCost)
{
  _vars = vars;
  _arity = _vars.size();
  _scope = new int[_arity];
  string seman(semantics);
  _semantics = seman; 
    string tempType(type);
  _type = tempType; 
  _baseCost = string2Cost(baseCost);  
}

Toulbar2_Expression* Toulbar2_Same::add(Toulbar2Solver *solver, bool top_level)
{
  if(!has_been_added()) {
    _solver = solver;
    if(_vars.size()%2 == 0) {
      for(int i = 0; i < _vars.size()/2; i++) {
        _vars.get_item(i)->add(solver, false);  
        _scope[i] = _vars.get_item(i)->_wcspIndex;
      }
      for(int i = _vars.size()/2; i < _vars.size(); i++) {
        _vars.get_item(i)->add(solver, false);  
        _scope[i] = _vars.get_item(i)->_wcspIndex;
      }
      if (top_level) {
		_solver->wcsp->postWSame(_scope, _arity, _semantics, _baseCost);    
		// string same = "ssame";
		// string parameters_ = "";
		// if(_type == "n")
		//   parameters_ = to_string(MAXCOST) + " " + to_string(_vars.size()/2) + " " + to_string(_vars.size()/2) + " ";
		// else
		//   parameters_ = _costValue + " " + to_string(_vars.size()/2) + " " + to_string(_vars.size()/2) + " ";
		// for(int i = 0; i < _vars.size();i++) 
        //   {
  		// 	parameters_ = parameters_ + to_string(_vars.get_item(i)->_wcspIndex) + " "; 
        //   }  
		// istringstream parameters(parameters_);
		// _solver->wcsp->postGlobalConstraint(_scope, _arity, same, parameters);    
	  }
    } else {
      cout << "Error on Same constraint: the two lists must have the same size (Same constraint not added)!" << endl; 
    }
  }
  return this;
}

Toulbar2_PostWSameGcc::Toulbar2_PostWSameGcc(Toulbar2ExpArray& vars, Toulbar2IntArray& vals, Toulbar2IntArray& lb_card, Toulbar2IntArray& ub_card, const char* type, const char* semantics, const char* baseCost)
{
  _vars = vars;
  _scope = new int[_vars.size()];
  _nbValues = vals.size();
  _values = new int[_nbValues];
  for(int i = 0; i < _nbValues; i++)
  {
    _values[i] = vals.get_item(i);  
  }
  _lb = new int[_nbValues];
  for(int i = 0; i < _nbValues; i++)
  {
    _lb[i] = lb_card.get_item(i);  
  }
  _ub = new int[_nbValues];
  for(int i = 0; i < _nbValues; i++)
  {
    _ub[i] = ub_card.get_item(i);  
  }
  string temp(type);
  _type = temp; 
  string tempSeman(semantics);
  _semantics = tempSeman;   
  _baseCost = string2Cost(baseCost);  
}

Toulbar2_Expression* Toulbar2_PostWSameGcc::add(Toulbar2Solver *solver, bool top_level)
{
  if(!has_been_added()) {
    _solver = solver;  
    for(int i = 0; i < _vars.size(); i++) {  
      _vars.get_item(i)->add(solver, false);  
      _scope[i] = _vars.get_item(i)->_wcspIndex;
    }
    if(top_level) {
	  _solver->wcsp->postWSameGcc(_scope, _vars.size(), _semantics, _baseCost, _values, _nbValues, _lb, _ub);
	}
  }
  return this;
} 

Toulbar2_PostWOverlap::Toulbar2_PostWOverlap(Toulbar2ExpArray& vars, int arity, const char* semantics, const char* baseCost, const char* comparator, int rightRes): Toulbar2_Expression() 
{
  _vars = vars;
  _arity = arity; 
  string seman(semantics);
  _semantics = seman; 
  _baseCost = string2Cost(baseCost);
  string comp(comparator); 
  _comparator = comp;
  _rightRes = rightRes;
  _scope = new int[_vars.size()];
}

Toulbar2_Expression* Toulbar2_PostWOverlap::add(Toulbar2Solver *solver, bool top_level)
{
  if(!has_been_added()) {
    _solver = solver;
    for(int i = 0; i < _vars.size(); i++) {
      _vars.get_item(i)->add(solver, false);
      _scope[i] = _vars.get_item(i)->_wcspIndex;
    }
    if (top_level) { 
      _solver->wcsp->postWOverlap(_scope, _arity, _semantics, _baseCost, _comparator, _rightRes);
    } 
  }
  return this;
}

/******************** Table constraint *************************/

Toulbar2_Table::Toulbar2_Table( Toulbar2ExpArray& vars, Toulbar2IntArray& tuples, const char* type ) 
  : Toulbar2_Expression() 
{
#ifdef _DEBUGWRAP
  cout << "Table scope:";
  for(int i=0; i<vars.size(); i++) cout << " " << vars.get_item(i)->_wcspIndex;
  cout << endl;
#endif
#ifndef NDEBUG
  for(int i=0; i<vars.size(); i++) for (int j=i+1; j<vars.size(); j++) assert(vars.get_item(i) != vars.get_item(j));
#endif
  _vars = vars;
  _tuples = tuples;
  _arity = _vars.size();
  _scope = new int[_arity]; 
  if(strcmp(type,"support") == 0) {
    _spin = 0;
  } else if (strcmp(type,"conflict") == 0) {
    _spin = 1;
  } else {
    _spin = 0;
    cout << "Warning: only \"support\" or \"conflict\" are allowed as Table constraint type (set to \"support\" by default)!" << endl;     
  }
}

Toulbar2_Table::Toulbar2_Table( Toulbar2_Expression *var1, Toulbar2_Expression *var2, Toulbar2IntArray& tuples, const char* type ) : Toulbar2_Expression() 
{
  _vars.add(var1);
  _vars.add(var2); 
  _tuples = tuples;
  _scope = new int[2];
  _arity = 2;
  if(strcmp(type,"support") == 0) {
    _spin = 0;
  } else if (strcmp(type,"conflict") == 0) {
    _spin = 1;
  } else {
    _spin = 0;
    cout << "Warning: only \"support\" or \"conflict\" are allowed as Table constraint type (set to \"support\" by default)!" << endl;     
  }
}

Toulbar2_Expression* Toulbar2_Table::add(Toulbar2Solver *solver, bool top_level) 
{
  if (!has_been_added()) {
    _solver = solver;
    for(int i = 0; i < _arity; i++) {
      _vars.get_item(i)->add(solver, false);
#ifdef _DEBUGWRAP    
		cout << "add Table var: " << _vars.get_item(i)->_wcspIndex << " " << ((_vars.get_item(i)->_wcspIndex>=0)?_solver->wcsp->getName(_vars.get_item(i)->_wcspIndex):".") << " " << ((_vars.get_item(i)->_wcspIndex>=0)?_solver->wcsp->getDomainInitSize(_vars.get_item(i)->_wcspIndex):0) << endl;
#endif
      _scope[i] = _vars.get_item(i)->_wcspIndex; 
    }
    if(top_level) { 
	  Cost top = solver->wcsp->getUb();
      if(_arity == 2) {
        int domainSizeV0 = _vars.get_item(0)->get_size();
        int domainSizeV1 = _vars.get_item(1)->get_size();
        vector<Cost> _costs(domainSizeV0*domainSizeV1, (_spin == 0)?top:0);
        for(int l = 0; l < _tuples.size(); l += 2) {
		  int idx = _solver->wcsp->toIndex(_vars.get_item(0)->_wcspIndex, _tuples.get_item(l)) * domainSizeV1 + _solver->wcsp->toIndex(_vars.get_item(1)->_wcspIndex, _tuples.get_item(l+1));
		  assert(idx >= 0);
		  assert(idx < _costs.size());
		  _costs[ idx ] = ((_spin == 0)?0:top);
		}		
        _solver->wcsp->postBinaryConstraint(_vars.get_item(0)->_wcspIndex, _vars.get_item(1)->_wcspIndex, _costs);
      } else if (_arity == 3) {
        int domainSizeV0 = _vars.get_item(0)->get_size();
        int domainSizeV1 = _vars.get_item(1)->get_size();
        int domainSizeV2 = _vars.get_item(2)->get_size();
        vector<Cost> _costs(domainSizeV0*domainSizeV1*domainSizeV2 , (_spin == 0)?top:0);
        for(int l = 0; l < _tuples.size(); l += 3) {
		  int idx =  _solver->wcsp->toIndex(_vars.get_item(0)->_wcspIndex, _tuples.get_item(l)) * domainSizeV1 * domainSizeV2 + _solver->wcsp->toIndex(_vars.get_item(1)->_wcspIndex, _tuples.get_item(l+1)) * domainSizeV2 + _solver->wcsp->toIndex(_vars.get_item(2)->_wcspIndex, _tuples.get_item(l+2));
		  assert(idx >= 0);
		  assert(idx < _costs.size());
		  _costs[ idx ] = ((_spin == 0)?0:top);
		}
        _solver->wcsp->postTernaryConstraint(_vars.get_item(0)->_wcspIndex, _vars.get_item(1)->_wcspIndex, _vars.get_item(2)->_wcspIndex, _costs);          
      } else if(_arity > 3) {
        int ctrIndex = _solver->wcsp->postNaryConstraintBegin(_scope,_arity, (_spin == 0)?top:0);
        if(_tuples.size() != 0) { 
          int i = 0;  
          while(i < _tuples.size()) {
            Value val[_arity];  
            for(int j = 0; j < _arity; j++) {
              val[j] = _tuples.get_item(i);
			  assert(_solver->wcsp->toIndex(_scope[j], val[j]) >= 0);
			  assert(_solver->wcsp->toIndex(_scope[j], val[j]) < _solver->wcsp->getDomainInitSize(_scope[j]));
              i++;
            }
            _solver->wcsp->postNaryConstraintTuple(ctrIndex,val,_arity, (_spin == 0)?0:top);
          }
        }
        _solver->wcsp->postNaryConstraintEnd(ctrIndex);
      } else {
        cout << "Table arity must be > 1, Table constraint not added..." << endl; 
      }
    } 
  }
  return this;
}


/* Minimise object */

// Toulbar2_Minimise::Toulbar2_Minimise(Toulbar2_Expression *var)
//   : Toulbar2_Expression()
// {
//   _exp = var;
// }

// Toulbar2_Expression* Toulbar2_Minimise::add(Toulbar2Solver *solver, bool top_level)
// {
//   if(!has_been_added()) {
//       _solver = solver;
//       // This will be the objective function
//       _exp->add(solver, false);
// #ifdef _DEBUGWRAP    
// 		cout << "add Minimise var: " << _exp->_wcspIndex << " " << ((_exp->_wcspIndex>=0)?_solver->wcsp->getName(_exp->_wcspIndex):".") << " " << ((_exp->_wcspIndex>=0)?_solver->wcsp->getDomainInitSize(_exp->_wcspIndex):0) << endl;
// #endif
//       if(top_level) {
// 		_solver->costshift += _solver->wcsp->getInf(_exp->_wcspIndex);
// 		vector<Cost> costs;
// 		for(unsigned int i=0; i < _solver->wcsp->getDomainInitSize(_exp->_wcspIndex); i++) {
// 		  costs.push_back(_solver->wcsp->toValue(_exp->_wcspIndex,i) - _solver->wcsp->getInf(_exp->_wcspIndex));    
// 		}
// 		_solver->wcsp->postUnary(_exp->_wcspIndex,costs);
// 	  }
//   }
//   return this;
// }

/**************************************************************
 ********************     Solver        ***********************
 **************************************************************/

Toulbar2Solver::Toulbar2Solver()
{
  tb2init();
  tb2checkOptions(MAXCOST);
  ToulBar2::verbose = -1;
  ToulBar2::startCpuTime = cpuTime();
  initCosts(MAXCOST);
  solver = WeightedCSPSolver::makeWeightedCSPSolver(MAXCOST);
  wcsp = solver->getWCSP();
  upperbound = MAXCOST;
  optimum = MAXCOST;
  costshift = MIN_COST;
  unsatisfiable = false;
  interrupted = false;
  MyWCSPSolver = this;
}

void Toulbar2Solver::add(Toulbar2_Expression* arg)
{
  try {
	arg->add(this, true);
  } catch (Contradiction) {
	wcsp->whenContradiction();
	unsatisfiable = true;
  }
#ifdef _DEBUGWRAP    
  cout << "add top-level var: " << arg->_wcspIndex << " " << ((arg->_wcspIndex>=0)?wcsp->getName(arg->_wcspIndex):".") << " " << ((arg->_wcspIndex>=0)?wcsp->getDomainInitSize(arg->_wcspIndex):0) << endl;
#endif
}

void Toulbar2Solver::initialise()
{
  ToulBar2::setvalue = NULL;
}

bool Toulbar2Solver::propagate()
{
  if (!unsatisfiable) {
	try {
	  wcsp->propagate();
	} catch (Contradiction) {
	  wcsp->whenContradiction();
	  unsatisfiable = true;
	}
  }
  return (!unsatisfiable);
}

int Toulbar2Solver::solve()
{
  if (!unsatisfiable) {
	wcsp->sortConstraints();
	if (ToulBar2::dumpWCSP == 1) {
	  solver->dump_wcsp(ToulBar2::problemsaved_filename.c_str());
	} else {
	  ToulBar2::timeOut = timeout;
	  signal(SIGINT,  timeOut);
	  try {
		try {
		  upperbound = wcsp->getUb();
		  unsatisfiable = !solver->solve();
		} catch (Contradiction) {
		  wcsp->whenContradiction();
		  unsatisfiable = true;
		}
	  } catch (TimeOut) {
		interrupted = true;
	  }
	  timerStop();
	  if(wcsp->getUb() < upperbound) {
		optimum = solver->getSolution(solution);
	  }
	}
  }
  return is_sat();
}

int Toulbar2Solver::solveAndRestart(const int policy, 
           const unsigned int base, 
           const double factor,
           const double decay,
           const int reinit)
{
  restart(base);
  return Toulbar2Solver::solve();
}

void Toulbar2Solver::printPython()
{
  if (optimum < upperbound) {
	if (interrupted) cout << "Toulbar2 best solution found: " << getOptimum() << endl;
	else cout << "Toulbar2 optimum found: " << getOptimum() << endl;
	for (unsigned int i=0; i < solution.size(); i++) {
	  cout << "Variable with wcspIndex " << i << " has value " << solution[i] << endl;  
	} 
  } else if (unsatisfiable) {
	cout << "Problem is unsatisfiable." << endl;
  } else {
	assert(interrupted);
	cout << "No solution found before being interrupted." << endl;
  }
}

/****************************************************************/
/* toulbar2 options                                             */
/****************************************************************/

void Toulbar2Solver::debug(const bool debug)
{
  ToulBar2::debug = debug;
}

void Toulbar2Solver::dumpWCSP(const int level, const char* problem)
{
  ToulBar2::dumpWCSP = level;
  ToulBar2::problemsaved_filename = to_string(problem);
}

void Toulbar2Solver::updateUb(const char* newUb)
{
    Cost ctUb = string2Cost(newUb);
    wcsp->updateUb(ctUb);
}

void Toulbar2Solver::nopre()
{
    ToulBar2::elimDegree = -1;
	ToulBar2::elimDegree_preprocessing = -1;
	ToulBar2::preprocessTernaryRPC = 0;
	ToulBar2::preprocessFunctional  = 0;
	ToulBar2::preprocessNary  = 0;
	ToulBar2::costfuncSeparate = false;
	ToulBar2::MSTDAC = false;
	ToulBar2::DEE = 0;
}

void Toulbar2Solver::btdMode(const int mode)
{
    ToulBar2::btdMode = mode;
}

void Toulbar2Solver::btdRootCluster(const int rCluster)
{
    ToulBar2::btdRootCluster = rCluster;
}

void Toulbar2Solver::btdSubTree(const int sTree)
{
    ToulBar2::btdSubTree = sTree;
}

void Toulbar2Solver::splitClusterMaxSize(const int size)
{
    ToulBar2::splitClusterMaxSize = size;
}

void Toulbar2Solver::boostingBTD(const bool boost)
{
    ToulBar2::boostingBTD = boost;
}
void Toulbar2Solver::maxSeparatorSize(const int size)
{
    ToulBar2::maxSeparatorSize = size;
}

void Toulbar2Solver::minProperVarSize(const int size)
{
    ToulBar2::minProperVarSize = size;
}

void Toulbar2Solver::writeSolution(const char* write)
{
    char* filename = new char[strlen(write)+1];
    strcpy(filename, write);
    ToulBar2::writeSolution = filename;
}

void Toulbar2Solver::showSolutions(const bool show)
{
    ToulBar2::showSolutions = show;
}

void Toulbar2Solver::allSolutions(const bool sol)
{
    ToulBar2::allSolutions = sol;
}

void Toulbar2Solver::approximateCountingBTD(const bool aproxim)
{
    ToulBar2::approximateCountingBTD = aproxim;
}

void Toulbar2Solver::binaryBranching(const bool boost)
{
    ToulBar2::binaryBranching = boost;
}

void Toulbar2Solver::staticVariableOrdering(const bool staticOrdering)
{
    ToulBar2::Static_variable_ordering = staticOrdering;
}

void Toulbar2Solver::lastConflict(const bool last)
{
    ToulBar2::lastConflict = last;
}

void Toulbar2Solver::dichotomicBranching(const int dicho)
{
    ToulBar2::dichotomicBranching = dicho;
}

void Toulbar2Solver::sortDomains(const bool sort)
{
    ToulBar2::sortDomains = sort;
}

void Toulbar2Solver::weightedDegree(const int wDegree)
{
    ToulBar2::weightedDegree = wDegree;
}

void Toulbar2Solver::weightedTightness(const int wTight)
{
    ToulBar2::weightedTightness = wTight;
}

void Toulbar2Solver::nbDecisionVars(const int nbDecision)
{
    ToulBar2::nbDecisionVars = nbDecision;
}

void Toulbar2Solver::elimDegree(const int degree)
{
    ToulBar2::elimDegree = degree;
}

void Toulbar2Solver::elimDegree_preprocessing(const int degree_prepoc)
{
    ToulBar2::elimDegree_preprocessing = degree_prepoc;
}

void Toulbar2Solver::elimSpaceMaxMB(const int size)
{
    ToulBar2::elimSpaceMaxMB = size;
}

void Toulbar2Solver::costfuncSeparate(const bool separate)
{
    ToulBar2::costfuncSeparate = separate;
}

void Toulbar2Solver::partialAssign(const char* certificate)
{
    solver->parse_solution(certificate);
}

void Toulbar2Solver::deadEndElimination(const int level)
{
  ToulBar2::DEE = level;
}

void Toulbar2Solver::vac(const int depth)
{
    ToulBar2::vac = depth;
}

void Toulbar2Solver::minsumDiffusion(const int min)
{
    ToulBar2::minsumDiffusion = min;
}

void Toulbar2Solver::costThreshold(const char* cost) 
{
    Cost ct = string2Cost(cost);
    ToulBar2::costThreshold = ct;
}

void Toulbar2Solver::costThresholdPre(const char* cost)
{
    Cost ct = string2Cost(cost);
    ToulBar2::costThresholdPre = ct;
}

void Toulbar2Solver::costMultiplier(const char* cost)
{
    double cm = atof(cost);
    ToulBar2::costMultiplier = cm;
}

void Toulbar2Solver::singletonConsistency(const bool singleCons)
{
    ToulBar2::singletonConsistency = singleCons;
}

void Toulbar2Solver::vacValueHeuristic(const bool vacVal)
{
    ToulBar2::vacValueHeuristic = vacVal;
}

void Toulbar2Solver::preprocessTernaryRPC(const int size)
{
    ToulBar2::preprocessTernaryRPC = size;
}

void Toulbar2Solver::preprocessFunctional(const int func)
{
    ToulBar2::preprocessFunctional = func;
}

void Toulbar2Solver::preprocessNary(const int maxnary)
{
    ToulBar2::preprocessNary = maxnary;
}

void Toulbar2Solver::QueueComplexity(const bool queue)
{
    ToulBar2::QueueComplexity = queue;
}

void Toulbar2Solver::lds(const int maxlds)
{
    ToulBar2::lds = maxlds;
}

void Toulbar2Solver::restart(const long maxrestarts)
{
    ToulBar2::restart = maxrestarts;
}

void Toulbar2Solver::lcLevel(const int level)
{
    LcLevelType lclevel = (LcLevelType) level;
    ToulBar2::LcLevel = lclevel;
}

void Toulbar2Solver::hbfs(const long hbfsgloballimit)
{
  ToulBar2::hbfs = 1;
  ToulBar2::hbfsGlobalLimit = hbfsgloballimit;
}

void Toulbar2Solver::hbfsAlpha(const long hbfsalpha)
{
  ToulBar2::hbfsAlpha = hbfsalpha;
}

void Toulbar2Solver::hbfsBeta(const long hbfsbeta)
{
  ToulBar2::hbfsBeta = hbfsbeta;
}

void Toulbar2Solver::hbfsOpenNodeLimit(const long openlimit)
{
  ToulBar2::hbfs = 1;
  if (ToulBar2::hbfsGlobalLimit==0) ToulBar2::hbfsGlobalLimit = 10000;
  ToulBar2::hbfsOpenNodeLimit = openlimit;
}

void Toulbar2Solver::variableEliminationOrdering(const int order)
{
  ToulBar2::varOrder = reinterpret_cast<char *>(abs(order));
}

void Toulbar2Solver::incop(const char* cmd)
{
  if (cmd == NULL || strlen(cmd) == 0) {
	ToulBar2::incop_cmd = "0 1 3 idwa 100000 cv v 0 200 1 0 0";
  } else {
    char* cmd_ = new char[strlen(cmd)+1];
    strcpy(cmd_, cmd);
	ToulBar2::incop_cmd = cmd_;
  }
}
