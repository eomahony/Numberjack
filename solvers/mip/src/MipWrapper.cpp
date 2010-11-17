
#include <iostream>

#include "MipWrapper.hpp"
#include <math.h>

using namespace std;

/**************************************************************
 ********************     EXPRESSION        *******************
 **************************************************************/

void MipWrapper_Expression::initialise(bool c){ 
  _lower = _upper = _coef = 0;
  _ident = nbj_ident = -1;
  _expr_encoding = NULL;
  _solver = NULL;
  _var = NULL;
  _continuous = c;
}

MipWrapper_Expression::MipWrapper_Expression(){
  initialise(true);
}

int MipWrapper_Expression::get_size(){ return _upper - _lower +1; }

int MipWrapper_Expression::get_max(){ return _upper; }

int MipWrapper_Expression::get_min(){ return _lower; }

void MipWrapper_Expression::leq(double value, MipWrapperSolver* solver){
  if(_expr_encoding != NULL){
    for(int i = _upper; i >= value && i >= _lower; ++i)
      if(_expr_encoding[i] != NULL){
	LinearConstraint *con = new LinearConstraint(0, 0);
	con->add_coef(_expr_encoding[i], 1);
	solver->_constraints.push_back(con);
      }
  } else if(value < _upper){
    LinearConstraint *con = new LinearConstraint(-INFINITY, value);
    con->add_coef(this, 1);
    solver->_constraints.push_back(con);
  }
}

void MipWrapper_Expression::lt(double value, MipWrapperSolver* solver){
  leq(value-1, solver);
}

void MipWrapper_Expression::geq(double value, MipWrapperSolver* solver){
  if(_expr_encoding != NULL){
    for(int i = _lower; i <= value && i <= _upper; ++i)
      if(_expr_encoding[i] != NULL){
	LinearConstraint *con = new LinearConstraint(0, 0);
	con->add_coef(_expr_encoding[i], 1);
	solver->_constraints.push_back(con);
      }
  } else if(value > _lower){
    LinearConstraint *con = new LinearConstraint(value, INFINITY);
    con->add_coef(this, 1);
    solver->_constraints.push_back(con);
  }
}

void MipWrapper_Expression::gt(double value, MipWrapperSolver* solver){
  geq(value+1, solver);
}

void MipWrapper_Expression::eq(double value, MipWrapperSolver* solver){
  // TODO: Check this
  if(_expr_encoding != NULL){
    if(value >= _lower && value <= _upper &&
       _expr_encoding[(int)value] != NULL){
      LinearConstraint *con = new LinearConstraint(1, 1);
      con->add_coef(_expr_encoding[(int) value], 1);
      solver->_constraints.push_back(con);
    }
  }  else {
    LinearConstraint *con = new LinearConstraint(value, value);
    con->add_coef(this, 1);
    solver->_constraints.push_back(con);
  }
}

void MipWrapper_Expression::neq(double value, MipWrapperSolver* solver){
  if(value >= _lower && value <= _upper)
  encode(solver);
  if(value >= _lower && value <= _upper && _expr_encoding[(int)value] != NULL){
    LinearConstraint *con = new LinearConstraint(0, 0);
    con->add_coef(_expr_encoding[(int)value], 1);
    solver->_constraints.push_back(con);
  }
}

void MipWrapper_Expression::display(){
  std::cout << _coef << "*X" << _ident << " ";
}

MipWrapper_Expression::~MipWrapper_Expression(){}

void MipWrapper_Expression::encode(MipWrapperSolver *solver){
  if(!_expr_encoding) {
    DBG("Creating an encoding %s\n", "");
    
    _expr_encoding = new MipWrapper_Expression*[get_size()];
    _expr_encoding -= (int)_lower;
    LinearConstraint *con = new LinearConstraint(1, 1);
    for(int i = _lower; i <= _upper; ++i){
      _expr_encoding[i] = new MipWrapper_IntVar(0, 1);
      _expr_encoding[i]->add(solver, false);
      con->add_coef(_expr_encoding[i], 1);
    }
    solver->_constraints.push_back(con);
  }
}

bool MipWrapper_Expression::has_been_added() const{
  return _solver != NULL;
}

MipWrapper_Expression* MipWrapper_Expression::add(MipWrapperSolver *solver,
						  bool top_level){
  if(!this->has_been_added()){
    _solver = solver;
    DBG("add variable [%f..%f]\n", _lower, _upper);
    this->_ident = solver->var_counter++;
  }
  return this;
}

double MipWrapper_Expression::get_whatever_value(){ return 0; }

LINEAR_ARG* MipWrapper_Expression::for_linear(){
  LINEAR_ARG *larg = (LINEAR_ARG*) calloc(1, sizeof(LINEAR_ARG));
  larg->expr = this;
  larg->coef = 1;
  larg->offset = 0;
  return larg;
}

int MipWrapper_Expression::for_linear_size(){ return 1; }

/**
 * IntVar stuff
 */

MipWrapper_IntVar::MipWrapper_IntVar(const int lb,
				     const int ub){
  initialise(false);
  _upper = (double)ub;
  _lower = (double)lb;
  _has_holes_in_domain = false;
  DBG("creating integer variable [%d..%d]\n", lb, ub);
}

MipWrapper_IntVar::MipWrapper_IntVar(const int lb,
				     const int ub,
				     const int ident){
  initialise(false);
  _upper = (double)ub;
  _lower = (double)lb;
  _has_holes_in_domain = false;
  nbj_ident = ident;
  DBG("creating integer variable [%d..%d]\n", lb, ub);
}

MipWrapper_IntVar::MipWrapper_IntVar(MipWrapperIntArray& values,
				     const int ident){
  initialise(false);
  
  _values = values;
  _has_holes_in_domain = true;
  nbj_ident = ident;
  
  // Ensure that bounds are set up properly  
  _upper = _values.get_item(0);
  _lower = _values.get_item(0);
  for(int i = 0; i < _values.size(); ++i){
    if(_values.get_item(i) > _upper) _upper = _values.get_item(i);
    if(_values.get_item(i) < _lower) _lower = _values.get_item(i);
  }
  DBG("Creating variable with holes in the domain, size:%d \n", _values.size()); 
}

MipWrapper_Expression* MipWrapper_IntVar::add(MipWrapperSolver* solver,
					      bool top_level){
  MipWrapper_Expression::add(solver, top_level);
  if(_has_holes_in_domain) encode(solver);
  return this;
}

void MipWrapper_IntVar::encode(MipWrapperSolver* solver){
  if(!_expr_encoding) {
    if(_has_holes_in_domain){
      DBG("Making an encoding for a variable with holes %s\n", "");

      _expr_encoding = new MipWrapper_Expression*[get_size()];
      for(int i = 0; i < get_size(); ++i) _expr_encoding[i] = NULL;
      _expr_encoding -= (int)_lower;
      
      LinearConstraint *con = new LinearConstraint(1, 1);
      for(int i = 0; i < _values.size(); ++i){
	_expr_encoding[_values.get_item(i)] = new MipWrapper_IntVar(0, 1);
	_expr_encoding[_values.get_item(i)]->add(solver, false);
	con->add_coef(_expr_encoding[i], 1);
      }
      _solver->_constraints.push_back(con);

    } else MipWrapper_Expression::encode(solver);
  }
}

int MipWrapper_IntVar::get_value(){
  if(_expr_encoding == NULL) return (int)(round(_solver->get_value(_var)));
  double res = 0;
  for(int i = _lower; i <= _upper; ++i)
    if(_expr_encoding[i] != NULL)
      res += i * _expr_encoding[i]->get_whatever_value();
  return res;
}

double MipWrapper_IntVar::get_whatever_value(){return get_value();}

/**
 * FloatVar stuff
 */

MipWrapper_FloatVar::MipWrapper_FloatVar(){
  initialise(true);
  _upper = 0;
  _lower = 0;
  DBG("Crearing continuous variable (%f .. %f)\n", _upper, _lower);
}

MipWrapper_FloatVar::MipWrapper_FloatVar(const double lb, const double ub){
  initialise(true);
  _upper = ub;
  _lower = lb;
  DBG("Crearing continuous variable (%f .. %f)\n", _upper, _lower);
}

MipWrapper_FloatVar::MipWrapper_FloatVar(const double lb,
					 const double ub,
					 const int ident){
  initialise(true);
  _upper = ub;
  _lower = lb;
  nbj_ident = ident;
  DBG("Crearing continuous variable (%f .. %f)\n", _upper, _lower);
}

double MipWrapper_FloatVar::get_value(){
  if(_expr_encoding == NULL) return _solver->get_value(_var);  
  double res = 0;
  for(int i = _lower; i <= _upper; ++i)
    if(_expr_encoding[i] != NULL)
      res += i * _expr_encoding[i]->get_whatever_value();
  return res;
}

double MipWrapper_FloatVar::get_whatever_value(){ return get_value(); }

/**
 * Constraint stuff
 */

MipWrapper_add::MipWrapper_add(MipWrapper_Expression *arg1,
			       MipWrapper_Expression *arg2)
  : MipWrapper_Sum(){
  addVar(arg1);
  addVar(arg2);
  addWeight(1);
  addWeight(1);
  initialise();
}

MipWrapper_add::MipWrapper_add(MipWrapper_Expression *arg1,
			       const int arg2)
  : MipWrapper_Sum(){
  addVar(arg1);
  addWeight(1);
  set_rhs(-arg2);
  initialise();
}

MipWrapper_add::~MipWrapper_add(){}

MipWrapper_sub::MipWrapper_sub(MipWrapper_Expression *arg1, const int arg2)
  : MipWrapper_Sum(){
  addVar(arg1);
  addWeight(1);
  set_rhs(arg2);
  initialise();
}

MipWrapper_sub::MipWrapper_sub(MipWrapper_Expression *arg1,
			       MipWrapper_Expression *arg2)
  : MipWrapper_Sum(){
  addVar(arg1);
  addVar(arg2);
  addWeight(1);
  addWeight(-1);
  initialise();
}

MipWrapper_sub::~MipWrapper_sub() {}

MipWrapper_AllDiff::MipWrapper_AllDiff( MipWrapper_Expression* arg1,
				       MipWrapper_Expression* arg2 ) 
  : MipWrapper_Flow(){
  addVar(arg1);
  addVar(arg2);
  initbounds();
  DBG("Creating an AllDiff constraint %s\n", "");
  std::fill(this->card_ub+this->min_val, this->card_ub+this->max_val+1, 1);
}

MipWrapper_AllDiff::MipWrapper_AllDiff( MipWrapperExpArray& vars ) 
  : MipWrapper_Flow(vars) {
  DBG("Creating an AllDiff constraint %s\n", "");
  std::fill(this->card_ub+this->min_val, this->card_ub+this->max_val+1, 1);
}

MipWrapper_AllDiff::~MipWrapper_AllDiff(){}

MipWrapper_Gcc::MipWrapper_Gcc(MipWrapperExpArray& vars,
			       MipWrapperIntArray& vals,
			       MipWrapperIntArray& lb,
			       MipWrapperIntArray& ub ) 
  : MipWrapper_Flow(vars){
  DBG("Creating an GCC constraint %s\n", "");
  int j, n=vals.size();
  for(int i=0; i<n; ++i){
      j = vals.get_item(i);
      this->card_lb[j] = lb.get_item(i);
      this->card_ub[j] = ub.get_item(i);
  }
}

MipWrapper_Gcc::~MipWrapper_Gcc(){}

MipWrapper_Flow::MipWrapper_Flow(MipWrapperExpArray& vars)
  : MipWrapper_FloatVar() {
  _vars = vars;
  initbounds();
}

MipWrapper_Flow::MipWrapper_Flow()
  : MipWrapper_FloatVar() {}

void MipWrapper_Flow::initbounds(){
  int i, lb, ub;
  max_val = -INT_MAX;
  min_val = INT_MAX;
  
  for(i = 0; i < _vars.size(); ++i) {
    lb = (int)(this->_vars.get_item(i)->_lower);
    ub = (int)(this->_vars.get_item(i)->_upper);
    
    if(lb < min_val) min_val = lb;
    if(ub > max_val) max_val = ub;
  }

  card_lb = new int[max_val-min_val+1];
  card_ub = new int[max_val-min_val+1];

  card_lb -= min_val;
  card_ub -= min_val;

  std::fill(card_lb+min_val, card_lb+max_val+1, 0);
  std::fill(card_ub+min_val, card_ub+max_val+1, _vars.size());
}

MipWrapper_Flow::~MipWrapper_Flow(){}

void MipWrapper_Flow::addVar(MipWrapper_Expression* v){ _vars.add(v); }

MipWrapper_Expression* MipWrapper_Flow::add(MipWrapperSolver *solver,
					    bool top_level){
  if(!has_been_added()){
    _solver = solver;
    DBG("Adding in Flow %s\n", "");

    if(top_level){
      
      for(int i = 0; i < _vars.size(); ++i) {
	_vars.get_item(i)->add(solver, false);
	_vars.get_item(i)->encode(solver);
      }
      
      for(int i = min_val; i <= max_val; ++i) {
	LinearConstraint *con = new LinearConstraint(card_lb[i], card_ub[i]);
	for(int j = 0; j < _vars.size(); ++j)
	  if( i <= (int)(_vars.get_item(j)->_upper) &&
	      i >= (int)(_vars.get_item(j)->_lower) &&
	      _vars.get_item(j)->_expr_encoding[i] != NULL)
	      con->add_coef(_vars.get_item(j)->_expr_encoding[i], 1);
	solver->_constraints.push_back(con);     
      }
      
    } else {
      std::cout << "Warning Flow constraint supported only on top level"
	        << std::endl;
      exit(1);
    }
  }
  return this;
}


MipWrapper_Sum::MipWrapper_Sum(MipWrapperExpArray& vars, 
		   MipWrapperIntArray& weights, 
		   const int offset)
  : MipWrapper_FloatVar(){
  _offset = offset;
  _vars = vars;
  _weights = weights;
  initialise();
}

MipWrapper_Sum::MipWrapper_Sum(MipWrapper_Expression *arg1, 
		   MipWrapper_Expression *arg2, 
		   MipWrapperIntArray& w, 
		   const int offset)
  : MipWrapper_FloatVar() {
  _offset = offset;
  _vars.add(arg1);
  _vars.add(arg2);
  _weights = w;
  initialise();
}

MipWrapper_Sum::MipWrapper_Sum(MipWrapper_Expression *arg, 
		   MipWrapperIntArray& w, 
		   const int offset)
  : MipWrapper_FloatVar() {
  _offset = offset;
  _vars.add(arg);
  _weights = w;
  initialise();
}

MipWrapper_Sum::MipWrapper_Sum()
  : MipWrapper_FloatVar(){
  _offset = 0;
}

void MipWrapper_Sum::initialise() {
  DBG("creating sum: size of params: %d\n", _vars.size());
  for(int i = 0; i < this->_vars.size(); ++i){
    int weight = this->_weights.get_item(i);
    if( weight > 0 ) {
      _lower += (weight * _vars.get_item(i)->_lower);
      _upper += (weight * _vars.get_item(i)->_upper);
    } else {
      _upper += (weight * _vars.get_item(i)->_lower);
      _lower += (weight * _vars.get_item(i)->_upper);
    }
  }
  _lower += _offset;
  _upper += _offset;
  DBG("Sum has domain [%f..%f]\n", _lower, _upper);    
}

MipWrapper_Sum::~MipWrapper_Sum(){}

void MipWrapper_Sum::addVar(MipWrapper_Expression* v){ _vars.add(v); }

void MipWrapper_Sum::addWeight(const int w){ _weights.add(w);}

void MipWrapper_Sum::set_rhs(const int k){ _offset = k; }

MipWrapper_Expression* MipWrapper_Sum::add(MipWrapperSolver *solver,
					   bool top_level){
  if(!this->has_been_added()){
    _solver = solver;
    DBG("Adding sum constraint %s\n", "");
    for(int i = 0; i < _vars.size(); ++i)
      _vars.set_item(i, _vars.get_item(i)->add(solver, false));
    
    if(top_level){
      std::cout << "Warning SUM constraint on top level not supported"
		<< std::endl;
      exit(1);
    }
  }
  return this;
}

double MipWrapper_Sum::get_value(){
    double res = 0;
    for(int i = 0; i < _vars.size(); ++i)
        res += _vars.get_item(i)->get_whatever_value() * _weights.get_item(i);
    return res;
}

LINEAR_ARG* MipWrapper_Sum::for_linear(){
  LINEAR_ARG* largs = (LINEAR_ARG*) calloc(_vars.size(), sizeof(LINEAR_ARG));
  for(int i = 0; i < _vars.size(); ++i){
    largs[i].expr = _vars.get_item(i);
    largs[i].coef = _weights.get_item(i);
    largs[i].offset = _offset;
    //std::cout << "Weight:" << largs[i].coef << std::endl;
  }
  return largs;
}

int MipWrapper_Sum::for_linear_size(){ return _vars.size(); }

void MipWrapper_Sum::leq(double value, MipWrapperSolver* solver){
  LinearConstraint *con = new LinearConstraint(-INFINITY, value);
  con->add_coef(this, 1);
  solver->_constraints.push_back(con);
}

void MipWrapper_Sum::geq(double value, MipWrapperSolver* solver){
  LinearConstraint *con = new LinearConstraint(value, INFINITY);
  con->add_coef(this, 1);
  solver->_constraints.push_back(con);
}

void MipWrapper_Sum::eq( double value, MipWrapperSolver* solver){
  LinearConstraint *con = new LinearConstraint(value, value);
  con->add_coef(this, 1);
  solver->_constraints.push_back(con);
}

void MipWrapper_Sum::neq(double value, MipWrapperSolver* solver){
  MipWrapper_FloatVar *var = new MipWrapper_FloatVar(_lower, _upper);
  var->add(solver, false);
  
  LinearConstraint *con = new LinearConstraint(0, 0);
  con->add_coef(this, 1);
  con->add_coef(var, -1);
  solver->_constraints.push_back(con);
  
  var->neq(value, solver);
}

void MipWrapper_Sum::encode(MipWrapperSolver* solver){
  MipWrapper_Expression::encode(solver);
  LinearConstraint *con = new LinearConstraint(0, 0);
  for(int i = _lower; i <= _upper; ++i){
    con->add_coef(_expr_encoding[i], -i);
  }
  con->add_coef(this, 1, false);
  solver->_constraints.push_back(con);
}

/* Binary operators */

MipWrapper_binop::MipWrapper_binop(MipWrapper_Expression *var1,
				   MipWrapper_Expression *var2)
  : MipWrapper_FloatVar(){
  _vars[0] = var1;
  _vars[1] = var2;
  _is_proper_coef = false;
  _upper = 1.0;
  _lower = 0.0;
  DBG("Crearing binary operator %s\n", "");
}

MipWrapper_binop::MipWrapper_binop(MipWrapper_Expression *var1, double rhs)
  : MipWrapper_FloatVar(){
  _vars[0] = var1;
  _vars[1] = NULL;
  _rhs = rhs;
  _is_proper_coef = true;
  _upper = 1.0;
  _lower = 0.0;
  DBG("Crearing binary operator %s\n", "");
}

MipWrapper_binop::~MipWrapper_binop(){}

MipWrapper_NoOverlap::MipWrapper_NoOverlap(MipWrapper_Expression *var1,
		      MipWrapper_Expression *var2,
		      int d1, int d2):
  MipWrapper_binop(var1, var2){
  _coefs = *(new MipWrapperIntArray());
  _coefs.add(d1);
  _coefs.add(d2);  
}

MipWrapper_NoOverlap::~MipWrapper_NoOverlap(){}

MipWrapper_Expression* MipWrapper_NoOverlap::add(MipWrapperSolver *solver,
						 bool top_level){
  if(!has_been_added()){
    _solver = solver;
   _vars[0] = _vars[0]->add(solver, false);
   _vars[1] = _vars[1]->add(solver, false);
   
   if(top_level){
     
     MipWrapperIntArray *arr1 = new MipWrapperIntArray();
     arr1->add(_coefs.get_item(0));
     solver->add_int_array(arr1);
     
     MipWrapperIntArray *arr2 = new MipWrapperIntArray();
     arr2->add(_coefs.get_item(1));
     solver->add_int_array(arr2);
     
     MipWrapper_Expression *prec1 = new MipWrapper_Precedence(_vars[0],
  							     _vars[1], *arr1);
     solver->add_expr(prec1);
     
     MipWrapper_Expression *prec2 = new MipWrapper_Precedence(_vars[1],
  							     _vars[0], *arr2);
     solver->add_expr(prec2);
     
     MipWrapper_Expression *orexp = new MipWrapper_or(prec1, prec2);
     solver->add_expr(orexp);
     orexp->add(solver, true);
     
   } else {
     std::cout << "No support for NoOverLap operator not in top level"
	       << std::endl;
     exit(1);
   }
  }
  return this;
}

MipWrapper_Precedence::MipWrapper_Precedence(MipWrapper_Expression *var1,
					     MipWrapper_Expression *var2,
					     MipWrapperIntArray &coefs):
  MipWrapper_binop(var1, var2), _coefs(coefs){}

MipWrapper_Precedence::~MipWrapper_Precedence(){}

MipWrapper_Expression* MipWrapper_Precedence::add(MipWrapperSolver *solver,
						  bool top_level){
  
  if(!has_been_added()){
    _solver = solver;
   _vars[0] = _vars[0]->add(solver, false);
   _vars[1] = _vars[1]->add(solver, false);
   
   if(top_level){
    LinearConstraint *con = new LinearConstraint(-_coefs.get_item(0),INFINITY);
     con->add_coef(_vars[0], 1);
     con->add_coef(_vars[1], -1);
     solver->_constraints.push_back(con);
     return this;
   } else {
    
    
    
     std::cout << "Precedence not at top level not supported yet" << std::endl;
     exit(1);
   }
  }
  return this;
}

MipWrapper_eq::MipWrapper_eq(MipWrapper_Expression *var1,
			     MipWrapper_Expression *var2)
  : MipWrapper_binop(var1,var2){ DBG("creating equality constraint", NULL); }

MipWrapper_eq::MipWrapper_eq(MipWrapper_Expression *var1, double rhs)
  : MipWrapper_binop(var1,rhs){ DBG("creating equality constraint", NULL); }

MipWrapper_eq::~MipWrapper_eq(){}

MipWrapper_Expression* MipWrapper_eq::add(MipWrapperSolver *solver,
					  bool top_level){
  if(!has_been_added()){
    _solver = solver;
    DBG("adding equality constraint", NULL);
    _vars[0] = _vars[0]->add(solver, false);
    
    if(top_level){
      if(_is_proper_coef){
	_vars[0]->eq(_rhs, solver);
	//LinearConstraint *con = new LinearConstraint(_rhs, _rhs);
	//con->add_coef(_vars[0], 1);
	//solver->_constraints.push_back(con);
      } else {
	_vars[1] = _vars[1]->add(solver, false);
	LinearConstraint *con = new LinearConstraint(0, 0);
	con->add_coef(_vars[0], 1);
	con->add_coef(_vars[1], -1);
	solver->_constraints.push_back(con);
      }
    } else {
      _vars[0]->encode(solver);
      if(_is_proper_coef){
	if(_vars[0]->_expr_encoding[(int)_rhs] != NULL){
	  _repr = _vars[0]->_expr_encoding[(int)_rhs];
	  return _vars[0]->_expr_encoding[(int)_rhs];
	}
	else return new MipWrapper_IntVar(0, 0);
      } else {
	MipWrapper_Expression *bvar = new MipWrapper_IntVar(0, 1);
	bvar = bvar->add(solver, false);
	_vars[1] = _vars[1]->add(solver, false);
	_vars[1]->encode(solver);
	
	int lb = (int) std::max(_vars[0]->_lower, _vars[1]->_lower);
	int ub = (int) std::min(_vars[0]->_upper, _vars[1]->_upper);
	
	// \forall_{i \in D_a \cup D_b} -1 <= A_i + B_i - C <= 1
	for(int i = lb; i <= ub; ++lb){
	  if(_vars[0]->_expr_encoding[i] != NULL &&
	     _vars[1]->_expr_encoding[i]){
	     LinearConstraint *con = new LinearConstraint(-1, 1);
	     con->add_coef(bvar, -1);
	     con->add_coef(_vars[0]->_expr_encoding[i], 1);
	     con->add_coef(_vars[1]->_expr_encoding[i], 1);
	     solver->_constraints.push_back(con);
	  }
	}
	
	double range_max = std::max( _vars[0]->_upper - _vars[1]->_lower,
				  _vars[1]->_upper - _vars[0]->_lower);
	MipWrapper_FloatVar *y = new MipWrapper_FloatVar(0, range_max);
	
	LinearConstraint *con = new LinearConstraint(1, INFINITY);
	con->add_coef(_vars[0], 1);
	con->add_coef(_vars[1], -1);
	con->add_coef(y, 1);
	solver->_constraints.push_back(con);
	
	con = new LinearConstraint(1, INFINITY);
	con->add_coef(_vars[0], -1);
	con->add_coef(_vars[1], 1);
	con->add_coef(y, 1);
	solver->_constraints.push_back(con);
	
	// RM <= C*RM + Y <= RM
	con = new LinearConstraint(range_max, range_max);
	con->add_coef(bvar, range_max);
	con->add_coef(y, 1);
	solver->_constraints.push_back(con);
	
	_repr = bvar;
	return bvar;
      }
    }
  }
  if(top_level) return this;
  return _repr;
}

/* In equality operator */

MipWrapper_ne::MipWrapper_ne(MipWrapper_Expression *var1,
			     MipWrapper_Expression *var2)
  : MipWrapper_binop(var1,var2){}

MipWrapper_ne::MipWrapper_ne(MipWrapper_Expression *var1, double rhs)
  : MipWrapper_binop(var1,rhs){}

MipWrapper_ne::~MipWrapper_ne(){}

MipWrapper_Expression* MipWrapper_ne::add(MipWrapperSolver *solver,
					  bool top_level){
  if(!this->has_been_added()){
    _solver = solver;
    DBG("adding not equal constraint %s\n", "");
    _vars[0] = _vars[0]->add(solver, false);
    _vars[0]->encode(solver);
    if(top_level){
      if(_is_proper_coef){
	_vars[0]->neq(_rhs, solver);
      } else {
        _vars[1] = _vars[1]->add(solver, false);
        _vars[1]->encode(solver);
	
        int lb = std::max(this->_vars[0]->_lower, this->_vars[1]->_lower);
	int ub = std::min(this->_vars[0]->_upper, this->_vars[1]->_upper);
	
	for(int i = lb; i <= ub; ++i)
	  if(_vars[0]->_expr_encoding[i] != NULL &&
	     _vars[1]->_expr_encoding[i] != NULL){
	    LinearConstraint *con = new LinearConstraint(0, 1);
	    con->add_coef(_vars[0]->_expr_encoding[i], 1);
	    con->add_coef(_vars[1]->_expr_encoding[i], 1);
	    solver->_constraints.push_back(con);
	  }
      }
      return this;
    } else {
      DBG("Reifing not equal constraint %s\n", "");
      
      if(_is_proper_coef){
	
	//TODO:
	
      } else {	
	_vars[1] = _vars[1]->add(solver, false);
	_vars[0]->encode(solver);
	_vars[1]->encode(solver);
      
	MipWrapper_Expression *c = new MipWrapper_IntVar(0, 1);
	c->add(solver, false);
	
	//TODO:
	std::cout << "Warning ne reif not supported at the moment"
		  << std::endl;
	exit(1);
      }
    }
  }
  if(top_level) return this;
  return _repr;
}

/* Leq operator */

MipWrapper_le::MipWrapper_le(MipWrapper_Expression *var1,
			     MipWrapper_Expression *var2)
  : MipWrapper_binop(var1,var2){ DBG("Crearing le constraint  %s\n", ""); }

MipWrapper_le::MipWrapper_le(MipWrapper_Expression *var1, double rhs)
  : MipWrapper_binop(var1,rhs){ DBG("Crearing le constraint  %s\n", ""); }

MipWrapper_le::~MipWrapper_le(){}

MipWrapper_Expression* MipWrapper_le::add(MipWrapperSolver *solver,
					  bool top_level){
  if(!has_been_added()){
    _solver = solver;
    DBG("adding le constraint %s\n", "");
    _vars[0] = _vars[0]->add(solver, false);
    if(top_level){
      if(_is_proper_coef){
	_vars[0]->leq(_rhs, solver);
      } else {
	_vars[1] = _vars[1]->add(solver, false);
	LinearConstraint *con = new LinearConstraint(0, INFINITY);
	con->add_coef(_vars[0], -1);
	con->add_coef(_vars[1], 1);
	solver->_constraints.push_back(con);
      }
    } else {
      MipWrapper_Expression *C = new MipWrapper_IntVar(0, 1);
      C = C->add(solver, false);
      
      if(_is_proper_coef){
	double M = std::max( _vars[0]->_upper - _rhs,
			     _rhs - _vars[0]->_lower);
	
	//TODO: Check this
	LinearConstraint *con1 = new LinearConstraint(-M-_rhs, INFINITY);
	con1->add_coef(_vars[0], -1);
	con1->add_coef(C, -M);
	solver->_constraints.push_back(con1);
	
	LinearConstraint *con2 = new LinearConstraint(1+_rhs, INFINITY);
	con2->add_coef(_vars[0], 1);
	con2->add_coef(C, M+1);
	solver->_constraints.push_back(con2);
	
      } else {
	_vars[1] = _vars[1]->add(solver, false);
	double M = std::max( _vars[0]->_upper - _vars[1]->_lower,
			     _vars[1]->_upper - _vars[0]->_lower);
	
	LinearConstraint *con1 = new LinearConstraint(-M, INFINITY);
	con1->add_coef(_vars[1], 1);
	con1->add_coef(_vars[0], -1);
	con1->add_coef(C, -M);
	solver->_constraints.push_back(con1);
	
	LinearConstraint *con2 = new LinearConstraint(1, INFINITY);
	con2->add_coef(_vars[0], 1);
	con2->add_coef(_vars[1], -1);
	con2->add_coef(C, M+1);
	solver->_constraints.push_back(con2);
      }
      _repr = C;
      return C;
    }
  }
  if(top_level) return this;
  return _repr;
}

/* Geq operator */
MipWrapper_ge::MipWrapper_ge(MipWrapper_Expression *var1,
			     MipWrapper_Expression *var2)
  : MipWrapper_binop(var1,var2){ DBG("Crearing ge constraint  %s\n", ""); }

MipWrapper_ge::MipWrapper_ge(MipWrapper_Expression *var1,
			     double rhs)
  : MipWrapper_binop(var1,rhs){ DBG("Crearing ge constraint  %s\n", ""); }

MipWrapper_ge::~MipWrapper_ge(){}

MipWrapper_Expression* MipWrapper_ge::add(MipWrapperSolver *solver,
					  bool top_level){ 
  if(!has_been_added()){
    _solver = solver;
    DBG("Adding ge constraint  %s\n", "");
    _vars[0] = _vars[0]->add(solver, false);
    if(top_level){
      if(this->_is_proper_coef){
	_vars[0]->geq(_rhs, solver);
      } else {
	_vars[1] = _vars[1]->add(solver, false);
	LinearConstraint *con = new LinearConstraint(0, INFINITY);
	con->add_coef(_vars[0], 1);
	con->add_coef(_vars[1], -1);
	solver->_constraints.push_back(con);
      }
      return this;
    } else {
      MipWrapper_Expression *C = new MipWrapper_IntVar(0, 1);
      C = C->add(solver, false);
      
      if(_is_proper_coef){
	double M = std::max( _vars[0]->_upper - _rhs,
			     _rhs - _vars[0]->_lower);
	
	//TODO: Check this
	LinearConstraint *con1 = new LinearConstraint(-M+_rhs, INFINITY);
	con1->add_coef(_vars[0], 1);
	con1->add_coef(C, -M);
	solver->_constraints.push_back(con1);
	
	LinearConstraint *con2 = new LinearConstraint(1-_rhs, INFINITY);
	con2->add_coef(_vars[0], -1);
	con2->add_coef(C, M+1);
	solver->_constraints.push_back(con2);
	
      } else {
	_vars[1] = _vars[1]->add(solver, false);
	double M = std::max( _vars[0]->_upper - _vars[1]->_lower,
			     _vars[1]->_upper - _vars[0]->_lower);
	
	LinearConstraint *con1 = new LinearConstraint(-M, INFINITY);
	con1->add_coef(_vars[0], 1);
	con1->add_coef(_vars[1], -1);
	con1->add_coef(C, -M);
	solver->_constraints.push_back(con1);
	
	LinearConstraint *con2 = new LinearConstraint(1, INFINITY);
	con2->add_coef(_vars[1], 1);
	con2->add_coef(_vars[0], -1);
	con2->add_coef(C, M+1);
	solver->_constraints.push_back(con2);
      }
      _repr = C;
      return C;
    }
  }
  if(top_level) return this;
  return _repr;
}

/* Lt object */
MipWrapper_lt::MipWrapper_lt(MipWrapper_Expression *var1,
			     MipWrapper_Expression *var2)
  : MipWrapper_binop(var1,var2){ DBG("Crearing lt constraint %s\n", ""); }

MipWrapper_lt::MipWrapper_lt(MipWrapper_Expression *var1,
			     double rhs)
  : MipWrapper_binop(var1,rhs){ DBG("Crearing lt constraint %s\n", ""); }

MipWrapper_lt::~MipWrapper_lt(){}

MipWrapper_Expression* MipWrapper_lt::add(MipWrapperSolver *solver,
					  bool top_level){
  if(!has_been_added()){
    _solver = solver;
    DBG("Adding lt constraint %s\n", "");
    _vars[0] = _vars[0]->add(solver, false);
    if(top_level){
      if(_is_proper_coef){
	_vars[0]->lt(_rhs, solver);
      } else {
	_vars[1] = _vars[1]->add(solver, false);
	LinearConstraint *con = new LinearConstraint(1, INFINITY);
	con->add_coef(_vars[0], -1);
	con->add_coef(_vars[1], 1);
	solver->_constraints.push_back(con);
      }
    } else {
      MipWrapper_Expression *C = new MipWrapper_IntVar(0, 1);
      C = C->add(solver, false);
      
      if(_is_proper_coef){
	double M = std::max( _vars[0]->_upper - _rhs,
			     _rhs - _vars[0]->_lower);
	
	//TODO: Check this
	LinearConstraint *con1 = new LinearConstraint(-M-_rhs, INFINITY);
	con1->add_coef(_vars[0], -1);
	con1->add_coef(C, -(M+1));
	solver->_constraints.push_back(con1);
	
	LinearConstraint *con2 = new LinearConstraint(0+_rhs, INFINITY);
	con2->add_coef(_vars[0], 1);
	con2->add_coef(C, M+1);
	solver->_constraints.push_back(con2);
	
      } else {
	_vars[1] = _vars[1]->add(solver, false);
	double M = std::max( _vars[0]->_upper - _vars[1]->_lower,
			     _vars[1]->_upper - _vars[0]->_lower);
	
	LinearConstraint *con1 = new LinearConstraint(-M, INFINITY);
	con1->add_coef(_vars[1], 1);
	con1->add_coef(_vars[0], -1);
	con1->add_coef(C, -(M+1));
	solver->_constraints.push_back(con1);
	
	LinearConstraint *con2 = new LinearConstraint(0, INFINITY);
	con2->add_coef(_vars[1], -1);
	con2->add_coef(_vars[0], 1);
	con2->add_coef(C, M+1);
	solver->_constraints.push_back(con2);
      }
      _repr = C;
      return C;
    }
  }
  if(top_level) return this;
  return _repr;
}

/* Gt object */
MipWrapper_gt::MipWrapper_gt(MipWrapper_Expression *var1,
			     MipWrapper_Expression *var2)
  : MipWrapper_binop(var1,var2){ DBG("Crearing gt constraint %s\n", ""); }

MipWrapper_gt::MipWrapper_gt(MipWrapper_Expression *var1,
			     double rhs)
  : MipWrapper_binop(var1,rhs){ DBG("Crearing gt constraint %s\n", ""); }

MipWrapper_gt::~MipWrapper_gt(){}

MipWrapper_Expression* MipWrapper_gt::add(MipWrapperSolver *solver,
					  bool top_level){
  if(!has_been_added()){
    _solver = solver;
    DBG("Adding gt constraint %s\n", ""); 
    _vars[0] = _vars[0]->add(solver, false);
    if(top_level){
      if(_is_proper_coef){
	_vars[0]->gt(_rhs, solver);
      } else {
	_vars[1] = _vars[1]->add(solver, false);
	LinearConstraint *con = new LinearConstraint(1, INFINITY);
	con->add_coef(_vars[0], 1);
	con->add_coef(_vars[1], -1);
	solver->_constraints.push_back(con);
      }
    } else {
      MipWrapper_Expression *C = new MipWrapper_IntVar(0, 1);
      C = C->add(solver, false);
      
      if(_is_proper_coef){
	double M = std::max( _vars[0]->_upper - _rhs,
			     _rhs - _vars[0]->_lower);
	
	//TODO: Check this
	LinearConstraint *con1 = new LinearConstraint(-M+_rhs, INFINITY);
	con1->add_coef(_vars[0], 1);
	con1->add_coef(C, -(M+1));
	solver->_constraints.push_back(con1);
	
	LinearConstraint *con2 = new LinearConstraint(0-_rhs, INFINITY);
	con2->add_coef(_vars[0], -1);
	con2->add_coef(C, M+1);
	solver->_constraints.push_back(con2);
	
      } else {
	_vars[1] = _vars[1]->add(solver, false);
	double M = std::max( _vars[0]->_upper - _vars[1]->_lower,
			     _vars[1]->_upper - _vars[0]->_lower);
	
	LinearConstraint *con1 = new LinearConstraint(-M, INFINITY);
	con1->add_coef(_vars[0], 1);
	con1->add_coef(_vars[1], -1);
	con1->add_coef(C, -(M+1));
	solver->_constraints.push_back(con1);
	
	LinearConstraint *con2 = new LinearConstraint(0, INFINITY);
	con2->add_coef(_vars[0], -1);
	con2->add_coef(_vars[1], 1);
	con2->add_coef(C, M+1);
	solver->_constraints.push_back(con2);
      }
      _repr = C;
      return C;
    }
  }
  if(top_level) return this;
  return _repr;
}

/*
 *
 * Boolean operator stuff
 *
 */
MipWrapper_and::MipWrapper_and(MipWrapper_Expression *var1,
		    MipWrapper_Expression *var2):
  MipWrapper_binop(var1, var2){ DBG("Creating and constraint %s\n", ""); }

MipWrapper_and::~MipWrapper_and(){}

MipWrapper_Expression* MipWrapper_and::add(MipWrapperSolver *solver,
					   bool top_level){

  if(!has_been_added()){
    _solver = solver;
    DBG("Adding in and constraint  %s\n", "");
    _vars[0] = _vars[0]->add(solver, false);
    _vars[1] = _vars[1]->add(solver, false);
    if(top_level){
      _vars[0]->_lower = 1;
      _vars[1]->_lower = 1;
    } else {
      MipWrapper_Expression *C = new MipWrapper_IntVar(0, 1);
      C = C->add(solver, false);
      (new MipWrapper_le(C, _vars[0]))->add(solver, true);
      (new MipWrapper_le(C, _vars[1]))->add(solver, true);
       
      LinearConstraint *con = new LinearConstraint(1, INFINITY);
      con->add_coef(_vars[0], 1);
      con->add_coef(_vars[1], 1);
      con->add_coef(C, -1);
      solver->_constraints.push_back(con);
      
      _repr = C;
      return C;
    }
  }
  if(top_level) return this;
  return _repr;
}

MipWrapper_or::MipWrapper_or(MipWrapper_Expression *var1,
			     MipWrapper_Expression *var2):
  MipWrapper_binop(var1, var2){ DBG("Creating or constraint %s\n", ""); }

MipWrapper_or::~MipWrapper_or(){}

MipWrapper_Expression* MipWrapper_or::add(MipWrapperSolver *solver,
					  bool top_level){

  if(!has_been_added()){
    _solver = solver;
    DBG("Adding or constraint %s\n", "");
    _vars[0] = _vars[0]->add(solver, false);
    _vars[1] = _vars[1]->add(solver, false);
    if(top_level){
      LinearConstraint *con = new LinearConstraint(1, 2);
      con->add_coef(_vars[0], 1);
      con->add_coef(_vars[1], 1);
      solver->_constraints.push_back(con);
    } else {
      MipWrapper_Expression *C = new MipWrapper_IntVar(0, 1);
      C = C->add(solver, false);
      (new MipWrapper_ge(C, _vars[0]))->add(solver, true);
      (new MipWrapper_ge(C, _vars[1]))->add(solver, true);
      
      LinearConstraint *con = new LinearConstraint(0, INFINITY);
      con->add_coef(_vars[0], 1);
      con->add_coef(_vars[1], 1);
      con->add_coef(C, -1);
      solver->_constraints.push_back(con);
      
      _repr = C;
      return C;
    }
  }
  if(top_level) return this;
  return _repr;
}

MipWrapper_not::MipWrapper_not(MipWrapper_Expression *var1):
  MipWrapper_binop(var1, 0.0){ DBG("Creating not constraint %s\n", ""); }

MipWrapper_not::~MipWrapper_not(){}

MipWrapper_Expression* MipWrapper_not::add(MipWrapperSolver *solver,
					   bool top_level){
  if(!has_been_added()){
    _solver = solver;
    _vars[0] = _vars[0]->add(solver, false);
    if(top_level){
      _vars[0]->_upper = 0;
    } else {
      MipWrapper_Expression *c = new MipWrapper_IntVar(0, 1);
      c->add(solver, false);
      LinearConstraint *con = new LinearConstraint(1, 1);
      con->add_coef(_vars[0], 1);
      con->add_coef(c, 1);
      solver->_constraints.push_back(con);
      _repr = c;
      return c;
    }
  }
  if(top_level) return this;
  return _repr;
} 

/**
 * Optimisation stuff
 */

// Minimise Class
MipWrapper_Minimise::MipWrapper_Minimise(MipWrapper_Expression *arg1):
  _obj(arg1){}
MipWrapper_Minimise::~MipWrapper_Minimise(){}
MipWrapper_Expression* MipWrapper_Minimise::add(MipWrapperSolver *solver,
						bool top_level){
  _obj = _obj->add(solver, false);
  LinearConstraint *obj_con = new LinearConstraint(-INFINITY, INFINITY);
  obj_con->add_coef(_obj, 1);
  solver->_obj = obj_con;
  solver->_obj_coef = -1;
  return this;
}

// Maximise Class
MipWrapper_Maximise::MipWrapper_Maximise(MipWrapper_Expression *arg1):
    _obj(arg1){}
MipWrapper_Maximise::~MipWrapper_Maximise(){}

MipWrapper_Expression* MipWrapper_Maximise::add(MipWrapperSolver *solver,
						bool top_level){
  _obj = _obj->add(solver, false);
  LinearConstraint *obj_con = new LinearConstraint(-INFINITY, INFINITY);
  obj_con->add_coef(_obj, 1);
  solver->_obj = obj_con;
   solver->_obj_coef = 1;
  return this;
}

/**************************************************************
 ********************     Solver        ***********************
 **************************************************************/

MipWrapperSolver::MipWrapperSolver(){
  DBG("Create a MIP solver %s\n", "");
  var_counter = 0;
  _verbosity = 0;
  _obj = NULL;
}

MipWrapperSolver::~MipWrapperSolver(){ DBG("Delete a MIP Solver %s\n", "");}

void MipWrapperSolver::add(MipWrapper_Expression* arg){
   DBG("Adding an expression %s\n", "");
  if(arg != NULL) arg->add(this, true);
}

void MipWrapperSolver::add_expr(MipWrapper_Expression *expr){
  _exprs.push_back(expr);
}

void MipWrapperSolver::add_int_array(MipWrapperIntArray *arr){
  _int_array.push_back(arr);
}

void MipWrapperSolver::add_var_array(MipWrapperExpArray *arr){
  _var_array.push_back(arr);
}

void MipWrapperSolver::initialise(MipWrapperExpArray& arg){
  DBG("Initialise the Solver with search vars %s\n", "");

  initialise();
}

void MipWrapperSolver::initialise(){
  DBG("Initialise the Solver %s\n", "");
  
   // This is where I need to filter through the the linear constraints
  // and deal with anything that
  
  for(unsigned int i = 0; i < _constraints.size(); ++i){
    // For all constraints
    for(unsigned int j = 0; j < _constraints[i]->_variables.size(); ++j){
      if(_constraints[i]->_variables[j]->_expr_encoding != NULL){
	double coef = _constraints[i]->_coefficients[j];
	_constraints[i]->add_coef(_constraints[i]->_variables[j], coef);
	
	_constraints[i]->_variables.erase(
			  _constraints[i]->_variables.begin()+j);
	_constraints[i]->_coefficients.erase(
			  _constraints[i]->_coefficients.begin()+j);
      }
    }
  }
  
}

int MipWrapperSolver::solve(){
  DBG("Solve %s\n", "");
  
  initialise();
  
  if(_verbosity) std::cout << "c solve with mip wrapper " << std::endl;
  for(unsigned int i = 0; i < _constraints.size(); ++i)
    _constraints[i]->display();
  
  return SAT;
}

int MipWrapperSolver::solveAndRestart(const int policy,
				const unsigned int base,
				const double factor,
				const double decay ){
  return solve();
}

int MipWrapperSolver::startNewSearch(){
  if(_verbosity) std::cout << "c start a new interuptable search" << std::endl;
  return 0;
}

int MipWrapperSolver::getNextSolution(){
  if(_verbosity) std::cout << "c seek next solution" << std::endl;
  return 0;
}

int MipWrapperSolver::sacPreprocess(const int type){
  if(_verbosity) std::cout << "c enforces singleton arc consistency" << std::endl;
  return 0;
}

void MipWrapperSolver::setHeuristic(const char* var_heuristic,
				    const char* val_heuristic,
				    const int rand){
  if(_verbosity) std::cout << "c set the variable/value ordering (ignored)" << std::endl;
}

void MipWrapperSolver::setFailureLimit(const int cutoff){
  if(_verbosity) std::cout << "c set a cutoff on failures" << std::endl;
}

void MipWrapperSolver::setTimeLimit(const int cutoff){}

void MipWrapperSolver::setNodeLimit(const int cutoff){}

void MipWrapperSolver::setVerbosity(const int degree){
  _verbosity = degree;
}

void MipWrapperSolver::setRandomized(const int degree){
  if(_verbosity) std::cout << "c set the type of randomization" << std::endl;
}

void MipWrapperSolver::setRandomSeed(const int seed){
  if(_verbosity) std::cout << "c set the random seed" << std::endl;
}

bool MipWrapperSolver::is_sat(){
  return true;
}

int MipWrapperSolver::get_lin_cons_size(){
  return _constraints.size();
}

bool MipWrapperSolver::is_unsat(){
  return ! is_sat();
}

void MipWrapperSolver::printStatistics(){
  if(_verbosity) std::cout << "\td Time: " << getTime() << "\tNodes:" << getNodes()
	    << std::endl;
}

int MipWrapperSolver::getBacktracks(){
  if(_verbosity) std::cout << "c print the number of backtracks" << std::endl;
  return 0;
}

int MipWrapperSolver::getNodes(){
  return 0;
}

int MipWrapperSolver::getFailures(){
  if(_verbosity) std::cout << "c print the number of failures" << std::endl;
  return 0;
}

int MipWrapperSolver::getChecks(){
  if(_verbosity) std::cout << "c print the number of checks" << std::endl;
  return 0;
}

int MipWrapperSolver::getPropags(){
  if(_verbosity) std::cout << "c print the number of propags" << std::endl;
  return 0;
}

double MipWrapperSolver::getTime(){
  return 0;
}

double MipWrapperSolver::get_value(void *ptr){ return 0; }

/**
 * Creates an empty linear constraint object
 */
LinearConstraint::LinearConstraint(double lhs, double rhs):
  _lhs(lhs), _rhs(rhs){}

LinearConstraint::~LinearConstraint(){}
    
void LinearConstraint::add_coef(MipWrapper_Expression* expr,
				double coef,
				bool use_encoding){
  DBG("Beginning to add in linear arguments\n%s", "");
  
  LINEAR_ARG* larg = expr->for_linear();
  _lhs -= larg->offset * coef;
  _rhs -= larg->offset * coef;
  for(int i = 0; i < expr->for_linear_size(); ++i)
    add_coef(larg+i, coef, use_encoding);
}

void LinearConstraint::add_coef(LINEAR_ARG* arg_struct,
				double coef,
				bool use_encoding){
  DBG("Adding in expression to linear constraint %s\n", "");
  if( use_encoding && arg_struct->expr->_expr_encoding != NULL){
    DBG("Adding in expression encoding %s\n", "");
    for(int i = arg_struct->expr->_lower; i <= arg_struct->expr->_upper; ++i )
      if(arg_struct->expr->_expr_encoding[i] != NULL)
	add_coef(arg_struct->expr->_expr_encoding[i], arg_struct->coef*i*coef);
  } else {
    _variables.push_back(arg_struct->expr);
    _coefficients.push_back(arg_struct->coef*coef);
  }			  
}
    
/**
 * Prints the linear constraint
 */
void LinearConstraint::display(){
  std::cout << _lhs << " <= ";
  for(unsigned int i = 0; i < _variables.size(); ++i){
    std::cout << _coefficients[i] << "*X" << _variables[i]->_ident;
    if(i + 1 < _variables.size()) std::cout << " + ";
  }
  std::cout << " <= " << _rhs << std::endl;
}
