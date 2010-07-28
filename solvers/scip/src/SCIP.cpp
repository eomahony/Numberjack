
#include <iostream>

#include "SCIP.hpp"
#include "scip_exception.hpp"
#include <math.h>

/**************************************************************
 ********************     EXPRESSION        *******************
 **************************************************************/

void Expression::initialise(bool c)
{ 
  _ident = -1;
  _scip = NULL;
  _var = NULL;
  _encoding = NULL;
  _lower = 0;
  _upper = 0; 
  _coef = 0; 
  _continuous = c;
}

Expression::Expression()
{
  initialise(true);

#ifdef _DEBUGWRAP
  std::cout << "creating an empty expression" << std::endl;
#endif

}

int Expression::get_size(){
  return _upper - _lower +1;
}

int Expression::get_max(){
  return _upper;
}

int Expression::get_min(){
  return _lower;
}

void Expression::display(){
  std::cout << _coef << "*X" << _ident << " ";
}

// SCIP_IntVar::SCIP_IntVar(const int nval)
// {
//   initialise(false);
//   _upper = (double)(nval-1);
//   _has_holes_in_domain = false;

// #ifdef _DEBUGWRAP
//   std::cout << "creating integer variable [" << 0 << ".." << (nval-1) << "]" << std::endl;
// #endif

// }

SCIP_IntVar::SCIP_IntVar(const int lb, const int ub)
{
  initialise(false);
  _upper = (double)ub;
  _lower = (double)lb;
  _has_holes_in_domain = false;
  nbj_ident = -1;

#ifdef _DEBUGWRAP
  std::cout << "creating integer variable [" << lb << ".." << ub << "]" << std::endl;
#endif 

}

SCIP_IntVar::SCIP_IntVar(const int lb, const int ub, const int ident)
{
  initialise(false);
  _upper = (double)ub;
  _lower = (double)lb;
  _has_holes_in_domain = false;
  nbj_ident = ident;

#ifdef _DEBUGWRAP
  std::cout << "creating integer variable [" << lb << ".." << ub << "]" << std::endl;
#endif 

}

SCIP_IntVar::SCIP_IntVar(SCIPIntArray& values, const int ident){
  _values = values;
  _has_holes_in_domain = true;
  
  _upper = _values.get_item(0);
  _lower = _values.get_item(0);
  
  nbj_ident = ident;

  for(int i = 0; i < _values.size(); ++i){
    if(_values.get_item(i) > _upper) _upper = _values.get_item(i);
    if(_values.get_item(i) < _lower) _lower = _values.get_item(i);
  }
}

void SCIP_IntVar::add(SCIPSolver* solver, bool top_level){
  Expression::add(solver, top_level);
  if(_has_holes_in_domain)
    encode(solver);
}

void SCIP_IntVar::encode(SCIPSolver* solver){
  
  if(_has_holes_in_domain){
    
#ifdef _DEBUGWRAP
  std::cout << "Making an encoding for a variable with holes in the domain" << std::endl;
#endif

    _scip = solver->get_scip();	

    if(!_encoding) {
      
      int k, lb=(int)_lower, ub = (int)_upper;
      int m = (ub-lb) + 1;
      int actual_var_idx = 0;

      _encoding = new SCIP_VAR*[m+1];
      _encoding[m] = _var;
      _encoding -= lb;
      
      if(_var == NULL)
	std::cout << "Variable is NULL" << std::endl;

      char res[100];

      for(k= lb ; k <= ub; ++k) {
	if(_values.get_item(actual_var_idx) == k){
	  actual_var_idx++;
	  sprintf(res, "x_%d_%d", _ident, k);

	  SCIP_CALL_EXC(SCIPcreateVar(_scip, _encoding+k,
				  res,
				  0,
				  1,
				  0, 
				  SCIP_VARTYPE_BINARY,
				  TRUE, FALSE, NULL, NULL, NULL, NULL) );      
	  SCIP_CALL_EXC( SCIPaddVar(_scip, _encoding[k]) );
      
	  solver->add_scip_var(_encoding[k]);
	} else {
	  // It is not in the domain
	  sprintf(res, "x_blank_%d_%d", _ident, k);

	  SCIP_CALL_EXC(SCIPcreateVar(_scip, _encoding+k,
				  res,
				  0,
				  0,
				  0, 
				  SCIP_VARTYPE_BINARY,
				  TRUE, FALSE, NULL, NULL, NULL, NULL) );      
	  SCIP_CALL_EXC( SCIPaddVar(_scip, _encoding[k]) );
      
	  solver->add_scip_var(_encoding[k]);
	}

#ifdef _DEBUGWRAP
	std::cout << "add encoding var to scip" << std::endl;
#endif
      }

      double var_vals[m+1];
      for(k=0; k<m; ++k)
	var_vals[k] = 1.0;

      SCIP_CONS *domain;
      SCIP_CALL_EXC( SCIPcreateConsLinear(_scip, &domain, "domain",
					  m, // # vars
					  _encoding+lb, // variables
					  var_vals, // values
					  1, // LHS
					  1, // RHS
					  TRUE, TRUE, TRUE,
					  TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
      SCIP_CALL_EXC( SCIPaddCons(_scip, domain) );

#ifdef _DEBUGWRAP
      std::cout << "add domain constraint to scip" << std::endl;
#endif

      solver->add_scip_cons(domain);

      int index = 0;
      for(k=lb; k<=ub; ++k)
	var_vals[index++] = k;
      var_vals[m] = -1.0;

      SCIP_CONS *channel;
      SCIP_CALL_EXC( SCIPcreateConsLinear(_scip, &channel,"channel",
					m+1, // # vars
					_encoding+lb, // variables
					var_vals, // values
					0, // LHS
					0, // RHS
					TRUE, TRUE, TRUE,
					TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
      SCIP_CALL_EXC( SCIPaddCons(_scip, channel) );

#ifdef _DEBUGWRAP
      std::cout << "add channel constraint to scip" << std::endl;
#endif

      solver->add_scip_cons(channel);
    
#ifdef _DEBUGWRAP
    std::cout << "Variable with holes in it's domain encoded" << std::endl;
#endif
    
    }
  }else {
    Expression::encode(solver);
  }
}

int SCIP_IntVar::get_value() const
{
  SCIP_SOL * sol = SCIPgetBestSol(_scip);
  return (int)(round(SCIPgetSolVal(_scip, sol, _var)));
}

SCIP_FloatVar::SCIP_FloatVar()
{
  initialise(true);
  _upper = 0;
  _lower = 0;

  nbj_ident = -1;

#ifdef _DEBUGWRAP
  std::cout << "creating continuous variable" << std::endl;
#endif

}

SCIP_FloatVar::SCIP_FloatVar(const double lb, const double ub)
{
  initialise(true);
  _upper = ub;
  _lower = lb;

  nbj_ident = -1;

#ifdef _DEBUGWRAP
  std::cout << "creating continuous variable [" << lb << ".." << ub << "]" << std::endl;
#endif

}

SCIP_FloatVar::SCIP_FloatVar(const double lb, const double ub, const int ident)
{
  initialise(true);
  _upper = ub;
  _lower = lb;

  nbj_ident = ident;

#ifdef _DEBUGWRAP
  std::cout << "creating continuous variable [" << lb << ".." << ub << "]" << std::endl;
#endif

}

double SCIP_FloatVar::get_value() const
{
  SCIP_SOL * sol = SCIPgetBestSol(_scip);
  return SCIPgetSolVal(_scip, sol, _var);
}

Expression::~Expression()
{
  if(_encoding) {
    _encoding += (int)_lower;
    delete [] _encoding;
  }
}

void Expression::encode(SCIPSolver *solver)
{
  _scip = solver->get_scip();

  if(!_encoding) {
    int k, lb=(int)_lower, ub = (int)_upper;
    int m = (ub-lb+1);

    _encoding = new SCIP_VAR*[m+1];
    _encoding[m] = _var;
    _encoding -= lb;

#ifdef _DEBUGWRAP
    std::cout << "Creating an encoding with " << (m+1) << " vars" << std::endl;
#endif

    char res[100];

    for(k=lb; k<=ub; ++k) {

      sprintf(res, "x_%d_%d", _ident, k);

      SCIP_CALL_EXC(SCIPcreateVar(_scip, _encoding+k,
				  res,
				  0,
				  1,
				  0, 
				  SCIP_VARTYPE_BINARY,
				  TRUE, FALSE, NULL, NULL, NULL, NULL) );      
      SCIP_CALL_EXC( SCIPaddVar(_scip, _encoding[k]) );
      
      solver->add_scip_var(_encoding[k]);

#ifdef _DEBUGWRAP
      std::cout << "add encoding var to scip" << std::endl;
#endif

    }

    double var_vals[m+1];
    for(k=0; k<m; ++k)
      var_vals[k] = 1.0;

    SCIP_CONS *domain;
    SCIP_CALL_EXC( SCIPcreateConsLinear(_scip, &domain,"domain",
					m, // # vars
					_encoding+lb, // variables
					var_vals, // values
					1, // LHS
					1, // RHS
					TRUE, TRUE, TRUE,
					TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
    SCIP_CALL_EXC( SCIPaddCons(_scip, domain) );

#ifdef _DEBUGWRAP
    std::cout << "add domain constraint to scip" << std::endl;
#endif

    solver->add_scip_cons(domain);

    for(k=0; k<m; ++k)
      var_vals[k] = (double)(k+lb);
    var_vals[m] = -1.0;

    SCIP_CONS *channel;
    SCIP_CALL_EXC( SCIPcreateConsLinear(_scip, &channel,"channel",
					m+1, // # vars
					_encoding+lb, // variables
					var_vals, // values
					0, // LHS
					0, // RHS
					TRUE, TRUE, TRUE,
					TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
    SCIP_CALL_EXC( SCIPaddCons(_scip, channel) );

#ifdef _DEBUGWRAP
    std::cout << "add channel constraint to scip" << std::endl;
#endif

    solver->add_scip_cons(channel);
    
#ifdef _DEBUGWRAP
    std::cout << "Variable encoded" << std::endl;
#endif

  }
}

bool Expression::has_been_added() const
{
  return (_scip != NULL);
}

void Expression::add(SCIPSolver *solver, bool top_level){

#ifdef _DEBUGWRAP
  std::cout << "add variable [" << _lower << ".." << _upper << "]" << std::endl;
#endif

  if(!has_been_added()){
  
    _scip = solver->get_scip();

    _ident = solver->var_counter++;
  
    SCIP_VAR* var_ptr;
    
    SCIP_VARTYPE type;
    
    if(_continuous) type = SCIP_VARTYPE_CONTINUOUS;
    else type = SCIP_VARTYPE_INTEGER;
    
    char res[100];
    sprintf(res, "X_%d", _ident);
    
    SCIP_CALL_EXC(
		  SCIPcreateVar(_scip, &var_ptr,
				res,
				_lower, // LB
				_upper, // UB
				_coef, // ective
				type,
				TRUE, FALSE, NULL, NULL, NULL, NULL) );
    
    SCIP_CALL_EXC( SCIPaddVar(_scip, var_ptr) );
    _var = var_ptr;

#ifdef _DEBUGWRAP
    std::cout << "add variable to scip" << std::endl;    
#endif

    solver->add_scip_var(_var);
  }
}

SCIP_add::SCIP_add(Expression *arg1, Expression *arg2)
  : SCIP_Sum()
{
  addVar(arg1);
  addVar(arg2);
  addWeight(1);
  addWeight(1);
  initialise();
}

SCIP_add::SCIP_add(Expression *arg1, const int arg2)
  : SCIP_Sum()
{
  addVar(arg1);
  addWeight(1);
  set_rhs(-arg2);
  initialise();
}

SCIP_add::~SCIP_add() {}


SCIP_sub::SCIP_sub(Expression *arg1, const int arg2)
  : SCIP_Sum()
{
  addVar(arg1);
  addWeight(1);
  set_rhs(arg2);
  initialise();
}

SCIP_sub::SCIP_sub(Expression *arg1, Expression *arg2)
  : SCIP_Sum()
{
  addVar(arg1);
  addVar(arg2);
  addWeight(1);
  addWeight(-1);
  initialise();
}

SCIP_sub::~SCIP_sub() {}


SCIP_AllDiff::SCIP_AllDiff( Expression* arg1, Expression* arg2 ) 
  : SCIP_Flow()
{
  addVar(arg1);
  addVar(arg2);
  initbounds();

#ifdef _DEBUGWRAP
  std::cout << "creating alldiff" << std::endl;
#endif

  std::fill(card_ub+min_val, card_ub+max_val+1, 1);

}

SCIP_AllDiff::SCIP_AllDiff( SCIPExpArray& vars ) 
  : SCIP_Flow(vars) 
{

#ifdef _DEBUGWRAP
  std::cout << "creating alldiff" << std::endl;
#endif

  std::fill(card_ub+min_val, card_ub+max_val+1, 1);

}

SCIP_AllDiff::~SCIP_AllDiff()
{
}

SCIP_Gcc::SCIP_Gcc( SCIPExpArray& vars, SCIPIntArray& vals, SCIPIntArray& lb, SCIPIntArray& ub ) 
  : SCIP_Flow(vars) 
{

#ifdef _DEBUGWRAP
  std::cout << "creating alldiff" << std::endl;
#endif

  int i, j, n=vals.size();
  for(i=0; i<n; ++i)
    {
      j = vals.get_item(i);
      card_lb[j] = lb.get_item(i);
      card_ub[j] = ub.get_item(i);
    }

}

SCIP_Gcc::~SCIP_Gcc()
{
}


SCIP_Flow::SCIP_Flow(SCIPExpArray& vars)
  : SCIP_FloatVar() 
{
  _vars = vars;
  initbounds();
}

SCIP_Flow::SCIP_Flow()
  : SCIP_FloatVar() 
{
}


void SCIP_Flow::initbounds() {
  int i, lb, ub;
  max_val = -INT_MAX;
  min_val = INT_MAX;
  
  for(i = 0; i < _vars.size(); ++i) {
    
    lb = (int)(_vars.get_item(i)->_lower);
    ub = (int)(_vars.get_item(i)->_upper);
    
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

SCIP_Flow::~SCIP_Flow(){
}

void SCIP_Flow::addVar(Expression* v){
  _vars.add(v);
}

void SCIP_Flow::add(SCIPSolver *solver, bool top_level){
  if(!has_been_added()){

    _scip = solver->get_scip();

#ifdef _DEBUGWRAP
    std::cout << "add flow" << std::endl;
#endif

    if(top_level){

      int i, j, k;
//       max_val = -INT_MAX;
//       min_val = INT_MAX;
  
      for(i = 0; i < _vars.size(); ++i) {
	_vars.get_item(i)->add(solver, false);
	_vars.get_item(i)->encode(solver);

// 	lb = (int)(_vars.get_item(i)->_lower);
// 	ub = (int)(_vars.get_item(i)->_upper);

// 	if(lb < min_val) min_val = lb;
// 	if(ub > max_val) max_val = ub;
      }

      SCIP_CONS *cons;

      SCIP_VAR *var_array[_vars.size()];
      double var_vals[_vars.size()];
      std::fill(var_vals, var_vals+_vars.size(), 1.0);

      for(i = min_val; i <= max_val; ++i) {
	k=0;
	for(j = 0; j < _vars.size(); ++j)
	  if( i <= (int)(_vars.get_item(j)->_upper) && i >= (int)(_vars.get_item(j)->_lower) )
	    var_array[k++] = _vars.get_item(j)->_encoding[i];
      
	SCIP_CALL_EXC( SCIPcreateConsLinear(_scip, & cons,"Flow value constraint",
					    k, // # vars
					    var_array, // variables
					    var_vals, // values
					    card_lb[i], // LHS
					    card_ub[i], // RHS
					    TRUE, TRUE, TRUE,
					    TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
	SCIP_CALL_EXC( SCIPaddCons(_scip, cons) );
      
#ifdef _DEBUGWRAP
	std::cout << "add flow constraint to scip" << std::endl;    
#endif    

	solver->add_scip_cons(cons);
      }
  
    } else {
      std::cout << "Warning Flow constraint supported only on top level" << std::endl;
    }
  }
}


SCIP_Sum::SCIP_Sum(SCIPExpArray& vars, 
		   SCIPIntArray& weights, 
		   const int offset)
  : SCIP_FloatVar() 
{
  _offset = offset;
  _vars = vars;
  _weights = weights;
  initialise();
}

SCIP_Sum::SCIP_Sum(Expression *arg1, 
		   Expression *arg2, 
		   SCIPIntArray& w, 
		   const int offset)
  : SCIP_FloatVar() 
{
  _offset = offset;
  _vars.add(arg1);
  _vars.add(arg2);
  _weights = w;
  initialise();
}

SCIP_Sum::SCIP_Sum(Expression *arg, 
		   SCIPIntArray& w, 
		   const int offset)
  : SCIP_FloatVar() 
{
  _offset = offset;
  _vars.add(arg);
  _weights = w;
  initialise();
}

SCIP_Sum::SCIP_Sum()
  : SCIP_FloatVar()
{
  _offset = 0;
}

void SCIP_Sum::initialise() {

#ifdef _DEBUGWRAP
  std::cout << "creating sum: Size of parameters is " << _vars.size() 
	    << std::endl; 
#endif
  
  //_lower = 0;
  //_upper = 0;

  int weight;
  for(int i = 0; i < _vars.size(); ++i){
    weight = _weights.get_item(i);
    
#ifdef _DEBUGWRAP
    std::cout << "bounds of sum variable:" << _vars.get_item(i)->_lower 
	      << " " << _vars.get_item(i)->_upper << std::endl;
#endif    
    
	//assert(weight);
    
    if( weight > 0 ) {
      _lower += (weight * _vars.get_item(i)->_lower);
      _upper += (weight * _vars.get_item(i)->_upper);
    } else {
      _upper += (weight * _vars.get_item(i)->_lower);
      _lower += (weight * _vars.get_item(i)->_upper);
    }
  }
  
  //std::cout << _offset << std::endl;
  
  _lower += _offset;
  _upper += _offset;
      
  // #ifdef _DEBUGWRAP
//       std::cout << "bounds right before getting added: " <<
// 	_lower << ":" << _upper << std::endl;
// #endif

#ifdef _DEBUGWRAP
  std::cout << "Intermediate variable has values: " << _lower << " to " 
	    << _upper << std::endl;
#endif

}

SCIP_Sum::~SCIP_Sum(){

#ifdef _DEBUGWRAP
  std::cout << "delete sum" << std::endl;
#endif

}

void SCIP_Sum::addVar(Expression* v){
  _vars.add(v);
}

void SCIP_Sum::addWeight(const int w){
  _weights.add(w);
}

void SCIP_Sum::set_rhs(const int k){
  _offset = k;
}

void SCIP_Sum::add(SCIPSolver *solver, bool top_level){
  if(!has_been_added()){

#ifdef _DEBUGWRAP    
    std::cout << "add sum" << std::endl;
#endif

    _scip = solver->get_scip();

    if(!top_level){
      
      for(int i = 0; i < _vars.size(); ++i)
	_vars.get_item(i)->add(solver, false);
      _scip = NULL;
      
//       int weight;
//       for(int i = 0; i < _vars.size(); ++i){
// 	weight = _weights.get_item(i);
	
// #ifdef _DEBUGWRAP
// 	std::cout << "bounds of sum variable:" << _vars.get_item(i)->_lower 
// 		  << " " << _vars.get_item(i)->_upper << std::endl;
// #endif    
	
// 	//assert(weight);
	
// 	if( weight > 0 ) {
// 	  _lower += (weight * _vars.get_item(i)->_lower);
// 	  _upper += (weight * _vars.get_item(i)->_upper);
// 	} else {
// 	  _upper += (weight * _vars.get_item(i)->_lower);
// 	  _lower += (weight * _vars.get_item(i)->_upper);
// 	}
//       }
      
//       //std::cout << _offset << std::endl;
      
//       _lower -= _offset;
//       _upper -= _offset;
      
// #ifdef _DEBUGWRAP
//       std::cout << "bounds right before getting added: " <<
// 	_lower << ":" << _upper << std::endl;
// #endif
      
     
      Expression::add(solver, false);
      
      
      SCIP_VAR *var_array[_vars.size() + 1];
      double var_vals[_vars.size() + 1];
      
      for(int i = 0; i < _vars.size(); ++i){
	var_array[i] = _vars.get_item(i)->_var;
	var_vals[i] = (double)(_weights.get_item(i));
      }
      
      var_array[_vars.size()] = _var;
      var_vals[_vars.size()] = -1;
      
      SCIP_CONS *cons;
      
      
      SCIP_CALL_EXC( SCIPcreateConsLinear(_scip, & cons,"Sum eq intervar con",
					  _vars.size() + 1, // # vars
					  var_array, // variables
					  var_vals, // values
					  -_offset, // LHS
					  -_offset, // RHS
					  TRUE, TRUE, TRUE,
					  TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
      SCIP_CALL_EXC( SCIPaddCons(_scip, cons) );
      
      solver->add_scip_cons(cons);

#ifdef _DEBUGWRAP      
      std::cout << "add sum constraint to scip" << std::endl;
#endif      
      
    } else {
      std::cout << "Warning SUM constraint on top level not supporeted" << std::endl;
    }
  }
}

/* Binary operators */

SCIP_binop::SCIP_binop(Expression *var1, Expression *var2)
  : SCIP_FloatVar()
{
  _vars[0] = var1;
  _vars[1] = var2;
  _is_proper_coef = false;
  _upper = 1.0;
  _lower = 0.0;


#ifdef _DEBUGWRAP
  std::cout << "creating binary operator" << std::endl;
#endif

}

SCIP_binop::SCIP_binop(Expression *var1, double rhs)
  : SCIP_FloatVar()
{
  _vars[0] = var1;
  _vars[1] = NULL;
  _rhs = rhs;
  _is_proper_coef = true;
  _upper = 1.0;
  _lower = 0.0;


#ifdef _DEBUGWRAP
  std::cout << "creating binary operator" << std::endl;
#endif

}


SCIP_binop::~SCIP_binop(){

#ifdef _DEBUGWRAP
  std::cout << "delete binary operator" << std::endl;
#endif

}

SCIP_NoOverlap::SCIP_NoOverlap(Expression *var1, Expression *var2, SCIPIntArray &coefs):
  SCIP_binop(var1, var2), _coefs(coefs){
  // Nothing else to do
}

SCIP_NoOverlap::~SCIP_NoOverlap(){
  // Nothing to do here I think
}

void SCIP_NoOverlap::add(SCIPSolver *solver, bool top_level){
  
  _vars[0]->add(solver, false);
  _vars[1]->add(solver, false);
  
  if(top_level){
    
    SCIPIntArray *arr1 = new SCIPIntArray();
    arr1->add(_coefs.get_item(0));
    solver->add_scip_int_array(arr1);
    
    SCIPIntArray *arr2 = new SCIPIntArray();
    arr2->add(_coefs.get_item(1));
    solver->add_scip_int_array(arr2);
    
    Expression *prec1 = new SCIP_Precedence(_vars[0], _vars[1], *arr1);
    solver->add_scip_expr(prec1);
    
    Expression *prec2 = new SCIP_Precedence(_vars[1], _vars[0], *arr2);
    solver->add_scip_expr(prec2);
    
    Expression *orexp = new SCIP_or(prec1, prec2);
    solver->add_scip_expr(orexp);
    orexp->add(solver, true);
    
    /*
    
    // Sort out b1 + b2 = 1
    
    Expression *b1 = new SCIP_IntVar(0, 1);
    solver->add_scip_expr(b1);
    
    Expression *b2 = new SCIP_IntVar(0, 1);
    solver->add_scip_expr(b2);
    
    SCIPExpArray *bs = new SCIPExpArray();
    bs->add(b1);
    bs->add(b2);
    solver->add_scip_var_array(bs);
    
    SCIPIntArray *bcoefs = new SCIPIntArray();
    bcoefs->add(1);
    bcoefs->add(1);
    solver->add_scip_int_array(bcoefs);
    
    Expression *bequn = new SCIP_Sum(*bs, *bcoefs);
    solver->add_scip_expr(bequn);
    
    Expression *eq0 = new SCIP_eq(bequn, new SCIP_IntVar(1,1));
    solver->add_scip_expr(eq0);
    
    eq0->add(solver, true);
    
    //std::cout << "Done b+b" << std::endl;
    
    //TODO: change me
    
    SCIP *scip = solver->get_scip();
    
    // X = x1 + y2
    Expression *x1 = new SCIP_FloatVar(-SCIPinfinity(scip), SCIPinfinity(scip));
    solver->add_scip_expr(x1);
    
    Expression *y2 = new SCIP_FloatVar(-SCIPinfinity(scip), SCIPinfinity(scip));
    solver->add_scip_expr(y2);
    
    Expression *addx1y2 = new SCIP_add(x1, y2);
    solver->add_scip_expr(addx1y2);
    
    Expression *equalX = new SCIP_eq(_vars[0], addx1y2);
    solver->add_scip_expr(equalX);
    
    equalX->add(solver, true);
    
    //std::cout << "X = x1 + y1" << std::endl;
    
    // Y = x2 + y1
    Expression *x2 = new SCIP_FloatVar(-SCIPinfinity(scip), SCIPinfinity(scip));
    solver->add_scip_expr(x2);
    
    Expression *y1 = new SCIP_FloatVar(-SCIPinfinity(scip), SCIPinfinity(scip));
    solver->add_scip_expr(y1);
    
    Expression *addx2y1 = new SCIP_add(x2, y1);
    solver->add_scip_expr(addx2y1);
    
    Expression *equalY = new SCIP_eq(_vars[1], addx2y1);
    solver->add_scip_expr(equalY);
    
    equalY->add(solver, true);
    
    //std::cout << "Y = x2 + y2" << std::endl;
    
    // x1 - y1 >= By x b1
    // x1 - y1 - Byb1 >= 0
    
    SCIPExpArray *vars1 = new SCIPExpArray();
    vars1->add(x1);
    vars1->add(y1);
    vars1->add(b1);
    solver->add_scip_var_array(vars1);
    
    SCIPIntArray *coefs1 = new SCIPIntArray();
    coefs1->add(1);
    coefs1->add(-1);
    coefs1->add(-_coefs.get_item(1));
    solver->add_scip_int_array(coefs1);
    
    Expression *sum1 = new SCIP_Sum(*vars1, *coefs1);
    solver->add_scip_expr(sum1);
    
    Expression *equn1 = new SCIP_ge(sum1, 0.0);
    solver->add_scip_expr(equn1);
    
    equn1->add(solver, true);
    
    //std::cout << "x1 - x1 >= By x b1" << std::endl;
    
    // y2 - x2 >= Bx x b2
    
    SCIPExpArray *vars2 = new SCIPExpArray();
    vars2->add(y2);
    vars2->add(x2);
    vars2->add(b2);
    solver->add_scip_var_array(vars2);
    
    SCIPIntArray *coefs2 = new SCIPIntArray();
    coefs2->add(1);
    coefs2->add(-1);
    coefs2->add(-_coefs.get_item(0));
    solver->add_scip_int_array(coefs2);
    
    Expression *sum2 = new SCIP_Sum(*vars2, *coefs2);
    solver->add_scip_expr(sum2);
    
    Expression *equn2 = new SCIP_ge(sum2, 0.0);
    solver->add_scip_expr(equn2);
    
    equn2->add(solver, true);
    
    //std::cout << "y2 - x2 >= Bx x b2" << std::endl;
    
    */
    
  } else {
    
    std::cout << "No support for NoOverLap operator not in top level" << std::endl;
    
  }
  
}

SCIP_Precedence::SCIP_Precedence(Expression *var1, Expression *var2, SCIPIntArray &coefs):
  SCIP_binop(var1, var2), _coefs(coefs){
  
}

SCIP_Precedence::~SCIP_Precedence(){
  
}

void SCIP_Precedence::add(SCIPSolver *solver, bool top_level){
  _vars[0]->add(solver, false);
  _vars[1]->add(solver, false);
  
  if(top_level){
    
    Expression *addlhs = new SCIP_add(_vars[0], _coefs.get_item(0));
    solver->add_scip_expr(addlhs);
    
    Expression *equn = new SCIP_ge(_vars[1], addlhs);
    solver->add_scip_expr(equn);
    
    equn->add(solver, true);
    
  } else {
    std::cout << "Precedence at top level not supported yet" << std::endl;
    
    Expression *addlhs = new SCIP_add(_vars[0], _coefs.get_item(0));
    solver->add_scip_expr(addlhs);
    
    Expression *equn = new SCIP_ge(_vars[1], addlhs);
    solver->add_scip_expr(equn);
    
    equn->add(solver, false);
    
    // Link the variables
    _var = equn->_var;
    
  }
}

SCIP_eq::SCIP_eq(Expression *var1, Expression *var2)
  : SCIP_binop(var1,var2)
{

#ifdef _DEBUGWRAP
  std::cout << "creating equality" << std::endl;
#endif

}

SCIP_eq::SCIP_eq(Expression *var1, double rhs)
  : SCIP_binop(var1,rhs)
{

#ifdef _DEBUGWRAP
  std::cout << "creating equality" << std::endl;
#endif

}

SCIP_eq::~SCIP_eq(){

#ifdef _DEBUGWRAP
  std::cout << "delete eq" << std::endl;
#endif

}

void SCIP_eq::add(SCIPSolver *solver, bool top_level){
  if(!has_been_added()){

    _scip = solver->get_scip();

#ifdef _DEBUGWRAP
    std::cout << "add equality" << std::endl;
#endif
    
    // Add children
    _vars[0]->add(solver, false);
    
    if(top_level){
      
      if(_is_proper_coef){
      
	SCIP_VAR *vars[1];
	
	vars[0] = _vars[0]->_var;
	
	SCIP_CONS *cons;
	
	double var_vals[1];
	var_vals[0] = 1.0;
	
	SCIP_CALL_EXC( SCIPcreateConsLinear(_scip, &cons, "Equality constraint",
					    2, // # vars
					    vars, // variables
					    var_vals, // values
					    _rhs, // LHS
					    _rhs, // RHS
					    TRUE, TRUE, TRUE,
					    TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
	SCIP_CALL_EXC( SCIPaddCons(_scip, cons) );
	
	solver->add_scip_cons(cons);

#ifdef _DEBUGWRAP	
	std::cout << "add equality constraint to scip" << std::endl;    
#endif
	
      } else {
	_vars[1]->add(solver, false);
	
	SCIP_VAR *vars[2];
      
	vars[0] = _vars[0]->_var;
	vars[1] = _vars[1]->_var;
	
	SCIP_CONS *cons;
	
	double var_vals[2];
	var_vals[0] = 1.0;
	var_vals[1] = -1.0;
	
	SCIP_CALL_EXC( SCIPcreateConsLinear(_scip, &cons, "Equality constraint",
					    2, // # vars
					    vars, // variables
					    var_vals, // values
					    0, // LHS
					    0, // RHS
					    TRUE, TRUE, TRUE,
					    TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
	SCIP_CALL_EXC( SCIPaddCons(_scip, cons) );

	solver->add_scip_cons(cons);

#ifdef _DEBUGWRAP	
	std::cout << "add equality constraint to scip" << std::endl;    
#endif

      }
    
#ifdef _DEBUGWRAP
      std::cout << "create equal predicate" << std::endl;
#endif

    } else {
     
      // First create a variable y = |a-b|
      // Get the bounds
      // Create the variable
      _vars[1]->add(solver, false);

      Expression *c = new SCIP_IntVar(0, 1);
      solver->add_scip_expr(c);
      
      SCIPIntArray *coef_array = new SCIPIntArray();
      coef_array->add(1);
      coef_array->add(1);
      solver->add_scip_int_array(coef_array);

      Expression *neq_reif = new SCIP_ne(_vars[0], _vars[1]);
      neq_reif->add(solver, false);
      solver->add_scip_expr(neq_reif);

      SCIPExpArray *var_array = new SCIPExpArray();
      var_array->add(c);
      var_array->add(neq_reif);
      solver->add_scip_var_array(var_array);

      Expression *eqcon = new SCIP_eq(new SCIP_Sum(*var_array, 
							*coef_array), 1);
      eqcon->add(solver, true); // Cause it's gotta be true
      solver->add_scip_expr(eqcon);

      _var = c->_var;
    }
  }
}

/* Disequality operator */

SCIP_ne::SCIP_ne(Expression *var1, Expression *var2)
  : SCIP_binop(var1,var2)
{
}

SCIP_ne::SCIP_ne(Expression *var1, double rhs)
  : SCIP_binop(var1,rhs)
{
}

SCIP_ne::~SCIP_ne(){

#ifdef _DEBUGWRAP
  std::cout << "delete notequal" << std::endl;
#endif

}

void SCIP_ne::add(SCIPSolver *solver, bool top_level){
  if(!has_been_added()){
    _scip = solver->get_scip();
    
    // Add children
    _vars[0]->add(solver, false);
    
    if(top_level){
      
      //TODO: Add in support for constants
      if(_is_proper_coef){
	
	_vars[0]->encode(solver);
	
	double var_vals[1];
        var_vals[0] = 1.0;
	
	SCIP_VAR *vars[1];
	vars[0] = _vars[0]->_encoding[(int)_rhs];
	
	SCIP_CONS *cons;
	
	SCIP_CALL_EXC( SCIPcreateConsLinear(_scip, &cons, "Not equal constraint",
					    1, // # vars
					    vars, // variables
					    var_vals, // values
					    0, // LHS
					    0, // RHS
					    TRUE, TRUE, TRUE,
					    TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
	SCIP_CALL_EXC( SCIPaddCons(_scip, cons) );
	solver->add_scip_cons(cons);
	
      } else {
      
        _vars[1]->add(solver, false);
        
        _vars[0]->encode(solver);
        _vars[1]->encode(solver);  
        
        SCIP_VAR *vars[2];
        
        double var_vals[2];
        var_vals[0] = 1.0;
        var_vals[1] = 1.0;
        
        int ub = std::min(_vars[0]->_upper, _vars[1]->_upper);
        int lb = std::max(_vars[0]->_lower, _vars[1]->_lower);
        
        for(int i=lb; i<=ub; ++i) {
  	
	  vars[0] = _vars[0]->_encoding[i];
	  vars[1] = _vars[1]->_encoding[i];
  	
#ifdef _DEBUGWRAP
	  if(vars[0] == NULL)
	    std::cout << "First one is NULL" << std::endl;
	  if(vars[1] == NULL)
	    std::cout << "Second one is NULL" << std::endl;
#endif
	
	  SCIP_CONS *cons;
	
	  SCIP_CALL_EXC( SCIPcreateConsLinear(_scip, &cons, "Not equal constraint",
					    2, // # vars
					    vars, // variables
					    var_vals, // values
					    0, // LHS
					    1, // RHS
					    TRUE, TRUE, TRUE,
					    TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
	  SCIP_CALL_EXC( SCIPaddCons(_scip, cons) );
	
	  solver->add_scip_cons(cons);
      }
	
    }
      
#ifdef _DEBUGWRAP
      std::cout << "create notequal predicate" << std::endl;
#endif
      
    } else {
      
      // This is where we need to add the reification
      _vars[1]->add(solver, false);
      
      Expression *c = new SCIP_IntVar(0, 1);
      c->add(solver, false);
      
      _vars[0]->encode(solver);
      _vars[1]->encode(solver);
      
      SCIP_VAR **var0encoding = _vars[0]->_encoding;
      SCIP_VAR **var1encoding = _vars[1]->_encoding;
      
      double lower = _vars[0]->_lower > _vars[1]->_lower ? 
					  _vars[0]->_lower : _vars[1]->_lower;
      double upper = _vars[0]->_upper < _vars[1]->_upper ?
	_vars[0]->_upper : _vars[1]->_upper;
    
      _scip = solver->get_scip();
      
      for(int i = lower; i <= upper; ++i){
	
	// For every pair of variables
	// -inf <= Xi + Yi <= 2 - C
	// -inf <= Xi + Yi + C <= 2
	
	SCIP_VAR *temp_vars[3];
	
	temp_vars[0] = NULL;
	
	temp_vars[0] = var0encoding[(int)(i)];
	
	temp_vars[1] = var1encoding[(int)(i)];
	
	temp_vars[2] = c->_var;
	
	double var_vals[3];
	var_vals[0] = 1.0;
	var_vals[1] = 1.0;
	var_vals[2] = 1.0;

	SCIP_CONS *cons;
	SCIP_CALL_EXC( SCIPcreateConsLinear(_scip, &cons, 
					    "Not equal reified constraint",
					    3, // # vars
					    temp_vars, // variables
					    var_vals, // values
					    -SCIPinfinity(_scip), // LHS
					    2, // RHS
					    TRUE, TRUE, TRUE,
					    TRUE, TRUE, FALSE, FALSE, FALSE, 
					    FALSE, FALSE) );
	SCIP_CALL_EXC( SCIPaddCons(_scip, cons) );
	solver->add_scip_cons(cons);
      }

      // This is what to do when c is 0

      double difference_a = _vars[0]->_upper - _vars[1]->_lower;
      double difference_b = _vars[1]->_upper - _vars[0]->_lower;
      double difference = difference_a;
      if(difference_a < difference_b)
	difference = difference_b;
	
      // Bounds for y = {0, difference}
      // y is continuous
      
      SCIP_FloatVar *y = new SCIP_FloatVar(0, difference);
      solver->add_scip_expr(y);
      
      // a-b <= y
      SCIPIntArray *intarr1 = new SCIPIntArray();
      intarr1->add(1);
      intarr1->add(-1);
      solver->add_scip_int_array(intarr1);

      SCIP_le *equation1 = new SCIP_le(new SCIP_Sum(_vars[0], _vars[1], *intarr1), y);
      solver->add_scip_expr(equation1);

      // b-a <= y
      
      SCIPIntArray *intarr2 = new SCIPIntArray();
      intarr2->add(-1);
      intarr2->add(1);
      solver->add_scip_int_array(intarr2);

      SCIP_le *equation2 = new SCIP_le(new SCIP_Sum(_vars[0], _vars[1], *intarr2), y);
      solver->add_scip_expr(equation2);

      // add both on top level
      equation1->add(solver, true);
      equation2->add(solver, true);
      
      // Create constraint c.M <= y

      // M = upper bound of y + 1
      double M = difference+1;

      SCIPIntArray *intarr3 = new SCIPIntArray();
      intarr3->add(M);
      solver->add_scip_int_array(intarr3);

      SCIPExpArray *exparr1 = new SCIPExpArray();
      exparr1->add(c);
      solver->add_scip_var_array(exparr1);

      SCIP_le *equation3 = new SCIP_le(y, new SCIP_Sum(*exparr1, *intarr3));
      
      // add y <= c.M on top level
      equation3->add(solver, true);
      solver->add_scip_expr(equation3);
      
      // Use c as the variable for this expression
      _upper = 1;
      _lower = 0;
      _var = c->_var;

  }
}
}

/* Leq operator */

SCIP_le::SCIP_le(Expression *var1, Expression *var2)
  : SCIP_binop(var1,var2)
{ 

#ifdef _DEBUGWRAP
  std::cout << "Second one is ok" << std::endl;
#endif

}

SCIP_le::SCIP_le(Expression *var1, double rhs)
  : SCIP_binop(var1,rhs)
{
}


SCIP_le::~SCIP_le(){

#ifdef _DEBUGWRAP
  std::cout << "delete lessequal" << std::endl;
#endif

}

void SCIP_le::add(SCIPSolver *solver, bool top_level){
  if(!has_been_added()){

    _scip = solver->get_scip();

#ifdef _DEBUGWRAP
    std::cout << "Adding le " << std::endl;
#endif
  
    // Add children
    _vars[0]->add(solver, false);
    if(_vars[1]) _vars[1]->add(solver, false);
  
    if(top_level){
    
      if(_is_proper_coef){
      
	SCIP_VAR *vars[1];
      
	vars[0] = _vars[0]->_var;
      
#ifdef _DEBUGWRAP
	if(vars[0] == NULL)
	  std::cout << "First one is NULL" << std::endl;
	if(vars[1] == NULL)
	  std::cout << "Second one is NULL" << std::endl;
#endif
      
	SCIP_CONS *cons;
      
	double var_vals[1];
	var_vals[0] = 1.0;
    
	SCIP_CALL_EXC( SCIPcreateConsLinear(_scip, &cons, "Less or equal constraint",
					    1, // # vars
					    vars, // variables
					    var_vals, // values
					    -SCIPinfinity(_scip), // LHS
					    _rhs, // RHS
					    TRUE, TRUE, TRUE,
					    TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
	SCIP_CALL_EXC( SCIPaddCons(_scip, cons) );
      
	solver->add_scip_cons(cons);

      } else {
  
	_vars[1]->add(solver, false);
  
	SCIP_VAR *vars[2];
      
	vars[0] = _vars[0]->_var;
	vars[1] = _vars[1]->_var;

#ifdef _DEBUGWRAP      
	if(vars[0] == NULL)
	  std::cout << "First one is NULL" << std::endl;
	if(vars[1] == NULL)
	  std::cout << "Second one is NULL" << std::endl;
#endif
      
	SCIP_CONS *cons;
      
	double var_vals[2];
	var_vals[0] = 1.0;
	var_vals[1] = -1.0;
    
	SCIP_CALL_EXC( SCIPcreateConsLinear(_scip, &cons, "Less or equal constraint",
					    2, // # vars
					    vars, // variables
					    var_vals, // values
					    -SCIPinfinity(_scip), // LHS
					    0, // RHS
					    TRUE, TRUE, TRUE,
					    TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
	SCIP_CALL_EXC( SCIPaddCons(_scip, cons) );

	solver->add_scip_cons(cons);
      }

#ifdef _DEBUGWRAP    
      std::cout << "create Leq predicate" << std::endl;
#endif

    } else {
      std::cout << "Leq operator not in top level insupported at the moment" << std::endl;
      
    }
  }
}

/* Geq operator */

SCIP_ge::SCIP_ge(Expression *var1, Expression *var2)
  : SCIP_binop(var1,var2)
{
}

SCIP_ge::SCIP_ge(Expression *var1, double rhs)
  : SCIP_binop(var1,rhs)
{

#ifdef _DEBUGWRAP
  std::cout << "Geq on constant constructor" << std::endl;
#endif

}

SCIP_ge::~SCIP_ge(){

#ifdef _DEBUGWRAP
  std::cout << "delete greaterequal" << std::endl;
#endif

}

void SCIP_ge::add(SCIPSolver *solver, bool top_level){ 
  if(!has_been_added()){
    
    _scip = solver->get_scip(); 

    // Add children
    _vars[0]->add(solver, false);
  
  
    if(top_level){
  
      //std::cout << "Top level" << std::endl;
  
      if(_is_proper_coef){

	//std::cout << "Proper coef" << std::endl;

#ifdef _DEBUGWRAP      
	std::cout << "Geq on constant" << std::endl;
#endif
      
	SCIP_VAR *vars[1];
      
	vars[0] = _vars[0]->_var;
      
	SCIP_CONS *cons;
      
	double var_vals[1];
	var_vals[0] = 1.0;
    
	SCIP_CALL_EXC( SCIPcreateConsLinear(_scip, &cons, "Greater or equal constraint",
					    1, // # vars
					    vars, // variables
					    var_vals, // values
					    _rhs, // LHS
					    SCIPinfinity(_scip), // RHS
					    TRUE, TRUE, TRUE,
					    TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
      
	SCIP_CALL_EXC( SCIPaddCons(_scip, cons) );

	solver->add_scip_cons(cons);
      
      } else {
	
	//std::cout << "Not proper coef" << std::endl;
	
	_vars[1]->add(solver, false);
  
	SCIP_VAR *vars[2];
      
	vars[0] = _vars[0]->_var;
	vars[1] = _vars[1]->_var;

#ifdef _DEBUGWRAP      
	if(vars[0] == NULL)
	  std::cout << "First one is NULL" << std::endl;
	if(vars[1] == NULL)
	  std::cout << "Second one is NULL" << std::endl;
#endif
      
	SCIP_CONS *cons;
      
	double var_vals[2];
	var_vals[0] = 1.0;
	var_vals[1] = -1.0;
    
	SCIP_CALL_EXC( SCIPcreateConsLinear(_scip, &cons, "Greater or equal constraint",
					    2, // # vars
					    vars, // variables
					    var_vals, // values
					    0, // LHS
					    SCIPinfinity(_scip), // RHS
					    TRUE, TRUE, TRUE,
					    TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
	SCIP_CALL_EXC( SCIPaddCons(_scip, cons) );

	solver->add_scip_cons(cons);

      }
    
#ifdef _DEBUGWRAP
      std::cout << "create Geq predicate" << std::endl;
#endif

    } else {
      
      
      //std::cout << "Not top level" << std::endl;
      
      // Dirty dirty hack!
      if(_is_proper_coef){
	_vars[1] = new SCIP_FloatVar(_rhs, _rhs);
	_vars[1]->add(solver, false);
      } else {
	_vars[1]->add(solver, false);
      }
      
      //std::cout << "Added in cvar" << std::endl;
      
      // Create variable to be reified
      Expression *b = new SCIP_IntVar(0, 1);
      b->add(solver, false);
      solver->add_scip_expr(b);
      
      // calculate M, must be bigger than
      double M = _vars[0]->_upper - _vars[1]->_lower;
      if (_vars[1]->_upper - _vars[0]->_lower > M)
	M = _vars[1]->_upper - _vars[0]->_lower;
	
      //std::cout << "Starting equation 1" << std::endl;
        
      //std::cout << "M is " << M << std::endl;    
	
      // equation 1
      SCIPIntArray *intarr1 = new SCIPIntArray();
      intarr1->add(1.0);
      intarr1->add(-1.0);
      intarr1->add(-M);
      solver->add_scip_int_array(intarr1);

      SCIPExpArray *exparr1 = new SCIPExpArray();
      exparr1->add(_vars[0]);
      exparr1->add(_vars[1]);
      exparr1->add(b);
      solver->add_scip_var_array(exparr1);
      
      Expression *sum1 = new SCIP_Sum(*exparr1, *intarr1);
      solver->add_scip_expr(sum1);
      
      // Put M as a variable here as SCIP did not like the neg coef on the RHS
      Expression *expr1 = new SCIP_ge(sum1, new SCIP_IntVar(-M, -M));
      solver->add_scip_expr(expr1);
      
      //std::cout << "About to add constraints" << std::endl;
      
      expr1->add(solver, true);
      
      //std::cout << "Starting equation 2" << std::endl;
      
      // equation 2
      Expression *expr2 = new SCIP_lt(sum1, 0.0);
      solver->add_scip_expr(expr2);
      
      expr2->add(solver, true);
      
      // Link the variables
      _var = b->_var;
    }
  }
}

/* Lt object */

SCIP_lt::SCIP_lt(Expression *var1, Expression *var2)
  : SCIP_binop(var1,var2)
{
}

SCIP_lt::SCIP_lt(Expression *var1, double rhs)
  : SCIP_binop(var1,rhs)
{
}

SCIP_lt::~SCIP_lt(){
}

void SCIP_lt::add(SCIPSolver *solver, bool top_level){
  if(!has_been_added()){

    _scip = solver->get_scip();
      
    // Add children
    _vars[0]->add(solver, false);
  
    if(top_level){
    
      if(_is_proper_coef){
      
	SCIP_VAR *vars[1];
      
	vars[0] = _vars[0]->_var;

	SCIP_CONS *cons;
      
	double var_vals[1];
	var_vals[0] = 1.0;
    
	SCIP_CALL_EXC( SCIPcreateConsLinear(_scip, &cons, "Less than constraint",
					    1, // # vars
					    vars, // variables
					    var_vals, // values
					    -SCIPinfinity(_scip), // LHS
					    _rhs - 1.0, // RHS
					    TRUE, TRUE, TRUE,
					    TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
	SCIP_CALL_EXC( SCIPaddCons(_scip, cons) );

	solver->add_scip_cons(cons);
      
      } else {
  
	_vars[1]->add(solver, false);
	SCIP_VAR *vars[2];
      
	vars[0] = _vars[0]->_var;
	vars[1] = _vars[1]->_var;

#ifdef _DEBUGWRAP      
	if(vars[0] == NULL)
	  std::cout << "First one is NULL" << std::endl;
	if(vars[1] == NULL)
	  std::cout << "Second one is NULL" << std::endl;
#endif
      
	SCIP_CONS *cons;
      
	double var_vals[2];
	var_vals[0] = 1.0;
	var_vals[1] = -1.0;
    
	SCIP_CALL_EXC( SCIPcreateConsLinear(_scip, &cons, "Less than constraint",
					    2, // # vars
					    vars, // variables
					    var_vals, // values
					    -SCIPinfinity(_scip), // LHS
					    -1.0, // RHS
					    TRUE, TRUE, TRUE,
					    TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
	SCIP_CALL_EXC( SCIPaddCons(_scip, cons) );

	solver->add_scip_cons(cons);

      }

#ifdef _DEBUGWRAP    
      std::cout << "create Lt predicate" << std::endl;
#endif

    } else {
      std::cout << "Lt operator not in top level insupported at the moment" << std::endl;
    }
  }
}

/* Gt object */

SCIP_gt::SCIP_gt(Expression *var1, Expression *var2)
  : SCIP_binop(var1,var2)
{
}

SCIP_gt::SCIP_gt(Expression *var1, double rhs)
  : SCIP_binop(var1,rhs)
{
}

SCIP_gt::~SCIP_gt(){
}

void SCIP_gt::add(SCIPSolver *solver, bool top_level){
  if(!has_been_added()){

#ifdef _DEBUGWRAP
    std::cout << "\t\t\t In add method of gt now" << std::endl;

#endif


    _scip = solver->get_scip();

    if(top_level){
    
      if(_is_proper_coef){
      
	 _vars[0]->add(solver, false);
      
	SCIP_VAR *vars[1];
      
	vars[0] = _vars[0]->_var;
      
	SCIP_CONS *cons;
      
	double var_vals[1];
	var_vals[0] = 1.0;
    
	SCIP_CALL_EXC( SCIPcreateConsLinear(_scip, &cons, "Greater than constraint",
					    1, // # vars
					    vars, // variables
					    var_vals, // values
					    _rhs+1, // LHS
					    SCIPinfinity(_scip), // RHS
					    TRUE, TRUE, TRUE,
					    TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
	SCIP_CALL_EXC( SCIPaddCons(_scip, cons) );

	solver->add_scip_cons(cons);
      
      } else {
  
	SCIP_VAR *vars[2];
      
      
       _vars[0]->add(solver, false);
	 _vars[1]->add(solver, false);
    
	vars[0] = _vars[0]->_var;
	vars[1] = _vars[1]->_var;

#ifdef _DEBUGWRAP      
	if(vars[0] == NULL)
	  std::cout << "First one is NULL" << std::endl;
	if(vars[1] == NULL)
	  std::cout << "Second one is NULL" << std::endl;
#endif
      
	SCIP_CONS *cons;
      
	double var_vals[2];
	var_vals[0] = 1.0;
	var_vals[1] = -1.0;
    
	SCIP_CALL_EXC( SCIPcreateConsLinear(_scip, &cons, "Greater than constraint",
					    2, // # vars
					    vars, // variables
					    var_vals, // values
					    1, // LHS
					    SCIPinfinity(_scip), // RHS
					    TRUE, TRUE, TRUE,
					    TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
	SCIP_CALL_EXC( SCIPaddCons(_scip, cons) );

	solver->add_scip_cons(cons);
    
      }

#ifdef _DEBUGWRAP
      std::cout << "create Gt predicate" << std::endl;
#endif

    } else {
      std::cout << "Gt operator not in top level insupported at the moment" << std::endl;
    }
  }
}

/*
 *
 * Boolean operator stuff
 *
 */

SCIP_and::SCIP_and(Expression *var1, Expression *var2):
  SCIP_binop(var1, var2){
  // Nothing else to do?
}

SCIP_and::~SCIP_and(){
  // Getting destroyed
}

void SCIP_and::add(SCIPSolver *solver, bool top_level){

  _vars[0]->add(solver, false);
  _vars[1]->add(solver, false);

  if(top_level){
    // We are in the top level
    SCIPExpArray *exprs = new SCIPExpArray();
    exprs->add(_vars[0]);
    exprs->add(_vars[1]);
    solver->add_scip_var_array(exprs);

    SCIPIntArray *coefs = new SCIPIntArray();
    coefs->add(1);
    coefs->add(1);
    solver->add_scip_int_array(coefs);

    Expression *sum = new SCIP_Sum(*exprs, *coefs);
    solver->add_scip_expr(sum);

    Expression *constraint = new SCIP_eq(sum, 2.0);
    constraint->add(solver, true);// Cause it's gotta be true
    solver->add_scip_expr(constraint);

  } else {

    Expression *c = new SCIP_IntVar(0, 1);
    solver->add_scip_expr(c);
    
    SCIPExpArray *exprs = new SCIPExpArray();
    exprs->add(_vars[0]);
    exprs->add(_vars[1]);
    solver->add_scip_var_array(exprs);

    SCIPIntArray *coefs = new SCIPIntArray();
    coefs->add(1);
    coefs->add(1);
    solver->add_scip_int_array(coefs);

    SCIPExpArray *exprs2 = new SCIPExpArray();
    exprs2->add(c);
    solver->add_scip_var_array(exprs2);

    SCIPIntArray *coefs2 = new SCIPIntArray();
    coefs2->add(2);
    solver->add_scip_int_array(coefs2);

    Expression *constraint2 = new SCIP_ge(new SCIP_Sum(*exprs, *coefs),
					       new SCIP_Sum(*exprs2, *coefs2));
    constraint2->add(solver, true);
    solver->add_scip_expr(constraint2);

    Expression *constraint3 = new SCIP_ge(new SCIP_add(c, 1.0),
					       new SCIP_Sum(*exprs, *coefs));
    constraint3->add(solver, true);
    solver->add_scip_expr(constraint3);
      
    _var = c->_var;
    _upper = 1.0;
    _lower = 0.0;
  }
}

SCIP_or::SCIP_or(Expression *var1, Expression *var2):
  SCIP_binop(var1, var2){
  // Nothing else to do?
}

SCIP_or::~SCIP_or(){
  // Getting destroyed
}

void SCIP_or::add(SCIPSolver *solver, bool top_level){

  _vars[0]->add(solver, false);
  _vars[1]->add(solver, false);

  if(top_level){
    
    // We are in the top level
    SCIPExpArray *exprs = new SCIPExpArray();
    exprs->add(_vars[0]);
    exprs->add(_vars[1]);
    solver->add_scip_var_array(exprs);

    SCIPIntArray *coefs = new SCIPIntArray();
    coefs->add(1);
    coefs->add(1);
    solver->add_scip_int_array(coefs);

    Expression *constraint = new SCIP_ge(new SCIP_Sum(*exprs, *coefs), 1.0);
    constraint->add(solver, true);
    solver->add_scip_expr(constraint);

  } else {

    Expression *c = new SCIP_IntVar(0, 1);
    solver->add_scip_expr(c);
    
    SCIPExpArray *exprs = new SCIPExpArray();
    exprs->add(_vars[0]);
    exprs->add(_vars[1]);
    exprs->add(c);
    solver->add_scip_var_array(exprs);

    SCIPIntArray *coefs = new SCIPIntArray();
    coefs->add(1);
    coefs->add(1);
    coefs->add(-1);
    solver->add_scip_int_array(coefs);
        
    Expression *constraint = new SCIP_ge(new SCIP_Sum(*exprs, *coefs), 0.0);
    constraint->add(solver, true);
    solver->add_scip_expr(constraint);

    _var = c->_var;
    _upper = 1.0;
    _lower = 0.0;
  }
}

SCIP_not::SCIP_not(Expression *var1):
  SCIP_binop(var1, 0.0){
  // Nothing else to do?
}

SCIP_not::~SCIP_not(){
  // getting destroyed
}

void SCIP_not::add(SCIPSolver *solver, bool top_level){
  _vars[0]->add(solver, false);

  if(top_level){
    
    Expression *constraint = new SCIP_eq(_vars[0], 0.0);
    constraint->add(solver, true);

  } else {

    Expression *c = new SCIP_IntVar(0, 1);

    SCIPExpArray *exprs = new SCIPExpArray();
    exprs->add(_vars[0]);
    exprs->add(c);
    solver->add_scip_var_array(exprs);
    
    SCIPIntArray *coefs = new SCIPIntArray();
    coefs->add(1);
    coefs->add(1);
    solver->add_scip_int_array(coefs);
    

    Expression *constraint = new SCIP_eq(new SCIP_Sum(*exprs, *coefs), 1);
    constraint->add(solver, true);
    solver->add_scip_expr(constraint);

    _var = c->_var;
    _upper = 1.0;
    _lower = 0.0;
  }
} 


// Minimise Class

SCIP_Minimise::SCIP_Minimise(Expression *arg1):
SCIP_eq(arg1, new SCIP_FloatVar(arg1->_lower, arg1->_upper)){
  _vars[1]->_coef = -1;
}

SCIP_Minimise::~SCIP_Minimise(){
}

// Maximise Class

SCIP_Maximise::SCIP_Maximise(Expression *arg1):
SCIP_eq(arg1, new SCIP_FloatVar(arg1->_lower, arg1->_upper)){
  _vars[1]->_coef = 1;
}

SCIP_Maximise::~SCIP_Maximise(){
}

// Expression* SCIP_Minimise(Expression* arg)
// {
//   std::cout << "create a minimisation ective" << std::endl;
  
//   arg->add(solver, false);
  
//   SCIP_Variable *var = new SCIP_Variable(arg->_lower, arg->_upper);
//   var->_coef = -1;

//   return new SCIP_eq(arg, var);
// }

// Expression* SCIP_Maximise(Expression* arg)
// {
//   std::cout << "create a maximisation objective" << std::endl;
  
//   arg->add(solver, false);
  
//   SCIP_Variable *var = new SCIP_Variable(arg->_lower, arg->_upper);
//   var->_coef = 1;
  
//   return new SCIP_eq(arg, var);
// }

/**************************************************************
 ********************     Solver        ***********************
 **************************************************************/

SCIPSolver::SCIPSolver()
{

#ifdef _DEBUGWRAP
  std::cout << "create a scip solver" << std::endl;
#endif

  var_counter = 0;

  // Create the scip object
  
  // Load up SCIP
  SCIP_CALL_EXC( SCIPcreate(& _scip) );

  // load default plugins linke separators, heuristics, etc.
  SCIP_CALL_EXC( SCIPincludeDefaultPlugins(_scip) );

  // create an empty problem
  SCIP_CALL_EXC( SCIPcreateProb(_scip, "PPP", NULL, NULL, NULL, NULL, NULL, NULL) );

  // set the objective sense to maximize, default is minimize
  SCIP_CALL_EXC( SCIPsetObjsense(_scip, SCIP_OBJSENSE_MAXIMIZE) );
}

SCIPSolver::~SCIPSolver()
{

#ifdef _DEBUGWRAP
  std::cout << "delete wrapped solver" << std::endl;
#endif

  for(unsigned int i=0; i<_scip_vars.size(); ++i)
    SCIP_CALL_EXC( SCIPreleaseVar(_scip, &_scip_vars[i]) );

  for(unsigned int i=0; i<_scip_cons.size(); ++i)
    SCIP_CALL_EXC( SCIPreleaseCons(_scip, &_scip_cons[i]) );
    
  for(unsigned int i=0; i < _scip_exprs.size(); ++i)
    delete(_scip_exprs[i]);
    
  for(unsigned int i=0; i < _int_array.size(); ++i)
    delete(_int_array[i]);
    
  for(unsigned int i=0; i < _var_array.size(); ++i)
    delete(_var_array[i]);

  SCIP_CALL_EXC( SCIPfree( &_scip) );

  _verbosity = 0;
}

void SCIPSolver::add(Expression* arg)
{

#ifdef _DEBUGWRAP
  std::cout << "add an expression" << std::endl;
#endif

  if(arg != NULL)
    arg->add(this, true);
}

void SCIPSolver::add_scip_var (SCIP_VAR * v) {_scip_vars.push_back(v);}
void SCIPSolver::add_scip_cons(SCIP_CONS* c) {_scip_cons.push_back(c);}

void SCIPSolver::add_scip_expr(Expression *expr){
  _scip_exprs.push_back(expr);
}

void SCIPSolver::add_scip_int_array(SCIPIntArray *arr){
  _int_array.push_back(arr);
}

void SCIPSolver::add_scip_var_array(SCIPExpArray *arr){
  _var_array.push_back(arr);
}

SCIP* SCIPSolver::get_scip() {return _scip;}

void SCIPSolver::initialise(SCIPExpArray& arg)
{

#ifdef _DEBUGWRAP
  std::cout << "initialise the solver" << std::endl;
#endif

}

void SCIPSolver::initialise()
{

#ifdef _DEBUGWRAP
  std::cout << "initialise the solver" << std::endl;
#endif

}

int SCIPSolver::solve()
{

#ifdef _DEBUGWRAP
  std::cout << "solve!" << std::endl;  
  SCIP_CALL_EXC(SCIPprintOrigProblem(_scip, NULL, NULL, FALSE));
#endif

  //std::cout << "c solve using SCIP " << std::endl;

  if(_verbosity == 1){
    // Do nothing extra
  } else if(_verbosity == 2){
    SCIP_CALL_EXC(SCIPprintOrigProblem(_scip, NULL, NULL, FALSE));
  } else {
      // disable scip output to stdout
    SCIP_CALL_EXC( SCIPsetMessagehdlr(NULL) );
  }

  SCIP_CALL_EXC( SCIPsolve(_scip) );
 
  SCIP_STATUS status = SCIPgetStatus(_scip);
 
  if( status == SCIP_STATUS_OPTIMAL ) return SAT;
  else if( status == SCIP_STATUS_INFEASIBLE ) return UNSAT;
  else return UNKNOWN;
}

int SCIPSolver::solveAndRestart(const int policy,
				const unsigned int base,
				const double factor,
				const double decay )
{
  return solve();
}

int SCIPSolver::startNewSearch()
{
  std::cout << "c start a new interuptable search" << std::endl;
  return 0;
}

int SCIPSolver::getNextSolution()
{
  std::cout << "c seek next solution" << std::endl;
  return 0;
}

int SCIPSolver::sacPreprocess(const int type)
{
  std::cout << "c enforces singleton arc consistency" << std::endl;
  return 0;
}

void SCIPSolver::setHeuristic(const char* var_heuristic, const char* val_heuristic, const int rand)
{
  std::cout << "c set the variable/value ordering (ignored)" << std::endl;
}

void SCIPSolver::setFailureLimit(const int cutoff)
{
  std::cout << "c set a cutoff on failures" << std::endl;
}

void SCIPSolver::setTimeLimit(const int cutoff)
{
  //std::cout << "set a cutoff on time" << std::endl;
  //SCIPparamSetReal (_scip->, _scip, (double)cutoff)
  //SCIP_PARAM* param = (SCIP_PARAM*)SCIPdialogGetData(time);
  SCIPsetRealParam(_scip, "limits/time", (double)cutoff);
}

void SCIPSolver::setVerbosity(const int degree)
{
#ifdef _DEBUGWRAP
  std::cout << "c set the verbosity" << std::endl;
#endif
  _verbosity = degree;
}

void SCIPSolver::setRandomized(const int degree)
{
  std::cout << "c set the type of randomization" << std::endl;
}

void SCIPSolver::setRandomSeed(const int seed)
{
  std::cout << "c set the random seed" << std::endl;
}

bool SCIPSolver::is_sat()
{
  
  SCIP_STATUS status = SCIPgetStatus(_scip);
  return !( status == SCIP_STATUS_INFEASIBLE );

  //std::cout << "whether the problem was satisfiable" << std::endl;

}

bool SCIPSolver::is_unsat()
{
  //std::cout << "whether the problem was unsatisfiable" << std::endl;
  return ! is_sat();
}

void SCIPSolver::printStatistics()
{
#ifdef _DEBUGWRAP
  std::cout << "c print a bunch of statistics" << std::endl;
#endif
  std::cout << "\td Time: " << getTime() << "\tNodes:" << getNodes() << std::endl;
  
  if(_verbosity == 2){
    SCIP_CALL_EXC(SCIPprintStatistics(_scip, NULL));
  }
}

int SCIPSolver::getBacktracks()
{
  std::cout << "c print the number of backtracks" << std::endl;
  return 0;
}

int SCIPSolver::getNodes()
{
#ifdef _DEBUGWRAP
  std::cout << "c Print the number of nodes" << std::endl;
#endif
  
  return _scip->stat->nnodes;
}

int SCIPSolver::getFailures()
{
  std::cout << "c print the number of failures" << std::endl;
  return 0;
}

int SCIPSolver::getChecks()
{
  std::cout << "c print the number of checks" << std::endl;
  return 0;
}

int SCIPSolver::getPropags()
{
  std::cout << "c print the number of propags" << std::endl;
  return 0;
}

double SCIPSolver::getTime()
{
#ifdef _DEBUGWRAP
  std::cout << "c print the cpu time" << std::endl;
#endif
  return SCIPclockGetTime(_scip->stat->solvingtime);
}

/**
 *
 *
 *
 * Stuff to improve SCIP
 *
 *
 *
 */


/**
 * Creates an empty linear constraint object
 */
LinearConstraint:: LinearConstraint(double lhs, double rhs):
  _lhs(lhs), _rhs(rhs){
}
    
/**
 * Destructor
 */
LinearConstraint::~LinearConstraint(){}
    
/**
 * Add in a normal coefficient
 */
void LinearConstraint::add_coef(Expression* expr){
  _variables.push_back(expr);
  _coefficients.push_back(expr->_coef);
}
    
/**
 * Add in a coefficient that is a view (some reformulation required)
 */
void LinearConstraint::add_coef(LinearView* expr){
  
}
    
/**
 * Prints the linear constraint
 */
void LinearConstraint::display(){
  std::cout << _lhs << " <= ";
  for(int i = 0; i < _variables.size(); ++i){
    _variables[i]->display();
    if(i + 1 < _variables.size()) std::cout << "+ ";
  }
  std::cout << " <= " << _rhs << std::endl;
}