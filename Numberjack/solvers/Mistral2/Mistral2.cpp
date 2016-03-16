
/** \file Mistral2.cpp
    \brief Solver interface for PYTHON Wrapper.
*/

#include "Mistral2.hpp"

/**************************************************************
 ********************     EXPRESSION        *******************
 **************************************************************/

//#define _DEBUGWRAP true

Mistral2_Expression::Mistral2_Expression()
{

#ifdef _DEBUGWRAP
  std::cout << "creating a Boolean expression" << std::endl;
#endif 

  Mistral::Variable x(0,1);
  _self = x;
}

Mistral2_Expression::Mistral2_Expression(const int nval)
{

#ifdef _DEBUGWRAP
  std::cout << "creating a variable [" << 0 << ".." << (nval-1) << "]" << std::endl;
#endif 

  Mistral::Variable x(0,nval-1);
  _self = x;
}

Mistral2_Expression::Mistral2_Expression(const int lb, const int ub)
{

#ifdef _DEBUGWRAP
  std::cout << "creating a variable [" << lb << ".." << ub << "]" << std::endl;
#endif 

  Mistral::Variable x(lb,ub);
  _self = x;
}

Mistral2_Expression::Mistral2_Expression(Mistral2IntArray& vals)
{

#ifdef _DEBUGWRAP
  std::cout << "creating a domain variable" << std::endl;
#endif 

  Mistral::Vector<int> util;
  for(int i=0; i<vals.size(); ++i) {
    util.add(vals.get_item(i));
  }
  Mistral::Variable x(util);
  _self = x;
}


int Mistral2_Expression::getVariableId() const
{

#ifdef _DEBUGWRAP
  std::cout << "return identity of expression" << std::endl;
#endif 

  return _self.id() ;
}

int Mistral2_Expression::get_value() const
{

#ifdef _DEBUGWRAP
  std::cout << "return value of expression" << std::endl;
#endif 

  return _self.get_solution_int_value() ;
}

int Mistral2_Expression::get_size() const
{

#ifdef _DEBUGWRAP
  std::cout << "return size of expression" << std::endl;
#endif 

  return _self.get_size();
}

int Mistral2_Expression::get_max() const
{

#ifdef _DEBUGWRAP
  std::cout << "return max of expression" << std::endl;
#endif
  return _self.get_max();
}

int Mistral2_Expression::get_min() const
{

#ifdef _DEBUGWRAP
  std::cout << "return min of expression" << std::endl;
#endif

  return _self.get_min();
}

bool Mistral2_Expression::contain(const int v) const
{

#ifdef _DEBUGWRAP
  std::cout << "return min of expression" << std::endl;
#endif
  return _self.contain(v);
}

Mistral2_Expression::~Mistral2_Expression()
{

#ifdef _DEBUGWRAP
  std::cout << "delete expression" << std::endl;
#endif
}

bool Mistral2_Expression::has_been_added() const
{
  //std::cout << "has this expression already been added?" << std::endl;
  return _self.is_initialised();
}

Mistral2_Expression* Mistral2_Expression::add(Mistral2Solver *solver, bool top_level)
{

#ifdef _DEBUGWRAP
  std::cout << "add expression to model" << std::endl;
#endif
  
  if(top_level){
#ifdef _DEBUGWRAP
    std::cout << "\tAdding at top level" << std::endl;
#endif
  } else {
#ifdef _DEBUGWRAP
    std::cout << "\tAdding within tree" <<std::endl;
#endif
  }

  if(!has_been_added()) {
    _solver = solver;
    if(top_level) {
      _solver->add(this);
    }
  }

  return this;
}

// /* Binary operators */

Mistral2_binop::Mistral2_binop(Mistral2_Expression *var1,
			       Mistral2_Expression *var2)
  : Mistral2_Expression()
{

#ifdef _DEBUGWRAP
  std::cout << "creating a binary operator" << std::endl;
#endif
  _vars[0] = var1;
  _vars[1] = var2;
}

Mistral2_binop::Mistral2_binop(Mistral2_Expression *var1, int constant)
  : Mistral2_Expression()
{

#ifdef _DEBUGWRAP
  std::cout << "creating a binary (constant) operator" << std::endl;
#endif
  _vars[0] = var1;
  _vars[1] = NULL;
  _constant = constant;
}

// Mistral2_binop::Mistral2_binop(int constant, Mistral2_Expression *var1)
//   : Mistral2_Expression()
// {

// #ifdef _DEBUGWRAP
//   std::cout << "creating a binary (constant) operator" << std::endl;
// #endif
//   _vars[0] = NULL;
//   _vars[1] = var1;
//   _constant = constant;
// }


Mistral2_binop::~Mistral2_binop(){

#ifdef _DEBUGWRAP
  std::cout << "delete binary operator" << std::endl;
#endif
}

/**
 * Constraints 
 */

Mistral2_Min::Mistral2_Min( Mistral2ExpArray& vars ) 
  : Mistral2_Expression() 
{

#ifdef _DEBUGWRAP
  std::cout << "creating an Min constraint" << std::endl;
#endif

_vars = vars;

}

Mistral2_Min::Mistral2_Min( Mistral2_Expression *var1, Mistral2_Expression *var2 ) 
  : Mistral2_Expression() 
{

#ifdef _DEBUGWRAP
  std::cout << "creating a binary Min constraint" << std::endl;
#endif

  _vars.add(var1);
  _vars.add(var2); 
}

Mistral2_Min::~Mistral2_Min()
{

#ifdef _DEBUGWRAP
  std::cout << "delete Min" << std::endl;
#endif
}

Mistral2_Expression* Mistral2_Min::add(Mistral2Solver *solver, bool top_level)
{
  if(!has_been_added()) {

#ifdef _DEBUGWRAP
    std::cout << "add a Min constraint to solver" << std::endl;
#endif

    _solver = solver;
    
    int i, n=_vars.size();  
    for(i=0; i<n; ++i) 
      _vars.get_item(i)->add(_solver,false);

    if(n == 2) {
      _self = Min(_vars.get_item(0)->_self,
		  _vars.get_item(1)->_self);
    } else if(n > 2) {
      Mistral::VarArray scope(n);
      for(i=0; i<n; ++i) scope[i] = _vars.get_item(i)->_self;
      _self = Min(scope);
    }

    if( top_level )
      _solver->solver->add( _self );    
  
  }
  return this;
}


Mistral2_Max::Mistral2_Max( Mistral2ExpArray& vars ) 
  : Mistral2_Expression() 
{

#ifdef _DEBUGWRAP
  std::cout << "creating a Max constraint" << std::endl;
#endif

  _vars = vars;
}

Mistral2_Max::Mistral2_Max( Mistral2_Expression *var1, Mistral2_Expression *var2 ) 
  : Mistral2_Expression() 
{

#ifdef _DEBUGWRAP
  std::cout << "creating a binary Max constraint" << std::endl;
#endif

  _vars.add(var1);
  _vars.add(var2);
}

Mistral2_Max::~Mistral2_Max()
{

#ifdef _DEBUGWRAP
  std::cout << "delete Max" << std::endl;
#endif
}

Mistral2_Expression* Mistral2_Max::add(Mistral2Solver *solver, bool top_level)
{
  if(!has_been_added()) {

#ifdef _DEBUGWRAP
    std::cout << "add a Max constraint" << std::endl;
#endif

    _solver = solver;    
    int i, n=_vars.size(); 
		 
    for(i=0; i<n; ++i) {
      _vars.get_item(i)->add(_solver,false);
		}

    if(n == 2) {
      _self = Max(_vars.get_item(0)->_self,
		  _vars.get_item(1)->_self);
    } else if(n > 2) {
      Mistral::VarArray scope(n);
      for(i=0; i<n; ++i) scope[i] = _vars.get_item(i)->_self;
      _self = Max(scope);
    }

    if( top_level )
      _solver->solver->add( _self );
    
  }
  return this;
}

Mistral2_AllDiff::Mistral2_AllDiff( Mistral2ExpArray& vars ) 
  : Mistral2_Expression() 
{
#ifdef _DEBUGWRAP
  std::cout << "creating an alldiff constraint" << std::endl;
#endif
  _vars = vars;
}

Mistral2_AllDiff::Mistral2_AllDiff( Mistral2_Expression *var1, Mistral2_Expression *var2 ) 
  : Mistral2_Expression() 
{
#ifdef _DEBUGWRAP
  std::cout << "creating a binary alldiff constraint" << std::endl;
#endif
  _vars.add(var1);
  _vars.add(var2); 
}

Mistral2_AllDiff::~Mistral2_AllDiff()
{
#ifdef _DEBUGWRAP
  std::cout << "delete alldiff" << std::endl;
#endif
}

Mistral2_Expression* Mistral2_AllDiff::add(Mistral2Solver *solver, bool top_level)
{
  if(!has_been_added()) {
#ifdef _DEBUGWRAP
    std::cout << "add an alldiff constraint" << std::endl;
#endif
    
    _solver = solver;

    int i, n=_vars.size();  
    for(int i = 0; i < n; ++i)
      _vars.set_item(i, _vars.get_item(i)->add(solver, false));
    
    if(n==2) {
      _self = (_vars.get_item(0)->_self != _vars.get_item(1)->_self);
    } else {
      Mistral::VarArray scope(n);
      for(i=0; i<n; ++i) scope[i] = _vars.get_item(i)->_self;
      _self = AllDiff(scope);
    }

    if(top_level){
#ifdef _DEBUGWRAP
      std::cout << "\tAdding at top level" << std::endl;
#endif
      _solver->solver->add( _self );
    } else {
#ifdef _DEBUGWRAP
      std::cout << "\tAdding within tree AllDiff constraint NOT A GOOD IDEA" <<std::endl;
#endif
      exit(1);
    }
    
  }
  return this;
}

Mistral2_Gcc::Mistral2_Gcc(Mistral2ExpArray& vars,
			   Mistral2IntArray& vals,
			   Mistral2IntArray& lb_card,
			   Mistral2IntArray& ub_card)
  : Mistral2_Expression() 
{
#ifdef _DEBUGWRAP
  std::cout << "creating a gcc constraint" << std::endl;
#endif

  _vars = vars;
  _vals = vals;
  _lb_card = lb_card;
  _ub_card = ub_card;

}

Mistral2_Gcc::~Mistral2_Gcc()
{
#ifdef _DEBUGWRAP
  std::cout << "delete gcc" << std::endl;
#endif
}

Mistral2_Expression* Mistral2_Gcc::add(Mistral2Solver *solver, bool top_level)
{
  if(!has_been_added()) {
#ifdef _DEBUGWRAP
    std::cout << "add gcc constraint" << std::endl;
#endif
      
    _solver = solver;
    
    int i, n=_vars.size(), m=_vals.size();
    int min_val=_vals.get_item(0);
    int max_val=_vals.get_item(m-1);
    int M = (max_val - min_val + 1);
    Mistral::VarArray scope(n);

    int *tmp_lb = new int[M];
    int *tmp_ub = new int[M];

    for(i=0; i<n; ++i) {
      _vars.get_item(i)->add(_solver,false);
      scope[i] = _vars.get_item(i)->_self;
    }

    //std::cout << min_val << " to " << max_val << std::endl;

    for(i=0; i<M; ++i) {
      tmp_lb[i] = _lb_card.get_item(i);
      tmp_ub[i] = _ub_card.get_item(i);

      //std::cout << " " << i+min_val << ": in [" << tmp_lb[i] << ".." << tmp_ub[i] << "]\n";
    }

    _self = Occurrences(scope, min_val, max_val, tmp_lb, tmp_ub);

    if( top_level )
      _solver->solver->add( _self );

    delete [] tmp_lb;
    delete [] tmp_ub;

  }

  return this;
}

Mistral2_Element::Mistral2_Element( Mistral2ExpArray& vars ) 
  : Mistral2_Expression() 
{
#ifdef _DEBUGWRAP
  std::cout << "creating element" << std::endl;
#endif

  _vars = vars;

}

Mistral2_Element::~Mistral2_Element()
{
#ifdef _DEBUGWRAP
  std::cout << "delete element" << std::endl;
#endif
}

Mistral2_Expression* Mistral2_Element::add(Mistral2Solver *solver, bool top_level)
{
  if(!has_been_added()) {
#ifdef _DEBUGWRAP
      std::cout << "add element constraint" << std::endl;
#endif
    _solver = solver;
    
    int i, n=_vars.size();  
    Mistral::VarArray scope(n-1);

    for(i=0; i<n-1; ++i) {
      _vars.get_item(i)->add(_solver,false);
      scope[i] = _vars.get_item(i)->_self;
    }	
    
    _vars.get_item(n-1)->add(_solver,false);
    Mistral::Variable index = _vars.get_item(n-1)->_self;
		
    _self = scope[index];

  }

  return this;
}

Mistral2_LeqLex::Mistral2_LeqLex( Mistral2ExpArray& vars ) 
  : Mistral2_Expression() 
{
#ifdef _DEBUGWRAP
  std::cout << "creating lexleq" << std::endl;
#endif
}

Mistral2_LeqLex::~Mistral2_LeqLex()
{
#ifdef _DEBUGWRAP
  std::cout << "delete leqlex" << std::endl;
#endif
}

Mistral2_Expression* Mistral2_LeqLex::add(Mistral2Solver *solver, bool top_level) {
  if(!has_been_added()) {
#ifdef _DEBUGWRAP
    std::cout << "add leqlex constraint" << std::endl;
#endif

    _solver = solver;
      
    for(int i = 0; i < _vars.size(); ++i)
      _vars.set_item(i, _vars.get_item(i)->add(solver, false));
      
    if(top_level){
#ifdef _DEBUGWRAP
      std::cout << "\tAdding at top level" << std::endl;
#endif
    } else {
#ifdef _DEBUGWRAP
      std::cout << "\tAdding within tree" <<std::endl;
#endif
    }
      
  }
  return this;
}


Mistral2_LessLex::Mistral2_LessLex( Mistral2ExpArray& vars ) 
  : Mistral2_Expression() 
{
#ifdef _DEBUGWRAP
  std::cout << "creating lexless" << std::endl;
#endif
}

Mistral2_LessLex::~Mistral2_LessLex()
{
#ifdef _DEBUGWRAP
  std::cout << "delete leslex" << std::endl;
#endif
}

Mistral2_Expression* Mistral2_LessLex::add(Mistral2Solver *solver, bool top_level) {
  if(!has_been_added()) {
#ifdef _DEBUGWRAP
    std::cout << "add lesslex constraint" << std::endl;
#endif

    _solver = solver;
      
    for(int i = 0; i < _vars.size(); ++i)
      _vars.set_item(i, _vars.get_item(i)->add(solver, false));
      
    if(top_level){
#ifdef _DEBUGWRAP
      std::cout << "\tAdding at top level" << std::endl;
#endif
    } else {
#ifdef _DEBUGWRAP
      std::cout << "\tAdding within tree" <<std::endl;
#endif
    }
      
  }
  return this;
}
 

Mistral2_Sum::Mistral2_Sum(Mistral2ExpArray& vars, 
			   Mistral2IntArray& weights, 
			   const int offset)
  : Mistral2_Expression() 
{
#ifdef _DEBUGWRAP
  std::cout << "creating sum" << std::endl;
#endif
  _vars = vars;
  _weights = weights;
  _offset = offset;
}

Mistral2_Sum::Mistral2_Sum(Mistral2_Expression *arg1, 
			   Mistral2_Expression *arg2, 
			   Mistral2IntArray& w, 
			   const int offset)
  : Mistral2_Expression() 
{
#ifdef _DEBUGWRAP
  std::cout << "creating sum" << std::endl;
#endif
  _vars.add(arg1);
  _vars.add(arg2); 
  _weights = w;
  _offset = offset;
}

Mistral2_Sum::Mistral2_Sum(Mistral2_Expression *arg, 
			   Mistral2IntArray& w, 
			   const int offset)
  : Mistral2_Expression() 
{
#ifdef _DEBUGWRAP
  std::cout << "creating sum" << std::endl;
#endif
  _vars.add(arg);
  _weights = w;
  _offset = offset;
}

Mistral2_Sum::Mistral2_Sum()
  : Mistral2_Expression()
{
  _offset = 0;
}

Mistral2_Sum::~Mistral2_Sum(){
#ifdef _DEBUGWRAP
  std::cout << "delete sum" << std::endl;
#endif
}

void Mistral2_Sum::addVar(Mistral2_Expression* v){
#ifdef _DEBUGWRAP
  std::cout << "adding variable" << std::endl;
#endif
}

void Mistral2_Sum::addWeight(const int w){
#ifdef _DEBUGWRAP
  std::cout << "adding weight" << std::endl;
#endif
}

Mistral2_Expression* Mistral2_Sum::add(Mistral2Solver *solver, bool top_level){
  if(!has_been_added()) {
#ifdef _DEBUGWRAP
    std::cout << "add sum constraint" << std::endl;
#endif
    _solver = solver;
    
    // _v * _w + o = _x
    // _v * _w - _x = -o 

    
    int i, n=_vars.size();  
    Mistral::VarArray scope(n);
    Mistral::Vector<int> w;
    bool weighted = false;
    bool boolean = true;
    w.initialise(n, n);

    for(i=0; i<n; ++i) {

      _vars.get_item(i)->add(_solver,false);
      if(_weights.size() > i)
	w[i] = _weights.get_item(i);
      else
	w[i] = 1;
      scope[i] = _vars.get_item(i)->_self;

      if(w[i] != 1) weighted = true;
      if(!scope[i].is_boolean()) boolean = false;
    }
    //w[n] = _weights.get_item(n);
 

    if(boolean) {
      if(weighted)
	_self = BoolSum(scope, w);
      else
	_self = BoolSum(scope, w);
    } else {
      if(weighted)
	_self = Sum(scope, w, -Mistral::INFTY, Mistral::INFTY, _offset);
      else
	_self = Sum(scope, w, -Mistral::INFTY, Mistral::INFTY, _offset);
    }
    
    if( top_level ) {
      _solver->solver->add( _self );
    } 
  }

  return this;
}



Mistral2_OrderedSum::Mistral2_OrderedSum(Mistral2ExpArray& vars, 
			   const int l, const int u)
  : Mistral2_Expression() 
{
#ifdef _DEBUGWRAP
  std::cout << "creating ordered sum" << std::endl;
#endif
  _vars = vars;
  _lb = l;
  _ub = u;
}


Mistral2_OrderedSum::~Mistral2_OrderedSum(){
#ifdef _DEBUGWRAP
  std::cout << "delete ordered sum" << std::endl;
#endif
}


Mistral2_Expression* Mistral2_OrderedSum::add(Mistral2Solver *solver, bool top_level){
	if(!has_been_added()) {
#ifdef _DEBUGWRAP
		std::cout << "add ordered sum constraint" << std::endl;
#endif
		_solver = solver;
    
		// _v * _w + o = _x
		// _v * _w - _x = -o 
    
		int i, n=_vars.size();  
		Mistral::VarArray scope(n);

		for(i=0; i<n; ++i) {
			_vars.get_item(i)->add(_solver,false);
			scope[i] = _vars.get_item(i)->_self;
		}

		_self = OSum(scope, _lb, _ub);
    
		if( top_level ) {
			_solver->solver->add( _self );
		} 
	}

  return this;
}

/**
 * Unary constraints
 */
Mistral2_Abs::Mistral2_Abs(Mistral2_Expression *var1)
  : Mistral2_binop(var1,0)
{
#ifdef _DEBUGWRAP
  std::cout << "creating Abs constraint" << std::endl;
#endif
}

Mistral2_Abs::~Mistral2_Abs(){
#ifdef _DEBUGWRAP
  std::cout << "delete Abs constraint" << std::endl;
#endif
}

Mistral2_Expression* Mistral2_Abs::add(Mistral2Solver *solver, bool top_level){
  if(!has_been_added()) {
#ifdef _DEBUGWRAP
    std::cout << "add Abs constraint" << std::endl;
#endif
      
    _solver = solver;
    
    _vars[0]->add(_solver,false);
    _self = Abs(_vars[0]->_self);

    if( top_level )
      _solver->solver->add( _self );
      
  }
  return this;
}

Mistral2_neg::Mistral2_neg(Mistral2_Expression *var1)
  : Mistral2_binop(var1,0)
{
#ifdef _DEBUGWRAP
  std::cout << "creating neg constraint" << std::endl;
#endif
}

Mistral2_neg::~Mistral2_neg(){
#ifdef _DEBUGWRAP
  std::cout << "delete neg constraint" << std::endl;
#endif
}

Mistral2_Expression* Mistral2_neg::add(Mistral2Solver *solver, bool top_level){
  if(!has_been_added()) {
#ifdef _DEBUGWRAP
    std::cout << "add neg constraint" << std::endl;
#endif
      
    _solver = solver;
    
    _vars[0]->add(_solver,false);
    _self = -(_vars[0]->_self);

    if( top_level )
      _solver->solver->add( _self );
      
  }
  return this;
}


/**
 * Binary constraints
 */
Mistral2_mul::Mistral2_mul(Mistral2_Expression *var1, Mistral2_Expression *var2)
  : Mistral2_binop(var1,var2)
{
#ifdef _DEBUGWRAP
  std::cout << "creating mul constraint between two variables" << std::endl;
#endif
}

Mistral2_mul::Mistral2_mul(Mistral2_Expression *var1, int constant)
  : Mistral2_binop(var1,constant)
{
#ifdef _DEBUGWRAP
  std::cout << "creating mul constraint between variable and constant" << std::endl;
#endif
}


Mistral2_mul::~Mistral2_mul(){
#ifdef _DEBUGWRAP
  std::cout << "delete mul constraint" << std::endl;
#endif
}

Mistral2_Expression* Mistral2_mul::add(Mistral2Solver *solver, bool top_level){
  if(!has_been_added()) {
#ifdef _DEBUGWRAP
    std::cout << "add mul constraint" << std::endl;
#endif
      
    _solver = solver;
    
    _vars[0]->add(_solver,false);
    if(_vars[1]) {
      _vars[1]->add(_solver,false);   
      _self = (_vars[0]->_self *
	       _vars[1]->_self);
    } else {
      _self = (_vars[0]->_self *
	       _constant);
    }
    if( top_level )
      _solver->solver->add( _self );
      
  }
  return this;
}

Mistral2_div::Mistral2_div(Mistral2_Expression *var1, Mistral2_Expression *var2)
  : Mistral2_binop(var1,var2)
{
#ifdef _DEBUGWRAP
  std::cout << "creating div constraint" << std::endl;
#endif
}

Mistral2_div::Mistral2_div(Mistral2_Expression *var1, int constant)
  : Mistral2_binop(var1,constant)
{
#ifdef _DEBUGWRAP
  std::cout << "creating div constraint" << std::endl;
#endif
}


Mistral2_div::~Mistral2_div(){
#ifdef _DEBUGWRAP
  std::cout << "delete div" << std::endl;
#endif
}

Mistral2_Expression* Mistral2_div::add(Mistral2Solver *solver, bool top_level){
  if(!has_been_added()) {
#ifdef _DEBUGWRAP
    std::cout << "add div constraint" << std::endl;
#endif
      
    _solver = solver;
    
    if(_vars[0]) {
      _vars[0]->add(_solver,false);   
      
      if(_vars[1]) {
	_vars[1]->add(_solver,false);   
	_self = (_vars[0]->_self /
		 _vars[1]->_self);
      } else {
	_self = (_vars[0]->_self /
		 _constant);
      }
    } else {

      std::cout << "Cannot divide constants" << std::endl;
      exit(1);

    }
    if( top_level )
      _solver->solver->add( _self );

  }
  return this;
}

Mistral2_mod::Mistral2_mod(Mistral2_Expression *var1, Mistral2_Expression *var2)
  : Mistral2_binop(var1,var2)
{
#ifdef _DEBUGWRAP
  std::cout << "creating mod predicate" << std::endl;
#endif
}

Mistral2_mod::Mistral2_mod(Mistral2_Expression *var1, int constant)
  : Mistral2_binop(var1,constant)
{
#ifdef _DEBUGWRAP
  std::cout << "creating mod predicate" << std::endl;
#endif
}


Mistral2_mod::~Mistral2_mod(){
#ifdef _DEBUGWRAP
  std::cout << "delete mod" << std::endl;
#endif
}

Mistral2_Expression* Mistral2_mod::add(Mistral2Solver *solver, bool top_level){
  if(!has_been_added()) {
#ifdef _DEBUGWRAP
    std::cout << "add mod predicate" << std::endl;
#endif
 
    _solver = solver;
    
    _vars[0]->add(_solver,false);
    if(_vars[1]) {
      _vars[1]->add(_solver,false);   
      _self = (_vars[0]->_self % _vars[1]->_self);
    } else {
      _self = (_vars[0]->_self % _constant);
    }
    if( top_level )
      _solver->solver->add( _self );

  }
  return this;
}

Mistral2_and::Mistral2_and(Mistral2_Expression *var1,
			   Mistral2_Expression *var2)
  : Mistral2_binop(var1,var2)
{
#ifdef _DEBUGWRAP
  std::cout << "creating and predicate" << std::endl;
#endif
}

Mistral2_and::Mistral2_and(Mistral2_Expression *var1,
			   int constant)
  : Mistral2_binop(var1,constant)
{
#ifdef _DEBUGWRAP
  std::cout << "creating and predicate" << std::endl;
  std::cout << "I DON'T THINK I SHOULD BE HERE" << std::endl;
#endif
  
  /**
   * Should never be in this constructor???
   */
  
}

Mistral2_and::~Mistral2_and(){
#ifdef _DEBUGWRAP
  std::cout << "delete and" << std::endl;
#endif
}

Mistral2_Expression* Mistral2_and::add(Mistral2Solver *solver,
				       bool top_level){
  if(!has_been_added()) {
#ifdef _DEBUGWRAP
    std::cout << "add and constraint" << std::endl;
#endif
      
    _solver = solver;
    
    _vars[0]->add(_solver,false);
    if(_vars[1]) {
      _vars[1]->add(_solver,false);   
      _self = (_vars[0]->_self && _vars[1]->_self);
    } else if(_constant) {
      _self = (_vars[0]->_self == 1);
    } else {
      _solver->solver->fail();
    }
    if( top_level )
      _solver->solver->add( _self );
      
  }
  return this;
}

Mistral2_or::Mistral2_or(Mistral2_Expression *var1, Mistral2_Expression *var2)
  : Mistral2_binop(var1,var2)
{
#ifdef _DEBUGWRAP
  std::cout << "creating or predicate" << std::endl;
#endif
}

Mistral2_or::Mistral2_or(Mistral2_Expression *var1, int constant)
  : Mistral2_binop(var1,constant)
{
#ifdef _DEBUGWRAP
  std::cout << "creating or predicate" << std::endl;
#endif
}

Mistral2_or::~Mistral2_or(){
#ifdef _DEBUGWRAP
  std::cout << "delete or" << std::endl;
#endif
}

Mistral2_Expression* Mistral2_or::add(Mistral2Solver *solver, bool top_level){
  if(!has_been_added()) {
#ifdef _DEBUGWRAP
    std::cout << "add or predicate" << std::endl;
#endif

    _solver = solver;

    bool used = false;
    _vars[0]->add(_solver,false);
    if(_vars[1]) {
      _vars[1]->add(_solver,false);   
      _self = (_vars[0]->_self || _vars[1]->_self);
      used = true;
    } else if(_constant == 0) {
      _self = (_vars[0]->_self == 1);
      used = true;
    } 
    if( top_level && used )
      _solver->solver->add( _self );
      
  }
  return this;
}

Mistral2_eq::Mistral2_eq(Mistral2_Expression *var1, Mistral2_Expression *var2)
  : Mistral2_binop(var1,var2)
{
#ifdef _DEBUGWRAP
  std::cout << "creating equality" << std::endl;
#endif
}

Mistral2_eq::Mistral2_eq(Mistral2_Expression *var1, int constant)
  : Mistral2_binop(var1,constant)
{
#ifdef _DEBUGWRAP
  std::cout << "creating equality" << std::endl;
#endif
}

Mistral2_eq::~Mistral2_eq(){
#ifdef _DEBUGWRAP
  std::cout << "delete eq" << std::endl;
#endif
}

Mistral2_Expression* Mistral2_eq::add(Mistral2Solver *solver, bool top_level){
  if(!has_been_added()) {
#ifdef _DEBUGWRAP
    std::cout << "add equality constraint (" << _constant << ")" << std::endl;
#endif
    
    _solver = solver;
    
    _vars[0]->add(_solver,false);
    if(_vars[1]) {
      _vars[1]->add(_solver,false);   
      _self = (_vars[0]->_self == _vars[1]->_self);
    } else {
      _self = (_vars[0]->_self == _constant);
    }

    if( top_level ) {
      _solver->solver->add( _self );
    }
  }

  return this;
}

/* Disequality operator */

Mistral2_ne::Mistral2_ne(Mistral2_Expression *var1, Mistral2_Expression *var2)
  : Mistral2_binop(var1,var2)
{
#ifdef _DEBUGWRAP
  std::cout << "creating notequal" << std::endl;
#endif
}

Mistral2_ne::Mistral2_ne(Mistral2_Expression *var1, int constant)
  : Mistral2_binop(var1,constant)
{
#ifdef _DEBUGWRAP
  std::cout << "creating notequal" << std::endl;
#endif
}

Mistral2_ne::~Mistral2_ne(){
#ifdef _DEBUGWRAP
  std::cout << "delete notequal" << std::endl;
#endif
}

Mistral2_Expression* Mistral2_ne::add(Mistral2Solver *solver, bool top_level){
  if(!has_been_added()) {
#ifdef _DEBUGWRAP
    std::cout << "add notequal constraint" << std::endl;
#endif

    _solver = solver;		
    _vars[0]->add(_solver,false);
				
    if(_vars[1]) {
			_vars[1]->add(_solver,false);  
      _self = (_vars[0]->_self != _vars[1]->_self);
    } else {
      _self = (_vars[0]->_self != _constant);
    }
    if( top_level )
      _solver->solver->add( _self );
      
  }
  return this;
}

/* Disjunctive constraint operator */

Mistral2_NoOverlap::Mistral2_NoOverlap(Mistral2_Expression *var1, Mistral2_Expression *var2, int constant, int bonstant)
  : Mistral2_binop(var1,var2)
{
#ifdef _DEBUGWRAP
  std::cout << "creating nooverlap" << std::endl;
#endif
  _constant = constant;
  _bonstant = bonstant;
}

Mistral2_NoOverlap::~Mistral2_NoOverlap()
{
#ifdef _DEBUGWRAP
  std::cout << "delete nooverlap" << std::endl;
#endif
}

Mistral2_Expression* Mistral2_NoOverlap::add(Mistral2Solver *solver, bool top_level) {

  if(!has_been_added()) {
#ifdef _DEBUGWRAP
      std::cout << "add nooverlap constraint" << std::endl;
#endif
    _solver = solver;
    
    _vars[0]->add(_solver,false);
    _vars[1]->add(_solver,false);    
    
    _self = Disjunctive(_vars[0]->_self, _vars[1]->_self, _constant, _bonstant);
    
    if( top_level )
      _solver->solver->add( _self );
  }

  return this;  
}


/* Leq operator */

Mistral2_le::Mistral2_le(Mistral2_Expression *var1, Mistral2_Expression *var2)
  : Mistral2_binop(var1,var2)
{ 
#ifdef _DEBUGWRAP
  std::cout << "Creating less than or equal" << std::endl;
#endif
}

Mistral2_le::Mistral2_le(Mistral2_Expression *var1, int constant)
  : Mistral2_binop(var1,constant)
{
#ifdef _DEBUGWRAP
  std::cout << "Creating less than or equal" << std::endl;
#endif
}


Mistral2_le::~Mistral2_le(){
#ifdef _DEBUGWRAP
  std::cout << "delete lessequal" << std::endl;
#endif
}

Mistral2_Expression* Mistral2_le::add(Mistral2Solver *solver, bool top_level) {
  if(!has_been_added()) {
#ifdef _DEBUGWRAP
    std::cout << "add lessequal constraint" << std::endl;
#endif

    _solver = solver;
    
    _vars[0]->add(_solver,false);
    if(_vars[1]) {
      _vars[1]->add(_solver,false);   
      _self = (_vars[0]->_self <= _vars[1]->_self);
    } else {
      _self = (_vars[0]->_self <= _constant);
    }
    if( top_level )
      _solver->solver->add( _self );
  
  }
  return this;
}

/* Geq operator */

Mistral2_ge::Mistral2_ge(Mistral2_Expression *var1, Mistral2_Expression *var2)
  : Mistral2_binop(var1,var2)
{
#ifdef _DEBUGWRAP
  std::cout << "Creating a greater or equal constraint" << std::endl;
#endif
}

Mistral2_ge::Mistral2_ge(Mistral2_Expression *var1, int constant)
  : Mistral2_binop(var1,constant)
{
#ifdef _DEBUGWRAP
  std::cout << "Creating a greater or equal constraint" << std::endl;
#endif
}

Mistral2_ge::~Mistral2_ge(){
#ifdef _DEBUGWRAP
  std::cout << "delete greaterequal" << std::endl;
#endif
}

Mistral2_Expression* Mistral2_ge::add(Mistral2Solver *solver, bool top_level){ 
  if(!has_been_added()) {
#ifdef _DEBUGWRAP
    std::cout << "add greaterequal constraint" << std::endl;
#endif

    _solver = solver;
    
    _vars[0]->add(_solver,false);
    if(_vars[1]) {
      _vars[1]->add(_solver,false);   
      _self = (_vars[0]->_self >= _vars[1]->_self);
    } else {
      _self = (_vars[0]->_self >= _constant);
    }
    if( top_level )
      _solver->solver->add( _self );
      
  }
 
  return this;
}

/* Lt object */

Mistral2_lt::Mistral2_lt(Mistral2_Expression *var1, Mistral2_Expression *var2)
  : Mistral2_binop(var1,var2)
{
#ifdef _DEBUGWRAP
  std::cout << "Creating a less than constraint" << std::endl;
#endif
}

Mistral2_lt::Mistral2_lt(Mistral2_Expression *var1, int constant)
  : Mistral2_binop(var1,constant)
{
#ifdef _DEBUGWRAP
  std::cout << "Creating a less than constraint" << std::endl;
#endif
}

Mistral2_lt::~Mistral2_lt()
{
#ifdef _DEBUGWRAP
  std::cout << "delete lessthan" << std::endl;
#endif
}

Mistral2_Expression* Mistral2_lt::add(Mistral2Solver *solver, bool top_level) {
  if(!has_been_added()) {
#ifdef _DEBUGWRAP
    std::cout << "add lessthan constraint" << std::endl;
#endif

    _solver = solver;
    
    _vars[0]->add(_solver,false);
    if(_vars[1]) {
      _vars[1]->add(_solver,false);   
      _self = (_vars[0]->_self < _vars[1]->_self);
    } else {
      _self = (_vars[0]->_self < _constant);
    }
    if( top_level )
      _solver->solver->add( _self );
      
  }
  return this;  
}

/* Gt object */

Mistral2_gt::Mistral2_gt(Mistral2_Expression *var1, Mistral2_Expression *var2)
  : Mistral2_binop(var1,var2)
{
#ifdef _DEBUGWRAP
  std::cout << "Creating a greater than constraint" << std::endl;
#endif
}

Mistral2_gt::Mistral2_gt(Mistral2_Expression *var1, int constant)
  : Mistral2_binop(var1,constant)
{
#ifdef _DEBUGWRAP
  std::cout << "Creating a greater than constraint" << std::endl;
#endif
}

Mistral2_gt::~Mistral2_gt()
{
#ifdef _DEBUGWRAP
  std::cout << "delete greaterthan" << std::endl;
#endif
}

Mistral2_Expression* Mistral2_gt::add(Mistral2Solver *solver, bool top_level) {
  if(!has_been_added()) {
#ifdef _DEBUGWRAP
    std::cout << "add greaterthan constraint" << std::endl;
#endif

    _solver = solver;
    
    _vars[0]->add(_solver,false);
    if(_vars[1]) {
      _vars[1]->add(_solver,false);   
      _self = (_vars[0]->_self > _vars[1]->_self);
    } else {
      _self = (_vars[0]->_self > _constant);
    }
    if( top_level )
      _solver->solver->add( _self );
            
  }
  return this;
}


/* Minimise object */

Mistral2_Minimise::Mistral2_Minimise(Mistral2_Expression *var)
  : Mistral2_Expression()
{
#ifdef _DEBUGWRAP
  std::cout << "Creating a minimise objective" << std::endl;
#endif
  _exp = var;
}

Mistral2_Minimise::~Mistral2_Minimise()
{
#ifdef _DEBUGWRAP
  std::cout << "delete minimise" << std::endl;
#endif
}

Mistral2_Expression* Mistral2_Minimise::add(Mistral2Solver *solver, bool top_level) {
  if(!has_been_added()) {
#ifdef _DEBUGWRAP
    std::cout << "add minimise objective" << std::endl;
#endif
      
    _solver = solver;

    //std::cout << 11 << std::endl;

    // This will be the objective function
    _exp = _exp->add(solver, false);

    //std::cout << 22 << std::endl;

    //std::cout << _solver->solver << std::endl;
      
    if(top_level){
      _solver->solver->minimize(_exp->_self);
      //_search_goal = new Mistral::Goal(Mistral::Goal::MINIMIZATION, _exp);
    }

    //std::cout << 33 << std::endl;
    
  }
  return this;
}

/* Maximise object */

Mistral2_Maximise::Mistral2_Maximise(Mistral2_Expression *var)
  : Mistral2_Expression()
{
#ifdef _DEBUGWRAP
  std::cout << "Creating a maximise objective" << std::endl;
#endif

  _exp = var;
}

Mistral2_Maximise::~Mistral2_Maximise()
{
#ifdef _DEBUGWRAP
  std::cout << "delete maximise" << std::endl;
#endif
}

Mistral2_Expression* Mistral2_Maximise::add(Mistral2Solver *solver, bool top_level) {
  if(!has_been_added()) {
#ifdef _DEBUGWRAP
    std::cout << "add maximise objective" << std::endl;
#endif

    _solver = solver;

    // This will be the objective function
    _exp = _exp->add(solver, false);
      
    if(top_level){
      _solver->solver->maximize(_exp->_self);
      //_search_goal = new Mistral::Goal(Mistral::Goal::MAXIMIZATION, _exp);
    }
      
  }
  return this;
}

/**************************************************************
 ********************     Solver        ***********************
 **************************************************************/

Mistral2Solver::Mistral2Solver()
{
#ifdef _DEBUGWRAP
  std::cout << "c (wrapper) creating solver" << std::endl;
#endif

  solver = new Mistral::Solver();
  solver->parameters.verbosity = 2;
  solver->consolidate();

  //_search_goal = NULL;
  
  _restart_policy_str = "geom";
  _heuristic_randomization = 2;
  _var_heuristic_str = "dom/wdeg";
  _val_heuristic_str = "minval+guided";
}

Mistral2Solver::~Mistral2Solver()
{
#ifdef _DEBUGWRAP
  std::cout << "c (wrapper) delete solver" << std::endl;
#endif

  delete solver;
}

void Mistral2Solver::add(Mistral2_Expression* arg)
{
#ifdef _DEBUGWRAP
  std::cout << "adding expression" << std::endl;
#endif
  
  arg->add(this, true);
  
}

void Mistral2Solver::initialise(Mistral2ExpArray& arg)
{
#ifdef _DEBUGWRAP
  std::cout << "Initialising solver with array of expressions" << std::endl;
#endif
}

void Mistral2Solver::initialise()
{
#ifdef _DEBUGWRAP
  std::cout << "Initialising solver with no expressions" << std::endl;
#endif
}

int Mistral2Solver::num_vars() {
#ifdef _DEBUGWRAP
  std::cout << "return number of variables" << std::endl;
#endif
  return 0;
}

int Mistral2Solver::get_degree(int i) {
#ifdef _DEBUGWRAP
  std::cout << "return degree of expression i" << std::endl;
#endif
  return 0;
}

int Mistral2Solver::solveAndRestart(const int policy, 
				    const unsigned int base, 
				    const double factor,
				    const double decay,
				    const int reinit)
{
#ifdef _DEBUGWRAP
  std::cout << " SOLVE!! at level: ..." << std::endl;
#endif

  return solve();
  //return 0;
}

int Mistral2Solver::solve()
{

#ifdef _DEBUGWRAP
  std::cout << " SOLVE!! " << std::endl;
#endif

  solver->consolidate();

  
#ifdef _DEBUGWRAP
  std::cout << solver << std::endl;
  solver->parameters.verbosity = 2; 
#endif

  _branching_heuristic = solver->heuristic_factory(_var_heuristic_str, _val_heuristic_str, _heuristic_randomization);

  //if(!_search_goal)
  //_search_goal = new Mistral::Goal(Mistral::Goal::SATISFACTION);
  _restart_policy = solver->restart_factory(_restart_policy_str); 

  //Mistral::Outcome result = 
  solver->depth_first_search(solver->variables, _branching_heuristic, _restart_policy, NULL, false); //, _search_goal);

  return (is_sat());
}

int Mistral2Solver::startNewSearch()
{
#ifdef _DEBUGWRAP
  std::cout << "starting new search" << std::endl;
#endif
	
  solver->consolidate();

  
#ifdef _DEBUGWRAP
  std::cout << solver << std::endl;
  solver->parameters.verbosity = 2; 
#endif

  _branching_heuristic = solver->heuristic_factory(_var_heuristic_str, _val_heuristic_str, _heuristic_randomization);

  _restart_policy = solver->restart_factory(_restart_policy_str); 

  solver->initialise_search(solver->variables, _branching_heuristic, _restart_policy, NULL);

  solver->statistics.start_time = Mistral::get_run_time();

	
  return 0;
}

int Mistral2Solver::getNextSolution()
{
#ifdef _DEBUGWRAP
  std::cout << "getting next solution" << std::endl;
#endif
	
  return solver->get_next_solution();
}

int Mistral2Solver::sacPreprocess(const int type)
{
#ifdef _DEBUGWRAP
  std::cout << "running sac preprocessing" << std::endl;
#endif
  return 0;
}

bool Mistral2Solver::propagate()
{
#ifdef _DEBUGWRAP
  std::cout << "propagating" << std::endl;
#endif
  return 0;
}

int Mistral2Solver::get_level()
{
#ifdef _DEBUGWRAP
  std::cout << "return current level of solver" << std::endl;
#endif
  return 0;
}

bool Mistral2Solver::undo(const int nlevel) {
#ifdef _DEBUGWRAP
  std::cout << "Going back nlevel levels" << std::endl;
#endif
  return 0;
} 

void Mistral2Solver::save() 
{
#ifdef _DEBUGWRAP
  std::cout << "saving state" << std::endl;
#endif
}

int Mistral2Solver::next(Mistral2_Expression* x, int v)
{
#ifdef _DEBUGWRAP
  std::cout << "get next" << std::endl;
#endif
  return 0;
}

void Mistral2Solver::post(const char* op, Mistral2_Expression* x, int v)
{
#ifdef _DEBUGWRAP
  std::cout << "Posting a constraint" << std::endl;
#endif


}

void Mistral2Solver::deduce(const char* op, Mistral2_Expression* x, int v) {
#ifdef _DEBUGWRAP
  std::cout << "deducing" << std::endl;
#endif
}

void Mistral2Solver::deduce() {
#ifdef _DEBUGWRAP
  std::cout << "deducing" << std::endl;
#endif
}

bool Mistral2Solver::branch_right() {
#ifdef _DEBUGWRAP
  std::cout << "branching right" << std::endl;
#endif
  return 0;
}

void Mistral2Solver::store_solution() {
#ifdef _DEBUGWRAP
  std::cout << "storing solution" << std::endl;
#endif
}

void Mistral2Solver::setHeuristic(const char* var_heuristic, const char* val_heuristic, const int rand)
{
#ifdef _DEBUGWRAP
  std::cout << "Setting heuristics" << std::endl;
#endif

  _heuristic_randomization = rand;
  _var_heuristic_str = var_heuristic;
  _val_heuristic_str = val_heuristic;

}

void Mistral2Solver::addNogood(Mistral2ExpArray& vars, 
			       Mistral2IntArray& vals)
{
#ifdef _DEBUGWRAP
  std::cout << "Adding nogood constraint" <<std::endl;
#endif
}

void Mistral2Solver::guide(Mistral2ExpArray& vars, 
			   Mistral2IntArray& vals,
			   Mistral2DoubleArray& probs)
{
#ifdef _DEBUGWRAP
  std::cout << "Adding nogood constraint" <<std::endl;
#endif
}

void Mistral2Solver::forceFiniteDomain(Mistral2ExpArray& vars)
{
#ifdef _DEBUGWRAP
  std::cout << "Forcing finite domain" <<std::endl;
#endif
}

void Mistral2Solver::backtrackTo(const int level) {
#ifdef _DEBUGWRAP
  std::cout << "backtracking to level" <<std::endl;
#endif
}

void Mistral2Solver::presolve()
{
#ifdef _DEBUGWRAP
  std::cout << "presolving" <<std::endl;
#endif
}

void Mistral2Solver::increase_init_level(const int i)
{
#ifdef _DEBUGWRAP
  std::cout << "increasing initial level" <<std::endl;
#endif
}

void Mistral2Solver::decrease_init_level(const int i)
{
#ifdef _DEBUGWRAP
  std::cout << "decreasing initial level" <<std::endl;
#endif
}

void Mistral2Solver::assign(Mistral2_Expression *X, const int v) {
#ifdef _DEBUGWRAP
  std::cout << "Setting domain of expression X to v" <<std::endl;
#endif
}

void Mistral2Solver::upOneLevel() {
#ifdef _DEBUGWRAP
  std::cout << "stepping up one level" <<std::endl;
#endif
}

void Mistral2Solver::reset(bool full) {
#ifdef _DEBUGWRAP
  std::cout << "resetting solver" <<std::endl;
#endif

  solver->restore(solver->search_root);
}

void Mistral2Solver::setLowerBounds(Mistral2ExpArray& vars, 
				    Mistral2IntArray& vals) {
#ifdef _DEBUGWRAP
  std::cout << "stepping up one level" <<std::endl;
#endif
}

void Mistral2Solver::setUpperBounds(Mistral2ExpArray& vars, 
				    Mistral2IntArray& vals) {
#ifdef _DEBUGWRAP
  std::cout << "stetting upper bounds" <<std::endl;
#endif
}

void Mistral2Solver::setRestartNogood() 
{
#ifdef _DEBUGWRAP
  std::cout << "setting restart no good" <<std::endl;
#endif
}

void Mistral2Solver::setFailureLimit(const int cutoff)
{
#ifdef _DEBUGWRAP
  std::cout << "setting failure limit" <<std::endl;
#endif

  solver->parameters.fail_limit = cutoff;

}

void Mistral2Solver::setNodeLimit(const int cutoff)
{
#ifdef _DEBUGWRAP
  std::cout << "setting node limit" <<std::endl;
#endif

  solver->parameters.node_limit = cutoff;

}

void Mistral2Solver::setTimeLimit(const int cutoff)
{
#ifdef _DEBUGWRAP
  std::cout << "setting time limit" <<std::endl;
#endif

  solver->parameters.time_limit = cutoff;

}

void Mistral2Solver::setVerbosity(const int degree)
{
#ifdef _DEBUGWRAP
  std::cout << "setting verbosity to " << degree <<std::endl;
#endif

  solver->parameters.verbosity = degree;

}

void Mistral2Solver::setRandomized(const int degree)
{
#ifdef _DEBUGWRAP
  std::cout << "setting randomised" <<std::endl;
#endif

  solver->parameters.randomization = degree;

}

void Mistral2Solver::setRandomSeed(const int seed)
{
#ifdef _DEBUGWRAP
  std::cout << "setting random seed" <<std::endl;
#endif

  solver->parameters.seed = seed;

}

bool Mistral2Solver::is_opt()
{
#ifdef _DEBUGWRAP
  std::cout << "returning is_opt()" <<std::endl;
#endif
  return solver->statistics.outcome == 3;
}

bool Mistral2Solver::is_sat()
{
#ifdef _DEBUGWRAP
  std::cout << "returning is satisfied?" <<std::endl;
#endif
  return solver->statistics.num_solutions > 0; //solver->statistics.outcome == Mistral::SAT || solver->statistics.outcome == Mistral::OPT;
}

bool Mistral2Solver::is_unsat()
{
#ifdef _DEBUGWRAP
  std::cout << "returning is NOT satisfied?" <<std::endl;
#endif
  return solver->statistics.outcome == 0;
}

void Mistral2Solver::printStatistics()
{
#ifdef _DEBUGWRAP
  std::cout << "printing statistics" <<std::endl;
#endif

  solver->statistics.print_full(std::cout);

}

int Mistral2Solver::getBacktracks()
{
#ifdef _DEBUGWRAP
  std::cout << "return number of backtracks" <<std::endl;
#endif
  return solver->statistics.num_backtracks;
}

int Mistral2Solver::getNodes()
{
#ifdef _DEBUGWRAP
  std::cout << "return number of nodes" <<std::endl;
#endif
  return solver->statistics.num_nodes;
}

int Mistral2Solver::getFailures()
{
#ifdef _DEBUGWRAP
  std::cout << "return number of failures" <<std::endl;
#endif
  return solver->statistics.num_failures;
}

int Mistral2Solver::getChecks()
{
#ifdef _DEBUGWRAP
  std::cout << "return number of checks" <<std::endl;
#endif
  return solver->statistics.num_propagations;
}

int Mistral2Solver::getPropags()
{
#ifdef _DEBUGWRAP
  std::cout << "return amount of propagation" <<std::endl;
#endif
  return solver->statistics.num_propagations;
}

double Mistral2Solver::getTime()
{
#ifdef _DEBUGWRAP
  std::cout << "return duration" <<std::endl;
#endif
  return solver->statistics.end_time - solver->statistics.start_time;
}

int Mistral2Solver::getNumVariables()
{
#ifdef _DEBUGWRAP
  std::cout << "return number of variables" <<std::endl;
#endif
  return solver->variables.size;
}

int Mistral2Solver::getNumConstraints()
{
#ifdef _DEBUGWRAP
  std::cout << "return number of constraints" <<std::endl;
#endif 
  return solver->posted_constraints.size;
}

int Mistral2Solver::getRandomNumber()
{
  return Mistral::randint(0xffffffff);
}

void Mistral2Solver::printPython()
{
	solver->consolidate();
	solver->display(std::cout);
}

