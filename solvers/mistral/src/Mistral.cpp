
#include "Mistral.hpp"

#include "../models/src/xcsp/XMLParser_libxml2.hh"
#include "../models/src/xcsp/MistralCallback.hh"
////#include "../models/src/xcsp/CSPParserCallback.hh"

#include <mistral_sol.h>

using namespace CSPXMLParser;
using namespace Mistral;

MistralCallback *cb; 
XMLParser_libxml2<MistralCallback> *parser;

const int num_features = 36;

/**************************************************************
 ********************     EXPRESSION        *******************
 **************************************************************/

void Mistral_Expression::initialise()
{ 
  _self = NULL;
  _solver = NULL;
}

void Mistral_Expression::print_python() const 
{
  _self->print_python();
}

Mistral_Expression::Mistral_Expression(BuildObject *x)
{

  initialise();
  nbj_ident = -1;
  _self = x;

#ifdef _DEBUGWRAP
  std::cout << "creating a variable [" << 0 << ".." << (nval-1) << "]" << std::endl;
#endif

}

Mistral_Expression::Mistral_Expression()
{
  initialise();

#ifdef _DEBUGWRAP
  std::cout << "creating an empty expression" << std::endl;
#endif

}

Mistral_Expression::Mistral_Expression(const int nval)
{
  initialise();
  _self = CSP::_Variable(nval);

#ifdef _DEBUGWRAP
  std::cout << "creating a variable [" << 0 << ".." << (nval-1) << "]" << std::endl;
#endif

}

Mistral_Expression::Mistral_Expression(const int lb, const int ub)
{
  initialise();
  _self = CSP::_Variable(lb, ub);

#ifdef _DEBUGWRAP
  std::cout << "creating a variable [" << lb << ".." << ub << "]" << std::endl;
#endif

}

Mistral_Expression::Mistral_Expression(MistralIntArray& vals)
{
  initialise();

  int _size   = vals.size();
  int *_values = new int[_size];
  for(int i=0; i<_size; ++i)
    _values[i] = vals.get_item(i);
  _self = CSP::_Variable(_values, _size);

  delete[] _values;

#ifdef _DEBUGWRAP
  std::cout << "creating a variable [" << _self->min() << ".." << _self->max() << "]" << std::endl;
#endif

}

const char* Mistral_Expression::get_type() const
{
  return _self->get_type();
}

int Mistral_Expression::get_arity() const
{
  return _self->get_arity();
}
 
Mistral_Expression* Mistral_Expression::get_child(const int i)
{
  return new Mistral_Expression(((BuildObjectPredicate*)_self)->scope[i]);
}

int Mistral_Expression::get_id() const
{
  return _self->getVarId();
}

int Mistral_Expression::getVariableId() const
{
  return _self->getVariableId();
}

int Mistral_Expression::next(int v)
{
  int nxt = v;

  if(_self->varptr_) 
    nxt = _self->varptr_->getNext(v);
  if(nxt == NOVAL) nxt = v;

  return nxt;
}

int Mistral_Expression::get_value() const
{
  return _self->getValue();
}

int Mistral_Expression::get_size() const
{
  return _self->getSize();
}

int Mistral_Expression::get_max() const
{
  return _self->getMax();
}

int Mistral_Expression::get_min() const
{
  return _self->getMin();
}

bool Mistral_Expression::contain(const int v) const
{
  return _self->isIn(v);
}


Mistral_Expression::~Mistral_Expression()
{

#ifdef _DEBUGWRAP
  std::cout << "delete expression" << std::endl;
#endif

  //_self = NULL;
  //_solver = NULL;
}


bool Mistral_Expression::has_been_added() const
{
  return (_solver != NULL);
}

Mistral_Expression* Mistral_Expression::add(MistralSolver *solver, bool top_level){

  if(!has_been_added()) {
#ifdef _DEBUGWRAP
    std::cout << "add variable [" << _self->min() << ".." << _self->max() << "]" << std::endl;
#endif
    _solver = solver;
    if(top_level)
      solver->model->add( _self );    
  }

  return this;
}

// /* Binary operators */

Mistral_binop::Mistral_binop(Mistral_Expression *var1, Mistral_Expression *var2)
  : Mistral_Expression()
{
#ifdef _DEBUGWRAP
  std::cout << "creating a binary operator" << std::endl;
#endif

  _vars[0] = var1;
  _vars[1] = var2;

}

Mistral_binop::Mistral_binop(Mistral_Expression *var1, int constant)
  : Mistral_Expression()
{
#ifdef _DEBUGWRAP
  std::cout << "creating a binary (constant) operator" << std::endl;
#endif

  _vars[0] = var1;
  _vars[1] = NULL;
  _constant = constant;

}


Mistral_binop::~Mistral_binop(){

#ifdef _DEBUGWRAP
  std::cout << "delete binary operator" << std::endl;
#endif

}


Mistral_Min::Mistral_Min( MistralExpArray& vars ) 
  : Mistral_Expression() 
{

  _vars = vars;

#ifdef _DEBUGWRAP
  std::cout << "creating an alldiff constraint" << std::endl;
#endif

}

Mistral_Min::Mistral_Min( Mistral_Expression *var1, Mistral_Expression *var2 ) 
  : Mistral_Expression() 
{
  _vars.add(var1);
  _vars.add(var2); 

#ifdef _DEBUGWRAP
  std::cout << "creating a binary alldiff constraint" << std::endl;
#endif

}

Mistral_Min::~Mistral_Min()
{

#ifdef _DEBUGWRAP
  std::cout << "delete alldiff" << std::endl;
#endif

}

Mistral_Expression* Mistral_Min::add(MistralSolver *solver, bool top_level) {

  if(!has_been_added()) {

#ifdef _DEBUGWRAP
    std::cout << "add an alldiff constraint" << std::endl;
#endif

    _solver = solver;
    
    int i, n=_vars.size();  
    for(i=0; i<n; ++i) 
      _vars.get_item(i)->add(_solver,false);

    if(n == 2) {
      _self = CSP::_Min(_vars.get_item(0)->_self,
			_vars.get_item(1)->_self);
    } else if(n > 2) {
      BuildObject **scope = new BuildObject*[n+1];
      for(i=0; i<n; ++i) scope[i] = _vars.get_item(i)->_self;
      _self = CSP::_Min(scope, n);
    }

    if( top_level )
      _solver->model->add( _self );
    
  }

  return this;
}


Mistral_Max::Mistral_Max( MistralExpArray& vars ) 
  : Mistral_Expression() 
{

  _vars = vars;

#ifdef _DEBUGWRAP
  std::cout << "creating an alldiff constraint" << std::endl;
#endif

}

Mistral_Max::Mistral_Max( Mistral_Expression *var1, Mistral_Expression *var2 ) 
  : Mistral_Expression() 
{
  _vars.add(var1);
  _vars.add(var2); 

#ifdef _DEBUGWRAP
  std::cout << "creating a binary alldiff constraint" << std::endl;
#endif

}

Mistral_Max::~Mistral_Max()
{

#ifdef _DEBUGWRAP
  std::cout << "delete alldiff" << std::endl;
#endif

}

Mistral_Expression* Mistral_Max::add(MistralSolver *solver, bool top_level) {

  if(!has_been_added()) {

#ifdef _DEBUGWRAP
    std::cout << "add an alldiff constraint" << std::endl;
#endif

    _solver = solver;
    
    int i, n=_vars.size();  
    for(i=0; i<n; ++i) 
      _vars.get_item(i)->add(_solver,false);

    if(n == 2) {
      _self = CSP::_Max(_vars.get_item(0)->_self,
			_vars.get_item(1)->_self);
    } else if(n > 2) {
      BuildObject **scope = new BuildObject*[n+1];
      for(i=0; i<n; ++i) scope[i] = _vars.get_item(i)->_self;
      _self = CSP::_Max(scope, n);
    }

    if( top_level )
      _solver->model->add( _self );
    
  }

  return this;
}


Mistral_Table::Mistral_Table( MistralExpArray& vars, MistralIntArray& tuples, const char* type ) 
  : Mistral_Expression() 
{
  _vars = vars;
  _tuples = tuples;

#ifdef _DEBUGWRAP
  std::cout << "creating a table constraint" << std::endl;
#endif
  
  spin = (!strcmp(type,"support"));

}

Mistral_Table::Mistral_Table( Mistral_Expression *var1, Mistral_Expression *var2, 
			      MistralIntArray& tuples, const char* type ) 
  : Mistral_Expression() 
{
  _vars.add(var1);
  _vars.add(var2); 
  _tuples = tuples;

#ifdef _DEBUGWRAP
  std::cout << "creating a binary table constraint" << std::endl;
#endif

  spin = (!strcmp(type,"support"));

}

Mistral_Table::~Mistral_Table()
{

#ifdef _DEBUGWRAP
  std::cout << "delete alldiff" << std::endl;
#endif

}

Mistral_Expression* Mistral_Table::add(MistralSolver *solver, bool top_level) {

  if(!has_been_added()) {

#ifdef _DEBUGWRAP
    std::cout << "add a table constraint" << std::endl;
#endif

    _solver = solver;
    
    int i, j, n=_vars.size(), m=_tuples.size();  
    for(i=0; i<n; ++i) 
      _vars.get_item(i)->add(_solver,false);

    BuildObject **scope = new BuildObject*[n+1];
    for(i=0; i<n; ++i) scope[i] = _vars.get_item(i)->_self;
    _self = CSP::_Table(scope, n, spin);

    BuildObjectTable* tabptr_ = (BuildObjectTable*)(((BuildObjectPredicate*)_self)->relation);
						    
    int tuple[n];
    for(int i=0; i<m; ++i) {
      j = (i%n);
      tuple[j] = _tuples.get_item(i);
      if(j == n-1) tabptr_->add(tuple);
    }
    
    if( top_level )
      _solver->model->add( _self );
  }

  return this;
}

void Mistral_Table::add(MistralIntArray& tuple) {
  for(int i=0; i<tuple.size(); ++i)
    _tuples.add(tuple.get_item(i));
}


Mistral_AllDiff::Mistral_AllDiff( MistralExpArray& vars ) 
  : Mistral_Expression() 
{

  _vars = vars;

#ifdef _DEBUGWRAP
  std::cout << "creating an alldiff constraint" << std::endl;
#endif

}

Mistral_AllDiff::Mistral_AllDiff( Mistral_Expression *var1, Mistral_Expression *var2 ) 
  : Mistral_Expression() 
{
  _vars.add(var1);
  _vars.add(var2); 

#ifdef _DEBUGWRAP
  std::cout << "creating a binary alldiff constraint" << std::endl;
#endif

}

Mistral_AllDiff::~Mistral_AllDiff()
{

#ifdef _DEBUGWRAP
  std::cout << "delete alldiff" << std::endl;
#endif

}

Mistral_Expression* Mistral_AllDiff::add(MistralSolver *solver, bool top_level) {

  if(!has_been_added()) {

#ifdef _DEBUGWRAP
    std::cout << "add an alldiff constraint" << std::endl;
#endif

    _solver = solver;
    
    int i, n=_vars.size();  
    for(i=0; i<n; ++i) 
      _vars.get_item(i)->add(_solver,false);

    if(n == 2) {
      _self = CSP::_Equal(_vars.get_item(0)->_self,
			  _vars.get_item(1)->_self, 0);      
    } else if(n > 2) {
      BuildObject **scope = new BuildObject*[n+1];
      for(i=0; i<n; ++i) scope[i] = _vars.get_item(i)->_self;
      _self = CSP::_AllDifferent(scope, n, 3, NOVAL/2);
    }

    if( top_level )
      _solver->model->add( _self );

  }

  return this;
}


Mistral_Gcc::Mistral_Gcc(MistralExpArray& vars,
			 MistralIntArray& vals,
			 MistralIntArray& lb_card,
			 MistralIntArray& ub_card)
  : Mistral_Expression() 
{

  _vars = vars;
  _vals = vals;
  _lb_card = lb_card;
  _ub_card = ub_card;

#ifdef _DEBUGWRAP
  std::cout << "creating a gcc constraint" << std::endl;
#endif

}

Mistral_Gcc::~Mistral_Gcc()
{

#ifdef _DEBUGWRAP
  std::cout << "delete gcc" << std::endl;
#endif

}

Mistral_Expression* Mistral_Gcc::add(MistralSolver *solver, bool top_level) {

  if(!has_been_added()) {
#ifdef _DEBUGWRAP
      std::cout << "add gcc constraint" << std::endl;
#endif

    _solver = solver;
    
    int i, n=_vars.size(), m=_vals.size();
    int min_val=_vals.get_item(0);
    int max_val=_vals.get_item(m-1);
    int M = (max_val - min_val + 1);
    BuildObject **scope = new BuildObject*[n+1];

    int *tmp_lb = new int[M];
    int *tmp_ub = new int[M];

    for(i=0; i<n; ++i) {
      _vars.get_item(i)->add(_solver,false);
      scope[i] = _vars.get_item(i)->_self;
    }
    for(i=0; i<M; ++i) {
      tmp_lb[i] = _lb_card.get_item(i);
      tmp_ub[i] = _ub_card.get_item(i);
    }

    _self = CSP::_Gcc(scope, n, min_val, max_val, tmp_lb, tmp_ub);

    delete [] tmp_lb;
    delete [] tmp_ub;

    if( top_level )
      _solver->model->add( _self );

  }

  return this;
}

Mistral_Element::Mistral_Element( MistralExpArray& vars ) 
  : Mistral_Expression() 
{

  _vars = vars;

#ifdef _DEBUGWRAP
  std::cout << "creating element" << std::endl;
#endif

}

Mistral_Element::~Mistral_Element()
{
#ifdef _DEBUGWRAP
  std::cout << "delete element" << std::endl;
#endif
}

Mistral_Expression* Mistral_Element::add(MistralSolver *solver, bool top_level) {

  if(!has_been_added()) {
#ifdef _DEBUGWRAP
      std::cout << "add element constraint" << std::endl;
#endif
    _solver = solver;
    
    int i, n=_vars.size();  
    BuildObject **scope = new BuildObject*[n+1];
      
    for(i=0; i<n; ++i) {
      _vars.get_item(i)->add(_solver,false);
      scope[i] = _vars.get_item(i)->_self;
    }

    _self = CSP::_Element(scope, n, 0);
  }

  return this;
}


Mistral_LeqLex::Mistral_LeqLex( MistralExpArray& vars ) 
  : Mistral_Expression() 
{

  _vars = vars;

#ifdef _DEBUGWRAP
  std::cout << "creating lexleq" << std::endl;
#endif

}

Mistral_LeqLex::~Mistral_LeqLex()
{
#ifdef _DEBUGWRAP
  std::cout << "delete leqlex" << std::endl;
#endif
}

Mistral_Expression* Mistral_LeqLex::add(MistralSolver *solver, bool top_level) {

  if(!has_been_added()) {
#ifdef _DEBUGWRAP
      std::cout << "add leqlex constraint" << std::endl;
#endif
    _solver = solver;
    
    int i, n=_vars.size();  
    BuildObject **scope = new BuildObject*[n+1];
    
    for(i=0; i<n; ++i) {
      _vars.get_item(i)->add(_solver,false);
      scope[i] = _vars.get_item(i)->_self;
    }

    _self = CSP::_LexOrder(scope, n, 1);

    if( top_level )
      _solver->model->add( _self );
  }
  
  return this;
}


Mistral_LessLex::Mistral_LessLex( MistralExpArray& vars ) 
  : Mistral_Expression() 
{

  _vars = vars;

#ifdef _DEBUGWRAP
  std::cout << "creating lexless" << std::endl;
#endif

}

Mistral_LessLex::~Mistral_LessLex()
{
#ifdef _DEBUGWRAP
  std::cout << "delete leslex" << std::endl;
#endif
}

Mistral_Expression* Mistral_LessLex::add(MistralSolver *solver, bool top_level) {

  if(!has_been_added()) {
#ifdef _DEBUGWRAP
      std::cout << "add lesslex constraint" << std::endl;
#endif
    _solver = solver;
    
    int i, n=_vars.size();  
    BuildObject **scope = new BuildObject*[n+1];
    
    for(i=0; i<n; ++i) {
      _vars.get_item(i)->add(_solver,false);
      scope[i] = _vars.get_item(i)->_self;
    }

    _self = CSP::_LexOrder(scope, n, 0);

    if( top_level )
      _solver->model->add( _self );

  }

  return this;
}
 

Mistral_Sum::Mistral_Sum(MistralExpArray& vars, 
			 MistralIntArray& weights, 
			 const int offset)
  : Mistral_Expression() 
{
#ifdef _DEBUGWRAP
  std::cout << "creating sum" << std::endl;
#endif
  //_offset = offset;
  _vars = vars;
  _weights = weights;
  _weights.add(offset);
}

Mistral_Sum::Mistral_Sum(Mistral_Expression *arg1, 
			 Mistral_Expression *arg2, 
			 MistralIntArray& w, 
			 const int offset)
  : Mistral_Expression() 
{
#ifdef _DEBUGWRAP
  std::cout << "creating sum" << std::endl;
#endif
  //_offset = offset;
  _vars.add(arg1);
  _vars.add(arg2);
  _weights = w;
  _weights.add(offset);
}

Mistral_Sum::Mistral_Sum(Mistral_Expression *arg, 
			 MistralIntArray& w, 
			 const int offset)
  : Mistral_Expression() 
{
#ifdef _DEBUGWRAP
  std::cout << "creating sum" << std::endl;
#endif
  //_offset = offset;
  _vars.add(arg);
  _weights = w;
  _weights.add(offset);
}

Mistral_Sum::Mistral_Sum()
  : Mistral_Expression()
{
  //_offset = 0;
}

Mistral_Sum::~Mistral_Sum(){
  
#ifdef _DEBUGWRAP
  std::cout << "delete sum" << std::endl;
#endif

}

void Mistral_Sum::addVar(Mistral_Expression* v){
  _vars.add(v);
}

void Mistral_Sum::addWeight(const int w){
  _weights.add(w);
}

Mistral_Expression* Mistral_Sum::add(MistralSolver *solver, bool top_level){

  if(!has_been_added()) {
#ifdef _DEBUGWRAP
      std::cout << "add sum constraint" << std::endl;
#endif
    _solver = solver;
    
    int i, n=_vars.size();  
    BuildObject **scope = new BuildObject*[n+1];
    int *w = new int[n+1];
    for(i=0; i<n; ++i) {
      _vars.get_item(i)->add(_solver,false);
      scope[i] = _vars.get_item(i)->_self;
      w[i] = _weights.get_item(i);
    }
    w[n] = _weights.get_item(n);
    
    _self = CSP::_Sum(scope, n, w);

    if( top_level )
      _solver->model->add( _self );

  }

  return this;
}

Mistral_mul::Mistral_mul(Mistral_Expression *var1, Mistral_Expression *var2)
  : Mistral_binop(var1,var2)
{

#ifdef _DEBUGWRAP
  std::cout << "creating muluality" << std::endl;
#endif

}

Mistral_mul::Mistral_mul(Mistral_Expression *var1, int constant)
  : Mistral_binop(var1,constant)
{

#ifdef _DEBUGWRAP
  std::cout << "creating muluality" << std::endl;
#endif

}


Mistral_mul::~Mistral_mul(){

#ifdef _DEBUGWRAP
  std::cout << "delete mul" << std::endl;
#endif

}

Mistral_Expression* Mistral_mul::add(MistralSolver *solver, bool top_level){

  if(!has_been_added()) {
#ifdef _DEBUGWRAP
      std::cout << "add mul constraint" << std::endl;
#endif
    _solver = solver;
    
    _vars[0]->add(_solver,false);
    if(_vars[1]) {
      _vars[1]->add(_solver,false);   
      _self = CSP::_Mul(_vars[0]->_self,
			_vars[1]->_self);
    } else {
      _self = CSP::_Mul(_vars[0]->_self,
			_constant);
    }
    if( top_level )
      _solver->model->add( _self );
  }
  
  return this;
}

Mistral_div::Mistral_div(Mistral_Expression *var1, Mistral_Expression *var2)
  : Mistral_binop(var1,var2)
{

#ifdef _DEBUGWRAP
  std::cout << "creating divuality" << std::endl;
#endif

}

Mistral_div::Mistral_div(Mistral_Expression *var1, int constant)
  : Mistral_binop(var1,constant)
{

#ifdef _DEBUGWRAP
  std::cout << "creating divuality" << std::endl;
#endif

}


Mistral_div::~Mistral_div(){

#ifdef _DEBUGWRAP
  std::cout << "delete div" << std::endl;
#endif

}

Mistral_Expression* Mistral_div::add(MistralSolver *solver, bool top_level){

  if(!has_been_added()) {
#ifdef _DEBUGWRAP
      std::cout << "add div constraint" << std::endl;
#endif
    _solver = solver;
    
    _vars[0]->add(_solver,false);
    if(_vars[1]) {
      _vars[1]->add(_solver,false);   
      _self = CSP::_Div(_vars[0]->_self,
			_vars[1]->_self);
    } else {
      std::cout << "s Not supported \nc (division by integer) - exiting" << std::endl;
      exit(1);
      //_self = CSP::_Div(_vars[0]->_self,
      //_constant);
    }
    if( top_level )
      _solver->model->add( _self );
    
  }
  
  return this;
}

Mistral_mod::Mistral_mod(Mistral_Expression *var1, Mistral_Expression *var2)
  : Mistral_binop(var1,var2)
{

#ifdef _DEBUGWRAP
  std::cout << "creating moduality" << std::endl;
#endif

}

Mistral_mod::Mistral_mod(Mistral_Expression *var1, int constant)
  : Mistral_binop(var1,constant)
{

#ifdef _DEBUGWRAP
  std::cout << "creating moduality" << std::endl;
#endif

}


Mistral_mod::~Mistral_mod(){

#ifdef _DEBUGWRAP
  std::cout << "delete mod" << std::endl;
#endif

}

Mistral_Expression* Mistral_mod::add(MistralSolver *solver, bool top_level){

  if(!has_been_added()) {
#ifdef _DEBUGWRAP
      std::cout << "add mod constraint" << std::endl;
#endif
    _solver = solver;
    
    _vars[0]->add(_solver,false);
    if(_vars[1]) {
      _vars[1]->add(_solver,false);   
      _self = CSP::_Mod(_vars[0]->_self,
			_vars[1]->_self);
    } else {
      _self = CSP::_Mod(_vars[0]->_self,
			_constant);
    }
    if( top_level )
      _solver->model->add( _self );
   
  }
  
  return this;
}

Mistral_and::Mistral_and(Mistral_Expression *var1, Mistral_Expression *var2)
  : Mistral_binop(var1,var2)
{

#ifdef _DEBUGWRAP
  std::cout << "creating anduality" << std::endl;
#endif

}

Mistral_and::Mistral_and(Mistral_Expression *var1, int constant)
  : Mistral_binop(var1,constant)
{

#ifdef _DEBUGWRAP
  std::cout << "creating anduality" << std::endl;
#endif

}


Mistral_and::~Mistral_and(){

#ifdef _DEBUGWRAP
  std::cout << "delete and" << std::endl;
#endif

}

Mistral_Expression* Mistral_and::add(MistralSolver *solver, bool top_level){

  if(!has_been_added()) {
#ifdef _DEBUGWRAP
      std::cout << "add and constraint" << std::endl;
#endif
    _solver = solver;
    
    _vars[0]->add(_solver,false);
    if(_vars[1]) {
      _vars[1]->add(_solver,false);   
      _self = CSP::_And(_vars[0]->_self,
			_vars[1]->_self);
    } else if(_constant) {
      _self = CSP::_Equal(_vars[0]->_self, 1, 1);
    } else {
      _self = CSP::_Equal(_vars[0]->_self, NOVAL, 1);
      //std::cerr << "inconsistent constraint, exiting" <<std::endl;
    }
    if( top_level )
      _solver->model->add( _self );
  }
  
  return this;
}

Mistral_or::Mistral_or(Mistral_Expression *var1, Mistral_Expression *var2)
  : Mistral_binop(var1,var2)
{

#ifdef _DEBUGWRAP
  std::cout << "creating oruality" << std::endl;
#endif

}

Mistral_or::Mistral_or(Mistral_Expression *var1, int constant)
  : Mistral_binop(var1,constant)
{

#ifdef _DEBUGWRAP
  std::cout << "creating oruality" << std::endl;
#endif

}


Mistral_or::~Mistral_or(){

#ifdef _DEBUGWRAP
  std::cout << "delete or" << std::endl;
#endif

}

Mistral_Expression* Mistral_or::add(MistralSolver *solver, bool top_level){

  if(!has_been_added()) {
#ifdef _DEBUGWRAP
      std::cout << "add or constraint" << std::endl;
#endif
    _solver = solver;
    
    _vars[0]->add(_solver,false);
    if(_vars[1]) {
      _vars[1]->add(_solver,false);   
      _self = CSP::_Or(_vars[0]->_self,
		       _vars[1]->_self);
    } else if(_constant == 0) {
      _self = CSP::_Equal(_vars[0]->_self, 1, 1);
    } else {
      _self = CSP::_Equal(_vars[0]->_self, NOVAL, 0);
    }
    if( top_level )
      _solver->model->add( _self );
  }
  
  return this;
}

Mistral_eq::Mistral_eq(Mistral_Expression *var1, Mistral_Expression *var2)
  : Mistral_binop(var1,var2)
{

#ifdef _DEBUGWRAP
  std::cout << "creating equality" << std::endl;
#endif

}

Mistral_eq::Mistral_eq(Mistral_Expression *var1, int constant)
  : Mistral_binop(var1,constant)
{

#ifdef _DEBUGWRAP
  std::cout << "creating equality" << std::endl;
#endif

}

Mistral_eq::~Mistral_eq(){

#ifdef _DEBUGWRAP
  std::cout << "delete eq" << std::endl;
#endif


}

Mistral_Expression* Mistral_eq::add(MistralSolver *solver, bool top_level){

  if(!has_been_added()) {
#ifdef _DEBUGWRAP
      std::cout << "add equality constraint" << std::endl;
#endif
    _solver = solver;
    
    _vars[0]->add(_solver,false);
    if(_vars[1]) {
      _vars[1]->add(_solver,false);   
      _self = CSP::_Equal(_vars[0]->_self,
			  _vars[1]->_self, 1);
    } else {
      _self = CSP::_Equal(_vars[0]->_self,
			  _constant, 1);
    }
    if( top_level )
      _solver->model->add( _self );
    
  }
  
  return this;
}

/* Disequality operator */

Mistral_ne::Mistral_ne(Mistral_Expression *var1, Mistral_Expression *var2)
  : Mistral_binop(var1,var2)
{

#ifdef _DEBUGWRAP
  std::cout << "creating notequal" << std::endl;
#endif

}

Mistral_ne::Mistral_ne(Mistral_Expression *var1, int constant)
  : Mistral_binop(var1,constant)
{

#ifdef _DEBUGWRAP
  std::cout << "creating notequal" << std::endl;
#endif

}

Mistral_ne::~Mistral_ne(){

#ifdef _DEBUGWRAP
  std::cout << "delete notequal" << std::endl;
#endif

}

Mistral_Expression* Mistral_ne::add(MistralSolver *solver, bool top_level){

  if(!has_been_added()) {
#ifdef _DEBUGWRAP
      std::cout << "add notequal constraint" << std::endl;
#endif
    _solver = solver;
    
    _vars[0]->add(_solver,false);
    if(_vars[1]) {
      _vars[1]->add(_solver,false);    
      _self = CSP::_Equal(_vars[0]->_self,
			  _vars[1]->_self, 0);
    } else {
      _self = CSP::_Equal(_vars[0]->_self, 
			  _constant, 0);
    }
    if( top_level )
      _solver->model->add( _self );
    
  }
  
  return this;
}

Mistral_NoOverlap::Mistral_NoOverlap(Mistral_Expression *var1, Mistral_Expression *var2, int constant, int bonstant)
  : Mistral_binop(var1,var2)
{
#ifdef _DEBUGWRAP
  std::cout << "creating nooverlap" << std::endl;
#endif
  _constant = constant;
  _bonstant = bonstant;
}

Mistral_NoOverlap::~Mistral_NoOverlap()
{
#ifdef _DEBUGWRAP
  std::cout << "delete nooverlap" << std::endl;
#endif
}

Mistral_Expression* Mistral_NoOverlap::add(MistralSolver *solver, bool top_level) {

  if(!has_been_added()) {
#ifdef _DEBUGWRAP
      std::cout << "add nooverlap constraint" << std::endl;
#endif
    _solver = solver;
    
    _vars[0]->add(_solver,false);
    _vars[1]->add(_solver,false);    
    
    _self = CSP::_Disjunctive(_vars[0]->_self, _constant, _vars[1]->_self, _bonstant, 0);
    
    if( top_level )
      _solver->model->add( _self );
  }

  return this;  
}

Mistral_Precedence::Mistral_Precedence(Mistral_Expression *var1, Mistral_Expression *var2, int constant)
  : Mistral_binop(var1,var2)
{
#ifdef _DEBUGWRAP
  std::cout << "creating precedence" << std::endl;
#endif
  _constant = constant;
}

Mistral_Precedence::~Mistral_Precedence()
{
#ifdef _DEBUGWRAP
  std::cout << "delete precedence" << std::endl;
#endif
}

Mistral_Expression* Mistral_Precedence::add(MistralSolver *solver, bool top_level) {

  if(!has_been_added()) {
#ifdef _DEBUGWRAP
      std::cout << "add precedence constraint" << std::endl;
#endif
    _solver = solver;
    
    _vars[0]->add(_solver,false);
    _vars[1]->add(_solver,false);    
    
    _self = CSP::_Precedence(_vars[0]->_self, _constant, _vars[1]->_self);
    
    if( top_level )
      _solver->model->add( _self );

  }

  return this;  
}

/* Leq operator */

Mistral_le::Mistral_le(Mistral_Expression *var1, Mistral_Expression *var2)
  : Mistral_binop(var1,var2)
{ 

#ifdef _DEBUGWRAP
  std::cout << "Creating less than or equal" << std::endl;
#endif

}

Mistral_le::Mistral_le(Mistral_Expression *var1, int constant)
  : Mistral_binop(var1,constant)
{

#ifdef _DEBUGWRAP
  std::cout << "Creating less than or equal" << std::endl;
#endif

}


Mistral_le::~Mistral_le(){

#ifdef _DEBUGWRAP
  std::cout << "delete lessequal" << std::endl;
#endif

}

Mistral_Expression* Mistral_le::add(MistralSolver *solver, bool top_level) {

  if(!has_been_added()) {
#ifdef _DEBUGWRAP
      std::cout << "add lessequal constraint" << std::endl;
#endif
    _solver = solver;
    
    _vars[0]->add(_solver,false);
    if(_vars[1])  {
      _vars[1]->add(_solver,false);    
      _self = CSP::_Precedence(_vars[0]->_self, 0, 
			       _vars[1]->_self);
    } else {
      _self = CSP::_Precedence(_vars[0]->_self, _constant); 
    }
    if( top_level )
      _solver->model->add( _self );

  }
  
  return this;
}

/* Geq operator */

Mistral_ge::Mistral_ge(Mistral_Expression *var1, Mistral_Expression *var2)
  : Mistral_binop(var1,var2)
{
#ifdef _DEBUGWRAP
  std::cout << "Creating a greater or equal constraint" << std::endl;
#endif
}

Mistral_ge::Mistral_ge(Mistral_Expression *var1, int constant)
  : Mistral_binop(var1,constant)
{

#ifdef _DEBUGWRAP
  std::cout << "Creating a greater or equal constraint" << std::endl;
#endif

}

Mistral_ge::~Mistral_ge(){

#ifdef _DEBUGWRAP
  std::cout << "delete greaterequal" << std::endl;
#endif

}

Mistral_Expression* Mistral_ge::add(MistralSolver *solver, bool top_level){ 

  if(!has_been_added()) {
#ifdef _DEBUGWRAP
      std::cout << "add greaterequal constraint" << std::endl;
#endif
    _solver = solver;
    
    _vars[0]->add(_solver,false);
    if(_vars[1])  {
      _vars[1]->add(_solver,false);    
      _self = CSP::_Precedence(_vars[1]->_self, 0, 
			       _vars[0]->_self);
    } else {
      _self = CSP::_Precedence(_constant, _vars[0]->_self); 
    }
    if( top_level )
      _solver->model->add( _self );
    
  }
  
  return this;
}

/* Lt object */

Mistral_lt::Mistral_lt(Mistral_Expression *var1, Mistral_Expression *var2)
  : Mistral_binop(var1,var2)
{
#ifdef _DEBUGWRAP
  std::cout << "Creating a less than constraint" << std::endl;
#endif
}

Mistral_lt::Mistral_lt(Mistral_Expression *var1, int constant)
  : Mistral_binop(var1,constant)
{
#ifdef _DEBUGWRAP
  std::cout << "Creating a less than constraint" << std::endl;
#endif
}

Mistral_lt::~Mistral_lt()
{
#ifdef _DEBUGWRAP
  std::cout << "delete lessthan" << std::endl;
#endif
}

Mistral_Expression* Mistral_lt::add(MistralSolver *solver, bool top_level) {

  if(!has_been_added()) {
#ifdef _DEBUGWRAP
      std::cout << "add lessthan constraint" << std::endl;
#endif
    _solver = solver;
    
    _vars[0]->add(_solver,false);
    if(_vars[1])  {
      _vars[1]->add(_solver,false);    
      _self = CSP::_Precedence(_vars[0]->_self, 1, 
			       _vars[1]->_self);
    } else {
      _self = CSP::_Precedence(_vars[0]->_self, _constant-1); 
    }
    if( top_level )
      _solver->model->add( _self );
    
  }

  return this;  
}

/* Gt object */

Mistral_gt::Mistral_gt(Mistral_Expression *var1, Mistral_Expression *var2)
  : Mistral_binop(var1,var2)
{
#ifdef _DEBUGWRAP
  std::cout << "Creating a greater than constraint" << std::endl;
#endif
}

Mistral_gt::Mistral_gt(Mistral_Expression *var1, int constant)
  : Mistral_binop(var1,constant)
{
#ifdef _DEBUGWRAP
  std::cout << "Creating a greater than constraint" << std::endl;
#endif
}

Mistral_gt::~Mistral_gt()
{
#ifdef _DEBUGWRAP
  std::cout << "delete greaterthan" << std::endl;
#endif
}

Mistral_Expression* Mistral_gt::add(MistralSolver *solver, bool top_level) {

  //std::cout << "ADD GT" << std::endl;

  if(!has_been_added()) {
#ifdef _DEBUGWRAP
      std::cout << "add greaterthan constraint" << std::endl;
#endif
    _solver = solver;
    
    _vars[0]->add(_solver,false);

    if(_vars[1])  {
      _vars[1]->add(_solver,false);    
      _self = CSP::_Precedence(_vars[1]->_self, 1, 
			       _vars[0]->_self);
    } else {
      _self = CSP::_Precedence(_constant+1, _vars[0]->_self); 
    }
    if( top_level ) 
      _solver->model->add( _self );
    
  }
  
  return this;
}


/* Minimise object */

Mistral_Minimise::Mistral_Minimise(Mistral_Expression *var)
  : Mistral_Expression()
{
#ifdef _DEBUGWRAP
  std::cout << "Creating a minimise objective" << std::endl;
#endif
  _exp = var;
  _obj = NULL;
}

Mistral_Minimise::~Mistral_Minimise()
{
#ifdef _DEBUGWRAP
  std::cout << "delete minimise" << std::endl;
#endif
}

Mistral_Expression* Mistral_Minimise::add(MistralSolver *solver, bool top_level) {

  if(!has_been_added()) {
#ifdef _DEBUGWRAP
      std::cout << "add minimise objective" << std::endl;
#endif
    _solver = solver;
    
    _exp->add(_solver,false);
    _obj = CSP::_Minimise(_exp->_self);

    if( top_level )
      _solver->model->add( _obj );

  }
  
  return this;
}

/* Maximise object */

Mistral_Maximise::Mistral_Maximise(Mistral_Expression *var)
  : Mistral_Expression()
{
#ifdef _DEBUGWRAP
  std::cout << "Creating a maximise objective" << std::endl;
#endif
  _exp = var;
  _obj = NULL;
}

Mistral_Maximise::~Mistral_Maximise()
{
#ifdef _DEBUGWRAP
  std::cout << "delete maximise" << std::endl;
#endif
}

Mistral_Expression* Mistral_Maximise::add(MistralSolver *solver, bool top_level) {

  if(!has_been_added()) {
#ifdef _DEBUGWRAP
      std::cout << "add maximise objective" << std::endl;
#endif
    _solver = solver;
    
    _exp->add(_solver,false);
    _obj = CSP::_Maximise(_exp->_self);

    if( top_level )
      _solver->model->add( _obj );
    
  }
  
  return this;
}

/**************************************************************
 ********************     Solver        ***********************
 **************************************************************/

MistralSolver::MistralSolver()
{
  solver = NULL;
  model = new CSP(); 
  //cb = NULL;
  //parser = NULL;
  feature_ready = false;
  
  nogood_base = NULL;

  first_decision_level = -1;

  is_copy = false;
}

MistralSolver::~MistralSolver()
{
#ifdef _DEBUGWRAP
  std::cout << "c (wrapper) delete solver" << std::endl;
#endif

  if(!is_copy)
    delete model;
  else
    model = NULL;
  delete solver;
  //delete parser;
  //delete cb;
  
}

void MistralSolver::add(Mistral_Expression* arg)
{
  if(arg != NULL)
    arg->add(this, true);
}

void MistralSolver::initialise(MistralExpArray& arg)
{
  BuildObject **args = new BuildObject*[arg.size()];
  for(int i=0; i<arg.size(); ++i)
    args[i] = arg.get_item(i)->_self;
  solver = new Solver(*model, args, arg.size());
  delete [] args;

  saved_level = solver->init_level;
}

void MistralSolver::initialise()
{
  solver = new Solver(*model);

  saved_level = solver->init_level;
}


void MistralSolver::load_xml(const char* filename, const int type) 
{

  cb = new MistralCallback();
  parser = new XMLParser_libxml2<MistralCallback>(*cb); 
  parser->setPreferredExpressionRepresentation(TREE);

  try
    {
      cb->setName         ( "default" );      
      cb->setModel        ( type );
      cb->setInferCliques ( 2 );
      cb->setCliquesLimit ( 10000 );
      cb->setEVarLimit    ( 12000 );
      cb->setEVarSize     ( 150 );
      cb->setRevertToSAT  ( 0 );
      cb->setVerbosity    ( 0 );
      //cb->setFeatureExt   ( 2 );

      parser->parse( filename ); // parse the input file

      model = &(cb->model);
      is_copy = true;

    }
  catch (exception &e)
    {
      cout.flush();
      cerr << "\n\tUnexpected exception :\n";
      cerr << "\t" << e.what() << endl;
      exit(1);
    }

}
 

const char*feature_names[num_features] = {"log_constants", "log_booleans", "log_ranges", 
				"log_bits", "log_lists", "log_values", 
				"log_extra_booleans", "log_extra_ranges", 
				"log_extra_bits", "log_extra_values", 
				"log_constraints", "log_search_vars", 
				"dyn_log_avg_weight", "dyn_log_stdev_weight", 
				"dyn_log_nodes", "dyn_log_propags", "max_arity", 
				"percent_ext", "percent_global", 
				"percent_gac_predicate", "percent_dec_predicate", 
				"sqrt_avg_domsize", "sqrt_max_domsize", 
				"percent_avg_continuity", "percent_min_continuity", 
				"num_alldiff", "perten_binext", "perten_naryext", 
				"perten_largeext", "percent_alldiff", 
				"percent_wsum", "percent_element", 
				"percent_cumulative", "perten_avg_predshape", 
				"perten_avg_predsize", "perten_avg_predarity"};

const char* MistralSolver::get_feature_name(int i) {
  //   std::cout << "c (wrapper) get feature name " << feature_names[i] << std::endl;
  return feature_names[i];
}

double MistralSolver::get_feature(int i) {
  if(!feature_ready) get_features();

//   std::cout << "c (wrapper) get feature value " << feature_names[i] ;
//   std::cout.flush();
//   std::cout << " = " << cb->features_vec[i] << std::endl;
  return cb->features_vec[i];
}

int MistralSolver::num_vars() {
  return cb->numVariables;
}

int MistralSolver::get_degree(int i) {
  BuildObject *x = cb->X[i].var_ptr_;
  if(x)
    return x->getBuildObject()->parent.size;
  return 0;
}


void MistralSolver::extract_graph() {
  VariableInt *x, *y;
  Constraint *con = NULL;
  MistralNode<Constraint*> *nd;
  
  if(!feature_ready) get_features();
  std::cout << "C++ extract graph: " << (solver->variables.size) << std::endl;

  graph.resize(solver->variables.size);
  for(int i=0; i<solver->variables.size; ++i) {
    
    x = solver->variables[i];
    nd = x->constraintsOnValue();
    while( nextNode(nd) ) {
      con = nd->elt;
      for(int i=0; i<con->arity; ++i) {
	y = con->_scope[i];
	graph[i].push_back(y->id);
      }
    }
  }
}


int MistralSolver::numNodes() { return graph.size(); }
int MistralSolver::degree(const int x) { return graph[x].size(); }
int MistralSolver::get_neighbor(const int x, const int y) { return graph[x][y]; }

void MistralSolver::get_features(MistralDoubleArray& features)
{
  if(!feature_ready) get_features();
  for(int i=0; i<num_features; ++i)
    features.add(cb->features_vec[i]);
}

void MistralSolver::get_features()
{
  
  //std::cout << "c (wrapper) extract features " << std::endl; 

  cb->setHeuristic    ( "dom/wldeg" );      
  cb->setRestartPolicy( "dyn" );
  cb->setRestartFactor( (1+1.0/3.0) );
  cb->setRestartBase  ( -1 );
  cb->setRandomize    ( 1 );
  cb->setRandomSeed   ( 11041979 );      
  cb->setTimeLimit    ( 2.0 );
  cb->setNodeLimit    ( -1 );
  
  cb->setProbing      ( 0 );
  cb->setPIteration   ( 100 );
  cb->setPLimit       ( 30 );
  cb->setProbeHeuris  ( 0 );
  cb->setUpdateLW     ( 0 );
  cb->setUpdateIP     ( 0 );
  
  cb->setSAC          ( 0 );
  cb->setLDS          ( 0 );
  
  cb->setDomainSplit  ( 0 );
  
  cb->setFeatureExt   ( 2 );
  cb->setAllSolution  ( 0 );
  
  cb->solve();

  solver = cb->cp_solver;

  feature_ready = true;
  //
}

int MistralSolver::solveAndRestart(const int policy, 
				   const unsigned int base, 
				   const double factor,
				   const double decay,
				   const int reinit)
{
  saved_level = solver->init_level;
  if(solver->init_level < solver->level)
    solver->init_level = solver->level;

  solver->solve_and_restart(policy, base, factor, decay, reinit);
  //solver->reset_trail(false);
  //solver->init_level = saved_level;
  return (is_sat());
}

int MistralSolver::solve()
{
  //std::cout << " SOLVE!! " << solver->init_level << " " << solver->level << std::endl;
  
  saved_level = solver->init_level;
  if(solver->init_level < solver->level)
    solver->init_level = solver->level;
  solver->solve();
  //solver->reset_trail(false);
  //solver->init_level = saved_level;
  return (is_sat());
}

int MistralSolver::startNewSearch()
{
  return solver->startNewSearch();
}

int MistralSolver::getNextSolution()
{
  return solver->getNextSolution();
  //return (solver->getNextSolution() == SAT);
}

int MistralSolver::sacPreprocess(const int type)
{
  return solver->sacPreprocess(type);
}

bool MistralSolver::propagate()
{
  if(solver->level < solver->init_level) {
    solver->STARTTIME = getRunTime();
    return (solver->presolve() != UNSAT);
  }
  return solver->filtering();
}

// void MistralSolver::branch_on(const char* op, Mistral_Expression* x, int v) 
// {
//   if(solver->level < solver->init_level) {
//     solver->STARTTIME = getRunTime();
//     solver->presolve();
//   }
//   if(op[1] == 't') 
//     if(op[0] == 'g') ++v;
//     else --v;
//   solver->newNode(x->_self->varptr_, op[0], v);
// }

int MistralSolver::get_level() { return solver->level; }

int MistralSolver::get_decision_id() { return decisions[decisions.size]->nbj_ident; }

bool MistralSolver::undo(const int nlevel) {
  int lvl = solver->level, okay = true;
  //std::cout << lvl << " -> " << (lvl-nlevel) << " / " << first_decision_level << std::endl;
  if(lvl-nlevel < first_decision_level) okay = false;
  
  lvl -= nlevel;
  if(lvl < 0) lvl = 0;
  solver->backtrackTo(lvl);

  int i=nlevel;
  if(i > valuation_size.size) i = valuation_size.size;
  while(i--) {
    //solver->undo();
    valuation_size.pop(valuation.size);
    decisions.pop();
  }
  return okay;
} 

void MistralSolver::save() 
{
  if(solver->level < solver->init_level) {
    solver->STARTTIME = getRunTime();
    solver->presolve();
  }
  solver->save();
  valuation_size.push(valuation.size);
}

void MistralSolver::post(const char* op, Mistral_Expression* x, int v)
{
  if(op[1] == 't') {
    if(op[0] == 'g') ++v;
    else --v;
  }

  int lvl = solver->level;
  solver->decision[lvl] = x->_self->varptr_; //currentDecision;
  decisions.push(x);
  SimpleUnaryConstraint dec(x->_self->varptr_);
  

  if(first_decision_level < 0) {
    first_decision_level = lvl-1;
    //std::cout << "first = " << first_decision_level << std::endl;
  }

  int vald = v;
  if(op[0] == 'e') {
    solver->decision[lvl]->setDomain(v);
  } else if(op[0] == 'n') {
    solver->decision[lvl]->remove(v);
  } else if(op[0] == 'g') {
    solver->decision[lvl]->setMin(v);
    if(op[1] == 'e') ++vald;
  } else if(op[0] == 'l') {
    solver->decision[lvl]->setMax(v);
    if(op[1] == 't') --vald;
  }

  dec.init_data(op[0], vald);
  solver->branching_decision[lvl] = dec;

  VarValuation vv(v, op[0]);
  valuation.push(vv);

  
  int i=solver->learners.size;
  while( i-- )
    solver->learners[i]->notifyChoice( );

  //solver->post(x->_self->varptr_, op[0], v);
}

void MistralSolver::deduce(const char* op, Mistral_Expression* x, int v) {

  if(op[1] == 't') {
    if(op[0] == 'g') ++v;
    else --v;
  }

  VariableInt *lastDecision = x->_self->varptr_;

  if(op[0] == 'e') {
    lastDecision->setDomain(v);
  } else if(op[0] == 'n') {
    lastDecision->remove(v);
  } else if(op[0] == 'g') {
    lastDecision->setMin(v);
  } else if(op[0] == 'l') {
    lastDecision->setMax(v);
  }

  //solver->deduce();
}

void MistralSolver::deduce() {

  VariableInt *lastDecision = solver->decision[solver->level+1];
  //VarValuation vv = valuation.pop();
  VarValuation vv = valuation[valuation.size];

  if(vv.type == 'e') {
    lastDecision->remove(vv.value);
  } else if(vv.type == 'n') {
    lastDecision->setDomain(vv.value);
  } else if(vv.type == 'g') {
    lastDecision->setMax(vv.value-1);
  } else if(vv.type == 'l') {
    lastDecision->setMin(vv.value+1);
  }


  //solver->deduce();
}

bool MistralSolver::branch_right() {
  //int lvl = solver->level;
  //if(lvl-nlevel < first_decision_level) return false;
  int lvl = solver->level, okay = true;
  //std::cout << lvl << " -> " << (lvl-1) << " / " << first_decision_level << std::endl;
  if(lvl-1 < first_decision_level) okay = false;
  // std::cout << solver->level << " " ;
  //   if(!(solver->undo())) {
  //     std::cout << false << std::endl;
  //     return false;
  //   }
  if(lvl) {
    solver->backtrackTo(lvl-1);
    valuation_size.pop(valuation.size);
    decisions.pop();
    //undo();
    if(okay)
      deduce();
  }
  // std::cout << true << std::endl;
  return okay;
}

// bool MistralSolver::deduce() {
//   if( !solver->undo() ) return false;
//   else solver->deduce();
//   return true;
// }

// bool MistralSolver::deduce() {
//   if( !solver->undo() ) return false;
//   else solver->deduce();
//   return true;
// }

void MistralSolver::store_solution() {
  //std::cout << solver->status << std::endl;
  if(solver->SOLUTIONS == 0 || solver->level != solver->init_level) {
    solver->store_solution();
    if( solver->goal ) solver->status = solver->goal->update();
    solver->closeSearch();
  }
}

void MistralSolver::setHeuristic(const char* var_heuristic, const char* val_heuristic, const int rand)
{
  solver->setHeuristic(var_heuristic, val_heuristic, rand);
}

void MistralSolver::addNogood(MistralExpArray& vars, 
			      MistralIntArray& vals)
{
  int i, n = vars.size();
  if(n != vals.size()) {
    std::cerr << "Warning argument do not align in call  to guide()" <<std::endl;
  }

  if(!nogood_base) {
    nogood_base = new ConstraintNogoodBase(solver);
  }
  //std::cout << std::endl;

  VariableInt* aux;

  bool is_virtual = false;

  Vector< GenLiteral > nogood;
  for(i=0; i<n; ++i) {

    //std::cout << " -- " << i << std::endl;
    //std::cout << " -- " ;
    //vars.get_item(i)->_self->print(std::cout);
    //std::cout << std::endl;
    //std::cout << " -- " << vars.get_item(i)->_self->getVariable() << std::endl;
    
    aux = vars.get_item(i)->_self->getVariable();
    if(aux && aux->getType() != VariableInt::CONST) {
      GenLiteral p(aux, vals.get_item(i));

      if(aux->getType() == VariableInt::VIRTUAL) {
	is_virtual = true;
	p.devirtualize();
      }
      nogood.push(p);
    }
  }

  if(!is_virtual)
    nogood_base->add(nogood);

  // nogood_base->print(std::cout);
//   std::cout << std::endl;


}

void MistralSolver::setAntiLex(MistralExpArray& vars) {
  int i, n = vars.size();
  BuildObject **scope = new BuildObject*[n];
  for(i=0; i<n; ++i) {
    scope[i] = vars.get_item(i)->_self;
  }

  solver->setAntiLex(scope, n);

  delete [] scope;
}

void MistralSolver::guide(MistralExpArray& vars, 
			  MistralIntArray& vals,
			  MistralDoubleArray& probs)
{
  
  int i, n = vars.size();
  if(n != vals.size()) {
    std::cerr << "Warning argument do not align in call  to guide()" <<std::endl;
  }

  BuildObject **scope = new BuildObject*[n];
  int *ideal = new int[n];
  
  for(i=0; i<n; ++i) {
    scope[i] = vars.get_item(i)->_self;
    ideal[i] = vals.get_item(i);
  }

  if(probs.size() == n) {
    int *probabilities = new int[n];
    for(i=0; i<n; ++i) 
      probabilities[i] = (int)(1000.0 * probs.get_item(i));
    solver->setGuidedSplitOrdering(scope, n, ideal, probabilities);
    delete [] probabilities;
  } else {
    solver->setGuidedOrdering(scope, n, ideal);
  }

  delete [] scope;
  delete [] ideal;
}

void MistralSolver::forceFiniteDomain(MistralExpArray& vars)
{
  
  int i, n = vars.size();
  for(i=0; i<n; ++i) 
    vars.get_item(i)->_self->unsetRange();

}

void MistralSolver::backtrackTo(const int level) {
  solver->backtrackTo(level);
}

void MistralSolver::presolve() { solver->presolve(); }
void MistralSolver::increase_init_level(const int i) { solver->init_level += i; }
void MistralSolver::decrease_init_level(const int i) { solver->init_level -= i; }
void MistralSolver::assign(Mistral_Expression *X, const int v) { X->_self->getVariable()->setDomain(v); }

void MistralSolver::upOneLevel() {
  solver->upOneLevel();
}

void MistralSolver::reset(bool full) {
  solver->reset(full);
  solver->init_level = saved_level;
  //if(full) --solver->level;
}

void MistralSolver::setLowerBounds(MistralExpArray& vars, 
				   MistralIntArray& vals) {
  int i, n = vars.size();
  if(n != vals.size()) {
    std::cerr << "Warning argument do not align in call  to guide()" <<std::endl;
  }

  BuildObject **scope = new BuildObject*[n];
  int *lower = new int[n];
  
  for(i=0; i<n; ++i) {
    scope[i] = vars.get_item(i)->_self;
    lower[i] = vals.get_item(i);
  }

  solver->setLowerBounds(scope, n, lower);

  delete [] scope;
  delete [] lower;
}

void MistralSolver::setUpperBounds(MistralExpArray& vars, 
				   MistralIntArray& vals) {
  int i, n = vars.size();
  if(n != vals.size()) {
    std::cerr << "Warning argument do not align in call  to guide()" <<std::endl;
  }

  BuildObject **scope = new BuildObject*[n];
  int *upper = new int[n];
  
  for(i=0; i<n; ++i) {
    scope[i] = vars.get_item(i)->_self;
    upper[i] = vals.get_item(i);
  }

  solver->setUpperBounds(scope, n, upper);

  delete [] scope;
  delete [] upper;
}

void MistralSolver::setRestartNogood() 
{
  solver->setRestartNogood();
}

void MistralSolver::setFailureLimit(const int cutoff)
{
  solver->setFailureLimit(cutoff);
}

void MistralSolver::setNodeLimit(const int cutoff)
{
  solver->setNodeLimit(cutoff);
}

void MistralSolver::setTimeLimit(const int cutoff)
{
  solver->setTimeLimit((double)cutoff);
}

void MistralSolver::setVerbosity(const int degree)
{
  solver->setVerbosity(degree);
}

void MistralSolver::setRandomized(const int degree)
{
  solver->setRandomized(degree);
}

void MistralSolver::setRandomSeed(const int seed)
{
  solver->setRandomSeed(seed);
}

bool MistralSolver::is_opt()
{
  return (solver && solver->status == OPT);
}

bool MistralSolver::is_sat()
{
  return (solver && (solver->status == SAT || solver->status == OPT));
}

bool MistralSolver::is_unsat()
{
  return (solver && solver->status == UNSAT);
}

void MistralSolver::printStatistics()
{
}

int MistralSolver::getBacktracks()
{
  return solver->getBacktracks();
}

int MistralSolver::getNodes()
{
  return solver->getNodes();
}

int MistralSolver::getFailures()
{
  return solver->getFailures();
}

int MistralSolver::getChecks()
{
  return solver->getChecks();
}

int MistralSolver::getPropags()
{
  return solver->getPropags();
}

double MistralSolver::getTime()
{
  return solver->getTime();
}

int MistralSolver::getNumVariables()
{
  return solver->variables.size;
}

int MistralSolver::getNumConstraints()
{
  return solver->constraints.size;
}

void MistralSolver::printPython()
{
  solver->printPython();
}

int MistralSolver::getRandomNumber()
{
  return randint(100);
}

void MistralSolver::test_x60() { solver->test_x60(); }


Mistral_Expression* MistralSolver::get_expression(const int i) 
{
  return new Mistral_Expression(model->toplevel_expressions[i]->getBuildObject());
}
int MistralSolver::num_expression() 
{
  return model->toplevel_expressions.size;
}
int MistralSolver::max_expression_id() 
{
  return model->declarations.size;
}
