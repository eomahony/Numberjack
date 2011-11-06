
/** \file Osi.cpp
    \brief Solver interface for PYTHON Wrapper.
*/

#include "Osi.hpp"
#include <coin/CoinPackedMatrix.hpp>
#include <coin/CoinPackedVector.hpp>
#include <coin/OsiClpSolverInterface.hpp>

/**************************************************************
********************     EXPRESSION        *******************
**************************************************************/

void Osi_Expression::initialise() {
}

Osi_Expression::Osi_Expression() {
#ifdef _DEBUGWRAP
    std::cout << "creating an empty expression" << std::endl;
#endif
    initialise();
}

Osi_Expression::Osi_Expression(const int nval) {
#ifdef _DEBUGWRAP
    std::cout   << "creating a variable [" << 0 << ".." << (nval-1) << "]"
                << std::endl;
#endif
    nbj_ident++;
}

Osi_Expression::Osi_Expression(const int lb, const int ub)
{
  std::cout << "creating a variable [" << lb << ".." << ub << "]" << std::endl;
}

Osi_Expression::Osi_Expression(OsiIntArray& vals)
{
  std::cout << "creating a variable [lb..ub]" << std::endl;
}


int Osi_Expression::getVariableId() const
{
  std::cout << "return identity of expression" << std::endl;
  return 0 ;
}

int Osi_Expression::get_value() const
{
  std::cout << "return value of expression" << std::endl;
  return 0 ;
}

int Osi_Expression::get_size() const
{
  std::cout << "return size of expression" << std::endl;
  return 0;
}

int Osi_Expression::get_max() const
{
  std::cout << "return max of expression" << std::endl;
  return 0;
}

int Osi_Expression::get_min() const
{
  std::cout << "return min of expression" << std::endl;
  return 0;
}

bool Osi_Expression::contain(const int v) const
{
  std::cout << "return min of expression" << std::endl;
  return 0;
}

Osi_Expression::~Osi_Expression()
{
  std::cout << "delete expression" << std::endl;
}

bool Osi_Expression::has_been_added() const
{
  std::cout << "has this expression already been added?" << std::endl;
  return false;
}

Osi_Expression* Osi_Expression::add(OsiSolver *solver, bool top_level)
{
  std::cout << "add expression to model" << std::endl;
  
  if(top_level){
    std::cout << "\tAdding at top level" << std::endl;
  } else {
    std::cout << "\tAdding within tree" <<std::endl;
  }
  
  return this;
}

// /* Binary operators */

Osi_binop::Osi_binop(Osi_Expression *var1,
                          Osi_Expression *var2)
  : Osi_Expression()
{
  std::cout << "creating a binary operator" << std::endl;
  _vars[0] = var1;
  _vars[1] = var2;
}

Osi_binop::Osi_binop(Osi_Expression *var1, int constant)
  : Osi_Expression()
{
  std::cout << "creating a binary (constant) operator" << std::endl;
  _vars[0] = var1;
  _vars[1] = NULL;
}


Osi_binop::~Osi_binop(){
  std::cout << "delete binary operator" << std::endl;
}

/**
* Constraints 
*/

Osi_Min::Osi_Min( OsiExpArray& vars ) 
  : Osi_Expression() 
{
  std::cout << "creating an Min constraint" << std::endl;
}

Osi_Min::Osi_Min( Osi_Expression *var1, Osi_Expression *var2 ) 
  : Osi_Expression() 
{
  std::cout << "creating a binary Min constraint" << std::endl;
}

Osi_Min::~Osi_Min()
{
  std::cout << "delete Min" << std::endl;
}

Osi_Expression* Osi_Min::add(OsiSolver *solver, bool top_level)
{
  if(!has_been_added()) {
    std::cout << "add a Min constraint to solver" << std::endl;
    
    for(int i = 0; i < _vars.size(); ++i)
      _vars.set_item(i, _vars.get_item(i)->add(solver, false));
    
    if(top_level){
      std::cout << "\tAdding at top level" << std::endl;
    } else {
      std::cout << "\tAdding within tree" <<std::endl;
    }
  
  }
  return this;
}


Osi_Max::Osi_Max( OsiExpArray& vars ) 
  : Osi_Expression() 
{
  std::cout << "creating an alldiff constraint" << std::endl;
}

Osi_Max::Osi_Max( Osi_Expression *var1, Osi_Expression *var2 ) 
  : Osi_Expression() 
{
  std::cout << "creating a binary alldiff constraint" << std::endl;
}

Osi_Max::~Osi_Max()
{
  std::cout << "delete alldiff" << std::endl;
}

Osi_Expression* Osi_Max::add(OsiSolver *solver, bool top_level)
{
  if(!has_been_added()) {
    std::cout << "add an alldiff constraint" << std::endl;
    
    for(int i = 0; i < _vars.size(); ++i)
      _vars.set_item(i, _vars.get_item(i)->add(solver, false));
    
    if(top_level){
      std::cout << "\tAdding at top level" << std::endl;
    } else {
      std::cout << "\tAdding within tree" <<std::endl;
    }
    
  }
  return this;
}

Osi_AllDiff::Osi_AllDiff( OsiExpArray& vars ) 
  : Osi_Expression() 
{
  std::cout << "creating an alldiff constraint" << std::endl;
}

Osi_AllDiff::Osi_AllDiff( Osi_Expression *var1, Osi_Expression *var2 ) 
  : Osi_Expression() 
{
  std::cout << "creating a binary alldiff constraint" << std::endl;
}

Osi_AllDiff::~Osi_AllDiff()
{
  std::cout << "delete alldiff" << std::endl;
}

Osi_Expression* Osi_AllDiff::add(OsiSolver *solver, bool top_level)
{
  if(!has_been_added()) {
    std::cout << "add an alldiff constraint" << std::endl;
    
    for(int i = 0; i < _vars.size(); ++i)
      _vars.set_item(i, _vars.get_item(i)->add(solver, false));
    
    if(top_level){
      std::cout << "\tAdding at top level" << std::endl;
    } else {
      std::cout << "\tAdding within tree AllDiff constraint NOT A GOOD IDEA" <<std::endl;
    }
    
  }
  return this;
}

Osi_Gcc::Osi_Gcc(OsiExpArray& vars,
      OsiIntArray& vals,
      OsiIntArray& lb_card,
      OsiIntArray& ub_card)
  : Osi_Expression() 
{
  std::cout << "creating a gcc constraint" << std::endl;
}

Osi_Gcc::~Osi_Gcc()
{
  std::cout << "delete gcc" << std::endl;
}

Osi_Expression* Osi_Gcc::add(OsiSolver *solver, bool top_level)
{
  if(!has_been_added()) {
    std::cout << "add gcc constraint" << std::endl;
      
    for(int i = 0; i < _vars.size(); ++i)
      _vars.set_item(i, _vars.get_item(i)->add(solver, false));
      
      
    if(top_level){
      std::cout << "\tAdding at top level" << std::endl;
    } else {
      std::cout << "\tAdding within tree GCC constraint NOT A GOOD IDEA"
        << std::endl;
    }
      
  }
  return this;
}

Osi_Element::Osi_Element( OsiExpArray& vars ) 
  : Osi_Expression() 
{
  std::cout << "creating element" << std::endl;
}

Osi_Element::~Osi_Element()
{
  std::cout << "delete element" << std::endl;
}

Osi_Expression* Osi_Element::add(OsiSolver *solver, bool top_level)
{
  if(!has_been_added()) {
    std::cout << "add element constraint" << std::endl;
      
    for(int i = 0; i < _vars.size(); ++i)
      _vars.set_item(i, _vars.get_item(i)->add(solver, false));
      
    if(top_level){
      std::cout << "\tAdding at top level" << std::endl;
    } else {
      std::cout << "\tAdding within tree Element NOT A GOOD IDEA"
        << std::endl;
    }
      
  }
  return this;
}

Osi_LeqLex::Osi_LeqLex( OsiExpArray& vars ) 
  : Osi_Expression() 
{
  std::cout << "creating lexleq" << std::endl;
}

Osi_LeqLex::~Osi_LeqLex()
{
  std::cout << "delete leqlex" << std::endl;
}

Osi_Expression* Osi_LeqLex::add(OsiSolver *solver, bool top_level) {
  if(!has_been_added()) {
    std::cout << "add leqlex constraint" << std::endl;
      
    for(int i = 0; i < _vars.size(); ++i)
      _vars.set_item(i, _vars.get_item(i)->add(solver, false));
      
    if(top_level){
      std::cout << "\tAdding at top level" << std::endl;
    } else {
      std::cout << "\tAdding within tree" <<std::endl;
    }
      
  }
  return this;
}


Osi_LessLex::Osi_LessLex( OsiExpArray& vars ) 
  : Osi_Expression() 
{
  std::cout << "creating lexless" << std::endl;
}

Osi_LessLex::~Osi_LessLex()
{
  std::cout << "delete leslex" << std::endl;
}

Osi_Expression* Osi_LessLex::add(OsiSolver *solver, bool top_level) {
  if(!has_been_added()) {
    std::cout << "add lesslex constraint" << std::endl;
      
    for(int i = 0; i < _vars.size(); ++i)
      _vars.set_item(i, _vars.get_item(i)->add(solver, false));
      
    if(top_level){
      std::cout << "\tAdding at top level" << std::endl;
    } else {
      std::cout << "\tAdding within tree" <<std::endl;
    }
      
  }
  return this;
}


Osi_Sum::Osi_Sum(OsiExpArray& vars, 
      OsiIntArray& weights, 
      const int offset)
  : Osi_Expression() 
{
  std::cout << "creating sum" << std::endl;
}

Osi_Sum::Osi_Sum(Osi_Expression *arg1, 
      Osi_Expression *arg2, 
      OsiIntArray& w, 
      const int offset)
  : Osi_Expression() 
{
  std::cout << "creating sum" << std::endl;
}

Osi_Sum::Osi_Sum(Osi_Expression *arg, 
      OsiIntArray& w, 
      const int offset)
  : Osi_Expression() 
{
  std::cout << "creating sum" << std::endl;
}

Osi_Sum::Osi_Sum()
  : Osi_Expression()
{
  //_offset = 0;
}

Osi_Sum::~Osi_Sum(){
  std::cout << "delete sum" << std::endl;
}

void Osi_Sum::addVar(Osi_Expression* v){
  std::cout << "adding variable" << std::endl;
}

void Osi_Sum::addWeight(const int w){
  std::cout << "adding weight" << std::endl;
}

Osi_Expression* Osi_Sum::add(OsiSolver *solver, bool top_level){
  if(!has_been_added()) {
    std::cout << "add sum constraint" << std::endl;
      
    for(int i = 0; i < _vars.size(); ++i)
      _vars.set_item(i, _vars.get_item(i)->add(solver, false));
      
    if(top_level){
      std::cout << "\tAdding at top level" << std::endl;
    } else {
      std::cout << "\tAdding within tree" <<std::endl;
    }
      
  }
  return this;
}

/**
* Binary constraints
*/
Osi_mul::Osi_mul(Osi_Expression *var1, Osi_Expression *var2)
  : Osi_binop(var1,var2)
{
  std::cerr << "s Not supported \nc (mulitply by variable) - exiting" << std::endl;
  exit(1);
}

Osi_mul::Osi_mul(Osi_Expression *var1, int constant)
  : Osi_binop(var1,constant)
{
  std::cout << "creating mul constraint between variable and constant" << std::endl;
}


Osi_mul::~Osi_mul(){
  std::cout << "delete mul constraint" << std::endl;
}

Osi_Expression* Osi_mul::add(OsiSolver *solver, bool top_level){
  if(!has_been_added()) {
      std::cout << "add mul constraint" << std::endl;
      
      _vars[0] = _vars[0]->add(solver, false);
      
      if(top_level){
    std::cout << "\tAdding at top level Mul constraint NOT A GOOD IDEA"
          << std::endl;
      } else {
    std::cout << "\tAdding within tree" <<std::endl;
    
    if(_vars[1]){
      std::cout << "\t\t Adding with two variables" << std::endl;
      
      _vars[1] = _vars[1]->add(solver, false);
      
      /**
      * This is where you put your code
      */
      
    } else {
      
      std::cout << "\t\t Adding with variable and constant" << std::endl;
      
      /**
      * This is where you put your code
      */
      
    }
    
      }
      
  }
  return this;
}

Osi_div::Osi_div(Osi_Expression *var1, Osi_Expression *var2)
  : Osi_binop(var1,var2)
{
  std::cout << "creating div constraint" << std::endl;
}

Osi_div::Osi_div(Osi_Expression *var1, int constant)
  : Osi_binop(var1,constant)
{
  std::cout << "creating div constraint" << std::endl;
}


Osi_div::~Osi_div(){
  std::cout << "delete div" << std::endl;
}

Osi_Expression* Osi_div::add(OsiSolver *solver, bool top_level){
  if(!has_been_added()) {
      std::cout << "add div constraint" << std::endl;
      
      _vars[0] = _vars[0]->add(solver, false);
      
      if(top_level){
    std::cout << "\tAdding at top level Div constraint NOT A GOOD IDEA"
          << std::endl;
    
      } else {
    std::cout << "\tAdding within tree" <<std::endl;
    
    if(_vars[1]){
      
      std::cout << "\t\tAdding two variable binary expression" << std::endl;
      
      _vars[1] = _vars[1]->add(solver, false);
      
      /**
      *  Your code code here
      */
      
    } else {
      
      std::cout << "\t\tAdding in variable and constraint expression" << std::endl;
      
      /**
      *  Your code code here
      */
      
    }
    
      }
  }
  return this;
}

Osi_mod::Osi_mod(Osi_Expression *var1, Osi_Expression *var2)
  : Osi_binop(var1,var2)
{
  std::cout << "creating mod predicate" << std::endl;
}

Osi_mod::Osi_mod(Osi_Expression *var1, int constant)
  : Osi_binop(var1,constant)
{
  std::cout << "creating mod predicate" << std::endl;
}


Osi_mod::~Osi_mod(){
  std::cout << "delete mod" << std::endl;
}

Osi_Expression* Osi_mod::add(OsiSolver *solver, bool top_level){
  if(!has_been_added()) {
      std::cout << "add mod predicate" << std::endl;
      
      if(top_level){
    std::cout << "\tAdding at top level NOT A GOOD IDEA!" << std::endl;
    
    /**
      * Your code goes here. This should probably fail
      */
    
      } else {
    std::cout << "\tAdding within tree" <<std::endl;
    
    if(_vars[1]){
      std::cout << "Adding with two expressions" << std::endl;
      
      _vars[1] = _vars[1]->add(solver, false);
      
      /**
      *  Your code code here
      */
      
    } else {
      
      std::cout << "Adding with an expression and a constant" << std::endl;
      
      /**
      * Your code goes here
      */
      
    }
    
      }
  }
  return this;
}

Osi_and::Osi_and(Osi_Expression *var1,
                      Osi_Expression *var2)
  : Osi_binop(var1,var2)
{
  std::cout << "creating and predicate" << std::endl;
}

Osi_and::Osi_and(Osi_Expression *var1,
                      int constant)
  : Osi_binop(var1,constant)
{
  std::cout << "creating and predicate" << std::endl;
  std::cout << "I DON'T THINK I SHOULD BE HERE" << std::endl;
  
  /**
  * Should never be in this constructor???
  */
  
}

Osi_and::~Osi_and(){
  std::cout << "delete and" << std::endl;
}

Osi_Expression* Osi_and::add(OsiSolver *solver,
                              bool top_level){
  if(!has_been_added()) {
      std::cout << "add and constraint" << std::endl;
      
      _vars[0] = _vars[0]->add(solver, false);
      
      if(top_level){
    std::cout << "\tAdding at top level" << std::endl;
    
    if(_vars[1]){
      
      std::cout << "\tCreating with two variables\n" << std::endl;
      _vars[1] = _vars[1]->add(solver, false);
      
      /**
      * Your code goes here
      */
      
    } else {
      
      std::cout << "Adding with a constant" << std::endl;
      
      /**
      * Your code goes here
      */
      
    }
    
      } else {
    std::cout << "\tAdding within tree" <<std::endl;
    
    if(_vars[1]){
      
      std::cout << "\tAdding with two variables" << std::endl;
      _vars[1] = _vars[1]->add(solver, false);
      
      /**
      * Your code goes here
      */
      
    } else {
      
      std::cout << "\tAdding with a constant" << std::endl;
      
      /**
      * Your code goes here
      */
      
    }
    
      }
      
  }
  return this;
}

Osi_or::Osi_or(Osi_Expression *var1, Osi_Expression *var2)
  : Osi_binop(var1,var2)
{
  std::cout << "creating or predicate" << std::endl;
}

Osi_or::Osi_or(Osi_Expression *var1, int constant)
  : Osi_binop(var1,constant)
{
  std::cout << "creating or predicate" << std::endl;
}

Osi_or::~Osi_or(){
  std::cout << "delete or" << std::endl;
}

Osi_Expression* Osi_or::add(OsiSolver *solver, bool top_level){
  if(!has_been_added()) {
      std::cout << "add or predicate" << std::endl;
      
      _vars[0] = _vars[0]->add(solver, false);
      
      if(top_level){
    std::cout << "\tAdding at top level" << std::endl;
    
    if(_vars[1]){
      std::cout << "\tAdding with two predicates" << std::endl;
      _vars[1] = _vars[1]->add(solver, false);
      
      /**
      * Your code goes here
      */
      
      
    } else {
      
      std::cout << "\tAdding with a constant" << std::endl;
      
      /**
      * Your code goes here
      */
      
    }
    
      } else {
    std::cout << "\tAdding within tree" <<std::endl;
    
    if(_vars[1]){
      std::cout << "\tAdding with two predicates" << std::endl;
      _vars[1] = _vars[1]->add(solver, false);
      
      /**
      * Your code goes here
      */
      
      
    } else {
      
      std::cout << "\tAdding with a constant" << std::endl;
      
      /**
      * Your code goes here
      */
      
    }
    
      }
      
  }
  return this;
}

Osi_eq::Osi_eq(Osi_Expression *var1, Osi_Expression *var2)
  : Osi_binop(var1,var2)
{
  std::cout << "creating equality" << std::endl;
}

Osi_eq::Osi_eq(Osi_Expression *var1, int constant)
  : Osi_binop(var1,constant)
{
  std::cout << "creating equality" << std::endl;
}

Osi_eq::~Osi_eq(){
  std::cout << "delete eq" << std::endl;
}

Osi_Expression* Osi_eq::add(OsiSolver *solver, bool top_level){
  if(!has_been_added()) {
      std::cout << "add equality constraint" << std::endl;
      
      _vars[0] = _vars[0]->add(solver, false);  
    
      if(top_level){
    std::cout << "\tAdding at top level" << std::endl;
      
    if(_vars[1]){
      std::cout << "\tAdding with two predicates" << std::endl;
      _vars[1] = _vars[1]->add(solver, false);
      
      /**
      * Your code goes here
      */
      
      
    } else {
      
      std::cout << "\tAdding with a constant" << std::endl;
      
      /**
      * Your code goes here
      */
      
    }
      
      } else {
    std::cout << "\tAdding within tree" <<std::endl;
      
    if(_vars[1]){
      std::cout << "\tAdding with two predicates" << std::endl;
      _vars[1] = _vars[1]->add(solver, false);
      
      /**
      * Your code goes here
      */
      
      
    } else {
      
      std::cout << "\tAdding with a constant" << std::endl;
      
      /**
      * Your code goes here
      */
      
    }
      
      }
      
  }
  return this;
}

/* Disequality operator */

Osi_ne::Osi_ne(Osi_Expression *var1, Osi_Expression *var2)
  : Osi_binop(var1,var2)
{
  std::cout << "creating notequal" << std::endl;
}

Osi_ne::Osi_ne(Osi_Expression *var1, int constant)
  : Osi_binop(var1,constant)
{
  std::cout << "creating notequal" << std::endl;
}

Osi_ne::~Osi_ne(){
  std::cout << "delete notequal" << std::endl;
}

Osi_Expression* Osi_ne::add(OsiSolver *solver, bool top_level){
  if(!has_been_added()) {
      std::cout << "add notequal constraint" << std::endl;
      
      _vars[0] = _vars[0]->add(solver, false);
      
      if(top_level){
    std::cout << "\tAdding at top level" << std::endl;
      
    if(_vars[1]){
      std::cout << "\tAdding with two predicates" << std::endl;
      _vars[1] = _vars[1]->add(solver, false);
      
      /**
      * Your code goes here
      */
      
      
    } else {
      
      std::cout << "\tAdding with a constant" << std::endl;
      
      /**
      * Your code goes here
      */
      
    }
      
      } else {
    std::cout << "\tAdding within tree" <<std::endl;
    if(_vars[1]){
      std::cout << "\tAdding with two predicates" << std::endl;
      _vars[1] = _vars[1]->add(solver, false);
      
      /**
      * Your code goes here
      */
      
      
    } else {
      
      std::cout << "\tAdding with a constant" << std::endl;
      
      /**
      * Your code goes here
      */
      
    }
      
      }
      
  }
  return this;
}

/* Leq operator */

Osi_le::Osi_le(Osi_Expression *var1, Osi_Expression *var2)
  : Osi_binop(var1,var2)
{ 
  std::cout << "Creating less than or equal" << std::endl;
}

Osi_le::Osi_le(Osi_Expression *var1, int constant)
  : Osi_binop(var1,constant)
{
  std::cout << "Creating less than or equal" << std::endl;
}


Osi_le::~Osi_le(){
  std::cout << "delete lessequal" << std::endl;
}

Osi_Expression* Osi_le::add(OsiSolver *solver, bool top_level) {
  if(!has_been_added()) {
      std::cout << "add lessequal constraint" << std::endl;
      
      if(top_level){
    std::cout << "\tAdding at top level" << std::endl;
      
    if(_vars[1]){
      std::cout << "\tAdding with two predicates" << std::endl;
      _vars[1] = _vars[1]->add(solver, false);
      
      /**
      * Your code goes here
      */
      
    } else {
      
      std::cout << "\tAdding with a constant" << std::endl;
      
      /**
      * Your code goes here
      */
      
    }
      
      } else {
    std::cout << "\tAdding within tree" <<std::endl;
    
    if(_vars[1]){
      std::cout << "\tAdding with two predicates" << std::endl;
      _vars[1] = _vars[1]->add(solver, false);
      
      /**
      * Your code goes here
      */
      
      
    } else {
      
      std::cout << "\tAdding with a constant" << std::endl;
      
      /**
      * Your code goes here
      */
      
    }
    
      }
  }
  return this;
}

/* Geq operator */

Osi_ge::Osi_ge(Osi_Expression *var1, Osi_Expression *var2)
  : Osi_binop(var1,var2)
{
  std::cout << "Creating a greater or equal constraint" << std::endl;
}

Osi_ge::Osi_ge(Osi_Expression *var1, int constant)
  : Osi_binop(var1,constant)
{
  std::cout << "Creating a greater or equal constraint" << std::endl;
}

Osi_ge::~Osi_ge(){
  std::cout << "delete greaterequal" << std::endl;
}

Osi_Expression* Osi_ge::add(OsiSolver *solver, bool top_level){ 
  if(!has_been_added()) {
      std::cout << "add greaterequal constraint" << std::endl;
      
      _vars[0] = _vars[0]->add(solver, false);
      
      if(top_level){
    std::cout << "\tAdding at top level" << std::endl;
    if(_vars[1]){
      std::cout << "\tAdding with two predicates" << std::endl;
      _vars[1] = _vars[1]->add(solver, false);
      
      /**
      * Your code goes here
      */
      
      
    } else {
      
      std::cout << "\tAdding with a constant" << std::endl;
      
      /**
      * Your code goes here
      */
      
    }
    
      } else {
    std::cout << "\tAdding within tree" <<std::endl;
    if(_vars[1]){
      std::cout << "\tAdding with two predicates" << std::endl;
      _vars[1] = _vars[1]->add(solver, false);
      
      /**
      * Your code goes here
      */
      
      
    } else {
      
      std::cout << "\tAdding with a constant" << std::endl;
      
      /**
      * Your code goes here
      */
      
    }
    
      }
      
  }

  return this;
}

/* Lt object */

Osi_lt::Osi_lt(Osi_Expression *var1, Osi_Expression *var2)
  : Osi_binop(var1,var2)
{
  std::cout << "Creating a less than constraint" << std::endl;
}

Osi_lt::Osi_lt(Osi_Expression *var1, int constant)
  : Osi_binop(var1,constant)
{
  std::cout << "Creating a less than constraint" << std::endl;
}

Osi_lt::~Osi_lt()
{
  std::cout << "delete lessthan" << std::endl;
}

Osi_Expression* Osi_lt::add(OsiSolver *solver, bool top_level) {
  if(!has_been_added()) {
      std::cout << "add lessthan constraint" << std::endl;
      
      _vars[0] = _vars[0]->add(solver, false);
      
      if(top_level){
    std::cout << "\tAdding at top level" << std::endl;
    
    if(_vars[1]){
      std::cout << "\tAdding with two predicates" << std::endl;
      _vars[1] = _vars[1]->add(solver, false);
      
      /**
      * Your code goes here
      */
      
      
    } else {
      
      std::cout << "\tAdding with a constant" << std::endl;
      
      /**
      * Your code goes here
      */
      
    }
    
      } else {
    std::cout << "\tAdding within tree" <<std::endl;
    
    if(_vars[1]){
      std::cout << "\tAdding with two predicates" << std::endl;
      _vars[1] = _vars[1]->add(solver, false);
      
      /**
      * Your code goes here
      */
      
      
    } else {
      
      std::cout << "\tAdding with a constant" << std::endl;
      
      /**
      * Your code goes here
      */
      
    }
      }
      
  }
  return this;  
}

/* Gt object */

Osi_gt::Osi_gt(Osi_Expression *var1, Osi_Expression *var2)
  : Osi_binop(var1,var2)
{
  std::cout << "Creating a greater than constraint" << std::endl;
}

Osi_gt::Osi_gt(Osi_Expression *var1, int constant)
  : Osi_binop(var1,constant)
{
  std::cout << "Creating a greater than constraint" << std::endl;
}

Osi_gt::~Osi_gt()
{
  std::cout << "delete greaterthan" << std::endl;
}

Osi_Expression* Osi_gt::add(OsiSolver *solver, bool top_level) {
  if(!has_been_added()) {
      std::cout << "add greaterthan constraint" << std::endl;
      
      _vars[0] = _vars[0]->add(solver, false);
      
      if(top_level){
    std::cout << "\tAdding at top level" << std::endl;
    
    if(_vars[1]){
      std::cout << "\tAdding with two predicates" << std::endl;
      _vars[1] = _vars[1]->add(solver, false);
      
      /**
      * Your code goes here
      */
      
    } else {
      
      std::cout << "\tAdding with a constant" << std::endl;
      
      /**
      * Your code goes here
      */
      
    }
    
      } else {
    std::cout << "\tAdding within tree" <<std::endl;
    if(_vars[1]){
      std::cout << "\tAdding with two predicates" << std::endl;
      _vars[1] = _vars[1]->add(solver, false);
      
      /**
      * Your code goes here
      */
      
    } else {
      
      std::cout << "\tAdding with a constant" << std::endl;
      
      /**
      * Your code goes here
      s*/ 
      
    }
    
      }
      
  }
  return this;
}


/* Minimise object */

Osi_Minimise::Osi_Minimise(Osi_Expression *var)
  : Osi_Expression()
{
  std::cout << "Creating a minimise objective" << std::endl;
}

Osi_Minimise::~Osi_Minimise()
{
  std::cout << "delete minimise" << std::endl;
}

Osi_Expression* Osi_Minimise::add(OsiSolver *solver, bool top_level) {
  if(!has_been_added()) {
      std::cout << "add minimise objective" << std::endl;
      
      // This will be the objective function
      _exp = _exp->add(solver, false);
      
      if(top_level){
    std::cout << "\tAdding at top level" << std::endl;
    
    /**
    * Your code goes here
    */
    
      } else {
    std::cout << "\tAdding within tree" <<std::endl;
    std::cout << "SHOULD NOT BE HERE!!!" <<std::endl;
      }
      
  }
  return this;
}

/* Maximise object */

Osi_Maximise::Osi_Maximise(Osi_Expression *var)
  : Osi_Expression()
{
  std::cout << "Creating a maximise objective" << std::endl;
}

Osi_Maximise::~Osi_Maximise()
{
  std::cout << "delete maximise" << std::endl;
}

Osi_Expression* Osi_Maximise::add(OsiSolver *solver, bool top_level) {
  if(!has_been_added()) {
      std::cout << "add maximise objective" << std::endl;
      
      // This will be the objective function
      _exp = _exp->add(solver, false);
      
      if(top_level){
    std::cout << "\tAdding at top level" << std::endl;
    
    /**
    * Your code goes here
    */
    
      } else {
    std::cout << "\tAdding within tree" <<std::endl;
    std::cout << "SHOULD NOT BE HERER!!!" << std::endl;
      }
      
  }
  return this;
}

/**************************************************************
********************     Solver        ***********************
**************************************************************/

OsiSolver::OsiSolver()
{
  std::cout << "c (wrapper) creating solver" << std::endl;
}

OsiSolver::~OsiSolver()
{
  std::cout << "c (wrapper) delete solver" << std::endl;
}

void OsiSolver::add(Osi_Expression* arg)
{
  std::cout << "adding expression" << std::endl;
  
  arg->add(this, true);
  
}

void OsiSolver::initialise(OsiExpArray& arg)
{
  std::cout << "Initialising solver with array of expressions" << std::endl;
}

void OsiSolver::initialise()
{
  std::cout << "Initialising solver with no expressions" << std::endl;
}

int OsiSolver::num_vars() {
  std::cout << "return number of variables" << std::endl;
  return 0;
}

int OsiSolver::get_degree(int i) {
  std::cout << "return degree of expression i" << std::endl;
  return 0;
}

int OsiSolver::solveAndRestart(const int policy, 
      const unsigned int base, 
      const double factor,
      const double decay,
      const int reinit)
{
  std::cout << " SOLVE!! at level: ..." << std::endl;
  return 0;
}

int OsiSolver::solve()
{
  std::cout << " SOLVE!! " << std::endl;
  return 0;
}

int OsiSolver::startNewSearch()
{
  std::cout << "starting new search" << std::endl;
  return 0;
}

int OsiSolver::getNextSolution()
{
  std::cout << "getting next solution" << std::endl;
  return 0;
}

int OsiSolver::sacPreprocess(const int type)
{
  std::cout << "running sac preprocessing" << std::endl;
  return 0;
}

bool OsiSolver::propagate()
{
  std::cout << "propagating" << std::endl;
  return 0;
}

int OsiSolver::get_level()
{
  std::cout << "return current level of solver" << std::endl;
  return 0;
}

bool OsiSolver::undo(const int nlevel) {
  std::cout << "Going back nlevel levels" << std::endl;
  return 0;
} 

void OsiSolver::save() 
{
  std::cout << "saving state" << std::endl;
}

int OsiSolver::next(Osi_Expression* x, int v)
{
  std::cout << "get next" << std::endl;
  return 0;
}

void OsiSolver::post(const char* op, Osi_Expression* x, int v)
{
  std::cout << "Posting a constraint" << std::endl;
}

void OsiSolver::deduce(const char* op, Osi_Expression* x, int v) {
  std::cout << "deducing" << std::endl;
}

void OsiSolver::deduce() {
  std::cout << "deducing" << std::endl;
}

bool OsiSolver::branch_right() {
  std::cout << "branching right" << std::endl;
  return 0;
}

void OsiSolver::store_solution() {
  std::cout << "storing solution" << std::endl;
}

void OsiSolver::setHeuristic(const char* var_heuristic, const char* val_heuristic, const int rand)
{
  std::cout << "Setting heuristics" << std::endl;
}

void OsiSolver::addNogood(OsiExpArray& vars, 
        OsiIntArray& vals)
{
  std::cout << "Adding nogood constraint" <<std::endl;
}

void OsiSolver::guide(OsiExpArray& vars, 
    OsiIntArray& vals,
    OsiDoubleArray& probs)
{
  std::cout << "Adding nogood constraint" <<std::endl;
}

void OsiSolver::forceFiniteDomain(OsiExpArray& vars)
{
  std::cout << "Forcing finite domain" <<std::endl;
}

void OsiSolver::backtrackTo(const int level) {
  std::cout << "backtracking to level" <<std::endl;
}

void OsiSolver::presolve()
{
  std::cout << "presolving" <<std::endl;
}

void OsiSolver::increase_init_level(const int i)
{
  std::cout << "increasing initial level" <<std::endl;
}

void OsiSolver::decrease_init_level(const int i)
{
  std::cout << "decreasing initial level" <<std::endl;
}

void OsiSolver::assign(Osi_Expression *X, const int v) {
  std::cout << "Setting domain of expression X to v" <<std::endl;
}

void OsiSolver::upOneLevel() {
  std::cout << "stepping up one level" <<std::endl;
}

void OsiSolver::reset(bool full) {
  std::cout << "resetting solver" <<std::endl;
}

void OsiSolver::setLowerBounds(OsiExpArray& vars, 
      OsiIntArray& vals) {
  std::cout << "stepping up one level" <<std::endl;
}

void OsiSolver::setUpperBounds(OsiExpArray& vars, 
      OsiIntArray& vals) {
  std::cout << "stetting upper bounds" <<std::endl;
}

void OsiSolver::setRestartNogood() 
{
  std::cout << "setting restart no good" <<std::endl;
}

void OsiSolver::setFailureLimit(const int cutoff)
{
  std::cout << "setting failure limit" <<std::endl;
}

void OsiSolver::setNodeLimit(const int cutoff)
{
  std::cout << "setting node limit" <<std::endl;
}

void OsiSolver::setTimeLimit(const int cutoff)
{
  std::cout << "setting time limit" <<std::endl;
}

void OsiSolver::setVerbosity(const int degree)
{
  std::cout << "setting verbosity" <<std::endl;
}

void OsiSolver::setRandomized(const int degree)
{
  std::cout << "setting randomised" <<std::endl;
}

void OsiSolver::setRandomSeed(const int seed)
{
  std::cout << "setting random seed" <<std::endl;
}

bool OsiSolver::is_opt()
{
  std::cout << "returning is_opt()" <<std::endl;
  return true;
}

bool OsiSolver::is_sat()
{
  std::cout << "returning is satisfied?" <<std::endl;
  return true;
}

bool OsiSolver::is_unsat()
{
  std::cout << "returning is NOT satisfied?" <<std::endl;
  return true;
}

void OsiSolver::printStatistics()
{
  std::cout << "printing statistics" <<std::endl;
}

int OsiSolver::getBacktracks()
{
  std::cout << "return number of backtracks" <<std::endl;
  return 0;
}

int OsiSolver::getNodes()
{
  std::cout << "return number of nodes" <<std::endl;
  return 0;
}

int OsiSolver::getFailures()
{
  std::cout << "return number of failures" <<std::endl;
  return 0;
}

int OsiSolver::getChecks()
{
  std::cout << "return number of checks" <<std::endl;
  return 0;
}

int OsiSolver::getPropags()
{
  std::cout << "return amount of propagation" <<std::endl;
  return 0;
}

double OsiSolver::getTime()
{
  std::cout << "return duration" <<std::endl;
  return 0;
}

int OsiSolver::getNumVariables()
{
  std::cout << "return number of variables" <<std::endl;
  return 0;
}

int OsiSolver::getNumConstraints()
{
  std::cout << "return number of constraints" <<std::endl;
  return 0;
}

int OsiSolver::getRandomNumber()
{
  return 8;
}
