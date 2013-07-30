
/** \file ExampleInterface.cpp
    \brief Solver interface for PYTHON Wrapper.
*/

#include "ExampleInterface.hpp"

/**************************************************************
 ********************     EXPRESSION        *******************
 **************************************************************/

ExampleInterface_Expression::ExampleInterface_Expression()
{
  std::cout << "creating an boolean expression" << std::endl;
}

ExampleInterface_Expression::ExampleInterface_Expression(const int nval)
{
  std::cout << "creating a variable [" << 0 << ".." << (nval-1) << "]" << std::endl;
}

ExampleInterface_Expression::ExampleInterface_Expression(const int lb, const int ub)
{
  std::cout << "creating a variable [" << lb << ".." << ub << "]" << std::endl;
}

ExampleInterface_Expression::ExampleInterface_Expression(ExampleInterfaceIntArray& vals)
{
  std::cout << "creating a variable [lb..ub]" << std::endl;
}


int ExampleInterface_Expression::getVariableId() const
{
  std::cout << "return identity of expression" << std::endl;
  return 0 ;
}

int ExampleInterface_Expression::get_value() const
{
  std::cout << "return value of expression" << std::endl;
  return 0 ;
}

int ExampleInterface_Expression::get_size() const
{
  std::cout << "return size of expression" << std::endl;
  return 0;
}

int ExampleInterface_Expression::get_max() const
{
  std::cout << "return max of expression" << std::endl;
  return 0;
}

int ExampleInterface_Expression::get_min() const
{
  std::cout << "return min of expression" << std::endl;
  return 0;
}

bool ExampleInterface_Expression::contain(const int v) const
{
  std::cout << "return min of expression" << std::endl;
  return 0;
}

ExampleInterface_Expression::~ExampleInterface_Expression()
{
  std::cout << "delete expression" << std::endl;
}

bool ExampleInterface_Expression::has_been_added() const
{
  std::cout << "has this expression already been added?" << std::endl;
  return false;
}

ExampleInterface_Expression* ExampleInterface_Expression::add(ExampleInterfaceSolver *solver, bool top_level)
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

ExampleInterface_binop::ExampleInterface_binop(ExampleInterface_Expression *var1,
                                               ExampleInterface_Expression *var2)
  : ExampleInterface_Expression()
{
  std::cout << "creating a binary operator" << std::endl;
  _vars[0] = var1;
  _vars[1] = var2;
}

ExampleInterface_binop::ExampleInterface_binop(ExampleInterface_Expression *var1, int constant)
  : ExampleInterface_Expression()
{
  std::cout << "creating a binary (constant) operator" << std::endl;
  _vars[0] = var1;
  _vars[1] = NULL;
}


ExampleInterface_binop::~ExampleInterface_binop(){
  std::cout << "delete binary operator" << std::endl;
}

/**
 * Constraints 
 */

ExampleInterface_Min::ExampleInterface_Min( ExampleInterfaceExpArray& vars ) 
  : ExampleInterface_Expression() 
{
  std::cout << "creating an Min constraint" << std::endl;
}

ExampleInterface_Min::ExampleInterface_Min( ExampleInterface_Expression *var1, ExampleInterface_Expression *var2 ) 
  : ExampleInterface_Expression() 
{
  std::cout << "creating a binary Min constraint" << std::endl;
}

ExampleInterface_Min::~ExampleInterface_Min()
{
  std::cout << "delete Min" << std::endl;
}

ExampleInterface_Expression* ExampleInterface_Min::add(ExampleInterfaceSolver *solver, bool top_level)
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


ExampleInterface_Max::ExampleInterface_Max( ExampleInterfaceExpArray& vars ) 
  : ExampleInterface_Expression() 
{
  std::cout << "creating an alldiff constraint" << std::endl;
}

ExampleInterface_Max::ExampleInterface_Max( ExampleInterface_Expression *var1, ExampleInterface_Expression *var2 ) 
  : ExampleInterface_Expression() 
{
  std::cout << "creating a binary alldiff constraint" << std::endl;
}

ExampleInterface_Max::~ExampleInterface_Max()
{
  std::cout << "delete alldiff" << std::endl;
}

ExampleInterface_Expression* ExampleInterface_Max::add(ExampleInterfaceSolver *solver, bool top_level)
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

ExampleInterface_AllDiff::ExampleInterface_AllDiff( ExampleInterfaceExpArray& vars ) 
  : ExampleInterface_Expression() 
{
  std::cout << "creating an alldiff constraint" << std::endl;
}

ExampleInterface_AllDiff::ExampleInterface_AllDiff( ExampleInterface_Expression *var1, ExampleInterface_Expression *var2 ) 
  : ExampleInterface_Expression() 
{
  std::cout << "creating a binary alldiff constraint" << std::endl;
}

ExampleInterface_AllDiff::~ExampleInterface_AllDiff()
{
  std::cout << "delete alldiff" << std::endl;
}

ExampleInterface_Expression* ExampleInterface_AllDiff::add(ExampleInterfaceSolver *solver, bool top_level)
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

ExampleInterface_Gcc::ExampleInterface_Gcc(ExampleInterfaceExpArray& vars,
       ExampleInterfaceIntArray& vals,
       ExampleInterfaceIntArray& lb_card,
       ExampleInterfaceIntArray& ub_card)
  : ExampleInterface_Expression() 
{
  std::cout << "creating a gcc constraint" << std::endl;
}

ExampleInterface_Gcc::~ExampleInterface_Gcc()
{
  std::cout << "delete gcc" << std::endl;
}

ExampleInterface_Expression* ExampleInterface_Gcc::add(ExampleInterfaceSolver *solver, bool top_level)
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

ExampleInterface_Element::ExampleInterface_Element( ExampleInterfaceExpArray& vars ) 
  : ExampleInterface_Expression() 
{
  std::cout << "creating element" << std::endl;
}

ExampleInterface_Element::~ExampleInterface_Element()
{
  std::cout << "delete element" << std::endl;
}

ExampleInterface_Expression* ExampleInterface_Element::add(ExampleInterfaceSolver *solver, bool top_level)
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

ExampleInterface_LeqLex::ExampleInterface_LeqLex( ExampleInterfaceExpArray& vars ) 
  : ExampleInterface_Expression() 
{
  std::cout << "creating lexleq" << std::endl;
}

ExampleInterface_LeqLex::~ExampleInterface_LeqLex()
{
  std::cout << "delete leqlex" << std::endl;
}

ExampleInterface_Expression* ExampleInterface_LeqLex::add(ExampleInterfaceSolver *solver, bool top_level) {
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


ExampleInterface_LessLex::ExampleInterface_LessLex( ExampleInterfaceExpArray& vars ) 
  : ExampleInterface_Expression() 
{
  std::cout << "creating lexless" << std::endl;
}

ExampleInterface_LessLex::~ExampleInterface_LessLex()
{
  std::cout << "delete leslex" << std::endl;
}

ExampleInterface_Expression* ExampleInterface_LessLex::add(ExampleInterfaceSolver *solver, bool top_level) {
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
 

ExampleInterface_Sum::ExampleInterface_Sum(ExampleInterfaceExpArray& vars, 
       ExampleInterfaceIntArray& weights, 
       const int offset)
  : ExampleInterface_Expression() 
{
  std::cout << "creating sum" << std::endl;
}

ExampleInterface_Sum::ExampleInterface_Sum(ExampleInterface_Expression *arg1, 
       ExampleInterface_Expression *arg2, 
       ExampleInterfaceIntArray& w, 
       const int offset)
  : ExampleInterface_Expression() 
{
  std::cout << "creating sum" << std::endl;
}

ExampleInterface_Sum::ExampleInterface_Sum(ExampleInterface_Expression *arg, 
       ExampleInterfaceIntArray& w, 
       const int offset)
  : ExampleInterface_Expression() 
{
  std::cout << "creating sum" << std::endl;
}

ExampleInterface_Sum::ExampleInterface_Sum()
  : ExampleInterface_Expression()
{
  //_offset = 0;
}

ExampleInterface_Sum::~ExampleInterface_Sum(){
  std::cout << "delete sum" << std::endl;
}

void ExampleInterface_Sum::addVar(ExampleInterface_Expression* v){
  std::cout << "adding variable" << std::endl;
}

void ExampleInterface_Sum::addWeight(const int w){
  std::cout << "adding weight" << std::endl;
}

ExampleInterface_Expression* ExampleInterface_Sum::add(ExampleInterfaceSolver *solver, bool top_level){
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
ExampleInterface_mul::ExampleInterface_mul(ExampleInterface_Expression *var1, ExampleInterface_Expression *var2)
  : ExampleInterface_binop(var1,var2)
{
  std::cout << "creating mul constraint between two variables" << std::endl;
}

ExampleInterface_mul::ExampleInterface_mul(ExampleInterface_Expression *var1, int constant)
  : ExampleInterface_binop(var1,constant)
{
  std::cout << "creating mul constraint between variable and constant" << std::endl;
}


ExampleInterface_mul::~ExampleInterface_mul(){
  std::cout << "delete mul constraint" << std::endl;
}

ExampleInterface_Expression* ExampleInterface_mul::add(ExampleInterfaceSolver *solver, bool top_level){
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

ExampleInterface_div::ExampleInterface_div(ExampleInterface_Expression *var1, ExampleInterface_Expression *var2)
  : ExampleInterface_binop(var1,var2)
{
  std::cout << "creating div constraint" << std::endl;
}

ExampleInterface_div::ExampleInterface_div(ExampleInterface_Expression *var1, int constant)
  : ExampleInterface_binop(var1,constant)
{
  std::cout << "creating div constraint" << std::endl;
}


ExampleInterface_div::~ExampleInterface_div(){
  std::cout << "delete div" << std::endl;
}

ExampleInterface_Expression* ExampleInterface_div::add(ExampleInterfaceSolver *solver, bool top_level){
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

ExampleInterface_mod::ExampleInterface_mod(ExampleInterface_Expression *var1, ExampleInterface_Expression *var2)
  : ExampleInterface_binop(var1,var2)
{
  std::cout << "creating mod predicate" << std::endl;
}

ExampleInterface_mod::ExampleInterface_mod(ExampleInterface_Expression *var1, int constant)
  : ExampleInterface_binop(var1,constant)
{
  std::cout << "creating mod predicate" << std::endl;
}


ExampleInterface_mod::~ExampleInterface_mod(){
  std::cout << "delete mod" << std::endl;
}

ExampleInterface_Expression* ExampleInterface_mod::add(ExampleInterfaceSolver *solver, bool top_level){
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

ExampleInterface_and::ExampleInterface_and(ExampleInterface_Expression *var1,
                                           ExampleInterface_Expression *var2)
  : ExampleInterface_binop(var1,var2)
{
  std::cout << "creating and predicate" << std::endl;
}

ExampleInterface_and::ExampleInterface_and(ExampleInterface_Expression *var1,
                                           int constant)
  : ExampleInterface_binop(var1,constant)
{
  std::cout << "creating and predicate" << std::endl;
  std::cout << "I DON'T THINK I SHOULD BE HERE" << std::endl;
  
  /**
   * Should never be in this constructor???
   */
  
}

ExampleInterface_and::~ExampleInterface_and(){
  std::cout << "delete and" << std::endl;
}

ExampleInterface_Expression* ExampleInterface_and::add(ExampleInterfaceSolver *solver,
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

ExampleInterface_or::ExampleInterface_or(ExampleInterface_Expression *var1, ExampleInterface_Expression *var2)
  : ExampleInterface_binop(var1,var2)
{
  std::cout << "creating or predicate" << std::endl;
}

ExampleInterface_or::ExampleInterface_or(ExampleInterface_Expression *var1, int constant)
  : ExampleInterface_binop(var1,constant)
{
  std::cout << "creating or predicate" << std::endl;
}

ExampleInterface_or::~ExampleInterface_or(){
  std::cout << "delete or" << std::endl;
}

ExampleInterface_Expression* ExampleInterface_or::add(ExampleInterfaceSolver *solver, bool top_level){
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

ExampleInterface_eq::ExampleInterface_eq(ExampleInterface_Expression *var1, ExampleInterface_Expression *var2)
  : ExampleInterface_binop(var1,var2)
{
  std::cout << "creating equality" << std::endl;
}

ExampleInterface_eq::ExampleInterface_eq(ExampleInterface_Expression *var1, int constant)
  : ExampleInterface_binop(var1,constant)
{
  std::cout << "creating equality" << std::endl;
}

ExampleInterface_eq::~ExampleInterface_eq(){
  std::cout << "delete eq" << std::endl;
}

ExampleInterface_Expression* ExampleInterface_eq::add(ExampleInterfaceSolver *solver, bool top_level){
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

ExampleInterface_ne::ExampleInterface_ne(ExampleInterface_Expression *var1, ExampleInterface_Expression *var2)
  : ExampleInterface_binop(var1,var2)
{
  std::cout << "creating notequal" << std::endl;
}

ExampleInterface_ne::ExampleInterface_ne(ExampleInterface_Expression *var1, int constant)
  : ExampleInterface_binop(var1,constant)
{
  std::cout << "creating notequal" << std::endl;
}

ExampleInterface_ne::~ExampleInterface_ne(){
  std::cout << "delete notequal" << std::endl;
}

ExampleInterface_Expression* ExampleInterface_ne::add(ExampleInterfaceSolver *solver, bool top_level){
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

ExampleInterface_le::ExampleInterface_le(ExampleInterface_Expression *var1, ExampleInterface_Expression *var2)
  : ExampleInterface_binop(var1,var2)
{ 
  std::cout << "Creating less than or equal" << std::endl;
}

ExampleInterface_le::ExampleInterface_le(ExampleInterface_Expression *var1, int constant)
  : ExampleInterface_binop(var1,constant)
{
  std::cout << "Creating less than or equal" << std::endl;
}


ExampleInterface_le::~ExampleInterface_le(){
  std::cout << "delete lessequal" << std::endl;
}

ExampleInterface_Expression* ExampleInterface_le::add(ExampleInterfaceSolver *solver, bool top_level) {
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

ExampleInterface_ge::ExampleInterface_ge(ExampleInterface_Expression *var1, ExampleInterface_Expression *var2)
  : ExampleInterface_binop(var1,var2)
{
  std::cout << "Creating a greater or equal constraint" << std::endl;
}

ExampleInterface_ge::ExampleInterface_ge(ExampleInterface_Expression *var1, int constant)
  : ExampleInterface_binop(var1,constant)
{
  std::cout << "Creating a greater or equal constraint" << std::endl;
}

ExampleInterface_ge::~ExampleInterface_ge(){
  std::cout << "delete greaterequal" << std::endl;
}

ExampleInterface_Expression* ExampleInterface_ge::add(ExampleInterfaceSolver *solver, bool top_level){ 
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

ExampleInterface_lt::ExampleInterface_lt(ExampleInterface_Expression *var1, ExampleInterface_Expression *var2)
  : ExampleInterface_binop(var1,var2)
{
  std::cout << "Creating a less than constraint" << std::endl;
}

ExampleInterface_lt::ExampleInterface_lt(ExampleInterface_Expression *var1, int constant)
  : ExampleInterface_binop(var1,constant)
{
  std::cout << "Creating a less than constraint" << std::endl;
}

ExampleInterface_lt::~ExampleInterface_lt()
{
  std::cout << "delete lessthan" << std::endl;
}

ExampleInterface_Expression* ExampleInterface_lt::add(ExampleInterfaceSolver *solver, bool top_level) {
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

ExampleInterface_gt::ExampleInterface_gt(ExampleInterface_Expression *var1, ExampleInterface_Expression *var2)
  : ExampleInterface_binop(var1,var2)
{
  std::cout << "Creating a greater than constraint" << std::endl;
}

ExampleInterface_gt::ExampleInterface_gt(ExampleInterface_Expression *var1, int constant)
  : ExampleInterface_binop(var1,constant)
{
  std::cout << "Creating a greater than constraint" << std::endl;
}

ExampleInterface_gt::~ExampleInterface_gt()
{
  std::cout << "delete greaterthan" << std::endl;
}

ExampleInterface_Expression* ExampleInterface_gt::add(ExampleInterfaceSolver *solver, bool top_level) {
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

ExampleInterface_Minimise::ExampleInterface_Minimise(ExampleInterface_Expression *var)
  : ExampleInterface_Expression()
{
  std::cout << "Creating a minimise objective" << std::endl;
}

ExampleInterface_Minimise::~ExampleInterface_Minimise()
{
  std::cout << "delete minimise" << std::endl;
}

ExampleInterface_Expression* ExampleInterface_Minimise::add(ExampleInterfaceSolver *solver, bool top_level) {
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

ExampleInterface_Maximise::ExampleInterface_Maximise(ExampleInterface_Expression *var)
  : ExampleInterface_Expression()
{
  std::cout << "Creating a maximise objective" << std::endl;
}

ExampleInterface_Maximise::~ExampleInterface_Maximise()
{
  std::cout << "delete maximise" << std::endl;
}

ExampleInterface_Expression* ExampleInterface_Maximise::add(ExampleInterfaceSolver *solver, bool top_level) {
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

ExampleInterfaceSolver::ExampleInterfaceSolver()
{
  std::cout << "c (wrapper) creating solver" << std::endl;
}

ExampleInterfaceSolver::~ExampleInterfaceSolver()
{
  std::cout << "c (wrapper) delete solver" << std::endl;
}

void ExampleInterfaceSolver::add(ExampleInterface_Expression* arg)
{
  std::cout << "adding expression" << std::endl;
  
  arg->add(this, true);
  
}

void ExampleInterfaceSolver::initialise(ExampleInterfaceExpArray& arg)
{
  std::cout << "Initialising solver with array of expressions" << std::endl;
}

void ExampleInterfaceSolver::initialise()
{
  std::cout << "Initialising solver with no expressions" << std::endl;
}

int ExampleInterfaceSolver::num_vars() {
  std::cout << "return number of variables" << std::endl;
  return 0;
}

int ExampleInterfaceSolver::get_degree(int i) {
  std::cout << "return degree of expression i" << std::endl;
  return 0;
}

int ExampleInterfaceSolver::solveAndRestart(const int policy, 
           const unsigned int base, 
           const double factor,
           const double decay,
           const int reinit)
{
  std::cout << " SOLVE!! at level: ..." << std::endl;
  return 0;
}

int ExampleInterfaceSolver::solve()
{
  std::cout << " SOLVE!! " << std::endl;
  return 0;
}

int ExampleInterfaceSolver::startNewSearch()
{
  std::cout << "starting new search" << std::endl;
  return 0;
}

int ExampleInterfaceSolver::getNextSolution()
{
  std::cout << "getting next solution" << std::endl;
  return 0;
}

int ExampleInterfaceSolver::sacPreprocess(const int type)
{
  std::cout << "running sac preprocessing" << std::endl;
  return 0;
}

bool ExampleInterfaceSolver::propagate()
{
  std::cout << "propagating" << std::endl;
  return 0;
}

int ExampleInterfaceSolver::get_level()
{
  std::cout << "return current level of solver" << std::endl;
  return 0;
}

bool ExampleInterfaceSolver::undo(const int nlevel) {
  std::cout << "Going back nlevel levels" << std::endl;
  return 0;
} 

void ExampleInterfaceSolver::save() 
{
  std::cout << "saving state" << std::endl;
}

int ExampleInterfaceSolver::next(ExampleInterface_Expression* x, int v)
{
  std::cout << "get next" << std::endl;
  return 0;
}

void ExampleInterfaceSolver::post(const char* op, ExampleInterface_Expression* x, int v)
{
  std::cout << "Posting a constraint" << std::endl;
}

void ExampleInterfaceSolver::deduce(const char* op, ExampleInterface_Expression* x, int v) {
  std::cout << "deducing" << std::endl;
}

void ExampleInterfaceSolver::deduce() {
  std::cout << "deducing" << std::endl;
}

bool ExampleInterfaceSolver::branch_right() {
  std::cout << "branching right" << std::endl;
  return 0;
}

void ExampleInterfaceSolver::store_solution() {
  std::cout << "storing solution" << std::endl;
}

void ExampleInterfaceSolver::setHeuristic(const char* var_heuristic, const char* val_heuristic, const int rand)
{
  std::cout << "Setting heuristics" << std::endl;
}

void ExampleInterfaceSolver::addNogood(ExampleInterfaceExpArray& vars, 
            ExampleInterfaceIntArray& vals)
{
  std::cout << "Adding nogood constraint" <<std::endl;
}

void ExampleInterfaceSolver::guide(ExampleInterfaceExpArray& vars, 
        ExampleInterfaceIntArray& vals,
        ExampleInterfaceDoubleArray& probs)
{
  std::cout << "Adding nogood constraint" <<std::endl;
}

void ExampleInterfaceSolver::forceFiniteDomain(ExampleInterfaceExpArray& vars)
{
  std::cout << "Forcing finite domain" <<std::endl;
}

void ExampleInterfaceSolver::backtrackTo(const int level) {
  std::cout << "backtracking to level" <<std::endl;
}

void ExampleInterfaceSolver::presolve()
{
  std::cout << "presolving" <<std::endl;
}

void ExampleInterfaceSolver::increase_init_level(const int i)
{
  std::cout << "increasing initial level" <<std::endl;
}

void ExampleInterfaceSolver::decrease_init_level(const int i)
{
  std::cout << "decreasing initial level" <<std::endl;
}

void ExampleInterfaceSolver::assign(ExampleInterface_Expression *X, const int v) {
  std::cout << "Setting domain of expression X to v" <<std::endl;
}

void ExampleInterfaceSolver::upOneLevel() {
  std::cout << "stepping up one level" <<std::endl;
}

void ExampleInterfaceSolver::reset(bool full) {
  std::cout << "resetting solver" <<std::endl;
}

void ExampleInterfaceSolver::setLowerBounds(ExampleInterfaceExpArray& vars, 
           ExampleInterfaceIntArray& vals) {
  std::cout << "stepping up one level" <<std::endl;
}

void ExampleInterfaceSolver::setUpperBounds(ExampleInterfaceExpArray& vars, 
           ExampleInterfaceIntArray& vals) {
  std::cout << "stetting upper bounds" <<std::endl;
}

void ExampleInterfaceSolver::setRestartNogood() 
{
  std::cout << "setting restart no good" <<std::endl;
}

void ExampleInterfaceSolver::setFailureLimit(const int cutoff)
{
  std::cout << "setting failure limit" <<std::endl;
}

void ExampleInterfaceSolver::setNodeLimit(const int cutoff)
{
  std::cout << "setting node limit" <<std::endl;
}

void ExampleInterfaceSolver::setTimeLimit(const int cutoff)
{
  std::cout << "setting time limit" <<std::endl;
}

void ExampleInterfaceSolver::setVerbosity(const int degree)
{
  std::cout << "setting verbosity" <<std::endl;
}

void ExampleInterfaceSolver::setRandomized(const int degree)
{
  std::cout << "setting randomised" <<std::endl;
}

void ExampleInterfaceSolver::setRandomSeed(const int seed)
{
  std::cout << "setting random seed" <<std::endl;
}

bool ExampleInterfaceSolver::is_opt()
{
  std::cout << "returning is_opt()" <<std::endl;
  return true;
}

bool ExampleInterfaceSolver::is_sat()
{
  std::cout << "returning is satisfied?" <<std::endl;
  return true;
}

bool ExampleInterfaceSolver::is_unsat()
{
  std::cout << "returning is NOT satisfied?" <<std::endl;
  return true;
}

void ExampleInterfaceSolver::printStatistics()
{
  std::cout << "printing statistics" <<std::endl;
}

int ExampleInterfaceSolver::getBacktracks()
{
  std::cout << "return number of backtracks" <<std::endl;
  return 0;
}

int ExampleInterfaceSolver::getNodes()
{
  std::cout << "return number of nodes" <<std::endl;
  return 0;
}

int ExampleInterfaceSolver::getFailures()
{
  std::cout << "return number of failures" <<std::endl;
  return 0;
}

int ExampleInterfaceSolver::getChecks()
{
  std::cout << "return number of checks" <<std::endl;
  return 0;
}

int ExampleInterfaceSolver::getPropags()
{
  std::cout << "return amount of propagation" <<std::endl;
  return 0;
}

double ExampleInterfaceSolver::getTime()
{
  std::cout << "return duration" <<std::endl;
  return 0;
}

int ExampleInterfaceSolver::getNumVariables()
{
  std::cout << "return number of variables" <<std::endl;
  return 0;
}

int ExampleInterfaceSolver::getNumConstraints()
{
  std::cout << "return number of constraints" <<std::endl;
  return 0;
}

int ExampleInterfaceSolver::getRandomNumber()
{
  return 8;
}
