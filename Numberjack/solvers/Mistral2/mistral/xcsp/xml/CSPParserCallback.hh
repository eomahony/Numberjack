/*=============================================================================
 * parser for CSP instances represented in XML format
 * 
 * Copyright (c) 2008 Olivier ROUSSEL (olivier.roussel <at> cril.univ-artois.fr)
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *=============================================================================
 */
#ifndef _CSPParserCallback_h_
#define _CSPParserCallback_h_

#include <string>
#include "XMLParser_constants.h"

namespace CSPXMLParser
{

using namespace std;

enum CSPDefinitionType 
  {RelationType,PredicateType,GlobalConstraintType, // don't change these ones
   DomainType,VariableType,ConstraintType,
   UndefinedType}; // UndefinedType must come last

class AST; // Abstract Syntax Tree representing an expression

class CSPParserCallback
{
public:
  virtual ~CSPParserCallback() {}

  /**
   * signal the beginning of parsing
   *
   * @param name: name of the instance
   */
  virtual void beginInstance(const string & name) {}

  /********************************************************************/

  /**
   * callback called at the beginning of the domains declarations
   *
   * @param nbDomains: number of domains that will be declared
   */
  virtual void beginDomainsSection(int nbDomains) {}
  
  /**
   * callback called at the beginning of the declaration of one domain
   *
   * @param name: identifier of the domain
   * @param idDomain: identifier assigned to the domain name (starting from 0)
   * @param nbValue: number of values in the domain
   */
  virtual void beginDomain(const string & name, int idDomain, int nbValue) {}

  /**
   * add a single value to the current domain
   *
   * @param v: value to add to the domain
   */
  virtual void addDomainValue(int v) {}

  /**
   * add the range of values [first..last] to the current domain
   *
   * @param first: first value to add to the domain
   * @param last: last value to add to the domain
   */
  virtual void addDomainValue(int first,int last) {}

  /**
   * ends the definition of the current domain
   */
  virtual void endDomain() {}

  /**
   * end the definition of all domains
   */
  virtual void endDomainsSection() {}

  /********************************************************************/

  /**
   * callback called at the beginning of the variables declarations
   *
   * @param nbVariables: number of variables that will be declared
   */
  virtual void beginVariablesSection(int nbVariables) {}
  
  /**
   * callback called to define a new variable
   *
   * @param name: identifier of the variable
   * @param idVar: identifier assigned to the variable name (starting from 0)
   * @param domain: identifier of the variable domain
   * @param idDomain: identifier assigned to the domain name (starting from 0)
   */
  virtual void addVariable(const string & name, int idVar, 
			   const string & domain, int idDomain) {}

  /**
   * end the definition of all variables
   */
  virtual void endVariablesSection() {}


  /********************************************************************/

  /**
   * callback called at the beginning of the relations declarations
   *
   * @param nbRelations: number of relations that will be declared
   */
  virtual void beginRelationsSection(int nbRelations) {}
  
  /**
   * callback called at the beginning of the declaration of one relation
   *
   * @param name: identifier of the relation
   * @param idRel: identifier assigned to the relation name (starting from 0)
   * @param arity: arity of the relation
   * @param nbTuples: number of tuples in the relation
   * @param isSupport: true if tuples represent support, false if
   *                  tuples represent conflicts
   */
  virtual void beginRelation(const string & name, int idRel,
			     int arity, int nbTuples, RelType relType) {}

  /**
   * add a single tuple to the current relation
   *
   * @param arity: the tuple arity
   * @param tuple: tuple to add to the relation (contains arity elements)
   */
  virtual void addRelationTuple(int arity, int tuple[]) {}

  /**
   * add a single weighted tuple to the current relation
   *
   * @param arity: the tuple arity
   * @param tuple: tuple to add to the relation (contains arity elements)
   * @param cost: the cost of this tuple
   */
  virtual void addRelationTuple(int arity, int tuple[], int cost) {}

  /**
   * ends the definition of the current relation
   */
  virtual void endRelation() {}

  /**
   * end the definition of all relations
   */
  virtual void endRelationsSection() {}

  /********************************************************************/

  /**
   * callback called at the beginning of the predicates declarations
   *
   * @param nbPredicates: number of predicates that will be declared
   */
  virtual void beginPredicatesSection(int nbPredicates) {}
  
  /**
   * callback called at the beginning of the declaration of one predicate
   *
   * @param name: identifier of the predicate
   * @param idPred: identifier assigned to the predicate name (starting from 0)
   */
  virtual void beginPredicate(const string & name, int idPred) {}

  /**
   * add a formal parameter to the current predicate
   *
   * @param pos: position of the formal parameter (0=first)
   * @param name: name of the parameter
   * @param type: type of the parameter
   */
  virtual void addFormalParameter(int pos, const string & name, const string & type) {}

  /**
   * provide the expression of the current predicate
   *
   * @param tree: the abstract syntax tree representing the expression
   */
  virtual void predicateExpression(AST *tree) {}


  /**
   * provide the expression of the current predicate
   *
   * @param expr: the string representing the expression
   */
  virtual void predicateExpression(const string &expr) {}

  /**
   * ends the definition of the current predicate
   */
  virtual void endPredicate() {}

  /**
   * end the definition of all predicates
   */
  virtual void endPredicatesSection() {}

  /********************************************************************/

  /**
   * callback called at the beginning of the constraints declarations
   *
   * @param nbConstraints: number of constraints that will be declared
   */
  virtual void beginConstraintsSection(int nbConstraints) {}
  
  /**
   * callback called at the beginning of the declaration of one constraint
   *
   * @param name: identifier of the constraint
   * @param idConstr: identifier assigned to the constraint name (starting from 0)
   * @param arity: arity of the constraint
   * @param name: the refererence to the definition of this constraint. May be a relation, a predicate or the name of a global constraint
   * @param reference: the name of the relation/predicate or global constraint which defines the support/conflict tuples of this constraint
   * @param type: type of reference (RelationType,PredicateType or
   *             GlobalConstraintType)
   * @param id: identifier associated to the reference
   */
  virtual void beginConstraint(const string & name, int idConstr,
			       int arity, 
			       const string & reference, 
			       CSPDefinitionType type, int id,
			       const ASTList &scope) {}

  /**
   * provides the list of parameters of the constraint
   *
   * @param args: the list of effective parameters of the constraint
   */
  virtual void constraintParameters(const ASTList &args) {}

  /**
   * ends the definition of the current constraint
   */
  virtual void endConstraint() {}

  /**
   * end the definition of all constraints
   */
  virtual void endConstraintsSection() {}


  /********************************************************************/


  /**
   * signal the end of parsing
   */
  virtual void endInstance() {}
};

} // namespace
#endif
