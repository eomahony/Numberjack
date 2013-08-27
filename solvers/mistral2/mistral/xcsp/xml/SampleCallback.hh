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
#ifndef _SampleCallback_h_
#define _SampleCallback_h_

#include "CSPParserCallback.hh"

#define ShowNumericalIDS

namespace CSPXMLParser
{
  using namespace std;

// definition of the functions which are called by the parser when it
// reads new information in the file. These functions are in charge of
// creating the data structure that the solver needs to do his job
class SampleCallback : public CSPParserCallback
{
public:
  virtual void beginInstance(const string & name)
  {
    cout << "callback: <instance>" << endl;
    cout << "callback: <presentation name='"
	 << name << "'/>" << endl;
  }

  virtual void beginDomainsSection(int nbDomains) 
  {  
    cout << "callback: <domains nbDomains='" 
	 << nbDomains << "'>" <<endl;
  }
  
  virtual void beginDomain(const string & name, int idDomain, int nbValue) 
  {
    cout << "callback:   <domain name='" << name 
#ifdef ShowNumericalIDS
	 << '/' << idDomain
#endif
	 << "' nbValue='" << nbValue << "'>" << endl;
  }

  void addDomainValue(int v) 
  {
    cout << "callback:     " << v << endl;
  }

  virtual void addDomainValue(int first,int last) 
  {
    cout << "callback:     " << first << ".." << last << endl;
  }

  virtual void endDomain() 
  {
    cout << "callback:   </domain>" <<endl;
  }

  /**
   * end the definition of all domains
   */
  virtual void endDomainsSection() 
  {
    cout << "callback: </domains>" <<endl;
  }


  virtual void beginVariablesSection(int nbVariables) 
  {
    cout << "callback: <variables nbVariables='"
	 << nbVariables << "'>" <<endl;
  }
  
  virtual void addVariable(const string & name, int idVar,
			   const string & domain, int idDomain) 
  {
    cout << "callback:   <variable name='" << name
#ifdef ShowNumericalIDS
	 << '/' << idVar
#endif
	 << "' domain='" << domain 
#ifdef ShowNumericalIDS
	 << '/' << idDomain
#endif
	 << "'/>" <<endl;
  }

  virtual void endVariablesSection() 
  {
    cout << "callback: </variables>" <<endl;
  }


  virtual void beginRelationsSection(int nbRelations) 
  {
    cout << "callback: <relations nbRelations='" << nbRelations 
	 << "'>" <<endl;
  }

  
  virtual void beginRelation(const string & name, int idRel,
			     int arity, int nbTuples, bool isSupport) 
  {
    cout << "callback:   <relation name='" << name 
#ifdef ShowNumericalIDS
	 << '/' << idRel
#endif
	 << "' arity='" << arity << "' nbTuples='"
	 << nbTuples << "' semantics='";

    if (isSupport)
      cout << "supports";
    else
      cout << "conflicts";

    cout << "'>" << endl;
  }

  virtual void addRelationTuple(int arity, int tuple[]) 
  {
    cout << "callback:     ";
    for(int i=0;i<arity;++i)
      cout << tuple[i] << " ";
    cout << endl;
  }

  virtual void endRelation() 
  {
    cout << "callback:   </relation>" <<endl;
  }

  virtual void endRelationsSection() 
  {
    cout << "callback: </relations>" <<endl;
  }

  virtual void beginPredicatesSection(int nbPredicates) 
  {
    cout << "callback: <predicates nbPredicates='" 
	 << nbPredicates << "'>" << endl;
  }
  
  virtual void beginPredicate(const string & name, int idPred) 
  {
    cout << "callback:   <predicate name='" << name 
#ifdef ShowNumericalIDS
	 << '/' << idPred
#endif
	 << "'>" << endl;
  }

  virtual void addFormalParameter(int pos, const string & name, 
				  const string & type) 
  {
    cout << "callback:    formal parameter " << pos << ": "
	 << type << " " << name << endl;
  }

  virtual void predicateExpression(AST *tree) 
  {
    cout << "callback:    predicate definition (AST) = ";
    tree->prefixExpression(cout);
    cout << endl;
  }

  virtual void predicateExpression(const string &expr) 
  {
    cout << "callback:    predicate definition=" << expr << endl;
  }

  virtual void endPredicate() 
  {
    cout << "callback: </predicate>" << endl;
  }

  virtual void endPredicatesSection() 
  {
    cout << "callback: </predicates>" << endl;
  }

  virtual void beginConstraintsSection(int nbConstraints) 
  {
    cout << "callback: <constraints nbConstraints='"
	 << nbConstraints << "'>" <<endl;
  }
  
  virtual void beginConstraint(const string & name, int idConstr,
			       int arity, 
			       const string & reference, 
			       CSPDefinitionType type, int id)
  {
    cout << "callback:   <constraint name='" << name
#ifdef ShowNumericalIDS
	 << '/' << idConstr
#endif
	 << "' arity='" << arity 
	 << "' reference='" << reference 
#ifdef ShowNumericalIDS
	 << "/type=" << type << "/id=" << id
#endif
	 << "'>" << endl;
  }

  virtual void addVariableToConstraint(const string & name, int idVar) 
  {
    cout << "callback:     scope -> " << name 
#ifdef ShowNumericalIDS
	 << '/' << idVar
#endif
	 << endl;
  }

  virtual void addEffectiveParameter(int pos, const string & expr) 
  {
    cout << "callback: effective parameter " << pos << " = "
	 << expr << endl;
  }

  virtual void addEffectiveParameter(int pos, AST *tree) 
  {
    cout << "callback: effective parameter AST" << endl;
  }

  virtual void addEffectiveParameterVariable(int pos, 
					     const string & name, 
					     int idVar) 
  {
    cout << "callback:     effective parameter " << pos << " = variable "
	 << name
#ifdef ShowNumericalIDS
	 << '/' << idVar
#endif
	 << endl;
  }

  virtual void addEffectiveParameterInteger(int pos, 
					    int value)
  {
    cout << "callback:     effective parameter " << pos << " = integer "
	 << value << endl;
  }

  virtual void addEffectiveParameterList(int pos) 
  {
    cout << "callback:     effective parameter " << pos 
	 << " = list (see below) "
	 << endl;
  }

  /**
   * indicate the start of a list
   */
  virtual void beginList() 
  {
    cout << "callback:     <list>" << endl;
  }

  virtual void addListVariable(int pos, const string & name, int idVar) 
  {
    cout << "callback:       element " << pos << " = variable "
	 << name
#ifdef ShowNumericalIDS
	 << '/' << idVar
#endif
	 << endl;
  }

  virtual void addListInteger(int pos, int value) 
  {
    cout << "callback:       element " << pos << " = integer "
	 << value << endl;
  }

  virtual void endList() 
  {
    cout << "callback:     </list>" << endl;
  }

  virtual void endConstraint() 
  {
    cout << "callback:   </constraint>" << endl;
  }

  /**
   * end the definition of all constraints
   */
  virtual void endConstraintsSection() 
  {
    cout << "callback: </constraints>" <<endl;
  }

  /********************************************************************/


  /**
   * signal the end of parsing
   */
  virtual void endInstance() 
  {
    cout << "callback: </instance>" <<endl;
  }

};

} // namespace

#endif

