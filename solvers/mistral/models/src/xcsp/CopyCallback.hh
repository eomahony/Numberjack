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
#ifndef _CopyCallback_h_
#define _CopyCallback_h_

#include "CSPParserCallback.hh"

namespace CSPXMLParser
{
  using namespace std;

/**
 * this callback only only tries to regenerate an XML instance (it
 * generates a Copy)
 *
 * THE IMPLEMENTATION IS CURRENTLY INCOMPLETE !!!
 *
 * it's a sample callback that may be used to test the parser
 */
class CopyCallback : public CSPParserCallback
{
private:
  bool firstDomainValue;
  bool firstTuple;

public:
  virtual void beginInstance(const string & name)
  {
    cout << "<instance>" << endl;
    cout << "<presentation name='"
	 << name << "'/>" << endl;
  }

  virtual void beginDomainsSection(int nbDomains) 
  {  
    cout << "<domains nbDomains='" 
	 << nbDomains << "'>" <<endl;
  }
  
  virtual void beginDomain(const string & name, int idDomain, int nbValue) 
  {
    cout << " <domain name='" << name 
	 << "' nbValues='" << nbValue << "'>";
    firstDomainValue=true;
  }

  virtual void addDomainValue(int v) 
  {
    if (!firstDomainValue)
      cout << ' ';
    firstDomainValue=false;

    cout << v;
  }

  virtual void addDomainValue(int first,int last) 
  {
    if (!firstDomainValue)
      cout << ' ';
    firstDomainValue=false;

    cout << first << ".." << last;
  }

  virtual void endDomain() 
  {
    cout << "</domain>" <<endl;
  }

  /**
   * end the definition of all domains
   */
  virtual void endDomainsSection() 
  {
    cout << "</domains>" <<endl;
  }


  virtual void beginVariablesSection(int nbVariables) 
  {
    cout << "<variables nbVariables='"
	 << nbVariables << "'>" <<endl;
  }
  
  virtual void addVariable(const string & name, int idVar,
			   const string & domain, int idDomain) 
  {
    cout << "<variable name='" << name
	 << "' domain='" << domain 
	 << "'/>" <<endl;
  }

  virtual void endVariablesSection() 
  {
    cout << "</variables>" <<endl;
  }


  virtual void beginRelationsSection(int nbRelations) 
  {
    cout << "<relations nbRelations='" << nbRelations 
	 << "'>" <<endl;
  }

  
  virtual void beginRelation(const string & name, int idRel,
			     int arity, int nbTuples, RelType relType) 
  {
    cout << "<relation name='" << name 
	 << "' arity='" << arity << "' nbTuples='"
	 << nbTuples << "' semantics='";

    switch(relType)
    {
    case REL_SUPPORT:
      cout << "supports";
      break;
    case REL_CONFLICT:
      cout << "conflicts";
      break;
    case REL_SOFT:
      cout << "soft";
      break;
    default:
      throw runtime_error("unknown relation type");
    }

    cout << "'>";
    firstTuple=true;
  }

  virtual void addRelationTuple(int arity, int tuple[]) 
  {      
    if (!firstTuple)
      cout << '|';

    firstTuple=false;

    cout << tuple[0];
    for(int i=1;i<arity;++i)
      cout << ' ' << tuple[i];
  }

  virtual void addRelationTuple(int arity, int tuple[], int cost) 
  {      
    if (!firstTuple)
      cout << '|';

    firstTuple=false;

    cout << cost << ":" << tuple[0];
    for(int i=1;i<arity;++i)
      cout << ' ' << tuple[i];
  }

  virtual void endRelation() 
  {
    cout << "</relation>" <<endl;
  }

  virtual void endRelationsSection() 
  {
    cout << "</relations>" <<endl;
  }

  virtual void beginPredicatesSection(int nbPredicates) 
  {
    cout << "<predicates nbPredicates='" 
	 << nbPredicates << "'>" << endl;
  }
  
  virtual void beginPredicate(const string & name, int idPred) 
  {
    cout << "<predicate name='" << name 
	 << "'>" << endl;
  }

  virtual void addFormalParameter(int pos, const string & name, 
				  const string & type) 
  {
    cout << "   formal parameter " << pos << ": "
	 << type << " " << name << endl;
  }

  virtual void predicateExpression(AST *tree) 
  {
    cout << "   predicate definition (AST) = ";
    tree->prefixExpression(cout);
    cout << endl;
  }

  virtual void predicateExpression(const string &expr) 
  {
    cout << "   predicate definition=" << expr << endl;
  }

  virtual void endPredicate() 
  {
    cout << "</predicate>" << endl;
  }

  virtual void endPredicatesSection() 
  {
    cout << "</predicates>" << endl;
  }

  virtual void beginConstraintsSection(int nbConstraints) 
  {
    cout << "<constraints nbConstraints='"
	 << nbConstraints << "'>" <<endl;
  }
  
  virtual void beginConstraint(const string & name, int idConstr,
			       int arity, 
			       const string & reference, 
			       CSPDefinitionType type, int id,
			       const ASTList &scope)
  {
    cout << "<constraint name='" << name
	 << "' arity='" << arity 
	 << "' scope='";

    if (scope.size())
      cout << scope[0].getVarName();

    for(int i=1;i<scope.size();++i)
      cout << ' ' << scope[i].getVarName();

    cout << "' reference='" << reference 
	 << "'>" << endl;
  }

  virtual void constraintParameters(const ASTList &args)
  {
    cout << "constraint parameters=";
    args.postfixExpression(cout);
    cout << endl;

    //cout << "duration=";
    //args[0].list()[0].dict()["duration"].postfixExpression(cout);
    //args[0][1]["duration"].postfixExpression(cout);
  } 


  virtual void endConstraint() 
  {
    cout << "  </constraint>" << endl;
  }

  /**
   * end the definition of all constraints
   */
  virtual void endConstraintsSection() 
  {
    cout << "</constraints>" <<endl;
  }

  /********************************************************************/


  /**
   * signal the end of parsing
   */
  virtual void endInstance() 
  {
    cout << "</instance>" <<endl;
  }

};

} // namespace

#endif

