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
#ifndef _XMLParser_h_
#define _XMLParser_h_

/**
 * @file XMLParser.h
 * @brief defines the XMLParser class to parse a CSP instance in XML format
 */

/*
 * TODO
 *
 * - check that arity and scope of constraints correspond to the parameters
 *
 * - handle the formal parameters via the operandStack
 *
 * - detect unused attributes
 * - add the necessary syntax checks
 * - check the number of elements per section
 * - check arity and type of formal/effective parameters
 * - be more explicit on errors
 */


#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <deque>
#include <vector>
#include <sstream>
#include <stdexcept>
#include <cassert>

#include "ExpressionParser.hh"
#include "CSPParserCallback.hh"

/**
 * @namespace CSPXMLParser
 * @brief this namespace encloses all definitions relative to the
 * CSP XML format parser.
 */
namespace CSPXMLParser
{

  using namespace std;

  /**
   * @brief contains a parser for the CSP XML format.
   *
   * This class proposes a parser in SAX mode (to use less memory than
   * in DOM mode)
   *
   * @param Callback the type of the class that will be used to transmit
   * information from the parser to the solver
   *
   * @param XMLLibraryString a class that encapsulates the strings
   * used by the XML library
   *
   * @param XMLLibraryAttributeList a class that encapsulates the list
   * of attributes of a tag
   *
   */
  template<class Callback, 
	   class XMLLibraryString, 
	   class XMLLibraryAttributeList>
  class XMLParser
  {
  public:
    // list of attributes and values for a tag
    typedef XMLLibraryAttributeList AttributeList;
    typedef XMLLibraryString XMLString;

  private:
    Callback *cb;
    Syntax preferredSyntax; // preferred syntax for
    // transmitting an expression to the solver

    ExpressionParser exprParser;

    /// stores some information on each symbol defined in the CSP instance
    struct SymbolInfo
    {
      CSPDefinitionType type; ///< kind of definition
      int id; ///< identifier

      // temporary fix ???
      vector<string> *convOrder;

      SymbolInfo()
      {
	id=-1;
	type=UndefinedType;
	convOrder=NULL;
      }
      
      SymbolInfo(int id, CSPDefinitionType type)
      {
	this->id=id;
	this->type=type;
	convOrder=NULL;
      }

      SymbolInfo(int id, CSPDefinitionType type, vector<string> *convOrder)
      {
	this->id=id;
	this->type=type;
	this->convOrder=convOrder;
      }
    };

    typedef map<string,SymbolInfo> SymbolDictionnary;

    SymbolDictionnary symbolDict;

    // stores the next available ID in each category
    int nextId[UndefinedType];

    /**
     * register a new symbol in its category.
     *
     * @return the identifier assigned to this symbol
     */
    int defineSymbol(const string &name, CSPDefinitionType type)
    {
      if (type==UndefinedType)
	throw runtime_error("Cannot register a symbol without a type");

      if (type<0 || type>UndefinedType)
	throw runtime_error("internal error: illegal CSPDefinitionType");

      SymbolInfo &info=symbolDict[name];

      if (info.type!=UndefinedType)
	throw runtime_error("symbol "+name+" is defined twice");

      info.type=type;
      info.id=nextId[type]++;

      return info.id;
    }


    /**
     * get the identifier assigned to an already registered symbol
     *
     * @param name symbol to look up
     * @param expectedType if different from UndefinedType, check that
     * the symbol has the expected type given by this parameter
     * @return the identifier assigned to this symbol
     */
    int getSymbolId(const string &name, 
		    CSPDefinitionType expectedType=UndefinedType)
    {
      SymbolInfo &info=symbolDict[name];

      if (info.type==UndefinedType)
	throw runtime_error("symbol "+name+" is undefined");

      if (expectedType!=UndefinedType && info.type!=expectedType)
	throw runtime_error("symbol "+name+" doesn't have the expected type");

      return info.id;
    }


    /**
     * get the information about an already registered symbol
     *
     * @param name symbol to look up
     * @return type and identifier associated to this symbol
     */
    SymbolInfo getSymbolInfo(const string &name)
    {
      typename SymbolDictionnary::iterator it=symbolDict.find(name);

      if (it==symbolDict.end())
	throw runtime_error("symbol "+name+" is undefined");

      return (*it).second;
    }

    /**
     * transmit an expression to the solver
     *
     */
    void xmitExpressionToSolver(AST *tree)
    {
      if (preferredSyntax==TREE)
      {
	getCallback()->
	  predicateExpression(tree);
      }
      else
      {
	ostringstream tmp;

	tree->expression(tmp,preferredSyntax);

	getCallback()->
	  predicateExpression(tmp.str());

	delete tree;
      }
    }


    /**
     * @brief the action to be performed (callback) when a tag is read in the
     * XML file
     */
    class TagAction
    {
    protected:
      bool activated;
      XMLParser *parser;

    public:
      TagAction(XMLParser *parser): parser(parser) {activated=false;}
      virtual ~TagAction() {}

      void activate() {activated=true;}
      void deactivate() {activated=false;}

      bool isActivated() {return activated;}

      /**
       * return the name of this tag (constant)
       */
      virtual const char *getTagName()=0;
      virtual void beginTag(const AttributeList &attributes)=0;

      /**
       *
       * @parm last must be true if this is the last chunk of text
       */
      virtual void text(const XMLString txt, bool last)=0;
      virtual void endTag()=0;

    protected:
      /**
       * check that the parent tag in the XML file has the indicated name
       *
       * when name is NULL, check that this tag has no parent in the XML file
       *
       * @param parentTag: the expected name of the parent tag or NULL to
       * check that we have no parent
       * @param n: 1 for parent, 2 for grand-parent and so on
       */
      void checkParentTag(const char *parentTag, int n=1)
      {
	if (this->parser->getParentTagAction(n)==NULL)
	{
	  if (parentTag!=NULL)
	    throw runtime_error("tag has no parent but it should have one");
	  else
	    return;
	}

	if (parentTag==NULL)
	  throw runtime_error("tag should have no parent");

#ifdef debug
	cout << "parent tag=" 
	     << this->parser->getParentTagAction(n)->getTagName() << endl;
#endif
	if (strcmp(this->parser->getParentTagAction(n)->getTagName(),
		   parentTag)!=0)
	  throw runtime_error("wrong parent for tag");
      }
    };

    typedef 
    map<XMLString,TagAction *> TagActionList;

    TagActionList tagList;

    struct State
    {
      // may we used abridged notations
      bool abridgedListAllowed,abridgedDictAllowed,
	abridgedVariableAllowed,abridgedIntegerAllowed;

      bool subtagAllowed;

      State()
      {
	subtagAllowed=true;

	setAllAbridgedAllowed(false);
      }

      void setAllAbridgedAllowed(bool v) 
      {
	abridgedListAllowed=abridgedDictAllowed=
	  abridgedVariableAllowed=abridgedIntegerAllowed=v;
      }
    };

    // the top of the following stacks is at the FRONT
    deque<TagAction *> actionStack;
    deque<State> stateStack;


    void registerTagAction(TagActionList &tagList, TagAction *action)
    {
      tagList[action->getTagName()]=action;
    }

    /**
     * a handler to silently ignore unkown tags
     */
    class UnknownTagAction : public TagAction
    {
    public:
      UnknownTagAction(XMLParser *parser): TagAction(parser) {}

      virtual const char *getTagName() {return "";}

      virtual void beginTag(const AttributeList &attributes)
      {
	// ignore
      }

      virtual void text(const XMLString txt, bool last)
      {
	// ignore
      }

      virtual void endTag()
      {
	// ignore
      }
    };

    class InstanceTagAction : public TagAction
    {
    public:
      InstanceTagAction(XMLParser *parser): TagAction(parser) {}

      virtual const char *getTagName() {return "instance";}

      virtual void beginTag(const AttributeList &attributes)
      {
      }

      virtual void text(const XMLString txt, bool last)
      {
	if (!txt.isWhiteSpace())
	  throw runtime_error("<instance> tag should not have meaningful text");
      }

      virtual void endTag()
      {
	this->parser->getCallback()->endInstance();
      }
    };

    class PresentationTagAction : public TagAction
    {
    public:
      PresentationTagAction(XMLParser *parser): TagAction(parser) {}

      virtual const char *getTagName() {return "presentation";}

      virtual void beginTag(const AttributeList &attributes)
      {
	this->parser->stateStack.front().setAllAbridgedAllowed(false);
	this->parser->stateStack.front().subtagAllowed=false;

	string name;

	this->checkParentTag("instance");

	if (!attributes["name"].to(name))
	  name="?";

	this->parser->getCallback()->beginInstance(name);
      }

      virtual void text(const XMLString txt, bool last)
      {
	// ignore
      }

      virtual void endTag()
      {
      }
    };

    class DomainsTagAction : public TagAction
    {
    public:
      DomainsTagAction(XMLParser *parser): TagAction(parser) {}

      virtual const char *getTagName() {return "domains";}

      virtual void beginTag(const AttributeList &attributes)
      {
	int nbDomains;

	this->checkParentTag("instance");

	if (!attributes["nbDomains"].to(nbDomains))
	  throw runtime_error("expected attribute nbDomains for tag "
			      "<domains>");

	this->parser->getCallback()->beginDomainsSection(nbDomains);
      }

      virtual void text(const XMLString txt, bool last)
      {
	if (!txt.isWhiteSpace())
	  throw runtime_error("<domains> tag should not have meaningful text");
      }

      virtual void endTag()
      {
	this->parser->getCallback()->endDomainsSection();
      }
    };

    class DomainTagAction : public TagAction
    {
    private:
      XMLString dotdot;

    public:
      DomainTagAction(XMLParser *parser): TagAction(parser), dotdot("..") 
      {
	this->parser->stateStack.front().setAllAbridgedAllowed(false);
	this->parser->stateStack.front().subtagAllowed=false;
      }

      virtual const char *getTagName() {return "domain";}

      virtual void beginTag(const AttributeList &attributes)
      {
	string name;
	int nbValues;

	this->checkParentTag("domains");

	if (!attributes["name"].to(name))
	  throw runtime_error("expected attribute name for tag "
			      "<domain>");

	if (!attributes["nbValues"].to(nbValues))
	  throw runtime_error("expected attribute nbValue for tag "
			      "<domain>");

	int idDomain=this->parser->defineSymbol(name,DomainType);

	this->parser->getCallback()->beginDomain(name,idDomain,nbValues);
      }

      virtual void text(const XMLString txt, bool last)
      {
	typename XMLString::Tokenizer tokenizer(txt);

	while(tokenizer.hasMoreTokens())
	{
	  XMLString token=tokenizer.nextToken();

	  int pos=token.find(dotdot);

	  if (pos==XMLString::npos)
	  {
	    int val;
	    token.to(val);
	    this->parser->getCallback()->addDomainValue(val);
	  }
	  else
	  {
	    int first,last;

	    token.substr(0,pos).to(first);
	    token.substr(pos+2).to(last);

	    this->parser->getCallback()->addDomainValue(first,last);
	  }
	}
      }

      virtual void endTag()
      {
	this->parser->getCallback()->endDomain();
      }
    };

    class VariablesTagAction : public TagAction
    {
    public:
      VariablesTagAction(XMLParser *parser): TagAction(parser) {}

      virtual const char *getTagName() {return "variables";}

      virtual void beginTag(const AttributeList &attributes)
      {
	int nbVariables;

	this->checkParentTag("instance");

	if (!attributes["nbVariables"].to(nbVariables))
	  throw runtime_error("expected attribute nbVariables for tag "
			      "<variables>");

	this->parser->getCallback()->beginVariablesSection(nbVariables);
      }

      virtual void text(const XMLString txt, bool last)
      {
	if (!txt.isWhiteSpace())
	  throw runtime_error("<variables> tag should not have meaningful text");
      }

      virtual void endTag()
      {
	this->parser->getCallback()->endVariablesSection();
      }
    };

    class VariableTagAction : public TagAction
    {
    public:
      VariableTagAction(XMLParser *parser): TagAction(parser) {}

      virtual const char *getTagName() {return "variable";}

      virtual void beginTag(const AttributeList &attributes)
      {
	this->parser->stateStack.front().setAllAbridgedAllowed(false);
	this->parser->stateStack.front().subtagAllowed=false;

	string name;
	string domain;

	this->checkParentTag("variables");

	if (!attributes["name"].to(name))
	  throw runtime_error("expected attribute name for tag "
			      "<variable>");

	if (!attributes["domain"].to(domain))
	  throw runtime_error("expected attribute domain for tag "
			      "<variable>");

	int idVar=this->parser->defineSymbol(name,VariableType);
	int idDomain=this->parser->getSymbolId(domain,DomainType);

	this->parser->getCallback()->addVariable(name,idVar,domain,idDomain);
      }

      virtual void text(const XMLString txt, bool last)
      {
	if (!txt.isWhiteSpace())
	  throw runtime_error("<variable> tag should not have meaningful text");
      }

      virtual void endTag()
      {
      }
    };

    class RelationsTagAction : public TagAction
    {
    public:
      RelationsTagAction(XMLParser *parser): TagAction(parser) {}

      virtual const char *getTagName() {return "relations";}

      virtual void beginTag(const AttributeList &attributes)
      {
	int nbRelations;

	this->checkParentTag("instance");

	if (!attributes["nbRelations"].to(nbRelations))
	  throw runtime_error("expected attribute nbRelation for tag "
			      "<relations>");

	this->parser->getCallback()->beginRelationsSection(nbRelations);
      }

      virtual void text(const XMLString txt, bool last)
      {
	if (!txt.isWhiteSpace())
	  throw runtime_error("<relations> tag should not have meaningful text");
      }

      virtual void endTag()
      {
	this->parser->getCallback()->endRelationsSection();
      }
    };

    class RelationTagAction : public TagAction
    {
    public:
      RelationTagAction(XMLParser *parser): TagAction(parser) 
      {
	tuple=NULL;
	i=0;
	nb=0;
      }

      virtual ~RelationTagAction()
      {
	if (tuple)
	  delete[] tuple;
      }

      virtual const char *getTagName() {return "relation";}

      virtual void beginTag(const AttributeList &attributes)
      {
	this->parser->stateStack.front().setAllAbridgedAllowed(false);
	this->parser->stateStack.front().subtagAllowed=false;

	string name,semantics;

	this->checkParentTag("relations");

	if (!attributes["arity"].to(arity))
	  throw runtime_error("expected attribute arity for tag "
			      "<relation>");

	if (arity<=0)
	  throw runtime_error("arity must be >0");

	if (!attributes["nbTuples"].to(nbTuples))
	  throw runtime_error("expected attribute nbTuples for tag "
			      "<relation>");

	if (!attributes["name"].to(name))
	  throw runtime_error("expected attribute name for tag "
			      "<relation>");

	if (!attributes["semantics"].to(semantics))
	  throw runtime_error("expected attribute semantics for tag "
			      "<relation>");

	RelType relType;

	isSoft=false;
	defaultCost=1;

	if (semantics=="supports")
	  relType=REL_SUPPORT;
	else
	  if (semantics=="conflicts")
	    relType=REL_CONFLICT;
	  else
	    if (semantics=="soft")
	    {
	      isSoft=true;
	      relType=REL_SOFT;
	    }
	    else
	      throw runtime_error("unknown semantics for tag <relation>");

	if (isSoft)
	{
	  if (!attributes["defaultCost"].to(defaultCost))
	    throw runtime_error("expected attribute defaultCost for tag "
				"<relation>");
	}


	int idRel=this->parser->defineSymbol(name,RelationType);
	this->parser->getCallback()->beginRelation(name,idRel,
						   arity,nbTuples,relType);

	if (tuple)
	  delete[] tuple;

	tuple=new int[arity];

	i=0;
	nb=0;
	currentCost=defaultCost;
      }

      virtual void text(const XMLString txt, bool last)
      {
	typename XMLString::Tokenizer tokenizer(txt);

	tokenizer.addSeparator('|');
	tokenizer.addSeparator(':');

	while(tokenizer.hasMoreTokens())
	{
	  XMLString token=tokenizer.nextToken();

	  //cout << "token=\"" << token << "\"" << endl;

	  if (token==XMLString(":"))
	  {
	    if (isSoft && i==1)
	    {
	      currentCost=tuple[0];
	      i=0;
	    }
	    else 
	      throw runtime_error("unexpected colon (':') in relation");
	  }
	  else
	    if (token==XMLString("|"))
	    {
	      if (i==arity)
		sendTuple();
	      else
		throw runtime_error("unexpected pipe ('|') in relation");
	    }
	    else
	    {
	      if (!token.to(tuple[i++]))
		throw runtime_error("failed to read integer");
	    }
	  
	} // while
      }

      virtual void endTag()
      {
	if (i==arity)
	  sendTuple();

	if (i!=0)
	  throw runtime_error("last tuple is incomplete !");

	if (nb!=nbTuples)
	  throw runtime_error("effective number of tuples doesn't match the announced number");

	this->parser->getCallback()->endRelation();
      }

    protected:
      void sendTuple()
      {
	if (isSoft)
	  this->parser->getCallback()->addRelationTuple(arity,tuple,
							currentCost);
	else
	  this->parser->getCallback()->addRelationTuple(arity,tuple);

	nb++;
	i=0;
      }

    protected:
      int arity;
      int nbTuples;

      
      /*
       * the text() method may not get the whole text in one single call
       * so we must be extra careful
       */

      int *tuple; // the tuple array
      int i; // position in the tuple of the next value 
      int nb; // number of tuples read so far

      bool isSoft;
      int currentCost,defaultCost;

    };

    class PredicatesTagAction : public TagAction
    {
    public:
      PredicatesTagAction(XMLParser *parser): TagAction(parser) {}

      virtual const char *getTagName() {return "predicates";}

      virtual void beginTag(const AttributeList &attributes)
      {
	int nbPredicates;

	this->checkParentTag("instance");

	if (!attributes["nbPredicates"].to(nbPredicates))
	  throw runtime_error("expected attribute nbPredicate for tag "
			      "<predicates>");

	this->parser->getCallback()->beginPredicatesSection(nbPredicates);
      }

      virtual void text(const XMLString txt, bool last)
      {
	if (!txt.isWhiteSpace())
	  throw runtime_error("<predicates> tag should not have meaningful text");
      }

      virtual void endTag()
      {
	this->parser->getCallback()->endPredicatesSection();
      }
    };

    class PredicateTagAction : public TagAction
    {
    public:
      PredicateTagAction(XMLParser *parser): TagAction(parser) {}

      virtual const char *getTagName() {return "predicate";}

      virtual void beginTag(const AttributeList &attributes)
      {
	string name;

	this->checkParentTag("predicates");

	if (!attributes["name"].to(name))
	  throw runtime_error("expected attribute name for tag "
			      "<predicate>");

	int idPred=this->parser->defineSymbol(name,PredicateType);
	this->parser->getCallback()->beginPredicate(name,idPred);
      }

      virtual void text(const XMLString txt, bool last)
      {
      }

      virtual void endTag()
      {
	this->parser->getCallback()->endPredicate();
      }


      const VariableInfo *getParmInfo() 
      {
	return &parametersInfo;
      }

    private:
      VariableInfo parametersInfo; // information collected about the
				   // formal parameters

      // this class may alter parametersInfo
      friend class XMLParser<Callback,XMLString,XMLLibraryAttributeList>::ParametersTagAction;  
    };


    class ParametersTagAction : public TagAction
    {
    private:
      bool formalParameters; // formal or effective parameters ?
      int pos; // position of a parameter
      XMLLibraryString formalparms;

    public:
      ParametersTagAction(XMLParser *parser): TagAction(parser) {}

      virtual const char *getTagName() {return "parameters";}

      virtual void beginTag(const AttributeList &attributes)
      {
	pos=0;
	formalparms.clear();

	const char *parentTag=this->parser->getParentTagAction()->getTagName();

	if (strcmp(parentTag,"predicate")==0)
	{
	  formalParameters=true;
	}
	else
	  if (strcmp(parentTag,"constraint")==0)
	  {
	    formalParameters=false;
	  }
	  else
	    throw runtime_error("unexpected <parameters> tag");

	if (!formalParameters)
	{
	  // push a mark to denote the start of the parameters list
	  this->parser->operandStack.push_back(NULL);

	  this->parser->stateStack.front().setAllAbridgedAllowed(true);
	}
	else
	{
	  this->parser->stateStack.front().subtagAllowed=false;
	}
      }

      virtual void text(const XMLString txt, bool last)
      {
	if (formalParameters)
	  formalparms.append(txt);
	else
	  if (!txt.isWhiteSpace())
	    throw runtime_error("<parameters> tag should not have meaningful text");
      }

      virtual void endTag()
      {
	if (formalParameters)
	  parseFormalParameters(formalparms);
	else
	{
	  // create a new list and gather all parameters stored on the
	  // operandStack

	  ASTList *args=new ASTList;

	  while(!this->parser->operandStack.empty())
	  {
	    AST *val=this->parser->operandStack.back();
	    this->parser->operandStack.pop_back();

	    if (val==NULL)
	      break; // end of parameters list

	    args->push_front(val);
	  }

	  // push the list
	  this->parser->operandStack.push_back(args);

#ifdef debug
	  cout << "args=";
	  args->postfixExpression(cout);
	  cout << endl;
#endif
	}
      }

    private:
      void parseFormalParameters(const XMLString txt)
      {
	typename XMLString::Tokenizer tokenizer(txt);

	PredicateTagAction *parent=
	  dynamic_cast<PredicateTagAction *>(this->parser
					     ->getParentTagAction());
	if (parent==NULL)
	  throw runtime_error("internal error: couldn't get my parent or wasn't a <predicate>");

	VariableInfo &parmInfo=parent->parametersInfo;

	while(tokenizer.hasMoreTokens())
	{
	  XMLString thetype=tokenizer.nextToken();

	  if (!tokenizer.hasMoreTokens())
	    throw runtime_error("missing formal parameter name after type");

	  XMLString thename=tokenizer.nextToken();

	  string name,type;
	  thename.to(name);
	  thetype.to(type);

	  parmInfo[name].id=pos;
	  // ??? parmInfo[name].type=...

	  this->parser->getCallback()->addFormalParameter(pos,name,type);

	  ++pos;
	}
      }

    };

    class ConstraintTagAction;

    class SymbTagAction : public TagAction
    {
    public:
      SymbTagAction(XMLParser *parser, NodeType symb, const char *tag): 
	TagAction(parser),tag(tag),symb(symb) {}

      virtual const char *getTagName() {return tag.c_str();}

      virtual void beginTag(const AttributeList &attributes)
      {
	this->parser->operandStack.push_back(new ASTSymb(symb));

	this->parser->stateStack.front().setAllAbridgedAllowed(false);
	this->parser->stateStack.front().subtagAllowed=false;
      }

      virtual void text(const XMLString txt, bool last)
      {
	throw runtime_error("<"+tag+"/> tag must not have any text");
      }

      virtual void endTag()
      {
      }

    private:
      const string tag;
      NodeType symb;
    };

    class IntegerConstantTagAction : public TagAction
    {
    public:
      IntegerConstantTagAction(XMLParser *parser): TagAction(parser) {}

      virtual const char *getTagName() {return "i";}

      virtual void beginTag(const AttributeList &attributes)
      {
	this->parser->stateStack.front().setAllAbridgedAllowed(false);
	this->parser->stateStack.front().subtagAllowed=false;
      }

      virtual void text(const XMLString txt, bool last)
      {
	int val;

	if (!txt.to(val))
	  throw runtime_error("invalid integer in <i> element");

	this->parser->operandStack.push_back(new ASTInteger(val));
      }

      virtual void endTag()
      {
      }
    };

    class VarTagAction : public TagAction
    {
    public:
      VarTagAction(XMLParser *parser): TagAction(parser) {}

      virtual const char *getTagName() {return "var";}

      virtual void beginTag(const AttributeList &attributes)
      {
	string name;

	this->parser->stateStack.front().setAllAbridgedAllowed(false);
	this->parser->stateStack.front().subtagAllowed=false;

	if (!attributes["name"].to(name))
	  throw runtime_error("expected attribute name for tag "
			      "<var>");

	pushVariable(name);
      }

      virtual void text(const XMLString txt, bool last)
      {
	if (!txt.isWhiteSpace())
	  throw runtime_error("<var> tag should not have meaningful text");
      }

      virtual void endTag()
      {
      }

      void pushVariable(const string &name)
      {
	int id=this->parser->getSymbolId(name,VariableType);

	this->parser->operandStack.push_back(new ASTVar(name,id));
      }
    };

    class ListTagAction : public TagAction
    {
    public:
      ListTagAction(XMLParser *parser): TagAction(parser) {}

      virtual const char *getTagName() {return "list";}

      virtual void beginTag(const AttributeList &attributes)
      {
	// push a mark to denote the start of the list
	this->parser->operandStack.push_back(NULL);

	this->parser->stateStack.front().setAllAbridgedAllowed(true);
      }

      virtual void text(const XMLString txt, bool last)
      {
	if (!txt.isWhiteSpace())
	  throw runtime_error("<list> tag should not have meaningful text");
      }

      virtual void endTag()
      {
	// create a new list and gather all items stored on the
	// operandStack (but keep in mind that the

	ASTList *list=new ASTList;

	while(!this->parser->operandStack.empty())
	{
	  AST *val=this->parser->operandStack.back();
	  this->parser->operandStack.pop_back();

	  if (val==NULL)
	    break; // end of list

	  list->push_front(val);
	}

	// push the list
	this->parser->operandStack.push_back(list);

#ifdef debug
	cout << "list=";
	list->postfixExpression(cout);
	cout << endl;
#endif
      }
    };

    class DictTagAction : public TagAction
    {
    public:
      DictTagAction(XMLParser *parser): TagAction(parser) 
      {
	conventionalKeyOrder=NULL;
      }

      virtual const char *getTagName() {return "dict";}

      virtual void beginTag(const AttributeList &attributes)
      {
	// push a mark to denote the start of the dictionary
	this->parser->operandStack.push_back(NULL);

	this->parser->stateStack.front().setAllAbridgedAllowed(true);
      }

      virtual void text(const XMLString txt, bool last)
      {
	if (!txt.isWhiteSpace())
	  throw runtime_error("<dict> tag should not have meaningful text");
      }

      virtual void endTag()
      {
	// create a new dictionary and gather all <key,values> stored
	// on the operandStack

	ASTDict *dict=new ASTDict;

	deque<AST *> &stack=this->parser->operandStack;

	// check if the keys were provided or if we should use the
	// conventional order.
	int nbKey=0,nbVal=0,base;
	for(base=stack.size()-1;base>=0 && stack[base]!=NULL;--base)
	  if (stack[base]->getType()==_DICTKEY)
	    ++nbKey;
	  else
	    ++nbVal;

	if (base<0)
	  throw runtime_error("internal error: lost dict mark");
	else
	  ++base;

	if (nbKey==0)
	{
	  // use conventional order
	  if (!conventionalKeyOrder)
	    throw runtime_error("no conventional key order defined for this dictionary");

	  const vector<string> &keyOrder=*conventionalKeyOrder;

	  unsigned int i;
	  for(i=0;i<stack.size()-base && i<keyOrder.size();++i)
	  {
	    if (stack[base+i]->getType()!=SYMB_NIL)
	      dict->setArg(keyOrder[i],stack[base+i]);
	    else
	    {
	      // delete the <nil> objects
	      delete stack[base+i];
	    }
	  }

	  if (i!=keyOrder.size())
	    throw runtime_error("not enough values in dictionary defined with conventional order");

	  assert(stack[base-1]==NULL);
	  stack.resize(base-1); // pop values from the stack
	}
	else
	{
	  while(!this->parser->operandStack.empty())
	  {
	    AST *val=this->parser->operandStack.back();
	    this->parser->operandStack.pop_back();

	    if (val==NULL)
	      break; // end of dictionary

	    if (this->parser->operandStack.empty())
	      throw runtime_error("internal error");

	    ASTDictKey *key=dynamic_cast<ASTDictKey *>
	      (this->parser->operandStack.back());

	    this->parser->operandStack.pop_back();

	    if (key==NULL)
	      throw runtime_error("missing key in dictionary");

	    dict->setArg(key->getKey(),val);
	    delete key;
	  }
	}
	// push the dictionary
	stack.push_back(dict);

#ifdef debug
	cout << "dict=";
	dict->postfixExpression(cout);
	cout << endl;
	//this->parser->dumpOperandStack();
#endif
      }

      void setConventionalKeyOrder(const vector<string> *order=NULL)
      {
	conventionalKeyOrder=order;
      }

    private:
      const vector<string> *conventionalKeyOrder;
    };

    class DictEntryTagAction : public TagAction
    {
    public:
      DictEntryTagAction(XMLParser *parser): TagAction(parser) {}

      virtual const char *getTagName() {return "entry";}

      virtual void beginTag(const AttributeList &attributes)
      {
	string key;

	if (!attributes["key"].to(key))
	  throw runtime_error("expected attribute key for tag "
			      "<entry>");

	this->parser->operandStack.push_back(new ASTDictKey(key));

	this->parser->stateStack.front().setAllAbridgedAllowed(true);
      }

      virtual void text(const XMLString txt, bool last)
      {
	if (!txt.isWhiteSpace())
	  throw runtime_error("<entry> tag should not have meaningful text");
      }

      virtual void endTag()
      {
      }
    };

    class ExpressionTagAction : public TagAction
    {
    private:
      bool inPredicateDefinition;

    public:
      ExpressionTagAction(XMLParser *parser): TagAction(parser) {}

      virtual const char *getTagName() {return "expression";}

      virtual void beginTag(const AttributeList &attributes)
      {
	const char *parentTag=this->parser->getParentTagAction()->getTagName();

	if (strcmp(parentTag,"predicate")==0)
	{
	  inPredicateDefinition=true;
	}
	else
	  if (strcmp(parentTag,"constraint")==0)
	  {
	    inPredicateDefinition=false;
	  }
	  else
	    throw runtime_error("unexpected <expression> tag");
      }

      virtual void text(const XMLString txt, bool last)
      {
	if (!txt.isWhiteSpace())
	  throw runtime_error("<expression> tag should not have meaningful text");
      }

      virtual void endTag()
      {
      }

      bool isPredicateDefinition() {return inPredicateDefinition;}
    };

    class FunctionalTagAction : public TagAction
    {
    public:
      FunctionalTagAction(XMLParser *parser): TagAction(parser) {}

      virtual const char *getTagName() {return "functional";}

      virtual void beginTag(const AttributeList &attributes)
      {
	this->checkParentTag("expression");
	expr.clear();
      }

      virtual void text(const XMLString txt, bool last)
      {
	// gather all chunks together
	txt.appendTo(expr);
      }

      virtual void endTag()
      {
	ExpressionTagAction *parentExpressionTag=
	  dynamic_cast<ExpressionTagAction *>(this->parser->getParentTagAction(1));
	if (parentExpressionTag==NULL)
	  throw runtime_error("internal error: couldn't get my parent or wasn't an <expression>");

	bool predicateDefinition=parentExpressionTag->isPredicateDefinition();

	if (predicateDefinition)
	{
	  PredicateTagAction *parentPredicateTag=
	    dynamic_cast<PredicateTagAction *>(this->parser
					       ->getParentTagAction(2));
	  if (parentPredicateTag==NULL)
	    throw runtime_error("internal error: couldn't get my parent or wasn't a <predicate>");

	  this->parser->exprParser.setVarInfo(parentPredicateTag->getParmInfo());

	  AST *tree=this->parser->exprParser.prefixParser(expr);

	  this->parser->exprParser.unsetVarInfo();

	  this->parser->xmitExpressionToSolver(tree);
	}
	else
	  throw runtime_error("expressions not yet suported as effective parameters"); // ???
      }

    private:
      string expr;
    };

    class InfixTagAction : public TagAction
    {
    public:
      InfixTagAction(XMLParser *parser): TagAction(parser) {}

      virtual const char *getTagName() {return "infix";}

      virtual void beginTag(const AttributeList &attributes)
      {
	this->checkParentTag("expression");
	expr.clear();
      }

      virtual void text(const XMLString txt, bool last)
      {
	// gather all chunks together
	txt.appendTo(expr);
      }

      virtual void endTag()
      {

	cout << "WARNING: ignoring infix expression: " <<  expr << endl;
	return; // ??? currently unsupported

	ExpressionTagAction *parentExpressionTag=
	  dynamic_cast<ExpressionTagAction *>(this->parser
					      ->getParentTagAction(1));
	if (parentExpressionTag==NULL)
	  throw runtime_error("internal error: couldn't get my parent or wasn't an <expression>");

	bool predicateDefinition=parentExpressionTag->isPredicateDefinition();

	if (predicateDefinition)
	{
	  PredicateTagAction *parentPredicateTag=
	    dynamic_cast<PredicateTagAction *>(this->parser
					       ->getParentTagAction(2));
	  if (parentPredicateTag==NULL)
	    throw runtime_error("internal error: couldn't get my parent or wasn't a <predicate>");

	  this->parser->exprParser.setVarInfo(parentPredicateTag->getParmInfo());

	  AST *tree=this->parser->exprParser.infixParser(expr);

	  this->parser->exprParser.unsetVarInfo();

	  ostringstream tmp;

	  // prefix
	  tree->expression(tmp);
	  this->parser->getCallback()->predicateExpression(tmp.str());

	  delete tree;
	}
	else
	  throw runtime_error("expressions not yet supported as effective parameters"); // ???
      }

    private:
      string expr;
    };

    class PostfixTagAction : public TagAction
    {
    public:
      PostfixTagAction(XMLParser *parser): TagAction(parser) {}

      virtual const char *getTagName() {return "postfix";}

      virtual void beginTag(const AttributeList &attributes)
      {
	this->checkParentTag("expression");
	expr.clear();
      }

      virtual void text(const XMLString txt, bool last)
      {
	// gather all chunks together
	txt.appendTo(expr);
      }

      virtual void endTag()
      {
	cout << "WARNING: ignoring postfix expression: " <<  expr << endl;
	return; // ??? currently unsupported

	ExpressionTagAction *parentExpressionTag=
	  dynamic_cast<ExpressionTagAction *>(this->parser
					      ->getParentTagAction(1));
	if (parentExpressionTag==NULL)
	  throw runtime_error("internal error: couldn't get my parent or wasn't an <expression>");

	bool predicateDefinition=parentExpressionTag->isPredicateDefinition();

	if (predicateDefinition)
	{
	  PredicateTagAction *parentPredicateTag=
	    dynamic_cast<PredicateTagAction *>(this->parser
					       ->getParentTagAction(2));
	  if (parentPredicateTag==NULL)
	    throw runtime_error("internal error: couldn't get my parent or wasn't a <predicate>");

	  this->parser->exprParser.setVarInfo(parentPredicateTag->getParmInfo());

	  AST *tree=this->parser->exprParser.postfixParser(expr);

	  this->parser->exprParser.unsetVarInfo();

	  ostringstream tmp;

	  // prefix
	  tree->expression(tmp);
	  this->parser->getCallback()->predicateExpression(tmp.str());

	  delete tree;
	}
	else
	  throw runtime_error("expressions not yet suported as effective parameters"); // ???
      }

    private:
      string expr;
    };



    class ConstraintsTagAction : public TagAction
    {
    public:
      ConstraintsTagAction(XMLParser *parser): TagAction(parser) {}

      virtual const char *getTagName() {return "constraints";}

      virtual void beginTag(const AttributeList &attributes)
      {
	int nbConstraints;

	this->checkParentTag("instance");

	if (!attributes["nbConstraints"].to(nbConstraints))
	  throw runtime_error("expected attribute nbConstraints for tag "
			      "<constraints>");

	this->parser->getCallback()->beginConstraintsSection(nbConstraints);
      }

      virtual void text(const XMLString txt, bool last)
      {
	if (!txt.isWhiteSpace())
	  throw runtime_error("<constraints> tag should not have meaningful text");
      }

      virtual void endTag()
      {
	this->parser->getCallback()->endConstraintsSection();
      }
    };

    class ConstraintTagAction : public TagAction
    {
    public:
      ConstraintTagAction(XMLParser *parser): TagAction(parser) {}

      virtual const char *getTagName() {return "constraint";}

      virtual void beginTag(const AttributeList &attributes)
      {
	string name,reference;
	int arity;

	initialOperandStackSize=this->parser->operandStack.size();

	this->checkParentTag("constraints");

	if (!attributes["arity"].to(arity))
	  throw runtime_error("expected attribute arity for tag "
			      "<constraint>");

	if (!attributes["name"].to(name))
	  throw runtime_error("expected attribute name for tag "
			      "<constraint>");

	if (!attributes["reference"].to(reference))
	  throw runtime_error("expected attribute reference for tag "
			      "<constraint>");

	if (reference.length()>7 && reference.substr(0,7)=="global:")
	{
	  // global constraints are case insensitive
	  for(unsigned int i=0;i<reference.size();++i)
	    reference[i]=tolower(reference[i]);
	}

	SymbolInfo info=this->parser->getSymbolInfo(reference);

	// temporary way to set up conventional dict. order ???
	if (info.convOrder)
	  this->parser->dictTagAction->setConventionalKeyOrder(info.convOrder);

	XMLString scope=attributes["scope"];
	if (scope.isNull())
	  throw runtime_error("expected attribute scope for tag "
			      "<constraint>");

	// clear scopeList
	scopeList.clear();

	typename XMLString::Tokenizer tokenizer(scope);

	while(tokenizer.hasMoreTokens())
	{
	  string var;

	  tokenizer.nextToken().to(var);
	  int idVar=this->parser->getSymbolId(var,VariableType);
	  scopeList.push_back(new ASTVar(var,idVar));

	  scopeInfo[var].id=idVar;
	  //scopeInfo[var].type=...; ???
	}

	int idConstr=this->parser->defineSymbol(name,ConstraintType);
	this->parser->getCallback()->beginConstraint(name,idConstr,arity,
						     reference,
						     info.type,
						     info.id,
						     scopeList);
      }

      virtual void text(const XMLString txt, bool last)
      {
	if (!txt.isWhiteSpace())
	  throw runtime_error("<constraint> tag should not have meaningful text");
      }

      virtual void endTag()
      {
	if (this->parser->operandStack.size()-initialOperandStackSize==1)
	{
	  ASTList *args=dynamic_cast<ASTList *>
	    (this->parser->operandStack.back());

	  this->parser->operandStack.pop_back();

	  if (args==NULL)
	    throw runtime_error("constraint parameters doesn't appear to be a list");
	  this->parser->getCallback()->constraintParameters(*args);

	  delete args;
	}
	else
	  if (this->parser->operandStack.size()==(unsigned int)initialOperandStackSize)
	  {
	    // no parameter given to the constraint, pass the scope list
	    this->parser->getCallback()->constraintParameters(scopeList);
	  }
	  else
	    throw runtime_error("internal error while handling constraint parameters");

	this->parser->getCallback()->endConstraint();

	this->parser->dictTagAction->setConventionalKeyOrder();
      }

      const VariableInfo &getScopeInfo()
      {
	return scopeInfo;
      }

    private:
      ASTList scopeList;
      VariableInfo scopeInfo;
      int initialOperandStackSize;
    };

    class ExtensionTagAction : public TagAction
    {
    public:
      ExtensionTagAction(XMLParser *parser): TagAction(parser) {}

      virtual const char *getTagName() {return "extension";}

      virtual void beginTag(const AttributeList &attributes)
      {
	// ignore
      }

      virtual void text(const XMLString txt, bool last)
      {
	// ignore
      }

      virtual void endTag()
      {
	// ignore
      }
    };

    /*************************************************************************
     *
     * MathML tags
     *
     *************************************************************************/

    class MathMLFunctionTagAction : public TagAction
    {
    public:
      MathMLFunctionTagAction(XMLParser *parser,
			      ASTAbstractFunction *ast): TagAction(parser) {}

      virtual const char *getTagName()=0;

      virtual void beginTag(const AttributeList &attributes)
      {
	this->checkParentTag("apply");
      }

      virtual void text(const XMLString txt, bool last)
      {
	if (!txt.isWhiteSpace())
	  throw runtime_error("MathML function tag should not have meaningful text");
      }

      virtual void endTag()
      {
      }
    };

    class MathTagAction : public TagAction
    {
    public:
      MathTagAction(XMLParser *parser): TagAction(parser) {}

      virtual const char *getTagName() {return "math";}

      virtual void beginTag(const AttributeList &attributes)
      {
	this->checkParentTag("expression");
      }

      virtual void text(const XMLString txt, bool last)
      {
      }

      virtual void endTag()
      {
      }
    };

    class ApplyTagAction : public TagAction
    {
    public:
      ApplyTagAction(XMLParser *parser): TagAction(parser) {}

      virtual const char *getTagName() {return "apply";}

      virtual void beginTag(const AttributeList &attributes)
      {
	//this->checkParentTag("expression");
      }

      virtual void text(const XMLString txt, bool last)
      {
      }

      virtual void endTag()
      {
      }
    };

    class CITagAction : public TagAction
    {
    public:
      CITagAction(XMLParser *parser): TagAction(parser) {}

      virtual const char *getTagName() {return "ci";}

      virtual void beginTag(const AttributeList &attributes)
      {
	//this->checkParentTag("expression");

	// ??? don't allow subtags
      }

      virtual void text(const XMLString txt, bool last)
      {
	string var;
	txt.to(var);
	new ASTVar(var);
      }

      virtual void endTag()
      {
      }
    };


    class CNTagAction : public TagAction
    {
    public:
      CNTagAction(XMLParser *parser): TagAction(parser) {}

      virtual const char *getTagName() {return "cn";}

      virtual void beginTag(const AttributeList &attributes)
      {
	//this->checkParentTag("expression");

	// ??? don't allow subtags
      }

      virtual void text(const XMLString txt, bool last)
      {
	int val;
	txt.to(val);
	new ASTInteger(val);
      }

      virtual void endTag()
      {
      }
    };

    class TrueTagAction : public TagAction
    {
    public:
      TrueTagAction(XMLParser *parser): TagAction(parser) {}

      virtual const char *getTagName() {return "true";}

      virtual void beginTag(const AttributeList &attributes)
      {
	this->parser->stateStack.front().setAllAbridgedAllowed(false);
	this->parser->stateStack.front().subtagAllowed=false;

	this->parser->operandStack.push_back(new ASTBoolean(true));
      }

      virtual void text(const XMLString txt, bool last)
      {
	throw runtime_error("<true/> tag should not contain text");
      }

      virtual void endTag()
      {
      }
    };

    class FalseTagAction : public TagAction
    {
    public:
      FalseTagAction(XMLParser *parser): TagAction(parser) {}

      virtual const char *getTagName() {return "false";}

      virtual void beginTag(const AttributeList &attributes)
      {
	this->parser->stateStack.front().setAllAbridgedAllowed(false);
	this->parser->stateStack.front().subtagAllowed=false;

	this->parser->operandStack.push_back(new ASTBoolean(false));
      }

      virtual void text(const XMLString txt, bool last)
      {
	throw runtime_error("<false/> tag should not contain text");
      }

      virtual void endTag()
      {
      }
    };


  public:
    XMLParser(Callback &cb)
    {
      preferredSyntax=TREE;

      this->cb=&cb;

      unknownTagHandler=new UnknownTagAction(this);

      // symbolic constants
      registerTagAction(tagList,new SymbTagAction(this,SYMB_NIL,"nil"));
      registerTagAction(tagList,new SymbTagAction(this,SYMB_EQ,"eq"));
      registerTagAction(tagList,new SymbTagAction(this,SYMB_NE,"ne"));
      registerTagAction(tagList,new SymbTagAction(this,SYMB_GE,"ge"));
      registerTagAction(tagList,new SymbTagAction(this,SYMB_GT,"gt"));
      registerTagAction(tagList,new SymbTagAction(this,SYMB_LE,"le"));
      registerTagAction(tagList,new SymbTagAction(this,SYMB_LT,"lt"));

      // base constructions
      registerTagAction(tagList,varTagAction=new VarTagAction(this));
      registerTagAction(tagList,new IntegerConstantTagAction(this));
      registerTagAction(tagList,new ListTagAction(this));
      registerTagAction(tagList,dictTagAction=new DictTagAction(this));
      registerTagAction(tagList,new DictEntryTagAction(this));


      registerTagAction(tagList,new InstanceTagAction(this));
      registerTagAction(tagList,new PresentationTagAction(this));
      registerTagAction(tagList,new DomainsTagAction(this));
      registerTagAction(tagList,new DomainTagAction(this));
      registerTagAction(tagList,new VariablesTagAction(this));
      registerTagAction(tagList,new VariableTagAction(this));
      registerTagAction(tagList,new RelationsTagAction(this));
      registerTagAction(tagList,new RelationTagAction(this));
      registerTagAction(tagList,new PredicatesTagAction(this));
      registerTagAction(tagList,new PredicateTagAction(this));

      registerTagAction(tagList,new ParametersTagAction(this));
      registerTagAction(tagList,new ExpressionTagAction(this));

      registerTagAction(tagList,new FunctionalTagAction(this));
      registerTagAction(tagList,new InfixTagAction(this));
      registerTagAction(tagList,new PostfixTagAction(this));

      registerTagAction(tagList,new ConstraintsTagAction(this));
      registerTagAction(tagList,new ConstraintTagAction(this));

      registerTagAction(tagList,new ExtensionTagAction(this));

      // mathml
      //registerTagAction(tagList,new MathTagAction(this));
      //registerTagAction(tagList,new ApplyTagAction(this));
      //registerTagAction(tagList,new CITagAction(this));
      //registerTagAction(tagList,new TrueTagAction(this));
      //registerTagAction(tagList,new FalseTagAction(this));

      for(int i=0;i<UndefinedType;++i)
	nextId[i]=0;

      symbolDict["global:alldifferent"]=SymbolInfo(0,GlobalConstraintType);


      cumulative_dictOrder.push_back("origin");
      cumulative_dictOrder.push_back("duration");
      cumulative_dictOrder.push_back("end");
      cumulative_dictOrder.push_back("height");

      symbolDict["global:cumulative"]=
	SymbolInfo(4,GlobalConstraintType,&cumulative_dictOrder);

      weightedSum_dictOrder.push_back("coef");
      weightedSum_dictOrder.push_back("var");

      symbolDict["global:weightedsum"]=
	SymbolInfo(5,GlobalConstraintType,&weightedSum_dictOrder);

      symbolDict["global:element"]=
	SymbolInfo(6,GlobalConstraintType);

      globalCardinality_dictOrder.push_back("value");
      globalCardinality_dictOrder.push_back("card");

      symbolDict["global:global_cardinality"]=
	SymbolInfo(7,GlobalConstraintType,&globalCardinality_dictOrder);
    }

    ~XMLParser()
    {
      delete unknownTagHandler;
      for(typename TagActionList::iterator it=tagList.begin();
	  it!=tagList.end();++it)
	delete (*it).second;
    }

    Callback *getCallback()
    {
      return cb;
    }

    /**
     * get the parent tag action that is n levels higher in the current
     * branch of the XML parse tree
     */
    TagAction *getParentTagAction(int n=1)
    {
      if (n<0 || n>=(int)actionStack.size())
	return NULL;

      return actionStack[n];
    }


    /**
     * choose the representation of an expression that the parser must
     * use to transmit an expression to the solver
     *
     * @param syntax 
     *   @arg TREE    transmit as an Abstract Syntax Tree (AST)
     *   @arg PREFIX  transmit as a string (prefix notation)
     *   @arg INFIX_C transmit as a string (C infix notation)
     *   @arg POSTFIX transmit as a string (postfix notation)
     *   @arg MATHML  transmit as a string (MathML expression) 
     *
     */
    void setPreferredExpressionRepresentation(Syntax syntax)
    {
      preferredSyntax=syntax;
    }


    /**
     * callbacks from the XML parser
     */
    void startElement(XMLString name, const AttributeList &attributes)
    {
      // consume the last tokens before we switch to the next element
      if (!textLeft.empty())
      {
	handleAbridgedNotation(textLeft,true);
	textLeft.clear();
      }

      if (!stateStack.empty() && !stateStack.front().subtagAllowed)
	throw runtime_error("this element must not contain any element");

      typename TagActionList::iterator iAction=tagList.find(name);
      TagAction *action;

      if (iAction!=tagList.end())
      {
	action=(*iAction).second;

	// ???
	//if (!action->isActivated())
	//  throw runtime_error("unexpected tag");
      }
      else
      {
	// add a handler to ignore the text and end element
	action=unknownTagHandler;

	cerr << "WARNING: unknown tag <" << name 
	     << "> will be ignored"
	     << endl;
      }

      stateStack.push_front(State());
      actionStack.push_front(action);
      action->beginTag(attributes);
    }

    void endElement(XMLString name)
    {
      // consume the last tokens 
      if (!textLeft.empty())
      {
	handleAbridgedNotation(textLeft,true);
	textLeft.clear();
      }

      typename TagActionList::iterator iAction=tagList.find(name);
      
      if (iAction!=tagList.end())
	(*iAction).second->endTag();

      actionStack.pop_front();
      stateStack.pop_front();
    }

    void characters(XMLString chars)
    {
      //cout << "chars=" << chars << "#" << endl;

      if (actionStack.empty())
      {
	if (chars.isWhiteSpace())
	  return;
	else
	  throw runtime_error("Text found outside any tag");
      }

      if (!textLeft.empty())
      {
	// break at first space, concatenate with textLeft and call
	// text()
	typename XMLString::iterator it=chars.begin(),end=chars.end();

	while(it!=end && !it.isWhiteSpace())
	{
	  textLeft.append(*it);
	  ++it;
	}

	while(it!=end && it.isWhiteSpace())
	{
	  textLeft.append(*it);
	  ++it;
	}

	handleAbridgedNotation(textLeft,false);
	textLeft.clear();

	chars=chars.substr(it,chars.end());
      }

      // break after last space, call text() with the first part and
      // store the last part in textLeft

      typename XMLString::iterator it,brk;

      brk=chars.end();
      while(brk!=chars.begin())
      {
	--brk;
	if (brk.isWhiteSpace())
	{
	  ++brk;
	  break;
	}
      }

      for(it=brk;it!=chars.end();++it)
	textLeft.append(*it);

      chars=chars.substr(chars.begin(),brk);

      if (!chars.empty())
	handleAbridgedNotation(chars,false);
    }

    void handleAbridgedNotation(XMLString chars, bool lastChunk)
    {
      typename XMLString::iterator it,beg,end;

      it=chars.begin();
      beg=chars.begin();
      end=chars.end();

      while(it!=end)
      {
	// skip spaces
	while(it!=end && it.isWhiteSpace())
	  ++it;

	if (it==end)
	  break;

	// does it look like a special token?
	int c=*it;
	switch(c)
	{
	case '[':
	  //cout << "got [" << endl;
	  if (stateStack.front().abridgedListAllowed)
	  {
	    if (beg!=it)
	      actionStack.front()->text(chars.substr(beg,it),true);


	    stateStack.push_front(State());
	    actionStack.push_front(tagList["list"]);
	    actionStack.front()->beginTag(AttributeList());
	    ++it;
	    beg=it;
	    continue;
	  }
	  break;
	case ']':
	  //cout << "got ]" << endl;
	  if (stateStack.front().abridgedListAllowed)
	  {
	    if (beg!=it)
	      actionStack.front()->text(chars.substr(beg,it),true);

	    actionStack.front()->endTag();
	    actionStack.pop_front();
	    stateStack.pop_front();
	    ++it;
	    beg=it;
	    continue;
	  }
	  break;

	case '{':
	  //cout << "got {" << endl;
	  if (stateStack.front().abridgedDictAllowed)
	  {
	    if (beg!=it)
	      actionStack.front()->text(chars.substr(beg,it),true);

	    stateStack.push_front(State());
	    actionStack.push_front(tagList["dict"]);
	    actionStack.front()->beginTag(AttributeList());
	    ++it;
	    beg=it;
	    continue;
	  }
	  break;
	case '}':
	  //cout << "got }" << endl;
	  if (stateStack.front().abridgedDictAllowed)
	  {
	    if (beg!=it)
	      actionStack.front()->text(chars.substr(beg,it),true);

	    actionStack.front()->endTag();
	    actionStack.pop_front();
	    stateStack.pop_front();
	    ++it;
	    beg=it;
	    continue;
	  }
	  break;

	case '/':
	  //cout << "got /" << endl;
	  if (stateStack.front().abridgedDictAllowed)
	  {
	    if (beg!=it)
	      actionStack.front()->text(chars.substr(beg,it),true);

	    typename XMLString::iterator begKey,endKey;

	    ++it;
	    begKey=it;
	    endKey=it;
	    while (endKey!=end && isalpha(*endKey))
	      ++endKey;

	    string key;

	    chars.substr(begKey,endKey).to(key);

	    operandStack.push_back(new ASTDictKey(key));

	    it=endKey;
	    beg=endKey;
	    continue;
	  }
	  break;

	default:
	  if (stateStack.front().abridgedVariableAllowed && 
	      (c=='_'||isalpha(c)))
	  {
	    typename XMLString::iterator endVar=it;
	    ++endVar;
	    while (endVar!=end && 
		   (isdigit(*endVar)||isalpha(*endVar)||*endVar=='_'))
	      ++endVar;

	    if (beg!=it)
	      actionStack.front()->text(chars.substr(beg,it),true);

	    string name;
	    chars.substr(it,endVar).to(name);
	    varTagAction->pushVariable(name);

	    it=endVar;
	    beg=it;
	    continue;
	  }
	  else
	    if (stateStack.front().abridgedIntegerAllowed && 
		(isdigit(c)||c=='+'||c=='-'))
	    {
	      typename XMLString::iterator endNum=it;
	      ++endNum;
	      while (endNum!=end && isdigit(*endNum))
		++endNum;

	      // ??? test for . (real) and .. (interval)

	      if (beg!=it)
		actionStack.front()->text(chars.substr(beg,it),true);

	      stateStack.push_front(State());
	      actionStack.push_front(tagList["i"]);
	      actionStack.front()->beginTag(AttributeList());
	      actionStack.front()->text(chars.substr(it,endNum),true);
	      actionStack.front()->endTag();
	      actionStack.pop_front();
	      stateStack.pop_front();
	      it=endNum;
	      beg=it;
	      continue;
	    }
	  break;
	}

	// no special token
	while(it!=end && !it.isWhiteSpace())
	  ++it;
      }

      if (beg!=end)
	actionStack.front()->text(chars.substr(beg,end),lastChunk);
    }

    void startDocument()
    {
      clearStacks();
    }


    void endDocument()
    {
    }

  protected:
    void clearStacks()
    {
      clearOperandStack();
      stateStack.clear();
    }

    void clearOperandStack()
    {
      for(unsigned int i=0;i<operandStack.size();++i)
	delete operandStack[i];

      operandStack.clear();
    }

    void dumpOperandStack()
    {
      cout << "operand stack:\n";

      for(int i=0;i<operandStack.size();++i)
      {
	cout << "  ";
	if (operandStack[i]==NULL)
	  cout << "NULL";
	else
	  operandStack[i]->postfixExpression(cout);

	cout << endl;
      }
    }

  protected:
    // text which is left for the next call to characters() because it
    // may not be a complete token
    XMLLibraryString textLeft;

    // specific actions
    VarTagAction *varTagAction;
    DictTagAction *dictTagAction;
    TagAction *unknownTagHandler; // handler to help ignore all unknown tags

    // stack of operands to construct list, dictionaries, predicate
    // parameters and so on
    deque<AST *> operandStack;


    // ??? temporary solution
    // conventional order for some global constraints
    vector<string> cumulative_dictOrder,weightedSum_dictOrder,globalCardinality_dictOrder;
  }; // class XMLParser

} // namespace


// Local Variables:
// mode: C++
// End:

#endif

