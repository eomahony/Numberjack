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
#ifndef _XMLParser_Xerces_h_
#define _XMLParser_Xerces_h_

#error "The Xerces-C library is currently unsupported"

#include <iostream>
#include <stdexcept>

#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/util/TransService.hpp>
#include <xercesc/parsers/SAXParser.hpp>
#include <xercesc/util/OutOfMemoryException.hpp>
#include <xercesc/sax/HandlerBase.hpp>
#include <xercesc/framework/XMLFormatter.hpp>
#include <xercesc/util/XMLChar.hpp>
#include <xercesc/util/XMLString.hpp>
#include <xercesc/util/XMLStringTokenizer.hpp>
#include <xercesc/sax/InputSource.hpp>
#include <xercesc/util/BinInputStream.hpp>

#include "XMLParser.h"

namespace CSPXMLParser
{
  using namespace std;
  using namespace xercesc;

  /**
   * this class is meant to be passed by value so it must only
   * contain reference to the actual string
   */
  class XercesString
  {
  private:
    const XMLCh *t; // the string
    int *refCounter; // counter of references to the string
    // if refCounter is NULL, we don't have to release the
    // string. Otherwise, we must release the string as soon as it is
    // not referenced any more

  public:
    static const int npos=-1;

    class Tokenizer
    {
    private:
      XMLStringTokenizer tokenizer;
    public:
      Tokenizer(const XercesString s): tokenizer(s.t) {}

      bool hasMoreTokens()
      {
	return tokenizer.hasMoreTokens();
      }

      XercesString nextToken()
      {
	return tokenizer.nextToken();
      }
    };


    /**
     * initialize to null
     *
     * isNull() must return true on a string initialized by this
     * default constructor
     */
    XercesString()
    {
      t=NULL;
      refCounter=NULL;
    }

    XercesString(int size)
    {
      // how does xerces allocates its strings ???
      t=new XMLCh[size];
      refCounter=new int;
      *refCounter=1;

      XMLCh *p=const_cast<XMLCh *>(t);
      *p=0;
    }

    XercesString(const XMLCh *s)
    {
      t=s;
      refCounter=NULL;
    }

    XercesString(const char *s)
    {
      t=XMLString::transcode(s);
      refCounter=new int;
      *refCounter=1;
    }

    XercesString(const XercesString &s)
    {
      t=s.t;
      refCounter=s.refCounter;

      incRefCounter();
    }

    ~XercesString()
    {
      decRefCounter();
    }

    XercesString &operator =(const XercesString &s)
    {
      throw runtime_error("unimplemented");
      return *this;
    }

    bool isNull() const
    {
      return t==NULL;
    }

    /**
     * returns true iff the string contains only white space
     */
    bool isWhiteSpace() const
    {
      return XMLString::isAllWhiteSpace(t);
    }

    /**
     * returns the first character of the string
     */
    int firstChar() const
    {
      return *t;//??? how to get the first character in the UTF-16 string ?
    }

    // return the position of sub in this string or npos it not found
    int find(XercesString sub) const
    {
      return XMLString::patternMatch(t,sub.t);
    }

    XercesString substr(int pos, int count=npos)
    {
      if (count==npos)
	count=XMLString::stringLen(t)-pos;

      XercesString sub(count+1);

      XMLCh *p=const_cast<XMLCh *>(sub.t);
      for(int i=0;i<count;++i)
	p[i]=t[pos+i];

      p[count]=0; // terminator

      return sub;
    }

    bool operator ==(const XercesString &s) const
    {
      return XMLString::equals(t,s.t); 
    }

    bool operator !=(const XercesString &s) const
    {
      return !XMLString::equals(t,s.t); 
    }

    bool operator <(const XercesString &s) const
    {
      return XMLString::compareString(t,s.t)<0; 
    }

    /**
     * returns true iff value was successfully transfered in the variable
     */
    template<typename T>
    bool to(T &v) const
    {
      throw runtime_error("don't know how to convert to this type");
    }

    friend ostream &operator <<(ostream &s, XercesString str);

  protected:
    void decRefCounter()
    {
      if (refCounter && --*refCounter==0)
      {
	XMLString::release(const_cast<XMLCh**>(&t));
	delete refCounter;
      }
    }

    void incRefCounter()
    {
      if (refCounter)
	++*refCounter;
    }
  };

  template<>
  bool XercesString::to<int>(int &v) const
  {
    v=XMLString::parseInt(t);
    // how to check for errors ??? and empty string ???
    return true;
  }

  template<>
  bool XercesString::to<string>(string &v) const
  {
    char *tmp=XMLString::transcode(t);
    v=tmp;
    XMLString::release(&tmp);

    return true;
  }

  ostream &operator <<(ostream &s, XercesString str)
  {
    char *txt=XMLString::transcode(str.t);
    s << txt;
    XMLString::release(&txt);

    return s;
  }

  ostream &operator <<(ostream &s, const XMLCh *t)
  {
    char *txt=XMLString::transcode(t);
    s << txt;
    XMLString::release(&txt);

    return s;
  }

  /**
   * an interface to the liste of attributes.
   *
   * This[name] must return the value of attribute
   * if attribute doesn't exist This[name].isNull() must be true 
   */
  class XercesAttributeList
  {
  private:
    AttributeList &list;
  public:
    XercesAttributeList(AttributeList &l): list(l) {}

    const XercesString operator[] (const char *name) const
    {
      return XercesString(list.getValue(name));
    }
  };



  /**
   * @brief the SAX handler for the Xerces library
   */
  template<class Callback=CSPParserCallback>
  class XMLParser_Xerces: public HandlerBase
  {
  public:
    XMLParser_Xerces(Callback &callback): xercesLibrary(), cspParser(callback)
    {
      // the library must be initialized before the cspParser
    }

    int parse(istream &in)
    {
      return libraryParse(in);
    }

    int parse(const char *filename)
    {
      ifstream in(filename);
      return parse(in);
    }

    void setPreferredExpressionRepresentation(AST::Syntax syntax)
    {
      cspParser.setPreferredExpressionRepresentation(syntax);
    }
  
  protected:

    /**
     * so far, this is a simple copy from a SAX parser example
     */
    int libraryParse(istream &in)
    {
      //
      //  Create a SAX parser object. Then, according to what we were told on
      //  the command line, set it to validate or not.
      //
      SAXParser* parser = new SAXParser;
      //parser->setValidationScheme(false);
      parser->setDoNamespaces(false);
      parser->setDoSchema(false);
      parser->setValidationSchemaFullChecking(false);

      //
      //  Create the handler object and install it as the document and error
      //  handler for the parser-> Then parse the file and catch any exceptions
      //  that propagate out
      //
      int errorCode = 0;
      int errorCount = 0;
      try
      {
	StreamInputSource inputStream(in);

	parser->setDocumentHandler(this);
	parser->setErrorHandler(this);
	parser->parse(inputStream);
	errorCount = parser->getErrorCount();
      }
      catch (const SAXParseException &e)
      {
	cerr << "On line " << e.getLineNumber() 
	     << ", column " << e.getColumnNumber() << endl;
	cerr << e.getMessage() << endl;
	errorCode = 5;
      }
      catch (const OutOfMemoryException&)
      {
	cerr << "OutOfMemoryException" << endl;
	errorCode = 5;
      }
      catch (const XMLException& toCatch)
      {
	cerr << "\nAn error occurred\n  Error: "
	     << toCatch.getMessage()
	     << "\n" << endl;
	errorCode = 4;
      }

      if(errorCode) {
	XMLPlatformUtils::Terminate();
	return errorCode;
      }

      //
      //  Delete the parser itself.  Must be done prior to calling
      //  Terminate, below.
      //
      delete parser;

      // And call the termination method
      XMLPlatformUtils::Terminate();

      //cout << "Error count=" << errorCount << endl;

      return errorCount;
    }



    /*************************************************************************
     *
     * SAX Handler
     *
     *************************************************************************/

    void startDocument()
    {
#ifdef debug
      cout << "Parsing begins" << endl;
#endif
    }

    void endDocument()
    {
#ifdef debug
      cout << "Parsing ends" << endl;
#endif
    }

    void ignorableWhitespace(const XMLCh* const chars,
			     const unsigned int length) {}

    void processingInstruction(const XMLCh* const target, 
			       const XMLCh* const data) {}

    void startElement(const XMLCh* const name, AttributeList& attributes) 
    {
#ifdef debug
      cout << "  begin element " << name << endl;
      for(int i=0;i<attributes.getLength();++i)
      {
	cout << "    attribute " << attributes.getName(i)
	     << " = " << attributes.getValue(i) << endl;
      }
#endif

      cspParser.startElement(name,attributes);
    }

    void characters(const XMLCh* const chars, const unsigned int length) 
    {
#ifdef debug
      char *n=XMLString::transcode(chars);
      XMLString::trim(n);
      cout << "    chars '" << n << "'" << endl;
      XMLString::release(&n);
#endif

      // temporarily modify the const string to cut it at length (this
      // is probably not very nice but it avoids useless allocation)
      XMLCh *s=const_cast<XMLCh *>(chars);
      XMLCh backup=s[length];

      s[length]=0;
      cspParser.characters(s);

      s[length]=backup;
    }

    void endElement(const XMLCh* const name) 
    {
#ifdef debug
      cout << "  end element " << name << endl;
#endif

      cspParser.endElement(name);
    }


    // -----------------------------------------------------------------------
    //  Implementations of the SAX ErrorHandler interface
    // -----------------------------------------------------------------------
    void warning(const SAXParseException& exc) 
    {
      throw exc;
    }
    
    void error(const SAXParseException& exc) 
    {
      throw exc;
    }

    void fatalError(const SAXParseException& exc) 
    {
      throw exc;
    }

  private:

    /**
     * @brief an input source to read from a stream
     */
    class StreamInputSource : public InputSource
    {
    private:
      istream &in;
    public:
      class StreamBinInputStream : public BinInputStream
      {
      private:
	unsigned int curpos;
	istream &in;
      public:
	StreamBinInputStream(istream &in) : in(in)
	{
	  curpos=0;
	}

	virtual unsigned int curPos() const
	{
	  return curpos;
	}

	virtual unsigned int readBytes(XMLByte *const toFill, 
				       const unsigned int maxToRead)
	{
	  in.read((char *)toFill,maxToRead);

	  int n=in.gcount();

	  if (n>=0)
	    curpos+=n;

	  return n;
	}
      };

      StreamInputSource(istream &in) : in(in) {}

      virtual BinInputStream *makeStream () const
      {
	return new StreamBinInputStream(in);
      }
    }; // class StreamInputSource


  private:
    /**
     * a class to initialize and finalize the library
     */
    class XercesLibrary
    {
    public:
      XercesLibrary()
      {
	// Initialize the XML4C2 system
	try
	{
	  XMLPlatformUtils::Initialize();
	}
	catch (const XMLException& toCatch)
	{
	  cerr << "Error during initialization! :\n"
	       << toCatch.getMessage() << endl;
	  throw;
	}
      }
    };


  private:
    XercesLibrary xercesLibrary; // only used to initialize and
				 // finalize the library
    XMLParser<Callback,XercesString,XercesAttributeList> cspParser;
  };

}

#endif
