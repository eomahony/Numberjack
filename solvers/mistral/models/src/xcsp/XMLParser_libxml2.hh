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
#ifndef _XMLParser_libxml2_h_
#define _XMLParser_libxml2_h_

#include <iostream>
#include <stdexcept>
#include <cerrno>
#include <limits>
#include <libxml/parser.h>

#include <mistral_glo.h>

#include "XMLParser.hh"



namespace CSPXMLParser
{
  using namespace std;

  /**
   * A class to represent a UTF8 encoded string passed by the XML
   * parser. This string is normally read-only and the memory is not
   * owned by this class. The end of the string is represented
   * either by the first NUL character, or by an end pointer
   * (whichever comes first). The end pointer is used to represent
   * substring without having to copy the substring.
   *
   * The append() method provides a limited interface to a writeable
   * string. When this method is used, memory is allocated (and owned
   * by this class) and is automatically resized when necessary. This
   * is a kind of copy-on-write mechanism.
   *
   * This class can also represent a null string (in the SQL sense).
   */
  class UTF8String
  {
  public:
    typedef unsigned char Byte;

    static const int npos=NOVAL;

    UTF8String()
    {
      _beg=_end=NULL;
      allocated=0;
    }

    UTF8String(const UTF8String &s)
    {
      _beg=s._beg;
      _end=s._end;
      allocated=0;
    }

    UTF8String(const Byte *b, const Byte *e)
    {
      _beg=b;
      _end=e;
      allocated=0;
    }

    // nul terminated string
    UTF8String(const Byte *s)
    {
      _beg=s;
      _end=NULL;
      allocated=0;
    }

    // nul terminated string
    UTF8String(const char *s)
    {
      _beg=reinterpret_cast<const Byte *>(s);
      _end=NULL;
      allocated=0;
    }

    ~UTF8String()
    {
      if (allocated)
	delete[] _beg;
    }

    bool isNull() const {return _beg==NULL;}

    bool empty() const {return _beg==_end || (_beg && *_beg==0);}

    int byteLength() const
    {
      if (!_end && _beg)
      {
	// identify the end of the string
	_end=_beg;
	while (*_end)
	  ++_end;
      }

      return _end-_beg;
    }

    void clear()
    {
      _end=_beg;
    }


    /**
     * an iterator on characters
     */
    class iterator
    {
    private:
      const Byte *p;
    public:
      iterator() {p=NULL;}

      explicit iterator(const Byte *s) {p=s;}

      int operator *()
      {
	int ch=*p;

	switch(codeLength(ch))
	{
	case 1:
	  break;
	  //return ch;
	case 2:
	  ch&=0x1F;
	  addNextByte(ch);
	  break;
	  //return ch;
	case 3:
	  if ((ch==0xE0 && p[1]<0xA0)
	      || (ch==0xED && p[1]>0x9F))
	    throw runtime_error("invalid UTF8 character");

	  ch &=0x0F;
	  addNextByte(ch);
	  addNextByte(ch);
	  break;
	  //return ch;
	case 4:
	  if ((ch==0xF0 && p[1]<0x90)
	      || (ch==0xF4 && p[1]>0x8F))
	    throw runtime_error("invalid UTF8 character");

	  ch&=0x07;
	  addNextByte(ch);
	  addNextByte(ch);
	  addNextByte(ch);
	  break;
	  //return ch;
	}
	
	return ch;
      }

      iterator &operator ++() // prefix
      {
	p+=codeLength(*p);
	return *this;
      }

      iterator operator ++(int) // postfix
      {
	iterator tmp=*this;
	p+=codeLength(*p);
	return tmp;
      }

      iterator &operator --() // prefix
      {
	const Byte *q=p;

	do
	{
	  --p;
	} while ((*p & 0xC0)==0x80);

	if (p+codeLength(*p)!=q)
	  throw runtime_error("invalid UTF8 sequence");

	return *this;
      }

      iterator operator --(int) // postfix
      {
	iterator tmp=*this;
	--(*this);
	return tmp;
      }

      iterator &operator =(iterator it)
      {
	p=it.p;
	return *this;
      }

      bool operator !=(iterator it)
      {
	return p!=it.p;
      }

      bool operator ==(iterator it)
      {
	return p==it.p;
      }

      inline bool isWhiteSpace() const
      {
	return *p<128 && (*p==' ' || *p=='\n' || *p=='\r' || *p=='\t' 
			  || *p== '\v' || *p== '\f');
      }

      const Byte *getPointer() const {return p;}

      Byte firstByte() const {return *p;}

    protected:

      /**
       * return the number of bytes of the current code point
       */
      inline int codeLength(int ch)
      {
	if (ch<0x80)
	  return 1; // only one byte
	else
	  if (ch<0xC2)
	    throw runtime_error("invalid UTF8 character");
	  else
	    if (ch<0xD0)
	      return 2; // 2 bytes
	    else
	      if (ch<0xF0)
		return 3; // 3 bytes
	      else
		if (ch<0xF5)
		  return 4; // 4 bytes
		else
		  throw runtime_error("invalid UTF8 character"); 
      }

      inline void addNextByte(int &ch)
      {
	ch<<=6;
	++p;
	if (*p<0x80 || *p>=0xC0)
	  throw runtime_error("invalid UTF8 character");
	ch|=*p & 0x3F;
      }

    };

    iterator begin() const {return iterator(_beg);}
    iterator end() const 
    {
      if (!_end && _beg)
      {
	// identify the end of the string
	_end=_beg;
	while (*_end)
	  ++_end;
      }
      return iterator(_end);
    }



    void append(int ch)
    {
      if (_end+4>_beg+allocated)
	resize();

      // when it's allocated, it's writeable
      Byte *p=const_cast<Byte *>(_end);
      write(p,ch);
      _end=p;
    }


    void append(UTF8String s)
    {
      for(iterator it=s.begin();it!=s.end() && *it;++it)
	append(*it);
    }
    
    /**
     * returns true iff the string contains only white space
     */
    bool isWhiteSpace() const
    {
      iterator it(_beg),end(_end);

      while(it!=end && it.isWhiteSpace())
	++it;

      return it==end || it.firstByte()==0;
    }

    int firstChar() const
    {
      return *iterator(_beg);
    }

    // return the position of sub in this string or npos if not found
    int find(UTF8String sub) const
    {
      int pos=0;
      int firstChar=sub.firstChar();

      for(iterator it(_beg),end(_end);it!=end && it.firstByte();++it)
      {
	if (*it!=firstChar)
	  ++pos;
	else
	{
	  // try to find a match
	  iterator itSub(sub._beg),endSub(sub._end),it2(it);

	  ++itSub;
	  ++it2;

	  while(it2!=end && it2.firstByte() && 
		itSub!=endSub && itSub.firstByte() && 
		*it2==*itSub)
	  {
	    ++it;
	    ++itSub;
	  }

	  if (itSub==endSub||itSub.firstByte()==0)
	    return pos; // found it

	  if (it2==end||it2.firstByte()==0)
	    return npos; // can't find it any more
	}
      }

      return npos; // not found 
    }

    UTF8String substr(int pos, int count=npos)
    {
      const Byte *beg;
      iterator it(_beg),end(_end);

      for(int i=0;i<pos && it!=end && it.firstByte();++i)
	++it;

      beg=it.getPointer();

      for(int i=0;i<count && it!=end && it.firstByte();++i)
	++it;

      return UTF8String(beg,it.getPointer());
    }

    UTF8String substr(iterator beg, iterator end)
    {
      return UTF8String(beg.getPointer(),end.getPointer());
    }

    bool operator ==(const UTF8String s) const
    {
      const Byte *p=_beg,*q=s._beg;

      while(p!=_end && q!=s._end && *p && *q && *p==*q)
      {
	++p;
	++q;
      }
     
      return (*p==0 || p==_end) && (*q==0 || q==s._end);
    }

    bool operator !=(const UTF8String s) const
    {
      return !operator==(s);
    }

    bool operator <(const UTF8String s) const
    {
      const Byte *p=_beg,*q=s._beg;

      while(p!=_end && q!=s._end && *p && *q && *p==*q)
      {
	++p;
	++q;
      }

      return (q!=s._end && *q!=0 ) && ((p==_end || *p==0) || *p<*q);
    }

    /**
     * returns true iff value was successfully transfered in the variable
     */
    template<typename T>
    bool to(T &v) const
    {
      throw runtime_error("don't know how to convert to this type");
    }

    void appendTo(string &v) const
    {
      // fill v with the UTF8 encoding

      for(const Byte *p=_beg;p!=_end && *p;++p)
	v+=*p;
    }

    friend ostream &operator <<(ostream &f, const UTF8String s);

    class Tokenizer
    {
    private:
      iterator it,end;
      vector<int> separators;
    public:
      Tokenizer(const UTF8String s): it(s._beg),end(s._end)
      {
	skipWhiteSpace();
      }

      /**
       * Character ch will be returned as one token
       */ 
      void addSeparator(int ch)
      {
	separators.push_back(ch);
      }

      bool hasMoreTokens()
      {
	return it!=end && *it;
      }

      UTF8String nextToken()
      {
	const Byte *b=it.getPointer(),*e;

	if (it==end || it.firstByte()==0)
	  return UTF8String();

	if (isSeparator(*it))
	  ++it;
	else
	  while(it!=end && *it && !it.isWhiteSpace() && !isSeparator(*it))
	    ++it;

	e=it.getPointer();

	skipWhiteSpace();

	return UTF8String(b,e);
      }

    protected:
      inline bool isSeparator(int ch)
      {
	for(vector<int>::const_iterator it=separators.begin();
	    it!=separators.end();++it)
	  if (*it==ch)
	    return true;

	return false;
      }

      inline void skipWhiteSpace()
      {
	while(it!=end && *it && it.isWhiteSpace())
	  ++it;
      }
    };

  protected:
    void resize()
    {
      if (allocated)
      {
	const Byte *q;
	Byte *tmp,*p;

	tmp=new Byte[2*allocated];

	for(p=tmp,q=_beg;q!=_end;++p,++q)
	  *p=*q;

	delete[] _beg;
	_beg=tmp;
	_end=p;
	allocated*=2;
      }
      else
      {
	allocated=64;
	_beg=_end=new Byte[allocated];
      }
    }

    void write(Byte *&p, int ch)
    {
      if (ch<0x80)
	*p++=ch;
      else
	if (ch<0x800)
	{
	  *p++=0xC0|(ch>>6);
	  *p++=0x80|(ch&0x3F);
	}
	else
	  if (ch<0x10000)
	  {
	    if (ch>=0xD800 && ch<=0xDFFF)
	      throw runtime_error("invalid UTF8 character");

	    *p++=0xE0|(ch>>12);
	    *p++=0x80|((ch>>6)&0x3F);
	    *p++=0x80|(ch&0x3F);
	  }
	  else
	    if (ch<0x110000)
	    {
	      *p++=0xF0|(ch>>18);
	      *p++=0x80|((ch>>12)&0x3F);
	      *p++=0x80|((ch>>6)&0x3F);
	      *p++=0x80|(ch&0x3F);
	    }
	    else
	      throw runtime_error("invalid UTF8 character");
    }

  protected:
    const Byte *_beg;
    mutable const Byte *_end;
    int allocated; // size of allocated array (or 0 if don't own the memory) 
  };

  ostream &operator <<(ostream &f, const UTF8String s)
  {
    // directly output UTF8
    f.write(reinterpret_cast<const char *>(s._beg),s.byteLength());
    return f;
  }

  template<>
  bool UTF8String::to<string>(string &v) const
  {
    // fill v with the UTF8 encoding
    v.clear();

    for(const Byte *p=_beg;p!=_end && *p;++p)
      v+=*p;

    return true;
  }

  template<>
  bool UTF8String::to<int>(int &v) const
  {
    iterator it(_beg),end(_end);
    bool neg=false;

    while(it!=end && it.isWhiteSpace())
      ++it;

    if (it==end || it.firstByte()==0) // end of string?
      return false;

    if (*it=='+')
      ++it;
    else
      if (*it=='-')
      {
	++it;
	neg=true;
      }

    v=0;
    while(it!=end && *it>='0' && *it<='9')
    {
      //??? #warning "must check for overflows" 
      v=v*10+(*it-'0');
      ++it;
    }

    if (neg)
      v=-v;

    while(it!=end && it.isWhiteSpace())
      ++it;

    return it==end || it.firstByte()==0;
  }


  /**
   * represents the attribute list of a XML tag
   */
  class AttributeList
  {
  public:
    typedef unsigned char Byte;

    /**
     * an empty list of attributes
     */
    AttributeList()
    {
      n=0;
      list=NULL;
    }

    AttributeList(const Byte **attr)
    {
      list=attr;

      n=0;
      if (list==NULL)
	return;

      while(list[n]!=NULL)
	n+=2;

      n/=2;
    }

    inline int size() const
    {
      return n;
    }

    UTF8String operator[](const char *name) const
    {
      for(int i=0;i<n;++i)
	if (xmlStrEqual(list[2*i],reinterpret_cast<const Byte *>(name)))
	  return UTF8String(list[2*i+1]);

      return UTF8String();
    }

    inline UTF8String getName(int i) const
    {
      return UTF8String(list[2*i]);
    }

    inline UTF8String getValue(int i) const
    {
      return UTF8String(list[2*i+1]);
    }

  private:
    int n; // number of attributes
    const Byte **list; // list[2*i] is the name of the i-th attribute,
    // list[2*i+1] is its value
  };


  /**
   * @brief the parser using the libxml2 library
   */
  template<class Callback=CSPParserCallback>
  class XMLParser_libxml2
  {
  public:
    XMLParser_libxml2(Callback &callback): cspParser(callback)
    {
      LIBXML_TEST_VERSION
	}

    int parse(istream &in)
    {
      /**
       * We don't use the DOM interface because it reads the document as
       * a whole and it is too memory consuming. The TextReader
       * interface is better because it reads only one node at a
       * time,but the text of a node is read as a whole and this can
       * still be large for the definition of some
       * relations. Therefore, we use the SAX interface.
       *
       * We also use the push mode to be able to read from any C++
       * stream.
       */
      const char *filename=NULL; // name of the input file
      xmlSAXHandler handler;
      xmlParserCtxtPtr parserCtxt=NULL;
    
      const int bufSize=4096;
      char *buffer=new char[bufSize];

      int size;

      xmlSAXVersion(&handler,1); // use SAX1 for now ???

      handler.startDocument=startDocument;
      handler.endDocument=endDocument;
      handler.characters=characters;
      handler.startElement=startElement;
      handler.endElement=endElement;
      handler.comment=comment;

      try
      {
	xmlSubstituteEntitiesDefault(1);

	in.read(buffer,bufSize);
	size=in.gcount();

	if (size>0)
	{
	  parserCtxt=xmlCreatePushParserCtxt(&handler,&cspParser,
					     buffer,size,filename);

	  while(in.good())
	  {
	    in.read(buffer,bufSize);
	    size=in.gcount();

	    if (size>0)
	      xmlParseChunk(parserCtxt,buffer,size,0);
	  }

	  xmlParseChunk(parserCtxt,buffer,0,1);

	  xmlFreeParserCtxt(parserCtxt);

	  xmlCleanupParser();
	}
      }
      catch(runtime_error &e)
      {
	// ???
	cout << "Exception at line " << parserCtxt->input->line << endl;
	throw(e);
      }

      delete[] buffer;

      return 0;
    }

    int parse(const char *filename)
    {
      ifstream in(filename);
      return parse(in);
    }

    void setPreferredExpressionRepresentation(Syntax syntax)
    {
      cspParser.setPreferredExpressionRepresentation(syntax);
    }
  
  protected:
    typedef XMLParser<Callback,UTF8String,AttributeList> Parser;

    /*************************************************************************
     *
     * SAX Handler
     *
     *************************************************************************/

    static void comment(void *parser, const xmlChar * value)
    {
    }

    static void startDocument(void *parser)
    {
#ifdef debug
      cout << "Parsing begins" << endl;
#endif
      static_cast<Parser *>(parser)->startDocument();
    }

    static void endDocument(void *parser)
    {
#ifdef debug
      cout << "Parsing ends" << endl;
#endif
      static_cast<Parser *>(parser)->endDocument();
    }

    static void characters(void *parser,const xmlChar *ch,int len)
    {
#ifdef debug
      cout << "    chars '" << UTF8String(ch,ch+len) << "'" << endl;
#endif
      static_cast<Parser *>(parser)->characters(UTF8String(ch,ch+len));
    }

    static void startElement(void * parser, 
			     const xmlChar *name, 
			     const xmlChar **attr)
    {
      AttributeList attributes(attr);

#ifdef debug
      cout << "  begin element " << UTF8String(name) << endl;
      for(int i=0;i<attributes.size();++i)
      {
	cout << "    attribute " << attributes.getName(i)
	     << " = " << attributes.getValue(i) << endl;
      }
#endif
      
      static_cast<Parser *>(parser)->startElement(UTF8String(name),
						  attributes);
    }

    static void endElement(void *parser, const xmlChar *name) 
    {
#ifdef debug
      cout << "  end element " << UTF8String(name) << endl;
#endif
      static_cast<Parser *>(parser)->endElement(UTF8String(name));
    }

  protected:
    Parser cspParser;
  };

}

#endif
