
/*
  Mistral is a constraint satisfaction and optimisation library
  Copyright (C) 2003-2005  Emmanuel Hebrard
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

  The author can be contacted electronically at 
  ehebrard@cse.unsw.edu.au.
*/

/** \file set.h
    \brief Header for the class MistralSet.
*/

//-------------Representation of sets as vectors of bits--------------------------
//
//    the set has to be bounded in size (N)
//
//--------------------------------------------------------------------------------
// Emmanuel Hebrard 29st Mar 2005


#ifndef _MistralSet_h
#define _MistralSet_h

#include <stdlib.h> 
#include <iostream>
#include <sstream>
#include <cstring>
#include <mistral_glo.h>

//using namespace std;
typedef unsigned int uint;
void showUint(uint n, std::ostream & os);


namespace Mistral {

  /**********************************************
   * Vector
   **********************************************/

  /// An extendable array.
  /**
     Could be represented as a vector, though
     this class is simpler and faster.
  */
  template <class T>
  class Vector {
  public:

    //static const int MIN_CAPACITY = 0;

    T *   stack_;
    int capacity;
    int    size;
  
    Vector()
    {

      //std::cout << "create an empty vector" << std::endl;

      capacity = 0;
      size = 0;
      stack_ = NULL;
    }

    Vector( const Vector<T>& s )
    {
      capacity = s.size;
      stack_ = (T*) malloc(capacity*sizeof(T));
      for(size = 0; size<capacity; ++size)
	stack_[size] = s[size];
    }

    Vector(const int n)
    {
      capacity = n;
      size = 0;
      if( capacity )
	stack_ = (T*) malloc(capacity*sizeof(T));
      else stack_ = NULL;
    }
  
    virtual ~Vector()
    {

      //std::cout << "delete vector" << std::endl;

      free( stack_ );
    }


    void init(const int s, const int c)
    {

      //std::cout << "init vector" << std::endl;

      size = s;
      if( capacity < c ) {
	capacity = c;
	T* new_stack = (T*) realloc(stack_, capacity*sizeof(T));
	stack_ = new_stack;
      }
    }

    void extendStack()
    {
      capacity = ((capacity+1) << 1);
      T* new_stack = (T*) realloc(stack_, capacity*sizeof(T));
      stack_ = new_stack;
    }


    void extendStack( int l )
    {
      capacity += l;
      T* new_stack = (T*) realloc(stack_, capacity*sizeof(T));
      stack_ = new_stack;
    }

    void resize( int l )
    {

      //std::cout << "resize vector" << std::endl;

      if( capacity < l ) {
	capacity = l;
	T* new_stack = (T*) realloc(stack_, capacity*sizeof(T));
	stack_ = new_stack;
      }
      size = l;
    }

    inline int empty() const
    {
      return !size;
    }
  
    inline void push(T x)
    {
      if( capacity == size ) 
	extendStack();
      stack_[size++] = x;
    }
    inline void push_back(T x)
    {
      if( capacity == size )
	extendStack();
      stack_[size++] = x;
    }
    inline T popUntil( const int level )
    {
      size = level;
      return stack_[size];
    }
    inline T pop()
    {
      return stack_[--size];
    }
    inline void pop(T& x)
    {
      x = stack_[--size];
    }
    inline void clear()
    {
      size = 0;
    }

    inline void remove(T& elt)
    {
      int j=size;
      while(j && stack_[--j] != elt);
      //{std::cout << j << " " << stack_[j] << " " << elt << std::endl;};
      //std::cout << j << std::endl;
      //assert( stack_[j] == elt );

      stack_[j] = stack_[--size];
    }
  
    inline void erase(const int i)
    {  
      stack_[i] = stack_[--size];
    }

    inline void setBack(const T& x)
    {
      stack_[size-1] = x;
    }

    inline T& back()
    {
      return stack_[size-1];
    }

    inline const T back() const
    {
      return stack_[size-1];
    }

    //   inline const T out() const
    //     {
    //       return stack_[size];
    //     }
    inline T& operator[](const int i)
    {
      return stack_[i];
    }
    inline T operator[](const int i) const
    {
      return stack_[i];
    }

    //   Vector<T>& operator=( const Vector<T>& s )
    //     {
    //       std::cout << 11 << std::endl;
    //       if( capacity < s.capacity ) {
    // 	capacity = s.capacity;
    // 	T *new_stack_ = (T*) malloc(capacity*sizeof(T));
    // 	free( stack_ );
    // 	stack_ = new_stack_;
    //       }
    //       for(size = 0; size<s.size; ++size)
    // 	stack_[size] = s[size];
    //       std::cout << 22 << std::endl;
    //       return *this;
    //     }

    void print(std::ostream& o) const
    {
      for(int i=0; i<size; ++i)
	o << stack_[i] << " ";
      o << "(" << size << ")" << std::endl;
    }

  };


  template <class T>
  class Array {

  public: 
    int size;
    T data[0];

    Array(const Vector<T>& ps) 
    {
      size = ps.size;
      for (int i=0; i<ps.size; ++i) 
	data[i] = ps[i];
    }

    static Array<T>* Array_new(const Vector<T>& ps)
    {
      void* mem = malloc(sizeof(Array<T>) + sizeof(T)*(ps.size));
      return new (mem) Array<T>(ps); 
    }

//     virtual ~Array() {
      
//       std::cout << ((void*)(data)) << std::endl;

//       free((void*)(data));
//     }

    inline T& operator [] (int i) { return data[i]; }
    inline T operator [] (int i) const { return data[i]; }
  };

  template <class T>
  class MistralNode {
  public:

    int index;

    T elt;  
    MistralNode< T > *next;
    MistralNode< T > *pred;

    MistralNode( )
    {
      index = -1;

      next = NULL;
      pred = NULL;
      elt = NULL;
    }

    ~MistralNode() {}

    inline void erase( )
    {
      if( pred ) {
	pred->next = next;
	if( next )
	  next->pred = pred;
	pred = NULL;
      }
    }

    inline bool isIn()
    {
      return pred;
    }

  };


  template <class T>
  inline MistralNode< T >* nextNode( MistralNode< T >*& nd )
  {
    do nd = nd->next;
    while( nd && !nd->elt );
    return nd;
  }





  template <class T>
  class MistralList {
  public:
    MistralNode< T > head;

    MistralList( ) {}
    ~MistralList( ) {}

    void checkConsistency()
    {
      MistralNode< T > *nd;
      nd = &head;
      //     while( nextNode(nd) ) 
      //       assert( nd->pred->next == nd );
    }

    inline void insert( MistralNode< T >& x )
    {
      if( !x.pred ) {
	x.next = head.next;
	x.pred = &head;
	if( head.next ) head.next->pred = &x;
	head.next = &x;
      }
    }

    inline void clear( )
    {
      head.next = NULL;
    }

    void print(std::ostream & os) 
    {
      MistralNode<T> *nd = &head;
      while( nextNode(nd) ) {
	os << " <";
	nd->elt->printshort( os );
	os << ">";
      }
      os << std::endl;
    }

  };


  template<typename T, typename F> T StrmConvert( F from )
  {
    std::stringstream temp;
    temp << from;
    T to = T();
    temp >to;
    return to;
  }

  template<typename F> std::string StrmConvert( F from )
  {
    return StrmConvert<std::string>( from );
  }

  inline std::string int2string(int x)
  {
    std::ostringstream o;
    o << x;
    return o.str();
  }


  /**********************************************
   * Set
   **********************************************/

  ///   A representation of sets using a vector of bits.
  /*!
    The sets are bounded. Beware that the space complexity
    depends only on the bound, and not on the current
    number of elements. As 32 elements can be 
    handled at once, the usual set operations (union,
    inetersectiuon, difference...) are very efficient.
    However iterating on the elements or counting
    them may not be so.
  */
  template< class T >
  class Bitset
  { 
  public:

    static const uint idx1 = 0xff;
    static const uint idx2 = 0xff00;
    static const uint idx3 = 0xff0000;
    static const uint idx4 = 0xff000000;

    /// 32 bits words
    static const uint EXP = 5;

    /*   /// 64 bits words */
    /*   static const uint EXP = 6; */

    static const uint size_word_bit = (1 << EXP);
    static const uint size_word_byte = (size_word_bit >> 3);
    static const uint CACHE = (size_word_bit - 1);
    static const uint LASTCHAR = (size_word_bit - 8);

  public:
    static const uint empt = 0;
    static const uint full = ~0;//0xffffffff;

    /// A vector of bits
    T* table;
    /// Number of 4-bytes words used to represent the set
    int pos_words;
    int neg_words;
    /*!
      Return the cardinality of the ith word of the set [O(N/32)].
    */
    inline  uint size( const int i )const
    {  
      int thecount = 0, tab;
      if( (tab = table[i]) ) {
	tab = tab - ((tab >> 1) & 0x55555555);
	tab = (tab & 0x33333333) + ((tab >> 2) & 0x33333333);
	tab = (tab + (tab >> 4)) & 0x0F0F0F0F;
	thecount += (tab*0x01010101 >> 24);
      }
      return thecount;  
    }



    inline int random() const
    {
      // get a random element of the bitset    
      int j = pos_words - neg_words;
      //int idx = (rand() % j);
      int idx = (randint(j));

      while( --j && !table[neg_words+idx] )
	//idx = (idx + (rand() % j)) % (pos_words - neg_words);
	idx = (idx + (randint(j))) % (pos_words - neg_words);

      idx += neg_words;

      int bs = table[idx];
      int vmn = 0;
      int vmx = size_word_bit-1;
      int side[2], dir, mask = 16;
    
      while(mask)
	{
	  side[0] = bs & mask_inf[(vmn+vmx)/2+1];
	  side[1] = bs & mask_sup[(vmn+vmx)/2+1];

	  // are there values on the right?
	  if(!(side[1])) {
	    dir = 0;
	  }
	  else {
	    // are there values on the left?
	    if(!(side[0])) {
	      //std::cout << "go right because no values on the left" << std::endl; 
	      dir = 1;
	    }
	    // if there are values both ways, we select randomly
	    else {
	      //dir=rand()%2;	    
	      dir=randint(2);
	      //std::cout << "go " << ( dir ? "right" : "left") << " at random" << std::endl;
	    }
	  }
	  if(dir) {
	    // go right
	    vmn += mask;
	    mask /= 2; 
	  } else {
	    // go left
	    vmx -= mask;
	    mask /= 2; 
	  }
	  bs = side[dir];
	}

      return idx*size_word_bit+vmn;
    }

    /*!
      Minimum element in the set [O(N/8)]. 
    */
    inline  int min()const
    { 
      uint tab, j;
      int i=neg_words, k;
      

      while( i < pos_words ) {
	if( (tab = table[i]) )
	  for( j=0; j<size_word_byte; ++j ) {
	    if( (k = getfirst[tab & idx1]) >= 0 )
	      return ( (i<<EXP)+(j<<3)+k ); 
	    tab = (tab >> 8);
	  }
	++i;
      }
      return NOVAL;
    }

  
    /*!
      Maximum element in the set [O(N/32)]. (however, not very efficient)
    */
    inline int max() const
    { 
      uint tab;
      int i=pos_words, j, k;
    
      while( i-- > neg_words )
	if( (tab = table[i]) ) {
	  j = size_word_byte;
	  while( j-- ) {
	    if( (k = getlast[(tab & idx4) >> LASTCHAR]) >= 0 ) 
	      return ( (i<<EXP)+(j<<3)+k );	
	    tab = (tab << 8);
	  }
	}
      return NOVAL;
    }

    Bitset()
    {
      init();
    }

    void init()
    {
      pos_words = 0;
      neg_words = 0;
      table = NULL;
    }


    Bitset(const int sz, const int* elt) 
    {
      int lb = MAXINT;
      int ub = MININT;
      for(int i=0; i<sz; ++i) {
	if(elt[i] > ub) ub = elt[i];
	if(elt[i] < lb) lb = elt[i];
      }

      init(lb,ub,empt);

      for(int i=0; i<sz; ++i) 
	insert( elt[i] );
    }

    Bitset(const int lb, const int ub, const uint p) 
    {
      init(lb,ub,p);
    }



    void reinit(const int lb, const int ub, const uint p) 
    {
      table += neg_words;
      delete [] table;
      init(lb, ub, p);
    }


    void init(const int sz, const uint p) 
    {
      pos_words = sz;
      neg_words = 0;

      if( sz ) {
	table = new T[pos_words];
	for(int i=0; i<pos_words; ++i) 
	  table[i]=p;
      } else table = NULL;
    }

    void init(const int lb, const int ub, const uint p) 
    {
      neg_words = (lb >> EXP);
      pos_words = (ub >> EXP)+1;

      table = new T[pos_words-neg_words];
      for(int i=0; i<pos_words-neg_words; ++i) 
	table[i]=p;
      table[pos_words-neg_words-1] &= 
	(p >> (size_word_bit-1-(ub & CACHE)));
      table[0] &= (p << (lb & CACHE));
      table -= neg_words;

      //std::cout << table << std::endl;

    }

    void init(Bitset<unsigned int>& s) 
    {
      //std::cout << "init " << s.table 

      pos_words = s.pos_words;
      neg_words = s.neg_words;
  
      table = new T[pos_words-neg_words];
      table -= neg_words;
      for(int i=neg_words; i<pos_words; ++i) {
	//std::cout << i << "  *" << (s.table+i) << std::endl;
	table[i].init(s.table+i, s.size(i));
      }
    }


    Bitset(const Bitset<T>& s) 
    {
      clone( s );
    }

    void clone(const Bitset<T>& s)
    {
      neg_words = s.neg_words;
      pos_words = s.pos_words;
      table = new T[pos_words-neg_words];
      memcpy(table, s.table+neg_words,
	     size_word_byte*(pos_words-neg_words));
      table -= neg_words;
    }

    void pointTo(Bitset<T>& s)
    {
      neg_words = s.neg_words;
      pos_words = s.pos_words;
      table = s.table;
    }

    void pointTo(unsigned int *t)
    {
      neg_words = 0;
      pos_words = 1;
      table = t;
    }

    void copy(const Bitset<T>& s) 
    {
      int i = (pos_words > s.pos_words ? s.pos_words : pos_words);
      int k, j = (neg_words < s.neg_words ? s.neg_words : neg_words);
      for( k=neg_words; k<j; ++k)
	table[k] = empt;
      for( k=i; k<pos_words; ++k)
	table[k] = empt;
      if( i>j )
	memcpy(table+j,s.table+j,size_word_byte*(i-j));
    }

    // void copyTo(Bitset<T>& s) const
    // {
    //   int i = (pos_words > s.pos_words ? s.pos_words : pos_words);
    //   int j = (neg_words < s.neg_words ? s.neg_words : neg_words);
    //   //s.clear();
  
    //   if( i>j )
    //     memcpy(s.table+j,table+j,size_word_byte*(i-j));
    // }

    virtual ~Bitset() 
    {
      table += neg_words;
      delete [] table; 
    }

    void destroy() 
    {
      table += neg_words;
      neg_words = 0;
      delete [] table; 
      table = NULL;
    }

    bool isBuilt()
    {
      return (table != NULL);
    }


    inline void swap(Bitset<T>& s)
    {
      T *aux = s.table;
      s.table = table;
      table = aux;
    }

    inline void  erase(const int elt) 
    {
      int i = (elt >> EXP);
      if( (i >= neg_words) && 
	  (i <  pos_words) )
	table[i] &= (full ^ (1 << (elt & CACHE)));
    }

    inline void  fastErase(const int elt) 
    {
      table[(elt >> EXP)] &= (full ^ (1 << (elt & CACHE)));
    }

    inline void  wordErase(const int elt) 
    {
      table[neg_words] &= (full ^ (1 << (elt & CACHE)));
    }

    inline int next(const int elt) const {

      uint j, tab;
      int e, k, i = ((1+elt) >> EXP);
      if( i < pos_words ) {
	e = ((1+elt) & CACHE);
	j = (e >> 3);

	if( (tab = ((table[i] & (full << e)) >> (j<<3))) )
	  while( j<size_word_byte ) {
	    if( (k = getfirst[tab & idx1]) >= 0 )
	      return ( (i<<EXP)+(j<<3)+k );
	    tab = (tab >> 8);
	    ++j;
	  }
	while( ++i < pos_words ) 
	  if( (tab = table[i]) ) 
	    for(j=0; j<size_word_byte; ++j) {
	      if( (k = getfirst[tab & idx1]) >= 0 )
		return ( (i<<EXP)+(j<<3)+k );
	      tab = (tab >> 8);
	    }
      }
      return NOVAL;
    }


    inline int prev(const int elt) const {
      uint tab;
      int i = ((elt-1) >> EXP);

      if( i >= neg_words ) {
	int e = ((elt-1) & CACHE), k;
	int j = 1+(e >> 3);

	if( (tab = ((table[i] & (full >> (CACHE - e))) << ((4-j) << 3))) ) 
	  while( j-- ) {
	    if( (k = getlast[(tab & idx4) >> LASTCHAR]) >= 0 )
	      return ( (i<<EXP)+(j<<3)+k );
	    tab = (tab << 8);
	  }
	while( i-- > neg_words ) 
	  if( (tab = table[i]) ) {
	    j = size_word_byte;
	    while( j-- ) {
	      if( (k = getlast[(tab & idx4) >> LASTCHAR]) >= 0 )
		return ( (i<<EXP)+(j<<3)+k );
	      tab = (tab << 8);
	    }
	  }
      }
      return NOVAL;
    }



    inline void xorTo(const Bitset<T>& s) 
    {
      int i = (pos_words > s.pos_words ? s.pos_words : pos_words);
      int j = (neg_words < s.neg_words ? s.neg_words : neg_words);
      while( i-- > j )
	s.table[i] ^= table[i];
    }

    inline void fastXorTo(const Bitset<T>& s) 
    {
      int i = pos_words;
      while( i-- > neg_words )
	s.table[i] ^= table[i];
    }

    inline void xorWith(const Bitset<T>& s) 
    {
      int i = (pos_words > s.pos_words ? s.pos_words : pos_words);
      int j = (neg_words < s.neg_words ? s.neg_words : neg_words);
      while( i-- > j )
	table[i] ^= s.table[i];
    }

    inline void fastXorWith(const Bitset<T>& s) 
    {
      int i = pos_words;
      while( i-- > neg_words )
	table[i] ^= s.table[i];
    }

    inline void unionWith(const Bitset<T>& s) 
    {
      int i = (pos_words > s.pos_words ? s.pos_words : pos_words);
      int j = (neg_words < s.neg_words ? s.neg_words : neg_words);
      while( i-- > j )
	table[i] |= s.table[i];
    }

    inline void unionTo(Bitset<T>& s) const 
    {
      s.unionWith( *this );
    }

    inline void intersectWith(const Bitset<unsigned int>& s) 
    {
      int i = (pos_words > s.pos_words ? s.pos_words : pos_words);
      int j = (neg_words < s.neg_words ? s.neg_words : neg_words);
      int k = pos_words;
      while( k > i ) {
	--k;
	table[k] = empt;
      }
      while( k > j ) {
	--k;
	table[k] &= s.table[k];
      }
      while( k-- > neg_words )
	table[k] = empt;
    }

    // inline void intersectWithInt(const Bitset<unsigned int>& s) 
    // {
    //   int i = (pos_words > s.pos_words ? s.pos_words : pos_words);
    //   int j = (neg_words < s.neg_words ? s.neg_words : neg_words);
    //   int k = pos_words;
    //   while( k > i ) {
    //     --k;
    //     table[k] = empt;
    //   }
    //   while( k > j ) {
    //     --k;
    //     table[k] &= s.table[k];
    //   }
    //   while( k-- > neg_words )
    //     table[k] = empt;
    // }

    inline void intersectTo(Bitset<T>& s) const
    {
      s.intersectWith( *this );
    }

    inline void setminusWith (const Bitset<unsigned int>& s) 
    {
      int i = (pos_words > s.pos_words ? s.pos_words : pos_words);
      int j = (neg_words < s.neg_words ? s.neg_words : neg_words);
      while( i-- > j )
	table[i] &= (~(s.table[i]));
    }

    inline void setminusTo (Bitset<T>& s) const
    {
      s.setminusWith( *this );
    }

    inline void xorWith(const Bitset<T>* s) 
    {
      xorWith(*s);
    }

    inline void unionWith  (const Bitset<T>* s) 
    {
      unionWith(*s);
    }

    inline void intersectWith(const Bitset<T>* s) 
    {
      intersectWith(*s);
    }

    inline void setminusWith (const Bitset<T>* s) 
    {
      setminusWith(*s);
    }
    inline void unionTo  (Bitset<T>* s) const
    {
      s->unionWith(*this);
    }

    inline void intersectTo(Bitset<T>* s) const
    {
      s->intersectWith(*this);
    }

    inline void setminusTo (Bitset<T>* s) const
    {
      s->setminusWith(*this);
    }

    inline bool equal(const Bitset<T>& s) const 
    {

      //   print( std::cout );
      //   std::cout << " =?= " ;
      //   s.print( std::cout );
      //   std::cout << std::endl;

      int i=pos_words;
      int j=s.pos_words;
      int k;
  

      //   if( i==3 )
      //     std::cout << "ahah" << std::endl;
      //   if( j==3 )
      //     std::cout << "ahah" << std::endl;


      while( j > i )
	if( s.table[--j] ) return false; 
      while( i > j )
	if( table[--i] ) return false;

      j=neg_words;
      k=s.neg_words;

      //   if( j==3 )
      //     std::cout << "ahah" << std::endl;
      //   if( k==3 )
      //     std::cout << "ahah" << std::endl;

      while( j > k )
	if( s.table[k++] ) return false;
      while( k > j )
	if( table[j++] ) return false;

      while( i-- > j )
	if( table[i] != s.table[i] ) return false;
  
      return true;
    }

    inline bool included(const Bitset<T>& s) const 
    {
      int i = (pos_words > s.pos_words ? s.pos_words : pos_words);
      int j = (neg_words < s.neg_words ? s.neg_words : neg_words);
      int k = pos_words;
      while( k > i ) {
	--k;
	if( table[k] ) return false;
      }
      while( k > j ) {
	--k;
	if( table[k] != (table[k] & s.table[k]) ) return false;
      }
      while( k-- > neg_words )
	if( table[k] ) return false;
      return true;
    }

    inline bool included(const Bitset<T>& s, int ww) const 
    {

      //std::cout << 11 << std::endl;

      int i = (pos_words > s.pos_words ? s.pos_words : pos_words);
      int j = (neg_words < s.neg_words ? s.neg_words : neg_words);
      int k = pos_words;
      while( k > i ) {
	--k;
	if( table[k] ) {
	  //std::cout << "false" << std::endl;
	  return false;
	}
      }
      while( k > j ) {
	--k;
	if( table[k] != (table[k] & s.table[k]) ) {
	  //std::cout << "false" << std::endl;
	  return false;
	}
      }
      while( k-- > neg_words )
	if( table[k] ) {
	  //std::cout << "false" << std::endl;
	  return false;
	}

      //std::cout << "true" << std::endl;

      return true;
    }

    inline bool included(const Bitset<T>* s) const 
    {
      return included( *s );
    }

    inline bool intersect(const Bitset<T>* s) const 
    {
      return intersect( *s );
    }


	
    /**
     * Returns the number of bits set in val.
     * For a derivation of this algorithm, see
     * "Algorithms and data structures with applications to 
     *  graphics and geometry", by Jurg Nievergelt and Klaus Hinrichs,
     *  Prentice Hall, 1993.
     */
    inline uint size()const 
    {  
      int i=pos_words;
      uint tab, thecount=0;
      while( i-- > neg_words ) 
	if( (tab = table[i]) ) {
	  tab = tab - ((tab >> 1) & 0x55555555);
	  tab = (tab & 0x33333333) + ((tab >> 2) & 0x33333333);
	  tab = (tab + (tab >> 4)) & 0x0F0F0F0F;
	  thecount += (tab*0x01010101 >> 24);
	}
      return thecount;  
    }

    inline uint wordSize()const 
    {  
      uint tab, thecount=0;
      if( (tab = table[neg_words]) ) {
	tab = tab - ((tab >> 1) & 0x55555555);
	tab = (tab & 0x33333333) + ((tab >> 2) & 0x33333333);
	tab = (tab + (tab >> 4)) & 0x0F0F0F0F;
	thecount += (tab*0x01010101 >> 24);
      }
      return thecount;  
    }

    /*!
      Check if element elt belong to the set [O(1)]
    */
    inline  bool member(const int elt)const 
    {
      int i = (elt >> EXP);
      return ( (i >= neg_words) && 
	       (i <  pos_words) && 
	       (table[i] & (1 << (elt & CACHE))) );
    }

    inline  bool fastMember(const int elt)const 
    {
      return ( (table[(elt >> EXP)] & (1 << (elt & CACHE))) );
    }
    inline  bool wordMember(const int elt)const 
    {
      return ( (table[neg_words] & (1 << (elt & CACHE))) );
    }
    /*!
      Insert element elt into the set [O(1)]
    */

    inline  void insert(const int elt)
    {
      int i = (elt >> EXP);
      if( (i >= neg_words) && 
	  (i <  pos_words) )
	table[i] |= (1 << (elt & CACHE));
    }
    inline  void fastInsert(const int elt)
    {
      table[(elt >> EXP)] |= (1 << (elt & CACHE));
    }
    inline   void wordInsert(const int elt)
    {
      table[neg_words] |= (1 << (elt & CACHE));
    }

    /*!
      Insert element elt into the set or erase it if it is already member [O(1)]
    */

    inline  void invert(const int elt)
    {
      int i = (elt >> EXP);
      if( (i >= neg_words) && 
	  (i <  pos_words) )
	table[i] ^= (1 << (elt & CACHE));
    }

    inline  void fastInvert(const int elt)
    {
      table[(elt >> EXP)] ^= (1 << (elt & CACHE));
    }

    inline  void wordInvert(const int elt)
    {
      table[neg_words] ^= (1 << (elt & CACHE));
    }

    /*!
      Return true iff the set is empty [O(N/32)]
    */
    inline  bool empty()const
    { 
      int i = pos_words;
      while( i-- > neg_words ) 
	if(table[i]) return false;
      return true;  
    }



    /*!
      Return true iff the calling object intersect s [O(N/32)]
    */

    inline bool intersect(const Bitset<T>& s)const
    {
      int i = (pos_words > s.pos_words ? s.pos_words : pos_words);
      int j = (neg_words < s.neg_words ? s.neg_words : neg_words);
      while( i-- > j )
	if(table[i] & s.table[i]) return true;
      return false;
    }

    // bool wordIntersect(const Bitset<T>& s)const;
    inline bool wordIntersect(const Bitset<T>& s) const 
    {
      return ( table[neg_words] & s.table[neg_words] ) ;
    }
    inline  bool fastIntersect(const Bitset<T>& s, int& idx)const
    {
      if( table[idx] & s.table[idx] ) return true; 
      if( pos_words > neg_words ) {
	idx = pos_words;
	while( idx > neg_words ) {
	  --idx;
	  if(table[idx] & s.table[idx]) return true;
	}
      }
      //   int x = pos_words;
      //   while( x > neg_words ) {
      //     --x;
      //     if(table[x] & s.table[x]) return true;
      //   }
      return false;
    }
    /*!
      Return true iff the calling object intersect [lo..up] [O(N/32)]
    */
    inline   bool intersect(const int lb, const int ub)const
    {
      int i = (ub >> EXP);
      int j = (lb >> EXP);

      if( i < neg_words || j >= pos_words ) 
	return false;

      uint masked_lb = (full << (lb & CACHE));
      uint masked_ub = (full >> (CACHE - (ub & CACHE)));

      if( i == j ) {
	if( table[i] & (masked_lb & masked_ub) ) return true;
	else return false;
      }

      if( i >= pos_words )
	i = pos_words-1;
      else if( table[i--] & masked_ub ) return true;
  
      if( j < neg_words ) 
	j = neg_words;
      else if( table[j++] & masked_lb ) return true;

      while( i >= j )
	if(table[i--]) return true;
      return false;
    }

    /*!
      Increment by x all elements in the set.
      Any element greater or equal than the capacity 
      of the set minus x is erased [O(N/32)]
    */
    inline void increment(const int v)
    {
      int step = (v >> EXP); 
      int i = pos_words;
      int e = (v & CACHE);
      int f = 32-e;
      int j = neg_words+step;
      while( --i > j ) 
	table[i] = ((table[i-step] << e) | ((table[i-step-1] & mask_sup[f]) >> f));
      if( i >= neg_words+step ) table[i] = (table[i-step] << e);
      while( i > neg_words ) table[--i] = 0;
    }

    /*!
      Decrement by x all elements in the set.
      Any element lower than x is erased [O(N/32)]
    */
    inline  void decrement(const int v)
    {
      int step = (v >> EXP); 
      int i = neg_words-1;
      int e = (v & CACHE);
      int f = 32-e;
      int j = pos_words-step-1;
      while( ++i < j )
	table[i] = ((table[i+step] >> e) | ((table[i+step+1] & mask_inf[e]) << f));
      if( i < pos_words-step ) table[i] = (table[i+step] >> e);
      while( ++i < pos_words ) table[i] = 0;
    }

    /*!
      Changes every value to its arythmetic negation
    */
    inline void negate( Bitset<T>& s )
    {
      int i = (pos_words > -s.neg_words ? -s.neg_words : pos_words);
      int j = (neg_words < -s.pos_words ? -s.pos_words : neg_words);
      unsigned int x, aux, rest = ( i < pos_words && (table[i] & 1) );
      while( i-- > j ) {
	aux = (table[i] & 1);
	x = (table[i] >> 1);
	x = ((x & 0x55555555) <<  1) | ((x >>  1) & 0x55555555);
	x = ((x & 0x33333333) <<  2) | ((x >>  2) & 0x33333333);
	x = ((x & 0x0F0F0F0F) <<  4) | ((x >>  4) & 0x0F0F0F0F);
	x = (x << 24) | ((x & 0xFF00) << 8) |
	  ((x >> 8) & 0xFF00) | (x >> 24);
	s.table[-i-1] = (x | rest);
	rest = aux;
      }
      if(rest)
	s.table[i+1] |= rest;
    }


    /*!
      Insert all elements between 0 to capacity [O(N/32)]
    */
    inline  void fill()
    {
      int i = pos_words;
      while( i > neg_words )
	table[--i] = full;
    }
    /*!
      Erase all elements [O(N/32)]
    */
    inline  void clear()
    {
      int i = pos_words;
      while( i > neg_words )
	table[--i] = empt;
    }
    /*!
      Erase all elements but v [O(N/32)]
    */
    inline  void setTo( const unsigned int v )
    {
      int i, j = (v >> EXP);
      for(i=neg_words; i<j; ++i)
	table[i] = empt;
      table[j] = (1 << v);
      for(i=j+1; i<pos_words; ++i) 
	table[i] = empt;
    }

    /*!
      flip all elements [O(N/32)]
    */
    inline  void flip()
    {
      int i = pos_words;
      while( i > neg_words )
	table[--i] ^= full;
    }

    /*!
      Erase all elements strictly lower than l [O(N/32)]
    */
    inline void setMin(const int bound)
    {
      int ith_word=(bound >> EXP);
      if( ith_word >= neg_words ) {
	if( ith_word <  pos_words ) {
	  int i=ith_word;
	  while( i-- > neg_words ) table[i]=0;
	  table[ith_word] &= (full << (bound & CACHE));
	} else clear();
      }
    }

    /*!
      Erase all elements strictly greater than u [O(N/32)]
    */
    inline void setMax(const int bound)
    {
      int ith_word=(bound >> EXP);
      if( ith_word <  pos_words ) {
	if( ith_word >= neg_words ) {
	  int i=pos_words;
	  while( --i > ith_word ) table[i]=0;
	  table[ith_word] &= (full >> (CACHE - (bound & CACHE)));
	} else clear();
      }
    }

    /*!
      Erase all elements in the interval [l..u] [O(N/32)]
    */
    inline  void removeInterval(const int lb, const int ub)
    {
      if( lb <= ub ) {
	int lb_word = lb >> EXP;
	int ub_word = ub >> EXP;
    
	lb_word = ( lb_word < neg_words ? neg_words : lb_word );
	ub_word = ( ub_word >= pos_words ? pos_words-1 : ub_word );

	uint masked_lb = 0;
	uint masked_ub = 0;
	if( lb_word >= neg_words ) 
	  masked_lb = (NOVAL >> (CACHE - (lb & CACHE)));
	if( ub_word < pos_words ) 
	  masked_ub = ((full-1) << (ub & CACHE) );
    
	if( lb_word == ub_word ) 
	  table[lb_word] &= (masked_lb | masked_ub);
	else {
	  table[lb_word] &= masked_lb;
	  table[ub_word] &= masked_ub;
	  while( --ub_word > lb_word )
	    table[ub_word] = 0;
	}
      }
    }

    /*!
      Add all elements in the interval [l..u] [O(N/32)]
    */
    inline void addInterval(int lb, int ub)
    {
      if( lb <= ub ) {
	int lb_word = lb >> EXP;
	int ub_word = ub >> EXP;
    
	lb_word = ( lb_word < neg_words ? neg_words : lb_word );
	ub_word = ( ub_word >= pos_words ? pos_words-1 : ub_word );

	uint masked_lb = full;
	uint masked_ub = full;
	if( lb_word >= neg_words ) 
	  masked_lb ^= (NOVAL >> (CACHE - (lb & CACHE)));
	if( ub_word < pos_words ) 
	  masked_ub ^= ((full-1) << (ub & CACHE) );
    
	if( lb_word == ub_word ) 
	  table[lb_word] |= (masked_lb & masked_ub);
	else {
	  table[lb_word] |= masked_lb;
	  table[ub_word] |= masked_ub;
	  while( --ub_word > lb_word )
	    table[ub_word] = full;
	}
      }
    }


    inline bool operator[](const int i)
    {
      return fastMember(i);
    }


    std::string toString()const
    {
      std::string s;
      if( !empty() ) {
	int last = NOVAL, cur=min(), aft;
	bool flag=false;
	do{
	  aft = next(cur);

	  if(aft != cur+1 || cur != last+1) {
	    if( flag )
	      s += " ";
	    s += int2string( cur );
	    flag = true;
	  } else if(flag) {
	    s += "..";
	    flag = false;
	  }
	  last = cur;
	  cur = aft;
	} while( cur != NOVAL && cur != last );
      }
      return s;
    }

    void  print(std::ostream & os) const 
    {
      if( table )
	print( os, "{,}" );
    }

    void  print(std::ostream & os, 
		const char *delim) const
    {
      os << delim[0];
      if( !empty() ) {
	int last = NOVAL, cur=min(), aft;
	bool flag=false;
	do{
	  aft = next(cur);

	  if(aft != cur+1 || cur != last+1) {
	    if( flag )
	      os << delim[1];
	    os << cur;
	    flag = true;
	  } else if(flag) {
	    os << "..";
	    flag = false;
	  }
	  last = cur;
	  cur = aft;
	} while( cur != NOVAL && cur != last );
      }
      os << delim[2];
    }

    void  printBits(std::ostream & os) const 
    {
      os << "[";
      for(int i=neg_words; i<pos_words; ++i)
	showUint( table[i], os );
      os << "]";
    }

  };




  typedef Bitset< unsigned int > BitSet;

  /**********************************************
   * List
   **********************************************/

  /// Sparse set representation
  class List 
  {
  public:

    /*!@name Parameters*/
    //@{
    unsigned int *absidx;
    unsigned int capacity;

    /// current size
    unsigned int size;
    /// list of values
    int *list_;
    /// values' indices
    unsigned int *index_;
    //@}

    /*!@name Constructors*/
    //@{
    List();
    virtual ~List();
    void init(const int, const int, const int=0);
    void init(Bitset<unsigned int>&, const int=0);
    //void pointTo(List& l);
    //@}    

    /*!@name Accessors*/
    //@{  
    inline bool member(const int elt)const 
    {
      return index_[elt]<size;
    } 
  
    inline bool empty()const 
    {
      return !size;
    } 

    inline int next(const int elt) const
    {
      unsigned int idx = index_[elt]+1;
      return (idx < size ? list_[idx] : elt);
    }

    inline int operator[](const int idx) const
    {
      return list_[idx];
    }

    inline int& operator[](const int idx)
    {
      return list_[idx];
    }
    //@}

    /*!@name List Manipulation*/
    //@{
    inline void fill()
    {
      size = capacity;
    }

    inline void clear()
    {
      size = 0;
    }
  
    inline void setTo(const int elt)
    {
      size=1;
      index_[*list_] = index_[elt];
      list_[index_[elt]] = *list_;
      *list_ = elt;
      index_[elt] = 0;
    }

    inline void erase(const int elt)
    {
      --size;
      index_[list_[size]] = index_[elt];
      list_[index_[elt]] = list_[size];
      list_[size] = elt;
      index_[elt] = size;
    }

    inline int next()
    {
      return list_[size];
    }

    inline int pop()
    {
      return list_[--size];
    }

    inline int popHead()
    {
      --size;
      index_[list_[size]] = 0;
      const int elt = *list_;
      *list_ = list_[size];
      list_[size] = elt;
      index_[elt] = size;
      return elt;
    }

    inline int head()
    {
      return *list_;
    }
    
    inline int back()
    {
      return list_[size-1];
    }

    inline void insert(const int elt)
    {

      //assert( index_[elt] >= size );

      //assert(check());

      index_[list_[size]] = index_[elt];
      list_[index_[elt]] = list_[size];
      list_[size] = elt;
      index_[elt] = size;
      ++size;

      //assert(check());
    }

    inline void revertTo(const int level)
    {
      size = level;
    }

    inline void index()
    {
      for(unsigned int i=0; i<capacity; ++i)
	index_[list_[i]] = i;
    }
    //@}

    /*!@name Miscellaneous*/
    //@{
    /// Print on out stream
    virtual void print(std::ostream& o) const ;
    //@}
  };

};
     
#endif
