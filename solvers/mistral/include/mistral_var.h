
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

/*! \file var.h
    \brief Header for alternative implementations of Variable.
*/

/**********************************************
 * Some alternative variable implementations 
 **********************************************/

#ifndef _VAR_H
#define _VAR_H

#include <mistral_csp.h>
#include <mistral_glo.h>


namespace Mistral {

  /********************************************
   * Reversible Boolean
   ********************************************/
  /*! \class ReversibleBool
    \brief Backtrackable Boolean
  */
  class ReversibleBool : public ReversibleObj
  {
  public:

    /*!@name Parameters*/
    //@{
    /// 0 -> {}, 1 -> {0}, 2 -> {1}, 3 -> {0,1}
    int state;
    //@}

    /*!@name Constructors*/
    //@{
    ReversibleBool() : ReversibleObj()
    {
      state = 3;
    }
    ReversibleBool(const int s) : ReversibleObj()
    {
      state = s;
    }
    //@}
  
    /*!@name Accessors*/
    //@{
    inline operator const int() const
    {
      return state;
    }
  
    inline void operator= ( const int nstat )
    {
      state = nstat;
      store->push( this );
    }
    //@}

    /*!@name Backtrack*/
    //@{
    inline void restore()
    {
      state = 3;
    }
    //@}

    /*!@name Miscellaneous*/
    //@{
    /// Print on out stream
    virtual void print(std::ostream& o) const {}
    //@}
  };


  /********************************************
   * Reversible Integer List + Bitset
   ********************************************/
  /*! \class ReversibleBitList
    \brief Backtrackable List+Set
  */
  class ReversibleBitList : public ReversibleObj
  {
  private:
    /*!@name Parameters*/
    //@{
    unsigned int *absidx;
  
  public:
    /// size trail
    Vector<unsigned int> size_;
    /// level trail
    Vector<int> lvl_;
    /// current size
    unsigned int size;
    /// list of values
    int *list_;
    /// values' indices
    unsigned int *index_;
    /// Bitset of values
    BitSet value_;
    //@}

    /*!@name Constructors*/
    //@{
    ReversibleBitList() : ReversibleObj() {}
    virtual ~ReversibleBitList();
    virtual void setValue(BitSet&, const int, const int, const unsigned int);
    //@}

    /*!@name Accessors*/
    //@{
    inline operator const BitSet() const
    {
      return value_;
    }

    inline int next(const int v) const
    {
      unsigned int idx = index_[v]+1;
      return (idx < size ? list_[idx] : NOVAL);
    }

    inline int min()const {return value_.min();}
    inline int max()const {return value_.max();}
    inline void unionTo    (BitSet* s) const {value_.unionTo(s);} 
    inline void copyTo     (BitSet* s) const {s->copy(value_);} 
    inline void intersectTo(BitSet* s) const {value_.intersectTo(s);} 
    inline void setminusTo (BitSet* s) const {value_.setminusTo(s);} 
    inline bool member(const int elt)const {return index_[elt]<size;} 
    inline bool empty()const {return !size;} 
    inline bool included (const BitSet& s)const {return value_.included(s);} 
    inline bool include (const BitSet& s)const {return s.included(value_);} 
    inline bool intersect(const BitSet& s)const {return value_.intersect(s);} 
    inline bool intersect(const int lb, const int ub)const {return value_.intersect(lb,ub);} 
    inline bool included (const BitSet* s)const {return value_.included(s);} 
    inline bool include (const BitSet* s)const {return s->included(value_);} 
    inline bool intersect(const BitSet* s)const {return value_.intersect(s);}
    //@}

    /*!@name Backtrack*/
    //@{
    inline void saveup()
    {
      if( *level != lvl_.back() ) {
	lvl_.push( *level );
	size_.push( size );
	store->push( this );
      }
    }

    inline void restore()
    {
      unsigned int i=size;
      size_.pop( size );
      lvl_.pop();
    
      while( i < size ) {
	value_.fastInvert( list_[i] );
	++i;
      }
    }
    //@}

    /*!@name List Manipulation*/
    //@{
    inline void setTo(int elt)
    {
      saveup();
      size=1;
      index_[*list_] = index_[elt];
      list_[index_[elt]] = *list_;
      *list_ = elt;
      index_[elt] = 0;
      value_.clear();
      value_.fastInvert( elt );
    }
  
    inline void erase(const int elt)
    {
      saveup();
      --size;
      index_[list_[size]] = index_[elt];
      list_[index_[elt]] = list_[size];
      list_[size] = elt;
      index_[elt] = size;
      value_.fastInvert(elt);
    }
    //@}
  
    /*!@name Miscellaneous*/
    //@{
    /// Print on out stream
    virtual void print(std::ostream& o) const ;
    //@}
  };


  /********************************************
   * Reversible Integer List 
   ********************************************/
  /*! \class ReversibleIntList
    \brief Backtrackable List
  */
  class ReversibleIntList : public ReversibleObj
  {
  private:
    /*!@name Parameters*/
    //@{
    unsigned int *absidx;
    unsigned int capacity;
  
    /// size trail
    Vector<unsigned int> size_;
    /// level trail
    Vector<int> lvl_;

  public:

    /// current size
    unsigned int size;
    /// list of values
    int *list_;
    /// values' indices
    unsigned int *index_;
    //@}

    /*!@name Constructors*/
    //@{
    ReversibleIntList() : ReversibleObj() {}
    virtual ~ReversibleIntList();
    virtual void setValue(const int, const int, const int);
    virtual void setValue(BitSet&, const int, const int, const int);
    //@}    

    /*!@name Backtrack*/
    //@{
    inline void saveup()
    {
      if( *level != lvl_.back() ) {
	lvl_.push( *level );
	size_.push( size );
	store->push( this );
      }
    }

    inline void restore()
    {
      size_.pop( size );
      lvl_.pop();
    }
    //@}
  
    /*!@name Accessors*/
    //@{  
    inline bool member(const int elt)const {
      return index_[elt]<size;
    } 
    inline bool empty()const {return !size;} 

    inline int next(const int elt) const
    {
      unsigned int idx = index_[elt]+1;
      return (idx < size ? list_[idx] : NOVAL);
    }

    inline int operator[](const int idx) const
    {
      return list_[idx];
    }
    //@}

    /*!@name List Manipulation*/
    //@{
    inline void fill()
    {
      saveup();
      size = capacity;
    }
  
    inline void setTo(const int elt)
    {
      saveup();
      size=1;
      index_[*list_] = index_[elt];
      list_[index_[elt]] = *list_;
      *list_ = elt;
      index_[elt] = 0;
    }

    inline void erase(const int elt)
    {
      saveup();
      --size;
      index_[list_[size]] = index_[elt];
      list_[index_[elt]] = list_[size];
      list_[size] = elt;
      index_[elt] = size;
    }

    inline void insert(const int elt)
    {
      saveup();
      index_[list_[size]] = index_[elt];
      list_[index_[elt]] = list_[size];
      list_[size] = elt;
      index_[elt] = size;
      ++size;
    }
    //@}

    /*!@name Miscellaneous*/
    //@{
    /// Print on out stream
    virtual void print(std::ostream& o) const ;

    virtual void printDebug(std::ostream& o) const 
    {
      unsigned int osize = size_[0];
      for(unsigned int i=0; i<osize; ++i) {
	if(i==size) std::cout << "|" ;
	else std::cout << " ";
	std::cout << std::setw(2) << list_[i];
      }
      std::cout << std::endl;
      for(int unsigned i=0; i<osize; ++i)
	std::cout << std::setw(3) << absidx[i];
      std::cout << std::endl;
    }
    //@}
  };


  // /********************************************
  //  * Reversible Integer MultiList
  //  ********************************************/
  // /*! \class ReversibleMultiList
  //   \brief Backtrackable MultiList
  // */
  // class ReversibleMultiList : public ReversibleObj
  // {
  //  private:
  //   /*!@name Parameters*/
  //   //@{
  //   int *next;
  //   int *prev;
  //   int *head;
  //   int *tail;

  //   /// size trail
  //   Vector<int> *delta_;
  //   /// level trail
  //   Vector<int> lvl_;
  //   Vector<int> size_;

  //  public:

  //   /*!@name Constructors*/
  //   //@{
  //   ReversibleMultiList() : ReversibleObj() {}
  //   virtual ~ReversibleMultiList();
  //   virtual void setValue(const int, const int, const int, const int*);
  //   //@}    

  //   /*!@name Backtrack*/
  //   //@{
  //   inline void saveup()
  //   {
  //     if( *level != lvl_.back() ) {
  //       lvl_.push( *level );
  //       size_.push( size );
  //       store->push( this );
  //     }
  //   }

  //   inline void restore()
  //   {
  //     size_.pop( size );
  //     lvl_.pop();
  //   }
  //   //@}
  
  //   /*!@name Accessors*/
  //   //@{  
  //   inline bool member(const int elt)const {return index_[elt]<size;} 
  //   inline bool empty()const {return !size;} 

  //   inline int next(const int elt) const
  //   {
  //     int idx = index_[elt]+1;
  //     return (idx < size ? list_[idx] : NOVAL);
  //   }

  //   inline int operator[](const int idx) const
  //   {
  //     return list_[idx];
  //   }
  //   //@}

  //   /*!@name List Manipulation*/
  //   //@{
  //   inline void fill()
  //   {
  //     saveup();
  //     size = capacity;
  //   }
  
  //   inline void setTo(const int elt)
  //   {
  //     saveup();
  //     size=1;
  //     index_[*list_] = index_[elt];
  //     list_[index_[elt]] = *list_;
  //     *list_ = elt;
  //     index_[elt] = 0;
  //   }

  //   inline void erase(const int elt)
  //   {
  //     saveup();
  //     --size;
  //     index_[list_[size]] = index_[elt];
  //     list_[index_[elt]] = list_[size];
  //     list_[size] = elt;
  //     index_[elt] = size;
  //   }

  //   inline void insert(const int elt)
  //   {
  //     saveup();
  //     index_[list_[size]] = index_[elt];
  //     list_[index_[elt]] = list_[size];
  //     list_[size] = elt;
  //     index_[elt] = size;
  //     ++size;
  //   }
  //   //@}

  //   /*!@name Miscellaneous*/
  //   //@{
  //   /// Print on out stream
  //   virtual void print(std::ostream& o) const ;

  //   virtual void printDebug(std::ostream& o) const 
  //   {
  //     int osize = size_[0];
  //     for(int i=0; i<osize; ++i) {
  //       if(i==size) std::cout << "|" ;
  //       else std::cout << " ";
  //       std::cout << std::setw(2) << list_[i];
  //     }
  //     std::cout << std::endl;
  //     for(int i=0; i<osize; ++i)
  //       std::cout << std::setw(3) << absidx[i];
  //     std::cout << std::endl;
  //   }
  //   //@}
  // };


  /********************************************
   * Reversible Set // TO RENAME IntSet?
   ********************************************/
  /*! \class ReversibleSet
    \brief Backtrackable Set
  */
  class ReversibleSet : public ReversibleObj
  {
    //private:
  public:
    /*!@name Parameters*/
    //@{
    /// bitset trail
    BitSet *val_;
    /// level trail
    Vector<int> *lvl_;
    /// current level
    int glvl;
    bool clone;

  public:
    /// current bitset
    BitSet value;
    //@}

    /*!@name Constructors*/
    //@{
    ReversibleSet() : ReversibleObj() 
    {
      lvl_ = NULL;
    }
  
    virtual void setValue(BitSet& e, const int l);
    virtual void setValue(const int lb, const int ub, const unsigned int val=BitSet::full);

    virtual ~ReversibleSet() 
    {
      if( lvl_ ) {
	lvl_ += value.neg_words;
	delete [] lvl_;

	val_ += value.neg_words;
	delete [] val_;
      }

      if( clone ) {
	value.neg_words = 0;
	value.table = NULL;
      }
    }
    //@}

    /*!@name Accessors*/
    //@{  
    inline operator const BitSet() const
    {
      return value;
    }

    inline uint size()const {return value.size();}
    inline int min()const {return value.min();}
    inline int max()const {return value.max();}
    inline void unionTo    (BitSet* s) const {value.unionTo(s);} 
    inline void copyTo     (BitSet* s) const {s->copy(value);} 
    inline void intersectTo(BitSet* s) const {value.intersectTo(s);} 
    inline void setminusTo (BitSet* s) const {value.setminusTo(s);} 
    inline bool member(const int elt)const {return value.member(elt);} 
    inline int next(const int elt) const {return value.next(elt);} 
    inline int prev(const int elt) const {return value.prev(elt);} 
    inline bool empty()const {return value.empty();} 
    inline bool included (const BitSet& s)const {return value.included(s);} 
    inline bool include (const BitSet& s)const {return s.included(value);} 
    inline bool intersect(const BitSet& s)const {return value.intersect(s);} 
    inline bool intersect(const int lb, const int ub)const {return value.intersect(lb,ub);} 
    inline bool included (const BitSet* s)const {return value.included(s);} 
    inline bool include (const BitSet* s)const {return s->included(value);} 
    inline bool intersect(const BitSet* s)const {return value.intersect(s);}
    //@}

    /*!@name Backtrack*/
    //@{
    inline void saveup()
    {
      int l=*level;
      if( glvl < l ) 
	{
	  glvl = l;
	  int i = value.pos_words;
	  int j = value.neg_words;
	  unsigned int w;
	  while( i-- > j ) 
	    {
	      w = value.table[i];
	      if( w != val_[i].table[lvl_[i].size-1] ) 
		{
		  val_[i].table[lvl_[i].size] = w;
		  lvl_[i].push(glvl);
		}
	    }
	  store->push( this );
	}
    }
  
    inline void restore()
    {
      int i = value.pos_words;
      int j = value.neg_words;
      int l = *level;
      glvl = 0;
      while( i-- > j ) 
	{
	  value.table[i] = val_[i].table[lvl_[i].size-1];
	  if( lvl_[i].back() == l )
	    lvl_[i].pop();
	  if( lvl_[i].back() > glvl )
	    glvl = lvl_[i].back();
	}
    }
    //@}

    /*!@name List Manipulation*/
    //@{  
    inline void setTo(int v)
    {
      saveup(); // saveup( v >> BitSet::EXP );
      value.clear();
      value.insert(v);
    }
    inline void unionWith(const BitSet& s)
    {
      saveup(); // saveup( v >> BitSet::EXP );
      value.unionWith(s);
    }
    inline void intersectWith(const BitSet& s)
    {
      saveup();
      value.intersectWith(s);
    }
    inline void setminusWith(const BitSet& s)
    {
      saveup();
      value.setminusWith(s);
    }
    inline void erase(const int elt)
    {
      saveup();
      value.fastInvert(elt);
    }
    inline void insert(const int elt)
    {
      saveup();
      value.insert(elt);
    }
    inline void setMin(int l)
    {
      saveup();
      value.setMin(l);
    }
    inline void setMax(int u)
    {
      saveup();
      value.setMax(u);
    }
    inline void removeInterval(int l, int u)
    {
      saveup();
      value.removeInterval(l,u);
    }
    //@}

    /*!@name Miscellaneous*/
    //@{
    /// Print on out stream
    virtual void print(std::ostream& o) const 
    {
      value.print( o );
    }
    //@}
  };




  /********************************************
   * Reversible Sparse Set
   ********************************************/
  /*! \class ReversibleSparseSet
    \brief Backtrackable Sparse Set 
  */
  class ReversibleSparseSet //: public ReversibleObj
  {
  
  private:
  
    /*!@name Parameters*/
    //@{
    Solver *solver;
    /// List of variables whose truth value are known
    List decided;
    /// level trail
    Vector<int> btlvl_;
    /// Truth values
    BitSet value;
    //@}

  public:

    /*!@name Parameters*/
    //@{
    /// number of variables
    int size;
    //@}

    /*!@name Constructors*/
    //@{
    ReversibleSparseSet();
    void initialise( Solver* );
    virtual ~ReversibleSparseSet();
    //@}

    /*!@name Accessors*/
    //@{
    inline bool setValue(const int idx, const int var, const int val);
    inline int getValue(const int var)
    {
      if( decided.member(var) ) 
	return value.member(var);
      return 2;
    }
    inline bool contain(const int var, const int val)
    {
      if( val >> 1 ) return false;
      if( decided.member(var) ) 
	return (value.member(var) == val);
      return true;
    }
    //@}

    /*!@name Backtrack*/
    //@{
    inline void save()
    {
      btlvl_.push( decided.size );
    }
    inline void restore()
    {
      decided.revertTo( btlvl_.pop() );
    }
    //@}
  
    /*!@name Miscellaneous*/
    //@{
    /// Print on out stream
    virtual void print(std::ostream& o) const ;
    //@}
  };



  /********************************************
   * Reversible Primitive Type
   ********************************************/
  /*! \class ReversibleNum
    \brief Backtrackable Primitive
  */
  template < class T >
  class ReversibleNum : public ReversibleObj
  {
  
  public:
    /*!@name Parameters*/
    //@{  
    /// value trail
    Vector<T> val_;
    /// level trail
    Vector<int> lvl_;
    /// current value
    T value;
    //@}

    /*!@name Constructors*/
    //@{  
    ReversibleNum() : ReversibleObj() {};
    virtual ~ReversibleNum() {};
    inline void setValue( const T x, const int t=16 )
    {
      if( !(val_.size || lvl_.size) ) {
	value = x;
	val_.push(value);
	lvl_.push(0);
      }
    }
    //@}

    /*!@name Accessors*/
    //@{  
    inline operator const T() const
    {
      return value;
    }
    //   inline const T out() const
    //   {
    //     return val_[val_.size+1];
    //   }
    //@}  

    /*!@name Backtrack*/
    //@{
    inline void saveup()
    {
      if( *level != lvl_.back() ) {
	lvl_.push( *level );
	val_.push( value );
	store->push( this );
      }
    }

    inline void restore()
    {
      val_.pop( value );
      lvl_.pop();
    }
    //@}

    /*!@name Manipulation*/
    //@{  
    inline void operator= ( const T x )
    {
      saveup();
      value = x;
    }
  
    inline void operator+= ( const T x )
    {
      saveup();
      value += x;
    }
  
    inline void operator-= ( const T x )
    {
      saveup();
      value -= x;
    }

    inline void operator*= ( const T x )
    {
      saveup();
      value *= x;
    }
  
    inline void operator/= ( const T x )
    {
      saveup();
      value /= x;
    }

    inline void operator|= ( const T x )
    {
      saveup();
      value |= x;
    }
  
    inline void operator&= ( const T x )
    {
      saveup();
      value &= x;
    }
  
    inline void operator^= ( const T x )
    {
      saveup();
      value ^= x;
    }
  
    inline void operator++ ()
    {
      saveup();
      ++value;
    }
  
    inline void operator-- ()
    {
      saveup();
      --value;
    } 
    //@}
  };


  /********************************************
   * Reversible Word
   ********************************************/
  /*! \class ReversibleWord
    \brief Backtrackable Word, used for bitsets
  */
  class ReversibleWord : public ReversibleObj
  {
  
  public:
    /*!@name Parameters*/
    //@{  
    /// current value
    unsigned int *word_;
    /// value trail
    unsigned int *val_;
    /// level trail
    int *lvl_;
    int last_;

    //@}

    /*!@name Constructors*/
    //@{  
    ReversibleWord() : ReversibleObj() {} 
    void init( unsigned int* w, const unsigned int count ) 
    { 
      word_ = w;
      val_ = new unsigned int[count+1];
      lvl_ = new int[count+1];

      last_ = 0;
      val_[last_] = *word_;
      lvl_[last_] = 0;
    };
    virtual ~ReversibleWord() {
      delete [] val_;
      delete [] lvl_;
    };
    //@}

    /*!@name Backtrack*/
    //@{  
    inline void saveup()
    {
      if( *level != *lvl_ ) {
	++last_;
	lvl_[last_] = *level;
	val_[last_] = *word_;
	store->push( this );
      }
    }

    inline void restore()
    {
      *word_ = val_[last_];
      --last_;
    }
    //@}

    /*!@name Manipulation*/
    //@{  
    inline void operator= ( const unsigned int x )
    {
      saveup();
      *word_ = x;
    }
  
    inline void operator+= ( const unsigned int x )
    {  
      saveup();
      *word_ += x;
    }
  
    inline void operator-= ( const unsigned int x )
    {
      saveup();
      *word_ -= x;
    }

    inline void operator*= ( const unsigned int x )
    {
      saveup();
      *word_ *= x;
    }
  
    inline void operator/= ( const unsigned int x )
    {
      saveup();
      *word_ /= x;
    }

    inline void operator|= ( const unsigned int x )
    {
      saveup();
      *word_ |= x;
    }
  
    inline void operator&= ( const unsigned int x )
    {
      saveup();
      *word_ &= x;
    }
  
    inline void operator^= ( const unsigned int x )
    {
      saveup();
      *word_ ^= x;
    }
  
    inline void operator++ ()
    {
      saveup();
      ++*word_;
    }
  
    inline void operator-- ()
    {
      saveup();
      --*word_;
    } 

    inline operator const unsigned int() const { return *word_; }
    //@}
  };


  /********************************************
   * Finite Domain Integer Variable
   ********************************************/
  /*! \class VariableDomain
    \brief Finite domain integer variable
  */
  class VariableDomain : public VariableInt {

  protected:

    inline bool boundchanged(const int v) 
    {
      if( v == vmin ) {
	vmin = values.next( v );
	return true;
      }
      if( v == vmax ) {
	vmax = values.prev( v );
	return true;
      }
      return false;
    }

    inline bool boundchanged(const BitSet& s, const bool spin) 
    {
      bool bc = false;
      if( spin ^ s.member(vmin) ) {
	vmin = values.min();
	bc = true;
      }
      if( spin ^ s.member(vmax) ) {
	vmax = values.max();
	bc = true;
      }

      return bc;
    }

    inline void triggerEvent( const bool bc )
    {
      if( vmin == vmax ) 
	//ACQueue->push(this, Constraint::VALUETRIGGER);
	ACQueue->push(id, Constraint::VALUETRIGGER);
      else if( bc )
	//ACQueue->push(this, Constraint::RANGETRIGGER);
	ACQueue->push(id, Constraint::RANGETRIGGER);
      else //ACQueue->push(this, Constraint::DOMAINTRIGGER);
	ACQueue->push(id, Constraint::DOMAINTRIGGER);
    }


  public:

    /*!@name Parameters*/
    //@{
    /// Minimum value in the domain
    ReversibleNum<int> vmin;
    /// Maximum value in the domain
    ReversibleNum<int> vmax;
    /// Minimum value
    int negval;
    /// Maximum value + 1
    int posval;
    /// The current domain
    BitSet values;
    //@}


    /*!@name Constructors*/
    //@{
    VariableDomain( Solver *s ) ;
    virtual ~VariableDomain() {}
    virtual int& getIntDomain();

    void setVariable(const int low, const int up);
    void setVariable(BitSet *s, const int length, 
		     const int low, const int up);
    //@}

    /*!@name Accessors and Iterators*/
    //@{
    /// Returns the assigned value if it exists
    inline int value() const 
    {
      return ( vmin );    
    }

    /// Returns the minimum value in the domain
    inline int min() const { return vmin; }

    /// Returns the maximum value in the domain
    inline int max() const { return vmax; }

    /// Returns the minimum value that could belong to the domain
    inline int minCapacity() const {return negval;}

    /// Returns 1 + the maximum value that could belong to the domain
    inline int maxCapacity() const {return posval;}
    //@}

    /*!@name Query methods*/
    //@{
    /// Whether or not the Variable is bound to a ground value
    inline bool isGround() const { return (vmin == vmax); }
    /// Whether or not the Variable is bound to a given ground value
    inline bool equal(const int v) const { return (vmin == vmax && vmin == v); }
    /*!
      Whether "v" is currently contained in the domain
    */
    inline bool contain(const int v) const 
    { 
      return ( values.member(v) ); 
    }
    inline bool fastContain(const int v) const 
    { 
      return ( values.fastMember(v) ); 
    }

    /*!
      Whether the domain has a nonempty intersection
      with the interval [l..u]
    */
    inline bool intersect(const int lo, const int up) const
    {
      return ( vmin <= up && vmax >= lo );
    }
    /*!
      Whether the domain is included in
      the interval [l..u]
    */
    inline bool included(const int lo, const int up) const
    {
      return ( vmin >= lo && vmax <= up );
    }
    /*!
      Whether the domain has a nonempty intersection
      with the Variable x
    */
    inline bool intersect(const VariableInt* x) const 
    {
      return x->intersect( values );
    }
    /*!
      Whether the domain is included
      in the Variable x 
    */
    inline bool included(const VariableInt* x) const 
    {
      return x->included( values );
    }
    /*!
      Whether the domain is included
      in the Variable x 
    */
    inline bool include(const VariableInt* x) const 
    {
      return x->included( this );
    }
    /*!
      Whether the domain has a nonempty intersection
      with the set s 
    */
    inline bool intersect(const BitSet* s) const
    {
      return values.intersect(s);
    }
    /*!
      Whether the domain is included
      in the set s
    */
    inline bool included(const BitSet* s) const
    {
      return values.included(s);
    }
    inline bool include(const BitSet* s) const
    {
      return s->included(values);
    }
    /*!
      Whether the domain has a nonempty intersection
      with the set s 
    */
    inline bool intersect(const BitSet& s) const
    {
      return values.intersect(s);
      //return values.fastIntersect(s);
    }
    /*!
      Whether the domain has a nonempty intersection
      with the set s
    */
    inline bool intersect(const BitSet& s, int& r) const
    {
      return values.fastIntersect(s,r);
      //return values.wordIntersect(s);
    }
    inline bool wordIntersect(const BitSet& s) const
    {

      //    s.print( std::cout );
      //    std::cout << " i ";
      //    values.print( std::cout );
    
      return values.wordIntersect(s);
    }
    /*!
      Whether the domain is included
      in the set s
    */
    inline bool included(const BitSet& s) const
    {
      return values.included(s);
    }
    inline bool include(const BitSet& s) const
    {
      return s.included(values);
    }
    /*!
      Intersect its domain with a set s
    */
    inline void intersectTo( BitSet* s ) const
    {
      s->intersectWith(values);
    }
    /*!
      Do the union of its domain with a set s
    */
    inline void unionTo( BitSet* s ) const
    {
      s->unionWith(values);
    }
    /*!
      Do the union of its domain with a set s
    */
    inline void copyTo( BitSet* s ) const
    {
      s->copy(values);
    }
    /*!
      Intersect its domain with a set s
    */
    inline void intersectTo( BitSet& s ) const
    {
      s.intersectWith(values);
    }
    /*!
      Do the union of its domain with a set s
    */
    inline void unionTo( BitSet& s ) const
    {
      s.unionWith(values);
    }
    /*!
      Do the union of its domain with a set s
    */
    inline void copyTo( BitSet& s ) const
    {
      s.copy(values);
    }
    //@}

  };




  // class BitsetIterator : public DomainIterator {

  //  public: 

  //   BitsetIterator( BitSet& vals ) 
  //     {
  //       domain.pointTo( vals );
  //     }

  // //   void BitsetIterator::init( int vmin, BitSet& vals ) 
  // //     {
  // //       domain.pointTo( vals );
      
  // //     }

  //   virtual ~BitsetIterator() 
  //     {
  //       domain.table = NULL;
  //       domain.neg_words = 0;
  //     }

  //   BitSet domain;
  //   int curval;

  //   inline bool next() { curval=domain.next(curval); return curval != NOVAL; }
  //   inline operator const int() const { return curval; }
  // };


  class BitsetIterator : public DomainIterator {

  public: 

    BitsetIterator( const int vmin, const BitSet& vals ) 
    {
      init( vmin, vals );
    }

    BitsetIterator( const BitSet& vals ) 
    {
      init( vals.min(), vals );
      //       //domain.pointTo( vals );
      //       table = vals.table;
      //       pos_words = vals.pos_words;
      //       neg_words = vals.neg_words;
      //       curval = vals.min();
      //       wordIdx = (curval >> BitSet::EXP);
      //       remainingVals = (table[wordIdx] & ((~0) << ((curval & BitSet::CACHE)+1)));
    }

    void init( const int vmin, const BitSet& vals ) 
    {
      table = vals.table;
      pos_words = vals.pos_words;
      neg_words = vals.neg_words;
      curval = vmin;
      wordIdx = (curval >> BitSet::EXP);
      remainingVals = (table[wordIdx] & ((~0) << ((curval & BitSet::CACHE)+1)));
      //offset = 32*wordIdx;
    }

    virtual ~BitsetIterator() 
    {
      //domain.table = NULL;
      //domain.neg_words = 0;
    }


    //       SAT    14025 NDS     2830 BTS/s      612714 PPGS/s       2990045 PPGS    5.01 s
    //       SAT      251 NDS     2188 BTS/s      527344 PPGS/s         47461 PPGS    0.23 s
    //     UNSAT    34988 NDS     3387 BTS/s      611582 PPGS/s       6317651 PPGS   10.47 s
    //     UNSAT     4197 NDS     2835 BTS/s      412213 PPGS/s        610076 PPGS    1.62 s
    //       SAT    15276 NDS     3505 BTS/s      655034 PPGS/s       2842851 PPGS    4.48 s

    //    1/ 0.6    13747 NDS     3239 BTS/s       2561616 PPGS   4.224 s


    //BitSet domain;
    unsigned int *table;
    int pos_words;
    int neg_words;
    int wordIdx;
    int remainingVals;
    int curval;
    //int offset;




    inline bool next() { 
      //return next1(); 
      //return next2(); 
      return next3(); 
    }
    inline bool next1() {

      /* While all bits in the word were processed */
      while (remainingVals == 0) {
	if (++wordIdx == pos_words) /* If the last word was processed */
	  return false;			 /* This is the end of the iteration */
	remainingVals = table[wordIdx]; /* Move to the next word */
	//offset += 32;
      }

      int b = remainingVals & -remainingVals; /* b contains a mask with only the least significant bit of remaingVals set */
      remainingVals ^= b;			    /* Unset the least significant bit */

      /* Set curval to the index of the least significant bit */
      curval = ((b & 0xFFFF0000) != 0) << 4;
      curval |= ((b & 0xFF00FF00) != 0) << 3;
      curval |= ((b & 0xF0F0F0F0) != 0) << 2;
      curval |= ((b & 0xCCCCCCCC) != 0) << 1;
      curval |= ((b & 0xAAAAAAAA) != 0);

      curval += (wordIdx - neg_words)*32; /* Add to curval its offset */
      return true;
    }


    inline bool next2() {
      const unsigned int T1 = 0x34572610;    // Lookup table 1
      const unsigned int T2 = 0xBCDFAE98;    // Lookup table 2 (T2 = T1 + 0x88888888)
      const unsigned int seq = 0x1D;         // A Bruijn sequence of length 8

      /* While all bits in the word were processed */
      while (remainingVals == 0) {
	if (++wordIdx == pos_words) /* If the last word was processed */
	  return false;			 /* This is the end of the iteration */
	remainingVals = table[wordIdx]; /* Move to the next word */
      }

      curval = (wordIdx - neg_words)*32; /* Add to curval its offset */
    
      unsigned int b = remainingVals & -remainingVals; // Create a mask with the least significant bit set to one
      remainingVals ^= b;
    
      const unsigned int y = b * seq;        // Use the multiplication by b as a shift operator
      if ((b & 0xFFFF0000) == 0) {		 // First step of the binary search
	if ((b & 0xFF00) == 0) {		 // Second step of the binary search
	  curval += (T1 >> ((y >> 3) & 0xFC)) & 0x0F;  // Look up the first table
	  return true;
	}
	else {
	  curval += (T2 >> ((y >> 11) & 0xFC)) & 0x0F; // Look up the second table
	  return true;
	}
      }
      else {
	if ((b & 0xFF000000) == 0) {	         // Second step of the binary search
	  curval += 16 | ((T1 >> ((y >> 19) & 0xFC)) & 0x0F); // Look up the first table
	  return true;
	}
	else {
	  curval += 16 | ((T2 >> ((y >> 27) & 0xFC)) & 0x0F); // Look up the second table
	  return true;
	}
      }
    }


    inline bool next3() {

      //std::cout << "bnext: " << curval << " " << remainingVals << std::endl;

      /* While all bits in the word were processed */
      while (remainingVals == 0) {
	if (++wordIdx == pos_words) /* If the last word was processed */
	  return false;			 /* This is the end of the iteration */
	remainingVals = table[wordIdx]; /* Move to the next word */
      }

      union {float f; unsigned int i; } t;
      unsigned int b = remainingVals & -remainingVals;
      remainingVals ^= b;
      t.f = (float)b; // cast the least significant bit in v to a float
      b = t.i >> 23;

      //curval = b + (wordIdx - neg_words) * 32 - 0x7F;
      curval = b + wordIdx * 32 - 0x7F;

      return true;
    }

    inline bool operator++() { 

      /* While all bits in the word were processed */
      while (remainingVals == 0) {
	if (++wordIdx == pos_words) /* If the last word was processed */
	  return false;			 /* This is the end of the iteration */
	remainingVals = table[wordIdx]; /* Move to the next word */
      }

      union {float f; unsigned int i; } t;
      unsigned int b = remainingVals & -remainingVals;
      remainingVals ^= b;
      t.f = (float)b; // cast the least significant bit in v to a float
      b = t.i >> 23;

      curval = b + wordIdx * 32 - 0x7F;
      return true;
    }
    inline operator const int() const { return curval; }
  };


  /*********************************************************
   * VariableBit || An implementation of Variable using a bit vector
   *********************************************************/
  /*! \class VariableBit
    \brief Implementation of integer domain Variable

    This implementation uses a revesible vector of bits for 
    representing the domain.
  */
  class VariableBit : public VariableDomain {

  public:

    /*!@name Parameters*/
    //@{
    /// Current domain size
    ReversibleNum<int> size;  
    /// Domain
    ReversibleSet domain;
    //@}

    /*!@name Constructors*/
    //@{
    /*! 
      Constructor, create a domain containing all values
      in the interval [0,-,s-1]
    */
    VariableBit(Solver *s, const int v);

    /*! 
      Constructor, create a domain containing all values
      in the interval [l,-,u]
    */
    VariableBit(Solver *s, const int l, const int u);

    /*! 
      Constructor, create a domain containing all values 
      in the set e
    */
    VariableBit(Solver *s, BitSet *e, 
		const int lg, 
		const int lo, 
		const int up);

    VariableBit(Solver *s) : VariableDomain(s) {}
    virtual ~VariableBit();
    //@}

    /*!@name Accessors and Iterators*/
    //@{
    /// access to the domain iterator
    DomainIterator *begin();

    /// Returns a random value in the domain
    int random() const;
//     {
//       int val;
//       if((1 + vmax - vmin ) == size) val = (randint(size) + vmin);
//       else {
// 	val = randint(size);      
// 	static_bit_domain_it->init( vmin, values );

// 	//DomainIterator *iter = begin();
// 	while(val--) static_bit_domain_it->next();
// 	  //iter->next();   
// 	val = (*static_bit_domain_it);
//       }
//       return val;
      
//       //if((1 + vmax - vmin ) != size) return values.random();
//       //return (randint(size) + vmin);
//     }

    /// Return the domain size
    inline int domsize() const 
    {
      return size;
    }
    /// Whether or not the Variable is currently an interval
    inline bool isRange() const 
    { 
      return (size == (vmax - vmin + 1)); 
    }
    /// Return the first (equal to min) value in the domain
    inline int first() const 
    { 
      return vmin; 
    }
    /// Return the last (equal to max) value in the domain
    inline int  last() const 
    { 
      return vmax; 
    }
    /// Return the smallest value currently in the domain that is strictly greater than "v"
    inline int getNext(const int v) const 
    { 
      return ( values.next(v) ); 
    }
    /*! Used when looping on all values in the domain, return    
      true if "v" is valid value 
    */
    inline bool good(const int v) const 
    { 
      return ( v != NOVAL ); 
    }
    /*! Set "v" to the smallest value currently in the domain that is strictly greater than "v"
     */
    inline bool setNext(int& v) const 
    { 
      return ( (v = values.next(v)) != NOVAL ); 
    }
  
    /// Returns the minimum absolute value in [0..infty] \\inter domain, NOVAL if there are none
    inline int minPosAbs() const 
    {
      if( vmax < 0 ) return NOVAL;
      if( vmin > 0 ) return vmin;
      if( !values.table || values.member( 0 ) ) return 0;
      return values.next( 0 );
    }
  
    /// Returns the minimum absolute value in [-infty..0] \\inter domain, NOVAL if there are none
    inline int minNegAbs() const  
    {
      if( vmin > 0 ) return NOVAL;
      if( vmax < 0 ) return vmax;
      if( !values.table || values.member( 0 ) ) return 0;
      return values.prev( 0 );
    }
    //@}
  

    /*!@name Domain handling*/
    //@{
    /*! 
      Remove the value "v"
    */
    inline bool remove(const int v)
    {
      if(!contain(v)) return true;
      if( 1 == size ) return false;

      domain.erase( v );
      --size;
  
      triggerEvent( boundchanged( v ) );

      return true;
    }

    /*! 
      Remove all values but "v"
    */
    inline bool setDomain(const int v)
    {
      if( !contain( v ) ) return false;
   
      if( 1 < size ) {

	size = 1;

	domain.setTo( v );

	vmin = v;
	vmax = v;

	//ACQueue->push(this, Constraint::VALUETRIGGER)
	ACQueue->push(id, Constraint::VALUETRIGGER);
      }

      return true;
    }

    /*! 
      Remove all values strictly lower than low
    */
    inline bool setMin(const int low)
    {
      if(vmax < low)  return false;  
      if(vmin >= low) return true;  

      domain.setMin(low);
      size = values.size();

      if( values.fastMember(low) ) vmin = low;
      else vmin = values.min();

      triggerEvent( true );

      return true;
    }

    /*! 
      Remove all values strictly greater than up
    */
    inline bool setMax(const int up)
    {
      if(vmin > up)  return false;  
      if(vmax <= up) return true;   

      domain.setMax(up);
      size = values.size();

      if( values.fastMember(up) ) vmax = up;
      else vmax = values.max();

      triggerEvent( true );

      return true;
    }

    /*! 
      Remove all values that do not appear in the array "d" of length "n"
    */
    inline bool setDomain(const int* d, const int n)
    {
      BitSet aux(vmin, vmax, BitSet::full);
      int i=n;
      while( i-- ) aux.insert( d[i] );
      return setDomain( aux );
    }

    /*! 
      Remove all values that do not appear in the set "s"
    */
    inline bool setDomain(const BitSet& s)
    {
      if( !values.intersect(s) ) return false;
      if( values.included(s) ) return true;

      domain.intersectWith(s);
      size = values.size();

      triggerEvent( boundchanged( s, true ) );

      return true;
    }

    /*! 
      Remove all values that do not appear in the current domain of the Variable "x"
    */
    inline bool setDomain( VariableInt* x ) const
    {
      if( 1 == size ) 
	return x->setDomain( vmin );
      return x->setDomain( values );
    }

    /*! 
      Remove all values that belong to the set "s"
    */
    inline bool removeSet(const BitSet& s)
    {  
      if( values.included(s) ) return false;
      if( !values.intersect(s) ) return true;

      domain.setminusWith(s);    
      size = values.size();

      triggerEvent( boundchanged( s, false ) );

      return true;
    }

    /*! 
      Remove all values in the interval [lo..up]
    */
    inline bool removeRange(const int lo, const int up)
    {
      if( !(values.intersect( lo, up )) )
	return true;

      if( lo <= vmin ) 
	return setMin(up+1);
      if( up >= vmax ) 
	return setMax(lo-1);  

      domain.removeInterval(lo, up);
      size = values.size();

      //ACQueue->push(this, Constraint::DOMAINTRIGGER);  
      ACQueue->push(id, Constraint::DOMAINTRIGGER);  

      return true;
    }
    //@}

    /*!@name Propagation method*/
    //@{
    inline bool revise( BitSet* supports, 
			int* residues, 
			VariableInt* X)
    {
      //int vali = vmin;
      bool consistent;
      DomainIterator *valit = begin();
      if( X->isWord ) 
	do consistent = (X->wordIntersect( supports[*valit] ) || remove( *valit ));
	while( consistent && valit->next() );
      else 
	do consistent = (X->intersect( supports[*valit], residues[*valit] ) || remove( *valit ));
	while( consistent && valit->next() );
      return consistent;
    }

    inline bool revise( Constraint* c, const int i )
    {

      int vali = vmin;
      do {
	// This procedure does two things:
	// * first it checks if we can say that the value is supported
	//   without constraint check
	// * second it ready the solution in order to be check
	//   by findSupport	

	if( !( c->firstSupport(i, vali) || 
	       c->findSupport(i, vali) ||
	       remove(vali) ) ) 
	  return false;	      
	
      } while( setNext(vali) );
  
      return true;
    }
    //@}
//     /*!@name Propagation method*/
//     //@{
//     inline bool revise( BitSet* supports, 
// 			int* residues, 
// 			VariableInt* X)
//     {
//       //int vali = vmin;
//       VariableInt* not_consistent; 
// //	  bool consistent;	// DG Change
//       DomainIterator *valit = begin();
//       if( X->isWord ) 
// 	do not_consistent = (!(X->wordIntersect( supports[*valit] )) && remove( *valit ));
// 	while( not_consistent && valit->next() );
//       else 
// 	do not_consistent = (!(X->intersect( supports[*valit], residues[*valit] )) && remove( *valit ));
// 	while( !not_consistent && valit->next() );
//       return not_consistent;
//     }

//     inline bool revise( Constraint* c, const int i )
//     {

//       int vali = vmin;
//       do {
// 	// This procedure does two things:
// 	// * first it checks if we can say that the value is supported
// 	//   without constraint check
// 	// * second it ready the solution in order to be check
// 	//   by findSupport	

// 	if( !( c->firstSupport(i, vali) || 
// 	       c->findSupport(i, vali) ||
// 	       remove(vali) ) ) 
// 	  return false;	      
	
//       } while( setNext(vali) );
  
//       return true;
//     }
//     //@}

    /*!@name Miscellaneous*/
    //@{
    void check_domain(int x) const ;
    /// Print on out stream
    virtual void print(std::ostream& o) const;
    virtual void printshort(std::ostream&) const ;
    virtual void printDomain(std::ostream&) const ;
    virtual int getType() const { return VariableInt::BIT; };
    //@}
  };



  /*********************************************************
   * VariableBitset || An implementation of Variable using a bit vector
   *********************************************************/
  /*! \class VariableBitset
    \brief Implementation of integer domain Variable

    This implementation uses a revesible vector of bits for 
    representing the domain.
  */
  class VariableBitset : public VariableDomain {

  public:

    void checkDomain( const char* message ) const {
      if( domain.size() != values.size() ) {
	std::cerr << message << " dom/val" << std::endl; 
	exit( 0 );
      }
      if( values.size() != (unsigned int)size ) {
	std::cerr << message << " val/size" << std::endl; 
	exit( 0 );
      }
      if( vmin != values.min() ) {
	std::cerr << message << " min/val" << std::endl; 
	exit( 0 );
      }
      if( vmax != values.max() ) {
	std::cerr << message << " max/val" << std::endl; 
	exit( 0 );
      }
    }

    /*!@name Parameters*/
    //@{
    /// Current domain size
    ReversibleNum<int> size;  
    /// Domain
    Bitset< ReversibleWord > domain;
    //@}

    /*!@name Constructors*/
    //@{
    /*! 
      Constructor, create a domain containing all values
      in the interval [0,-,s-1]
    */
    VariableBitset(Solver *s, const int v);

    /*! 
      Constructor, create a domain containing all values
      in the interval [l,-,u]
    */
    VariableBitset(Solver *s, const int l, const int u);

    /*! 
      Constructor, create a domain containing all values 
      in the set e
    */
    VariableBitset(Solver *s, BitSet *e, 
		   const int lg, 
		   const int lo, 
		   const int up);

    VariableBitset(Solver *s) : VariableDomain(s) {}
    virtual ~VariableBitset();
    //@}

    /*!@name Accessors and Iterators*/
    //@{
    /// access to the domain iterator
    DomainIterator *begin();

    /// Returns a random value
    int random() const;
//     {
//       int val;
//       if((1 + vmax - vmin ) == size) val = (randint(size) + vmin);
//       else {
// 	val = randint(size);      
// 	DomainIterator *iter = begin();
// 	while(val--) iter->next();   
// 	val = (*iter);
//       }

//       return val;
//       //return values.random();
//     }
    /// Return the domain size
    inline int domsize() const 
    {
      return size;
    }
    /// Whether or not the Variable is currently an interval
    inline bool isRange() const 
    { 
      return (size == (vmax - vmin + 1)); 
    }
    /// Return the first (equal to min) value in the domain
    inline int first() const 
    { 
      return vmin; 
    }
    /// Return the last (equal to max) value in the domain
    inline int  last() const 
    { 
      return vmax; 
    }
    /// Return the smallest value currently in the domain that is strictly greater than "v"
    inline int getNext(const int v) const 
    { 
      return ( values.next(v) ); 
    }
    /*! Used when looping on all values in the domain, return    
      true if "v" is valid value 
    */
    inline bool good(const int v) const 
    { 
      return ( v != NOVAL ); 
    }
    /*! Set "v" to the smallest value currently in the domain that is strictly greater than "v"
     */
    inline bool setNext(int& v) const 
    { 
      return ( (v = values.next(v)) != NOVAL ); 
    }
  
    /// Returns the minimum absolute value in [0..infty] \\inter domain, NOVAL if there are none
    inline int minPosAbs() const 
    {
      if( vmax < 0 ) return NOVAL;
      if( vmin > 0 ) return vmin;
      if( !values.table || values.member( 0 ) ) return 0;
      return values.next( 0 );
    }
  
    /// Returns the minimum absolute value in [-infty..0] \\inter domain, NOVAL if there are none
    inline int minNegAbs() const  
    {
      if( vmin > 0 ) return NOVAL;
      if( vmax < 0 ) return vmax;
      if( !values.table || values.member( 0 ) ) return 0;
      return values.prev( 0 );
    }
    //@}
  

    /*!@name Domain handling*/
    //@{
    /*! 
      Remove the value "v"
    */
    inline bool remove(const int v)
    {
      if(!contain(v)) return true;
      if( 1 == size ) return false;
    
      domain.erase( v );
      --size;
  
      triggerEvent( boundchanged( v ) );

      return true;
    }

    /*! 
      Remove all values but "v"
    */
    inline bool setDomain(const int v)
    {
      if( !contain( v ) ) return false;
      if( 1 < size ) {
	size = 1;

	domain.setTo( v );

	vmin = v;
	vmax = v;

	//ACQueue->push(this, Constraint::VALUETRIGGER);
	ACQueue->push(id, Constraint::VALUETRIGGER);
      }
      return true;
    }

    /*! 
      Remove all values strictly lower than low
    */
    inline bool setMin(const int low)
    {
      if(vmax < low)  return false;  
      if(vmin >= low) return true;  

      domain.setMin(low);
      size = values.size();

      if( values.fastMember(low) ) vmin = low;
      else vmin = values.min();

      triggerEvent( true );

      return true;
    }

    /*! 
      Remove all values strictly greater than up
    */
    inline bool setMax(const int up)
    {
      if(vmin > up)  return false;  
      if(vmax <= up) return true;   

      domain.setMax(up);
      size = values.size();

      if( values.fastMember(up) ) vmax = up;
      else vmax = values.max();

      triggerEvent( true );

      return true;
    }

    /*! 
      Remove all values that do not appear in the array "d" of length "n"
    */
    inline bool setDomain(const int* d, const int n)
    {
      BitSet aux(vmin, vmax, BitSet::full);
      int i=n;
      while( i-- ) aux.insert( d[i] );
      return setDomain( aux );
    }

    /*! 
      Remove all values that do not appear in the set "s"
    */
    inline bool setDomain(const BitSet& s)
    {
      if( !values.intersect(s) ) return false;
      if( values.included(s) ) return true;

      domain.intersectWith(s); // TO REPAIR
      size = values.size();

      triggerEvent( boundchanged( s, true ) );

      return true;
    }

    /*! 
      Remove all values that do not appear in the current domain of the Variable "x"
    */
    inline bool setDomain( VariableInt* x ) const
    {
      if( 1 == size ) 
	return x->setDomain( vmin );
      return x->setDomain( values );
    }

    /*! 
      Remove all values that belong to the set "s"
    */
    inline bool removeSet(const BitSet& s)
    {  
      if( values.included(s) ) return false;
      if( !values.intersect(s) ) return true;

      domain.setminusWith(s);    // TO REPAIR
      size = values.size();

      triggerEvent( boundchanged( s, false ) );

      return true;
    }

    /*! 
      Remove all values in the interval [lo..up]
    */
    inline bool removeRange(const int lo, const int up)
    {
      if( !(values.intersect( lo, up )) )
	return true;

      if( lo <= vmin ) 
	return setMin(up+1);
      if( up >= vmax ) 
	return setMax(lo-1);  

      domain.removeInterval(lo, up);
      size = values.size();

      //ACQueue->push(this, Constraint::DOMAINTRIGGER);  
      ACQueue->push(id, Constraint::DOMAINTRIGGER);  

      return true;
    }
    //@}

    /*!@name Propagation method*/
    //@{
    inline bool revise( BitSet* supports, 
			int* residues, 
			VariableInt* X)
    {
      //int vali = vmin;
      bool consistent;
      DomainIterator *valit = begin();
      if( X->isWord )
	do consistent = (X->wordIntersect( supports[*valit] ) || remove( *valit ));
	while( consistent && valit->next() );
      else 
	do consistent = (X->intersect( supports[*valit], residues[*valit] ) || remove( *valit ));
	while( consistent && valit->next() );
      return consistent;
    }

    inline bool revise( Constraint* c, const int i )
    {
 
      int vali = vmin;
      do {
	// This procedure does two things:
	// * first it checks if we can say that the value is supported
	//   without constraint check
	// * second it ready the solution in order to be check
	//   by findSupport	
 
	if( !( c->firstSupport(i, vali) || 
	       c->findSupport(i, vali) ||
	       remove(vali) ) )
	  return false;	
	
      } while( setNext(vali) );
  
      return true;
    }
    //@}
//     /*!@name Propagation method*/
//     //@{
//     inline bool revise( BitSet* supports, 
// 			int* residues, 
// 			VariableInt* X)
//     {
//       //int vali = vmin;
//       VariableInt* not_consistent; // DG Change
//       DomainIterator *valit = begin();
//       if( X->isWord )
// 	do not_consistent = (!(X->wordIntersect( supports[*valit] )) && remove( *valit ));
// 	while( !not_consistent && valit->next() );
//       else 
// 	do not_consistent = (!(X->intersect( supports[*valit], residues[*valit] )) && remove( *valit ));
// 	while( not_consistent && valit->next() );
//       return not_consistent;
//     }

//     inline bool revise( Constraint* c, const int i )
//     {
 
//       int vali = vmin;
//       do {
// 	// This procedure does two things:
// 	// * first it checks if we can say that the value is supported
// 	//   without constraint check
// 	// * second it ready the solution in order to be check
// 	//   by findSupport	
 
// 	if( !( c->firstSupport(i, vali) || 
// 	       c->findSupport(i, vali) ||
// 	       remove(vali) ) )
// 	  return false;	
	
//       } while( setNext(vali) );
  
//       return true;
//     }
//     //@}

    /*!@name Miscellaneous*/
    //@{
    void check_domain(int x) const ;
    /// Print on out stream
    virtual void print(std::ostream& o) const;
    virtual void printshort(std::ostream&) const ;
    virtual void printDomain(std::ostream&) const ;
    virtual int getType() const { return VariableInt::BIT; };
    //@}
  };


  class ListIterator : public DomainIterator {

  public: 

    virtual ~ListIterator() {}

    int* vals;
    int  size;

    inline bool next() { return size--; }
    inline operator const int() const { return vals[size]; }
  };


  /*********************************************************
   * VariableBitList || An implementation of Variable using a list
   *********************************************************/
  /*! \class VariableBitList
    \brief Implementation of integer domain Variable

    This implementation uses a revesible list
    representing the domain.
  */
  class VariableBitList : public VariableDomain {

  public:

    /*!@name Parameters*/
    //@{
    /// Domain
    ReversibleBitList domain;
    //@}

    /*!@name Constructors*/
    //@{
    /*! 
      Constructor, create a domain containing all values
      in the interval [0,-,s-1]
    */
    VariableBitList(Solver *s, const int v); 

    /*! 
      Constructor, create a domain containing all values
      in the interval [l,-,u]
    */
    VariableBitList(Solver *s, const int l, const int u); 

    /*! 
      Constructor, create a domain containing all values 
      in the set e
    */
    VariableBitList(Solver *s, 
		    BitSet *e, 
		    const int lg, 
		    const int lo, 
		    const int up); 
  

    VariableBitList(Solver *s) : VariableDomain(s) {}
    virtual ~VariableBitList();
    //@}

    /*!@name Accessors and Iterators*/
    //@{
    /// access to the domain iterator
    DomainIterator *begin();

    /// Returns a random value
    inline int random() const
    {
      //return domain.list_[rand()%domain.size];
      return domain.list_[randint(domain.size)];
    }
    /// Return the domain size
    inline int domsize() const 
    {
      return domain.size;
    }
    /// Whether or not the Variable is currently an interval
    inline bool isRange() const 
    { 
      return (domain.size == (unsigned int)(vmax - vmin + 1)); 
    }
    /// Return the first (equal to min) value in the domain
    inline int first() const 
    { 
      return *(domain.list_);
    }
    /// Return the last (equal to max) value in the domain
    inline int  last() const 
    { 
      return domain.list_[domain.size-1];
    }
    /// Return the smallest value currently in the domain that is strictly greater than "v"
    inline int getNext(const int v) const 
    { 
      return ( domain.next(v) ); 
    }
    /*! Used when looping on all values in the domain, return    
      true if "v" is valid value 
    */
    inline bool good(const int v) const 
    { 
      return ( v != NOVAL ); 
    }
    /*! Set "v" to the smallest value currently in the domain that is strictly greater than "v"
     */
    inline bool setNext(int& v) const 
    { 
      return ( (v = domain.next(v)) != NOVAL ); 
    }
  
    /// Returns the minimum absolute value in [0..infty] \\inter domain, NOVAL if there are none
    inline int minPosAbs() const 
    {
      if( vmax < 0 ) return NOVAL;
      if( vmin > 0 ) return vmin;
      if( !values.table || values.member( 0 ) ) return 0;
      return values.next( 0 );
    }
  
    /// Returns the minimum absolute value in [-infty..0] \\inter domain, NOVAL if there are none
    inline int minNegAbs() const  
    {
      if( vmin > 0 ) return NOVAL;
      if( vmax < 0 ) return vmax;
      if( !values.table || values.member( 0 ) ) return 0;
      return values.prev( 0 );
    }
    //@}
  

    /*!@name Domain handling*/
    //@{

    /*! 
      Remove the value "v"
    */
    inline bool remove(const int v)
    {
      if(!contain(v)) return true;
      if( 1 == domain.size ) return false;

      domain.erase( v );
  
      triggerEvent( boundchanged( v ) );

      return true;
    }

    /*! 
      Remove all values but "v"
    */
    inline bool setDomain(const int v)
    {
      if( !contain( v ) ) return false;
   
      if( 1 < domain.size ) {

	domain.setTo( v );

	vmin = v;
	vmax = v;

	//ACQueue->push(this, Constraint::VALUETRIGGER);
	ACQueue->push(id, Constraint::VALUETRIGGER);
      }

      return true;
    }

    /*! 
      Remove all values strictly lower than l
    */
    inline bool setMin(const int l)
    {  
      return true;
    }

    /*! 
      Remove all values strictly greater than u
    */
    inline bool setMax(const int u)
    {  
      return true;
    }

    /*! 
      Remove all values that do not appear in the array "a" of length "l"
    */
    inline bool setDomain(const int* a, const int l)
    {  
      return true;
    }

    /*! 
      Remove all values that do not appear in the set "s"
    */
    inline bool setDomain(const BitSet& s)
    {  
      return true;
    }

    /*! 
      Remove all values that do not appear in the current domain of the Variable "x"
    */
    inline bool setDomain( VariableInt* x ) const
    {  
      if( 1 == domain.size ) 
	return x->setDomain( vmin );
      return x->setDomain( values );
    }

    /*! 
      Remove all values that belong to the set "s"
    */
    inline bool removeSet(const BitSet& s)
    {  
      return true;
    }

    /*! 
      Remove all values in the interval [l..u]
    */
    inline bool removeRange(const int l, const int u)
    {  
      return true;
    }
    //@}

    /*!@name Propagation method*/
    //@{
    inline bool revise( BitSet* supports, 
			int* residues, 
			VariableInt* X)
    {
      bool not_consistent = true;
      int *val = domain.list_, i=domain.size;
      if( X->isWord )
	while( not_consistent && i-- ) 
	  not_consistent = (X->wordIntersect( supports[val[i]] ) || remove( val[i] ));
      else 
	while( not_consistent && i-- )
	  not_consistent = (X->intersect( supports[val[i]], residues[val[i]] ) || remove( val[i] ));
      return not_consistent;
    }
    inline bool revise( Constraint* c, const int i )
    {

      bool not_consistent = true;
      int *val = domain.list_, j=domain.size;
      while( not_consistent && j-- ) {
	not_consistent = ( c->firstSupport(i, val[j]) || 
		       c->findSupport (i, val[j]) ||
		       remove(val[j]) );
      }
      return not_consistent;
    }
    //@}

    /*!@name Miscellaneous*/
    //@{
    void check_domain(int x) const ;
    /// Print on out stream
    virtual void print(std::ostream& o) const;
    virtual void printshort(std::ostream&) const ;
    virtual void printDomain(std::ostream&) const ;
    virtual int getType() const { return VariableInt::LIST; };
    //@}
  };


  /*********************************************************
   * VariableList || An implementation of Variable using a list
   *********************************************************/
  /*! \class VariableList
    \brief Implementation of integer domain Variable

    This implementation uses a revesible list
    representing the domain.
  */
  class VariableList : public VariableInt {

  public:

    /*!@name Parameters*/
    //@{
    // WARNING! min and max are not maintained
    /// Domain
    ReversibleIntList domain;
    //@}

    /*!@name Constructors*/
    //@{
    /*! 
      Constructor, create a domain containing all values
      in the interval [0,-,s-1]
    */
    VariableList(Solver *s, const int v); 

    /*! 
      Constructor, create a domain containing all values
      in the interval [l,-,u]
    */
    VariableList(Solver *s, const int l, const int u); 

    /*! 
      Constructor, create a domain containing all values 
      in the set e
    */
    VariableList(Solver *s, 
		 BitSet* e, 
		 const int lg, 
		 const int lo, 
		 const int up); 
  

    VariableList(Solver *s) : VariableInt(s) {}
    virtual ~VariableList();
    //@}

    /*!@name Accessors and Iterators*/
    //@{
    /// Returns the assigned value if it exists
    inline int value() const 
    {
      return *(domain.list_);
    }

    /// Returns the minimum value in the domain
    inline int min() const { return *(domain.list_); } //TODO

    /// Returns the maximum value in the domain
    inline int max() const { return domain.list_[domain.size-1]; } //TODO

    /// Returns the minimum value that could belong to the domain
    inline int minCapacity() const {return *(domain.list_);} //TODO

    /// Returns 1 + the maximum value that could belong to the domain
    inline int maxCapacity() const {return domain.list_[domain.size-1];} //TODO
    //@}

    /*!@name Query methods*/
    //@{
    /// Whether or not the Variable is bound to a ground value
    inline bool isGround() const { return (domain.size == 1); }

    /// Whether or not the Variable is bound to a given ground value
    inline bool equal(const int v) const { return (domain.size == 1 && *(domain.list_) == v); }

    /*!
      Whether "v" is currently contained in the domain
    */
    inline bool contain(const int v) const 
    { 
      return ( domain.member(v) ); 
    }
    inline bool fastContain(const int v) const 
    { 
      return ( domain.member(v) ); 
    }

    /*!
      Whether the domain has a nonempty intersection
      with the interval [l..u]
    */
    inline bool intersect(const int lo, const int up) const
    {
      return false; //TODO
    }
    /*!
      Whether the domain is included in
      the interval [l..u]
    */
    inline bool included(const int lo, const int up) const
    {
      return false; //TODO
    }
    /*!
      Whether the domain has a nonempty intersection
      with the Variable x
    */
    inline bool intersect(const VariableInt* x) const 
    {
      return false; //TODO
    }
    /*!
      Whether the domain is included
      in the Variable x 
    */
    inline bool included(const VariableInt* x) const 
    {
      return false; //TODO
    }
    /*!
      Whether the domain is included
      in the Variable x 
    */
    inline bool include(const VariableInt* x) const 
    {
      return false; //TODO
    }
    /*!
      Whether the domain has a nonempty intersection
      with the set s 
    */
    inline bool intersect(const BitSet* s) const
    {
      return false; //TODO
    }
    /*!
      Whether the domain is included
      in the set s
    */
    inline bool included(const BitSet* s) const
    {
      return false; //TODO
    }
    inline bool include(const BitSet* s) const
    {
      return false; //TODO
    }
    /*!
      Whether the domain has a nonempty intersection
      with the set s 
    */
    inline bool intersect(const BitSet& s) const
    {
      return false; //TODO
      //return values.fastIntersect(s);
    }
    /*!
      Whether the domain has a nonempty intersection
      with the set s
    */
    inline bool intersect(const BitSet& s, int& r) const
    {
      return false; //TODO
      //return values.wordIntersect(s);
    }
    inline bool wordIntersect(const BitSet& s) const
    {
      return false; //TODO
    }
    /*!
      Whether the domain is included
      in the set s
    */
    inline bool included(const BitSet& s) const
    {
      return false; //TODO
    }
    inline bool include(const BitSet& s) const
    {
      return false; //TODO
    }
    /*!
      Intersect its domain with a set s
    */
    inline void intersectTo( BitSet* s ) const
    {
      //TODO
    }
    /*!
      Do the union of its domain with a set s
    */
    inline void unionTo( BitSet* s ) const
    {
      //TODO
    }
    /*!
      Do the union of its domain with a set s
    */
    inline void copyTo( BitSet* s ) const
    {
      //TODO
    }
    /*!
      Intersect its domain with a set s
    */
    inline void intersectTo( BitSet& s ) const
    {
      //TODO
    }
    /*!
      Do the union of its domain with a set s
    */
    inline void unionTo( BitSet& s ) const
    {
      //TODO
    }
    /*!
      Do the union of its domain with a set s
    */
    inline void copyTo( BitSet& s ) const
    {
      //TODO
    }


    /// access to the domain iterator
    DomainIterator *begin();

    /// Returns a random value
    inline int random() const
    {
      //return domain.list_[rand()%domain.size];
      return domain.list_[randint(domain.size)];
    }
    /// Return the domain size
    inline int domsize() const 
    {
      return domain.size;
    }
    /// Whether or not the Variable is currently an interval
    inline bool isRange() const 
    { 
      return false; 
    }
    /// Return the first (equal to min) value in the domain
    inline int first() const 
    { 
      return *(domain.list_);
    }
    /// Return the last (equal to max) value in the domain
    inline int  last() const 
    { 
      return domain.list_[domain.size-1];
    }
    /// Return the smallest value currently in the domain that is strictly greater than "v"
    inline int getNext(const int v) const 
    { 
      return ( domain.next(v) ); 
    }
    /*! Used when looping on all values in the domain, return    
      true if "v" is valid value 
    */
    inline bool good(const int v) const 
    { 
      return ( v != NOVAL ); 
    }
    /*! Set "v" to the smallest value currently in the domain that is strictly greater than "v"
     */
    inline bool setNext(int& v) const 
    { 
      return ( (v = domain.next(v)) != NOVAL ); 
    }
  
    /// Returns the minimum absolute value in [0..infty] \\inter domain, NOVAL if there are none
    inline int minPosAbs() const 
    {
      return 0;//TODO
    }
  
    /// Returns the minimum absolute value in [-infty..0] \\inter domain, NOVAL if there are none
    inline int minNegAbs() const  
    {
      return 0;//TODO
    }
    //@}
  

    /*!@name Domain handling*/
    //@{

    /*! 
      Remove the value "v"
    */
    inline bool remove(const int v)
    {
      if(!contain(v)) return true;
      if( 1 == domain.size ) return false;

      domain.erase( v );

      if( domain.size == 1 )
	ACQueue->push(id, Constraint::VALUETRIGGER);  
      else
	ACQueue->push(id, Constraint::DOMAINTRIGGER);  

      return true;
    }

    /*! 
      Remove all values but "v"
    */
    inline bool setDomain(const int v)
    {
      if( !contain( v ) ) return false;
   
      if( 1 < domain.size ) {

	domain.setTo( v );

	ACQueue->push(id, Constraint::VALUETRIGGER);
      }

      return true;
    }

    /*! 
      Remove all values strictly lower than l
    */
    inline bool setMin(const int l)
    {  
      return true; //TODO
    }

    /*! 
      Remove all values strictly greater than u
    */
    inline bool setMax(const int u)
    {  
      return true; //TODO
    }

    /*! 
      Remove all values that do not appear in the array "a" of length "l"
    */
    inline bool setDomain(const int* a, const int l)
    {  
      return true; //TODO
    }

    /*! 
      Remove all values that do not appear in the set "s"
    */
    inline bool setDomain(const BitSet& s)
    {  
      return true; //TODO
    }

    /*! 
      Remove all values that do not appear in the current domain of the Variable "x"
    */
    inline bool setDomain( VariableInt* x ) const
    {  
      return true; //TODO
    }

    /*! 
      Remove all values that belong to the set "s"
    */
    inline bool removeSet(const BitSet& s)
    {  
      return true; //TODO
    }

    /*! 
      Remove all values in the interval [l..u]
    */
    inline bool removeRange(const int l, const int u)
    {  
      return true; //TODO
    }
    //@}

    /*!@name Propagation method*/
    //@{
    inline bool revise( BitSet* supports, 
			int* residues, 
			VariableInt* X)
    {
      return false; //TODO
    }
    inline bool revise( Constraint* c, const int i )
    {
      bool not_consistent = true;
      int *val = domain.list_, j=domain.size;
      while( not_consistent && j-- ) {
	not_consistent = ( c->firstSupport(i, val[j]) || 
		       c->findSupport (i, val[j]) ||
		       remove(val[j]) );
      }
      return not_consistent;
    }
    //@}

    /*!@name Miscellaneous*/
    //@{
    void check_domain(int x) const ;
    /// Print on out stream
    virtual void print(std::ostream& o) const;
    virtual void printshort(std::ostream&) const ;
    virtual void printDomain(std::ostream&) const ;
    virtual int getType() const { return VariableInt::INTLIST; };
    //@}
  };



  class BoolIterator : public DomainIterator {

  public: 

    virtual ~BoolIterator() {}

    int curval;

    inline bool next() { 
      if(curval < 3) return false;  
      curval >>= 1;
      return true;
    }
    inline operator const int() const { return (curval >> 1); }
  };



  /*********************************************************
   * VariableBool || Boolean variable
   *********************************************************/
  /*! \class VariableBool
    \brief Implementation of Boolean Variable
  */
  class VariableBool : public VariableInt {

  protected:
    inline bool setState( const int nstat ) ;

  public:

    /*!@name Parameters*/
    //@{
    ReversibleBool domain;
    //@}

    /*!@name Constructors*/
    //@{
    VariableBool(Solver *s); 
    virtual int& getIntDomain();
    virtual ~VariableBool();
    //@}

    /*!@name Accessors and Iterators*/
    //@{
    /// access to the domain iterator
    DomainIterator *begin();

    /// Return the first (equal to min) value in the domain
    inline int first() const
    {
      return !(domain.state & 1);
    }
    /// Return the last (equal to max) value in the domain
    inline int  last() const
    {
      return (domain.state >> 1);
    }
    /// Return the smallest value currently in the domain that is strictly greater than "v"
    inline int getNext(const int v) const
    {
      return ((domain.state >> 1) & (1-v));
    }
    /*! 
      Used when looping on all values in the domain, return
      true if "v" is valid value
    */
    inline bool good(const int v) const
    {
      return ( v );
    }
    /*! 
      Set "v" to the smallest value currently in the domain that is strictly greater than "v"
    */
    inline bool setNext(int& v) const
    {
      return ( v = ((domain.state >> 1) & (1-v)) );
    }

    /// Returns the minimum absolute value in [0..infty] \\inter domain, NOVAL if there are none
    inline int minPosAbs() const
    {
      if( domain.state & 1 ) return 0;
      return 1;
    }

    /// Returns the minimum absolute value in [-infty..0] \\inter domain, NOVAL if there are none
    inline int minNegAbs() const
    {
      if( domain.state & 1 ) return 0;
      return NOVAL;
    }

    /*!
      Whether "v" is currently contained in the domain
    */
    inline bool contain(const int v) const
    {
      return ( !(v >> 1) && (domain.state & (v+1)) );
    }
    inline bool fastContain(const int v) const
    {
      return ( domain.state & (v+1) );
    }

    /// Returns the assigned value if it exists
    inline int value() const
    {
      return ( domain.state == 3 ? NOVAL : domain.state-1 );
    }

    /// Returns the minimum value in the domain
    inline int min() const
    {
      return !(domain.state & 1);
    }

    /// Returns the maximum value in the domain
    inline int max() const
    {
      return ( domain.state >> 1 );
    }

    /// Returns the minimum value that could belong to the domain
    inline int minCapacity() const
    {
      return 0;
    }

    /// Returns 1 + the maximum value that could belong to the domain
    inline int maxCapacity() const
    {
      return 2;
    }

    /// Returns a random value
    inline int random() const
    {
      if( domain.state == 3 ) //return (rand()%2);
	return (randint(2));
      else return domain.state-1;
    }
    /// Returns the domain size
    inline int domsize() const
    {
      return ( (domain.state & 1) + ((domain.state & 2) >> 1) );
    }
    //@}

    /*!@name Query methods*/
    //@{
    /// Whether or not the Variable is currently an interval
    inline bool isRange() const
    {
      return true;
    }
    /// Whether or not the Variable is bound to a ground value
    inline bool isGround() const
    {
      return ( domain.state != 3 );
    }
    /// Whether or not the Variable is bound to a given ground value
    inline bool equal(const int v) const { return (domain.state-1 == v); }
    /*!
      Whether the domain has a nonempty intersection
      with the interval [l..u]
    */
    inline bool intersect(int lo, int up) const
    {
      if( lo > up || lo > 1 || up < 0 ) return false;
      int nstat = 3;
      if(lo < 1) nstat ^= 1;
      if(up > 0) nstat ^= 2;
      return ( nstat & domain.state );
    }
    inline bool included(int lo, int up) const
    {
      //return ( !(domain.state ^ lo) && !(domain.state ^ up) );
      return ( !(domain.state ^ lo) && !(domain.state ^ up) );
    }
    /*!
      Whether the domain has a nonempty intersection
      with the Variable x
    */
    inline bool intersect(const VariableInt* v) const
    {
      return (v->contain( !(domain.state & 1) ) || v->contain( (domain.state >> 1) )) ;
    }
    /*!
      Whether the domain is included
      in the Variable x
    */
    inline bool included(const VariableInt* v) const
    {
      return (v->min() >= !(domain.state & 1) && v->max() <= (domain.state >> 1));
    }
    inline bool include(const VariableInt* v) const
    {
      return v->contain( min() ) && v->contain( max() );
    }
    /*!
      Whether the domain has a nonempty intersection
      with the set s
    */
    inline bool intersect(const BitSet* s) const
    {
      return ( s->pos_words > 0 && s->neg_words < 1 &&
	       (domain.state & s->table[0]) );
    }
    /*!
      Whether the domain is included
      in the set s
    */
    inline bool included(const BitSet* s) const
    {
      return (
	      s->pos_words > 0 &&
	      s->neg_words <= 0 &&
	      (s->table[0] & domain.state) == (unsigned int)(domain.state)
	      );
    }
    inline bool include(const BitSet* s) const
    {
      return (
	      s->size() <= 2 &&//change that
	      s->pos_words > 0 &&
	      s->neg_words <= 0 &&
	      (s->table[0] & domain.state) == s->table[0]
	      );
    }
    /*!
      Whether the domain has a nonempty intersection
      with the set s
    */
    inline bool intersect(const BitSet& s) const
    {
      return ( s.pos_words > 0 && s.neg_words < 1 &&
	       (domain.state & s.table[0]) );
    }
    /*!
      Whether the domain has a nonempty intersection
      with the set s
    */
    inline bool intersect(const BitSet& s, int& r) const
    {
      return ( s.pos_words > 0 && s.neg_words < 1 &&
	       (domain.state & s.table[0]) );
    }
    inline bool wordIntersect(const BitSet& s) const
    {
      return ( (domain.state & s.table[0]) );
    }
    /*!
      Whether the domain is included
      in the set s
    */
    inline bool included(const BitSet& s) const
    {
      return (
	      s.pos_words > 0 &&
	      s.neg_words <= 0 &&
	      (s.table[0] & domain.state) == (unsigned int)(domain.state)
	      );
    }
    inline bool include(const BitSet& s) const
    {
      return (
	      s.size() <= 2 &&//change that
	      s.pos_words > 0 &&
	      s.neg_words <= 0 &&
	      (s.table[0] & domain.state) == s.table[0]
	      );
    }
    /*!
      Intersect its domain with a set s
    */
    inline void intersectTo( BitSet* s ) const
    {
      int i=s->pos_words;
      int j=s->neg_words;
      while( i-- > j )
	if( !i ) s->table[i] &= domain.state;
	else s->table[i] = 0;
    }
    /*!
      Do the union of its domain with a set s
    */
    inline void unionTo( BitSet* s ) const
    {
      if( s->pos_words > 0 && s->neg_words < 1 )
	s->table[0] |= domain.state;
    }
    /*!
      Do the union of its domain with a set s
    */
    inline void copyTo( BitSet* s ) const
    {
      s->clear();
      if( s->pos_words > 0 && s->neg_words < 1 )
	s->table[0] |= domain.state;
    }
    /*!
      Intersect its domain with a set s
    */
    inline void intersectTo( BitSet& s ) const
    {
      int i=s.pos_words;
      int j=s.neg_words;
      while( i-- > j )
	if( !i ) s.table[i] &= domain.state;
	else s.table[i] = 0;
    }
    /*!
      Do the union of its domain with a set s
    */
    inline void unionTo( BitSet& s ) const
    {
      if( s.pos_words > 0 && s.neg_words < 1 )
	s.table[0] |= domain.state;
    }
    /*!
      Do the union of its domain with a set s
    */
    inline void copyTo( BitSet& s ) const
    {
      s.clear();
      if( s.pos_words > 0 && s.neg_words < 1 )
	s.table[0] |= domain.state;
    }
    //@}
    

    /*!@name Domain handling*/
    //@{
    /*!
      Remove the value "v"
    */
    bool remove(const int v);
    /*!
      Remove all values but "v"
    */
    bool setDomain(const int v);
    /*!
      Remove all values strictly lower than l
    */
    bool setMin(const int l);
    /*!
      Remove all values strictly greater than u
    */
    bool setMax(const int u);
    /*!
      Remove all values that do not appear in the array "a" of length "l"
    */
    bool setDomain(const int* a, const int l);
    /*!
      Remove all values that do not appear in the {@link BitSet set} "s"
    */
    bool setDomain(const BitSet& s);
    /*!
      Remove all values that do not appear in the current domain of the variable "x"
    */
    bool setDomain(VariableInt*) const;
    /*!
      Remove all values that belong to the {@link BitSet set} "s"
    */
    bool removeSet(const BitSet& s);
    /*!
      Remove all values in the interval [l..u]
    */
    bool removeRange(const int l, const int u);
    //@}

    /// The mandatory methods
    /*!@name Decision methods*/
    //@{
    bool revise( BitSet*, int*, VariableInt* );
    bool revise( Constraint*, const int );
    //@}

    /*!@name Miscellaneous*/
    //@{
    /// Print on out stream
    virtual void print(std::ostream& o) const;
    virtual void printshort(std::ostream&) const ;
    virtual void printDomain(std::ostream&) const ;
    virtual int getType() const { return VariableInt::BOOL; };
    //@}
  };



  // /*********************************************************
  //  * VariableBool2 || Boolean variable
  //  *********************************************************/
  // /*! \class VariableBool2
  //     \brief Implementation of Boolean Variable
  // */
  // class VariableBool2 : public VariableInt {

  // protected:
  
  //   /// Old utility function
  //   //inline bool setState( const int nstat ) ;
  //   //inline bool setVal( const int nstat ) ;

  //  public:

  //   /*!@name Parameters*/
  //   //@{
  //   /// Old Domain
  //   //ReversibleBool domain;

  //   /// New Domain
  //   ReversibleSparseSet& domainPool;
  //   int index;
  //   //@}

  //   /*!@name Constructors*/
  //   //@{
  //   VariableBool2(Solver *s); 
  //   virtual ~VariableBool2();
  //   //@}

  //   /*!@name Accessors and Iterators*/
  //   //@{
  //   /// access to the domain iterator
  //   DomainIterator *begin();

  //   /// Return the first (equal to min) value in the domain
  //   inline int first() const
  //   {

  //     //assert( !(domain.state & 1) == (domainPool.getValue(index) % 2) );

  //     return (domainPool.getValue(index) % 2);
  //   }
  //   /// Return the last (equal to max) value in the domain
  //   inline int  last() const
  //   {

  //     //assert( (domain.state >> 1) == (domainPool.getValue(index) > 0) );

  //     return (domainPool.getValue(index) > 0);
  //   }
  //   /// Return the smallest value currently in the domain that is strictly greater than "v"
  //   inline int getNext(const int v) const
  //   {    

  //     //assert( ((domain.state >> 1) & (1-v)) == (!v && (domainPool.getValue(index) > 0)) );

  //     return (!v && (domainPool.getValue(index) > 0));
  //   }
  //   /*! 
  //     Used when looping on all values in the domain, return
  //     true if "v" is valid value
  //   */
  //   inline bool good(const int v) const
  //     {
  //       return (v);
  //     }
  //   /*! 
  //     Set "v" to the smallest value currently in the domain that is strictly greater than "v"
  //   */
  //   inline bool setNext(int& v) const
  //     {

  //       //assert( ((domain.state >> 1) & (1-v)) == (!v && (domainPool.getValue(index) > 0)) );

  //       return ( v = (!v && (domainPool.getValue(index) > 0)) );
  //     }

  //   /// Returns the minimum absolute value in [0..infty] \\inter domain, NOVAL if there are none
  //   inline int minPosAbs() const
  //   {
  //     return (domainPool.getValue(index) % 2);
  //   }

  //   /// Returns the minimum absolute value in [-infty..0] \\inter domain, NOVAL if there are none
  //   inline int minNegAbs() const
  //   {
  //     if(domainPool.getValue(index) % 2)
  //       return 0;
  //     else return NOVAL;
  //   }

  //   /*!
  //      Whether "v" is currently contained in the domain
  //   */
  //   inline bool contain(const int v) const
  //   {
  //     return ( domainPool.contain(index,v) );
  //   }
  //   inline bool fastContain(const int v) const
  //   {
  //     return ( domainPool.contain(index,v) );
  //   }

  //   /// Returns the assigned value if it exists
  //   inline int value() const
  //   {
    
  //     int val = domainPool.getValue(index);
  //     if(val == 2) val = NOVAL;

  //     //assert( (domain.state == 3 ? NOVAL : domain.state-1) == val );

  //     return ( val );
  //   }

  //   /// Returns the minimum value in the domain
  //   inline int min() const
  //   {

  //     //assert( !(domain.state & 1) == (domainPool.getValue(index) % 2) );

  //     return (domainPool.getValue(index) % 2);
  //   }

  //   /// Returns the maximum value in the domain
  //   inline int max() const
  //   {

  //     //std::cout << (domain.state) << " " << (domainPool.getValue(index)) << std::endl;

  //     //assert( (domain.state >> 1) == (domainPool.getValue(index) > 0) );

  //     return (domainPool.getValue(index) > 0);
  //   }

  //   /// Returns the minimum value that could belong to the domain
  //   inline int minCapacity() const
  //   {
  //     return 0;
  //   }

  //   /// Returns 1 + the maximum value that could belong to the domain
  //   inline int maxCapacity() const
  //   {
  //     return 2;
  //   }

  //   /// Returns a random value
  //   inline int random() const
  //   {
  //     int val = domainPool.getValue(index);
  //     if( val == 2 ) return (rand()%2);
  //     return val;
  //   }
  //   /// Returns the domain size
  //   inline int domsize() const
  //   {

  //     //assert( ((domain.state & 1) + ((domain.state & 2) >> 1)) == (1 + ( domainPool.getValue(index) > 1 )) );

  //     return 1 + ( domainPool.getValue(index) > 1 );
  //   }
  //   //@}

  //   /*!@name Query methods*/
  //   //@{
  //   /// Whether or not the Variable is currently an interval
  //   inline bool isRange() const
  //   {
  //     return true;
  //   }
  //   /// Whether or not the Variable is bound to a ground value
  //   inline bool isGround() const
  //   {
  //     return ( domainPool.getValue(index) != 2 );
  //   }
  //   /// Whether or not the Variable is bound to a given ground value
  //   inline bool equal(const int v) const 
  //   { 
  //     return (domainPool.getValue(index) == v); 
  //   }
  //   /*!
  //      Whether the domain has a nonempty intersection
  //      with the interval [l..u]
  //   */
  //   inline bool intersect(int lo, int up) const
  //   {
  //     if( lo > up || lo > 1 || up < 0 ) return false;
  //     int nstat = 3;
  //     if(lo < 1) nstat ^= 1;
  //     if(up > 0) nstat ^= 2;
  //     return ( nstat & (domainPool.getValue(index)+1) );
  //   }
  //   inline bool included(int lo, int up) const
  //   {
  //     //return ( !(domain.state ^ lo) && !(domain.state ^ up) );
  //     int val = domainPool.getValue(index)+1;
  //     return ( !(val ^ lo) && !(val ^ up) );
  //   }
  //   /*!
  //      Whether the domain has a nonempty intersection
  //      with the Variable x
  //   */
  //   inline bool intersect(const VariableInt* v) const
  //   {
  //     int val = domainPool.getValue(index)+1;
  //     return (v->contain( !(val & 1) ) || v->contain( (val >> 1) )) ;
  //   }
  //   /*!
  //      Whether the domain is included
  //      in the Variable x
  //   */
  //   inline bool included(const VariableInt* v) const
  //   {
  //     int val = domainPool.getValue(index)+1;
  //     return (v->min() >= !(val & 1) && v->max() <= (val >> 1));
  //   }
  //   inline bool include(const VariableInt* v) const
  //   {
  //     return v->contain( min() ) && v->contain( max() );
  //   }
  //   /*!
  //      Whether the domain has a nonempty intersection
  //      with the set s
  //   */
  //   inline bool intersect(const BitSet* s) const
  //   {
  //     return ( s->pos_words > 0 && s->neg_words < 1 &&
  // 	     ((domainPool.getValue(index)+1) & s->table[0]) );
  //   }
  //   /*!
  //      Whether the domain is included
  //      in the set s
  //   */
  //   inline bool included(const BitSet* s) const
  //   {
  //     int val = domainPool.getValue(index)+1;
  //     return (
  // 	    s->pos_words > 0 &&
  // 	    s->neg_words <= 0 &&
  // 	    (s->table[0] & val) == val
  // 	    );
  //   }
  //   inline bool include(const BitSet* s) const
  //   {
  //     return (
  // 	    s->size() <= 2 &&//change that
  // 	    s->pos_words > 0 &&
  // 	    s->neg_words <= 0 &&
  // 	    (s->table[0] & (domainPool.getValue(index)+1)) == s->table[0]
  // 	    );
  //   }
  //   /*!
  //      Whether the domain has a nonempty intersection
  //      with the set s
  //   */
  //   inline bool intersect(const BitSet& s) const
  //   {
  //     return ( s.pos_words > 0 && s.neg_words < 1 &&
  // 	     ((domainPool.getValue(index)+1) & s.table[0]) );
  //   }
  //   /*!
  //      Whether the domain has a nonempty intersection
  //      with the set s
  //   */
  //   inline bool intersect(const BitSet& s, int& r) const
  //   {
  //     return ( s.pos_words > 0 && s.neg_words < 1 &&
  // 	     ((domainPool.getValue(index)+1) & s.table[0]) );
  //   }
  //   inline bool wordIntersect(const BitSet& s) const
  //   {
  //     return ( ((domainPool.getValue(index)+1) & s.table[0]) );
  //   }
  //   /*!
  //      Whether the domain is included
  //      in the set s
  //   */
  //   inline bool included(const BitSet& s) const
  //   {
  //     int val = domainPool.getValue(index)+1;
  //     return (
  // 	    s.pos_words > 0 &&
  // 	    s.neg_words <= 0 &&
  // 	    (s.table[0] & val) == val
  // 	    );
  //   }
  //   inline bool include(const BitSet& s) const
  //   {
  //     return (
  // 	    s.size() <= 2 &&//change that
  // 	    s.pos_words > 0 &&
  // 	    s.neg_words <= 0 &&
  // 	    (s.table[0] & (domainPool.getValue(index)+1)) == s.table[0]
  // 	    );
  //   }
  //   /*!
  //      Intersect its domain with a set s
  //   */
  //   inline void intersectTo( BitSet* s ) const
  //   {
  //     int i=s->pos_words;
  //     int j=s->neg_words;
  //     while( i-- > j )
  //       if( !i ) s->table[i] &= (domainPool.getValue(index)+1);
  //       else s->table[i] = 0;
  //   }
  //   /*!
  //      Do the union of its domain with a set s
  //   */
  //   inline void unionTo( BitSet* s ) const
  //   {
  //     if( s->pos_words > 0 && s->neg_words < 1 )
  //       s->table[0] |= (domainPool.getValue(index)+1);
  //   }
  //   /*!
  //      Intersect its domain with a set s
  //   */
  //   inline void intersectTo( BitSet& s ) const
  //   {
  //     int i=s.pos_words;
  //     int j=s.neg_words;
  //     while( i-- > j )
  //       if( !i ) s.table[i] &= (domainPool.getValue(index)+1);
  //       else s.table[i] = 0;
  //   }
  //   /*!
  //      Do the union of its domain with a set s
  //   */
  //   inline void unionTo( BitSet& s ) const
  //   {
  //     if( s.pos_words > 0 && s.neg_words < 1 )
  //       s.table[0] |= (domainPool.getValue(index)+1);
  //   }
  //   //@}
    

  //   /*!@name Domain handling*/
  //   //@{
  //   /*!
  //       Remove the value "v"
  //   */
  //   bool remove(const int v);
  //   /*!
  //       Remove all values but "v"
  //   */
  //   bool setDomain(const int v);
  //   /*!
  //       Remove all values strictly lower than l
  //   */
  //   bool setMin(const int l);
  //   /*!
  //       Remove all values strictly greater than u
  //   */
  //   bool setMax(const int u);
  //   /*!
  //       Remove all values that do not appear in the array "a" of length "l"
  //   */
  //   bool setDomain(const int* a, const int l);
  //   /*!
  //       Remove all values that do not appear in the {@link BitSet set} "s"
  //   */
  //   bool setDomain(const BitSet& s);
  //   /*!
  //       Remove all values that do not appear in the current domain of the variable "x"
  //   */
  //   bool setDomain(VariableInt*) const;
  //   /*!
  //       Remove all values that belong to the {@link BitSet set} "s"
  //   */
  //   bool removeSet(const BitSet& s);
  //   /*!
  //       Remove all values in the interval [l..u]
  //   */
  //   bool removeRange(const int l, const int u);
  //   //@}

  //   /// The mandatory methods
  //   /*!@name Decision methods*/
  //   //@{
  //   bool revise( BitSet*, int*, VariableInt* );
  //   bool revise( Constraint*, const int );
  //   //@}

  //   /*!@name Miscellaneous*/
  //   //@{
  //   /// Print on out stream
  //   virtual void print(std::ostream& o) const;
  //   virtual void printshort(std::ostream&) const ;
  //   virtual void printDomain(std::ostream&) const ;
  //   virtual int getType() const { return VariableInt::BOOL; };
  //   //@}
  // };




  class Convertor {

  public:

    Convertor() {}
    virtual ~Convertor() {}

    virtual int forward( const int v ) = 0;
    virtual int backward( const int v ) = 0;
    virtual void forward( BitSet& s ) = 0;
    virtual void backward( BitSet& s ) = 0;

    virtual int minPosAbs( VariableInt *x ) const = 0;
    virtual int minNegAbs( VariableInt *x ) const = 0;
    virtual int min( VariableInt *x ) const = 0;
    virtual int max( VariableInt *x ) const = 0;
    virtual int minCapacity( VariableInt *x ) const = 0;
    virtual int maxCapacity( VariableInt *x ) const = 0;
    virtual int random( VariableInt *x ) const = 0;
    virtual int domsize( VariableInt *x ) const = 0;
    virtual bool isRange( VariableInt *x ) const = 0;
    virtual bool included( VariableInt *x, const int lo, const int up ) const = 0;
    virtual bool intersect( VariableInt *x, const int lo, const int up ) const = 0;
    virtual bool include( VariableInt *x, const VariableInt *y ) const = 0;
    virtual bool included( VariableInt *x, const VariableInt *y ) const = 0;
    virtual bool intersect( VariableInt *x, const VariableInt *y ) const = 0;
    virtual bool include( VariableInt *x, const BitSet& s ) const = 0;
    virtual bool included( VariableInt *x, const BitSet& s ) const = 0;
    virtual bool intersect( VariableInt *x, const BitSet& s ) const = 0;
    virtual void unionTo( VariableInt *x, BitSet& s ) = 0;
    virtual void copyTo( VariableInt *x, BitSet& s ) = 0;
    virtual void intersectTo( VariableInt *x, BitSet& s ) = 0;
    virtual bool setMin( VariableInt *x, const int b ) = 0;
    virtual bool setMax( VariableInt *x, const int b ) = 0;
    virtual bool setDomain( VariableInt *x, const VariableInt *y ) = 0;
    virtual bool removeRange( VariableInt *x, const int lo, const int up ) = 0;

    virtual void print(std::ostream& o, const VariableInt* x) const = 0;
  };



  /********************************************
   * Virtual Iterator
   ********************************************/
  /*! \class VirtualIterator
    \brief Values Iteration Accessors.
  */
  class VirtualIterator : public DomainIterator {

  public:

    Convertor *conversion;
    //DomainIterator *reference;
    VariableInt *reference;
    int curval;

    /*!@name Constructors*/
    //@{
    VirtualIterator() {}
    void init(Convertor *c, VariableInt *r)
    {
      conversion = c;
      reference = r;
      curval = r->first();
    }
    virtual ~VirtualIterator() 
    {
    }
    //@}

    /*!@name Accessors*/
    //@{
    virtual bool next()
    {
      //reference->next();
      reference->setNext( curval );
      return curval != NOVAL;
    }
    virtual operator const int() const
    {
      //return conversion->forward( (int)(*reference) );
      return conversion->forward( curval );
    }
    //@}
  };


  /*********************************************************
   * VariableVirtual || Virtual variable
   *********************************************************/
  /*! \class VariableVirtual
    \brief Implementation of Virtual Variable ('View')
  */
  class VariableVirtual : public VariableInt {
  public:

    /*!@name Parameters*/
    //@{
    int lb_;
    int ub_;
    BitSet util;
    VariableInt *reference;
    Convertor *conversion;
    //@}

    /*!@name Constructors*/
    //@{
    VariableVirtual(Solver *s); 
    virtual ~VariableVirtual();
    virtual int& getIntDomain();
    //@}

    /*!@name Accessors and Iterators*/
    //@{

    virtual MistralList<Constraint*>* triggerOnValue() {return reference->triggerOnValue();}
    virtual MistralList<Constraint*>* triggerOnRange() {return reference->triggerOnRange();}
    virtual MistralList<Constraint*>* triggerOnDomain() {return reference->triggerOnDomain();}

    /// access to the domain iterator
    DomainIterator *begin();
    virtual VariableInt* getVar() { return reference; }

    /// Return the first (equal to min) value in the domain
    inline int first() const
    {
      return conversion->forward( reference->first() );
    }
    //   /// Return the last (equal to max) value in the domain
    //   inline int  last() const
    //   {
    //     return conversion->forward( reference->last() );
    //   }
    /// Return the smallest value currently in the domain that is strictly greater than "v"
    inline int getNext(const int v) const
    {
      return conversion->forward( reference->getNext( conversion->backward( v ) ) );
    }
    ///*! 
    //     Used when looping on all values in the domain, return
    //     true if "v" is valid value
    //   */
    //   inline bool good(const int v) const
    //     {
    //       return reference->good( conversion->backward( v ) );
    //     }
    /*! 
      Set "v" to the smallest value currently in the domain that is strictly greater than "v"
    */
    inline bool setNext(int& v) const
    {
      int v_ref = conversion->backward( v );
      bool res = reference->setNext( v_ref );
      v = conversion->forward( v_ref );
      return res;
    }

    /// Returns the minimum absolute value in [0..infty] \\inter domain, NOVAL if there are none
    inline int minPosAbs() const
    {
      //std::cerr << "not implemented" << std::endl;
      //exit(0);
      return conversion->minPosAbs( reference );
    }

    /// Returns the minimum absolute value in [-infty..0] \\inter domain
    inline int minNegAbs() const
    {
      //std::cerr << "not implemented" << std::endl;
      //exit(0);
      return conversion->minNegAbs( reference );
    }

    /*!
      Whether "v" is currently contained in the domain
    */
    inline bool contain(const int v) const
    {
      return reference->contain( conversion->backward( v ) );
    }
    inline bool fastContain(const int v) const
    {
      return reference->fastContain( conversion->backward( v ) );
    }

    /// Returns the assigned value if it exists
    inline int value() const
    {
      int val = reference->value();
      if(val != NOVAL)
	val = conversion->forward( val );
      return val;
    }

    /// Returns the minimum value in the domain
    inline int min() const
    {
      return conversion->min( reference );
    }

    /// Returns the maximum value in the domain
    inline int max() const
    {
      return conversion->max( reference );
    }

    /// Returns the minimum value that could belong to the domain
    inline int minCapacity() const
    {
      return conversion->minCapacity( reference );
    }

    /// Returns 1 + the maximum value that could belong to the domain
    inline int maxCapacity() const
    {
      return conversion->maxCapacity( reference );
    }

    // Returns a random value in the domain
    inline int random() const
    {
      return conversion->random( reference );
    }

    /// Returns the domain size
    inline int domsize() const
    {
      return conversion->domsize( reference );
    }
    //@}

    /*!@name Query methods*/
    //@{
    /// Whether or not the Variable is currently an interval
    inline bool isRange() const
    {
      return conversion->isRange( reference );
    }
    /// Whether or not the Variable is bound to a ground value
    inline bool isGround() const
    {
      return reference->isGround();
    }
    /// Whether or not the Variable is bound to a given ground value
    inline bool equal(const int v) const 
    { 
      return reference->equal( conversion->backward( v ) );
    }
    /*!
      Whether the domain has a nonempty intersection
      with the interval [l..u]
    */
    inline bool intersect(int lo, int up) const
    {
      return conversion->intersect( reference, lo, up );
    }
    inline bool included(int lo, int up) const
    {
      return conversion->included( reference, lo, up );
    }
    /*!
      Whether the domain has a nonempty intersection
      with the Variable x
    */
    inline bool intersect(const VariableInt* v) const
    {
      BitSet local_util(lb_, ub_, BitSet::empt);

      reference->unionTo( local_util );
      //local_util.print( std::cout );

      conversion->forward( local_util );
      //conversion->backward( local_util );
      //local_util.print( std::cout );

      return v->intersect( local_util );

      //reference->unionTo( util );
      //return true;

      //return conversion->intersect( reference, v );
    }
    /*!
      Whether the domain is included
      in the Variable x
    */
    inline bool included(const VariableInt* v) const
    {
      return conversion->included( reference, v );
    }
    inline bool include(const VariableInt* v) const
    {
      return conversion->include( reference, v );
    }
    /*!
      Whether the domain has a nonempty intersection
      with the set s
    */
    inline bool intersect(const BitSet* s) const
    {
      return conversion->intersect( reference, *s );
    }
    /*!
      Whether the domain is included
      in the set s
    */
    inline bool included(const BitSet* s) const
    {
      return conversion->included( reference, *s );
    }
    inline bool include(const BitSet* s) const
    {
      return conversion->include( reference, *s );
    }
    /*!
      Whether the domain has a nonempty intersection
      with the set s
    */
    inline bool intersect(const BitSet& s) const
    {
      return conversion->intersect( reference, s );
    }
    /*!
      Whether the domain has a nonempty intersection
      with the set s
    */
    inline bool intersect(const BitSet& s, int& r) const
    {
      std::cerr << "TODO" << std::endl;
      exit(0);
    }
    inline bool wordIntersect(const BitSet& s) const
    {
      return conversion->intersect( reference, s );
    }
    /*!
      Whether the domain is included
      in the set s
    */
    inline bool included(const BitSet& s) const
    {
      return conversion->included( reference, s );
    }
    inline bool include(const BitSet& s) const
    {
      return conversion->include( reference, s );
    }
    /*!
      Intersect its domain with a set s
    */
    inline void intersectTo( BitSet* s ) const
    {
      conversion->intersectTo( reference, *s );
    }
    /*!
      Do the union of its domain with a set s
    */
    inline void unionTo( BitSet* s ) const
    {
      conversion->unionTo( reference, *s );
    }
    /*!
      Do the union of its domain with a set s
    */
    inline void copyTo( BitSet* s ) const
    {
      conversion->copyTo( reference, *s );
    }
    /*!
      Intersect its domain with a set s
    */
    inline void intersectTo( BitSet& s ) const
    {
      conversion->intersectTo( reference, s );
    }
    /*!
      Do the union of its domain with a set s
    */
    inline void unionTo( BitSet& s ) const
    {
      conversion->unionTo( reference, s );
    }
    /*!
      Do the union of its domain with a set s
    */
    inline void copyTo( BitSet& s ) const
    {
      conversion->copyTo( reference, s );
    }
    //@}
    

    /*!@name Domain handling*/
    //@{
    /*!
      Remove the value "v"
    */
    bool remove(const int v)
    {
      return reference->remove( conversion->backward( v ) );
    }
    /*!
      Remove all values but "v"
    */
    bool setDomain(const int v)
    {
      int v_prime = conversion->backward( v );

      //print( std::cout );
      //std::cout << " <- " << v << std::endl;
      //reference->print( std::cout );
      //std::cout << " <- " << v_prime << std::endl;

      bool result = reference->setDomain( v_prime );
    
      //     print( std::cout );
      //     std::cout << std::endl;
      //     reference->print( std::cout );
      //     std::cout << std::endl;
      //     std::cout << result << std::endl;
      //     //exit(0);
    
      return result;
    }
    /*!
      Remove all values strictly lower than l
    */
    bool setMin(const int l)
    {
      return conversion->setMin( reference, l );
    }
    /*!
      Remove all values strictly greater than u
    */
    bool setMax(const int u)
    {
      return conversion->setMax( reference, u );
    }
    /*!
      Remove all values that do not appear in the array "a" of length "l"
    */
    bool setDomain(const int* a, const int l)
    {
      int i, tmp[l];
      for(i=0; i<l; ++i)
	tmp[i] = conversion->backward( a[i] );
      return reference->setDomain( tmp, l );
    }
    /*!
      Remove all values that do not appear in the {@link BitSet set} "s"
    */
    bool setDomain(const BitSet& s)
    {
      if(util.table) {
	util.copy(s);
	conversion->backward( util );
	return reference->setDomain( util );
      } else {
	int lb = s.next( min()-1 );
	int ub = s.prev( max()+1 );
	
	return(setMin(lb) && setMax(ub));
      }
    }
    /*!
      Remove all values that do not appear in the current domain of the variable "x"
    */
    bool setDomain(VariableInt* x) const 
    {

      //BitSet local_util(lb_, ub_, BitSet::full);
      //local_util.print( std::cout );
      //local_util.clear();
      //local_util.print( std::cout );
      BitSet local_util(lb_, ub_, BitSet::empt);

      reference->unionTo( local_util );
      //local_util.print( std::cout );

      conversion->forward( local_util );
      //conversion->backward( local_util );
      //local_util.print( std::cout );

      return x->setDomain( local_util );
      //x->print( std::cout );
    
      //reference->setDomain( util );
      //conversion->setDomain( reference, x );
    }
    /*!
      Remove all values that belong to the {@link BitSet set} "s"
    */
    bool removeSet(const BitSet& s)
    {
      if(util.table) {
	util.copy(s);
	conversion->backward( util );
	return reference->removeSet( util );
      } else {
	int lb = min();
	int ub = max();
	while( s.member(lb) ) ++lb; 
	while( s.member(ub) ) --ub; 
	return(setMin(lb) && setMax(ub));
      }
    }
    /*!
      Remove all values in the interval [l..u]
    */
    bool removeRange(const int l, const int u)
    {
      return conversion->removeRange( reference, l, u );
    }
    //@}

    /// The mandatory methods
    /*!@name Decision methods*/
    //@{
    bool revise( BitSet* s, int* i, VariableInt* v )
    {
      std::cerr << "TODO" << std::endl;
      exit(0);
    }
    bool revise( Constraint*, const int )
    {
      std::cerr << "TODO" << std::endl;
      exit(0);
    }
    //@}

    /*!@name Miscellaneous*/
    //@{
    /// Print on out stream
    virtual void print(std::ostream& o) const
    {
      printshort(o);
      if( isGround() )
	o << " = " << value();
      else {
	o << " = ";
	conversion->print( o, reference );
      }
    }
    virtual void printshort(std::ostream& o) const 
    {
      o << "v" << id+1 ;
    }
    virtual void printDomain(std::ostream& o) const 
    {
      //     DomainIterator *di = begin();
      //     o << "{";
      //     do {
      //       o << (*di) << ",";
      //     } while( ++(*di) );
      //     o << "}";
    }
    virtual int getType() const { return VariableInt::VIRTUAL; }
    //@}
  };


  class MapOffset : public Convertor {

  public:

    int K;

    MapOffset( const int o, VariableVirtual* x ) { 
      K = o;
      //x->reference->print( std::cout );    
      //std::cout << std::endl ;//<< (x->reference->min()+(K<0 ? K : 0)) << " "
      x->lb_ = x->reference->min()+(K<0 ? K : 0);
      x->ub_ = x->reference->max()+(K>0 ? K : 0);
      x->util.init( x->lb_, x->ub_, BitSet::empt );
      //x->util.init( x->reference->min()+(K<0 ? K : 0), x->reference->max()+(K>0 ? K : 0), BitSet::empt );
    }
    virtual ~MapOffset() {}

    virtual int forward( const int v ) { return v+K; }
    virtual int backward( const int v ) { return v-K; }
    virtual void forward( BitSet& s ) { s.increment(K); }
    virtual void backward( BitSet& s ) { s.decrement(K); }

    virtual int minPosAbs( VariableInt *x ) const
    {
      int ub = x->max(), lb = x->min();
      if(ub+K >= 0) {
	if(lb < -K) lb = -K;
	for(int v=lb; v<=ub; ++v)
	  if( x->contain(v) ) return (v+K);
      }
      return NOVAL;
    }
    virtual int minNegAbs( VariableInt *x ) const
    {
      int ub = x->max(), lb = x->min();
      if(lb-K <= 0) {
	if(ub > -K) ub = -K;
	for(int v=ub; v>=lb; --v)
	  if( x->contain(v) ) return (v+K);
      }
      return NOVAL;
    }
    virtual int min( VariableInt *x ) const { return x->min()+K; }
    virtual int max( VariableInt *x ) const { return x->max()+K; }
    virtual int minCapacity( VariableInt *x ) const { return x->minCapacity()+K; }
    virtual int maxCapacity( VariableInt *x ) const { return x->maxCapacity()+K; }
    virtual int random( VariableInt *x ) const { return x->random()+K; }
    virtual int domsize( VariableInt *x ) const { return x->domsize(); }
    virtual bool isRange( VariableInt *x ) const { return true; } //return x->isRange(); }
    virtual bool included( VariableInt *x, const int lo, const int up ) const { return x->included(lo-K, up-K); }
    virtual bool intersect( VariableInt *x, const int lo, const int up ) const { return x->intersect(lo-K, up-K); }
    virtual bool include( VariableInt *x, const VariableInt *y ) const { std::cerr << "TODO 1" << std::endl; exit(0); }
    virtual bool included( VariableInt *x, const VariableInt *y ) const { std::cerr << "TODO 2" << std::endl; exit(0); }
    virtual bool intersect( VariableInt *x, const VariableInt *y ) const { 
      //x->util.copy( )
      std::cerr << "TODO 3" << std::endl; exit(0); 
    }
    virtual bool include( VariableInt *x, const BitSet& s ) const { std::cerr << "TODO 4" << std::endl; exit(0); }
    virtual bool included( VariableInt *x, const BitSet& s ) const { std::cerr << "TODO 5" << std::endl; exit(0); }
    virtual bool intersect( VariableInt *x, const BitSet& s ) const { std::cerr << "TODO 6" << std::endl; exit(0); }
    virtual void unionTo( VariableInt *x, BitSet& s ) {  

      //  x->print( std::cout );
      // std::cout << std::endl;
      DomainIterator *iter = x->begin();
      do {
	//  std::cout << "* " << (*iter) << std::endl;
	// std::cout << "  " << forward(*iter) << std::endl;
	s.insert( forward(*iter) );
      } while( iter->next() );

      //     s.print( std::cout );
      //  std::cout << std::endl;
      //     //exit(0);
    }
    virtual void copyTo( VariableInt *x, BitSet& s ) { std::cerr << "TODO 7.5" << std::endl; exit(0); }
    virtual void intersectTo( VariableInt *x, BitSet& s ) { std::cerr << "TODO 8" << std::endl; exit(0); }
    virtual bool setMin( VariableInt *x, const int b ) { return x->setMin(b-K); }
    virtual bool setMax( VariableInt *x, const int b ) { return x->setMax(b-K); }
    virtual bool setDomain( VariableInt *x, const VariableInt *y ) { std::cerr << "TODO 9" << std::endl; exit(0); }
    virtual bool removeRange( VariableInt *x, const int lo, const int up ) { return x->removeRange(lo-K, up-K); }

    virtual void print(std::ostream& o, const VariableInt* x) const;
  };



  class MapFactor : public Convertor {

  public:

    int K;

    MapFactor( const int o, VariableVirtual* x ) { 
      K = o;
    
      //     x->reference->print( std::cout );
      //     std::cout << " * " << K << std::endl;

      int bound[2], i=(K<0);
      bound[i] = x->reference->min();
      if( bound[i] < 0 )
	bound[i] *= K;
      bound[1-i] = x->reference->max();
      if( bound[1-i] > 0 )
	bound[1-i] *= K;


      //    std::cout << "[" << bound[0] << ".." << bound[1] << "]" << std::endl;

      x->lb_ = bound[0];
      x->ub_ = bound[1];
      //x->util.init( bound[0], bound[1], BitSet::empt );
    }
    virtual ~MapFactor() {}

    // from the pointed variable to the view
    virtual int forward( const int v ) { return v*K; }
    // from the view to the pointed variable
    virtual int backward( const int v ) { return v/K; }
    virtual void forward( BitSet& s ) { /*TODO*/ }
    virtual void backward( BitSet& s ) { /*TODO*/ }
  
    virtual int minPosAbs( VariableInt *x ) const
    {
      //     int ub = x->max(), lb = x->min();
      //     if(ub+K >= 0) {
      //       if(lb < -K) lb = -K;
      //       for(int v=lb; v<=ub; ++v)
      // 	if( x->contain(v) ) return (v+K);
      //     }
      //     return NOVAL;
      std::cerr << "TODO" << std::endl;
      exit(0);
    }
    virtual int minNegAbs( VariableInt *x ) const
    {
      //     int ub = x->max(), lb = x->min();
      //     if(lb-K <= 0) {
      //       if(ub > -K) ub = -K;
      //       for(int v=ub; v>=lb; --v)
      // 	if( x->contain(v) ) return (v+K);
      //     }
      //     return NOVAL;
      std::cerr << "TODO" << std::endl;
      exit(0);
    }
    virtual int min( VariableInt *x ) const {     
      if( K>0 )
	return x->min()*K;
      else
	return x->max()*K;
    }
    virtual int max( VariableInt *x ) const { 
      if( K>0 )
	return x->max()*K;
      else
	return x->min()*K;
    }
    virtual int minCapacity( VariableInt *x ) const { 
      if( K>0 )
	return x->minCapacity()*K;
      else
	return x->maxCapacity()*K;
    }
    virtual int maxCapacity( VariableInt *x ) const { 
      if( K>0 )
	return x->maxCapacity()*K;
      else
	return x->minCapacity()*K;
    }
    virtual int random( VariableInt *x ) const { return x->random()*K; }
    virtual int domsize( VariableInt *x ) const { return x->domsize(); }
    virtual bool isRange( VariableInt *x ) const { return true; }
    virtual bool included( VariableInt *x, const int lo, const int up ) const { 
      if( K>0 ) {
	int lo_prime = lo/K;
	if(lo%K)
	  ++lo_prime;
	return x->included(lo_prime, up/K); 
      } else {
	int up_prime = up/K;
	if(up%K)
	  ++up_prime;
	return x->included(up_prime, lo/K);
      }
    }
    virtual bool intersect( VariableInt *x, const int lo, const int up ) const { 
      if( K>0 ) {
	int lo_prime = lo/K;
	if(lo%K)
	  ++lo_prime;
	return x->intersect(lo_prime, up/K); 
      } else {
	int up_prime = up/K;
	if(up%K)
	  ++up_prime;
	return x->intersect(up_prime, lo/K);
      }
    }
    virtual bool include( VariableInt *x, const VariableInt *y ) const { std::cerr << "TODO" << std::endl; exit(0); }
    virtual bool included( VariableInt *x, const VariableInt *y ) const { std::cerr << "TODO" << std::endl; exit(0); }
    virtual bool intersect( VariableInt *x, const VariableInt *y ) const { std::cerr << "TODO" << std::endl; exit(0); }
    virtual bool include( VariableInt *x, const BitSet& s ) const { std::cerr << "TODO" << std::endl; exit(0); }
    virtual bool included( VariableInt *x, const BitSet& s ) const { std::cerr << "TODO" << std::endl; exit(0); }
    virtual bool intersect( VariableInt *x, const BitSet& s ) const { std::cerr << "TODO" << std::endl; exit(0); }
    virtual void unionTo( VariableInt *x, BitSet& s ) { 

      //     x->print( std::cout );
      //     std::cout << std::endl;
      DomainIterator *iter = x->begin();
      do {
	//       std::cout << "* " << (*iter) << std::endl;
	//       std::cout << "  " << forward(*iter) << std::endl;
	s.insert( forward(*iter) );
      } while( iter->next() );

      //     s.print( std::cout );
      //     std::cout << std::endl;
      //     exit(0);
    }
    virtual void copyTo( VariableInt *x, BitSet& s ) { std::cerr << "TODO" << std::endl; exit(0); }
    virtual void intersectTo( VariableInt *x, BitSet& s ) { std::cerr << "TODO" << std::endl; exit(0); }
    virtual bool setMin( VariableInt *x, const int b ) 
    { 
      int b_prime = b/K;
      if( K>0 ) {      
	if( b>0 && b%K )
	  ++b_prime;

	//         std::cout << "(1) i.e. ";
	//         x->print( std::cout );
	//         std::cout << " >= " << ((double)b/(double)K) << " - " << b_prime << std::endl;

	bool result = x->setMin(b_prime);
	//std::cout << result << std::endl;
	return result;
      } else {
	if( b>0 && b%K )
	  --b_prime;
      
	//       std::cout << "(2) i.e. ";
	//       x->print( std::cout );
	//       std::cout << " <= " << ((double)b/(double)K) << " - " << b_prime << std::endl;

	bool result = x->setMax(b_prime);
	//       std::cout << result << std::endl;      
	return result;
      }
    }
    virtual bool setMax( VariableInt *x, const int b ) 
    { 
      int b_prime = b/K;
      if( K>0 ) {
	if( b<0 && b%K )
	  --b_prime;

	//       std::cout << "(3) i.e. ";
	//       x->print( std::cout );
	//       std::cout << " <= " << ((double)b/(double)K) << " - " << b_prime << std::endl;

	bool result = x->setMax(b_prime);
	//      std::cout << result << std::endl;
	return result;
      } else {
	if( b<0 && b%K )
	  ++b_prime;

	//       std::cout << "(4) i.e. ";
	//       x->print( std::cout );
	//       std::cout << " >= " << ((double)b/(double)K) << " - " << b_prime << std::endl;

	bool result = x->setMin(b_prime);
	//      std::cout  << result << std::endl;
	return result;
      }
    }
    virtual bool setDomain( VariableInt *x, const VariableInt *y ) { std::cerr << "TODO" << std::endl; exit(0); }
    virtual bool removeRange( VariableInt *x, const int lo, const int up ) { 
      int lo_prime = lo/K;
      int up_prime = up/K;

      if( K>0 ) {
	if(lo>0 && lo%K)
	  ++lo_prime;
	if(up<0 && up%K)
	  --up_prime;
	return x->removeRange(lo_prime, up_prime); 
      } else {
	if(up<0 && up%K)
	  ++up_prime;
	if(lo>0 && lo%K)
	  --lo_prime;
	return x->removeRange(up_prime, lo_prime);
      }
    }

    virtual void print(std::ostream& o, const VariableInt* x) const;
  };



  class MapQuotient : public Convertor {

  public:

    int K;

    MapQuotient( const int o, VariableVirtual* x ) { 
      K = o;
      int bound[2], i=(K<0);
      bound[i] = x->reference->min();
      if( bound[i] > 0 )
	bound[i] /= K;
      bound[1-i] = x->reference->max();
      if( bound[1-i] < 0 )
	bound[1-i] /= K;

      x->lb_ = bound[0];
      x->ub_ = bound[1];
      x->util.init( bound[0], bound[1], BitSet::empt );
    }
    virtual ~MapQuotient() {}

    virtual int forward( const int v ) { return v/K; }
    virtual int backward( const int v ) { return v*K; }
    virtual void forward( BitSet& s ) { /*TODO*/ }
    virtual void backward( BitSet& s ) { /*TODO*/ }
  
    virtual int minPosAbs( VariableInt *x ) const
    {
      //     int ub = x->max(), lb = x->min();
      //     if(ub+K >= 0) {
      //       if(lb < -K) lb = -K;
      //       for(int v=lb; v<=ub; ++v)
      // 	if( x->contain(v) ) return (v+K);
      //     }
      return NOVAL;
    }
    virtual int minNegAbs( VariableInt *x ) const
    {
      //     int ub = x->max(), lb = x->min();
      //     if(lb-K <= 0) {
      //       if(ub > -K) ub = -K;
      //       for(int v=ub; v>=lb; --v)
      // 	if( x->contain(v) ) return (v+K);
      //     }
      return NOVAL;
    }
    virtual int min( VariableInt *x ) const { return x->min()/K; }
    virtual int max( VariableInt *x ) const { return x->max()/K; }
    virtual int minCapacity( VariableInt *x ) const { return x->minCapacity()/K; }
    virtual int maxCapacity( VariableInt *x ) const { return x->maxCapacity()/K; }
    virtual int random( VariableInt *x ) const { return x->random()/K; }
    virtual int domsize( VariableInt *x ) const { return x->domsize(); }
    virtual bool isRange( VariableInt *x ) const { return true; }
    virtual bool included( VariableInt *x, const int lo, const int up ) const { 
      if( K>0 )
	return x->included(lo*K, up*K); 
      else
	return x->included(up*K, lo*K);
    }
    virtual bool intersect( VariableInt *x, const int lo, const int up ) const { 
      if( K>0 )
	return x->intersect(lo*K, up*K); 
      else
	return x->intersect(up*K, lo*K);
    }
    virtual bool include( VariableInt *x, const VariableInt *y ) const { std::cerr << "TODO" << std::endl; exit(0); }
    virtual bool included( VariableInt *x, const VariableInt *y ) const { std::cerr << "TODO" << std::endl; exit(0); }
    virtual bool intersect( VariableInt *x, const VariableInt *y ) const { std::cerr << "TODO" << std::endl; exit(0); }
    virtual bool include( VariableInt *x, const BitSet& s ) const { std::cerr << "TODO" << std::endl; exit(0); }
    virtual bool included( VariableInt *x, const BitSet& s ) const { std::cerr << "TODO" << std::endl; exit(0); }
    virtual bool intersect( VariableInt *x, const BitSet& s ) const { std::cerr << "TODO" << std::endl; exit(0); }
    virtual void unionTo( VariableInt *x, BitSet& s ) { std::cerr << "TODO" << std::endl; exit(0); }
    virtual void copyTo( VariableInt *x, BitSet& s ) { std::cerr << "TODO" << std::endl; exit(0); }
    virtual void intersectTo( VariableInt *x, BitSet& s ) { std::cerr << "TODO" << std::endl; exit(0); }
    virtual bool setMin( VariableInt *x, const int b ) 
    { 
      if( K>0 )
	return x->setMin(b*K);
      else
	return x->setMax(b*K);
    }
    virtual bool setMax( VariableInt *x, const int b ) 
    { 
      if( K>0 )
	return x->setMax(b*K);
      else
	return x->setMin(b*K);
    }
    virtual bool setDomain( VariableInt *x, const VariableInt *y ) { std::cerr << "TODO" << std::endl; exit(0); }
    virtual bool removeRange( VariableInt *x, const int lo, const int up ) { 
      if( K>0 )
	return x->removeRange(lo*K, up*K); 
      else
	return x->removeRange(up*K, lo*K);
    }

    virtual void print(std::ostream& o, const VariableInt* x) const;
  };


  class ConstantIterator : public DomainIterator {

  public: 

    virtual ~ConstantIterator() {}

    ConstantIterator( const int v )
    {
      curval = v;
    }
    ConstantIterator( )
    {
    }

    int curval;

    inline bool next() { return false; }
    inline operator const int() const { return curval; }
  };

  /*********************************************************
   * Constant 
   *********************************************************/
  /*! \class Constant
    \brief Implementation of a Constant
  */
  class Constant : public VariableInt {

  public:
  
    int val;

    Constant(Solver *s, const int v );
    virtual ~Constant(); 
    virtual int& getIntDomain();

    /*!@name Accessors and Iterators*/
    //@{
    /// access to the domain iterator
    DomainIterator *begin();

    /// Return the first (equal to min) value in the domain
    inline int first() const 
    {
      return val;
    }
    /// Return the last (equal to max) value in the domain
    inline int  last() const 
    { 
      return val; 
    }
    /// Return the smallest value currently in the domain that is strictly greater than "v"
    inline int getNext(const int v) const 
    {
      return NOVAL;
    }
    /*! Used when looping on all values in the domain, return    
      true if "v" is valid value 
    */
    inline bool good(const int v) const 
    { 
      return ( v != NOVAL );
    }
    /*! Set "v" to the smallest value currently in the domain that is strictly greater than "v"
     */
    inline bool setNext(int& v) const 
    { 
      return false; 
    }
    /// Returns the minimum absolute value in [0..infty] \\inter domain, NOVAL if there are none
    inline int minPosAbs() const 
    {
      return (val >= 0 ? val : NOVAL);
    }
    /// Returns the minimum absolute value in [-infty..0] \\inter domain, NOVAL if there are none
    inline int minNegAbs() const  
    {
      return (val <= 0 ? val : NOVAL);
    }
    //@}


    /*!@name Utilities*/
    //@{
    /*!
      Whether "v" is currently contained in the domain
    */
    inline bool contain(const int v) const 
    { 
      return ( v == val ); 
    }
    inline bool fastContain(const int v) const 
    { 
      return ( v == val ); 
    }
    /// Returns the assigned value if it exists
    inline int value() const
    {
      return val;
    }    
    /// Returns the minimum value in the domain
    inline int min() const
    {
      return val;
    }
    /// Returns the maximum value in the domain
    inline int max() const
    {
      return val;
    }  
    /// Returns the minimum value that could belong to the domain
    inline int minCapacity() const 
    {
      return val;
    }
    /// Returns 1 + the maximum value that could belong to the domain
    inline int maxCapacity() const 
    {
      return val+1;
    }
    /// Returns a random value in the domain
    inline int random() const
    {
      return val;
    }
    /// Returns the domain size
    inline int domsize() const
    {
      return 1;
    }
    /// Whether or not the Variable is currently an interval
    inline bool isRange() const 
    { 
      return true; 
    }
    /// Whether or not the Variable is bound to a ground value
    inline bool isGround() const
    {
      return true;   
    }
    /// Whether or not the Variable is bound to a given ground value
    inline bool equal(const int v) const { return (val == v); }
    /*!
      Whether the domain has a nonempty intersection
      with the interval [l..u]
    */
    inline bool intersect(int lo, int up) const
    {
      return ( lo <= val && val <= up );
    }
    /*!
      Whether the domain is included in
      the interval [l..u]
    */
    inline bool included(int lo, int up) const
    {
      return ( lo == up && up == val );
    }
    /*!
      Whether the domain has a nonempty intersection
      with the Variable x 
    */
    inline bool intersect(const VariableInt* v) const 
    {  
      return v->contain( val ) ;
    }
    /*!
      Whether the domain is included
      in the Variable x 
    */
    inline bool included(const VariableInt* v) const 
    {
      return (v->min() == val && v->max() == val);
    }
    inline bool include(const VariableInt* v) const 
    {
      return (v->contain( val));
    }
    /*!
      Whether the domain has a nonempty intersection
      with the set s 
    */
    inline bool intersect(const BitSet* s) const
    {
      return s->member( val );
    }
    /*!
      Whether the domain is included
      in the set s
    */
    inline bool included(const BitSet* s) const
    {
      return (s->member(val));
    }
    inline bool include(const BitSet* s) const
    {
      return (s->size() == 1 && s->member(val));
    }
    /*!
      Whether the domain has a nonempty intersection
      with the set s 
    */
    inline bool intersect(const BitSet& s) const
    {
      return s.member( val );
    }
    /*!
      Whether the domain has a nonempty intersection
      with the set s
    */
    inline bool intersect(const BitSet& s, int& r) const
    {
      return s.member( val );
    }
    inline bool wordIntersect(const BitSet& s) const
    {
      return s.wordMember( val );
    }
    /*!
      Whether the domain is included
      in the set s
    */
    inline bool included(const BitSet& s) const
    {
      return (s.member(val));
    }
    inline bool include(const BitSet& s) const
    {
      return (s.size() == 1 && s.member(val));
    }
    /*!
      Intersect its domain with a set s
    */
    inline void intersectTo( BitSet* s ) const
    {
      if( s->member( val ) ) {
	s->clear();
	s->insert( val );
      } else s->clear();
    }
    /*!
      Do the union of its domain with a set s
    */
    inline void unionTo( BitSet* s ) const
    {
      s->insert( val );
    }
    /*!
      Do the union of its domain with a set s
    */
    inline void copyTo( BitSet* s ) const
    {
      s->clear();
      s->insert( val );
    }
    /*!
      Intersect its domain with a set s
    */
    inline void intersectTo( BitSet& s ) const
    {
      if( s.member( val ) ) {
	s.clear();
	s.insert( val );
      } else s.clear();
    }
    /*!
      Do the union of its domain with a set s
    */
    inline void unionTo( BitSet& s ) const
    {
      s.insert( val );
    }
    /*!
      Do the union of its domain with a set s
    */
    inline void copyTo( BitSet& s ) const
    {
      s.clear();
      s.insert( val );
    }
    //@}
    

    /*!@name Domain handling*/
    //@{
    /*! 
      Remove the val "v"
    */
    bool remove(const int v) 
    {
      return ( v != val );
    }
    /*! 
      Remove all vals but "v"
    */
    bool setDomain(const int v)
    {
      return ( v == val );
    }
    /*! 
      Remove all vals strictly lower than l
    */
    bool setMin(const int l) 
    {
      return ( l <= val );
    }
    /*! 
      Remove all values strictly greater than u
    */
    bool setMax(const int u)
    {
      return ( u >= val );
    }
    /*! 
      Remove all values that do not appear in the array "a" of length "l"
    */
    bool setDomain(const int* a, const int l) 
    {
      int i = 0;
      while( i < l )
	if (a[i++] == val) return true;
      return false;
    }
    /*! 
      Remove all values that do not appear in the {@link BitSet set} "s"
    */
    bool setDomain(const BitSet& s)
    {
      return s.member( val ) ;
    }
    /*! 
      Remove all values that do not appear in the current domain of the variable "x"
    */
    bool setDomain(VariableInt* x) const
    {
      //return x->contain( val );
      return x->setDomain( val );
    }
    /*! 
      Remove all values that belong to the {@link BitSet set} "s"
    */
    bool removeSet(const BitSet& s)
    {
      return !s.member( val ) ;
    }
    /*! 
      Remove all values in the interval [l..u]
    */
    bool removeRange(const int l, const int u)
    {
      return (l > val || u < val );
    }
    //@}

    inline bool revise( BitSet* supports, 
			int* residues, 
			VariableInt* X)
    {
      return X->intersect( supports[val], residues[val] );
    }
    bool revise( Constraint* c, const int i)
    {
      return( ( c->firstSupport(i, val) || 
		c->findSupport(i, val) ) );
    }

    /*!@name Miscellaneous*/
    //@{
    /// Print on out stream
    virtual void print(std::ostream& o) const
    {
      o << val;
    }
    virtual void printshort(std::ostream& o) const 
    {
      o << val;
    }
    virtual void printDomain(std::ostream& o) const 
    {
      o << val;
    }
    virtual int getType() const { return VariableInt::CONST; };
    //@}
  };


  class RangeIterator : public DomainIterator {

  public: 

    virtual ~RangeIterator() {}

    int curval;
    int ubound;

    inline bool next() { return ++curval <= ubound; }
    inline operator const int() const { return curval; }
  };

  /********************************************************
   * VariableRange || An implementation of Variable using an interval
   ********************************************************/
  /*! \class VariableRange
    \brief Implementation of interval variable

    Only the bounds are stored, as reversible int
  */
  class VariableRange : public VariableInt {

  public:
    /*!@name Parameters*/
    //@{
    ReversibleNum<int> vmin;
    ReversibleNum<int> vmax;
    //@}

    /*!@name Constructors*/
    //@{
    VariableRange(Solver *s, const int v); 

    VariableRange(Solver *s, const int l, const int u); 

    //Initialize the variable to have id num and domain size "size"
    void setVariable(const int low, const int up);
    VariableRange(Solver *s); 

    virtual ~VariableRange();

    virtual int& getIntDomain();
    //@}

    /*!@name Accessors and Iterators*/
    //@{
    /// access to the domain iterator
    DomainIterator *begin();

    ///
    inline int whichBound()
    {
      //return( vmin.lvl_.back() > vmax.lvl_.back() );

      int curmean = (vmin + vmax);
      int oldmean = (vmin.val_[0] + vmax.val_[0]);

      return( curmean > oldmean );

    }
    /// Return the first (equal to min) value in the domain
    inline int first() const
    {
      return vmin;
    }
    /// Return the last (equal to max) value in the domain
    inline int  last() const
    {
      return vmax;
    }
    /// Return the smallest value currently in the domain that is strictly greater than "v"
    inline int getNext(const int v) const
    {
      return ( v+1 );
    }
    /*! Used when looping on all values in the domain, return
      true if "v" is valid value
    */
    inline bool good(const int v) const
    {
      return ( v <= vmax );
    }
    /*! Set "v" to the smallest value currently in the domain that is strictly greater than "v"
     */
    inline bool setNext(int& v) const
    {
      return ( ++v <= vmax );
    }
    /// Returns the minimum absolute value in [0..infty] \\inter domain, NOVAL if there are none
    inline int minPosAbs() const 
    {
      if( 0 > vmax ) return NOVAL;
      return ( 0 >= vmin ? 0 : (int)vmin );
    }
    /// Returns the minimum absolute value in [-infty..0] \\inter domain, NOVAL if there are none
    inline int minNegAbs() const  
    {
      if( 0 < vmin ) return NOVAL;
      return ( 0 <= vmax ? 0 : (int)vmax );
    }
    //@}


    /*!@name Utilities*/
    //@{
    inline bool contain(const int v) const
    {
      return ( v >= vmin && v <= vmax );
    }
    inline bool fastContain(const int v) const
    {
      return ( v >= vmin && v <= vmax );
    }
    /// Returns the assigned value if it exists
    inline int value() const
    {
      //return ( (int)vmin == (int)vmax ? vmin : NOVAL );    
      return (int)vmin;    
    }    
    /// Returns the minimum value in the domain
    inline int min() const
    {
      return vmin;
    }
    /// Returns the maximum value in the domain
    inline int max() const
    {
      return vmax;
    }  
    /// Returns the minimum value that could belong to the domain
    inline int minCapacity() const 
    {
      return vmin.val_[0];
    }
    /// Returns 1 + the maximum value that could belong to the domain
    inline int maxCapacity() const 
    {
      return vmax.val_[0]+1;
    }
    /// Returns a random value in the domain
    inline int random() const
    {
      //return vmin + (rand()%(vmax-vmin));
      return vmin + (randint(vmax-vmin));
    }
    /// Returns the domain size
    inline int domsize() const
    {
      return ( 1 + vmax - vmin );
    }
    /// Whether or not the Variable is currently an interval
    inline bool isRange() const 
    { 
      return true; 
    }
    /// Whether or not the Variable is bound to a ground value
    inline bool isGround() const
    {
      return ( (int)vmin == (int)vmax );   
    }
    /// Whether or not the Variable is bound to a given ground value
    inline bool equal(const int v) const { return ((int)vmin == (int)vmax && (int)vmin == v); }
    /*!
      Whether the domain has a nonempty intersection
      with the interval [l..u]
    */
    inline bool intersect(const int lo, const int up) const
    {
      return ( up >= vmin && lo <= vmax );
    }
    /*!
      Whether the domain is included in
      the interval [l..u]
    */
    inline bool included(const int lo, const int up) const
    {
      return ( lo >= vmin && up <= vmax );
    }
    /*!
      Whether the domain has a nonempty intersection
      with the Variable x
    */
    inline bool intersect(const VariableInt* v) const
    {
      return intersect(v->min(), v->max());
    }
    inline bool included(const VariableInt* v) const
    {
      return included(v->min(), v->max());
    }
    inline bool include(const VariableInt* v) const
    {
      return v->min() >= vmin && v->max() <= vmax;
    }
    inline bool intersect(const BitSet* s) const
    {
      return ( s->max() >= vmin && s->min() <= vmax );
    }
    inline bool included(const BitSet* s) const
    {
      return ( s->min() >= vmin && s->max() <= vmax );
    }
    inline bool include(const BitSet* s) const
    {
      return ( true );
    }
    inline bool intersect(const BitSet& s) const
    {
      return ( s.max() >= vmin && s.min() <= vmax );
    }
    inline bool intersect(const BitSet& s, int& r) const
    {
      return ( s.max() >= vmin && s.min() <= vmax );
    }
    inline bool wordIntersect(const BitSet& s) const
    {
      return ( s.max() >= vmin && s.min() <= vmax );
    }
    inline bool included(const BitSet& s) const
    {
      return ( s.min() >= vmin && s.max() <= vmax );
    }
    inline bool include(const BitSet& s) const
    {
      return ( true );
    }
    /*!
      Intersect its domain with a set s
    */
    inline void intersectTo( BitSet* s ) const
    {
      s->setMin( vmin );
      s->setMax( vmax );
    }
    /*!
      Do the union of its domain with a set s
    */
    inline void unionTo( BitSet* s ) const
    {
      s->addInterval(vmin, vmax);
    }
    /*!
      Do the union of its domain with a set s
    */
    inline void copyTo( BitSet* s ) const
    {
      s->clear();
      s->addInterval(vmin, vmax);
    }
    /*!
      Intersect its domain with a set s
    */
    inline void intersectTo( BitSet& s ) const
    {
      s.setMin( vmin );
      s.setMax( vmax );
    }
    /*!
      Do the union of its domain with a set s
    */
    inline void unionTo( BitSet& s ) const
    {
      s.addInterval(vmin, vmax);
    }
    /*!
      Do the union of its domain with a set s
    */
    inline void copyTo( BitSet& s ) const
    {
      s.clear();
      s.addInterval(vmin, vmax);
    }
    //@}

    /*!@name Domain handling*/
    //@{
    /*!
      Remove the value "v"
    */
    bool remove(const int v);
    /*!
      Remove all values but "v"
    */
    bool setDomain(const int v);
    /*!
      Remove all values strictly lower than l
    */
    bool setMin(const int l);
    /*!
      Remove all values strictly greater than u
    */
    bool setMax(const int u);
    /*!
      Remove all values that do not appear in the array "a" of length "l"
    */
    bool setDomain(const int* a, const int l);
    /*!
      Remove all values that do not appear in the {@link BitSet set} "s"
    */
    bool setDomain(const BitSet& s);
    /*!
      Remove all values that do not appear in the current domain of the variable "x"
    */
    bool setDomain(VariableInt*) const;
    /*!
      Remove all values that belong to the {@link BitSet set} "s"
    */
    bool removeSet(const BitSet& s);
    /*!
      Remove all values in the interval [l..u]
    */
    bool removeRange(const int l, const int u);
    //@}

    /// The mandatory methods
    /*!@name Decision methods*/
    //@{
    bool revise( BitSet*, int*, VariableInt* );
    bool revise( Constraint*, const int );
    //@}

    /*!@name Miscellaneous*/
    //@{
    /// Print on out stream
    virtual void print(std::ostream& o) const;
    virtual void printshort(std::ostream&) const ;
    virtual void printDomain(std::ostream&) const ;
    virtual int getType() const { return VariableInt::RANGE; };
    //@}
  };

};

#endif // __VAR_H
