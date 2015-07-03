
/*
  Mistral 2.0 is a constraint satisfaction and optimisation library
  Copyright (C) 2009  Emmanuel Hebrard
  
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

  The author can be contacted electronically at emmanuel.hebrard@gmail.com.
*/


/*! \file mistral_structure.hpp
  \brief Header file for the basic data-structures.
*/


#include <iostream>
#include <iomanip>
#include <stdlib.h> 
#include <sstream> 
#include <string.h>
#include <limits.h>

#include <mistral_global.hpp>


#ifndef __STRUCTURE_HPP
#define __STRUCTURE_HPP


namespace Mistral {


template <class WORD_TYPE>
void showUint(WORD_TYPE n, std::ostream& os) {
  WORD_TYPE mask=1;
  while(mask){
    if(mask & n) os<<1;
    else os<<0;
    mask = mask << 1;
  }
}

const int getlast[256] = {-1, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7};

const int NOVAL = (int)((~(unsigned int)0)/2);
#define INFTY  NOVAL/2


//   /**********************************************
//    * Vector
//    **********************************************/
//   /*! \class Vector
//     \brief Simple vector class     
//   */
  
//   template <class DATA_TYPE>
//   class Vector {
//   public:

//     /*!@name Parameters*/
//     //@{
//     DATA_TYPE* stack_;
//     unsigned int capacity;
//     unsigned int size;
//     //@}

//     /*!@name Constructor*/
//     //@{
//     Vector()
//     {
//       capacity = 0;
//       size = 0;
//       stack_ = NULL;
//     }

//     Vector(const Vector<DATA_TYPE>& s)
//     {
//       capacity = s.size;
//       stack_ = (DATA_TYPE*) malloc(capacity*sizeof(DATA_TYPE));
//       //stack_ = new DATA_TYPE[capacity];
//       for(size = 0; size<capacity; ++size)
// 	stack_[size] = s[size];
//     }

//     Vector(const int n)
//     {
//       capacity = n;
//       size = n;
//       if( capacity ) {
// 	stack_ = (DATA_TYPE*) malloc(capacity*sizeof(DATA_TYPE));
// 	//stack_ = new DATA_TYPE[capacity];
// 	int f = sizeof(DATA_TYPE)/sizeof(int);
// 	std::fill((int*)stack_, (int*)stack_+(capacity*f), 0);
//       }
//       else stack_ = NULL;
//     }
//     //@}

//     /*!@name Destructor*/
//     //@{  
//     virtual ~Vector()
//     {
// #ifdef _DEBUG_MEMORY
//       std::cout << "c delete vector: " << size << " " << capacity << " " // ;
//       // display(std::cout);
//       // std::cout 
// 	<< std::endl;
// #endif
//       free( stack_ );
//       //delete [] stack_;
//     }
//     //@}

//     /*!@name Initialisation*/
//     //@{
//     void initialise(const unsigned int c)
//     {
//       size = 0;
//       capacity = c;

//       //std::cout << "init(c): " << capacity << std::endl;

//       stack_ = (DATA_TYPE*) malloc(capacity*sizeof(DATA_TYPE));
      
//       // stack_ = new DATA_TYPE[capacity];
      
//       int f = sizeof(DATA_TYPE)/sizeof(int);
//       std::fill((int*)stack_, (int*)stack_+(capacity*f), 0);

//       //DATA_TYPE x(0);
//       //std::fill(stack_, stack_+capacity, x);
//     }

//     void initialise(const unsigned int s, const unsigned int c)
//     {
//       size = s;
//       capacity = c;
//       stack_ = (DATA_TYPE*) malloc(capacity*sizeof(DATA_TYPE));

//       //std::cout << "init(s,c): " << capacity << std::endl;

//       // stack_ = new DATA_TYPE[capacity];
      
//       int f = sizeof(DATA_TYPE)/sizeof(int);
//       std::fill((int*)stack_, (int*)stack_+(capacity*f), 0);
      
//       //DATA_TYPE x(0);
//       //std::fill(stack_, stack_+capacity, x);
//     }

//     void extendStack( const unsigned int l=0 )
//     {

//       //std::cout << "extend stack!! " << this << std::endl;

//       unsigned int increment = (l ? l : (capacity+1) << 1);
//       capacity += increment;

//       // DATA_TYPE* new_stack = (DATA_TYPE*) realloc(stack_, capacity*sizeof(DATA_TYPE));
//       // std::cout << (int*)new_stack << " " << (int*)stack_ << std::endl;
//       // stack_ = new_stack;


//       DATA_TYPE* new_stack = (DATA_TYPE*) malloc(capacity*sizeof(DATA_TYPE));
//       // for(int i=0; i<capacity-increment; ++i)
//       // 	new_stack[i] = stack_[i];

//       memcpy(new_stack, stack_, (capacity-increment)*sizeof(DATA_TYPE));

//       //memcpy((int*)new_stack, (int*)stack_, (capacity-increment)*f);
//       free(stack_); 
//       stack_ = new_stack;

//       // std::cout << "extend: " << capacity << std::endl;
     
//       // DATA_TYPE* new_stack_ = new DATA_TYPE[capacity];
//       // memcpy(new_stack_, stack_, capacity-increment);
//       // delete [] stack_; 
//       // stack_ = new_stack_;
      
//       int f = sizeof(DATA_TYPE)/sizeof(int);
//       std::fill((int*)stack_+(capacity-increment)*f, (int*)stack_+(capacity*f), 0);


//       //std::cout << " ==> " << this << std::endl;

//       //DATA_TYPE x(0);
//       //std::fill(stack_+capacity-increment, stack_+capacity, x);
//     }

//     void resize( const unsigned int l )
//     {
//       if( capacity < l ) {
// 	capacity = l;
// 	DATA_TYPE* new_stack = (DATA_TYPE*) realloc(stack_, capacity*sizeof(DATA_TYPE));
// 	stack_ = new_stack;

// 	// std::cout << "resize: " << l << std::endl;
	
// 	// DATA_TYPE* new_stack_ = new DATA_TYPE[l];
// 	// memcpy(new_stack_, stack_, capacity);
// 	// delete [] stack_;
// 	// stack_ = new_stack_;

// 	// capacity = l;
//       }
//       size = l;
//     }
//     //@}

//     /*!@name Accessors*/
//     //@{
//     inline int empty() const
//     {
//       return !size;
//     }
  
//     inline void add(DATA_TYPE x)
//     {
//       if( capacity == size ) 
// 	extendStack();
//       stack_[size++] = x;
//     }

//     inline void secure(const DATA_TYPE x) 
//     {
//       if( capacity == size ) 
// 	extendStack();
//     }

//     inline void fast_add(DATA_TYPE x)
//     {
//       stack_[size++] = x;
//     }

//     inline void push_back(DATA_TYPE x)
//     {
//       if( capacity == size ) 
// 	extendStack();
//       stack_[size++] = x;
//     }

//     inline DATA_TYPE pop_until(const unsigned int level)
//     {
//       size = level;
//       return stack_[size];
//     }

//     inline DATA_TYPE pop()
//     {
//       return stack_[--size];
//     }

//     inline void pop(DATA_TYPE& x)
//     {
//       x = stack_[--size];
//     }

//     inline void clear()
//     {
//       size = 0;
//     }

//     inline void remove(const unsigned int i)
//     {  
//       stack_[i] = stack_[--size];
//     }

//     inline void remove_elt(DATA_TYPE& elt)
//     {
//       unsigned int j=size;
//       while(j && stack_[--j] != elt);
//       stack_[j] = stack_[--size];
//     }

//     inline void setBack(const DATA_TYPE& x, const int k=1)
//     {
//       stack_[size-k] = x;
//     }

//     inline DATA_TYPE& front(const int k=0)
//     {
//       return stack_[k];
//     }

//     inline const DATA_TYPE front(const int k=0) const
//     {
//       return stack_[k];
//     }

//     inline DATA_TYPE& back(const int k=1)
//     {
//       return stack_[size-k];
//     }

//     inline const DATA_TYPE back(const int k=1) const
//     {
//       return stack_[size-k];
//     }

//     inline DATA_TYPE& operator[](const int i)
//     {
//       return stack_[i];
//     }

//     inline const DATA_TYPE operator[](const int i) const
//     {
//       return stack_[i];
//     }

//     inline Vector< DATA_TYPE >& operator=(const Vector< DATA_TYPE >& x)
//     {
//       initialise(0, x.capacity);
//       for(unsigned int i=0; i<x.size; ++i)
// 	add(x[i]);
//       return *this;
//     }
//     //@}

//     /*!@name Printing*/
//     //@{
//     std::ostream& display(std::ostream& os) const {
//       os << "[";
//       if(size) os << stack_[0] ;
//       for(unsigned int i=1; i<size; ++i)
// 	os << " " << stack_[i];
//       os << "]";
//       return os;
//     }
//     //@}

//   };



  /**********************************************
   * Vector
   **********************************************/
  /*! \class Vector
    \brief Simple vector class     
  */
  //int global = 0;
  
  template <class DATA_TYPE>
  int increasing_order(const void *x, const void *y) {
    DATA_TYPE& x_ = *((DATA_TYPE*)x);
    DATA_TYPE& y_ = *((DATA_TYPE*)y);
    return(x_ < y_ ? -1 : (x_ > y_ ? 1 : 0));
  }
  
  template <class DATA_TYPE>
  class Vector {
  public:

    typedef DATA_TYPE* iterator;
    
    /*!@name Parameters*/
    //@{
    DATA_TYPE* stack_;
    unsigned int capacity;
    unsigned int size;
    //@}

    /*!@name Constructor*/
    //@{
    Vector()
    {
      capacity = 0;
      size = 0;
      stack_ = NULL;
    }

    Vector(const Vector<DATA_TYPE>& s)
    {
      capacity = s.size;
      //stack_ = (DATA_TYPE*) malloc(capacity*sizeof(DATA_TYPE));
      stack_ = new DATA_TYPE[capacity];
      for(size = 0; size<capacity; ++size)
	stack_[size] = s[size];
    }

    Vector(const int n)
    {
		initialise(n);
	//       capacity = n;
	//       size = n;
	//       if( capacity ) {
	// //stack_ = (DATA_TYPE*) malloc(capacity*sizeof(DATA_TYPE));
	// stack_ = new DATA_TYPE[capacity];
	// //int f = sizeof(DATA_TYPE)/sizeof(int);
	// //std::fill((int*)stack_, (int*)stack_+(capacity*f), 0);
	// std::fill(stack_, stack_+capacity, (DATA_TYPE)0);
	//       }
	//       else stack_ = NULL;
    }
    //@}

    /*!@name Destructor*/
    //@{  
    virtual ~Vector()
    {
#ifdef _DEBUG_MEMORY
      std::cout << "c delete vector: " << size << " " << capacity << " " // ;
      // display(std::cout);
      // std::cout 
	<< std::endl;
#endif
      //free( stack_ );
      delete [] stack_;
    }
    //@}

    /*!@name Initialisation*/
    //@{
    void initialise(const unsigned int c)
    {
      size = 0;
      capacity = c;
      stack_ = new DATA_TYPE[capacity];
      std::fill(stack_, stack_+capacity, DATA_TYPE());
    }

    void initialise(const unsigned int s, const unsigned int c)
    {
      size = s;
      capacity = c;
      stack_ = new DATA_TYPE[capacity];
      std::fill(stack_, stack_+capacity, DATA_TYPE());
    }

    void initialise(const unsigned int s, const unsigned int c, DATA_TYPE x)
    {
      size = s;
      capacity = c;
      stack_ = new DATA_TYPE[capacity];
      std::fill(stack_, stack_+capacity, x);
    }

    void copy(Vector<DATA_TYPE>& vec) {
      size = vec.size;
      capacity = vec.capacity;
      stack_ = vec.stack_;
    }

    void neutralise() {
      stack_ = NULL;
    }

    void extendStack( const unsigned int l=0 )
    {

      // if(size > 1400) {
      // std::cerr << "extend " << size << std::endl;
      // //display(std::cout);
      // std::cout << std::endl;

      unsigned int increment = (l ? l : (capacity+1) << 1);
      capacity += increment;

      DATA_TYPE* new_stack = new DATA_TYPE[capacity];
      for(unsigned int i=0; i<capacity-increment; ++i)
       	new_stack[i] = stack_[i];

      delete [] stack_;
      stack_ = new_stack;
      
      std::fill(stack_+capacity-increment, stack_+capacity, DATA_TYPE());


      // std::cout << "extended " ;
      // //display(std::cout);
      // std::cout << std::endl << std::endl;


    }

    void resize( const unsigned int l )
    {
      if( capacity < l ) {
	DATA_TYPE* new_stack_ = new DATA_TYPE[l];
	memcpy(new_stack_, stack_, capacity);
	delete [] stack_;
	stack_ = new_stack_;
	std::fill(stack_+capacity, stack_+l, (DATA_TYPE)0);

	capacity = l;
      }
      size = l;
    }
    //@}

    /*!@name Accessors*/
    //@{

    void sort() {
      qsort(stack_, size, sizeof(DATA_TYPE), increasing_order<DATA_TYPE>);
    }

    inline iterator begin() const {
      return stack_;
    }
    inline iterator end() const {
      return stack_+size;
    }

    inline int empty() const
    {
      return !size;
    }
  
    inline void add(DATA_TYPE x)
    {
      if( capacity == size ) 
	extendStack();
      stack_[size++] = x;
    }

    inline void secure(const DATA_TYPE x) 
    {
      if( capacity == size ) 
	extendStack();
    }

    inline void fast_add(DATA_TYPE x)
    {
      stack_[size++] = x;
    }

    inline void push_back(DATA_TYPE x)
    {
      if( capacity == size ) 
	extendStack();
      stack_[size++] = x;
    }

    inline DATA_TYPE pop_until(const unsigned int level)
    {
      size = level;
      return stack_[size];
    }

    inline DATA_TYPE pop()
    {
      return stack_[--size];
    }

    inline void pop(DATA_TYPE& x)
    {
      x = stack_[--size];
    }

    inline void clear()
    {
      size = 0;
    }

    inline void remove(const unsigned int i)
    {  
      stack_[i] = stack_[--size];
    }

    inline void remove_elt(DATA_TYPE& elt)
    {
      unsigned int j=size;
      while(j && stack_[--j] != elt);
      stack_[j] = stack_[--size];
    }

    inline void setBack(const DATA_TYPE& x, const int k=1)
    {
      stack_[size-k] = x;
    }

    inline DATA_TYPE& front(const int k=0)
    {
      return stack_[k];
    }

    inline const DATA_TYPE front(const int k=0) const
    {
      return stack_[k];
    }

    inline DATA_TYPE& back(const int k=1)
    {
      return stack_[size-k];
    }

    inline const DATA_TYPE back(const int k=1) const
    {
      return stack_[size-k];
    }

    inline DATA_TYPE& operator[](const int i)
    {
      return stack_[i];
    }

    inline const DATA_TYPE operator[](const int i) const
    {
      return stack_[i];
    }

    inline Vector< DATA_TYPE >& operator=(const Vector< DATA_TYPE >& x)
    {
      if(x.capacity && !stack_) {
	initialise(0, x.capacity);
      } else if(capacity<x.size) {
	extendStack(x.capacity-capacity);
      }

      clear();

      for(unsigned int i=0; i<x.size; ++i)
	add(x[i]);
      
      return *this;
    }
    //@}

    /*!@name Printing*/
    //@{
    std::ostream& display(std::ostream& os) const {
      os << "[";
      if(size) os << stack_[0] ;
      for(unsigned int i=1; i<size; ++i)
	os << " " << stack_[i];
      os << "]";
      return os;
    }
    //@}

  };


template < int N, class T >
  
  class Tuple {
  
 public:
  
  T data[N];
  
  Tuple() {}
  Tuple(T a) {data[0] = a;}
  Tuple(T a, T b) {data[0] = a; data[1] = b;}
  Tuple(T a, T b, T c) {data[0] = a; data[1] = b; data[2] = c;}
  Tuple(T a, T b, T c, T d) {data[0] = a; data[1] = b; data[2] = c; data[3] = d;}
  Tuple(T a, T b, T c, T d, T e) {data[0] = a; data[1] = b; data[2] = c; data[3] = d; data[4] = e;}

  T operator[](const int i) const { return data[i]; }
  T& operator[](const int i) { return data[i]; }
  
  /*!@name Printing*/
  //@{
  std::ostream& display(std::ostream& os) const {
    os << "<";
    if(N) os << data[0] ;
    for(unsigned int i=1; i<N; ++i)
      os << " " << data[i];
    os << ">";
    return os;
  }
  //@}

};


  typedef Tuple< 2, int > Pair;
  
  template <class T1, class T2, class T3>
  class Triplet {
  public:
    T1 first;
    T2 second;
    T3 third;

    Triplet() {}
    Triplet(T1 t1, T2 t2, T3 t3) : first(t1), second(t2), third(t3) {}
    ~Triplet() {}

    operator T1() { return first; }

   std::ostream& display(std::ostream& os) const {
      os << "<" << first << ", " << second << ", " << third << ">";
      return os;
    }    

  };


  template <class T>
  class TwoWayStack {
  public:

    /*!@name Parameters*/
    //@{
    /** <stack_> is an array of object T such that:
	- The first element of the stack is the one at index <start>
	- There are <size> elements;
	- The current size of the stack is <capacity>
	Moreover, we assume that elements T can be casted into int, and that
	this integer is between 0 and <capacity>, then we have:
	- <index_[i]> is the index_ of element i in <stack_>
     */

    T   *stack_;
    int *index_;

    // index
    unsigned int start;
    unsigned int capacity;
    unsigned int size;
    //@}

    /*!@name Constructor*/
    //@{
    TwoWayStack()
    {
      start = 0;
      capacity = 0;
      size = 0;
      stack_ = NULL;
      index_ = NULL;
    }
    //@}

    /*!@name Destructor*/
    //@{  
    virtual ~TwoWayStack()
    {
      free( stack_ );
      free( index_ );
    }
    //@}

    /*!@name Initialisation*/
    //@{
    void initialise(const unsigned int c)
    {
      start = 0;
      size = 0;
      capacity = c;

      stack_ = (T  *) malloc(capacity*sizeof(T  ));
      index_ = (int*) malloc(capacity*sizeof(int));
      
      int f = sizeof(T)/sizeof(int);

      //int out = (start+size+1)%capacity;

      std::fill((int*)stack_, (int*)stack_+(capacity*f), 0);
      std::fill(index_, index_+capacity, -1);



      // std::cout << "CREATE " << std::endl;
      // for(int i=0; i<capacity; ++i)
      //  	std::cout << index_[i] << " ";
      // std::cout << std::endl;
    }


    // void check_intergrity() {
    //   // check that elements 
    // }


    void extend_stack( const unsigned int l=0 )
    {
      unsigned int increment = (l ? l : (capacity+1) << 1);
      int f = sizeof(T)/sizeof(int);
      //capacity += increment;

      int* new_index = (int*) malloc((capacity+increment)*sizeof(int));
      for(unsigned int i=0; i<capacity; ++i) {
	new_index[i] = (index_[i]>=0 ? (index_[i]+capacity-start)%capacity : -1);
      }
      // for(int i=start; i<capacity; ++i) {
      // 	new_index[i] = index_[i]-start;
      // }
      //int out = (start+size+1)%(capacity+increment);
      std::fill(new_index+capacity, new_index+capacity+increment, -1);

      T* new_stack = (T*) malloc((capacity+increment)*sizeof(T));
      memcpy(new_stack, stack_+start, (capacity-start)*sizeof(T));
      memcpy(new_stack+capacity-start, stack_, (start)*sizeof(T));
      std::fill((int*)new_stack+(capacity*f), (int*)new_stack+(capacity+increment)*f, 0);

      free(stack_); 
      stack_ = new_stack;

      free(index_); 
      index_ = new_index;

      capacity += increment;
      start = 0;

      // std::cout << "EXTEND" << std::endl;
      // for(int i=0; i<capacity; ++i)
      //  	std::cout << index_[i] << " ";
      // std::cout << std::endl;

    }

    void declare( const int x ) {
      if(x >= (int)capacity) extend_stack();
    }
    //@}

    /*!@name Accessors*/
    //@{
    inline int empty() const
    {
      return !size;
    }

    inline void probe() { if(size == capacity) extend_stack(); }
  
    inline void push_back(T x)
    {

      //if((int)x == 355) std::cout << "PUSH B x355" << std::endl; 

      int idx = (start+(size++))%capacity;
      stack_[idx] = x;
      index_[(int)x] = idx;
    }

    inline void push_front(T x)
    {

      //if((int)x == 355) std::cout << "PUSH F x355" << std::endl; 

      ++size;
      start = ((start+capacity-1)%capacity);
      stack_[start] = x;

      index_[(int)x] = start;
    }

    inline bool contain(const int x) {

      // if(x == 355) {

      // 	std::cout << "contain 355? " << index_[x] << std::endl;

      // }


      return index_[x] >= 0;

      //return (int)(stack_[index_[x]]) == x;

      // return (( capacity+
      // 		index_[x]-start)%capacity) < size;
    }

    inline T& operator[](const int x) {
      return stack_[index_[x]];
    }

    inline T pop_back()
    {
      int idx = (start+--size)%capacity;
      index_[(int)stack_[idx]] = -1;
      return stack_[idx];
    }

    inline T pop_front()
    {
      --size;
      T rval = stack_[start];
      index_[(int)stack_[start]] = -1;
      start = (start+1)%capacity;
      return rval;
    }

    inline void pop_back(T& x)
    {
      int idx = (start+--size)%capacity;
      index_[(int)stack_[idx]] = -1;
      //index_[idx] = -1;
      x = stack_[idx];
    }

    inline void pop_front(T& x)
    {
      --size;
      x = stack_[start];
      index_[(int)stack_[start]] = -1;
      //index_[start] = -1;
      start = (start+1)%capacity;
    }

    inline void clear()
    {
      for(unsigned int i=0; i<size; ++i) {
       	index_[stack_[(start+i)%capacity]] = -1;
      }
      size = 0;
    }

    //@}

    /*!@name Printing*/
    //@{
    std::ostream& display(std::ostream& os) const {
      os << "[";
      if(size) os << stack_[start]  ;
      for(unsigned int i=1; i<size; ++i)
	os << " " << stack_[(start+i)%capacity] ;
      os << "]";
      return os;
    }
    //@}

  };




  template <class MAIN_TYPE, class AUX_TYPE>
  class BiStack {
  public:

    /*!@name Parameters*/
    //@{
    MAIN_TYPE* main_stack_;
    AUX_TYPE* aux_stack_;
    unsigned int start;
    unsigned int capacity;
    unsigned int size;
    //@}

    /*!@name Constructor*/
    //@{
    BiStack()
    {
      start = 0;
      capacity = 0;
      size = 0;
      main_stack_ = NULL;
      aux_stack_ = NULL;
    }
    //@}

    /*!@name Destructor*/
    //@{  
    virtual ~BiStack()
    {
      free( main_stack_ );
      free( aux_stack_ );
    }
    //@}

    /*!@name Initialisation*/
    //@{
    void initialise(const unsigned int c)
    {
      start = 0;
      size = 0;
      capacity = c;

      main_stack_ = (MAIN_TYPE*) malloc(capacity*sizeof(MAIN_TYPE));
      aux_stack_ = (AUX_TYPE*) malloc(capacity*sizeof(AUX_TYPE));
      
      int f = sizeof(MAIN_TYPE)/sizeof(int);
      std::fill((int*)main_stack_, (int*)main_stack_+(capacity*f), 0);
      f = sizeof(AUX_TYPE)/sizeof(int);
      std::fill((int*)aux_stack_, (int*)aux_stack_+(capacity*f), 0);
    }


    void extend_stack( const unsigned int l=0 )
    {
      unsigned int increment = (l ? l : (capacity+1) << 1);
      capacity += increment;

      MAIN_TYPE* new_main_stack = (MAIN_TYPE*) malloc(capacity*sizeof(MAIN_TYPE));
      memcpy(new_main_stack, main_stack_, (capacity-increment)*sizeof(MAIN_TYPE));

      AUX_TYPE* new_aux_stack = (AUX_TYPE*) malloc(capacity*sizeof(AUX_TYPE));
      memcpy(new_aux_stack, aux_stack_, (capacity-increment)*sizeof(AUX_TYPE));

      free(main_stack_); 
      main_stack_ = new_main_stack;

      free(aux_stack_); 
      aux_stack_ = new_aux_stack;

      int f = sizeof(MAIN_TYPE)/sizeof(int);
      std::fill((int*)main_stack_, (int*)main_stack_+(capacity*f), 0);
      f = sizeof(AUX_TYPE)/sizeof(int);
      std::fill((int*)aux_stack_, (int*)aux_stack_+(capacity*f), 0);
    }
    //@}

    /*!@name Accessors*/
    //@{
    inline int empty() const
    {
      return !size;
    }

    inline void probe() { if(size == capacity) extend_stack(); }
  
    inline void push_back(MAIN_TYPE x, AUX_TYPE y)
    {
      int idx = (start+(size++))%capacity;
      main_stack_[idx] = x;
      aux_stack_[idx] = y;
    }

    inline void push_front(MAIN_TYPE x, AUX_TYPE y)
    {
      ++size;
      start = (start+capacity-1)%capacity;
      main_stack_[start] = x;
      aux_stack_[start] = y;
    }


    inline void pop_back(MAIN_TYPE& x, AUX_TYPE& y)
    {
      int idx = (start+--size)%capacity;
      x = main_stack_[idx];
      y = aux_stack_[idx];
    }

    inline void pop_front(MAIN_TYPE& x, AUX_TYPE& y)
    {
      --size;
      x = main_stack_[start];
      y = aux_stack_[start];
      start = (start+1)%capacity;
    }

    inline void clear()
    {
      size = 0;
    }

    //@}

    /*!@name Printing*/
    //@{
    std::ostream& display(std::ostream& os) const {
      for(unsigned int i=0; i<capacity; ++i)
	os << main_stack_[i] << " ";
      os << std::endl;

      os << "[";
      if(size) os << main_stack_[start] << "/" << aux_stack_[start] ;
      for(unsigned int i=1; i<size; ++i)
	os << " " << main_stack_[(start+i)%capacity] << "/" << aux_stack_[(start+i)%capacity];
      os << "]";
      return os;
    }
    //@}

  };





  /**********************************************
   * Vector2
   **********************************************/
  /*! \class Vector2
    \brief Simple vector class     
  */
  
  template <class DATA_TYPE>
  class Vector2 {
  public:

    /*!@name Parameters*/
    //@{
    DATA_TYPE* stack_;
    unsigned int capacity;
    unsigned int size;
    //@}

    /*!@name Constructor*/
    //@{
    Vector2()
    {
      capacity = 0;
      size = 0;
      stack_ = NULL;
    }

    Vector2(const Vector2<DATA_TYPE>& s)
    {
      capacity = s.size;
      stack_ = new DATA_TYPE[capacity];
      //(DATA_TYPE*) malloc(capacity*sizeof(DATA_TYPE));
      for(size = 0; size<capacity; ++size)
  	stack_[size] = s[size];
    }

    Vector2(const int n)
    {
      capacity = n;
      size = n;
      if( capacity )
  	//stack_ = (DATA_TYPE*) malloc(capacity*sizeof(DATA_TYPE));
  	stack_ = new DATA_TYPE[capacity]; //(DATA_TYPE*) malloc(capacity*sizeof(DATA_TYPE));
      else stack_ = NULL;
    }
    //@}

    /*!@name Destructor*/
    //@{  
    virtual ~Vector2()
    {

#ifdef _DEBUG_MEMORY
      std::cout << "c delete vector: " << size << " " << capacity << " " // ;
      // display(std::cout);
      // std::cout 
	<< std::endl;
#endif

      delete [] stack_;
      //free( stack_ );
    }
    //@}

    /*!@name Initialisation*/
    //@{
    void initialise(const unsigned int c)
    {
      size = 0;
      capacity = c;
      //stack_ = (DATA_TYPE*) malloc(capacity*sizeof(DATA_TYPE));
      stack_ = new DATA_TYPE[capacity];

      DATA_TYPE x(0);
      std::fill(stack_, stack_+capacity, x);
    }

    void initialise(const unsigned int s, const unsigned int c)
    {
      size = s;
      capacity = c;
      //stack_ = (DATA_TYPE*) malloc(capacity*sizeof(DATA_TYPE));
      stack_ = new DATA_TYPE[capacity];

      DATA_TYPE x(0);
      std::fill(stack_, stack_+capacity, x);
    }

    //     void extendStack()
    //     {
    //       unsigned int new_capacity = ((capacity+1) << 1);
    //       DATA_TYPE* new_stack = (DATA_TYPE*) realloc(stack_, new_capacity*sizeof(DATA_TYPE));
    //       stack_ = new_stack;

    //       DATA_TYPE x;
    //       std::fill(stack_+capacity, stack_+new_capacity, x);
      
    //       capacity = new_capacity;
    //     }

    void extendStack( const unsigned int l=0 )
    {
      unsigned int increment = (l ? l : (capacity+1) << 1);
      capacity += increment;
      //DATA_TYPE* new_stack = (DATA_TYPE*) realloc(stack_, capacity*sizeof(DATA_TYPE));
      DATA_TYPE* new_stack_ = new DATA_TYPE[capacity];
      memcpy(new_stack_, stack_, capacity-increment);
      delete [] stack_;
      stack_ = new_stack_;

      DATA_TYPE x(0);
      std::fill(stack_+capacity-increment, stack_+capacity, x);
    }

    void resize( const unsigned int l )
    {
      if( capacity < l ) {
	int old_capacity = capacity;
  	capacity = l;
  	// DATA_TYPE* new_stack = (DATA_TYPE*) realloc(stack_, capacity*sizeof(DATA_TYPE));
  	// stack_ = new_stack;

	DATA_TYPE* new_stack_ = new DATA_TYPE[capacity];
	memcpy(new_stack_, stack_, old_capacity);
	delete [] stack_;
	stack_ = new_stack_;
	
	DATA_TYPE x(0);
	std::fill(stack_+old_capacity, stack_+capacity, x);
      }
      size = l;
    }
    //@}

    /*!@name Accessors*/
    //@{
    inline int empty() const
    {
      return !size;
    }
  
    inline void add(DATA_TYPE x)
    {
      if( capacity == size ) 
  	extendStack();
      stack_[size++] = x;
    }

    inline void secure(const DATA_TYPE x) 
    {
      if( capacity == size ) 
  	extendStack();
    }

    inline void fast_add(DATA_TYPE x)
    {
      stack_[size++] = x;
    }

    inline void push_back(DATA_TYPE x)
    {
      if( capacity == size ) 
  	extendStack();
      stack_[size++] = x;
    }

    inline DATA_TYPE pop_until(const unsigned int level)
    {
      size = level;
      return stack_[size];
    }

    inline DATA_TYPE pop()
    {
      return stack_[--size];
    }

    inline void pop(DATA_TYPE& x)
    {
      x = stack_[--size];
    }

    inline void clear()
    {
      size = 0;
    }

    inline void remove(const unsigned int i)
    {  
      stack_[i] = stack_[--size];
    }

    inline void remove_elt(DATA_TYPE& elt)
    {
      unsigned int j=size;
      while(j && stack_[--j] != elt);
      stack_[j] = stack_[--size];
    }

    inline void setBack(const DATA_TYPE& x, const int k=1)
    {
      stack_[size-k] = x;
    }

    inline DATA_TYPE& front(const int k=0)
    {
      return stack_[k];
    }

    inline const DATA_TYPE front(const int k=0) const
    {
      return stack_[k];
    }

    inline DATA_TYPE& back(const int k=1)
    {
      return stack_[size-k];
    }

    inline const DATA_TYPE back(const int k=1) const
    {
      return stack_[size-k];
    }

    inline DATA_TYPE& operator[](const int i)
    {
      return stack_[i];
    }

    inline const DATA_TYPE operator[](const int i) const
    {
      return stack_[i];
    }

    inline Vector2< DATA_TYPE >& operator=(const Vector2< DATA_TYPE >& x)
    {
      initialise(0, x.capacity);
      for(unsigned int i=0; i<x.size; ++i)
  	add(x[i]);
      return *this;
    }
    //@}

    /*!@name Printing*/
    //@{
    std::ostream& display(std::ostream& os) const {
      os << "[";
      if(size) os << stack_[0] ;
      for(unsigned int i=1; i<size; ++i)
  	os << " " << stack_[i];
      os << "]";
      return os;
    }
    //@}

  };


  /**********************************************
   * Array
   **********************************************/    
  /*! \class Array
    \brief Simple array class 
  */

  typedef unsigned int Atom;
  typedef unsigned int Literal;

  //class Decision;
  class Explanation {

  public:
    
    typedef Literal* iterator;

 
    //typedef Decision* iterator;

    virtual iterator get_reason_for(const Atom a, const int lvl, iterator& end) = 0;

    virtual std::ostream& display(std::ostream& os) const = 0;

    virtual bool is_clause() {return true;}
    
  };


  
  template <class DATA_TYPE>
  class Array : public Explanation {

  public: 

    /*!@name Parameters*/
    //@{
    unsigned int size;
    DATA_TYPE data[0];
    //@}

    Array(const Vector<DATA_TYPE>& ps) 
    {
      size = ps.size;
      for (unsigned int i=0; i<ps.size; ++i) 
	data[i] = ps[i];
    }

    virtual ~Array() {}

    // virtual Explanation::iterator begin(Atom a) { return &(data[0]); }
    // virtual Explanation::iterator end  (Atom a) { return &(data[size]); }
    
    virtual iterator get_reason_for(const Atom a, const int lvl, iterator& end) { end = &(data[size]); return &(data[0]); }

    static Array<DATA_TYPE>* Array_new(const Vector<DATA_TYPE>& ps)
    {
      void* mem = malloc(sizeof(Array<DATA_TYPE>) + sizeof(DATA_TYPE)*(ps.size));
      return new (mem) Array<DATA_TYPE>(ps); 
    }

    inline DATA_TYPE& operator [] (const int i) { return data[i]; }
    inline DATA_TYPE operator [] (const int i) const { return data[i]; }

    std::ostream& display(std::ostream& os) const {
      os << "[";
      if(size) os << data[0];
      for(unsigned int i=1; i<size; ++i)
	os << " " << data[i];
      os << "]";
      return os;
    }

  };


  /**********************************************
   * MultiSet
   **********************************************/
  /// Sparse multi-set representation

  class MultiSet 
  {
  public:

    /*!@name Parameters*/
    //@{
    /// list of values
    int *values_; // contiguous list of values (repeated)
    unsigned int **index_; // index[v][o] is the index of occurrence 'o' of value 'v' in 'values_'
    unsigned int *occurrence_; // self explanatory
    
    /// current max capacity
    unsigned int capacity;
    /// current size (total number of occurrences)
    unsigned int size;

    /// minimum and maximum value 
    int lb_;
    int ub_;
    //@}

    /*!@name Constructors*/
    //@{
    MultiSet()
    {
      size = 0;
      capacity = 0;
      lb_ = 0;
      ub_ = -1;

      values_ = NULL;
      index_ = NULL;
      occurrence_ = NULL;
    }

    MultiSet(const int lb, const int ub, const int* occ=NULL, const int* mocc=NULL)
    {
      if(!occ) {
	int tmp_occ[ub-lb+1];
	int tmp_mocc[ub-lb+1];
	for(int i=0; i<ub-lb+1; ++i) {
	  tmp_occ[i] = 0;
	  tmp_mocc[i] = 4;
	}
	initialise(lb, ub, tmp_occ, tmp_mocc);
      } else if(!mocc) {
	int tmp_mocc[ub-lb+1];
	for(int i=0; i<ub-lb+1; ++i) 
	  tmp_mocc[i] = occ[i];
	initialise(lb, ub, occ, tmp_mocc);
      } else {
	initialise(lb, ub, occ, mocc);
      }
    }

    virtual ~MultiSet()
    {
      delete [] values_;
      occurrence_ += lb_;
      delete [] occurrence_;
      for(int v=lb_; v<=ub_; ++v)
	delete [] index_[v];
      index_ += lb_;
      delete [] index_;
    }

    void initialise(const int lb, const int ub, const int* occ, const int* mocc)
    {
      int span = ub-lb+1, i=0, j=0, v, o;

      lb_ = lb;
      ub_ = ub;

      index_ = new unsigned int*[span];
      occurrence_ = new unsigned int[span];
      capacity = 0;
      size = 0;

      for(v=0; v<span; ++v) {
	capacity += mocc[v];
	size += occ[v];
	occurrence_[v] = occ[v];
	index_[v] = new unsigned int[mocc[v]];
      }
      
      occurrence_ -= lb_;
      index_ -= lb_;
      values_ = new int[capacity];
      
      for(v=lb_; v<=ub_; ++v) {
	for(o=0; o<occ[v-lb_]; ++o) {
	  values_[i] = v;
	  index_[v][o] = i++;
	}
	for(o=occ[v-lb_]; o<mocc[v-lb_]; ++o) {
	  values_[size+j] = v;
	  index_[v][o] = size+j++;
	}
      }
    }
    //@}    


    inline void remove(const int elt)
    {
      --size;
      --occurrence_[elt];

      int v = values_[size];
      int i = index_[elt][occurrence_[elt]];

      index_[v][occurrence_[v]-1] = i;
      index_[elt][occurrence_[elt]] = size;

      values_[i] = values_[size];
      values_[size] = elt;
    }

    inline int pop()
    {
      int v = values_[--size];
      --occurrence_[v];
      return v;
    }

    inline int pop_head()
    {
      const int elt = values_[0];

      --size;
      --occurrence_[elt];

      index_[values_[size]][occurrence_[values_[size]]-1] = 0;
      index_[elt][occurrence_[elt]] = size;

      values_[0] = values_[size];
      values_[size] = elt;

      return elt;
    }

    inline int head()
    {
      return *values_;
    }
    
    inline int back()
    {
      return values_[size-1];
    }

    inline void add(const int elt)
    {
      int i = index_[elt][occurrence_[elt]];
      int v = values_[size]; 

      index_[v][occurrence_[v]] = i;
      index_[elt][occurrence_[elt]] = size;

      values_[i] = v;
      values_[size] = elt;

      ++occurrence_[elt];
      ++size;
    }
    //@}

    /*!@name Miscellaneous*/
    //@{
    //     std::string getString() const {
    //       std::string return_str = "{{";
    //       bool first = true;
    //       for(int v=lb_; v<=ub_; ++v) {
    // 	if(occurrence_[v] > 0) {
    // 	  if(first) first = false;
    // 	  else return_str += ",";
    // 	  return_str += ("("+toString((int)(occurrence_[v]))+"*"+toString(v)+")");
    // 	}
    //       }
    //       return_str += "}}";
    //       return return_str;
    //     }

    std::ostream& display(std::ostream& os) const {
      os << "{{";
      bool first = true;
      for(int v=lb_; v<=ub_; ++v) {
	if(occurrence_[v] > 0) {
	  if(first) first = false;
	  else os << ",";
	  os << "(" << occurrence_[v] << "*" << v << ")" ;
	}
      }
      os << "}}";
      return os;
    }
    //@}
  };




  /**********************************************
   * IntStack
   **********************************************/
  /// Sparse set representation

  class IntStack 
  {
  public:

    //typedef int* iterator;

    /*!@name Parameters*/
    //@{
    /// list of values
    int *list_;
    /// current max capacity
    unsigned int index_capacity;
    unsigned int list_capacity;
    /// current size
    unsigned int size;
    /// values' indices
    unsigned int *index_;
    unsigned int *start_;
    //@}

    /*!@name Constructors*/
    //@{
    IntStack()
    {
      size = 0;
      index_capacity = 0;
      list_capacity = 0;

      list_ = NULL;
      index_ = NULL;
      start_ = NULL;
    }

    IntStack(const int lb, const int ub, const bool full=true)
    {
      initialise(lb, ub, ub-lb+1, full);
    }

    IntStack(const int lb, const int ub, const int sz, const bool full)
    {
      initialise(lb, ub, sz, full);
    }

    IntStack(IntStack& shared, const int sz)
    {
      initialise(shared, sz);
    }

    void  neutralise() {
      start_ = NULL;
    }
    
    virtual ~IntStack()
    {
      delete [] list_;
      delete [] start_;
    }

    virtual void initialise(IntStack& shared, const int sz);
    virtual void initialise(const int lb, const int ub, const int sz, const bool full);
    // {

    //   std::cout << sz << std::endl;

    //   capacity = ub-lb+1;
    //   list_ = new int[sz];
    //   start_ = new unsigned int[capacity];
    //   index_ = start_ - lb;      

    //   for(int i=lb; i<=ub; ++i) 
    // 	{
    // 	  index_[i] = i-lb;
    // 	  list_[i-lb] = i;
    // 	}
      
    //   size = (full ? sz : 0);
    // }




    // void extend_list()
    // {
    //   unsigned int increment = ((list_capacity+1) << 1);
    //   list_capacity += increment;

    //   int* new_list = new int[list_capacity];
    //   memcpy(new_list, list_, (list_capacity-increment)*sizeof(int));

    //   // for(unsigned int i=0; i<list_capacity-increment; ++i)
    //   //  	new_list[i] = list_[i];

    //   delete [] list_;
    //   list_ = new_list;
    // }

    // void extend(const int new_elt)
    // {
    //   int lb = (int)(start_-index_), new_lb = lb;
    //   int ub = index_capacity+lb-1, new_ub = ub;
    //   if(new_elt < lb) {
    // 	new_lb = new_elt;
    //   } else if(new_elt > ub) {
    // 	new_ub = new_elt;
    //   } else {
    // 	return;
    //   }
      
    //   unsigned int new_index_capacity = new_ub-new_lb+1;
    //   if(new_index_capacity < index_capacity*2) new_index_capacity = index_capacity*2;
    //   if(new_lb < lb) {
    // 	new_lb = ub-new_index_capacity+1;
    //   } else {
    // 	new_ub = lb+new_index_capacity-1;
    //   }



    //   unsigned int *aux_start = start_;
    //   start_ = new unsigned int[new_index_capacity];
    //   memcpy(start_+(lb-new_lb), aux_start, index_capacity*sizeof(unsigned int));
    //   delete [] aux_start;

    //   index_ = start_ - new_lb;
    //   int k = 0;
    //   for(int i=new_lb; i<lb; ++i) {
    // 	index_[i] = size+k;
    // 	//list_[index_capacity+k++] = i;
    //   }
    //   for(int i=ub+1; i<=new_ub; ++i) {
    // 	index_[i] = size+k;
    // 	//list_[index_capacity+k++] = i;
    //   }

    //   index_capacity = new_index_capacity;
    // }
    // //@}    

    // /*!@name Accessors*/
    // //@{ 
    // inline int get_min() const 
    // {
    //   int the_min = INFTY;
    //   int offset = (index_ - start_);
    //   unsigned int explored;
    //   int stop_crit = size+offset;

    //   for(explored = 0; explored < size 
    // 	    //&& size - explored <= (the_min - offset); 
    // 	    && the_min + (int)explored >= stop_crit; 
    // 	  ++explored) {
    // 	if(list_[explored] < the_min) the_min = list_[explored];
    //   }
    //   if(explored < size) {
    // 	while(offset < the_min && index_[offset]>=size) ++offset;
    // 	the_min = offset;
    //   }



    //   // int the_min = INFTY;
    //   // if(size) {
    //   // 	the_min = list_[0];
    //   // 	// ratio size / index_capacity
    //   // 	if(4*size < index_capacity) {
    //   // 	  for(unsigned int i=1; i<size; ++i)
    //   // 	    if(list_[i] < the_min) the_min = list_[i];
    //   // 	} else {
    //   // 	  int val=(index_ - start_);
    //   // 	  while( val<the_min && index_[val] >= size )
    //   // 	    ++val;
    //   // 	  the_min = val;
    //   // 	}
    //   // }


    //   return the_min;
    // }

    // inline int get_max() const 
    // {
    //   int the_max = -INFTY;
    //   if(size) {
    // 	the_max = list_[0];
    // 	// ratio size / index_capacity
    // 	if(4*size < index_capacity) {
    // 	  for(unsigned int i=1; i<size; ++i)
    // 	    if(list_[i] > the_max) the_max = list_[i];
    // 	} else {
    // 	  int val=(index_capacity + index_ - start_);
    // 	  while( val>the_max && index_[val] >= size )
    // 	    --val;
    // 	  the_max = val;
    // 	}
    //   }
    //   return the_max;
    // }


    // inline bool safe_contain(const int elt) const 
    // {
    //   return ((unsigned)(elt-(int)(start_-index_))<index_capacity && index_[elt]<size);
    // } 
 
    // inline bool contain(const int elt) const 
    // {
    //   return index_[elt]<size;
    // } 
  
    // inline bool empty()const 
    // {
    //   return !size;
    // } 

    // inline int next(const int elt) const
    // {
    //   unsigned int idx = index_[elt]+1;
    //   return (idx < size ? list_[idx] : elt);
    // }
    // inline int prev(const int elt) const
    // {
    //   unsigned int idx = index_[elt]-1;
    //   return (idx >= 0 ? list_[idx] : elt);
    // }
    
    // inline int operator[](const unsigned int idx) const
    // {
    //   return list_[idx];
    // }

    // inline int& operator[](const unsigned int idx)
    // {
    //   return list_[idx];
    // }
    // //@}

    // /*!@name List Manipulation*/
    // //@{
    // inline int* begin() 
    // {
    //   return list_;
    // }

    // inline int* end() 
    // {
    //   return &(list_[size]);
    // }

    // inline int* end_mem() 
    // {
    //   return list_+list_capacity;
    // }

    // inline void fill()
    // {
    //   size = list_capacity;
    // }

    // inline void clear()
    // {
    //   size = 0;
    // }
  
    // inline void set_to(const int elt)
    // {
    //   size=1;
    //   index_[*list_] = index_[elt];
    //   list_[index_[elt]] = *list_;
    //   *list_ = elt;
    //   index_[elt] = 0;
    // }

    // inline void remove(const int elt)
    // {
    //   --size;
    //   index_[list_[size]] = index_[elt];
    //   list_[index_[elt]] = list_[size];
    //   list_[size] = elt;
    //   index_[elt] = size;
    // }

    // inline int next()
    // {
    //   return list_[size];
    // }

    // inline int pop()
    // {
    //   return list_[--size];
    // }

    // inline int pop_head()
    // {
    //   --size;
    //   index_[list_[size]] = 0;
    //   const int elt = *list_;
    //   *list_ = list_[size];
    //   list_[size] = elt;
    //   index_[elt] = size;
    //   return elt;
    // }

    // inline int head() const
    // {
    //   return *list_;
    // }
    
    // inline int back() const
    // {
    //   return list_[size-1];
    // }

    // inline void init_add(const int elt)
    // {
    //   if(size == list_capacity) extend_list();
    //   if((unsigned)(elt-(int)(start_-index_)) >= index_capacity) {
    // 	extend(elt);
    //   }
    //   list_[size] = elt;
    //   index_[elt] = size;
    //   ++size;
    // }

    // inline void add(const int elt)
    // {
    //   index_[list_[size]] = index_[elt];
    //   list_[index_[elt]] = list_[size];
    //   list_[size] = elt;
    //   index_[elt] = size;
    //   ++size;
    // }

    // inline void safe_add(const int elt)
    // {
    //   if(!safe_contain(elt)) extend(elt);
    //   add(elt);
    // }

    // // create a new element that can potentially be outside the bounds
    // inline void create(const int elt)
    // {
    //   extend(elt);
    //   add(elt);
    // }

    // inline void ordered_add(const int elt)
    // {
    //   // the first non-element goes where elt was
    //   index_[list_[size]] = index_[elt];
    //   list_[index_[elt]] = list_[size];

    //   int idx = size;
    //   while(idx && list_[idx-1] > elt) { // push every values greater than elt above elt
    // 	list_[idx] = list_[idx-1];
    // 	index_[list_[idx-1]] = idx;
    // 	--idx;
    //   }

    //   list_[idx] = elt;
    //   index_[elt] = idx;
    //   ++size;
    // }


    // inline void revert_to(const int level)
    // {
    //   size = level;
    // }

    // inline void index()
    // {
    //   for(unsigned int i=0; i<list_capacity; ++i)
    // 	index_[list_[i]] = i;
    // }
    // //@}



    void extend_list();
    void extend(const int new_elt);


    /*!@name Accessors*/
    //@{ 
     int get_min() const;

     int get_max() const;

     bool safe_contain(const int elt) const ;

     bool contain(const int elt) const ;
	 
	 int get_index(const int elt) const ;
  
     bool empty()const ;

     int next(const int elt) const ;

     int prev(const int elt) const ;
    
     int operator[](const unsigned int idx) const ;

     int& operator[](const unsigned int idx) ;
    //@}

    /*!@name List Manipulation*/
    //@{
     int* begin() ;

     int* end() ;

     int* end_mem() ;

     void fill() ;

     void clear() ;
  
     void set_to(const int elt) ;

     void remove(const int elt) ;

     int next() ;

     int pop() ;

     int pop_head() ;

     int head() const ;
    
     int back() const ;

     void init_add(const int elt) ;

     void add(const int elt) ;
	 
	 void set(const int elt, const int idx) ;

     void safe_add(const int elt) ;

    // create a new element that can potentially be outside the bounds
     void create(const int elt) ;

     void ordered_add(const int elt) ;

     void revert_to(const int level) ;

     void index() ;
    //@}


    /*!@name Miscellaneous*/
    //@{
    std::ostream& display(std::ostream& os) const;
	std::string to_str() const;
    // std::ostream& display(std::ostream& os) const {
    //   int min =  INFTY;
    //   int max = -INFTY;

    //   for(unsigned int i=0; i<size; ++i) {
    // 	if(list_[i] < min) min = list_[i];
    // 	if(list_[i] > max) max = list_[i];
    //   }

    //   Bitset< unsigned long int > elts(min, max, BitSet::empt);
    //   for(unsigned int i=0; i<size; ++i) {
    // 	elts.add(list_[i]);
    //   }

    //   os << elts;
    //   // os << "(";
    //   // bool not_empty = (size>0);

    //   // if(not_empty) os << list_[0];
    //   // for(unsigned int i=1; i<size; ++i)
    //   // 	os << " " << list_[i];
    //   // os << ")";
    //   return os;
    // }
    //@}
  };
  
  
  
  class Graph {
	
  public:

  	int       capacity;
  	IntStack      node;
  	IntStack *neighbor;
	
  	Graph() {}
	
  	Graph(const int n) {
  		initialise(n);
  	}
	
  	Graph(const Graph& g) {
				
  		initialise(g.capacity,g.size() == g.capacity);
		if(g.size() != g.capacity)
			for(int i=0; i<g.node.size; ++i) { node.add(g.node[i]); }
		for(int x=0; x<capacity; ++x) {
			for(int i=0; i<g.neighbor[x].size; ++i) {
				neighbor[x].add(g.neighbor[x][i]);
			}
		}
		
  	}
	
  	virtual void initialise(const int n, bool full=true) {
  		capacity = n;
  		node.initialise(0,capacity-1,capacity,full);
  		neighbor = new IntStack[n];
  		for(int i=0; i<n; ++i) {
  			neighbor[i].initialise(0, capacity-1, capacity, false);
  		}
  	}
	
  	virtual ~Graph() {
  		delete [] neighbor;
  	}
	
  	int size() const { return node.size; }
	
  	void clear() { 
  		while(!node.empty()) { 
  			neighbor[node.pop()].clear(); 
  		} 
  	}
	
  	void add_node(const int x) {
  		node.add(x);
  	}
	
  	void add_directed(const int x, const int y) {
  		neighbor[x].add(y);
  	}
	
  	void add_undirected(const int x, const int y) {
  		neighbor[x].add(y);
  		neighbor[y].add(x);
  	}
	
  	bool exist_arc(const int x, const int y) const {
  		return neighbor[x].contain(y);
  	}
	
  	int degree(const int x) const {
  		return neighbor[x].size;
  	}
	
  	int get_neighbor(const int x, const int i) const {
  		return neighbor[x][i];
  	}
	
  	int get_next_neighbor(const int x, const int y) const {
  		return neighbor[x].next(y);
  	}
	
  	int get_prev_neighbor(const int x, const int y) const {
  		return neighbor[x].prev(y);
  	}
	
  	std::ostream& display(std::ostream& os) const {
  		for(int i=0; i<node.size; ++i) {
  			os << node[i] << ": ";
  			neighbor[node[i]].display(os);
  			os << std::endl;
  		}
  		return os;
  	}
	
  };


  std::ostream& operator<< (std::ostream& os, const Graph& x);

  std::ostream& operator<< (std::ostream& os, const Graph* x);
  




//   //#define TO_ELEMENT(x) ((x) << 3);
//   //#define TO_INDEX(x) ((x) << 11);
//   //#define SIZE 0xfffffff8

// #define SIZE 0xfffffffc
// #define CURRENT 0x3fffffff

  
//   const int ELEMENT = {
//     /*    
// 	  0xffffffe7,
// 	  0xffffff9f,
// 	  0xfffffe7f,
// 	  0xfffff9ff
//     */

//     0xfffffff3,
//     0xffffffcf,
//     0xffffff3f,
//     0xfffffcff

//   };

//   const int INDEX = {
//     /*
//     0xffffe7ff,
//     0xffff9fff,
//     0xfffe7fff,
//     0xfff9ffff
//     */

//     0xfffff3ff,
//     0xffffcfff,
//     0xffff3fff,
//     0xfffcffff

//   };


//   const int HISTORY = {

//     0xfff3ffff,
//     0xffcfffff,
//     0xff3fffff,
//     0xfcffffff

//   };





//   /*
// #define ELT0 0xffffffe7
// #define ELT1 0xffffff9f
// #define ELT2 0xfffffe7f
// #define ELT3 0xfffff9ff
// #define IDX0 0xffffe7ff
// #define IDX1 0xffff9fff
// #define IDX2 0xfffe7fff
// #define IDX3 0xfff9ffff
// */

//   /**********************************************
//    * HexStack
//    **********************************************/
//   /// Sparse set representation

//   class HexStack 
//   {
//   public:

//     /*!@name Parameters*/
//     //@{
//     int _data_;

//     // the first 2 bits stand for the current size (1-4)
//     // the following 8 bits store the values (0-3)
//     // the following 8 bits store the index (0-3)
//     // the following 8 bits store the suucessive sizes (1-4)
//     // the following 4 bits don't do anything
//     // the last 2 bits store the number of changes
//     //@}

//     /*!@name Constructors*/
//     //@{
//     HexStack()
//     {
//       _data_ = 0;
//     }

//     HexStack(const int n=1)
//     {
//       initialise(n);
//     }

//     virtual ~HexStack()
//     {
//     }

//     inline int size() { return (_data_ & 3)+1; }
//     inline int element(const int i) { return (_data_ >> (2+i)) & 3; }
//     inline int index(const int e) { return (_data_ >> (10+e)) & 3; }

//     virtual void initialise(const int n)
//     {
//       _data_ = n-1;

//       for(int i=0, int e=0; e<4; ++i, ++e) 
//   	{
//   	  // store the value
//   	  _data_ |= (e << (i+2));
//   	  // store the index
//   	  _data_ |= (i << (e+10));
//   	}
//     }
//     //@}    

//     /*!@name Accessors*/
//     //@{ 
//     inline bool contain(const int elt) const 
//     {
//       return ((_data_ >> (elt+10)) & 3)<=(_data_ & 3);
//     } 
  
//     inline bool singleton() const 
//     {
//       return !(_data_ & 3);
//     } 

//     inline int next(const int elt) const
//     {
//       unsigned int idx = ((_data_ >> (elt+10)) & 3);
//       return ((_data_ >> (2+idx+(idx < (_data_ & 3)))) & 3);
//     }
    
//     inline int operator[](const unsigned int idx) const
//     {
//       return (_data_ >> (2+idx)) & 3;
//     }
//     //@}

//     /*!@name List Manipulation*/
//     //@{
//     inline void fill()
//     {
//       _data_ &= SIZE;
//       _data_ |= 3;
//     }

//     inline void clear()
//     {
//       _data_ &= SIZE;
//     }
  
//     inline void remove(const int elt)
//     {
//       // get the index of the last element
//       unsigned int last_idx = (_data_ & 3);

//       // decrease the size
//       --_data_;

//       // idx of elt: 
//       unsigned int elt_idx = ((_data_ >> (10+elt)) & 3);

//       // idx of last: 
//       unsigned int last_elt = ((_data_ >> (2+last_idx)) & 3);
      
//       // set the index of the last element to the index of elt
//       _data_ &= INDEX[last_elt];
//       _data_ |= (elt_idx << (10+last_elt));


//       index_[list_[size]] = index_[elt];
//       list_[index_[elt]] = list_[size];
//       list_[size] = elt;
//       index_[elt] = size;
//     }

//     inline int next()
//     {
//       return list_[size];
//     }

//     inline int pop()
//     {
//       return list_[--size];
//     }

//     inline int pop_head()
//     {
//       --size;
//       index_[list_[size]] = 0;
//       const int elt = *list_;
//       *list_ = list_[size];
//       list_[size] = elt;
//       index_[elt] = size;
//       return elt;
//     }

//     inline int head()
//     {
//       return *list_;
//     }
    
//     inline int back()
//     {
//       return list_[size-1];
//     }

//     inline void add(const int elt)
//     {
//       //       std::cout << elt << ", " << size << " <= " << capacity << std::endl; 
//       //       std::cout << index_[elt] << " <= " << capacity << std::endl; 

//       index_[list_[size]] = index_[elt];
//       list_[index_[elt]] = list_[size];
//       list_[size] = elt;
//       index_[elt] = size;
//       ++size;
//     }

//     // create a new element that can potentially be outside the bounds
//     inline void create(const int elt)
//     {
//       extend(elt);
//       add(elt);
//     }

//     inline void ordered_add(const int elt)
//     {
//       // the first non-element goes where elt was
//       index_[list_[size]] = index_[elt];
//       list_[index_[elt]] = list_[size];

//       int idx = size;
//       while(idx && list_[idx-1] > elt) { // push every values greater than elt above elt
//   	list_[idx] = list_[idx-1];
//   	index_[list_[idx-1]] = idx;
//   	--idx;
//       }

//       list_[idx] = elt;
//       index_[elt] = idx;
//       ++size;
//     }


//     inline void revert_to(const int level)
//     {
//       size = level;
//     }

//     inline void index()
//     {
//       for(unsigned int i=0; i<capacity; ++i)
//   	index_[list_[i]] = i;
//     }
//     //@}

//     /*!@name Miscellaneous*/
//     //@{
//     //     std::string getString() const {
//     //       std::string return_str = "(";
//     //       if(size) return_str += toString(list_[0]);
//     //       for(unsigned int i=1; i<size; ++i)
//     // 	return_str += (" "+toString(list_[i]));
//     //       return_str += ")";
      
//     //       return return_str;
//     //     }

//     std::ostream& display(std::ostream& os) const {
//       os << "(";
//       if(size) os << list_[0];
//       for(unsigned int i=1; i<size; ++i)
//   	os << " " << list_[i];
//       os << ")";
//       return os;
//     }
//     //@}
//   };



  // template < int ARITY >
  // class ActiveStack {

  //   int _data_;
  //   //int _level_;

  //   ActiveStack(const int n) {
  //     _data_ = (1 << n)-1;
  //   }

  //   ~ActiveStack() {}

  //   void remove(const int x, const int lvl) 
  //   { 
  //     if(lvl == _level_) {
  // 	_data_ ^= (1 << x);
  //     } else {
  // 	_data_ = (_data_ << ARITY) | (_data_ ^ (1 << x)); }
  //   }


  // };


  // /**********************************************
  //  * Stack
  //  **********************************************/
  // /// Sparse set representation

  // template< class PTR_TYPE >
  // class Stack 
  // {
  // public:

  //   /*!@name Parameters*/
  //   //@{
  //   /// list of values
  //   PTR_TYPE *list_;
  //   /// current max capacity
  //   unsigned int capacity;
  //   /// current size
  //   unsigned int size;
  //   /// values' indices
  //   unsigned int *index_;
  //   int offset;
  //   //unsigned int *start_;
  //   //@}

  //   /*!@name Constructors*/
  //   //@{
  //   Stack()
  //   {
  //     size = 0;
  //     capacity = 0;
  //     list_ = NULL;
  //     offset = 0;
  //     index_ = NULL;
  //   }

  //   Stack(Vector< PTR_TYPE >& obj, bool full=true)
  //   {
  //     initialise(obj, full);
  //   }

  //   virtual ~Stack()
  //   {
  //     delete [] list_;
  //     index_  += offset;
  //     delete [] index_;
  //   }

  //   void initialise(Vector< PTR_TYPE >& obj, const bool full=true)
  //   {
  //     assert((obj.size == 0) || ((unsigned int)(obj.back()->id - obj[0]->id + 1) == obj.size));

  //     capacity = (obj.size ? obj.size : obj.capacity);
  //     list_ = new PTR_TYPE[capacity];
  //     offset = (obj.size ? obj[0]->id : 0);
  //     index_ = new unsigned int[capacity];
  //     index_ -= offset;
  //     for(unsigned int i=0; i<capacity; ++i) 
  // 	{
  // 	  index_[i+offset] = i;
  // 	  list_[i] = obj[i];
  // 	}
      
  //     size = (full ? capacity : 0);
  //   }
  //   //@}    

  //   /*!@name Accessors*/
  //   //@{  
  //   inline bool contain(const PTR_TYPE elt) const 
  //   {
  //     return index_[elt->id]<size;
  //   } 
  //   inline bool contain(const int elt) const 
  //   {
  //     return index_[elt]<size;
  //   } 
  
  //   inline bool empty()const 
  //   {
  //     return !size;
  //   } 

  //   inline PTR_TYPE next(const PTR_TYPE elt) const
  //   {
  //     unsigned int idx = index_[elt->id]+1;
  //     return (idx < size ? list_[idx] : elt);
  //   }
  //   inline PTR_TYPE next(const int elt) const
  //   {
  //     unsigned int idx = index_[elt]+1;
  //     return (idx < size ? list_[idx] : elt);
  //   }
    
  //   inline PTR_TYPE operator[](const unsigned int idx) const
  //   {
  //     return list_[idx];
  //   }

  //   inline PTR_TYPE& operator[](const unsigned int idx)
  //   {
  //     return list_[idx];
  //   }
  //   //@}

  //   /*!@name List Manipulation*/
  //   //@{
  //   inline void fill()
  //   {
  //     size = capacity;
  //   }

  //   inline void clear()
  //   {
  //     size = 0;
  //   }
  
  //   inline void set_to(const PTR_TYPE elt)
  //   {
  //     int idx = elt->id;
  //     size=1;
  //     index_[(*list_)->id] = index_[idx];
  //     list_[index_[idx]] = *list_;
  //     *list_ = elt;
  //     index_[idx] = 0;
  //   }

  //   inline void remove(const PTR_TYPE elt)
  //   {
  //     int idx = elt->id;
  //     --size;
  //     index_[list_[size]->id] = index_[idx];
  //     list_[index_[idx]] = list_[size];
  //     list_[size] = elt;
  //     index_[idx] = size;
  //   }

  //   inline void remove(const int idx)
  //   {
  //     PTR_TYPE elt = list_[index_[idx]];
  //     --size;
  //     index_[list_[size]->id] = index_[idx];
  //     list_[index_[idx]] = list_[size];
  //     list_[size] = elt;
  //     index_[idx] = size;
  //   }

  //   inline PTR_TYPE next()
  //   {
  //     return list_[size];
  //   }

  //   inline PTR_TYPE pop()
  //   {
  //     return list_[--size];
  //   }

  //   inline PTR_TYPE pop_head()
  //   {
  //     --size;
  //     index_[list_[size]->id] = 0;
  //     const PTR_TYPE elt = *list_;
  //     *list_ = list_[size];
  //     list_[size] = elt;
  //     index_[elt->id] = size;
  //     return elt;
  //   }

  //   inline PTR_TYPE head()
  //   {
  //     return *list_;
  //   }
    
  //   inline PTR_TYPE back()
  //   {
  //     return list_[size-1];
  //   }

  //   inline void add(const PTR_TYPE elt)
  //   {
  //     int idx = elt->id;
  //     index_[list_[size]->id] = index_[idx];
  //     list_[index_[idx]] = list_[size];
  //     list_[size] = elt;
  //     index_[idx] = size;
  //     ++size;
  //   }

  //   inline void ordered_add(const PTR_TYPE elt)
  //   {
  //     int idx = elt->id;

  //     // the first non-element goes where elt was
  //     index_[list_[size]->id] = index_[idx];
  //     list_[index_[idx]] = list_[size];

  //     int rank = size;
  //     while(idx && list_[rank-1]->id > elt->id) { // push every values greater than elt above elt
  // 	list_[rank] = list_[rank-1];
  // 	index_[list_[rank-1]->id] = rank;
  // 	--rank;
  //     }

  //     list_[rank] = elt;
  //     index_[idx] = rank;
  //     ++size;
  //   }

  //   inline void revert_to(const int level)
  //   {
  //     size = level;
  //   }

  //   inline void index()
  //   {
  //     for(unsigned int i=0; i<capacity; ++i)
  // 	index_[list_[i]->id] = list_[i]->id;
  //   }
  //   //@}

  //   /*!@name Miscellaneous*/
  //   //@{
  //   //     std::string getString() const {
  //   //       std::string return_str = "(";
  //   //       if(size) return_str += toString(list_[0]);
  //   //       for(unsigned int i=1; i<size; ++i)
  //   // 	return_str += (" "+toString(list_[i]));
  //   //       return_str += ")";
      
  //   //       return return_str;
  //   //     }

  //   std::ostream& display(std::ostream& os) const {
  //     os << "(";
  //     if(size) os << list_[0];
  //     for(unsigned int i=1; i<size; ++i)
  // 	os << " " << list_[i];
  //     os << ")";
  //     return os;
  //   }
  //   //@}
  // };


  /**********************************************
   * VarStack
   **********************************************/
  /// Sparse set representation

  class Solver;
  template< class VAR_TYPE, class SIZE_TYPE >
  class VarStack 
  {
  public:

    /*!@name Parameters*/
    //@{
    /// list of values
    VAR_TYPE *list_;
    /// current max capacity
    unsigned int capacity;
    /// current size
    SIZE_TYPE size;
    /// values' indices
    unsigned int *index_;
    int offset;
    //unsigned int *start_;
    //@}

    /*!@name Constructors*/
    //@{
    VarStack() : size(0)
    {
      //size = 0;
      capacity = 0;
      list_ = NULL;
      offset = 0;
      index_ = NULL;
    }

    VarStack(Vector< VAR_TYPE >& obj, bool full=true)
    {
      initialise(obj, full);
    }

    virtual ~VarStack()
    {
      delete [] list_;
      index_  += offset;
      delete [] index_;
    }

    void point_to(VarStack< VAR_TYPE, SIZE_TYPE >& vs) {
      capacity = vs.capacity;
      size = vs.size;
      list_ = vs.list_;
      offset = vs.offset;
      index_ = vs.index_;
    }

    void initialise(Solver *s) {
      size.initialise(s);
    }

    void initialise(const int n)
    {
      capacity = n;
      list_ = new VAR_TYPE[capacity];
      VAR_TYPE x;
      std::fill(list_, list_+capacity, x);
      
      offset = 0;
      index_ = new unsigned int[capacity];
      for(unsigned int i=0; i<capacity; ++i) 
	{
	  index_[i] = i;
	}
      index_ -= offset;
      size = 0;
    }

    void extend(const int idx) {
      //int idx = elt.id();
      int new_lb = offset;
      int new_ub = capacity-offset-1;
      
      if(idx < new_lb || idx > new_ub) {

	if(idx < new_lb) new_lb = idx;
	else if(idx > new_ub) new_ub = idx;
	
	unsigned int new_capacity = new_ub-new_lb+1;
	if(new_capacity < 2*capacity) new_capacity = 2*capacity;
	
	unsigned int *aux_index = index_+offset;
	index_ = new unsigned int[new_capacity];
	memcpy(index_, aux_index, capacity*sizeof(unsigned int));
	for(unsigned int i=capacity; i<new_capacity; ++i) 
	  {
	    index_[i] = i;
	  }
	delete [] aux_index;
	
	VAR_TYPE *aux_list = list_;
	list_ = new VAR_TYPE[new_capacity];
	memcpy(list_, aux_list, capacity*sizeof(VAR_TYPE));

	VAR_TYPE x;
	std::fill(list_+capacity, list_+new_capacity, x);

	delete [] aux_list;
	
	index_ -= new_lb;
	capacity = new_capacity;
	offset = new_lb;
      }
    }

    void declare(VAR_TYPE elt) {
      int idx = elt.id();
      extend(idx);
      
      if(idx < offset || idx >= offset+(int)size) {
	list_[index_[idx+offset]] = list_[size];
	list_[size] = elt;
	index_[idx+offset] = size;
	++size;
      } else {
	add(elt);
      }
    }

    void initialise(Vector< VAR_TYPE >& obj, const bool full=true)
    {
      assert((obj.size == 0) || ((unsigned int)(obj.back().id() - obj[0].id() + 1) == obj.size));

      capacity = (obj.size ? obj.size : obj.capacity);
      list_ = new VAR_TYPE[capacity];
      offset = (obj.size ? obj[0].id() : 0);
      index_ = new unsigned int[capacity];
      index_ -= offset;
      for(unsigned int i=0; i<capacity; ++i) 
	{
	  index_[i+offset] = i;
	  list_[i] = obj[i];
	}
      
      size = (full ? capacity : 0);
    }
    //@}    

    /*!@name Accessors*/
    //@{  
    inline bool contain(const VAR_TYPE elt) const 
    {
      return (int)(index_[elt.id()])<(int)size;
    } 
    inline bool safe_contain(const VAR_TYPE elt) const 
    {
      return (elt.id()-offset <= (int)capacity && (int)(index_[elt.id()])<(int)size);
    } 
    inline bool contain(const int elt) const 
    {
      return index_[elt]<(unsigned int)size;
    } 
    inline bool safe_contain(const int elt) const 
    {
      return (elt-offset <= capacity && index_[elt]<(unsigned int)size);
    } 
  
    inline bool empty()const 
    {
      return !(int)size;
    } 

    inline VAR_TYPE next(const VAR_TYPE elt) const
    {
      unsigned int idx = index_[elt.id()]+1;
      return (idx < size ? list_[idx] : elt);
    }
    inline VAR_TYPE next(const int elt) const
    {
      unsigned int idx = index_[elt]+1;
      return (idx < size ? list_[idx] : elt);
    }
    
    inline VAR_TYPE operator[](const unsigned int idx) const
    {
      return list_[idx];
    }

    inline VAR_TYPE& operator[](const unsigned int idx)
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
  
    inline void set_to(const VAR_TYPE elt)
    {
      int idx = elt.id();
      size=1;
      index_[(*list_).id()] = index_[idx];
      list_[index_[idx]] = *list_;
      *list_ = elt;
      index_[idx] = 0;
    }

    inline void remove(const VAR_TYPE elt)
    {
      int idx = elt.id();
      --size;
      index_[list_[size].id()] = index_[idx];
      list_[index_[idx]] = list_[size];
      list_[size] = elt;
      index_[idx] = size;
    }

    inline void remove(const int idx)
    {
      VAR_TYPE elt = list_[index_[idx]];
      --size;
      index_[list_[size].id()] = index_[idx];
      list_[index_[idx]] = list_[size];
      list_[size] = elt;
      index_[idx] = size;
    }

    inline VAR_TYPE next()
    {
      return list_[size];
    }

    inline VAR_TYPE pop()
    {
      return list_[--size];
    }

    inline VAR_TYPE pop_head()
    {
      --size;
      index_[list_[size].id()] = 0;
      const VAR_TYPE elt = *list_;
      *list_ = list_[size];
      list_[size] = elt;
      index_[elt.id()] = size;
      return elt;
    }

    inline VAR_TYPE head()
    {
      return *list_;
    }
    
    inline VAR_TYPE back()
    {
      return list_[size-1];
    }

    inline void add(const VAR_TYPE elt)
    {
      
      // std::cout << "add " << elt <<  " to " ;
      // display(std::cout);
      // std::cout << " " << (int*)list_ << " " << index_ << std::endl;

      int idx = elt.id();
      index_[list_[size].id()] = index_[idx];
      list_[index_[idx]] = list_[size];
      list_[size] = elt;
      index_[idx] = size;
      ++size;
    }

    inline void ordered_add(const VAR_TYPE elt)
    {
      int idx = elt.id();

      // the first non-element goes where elt was
      index_[list_[size].id()] = index_[idx];
      list_[index_[idx]] = list_[size];

      int rank = size;
      while(idx && list_[rank-1].id() > elt.id()) { // push every values greater than elt above elt
	list_[rank] = list_[rank-1];
	index_[list_[rank-1].id()] = rank;
	--rank;
      }

      list_[rank] = elt;
      index_[idx] = rank;
      ++size;
    }

    inline void revert_to(const int level)
    {
      size = level;
    }

    inline int index(const int elt_idx) {
      if(elt_idx>((int)capacity-offset)) return -1;
      return index_[elt_idx];
    }

    inline void index()
    {
      for(unsigned int i=0; i<capacity; ++i)
	index_[list_[i].id()] = list_[i].id();
    }
    //@}

    /*!@name Miscellaneous*/
    //@{
    //     std::string getString() const {
    //       std::string return_str = "(";
    //       if(size) return_str += toString(list_[0]);
    //       for(unsigned int i=1; i<size; ++i)
    // 	return_str += (" "+toString(list_[i]));
    //       return_str += ")";
      
    //       return return_str;
    //     }

    std::ostream& display(std::ostream& os) const {
      os << "(";
      if(size) os << list_[0];
      for(int i=1; i<size; ++i)
	os << " " << list_[i];
      os << ")";
      return os;
    }
    //@}
  };



  /**********************************************
   * ConStack
   **********************************************/
  /// Sparse set representation

  template< class CON_TYPE >
  class ConStack 
  {
  public:

    /*!@name Parameters*/
    //@{
    /// list of values
    CON_TYPE *list_;
    /// current max capacity
    unsigned int capacity;
    /// current size
    unsigned int size;
    /// values' indices
    unsigned int *index_;
    int offset;
    //unsigned int *start_;
    //@}

    /*!@name Constructors*/
    //@{
    ConStack()
    {
      size = 0;
      capacity = 0;
      list_ = NULL;
      offset = 0;
      index_ = NULL;
    }

    ConStack(Vector< CON_TYPE >& obj, bool full=true)
    {
      initialise(obj, full);
    }

    virtual ~ConStack()
    {
      delete [] list_;
      index_  += offset;
      delete [] index_;
    }

    void initialise(const int n)
    {
      capacity = n;
      list_ = new CON_TYPE[capacity];
      offset = 0;
      index_ = new unsigned int[capacity];
      for(unsigned int i=0; i<capacity; ++i) 
	{
	  index_[i] = i;
	}
      index_ -= offset;
      size = 0;
    }

    void extend(const int idx) {
      //int idx = elt.id();
      int new_lb = offset;
      int new_ub = capacity-offset-1;
      
      if(idx < new_lb || idx > new_ub) {

	if(idx < new_lb) new_lb = idx;
	else if(idx > new_ub) new_ub = idx;
	
	unsigned int new_capacity = new_ub=new_lb+1;
	if(new_capacity < 2*capacity) new_capacity = capacity;
	
	unsigned int *aux_index = index_+offset;
	index_ = new unsigned int[new_capacity];
	memcpy(index_, aux_index, capacity*sizeof(unsigned int));
	delete [] aux_index;
	
	CON_TYPE *aux_list = list_;
	list_ = new CON_TYPE[new_capacity];
	memcpy(list_, aux_list, capacity*sizeof(CON_TYPE));
	delete [] aux_list;

	index_ -= new_lb;
	capacity = new_capacity;
      }
    }

    void declare(CON_TYPE elt) {
      int idx = elt->id;
      extend(idx);
      
      if(idx < offset || idx >= offset+(int)size) {
	list_[index_[idx+offset]] = list_[size];
	list_[size] = elt;
	index_[idx+offset] = size;
	++size;
      } else {
	add(elt);
      }
    }

    void initialise(Vector< CON_TYPE >& obj, const bool full=true)
    {
      assert((obj.size == 0) || ((unsigned int)(obj.back()->id - obj[0]->id + 1) == obj.size));

      capacity = (obj.size ? obj.size : obj.capacity);
      list_ = new CON_TYPE[capacity];
      offset = (obj.size ? obj[0]->id : 0);
      index_ = new unsigned int[capacity];
      index_ -= offset;
      for(unsigned int i=0; i<capacity; ++i) 
	{
	  index_[i+offset] = i;
	  list_[i] = obj[i];
	}
      
      size = (full ? capacity : 0);
    }
    //@}    

    /*!@name Accessors*/
    //@{  
    inline bool contain(const CON_TYPE elt) const 
    {
      return index_[elt->id]<size;
    } 
    inline bool contain(const int elt) const 
    {
      return index_[elt]<size;
    } 
  
    inline bool empty()const 
    {
      return !size;
    } 

    inline CON_TYPE next(const CON_TYPE elt) const
    {
      unsigned int idx = index_[elt->id]+1;
      return (idx < size ? list_[idx] : elt);
    }
    inline CON_TYPE next(const int elt) const
    {
      unsigned int idx = index_[elt]+1;
      return (idx < size ? list_[idx] : elt);
    }
    
    inline CON_TYPE operator[](const unsigned int idx) const
    {
      return list_[idx];
    }

    inline CON_TYPE& operator[](const unsigned int idx)
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
  
    inline void set_to(const CON_TYPE elt)
    {
      int idx = elt->id;
      size=1;
      index_[(*list_)->id] = index_[idx];
      list_[index_[idx]] = *list_;
      *list_ = elt;
      index_[idx] = 0;
    }

    inline void remove(const CON_TYPE elt)
    {
      int idx = elt->id;
      --size;
      index_[list_[size]->id] = index_[idx];
      list_[index_[idx]] = list_[size];
      list_[size] = elt;
      index_[idx] = size;
    }

    inline void remove(const int idx)
    {
      CON_TYPE elt = list_[index_[idx]];
      --size;
      index_[list_[size]->id] = index_[idx];
      list_[index_[idx]] = list_[size];
      list_[size] = elt;
      index_[idx] = size;
    }

    inline CON_TYPE next()
    {
      return list_[size];
    }

    inline CON_TYPE pop()
    {
      return list_[--size];
    }

    inline CON_TYPE pop_head()
    {
      --size;
      index_[list_[size]->id] = 0;
      const CON_TYPE elt = *list_;
      *list_ = list_[size];
      list_[size] = elt;
      index_[elt->id] = size;
      return elt;
    }

    inline CON_TYPE head()
    {
      return *list_;
    }
    
    inline CON_TYPE back()
    {
      return list_[size-1];
    }

    inline void add(const CON_TYPE elt)
    {
      int idx = elt->id;
      index_[list_[size]->id] = index_[idx];
      list_[index_[idx]] = list_[size];
      list_[size] = elt;
      index_[idx] = size;
      ++size;
    }

    inline void ordered_add(const CON_TYPE elt)
    {
      int idx = elt->id;

      // the first non-element goes where elt was
      index_[list_[size]->id] = index_[idx];
      list_[index_[idx]] = list_[size];

      int rank = size;
      while(idx && list_[rank-1]->id > elt->id) { // push every values greater than elt above elt
	list_[rank] = list_[rank-1];
	index_[list_[rank-1]->id] = rank;
	--rank;
      }

      list_[rank] = elt;
      index_[idx] = rank;
      ++size;
    }

    inline void revert_to(const int level)
    {
      size = level;
    }

    inline void index()
    {
      for(unsigned int i=0; i<capacity; ++i)
	index_[list_[i]->id] = list_[i]->id;
    }
    //@}

    /*!@name Miscellaneous*/
    //@{
    //     std::string getString() const {
    //       std::string return_str = "(";
    //       if(size) return_str += toString(list_[0]);
    //       for(unsigned int i=1; i<size; ++i)
    // 	return_str += (" "+toString(list_[i]));
    //       return_str += ")";
      
    //       return return_str;
    //     }

    std::ostream& display(std::ostream& os) const {
      os << "(";
      if(size) os << list_[0];
      for(unsigned int i=1; i<size; ++i)
	os << " " << list_[i];
      os << ")";
      return os;
    }
    //@}
  };

 
  // /**********************************************
  //  * IdStack
  //  **********************************************/
  // /// Sparse set representation
  // template< class DATA_TYPE >
  // class IdStack 
  // {
  // public:

  //   /*!@name Parameters*/
  //   //@{
  //   /// list of values
  //   DATA_TYPE *list_;
  //   /// current max capacity
  //   unsigned int capacity;
  //   /// current size
  //   unsigned int size;
  //   //@}

  //   /*!@name Constructors*/
  //   //@{
  //   IdStack()
  //   {
  //     size = 0;
  //     capacity = 0;
  //     list_ = NULL;
  //   }

  //   IdStack(Vector< DATA_TYPE >& obj, bool full=true)
  //   {
  //     initialise(obj, full);
  //   }

  //   virtual ~IdStack()
  //   {
  //     delete [] list_;
  //   }

  //   void initialise(Vector< DATA_TYPE >& obj, const bool full=true)
  //   {
  //     capacity = obj.size;
  //     size = (full ? capacity : 0);
  //     list_ = new DATA_TYPE[capacity];
  //     for(unsigned int i=0; i<capacity; ++i) 
  // 	{
  // 	  list_[i] = obj[i];
  // 	  list_[i].set_stack_id(i);
  // 	}
  //   }
  //   //@}    

  //   /*!@name Accessors*/
  //   //@{  
  //   inline bool contain(const DATA_TYPE elt) const 
  //   {
  //     return elt.get_stack_id()<size;
  //   } 
  
  //   inline bool empty()const 
  //   {
  //     return !size;
  //   } 

  //   inline DATA_TYPE next(const DATA_TYPE elt) const
  //   {
  //     unsigned int idx = elt.get_stack_id()+1;
  //     return (idx < size ? list_[idx] : elt);
  //   }
    
  //   inline DATA_TYPE operator[](const unsigned int idx) const
  //   {
  //     return list_[idx];
  //   }

  //   inline DATA_TYPE& operator[](const unsigned int idx)
  //   {
  //     return list_[idx];
  //   }
  //   //@}

  //   /*!@name List Manipulation*/
  //   //@{
  //   inline void fill()
  //   {
  //     size = capacity;
  //   }

  //   inline void clear()
  //   {
  //     size = 0;
  //   }
  
  //   inline void set_to(DATA_TYPE elt)
  //   {
  //     int idx = elt.get_stack_id();
  //     size=1;
      
  //     (*list_).set_stack_id(idx);
  //     elt.set_stack_id(0);

  //     list_[idx] = *list_;
  //     *list_ = elt;
  //   }

  //   inline void remove(DATA_TYPE elt)
  //   {
  //     int idx = elt.get_stack_id();

  //     list_[--size].set_stack_id(idx);
  //     elt.set_stack_id(size);
      
  //     list_[idx] = list_[size];
  //     list_[size] = elt;
  //   }

  //   inline DATA_TYPE pop()
  //   {
  //     return list_[--size];
  //   }

  //   inline DATA_TYPE pop_head()
  //   {
  //     --size;
  //     const DATA_TYPE elt = *list_;

  //     list_[size].set_stack_id(0);
  //     elt.set_stack_id(size);

  //     *list_ = list_[size];
  //     list_[size] = elt;
  //   }

  //   inline DATA_TYPE head()
  //   {
  //     return *list_;
  //   }
    
  //   inline DATA_TYPE back(const int offset=0)
  //   {
  //     return list_[size-1+offset];
  //   }

  //   inline void add(DATA_TYPE elt)
  //   {
  //     int idx = elt.get_stack_id();
      
  //     elt.set_stack_id(size);
  //     list_[size].set_stack_id(idx);

  //     list_[idx] = list_[size];
  //     list_[size] = elt;
      
  //     ++size;
  //   }

  //   inline void revert_to(const int level)
  //   {
  //     size = level;
  //   }

  //   std::ostream& display(std::ostream& os) const {
  //     os << "(";
  //     if(size) os << list_[0];
  //     for(unsigned int i=1; i<size; ++i)
  // 	os << " " << list_[i];
  //     os << ")";
  //     return os;
  //   }
  //   //@}
  // };


  template< class DATA_TYPE >
  class Indexed 
  {
    
  public:
    
    DATA_TYPE element;
    int index;

    Indexed() {};
    Indexed(DATA_TYPE x) { element = x; index = INFTY; };
    virtual ~Indexed() {};
  };

  /**********************************************
   * Stack
   **********************************************/
  /// Sparse set representation
  template< class DATA_TYPE >
  class Stack 
  {
  public:

    /*!@name Parameters*/
    //@{
    /// list of values
    Indexed< DATA_TYPE > *list_;
    /// current max capacity
    unsigned int capacity;
    /// current size if full
    unsigned int fsize;
    /// current size
    unsigned int size;
    //@}

    /*!@name Constructors*/
    //@{
    Stack()
    {
      size = 0;
      fsize = 0;
      capacity = 0;
      list_ = NULL;
    }

    Stack(Vector< DATA_TYPE >& obj, bool full=true)
    {
      initialise(obj, full);
    }

    virtual ~Stack()
    {
      delete [] list_;
    }

    void extend()
    {
      capacity = 2*fsize;
      Indexed<DATA_TYPE> *new_list_ = new Indexed<DATA_TYPE>[capacity];
      memcpy(new_list_, list_, sizeof(Indexed<DATA_TYPE>)*fsize);
      delete [] list_;
      list_ = new_list_;
    }

    void initialise(Vector< DATA_TYPE >& obj, const bool full=true)
    {
      fsize = obj.size;
      size = (full ? fsize : 0);
      capacity = 2*fsize;
      list_ = new Indexed<DATA_TYPE>[capacity];
      for(unsigned int i=0; i<fsize; ++i) 
	{
	  list_[i] = Indexed<DATA_TYPE>(obj[i]);
	  list_[i].index = i;
	}
    }
    //@}    

    /*!@name Accessors*/
    //@{  
    inline bool contain(const Indexed<DATA_TYPE> elt) const 
    {
      return elt.index<size;
    } 
  
    inline bool empty()const 
    {
      return !size;
    } 

    inline DATA_TYPE next(const Indexed<DATA_TYPE> elt) const
    {
      unsigned int idx = elt.index+1;
      return (idx < size ? list_[idx] : elt);
    }
    
    inline DATA_TYPE operator[](const unsigned int idx) const
    {
      return list_[idx].element;
    }

    inline DATA_TYPE& operator[](const unsigned int idx)
    {
      return list_[idx].element;
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
  
    inline void set_to(Indexed<DATA_TYPE> elt)
    {
      int idx = elt.index;
      size=1;
      
      list_->index = idx;
      elt.index = 0;

      list_[idx] = *list_;
      *list_ = elt;
    }

    inline void remove(Indexed<DATA_TYPE> elt)
    {
      int idx = elt.index;

      list_[--size].index = idx;
      elt.index = size;
      
      list_[idx] = list_[size];
      list_[size] = elt;
    }

    inline DATA_TYPE pop()
    {
      return list_[--size].element;
    }

    inline DATA_TYPE pop_head()
    {
      --size;
      const DATA_TYPE elt = list_->element;

      list_[size].index = 0;
      elt.index = size;

      *list_ = list_[size];
      list_[size] = elt;
    }

    inline DATA_TYPE head()
    {
      return list_->element;
    }
    
    inline DATA_TYPE back(const int offset=0)
    {
      return list_[size-1+offset].element;
    }

    inline void add(Indexed<DATA_TYPE> elt)
    {
      int idx = elt.index;
      
      elt.index = size;
      list_[size].index = idx;

      list_[idx] = list_[size];
      list_[size] = elt;
      
      ++size;
    }

    inline void declare(DATA_TYPE elt)
    {
      Indexed<DATA_TYPE> x(elt, fsize);
      if(fsize == capacity) extend();
      list_[fsize++] = x;
      
      add(x);
    }

    inline void revert_to(const int level)
    {
      size = level;
    }

    std::ostream& display(std::ostream& os) const {
      os << "(";
      if(size) os << list_[0].element;
      for(unsigned int i=1; i<size; ++i)
	os << " " << list_[i].element;
      os << ")";
      return os;
    }
    //@}
  };


  /**********************************************
   * Node
   **********************************************/
  /// Node of Multilist

  template < class DATA_TYPE >
  class Node {
  public:

    /*!@name Parameters*/
    //@{
    int prev;
    int next;
    DATA_TYPE elt;
    //@}

    Node(DATA_TYPE e=(DATA_TYPE(0)), int p=-1, int n=-1) {elt = e; prev=p; next=n;}
    inline operator const DATA_TYPE() const { return elt; }

    //inline Node& operator=(const Node& nd) const { prev = nd.prev; return elt; }

    /*!@name Miscellaneous*/
    //@{
    /// Print on out stream
    //     std::string getString() const {
    //       std::string return_str = "("+toString(prev)+"|"+toString(elt)+"|"+toString(next)+")";
    //       return return_str;
    //     }
    std::ostream& display(std::ostream& os) const {
      //os << "(" << prev << "|" << elt << "|" << next << ")";
      os << elt;
      return os;
    }
    //@}
  };


  /**********************************************
   * Multilist
   **********************************************/
  /// List with multiple entry points

  template < class DATA_TYPE, int NUM_HEAD >
  class MultiList 
  {

  public:
    /*!@name Parameters*/
    //@{
    Vector< Node< DATA_TYPE > > data;
    int head[NUM_HEAD+1]; // index of the first element of the list
    unsigned int degree; // number of element in the entire list
    //@}

    /*!@name Constructors*/
    //@{
    MultiList() {
      data.initialise(0, std::max(2*NUM_HEAD,16));
      for(int k=0; k<=NUM_HEAD; ++k) {
	head[k] = k; 
	Node< DATA_TYPE > x((DATA_TYPE)0, k-1, k+1); 
	data.add(x);
      }
      degree=0;
    }
    //@}

    /*!@name List Manipulation*/
    //@{
    inline int create(DATA_TYPE elt, const int k=0) {
      Node<DATA_TYPE> x(elt);
      data.add(x);
      add(data.size-1, k);
      return data.size-1;
    }

    inline void add(const int idx, const int k=0) {
      int succ = data.stack_[head[k]].next;
      data.stack_[succ].prev = idx;
      data.stack_[head[k]].next = idx;
      data.stack_[idx].next = succ;
      data.stack_[idx].prev = head[k];
      ++degree;
    }

    inline void remove(const int idx, const int k=0) {
      int succ = data.stack_[idx].next;
      int prec = data.stack_[idx].prev;   

      data.stack_[succ].prev = prec;
      data.stack_[prec].next = succ;

      --degree;
    }

    inline Node<DATA_TYPE>& first(const int h) {
      return data.stack_[head[h]];
    }

    inline Node<DATA_TYPE>& last(const int h) {
      return data.stack_[data.stack_[head[h+1]].prev];
    }

    inline DATA_TYPE pop(const int h) {
      int idx = data.stack_[head[h+1]].prev;
      remove(idx, h);
      return data.stack_[idx].elt;
    }

    inline bool empty(const int h) {
      return (data.stack_[head[h]].next > NUM_HEAD);
    }

    inline bool next(Node< DATA_TYPE >& node) const {
      unsigned int idx_next;
      do { 
	idx_next = node.next;
	if(idx_next == NUM_HEAD) return false;
	node = data.stack_[idx_next];
      } while(!node.elt);
      return true;
    }

    inline bool prev(Node< DATA_TYPE >& node) const {
      unsigned int idx_prev;
      do { 
	idx_prev = node.prev;
	if(!idx_prev) return false;
	node = data.stack_[idx_prev];
      } while(!node.elt);
      return true;
    }
    //@}

    /*!@name Miscellaneous*/
    //@{
    /// Print on out stream
    //     std::string getString() const {
    //       std::string return_str = "";
    //       for(unsigned int i=0; i<data.size; ++i)
    // 	return_str += toString(data[i]);
    //       return return_str;
    //     }
    std::ostream& display(std::ostream& os) const;//  {

    void debug_print(std::ostream& os) const 
    {
      for(unsigned int i=head[0]; i<data.size; ++i)
	os << "(" << data.stack_[i].prev << "|" 
	   << data.stack_[i].elt << "|" 
	   << data.stack_[i].next << ")";
      //o << data.stack_[i].print(o);
    }
    void flat_print(std::ostream& o) const 
    {
      o << "[ ";
      for(unsigned int i=NUM_HEAD+1; i<data.size; ++i)
	o << data.stack_[i].elt << " ";
      o << "] ";
    }
    void print(std::ostream& o) const 
    {
      for(int k=0; k<NUM_HEAD; ++k) {
	o << "[ ";

	Node< DATA_TYPE > nd = data.stack_[head[k]];
	while(next(nd)) {
	  o << (DATA_TYPE)nd << " ";
	}
	o << "] ";
      }
    }
    //@}

  };


  /**********************************************
   * BitSet
   **********************************************/

  /*! \class BitSet
    \brief A representation of sets using a vector of bits.

    The sets have a static capacity. 
    Template are used so that ReversibleWords can used instead of unsigned int
  */
  template< class WORD_TYPE, class FLOAT_TYPE >
  class Bitset
  { 

  public:
    
    /*!@name Class attributes*/
    //@{
    static const WORD_TYPE empt = 0;
    static const WORD_TYPE full = ~0;
    static const unsigned int EXP = (sizeof(empt) == 4 ? 5 /*32 bits*/ : 6 /*64 bits*/);
    static const unsigned int size_word_bit = (1 << EXP);
    static const unsigned int size_word_byte = (size_word_bit >> 3);
    static const unsigned int CACHE = (size_word_bit - 1);
    static const unsigned int LASTCHAR = (size_word_bit - 8);
    static const unsigned int mantissa = (sizeof(empt) == 4 ? 23 /*32 bits*/ : 52 /*64 bits*/);
    static const unsigned int float_offset = (sizeof(empt) == 4 ? (0x7f) /*32 bits*/ : (0x3ff) /*64 bits*/);
    static const WORD_TYPE mask_first_char = 0xff;
    static const WORD_TYPE mask_last_char = (mask_first_char << ((size_word_bit) - 8));
    //@}

    /*!@name Parameters*/
    //@{
    /// index of the first word used to represent the set
    int pos_words;
    /// 1 + index of the last word used to represent the set
    int neg_words;
    /// A vector of bits 
    WORD_TYPE* table;
    //@}

    Bitset()
    {
      initialise();
    }

    Bitset(int sz)
    {
      if(sz>0) {
	initialise(sz,0);
      } else {
      	initialise();
      }
    }

    void initialise()
    {
      pos_words = 0;
      neg_words = 0;
      table = NULL;
    }

    Bitset(const int sz, const int* elt) 
    {
      int lb =  NOVAL;
      int ub = -NOVAL;
      for(int i=0; i<sz; ++i) {
	if(elt[i] > ub) ub = elt[i];
	if(elt[i] < lb) lb = elt[i];
      }

      initialise(lb,ub,empt);

      for(int i=0; i<sz; ++i) 
	add( elt[i] );
    }

    void initialise(const Vector< int >& elt) 
    {
      int min = elt.front();
      // int max = elt.front();
      // for(unsigned int i=1; i<elt.size; ++i) {
      // 	if(elt[i] < min) min = elt[i];
      // 	if(elt[i] > max) max = elt[i];
      // }
      int max = elt.back();

      initialise(min,max,empt);
      
      for(unsigned int i=0; i<elt.size; ++i) 
	add( elt[i] );
    }

    void initialise(const int lb, const int ub, const Vector< int >& elt) 
    {

      // std::cout << "initialise bitset with " << lb << ".." << ub << ": " ;
      // elt.display(std::cout);
      // std::cout << std::endl; 

      initialise(lb,ub,empt);

      // display(std::cout);
      // std::cout << std::endl;

      for(unsigned int i=0; i<elt.size; ++i) {
	add( elt[i] );
	
	// display(std::cout);
	// std::cout << std::endl;
      }
      
    }

    Bitset(const int lb, const int ub, const WORD_TYPE p)
    {
      initialise(lb,ub,p,NULL);
    }

    inline int word_index(const int elt) const
    {
      return (elt >> EXP);
    }

    bool operator==(const Bitset<WORD_TYPE,FLOAT_TYPE>& s) {
      return equal(s);
    }

    bool operator!=(const Bitset<WORD_TYPE,FLOAT_TYPE>& s) {
      return !equal(s);
    }

    Bitset<WORD_TYPE,FLOAT_TYPE>& operator=(const Bitset<WORD_TYPE,FLOAT_TYPE>& q) 
    {
      if(!table)
	clone(q);
      else
	copy(q);
      return *this;
    }

    void reinitialise(const int lb, const int ub, const WORD_TYPE p) 
    {
      table += neg_words;
      delete [] table;
      initialise(lb, ub, p, NULL);
    }

    void initialise(const int sz, const WORD_TYPE p) 
    {
      pos_words = sz;
      neg_words = 0;

      if( sz>=0 ) {
	table = new WORD_TYPE[pos_words];
	for(int i=0; i<pos_words; ++i) 
	  table[i]=p;
      } else table = NULL;
    }

    void initialise(const int lb, const int ub, const WORD_TYPE p, WORD_TYPE *pool=NULL) 
    {
      neg_words = (lb >> EXP);
      pos_words = (ub >> EXP)+1;
      if(pool==NULL) table = new WORD_TYPE[pos_words-neg_words];
      else table = pool;
      for(int i=0; i<pos_words-neg_words; ++i) 
	table[i]=p;
      table[pos_words-neg_words-1] &= 
	(p >> (size_word_bit-1-(ub & CACHE)));
      table[0] &= (p << (lb & CACHE));
      table -= neg_words;
    }

    void initialise(const Bitset<WORD_TYPE,FLOAT_TYPE>& s) 
    {
      pos_words = s.pos_words;
      neg_words = s.neg_words;
  
      table = new WORD_TYPE[pos_words-neg_words];
      table -= neg_words;
      for(int i=neg_words; i<pos_words; ++i) 
	//table[i].initialise(s.table+i, s.size(i));
	table[i] = s.table[i];
    }

    inline void declare(const int elt)
    {
      int i = (elt >> EXP);
      if( (i < neg_words) ||
	  (i >=  pos_words) ) 
	{
	  extend(elt);
	}
      fast_add(elt);
    }

    void extend(const int elt) 
    {
      int nval = (elt >> EXP);
      if( (nval < neg_words) ||
	  (nval >=  pos_words) ) 
	{
	  int new_neg_words = neg_words;
	  //nval;
	  int new_pos_words = pos_words;
	  //nval+1;
	  bool need_to_extend = false;
	  if(nval < new_neg_words) {
	    new_neg_words = nval;
	    need_to_extend = true;
	  }
	  if(nval >= new_pos_words) {
	    new_pos_words = nval+1;
	    need_to_extend = true;
	  }

	  if(need_to_extend) {
	    WORD_TYPE *aux = table;
	    table = new WORD_TYPE[new_pos_words-new_neg_words];
	    table -= new_neg_words;
	    
	    memcpy(table+neg_words, aux+neg_words, 
		   (pos_words-neg_words)*sizeof(WORD_TYPE));

	    if(new_neg_words < neg_words)
	      std::fill(table+new_neg_words, table+neg_words, 0);

	    if(new_pos_words > pos_words)
	      std::fill(table+pos_words, table+new_pos_words, 0);
	    
	    aux += neg_words;
	    delete [] aux;
	    
	    pos_words = new_pos_words; 
	    neg_words = new_neg_words; 
	  }
	}
    }

    Bitset(const Bitset<WORD_TYPE,FLOAT_TYPE>& s) 
    {      
      initialise();
      clone( s );
    }

    void clone(const Bitset<WORD_TYPE,FLOAT_TYPE>& s)
    {
      if(table) {
	table += neg_words;
	delete [] table;
      }
      neg_words = s.neg_words;
      pos_words = s.pos_words;
      table = new WORD_TYPE[pos_words-neg_words];
      memcpy(table, s.table+neg_words,
	     size_word_byte*(pos_words-neg_words));
      table -= neg_words;
    }

    void point_to(Bitset<WORD_TYPE,FLOAT_TYPE>& s)
    {
      neg_words = s.neg_words;
      pos_words = s.pos_words;
      table = s.table;
    }

    void point_to(WORD_TYPE *t)
    {
      neg_words = 0;
      pos_words = 1;
      table = t;
    }

    void copy(const Bitset<WORD_TYPE,FLOAT_TYPE>& s) 
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

    bool is_built()
    {
      return (table != NULL);
    }

    inline void swap(Bitset<WORD_TYPE,FLOAT_TYPE>& s)
    {
      WORD_TYPE *aux = s.table;
      s.table = table;
      table = aux;
    }


    void iterate_into_b(const int size, int *buffer) {
      int elt;
      int idx;
      int nval;

      union {FLOAT_TYPE f; WORD_TYPE i; } t;
      WORD_TYPE b;
      WORD_TYPE v; 

      nval = 1;
      elt = buffer[0];
      idx = ((elt+1) >> EXP);
      v = (table[idx] & (full << ((elt+1) & CACHE)));  
      
      while(nval < size) {
	// find the next word that is not null
	while(!v) v = table[++idx];

	// find the first element in the set:
	// remove all other element
	b = v & -v;

	// cast into float, which will be coded as 1*2^exp, and 'exp' is precisely the index of the first element 
	t.f = (FLOAT_TYPE)b; // cast the least significant bit in v to a float
	// keep only the exponant part
	elt = t.i >> mantissa;

	elt += idx * size_word_bit - float_offset;

	do {
	  // put the element in the buffer
	  buffer[nval] = elt;
	  ++nval;
	
	  // remove it from v
	  v ^= b;
	  
	  do {
	    
	    // try the next element
	    b <<= 1;
	    ++elt;

	  } while( b && !(v & b) );

	} while(v);
	
      }
      
    }


    void iterate_into(const int size, int *buffer) {
      int elt;
      int idx;
      int nval;

      union {FLOAT_TYPE f; WORD_TYPE i; } t;
      WORD_TYPE b;
      WORD_TYPE v; 

      nval = 1;
      elt = buffer[0];
      idx = ((elt+1) >> EXP);
      v = (table[idx] & (full << ((elt+1) & CACHE)));  
      
      elt = (idx * size_word_bit - float_offset);

      while(nval < size) {
	// find the next word that is not null
	while(!v) {
	  v = table[++idx];
	  elt += size_word_bit;
	}

	// find the first element in the set:
	do {
	  // remove all other element
	  b = v & -v;
	  
	  // cast into float, which will be coded as 1*2^exp, and 'exp' is precisely the index of the first element 
	  t.f = (FLOAT_TYPE)b; // cast the least significant bit in v to a float
	  	  	  
	  // put the element in the buffer
	  buffer[nval++] = (t.i >> mantissa)+elt;
	  
	  // remove it from v
	  v ^= b;
	  
	} while(v);      
      }
    }


    inline int lsb_mantissa(const WORD_TYPE v) const {
      union {FLOAT_TYPE f; WORD_TYPE i; } t;
      WORD_TYPE b = v & -v;

      t.f = (FLOAT_TYPE)b; // cast the least significant bit in v to a float
      b = t.i >> mantissa;
      return b - float_offset;
    }

    inline int lsb_gcc(const WORD_TYPE v) const {
      return __builtin_ctz(v);
    }

    inline int msb_gcc(const WORD_TYPE v) const {
      return __builtin_clz(v);
    }

      
    inline int minimum_element(int idx, WORD_TYPE v, const int def=NOVAL) const
    {
           
      while(v == 0) {
	if( ++idx >= pos_words )
	  return def;
	v = table[idx];
      }

      
      /*
      
      union {FLOAT_TYPE f; WORD_TYPE i; } t;
      WORD_TYPE b = v & -v;

      t.f = (FLOAT_TYPE)b; // cast the least significant bit in v to a float
      b = t.i >> mantissa;
      */
      //return b + idx * size_word_bit - float_offset;


     return __builtin_ctz(v) //__builtin_ffs(v) - 1 
       + (idx * size_word_bit);
      
    }

    /*!
      Minimum element in the set [O(N/32)]. 
    */
    inline int min() const
    { 
      int idx = neg_words;
      WORD_TYPE v = table[idx];
      return minimum_element(idx,v);
    }

    /*!
      Maximum element in the set [O(N/8)]. 
    */
    inline int max() const
    { 
      WORD_TYPE tab;
      int i=pos_words, j, k;
    
      while( i-- > neg_words )
	if( (tab = table[i]) ) {
	  j = size_word_byte;
	  while( j-- ) {
	    if( (k = getlast[(tab & mask_last_char) >> LASTCHAR]) >= 0 ) 
	      return ( (i<<EXP)+(j<<3)+k );	
	    tab = (tab << 8);
	  }
	}
      return NOVAL;
    }

    inline void  remove(const int elt) 
    {
      int i = (elt >> EXP);
      if( (i >= neg_words) && 
	  (i <  pos_words) )
	table[i] &= (full ^ ((WORD_TYPE)1 << (elt & CACHE)));
    }

    inline void  fast_remove(const int elt) 
    {
      table[(elt >> EXP)] &= (full ^ ((WORD_TYPE)1 << (elt & CACHE)));
    }

    inline void  word_remove(const int elt) 
    {
      table[neg_words] &= (full ^ ((WORD_TYPE)1 << (elt & CACHE)));
    }


   inline int next(const int elt) const {
      int idx = ((elt+1) >> EXP);
      if(idx >= pos_words) return elt;
      WORD_TYPE v = (table[idx] & (full << ((elt+1) & CACHE)));
      while(v == 0) {
	if(++idx >= pos_words) return elt;
	v = table[idx];
      }
      return lsb_gcc(v) + (idx * size_word_bit);
    }


    /*
    inline int next(const int elt) const {
      int idx = ((elt+1) >> EXP);
      if(idx >= pos_words) return elt;
      WORD_TYPE v = (table[idx] & (full << ((elt+1) & CACHE)));
      return minimum_element(idx,v,elt);
    }
    */

    inline int prev(const int elt) const {

      WORD_TYPE tab;
      int i = ((elt-1) >> EXP);
      int SHFT = size_word_byte;

      if( i >= neg_words ) {
	int e = ((elt-1) & CACHE), k;
	int j = 1+(e >> 3);

	if( (tab = ((table[i] & (full >> (CACHE - e))) << ((SHFT-j) << 3))) ) 
	  while( j-- ) {
	    if( (k = getlast[(tab & mask_last_char) >> LASTCHAR]) >= 0 )
	      return ( (i<<EXP)+(j<<3)+k );
	    tab = (tab << 8);
	  }
	while( i-- > neg_words ) 
	  if( (tab = table[i]) ) {
	    j = size_word_byte;
	    while( j-- ) {
	      if( (k = getlast[(tab & mask_last_char) >> LASTCHAR]) >= 0 )
		return ( (i<<EXP)+(j<<3)+k );
	      tab = (tab << 8);
	    }
	  }
      }

      return elt;
    }

    inline void xor_to(const Bitset<WORD_TYPE,FLOAT_TYPE>& s) 
    {
      int i = (pos_words > s.pos_words ? s.pos_words : pos_words);
      int j = (neg_words < s.neg_words ? s.neg_words : neg_words);
      while( i-- > j )
	s.table[i] ^= table[i];
    }

    inline void fast_xor_to(const Bitset<WORD_TYPE,FLOAT_TYPE>& s) 
    {
      int i = pos_words;
      while( i-- > neg_words )
	s.table[i] ^= table[i];
    }

    inline void xor_with(const Bitset<WORD_TYPE,FLOAT_TYPE>& s) 
    {
      int i = (pos_words > s.pos_words ? s.pos_words : pos_words);
      int j = (neg_words < s.neg_words ? s.neg_words : neg_words);
      while( i-- > j )
	table[i] ^= s.table[i];
    }

    inline void fast_xor_with(const Bitset<WORD_TYPE,FLOAT_TYPE>& s) 
    {
      int i = pos_words;
      while( i-- > neg_words )
	table[i] ^= s.table[i];
    }

    inline void union_with(const Bitset<WORD_TYPE,FLOAT_TYPE>& s) 
    {
      int i = (pos_words > s.pos_words ? s.pos_words : pos_words);
      int j = (neg_words < s.neg_words ? s.neg_words : neg_words);
      while( i-- > j )
	table[i] |= s.table[i];
    }

    inline void union_with(const int s) 
    {
      if(pos_words>0 && neg_words<=0) table[0] |= s;
    }

    inline void union_to(Bitset<WORD_TYPE,FLOAT_TYPE>& s) const 
    {
      s.union_with( *this );
    }

    inline void intersect_with(const Bitset<WORD_TYPE,FLOAT_TYPE>& s) 
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

    inline void intersect_with(const int s) 
    {
      int i=pos_words;
      while(i-- > 1) table[i] = empt;
      i = 0;
      while(i-- > neg_words) table[i] = empt;
      if(pos_words>0 && neg_words<=0) table[0] &= s;
    }

    inline void intersect_to(Bitset<WORD_TYPE,FLOAT_TYPE>& s) const
    {
      s.intersect_with( *this );
    }

    inline void setminus_with (const Bitset<WORD_TYPE,FLOAT_TYPE>& s) 
    {
      int i = (pos_words > s.pos_words ? s.pos_words : pos_words);
      int j = (neg_words < s.neg_words ? s.neg_words : neg_words);
      while( i-- > j )
	table[i] &= (~(s.table[i]));
    }

    inline void setminus_to (Bitset<WORD_TYPE,FLOAT_TYPE>& s) const
    {
      s.setminus_with( *this );
    }

    inline void xor_with(const Bitset<WORD_TYPE,FLOAT_TYPE>* s) 
    {
      xor_with(*s);
    }

    inline void union_with  (const Bitset<WORD_TYPE,FLOAT_TYPE>* s) 
    {
      union_with(*s);
    }

    inline void intersect_with(const Bitset<WORD_TYPE,FLOAT_TYPE>* s) 
    {
      intersect_with(*s);
    }

    inline void setminus_with (const Bitset<WORD_TYPE,FLOAT_TYPE>* s) 
    {
      setminus_with(*s);
    }
    inline void union_to  (Bitset<WORD_TYPE,FLOAT_TYPE>* s) const
    {
      s->union_with(*this);
    }

    inline void intersect_to(Bitset<WORD_TYPE,FLOAT_TYPE>* s) const
    {
      s->intersect_with(*this);
    }

    inline void setminus_to (Bitset<WORD_TYPE,FLOAT_TYPE>* s) const
    {
      s->setminus_with(*this);
    }

    inline bool equal(const Bitset<WORD_TYPE,FLOAT_TYPE>& s) const 
    {
      int i=pos_words;
      int j=s.pos_words;
      int k;
  
      while( j > i )
	if( s.table[--j] ) return false; 
      while( i > j )
	if( table[--i] ) return false;

      j=neg_words;
      k=s.neg_words;

      while( j > k )
	if( s.table[k++] ) return false;
      while( k > j )
	if( table[j++] ) return false;

      while( i-- > j )
	if( table[i] != s.table[i] ) return false;
  
      return true;
    }
    
    inline bool includes(const WORD_TYPE s) const 
    {
      return( pos_words && neg_words<1 && (table[0] & s) == s );
    }

    inline bool included(const WORD_TYPE s) const 
    {
      bool inc = true;
      int k = pos_words;
      if(neg_words>0 || pos_words<1) {
	while(k>neg_words && inc) inc = !(table[--k]);
      } else {
	while(k>1 && inc) inc = !(table[--k]);
	inc = ((table[--k] & s) == table[0]);
	while(k>neg_words && inc) inc = !(table[--k]);
      }
      return inc;
    }

    inline bool included(const Bitset<WORD_TYPE,FLOAT_TYPE>& s) const 
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

    inline bool includes(const Bitset<WORD_TYPE,FLOAT_TYPE>& s) const 
    {
      int i = (pos_words > s.pos_words ? s.pos_words : pos_words);
      int j = (neg_words < s.neg_words ? s.neg_words : neg_words);
      int k = s.pos_words;
      while( k > i ) {
	--k;
	if( s.table[k] ) return false;
      }
      while( k > j ) {
	--k;
	if( s.table[k] != (table[k] & s.table[k]) ) {
	  return false;
	}
      }
      while( k-- > s.neg_words ) {
	if( s.table[k] ) return false;
      }

      return true;
    }

    inline bool included(const Bitset<WORD_TYPE,FLOAT_TYPE>* s) const 
    {
      return included( *s );
    }

    inline bool includes(const Bitset<WORD_TYPE,FLOAT_TYPE>* s) const 
    {
      return includes( *s );
    }

    inline bool intersect(const Bitset<WORD_TYPE,FLOAT_TYPE>* s) const 
    {
      return intersect( *s );
    }

    inline bool included(const int lb, const int ub) const 
    {
      int neg_int = lb >> EXP;
      int pos_int = ub >> EXP;
      int k = pos_words;
      while( k > pos_int ) 
	if( table[--k] ) return false;
      k = neg_words;
      while( k < neg_int )
	if( table[k++] ) return false;
      if(neg_int == pos_int) {
	k = ((full << (lb & CACHE)) & (full >> (CACHE - (ub & CACHE))));
	return (k & table[neg_int]) == table[neg_int];
      } else {
	return ((((full << (lb & CACHE)) & table[neg_int]) == table[neg_int]) &&
		(((full >> (CACHE - (ub & CACHE))) & table[pos_int]) == table[pos_int]));
      }
    }

    inline bool includes(const int lb, const int ub) const
    {
      int neg_int = lb >> EXP;
      int pos_int = ub >> EXP;
      int k = pos_int-1;
      unsigned int u, l;
      while( k > neg_int ) {
	if( table[k] != full ) return false;
	--k;
      }
      if(neg_int == pos_int) {
	u = ((full << (lb & CACHE)) & (full >> (CACHE - (ub & CACHE))));
	return (u & table[neg_int]) == u;
      } else {
	u = (full >> (CACHE - (ub & CACHE)));
	l = (full << (lb & CACHE));
	return (((l & table[neg_int]) == l) &&
		((u & table[pos_int]) == u));
      }
    }



    /*!
     * Returns the number of bits set in v.
     * For a derivation of this algorithm, see
     * "Algorithms and data structures with applications to 
     *  graphics and geometry", by Jurg Nievergelt and Klaus Hinrichs,
     *  Prentice Hall, 1993.
     */
    inline unsigned int word_size(WORD_TYPE v) const 
    {
      return __builtin_popcount(v);
      /*
      v = v - ((v >> 1) & (WORD_TYPE)~(WORD_TYPE)0/3);                           // temp
      v = (v & (WORD_TYPE)~(WORD_TYPE)0/15*3) + ((v >> 2) & (WORD_TYPE)~(WORD_TYPE)0/15*3);      // temp
      v = (v + (v >> 4)) & (WORD_TYPE)~(WORD_TYPE)0/255*15;                      // temp
      return (WORD_TYPE)(v * ((WORD_TYPE)~(WORD_TYPE)0/255)) >> (sizeof(v) - 1) * CHAR_BIT; // count
      */
    }
	
    inline unsigned int size() const 
    {  
      int i=pos_words;
      unsigned int c=0;
      WORD_TYPE v;
      while( i-- > neg_words ) 
	if( (v = table[i]) ) 
	  c += word_size(v);
      return c;  
    }

    inline unsigned int word_size() const 
    {  
      unsigned int v, c=0;
      if( (v = table[neg_words]) ) 
	c = word_size(v);
      return c;  
    }

    inline unsigned int size( const int i ) const
    {  
      WORD_TYPE v;
      unsigned int c=0;
      if( (v = table[i]) ) 
	c = word_size(v);
      return c;  
    }



    /*!
      Check if element elt belong to the set [O(1)]
    */
    inline  bool contain(const int elt)const 
    {
      int i = (elt >> EXP);
      return ( (i >= neg_words) && 
	       (i <  pos_words) && 
	       (table[i] & ((WORD_TYPE)1 << (elt & CACHE))) );
    }

    inline  bool fast_contain(const int elt)const 
    {
      return ( (table[(elt >> EXP)] & ((WORD_TYPE)1 << (elt & CACHE))) );
    }

    inline  bool word_contain(const int elt)const 
    {
      return ( (table[neg_words] & ((WORD_TYPE)1 << (elt & CACHE))) );
    }
    /*!
      Add element elt into the set [O(1)]
    */

    inline  void add(const int elt)
    {
      int i = (elt >> EXP);
      if( (i >= neg_words) && 
	  (i <  pos_words) ) 
	table[i] |= ((WORD_TYPE)1 << (elt & CACHE));
    }

    inline  void fast_add(const int elt)
    {
      table[(elt >> EXP)] |= ((WORD_TYPE)1 << (elt & CACHE));
    }

    inline   void word_add(const int elt)
    {
      table[neg_words] |= ((WORD_TYPE)1 << (elt & CACHE));
    }

    /*!
      Add element elt into the set or remove it if it is already contain [O(1)]
    */
    inline  void invert(const int elt)
    {
      int i = (elt >> EXP);
      if( (i >= neg_words) && 
	  (i <  pos_words) )
	table[i] ^= ((WORD_TYPE)1 << (elt & CACHE));
    }

    inline  void fast_invert(const int elt)
    {
      table[(elt >> EXP)] ^= ((WORD_TYPE)1 << (elt & CACHE));
    }

    inline  void word_invert(const int elt)
    {
      table[neg_words] ^= ((WORD_TYPE)1 << (elt & CACHE));
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
    inline bool intersect(const Bitset<WORD_TYPE,FLOAT_TYPE>& s)const
    {
      int i = (pos_words > s.pos_words ? s.pos_words : pos_words);
      int j = (neg_words < s.neg_words ? s.neg_words : neg_words);
      while( i-- > j )
	if(table[i] & s.table[i]) return true;
      return false;
    }

    /*!
      Return true iff the calling object intersect s (s is assumed to be a bitset in {0,..,31}) [O(N/32)]
    */
    inline bool intersect(const int s) const
    {
      return(pos_words && neg_words<1 && (table[0]&s));
    }

    inline bool word_intersect(const Bitset<WORD_TYPE,FLOAT_TYPE>& s) const 
    {
      return ( table[neg_words] & s.table[neg_words] ) ;
    }

    inline  bool fast_intersect(const Bitset<WORD_TYPE,FLOAT_TYPE>& s, int& idx)const
    {
      if( table[idx] & s.table[idx] ) return true; 
      if( pos_words > neg_words ) {
	idx = pos_words;
	while( idx > neg_words ) {
	  --idx;
	  if(table[idx] & s.table[idx]) return true;
	}
      }
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

      WORD_TYPE masked_lb = (full << (lb & CACHE));
      WORD_TYPE masked_ub = (full >> (CACHE - (ub & CACHE)));

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
      of the set minus x is removed [O(N/32)]
    */
    inline void increment(const int v)
    {
      int step = (v >> EXP); 
      int i = pos_words;
      int e = (v & CACHE);
      int f = size_word_bit-e;
      int j = neg_words+step;
      WORD_TYPE mask = ((WORD_TYPE)~0 << f);
      while( --i > j ) 
	table[i] = ((table[i-step] << e) | ((table[i-step-1] & mask) >> f));
      if( i >= neg_words+step ) table[i] = (table[i-step] << e);
      while( i > neg_words ) table[--i] = 0;
    }

    /*!
      Decrement by x all elements in the set.
      Any element lower than x is removed [O(N/32)]
    */
    inline  void decrement(const int v)
    {
      int step = (v >> EXP); 
      int i = neg_words-1;
      int e = (v & CACHE);
      int f = size_word_bit-e;
      int j = pos_words-step-1;
      WORD_TYPE mask = ((WORD_TYPE)~0 >> e);
      while( ++i < j )
	table[i] = ((table[i+step] >> e) | ((table[i+step+1] & mask) << f));
      if( i < pos_words-step ) table[i] = (table[i+step] >> e);
      while( ++i < pos_words ) table[i] = 0;
    }

    /*!
      Changes every value to its arythmetic negation
    */
    inline void negate( Bitset<WORD_TYPE,FLOAT_TYPE>& s ) const
    {
      int i = (pos_words > -s.neg_words ? -s.neg_words : pos_words);
      int j = (neg_words < -s.pos_words ? -s.pos_words : neg_words);

      unsigned int a; 
      WORD_TYPE mask, v, aux, rest = ( i < pos_words && (table[i] & 1) );

      while( i-- > j ) {
	aux = (table[i] & 1);
	v = (table[i] >> 1);
	mask = ~0;         
	a = sizeof(v) * CHAR_BIT; // bit size; must be power of 2 
	while ((a >>= 1) > 0) 
	  {
	    mask ^= (mask << a);
	    v = ((v >> a) & mask) | ((v << a) & ~mask);
	  }	
	s.table[-i-1] = (v | rest);
	rest = aux;
      }
      if(rest)
	s.table[i+1] |= rest;
    }

    /*!
      Add all elements between 0 to capacity [O(N/32)]
    */
    inline  void fill()
    {
      int i = pos_words;
      while( i > neg_words )
	table[--i] = full;
    }

    inline  void fill(const int lb, const int ub)
    {
      int i = (ub >> EXP);
      int j = (lb >> EXP);
      
      if( i >= neg_words || j < pos_words ) {
		
	WORD_TYPE masked_lb = (full << (lb & CACHE));
	WORD_TYPE masked_ub = (full >> (CACHE - (ub & CACHE)));
	
	if( i == j ) {

	  table[i] |= (masked_lb & masked_ub);
	  
	} else {
	  
	  if( i >= pos_words ) {
	    i = pos_words-1;
	  } else {
	    table[i--] |= masked_ub;
	  }
	  
	  if( j < neg_words ) {
	    j = neg_words;
	  } else {
	    table[j++] |= masked_lb;
	  }
	  
	  while( i >= j )
	    table[i--] |= full;
	}
      }
    }

    /*!
      Remove all elements [O(N/32)]
    */
    inline  void clear()
    {
      int i = pos_words;
      while( i > neg_words )
	table[--i] = empt;
    }

    /*!
      Remove all elements but v [O(N/32)]
    */
    inline  void set_to( const int v )
    {
      int i, j = (v >> EXP);
      for(i=neg_words; i<j; ++i)
	table[i] = empt;
      table[j] = ((WORD_TYPE)1 << v);
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
      Remove all elements strictly lower than l [O(N/32)]
    */
    inline void set_min(const int bound)
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
      Remove all elements strictly greater than u [O(N/32)]
    */
    inline void set_max(const int bound)
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
      Remove all elements in the interval [l..u] [O(N/32)]
    */
    inline  void remove_interval(const int lb, const int ub)
    {
      if( lb <= ub ) {
	int lb_word = lb >> EXP;
	int ub_word = ub >> EXP;

	lb_word = ( lb_word < neg_words ? neg_words : lb_word );
	ub_word = ( ub_word >= pos_words ? pos_words-1 : ub_word );

	WORD_TYPE masked_lb = 0;
	WORD_TYPE masked_ub = 0;

	if( lb_word >= neg_words ) 
	  // add a '0' on the 32nd bit, because >> 32 does nothing
	  masked_lb = ((full/2) >> (CACHE - (lb & CACHE)));
	if( ub_word < pos_words ) 
	  masked_ub = ((full-1) << (ub & CACHE));

	if( lb_word == ub_word ) {
	  table[lb_word] &= (masked_lb | masked_ub);
	} else {
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
    inline void add_interval(int lb, int ub)
    {

      if( lb <= ub ) {
	int lb_word = lb >> EXP;
	int ub_word = ub >> EXP;
    
	lb_word = ( lb_word < neg_words ? neg_words : lb_word );
	ub_word = ( ub_word >= pos_words ? pos_words-1 : ub_word );

	WORD_TYPE masked_lb = full;
	WORD_TYPE masked_ub = full;
	if( lb_word >= neg_words ) 
	  //masked_lb ^= (full >> (CACHE - (lb & CACHE) + 1));
	  masked_lb ^= ((full/2) >> (CACHE - (lb & CACHE)));
	if( ub_word < pos_words ) 
	  //masked_ub ^= ((full-1) << (ub & CACHE) );
	  masked_ub ^= ((full-1) << (ub & CACHE));

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
      return fast_contain(i);
    }

    std::ostream& display(std::ostream& os) const {
      os << "{";
      if( !empty() ) {
	int last = NOVAL, cur=min(), aft;

	bool flag=false;
	do{
	  aft = next(cur);

	  if(aft != cur+1 || cur != last+1) {
	    if( flag )
	      os << ",";
	    os << (int)cur;
	    flag = true;
	  } else if(flag) {
	    os << "..";
	    flag = false;
	  }
	  last = cur;
	  cur = aft;
	} while( cur != NOVAL && cur != last );
      }
      os << "}";
      return os;
    }
	
    std::string to_str() const {
			    std::ostringstream oss;
				oss << "{";
				//std::string rstr = std::string("{");
      if( !empty() ) {
	int last = NOVAL, cur=min(), aft;

	bool flag=false;
	do{
	  aft = next(cur);

	  if(aft != cur+1 || cur != last+1) {
	    if( flag )
	      //rstr += std::string(",");
			oss << ",";

  	    oss << (int)cur;

	    flag = true;
	  } else if(flag) {
	    //rstr += std::string("..");
		  oss << "..";
	    flag = false;
	  }
	  last = cur;
	  cur = aft;
	} while( cur != NOVAL && cur != last );
      }
      //rstr += std::string("}");
	  oss << "}";
      return oss.str();
    }

    void  print_bits(std::ostream & os) const 
    {
      os << "[";
      for(int i=neg_words; i<pos_words; ++i)
	showUint( table[i], os );
      os << "]";
    }

  };


  /**********************************************
   * InlinedBitSet
   **********************************************/
  /// The data array is stored directly in the object

  template< class WORD_TYPE, class FLOAT_TYPE >
  class InlinedBitset : public Bitset<WORD_TYPE,FLOAT_TYPE> 
  {

  public:

    /**@name Parameters*/
    //@{    
    WORD_TYPE data[0];
    //@}
    
    static InlinedBitset<WORD_TYPE,FLOAT_TYPE>* make_new(const int lb, const int ub, const WORD_TYPE p)
    {
      int n_w = (lb >> Bitset<WORD_TYPE,FLOAT_TYPE>::EXP);
      int p_w = (ub >> Bitset<WORD_TYPE,FLOAT_TYPE>::EXP)+1;
      int size = p_w - n_w;
      
      void* mem = malloc(sizeof(InlinedBitset<WORD_TYPE,FLOAT_TYPE>) + sizeof(WORD_TYPE)*(size));
      return new (mem) InlinedBitset<WORD_TYPE,FLOAT_TYPE>(lb, ub, p);
    }

    InlinedBitset(const int lb, const int ub, const WORD_TYPE p)
    {
      initialise(lb,ub,p,&(data[0]));
    }
    
  };

  typedef Bitset< unsigned long long int, double > Bitset64;
  typedef Bitset< unsigned int, float > Bitset32;
  typedef InlinedBitset< unsigned long long int, double > iBitset64;
  typedef InlinedBitset< unsigned int, float > iBitset32;

#ifdef _BIT64
  
  typedef Bitset64 BitSet;

#else

  typedef Bitset32 BitSet;

#endif


  /**********************************************
   * Queue
   **********************************************/
  /// List of constraints for computing the GAC closure (fifo).
  
  class Queue {
  public:
    
    /**@name Parameters*/
    //@{
    int *next; // next element in the list (one for each elt + one for the head)
    int _head; // last_index+1
    int _tail; // last element
    int offset; // first index
    //@}

    /**@name Constructors*/
    //@{
    Queue() {
      next = NULL;
      _head = NOVAL;
      _tail = NOVAL;
      offset = 0;
    }
    ~Queue()
    {
      next += offset;
      delete [] next;
    }
    void cancel()
    {
      next = NULL;
      offset = 0;
    }
    bool is_initialised() 
    {
      return (next != NULL);
    }
    void initialise(const int n)
    {
      next = new int[n+1];
      std::fill( next, next+n, NOVAL );
      _head = _tail = n;
      next[n] = n;
      offset = 0;
    }
    void initialise(const int lb, const int ub)
    {
      next = new int[ub-lb+2];
      next-=lb;
      std::fill( next+lb, next+ub+1, NOVAL );
      _head = _tail = ub+1;
      next[ub+1] = ub+1;
      offset = lb;
    }
    void extend(const int elt)
    {
      int new_lb = (elt < offset ? elt : offset);
      int new_ub = (elt >= _head ? elt : _head-1);

      if(new_lb < offset || new_ub >= _head) {

	// 	std::cout << "EXTEND QUEUE" << std::endl;
	
	int *old_next = next;
	int old_head = _head;
	int old_tail = _tail;
	int old_offset = offset;

	// 	std::cout << "before: "
	// 		  << _head << " " 
	// 		  << _tail << " "
	// 		  << next[_head] 
	// 		  << std::endl;

	initialise(new_lb, new_ub);

	// 	std::cout << "afteri: "
	// 		  << _head << " " 
	// 		  << _tail << " "
	// 		  << next[_head] 
	// 		  << std::endl;


	// 	std::cout << " saved: "
	// 		  << old_head << " " 
	// 		  << old_tail << " "
	// 		  << old_next[old_head] 
	// 		  << std::endl;
	
	if(old_tail != old_head) {
	  int elt = old_next[old_head];
	  while(elt != old_head) {
	    //std::cout << "ADD " << elt << std::endl;
	    add(elt);
	    elt = old_next[elt];
	  }
	}

	// 	std::cout << " after: "
	// 		  << _head << " " 
	// 		  << _tail << " "
	// 		  << next[_head] 
	// 		  << std::endl;


	// 	int *aux = next;
	// 	next = new int[new_ub-new_lb+2];
	// 	std::fill(next, next+new_ub-new_lb+2, NOVAL);
	// 	next-=new_lb;
	// 	//memcpy(next+offset, aux+offset, (_head-offset+1)*sizeof(int));

	
	

	// // 	for(int i=new_lb; i<=new_ub+1; ++i)
	// // 	  if(next[i] != NOVAL)
	// // 	    std::cout << std::setw(3) << next[i] ;
	// // 	  else 
	// // 	    std::cout << " . ";
	// // 	std::cout << std::endl;
	// // 	for(int i=new_lb; i<=new_ub+1; ++i)	  
	// // 	  if(i>=offset && i<=_head && aux[i] != NOVAL)
	// // 	    std::cout << std::setw(3) << next[i] ;
	// // 	  else 
	// // 	    std::cout << " . ";
	// // 	std::cout << std::endl;


	// 	for(int i=new_lb; i<offset; ++i)
	// 	  next[i] = NOVAL;

	// 	for(int i=_head+1; i<new_ub+1; ++i)
	// 	  next[i] = NOVAL;

	// 	if(_tail == _head) {
	// 	  _tail = new_ub+1;
	// 	}
	
	// 	next[_tail] = new_ub+1;
	// 	next[new_ub+1] = next[_head];
	

	// 	for(int i=new_lb; i<=new_ub+1; ++i)
	// 	  if(next[i] != NOVAL)
	// 	    std::cout << std::setw(3) << next[i] ;
	// 	  else 
	// 	    std::cout << " . ";
	// 	std::cout << std::endl;
	// 	for(int i=new_lb; i<=new_ub+1; ++i)	  
	// 	  if(i>=offset && i<=_head && aux[i] != NOVAL)
	// 	    std::cout << std::setw(3) << next[i] ;
	// 	  else 
	// 	    std::cout << " . ";
	// 	std::cout << std::endl;

	// 	//exit(1);

	
	// // 	if(_head <= new_ub) {
	// // 	  next[new_ub+1] = aux[_head];
	// // 	  next[_head] = NOVAL;
	// // 	}

	// 	for(int i=new_lb; i<=new_ub+1; ++i)
	// 	  if(next[i] != NOVAL)
	// 	    std::cout << std::setw(3) << next[i] ;
	// 	  else 
	// 	    std::cout << " . ";
	// 	std::cout << std::endl;
	// 	for(int i=new_lb; i<=new_ub+1; ++i)	  
	// 	  if(i>=offset && i<=_head && aux[i] != NOVAL)
	// 	    std::cout << std::setw(3) << aux[i] ;
	// 	  else 
	// 	    std::cout << " . ";
	// 	std::cout << std::endl;

	// // 	//exit(1);

	// 	aux += offset;
	// 	delete [] aux;

	// 	offset = new_lb;
	// 	_head = new_ub+1;

	old_next += old_offset;
	delete [] old_next;

	//	display(std::cout);
      }
    }
    void declare(const int elt) {
      extend(elt);
      add(elt);
    }
    //@}
  
    /*!@name Accessors*/
    //@{
    /// Add an element to the queue
    inline int first() {
      return next[_head];
    }

    inline void add(const int elt)
    {
      next[_tail] = elt;
      _tail = elt;
      next[elt] = _head;
    }

    /// Pop the first event of the queue
    inline int pop()
    {
      int elt=next[_head];
      next[_head] = next[elt];
      if( next[_head] == _head ) _tail = _head;

      return elt;
    }

    /// Is the queue empty?
    inline bool empty() const
    {
      return _head == _tail;
    }

    /// Clear all elements of the queue
    inline void clear()
    {
      _tail = _head;
    }
    //@}
    
    Queue& operator=(const Queue& q) {
      _head = q._head;
      _tail = q._tail;
      next = q.next;
      offset = q.offset;
      return *this;
    }

    /*!@name Miscellanous*/
    //@{ 
    std::ostream& display(std::ostream& os) const {

      os << _head << " " << _tail << " " << offset << ": ";
      if(_head != NOVAL){
	for(int i=offset; i<=_head; ++i) {
	  if(next[i] != NOVAL)
	    os << next[i] << " ";
	  else 
	    os << ". ";
	}
      }
      os << "\n" ;

      if(_head != NOVAL){
	int i;
	os << _head << "(";
	os.flush();
	i = next[_head];
	if(i != _head) {
	  os << i;
	  i = next[i];
	  while( i != _head ) {
	    os << " " << i;
	    
	    if(i == NOVAL) break;
	    
	    i = next[i];
	  }
	}
	os << ")";
      }


      return os;
    }
    //@}

  };



  template <class T>
  class BinaryMinHeap {

  private:
    
    inline void sift_up(unsigned int index, const T x) {
      ++num_operations;
      
      unsigned int ascendant = (index-1)/2;    
      while(index && x < data[ascendant]) {
	++num_operations;
	
	data[index] = data[ascendant];
	data[ascendant] = x;
	index = ascendant;
	ascendant /= 2;
      }
    }
    
    inline void sift_down(unsigned int index, const T x) {
      unsigned int descendant = index*2+1;
      unsigned int smallest = index;
      ++num_operations;
      
      while(descendant < data.size) {
	++num_operations;
	
	if(data[descendant] < data[smallest]) {
	  smallest = descendant;
	}
	++descendant;
	if(descendant < data.size && data[descendant] < data[smallest]) {
	  smallest = descendant;
	}
	if(smallest == index) break;
	
	data[index] = data[smallest];
	data[smallest] = x;
	
	index = smallest;
	descendant = smallest*2+1;
      }
    }
    
  public: 
    
    int num_operations;
    Vector< T > data;
    
    BinaryMinHeap() { num_operations=0; }
    virtual ~BinaryMinHeap() { }
    
    void add(const T x) {
      data.add(x);
      sift_up(data.size-1, x);
    }
    
    T pop_min() {
      T the_min = data[0];
      data[0] = data[--data.size];
      data[data.size] = the_min;
      sift_down(0, data[0]);
      return the_min;
    }
    
    /*!@name Printing*/
    //@{
    std::ostream& display(std::ostream& os) const {
      

      unsigned int n_level = 0;
      while((unsigned int)(1<<n_level)<=data.size) ++n_level;
      
      unsigned int offset = 0;
      unsigned int step = 1;

      unsigned int i, j, l, first, last;

      for(i=0; i<30; ++i) {
	os << " " << data[i];
      }
      os << std::endl;
	
      for(l=0; l<n_level; ++l) {
	first = (1<<(n_level-l-1))-1;
	last = (1<<(n_level-l))-1;
	if(last>data.size) last = data.size;
	for(j=0; j<offset; ++j) os << "   ";
	for(i=first; i<last; ++i) {
	  os << std::setw(3) << data[i];
	  if(i<last-1) for(j=0; j<step; ++j) os << "   ";
	}
	os << std::endl;
	offset = offset+(1+step)/2;
	step = 2*step+1;
      }
      return os;
    }
    //@}
    
  };





  template <class T>
  class IndexedBinaryMaxHeap {

  private:
    
    // starting at position 'rank' in the tree, move up to the root swapping elements until this branch is sorted
    inline void sift_up(unsigned int rank) {
      unsigned int ascendant = (rank-1)/2;
      int x = heap[rank];
      while(rank && value[x] > value[heap[ascendant]]) {
	index[x] = ascendant;
	index[heap[ascendant]] = rank;
	
	heap[rank] = heap[ascendant];
	heap[ascendant] = x;

	rank = ascendant--;
	ascendant /= 2;
      }
    }
    
    // starting at position 'rank' in the tree, move down to a leaf swapping elements until this branch is sorted
    inline void sift_down(unsigned int rank) {
      unsigned int descendant = rank*2+1;
      unsigned int largest = rank;
      int x = heap[rank];
      while(descendant < heap.size) {
	if(value[heap[descendant]] > value[heap[largest]]) {
	  largest = descendant;
	}
	++descendant;
	if(descendant < heap.size && value[heap[descendant]] > value[heap[largest]]) {
	  largest = descendant;
	}
	if(largest == rank) break;
	
	index[heap[largest]] = rank;
	index[x] = largest;

	heap[rank] = heap[largest];
	heap[largest] = x;
	
	rank = largest;
	descendant = largest*2+1;
      }
    }
    
  public: 

    // the value of each element (indexed by element names)
    Vector< T >   value;

    // the data structure of the heap (populated with element names)
    Vector< int >  heap;

    // the current position of each element in the heap (indexed by element names)
    Vector< int > index;

    
    IndexedBinaryMaxHeap() {}
    virtual ~IndexedBinaryMaxHeap() { }
    
    int add(const T x) {
      int i = value.size;
      value.add(x);
      index.add(heap.size);
      heap.add(i);
      sift_up(heap.size-1);
      
      return i;
    }

    void change(const int i, const T new_value) {
      if(new_value > value[i]) {
	value[i] = new_value;
	sift_up(index[i]);
      } else {
	value[i] = new_value;
	sift_down(index[i]);
      }
    }

    T pop_max() {
      int the_root = heap[0];

      heap[0] = heap[--heap.size];
      heap[heap.size] = the_root;

      index[the_root] = heap.size;
      index[heap[0]] = 0;

      sift_down(0);
      return value[the_root];
    }



    
    /*!@name Printing*/
    //@{
    std::ostream& display(std::ostream& os) const {
      

      unsigned int n_level = 0;
      while((unsigned int)(1<<n_level)<=heap.size) ++n_level;
      
      unsigned int offset = 0;
      unsigned int step = 1;

      unsigned int i, j, l, first, last;

   
	
      for(l=0; l<n_level; ++l) {
	first = (1<<(n_level-l-1))-1;
	last = (1<<(n_level-l))-1;
	if(last>heap.size) last = heap.size;
	for(j=0; j<offset; ++j) os << "   ";
	for(i=first; i<last; ++i) {
	  if(index[heap[i]] != i) {
	    os << "!!!!" << std::endl;
	    exit(1);
	  }
	  os << std::setw(3) << value[heap[i]] ;
	  if(i<last-1) for(j=0; j<step; ++j) os << "   ";
	}
	os << std::endl;
	offset = offset+(1+step)/2;
	step = 2*step+1;
      }
      return os;
    }
    //@}
    
  };


  int __modulo_fct__(const int x, const int m) ;

  //class Variable;
  class BiInterval;
  class Interval {

  public:
    
    int min;
    int max;
    
    Interval();
    Interval(const BiInterval b);
    Interval(const int _min, const int _max);
    virtual ~Interval();
    
    Interval get_union(Interval arg);
    
    bool contain(const int x) const;
    bool empty() const;
    
    Interval operator*(const Interval);
    Interval operator/(const Interval);
    Interval anti_mul(const Interval);
    Interval operator-() const;
    Interval operator%(const int);
    Interval positive_modulo(const int);
    Interval operator_modulo(const int);

    Interval target_positive_modulo(const int, const Interval);
    Interval target_modulo(const int, const Interval);
    Interval target_c_modulo(const int, const Interval);

    void operator+=(const int x);
    void operator-=(const int x);

    std::ostream& display(std::ostream& os) const;
    
  };




  class Variable;
  class PositiveHalfDomain;
  class NegativeHalfDomain : public Interval {

public:
  
  NegativeHalfDomain();
  NegativeHalfDomain(const int _min, const int _max);
  virtual ~NegativeHalfDomain();


  Interval operator*(const PositiveHalfDomain arg);
  Interval operator*(const NegativeHalfDomain arg);
  Interval operator*(const int arg);

  Interval anti_mul(const PositiveHalfDomain arg);
  Interval anti_mul(const NegativeHalfDomain arg);
  Interval anti_mul(const int arg);

  Interval operator/(const PositiveHalfDomain arg);
  Interval operator/(const NegativeHalfDomain arg);
  Interval operator/(const int arg);
  Interval divided_by(const PositiveHalfDomain arg, const Variable target);
  Interval divided_by(const NegativeHalfDomain arg, const Variable target);

  Interval anti_div_X(const PositiveHalfDomain arg);
  Interval anti_div_X(const NegativeHalfDomain arg);
  Interval anti_div_X_pos(const int arg);
  Interval anti_div_X_neg(const int arg);

  Interval anti_div_Y(const PositiveHalfDomain arg);
  Interval anti_div_Y(const NegativeHalfDomain arg);
  Interval anti_div_Y_pos(const int arg);
  Interval anti_div_Y_neg(const int arg);


  // Interval operator%(const int mod);
  // Interval operator%(const int mod);

  // // the return value is the interval I such that I%mod = this
  // Interval anti_modulo(const int mod);
  // Interval anti_modulo(const int mod);

  //void operator=(const Interval I) {min = I.min; max = I.max;}
    //Interval operator-() { return Interval(-max, -min); }
  
};


class PositiveHalfDomain : public Interval {

public:
  
  PositiveHalfDomain();
  PositiveHalfDomain(const int _min, const int _max);
  virtual ~PositiveHalfDomain();

  Interval operator*(const PositiveHalfDomain arg);
  Interval operator*(const NegativeHalfDomain arg);
  Interval operator*(const int arg);

  Interval anti_mul(const PositiveHalfDomain arg);
  Interval anti_mul(const NegativeHalfDomain arg);
  Interval anti_mul(const int arg);

  Interval operator/(const PositiveHalfDomain arg);
  Interval operator/(const NegativeHalfDomain arg);
  Interval divided_by(const PositiveHalfDomain arg, const Variable target, const bool pos=true);
  Interval divided_by(const NegativeHalfDomain arg, const Variable target);
  Interval operator/(const int arg);

  Interval anti_div_X(const PositiveHalfDomain arg);
  Interval anti_div_X(const NegativeHalfDomain arg);
  Interval anti_div_X_pos(const int arg);
  Interval anti_div_X_neg(const int arg);

  Interval anti_div_Y(const PositiveHalfDomain arg);
  Interval anti_div_Y(const NegativeHalfDomain arg);
  Interval anti_div_Y_pos(const int arg);
  Interval anti_div_Y_neg(const int arg);
  // Interval operator%(const int mod);
  // Interval operator%(const int mod);

  // // the return value is the interval I such that I%mod = this
  // Interval anti_modulo(const int mod);
  // Interval anti_modulo(const int mod);


  //void operator=(const Interval I) {min = I.min; max = I.max;}
  //Interval operator-() { return Interval(-max, -min); }
};


  class Variable;
  class BiInterval {

  public:
    PositiveHalfDomain positive;
    NegativeHalfDomain negative;
    bool zero;

    BiInterval();
    BiInterval(const Variable x);
    BiInterval(const int _min, const int _max);
    BiInterval(const Interval I);
    BiInterval(const Interval neg, const Interval pos, const bool z);
    BiInterval(const Interval neg, const Interval pos, const Interval z);
    BiInterval(const int n_min, const int n_max, const int p_min, const int p_max, const bool z);
    void initialise(const int _min, const int _max);
    virtual ~BiInterval();
    

    int get_min() const;
    int get_max() const;

    int get_min_abs() const;
    int get_max_abs() const;

    int is_hollow() const;
    //Interval get_hole() const;

    BiInterval operator*(const BiInterval arg);
    BiInterval anti_mul(const BiInterval arg);
    BiInterval operator/(const BiInterval arg);
    BiInterval divided_by(const BiInterval arg, const Variable target);
    BiInterval anti_div_X(const BiInterval arg);
    BiInterval anti_div_Y(const BiInterval arg);

    BiInterval operator*(const int arg);
    BiInterval operator/(const int arg);
    BiInterval anti_div_X(const int arg);
    BiInterval anti_div_Y(const int arg);
    // BiInterval operator%(const int mod);
    // BiInterval anti_modulo(const int mod);

    bool operator==(const int x) const;
    void operator=(const int x);

    std::ostream& display(std::ostream& os) const;
    
  };


  class IntervalList : public Vector<Interval> {
    
  public:
    IntervalList();
    virtual ~IntervalList();

    void union_with(const IntervalList& with, IntervalList& into) const;
    void intersect_with(const IntervalList& with, IntervalList& into) const;

    void operator=(const IntervalList& l);

    void push(const int lb, const int ub);
    void push(const Interval& I);

  };

 


  // /**********************************************
  //  * VariableQueue
  //  **********************************************/
  // /// List of integer, used because we can do fifo push/pop
  // class VariableQueue {
  // public:

  //   /**@name Constructors*/
  //   //@{
  //   // index of the following element 
  //   int *next;
  //   // type of event
  //   int *trigger;
  //   // index of the first element
  //   int head;
  //   // last element
  //   int tail;
  //   // size of the structure
  //   int capacity;
  //   //@}

  //   /**@name Constructors*/
  //   //@{
  //   VariableQueue() {
  //     trigger = NULL;
  //     next = NULL;
  //     capacity = 0;
  //   }
  //   ~VariableQueue()
  //   {
  //     delete [] trigger;
  //     delete [] next;
  //   }
  //   void initialise(int n)
  //   {
  //     //      X = x;
  //     trigger = new int[n];
  //     std::fill( trigger, trigger+n, 0 );
  //     next = new int[n+2];
  //     std::fill( next, next+n, NOVAL );
  //     head = tail = n;
  //     next[n] = n;
  //     capacity = n;
  //   }

  //   void declare(int x) {
  //     if(x>=capacity) {
  // 	int new_capacity = 2*capacity;

  // 	int *new_trigger = new int[new_capacity];
  // 	memcpy(new_trigger, trigger, capacity*sizeof(int));
  // 	std::fill( new_trigger+capacity, new_trigger+new_capacity, 0 );

  // 	int *new_next = new int[new_capacity+2];
  // 	memcpy(new_next, next, capacity*sizeof(int));
  // 	std::fill( new_next+capacity, new_next+new_capacity, 0 );
  // 	std::fill( next+capacity, next+new_capacity, NOVAL );

  // 	if(head == tail) {
  // 	  head = tail = new_capacity;
  // 	  next[new_capacity] = new_capacity;
  // 	}

  // 	delete [] trigger;
  // 	delete [] next;
	
  // 	trigger = new_trigger;
  // 	next = new_trigger;
  // 	capacity = new_capacity;
  //     }
  //   }
  //   //@}
  
  //   /*!@name Accessors*/
  //   //@{
  //   /// Add an element to the queue
  //   inline void add(const int i, const int event)
  //   {
  //     //int i=v->id;
  //     if( !trigger[i] ) {
  // 	next[tail] = i;
  // 	tail = i;
  // 	next[i] = head;
  //     }
  //     trigger[i] |= event;
  //   }

  //   /// Pop the first event of the queue
  //   inline int pop( int& event )
  //   {
  //     int i=next[head];
  //     next[head] = next[i];
  //     event = trigger[i];
  //     trigger[i] = 0;

  //     if( next[head] == head )
  // 	tail = head;
  //     return i;
  //   }

  //   /// Is the queue empty?
  //   inline bool empty() const
  //   {
  //     return (next[head] == head);
  //   }

  //   /// Clear all elements of the queue
  //   inline void clear()
  //   {
  //     while( next[head] != head ) {
  // 	trigger[next[head]] = 0;
  // 	next[head] = next[next[head]];
  //     }
  //     tail = head;
  //   }
  //   //@}

  //   /*!@name Miscellanous*/
  //   //@{ 
  //   std::ostream& display(std::ostream& os) const
  //   {
  //     int i;
  //     os << "{";
  //     i = next[head];
  //     while( i != head ) {
  // 	os << i << " ";
  // 	i = next[i];
  //     }
  //     os << "}";
  //     return os;
  //   }
  //   //@}

  // };



  //   template < class DATA_TYPE > 
  //   std::string toString(const Vector< DATA_TYPE >& x) {
  //     std::ostringstream os;
  //     os << x;
  //     return os.str();
  //   }
  
  //   template < class DATA_TYPE > 
  //   std::string toString(const Stack< DATA_TYPE >& x) {
  //     std::ostringstream os;
  //     os << x;
  //     return os.str();
  //   }
  
  //   template < class DATA_TYPE > 
  //   std::string toString(const Array< DATA_TYPE >& x) {
  //     return x.getString();
  //   }

  //   template < class DATA_TYPE > 
  //   std::string toString(const Node< DATA_TYPE >& x) {
  //     return x.getString();
  //   }

  //   template < class DATA_TYPE, int NUM_HEAD > 
  //   std::string toString(const MultiList< DATA_TYPE, NUM_HEAD >& x) {
  //     return x.getString();
  //   }

  //   template< class WORD_TYPE, class FLOAT_TYPE >
  //   std::string toString(const Bitset< WORD_TYPE, FLOAT_TYPE >& x) {
  //     return x.getString();
  //   }



  std::ostream& operator<< (std::ostream& os, const Explanation& x);  

  std::ostream& operator<< (std::ostream& os, const Interval& x);  

  std::ostream& operator<< (std::ostream& os, const BiInterval& x);

  std::ostream& operator<< (std::ostream& os, const MultiSet& x);

  std::ostream& operator<< (std::ostream& os, const IntStack& x);

  std::ostream& operator<< (std::ostream& os, const Queue& x);

  template < int N, class T > 
  std::ostream& operator<< (std::ostream& os, const Tuple< N, T >& x) {
    return x.display(os);
  }

  template < class DATA_TYPE > 
  std::ostream& operator<< (std::ostream& os, const TwoWayStack< DATA_TYPE >& x) {
    return x.display(os);
  }

  template <class T1, class T2, class T3>
  std::ostream& operator<< (std::ostream& os, const Triplet<T1, T2, T3>& x) {
    return x.display(os);
  }


  template < class DATA_TYPE > 
  std::ostream& operator<< (std::ostream& os, const BinaryMinHeap< DATA_TYPE >& x) {
    return  x.display(os);
  }

  template < class DATA_TYPE > 
  std::ostream& operator<< (std::ostream& os, const IndexedBinaryMaxHeap< DATA_TYPE >& x) {
    return  x.display(os);
  }

  template < class DATA_TYPE > 
  std::ostream& operator<< (std::ostream& os, const Vector< DATA_TYPE >& x) {
    return x.display(os);
  }

  // template < class DATA_TYPE > 
  // std::ostream& operator<< (std::ostream& os, const Vector< DATA_TYPE >& x) {
  //   return x.display(os);
  // }

  template < class MAIN_TYPE, class AUX_TYPE > 
  std::ostream& operator<< (std::ostream& os, const BiStack< MAIN_TYPE, AUX_TYPE >& x) {
    return x.display(os);
  }

  // template < class DATA_TYPE > 
  // std::ostream& operator<< (std::ostream& os, const Stack< DATA_TYPE >& x) {
  //   return x.display(os);
  // }

  template < class DATA_TYPE, class SIZE_TYPE > 
  std::ostream& operator<< (std::ostream& os, const VarStack< DATA_TYPE, SIZE_TYPE >& x) {
    return x.display(os);
  }

  template < class DATA_TYPE > 
  std::ostream& operator<< (std::ostream& os, const ConStack< DATA_TYPE >& x) {
    return x.display(os);
  }

  template < class DATA_TYPE > 
  std::ostream& operator<< (std::ostream& os, const Stack< DATA_TYPE >& x) {
    return x.display(os);
  }

  template < class DATA_TYPE > 
  std::ostream& operator<< (std::ostream& os, const Array< DATA_TYPE >& x) {
    return x.display(os);
  }

  template < class DATA_TYPE > 
  std::ostream& operator<< (std::ostream& os, const Node< DATA_TYPE >& x) {
    return x.display(os);
  }

  template < class DATA_TYPE, int NUM_HEAD > 
  std::ostream& operator<< (std::ostream& os, const MultiList< DATA_TYPE, NUM_HEAD >& x) {
    return x.display(os);
  }

  template< class WORD_TYPE, class FLOAT_TYPE >
  std::ostream& operator<< (std::ostream& os, const Bitset< WORD_TYPE, FLOAT_TYPE >& x) {
    return x.display(os);
  }


  std::ostream& operator<< (std::ostream& os, const Explanation* x);  

  std::ostream& operator<< (std::ostream& os, const IntStack* x);

  std::ostream& operator<< (std::ostream& os, const Queue* x);

  std::ostream& operator<< (std::ostream& os, const MultiSet* x);

  template < int N, class T > 
  std::ostream& operator<< (std::ostream& os, const Tuple< N, T >* x) {
    return x->display(os);
  }

  template < class DATA_TYPE > 
  std::ostream& operator<< (std::ostream& os, const TwoWayStack< DATA_TYPE >* x) {
    return x->display(os);
  }

  template <class T1, class T2, class T3>
  std::ostream& operator<< (std::ostream& os, const Triplet<T1, T2, T3>* x) {
    return x->display(os);
  }


  template < class DATA_TYPE > 
  std::ostream& operator<< (std::ostream& os, const BinaryMinHeap< DATA_TYPE >* x) {
    return (x ? x->display(os) : os);
  }

 template < class DATA_TYPE > 
  std::ostream& operator<< (std::ostream& os, const IndexedBinaryMaxHeap< DATA_TYPE >* x) {
    return (x ? x->display(os) : os);
  }

  template < class DATA_TYPE > 
  std::ostream& operator<< (std::ostream& os, const Vector< DATA_TYPE >* x) {
    return (x ? x->display(os) : (os << "nill"));
  }

  // template < class DATA_TYPE > 
  // std::ostream& operator<< (std::ostream& os, const Vector< DATA_TYPE >* x) {
  //   return (x ? x->display(os) : (os << "nill"));
  // }

  template < class MAIN_TYPE, class AUX_TYPE > 
  std::ostream& operator<< (std::ostream& os, const BiStack< MAIN_TYPE, AUX_TYPE >* x) {
    return x->display(os);
  }

  // template < class DATA_TYPE > 
  // std::ostream& operator<< (std::ostream& os, const Stack< DATA_TYPE >* x) {
  //   return (x ? x->display(os) : os);
  // }

  template < class DATA_TYPE > 
  std::ostream& operator<< (std::ostream& os, const Array< DATA_TYPE >* x) {
    return (x ? x->display(os) : os);
  }

  template < class DATA_TYPE > 
  std::ostream& operator<< (std::ostream& os, const Node< DATA_TYPE >* x) {
    return (x ? x->display(os) : os);
  }

  //   template < class DATA_TYPE, int NUM_HEAD > 
  //   std::ostream& operator<< (std::ostream& os, const MultiList< DATA_TYPE, NUM_HEAD >* x) {
  //     return x->display(os);
  //   }

  template< class WORD_TYPE, class FLOAT_TYPE >
  std::ostream& operator<< (std::ostream& os, const Bitset< WORD_TYPE, FLOAT_TYPE >* x) {
    return (x ? x->display(os) : os);
  }

  //   template < class DATA_TYPE, int NUM_HEAD >
  //   std::ostream& MultiList< DATA_TYPE, NUM_HEAD >::display(std::ostream& os) const {
  //     //for(unsigned int i=0; i<data.size; ++i)
  //     //os << data[i];
    
  //     for(int k=0; k<NUM_HEAD; ++k) {
      
  //       int min_elt =  NOVAL;
  //       int max_elt = -NOVAL;
      
  //       Node< DATA_TYPE > nd = data.stack_[head[k]];
  //       while(next(nd)) {
  // 	if(min_elt > (int)nd) min_elt = (int)nd;
  // 	if(max_elt < (int)nd) max_elt = (int)nd;
  //       }
      
  //       if(min_elt == NOVAL) {
  // 	min_elt = 0;
  // 	max_elt = 1;
  //       }
      
  //       BitSet tmp(min_elt, max_elt, BitSet::empt);
      
  //       nd = data.stack_[head[k]];
  //       while(next(nd)) {
  // 	tmp.add((int)nd);
  //       }
      
  //       os << tmp;
  //     }
  //     return os;
  //   }

}


// class IntervalList {

// public:

//   unsigned int size;
//   Vector<int> min;
//   Vector<int> max;

//   bool contain(const int v) {
//     int lmin=0;
//     int lmax=0;

//     bool is_in = true;

//     while(is_in) {
//       if(v < min[lmin] || v > max[lmax]) {
// 	is_in = false;
//       }
//     }
//   }

// };

// // [0,10][23,25][40,100][120,150]
// //    min   max
// /*   -inf   inf
// 	0   150
//        23    10
//       120   100
//        40    25
//       inf  -inf
//  */

#endif // __STRUCTURE_HPP
