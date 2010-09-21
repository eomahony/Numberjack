
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

#include <iostream>
#include <iomanip>
#include <assert.h>

#include <mistral_csp.h>
#include <mistral_con.h>
#include <mistral_sol.h>
#include <mistral_var.h>


using namespace Mistral;
using namespace std;

BitsetIterator* static_bit_domain_it = NULL;
ListIterator* static_list_domain_it = NULL;
BoolIterator* static_bool_domain_it = NULL;
RangeIterator* static_range_domain_it = NULL;
ConstantIterator* static_constant_domain_it = NULL;
VirtualIterator* static_reference_domain_it = NULL;



ReversibleSparseSet::ReversibleSparseSet() //: ReversibleObj() 
{
  size = 0;
}

void ReversibleSparseSet::initialise( Solver* s )
{  
  solver = s;
  value.init(0, size-1, BitSet::empt);
  btlvl_.init(0, size);
  decided.init(0, size-1);
}
// void ReversibleSparseSet::initialise( const int n )
// {
//   solver = s;
//   value.init(0, size-1, BitSet::empt);
//   btlvl_.init(0, size);
//   decided.init(0, size-1);
// }
// void ReversibleSparseSet::pointTo( ReversibleDomainPool& pool )
// {
//   size = pool.size;
//   solver = pool.solver;
//   decided 
// }
ReversibleSparseSet::~ReversibleSparseSet() 
{
}
void ReversibleSparseSet::print(std::ostream& o) const 
{
}
bool ReversibleSparseSet::setValue(const int id, const int var, const int val)
  {
    if( decided.member(var) ) 
      return (value.member(var) == val);
    //saveup();

    solver->triggerEvent(id, Constraint::VALUETRIGGER);

    decided.insert(var);
    if(val)
      value.insert(var);
    else
      value.erase(var);
    return true;
  }

ReversibleBitList::~ReversibleBitList()
{
  delete [] list_;
  delete [] absidx;

  value_.neg_words = 0;
  value_.table = NULL;
}

void ReversibleBitList::setValue(BitSet& e, 
			       const int lb, 
			       const int ub,
			       const unsigned int l)
{
  size = l;
  size_.push(size);
  lvl_.push(0);

  value_.pos_words = e.pos_words;
  value_.neg_words = e.neg_words;
  value_.table = e.table;

  list_ = new int[size+1];
  list_[size] = NOVAL;
  absidx = new unsigned int[ub-lb+1];
  index_ = (absidx - lb);

  unsigned int i;
  if( (int)size == (ub-lb+1) ) 
    {
      for(i=0; i<size; ++i) 
	{
	  index_[lb+i] = i;
	  list_[i] = (lb+i);
	}
    }
  else
    {
      i=0;
      int j=value_.min();
      do {
	index_[j] = i;
	list_[i] = j;
	++i;
      } while( (j = value_.next(j)) != NOVAL );
    }
}

void ReversibleBitList::print(std::ostream& o) const 
{
  value_.print( o );
}


ReversibleIntList::~ReversibleIntList()
{
  delete [] list_;
  delete [] absidx;
}

void ReversibleIntList::setValue(const int lb, 
				 const int ub,
				 const int l)
{
  size = ub-lb+1;
  size_.push(l);
  lvl_.push(0);

  list_ = new int[size+1];
  list_[size] = NOVAL;
  absidx = new unsigned int[size];
  index_ = (absidx - lb);

  for(unsigned int i=0; i<size; ++i) 
    {
      index_[lb+i] = i;
      list_[i] = (lb+i);
    }
  
  size = l;
}

void ReversibleIntList::setValue(BitSet& e,
				 const int lb, 
				 const int ub,
				 const int l)
{
  size = l;
  size_.push(size);
  lvl_.push(0);

  list_ = new int[size+1];
  list_[size] = NOVAL;
  absidx = new unsigned int[ub-lb+1];
  index_ = (absidx - lb);

  unsigned int i;
  if( (int)size == (ub-lb+1) ) 
    {
      for(i=0; i<size; ++i) 
	{
	  index_[lb+i] = i;
	  list_[i] = (lb+i);
	}
    }
  else
    {
      i=0;
      int j=e.min();
      do {
	index_[j] = i;
	list_[i] = j;
	++i;
      } while( (j = e.next(j)) != NOVAL );
    }
}

void ReversibleIntList::print(std::ostream& o) const 
{
  for(unsigned int i=0; i<size; ++i)
    o << " " << list_[i];
  o << " .";
  for(unsigned int i=size; i<size_[0]; ++i)
    o << " " << list_[i];
}

void ReversibleSet::setValue(const int lb, const int ub, const unsigned int val)
{
  clone = false;
  glvl = 0;
  value.init(lb, ub, val);

  int i = value.pos_words;
  int j = value.neg_words;
  lvl_ = new Vector<int>[i-j];
  val_ = new BitSet[i-j];
  lvl_ -= j;
  val_ -= j;

  int lmax=32;
  
//   ////////////////
//  int k;
//   nChanges = new int[BitSet::size_word_bit];
//   words_ = new int*[BitSet::size_word_bit];
//   for(k=0; k<BitSet::size_word_bit; ++k)
//     words_[k] = new int[i-j];
//   ////////////////

  while( i-- > j ) { 
    if( val == BitSet::full )
      lmax = value.size(i);
    val_[i].init( lmax+1, value.table[i] );
    lvl_[i].init( 0, lmax+1 );
    lvl_[i].push(0);
  }
}

void ReversibleSet::setValue(BitSet& e, const int l)
{
  clone = true;
  glvl = 0;
  value.pos_words = e.pos_words;
  value.neg_words = e.neg_words;
  value.table = e.table;
  int i = value.pos_words;
  int j = value.neg_words;
  lvl_ = new Vector<int>[i-j];
  val_ = new BitSet[i-j];
  lvl_ -= e.neg_words;
  val_ -= e.neg_words;

  int lmax;  

//   ////////////////
//   int k;
//   nChanges = new int[BitSet::size_word_bit];
//   words_ = new int*[BitSet::size_word_bit];
//   for(k=0; k<BitSet::size_word_bit; ++k)
//     words_[k] = new int[i-j];
//   ////////////////

  while( i-- > j ) { 
    lmax = value.size(i);
    val_[i].init( lmax+1, value.table[i] );
    lvl_[i].init( 0, lmax+1 );
    lvl_[i].push(0);
  }
}


/**********************************************
 * VariableDomain
 **********************************************/

VariableDomain::VariableDomain( Solver *s ) : VariableInt(s) 
{
  s->binds( vmin );
  s->binds( vmax );
}

// void VariableDomain::triggerEvent( const bool bc )
// {
//   if( vmin == vmax ) {
//     solver->triggerEvent(this, Constraint::VALUETRIGGER);
//   }
//   else if( bc )
//     solver->triggerEvent(this, Constraint::RANGETRIGGER);
//   else solver->triggerEvent(this, Constraint::DOMAINTRIGGER);
// }

void VariableDomain::setVariable(BitSet *s, const int lg, 
				 const int lo, const int up)
{
  isWord = (lo >= 0 && up < 32);
  vmax.setValue( up );
  vmin.setValue( lo );
  posval = up+1;
  negval = lo;
  values.pos_words = s->pos_words;
  values.neg_words = s->neg_words;
  values.table = s->table;
  s->table = NULL;
}

void VariableDomain::setVariable(const int lo, const int up) 
{
  isWord = (lo >= 0 && up < 32);
  vmin.setValue( lo );
  vmax.setValue( up );

  negval = lo;
  posval = up+1;

  values.init(lo, up, BitSet::full);
}

int& VariableDomain::getIntDomain()
{
  return ((int*)values.table)[0];
}

/**********************************************
 * VariableBit
 **********************************************/

void VariableBit::check_domain(int x) const
{
  cout << "check domain " << x << endl;
  assert( values.size() == (unsigned int)size );
  assert( values.max() == vmax );
  assert( values.min() == vmin );

  int real_size = 0;
  for(int i=negval; i<posval; ++i)
    real_size += values.member( i );

  assert( real_size == size );
}

DomainIterator *VariableBit::begin()
{
  //static_bit_domain_it->curval = vmin;
  //static_bit_domain_it->domain.pointTo( values );

  //cout << "here" << endl;

//   int test = vmin;
//   if( test > 1000 ){
//     cout << "here" << endl;
//   }


  static_bit_domain_it->init( vmin, values );
  return static_bit_domain_it;
}

VariableBit::VariableBit(Solver *s, const int v) : VariableDomain(s) 
  { 
    s->binds( size );
    size.setValue(v);
    setVariable(0, v-1); 
    //domain.init(s, values, size);
    s->binds( domain );
    domain.setValue(values, size);
    if( !static_bit_domain_it ) 
      static_bit_domain_it = new BitsetIterator( values );
  }
  /*! 
    Constructor, create a domain containing all values
    in the interval [l,-,u]
  */
VariableBit::VariableBit(Solver *s, const int l, const int u) : VariableDomain(s) 
  { 
    s->binds( size );
    size.setValue( u-l+1 );
    setVariable(l, u); 
    //domain.init(s, values, size);
    s->binds( domain );
    domain.setValue(values, size);
    if( !static_bit_domain_it ) 
      static_bit_domain_it = new BitsetIterator( values );
  }
  /*! 
    Constructor, create a domain containing all values 
    in the set e
  */
 VariableBit::VariableBit(Solver *s, BitSet *e, 
	     const int lg, 
	     const int lo, 
	     const int up) : VariableDomain(s)
  { 
    s->binds( size );
    size.setValue( lg );
    setVariable(e, lg, lo, up); 
    //domain.init(s, values, size);
    s->binds( domain );
    domain.setValue(values, size);
    if( !static_bit_domain_it ) 
      static_bit_domain_it = new BitsetIterator( values );
  }

VariableBit::~VariableBit() 
{
  if( static_bit_domain_it )
    delete static_bit_domain_it;
  static_bit_domain_it = NULL;
}

int VariableBit::random() const
{
  int val;
  if((1 + vmax - vmin ) == size) val = (randint(size) + vmin);
  else {
    val = randint(size);      
    static_bit_domain_it->init( vmin, values );
    
    //DomainIterator *iter = begin();
    while(val--) static_bit_domain_it->next();
    //iter->next();   
    val = (*static_bit_domain_it);
  }
  return val;
  
  //if((1 + vmax - vmin ) != size) return values.random();
  //return (randint(size) + vmin);
}

// bool VariableBit::remove(const int v) 
// {

// #ifdef _CONSISTENCYCHECK
//   check_domain(1);
// #endif

//   if(!contain(v)) return true;
//   if( 1 == size ) return false;

//   domain.erase( v );
//   --size;
  
//   triggerEvent( boundchanged( v ) );

// #ifdef _CONSISTENCYCHECK
//   check_domain(2);
// #endif

//   return true;
// }

// bool VariableBit::setDomain(const int v) 
// {

// #ifdef _CONSISTENCYCHECK
//   check_domain(3);
// #endif

//   if( !contain( v ) ) return false;
   
//   if( 1 < size ) {

//     size = 1;

//     domain.setTo( v );

//     vmin = v;
//     vmax = v;

//     solver->triggerEvent(this, Constraint::VALUETRIGGER);
//   }

// #ifdef _CONSISTENCYCHECK
//   check_domain(4);
// #endif

//   return true;
// }

// bool VariableBit::setMax(const int up) 
// {

// #ifdef _CONSISTENCYCHECK
//   check_domain(5);
// #endif

//   if(vmin > up)  return false;  
//   if(vmax <= up) return true;   

//   domain.setMax(up);
//   size = values.size();

//   if( values.fastMember(up) ) vmax = up;
//   else vmax = values.max();

//   triggerEvent( true );

// #ifdef _CONSISTENCYCHECK
//   check_domain(6);
// #endif

//   return true;
// }

// bool VariableBit::setMin(const int low) 
// {

// #ifdef _CONSISTENCYCHECK
//   check_domain(7);
// #endif

//   if(vmax < low)  return false;  
//   if(vmin >= low) return true;  

//   domain.setMin(low);
//   size = values.size();

//   if( values.fastMember(low) ) vmin = low;
//   else vmin = values.min();

//   triggerEvent( true );

// #ifdef _CONSISTENCYCHECK
//   check_domain(8);
// #endif

//   return true;
// }

// bool VariableBit::setDomain(const int *d, const int n) 
// {
//   BitSet aux(vmin, vmax, BitSet::full);
//   int i=n;
//   while( i-- ) aux.insert( d[i] );
//   return setDomain( aux );
// }

// bool VariableBit::setDomain(VariableInt *x) const 
// {
//   if( 1 == size ) 
//     return x->setDomain( vmin );
//   return x->setDomain( values );
// }

// bool VariableBit::setDomain(const BitSet& s) 
// {

// #ifdef _CONSISTENCYCHECK
//   check_domain(9);
// #endif

//   if( !values.intersect(s) ) return false;
//   if( values.included(s) ) return true;

//   domain.intersectWith(s);
//   size = values.size();

//   triggerEvent( boundchanged( s, true ) );

// #ifdef _CONSISTENCYCHECK
//   check_domain(10);
// #endif

//   return true;
// }

// bool VariableBit::removeRange(const int lo, const int up) 
// {

// #ifdef _CONSISTENCYCHECK
//   check_domain(11);
// #endif

//   if( !(values.intersect( lo, up )) )
//     return true;

//   if( lo <= vmin ) 
//     return setMin(up+1);
//   if( up >= vmax ) 
//     return setMax(lo-1);  

//   domain.removeInterval(lo, up);
//   size = values.size();

//   solver->triggerEvent(this, Constraint::DOMAINTRIGGER);  

// #ifdef _CONSISTENCYCHECK
//   check_domain(12);
// #endif

//   return true;
// }

// bool VariableBit::removeSet(const BitSet& s) 
// {  

// #ifdef _CONSISTENCYCHECK
//   check_domain(13);
// #endif

//   if( values.included(s) ) return false;
//   if( !values.intersect(s) ) return true;

//   domain.setminusWith(s);    
//   size = values.size();

//   triggerEvent( boundchanged( s, false ) );

// #ifdef _CONSISTENCYCHECK
//   check_domain(14);
// #endif

//   return true;
// }

// inline bool VariableBit::revise( BitSet* supports, 
// 				 int* residues, 
// 				 VariableInt* X)
// {
//   int vali = vmin;
//   bool consistent;
//   if( X->isWord )
//     do consistent = (X->wordIntersect( supports[vali] ) || remove( vali ));
//     while( consistent && setNext(vali) );
//   else 
//     do consistent = (X->intersect( supports[vali], residues[vali] ) || remove( vali ));
//     while( consistent && setNext(vali) );
//   return consistent;
// }

// inline bool VariableBit::revise( Constraint* c, const int i )
// {
//   int vali = vmin;
//   do {
//     // This procedure does two things:
//     // * first it checks if we can say that the value is supported
//     //   without constraint check
//     // * second it ready the solution in order to be check
//     //   by findSupport	
    
//     if( !( c->firstSupport(i, vali) || 
// 	   c->findSupport(i, vali) ||
// 	   remove(vali) ) )
//       return false;	
	
//   } while( setNext(vali) );
  
//   return true;
// }

void VariableBit::print(std::ostream& o) const 
{
  printshort(o);
  if( isGround() ) {
    o << " = " << min();
  } else {
    o << " in ";
    printDomain(o);
    //o << " w=" << weight;
  }
}

void VariableBit::printshort(std::ostream& o) const 
{
  o << "x" << id+1 ;
}

void VariableBit::printDomain(std::ostream& o) const 
{
  values.print( o );
}


/**********************************************
 * VariableBitset
 **********************************************/

void VariableBitset::check_domain(int x) const
{
  cout << "check domain " << x << endl;
  assert( values.size() == (unsigned int)size );
  assert( values.max() == vmax );
  assert( values.min() == vmin );

  int real_size = 0;
  for(int i=negval; i<posval; ++i)
    real_size += values.member( i );

  assert( real_size == size );
}

DomainIterator *VariableBitset::begin()
{
  //static_bit_domain_it->curval = vmin;
  //static_bit_domain_it->domain.pointTo( values );

  static_bit_domain_it->init( vmin, values );
  return static_bit_domain_it;
}

VariableBitset::VariableBitset(Solver *s, const int v) 
  : VariableDomain(s) 
  { 
    s->binds( size );
    size.setValue(v);
    setVariable(0, v-1); 
    domain.init(values);
    for(int i=domain.neg_words; i<domain.pos_words; ++i) {
      s->binds( domain.table[i] );
    }
    if( !static_bit_domain_it ) 
      static_bit_domain_it = new BitsetIterator( values );
  }
  /*! 
    Constructor, create a domain containing all values
    in the interval [l,-,u]
  */
VariableBitset::VariableBitset(Solver *s, const int l, const int u) 
  : VariableDomain(s) 
  { 
    s->binds( size );
    size.setValue( u-l+1 );
    setVariable(l, u); 
    domain.init(values);
    for(int i=domain.neg_words; i<domain.pos_words; ++i) {
      s->binds( domain.table[i] );
    }
    if( !static_bit_domain_it ) 
      static_bit_domain_it = new BitsetIterator( values );
  }
  /*! 
    Constructor, create a domain containing all values 
    in the set e
  */
VariableBitset::VariableBitset(Solver *s, BitSet *e, 
			       const int lg, 
			       const int lo, 
			       const int up) 
  : VariableDomain(s)
{ 
  s->binds( size );
  size.setValue( lg );
  setVariable(e, lg, lo, up); 
  domain.init(values);
  for(int i=domain.neg_words; i<domain.pos_words; ++i) {
    s->binds( domain.table[i] );
  }
  if( !static_bit_domain_it ) 
    static_bit_domain_it = new BitsetIterator( values );
}

VariableBitset::~VariableBitset() 
{
  if( static_bit_domain_it )
    delete static_bit_domain_it;
  static_bit_domain_it = NULL;
}

void VariableBitset::print(std::ostream& o) const 
{
  //o << "x" << id ;
  printshort( o );
  if( isGround() )
    o << " = " << min();
  else {
    o << " in ";
    values.print( o );
    cout << " / " ;
    domain.print( o );
    cout << " / [" << (vmin) << ".." << (vmax) << "] " << (size);
  }
  checkDomain( "print" );
}

void VariableBitset::printshort(std::ostream& o) const 
{
  o << "x" << id+1 ;
}

void VariableBitset::printDomain(std::ostream& o) const 
{
  values.print( o );
}

int VariableBitset::random() const
{
  int val;
  if((1 + vmax - vmin ) == size) val = (randint(size) + vmin);
  else {
    val = randint(size);      
    static_bit_domain_it->init( vmin, values );
    
    //DomainIterator *iter = begin();
    while(val--) static_bit_domain_it->next();
    //iter->next();   
    val = (*static_bit_domain_it);
  }
  return val;
}



/**********************************************
 * VariableBitList
 **********************************************/

void VariableBitList::check_domain(int x) const
{
  cout << "check domain " << x << endl;
  assert( values.size() == domain.size );
  assert( values.max() == vmax );
  assert( values.min() == vmin );

  unsigned int real_size = 0;
  for(int i=negval; i<posval; ++i)
    real_size += values.member( i );

  assert( real_size == domain.size );
}

VariableBitList::VariableBitList(Solver *s, const int v) : VariableDomain(s) 
{ 
  setVariable(0, v-1); 
  //domain.init(s, values, 0, v-1, v);
  s->binds( domain );
  domain.setValue(values, 0, v-1, v);
  if( !static_list_domain_it ) static_list_domain_it = new ListIterator();
}

VariableBitList::VariableBitList(Solver *s, const int l, const int u) : VariableDomain(s) 
{ 
  setVariable(l, u); 
  //domain.init(s, values, l, u, u-l+1);
  s->binds( domain );
  domain.setValue(values, l, u, u-l+1);
  if( !static_list_domain_it ) static_list_domain_it = new ListIterator();
}

VariableBitList::VariableBitList(Solver *s, 
			   BitSet *e, 
			   const int lg, 
			   const int lo, 
			   const int up) : VariableDomain(s)
{ 
  setVariable(e, lg, lo, up); 
  //domain.init(s, values, lo, up, lg);
  s->binds( domain );
  domain.setValue(values, lo, up, lg);
  if( !static_list_domain_it ) static_list_domain_it = new ListIterator();
}

DomainIterator *VariableBitList::begin()
{
  static_list_domain_it->size = domain.size-1;
  static_list_domain_it->vals = domain.list_;
  return static_list_domain_it;
}

VariableBitList::~VariableBitList() 
{
  if( static_list_domain_it )
    delete static_list_domain_it;
  static_list_domain_it = NULL;
}

void VariableBitList::print(std::ostream& o) const 
{
  //o << "y" << id ;
  printshort( o );
  if( isGround() )
    o << " = " << min();
  else {
    o << " in ";
    values.print( o );
  }
}

void VariableBitList::printshort(std::ostream& o) const 
{
  o << "y" << id+1 ;
}

void VariableBitList::printDomain(std::ostream& o) const 
{
  values.print( o );
}



// bool VariableBitList::remove(const int v) 
// {

// #ifdef _CONSISTENCYCHECK
//   check_domain(1);
// #endif

//   if(!contain(v)) return true;
//   if( 1 == domain.size ) return false;

//   domain.erase( v );
  
//   triggerEvent( boundchanged( v ) );

// #ifdef _CONSISTENCYCHECK
//   check_domain(2);
// #endif

//   return true;
// }

// bool VariableBitList::setDomain(const int v) 
// {

// #ifdef _CONSISTENCYCHECK
//   check_domain(3);
// #endif

//   if( !contain( v ) ) return false;
   
//   if( 1 < domain.size ) {

//     domain.setTo( v );

//     vmin = v;
//     vmax = v;

//     solver->triggerEvent(this, Constraint::VALUETRIGGER);
//   }

// #ifdef _CONSISTENCYCHECK
//   check_domain(4);
// #endif

//   return true;
// }

// bool VariableBitList::setMax(const int up) 
// {

// // #ifdef _CONSISTENCYCHECK
// //   check_domain(5);
// // #endif

// //   if(vmin > up)  return false;  
// //   if(vmax <= up) return true;   

// //   domain.setMax(up);

// //   if( values.member(up) ) vmax = up;
// //   else vmax = values.max();

// //   triggerEvent( true );

// // #ifdef _CONSISTENCYCHECK
// //   check_domain(6);
// // #endif

// //TODO

//   return true;
// }

// bool VariableBitList::setMin(const int low) 
// {

// // #ifdef _CONSISTENCYCHECK
// //   check_domain(7);
// // #endif

// //   if(vmax < low)  return false;  
// //   if(vmin >= low) return true;  

// //   domain.setMin(low);
// //   size = values.size();

// //   if( values.member(low) ) vmin = low;
// //   else vmin = values.min();

// //   triggerEvent( true );

// // #ifdef _CONSISTENCYCHECK
// //   check_domain(8);
// // #endif

// //TODO

//   return true;
// }

// bool VariableBitList::setDomain(const int *d, const int n) 
// {
// //   BitSet aux(vmin, vmax);
// //   int i=n;
// //   while( i-- ) aux.insert( d[i] );
// //   return setDomain( aux );

// //TODO

// }

// bool VariableBitList::setDomain(VariableInt *x) const 
// {
//   if( 1 == domain.size ) 
//     return x->setDomain( vmin );
//   return x->setDomain( values );

// //TODO

// }

// bool VariableBitList::setDomain(const BitSet& s) 
// {

// // #ifdef _CONSISTENCYCHECK
// //   check_domain(9);
// // #endif

// //   if( !values.intersect(s) ) return false;
// //   if( values.included(s) ) return true;

// //   domain.intersectWith(s);
// //   size = values.size();

// //   triggerEvent( boundchanged( s, true ) );

// // #ifdef _CONSISTENCYCHECK
// //   check_domain(10);
// // #endif

// //TODO

//   return true;
// }

// bool VariableBitList::removeRange(const int lo, const int up) 
// {

// // #ifdef _CONSISTENCYCHECK
// //   check_domain(11);
// // #endif

// //   if( !(values.intersect( lo, up )) ) return true;
// //   if( lo <= vmin ) return setMin(up+1);
// //   if( up >= vmax ) return setMax(lo-1);  

// //   domain.removeInterval(lo, up);
// //   size = values.size();

// //   solver->triggerEvent(this, Constraint::DOMAINTRIGGER);  

// // #ifdef _CONSISTENCYCHECK
// //   check_domain(12);
// // #endif

// //TODO

//   return true;
// }

// bool VariableBitList::removeSet(const BitSet& s) 
// {  

// // #ifdef _CONSISTENCYCHECK
// //   check_domain(13);
// // #endif

// //   if( values.included(s) ) return false;
// //   if( !values.intersect(s) ) return true;

// //   domain.setminusWith(s);    
// //   size = values.size();

// //   triggerEvent( boundchanged( s, false ) );

// // #ifdef _CONSISTENCYCHECK
// //   check_domain(14);
// // #endif

// //TODO

//   return true;
// }

// inline bool VariableBitList::revise( BitSet* supports, 
// 				  int* residues, 
// 				  VariableInt* X)
// {
//   bool consistent = true;
//   int *val = domain.list_, i=domain.size;
//   if( X->isWord )
//     while( consistent && i-- ) 
//       consistent = (X->wordIntersect( supports[val[i]] ) || remove( val[i] ));
//   else 
//     while( consistent && i-- )
//       consistent = (X->intersect( supports[val[i]], residues[val[i]] ) || remove( val[i] ));
//   return consistent;

// //   int *val = domain.list_, i=domain.size;
// //   if( X->isWord )
// //     while( i-- && (X->wordIntersect( supports[val[i]] ) || remove( val[i] )) ) ;
// //   else 
// //     while( i-- && (X->intersect( supports[val[i]], residues[val[i]] ) || remove( val[i] )) ) ;
// //   return (i < 0);


// //   //bool consistent = true;
// //   int *val = domain.list_, i=domain.size;
// //   if( isWord )
// //     while( i-- ) {
// //       if( !(X->wordIntersect( supports[val[i]] ) || remove( val[i] )) ) break;
// //     } 
// //   else 
// //     while( i-- ) {
// //       if( !(X->intersect( supports[val[i]], 
// // 			  residues[val[i]] ) || remove( val[i] )) ) break;  
// //     }
// //   return ( i < 0 );


// //   int *val = domain.list_, i=domain.size;
// //   if( isWord )
// //     while( i-- )
// //       if(!X->wordIntersect( supports[val[i]] ) && !remove( val[i] )) break;
// //   else 
// //     while( i-- )
// //       if(!X->intersect( supports[val[i]], 
// // 			residues[val[i]] ) && !remove( val[i] )) break;
// //   return (i < 0);

// }

// inline bool VariableBitList::revise( Constraint* c, const int i )
// {
//   bool consistent = true;
//   int *val = domain.list_, j=domain.size;
//   while( consistent && j-- ) {
//     consistent = ( c->firstSupport(i, val[j]) || 
// 		   c->findSupport (i, val[j]) ||
// 		   remove(val[j]) );
//   }
//   return consistent;
// }



/**********************************************
 * VariableList
 **********************************************/

void VariableList::check_domain(int x) const
{
  cout << "check domain " << x << endl;
}

VariableList::VariableList(Solver *s, const int v) : VariableInt(s) 
{ 
  s->binds( domain );
  domain.setValue(0, v-1, v);
  if( !static_list_domain_it ) static_list_domain_it = new ListIterator();
}

VariableList::VariableList(Solver *s, const int l, const int u) : VariableInt(s) 
{ 

  //cout << l << " - " << u << endl;

  s->binds( domain );
  domain.setValue(l, u, u-l+1);
  

  if( !static_list_domain_it ) static_list_domain_it = new ListIterator();
}

VariableList::VariableList(Solver *s, 
			   BitSet* e, 
			   const int lg, 
			   const int lo, 
			   const int up) : VariableInt(s)
{ 
  s->binds( domain );
  domain.setValue(*e, lo, up, lg);
  if( !static_list_domain_it ) static_list_domain_it = new ListIterator();
}

DomainIterator *VariableList::begin()
{
  static_list_domain_it->size = domain.size-1;
  static_list_domain_it->vals = domain.list_;
  return static_list_domain_it;
}

VariableList::~VariableList() 
{
  if( static_list_domain_it )
    delete static_list_domain_it;
  static_list_domain_it = NULL;
}

void VariableList::print(std::ostream& o) const 
{
  //o << "z" << id ;
  printshort( o );
  if( isGround() )
    o << " = " << first();
  else {
    o << " in {";
    for(unsigned int i=0; i<domain.size-1; ++i)
      cout << domain.list_[i] << ", ";
    cout << domain.list_[domain.size-1] << "}";
    //values.print( o );
  }
}

void VariableList::printshort(std::ostream& o) const 
{
  o << "z" << id+1 ;
}

void VariableList::printDomain(std::ostream& o) const 
{
  //values.print( o );
}



/**********************************************
 * VariableBool
 **********************************************/

VariableBool::VariableBool(Solver *s) : VariableInt(s) 
{
  //domain.init(s);
  if(s) s->binds( domain );
  if( !static_bool_domain_it ) static_bool_domain_it = new BoolIterator();
}

int& VariableBool::getIntDomain()
{
  return domain.state;
}

DomainIterator *VariableBool::begin()
{
  static_bool_domain_it->curval = domain.state;//!(domain.state & 1);
  return static_bool_domain_it;
}

VariableBool::~VariableBool() 
{
  if( static_bool_domain_it )
    delete static_bool_domain_it;
  static_bool_domain_it = NULL;
}

inline bool VariableBool::setState( const int nstat ) 
{
  if( !nstat ) return false;
  if( nstat == domain.state ) return true;
  domain = nstat;
  //solver->triggerEvent(this, Constraint::VALUETRIGGER);
  solver->triggerEvent(id, Constraint::VALUETRIGGER);
  return true;
}

bool VariableBool::remove(const int v) 
{
//   if(!((1 << v) & 3)) return true;
  if(v>1 || v<0) return true;
  return setState( domain.state & (2-v) );
}

bool VariableBool::setDomain(const int v) 
{
  if((1 << v) & 3) 
    return setState( domain.state & (v+1) );
  else return false;
}

bool VariableBool::setMax(const int up) 
{
  if( up > 0 ) return true;
  if( up < 0 ) return false;
  return setState( domain.state & 1 );
}

bool VariableBool::setMin(const int low) 
{
  if( low < 1 ) return true;
  if( low > 1 ) return false;
  return setState( domain.state & 2 );
}

bool VariableBool::setDomain(const int *d, const int n) 
{
  int nstat = 0, i=n;
  while( i-- ) 
    if( (nstat |= ((1 << d[i]) & 3) ) == 3 )  break;
  return setState( nstat );
}

bool VariableBool::setDomain(VariableInt* x) const
{ 
  if( domain.state == 3 )
    return ( x->setMin(0) && x->setMax(1) );
  else return ( x->setDomain(domain.state-1) );
}

bool VariableBool::setDomain(const BitSet& s) 
{
  if( s.pos_words < 1 || s.neg_words > 0 ) return false;
  return setState( domain.state & s.table[0] );
}

bool VariableBool::removeRange(const int lo, const int up) 
{
  if( lo > up || lo > 1 || up < 0 ) return true;
  int nstat = 3;
  if(lo < 1) nstat ^= 1;
  if(up > 0) nstat ^= 2;
  return setState( nstat );
}

bool VariableBool::removeSet(const BitSet& s) 
{
  if( s.pos_words < 1 || s.neg_words > 0 ) return true;
  return setState( (domain.state ^ s.table[0]) & 3 );
}

inline bool VariableBool::revise( BitSet* supports, 
				  int* residues, 
				  VariableInt* X)
{
  int i=2;
  while( i-- )
    if( !(X->intersect( supports[i], 
			residues[i] ) ||
	  remove( i )) ) return false;
  return true;
}

inline bool VariableBool::revise( Constraint* c, const int i )
{
  int vali=2;
  while( vali-- )    
    if( !( c->firstSupport(i, vali) || 
	   c->findSupport(i, vali) ||
	   remove(vali) ) )
      return false;
  return true;
}

void VariableBool::print(std::ostream& o) const 
{
  printshort( o );  
  if( domain.state == 3 )
    o << " in {0,1}";
  else o << " = " << (domain.state-1) ;
}

void VariableBool::printshort(std::ostream& o) const 
{
  o << "b" << id+1 ;
}

void VariableBool::printDomain(std::ostream& o) const 
{
  if( domain.state == 3 )
    o << "{0,1}";
  else o << "{" << (domain.state-1) << "}";
}



// /**********************************************
//  * VariableBool2
//  **********************************************/

// VariableBool2::VariableBool2(Solver *s) 
//   : VariableInt(s), domainPool(s->getDomainPool())
// {
//   //domain.init(s);
//   //s->binds( domain );
//   //domainPool = s->getDomainPool( index );

//   index = s->suscribeToPool();
// //   cout << index ;

// //   cout << " " ;

// //   print( cout );

// //   cout << endl;

//   if( !static_bool_domain_it ) static_bool_domain_it = new BoolIterator();
// }

// DomainIterator *VariableBool2::begin()
// {
//   static_bool_domain_it->curval = domainPool.getValue(index)+1;//!(domain.state & 1);
//   return static_bool_domain_it;
// }

// VariableBool2::~VariableBool2() 
// {
//   if( static_bool_domain_it )
//     delete static_bool_domain_it;
//   static_bool_domain_it = NULL;
// }

// // inline bool VariableBool2::setState( const int nstat ) 
// // {
// //   if( !nstat ) return false;
  
// //   if( nstat == domain.state ) return true;
// //   domain = nstat;
  
// //   //if( nstat == domainPool.getValue(index) ) return true;
// //   //domainPool.setValue(index, nstat);

// //   //solver->triggerEvent(this, Constraint::VALUETRIGGER);
// //   solver->triggerEvent(id, Constraint::VALUETRIGGER);
// //   return true;
// // }


// // inline bool VariableBool2::setVal( const int nstat ) 
// // {
// //   if( nstat < 0 ) return false;
  

// // //   if( nstat == domain.state ) return true;
// // //   domain = nstat;
  
// //   //if( nstat == domainPool.getValue(index) ) return true;
// //   domainPool.setValue(id, index, nstat);

// //   //solver->triggerEvent(this, Constraint::VALUETRIGGER);
// //   //solver->triggerEvent(id, Constraint::VALUETRIGGER);
// //   return true;
// // }

// bool VariableBool2::remove(const int v) 
// {
// //   if(!((1 << v) & 3)) return true;
//   if(v>1 || v<0) return true;
  
//   //int oval = (domain.state & (2-v));
//   //int nval = !v

//   //bool old = setState( oval );
  
//   bool cur = domainPool.setValue( id, index, !v );

//   //bool cur = setState( !v );
//   //if( cur )
//     //solver->triggerEvent(id, Constraint::VALUETRIGGER);

//   //assert( old == cur );
  
//   //assert( domainPool.getValue(index)+1 == domain.state );


//   return cur;

// }

// bool VariableBool2::setDomain(const int v) 
// {
//   if((1 << v) & 3) {

//     //bool old = setState( domain.state & (v+1) );

//     bool cur = domainPool.setValue( id, index, v );
//     //if( cur )
//     //solver->triggerEvent(id, Constraint::VALUETRIGGER);
    
//     //bool cur = setState( v );

//     //assert( old == cur );
    
//     //assert( domainPool.getValue(index)+1 == domain.state );

//     return cur;

//   } else return false;
// }

// bool VariableBool2::setMax(const int up) 
// {
//   if( up > 0 ) return true;
//   if( up < 0 ) return false;

//   //bool old = setState( domain.state & 1 );

//   bool cur = domainPool.setValue( id, index, 0 );
//   //if( cur )
//   //solver->triggerEvent(id, Constraint::VALUETRIGGER);
  
//   //bool cur = setState( 0 );

//   //assert( cur == old );
  
//   //assert( domainPool.getValue(index)+1 == domain.state );
  
//   return cur;

// }

// bool VariableBool2::setMin(const int low) 
// {
//   if( low < 1 ) return true;
//   if( low > 1 ) return false;

//   //bool old = setState( domain.state & 2 );

//   bool cur = domainPool.setValue( id, index, 1 );
//   //if( cur )
//   //solver->triggerEvent(id, Constraint::VALUETRIGGER);

//   //bool cur = setState( 1 );

//   //assert( cur == old );

//   //assert( domainPool.getValue(index)+1 == domain.state );
  
//   return cur;

// }

// bool VariableBool2::setDomain(const int *d, const int n) 
// {
//   int nstat = 0, i=n;
//   while( i-- ) 
//     if( (nstat |= ((1 << d[i]) & 3) ) == 3 )  break;

//   //bool old = setState( nstat );

//   --nstat;

//   bool cur = ( nstat>1 || domainPool.setValue( id, index, nstat ) );
//   //if( cur )
//   //solver->triggerEvent(id, Constraint::VALUETRIGGER);
    
//   //assert( cur == old );

//   //assert( domainPool.getValue(index)+1 == domain.state );
  
//   return cur;

// }

// bool VariableBool2::setDomain(VariableInt* x) const
// { 
//   int val = domainPool.getValue(index);
//   if( val == 2 )
//     return ( x->setMin(0) && x->setMax(1) );
//   else return ( x->setDomain(val) );
// }

// bool VariableBool2::setDomain(const BitSet& s) 
// {
//   if( s.pos_words < 1 || s.neg_words > 0 ) return false;

//   //bool old = setState( domain.state & s.table[0] );
    
//   int vals = ((domainPool.getValue(index)+1) & s.table[0]);
//   bool cur =( (--vals == 2) || domainPool.setValue(id, index, vals) );
//   //if( cur )
//   //solver->triggerEvent(id, Constraint::VALUETRIGGER);
  
//   //assert( cur == old );

//   //assert( domainPool.getValue(index)+1 == domain.state );
  
//   return cur;

// }

// bool VariableBool2::removeRange(const int lo, const int up) 
// {
//   if( lo > up || lo > 1 || up < 0 ) return true;
//   int nstat = 3;
//   if(lo < 1) nstat ^= 1;
//   if(up > 0) nstat ^= 2;

//   //bool old = setState( nstat );
   
//   bool cur =( (--nstat == 2) || domainPool.setValue(id, index, nstat) );
//   //if( cur )
//   //solver->triggerEvent(id, Constraint::VALUETRIGGER);
  
//   //assert( cur == old );
    
//   //assert( domainPool.getValue(index)+1 == domain.state );
  
//   return cur;


// }

// bool VariableBool2::removeSet(const BitSet& s) 
// {
//   if( s.pos_words < 1 || s.neg_words > 0 ) return true;
  
//   //bool old = setState( (domain.state ^ s.table[0]) & 3 );

//   int nstat = (((domainPool.getValue(index)+1) ^ s.table[0]) & 3);
//   bool cur =( (--nstat == 2) || domainPool.setValue(id, index, nstat) );
//   //if( cur )
//   //solver->triggerEvent(id, Constraint::VALUETRIGGER);
  
//   //assert( cur == old );
  
//   //assert( domainPool.getValue(index)+1 == domain.state );

//   return cur;

// }

// inline bool VariableBool2::revise( BitSet* supports, 
// 				   int* residues, 
// 				   VariableInt* X)
// {
//   int i=2;
//   while( i-- )
//     if( !(X->intersect( supports[i], 
// 			residues[i] ) ||
// 	  remove( i )) ) return false;
//   return true;
// }

// inline bool VariableBool2::revise( Constraint* c, const int i )
// {
//   int vali=2;
//   while( vali-- )    
//     if( !( c->firstSupport(i, vali) || 
// 	   c->findSupport(i, vali) ||
// 	   remove(vali) ) )
//       return false;
//   return true;
// }

// void VariableBool2::print(std::ostream& o) const 
// {
//   o << "b" << id ;
//   if( domainPool.getValue(index) == 2 )
//     o << " in {0,1}";
//   else o << " = " << (domainPool.getValue(index)) ;
// }

// void VariableBool2::printshort(std::ostream& o) const 
// {
//   o << "b" << id ;
// }

// void VariableBool2::printDomain(std::ostream& o) const 
// {
//   if( domainPool.getValue(index) == 2 )
//     o << "{0,1}";
//   else o << "{" << domainPool.getValue(index) << "}";
// }


/**********************************************
 * VariableRange
 **********************************************/

VariableRange::VariableRange(Solver *s, const int v) : VariableInt(s) 
{ 
  setVariable(0, v-1); 
  s->binds( vmin );
  s->binds( vmax );
  if( !static_range_domain_it ) 
    static_range_domain_it = new RangeIterator();
}

VariableRange::VariableRange(Solver *s, const int l, const int u) : VariableInt(s) 
{ 
  setVariable(l, u);
  s->binds( vmin );
  s->binds( vmax ); 
  if( !static_range_domain_it ) 
    static_range_domain_it = new RangeIterator();
}

VariableRange::VariableRange(Solver *s) : VariableInt(s) 
{  
  if( !static_range_domain_it ) 
    static_range_domain_it = new RangeIterator();
}

DomainIterator *VariableRange::begin()
{
  static_range_domain_it->curval = vmin;
  static_range_domain_it->ubound = vmax;
  return static_range_domain_it;
}

VariableRange::~VariableRange() 
{
  if( static_range_domain_it ) {
    delete static_range_domain_it;
  }
  static_range_domain_it = NULL;
}

int& VariableRange::getIntDomain()
{
  // BUGGY!
  return vmin.value;
}

void VariableRange::setVariable(const int low, const int up) 
{
  vmax.setValue(up);
  vmin.setValue(low);
}


bool VariableRange::remove(const int v) 
{
  if( v > vmax || v < vmin ) return true;
  if( (int)vmin == (int)vmax ) return false;

  if( v == vmin ) 
    ++vmin;
  else if( v == vmax )
    --vmax;
  else return true;

  if((int)vmin < (int)vmax)
    //solver->triggerEvent(this, Constraint::RANGETRIGGER);
    solver->triggerEvent(id, Constraint::RANGETRIGGER);
  else //solver->triggerEvent(this, Constraint::VALUETRIGGER);
    solver->triggerEvent(id, Constraint::VALUETRIGGER);

  return true;
}

bool VariableRange::setDomain(VariableInt* x) const
{ 
  if( (int)vmin == (int)vmax )
    return x->setDomain( vmin );
  return ( x->setMin(vmin) && x->setMax(vmax) );
}

bool VariableRange::setDomain(const int v) 
{
  if( v < vmin || v > vmax )
    return false;
  if( (int)vmin < (int)vmax ) {
    //solver->triggerEvent(this, Constraint::VALUETRIGGER);
    solver->triggerEvent(id, Constraint::VALUETRIGGER);
    vmin = v;
    vmax = v;
  }
  return true;
}

bool VariableRange::setMax(const int up) 
{
  if(vmin > up)
    return false;
  if(vmax <= up)
    return true;  

  vmax = up;
  if( (int)vmin == (int)vmax )
    //solver->triggerEvent(this, Constraint::VALUETRIGGER);
    solver->triggerEvent(id, Constraint::VALUETRIGGER);
  else //solver->triggerEvent(this, Constraint::RANGETRIGGER);
    solver->triggerEvent(id, Constraint::RANGETRIGGER);

  return true;
}

bool VariableRange::setMin(const int low) 
{
  if(vmax < low)
    return false;
  if(vmin >= low)
    return true;  

  vmin = low;
  if( (int)vmin == (int)vmax )
    //solver->triggerEvent(this, Constraint::VALUETRIGGER);
    solver->triggerEvent(id, Constraint::VALUETRIGGER);
  else //solver->triggerEvent(this, Constraint::RANGETRIGGER);
    solver->triggerEvent(id, Constraint::RANGETRIGGER);

  return true;
}

bool VariableRange::setDomain (const int *d, const int n) 
{
  int lb = NOVAL;
  int ub = -1*NOVAL;
  int i=n;
  while( i-- ) {
    if(d[i] > ub) ub=d[i];
    if(d[i] < lb) lb=d[i];
  }    
  return(setMin(lb) && setMax(ub));
}

bool VariableRange::setDomain (const BitSet& s) 
{
  int lb = s.next( vmin-1 );
  int ub = s.prev( vmax+1 );

  return(setMin(lb) && setMax(ub));
}

bool VariableRange::removeRange(const int lo, const int up) 
{
  if( lo > vmax || up < vmin ) return true;
  if( lo <= vmin )
    return setMin(up+1);
  if( up >= vmax ) 
    return setMax(lo-1);  
  return true;
}

bool VariableRange::removeSet(const BitSet& s) 
{
  int lb = vmin;
  int ub = vmax;
  while( s.member(lb) ) ++lb; 
  while( s.member(ub) ) --ub; 
  return(setMin(lb) && setMax(ub));
}

inline bool VariableRange::revise( BitSet* supports, 
				   int* residues, 
				   VariableInt* X)
{
  // compute ub
  int i=vmax;
  while( !(X->intersect( supports[i], 
			 residues[i] )) ) --i;
  if( !setMax(i) ) return false;

  // compute lb
  i=vmin;
  while( !(X->intersect( supports[i], 
			 residues[i] )) ) ++i;
  if( !setMin(i) ) return false;

  return true;
}


inline bool VariableRange::revise( Constraint *c, const int i )
{
  // compute ub
  int vali=vmax;
  while ( !( c->firstSupport(i, vali) || 
	     c->findSupport(i, vali) ) ) --vali;
  if( !setMax(vali) ) return false;

  // compute lb
  vali=vmin;
  while ( !( c->firstSupport(i, vali) || 
	     c->findSupport(i, vali) ) ) ++vali;
  if( !setMin(vali) ) return false;

  return true;
}

void VariableRange::print(std::ostream& o) const
{
  printshort( o );
  if( isGround() )
    o << " = " << vmin;
  else 
    o << "[" << vmin << ".." << vmax << "]";
}

void VariableRange::printshort(std::ostream& o) const
{
  o << "r" << id+1;
}

void VariableRange::printDomain(std::ostream& o) const 
{
  o << "{" << vmin << ".." << vmax << "}";
}



Constant::Constant(Solver *s, const int v ) : VariableInt(s) { 
  val = v; 
  if( !static_constant_domain_it ) static_constant_domain_it = new ConstantIterator( );
}

Constant::~Constant() {
  if( static_constant_domain_it )
    delete static_constant_domain_it;
  static_constant_domain_it = NULL;
}

int& Constant::getIntDomain()
{
  return val;
}


DomainIterator *Constant::begin()
{ 
  static_constant_domain_it->curval = val;
  return static_constant_domain_it;
}


/**********************************************
 * VariableVirtual
 **********************************************/

VariableVirtual::VariableVirtual(Solver *s) : VariableInt(s) 
{ 
  if( !static_reference_domain_it ) 
    static_reference_domain_it = new VirtualIterator();
}

VariableVirtual::~VariableVirtual() 
{ 
  if( static_reference_domain_it )     
    delete static_reference_domain_it;
  static_reference_domain_it = NULL;
  delete conversion;
}

int& VariableVirtual::getIntDomain()
{
  // BUGGY!
  return reference->getIntDomain();
}

DomainIterator *VariableVirtual::begin()
{ 
  static_reference_domain_it->init(conversion, reference);
  return static_reference_domain_it;
}


void MapOffset::print(ostream& o, const VariableInt* x) const
  {
    //x->printshort(o);
    x->print(o);
    if( K > 0 )
      o << " + " << K;
    else if( K < 0 )
      o << " - " << (-K);
  };


void MapFactor::print(ostream& o, const VariableInt* x) const
  {
    
    //x->printshort(o);
    x->print(o);
    o << " * " << K;
  };


void MapQuotient::print(ostream& o, const VariableInt* x) const
  {
    x->printshort(o);
    o << " / " << K;
  };

