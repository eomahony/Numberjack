
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

#include <mistral_mod.h>
#include <mistral_sol.h>
#include <cstring>
//#include "../../src/Mistral.hpp"

using namespace Mistral;
using namespace std;


static ConstraintStore ENVIRONMENT;

int P_GAC;

// Mistral_Expression* BuildObject::ExpressionWrap() 
// {
//   Mistral_Expression *x = new Mistral_Expression(getMin(), getMax());
//   return x;
// }

void BuildObject::add( CSP *m, const int l )
{ 
  //if(l) 
  if(l || !isRel()) // reference "artificially" the variables added at the top level 
    reference();
  if( id == NOVAL ) {
    model = m;
    id = model->declarations.size ;
    model->declarations.push( this ) ;
  }
}

// void BuildObject::gather_variables()
// {
//   model->variables.insert(id);
// }

void BuildObject::close()
{  
  //cout << (isMaperenced()) << " " << (!veqptr_) << " " << (min_ < 0 || max_ > 1) << endl;
  if( isReferenced() && !veqptr_ && (min_ < 0 || max_ > 1) ) model->unsetSat();
}

/**********************************************
 * Variable BuildObject
 **********************************************/

void BuildObject::initVar(const int lb, const int ub) 
{
  refs = 0;
  //id = -1;
  id = NOVAL;
  min_ = lb;
  max_ = ub;
  size_ = (ub-lb+1);
  state = ( BLIST | ILIST | RANGE );
  veqptr_ = NULL;
  varptr_ = NULL;
  domain_ = NULL;
  model = NULL;
}

BuildObject::BuildObject() 
{
  initVar(0,-1);
}

BuildObject::BuildObject( const int lo, const int up ) 
{
  initVar(lo, up);
}

BuildObject::BuildObject( const int* array, const int length ) 
{
  min_ = NOVAL;
  max_ = -NOVAL;
  for(int i=0; i<length; ++i) {
    if(min_ > array[i]) min_ = array[i];
    if(max_ < array[i]) max_ = array[i];
  }
  initVar(min_, max_);
  size_ = length;
  if( length < (max_ - min_ + 1) ) {
    domain_ = new BitSet(min_, max_, BitSet::empt);
    for(int i=0; i<length; ++i) domain_->insert( array[i] );
  } else state |= RANGE ;
}

BuildObject::BuildObject( const int* array, 
			  const int length, 
			  const int lo, 
			  const int up ) 
{
  initVar(lo, up);
  size_ = length;
  if( length < (max_ - min_ + 1) ) {
    domain_ = new BitSet(min_, max_, BitSet::empt);
    for(int i=0; i<length; ++i) domain_->insert( array[i] );
  } else state |= RANGE ;
}

BuildObject::BuildObject( const BitSet& d, 
			  const int length, 
			  const int lo, 
			  const int up ) 
{
  initVar(lo, up);
  size_ = length;
  if( length < (max_ - min_ + 1) ) {
    domain_ = new BitSet(min_, max_, BitSet::empt);
    domain_->unionWith( d );
  } else state |= RANGE ;
}

BuildObject* BuildObject::getBuildObject()
{
  return (veqptr_ ? veqptr_->getBuildObject() : this);
}

const BuildObject* BuildObject::getConstBuildObject() const
{

  //std::cout << veqptr_ << std::endl;

  return (veqptr_ ? veqptr_->getConstBuildObject() : this);
}

void BuildObject::build(SatSolver *s, const int *idx)
{
  if( !veqptr_ && isReferenced() ) {
    if( min_ < 0 || max_ > 1 ) {
      cerr << " cannot be converted to a sat problem: non-Boolean variables" << endl;
      exit(1);
    } else if( min_ == 1 ) {
      Vector<Literal> new_clause;
      new_clause.push( idx[id] );
      s->addClause( s->base, new_clause, s->stats.base_avg_size );
      s->addOriginalClause( new_clause );
      new_clause.clear();
    } else if( max_ == 0 ) {
      Vector<Literal> new_clause;
      new_clause.push( -idx[id] );
      s->addClause( s->base, new_clause, s->stats.base_avg_size );
      s->addOriginalClause( new_clause );
      new_clause.clear();
    }
  }
}

void BuildObject::build(Solver *s)
{
  if( !veqptr_ && isReferenced() ) {
    if( size_ == 1) {
      map<int,VariableInt*>::iterator iter = s->constants.find( (int)min_ );
      if( iter != s->constants.end() ) {
	varptr_ = iter->second;
      } else {// this constant as not yet been created
	varptr_ = new Constant( s, min_ );
	s->constants[(int)min_] = varptr_;
      }
    } else if( size_ == 2 && min_ == 0 && max_ == 1 ) {
      //cout << "bool" << endl;
      varptr_ = new VariableBool(s);
    } else if( size_ == (max_ - min_ + 1) ) {
      if( isRange() && size_ > 64 ) {
	//cout << "range" << endl;
	varptr_ = new VariableRange( s, min_, max_ );
      } else {
	if( isIList() ) {
	  //cout << "integer list" << endl;
	  varptr_ = new VariableBitList( s, min_, max_ );
	} else if( isBList() ) {
	  //cout << "integer list + bitset" << endl;
	  varptr_ = new VariableBitList( s, min_, max_ );
	} else {
	  //cout << "bitset" << endl;
	  varptr_ = new VariableBit( s, min_, max_ );
	}
      }
    } else {
      if( isIList() ) {
	//cout << "integer list" << endl;
	varptr_ = new VariableBitList( s, domain_, size_, min_, max_ );
      } else if( isBList() ) {
	//cout << "integer list + bitset" << endl;
	varptr_ = new VariableBitList( s, domain_, size_, min_, max_ );
      } else {
	//cout << "bitset" << endl;
	varptr_ = new VariableBit( s, domain_, size_, min_, max_ );
      }
    }
  }
}

// void BuildObject::releaseEvent( const int sb )
// {
//   int sn = size_;
//   if( model && sb != sn ) {
//     model->triggerEvent( id, 1 );
//     //model->triggerEvent( getVarId(), 1 );
    
// #ifdef _DEBUGMODEL
//      cout << "\t\ttrigger event for ";
//      print( cout );
//      cout << endl;
// #endif
     
//   }
// }


void BuildObject::releaseEvent( const int sb )
{
  int sn = size();
  if( model && sb != sn ) {
    model->triggerEvent( id, 1 );
    //model->triggerEvent( getVarId(), 1 );
    
#ifdef _DEBUGMODEL
     cout << "\t\ttrigger event for ";
     print( cout );
     cout << endl;
#endif
     
  }
}

int BuildObject::_size() const
{
  return size_;
}

int BuildObject::size() const
{
  return getConstBuildObject()->_size();
}

int BuildObject::_min() const
{
  return min_;
}

int BuildObject::min() const
{
  return getConstBuildObject()->_min();
}

int BuildObject::_max() const
{
  return max_;
}

int BuildObject::max() const
{
  return getConstBuildObject()->_max();
}

int BuildObject::_domsize() const
{
  return size_;
}

int BuildObject::domsize() const
{
  return getConstBuildObject()->_domsize();
}

int BuildObject::_value() const
{
  return min_;
}

int BuildObject::value() const
{
  //if( varptr_ )
  //return varptr_->getValue();
  return getConstBuildObject()->_value();
}

int BuildObject::getValue() const 
{
  const BuildObject *self = getConstBuildObject();

  if( self->varptr_ )
    return self->varptr_->getValue();
  else
    return self->_value();
}

int BuildObject::getSize() const 
{
  const BuildObject *self = getConstBuildObject();

  if( self->varptr_ )
    return self->varptr_->domsize();
  else
    return self->size_;
}

int BuildObject::getVariableId() const 
{
  const BuildObject *self = getConstBuildObject();

  if( self->varptr_ )
    return self->varptr_->id;
  else
    return -1;
}

int BuildObject::getMin() const 
{
  const BuildObject *self = getConstBuildObject();

  if( self->varptr_ )
    //return self->varptr_->getMin();
    return self->varptr_->min();
  else
    return self->_min();
}

int BuildObject::getMax() const 
{
  const BuildObject *self = getConstBuildObject();

  if( self->varptr_ )
    //return self->varptr_->getMax();
    return self->varptr_->max();
  else
    return self->_max();
}

bool BuildObject::isIn(const int v) const 
{
  const BuildObject *self = getConstBuildObject();

  if( self->varptr_ )
    return self->varptr_->contain(v);
  else
    return self->_contain(v);
}

int BuildObject::_first() const
{
  return min_;
}

int BuildObject::first() const
{
  return getConstBuildObject()->_first();
}

int BuildObject::_last() const
{
  return max_;
}

int BuildObject::last() const
{
  return getConstBuildObject()->_last();
}

int BuildObject::_minCapacity() const
{
  return min_;
}

int BuildObject::minCapacity() const
{
  return getConstBuildObject()->_minCapacity();
}

int BuildObject::_maxCapacity() const
{
  return max_;
}

int BuildObject::maxCapacity() const
{
  return getConstBuildObject()->_maxCapacity();
}

int BuildObject::_minPosAbs() const 
{
  if( max_ < 0 ) return NOVAL;
  if( min_ > 0 ) return min_;
  if( !domain_ || domain_->member( 0 ) ) return 0;
  return domain_->next( 0 );
}

int BuildObject::minPosAbs() const
{
  return getConstBuildObject()->_minPosAbs();
}

int BuildObject::_minNegAbs() const 
{
  if( min_ > 0 ) return NOVAL;
  if( max_ < 0 ) return max_;
  if( !domain_ || domain_->member( 0 ) ) return 0;
  return domain_->prev( 0 );
}

int BuildObject::minNegAbs() const
{
  return getConstBuildObject()->_minNegAbs();
}

bool BuildObject::_remove(const int v) 
{

  int sb = size_;

  if( min_ <= v && max_ >= v ) {
    if( domain_ != NULL ) {
      if(domain_->member(v)) {
	domain_->erase( v );
	--size_;
	if( v == min_ )
	  min_ = domain_->min();
	if( v == max_ )
	  max_ = domain_->max();
      }
    } else {
      if( v == min_ ) {
	--size_;
	++min_;
      } else if( v == max_ ) {
	--size_;
	--max_;
      } else {
	domain_ = new BitSet(min_, max_, BitSet::full);
	domain_->erase( v );
	--size_;
      }
    }
  }
 
  releaseEvent( sb );

  return size_ > 0;
}

bool BuildObject::remove(const int v) 
{
  int sb = size();
  bool res = getBuildObject()->_remove( v );
  releaseEvent( sb );
  return res;
}

bool BuildObject::_setDomain(const int v) 
{
  int sb = size_;

  if( domain_ != NULL ) {
    if( !domain_->member(v) )
      return false;
    else {
      delete domain_;
      domain_ = NULL;
    }
  } else if( min_ > v || v > max_ ) return false;

  size_ = 1;
  min_ = v;
  max_ = v;

  releaseEvent( sb );

  return true;
}

bool BuildObject::setDomain(const int v) 
{
  int sb = size();
  bool res = getBuildObject()->_setDomain( v );
  releaseEvent( sb );
  return res;
}

bool BuildObject::_setEqual(BuildObject *v) 
{
  bool rval = true;
  while( !v->parent.empty() ) {
    parent.push( v->parent.pop() );
  }
  refs += v->refs;
  if( !v->isRange() )
    unsetRange();
  rval = (_setMin( v->min_ ) &&
	  _setMax( v->max_ ));
  if( v->domain_ ) 
    rval = _setDomain( *(v->domain_) );
  v->veqptr_ = this;
  
  return rval;
}

bool BuildObject::setEqual(BuildObject *v) 
{
  //getVarId() < v->getVarId()
  if( getConstBuildObject() == v->getConstBuildObject() ) return true;
  if( getVarId() < v->getVarId() )
    return getBuildObject()->_setEqual( v->getBuildObject() );
  else
    return v->getBuildObject()->_setEqual( getBuildObject() );
}

bool BuildObject::_setMin(const int l) 
{
  int sb = size_;

  if( l > min_ ) { 
    if( domain_ != NULL ) {
      domain_->setMin(l);
      size_ = domain_->size();
      min_ = domain_->min();
    } else {
      size_ -= ( l - min_ );
      min_ = l;
    }
  }

  releaseEvent( sb );

  return size_ > 0;
}

bool BuildObject::setMin(const int l) 
{
  int sb = size();
  bool res = getBuildObject()->_setMin( l );
  releaseEvent( sb );
  return res;
}

bool BuildObject::_setMax(const int u)
{

  int sb = size_;

  if( u < max_ ) { 
    if( domain_ != NULL ) {
      domain_->setMax(u);
      size_ = domain_->size();
      max_ = domain_->max();
    } else {
      size_ -= ( max_ - u);
      max_ = u;
    }
  }

  releaseEvent( sb );

  return size_ > 0;
}

bool BuildObject::setMax(const int u) 
{
  int sb = size();
  bool res = getBuildObject()->_setMax( u );
  releaseEvent( sb );
  return res;
}

bool BuildObject::_setDomain(const BitSet& s) 
{

  int sb = size_;

  if( domain_ == NULL )  
    domain_ = new BitSet(min_, max_, BitSet::full);
  domain_->intersectWith( s );
  size_ = domain_->size();
  min_ = domain_->min();
  max_ = domain_->max();

  releaseEvent( sb );

  return size_ > 0;
}

bool BuildObject::setDomain(const BitSet& s) 
{
  int sb = size();
  bool res = getBuildObject()->_setDomain( s );
  releaseEvent( sb );
  return res;
}

bool BuildObject::_removeSet(const BitSet& s) 
{
  int sb = size_;

  if( domain_ == NULL ) 
    domain_ = new BitSet(min_, max_, BitSet::full);
  domain_->setminusWith( s );
  size_ = domain_->size();
  min_ = domain_->min();
  max_ = domain_->max();

  releaseEvent( sb );

  return size_ > 0;
}

bool BuildObject::removeSet(const BitSet& s) 
{
  int sb = size();
  bool res = getBuildObject()->_removeSet( s );
  releaseEvent( sb );
  return res;
}

bool BuildObject::_removeSet(const int* vals, const int n) 
{
  for(int i=0; i<n; ++i)
    if( !_remove(vals[i]) ) return false;
  return true;
}

bool BuildObject::removeSet(const int* vals, const int n) 
{
  int sb = size();
  bool res = getBuildObject()->_removeSet( vals, n );
  releaseEvent( sb );
  return res;
}

bool BuildObject::_removeRange(const int l, const int u) 
{
  
  int sb = size_;
 
  if( l <= min_ )
    setMin( u+1 );
  else if( u >= max_ )
    setMax( l-1 );
  else {
    if( domain_ == NULL ) 
      domain_ = new BitSet( min_, max_, BitSet::full );
    domain_->removeInterval( l, u );
    size_ = domain_->size();
  }

  releaseEvent( sb );

  return size_ > 0;
}

bool BuildObject::removeRange(const int l, const int u) 
{
  int sb = size();
  bool res = getBuildObject()->_removeRange(l,u);
  releaseEvent( sb );
  return res;
}

bool BuildObject::_contain(const int v) const 
{ 
  return (v >= min_ &&
	  v <= max_ &&
	  ( !domain_ || domain_->member(v) ));
}

bool BuildObject::contain(const int v) const
{
  return getConstBuildObject()->_contain( v );
}

bool BuildObject::_equal(const int v) const 
{ 
  return (v == min_ && v == max_ );
}

bool BuildObject::equal(const int v) const
{
  return getConstBuildObject()->_equal( v );
}

bool BuildObject::_intersect(const BuildObject* x) const
{
  if( x->min_ > max_ || x->max_ < min_ ) return false;
  if( domain_ ) {
    if( x->domain_ )
      return domain_->intersect( x->domain_ );
    else {
      BitsetIterator bit( *domain_ );
      do if( bit >= x->min_ && bit <= x->max_ ) return true;
      while( ++bit );
      return false;
    }
  } else if( x->domain_ ) {
    BitsetIterator bit( *(x->domain_) );
    do if( bit >= min_ && bit <= max_ ) return true;
      while( ++bit );
    return false;
  } 
  return true;
}

bool BuildObject::intersect(const BuildObject* x) const
{
  return getConstBuildObject()->_intersect( x->getConstBuildObject() );
}

bool BuildObject::_included(const BuildObject* x) const
{
  if( x->min_ > min_ || x->max_ < max_ ) return false;
  if( x->domain_ ) {
    if( domain_ )
      return domain_->included( x->domain_ );
    else
      return false;
  }  
  return true;
}

bool BuildObject::included(const BuildObject* x) const
{
  return getConstBuildObject()->_included( x->getConstBuildObject() );
}

int BuildObject::_getNext(const int v) const 
{ 
  if( domain_ != NULL )
    return ( domain_->next(v) ); 
  else return ( v < max_ ? v+1 : NOVAL );
}

int BuildObject::getNext(const int v) const
{
  return getConstBuildObject()->_getNext( v );
}

bool BuildObject::_good(const int v) const 
{ 
  return ( v != NOVAL ); 
}

bool BuildObject::good(const int v) const
{
  return getConstBuildObject()->_good(v);
}

bool BuildObject::_setNext(int& v) const 
{ 
  if( domain_ != NULL )
    return ( (v = domain_->next(v)) != NOVAL ); 
  else return ( ++v <= max_ );
}

bool BuildObject::setNext(int& v) const
{
  return getConstBuildObject()->_setNext(v);
}

bool BuildObject::_isInterval() const
{
  return (domain_ == NULL);
}

bool BuildObject::_isGround() const
{
  return (size_ <= 1);
}

bool BuildObject::isInterval() const
{
  return getConstBuildObject()->_isInterval();
}

bool BuildObject::isGround() const
{
  return getConstBuildObject()->_isGround();
}

int BuildObject::getVarId() const { return getConstBuildObject()->id; } 

int BuildObject::getState() const
{
  return getConstBuildObject()->state;
}
void BuildObject::setState( const int s ) 
{
  getBuildObject()->state = s;
}

// bool BuildObject::isRange() const { return (state & RANGE); }
// void BuildObject::setRange() { if( !(state & RANGE) ) state |= RANGE; }
// void BuildObject::unsetRange() { if(state & RANGE) state &= (0xffffffff ^ RANGE); }

// bool BuildObject::isBList() const { return (state & BLIST); }
// void BuildObject::setBList() { if( !(state & BLIST) ) state |= BLIST; }
// void BuildObject::unsetBList() { if( state & BLIST ) state &= (0xffffffff ^ BLIST); }

// bool BuildObject::isIList() const { return (state & ILIST); }
// void BuildObject::setIList() { if( !(state & ILIST) ) state |= ILIST; }
// void BuildObject::unsetIList() { if( state & ILIST ) state &= (0xffffffff ^ ILIST); }

// bool BuildObject::isRel() const { return false; }
// void BuildObject::unsetRel() { }

// bool BuildObject::isNeq() const { return (state & NEQ); }
// void BuildObject::setNeq() { if( !(state & NEQ) ) state |= NEQ; }
// void BuildObject::unsetNeq() { if( state & NEQ ) state &= (0xffffffff ^ NEQ); }

bool BuildObject::isRange() const { return (getState() & RANGE); }
void BuildObject::setRange() { int s = getState(); if( !(s & RANGE) ) setState( s | RANGE ); }
void BuildObject::unsetRange() { int s = getState(); if(s & RANGE) setState( s & (0xffffffff ^ RANGE) ); }

bool BuildObject::isBList() const { return (getState() & BLIST); }
void BuildObject::setBList() { int s = getState(); if( !(s & BLIST) ) setState ( s | BLIST ); }
void BuildObject::unsetBList() { int s = getState(); if( s & BLIST ) setState( s & (0xffffffff ^ BLIST) ); }

bool BuildObject::isIList() const { return (getState() & ILIST); }
void BuildObject::setIList() { int s = getState(); if( !(s & ILIST) ) setState( s | ILIST ); }
void BuildObject::unsetIList() { int s = getState(); if( s & ILIST ) setState( s & (0xffffffff ^ ILIST) ); }

bool BuildObject::isRel() const { return false; }
void BuildObject::unsetRel() { }

bool BuildObject::isNeq() const { return (getState() & NEQ); }
void BuildObject::setNeq() { int s = getState(); if( !(s & NEQ) ) setState( s | NEQ ); }
void BuildObject::unsetNeq() { int s = getState(); if( s & NEQ ) setState( s & (0xffffffff ^ NEQ) ); }

bool BuildObject::isVirtual() const { return (getState() & VIRTUAL); }
void BuildObject::setVirtual() { int s = getState(); if( !(s & VIRTUAL) ) setState( s | VIRTUAL ); }
void BuildObject::unsetVirtual() { int s = getState(); if( s & VIRTUAL ) setState( s & (0xffffffff ^ VIRTUAL) ); }


bool BuildObject::isTop() const { return (isReferenced()-isRel() <= 0); }


void BuildObject::_reference() { ++refs; }
void BuildObject::reference() { getBuildObject()->_reference(); }
void BuildObject::_deReference() { --refs; }
void BuildObject::deReference() { getBuildObject()->_deReference(); }
int BuildObject::_isReferenced() const { return refs; }
int BuildObject::isReferenced() const { return getConstBuildObject()->_isReferenced(); }

bool BuildObject::isVar()
{
  return (isReferenced() && !veqptr_);
}

bool BuildObject::isBoolean()
{
  return (size() == 2 && min() == 0 && max() == 1);
}

bool BuildObject::isFDomain()
{
  return (size() > 1 && !isBoolean() && !isInterval());
}

bool BuildObject::isInterval()
{
  return (size() == (max() - min() + 1) && isRange());
}

bool BuildObject::isConstant()
{
  return (size() == 1);
}

bool BuildObject::isSearchable()
{
  
  //std::cout << "<" << (size_ > 1) << " " << (isVirtual()) << ">=" << ();

  return (size_ > 1 // && veqptr_ == NULL
	  && !isVirtual()
	  );
}

bool BuildObject::isConstraint()
{
  return false;
}

VariableInt* BuildObject::_getVariable()
{
  return varptr_;
}

VariableInt* BuildObject::getVariable()
{
  return getBuildObject()->_getVariable();
}


void BuildObject::printshort(std::ostream& o) const   
{
  if( veqptr_ ) {
    o << "*";
    veqptr_->BuildObject::printshort( o );
  } else if( size_ == 1 ) {
    if( getVarId() < NOVAL/2 )
      o << "c" << getVarId() << "(" << min_ << ")";
    else 
      o << min_;
  } else {
    
    if( size_ == 2 && min_ == 0 && max_ == 1 ) {
      o << "b";
    } else if( size_ == (max_ - min_ + 1) ) {
      if( isRange() && size_ > 64 ) {
	o << "r";
      } else {
	if( isIList() ) {
	  o << "z";
	} else if( isBList() ) {
	  o << "y";
	} else {
	  o << "x";
	}
      }
    } else {
      if( isIList() ) {
	o << "z";
      } else if( isBList() ) {
	o << "y";
      } else {
	o << "x";
      }
    }
    if( getVarId() < NOVAL/2 )
      o //<< "x" 
	<< getVarId();
  } 
}

void BuildObject::print_python() 
{
  VariableInt *x = getVariable();
  if(x) {
    if(x->getType() == VariableInt::VIRTUAL)
      std::cout << "v" << ((VariableVirtual*)x)->reference->id ;
    else
      std::cout << "x" << x->id ;
  }else
    std::cout << "c" << getMin() ;
  //this->getConstBuildObject()->BuildObject::printshort(std::cout);
}

void BuildObject::print(std::ostream& o) const 
{
  if( varptr_ )
    varptr_->print(o);
  else {
    BuildObject::printshort( o );
    const BuildObject *self = getConstBuildObject();
    if( self->domain_ != NULL ) {
      o << "in" ;
      self->domain_->print( o );
    } else if( self->size_ != 1 ) {
      o << "in{" << self->min_ << ",.." << self->max_ << "}";
    } 
  }
}

string BuildObject::toString() const 
{
  string s = ("x" + int2string(id) + " [" + int2string(getMin()) + ".." + int2string(getMax()) + "]"); 
  return s;
}

string BuildObject::xmlPred( int& idx, int level, int *sub ) const 
{
  string s = ("X" + int2string(idx++)); 
  return s;
}

string BuildObject::xmlCon( BitSet& scope, int *index ) const 
{
  int i=getVarId();
  if( size() > 1 ) scope.insert( index[i] );
  string s = (size() > 1 ? ("V" + int2string(index[i])) : (int2string(min_)) ); 
  return s;
}

string BuildObject::xmlDom( ) const 
{
  
  string s;
  if( domain_ != NULL ) {
    if( getConstBuildObject()->varptr_ )
      s += ((VariableDomain*)(getConstBuildObject()->varptr_))->values.toString();
    else
      s += getConstBuildObject()->domain_->toString();
  } else if( size() > 2 ) {
    s += ((int2string(min())) + ".." + (int2string(max()))); 
  } else if( size() == 2 ) {
    s += ((int2string(min())) + " " + (int2string(max())));
  } else {
    s += (int2string(min()));
  }
  
  return s;
}

void BuildObject::printstate( std::ostream& o ) const 
{
}

BuildObject::~BuildObject() 
{
  if( domain_ ) {
    domain_->neg_words = 0;
    domain_->pos_words = 0;
    domain_->table = NULL;
    delete domain_;
  }
}


bool BuildObject::operator==( BuildObject& x ) const
{
  return (getConstBuildObject() == x.getConstBuildObject());
}

/**********************************************
 * Constraint BuildObject
 **********************************************/ 

BuildObjectConstraint::BuildObjectConstraint() {}

BuildObjectConstraint::~BuildObjectConstraint() {}

void BuildObjectConstraint::add( CSP *m, const int l, 
				 BuildObjectPredicate* pred ) 
{
  int i, n=pred->arity;
  for(i=0; i<n; ++i) {
    pred->scope[i]->add( m, l+1 );
  }
  //++pred->refs;
  pred->reference();

//   cout << "ADDING ";
//   pred->print( cout );
//   cout << endl;

}

void BuildObjectConstraint::build(SatSolver *, const int *, BuildObjectPredicate *)
{
  cerr << "Cannot convert constraint to clausal form" << endl;
  exit(1);
}

void BuildObjectConstraint::self_build( Solver *s, VariableInt **tmp, BuildObjectPredicate *pred )
{
  pred->BuildObject::build( s );
}

void BuildObjectConstraint::print(std::ostream&, const BuildObjectPredicate*) const {}

string BuildObjectConstraint::toString(const BuildObjectPredicate*) const 
{
  string s;
  return s; 
}

string BuildObjectConstraint::xmlPred(int&, int, const BuildObjectPredicate*, int*) const
{
  string s; 
  return s;
}


/**********************************************
 * Predicate BuildObject
 **********************************************/
BuildObjectPredicate::BuildObjectPredicate(BuildObject **x, const int n, 
					   const int lb, const int ub, 
					   BuildObjectConstraint *r, int *p)
  : BuildObject(lb, ub)
{
  arity = n;
  params = p;
  scope = x;
  relation = r;
  scope[n] = this;
}

BuildObjectPredicate::~BuildObjectPredicate() 
{
  delete [] scope;
  delete [] params;
}

// void BuildObject::gather_variables()
// {
// }

bool BuildObjectPredicate::isRel() const 
{
  return relation;
}

void BuildObjectPredicate::unsetRel() 
{
  relation = NULL;
  --refs;
}

void BuildObjectPredicate::add( CSP *m, const int l )
{
  for(int i=0; i<arity; ++i) {
    scope[i]->parent.push( this );
  }
  relation->add( m, l, this );
  BuildObject::add( m, l );
}

void BuildObjectPredicate::close()
{
  relation->close( this );
}

int BuildObjectPredicate::propagateUpward( ) 
{
  return relation->propagateUpward( this );
}

int BuildObjectPredicate::propagateDownward( ) 
{
  return relation->propagateDownward( this );
}

void BuildObjectPredicate::build(SatSolver *s, const int *idx) 
{
  if( isReferenced() ) 
    BuildObject::build( s, idx );
  
  if( isRel() ) 
    relation->build( s, idx, this );
}

void BuildObjectPredicate::build(Solver *s) 
{

  if( !isRel() ) {
    if( isReferenced() )
      BuildObject::build( s );
  } else {
    int n=arity+(isReferenced() > 0);
    VariableInt *tmpvars[n];
    for(int i=0; i<n; ++i) 
      tmpvars[i] = scope[i]->getVariable();
    relation->self_build( s, tmpvars, this );
    if( isReferenced() > 0 )
      tmpvars[n-1] = scope[n-1]->getVariable();
    relation->build( s, tmpvars, this );
  }
  
}
  
string BuildObjectPredicate::xmlCon( BitSet& scp, int *index ) const
{
  string s;
  if( arity > 3 )
    s = BuildObject::xmlCon( scp, index );
  else {
    int n=arity-1;
    for(int i=0; i<n; ++i)
      s += (scope[i]->xmlCon( scp, index ) + " ");
    s += scope[n]->xmlCon( scp, index );
  }
  return s;
}

void BuildObjectPredicate::printshort(std::ostream& o) const
{
  if( isRel() )
    relation->print( o, this );
  else if( isReferenced() )
    BuildObject::printshort( o );
}
 
void BuildObjectPredicate::print(std::ostream& o) const
{
  if( getConstBuildObject()->varptr_ ) {
    getConstBuildObject()->varptr_->print( o );
    o << " = "; 
  } else {
    if( isReferenced() ) {
      BuildObject::print( o );
      o << " = "; 
    } 
  }
  if( isRel() ) {
    relation->print( o, this );
  }
}

string BuildObjectPredicate::toString() const 
{
  return relation->toString( this );
}
string BuildObjectPredicate::xmlPred( int& idx, int level, int *sub ) const
{
  return relation->xmlPred( idx, level, this, sub );
}


/**********************************************
 * Variable Wrapper
 **********************************************/

Variable::Variable() { var_ptr_=NULL; }

Variable::~Variable() 
{}

int Variable::value()
{
  int v = NOVAL;
//   if( var_ptr_->getVariable() != NULL ) {
//     v = var_ptr_->getVariable()->value();
//   } else if( var_ptr_->getConstBuildObject() != NULL ) {
//     v = var_ptr_->getConstBuildObject()->value();
//   } 
  return var_ptr_->getValue();
  return v;
}

int Variable::size()
{
  int s = NOVAL;
  if( var_ptr_->getVariable() != NULL )
    s = var_ptr_->getVariable()->domsize();
  else s = var_ptr_->size();
  return s;
}

int Variable::min()
{
  int m = var_ptr_->min();
  if( var_ptr_->getVariable() != NULL )
    m = var_ptr_->getVariable()->min();
  else m = var_ptr_->min();
  return m;
}

int Variable::max()
{
  int m = var_ptr_->max();
  if( var_ptr_->getVariable() != NULL )
    m = var_ptr_->getVariable()->max();
  else m = var_ptr_->max();
  return m;
}

Variable& Variable::operator= ( const Variable& x )
{
  var_ptr_ = x.var_ptr_;
  return *this;
}

VarArray Variable::operator,( Variable x )
{
  VarArray table(2);
  table[0] = *this;
  table[1] = x;
  return table;
}

Variable::Variable( const int nvals ) 
{
  var_ptr_ = CSP::_Variable( nvals );
}

Variable::Variable( const Variable& v )
{
  var_ptr_ = v.var_ptr_;
}

Variable::Variable( BuildObject *v )
{
  var_ptr_ = v;
}

Variable::Variable( const int lb, const int ub ) 
{
  var_ptr_ = CSP::_Variable( lb, ub );
}

Variable::Variable( const int* array, const int length ) 
{
  var_ptr_ = CSP::_Variable( array, length );
}

Variable::Variable( const int* array, const int length, 
		    const int lb, const int ub ) 
{
  var_ptr_ = CSP::_Variable( array, length, lb, ub );
}

Variable::Variable( const BitSet& d ) 
{
  var_ptr_ = CSP::_Variable( d, d.size(), d.min(), d.max() );
}

Variable::Variable( const BitSet& d, const int length, 
		    const int lb, const int ub ) 
{
  var_ptr_ = CSP::_Variable( d, length, lb, ub );
}

VariableInt *Variable::getVariable() 
{
  if(var_ptr_) return var_ptr_->getVariable();
  else return NULL;
}

void Variable::printshort(std::ostream& o) const 
{
  if( var_ptr_->getVariable() != NULL )
    var_ptr_->getVariable()->printshort( o );
  else var_ptr_->printshort( o );
}

void Variable::print(std::ostream& o) const 
{
  if( var_ptr_->getVariable() != NULL )
    var_ptr_->getVariable()->print( o );
  else var_ptr_->print( o );
}


/**********************************************
 * VarArray Wrapper
 **********************************************/

VarArray::VarArray( ) {

  //cout << "create an empty var array" << endl;

}

VarArray::VarArray( const int nvars )
{
  array_.resize( nvars );
}

VarArray::VarArray( const int nvars, const int nvals )
{
#ifdef _NARRAY
  array_.resize( nvars );
  int i=nvars;
  while( i-- ) {
    array_[i] = CSP::_Variable( nvals );
  }
#else
  array_.resize( nvars );
  tmp_.resize( nvars );
  int i=nvars;
  while( i-- ) {
//     tmp_[i] = CSP::_Variable( nvals );
//     array_[i] = Variable( tmp_[i] );
    array_[i] = Variable( nvals );
  }
#endif
}

VarArray::VarArray( const int nvars, const int lo, const int up ) 
{

  //cout << "create a var array" << endl;

#ifdef _NARRAY
  array_.resize( nvars );
  int i=nvars;
  while( i-- ) {
    array_[i] = CSP::_Variable( lo, up );
  }
#else
  array_.resize( nvars );
  int i=nvars;
  while( i-- ) {
//     tmp_[i] = CSP::_Variable( lo, up );
//     array_[i] = Variable( tmp_[i] );
    array_[i] = Variable( lo, up );
  }
#endif
}

VarArray::VarArray( const int nvars, const int* array, const int length ) 
{
#ifdef _NARRAY
  array_.resize( nvars );
  int i=nvars;
  while( i-- ) {
    array_[i] = CSP::_Variable( array, length );
  }
#else
  array_.resize( nvars );
  int i=nvars;
  while( i-- ) {
//     tmp_[i] = CSP::_Variable( array, length );
//     array_[i] = Variable( tmp_[i] );
    array_[i] = Variable( array, length );
  }
#endif
}

VarArray::VarArray( const int nvars, const int* array, const int length, const int lo, const int up )
{
#ifdef _NARRAY
  array_.resize( nvars );
  int i=nvars;
  while( i-- ) {
    array_[i] = CSP::_Variable( array, length, lo, up );
  }
#else
  array_.resize( nvars );
  int i=nvars;
  while( i-- ) {
//     tmp_[i] = CSP::_Variable( array, length, lo, up );
//     array_[i] = Variable( tmp_[i] );
    array_[i] = Variable( array, length, lo, up );
  }
#endif
}

VarArray::~VarArray() 
{
#ifdef _NARRAY
#else
//   int i=array_.size();
//   while( i-- ) {
//     if( array_[i].var_ptr_->id == NOVAL )
//       delete array_[i].var_ptr_;
//   }
#endif
  //cout << "delete var array" << endl;
}

BuildObject **VarArray::getArgs() 
{
#ifdef _NARRAY
  return array_.stack_;
#else
//   int n = array_.size();
//   tmp_.clear();
//   tmp_.resize( n );
//   for( int i=0; i<n; ++i ) 
//     tmp_[i] = array_[i].var_ptr_;
//     //tmp_.push( array_[i].var_ptr_ );
//   return &(tmp_[0]);//tmp_.stack_;
// //   int n = array_.size();
// //   BuildObject **tmp = new BuildObject*[n+1];
// //   while( n-- ) tmp[n] = array_[n].var_ptr_;  
// //   return tmp;
  int n = array_.size();
  BuildObject **args = new BuildObject*[n+1];
  for( int i=0; i<n; ++i ) 
    args[i] = array_[i].var_ptr_;
  return args;
#endif
}


VariableInt **VarArray::getVarArray()
{
  int n = array_.size();
  VariableInt **x_array = new VariableInt*[n];
  while( n-- )
    x_array[n] = array_[n].var_ptr_->getVariable();
  return x_array;
}

// VariableInt *Variable::getVariable()
// {
//   return var_ptr_->varptr_;
// }


void VarArray::add( const int k )
{
#ifdef _NARRAY
  array_.push_back( CSP::_Variable(k, k) );
#else
//   tmp_.push_back( CSP::_Variable(k, k) );
//   array_.push_back( Variable( tmp_.back() ) );
  array_.push_back( Variable(k, k) );
#endif
}

void VarArray::add( Variable x )
{
#ifdef _NARRAY
  array_.push_back( x.var_ptr_ );
#else
  //tmp_.push_back( x.var_ptr_ );
  array_.push_back( x );
#endif
}

void VarArray::add( VarArray& a )
{
#ifdef _NARRAY
  int i=0, n=a.size();
  for(i=0; i<n; ++i)
    add( a.array_[i] );
#else
  int i=0, n=a.size();
  for(i=0; i<n; ++i) {
    //tmp_.push_back( a.tmp_[i] );
    add( a[i] );
  }
#endif
}

void VarArray::build( Solver *s ) 
{
#ifdef _NARRAY
  int n=array_.size;
  for(int i=0; i<n; ++i)
    array_[i]->build( s );
#else
  int n=array_.size();
  for(int i=0; i<n; ++i)
    array_[i].var_ptr_->build( s );
  //tmp_[i]->build( s );
#endif
}  

void VarArray::resize( const int n )
{
#ifdef _NARRAY
  array_.resize( n );
#else
  array_.resize( n );
  tmp_.resize( n );
#endif
}

void VarArray::clear()
{
#ifdef _NARRAY
  array_.clear();
#else
  array_.clear();
  //tmp_.clear();
#endif
}

int VarArray::size() const
{
  return array_.size();
}


/*!@name Operators*/
Variable& VarArray::operator[]( const int i ) 
{
#ifdef _NARRAY
  buffer.var_ptr_ = array_[i];
  return buffer; 
#else
  return array_[i];
#endif
}

VarArray VarArray::operator,( Variable x )
{
#ifdef _NARRAY
  int i=size();
  VarArray table(i+1);
  table.array_[i] = x.var_ptr_;
  while( i-- ) 
    table.array_[i] = array_[i];
  return table;
#else
  int i=size();
  VarArray table(i+1);
  table[i] = x;
  //table.tmp_[i] = x.var_ptr_;
  while( i-- ) {
    table[i] = array_[i];
    //table.tmp_[i] = tmp_[i];
  }
  return table;
#endif
}

VarArray VarArray::operator,( VarArray& x )
{
#ifdef _NARRAY
  int i=size(), j=x.size();
  VarArray table(i+j);
  while( j-- )
    table.array_[i+j] = x.array_[j];
  while( i-- )
    table.array_[i] = array_[i];
  return table;
#else
  int i=size(), j=x.size();
  VarArray table(i+j);
  while( j-- ) {
    //table.tmp_[i+j] = tmp_[i+j];
    table[i+j] = x.array_[j];
  }
  while( i-- ) {
    //table.tmp_[i] = tmp_[i];
    table[i] = array_[i];
  }
  return table;
#endif
}

VarArray& VarArray::operator=( VarArray& x )
{
#ifdef _NARRAY
  array_ = x.array_;
#else
  array_ = x.array_;
  //tmp_ = x.tmp_;
#endif

  return *this;
}

void VarArray::printshort(std::ostream& o) const 
{
  int i=0, n=size()-1;
  if( n >= 0 ) {
    o << "[";
    o.flush();
    for(i=0; i<n; ++i)
      {
#ifdef _NARRAY
	array_[i]->printshort( o );
#else
	array_[i].printshort( o );
#endif
	o << ", ";
	o.flush();
      }
#ifdef _NARRAY
	array_[n]->printshort( o );
#else
	array_[n].printshort( o );
#endif
    o << "]";
    o.flush();
  }
}
  
void VarArray::print(std::ostream& o) const 
{
  int i=0, n=size()-1;
  if( n >= 0 ) {
    o << "[";
    o.flush();
    for(i=0; i<n; ++i)
      {
#ifdef _NARRAY
	array_[i]->print( o );
#else
	array_[i].print( o );
#endif
	o << ", ";
	o.flush();
      }
#ifdef _NARRAY
	array_[n]->print( o );
#else
	array_[n].print( o );
#endif
    o << "]";
    o.flush();
  }
}


/**********************************************
 * Matrix Wrapper
 **********************************************/

Matrix::Matrix( const int nr, const int nc ) 
{
  nCols=nc;
  nRows=nr;
  array_.resize( nr*nc );
  tmp_.resize( std::max(nr, nc) );
}

Matrix::Matrix( const int nr, const int nc, const int nv ) 
{
  nCols=nc;
  nRows=nr;
  array_.resize( nr*nc );
  tmp_.resize( std::max(nr, nc) );
  int i=nc*nr;
  while( i-- ) 
    array_[i] = CSP::_Variable( nv );
}

Matrix::Matrix( const int nr, const int nc, const int lo, const int up ) 
{
  nCols=nc;
  nRows=nr;
  array_.resize( nr*nc );
  tmp_.resize( std::max(nr, nc) );
  int i=nc*nr;
  while( i-- ) 
    array_[i] = CSP::_Variable( lo, up );
}

Matrix::Matrix( const int nr, const int nc, const int* matrix, const int length ) 
{
  nCols=nc;
  nRows=nr;
  array_.resize( nr*nc );
  tmp_.resize( std::max(nr, nc) );
  int i=nc*nr;
  while( i-- ) 
    array_[i] = CSP::_Variable( matrix, length );
}

Matrix::Matrix( const int nr, const int nc, const int* matrix, const int length, const int lo, const int up ) 
{
  nCols=nc;
  nRows=nr;
  array_.resize( nr*nc );
  tmp_.resize( std::max(nr, nc) );
  int i=nc*nr;
  while( i-- ) 
    array_[i] = CSP::_Variable( matrix, length, lo, up );
}


Matrix::~Matrix() {}

Variable Matrix::cell( const int i, const int j ) 
{
  tmp_[0] = Variable( array_[i*nCols+j] );
  return tmp_[0];
}

VarArray& Matrix::row( const int i ) 
{
  tmp_.resize( nCols );
  for( int j=0; j<nCols; ++j )
    tmp_[j] = Variable( array_[i*nCols+j] );
  return tmp_;
}
VarArray& Matrix::column( const int j ) 
{
  tmp_.resize( nRows );
  for( int i=0; i<nRows; ++i )
    tmp_[i] = Variable( array_[i*nCols+j] );
  return tmp_;
}

void Matrix::printshort(std::ostream& o) const {}
void Matrix::print(std::ostream& o) const 
{
  for( int i=0; i<nRows; ++i ) {
    for( int j=0; j<nCols; ++j ) {
      array_[i*nCols+j]->print( cout );
      cout << " ";
    }
    cout << endl;
  }
}


/**********************************************
 *  Negation Predicate BuildObject
 **********************************************/ 
void BuildObjectNegation::build( Solver *s, VariableInt **tmp, BuildObjectPredicate *pred )
{
  new PredicateNegation( s, tmp );
} 

//int getMin( BuildObjectPredicate *pred ) const ;
//int getMax( BuildObjectPredicate *pred ) const ;

int BuildObjectNegation::propagateUpward( BuildObjectPredicate *pred ) const 
{
  return ( pred->setMin( -(pred->scope[0]->max()) ) && pred->setMax( -(pred->scope[0]->min()) ) );
}

int BuildObjectNegation::propagateDownward( BuildObjectPredicate *pred ) const 
{  
  return ( pred->scope[0]->setMin( -(pred->max()) ) && pred->scope[0]->setMax( -(pred->min()) ) );
}

void BuildObjectNegation::close( BuildObjectPredicate *pred )
{  
  if( pred->isConstant() ) {
    pred->unsetRel();
    pred->scope[0]->deReference();
  } else {
    pred->model->unsetSat();
    pred->scope[0]->unsetBList();
    pred->scope[0]->unsetIList();
    pred->unsetBList();
    pred->unsetIList();
  }
}

void BuildObjectNegation::print(std::ostream& o, const BuildObjectPredicate *pred) const
{
  o << "-";
  pred->scope[0]->print( o );
}

string BuildObjectNegation::toString(const BuildObjectPredicate* pred) const 
{
  string x1 = pred->scope[0]->toString();
  string s = ( "-" + x1 );
  return s;
}

string BuildObjectNegation::xmlPred(int& idx, int level, const BuildObjectPredicate *pred, int *sub) const 
{
  string x1 = pred->scope[0]->xmlPred(idx, level+1, sub);
  string s = ( "neg(" + x1 + ")" );
  return s;
}


/**********************************************
 *  Not Predicate BuildObject
 **********************************************/ 
void BuildObjectNot::build( Solver *s, VariableInt **tmp, BuildObjectPredicate *pred ) 
{
  CSP *model = pred->model;
  if( model->satCompatible ) {
    if( !model->clauseBase )
      model->clauseBase = new ConstraintClauseBase(s);
    Vector<Literal> new_clause;
    for(int i=0; i<2; ++i) {
      for(int j=0; j<2; ++j) 
	new_clause.push( (i == j ? 1 : -1) * (tmp[j]->id+1) );
      model->clauseBase->addClause( new_clause );
      new_clause.clear();
    }
  } else {
    if( !pred->isConstant() )
      new PredicateNot( s, tmp );
  }
}

void BuildObjectNot::build( SatSolver *s, const int *idx, BuildObjectPredicate *pred )
{
  Vector<Literal> new_clause;
  for(int i=0; i<2; ++i) {
    for(int j=0; j<2; ++j) 
      new_clause.push( (i == j ? 1 : -1) * idx[pred->scope[j]->getVarId()] );
    s->addClause( s->base, new_clause, s->stats.base_avg_size );
    s->addOriginalClause( new_clause );
    new_clause.clear();
  }
} 

int BuildObjectNot::propagateUpward( BuildObjectPredicate *pred ) const 
{
  BuildObject *x_ = pred->scope[0];
  int ub = (x_->contain(0));
  int lb = (ub && x_->size() == 1);
  int consistent = ( pred->setMin( lb ) && pred->setMax( ub ) );
  if( consistent && pred->isTop() ) 
    consistent = pred->setDomain(1);
  return consistent;
}

int BuildObjectNot::propagateDownward( BuildObjectPredicate *pred ) const 
{
  int consistent  = true;
  BuildObject *x_ = pred->scope[0];
  if( pred->equal(1) )
    consistent = x_->setDomain(0);
  if( pred->equal(0) )
    consistent = x_->remove(0);
  return consistent;
}

void BuildObjectNot::close( BuildObjectPredicate *pred ) 
{  
  if( pred->isConstant() ) {
    pred->unsetRel();
    pred->scope[0]->deReference();
  } else {
    pred->scope[0]->unsetBList();
    pred->scope[0]->unsetIList();
    pred->unsetBList();
    pred->unsetIList();
  }
}

void BuildObjectNot::print(std::ostream& o, const BuildObjectPredicate *pred) const
{
  o << "!";
  pred->scope[0]->print( o );
}

string BuildObjectNot::toString(const BuildObjectPredicate* pred) const 
{
  string x1 = pred->scope[0]->toString();
  string s = ( "not " + x1 );
  return s;
}

string BuildObjectNot::xmlPred(int& idx, int level, const BuildObjectPredicate *pred, int *sub) const 
{
  string x1 = pred->scope[0]->xmlPred(idx, level+1, sub);
  string s = ( "not(" + x1 + ")" );
  return s;
}


/**********************************************
 *  Abs Predicate BuildObject
 **********************************************/ 
void BuildObjectAbs::build( Solver *s, VariableInt **tmp, BuildObjectPredicate *pred ) 
{
  new PredicateAbs( s, tmp );
}

int BuildObjectAbs::propagateUpward( BuildObjectPredicate *pred ) const 
{
  BuildObject *x = pred->scope[0];
  int lb, ub, aux;     
  if( x->min() >= 0 ) {
    lb = x->min();
    ub = x->max();
  } else if( x->max() <= 0 ) {
    ub = -1*x->min();
    lb = -1*x->max();
  } else {
    aux = -1*(x->minNegAbs());
    lb = x->minPosAbs();
    if( aux < lb ) lb = aux;    
    aux = -1*(x->min());
    ub = x->max();
    if( aux > ub ) ub = aux;
  }

  int consistent = ( pred->setMin( lb ) && pred->setMax( ub ) );
  if( consistent && pred->isTop() ) 
    consistent = pred->setDomain(1);
  return consistent;
}

int BuildObjectAbs::propagateDownward( BuildObjectPredicate *pred ) const 
{
  BuildObject *x = pred->scope[0];
  int mn = pred->max();
  int mx = pred->min()-1;  
  return (x->setMin(-mn) && x->setMax(mn) && (-mx >= mx  || x->removeRange(-mx, mx)));
}

void BuildObjectAbs::close( BuildObjectPredicate *pred ) 
{  
  if( pred->isConstant() ) {
    pred->unsetRel();
    pred->scope[0]->deReference();
  } else {
    pred->model->unsetSat();
    pred->scope[0]->unsetRange();
    pred->scope[0]->unsetBList();
    pred->scope[0]->unsetIList();
    pred->unsetBList();
    pred->unsetIList();
  }
}

void BuildObjectAbs::print(std::ostream& o, const BuildObjectPredicate *pred) const
{
  o << "|";
  pred->scope[0]->print( o );
  o << "|";
}

string BuildObjectAbs::toString(const BuildObjectPredicate* pred) const 
{
  string x1 = pred->scope[0]->toString();
  string s = ( "abs(" + x1 + ")" );
  return s;
}

string BuildObjectAbs::xmlPred(int& idx, int level, const BuildObjectPredicate *pred, int *sub) const 
{
  string x1 = pred->scope[0]->xmlPred(idx, level+1, sub);
  string s = ( "abs(" + x1 + ")" );
  return s;
}

/**********************************************
 * And Constraint BuildObject
 **********************************************/
void BuildObjectAnd::build( Solver *s, VariableInt **tmp, BuildObjectPredicate *pred ) 
{
  if( !pred->isConstant() || pred->equal(0) )
    new PredicateAnd( s, tmp );
}

void BuildObjectAnd::build( SatSolver *s, const int *idx, BuildObjectPredicate *pred )
{
  Vector<Literal> new_clause;

  if( !pred->isConstant()  ) {

    new_clause.push( -(idx[pred->scope[0]->getVarId()]) );
    new_clause.push( -(idx[pred->scope[1]->getVarId()]) );
    new_clause.push( idx[pred->scope[2]->getVarId()] );
    s->addClause( s->base, new_clause, s->stats.base_avg_size );
    s->addOriginalClause( new_clause );
    new_clause.clear();

    new_clause.push( idx[pred->scope[0]->getVarId()] );
    new_clause.push( -(idx[pred->scope[2]->getVarId()]) );
    s->addClause( s->base, new_clause, s->stats.base_avg_size );
    s->addOriginalClause( new_clause );
    new_clause.clear();

    new_clause.push( idx[pred->scope[1]->getVarId()] );
    new_clause.push( -(idx[pred->scope[2]->getVarId()]) );
    s->addClause( s->base, new_clause, s->stats.base_avg_size );
    s->addOriginalClause( new_clause );

  } else if( pred->equal(0) ) {

    new_clause.push( -(idx[pred->scope[0]->getVarId()]) );
    new_clause.push( -(idx[pred->scope[1]->getVarId()]) );
    s->addClause( s->base, new_clause, s->stats.base_avg_size );
    s->addOriginalClause( new_clause );

  }
} 

int BuildObjectAnd::propagateUpward( BuildObjectPredicate *pred ) const 
{
  BuildObject **x = pred->scope;
  int lb = (!x[0]->contain(0) && !x[1]->contain(0));
  int ub = (!x[0]->equal(0) && !x[1]->equal(0));  
  int consistent = ( pred->setMin( lb ) && pred->setMax( ub ) );
  if( consistent ) {
    if ( pred->isTop() ) 
      consistent = pred->setDomain(1);
    else if(lb < ub) {
      if(!x[0]->contain(0))
	pred->setEqual(x[1]);
      if(!x[1]->contain(0))
	pred->setEqual(x[0]);
    }
  }
  return consistent;
}

int BuildObjectAnd::propagateDownward( BuildObjectPredicate *pred ) const 
{
  bool consistent = true;
  BuildObject **x = pred->scope;
  if( pred->equal(0) ) {
    if( !x[0]->contain(0) )
      consistent = x[1]->setDomain(0);
    else if( !x[1]->contain(0) )
      consistent = x[0]->setDomain(0);
  } else if( !pred->contain(0) ) 
    consistent = (x[0]->remove(0) && x[1]->remove(0));

  return consistent;
}

void BuildObjectAnd::close( BuildObjectPredicate *pred ) 
{  
  //bool valid = !pred->contain(0);
  BuildObject **x = pred->scope;
  bool valid = (x[0]->isConstant() || x[1]->isConstant());
  //if( !valid && pred->equal(0) ) 
  //valid = ( x[0]->equal(0) || x[1]->equal(0) );
  if( valid ) {
    pred->unsetRel();
    x[0]->deReference();
    x[1]->deReference();
  } else {
    x[0]->unsetBList();
    x[0]->unsetIList();
    x[1]->unsetBList();
    x[1]->unsetIList();
    pred->unsetBList();
    pred->unsetIList();
  }
}

void BuildObjectAnd::print(std::ostream& o, const BuildObjectPredicate *pred) const
{
  o << "(";
  pred->scope[0]->print( o );
  o << " & ";
  pred->scope[1]->print( o );
  o << ")";
}

string BuildObjectAnd::toString(const BuildObjectPredicate* pred) const 
{
  string x1 = pred->scope[0]->toString();
  string x2 = pred->scope[1]->toString();
  string s = ( "(" + x1 + " and " + x2 + ")" );
  return s;
}

string BuildObjectAnd::xmlPred(int& idx, int level, const BuildObjectPredicate *pred, int *sub) const 
{
  string x1 = pred->scope[0]->xmlPred(idx, level+1, sub);
  string x2 = pred->scope[1]->xmlPred(idx, level+1, sub);
  string s = ( "and(" + x1 + "," + x2 + ")" );
  return s;
}


/**********************************************
 * Or Constraint BuildObject
 **********************************************/
void BuildObjectOr::build( Solver *s, VariableInt **tmp, BuildObjectPredicate *pred ) 
{
  if( !pred->isConstant() ) {
    new PredicateOr( s, tmp );
  } else if( pred->equal(1) ) {
    new ConstraintOr( s, tmp );
  } else {
    cout << "s NOT SUPPORTED" << endl;
    exit(0);
  }
}

void BuildObjectOr::build( SatSolver *s, const int *idx, BuildObjectPredicate *pred )
{
  Vector<Literal> new_clause;

  if( !pred->isConstant()  ) {

    new_clause.push( idx[pred->scope[0]->getVarId()] );
    new_clause.push( idx[pred->scope[1]->getVarId()] );
    new_clause.push( -(idx[pred->scope[2]->getVarId()]) );
    s->addClause( s->base, new_clause, s->stats.base_avg_size );
    s->addOriginalClause( new_clause );
    new_clause.clear();

    new_clause.push( -(idx[pred->scope[0]->getVarId()]) );
    new_clause.push( idx[pred->scope[2]->getVarId()] );
    s->addClause( s->base, new_clause, s->stats.base_avg_size );
    s->addOriginalClause( new_clause );
    new_clause.clear();

    new_clause.push( -(idx[pred->scope[1]->getVarId()]) );
    new_clause.push( idx[pred->scope[2]->getVarId()] );
    s->addClause( s->base, new_clause, s->stats.base_avg_size );
    s->addOriginalClause( new_clause );

  } else if( pred->equal(1) ) {

    new_clause.push( idx[pred->scope[0]->getVarId()] );
    new_clause.push( idx[pred->scope[1]->getVarId()] );
    s->addClause( s->base, new_clause, s->stats.base_avg_size );
    s->addOriginalClause( new_clause );

  }
} 

int BuildObjectOr::propagateUpward( BuildObjectPredicate *pred ) const 
{
  BuildObject **x = pred->scope;
  int lb = (!x[0]->contain(0) || !x[1]->contain(0));
  int ub = (!x[0]->equal(0) || !x[1]->equal(0));
  int consistent = ( pred->setMin( lb ) && pred->setMax( ub ) );
  if( consistent && pred->isTop() ) 
    consistent = pred->setDomain(1);
  return consistent;
}

int BuildObjectOr::propagateDownward( BuildObjectPredicate *pred ) const 
{
  bool consistent = true;
  BuildObject **x = pred->scope;
  if( !pred->contain(0) ) {
    if( x[0]->equal(0) )
      consistent = x[1]->remove(0);
    else if( x[1]->equal(0) )
      consistent = x[0]->remove(0);
  } else if( pred->equal(0) ) 
    consistent = (x[0]->setDomain(0) && x[1]->setDomain(0));
  return consistent;
}

void BuildObjectOr::close( BuildObjectPredicate *pred ) 
{  
  bool valid = pred->equal(0);
  BuildObject **x = pred->scope;
  if( !valid && !pred->contain(0) ) 
    valid = ( !x[0]->contain(0) || !x[1]->contain(0) );
  else if(pred->isTop() && pred->isReferenced() == 1) {
    valid = true;
  }
  if( valid ) {
    pred->unsetRel();
    pred->scope[0]->deReference();
    pred->scope[1]->deReference();
  } else {
//     if(pred->isTop()) {
//       cout << "deref" << endl;
//       pred->deReference();
//     }
    x[0]->unsetBList();
    x[0]->unsetIList();
    x[1]->unsetBList();
    x[1]->unsetIList();
    pred->unsetBList();
    pred->unsetIList();
  }
}

void BuildObjectOr::print(std::ostream& o, const BuildObjectPredicate *pred) const
{
  o << "(";
  pred->scope[0]->print( o );
  o << " || ";
  pred->scope[1]->print( o );
  o << ")";
}

string BuildObjectOr::toString(const BuildObjectPredicate* pred) const 
{
  string x1 = pred->scope[0]->toString();
  string x2 = pred->scope[1]->toString();
  string s = ( "(" + x1 + " or " + x2 + ")" );
  return s;
}

string BuildObjectOr::xmlPred(int& idx, int level, const BuildObjectPredicate *pred, int *sub) const 
{
  string x1 = pred->scope[0]->xmlPred(idx, level+1, sub);
  string x2 = pred->scope[1]->xmlPred(idx, level+1, sub);
  string s = ( "or(" + x1 + "," + x2 + ")" );
  return s;
}


/**********************************************
 *  Equal Predicate BuildObject
 **********************************************/ 
void BuildObjectEqual::build( Solver *s, VariableInt **tmp, BuildObjectPredicate *pred ) 
{
  if( !pred->isConstant() ) {
    if( pred->arity == 2 )
      new PredicateEqual( s, tmp, pred->params[0] );
    else
      new PredicateEqualConstant( s, tmp, pred->params[1], pred->params[0] );
  } else {
    if( pred->equal(0) == pred->params[0] )
      new ConstraintNotEqual( s, tmp );
  }
}

void BuildObjectEqual::build( SatSolver *s, const int *idx, BuildObjectPredicate *pred )
{
  Vector<Literal> new_clause;
  int i, spin = pred->params[0];

  if( !pred->isConstant() ) {
    if( pred->arity == 2 ) {
      for(i=0; i<2; ++i) {
	new_clause.push( idx[pred->scope[i]->getVarId()] );
	new_clause.push( -(idx[pred->scope[1-i]->getVarId()]) );
	new_clause.push( (spin ? -1 : 1) * idx[pred->scope[2]->getVarId()] );
	s->addClause( s->base, new_clause, s->stats.base_avg_size );
	s->addOriginalClause( new_clause );
	new_clause.clear();
      }

      for(i=0; i<2; ++i) {
	new_clause.push( (i ? -1 : 1) * idx[pred->scope[0]->getVarId()] );
	new_clause.push( (i ? -1 : 1) * idx[pred->scope[1]->getVarId()] );
	new_clause.push( (spin ? 1 : -1) * idx[pred->scope[2]->getVarId()] );
	s->addClause( s->base, new_clause, s->stats.base_avg_size );
	s->addOriginalClause( new_clause );
	new_clause.clear();
      }
    } else {      
      cerr << "NOT IMPLEMENTED!" << endl;
      exit(1);
    }
  } else {
    if( pred->equal(0) == pred->params[0] ) {
      for(i=0; i<2; ++i) {
	new_clause.push( (i ? -1 : 1) * idx[pred->scope[0]->getVarId()] );
	new_clause.push( (i ? -1 : 1) * idx[pred->scope[1]->getVarId()] );
	s->addClause( s->base, new_clause, s->stats.base_avg_size );
	s->addOriginalClause( new_clause );
	new_clause.clear();
      }
    }
  }
} 

int BuildObjectEqual::propagateUpward( BuildObjectPredicate *pred ) const 
{

  //  cout << "up" << endl;

  BuildObject **x = pred->scope;
  int n=pred->arity, spin=pred->params[0], b[2], k=pred->params[1];
  int consistent = true;
  //////
  if(pred->isConstant() && n==2 && (spin != pred->value())) {
    if(x[1]->isConstant()) {
      pred->params[1] = x[1]->value();
      pred->arity = --n;
      x[1]->deReference();
      x[1] = x[2];
    } else if(x[0]->isConstant()) {
      pred->params[1] = x[0]->value();
      pred->arity = --n;
      x[0]->deReference();
      x[0] = x[1];
      x[1] = x[2];
    }
    consistent = propagateDownward(pred);
  }
  ///////

//   pred->print(cout);
//   cout << " " << consistent << endl;

  if( consistent ) {
    if( n == 2 ) {
      b[spin] = (spin ^ !x[0]->intersect(x[1]));
      b[1-spin] = (spin ^ (x[0]->size() != 1 || x[1]->size() != 1 || !x[0]->equal(x[1]->value())));
    } else {
      b[spin] = (spin ^ !x[0]->contain(k));
      b[1-spin] = (spin ^ (x[0]->size() != 1 || !x[0]->equal(k)));
    }

    consistent &= ( pred->setMin(b[0]) && pred->setMax(b[1]) );
    
    //    cout << "---- " << (pred->isTop()) << endl;

    if( consistent && pred->isTop() ) {
      consistent = pred->setDomain(1);

//       pred->print(cout);
//       cout << " " << consistent << endl;


    }
  }
  
//   pred->print(std::cout);
//   std::cout << " " << consistent << std::endl << std::endl;

  return consistent;
}

int BuildObjectEqual::propagateDownward( BuildObjectPredicate *pred ) const 
{
  //cout << "down" << endl;
  //exit(0);

  int consistent = true;

  if( pred->isConstant() ) {
    int n = pred->arity, spin = pred->params[0], k = pred->params[1];
    BuildObject **x = pred->scope;
    if( spin == pred->equal(0) ) {
      if( n == 2 ) 
	consistent = ((!x[0]->isConstant() || x[1]->remove(x[0]->value()))
		      &&
		      (!x[1]->isConstant() || x[0]->remove(x[1]->value())));
      else
	consistent = ((!x[0]->isConstant() || x[0]->value() != k)
		      &&
		      (x[0]->remove(k)));
    } else {
      if( n == 2 ) 
	consistent = x[0]->setEqual(x[1]);
      else
	consistent = x[0]->setDomain(k);
    }
  }

//   pred->print(cout);
//   cout << endl;

  return consistent;
}

void BuildObjectEqual::close( BuildObjectPredicate *pred ) 
{  
  BuildObject **x = pred->scope;
  int spin = pred->params[0], valid, n=pred->arity;
  if( spin )
    valid  = (!pred->contain(0) ||
	      (pred->equal(0) && (x[0]->isConstant() || x[1]->isConstant())));
  else
    valid = (pred->equal(0) ||
	     (!pred->contain(0) && (x[0]->isConstant() || x[1]->isConstant())));
  
  if( valid ) {
    pred->unsetRel();
    for(int i=0; i<n; ++i)
      x[i]->deReference();
    //x[1]->deReference();
  } else {
    for(int i=0; i<n; ++i) {
      x[i]->unsetBList();
      x[i]->unsetIList();
    }
    //x[1]->unsetBList();
    //x[1]->unsetIList();
    pred->unsetBList();
    pred->unsetIList();

    if( pred->isConstant() && spin == pred->equal(0) ) 
      pred->setNeq();

    if( !pred->isConstant() || spin == pred->equal(0) ) {
      for(int i=0; i<n; ++i)
	x[i]->unsetRange();
      //x[1]->unsetRange();
    }
    if(pred->isTop() && pred->isConstant())
      pred->deReference();
    else if( n == 2 ) {
      if( x[0]->isConstant() ) {
	x[0]->deReference();
	pred->params[1] = x[0]->value();
	x[0] = x[1];
	x[1] = pred;
	--pred->arity;
      } else if( x[1]->isConstant() ) {
	x[1]->deReference();
	pred->params[1] = x[1]->value();
	x[1] = pred;
	--pred->arity;
      }
    }
  }
}

void BuildObjectEqual::print(std::ostream& o, const BuildObjectPredicate *pred) const
{
  o << "(";
  pred->scope[0]->print( o );
  if( pred->params[0] )
    o << " == " ;
  else
    o << " =/= " ;
  if(pred->arity == 2)
    pred->scope[1]->print( o );
  else 
    cout << pred->params[1];
  o << ")";
}

string BuildObjectEqual::toString(const BuildObjectPredicate* pred) const 
{
  string x1 = pred->scope[0]->toString();
  string x2 = pred->scope[1]->toString();
  string s = ( "(" + x1 + (pred->params[0] ? " == " : " =/= ") + x2 + ")" );
  return s;
}

string BuildObjectEqual::xmlPred(int& idx, int level, const BuildObjectPredicate *pred, int *sub) const 
{
  string x1 = pred->scope[0]->xmlPred(idx, level+1, sub);
  string x2 = pred->scope[1]->xmlPred(idx, level+1, sub);
  string s = ( "eq(" + x1 + "," + x2 + ")" );
  return s;
}


/**********************************************
 * IfThenElse Constraint BuildObject
 **********************************************/
void BuildObjectIfThenElse::build( Solver *s, VariableInt **tmp, BuildObjectPredicate *pred ) 
{
  new PredicateIfThenElse( s, tmp );
}

int BuildObjectIfThenElse::propagateUpward( BuildObjectPredicate *pred ) const 
{
  BuildObject **x = pred->scope;
  int lb, ub, consistent = 1;
  if(x[0]->equal(0)) {
    consistent = pred->setEqual( x[2] );
  } else if(!x[0]->contain(0)) {
    consistent = pred->setEqual( x[1] );
  } else {
    lb = std::min( x[1]->min(), x[2]->min() );
    ub = std::max( x[1]->max(), x[2]->max() );
    consistent = ( pred->setMin( lb ) && pred->setMax( ub ) );
  }
  return consistent;
}

int BuildObjectIfThenElse::propagateDownward( BuildObjectPredicate *pred ) const 
{
  bool consistent = true;
  BuildObject **x = pred->scope;

  if(x[0]->equal(0)) {
    consistent = pred->setEqual( x[2] );
  } else if(!x[0]->contain(0)) {
    consistent = pred->setEqual( x[2] );
  } 
  else {
    if( !pred->intersect(x[1]) ) {
      consistent = ( x[0]->setDomain(0) && pred->setEqual(x[2]) );
    } else if( !pred->intersect(x[2]) ) {
      consistent = ( x[0]->remove(0) && pred->setEqual(x[1]) );
    }
  }
  return consistent;
}

void BuildObjectIfThenElse::close( BuildObjectPredicate *pred ) 
{  
  BuildObject **x = pred->scope;
  bool valid = (x[0]->equal(0) || !x[0]->contain(0));
  if( valid ) {
    pred->unsetRel();
    x[0]->deReference();
    x[1]->deReference();
    x[2]->deReference();
  } else {
    pred->model->unsetSat();
    x[0]->unsetBList();
    x[0]->unsetIList();
    x[1]->unsetBList();
    x[1]->unsetIList();
    x[2]->unsetBList();
    x[2]->unsetIList();
    pred->unsetBList();
    pred->unsetIList();
  }
}

void BuildObjectIfThenElse::print(std::ostream& o, const BuildObjectPredicate *pred) const
{
  o << "(";
  pred->scope[0]->print( o );
  o << " ? ";
  pred->scope[1]->print( o );
  o << " : ";
  pred->scope[2]->print( o );
  o << ")";
}

string BuildObjectIfThenElse::toString(const BuildObjectPredicate* pred) const 
{
  string x1 = pred->scope[0]->toString();
  string x2 = pred->scope[1]->toString();
  string x3 = pred->scope[2]->toString();
  string s = ( "if " + x1 + " then " + x2 + " else " + x3 );
  return s;
}

string BuildObjectIfThenElse::xmlPred(int& idx, int level, const BuildObjectPredicate *pred, int *sub) const 
{
  string x1 = pred->scope[0]->xmlPred(idx, level+1, sub);
  string x2 = pred->scope[1]->xmlPred(idx, level+1, sub);
  string x3 = pred->scope[2]->xmlPred(idx, level+1, sub);
  string s = ( "if(" + x1 + "," + x2 + "," + x3 + ")" );
  return s;
}


/**********************************************
 *  Disjunctive Predicate BuildObject
 **********************************************/ 
void BuildObjectDisjunctive::build( Solver *s, VariableInt **tmp, BuildObjectPredicate *pred ) 
{
  if(pred->params[2]) {
    if( !pred->model->clauseBase )
      pred->model->clauseBase = new ConstraintClauseBase(s);
    new PredicateDisjunct( s, tmp, pred->params, pred->model->clauseBase );
  } else {
    if(pred->equal(0)) {
      new ConstraintLess( s, tmp, pred->params[0] );
    } else if(pred->equal(1)) {
      pred->params[0] = pred->params[1];
      VariableInt *aux = tmp[0];
      tmp[0] = tmp[1];
      tmp[1] = aux;
      new ConstraintLess( s, tmp, pred->params[0] );
    } else {
#ifdef _SAFE
      new PredicateSafeDisjunct( s, tmp, pred->params );
#else
      if((tmp[0]->getType() == VariableInt::RANGE ||
	  tmp[0]->getType() == VariableInt::BIT)
	 && 
	 (tmp[1]->getType() == VariableInt::RANGE || 
	  tmp[1]->getType() == VariableInt::BIT)
	 ) {
	//std::cout << "optimised" << std::endl;
	new PredicateDisjunctive( s, tmp, pred->params );
      } else {
	//std::cout << "safe" << std::endl;
	new PredicateSafeDisjunct( s, tmp, pred->params );
      }
#endif
    }
  }
}

int BuildObjectDisjunctive::propagateUpward( BuildObjectPredicate *pred ) const 
{
//   BuildObject **x = pred->scope;
//   int *d = pred->params;
//   int lb = (x[0]->min() + d[0] > x[1]->max());
//   int ub = (x[1]->min() + d[1] <= x[0]->max());
//   return ( lb <= ub && pred->setMin( lb ) && pred->setMax( ub ) && propagateDownward(pred) );

  return propagateDownward(pred);
}

int BuildObjectDisjunctive::propagateDownward( BuildObjectPredicate *pred ) const 
{
  bool consistent = true;
  BuildObject **x = pred->scope;
  int *d = pred->params;
  int lb = (x[0]->min() + d[0] > x[1]->max());
  int ub = (x[1]->min() + d[1] <= x[0]->max());
  consistent = ( lb <= ub && pred->setMin( lb ) && pred->setMax( ub ) );

  
  //BuildObject **x = pred->scope;
  //int *d = pred->params;

  if(consistent) {
    int lb, ub;
    if(pred->equal(0)) {
      ub = x[1]->max()-d[0];
      lb = x[0]->min()+d[0];
      consistent = ( x[0]->min() <= ub && x[1]->max() >= lb );
      if( consistent ) {
	x[0]->setMax( ub );
	x[1]->setMin( lb );
      }
    } else if(pred->equal(1)) {
      ub = x[0]->max()-d[1];
      lb = x[1]->min()+d[1];
      consistent = ( x[1]->min() <= ub && x[0]->max() >= lb );
      if( consistent ) {
	x[1]->setMax( ub );
	x[0]->setMin( lb );
      }
    } 
  }

  return consistent;
}

void BuildObjectDisjunctive::close( BuildObjectPredicate *pred ) 
{  
  BuildObject **x = pred->scope;
  bool valid = (pred->isConstant() && (x[0]->isConstant() || x[1]->isConstant()));
  if( valid ) {
    pred->unsetRel();
    x[0]->deReference();
    x[1]->deReference();
  } else {
    pred->model->unsetSat();
    x[0]->unsetBList();
    x[0]->unsetIList();
    x[1]->unsetBList();
    x[1]->unsetIList();
    pred->unsetBList();
    pred->unsetIList();
  }
}

void BuildObjectDisjunctive::print(std::ostream& o, const BuildObjectPredicate *pred) const
{
  o << "(";
  pred->scope[0]->print( o );
  o << " + " << pred->params[0] << " <= ";
  pred->scope[1]->print( o );
  o << " || ";
  pred->scope[1]->print( o );
  o << " + " << pred->params[1] << " <= ";
  pred->scope[0]->print( o );
  o << ") ";
  o << (pred->isReferenced()) << " ";
  pred->BuildObject::print(cout);
}

string BuildObjectDisjunctive::toString(const BuildObjectPredicate* pred) const 
{
  string x1 = pred->scope[0]->toString();
  string x2 = pred->scope[1]->toString();
  string s = ( x1 + " + " + int2string(pred->params[0]) + " <= " + x2 + " or " + x2 + " + " + int2string(pred->params[1]) + " <= " + x1 );
  return s;
}

string BuildObjectDisjunctive::xmlPred(int& idx, int level, const BuildObjectPredicate *pred, int *sub) const 
{
  string x1 = pred->scope[0]->xmlPred(idx, level+1, sub);
  string x2 = pred->scope[1]->xmlPred(idx, level+1, sub);
  string s = ( "or(le(add(" + x1 + "," + int2string(pred->params[0]) + ")," + x2 + ")," +
	       "le(add(" + x2 + "," + int2string(pred->params[1]) + ")," + x1 + "))" );
  return s;
}


/**********************************************
 *  Overlap Predicate BuildObject
 **********************************************/ 
void BuildObjectOverlap::build( Solver *s, VariableInt **tmp, BuildObjectPredicate *pred ) 
{
  if(pred->equal(0)) // x0 << x1
    new ConstraintLess( s, tmp, pred->params[0] );
  else if(pred->equal(1)) { // x1 << x0
    pred->params[0] = pred->params[1];
    VariableInt *aux = tmp[0];
    tmp[0] = tmp[1];
    tmp[1] = aux;
    new ConstraintLess( s, tmp, pred->params[0] );
  } 
  else {
    if(pred->size() > 2) {
      new PredicateOverlap(s, tmp, pred->params );
    } else if(!pred->contain(2)) { // regular disjunctive
      new PredicateDisjunctive( s, tmp, pred->params );
    } else {
      if(!pred->contain(0)) { // not( x0 << x1 ) -> x1 + (1 - d0) <= x0
	int d = (1-pred->params[0]);
	VariableInt *aux = tmp[0];
	tmp[0] = tmp[1];
	tmp[1] = aux;
	new ConstraintLess( s, tmp, d );
      }
      if(!pred->contain(1)) { // not( x1 << x0 ) -> x0 + (1 - d1) <= x1
	int d = (1-pred->params[1]);
	new ConstraintLess( s, tmp, d );
      }
    }
  }
}

int BuildObjectOverlap::propagateUpward( BuildObjectPredicate *pred ) const 
{
  return propagateDownward(pred);
}

int BuildObjectOverlap::propagateDownward( BuildObjectPredicate *pred ) const 
{
  bool consistent = (pred->setMin(0) && pred->setMax(2));

  if(consistent) {
    BuildObject **x = pred->scope;
    int *d = pred->params;
    int x1_is_before_x0 = (x[1]->max() + d[1] <= x[0]->min());
    int x1_canbe_before_x0 = (x[1]->min() + d[1] <= x[0]->max());
    int x0_is_before_x1 = (x[0]->max() + d[0] <= x[1]->min());
    int x0_canbe_before_x1 = (x[0]->min() + d[0] <= x[1]->max());
    
    if(consistent && x0_is_before_x1) consistent &= pred->setDomain(0);
    if(consistent && x1_is_before_x0) consistent &= pred->setDomain(1);
    if(consistent && !x0_canbe_before_x1) consistent &= pred->remove(0);
    if(consistent && !x1_canbe_before_x0) consistent &= pred->remove(1);
  }
//   if(consistent) {
//     // to do
    
//   }

  return consistent;
}

void BuildObjectOverlap::close( BuildObjectPredicate *pred ) 
{  
  BuildObject **x = pred->scope;
  bool valid = (pred->isConstant() && (x[0]->isConstant() || x[1]->isConstant()));
  if( valid ) {
    pred->unsetRel();
    x[0]->deReference();
    x[1]->deReference();
  } else {
    pred->model->unsetSat();
    x[0]->unsetBList();
    x[0]->unsetIList();
    x[1]->unsetBList();
    x[1]->unsetIList();
    pred->unsetBList();
    pred->unsetIList();
  }
}

void BuildObjectOverlap::print(std::ostream& o, const BuildObjectPredicate *pred) const
{
  o << "(";
  pred->scope[0]->print( o );
  o << " + " << pred->params[0] << " <= ";
  pred->scope[1]->print( o );
  o << " || ";
  pred->scope[1]->print( o );
  o << " + " << pred->params[1] << " <= ";
  pred->scope[0]->print( o );
  o << ") ";
  o << (pred->isReferenced()) << " ";
  pred->BuildObject::print(cout);
}

string BuildObjectOverlap::toString(const BuildObjectPredicate* pred) const 
{
  string x1 = pred->scope[0]->toString();
  string x2 = pred->scope[1]->toString();
  string s = ( x1 + " + " + int2string(pred->params[0]) + " <= " + x2 + " or " + x2 + " + " + int2string(pred->params[1]) + " <= " + x1 );
  return s;
}

string BuildObjectOverlap::xmlPred(int& idx, int level, const BuildObjectPredicate *pred, int *sub) const 
{
  string x1 = pred->scope[0]->xmlPred(idx, level+1, sub);
  string x2 = pred->scope[1]->xmlPred(idx, level+1, sub);
  string s = ( "or(le(add(" + x1 + "," + int2string(pred->params[0]) + ")," + x2 + ")," +
	       "le(add(" + x2 + "," + int2string(pred->params[1]) + ")," + x1 + "))" );
  return s;
}


/**********************************************
 *  Precedence Predicate BuildObject
 **********************************************/ 
void BuildObjectPrecedence::build( Solver *s, VariableInt **tmp, BuildObjectPredicate *pred ) 
{
  BuildObject **x = pred->scope;
  if( pred->isConstant() ) { 
    if( !x[0]->isConstant() && !x[1]->isConstant() ) {
      if( pred->equal(0) ) {
	VariableInt *aux = tmp[0];
	tmp[0] = tmp[1];
	tmp[1] = aux;
	pred->params[0] = (1-pred->params[0]);
      }    
      new ConstraintLess( s, tmp, pred->params[0] );
    }
  } else {
    int b = pred->params[0];
    if( x[0]->isConstant() ) {
      b += x[0]->value();
      new PredicateLowerBound( s, &(tmp[1]), b );      
    } else if( x[1]->isConstant() ) {
      b = x[1]->value() - b;
      tmp[1] = tmp[0];
      new PredicateUpperBound( s, &(tmp[1]), b );
    } else new PredicateLess( s, tmp, pred->params[0] );
  }
}

int BuildObjectPrecedence::propagateUpward( BuildObjectPredicate *pred ) const 
{
//   BuildObject **x = pred->scope;
//   int d = pred->params[0];
//   int lb = (x[0]->max() + d <= x[1]->min());
//   int ub = (x[0]->min() + d <= x[1]->max());
//   int consistent = ( pred->setMin( lb ) && pred->setMax( ub ) );
//   if( consistent && pred->isTop() ) 
//     consistent = pred->setDomain(1);
//  return (consistent && propagateDownward(pred));
  return propagateDownward(pred);
}

int BuildObjectPrecedence::propagateDownward( BuildObjectPredicate *pred ) const 
{

//    pred->print(cout);
//    cout << endl;

  BuildObject **x = pred->scope;
  int d = pred->params[0];
  int lb = (x[0]->max() + d <= x[1]->min());
  int ub = (x[0]->min() + d <= x[1]->max());

//    cout << (x[0]->min()) << ".." << (x[0]->max()) << " | "
//         << (x[1]->min()) << ".." << (x[1]->max()) << endl;

  int consistent = ( pred->setMin( lb ) && pred->setMax( ub ) );


  if( consistent && pred->isTop() ) {
    if( pred->size() > 1 ) 
      consistent = pred->setDomain(1);
    else
      consistent = pred->equal(1);
  }
  //  cout << lb << " " << ub << " " << consistent << endl;

  //bool consistent = true;
  //BuildObject **x = pred->scope;
  //int d = pred->params[0];

  if(consistent) {

//     pred->print(cout);
//     cout << endl;
    
    
    if(pred->equal(0)) {
      
      //      cout << "propagate opposite" << endl;
      
      consistent = ( x[0]->setMin( x[1]->min()-d+1 ) && x[1]->setMax( x[0]->max()+d-1 ) );
    } else if(pred->equal(1)) {
      
      //      cout << "propagate normal" << endl;
      
      consistent = ( x[0]->setMax( x[1]->max()-d ) && x[1]->setMin( x[0]->min()+d ) );
    } 
    

//     pred->print(cout);
//     cout << ": " << consistent << endl;
  }

  return consistent;
}

void BuildObjectPrecedence::close( BuildObjectPredicate *pred ) 
{  
//   cout << "close ";
//   cout << endl;
//   pred->print( cout );
//   cout << endl;

  BuildObject **x = pred->scope;
  bool valid = (pred->isConstant() && (x[0]->isConstant() || x[1]->isConstant()));
  if( valid ) {
    pred->unsetRel();
    x[0]->deReference();
    x[1]->deReference();
  } else {
    pred->model->unsetSat();
    if(pred->isTop() && pred->isConstant()) {
      pred->deReference();
      if( (pred->equal(1) == pred->params[0]) > 0 )
	pred->setNeq();
    }
    for(int i=0; i<2; ++i) {
      if(x[i]->isConstant())
	x[i]->deReference();
      x[i]->unsetBList();
      x[i]->unsetIList();
    }
    pred->unsetBList();
    pred->unsetIList();
  }
}

void BuildObjectPrecedence::print(std::ostream& o, const BuildObjectPredicate *pred) const
{
  o << "(";
  pred->scope[0]->print( o );
  if( pred->params[0] == 1 )
    o << " < ";
  else {
    if( pred->params[0] )
      o << " + " << pred->params[0] ;
    o << " <= ";
  }
  pred->scope[1]->print( o );
  o << ")";
}    

string BuildObjectPrecedence::toString(const BuildObjectPredicate* pred) const 
{
  string x1 = pred->scope[0]->toString();
  string x2 = pred->scope[1]->toString();
  string s = ( x1 + " <= " + x2 );
  return s;
}

string BuildObjectPrecedence::xmlPred(int& idx, int level, const BuildObjectPredicate *pred, int *sub) const 
{
  string x1 = pred->scope[0]->xmlPred(idx, level+1, sub);
  string x2 = pred->scope[1]->xmlPred(idx, level+1, sub);
  string s;
  if(pred->params[0] == 0 )
    s = ( "le(" + x1 + "," + x2 + ")" );
  else if(pred->params[0] == 1)
    s = ( "lt(" + x1 + "," + x2 + ")" );
  else
    s = ( "le(add(" + x1 + "," + int2string(pred->params[0]) + ")," + x2 + ")" );
  return s;
}


/**********************************************
 *  Mul Predicate BuildObject
 **********************************************/
void BuildObjectMul::self_build( Solver *s, VariableInt **tmp, BuildObjectPredicate *pred )
{
  BuildObject **x = pred->scope;
  
  assert( !x[0]->isConstant() || !x[1]->isConstant() );

  if( x[0]->isConstant() )
    {
      VariableVirtual *vv;
      s->MISC += 1000;
      vv = new VariableVirtual( s );
      if(pred->veqptr_) {
	pred->veqptr_->varptr_ = vv;
	//pred->veqptr_->setVirtual();
      } else {
	pred->varptr_ = vv; //new VariableVirtual( s );
      }
// 	((VariableVirtual*)(pred->varptr_))->reference = tmp[1];
// 	((VariableVirtual*)(pred->varptr_))->conversion = new MapFactor( x[0]->min(), ((VariableVirtual*)(pred->varptr_)) );
      vv->reference = tmp[1];
      vv->conversion = new MapFactor( x[0]->min(), vv ); //((VariableVirtual*)(pred->varptr_)) );
      pred->setVirtual();
      //pred->veqptr_ = x[1];
    }
  else if( x[1]->isConstant() )
    {
      VariableVirtual *vv;
      s->MISC += 1000;
      vv = new VariableVirtual( s );
      if(pred->veqptr_) {
	pred->veqptr_->varptr_ = vv;
	//pred->veqptr_->setVirtual();
      } else {
	pred->varptr_ = vv; //new VariableVirtual( s );
      }
// 	((VariableVirtual*)(pred->varptr_))->reference = tmp[1];
// 	((VariableVirtual*)(pred->varptr_))->conversion = new MapFactor( x[0]->min(), ((VariableVirtual*)(pred->varptr_)) );
      vv->reference = tmp[0];
      vv->conversion = new MapFactor( x[1]->min(), vv ); //((VariableVirtual*)(pred->varptr_)) );
      pred->setVirtual();
      //pred->veqptr_ = x[1];

//       s->MISC += 1000;
//       pred->varptr_ = new VariableVirtual( s );
//       pred->setVirtual();
//       //pred->veqptr_ = x[0];
//       ((VariableVirtual*)(pred->varptr_))->reference = tmp[0];
//       ((VariableVirtual*)(pred->varptr_))->conversion = new MapFactor( x[1]->min(), ((VariableVirtual*)(pred->varptr_)) );
    }
  else
    {
      pred->BuildObject::build(s);
    }
}

void BuildObjectMul::build( Solver *s, VariableInt **tmp, BuildObjectPredicate *pred ) 
{  
  BuildObject **x = pred->scope;
  if( x[0]->isConstant() )
    {
//       int f = x[0]->value();
//       new PredicateFactor( s, &tmp[1], f );
    }
  else if( x[1]->isConstant() )
    {
//       int f = x[1]->value();
//       tmp[1] = tmp[0];
//       new PredicateFactor( s, &tmp[1], f );
    }
  else new PredicateMul( s, tmp );
}

int BuildObjectMul::propagateUpward( BuildObjectPredicate *pred ) const 
{
  BuildObject **x = pred->scope;
  int z[4], lb = NOVAL, ub = -1*NOVAL, i = 4;
  z[0] = x[0]->min() * x[1]->min();
  z[1] = x[0]->min() * x[1]->max();
  z[2] = x[0]->max() * x[1]->min();
  z[3] = x[0]->max() * x[1]->max();  
  while( i-- ) {
    if( z[i] > ub ) ub = z[i];
    if( z[i] < lb ) lb = z[i];
  }
  return ( pred->setMin( lb ) && pred->setMax( ub ) && (x[0]->contain(0) || x[1]->contain(0) || pred->remove(0)) );
}

int BuildObjectMul::propagateDownward( BuildObjectPredicate *pred ) const 
{
  bool consistent = true;
  BuildObject **x = pred->scope;

  if( pred->equal(0) ) {
    if( !x[0]->contain(0) )
      consistent = x[1]->setDomain(0);
    if( !x[1]->contain(0) )
      consistent = x[0]->setDomain(0);
  } else {
    int k, lb, ub, i;

    if( !pred->contain(0) ) 
      consistent = ( x[0]->remove(0) && x[1]->remove(0) );
    
    if( consistent && pred->isConstant() ) {
      k = pred->value();
      
      for(i=0; i<2; ++i) 
	if( !x[i]->contain(0) ) {
	  
	  ub = x[i]->minPosAbs();
	  if( ub != NOVAL )
	    ub = k/ub;
	  else
	    ub = k/x[i]->min();
	  
	  lb = x[i]->minNegAbs();
	  if( lb != NOVAL )
	    lb = k/lb;
	  else
	    lb = k/x[i]->max();
	  
	  if( k>=0 )
	    consistent = (x[1-i]->setMin(lb) && x[1-i]->setMax(ub));
	  else
	    consistent = (x[1-i]->setMin(ub) && x[1-i]->setMax(lb));
	}
    } else {

      int factor, f, rest;
      
      for(i=0; i<2; ++i)
	if( x[i]->isConstant() ) {
	  
	  factor = x[i]->value(); 
	  if( factor ) {
	    f = (factor > 0 ? factor : -factor);
	    
	    lb = pred->min();
	    rest = (lb % factor);   
	    
	    if( rest ) {
	      if( lb > 0 ) lb += f;
	      lb -= rest;
	      
	      consistent = pred->setMin( lb );
	      lb = pred->min();
	    }
	    
	    if( consistent ) {
	      ub = pred->max();
	      rest = (ub % factor);
	      
	      if( rest ) {
		if( ub < 0 ) ub -= f;
		ub -= rest;
		
		consistent = pred->setMax( ub );
		ub = pred->max();
	      }
	      
	      if ( factor > 0 ) {
		consistent &= ( x[1-i]->setMin( lb / factor ) &&
				x[1-i]->setMax( ub / factor ) );	
	      } else {
		consistent &= ( x[1-i]->setMin( ub / factor ) &&
				x[1-i]->setMax( lb / factor ) );
	      }            
	    }
	  }
	}
    }
  }

  return consistent;
}

void BuildObjectMul::close( BuildObjectPredicate *pred ) 
{  
  BuildObject **x = pred->scope;
  bool valid = ( (pred->isConstant() && x[0]->isConstant() && x[1]->isConstant())
		 || (pred->equal(0) && (x[0]->equal(0) || x[1]->equal(0)))
		 );
  if( valid ) {
    pred->unsetRel();
    x[0]->deReference();
    x[1]->deReference();
  } else {
    pred->model->unsetSat();
    for(int i=0; i<2; ++i) {
      if( x[i]->isConstant() )
	x[i]->deReference();

      x[i]->unsetBList();
      x[i]->unsetIList();
    }
    pred->unsetBList();
    pred->unsetIList();
  }
}

void BuildObjectMul::print(std::ostream& o, const BuildObjectPredicate *pred) const
{
  o << "(";
  pred->scope[0]->print( o );
  o << " * ";
  pred->scope[1]->print( o );
  o << ")";
}

string BuildObjectMul::toString(const BuildObjectPredicate* pred) const 
{
  string x1 = pred->scope[0]->toString();
  string x2 = pred->scope[1]->toString();
  string s = ( "(" + x1 + " * " + x2 + ")" );
  return s;
}

string BuildObjectMul::xmlPred(int& idx, int level, const BuildObjectPredicate *pred, int *sub) const 
{
  string x1 = pred->scope[0]->xmlPred(idx, level+1, sub);
  string x2 = pred->scope[1]->xmlPred(idx, level+1, sub);
  string s = ( "mul(" + x1 + "," + x2 + ")" );
  return s;
}


/**********************************************
 *  Div Predicate BuildObject
 **********************************************/
void BuildObjectDiv::self_build( Solver *s, VariableInt **tmp, BuildObjectPredicate *pred )
{
  BuildObject **x = pred->scope;
  if( x[1]->isConstant() )
    {
      VariableVirtual *vv;
      s->MISC += 1000000;
      vv = new VariableVirtual( s );
      if(pred->veqptr_) {
	pred->veqptr_->varptr_ = vv;
	//pred->veqptr_->setVirtual();
      } else {
	pred->varptr_ = vv; //new VariableVirtual( s );
      }
// 	((VariableVirtual*)(pred->varptr_))->reference = tmp[1];
// 	((VariableVirtual*)(pred->varptr_))->conversion = new MapFactor( x[0]->min(), ((VariableVirtual*)(pred->varptr_)) );
      vv->reference = tmp[0];
      vv->conversion = new MapQuotient( x[1]->min(), vv ); //((VariableVirtual*)(pred->varptr_)) );
      pred->setVirtual();
      //pred->veqptr_ = x[1];

//       s->MISC += 1000000;
//       pred->varptr_ = new VariableVirtual( s );
//       pred->setVirtual();
//       //pred->veqptr_ = x[0];
//       ((VariableVirtual*)(pred->varptr_))->reference = tmp[0];
//       ((VariableVirtual*)(pred->varptr_))->conversion = new MapQuotient( tmp[1]->min(), ((VariableVirtual*)(pred->varptr_)) );
    }
  else
    {
      pred->BuildObject::build(s);
    }
}
void BuildObjectDiv::build( Solver *s, VariableInt **tmp, BuildObjectPredicate *pred ) 
{
  if( !pred->scope[1]->isConstant() )
    new PredicateDiv( s, tmp );
}

int BuildObjectDiv::propagateUpward( BuildObjectPredicate *pred ) const 
{
  BuildObject **x = pred->scope;
  int lb, ub, aux;
      
  int mp1 = x[1]->minPosAbs();
  int mn1 = x[1]->minNegAbs();
  int mp0 = x[0]->minPosAbs();
  int mn0 = x[0]->minNegAbs();
      
  //int m1 = x[1]->min();
  //int x1 = x[1]->max();
  int m0 = x[0]->min();
  int x0 = x[0]->max();
      
  if( !mp1 ) {
    if( x[1]->max() )
      mp1 = 1;
    else mp1 = -1;
  }
  if( !mn1 ) {
    if( x[1]->min() )
      mn1 = -1;
    else mn1 = 1;
  }
  if( !mp0 ) {
    if( x[0]->max() )
      mp0 = 1;
    else mp0 = -1;
  }
  if( !mn0 ) {
    if( x[0]->min() )
      mn0 = -1;
    else mn0 = 1;
  }
      
  lb = x0/mn1;
  aux = m0/mp1;
  if( lb > aux ) lb = aux;
      
  ub = x0/mp1;
  aux = m0/mn1;
  if( ub < aux ) ub = aux;

  return ( pred->setMin( lb )
	   &&
	   pred->setMax( ub ) );
}

int BuildObjectDiv::propagateDownward( BuildObjectPredicate *pred ) const 
{
  BuildObject **x = pred->scope;
  bool consistent = x[1]->remove(0);

  if( consistent ) { 
    if( pred->equal(0) ) {
      consistent = x[0]->setDomain(0);
    } else {
      int k, //z[4], 
	lb, ub;//, i;
      
      if( !pred->contain(0) ) {
	
	consistent = x[0]->remove(0);
	
	if( consistent ) {
	  lb = std::abs(x[0]->min());
	  ub = std::abs(x[0]->max());
	  if(ub<lb) ub=lb;
	  consistent = (x[1]->setMin(-ub) && x[1]->setMax(-ub));
	}

      }
      if( consistent && pred->isConstant() ) {
	k = pred->value();
	
	ub = k * x[1]->max();
	lb = k * x[1]->min();	    
	if( k>=0 )
	  consistent = (x[0]->setMin(lb) && x[0]->setMax(ub));
	else
	  consistent = (x[0]->setMin(ub) && x[0]->setMax(lb));

	ub = x[0]->max()/k;
	lb = x[0]->min()/k;
	if( k>=0 )
	  consistent = (x[1]->setMin(lb) && x[1]->setMax(ub));
	else
	  consistent = (x[1]->setMin(ub) && x[1]->setMax(lb));	
      }
    }
  }
  return consistent;
}

void BuildObjectDiv::close( BuildObjectPredicate *pred ) 
{  
  BuildObject **x = pred->scope;
  bool valid = ( (pred->isConstant() && x[0]->isConstant() && x[1]->isConstant())
		 || (pred->equal(0) && x[0]->equal(0))
		 );
  if( valid ) {
    pred->unsetRel();
    x[0]->deReference();
    x[1]->deReference();
  } else {
    pred->model->unsetSat();
    x[0]->unsetBList();
    x[0]->unsetIList();
    x[1]->unsetBList();
    x[1]->unsetIList();
    pred->unsetBList();
    pred->unsetIList();
  }
}

void BuildObjectDiv::print(std::ostream& o, const BuildObjectPredicate *pred) const
{
  o << "(";
  pred->scope[0]->print( o );
  o << " / ";
  pred->scope[1]->print( o );
  o << ")";
}

string BuildObjectDiv::toString(const BuildObjectPredicate* pred) const 
{
  string x1 = pred->scope[0]->toString();
  string x2 = pred->scope[1]->toString();
  string s = ( "(" + x1 + " / " + x2 + ")" );
  return s;
}

string BuildObjectDiv::xmlPred(int& idx, int level, const BuildObjectPredicate *pred, int *sub) const 
{
  string x1 = pred->scope[0]->xmlPred(idx, level+1, sub);
  string x2 = pred->scope[1]->xmlPred(idx, level+1, sub);
  string s = ( "div(" + x1 + "," + x2 + ")" );
  return s;
}


/**********************************************
 *  Mod Predicate BuildObject
 **********************************************/
void BuildObjectMod::build( Solver *s, VariableInt **tmp, BuildObjectPredicate *pred ) 
{
  new PredicateMod( s, tmp );
}

int BuildObjectMod::propagateUpward( BuildObjectPredicate *pred ) const 
{
  BuildObject **x = pred->scope;
  int lb = ( x[0]->min() < 0 ? -x[1]->max()+1 : 0 );
  int ub = ( x[0]->max() > 0 ?  x[1]->max()-1 : 0 );
  return ( pred->setMin( lb ) && pred->setMax( ub ) );
}

int BuildObjectMod::propagateDownward( BuildObjectPredicate *pred ) const 
{
  //BuildObject **x = pred->scope;
  bool consistent = true;
  return consistent;
}

void BuildObjectMod::close( BuildObjectPredicate *pred ) 
{  
  BuildObject **x = pred->scope;
  bool valid = ( (pred->isConstant() && x[0]->isConstant() && x[1]->isConstant()) );
  if( valid ) {
    pred->unsetRel();
    x[0]->deReference();
    x[1]->deReference();
  } else {
    pred->model->unsetSat();
    x[0]->unsetBList();
    x[0]->unsetIList();
    x[0]->unsetRange();
    x[1]->unsetBList();
    x[1]->unsetIList();
    x[1]->unsetRange();
    pred->unsetBList();
    pred->unsetIList();
    pred->unsetRange();
  }
}

void BuildObjectMod::print(std::ostream& o, const BuildObjectPredicate *pred) const
{
  o << "(";
  pred->scope[0]->print( o );
  o << " % ";
  pred->scope[1]->print( o );
  o << ")";
}

string BuildObjectMod::toString(const BuildObjectPredicate* pred) const 
{
  string x1 = pred->scope[0]->toString();
  string x2 = pred->scope[1]->toString();
  string s = ( "(" + x1 + " % " + x2 + ")" );
  return s;
}

string BuildObjectMod::xmlPred(int& idx, int level, const BuildObjectPredicate *pred, int *sub) const 
{
  string x1 = pred->scope[0]->xmlPred(idx, level+1, sub);
  string x2 = pred->scope[1]->xmlPred(idx, level+1, sub);
  string s = ( "mod(" + x1 + "," + x2 + ")" );
  return s;
}


/**********************************************
 *  Member Predicate BuildObject
 **********************************************/ 
void BuildObjectMember::build( Solver *s, VariableInt **tmp, BuildObjectPredicate *pred ) 
{
  if( !pred->isConstant() )
    new PredicateMember( s, tmp, &(pred->params[1]), pred->params[0], 1 );
}

int BuildObjectMember::propagateUpward( BuildObjectPredicate *pred ) const 
{
  BuildObject *x = pred->scope[0];
  int i, n = pred->params[0], inter_size=0;
  int *vals = &(pred->params[1]);
  for(i=0; i<n; ++i)
    if( x->contain(vals[i]) ) 
      ++inter_size;
  int lb = inter_size == x->size();
  int ub = inter_size > 0;
  int consistent = ( pred->setMin( lb ) && pred->setMax( ub ) );
  if( consistent && pred->isTop() ) 
    consistent = pred->setDomain(1);
  return consistent;
}

int BuildObjectMember::propagateDownward( BuildObjectPredicate *pred ) const 
{
  bool consistent = true;
  BuildObject *x = pred->scope[0];
  int i, inter_size=0, n=pred->params[0], *vals=&(pred->params[1]);
  if( pred->equal(0) ) {
    for(i=0; i<n; ++i)
      if( x->contain(vals[i]) ) 
	++inter_size;
    if( inter_size == x->size() )
      consistent = false;
    else if( inter_size == x->size()-1 ) 
      x->removeSet( vals, n );
  } else if( pred->equal(1) ) {
    consistent = x->setDomain( vals, n );
  }
  return consistent;
}

void BuildObjectMember::close( BuildObjectPredicate *pred ) 
{  
  BuildObject *x = pred->scope[0];
  bool valid = ( pred->equal(1) );
  if( pred->equal(0) ) { 
    int i, inter_size=0, n=pred->params[0], *vals=&(pred->params[1]);
    for(i=0; i<n; ++i)
      if( x->contain(vals[i]) ) 
	++inter_size;
    valid = !inter_size;
  }
  if( valid ) {
    pred->unsetRel();
    x->deReference();
    if(x->size() != (x->max()-x->min()+1))
      x->unsetRange();
  } else {
    pred->model->unsetSat();
    x->unsetBList();
    x->unsetIList();
    x->unsetRange();
    pred->unsetBList();
    pred->unsetIList();
  }

}

void BuildObjectMember::print(std::ostream& o, const BuildObjectPredicate *pred) const
{
  o << "(";
  pred->scope[0]->print( o );
  o << " in {" ;
  for(int i=0; i<pred->params[0]-1; ++i)
    o << (pred->params[1+i]) << ",";
  o << pred->params[pred->params[0]] << "}";
}

string BuildObjectMember::toString(const BuildObjectPredicate* pred) const 
{
  string x1 = pred->scope[0]->toString();
  string s = ( "(" + x1 + " in " + int2string(pred->params[0]) + ")" );
  return s;
}

 string BuildObjectMember::xmlPred(int& idx, int level, const BuildObjectPredicate *pred, int *sub) const 
 {
   string x1 = pred->scope[0]->xmlPred(idx, level+1, sub);
   string s = ( "member(" + x1 + "," + int2string(pred->params[0]) + ")" );
   return s;
 }


/**********************************************
 * Element Constraint BuildObject
 **********************************************/
void BuildObjectElement::build( Solver *s, VariableInt **tmp, BuildObjectPredicate *pred ) 
{
  new PredicateElement( s, tmp, pred->arity+1, pred->params[0] );
}

int BuildObjectElement::propagateUpward( BuildObjectPredicate *pred ) const 
{
  BuildObject **x = pred->scope;
  int n=pred->arity-1, k=pred->params[0];
  int ub = -NOVAL;
  int lb = NOVAL;
  for(int i=0; i<n; ++i) {
    if( x[n]->contain(i+k) ) {
      if( ub < x[i]->max() ) ub = x[i]->max(); 
      if( lb > x[i]->min() ) lb = x[i]->min(); 
    }
  }
  return ( pred->setMin( lb ) && pred->setMax( ub ) );
}

int BuildObjectElement::propagateDownward( BuildObjectPredicate *pred ) const 
{
  BuildObject **x = pred->scope;
  int // lb=NOVAL, ub=-NOVAL,
    i, offset=pred->params[0], n=pred->arity-1;  
  bool consistent = true;
  for(i=0; consistent && i<n; ++i)
    if( !pred->intersect( x[i] ) )
      consistent = x[n]->remove( i+offset );
  if( consistent && x[n]->isConstant() ) {
    i = x[n]->value();
    consistent = pred->setEqual( x[i-offset] );
  }
  return consistent;
}

void BuildObjectElement::close( BuildObjectPredicate *pred ) 
{  
  BuildObject **x = pred->scope;
  int i, n=pred->arity-1;  
  if( x[n]->isConstant() ) {
    pred->unsetRel();
    for(i=0; i<=n; ++i) 
      x[i]->deReference();
  } else {
    pred->model->unsetSat();
    for(i=0; i<=n; ++i) {
      x[i]->unsetBList();
      x[i]->unsetIList();
      x[i]->unsetRange();
    }
    pred->unsetBList();
    pred->unsetIList();
    pred->unsetRange();
  }
}

void BuildObjectElement::print(std::ostream& o, const BuildObjectPredicate *pred) const
{
  o << "T[";
  pred->scope[pred->arity-1]->print( o );
  o << "], T={ ";
  for(int i=0; i<pred->arity-1; ++i) {
    pred->scope[i]->print( cout );
    cout << " ";
  }
  o << "}";
}

string BuildObjectElement::toString(const BuildObjectPredicate* pred) const 
{
  string x1 = pred->scope[pred->arity-1]->toString();
  string s = ( "X[" + x1 + "]" );
  return s;
}

string BuildObjectElement::xmlPred(int& idx, int level, const BuildObjectPredicate *pred, int *sub) const 
{
  string s;
  if( !level ) 
    s = "global:element";
  else
    s = ("X" + int2string(idx++)); 
  return s;
}


/**********************************************
 * Sum Constraint BuildObject
 **********************************************/
void BuildObjectSum::self_build( Solver *s, VariableInt **tmp, BuildObjectPredicate *pred )
{
  //if( n == 1 && !weighted ) {

  //BuildObject **x = pred->scope;
  int *w=pred->params+1;
  int n = w[-1];
  //int wpos = w[n+2], wneg = w[n+3];
 
  /// check if we can use an offset predicate, or a less/more constraint:
  if( // false &&
     pred->size() > 1 && n == 2 && w[0] == 1 && w[1] == -1 && w[2] == w[3] ) {
    //cout << "SELF BUILD" << endl;    
    
    if( !w[2] )
      {
	cout << "equality!" << endl;
	exit(0);
      }

      VariableVirtual *vv;
      s->MISC += 1;
      vv = new VariableVirtual( s );
      if(pred->veqptr_) {
	pred->veqptr_->varptr_ = vv;
	pred->veqptr_->setVirtual();
      } else {
	pred->varptr_ = vv;
      }
// 	((VariableVirtual*)(pred->varptr_))->reference = tmp[1];
// 	((VariableVirtual*)(pred->varptr_))->conversion = new MapFactor( x[0]->min(), ((VariableVirtual*)(pred->varptr_)) );
      vv->reference = tmp[0];
      vv->conversion = new MapOffset( -w[2], vv );
      pred->setVirtual();
      //pred->veqptr_ = x[1];

//     ++s->MISC;
//     //pred->veqptr_ = x[0];
//     pred->setVirtual();
//     pred->varptr_ = new VariableVirtual( s );    
//     ((VariableVirtual*)(pred->varptr_))->reference = tmp[0];
//     ((VariableVirtual*)(pred->varptr_))->conversion = new MapOffset( -w[2], ((VariableVirtual*)(pred->varptr_)) );
  } else {  
    //cout << "build (standard)" << endl;
    pred->BuildObject::build( s );
  }
  //return X;
}

void BuildObjectSum::build( Solver *s, VariableInt **tmp, BuildObjectPredicate *pred ) 
{
  //BuildObject **x = pred->scope;
  int *w=pred->params+1;
  int n = w[-1];
  int wpos = w[n+2], wneg = w[n+3];
 
 Constraint *cons;

 /// check if we can use an offset predicate, or a less/more constraint:
 if( n == 2 && w[0] == 1 && w[1] == -1 )
   {
     if( w[2] == w[3] ) {
       if( pred->size() == 1 )
	 cons = new PredicateOffset( s, tmp, -w[2] );
     } else {
       if( w[2] != NOVAL/2 ) {
	 cons = new ConstraintLess( s, tmp, -w[2] );
       } if( w[3] != -NOVAL/2 ) {
	 VariableInt *v_aux = tmp[0];
	 tmp[0] = tmp[1];
	 tmp[1] = v_aux;	    
	 cons = new ConstraintLess( s, tmp, w[3] );
       }
     }
   }
 else if( n == 3 && w[0] == 1 && w[2] == -1 && !w[3] && !w[4] && (w[1] * w[1] == 1) )
   {
     if( w[1] == 1 )
       cons = new PredicateAdd( s, tmp );
     else
       cons = new PredicateSub( s, tmp );
   }
 else if( w[pred->arity+5] && wpos == n ) {
   //cout << "ConstraintBool" << endl;
   if( w[n] == w[n+1] ) {
     cons = new ConstraintBoolSumEqual( s, tmp, n, w[n] );
   } else {
     if( w[n] != NOVAL/2 ) {
       //cons = new ConstraintWeightedSum( s, tmp, n, w[n], w[n+1], w, wpos, wneg );
       cons = new ConstraintBoolSumLess( s, tmp, n, w[n] );
     } if( w[n+1] != -NOVAL/2 ) {
       //cons = new ConstraintWeightedSum( s, tmp, n, w[n], w[n+1], w, wpos, wneg );
       cons = new ConstraintBoolSumMore( s, tmp, n, w[n+1] );
     }
   }
  } else {
   //cout << w[n] << ".." << w[n+1] << endl;
   //cout << "ConstraintWeightedSum" << endl;
   if(w[n] != NOVAL/2 || w[n+1] != -NOVAL/2)
     cons = new ConstraintWeightedSum( s, tmp, n, w[n], w[n+1], w, wpos, wneg );
 }
//   cons->print( cout );
//   cout << endl << endl;

}

int BuildObjectSum::propagateUpward( BuildObjectPredicate *pred ) const 
{

  BuildObject **x = pred->scope;
  int *w = pred->params+1, n=pred->arity;

  //std::cout << n << std::endl;

  w[n+2] = 0;
  w[n+3] = 0;
  for(int i=0; i<n; ++i) {
    //std::cout << (w[i]) << " ";
    //x[i]->print( std::cout );
    //std::cout << " ";

    w[n+2] += (w[i] < 0 ? w[i]*x[i]->min() : w[i]*x[i]->max() );
    w[n+3] += (w[i] < 0 ? w[i]*x[i]->max() : w[i]*x[i]->min() );

    //std::cout << "[" << w[n+3] << ".." << w[n+2] << "]" << std::endl;

  }
  w[n+3] -= w[n+1];
  w[n+2] -= w[n+1];

  //std::cout << "[" << w[n+3] << ".." << w[n+2] << "]" << std::endl;

  return ( pred->setMin( w[n+3] ) && pred->setMax( w[n+2] ) );
}

int BuildObjectSum::propagateDownward( BuildObjectPredicate *pred ) const
{


//   w[n+2] = 0;
//   w[n+3] = 0;
//   for(int i=0; i<n; ++i) {
//     std::cout << (w[i]) << " ";
//     x[i]->print( std::cout );
//     std::cout << " ";

//     w[n+2] += (w[i] < 0 ? w[i]*x[i]->min() : w[i]*x[i]->max() );
//     w[n+3] += (w[i] < 0 ? w[i]*x[i]->max() : w[i]*x[i]->min() );

//     std::cout << "[" << w[n+3] << ".." << w[n+2] << "]" << std::endl;

//   }
//   w[n+3] -= w[n+1];
//   w[n+2] -= w[n+1];

//   

  //if( propagateUpward(pred) ) {

    BuildObject **x = pred->scope;
    int consistent = true, *w = pred->params+1, n = pred->arity, i;

    int lb = w[n+3], ub = w[n+2];

    //std::cout << (pred->min()) << " " << (pred->max()) << std::endl;
    
    //std::cout << lb << " " << ub << std::endl;
    
    if( lb != pred->min() || ub != pred->max() ) {
      
      
      //lb -= pred->max();
      
      //    cout << "[" << lb << "," << ub << "]  [" << (pred->min()) << "," << (pred->max()) << "]" << endl; 
      
      lb = pred->max() - lb;
      ub -= pred->min();
      int aux;
      for(i=0; consistent && i<n; ++i) {
	if( w[i] > 0 ) {
	  aux = x[i]->max();
	  
	  // 	cout << (aux - ub/w[i]) << " <= " ;
	  // 	x[i]->print( cout );
	  // 	cout << " <= " << (x[i]->min() + lb/w[i]) << endl;
	  
	  
	  consistent = ( (x[i]->setMax(x[i]->min() + lb/w[i])) 
			 &&
			 (x[i]->setMin(aux - ub/w[i])) );
	  
	} else if( w[i] < 0 ) {
	  aux = x[i]->min();
	  
	  
	  // 	cout << (x[i]->max() + lb/w[i]) << " <= " ;
	  // 	x[i]->print( cout );
	  // 	cout << " <= " << (aux - ub/w[i]) << endl;
	  
	  
	  consistent = ( (x[i]->setMin(x[i]->max() + lb/w[i]))
			 &&
			 (x[i]->setMax(aux - ub/w[i])) );
	}
      }
    }
    return consistent;
    //} 
    //return false;
}

void BuildObjectSum::close( BuildObjectPredicate *pred ) 
{  
  // before close:
  // w[0]..w[n-1] -> weights



  BuildObject **x = pred->scope, *v_aux;
  int n=pred->arity, i=n, k, *w=pred->params+1, i_aux, j, m=n;
  w[n+5]=1;

//   cout << pred->id << " " << endl;
//   for(int y=0; y<=m; ++y) {
//     cout << "  ";
//     x[y]->BuildObject::printshort( cout );
//     cout << " ";
//   }
//   cout << endl;
//   for(int y=0; y<m+6; ++y)
//     cout << setw(4) << w[y] << " ";
//   cout << endl;

  while( i-- ) {
    if( x[i]->isConstant() || !w[i] ) {
      /// remove constant or null-weighted vars
      if( w[i] )
	w[m+1] -= (w[i] * x[i]->value());
      x[i]->deReference();
      x[i] = x[--n];
      w[i] = w[n];
    } else {
      /// aggregate multiple variables
      for( k=i+1; k<n; ++k)
	if(x[i]->getVarId() == x[k]->getVarId()) {
	  w[k] += w[i];
	  if( !w[k] ) {
	    x[k] = x[--n];
	    w[k] = w[n];
	  }
	  x[i] = x[--n];
	  w[i] = w[n];
	}
    }
  }

  /// find out if the sum is weighted and/or Boolean
  for(i=0; w[m+5] && i<n; ++i) 
    w[m+5] &= x[i]->isBoolean();

  /// reset the arity
  x[n] = x[m];
  for(i=0; i<6; ++i)
    w[n+i] = w[m+i];
  pred->arity = n; 


//   for(int y=0; y<=m; ++y) {
//     cout << "  ";
//     x[y]->BuildObject::printshort( cout );
//     cout << " ";
//   }
//   cout << endl;
//   for(int y=0; y<m+6; ++y)
//     cout << setw(4) << w[y] << " ";
//   cout << endl;


  /// equality
  if( n == 1 && !w[n+1] && w[0] == 1 ) {

    //cout << "equality!" << endl;
    //exit(0);

    pred->setEqual( x[0] );
    x[0]->deReference();
    n=0;
  }

  //  assert( w[n] == -1 );

  /// classic closure
  if(!n) pred->unsetRel();  
  else {
    pred->model->unsetSat();
    if(pred->isConstant()) pred->deReference();
    for(i=0; i<n; ++i) {
      x[i]->unsetBList();
      x[i]->unsetIList();
    }
    pred->unsetBList();
    pred->unsetIList();

    //cout << n << endl;
    // assert( w[n] == -1 );
    if( pred->isReferenced() < 1 || pred->isConstant() ) {
      // The reified var is not referenced by another constraint, or it is constant.
      // We can get rid of it
      if( w[n+2] > pred->max() ) // need an upper bound;
	w[n] = pred->max()+w[n+1];
      else w[n] = NOVAL/2;
      if( w[n+3] < pred->min() ) // need an upper bound;
	w[n+1] = pred->min()+w[n+1];
      else w[n+1] = -NOVAL/2;
      if( pred->isReferenced() == 1 )
	pred->deReference();

    } else {
      ++n; 
      w[n+1] = w[n];

      //for(i=4; i>=0; --i)
      //w[n+i+1] = w[n+i];

    }
    /// set the "real" arity, that is, counting the predicate itself if needed.
    w[-1] = n;


//     for(int y=0; y<=m; ++y) {
//       cout << "  ";
//       x[y]->BuildObject::printshort( cout );
//       cout << " ";
//     }
//     cout << endl;
//     for(int y=0; y<m+6; ++y)
//       cout << setw(4) << w[y] << " ";
//     cout << endl;
    

    i=0;
    j=n-1;
    // we put it in a normal form:
    // 1/ all vars with unit weight, 2/ all vars with weight>1, 3/ all vars with weight < 0 
    // w[0]..w[n-1] keep the weights, w[n]/w[n+1] keep ub/lb, w[n+2]/w[n+3] keep wpos/wneg, w[n+4] keeps 'bool'
    while( true )
      {
	while( i<n && w[i] == 1 ) ++i;
	while( j && w[j] != 1 ) --j;
	if( i >= j ) break;
	else {
	  i_aux  = w[i];
	  v_aux  = x[i];
	  w[i]   = w[j];
	  x[i] = x[j];
	  w[j]   = i_aux;
	  x[j] = v_aux;
	}
      } 
    
    w[n+2] = i;
    j = n-1;
    while( true )
      {
	while( i<n && w[i] > 1 ) ++i;
	while( j>=w[n+2] && w[j] < 0 ) --j;

	if( i >= j ) break;
	else {
	  i_aux  = w[i];
	  v_aux  = x[i];
	  w[i]   = w[j];
	  x[i] = x[j];
	  w[j]   = i_aux;
	  x[j] = v_aux;
	}
      }
    w[n+3] = i;
  }


//   for(int y=0; y<=m; ++y) {
//     cout << "  ";
//     x[y]->BuildObject::printshort( cout );
//     cout << " ";
//   }
//   cout << endl;
//   for(int y=0; y<m+6; ++y)
//     cout << setw(4) << w[y] << " ";
//   cout << endl << (pred->isReferenced()) << endl << endl;

}

void BuildObjectSum::print(std::ostream& o, const BuildObjectPredicate *pred) const
{
  o << "(";
  o.flush();

  //o << pred->arity << std::endl;

  int * w = pred->params+1;
  if( pred->params ){
    if( w[0] == -1 )
      o << "-";
    else if( w[0] != 1 )
      o << w[0] << "*";
  }
  for(int i=0; i<pred->arity-1; ++i) {
    pred->scope[i]->print( o );

    if( w ){
      if( w[i+1] == -1 )
	o << " - ";
      else if( w[i+1] < 0 )
	o << " - " << (-1 * w[i+1]) << "*";
      else if( w[i+1] > 1 )
	o << " + " << w[i+1] << "*";
      else o << " + " ;
    } else 
      o << " + " ;
  }
  pred->scope[pred->arity-1]->print( o );
  o << " + " << ( !w ? 0 : w[pred->arity+1]) << ")";
    //pred->BuildObject::print(o);
}

string BuildObjectSum::toString(const BuildObjectPredicate* pred) const 
{
  //int * w = pred->params+1;
  string s;
  int * w = pred->params+1;
  if( pred->params ){
    if( w[0] == -1 )
      s += "-";
    else if( w[0] != 1 )
      s += (int2string(w[0]) + "*");
  }
  for(int i=0; i<pred->arity-1; ++i) {
    s += pred->scope[i]->toString();
    if( w ){
      if( w[i+1] == -1 )
	s += " - ";
      else if( w[i+1] < 0 )
	s += ( " - " + int2string(-1 * w[i+1]) + "*" );
      else if( w[i+1] > 1 )
	s += ( " + " + int2string(w[i+1]) + "*" );
      else s += " + " ;
    } else 
      s += " + " ;
  }
  s += pred->scope[pred->arity-1]->toString();
  return s;
}

string BuildObjectSum::xmlPred(int& idx, int level, const BuildObjectPredicate *pred, int *sub) const 
{
  string s;
  if( !level ) 
    s = "global:weightedSum";
  else
    s = ("X" + int2string(idx++)); 
  return s;
}


/**********************************************
 * Max Constraint BuildObject
 **********************************************/
void BuildObjectMax::build( Solver *s, VariableInt **tmp, BuildObjectPredicate *pred ) 
{
  new PredicateMax( s, tmp, pred->arity+1 );
}

int BuildObjectMax::propagateUpward( BuildObjectPredicate *pred ) const 
{
  BuildObject **x = pred->scope;
  int i, j, n=pred->arity;
  int ub = x[0]->max();
  int lb = x[0]->min();
  for(i=1; i<n; ++i) {
    j = x[i]->min();
    if(lb < j) lb = j;
    j = x[i]->max();
    if(ub < j) ub = j;
  }
  return ( pred->setMin( lb ) && pred->setMax( ub ) );
}

int BuildObjectMax::propagateDownward( BuildObjectPredicate *pred ) const
{
  int consistent = true;
  BuildObject **x = pred->scope;
  int i, n=pred->arity;
  int ub = pred->max();
  for(i=0; consistent && i<n; ++i) {
    consistent = x[i]->setMax(ub);
  }
  return consistent;
}

void BuildObjectMax::close( BuildObjectPredicate *pred ) 
{  
  BuildObject **x = pred->scope;
  int n=pred->arity, i=n;
  while( i-- )
    if( x[i]->max() < pred->min() ) {
      x[i]->deReference();
      x[i] = x[--n];
    }
  pred->arity = n;
  if(!n) pred->unsetRel();
  else {
    pred->model->unsetSat();
    for(i=0; i<n; ++i) {
      x[i]->unsetBList();
      x[i]->unsetIList();
    }
    pred->unsetBList();
    pred->unsetIList();
  }
}

void BuildObjectMax::print(std::ostream& o, const BuildObjectPredicate *pred) const
{
  o << pred->arity;
  o << " max(";
  for(int i=0; i<pred->arity-1; ++i) {
    pred->scope[i]->print( o );
    o << ", ";
  }
  pred->scope[pred->arity-1]->print( o );
  o << ") " ;
}

string BuildObjectMax::toString(const BuildObjectPredicate* pred) const 
{
  //int * w = pred->params+1;
  string s = "Max(";
  for(int i=0; i<pred->arity-1; ++i) {
    s += pred->scope[i]->toString();
    s += ", " ;
  }
  s += pred->scope[pred->arity-1]->toString();
  s += ")";
  return s;
}

string BuildObjectMax::xmlPred(int& idx, int level, const BuildObjectPredicate *pred, int *sub) const 
{
  string s;
  if( !level ) 
    s = "global:max";
  else
    s = ("X" + int2string(idx++)); 
  return s;
}


/**********************************************
 * Min Constraint BuildObject
 **********************************************/
void BuildObjectMin::build( Solver *s, VariableInt **tmp, BuildObjectPredicate *pred ) 
{
  new PredicateMin( s, tmp, pred->arity+1 );
}

int BuildObjectMin::propagateUpward( BuildObjectPredicate *pred ) const 
{
  BuildObject **x = pred->scope;
  int i, j, n=pred->arity;
  int ub = x[0]->max();
  int lb = x[0]->min();
  for(i=1; i<n; ++i) {
    j = x[i]->min();
    if(lb > j) lb = j;
    j = x[i]->max();
    if(ub > j) ub = j;
  }
  return ( pred->setMin( lb ) && pred->setMax( ub ) );
}

int BuildObjectMin::propagateDownward( BuildObjectPredicate *pred ) const
{
  int consistent = true;
  BuildObject **x = pred->scope;
  int i, n=pred->arity;
  int lb = pred->min();
  for(i=0; consistent && i<n; ++i) {
    consistent = x[i]->setMin(lb);
  }
  return consistent;
}

void BuildObjectMin::close( BuildObjectPredicate *pred ) 
{  
  BuildObject **x = pred->scope;
  int n=pred->arity, i=n;
  while( i-- )
    if( x[i]->min() > pred->max() ) {
      x[i]->deReference();
      x[i] = x[--n];
    }
  pred->arity = n;
  if(!n) pred->unsetRel();
  else {
    pred->model->unsetSat();
    for(i=0; i<n; ++i) {
      x[i]->unsetBList();
      x[i]->unsetIList();
    }
    pred->unsetBList();
    pred->unsetIList();
  }
}

void BuildObjectMin::print(std::ostream& o, const BuildObjectPredicate *pred) const
{
  o << pred->arity;
  o << " min(";
  for(int i=0; i<pred->arity-1; ++i) {
    pred->scope[i]->print( o );
    o << ", ";
  }
  pred->scope[pred->arity-1]->print( o );
  o << ") " ;
}

string BuildObjectMin::toString(const BuildObjectPredicate* pred) const 
{
  //int * w = pred->params+1;
  string s = "Min(";
  for(int i=0; i<pred->arity-1; ++i) {
    s += pred->scope[i]->toString();
    s += ", " ;
  }
  s += pred->scope[pred->arity-1]->toString();
  s += ")";
  return s;
}

string BuildObjectMin::xmlPred(int& idx, int level, const BuildObjectPredicate *pred, int *sub) const 
{
  string s;
  if( !level ) 
    s = "global:min";
  else
    s = ("X" + int2string(idx++)); 
  return s;
}

/**********************************************
 * AllDiff Constraint BuildObject
 **********************************************/
void BuildObjectAllDiff::build( Solver *s, VariableInt **tmp, BuildObjectPredicate *pred ) 
{
  int c = pred->params[0];
  if((!(c & 1) || pred->arity > 2) && (c & 2))
    new ConstraintAllDiff( s, tmp, pred->arity );
}

void BuildObjectAllDiff::add( CSP *mod, const int l, 
			      BuildObjectPredicate* pred ) 
{
  int i, j, n=pred->arity, m;
  BuildObject **x = pred->scope;

  /// clique of notEqual
  if(pred->params[0] & 1) {
    for(i=0; i<n; ++i)
      for(j=i+1; j<n; ++j)
	{	    
	  //mod->add( CSP::_Equal( x[i], x[j], 0 ) );
	  BuildObject* p = CSP::_Equal( x[i], x[j], 0 );
	  p->add( mod, 0 );
	}
  }

  /// Bound Consistency
  if( pred->params[0] & 2 ) {
    BuildObjectConstraint::add(mod, l, pred);
  } 


  /// Bound Consistency decomposition
  if( pred->params[0] & 4 ) {    
    int MaxHall = pred->params[1], lim;
    int k, minb=NOVAL, maxb=-NOVAL, itv;
    int lb[n], ub[n];
    BuildObject  ***b = new BuildObject**[n];
    BuildObject ****a = new BuildObject***[n];
    for(i=0; i<n; ++i) {
      lb[i] = x[i]->min();
      ub[i] = x[i]->max();
      if(lb[i] < minb ) minb = lb[i];
      if(ub[i] > maxb ) maxb = ub[i];
      m = (ub[i] - lb[i] + 1);
      b[i] = new BuildObject*[m];
      a[i] = new BuildObject**[m];
      b[i]-=lb[i];
      a[i]-=lb[i];

      for(j=lb[i]; j<=ub[i]; ++j) 
	b[i][j] = CSP::_Precedence( x[i], j );

      for(j=lb[i]; j<=ub[i]; ++j) {
	a[i][j] = new BuildObject*[m];
	a[i][j]-=lb[i];

	lim = std::min(ub[i], j+MaxHall-1);
	//for(k=j; k<=ub[i]; ++k) {
	for(k=j; k<=lim; ++k) {

#ifdef _DEBUGDECOMP
	  cout << "a[" << i << "][" << j << "][" << k << "]" << " ";
#endif

	  if(j==lb[i] && k==ub[i])
	    a[i][j][k] = CSP::_Variable(1,1);
	  else if(j==lb[i])
	    a[i][j][k] = b[i][k];
	  else if(k==ub[i]) {
	    // beware, we use in fact 1-b[i][j-1];
	    a[i][j][k] = b[i][j-1];

#ifdef _DEBUGDECOMP
	    cout << "1 - ";
#endif

	  }
	  else
	    a[i][j][k] = CSP::_Precedence( b[i][j-1], 1, b[i][k] );
	  
#ifdef _DEBUGDECOMP
	  a[i][j][k]->print( cout );
	  cout << endl;
#endif

	}
      }
    }

#ifdef _DEBUGDECOMP
    cout << endl;
    for(i=0; i<n; ++i) {
      cout << "a[" << i << "][" << maxb << "][" << maxb << "]" << " ";
      a[i][maxb][maxb]->print( cout );
      cout << endl;
    }
    cout << endl;
    cout << ((maxb-minb+1) * (maxb - minb + 2) / 2) << endl;
#endif

    for(j=minb; j<maxb; ++j) {
      lim = std::min(maxb, j+MaxHall);

      for(k=j; k<lim; ++k) {
	if(k-j <= MaxHall) {
	  //for(k=j; k<maxb; ++k) 
	  //	if(k-j+1 < 4) {
	  itv = 0;
	  
#ifdef _DEBUGDECOMP
	  cout << endl << " A_" << j << "." << k << " ";
#endif
	  
	  BuildObject **alu = new BuildObject*[n+1];
	
	  for(i=0; i<n; ++i)
	    if(lb[i] <= j && ub[i] >= k) {
	    
#ifdef _DEBUGDECOMP
	      a[i][j][k]->print( cout );
	      cout << " + " ;
#endif
	    
	      alu[itv++] = a[i][j][k];
	    }  
	  
#ifdef _DEBUGDECOMP
	  cout << " <= " << (k-j+1) << endl;
#endif
	
	  mod->add( CSP::_Precedence( CSP::_Sum( alu, itv ), (k-j+1) ) );
	}
      }
      
      if( maxb-j < MaxHall ) {
#ifdef _DEBUGDECOMP
	cout << endl << " A_" << j << "x" << maxb << " ";
#endif
	
	BuildObject **alu = new BuildObject*[n+1];
	itv = 0;
	for(i=0; i<n; ++i)
	  if(lb[i] <= j && ub[i] >= maxb) {
	    
#ifdef _DEBUGDECOMP
	    a[i][j][maxb]->print( cout );
	    cout << " + " ;
#endif
	    
	    alu[itv++] = a[i][j][maxb];
	  }  
	
#ifdef _DEBUGDECOMP
	cout << " >= " << (n-k+j) << endl;
#endif
	
	mod->add( CSP::_Precedence( (n-k+j), CSP::_Sum( alu, itv ) ) );
      }
    }

#ifdef _DEBUGDECOMP
    cout << endl << " A_" << maxb << "," << maxb << " ";
#endif

    BuildObject **alu = new BuildObject*[n+1];
    itv = 0;
    for(i=0; i<n; ++i)
      if(lb[i] <= maxb && ub[i] >= maxb) {

#ifdef _DEBUGDECOMP
	a[i][maxb][maxb]->print( cout );
	cout << " + " ;
#endif

	alu[itv++] = a[i][maxb][maxb];
      }  

#ifdef _DEBUGDECOMP
    cout << " >= " << (n-1) << endl;
#endif

    mod->add( CSP::_Precedence( (n-1), CSP::_Sum( alu, itv ) ) );
  }
}

void BuildObjectAllDiff::close( BuildObjectPredicate *pred ) 
{  
  BuildObject **x = pred->scope;
  int n=pred->arity, i=n;  
  if( pred->params[0] & 2 ) {
    while( i-- )
      if( x[i]->isConstant() ) {
	x[i]->deReference();
	x[i] = x[--n];
      }
    pred->arity = n;
    if(!n) pred->unsetRel();
    else {
      pred->model->unsetSat();
      pred->deReference();
      for(i=0; i<n; ++i) {
	x[i]->unsetBList();
	x[i]->unsetIList();
      }
      pred->unsetBList();
      pred->unsetIList();
    }
  }
//  else {
//     pred->unsetRel();
//     pred->deReference();
//   }
}

void BuildObjectAllDiff::print(std::ostream& o, const BuildObjectPredicate *pred) const
{
  o << "AllDifferent(";
  for(int i=0; i<pred->arity-1; ++i) {
    pred->scope[i]->print( o );
    o << " " ;
  }
  pred->scope[pred->arity-1]->print( o );
  o << ")";
}

string BuildObjectAllDiff::toString(const BuildObjectPredicate* pred) const 
{
  //int * w = pred->params+1;
  string s = "AllDiff(";
  for(int i=0; i<pred->arity-1; ++i) {
    s += pred->scope[i]->toString();
    s += ", " ;
  }
  s += pred->scope[pred->arity-1]->toString();
  s += ")";
  return s;
}

string BuildObjectAllDiff::xmlPred(int& idx, int level, const BuildObjectPredicate *pred, int *sub) const 
{
  string s;
  if( !level ) 
    s = "global:allDifferent";
  else
    s = ("X" + int2string(idx++)); 
  return s;
}


/**********************************************
 * Clause Constraint BuildObject
 **********************************************/
void BuildObjectClause::build( Solver *s, VariableInt **tmp, BuildObjectPredicate *pred ) 
{
  if( !s->sat )
    s->sat = new ConstraintClauseBase(s);
  s->sat->addDelayedClause(tmp, pred->arity, pred->params);
}

void BuildObjectClause::close( BuildObjectPredicate *pred ) 
{
  pred->deReference();
}

void BuildObjectClause::print(std::ostream& o, const BuildObjectPredicate *pred) const
{
  o << "Clause(";
  for(int i=0; i<pred->arity-1; ++i) {
    pred->scope[i]->print( o );
    o << " " ;
  }
  pred->scope[pred->arity-1]->print( o );
  o << ")";
}

string BuildObjectClause::toString(const BuildObjectPredicate* pred) const 
{
  //int * w = pred->params+1;
  string s = "Clause(";
  for(int i=0; i<pred->arity-1; ++i) {
    s += pred->scope[i]->toString();
    s += ", " ;
  }
  s += pred->scope[pred->arity-1]->toString();
  s += ")";
  return s;
}

string BuildObjectClause::xmlPred(int& idx, int level, const BuildObjectPredicate *pred, int *sub) const 
{
  string s;
  if( !level ) 
    s = "global:clause";
  else
    s = ("X" + int2string(idx++)); 
  return s;
}


/**********************************************
 * TDAG Constraint BuildObject
 **********************************************/
void BuildObjectTDAG::build( Solver *s, VariableInt **tmp, BuildObjectPredicate *pred ) 
{  
//   VariableInt **t1 = new VariableInt[pred->arity];
//   VariableInt **t2 = new VariableInt[pred->arity];
  VariableInt *t1[pred->arity];
  VariableInt *t2[pred->arity];

  for(int i=0; i<pred->arity; ++i) {
    t1[i] = ((BuildObjectPredicate*)(pred->scope[i]))->scope[0]->getVariable();
    t2[i] = ((BuildObjectPredicate*)(pred->scope[i]))->scope[1]->getVariable();
  }

  new ConstraintTDAG( s, tmp, pred->arity, t1, t2 );
}

void BuildObjectTDAG::close( BuildObjectPredicate *pred ) 
{  
  pred->deReference();
}

void BuildObjectTDAG::print(std::ostream& o, const BuildObjectPredicate *pred) const
{
  o << "TDAG(";
  for(int i=0; i<pred->arity-1; ++i) {
    pred->scope[i]->print( o );
    o << " " ;
  }
  pred->scope[pred->arity-1]->print( o );
  o << ")";
}

string BuildObjectTDAG::toString(const BuildObjectPredicate* pred) const 
{
  //int * w = pred->params+1;
  string s = "TDAG(";
  for(int i=0; i<pred->arity-1; ++i) {
    s += pred->scope[i]->toString();
    s += ", " ;
  }
  s += pred->scope[pred->arity-1]->toString();
  s += ")";
  return s;
}

string BuildObjectTDAG::xmlPred(int& idx, int level, const BuildObjectPredicate *pred, int *sub) const 
{
  string s = "global:tdag";
  return s;
}

/**********************************************
 * Tree Constraint BuildObject
 **********************************************/
void BuildObjectTree::build( Solver *s, VariableInt **tmp, BuildObjectPredicate *pred ) 
{  
  int n=pred->arity;
  if(pred->params[0]) {
    new ConstraintTreeExplicit( s, tmp, pred->params[0] );
  } else
    new ConstraintTree( s, tmp, n );
}

void BuildObjectTree::print(std::ostream& o, const BuildObjectPredicate *pred) const
{
  o << "Tree(";
  for(int i=0; i<pred->arity-1; ++i) {
    pred->scope[i]->print( o );
    o << " " ;
  }
  pred->scope[pred->arity-1]->print( o );
  o << ")";
}

string BuildObjectTree::toString(const BuildObjectPredicate* pred) const 
{
  //int * w = pred->params+1;
  string s = "Tree(";
  for(int i=0; i<pred->arity-1; ++i) {
    s += pred->scope[i]->toString();
    s += ", " ;
  }
  s += pred->scope[pred->arity-1]->toString();
  s += ")";
  return s;
}

string BuildObjectTree::xmlPred(int& idx, int level, const BuildObjectPredicate *pred, int *sub) const 
{
  string s = "global:tree";
  return s;
}

/**********************************************
 * Gcc Constraint BuildObject
 **********************************************/
void BuildObjectGcc::build( Solver *s, VariableInt **tmp, 
			    BuildObjectPredicate *pred )
{
  int n=pred->arity, m = (pred->params[1]-pred->params[0]+1);
  
    new ConstraintGlobalCardinality( s, tmp, n,
				     pred->params[0],
				     pred->params[1],
				     &(pred->params[2]),
				     &(pred->params[2+m]) );
}

void BuildObjectGcc::add( CSP *mod, const int l, 
			  BuildObjectPredicate* pred ) 
{
  assert( l == 0 );
  BuildObjectConstraint::add(mod, l, pred);
  BuildObject **x = pred->scope;
  int i, j, n=pred->arity, m = (pred->params[1]-pred->params[0]+1);
  int max_occ, min_occ, value;
  for(j=0; j<m; ++j) {
    min_occ = pred->params[2+j];
    max_occ = pred->params[2+m+j];
    value = pred->params[0]+j;
    //BuildObject *xval = CSP::_Variable(value, value);
    if( min_occ > 0 || max_occ < n ) {     
      BuildObject **B = new BuildObject*[n+1];
      for(i=0; i<n; ++i) {
// 	BuildObject **scp = new BuildObject*[3];
// 	scp[0] = x[i];
// 	scp[1] = xval;
// 	int *p = new int[1];
// 	p[0] = 1;
// 	B[i] =  new BuildObjectPredicate( scp, 2, -NOVAL/2, NOVAL/2, ENVIRONMENT[ConstraintStore::EQUAL], p );
	B[i] = CSP::_Equal( x[i], value, 1 );
      }
//       int *os = new int[n+1];
//       std::fill(os, os+n, 1);
//       os[n] = 0;
//       BuildObject *sum = new BuildObjectPredicate( B, n, -NOVAL/2, NOVAL/2, ENVIRONMENT[ConstraintStore::SUM], os );
      BuildObject *sum = CSP::_Sum( B, n );
      if( min_occ > 0 ) {
// 	int *d = new int[1];
// 	d[0] = 0;
// 	BuildObject **scp = new BuildObject*[3];
// 	scp[1] = sum;
// 	scp[0] = CSP::_Variable(min_occ, min_occ);
// 	mod->add( new BuildObjectPredicate( scp, 2, -NOVAL/2, NOVAL/2, ENVIRONMENT[ConstraintStore::PRECEDENCE], d ) );
	mod->add( CSP::_Precedence( min_occ, sum ) );
      }
      if( max_occ < n ) {
// 	int *d = new int[1];
// 	d[0] = 0;
// 	BuildObject **scp = new BuildObject*[3];
// 	scp[0] = sum;
// 	scp[1] = CSP::_Variable(min_occ, min_occ);
// 	mod->add( new BuildObjectPredicate( scp, 2, -NOVAL/2, NOVAL/2, ENVIRONMENT[ConstraintStore::PRECEDENCE], d ) );
	mod->add( CSP::_Precedence( sum, max_occ ) );
      }
    }
  }
}

void BuildObjectGcc::close( BuildObjectPredicate *pred ) 
{  
  BuildObject **x = pred->scope;
  int j, n=pred->arity, i=n, m=(pred->params[1]-pred->params[0]+1);  
  while( i-- )
    if( x[i]->isConstant() ) {
      j = x[i]->value()-pred->params[0];
      --pred->params[2+j];
      --pred->params[2+m+j];
      x[i]->deReference();
      x[i] = x[--n];
    }
  pred->arity = n;
  if(!n) pred->unsetRel();
  else {
    pred->model->unsetSat();
    pred->deReference();
    for(i=0; i<n; ++i) {
      x[i]->unsetBList();
      x[i]->unsetIList();
    }
    pred->unsetBList();
    pred->unsetIList();
  }
}


void BuildObjectGcc::print(std::ostream& o, const BuildObjectPredicate *pred) const
{
  o << "Gcc(";
  for(int i=0; i<pred->arity-1; ++i) {
    pred->scope[i]->print( o );
    o << " " ;
  }
  pred->scope[pred->arity-1]->print( o );
  o << ", [" << pred->params[0] << ".." 
    << pred->params[1] << "] >={";
  for(int i=0; i<(pred->params[1]-pred->params[0]+1); ++i) {
    o << " " << pred->params[2+i];
  }
  o << " }, <={";
  for(int i=0; i<(pred->params[1]-pred->params[0]+1); ++i) {
    o << " " << pred->params[3+pred->params[1]-pred->params[0]+i];
  }
  o << " })";
}

string BuildObjectGcc::toString(const BuildObjectPredicate* pred) const 
{
  //int * w = pred->params+1;
  string s = "Gcc(";
  for(int i=0; i<pred->arity-1; ++i) {
    s += pred->scope[i]->toString();
    s += ", " ;
  }
  s += pred->scope[pred->arity-1]->toString();
  s += ")";
  return s;
}

string BuildObjectGcc::xmlPred(int& idx, int level, const BuildObjectPredicate *pred, int *sub) const 
{
  string s = "global:gcc";
  return s;
}


/**********************************************
 * Lex Constraint BuildObject
 **********************************************/
void BuildObjectLex::add( CSP *m, const int l, 
			  BuildObjectPredicate* pred ) 
{
  assert( l == 0 );
  //BuildObjectConstraint::add(m, l, pred);

  int i, n=(pred->arity >> 1);
  BuildObject **x = pred->scope;
  BuildObject *B[n+1];
  
  B[0] = CSP::_Variable(0,0);
  for(i=1; i<n; ++i) 
    B[i] = CSP::_Variable(0,1);
  B[n] = CSP::_Variable(pred->params[0],1);

  for(i=0; i<n; ++i) {
    BuildObject **scp = new BuildObject*[5];
    scp[0] = x[i];
    scp[1] = x[n+i];
    scp[2] = B[i];
    scp[3] = B[i+1];
    m->add( new BuildObjectPredicate( scp, 4, -NOVAL/2, NOVAL/2, ENVIRONMENT[ConstraintStore::DLEX], NULL ) );
  }

  //pred->deReference();
}

void BuildObjectLex::build( Solver *s, VariableInt **tmp, 
			    BuildObjectPredicate *pred ) 
{
}

void BuildObjectLex::print(std::ostream& o, const BuildObjectPredicate *pred) const
{
  int i;
  o << "[ ";
  for(i=0; i<(pred->arity/2-1); ++i) {
    pred->scope[i]->print( o );
    o << " " ;
  }
  pred->scope[i]->print( o );
  o << "]"; 
  if(pred->params[0])
    o << " <= ";
  else
    o << " < ";
  o << "[ ";
  for(i=(pred->arity/2); i<pred->arity-1; ++i) {
    pred->scope[i]->print( o );
    o << " " ;
  }
  pred->scope[i]->print( o );
  o << "]";
}

string BuildObjectLex::toString(const BuildObjectPredicate* pred) const 
{
  //int * w = pred->params+1;
  string s = "Lex(";
  for(int i=0; i<pred->arity-1; ++i) {
    s += pred->scope[i]->toString();
    s += ", " ;
  }
  s += pred->scope[pred->arity-1]->toString();
  s += ")";
  return s;
}

string BuildObjectLex::xmlPred(int& idx, int level, const BuildObjectPredicate *pred, int *sub) const 
{
  string s = "global:lex";
  return s;
}


/**********************************************
 * DLex Constraint BuildObject
 **********************************************/
void BuildObjectDLex::build( Solver *s, VariableInt **tmp, BuildObjectPredicate *pred ) 
{
  new ConstraintLex( s, tmp );
}

int BuildObjectDLex::propagateUpward( BuildObjectPredicate *pred ) const 
{
  pred->releaseEvent(-1);
  return true;
}

int BuildObjectDLex::propagateDownward( BuildObjectPredicate *pred ) const 
{
  bool consistent = true;
  BuildObject **x = pred->scope;

//   int domains_before = (x[0]->size() + x[1]->size() + x[2]->size() + x[3]->size());
  

  consistent = ( x[2]->setMax( x[3]->max() ) && x[3]->setMin( x[2]->min() ) );
  if( consistent && (x[2]->max() < x[3]->min()) )
    consistent = ( x[0]->setMax( x[1]->max()-1 ) && x[1]->setMin( x[0]->min()+1 ) );
  if( consistent && (x[2]->equal(0) && x[3]->equal(0)) )
    consistent = x[0]->setEqual( x[1] );
  //if( consistent && (x[0]->min() >= x[1]->max()) )
    //consistent =  ( x[2]->setMax( x[3]->max()-1 ) && x[3]->setMin( x[2]->min()+1 ) );  
  if( consistent && (x[0]->min() > x[1]->max()) )
    consistent =  ( x[3]->setMin( 1 ) );  


//   //if( domains_before > (x[0]->size() + x[1]->size() + x[2]->size() + x[3]->size()) ) {
//     pred->print( cout );
//     cout << endl;
//     //}


  return consistent;
}

void BuildObjectDLex::close( BuildObjectPredicate *pred ) 
{  
  BuildObject **x = pred->scope;
  bool valid = ( (x[2]->equal(0) && x[3]->equal(0)) || (x[0]->max() < x[1]->min()) );
  pred->deReference();
  if( valid ) {
    pred->unsetRel();
    x[0]->deReference();
    x[1]->deReference();
    x[2]->deReference();
    x[3]->deReference();
  } else {
    pred->model->unsetSat();
    for(int i=0; i<4; ++i) {
      x[i]->unsetBList();
      x[i]->unsetIList();
    }
  }
}

void BuildObjectDLex::print(std::ostream& o, const BuildObjectPredicate *pred) const
{
  o << "(";
  pred->scope[0]->print( o );
  o << " <=lex ";
  pred->scope[1]->print( o );
  o << " / ";
  pred->scope[2]->print( o );
  o << " ";
  pred->scope[3]->print( o );
  o << ")";

}

string BuildObjectDLex::toString(const BuildObjectPredicate* pred) const 
{
  //int * w = pred->params+1;
  string s = "DLex(";
  for(int i=0; i<pred->arity-1; ++i) {
    s += pred->scope[i]->toString();
    s += ", " ;
  }
  s += pred->scope[pred->arity-1]->toString();
  s += ")";
  return s;
}

string BuildObjectDLex::xmlPred(int& idx, int level, const BuildObjectPredicate *pred, int *sub) const 
{
  string s = " ";
  return s;
}


/**********************************************
 * Cardinality Constraint BuildObject
 **********************************************/
void BuildObjectCardinality::build( Solver *s, VariableInt **tmp, 
				    BuildObjectPredicate *pred ) 
{
}

void BuildObjectCardinality::add( CSP *m, const int l, 
				  BuildObjectPredicate* pred ) 
{
  assert( l == 0 );
  //BuildObjectConstraint::add(m, l, pred);

  BuildObject **x = pred->scope;
  int i, n=pred->arity, k = pred->params[0];

  //BuildObject *kval = CSP::_Variable(k, k);
  BuildObject **B = new BuildObject*[n+1];

  for(i=0; i<n; ++i) {
//     BuildObject **scp = new BuildObject*[3];
//     scp[0] = x[i];
//     scp[1] = kval;
//     int *p = new int[1];
//     p[0] = 1;
//     B[i] =  new BuildObjectPredicate( scp, 2, -NOVAL/2, NOVAL/2, ENVIRONMENT[ConstraintStore::EQUAL], p );
    B[i] = CSP::_Equal( x[i], k, 1 );
  }
//   int *os = new int[n+1];
//   std::fill(os, os+n, 1);
//   os[n] = 0;
//   BuildObject *sum = new BuildObjectPredicate( B, n, -NOVAL/2, NOVAL/2, ENVIRONMENT[ConstraintStore::SUM], os );
  BuildObject *sum = CSP::_Sum( B, n );

//   BuildObject **escp = new BuildObject*[3];
//   int *p = new int[1];
//   p[0] = 1;
//   escp[0] = sum;
//   escp[1] = pred;
//   BuildObject *equ = new BuildObjectPredicate( escp, 2, -NOVAL/2, NOVAL/2, ENVIRONMENT[ConstraintStore::EQUAL], p );

  m->add( CSP::_Equal( sum, pred, 1 ) );
}

void BuildObjectCardinality::print(std::ostream& o, const BuildObjectPredicate *pred) const
{
  int i;
  o << "occ(" << pred->params[0] << ") <= " << pred->params[1] << " in [";
  for(i=0; i<pred->arity-1; ++i) {
    pred->scope[i]->print( o );
    o << "," ;
  }
  pred->scope[pred->arity-1]->print( o );
  o << "]";
}

string BuildObjectCardinality::toString(const BuildObjectPredicate* pred) const 
{
  //int * w = pred->params+1;
  string s = "Cardinality(";
  for(int i=0; i<pred->arity-1; ++i) {
    s += pred->scope[i]->toString();
    s += ", " ;
  }
  s += pred->scope[pred->arity-1]->toString();
  s += ")";
  return s;
}

string BuildObjectCardinality::xmlPred(int& idx, int level, const BuildObjectPredicate *pred, int *sub) const 
{
  string s = "global:maxOccurence";
  return s;
}


/**********************************************
 *  Extensional Constraint BuildObject
 **********************************************/ 
BuildObjectTable::BuildObjectTable(const int a) 
{ 
  arity = a; 
  spin  = 0;
  //table = NULL;
  built_count = 0;
}
BuildObjectTable::BuildObjectTable(const int a, const int nbt) 
{ 
  arity = a; 
  tuples.init(0, nbt);
  spin=0;
  //table = NULL;
  built_count = 0;
}
void BuildObjectTable::addInstance( BuildObjectPredicate* pred )
{
  // params[0] == type (of the leader) 
  // params[1] == index in instance
  // params[2] == next predicate in the same group
  // params[3] == index in leader

  pred->params[0] = -1;
  int idx = instance.size(), settled=false, i=leader.size(), j, included;
  BuildObjectPredicate *L_i;

  pred->params[1] = idx;
  instance.push_back( pred );
  tables.push_back( NULL );
  ++built_count;
  
  assert( instance.size() == (unsigned int)built_count );

  while( !settled && i-- )
    {
      L_i = instance[leader[i]];

      //cout << "compare ";
      //pred->print( cout );
      //cout << " with ";
      //L_i->print( cout );
      //cout << endl;
      
      if(pred->arity == L_i->arity) {
	included = true;
	// L_i is a vaid leader for pred
	for(j=0; included && j<pred->arity; ++j)
	  included = pred->scope[j]->included( L_i->scope[j] );
	
	if( included ) {
	  pred->params[3] = L_i->params[3];
	  pred->params[2] = L_i->params[2];
	  L_i->params[2] = idx;
	  settled = true;

// 	  pred->print( cout );
// 	  cout << " <= ";
// 	  L_i->print( cout );
// 	  cout << endl;

	  // pred is a new leader for L_i's group
	} else {

	  included = true;
	  for(j=0; included && j<pred->arity; ++j)
	    included = L_i->scope[j]->included( pred->scope[j] );
	  if( included ) {
	    selectType( pred );
	    pred->params[3] = L_i->params[3];
	    leader[L_i->params[3]] = idx;
	    pred->params[2] = L_i->params[1];	
	    settled = true;

// 	    L_i->print( cout );
// 	    cout << " <= ";
// 	    pred->print( cout );
// 	    cout << endl;

	  }
	}
      }
    }

  if( !settled ) {

    //cout << "new leader: " << idx << endl;
    selectType( pred );
    pred->params[2] = -1;
    pred->params[3] = leader.size();
    leader.push_back( idx );
  }
}
void BuildObjectTable::printInstances()
{
  unsigned int i, j;
  cout << "constraint "<< arity << " " << tuples.size << endl;
  for(i=0; i<leader.size(); ++i) {
    cout << "Group_" << (i+1) << ": ";
    for( j=leader[i]; j>=0; j=instance[j]->params[2] )
      instance[j]->print( cout );
    cout << endl;
  }    
}
int BuildObjectTable::propagateUpward( BuildObjectPredicate *pred ) const
{
  int consistent=true, n=pred->arity, m=tuples.size;
  if( n == 1 ) {
    if( spin ) {
      int vals[m];
      for(int i=0; i<m; ++i)
	vals[i] = tuples[i][0];
      consistent = pred->scope[0]->setDomain( vals, m );
    } else {
      for(int i=0; consistent && i<m; ++i)	
	consistent = pred->scope[0]->remove( tuples[i][0] );
    }
  } else if( n == 2 ) {
    BuildObject **x = pred->scope;
    int j=-1, k;
    while( ++j<2 && !x[j]->isConstant() );
    if( j<2 ) {
      k = x[j]->value();
      if( spin ) {
	BitSet supports( x[1-j]->min(), x[1-j]->max(), BitSet::empt );
	for(int i=0; consistent && i<m; ++i)
	  if( tuples[i][j] == k )
	    supports.insert( tuples[i][1-j] );
	consistent = x[1-j]->setDomain( supports );
      } else {
	for(int i=0; consistent && i<m; ++i)	
	  consistent = (tuples[i][j] != k || pred->scope[1-j]->remove(tuples[i][1-j]));
      }
    }
  }
  return consistent;
}
 
void BuildObjectTable::close( BuildObjectPredicate *pred ) 
{
  BuildObject **x = pred->scope;
  int j, n=pred->arity, i;
  
  for(i=0; i<n; ++i)
    x[i]->unsetRange();
  
  i=tuples.size;
  if( n == 1 || (!spin && !tuples.size) ) {
    pred->unsetRel();
    x[0]->deReference();
  } else if( n == 2 ) {
    j = -1;
    while( ++j<2 && !x[j]->isConstant() );
    if( j<2 ) {
      pred->unsetRel();
      x[0]->deReference();
      x[1]->deReference();
    } else {
      pred->deReference();
      setNotEqual( pred );
    }
  } else pred->deReference();  


  // if( (spin || tuples.size != 1) && pred->isRel() )
  if( pred->isRel() ) 
    {
 //      //selectType( pred );
//       cout << endl << "add instance: ";
//       pred->print( cout );
//       cout << "  ";
//       //cout << endl;
//       //printInstances();
      
      addInstance( pred );
//       cout << built_count << endl;
//       //printInstances();
    }
}

void BuildObjectTable::build( SatSolver *s, const int *idx, BuildObjectPredicate *pred )
{

  // cout << "SAT" << endl;

  Vector<Literal> new_clause;
  int i, j, k, n=pred->arity, m=tuples.size;

  if( !spin ) {
    for(i=0; i<m; ++i) {
      for(j=0; j<n; ++j) {
	if(tuples[i][j] == 0) k=1;
	else if(tuples[i][j] == 1) k=-1;
	else {
	  cerr << "NOT IMPLEMENTED!" << endl;
	  exit(1);
	}
	new_clause.push( k * idx[pred->scope[j]->getVarId()] );
      }
      s->addClause( s->base, new_clause, s->stats.base_avg_size );
      s->addOriginalClause( new_clause );
      new_clause.clear();
    } 
  } else {
    BitSet matrix(0, ((1 << n)-1), BitSet::full);
    for(i=0; i<m; ++i) {
      k = 0;
      for(j=0; j<n; ++j)
	if( tuples[i][j] == 1 )
	  k += (1 << j);
      matrix.erase( k );
    }
    BitsetIterator bit( matrix );
    do {
      k = bit;
      for(j=0; j<n; ++j)
	if( k & (1 << j) )
	  new_clause.push( -idx[pred->scope[j]->getVarId()] );
	else
	  new_clause.push( (idx[pred->scope[j]->getVarId()]) );
      s->addClause( s->base, new_clause, s->stats.base_avg_size );
      s->addOriginalClause( new_clause );
      new_clause.clear();
    } while( ++bit );
  }
}


void BuildObjectTable::build( Solver *s, VariableInt **tmp, 
			      BuildObjectPredicate *pred ) 
{    
  CSP *model = pred->model;
  if( model->satCompatible ) {
    if( !model->clauseBase )
      model->clauseBase = new ConstraintClauseBase(s);

    Vector<Literal> new_clause;
    int i, j, k, n=pred->arity, m=tuples.size;

    if( !spin ) {
      for(i=0; i<m; ++i) {
	for(j=0; j<n; ++j) {
	  if(tuples[i][j] == 0) k=1;
	  else if(tuples[i][j] == 1) k=-1;
	  else {
	    cout << endl << "s NOT SUPPORTED!" << endl;
	    exit(1);
	  }
	  new_clause.push( k * (tmp[j]->id+1) );
	}
	model->clauseBase->addClause( new_clause );
	new_clause.clear();
      } 
    } else {
      BitSet matrix(0, ((1 << n)-1), BitSet::full);
      for(i=0; i<m; ++i) {
	k = 0;
	for(j=0; j<n; ++j)
	  if( tuples[i][j] == 1 )
	    k += (1 << j);
	matrix.erase( k );
      }      
      BitsetIterator bit( matrix );
      do {
	k = bit;
	for(j=0; j<n; ++j)
	  if( k & (1 << j) )
	    new_clause.push( -(tmp[j]->id+1) );
	  else
	    new_clause.push( tmp[j]->id+1 );
	model->clauseBase->addClause( new_clause );
	new_clause.clear();
      } while( ++bit );
    }    
  } else {

    unsigned int i;
    int j, n=pred->arity, type=instance[leader[pred->params[3]]]->params[0];
//     if( !spin && tuples.size == 1 )
//       {
	
//       }
//     else
//       {    
	switch(type) {
	case 0: {
	  tables[pred->params[1]] = new ConstraintClause( s, tmp, n, tuples[0] );
	  //assert( false );
	} break;
	case 1: {
	  tables[pred->params[1]] = new ConstraintAC3Bitset( s, tmp, spin );
	} break;      
	case 2: {
	  tables[pred->params[1]] = new ConstraintGAC3Valid( s, tmp, n );
	} break;
	default: {
	  tables[pred->params[1]] = new ConstraintGAC2001Allowed( s, tmp, n );
	}
	}
	//   }
    
    if( !--built_count ) 
      {
	// if the tables are not yet initialised, do it!
	for(i=0; i<leader.size(); ++i) {
	  switch( instance[leader[i]]->params[0] ) {
	  case 0: {
	    //assert( false );
	  } break;
	  case 1: {
	    ConstraintAC3Bitset* tab = (ConstraintAC3Bitset*)(tables[leader[i]]);
	    for( j=leader[i]; j>=0; j=instance[j]->params[2] )
	      ((ConstraintAC3Bitset*)(tables[j]))->init( tuples, spin, tab );
	  } break;      
	  case 2: {
	    ConstraintGAC3Valid* tab = (ConstraintGAC3Valid*)(tables[leader[i]]);
	    for( j=leader[i]; j>=0; j=instance[j]->params[2] ) {
	      ((ConstraintGAC3Valid*)(tables[j]))->init( tuples, spin, tab );
	      ((ConstraintGAC3Valid*)(tables[j]))->useResidualSupports();
	    }
	  } break;
	  default: {
	    ConstraintGAC2001Allowed* tab = (ConstraintGAC2001Allowed*)(tables[leader[i]]);
	    for( j=leader[i]; j>=0; j=instance[j]->params[2] )
	      ((ConstraintGAC2001Allowed*)(tables[j]))->init( tuples, spin, tab );
	  }
	  }
	}    
      }
  } 
}


void BuildObjectTable::selectType( BuildObjectPredicate *pred )
{
  BuildObject **x = pred->scope;
  int n=pred->arity, i;
  long int matrixsize = 1;
  double rtightness;
  for(i=0; matrixsize > 0 && i<n; ++i)
    matrixsize *= (pred->scope[i]->max() - pred->scope[i]->min() + 1);
  if( matrixsize < 0 ) {
    rtightness = 1;
  } else {
    if( spin ) {
      rtightness = 1-((double)(tuples.size) / (double)(matrixsize));
    } else {
      rtightness = ((double)(tuples.size) / (double)(matrixsize));
    }
  }

  // cout << spin << " " << tuples.size << " " << rtightness << " " << matrixsize << endl;

  
  if( spin || tuples.size ) {
    if( !spin && tuples.size == 1 ) { // CLAUSE 
      //assert( false );
      pred->params[0] = 0;
    } else if( rtightness >= .001 && n == 2 ) { // SUPPORT BITSETS
      pred->params[0] = 1;
      for(i=0; i<n; ++i)
	x[i]->unsetIList();
    } else if( !spin || (matrixsize < 32000000 && rtightness < 0.9) ) { // N-ARY MATRIX
      pred->params[0] = 2;
      for(i=0; i<n; ++i)
	x[i]->unsetIList();
    } else { // LIST OF TUPLES
      pred->params[0] = 3;
    }
  }
  if( spin && matrixsize > 512 )
    pred->model->unsetSat();   
}

void BuildObjectTable::setNotEqual( BuildObjectPredicate *pred )
{
  BuildObject **x = pred->scope;
  int k, j, i=tuples.size;
  if( spin ) {
    while( i-- )
      if(tuples[i][0] == tuples[i][1]) break;
    if( i < 0 )
      pred->setNeq();
  } else {
    int borneinf = std::max( x[0]->min(), x[1]->min() );
    int bornesup = std::min( x[0]->max(), x[1]->max() );
    k=(bornesup - borneinf + 1);
    int values[k];
    int *vals = values-borneinf;
    for(j=borneinf; j<=bornesup; ++j)
      if( x[0]->contain(j) && x[1]->contain(j) )
	vals[j] = 1;
      else {
	vals[j] = 0;
	--k;
      }
    while( k && i-- ) {
      j = tuples[i][0];
      bool flag1 = (j == tuples[i][1]);
      bool flag2 = (j >= borneinf);
      bool flag3 = (j <= bornesup);
      bool flag4 = (vals[j]);

      if( flag1 && flag2 && flag3 && flag4 )
      //if(j == tuples[i][1] && j >= borneinf && j <= bornesup && vals[j]) 
	{
	  --k;
	  vals[j] = 0;
	}
    }
    if( !k )
      pred->setNeq();
  }
}

BuildObjectTable::~BuildObjectTable() 
{
  int i;
  i=tuples.size;
  while( i-- ) {
    delete [] tuples[i];
  }
}

void BuildObjectTable::print(std::ostream& o, const BuildObjectPredicate *pred) const
{
  o << "table(";
  for(int i=0; i<pred->arity-1; ++i) {
    pred->scope[i]->print( o );
    o << " " ;
  }
  pred->scope[pred->arity-1]->print( o );
  o << ") " << tuples.size;
}

string BuildObjectTable::toString(const BuildObjectPredicate* pred) const 
{
  //int * w = pred->params+1;
  string s = "Table(";
  for(int i=0; i<pred->arity-1; ++i) {
    s += pred->scope[i]->toString();
    s += ", " ;
  }
  s += pred->scope[pred->arity-1]->toString();
  s += ")";
  return s;
}

string BuildObjectTable::xmlPred(int& idx, int level, const BuildObjectPredicate *pred, int *sub) const 
{
  string s = "table";
  return s;
}

void BuildObjectTable::add( const int* sol )
{
  int i = arity, *ngd = new int[arity];
  while( i-- ) 
    ngd[i] = sol[i];
  tuples.push_back( ngd );
}







/**********************************************
 * Adding variables/constraints to the model
 **********************************************/

int CSP::numVars() const
{
  return (nList + nBit + nBoolean + nRange) ;
}

int CSP::numEVars() const
{
  return (eList + eBit + eBoolean + eRange) ;
}

inline void CSP::triggerEvent(int id, int e) 
{
  buildqueue.push( id, e );
}


void CSP::add( BuildObject *x ) 
{
  //  #ifdef _DEBUGMODEL
//     cout << "add ";
//     x->print( cout );
//     cout << endl;
//  #endif

// // //   if( x->id >= 98 ) {
//      cout << "add ";
//      x->print( cout );
//      cout << endl;
// // //   }

     
     if(x->isRel()) toplevel_expressions.push(x);
     x->add( this, 0 );

// #ifdef _DEBUGMODEL
//   cout << " : " << (x->id) << endl;
  
//   print( cout );
//   cout << endl << endl;
// #endif
// //   if( x->id >= 98 ) {

// //      cout << declarations[98]->isReferenced() << endl;
// //      declarations[98]->print( cout );
// //      cout << endl;
    
// //     exit( 0 );

// //   }
}

void CSP::add( BuildObjectObjective *x ) 
{
  x->X->add( this, 1 );
  goald = x;
}

void CSP::add( BuildObject **x, const int n ) 
{
  for(int i=0; i<n; ++i) 
    x[i]->add( this, 0 );
}

void CSP::add( Variable x )
{
  x.var_ptr_->add( this, 0 );
}

void CSP::add( VarArray& x )
{
  int n=x.size();
  for(int i=0; i<n; ++i) {
    add( x[i].var_ptr_ );
    //x[i].var_ptr_->add( this, 0 );
  }
}

void CSP::add( Objective g )
{
  g.X->add( this, 1 );
  goald = g.obj;
}

// void CSP::close()
// {
//   //buildqueue.initStack( declarations.size ) ;
// }

int CSP::inferAllDiffs( const int all, const int limit )
{
  int i, n=declarations.size, // m=notequals.size,
    x, y, ncliques=0, maxsize;

  // compute nodes and edges of the graph
  vector< int > edges;
  vector< int > vertice;
  int vertice_index[n];
  std::fill(vertice_index, vertice_index+n, -1);
  
  i = notequals.size;

  while( i-- )
    {
      x = ((BuildObjectPredicate*)notequals[i])->scope[0]->getVarId();
      y = ((BuildObjectPredicate*)notequals[i])->scope[1]->getVarId();
      
      if( vertice_index[x] < 0  ) {
	vertice_index[x] = vertice.size();
	vertice.push_back( x );
      }
      if( vertice_index[y] < 0  ) {
	vertice_index[y] = vertice.size();
	vertice.push_back( y );
      }
      
      edges.push_back( vertice_index[x] );
      edges.push_back( vertice_index[y] );
    }
    
  if( vertice.size() < 4096 && edges.size() > (vertice.size() * 2) )
    {
      MistralGraph G( vertice.size(), edges );      
      vector< vector< int > > cliques;      
      if( all ) G.getAllCliques( cliques, limit );
      else G.getSomeCliques( cliques );
      
      x = cliques.size();
      while( x-- ) {
	i = n = cliques[x].size();
	maxsize = 0;
// 	while( i-- ) {
// 	  y = declarations[vertice[cliques[x][i]]]->size();
// 	  if(maxsize < y) maxsize = y;
// 	}
	
// 	cout << '\t' << n << " " << maxsize << endl;

	if( n > 3 && maxsize <= n+1 ) {
	  i=n;
	  ++ncliques;
	  BuildObject **tmp = new BuildObject*[n+1];
	  while( i-- ) tmp[n-i-1] = declarations[vertice[cliques[x][i]]];
	  int *par = new int[1];
	  par[0] = 2;
	  BuildObjectPredicate *cliqueAllDiff = 
	    new BuildObjectPredicate( tmp, n, -NOVAL/2, NOVAL/2, 
				      ENVIRONMENT[ConstraintStore::ALLDIFF],
				      par ) ;
	  cliqueAllDiff->add( this, 0 );
	  cliqueAllDiff->close();
	}
      }
    }

  return ncliques;
}
  
bool CSP::preBuild()
{
  if(closed) return true;

  //close();
  BuildObject* x;
  int i, n=declarations.size, event;
  bool consistent = true;

  for(i=0; i<n; ++i) {
    declarations[i]->varptr_ = NULL;
  }

#ifdef _DEBUGMODEL
  cout << "pre-build" << endl;
  print(cout);
  cout << endl;
#endif


  buildqueue.initList( n, declarations.stack_ ) ;
  for(i=0; i<n; ++i) {
    if( declarations[i]->isRel() )  {
	consistent = declarations[i]->propagateUpward( );
	if( !consistent ) {
	  //cout << "c found inconsistent while building" << endl;
	  return false;
	} 

	//      declarations[i]->propagateUpward();
    } 
    //if( !declarations[i]->isRel() ) 
    //if( declarations[i]->isReferenced() ) 
    //declarations[i]->releaseEvent( -1 );
  }


//   for(i=0; i<n; ++i) 
//     if( declarations[i]->isConstant() ) 
//       triggerEvent( i, 1 );      
  

  while( !buildqueue.empty() ) {

    x = buildqueue.pop( event );
   
#ifdef _DEBUGMODEL
    cout << endl << "event on ";     
    x->print( cout );
    cout << endl;
#endif    

    i=x->parent.size;

    while( i-- )
      if( x->parent[i]->isRel() ) {

#ifdef _DEBUGMODEL	
	cout << "\tpropagate up: ";
	x->parent[i]->print( cout );
	cout << endl ;
#endif
	
	consistent = x->parent[i]->propagateUpward( );
	if( !consistent ) {
	  //cout << "c found inconsistent while building" << endl;
	  return false;
	} 

#ifdef _DEBUGMODEL
	cout << "\tresult:       ";
	x->parent[i]->print( cout );       
	cout << endl;
#endif

      }
    if( x->isRel() ) {

#ifdef _DEBUGMODEL
      cout << "\tpropagate down: ";
      x->print( cout );
      cout << endl ;
#endif

      consistent = x->propagateDownward( );
      if( !consistent ) {
	//cout << "c found inconsistent while building" << endl;
	return false;
      } 
      

#ifdef _DEBUGMODEL
      cout << "\tresult:         ";
      x->print( cout );
      cout << endl;
#endif

    }    
  }

#ifdef _DEBUGMODEL
  cout << endl;
  print(cout);
  cout << endl;
  cout << "END PROPAG" << endl;
#endif


  //for(i=0; i<n; ++i)
  i=n;
  while( i-- )
    //if( declarations[i]->isRel() ) 
    {

#ifdef _DEBUGMODEL
      cout << "\tclose:  _" << (declarations[i]->isReferenced()) << "_ ";
      declarations[i]->print( cout );
      cout << endl ;
#endif

      declarations[i]->close();

#ifdef _DEBUGMODEL
      cout << "\tresult: _" << (declarations[i]->isReferenced()) << "_ ";
      declarations[i]->print( cout );
      cout << endl;
#endif

    }

#ifdef _DEBUGMODEL
  cout << endl;
  print(cout);
  cout << endl;
  cout << "END" << endl;
#endif

  close();

  return true;
}

void CSP::close()
{
  if(closed) return;

  int i=declarations.size;
  while( i-- ) {

    if( declarations[i]->isNeq() ) {
      ++nNeq;
      notequals.push( declarations[i] );
    }

    if( declarations[i]->isReferenced() ) {
      if( declarations[i]->isSearchable() ) {
	nValues += declarations[i]->size();
	if( declarations[i]->isBoolean() ) {
	  if( declarations[i]->isRel() )
	    ++eBoolean;
	  else
	    ++nBoolean;
	  //cout << "Boolean " ;
	} else if( declarations[i]->isInterval() ) {
	  if( declarations[i]->isRel() )
	    ++eRange;
	  else
	    ++nRange;
	  //cout << "Range " ;
	} else if( declarations[i]->isBList() || declarations[i]->isIList() ) {
	  if( declarations[i]->isRel() )
	    ++eList;
	  else
	    ++nList;
	  //cout << "Domain " ;
	} else {
	  if( declarations[i]->isRel() )
	    ++eBit;
	  else
	    ++nBit;
	  //cout << "Domain " ;
	}
      } else {
	//cout << "Constant " ;
	if( declarations[i]->isRel() )
	  ++eConstant;
	else
	  ++nConstant;
      }
    }
    //if( declarations[i]->isConstraint() ) {
    if( declarations[i]->isRel() ) {

      //cout << "Constraint " ;
      if( declarations[i]->isReferenced() )
	++eConstraint;
      else
	++nConstraint;
    }
    //    cout << endl;
  }

  closed = true;
}


void CSP::build( Solver *s )
{

//   cout << nBoolean << endl
//        << eBoolean << endl
//        << nRange << endl
//        << eRange << endl
//        << nList << endl
//        << eList << endl
//        << nBit << endl
//        << eBit << endl
//        << nConstant << endl
//        << eConstant << endl;
  
//   exit(0);


  //cout << "BUILD" << endl;

  //s->initDomainPool( nBoolean+eBoolean );

  int i, n=declarations.size;

  for(i=0; i<n; ++i) {
#ifdef _DEBUGMODEL
    cout << "\tbuild:   " ;
    declarations[i]->print( cout );
    cout << endl ;
#endif

    declarations[i]->build( s );
    
#ifdef _DEBUGMODEL
    cout << "\tresult: " ;
    declarations[i]->print( cout );
    cout << endl;
#endif

  }


  if( clauseBase )
    clauseBase->initialise();

  if( goald )
    goald->build( s );
}


void CSP::build( SatSolver *s ) 
{

  int i, n=declarations.size, k=0;

  // get variables
  int index[n];
  for(i=0; i<n; ++i) 
    if( declarations[i]->isVar() )
      index[i] = ++k;

  s->init(k, 32768);

  for(i=0; i<n; ++i) {

#ifdef _DEBUGMODEL
    cout << "\tbuild:   " ;
    declarations[i]->print( cout );
    cout << endl ;
#endif

    declarations[i]->build( s, index );
    
#ifdef _DEBUGMODEL
    cout << "\tresult: " ;
    declarations[i]->print( cout );
    cout << endl;
#endif
  }

}



/**********************************************
 * Operator overloading
 **********************************************/

Variable VarArray::operator[]( Variable x ) 
{
  return Element( *this, x );
}

Variable VarArray::operator< ( VarArray& x ) 
{
  return LexOrder( *this, x, 1 );
}
Variable VarArray::operator> ( VarArray& x ) 
{
  return LexOrder( x, *this, 1 );
}
Variable VarArray::operator<=( VarArray& x ) 
{
  return LexOrder( *this, x, 0 );
}
Variable VarArray::operator>=( VarArray& x )
{
  return LexOrder( x, *this, 0 );
}

Variable Variable::operator==( const BitSet& s ) 
{
  return Member( *this, s, 1 );
}

Variable Variable::operator!=( const BitSet& s ) 
{
  return Member( *this, s, 0 );
}

Variable Variable::operator==( const int k ) 
{
  return Equal( *this, k, 1 );
}

Variable Variable::operator%( const int k ) 
{
  return Mod( *this, k );
}

Variable Variable::operator!=( const int k ) 
{
  return Equal( *this, k, 0 );
}

Variable Variable::operator<=( const int k ) 
{
  return Precedence( *this, k );
}

Variable Variable::operator< ( const int k ) 
{
  return Precedence( *this, k-1 );
}

Variable Variable::operator>=( const int k ) 
{
  return Precedence( k, *this );
}

Variable Variable::operator> ( const int k ) 
{
  return Precedence( k+1, *this );
}

Variable Variable::operator+ ( Variable x )
{
  return Sum( *this, x, 1 );
}

Variable Variable::operator- ( Variable x )
{
  return Sum( *this, x, -1 );
}

Variable Variable::operator* ( Variable x )
{
  return Mul( *this, x );
}

Variable Variable::operator/ ( Variable x )
{
  return Div( *this, x );
}

Variable Variable::operator^ ( Variable x )
{
  cerr << "The operator \"^\" is not implemented" << endl;
  exit(0);
}

Variable Variable::operator%( Variable x )
{
  return Mod( *this, x );
}

Variable Variable::operator&&( Variable x )
{
  return And( *this, x );
}

Variable Variable::operator||( Variable x )
{
  return Or( *this, x );
}

Variable Variable::operator==( Variable x ) 
{
  return Equal( *this, x, 1 );
}
 
Variable Variable::operator!=( Variable x ) 
{
  return Equal( *this, x, 0 );
}

Variable Variable::operator<=( Variable x ) 
{
  return Precedence( *this, 0, x );
}

Variable Variable::operator>=( Variable x ) 
{
  return Precedence( x, 0, *this );
}

Variable Variable::operator< ( Variable x ) 
{
  return Precedence( *this, 1, x );
}

Variable Variable::operator> ( Variable x ) 
{
  return Precedence( x, 1, *this );
}

Variable Variable::operator-( )
{
  return Negation( *this );
}

Variable Variable::operator!( )
{
  return Not( *this );
}

Variable Variable::operator*( const int k )
{
  if( k != 1 ) 
    return Mul( *this, k );
  else return *this;
}

Variable Variable::operator-( const int k )
{
  if( k ) 
    return Sum( *this, -k );
    //return 
  else return *this;
}

Variable Variable::operator+( const int k )
{
  if( k ) 
    return Sum( *this, k );
  else return *this;
}


/**********************************************
 *  Negation Predicate Wrapper
 **********************************************/ 
Negation::Negation( Variable& x_ )
{
  var_ptr_ = CSP::_Negation( x_.var_ptr_ );
}

/**********************************************
 *  Not Predicate Wrapper
 **********************************************/ 
Not::Not( Variable& x_ )
{
  var_ptr_ = CSP::_Not( x_.var_ptr_ );
}

/**********************************************
 *  Abs Predicate Wrapper
 **********************************************/ 
Abs::Abs( Variable x_ )
{
  var_ptr_ = CSP::_Abs( x_.var_ptr_ );
}

/**********************************************
 *  And Predicate Wrapper
 **********************************************/ 
And::And( Variable& x_, Variable& y_ )
{
  var_ptr_ = CSP::_And( x_.var_ptr_, y_.var_ptr_ );
}

/**********************************************
 *  Or Predicate Wrapper
 **********************************************/ 
Or::Or( Variable& x_, Variable& y_ )
{
  var_ptr_ = CSP::_Or( x_.var_ptr_, y_.var_ptr_ );
}

/**********************************************
 *  IfThenElse Predicate Wrapper
 **********************************************/ 
IfThenElse::IfThenElse( Variable& x_, Variable& y_, Variable& z_ )
{
  var_ptr_ = CSP::_IfThenElse( x_.var_ptr_, y_.var_ptr_, z_.var_ptr_ );
}

/**********************************************
 *  Disjunctive Predicate Wrapper
 **********************************************/ 
Disjunctive::Disjunctive( Variable& x_, int dx, Variable& y_, int dy, const int t )
{
  var_ptr_ = CSP::_Disjunctive( x_.var_ptr_, dx, y_.var_ptr_, dy, t );
}

/**********************************************
 *  Overlap Predicate Wrapper
 **********************************************/ 
Overlap::Overlap( Variable& x_, int dx, Variable& y_, int dy )
{
  var_ptr_ = CSP::_Overlap( x_.var_ptr_, dx, y_.var_ptr_, dy );
}

// Disjunctive::Disjunct( Variable& x_, int dx, Variable& y_, int dy )
// {
//   var_ptr_ = CSP::_Disjunctive( x_.var_ptr_, dx, y_.var_ptr_, dy, 1 );
// }

/**********************************************
 *  Precedence Predicate Wrapper
 **********************************************/ 
Precedence::Precedence( Variable& x_, const int d, Variable& y_, const int tp )
{
  var_ptr_ = CSP::_Precedence( x_.var_ptr_, d, y_.var_ptr_ );
}

Precedence::Precedence( Variable& x_, const int v )
{
  var_ptr_ = CSP::_Precedence( x_.var_ptr_, v );
}

Precedence::Precedence( const int v, Variable& x_ )
{
  var_ptr_ = CSP::_Precedence( v, x_.var_ptr_ );
}

/**********************************************
 *  Mul Predicate Wrapper
 **********************************************/
Mul::Mul( Variable& x_, Variable& y_ )
{
  var_ptr_ = CSP::_Mul( x_.var_ptr_, y_.var_ptr_ );
}

Mul::Mul( Variable& x_, const int v )
{
  var_ptr_ = CSP::_Mul( x_.var_ptr_, v );
}

/**********************************************
 *  Div Predicate Wrapper
 **********************************************/
Div::Div( Variable& x_, Variable& y_ )
{
  var_ptr_ = CSP::_Div( x_.var_ptr_, y_.var_ptr_ );
}

/**********************************************
 *  Mod Predicate Wrapper
 **********************************************/
Mod::Mod( Variable& x_, Variable& y_ )
{
  var_ptr_ = CSP::_Mod( x_.var_ptr_, y_.var_ptr_ );
}

Mod::Mod( Variable& x_, const int v )
{
  var_ptr_ = CSP::_Mod( x_.var_ptr_, v );
}

/**********************************************
 *  Equal Predicate Wrapper
 **********************************************/ 
Equal::Equal( Variable& x_, Variable& y_, const int spin )
{
  var_ptr_ = CSP::_Equal( x_.var_ptr_, y_.var_ptr_, spin );
}

Equal::Equal( Variable& x_, const int v, const int spin )
{
  var_ptr_ = CSP::_Equal( x_.var_ptr_, v, spin );
}

/**********************************************
 *  Member Predicate Wrapper
 **********************************************/ 
Member::Member( Variable& x_, const BitSet& k, const int spin )
{
  var_ptr_ = CSP::_Member( x_.var_ptr_, k, spin );
}

/**********************************************
 * Sum Predicate Wrapper
 **********************************************/
Sum::Sum( Variable& x_, const int k )
{
  //  if( k )
  var_ptr_ = CSP::_Sum( x_.var_ptr_, k );
  //  else var_ptr_ = x_.var_ptr_;
}
Sum::Sum( VarArray& x_, const int *os )
{
  var_ptr_ = CSP::_Sum( x_.getArgs(), x_.size(), os );
}
Sum::Sum( Variable& x_, Variable& y_, const int t )
{
  var_ptr_ = CSP::_Sum( x_.var_ptr_, y_.var_ptr_, t );
}

/**********************************************
 * Min Predicate Wrapper
 **********************************************/
Min::Min( VarArray& x_ )
{
  var_ptr_ = CSP::_Min( x_.getArgs(), x_.size() );
}

/**********************************************
 * Max Predicate Wrapper
 **********************************************/
Max::Max( VarArray& x_ )
{
  var_ptr_ = CSP::_Max( x_.getArgs(), x_.size() );
}

/**********************************************
 * Element Predicate Wrapper
 **********************************************/
Element::Element( VarArray& x_, Variable& y_ )
{
  int n = x_.size();//, *p = new int[1];
  //p[0] = 0;
  BuildObject **x = new BuildObject*[n+2];
  for(int i=0; i<n; ++i)
    x[i] = x_[i].var_ptr_;
  x[n] = y_.var_ptr_;  
  //var_ptr_ = CSP::_Constraint(ConstraintStore::ELEMENT, x, n+1, p );
  var_ptr_ = CSP::_Element( x, n+1, 0 );
}

Element::Element( VarArray& x_ )
{
  int n = x_.size(), *p = new int[1];
  p[0] = 1;
  BuildObject **x = new BuildObject*[n];
  for(int i=0; i<n-2; ++i)
    x[i] = x_[i+1].var_ptr_;
  x[n-2] = x_[0].var_ptr_; 
  BuildObject *y = x_[n-1].var_ptr_;
  var_ptr_ = CSP::_Constraint(ConstraintStore::ELEMENT, x, n-1, p );
  y->setEqual( var_ptr_ ) ;
  var_ptr_->id = NOVAL;
}

/**********************************************
 * Cardinality < Constraint Wrapper
 **********************************************/
Cardinality::Cardinality( VarArray& x_, const int k )
{
  var_ptr_ = CSP::_Cardinality( x_.getArgs(), x_.size(), k );
}

/**********************************************
 * AllDifferent Constraint Wrapper
 **********************************************/
AllDifferent::AllDifferent( VarArray& x_, const int c, const int k )
{
  var_ptr_ = CSP::_AllDifferent( x_.getArgs(), x_.size(), c, k );
}


/**********************************************
 * Gcc Constraint Wrapper
 **********************************************/
Gcc::Gcc( VarArray& x_,
	  const int firstDomainValue, 
	  const int lastDomainValue,
	  const int* minOccurrences,
	  const int* maxOccurrences )
{
  int n = (lastDomainValue-firstDomainValue+1);
  int *p = new int[2+2*n];
  p[0] = firstDomainValue;
  p[1] = lastDomainValue;  
  for(int i=0; i<n; ++i) {
    p[2+i] = minOccurrences[i];
    p[2+n+i] = maxOccurrences[i];
  }
  var_ptr_ = CSP::_Constraint(ConstraintStore::GCC, x_.getArgs(), x_.size(), p );
}

/**********************************************
 * TransitiveDAG Constraint Wrapper
 **********************************************/
TransitiveDAG::TransitiveDAG( VarArray& x_ )
{
  var_ptr_ = CSP::_Constraint(ConstraintStore::TDAG, x_.getArgs(), x_.size(), NULL);
}

/**********************************************
 * Tree Constraint Wrapper
 **********************************************/
Tree::Tree( VarArray& x_ )
{
  int *p = new int[1];
  p[0] = 0;
  var_ptr_ = CSP::_Constraint(ConstraintStore::TREE, x_.getArgs(), x_.size(), p);
}

Tree::Tree( VarArray& x_, VarArray& y_ )
{
  int n = x_.size(), m = y_.size();
  BuildObject **x = new BuildObject*[n+m+1];
  for(int i=0; i<n; ++i) 
    x[i] = x_[i].var_ptr_;
  for(int i=0; i<m; ++i) 
    x[i] = y_[i].var_ptr_;
  int *p = new int[1];
  p[0] = n;
  var_ptr_ = CSP::_Constraint(ConstraintStore::TREE, x, n+m, p); 
}

/**********************************************
 * Lex < Constraint Wrapper
 **********************************************/
LexOrder::LexOrder( VarArray& x_, VarArray& y_, int eq )
{
  int n = x_.size();
  BuildObject **x = new BuildObject*[2*n+1];
  for(int i=0; i<n; ++i) {
    x[i] = x_[i].var_ptr_;
    x[i+n] = y_[i].var_ptr_;
  }
  int *p = new int[1];
  p[0] = eq;
  //   var_ptr_ = new BuildObjectPredicate( x, 2*n, -NOVAL/2, NOVAL/2, 
  // 				       ENVIRONMENT[ConstraintStore::LEX], p );
  var_ptr_ = CSP::_Constraint(ConstraintStore::LEX, x, 2*n, p );
}


/**********************************************
 *  Extensional Constraint Wrapper
 **********************************************/ 
int CSP::addTable( const int arity, const int nbTuples )
{
  return ENVIRONMENT.addConstraint( new BuildObjectTable(arity, nbTuples) );
}
void CSP::addTuple( const int idx, int *tuple )
{
  ((BuildObjectTable*)(ENVIRONMENT[idx]))->add( tuple );
}
Table::Table() {}
Table::Table( VarArray& x_, const bool sp )
{
  var_ptr_ = CSP::_Table( x_.getArgs(), x_.size(), sp ); 
  tabptr_ = (BuildObjectTable*)(((BuildObjectPredicate*)var_ptr_)->relation);
}
Table::Table( VarArray& x_, const int tidx, const bool sp )
{
  var_ptr_ = CSP::_Table( x_.getArgs(), x_.size(), tidx, sp );
  tabptr_ = (BuildObjectTable*)(((BuildObjectPredicate*)var_ptr_)->relation);
}

void Table::add( const int* sol )
{
  tabptr_->add( sol ) ;
}



/**********************************************
 * Printing routines
 **********************************************/

std::ostream& operator<< (std::ostream& os, const Mistral::Matrix& x) 
{
  x.print( os );
  return os;
}

std::ostream& operator<< (std::ostream& os, const Mistral::Variable& x)
{
  x.print( os );
  return os;
}

std::ostream& operator<< (std::ostream& os, const Mistral::VarArray& x)
{
  //  x.print( os );
  x.print( os );
  return os;
}

std::ostream& operator<< (std::ostream& os, const Mistral::CSP& p)
{
  p.print( os );
  return os;
}



/**********************************************
 * Objective Function BuildObjects
 **********************************************/
BuildObjectObjective::BuildObjectObjective(BuildObject *x) : X(x) {}

void BuildObjectObjective::print( std::ostream& o ) const {}


/// Maximisation declaration
BuildObjectMaximiseVar::BuildObjectMaximiseVar(BuildObject *x) 
  : BuildObjectObjective(x)
{
}

void BuildObjectMaximiseVar::print( std::ostream& o ) const 
{
  o << "Maximise ";
  X->print( o ) ;
}

void BuildObjectMaximiseVar::build( Solver *s )
{
  s->goal = new MaximiseVar( s, X->getVariable() );
}

/// Minimisation declaration
BuildObjectMinimiseVar::BuildObjectMinimiseVar(BuildObject *x) 
  : BuildObjectObjective(x)
{
}

void BuildObjectMinimiseVar::print( std::ostream& o ) const 
{
  o << "Minimise ";
  X->print( o ) ;
}

void BuildObjectMinimiseVar::build( Solver *s )
{
  s->goal = new MinimiseVar( s, X->getVariable() );
}


/**********************************************
 * CSP
 **********************************************/

void CSP::setGAC( const int g ) {P_GAC = g;}

CSP::CSP() 
{
  satCompatible = true;
  clauseBase = NULL;
  closed = false;

  P_GAC = DYNAMIC;
  //P_BOUND = 0;
  //P_DOMAIN = DYNAMIC;

  nValues     = 0;
  nList       = 0;
  nBit        = 0;
  nBoolean    = 0;
  nRange      = 0;
  nConstant   = 0;
  nConstraint = 0;
  eValues     = 0;
  eList       = 0;
  eBit        = 0;
  eBoolean    = 0;
  eRange      = 0;
  eConstant   = 0;
  eConstraint = 0;
  nNeq        = 0;
  goald = NULL;
}

CSP::~CSP() 
{
  
  //std::cout << "c (mistral) delete model" << std::endl;
  //std::cout  << "*** DELETE MODEL ***" << endl;

  ENVIRONMENT.flush();
  
  cleanUp();
  
  delete goald;
  
  declarations.clear();
  goald = NULL;
}

void CSP::cleanUp()
{
  nValues     = 0;
  nList       = 0;
  nBit        = 0;
  nBoolean    = 0;
  nRange      = 0;
  nConstant   = 0;
  nConstraint = 0;
  eValues     = 0;
  eList       = 0;
  eBit        = 0;
  eBoolean    = 0;
  eRange      = 0;
  eConstant   = 0;
  eConstraint = 0;
  int i = declarations.size;
  while( i-- ) {
    delete declarations[i];
  }
}

bool CSP::isSatCompatible() const {
  return satCompatible;
}
void CSP::unsetSat() {
  satCompatible = false;
}


BuildObject* CSP::_Variable( const int nvals ) 
{
  BuildObject *x = new BuildObject( 0, nvals-1 );
  //x->setVar();
  return x;
}

BuildObject* CSP::_Variable( const int lb, const int ub ) 
{
  BuildObject *x = new BuildObject( lb, ub );
  //x->setVar();
  return x;
}

BuildObject* CSP::_Variable( const int* array, const int length ) 
{
  BuildObject *x = new BuildObject(array, length);
  //x->setVar();
  return x;
}

BuildObject* CSP::_Variable( const int* array, const int length, 
			     const int lb, const int ub ) 
{
  BuildObject *x = new BuildObject(array, length, lb, ub);
  //x->setVar();
  return x;
}

BuildObject* CSP::_Variable( const BitSet& d, const int length, 
			     const int lb, const int ub ) 
{
  BuildObject *x = new BuildObject(d, length, lb, ub);
  //x->setVar();
  return x;
}

BuildObject** CSP::_VarArray( const int nvars, const int nvals )
{
  BuildObject **array_ = new BuildObject*[nvars];
  int i=nvars;
  while( i-- )
    array_[i] = CSP::_Variable( nvals );
  return array_;
}

BuildObject** CSP::_VarArray( const int nvars, const int lb, const int ub ) 
{
  BuildObject **array_ = new BuildObject*[nvars];
  int i=nvars;
  while( i-- ) 
    array_[i] = CSP::_Variable( lb, ub );
  return array_;
}

BuildObject** CSP::_VarArray( const int nvars, const int* array, const int length ) 
{
  BuildObject **array_ = new BuildObject*[nvars];
  int i=nvars;
  while( i-- ) 
    array_[i] = CSP::_Variable( array, length );
  return array_;
}

BuildObject** CSP::_VarArray( const int nvars, const int* array, 
			     const int length, const int lb, const int ub )
{
  BuildObject **array_ = new BuildObject*[nvars];
  int i=nvars;
  while( i-- ) 
    array_[i] = CSP::_Variable( array, length, lb, ub );
  return array_;
}


BuildObject* CSP::_Constraint( const int t, BuildObject **x, const int n, int*p )
{
  return new BuildObjectPredicate( x, n, -NOVAL/2, NOVAL/2, ENVIRONMENT[t], p );
}

BuildObject* CSP::_Negation( BuildObject *x_ )
{
//   BuildObject **x = new BuildObject*[2];
//   x[0] = x_.var_ptr_;
//   var_ptr_ = CSP::_Constraint(ConstraintStore::NEG, x, 1, NULL);
  BuildObject **x = new BuildObject*[2];
  x[0] = x_;
  return CSP::_Constraint(ConstraintStore::NEG, x, 1, NULL);
}
BuildObject* CSP::_Not( BuildObject *x_ )
{
//   BuildObject **x = new BuildObject*[2];  
//   x[0] = x_.var_ptr_;
//   var_ptr_ = CSP::_Constraint(ConstraintStore::NOT, x, 1, NULL);
  BuildObject **x = new BuildObject*[2];  
  x[0] = x_;
  return CSP::_Constraint(ConstraintStore::NOT, x, 1, NULL);
}
BuildObject* CSP::_Abs( BuildObject *x_ )
{
//   BuildObject **x = new BuildObject*[2];  
//   x[0] = x_.var_ptr_;
//   var_ptr_ = CSP::_Constraint(ConstraintStore::ABS, x, 1, NULL);
  BuildObject **x = new BuildObject*[2];  
  x[0] = x_;
  return CSP::_Constraint(ConstraintStore::ABS, x, 1, NULL);
}
BuildObject* CSP::_And( BuildObject *x_, BuildObject *y_ )
{
//   BuildObject **x = new BuildObject*[3];  
//   x[0] = x_.var_ptr_;
//   x[1] = y_.var_ptr_;
//   var_ptr_ = CSP::_Constraint(ConstraintStore::AND, x, 2, NULL);
  BuildObject **x = new BuildObject*[3];  
  x[0] = x_;
  x[1] = y_;
  return CSP::_Constraint(ConstraintStore::AND, x, 2, NULL);
}
BuildObject* CSP::_Or( BuildObject *x_, BuildObject *y_ )
{
//   BuildObject **x = new BuildObject*[3];  
//   x[0] = x_.var_ptr_;
//   x[1] = y_.var_ptr_;
//   var_ptr_ = CSP::_Constraint(ConstraintStore::OR, x, 2, NULL);
  BuildObject **x = new BuildObject*[3];  
  x[0] = x_;
  x[1] = y_;
  return CSP::_Constraint(ConstraintStore::OR, x, 2, NULL);
}
BuildObject* CSP::_IfThenElse( BuildObject *x_, BuildObject *y_, BuildObject *z_ )
{
  BuildObject **x = new BuildObject*[4];
  x[0] = x_;
  x[1] = y_;
  x[2] = z_;
  return CSP::_Constraint(ConstraintStore::IFTHENELSE, x, 3, NULL);
}
BuildObject* CSP::_Disjunctive( BuildObject *x_, const int dx, BuildObject *y_, const int dy, const int t )
{
  BuildObject **x = new BuildObject*[3];  
  x[0] = x_;
  x[1] = y_;
  int *p = new int[3];
  p[0] = dx;
  p[1] = dy;
  p[2] = t;
  return CSP::_Constraint(ConstraintStore::DISJUNCTIVE, x, 2, p);
}
BuildObject* CSP::_Overlap( BuildObject *x_, const int dx, BuildObject *y_, const int dy )
{
  BuildObject **x = new BuildObject*[3];  
  x[0] = x_;
  x[1] = y_;
  int *p = new int[2];
  p[0] = dx;
  p[1] = dy;
  return CSP::_Constraint(ConstraintStore::OVERLAP, x, 2, p);
}
BuildObject* CSP::_Precedence( BuildObject *x_, const int d, BuildObject *y_ )
{
//   BuildObject **x = new BuildObject*[3];  
//   x[0] = x_.var_ptr_;
//   x[1] = y_.var_ptr_;
//   int *p = new int[1];
//   p[0] = d;
//   var_ptr_ = CSP::_Constraint(ConstraintStore::PRECEDENCE, x, 2, p);
  BuildObject **x = new BuildObject*[3];  
  x[0] = x_;
  x[1] = y_;
  int *p = new int[1];
  p[0] = d;
  return CSP::_Constraint(ConstraintStore::PRECEDENCE, x, 2, p);
}
BuildObject* CSP::_Precedence( BuildObject *x_, const int v )
{
  BuildObject **x = new BuildObject*[3];
  x[0] = x_;
  x[1] = CSP::_Variable(v,v);
  int *p = new int[1];
  p[0] = 0;
  return CSP::_Constraint(ConstraintStore::PRECEDENCE, x, 2, p);
}
BuildObject* CSP::_Precedence( const int v, BuildObject *x_ )
{
  BuildObject **x = new BuildObject*[3];
  x[1] = x_;
  x[0] = CSP::_Variable(v,v);
  int *p = new int[1];
  p[0] = 0;
  return CSP::_Constraint(ConstraintStore::PRECEDENCE, x, 2, p);
}
BuildObject* CSP::_Add( BuildObject *x_, BuildObject *y_ )
{
//   BuildObject **x = new BuildObject*[3];
//   x[0] = x_;
//   x[1] = y_;
//   int *p = new int[3];
//   p[0] = 1;
//   p[1] = 1;
//   p[2] = 0;
//   return CSP::_Constraint(ConstraintStore::SUM, x, 2, p);
  //cout << "here" << endl;
  //exit(0);
  //  if( y_->size() == 1 && )
  return CSP::_Sum( x_, y_, 1 );
}
BuildObject* CSP::_Sub( BuildObject *x_, BuildObject *y_ )
{
//   BuildObject **x = new BuildObject*[3];
//   x[0] = x_;
//   x[1] = y_;
//   int *p = new int[3];
//   p[0] = 1;
//   p[1] = -1;
//   p[2] = 0;
//   return CSP::_Constraint(ConstraintStore::SUM, x, 2, p);
  return CSP::_Sum( x_, y_, -1 );
}
BuildObject* CSP::_Mul( BuildObject *x_, BuildObject *y_ )
{
//   BuildObject **x = new BuildObject*[3];  
//   x[0] = x_.var_ptr_;
//   x[1] = y_.var_ptr_;
//   var_ptr_ = CSP::_Constraint(ConstraintStore::MUL, x, 2, NULL);
  BuildObject **x = new BuildObject*[3];  
  x[0] = x_;
  x[1] = y_;
  return CSP::_Constraint(ConstraintStore::MUL, x, 2, NULL);
}
BuildObject* CSP::_Mul( BuildObject *x_, const int v )
{
  BuildObject **x = new BuildObject*[3];  
  x[0] = x_;
  x[1] = CSP::_Variable(v,v);
  return CSP::_Constraint(ConstraintStore::MUL, x, 2, NULL);
}
BuildObject* CSP::_Div( BuildObject *x_, BuildObject *y_ )
{
//   BuildObject **x = new BuildObject*[3];  
//   x[0] = x_.var_ptr_;
//   x[1] = y_.var_ptr_;
//   var_ptr_ = CSP::_Constraint(ConstraintStore::DIV, x, 2, NULL);
  BuildObject **x = new BuildObject*[3];  
  x[0] = x_;
  x[1] = y_;
  return CSP::_Constraint(ConstraintStore::DIV, x, 2, NULL);
}
BuildObject* CSP::_Mod( BuildObject *x_, BuildObject *y_ )
{
//   BuildObject **x = new BuildObject*[3];  
//   x[0] = x_.var_ptr_;
//   x[1] = y_.var_ptr_;
//   var_ptr_ = CSP::_Constraint(ConstraintStore::MOD, x, 2, NULL);
  BuildObject **x = new BuildObject*[3];  
  x[0] = x_;
  x[1] = y_;
  return CSP::_Constraint(ConstraintStore::MOD, x, 2, NULL);
}
BuildObject* CSP::_Mod( BuildObject *x_, const int v )
{
//   BuildObject **x = new BuildObject*[3];  
//   x[0] = x_.var_ptr_;
//   x[1] = CSP::_Variable(v,v);
//   var_ptr_ = CSP::_Constraint(ConstraintStore::MOD, x, 2, NULL);
  BuildObject **x = new BuildObject*[3];  
  x[0] = x_;
  x[1] = CSP::_Variable(v,v);
  return CSP::_Constraint(ConstraintStore::MOD, x, 2, NULL);
}
BuildObject* CSP::_Equal( BuildObject *x_, BuildObject *y_, const int spin )
{
//   BuildObject **x = new BuildObject*[3];  
//   x[0] = x_.var_ptr_;
//   x[1] = y_.var_ptr_;
//   int *p = new int[2];
//   p[0] = spin;
//   p[1] = NOVAL;
//   var_ptr_ = CSP::_Constraint(ConstraintStore::EQUAL, x, 2, p);
  BuildObject **x = new BuildObject*[3];  
  x[0] = x_;
  x[1] = y_;
  int *p = new int[2];
  p[0] = spin;
  p[1] = NOVAL;
  return CSP::_Constraint(ConstraintStore::EQUAL, x, 2, p);
}
BuildObject* CSP::_Equal( BuildObject *x_, const int v, const int spin )
{
//   BuildObject **x = new BuildObject*[2];  
//   x[0] = x_.var_ptr_;
//   int *p = new int[2];
//   p[0] = spin;
//   p[1] = v;
//   var_ptr_ = CSP::_Constraint(ConstraintStore::EQUAL, x, 1, p);
  BuildObject **x = new BuildObject*[2];  
  x[0] = x_;
  int *p = new int[2];
  p[0] = spin;
  p[1] = v;
  return CSP::_Constraint(ConstraintStore::EQUAL, x, 1, p);
}
BuildObject* CSP::_Member( BuildObject *x_, const BitSet& k, const int spin )
{
//   BuildObject **x = new BuildObject*[2];
//   x[0] = x_.var_ptr_;
//   int *p = new int[2+k.size()];
//   p[0] = spin;
//   p[1] = 2;
//   BitsetIterator bit(k);
//   do p[p[1]++] = bit;
//   while( ++bit );
//   var_ptr_ = CSP::_Constraint(ConstraintStore::MEMBER, x, 1, p);
  BuildObject **x = new BuildObject*[2];
  x[0] = x_;
  int *p = new int[2+k.size()];
  p[0] = spin;
  p[1] = 2;
  BitsetIterator bit(k);
  do p[p[1]++] = bit;
  while( ++bit );
  return CSP::_Constraint(ConstraintStore::MEMBER, x, 1, p);
}
// BuildObject* CSP::_Sum( BuildObject *x_, const int o )
// {
// //   BuildObject **x = new BuildObject*[2];
// //   int *p = new int[7];
// //   p[0] = 1;
// //   p[1] = k;
// //   x[0] = x_.var_ptr_;
// //   var_ptr_ = CSP::_Constraint(ConstraintStore::SUM, x, 1, p );
//   BuildObject **x = new BuildObject*[2];
//   int *p = new int[7];
//   p[0] = 1;
//   p[1] = o;
//   x[0] = x_;
//   return CSP::_Constraint(ConstraintStore::SUM, x, 1, p );
// }
// BuildObject* CSP::_Sum( BuildObject **x_, const int n, const int *os )
// {
// //   int n = x_.size();
// //   int *p = new int[n+6]; //w[n]->offset, w[n+1]/w[n+2]->lb/ub, w[n+3]->weighted?, w[n+4]->Boolean?
// //   for(int i=0; i<n; ++i)
// //     p[i] = (os ? os[i] : 1);
// //   p[n] = (os ? os[n] : 0);
// //   BuildObject **x = x_.getArgs();
// //   var_ptr_ = CSP::_Constraint(ConstraintStore::SUM, x, n, p);
//   int i, *p = new int[n+6];
//   for(i=0; i<n; ++i)
//     p[i] = (os ? os[i] : 1);
//   p[n] = (os ? os[n] : 0);
//   return CSP::_Constraint(ConstraintStore::SUM, x_, n, p);
// }
// BuildObject* CSP::_Sum( BuildObject *x_, BuildObject *y_, const int t )
// {
// //   BuildObject **x = new BuildObject*[3];
// //   x[0] = x_.var_ptr_;
// //   x[1] = y_.var_ptr_;
// //   int *p = new int[8];
// //   p[0] = 1;
// //   p[1] = t;
// //   p[2] = 0;
// //   var_ptr_ = CSP::_Constraint(ConstraintStore::SUM, x, 2, p);
//    BuildObject **x = new BuildObject*[3];
//    x[0] = x_;
//    x[1] = y_;
//    int *p = new int[8];
//    p[0] = 1;
//    p[1] = t;
//    p[2] = 0;
//    return CSP::_Constraint(ConstraintStore::SUM, x, 2, p);
// }
BuildObject* CSP::_Sum( BuildObject *x_, const int o )
{
//   BuildObject **x = new BuildObject*[2];
//   int *p = new int[7];
//   p[0] = 1;
//   p[1] = k;
//   x[0] = x_.var_ptr_;
//   var_ptr_ = CSP::_Constraint(ConstraintStore::SUM, x, 1, p );
  BuildObject **x = new BuildObject*[2];
  int *p = new int[8];
  p[1] = 1;
  p[2] = -1;
  p[3] = -o;
  //p[3] = -o;
  x[0] = x_;
  return CSP::_Constraint(ConstraintStore::SUM, x, 1, p );
}
BuildObject* CSP::_Sum( BuildObject **x_, const int n, const int *os )
{
//   int n = x_.size();
//   int *p = new int[n+6]; //w[n]->offset, w[n+1]/w[n+2]->lb/ub, w[n+3]->weighted?, w[n+4]->Boolean?
//   for(int i=0; i<n; ++i)
//     p[i] = (os ? os[i] : 1);
//   p[n] = (os ? os[n] : 0);
//   BuildObject **x = x_.getArgs();
//   var_ptr_ = CSP::_Constraint(ConstraintStore::SUM, x, n, p);

  int i, *p = new int[n+7];
  for(i=0; i<n; ++i) {
    p[i+1] = (os ? os[i] : 1);
  }

  p[n+1] = -1;
  p[n+2] = (os ? -os[n] : 0);

  //p[n+2] = (os ? -os[n] : 0);
  BuildObject *res = CSP::_Constraint(ConstraintStore::SUM, x_, n, p);

  return res;
}
BuildObject* CSP::_Sum( BuildObject *x_, BuildObject *y_, const int t )
{
//   BuildObject **x = new BuildObject*[3];
//   x[0] = x_.var_ptr_;
//   x[1] = y_.var_ptr_;
//   int *p = new int[8];
//   p[0] = 1;
//   p[1] = t;
//   p[2] = 0;
//   var_ptr_ = CSP::_Constraint(ConstraintStore::SUM, x, 2, p);
   BuildObject **x = new BuildObject*[3];
   x[0] = x_;
   x[1] = y_;
   int *p = new int[9];
   p[1] = 1;
   p[2] = t;
   p[3] = -1;
   p[4] = 0;
   //p[4] = 0;
   return CSP::_Constraint(ConstraintStore::SUM, x, 2, p);
}
BuildObject* CSP::_Clause( BuildObject **x_, const int n, const int *pol )
{
  BuildObject **x = new BuildObject*[n];
  int *lits = new int[n];
  for(int i=0; i<n; ++i) {
    x[i] = x_[i];
    lits[i] = pol[i];
  }
  // BuildObjectPredicate *pred = 
//   pred->deReference();
//   return pred;
  return CSP::_Constraint(ConstraintStore::CLAUSE, x, n, lits);
  //cl_scope.push(x);
  //cl_lit.push(lits);
  //cl_size.push(n);
}
BuildObject* CSP::_Min( BuildObject *x_, BuildObject *y_ )
{
  BuildObject **x = new BuildObject*[3];
  x[0] = x_;
  x[1] = y_;
  return CSP::_Constraint(ConstraintStore::MIN, x, 2, NULL );
}
BuildObject* CSP::_Max( BuildObject *x_, BuildObject *y_ )
{
  BuildObject **x = new BuildObject*[3];
  x[0] = x_;
  x[1] = y_;
  return CSP::_Constraint(ConstraintStore::MAX, x, 2, NULL );
}
BuildObject* CSP::_Min( BuildObject **x_, const int n )
{
//   int n = x_.size();
//   BuildObject **x = x_.getArgs();
//   var_ptr_ = CSP::_Constraint(ConstraintStore::MIN, x, n, NULL );
  return CSP::_Constraint(ConstraintStore::MIN, x_, n, NULL );
}
BuildObject* CSP::_Max( BuildObject **x_, const int n )
{
//   int n = x_.size();
//   BuildObject **x = x_.getArgs();
//   var_ptr_ = CSP::_Constraint(ConstraintStore::MAX, x, n, NULL );
  return CSP::_Constraint(ConstraintStore::MAX, x_, n, NULL );
}
BuildObject* CSP::_AllDifferent( BuildObject **x_, const int n )
{
  int *p = new int[2];
  p[0] = 2;
  p[1] = 0;
  return CSP::_Constraint(ConstraintStore::ALLDIFF, x_, n, p );
}
BuildObject* CSP::_AllDifferent( BuildObject **x_, const int n, const int c, const int k )
{
//   int *p = new int[1];
//   p[0] = c;
//   var_ptr_ = CSP::_Constraint(ConstraintStore::ALLDIFF, x_.getArgs(), x_.size(), p );
  int *p = new int[2];
  p[0] = c;
  p[1] = k;
  return CSP::_Constraint(ConstraintStore::ALLDIFF, x_, n, p );
}
BuildObject* CSP::_Cardinality( BuildObject **x_, const int n, const int k )
{
//   int *p = new int[1];
//   p[0] = k;
//   var_ptr_ = CSP::_Constraint(ConstraintStore::CARD, x_.getArgs(), x_.size(), p );
  int *p = new int[1];
  p[0] = k;
  return CSP::_Constraint(ConstraintStore::CARD, x_, n, p );
}
BuildObject* CSP::_Element( BuildObject **x_, const int n, const int k )
{
  int *p = new int[1];
  p[0] = k;
  return CSP::_Constraint(ConstraintStore::ELEMENT, x_, n, p );
}
BuildObject* CSP::_Gcc( BuildObject **x_, const int n,
			const int firstDomainValue, 
			const int lastDomainValue,
			const int* minOccurrences,
			const int* maxOccurrences )
{
  int m = (lastDomainValue-firstDomainValue+1);
  int *p = new int[2+2*m];
  p[0] = firstDomainValue;
  p[1] = lastDomainValue;  
  for(int i=0; i<m; ++i) {
    p[2+i] = minOccurrences[i];
    p[2+m+i] = maxOccurrences[i];
  }
  return CSP::_Constraint(ConstraintStore::GCC, x_, n, p );
}
BuildObject* CSP::_LexOrder( BuildObject **x_, const int n, const int eq )
{
  int *p = new int[1];
  p[0] = eq;
  return CSP::_Constraint(ConstraintStore::LEX, x_, n, p );
}

BuildObject* CSP::_Table( BuildObject **x_, const int nx, const bool sp )
{
  int cidx = ENVIRONMENT.addConstraint( new BuildObjectTable(nx) ) ;
  BuildObjectTable* tab = (BuildObjectTable*)(ENVIRONMENT[cidx]) ;
  tab->spin = sp;

  int *p = new int[5];
  return new BuildObjectPredicate( x_, nx, -NOVAL/2, NOVAL/2, tab, p ) ;
}
BuildObject* CSP::_Table( BuildObject **x_, const int nx, 
			 const int tidx, const bool sp )
{
  BuildObjectTable* tab = (BuildObjectTable*)(ENVIRONMENT[tidx]) ;
  tab->spin = sp;

  int *p = new int[5];
  return new BuildObjectPredicate( x_, nx, -NOVAL/2, NOVAL/2, tab, p ) ;
}

void CSP::_addTuple( BuildObject *tab, const int* sol )
{
  ((BuildObjectTable*)(((BuildObjectPredicate*)tab)->relation))->add(sol);
}

BuildObjectObjective* CSP::_Maximise( BuildObject *x )
{
  return new BuildObjectMaximiseVar( x ) ; 
}

BuildObjectObjective* CSP::_Minimise( BuildObject *x )
{
  return new BuildObjectMinimiseVar( x ) ; 
}


void CSP::printPython() const
{
  print( cout );
}

void CSP::print(std::ostream& o) const
{

  int i, n=declarations.size;
  for(i=0; i<n; ++i) {
    
    //cout << declarations[i]->id << " ";
    //cout << declarations[i]->getVarId() << " ";
    if( declarations[i]->isReferenced() || declarations[i]->isRel() ) {
      cout << " _" << declarations[i]->isReferenced() << "_ ";
      declarations[i]->print( cout );
      cout << endl;
    }
  }
  
}

int CSP::size() 
{
  return declarations.size ;
}

void CSP::printXML(std::ostream& o) const
{
  int nVariables = 0;
  int nConstraints = 0;
  int skippedConstraints = 0;
  int skippedVariables = 0;
  int maxarity = 0;

  int i, j=declarations.size, k, l;

  int var2dec[j], con2dec[j], dec2var[j], subconstraint[j];
  for( i=0; i<j; ++i )
    { 
      subconstraint[i] = 0;
      if( declarations[i]->isRel() ) {
	if( declarations[i]->isConstraint() ) {
	  con2dec[nConstraints++] = i;	
	  if( declarations[i]->isReferenced() ) {
	    dec2var[declarations[i]->getVarId()] = nVariables;
	    var2dec[nVariables++] = i;
	  }
	}
	k = ((BuildObjectPredicate*)(declarations[i]))->arity;
	if( k > maxarity )
	  maxarity = k;
      } else {
	cerr << "MUT FIND OUT WHAT WAS THAT!" << endl;
	exit(0);
	//if( declarations[i]->isReferenced() && (!declarations[i]->isPar() || declarations[i]->size() > 1) ) {
	//dec2var[declarations[i]->getVarId()] = nVariables;
	//var2dec[nVariables++] = i;
	//}
      }
    }

  int var2dom[nVariables];
  int con2pre[nConstraints];
  int con_type[nConstraints];
  map<string, int> domains_map;
  map<string, int> predicates_map;
  vector<string> domains;
  vector<string> predicates;
  vector<int> domainsize;
  vector<int> predicatesize;

  vector<string> global;
  vector<int> table;

  int nDomains=0;
  int nPredicates=0;
  int nparams=0;

  i=nConstraints;
  while( i-- ) {
    if( !subconstraint[con2dec[i]] ) {

      nparams=0;
      string pred = declarations[con2dec[i]]->xmlPred( nparams, 0, subconstraint );           
      subconstraint[con2dec[i]] = 0;
      
      if( pred.substr(0,6) == "global" ) {
	//cout << "GLOBAL: " << pred << endl;
	global.push_back( pred );
	con_type[i] = 2;
      } else if( pred == "table" ) {
	//cout << "TABLE: " << pred << endl;
	table.push_back( i );
	con_type[i] = 1;
      } else {
	//cout << "PREDICATE: " << pred << endl;
	if( declarations[con2dec[i]]->isReferenced() )
	  pred = "eq(" + pred + ",X" + int2string(nparams++) + ")";
	con_type[i] = 0;
	if( predicates_map.find(pred) != predicates_map.end() ) {
	  con2pre[i] = predicates_map[pred];
	} else {
	  predicates_map[pred] = nPredicates;
	  predicates.push_back(pred);
	  predicatesize.push_back(nparams);
	  con2pre[i] = nPredicates++;
	}
      }
    } else {
      ++skippedConstraints;
    }
  }



  for(i=0; i<nVariables; ++i) {
    string dom = declarations[var2dec[i]]->xmlDom();
    if( domains_map.find( dom ) != domains_map.end() ) {
      var2dom[i] = domains_map[dom];      
    } else {
      domains_map[dom] = nDomains;
      domains.push_back( dom );
      domainsize.push_back( declarations[var2dec[i]]->size() );

      var2dom[i] = nDomains++;
    }
  }


  // header	  
  o << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" 
    << endl << endl << "<instance>" << endl 
    << "<presentation maxConstraintArity=\""
    << maxarity <<"\" format=\"XCSP 2.1\"/>" 
    << endl << endl;
	
  // domains
  o << "  <domains nbDomains=\"" << nDomains << "\">" << endl;
  for(i=0; i<nDomains; ++i)
    {
      o << "    <domain name=\"D" << i 
	<< "\" nbValues=\"" << domainsize[i] << "\">" 
	<< domains[i]
	<< "</domain>" << endl; 
    }
  o << "  </domains>" << endl << endl;

  // variables
  i=nVariables;
  while( i-- )
    if( subconstraint[var2dec[i]] ) {
      var2dec[i] = var2dec[--nVariables];
      var2dom[i] = var2dom[nVariables];
      dec2var[var2dec[i]] = i;
      //++skippedVariables;
    }
  k=0;
  o << "  <variables nbVariables=\"" << (nVariables-skippedVariables) << "\">" << endl;
  for(i=0; i<nVariables; ++i)
    {
      //if( !subconstraint[var2dec[i]] )
      o << "    <variable name=\"V" << i
	<< "\" domain=\"D" << var2dom[i] << "\"/>" << endl;
    }
  o << "  </variables>" << endl << endl;
  
  // predicates
  if( nPredicates ) {
    o << "  <predicates nbPredicates=\"" 
      << (nPredicates) << "\">" << endl;
    for(i=0; i<nPredicates; ++i)
      {
	o << "    <predicate name=\"P"<< i << "\">" << endl
	  << "      <parameters>";
	for(j=0; j<predicatesize[i]; ++j)
	  o << " int X" << j;
	o << "</parameters>" << endl
	  << "      <expression> " << endl
	  << "        <functional>" 
	  << predicates[i] 
	  << "</functional>" << endl
	  << "      </expression>" << endl
	  << "    </predicate>" << endl;
      }
    o << "  </predicates>" << endl << endl;
  }

  // constraints
  int arity;
  BuildObjectPredicate *bcon;
  BitSet scope(0, nVariables-1, BitSet::empt);
  o << "  <constraints nbConstraints=\"" 
    << nConstraints-skippedConstraints << "\">" << endl;
  k=0;
  l=0;

  i=nConstraints;
  while( i-- ) 
    if( !subconstraint[con2dec[i]] )
      {
	bcon = ((BuildObjectPredicate*)(declarations[con2dec[i]]));
	
	/// PREDICATES
	if( con_type[i] == 0 ) {
	  scope.clear();
	  
	  string params = declarations[con2dec[i]]->xmlCon( scope, dec2var ); 
	  if( declarations[con2dec[i]]->isReferenced() )
	    {
	      if( declarations[con2dec[i]]->size() > 1 ) 
		{
		  j = dec2var[declarations[con2dec[i]]->getVarId()];
		  scope.insert(j);
		  params += (" V" + int2string(j));
		}
	      else
		{
		  params += (" " + int2string(declarations[con2dec[i]]->min()));
		}
	    }
	  arity = scope.size();
	  if( arity ) {
	    o << "    <constraint name=\"C" << l++ <<"\" arity=\""
	      << arity << "\" scope=\"";
	    j=scope.min();
	    while(true) {
	      o << "V" << j;//dec2var[j] ;
	      j = scope.next( j );
	      if( j != NOVAL )
		o << " ";
	      else break;
	    } 
	  o << "\" reference=\"P" << con2pre[i] << "\">" << endl
	    << "      <parameters>" 
	    << params
	    << "</parameters> " << endl
	    << "    </constraint>" << endl; 
	  }
	}
	
	/// RELATIONS
	else if( con_type[i] == 1 ) {
	  ;
	}
      
	/// GLOBAL CONSTRAINTS
	else if( con_type[i] == 2 ) {

// 	  cout << i << " -- " << bcon->arity << " ";
// 	  bcon->print( cout );
// 	  cout << " " << global[k] << endl;

	  scope.clear();
	  arity = bcon->arity;
	  for(j=0; j<arity; ++j) {
	    scope.insert(dec2var[bcon->scope[j]->getVarId()]);
	  }
	  if( bcon->isReferenced() ) scope.insert( dec2var[bcon->scope[arity]->getVarId()] );
	  arity = scope.size();
	  o << "    <constraint name=\"C" << l++ <<"\" arity=\""
	    << arity << "\" scope=\"";
	  j = scope.min();
	  do { 
	    o << "V" << j << " ";
	    j = scope.next( j );
	  } while(j != NOVAL);

	  o << "\" reference=\"" << global[k] ;
	  
	  if( global[k] == "global:element" )
	    {

	      o << "\">" << endl
		<< "      <parameters>" << endl
		<< "        V" << dec2var[bcon->scope[bcon->arity-1]->getVarId()] 
		<< endl << "        [ ";
	      for(j=0; j<bcon->arity-1; ++j)
		o << "V" << dec2var[bcon->scope[j]->getVarId()] << " ";
	      o << "]\n        V" 
		<< dec2var[bcon->scope[bcon->arity]->getVarId()] 
		<< endl << "      </parameters>"
		<< endl << "    </constraint>" << endl;
	    }
	  else if( global[k] == "global:cumulative" )
	    {
	      o << "\"/>" << endl;
	    }	  
	  else if( global[k] == "global:weightedSum" )
	    {
	      o << "\"/>" << endl;
	    }
	  else if( global[k] == "global:allDifferent" )
	    {
	      o << "\"/>" << endl;
	    }
	  ++k;
	}
      }
  o << "  </constraints>" 
    << endl << endl
    << "</instance>" << endl;

}



int ConstraintStore::addConstraint(BuildObjectConstraint *c)
{
  store.push( c );
  return (store.size-1);
}

void ConstraintStore::remConstraint(int idx)
{
  delete store[idx];
  store[idx] = NULL;
  while( !store.back() )
    --store.size;
}

void ConstraintStore::flush()
{
  int i = store.size;
  while( i-- > NUMCONS ) {    
    if( store[i] ) 
      delete store[i]; 
    store.pop();
  }
}

BuildObjectConstraint* ConstraintStore::getConstraint(const int c)
{
  if( !store[c] ) {
    switch( c ) {
    case ABS         : store[c] = new BuildObjectAbs();              break;
    case ALLDIFF     : store[c] = new BuildObjectAllDiff();          break;
    case AND         : store[c] = new BuildObjectAnd();              break;
    case DISJUNCTIVE : store[c] = new BuildObjectDisjunctive();      break;
    case OVERLAP     : store[c] = new BuildObjectOverlap();          break;
    case DIV         : store[c] = new BuildObjectDiv();              break;
    case ELEMENT     : store[c] = new BuildObjectElement();          break;
    case EQUAL       : store[c] = new BuildObjectEqual();            break;
    case MEMBER      : store[c] = new BuildObjectMember();           break;
    case GCC         : store[c] = new BuildObjectGcc();              break;
    case LEX         : store[c] = new BuildObjectLex();              break;
    case DLEX        : store[c] = new BuildObjectDLex();             break;
    case MAX         : store[c] = new BuildObjectMax();              break;
    case MIN         : store[c] = new BuildObjectMin();              break;
    case CARD        : store[c] = new BuildObjectCardinality();      break;
    case MOD         : store[c] = new BuildObjectMod();              break;
    case MUL         : store[c] = new BuildObjectMul();              break;
    case NEG         : store[c] = new BuildObjectNegation();         break;
    case NOT         : store[c] = new BuildObjectNot();              break;
    case OR          : store[c] = new BuildObjectOr();               break;
    case PRECEDENCE  : store[c] = new BuildObjectPrecedence();       break;
    case SUM         : store[c] = new BuildObjectSum();              break;
    case IFTHENELSE  : store[c] = new BuildObjectIfThenElse();       break;
    case TREE        : store[c] = new BuildObjectTree();             break;
    case TDAG        : store[c] = new BuildObjectTDAG();             break;
    case CLAUSE      : store[c] = new BuildObjectClause();           break;
      //case REGULAR     : store[c] = new BuildObjectRegular();          break;

    }
  }
  return store[c];
}


ConstraintStore::ConstraintStore()
{
  store.init(NUMCONS, 1024);
  for(int i=0; i<NUMCONS; ++i)
    store[i] = NULL;
}

ConstraintStore::~ConstraintStore()
{
  int i = store.size;
  while( i-- )
    if( store[i] )
      delete store[i];
}

BuildObjectConstraint *ConstraintStore::operator[]( const int c ) 
{
  return getConstraint( c );
}

Objective::Objective(BuildObject *x) : X(x) 
{ 
  obj = NULL; 
}

Maximise::Maximise(Variable x) : Objective( x.var_ptr_ ) 
{
  obj = new BuildObjectMaximiseVar( X ) ; 
}

Minimise::Minimise(Variable x) : Objective( x.var_ptr_ ) 
{
  obj = new BuildObjectMinimiseVar( X ) ;
}


int registerConstraint(BuildObjectConstraint*c)
{
  return ENVIRONMENT.addConstraint( c );
}

BuildObjectConstraint *getConstraint(const int idx)
{
  return ENVIRONMENT[idx];
}
