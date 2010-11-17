
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

#include <mistral_sol.h>
#include <mistral_con.h>
#include <mistral_csp.h>
#include <mistral_mod.h>
//#include <sat.h>

using namespace Mistral;


//VariableBool is_true(NULL,1);
//VariableBool is_false(NULL,0);

/**********************************************
 * Relations described in extension
 **********************************************/

/**********************************************
 * Unary Constraint
 **********************************************/

UnaryConstraint::UnaryConstraint(Solver *s, VariableInt* x, const char* n)
  : X(x)
{
  solver = s;
  name = n;
}

UnaryConstraint::~UnaryConstraint()
{
}

void UnaryConstraint::activate()
{
  solver->unaryCons.push( this );
}


/**********************************************
 * Unary More Constraint
 **********************************************/

UnaryConstraintMore::UnaryConstraintMore(Solver *s, VariableInt* x, int lb)
  : UnaryConstraint( s, x, "lower_bound" )
{
  bound = lb;
}

bool UnaryConstraintMore::propagate()
{
  return X->setMin( bound );
}

void UnaryConstraintMore::print(std::ostream& o) const
{
  X->print(o);
  o << " >= " << bound;
}

/**********************************************
 * Unary Less Constraint
 **********************************************/

UnaryConstraintLess::UnaryConstraintLess(Solver *s, VariableInt* x, int ub)
  : UnaryConstraint( s, x, "upper_bound" )
{
  bound = ub;
}

bool UnaryConstraintLess::propagate()
{
  return X->setMax( bound );
}

void UnaryConstraintLess::print(std::ostream& o) const
{
  X->print(o);
  o << " <= " << bound;
}


/**********************************************
 * GAC2001 Extensional Constraint
 **********************************************/

ConstraintGAC2001Allowed::ConstraintGAC2001Allowed(Solver* s, VariableInt** v, const int n)
: Constraint(s, v, n, DOMAINTRIGGER) 
{
  int i, j, nval;

  firstSupport = new ReversibleNum<int>*[arity];
  themins = new int[arity];
  D_X = new BitSet[arity];
  order = new int[arity];

  supportList = new Vector<int*>*[arity];

  for(i=0; i<arity; ++i) {
    order[i] = i;
    themins[i] = scope[i]->min();
    nval=(scope[i]->max() - themins[i] + 1);

//     /// init the list of supports for each value
//     supportList[i] = new Vector<int*> [scope[i]->max() - themins[i] + 1];
//     supportList[i] -= themins[i];
//     // reserve one index for the residual support
//     DomainIterator *valit = scope[i]->begin();
//     do {
//       j = *valit;
//       supportList[i][j].push( NULL );
//     } while( valit->next() );

    /// init the reversible data structure 'firstSupport'
    firstSupport[i] = new ReversibleNum<int>[nval];
    firstSupport[i] -= themins[i];
    DomainIterator *valit = scope[i]->begin();
    do {
      j = *valit;
      s->binds( firstSupport[i][j] );
      firstSupport[i][j].setValue(1);
    } while( valit->next() );

    /// init the domain shortcuts
    switch( scope[i]->getType() ) {
    case VariableInt::CONST : { 
      D_X[i].init( scope[i]->min(), scope[i]->max(), BitSet::full); 
    } break;
    case VariableInt::RANGE : { 
      std::cerr << "not handled" << std::endl; exit(0); 
    } break;
    case VariableInt::INTLIST : { 
      std::cerr << "not handled" << std::endl; exit(0); 
    } break;
    case VariableInt::BOOL : { 
      D_X[i].pointTo( (unsigned int*)(&(((VariableBool*)scope[i])->domain.state)) ); 
    } break;
    default : { 
      D_X[i].pointTo( ((VariableDomain*)scope[i])->values ); 
    }
    }
  }

  isClone = false;
}


void ConstraintGAC2001Allowed::init( Vector<int*>& tuples, const int spin, ConstraintGAC2001Allowed* con )
{

  int relevant, i, j;
  DomainIterator *valit;

  for(i=0; i<arity; ++i) {
    /// init the list of supports for each value
    supportList[i] = new Vector<int*> [scope[i]->max() - themins[i] + 1];
    supportList[i] -= themins[i];
    // reserve one index for the residual support
    valit = scope[i]->begin();
    do {
      j = *valit;
      supportList[i][j].push( NULL );
    } while( valit->next() );
  }

  i=tuples.size;
  while( i-- ) {
    relevant = true;
    for(j=0; relevant && j<arity; ++j) {
      relevant = D_X[j].member( tuples[i][j] );
    }
    if( relevant )
      for(j=0; j<arity; ++j)
	supportList[j][tuples[i][j]].push( tuples[i] );	   
  }
}


// ConstraintGAC2001Allowed::ConstraintGAC2001Allowed(Solver *s, VariableInt **v,
// 						   const int n, Vector<int*>** m)
//   : Constraint(s, v, n, DOMAINTRIGGER) 
// { 
//   //std::cout << "arity: " << arity << std::endl;

//   firstSupport = new ReversibleNum<int>*[arity];
//   supportList = m;
//   themins = new int[arity];
//   D_X = new BitSet[arity];
//   order = new int[arity];
//   int i, j, k, l, nTuples, nval;
//   for(i=0; i<arity; ++i) {
//     order[i] = i;
//     themins[i] = scope[i]->min();
//     nval=(scope[i]->max() - themins[i] + 1);
    
//     //std::cout << (scope[i]->max()) << " " << themins[i] << " " << nval << std::endl; 

//     firstSupport[i] = new ReversibleNum<int>[nval];
//     firstSupport[i] -= themins[i];
//     DomainIterator *valit = scope[i]->begin();
//     do {
//       j = *valit;

//       //firstSupport[i][j].init(s);
//       s->binds( firstSupport[i][j] );
//       firstSupport[i][j].setValue(1);
//     } while( valit->next() );

//     switch( scope[i]->getType() ) {
//     case VariableInt::CONST : { 
//       D_X[i].init( scope[i]->min(), scope[i]->max(), BitSet::full); 
//     } break;
//     case VariableInt::RANGE : { 
//       std::cerr << "not handled" << std::endl; exit(0); 
//     } break;
//     case VariableInt::INTLIST : { 
//       std::cerr << "not handled" << std::endl; exit(0); 
//     } break;
//     case VariableInt::BOOL : { 
//       D_X[i].pointTo( (unsigned int*)(&(((VariableBool*)scope[i])->domain.state)) ); 
//     } break;
//     default : { 
//       D_X[i].pointTo( ((VariableDomain*)scope[i])->values ); 
//     }
//     }
//   }
// }


ConstraintGAC2001Allowed::~ConstraintGAC2001Allowed() 
{ 
  int i;
  for(i=0; i<arity; ++i) {
    if( scope[i]->getType() != VariableInt::CONST ) {
      D_X[i].table = NULL;
      D_X[i].neg_words = 0;
    }
  }
  for(i=0; i<arity; ++i) {
    firstSupport[i] += themins[i];
    delete [] firstSupport[i];
  }
  delete [] firstSupport;
  delete [] D_X;
  delete [] order;

  if( !isClone ) {
    for(i=0; i<arity; ++i) {
      supportList[i] += themins[i];
      delete [] supportList[i];
    }
    delete [] supportList;
  }

  delete [] themins;
}


int ConstraintGAC2001Allowed::getNumSupports(const int var, const int val) const
{
  return (supportList[var][val].size - firstSupport[var][val]);
}


bool ConstraintGAC2001Allowed::propagate(const int changedIdx, const int e) 
{
  int consistent = true, i, oi, j, k, ok, valid, index, n;
  //VariableInt* not_consistent = NULL;
  //int i, oi, j, k, ok, valid, index, n;

  for(i=1; i<arity; ++i) {   // order variables by increasing domain size
    j=i;
    while( j && scope[order[j]]->domsize() < scope[order[j-1]]->domsize() ) {      
      ok = order[j];
      order[j] = order[j-1];
      order[--j] = ok;
    }
  }

  DomainIterator *valit;
  for( i=0; consistent && i<arity; ++i ) {
    oi = order[i];
    if( oi != changedIdx && scope[oi]->isLinked() ) 
      {
	valit = scope[oi]->begin(); // revise the domain of scope[oi]
	do {
	  j = *valit; // for each value j	  

	  n = supportList[oi][j].size;
 	  supports_X = supportList[oi][j].stack_; // init the list of supports
	  valid = false;
	  for(index=firstSupport[oi][j]; !valid && index<n; ++index) {
	    valid = true;
	    for(k=0; valid && k<arity; ++k) {
	      ok = order[k];
	      valid = ( oi == ok || D_X[ok].fastMember( supports_X[index][ok] ) );
	      //valid = ( oi == ok || D_X[ok].member( supports_X[index][ok] ) );
	    }
	  }
	  if( valid ) firstSupport[oi][j] = index-1;
	  else consistent = scope[oi]->remove( j );
	} while( consistent && valit->next() );
      }
  }

  return consistent;
}


int ConstraintGAC2001Allowed::check( const int* s ) const 
{
  //  int i, equal=0;
//   initSupport( 0, s[0] );
//   while( !equal && index-- )
//     {
//       i=arity;
//       equal=1;
//       while( equal && --i )
// 	equal = (s[i] == support[index][i]);
//     }
//  return !equal;
  return 0;
}

void ConstraintGAC2001Allowed::print(std::ostream& o) const
{
//   if( id == 6 ) {
//     //printlong(o);
//     std::cout << supportList[7][23][4][0] << std::endl;
//   }

  o << "GAC2001Ext";
  Constraint::print( o );

  std::cout << std::endl;
  int j,n,oi,index,k;
  DomainIterator *valit;

  for( oi=0; oi<arity; ++oi ) {
    valit = scope[oi]->begin(); // revise the domain of scope[oi]
    do {
      j = *valit; // for each value j 
      scope[oi]->printshort( std::cout );      
      std::cout << " = " << j << "\n";
      n = supportList[oi][j].size;
      for(index=//0;
	  firstSupport[oi][j]; 
	  index<n; ++index) {
	std::cout << "\t";
	//if( index == firstSupport[oi][j] ) std::cout << "* ";
	for(k=0; k<arity; ++k) {
	  std::cout << ((char)(supportList[oi][j][index][k]+65)) ;
	}
	std::cout << std::endl;
      }
    } while( valit->next() );
  }  
}


void ConstraintGAC2001Allowed::printlong(std::ostream& o) const
{
  o << id << "GAC2001Ext";
  Constraint::print( o );

  //int j,n,oi,index,k;
  //DomainIterator *valit;
  int oi;

  for(oi=0; oi<arity; ++oi) {
    scope[oi]->print( std::cout );
    std::cout << "  ";
  }
  std::cout << std::endl;


//   for( oi=0; oi<arity; ++oi ) {
//     valit = scope[oi]->begin(); // revise the domain of scope[oi]
//     do {
//       j = *valit; // for each value j 
//       scope[oi]->printshort( std::cout );      
//       std::cout << " = " << j << "\n";
//       n = supportList[oi][j].size;
//       for(index=//0;
// 	  firstSupport[oi][j]; 
// 	  index<n; ++index) {
// 	std::cout << "\t";
// 	//if( index == firstSupport[oi][j] ) std::cout << "* ";
// 	for(k=0; k<arity; ++k) {
// 	  std::cout << ((char)(supportList[oi][j][index][k]+65)) ;
// 	}
// 	std::cout << std::endl;
//       }
//     } while( valit->next() );
//   }  

}



/**********************************************
 * Extensional Constraint
 **********************************************/

int ConstraintGAC3Valid::getpos(const int *sol) const 
{
  int pos = sol[arity-1], i = arity-1;
  while( i-- ) pos += (sol[i])*var_sizes[i]; 
  return pos;
}

ConstraintGAC3Valid::ConstraintGAC3Valid(Solver* s, VariableInt** v, const int n)
  : Constraint(s, v, n, DOMAINTRIGGER) 
{
//   int i;
//   var_sizes = new int[n-1];
//   var_sizes[n-2] = (scope[n-1]->max() - scope[n-1]->min() + 1);
//   for(i=n-3; i >= 0; --i) {
//     var_sizes[i] = var_sizes[i+1] * 
//       (scope[i+1]->max() - scope[i+1]->min() + 1);
//     }
//   int matrixmin = scope[n-1]->min();
//   int matrixmax = scope[n-1]->max();
    
//   i = n-1;
//   while( i-- ) {
//     matrixmax += (scope[i]->max())*var_sizes[i]; 
//     matrixmin += (scope[i]->min())*var_sizes[i]; 
//   }
	
//   matrix.init(matrixmin, matrixmax, BitSet::empt);

  isClone = false;
}

bool ConstraintGAC3Valid::isValid(const int *tuple) const
{
  int i=arity, valid=true;
  while(valid && i--)
    valid = scope[i]->contain(tuple[i]);
  return valid;
}

void ConstraintGAC3Valid::init( Vector<int*>& tuples, const int spin, ConstraintGAC3Valid *con )
{
  if( con == this ) {
    int i;
    int matrixmin = scope[arity-1]->min();
    int matrixmax = scope[arity-1]->max();
    
    var_sizes = new int[arity-1];
    var_sizes[arity-2] = (scope[arity-1]->max() - scope[arity-1]->min() + 1);
    for(i=arity-3; i >= 0; --i) {
      var_sizes[i] = var_sizes[i+1] * 
	(scope[i+1]->max() - scope[i+1]->min() + 1);
    }
    
    i = arity-1;
    while( i-- ) {
      matrixmax += (scope[i]->max())*var_sizes[i]; 
      matrixmin += (scope[i]->min())*var_sizes[i]; 
    }
    
    matrix.init(matrixmin, matrixmax, BitSet::empt);
    
    if( spin ) {
      matrix.fill();
      i=tuples.size;
      while( i-- ) {
	if(isValid(tuples[i]))
	  matrix.erase(getpos(tuples[i]));
      }
      
    } else {
      matrix.clear();
      i=tuples.size;
      while( i-- )
	if(isValid(tuples[i]))
	  matrix.insert(getpos(tuples[i]));
    }
    isClone = false;
  } else {
    isClone = true;
    matrix.pos_words = con->matrix.pos_words;
    matrix.neg_words = con->matrix.neg_words;
    matrix.table = con->matrix.table;
    var_sizes = con->var_sizes;
  }
}


// ConstraintGAC3Valid::ConstraintGAC3Valid(Solver *s, VariableInt **v,
// 					     const int n, BitSet& m,
// 					     int *vs)
//   : Constraint(s, v, n, DOMAINTRIGGER) 
// {

//   //std::cout << "VALID" << std::endl;

//   matrix.pos_words = m.pos_words;
//   matrix.neg_words = m.neg_words;
//   matrix.table = m.table;
//   var_sizes = vs;
// }

ConstraintGAC3Valid::~ConstraintGAC3Valid() 
{ 
  if( isClone ) {
    matrix.table = NULL;
    matrix.neg_words = 0;
  } else {
    delete [] var_sizes;
  }
}

int ConstraintGAC3Valid::check( const int* s ) const 
{
  return matrix.member(getpos(s));
}

void ConstraintGAC3Valid::fillNogood() 
{
  matrix.fill();
}

void ConstraintGAC3Valid::fillSupport() 
{
  matrix.clear();
}

void ConstraintGAC3Valid::addNogood(const int* sol) 
{
  matrix.insert(getpos(sol));
}

void ConstraintGAC3Valid::addSupport(const int* sol) 
{
  matrix.erase(getpos(sol));
}

void ConstraintGAC3Valid::print(std::ostream& o) const
{
  o << "Ext";
  Constraint::print( o );
}


/**********************************************
 * Extensional Binary Constraint
 **********************************************/

void ConstraintAC3Bitset::init( )
{
  int i=2, j;
  BitSet aux[2];
  while( i-- ) 
    aux[i].init( scope[1-i]->minCapacity(), 
		 scope[1-i]->maxCapacity(),
		 BitSet::empt );

  i=2;
  while( i-- ) {
    size[i] = scope[i]->domsize();

    residualSupport[i] = new int[scope[i]->maxCapacity() - 
				 scope[i]->minCapacity()];
    residualSupport[i] -= scope[i]->minCapacity();

    minSup[i] = NOVAL;

    DomainIterator *valit = scope[i]->begin();
    do {
      residualSupport[i][*valit] = supportList[i][*valit].pos_words-1;	

      aux[i].clear();
      scope[1-i]->unionTo( aux[i] );
      aux[i].intersectWith( supportList[i][*valit] );

      j = aux[i].size();

      if( minSup[i] > j ) minSup[i] = j;
    } while( valit->next() );
  }
}

ConstraintAC3Bitset::ConstraintAC3Bitset(Solver* s, VariableInt** v, const int spin)
  : Constraint(s,v,2,DOMAINTRIGGER)
{
//   int i, j;
//   i=2;
//   while( i-- ) { 
//     supportList[i] = new BitSet[scope[i]->maxCapacity() - scope[i]->minCapacity()];
//     supportList[i]-= scope[i]->minCapacity();
//     j = scope[i]->maxCapacity();
//     while( j-- > scope[i]->minCapacity() )
//       supportList[i][j].init( scope[1-i]->min(),
// 			      scope[1-i]->max(),
// 			      (spin ? 0 : 0xffffffff) );
//   }
  isClone = false;
}

void ConstraintAC3Bitset::init( Vector<int*>& tuples, const int spin, ConstraintAC3Bitset* con )
{
  if( con == this ) {
    int i=2, j;
    while( i-- ) { 
      supportList[i] = new BitSet[scope[i]->maxCapacity() - scope[i]->minCapacity()];
      supportList[i]-= scope[i]->minCapacity();
      j = scope[i]->maxCapacity();
      while( j-- > scope[i]->minCapacity() )
	supportList[i][j].init( scope[1-i]->min(),
				scope[1-i]->max(),
				(spin ? 0 : 0xffffffff) );
    }
    
    i=tuples.size;
    while( i-- ) {
      j=2;
      while( j-- ) 
	if( scope[j]->contain(tuples[i][j]) ) {
	  if( spin ) 
	    supportList[j][tuples[i][j]].insert( tuples[i][1-j] );
	  else 
	    supportList[j][tuples[i][j]].erase( tuples[i][1-j] );
	}
    }
    isClone = false;
  } else {
    isClone = true;
    supportList[0] = con->supportList[0];
    supportList[1] = con->supportList[1];
  }
  init();
}

// ConstraintAC3Bitset::ConstraintAC3Bitset(Solver *s, VariableInt **v, BitSet **sl)
//     : Constraint(s,v,2,DOMAINTRIGGER)
// {
//   supportList[0] = sl[0];
//   supportList[1] = sl[1];
// }

ConstraintAC3Bitset::~ConstraintAC3Bitset() 
{ 
  int i=2;
  while( i-- ) {
    residualSupport[i]+= scope[i]->minCapacity();
    delete [] residualSupport[i];
  }
  if( !isClone ) {
    i=2;
    while( i-- ) {
      supportList[i] += scope[i]->minCapacity();
      delete [] supportList[i];
    }
  }
}

int ConstraintAC3Bitset::check( const int* s ) const 
{
  return !(supportList[0][s[0]].member( s[1] ));
}

void ConstraintAC3Bitset::fillNogood() 
{
  int i=2;
  while( i-- ) {
    DomainIterator *valit = scope[i]->begin();
    do supportList[i][*valit].clear();
    while( valit->next() );
  }  
}

void ConstraintAC3Bitset::fillSupport() 
{
  int i=2;
  while( i-- ) {
    DomainIterator *valit = scope[i]->begin();
    do supportList[i][*valit].fill();
    while( valit->next() );
  }
}

void ConstraintAC3Bitset::addNogood(const int* sol) 
{
  int i=2;
  while( i-- ) {
    supportList[i][sol[i]].erase( sol[1-i] );
  }
}

void ConstraintAC3Bitset::addSupport(const int* sol) 
{
  int i=2;
  while( i-- ) {
    supportList[i][sol[i]].insert( sol[1-i] );
  }
}

bool ConstraintAC3Bitset::propagate(const int changedIdx, const int e) 
{
  int i = 1-changedIdx;
  bool consistent = true;
  if( (size[changedIdx] - scope[changedIdx]->domsize()) >= minSup[i] ) {
    consistent = scope[i]->revise( supportList[i], 
				   residualSupport[i], 
				   scope[changedIdx] );
  }
  return consistent;
}

void ConstraintAC3Bitset::print(std::ostream& o) const
{
  o << "BinExt";
  Constraint::print( o );
  int i=2;
  while( i-- ) {
    //o << " |";
    o << std::endl;
    DomainIterator *valit = scope[i]->begin();
    do {
      //o << " " << supportList[i][*valit].size();
      supportList[i][*valit].print( o );
      o << std::endl;
    } while( valit->next() );
  }
  o << std::endl;
}



/**********************************************
 * Clause Constraint
 **********************************************/

ConstraintClause::ConstraintClause(Solver *s, VariableInt** v, const int n, const int* c) :
  Constraint(s, v, n, VALUETRIGGER) 
{
  conflict = new int[n];
  int i = n;

  while( i-- ) {
    conflict[i] = c[i];
  }
}

ConstraintClause::~ConstraintClause() 
{ 
  delete [] conflict;
}

int ConstraintClause::check( const int* s ) const 
{
  int i = arity;
  while( i-- )
    if( s[i] != conflict[i] ) return 0;
  return 1;
}

bool ConstraintClause::propagate()
{
  bool consistent = true;
//   for(int i=0; i<arity; ++i)
//     if( scope[i]->isGround() )
//       consistent = propagate( scope[i], VALUETRIGGER );
  return consistent;
} 

bool ConstraintClause::propagate(const int changedIdx, const int e) 
{

  //if( scope[changedIdx] == watch[0] )

  int i=arity, j=-1;
  while( i-- ) {
    if( scope[i]->isLinked() ) {
      if( j < 0 ) j = i;
      else {
	return true;
      }
    } else if( scope[i]->value() != conflict[i] ) {
      return true;
    }
  }

  return scope[j]->remove( conflict[j] );
}

void ConstraintClause::print(std::ostream& o) const
{

  bool sat = false;
  int i=arity;

  while( !sat && i-- ) {
    if( scope[i]->isGround() && scope[i]->value() != conflict[i] )
      sat = true;
  }
  if( sat ) o << "T";
  else {
    o << "{ ";
    for(i=0; i<arity; ++i) {
      if( scope[i]->isGround() ) {
	if( scope[i]->value() !=  conflict[i] )
	  o << "T ";
      } else {
	scope[i]->printshort( o );
	o << "~" << conflict[i] << " ";
      }
    }
    o << "}";
  }
}


/**********************************************
 * =/= Constraint
 **********************************************/

ConstraintNotEqual::ConstraintNotEqual(Solver *s, VariableInt** v)
  : Constraint(s, v, 2, VALUETRIGGER)
{
}

ConstraintNotEqual::~ConstraintNotEqual() {}

int ConstraintNotEqual::check( const int* s ) const 
{
  return ( s[0] == s[1] );
}

bool ConstraintNotEqual::propagate()
{
  return ((!scope[1]->isGround() || scope[0]->remove(scope[1]->min())) &&
	  (!scope[0]->isGround() || scope[1]->remove(scope[0]->min())));
}

bool ConstraintNotEqual::propagate(const int changedIdx, const int e)
{
  return scope[1-changedIdx]->remove(scope[changedIdx]->min());  
}

void ConstraintNotEqual::print(std::ostream& o) const
{
  scope[0]->printshort( o ) ;
  o << " =/= " ;
  scope[1]->printshort( o ) ;
}


/**********************************************
 * =/= Constraint with wildcard
 **********************************************/

ConstraintNotEqualWildcard::ConstraintNotEqualWildcard(Solver *s, VariableInt** v, const int w)
  : Constraint(s, v, 2, VALUETRIGGER)
{
  wildcard = w;
}

ConstraintNotEqualWildcard::~ConstraintNotEqualWildcard() {}

int ConstraintNotEqualWildcard::check( const int* s ) const 
{
  return ( s[0] == s[1] && s[0] != wildcard );
}

bool ConstraintNotEqualWildcard::propagate()
{
  return ((!scope[0]->isGround() || propagate(0, VALUETRIGGER)) &&
  	  (!scope[1]->isGround() || propagate(1, VALUETRIGGER)));
}

bool ConstraintNotEqualWildcard::propagate(const int changedIdx, const int e)
{
  bool consistent = true;
  int idx = 1-changedIdx, w = scope[changedIdx]->min();
  if( w != wildcard )
    consistent = scope[idx]->remove(w);  
  else if( idx )
    consistent = scope[idx]->setDomain(wildcard);

  return consistent;
}

void ConstraintNotEqualWildcard::print(std::ostream& o) const
{
  scope[0]->print( o ) ;
  o << " =/= " ;
  scope[1]->print( o ) ;
  o << " or " ;
  scope[0]->print( o ) ;
  o << " = " ;
  o << wildcard ;
}


/**********************************************
 * Or Constraint
 **********************************************/

ConstraintOr::ConstraintOr(Solver *s, VariableInt** v)
  : Constraint(s, v, 2, DOMAINTRIGGER)
{
}

ConstraintOr::~ConstraintOr() {}

int ConstraintOr::check( const int* s ) const 
{
  return ( !s[0] && !s[1] );
}

bool ConstraintOr::propagate(const int changedIdx, const int e)
{
  return( scope[changedIdx]->max() || scope[1-changedIdx]->setDomain(1) );  
}

void ConstraintOr::print(std::ostream& o) const
{
  scope[0]->printshort( o ) ;
  o << " || " ;
  scope[1]->printshort( o ) ;
}

/**********************************************
 * NogoodBase Constraint
 **********************************************/
ConstraintNogoodBase::ConstraintNogoodBase(Solver* s) 
  : Constraint(s, s->variables.stack_, s->length, VALUETRIGGER, 0)
{

//   //std::cout << s->length << std::endl;

 
  int i, n; // = s->length;
//   scp = new VariableInt*[n];
//   for(i=0; i<n; ++i) {
//     scp[i] = s->variables[i];
//   }
//   initScope(scp, n, VALUETRIGGER, 0);

  weight = 0;
  minimums = new int[arity];
  watched = new Vector< Array < GenLiteral >* > *[arity];
  for(i=0; i<arity; ++i) {
    minimums[i] = scope[i]->minCapacity();
    n = (scope[i]->maxCapacity() - minimums[i] + 1);
    watched[i] = new Vector< Array < GenLiteral > * >[n];
    watched[i] -= minimums[i];
  }

  //std::cout << 11 << std::endl;
}

ConstraintNogoodBase::~ConstraintNogoodBase() {
  for(int i=0; i<arity; ++i) {
    watched[i] += minimums[i];
    delete [] watched[i];
  }
  /*
  Array<GenLiteral> *ngd;
  while(!nogood.empty()) {
    nogood.pop(ngd);
    delete ngd;
  }
  */
  delete [] minimums;
}

void ConstraintNogoodBase::add( Vector<GenLiteral>& conflict ) {
  //std::cout << "beg add" << std::endl;
  if(conflict.size > 1) {

//      std::cout << (conflict.size) << std::endl;
	      
//      conflict[0].print(std::cout);
//      std::cout << std::endl;

//      conflict[1].print(std::cout);
//      std::cout << std::endl;

//      std::cout << arity << std::endl;
      

    Array<GenLiteral> *cl = Array<GenLiteral>::Array_new(conflict);
    nogood.push( cl );
    watched[conflict[0].var->id][conflict[0].val].push(cl);
    watched[conflict[1].var->id][conflict[1].val].push(cl);
  } else {
    SimpleUnaryConstraint cons('r', conflict[0].val, conflict[0].var);
    scope[0]->solver->sUnaryCons.push(cons);
    //std::cout << cons.var->id << "," << cons.val << ":" ;
    //std::cerr << "unary nogood!!" << std::endl;
    //exit(0);
  }
  //std::cout << "end add" << std::endl;
}

int ConstraintNogoodBase::check( const int* ) const {
  return 0;
}

bool ConstraintNogoodBase::propagate() {
  return true;
}

bool ConstraintNogoodBase::propagate(const int changedIdx, const int e) {
  //   std::cout << "unit propag ";
   //   scope[changedIdx]->printshort(std::cout);
//   std::cout << " == " << scope[changedIdx]->value() << std::endl;

  int j, rank, consistent=true, v=scope[changedIdx]->value();
  int i=watched[changedIdx][v].size;

  //std::cout << i << std::endl;

  GenLiteral p;
  //Array< GenLiteral >  *clause = watched[changedIdx][v];
  while(consistent && i--) {

//     std::cout << "\t" << i << std::endl;
//     std::cout << "\t" << watched[changedIdx] << std::endl ;
//     //std::cout << "\t" << watched[changedIdx][v] << std::endl ;
//     watched[changedIdx][v].print(std::cout) ;
//     std::cout << std::endl;
//     std::cout << "\t" << watched[changedIdx][v][i] << std::endl << std::endl;

    //for(i=0; consistent && i<watched[changedIdx][v].size; ++i) {
    Array< GenLiteral >& cl = *(watched[changedIdx][v][i]);

//     std::cout << "watched by:" << std::endl;
//     for(j=0; j<cl.size; ++j) {
//       cl[j].print(std::cout);
//       std::cout << " ";
//     }
//     std::cout << std::endl;

    // first check if this nogood is already satisfied (first literal is satisfied)
    if(cl[0].satisfied() || cl[1].satisfied()) {
      //std::cout << " => satisfied by watchee" << std::endl;
      continue;
    }

    // otherwise, find out if the literal is the first, or second 
    rank = (cl[1].var == scope[changedIdx]);

    //std::cout << " the " << (rank ? "second" : "first") << " literal is the trigger " << std::endl;

//     if(rank && cl[1].satisfied()) {
//       // it is the second, and it is satisfied, we put it first.
//       p = cl[0];
//       cl[0] = cl[1];
//       cl[1] = p;
//       continue;
//     } else {
      // cl[rank] is violated, we need to replace it
    assert(scope[changedIdx]->id == cl[rank].var->id);
    assert(scope[changedIdx]->equal(cl[rank].val));
    
    for(j=2; j<cl.size; ++j) {
      // for each subsequent literal, if it non-violated, it is a possible watcher
      if(!cl[j].violated()) {
	// in which case we swap it with cl[rank]
	p = cl[j];
	cl[j] = cl[rank];
	cl[rank] = p;
	
	watched[p.var->id][p.val].push(watched[changedIdx][v][i]);
	watched[changedIdx][v].erase(i);
	
	break;
      } 
    }

    if(j>=cl.size) {
      // there wasn't any possible replacement, so we need to propagate:
      consistent = cl[1-rank].var->remove(cl[1-rank].val);
    }
  }
  //}

  //  std::cout << "end propag" << std::endl;
  
  return consistent;
}

void ConstraintNogoodBase::print(std::ostream& o) const {
  for(int i=0; i<nogood.size; ++i) {
    Array< GenLiteral >& cl = *(nogood[i]);
    for(int j=0; j<cl.size; ++j) {
      cl[j].print(o);
      o << " ";
    }
    o << std::endl;
  }
//   for(int i=0; i<arity; ++i) {
//     for(int j=minimums[i]; j<)
//   }
}


/**********************************************
 * GenNogoodBase Constraint
 **********************************************/
ConstraintGenNogoodBase::ConstraintGenNogoodBase(Solver* s) 
  : Constraint(s, s->variables.stack_, s->length, VALUETRIGGER, 0)
{
  int i, n; // = s->length;
  weight = 0;
  minimums = new int[arity];
  watched_values  = new Vector< Array < GeneralLiteral >* > *[arity];
  watched_domains = new Vector< Array < GeneralLiteral >* > [arity];
  watched_bounds  = new Vector< Array < GeneralLiteral >* > [arity];

  
  for(i=0; i<arity; ++i) {
    minimums[i] = 0;
    watched_values[i] = NULL;
  }

  for(i=0; i<arity; ++i) {
    minimums[i] = scope[i]->minCapacity();
    n = (scope[i]->maxCapacity() - minimums[i] + 1);

    if(n>10000) {
      std::cout << n <<std::endl;
      exit(1);
    }

    watched_values[i] = new Vector< Array < GeneralLiteral > * >[n];
    watched_values[i] -= minimums[i];
  }
}

ConstraintGenNogoodBase::~ConstraintGenNogoodBase() {
  for(int i=0; i<arity; ++i) {
    watched_values[i] += minimums[i];
    delete [] watched_values[i];
  }
  delete [] watched_values;
  delete [] watched_domains;
  delete [] watched_bounds;
  delete [] minimums;

  for(int i=0; i<nogood.size; ++i) {
    free(nogood[i]);
  }
}

void ConstraintGenNogoodBase::add( Vector<GeneralLiteral>& conflict ) {
  if(conflict.size > 1) {
    Array<GeneralLiteral> *cl = Array<GeneralLiteral>::Array_new(conflict);
    nogood.push( cl );
    
    int var, val;

    for(int i=0; i<2; ++i) {
      var = conflict[i].var->id;
      val = conflict[i].value();

      if(conflict[i].type() == GeneralLiteral::REMOVAL) {
// 	if(!watched_values[var]) {
// 	  minimums[var] = scope[var]->minCapacity();
// 	  n = (scope[var]->maxCapacity() - minimums[var] + 1);
// 	  watched_values[var] = new Vector< Array < GeneralLiteral > * >[n];
// 	  watched_values[var] -= minimums[var];
// 	}
	watched_values[var][val].push(cl);
      } else if(conflict[i].type() == GeneralLiteral::ASSIGNMENT) {
	watched_domains[var].push(cl);
      } else {
      	watched_bounds[var].push(cl);
      }
    }
  } else {
    //SimpleUnaryConstraint cons('r', conflict[0].value(), conflict[0].var);
    scope[0]->solver->sUnaryCons.push(conflict[0]);
    //std::cout << conflict[0].var->id << "," << conflict[0].value() << ":" ;
    //conflict[0].print(std::cout);
    //std::cerr << " is a unary nogood!!" << std::endl;
    //exit(0);
  }
  //std::cout << "end add" << std::endl;
}

void ConstraintGenNogoodBase::removeClause( const int idx ) {

  int var, val;

//   std::cout << "remove";
//   for(int i=0; i<nogood[idx]->size; ++i) {
//     std::cout << " ";
//     (*(nogood[idx]))[i].print(std::cout);
//   }
//   std::cout << std::endl;

  for(int i=0; i<2; ++i) {
    var = (*(nogood[idx]))[i].var->id;
    val = (*(nogood[idx]))[i].value();

    if((*(nogood[idx]))[i].type() == GeneralLiteral::REMOVAL) {
     
      
      //(*(nogood[idx]))[i].print(std::cout);

      //std::cout << std::endl; 

      watched_values[var][val].remove(nogood[idx]);

    } else if((*(nogood[idx]))[i].type() == GeneralLiteral::ASSIGNMENT) {

      //std::cout << "this should not happen" << std::endl; 

      watched_domains[var].remove(nogood[idx]);
    } else {

      //std::cout << "this should not happen" << std::endl; 

      watched_bounds[var].remove(nogood[idx]);
    }
  } 

  nogood.erase(idx);
  //delete nogood[idx];
  free(nogood[idx]);
}

void ConstraintGenNogoodBase::forget( const int lim ) {
  //std::cout << " FORGET! " << std::endl;
  //std::cout << nogood.size << " -> " << lim << std::endl;

  while(nogood.size>lim) {
    //std::cout << nogood.size << " -> " << lim << std::endl;
    removeClause(nogood.size-1);
  }

  //std::cout << nogood.size << " -> " << lim << std::endl;
}

int ConstraintGenNogoodBase::check( const int* ) const {
  return 0;
}

bool ConstraintGenNogoodBase::propagate() {
  return true;
}

bool ConstraintGenNogoodBase::propagate(const int changedIdx, const int e) {
  int j, rank, consistent=true, v=scope[changedIdx]->value();
  int i;

  if(e == VALUETRIGGER && watched_values[changedIdx]) {
    i=watched_values[changedIdx][v].size;
    GeneralLiteral p;
    while(consistent && i--) {
      Array< GeneralLiteral >& cl = *(watched_values[changedIdx][v][i]);
      if(cl[0].satisfied() || cl[1].satisfied()) {
	continue;
      }
      
      // otherwise, find out if the literal is the first, or second 
      rank = (cl[1].var == scope[changedIdx]);
      
      assert(scope[changedIdx]->id == cl[rank].var->id);
      assert(scope[changedIdx]->equal(cl[rank].value()));
      
      for(j=2; j<cl.size; ++j) {
	// for each subsequent literal, if it non-violated, it is a possible watcher
	if(!cl[j].violated()) {
	  // in which case we swap it with cl[rank]
	  p = cl[j];
	  cl[j] = cl[rank];
	  cl[rank] = p;
	  
	  watched_values[p.var->id][p.value()].push(watched_values[changedIdx][v][i]);
	  watched_values[changedIdx][v].erase(i);
	  
	  break;
	} 
      }
      
      if(j>=cl.size) {
	// there wasn't any possible replacement, so we need to propagate:
	consistent = cl[1-rank].var->remove(cl[1-rank].value());
      }
    }
  } else if(e == RANGETRIGGER) {
    std::cout << "this should not happen" << std::endl; 
  } else if(e == DOMAINTRIGGER) {
    std::cout << "this should not happen" << std::endl; 
  }

  //  std::cout << "end propag" << std::endl;
  
  return consistent;
}

void ConstraintGenNogoodBase::print(std::ostream& o) const {
  for(int i=0; i<nogood.size; ++i) {
    Array< GeneralLiteral >& cl = *(nogood[i]);
    for(int j=0; j<cl.size; ++j) {
      cl[j].print(o);
      o << " ";
    }
    o << std::endl;
  }
//   for(int i=0; i<arity; ++i) {
//     for(int j=minimums[i]; j<)
//   }
}




/**********************************************
 * inverse Constraint
 **********************************************/
ConstraintInverse::ConstraintInverse(Solver *s, VariableInt** v, int i, int j)
  : Constraint(s, v, 2, DOMAINTRIGGER)
{
  i_[0] = i;
  i_[1] = j;
}

ConstraintInverse::~ConstraintInverse() {}

int ConstraintInverse::check( const int* s ) const 
{
  return ( (s[0] == i_[0]) != (s[1] == i_[1]) );
}

bool ConstraintInverse::propagate(const int changedIdx, const int e)
{
  int idx = 1-changedIdx;
  if( !scope[changedIdx]->contain( i_[changedIdx] ) ) {
    if( !scope[idx]->remove( i_[idx] ) ) return false;
  } else {
    if( scope[changedIdx]->isGround() && 
	!scope[idx]->setDomain( i_[idx] ) ) return false;
  }
  return true;  
}

void ConstraintInverse::print(std::ostream& o) const
{
  o << "(" ;
  scope[0]->printshort( o ) ;
  o << " = " << i_[0] << ") <=> (";
  scope[1]->printshort( o ) ;
  o << " = " << i_[1] << ")" ;
}


/**********************************************
 * Lex Constraint
 **********************************************/
// this constraints checks one coordinate i of the sequences x[] and y[]
// scope[0] = x[i]
// scope[1] = y[i]
// scope[2] = b[i]
// scope[3] = b[i+1]
// the constraint enforces 3 things:
// 1/ b[i] <= b[i+1]
// 2/ (b[i] < b[i+1]) => (x[i] < y[i])
// 3/ (b[i] = b[i+1] = 0) => (x[i] = y[i])
ConstraintLex::ConstraintLex(Solver *s, VariableInt** v)
  : Constraint(s, v, 4, DOMAINTRIGGER),
    //domain_b1( ((VariableBool*)(v[2]))->domain ),
    //domain_b2( ((VariableBool*)(v[3]))->domain )
    domain_b1(v[2]->getIntDomain()),
    domain_b2(v[3]->getIntDomain())
{  
}

ConstraintLex::~ConstraintLex() 
{
}
		   
int ConstraintLex::check( const int* s ) const 
{
  return !( ( (s[2] < s[3]) && (s[0] < s[1]) )
	    ||
	    ( (s[2] + s[3] == 0) && (s[0] == s[1]) )
	    ||
	    (s[2] + s[3] == 2) );
}

// 1/ b[i] <= b[i+1], b[i] has changed
inline bool ConstraintLex::rule1forward()
{
  bool consistent = scope[3]->setMin( scope[2]->min() );
//   std::cout << "rule1f" << std::endl;
//   print( std::cout );
//   std::cout << std::endl;
  return consistent;
}
// 1/ b[i] <= b[i+1], b[i+1] has changed
inline bool ConstraintLex::rule1backward()
{
  bool consistent = scope[2]->setMax( scope[3]->max() );
//   std::cout << "rule1b" << std::endl;
//   print( std::cout );
//   std::cout << std::endl;
  return consistent;
}

// 2/ (b[i] = b[i+1] = 0) => (x[i] = y[i])
//    (b[i] < b[i+1]) => (x[i] < y[i]), b[i] or b[i+1] have changed
inline bool ConstraintLex::rule2forward()
{
  bool consistent = true;
  //if( scope[2]->max() ) return true;

  if( !scope[2]->max() ) {
    int x=scope[3]->min();
    consistent = ( scope[1]->setMin( scope[0]->min()+x ) &&
		   scope[0]->setMax( scope[1]->max()-x ) );
  }

//   std::cout << "rule2f" << std::endl;
//   print( std::cout );
//   std::cout << std::endl;
  return consistent;
}
// 2/ (b[i] = b[i+1] = 0) => (x[i] = y[i])
//    (b[i] < b[i+1]) => (x[i] < y[i]), x[i] or y[i] have changed
inline bool ConstraintLex::rule2backward()
{
  bool consistent = true;
  int diff = (scope[0]->min() - scope[1]->max());
  if( !diff ) // we therefore have x[i] >= y[i], hence b[i] = b[i+1]
    consistent = ( scope[2]->setMin( scope[3]->min() ) &&
		   scope[3]->setMax( scope[2]->max() ) );
  else if( diff > 0 ) // if x[i] > y[i] then we should have: b[i] = b[i+1] = 1
    consistent = (scope[2]->setDomain(1) && scope[3]->setDomain(1));
  // else, we cannot conclude

//   std::cout << "rule2b" << std::endl;
//   print( std::cout );
//   std::cout << std::endl;
  return consistent;
}



// bool ConstraintLex::propagate(const int changedIdx, const int e)
// { 
//   bool consistent = true;

// //   bool print_it = (scope[0]->id == 3 && scope[1]->id == 93 && !(scope[2]->min()) && !(scope[3]->min()) && (scope[0]->min() != scope[1]->min() || scope[0]->max() != scope[1]->max()));
// //   if(print_it)
// //     {
// //       scope[0]->print(std::cout);
// //       std::cout << " < ";
// //       scope[1]->print(std::cout);
// //       std::cout << std::endl;

// //       scope[2]->print(std::cout);
// //       std::cout << " | ";
// //       scope[3]->print(std::cout);
// //       std::cout << std::endl;
// //     }

//   // prunes B[i+1]
//   // B[i] <= B[i+1]
//   if( (scope[2]->min() && !(scope[3]->remove(0)))
//       ||
//       (!scope[3]->max() && !(scope[2]->remove(1)))
//       ) 
//     consistent = false;


//   if(consistent && !scope[2]->max() && !scope[3]->max()) {
//     consistent = (scope[0]->setDomain(scope[1]) &&
// 		  scope[1]->setDomain(scope[0]));
//   }



//   // irrelevant part of the vectors
//   if( consistent && !scope[2]->min() /*&& scope[3]->min()*/ ) {

//     // x[i] == y[i] -> b[i]+b[i+1] > 0 -> b[i+1] = 0
//     if( scope[0]->isGround() && scope[1]->isGround() && 
// 	!(scope[2]->max()) &&
// 	(scope[0]->value() == scope[1]->value()) &&
// 	!(scope[3]->remove(1)) ) 
//       consistent = false;

//     if( consistent ) {
//       // prunes B[i]
//       // x[i] > y[i] -> b[i] > 0
//       if( scope[0]->min() > scope[1]->max() &&
// 	  ( !(scope[2]->setDomain(1)) ||
// 	    !(scope[3]->setDomain(1)) ) ) 
// 	consistent = false;
      
//       if( consistent ) {
// 	// x[i] >= y[i] -> b[i] >= b[i+1]
// 	if( scope[0]->min() >= scope[1]->max() &&
// 	    ( !(scope[2]->setMin( scope[3]->min() )) || 
// 	      !(scope[3]->setMax( scope[2]->max() )) )
// 	    ) 
// 	  consistent = false;
	
// 	if( consistent ) {
// 	  // x[i] < y[i] -> b[i+1] = 1
// 	  if( scope[0]->max() < scope[1]->min() )
// 	    if( !scope[3]->setDomain(1) ) 
// 	      consistent = false;
	  
// 	  if( consistent ) {
// 	    // prunes X[i] and Y[i]
// 	    if( !(scope[2]->max()) ) {
// 	      if( scope[3]->min() ) {
// 		if( !(scope[0]->setMax(scope[1]->max()-1)) ) 
// 		  consistent = false;
// 		else if( !(scope[1]->setMin(scope[0]->min()+1)) ) 
// 		  consistent = false;
// 	      }  else {
// 		if( !(scope[0]->setMax(scope[1]->max())) ) 
// 		  consistent = false;
// 		else if( !(scope[1]->setMin(scope[0]->min())) ) 
// 		  consistent = false;
// 	      }
// 	    }
// 	  }
// 	}
//       }
//     }
//   }


// //   if(print_it)
// //     {
// //       scope[0]->print(std::cout);
// //       std::cout << " < ";
// //       scope[1]->print(std::cout);
// //       std::cout << std::endl;

// //       scope[2]->print(std::cout);
// //       std::cout << " | ";
// //       scope[3]->print(std::cout);
// //       std::cout << std::endl << consistent << std::endl << std::endl;
// //     }


// //   bool consistent = true;
// //   //   if( v == scope[2] )
// //   //     consistent = (rule1forward() && rule2forward());
// //   //   else if( v == scope[3] )
// //   //     consistent = (rule1backward() && rule2forward());
// //   //   else
// //   //     consistent = rule2backward();
// //   //   return consistent;


// // //   int before = (scope[0]->domsize() +
// // // 		scope[1]->domsize() +
// // // 		scope[2]->domsize() +
// // // 		scope[3]->domsize());
// // //   print( std::cout );
// // //   std::cout << std::endl;


// //   // prunes B[i+1]
// //   // B[i] <= B[i+1]
// // //   if( scope[2]->min() && !(scope[3]->remove(0)) 
// // //       ||
// // //       !scope[3]->max() && !(scope[2]->remove(1))
// // //       ) 
// // //     consistent = false;
// //   if( domain_b1.state == 2 && !(scope[3]->remove(0)) 
// //       ||
// //       domain_b2.state == 1 && !(scope[2]->remove(1))
// //       ) 
// //     consistent = false;



// // //return false;

// //   // irrelevant part of the vectors
// //   //  if( consistent && !scope[2]->min() /*&& scope[3]->min()*/ ) {
// //   if( consistent && domain_b1.state != 2 /*&& scope[3]->min()*/ ) {
// // //     //std::cout << "relax, relax, ohoh, ohoh" << std::endl;
// // //     //relax();
// // //     return true;
// // //   }
  
// //   // prunes B[i+1]
// //   //   // x[i] =/= y[i] -> b[i]+b[i+1] > 0 -> b[i+1] > 0
// //   //   if( !(scope[0]->intersect( scope[1] )) &&
// //   //       !(scope[3]->remove(0)) ) return false;
// //   // x[i] == y[i] -> b[i]+b[i+1] > 0 -> b[i+1] = 0
// //     if( scope[0]->isGround() && scope[1]->isGround() && 
// // 	//!(scope[2]->max()) &&
// // 	domain_b1.state == 1 &&
// // 	(scope[0]->value() == scope[1]->value()) &&
// // 	!(scope[3]->remove(1)) ) 
// //       consistent = false;
// //       //return false;

// //     if( consistent ) {
// //       // prunes B[i]
// //       // x[i] > y[i] -> b[i] > 0
// //       if( scope[0]->min() > scope[1]->max() &&
// // 	  ( !(scope[2]->setDomain(1)) ||
// // 	    !(scope[3]->setDomain(1)) ) ) 
// // 	consistent = false;
// // 	//return false;
      
// //       if( consistent ) {
// // 	// x[i] >= y[i] -> b[i] >= b[i+1]
// // 	if( scope[0]->min() >= scope[1]->max() &&
// // 	    // 	    ( !(scope[2]->setMin( scope[3]->min() )) || 
// // 	    // 	      !(scope[3]->setMax( scope[2]->max() )) )
// // 	    // 	    ) 
// // 	    ( !(scope[2]->setMin( !(domain_b2.state & 1)  )) || 
// // 	      !(scope[3]->setMax( (domain_b1.state >> 1) )) )
// // 	    ) 
// // 	  consistent = false;
// // 	//return false;
	
// // 	if( consistent ) {
// // 	  // x[i] < y[i] -> b[i+1] = 1
// // 	  if( scope[0]->max() < scope[1]->min() )
// // 	    if( !scope[3]->setDomain(1) ) 
// // 	      consistent = false;
// // 	      //return false;
	  
// // 	  if( consistent ) {
// // 	    // prunes X[i] and Y[i]
// // // 	    if( !(scope[2]->max()) )
// // // 	      if( scope[3]->min() ) {
// // 	    if( domain_b1.state == 1 )
// // 	      if( domain_b2.state == 2 ) {
// // 		if( !(scope[0]->setMax(scope[1]->max()-1)) ) 
// // 		  consistent = false;
// // 		//return false;
// // 		else if( !(scope[1]->setMin(scope[0]->min()+1)) ) 
// // 		  consistent = false;
// // 		//return false;
// // 	      }  else {
// // 		if( !(scope[0]->setMax(scope[1]->max())) ) 
// // 		  consistent = false;
// // 		//return false;
// // 		else if( !(scope[1]->setMin(scope[0]->min())) ) 
// // 		  consistent = false;
// // 		//return false;
// // 	      }
// // 	  }
// // 	}
// //       }
// //     }
// //   }


// // //   print( std::cout );
// // //   std::cout << " " << (consistent ? "ok" : "inconsistent") << std::endl;

// // //   if( before > (scope[0]->domsize() +
// // // 		scope[1]->domsize() +
// // // 		scope[2]->domsize() +
// // // 		scope[3]->domsize()) )
// // //     std::cout << "***" ;

// // //   std::cout << std::endl;

//   return consistent;
//   //return true;
// }





bool ConstraintLex::propagate(const int changedIdx, const int e)
{ 
  bool consistent = true;

  // prunes B[i+1]
  // B[i] <= B[i+1]
  if( (scope[2]->min() && !(scope[3]->remove(0)))
      ||
      (!scope[3]->max() && !(scope[2]->remove(1)))
      ) 
    consistent = false;

  //if(consistent && !scope[2]->max() && !scope[3]->max()) {
  if(consistent && !(domain_b1 >> 1) && !(domain_b2 >> 1)) {
    consistent = (scope[0]->setDomain(scope[1]) &&
		  scope[1]->setDomain(scope[0]));
  }

  // irrelevant part of the vectors
  if( consistent && !scope[2]->min() /*&& scope[3]->min()*/ ) {
    //if( consistent && !(domain_b1 & 1) /*&& scope[3]->min()*/ ) {

    // x[i] == y[i] -> b[i]+b[i+1] > 0 -> b[i+1] = 0
    if( scope[0]->isGround() && scope[1]->isGround() && 
	//if( (domain_b1 != 3) && (domain_b2 != 3) && 
	!(scope[2]->max()) &&
	(scope[0]->value() == scope[1]->value()) &&
	//(domain_b1 == domain_b2) &&
	!(scope[3]->remove(1)) ) 
      consistent = false;

    if( consistent ) {
      // prunes B[i]
      // x[i] > y[i] -> b[i] > 0
      if( scope[0]->min() > scope[1]->max() &&
	  ( !(scope[2]->setDomain(1)) ||
	    !(scope[3]->setDomain(1)) ) ) 
	consistent = false;
      
      if( consistent ) {
	// x[i] >= y[i] -> b[i] >= b[i+1]
	if( scope[0]->min() >= scope[1]->max() &&
	    ( !(scope[2]->setMin( scope[3]->min() )) || 
	      !(scope[3]->setMax( scope[2]->max() )) )
	    ) 
	  consistent = false;
	
	if( consistent ) {
	  // x[i] < y[i] -> b[i+1] = 1
	  if( scope[0]->max() < scope[1]->min() )
	    if( !scope[3]->setDomain(1) ) 
	      consistent = false;
	  
	  if( consistent ) {
	    // prunes X[i] and Y[i]
	    if( !(scope[2]->max()) ) {
	      if( scope[3]->min() ) {
		if( !(scope[0]->setMax(scope[1]->max()-1)) ) 
		  consistent = false;
		else if( !(scope[1]->setMin(scope[0]->min()+1)) ) 
		  consistent = false;
	      }  else {
		if( !(scope[0]->setMax(scope[1]->max())) ) 
		  consistent = false;
		else if( !(scope[1]->setMin(scope[0]->min())) ) 
		  consistent = false;
	      }
	    }
	  }
	}
      }
    }
  }


  return consistent;
}

void ConstraintLex::print(std::ostream& o) const 
{

  
  

//   o << "[" << scope[0]->min() << " "
//     << scope[0]->max() << "] < ["
//     << scope[1]->min() << " "
//     << scope[1]->max() << "] subject to b[i]={"
//     << scope[2]->min();
//   if( !scope[2]->isGround() )
//     o << "," << scope[2]->max();
//   o << "} and b[i+1]={"
//     << scope[3]->min();
//   if( !scope[3]->isGround() )
//     o << "," << scope[3]->max();
//   o << "}";


  scope[0]->printshort( o ); 
  o << " < ";
  scope[1]->printshort( o ); 
  o << " subject to ";
  scope[2]->printshort( o ); 
  o << " < ";
  scope[3]->printshort( o ); 


//   o << id << " ( (" ;
//   scope[2]->printshort( o ) ; 
//   o << " < " ;
//   scope[3]->printshort( o ) ;
//   o << " -> (" ;
//   scope[0]->printshort( o ) ;
//   o << " < " ; 
//   scope[1]->printshort( o ) ; 
//   o << ") ) & ( (" ;
//   scope[2]->printshort( o ) ;
//   o << " + " ;
//   scope[3]->printshort( o ) ;
//   o << " = 0) -> " ;
//   scope[0]->printshort( o ) ; 
//   o << " = " ; 
//   scope[1]->printshort( o ) ;
//   o << ")";
}


/**********************************************
 * <= Constraint
 **********************************************/
ConstraintLess::ConstraintLess(Solver *s, VariableInt** v)
  : Constraint(s, v, 2, RANGETRIGGER, WEIGHTTWO)
{
  offset = 0;
}

ConstraintLess::ConstraintLess(Solver *s, VariableInt** v, int o)
  : Constraint(s, v, 2, RANGETRIGGER, WEIGHTTWO)
{
  offset = o;
}

ConstraintLess::~ConstraintLess() {}

int ConstraintLess::check( const int* s ) const 
{
  return ( s[0] + offset > s[1] );
}

bool ConstraintLess::propagate(const int changedIdx, const int e)
{
  return scope[1]->setMin(scope[0]->min() + offset) &&
    scope[0]->setMax(scope[1]->max() - offset);
}

void ConstraintLess::print(std::ostream& o) const
{
  scope[0]->printshort( o ) ;
  if( !offset )
    o << " <= ";
  else if( offset == 1 )
    o << " < ";
  else 
    o << " + " << offset-1 << " < " ;
  scope[1]->printshort( o ) ;
}



/**********************************************
 * TDAG Constraint 
 **********************************************/
// v contains n*n Boolean variables a
ConstraintTDAG::ConstraintTDAG(Solver* s, VariableInt** v, const int n,
			       VariableInt** interval1, 
			       VariableInt** interval2)
  : Constraint( s, v, n, VALUETRIGGER )
{
  int i, k, x, y;

  num_nodes = (int)(floor(sqrt(n*2)))+1;
  successor = new ReversibleIntList[num_nodes];
  predecessor = new ReversibleIntList[num_nodes];
  pair_index = new int*[num_nodes];
  //node_name = new int[num_nodes];
  intervals = new VariableInt*[num_nodes];
  node_index = new int[s->variables.size];
  std::fill(node_index, node_index+s->variables.size, -1);

  first_index = new int[arity];
  second_index = new int[arity];
  for(i=0; i<num_nodes; ++i) {
    pair_index[i] = new int[num_nodes];
   
    s->binds( successor[i] );
    successor[i].setValue(0, num_nodes-1, 0);
    s->binds( predecessor[i] );
    predecessor[i].setValue(0, num_nodes-1, 0);

    successor[i].insert(i);
    predecessor[i].insert(i);
  }
  
  //PredicateDisjunctive* d;
  i=0;
  for(k=0; k<arity; ++k) {
    
    x = interval1[k]->id;
    y = interval2[k]->id;
    
    if(node_index[x] < 0) {
      node_index[x] = i;
      //node_name[i] = x;
      intervals[i] = interval1[k];
      ++i;
    }
    if(node_index[y] < 0) {
      node_index[y] = i;
      intervals[i] = interval2[k];
      //node_name[i] = y;
      ++i;
    }
    
    pair_index[node_index[x]][node_index[y]] = k;
    pair_index[node_index[y]][node_index[x]] = k;
    first_index[k] = node_index[x];
    second_index[k] = node_index[y]; 

  }



//   for(k=0,i=0; i<num_nodes; ++i) {
//     for(j=i+1; j<num_nodes; ++j, ++k) {
//       pair_index[i][j] = k;
//       pair_index[j][i] = k;
//       first_index[k] = i;
//       second_index[k] = j; 
//     }

//     s->binds( successor[i] );
//     successor[i].setValue(0, num_nodes-1, 0);
//     s->binds( predecessor[i] );
//     predecessor[i].setValue(0, num_nodes-1, 0);
//   }

//   for(i=0; i<arity; ++i) {
//     std::cout << std::setw(3) << first_index[i]  
// 	      << std::setw(3) << second_index[i] << std::endl;  
//   }
//   exit(1);
}

ConstraintTDAG::~ConstraintTDAG()
{
  delete [] successor;
  delete [] predecessor;
  for(int i=0; i<num_nodes; ++i) 
    delete [] pair_index[i];
  delete [] pair_index;
  delete [] first_index;
  delete [] second_index;
}

inline int ConstraintTDAG::check( const int* ) const 
{
  return 0; // TODO
}

inline bool ConstraintTDAG::propagate()
{
  return true;
}

inline bool ConstraintTDAG::propagate(const int changedIdx, const int e)
{
  bool consistent = true;
  unsigned int i, j, x, y;
  int node_i = -1;
  int node_j = -1;

  if(e == VALUETRIGGER) {
    node_i = first_index[changedIdx];
    node_j = second_index[changedIdx];
    
    if(scope[changedIdx]->value() == 0) {
      // node_i has to be before node_j
      //successor[node_i].insert(node_j);

      // node_j has to be after node_i
      //predecessor[node_j].insert(node_i);

      if(!successor[node_i].member(node_j)) {

	std::cout << "register " << intervals[node_i]->id+1 << " -> " 
		  << intervals[node_j]->id+1 << std::endl;

	// every successor of node_j (included) 
	// is also successor of 
	// every predecessor of node_i (included)
	for(i=0; consistent && i<predecessor[node_i].size; ++i)
	  for(j=0; consistent && j<successor[node_j].size; ++j) {
	    x = predecessor[node_i][i];
	    y = successor[node_j][j];
	    //successor[x].insert(y);
	    //predecessor[y].insert(x);
	    
	    std::cout << "\ttransitivity: " << intervals[x]->id+1 << " -> " 
		      << intervals[y]->id+1 << std::endl;
	    if(!successor[x].member(y)) {
	      consistent &= scope[pair_index[x][y]]->setDomain(x>y);
	      successor[x].insert(y);
	      predecessor[y].insert(x);
	    }
	  }
      } else {

	std::cout << "already know that " << intervals[node_i]->id+1 << " -> " 
		  << intervals[node_j]->id+1 << std::endl;


      }
    } else {
      // node_j has to be before node_i
      //successor[node_j].insert(node_i);

      // node_i has to be after node_j
      //predecessor[node_i].insert(node_j);


      if(!successor[node_j].member(node_i)) {

	std::cout << "register " << intervals[node_j]->id+1 << " -> " 
		  << intervals[node_i]->id+1 << std::endl;

	// every successor of node_i (included) 
	// is also successor of 
	// every predecessor of node_j (included)
	for(j=0; consistent && j<predecessor[node_j].size; ++j)
	  for(i=0; consistent && i<successor[node_i].size; ++i) {
	    x = predecessor[node_j][j];
	    y = successor[node_i][i];
	    //successor[x].insert(y);
	    //predecessor[y].insert(x);
	    
	    std::cout << "\ttransitivity: " << intervals[x]->id+1 << " -> " 
		      << intervals[y]->id+1 << std::endl;
	    
	    if(!successor[x].member(y)) {
	      consistent &= scope[pair_index[x][y]]->setDomain(x>y);
	      successor[x].insert(y);
	      predecessor[y].insert(x);
	    }
	    //consistent &= scope[pair_index[x][y]]->setDomain(x>y);
	  }
      } else {

	std::cout << "already know that " << intervals[node_j]->id+1 << " -> " 
		  << intervals[node_i]->id+1 << std::endl;


      }
    }
  }


//   if(e == VALUETRIGGER) {
//     if(scope[changedIdx]->value() == 0) {
//       first = get_first(changedIdx);
//       second = get_second(changedIdx);
//     } else if(scope[changedIdx]->value() == 1) {
//       first = get_second(changedIdx);
//       second = get_first(changedIdx);
//     }

//     // first come before second, hence we insert second
//     // to the successor of first
//     if(!successor[first].member(second)) successor[first].insert(second);

//     // make sure that every 
    

//   }

  return consistent;
}

void ConstraintTDAG::print(std::ostream& o) const 
{
  o << "TDAG ";
  Constraint::print( o );
}


/**********************************************
 * Tree Constraint 
 **********************************************/

ConstraintTree::ConstraintTree(Solver* s, VariableInt** v, const int n)
  : Constraint( s, v, n, VALUETRIGGER )
{
  root = new ReversibleNum<int>[arity];
  leftChild = new ReversibleNum<int>[arity];
  rightChild = new ReversibleNum<int>[arity];
  descendant = new ReversibleSet[arity];
  
  for(int i=0; i<arity; ++i)
    {
      s->binds( root[i] );
      root[i].setValue( i );
      s->binds( leftChild[i] );
      leftChild[i].setValue( -1 );
      s->binds( rightChild[i] );
      rightChild[i].setValue( -1 );

      //root[i].init( s, i );
      //leftChild[i].init( s, -1 );
      //rightChild[i].init( s, -1 );
      s->binds( descendant[i] );
      descendant[i].setValue( 0, arity-1, BitSet::empt );
      descendant[i].insert( i );
    }
}

ConstraintTree::~ConstraintTree()
{
  delete [] root;
  delete [] descendant;
  delete [] leftChild;
  delete [] rightChild;
}

inline int ConstraintTree::check( const int* ) const 
{
  return 0; // TODO
}

inline bool ConstraintTree::propagate()
{
  return true;
}

inline bool ConstraintTree::propagate(const int changedIdx, const int e)
{

//   std::cout << "propagate ";
//   scope[changedIdx]->print( std::cout );
//   std::cout << std::endl;
//   for(int x=0; x<arity; ++x)
//     {
//       scope[x]->print( std::cout );
//       std::cout << " root: " << (root[x]) << "  descendants: ";
//       descendant[x].print( std::cout );
//       std::cout << " left: " << (leftChild[x]) << " right: " << (rightChild[x]) << std::endl;
//     }
//   std::cout << std::endl;
  

  int i, j = scope[changedIdx]->min(), consistent = true;
  // rule 1: only one root, its parent is itself
  if( j == changedIdx )
    {
      consistent = (root[changedIdx] == changedIdx);
      for(i=0; consistent && i<arity; ++i)
	consistent &= ( i == changedIdx || scope[i]->remove( i ) );
    }
  else 
    {
 //      int lc = leftChild[j];
//       int rc = rightChild[j];
//       // binary tree:
//       if( lc < 0 )
// 	leftChild[j] = changedIdx;
//       else if( rc < 0 ) {
// 	rightChild[j] = changedIdx;
// 	for(i=0; consistent && i<arity; ++i) 
// 	  consistent &= ( i==lc || i==changedIdx || i==j || scope[i]->remove( j ) );
//       } else consistent = false;

      for(i=0; i<arity; ++i)
	if(descendant[i].member(j))
	  descendant[i].insert(changedIdx);

      if( consistent ) {
	// point 2:
	int r=root[j];
	
	// rule 3:
	descendant[r].unionWith( descendant[changedIdx] );
	
	// rule 4:
	consistent = scope[r]->removeSet( descendant[changedIdx] );
	
	// rule 5:
	i = descendant[changedIdx].min();
	do { 
	  root[i] = r;
	  i = descendant[changedIdx].next( i );
	} while( i != NOVAL );	
      }
    }
//   for(int x=0; x<arity; ++x)
//     {
//       scope[x]->print( std::cout );
//       std::cout << " root: " << (root[x]) << "  descendants: ";
//       descendant[x].print( std::cout );
//       std::cout << " left: " << (leftChild[x]) << " right: " << (rightChild[x]) << std::endl;
//     }
//   std::cout << consistent << std::endl << std::endl;

  return consistent;
}

void ConstraintTree::print(std::ostream& o) const 
{
  o << "Tree ";
  Constraint::print( o );
}



/**********************************************
 * Explicit Tree Constraint 
 **********************************************/

ConstraintTreeExplicit::ConstraintTreeExplicit(Solver* s, VariableInt** v, const int n)
  : Constraint( s, v, n+(n*n), VALUETRIGGER )
{
  
  //std::cout << n << std::endl;
  descendant = new VariableInt*[n*n];
  memcpy(descendant, v+n, n*n*sizeof(VariableInt*));
  N = n;
//   //(scope[n]);

//   std::cout << "here" << std::endl;

//   for(int i=0; i<(n*n); ++i) {
//     descendant[i]->print(std::cout);
//     std::cout << std::endl;
//   }

  root = new ReversibleNum<int>[N];
//   leftChild = new ReversibleNum<int>[N];
//   rightChild = new ReversibleNum<int>[N];
  
  for(int i=0; i<N; ++i)
    {
      s->binds( root[i] );
      root[i].setValue( i );
 //      s->binds( leftChild[i] );
//       leftChild[i].setValue( -1 );
//       s->binds( rightChild[i] );
//       rightChild[i].setValue( -1 );
    }
}

ConstraintTreeExplicit::~ConstraintTreeExplicit()
{
  delete [] root;
  delete [] descendant;
  //delete [] leftChild;
  //delete [] rightChild;
}

inline int ConstraintTreeExplicit::check( const int* ) const 
{
  return 0; // TODO
}

inline bool ConstraintTreeExplicit::propagate()
{
//   std::cout << "propagate ";
//   std::cout << std::endl;
//   for(int x=0; x<N; ++x)
//     {
//       scope[x]->print( std::cout );
//       std::cout << " root: " << (root[x]) << "  descendants: {";
//       for(int y=0; y<N-1; ++y)
// 	if(descendant[x*N+y]->min())
// 	  std::cout << y << " ";
//       if(descendant[(x+1)*N-1]->min())
// 	std::cout << (N-1) ;
//       std::cout << "} " << std::endl;
//     }
//   std::cout << std::endl;

  bool consistent = true;
  for(int i=0; consistent && i<N; ++i)
    consistent = descendant[i*N+i]->setDomain(1);

//   for(int x=0; x<N; ++x)
//     {
//       scope[x]->print( std::cout );
//       std::cout << " root: " << (root[x]) << "  descendants: {";
//       for(int y=0; y<N-1; ++y)
// 	if(descendant[x*N+y]->min())
// 	  std::cout << y << " ";
//       if(descendant[(x+1)*N-1]->min())
// 	std::cout << (N-1) ;
//       std::cout << "} left: " << (leftChild[x]) << " right: " << (rightChild[x]) << std::endl;
//     }
//   std::cout << consistent << std::endl << std::endl;
  return consistent;
}

inline bool ConstraintTreeExplicit::propagate(const int changedIdx, const int e)
{

//   std::cout << "propagate ";
//   scope[changedIdx]->print( std::cout );
//   std::cout << std::endl;
//   for(int x=0; x<N; ++x)
//     {
//       scope[x]->print( std::cout );
//       std::cout << " root: " << (root[x]) << "  descendants: {";
//       for(int y=0; y<N-1; ++y)
// 	if(descendant[x*N+y]->min())
// 	  std::cout << y << " ";
//       if(descendant[(x+1)*N-1]->min())
// 	std::cout << (N-1) ;
//       std::cout << "} "<< std::endl;
//     }
//   std::cout << std::endl;

  

//  for(int x=0; x<N; ++x)
//     {
//       std::cout << "descendant of " << x << ": " ;
//       for(int y=0; y<N; ++y)
// 	if(descendant[x*N+y]->min())
// 	  std::cout << y << " ";
//       std::cout << std::endl;
//     }
//   std::cout << std::endl;  


  int i, j = scope[changedIdx]->min(), consistent = true;

  if( changedIdx < N ) {
    //rule 1: only one root, its parent is itself
    if( j == changedIdx )
      {
// 	consistent = (root[changedIdx] == changedIdx);
// 	for(i=0; consistent && i<N; ++i)
// 	  consistent &= ( i == changedIdx || scope[i]->remove( i ) );
      }
    else 
      {
	for(i=0; i<N; ++i) 
	  if(descendant[i*N+j]->min())
	    consistent = descendant[i*N+changedIdx]->setDomain(1);
	
	if( consistent ) {
	  //point 2:
	  int r=root[j];
	  
	  //rule 3:
	  for(i=0; consistent && i<N; ++i)
	  consistent = descendant[r*N+i]->setMin( descendant[changedIdx*N+i]->min() );
	
	  //rule 4:
	  for(i=0; consistent && i<N; ++i)
	    consistent = ( !descendant[changedIdx*N+i]->min() ||
			   scope[r]->remove( i ) );
	  
	  //rule 5:
	  for(i=0; i<N; ++i)
	    if( descendant[changedIdx*N+i]->min() )
	      root[i] = r;
	}
      }
//     for(int x=0; x<N; ++x)
//       {
// 	scope[x]->print( std::cout );
// 	std::cout << " root: " << (root[x]) << "  descendants: {";
// 	for(int y=0; y<N-1; ++y)
// 	  if(descendant[x*N+y]->min())
// 	    std::cout << y << " ";
// 	if(descendant[(x+1)*N-1]->min())
// 	  std::cout << (N-1) ;
// 	std::cout << "} "<< std::endl;
//       }
//     std::cout << consistent << std::endl << std::endl;
  }

//   std::cout << std::endl;
//  for(int x=0; x<N; ++x)
//     {
//       std::cout << "descendant of " << x << ": " ;
//       for(int y=0; y<N; ++y)
// 	if(descendant[x*N+y]->min())
// 	  std::cout << y << " ";
//       std::cout << " <= D <= ";
//       for(int y=0; y<N; ++y)
// 	if(descendant[x*N+y]->max())
// 	  std::cout << y << " ";
//       std::cout << std::endl;
//     }
//  std::cout << std::endl << std::endl;  

  return consistent;
}

void ConstraintTreeExplicit::print(std::ostream& o) const 
{
  o << "TreeExplicit ";
  Constraint::print( o );
}


// /**********************************************
//  * NotAllEqual Constraint
//  **********************************************/ 

// ConstraintNotAllEqual::ConstraintNotAllEqual(Solver*, VariableInt** v, const int n)
// {

// }

// ConstraintNotAllEqual::~ConstraintNotAllEqual()
// {

// }

// inline int ConstraintNotAllEqual::check( const int* ) const 
// {

// }

// inline bool ConstraintNotAllEqual::propagate()
// {

// }

// inline bool ConstraintNotAllEqual::propagate(const int changedIdx, const int e)
// {

// }

// void ConstraintNotAllEqual::print(std::ostream& o) const
// {
//   for(int i=0; i<arity-1; ++i)
// }


/**********************************************
 * AllDiff Constraint 
 **********************************************/

const int INCONSISTENT = 0;
const int CHANGES      = 1;
const int NO_CHANGES   = 2;

ConstraintAllDiff::ConstraintAllDiff(Solver *s, VariableInt** v, const int n) 
  : Constraint(s, v, n, RANGETRIGGER) 
{
  delayed = true;
  init();
}

void ConstraintAllDiff::init() 
{
  int i;
  lastLevel = -1;
  nb = 0;

  iv        = new Interval[arity];
  minsorted = new Interval*[arity];
  maxsorted = new Interval*[arity];
  bounds    = new int[2*arity+2];
  std::fill(bounds, bounds+2*arity+2, 0);

  for( i=0; i<arity; ++i ) {
    minsorted[i] = maxsorted[i] = &iv[i];  
    iv[i].min = iv[i].max = NOVAL;
  }

  t = new int[2*arity+2];
  d = new int[2*arity+2];
  h = new int[2*arity+2];
}


ConstraintAllDiff::~ConstraintAllDiff() 
{
  delete [] bounds;
  delete [] maxsorted;
  delete [] minsorted;
  delete [] iv;
  delete [] h;
  delete [] d;
  delete [] t;
}

void sortmin( Interval *v[], int n ) 
{
  int i, current;
  bool sorted;
  Interval *t;

  current = n-1;
  sorted = false;
  while( !sorted ) {
    sorted = true;
    for( i = 0; i < current; i++ ) {
      if( v[i]->min > v[i+1]->min ) {
        t = v[i];
        v[i] = v[i+1];
        v[i+1] = t;
        sorted = false;
      }
    }
    current--;
  }
}

void sortmax( Interval *v[], int n ) 
{
  int i, current;
  bool sorted;
  Interval *t;

  current = 0;
  sorted = false;
  while( !sorted ) {
    sorted = true;
    for( i = n-1; i > current; i-- ) {
      if( v[i]->max < v[i-1]->max ) {
        t = v[i];
        v[i] = v[i-1];
        v[i-1] = t;
        sorted = false;
      }
    }
    current++;
  }
}

void ConstraintAllDiff::sortit() 
{  
  int i,j,nb,min,max,last;

  sortmin(minsorted, arity);
  sortmax(maxsorted, arity);

  min = minsorted[0]->min;
  max = maxsorted[0]->max + 1;
  bounds[0] = last = min-2;

  for (i=j=nb=0;;) { // merge minsorted[] and maxsorted[] into bounds[]
    if (i<arity && min<=max) {	// make sure minsorted exhausted first
      if (min != last)
        bounds[++nb] = last = min;
      minsorted[i]->minrank = nb;
      if (++i < arity)
        min = minsorted[i]->min;
    } else {
      if (max != last)
	bounds[++nb] = last = max;
      maxsorted[j]->maxrank = nb;
      if (++j == arity) break;
      max = maxsorted[j]->max + 1;
    }
  }
  ConstraintAllDiff::nb = nb;
  bounds[nb+1] = bounds[nb] + 2;
}


void pathset(int *t, int start, int end, int to) 
{
  int k, l;
  for (l=start; (k=l) != end; t[k]=to) 
    l = t[k];  
}

int pathmin(int *t, int i) 
{
  for (; t[i] < i; i=t[i]) ;
  return i;
}

int pathmax(int *t, int i) 
{
  for (; t[i] > i; i=t[i]) ;  
  return i;
}


int ConstraintAllDiff::filterlower() 
{
  int i,j,w,x,y,z;
  int changes = 0;

  for (i=1; i<=nb+1; i++)
    d[i] = bounds[i] - bounds[t[i]=h[i]=i-1];
  for (i=0; i<arity; i++) { // visit Intervals in increasing max order
    x = maxsorted[i]->minrank; y = maxsorted[i]->maxrank;
    j = t[z = pathmax(t, x+1)];
    if (--d[z] == 0)
      t[z = pathmax(t, t[z]=z+1)] = j;
    pathset(t, x+1, z, z); // path compression
    if (d[z] < bounds[z]-bounds[y]) return INCONSISTENT; // no solution
    if (h[x] > x) {
      maxsorted[i]->min = bounds[w = pathmax(h, h[x])];
      pathset(h, x, w, w); // path compression
      changes = 1;
    }
    if (d[z] == bounds[z]-bounds[y]) {
      pathset(h, h[y], j-1, y); // mark hall Interval
      h[y] = j-1; //("hall Interval [%d,%d)\n",bounds[j],bounds[y]);
    }
  }
  if( changes )
    return CHANGES;
  else
    return NO_CHANGES;
}


int ConstraintAllDiff::filterupper()
{
  int i,j,w,x,y,z;
  int changes = 0;

  for (i=0; i<=nb; i++)
    d[i] = bounds[t[i]=h[i]=i+1] - bounds[i];
  for (i=arity; --i>=0; ) { // visit Intervals in decreasing min order
    x = minsorted[i]->maxrank; y = minsorted[i]->minrank;
    j = t[z = pathmin(t, x-1)];
    if (--d[z] == 0)
      t[z = pathmin(t, t[z]=z-1)] = j;
    pathset(t, x-1, z, z);
    if (d[z] < bounds[y]-bounds[z]) return INCONSISTENT; // no solution
    if (h[x] < x) {
      minsorted[i]->max = bounds[w = pathmin(h, h[x])] - 1;
      pathset(h, x, w, w);
      changes = 1;
    }
    if (d[z] == bounds[y]-bounds[z]) {
      pathset(h, h[y], j+1, y);
      h[y] = j+1;
    }
  }
  if( changes )
    return CHANGES;
  else
    return NO_CHANGES;
}


bool ConstraintAllDiff::propagate(const int changedIdx, const int e) 
{
  return propagate();
}

bool ConstraintAllDiff::propagate() 
{
  int i, a, b;

  int status_lower, status_upper;
  int l, u;

  a = 0;
  b = arity;

  //if( lastLevel != ((solver->level) - 1) ) {
  if( lastLevel != ((*level) - 1) ) {
    // not incremental
    status_lower = CHANGES;
    status_upper = CHANGES;
    i = 0;
    while (i < arity) {
      iv[i].min = scope[i]->min();
      iv[i].max = scope[i]->max();
      i++;
    }
  }
  else {
    // incremental
    status_lower = NO_CHANGES;
    status_upper = NO_CHANGES;
    for( i = a; i < b; i++ ) {
      l = iv[i].min;
      u = iv[i].max;
      iv[i].min = scope[i]->min();
      iv[i].max = scope[i]->max();
      if( l != iv[i].min ) status_lower = CHANGES;
      if( u != iv[i].max ) status_upper = CHANGES;
    }
  }

  lastLevel = *level;//(solver->level);

  if( status_lower == NO_CHANGES && status_upper == NO_CHANGES ) 
    return true;

  sortit();

  status_lower = filterlower();
  if( status_lower != INCONSISTENT )
    status_upper = filterupper();  

  if( (status_lower == INCONSISTENT) || (status_upper == INCONSISTENT) ) 
    { return false; }
  else
    if( (status_lower == CHANGES) || (status_upper == CHANGES) ) {
      i = 0;
      while (i < arity) {
	if( !scope[i]->setMin( iv[i].min ) )  { return false; }
	if( !scope[i]->setMax( iv[i].max ) )  { return false; }
	i++;
      }
    }  
  return true;
}

int ConstraintAllDiff::check( const int* s ) const 
{
  int i=arity, j;
  while(--i) {
    j=i;
    while(j--)
      if( s[i] == s[j] ) return 1;
  }
  return 0; 
}

void ConstraintAllDiff::print(std::ostream& o) const 
{    
  o << "alldiff(";
  for(int i=0; i<arity-1; ++i) {
    scope[i]->printshort( o );
    o << ", ";
  }
  scope[arity-1]->printshort( o );
  o << ")";
}

/**********************************************
 * EdgeFinder Constraint
 **********************************************/ 

bool ConstraintEdgeFinder::jacksonPSchedule()
{  
  int residue[arity], rest=arity;
  int min_rd = NOVAL, max_dd = 0;
  int i=arity, j, x, y, tmp;
  while( i-- ) {
    sorted[i] = i;
    if( min_rd > scope[i]->min() ) min_rd = scope[i]->min();
    if( max_dd < (scope[i]->max() + durations[i]) ) 
      max_dd = (scope[i]->max() + durations[i]);
  }
  int jps[max_dd];
  //int jps_length = 0;

  // Bubble Sort
  for( x=0; x<arity; x++ ) {
    for( y=0; y<arity-1-x; y++ ) {
      if( scope[sorted[y  ]]->max()+durations[sorted[y  ]] < 
	  scope[sorted[y+1]]->max()+durations[sorted[y+1]] ) {
	tmp         = sorted[y+1];
	sorted[y+1] = sorted  [y];
	sorted  [y] =         tmp;
      }
    }
    // initialisation of data-structures
    residue[arity-1-x] = durations[sorted[arity-1-x]];
    Pk[x][0] = -1; Pk[x][arity] = 0;
  }

  int dur, cur = min_rd;
  while( rest ) { // there is still something to schedule
    i = arity;
    // compute Pk
    while( i--) 
      if( Pk[i][0] < 0 && scope[i]->min() <= cur ) 
	for(j=0; j<arity; ++j) {
	  Pk[i][sorted[j]] = residue[j];
	  Pk[i][arity] += residue[j];	
	}
    // jps
    i = arity; dur = NOVAL; tmp = NOVAL;
    while( tmp == NOVAL && i-- ) 
      if( residue[i] ) {
	if( scope[sorted[i]]->min() <= cur ) {
	  // this is the task to schedule next in the jps.
	  if( residue[i] < dur ) dur = residue[i];
	  tmp = sorted[i];
	  if( !(residue[i] -= dur)) --rest;	  
	} else 
	  // this activity has an earlier due date but its release date is later.
	  // we change dur so as not to miss it
	  if( scope[sorted[i]]->min() - cur < dur ) dur = (scope[sorted[i]]->min() - cur);
      }
    if( cur+dur <= max_dd )
      std::fill( jps+cur, jps+cur+dur, tmp );
    cur += dur;
  }
 
  return (cur <= max_dd);
}

ConstraintEdgeFinder::ConstraintEdgeFinder(Solver *s, VariableInt** v, 
					   const int n, const int *d)
  : Constraint(s, v, n, RANGETRIGGER)
{
  sorted = new int[arity];
  durations = new int[arity];
  Pk = new int*[arity];
  for(int i=0; i<arity; ++i) {
    durations[i] = d[i];
    Pk[i] = new int[arity+1];  
  }
  //int i=arity;
  //while( i-- )
  //  triggers[i] = scope[i]->triggerOnRange();//&(scope[i]->rangeTrigger);
}

ConstraintEdgeFinder::~ConstraintEdgeFinder()
{
  delete [] sorted;
  delete [] durations;
  for(int i=0; i<arity; ++i) 
    delete [] Pk[i];  
  delete [] Pk;
}

int ConstraintEdgeFinder::check( const int* s ) const 
{
  // todo
  return 0;
}

bool ConstraintEdgeFinder::propagate(const int changedIdx, const int e)
{
  return propagate();
}

bool ConstraintEdgeFinder::propagate()
{
  if( !jacksonPSchedule() ) return false;

  int i, j, k;
  for( i=0; i<arity; ++i ) {
    Pk[i][arity] += scope[i]->min();
    for( j=0; j<arity; ++j ) {
      k = sorted[j];
      if ( i != k && Pk[i][k] ) {
	if(  Pk[i][arity] > scope[k]->max()+durations[k] ) {
	  if( !scope[i]->setMin(scope[k]->min() + durations[k]) ) return false;
	  if( !scope[k]->setMax(scope[i]->max() - durations[k]) ) return false;
	}
	Pk[i][arity] -= Pk[i][k];
      }
    }
  }
  return true;
}

void ConstraintEdgeFinder::print(std::ostream& o) const
{
  Constraint::print( o );
  o << std::endl;
  for(int i=0; i<arity-1; ++i) {
    o << "[" << scope[i]->min() << "..|" << durations[i]
      << "|.." << scope[i]->max()+durations[i] << "]," ;
  }
  o << "[" << scope[arity-1]->min() << "..|" << durations[arity-1]
    << "|.." << scope[arity-1]->max()+durations[arity-1] << "]" ;
}


/**********************************************
 * Boolean Sum Equality Constraint  
 **********************************************/

ConstraintBoolSumEqual::ConstraintBoolSumEqual(Solver *s, 
					       VariableInt **v, 
					       const int n,
					       const int t)
  : Constraint(s, v, n, VALUETRIGGER)
{
  total = t;
  s->binds( min_ );
  s->binds( max_ );
  min_.setValue(0);
  max_.setValue(n);
  s->binds( unassigned );
  unassigned.setValue(0, n-1, n);
}

ConstraintBoolSumEqual::~ConstraintBoolSumEqual()
{
}

int ConstraintBoolSumEqual::check( const int* s ) const 
{
  int i = arity, tot = 0;
   while( i-- )
    tot += s[i];

  return ( tot > total ? tot - total : total - tot );
}

bool ConstraintBoolSumEqual::propagate()
{
  int _min_ = 0;
  int _max_ = 0;

  for( int i=0; i<arity; ++i ) {
    _min_ += scope[i]->min();
    _max_ += scope[i]->max();
  }
  if(_min_ > total || _max_ < total) return false;
  else if(_min_ == total ) {
    for( int i=0; i<arity; ++i ) 
      if( !scope[i]->isGround() ) scope[i]->setDomain(0);
  } else if(_max_ == total ) {
    for( int i=0; i<arity; ++i ) 
      if( !scope[i]->isGround() ) scope[i]->setDomain(1);
  }
  return true;
} 
 

bool ConstraintBoolSumEqual::propagate(const int changedIdx, const int event) 
{
  int i, bound[2], polarity=scope[changedIdx]->min();

  if( polarity ) ++min_;
  else --max_;
  
  bound[0] = min_;
  bound[1] = max_;

  if( bound[polarity] == total ) {
    return true;
  } else if( bound[1-polarity] == total ) {

    bound[0] = 0;
    bound[1] = arity;
    
    if( polarity ) {
      for(i=0; i<arity; ++i) {
	if( !scope[i]->min() ) {
	  if( --bound[1] < total ) return false;
	  scope[i]->setDomain(0);
	} else if( ++bound[0] > total ) return false;
      }
    } else {
      for(i=0; i<arity; ++i) {
	if( scope[i]->max() ) {
	  if( ++bound[0] > total ) return false;
	  scope[i]->setDomain(1);
	} else if( --bound[1] < total ) return false;
      }
    }
  }

  return true;
}



// bool ConstraintBoolSumEqual::propagate(const int changedIdx, const int event) 
// {
//   int i, j, bound[2], polarity=scope[changedIdx]->min();
//   unassigned.erase( changedIdx );
//   if( polarity ) ++min_;
  
//   bound[0] = min_;
//   bound[1] = (min_ + unassigned.size);

//   if( bound[polarity] == total ) {
//     return true;
//   } else if( bound[1-polarity] == total ) {
//     i = unassigned.size;
//     if( polarity ) {
//       while( i-- ) {
// 	j = unassigned[i];      
// 	if( !scope[j]->min() ) {
// 	  if( --bound[1] < total ) return false;
// 	  scope[j]->setDomain(0);
// 	} else if( ++bound[0] > total ) return false;
//       }
//     } else {
//       while( i-- ) {
// 	j = unassigned[i];
// 	if( scope[j]->max() ) {
// 	  if( ++bound[0] > total ) return false;
// 	  scope[j]->setDomain(1);
// 	} else if( --bound[1] < total ) return false;
//       }
//     }
//   }

//   return true;
// }

void ConstraintBoolSumEqual::print(std::ostream& o) const
{ 
  //std::cout << "BOOLSUMEQ" << std::endl;

  o<< "(" ;
  for(int i=0; i<arity-1; ++i) {
    scope[i]->printshort( o );
    o << " + ";
  }
  scope[arity-1]->printshort( o );
  o << ") = " << total ;
}


/**********************************************
 * Boolean Sum Less Constraint  
 **********************************************/

ConstraintBoolSumLess::ConstraintBoolSumLess(Solver *s, 
					     VariableInt **v, 
					     const int n,
					     const int t)
  : Constraint(s, v, n, VALUETRIGGER)
{
  total = t;
  s->binds( min_ );
  min_.setValue(0);
  //min_.init( s );
}

ConstraintBoolSumLess::~ConstraintBoolSumLess()
{
}

int ConstraintBoolSumLess::check( const int* s ) const 
{
  int i = arity, tot = 0;
   while( i-- )
    tot += s[i];

  return ( tot > total ? tot - total : 0 );
}

bool ConstraintBoolSumLess::propagate()
{
  int _min_ = 0;
  
  for( int i=0; i<arity; ++i ) 
    _min_ += scope[i]->min();
  if(_min_ > total) return false;
  else if(_min_ == total ) 
    for( int i=0; i<arity; ++i ) 
      if( !scope[i]->isGround() ) scope[i]->setDomain(0);
  return true;

//   bool consistent = true;
//   for( int i=0; consistent && i<arity; ++i )
//     if( scope[i]->isGround() )
//       consistent = propagate( i, VALUETRIGGER );
//   return consistent;
} 

bool ConstraintBoolSumLess::propagate(const int changedIdx, const int event) 
{
  int i, smin;
  if( scope[changedIdx]->min() ) {
    ++min_;
    if( min_ == total ) {
      smin = 0;
      for(i=0; i<arity; ++i) {
	if( !scope[i]->min() ) {
	  scope[i]->setDomain(0);
	} else if( ++smin > total ) return false;
      }
    }
  }
  return true;
}

// bool ConstraintBoolSumLess::propagate(const int changedIdx, const int e) 
// {
//   int i, smin = 0;

//   if( scope[changedIdx]->min() ) ++min_;
//   smin = min_;

//   if( !isLinked() ) {
//     int idx = (watch[0] == scope[changedIdx]);
//     if( watch[idx]->min() ) ++smin;
//   }
  
//   if( total < smin ) return false;

//   if( smin == total ) {
//     i = arity;
//     while( i-- )
//       scope[i]->setMax( scope[i]->min() );
//   }

//   return true;
// }

void ConstraintBoolSumLess::print(std::ostream& o) const
{ 
  o<< "bsl(" ;
  for(int i=0; i<arity-1; ++i) {
    scope[i]->printshort( o );
    o << " + ";
  }
  scope[arity-1]->printshort( o );
  o << ") <= " << total ;
}


/**********************************************
 * Boolean Sum More Constraint  
 **********************************************/

ConstraintBoolSumMore::ConstraintBoolSumMore(Solver *s, 
					     VariableInt **v, 
					     const int n,
					     const int t)
  : Constraint(s, v, n, VALUETRIGGER)
{
  total = t;
  s->binds( max_ );
  max_.setValue(arity);
  //max_.init( s );
}

ConstraintBoolSumMore::~ConstraintBoolSumMore()
{
}

int ConstraintBoolSumMore::check( const int* s ) const 
{
  int i = arity, tot = 0;
   while( i-- )
    tot += s[i];

  return ( tot < total ? tot - total : 0 );
}

bool ConstraintBoolSumMore::propagate()
{
  int _max_ = 0;
  
  for( int i=0; i<arity; ++i ) 
    _max_ += scope[i]->max();
  if(_max_ < total) return false;
  else if(_max_ == total ) 
    for( int i=0; i<arity; ++i ) 
      if( !scope[i]->isGround() ) scope[i]->setDomain(1);
  return true;

//   bool consistent = true;
//   for( int i=0; consistent && i<arity; ++i )
//     if( scope[i]->isGround() )
//       consistent = propagate( i, VALUETRIGGER );
//   return consistent;
} 

bool ConstraintBoolSumMore::propagate(const int changedIdx, const int event) 
{
  int i, smax;
  if( scope[changedIdx]->equal(0) ) {
    --max_;
    if( max_ == total ) {
      smax = arity;   
      for(i=0; i<arity; ++i) {
	if( scope[i]->max() ) {
	  scope[i]->setDomain(1);
	} else if( --smax < total ) return false;
      }
    }
  }
  return true;
}

// bool ConstraintBoolSumMore::propagate(const int changedIdx, const int e) 
// {
//   int i, smax = 0;

//   if( !scope[changedIdx]->max() ) --max_;
//   smax = max_;

//   if( !isLinked() ) {
//     int idx = (watch[0] == scope[changedIdx]);
//     if( !watch[idx]->max() ) --smax;
//   }
  
//   if( smax < total ) return false;

//   if( smax == total ) {
//     i = arity;
//     while( i-- )
//       scope[i]->setMin( scope[i]->max() );
//   }

//   return true;
// }

void ConstraintBoolSumMore::print(std::ostream& o) const
{ 
  o<< "(" ;
  for(int i=0; i<arity-1; ++i) {
    scope[i]->printshort( o );
    o << " + ";
  }
  scope[arity-1]->printshort( o );
  o << ") >= " << total ;
}


// /**********************************************
//  * SimpleSum Predicate 
//  **********************************************/

// PredicateSimpleSum::PredicateSimpleSum(Solver *s, VariableInt** v, const int n)
//   : Constraint(s, v, n, RANGETRIGGER)
// {
//   //  delayed = true;
// }

// PredicateSimpleSum::~PredicateSimpleSum()
// {
// }

// int PredicateSimpleSum::check( const int* s ) const 
// {
//   int i = arity-1, total = 0;
//    while( i-- )
//     total += s[i];

//   return ( total > s[arity-1] ? total - s[arity-1] : s[arity-1] - total );
// }

// bool PredicateSimpleSum::propagate(const int changedIdx, const int e) 
// {
//   return propagate();
// }

// bool PredicateSimpleSum::propagate() 
// {

//   int i, n = arity-1;
//   // compute the max and th min
//   int smin = 0, smax = 0;
//   for(i=0; i<n; ++i) {
//     smax += scope[i]->max();
//     smin += scope[i]->min();
//   }

//   if( !scope[n]->setMax( smax ) || 
//       !scope[n]->setMin( smin ) ) return false;

//   smax -= scope[n]->min();
//   smin -= scope[n]->max();

//   // prune domains
//   int aux;
//   for(i=0; i<n; ++i) {
//     aux = scope[i]->max();
//     if( (!scope[i]->setMax(scope[i]->min() - smin)) 
// 	||
// 	(!scope[i]->setMin(aux - smax)) ) 
//       return false;
//   }

//   return true;
// }

// void PredicateSimpleSum::print(std::ostream& o) const
// { 
//   o<< "(" ;
//   for(int i=0; i<arity-2; ++i) {
//     scope[i]->printshort( o );
//     o << " + ";
//   }
//   scope[arity-2]->printshort( o );
//   o << ") = " ;
//   scope[arity-1]->printshort( o );
// }



// /**********************************************
//  * SumEqual Constraint 
//  **********************************************/

// ConstraintSumEqual::ConstraintSumEqual(Solver *s, VariableInt** v, 
// 				       const int n,
// 				       const int t)
//   : Constraint(s, v, n, RANGETRIGGER)
// {
//   total = t;
// }

// ConstraintSumEqual::~ConstraintSumEqual()
// {
// }

// int ConstraintSumEqual::check( const int* s ) const 
// {
//   int i = arity, tot = 0;
//    while( i-- )
//     tot += s[i];

//   return ( tot > total ? tot - total : total - tot );
// }


// bool ConstraintSumEqual::propagate(const int changedIdx, const int e) 
// {
//   return propagate();
// } 

// bool ConstraintSumEqual::propagate()
// {
//   int i;
//   // compute the max and th min
//   int smin = 0, smax = 0;
//   for(i=0; i<arity; ++i) {
//     smax += scope[i]->max();
//     smin += scope[i]->min();
//   }

//   if( smax < total || smin > total ) return false;

//   smax -= total;
//   smin -= total;

//   // prune domains
//   int aux;
//   for(i=0; i<arity; ++i) {
//     aux = scope[i]->max();
//     if( (!scope[i]->setMax(scope[i]->min() - smin)) 
// 	||
// 	(!scope[i]->setMin(aux - smax)) ) 
//       return false;
//   }

//   return true;
// }

// void ConstraintSumEqual::print(std::ostream& o) const
// { 
//   o<< "(" ;
//   for(int i=0; i<arity-1; ++i) {
//     scope[i]->printshort( o );
//     o << " + ";
//   }
//   scope[arity-1]->printshort( o );
//   o << ") = " << total;
// }


// /**********************************************
//  * WeightedSumEqual Constraint 
//  **********************************************/

// ConstraintWeightedSumEqual::ConstraintWeightedSumEqual(Solver *s, VariableInt** v, 
// 						       const int n,
// 						       const int *w)
//   : Constraint(s, v, n, RANGETRIGGER)
// {
//   weights = new int[n+1];
//   int i;
//   for(i=0; i<=n; ++i) 
//     weights[i] = w[i];
// }

// ConstraintWeightedSumEqual::~ConstraintWeightedSumEqual()
// {
//   delete [] weights;
// }

// int ConstraintWeightedSumEqual::check( const int* s ) const 
// {
//   int i = arity, tot = 0;
//    while( i-- )
//      tot += (weights[i] * s[i]);

//   return ( tot > weights[arity] ? tot - weights[arity] : weights[arity] - tot );
// }

// bool ConstraintWeightedSumEqual::propagate(const int changedIdx, const int e) 
// {
//   return propagate();
// }

// bool ConstraintWeightedSumEqual::propagate() 
// {
//   int i;
//   // compute the max and th min
//   int smin=0, smax=0;
//   for(i=0; i<arity; ++i) {
//     if( weights[i] > 0 ) {
//       smax += weights[i] * scope[i]->max();
//       smin += weights[i] * scope[i]->min();
//     } else {
//       smin += weights[i] * scope[i]->max();
//       smax += weights[i] * scope[i]->min();
//     }
//   }

//   if( smax < weights[arity] || smin > weights[arity] ) return false;

//   smax -= weights[arity];
//   smin -= weights[arity];

//   // prune domains
//   int aux;
//   for(i=0; i<arity; ++i) {
//     if( weights[i] > 0 ) {
//       aux = scope[i]->max();
//       if( (!scope[i]->setMax(scope[i]->min() - smin/weights[i])) 
// 	  ||
// 	  (!scope[i]->setMin(aux - smax/weights[i])) ) 
// 	return false;
//     } else {
//       aux = scope[i]->min();
//       if( (!scope[i]->setMin(scope[i]->max() - smin/weights[i]))
// 	  ||
// 	  (!scope[i]->setMax(aux - smax/weights[i])) )
// 	return false;
//     }
//   }

//   return true;
// }

// void ConstraintWeightedSumEqual::print(std::ostream& o) const
// { 
//   o<< "(" ;
//   if( weights[0] < 0 ) {
//     o << " - ";
//     o.flush() ; 
//   }
//   if( weights[0] != 1 && weights[0] != -1 ) {
//     o << weights[0];
//     o.flush();
//   }
//   scope[0]->printshort( o );
//   o.flush();
//   for(int i=1; i<arity; ++i) {
//     if( weights[i] < 0 ) {
//       o << " - ";
//       o.flush();
//     }
//     else {
//       o << " + ";
//       o.flush();
//     }
//     if( weights[i] != 1 && weights[i] != -1 ) {
//       o << abs(weights[i]);
//       o.flush();
//     }
//     scope[i]->printshort( o );
//     o.flush();
//   }
//   o << ") = " << weights[arity];
// }




/**********************************************
 * WeightedSum Constraint 
 **********************************************/

ConstraintWeightedSum::ConstraintWeightedSum(Solver *s, VariableInt** v, 
					     const int n, 
					     const int ub, const int lb,
					     const int *w,
					     const int wp, const int wn
					     )
  : Constraint(s, v, n, RANGETRIGGER)
{
  weights = new int[n+1];
  int i;
  for(i=0; i<=n; ++i) 
    weights[i] = w[i];
  wpos = wp;
  wneg = wn;
  for(i=wneg; i<arity; ++i) weights[i] = -weights[i];
  lo_bound = new int[n];
  up_bound = new int[n];
  UB = ub;
  LB = lb;
}

ConstraintWeightedSum::~ConstraintWeightedSum()
{
  delete [] weights;
  delete [] lo_bound;
  delete [] up_bound;
}

int ConstraintWeightedSum::check( const int* s ) const 
{
  int i, total=0;
  for(i=0; i<wpos; ++i)
    total += s[i];
  for(i=wpos; i<wneg; ++i)
    total += (weights[i] * s[i]);
  for(i=wneg; i<arity; ++i)
    total -= (weights[i] * s[i]);

  if( LB <= total && total <= UB ) return 0;
  else return ( LB <= total ? total - UB : LB - total );
}


void ConstraintWeightedSum::checkBounds()
{
//   int j;
//   int tot_max, tot_min;

//   tot_max = 0;
//   tot_min = 0;
//   for(j=0; j<arity-1; ++j) {
//     tot_min += (weights[j] * scope[j]->min());
//     tot_max += (weights[j] * scope[j]->max());
//   }
//   tot_min -= (weights[j] * scope[j]->max());
//   tot_max -= (weights[j] * scope[j]->min());

//   if( tot_max < LB || tot_min > UB ) {
//     print(std::cout);
//     std::cout << std::endl << tot_min << " " << LB << " " << UB << " " << tot_max  << std::endl;
//     exit(0);
//   }

//   for(j=0; j<arity; ++j) {
//     tot_min -= (weights[j] * scope[j]->min());
//     tot_min += (weights[j] * scope[j]->max());
    
//     if( UB < tot_min ) {
//       print(std::cout);
//       std::cout << std::endl << tot_min << " " << LB << " " << UB << " " << tot_max << " " << (scope[j]->max()) << std::endl;
//       exit(0);
//     }

//     tot_min += (weights[j] * scope[j]->min());
//     tot_min -= (weights[j] * scope[j]->max());



//     tot_max -= (weights[j] * scope[j]->max());
//     tot_max += (weights[j] * scope[j]->min());

//     if( tot_max < LB ) {
//       print(std::cout);
//       std::cout << std::endl << tot_min << " " << LB << " " << UB << " " << tot_max << " " << (scope[j]->min()) << std::endl;
//       exit(0);
//     }

//     tot_max += (weights[j] * scope[j]->max());
//     tot_max -= (weights[j] * scope[j]->min());
//   }


}

bool ConstraintWeightedSum::propagate(const int changedIdx, const int e) 
{
  return propagate();
}

bool ConstraintWeightedSum::propagate() 
{
  int i;
  // compute the max and th min
  int smin=0, smax=0, maxspan=0, span, consistent = true;
  for(i=0; i<wpos; ++i) {
    smax += (up_bound[i] = scope[i]->max());
    smin += (lo_bound[i] = scope[i]->min());
    span = (up_bound[i]-lo_bound[i]);
    if( span > maxspan ) maxspan = span;
  }
  for(i=wpos; i<wneg; ++i) {
    smax += weights[i] * (up_bound[i] = scope[i]->max());
    smin += weights[i] * (lo_bound[i] = scope[i]->min());
    span = weights[i] * (up_bound[i]-lo_bound[i]);
    if( span > maxspan ) maxspan = span;
  }
  for(i=wneg; i<arity; ++i) {
    smax -= weights[i] * (lo_bound[i] = scope[i]->min());
    smin -= weights[i] * (up_bound[i] = scope[i]->max());
    span = weights[i] * (up_bound[i]-lo_bound[i]);
    if( span > maxspan ) maxspan = span;
  }

  /// check inconsistency
  if( smax < LB || smin > UB ) consistent = false;
  else {
    smax -= LB;
    smin = (UB - smin);
  
    /// prune with respect to the lower bound?
    if( smin < maxspan ) {

      for(i=0; consistent && i<wpos; ++i) 
	if( smin < (up_bound[i] - lo_bound[i]) ) 
	  consistent = scope[i]->setMax( lo_bound[i] + smin );

      for(i=wpos; consistent && i<wneg; ++i) 
	if( smin < (up_bound[i] - lo_bound[i]) * weights[i] ) 
	  consistent = scope[i]->setMax( lo_bound[i] + smin/weights[i] );
      
      for(i=wneg; consistent && i<arity; ++i) 
	if( smin < (up_bound[i] - lo_bound[i]) * weights[i] )
	  consistent = scope[i]->setMin( up_bound[i] - smin/weights[i] );

    }
  
    /// prune with respect to the upwer bound?
    if( smax < maxspan ) {

      for(i=0; consistent && i<wpos; ++i) 
	if( smax < (up_bound[i] - lo_bound[i]) ) 
	  consistent = scope[i]->setMin( up_bound[i] - smax );

      for(i=wpos; consistent && i<wneg; ++i) 
	if( smax < (up_bound[i] - lo_bound[i]) * weights[i] ) 
	  consistent = scope[i]->setMin( up_bound[i] - smax/weights[i] );

      for(i=wneg; consistent && i<arity; ++i) 
	if( smax < (up_bound[i] - lo_bound[i]) * weights[i] ) 	  
	  consistent = scope[i]->setMax( lo_bound[i] + smax/weights[i] );

      }

  }
  
  return consistent;
}

// bool ConstraintWeightedSum::propagate() 
// {
//   int i;
//   // compute the max and th min
//   int smin=0, smax=0, maxspan=0, span, consistent = true, aux;
//   for(i=0; i<wpos; ++i) {
//     smax += (up_bound[i] = scope[i]->max());
//     smin += (lo_bound[i] = scope[i]->min());
//     span = (up_bound[i]-lo_bound[i]);
//     if( span > maxspan ) maxspan = span;
//   }
//   for(i=wpos; i<wneg; ++i) {
//     smax += weights[i] * (up_bound[i] = scope[i]->max());
//     smin += weights[i] * (lo_bound[i] = scope[i]->min());
//     span = weights[i] * (up_bound[i]-lo_bound[i]);
//     if( span > maxspan ) maxspan = span;
//   }
//   for(i=wneg; i<arity; ++i) {
//     smax -= weights[i] * (lo_bound[i] = scope[i]->min());
//     smin -= weights[i] * (up_bound[i] = scope[i]->max());
//     span = weights[i] * (up_bound[i]-lo_bound[i]);
//     if( span > maxspan ) maxspan = span;
//   }

// //   std::cout <<"max: " ;
// //   for(i=0; i<wpos; ++i) {
// //     std::cout << std::setw(4) << up_bound[i] << " ";
// //   }
// //   std::cout << std::endl;

// //   std::cout <<"min: " ;
// //   for(i=0; i<wpos; ++i) {
// //     std::cout << std::setw(4) << lo_bound[i] << " ";
// //   }
// //   std::cout << std::endl;


// //   std::cout << "[" << smin << " " << smax << "]" << std::endl;

//   /// check inconsistency
//   if( smax < LB || smin > UB ) return false;
//   smax -= LB;
//   smin = (UB - smin);
  

// //    std::cout <<  "SMIN, SMAX, MAXSPAN: " << smin << " " << smax << " " << maxspan << std::endl;
// //    print( std::cout );
// //    std::cout << std::endl;


//   /// prune with respect to the lower bound?
//   if( smin < maxspan ) {
//     for(i=0; i<wpos; ++i) {

// //       scope[i]->print( std::cout );
// //       std::cout << " <= " << (lo_bound[i] + smin) << std::endl;

//       if( smin < (up_bound[i] - lo_bound[i]) )
// 	if( !scope[i]->setMax( lo_bound[i] + smin ) ) return false;
//     }
//     for(i=wpos; i<wneg; ++i) 
//       if( smin < (up_bound[i] - lo_bound[i]) * weights[i] ) {
// 	aux = smin/weights[i];
		
// 	if( !scope[i]->setMax( lo_bound[i] + aux ) ) return false;
	
// 	//	assert( smin >= (scope[i]->max() - lo_bound[i]) * weights[i] );
	
//       }
//     for(i=wneg; i<arity; ++i) 
//       if( smin < (up_bound[i] - lo_bound[i]) * weights[i] ) {
// 	aux = smin/weights[i];
	
	
// 	if( !scope[i]->setMin( up_bound[i] - aux ) ) return false;
	
// 	//	assert( smin >= (up_bound[i] - scope[i]->min()) * weights[i] );
	
//       }
//   }      
  
//   /// prune with respect to the upwer bound?
//   if( smax < maxspan ) {
//     for(i=0; i<wpos; ++i) {

// //       scope[i]->print( std::cout );
// //       std::cout << " >= " << (up_bound[i] - smax) << std::endl;

//       if( smax < (up_bound[i] - lo_bound[i]) )
// 	if( !scope[i]->setMin( up_bound[i] - smax ) ) return false;
//     }
//     for(i=wpos; i<wneg; ++i) 
//       if( smax < (up_bound[i] - lo_bound[i]) * weights[i] ) {

// // 	std::cout << std::endl;
// // 	scope[i]->print( std::cout );
// // 	std::cout << std::endl;

// 	aux = smax/weights[i];


// 	//std::cout << smax << " " << weights[i] << " " << aux << " " << (up_bound[i] - aux) << std::endl;

// 	if( !scope[i]->setMin( up_bound[i] - aux ) ) return false;

// // 	if( smax < (up_bound[i] - scope[i]->min()) * weights[i] )  {
// // 	  print( std::cout );
// // 	  std::cout << std::endl;
// // 	  //scope[i]->print( std::cout );
// // 	  //std::cout << std::endl;
// // 	}

// // 	assert( smax >= (up_bound[i] - scope[i]->min()) * weights[i] );

//       }
//     for(i=wneg; i<arity; ++i) 
//       if( smax < (up_bound[i] - lo_bound[i]) * weights[i] ) {
// 	aux = smax/weights[i]; 
// 	if( !scope[i]->setMax( lo_bound[i] + aux ) ) return false;

// 	//	assert( smax >= (scope[i]->max() - lo_bound[i]) * weights[i] );
//       }
//   }

// //   if(arity == 6 && weights[0] > 1000)
// //     checkBounds();
  
//   return consistent;
// }

void ConstraintWeightedSum::print(std::ostream& o) const
{ 
  int i;
  o<< "(" ;
  for(i=0; i<wpos; ++i) {
    if(i) o << " + ";
    scope[i]->printshort( o );
  }
  for(i=wpos; i<wneg; ++i) {
    if(i) o << " + " ;
    o << weights[i];
    scope[i]->printshort( o );
  }
  for(i=wneg; i<arity; ++i) {
    if(i) o << " - " ;
    if(weights[i] != 1) o << weights[i];
    scope[i]->printshort( o );
  }
  if( LB == UB )
    o << ") = " << UB ;
  else if( LB == -NOVAL/2 )
    o << ") < " << (UB+1) ;
  else if( UB == NOVAL/2 )
    o << ") > " << (LB-1) ;
  else 
    o << ") in [" << LB << ".." << UB << "]" ;

//   std::cout << std::endl;
//   for(i=0; i<arity; ++i) {
//     scope[i]->print( std::cout );
//     std::cout << std::endl;
//   }
}

void ConstraintWeightedSum::printFull(std::ostream& o) const
{ 
  int i;
  o<< "(" ;
  for(i=0; i<wpos; ++i) {
    if(i) o << " + ";
    scope[i]->print( o );
  }
  for(i=wpos; i<wneg; ++i) {
    if(i) o << " + " ;
    o << weights[i];
    scope[i]->print( o );
  }
  for(i=wneg; i<arity; ++i) {
    if(i) o << " - " ;
    o << weights[i];
    scope[i]->print( o );
  }
  o << ") " ;
  //o << ") " << = " << weights[arity];
  if( LB == UB )
    o << "= " << UB ;
  else if( LB == -NOVAL/2 )
    o << "< " << (UB+1) ;
  else if( UB == NOVAL/2 )
    o << "> " << (LB-1) ;
  else 
    o << " in [" << LB << ".." << UB << "]" ;
}



/**********************************************
 * LinearCoef Constraint 
 **********************************************/

ConstraintLinearCoef::ConstraintLinearCoef(Solver *s, VariableInt** v, 
					   const int n, 
					   const double ub, const double lb,
					   const double *w,
					   const int wp, const int wn
					   )
  : Constraint(s, v, n, RANGETRIGGER)
{
  weights = new double[n+1];
  int i;
  for(i=0; i<=n; ++i) 
    weights[i] = w[i];
  wpos = wp;
  wneg = wn;
  for(i=wneg; i<arity; ++i) weights[i] = -weights[i];
  lo_bound = new double[n];
  up_bound = new double[n];
  UB = ub;
  LB = lb;
}

ConstraintLinearCoef::~ConstraintLinearCoef()
{
  delete [] weights;
  delete [] lo_bound;
  delete [] up_bound;
}

int ConstraintLinearCoef::check( const int* s ) const 
{
  int i, total=0;
  for(i=0; i<wpos; ++i)
    total += s[i];
  for(i=wpos; i<wneg; ++i)
    total += (weights[i] * s[i]);
  for(i=wneg; i<arity; ++i)
    total -= (weights[i] * s[i]);

  if( LB <= total && total <= UB ) return 0;
  else return ( LB <= total ? total - UB : LB - total );
}


void ConstraintLinearCoef::checkBounds()
{
}

bool ConstraintLinearCoef::propagate(const int changedIdx, const int e) 
{
  return propagate();
}

bool ConstraintLinearCoef::propagate() 
{
  int i;
  // compute the max and th min
  double smin=0, smax=0, maxspan=0, span;
  int consistent = true;

  for(i=0; i<wpos; ++i) {
    smax += (up_bound[i] = (double)(scope[i]->max()));
    smin += (lo_bound[i] = (double)(scope[i]->min()));
    span = (up_bound[i]-lo_bound[i]);
    if( span > maxspan ) maxspan = span;
  }
  for(i=wpos; i<wneg; ++i) {
    smax += weights[i] * (up_bound[i] = (double)(scope[i]->max()));
    smin += weights[i] * (lo_bound[i] = (double)(scope[i]->min()));
    span = weights[i] * (up_bound[i]-lo_bound[i]);
    if( span > maxspan ) maxspan = span;
  }
  for(i=wneg; i<arity; ++i) {
    smax -= weights[i] * (lo_bound[i] = (double)(scope[i]->min()));
    smin -= weights[i] * (up_bound[i] = (double)(scope[i]->max()));
    span = weights[i] * (up_bound[i]-lo_bound[i]);
    if( span > maxspan ) maxspan = span;
  }

  /// check inconsistency
  if( smax < LB || smin > UB ) consistent = false;
  else {
    smax -= LB;
    smin = (UB - smin);
  
    /// prune with respect to the lower bound?
    if( smin < maxspan ) {

      for(i=0; consistent && i<wpos; ++i) 
	if( smin < (up_bound[i] - lo_bound[i]) ) 
	  consistent = scope[i]->setMax( (int)(lo_bound[i] + smin) );

      for(i=wpos; consistent && i<wneg; ++i) 
	if( smin < (up_bound[i] - lo_bound[i]) * weights[i] ) 
	  consistent = scope[i]->setMax( (int)(lo_bound[i] + smin/weights[i]) );
      
      for(i=wneg; consistent && i<arity; ++i) 
	if( smin < (up_bound[i] - lo_bound[i]) * weights[i] )
	  consistent = scope[i]->setMin( (int)(up_bound[i] - smin/weights[i]) );

    }
  
    /// prune with respect to the upwer bound?
    if( smax < maxspan ) {

      for(i=0; consistent && i<wpos; ++i) 
	if( smax < (up_bound[i] - lo_bound[i]) ) 
	  consistent = scope[i]->setMin( (int)(up_bound[i] - smax) );

      for(i=wpos; consistent && i<wneg; ++i) 
	if( smax < (up_bound[i] - lo_bound[i]) * weights[i] ) 
	  consistent = scope[i]->setMin( (int)(up_bound[i] - smax/weights[i]) );

      for(i=wneg; consistent && i<arity; ++i) 
	if( smax < (up_bound[i] - lo_bound[i]) * weights[i] ) 	  
	  consistent = scope[i]->setMax( (int)(lo_bound[i] + smax/weights[i]) );

      }

  }
  
  return consistent;
}

void ConstraintLinearCoef::print(std::ostream& o) const
{ 
  int i;
  o<< "(" ;
  for(i=0; i<wpos; ++i) {
    if(i) o << " + ";
    scope[i]->printshort( o );
  }
  for(i=wpos; i<wneg; ++i) {
    if(i) o << " + " ;
    o << weights[i];
    scope[i]->printshort( o );
  }
  for(i=wneg; i<arity; ++i) {
    if(i) o << " - " ;
    if(weights[i] != 1) o << weights[i];
    scope[i]->printshort( o );
  }
  if( LB == UB )
    o << ") = " << UB ;
  else if( LB == -NOVAL/2 )
    o << ") < " << (UB+1) ;
  else if( UB == NOVAL/2 )
    o << ") > " << (LB-1) ;
  else 
    o << ") in [" << LB << ".." << UB << "]" ;

}

void ConstraintLinearCoef::printFull(std::ostream& o) const
{ 
  int i;
  o<< "(" ;
  for(i=0; i<wpos; ++i) {
    if(i) o << " + ";
    scope[i]->print( o );
  }
  for(i=wpos; i<wneg; ++i) {
    if(i) o << " + " ;
    o << weights[i];
    scope[i]->print( o );
  }
  for(i=wneg; i<arity; ++i) {
    if(i) o << " - " ;
    o << weights[i];
    scope[i]->print( o );
  }
  o << ") " ;
  //o << ") " << = " << weights[arity];
  if( LB == UB )
    o << "= " << UB ;
  else if( LB == -NOVAL/2 )
    o << "< " << (UB+1) ;
  else if( UB == NOVAL/2 )
    o << "> " << (LB-1) ;
  else 
    o << " in [" << LB << ".." << UB << "]" ;
}



// /**********************************************
//  * SumLess Constraint 
//  **********************************************/

// ConstraintSumLess::ConstraintSumLess(Solver *s, VariableInt** v, 
// 				     const int n,
// 				     const int t)
//   : Constraint(s, v, n, RANGETRIGGER)
// {
//   total = t;
// }

// ConstraintSumLess::~ConstraintSumLess()
// {
// }

// int ConstraintSumLess::check( const int* s ) const 
// {
//   int i = arity, tot = 0;
//   while( i-- )
//     tot += s[i];
  
//   return ( tot > total ? tot - total : 0 );
// }

// bool ConstraintSumLess::propagate(const int changedIdx, const int e) 
// {
//   return propagate();
// }

// bool ConstraintSumLess::propagate() 
// {
//   int i;
//   // compute the max and th min
//   int smin = 0;
//   for(i=0; i<arity; ++i) 
//     smin += scope[i]->min();
  
//   if( total < smin ) return false;
//   smin -= total;

//   // prune domains
//   for(i=0; i<arity; ++i) {
//     if( !scope[i]->setMax(scope[i]->min() - smin) ) 
//       return false;
//   }

//   return true;
// }

// void ConstraintSumLess::print(std::ostream& o) const
// { 
//   o<< "(" ;
//   for(int i=0; i<arity-1; ++i) {
//     scope[i]->printshort( o );
//     o << " + ";
//   }
//   scope[arity-1]->printshort( o );
//   o << ") <= " << total;
// }


// // /**********************************************
// //  * Tuple Predicate 
// //  **********************************************/

// // int exponent(int x, int e) {
// //   int res = 1;
// //   while( e-- )
// //     res = (res * x);
// //   return res;
// // }

// // void int2tuple(int x, int* tuple, int n, int* exposant) {
// //   int res = x;
// //   while( n-- ) {
// //     tuple[n] = (res/exposant[n]);
// //     res -= (tuple[n] * exposant[n]);
// //   }
// // }

// // void tuple2int(int& x, int* tuple, int n, int* exposant) {
// //   x = 0;
// //   while( n-- )
// //     x += (tuple[n] * exposant[n]);  
// // }

// // TuplePredicate::TuplePredicate(Solver *s, VariableInt** nvars, 
// // 			       const int n)
// //   : Constraint(s,), BitVar()
// // {
// //   name = "tuple";

// //   int lb=0, ub=0, i, size;
// //   VariableInt *scp[n+1];
// //   vals = new BitSet[n];

// //   exps = new int[n];
// //   for(i=0; i<n; ++i) {
// //     scp[i] = nvars[i];
// //     if(i) {
// //       exps[i] = (scp[i]->max() - scp[i]->min() + 1) * exps[i-1];
// //     } else exps[i] = 1;

// //     lb += (scp[i]->min() * exps[i]);
// //     ub += (scp[i]->max() * exps[i]);

// //     vals[i].init( scp[i]->min(), scp[i]->max(), BitSet::empt );
// //   }

// //   setVariable( lb, ub );
// //   scp[n] = this;
// //   initScope( scp, n+1 );
// // }

// // TuplePredicate::~TuplePredicate()
// // {
// //   delete [] exps;
// //   delete [] vals;
// // }

// // bool TuplePredicate::add( CSP* model, int toplevel )
// // {
// //   VariableInt::add( model, 0 );
// //   Constraint::add( model, 0 );
// //   return true;
// // }

// // int TuplePredicate::check( const int* s ) const 
// // {
// //   int x;
// //   tuple2int( x, s, arity-1, exps ); 
// //   return ( x != s[arity-1] );
// // }

// // bool TuplePredicate::propagate(const int changedIdx, const int e) 
// // {
// //   int vali, vnxt, i, tuple[arity-1];

// //   i = arity-1;
// //   while( i-- ) vals[i].clear();

// //   vali = first();
// //   do {
// //     vnxt = getNext( vali );
// //     int2tuple( vali, tuple, arity-1, exps );
// //     i=arity-2;
// //     while( i >= 0 && scope[i]->contain(tuple[i]) ) --i;
// //     if( i >= 0 && !remove( vali ) ) return false;
// //     vali = vnxt;

// //     if( contain(vali) ) {
// //       i = arity-1;
// //       while( i-- ) vals[i].insert( tuple[i] );
// //     }
// //   }
// //   while( good( vali ) );	

// //   i = arity-1;
// //   while( i-- ) scope[i]->setDomain( vals[i] );
  
// //   return true;
// // }

// // void TuplePredicate::print(std::ostream& o) const
// // { 
// //   BitVar::print( o );
// //   o << " := <";
// //   for(int i=0; i<arity-2; ++i) {
// //     scope[i]->printshort( o );
// //     o << ", ";
// //   }
// //   scope[arity-2]->printshort( o );
// //   o << ">";
// // }


/**********************************************
 * Max Predicate 
 **********************************************/

PredicateMax::PredicateMax(Solver *s, VariableInt** v, const int n)
  : Constraint(s, v, n, RANGETRIGGER)
{
}

PredicateMax::~PredicateMax()
{
}

int PredicateMax::check( const int* s ) const 
{
  int i = arity-2, themax = s[arity-2];
  while( i-- > 0 ) if( s[i] > themax ) themax = s[i];
  return (themax != s[arity-1] );
}

bool PredicateMax::propagate(const int changedIdx, const int e) 
{

  //std::cout << "beg max" << std::endl;

  // NOT GAC !!
  bool consistent = true;
  int i, j=-1, m;
  if( changedIdx == arity-1 ) {
    i = arity-1;
    while( consistent && i-- )
      consistent = scope[i]->setMax( scope[changedIdx]->max() );
  } else {
    consistent = scope[arity-1]->setMin( scope[changedIdx]->min() );
    if( consistent ) {
      int maxmax = scope[arity-2]->max(), maxmin = scope[arity-2]->min();
      i = arity-2;
      while( i-- > 0 ) {
	if( maxmax < scope[i]->max() )
	  maxmax = scope[i]->max();
	if( maxmin < scope[i]->min() )
	  maxmin = scope[i]->min();
      }
      consistent = (scope[arity-1]->setMax( maxmax ) && 
		    scope[arity-1]->setMin( maxmin ));
    }
  }
  if(consistent) {
    i = arity-1, m = scope[i]->min();
    while( i-- ) 
      if( scope[i]->max() >= m )
	{
	  if( j < 0 ) 
	    j = i;
	  else break;
	}
    if( i < 0 )
      consistent = (j >= 0 && scope[j]->setMin( m ));
  }

  //  std::cout << "end max" << std::endl;

  return consistent;
}

void PredicateMax::print(std::ostream& o) const
{ 
  o<< "max(" ;
  for(int i=0; i<arity-2; ++i) {
    scope[i]->printshort( o );
    o << ", ";
  }
  scope[arity-2]->printshort( o );
  o << ") = " ;
  scope[arity-1]->printshort( o );
}


/**********************************************
 * Min Predicate 
 **********************************************/

PredicateMin::PredicateMin(Solver *s, VariableInt** v, const int n)
  : Constraint(s, v, n, RANGETRIGGER)
{
}

PredicateMin::~PredicateMin()
{
}

int PredicateMin::check( const int* s ) const 
{
  int i = arity-2, themin = s[arity-2];
  while( i-- > 0 ) if( s[i] < themin ) themin = s[i];
  return (themin != s[arity-1] );
}

bool PredicateMin::propagate(const int changedIdx, const int e) 
{
  bool consistent = true;

  //std::cout << "beg min" << std::endl;


  // NOT GAC !!
  int i, j=-1, m;

//   if(id == 145) {
//   std::cout << id << " min(";
//   for(i=0; i<arity-1; ++i) {
//     std::cout << " ";
//     scope[i]->print(std::cout);
//   }
//   std::cout << ") = ";
//   scope[i]->print(std::cout);
//   std::cout << " " << (changedIdx) << std::endl;
//   }

  if( changedIdx == arity-1 ) {
    i = arity-1;
    while( consistent && i-- ) {
//       std::cout << "\tset min(";
//       scope[i]->print(std::cout);
//       std::cout << ") to " << scope[changedIdx]->min() << std::endl;
      consistent = scope[i]->setMin( scope[changedIdx]->min() );
      //      std::cout << "\t\t" << consistent << std::endl;
    }
  } else {
    consistent = scope[arity-1]->setMax( scope[changedIdx]->max() );
    if( consistent ) {
      int minmax = scope[arity-2]->max(), minmin = scope[arity-2]->min();
      i = arity-2;
      while( i-- > 0 ) {
	if( minmax > scope[i]->max() )
	  minmax = scope[i]->max();
	if( minmin > scope[i]->min() )
	  minmin = scope[i]->min();
      }
      consistent = (scope[arity-1]->setMax( minmax ) && 
		    scope[arity-1]->setMin( minmin ));
    }
  }
  if(consistent) {
    i = arity-1, m = scope[i]->max();
    while( i-- ) 
      if( scope[i]->min() <= m )
	{
	  if( j < 0 ) 
	    j = i;
	  else break;
	}
    if( i < 0 )
      consistent = (j >= 0 && scope[j]->setMax( m ));
  }

//   if(id == 145) {
//   std::cout << id << " min(";
//   for(i=0; i<arity-1; ++i) {
//     std::cout << " ";
//     scope[i]->print(std::cout);
//   }
//   std::cout << ") = ";
//   scope[i]->print(std::cout);
//   std::cout << " " << (consistent) << std::endl << std::endl;
//   }

//  std::cout << "end min" << std::endl;

  return consistent;
}

void PredicateMin::print(std::ostream& o) const
{ 
  o<< "min(" ;
  for(int i=0; i<arity-2; ++i) {
    scope[i]->printshort( o );
    o << ", ";
  }
  scope[arity-2]->printshort( o );
  o << ") = " ;
  scope[arity-1]->printshort( o );
}


// /**********************************************
//  * Arithmetic and Logic Constraints
//  **********************************************/



// /**********************************************
//  * WeightedSum Predicate 
//  **********************************************/

// PredicateWeightedSum::PredicateWeightedSum( Solver *s, VariableInt** v, 
// 					    const int n, 
// 					    const int* os )
//   : Constraint(s, v, n, RANGETRIGGER)
// {
//   weights = new int[n];
//   int i;
//   for(i=0; i<n; ++i) 
//     weights[i] = os[i];
// }

// PredicateWeightedSum::~PredicateWeightedSum()
// {
//   delete [] weights;
// }

// int PredicateWeightedSum::check( const int* s ) const 
// {
//   int i = arity-1, total = weights[i];
//    while( i-- )
//     total += weights[i] * s[i];

//   return ( total > s[arity-1] ? total - s[arity-1] : s[arity-1] - total );
// }

// bool PredicateWeightedSum::propagate(const int changedIdx, const int e) 
// {
//   return propagate();
// }

// bool PredicateWeightedSum::propagate() 
// {
//   int i, n = arity-1;
//   // compute the max and th min
//   int smin = weights[n], smax=weights[n];
//   for(i=0; i<n; ++i) {
//     if( weights[i] > 0 ) {
//       smax += weights[i] * scope[i]->max();
//       smin += weights[i] * scope[i]->min();
//     } else {
//       smin += weights[i] * scope[i]->max();
//       smax += weights[i] * scope[i]->min();
//     }
//   }

//   if( !scope[n]->setMax( smax ) || !scope[n]->setMin( smin ) ) return false;

//   smax -= scope[n]->min();
//   smin -= scope[n]->max();

//   // prune domains
//   int aux;
//   for(i=0; i<n; ++i) {
//     if( weights[i] > 0 ) {
//       aux = scope[i]->max();
//       if( (!scope[i]->setMax(scope[i]->min() - smin/weights[i])) 
// 	  ||
// 	  (!scope[i]->setMin(aux - smax/weights[i])) ) 
// 	return false;
//     } else {
//       aux = scope[i]->min();
//       if( (!scope[i]->setMin(scope[i]->max() - smin/weights[i]))
// 	  ||
// 	  (!scope[i]->setMax(aux - smax/weights[i])) )
// 	return false;
//     }
//   }

//   return true;
// }

// void PredicateWeightedSum::print(std::ostream& o) const
// { 
//   o<< "(" ;
//   if( weights[0] < 0 ) {
//     o << " - ";
//     o.flush() ; 
//   }
//   if( weights[0] != 1 && weights[0] != -1 ) {
//     o << weights[0];
//     o.flush();
//   }
//   scope[0]->printshort( o );
//   o.flush();
//   for(int i=1; i<arity-1; ++i) {
//     if( weights[i] < 0 ) {
//       o << " - ";
//       o.flush();
//     }
//     else {
//       o << " + ";
//       o.flush();
//     }
//     if( weights[i] != 1 && weights[i] != -1 ) {
//       o << abs(weights[i]);
//       o.flush();
//     }
//     scope[i]->printshort( o );
//     o.flush();
//   }
//   if( weights[arity-1] > 0 ) o << " + " << weights[arity-1];
//   if( weights[arity-1] < 0 ) o << " - " << -weights[arity-1];
//   o << ") = " ;
//   scope[arity-1]->printshort( o );
// }



/**********************************************
 * Nested Offset Predicate
 **********************************************/
/// x0 + k
PredicateOffset::PredicateOffset(Solver *s, VariableInt** v, const int k) 
  : Constraint(s, v, 2, DOMAINTRIGGER, WEIGHTTWO)
{
  offset = k;
  if( offset < 0 )
    xdom.init(scope[0]->min()+offset, scope[0]->max(), BitSet::empt);
  else xdom.init(scope[0]->min(), scope[0]->max()+offset, BitSet::empt);
}

PredicateOffset::~PredicateOffset()
{
}

int PredicateOffset::check( const int* s ) const 
{
  return( s[0]+offset != s[1] );
}

bool PredicateOffset::awakeOnDomain( VariableInt *v )
{
  xdom.clear();
  v->unionTo( xdom );
  
  if( offset < 0 ) {
    if( v == scope[0] )
      xdom.decrement( -offset );
    else xdom.increment( -offset );
  } else {
    if( v == scope[0] )
      xdom.increment( offset );
    else xdom.decrement( offset );
  }    
  return scope[(scope[0] == v)]->setDomain( xdom );     
}

bool PredicateOffset::awakeOnRange( VariableInt *v )
{
  return( v == scope[0] ?
	  (scope[1]->setMin( v->min()+offset ) &&
	   scope[1]->setMax( v->max()+offset )) :
	  (scope[0]->setMin( v->min()-offset ) &&
	   scope[0]->setMax( v->max()-offset )) );
}

bool PredicateOffset::propagate(const int changedIdx, const int e)
{
  return ( (e & Constraint::DOMAINTRIGGER) ?
	   awakeOnDomain( scope[changedIdx] ) :
	   awakeOnRange( scope[changedIdx] ) );
}

void PredicateOffset::print(std::ostream& o) const 
{
  scope[0]->printshort( o );
  if( offset < 0 )
    o << " - " << (-offset) ;
  else
    o << " + " << offset ;
  o << " = ";
  scope[1]->printshort( o );
}


/**********************************************
 * Nested Factor Predicate
 **********************************************/
/// x0 + k
PredicateFactor::PredicateFactor(Solver *s, VariableInt** v, const int k) 
  : Constraint(s, v, 2, DOMAINTRIGGER, WEIGHTTWO)
{
  factor = k;
}

PredicateFactor::~PredicateFactor()
{
}

int PredicateFactor::check( const int* s ) const 
{
  return( s[0] * factor != s[1] );
}

bool PredicateFactor::propagate(const int changedIdx, const int e)
{
  bool consistent = true;
  int lb, ub, rest, i = 1-changedIdx;

  if( i ) // v == scope[0] 
    {
      if( factor > 0 )
	consistent = ( scope[1]->setMin( scope[changedIdx]->min() * factor ) &&
		       scope[1]->setMax( scope[changedIdx]->max() * factor ) );
      else
	consistent = ( scope[1]->setMin( scope[changedIdx]->max() * factor ) &&
		       scope[1]->setMax( scope[changedIdx]->min() * factor ) );
    }
  else // v == scope[1]
    {
      int f = (factor > 0 ? factor : -factor);

      lb = scope[changedIdx]->min();
      rest = (lb % factor);   
   
      if( rest ) {
	if( lb > 0 ) lb += f;
	lb -= rest;

	consistent = scope[changedIdx]->setMin( lb );
	lb = scope[changedIdx]->min();
      }
      
      if( consistent ) {
	ub = scope[changedIdx]->max();
	rest = (ub % factor);

	if( rest ) {
	  if( ub < 0 ) ub -= f;
	  ub -= rest;

	  consistent = scope[changedIdx]->setMax( ub );
	  ub = scope[changedIdx]->max();
	}
	
	if ( factor > 0 ) {
	  consistent &= ( scope[0]->setMin( lb / factor ) &&
			  scope[0]->setMax( ub / factor ) );	
	} else {
	  consistent &= ( scope[0]->setMin( ub / factor ) &&
			  scope[0]->setMax( lb / factor ) );
	}            
      }
    }

  
  return consistent;
}

void PredicateFactor::print(std::ostream& o) const 
{
  o  << factor << "." ;
  scope[0]->printshort( o );
  o << " = ";
  scope[1]->printshort( o );
}


/**********************************************
 * Nested Abs Predicate
 **********************************************/
/// -x0
PredicateAbs::PredicateAbs(Solver *s, VariableInt** v) 
  : Constraint(s, v, 2)
{
}

PredicateAbs::~PredicateAbs()
{
}

int PredicateAbs::check( const int* s ) const 
{
  return( abs(s[0]) != s[1] );
}

bool PredicateAbs::propagate(const int changedIdx, const int e)
{  
  int lb, ub, aux;
  bool consistent = true;

  if( e == VALUETRIGGER ) {
    if( !changedIdx ) {
      consistent = scope[1]->setDomain( abs(scope[changedIdx]->value()) );
    } else {
      ub = scope[changedIdx]->value();
      lb = -ub;
      consistent = (scope[0]->setMin( lb )
		    &&
		    scope[0]->setMax( ub )
		    && 
		    scope[0]->removeRange( lb+1, ub-1 ));
    }
  } else {
    if( scope[changedIdx]->isRange() ) {
      if( !changedIdx ) {
	lb = abs(scope[changedIdx]->min());
	ub = abs(scope[changedIdx]->max());
	if(lb > ub) {
	  aux = ub;
	  ub = lb;
	} else aux = lb;
	consistent = scope[1]->setMax( ub );
      
	if( consistent ) {
	  lb = scope[1]->min();
	  if( !(scope[changedIdx]->contain(lb) || scope[changedIdx]->contain(-1*lb)) ) {
	    --lb;
	    ub = aux;
	    aux = ((ub + lb) >> 1);
	    while( ub > lb+1 ) {
	      if( scope[changedIdx]->intersect(-1*aux, aux) ) {
		ub = aux;
		aux = ((ub + lb) >> 1);
	      } else {
		lb = aux;
		aux = ((ub + lb) >> 1);
	      }
	    }  
	    consistent = scope[1]->setMin( aux+1 );
	  }
	}
      } else {
	consistent = ((scope[0]->setMax(scope[1]->max())) 
		      &&
		      (scope[0]->setMin(-1*scope[1]->max()))
		      &&
		      (!scope[1]->min() || 
		       scope[0]->removeRange( -1*scope[1]->min()+1, 
					      scope[1]->min()-1)));
      }
    } else {
      if( !changedIdx ) {
// 	lb = scope[1]->first();
// 	do {
// 	  aux = scope[1]->getNext( lb );
// 	  if( !(scope[0]->contain(lb) || scope[0]->contain(-lb)) )
// 	    consistent = scope[1]->remove( lb );
// 	  lb = aux;
// 	} while( consistent && scope[1]->good( lb ) );

	// std::cout << "NOT CHECKED 1" << std::endl;

	DomainIterator *valit = scope[1]->begin();
	do {
	  lb=*valit;
	  if( !(scope[0]->contain(lb) || scope[0]->contain(-lb)) )
	    consistent = scope[1]->remove( lb );
	} while( consistent && valit->next() );

      } else {

	// std::cout << "NOT CHECKED 2" << std::endl;

	DomainIterator *valit = scope[0]->begin();
	do {
	  lb=*valit;
	  if( !(scope[1]->contain(lb) || scope[1]->contain(-lb)) )
	    consistent = scope[0]->remove( lb );
	} while( consistent && valit->next() );

// 	lb = scope[0]->first();
// 	do {
// 	  aux = scope[0]->getNext( lb );
// 	  if( !scope[1]->contain(abs(lb)) )
// 	    consistent = scope[0]->remove( lb );
// 	  lb = aux;
// 	} while( consistent && scope[0]->good( lb ) );
      }
    }
  }

  return consistent;
}

void PredicateAbs::print(std::ostream& o) const 
{
  o << "abs(" ;
  scope[0]->printshort( o );
  o << ") = ";
  scope[1]->printshort( o ); 
}


/**********************************************
 * Nested Negation Predicate
 **********************************************/
/// -x0
PredicateNegation::PredicateNegation(Solver *s, VariableInt** v) 
  : Constraint(s, v, 2, DOMAINTRIGGER, WEIGHTTWO)
{
  dom[0].init(v[0]->min(), v[0]->max(), BitSet::empt);
  dom[1].init(v[1]->min(), v[1]->max(), BitSet::empt);
}

PredicateNegation::~PredicateNegation()
{
}

int PredicateNegation::check( const int* s ) const 
{
  return( s[0] != -s[1] );
}

bool PredicateNegation::propagate(const int changedIdx, const int e)
{
  dom[0].clear();
  dom[1].clear();

  scope[changedIdx]->unionTo( dom[changedIdx] );
  dom[changedIdx].negate( dom[1-changedIdx] );
  return scope[1-changedIdx]->setDomain( dom[1-changedIdx] );
}

void PredicateNegation::print(std::ostream& o) const 
{
  o << "-" ;
  scope[0]->printshort( o );
  o << " = ";
  scope[1]->printshort( o );
}


/**********************************************
 * Nested Not Predicate
 **********************************************/
/// not(x0)
PredicateNot::PredicateNot(Solver *s, VariableInt** v) 
  : Constraint(s, v, 2, DOMAINTRIGGER, WEIGHTTWO)
{
}

PredicateNot::~PredicateNot()
{
}

int PredicateNot::check( const int* s ) const 
{
  return( !(s[0]) != s[1] );
}

bool PredicateNot::propagate()
{
  bool consistent = true;
  for( int i=0; consistent && i<2; ++ i )
    consistent = ( (!scope[i]->min() || scope[1-i]->setDomain( 0 )) &&
		   (scope[i]->max() || scope[1-i]->remove( 0 ) ) );
  return consistent;
}

bool PredicateNot::propagate(const int changedIdx, const int e)
{
  return ( (!scope[changedIdx]->min() || scope[1-changedIdx]->setDomain( 0 )) &&
	   (scope[changedIdx]->max() || scope[1-changedIdx]->remove( 0 ) ) );
}

void PredicateNot::print(std::ostream& o) const 
{
  o << "!" ;
  scope[0]->printshort( o );
  o << " = ";
  scope[1]->printshort( o );
}


/**********************************************
 * Nested Equal Predicate
 **********************************************/
/// x0 =?= x1

PredicateEqual::PredicateEqual(Solver *s, VariableInt** v, int sp) 
  : Constraint(s, v, 3, DOMAINTRIGGER, WEIGHTTWO)
{
  spin = sp;
}

PredicateEqual::~PredicateEqual()
{
}

int PredicateEqual::check( const int* s ) const 
{
  return( (s[0] == s[1]) == (s[2] ^ spin) );
}

bool PredicateEqual::propagate(const int changedIdx, const int e)
{

  bool consistent = true;

  if( scope[2]->isGround() ) {
    if( (spin + scope[2]->min()) != 1 ) {
      consistent = ( scope[0]->setDomain(scope[1]) &&
		     scope[1]->setDomain(scope[0]) );
    } else {
      consistent = ( (!scope[0]->isGround() || scope[1]->remove(scope[0]->value()))
		     &&
		     (!scope[1]->isGround() || scope[0]->remove(scope[1]->value())) );
    }
  } else {
    if( !(scope[0]->intersect(scope[1])) ) {
      consistent = ( scope[2]->remove( spin ) );
    } else { 
      if( scope[0]->isGround() && scope[1]->isGround() ) {
	consistent = ( scope[2]->setDomain( spin ) );
      }
    }
  }

  return consistent;
}

void PredicateEqual::print(std::ostream& o) const 
{
  o << "(";
  scope[0]->printshort( o );
  if( spin )
    o << " == ";
  else
    o << " =/= ";
  scope[1]->printshort( o );
  o << ") = ";
  scope[2]->printshort( o );
}


/**********************************************
 * Nested EqualConstant Predicate
 **********************************************/
/// x0 =?= k

PredicateEqualConstant::PredicateEqualConstant(Solver *s, 
					       VariableInt** v, 
					       const int k,
					       const int sp) 
  : Constraint(s, v, 2, DOMAINTRIGGER, WEIGHTTHREE)
{
  K = k;
  spin = sp;
}

PredicateEqualConstant::~PredicateEqualConstant()
{
}

int PredicateEqualConstant::check( const int* s ) const 
{
  return( (s[0] == K) == (s[1] ^ spin) );
}

bool PredicateEqualConstant::propagate(const int changedIdx, const int e)
{
  bool consistent = true;
  if( scope[1]->isGround() ) {
    if( (spin + scope[1]->min()) != 1 ) {
      consistent = scope[0]->setDomain(K);
    } else {      
      consistent = scope[0]->remove(K);
    }
  } else {
    if( !(scope[0]->contain( K )) ) {
      consistent = scope[1]->remove( spin );
    } else { 
      if( scope[0]->isGround() ) {
	consistent = scope[1]->setDomain( spin );
      }
    }
  }
  return consistent;
}

void PredicateEqualConstant::print(std::ostream& o) const 
{
  o << "equalconstant  (";
  scope[0]->print( o );
  if( spin )
    o << " == ";
  else
    o << " =/= ";
  o << K << ") = ";
  scope[1]->print( o );
}



/**********************************************
 * Nested Member Predicate
 **********************************************/
/// x0 =?= k

PredicateMember::PredicateMember(Solver *s, 
				 VariableInt** v, 
				 const int *k,
				 const int l,
				 const int sp) 
  : Constraint(s, v, 2, DOMAINTRIGGER, WEIGHTTHREE)
{
  K.init( scope[0]->min(), scope[0]->max(), BitSet::empt );
  for(int i=0; i<l; ++i)
    K.insert(k[i]);
  spin = sp;
}

PredicateMember::~PredicateMember()
{
}

int PredicateMember::check( const int* s ) const 
{
  return( (K.fastMember(s[0])) == (s[1] ^ spin) );
}

bool PredicateMember::propagate(const int changedIdx, const int e)
{
  bool consistent = true;
  if( scope[1]->isGround() ) {
    if( (spin + scope[1]->min()) != 1 ) {
      consistent = scope[0]->setDomain(K);
    } else {      
      consistent = scope[0]->removeSet(K);
    }
  } else {
    if( !(scope[0]->intersect( K )) ) {
      consistent = scope[1]->remove( spin );
    } else { 
      if( scope[0]->isGround() ) {
	consistent = scope[1]->setDomain( spin );
      }
    }
  }
  return consistent;
}

void PredicateMember::print(std::ostream& o) const 
{
  o << "equalconstant  (";
  scope[0]->print( o );
  if( spin )
    o << " in ";
  else
    o << " not in ";
  K.print(o);
  o << ") = ";
  scope[1]->print( o );
}


/**********************************************
 * Nested And Predicate
 **********************************************/
/// x0 &?& x1
PredicateAnd::PredicateAnd(Solver *s, VariableInt** v) 
  : Constraint(s, v, 3, DOMAINTRIGGER, WEIGHTTHREE), 
    domain_x(v[0]->getIntDomain()),
    domain_y(v[1]->getIntDomain()),
    domain_z(v[2]->getIntDomain())

    /*
    domain_x(((VariableBool*)(v[0]))->domain),
    domain_y(((VariableBool*)(v[1]))->domain),
    domain_z(((VariableBool*)(v[2]))->domain)
    domain_z((v[2]->getType() == VariableInt::CONST ? (((Constant*)(v[2]))->val ? is_true : is_false) : ((VariableBool*)(v[2]))->domain ))
    */
{

  //std::cout << (domain_x) << " " << (domain_y) << " " << (domain_z) << " " << std::endl;

  //std::cout << sizeof(Constant) << " " << std::endl;
  //std::cout << sizeof(VariableBool) << " " << std::endl;

  //std::cout << (&domain_z) << " " << std::endl;
  //std::cout << (domain_z.state) << " " << std::endl;
  //if(v[2]->getType() == VariableInt::CONST) domain_z.state = (((Constant*)(v[2]))->val + 1);
  //std::cout << domain_z.state << std::endl;
  //exit(1);
}

PredicateAnd::~PredicateAnd()
{
}

int PredicateAnd::check( const int* s ) const 
{
  return( (s[0] && s[1]) != s[2] );
}

// bool PredicateAnd::propagate(const int changedIdx, const int e)
// {
// }


bool PredicateAnd::propagate(const int changedIdx, const int e)
{
//   if( scope[2]->min() ) {
//     if( !scope[0]->setDomain( 1 ) ||
// 	!scope[1]->setDomain( 1 ) )
//       return false;
//   } else if( !(scope[2]->max()) ) {
//     if( (scope[0]->min() && !scope[1]->setDomain( 0 ))
// 	||
// 	(scope[1]->min() && !scope[0]->setDomain( 0 )) )
//       return false;
//   }
//   if( (!scope[0]->max() || !scope[1]->max()) &&
//       !scope[2]->setDomain( 0 ) ) return false;
//   else if( scope[0]->min() && scope[1]->min() &&
// 	   !scope[2]->setDomain( 1 ) ) return false;
//   return true;


//   if(scope[0]->id == 32 && scope[1]->id == 33) {
//     scope[0]->print(std::cout);
//     std::cout << " ";
//     scope[1]->print(std::cout);
      
//   }


//   bool consistent = true;
//   if( domain_z.state == 2 ) {

// //     if( !(domain_x.state & 2) || 
// // 	!(domain_y.state & 2) ) consistent = false;
// //     else {
// //       scope[0]->setDomain(1);
// //       scope[1]->setDomain(1);
// //     }

//      if( !scope[0]->setDomain( 1 ) ||
// 	 !scope[1]->setDomain( 1 ) )
//        consistent = false;

//   } else if( domain_z.state == 1 ) {
//     if( (domain_x.state == 2 && !scope[1]->setDomain( 0 ))
// 	||
// 	(domain_y.state == 2 && !scope[0]->setDomain( 0 )) )
//       consistent = false;
//   }
//   if( consistent ) {
//     if( (domain_x.state == 1 || domain_y.state == 1) &&
// 	!scope[2]->setDomain( 0 ) ) consistent = false;
//     else if( domain_x.state == 2 && domain_y.state == 2 &&
// 	     !scope[2]->setDomain( 1 ) ) consistent = false;
//   } 
//   return consistent;


  bool consistent = true;
  if( domain_z == 2 ) {

//     if( !(domain_x & 2) || 
// 	!(domain_y & 2) ) consistent = false;
//     else {
//       scope[0]->setDomain(1);
//       scope[1]->setDomain(1);
//     }

     if( !scope[0]->setDomain( 1 ) ||
	 !scope[1]->setDomain( 1 ) )
       consistent = false;

  } else if( domain_z == 1 ) {
    if( (domain_x == 2 && !scope[1]->setDomain( 0 ))
	||
	(domain_y == 2 && !scope[0]->setDomain( 0 )) )
      consistent = false;
  }
  if( consistent ) {
    if( (domain_x == 1 || domain_y == 1) &&
	!scope[2]->setDomain( 0 ) ) consistent = false;
    else if( domain_x == 2 && domain_y == 2 &&
	     !scope[2]->setDomain( 1 ) ) consistent = false;
  } 
  return consistent;
}


void PredicateAnd::print(std::ostream& o) const 
{
  o << "(";
  scope[0]->printshort( o );
  o << " && " ;
  scope[1]->printshort( o );
  o << ") = ";
  scope[2]->printshort( o );
}


/**********************************************
 * Nested Or Predicate
 **********************************************/
/// x0 &?& x1
PredicateOr::PredicateOr(Solver *s, VariableInt** v) 
  : Constraint(s, v, 3)
{
}

PredicateOr::~PredicateOr()
{
}

int PredicateOr::check( const int* s ) const 
{
  return( (s[0] || s[1]) != s[2] );
}

bool PredicateOr::propagate(const int changedIdx, const int e)
{
  bool consistent = true;
  if( !scope[2]->max() ) {
    consistent = ( scope[0]->setDomain( 0 ) &&
		   scope[1]->setDomain( 0 ) );
  } else if( scope[2]->min() ) {
    consistent = ( ( scope[0]->max() || scope[1]->setDomain( 1 ) )
		   &&
		   ( scope[1]->max() || scope[0]->setDomain( 1 ) ) );
  }
  if( consistent ) {
    if( scope[0]->min() || scope[1]->min() )
      consistent = scope[2]->setDomain( 1 );
    else if( !scope[0]->max() && !scope[1]->max() )
      consistent = scope[2]->setDomain( 0 );
  }
  return consistent;
}

void PredicateOr::print(std::ostream& o) const 
{
  o << "(";
  scope[0]->printshort( o );
  o << " || " ;
  scope[1]->printshort( o );
  o << ") = ";
  scope[2]->printshort( o );
}


/**********************************************
 * IfThenElse Predicate
 **********************************************/

PredicateIfThenElse::PredicateIfThenElse(Solver *s, VariableInt** v)
  : Constraint(s, v, 4, DOMAINTRIGGER)
{
}

PredicateIfThenElse::~PredicateIfThenElse() {}

int PredicateIfThenElse::check( const int* s ) const 
{
  return ( (s[0] ? s[1] : s[2]) != s[3] );
}

bool PredicateIfThenElse::propagate(const int changedIdx, const int e)
{
  bool consistent = true;
  if( !scope[3]->intersect( scope[1] ) )
    consistent = scope[0]->remove(1);
  else if( !scope[3]->intersect( scope[2] ) )
    consistent = scope[0]->remove(0);
      
  if( scope[0]->min() )
    consistent = (scope[3]->setDomain( scope[1] ) && scope[1]->setDomain( scope[3] ));
  else if( !scope[0]->max() )
    consistent = (scope[3]->setDomain( scope[2] ) && scope[2]->setDomain( scope[3] ));
  return consistent;
}

void PredicateIfThenElse::print(std::ostream& o) const
{
  scope[3]->print( o ) ;
  o << " =";
  o << "(" ;
  scope[0]->print( o ) ;
  o << " ? " ;
  scope[1]->print( o ) ;
  o << " : " ;
  scope[2]->print( o ) ;
  o << ")";
}

/**********************************************
 * <= Predicate 
 **********************************************/
/// (x0 + k < x1)
PredicateLess::PredicateLess(Solver *s, VariableInt** v, const int offs) 
  : Constraint(s, v, 3, RANGETRIGGER, WEIGHTTWO)
{
  offset = offs;
}

PredicateLess::~PredicateLess()
{
}

int PredicateLess::check( const int* s ) const 
{
  return ( (s[0] + offset <= s[1]) != s[2] );
}

bool PredicateLess::propagate(const int changedIdx, const int e) 
{
  bool consistent = true;

  if( scope[2]->domsize() < 2 ) {
    if( scope[2]->min() ) {
      consistent = (scope[1]->setMin(scope[0]->min() + offset) &&
		    scope[0]->setMax(scope[1]->max() - offset));
    } else {
      consistent = (scope[0]->setMin(scope[1]->min() - offset + 1) &&
		    scope[1]->setMax(scope[0]->max() + offset - 1));
    }
  } else {
    if( scope[0]->max() + offset <= scope[1]->min() ) 
      consistent = ( scope[2]->remove( 0 ) );
    else if( scope[0]->min() + offset > scope[1]->max() )
      consistent = ( scope[2]->remove( 1 ) );
  }

  return consistent;
}

void PredicateLess::print(std::ostream& o) const
{  
  o << "(";
  scope[0]->printshort( o );
  if( !offset )
    o << " <= ";
  else if( offset == 1 )
    o << " < ";
  else 
    o << " + " << offset-1 << " < " ;
  scope[1]->printshort( o );
  o << ") = ";
  scope[2]->printshort( o );
}


/**********************************************
 * >= k Predicate 
 **********************************************/
/// (x0 + k < x1)
PredicateLowerBound::PredicateLowerBound(Solver *s, VariableInt** v, const int b) 
  : Constraint(s, v, 2, RANGETRIGGER, WEIGHTTHREE)
{
  bound = b;
}

PredicateLowerBound::~PredicateLowerBound()
{
}

int PredicateLowerBound::check( const int* s ) const 
{
  return ( (s[0] >= bound) != s[1] );
}

bool PredicateLowerBound::propagate(const int changedIdx, const int e) 
{
  bool consistent = true;

  if( scope[1]->domsize() < 2 ) {
    if( scope[1]->min() ) 
      consistent = scope[0]->setMin(bound);
    else 
      consistent = scope[0]->setMax(bound-1);    
  } else {
    if( scope[0]->max() < bound  ) 
      consistent =  scope[1]->remove( 1 );
    else if( scope[0]->min() >= bound )
      consistent = scope[1]->remove( 0 );
  }

  return consistent;
}

void PredicateLowerBound::print(std::ostream& o) const
{  
  o << "(";
  scope[0]->printshort( o );
  o << " > " << (bound-1) << ") = ";
  scope[1]->printshort( o );
}


/**********************************************
 * >= k Predicate 
 **********************************************/
/// (x0 + k < x1)
PredicateUpperBound::PredicateUpperBound(Solver *s, VariableInt** v, const int b) 
  : Constraint(s, v, 2, RANGETRIGGER, WEIGHTTHREE)
{
  bound = b;
}

PredicateUpperBound::~PredicateUpperBound()
{
}

int PredicateUpperBound::check( const int* s ) const 
{
  return ( (s[0] <= bound) != s[1] );
}

bool PredicateUpperBound::propagate(const int changedIdx, const int e) 
{
  bool consistent = true;
  if( scope[1]->domsize() < 2 ) {
    if( scope[1]->min() ) 
      consistent = scope[0]->setMax(bound);
    else 
      consistent = scope[0]->setMin(bound+1);    
  } else {
    if( scope[0]->max() <= bound  ) 
      consistent =  scope[1]->remove( 0 );
    else if( scope[0]->min() > bound )
      consistent = scope[1]->remove( 1 );
  }
  return consistent;
}

void PredicateUpperBound::print(std::ostream& o) const
{  
  o << "(";
  scope[0]->printshort( o );
  o << " < " << (bound+1) << ") = ";
  scope[1]->printshort( o );
}


/**********************************************
 * Add Predicate 
 **********************************************/
/// (x0 + x1)
PredicateAdd::PredicateAdd(Solver *s, VariableInt** v) 
  : Constraint(s, v, 3, DOMAINTRIGGER, WEIGHTTHREE)
{
  // int lb = std::min(v[2]->min(), std::min(v[0]->min(), v[1]->min()));
  // int ub = std::max(v[2]->max(), std::max(v[0]->max(), v[1]->max()));

//   if( ub - lb < 1024 && 
//       scope[0]->getType() != VariableInt::VIRTUAL &&
//       scope[1]->getType() != VariableInt::VIRTUAL &&
//       scope[2]->getType() != VariableInt::VIRTUAL 
//       ) {
//     xdom.init(lb, ub, BitSet::empt);
//     tdom.init(lb, ub, BitSet::empt);
//   }
}

PredicateAdd::~PredicateAdd()
{
}

int PredicateAdd::check( const int* s ) const 
{
  return ( (s[0] + s[1]) != s[2] );
}

bool PredicateAdd::propagate(const int changedIdx, const int e) 
{
  int i;
  bool consistent = true;//, hole =false;
  
//   if(xdom.table)
//     for(i=0; i<3 && !hole; ++i) 
//       if( !scope[i]->isRange() )
// 	hole = true;

//   if( hole ) {
//     //int vali, vnxt;
    
//     /// revise d(scope[2])
//     if( changedIdx != 2 && consistent ) {
//       tdom.clear();
//       i = ( scope[0]->domsize() > scope[1]->domsize() );

//       // std::cout << "NOT CHECKED 3" << std::endl;

//       DomainIterator *valit = scope[i]->begin();
//       do {
// 	xdom.clear();
// 	scope[1-i]->unionTo( xdom );
// 	if( *valit < 0 )
// 	  xdom.decrement( -(*valit) );
// 	else
// 	  xdom.increment( *valit );
// 	tdom.unionWith( xdom );
//       } while( valit->next() );

//  //      vali = scope[i]->first();
// //       do {
// // 	vnxt = scope[i]->getNext( vali );
// // 	xdom.clear();
// // 	scope[1-i]->unionTo( xdom );
// // 	if( vali < 0 )
// // 	  xdom.decrement( -vali );
// // 	else
// // 	  xdom.increment( vali );
// // 	tdom.unionWith( xdom );
// // 	vali = vnxt;
// //       } while( scope[i]->good( vali ) );
//       consistent = scope[2]->setDomain( tdom );
//     }
    
//     /// revise d(scope[0])
//     if( changedIdx && consistent ) {  
//       tdom.clear();
// //       vali = scope[1]->first();
// //       do {
// // 	vnxt = scope[1]->getNext( vali );
// // 	xdom.clear();
// // 	scope[2]->unionTo( xdom );
// // 	if( vali < 0 )
// // 	  xdom.increment( -vali );
// // 	else
// // 	  xdom.decrement( vali );
// // 	tdom.unionWith( xdom );
// // 	vali = vnxt;
// //       } while( scope[1]->good( vali ) );

//       // std::cout << "NOT CHECKED 4" << std::endl;

//       DomainIterator *valit = scope[1]->begin();
//       do {
// 	xdom.clear();
// 	scope[2]->unionTo( xdom );
// 	if( *valit < 0 )
// 	  xdom.increment( -(*valit) );
// 	else
// 	  xdom.decrement( *valit );
// 	tdom.unionWith( xdom );
//       } while( valit->next() );

//       consistent = scope[0]->setDomain( tdom );
//     }

//     /// revise d(scope[1])
//     if( changedIdx != 1 && consistent ) {  
//       tdom.clear();

//       // std::cout << "NOT CHECKED 5" << std::endl;

//  //      vali = scope[0]->first();
// //       do {
// // 	vnxt = scope[0]->getNext( vali );
// // 	xdom.clear();
// // 	scope[2]->unionTo( xdom );
// // 	if( vali < 0 )
// // 	  xdom.increment( -vali );
// // 	else
// // 	  xdom.decrement( vali );
// // 	tdom.unionWith( xdom );
// // 	vali = vnxt;
// //       } while( scope[0]->good( vali ) );

//       DomainIterator *valit = scope[0]->begin();
//       do {
// 	xdom.clear();
// 	scope[2]->unionTo( xdom );
// 	if( *valit < 0 )
// 	  xdom.increment( -(*valit) );
// 	else
// 	  xdom.decrement( *valit );
// 	tdom.unionWith( xdom );
//       } while( valit->next() );

//       consistent = scope[1]->setDomain( tdom );
//     }

//   } else {

    int lb[3], ub[3];
    for(i=0; i<2; ++i)
      {
	lb[i] = scope[i]->min();
	ub[i] = scope[i]->max();
      }
    

// //     std::cout << std::endl << lb[0] << " + " << lb[1] << " <= " << scope[2]->min() << std::endl;
// //     std::cout << ub[0] << " + " << ub[1] << " >= " << scope[2]->max() << std::endl;

//      scope[2]->print( std::cout ) ;
//      std::cout << " >= " << (lb[0]+lb[1]) << std::endl;
//      scope[2]->print( std::cout ) ;
//      std::cout << " <= " << (ub[0]+ub[1]) << std::endl << std::endl;


    consistent = (scope[2]->setMin( lb[0]+lb[1] ) &&
		  scope[2]->setMax( ub[0]+ub[1] ));


    if( consistent ) {
      lb[2] = scope[2]->min();
      ub[2] = scope[2]->max();

      
//        scope[0]->print( std::cout ) ;
//        std::cout << " >= " << (lb[2]-ub[1]) << std::endl;
//        scope[0]->print( std::cout ) ;
//        std::cout << " <= " << (ub[2]-lb[1]) << std::endl;

//        scope[1]->print( std::cout ) ;
//        std::cout << " >= " << (lb[2]-ub[0]) << std::endl;
//        scope[1]->print( std::cout ) ;
//        std::cout << " <= " << (ub[2]-lb[0]) << std::endl;


      consistent = (scope[0]->setMin( lb[2] - ub[1] ) &&
		    scope[0]->setMax( ub[2] - lb[1] ) &&
		 
		    scope[1]->setMin( lb[2] - ub[0] ) &&
		    scope[1]->setMax( ub[2] - lb[0] ));

      //std::cout << consistent << std::endl;


//       std::cout << lb[2] << " - " << ub[1] << " <= " << lb[0] << std::endl;
//       std::cout << ub[2] << " - " << lb[1] << " >= " << ub[0] << std::endl;
//       std::cout << lb[2] << " - " << ub[0] << " <= " << lb[1] << std::endl;
//       std::cout << ub[2] << " - " << lb[0] << " >= " << ub[1] << std::endl << std::endl;
    

      //     int lb = (scope[0]->min() + scope[1]->min());
      //     int ub = (scope[0]->max() + scope[1]->max());
      
      //     consistent = (scope[2]->setMin( lb ) &&
      // 		  scope[2]->setMax( ub ) &&
      
      // 		  scope[0]->setMin( scope[2]->min() - scope[1]->max() ) &&
      // 		  scope[0]->setMax( scope[2]->max() - scope[1]->min() ) &&
      
      // 		  scope[1]->setMin( scope[2]->min() - scope[0]->max() ) &&
      // 		  scope[1]->setMax( scope[2]->max() - scope[0]->min() ) );
      
    }
    //  }

  return consistent;
}

void PredicateAdd::print(std::ostream& o) const
{  
  scope[0]->printshort( o );
  o << " + ";
  scope[1]->printshort( o );
  o << " = ";
  scope[2]->printshort( o );
}


/**********************************************
 * Sub Predicate 
 **********************************************/
/// (x0 - x1)
PredicateSub::PredicateSub(Solver *s, VariableInt** v) 
  : Constraint(s, v, 3, DOMAINTRIGGER, WEIGHTTHREE)
{
  //int lb = std::min(scope[2]->min(), std::min(v[0]->min(), v[1]->min()));
  //int ub = std::max(scope[2]->max(), std::max(v[0]->max(), v[1]->max()));

//   if( !s->boundOptimised && ub - lb < 1024 ) {
//     xdom.init(lb, ub, BitSet::empt);
//     tdom.init(lb, ub, BitSet::empt);
//   }
}

PredicateSub::~PredicateSub()
{
}

int PredicateSub::check( const int* s ) const 
{
  return ( (s[0] - s[1]) != s[2] );
}

bool PredicateSub::propagate(const int changedIdx, const int e) 
{
  //int i;
  bool consistent = true;//, hole =false;

//   for(i=0; i<3 && !hole; ++i)
//     hole = !(scope[i]->isRange()); 
  
//   if( hole && xdom.table ) {
    
//     /// revise d(scope[2])
//     if( changedIdx != 2 && consistent ) {
//       tdom.clear();
//       DomainIterator *valit = scope[1]->begin();
//       tdom.clear();
//       do {
// 	xdom.clear();
// 	scope[0]->unionTo( xdom );
// 	if( *valit < 0 ) 
// 	  xdom.increment( -(*valit) );
// 	else 
// 	  xdom.decrement( *valit );
// 	tdom.unionWith( xdom );
//       } while( valit->next() );
//       consistent = scope[2]->setDomain( tdom );
//     }
    
//     /// revise d(scope[0])
//     if( changedIdx && consistent ) {  
//       tdom.clear();
//       DomainIterator *valit = scope[1]->begin();
//       do {
// 	xdom.clear();
// 	scope[2]->unionTo( xdom );
// 	if( *valit < 0 )
// 	  xdom.decrement( -(*valit) );
// 	else
// 	  xdom.increment( *valit );
// 	tdom.unionWith( xdom );
//       } while( valit->next() );
//       consistent = scope[0]->setDomain( tdom );
//     }

//     /// revise d(scope[1])
//     if( changedIdx != 1 && consistent ) {  
//       tdom.clear();
//       DomainIterator *valit = scope[2]->begin();
//       do {
// 	xdom.clear();
// 	scope[0]->unionTo( xdom );
// 	if( *valit < 0 )
// 	  xdom.increment( -(*valit) );
// 	else
// 	  xdom.decrement( *valit );
// 	tdom.unionWith( xdom );
//       } while( valit->next() );
//       consistent = scope[1]->setDomain( tdom );
//     }
//   } else {

  int lb = (scope[0]->min() - scope[1]->max());
  int ub = (scope[0]->max() - scope[1]->min());
  
  consistent = (scope[2]->setMin( lb ) &&
		scope[2]->setMax( ub ) &&
		
		scope[0]->setMin( scope[2]->min() + scope[1]->min() ) &&
		scope[0]->setMax( scope[2]->max() + scope[1]->max() ) &&
		
		scope[1]->setMin( scope[0]->min() - scope[2]->max() ) &&
		scope[1]->setMax( scope[0]->max() - scope[2]->min() ) );

  //}
  
  return consistent;
}

void PredicateSub::print(std::ostream& o) const
{  
  scope[0]->printshort( o );
  o << " - ";
  scope[1]->printshort( o );
  o << " = ";
  scope[2]->printshort( o );
}


/**********************************************
 * Mul Predicate 
 **********************************************/
/// (x0 * x1)
PredicateMul::PredicateMul(Solver *s, VariableInt** v) 
  : Constraint( s, v, 3 )
{
  int lb = std::min(v[2]->min(), std::min(v[0]->min(), v[1]->min()));
  int ub = std::max(v[2]->max(), std::max(v[0]->max(), v[1]->max()));

  if( ub - lb < 1024 ) {
    xdom.init(lb, ub, BitSet::empt);
  }
  //tdom.init(lb, ub);
}

PredicateMul::~PredicateMul()
{
}

int PredicateMul::check( const int* s ) const 
{
  return ( (s[0] * s[1]) != s[2] );
}


bool PredicateMul::propagate(const int changedIdx, const int e) 
{

//   if(scope[0]->id == 1079) {
//     std::cout << "propagate " ;
//     print(std::cout); 
//     std::cout << " because of ";
//     scope[changedIdx]->printshort(std::cout);
//     std::cout << std::endl;
//     print( std::cout );
//     std::cout << std::endl;
//   }

  int v;
  int i = (changedIdx+1)%3;
  int j = (changedIdx+2)%3;
  //int consistent = pruneZeros(changedIdx);
  int consistent, k=3;
  do consistent = pruneZeros(--k);
  while( consistent == 1 && k );


  if( consistent > 0 ) {
    /// first, we check the particular cases:
    /*************************************************
     * Case 1: two variables are ground 
     * Case 2: the calling variable is ground
     * Case 3: one variable is ground
     *************************************************/
    
    if( scope[changedIdx]->isGround() )
      {
	if( scope[i]->isGround() ) {
	  //std::cout << "FC" << std::endl;
	  consistent = pruneUnary(j);
	} else if( scope[j]->isGround() ) {
	  //std::cout << "FC" << std::endl;
	  consistent = pruneUnary(i);
	} else {
	  v = scope[changedIdx]->first();
	  //std::cout << "AC-1" << std::endl;
	  consistent = (pruneBinary(i, j, v) &&
			pruneBinary(j, i, v));
	}
      }
    else {
      switch( scope[i]->isGround() + 2*scope[j]->isGround() ) {
      case 0: {
// 	if(scope[0]->id == 1079) {
// 	  std::cout << "GAC" << std::endl;
// 	}
	consistent = (pruneTernary(i) && pruneTernary(j));
      } break;
      case 1: {
	v = scope[i]->first();
// 	if(scope[0]->id == 1079) {
// 	  std::cout << "AC-2" << std::endl;
// 	}
	consistent = pruneBinary(changedIdx, j, v);
      } break;
      case 2: {
	v = scope[j]->first();
// 	if(scope[0]->id == 1079) {
// 	  std::cout << "AC-2" << std::endl;
// 	}
	consistent = pruneBinary(changedIdx, i, v);
      } break;
      default: {
	//consistent = pruneUnary(changedIdx);
	//std::cout << "nothing" << std::endl;
      }
      }
    }
  }
//  else {
//     if(scope[0]->id == 1079) {
//       std::cout << "no need to check further" << std::endl;
//     }
//   }
  
//   if(scope[0]->id == 1079) {
//        if( consistent ) {
//     print( std::cout );
//      std::cout << " ok" << std::endl;
//    } else {
//      std::cout << "INCONSISTENT" << std::endl;
//    }
//   std::cout << std::endl;
//   }

  return consistent;
}


inline int xtimey( const int x, const int y, int& r )
{
  //std::cout << x << " * " << y << std::endl;
  r = 0 ;
  return (x*y);
}

inline int xovery( const int x, const int y, int& r )
{
  //  std::cout << x << " x/y " << y << " = " << (x/y) 
  //   << " r " << (x%y != 0) << std::endl;
  r = (x%y != 0);
  return (x/y);
}

inline int yoverx( const int x, const int y, int& r )
{
  //std::cout << y << " y/x " << x << " = "  << (y/x) 
  //   << " r "<< (x%y != 0) << std::endl;
  r = -(y%x != 0);
  return (y/x);
}


int PredicateMul::pruneZeros(const int changedIdx)
{
  ///std::cout << "check the zeros" << std::endl;

//   bool sprint = ( scope[0]->id == 5740 && scope[1]->id == 5741 );
    
//   if( sprint )
//     {
//       std::cout << "check the zeros" << std::endl;
//       print( std::cout );
//       std::cout << " " << changedIdx << std::endl;
//     }

  int consistent = 1;
  if( !scope[changedIdx]->contain(0) ) { // scope[changedIdx] cannot be 0
    if( changedIdx == 2 ) { // the product is not 0

 //      if( sprint )
// 	std::cout << "scope[0]->remove(0) && scope[1]->remove(0)" << std::endl;

      consistent = (scope[0]->remove(0) && scope[1]->remove(0));
    } else { // one factor is not zero     
      consistent = ( (scope[1-changedIdx]->contain(0) || scope[2]->remove(0)) &&
		     (!scope[2]->isGround() || scope[2]->first() || scope[1-changedIdx]->setDomain(0)) ); 
    }
  } else if( scope[changedIdx]->isGround() ) { // scope[changedIdx] must be 0
    if( changedIdx == 2 ) { // the product is 0
      consistent = -( (scope[0]->contain(0) || scope[1]->setDomain(0)) &&
		      (scope[1]->contain(0) || scope[0]->setDomain(0)) );
    } else { // one factor is 0
      consistent = -(scope[2]->setDomain(0));
    }
  }

  return consistent;
}

bool PredicateMul::pruneUnary(const int last) 
{
  if( last < 2 ) {
    int v = scope[1-last]->first();
    int w = scope[2]->first();
    if(v) return (!(w%v) && scope[last]->setDomain( w/v ));      
    else return !w;
  } else
    return scope[last]->setDomain( scope[0]->first() * scope[1]->first() );
}

bool PredicateMul::pruneBinary(const int otherIdx, const int reviseIdx, const int v)
{

//   std::cout << "revise " ;
//   scope[reviseIdx]->print(std::cout);
//   std::cout << " with respect to ";
//   scope[otherIdx]->print(std::cout);
//   std::cout << std::endl;

  // bool yox = false, xoy = false;

  int (*oper)(const int, const int, int&);
  if( otherIdx+reviseIdx != 1 ) {
    if( reviseIdx == 2 ) oper = xtimey;
    else {
      //yox = true;
      oper = yoverx;
    }
  } else {
    //xoy = true;
    oper = xovery;
  }


  int w1 = scope[otherIdx]->min(), w2 = scope[otherIdx]->max(), lb, ub, r;
  if(v<0) {
    lb = w1;
    w1 = w2;
    w2 =lb;
  }


 //   if( (xoy && !w1) || (yox && !v) ) {

//     print( std::cout );
//     std::cout << std::endl;

//     std::cout << "revise " ;
//     scope[reviseIdx]->print(std::cout);
//     std::cout << " with respect to ";
//     scope[otherIdx]->print(std::cout);
//     std::cout << " (" << v << ")" << std::endl;    
//   }


  lb = oper(v, w1, r);
//   if(r) {
//     std::cout << "**" << std::endl;
//     ++lb;
//   }
  switch( r ) {
  case -1: {
    //std::cout << "lb is not exact" << std::endl;
    //std::cout << w1 << "/" << v << " =/= " << lb << std::endl;    
    if(v>0) {
      //std::cout << v << "*" << lb << " < " << w1 << std::endl;
      if(v*lb < w1) ++lb;
    } else {
      //std::cout << v << "*" << lb << " > " << w1 << std::endl;
      if(v*lb > w1) ++lb;
    }
    //std::cout << " => " << lb << std::endl;
  } break; 
  case 1: {
    //std::cout << "lb is not exact" << std::endl;
    //std::cout << v << "/" << w1 << " =/= " << lb << std::endl;    
    if(v>0) {
      //std::cout << v << "/" << lb << " < " << w1 << "?" << std::endl;
      if(v/lb < w1) ++lb;
    } else {
      //std::cout << v << "/" << lb << " > " << w1 << "?" << std::endl;
      if(v/lb > w1) ++lb;
    }
    //std::cout << " => " << lb << std::endl;
  } break;
//   default : {
//     std::cout << "lb is exact" << std::endl;
//   }
  }

  ub = oper(v, w2, r);
//   if(r) {
//     std::cout << "gg" << std::endl;
//     --ub;
//   }
  switch( r ) {
  case -1: {
    //std::cout << "ub is not exact" << std::endl;
    //std::cout << w2 << "/" << v << " =/= " << ub << std::endl;
    if(v>0) {
      //std::cout << v << "*" << ub << " > " << w2 << "?" << std::endl;
      if(v*ub > w2) --ub;
    } else {
      //std::cout << v << "*" << ub << " < " << w2 << "?" << std::endl;
      if(v*ub < w2) --ub;
    }
    //std::cout << " => " << ub << std::endl;
  } break; 
  case 1: {
    //std::cout << "ub is not exact" << std::endl;
    //std::cout << v << "/" << w2 << " =/= " << ub << std::endl;    
    if(v>0) {
      //std::cout << v << "/" << ub << " > " << w2 << "?" << std::endl;
      if(v/ub > w2) --ub;
    } else {
      //std::cout << v << "/" << ub << " < " << w2 << "?" << std::endl;
      if(v/ub < w2) --ub;
    }
    //std::cout << " => " << ub << std::endl;
  } break;
//   default : {
//     std::cout << "ub is exact" << std::endl;
//   }
  }
  
  if( !(scope[reviseIdx]->setMin(lb) && scope[reviseIdx]->setMax(ub)) ) return false;
  

  if( !scope[otherIdx]->isRange() && !scope[reviseIdx]->isRange() && xdom.table )
    {
      //std::cout << "there are holes, we achieve AC" << std::endl; 
      xdom.clear();
      DomainIterator *valit = scope[otherIdx]->begin();
      do xdom.insert( oper(v, *valit, r) );
      while( valit->next() ); 
      if( !scope[reviseIdx]->setDomain( xdom ) ) return false;
    }
  return true;
}

bool PredicateMul::pruneTernary(const int reviseIdx)
{

//   if(scope[0]->id == 1079) {
//     std::cout << "prune ternary " ;
//     scope[reviseIdx]->print( std::cout );
//     std::cout << ": ";
//     print(std::cout);
//     std::cout << std::endl;
//   }


  if( reviseIdx == 2 || !scope[2]->contain(0) || !scope[1-reviseIdx]->contain(0) )
    {

      int (*oper)(const int, const int, int&);
      int x, y, bound[6], zero[3], lb, ub, v[4], i, r;
      
      for(i=0; i<3; ++i)
	{      
	  bound[2*i] = scope[i]->min(); 
	  bound[2*i+1] = scope[i]->max(); 

	  zero[i] = 0;
	  if(!bound[2*i]) { ++bound[2*i]; ++zero[i]; }
	  if(!bound[2*i+1]) { --bound[2*i+1]; zero[i]+=2; }

// 	  if(scope[0]->id == 1079) {
// 	    std::cout << "{" << bound[2*i] << "," << bound[2*i+1] << "} " << zero[i] << std::endl;
// 	  }
	}

      if(reviseIdx == 2) {
	oper = xtimey;
	x = 0;
	y = 1;
      } else {
	oper = xovery;
	x = 2;
	y = 1-reviseIdx;
      }

      for(i=0; i<4; ++i) {
	v[i] = oper(bound[2*x+i/2], bound[2*y+i%2], r);

// 	  if(scope[0]->id == 1079) {
// 	std::cout << bound[2*x+i/2] << " * " << bound[2*y+i%2] << " == " << v[i] << std::endl;
// 	  }
      }

      ub = 0;
      lb = 0;

      for(i=1; i<4; ++i) 
	if( v[i] > v[ub] ) ub = i;
	else if( v[i] < v[lb] ) lb = i;

      lb = v[lb];
      ub = v[ub];
  
// 	  if(scope[0]->id == 1079) {
//       std::cout << "[" << lb << "," << ub << "]" << std::endl;
// 	  }

	  if(zero[reviseIdx]) {
	    if(lb > 0) lb = 0;
	    if(ub < 0) ub = 0;
	  }
	  // 	    // WARNING CHANGE, MAY BE BUGGY!!
//       if( (zero[reviseIdx] & 1) && lb > 0 ) lb = 0;
//       if( (zero[reviseIdx] & 2) && ub < 0 ) ub = 0;

 // 	  if(scope[0]->id == 1079) {
//       std::cout << "[" << lb << "," << ub << "]" << std::endl;
// 	  }

      return ( (scope[reviseIdx]->setMin( lb )) &&
	       (scope[reviseIdx]->setMax( ub )) );
    }
  return true;
}


void PredicateMul::print(std::ostream& o) const
{  
  o << "(";
  scope[0]->print( o );
  o << " * ";
  scope[1]->print( o );
  o << ") = ";
  scope[2]->print( o );
}


/**********************************************
 * Div Predicate 
 **********************************************/
/// (x0 / x1)
PredicateDiv::PredicateDiv(Solver *s, VariableInt** v) 
  : Constraint(s, v, 3)
{
}

PredicateDiv::~PredicateDiv()
{
}

int PredicateDiv::check( const int* s ) const 
{
  return ( (s[0] / s[1]) != s[2] );
}

bool PredicateDiv::propagate(const int changedIdx, const int e) 
{
  // NOT GAC !!
  if( scope[0]->isGround() && 
      !(scope[0]->max()) ) {
    if( !(scope[2]->setDomain( 0 )) ) return false;
    else return true;
  }

  if( !(scope[1]->remove( 0 ) ) ) return false;

  if( scope[0]->isGround() && scope[1]->isGround() )
    {      
      if( !(scope[2]->setDomain( (scope[0]->min() / scope[1]->min()) )) )
	return false;
      else return true;
    }

  if( scope[1]->isGround() && scope[2]->isGround() )
    {
      if( !scope[2]->max() ) {
	if( scope[1]->min() < 0 ) {
	  if( !(scope[0]->setMin( scope[1]->min()+1 )) ||
	      !(scope[0]->setMax( -1*scope[1]->min()-1 )) )
	    return false;
	} else if( !(scope[0]->setMax( scope[1]->min()-1 )) ||
		   !(scope[0]->setMin( -1*scope[1]->min()+1 )) )
	  return false;
      } else {
	int aux, ub, lb = (scope[1]->min() * scope[2]->min());
	
	aux = ( scope[1]->min() > 0 ? scope[1]->min() : -1*scope[1]->min() );
	--aux;
	if( lb < 0 ) {
	  ub = lb;
	  lb = lb - aux;
	} else ub = lb + aux;
	if( !(scope[0]->setMin( lb )) ||
	    !(scope[0]->setMax( ub )) )
	  return false;
	else return true;
      }
    }

  if( scope[0]->isGround() && scope[2]->isGround() )
    {
      if( !scope[2]->max() ) {
	if( scope[0]->min() < 0 ) {
	  if( !(scope[1]->removeRange( scope[0]->min(),
				       -1*scope[0]->min() )) )
	    return false;
	} else if( !(scope[1]->removeRange( -1*scope[0]->min(),
					    scope[0]->min() )) )
	  return false;
      } else {
	int aux, lb = 0, ub = (scope[0]->min() / scope[2]->min());
	if( scope[2]->min() > 0 )
	  lb = (scope[0]->min() / (scope[2]->min()+1));
	else 
	  lb = (scope[0]->min() / (scope[2]->min()-1));
	if( ub > 0 ) ++lb;
	else --lb;
	
	if( lb > ub ) {
	  aux = lb;
	  lb = ub;
	  ub = aux;
	}
	
	if( !(scope[1]->setMin( lb )) ||
	    !(scope[1]->setMax( ub )) )
	  return false;
	else return true;
      }
    }
  return true;
}

void PredicateDiv::print(std::ostream& o) const
{  
  o << "(";
  scope[0]->printshort( o );
  o << " / ";
  scope[1]->printshort( o );
  o << ") = ";
  scope[2]->printshort( o );
}


/**********************************************
 * Mod Predicate 
 **********************************************/
/// (x0 % x1)
PredicateMod::PredicateMod(Solver *s, VariableInt** v) 
  : Constraint(s, v, 3)
{
}

PredicateMod::~PredicateMod()
{
}

int PredicateMod::check( const int* s ) const 
{
  return ( (s[0] % s[1]) != s[2] );
}

bool PredicateMod::propagate(const int changedIdx, const int e) 
{
  bool consistent = true;
  int ub, lb, modulo, target, k;  
  if( changedIdx ) { // prune x0
    if( scope[1]->isGround() ) { // modulo is known
      modulo = scope[1]->min();
      if( modulo == 1 ) // special case
	consistent = scope[0]->setDomain(0);
      else if( scope[2]->isGround() ) { // target is known
	target = scope[2]->min();
	consistent = (target < modulo);
      
	// positive/negative target
	if( target < 0 )
	  consistent = scope[0]->setMax(0);
	else if( target > 0 )
	  consistent = scope[0]->setMin(0);
	else 
	  consistent = scope[0]->setDomain(0);

	if( consistent ) {
	  // remove intervals [target+1+k*modulo..target+k*(modulo+1)-1]
	  k = (scope[0]->max()-target-1)/modulo;
	  while( consistent ) {
	    lb = (target+1+k*modulo);
	    ub = (target+(k+1)*modulo-1);
	    if( ub < scope[0]->min() ) break;
	    consistent = scope[0]->removeRange( lb, ub );
	    --k;
	  }
	}      
      } else {
	// prune x0 with respect to x2
	
	DomainIterator *valit = scope[0]->begin();
	do {
	k = (*valit % modulo);
	consistent = ( scope[2]->contain( k ) || scope[0]->remove( *valit ) );
	} while( consistent && valit->next() );
	
      }
    } else {
      // modulo is not known we want to prune x0

    }
  } else if( changedIdx != 2 ) {
    // prune x2

    if( scope[1]->isGround() ) { // modulo is known

      if( scope[0]->isGround() )
	
	consistent = scope[2]->setDomain( (scope[0]->min() % scope[1]->min()) );

      else {
      
	modulo = scope[1]->min();
	ub = scope[0]->max();
	if( ub > 0 && modulo <= ub )
	  ub = modulo-1;
	else ub = 0;
	lb = scope[0]->min();
	if( lb < 0 && 1-modulo > lb )
	  lb = 1-modulo;
	else lb = 0;
	
	consistent = scope[2]->setMax( ub ) && scope[2]->setMin( lb );
	
	DomainIterator *valit = scope[2]->begin();
	do {
	  k = *valit;
	  if( k > 0 ) {
	    lb = (scope[0]->min()/modulo)*modulo;
	    k = std::max( k, lb+k );
	    while( !scope[0]->contain(k) && k <= scope[0]->max() )  {
	      k+=modulo;
	    }
	    if(k && k > scope[0]->max()) {
	      consistent = scope[2]->remove( *valit );
	    } 
	  }
	  else {
	    ub = (scope[0]->max()/modulo)*modulo;
	    if(k) k = std::min( k, ub+k );
	    else k = ub;
	    while( !scope[0]->contain(k) && k >= scope[0]->min() ) {
	      k-=modulo;
	    }
	    if(k < scope[0]->min()) {
	      consistent = scope[2]->remove( *valit );
	    } 
	  }
	} while( consistent && valit->next() );
      } 
    } else {
      // modulo is not known we want to prune x2
    }
  } 

  return consistent;
}

void PredicateMod::print(std::ostream& o) const
{  
  o << "(";
  scope[0]->printshort( o );
  o << " % ";
  scope[1]->printshort( o );
  o << ") = ";
  scope[2]->printshort( o );
}


/**********************************************
 * ModConstant Predicate 
 **********************************************/
/// (x0 % x1)
PredicateModConstant::PredicateModConstant(Solver *s, VariableInt** v, const int mod) 
  : Constraint(s, v, 2)
{
  modulo = mod;
}

PredicateModConstant::~PredicateModConstant()
{
}

int PredicateModConstant::check( const int* s ) const 
{
  return ( (s[0] % modulo) != s[1] );
}

bool PredicateModConstant::propagate(const int changedIdx, const int e) 
{
  bool consistent = true;
  int lb, ub, target, k;  
  if( changedIdx ) { // prune x0

    if( modulo == 1 ) // special case
      consistent = scope[0]->setDomain(0);

    else if( scope[1]->isGround() ) { // target is known
      target = scope[1]->min();
      consistent = (target < modulo);
      
      // positive/negative target
      if( target < 0 )
	consistent = scope[0]->setMax(0);
      else if( target > 0 )
	consistent = scope[0]->setMin(0);
      else 
	consistent = scope[0]->setDomain(0);

      if( consistent ) {
	// remove intervals [target+1+k*modulo..target+k*(modulo+1)-1]
	k = (scope[0]->max()-target-1)/modulo;
	while( consistent ) {
	  lb = (target+1+k*modulo);
	  ub = (target+(k+1)*modulo-1);
	  if( ub < scope[0]->min() ) break;
	  consistent = scope[0]->removeRange( lb, ub );
	  --k;
	}
      }
    } else {
      // prune x0 with respect to x1

      DomainIterator *valit = scope[0]->begin();
      do {
	k = (*valit % modulo);
	consistent = ( scope[1]->contain( k ) || scope[0]->remove( *valit ) );
      } while( consistent && valit->next() );
      
    }
  } else {

    if( scope[0]->isGround() )
	consistent = scope[2]->setDomain( (scope[0]->min() % scope[1]->min()) );
    else {
      ub = scope[0]->max();
      if( ub > 0 && modulo <= ub )
	ub = modulo-1;
      else ub = 0;
      
      lb = scope[0]->min();
      if( lb < 0 && 1-modulo > lb )
	lb = 1-modulo;
      else lb = 0;
      
      consistent = scope[1]->setMax( ub ) && scope[1]->setMin( lb );
      
      
      DomainIterator *valit = scope[1]->begin();
      do {
	k = *valit;

	if( k > 0 ) {
	  
	  lb = (scope[0]->min()/modulo)*modulo;
	  k = std::max( k, lb+k );
	  
	  while( !scope[0]->contain(k) && k <= scope[0]->max() )  {
	    k+=modulo;
	  }
	  if(k && k > scope[0]->max()) {
	    consistent = scope[1]->remove( *valit );
	  } 
	}
	
	else {
	  ub = (scope[0]->max()/modulo)*modulo;
	  if(k) k = std::min( k, ub+k );
	  else k = ub;

	  while( !scope[0]->contain(k) && k >= scope[0]->min() ) {
	    k-=modulo;
	  }
	  if(k < scope[0]->min()) {
	    consistent = scope[1]->remove( *valit );
	  } 
	}
	
      } while( consistent && valit->next() );
    }

  }

  return consistent;
}

void PredicateModConstant::print(std::ostream& o) const
{  
  o << "(";
  scope[0]->printshort( o );
  o << " % " << modulo << ") = ";
  scope[1]->printshort( o );
}


/**********************************************
 * Pow Predicate 
 **********************************************/
/// (x0 ^ x1)
PredicatePow::PredicatePow(Solver *s, VariableInt** v) 
  : Constraint(s, v, 3)
{
  std::cerr << "Not implemented" << std::endl;
  exit( 0 );
}

PredicatePow::~PredicatePow()
{
}

int PredicatePow::check( const int* s ) const 
{
  return ( pow((double)(s[0]), (double)(s[1])) != (double)(s[2]) );
}

bool PredicatePow::propagate(const int changedIdx, const int e) 
{
  std::cerr << "Not implemented" << std::endl;
  exit( 0 );
  return true;
}

void PredicatePow::print(std::ostream& o) const
{  
  o << "(";
  scope[0]->printshort( o );
  o << " ^ ";
  scope[1]->printshort( o );
  o << ") = ";
  scope[2]->printshort( o );
}


/**********************************************
 * Disjunctive Predicate 
 **********************************************/
/// 0 :-> (x0 + d0 < x1), 1 :-> (x1 + d1 < x0) 
PredicateDisjunctive::PredicateDisjunctive(Solver *s, VariableInt** v, const int* dur) 
  : Constraint(s, v, 3, RANGETRIGGER), 
    state(((VariableBool*)(v[2]))->domain.state),
    min_0(((VariableRange*)(v[0]))->vmin.value),
    min_1(((VariableRange*)(v[1]))->vmin.value),
    max_0(((VariableRange*)(v[0]))->vmax.value),
    max_1(((VariableRange*)(v[1]))->vmax.value)
{
  duration[0] = dur[0];
  duration[1] = dur[1];
}

PredicateDisjunctive::~PredicateDisjunctive()
{
}

int PredicateDisjunctive::check( const int* s ) const 
{
//   return ( (s[2] && (s[0] + duration[0] <= s[1])) 
// 	   ||
// 	   (!s[2] && (s[1] + duration[1] <= s[0])) );

//   print(std::cout);
//   std::cout << std::endl;

//   std::cout << "s[0] = " << s[0] << " s[1] = " 
//  	   << s[1] << " s[2] = " << s[2] << std::endl;




  return ( (s[2] && (s[1] + duration[1] > s[0])) 
	   ||
	   (!s[2] && (s[0] + duration[0] > s[1])) );
}

bool PredicateDisjunctive::propagate(const int changedIdx, const int e) 
{

  //assert( state)

//   assert( min_0 == scope[0]->min() );
//   assert( min_1 == scope[1]->min() );
//   assert( max_0 == scope[0]->max() );
//   assert( max_1 == scope[1]->max() );


//   if( changedIdx != 2 )
//     for(int i=0; i<2; ++i)
//       if( scope[i]->min() + duration[i] > scope[1-i]->max() 
// 	  && !scope[2]->setDomain( 1-i ) ) return false;

//   if( !scope[2]->max() ) {
//     // x + a <= y
//     return ( scope[1]->setMin(scope[0]->min() + duration[0]) &&
// 	     scope[0]->setMax(scope[1]->max() - duration[0]) );

//   } else if( scope[2]->min() ) {
//     // y + b <= x
//     return ( scope[0]->setMin(scope[1]->min() + duration[1]) &&
// 	     scope[1]->setMax(scope[0]->max() - duration[1]) );    
//   } else {
//     int MaxIncX = scope[1]->min() + duration[1] -1;
//     int MinIncX = scope[1]->max() - duration[0] +1;

//     if( MaxIncX >= MinIncX && 
// 	!scope[0]->removeRange(MinIncX, MaxIncX) )
//       return false;
    
//     MaxIncX = scope[0]->min() + duration[0] -1;
//     MinIncX = scope[0]->max() - duration[1] +1;

//     if( MaxIncX >= MinIncX && 
// 	!scope[1]->removeRange(MinIncX, MaxIncX) )
//       return false;
//   }


  if( changedIdx != 2 ) {
    if( min_0 + duration[0] > max_1 ) {
      if( !scope[2]->setDomain(1) ) return false;
    } else if( min_1 + duration[1] > max_0 ) {
      if( !scope[2]->setDomain(0) ) return false;
    }
  }


  int bound;
  if( state == 1 ) {
    // x + a <= y
    bound = min_0 + duration[0];
    if( bound > min_1 ) {
      if( !scope[1]->setMin( bound ) ) return false;
    }
    bound = max_1 - duration[0];
    if( bound < max_0 ) {
      if( !scope[0]->setMax( bound ) ) return false;
    }

  } else if( state == 2 ) {
    // y + b <= x
    bound = min_1 + duration[1];
    if( bound > min_0 ) {
      if( !scope[0]->setMin( bound ) ) return false;
    }
    bound = max_0 - duration[1];
    if( bound < max_1 ) {
      if( !scope[1]->setMax( bound ) ) return false;
    }

  } 

//   if( state == 1 ) {
//     // x + a <= y
//     return ( scope[1]->setMin(min_0 + duration[0]) &&
// 	     scope[0]->setMax(max_1 - duration[0]) );

//   } else if( state == 2 ) {
//     // y + b <= x
//     return ( scope[0]->setMin(min_1 + duration[1]) &&
// 	     scope[1]->setMax(max_0 - duration[1]) );    
//   } 
  
  return true;
}

void PredicateDisjunctive::print(std::ostream& o) const
{  
  if( scope[2]->contain(0) ) {
    o << "(";
    scope[0]->printshort( o );
    if( !duration[0] )
      o << " <= ";
    else if( duration[0] == 1 )
      o << " < ";
    else 
      o << " + " << duration[0]-1 << " < " ;
    scope[1]->printshort( o );
    o << ")";
  }
  if( scope[2]->domsize() == 2 )
    o << " || ";
  if( scope[2]->contain(1) ) {
    o << "(";
    scope[1]->printshort( o );
    if( !duration[1] )
      o << " <= ";
    else if( duration[1] == 1 )
      o << " < ";
    else 
      o << " + " << duration[1]-1 << " < " ;
    scope[0]->printshort( o );
    o << ") = ";
  }
  //  o << " [" << vmin << ".." << vmax << "]";
  scope[2]->printshort( o );
}



/**********************************************
 * Overlap Predicate 
 **********************************************/
/// 0 :-> (x0 + d0 < x1)
/// 1 :-> (x1 + d1 < x0) 
/// 2 :-> (x0 + d0 > x1) & (x1 + d1 > x0)
PredicateOverlap::PredicateOverlap(Solver *s, VariableInt** v, const int* dur) 
  : Constraint(s, v, 3, RANGETRIGGER), 
    state(((VariableBit*)(v[2]))->values.table[0]),
    min_0(((VariableRange*)(v[0]))->vmin.value),
    min_1(((VariableRange*)(v[1]))->vmin.value),
    max_0(((VariableRange*)(v[0]))->vmax.value),
    max_1(((VariableRange*)(v[1]))->vmax.value)
{
  duration[0] = dur[0];
  duration[1] = dur[1];
  triggers[2] = _scope[2]->triggerOnDomain();
}

PredicateOverlap::~PredicateOverlap()
{
}

int PredicateOverlap::check( const int* s ) const 
{
  return (((s[2] != 2) || ((s[0] + duration[0] <= s[1]) || (s[1] + duration[1] <= s[0]))) &&
	  ((s[2] != 0) || (s[0] + duration[0] > s[1])) &&
	  ((s[2] != 1) || (s[1] + duration[1] > s[0])));
}

bool PredicateOverlap::propagate(const int changedIdx, const int e) 
{

  bool consistent = true;
  int bound;

  if(scope[0]->id == 386 && scope[1]->id == 392) {
    print( std::cout );
    std::cout << std::endl;
    for(int a=0; a<arity; ++a) {
      scope[a]->print( std::cout );
      std::cout << " ";
    }
    std::cout << std::endl;
  }

  if(changedIdx != 2) {
    if(min_0 + duration[0] > max_1) {
      if(!scope[2]->remove(0)) consistent = false;
    } else if(min_1 + duration[1] > max_0 ) {
      if(!scope[2]->remove(1)) consistent = false;
    } else if((max_0 + duration[0] <= min_1) || (max_1 + duration[1] <= min_0)) {
      if(!scope[2]->remove(2)) consistent = false;
    }
  }


  if(consistent) {

    switch(state) {
    case 1 : { // {0} -> x0 << x1
      // x0 + d0 <= x1
      bound = min_0 + duration[0];
      if( bound > min_1 ) consistent &= scope[1]->setMin( bound );      
      bound = max_1 - duration[0];
      if( bound < max_0 ) consistent &= scope[0]->setMax( bound );
    } break;
    case 2 : { // {1} -> x1 << x0
      // x1 + d1 <= x0
      bound = min_1 + duration[1];
      if( bound > min_0 ) consistent &= scope[0]->setMin( bound );      
      bound = max_0 - duration[1];
      if( bound < max_1 ) consistent &= scope[1]->setMax( bound );
    } break;
    case 3 : { // {0,1} -> no overlap
      if( min_0 + duration[0] > max_1 ) scope[2]->setDomain(1);
      else if( min_1 + duration[1] > max_0 ) !scope[2]->setDomain(0);
    } break;
    case 4 : { // {2} -> overlap
      // x1 + d1 - 1 >= x0
      // x0 - d1 + 1 <= x1
      bound = min_0 - duration[1] + 1;
      if( bound > min_1 ) consistent &= scope[1]->setMin( bound );      
      bound = max_1 + duration[1] - 1;
      if( bound < max_0 ) consistent &= scope[0]->setMax( bound );
      // x0 + d0 - 1 >= x1
      // x1 - d0 + 1 <= x0
      bound = min_1 - duration[0] + 1;
      if( bound > min_0 ) consistent &= scope[0]->setMin( bound );      
      bound = max_0 + duration[0] - 1;
      if( bound < max_1 ) consistent &= scope[1]->setMax( bound );
    } break;
    case 5 : { // {0,2} -> not( x1 << x0 )
      // x1 + d1 - 1 >= x0
      // x0 - d1 + 1 <= x1
      bound = min_0 - duration[1] + 1;
      if( bound > min_1 ) consistent &= scope[1]->setMin( bound );      
      bound = max_1 + duration[1] - 1;
      if( bound < max_0 ) consistent &= scope[0]->setMax( bound );
    } break;
    case 6 : { // {1,2} -> not( x0 << x1 )
      // x0 + d0 - 1 >= x1
      // x1 - d0 + 1 <= x0
      bound = min_1 - duration[0] + 1;
      if( bound > min_0 ) consistent &= scope[0]->setMin( bound );      
      bound = max_0 + duration[0] - 1;
      if( bound < max_1 ) consistent &= scope[1]->setMax( bound );
    } break;
    }
  }


  if(scope[0]->id == 386 && scope[1]->id == 392) {
    for(int a=0; a<arity; ++a) {
      scope[a]->print( std::cout );
      std::cout << " ";
    }
    std::cout << std::endl;
    print( std::cout );
    std::cout << (consistent ? " OK" : " FAIL") << std::endl << std::endl;
  }
  return consistent;
}

void PredicateOverlap::print(std::ostream& o) const
{  
  scope[2]->print(o);
  o << " <=> {(";
  scope[0]->printshort( o );
  o << " + " << duration[0]-1 << " < " ;
  scope[1]->printshort( o );
  o << ") | (";
  scope[1]->printshort( o );
  o << " + " << duration[1]-1 << " < " ;
  scope[0]->printshort( o );
  o << ") | (";
  scope[0]->printshort( o );
  o << " overlap with " ;
  scope[1]->printshort( o );
  o << ")}";
}


/**********************************************
 * SafeDisjunct Predicate 
 **********************************************/
/// 0 :-> (x0 + d0 < x1), 1 :-> (x1 + d1 < x0) 
PredicateSafeDisjunct::PredicateSafeDisjunct(Solver *s, VariableInt** v, const int* dur) 
  : Constraint(s, v, 3, RANGETRIGGER)
{
  duration[0] = dur[0];
  duration[1] = dur[1];
}

PredicateSafeDisjunct::~PredicateSafeDisjunct()
{
}

int PredicateSafeDisjunct::check( const int* s ) const 
{

  //std::cout << "s[0] = " << s[0] << " s[1] = " 
  //<< s[1] << " s[2] = " << s[2] << std::endl;



  return ( (s[2] && (s[1] + duration[1] > s[0])) 
	   ||
	   (!s[2] && (s[0] + duration[0] > s[1])) );


//   return ( (s[2] && (s[0] + duration[0] <= s[1])) 
// 	   ||
// 	   (!s[2] && (s[1] + duration[1] <= s[0])) );
}

bool PredicateSafeDisjunct::propagate(const int changedIdx, const int e) 
{

//   std::cout << "PROPAGATE " << scope[2]->id << " " << scope[0]->id << " " << scope[1]->id << std::endl;
//   if( scope[2]->id == 26 && scope[0]->id == 43 ) {
//     scope[changedIdx]->print(std::cout);
//     std::cout << std::endl;
//     print(std::cout);
//     std::cout << std::endl;
//   }

  if( changedIdx != 2 )
    for(int i=0; i<2; ++i)
      if( scope[i]->min() + duration[i] > scope[1-i]->max() 
	  && !scope[2]->setDomain( 1-i ) ) return false;

  if( !scope[2]->max() ) {
    // x + a <= y
    return ( scope[1]->setMin(scope[0]->min() + duration[0]) &&
	     scope[0]->setMax(scope[1]->max() - duration[0]) );

  } else if( scope[2]->min() ) {
    // y + b <= x
    return ( scope[0]->setMin(scope[1]->min() + duration[1]) &&
	     scope[1]->setMax(scope[0]->max() - duration[1]) );    
  } 

//   if( scope[2]->id == 26 && scope[0]->id == 43 ) {
//     print(std::cout);
//     std::cout << std::endl << std::endl;
//   }

  return true;
}

void PredicateSafeDisjunct::print(std::ostream& o) const
{
  o << "(safe) ";
  if( scope[2]->contain(0) ) {
    o << "(";
    scope[0]->print( o );
    if( !duration[0] )
      o << " <= ";
    else if( duration[0] == 1 )
      o << " < ";
    else 
      o << " + " << duration[0]-1 << " < " ;
    scope[1]->print( o );
    o << ")";
  }
  if( scope[2]->domsize() == 2 )
    o << " || ";
  if( scope[2]->contain(1) ) {
    o << "(";
    scope[1]->print( o );
    if( !duration[1] )
      o << " <= ";
    else if( duration[1] == 1 )
      o << " < ";
    else 
      o << " + " << duration[1]-1 << " < " ;
    scope[0]->print( o );
    o << ") = ";
  }
  //  o << " [" << vmin << ".." << vmax << "]";
  scope[2]->print( o );
}


/**********************************************
 * Hole Predicate 
 **********************************************/
/// 0 :-> (x0 + d0 < x1), 1 :-> (x1 + d1 < x0) 
PredicateHole::PredicateHole(Solver *s, VariableInt** v, const int* dur) 
  : Constraint(s, v, 2, RANGETRIGGER), 
    state(((VariableBool*)(v[1]))->domain.state),
    min_x(((VariableRange*)(v[0]))->vmin.value),
    max_x(((VariableRange*)(v[0]))->vmax.value)
{
  bounds[0] = dur[0];
  bounds[1] = dur[1];
}

PredicateHole::~PredicateHole()
{
}

int PredicateHole::check( const int* s ) const 
{
  return ( (s[1] && (s[0] > bounds[0]))
	   ||
	   (!s[1] && (s[0] < bounds[1])) );
}

bool PredicateHole::propagate(const int changedIdx, const int e) 
{

  if( changedIdx != 1 ) {
    if( min_x > bounds[0] ) {
      if( !scope[1]->setDomain(1) ) return false;
    } else if( bounds[1] > max_x ) {
      if( !scope[1]->setDomain(0) ) return false;
    }
  }

  if( state == 1 ) {
    if( !scope[0]->setMax( bounds[0] ) ) return false;
  } else if( state == 2 ) {
    if( !scope[0]->setMin( bounds[1] ) ) return false;
  } 
  
  return true;
}

void PredicateHole::print(std::ostream& o) const
{  
  if( scope[2]->contain(0) ) {
    o << "(";
    scope[0]->printshort( o );
    o << " <= " << bounds[0] << " || ";
    scope[0]->printshort( o );
    o << " => " << bounds[1] << ") ";
  }
  scope[2]->printshort( o );
}


/**********************************************
 * Element Predicate
 **********************************************/ 
/// x[N]
PredicateElement::PredicateElement(Solver *s, VariableInt** v, const int n, const int o)
  : Constraint(s, v, n)
{
//   aux_dom.init( std::min( 0, scope[n-1]->min() ), 
// 		std::max( n-1, scope[n-1]->max() ), 
// 		BitSet::empt );
  aux_dom.init( std::min( 0, scope[n-1]->min() ), 
		std::max( n-1, scope[n-1]->max() ), 
		BitSet::empt );
  offset = o;

}

PredicateElement::~PredicateElement()
{
}

int PredicateElement::check( const int* s ) const 
{ 
  //return ( s[s[arity-2]] != s[arity-1] );
  return ( s[s[arity-2]-offset] != s[arity-1] );
}

bool PredicateElement::propagate(const int changedIdx, const int e)
{
  int i = arity-2;
  VariableInt *N = scope[i], *V = scope[i+1];
  aux_dom.clear();
  bool consistent = true;

  /* propagation rule for N */ 
  if( changedIdx != i ) {
    while( i-- ) 
      if( V->intersect(scope[i]) )
	//aux_dom.insert(i);
	aux_dom.insert(i+offset);
    consistent = ( N->setDomain(aux_dom) );
  }

  /* propagation rule for V */ 
  if( changedIdx != arity-1 && consistent ) {
    aux_dom.clear();

//     i = N->first();
//     do scope[i-1]->unionTo(aux_dom); 
//     while( N->setNext(i) );

    DomainIterator *valit = N->begin();
    do scope[*valit-offset]->unionTo(aux_dom); 
    //do scope[*valit]->unionTo(aux_dom); 
    while( valit->next() );

    consistent = ( V->setDomain(aux_dom) );
  }

  /* propagation rule for x[N] */ 
  if( N->isGround() ) 
    //consistent = ( scope[N->min()]->setDomain(V) &&
    //V->setDomain(scope[N->min()]) );
    consistent = ( scope[N->min()-offset]->setDomain(V) &&
		   V->setDomain(scope[N->min()-offset]) );

  return consistent;
}

void PredicateElement::print(std::ostream& o) const
{
  o << "x"
    // << scope[0]->id << ".." << scope[arity-2]->id
    << "[" ;
  scope[arity-2]->printshort( o );
  o << "] = " ;
  scope[arity-1]->printshort( o );
}

ConstraintLDS::ConstraintLDS(Solver* s, VariableInt** v, const int n)
  : Constraint( s, v, n )
{
}

ConstraintHamming::ConstraintHamming(Solver* s)
  : ConstraintLDS( s, s->variables.stack_, s->length )
{
  int i;
  lb_threshold = ub_threshold = -1;
  s->binds(mindiff);
  s->binds(maxdiff);
  mindiff.setValue( 0 );
  maxdiff.setValue( arity );

  s->binds( maybeequal );
  s->binds( cannotequal );
  maybeequal.setValue( 0, arity-1 );
  cannotequal.setValue( 0, arity-1 );
  

  ideal = new int[arity];

  for(i=0; i<arity; ++i) {
    ideal[i] = scope[i]->branch->getBest();
    if( scope[i]->isGround() ) {
      --maxdiff;
      cannotequal.erase( i );
    }
  }
}

ConstraintHamming::ConstraintHamming(Solver *s, const int n, BuildObject** exp, int *sol)
  : ConstraintLDS( s, s->variables.stack_, s->length )
{
  int i=0, j=0;//, ideal[arity];
  BuildObject *bvar;
  VariableInt *temp;
  while( i<n ) {
    bvar = exp[i];

//     std::cout << sol[i] << "  -- ";
//     bvar->print(std::cout);
//     std::cout << std::endl;

    if( bvar->isSearchable() &&
	((temp = bvar->getVariable()) != NULL) ) {
      
      //ideal[j++] = sol[i];
      sol[j++] = sol[i];
      
    } 
    ++i;
  }
  
  assert(j==arity);
  
  initialise(s,sol);
  
}

ConstraintHamming::ConstraintHamming(Solver *s, const int *sol)
  : ConstraintLDS( s, s->variables.stack_, s->length )
{
  initialise(s,sol);
}

void ConstraintHamming::initialise(Solver *s, const int *sol)
{
  int i;
  lb_threshold = ub_threshold = -1;
  //threshold = -1;
  s->binds(mindiff);
  s->binds(maxdiff);
  mindiff.setValue( 0 );
  maxdiff.setValue( arity );

  s->binds( maybeequal );
  s->binds( cannotequal );
  maybeequal.setValue( 0, arity-1 );
  cannotequal.setValue( 0, arity-1 );
  

  ideal = new int[arity];

  //std::cout << "ideal: ";
  for(i=0; i<arity; ++i) {
    ideal[i] = sol[i];
    //std::cout << ideal[i] << " "; 
    if( !scope[i]->contain( ideal[i] ) ) {
      ++mindiff;
      maybeequal.erase( i );
    } else if( scope[i]->isGround() ) {
      --maxdiff;
      cannotequal.erase( i );
    }
  }
  //std::cout << std::endl;
}

ConstraintHamming::~ConstraintHamming()
{
  delete [] ideal;
}

int ConstraintHamming::check( const int* s ) const 
{
  int i, nd=0;
  for(i=0; i<arity; ++i)
    nd += ( s[i] != ideal[i] );
  return (nd < lb_threshold || nd > ub_threshold);
}

bool ConstraintHamming::propagate(const int changedIdx, const int e)
{
  bool consistent = true;
  int i;//, idx=changedIdx;//v->id;

//   for(i=0; i<arity; ++i) {
//     if(scope[i]->equal(ideal[i]))
//       std::cout << "0";
//     else if(!scope[i]->contain(ideal[i]))
//       std::cout << "X";
//     else
//       std::cout << "?";
//   }
//   std::cout << std::endl;


  if( maybeequal.member( changedIdx ) )
    {
      if( !scope[changedIdx]->contain( ideal[changedIdx] ) ) {
	maybeequal.erase( changedIdx );
	++mindiff;
      } else if( e & Constraint::VALUETRIGGER ) {
	--maxdiff;
	cannotequal.erase( changedIdx );
      }
    }

  if( mindiff > ub_threshold || maxdiff < lb_threshold ) consistent = false;
  else {
    if( maxdiff == lb_threshold )
      {
	for(i=0; consistent && i<arity; ++i)
	  if( cannotequal.member( i ) )
	    consistent = scope[i]->remove( ideal[i] );
      }
    else if( mindiff == ub_threshold )
      {
	for(i=0; consistent && i<arity; ++i)
	  if( maybeequal.member( i ) ) 
	    consistent = scope[i]->setDomain( ideal[i] );
      }
  }

//   for(i=0; i<arity; ++i) {
//     if(!cannotequal.member(i))
//       std::cout << "0";
//     else if(!maybeequal.member(i))
//       std::cout << "X";
//     else
//       std::cout << "?";
//   }
//   std::cout << "[" << (mindiff) << ","
// 	    << (maxdiff) << "] ["
// 	    << (lb_threshold) << "," 
// 	    << (ub_threshold) << "]"
// 	    << std::endl
// 	    << std::endl;

  return consistent;
}

void ConstraintHamming::print(std::ostream& o) const 
{
  o << "Hamming-" ;
  if(lb_threshold == ub_threshold)
    o << lb_threshold ;
  else 
    o << "[" << lb_threshold << ".." << ub_threshold << "]";
  o << ": ";
  for(int i=0; i<arity; ++i)
    o << std::setw(4) << ideal[i] << " ";

  o << std::endl;
  o << "maybeequal: ";
  maybeequal.print( o );
  o << std::endl;

  o << "cannotequal: ";
  cannotequal.print( o );
  o << std::endl;

  o << "[" << mindiff << ".." << maxdiff << "]" ;

}


ConstraintDistance::ConstraintDistance(Solver* s)
  : ConstraintLDS( s, s->sequence, s->length )
{
  int i, m;
  maxthreshold=0;

  mindiff = new ReversibleNum<int>[arity+1];
  maxdiff = new ReversibleNum<int>[arity+1]; 
  distance_to_val = new int*[arity];
  val_to_distance = new int*[arity];

  for(i=0; i<arity; ++i) {

    m = (scope[i]->maxCapacity() - scope[i]->minCapacity() + 1);
    distance_to_val[i] = new int[scope[i]->domsize()]; 
    val_to_distance[i] = new int[m]; 
    val_to_distance[i] -= scope[i]->minCapacity(); 
    scope[i]->branch->getOrder( distance_to_val[i], val_to_distance[i] );
    maxthreshold += scope[i]->domsize()-1;
    
    s->binds( mindiff[i] );
    s->binds( maxdiff[i] );
    mindiff[i].setValue( 0 );
    maxdiff[i].setValue( scope[i]->domsize()-1 );

    scope[i]->print( std::cout );
    for(int j=scope[i]->minCapacity(); j<scope[i]->maxCapacity(); ++j)
      std::cout << " " << distance_to_val[i][j] ;
    std::cout << " " << maxthreshold ;
    std::cout << std::endl;

  }

  lb_threshold = -1;
  ub_threshold = -1;
  s->binds( mindiff[arity] );
  s->binds( maxdiff[arity] );
  mindiff[arity].setValue( 0 );
  maxdiff[arity].setValue( maxthreshold );  

  std::cout << "max thresh: " << maxthreshold << std::endl;  

}

ConstraintDistance::~ConstraintDistance()
{
  for(int i=0; i<arity; ++i) {
    delete [] distance_to_val[i];
    val_to_distance[i] += scope[i]->minCapacity();
    delete [] val_to_distance[i];
  }
  delete [] distance_to_val;
  delete [] val_to_distance;
}

int ConstraintDistance::check( const int* s ) const {
  int i, nd=0;
  for(i=0; i<arity; ++i)
    nd += ( val_to_distance[i][s[i]] );
  return (nd < lb_threshold || nd > ub_threshold);
}

bool ConstraintDistance::propagate(const int changedIdx, const int e)
{
  bool consistent = true;
  int i = mindiff[changedIdx], j = maxdiff[changedIdx], k;

//   for(k=0; k<arity; ++k)
//     {
//       std::cout << "\t[" << mindiff[k] << ".." << maxdiff[k] << "]  " ;
//       scope[k]->print( std::cout );
//       if( k == idx )
// 	std::cout << "*";
//       std::cout << std::endl;
//     }
//   std::cout << "\t[" << mindiff[arity] << ".." << maxdiff[arity] << "]" << std::endl << std::endl;


  if( e & Constraint::VALUETRIGGER ) {
    k = val_to_distance[changedIdx][scope[changedIdx]->value()];

    mindiff[changedIdx] = k;
    maxdiff[changedIdx] = k;
    mindiff[arity] += (k-i);
    maxdiff[arity] += (k-j);

  } else {
    
    if( !scope[changedIdx]->fastContain( distance_to_val[changedIdx][i] ) ) {
      k=i;
      do { ++k;
      } while( !scope[changedIdx]->fastContain( distance_to_val[changedIdx][k] ) );

      mindiff[arity] += (k-i);
      mindiff[changedIdx] = k;
    }
    
    if( !scope[changedIdx]->fastContain( distance_to_val[changedIdx][j] ) ) {
      k=j;
      do --k;
      while( !scope[changedIdx]->fastContain( distance_to_val[changedIdx][k] ) );
      maxdiff[arity] += (k-j);
      maxdiff[changedIdx] = k;
    }
  }

  int dmax = maxdiff[arity];
  int dmin = mindiff[arity];

  if( dmin > ub_threshold || dmax < lb_threshold ) consistent = false;
  else {
    for(k=0; k<arity; ++k)
      {

	if( ( dmin - mindiff[k] + maxdiff[k]) > ub_threshold )
	  {
	    i = (ub_threshold - dmin + mindiff[k]);
	    j = maxdiff[k];
	    while( ++i <= j && consistent )
	      consistent = scope[k]->remove( distance_to_val[k][i] );
	  }

	if( ( dmax - maxdiff[k] + mindiff[k]) < lb_threshold )
	  {
	    i = (maxdiff[k] + lb_threshold - dmax);
	    j = mindiff[k];
	    while( --i >= j && consistent )
	      consistent = scope[k]->remove( distance_to_val[k][i] );
	  }

      }
  }

//   std::cout << "\t[" << mindiff[arity] << ".." << maxdiff[arity] << "]" << std::endl ;
//   for(k=0; k<arity; ++k)
//     {
//       std::cout << "\t[" << mindiff[k] << ".." << maxdiff[k] << "]  " ;
//       scope[k]->print( std::cout );
//       if( k == changedIdx )
// 	std::cout << "*";
//       std::cout << std::endl;
//     }
//   std::cout << std::endl << consistent << " ======================" << std::endl << std::endl;
    
  return consistent;
}

void ConstraintDistance::print(std::ostream& o) const 
{
  std::cout << "Distance-[" << lb_threshold  << ".." << ub_threshold << "]";
}




/******************************************************************************
 File: alldiff.cpp

 Implementation of the algorithm for bounds consistency of the
 generalized cardinality constraint described in:

	C.-G. Quimper, P. van Beek, A. Lopez-Ortiz, A. Golynski, and
	S.B. Sadjad. An efficient bounds consistency algorithm for the
	global cardinality constraint. CP-2003.

 By: Claude-Guy Quimper
 ******************************************************************************/

ConstraintGlobalCardinality::ConstraintGlobalCardinality(Solver *s,
							 VariableInt **v,
							 const int n,
							 const int firstDomainValue, 
							 const int lastDomainValue,
							 const int* minOccurrences,
							 const int* maxOccurrences )
  : Constraint(s, v, n, RANGETRIGGER )
{
  delayed = true;

  int i, range;

  range = lastDomainValue - firstDomainValue + 1;

  lastLevel = -1;

  iv        = (Interval  *)calloc(arity, sizeof(Interval  ));
  minsorted = (Interval **)calloc(arity, sizeof(Interval *));
  maxsorted = (Interval **)calloc(arity, sizeof(Interval *));
  bounds    = (int *)calloc(2*arity+2, sizeof(int));

  for( i = 0; i < arity; i++ ) {
    minsorted[i] = maxsorted[i] = &iv[i];
  }

  t = (int *)calloc(2*arity+2, sizeof(int));
  d = (int *)calloc(2*arity+2, sizeof(int));
  h = (int *)calloc(2*arity+2, sizeof(int));

  stableInterval      = (int *)calloc(2*arity+2, sizeof(int));
  potentialStableSets = (int *)calloc(2*arity+2, sizeof(int));
  newMin              = (int *)calloc(  arity,   sizeof(int));

  l = initializePartialSum(firstDomainValue, range, minOccurrences);
  u = initializePartialSum(firstDomainValue, range, maxOccurrences);
}


ConstraintGlobalCardinality::~ConstraintGlobalCardinality()
{
  free(bounds);
  free(maxsorted);
  free(minsorted);
  free(iv);
  free(h);
  free(d);
  free(t);
  free(newMin);
  free(potentialStableSets);
  free(stableInterval);
  destroyPartialSum(u);
  destroyPartialSum(l);
}

void
ConstraintGlobalCardinality::sortit()
{
  int i,j,nb,min,max,last;

  sortmin(minsorted, arity);
  sortmax(maxsorted, arity);

  min = minsorted[0]->min;
  max = maxsorted[0]->max + 1;
  //MODIFIED: bounds[0] = last = min-2;
  bounds[0] = last = l->firstValue + 1;

  for (i=j=nb=0;;) { // merge minsorted[] and maxsorted[] into bounds[]
    if (i<arity && min<=max) {	// make sure minsorted exhausted first
      if (min != last)
        bounds[++nb] = last = min;
      minsorted[i]->minrank = nb;
      if (++i < arity)
        min = minsorted[i]->min;
    } else {
      if (max != last)
         bounds[++nb] = last = max;
      maxsorted[j]->maxrank = nb;
      if (++j == arity) break;
      max = maxsorted[j]->max + 1;
    }
  }
  ConstraintGlobalCardinality::nb = nb;
  //MODIFIED: bounds[nb+1] = bounds[nb] + 2;
  bounds[nb+1] = u->lastValue + 1;
}


/* 
 * Shrink the lower bounds for the max occurrences problem.
 */
int
ConstraintGlobalCardinality::filterLowerMax()
{
  int i,j,w,x,y,z;
  int changes = 0;

  for (i=1; i<=nb+1; i++) {
    t[i] = h[i] = i-1;
    d[i] = sum(u, bounds[i-1], bounds[i]-1);
  }
  for (i=0; i<arity; i++) { // visit intervals in increasing max order
    // get interval bounds
    x = maxsorted[i]->minrank; y = maxsorted[i]->maxrank;
    j = t[z = pathmax(t, x+1)];
    if (--d[z] == 0) {
      t[z = pathmax(t, t[z]=z+1)] = j;
    }
    pathset(t, x+1, z, z);
    if (d[z] < sum(u, bounds[y], bounds[z] - 1)) {
      return INCONSISTENT; // no solution
    }
    if (h[x] > x) {
      maxsorted[i]->min = bounds[w = pathmax(h, h[x])];
      pathset(h, x, w, w);
      changes = 1;
    }
    if (d[z] == sum(u, bounds[y], bounds[z] - 1)) {
      pathset(h, h[y], j-1, y); // mark hall interval
      h[y] = j-1; //("hall interval [%d,%d)\n",bounds[j],bounds[y]);
    }
  }
  if( changes )
    return CHANGES;
  else
    return NO_CHANGES;
}

/*
 * Shrink the upper bounds for the max occurrences problem.
 */
int
ConstraintGlobalCardinality::filterUpperMax()
{
  // Assertion: filterLowerMax returns true
  int i,j,w,x,y,z;
  int changes = 0;

  for (i=0; i<=nb; i++) {
    d[i] = sum(u, bounds[i], bounds[t[i]=h[i]=i+1]-1);
  }
  for (i=arity; --i>=0; ) { // visit intervals in decreasing min order
    // get interval bounds
    x = minsorted[i]->maxrank; y = minsorted[i]->minrank;
    j = t[z = pathmin(t, x-1)];
    if (--d[z] == 0) {
      t[z = pathmin(t, t[z]=z-1)] = j;
    }
    pathset(t, x-1, z, z);
    if (d[z] < sum(u, bounds[z], bounds[y]-1)) {
      return INCONSISTENT; // no solution
    }
    if (h[x] < x) {
      minsorted[i]->max = bounds[w = pathmin(h, h[x])] - 1;
      pathset(h, x, w, w);
      changes = 1;
    }
    if (d[z] == sum(u, bounds[z], bounds[y]-1)) {
      pathset(h, h[y], j+1, y);
      h[y] = j+1;
    }
  }
  if( changes )
    return CHANGES;
  else
    return NO_CHANGES;
}


/*
 * Shrink the lower bounds for the min occurrences problem.
 * called as: filterLowerMin(t, d, h, stableInterval, potentialStableSets, newMin);
 */
int
ConstraintGlobalCardinality::filterLowerMin(int *tl, int *c,
					    int* stableAndUnstableSets,
					    int* stableInterval,
					    int* potentialStableSets,
					    int* newMin)
{
  int i,j,w,x,y,z,v;
  int changes = 0;

  for (w=i=nb+1; i>0; i--) {
    //c[i] = sum(l, bounds[potentialStableSets[i]=stableInterval[i]=i-1], bounds[i]-1);
    potentialStableSets[i]=stableInterval[i]=i-1;
    c[i] = sum(l, bounds[i-1], bounds[i]-1);
    // If the capacity between both bounds is zero, we have
    // an unstable set between these two bounds.
    if (c[i] == 0) {
      stableAndUnstableSets[i-1] = w;
    }
    else {
      w = stableAndUnstableSets[w] = i - 1;
    }
  }

  for (i = w = nb + 1; i >= 0; i--) {
    if (c[i] == 0)
      tl[i] = w;
    else
      w = tl[w] = i;
  }

  for (i = 0; i < arity; i++) { // visit intervals in increasing max order
    // Get interval bounds
    x = maxsorted[i]->minrank; y = maxsorted[i]->maxrank;
    j = tl[z = pathmax(tl, x+1)];
    if (z != x+1) {
      // if bounds[z] - 1 belongs to a stable set,
      // [bounds[x], bounds[z]) is a sub set of this stable set
      v = potentialStableSets[w = pathmax(potentialStableSets, x + 1)];
      pathset(potentialStableSets, x+1, w, w); // path compression
      w = y < z ? y : z;
      pathset(potentialStableSets, potentialStableSets[w], v , w);
      potentialStableSets[w] = v;
    }

    if (c[z] <= sum(l, bounds[y], bounds[z] - 1)) {
      // (potentialStableSets[y], y] is a stable set
      w = pathmax(stableInterval, potentialStableSets[y]);
      pathset(stableInterval, potentialStableSets[y], w, w); // Path compression
      pathset(stableInterval, stableInterval[y], v=stableInterval[w], y);
      stableInterval[y] = v;
    }
    else {
      // Decrease the capacity between the two bounds
      if (--c[z] == 0) {
	tl[z = pathmax(tl, tl[z]=z+1)] = j;
      }

      // If the lower bound belongs to an unstable or a stable set,
      // remind the new value we might assigned to the lower bound
      // in case the variable doesn't belong to a stable set.
      if (stableAndUnstableSets[x] > x) {
	w = newMin[i] = pathmax(stableAndUnstableSets, x);
	pathset(stableAndUnstableSets, x, w, w); // path compression
      }
      else {
	newMin[i] = x; // Do not shrink the variable
      }

      // If an unstable set is discovered
      if (c[z] == sum(l, bounds[y], bounds[z] - 1)) {
 	if (stableAndUnstableSets[y] > y) // Consider stable and unstable sets beyong y
 	  y = stableAndUnstableSets[y]; // Equivalent to pathmax since the path is fully compressed
	pathset(stableAndUnstableSets, stableAndUnstableSets[y], j-1, y); // mark the new unstable set
	stableAndUnstableSets[y] = j-1;
      }
    }
    pathset(tl, x+1, z, z); // path compression
  }

  // If there is a failure set
  if (stableAndUnstableSets[nb] != 0) {
    return INCONSISTENT; // no solution
  }

  // Perform path compression over all elements in
  // the stable interval data structure. This data
  // structure will no longer be modified and will be
  // accessed n or 2n times. Therefore, we can afford
  // a linear time compression.
  for (i = nb+1; i > 0; i--) {
    if (stableInterval[i] > i)
      stableInterval[i] = w;
    else
      w = i;
  }

  // For all variables that are not a subset of a stable set, shrink the lower bound
  for (i=arity-1; i>=0; i--) {
    x = maxsorted[i]->minrank; y = maxsorted[i]->maxrank;
    if ((stableInterval[x] <= x) || (y > stableInterval[x])) {
      maxsorted[i]->min = skipNonNullElementsRight(l, bounds[newMin[i]]);
      changes = 1;
    }
  }

  if( changes )
    return CHANGES;
  else
    return NO_CHANGES;
}


/*
 * Shrink the upper bounds for the min occurrences problem.
 * called as: filterUpperMin(t, d, h, stableInterval, newMin);
 */
int 
ConstraintGlobalCardinality::filterUpperMin(int *tl, int *c,
					    int* stableAndUnstableSets,
					    int* stableInterval,
					    int* newMax)
{
  // ASSERTION: filterLowerMin returns true
  int i,j,w,x,y,z;
  int changes = 0;

  for (w=i=0; i<=nb; i++) {
    //    d[i] = bounds[t[i]=h[i]=i+1] - bounds[i];
    c[i] = sum(l, bounds[i], bounds[i+1]-1);
    if (c[i] == 0)
      tl[i]=w;
    else
      w=tl[w]=i;
  }
  tl[w]=i;
  for (i = 1, w = 0; i<=nb; i++) {
    if (c[i-1] == 0)
      stableAndUnstableSets[i] = w;
    else
      w = stableAndUnstableSets[w] = i;
  }
  stableAndUnstableSets[w] = i;

  for (i=arity; --i>=0; ) { // visit intervals in decreasing min order
    // Get interval bounds
    x = minsorted[i]->maxrank; y = minsorted[i]->minrank;

    // Solve the lower bound problem
    j = tl[z = pathmin(tl, x-1)];

    // If the variable is not in a discovered stable set
    // Possible optimization: Use the array stableInterval to perform this test
    if (c[z] > sum(l, bounds[z], bounds[y]-1)) {
      if (--c[z] == 0) {
	tl[z = pathmin(tl, tl[z]=z-1)] = j;
      }
      if (stableAndUnstableSets[x] < x) {
	newMax[i] = w = pathmin(stableAndUnstableSets, stableAndUnstableSets[x]);
	pathset(stableAndUnstableSets, x, w, w); // path compression
      }
      else {
	newMax[i] = x;
      }
      if (c[z] == sum(l, bounds[z], bounds[y]-1)) {
	if (stableAndUnstableSets[y] < y) {
	  y = stableAndUnstableSets[y];
        }
	pathset(stableAndUnstableSets, stableAndUnstableSets[y], j+1, y);
	stableAndUnstableSets[y] = j+1;
      }
    }
    pathset(tl, x-1, z, z);
  }

  // For all variables that are not subsets of a stable set, shrink the lower bound
  for (i=arity-1; i>=0; i--) {
    x = minsorted[i]->minrank; y = minsorted[i]->maxrank;
    if ((stableInterval[x] <= x) || (y > stableInterval[x]))
      minsorted[i]->max = skipNonNullElementsLeft(l, bounds[newMax[i]]-1);
      changes = 1;
  }

  if( changes )
    return CHANGES;
  else
    return NO_CHANGES;
}


bool ConstraintGlobalCardinality::propagate(const int changedIdx, const int e)
{
  return propagate();
}

bool ConstraintGlobalCardinality::propagate()
{
  int statusLower, statusUpper;
  int statusLowerMin=NOVAL, statusUpperMin;
  int statusLowerMax=-NOVAL, statusUpperMax;
  int i, a, b;
  int dl, du;

  a = 0;
  b = arity;
  //  a = _vars.getRangeIndexMin();
  //  b = _vars.getRangeIndexMax();

/*
  if( _prop == WithValueRemoval && (a == (b-1)) && _vars[a].isBound() ) {
    return;
  }
*/


  //if( lastLevel != scope[0]->solver->level ) {
  //if( lastLevel != solver->level ) {
  if( lastLevel != *level ) {
    // not incremental
    statusLower = CHANGES;
    statusUpper = CHANGES;

    i = arity;
    while( i-- ) {
      iv[i].min = scope[i]->min();
      iv[i].max = scope[i]->max();
    }

  }
  else {
    // incremental
    statusLower = NO_CHANGES;
    statusUpper = NO_CHANGES;
    for( i = a; i < b; i++ ) {
      dl = iv[i].min;
      du = iv[i].max;
      iv[i].min = scope[i]->min();
      iv[i].max = scope[i]->max();
      if( dl != iv[i].min ) statusLower = CHANGES;
      if( du != iv[i].max ) statusUpper = CHANGES;
    }
  }

  //lastLevel = scope[0]->solver->level;
  //lastLevel = solver->level;
  lastLevel = *level;

  if( statusLower == NO_CHANGES && statusUpper == NO_CHANGES ) {
    return true;
  }

  sortit();

  // The variable domains must be inside the domain defined by
  // the lower bounds (l) and the upper bounds (u).
  //assert(minValue(l) == minValue(u));
  //assert(maxValue(l) == maxValue(u));
  //assert(minValue(l) <= minsorted[0]->min);
  //assert(maxsorted[n-1]->max <= maxValue(u));

  // Checks if there are values that must be assigned before the
  // smallest interval or after the last interval. If this is
  // the case, there is no solution to the problem
  // This is not an optimization since
  // filterLower{Min,Max} and
  // filterUpper{Min,Max} do not check for this case.
  if ((sum(l, minValue(l), minsorted[0]->min - 1) > 0) ||
      (sum(l, maxsorted[arity-1]->max + 1, maxValue(l)) > 0)) {
    //_vars.getManager().fail();
    return false;
  }

  statusLowerMax = filterLowerMax();
  if( statusLowerMax != INCONSISTENT ) {
    statusLowerMin = filterLowerMin(t, d, h,
			stableInterval, potentialStableSets, newMin);
  }

  if( (statusLowerMax == INCONSISTENT) || (statusLowerMin == INCONSISTENT) ) {
    //_vars.getManager().fail();
    return false;
  }
  else {

    statusUpperMax = filterUpperMax();
    statusUpperMin = filterUpperMin(t, d, h, stableInterval, newMin);

    if( (statusLowerMax == CHANGES) || (statusLowerMin == CHANGES) ||
        (statusUpperMax == CHANGES) || (statusUpperMin == CHANGES) ) {
      i = arity;
      while( i-- )
	if( !(scope[i]->setMin( iv[i].min ) && scope[i]->setMax( iv[i].max )) )
	  return false;
    } // if
  } // else

  return true;
}

// Create a partial sum data structure adapted to the
// filterLower{Min,Max} and filterUpper{Min,Max} functions.
// Two elements before and after the element list will be added
// with a weight of 1. 
partialSum*
ConstraintGlobalCardinality::initializePartialSum(const int firstValue, int count, 
						  const int* elements)
{
  int i,j;
  partialSum* res = new partialSum;
  res->sum = new int[count+1+2+2];
  res->firstValue = firstValue - 3; // We add three elements at the beginning
  res->lastValue = firstValue + count + 1;
  res->sum[0] = 0;
  res->sum[1] = 1;
  res->sum[2] = 2;
  for (i = 2; i < count+2; i++) {
    res->sum[i+1] = res->sum[i] + elements[i-2];
  }
  res->sum[i+1] = res->sum[i] + 1;
  res->sum[i+2] = res->sum[i+1] + 1;

  res->ds = new int[count+1+2+2];

  for (j=(i=count+3)+1; i > 0;) {
    while (res->sum[i] == res->sum[i-1])
      res->ds[i--]=j;
    j=res->ds[j]=i--;
  }
  res->ds[j]=0;
  return res;
}

void
ConstraintGlobalCardinality::destroyPartialSum(partialSum *p)
{
  delete p->ds;
  delete p->sum;
  delete p;
}

int
ConstraintGlobalCardinality::sum(partialSum *p, int from, int to)
{
  if (from <= to) {
    //assert((p->firstValue <= from) && (to <= p->lastValue));
    return p->sum[to - p->firstValue] - p->sum[from - p->firstValue - 1];
  }
  else {
    //assert((p->firstValue <= to) && (from <= p->lastValue));
    return p->sum[to - p->firstValue - 1] - p->sum[from - p->firstValue];
  }
}

int
ConstraintGlobalCardinality::minValue(partialSum *p) {
  return p->firstValue + 3;
}

int
ConstraintGlobalCardinality::maxValue(partialSum *p) {
  return p->lastValue - 2;
}

int
ConstraintGlobalCardinality::skipNonNullElementsRight(partialSum *p, int value) {
  value -= p->firstValue;
  return (p->ds[value] < value ? value : p->ds[value]) + p->firstValue;
}

int
ConstraintGlobalCardinality::skipNonNullElementsLeft(partialSum *p, int value) {
  value -= p->firstValue;
  return (p->ds[value] > value ? p->ds[p->ds[value]] : value) + p->firstValue;
}


int ConstraintGlobalCardinality::check( const int* s ) const 
{
//   int i=arity, j;
//   while(--i) {
//     j=i;
//     while(j--)
//       if( s[i] == s[j] ) return 1;
//   }
  return 0; 
}

void ConstraintGlobalCardinality::print(std::ostream& o) const 
{    
    o << "Gcc(";
    for(int i=0; i<arity-1; ++i) {
      scope[i]->printshort( o );
      o << " " ;
    }
    scope[arity-1]->printshort( o );
    o << ", [" << (l->firstValue+3) << ".." 
      << (l->lastValue-2) << "] >={";
    o << " }, <={";
    o << " })";
}

/*
 *  End of user defined propagator for enforcing bounds consistency
 *=================================================================*/


ConstraintSlideDecomp::ConstraintSlideDecomp(Solver *s, VariableInt** v, 
					     const int n, const int m,
					     const int nv, const int sw, const int K,
					     int* cft, int* t2d, int* d2t, int* exp) 
{
  nbVals = nv; 
  sequenceWidth = sw;
  conflicts = cft;
  int2ind = t2d;
  ind2int = d2t;
  exposant = exp;

  int i, j, k, nSlideVars = m-sequenceWidth+2;
  VariableInt *slideVars[n*nSlideVars];
  for(i=0; i<n; ++i)
    for(j=0; j<nSlideVars; ++j)
      slideVars[i*nSlideVars+j] = new VariableBit( s, 0, K-1 );
  
  VariableInt *scp[sequenceWidth];
  for(i=0; i<n; ++i) {
    for(j=0; j<nSlideVars; ++j) {
      for(k=0; k<sequenceWidth-1; ++k)
	scp[k] = v[i*m+j+k];
      scp[sequenceWidth-1] = slideVars[i*nSlideVars+j];
      new PredicateTuple( s, scp,
			  sequenceWidth,
			  nbVals,
			  int2ind,
			  ind2int,
			  exposant );
    }
    for(j=0; j<nSlideVars-1; ++j) {
      scp[0] = slideVars[i*nSlideVars+j];
      scp[1] = slideVars[i*nSlideVars+j+1];
      new ConstraintSlide( s, scp, this );
    }
  }
}

ConstraintSlideDecomp::~ConstraintSlideDecomp()
{
  delete [] conflicts;
  delete [] int2ind;
  delete [] ind2int;
  delete [] exposant;
}

ConstraintSlide::ConstraintSlide(Solver *s, VariableInt** v, ConstraintSlideDecomp* c) 
  : Constraint(s, v, 2) 
{
  nbVals = c->nbVals; 
  conflicts = c->conflicts;
  int2ind = c->int2ind;
  ind2int = c->ind2int;
  exposant = c->exposant;
  sequenceWidth = c->sequenceWidth;
}

ConstraintSlide::~ConstraintSlide()
{
}

int ConstraintSlide::check( const int* s ) const
{
  return conflicts[appendLeft( s[0], s[1] )];
}

bool ConstraintSlide::propagate(const int changedIdx, const int e) 
{
  int // vnxt, value,
    tuple, support;
  bool is_supported, consistent=true;;
  
  if( !changedIdx ) {

//     value = scope[1]->first();
//     do {
//       vnxt = scope[1]->getNext( value );

//       is_supported = false;
//       support = -1;
//       while( !is_supported && ++support < nbVals ) {
// 	tuple = appendLeft(value,support);
// 	is_supported = (
// 			!conflicts[tuple]
// 			&&
// 			scope[0]->contain( int2ind[tuple%exposant[sequenceWidth-1]] )
// 			);
//       }
//       if( !is_supported ) 
// 	consistent = scope[1]->remove( value );
//       vnxt = value;
//     }
//     while( consistent && scope[1]->good(value) );

    DomainIterator *valit = scope[1]->begin();
    do {
      is_supported = false;
      support = -1;
      while( !is_supported && ++support < nbVals ) {
	tuple = appendLeft(*valit, support);
	is_supported = (
			!conflicts[tuple]
			&&
			scope[0]->contain( int2ind[tuple%exposant[sequenceWidth-1]] )
			);
      }
      if( !is_supported ) 
	consistent = scope[1]->remove( *valit );
    }
    while( consistent && valit->next() );

  } else {
    
//     value = scope[0]->first();
//     do {
//       vnxt = scope[0]->getNext( value );

// 	is_supported = false;
// 	support = -1;
// 	while( !is_supported && ++support < nbVals ) {
// 	  tuple = appendRight(value,support);
// 	  is_supported = (
// 			  !conflicts[tuple]
// 			  &&
// 			  scope[1]->contain( int2ind[tuple/nbVals] )
// 			  );
// 	}
// 	if( !is_supported ) 
// 	  consistent = scope[0]->remove( value );
// 	vnxt = value;
//     }
//     while( consistent && scope[0]->good(value) );

    DomainIterator *valit = scope[0]->begin();
    do {
	is_supported = false;
	support = -1;
	while( !is_supported && ++support < nbVals ) {
	  tuple = appendRight(*valit,support);
	  is_supported = (
			  !conflicts[tuple]
			  &&
			  scope[1]->contain( int2ind[tuple/nbVals] )
			  );
	}
	if( !is_supported ) 
	  consistent = scope[0]->remove( *valit );
    }
    while( consistent && valit->next() );
    
  }

  return consistent;
}

void ConstraintSlide::print(std::ostream& o) const 
{    
  o << "(";
  scope[0]->printshort( o );
  o << " >> ";
  scope[1]->printshort( o );
  o << ")";
}


PredicateTuple::PredicateTuple(Solver *s, VariableInt** v, 
			       const int n, const int nv,
			       int* t2d, int* d2t, int* exp) 
  : Constraint(s, v, n) 
{
  nbVals = nv; 
  int2ind = t2d;
  ind2int = d2t;
  exposant = exp;
}

PredicateTuple::~PredicateTuple()
{
}

int PredicateTuple::check( const int* s ) const
{
  int i, tuple=0, n=arity-1;
  for(i=0; i<n; ++i)
    tuple += (s[i] * exposant[i]);
  return ( s[n] != tuple );
}

bool PredicateTuple::propagate(const int changedIdx, const int e) 
{
  int // value, vnxt,
    i, n=arity-1;
  bool consistent = true;

  if( changedIdx == n ) {

    for(i=0; i<n; ++i)
      vals[i].clear();

//     value = scope[n]->first();
//     do for(i=0; i<n; ++i) 
// 	 vals[i].insert((ind2int[value]/exposant[i])%nbVals);
//     while( scope[n]->setNext( value ) );

    DomainIterator *valit = scope[n]->begin();
    do for(i=0; i<n; ++i) 
	 vals[i].insert((ind2int[*valit]/exposant[i])%nbVals);
    while( valit->next() );


    for(i=0; consistent && i<n; ++i) 
      consistent = scope[i]->setDomain( vals[i] );
 
  } else {

//     value = scope[n]->first();
//     do {	
//       vnxt = scope[n]->getNext( value );
//       if( !(v->contain((ind2int[value]/exposant[i])%nbVals)) )
// 	consistent = scope[n]->remove( value );
//       value = vnxt;
//     }
//     while( consistent && scope[n]->good( value ) );

    DomainIterator *valit = scope[n]->begin();
    do if( !(scope[changedIdx]->contain((ind2int[*valit]/exposant[changedIdx])%nbVals)) )
	 consistent = scope[n]->remove( *valit );
    while( consistent && valit->next() );

  }

  return consistent;
}

void PredicateTuple::print(std::ostream& o) const 
{    
  o << "(";
  for(int i=0; i<arity-2; ++i) {
    scope[i]->printshort( o );
    o << ", ";
  }
  scope[arity-1]->printshort( o );
  o << ") = " ;
  scope[arity-1]->printshort( o );
}

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////


#ifdef _CPLEX

CPlexConstraint::CPlexConstraint(Solver *s, VariableInt** v, const int n) 
  : Constraint(s, v, n )
{
  delayed = true;

  int i, j, k=0, l;

  epsilon = 0.000001;
  lb=NOVAL;
  ub=-1*NOVAL;

  // the structs asgn2index, index2valasgn, index2varasgn
  // -- values2index and index2values are initialized
  //
  // asgn2index[i][j] is the num of the column corresponding to
  // scope[i] := j
  // index2valasgn[k] is the value involved in the column k;
  // index2varasgn[k] is the variable involved in the column k;
  //
  // values2index[i] is the rank of value i
  // index2values[j] is the vale at rank j

  nRows=0;
  nColumns=0;
  asgn2index = new int*[arity-1];
  lbs = new int[arity-1];

  for(i=0; i<arity-1; ++i) {
    lbs[i] = scope[i]->minCapacity();
    nColumns += scope[i]->domsize();
    if(lb > lbs[i]) lb = lbs[i];
    if(ub < scope[i]->max()) ub = scope[i]->max();

    asgn2index[i] = new int[scope[i]->maxCapacity() - lbs[i] + 1];
    asgn2index[i] -= lbs[i];
  }

  index2valasgn = new int[nColumns];
  index2varasgn = new int[nColumns];

  nvals = ub-lb+1;
  index2values = new int[nvals];
  values2index = new int[nvals];
  values2index-=lb;
  std::fill(values2index+lb, values2index+ub+1, NOVAL);

  k-0;
  nvals=0;
  for(i=0; i<arity-1; ++i) {

 //    j = scope[i]->first();
//     do {

//       asgn2index[i][j] = k;
//       index2valasgn[k] = j;
//       index2varasgn[k] = i;
//       ++k;

//       if( values2index[j] == NOVAL ) {
// 	values2index[j] = nvals;
// 	index2values[nvals] = j;
// 	++nvals;
//       }
//     } while( scope[i]->setNext(j) ); 

    DomainIterator *valit = scope[i]->begin();
    do {
      j=*valit;

      asgn2index[i][j] = k;
      index2valasgn[k] = j;
      index2varasgn[k] = i;
      ++k;

      if( values2index[j] == NOVAL ) {
	values2index[j] = nvals;
	index2values[nvals] = j;
	++nvals;
      }
    } while( valit->next() ); 

  }
  nRows = arity-1;

  /*************************
   * CPlex structs
   *************************/  

  obj     = (double *) malloc (nColumns * sizeof(double));
  rhs     = (double *) malloc (nRows    * sizeof(double));
  sense   = (char *)   malloc (nRows    * sizeof(char)); 
  matbeg  = (int *)    malloc (nColumns * sizeof(int));
  matcnt  = (int *)    malloc (nColumns * sizeof(int));   
  matind  = (int *)    malloc (nColumns * sizeof(int));   
  matval  = (double *) malloc (nColumns * sizeof(double));
  zlb     = (double *) malloc (nColumns * sizeof(double));
  zub     = (double *) malloc (nColumns * sizeof(double));
  qmatbeg = (int *)    malloc (nColumns * sizeof(int)); 
  qmatcnt = (int *)    malloc (nColumns * sizeof(int));   

  env = NULL;
  env = CPXopenCPLEX (&status);
  if ( env == NULL ) {
    char  errmsg[1024];
    fprintf (stderr, "Could not open CPLEX environment.\n");
    CPXgeterrorstring (env, status, errmsg);
    fprintf (stderr, "%s", errmsg);
  }

  status = CPXsetintparam (env, CPX_PARAM_SCRIND, CPX_ON);
  if ( status ) {
    fprintf (stderr, 
	     "Failure to turn on screen indicator, error %d.\n", status);
  }


  lp = CPXcreateprob (env, &status, "bidule");
  if ( lp == NULL ) {
    fprintf (stderr, "Failed to create problem.\n");
  }

  //////////////// CREATING THE LP /////////////////

  // maximise (not sure if it applies to the quad objective)
  objsen = -1;

  
  for(i=0; i<nColumns; ++i) { // for each columns
    // only one constraint on each var
    matcnt[i] = 1;

    // one coeff for each column
    matbeg[i] = i;

    // the index of the row is the CSP var corresponding to each CSP assignment
    j = index2varasgn[i];
    matind[i] = j;
    matval[i] = 1.0;

    // [0..1] variables
    zlb[i] = 0.0;
    zub[i] = 1.0;

    // no linear objective
    obj[i] = 0.0;
  }

  for(i=0; i<nRows; ++i) {
    // exactly one value per domain
    rhs[i] = 1.0;

    // all equalities
    sense[i] = 'E';
  }


  status = CPXcopylp (env, lp, nColumns, nRows, objsen, obj, rhs, 
		      sense, matbeg, matcnt, matind, matval,
		      zlb, zub, NULL);

  //////////////////////////////////////////////////

  x  = new double[nColumns];
  dj = new double[nColumns];
  pi    = new double[nRows];
  slack = new double[nRows];

/*************************
 *                       *
 *************************/


}

CPlexConstraint::~CPlexConstraint()
{
  values2index += lb;
  delete [] values2index;  
  delete [] index2values;  

  for(int i=0; i<arity-1; ++i) {
    asgn2index[i] += lbs[i];
    delete [] asgn2index[i];
  }
  delete [] asgn2index;
  delete [] index2valasgn;
  delete [] index2varasgn;
  delete [] lbs;
}

int CPlexConstraint::check(const int* sol) const 
{
  return 0;
}

void CPlexConstraint::updateModel()
{
}

bool CPlexConstraint::propagate(const int changedIdx, const int e)
{
}

void CPlexConstraint::print(std::ostream& o) const 
{
}

/**********************************************
 * CPlexMinDiff Predicate 
 **********************************************/

PredicateCPlexMinDiff::PredicateCPlexMinDiff(Solver *s, const int f, VariableInt** v, const int n) 
  : CPlexConstraint(s, v, n)
{ 
  filtering = f;

  int i, j, k, nRows=arity-1;

  K = arity-1;
  vars = new VariableInt*[K];
  occurrence = new int[ub-lb+1];
  occurrence -= lb;
  orderedvals = new int[ub-lb+1];

  std::fill( occurrence+lb, occurrence+ub+1, 0 );



  //////////////// CREATING THE QP /////////////////

  //std::cout << "here" << std::endl;

  int nAsgn = 0;
  DomainIterator *valit;
  for(i=0; i<nRows; ++i) {

    valit = scope[i]->begin();
    do ++occurrence[*valit];
    while( valit->next() ); 

//     j = scope[i]->first();
//     do {
//       ++occurrence[j];
//     } while( scope[i]->setNext(j) ); 
  } 

  std::cout << lb << ".." << ub << std::endl;

  for(i=lb; i<=ub; ++i) 
    if(occurrence[i])
      nAsgn += (occurrence[i] * (occurrence[i]-1));

  qmatind = (int *)    malloc (nAsgn * sizeof(int)); 
  qmatval = (double *) malloc (nAsgn * sizeof(double));   
  //qmatbeg = (int *)    malloc (nColumns * sizeof(int)); 
  //qmatcnt = (int *)    malloc (nColumns * sizeof(int)); 

  std::cout << nAsgn << std::endl;

  int str_ = 0;
  for(i=0; i<nColumns; ++i) { // for each columns
    qmatbeg[i] = str_;

    for(j=0; j<nRows; ++j) {
      if(j != index2varasgn[i] && 
	 scope[j]->contain(index2valasgn[i])) {

	qmatind[str_] = asgn2index[j][index2valasgn[i]];
	qmatval[str_++] = 0.5;
      }
    }
    
    qmatcnt[i] = (str_ - qmatbeg[i]);
    
//     std::cout << " => " << i << ": x" << index2varasgn[i] << "=" << index2valasgn[i] << " -> "
// 	 << (str_ - qmatbeg[i]) << " " << (occurrence[index2valasgn[i]]) 
// 	 << " " << qmatbeg[i] << std::endl << std::endl;
    //assert( (str - qmatbeg[i]) == (occurrence[index2valasgn[i]] - 1) );
  }

  status = CPXcopyquad (env, lp, qmatbeg, qmatcnt, qmatind, qmatval);
  if ( status ) {
    fprintf (stderr, "Failed to copy quadratic matrix.\n");
  }

  //////////////////////////////////////////////////

  //////////////// SOLVING /////////////////////////

  status = CPXqpopt (env, lp);
  if ( status ) {
    fprintf (stderr, "Failed to optimize QP.\n");
  }
  
  status = CPXsolution (env, lp, &solstat, &objval, x, pi, slack, dj);
  if ( status ) {
    fprintf (stderr, "Failed to obtain solution.\n");
  }

  status = CPXwriteprob (env, lp, "sae.lp", NULL);

  exit( 0 );

  //////////////////////////////////////////////////
  
  weight = 0;
}

PredicateCPlexMinDiff::~PredicateCPlexMinDiff()
{
  occurrence += lb;
  delete [] occurrence;
  delete [] vars;
  delete [] orderedvals;
}

int PredicateCPlexMinDiff::approximation( int& opt ) 
{
  int nbdiff = (((K)*(K-1)) >> 1);
  int i, j, nvals;//, vali;

  int nvars=K, aux;
  VariableInt *auxv;
  int vali, *firstval = orderedvals;
    
  /// count values
  std::fill( occurrence+lb, occurrence+ub+1, 0 );
  nvals=0;
  DomainIterator *valit;
  for(i=0; i<K; ++i) {
    vars[i] = scope[i];
    
    valit = scope[i]->begin();
    do if( ++occurrence[ *valit ] == 1 )
	orderedvals[nvals++] = *valit;
    while( valit->next() );

//     vali = scope[i]->first();
//     do {
//       if( ++occurrence[ vali ] == 1 )
// 	orderedvals[nvals++] = vali;
//     }
//     while( scope[i]->setNext(vali) );
  }
     
  int *lastval = orderedvals+nvals;

  while( nvars ) {
    /// order values
    for(i=1; firstval+i<lastval; ++i) {
      j = i;
      while( j && occurrence[firstval[j]] > occurrence[firstval[j-1]] )
	{
	  vali = firstval[j];
	  firstval[j] = firstval[j-1];
	  firstval[j-1] = vali;
	  --j;
	}
    }

    /// use the value with largest occurrences
    aux = occurrence[*firstval];
    aux = ((aux * (aux-1)) >> 1);
    nbdiff -= aux;

    // update the data-structures
    i = nvars;
    while( i-- ) 
      if( vars[i]->contain( *firstval ) )
	{
	  // update the array 'occurrence'
// 	  vali = vars[i]->first();
// 	  do --occurrence[ vali ];
// 	  while( vars[i]->setNext(vali) );

	  valit = vars[i]->begin();
	  do --occurrence[ *valit ];
	  while( valit->next() );
	  
	  // put all used variables at the end of the array 'vars'
	  --nvars;
	  auxv = vars[nvars];
	  vars[nvars] = vars[i];
	  vars[i] = auxv;
	}
    //--nvals;
    ++firstval;

  }

  return nbdiff;
}

bool PredicateCPlexMinDiff::propagate(const int changedIdx, const int e)
{
  bool consistent = true;
  return consistent;
}
  
void PredicateCPlexMinDiff::print(std::ostream& o) const 
{
  o << "#(";
  for(int i=0; i<K-1; ++i) {
    scope[i]->printshort(o);
    o << ", ";
  }
  scope[K-1]->printshort(o);
  o << ") = ";
  scope[K]->printshort( o );
}

#endif


#ifdef _COINLP

/**********************************************
 * Lp Constraint 
 **********************************************/

LpConstraint::LpConstraint( Solver *s ) 
  : Constraint(s)
{
  delayed = true;
}

LpConstraint::LpConstraint(Solver *s, VariableInt** v, const int n, const int t) 
  : Constraint(s, v, n )
{
  init( t );
  delayed = true;
}

void LpConstraint::init(const int t) 
{ 
  type = t;
  int i, j, k=0, l;
  epsilon = 0.000001;
  lb=NOVAL;
  ub=-1*NOVAL;

  imodel.setLogLevel( 0 );
  smodel.setLogLevel( 0 );

  // the structs asgn2index, index2valasgn, index2varasgn
  // -- values2index and index2values are initialized
  //
  // asgn2index[i][j] is the num of the column corresponding to
  // scope[i] := j
  // index2valasgn[k] is the value involved in the column k;
  // index2varasgn[k] is the variable involved in the column k;
  //
  // values2index[i] is the rank of value i
  // index2values[j] is the vale at rank j

  nRows=0;
  nColumns=0;
  asgn2index = new int*[arity-1];
  lbs = new int[arity-1];

  for(i=0; i<arity-1; ++i) {
    lbs[i] = scope[i]->minCapacity();
    nColumns += scope[i]->domsize();
    if(lb > lbs[i]) lb = lbs[i];
    if(ub < scope[i]->max()) ub = scope[i]->max();

    asgn2index[i] = new int[scope[i]->maxCapacity() - lbs[i] + 1];
    asgn2index[i] -= lbs[i];
  }

  index2valasgn = new int[nColumns];
  index2varasgn = new int[nColumns];

  nvals = ub-lb+1;
  index2values = new int[nvals];
  values2index = new int[nvals];
  values2index-=lb;
  std::fill(values2index+lb, values2index+ub+1, NOVAL);

  k-0;
  DomainIterator *valit;
  for(i=0; i<arity-1; ++i) {
    
    valit = scope[i]->begin();
    do {
      j = *valit;
      asgn2index[i][j] = k;
      index2valasgn[k] = j;
      index2varasgn[k] = i;
      ++k;

      if( values2index[j] == NOVAL ) {
	values2index[j] = nRows;
	index2values[nRows] = j;
	++nRows;
      }
    } while( valit->next() ); 

//     j = scope[i]->first();
//     do {

//       asgn2index[i][j] = k;
//       index2valasgn[k] = j;
//       index2varasgn[k] = i;
//       ++k;

//       if( values2index[j] == NOVAL ) {
// 	values2index[j] = nRows;
// 	index2values[nRows] = j;
// 	++nRows;
//       }
//     } while( scope[i]->setNext(j) ); 
  }
}

LpConstraint::~LpConstraint()
{
  values2index += lb;
  delete [] values2index;  
  delete [] index2values;  

  for(int i=0; i<arity-1; ++i) {
    asgn2index[i] += lbs[i];
    delete [] asgn2index[i];
  }
  delete [] asgn2index;
  delete [] index2valasgn;
  delete [] index2varasgn;
  delete [] lbs;
}

int LpConstraint::check(int* sol) const 
{
  return 0;
}

void LpConstraint::updateModel()
{
  int i, in;
  //double *colub;
  //  if( type == LpConstraint::INTERIOR ) {
  double *colub1 = imodel.columnUpper();
    // } else {
  double *colub2 = smodel.columnUpper();

 //  std::cout << colub1 << std::endl
//        << colub2 << std::endl << std::endl;


    // }
  for(i=0; i<nColumns; ++i) {

    //    assert( colub1[i] == colub2[i] );

    in = scope[index2varasgn[i]]->contain(index2valasgn[i]);
    if( colub1 && colub1[i] != 0 ) {
      if( !in ) {
	if( colub1 )
	  colub1[i] = 0.0;
	if( colub2 )
	  colub2[i] = 0.0;
      }
    } else if( in ) {
      if( colub1 )
	colub1[i] = 1.0;
      if( colub2 )
	colub2[i] = 1.0;
    }
  }
}

bool LpConstraint::propagate(const int changedIdx, const int e)
{

  switch( type ) {
  case INTERIOR : std::cout << "P Interior" << std::endl; break;
  case SIMPLEX_PRIMAL : std::cout << "P Simplex Primal" << std::endl; break;
  case SIMPLEX_DUAL : std::cout << "P Simplex Dual" << std::endl; break;
  case SIMPLEX_BARRIER : std::cout << "P Simplex Barrier" << std::endl; break;
  case SIMPLEX_PRIMAL_IT : std::cout << "P Simplex Primal*" << std::endl; break;
  case SIMPLEX_DUAL_IT : std::cout << "P Simplex Dual*" << std::endl; break;
  case SIMPLEX_BARRIER_IT : std::cout << "P Simplex Barrier*" << std::endl; break;
  }

  updateModel();

  int i, j, jnxt, k;
  bool consistent = true;
  DomainIterator *valit;

  switch( type ) {
  case INTERIOR : {
    
    imodel.primalDual();
    if( (consistent = imodel.primalFeasible()) ) {
      const double * solution = imodel.primalColumnSolution();
      for(i=0; consistent && i<arity; ++i) {

	valit = scope[i]->begin();
	do {
	  j = *valit;
	  k = asgn2index[i][j];
	  if( solution[k] < epsilon ) consistent = scope[i]->remove(j);
	} while( consistent && valit->next() ); 

// 	j = scope[i]->first();
// 	do {
// 	  jnxt = scope[i]->getNext(j);
// 	  k = asgn2index[i][j];
// 	  if( solution[k] < epsilon ) consistent = scope[i]->remove(j);
// 	  j = jnxt;
// 	} while( consistent && scope[i]->good(j) ); 
      }
    }
    
    break;
  }
    
  case SIMPLEX_PRIMAL : {
   
    consistent = smodel.primal();
    //consistent = ((ClpSimplex*)model)->primalFeasible();
    
    break;
  }

  case SIMPLEX_DUAL : {
   
    consistent = smodel.dual();
    
    break;
  }

  case SIMPLEX_BARRIER : {
   
    consistent = smodel.barrier();
    
    break;
  }

  case SIMPLEX_PRIMAL_IT : {

    double *collb  = smodel.columnLower();
    double *colub1 = smodel.columnUpper();
    double *colub2 = imodel.columnUpper();
    bool cons;

    for(i=0; consistent && i<arity; ++i) {
 //      j = scope[i]->first();
//       do {
// 	jnxt = scope[i]->getNext(j);
// 	k = asgn2index[i][j];
	
// 	collb[k] = 1.0;
// 	smodel.primal();
// 	cons = smodel.primalFeasible();

// 	if( !cons ) {
// 	  consistent = scope[i]->remove( j );
// 	  colub1[k] = 0.0;
// 	  colub2[k] = 0.0;
// 	}
// 	collb[k] = 0.0;

// 	j = jnxt;
//       } while( consistent && scope[i]->good(j) ); 

      valit = scope[i]->begin();
      do {
	j = *valit;

	k = asgn2index[i][j];	
	collb[k] = 1.0;
	smodel.primal();
	cons = smodel.primalFeasible();

	if( !cons ) {
	  consistent = scope[i]->remove( j );
	  colub1[k] = 0.0;
	  colub2[k] = 0.0;
	}
	collb[k] = 0.0;

      } while( consistent && valit->next() ); 

    }


    break;    
  }

  case SIMPLEX_DUAL_IT : {
   
    double *collb  = smodel.columnLower();
    double *colub1 = smodel.columnUpper();
    double *colub2 = imodel.columnUpper();    

    bool cons;
    for(i=0; consistent && i<arity; ++i) {
//       j = scope[i]->first();
//       do {
// 	jnxt = scope[i]->getNext(j);
// 	k = asgn2index[i][j];
	
// 	collb[k] = 1.0;
// 	smodel.dual();
// 	cons = smodel.primalFeasible();

// 	if( !cons ) {
// 	  consistent = scope[i]->remove( j );
// 	  colub1[k] = 0.0;
// 	  colub2[k] = 0.0;
// 	}
// 	collb[k] = 0.0;

// 	j = jnxt;
//       } while( consistent && scope[i]->good(j) ); 

      valit = scope[i]->begin();
      do {
	j = *valit;
	k = asgn2index[i][j];
	
	collb[k] = 1.0;
	smodel.dual();
	cons = smodel.primalFeasible();

	if( !cons ) {
	  consistent = scope[i]->remove( j );
	  colub1[k] = 0.0;
	  colub2[k] = 0.0;
	}
	collb[k] = 0.0;

      } while( consistent && valit->next() ); 
    }

    break;    
    
  }

  case SIMPLEX_BARRIER_IT : {

    double *collb = smodel.columnLower();
    double *colub1 = smodel.columnUpper();
    double *colub2 = imodel.columnUpper();
    bool cons;
    for(i=0; consistent && i<arity; ++i) {
//       j = scope[i]->first();
//       do {

// 	jnxt = scope[i]->getNext(j);
// 	k = asgn2index[i][j];
	
// 	collb[k] = 1.0;
// 	smodel.barrier();
// 	cons = smodel.primalFeasible();

// 	if( !cons ) {
// 	  consistent = scope[i]->remove( j );
// 	  colub1[k] = 0.0;
// 	  colub2[k] = 0.0;
// 	}
// 	collb[k] = 0.0;

// 	j = jnxt;
//       } while( consistent && scope[i]->good(j) ); 

      valit = scope[i]->begin();
      do {
	j = *valit;
	k = asgn2index[i][j];
	
	collb[k] = 1.0;
	smodel.barrier();
	cons = smodel.primalFeasible();

	if( !cons ) {
	  consistent = scope[i]->remove( j );
	  colub1[k] = 0.0;
	  colub2[k] = 0.0;
	}
	collb[k] = 0.0;

      } while( consistent && valit->next() ); 
    }

    break;
  }
    

  }

  return consistent;
}
  
void LpConstraint::print(std::ostream& o) const 
{
}


/**********************************************
 * MaxDiff Predicate 
 **********************************************/

PredicateMaxDiff::PredicateMaxDiff(Solver *s, VariableInt** v, const int n, const int t) 
  : LpConstraint(s, v, n)
{
  //delayed = true;
  filtering = LpConstraint::LPA; 
  initStruct(t);
}

PredicateMaxDiff::PredicateMaxDiff(Solver *s, const int ft, VariableInt** v, const int n, const int mt) 
  : LpConstraint(s, v, n)
{ 
  //delayed = true;
  filtering = ft;
  initStruct(mt);
}

void PredicateMaxDiff::initStruct(const int t) 
{ 
  int i, j, n=arity-1;

  if( filtering == LpConstraint::LPA ||
      filtering == LpConstraint::LPS ||
      filtering == LpConstraint::LPO ) {
        
    imodel.resize(n+nRows+1, 0);
    smodel.resize(n+nRows,   0);
    
    CoinBuild SbuildObject;
    CoinBuild IbuildObject;
    
    // the n first rows correspond to domain constraints
    // ( sum(x[i,j] | j in D(X[i])) == 1  )
    for (i=0;i<n;i++) {
      imodel.setRowLower(i, 1.0);
      imodel.setRowUpper(i, 1.0);
      smodel.setRowLower(i, 1.0);
      smodel.setRowUpper(i, 1.0);
    }

    // there are as many additional rows as values
    // for counting values
    nRows += n;
    for (i=n;i<nRows;i++) {
      imodel.setRowLower(i, 0.0);
      imodel.setRowUpper(i, 0.0);
      smodel.setRowLower(i, 0.0);
      smodel.setRowUpper(i, 0.0);
    }
    
    // ?    
    imodel.setRowLower(nRows, 0.0);
    imodel.setRowUpper(nRows, n*(n-1)/2);  
    
    DomainIterator *valit;
    for(i=0; i<n; ++i) {

      valit = scope[i]->begin();
      do {
	// assignment, => Xi takes j and j is used once
	int columnIndex[2] = {i, n+values2index[*valit]};
	double columnValue[2] = {1.0, 1.0};     
	// 0 -> Xi =/= j, 1 -> Xi == j
	SbuildObject.addColumn(2, columnIndex, columnValue, 0.0, 1.0);
	IbuildObject.addColumn(2, columnIndex, columnValue, 0.0, 1.0);
      } while( valit->next() ); 

//       j = scope[i]->first();
//       do {
// 	// assignment, => Xi takes j and j is used once
// 	int columnIndex[2] = {i, n+values2index[j]};
// 	double columnValue[2] = {1.0, 1.0};     
// 	// 0 -> Xi =/= j, 1 -> Xi == j
// 	SbuildObject.addColumn(2, columnIndex, columnValue, 0.0, 1.0);
// 	IbuildObject.addColumn(2, columnIndex, columnValue, 0.0, 1.0);
//       } while( scope[i]->setNext(j) ); 
    }
    
    int nValues[nRows-n];
    std::fill(nValues, nValues+nRows-n, 0);
    for(i=0; i<n; ++i) {
//       j = scope[i]->first();
//       do {
// 	// dummy vars used to express the objective function
// 	// they are symmetric, so the objective is to minimise
// 	// 0 * y1j + 1 * y2j + 2 * y3j + .... 
// 	int column1Index[1] = {n+values2index[j]};
// 	double column1Value[1] = {-1.0};
// 	SbuildObject.addColumn(1, column1Index, column1Value, 0.0, 1.0, 
// 			       (1.0 * nValues[values2index[j]]));
	
// 	int column2Index[2] = {n+values2index[j], nRows};
// 	double column2Value[2] = {-1.0, (1.0 * nValues[values2index[j]])};
// 	IbuildObject.addColumn(2, column2Index, column2Value, 0.0, 1.0);
	
// 	++nValues[values2index[j]];
//       } while( scope[i]->setNext(j) ); 

      valit = scope[i]->begin();
      do {
	// dummy vars used to express the objective function
	// they are symmetric, so the objective is to minimise
	// 0 * y1j + 1 * y2j + 2 * y3j + .... 
	j = *valit;
	int column1Index[1] = {n+values2index[j]};
	double column1Value[1] = {-1.0};
	SbuildObject.addColumn(1, column1Index, column1Value, 0.0, 1.0, 
			       (1.0 * nValues[values2index[j]]));
	
	int column2Index[2] = {n+values2index[j], nRows};
	double column2Value[2] = {-1.0, (1.0 * nValues[values2index[j]])};
	IbuildObject.addColumn(2, column2Index, column2Value, 0.0, 1.0);
	
	++nValues[values2index[j]];
      } while( valit->next() ); 

    }
    
    imodel.addColumns(IbuildObject);
    smodel.addColumns(SbuildObject);
  }

  if( filtering != LpConstraint::UPB ) {
    groundvals = new int[(ub - lb + 1)];
    groundvals -= lb;
    orderedvals =  new int[(ub - lb + 1)];
    //    toremove.init( lb, ub, BitSet::empt );
  }

  weight = 0;
  
  freevals.init( lb, ub, BitSet::empt );
}

PredicateMaxDiff::~PredicateMaxDiff()
{
  groundvals += lb;
  delete [] groundvals;
  delete [] orderedvals;
}

int PredicateMaxDiff::check( const int* s ) const 
{
  return 0;
}

int PredicateMaxDiff::approximation( int& opt ) 
{
  freevals.clear();

  int K = arity-1, i, j, k, aux, nbvals=0, nbvars=K, nbdiff=0, maxfree=0;
  std::fill( groundvals + lb, groundvals + ub + 1, 0 );
  for(i=lb; i<=ub; ++i)
    orderedvals[i] = i;

  // groundvals keeps the count of values already used by assigned variables
  // freevals is the union of unassigned variables' domains    
  for(i=0; i<K; ++i) 
    if(scope[i]->isGround()) {
      j = scope[i]->min();
      if( !groundvals[j] )
	++nbvals;
      nbdiff += (K - nbvars - groundvals[j]);
      --nbvars;
      ++groundvals[j];
    }
    else 
      scope[i]->unionTo( freevals );
  
  // separate freevals from non-freevals in groundvals (orderedval)
  // in [lb..k-1] are all pure ground values
  // in [k..ub] are other (useful) values
  i = ub;
  k = lb;
  while( i >= k ) {
    if( !freevals.member( orderedvals[i] ) ) {
      aux = orderedvals[k];
      orderedvals[k] = orderedvals[i];
      orderedvals[i] = aux;
      ++k;
    } else --i;
  }
  
  // we order values by increasing 'groundvals', that is,
  // the number of time there were used in assigned variables
  i = ub;
  while( i-- > k ) {
    j = i;
    while( j < ub && groundvals[orderedvals[j+1]] < groundvals[orderedvals[j]] ) {
      aux = orderedvals[j+1];
      orderedvals[j+1] = orderedvals[j];
      orderedvals[j] = aux;
      ++j;
    }
  }

  // now, for each unassigned var, we "assign" it to the least used
  // values, and then update both groundvals and orderedvals.
  while( nbvars ) {
    // we increment the current number of diffs by the number
    // of variables assigned to something else than orderedvals[k]
    nbdiff += (K - nbvars - groundvals[orderedvals[k]]);
    ++groundvals[orderedvals[k]];
    
    // maxfree is the maximum number of equalities induced by a choice so far
    if(maxfree < groundvals[orderedvals[k]])
      maxfree = groundvals[orderedvals[k]];
    
    // when the algorithm ends, freevals contains inconsistent values
    // therefore we remove orderedvals[k].
    freevals.erase( orderedvals[k] );

    // reorder the vals
    j = k;
    while( j < ub && 
	   groundvals[orderedvals[j+1]] < 
	   groundvals[orderedvals[j]] ) {
      aux = orderedvals[j+1];
      orderedvals[j+1] = orderedvals[j];
      orderedvals[j] = aux;
      ++j;
    }
    --nbvars;
  }
  
  if(opt < nbdiff) opt = nbdiff;


  // now we remove consistent values from freevals
  if( k <= ub && !(freevals.empty()) ) {

    // for each value j in freevals;
    j = freevals.min();
    do {
      // if we swap the worse choice so far with j, we obtain:
      // (opt + maxfree - groundvals[j])
      if( (opt + maxfree - groundvals[j]) > scope[K]->min() ) {
	// if this is less differences than allowed, j can be removed
	freevals.erase( j );
      }
      j = freevals.next(j);
    } while( j != NOVAL );
  }

  return nbdiff;
}

int PredicateMaxDiff::upperbound()
{
  int K = arity-1;
  freevals.clear();
  for(int i=0; i<K; ++i) {
    scope[i]->unionTo( freevals );
  }
  int c = freevals.size();
  
  int res = (K*(K-1)/2);
  if( c < K ) 
    res -= ((c + (K%c))*((int)(K/c)) - c);

  return res;
}


bool PredicateMaxDiff::propagate(const int changedIdx, const int e)
{
  int K = arity-1;
  bool consistent = true;
  DomainIterator *valit;

  /// only filter if a strictly positive number of differences is required
  if( scope[K]->min() ) {

    int i, j, nbdiff, jnxt, k = (((K) * (K-1)) >> 1);

    switch( filtering ) {

      // simple upper bound
    case LpConstraint::UPB : {
      nbdiff = upperbound();
      if( !scope[K]->setMax( nbdiff) ) consistent = false;
      break;
    }

      // approximation only
    case LpConstraint::APX : {
      nbdiff = scope[K]->max();
      nbdiff = approximation( nbdiff );

      if( !scope[K]->setMax( nbdiff ) ) consistent = false;


      if( consistent && !freevals.empty() ) {
	
	int isngrd[K];
	for(i=0; i<K; ++i)
	  isngrd[i] = !(scope[i]->isGround()); 
	
	for(i=0; consistent && i<K; ++i) 
	  if( isngrd[i] ) {
	    j = freevals.min();
	    do {
	      consistent = scope[i]->remove( j );
	      j = freevals.next(j);
	    } while( consistent && j != NOVAL );
	  }
      }
      
      break;
    }

    case LpConstraint::LPO : {
      
      updateModel();
      
      smodel.dual();
      nbdiff = (k - (int)(ceil(smodel.objectiveValue())));
       
      if( !scope[K]->setMax( nbdiff ) ) consistent = false;
      
      if( consistent ) {

	imodel.setRowLower(nRows, (double)( k - scope[K]->max() ));
	imodel.setRowUpper(nRows, (double)( k - scope[K]->min() ));
	
	imodel.primalDual();  
	
	if( (consistent = imodel.primalFeasible()) ) {
	  
	  const double * solution = imodel.primalColumnSolution();
	  for(i=0; consistent && i<K; ++i) {
// 	    j = scope[i]->first();
// 	    do {
// 	      jnxt = scope[i]->getNext(j);
// 	      k = asgn2index[i][j];	    
// 	      if( solution[k] < epsilon ) 
// 		consistent = scope[i]->remove(j);
// 	      j = jnxt;
// 	    } while( consistent && scope[i]->good(j) );

	    valit = scope[i]->begin();
	    do {
	      j = *valit;
	      k = asgn2index[i][j];	    
	      if( solution[k] < epsilon ) 
		consistent = scope[i]->remove(j);
	    } while( consistent && valit->next() );
	    
	  }
	}
      }
      
      break;
    }

    case LpConstraint::LPA : {

      nbdiff = scope[K]->max();
      nbdiff = approximation( nbdiff );

      if( !scope[K]->setMax( nbdiff ) ) consistent = false;
      
      if( consistent && !freevals.empty() ) {
	
	int isngrd[K];
	for(i=0; i<K; ++i)
	  isngrd[i] = !(scope[i]->isGround()); 
	
	for(i=0; consistent && i<K; ++i) 
	  if( isngrd[i] ) {

	    j = freevals.min();
	    do {
	      consistent = scope[i]->remove( j );
	      j = freevals.next(j);
	    } while( consistent && j != NOVAL );
	  }
      }

      if( consistent ) {
	
	updateModel();
	
	smodel.dual();
	nbdiff = (k - (int)(ceil(smodel.objectiveValue())));
	
	if( !scope[K]->setMax( nbdiff ) ) {
	  //return K+1;
	  //return scope[K];
	  //return NULL;
	  return false;
	}
	
	imodel.setRowLower(nRows, (double)( k - scope[K]->max() ));
	imodel.setRowUpper(nRows, (double)( k - scope[K]->min() ));
	
	imodel.primalDual();  
	
	if( (consistent = imodel.primalFeasible()) ) {
	  
	  const double * solution = imodel.primalColumnSolution();
	  for(i=0; consistent && i<K; ++i) {

// 	    j = scope[i]->first();
// 	    do {
// 	      jnxt = scope[i]->getNext(j);
// 	      k = asgn2index[i][j];	    
// 	      if( solution[k] < epsilon ) 
// 		consistent = scope[i]->remove(j);
// 	      j = jnxt;
// 	    } while( consistent && scope[i]->good(j) ); 

	    valit = scope[i]->begin();
	    do {
	      j = *valit;
	      k = asgn2index[i][j];	    
	      if( solution[k] < epsilon ) 
		consistent = scope[i]->remove(j);
	    } while( consistent && valit->next() ); 
	  }
	}
      }
      
      break;
    }

    case LpConstraint::LPS : {

      updateModel();
      
      smodel.dual();
      nbdiff = (k - (int)(ceil(smodel.objectiveValue())));
       
      if( !scope[K]->setMax( nbdiff ) ) consistent = false;

      if( consistent ) {
	
	nbdiff = scope[K]->max();
	nbdiff = approximation( nbdiff );
	
	if( !freevals.empty() ) {

	  if( !scope[K]->setMax( nbdiff ) ) consistent = false;
	  
	  if( consistent ) {
	    int isngrd[K];
	    for(i=0; i<K; ++i)
	      isngrd[i] = !(scope[i]->isGround()); 
	    
	    for(i=0; consistent && i<K; ++i) 
	      if( isngrd[i] ) {
		j = freevals.min();
		do {
		  consistent = scope[i]->remove( j );
		  j = freevals.next(j);
		} while( consistent && j != NOVAL );
	      }
	  }
	}
      }
	
      break;
    }

    }
  }
    
  return consistent;
}
  
void PredicateMaxDiff::print(std::ostream& o) const 
{
  int K=arity-1;
  o << "#(";
  for(int i=0; i<K-1; ++i) {
    scope[i]->printshort(o);
    o << ", ";
  }
  scope[K-1]->printshort(o);
  o << ") = ";
  scope[K]->printshort( o );
}


/**********************************************
 * MinDiff Predicate 
 **********************************************/

PredicateMinDiff::PredicateMinDiff(Solver *s, VariableInt** v, const int n, const int t) 
  : LpConstraint(s, v, n)
{
  filtering = LpConstraint::APX; 
  initStruct(t);
}

PredicateMinDiff::PredicateMinDiff(Solver *s, const int f, VariableInt** v, const int n, const int t) 
  : LpConstraint(s, v, n)
{ 
  filtering = f;
  initStruct(t);
}

void PredicateMinDiff::initStruct(const int t) 
{ 

  int i, j, k, nRows=arity-1;

  K = arity-1;
  vars = new VariableInt*[K];
  occurrence = new int[ub-lb+1];
  occurrence -= lb;
  orderedvals = new int[ub-lb+1];

  std::fill( occurrence+lb, occurrence+ub+1, 0 );
  
  weight = 0;


  if( filtering == LpConstraint::LPA ||
      filtering == LpConstraint::LPS ||
      filtering == LpConstraint::LPO ) 
    {
      
      imodel.resize(nRows, 0);
      smodel.resize(nRows, 0);
      
      CoinBuild IbuildObject;
      CoinBuild SbuildObject;
      
      for (i=0;i<nRows;i++) {
	imodel.setRowLower(i, 1.0);
	imodel.setRowUpper(i, 1.0);
	smodel.setRowLower(i, 1.0);
	smodel.setRowUpper(i, 1.0);
      }
      
      nColumns = 0;

      DomainIterator *valit;
      for(i=0; i<nRows; ++i) {
// 	j = scope[i]->first();
// 	do {
// 	  ++occurrence[j];
// 	  ++nColumns;
// 	  int columnIndex[1] = {i};
// 	  double columnValue[1] = {1.0};     
// 	  IbuildObject.addColumn(1, columnIndex, columnValue, 0.0, 1.0);
// 	  SbuildObject.addColumn(1, columnIndex, columnValue, 0.0, 1.0);
// 	} while( scope[i]->setNext(j) ); 

	valit = scope[i]->begin();
	do {
	  j = *valit;
	  ++occurrence[j];
	  ++nColumns;
	  int columnIndex[1] = {i};
	  double columnValue[1] = {1.0};     
	  IbuildObject.addColumn(1, columnIndex, columnValue, 0.0, 1.0);
	  SbuildObject.addColumn(1, columnIndex, columnValue, 0.0, 1.0);
	} while( valit->next() ); 
      } 
      imodel.addColumns(IbuildObject);
      smodel.addColumns(SbuildObject);

      int nProd = 0;
      for(i=lb; i<=ub; ++i)
	nProd += ((occurrence[i] - 1) * (occurrence[i]) / 2);
      
      int *start = new int[nColumns+1];
      int *column = new int[nProd];
      double *element = new double[nProd];
      
      int idx = 0, str = 0;
      for(j=0; j<nRows; ++j) {
	for(i=0; i<nvals; ++i) {
	  if(scope[j]->contain(index2values[i])) {
	    start[idx] = str;
	    for(k=j+1; k<nRows; ++k)
	      if(scope[k]->contain(index2values[i])) {
		

		column[str] = asgn2index[k][index2values[i]];
		
		assert( index2varasgn[column[str]] == k );
		assert( index2varasgn[column[str]] == k );
		assert( index2valasgn[column[str]] == index2values[i] );
		
		
		element[str] = -1.0;
		++str;
	      }
	    ++idx;
	  }
	}
      }
      start[idx] = str;

      imodel.loadQuadraticObjective(nColumns,start,column,element);
      smodel.loadQuadraticObjective(nColumns,start,column,element);

 //      ClpCholeskyBase * 
      cholesky = new ClpCholeskyBase();
      cholesky->setKKT(true);
      imodel.setCholesky(cholesky);

      delete [] start;
      delete [] column;
      delete [] element;

    }
}

PredicateMinDiff::~PredicateMinDiff()
{
  occurrence += lb;
  delete [] occurrence;
  delete [] vars;
  delete [] orderedvals;
}

int PredicateMinDiff::approximation( int& opt ) 
{
  int nbdiff = (((K)*(K-1)) >> 1);
  int i, j, nvals;//, vali;

  int nvars=K, aux;
  VariableInt *auxv;
  int *firstval = orderedvals;

  /// count values
  std::fill( occurrence+lb, occurrence+ub+1, 0 );
  nvals=0;
  DomainIterator *valit;
  for(i=0; i<K; ++i) {
    vars[i] = scope[i];

    valit = scope[i]->begin();
    do if( ++occurrence[ *valit ] == 1 )
	orderedvals[nvals++] = *valit;
    while( valit->next() );

//     vali = scope[i]->first();
//     do {
//       if( ++occurrence[ vali ] == 1 )
// 	orderedvals[nvals++] = vali;
//     }
//     while( scope[i]->setNext(vali) );
  }
     
  int *lastval = orderedvals+nvals;

  while( nvars ) {
    
    /// order values
    for(i=1; firstval+i<lastval; ++i) {
      j = i;
      while( j && occurrence[firstval[j]] > occurrence[firstval[j-1]] )
	{
	  vali = firstval[j];
	  firstval[j] = firstval[j-1];
	  firstval[j-1] = vali;
	  --j;
	}
    }
    
    /// use the value with largest occurrences
    aux = occurrence[*firstval];
    aux = ((aux * (aux-1)) >> 1);
    nbdiff -= aux;

    // update the data-structures
    i = nvars;
    while( i-- ) 
      if( vars[i]->contain( *firstval ) )
	{
	  // update the array 'occurrence'
// 	  vali = vars[i]->first();
// 	  do --occurrence[ vali ];
// 	  while( vars[i]->setNext(vali) );

	  valit = vars[i]->begin();
	  do --occurrence[ *valit ];
	  while( valit->next() );
	  
	  // put all used variables at the end of the array 'vars'
	  --nvars;
	  auxv = vars[nvars];
	  vars[nvars] = vars[i];
	  vars[i] = auxv;
	}
    //--nvals;
    ++firstval;
  }

  return nbdiff;
}

bool PredicateMinDiff::propagate(const int changedIdx, const int e)
{
  bool consistent = true;
  int i, j, nbdiff, jnxt, k = (((K) * (K-1)) >> 1);
  
  if( scope[K]->max() < k ) {

    switch( filtering ) {

    case LpConstraint::APX : {
      nbdiff = scope[K]->min();
      nbdiff = approximation( nbdiff );

      if( !scope[K]->setMin( nbdiff ) ) consistent = false;
     
      break;
    }

    case LpConstraint::LPO : {

      updateModel();

//       ClpInterior barrier;
//       barrier.borrowModel(imodel);
//       ClpCholeskyBase * cholesky = new ClpCholeskyBase();
//       cholesky->setKKT(true);
//       barrier.setCholesky(cholesky);
//       barrier.primalDual();
//       int barrierobj = ((K * (K-1) / 2) + (int)(barrier.objectiveValue()));
//       barrier.returnModel(imodel);


      smodel.primal();
      int simplexobj = (k + (int)round(smodel.objectiveValue()));
      
      if( !(smodel.primalFeasible()) ) simplexobj = -1;
      //ClpCholeskyBase * cholesky = new ClpCholeskyBase();
      //cholesky->setKKT(true);

      imodel.primalDual(); 

      nbdiff = scope[K]->min();
      nbdiff = approximation( nbdiff );

      int objective = (k + (int)round(imodel.objectiveValue()));

      if( !(imodel.primalFeasible()) ) objective = -1;

      if( simplexobj < 0 ) std::cout << " - ";
      else std::cout << std::setw(2) << simplexobj << " ";

      if( objective < 0 ) std::cout << " - ";
      else std::cout << std::setw(2) << objective << " ";

      std::cout << std::setw(2) << nbdiff << std::endl; 



      if( objective != nbdiff ) {//|| objective != simplexobj || simplexobj != nbdiff ) {

	std::cout << "I-model" << std::endl;

	const double * solution = imodel.primalColumnSolution();
	for(j=0; j<nColumns; ++j) 
	  if( solution[j] )
	    std::cout << std::setw(12) << solution[j] << ".x" << index2varasgn[j] << " = "
		 << index2valasgn[j] << std::endl;
	std::cout << "\t" << (imodel.objectiveValue()) << std::endl;

	std::cout << smodel.objectiveValue() << " " 
	     << (floor(smodel.objectiveValue())) << " " 
	     << (ceil(smodel.objectiveValue())) << " " 
	     << ((int)(smodel.objectiveValue())) << " " 
	     << std::endl
	     << imodel.objectiveValue() << " " 
	     << (floor(imodel.objectiveValue())) << " " 
	     << (ceil(imodel.objectiveValue())) << " " 
	     << ((int)(imodel.objectiveValue())) << " " 
	     << std::endl;

	for(i=0; i<K; ++i) {
	  scope[i]->print( std::cout );
	  std::cout << std::endl;
	}
	std::cout << std::endl;
	
	imodel.writeMps( "imodel.mps" );
	smodel.writeMps( "smodel.mps" );

	ifstream imod("imodel.mps");
	ifstream smod("smodel.mps");
	char ci, cs;
	while( imod.good() && smod.good() )
	  {
	    imod.get( ci );
	    smod.get( cs );
	    if( ci != cs ) {
	      std::cout << "DIFF" << std::endl;
	      exit(0);
	    }
	  }
	if( imod.good() || smod.good() ) {
	  std::cout << "DIFF" << std::endl;
	}
	//exit(0);
      }


     if( simplexobj != nbdiff ) {//|| objective != simplexobj || simplexobj != nbdiff ) {

	std::cout << "S-model" << std::endl;

	for(i=0; i<K; ++i) {
	  scope[i]->print( std::cout );
	  std::cout << std::endl;
	}
	std::cout << std::endl;
	

	const double * solution = smodel.primalColumnSolution();
	for(j=0; j<nColumns; ++j) 
	  if( solution[j] )
	    std::cout << std::setw(12) << solution[j] << ".x" << index2varasgn[j] << " = "
		 << index2valasgn[j] << std::endl;
	std::cout << "\t" << (smodel.objectiveValue()) << std::endl;


	smodel.writeMps( "smodel.mps" );
	
	ClpSimplex secondmodel;
	secondmodel.setLogLevel( 0 );
	
	secondmodel.readMps("smodel.mps");

	secondmodel.primal();

	const double * solution2 = secondmodel.primalColumnSolution();
	for(j=0; j<nColumns; ++j) 
	  if( solution2[j] )
	    std::cout << std::setw(12) << solution2[j] << ".x" << index2varasgn[j] << " = "
		 << index2valasgn[j] << std::endl;
	std::cout << "\t" << (secondmodel.objectiveValue()) << std::endl;

	exit(0);
      }
      
      consistent = scope[K]->setMin( nbdiff );
     
      break;
    }
    }
  }
    
  return consistent;
}
  
void PredicateMinDiff::print(std::ostream& o) const 
{
  o << "#(";
  for(int i=0; i<K-1; ++i) {
    scope[i]->printshort(o);
    o << ", ";
  }
  scope[K-1]->printshort(o);
  o << ") = ";
  scope[K]->printshort( o );
}

#endif



ConstraintSatLess::ConstraintSatLess(Solver* s, VariableInt** v, const int dur, 
				     ConstraintClauseBase *st)
  : Constraint(s, v, 2, RANGETRIGGER),
    min_0(((VariableRange*)(v[0]))->vmin.value),
    min_1(((VariableRange*)(v[1]))->vmin.value),
    max_0(((VariableRange*)(v[0]))->vmax.value),
    max_1(((VariableRange*)(v[1]))->vmax.value)
{
  sat = st;
  duration = dur;
  sat->addPrecedence( this );
  responsible_lb = &(sat->failure[0]);
  responsible_ub = &(sat->failure[1]);
  responsible_disjunct = &(sat->failure[2]);
}
ConstraintSatLess::~ConstraintSatLess()
{
}
int ConstraintSatLess::check( const int* s ) const 
{
  return ( s[0] + duration > s[1] );
}
bool ConstraintSatLess::propagate(const int changedIdx, const int e)
{
  bool consistent = false;
  if( min_0 + duration <= max_1)
    {
      if( changedIdx ) {
	if(  max_0 > max_1 - duration ) {
	  scope[0]->setMax(max_0 - duration);
	  *reason_for_ub_decrease = atom_id;
	}
      } else {
	if(  min_1 < min_0 + duration ) {
	  scope[1]->setMin(min_0 + duration);
	  *reason_for_lb_increase = atom_id;
	}
      }
      consistent = true;
    }
  else 
    {
      *responsible_lb = task_id[0];
      *responsible_ub = task_id[1];
      *responsible_disjunct = atom_id;
    }
  return consistent;
}
void ConstraintSatLess::print(std::ostream& o) const 
{
  scope[0]->printshort( o ) ;
  if( !duration )
    o << " <= ";
  else if( duration == 1 )
    o << " < ";
  else 
    o << " + " << duration-1 << " < " ;
  scope[1]->printshort( o ) ;
}

/**********************************************
 * Disjunct Predicate 
 **********************************************/
/// 0 :-> (x0 + d0 < x1), 1 :-> (x1 + d1 < x0) 
PredicateDisjunct::PredicateDisjunct(Solver *s, VariableInt** v, 
				     const int* dur, ConstraintClauseBase *st) 
  : Constraint(s, v, 3, RANGETRIGGER), 
    state(((VariableBool*)(v[2]))->domain.state),
    min_0(((VariableRange*)(v[0]))->vmin.value),
    min_1(((VariableRange*)(v[1]))->vmin.value),
    max_0(((VariableRange*)(v[0]))->vmax.value),
    max_1(((VariableRange*)(v[1]))->vmax.value)
    //failure(st->failure)
{
  sat = st;
  duration[0] = dur[0];
  duration[1] = dur[1];
  sat->addDisjunct( this );
  //atom_id = scope[2]->id+1;
  //failed_task[0] = sat->failed_task;
  //failed_task[1] = sat->failed_task+1;
  responsible_lb = &(sat->failure[0]);
  responsible_ub = &(sat->failure[1]);
  responsible_disjunct = &(sat->failure[2]);

  //std::cout << atom_id << " " << (&(sat->reason[atom_id])) << std::endl;

  //reason = &(sat->reason[atom_id]);
}

PredicateDisjunct::~PredicateDisjunct()
{
}

int PredicateDisjunct::check( const int* s ) const 
{
  return ( (s[2] && (s[0] + duration[0] <= s[1])) 
	   ||
	   (!s[2] && (s[1] + duration[1] <= s[0])) );
}

bool PredicateDisjunct::propagate(const int changedIdx, const int e) 
{
  int consistent=true;
  if( changedIdx != 2 ) {
    if( min_0 + duration[0] > max_1 ) {  
      if( state == 3 ) {
  
#ifdef _DEBUGNOGOOD
	std::cout << " prune disjunct b" << atom_id << " = 1  because ";
#endif
	
	scope[2]->setDomain(1);
	
	//sat->addLiteral(atom_id, atom_id);
	*reason = NULL;
	*reason_for_choice[0] = *reason_for_lb_increase[0];
	*reason_for_choice[1] = *reason_for_ub_decrease[1];
	
#ifdef _DEBUGNOGOOD
	std::cout << (*reason_for_choice[0]) << " and " << (*reason_for_choice[1]) << std::endl;
#endif

      } else if( state == 1 ) {

#ifdef _DEBUGNOGOOD
	std::cout << " should prune disjunct b" << atom_id << " = 1" << std::endl;
#endif

	consistent = false;
      }

    } else if( min_1 + duration[1] > max_0 ) {
      if( state == 3 ) {

#ifdef _DEBUGNOGOOD
	std::cout << " prune disjunct b" << atom_id << " = 0  because ";
#endif

	scope[2]->setDomain(0);

	//sat->addLiteral(atom_id, -atom_id);
	*reason = NULL;
// 	*reason_for_choice[0] = *reason_for_ub_decrease[0];
// 	*reason_for_choice[1] = *reason_for_lb_increase[1];
	*reason_for_choice[0] = *reason_for_lb_increase[1];
	*reason_for_choice[1] = *reason_for_ub_decrease[0];

#ifdef _DEBUGNOGOOD
	std::cout << (*reason_for_choice[0]) << " and " << (*reason_for_choice[1]) << std::endl;
#endif

      } else if( state == 2 ) {

#ifdef _DEBUGNOGOOD
	std::cout << " should prune disjunct b" << atom_id << " = 0" << std::endl;
#endif

	consistent = false;
      }

    }
  }

  int bound;
  if( state == 1 ) {

    // x + a <= y
    bound = min_0 + duration[0];
    if( bound > max_1 ) {

#ifdef _DEBUGNOGOOD
	std::cout << " should prune lb(t" << task_id[1] << ")" << std::endl;
#endif

      consistent = false;
    } else if( bound > min_1 ) {

      scope[1]->setMin( bound );
      *reason_for_lb_increase[1] = atom_id;

#ifdef _DEBUGNOGOOD
      std::cout << " prune lb(t" << task_id[1] << ") because b" << atom_id << " & lb(t" << task_id[0] << ")" << std::endl;
#endif

    }

    bound = max_1 - duration[0];
    if( bound < min_0 ) {

#ifdef _DEBUGNOGOOD
	std::cout << " should prune ub(t" << task_id[0] << ")" << std::endl;
#endif

      consistent = false;
    } else if( bound < max_0 ) {

      scope[0]->setMax( bound );
      *reason_for_ub_decrease[0] = atom_id;

#ifdef _DEBUGNOGOOD
      std::cout << " prune ub(t" << task_id[0] << ") because b" << atom_id << " & ub(t" << task_id[1] << ")" << std::endl;
#endif

    }

  } else if( state == 2 ) {

    // y + b <= x
    bound = min_1 + duration[1];
    if( bound > max_0 ) {

#ifdef _DEBUGNOGOOD
	std::cout << " should prune lb(t" << task_id[0] << ")" << std::endl;
#endif

      consistent = false;
    } else if( bound > min_0 ) {

      scope[0]->setMin( bound );
      *reason_for_lb_increase[0] = atom_id;

#ifdef _DEBUGNOGOOD
      std::cout << " prune lb(t" << task_id[0] << ") because b" << atom_id << " & lb(t" << task_id[1] << ")" << std::endl;
#endif

    }

    bound = max_0 - duration[1];
    if( bound < min_1 ) {

#ifdef _DEBUGNOGOOD
	std::cout << " should prune ub(t" << task_id[1] << ")" << std::endl;
#endif

      consistent = false;
    } else if( bound < max_1 ) {

      scope[1]->setMax( bound );
      *reason_for_ub_decrease[1] = atom_id;

#ifdef _DEBUGNOGOOD
      std::cout << " prune ub(t" << task_id[1] << ") because b" << atom_id << " & ub(t" << task_id[0] << ")" << std::endl;
#endif

    }

  } 


  if(!consistent) {

#ifdef _DEBUGNOGOOD  
    std::cout << "FAILURE on b" << atom_id << " = t" << task_id[0] << " . t" << task_id[1] << "   (";
    print( std::cout );
    std::cout << " -> " ;
    scope[0]->print( std::cout );
    std::cout << " / ";
    scope[1]->print( std::cout );
    std::cout << ")" << std::endl;
#endif

    *responsible_disjunct = atom_id;
    if( state == 1 ) { // x[0] < x[1]
	
	*responsible_lb = task_id[0];
	*responsible_ub = task_id[1];

    } else { // x[0] < x[1]
       
	*responsible_lb = task_id[1];
	*responsible_ub = task_id[0];

    }

    sat->update();
  }

  return consistent;
}

void PredicateDisjunct::print(std::ostream& o) const
{  
  
  //if( !scope[2]->isGround() ) {
     scope[2]->print( o );
     o << ": ";
     //}
  if( scope[2]->contain(0) ) {
    o << "(";
    scope[0]->print( o );
    if( !duration[0] )
      o << " <= ";
    else if( duration[0] == 1 )
      o << " < ";
    else 
      o << " + " << duration[0]-1 << " < " ;
    scope[1]->print( o );
    o << ")";
  }
  if( scope[2]->domsize() == 2 )
    o << " || ";
  if( scope[2]->contain(1) ) {
    o << "(";
    scope[1]->print( o );
    if( !duration[1] )
      o << " <= ";
    else if( duration[1] == 1 )
      o << " < ";
    else 
      o << " + " << duration[1]-1 << " < " ;
    scope[0]->print( o );
    o << ")";
  }
  //  o << " [" << vmin << ".." << vmax << "]";
 
}


/**********************************************
 * Clause Base
 **********************************************/ 

ConstraintClauseBase::ConstraintClauseBase(Solver* s)
  : Constraint(s), SatSolver()
{
  scope = NULL;
  solver = s;
  conflict = NULL;
  use_learning = true;
  use_backjump = true;
  is_choice = false;
  weight = 0;
}

ConstraintClauseBase::~ConstraintClauseBase()
{
  for(int i=0; i<cl_scope.size; ++i) {
    delete [] (cl_scope.stack_[i]);
    delete [] (cl_lits.stack_[i]);
  }
}

inline int ConstraintClauseBase::check( const int* s ) const 
{
  // TODO
  return 0;
}

inline bool ConstraintClauseBase::propagate()
{
  unitPropagate();
  int consistent = true;
  Atom a;
  for(unsigned int i=0; consistent && i<assumptions.size; ++i) {
    a = assumptions[i];
    consistent = X[a]->setDomain( a==polarity[a] );
  }
  return consistent;
}

void ConstraintClauseBase::checkConsistency()
{
  BitSet ground_vars(-arity, arity, BitSet::empt);
  for(int i=0; i<solver->length; ++i)
    if( solver->sequence[i]->isGround() ) {
      int idx = solver->sequence[i]->id+1;
      int val = solver->sequence[i]->value();
      ground_vars.insert( val ? idx : -idx );
    }
  
  BitSet literals(-arity, arity, BitSet::empt);
  for(unsigned int i=0; i<assumptions.size; ++i)
    literals.insert( polarity[assumptions[i]] );

  bool equal = ( literals.included(literals) && literals.included(ground_vars) );
  if( !equal ) {
    literals.print( std::cout );
    std::cout << std::endl;
    ground_vars.print( std::cout );
    std::cout << std::endl;
    printDecisions( std::cout );
    exit(0);
  }
}

inline bool ConstraintClauseBase::propagate(const int changedIdx, const int e)
{
  bool consistent = true;
  Atom a = changedIdx+1;
  Literal p = ( scope[changedIdx]->equal(1) ? a : -a );
  conflict = NULL;
  unsigned int i;
  int cw;

  //lvl[a] = *level;

  if( is_choice ) {

#ifdef _DEBUGNOGOOD
    std::cout << "c  notify choice of variable ";
    scope[changedIdx]->print( std::cout );
    std::cout << " -> " << p << std::endl;    
#endif

    if(mode == DTP) {
      reason_lit[0][a] = 0;
      reason_lit[1][a] = 0;
    }

    //std::cout << "decision at level " << (solver->level) << std::endl;

    makeDecision( p );
    
    store->push( this );
    is_choice = false;

  } else {

#ifdef _DEBUGNOGOOD
    std::cout << "c propagate after event on " ;
    scope[changedIdx]->print( std::cout );
    std::cout << " -> " << p << std::endl;
#endif    

    // channel from cp to sat
    if( !assumptions.member(a) ) {
      addLiteral(a, p);
    }
  }
  
  // unit propag
  i = assumptions.size;
  conflict = NULL;  
  cw = isWatchedBy[-p].size;
  while(cw-- && !conflict) {
    conflict = updateWatcher(cw, a, -p);
  }
  consistent = !conflict;
  
  if( consistent ) {
    // channel from sat to cp
    if( i < assumptions.size )       
      for(; consistent && i<assumptions.size; ++i) {
	a = assumptions[i];
	consistent &= X[a]->setDomain( a==polarity[a] );

	if( !consistent ) {

#ifdef _DEBUGNOGOOD      
 	  std::cout << "FAIL ON ";
	  X[a]->print( std::cout );
	  std::cout << " ";
	  printClause( std::cout, reason[a] );
	  std::cout << std::endl;
#endif

	  conflict = reason[a];
	  //reason[a] = NULL;
	  update();
	}

      }
  } else {

//     //std::cout << "sat-inconsistent " ;

//     if( i < assumptions.size )       
//       for(; consistent && i<assumptions.size; ++i) {
	
// 	a = assumptions[i];

// 	std::cout << a << " ";

// 	X[a]->setDomain( a==polarity[a] );
//       }
//     std::cout << std::endl;
    update();
  } 
  
  return consistent;
}


void ConstraintClauseBase::analyzeConflict() 
{

#ifdef _DEBUGNOGOOD
   std::cout << "analyze " << (*level) << " ";
   if(mode == SAT) std::cout << " SAT ";
   if(mode == CSP) std::cout << " CSP ";
   if(mode == DTP) std::cout << " DTP ";

   if( conflict )
     printClause( std::cout, conflict );
   std::cout << std::endl;
#endif

  Atom a;
  Literal p;
  int btLevel;
  

  //  std::cout << (*level) << std::endl;

  if(*level > solver->init_level) {
    switch(mode) {
      //case CSP: {
      //} break;
    case SAT: {
      
      //std::cout << *level << " " << solver->init_level << std::endl;
      
      //std::cout << conflict << std::endl;

      btLevel = *level-1;//analyze( conflict, p, use_learning );
      /*
      if( use_backjump || use_learning ) {
	a = atom(p);
	solver->decision[*level] = X[a];
	X[a]->postCut( p!=a );
	if( use_backjump )
	  solver->backtrackLevel = btLevel+1;
      }
      */   
    } break;
    case DTP: {

#ifdef _DEBUGNOGOOD      
      std::cout << "c  analyze conflict (Disjunctive Temporal Problem)" << std::endl;
      std::cout << "c failed on ";     
      if( conflict )
	{
	  std::cout << " cl";
	  printClause( std::cout, conflict);
	  std::cout << std::endl;
	}
      else
	{
	  std::cout << "b" << failure[2] << " = t" << failure[0] << " < t" << failure[1] << std::endl; 
	}
      int i;
      for(i=0; i<nTasks; ++i) {

	std::cout //<< (the_task[i]->id+1) << ": " 
	  << std::setw(2) << i << ": " 
	  << std::setw(2) << reason_lb[i] << " " 
	  << std::setw(2) << reason_ub[i] << "  -- " ;
	the_task[i]->print( std::cout );
	std::cout << std::endl;
      }

      int n = (solver->empty-solver->sequence);
      for(i=0; i<n; ++i)
	if(solver->sequence[i]->isGround())
	  {
	    if( solver->sequence[i]->id < numAtoms ) {
	      a = solver->sequence[i]->id+1;
	      p = (solver->sequence[i]->equal(1) ? a : -a);
	      std::cout << lvl[a] << std::setw(5) << p << "  (t" << std::setw(2) << disjunct[0][a] 
		 << ", t" << std::setw(2) << disjunct[1][a] << ")   ";
	    if( lvl[a] == 0 ) {
	      std::cout << "entailed";
	    } else if( reason[a] ) {    
	      printClause( std::cout, reason[a] );
	    } else if( (reason_lit[0][a] + reason_lit[1][a]) ) {
	      int lb, d_lb, ub, d_ub;
	      // therefore the lb of X0 and the ub of X1 are responsible
	      d_lb = reason_lit[0][a];
	      d_ub = reason_lit[1][a];
	      if( d_lb ) lb = disjunct[X[d_lb]->equal(1)][d_lb];
	      if( d_ub ) ub = disjunct[X[d_ub]->equal(0)][d_ub];
	      
	      if( d_lb )
		std::cout << ( X[d_lb]->equal(1) ? " b" : "-b" ) << d_lb 
		     << " lb(t" << lb << ") "; 
	      else 
		std::cout << " - ";
	      if( d_ub )
		std::cout << ( X[d_ub]->equal(1) ? " b" : "-b" ) << d_ub 
		     << " ub(t" << ub << ") ";
	      else 
		std::cout << " - ";	    
	    } else {
	      std::cout << "decision";
	    }
	    std::cout << std::endl;
	  }
	}
#endif      
      
      //exploreImplicationGraph();
      extractNogood();
      
      p = 0;
      btLevel = -1;
      int sLevel = -1;
      int cLevel = *level;

      if( learnt_clause.size > 1 ) {

	int j, idx_top=-1, idx_btk=-1;

#ifdef _DEBUGNOGOOD      
	int numSameLevel = 0;

	std::cout << cLevel << std::endl << "(";
	for(j=0; j<learnt_clause.size; ++j) {
	  a = atom(learnt_clause[j]);
	  std::cout << " " << std::setw(4) << learnt_clause[j] ;
	}
	std::cout << " ) " << std::endl;
	
	std::cout << "(";
	for(j=0; j<learnt_clause.size; ++j) {
	  a = atom(learnt_clause[j]);
	  std::cout << " " << std::setw(4) << lvl[a] ;
	}
	std::cout << " ) " << std::endl;
#endif
	


	Literal q;
	for(j=0; j<learnt_clause.size; ++j) {
	  q = learnt_clause[j];
	  a = atom(q);
	  activity[q] += params.activity_increment;
	  if( lvl[a] > sLevel ) {
	    btLevel = sLevel;
	    idx_btk = idx_top;
	    sLevel = lvl[a];
	    idx_top = j;
	    p = learnt_clause[j];
	  } else if( lvl[a] > btLevel ) {
	    btLevel = lvl[a];
	    idx_btk = j;
	  }
	}
	
	//std::cout << idx_top << " " << idx_btk << std::endl;
	if( !idx_btk ) idx_btk = idx_top;
	
	if(idx_top) {
	  // swap top-level literal with first literal
	  j = learnt_clause[0];
	  learnt_clause[0] = learnt_clause[idx_top];
	  learnt_clause[idx_top] = j;
	}

	if(idx_btk != 1) {
	  // swap backtrack-level literal with second literal
	  j = learnt_clause[idx_btk];
	  learnt_clause[idx_btk] = learnt_clause[1];
	  learnt_clause[1] = j;
	}

	//assert(lvl[atom(learnt_clause[0])] >= cLevel);



#ifdef _DEBUGNOGOOD      

 	std::cout << " conflict " << (learnt_clause.size);
 	for(j=0; j<learnt_clause.size; ++j) 
 	  std::cout << " " << learnt_clause[j] ;
 	std::cout << std::endl;


	std::cout << "(";
	for(j=0; j<learnt_clause.size; ++j) {
	  a = atom(learnt_clause[j]);
	  if( cLevel <= lvl[a] ) 
	    ++numSameLevel; 
	  std::cout << " " << learnt_clause[j] ;
	}
	std::cout << " ) ==> " << (cLevel - btLevel) << std::endl << std::endl;
	
	assert( lvl[atom(learnt_clause[0])] >= cLevel );
	assert( numSameLevel == 1 );
#endif

// 	if( learnt_clause.size == 6 &&
// 	    learnt_clause[0] == 24 && 
// 	    learnt_clause[1] == 18 && 
// 	    learnt_clause[2] == 16 && 
// 	    learnt_clause[3] == 19 && 
// 	    learnt_clause[4] == 23 && 
// 	    learnt_clause[5] == 12 )
// 	  exit( 0 );
	
	
	a = atom(p);
	addClause( learnt, learnt_clause, stats.learnt_avg_size );
	reason[a] = learnt.back();

      } else {
	btLevel = solver->init_level;
	p = learnt_clause[0];
	a = atom(p);
      }
      
      
#ifdef _DEBUGNOGOOD      
      std::cout << "c backtract to level " << btLevel << " and post b" << a << " =/= " << (p!=a) << std::endl;
      assert( a );      
#endif      
      
      //assert( p==a );
      
      //assert( solver->decision[cLevel] == X[a] );
      
      solver->decision[cLevel] = X[a];
      X[a]->postCut( p!=a );
      solver->backtrackLevel = btLevel;
      
      learnt_clause.clear();
      visitedUbs.clear();
      visitedLbs.clear();
      visitedDisjuncts.clear();
      disjunctToExplain.clear();
      
    }
    }
  }  
}








void ConstraintClauseBase::explainBounds()
{
  int not_finished=true, item;
  while( not_finished ) {
    if( !lbToExplain.empty() ) {
      item = lbToExplain.pop();
      explain_lowerbound( item );
    } else if( !ubToExplain.empty() ) {
      item = ubToExplain.pop();
      explain_upperbound( item );
    } else {
      not_finished=false;
    }
  }
}





void ConstraintClauseBase::extractNogood()
{

  int not_finished=true, culprit, lower_bound=-1, upper_bound=-1, 
    cLevel=*level, nbNotExplained=1, i, a, n;

  if( !conflict ) {

    //std::cout << "failed from disjunction" << std::endl;

    // get the 'failed' atom, lb and ub in the list of things to explain
    lower_bound = failure[0];
    upper_bound = failure[1];

    visitedDisjuncts.insert( failure[2] );
    if( lvl[failure[2]] < cLevel ) {
      //assert( X[failure[2]]->isGround() );
      learnt_clause.push( X[failure[2]]->equal(1) ? -failure[2] : failure[2] );
    }
    else
      disjunctToExplain.insert( failure[2] );
  } else {

    //std::cout << "failed from clause" << std::endl;

    Clause& clause = *conflict;
    n = clause.size;
    for(i=0; i<n; ++i) {
      a = atom(clause[i]);
      if( !visitedDisjuncts.member(a) ) {
	visitedDisjuncts.insert(a);
	if( lvl[a] < cLevel ) {
	  //assert( X[a]->isGround() );
	  learnt_clause.push(clause[i]);
	} else
	  disjunctToExplain.insert(a);
      }
    }
  }
   
  // now start exploring the graph from these nodes
  while( not_finished ) 
    {

#ifdef _DEBUGNOGOOD
      std::cout << std::endl << "start cycle, clause: (";
      for(int j=0; j<learnt_clause.size; ++j)
	std::cout << " " << learnt_clause[j];
      std::cout << " ). disjuncts to explain: " ;    
      disjunctToExplain.print( std::cout );
      std::cout << std::endl;
#endif 

      if( lower_bound >= 0 || upper_bound >= 0 ) {

#ifdef _DEBUGNOGOOD
	std::cout << "expand lower and upper bounds" << std::endl;
	
	if( lower_bound < 0 )
	  std::cout << "nothing to explain for the lower bound" << std::endl;
	else {
#endif 
	  
	  while( lower_bound >= 0 ) 
	    {
	      visitedLbs.insert( lower_bound );
#ifdef _DEBUGNOGOOD
	      std::cout << "lb(t" << lower_bound << ") " ;
#endif 
	      culprit = reason_lb[lower_bound];
	      if( culprit ) {
#ifdef _DEBUGNOGOOD   
		std::cout << "<- b" << culprit << " & ";
#endif
		if( !visitedDisjuncts.member(culprit) ) {
		  visitedDisjuncts.insert(culprit);
		  if( lvl[culprit] < cLevel ) {
		    //assert( X[culprit]->isGround() );
		    learnt_clause.push( X[culprit]->equal(1) ? -culprit : culprit );
		  } else
		    disjunctToExplain.insert(culprit);
		}
		lower_bound = disjunct[X[culprit]->equal(1)][culprit];
		if( visitedLbs.member( lower_bound ) ) {
		  lower_bound = -1;
		  break;
		}
	      } else lower_bound = -1;
	    }
#ifdef _DEBUGNOGOOD   
	  std::cout << std::endl;
	}
#endif 
	
#ifdef _DEBUGNOGOOD
	if( upper_bound < 0 )
	  std::cout << "nothing to explain for the upper bound" << std::endl;
	else {
#endif 
	  
	  while( upper_bound >= 0 ) {
	    visitedUbs.insert( upper_bound );
#ifdef _DEBUGNOGOOD
	    std::cout << "ub(t" << upper_bound << ") ";
#endif 
	    culprit = reason_ub[upper_bound];
	    if( culprit ) {
#ifdef _DEBUGNOGOOD   
	      std::cout << "<- b" << culprit << " & ";
#endif
	      if( !visitedDisjuncts.member(culprit) ) {
		visitedDisjuncts.insert(culprit);
		if( lvl[culprit] < cLevel )
		  learnt_clause.push( X[culprit]->equal(1) ? -culprit : culprit );
		else
		  disjunctToExplain.insert(culprit);
	      }
	      upper_bound = disjunct[X[culprit]->equal(0)][culprit];
	      if( visitedUbs.member( upper_bound ) ) {
		upper_bound = -1;
		break;
	      }
	    } else upper_bound = -1;
	  }
	  
#ifdef _DEBUGNOGOOD   
	  std::cout << std::endl;
	}
	std::cout << "after expanding the bounds: (";
	for(int j=0; j<learnt_clause.size; ++j)
	  std::cout << " " << learnt_clause[j];
	std::cout << " ) - " ;    
	disjunctToExplain.print( std::cout );
	std::cout << std::endl;
#endif
      }


      if( disjunctToExplain.size > (unsigned int)nbNotExplained ) {
	culprit = disjunctToExplain.pop();
#ifdef _DEBUGNOGOOD   
      	std::cout << "b" << culprit << " <- ";
#endif
	if( reason[culprit] ) {
#ifdef _DEBUGNOGOOD
	  printClause( std::cout, reason[culprit] );
	  std::cout << std::endl;
#endif
	  Clause& clause = *(reason[culprit]);
	  n = clause.size;
	  for(i=0; i<n; ++i) 
	    {
	      a = atom(clause[i]);
	      if( !visitedDisjuncts.member(a) ) 
		{
		  visitedDisjuncts.insert(a);
		  if( lvl[a] < cLevel )
		    learnt_clause.push(clause[i]);
		  else
		    disjunctToExplain.insert(a);
		}
	    }
	} else {
	  // why has this disjunct set to its value?
	  // 1/ find out its polarity
	  // 	int d_lb, d_ub, disjunct_polarity = //(d == polarity[d]);
	  // 	  X[culprit]->value();
	  // 	if( disjunct_polarity ) { // it is 1, hence X1 < X0
	  // 	  // therefore the lb of X0 and the ub of X1 are responsible
	  // 	  d_lb = reason_lit[0][culprit];
	  // 	  d_ub = reason_lit[1][culprit];
	  // 	  if( d_lb ) lower_bound = disjunct[X[d_lb]->equal(1)][d_lb];
	  // 	  if( d_ub ) upper_bound = disjunct[X[d_ub]->equal(0)][d_ub];
	  // 	} else { // it is 1, hence X0 < X1
	  // 	  // therefore the lb of X1 and the ub of X0 are responsible
	  // 	  d_lb = reason_lit[1][culprit];
	  // 	  d_ub = reason_lit[0][culprit];
	  // 	  if( d_lb ) lower_bound = disjunct[X[d_lb]->equal(1)][d_lb];
	  // 	  if( d_ub ) upper_bound = disjunct[X[d_ub]->equal(0)][d_ub];
	  // 	}
	  int d_lb, d_ub;
	  // therefore the lb of X0 and the ub of X1 are responsible
	  d_lb = reason_lit[0][culprit];
	  d_ub = reason_lit[1][culprit];
	  if( d_lb ) lower_bound = disjunct[X[d_lb]->equal(1)][d_lb];
	  if( d_ub ) upper_bound = disjunct[X[d_ub]->equal(0)][d_ub];
	  if( d_lb + d_ub ) {
#ifdef _DEBUGNOGOOD
	    if(d_lb) 
	      std::cout << "b" << d_lb << " & lb(t" << lower_bound << ")";
	    if(d_lb && d_ub)
	      std::cout << " & ";
	    if(d_ub)
	      std::cout << "b" << d_ub << " & ub(t" << upper_bound << ")" ;
	    std::cout << std::endl;
#endif
	    if( d_lb && !visitedDisjuncts.member(d_lb) ) {
	      visitedDisjuncts.insert(d_lb);
	      if( lvl[d_lb] < cLevel )
		learnt_clause.push( X[d_lb]->equal(1) ? -d_lb : d_lb );
	      else
		disjunctToExplain.insert(d_lb);
	    }
	    if( d_ub && !visitedDisjuncts.member(d_ub) ) {
	      visitedDisjuncts.insert(d_ub);      
	      if( lvl[d_ub] < cLevel )
		learnt_clause.push( X[d_ub]->equal(1) ? -d_ub : d_ub );
	      else
		disjunctToExplain.insert(d_ub);
	    }
	  } else {
#ifdef _DEBUGNOGOOD
	    std::cout << "that was the last decision" << std::endl;
#endif
	    --nbNotExplained;
	    learnt_clause.push( (X[culprit]->equal(1) ? -culprit : culprit) );
	  }
	}
      } else {

#ifdef _DEBUGNOGOOD
	std::cout << "the analysis is finished ";
#endif

	if(nbNotExplained) {
#ifdef _DEBUGNOGOOD
	  std::cout << ", still need to add the last unexplained disjunct" << std::endl;
#endif
	  culprit = disjunctToExplain.pop();
	  learnt_clause.push( (X[culprit]->equal(1) ? -culprit : culprit) );
	} 
#ifdef _DEBUGNOGOOD
	else std::cout << std::endl;
#endif
	not_finished = false;
      }
  

    }
    

#ifdef _DEBUGNOGOOD
  for(int i=0; i<learnt_clause.size; ++i) {
    int p = learnt_clause[i];
    int a = atom(p);
    
    //if( X[a]->isGround() ) {
    if( (p > 0) != (X[a]->equal(0)) ) {
      std::cout << p << " ";
      X[a]->print( std::cout );
      std::cout << std::endl;
    }
    assert( (p > 0) == (X[a]->equal(0)) );
    //}
  }
#endif
  //exit(0);
}




void ConstraintClauseBase::exploreImplicationGraph()
{

  int i=0, item, finished=false;

  if( !conflict ) {
    // get the 'failed' atom, lb and ub in the list of things to explain
    visitedLbs.insert( failure[0] );
    lbToExplain.insert( failure[0] );
    
    visitedUbs.insert( failure[1] );
    ubToExplain.insert( failure[1] );
    
    visitedDisjuncts.insert( failure[2] );
    disjunctToExplain.insert( failure[2] );
  } else {
    Clause& clause = *conflict;
    unsigned int i, n = clause.size;
    Atom a;
    for(i=0; i<n; ++i) {
      a = atom(clause[i]);
      if( !visitedDisjuncts.member(a) ) {
	visitedDisjuncts.insert(a);      
	disjunctToExplain.insert(a);
      }
    }
  }

#ifdef _DEBUGNOGOOD      
  std::cout << "UB: " ;
  visitedUbs.print( std::cout );
  std::cout << std::endl;

  std::cout << "LB: " ;
  visitedLbs.print( std::cout );
  std::cout << std::endl;

  std::cout << "DJ: " ;
  visitedDisjuncts.print( std::cout );
  std::cout << std::endl;
#endif

  i=0;
  int head[3];
  head[0] = 0;
  head[1] = 0;
  head[2] = 0;
  // now start exploring the graph from these nodes
  while( !finished ) {

    // pick something yet to be explained
    if( !lbToExplain.empty() ) {
      item = lbToExplain.pop();
      explain_lowerbound( item );
    } else if( !ubToExplain.empty() ) {
      item = ubToExplain.pop();
      explain_upperbound( item );
    } else if( !disjunctToExplain.empty() ) {
      item = disjunctToExplain.pop();
      explain_disjunct( item );
    } else {
      finished = true;
    }

//     // pick something yet to be explained
//     if( disjunctToExplain.size > head[0] ) {
//       item = disjunctToExplain[head[0]++];
//       explain_disjunct( item );
//     } else if( lbToExplain.size > head[1] ) {
//       item = lbToExplain[head[1]++];
//       explain_lowerbound( item );
//     } else if( ubToExplain.size > head[2] ) {
//       item = ubToExplain[head[2]++];
//       explain_upperbound( item );
//     } else {
//       finished = true;
//       ubToExplain.clear();
//       lbToExplain.clear();
//       disjunctToExplain.clear();
//     }

  }
}




void ConstraintClauseBase::explain_upperbound( const int task )
{
  // what are the other variables responsible for the upper bound of 'task' to be low?

#ifdef _DEBUGNOGOOD
  std::cout << "\tub(t" << task << ") <- ";
#endif

  // first, the disjunct corresponding to the constraint that last pruned task's ub
  int bound, culprit = reason_ub[task];
  if( culprit ) {

#ifdef _DEBUGNOGOOD
    std::cout << "b" << ( X[culprit]->equal(1) ? culprit : -culprit ) << " & ";
#endif

    if( !visitedDisjuncts.member(culprit) ) {
      visitedDisjuncts.insert(culprit);
      disjunctToExplain.insert(culprit);   
    }

    bound = disjunct[X[culprit]->equal(0)][culprit];

#ifdef _DEBUGNOGOOD
      std::cout << "ub(t" << bound << ")" << std::endl;
#endif 

    if( !visitedUbs.member(bound) ) {
      visitedUbs.insert(bound);      
      ubToExplain.insert(bound);
    }
  } 
#ifdef _DEBUGNOGOOD
  else {
    std::cout << "null" << std::endl; 
  }
#endif

}

void ConstraintClauseBase::explain_lowerbound( const int task )
{
#ifdef _DEBUGNOGOOD
  std::cout << "\tlb(t" << task << ") <- ";
#endif

  // what are the other variables responsible for the lower bound of 'task' to be high?

  // first, the disjunct corresponding to the constraint that last pruned task's lb
  int bound, culprit = reason_lb[task];
  if( culprit ) {

#ifdef _DEBUGNOGOOD   
      std::cout << "b" << ( X[culprit]->equal(1) ? culprit : -culprit ) << " & ";
#endif

    if( !visitedDisjuncts.member(culprit) ) {
      visitedDisjuncts.insert(culprit);
      disjunctToExplain.insert(culprit);
    }
    bound = disjunct[X[culprit]->equal(1)][culprit];

#ifdef _DEBUGNOGOOD
    std::cout << "lb(t" << bound << ")" << std::endl;
#endif 

    if( !visitedLbs.member(bound) ) {
      visitedLbs.insert(bound);      
      lbToExplain.insert(bound);
    }
  } 

#ifdef _DEBUGNOGOOD
  else {
    std::cout << "null" << std::endl; 
  }
#endif

}


void ConstraintClauseBase::explain_disjunct( const int d )
{

#ifdef _DEBUGNOGOOD
  std::cout << "\tb" << ( X[d]->equal(1) ? d : -d ) << " <- " ;
#endif

  if( isDecision(d) ) {

#ifdef _DEBUGNOGOOD
    std::cout << "null" << std::endl;
#endif

    learnt_clause.push(-polarity[d]);
  } else if( reason[d] ) {

#ifdef _DEBUGNOGOOD
    printClause( std::cout, reason[d] );
    std::cout << std::endl;
#endif

    Clause& clause = *(reason[d]);
    unsigned int i, n = clause.size;
    Atom a;
    for(i=0; i<n; ++i) {
      a = atom(clause[i]);
      if( !visitedDisjuncts.member(a) ) {
	visitedDisjuncts.insert(a);      
	disjunctToExplain.insert(a);
      }
    }

  } else {

    // why has this disjunct set to its value?
    
    // 1/ find out its polarity
    int lb=0, ub=0, d_lb, d_ub, disjunct_polarity = //(d == polarity[d]);
      X[d]->value();
      
    
    if( disjunct_polarity ) { // it is 1, hence X1 < X0
      // therefore the lb of X0 and the ub of X1 are responsible
      d_lb = reason_lit[0][d];
      d_ub = reason_lit[1][d];
      if( d_lb ) lb = disjunct[X[d_lb]->equal(1)][d_lb];
      if( d_ub ) ub = disjunct[X[d_ub]->equal(0)][d_ub];
    } else { // it is 1, hence X0 < X1
      // therefore the lb of X1 and the ub of X0 are responsible
      d_lb = reason_lit[1][d];
      d_ub = reason_lit[0][d];
      if( d_lb ) lb = disjunct[X[d_lb]->equal(1)][d_lb];
      if( d_ub ) ub = disjunct[X[d_ub]->equal(0)][d_ub];
    }

#ifdef _DEBUGNOGOOD
    if(d_lb) 
      std::cout << "b" << ( X[d_lb]->equal(1) ? d_lb : -d_lb ) << " & lb(t" << lb << ")";
    if(d_lb && d_ub)
      std::cout << " & ";
    if(d_ub)
      std::cout << "b" << ( X[d_ub]->equal(1) ? d_ub : -d_ub ) << " & ub(t" << ub << ")" ;
    std::cout << std::endl;
#endif

    
    if( d_lb && !visitedDisjuncts.member(d_lb) ) {
      visitedDisjuncts.insert(d_lb);      
      disjunctToExplain.insert(d_lb);
    }
    
    if( d_ub && !visitedDisjuncts.member(d_ub) ) {
      visitedDisjuncts.insert(d_ub);      
      disjunctToExplain.insert(d_ub);
    }
    
    if( d_lb && !visitedLbs.member(lb) ) {
      visitedLbs.insert(lb);      
      lbToExplain.insert(lb);
    }
    
    if( d_ub && !visitedUbs.member(ub) ) {
      visitedUbs.insert(ub);
      ubToExplain.insert(ub);
    }
  }
}



int ConstraintClauseBase::explain_disjunct2( const int d, int& count )
{

#ifdef _DEBUGNOGOOD
  std::cout << "\t(" << lvl[d] << ") b" << d << " <- " ;
#endif

  int res = 0;
  if( lvl[d] < *level ) {
#ifdef _DEBUGNOGOOD
    std::cout << "stop" << std::endl;
#endif

    learnt_clause.push(-polarity[d]);
  } else if( count > 0 ) {
    --count;
#ifdef _DEBUGNOGOOD
    std::cout << "quotat" << std::endl;
#endif

    learnt_clause.push(-polarity[d]);
  } else if( reason[d] ) {

#ifdef _DEBUGNOGOOD
    printClause( std::cout, reason[d] );
    std::cout << std::endl;
#endif

    Clause& clause = *(reason[d]);
    unsigned int i, n = clause.size;
    Atom a;
    for(i=0; i<n; ++i) {
      a = atom(clause[i]);
      if( !visitedDisjuncts.member(a) ) {
	visitedDisjuncts.insert(a);      
	disjunctToExplain.insert(a);
      }
    }

  } else {

    // why has this disjunct set to its value?
    
    // 1/ find out its polarity
    int lb=0, ub=0, d_lb, d_ub, disjunct_polarity = //(d == polarity[d]);
      X[d]->value();
      
    
    if( disjunct_polarity ) { // it is 1, hence X1 < X0
      // therefore the lb of X0 and the ub of X1 are responsible
      d_lb = reason_lit[0][d];
      d_ub = reason_lit[1][d];
      if( d_lb ) lb = disjunct[X[d_lb]->equal(1)][d_lb];
      if( d_ub ) ub = disjunct[X[d_ub]->equal(0)][d_ub];
    } else { // it is 1, hence X0 < X1
      // therefore the lb of X1 and the ub of X0 are responsible
      d_lb = reason_lit[1][d];
      d_ub = reason_lit[0][d];
      if( d_lb ) lb = disjunct[X[d_lb]->equal(1)][d_lb];
      if( d_ub ) ub = disjunct[X[d_ub]->equal(0)][d_ub];
    }

#ifdef _DEBUGNOGOOD
    if(d_lb) 
      std::cout << "b" << d_lb << " & lb(t" << lb << ")";
    if(d_lb && d_ub)
      std::cout << " & ";
    if(d_ub)
      std::cout << "b" << d_ub << " & ub(t" << ub << ")" ;
    std::cout << std::endl;
#endif

    
    if( d_lb && !visitedDisjuncts.member(d_lb) ) {
      visitedDisjuncts.insert(d_lb);      
      disjunctToExplain.insert(d_lb);
    }
    
    if( d_ub && !visitedDisjuncts.member(d_ub) ) {
      visitedDisjuncts.insert(d_ub);      
      disjunctToExplain.insert(d_ub);
    }
    
    if( d_lb && !visitedLbs.member(lb) ) {
      res = 1;
      visitedLbs.insert(lb);      
      lbToExplain.insert(lb);
    }
    
    if( d_ub && !visitedUbs.member(ub) ) {
      res = 1;
      visitedUbs.insert(ub);
      ubToExplain.insert(ub);
    }
  }

  return res;
}



void ConstraintClauseBase::update()
{

  // update the 'lvl' of the literals currently on the propag stack
  MistralGacList<VariableInt*>& gstack = solver->gacvarstack;
  if( !gstack.empty() ) {
    /// BUGGY!
#ifdef _DEBUGNOGOOD
    std::cout << "update " ;
#endif
    int l = *level, q = gstack.head;
    do {
      q = gstack.next[q];

#ifdef _DEBUGNOGOOD
      std::cout << " " << (q+1);
#endif

      if( q < arity ) {
	//std::cout << "(+) ";
	//addLiteral( (q+1), (X[q+1]->equal(0) ? -(q+1) : (q+1)) );
	lvl[q+1] = l;
      } 

#ifdef _DEBUGNOGOOD
      else
	std::cout << " end" ;
#endif


    } while( q != gstack.tail );
#ifdef _DEBUGNOGOOD
    std::cout << std::endl;
#endif
  }

  // update the valuation of literals in the assumption stack
  int i=assumptions.size, cLevel=*level;
  Atom a;
  while( i-- ) { 

    a = assumptions[i];

    //std::cout << "update " << a << " ";

    if( lvl[a] < cLevel ) break;
    if( !X[a]->isGround() ) {
      X[a]->setDomain( a==polarity[a] );
      //X[a]->print( std::cout );
    }
    
    //std::cout << std::endl;
  }
  //std::cout << std::endl;
}

void ConstraintClauseBase::addDelayedClause(VariableInt **x, const int n, int *l) {
  VariableInt** scp = new VariableInt*[n];
  int* lits = new int[n];
  for(int i=0; i<n; ++i) {
    scp[i] = x[i];
    lits[i] = l[i]; 
  }
  cl_scope.push(scp);
  cl_lits.push(lits);
  cl_size.push(n);
}

void ConstraintClauseBase::initialise()
{
  if(!scope) {
    params.forgetfulness = 0.9;
    solver->sat = this;

    int i, j, m = solver->variables.size;
    int atom_id[m];
    std::fill(atom_id, atom_id+m, -1);
  
    if( disjuncts.empty() ) {
      initScope( solver->variables.stack_, solver->variables.size, VALUETRIGGER, 0 );
      if( solver->constraints.size == 1 )
	mode = SAT;
	//mode = CSP;
      else
	mode = CSP;

      m = solver->length;
      for(i=0; i<m; ++i)
	atom_id[scope[i]->id] = i+1;
    
    } else {
   
      arity = disjuncts.size;
      VariableInt *bool_vars[arity];
      for(i=0; i<arity; ++i)
	bool_vars[i] = disjuncts[i]->scope[2];

      initScope( bool_vars, arity, VALUETRIGGER, 0 );

    }

    init( arity, base.size );
    X = scope-1;
    solver->learners.push( new WeighterSAT( solver, this ) );

    if( solver->constraints.size != 1 ) {
      //if( solver->constraints.size > m+1 )
      if(disjuncts.empty())
	mode = CSP;
      else {
	mode = DTP;
      
	nTasks=0;
	reason_lit[0]  = new int[arity+1];
	reason_lit[1]  = new int[arity+1];
	std::fill( reason_lit[0], reason_lit[0]+arity+1, 0 );
	std::fill( reason_lit[1], reason_lit[1]+arity+1, 0 );
	disjunct[0] = new int[arity+1];
	disjunct[1] = new int[arity+1];
	for(j=0; j<arity; ++j) {
	  disjuncts[j]->atom_id = j+1;
	  disjuncts[j]->reason = &(reason[j+1]);
	
	  i = disjuncts[j]->scope[0]->id;
	  if( atom_id[i] < 0 ) atom_id[i] = nTasks++;
	  disjunct[0][j+1] = atom_id[i];
	
	  i = disjuncts[j]->scope[1]->id;
	  if( atom_id[i] < 0 ) atom_id[i] = nTasks++;
	  disjunct[1][j+1] = atom_id[i];
       
	  disjuncts[j]->task_id[0] = disjunct[0][j+1];
	  disjuncts[j]->task_id[1] = disjunct[1][j+1];

#ifdef _DEBUGNOGOOD
	  std::cout << "b" << (j+1) ;
	  if( j+1 < 10 )
	    std::cout << " ";	
	  std::cout << ":  t" << disjunct[0][j+1] << ":" ;
	  disjuncts[j]->scope[0]->printDomain( std::cout );
	  std::cout << " (" << disjuncts[j]->duration[0] << ") ";
	  std::cout << " . t" << disjunct[1][j+1] << ":" ;
	  disjuncts[j]->scope[1]->printDomain( std::cout );
	  std::cout << " (" << disjuncts[j]->duration[1] << ") ";
	  std::cout << std::endl ;

#endif	

	}      

	the_task = new VariableInt*[nTasks];
	reason_lb = new ReversibleNum<int>[nTasks];
	reason_ub = new ReversibleNum<int>[nTasks];
	for(j=0; j<nTasks; ++j) {
	  solver->binds( reason_lb[j] );
	  solver->binds( reason_ub[j] );
	  reason_lb[j].setValue( 0 );
	  reason_ub[j].setValue( 0 );
	}
      
	j=arity;
	while( j-- ) {
	  disjuncts[j]->reason_for_lb_increase[0] = &(reason_lb[disjunct[0][j+1]]);
	  disjuncts[j]->reason_for_lb_increase[1] = &(reason_lb[disjunct[1][j+1]]);
	  disjuncts[j]->reason_for_ub_decrease[0] = &(reason_ub[disjunct[0][j+1]]);
	  disjuncts[j]->reason_for_ub_decrease[1] = &(reason_ub[disjunct[1][j+1]]);
	  disjuncts[j]->reason_for_choice[0] = &(reason_lit[0][j+1]);
	  disjuncts[j]->reason_for_choice[1] = &(reason_lit[1][j+1]);

	  the_task[disjunct[0][j+1]] = disjuncts[j]->scope[0];
	  the_task[disjunct[1][j+1]] = disjuncts[j]->scope[1];
	}

	visitedUbs.init(0, nTasks-1, BitSet:: empt);
	visitedLbs.init(0, nTasks-1, BitSet:: empt);
	visitedDisjuncts.init(1, arity, BitSet:: empt);
      
	lbToExplain.init(0, nTasks-1, 0);
	ubToExplain.init(0, nTasks-1, 0);
	disjunctToExplain.init(1, arity, 0);
      
      }
    }

      m = cl_scope.size;
      if(m) {
	Vector< Literal > tmp_cl;
	for(i=0; i<m; ++i) {
	  for(j=0; j<cl_size[i]; ++j) {
	    if(cl_lits[i][j])
	      tmp_cl.push(cl_scope[i][j]->id+1);
	    else
	      tmp_cl.push(-(cl_scope[i][j]->id+1));
	  }

	  addClause( base, tmp_cl, stats.base_avg_size );
	  addOriginalClause( tmp_cl );
	  tmp_cl.clear();
	}
      }

  }
}

void ConstraintClauseBase::restore()
{
  if(*level > 1) {
    //std::cout << (*level - 2) << std::endl;
    //backtrackTo( (*level)-(1+solver->init_level) );
    //backtrackTo( (*level)-2 );
    backtrackTo( (*level)-2 );
  }
}

void ConstraintClauseBase::print(std::ostream& o) const 
{

  o << "CLAUSE BASE" << std::endl;
  printClauses(o);

}


