
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

#include <mistral_csp.h>
#include <mistral_sol.h>
#include <mistral_con.h>
#include <mistral_set.h>

using namespace Mistral;

void checkCons( MistralList<Constraint*>& l, VariableInt* x )
{
  Constraint* con;
  int i = 100;
  MistralNode<Constraint*>* nd = &(l.head);
  while( nextNode(nd) && --i ) {
    con = nd->elt;
    bool isIn = false;
    for(int i=0; !isIn && i<con->arity; ++i)
      isIn = (con->_scope[i] == x);
    assert( isIn );
  }
  assert( i );
}


/**********************************************
 * Constraint
 **********************************************/

void Constraint::relax()
{
  //solver->save( this );
  store->push( this );
  for(int i=0; i<arity; ++i)
    elements[i].erase();
}

void Constraint::restore()
{
  for(int i=0; i<arity; ++i)
    triggers[i]->insert( elements[i] );
}

void Constraint::initScope(VariableInt ** scp, const int n, const int t, const int w) 
{
  weight = w;
  arity = n;
  _scope = new VariableInt*[arity];
  triggers = new MistralList<Constraint*>*[arity];
  elements = new MistralNode<Constraint*>[arity+2];
  
  //std::cout << "initScope " << arity << std::endl;


  int i = arity+2, refs=false;
  while( i-- ) {
    elements[i].elt = this;
    elements[i].index = i;

    //std::cout << " " << i ;
  }
  //std::cout << std::endl;

  i = arity;
  while( i-- ) {
    refs |= (scp[i]->getType() == VariableInt::VIRTUAL);
    _scope[i] = scp[i]->getVar();    
    switch(t) {
    case DOMAINTRIGGER:
      { 
	triggers[i] = _scope[i]->triggerOnDomain();
	break;
      }
    case RANGETRIGGER: 
      {
	triggers[i] = _scope[i]->triggerOnRange();
	break;
      }
    case VALUETRIGGER: 
      {
	triggers[i] = _scope[i]->triggerOnValue();
	break;
      }
    }
  }

  if(refs)
    {
      scope = new VariableInt*[arity];
      i = arity;
      while( i-- ) 
	scope[i] = scp[i];
    }
  else scope = _scope;

  /// this is in case some variables are duplicated in _scope
  /// the watches should always be different variables
  watch[1] = _scope[arity-1];
  i = arity-1;
  while( i-- && _scope[i] == watch[1] );
  if( i < 0 ) i=0;
  watch[0] = _scope[i];

  /// keep that test for now. If we couldn't find two different 
  /// variables to use as watchers, then we hope that the
  /// sole variable is a constant (apparently it is not a problem in that case)
  //assert( watch[0] != watch[1] || watch[0]->isGround() );


  lastlinked = arity-1;

  //elements[arity].type = 1;
  //elements[arity+1].type = 1;
  
  if( elements[arity].elt ) 
    watch[0]->watchedConstraints.insert( elements[arity] );
  if( elements[arity+1].elt ) 
    watch[1]->watchedConstraints.insert( elements[arity+1] );
  
  i = 0;
  while( i < arity ) {
    if(triggers[i]) {
      triggers[i]->insert( elements[i] );
      ++_scope[i]->degree;
      _scope[i]->weight += weight;
      
      //_scope[i]->print( std::cout );
      //std::cout << " " << elements[i].index << "  ";
    }    
    ++i;
  }
  //std::cout << std::endl;

  //printIndices();
}


// void Constraint::printIndices()
// {
//   for(int i=0; i<arity; ++i) {
//     _scope[i]->printshort( std::cout );
//     std::cout << ": " << i << " /";
//     MistralNode<Constraint*> *nd = _scope[i]->constraintsOnValue();
//     while( nextNode(nd) ) 
//       std::cout << " " << nd->index;
//     std::cout << std::endl;
//   }
// }

Constraint::Constraint(Solver *s,
		       VariableInt ** scp, 
		       const int n,
		       const int t,
		       const int w) 
//  : ReversibleObj( s )
{
  s->binds( *this );
  //solver = s;
  delayed = false;
  initScope( scp, n, t, w );
  supp = new int[n];
  supports = NULL;
  id = s->constraints.size;
  s->constraints.push( this );
}

Constraint::Constraint(Solver *s) 
//  : ReversibleObj( s )
{
  s->binds( *this );
  //solver = s;
  delayed = false;
  supp = NULL;
  supports = NULL;
  id = s->constraints.size;
  s->constraints.push( this );
}

Constraint::~Constraint() 
{
  if( supports ) {
    for(int x=0; x<arity; ++x) {
      for(int v=scope[x]->minCapacity(); v<scope[x]->maxCapacity(); ++v)
	delete [] supports[x][v];
      supports[x] += scope[x]->minCapacity();
      delete [] supports[x];
    }
    delete [] supports;
  }
  delete [] supp;
  if( scope != _scope )
    delete [] scope;
  delete [] _scope;
  delete [] elements;
  delete [] triggers;
}

void Constraint::useResidualSupports()
{  
  int x, v, vfill = NOVAL;
  supports = new int**[arity];
  for(x=0; x<arity; ++x) {
    supports[x] = new int*[scope[x]->maxCapacity()-scope[x]->minCapacity()];
    supports[x] -= scope[x]->minCapacity();
    for(v=scope[x]->minCapacity(); v<scope[x]->maxCapacity(); ++v) {
      if( scope[x]->contain(v) ) {
	supports[x][v] = new int[arity];
	std::fill(supports[x][v],
		  supports[x][v]+arity,
		  vfill);
      } else supports[x][v] = NULL;
    }
  }
}

bool Constraint::firstSupport(const int vri, 
			      const int vli) 
{
  int j;
  if( supports && supports[vri][vli][0] != NOVAL ) {
    j=arity;
    while( j-- ) 
      if( vri != j )
	if (!scope[j]->contain( supports[vri][vli][j] )) break;
    if( j < 0 ) 
      return true;
  } 
  j=arity;
  while( j-- )
    supp[j] = scope[j]->first();
  supp[vri] = vli; 

  return false;
}

bool Constraint::findSupport(const int vri, const int vli) 
{
  int i=arity, vali;
  bool found=false;
  // sol is initialized: either to the value 
  // a variable is already assigned to
  // or to the first value in its domain
  while(i >= 0) {
    // check this assignment
    if( !check( supp ) ) {
      found=true;
      if( supports ) {
	vali = arity;
	while( vali-- )
	  supports[vri][vli][vali] = supp[vali];
      }
      break;
    }
 
    // try to assign more things
    // find the last var whose domain we have not exhausted
    --i;
    while( i >= 0 ) {
      if( i == vri || scope[i]->isGround() ) {
	--i;
	continue;
      }
      if( scope[i]->setNext( supp[i] ) )
      	break;
      else
	supp[i] = scope[i]->first();
      --i;
    }
    if( i >= 0 )
      i = arity;
  } 
  return found;
}

bool Constraint::propagate()
{
  int i, consistent=1, evt = ( Constraint::RANGETRIGGER );
  //int i, consistent=1, evt = ( Constraint::VALUETRIGGER );
  for( i=0; consistent && i<arity; ++i )
    consistent = propagate( i, evt );
  return consistent;
} 

bool Constraint::propagate(const int changedIdx, const int e) 
{
  bool consistent = true;
  for( int i=0; consistent && i<arity; ++i )
    if( i != changedIdx && _scope[i]->isLinked() ) { 
      consistent = scope[i]->revise( this, i );
    }
  
  return consistent;
}

bool Constraint::isLinked()
{
  return ( watch[0]->isLinked() && watch[1]->isLinked() );
}

inline void Constraint::stopWatching( VariableInt* x ) 
{
  // x became unlinked.
  // we thus need to find another variable to replace x
  // lastlinked is set to the index of the new watcher if
  // one could be found
  // or to the index of the last remaining watcher otherwise

  // x == watch[idx]
  int i=arity, idx = (watch[0] != x);
  // find a new variable for replacing the first one (x)
  while( i-- ) {
    if( _scope[i] == watch[1-idx] )
      lastlinked = i;
    else if( _scope[i] != watch[idx] &&
	     _scope[i]->isLinked()
	     ) {

      // stop being watched by that variable and get a new watcher
      elements[arity+idx].erase();
      watch[idx] = _scope[i];
      _scope[i]->watchedConstraints.insert( elements[arity+idx] );

      lastlinked = i;
      break;
    }
  }
  // if a new variable could not be found we remove this constraint
  // from the list of the last linked variable
  if( i < 0 ) {
    oldweight = weight;
    _scope[lastlinked]->weight -= weight;
    --_scope[lastlinked]->degree;
    elements[lastlinked].erase();
    elements[arity+1-idx].erase();
  } 
}

inline void Constraint::resumeWatching( VariableInt* x ) 
{
  if( x != _scope[lastlinked] ) {
    int idx = (watch[0] != x);
    watch[1-idx]->watchedConstraints.insert( elements[arity+1-idx] );
    triggers[lastlinked]->insert( elements[lastlinked] );
    ++_scope[lastlinked]->degree;
    _scope[lastlinked]->weight += oldweight; 
  }
}

void Constraint::printshort(std::ostream& o) const 
{
  o << "c" << id;
}

void Constraint::print(std::ostream& o) const 
{
  o << "(" ;
  for(int i=0; i<arity-1; ++i) {
    scope[i]->printshort( o );
    o << ", ";
  }
  scope[arity-1]->printshort( o );
  o << ")";
}

bool Constraint::constrains( VariableInt *x ) const 
{
  bool cstr = false;
  for(int i=0; !cstr && i<arity; ++i) 
    cstr = ( _scope[i] == x ) ;
  return cstr;
}

/**********************************************
 * VariableInt
 **********************************************/

VariableInt::VariableInt() 
{
  isWord = false;
  id = -1;
  degree = 0;
  weight = 0;
  branch = NULL;
} 

VariableInt::VariableInt( Solver *s ) 
  : ACQueue( &(s->gacvarstack) )
{
  isWord = false;
  id = -1;
  degree = 0;
  weight = 0;
  branch = NULL;
  valueTrigger.head.next  = &(rangeTrigger.head);
  rangeTrigger.head.next  = &(domainTrigger.head);
  domainTrigger.head.pred = &(rangeTrigger.head);
  rangeTrigger.head.pred  = &(valueTrigger.head);
  solver = s;
  id = solver->variables.size;
  solver->variables.push( this );
}

VariableInt::~VariableInt() 
{
  delete branch;
}

void VariableInt::unLink( ) 
{
//   if(id == 59) {
//     (*seqIdx)->print(std::cout);
//     std::cout << std::endl;
//     (**seqBeg)->print(std::cout);
//     std::cout << std::endl;
//   }

//   std::cout << "unlink " ;
//   printshort(std::cout);
//   std::cout << std::endl;

  *seqIdx = **seqBeg;
  (*seqIdx)->seqIdx = seqIdx;
  seqIdx = *seqBeg;
  **seqBeg = this;
  ++(*seqBeg);

  MistralNode<Constraint*> *nd, *nxt = &(watchedConstraints.head);
  nd = nextNode(nxt);
  while( nd ) {
    nextNode(nxt);
    nd->elt->stopWatching( this );
    nd = nxt;
  }
}

void VariableInt::link( ) 
{
  
  --(*seqBeg);

  MistralNode<Constraint*> *nd = &(watchedConstraints.head);
  while( nextNode(nd) ) 
    nd->elt->resumeWatching( this );
  
}

void VariableInt::postCut( const int p )
{
  branch->postCut( p );
}

int VariableInt::getValue() const
{
  return solver->solution[id];
}

int VariableInt::getMin() const
{
  return solver->solution[id];
}

int VariableInt::getMax() const
{
  return solver->max_solution[id];
}


/**********************************************
 * MistralGraph
 **********************************************/

/// A directed graph implementation.
MistralGraph::MistralGraph( VariableInt **vars, const int l ) 
{
  int i=l, j, lb=NOVAL, ub=-1;
  while( i-- ) {
    if( vars[i]->min() < lb ) lb = vars[i]->min();
    if( vars[i]->max() > ub ) ub = vars[i]->max();
  }
  n = l+ub-lb+1;

  first = new int *[n];
  last  = new int *[n];    
  node  = new int**[n];

  for(i=0; i<n; ++i) {
    first[i] = new int[n];
    last[i]  = first[i];
    node [i] = new int*[n];
    j = n;
    while( j-- ) node[i][j] = NULL;
  }
  for(i=0; i<l; ++i) {
    j = vars[i]->first();
    do add(i, j);
    while( vars[i]->setNext(j) );
  }  
}

void MistralGraph::getAllCliques( std::vector< std::vector< int > >& cliques, const int limit )
{
  BitSet cand  (0, n-1, BitSet::full);
  BitSet K     (0, n-1, BitSet::empt);
  BitSet gammaK(0, n-1, BitSet::full);
  for(int i=0; i<n; ++i)
    {
      if( !degree(i) )
	{
	  cand.erase(i);
	  gammaK.erase(i);
	}
    }
  BronKerbosch( cand, K, gammaK, cliques, limit );
}

int MistralGraph::getAndRemoveMaxCand(BitSet& cand) {
  int index = -1; 
  if (!cand.empty()) {
    int i = cand.min(), max = 0; 
    do {
      if(degree(i)>max) {
	max = degree(i); 
	index = i; 
      } 
      i = cand.next(i);
    } while( i != NOVAL );
    cand.erase(index);
  }
  return index; 
}

void MistralGraph::removeNeighbours(const int x, BitSet& neighbours) 
{
  int *varit;
  for(varit=first[x]; varit!=last[x]; ++varit) 
    neighbours.erase( *varit );
}

void MistralGraph::updateGammaK(const int x, BitSet& res) 
{
  int i=res.min(), *varit, isntIn;
  do {
    isntIn=true;
    for(varit=first[x]; isntIn && varit!=last[x]; ++varit)
      isntIn = (*varit != i);
    if( isntIn ) res.erase(i);
    i = res.next(i);
  } while( i != NOVAL );
  res.erase(x);
}

void MistralGraph::BronKerbosch( BitSet& cand,
				 BitSet& K,
				 BitSet& gammaK,
				 std::vector< std::vector< int > >& cliques, 
				 const int limit )
{
  if( (int)(cliques.size()) < limit ) {
    BitSet updatedCand(cand);
    BitSet cloneK(K);
    BitSet cloneGammaK(gammaK);
    
    int x = getAndRemoveMaxCand(updatedCand);	
    removeNeighbours(x, updatedCand);
    if( !updatedCand.empty() )
      BronKerbosch(updatedCand, cloneK, cloneGammaK, cliques, limit);
    
    if( (int)(cliques.size()) < limit ) {
      BitSet updatedGammaK(gammaK);  
    
      K.insert(x);
      updateGammaK(x, updatedGammaK);
      if( !updatedGammaK.empty() )
	BronKerbosch(updatedGammaK, K, updatedGammaK, cliques, limit);
      else storeClique( K, cliques );
    }
  }
}

void MistralGraph::storeClique(BitSet& K, std::vector< std::vector< int > >& cliques) 
{
  std::vector<int> clique;
  int i=K.min();
  do { 
    //    std::cout << i << " ";
    clique.push_back(i);
    i = K.next(i);
  } while( i != NOVAL );
  //  std::cout << std::endl;
  cliques.push_back( clique );
}


void MistralGraph::getSomeCliques( std::vector< std::vector< int > >& cliques )
{
  int x, *varit, isNeighbor, i;

  BitSet visited( 0, n-1, BitSet::empt );
  Vector<int> vertice;
  for(x=0; x<n; ++x)
    vertice.push( x );

  while( !vertice.empty() )
    {
      // arbitrarily pick a node x
      vertice.pop( x );
      if( !visited.member( x ) ) {
	visited.insert( x );
	std::vector<int> neighbors;
	neighbors.push_back( x );

	// add all neighbors as long as they form a clique
	for(varit=first[x]; varit!=last[x]; ++varit) {
	  isNeighbor = true;
	  i = neighbors.size();
	  while( isNeighbor && --i )
	    isNeighbor = isIn( *varit, neighbors[i]);
	  if( isNeighbor ) {
	    neighbors.push_back( *varit );
	    visited.insert( *varit );
	  } 
	}
	
	cliques.push_back( neighbors );
      } 
    }
}

MistralGraph::MistralGraph( const int l, std::vector<int> E ) 
{
  n=l;
  int i, j;

  first = new int *[n];
  last  = new int *[n];    
  node  = new int**[n];

  for(i=0; i<n; ++i) {
    first[i] = new int[n];
    last[i]  = first[i];
    node [i] = new int*[n];
    j = n;
    while( j-- ) node[i][j] = NULL;
  }

  i=E.size();
  while( i ) {
    i-=2;
    add(E[i], E[i+1]);
    add(E[i+1], E[i]);
  }
}

MistralGraph::MistralGraph( const int l ) 
{
  init(l);
}

void MistralGraph::init( const int l )
{
  n = l;

  first = new int*[n];
  last  = new int*[n];
  node  = new int**[n];
  int i = n, j;
  while( i-- ) {
    first[i] = new int[n];
    last[i]  = first[i];
    node[i] = new int*[n];
    j = n;
    while( j-- ) node[i][j] = NULL;
  }
}

MistralGraph::~MistralGraph()
{
  cleanUp();
}

void MistralGraph::cleanUp()
{
  for(int i=0; i<n; ++i) {
    delete [] first[i];
    delete [] node[i];
  }
  delete [] first; 
  delete [] last; 
  delete [] node;
}

void MistralGraph::add(const int x, const int y)
{
  if( !node[x][y] ) {
    *last[x] = y;
    node[x][y] = last[x]++;
  }
}

void MistralGraph::remove(const int x, const int y)
{
  if( node[x][y] ) {
    --last[x];
    *node[x][y] = *last[x];
    node[x][*last[x]] = last[x];
    node[x][y] = NULL;
  }
}

int MistralGraph::degree(const int x) const
{
  return(last[x] - first[x]);
}

void MistralGraph::dfs( const int v ) const
{
  BitSet visited( 0, n-1, BitSet::empt );
  Vector<int> branches( n );
  int x, *varit;
  branches.push( v );

  while( branches.size ) {
    branches.pop( x );
    visited.insert( x );

    std::cout << x << std::endl;

    for(varit=first[x]; varit!=last[x]; ++varit) {
      if( !visited.member( *varit ) )
	branches.push( *varit );
    }
  }
}

/*!
 * return the maximum cardinality of a matching in a
 * bipartite graph
 */
int MistralGraph::AltBlumMehlorn( int nv )
{
//   int m=0, i=n, j, *k;
//   while( i-- )
//     m+=degree(i);

//   int LMAX = (int)(sqrt(n) * sqrt(m * log(n)) / n)+4;

//   int layer[n], isfree[n];
//   int next_free[n*(n+1)], pred_free[n*(n+1)];
//   int cardinality=0, l, vertex, path[n], length=0;

//   std::fill(isfree, isfree+n, 1);
//   std::fill(next_free, next_free+(n*(n+1)), NOVAL);
//   std::fill(pred_free, pred_free+(n*(n+1)), NOVAL);
 
//   // values are on layer 0
//   std::fill(layer+nv, layer+n, 0);
//   j = n;
//   for(i=nv; i<n; ++i) {
//     next_free[j] = i;
//     pred_free[i] = j;
//     j = i;
//   }

//   // variables are on layer 1
//   std::fill(layer, layer+nv, 1);
//   j = n;
//   for(i=0; i<nv; ++i) {
//     next_free[n+1+j] = i;
//     pred_free[n+1+i] = j;
//     j = i;
//   }

//   //g.invariant(layer);

//   l = 1;
//   while( l < n && cardinality < nv ) {

//     //std::cout << "layer " << l << std::endl;
//     //g.invariant(layer);

//     while( (vertex = next_free[l*(n+1) + n]) != NOVAL ) {

//       //g.invariant(layer);

//       path[length++] = vertex;
//       while( length ) {
// 	if( !layer[path[length-1]] && isfree[path[length-1]] ) {
// 	  // we found an augmenting path
// 	  ++cardinality;

// 	  // Both ends of the path are now matched 
// 	  isfree[path[0]] = 0;
// 	  isfree[path[length-1]] = 0;

// 	  for(j=0; j<length; j+=(length-1)) {
// 	      i = (layer[path[j]] * (n+1)) + path[j];
// 	      next_free[(layer[path[j]] * (n+1)) + pred_free[i]] = next_free[i];
// 	      pred_free[(layer[path[j]] * (n+1)) + next_free[i]] = pred_free[i];
// 	      pred_free[i] = NOVAL;
// 	      next_free[i] = NOVAL;
// 	      if( length < 2 ) break;
// 	  }

// 	  // reverse all edges in path;
// 	  i=0;
// 	  while( i < length-1 ) {
// 	    swap(path[i],path[i+1]);
// 	    ++i;
// 	  }

// 	  // clear the path
// 	  length = 0;

// 	  //g.invariant(layer);

// 	} else {
	  
// 	  j = path[length-1];
// 	  //do if( (j = next[path[length-1]][j]) == NOVAL ) break;
// 	  do {
// 	    k = node[path[length-1]][j]+1;
// 	    if( k = last[path[length-1]] ) break;
// 	    j = *k;
// 	  }
// 	  while( layer[path[length-1]] != layer[j]+1 );
	  
// 	  if( j != NOVAL ) {

// 	    // add jnxt to the path
// 	    path[length++] = j;
	    
// 	  }
// 	  else {

// 	    --length;
// 	    if( isfree[path[length]] ) {
// 	      i = (layer[path[length]] * (n+1)) + path[length];
// 	      next_free[(layer[path[length]] * (n+1)) + pred_free[i]] = next_free[i];
// 	      pred_free[(layer[path[length]] * (n+1)) + next_free[i]] = pred_free[i];
// 	      pred_free[i] = NOVAL;
// 	      next_free[i] = NOVAL;
// 	    }
// 	    layer[path[length]] += 2;

// 	    if( isfree[path[length]] ) {
// 	      i += 2*(n+1);
// 	      j = (layer[path[length]] * (n+1)) + n;
// 	      pred_free[(layer[path[length]] * (n+1)) + next_free[j]] = path[length];
// 	      next_free[i] = next_free[j];
// 	      next_free[j] = path[length];
// 	      pred_free[i] = n;
// 	    }

// 	    //g.invariant(layer);
// 	  }
// 	}
//       }
//     }
//     l += 2;
//   }

// //   g.print( std::cout );
// //   std::cout << std::endl;

// //   for(int x=0; x<n; ++x) {
// //     bool empty = true;
// //     for(int y=0; y<=n; ++y) 
// //       if( layer[y] == x ) empty = false;
// //     if( empty ) continue;

// //     std::cout << "layer" << x << ":";
// //     for(int y=0; y<=n; ++y) 
// //       if( layer[y] == x )
// // 	std::cout << std::setw(3) << y;
// //     std::cout << std::endl;
// //   }

// //   for(int x=0; x<n; ++x) {
// //     bool empty = ((j = next_free[x*(n+1)+n]) == NOVAL);
// //     if(empty) continue;
// //     std::cout << "free" << x << ":";
// //     do {
// //       std::cout << std::setw(3) << j ;
// //       std::cout.flush();
// //       assert(isfree[j]);
// //       assert(layer[j] == x);
// //       j = next_free[x*(n+1)+j];
// //     } while( j != NOVAL );
// //     std::cout << std::endl;
// //   }

// //   std::cout << cardinality << std::endl;

//   return cardinality; 
  return 0; 
}


void MistralGraph::print(std::ostream& o) const
{
  int i, *varit;
  for(i=0; i<n; ++i) {
    std::cout << i << ": ";
    for(varit=first[i]; varit!=last[i]; ++varit)
      std::cout << " ->" << std::setw(3) << (*varit);
    std::cout << std::endl;
  }
}

std::ostream& operator<< (std::ostream& os, const MistralGraph& g)
{
  g.print( os );
  return os;
}


