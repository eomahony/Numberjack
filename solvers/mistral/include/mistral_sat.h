
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

/** \file sat.h
    \brief Header for the SAT Solver.
*/


#ifndef _SAT_H
#define _SAT_H

#include <mistral_set.h>
#include <mistral_glo.h>
#include <iostream>
#include <iomanip>
#include <assert.h>

namespace Mistral {

  class CSP;

  // typedef int Literal;
  // typedef int Atom;
  // typedef Array<Literal> Clause;

  struct Statistics {
    // STATS
    double            start_time;
    long unsigned int nodes;
    long unsigned int unit_props;
    long unsigned int conflicts;
    unsigned int      literals;
    unsigned int      small;
    double            base_avg_size;
    double            learnt_avg_size;
  };

  struct Parameters {
    int               seed;
    double            time_limit;
    long unsigned int restart_limit;
    long unsigned int restart_base;
    double            restart_factor;
    double            activity_increment;
    int               policy;
    double            forgetfulness;
    unsigned int      randomization;
    bool              shuffle;
    int               verbosity;
    double            decay;
  };

  class RestartPolicy;



  class SatSolver 
  {

  public:
    /**@name Parameters*/
    //@{
    /// UNKNOWN - SAT - UNSAT - LIMITOUT
    int status;
    /// number of atoms
    int numAtoms;
    /// number of clauses
    int numClauses;

    /// List of atoms whose truth value are known, either through decisions or inference
    List assumptions;
    /// List of decision "levels", used to backtrack
    Vector<int> decisions;
    /// the truth value for each atom
    Literal *polarity;
    /// the clause that entailed this atom
    Clause **reason;
    /// the level at which this atom was entailed
    int *lvl;

    /// pointers to the scope of the clauses (base),
    Vector<Clause*> base;
    /// pointers to the scope of the clauses (learnt),
    Vector<Clause*> learnt;
    /// pointers to the scope of the clauses (original),
    Vector<Clause*> original;

    /// for each Literal, the list of clauses it is watched by
    Vector<Clause*> *isWatchedBy;
    //std::vector<Clause*> *isWatchedBy;

    /// Literal Activity 
    double *activity;
  
    /// utils
    BitSet visited;
    Vector<Literal> learnt_clause;
    int nextDeduction;

    /// search statistics
    Statistics stats;
    /// search parameters
    Parameters params;
    /// Restart
    RestartPolicy *restart_policy;
    //@}


    /**@name Constructors*/
    //@{
    SatSolver();
    SatSolver(const char* filename);
    SatSolver(CSP& model);
    virtual ~SatSolver();
    void initParams();
    void setPolicy(const int p);
    virtual void init(const int n, const int m);
    void parseDimacs(const char* filename);
    //@}

    /**@name Miscellanous*/
    //@{
    void printAll(std::ostream& o) const;
    void printWatchers(std::ostream& o, int beg=NOVAL, int end=NOVAL) const;
    void printDecisions(std::ostream& o, bool m=true) const;
    void printClauses(std::ostream& o) const;
    void printClause(std::ostream& o, Clause* cl) const;
    //@}

    /**@name Solving Methods*/
    //@{
    /// Solves the problem and set the status accordingly
    virtual int solve();
    /// Minisat style conflict directed search
    int iterative_search();
    /// Finds an explanation and uses it to backjump
    int analyze( Clause *conflict, int& lit, const bool learn );
    /// Chooses the next Literal to branch on
    Literal choice();
    /// Create a choice point and add l to the clause base
    void makeDecision(const Literal l);
    /// Backjump to the choice point where l was entailed
    void backtrackTo( const int backtrackLevel );
    /// Add l to the clause base
    void addLiteral(const Atom a, const Literal l);
    /// Returns the conflict index, or -1 if there was no conflict
    Clause* unitPropagate();
    /// Returns the conflict index, or -1 if there was no conflict
    Clause* updateWatcher(const int cw, const Atom x, const Literal p);
    //@}

    /**@name Clause Base Methods*/
    //@{
    /// Reduces every literal's activity by a constant factor
    void decayActivity();
    /// Returns NULL if the clause is always true, and a reduced version otherwise
    Clause* reduce(Clause* clause);
    /// Reduces all clauses that can be
    void simplifyDataBase();
    /// Add a clause to the base/learnt
    void addClause( Vector<Literal>& conflict );
    void addClause( Vector<Clause*>& clauseList, 
		    Vector<Literal>& conflict,
		    double& avgsize );
    /// Remove a clause from the base/learnt
    void removeClause( Vector<Clause*>& clauseList, 
		       const int cidx,
		       double& avgsize );
    /// Forget learnt clauses that do not meet a given criterion
    void forget();
    /// Add a clause to the original base
    void addOriginalClause( Vector<Literal>& conflict );
    //@}

    /**@name Utils*/
    //@{
    int checkSolution();
    //int atom( const Literal l ) const;
    void shuffle();
    bool rlimitExpired();
    bool limitExpired();
    //@}
  };




  class RestartPolicy {

  public:

    long unsigned int& restart_base;
    long unsigned int& restart_limit;
    double&            restart_factor;
  
    RestartPolicy(SatSolver*);
    virtual ~RestartPolicy();
  
    virtual void update() = 0;

  };



  class Geometric : public RestartPolicy {

  public:
  
    Geometric(SatSolver*);
    virtual ~Geometric();
  
    virtual void update();

  };


  class Luby : public RestartPolicy {

  private:

    unsigned int iteration;
    unsigned int binlog( unsigned int x );
    unsigned int luby_seq( unsigned int iteration );

  public:
  
    Luby(SatSolver*);
    virtual ~Luby();
  
    virtual void update();

  };

};


int atom( const Mistral::Literal l );
// int value( const Atom a );
// int literal( const Atom a );
// int ground( const Atom a );

int compar(const void *a, const void *b);
void initSort(double *sa);

using namespace Mistral;

inline void SatSolver::shuffle()
{
  int i, j, k;
  for(i=assumptions.size; i<numAtoms; ++i)
    {
      //j = i+rand()%(numAtoms-i);
      j = i+randint(numAtoms-i);
      k = assumptions[j];
      assumptions[j] = assumptions[i];
      assumptions[i] = k;
      assumptions.index_[k] = i;
    }
}


inline Clause* SatSolver::reduce(Clause* clause)
{
  Atom x;
  Literal* data = clause->data;
  int j, sz;

  // check the watchs first
  for(j=2; j;) {
    x = atom(data[--j]);
    if( !lvl[x] && polarity[x] == data[j] ) {
      return NULL;
    }
  }
  sz = clause->size;
  for(j=sz; j>2;) {
    x = atom(data[--j]);
    if( lvl[x] ) continue;
    if( polarity[x] == data[j] ) return NULL;
    data[j] = data[--sz];
  }

  clause->size = sz;
  return clause;
}

inline void SatSolver::simplifyDataBase()
{
  int i, oldsize;
  for(i=base.size; i;) {
    oldsize = base[i-1]->size;
    if(!reduce(base[--i])) 
      removeClause(base, i, stats.base_avg_size);
    else if(oldsize != base[i]->size) {
      oldsize -= base[i]->size;
      stats.base_avg_size -= oldsize/(double)(base.size);
    }
  }
  for(i=learnt.size; i;) {
    oldsize = learnt[i-1]->size;
    if(!reduce(learnt[--i])) 
      removeClause(learnt, i, stats.learnt_avg_size );
    else if(oldsize != learnt[i]->size) {
      oldsize -= learnt[i]->size;
      stats.learnt_avg_size -= oldsize/(double)(learnt.size);
    }
  }
}

inline void SatSolver::addClause( Vector<Literal>& conflict )
{
  if(conflict.size > 1) {
    Clause *cl = Clause::Array_new(conflict);
    base.push( cl );
    double size = base.size;
    stats.base_avg_size = (stats.base_avg_size*(size-1) + double(conflict.size))/size;
    if( conflict.size < 4 ) ++stats.small;
  } else {
    Literal p = conflict[0];
    Atom x = atom(p);
    if(!assumptions.member(x))
      addLiteral(x, p);
    else if(p != polarity[x]) 
      status = UNSAT;
  }
}

inline void SatSolver::addClause( Vector<Clause*>& clauseList, 
				  Vector<Literal>& conflict,
				  double& avgsize )
{
  if(conflict.size > 1) {
    Clause *cl = Clause::Array_new(conflict);
    clauseList.push( cl );
    isWatchedBy[conflict[0]].push(cl);
    isWatchedBy[conflict[1]].push(cl);
    
    //isWatchedBy[conflict[0]].push_back(cl);
    //isWatchedBy[conflict[1]].push_back(cl);
      
    double size = clauseList.size;
    avgsize = (avgsize*(size-1) + double(conflict.size))/size;
    if( conflict.size < 4 ) ++stats.small;

  } else {
    Literal p = conflict[0];
    Atom x = atom(p);
    if(!assumptions.member(x))
      addLiteral(x, p);
    else if(p != polarity[x]) 
      status = UNSAT;      
  }
}

inline void SatSolver::addOriginalClause( Vector<Literal>& conflict )
{
  Clause *cl = Clause::Array_new(conflict);
  original.push( cl );
}

inline void SatSolver::removeClause( Vector<Clause*>& clauseList, 
				     const int cidx,
				     double& avgsize )
{
  Clause *clause = clauseList[cidx];

  isWatchedBy[clause->data[0]].remove( clause );
  isWatchedBy[clause->data[1]].remove( clause );
  clauseList.erase( cidx );

  double size = clauseList.size;
  if(size > 0) avgsize = (avgsize*(size+1) - double(clause->size))/size;
  else         avgsize = 0;

  free(clause);
}

inline void SatSolver::forget()
{

  if( params.forgetfulness > 0.0 ) {
    int nlearnt = learnt.size;
    double sa[nlearnt];
    Clause *tmp[nlearnt];
    int i, j, order[nlearnt];
    initSort(&(sa[0]));
    for(i=0; i<nlearnt; ++i)
      {
	order[i] = i;
	sa[i] = 0.0;
	Clause& clause = *(learnt[i]);
	j=clause.size;
	while(j--)
	  //for(; j;) {
	  //--j;
	  sa[i] += activity[clause[j]];
	//}
	sa[i] /= (clause.size * clause.size);
      }
    qsort(order, nlearnt, sizeof(int), compar);
    for(i=0; i<nlearnt; ++i)
      tmp[i] = learnt[order[i]]; 
    for(i=0; i<nlearnt; ++i) {
      learnt[i] = tmp[i];
      //printClause( std::cout, learnt[i] );
      //    std::cout << " " << learnt[i]->size;
    }
    //  std::cout << std::endl;
    
    int keep = (int)((double)nlearnt * (1.0-params.forgetfulness));
    
    for(i=nlearnt; i>keep;)
      removeClause( learnt, --i, stats.learnt_avg_size );
    while(i>1) {
      --i;
      if(sa[order[i]] == sa[order[i-1]])
	removeClause( learnt, i, stats.learnt_avg_size );
    }
    
    //    for(i=0; i<learnt.size; ++i) {
    //      std::cout << " " << learnt[i]->size;
    //    }
    //    std::cout << std::endl;
    
    // //   //exit(0);
  }
}

inline bool SatSolver::rlimitExpired()
{
  return( params.restart_limit && (stats.conflicts > params.restart_limit) );
}

inline bool SatSolver::limitExpired()
{
  return( params.time_limit>0 && ((getRunTime() - stats.start_time) > params.time_limit) );
}

inline int SatSolver::iterative_search()
{
  Literal p = 0;
  Atom a;
  while(status == UNKNOWN) {
    Clause *conflict = unitPropagate();
    if(decisions.size == 0 && assumptions.size) simplifyDataBase();
    if( conflict ) {

#ifdef _DEBUGSEARCH
      for(int d=0; d<decisions.size; ++d) std::cout << " ";
      std::cout << "conflict: ";
      printClause(std::cout, conflict);
      std::cout << " " << (assumptions.size) << std::endl;
#endif

      if( decisions.size == 0 ) status = UNSAT;
      else if( limitExpired() ) return UNKNOWN;
      else if( rlimitExpired() ) {
	backtrackTo(0);
	status = LIMITOUT;
      } 
      else backtrackTo( analyze( conflict, p, true ) );
      if(status != LIMITOUT) {
	a = atom(p);

	//std::cout << decisions.size << " " << (assumptions.member(a)) << std::endl;
	if( !decisions.size && assumptions.member(a) )
	  status = UNSAT;
	else {	
	  nextDeduction = assumptions.size;
	  addLiteral(a, p);
	  //reason[a] = learnt.back();
	}
      }
    } else {
      if( (int)(assumptions.size) != numAtoms ) makeDecision(choice());	  
      else status = checkSolution();
    }
  }

  return status;
}

inline int SatSolver::analyze( Clause *conflict, int& lit, const bool learn )
{
#ifdef _DEBUGNOGOOD
  printDecisions( std::cout , 0 );
#endif

  learnt_clause.clear();

  int j, backtrackLevel = 0;
  int pathC = 0, index = assumptions.size;
  Literal p=0, q;
  Atom a;

  //if(conflict->size) {

    learnt_clause.push(p);
    
    do {

      //std::cout << conflict << std::endl;

      // add the parents of the conflict to the current set of visited atoms
      Clause& con = *conflict;

#ifdef _DEBUGNOGOOD
      printClause( std::cout, conflict );
      std::cout << std::endl;
#endif

      for(j=0; j<con.size; ++j) {
	q = con[j];
	a = atom(q);

#ifdef _DEBUGNOGOOD
	std::cout << "\t" << q << ": ";
#endif

	if( !visited.fastMember(a) ) {
	  activity[q] += params.activity_increment;
	  visited.fastInsert(a);
	  // we'll need to replace 'a' by its parents since its level is too high
	  if(lvl[a] >= decisions.size) {

#ifdef _DEBUGNOGOOD
	    std::cout << "expend" << std::endl;
#endif

	    ++pathC;
	  } else {
	    // q's level is below the current level, hence we are not expending it further

#ifdef _DEBUGNOGOOD
	    std::cout << "add to the clause" << std::endl;
#endif

	    learnt_clause.push(q);
	    if(lvl[a] > backtrackLevel)
	      backtrackLevel = lvl[a];
	  }
	}

#ifdef _DEBUGNOGOOD
	else {
	  std::cout << "visited" << std::endl;
	}
#endif

      }
      // jump to the next visited atom that need be further expended
      while(!visited.fastMember(assumptions[--index]));
      a = assumptions[index];
      p = polarity[a];

#ifdef _DEBUGNOGOOD
      std::cout << "explore " << p << " ";
      std::cout.flush();
#endif

      if( pathC > 1 ) {
	// there are still atoms to expand, we start with 'a'
	conflict = reason[a];
	visited.fastInsert(a);
      } 
#ifdef _DEBUGNOGOOD
      else {
	std::cout << std::endl;
      }
#endif

    } while( --pathC );
    // p is the last decision, since all atoms above it in the
    // assumption stack have been skipped or expended.
    learnt_clause[0] = -p;    

#ifdef _DEBUGNOGOOD
    std::cout << "add the negation of the last decision: " << -p << std::endl;
    std::cout << "(";
    for(int i=0; i<learnt_clause.size; ++i)
      std::cout << " " << learnt_clause[i];
    std::cout << " )" << std::endl;
#endif


      
    if( learn && learnt_clause.size != 1 ) {
      addClause( learnt, learnt_clause, stats.learnt_avg_size );
      reason[atom(p)] = learnt.back();
    }
    visited.clear();
    lit = -p;    

#ifdef _DEBUGNOGOOD
    std::cout << "backtraclLevel = " << backtrackLevel << "/" << (decisions.size) << std::endl;
#endif

    ++stats.conflicts;

    //std::cout << (decisions.size - backtrackLevel) << std::endl;
    //}

  return backtrackLevel;
}

inline void SatSolver::decayActivity()
{
  if(params.decay >= 0)
    for(int x=1; x<=numAtoms; ++x)
      {
	activity[ x] *= params.decay;
	activity[-x] *= params.decay;
      }
}

inline Literal SatSolver::choice()
{
  decayActivity();

  unsigned int crd = 0;
  double best[params.randomization], cur;
  int i=assumptions.size, j, k;
  Atom x[params.randomization], y;
  Literal p;
  while(i < numAtoms)
    {	       
      y = assumptions[i];

      cur = activity[y];
      cur += activity[-y]; 

      for(j=crd; j && cur>best[j-1]; --j);
      for(k=crd; k>j; --k) {
	x[k] = x[k-1];
	best[k] = best[k-1];
      }
      best[j] = cur;
      x[j] = y;

      if(crd<params.randomization) ++crd;
      ++i;
    }

  j = randint(crd);
  p = (activity[-x[j]] < activity[x[j]] ? -x[j] : x[j]);

#ifdef _DEBUGSEARCH
  for(int d=0; d<decisions.size; ++d) std::cout << " ";
  std::cout << "decide: " << p << " " << (assumptions.size) << std::endl;
#endif

  return p;    
}

inline void SatSolver::makeDecision(const Literal l)
{
  ++stats.nodes;
  nextDeduction = assumptions.size;
  decisions.push( nextDeduction );
  Atom a = atom(l);
  reason[a] = NULL;
  addLiteral(a, l);
}

inline void SatSolver::backtrackTo( const int backtrackLevel )
{
  assumptions.revertTo( decisions.popUntil(backtrackLevel) );
}

inline void SatSolver::addLiteral(const Atom a, const Literal l)
{    
  if(!decisions.size) ++stats.literals;
  assumptions.insert(a);
  polarity[a] = l;
  lvl[a] = decisions.size;
}

inline Clause* SatSolver::unitPropagate()
{
  Clause* conflict = NULL;
  int cw;
  Literal p;
  Atom a;
  while( nextDeduction < (int)(assumptions.size) && !conflict ) {
    ++stats.unit_props;
    a = assumptions[nextDeduction];
    p = -polarity[a];
    cw = isWatchedBy[p].size;
    while(cw-- && !conflict) {
      conflict = updateWatcher(cw, a, p);
    }
    ++nextDeduction;
  }
  return conflict;
}
  
inline Clause* SatSolver::updateWatcher(const int cw, const Atom x, const Literal p)
{
  Clause *cl = isWatchedBy[p][cw];
  Clause& clause = *cl;
  int j;
  Literal q, r;
  Atom y, z;

  //   std::cout << "\t";
  //   printClause( std::cout, isWatchedBy[p][cw] );
  //   //std::cout << std::endl;

  //ensure that p is the second watched lit
  if( clause[1] != p ) {
    q = clause[1];
    clause[0] = q;
    clause[1] = p;
  } else q = clause[0];
  y=atom(q);


  //   if( assumptions.member(y) && polarity[y] == q ) {
  //     std::cout << " -valid (" << y << ") " << std::endl;
  //   } 

  //check if the other watched lit is assigned
  if( !assumptions.member(y) || polarity[y] != q ) {

    //    std::cout << " search for a replacement " ;

    for(j=2; j<clause.size; ++j) {
      // for each literal q of the clause,
      r = clause[j];
      z = atom(r);
	
      if( !assumptions.member(z) ) { // this literal is not set
	// then it is a good candidate to replace p
	clause[1] = r;
	clause[j] = p;
	isWatchedBy[p].erase(cw);
	isWatchedBy[r].push(cl);

	//	std::cout << " " << r << std::endl;

	break;	
      }
	
      // if it is set true, then the clause is satisfied
      else if( polarity[z] == r ) {
	
	//	std::cout << " the clause is valid" << std::endl;
	
	break;
      }
    }
      
    if( j == clause.size ) // no replacement could be found
      { 
	if( !assumptions.member(y) ) {
	  // the last literal (other watched lit) is not set yet, we set it
	  addLiteral(y, q);
	  reason[y] = cl;

#ifdef _DEBUGNOGOOD
	  std::cout << "unit prune disjunct b" << y << " because ";
	  printClause( std::cout, cl );
	  std::cout << std::endl;
#endif

	  //	  std::cout << " prop " << q << std::endl;
	  
	} else 
	  // it is set to false already, we fail
	  if( polarity[y] == -q ) {

	    //	    std::cout << " fail" << std::endl;

	    return cl;
	  }
      }
  }

  //   std::cout << "\t";
  //   printClause( std::cout, isWatchedBy[p][cw] );
  //   std::cout << std::endl;

  return NULL;
}



#endif
