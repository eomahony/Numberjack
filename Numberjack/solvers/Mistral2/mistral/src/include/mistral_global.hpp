
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

  The author can be contacted electronically at emmanuel.hebrard@gmail.com.
*/


/** \file mistral_global.hpp
    \brief Global Definitions.
*/


#ifndef _MISTRAL_GLOBAL_HPP
#define _MISTRAL_GLOBAL_HPP

#include <stdint.h>

#include <string>
#include <iostream>

#include <mistral_structure.hpp>


// #define _DEBUG_AC true
// #define _DEBUG_PROPAGATE true
//#define _DEBUG_SEARCH true


namespace Mistral {


  typedef int* BoolDomain;
#define GET_MIN(d) (!(*d & 1))
#define GET_MAX(d)   (*d >> 1)
#define IS_GROUND(d) (*d != 3)
#define IS_FREE(d)   (*d == 3)
#define GET_VAL(d)   (*d-1)
#define GET_SIZE(d) ((*d+1)/2)



  typedef int Event;
  typedef int Outcome;
  typedef int PropagationOutcome;
  typedef int RewritingOutcome;


  //typedef unsigned int Literal;
  //typedef unsigned int Atom;
  typedef unsigned int Value;

  typedef Array<Literal> Clause;



#define _DIV_ARITY true

  // typedef int ExtLiteral;
  // typedef int LiteralType;
  // typedef int LiteralVar;

  // inline LiteralType get_type(ExtLiteral l) { return l&3; }
  // inline LiteralVar get_variable(ExtLiteral l) { return l/4; }
  // inline ExtLiteral get_literal(const LiteralType lt, const int LiteralVar lv) { return (lv*4) | lt; }


  /** event */
  




  // class Explanation {

  // public:
    
  //   typedef Literal* iterator;

  //   virtual iterator begin() = 0;
  //   virtual iterator end() = 0;
    
  // };

  // class Clause : virtual public Explanation, public Array<Literal> {
    
  // public:

  //   virtual Explanation::iterator begin() { return &(data[0]); }
  //   virtual Explanation::iterator end  () { return &(data[size]); }

  // //   inline Literal& operator [] (const int i)       { return Array<Literal>::data[i]; }
  // //   inline Literal  operator [] (const int i) const { return Array<Literal>::data[i]; }
  
  // };


  //#define _DEBUG_ACTIVITY true
  //#define _DEBUG_LESS (id == 84)
  //#define _DEBUG_PARITY ((id == 2))
  //#define _DEBUG_BUILD true
  //#define _DEBUG_BACKTRACK (id == 84)
  //#define _DEBUG_RESTORE true
  //#define _DEBUG_RELAX (id == 84)
  //#define _DEBUG_HISTORY true  
  //#define _DEBUG_GENPROPAG (id == 5)
  //#define _DEBUG_MUL (id==7660)
  //#define _DEBUG_MAX (id==288)
  //#define _DEBUG_REWRITE true
  //#define _DEBUG_DISJUNCTIVE true
  //#define _DEBUG_RDISJUNCTIVE true
  //#define _DEBUG_LEX true
  //#define _DEBUG_ADD (id == 48505)
  //#define _DEBUG_DIV true
  //#define _DEBUG_WEIGHTEDSUM (id == 0)
  //#define _DEBUG_WEIGHTEDBOOLSUM (id == 4)
  //#define _DEBUG_ELEMENT true
  //#define _DEBUG_ELEMENT ((id==0))
  //#define _DEBUG_CLIQUENOTEQUAL true
  //#define _DEBUG_MIN true
  //#define _DEBUG_MIN ((id==23))
  //#define _DEBUG_MIN ((scope[2].id()==801))
  //#define _DEBUG_MAX true
  //#define _DEBUG_MAX (get_solver()->parameters.verbosity == -7)
  //#define _DEBUG_MAX ((id==368))
  //#define _DEBUG_MEMORY true
  //#define _DEBUG_SEARCH true //(parameters.verbosity > 1)
  //#define _MONITOR true
  //#define _DEBUG_NOGOOD true
  //(statistics.num_filterings > 640)
  //((statistics.num_filterings == 48212) || (statistics.num_filterings == 46738))
  //#define _DEBUG_UNITPROP true
  //#define _DEBUG_WATCH true
  //#define _DEBUG_CHECKER (cur_iteration > 39950) 
  //#define _DEBUG_CHECKER true
  //#define _DEBUG_MOD true  
  //#define _DEBUG_INTERVALS true
  //#define _DEBUG_INTERVAL_DIV true
  //#define _DEBUG_CMOD true  
  //#define _DEBUG_ABS true  
  //#define _DEBUG_DIV true  
  //#define _DEBUG_OFFSET ((id==1309))
  //#define _DEBUG_AMSC  ((id==2912))
  //#define _DEBUG_AC true //((statistics.num_filterings >= 177))

  
#define NORESTART 0
#define GEOMETRIC 1
#define LUBY 2

#define LARGE_VALUE NOVAL/16384
#define SMALL_VALUE -NOVAL/16384
  //#define INFTY  NOVAL/2
#define MAXINT NOVAL
#define MININT -NOVAL
#define MIN_CAPACITY 16
#define NULL_ATOM 0xffffffff
  
#define BOUND_CONSISTENCY 1
#define FORWARD_CHECKING 0

#define CONST_VAR   1
#define BOOL_VAR    2
#define RANGE_VAR   4
#define BITSET_VAR  8
#define LIST_VAR    16
#define VIRTUAL_VAR 0
  //#define DYN_VAR     27
#define DYN_VAR     31
#define REMOVED_VAR 512
#define EXPRESSION  3

#define NUM_METHODS 30
#define NUM_VARTYPES 6


  // able to explain its pruning
#define EXPLAINED  0x00800000
  // binary, must be propagated on event
#define BINARY     0x80000000
  // ternary, must be propagated on event
#define TERNARY    0x40000000
  // all other constraints
#define GLOBAL     0x00000000
  // whether it should not be propagated on event
#define POSTPONED  0x20000000
  // whether it should be pushed on the constraint stack
#define PUSHED     0x10000000
  // whether it should be awakwen on its own changes
#define IDEMPOTENT 0x08000000
#define ACTIVITY   0x04000000
//   //#define CTYPE      0x07ffffff
// #define CTYPE      0x03ffffff
// #define ITYPE      0xf8000000
#define RELAXED    0x02000000
#define POSTED     0x01000000
#define CTYPE      0x007fffff
#define ITYPE      0xff800000

  static const int size_byte[8] = {0,1,1,2,1,2,2,3};

#ifdef _PROFILING

  static const int VARTYPE[17] = {0,1,2,-1,3,-1,-1,-1,4,-1,-1,-1,-1,-1,-1,-1,5};



    #define _m_get_size_          0
    #define _m_get_degree_        1
    #define _m_get_min_           2
    #define _m_get_max_           3
    #define _m_get_initial_min_   4
    #define _m_get_initial_max_   5
    #define _m_get_min_pos_       6
    #define _m_get_max_neg_       7
    #define _m_next_              8
    #define _m_prev_              9
    #define _m_is_range_          10
    #define _m_is_ground_         11
    #define _m_equal_             12
    #define _m_contain_           13
    #define _m_intersect_range_   14
    #define _m_included_range_    15
    #define _m_includes_range_    16
    #define _m_intersect_set_     17
    #define _m_included_set_      18
    #define _m_includes_set_      19
    #define _m_intersect_to_      20
    #define _m_union_to_          21
    #define _m_remove_            22
    #define _m_set_domain_value_  23
    #define _m_set_min_           24
    #define _m_set_max_           25
    #define _m_set_domain_set_    26
    #define _m_remove_set_        27
    #define _m_remove_interval_   28
    #define _m_restore_           29

  static const char* VAR_NAME[NUM_VARTYPES] = {"virtual", "constant", "boolean", "range", "bitset", "list"};
  static const char* METHOD_NAME[NUM_METHODS] = {
    "get_size"         ,
    "get_degree"       ,
    "get_min"          ,
    "get_max"          ,
    "get_initial_min"  ,
    "get_initial_max"  ,
    "get_min_pos"      ,
    "get_max_neg"      ,
    "next"             ,
    "prev"             ,
    "is_range"         ,
    "is_ground"        ,
    "equal"            ,
    "contain"          ,
    "intersect_range"  ,
    "included_range"   ,
    "includes_range"   ,
    "intersect_set"    ,
    "included_set"     ,
    "includes_set"     ,
    "intersect_to"     ,
    "union_to"         ,
    "remove"           ,
    "set_domain" ,
    "set_min"          ,
    "set_max"          ,
    "set_domain_set",
    "remove_set"       ,
    "remove_interval"  ,
    "restore"          
  };

#endif
  
#define DOMAIN_C     1
#define RANGE_C      2
#define UB_C         4
#define LB_C         8
#define VALUE_C     16

#define NO_EVENT     0
#define DOMAIN_EVENT 1
  // DOMAIN_C+RANGE_C
#define RANGE_EVENT  3
  // DOMAIN_C+RANGE_C+UB_C
#define UB_EVENT     7
  // DOMAIN_C+RANGE_C+LB_C
#define LB_EVENT     11
  // DOMAIN_C+RANGE_C+VALUE_C
#define VALUE_EVENT  31
#define FAIL_EVENT   32

#define SUPPRESSED   1
  //#define REPLACED     1

#define DOMAIN_CHANGED(e) (bool)((e)&DOMAIN_C)
#define RANGE_CHANGED(e)  (bool)((e)&RANGE_C)
#define LB_CHANGED(e)     (bool)((e)&LB_C)
#define UB_CHANGED(e)     (bool)((e)&UB_C)
#define ASSIGNED(e)       (bool)((e)&VALUE_C)
#define FAILED(e)         (bool)((e)&FAIL_EVENT)



  //   //#define __modulo_fct__(x,m)

  //int __modulo_fct__(const int x, const int m);
//   int mod = x%m;
//   if(mod && (mod<0) != (m<0))  mod += m;
//   return mod;
// }



  std::string int2str(const int x);

  std::string event2strc(Event e); //  {
  std::string event2str(Event e); //  {
  //   std::string str_evt;
  //   if(ASSIGNED(e)) {
  //     str_evt = "value";
  //     if(e != VALUE_C) {
  // 	str_evt += "("; 
  // 	if(LB_CHANGED(e)) str_evt += "+lb";
  // 	if(UB_CHANGED(e)) str_evt += "+ub";
  // 	str_evt += ")";
  //     } 
  //   } else if(LB_CHANGED(e)) {
  //     str_evt = "lb";
  //     if(UB_CHANGED(e)) str_evt += "/ub";
  //   } else if(UB_CHANGED(e)) {
  //     str_evt = "ub";
  //   } else {
  //     str_evt = "domain";
  //   }
  //   return str_evt;
  // }



  //#define IS_FAIL(e) ((e)&FAIL_EVENT)

#define SAT      1
#define OPT      3
#define UNSAT    0
#define UNKNOWN  2
#define LIMITOUT 4

#define IS_OK(o) (o<0)
#define FAILURE(x) x
#define CONSISTENT -1
  
#define _VALUE_ 0
#define _RANGE_ 1
#define _DOMAIN_ 2
#define _NEVER_ -1

#define EVENT_TYPE(e) (2-(RANGE_CHANGED(e))-(ASSIGNED(e)))

  /*
0.....1..2...........3
0.....1..2..........3.
.0....1..2..........3.
..0....1.2..........3.
..0....1..2..........3 XX
   */


  //std::string event2str(Event e);//  {

  std::string outcome2str(Outcome e);//  {
  //   std::string str_out;
  //   switch(e) {
  //   case SAT: str_out = "SAT"; break;
  //   case OPT: str_out = "OPT"; break;
  //   case UNSAT: str_out = "UNSAT"; break;
  //   case UNKNOWN: str_out = "UNKNOWN"; break;
  //   case LIMITOUT: str_out = "LIMITOUT"; break;
  //   }
  //   return str_out;
  // }


  std::string domain2str(int d);



/**
   This class stores information about the reason for pruning
   1/ in the case of Boolean variables, we keep a pointer to the constraint/clause that entailed the pruning
   2/ in the case of Range variables, we keep the last explanation for the upper bound 
      and the last explanation for the lower bound [WRONG]
   3/ otherwise [DONT KNOW]
 */
  class DomainExplanation {

  public:

    virtual void store_reason_for_change(const Event evt, Explanation *expl) = 0;
    virtual Explanation* get_explanation(const Event evt) = 0;

  };

  class BoolDomainExplanation : public DomainExplanation {

  public:

    BoolDomainExplanation() {reason = NULL;}
    virtual ~BoolDomainExplanation() {}

    Explanation *reason;
    
    virtual void store_reason_for_change(const Event evt, Explanation *expl) { reason = expl; };
    virtual Explanation* get_explanation(const Event evt) { return reason; };

  };



  //void print_trace(std::ostream& os, )

  /**********************************************
   * Timing Memory and Command line utilities 
   *********************************************/

  bool probe();
  double cpu_time(void);
#if defined(__linux__)
  int mem_read_stat(int );
#endif
  uint64_t mem_used();


  double get_run_time();
  unsigned long int get_memory();
  void get_command_line(const char**,int*,int,const char**,const char**,int,char**,int);

  template <class WORD_TYPE>
  void print_bitset(WORD_TYPE n, const int idx, std::ostream& os) {
    int offset = 8*sizeof(WORD_TYPE)*idx;
    WORD_TYPE mask=1;
    int last, cur, serie=0, k=0;
    os << "{";
    while(mask){
      if(mask & n) {
	last = (offset + k);
	os << last;
	break;
      }
      mask = mask << 1;
      ++k;
    }
    
    ++k;
    mask <<= 1;
    while(mask){
      if(mask & n) {
	cur = (k + offset);
	if(cur == last+1) ++serie;
	else {
	  if(serie > 1) os << ".." << last ;
	  else if(serie > 0) os << "," << last;
	  os << "," << cur;
	  serie = 0;
	}
	last = cur;
      }
      ++k;
      mask = mask << 1;
    }
    if(serie > 1) os << ".." << last ;
    else if(serie > 0) os << "," << last;
    
    os << "}";
  }

  /**********************************************
   * Knuth's Random number generator (code from sp-1.4)
   **********************************************/

  // class Random {

  // public:
  //   unsigned int mistral_rand_x[56];
  //   unsigned int mistral_rand_y[256];
  //   unsigned int mistral_rand_z;

  //   int mistral_rand_j;
  //   int mistral_rand_k;

  //   Random(unsigned int seed=12345) { usrand(seed); }
  //   virtual ~Random() {}

  bool random_generator_is_ready();
    void usrand (unsigned seed);
    unsigned urand0 (void);
    unsigned urand (void);
    int randint(int upto);
    double randreal();
  //  };

  /// Other utils:

  int log2_( const unsigned int v );

}

#define MAX_URAND 0xFFFFFFFFL

   
#endif // __GLOBAL_HPP


/**
   Boolean: 1 bit for the domain, 31 bits for the assignment level
   Range: 32 bits for the min, 32 bits for the max
   Word: 5 bits for the min, 5 bits for the max, 5 bits for the size, 32 bits for the domain, 17 bits for the offset
   Integer: 32 bits for the min, 32 bits for the max, 32 bits for a ptr to the rest of the domain 

   00
   01
   10
   11
*/

/**
   class Intvar {
   public:
   char type;
   void *_ptr;





   int min(, Variable *X) {
   if(type == BOOLEAN) {
   int x = (int)(*X);
   return (x&1 && x>>1 <= level); 
   } else {
   VariableFDomain *x = (VariableFDomain*)X;
   return X->min;
   }
   }
  
   }

   -- Variables are read only structs.
   -- WriteManager class takes care of the reversibility.
*/




