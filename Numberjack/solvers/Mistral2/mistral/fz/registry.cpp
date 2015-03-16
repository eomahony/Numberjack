/* -*- mode: C++; c-basic-offset: 2; indent-tabs-mode: nil -*- */
/*
 *  Main authors:
 *     Guido Tack <tack@gecode.org>
 *
 *  Contributing authors:
 *     Mikael Lagerkvist <lagerkvist@gmail.com>
 *
 *  Copyright:
 *     Guido Tack, 2007
 *     Mikael Lagerkvist, 2009
 *
 *  Last modified:
 *     $Date: 2010-07-21 11:42:47 +0200 (Wed, 21 Jul 2010) $ by $Author: tack $
 *     $Revision: 11243 $
 *
 *  This file is part of Gecode, the generic constraint
 *  development environment:
 *     http://www.gecode.org
 *
 *  Permission is hereby granted, free of charge, to any person obtaining
 *  a copy of this software and associated documentation files (the
 *  "Software"), to deal in the Software without restriction, including
 *  without limitation the rights to use, copy, modify, merge, publish,
 *  distribute, sublicense, and/or sell copies of the Software, and to
 *  permit persons to whom the Software is furnished to do so, subject to
 *  the following conditions:
 *
 *  The above copyright notice and this permission notice shall be
 *  included in all copies or substantial portions of the Software.
 *
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 *  EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 *  MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 *  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
 *  LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 *  OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 *  WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 */

#include "registry.hpp"
#include "flatzinc.hpp"

//#include <mistral_solver.hpp>
#include <mistral_variable.hpp>

#include <vector>
#include <set>

using namespace std;

namespace FlatZinc {

  void report_unsupported(const char* con) {
    std::cout << "% " << con << " is not yet supported!" << std::endl;
    exit(1);
  }


  Registry& registry(void) {
    static Registry r;
    return r;
  }

  void
  Registry::post(Solver& s, FlatZincModel &m,
                 const ConExpr& ce, AST::Node* ann) {
    std::map<std::string,poster>::iterator i = r.find(ce.id);
    if (i == r.end()) {
      throw FlatZinc::Error("Registry",
        std::string("Constraint ")+ce.id+" not found");
    }
    i->second(s, m, ce, ann);
  }

  void
  Registry::add(const std::string& id, poster p) {
    r[id] = p;
  }

  namespace {

    int ann2icl(AST::Node* ann) {
      // if (ann) {
      //   if (ann->hasAtom("val"))
      //     return ICL_VAL;
      //   if (ann->hasAtom("domain"))
      //     return ICL_DOM;
      //   if (ann->hasAtom("bounds") ||
      //       ann->hasAtom("boundsR") ||
      //       ann->hasAtom("boundsD") ||
      //       ann->hasAtom("boundsZ"))
      //     return ICL_BND;
      // }
      return 0;
    }

    inline Vector<int> arg2intargs(AST::Node* arg, int offset = 0) {
      AST::Array* a = arg->getArray();
      Vector<int> ia(a->a.size()+offset);
      for (int i=offset; i--;)
        ia[i] = 0;
      for (int i=a->a.size(); i--;)
        ia[i+offset] = a->a[i]->getInt();
      return ia;
    }

    inline Vector<int> arg2boolargs(AST::Node* arg, int offset = 0) {
      AST::Array* a = arg->getArray();
      Vector<int> ia(a->a.size()+offset);
      for (int i=offset; i--;)
        ia[i] = 0;
      for (int i=a->a.size(); i--;) {
        ia[i+offset] = a->a[i]->getBool();
      }
      return ia;
    }

    inline
    set<int> setrange(int min, int max)
    {
      set<int> rv;
      for(int i = min; i <= max; ++i)
        rv.insert(i);
      return rv;
    }

    inline set<int> arg2intset(Solver& s, AST::Node* n) {
      AST::SetLit* sl = n->getSet();
      set<int> d;
      if (sl->interval) {
        d = setrange(sl->min, sl->max);
      } else {
        for (int i=sl->s.size(); i--; )
          d.insert(sl->s[i]);
      }
      return d;
    }

    // inline set<int> arg2intvec(Solver& s, AST::Node* n) {
    //   AST::SetLit* sl = n->getSet();
    //   Vector<int> d(sl->size);
    //   for (int i=sl->s.size(); i--; )
    //     d.add(sl->s[i]);
    //   return d;
    // }

    inline Vector<int> arg2intvec(Solver& s, AST::Node* n) {
      AST::SetLit* sl = n->getSet();
      Vector<int> d;
      if (sl->interval) {
        for(int elt = sl->min; elt<=sl->max; ++elt)
          d.add(elt);
        //d = setrange(sl->min, sl->max);
      } else {
        for (unsigned int i=0; i<sl->s.size(); i++ )
          d.add(sl->s[i]);
      }
      return d;
    }

    // inline Vector<int> arg2intvec(Solver& s, AST::Node* n) {
    //   AST::SetLit* sl = n->getSet();
    //   Vector<int> d;
    //   if (sl->interval) {
    //     for(int elt = sl->min; elt<=sl->max; ++elt)
    //       d.add(elt);
    //   } else {
    //     for (unsigned int i=0; i<sl->s.size(); i++ )
    //       d.add(sl->s[i]);
    //   }
    //   return d;
    // } 

    inline vector<set<int> > arg2intsetargs(Solver& s,
                                            AST::Node* arg, int offset = 0) {
      AST::Array* a = arg->getArray();
      if (a->a.size() == 0) {
        vector<set<int> > emptyIa(0);
        return emptyIa;
      }
      vector<set<int> > ia(a->a.size()+offset);
      for (int i=a->a.size(); i--;) {
        ia[i+offset] = arg2intset(s, a->a[i]);
      }
      return ia;
    }

    inline Vector<Variable> arg2setvarargs(Solver& s,
                                           FlatZincModel& m,
                                           AST::Node* arg,
                                           int offset = 0) {
      AST::Array* a = arg->getArray();
      if (a->a.size() == 0) {
        Vector<Variable> emptyIa;
        return emptyIa;
      }
      Vector<Variable> ia(a->a.size()+offset);
      for (int i=a->a.size(); i--;) {
        if (a->a[i]->isSetVar()) {
          ia[i+offset] = m.sv[a->a[i]->getSetVar()];
        } else {
          Vector<int> e = arg2intvec(s, a->a[i]);


          //std::cout << e << std::endl;

          Variable ex = SetVariable(e,e,e.size,e.size);

          // std::cout << "TODO: set constants" << std::endl;
          // exit(1);

          //Variable constant(a->a[i]->getInt(), a->a[i]->getInt());
          ia[i+offset] = ex; //constant;
        }
      }
      return ia;
    }


    inline Vector<Variable> arg2intvarargs(Solver& s,
                                           FlatZincModel& m,
                                           AST::Node* arg,
                                           int offset = 0) {
      AST::Array* a = arg->getArray();
      if (a->a.size() == 0) {
        Vector<Variable> emptyIa;
        return emptyIa;
      }
      Vector<Variable> ia(a->a.size()+offset);
      // for (int i=offset; i--;) {
      //   if( m.constants.find(0) == m.constants.end() )
      //     m.constants[0] = s.newVariable(0, 0);
      //   ia[i] = m.constants[0];
      // }
      for (int i=a->a.size(); i--;) {
        if (a->a[i]->isIntVar()) {
          ia[i+offset] = m.iv[a->a[i]->getIntVar()];
        } else {
          Variable constant(a->a[i]->getInt(), a->a[i]->getInt());
          ia[i+offset] = constant;
        }
      }
      return ia;
    }

    inline Vector<Variable> arg2boolvarargs(Solver& s,
                                          FlatZincModel& m,
                                          AST::Node* arg,
                                          int offset = 0, int siv=-1) {
      AST::Array* a = arg->getArray();
      if (a->a.size() == 0) {
        Vector<Variable> emptyIa;
        return emptyIa;
      }
      Vector<Variable> ia(a->a.size()+offset-(siv==-1?0:1));
      for (int i=offset; i--;) {
        ia[i] = Variable(0);
      }
      for (int i=0; i<static_cast<int>(a->a.size()); i++) {
        if (i==siv)
          continue;
        if (a->a[i]->isBool()) {
          bool value = a->a[i]->getBool();
          Variable iv = Variable(value);
          ia[offset++] = iv;
        } else if (a->a[i]->isIntVar() &&
                   m.aliasBool2Int(a->a[i]->getIntVar()) != -1) {
          ia[offset++] = m.bv[m.aliasBool2Int(a->a[i]->getIntVar())];
        } else {
          ia[offset++] = m.bv[a->a[i]->getBoolVar()];
        }
      }
      return ia;
    }

    Variable getBoolVar(Solver& s,
                      FlatZincModel& m,
                      AST::Node* n) {
      Variable x0;
      if (n->isBool()) {
        Variable constvar(n->getBool());
        x0 = constvar;
      }
      else {
        x0 = m.bv[n->getBoolVar()];
      }
      return x0;
    }

    Variable getIntVar(Solver& s,
                     FlatZincModel& m,
                     AST::Node* n) {
      Variable x0;
      if (n->isIntVar()) {
        x0 = m.iv[n->getIntVar()];
      } else {
        x0 = Variable(n->getInt(), n->getInt());
      }
      return x0;
    }

    // void fz_set2mistral_set(set<int>& from, Vector<int>& to) {
    //   for(size_t i = 0; i < from.size(); ++i) {
    //     to.add(from[i]);
    //   }      
    // }
    
    Variable getSetVar(Solver& s,
                       FlatZincModel& m,
                       AST::Node* n) {
      Variable x0;
      if( n->isSetVar() ) {
        x0 = m.sv[n->getSetVar()];
      } else {
        AST::SetLit *sl = n->getSet();
        if( sl->interval ) {
          
          //std::cout << "build set from interval [" << sl->min << ".." << sl->max << "]\n";

          x0 = SetVariable( sl->min, sl->max, (sl->max-sl->min+1), (sl->max-sl->min+1) );
        } else {
          if( sl->s.empty() ) {
            x0 = SetVariable();
          } else {
            Vector<int> elts = arg2intvec(s, n);
            elts.sort();

            //std::cout << "build set from list " << elts << "\n";

            //std::cout << elts << std::endl;
            x0 = SetVariable(elts, elts, elts.size, elts.size);
          }
        }
      }

      //std::cout << x0 << std::endl;

      return x0;
    }

    bool isBoolArray(FlatZincModel& m, AST::Node* b) {
      AST::Array* a = b->getArray();
      if (a->a.size() == 0)
        return true;
      for (int i=a->a.size(); i--;) {
        if (a->a[i]->isBoolVar() || a->a[i]->isBool())
          continue;
        if ( !a->a[i]->isIntVar() )
          return false;
        if( m.aliasBool2Int(a->a[i]->getIntVar()) == -1)
          return false;
      }
      return true;
    }

    void p_alldifferent(Solver& s, FlatZincModel &m,
                        const ConExpr& ce, AST::Node* ann) {
      Vector<Variable> va = arg2intvarargs(s, m, ce[0]);
      if(va.size)
        s.add( AllDiff(va, BOUND_CONSISTENCY) );
    }

    void p_int_eq(Solver& s, FlatZincModel &m,
                  const ConExpr& ce, AST::Node* ann) {
      if (ce[0]->isIntVar()) {
        if (ce[1]->isIntVar()) {
          s.add( getIntVar(s, m, ce[0]) == getIntVar(s, m, ce[1]) );
        } else {
          s.add( getIntVar(s, m, ce[0]) == ce[1]->getInt() );
        }
      } else {
        s.add( getIntVar(s, m, ce[1]) == ce[0]->getInt() );
      }

      //cout << "post " << s.constraints.back() << endl;
    }

    void p_int_neq(Solver& s, FlatZincModel &m,
                   const ConExpr& ce, AST::Node* ann) {
      if (ce[0]->isIntVar()) {
        if (ce[1]->isIntVar()) {
          s.add( getIntVar(s, m, ce[0]) != getIntVar(s, m, ce[1]) );
        } else {
          s.add( getIntVar(s, m, ce[0]) != ce[1]->getInt() );
        }
      } else {
        s.add( getIntVar(s, m, ce[1]) != ce[0]->getInt() );
      }
    }

    // ce[0] <= ce[1] + c
    void p_int_leq_c(Solver& s, FlatZincModel &m,
                     AST::Node* ce0, AST::Node* ce1, int c,
                     AST::Node* ann) {
      if (ce0->isIntVar()) {
        if (ce1->isIntVar()) {
          s.add( Precedence( getIntVar(s, m, ce0), -c, getIntVar(s, m, ce1) ) );


#ifdef _DEBUG_FLATZINC
          Variable x = getIntVar(s, m, ce0);
          Variable y = getIntVar(s, m, ce1);          
          std::cout << x << " in " << x.get_domain() << " + " << (-c) << " <= " << y << std::endl; 
#endif

        } else {
          s.add( getIntVar(s, m, ce0) <= ce1->getInt()+c );

#ifdef _DEBUG_FLATZINC
          Variable x = getIntVar(s, m, ce0);
          std::cout << x << " in " << x.get_domain() << " <= " << (ce1->getInt()+c) << std::endl; 
#endif

        }
      } else {
        s.add( getIntVar(s, m, ce1) >= ce0->getInt()-c );

#ifdef _DEBUG_FLATZINC
         Variable x = getIntVar(s, m, ce1);
          std::cout << x << " in " << x.get_domain() << " >= " << (ce0->getInt()-c) << std::endl; 
#endif

      }
    }

    void p_int_geq(Solver& s, FlatZincModel &m,
                   const ConExpr& ce, AST::Node* ann) {
      p_int_leq_c(s, m, ce[1], ce[0], 0, ann);
    }
    void p_int_gt(Solver& s, FlatZincModel &m,
                  const ConExpr& ce, AST::Node* ann) {
      p_int_leq_c(s, m, ce[1], ce[0], -1, ann);
    }
    void p_int_leq(Solver& s, FlatZincModel &m,
                   const ConExpr& ce, AST::Node* ann) {
      p_int_leq_c(s, m, ce[0], ce[1], 0, ann);
    }
    void p_int_lt(Solver& s, FlatZincModel &m,
                  const ConExpr& ce, AST::Node* ann) {
      p_int_leq_c(s, m, ce[0], ce[1], -1, ann);
    }

    /* Comparisons */
    void p_int_eq_reif(Solver& s, FlatZincModel& m,
                       const ConExpr& ce, AST::Node* ann) {
     if (ce[0]->isIntVar()) {
        if (ce[1]->isIntVar()) {
          s.add( (getIntVar(s, m, ce[0]) == getIntVar(s, m, ce[1])) == getBoolVar(s, m, ce[2]) ); 
        } else {
          s.add( (getIntVar(s, m, ce[0]) == ce[1]->getInt()) == getBoolVar(s, m, ce[2]) ); 
        }
     } else {
       s.add( (getIntVar(s, m, ce[1]) == ce[0]->getInt()) == getBoolVar(s, m, ce[2]) ); 
     }
    }
    void p_int_ne_reif(Solver& s, FlatZincModel& m,
                       const ConExpr& ce, AST::Node* ann) {
     if (ce[0]->isIntVar()) {
        if (ce[1]->isIntVar()) {
          s.add( (getIntVar(s, m, ce[0]) != getIntVar(s, m, ce[1])) == getBoolVar(s, m, ce[2]) ); 
        } else {
          s.add( (getIntVar(s, m, ce[0]) != ce[1]->getInt()) == getBoolVar(s, m, ce[2]) ); 
        }
     } else {
       s.add( (getIntVar(s, m, ce[1]) != ce[0]->getInt()) == getBoolVar(s, m, ce[2]) ); 
     }
    }
    void p_int_ge_reif(Solver& s, FlatZincModel& m,
                       const ConExpr& ce, AST::Node* ann) {
      s.add( (getIntVar(s, m, ce[0]) >= getIntVar(s, m, ce[1])) == getBoolVar(s, m, ce[2]) ); 
    }
    void p_int_gt_reif(Solver& s, FlatZincModel& m,
                       const ConExpr& ce, AST::Node* ann) {
      s.add( (getIntVar(s, m, ce[0]) > getIntVar(s, m, ce[1])) == getBoolVar(s, m, ce[2]) ); 
    }
    void p_int_le_reif(Solver& s, FlatZincModel& m,
                       const ConExpr& ce, AST::Node* ann) {
      s.add( (getIntVar(s, m, ce[0]) <= getIntVar(s, m, ce[1])) == getBoolVar(s, m, ce[2]) ); 
    }
    void p_int_lt_reif(Solver& s, FlatZincModel& m,
                       const ConExpr& ce, AST::Node* ann) {
      s.add( (getIntVar(s, m, ce[0]) < getIntVar(s, m, ce[1])) == getBoolVar(s, m, ce[2]) ); 
    }

    /* linear (in-)equations */
    void p_int_lin(Solver& s, FlatZincModel& m,
                   event_type op, bool strict, bool reif,
                   const ConExpr& ce,
                   AST::Node* ann) {
      // vector<int> ia = arg2intargs(ce[0]);
      // if (isBoolArray(m,ce[1])) {
      //   int c = ce[2]->getInt();
      //   vector<Variable> iv = arg2boolvarargs(s, m, ce[1]);
      //   switch(op) {
      //   case LEQ_EXPR:
      //     for(size_t i = 0; i != ia.size(); ++i)
      //       ia[i] = -ia[i];
      //     c = -c;
      //     // continue on to GEQ_EXPR
      //   case GEQ_EXPR:
      //     if( strict ) {
      //       ++c;
      //     }
      //     if( !reif )
      //       post_pb(s, iv, ia, c);
      //     else {
      //       Variable b = getBoolVar(s, m, ce[3]);
      //       post_pb_iff_reif(s, iv, ia, c, b);
      //     }
      //     break;
      //   case EQ_EXPR:
      //     assert(!strict);
      //     // pseudo-boolean equality
      //     break;
      //   case NEQ_EXPR:
      //     assert(!strict);
      //     // pseudo-boolean inequality
      //     break;
      //   case NONE: assert(0); break;
      //   }
      // } else {
      //   vector<Variable> iv = arg2intvarargs(s, m, ce[1]);
      //   int c = ce[2]->getInt();
      //   switch(op) {
      //   case GEQ_EXPR:
      //     for(size_t i = 0; i != ia.size(); ++i)
      //       ia[i] = -ia[i];
      //     c = -c;
      //     // continue on to LEQ_EXPR
      //   case LEQ_EXPR:
      //     if( strict ) {
      //       --c;
      //     }
      //     if( !reif )
      //       post_lin_leq(s, iv, ia, -c);
      //     else {
      //       Variable b = getBoolVar(s, m, ce[3]);
      //       post_lin_leq_iff_reif(s, iv, ia, -c, b);
      //     }
      //     break;
      //   case EQ_EXPR:
      //     assert(!strict);
      //     if( !reif )
      //       post_lin_eq(s, iv, ia, -c);
      //     else {
      //       Variable b = getBoolVar(s, m, ce[3]);
      //       post_lin_eq_iff_reif(s, iv, ia, -c, b);
      //     }
      //     break;
      //   case NEQ_EXPR:
      //     assert(!strict);
      //     if( !reif )
      //       post_lin_neq(s, iv, ia, -c);
      //     else {
      //       Variable b = getBoolVar(s, m, ce[3]);
      //       post_lin_neq_iff_reif(s, iv, ia, -c, b);
      //     }
      //     break;
      //   case NONE: assert(0); break;
      //   }
      // }
    	  report_unsupported("p_int_lin");

    }

    void p_int_lin_eq(Solver& s, FlatZincModel& m,
                      const ConExpr& ce, AST::Node* ann) {
      Vector<int> ia = arg2intargs(ce[0]);
      Vector<Variable> iv = arg2intvarargs(s, m, ce[1]);
      if(ce[2]->isIntVar()) {
        Variable r = getIntVar(s, m, ce[2]);
        ia.add(-1);
        iv.add(r);
        s.add( Sum(iv, ia, 0, 0) );
      } else {
        int c = ce[2]->getInt();
        if(iv.size)
          s.add( Sum(iv, ia, c, c) );
        else if(c != 0)
          s.fail();
      }
    }



    void p_bool_lin_eq(Solver& s, FlatZincModel& m,
                      const ConExpr& ce, AST::Node* ann) {
      Vector<int> ia = arg2intargs(ce[0]);
      Vector<Variable> iv = arg2boolvarargs(s, m, ce[1]);

      if(ce[2]->isIntVar()) {
        Variable r = getIntVar(s, m, ce[2]);

        //std::cout << r << std::endl;

        ia.add(-1);
        iv.add(r);
        s.add( Sum(iv, ia, 0, 0) );
      } else {
        int c = ce[2]->getInt();
        if(iv.size)
          s.add( Sum(iv, ia, c, c) );
        else if(c != 0)
          s.fail();
        //s.add( Sum(iv, ia, c, c) );
      }

      //std::cout << s << std::endl;

    }
    void p_int_lin_eq_reif(Solver& s, FlatZincModel& m,
                           const ConExpr& ce, AST::Node* ann) {
      Vector<int> ia = arg2intargs(ce[0]);
      Vector<Variable> iv = arg2intvarargs(s, m, ce[1]);

      if(ce[2]->isIntVar()) {
        Variable c = getIntVar(s, m, ce[2]);
        if(iv.size)
          s.add( (Sum(iv, ia) == c) == getBoolVar(s, m, ce[3]) );
      } else {
        int c = ce[2]->getInt();
        //s.add( (Sum(iv, ia) == c) == getBoolVar(s, m, ce[3]) );
        if(iv.size)
          s.add( (Sum(iv, ia) == c) == getBoolVar(s, m, ce[3]) );
        // else if(c != 0)
        //   s.add( getBoolVar(s, m, ce[3]) == 0 );
      }
    }
    void p_int_lin_ne(Solver& s, FlatZincModel& m,
                      const ConExpr& ce, AST::Node* ann) {
      Vector<int> ia = arg2intargs(ce[0]);
      Vector<Variable> iv = arg2intvarargs(s, m, ce[1]);

      if(ce[2]->isIntVar()) {
        Variable c = getIntVar(s, m, ce[2]);
        if(iv.size)
          s.add( Sum(iv, ia) != c );
        else 
          s.add( c != 0 );
      } else {
        int c = ce[2]->getInt();
        if(iv.size)
          s.add( Sum(iv, ia) != c );
        else if(c == 0)
          s.fail();
      }
   //s.add( Sum(iv, ia) != c );
    }
    void p_int_lin_ne_reif(Solver& s, FlatZincModel& m,
                           const ConExpr& ce, AST::Node* ann) {
      Vector<int> ia = arg2intargs(ce[0]);
      Vector<Variable> iv = arg2intvarargs(s, m, ce[1]);

      if(ce[2]->isIntVar()) {
        Variable c = getIntVar(s, m, ce[2]);
        if(iv.size)
          s.add( (Sum(iv, ia) != c) == getBoolVar(s, m, ce[3]) );
      } else {
        int c = ce[2]->getInt();
        if(iv.size)
          s.add( (Sum(iv, ia) != c) == getBoolVar(s, m, ce[3]) );
      }
    }

    void p_int_lin_le(Solver& s, FlatZincModel& m,
                      const ConExpr& ce, AST::Node* ann) {
      Vector<int> ia = arg2intargs(ce[0]);
      Vector<Variable> iv = arg2intvarargs(s, m, ce[1]);

      if(ce[2]->isIntVar()) {
        Variable c = getIntVar(s, m, ce[2]);
        if(iv.size)
          s.add( Sum(iv, ia) <= c );
        else 
          s.add( c >= 0 );
      } else {
        int c = ce[2]->getInt();
        if(ia.size == iv.size 
           && iv.size == 2
           ) {
          if(ia[0] == 1 && ia[1] == -1) {
            s.add(Constraint(new ConstraintLess(iv, -c)));
          } else if(ia[0] == -1 && ia[1] == 1) {
            Variable x = iv[0];
            iv[0] = iv[1];
            iv[1] = x;
            s.add(Constraint(new ConstraintLess(iv, -c)));
          } else {
            s.add( Sum(iv, ia, -INFTY, c) );
          }
        } else {
          if(iv.size)
            s.add( Sum(iv, ia, -INFTY, c) );
          else if(c < 0)
            s.fail();
        }
      }
    }

    void p_bool_lin_le(Solver& s, FlatZincModel& m,
                      const ConExpr& ce, AST::Node* ann) {
      Vector<int> ia = arg2intargs(ce[0]);
      Vector<Variable> iv = arg2boolvarargs(s, m, ce[1]);
      if(ce[2]->isIntVar()) {
        Variable c = getIntVar(s, m, ce[2]);
        if(iv.size)
          s.add( Sum(iv, ia) <= c );
        else 
          s.add( c >= 0 );
      } else {
        int c = ce[2]->getInt();
        if(iv.size)
          s.add( Sum(iv, ia, -INFTY, c) );
        else if(c < 0)
          s.fail();
      }
    }
    void p_int_lin_le_reif(Solver& s, FlatZincModel& m,
                           const ConExpr& ce, AST::Node* ann) {
      Vector<int> ia = arg2intargs(ce[0]);
      Vector<Variable> iv = arg2intvarargs(s, m, ce[1]);
      Variable b = getBoolVar(s, m, ce[3]);

      if(ce[2]->isIntVar()) {
        Variable c = getIntVar(s, m, ce[2]);
        if(iv.size)
          s.add( (Sum(iv, ia) <= c) == b );
        else 
          s.add( (c >= 0) == b );
      } else {
        int c = ce[2]->getInt();
        // check if this is in fact a reified precedence
        if(ia.size == iv.size 
           && iv.size == 2
           ) {
          if(ia[0] == 1 && ia[1] == -1) {
            iv.add(b);
            s.add(Constraint(new PredicateLess(iv, -c)));
          } else if(ia[0] == -1 && ia[1] == 1) {
            Variable x = iv[0];
            iv[0] = iv[1];
            iv[1] = x;
            iv.add(b);
            s.add(Constraint(new PredicateLess(iv, -c)));
          } else {
            s.add( (Sum(iv, ia) <= c) == b );
          }
        } else {
          if(iv.size)
            s.add( (Sum(iv, ia) <= c) == b );
        }
      }
    }

    void p_int_lin_lt(Solver& s, FlatZincModel& m,
                      const ConExpr& ce, AST::Node* ann) {
      Vector<int> ia = arg2intargs(ce[0]);
      Vector<Variable> iv = arg2intvarargs(s, m, ce[1]);

      if(ce[2]->isIntVar()) {
        Variable c = getIntVar(s, m, ce[2]);
        if(iv.size)
          s.add( Sum(iv, ia) < c );
        else 
          s.add( c > 0 );
      } else {
        int c = ce[2]->getInt();
        if(ia.size == iv.size 
           && iv.size == 2
           ) {
          if(ia[0] == 1 && ia[1] == -1) {
            s.add(Constraint(new ConstraintLess(iv, 1-c)));
          } else if(ia[0] == -1 && ia[1] == 1) {
            Variable x = iv[0];
            iv[0] = iv[1];
            iv[1] = x;
            s.add(Constraint(new ConstraintLess(iv, 1-c)));
          } else {
            s.add( Sum(iv, ia, -INFTY, c-1) );
          }
        } else {
          if(iv.size)
            s.add( Sum(iv, ia, -INFTY, c-1) );
          else if(c <= 0)
            s.fail();
        }
      }

      // int c = ce[2]->getInt();
      // if(iv.size)
      //   s.add( Sum(iv, ia, -INFTY, c-1) );
      // else if(c <= 0)
      //   s.fail();
    }
    void p_int_lin_lt_reif(Solver& s, FlatZincModel& m,
                           const ConExpr& ce, AST::Node* ann) {
      Vector<int> ia = arg2intargs(ce[0]);
      Vector<Variable> iv = arg2intvarargs(s, m, ce[1]);
      Variable b = getBoolVar(s, m, ce[3]);

      if(ce[2]->isIntVar()) {
        Variable c = getIntVar(s, m, ce[2]);
        if(iv.size)
          s.add( (Sum(iv, ia) < c) == b );
        else 
          s.add( (c > 0) == b );
      } else {
        int c = ce[2]->getInt();
        // check if this is in fact a reified precedence
        if(ia.size == iv.size 
           && iv.size == 2
           ) {
          if(ia[0] == 1 && ia[1] == -1) {
            iv.add(b);
            s.add(Constraint(new PredicateLess(iv, 1-c)));
          } else if(ia[0] == -1 && ia[1] == 1) {
            Variable x = iv[0];
            iv[0] = iv[1];
            iv[1] = x;
            iv.add(b);
            s.add(Constraint(new PredicateLess(iv, 1-c)));
          } else {
            s.add( (Sum(iv, ia) < c) == b );
          }
        } else {
          if(iv.size)
            s.add( (Sum(iv, ia) < c) == b );
        }
      }

      // int c = ce[2]->getInt();
      // if(iv.size)
      //   s.add( (Sum(iv, ia) < c) == getBoolVar(s, m, ce[3]) );
    }
    void p_int_lin_ge(Solver& s, FlatZincModel& m,
                      const ConExpr& ce, AST::Node* ann) {
      Vector<int> ia = arg2intargs(ce[0]);
      Vector<Variable> iv = arg2intvarargs(s, m, ce[1]);

      if(ce[2]->isIntVar()) {
        Variable c = getIntVar(s, m, ce[2]);
        if(iv.size)
          s.add( Sum(iv, ia) >= c );
        else 
          s.add( c <= 0 );
      } else {
        int c = ce[2]->getInt();
        if(ia.size == iv.size 
           && iv.size == 2
           ) {
          if(ia[0] == -1 && ia[1] == 1) {
            s.add(Constraint(new ConstraintLess(iv, c)));
          } else if(ia[0] == 1 && ia[1] == -1) {
            Variable x = iv[0];
            iv[0] = iv[1];
            iv[1] = x;
            s.add(Constraint(new ConstraintLess(iv, c)));
          } else {
            s.add( Sum(iv, ia, c, INFTY) );
          }
        } else {
          if(iv.size)
            s.add( Sum(iv, ia, c, INFTY) );
          else if(c > 0)
            s.fail();
        }
      }

      // int c = ce[2]->getInt();
      // if(iv.size)
      //   s.add( Sum(iv, ia, c, INFTY) );
      // else if(c > 0)
      //   s.fail();
    }
    void p_int_lin_ge_reif(Solver& s, FlatZincModel& m,
                           const ConExpr& ce, AST::Node* ann) {
      Vector<int> ia = arg2intargs(ce[0]);
      Vector<Variable> iv = arg2intvarargs(s, m, ce[1]);
      Variable b = getBoolVar(s, m, ce[3]);

      if(ce[2]->isIntVar()) {
        Variable c = getIntVar(s, m, ce[2]);
        if(iv.size)
          s.add( (Sum(iv, ia) >= c) == b );
        else 
          s.add( (c <= 0) == b );
      } else {
        int c = ce[2]->getInt();
        // check if this is in fact a reified precedence
        if(ia.size == iv.size 
           && iv.size == 2
           ) {
          if(ia[0] == -1 && ia[1] == 1) {
            iv.add(b);
            s.add(Constraint(new PredicateLess(iv, c)));
          } else if(ia[0] == 1 && ia[1] == -1) {
            Variable x = iv[0];
            iv[0] = iv[1];
            iv[1] = x;
            iv.add(b);
            s.add(Constraint(new PredicateLess(iv, c)));
          } else {
            s.add( (Sum(iv, ia) >= c) == b );
          }
        } else {
          if(iv.size)
            s.add( (Sum(iv, ia) >= c) == b );
        }
      }

      // int c = ce[2]->getInt();
      // if(iv.size)
      //   s.add( (Sum(iv, ia) >= c) == getBoolVar(s, m, ce[3]) );
    }
    void p_int_lin_gt(Solver& s, FlatZincModel& m,
                      const ConExpr& ce, AST::Node* ann) {
      Vector<int> ia = arg2intargs(ce[0]);
      Vector<Variable> iv = arg2intvarargs(s, m, ce[1]);
      
      if(ce[2]->isIntVar()) {
        Variable c = getIntVar(s, m, ce[2]);
        if(iv.size)
          s.add( Sum(iv, ia) > c );
        else 
          s.add( c < 0 );
      } else {
        int c = ce[2]->getInt();
        if(ia.size == iv.size 
           && iv.size == 2
           ) {
          if(ia[0] == -1 && ia[1] == 1) {
            s.add(Constraint(new ConstraintLess(iv, c+1)));
          } else if(ia[0] == 1 && ia[1] == -1) {
            Variable x = iv[0];
            iv[0] = iv[1];
            iv[1] = x;
            s.add(Constraint(new ConstraintLess(iv, c+1)));
          } else {
            s.add( Sum(iv, ia, c+1, INFTY) );
          }
        } else {
          if(iv.size)
            s.add( Sum(iv, ia, c+1, INFTY) );
          else if(c >= 0)
            s.fail();
        }
      }

      // int c = ce[2]->getInt();      
      // if(iv.size)
      //   s.add( Sum(iv, ia, c+1, INFTY) );
      // else if(c >= 0)
      //   s.fail();
//s.add( Sum(iv, ia, c+1, INFTY) );
    }
    void p_int_lin_gt_reif(Solver& s, FlatZincModel& m,
                           const ConExpr& ce, AST::Node* ann) {
      Vector<int> ia = arg2intargs(ce[0]);
      Vector<Variable> iv = arg2intvarargs(s, m, ce[1]);
            Variable b = getBoolVar(s, m, ce[3]);

      if(ce[2]->isIntVar()) {
        Variable c = getIntVar(s, m, ce[2]);
        if(iv.size)
          s.add( (Sum(iv, ia) > c) == b );
        else 
          s.add( (c < 0) == b );
      } else {
        int c = ce[2]->getInt();
        // check if this is in fact a reified precedence
        if(ia.size == iv.size 
           && iv.size == 2
           ) {
          if(ia[0] == -1 && ia[1] == 1) {
            iv.add(b);
            s.add(Constraint(new PredicateLess(iv, c+1)));
          } else if(ia[0] == 1 && ia[1] == -1) {
            Variable x = iv[0];
            iv[0] = iv[1];
            iv[1] = x;
            iv.add(b);
            s.add(Constraint(new PredicateLess(iv, c+1)));
          } else {
            s.add( (Sum(iv, ia) > c) == b );
          }
        } else {
          if(iv.size)
            s.add( (Sum(iv, ia) > c) == b );
        }
      }

      // int c = ce[2]->getInt();      
      // if(iv.size)
      // s.add( (Sum(iv, ia) > c) == getBoolVar(s, m, ce[3]) );
    }

    /* arithmetic constraints */

    void p_int_plus(Solver& s, FlatZincModel& m,
                    const ConExpr& ce, AST::Node* ann) {

    	Vector<int> ia ;
    	Vector<Variable> iv ;
    	//s.add( Sum(iv, ia, 0, 0) );
   /* 	int c =0;
    	if (!ce[0]->isIntVar())
    	{
    		c-=ce[0]->getInt();
    	}
    	else
    	{
    		ia.add(1);
    		iv.add(getIntVar(s, m, ce[0]));
    	}

    	if (!ce[1]->isIntVar())
    	{
    		c-=ce[1]->getInt();
    	}
    	else
    	{
    		ia.add(1);
    		iv.add(getIntVar(s, m, ce[1]));
    	}

    	if (!ce[2]->isIntVar())
    	{
    		c+=ce[2]->getInt();
    	}
    	else
    	{
    		ia.add(-1);
    		iv.add(getIntVar(s, m, ce[2]));
    	}

    	if (ia.size)
    		s.add( Sum(iv, ia, c,c ));
    	else
    		if (c!=0)
    			s.fail();
*/

    	      if (!ce[0]->isIntVar()) {
        s.add((getIntVar(s, m, ce[1]) + ce[0]->getInt()) == getIntVar(s, m, ce[2]));
      } else if (!ce[1]->isIntVar()) {
        s.add((getIntVar(s, m, ce[0]) + ce[1]->getInt()) == getIntVar(s, m, ce[2]));
      } else if (!ce[2]->isIntVar()) {
        s.add((getIntVar(s, m, ce[0]) + getIntVar(s, m, ce[1])) == ce[2]->getInt());
      } else {
        s.add((getIntVar(s, m, ce[0]) + getIntVar(s, m, ce[1])) == getIntVar(s, m, ce[2]));
      }

    }

    void p_int_minus(Solver& s, FlatZincModel& m,
                     const ConExpr& ce, AST::Node* ann) {
      if (!ce[0]->isIntVar()) {
        s.add((getIntVar(s, m, ce[1]) - ce[0]->getInt()) == getIntVar(s, m, ce[2]));
      } else if (!ce[1]->isIntVar()) {
        s.add((getIntVar(s, m, ce[0]) - ce[1]->getInt()) == getIntVar(s, m, ce[2]));
      } else if (!ce[2]->isIntVar()) {
        s.add((getIntVar(s, m, ce[0]) - getIntVar(s, m, ce[1])) == ce[2]->getInt());
      } else {
        s.add((getIntVar(s, m, ce[0]) - getIntVar(s, m, ce[1])) == getIntVar(s, m, ce[2]));
      }
     }

    void p_int_abs(Solver& s, FlatZincModel& m,
                     const ConExpr& ce, AST::Node* ann) {

    	Variable b = getIntVar(s, m, ce[1]);
    	Variable a = getIntVar(s, m, ce[0]);

        s.add(Abs(a) == b);

    	// // s.add(((a <= 0) <= (b == -a)));  Unsat with (0,0) ??
    	// s.add(((a < 0) <= (b == -a)));
    	// s.add((a >= 0) <= (b == a));

    }


    void p_int_div(Solver& s, FlatZincModel& m,
                     const ConExpr& ce, AST::Node* ann) {

    	Variable a = getIntVar(s, m, ce[0]);
    	Variable b = getIntVar(s, m, ce[1]);
    	Variable d = getIntVar(s, m, ce[2]);

        if(b.is_ground()) {
          s.add(a / b.get_value() == d);
        } else {
          s.add(a / b == d);
        }
        
        /*
    	int bnd = b.get_max();
    	if (b.get_min() > bnd)
    		bnd = b.get_min();
    	if (-b.get_min() > bnd)
    		bnd = (-b.get_min());
    	if ((-b.get_max())> bnd)
    		bnd = (-b.get_max());

    	Variable r = Variable (-bnd, bnd);

    	Variable abs_b = Variable (-bnd, bnd);
    	Variable abs_r = Variable (-bnd, bnd);

        	s.add(((b < 0) <= (abs_b == -b)));
        	s.add((b >= 0) <= (abs_b == b));

        	s.add(((a*r) >= 0));
        	s.add(r== a - (b*d));

        	s.add(((r < 0) <= (abs_r == -r)));
        	s.add((r >= 0) <= (abs_r == r));

        	s.add(abs_r < abs_b);
        */
    }


    void p_int_mod(Solver& s, FlatZincModel& m,
                     const ConExpr& ce, AST::Node* ann) {

    	Variable a = getIntVar(s, m, ce[0]);
    	Variable b = getIntVar(s, m, ce[1]);
    	Variable r = getIntVar(s, m, ce[2]);

        if(b.is_ground()) {
          s.add(a % b.get_value() == r);
        } else {
          s.add(a % b == r);
        }
        

        // //Decomposition

    	// int bnd_a = a.get_max();
    	// if (a.get_min() > bnd_a)
    	// 	bnd_a = a.get_min();
    	// if (-a.get_min() > bnd_a)
    	// 	bnd_a = (-a.get_min());
    	// if ((-a.get_max())> bnd_a)
    	// 	bnd_a = (-a.get_max());


    	// int bnd_b = b.get_max();
    	// if (b.get_min() > bnd_b)
    	// 	bnd_b = b.get_min();
    	// if (-b.get_min() > bnd_b)
    	// 	bnd_b = (-b.get_min());
    	// if ((-b.get_max())> bnd_b)
    	// 	bnd_b = (-b.get_max());
    	// //cout << "bnd b" << bnd_b << "\n" ;

    	// Variable d = Variable (-bnd_a, bnd_a);

    	// Variable abs_b = Variable (-bnd_b, bnd_b);
    	// Variable abs_r = Variable (-bnd_b, bnd_b);

        // 	s.add(((b < 0) <= (abs_b == -b)));
        // 	s.add((b >= 0) <= (abs_b == b));

        // 	s.add(((a*r) >= 0));
        // 	s.add(r== a - (b*d));

        // 	s.add(((r < 0) <= (abs_r == -r)));
        // 	s.add((r >= 0) <= (abs_r == r));

        // 	s.add(abs_r < abs_b);

    }

    void p_int_times(Solver& s, FlatZincModel& m,
                     const ConExpr& ce, AST::Node* ann) {
      // cout << "CAN'T POST A MUL CONSTRAINT" << endl;
      // exit(1);

    	if (ce[0]->isIntVar() && ce[1]->isIntVar())
    	{
    		if (!ce[2]->isIntVar()) {
    			s.add((getIntVar(s, m, ce[0]) * getIntVar(s, m, ce[1])) == ce[2]->getInt());
    		} else {
    			s.add((getIntVar(s, m, ce[0]) * getIntVar(s, m, ce[1])) == getIntVar(s, m, ce[2]));
    		}
    	}
    	else
    	{
    		Vector<int> ia ;
    		Vector<Variable> iv ;

    		int c = 0;
    		if ((!ce[0]->isIntVar()) && (!ce[1]->isIntVar()))
    			c-= (ce[0]->getInt() * ce[1]->getInt());
    		else
    		{
    			if (!ce[0]->isIntVar())
    			{
    				ia.add(ce[0]->getInt());
    				iv.add(getIntVar(s, m, ce[1]));
    			}
    			else
    			{
    				ia.add(ce[1]->getInt());
    				iv.add(getIntVar(s, m, ce[0]));
    			}
    		}

    		if (!ce[2]->isIntVar())
    		{
    			c+=ce[2]->getInt();
    		}
    		else
    		{
    			ia.add(-1);
    			iv.add(getIntVar(s, m, ce[2]));
    		}

    		if (ia.size)
    			s.add( Sum(iv, ia, c,c ));
    		else  if(c != 0)
    			s.fail();
    	}


/*
      Variable x;

      if (!ce[0]->isIntVar()) {
        if(ce[0]->getInt()) {
          x = (getIntVar(s, m, ce[1]) * ce[0]->getInt()) == getIntVar(s, m, ce[2]);

          //std::cout << "add " << x << std::endl; 

          s.add(x);
        } else {
          x = (getIntVar(s, m, ce[2]) == 0);
          
          //std::cout << "add " << x << std::endl; 

          s.add(x);
        }
      } else if (!ce[1]->isIntVar()) {
        if(ce[1]->getInt()) {
          x = (getIntVar(s, m, ce[0]) * ce[1]->getInt()) == getIntVar(s, m, ce[2]);

          //std::cout << "add " << x << std::endl; 

          s.add(x);
        } else {
          x = (getIntVar(s, m, ce[2]) == 0);

          //std::cout << "add " << x << std::endl; 

          s.add(x);
        }
      } else if (!ce[2]->isIntVar()) {
        s.add((getIntVar(s, m, ce[0]) * getIntVar(s, m, ce[1])) == ce[2]->getInt());
      } else {
        s.add((getIntVar(s, m, ce[0]) * getIntVar(s, m, ce[1])) == getIntVar(s, m, ce[2]));
      }
*/
      //std::cout << s << std::endl;

      // Variable x0 = getIntVar(s, m, ce[0]);
      // Variable x1 = getIntVar(s, m, ce[1]);
      // Variable x2 = getIntVar(s, m, ce[2]);
      // s.add( x0*x1 == x2 );

      // //post_mult(s, x2, x0, x1); // note the order
    }

    void p_int_negate(Solver& s, FlatZincModel& m,
                      const ConExpr& ce, AST::Node* ann) {

  	  report_unsupported("p_int_negate");

      // if( !ce[0]->isIntVar() ) {
      //   if( !ce[1]->isIntVar() ) {
      //     if( ce[0]->getInt() != - ce[1]->getInt() )
      //       throw unsat();
      //     return;
      //   }
      //   Variable x1 = getIntVar(s, m, ce[1]);
      //   x1.assign(s, -ce[0]->getInt(), NO_REASON);
      // } else if( !ce[1]->isIntVar() ) {
      //   Variable x0 = getIntVar(s, m, ce[1]);
      //   x0.assign(s, -ce[1]->getInt(), NO_REASON);
      // } else {
      //   Variable x0 = getIntVar(s, m, ce[0]);
      //   Variable x1 = getIntVar(s, m, ce[1]);
      //   post_neg(s, x0, x1, 0);
      // }
    }

    void p_int_min(Solver& s, FlatZincModel& m,
                   const ConExpr& ce, AST::Node* ann) {
      // cout << "CAN'T POST A MIN CONSTRAINT" << endl;
      // exit(1);

      Variable x0 = getIntVar(s, m, ce[0]);
      Variable x1 = getIntVar(s, m, ce[1]);
      Variable x2 = getIntVar(s, m, ce[2]);

      s.add(Min(x0, x1) == x2);

      // post_min(s, x2, x0, x1);
    }
    void p_int_max(Solver& s, FlatZincModel& m,
                   const ConExpr& ce, AST::Node* ann) {
      // cout << "CAN'T POST A MAX CONSTRAINT" << endl;
      // exit(1);

      Variable x0 = getIntVar(s, m, ce[0]);
      Variable x1 = getIntVar(s, m, ce[1]);
      Variable x2 = getIntVar(s, m, ce[2]);

      s.add(Max(x0, x1) == x2);

      // post_max(s, x2, x0, x1);
    }

    /* element constraints */
    void p_array_int_element(Solver& s, FlatZincModel& m,
                             const ConExpr& ce, AST::Node* ann) {
      Variable selector = getIntVar(s, m, ce[0]);
      Variable result = getIntVar(s, m, ce[2]);
      Vector<Variable> iv = arg2intvarargs(s, m, ce[1]);

      s.add(selector > 0);
      s.add(selector <= iv.size);
      s.add(Element(iv, selector, 1) == result);
    }
    void p_array_set_element(Solver& s, FlatZincModel& m,
                             const ConExpr& ce, AST::Node* ann) {
      Variable selector = getIntVar(s, m, ce[0]);
      Variable result = getSetVar(s, m, ce[2]);
      Vector<Variable> sv = arg2setvarargs(s, m, ce[1]);

      report_unsupported("p_array_set_element");

      /*
      std::cout << selector.get_domain() << std::endl;
      std::cout << result << std::endl;
      std::cout << sv << std::endl;
      */

      s.add(selector > 0);
      s.add(selector <= sv.size);
      //s.add(ElementSet(sv, selector, 1) == result);
    }
    void p_array_var_bool_element(Solver& s, FlatZincModel& m,
                                  const ConExpr& ce, AST::Node* ann) {
      Variable selector = getIntVar(s, m, ce[0]);
      Variable result = getBoolVar(s, m, ce[2]);
      Vector<Variable> iv = arg2boolvarargs(s, m, ce[1]);

      s.add(selector > 0);
      s.add(selector <= iv.size);
      s.add(Element(iv, selector, 1) == result);
    }
    void p_array_bool_element(Solver& s, FlatZincModel& m,
                              const ConExpr& ce, AST::Node* ann) {
      Variable selector = getIntVar(s, m, ce[0]);
      Variable result = getBoolVar(s, m, ce[2]);
      Vector<int> iv = arg2boolargs(ce[1]);
      Vector<int> good_values;

      // std::cout << "possible values: " << selector.get_domain() << std::endl;
      // std::cout << "boolargs: " << iv << std::endl;
      
      for(unsigned int i=0; i<iv.size; ++i) {
        if(iv[i]) good_values.add(i+1);
      }

      //std::cout << "good values: " << good_values << std::endl;

      s.add(selector > 0);
      s.add(selector <= iv.size);
      s.add(Member(selector, good_values) == result);

      //exit(1);
      //s.add(Element(iv, selector, 1) == result);
    }

    /* alldiff */
    void p_all_different(Solver& s, FlatZincModel& m,
                         const ConExpr& ce, AST::Node* ann) {
      Vector< Variable > iv = arg2intvarargs(s, m, ce[0]);

      if(iv.size)
        s.add( AllDiff(iv) );
    }

    /*
    %-----------------------------------------------------------------------------%
    % Requires at least 'n' variables in 'x' to take the value 'v'.
    %-----------------------------------------------------------------------------%

    predicate at_least_int(int: n, array[int] of var int: x, int: v) =
        sum(i in index_set(x)) ( bool2int(x[i] == v) ) >= n;
     */
    void p_at_least_int(Solver& s, FlatZincModel& m,
    		const ConExpr& ce, AST::Node* ann) {

    	int n= ce[0]->getInt(), v = ce[2]->getInt();
    	Vector< Variable > iv = arg2intvarargs(s, m, ce[1]);

    	int size = iv.size;

    	if (size){

    		VarArray subsequence;
    		for (int i=0; i < size; ++i)
    		{
    			subsequence.add(iv[i]==v);
    			//s.add( Free(subsequence.back()));
    		}

    		s.add( BoolSum(subsequence) >= n);
    	}
    }


    /*    	%-----------------------------------------------------------------------------%
        	% Requires at most 'n' variables in 'x' to take the value 'v'.
        	%-----------------------------------------------------------------------------%

        	predicate at_most_int(int: n, array[int] of var int: x, int: v) %=
        	%    sum(i in index_set(x)) ( bool2int(x[i] == v) ) <= n;
     */

    void p_at_most_int(Solver& s, FlatZincModel& m,
    		const ConExpr& ce, AST::Node* ann) {

    	int n= ce[0]->getInt(), v = ce[2]->getInt();
    	Vector< Variable > iv = arg2intvarargs(s, m, ce[1]);

    	int size = iv.size;

    	if (size){

    		VarArray subsequence;
    		for (int i=0; i < size; ++i)
    		{
    			subsequence.add(iv[i]==v);
    			//s.add( Free(subsequence.back()));
    		}

    		s.add( BoolSum(subsequence) <= n);
    	}
    }


    /*
     *
     %-----------------------------------------------------------------------------%
    % Constrains 'c' to be the number of occurrences of 'y' in 'x'.
    %-----------------------------------------------------------------------------%

    predicate count_eq(array[int] of var int: x, var int: y, var int: c) =
        c = sum(i in index_set(x)) ( bool2int(x[i] == y) );

    %-----------------------------------------------------------------------------%
    %-----------------------------------------------------------------------------%
     */
    void p_count_eq(Solver& s, FlatZincModel& m,
    		const ConExpr& ce, AST::Node* ann) {

    	//std::cout << "count_eq" << std::endl;
    	Vector< Variable > x = arg2intvarargs(s, m, ce[0]);
    	int size = x.size;
    	if (size){
    		VarArray subsequence;
    		if (ce[1]->isInt() ){
    			int y = ce[1]->getInt();
    			for (int i=0; i < size; ++i)
    			{
    				subsequence.add(x[i]==y);
    				//s.add( Free(subsequence.back()));
    			}
    		}
    		else
    		{
    			//Variable y = ce[1]->getIntVar();
    			Variable y = getIntVar(s,m,ce[1]);

    			for (int i=0; i < size; ++i)
    			{
    				subsequence.add(y==x[i]);
    				//s.add( Free(subsequence.back()));
    			}
    		}


    		if (ce[2]->isInt() )
    			//s.add( BoolSum(subsequence, ce[2]->getInt(), ce[2]->getInt()));
    			s.add( BoolSum(subsequence) ==  ce[2]->getInt());
    		else
    			//s.add( BoolSum(subsequence) == ce[2]->getIntVar());
    			s.add( BoolSum(subsequence) == getIntVar(s,m,ce[2]));
    	}
    }


    //Reified count_eq
    //predicate count_eq_reif(array[int] of var int: x, var int: y, var int: c,  var bool: b);
    void p_count_eq_reif(Solver& s, FlatZincModel& m,
    		const ConExpr& ce, AST::Node* ann) {

    	//std::cout << "p_count_eq_reif" << std::endl;

    	Vector< Variable > x = arg2intvarargs(s, m, ce[0]);
    	int size = x.size;
    	if (size){
    		VarArray subsequence;
    		if (ce[1]->isInt() ){
    			int y = ce[1]->getInt();
    			for (int i=0; i < size; ++i)
    			{
    				subsequence.add(x[i]==y);
    				//s.add( Free(subsequence.back()));
    			}
    		}
    		else
    		{
    			//Variable y = ce[1]->getIntVar();
    			Variable y = getIntVar(s,m,ce[1]);

    			for (int i=0; i < size; ++i)
    			{
    				subsequence.add(y==x[i]);
    				//s.add( Free(subsequence.back()));
    			}
    		}

    		if (ce[3]->isBool() ){


    			int is_true = ce[3]->getBool();
    			if (is_true){
    				if (ce[2]->isInt() )
    					//s.add( BoolSum(subsequence, ce[2]->getInt(), ce[2]->getInt()));
    					s.add( BoolSum(subsequence) ==  ce[2]->getInt());
    				else
    					//s.add( BoolSum(subsequence) == ce[2]->getIntVar());
    					s.add( BoolSum(subsequence) == getIntVar(s,m,ce[2]));
    			}
    			else {

    				if (ce[2]->isInt() )
    					//s.add( BoolSum(subsequence, ce[2]->getInt(), ce[2]->getInt()));
    					s.add( BoolSum(subsequence) !=  ce[2]->getInt());
    				else
    					//s.add( BoolSum(subsequence) == ce[2]->getIntVar());
    					s.add( BoolSum(subsequence) != getIntVar(s,m,ce[2]));
    			}
    		}
    		else {
    			if (ce[2]->isInt() )
    				//s.add( BoolSum(subsequence, ce[2]->getInt(), ce[2]->getInt()));
    				s.add( (BoolSum(subsequence) ==  ce[2]->getInt() ) == getBoolVar(s,m,ce[3]) );
    			else
    				//s.add( BoolSum(subsequence) == ce[2]->getIntVar());
    				s.add( (BoolSum(subsequence) == getIntVar(s,m,ce[2])) == getBoolVar(s,m,ce[3]));
    		}

    	}
    }



    /*
    %-----------------------------------------------------------------------------%
    % Constrains 'c' to be greater than or equal to the number of occurrences of
    % 'y' in 'x'.
    %-----------------------------------------------------------------------------%

    predicate count_geq(array[int] of var int: x, var int: y, var int: c) =
        c >= sum(i in index_set(x)) ( bool2int(x[i] == y) );

    %-----------------------------------------------------------------------------%
    %-----------------------------------------------------------------------------%
     */
    void p_count_geq(Solver& s, FlatZincModel& m,
    		const ConExpr& ce, AST::Node* ann) {

    	//std::cout << "count_eq" << std::endl;
    	Vector< Variable > x = arg2intvarargs(s, m, ce[0]);
    	int size = x.size;
    	if (size){
    		VarArray subsequence;
    		if (ce[1]->isInt() ){
    			int y = ce[1]->getInt();
    			for (int i=0; i < size; ++i)
    			{
    				subsequence.add(x[i]==y);
    				//s.add( Free(subsequence.back()));
    			}
    		}
    		else
    		{
    			Variable y = getIntVar(s,m,ce[1]); //ce[1]->getIntVar();
    			for (int i=0; i < size; ++i)
    			{
    				subsequence.add(y==x[i]);
    				//s.add( Free(subsequence.back()));
    			}
    		}

    		if (ce[2]->isInt() )
    			s.add( BoolSum(subsequence) <= ce[2]->getInt());
    		else
    			s.add( BoolSum(subsequence) <= getIntVar(s,m,ce[2])); //ce[2]->getIntVar());
    	}
    }

    /*
    %-----------------------------------------------------------------------------%
    % Constrains 'c' to be strictly greater than the number of occurrences of 'y'
    % in 'x'.
    %-----------------------------------------------------------------------------%

    predicate count_gt(array[int] of var int: x, var int: y, var int: c) =
        c > sum(i in index_set(x)) ( bool2int(x[i] == y) );

    %-----------------------------------------------------------------------------%
    %-----------------------------------------------------------------------------%
     */
    void p_count_gt(Solver& s, FlatZincModel& m,
    		const ConExpr& ce, AST::Node* ann) {

    	//std::cout << "count_eq" << std::endl;
    	Vector< Variable > x = arg2intvarargs(s, m, ce[0]);
    	int size = x.size;
    	if (size){
    		VarArray subsequence;
    		if (ce[1]->isInt() ){
    			int y = ce[1]->getInt();
    			for (int i=0; i < size; ++i)
    			{
    				subsequence.add(x[i]==y);
    				//s.add( Free(subsequence.back()));
    			}
    		}
    		else
    		{
    			Variable y =getIntVar(s,m,ce[1]);  //ce[1]->getIntVar();
    			for (int i=0; i < size; ++i)
    			{
    				subsequence.add(y==x[i]);
    				//s.add( Free(subsequence.back()));
    			}
    		}



    		if (ce[2]->isInt() )
    			s.add( BoolSum(subsequence) < ce[2]->getInt());
    		else
    			s.add( BoolSum(subsequence) < getIntVar(s,m,ce[2])); //ce[2]->getIntVar());
    	}
    }

    /*
    %-----------------------------------------------------------------------------%
    % Constrains 'c' to be less than or equal to the number of occurrences of
    % 'y' in 'x'.
    %-----------------------------------------------------------------------------%

    predicate count_leq(array[int] of var int: x, var int: y, var int: c) =
        c <= sum(i in index_set(x)) ( bool2int(x[i] == y) );

    %-----------------------------------------------------------------------------%
    %-----------------------------------------------------------------------------%
     */

    void p_count_leq(Solver& s, FlatZincModel& m,
    		const ConExpr& ce, AST::Node* ann) {

    	//std::cout << "count_eq" << std::endl;
    	Vector< Variable > x = arg2intvarargs(s, m, ce[0]);
    	int size = x.size;
    	if (size){
    		VarArray subsequence;
    		if (ce[1]->isInt() ){
    			int y = ce[1]->getInt();
    			for (int i=0; i < size; ++i)
    			{
    				subsequence.add(x[i]==y);
    				s.add(subsequence.back());
    			}
    		}
    		else
    		{
    			Variable y =getIntVar(s,m,ce[1]); // ce[1]->getIntVar();
    			for (int i=0; i < size; ++i)
    			{
    				subsequence.add(y==x[i]);
    				//s.add( Free(subsequence.back()));
    			}
    		}



    		if (ce[2]->isInt() )
    			s.add( BoolSum(subsequence) >= ce[2]->getInt());
    		else
    			s.add( BoolSum(subsequence) >= getIntVar(s,m,ce[2])); //ce[2]->getIntVar());
    	}
    }

    /*
    %-----------------------------------------------------------------------------%
    % Constrains 'c' to be strictly less than the number of occurrences of 'y'
    % in 'x'.
    %-----------------------------------------------------------------------------%

    predicate count_lt(array[int] of var int: x, var int: y, var int: c) ;% =
    %    c < sum(i in index_set(x)) ( bool2int(x[i] == y) );

    %-----------------------------------------------------------------------------%
    %-----------------------------------------------------------------------------%
     */

    void p_count_lt(Solver& s, FlatZincModel& m,
    		const ConExpr& ce, AST::Node* ann) {

    	//std::cout << "count_eq" << std::endl;
    	Vector< Variable > x = arg2intvarargs(s, m, ce[0]);
    	int size = x.size;
    	if (size){
    		VarArray subsequence;
    		if (ce[1]->isInt() ){
    			int y = ce[1]->getInt();
    			for (int i=0; i < size; ++i)
    			{
    				subsequence.add(x[i]==y);
    				//s.add(Free(subsequence.back()));
    			}
    		}
    		else
    		{
    			Variable y = getIntVar(s,m,ce[1]); //ce[1]->getIntVar();
    			for (int i=0; i < size; ++i)
    			{
    				subsequence.add(y==x[i]);
    				//s.add(Free(subsequence.back()));
    			}
    		}



    		if (ce[2]->isInt() )
    			s.add( BoolSum(subsequence) > ce[2]->getInt());
    		else
    			s.add( BoolSum(subsequence) > getIntVar(s,m,ce[2]));//ce[2]->getIntVar());
    	}
    }

    /*
    %-----------------------------------------------------------------------------%
    % Constrains 'c' to not be the number of occurrences of 'y' in 'x'.
    %-----------------------------------------------------------------------------%

    predicate count_neq(array[int] of var int: x, var int: y, var int: c) ;%=
    %    c != sum(i in index_set(x)) ( bool2int(x[i] == y) );

    %-----------------------------------------------------------------------------%
    %-----------------------------------------------------------------------------%
     */

    void p_count_neq(Solver& s, FlatZincModel& m,
    		const ConExpr& ce, AST::Node* ann) {

    	//std::cout << "count_eq" << std::endl;
    	Vector< Variable > x = arg2intvarargs(s, m, ce[0]);
    	int size = x.size;
    	if (size){
    		VarArray subsequence;
    		if (ce[1]->isInt() ){
    			int y = ce[1]->getInt();
    			for (int i=0; i < size; ++i)
    			{
    				subsequence.add(x[i]==y);
    				//s.add( Free(subsequence.back()));
    			}
    		}
    		else
    		{
    			Variable y = getIntVar(s,m,ce[1]); //ce[1]->getIntVar();
    			for (int i=0; i < size; ++i)
    			{
    				subsequence.add(y==x[i]);
    				//s.add( Free(subsequence.back()));
    			}
    		}



    		if (ce[2]->isInt() )
    			s.add( BoolSum(subsequence) != ce[2]->getInt());
    		else
    			s.add( BoolSum(subsequence) != getIntVar(s,m,ce[2])); //ce[2]->getIntVar());
    	}
    }




    /* All variablmes should be equal!
     * predicate all_equal_int(array[int] of var int: x) =
    	forall(i, j in index_set(x) where i < j) ( x[i] = x[j] );
     */

    void p_all_equal_int(Solver& s, FlatZincModel& m,
    		const ConExpr& ce, AST::Node* ann) {

    	Vector< Variable > iv = arg2intvarargs(s, m, ce[0]);

    	int size = iv.size;

    	if (size){
    		for (int i=1; i < size; ++i)
    		{
    			s.add(iv[i]==iv[0]);
    		}
    	}
    }


    /* cumulative */
    //#define _DEBUG_CUMULATIVE
    
    Vector<Variable> demand;

    int largest_demand(const void *x_, const void *y_) {
      int x = *((int*)x_);
      int y = *((int*)y_);
      if(demand[x].get_min() > demand[y].get_min())
        return 1;
      else if(demand[x].get_min() < demand[y].get_min())
        return -1;
      return 0;
    }



    void p_disjunctive(Solver& s, 
                       Vector<Variable>& start, 
                       Vector<Variable>& dur) {
      
      int n = start.size;
      for(int i=0; i<n-1; ++i) {
    	  for(int j=i+1; j<n; ++j) {
    		  //Check if it's already a Precedence!
    		  if (start[i].get_min()+dur[i].get_min() >start[j].get_min()){
    			  s.add(Precedence(start[j], dur[j].get_min() , start[i]));
    		  }
    		  else
    			  if (start[j].get_min()+dur[j].get_min() >start[i].get_min()){
    				  s.add(Precedence(start[i], dur[i].get_min() , start[j]));
    			  }
    			  else{
    				  s.add(Free(ReifiedDisjunctive(start[i], start[j], dur[i].get_min(), dur[j].get_min())));
    			  }
    	  }
      }
    }
    
    
    void p_cumulative_flow(Solver& s, 
                           Vector<Variable>& start, 
                           Vector<Variable>& dur,
                           Vector<Variable> req,
                           Variable cap, 
                           int horizon) {

      int n = start.size;
      Variable **flow = new Variable*[n+1];
      for(int i=0; i<=n; ++i) {
        flow[i] = new Variable[n+1];
      }

      for(int i=0; i<n; ++i) {
        flow[i][n] = Variable(0, req[i].get_max());
        flow[n][i] = Variable(0, req[i].get_max());
        for(int j=i+1; j<n; ++j) {
          // for each pair of tasks, create a variable for the amount flow going from ti to tj
          int mf = (req[i].get_max()>req[j].get_max() ? req[j].get_max() : req[i].get_max());
          flow[i][j] = Variable(0, mf);
          flow[j][i] = Variable(0, mf);
          // the flow can go only in one direction 
          s.add((flow[i][j]<=0) || (flow[j][i]<=0));
        }
      }
      
      // flow conservation constraints
      Vector<Variable> outflow;
      Vector<Variable> inflow;
      for(int i=0; i<n; ++i) {
        for(int j=0; j<=n; ++j) {
          if(i != j) {
            outflow.add(flow[i][j]);
            inflow.add(flow[j][i]);
          }
        }
        if(req[i].is_ground()) {
          s.add(Sum(outflow, req[i].get_min(), req[i].get_max()));
          s.add(Sum(inflow, req[i].get_min(), req[i].get_max()));
        } else {
          s.add(Sum(outflow) == req[i]);
          s.add(Sum(inflow) == req[i]);
        }
        outflow.clear();
        inflow.clear();
      }
      for(int j=0; j<n; ++j) {
        outflow.add(flow[n][j]);
        inflow.add(flow[j][n]);
      }

      if(cap.is_ground()) {
        s.add(Sum(outflow, cap.get_min(), cap.get_max()));
        s.add(Sum(inflow, cap.get_min(), cap.get_max()));
      } else {
        s.add(Sum(outflow) == cap);
        s.add(Sum(inflow) == cap);
      }
    
      // There can be a flow going from i to j only if i precedes j
      for(int i=0; i<n; ++i) {
        for(int j=0; j<n; ++j) {
          if(i != j) {
            if(dur[i].is_ground())
              s.add( (flow[i][j]>0) <= (start[i]+dur[i].get_min() <= start[j]) );
            else
              s.add( (flow[i][j]>0) <= ((start[i]+dur[i]) <= start[j]) );
          }
        }
      }

      for(int i=0; i<=n; ++i) {
        delete [] flow[i];
      }
      delete [] flow;
    }

    void p_cumulative_discretization(Solver& s, 
                                     Vector<Variable>& start, 
                                     Vector<Variable>& dur,
                                     Vector<Variable> req,
                                     Variable cap, 
                                     int horizon) {

      Variable     in_process;
      Vector<Variable>  scope;
      //Vector<Variable> cstart;
      Vector<int>      demand;
      int        process_time;
      int               total;
      int             min_req;
      int             Boolean;
      // post a sum constraint for each time point
      for(int t=0; t<=horizon; ++t) {
        //std::cout << " - time=" << t << ":" ; //<< std::endl;

        demand.clear();
        scope.clear();
        total = 0;
        Boolean = true;
        min_req = Mistral::INFTY;
        for(int i=0; i<start.size; ++i) {
          if(start[i].get_min() <= t && start[i].get_max()+dur[i].get_max() >= t) {
            total += req[i].get_max();
            if(min_req > req[i].get_min()) {
              min_req = req[i].get_min();
            }

            // this tasks can be in process at time t
            if(dur[i].is_ground()) {
              process_time = dur[i].get_min();
              // constant duration
              if(start[i].get_max()>t) {
                if(start[i].get_min()+process_time<t) {

                  if(process_time==1) {
                    in_process = (start[i] == t);
                  } else {
                    in_process = (Member(start[i], t-process_time+1, t));
                  }

                  //in_process = ((start[i] <= t) && (start[i] > t-process_time));
                  //std::cout << "  -task t" << i << " is in process if it start before or at " 
                  //          << t << " and ends after or at " << t << ": " << in_process << " ";
                } else {
                  // must finish at t or after
                  in_process = (start[i] <= t);
                  //std::cout << "  -task t" << i << " is in process if it start before or at " 
                  //          << t << ": " << in_process << " ";
                }
              } else if(start[i].get_min()+process_time<t) {
                in_process = (start[i] >= t-process_time);
                //std::cout << "  -task t" << i << " is in process if it ends after or at " << t << ": " << in_process << " ";
              }
            } else {
              // TODO: tasks with variable durations
              report_unsupported("p_cumulative");
            }
            if(req[i].is_ground()) {
              //std::cout << "the demand is fixed: " << req[i].get_min() << std::endl;
              scope.add(in_process);
              demand.add(req[i].get_min());
            } else {
              //std::cout << "the demand is variable: " << req[i].get_domain() << std::endl;
              scope.add((in_process * req[i]));
              demand.add(1);
              if(req[i].get_max()>1) {
                Boolean = false;
              }
            }
          } else {
            //std::cout << "  -task t" << i << " cannot be in process\n";
          }
        }
        
        //std::cout << scope.size << " / " << start.size << std::endl;
        //std::cout << scope.size << " > 1 AND " << total << " >? " << cap.get_min() << " " << Boolean << std::endl;

        if(scope.size>1 && total>cap.get_min()) {
          if(cap.is_ground()) {
            if(Boolean) {
              // check if this is in fact a clause
              if(2*min_req > cap.get_max()) {
                s.add(BoolSum(scope,-Mistral::INFTY,1));
              } else {
                s.add(BoolSum(scope,demand,-Mistral::INFTY,cap.get_min()));
              }
            }
            else
              s.add(Sum(scope,demand,-Mistral::INFTY,cap.get_min()));
          } else {
            scope.add(cap);
            demand.add(-1);
            s.add(Sum(scope,demand,-Mistral::INFTY,0));
          }
        }
      }
    }

    void p_cumulative(Solver& s, FlatZincModel& m,
                      const ConExpr& ce, AST::Node* ann) {

      
      
      Vector<Variable> start = arg2intvarargs(s, m, ce[0]);
      Vector<Variable> dur = arg2intvarargs(s, m, ce[1]);
      Vector<Variable> req = arg2intvarargs(s, m, ce[2]);
      Variable cap = getIntVar(s, m, ce[3]);


 


      if(start.size != dur.size || start.size != req.size) {
        std::cout << " c Warning: wrong arguments for Cumulative - exiting\n";
        exit(0);
      }

 
      int orig=INFTY, horizon = 0, ddate, rdate, unit_req=true, unit_dur=true, disjunctive=false, same_req=true, same_dur=true, var_dur=false, var_req=false;
      //

      for(int i=0; i<start.size; ++i) {
        ddate = start[i].get_max() + dur[i].get_max();
        if(horizon<ddate) {
          horizon = ddate;
        }
        rdate = start[i].get_min();
        if(orig>rdate) {
          orig = rdate;
        }
        if(!req[i].is_ground())
          var_req = true;
        if(!dur[i].is_ground())
          var_dur = true;
        if(req[i].get_max() > 1)
          unit_req = false;
        if(dur[i].get_max() > 1)
          unit_dur = false;
        if(same_dur && (!dur[i].is_ground() || (i && dur[i-1].get_min() != dur[i].get_min())))
          same_dur = false;
        if(same_req && (!req[i].is_ground() || (i && req[i-1].get_min() != req[i].get_min())))
          same_req = false;


        // if(!dur[i].is_ground() || !req[i].is_ground()) {
        //   simple = false;
        // }
      }

      
#ifdef _DEBUG_CUMULATIVE
      std::cout << "\nhorizon  = " << horizon << std::endl
                << "capacity = " << cap << std::endl
                << "tasks    = " << start.size << std::endl;
#endif

      int *order = new int[start.size];
      for(int i=0; i<start.size; ++i) {
        order[i] = i;
      }
      demand.copy(req);
      qsort(order, start.size, sizeof(int), largest_demand);
      demand.neutralise();

      if(start.size>1 && req[order[0]].get_min()+req[order[1]].get_min()>cap.get_max()) {
        disjunctive = true;
      }


#ifdef _DEBUG_CUMULATIVE
      for(int i=0; i<start.size; ++i) {
        int oi = order[i];
        std::cout << "           [" << start[oi].get_min() <<  ".." 
                  <<  start[oi].get_max()+dur[oi].get_max() << "]";
        if(dur[oi].is_ground())
          std::cout << " p=" << dur[oi].get_max() ;
        else
          std::cout << " p=" << dur[oi].get_domain() ;
        if(req[oi].is_ground())
          std::cout << " r=" << req[oi].get_max() ;
        else
          std::cout << " r=" << req[oi].get_domain() ;
        std::cout << std::endl;
      }
      //exit(1);
#endif

      int size_discretization = horizon*start.size;
      int size_flow = start.size*start.size*cap.get_max();


      if(unit_dur && disjunctive) {

#ifdef _DEBUG_CUMULATIVE
        std::cout << "this is an alldiff" << std::endl;
#endif

        s.add(AllDiff(start));
      } else if(unit_dur && unit_req) {

#ifdef _DEBUG_CUMULATIVE
        std::cout << "this is a gcc" << std::endl;
#endif

        int *lb = new int[horizon-orig];
        int *ub = new int[horizon-orig];
        for(int i=orig; i<horizon; ++i) {
          lb[i-orig] = cap.get_min();
          ub[i-orig] = cap.get_max();
        }
        s.add(Occurrences(start, orig, horizon-1, lb, ub));
      } else if(disjunctive && !var_dur) {

#ifdef _DEBUG_CUMULATIVE
        if(same_dur) {
          std::cout << "this is an interdistance" << std::endl;
        } else {
          std::cout << "this is a disjunctive" << std::endl;
        }
#endif

        p_disjunctive(s, start, dur);
      } else {
        
#ifdef _DEBUG_CUMULATIVE
        if(same_dur && same_req) {
          std::cout << "this is a symmetric cumulative" << std::endl;
        } else if(same_dur) {
          std::cout << "this is a cumulative will equal processing times" << std::endl;
        } else if(same_req) {
          std::cout << "this is a cumulative will equal demands" << std::endl;
        } else {
          std::cout << "this is a general cumulative" << std::endl;
        }
        
        std::cout << "size of the time-discretization encoding = " << size_discretization << std::endl;
        std::cout << "size of the flow encoding = " << size_flow << std::endl;
#endif
  
        if(size_flow < size_discretization || var_dur) {
          p_cumulative_flow(s, start, dur, req, cap, horizon);
        } else {
          p_cumulative_discretization(s, start, dur, req, cap, horizon);
        }


      }

    }


    /* global cardinality constraint */
    void p_distribute(Solver& s, FlatZincModel& m,
                      const ConExpr& ce, AST::Node* ann) {
      
      std::cout << "GLOBAL CARDINALITY!!\n";
      exit(1);

    }

    /* coercion constraints */
    void p_bool2int(Solver& s, FlatZincModel& m,
                    const ConExpr& ce, AST::Node* ann) {
      Variable x0 = getBoolVar(s, m, ce[0]);
      Variable x1 = getIntVar(s, m, ce[1]);

      s.add( x0 == x1 );
    }

    void p_int_in(Solver& s, FlatZincModel& m,
                  const ConExpr& ce, AST::Node *) {
      
      if (ce[0]->isIntVar()) {

        Variable x = getIntVar(s, m, ce[0]);
        Vector< int > d = arg2intvec(s,ce[1]); // arg2intargs(ce[1]);
       
        BitSet sd(d[0], d.back(), BitSet::empt);
        for(unsigned int i=0; i<d.size; ++i) {
          sd.add(d[i]);
        }

        //x.set_domain(sd);
        s.add(Member(x,sd));
        //s.add()

    //    Variable y = getIntVar(s, m, ce[0]);
      } // else if (ce[0]->isSetVar()) {

      // } else if (ce[0]->isSetVar()) {

      // }
      //cout << "WARNING: DOES NOT TAKE DOMAINS INTO ACCOUNT" << endl;
      

      // set<int> d = arg2intset(s,ce[1]);
      // for(set<int>::iterator i = d.begin(); i != d.end(); ++i)
      //   cout << *i << " ";
      // cout << endl;
      // if (ce[0]->isBoolVar()) {
      //   Variable x = getBoolVar(s, m, ce[0]);
      //   if( d.find(0) == d.end() ) {
      //     x.setmin(s, 1, NO_REASON);
      //   }
      //   if( d.find(1) == d.end() ) {
      //     x.setmax(s, 0, NO_REASON);
      //   }
      // } else {
      //   Variable x = getIntVar(s, m, ce[0]);
      //   // FIXME: this can be more efficient by traversing the set
      //   for(int i = x.min(s), iend = x.max(s); i != iend; ++i) {
      //     if( d.find(i) == d.end() ) {
      //       x.remove(s, i, NO_REASON);
      //     }
      //   }
      // }
    }

    // /* Bool constraints */
    // Lit safeLit(Solver &s, Variable v) {
    //   assert(v.min(s) >= 0 && v.max(s) <= 1);
    //   if( v.max(s) == 0 )
    //     return ~Lit(v.eqi(s, 0));
    //   else
    //     return Lit(v.eqi(s,1));
    // }

    void p_array_bool_xor(Solver& s, FlatZincModel& m,
                          const ConExpr& ce, AST::Node* ann) {

      Vector<Variable> bv = arg2boolvarargs(s, m, ce[0]);
      if(!bv.empty())
      s.add(Parity(bv, 1));
      

      
      // Vector<Variable> bv = arg2boolvarargs(s, m, ce[0]);
      // Variable count(0, bv.size/2);
      // Vector<int> coefs(bv.size+1);
      // for(unsigned int i=0; i<bv.size; ++i) coefs[i] = 1;
      // coefs[bv.size] = -2;
      // bv.add(count);
        
      // if(!bv.empty())
      //   s.add(Sum(bv, coefs, 1, 1));
      
    }

    void p_array_bool_and(Solver& s, FlatZincModel& m,
                          const ConExpr& ce, AST::Node* ann) {
      
       Vector<Variable> bv = arg2boolvarargs(s, m, ce[0]);
       Variable r = getBoolVar(s, m, ce[1]);

       if(!bv.empty())
         s.add((BoolSum(bv) == bv.size) == r);
       else
         s.add(r == 1);

      // vec<Lit> up;
      // up.push( safeLit(s, r) );
      // for(size_t i = 0; i != bv.size(); ++i) {
      //   vec<Lit> down;
      //   down.push( ~safeLit(s, r) );
      //   down.push( safeLit(s, bv[i]) );
      //   s.addClause(down);

      //   up.push( ~safeLit(s, bv[i]) );
      // }
      // s.addClause(up);
    }

    void p_array_bool_or(Solver& s, FlatZincModel& m,
                         const ConExpr& ce, AST::Node* ann) {
      Vector<Variable> bv = arg2boolvarargs(s, m, ce[0]);
      Variable r = getBoolVar(s, m, ce[1]);


      // std::cout << bv << " >= 1 <-> " << r << std::endl;

      // for(int i=0; i<bv.size; ++i)
      //   std::cout << bv[i].get_domain() << std::endl;

      // std::cout << r.get_domain() << std::endl;

      if(!bv.empty()) {
        if(r.get_min() == 1 && bv.size == 2) {
          s.add( bv[0] || bv[1] );
        } else { 
          s.add((BoolSum(bv) >= 1) == r);
        }
      }
      else
        s.add(r == 0);



      // vec<Lit> up;
      // up.push( ~safeLit(s, r) );
      // for(size_t i = 0; i != bv.size(); ++i) {
      //   vec<Lit> down;
      //   down.push( safeLit(s, r) );
      //   down.push( ~safeLit(s, bv[i]) );
      //   s.addClause(down);

      //   up.push( safeLit(s, bv[i]) );
      // }
      // s.addClause(up);
    }

    void p_bool_and(Solver& s, FlatZincModel& m,
                    const ConExpr& ce, AST::Node* ann) {


      Variable x0 = getBoolVar(s, m, ce[0]);
      Variable x1 = getBoolVar(s, m, ce[1]);
      Variable r = getBoolVar(s, m, ce[2]);

      s.add( (x0 && x1) == r );

      // vec<Lit> ps1, ps2, ps3;
      // ps1.push( ~safeLit(s, r) );
      // ps1.push( safeLit(s, x0) );
      // ps2.push( ~safeLit(s, r) );
      // ps2.push( safeLit(s, x1) );
      // ps3.push( ~safeLit(s, x0) );
      // ps3.push( ~safeLit(s, x1) );
      // ps3.push( safeLit(s, r) );
      // s.addClause(ps1);
      // s.addClause(ps2);
      // s.addClause(ps3);
    }

    void p_bool_or(Solver& s, FlatZincModel& m,
                   const ConExpr& ce, AST::Node* ann) {
      Variable x0 = getBoolVar(s, m, ce[0]);
      Variable x1 = getBoolVar(s, m, ce[1]);
      Variable r = getBoolVar(s, m, ce[2]);

      s.add( (x0 || x1) == r );

      /*Can be decomposed into
       * not(r) \/ x0  \/x1
       * r  \/ not(x0)
       * r  \/ not(x1)
       */

/*

      Vector<Variable *> pos;
      Vector<Variable *> neg;
      //not(r) \/ x0  \/x1
      pos.clear();
      neg.clear();
      pos.add(&x0);
      pos.add(&x1);
      neg.add(&r);
      encode_clause(s,pos,neg );

    // r  \/ not(x0)
      pos.clear();
      neg.clear();
      pos.add(&r);
      neg.add(&x0);
      encode_clause(s,pos,neg );


     // r  \/ not(x1)
      pos.clear();
      neg.clear();
      pos.add(&r);
      neg.add(&x1);
      encode_clause(s,pos,neg );

*/

      // vec<Lit> ps1, ps2, ps3;
      // ps1.push( safeLit(s, r) );
      // ps1.push( ~safeLit(s, x0) );
      // ps2.push( safeLit(s, r) );
      // ps2.push( ~safeLit(s, x1) );
      // ps3.push( safeLit(s, x0) );
      // ps3.push( safeLit(s, x1) );
      // ps3.push( ~safeLit(s, r) );
      // s.addClause(ps1);
      // s.addClause(ps2);
      // s.addClause(ps3);
    }


    void p_bool_clause(Solver& s, FlatZincModel& m,
                       const ConExpr& ce, AST::Node* ann) {

      //report_unsupported("p_bool_clause");
      Vector<Variable> pos = arg2boolvarargs(s, m, ce[0]);
      Vector<Variable> neg = arg2boolvarargs(s, m, ce[1]);

      //
      //Posting a clause instead of boolean sums
      //m.add_clause(pos, neg);


      if(pos.empty())
        s.add( BoolSum(neg) < neg.size );
      else if(neg.empty())
        s.add( BoolSum(pos) > 0 );
      else 
        s.add( (BoolSum(pos) > 0) || (BoolSum(neg) < neg.size) );


      // vector<Variable> x1 = arg2boolvarargs(s, m, ce[1]);

      // vec<Lit> ps;
      // for(size_t i = 0; i != x0.size(); ++i)
      //   ps.push(safeLit(s, x0[i]));
      // for(size_t i = 0; i != x1.size(); ++i)
      //   ps.push(~safeLit(s, x1[i]));
      // s.addClause(ps);
    }

    void p_bool_clause_reif(Solver& s, FlatZincModel& m,
                            const ConExpr& ce, AST::Node* ann) {
    	report_unsupported("p_bool_clause_reif");

      // vector<Variable> x0 = arg2boolvarargs(s, m, ce[0]);
      // vector<Variable> x1 = arg2boolvarargs(s, m, ce[1]);
      // Variable r = getBoolVar(s, m, ce[2]);

      // vec<Lit> ps;
      // ps.push( ~safeLit(s, r) );
      // for(size_t i = 0; i != x0.size(); ++i) {
      //   ps.push(safeLit(s, x0[i]));
      //   vec<Lit> ps1;
      //   ps1.push(~safeLit(s, x0[i]));
      //   ps1.push(safeLit(s, r));
      //   s.addClause(ps1);
      // }
      // for(size_t i = 0; i != x1.size(); ++i) {
      //   ps.push(~safeLit(s, x1[i]));
      //   vec<Lit> ps1;
      //   ps1.push(safeLit(s, x1[i]));
      //   ps1.push(safeLit(s, r));
      //   s.addClause(ps1);
      // }
      // s.addClause(ps);
    }

    void p_bool_eq(Solver& s, FlatZincModel& m,
                   const ConExpr& ce, AST::Node* ann) {

       Variable a = getBoolVar(s, m, ce[0]);
       Variable b = getBoolVar(s, m, ce[1]);

       s.add(a == b);

      // vec<Lit> ps1, ps2;
      // ps1.push(~safeLit(s, a));
      // ps1.push(safeLit(s, b));

      // ps2.push(~safeLit(s, b));
      // ps2.push(safeLit(s, a));

      // s.addClause(ps1);
      // s.addClause(ps2);
    }


    void p_bool_eq_reif(Solver& s, FlatZincModel& m,
                        const ConExpr& ce, AST::Node* ann) {
      Variable a = getBoolVar(s, m, ce[0]);
      Variable b = getBoolVar(s, m, ce[1]);
      Variable r = getBoolVar(s, m, ce[2]);

      s.add((a == b) == r);

      //    p_bool_eq_reif can be decomposed as follows :
      /* r \/ a \/ b
       * r \/ not(a) \/ not(b)
       * not(r) \/ a \/ not(b)
       * not(r) \/ not(a) \/ b
       *
       * I HAD PROBLEMS WITH THE DECOMPOSITION! THE SOLVERS STOPS WITH THE INSTANCE ProjectPlannertest_12_8.mzn [MINIZINC'12 CHALLENGE]
       */

/*
      Vector<Variable *> pos;
      Vector<Variable *> neg;
      pos.clear();
      neg.clear();

      pos.add(&r);
      pos.add(&a);
      pos.add(&b);

      encode_clause(s,pos,neg );

      pos.clear();
      neg.clear();

      pos.add(&r);
      neg.add(&a);
      neg.add(&b);

      encode_clause(s,pos,neg );

      pos.clear();
      neg.clear();

      neg.add(&r);
      pos.add(&a);
      neg.add(&b);

      encode_clause(s,pos,neg );

      pos.clear();
      neg.clear();

      neg.add(&r);
      neg.add(&a);
      pos.add(&b);

      encode_clause(s,pos,neg );
*/


      // vec<Lit> ps1, ps2, ps3, ps4;
      // ps1.push(~safeLit(s, a));
      // ps1.push(~safeLit(s, b));
      // ps1.push(safeLit(s, r));

      // ps2.push(safeLit(s, a));
      // ps2.push(safeLit(s, b));
      // ps2.push(safeLit(s, r));

      // ps3.push(safeLit(s, a));
      // ps3.push(~safeLit(s, b));
      // ps3.push(~safeLit(s, r));

      // ps4.push(~safeLit(s, a));
      // ps4.push(safeLit(s, b));
      // ps4.push(~safeLit(s, r));

      // s.addClause(ps1);
      // s.addClause(ps2);
      // s.addClause(ps3);
      // s.addClause(ps4);
    }

    void p_bool_ge(Solver& s, FlatZincModel& m,
                   const ConExpr& ce, AST::Node* ann) {
      Variable a = getBoolVar(s, m, ce[0]);
      Variable b = getBoolVar(s, m, ce[1]);

      s.add(a >= b);

      // vec<Lit> ps;
      // ps.push( safeLit(s, a) );
      // ps.push( ~safeLit(s, b) );
      // s.addClause(ps);
    }

    void p_bool_ge_reif(Solver& s, FlatZincModel& m,
                        const ConExpr& ce, AST::Node* ann) {
      Variable a = getBoolVar(s, m, ce[0]);
      Variable b = getBoolVar(s, m, ce[1]);
      Variable r = getBoolVar(s, m, ce[2]);

      s.add((a >= b) == r);

      // vec<Lit> ps1, ps2, ps3;
      // ps1.push( safeLit(s, a) );
      // ps1.push( ~safeLit(s, b) );
      // ps1.push( ~safeLit(s, r) );
      // ps2.push( ~safeLit(s, a) );
      // ps2.push( safeLit(s, r) );
      // ps3.push( safeLit(s, b));
      // ps3.push( safeLit(s, r));
      // s.addClause(ps1);
      // s.addClause(ps2);
      // s.addClause(ps3);
    }

    void p_bool_gt(Solver& s, FlatZincModel& m,
                   const ConExpr& ce, AST::Node* ann) {
       Variable a = getBoolVar(s, m, ce[0]);
       Variable b = getBoolVar(s, m, ce[1]);
      
 //      s.add(a > b);
       s.add(a==1);
       s.add(b==0);
       // if( s.value(safeLit(s, a)) == l_False )
      //   throw unsat();
      // if( s.value(~safeLit(s, b)) == l_False )
      //   throw unsat();
      // s.enqueue(safeLit(s, a));
      // s.enqueue(~safeLit(s, b));
    }

    void p_bool_gt_reif(Solver& s, FlatZincModel& m,
                        const ConExpr& ce, AST::Node* ann) {
      Variable a = getBoolVar(s, m, ce[0]);
      Variable b = getBoolVar(s, m, ce[1]);
      Variable r = getBoolVar(s, m, ce[2]);

      s.add((a > b) == r);

      // vec<Lit> ps1, ps2, ps3;
      // ps1.push( ~safeLit(s, a) );
      // ps1.push( safeLit(s, b) );
      // ps1.push( safeLit(s, r) );

      // ps2.push( safeLit(s, a) );
      // ps2.push( ~safeLit(s, r) );
      // ps3.push( ~safeLit(s, b));
      // ps3.push( ~safeLit(s, r));
      // s.addClause(ps1);
      // s.addClause(ps2);
      // s.addClause(ps3);
    }

    void p_bool_le(Solver& s, FlatZincModel& m,
                   const ConExpr& ce, AST::Node* ann) {
    	Variable a = getBoolVar(s, m, ce[0]);
    	Variable b = getBoolVar(s, m, ce[1]);

    	s.add( a <= b );

    	//Add clause not(a) \/ b instead of a <= b
/*
    	Vector<Variable*> pos;
    	Vector<Variable*> neg;
    	pos.clear();
    	neg.clear();
    	pos.add(&b);
    	neg.add(&a);
    	encode_clause(s,pos,neg );
*/
    }
    /* bool_le_reif is : (Not(a) \/ b ) <--> r can be decomposed as follows
       * a \/ r
       * Not(b) \/ r
       * Not(r) \/ Not(a) \/ b
       * However, the solver stops when asking for all solutions! So we use the old constraint
      ( (a <= b) == r)
       */
    void p_bool_le_reif(Solver& s, FlatZincModel& m,
                        const ConExpr& ce, AST::Node* ann) {
      Variable a = getBoolVar(s, m, ce[0]); 
      Variable b = getBoolVar(s, m, ce[1]);
      Variable r = getBoolVar(s, m, ce[2]);

      s.add( (a <= b) == r );


      //cout << "(" << a << " <= " << b << ") = " << r << endl;

      /*
          Vector<Variable *> pos;
            Vector<Variable *> neg;

            pos.clear();
            neg.clear();

            pos.add(&a);
            pos.add(&r);

            encode_clause(s,pos,neg);

            pos.clear();
            neg.clear();

            pos.add(&r);
            neg.add(&b);

            encode_clause(s,pos,neg);

            pos.clear();
            neg.clear();

            pos.add(&b);
            neg.add(&a);
            neg.add(&r);

            encode_clause(s,pos,neg);
      */


    }

    void p_bool_lt(Solver& s, FlatZincModel& m,
                   const ConExpr& ce, AST::Node* ann) {
      Variable a = getBoolVar(s, m, ce[0]); 
      Variable b = getBoolVar(s, m, ce[1]);

      s.add (a==0);
      s.add (b==1);
      //s.add( a < b );
    }

    void p_bool_lt_reif(Solver& s, FlatZincModel& m,
                        const ConExpr& ce, AST::Node* ann) {
      Variable a = getBoolVar(s, m, ce[0]);
      Variable b = getBoolVar(s, m, ce[1]);
      Variable r = getBoolVar(s, m, ce[2]);

      s.add( (a < b) == r );
    }

    void p_bool_left_imp(Solver& s, FlatZincModel& m,
                         const ConExpr& ce, AST::Node* ann) {
      p_bool_ge_reif(s, m, ce, ann);
    }

    void p_bool_right_imp(Solver& s, FlatZincModel& m,
                          const ConExpr& ce, AST::Node* ann) {
      p_bool_le_reif(s, m, ce, ann);
    }

    void p_bool_ne(Solver& s, FlatZincModel& m,
                   const ConExpr& ce, AST::Node* ann) {
      Variable a = getBoolVar(s, m, ce[0]);
      Variable b = getBoolVar(s, m, ce[1]);

      s.add( a != b );
    }

    void p_bool_ne_reif(Solver& s, FlatZincModel& m,
                        const ConExpr& ce, AST::Node* ann) {
      Variable a = getBoolVar(s, m, ce[0]);
      Variable b = getBoolVar(s, m, ce[1]);
      Variable r = getBoolVar(s, m, ce[2]);

       s.add( (a != b) == r );

      /* Can be decomposed as follows
       * not(a) \/ b  \/ r
       * a \/ not(b)  \/ r
       * not(a) \/ not(b)  \/ not(r)
       * a \/ b  \/ not(r)
       */
/*
      Vector<Variable *> pos;
      Vector<Variable *> neg;
      // not(a) \/ b  \/ r
      pos.clear();
      neg.clear();

      pos.add(&b);
      pos.add(&r);
      neg.add(&a);

      encode_clause(s,pos,neg);
      // a \/ not(b)  \/ r
      pos.clear();
      neg.clear();

      pos.add(&a);
      pos.add(&r);
      neg.add(&b);

      encode_clause(s,pos,neg);

      //  not(a) \/ not(b)  \/ not(r)
      pos.clear();
      neg.clear();

      neg.add(&a);
      neg.add(&r);
      neg.add(&b);

      encode_clause(s,pos,neg);

      //  a \/ b  \/ not(r)
      pos.clear();
      neg.clear();

      pos.add(&a);
      pos.add(&b);
      neg.add(&r);

      encode_clause(s,pos,neg);
*/
    }

    void p_bool_xor(Solver& s, FlatZincModel& m,
                    const ConExpr& ce, AST::Node* ann) {
      p_bool_ne_reif(s, m, ce, ann);
    }

    void p_bool_not(Solver& s, FlatZincModel& m,
                    const ConExpr& ce, AST::Node* ann) {
      p_bool_ne(s, m, ce, ann);
    }

    /* ================================================== */
    /* Set constraints */

    void p_set_card(Solver& s, FlatZincModel& m,
                    const ConExpr& ce, AST::Node* ann) {
      Variable A = getSetVar(s, m, ce[0]);
      Variable card = getIntVar(s, m, ce[1]);

      s.add(Card(A) == card);
    }

    void p_set_diff(Solver& s, FlatZincModel& m,
                    const ConExpr& ce, AST::Node* ann) {
     
      Variable A = getSetVar(s, m, ce[0]);
      Variable B = getSetVar(s, m, ce[1]);
      Variable C = getSetVar(s, m, ce[2]);

      //std::cout << A << " - " << B << " = " << C << std::endl;
      
      s.add(SetDifference(A,B) == C);
    }

    void p_set_symdiff(Solver& s, FlatZincModel& m,
                       const ConExpr& ce, AST::Node* ann) {
      Variable A = getSetVar(s, m, ce[0]);
      Variable B = getSetVar(s, m, ce[1]);
      Variable C = getSetVar(s, m, ce[2]);
            
      s.add(SymmetricDifference(A,B) == C);
    }

    void p_set_eq(Solver& s, FlatZincModel& m,
                    const ConExpr& ce, AST::Node* ann) {
      Variable A = getSetVar(s, m, ce[0]);
      Variable B = getSetVar(s, m, ce[1]);
      
      s.add( A == B );
    }

    void p_set_ne(Solver& s, FlatZincModel& m,
                    const ConExpr& ce, AST::Node* ann) {
      Variable A = getSetVar(s, m, ce[0]);
      Variable B = getSetVar(s, m, ce[1]);
      
      s.add( A != B );
    }

    void p_set_le(Solver& s, FlatZincModel& m,
                  const ConExpr& ce, AST::Node* ann) {
      Variable A = getSetVar(s, m, ce[0]);
      Variable B = getSetVar(s, m, ce[1]);
      
      s.add( A <= B ); // Warning: this is lex leq (not set_le as defined in minizinc)
    }

    void p_set_lt(Solver& s, FlatZincModel& m,
                  const ConExpr& ce, AST::Node* ann) {
      Variable A = getSetVar(s, m, ce[0]);
      Variable B = getSetVar(s, m, ce[1]);
      
      s.add( A < B ); // Warning: this is lex leq (not set_le as defined in minizinc)
    }

   void p_set_le_reif(Solver& s, FlatZincModel& m,
                  const ConExpr& ce, AST::Node* ann) {
      Variable A = getSetVar(s, m, ce[0]);
      Variable B = getSetVar(s, m, ce[1]);
      Variable b = getBoolVar(s, m, ce[2]);
      report_unsupported("p_set_le_re");
      
      s.add( b == (A <= B) ); // Warning: this is lex leq (not set_le as defined in minizinc)
    }

    void p_set_lt_reif(Solver& s, FlatZincModel& m,
                  const ConExpr& ce, AST::Node* ann) {
      Variable A = getSetVar(s, m, ce[0]);
      Variable B = getSetVar(s, m, ce[1]);
      Variable b = getBoolVar(s, m, ce[2]);
      report_unsupported("p_set_lt_re");
      
      s.add( b == (A < B) ); // Warning: this is lex leq (not set_le as defined in minizinc)
    }

    void p_set_eq_reif(Solver& s, FlatZincModel& m,
                     const ConExpr& ce, AST::Node* ann) {
      Variable A = getSetVar(s, m, ce[0]);
      Variable B = getSetVar(s, m, ce[1]);
      Variable b = getBoolVar(s, m, ce[2]);
 
      s.add( b == (A == B) );
    }

    void p_set_ne_reif(Solver& s, FlatZincModel& m,
                     const ConExpr& ce, AST::Node* ann) {
      Variable A = getSetVar(s, m, ce[0]);
      Variable B = getSetVar(s, m, ce[1]);
      Variable b = getBoolVar(s, m, ce[2]);

      s.add( b == (A != B) );
    }

    void p_set_in(Solver &s, FlatZincModel& m,
                  const ConExpr& ce, AST::Node* ann) {
      Variable x = getIntVar(s, m, ce[0]);

      // std::cout << ce[1]->isSetVar() << std::endl;

      // std::cout << ce[1]->isSet() << std::endl;
      
      
      if(ce[1]->isSetVar()) {
        Variable A = getSetVar(s, m, ce[1]);
        
        // std::cout << x << " should be in " << A.get_domain() << std::endl;
        
        // exit(1);
        
        s.add( Member(x,A));
      } else if(ce[1]->isSet()) {
        
        Vector<int> elts = arg2intvec(s, ce[1]);


        //std::cout << "elts: " << elts << std::endl;
        //fz_set2mistral_set(ce[1], elts);

        s.add( Member(x, elts) );
      }
    }


    void p_set_in_reif(Solver &s, FlatZincModel& m,
                     const ConExpr& ce, AST::Node* ann) {

      

       Variable x = getIntVar(s, m, ce[0]);
       Variable b = getBoolVar(s, m, ce[2]);

      if(ce[1]->isSetVar()) {
        Variable A = getSetVar(s, m, ce[1]);
        s.add( b == Member(x,A));
      } else if(ce[1]->isSet()) {        
        Vector<int> elts = arg2intvec(s, ce[1]);
        s.add( b == Member(x, elts) );
      }

    }

    void p_set_isect(Solver &s, FlatZincModel& m,
                     const ConExpr& ce, AST::Node* ann) {
      Variable A = getSetVar(s, m, ce[0]);
      Variable B = getSetVar(s, m, ce[1]);
      Variable C = getSetVar(s, m, ce[2]);

      s.add( Intersection(A,B) == C );
    }

    void p_set_union(Solver &s, FlatZincModel& m,
                     const ConExpr& ce, AST::Node* ann) {
      Variable A = getSetVar(s, m, ce[0]);
      Variable B = getSetVar(s, m, ce[1]);
      Variable C = getSetVar(s, m, ce[2]);

      s.add( Union(A,B) == C );
    }

    /* Note that we subset in flatzinc is subseteq for us (and same
       for superset). Flatzinc does not have strict subset. */
    void p_set_subset(Solver &s, FlatZincModel& m,
                      const ConExpr& ce, AST::Node* ann) {
      Variable A = getSetVar(s, m, ce[0]);
      Variable B = getSetVar(s, m, ce[1]);

      s.add( Subset(A,B) );
    }

    void p_set_superset(Solver &s, FlatZincModel& m,
                        const ConExpr& ce, AST::Node* ann) {
      Variable A = getSetVar(s, m, ce[0]);
      Variable B = getSetVar(s, m, ce[1]);

      s.add( Subset(B,A) );
    }

    void p_set_subset_reif(Solver &s, FlatZincModel& m,
                      const ConExpr& ce, AST::Node* ann) {
      Variable A = getSetVar(s, m, ce[0]);
      Variable B = getSetVar(s, m, ce[1]);
      Variable p = getBoolVar(s, m, ce[2]);
      
      s.add( Subset(A,B) == p );
    }
 
    void p_set_superset_reif(Solver &s, FlatZincModel& m,
                           const ConExpr& ce, AST::Node* ann) {
      Variable A = getSetVar(s, m, ce[0]);
      Variable B = getSetVar(s, m, ce[1]);
      Variable p = getBoolVar(s, m, ce[2]);
      s.add( p == Subset(B,A) );
    }

    class IntPoster {
    public:
      IntPoster(void) {
        registry().add("int_eq", &p_int_eq);
        registry().add("int_ne", &p_int_neq);
        registry().add("int_ge", &p_int_geq);
        registry().add("int_gt", &p_int_gt);
        registry().add("int_le", &p_int_leq);
        registry().add("int_lt", &p_int_lt);
        registry().add("int_eq_reif", &p_int_eq_reif);
        registry().add("int_ne_reif", &p_int_ne_reif);
        registry().add("int_ge_reif", &p_int_ge_reif);
        registry().add("int_gt_reif", &p_int_gt_reif);
        registry().add("int_le_reif", &p_int_le_reif);
        registry().add("int_lt_reif", &p_int_lt_reif);
        registry().add("int_lin_le", &p_int_lin_le);
        registry().add("int_lin_le_reif", &p_int_lin_le_reif);
        registry().add("int_lin_lt", &p_int_lin_lt);
        registry().add("int_lin_lt_reif", &p_int_lin_lt_reif);
        registry().add("int_lin_ge", &p_int_lin_ge);
        registry().add("int_lin_ge_reif", &p_int_lin_ge_reif);
        registry().add("int_lin_gt", &p_int_lin_gt);
        registry().add("int_lin_gt_reif", &p_int_lin_gt_reif);
        registry().add("int_lin_eq", &p_int_lin_eq);
        registry().add("int_lin_eq_reif", &p_int_lin_eq_reif);
        registry().add("int_lin_ne", &p_int_lin_ne);
        registry().add("int_lin_ne_reif", &p_int_lin_ne_reif);
        registry().add("int_plus", &p_int_plus);
        registry().add("int_minus", &p_int_minus);
        registry().add("int_abs", &p_int_abs);
        registry().add("int_times", &p_int_times);
        registry().add("int_negate", &p_int_negate);
        registry().add("int_min", &p_int_min);
        registry().add("int_max", &p_int_max);
        registry().add("int_div", &p_int_div);
        registry().add("int_mod", &p_int_mod);

        registry().add("int_in", &p_int_in);

        registry().add("array_var_int_element", &p_array_int_element);
        registry().add("array_int_element", &p_array_int_element);
        //registry().add("array_set_element", &p_array_set_element);
        registry().add("array_var_set_element", &p_array_set_element);
        registry().add("array_var_bool_element", &p_array_var_bool_element);
        registry().add("array_bool_element", &p_array_bool_element);

        //Add here Mistral redefinitions of global constraints
        registry().add("distribute", &p_distribute);
        registry().add("all_different_int", &p_all_different);
        registry().add("all_equal_int", &p_all_equal_int);
        registry().add("at_most_int", &p_at_most_int);
        registry().add("at_least_int", &p_at_least_int);
        registry().add("count_eq", &p_count_eq);
        registry().add("count_eq_reif", &p_count_eq_reif);
        registry().add("count_geq", &p_count_geq);
        registry().add("count_gt", &p_count_gt);
        registry().add("count_leq", &p_count_leq);
        registry().add("count_lt", &p_count_lt);
        registry().add("count_neq", &p_count_neq);

        registry().add("cumulative", &p_cumulative);

        registry().add("bool2int", &p_bool2int);

        registry().add("array_bool_xor", &p_array_bool_xor);
        registry().add("array_bool_and", &p_array_bool_and);
        registry().add("array_bool_or", &p_array_bool_or);
        registry().add("bool_and", &p_bool_and);
        registry().add("bool_or", &p_bool_or);
        registry().add("bool_eq", &p_bool_eq);
        registry().add("bool_eq_reif", &p_bool_eq_reif);
        registry().add("bool_ge", &p_bool_ge);
        registry().add("bool_ge_reif", &p_bool_ge_reif);
        registry().add("bool_gt", &p_bool_gt);
        registry().add("bool_gt_reif", &p_bool_gt_reif);
        registry().add("bool_le", &p_bool_le);
        registry().add("bool_le_reif", &p_bool_le_reif);
        registry().add("bool_lt", &p_bool_lt);
        registry().add("bool_lt_reif", &p_bool_lt_reif);
        registry().add("bool_left_imp", &p_bool_left_imp);
        registry().add("bool_right_imp", &p_bool_right_imp);
        registry().add("bool_ne", &p_bool_ne);
        registry().add("bool_ne_reif", &p_bool_ne_reif);
        registry().add("bool_xor", &p_bool_xor);
        registry().add("bool_not", &p_bool_not);
        registry().add("bool_clause", &p_bool_clause);
        registry().add("bool_clause_reif", &p_bool_clause_reif);
        registry().add("bool_lin_eq", &p_bool_lin_eq);
        registry().add("bool_lin_le", &p_bool_lin_le);

        registry().add("set_card", &p_set_card);
        registry().add("set_eq", &p_set_eq);
        registry().add("set_ne", &p_set_ne);
        registry().add("set_eq_reif", &p_set_eq_reif);
        registry().add("set_ne_reif", &p_set_ne_reif);

        //registry().add("set_le", &p_set_le);
        //registry().add("set_le_reif", &p_set_le_reif);

        //registry().add("set_lt", &p_set_le);
        //registry().add("set_lt_reif", &p_set_le_reif);

        registry().add("set_in", &p_set_in);
        registry().add("set_in_reif", &p_set_in_reif);

        registry().add("set_diff", &p_set_diff);
        registry().add("set_symdiff", &p_set_symdiff);
        registry().add("set_intersect", &p_set_isect);
        registry().add("set_union", &p_set_union);

        registry().add("set_subset", &p_set_subset);
        registry().add("set_superset", &p_set_superset);
        registry().add("set_subset_reif", &p_set_subset_reif);
        registry().add("set_superset_reif", &p_set_superset_reif);
      }
    };
    IntPoster __int_poster;
  }
}
