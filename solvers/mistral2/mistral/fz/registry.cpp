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

    /* cumulative */
    void p_cumulative(Solver& s, FlatZincModel& m,
                      const ConExpr& ce, AST::Node* ann) {

      report_unsupported("p_cumulative");
      // vector<Variable> start = arg2intvarargs(s, m, ce[0]);
      // vector<Variable> dur = arg2intvarargs(s, m, ce[1]);
      // vector<Variable> req = arg2intvarargs(s, m, ce[2]);
      // Variable cap = getIntVar(s, m, ce[3]);
      // post_cumulative(s, start, dur, req, cap);
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

        registry().add("all_different_int", &p_all_different);
        //registry().add("cumulative", &p_cumulative);

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
