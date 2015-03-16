/* -*- mode: C++; c-basic-offset: 2; indent-tabs-mode: nil -*- */
/*
 *  Main authors:
 *     Guido Tack <tack@gecode.org>
 *
 *  Copyright:
 *     Guido Tack, 2007
 *
 *  Last modified:
 *     $Date: 2009-09-10 11:44:51 +0200 (Thu, 10 Sep 2009) $ by $Author: tack $
 *     $Revision: 9700 $
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

#ifndef __MINICSP_FLATZINC_REGISTRY_HH__
#define __MINICSP_FLATZINC_REGISTRY_HH__

#include "flatzinc.hpp"
#include <string>
#include <map>

// #define NEQ 1
// #define  EQ 2
// #define GEQ 3
// #define LEQ 4
// #define NONE 5



namespace FlatZinc {

  enum event_type { NEQ_EXPR, EQ_EXPR, GEQ_EXPR, LEQ_EXPR, NONE_EXPR };

  /// Map from constraint identifier to constraint posting functions
  class Registry {
  public:
    /// Type of constraint posting function
    typedef void (*poster) (Solver&,
                            FlatZincModel&,
                            const ConExpr&,
                            AST::Node*);
    /// Add posting function \a p with identifier \a id
    void add(const std::string& id, poster p);
    /// Post constraint specified by \a ce
    void post(Solver& s, FlatZincModel &m,
              const ConExpr& ce, AST::Node* ann);

  private:
    /// The actual registry
    std::map<std::string,poster> r;
  };

  /// Return global registry object
  Registry& registry(void);

}

#endif

// STATISTICS: flatzinc-any
