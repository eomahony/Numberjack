%module(package="Numberjack.solvers") Gecode

%{
#include <iostream>
#include <vector>
#include <Gecode.hpp>
#include <gecode/int.hh>
#include <gecode/search.hh>
#include <gecode/driver.hh>
%}

%include "Gecode.hpp"

%template(GecodeExpArray) GecodeArray< Gecode_Expression* >;
%template(GecodeIntArray) GecodeArray< int >;
%template(GecodeDoubleArray) GecodeArray< double >;



%pythoncode %{
import Numberjack

class Solver(Numberjack.NBJ_STD_Solver):
    def __init__(self, model=None, X=None, FD=False, clause_limit=-1, encoding=None):
        Numberjack.NBJ_STD_Solver.__init__(self, "Gecode", "Gecode", model, X, FD, clause_limit, encoding)
%}
