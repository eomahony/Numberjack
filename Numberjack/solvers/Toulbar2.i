%module(package="Numberjack.solvers") Toulbar2

%{
#include <Toulbar2.hpp>
#include <tb2solver.hpp>
%}

%include "Toulbar2.hpp"

%template(Toulbar2ExpArray) Toulbar2Array< Toulbar2_Expression* >;
%template(Toulbar2IntArray) Toulbar2Array< int >;
%template(Toulbar2DoubleArray) Toulbar2Array< double >;


%pythoncode %{
import Numberjack

class Solver(Numberjack.NBJ_STD_Solver):
    def __init__(self, model=None, X=None, FD=False, clause_limit=-1, encoding=None):
        Numberjack.NBJ_STD_Solver.__init__(self, "Toulbar2", "Toulbar2", model, X, FD, clause_limit, encoding)
%}
