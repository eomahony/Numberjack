%module(package="Numberjack.solvers") Mistral2

%{
#include <Mistral2.hpp>
#include <iostream>
#include <vector>
#include <mistral_solver.hpp>
#include <mistral_search.hpp>
%}

%include "Mistral2.hpp"

%template(Mistral2ExpArray) Mistral2Array< Mistral2_Expression* >;
%template(Mistral2IntArray) Mistral2Array< int >;
%template(Mistral2DoubleArray) Mistral2Array< double >;

%pythoncode %{
import Numberjack

class Solver(Numberjack.NBJ_STD_Solver):
    def __init__(self, model=None, X=None, FD=False, clause_limit=-1, encoding=None):
        Numberjack.NBJ_STD_Solver.__init__(self, "Mistral2", "Mistral2", model, X, FD, clause_limit, encoding)
%}
