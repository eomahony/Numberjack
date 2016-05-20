%module(package="Numberjack.solvers") Walksat
%import(module="Numberjack.solvers.SatWrapper") "SatWrapper.hpp"

%{
#include <vector>
#include <iostream>
#include "Walksat.hpp"
#include "cpp_walksat.hpp"
#include "SatWrapper.hpp"
%}

%include "Walksat.hpp"


%pythoncode %{
import Numberjack

class Solver(Numberjack.NBJ_STD_Solver):
    def __init__(self, model=None, X=None, FD=False, clause_limit=-1, encoding=None):
        Numberjack.NBJ_STD_Solver.__init__(self, "Walksat", "SatWrapper", model, X, FD, clause_limit, encoding)
%}

