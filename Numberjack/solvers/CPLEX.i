%module(package="Numberjack.solvers") CPLEX
%import(module="Numberjack.solvers.MipWrapper") "MipWrapper.hpp"

%{
#include <iostream>
#include <vector>
#include <ilcplex/ilocplex.h>
#include "CPLEX.hpp"
%}

%include "CPLEX.hpp"


%pythoncode %{
import Numberjack

class Solver(Numberjack.NBJ_STD_Solver):
    def __init__(self, model=None, X=None, FD=False, clause_limit=-1, encoding=None):
        Numberjack.NBJ_STD_Solver.__init__(self, "CPLEX", "MipWrapper", model, X, FD, clause_limit, encoding)
%}
