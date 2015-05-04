%module(package="Numberjack.solvers") MiniSat
%import(module="Numberjack.solvers.SatWrapper") "SatWrapper.hpp"

%{
#include "MiniSat.hpp"
%}

%include "MiniSat.hpp"


%pythoncode %{
import Numberjack

class Solver(Numberjack.NBJ_STD_Solver):
    def __init__(self, model=None, X=None, FD=False, clause_limit=-1, encoding=None):
        Numberjack.NBJ_STD_Solver.__init__(self, "MiniSat", "SatWrapper", model, X, FD, clause_limit, encoding)
%}
