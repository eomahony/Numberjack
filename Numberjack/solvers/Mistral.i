%module(package="Numberjack.solvers") Mistral

%{
#include <Mistral.hpp>
#include <mistral_mod.h>
#include <mistral_glo.h>
%}

%include "Mistral.hpp"

%template(MistralExpArray) MistralArray< Mistral_Expression* >;
%template(MistralIntArray) MistralArray< int >;
%template(MistralDoubleArray) MistralArray< double >;

%pythoncode %{
import Numberjack

class Solver(Numberjack.NBJ_STD_Solver):
    def __init__(self, model=None, X=None, FD=False, clause_limit=-1, encoding=None):
        Numberjack.NBJ_STD_Solver.__init__(self, "Mistral", "Mistral", model, X, FD, clause_limit, encoding)
%}

