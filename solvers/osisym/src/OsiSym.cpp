#include "OsiSym.hpp"

/**************************************************************
 ********************     Solver        ***********************
 **************************************************************/
OsiSymSolver::OsiSymSolver() : OsiSolver() {
	OsiSolver::setSolver(new OsiSymSolverInterface);
}

OsiSymSolver::~OsiSymSolver() {
}
