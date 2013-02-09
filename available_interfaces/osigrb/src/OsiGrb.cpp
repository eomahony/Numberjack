#include "OsiGrb.hpp"

/**************************************************************
 ********************     Solver        ***********************
 **************************************************************/
OsiGrbSolver::OsiGrbSolver() : OsiSolver() {
	clpSi = new OsiGrbSolverInterface;
	OsiSolver::setSolver(clpSi);
}

OsiGrbSolver::~OsiGrbSolver() {
}

