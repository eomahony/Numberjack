#include "OsiClp.hpp"

/**************************************************************
 ********************     Solver        ***********************
 **************************************************************/
OsiClpSolver::OsiClpSolver() : OsiSolver() {
	clpSi = new OsiClpSolverInterface;
	OsiSolver::setSolver(clpSi);
}

OsiClpSolver::~OsiClpSolver() {
}
