#include "OsiClp.hpp"

/**************************************************************
 ********************     Solver        ***********************
 **************************************************************/
OsiClpSolver::OsiClpSolver() : OsiSolver() {
	OsiSolver::setSolver(new OsiClpSolverInterface);
}

OsiClpSolver::~OsiClpSolver() {
}
