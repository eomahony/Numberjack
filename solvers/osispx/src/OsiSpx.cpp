#include "OsiSpx.hpp"

/**************************************************************
 ********************     Solver        ***********************
 **************************************************************/
OsiSpxSolver::OsiSpxSolver() : OsiSolver() {
	OsiSolver::setSolver(new OsiSpxSolverInterface);
}

OsiSpxSolver::~OsiSpxSolver() {
}
