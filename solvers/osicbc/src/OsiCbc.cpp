#include "OsiCbc.hpp"

/**************************************************************
 ********************     Solver        ***********************
 **************************************************************/
OsiCbcSolver::OsiCbcSolver() : OsiSolver() {
	OsiSolver::setSolver(new OsiCbcSolverInterface);
}

OsiCbcSolver::~OsiCbcSolver() {
}
