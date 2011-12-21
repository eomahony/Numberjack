#include "OsiCbc.hpp"

/**************************************************************
 ********************     Solver        ***********************
 **************************************************************/
OsiCbcSolver::OsiCbcSolver() : OsiSolver() {
	OsiSolver::setSolver(new OsiCbcSolverInterface);
}

OsiCbcSolver::~OsiCbcSolver() {
}

int OsiCbcSolver::getNodes() {
	OsiCbcSolverInterface s = (OsiCbcSolverInterface)OsiSolver::getSolver();
	printf("OSICBC GETNODES");
	return s.getNodeCount();
}
