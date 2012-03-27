#include "OsiCbc.hpp"

/**************************************************************
 ********************     Solver        ***********************
 **************************************************************/
OsiCbcSolver::OsiCbcSolver() : OsiSolver() {
	OsiSolver::setSolver(new OsiCbcSolverInterface);
}

OsiCbcSolver::~OsiCbcSolver() {
}

void OsiCbcSolver::setLPRelaxationSolver(OsiSolver njSolver) {
    OsiSolver::setSolver(new OsiCbcSolverInterface(njSolver.getSolver()));
}

int OsiCbcSolver::getNodes() {
	OsiCbcSolverInterface s = (OsiCbcSolverInterface)OsiSolver::getSolver();
	return s.getNodeCount();
}
