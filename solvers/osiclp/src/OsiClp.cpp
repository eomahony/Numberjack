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

void OsiClpSolver::setTimeLimit(int cutoff) {
	std::cout << "OsiClp does not support setTimeLimit, OsiClpSolverInterface stops after 3600 seconds..." << std::endl;
}
