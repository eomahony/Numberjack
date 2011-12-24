#include "OsiSpx.hpp"

/**************************************************************
 ********************     Solver        ***********************
 **************************************************************/
OsiSpxSolver::OsiSpxSolver() : OsiSolver() {
	spxS = new OsiSpxSolverInterface;
	OsiSolver::setSolver(spxS);
}

OsiSpxSolver::~OsiSpxSolver() {
}

void OsiSpxSolver::setTimeLimit(const int cutoff) {
	spxS->setTimeLimit(cutoff);
}
