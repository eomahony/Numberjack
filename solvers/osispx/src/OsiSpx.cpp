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

// Override OsiSolver::solve() because Soplex doesn't have branch & bound
int OsiSpxSolver::solve() {
	OsiSolverInterface* si = OsiSolver::getSolver();
	if(!OsiSolver::prepareSolve()) {
		return UNSAT;
	}

	try {
		si->initialSolve();
	} catch (CoinError err) {
		err.print(true);
		return UNSAT;
	}

	if (si->isProvenOptimal())
		return SAT;
	else if (si->isProvenPrimalInfeasible())
		return UNSAT;
	else
		return UNKNOWN;
}
