#include "OsiVol.hpp"

/**************************************************************
 ********************     Solver        ***********************
 **************************************************************/
OsiVolSolver::OsiVolSolver() :
		OsiSolver() {
	OsiSolver::setSolver(new OsiVolSolverInterface);
}

OsiVolSolver::~OsiVolSolver() {
}

// Override OsiSolver::solve() because Vol doesn't have a branch & bound
// method and is picky about the type of problem
int OsiVolSolver::solve() {
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
