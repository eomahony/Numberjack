#include "OsiDylp.hpp"

/**************************************************************
 ********************     Solver        ***********************
 **************************************************************/
OsiDylpSolver::OsiDylpSolver() :
		OsiSolver() {
	OsiSolver::setSolver(new OsiDylpSolverInterface);
}

OsiDylpSolver::~OsiDylpSolver() {
}

int OsiDylpSolver::solve() {
	OsiSolverInterface* si = OsiSolver::getSolver();

	if (!OsiSolver::prepareSolve()) {
		return UNKNOWN;
	}

	try {
		si->initialSolve();
	} catch (CoinError err) {
		err.print(true);
	}

	if (si->isProvenOptimal())
		return SAT;
	else if (si->isProvenPrimalInfeasible())
		return UNSAT;
	else
		return UNKNOWN;
}
