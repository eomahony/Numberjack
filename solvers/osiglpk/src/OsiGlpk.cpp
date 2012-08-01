#include "OsiGlpk.hpp"

/**************************************************************
 ********************     Solver        ***********************
 **************************************************************/
OsiGlpkSolver::OsiGlpkSolver() : OsiSolver() {
	OsiSolver::setSolver(new OsiGlpkSolverInterface);
}

OsiGlpkSolver::~OsiGlpkSolver() {
}
