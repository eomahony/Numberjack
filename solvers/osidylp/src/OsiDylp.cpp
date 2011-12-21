#include "OsiDylp.hpp"

/**************************************************************
 ********************     Solver        ***********************
 **************************************************************/
OsiDylpSolver::OsiDylpSolver() : OsiSolver() {
	OsiSolver::setSolver(new OsiDylpSolverInterface);
}

OsiDylpSolver::~OsiDylpSolver() {
}
