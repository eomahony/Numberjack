#include "OsiVol.hpp"

/**************************************************************
 ********************     Solver        ***********************
 **************************************************************/
OsiVolSolver::OsiVolSolver() : OsiSolver() {
	OsiSolver::setSolver(new OsiVolSolverInterface);
}

OsiVolSolver::~OsiVolSolver() {
}
