#include "Osi.hpp"

using namespace std;

/**************************************************************
 ********************     Solver        ***********************
 **************************************************************/
OsiSolver::OsiSolver() {
	_verbosity = 0;
	n_cols = 0;
	n_rows = 0;
	matrix = new CoinPackedMatrix(false, 0, 0);
	hasSolver = false;
}

void OsiSolver::setSolver(OsiSolverInterface* s) {
	si = s;
	hasSolver = true;
}

OsiSolverInterface* OsiSolver::getSolver() {
	return si;
}

int OsiSolver::load_gmpl(char *filename, char *dataname) {
	return dataname == NULL ?
			si->readGMPL(filename) : si->readGMPL(filename, dataname);
}

OsiSolver::~OsiSolver() {
}

void OsiSolver::initialise(MipWrapperExpArray& arg) {
	initialise();
}

void OsiSolver::initialise() {
	MipWrapperSolver::initialise();
	objective = new double[var_counter];

	col_lb = new double[var_counter];
	col_ub = new double[var_counter];
	row_lb = new double[_constraints.size()];row_ub
	= new double[_constraints.size()];

if(	_obj != NULL) {
		for (unsigned int i = 0; i < _obj->_coefficients.size(); i++) {
			objective[i] = -1.0 * _obj_coef * _obj->_coefficients[i];
		}
	}
	for (unsigned int i = 0; i < _constraints.size(); ++i)
		add_in_constraint(_constraints[i]);
}

inline double OsiSolver::manageInfinity(double value) {
	if (hasSolver)
		return (isinf(value) ? ((value > 0 ? 1. : -1.) * si->getInfinity()) :
		value)
		;
	else {
		printf("No OSI set\n");
		return value;
	}
}

void OsiSolver::add_in_constraint(LinearConstraint *con, double coef) {
	CoinPackedVector row;
	for (unsigned int i = 0; i < con->_variables.size(); ++i) {
		int *index = new int;
		MipWrapper_Expression *var = con->_variables[i];
		if (var->_var == NULL) {
			*index = n_cols++;
			var->_var = index;
			if (!var->_continuous) {
				integer_vars.push_back(*index);
			}
		} else {
			index = (int*) var->_var;
		}
		col_lb[*index] = manageInfinity(var->_lower);
		col_ub[*index] = manageInfinity(var->_upper);
		row.insert(*index, con->_coefficients[i]);
	}
	row_lb[n_rows] = manageInfinity(con->_lhs);
	row_ub[n_rows] = manageInfinity(con->_rhs);
	n_rows++;
	matrix->appendRow(row);
}

void OsiSolver::printModel() {
	printf("########\nMODEL:\nVARS :\n");
	printf("V%d(%.2lf, %.2lf)", 0, col_lb[0], col_ub[0]);
	for (int i = 1; i < n_cols; i++) {
		printf(", V%d(%.2lf,%.2lf)", i, col_lb[i], col_ub[i]);
	}

	printf("\n\nObjective : Maximise(%s): ",
			(_obj_coef == -1 ? "Minimise" : "Maximise"));
	printf("%.2lf * V%d", objective[0], 0);
	for (int i = 1; i < n_cols; i++) {
		printf(" + %.2lf * V%d", objective[i], i);
	}

	printf("\n\nMatrix:\n");
	for (int i = 0; i < matrix->getNumRows(); i++) {
		CoinShallowPackedVector row = matrix->getVector(i);
		if (row_lb[i] == -1 * si->getInfinity())
			printf("-infinity <= ");
		else
			printf("%.2lf <= ", row_lb[i]);
		printf("%.2lf * V%d", row.getElements()[0], row.getIndices()[0]);
		for (int j = 1; j < row.getNumElements(); j++)
			printf(" + %.2lf * V%d", row.getElements()[j], row.getIndices()[j]);
		if (row_ub[i] == si->getInfinity())
			printf(" <= infinity\n");
		else
			printf(" <= %.2lf\n", row_ub[i]);
	}
	printf("########\n");
}

/*
 * Loads the model and prepares the solver
 * returns true for success and false when there is no solver set.
 */
bool OsiSolver::prepareSolve() {
	if (n_cols == 0)
		initialise();

	if (!hasSolver) {
		printf("No OsiSolverInterface set\n");
		return false;
	}

	if (_verbosity == 0) {
		si->setHintParam(OsiDoReducePrint);
	} else {
		printModel();
	}
	si->messageHandler()->setLogLevel(_verbosity);

	si->loadProblem(*matrix, col_lb, col_ub, objective, row_lb, row_ub);
	si->setInteger(integer_vars.data(), integer_vars.size());
	return true;
}

int OsiSolver::solve() {
	if (!prepareSolve()) {
		return UNKNOWN;
	}

	try {
		si->initialSolve();
	} catch (CoinError err) {
		err.print(true);
	}

	try {
		si->branchAndBound();
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

int OsiSolver::getNextSolution() {
	return UNSAT;
}

void OsiSolver::setTimeLimit(const int cutoff) {
	std::cout << "Solver does not support time limit" << std::endl;
}

void OsiSolver::setVerbosity(const int degree) {
	_verbosity = degree;
}

bool OsiSolver::is_sat() {
	return si->isProvenOptimal();
}

bool OsiSolver::is_unsat() {
	return !is_sat();
}

void OsiSolver::printStatistics() {
	std::cout << "\td Time: " << getTime() << "\tNodes:" << getNodes()
			<< std::endl;
}

int OsiSolver::getNodes() {
	return 0;
}

double OsiSolver::getTime() {
	return si->getIterationCount();
}

double OsiSolver::get_value(void *ptr) {
	if (hasSolver)
		return si->getColSolution()[*(int*) ptr];
	else {
		printf("No OSI set\n");
		return -1;
	}
}
