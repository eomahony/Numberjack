#include "Osi.hpp"
#ifdef _NJ_CBC
#define DEFAULT_SOLVER "cbc"
#elif _NJ_GLPK
#define DEFAULT_SOLVER "glpk"
#else
#define DEFAULT_SOLVER "clp"
#endif
/**************************************************************
 ********************     Solver        ***********************
 **************************************************************/

OsiSolver::OsiSolver() {
	_verbosity = 0;
	n_cols = 0;
	n_rows = 0;
	matrix = new CoinPackedMatrix(false, 0, 0);
	//solver = DEFAULT_SOLVER;
	solver = "glpk";
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
			objective[i] = - _obj->_coefficients[i];
		}
	}
	for (unsigned int i = 0; i < _constraints.size(); ++i)
		add_in_constraint(_constraints[i]);
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
		col_lb[*index] = var->_lower;
		col_ub[*index] = var->_upper;
		row.insert(*index, con->_coefficients[i]);
	}
	int row_index = n_rows++;
	row_lb[row_index] = con->_lhs;
	row_ub[row_index] = con->_rhs;
	matrix->appendRow(row);
}

void OsiSolver::printModel() {
	printf("\n########\nMODEL:\nVARS :\n");
	printf("V%d(%.2lf, %.2lf)", 0, col_lb[0], col_ub[0]);
	for (int i = 1; i < n_cols; i++) {
		printf(", V%d(%.2lf,%.2lf)", i, col_lb[i], col_ub[i]);
	}

	printf("\n\nObjective : Maximise: ");
	printf("%.2lf * V%d", objective[0], 0);
	for (int i = 1; i < n_cols; i++) {
		printf(" + %.2lf * V%d", objective[i], i);
	}

	printf("\n\nMatrix:\n");
	for (int i = 0; i < matrix->getNumRows(); i++) {
		CoinShallowPackedVector row = matrix->getVector(i);
		printf("%.2lf <= ", row_lb[i]);
		printf("%.2lf * V%d", row.getElements()[0], row.getIndices()[0]);
		for (int j = 1; j < row.getNumElements(); j++)
			printf(" + %.2lf * V%d", row.getElements()[j], row.getIndices()[j]);
		printf(" <= %.2lf\n", row_ub[i]);
	}
	printf("########\n");
}

int OsiSolver::solve() {
	if (n_cols == 0)
		initialise();

	printModel();

	if(strcmp(solver.c_str(), "cbc") == 0) {
		si = new OsiCbcSolverInterface;
	} else if(strcmp(solver.c_str(), "glpk") == 0) {
		si = new OsiGlpkSolverInterface;
	} else {
		si = new OsiClpSolverInterface;
	}

	if (_verbosity == 0) {
		si->setHintParam(OsiDoReducePrint);
	}
	si->messageHandler()->setLogLevel(_verbosity);

	si->loadProblem(*matrix, col_lb, col_ub, objective, row_lb, row_ub);
	si->setInteger(integer_vars.data(), integer_vars.size());

	si->initialSolve();
	si->branchAndBound();

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
	return si->getColSolution()[*(int*) ptr];
}
