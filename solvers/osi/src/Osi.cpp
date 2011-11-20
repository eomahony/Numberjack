#include "Osi.hpp"

/**************************************************************
 ********************     Solver        ***********************
 **************************************************************/

OsiSolver::OsiSolver() {
	n_cols = 0;
	n_rows = 0;
	si = new OsiClpSolverInterface;
	matrix = new CoinPackedMatrix(false, 0, 0);
}

OsiSolver::~OsiSolver() {
}

void OsiSolver::initialise(MipWrapperExpArray& arg) {
	initialise();
}

void OsiSolver::initialise() {
	MipWrapperSolver::initialise();

	objective = new double[_var_counter];

	col_lb = new double[_var_counter];
	col_ub = new double[_var_counter];
	row_lb = new double[_constraints.size()];
	row_ub = new double[_constraints.size()];

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
			printf("Creating new variable:\n");
			*index = n_cols++;
			printf("\tVariable index: %d\n", *index);
			var->_var = index;
			if (!var->_continuous) {
				printf("\tInteger Var;\n");
				integer_vars.push_back(*index);
			}
		} else {
			printf("Found old variable:\n");
			index = (int*) var->_var;
			printf("\tVariable index: %d\n", *index);
		}
		col_lb[*index] = var->_lower;
		col_ub[*index] = var->_upper;
		printf("\tCoefficient: %.2lf\n", con->_coefficients[i]);
		row.insert(*index, con->_coefficients[i]);
	}
	int row_index = n_rows++;
	row_lb[row_index] = con->_lhs;
	row_ub[row_index] = con->_rhs;
	printf("n_rows: %d\n", n_rows);
	matrix->appendRow(row);
}

void OsiSolver::printModel() {
	printf("\n########\nMODEL:\nVARS :\n");
	printf("V%d(%.2lf, %.2lf)", 0, col_lb[0], col_ub[0]);
	for (int i = 1; i < n_cols; i++) {
		printf(", V%d(%.2lf,%.2lf)", i, col_lb[i], col_ub[i]);
	}

	printf("\n\nObjective : Minimise: ");
	printf("%.2lf * V%d", objective[0], 0);
	for (int i = 1; i < n_cols; i++) {
		printf(" + %.2lf * V%d", objective[i], i);
	}

	printf("\n\nMatrix:\n");
	for (int i = 0; i < matrix->getNumRows(); i++) {
		CoinShallowPackedVector row = matrix->getVector(i);
		printf("%.2lf <= ", row_lb[i]);
		for (int j = 0; j < row.getNumElements(); j++)
			printf(" + %.2lf * V%d", row.getElements()[j], row.getIndices()[j]);
		printf(" <= %.2lf\n", row_ub[i]);
	}
	printf("########\n");
}

int OsiSolver::solve() {
	if (n_cols == 0)
		initialise();

	printModel();

	si->loadProblem(*matrix, col_lb, col_ub, objective, row_lb, row_ub);
	si->setInteger(integer_vars.data(), integer_vars.size());
	si->branchAndBound();

	if (si->isProvenOptimal())
		return SAT;
	else if (si->isProvenPrimalInfeasible())
		return UNSAT;
	else
		return UNKNOWN;
}

void OsiSolver::setTimeLimit(const int cutoff) {

}

void OsiSolver::setVerbosity(const int degree) {
//_verbosity = degree;
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
	return 0.;
}

double OsiSolver::get_value(void *ptr) {
	return si->getColSolution()[*(int*) ptr];
}
