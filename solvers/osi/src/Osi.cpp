#include "Osi.hpp"

using namespace std;

Osi_Expression::Osi_Expression(const double lb, const double ub) {
	this->lb = lb;
	this->ub = ub;
	this->isReal = true;
	type = "Expression";
}
Osi_binop::Osi_binop(Osi_Expression *var, double constant) {
	_var = var;
	_constant = constant;
}
Osi_le::Osi_le(Osi_Expression *var, double constant) :
		Osi_binop(var, constant) {
	this->type = "le";
}
Osi_ge::Osi_ge(Osi_Expression *var, double constant) :
		Osi_binop(var, constant) {
	this->type = "ge";
}
Osi_Sum::Osi_Sum(OsiExpArray& vars, OsiDoubleArray& weights, const int offset) {
	this->type = "sum";
	_vars = vars;
	_weights = weights;
	_weights.add(offset);
}
Osi_Minimise::Osi_Minimise(Osi_Expression *var) {
	this->type = "minimise";
	_exp = var;
}

/**************************************************************
 ********************     Solver        ***********************
 **************************************************************/
OsiSolver::OsiSolver() {
	_verbosity = 0;
	n_cols = 0;
	n_rows = 0;
	matrix = new CoinPackedMatrix(false, 0, 0);
	hasSolver = false;
	hasIntegers = false;
}

void OsiSolver::setSolver(OsiSolverInterface* s) {
	si = s;
	hasSolver = true;
}

OsiSolverInterface* OsiSolver::getSolver() {
	return si;
}


void OsiSolver::initialise(MipWrapperExpArray& arg) {
	initialise();
}

void OsiSolver::initialise() {
	MipWrapperSolver::initialise();
	objective = new double[var_counter];

	col_lb = new double[var_counter];
	col_ub = new double[var_counter];
	row_lb = new double[_constraints.size()];
	row_ub = new double[_constraints.size()];

	if(_obj != NULL) {
		for (unsigned int i = 0; i < _obj->_coefficients.size(); i++) {
			objective[i] = -1.0 * _obj_coef * _obj->_coefficients[i];
		}
	}
	for (unsigned int i = 0; i < _constraints.size(); ++i)
		add_in_constraint(_constraints[i]);
}

inline double OsiSolver::manageInfinity(double value) {
	if (hasSolver) {
		return (isinf(value)
		? ((value > 0 ? 1. : -1.) * si->getInfinity())
		: value)
		;
	} else {
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
				hasIntegers = true;
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
		const double* elements = row.getElements();
		const int* indices = row.getIndices();
		int n = row.getNumElements();

		if (row_lb[i] == -1 * si->getInfinity())
			printf("-infinity <= ");
		else
			printf("%.2lf <= ", row_lb[i]);
		printf("%.2lf * V%d", elements[0], indices[0]);
		for (int j = 1; j < n; j++)
			printf(" + %.2lf * V%d", elements[j], indices[j]);
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
	//	si->setHintParam(OsiDoReducePrint);
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

	timer.reset();
	try {
		si->initialSolve();
	} catch (CoinError err) {
		err.print(true);
	}

	if(hasIntegers) {
		try {
			si->branchAndBound();
		} catch (CoinError err) {
			err.print(true);
		}
	}
	time = timer.timeElapsed();

	if (si->isProvenOptimal())
		return SAT;
	else if (si->isProvenPrimalInfeasible())
		return UNSAT;
	else
		return UNKNOWN;
}

int OsiSolver::getNextSolution() {
	timer.reset();
	try {
		si->resolve();
	} catch (CoinError err) {
		err.print(true);
	}
	time = timer.timeElapsed();

	if (si->isProvenOptimal())
		return SAT;
	else if (si->isProvenPrimalInfeasible())
		return UNSAT;
	else
		return UNKNOWN;
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
	std::cout << "\td Time: " << getTime()
			  << "\tNodes:" << getNodes()
			  << "\tIterations:" << si->getIterationCount()
			<< std::endl;
}

int OsiSolver::getNodes() {
	std::cout << "Solver does not support reporting of node count" << std::endl;
	return 0;
}

double OsiSolver::getTime() {
	return time;
}

double OsiSolver::get_value(void *ptr) {
	if (hasSolver)
		return si->getColSolution()[*(int*) ptr];
	else {
		printf("No OSI set\n");
		return -1;
	}
}

/**
 * File loading and model building functions.
 *
 * These must be used in conjunction with GMPLParse.py as get_value does not
 * have any variables to reference.
 */

int OsiSolver::load_gmpl(const char *filename, const char *dataname) {
	return dataname == NULL ?
			si->readGMPL(filename) : si->readGMPL(filename, dataname);
}

int OsiSolver::load_mps(const char * filename, const char * extension) {
	return si->readMps(filename, extension);
}

int OsiSolver::load_lp(const char *filename, const double epsilon) {
	return si->readLp(filename, epsilon);
}

void OsiSolver::build_expressions() {
	// build expressions*, first element is Sum() for the objective
	// every pair after that ((1,2),(3,4)) is a le and ge over a
	// common Sum
	const double*objCoef = si->getObjCoefficients();
	const int ncols = si->getNumCols();
	const int nrows = si->getNumRows();
	const double* col_ubs = si->getColUpper();
	const double* col_lbs = si->getColLower();
	const double* row_ubs = si->getRowUpper();
	const double* row_lbs = si->getRowLower();
	const CoinPackedMatrix* mtx = si->getMatrixByRow();

	// Build objective expression.
	OsiExpArray vars;
	OsiDoubleArray coefs;
    double sum = 0;
	for (int i = 0; i < ncols; i++) {
		Osi_Expression* var;

		if(si->isContinuous(i))
			var = new Osi_DoubleVar(col_lbs[i], col_ubs[i], i);
		else
			var = new Osi_IntVar(col_lbs[i], col_ubs[i], i);
		var->varname = si->getColName(i, 255);
		vars.add(var);
		coefs.add(objCoef[i]);
        sum += objCoef[i];
	}
	/* Only return an objective if there is one! (Empty objectives were not
	 * working with some solvers
	 */
    if(sum != 0) {
    	expressions.push_back(new Osi_Minimise(new Osi_Sum(vars, coefs, 0)));
    }

	// Build remaining expressions, expr <= upper, expr >= lower
	for (int i = 0; i < nrows; i++) {
		CoinShallowPackedVector row = mtx->getVector(i);
		const double* elements = row.getElements();
		const int* indices = row.getIndices();
		OsiExpArray sumvars;
		OsiDoubleArray coefs;
		for (int j = 0; j < row.getNumElements(); j++) {
			sumvars.add(vars.get_item(indices[j]));
			coefs.add(elements[j]);
		}
		Osi_Sum* expr = new Osi_Sum(sumvars, coefs, 0);
		expressions.push_back(new Osi_le(expr, row_ubs[i]));
		expressions.push_back(new Osi_ge(expr, row_lbs[i]));
	}
}
int OsiSolver::num_expression() {
	if (expressions.empty())
		build_expressions();
	return expressions.size();
}
Osi_Expression* OsiSolver::get_expression(const int index) {
	if (expressions.empty()) {
		build_expressions();
	}
	return expressions[index];
}
