#include "MipWrapper.hpp"

#include <iostream>
#include <vector>

#include <coin/OsiClpSolverInterface.hpp>
#include <coin/CoinPackedVector.hpp>
#include <coin/CoinShallowPackedVector.hpp>
#include <coin/OsiSolverInterface.hpp>

/**
 The solver itself
 */
class OsiSolver: public MipWrapperSolver {
private:

	OsiSolverInterface *si;

	int n_cols;
	double * col_lb; //the column lower bounds
	double * col_ub; //the column upper bounds

	int n_rows;
	double * row_lb;
	double * row_ub;

	double * objective; //the objective coefficients
	std::vector<int> integer_vars; //the indices of the integral vars

	CoinPackedMatrix* matrix;

	std::map<int, int> varIndices;

	int _var_counter;
	int _verbosity;
	void add_in_constraint(LinearConstraint *con, double coef = 0);

public:
	int var_counter;

	OsiSolver();
	virtual ~OsiSolver();

	// initialise the solver before solving (no more calls to add after this)
	void initialise(MipWrapperExpArray& arg);
	void initialise();

	// solving methods
	int solve();
	void printModel();

	// parameter tuning methods
	void setTimeLimit(const int cutoff);
	void setVerbosity(const int degree);

	// statistics methods
	int getNodes();
	bool is_sat();
	bool is_unsat();
	void printStatistics();
	double getTime();

	// Value stuff
	virtual double get_value(void *ptr);

};
