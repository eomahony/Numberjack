#include "MipWrapper.hpp"

#include <iostream>
#include <vector>

// solver interfaces
#ifdef _NJ_CBC
#include <coin/OsiCbcSolverInterface.hpp>
#include <coin/OsiClpSolverInterface.hpp>
#endif

#ifdef _NJ_CLP
#include <coin/OsiClpSolverInterface.hpp>
#endif

#ifdef _NJ_GLPK
#include <coin/OsiGlpkSolverInterface.hpp>
#endif

#include <coin/CoinPackedVector.hpp>
#include <coin/CoinShallowPackedVector.hpp>
#include <coin/OsiSolverInterface.hpp>

/**
 The solver itself
 */
class OsiSolver: public MipWrapperSolver {
private:

	std::string solver;
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

	int _verbosity;
	void add_in_constraint(LinearConstraint *con, double coef = 0);

public:
	OsiSolver();
	virtual ~OsiSolver();

	// initialise the solver before solving (no more calls to add after this)
	void initialise(MipWrapperExpArray& arg);
	void initialise();

	// solving methods
	int solve();
	int getNextSolution();

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
