#include "MipWrapper.hpp"

#include <iostream>
#include <vector>
#include <math.h>

// solver interfaces
#include <CoinPackedVector.hpp>
#include <CoinPackedMatrix.hpp>
#include <CoinShallowPackedVector.hpp>
#include <OsiSolverInterface.hpp>

template<class T>
class OsiArray {
private:
	std::vector<T> _array;
public:
	OsiArray() {}
	~OsiArray() {}
	int size() {return _array.size();}
	void add(T arg) {_array.push_back(arg);}
	T& get_item(const int i) {return _array[i];}
	void set_item(const int i, T item) {_array[i] = item;}
};
typedef OsiArray<double> OsiDoubleArray;

class Osi_Expression {
private:
	double lb;
	double ub;
public:
	std::string type;
	bool isReal;
	int nbj_ident;

	Osi_Expression() {}
	virtual ~Osi_Expression() {}
	Osi_Expression(const double lb, const double ub);
	bool is_continuous() {return isReal;}
	const char* get_type() {return type.c_str();}
	virtual int get_arity() {return 0.;}
	virtual Osi_Expression* get_child(const int i){return NULL;}
	virtual double get_parameter(const int i){return 0.;}
	double get_min() {return lb;}
	double get_max() {return ub;}
	int getVariableId() {return nbj_ident;}
};

typedef OsiArray<Osi_Expression*> OsiExpArray;

class Osi_DoubleVar: public Osi_Expression {
public:
	Osi_DoubleVar(const double lb, const double ub, const int ident) : Osi_Expression(lb, ub) {nbj_ident = ident; type="var";}
	int get_arity() {return 0;}
	Osi_Expression* get_child(const int i) {return NULL;}
};

class Osi_IntVar: public Osi_Expression {
public:
	Osi_IntVar(const double lb, const double ub, const int ident) : Osi_Expression(lb, ub) {nbj_ident = ident; type="var"; isReal=false;}
	int get_arity() {return 0;}
	Osi_Expression* get_child(const int i) {return NULL;}
};

class Osi_Sum: public Osi_Expression {
private:
	OsiExpArray _vars;
	OsiDoubleArray _weights;

public:
	Osi_Sum(OsiExpArray& vars, OsiDoubleArray& weights, const int offset = 0);
	int get_arity() {return _vars.size();}
	Osi_Expression* get_child(const int i) {return _vars.get_item(i);}
	double get_parameter(const int i) {return _weights.get_item(i);}
};

class Osi_binop: public Osi_Expression {
protected:
	Osi_Expression *_var;
	double _constant;

public:
	Osi_binop(Osi_Expression *var, double constant);
	int get_arity() {return 1;}
	Osi_Expression* get_child(const int i) {return _var;}
	double get_parameter(const int i) { return _constant; }
};

class Osi_le: public Osi_binop {
public:
	Osi_le(Osi_Expression *var, double constant);
};

class Osi_ge: public Osi_binop {
public:
	Osi_ge(Osi_Expression *var, double constant);
};

class Osi_Minimise: public Osi_Expression {
private:
	Osi_Expression* _exp;
public:
	Osi_Minimise(Osi_Expression *var);
	int get_arity() {return 1;}
	Osi_Expression* get_child(const int i) {return _exp;}
};

/**
 The solver itself
 */
class OsiSolver: public MipWrapperSolver {
private:
	std::string solver;
	bool hasSolver;
	OsiSolverInterface *si;
	std::vector<Osi_Expression*> expressions;

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
	inline double manageInfinity(double value);

public:
	bool prepareSolve();
	void setSolver(OsiSolverInterface* s);
	OsiSolverInterface* getSolver();
	OsiSolver();

	int load_gmpl(const char* filename, const char* dataname = NULL);
	int load_mps(const char* filename, const char* extension);
	int load_lp(const char* filename, const double epsilon);

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
	double get_value(void *ptr);
	void build_expressions();
	int num_expression();
	Osi_Expression* get_expression(const int i);
};
