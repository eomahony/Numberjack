/** \file Osi.cpp
 \brief Solver interface for PYTHON Wrapper.
 */

#include "Osi.hpp"

/**************************************************************
 ********************     EXPRESSION        *******************
 **************************************************************/

using namespace std;

void Osi_Expression::initialise() {
	hasBeenAdded = false;
}

Osi_Expression::Osi_Expression() {
	cout << "creating an empty expression" << endl;
	initialise();
}

Osi_Expression::Osi_Expression(const int nval) {
	cout << "creating a variable [" << 0 << ".." << (nval - 1) << "]" << endl;
}

Osi_Expression::Osi_Expression(const int lb, const int ub) {
	cout << "creating an integer variable [" << lb << ".." << ub << "]" << endl;
	this->lb = lb;
	this->ub = ub;
}

Osi_Expression::Osi_Expression(const double lb, const double ub) {
	cout << "creating a real variable [" << lb << ".." << ub << "]" << endl;
}

Osi_Expression::Osi_Expression(OsiIntArray& vals) {
	cout << "creating a variable [lb..ub]" << endl;
}

int Osi_Expression::getVariableId() const {
	cout << "return identity of expression" << endl;
	return 0;
}

int Osi_Expression::get_value() const {
	return ub == lb ? ub : 0;
}

int Osi_Expression::get_size() const {
	cout << "return size of expression" << endl;
	return 0;
}

int Osi_Expression::get_max() const {
	cout << "return max of expression" << endl;
	return 0;
}

int Osi_Expression::get_min() const {
	cout << "return min of expression" << endl;
	return 0;
}

bool Osi_Expression::contain(const int v) const {
	cout << "return min of expression" << endl;
	return 0;
}

Osi_Expression::~Osi_Expression() {
	cout << "delete expression" << endl;
}

bool Osi_Expression::has_been_added() const {
	return hasBeenAdded;
}

Osi_Expression* Osi_Expression::add(OsiSolver *solver, bool top_level) {
	cout << "add expression to model" << endl;
	if (top_level) {
		cout << "\tAdding at top level" << endl;
	}

	return this;
}

// /* Binary operators */

Osi_binop::Osi_binop(Osi_Expression *var1, Osi_Expression *var2) :
		Osi_Expression() {
	cout << "creating a binary operator" << endl;
	_vars[0] = var1;
	_vars[1] = var2;
}

Osi_binop::Osi_binop(Osi_Expression *var1, int constant) :
		Osi_Expression() {
	cout << "creating a binary (constant) operator" << endl;
	_vars[0] = var1;
	_vars[1] = NULL;
}

Osi_binop::~Osi_binop() {
	cout << "delete binary operator" << endl;
}

/**
 * Constraints
 */

Osi_Min::Osi_Min(OsiExpArray& vars) :
		Osi_Expression() {
	cout << "creating an Min constraint" << endl;
}

Osi_Min::Osi_Min(Osi_Expression *var1, Osi_Expression *var2) :
		Osi_Expression() {
	cout << "creating a binary Min constraint" << endl;
}

Osi_Min::~Osi_Min() {
	cout << "delete Min" << endl;
}

Osi_Expression* Osi_Min::add(OsiSolver *solver, bool top_level) {
	if (!has_been_added()) {
		solver->exprs.resize(solver->exprs.size() + 1, this);
		hasBeenAdded = true;
		cout << "add a Min constraint to solver" << endl;
		for (int i = 0; i < _vars.size(); ++i)
			_vars.set_item(i, _vars.get_item(i)->add(solver, false));
	}
	return this;
}

Osi_Max::Osi_Max(OsiExpArray& vars) :
		Osi_Expression() {
	cout << "creating an alldiff constraint" << endl;
}

Osi_Max::Osi_Max(Osi_Expression *var1, Osi_Expression *var2) :
		Osi_Expression() {
	cout << "creating a binary alldiff constraint" << endl;
}

Osi_Max::~Osi_Max() {
	cout << "delete alldiff" << endl;
}

Osi_Expression* Osi_Max::add(OsiSolver *solver, bool top_level) {
	if (!has_been_added()) {
		solver->exprs.resize(solver->exprs.size() + 1, this);
		hasBeenAdded = true;
		cout << "add an alldiff constraint" << endl;

		for (int i = 0; i < _vars.size(); ++i)
			_vars.set_item(i, _vars.get_item(i)->add(solver, false));
	}
	return this;
}

Osi_AllDiff::Osi_AllDiff(OsiExpArray& vars) :
		Osi_Expression() {
	cout << "creating an alldiff constraint" << endl;
}

Osi_AllDiff::Osi_AllDiff(Osi_Expression *var1, Osi_Expression *var2) :
		Osi_Expression() {
	cout << "creating a binary alldiff constraint" << endl;
}

Osi_AllDiff::~Osi_AllDiff() {
	cout << "delete alldiff" << endl;
}

Osi_Expression* Osi_AllDiff::add(OsiSolver *solver, bool top_level) {
	if (!has_been_added()) {
		solver->exprs.resize(solver->exprs.size() + 1, this);
		hasBeenAdded = true;
		cout << "add an alldiff constraint" << endl;

		for (int i = 0; i < _vars.size(); ++i)
			_vars.set_item(i, _vars.get_item(i)->add(solver, false));
	}
	return this;
}

Osi_Gcc::Osi_Gcc(OsiExpArray& vars, OsiIntArray& vals, OsiIntArray& lb_card,
		OsiIntArray& ub_card) :
		Osi_Expression() {
	cout << "creating a gcc constraint" << endl;
}

Osi_Gcc::~Osi_Gcc() {
	cout << "delete gcc" << endl;
}

Osi_Expression* Osi_Gcc::add(OsiSolver *solver, bool top_level) {

	if (!has_been_added()) {
		solver->exprs.resize(solver->exprs.size() + 1, this);
		hasBeenAdded = true;
		cout << "add gcc constraint" << endl;

		for (int i = 0; i < _vars.size(); ++i)
			_vars.set_item(i, _vars.get_item(i)->add(solver, false));

	}
	return this;
}

Osi_Element::Osi_Element(OsiExpArray& vars) :
		Osi_Expression() {
	cout << "creating element" << endl;
}

Osi_Element::~Osi_Element() {
	cout << "delete element" << endl;
}

Osi_Expression* Osi_Element::add(OsiSolver *solver, bool top_level) {

	if (!has_been_added()) {
		cout << "add element constraint" << endl;
		solver->exprs.resize(solver->exprs.size() + 1, this);
		hasBeenAdded = true;

		for (int i = 0; i < _vars.size(); ++i)
			_vars.set_item(i, _vars.get_item(i)->add(solver, false));
	}
	return this;
}

Osi_LeqLex::Osi_LeqLex(OsiExpArray& vars) :
		Osi_Expression() {
	cout << "creating lexleq" << endl;
}

Osi_LeqLex::~Osi_LeqLex() {
	cout << "delete leqlex" << endl;
}

Osi_Expression* Osi_LeqLex::add(OsiSolver *solver, bool top_level) {
	if (!has_been_added()) {
		solver->exprs.resize(solver->exprs.size() + 1, this);
		hasBeenAdded = true;
		cout << "add leqlex constraint" << endl;

		for (int i = 0; i < _vars.size(); ++i)
			_vars.set_item(i, _vars.get_item(i)->add(solver, false));
	}
	return this;
}

Osi_LessLex::Osi_LessLex(OsiExpArray& vars) :
		Osi_Expression() {
	cout << "creating lexless" << endl;
}

Osi_LessLex::~Osi_LessLex() {
	cout << "delete leslex" << endl;
}

Osi_Expression* Osi_LessLex::add(OsiSolver *solver, bool top_level) {
	if (!has_been_added()) {
		solver->exprs.resize(solver->exprs.size() + 1, this);
		hasBeenAdded = true;
		cout << "add lesslex constraint" << endl;

		for (int i = 0; i < _vars.size(); ++i)
			_vars.set_item(i, _vars.get_item(i)->add(solver, false));

	}
	return this;
}

Osi_Sum::Osi_Sum(OsiExpArray& vars, OsiIntArray& weights, const int offset) :
		Osi_Expression() {
	for (int i = 0; i < vars.size(); i++) {
		int weight = weights.get_item(i);
		Osi_Expression *expr = vars.get_item(i);
		cout << "creating sum with vars and weights" << endl;
		printf("%d: %d * %d(%d..%d)\n", i, weight, expr->nbj_ident, expr->lb,
				expr->ub);
	}
	_vars = vars;
	_weights = weights;
	_weights.add(offset);
}

Osi_Sum::Osi_Sum(Osi_Expression *arg1, Osi_Expression *arg2,
		OsiIntArray& weights, const int offset) :
		Osi_Expression() {
	cout << "creating sum from two expressions" << endl;
	printf("%d * %d(%d, %d) + %d * %d(%d..%d)\n", weights.get_item(0),
			arg1->nbj_ident, arg1->lb, arg1->ub, weights.get_item(1),
			arg2->nbj_ident, arg2->lb, arg2->ub);
	_vars.add(arg1);
	_vars.add(arg2);
	_weights = weights;
	_weights.add(offset);
}

Osi_Sum::Osi_Sum(Osi_Expression *arg, OsiIntArray& w, const int offset) :
		Osi_Expression() {
	cout << "creating sum" << endl;
	_vars.add(arg);
	_weights = w;
	_weights.add(offset);
}

Osi_Sum::Osi_Sum() :
		Osi_Expression() {
	//_offset = 0;
}

Osi_Sum::~Osi_Sum() {
	cout << "delete sum" << endl;
}

void Osi_Sum::addVar(Osi_Expression* v) {
	cout << "adding variable" << endl;
	_vars.add(v);
}

void Osi_Sum::addWeight(const int w) {
	cout << "adding weight" << endl;
	_weights.add(w);
}

Osi_Expression* Osi_Sum::add(OsiSolver *solver, bool top_level) {
	if (!has_been_added()) {
		solver->exprs.resize(solver->exprs.size() + 1, this);
		hasBeenAdded = true;
		if (solver->obj_id == nbj_ident) {
			int sign = solver->maximise == 1 ? -1 : 1;
			cout << "add " << (sign == 1 ? "Minimise" : "Maximise") << " ";
			for (int i = 0; i < _vars.size(); ++i) {
				Osi_Expression *expr = _vars.get_item(i);
				solver->addVar(expr);
				solver->objective[solver->objective.size() - 1] = sign
						* (double) _weights.get_item(i);
			}
			cout << endl;
		} else {
			cout << "add sum constraint" << endl;
			CoinPackedVector row;
			for (int i = 0; i < _vars.size(); ++i) {
				Osi_Expression *expr = _vars.get_item(i);
				solver->addVar(expr);
				row.insert(expr->nbj_ident, (double) _weights.get_item(i));
			}

			solver->matrix->appendRow(row);
			solver->rows.insert(
					make_pair(nbj_ident, solver->matrix->getNumRows() - 1));
		}

	}
	return this;
}

/**
 * Binary constraints
 */
Osi_mul::Osi_mul(Osi_Expression *var1, Osi_Expression *var2) :
		Osi_binop(var1, var2) {
	cerr << "s Not supported \nc (mulitply by variable) - exiting" << endl;
	exit(1);
}

Osi_mul::Osi_mul(Osi_Expression *var1, int constant) :
		Osi_binop(var1, constant) {
	cout << "creating mul constraint between variable and constant" << endl;
}

Osi_mul::~Osi_mul() {
	cout << "delete mul constraint" << endl;
}

Osi_Expression* Osi_mul::add(OsiSolver *solver, bool top_level) {
	if (!has_been_added()) {
		solver->exprs.resize(solver->exprs.size() + 1, this);
		hasBeenAdded = true;
		cout << "add mul constraint" << endl;
		_vars[0] = _vars[0]->add(solver, false);
		if (top_level) {
			cout << "\tAdding at top level Mul constraint NOT A GOOD IDEA"
					<< endl;
		} else {
			cout << "\tAdding within tree" << endl;
			if (_vars[1]) {
				cout << "\t\t Adding with two variables" << endl;
				_vars[1] = _vars[1]->add(solver, false);
			} else {
				cout << "\t\t Adding with variable and constant" << endl;
			}
		}
	}
	return this;
}

Osi_div::Osi_div(Osi_Expression *var1, Osi_Expression *var2) :
		Osi_binop(var1, var2) {
	cout << "creating div constraint" << endl;
}

Osi_div::Osi_div(Osi_Expression *var1, int constant) :
		Osi_binop(var1, constant) {
	cout << "creating div constraint" << endl;
}

Osi_div::~Osi_div() {
	cout << "delete div" << endl;
}

Osi_Expression* Osi_div::add(OsiSolver *solver, bool top_level) {
	if (!has_been_added()) {
		solver->exprs.resize(solver->exprs.size() + 1, this);
		hasBeenAdded = true;
		cout << "add div constraint" << endl;
		_vars[0] = _vars[0]->add(solver, false);
		if (top_level) {
			cout << "\tAdding at top level Div constraint NOT A GOOD IDEA"
					<< endl;
		} else {
			cout << "\tAdding within tree" << endl;
			if (_vars[1]) {
				cout << "\t\tAdding two variable binary expression" << endl;
				_vars[1] = _vars[1]->add(solver, false);
			} else {
				cout << "\t\tAdding in variable and constraint expression"
						<< endl;
			}
		}
	}
	return this;
}

Osi_mod::Osi_mod(Osi_Expression *var1, Osi_Expression *var2) :
		Osi_binop(var1, var2) {
	cout << "creating mod predicate" << endl;
}

Osi_mod::Osi_mod(Osi_Expression *var1, int constant) :
		Osi_binop(var1, constant) {
	cout << "creating mod predicate" << endl;
}

Osi_mod::~Osi_mod() {
	cout << "delete mod" << endl;
}

Osi_Expression* Osi_mod::add(OsiSolver *solver, bool top_level) {
	if (!has_been_added()) {
		solver->exprs.resize(solver->exprs.size() + 1, this);
		hasBeenAdded = true;
		cout << "add mod predicate" << endl;
		if (top_level) {
			cout << "\tAdding at top level NOT A GOOD IDEA!" << endl;
			/**
			 * Your code goes here. This should probably fail
			 */
		} else {
			cout << "\tAdding within tree" << endl;
			if (_vars[1]) {
				cout << "Adding with two expressions" << endl;
				_vars[1] = _vars[1]->add(solver, false);
				/**
				 *  Your code code here
				 */
			} else {
				cout << "Adding with an expression and a constant" << endl;
				/**
				 * Your code goes here
				 */
			}
		}
	}
	return this;
}

Osi_and::Osi_and(Osi_Expression *var1, Osi_Expression *var2) :
		Osi_binop(var1, var2) {
	cout << "creating and predicate" << endl;
}

Osi_and::Osi_and(Osi_Expression *var1, int constant) :
		Osi_binop(var1, constant) {
	cout << "creating and predicate" << endl;
	cout << "I DON'T THINK I SHOULD BE HERE" << endl;

	/**
	 * Should never be in this constructor???
	 */

}

Osi_and::~Osi_and() {
	cout << "delete and" << endl;
}

Osi_Expression* Osi_and::add(OsiSolver *solver, bool top_level) {
	if (!has_been_added()) {
		solver->exprs.resize(solver->exprs.size() + 1, this);
		hasBeenAdded = true;
		cout << "add and constraint" << endl;

		_vars[0] = _vars[0]->add(solver, false);

		if (top_level) {
			cout << "\tAdding at top level" << endl;

			if (_vars[1]) {

				cout << "\tCreating with two variables\n" << endl;
				_vars[1] = _vars[1]->add(solver, false);

				/**
				 * Your code goes here
				 */

			} else {

				cout << "Adding with a constant" << endl;

				/**
				 * Your code goes here
				 */

			}

		} else {
			cout << "\tAdding within tree" << endl;

			if (_vars[1]) {

				cout << "\tAdding with two variables" << endl;
				_vars[1] = _vars[1]->add(solver, false);

				/**
				 * Your code goes here
				 */

			} else {

				cout << "\tAdding with a constant" << endl;

				/**
				 * Your code goes here
				 */

			}

		}

	}
	return this;
}

Osi_or::Osi_or(Osi_Expression *var1, Osi_Expression *var2) :
		Osi_binop(var1, var2) {
	cout << "creating or predicate" << endl;
}

Osi_or::Osi_or(Osi_Expression *var1, int constant) :
		Osi_binop(var1, constant) {
	cout << "creating or predicate" << endl;
}

Osi_or::~Osi_or() {
	cout << "delete or" << endl;
}

Osi_Expression* Osi_or::add(OsiSolver *solver, bool top_level) {

	if (!has_been_added()) {
		solver->exprs.resize(solver->exprs.size() + 1, this);
		hasBeenAdded = true;
		cout << "add or predicate" << endl;

		_vars[0] = _vars[0]->add(solver, false);

		if (top_level) {
			cout << "\tAdding at top level" << endl;

			if (_vars[1]) {
				cout << "\tAdding with two predicates" << endl;
				_vars[1] = _vars[1]->add(solver, false);

				/**
				 * Your code goes here
				 */

			} else {

				cout << "\tAdding with a constant" << endl;

				/**
				 * Your code goes here
				 */

			}

		} else {
			cout << "\tAdding within tree" << endl;

			if (_vars[1]) {
				cout << "\tAdding with two predicates" << endl;
				_vars[1] = _vars[1]->add(solver, false);

				/**
				 * Your code goes here
				 */

			} else {

				cout << "\tAdding with a constant" << endl;

				/**
				 * Your code goes here
				 */

			}

		}

	}
	return this;
}

Osi_eq::Osi_eq(Osi_Expression *var1, Osi_Expression *var2) :
		Osi_binop(var1, var2) {
	cout << "creating equality" << endl;
	_vars[0] = var1;
	_vars[1] = var2;
}

Osi_eq::Osi_eq(Osi_Expression *var1, int constant) :
		Osi_binop(var1, constant) {
	cout << "creating equality" << endl;
	_vars[0] = var1;
	_constant = constant;
}

Osi_eq::~Osi_eq() {
	cout << "delete eq" << endl;
	_vars[0] = NULL;
	_vars[1] = NULL;
}

Osi_Expression* Osi_eq::add(OsiSolver *solver, bool top_level) {

	if (!has_been_added()) {
		solver->exprs.resize(solver->exprs.size() + 1, this);
		hasBeenAdded = true;
		cout << "add equality constraint" << endl;

		_vars[0] = _vars[0]->add(solver, false);
		CoinPackedVector row;
		row.insert(_vars[0]->nbj_ident, 1.);
		row.insert(_vars[1]->nbj_ident, -1.);

		if (solver->row_lb.size() == solver->row_lb.max_size()) {
			solver->row_lb.resize(solver->row_lb.size(), 0);
			solver->row_ub.resize(solver->row_ub.size(), 0);
		}

		if (top_level) {
			cout << "\tAdding at top level" << endl;

			if (_vars[1]) {
				cout << "\tAdding with two predicates" << endl;
				_vars[1] = _vars[1]->add(solver, false);
				solver->matrix->appendRow(row);

				/**
				 * Your code goes here
				 */

			} else {

				cout << "\tAdding with a constant" << endl;

				/**
				 * Your code goes here
				 */

			}

		} else {
			cout << "\tAdding within tree" << endl;

			if (_vars[1]) {
				cout << "\tAdding with two predicates" << endl;
				_vars[1] = _vars[1]->add(solver, false);

				/**
				 * Your code goes here
				 */

			} else {

				cout << "\tAdding with a constant" << endl;

				/**
				 * Your code goes here
				 */

			}

		}

	}
	return this;
}

/* Disequality operator */

Osi_ne::Osi_ne(Osi_Expression *var1, Osi_Expression *var2) :
		Osi_binop(var1, var2) {
	cout << "creating notequal" << endl;
}

Osi_ne::Osi_ne(Osi_Expression *var1, int constant) :
		Osi_binop(var1, constant) {
	cout << "creating notequal" << endl;
}

Osi_ne::~Osi_ne() {
	cout << "delete notequal" << endl;
}

Osi_Expression* Osi_ne::add(OsiSolver *solver, bool top_level) {

	if (!has_been_added()) {
		solver->exprs.resize(solver->exprs.size() + 1, this);
		hasBeenAdded = true;
		cout << "add notequal constraint" << endl;

		_vars[0] = _vars[0]->add(solver, false);

		if (top_level) {
			cout << "\tAdding at top level" << endl;

			if (_vars[1]) {
				cout << "\tAdding with two predicates" << endl;
				_vars[1] = _vars[1]->add(solver, false);

				/**
				 * Your code goes here
				 */

			} else {

				cout << "\tAdding with a constant" << endl;

				/**
				 * Your code goes here
				 */

			}

		} else {
			cout << "\tAdding within tree" << endl;
			if (_vars[1]) {
				cout << "\tAdding with two predicates" << endl;
				_vars[1] = _vars[1]->add(solver, false);

				/**
				 * Your code goes here
				 */

			} else {

				cout << "\tAdding with a constant" << endl;

				/**
				 * Your code goes here
				 */

			}

		}

	}
	return this;
}

/* Leq operator */

Osi_le::Osi_le(Osi_Expression *var1, Osi_Expression *var2) :
		Osi_binop(var1, var2) {
	cout << "Creating less than or equal" << endl;
}

Osi_le::Osi_le(Osi_Expression *var1, int constant) :
		Osi_binop(var1, constant) {
	cout << "Creating less than or equal" << endl;
}

Osi_le::~Osi_le() {
	cout << "delete lessequal" << endl;
}

Osi_Expression* Osi_le::add(OsiSolver *solver, bool top_level) {
	if (!has_been_added()) {
		solver->exprs.resize(solver->exprs.size() + 1, this);
		hasBeenAdded = true;
		cout << "add less equal constraint" << endl;

		if (top_level) {
			cout << "\tAdding at top level" << endl;

			if (_vars[1]) {
				cout << "\tAdding with two predicates" << endl;
				_vars[1] = _vars[1]->add(solver, false);

				/**
				 * Your code goes here
				 */

			} else {

				cout << "\tAdding with a constant" << endl;

				/**
				 * Your code goes here
				 */

			}

		} else {
			cout << "\tAdding within tree" << endl;

			if (_vars[1]) {
				cout << "\tAdding with two predicates" << endl;
				_vars[1] = _vars[1]->add(solver, false);

				/**
				 * Your code goes here
				 */

			} else {

				cout << "\tAdding with a constant" << endl;

				/**
				 * Your code goes here
				 */

			}

		}
	}
	return this;
}

/* Geq operator */

Osi_ge::Osi_ge(Osi_Expression *var1, Osi_Expression *var2) :
		Osi_binop(var1, var2) {
	cout << "Creating a greater or equal constraint" << endl;
}

Osi_ge::Osi_ge(Osi_Expression *var1, int constant) :
		Osi_binop(var1, constant) {
	cout << "Creating a greater or equal constraint" << endl;
}

Osi_ge::~Osi_ge() {
	cout << "delete greaterequal" << endl;
}

Osi_Expression* Osi_ge::add(OsiSolver *solver, bool top_level) {
	if (!has_been_added()) {
		solver->exprs.resize(solver->exprs.size() + 1, this);
		hasBeenAdded = true;
		cout << "add greaterequal constraint" << endl;

		_vars[0] = _vars[0]->add(solver, false);

		if (top_level) {
			cout << "\tAdding at top level" << endl;
			if (_vars[1]) {
				cout << "\tAdding with two predicates" << endl;
				_vars[1] = _vars[1]->add(solver, false);

				/**
				 * Your code goes here
				 */

			} else {

				cout << "\tAdding with a constant" << endl;

				/**
				 * Your code goes here
				 */

			}

		} else {
			cout << "\tAdding within tree" << endl;
			if (_vars[1]) {
				cout << "\tAdding with two predicates" << endl;
				_vars[1] = _vars[1]->add(solver, false);

				/**
				 * Your code goes here
				 */

			} else {

				cout << "\tAdding with a constant" << endl;

				/**
				 * Your code goes here
				 */

			}

		}

	}

	return this;
}

/* Lt object */

Osi_lt::Osi_lt(Osi_Expression *var1, Osi_Expression *var2) :
		Osi_binop(var1, var2) {
	cout << "Creating a less than constraint" << endl;
}

Osi_lt::Osi_lt(Osi_Expression *var1, int constant) :
		Osi_binop(var1, constant) {
	cout << "Creating a less than constraint" << endl;
	_vars[0] = var1;
	_constant = constant;
}

Osi_lt::~Osi_lt() {
	cout << "delete lessthan" << endl;
}

Osi_Expression* Osi_lt::add(OsiSolver *solver, bool top_level) {
	if (!has_been_added()) {
		solver->exprs.resize(solver->exprs.size() + 1, this);
		hasBeenAdded = true;

		_vars[0] = _vars[0]->add(solver, false);

		int row = solver->rows[_vars[0]->nbj_ident];
		if (solver->row_lb.size() <= (size_t) row) {
			solver->row_lb.resize(row + 1);
			solver->row_ub.resize(row + 1);
		}
		cout << "Adding LT to row :" << row << endl;
		solver->row_lb[row] = -1.0 * solver->si->getInfinity();
		solver->row_ub[row] = (double) (_constant - 1);

		if (top_level) {
			cout << "\tAdding at top level" << endl;

			if (_vars[1]) {
				cout << "\tAdding with two predicates" << endl;
				_vars[1] = _vars[1]->add(solver, false);

				/**
				 * Your code goes here
				 */

			} else {

				cout << "\tAdding with a constant" << endl;

				/**
				 * Your code goes here
				 */

			}

		} else {
			cout << "\tAdding within tree" << endl;

			if (_vars[1]) {
				cout << "\tAdding with two predicates" << endl;
				_vars[1] = _vars[1]->add(solver, false);

				/**
				 * Your code goes here
				 */

			} else {

				cout << "\tAdding with a constant" << endl;

				/**
				 * Your code goes here
				 */

			}
		}

	}
	return this;
}

/* Gt object */

Osi_gt::Osi_gt(Osi_Expression *var1, Osi_Expression *var2) :
		Osi_binop(var1, var2) {
	cout << "Creating a greater than constraint" << endl;
	_vars[0] = var1;
	_vars[1] = var2;
}

Osi_gt::Osi_gt(Osi_Expression *var1, int constant) :
		Osi_binop(var1, constant) {
	cout << "Creating a greater than constraint" << endl;
	_vars[0] = var1;
	_constant = constant;
}

Osi_gt::~Osi_gt() {
	cout << "delete greaterthan" << endl;
}

Osi_Expression* Osi_gt::add(OsiSolver *solver, bool top_level) {
	if (!has_been_added()) {
		solver->exprs.resize(solver->exprs.size() + 1, this);
		hasBeenAdded = true;
		cout << "add greaterthan constraint" << endl;

		_vars[0] = _vars[0]->add(solver, false);

		int row = solver->rows[_vars[0]->nbj_ident];
		if (solver->row_lb.size() <= (size_t) row) {
			solver->row_lb.resize(row + 1);
			solver->row_ub.resize(row + 1);
		}
		cout << "Adding GT to row :" << row << endl;
		solver->row_lb[row] = (double) (_constant - 1);
		solver->row_ub[row] = solver->si->getInfinity();

		if (top_level) {
			cout << "\tAdding at top level" << endl;

			if (_vars[1]) {
				cout << "\tAdding with two predicates" << endl;
				_vars[1] = _vars[1]->add(solver, false);

				/**
				 * Your code goes here
				 */

			} else {

				cout << "\tAdding with a constant" << endl;

				/**
				 * Your code goes here
				 */

			}

		} else {
			cout << "\tAdding within tree" << endl;
			if (_vars[1]) {
				cout << "\tAdding with two predicates" << endl;
				_vars[1] = _vars[1]->add(solver, false);

				/**
				 * Your code goes here
				 */

			} else {

				cout << "\tAdding with a constant" << endl;

				/**
				 * Your code goes here
				 s*/

			}

		}

	}
	return this;
}

/* Minimise object */

Osi_Minimise::Osi_Minimise(Osi_Expression *var) :
		Osi_Expression() {
	_exp = var;
	cout << "Creating a minimise objective" << endl;
}

Osi_Minimise::~Osi_Minimise() {
	cout << "delete minimise" << endl;
}

Osi_Expression* Osi_Minimise::add(OsiSolver *solver, bool top_level) {
	if (!has_been_added()) {
		solver->exprs.resize(solver->exprs.size() + 1, this);
		hasBeenAdded = true;
		cout << "add minimise objective" << endl;

		// This will be the objective function
		solver->maximise = false;
		solver->obj_id = _exp->nbj_ident;
		_exp = _exp->add(solver, false);

		if (top_level) {
			cout << "\tAdding at top level" << endl;

			/**
			 * Your code goes here
			 */

		} else {
			cout << "\tAdding within tree" << endl;
			cout << "SHOULD NOT BE HERE!!!" << endl;
		}

	}
	return this;
}

/* Maximise object */

Osi_Maximise::Osi_Maximise(Osi_Expression *var) :
		Osi_Expression() {
	cout << "Creating a maximise objective" << endl;
	_exp = var;
}

Osi_Maximise::~Osi_Maximise() {
	cout << "delete maximise" << endl;
	_exp = NULL;
}

Osi_Expression* Osi_Maximise::add(OsiSolver *solver, bool top_level) {
	if (!has_been_added()) {
		solver->exprs.resize(solver->exprs.size() + 1, this);
		hasBeenAdded = true;
		cout << "add maximise objective" << endl;

		solver->obj_id = _exp->nbj_ident;
		solver->maximise = true;
		// This will be the objective function
		_exp = _exp->add(solver, false);

		if (top_level) {
			cout << "\tAdding at top level" << endl;

			/**
			 * Your code goes here
			 */

		} else {
			cout << "\tAdding within tree" << endl;
			cout << "SHOULD NOT BE HERER!!!" << endl;
		}

	}
	return this;
}

/**************************************************************
 ********************     Solver        ***********************
 **************************************************************/

OsiSolver::OsiSolver() {
	cout << "c (wrapper) creating solver" << endl;
	this->si = new OsiXXXSolverInterface;
	this->matrix = new CoinPackedMatrix(false, 0, 0);
}

OsiSolver::~OsiSolver() {
	cout << "c (wrapper) delete solver" << endl;
}

void OsiSolver::add(Osi_Expression* arg) {
	cout << "adding expression" << endl;
	arg->add(this, true);
}

void OsiSolver::initialise(OsiExpArray& arg) {
	cout << "Initialising solver with array of expressions" << endl;
}

void OsiSolver::initialise() {
	cout << "Initialising solver with no expressions" << endl;
}

int OsiSolver::num_vars() {
	return col_lb.size();
}

int OsiSolver::get_degree(int i) {
	cout << "return degree of expression i" << endl;
	return 0;
}

int OsiSolver::solveAndRestart(const int policy, const unsigned int base,
		const double factor, const double decay, const int reinit) {
	cout << " SOLVE!! at level: ..." << endl;
	return 0;
}

int OsiSolver::solve() {
	printf("SOLVE!!\n");

	int nvars = col_lb.size();

	int *integer_vars = new int[nvars];
	for (int i = 0; i < nvars; i++)
		integer_vars[i] = i;

	printModel();

	si->loadProblem(*matrix, col_lb.data(), col_ub.data(), objective.data(),
			row_lb.data(), row_ub.data());
	si->setInteger(integer_vars, nvars);
	si->branchAndBound();
	if (si->isProvenOptimal()) {
		cout << "Found optimal solution!" << std::endl;
		cout << "Objective value is " << si->getObjValue() << std::endl;

		int n = si->getNumCols();
		const double *solution = si->getColSolution();

		for (map<Osi_Expression*, int>::iterator iter = vars.begin();
				iter != vars.end(); iter++) {
			iter->first->lb = solution[iter->second];
			iter->first->ub = solution[iter->second];
		}

		for (int i = 0; i < n; i++)
			cout << "solution[" << i << "] = " << solution[i] << endl;
		return true;
	} else {
		cout << "Didn't find optimal solution." << std::endl;
		// Could then check other status functions.
		return false;
	}
}

void OsiSolver::printModel() {
	int nvars = col_lb.size();
	int nrows = row_lb.size();
	printf("\n########\nMODEL :\nVARS :");
	printf("V%d(%.2lf, %.2lf)", 0, col_lb[0], col_ub[0]);
	for (int i = 1; i < nvars; i++) {
		printf(", V%d(%.2lf,%.2lf)", i, col_lb[i], col_ub[i]);
	}
	printf("\n\nObjective : Minimise: ");
	printf("%.2lf * V%d", objective[0], 0);
	for (int i = 1; i < nvars; i++) {
		printf(" + %.2lf * V%d", objective[i], i);
	}
	printf("\n\nMatrix:\n");
	for (int i = 0; i < nrows; i++) {
		CoinPackedVector row = matrix->getVector(i);
		printf("%.2lf <= ", row_lb[i]);
		double *element = row.getElements();
		int *indices = row.getIndices();
		printf("%.2lf * V%d", element[0], indices[0]);
		for (int j = 1; j < row.getNumElements(); j++) {
			printf(" + %.2lf * V%d", element[j], indices[j]);
		}
		printf(" <= %.2lf\n", row_ub[i]);
	}
	printf("########\n");
}

void OsiSolver::addVar(Osi_Expression *expr) {
	if (vars.find(expr) == vars.end()) {
		col_lb.resize(col_lb.size() + 1, expr->lb);
		col_ub.resize(col_ub.size() + 1, expr->ub);
		objective.resize(objective.size() + 1, 0);

		vars[expr] = col_lb.size() - 1;

		printf("VARS:");
		for (map<Osi_Expression*, int>::iterator iter = vars.begin();
				iter != vars.end(); iter++) {
			printf(", N%d V%d", iter->first->nbj_ident, iter->second);
		}
		printf("\n");
	}
}

int OsiSolver::startNewSearch() {
	cout << "starting new search" << endl;
	return 0;
}

int OsiSolver::getNextSolution() {
	cout << "getting next solution" << endl;
	return 0;
}

int OsiSolver::sacPreprocess(const int type) {
	cout << "running sac preprocessing" << endl;
	return 0;
}

bool OsiSolver::propagate() {
	cout << "propagating" << endl;
	return 0;
}

int OsiSolver::get_level() {
	cout << "return current level of solver" << endl;
	return 0;
}

bool OsiSolver::undo(const int nlevel) {
	cout << "Going back nlevel levels" << endl;
	return 0;
}

void OsiSolver::save() {
	cout << "saving state" << endl;
}

int OsiSolver::next(Osi_Expression* x, int v) {
	cout << "get next" << endl;
	return 0;
}

void OsiSolver::post(const char* op, Osi_Expression* x, int v) {
	cout << "Posting a constraint" << endl;
}

void OsiSolver::deduce(const char* op, Osi_Expression* x, int v) {
	cout << "deducing" << endl;
}

void OsiSolver::deduce() {
	cout << "deducing" << endl;
}

bool OsiSolver::branch_right() {
	cout << "branching right" << endl;
	return 0;
}

void OsiSolver::store_solution() {
	cout << "storing solution" << endl;
}

void OsiSolver::setHeuristic(const char* var_heuristic,
		const char* val_heuristic, const int rand) {
	cout << "Setting heuristics" << endl;
}

void OsiSolver::addNogood(OsiExpArray& vars, OsiIntArray& vals) {
	cout << "Adding nogood constraint" << endl;
}

void OsiSolver::guide(OsiExpArray& vars, OsiIntArray& vals,
		OsiDoubleArray& probs) {
	cout << "Adding nogood constraint" << endl;
}

void OsiSolver::forceFiniteDomain(OsiExpArray& vars) {
	cout << "Forcing finite domain" << endl;
}

void OsiSolver::backtrackTo(const int level) {
	cout << "backtracking to level" << endl;
}

void OsiSolver::presolve() {
	cout << "presolving" << endl;
	si->initialSolve();
}

void OsiSolver::increase_init_level(const int i) {
	cout << "increasing initial level" << endl;
}

void OsiSolver::decrease_init_level(const int i) {
	cout << "decreasing initial level" << endl;
}

void OsiSolver::assign(Osi_Expression *X, const int v) {
	cout << "Setting domain of expression X to v" << endl;
}

void OsiSolver::upOneLevel() {
	cout << "stepping up one level" << endl;
}

void OsiSolver::reset(bool full) {
	si->reset();
}

void OsiSolver::setLowerBounds(OsiExpArray& vars, OsiIntArray& vals) {
	cout << "stepping up one level" << endl;
}

void OsiSolver::setUpperBounds(OsiExpArray& vars, OsiIntArray& vals) {
	cout << "stetting upper bounds" << endl;
}

void OsiSolver::setRestartNogood() {
	cout << "setting restart no good" << endl;
}

void OsiSolver::setFailureLimit(const int cutoff) {
	cout << "setting failure limit" << endl;
}

void OsiSolver::setNodeLimit(const int cutoff) {
	cout << "setting node limit" << endl;
}

void OsiSolver::setTimeLimit(const int cutoff) {
	cout << "setting time limit" << endl;
}

void OsiSolver::setVerbosity(const int degree) {
	cout << "setting verbosity" << endl;
}

void OsiSolver::setRandomized(const int degree) {
	cout << "setting randomised" << endl;
}

void OsiSolver::setRandomSeed(const int seed) {
	cout << "setting random seed" << endl;
}

bool OsiSolver::is_opt() {
	return si->isProvenOptimal();
}

bool OsiSolver::is_sat() {
	return si->isProvenOptimal();
}

bool OsiSolver::is_unsat() {
	return !si->isProvenOptimal();
}

void OsiSolver::printStatistics() {
	cout << "printing statistics" << endl;
}

int OsiSolver::getBacktracks() {
	cout << "return number of backtracks" << endl;
	return 0;
}

int OsiSolver::getNodes() {
	cout << "return number of nodes" << endl;
	return 0;
}

int OsiSolver::getFailures() {
	cout << "return number of failures" << endl;
	return 0;
}

int OsiSolver::getChecks() {
	cout << "return number of checks" << endl;
	return 0;
}

int OsiSolver::getPropags() {
	cout << "return amount of propagation" << endl;
	return 0;
}

double OsiSolver::getTime() {
	cout << "return duration" << endl;
	return 0;
}

int OsiSolver::getNumVariables() {
	return col_lb.size();
}

int OsiSolver::getNumConstraints() {
	return matrix->getNumRows();
}

int OsiSolver::getRandomNumber() {
	return 8;
}
