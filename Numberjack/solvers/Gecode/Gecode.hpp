#include <iostream>
#include <vector>
#include <gecode/driver.hh>
#include <gecode/int.hh>
#include <gecode/search.hh>


#ifndef _PYTHON_H
#define _PYTHON_H


/**
   Array of expressions (for nary constraints)
   These are used to pass in a number of variables or expressions to constraints
*/
template<class T>
class GecodeArray {
private:
    std::vector< T > _array;

public:
    GecodeArray() {}
    virtual ~GecodeArray() {}
    int size() const {
        return _array.size();
    }
    void add(T arg) {
        _array.push_back(arg);
    }
    const T& get_item(const int i) const {
        return _array[i];
    }
    void set_item(const int i, T item) {
        _array[i]=item;
    }
};

typedef GecodeArray< int > GecodeIntArray;
typedef GecodeArray< double > GecodeDoubleArray;

// typedef enum {NOTVAR, BOOLVAR, INTVARNVAL, INTVARBOUNDS, INTVARVALUES} VariableType;
// typedef enum VariableType;
enum VariableType {
    NOTVAR=0,
    BOOLVAR=1,
    INTVARNVAL=2,
    INTVARBOUNDS=3,
    INTVARVALUES=4
};

enum SATType {
    UNKNOWN=0,
    SAT=1,
    UNSAT=2
};

/*
This class acts as a wrapper between the Numberjack and Gecode representations
of the problem. Gecode will return copies of this object for each
solution/during search, so we Numberjack does not maintain pointers to the
Gecode variables like in other solvers, we keep the index into 'variables'.
*/
namespace Gecode{
class Space;
};

class NJGecodeSpace : public Gecode::Space {

protected:
    Gecode::IntVarArray gcvariables;
    std::vector<Gecode::IntVar> variables;
    bool closed;

public:

    NJGecodeSpace(void) {
        std::cout << "In NJGecodeSpace constructor." << std::endl;
        closed = false;
    }

    NJGecodeSpace(bool share, NJGecodeSpace& s) : Gecode::Space(share, s), closed(s.closed) {
        // variables.update(*this, share, s.variables);
        std::cout << "In NJGecodeSpace copy constructor" << std::endl;
        std::cout << "copying:" << s.gcvariables << std::endl;
        this->gcvariables.update(*this, share, s.gcvariables);
        std::cout << "newcopy:" << gcvariables << std::endl;
        // variables.reserve(s.variables.size());
        // for(unsigned int i=0; i<s.variables.size(); i++){
        //     std::cout << "creating copy of var " << s.variables[i] << std::endl;
        //     // variables[i].update(*this, share, s.variables[i]);
        //     // update(*this, share, s.variables[i]);
        //     // Gecode::IntVar x(s.variables[i]);
        //     Gecode::IntVar x;
        //     x.update(*this, share, s.variables[i]);
        //     variables.push_back(x);
        // }
        std::cout << "Finished copy constructor" << std::endl;
    }

    virtual Gecode::Space* copy(bool share) {
        std::cout << "In NJGecodeSpace copy share:" << share << std::endl;
        return new NJGecodeSpace(share, *this);
    }

    // Returns the index into 'variables' of the new variable
    unsigned int createIntVar(int lb, int ub){
// #ifdef _DEBUGWRAP
        std::cout << "Creating gecode variable with bounds " << lb << ".." << ub << std::endl;
        std::cout.flush();
// #endif
        assertNotClosed();

        Gecode::IntVar x(*this, lb, ub);


        variables.push_back(x);
        return variables.size() - 1;
    }

    unsigned int createIntVar(GecodeIntArray vals){
// #ifdef _DEBUGWRAP
        std::cout << "Creating gecode variable with value set of size " << vals.size() << std::endl;
// #endif
        assertNotClosed();
        int *values = new int[vals.size()];
        for(int i=0; i<vals.size(); i++) values[i] = vals.get_item(i);

        Gecode::IntVar x(*this, Gecode::IntSet(*values, vals.size()));

        variables.push_back(x);
        return variables.size() - 1;
    }

    inline Gecode::IntVar getVar(unsigned int ind){
        std::cout << "space:" << (this) << " getVar with ind:" << ind << std::endl;
        if(closed){
            if(ind >= gcvariables.size()){
                std::cerr << "Error: index " << ind << " out of range for variable array in closed model of size " << gcvariables.size() << std::endl;
                exit(1);
            }
            return gcvariables[ind];
        } else {
            if(ind >= variables.size()){
                std::cerr << "Error: index " << ind << " out of range for variable array of size " << variables.size() << std::endl;
                exit(1);
            }
            return variables[ind];
        }
    }

    Gecode::IntVarArray& getGecodeSearchVariables() {
        if(!closed){
            std::cerr << "Error: search variables were requested from the space before the space was closed.";
            exit(1);
        }
        std::cout << "Returning vararry of size " << gcvariables.size() << std::endl;
        return gcvariables;
    }

    void assertNotClosed(){
        if(closed){
            std::cerr << "Error Numberjack space was closed and trying to add something new." << std::endl;
            exit(1);
        }
    }
    void closeNJ(){
        closed = true;
        // gcvariables = Gecode::IntVarArray(*this, variables.size());
        Gecode::IntVarArgs varargs(0);
        for(unsigned int i=0; i<variables.size(); i++){
            // gcvariables[i] = variables[i];
            // varargs << Gecode::IntVar(variables[i]);
            Gecode::IntVar v(variables[i]);
            varargs << v;
        }
        // gcvariables = Gecode::IntVarArray(*this, varargs);
        this->gcvariables = Gecode::IntVarArray(*this, varargs);

        Gecode::IntVarArray vars_labeling(*this, varargs);
        // Gecode::distinct(*this, gcvariables);
        // Gecode::branch(*this, gcvariables, Gecode::INT_VAR_SIZE_MIN(), Gecode::INT_VAL_MIN());
        // Gecode::branch(*this, gcvariables, Gecode::INT_VAR_MIN_MIN(), Gecode::INT_VAL_SPLIT_MIN());
        Gecode::branch(*this, vars_labeling, Gecode::INT_VAR_MIN_MIN(), Gecode::INT_VAL_SPLIT_MIN());

    }

    void print(void) const {
        for(unsigned int i=0; i<variables.size(); i++){
            std::cout << "Variable " << i << ": " << variables[i] << std::endl;
        }
    }
};



/**
   Expression (Used to encode variables & constraints)

   Everything is an expression.
*/
class GecodeSolver;
class Gecode_Expression {

public:

    /**
     * Unique indentifier
     */
    int nbj_ident, _nval, _lb, _ub;
    GecodeIntArray _vals;
    VariableType _vartype;
    int _gcvarid;
    GecodeSolver *_solver;



    /**
     * Returns true if the expression has been added to the underlying solvers
     * used to ensure that things are not added into the solver twice
     */
    bool has_been_added() const;

    /**
     * Creates an expression that has a binary domain
     */
    Gecode_Expression();

    /**
     * Creates an expression that has a domain [0, nval-1] or a domain with
     * nval values
     */
    Gecode_Expression(const int nval);

    /**
     * Creates an expression that has a domain lb to ub
     */
    Gecode_Expression(const int lb, const int ub);

    /**
     * Creates an expression whose domain is specified by a list of values
     */
    Gecode_Expression(GecodeIntArray& vals);

    /**
     * Destructor
     */
    virtual ~Gecode_Expression();

    /**
     * Returns the identifier of this expression
     */
    int getVariableId() const;

    /**
     * Returns the value of the variable
     */
    virtual int get_value() const;

    /**
     * Returns the current size of the domain
     */
    int get_size() const;

    /**
     * Returns the upper bound of the expression
     */
    int get_min() const;

    /**
     * Returs the lower bound of the expression
     */
    int get_max() const;

    /**
     * Returns true iff the value v is in the domain of the variable
     */
    bool contain(const int v) const;

    /**
     * VERY IMPORTANT FUNCTION
     *
     * This method adds the object to the underlying solver.
     *
     * @param solver :- Underlying solver wrapper to which the expression will be
     * added
     * @param top_level :- True iff the expresion is the root expression of an
     * expression tree, false otherwise (leaf or within an expression tree)
     *
     * @return The expression that represents the node in the expression tree.
     * Normally this is the object itself apart from cases when things are encoded
     * or decomposed.
     *
     * Exery class that inherits from expression overloads this method and uses it
     * to ass itself to the underlying solver. Variable does not overload this add
     * method so if Expression::add() is called a variable should be created
     *
     */
    virtual Gecode_Expression* add(GecodeSolver *solver, bool top_level);

    virtual Gecode::IntVar getGecodeVar() const;

};

/**
 * This is the Integer varibale class
 */
class Gecode_IntVar : public Gecode_Expression {

public:

    /**
     * Creates a binary integer variable
     */
    Gecode_IntVar() : Gecode_Expression() {}

    /**
     * Creates an integer variable object
     * whose domain is [lb, ub] with identifier 1
     */
    Gecode_IntVar(const int lb, const int ub, const int ident) : Gecode_Expression(lb, ub) {
        nbj_ident = ident;
    }

    /**
     * Creates an integer variable object whose domain is specified by
     * a list of integers and whose identifier is one
     */
    Gecode_IntVar(GecodeIntArray& vals, const int ident) : Gecode_Expression(vals) {
        nbj_ident = ident;
    }

};

typedef GecodeArray< Gecode_Expression* > GecodeExpArray;

/**
 * Expression that represents the minimum of a set of variables
 */
class Gecode_Min : public Gecode_Expression {
private:
    /**
     * The expressions that are the scope of this constraint
     */
    GecodeExpArray _vars;

public:
    /**
     * Creates a Min expression over an array of variables
     */
    Gecode_Min(GecodeExpArray& vars);

    /**
     * Creates a Min expression over two variabels
     */
    Gecode_Min(Gecode_Expression *var1, Gecode_Expression *var2);

    /**
     * Descructor
     */
    virtual ~Gecode_Min();

    /**
     * Add the Min expression to the underlying solver
     *
     * See: Expression::add()
     *
     */
    virtual Gecode_Expression* add(GecodeSolver *solver, bool top_level);
};

/**
 * Expression that represents the maximum of a set of variables
 */
class Gecode_Max : public Gecode_Expression {
private:

    /**
     * The expressions that are the scope of this constraint
     */
    GecodeExpArray _vars;

public:

    /**
     * Creates a Max expression over an array of variables
     */
    Gecode_Max(GecodeExpArray& vars);

    /**
     * Creates a Max expression over two variabels
     */
    Gecode_Max(Gecode_Expression *var1, Gecode_Expression *var2);

    /**
     * Descructor
     */
    virtual ~Gecode_Max();

    /**
     * Add the Min expression to the underlying solver
     *
     * See: Expression::add()
     *
     */
    virtual Gecode_Expression* add(GecodeSolver *solver, bool top_level);
};

/**
 * AllDifferent constraint
 */
class Gecode_AllDiff : public Gecode_Expression {
private:

    /**
     * Variables that are in the scope of the constraint
     */
    GecodeExpArray _vars;

public:

    /**
     * All Different constraint on an array of expressions
     */
    Gecode_AllDiff(GecodeExpArray& vars);

    /**
     * All Different constraint on two expressions, equivalent to a not equal
     */
    Gecode_AllDiff(Gecode_Expression *var1, Gecode_Expression *var2);

    /**
     * Destructor
     */
    virtual ~Gecode_AllDiff();

    /**
     * Add the all different constraint to the solver
     *
     * see Expression::add()
     */
    virtual Gecode_Expression* add(GecodeSolver *solver, bool top_level);
};

/**
 * GCC constraint
 */
class Gecode_Gcc : public Gecode_Expression {
private:

    /**
     * Variables in scope of constraint
     */
    GecodeExpArray _vars;

    /**
     *
     */
    GecodeIntArray _vals;
    GecodeIntArray _lb_card;
    GecodeIntArray _ub_card;

public:
    Gecode_Gcc(GecodeExpArray& vars,
               GecodeIntArray& vals,
               GecodeIntArray& lb_card,
               GecodeIntArray& ub_card);

    /**
     * Destructor
     */
    virtual ~Gecode_Gcc();

    /**
     * Adds the constraint into the solver
     *
     * see Expression::add()
     */
    virtual Gecode_Expression* add(GecodeSolver *solver, bool top_level);
};

/**
 * Element constraint
 */
class Gecode_Element : public Gecode_Expression {
private:

    /**
     * Variables in the scope of the constraint
     */
    GecodeExpArray _vars;

public:

    /**
     * Create an element constraint object
     *
     * @param vars :-
     */
    Gecode_Element(GecodeExpArray& vars);

    /**
     * Destructor
     */
    virtual ~Gecode_Element();

    /**
     * Initialise internal bounds.
     */
    virtual void initialise();

    /**
     * Adds the constraint into the underlying solver
     *
     * see Expression::add()
     */
    virtual Gecode_Expression* add(GecodeSolver *solver, bool top_level);
};

/**
 * LeqLex constraint
 */
class Gecode_LeqLex : public Gecode_Expression {
private:

    /**
     * Variables in the scope of the constriant
     */
    GecodeExpArray _vars;

public:

    /**
     * Create an LeqLex constraint object
     */
    Gecode_LeqLex(GecodeExpArray& vars);

    /**
     * Destructor
     */
    virtual ~Gecode_LeqLex();

    /**
     * Add the constraint into the solver
     *
     * see Expression::add()
     */
    virtual Gecode_Expression* add(GecodeSolver *solver, bool top_level);
};

/**
 * LessLex constraint
 */
class Gecode_LessLex : public Gecode_Expression {
private:

    /**
     * Variables in the scope of the constraint
     */
    GecodeExpArray _vars;

public:

    /**
     * Create a less lex constraint object
     */
    Gecode_LessLex(GecodeExpArray& vars);

    /**
     * Destructor
     */
    virtual ~Gecode_LessLex();

    /**
     * Add the constraint to the underlying solver
     *
     * see Expression::add()
     */
    virtual Gecode_Expression* add(GecodeSolver *solver, bool top_level);
};

/**
 * Weighted sum expression
 */
class Gecode_Sum: public Gecode_Expression {

private:

    /**
     * Variables in the scope of the constraint
     */
    GecodeExpArray _vars;

    /**
     * Weights of the linear equation
     */
    GecodeIntArray _weights;

    int _offset;


public:

    /**
     * Create a weighted sum expression object that is the weighted sum
     * of the expression array vars and the integer array weights with a optional
     * offset
     */
    Gecode_Sum( GecodeExpArray& vars,
                GecodeIntArray& weights,
                const int offset=0);

    /**
     * Create a weighted sum expression object that is the weighted sum
     * of the two expressions and the integer array weights with an offset
     */
    Gecode_Sum( Gecode_Expression* arg1,
                Gecode_Expression* arg2,
                GecodeIntArray& weights,
                const int offset );

    /**
     * Create a weighted sum expression object that is the weighted sum
     * of the one expression and one array of expressions
     * and the integer array weights with an offset
     */
    Gecode_Sum( Gecode_Expression* arg,
                GecodeIntArray& weights,
                const int offset );

    /**
     * Empty sum object addVar and addWeight methods can be used to add in
     * variables and weights
     */
    Gecode_Sum();

    /**
     * Add a variable to the weighted sum
     */
    void addVar( Gecode_Expression* v );

    /**
     * Add a weight to the weighted sum
     */
    void addWeight( const int w );

    /**
     * Destructor
     */
    virtual ~Gecode_Sum();

    /**
     * Initialise internal bounds.
     */
    virtual void initialise();

    /**
     * Get the value of the expression after solving.
     */
    virtual int get_value() const;

    /**
     * Add the weighted sum to the solver
     *
     * see Expression:add()
     *
     * @return This method needs to return an expression that will represent the
     * value of the weighted sum in the underlying solver
     * be it an intermediate variable or otherwise.
     */
    virtual Gecode_Expression* add(GecodeSolver *solver, bool top_level);
};

/**
 * Class that all binary and unary expressions inherit from
 */
class Gecode_binop: public Gecode_Expression {
protected:

    /**
     * Two expression array that will represent the scope of the expression
     */
    Gecode_Expression *_vars[2];

    /**
     *
     */
    int _constant;

public:

    /**
     * @return 1 iff the expression is unary, 2 otherwise
     */
    int arity() {
        return 1+(_vars[1]==NULL);
    }

    /**
     *
     */
    Gecode_binop(Gecode_Expression *var1,
                 Gecode_Expression *var2);

    /**
     *
     */
    Gecode_binop(Gecode_Expression *var1, int int_arg);

    /**
     * Destructor
     */
    virtual ~Gecode_binop();

    /**
     *
     */
    virtual Gecode_Expression* add(GecodeSolver *solver, bool top_level) = 0;
};

class Gecode_mul: public Gecode_binop {
public:
    Gecode_mul(Gecode_Expression *var1, Gecode_Expression *var2);
    Gecode_mul(Gecode_Expression *var1, int int_arg);
    virtual ~Gecode_mul();
    virtual Gecode_Expression* add(GecodeSolver *solver, bool top_level);
};

class Gecode_div: public Gecode_binop {
public:
    Gecode_div(Gecode_Expression *var1, Gecode_Expression *var2);
    Gecode_div(Gecode_Expression *var1, int int_arg);
    virtual ~Gecode_div();
    virtual Gecode_Expression* add(GecodeSolver *solver, bool top_level);
};

class Gecode_mod: public Gecode_binop {
public:
    Gecode_mod(Gecode_Expression *var1, Gecode_Expression *var2);
    Gecode_mod(Gecode_Expression *var1, int int_arg);
    virtual ~Gecode_mod();
    virtual Gecode_Expression* add(GecodeSolver *solver, bool top_level);
};

class Gecode_and: public Gecode_binop {
public:
    Gecode_and(Gecode_Expression *var1, Gecode_Expression *var2);
    Gecode_and(Gecode_Expression *var1, int int_arg);
    virtual ~Gecode_and();
    virtual Gecode_Expression* add(GecodeSolver *solver, bool top_level);
};

class Gecode_or: public Gecode_binop {
public:
    Gecode_or(Gecode_Expression *var1, Gecode_Expression *var2);
    Gecode_or(Gecode_Expression *var1, int int_arg);
    virtual ~Gecode_or();
    virtual Gecode_Expression* add(GecodeSolver *solver, bool top_level);
};

class Gecode_eq: public Gecode_binop {
public:
    Gecode_eq(Gecode_Expression *var1, Gecode_Expression *var2);
    Gecode_eq(Gecode_Expression *var1, int int_arg);
    virtual ~Gecode_eq();
    virtual Gecode_Expression* add(GecodeSolver *solver, bool top_level);
};

class Gecode_ne: public Gecode_binop {
public:
    Gecode_ne(Gecode_Expression *var1, Gecode_Expression *var2);
    Gecode_ne(Gecode_Expression *var1, int int_arg);
    virtual ~Gecode_ne();
    virtual Gecode_Expression* add(GecodeSolver *solver, bool top_level);
};

class Gecode_NoOverlap: public Gecode_binop {
private:
    int _bonstant;
public:
    Gecode_NoOverlap(Gecode_Expression *var1, Gecode_Expression *var2, int d1, int d2);
    virtual ~Gecode_NoOverlap();
    virtual Gecode_Expression* add(GecodeSolver *solver, bool top_level);
};

class Gecode_le: public Gecode_binop {
public:
    Gecode_le(Gecode_Expression *var1, Gecode_Expression *var2);
    Gecode_le(Gecode_Expression *var1, int int_arg);
    virtual ~Gecode_le();
    virtual Gecode_Expression* add(GecodeSolver *solver, bool top_level);
};

class Gecode_ge: public Gecode_binop {
public:
    Gecode_ge(Gecode_Expression *var1, Gecode_Expression *var2);
    Gecode_ge(Gecode_Expression *var1, int int_arg);
    virtual ~Gecode_ge();
    virtual Gecode_Expression* add(GecodeSolver *solver, bool top_level);
};

class Gecode_lt: public Gecode_binop {
public:
    Gecode_lt(Gecode_Expression *var1, Gecode_Expression *var2);
    Gecode_lt(Gecode_Expression *var1, int int_arg);
    virtual ~Gecode_lt();
    virtual Gecode_Expression* add(GecodeSolver *solver, bool top_level);
};

class Gecode_gt: public Gecode_binop {
public:
    Gecode_gt(Gecode_Expression *var1, Gecode_Expression *var2);
    Gecode_gt(Gecode_Expression *var1, int int_arg);
    virtual ~Gecode_gt();
    virtual Gecode_Expression* add(GecodeSolver *solver, bool top_level);
};


class Gecode_Minimise : public Gecode_Expression {
private:
    Gecode_Expression* _exp;

public:
    Gecode_Minimise(Gecode_Expression *var);
    virtual ~Gecode_Minimise();
    virtual Gecode_Expression* add(GecodeSolver *solver, bool top_level);
};


class Gecode_Maximise : public Gecode_Expression {
private:
    Gecode_Expression* _exp;

public:
    Gecode_Maximise(Gecode_Expression *var);
    virtual ~Gecode_Maximise();
    virtual Gecode_Expression* add(GecodeSolver *solver, bool top_level);
};


/**
   The solver itself
*/
class GecodeSolver {

private:
    SATType lastsolvestatus;

public:

    NJGecodeSpace *gecodespace, *originalspace;

    // int first_decision_level;
    // int saved_level;

    // int _heuristic_randomization;
    // std::string _var_heuristic_str;
    // std::string _val_heuristic_str;
    // std::string _restart_policy_str;

    // Mistral::BranchingHeuristic *_branching_heuristic;
    // Mistral::RestartPolicy *_restart_policy;
    // Mistral::Goal *_search_goal;

    GecodeSolver();
    virtual ~GecodeSolver();

    // add an expression, in the case of a tree of expressions,
    // each node of the tree is added separately, depth first.
    void add(Gecode_Expression* arg);

    // used to initialise search on a given subset of variables
    void initialise(GecodeExpArray& arg);
    // initialise the solver before solving (no more calls to add after this)
    void initialise();

    // solving methods
    int solve();
    int solveAndRestart(const int policy = 0,
                        const unsigned int base = 32,
                        const double factor = 1.3333333,
                        const double decay = 0.0,
                        const int reinit = -1);
    int startNewSearch();
    int getNextSolution();
    // int sacPreprocess(const int type);

    int next(Gecode_Expression* x, int v);

    int get_level();

    bool propagate();
    void save();
    void post(const char* op, Gecode_Expression* x, int v);
    bool undo(const int nlevel);
    void deduce(const char* op, Gecode_Expression* x, int v);
    void deduce();
    bool branch_right();
    void store_solution();

    // parameter tuning methods
    void setHeuristic(const char* var_heuristic, const char* val_heuristic, const int rand);
    void setFailureLimit(const int cutoff);
    void setNodeLimit(const int cutoff);
    void setTimeLimit(const int cutoff);
    void setVerbosity(const int degree);
    void setRandomized(const int degree);
    void setRandomSeed(const int seed);
    void forceFiniteDomain(GecodeExpArray& vars);
    void addNogood(GecodeExpArray& vars,
                   GecodeIntArray& vals);
    void guide(GecodeExpArray& vars,
               GecodeIntArray& vals,
               GecodeDoubleArray& probs);
    void backtrackTo(const int level);
    void upOneLevel();
    void presolve();
    void assign(Gecode_Expression *X, const int v);
    void increase_init_level(const int i);
    void decrease_init_level(const int i);
    void reset(bool full);
    void setLowerBounds(GecodeExpArray& vars,
                        GecodeIntArray& vals);
    void setUpperBounds(GecodeExpArray& vars,
                        GecodeIntArray& vals);
    void setRestartNogood();

    // statistics methods
    bool is_opt();
    bool is_sat();
    bool is_unsat();
    void printStatistics();
    int getBacktracks();
    int getNodes();
    int getFailures();
    int getChecks();
    int getPropags();
    double getTime();

    int getRandomNumber();

    int getNumVariables();
    int getNumConstraints();

    int num_vars();
};

#endif // _PYTHON_H
