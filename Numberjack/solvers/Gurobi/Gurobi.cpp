#include "Gurobi.hpp"
#include <iostream>


/**************************************************************
 ********************     Solver        ***********************
 **************************************************************/

GurobiSolver::GurobiSolver(){
    _verbosity = 0;
    has_been_added = false;
    optimstatus = -1;

    /*  GRBModel takes a copy of env, so to change parameters of the env we
        need to use model->getEnv() and modify that. */
    try {
        env = new GRBEnv();
        model = new GRBModel(*env);
        variables = new vector<GRBVar>;

        setVerbosity(_verbosity);  // Default to no output
    } catch (GRBException e) {
        std::cout << "Gurobi execption while initialising number: " << e.getErrorCode() << std::endl;
        std::cout << e.getMessage() << std::endl;
        exit(1);
    }
}

GurobiSolver::~GurobiSolver(){
    DBG("Delete wrapped solver %s\n", "");

    variables->clear();
    delete variables;
    delete model;
    delete env;
}

void GurobiSolver::initialise(MipWrapperExpArray& arg){
    DBG("Initialise the Solver with search vars %s\n", "");
    // TODO check if we can specify the search variables in Gurobi.
    initialise();
}

void GurobiSolver::initialise(){
    DBG("initialise the solver %s\n", "");
    MipWrapperSolver::initialise();

    // Add all variables
    if (_obj != NULL){
        add_variables_from(_obj, _obj_coef);
    }
    for(unsigned int i = 0; i < _constraints.size(); ++i)
        add_variables_from(_constraints[i]);

    // Update model so we can use variables in constraints
    model->update();

    // Add constraints
    if(_obj != NULL){
        add_in_constraint(_obj, _obj_coef);
    }

    for(unsigned int i = 0; i < _constraints.size(); ++i)
        add_in_constraint(_constraints[i]);
    has_been_added = true;
}

void GurobiSolver::add_variables_from(LinearConstraint *con, double coef) {
    for(unsigned int i = 0; i < con->_variables.size(); ++i){
        if(con->_variables[i]->_var == NULL){
          DBG("\tAdding variable to Gurobi\n%s", "");
        
            GRBVar var_ptr;
            char type;
        
            // GRB_CONTINUOUS, GRB_BINARY, GRB_INTEGER, GRB_SEMICONT, or GRB_SEMIINT
            if(con->_variables[i]->_continuous) type = GRB_CONTINUOUS;
            else type = GRB_INTEGER;

            var_ptr = model->addVar(
                con->_variables[i]->_lower, // LB
                con->_variables[i]->_upper, // UB
                coef,
                type);
            int *var_id = new int;
            *var_id = variables->size();
            variables->push_back(var_ptr);
            con->_variables[i]->_var = (void*) var_id;
        }
    }
}

void GurobiSolver::add_in_constraint(LinearConstraint *con, double coef){
    DBG("Creating a Gurobi representation of a constriant %s\n", "");

    double *weights = new double[con->_coefficients.size()];
    GRBVar *vars = new GRBVar[con->_variables.size()];

    for(unsigned int i = 0; i < con->_variables.size(); ++i){
        assert(con->_variables[i]->_var != NULL);

        vars[i] = variables->at(*(int*)(con->_variables[i]->_var));
        weights[i] = con->_coefficients[i];
    }

    GRBLinExpr lin_expr = GRBLinExpr();
    lin_expr.addTerms(weights, vars, con->_variables.size());

    if(coef < -0.1 || coef > 0.1){
        // Gurobi's objective sense is the opposite of Numberjack's. Gurobi: -1 Minimize, +1 Maximize
        model->set(GRB_IntAttr_ModelSense, -1 * (int)coef);
        model->setObjective(lin_expr);
    } else {
        if(con->_lhs > -INFINITY && con->_rhs < INFINITY){
            if(con->_lhs == con->_rhs) {
                model->addConstr(lin_expr, GRB_EQUAL, con->_lhs);
            } else {
                model->addRange(lin_expr, con->_lhs, con->_rhs);
            }
        } else if(con->_lhs > -INFINITY) model->addConstr(lin_expr, GRB_GREATER_EQUAL, con->_lhs);
        else if(con->_rhs < INFINITY) model->addConstr(lin_expr, GRB_LESS_EQUAL, con->_rhs);
    }
    delete[] weights;
    delete[] vars;
}

int GurobiSolver::solve(){
    DBG("solve!%s\n", "");
    
    if(!has_been_added) initialise();
    model->optimize();
    optimstatus = model->get(GRB_IntAttr_Status);
 
    if( optimstatus == GRB_OPTIMAL || optimstatus == GRB_SUBOPTIMAL ) return SAT;
    else if( optimstatus == GRB_INFEASIBLE ) return UNSAT;
    else return UNKNOWN;
}

void GurobiSolver::setTimeLimit(const int cutoff){
    model->getEnv().set(GRB_DoubleParam_TimeLimit, cutoff);
}

void GurobiSolver::setNodeLimit(const int cutoff){
    model->getEnv().set(GRB_DoubleParam_NodeLimit, cutoff);
}

void GurobiSolver::setThreadCount(const int nr_threads){
    if(nr_threads < 0){
        std::cerr << "Warning: cannot specify a negative thread count, ignoring." << std::endl;
    } else {
        model->getEnv().set(GRB_IntParam_Threads, nr_threads);
    }
}
 
void GurobiSolver::setVerbosity(const int degree){
    // Gurobi's verbosity should be either 0/1.
    if(degree >= 0)
        _verbosity = degree <= 1 ? degree : 1;
    model->getEnv().set(GRB_IntParam_OutputFlag, _verbosity);
}

bool GurobiSolver::is_sat(){
    return (model->get(GRB_IntAttr_Status) == GRB_OPTIMAL || (model->get(GRB_IntAttr_Status) == GRB_SUBOPTIMAL) );
}

bool GurobiSolver::is_unsat(){
    return (model->get(GRB_IntAttr_Status) == GRB_INFEASIBLE);
    return ! is_sat();
}

bool GurobiSolver::is_opt(){
    return optimstatus == GRB_OPTIMAL;
}

void GurobiSolver::output_lp(const char *filename){
    model->write(filename);
}

void GurobiSolver::output_mps(const char *filename){
    model->write(filename);
}

void GurobiSolver::printStatistics(){
    std::cout << "\td Time: " << getTime() << "\tNodes:" << getNodes() << std::endl;
}

int GurobiSolver::getNodes(){
    return (int)(model->get(GRB_DoubleAttr_NodeCount));
}

double GurobiSolver::getTime(){
    return model->get(GRB_DoubleAttr_Runtime);
}

double GurobiSolver::getOptimalityGap(){
    if ((model->get(GRB_IntAttr_SolCount) == 0) || (fabs(model->get(GRB_DoubleAttr_ObjVal)) < 1e-6)) {
        return GRB_INFINITY;
    }
    return fabs(model->get(GRB_DoubleAttr_ObjBound) - model->get(GRB_DoubleAttr_ObjVal)) / fabs(model->get(GRB_DoubleAttr_ObjVal));
}

double GurobiSolver::get_value(void *ptr){
    int index = *(int*)(ptr);
    if(index >= 0 && index < variables->size()){
        try {
            GRBVar v = variables->at(index);
            return v.get(GRB_DoubleAttr_X);
        } catch (GRBException e) {
            std::cout << "Gurobi exception getting value: " << e.getErrorCode() << std::endl;
            std::cout << e.getMessage() << std::endl;
        }
    }
    return 0.0;
}
