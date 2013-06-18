#include <iostream>
#include <math.h>
#include "Gurobi.hpp"

using namespace std;

/**************************************************************
 ********************     Solver        ***********************
 **************************************************************/

GurobiSolver::GurobiSolver(){
    // DBG("Create a Gurobi solver %s\n", "");

    var_counter = 0;
    _verbosity = 0;
    has_been_added = false;

    env = new GRBEnv();
    model = new GRBModel(*env);
    variables = new vector<GRBVar>;

    // set the objective sense to maximize, default is minimize
    // SCIP_CALL_EXC( SCIPsetObjsense(_scip, SCIP_OBJSENSE_MAXIMIZE) );
}

GurobiSolver::~GurobiSolver(){
    DBG("Delete wrapped solver %s\n", "");
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

    if(_obj != NULL){
        add_in_constraint(_obj, _obj_coef);
    }

    for(unsigned int i = 0; i < _constraints.size(); ++i)
        add_in_constraint(_constraints[i]);
    has_been_added = true;
}

void GurobiSolver::add_in_constraint(LinearConstraint *con, double coef){
    DBG("Creating a Gurobi representation of a constriant %s\n", "");
    
    double *weights = new double[con->_coefficients.size()];
    // GRBVar** vars = new GRBVar*[con->_variables.size()];
    GRBVar *vars = new GRBVar[con->_variables.size()];
    
    for(unsigned int i = 0; i < con->_variables.size(); ++i){
        
        DBG("\tAdding variable to Gurobi\n%s", "");
        
        if(con->_variables[i]->_var == NULL){
        
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
            vars[i] = var_ptr;
            weights[i] = con->_coefficients[i];
            DBG("Created new variable with id %d. type:%c lb:%f ub:%f coef:%f\n", *var_id, type, con->_variables[i]->_lower, con->_variables[i]->_upper, coef);

            // Need to tell Gurobi to update the model to integrate new variables.
            model->update();
        } else {
            vars[i] = variables->at(*(int*)(con->_variables[i]->_var));
            weights[i] = con->_coefficients[i];
        }
    }

    GRBLinExpr lin_expr = GRBLinExpr();
    lin_expr.addTerms(weights, vars, con->_variables.size());

    if(coef < -0.1 || coef > 0.1){
        // Gurobi's objective sense is the opposite of Numberjack's. Gurobi: -1 Minimize, +1 Maximize
        model->set(GRB_IntAttr_ModelSense, -1 * (int)coef);
        model->setObjective(lin_expr);
    } else {
        model->addRange(lin_expr, con->_lhs, con->_rhs);
    }
}

int GurobiSolver::solve(){
    DBG("solve!%s\n", "");
    
    if(!has_been_added) initialise();

    model->optimize();
    optimstatus = model->get(GRB_IntAttr_Status);
 
    if( optimstatus == GRB_OPTIMAL ) return SAT;
    else if( optimstatus == GRB_INFEASIBLE) return UNSAT;
    else return UNKNOWN;
}

void GurobiSolver::setTimeLimit(const int cutoff){
    env->set(GRB_DoubleParam_TimeLimit, cutoff);
}

void GurobiSolver::setNodeLimit(const int cutoff){
    env->set(GRB_DoubleParam_NodeLimit, cutoff);
}
 
void GurobiSolver::setVerbosity(const int degree){
    // Gurobi's verbosity should be either 0/1.
    if(degree >= 0)
        _verbosity = degree <= 1 ? degree : 1;
    env->set(GRB_IntParam_OutputFlag, _verbosity);
}

bool GurobiSolver::is_sat(){
    return !(model->get(GRB_IntAttr_Status) == GRB_INFEASIBLE);
}

bool GurobiSolver::is_unsat(){
    return ! is_sat();
}

bool GurobiSolver::is_opt(){
    return optimstatus == GRB_OPTIMAL;
}

void GurobiSolver::printStatistics(){
    // std::cout << "\td Time: " << getTime() << "\tNodes:" << getNodes()
	   //      << std::endl;
    // if(_verbosity > 1){
    //     SCIP_CALL_EXC(SCIPprintStatistics(_scip, NULL));
    // }
}

int GurobiSolver::getNodes(){
    return (int)(model->get(GRB_DoubleAttr_NodeCount));
}

double GurobiSolver::getTime(){
    return model->get(GRB_DoubleAttr_Runtime);
}

double GurobiSolver::get_value(void *ptr){
    int index = *(int*)(ptr);
    if(index >= 0 and index < variables->size()){
        GRBVar v = variables->at(index);
        return v.get(GRB_DoubleAttr_X);
    }
    return 0.0;
}


// Simple function to determine the MIP gap
// double gap(GRBModel *model) throw(GRBException)
// {
//   if ((model->get(GRB_IntAttr_SolCount) == 0) ||
//       (fabs(model->get(GRB_DoubleAttr_ObjVal)) < 1e-6))
//   {
//     return GRB_INFINITY;
//   }
//   return fabs(model->get(GRB_DoubleAttr_ObjBound) -
//               model->get(GRB_DoubleAttr_ObjVal)) /
//          fabs(model->get(GRB_DoubleAttr_ObjVal));
// }
