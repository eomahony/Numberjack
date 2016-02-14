#include "CPLEX.hpp"
ILOSTLBEGIN

CPLEXSolver::CPLEXSolver() {
    DBG("Create a CPLEX solver %s\n", "");

    var_counter = 0;
    _verbosity = 0;
    cplextime = 0.0;
    has_been_added = false;
    optimstatus = IloAlgorithm::Unknown;
    env = new IloEnv();
    model = new IloModel(*env);
    cplex = new IloCplex(*model);
    variables = new IloNumVarArray(*env);
    setVerbosity(_verbosity);  // Default to no verbose output
}

CPLEXSolver::~CPLEXSolver() {
    DBG("Delete wrapped solver%s\n", "");

    for(int i=0; i<_constraints.size(); i++) delete _constraints[i];
    _constraints.clear();

    for(int i=0; i<variableptrs.size(); i++) delete variableptrs[i];
    variableptrs.clear();

    env->end();
    delete variables;
    delete cplex;
    delete model;
    delete env;
}

void CPLEXSolver::initialise(MipWrapperExpArray &arg) {
    DBG("Initialising solver with array of expressions%s\n", "");
}

void CPLEXSolver::initialise() {
    DBG("Initialising solver with no expressions%s\n", "");
    MipWrapperSolver::initialise();

    if(_obj != NULL){
        add_in_constraint(_obj, _obj_coef);
    }

    for(unsigned int i = 0; i < _constraints.size(); ++i)
        add_in_constraint(_constraints[i]);
    has_been_added = true;
}

void CPLEXSolver::add_in_constraint(LinearConstraint *con, double coef){
    DBG("Creating a CPLEX representation of a constriant %s\n", "");
    IloNumArray weights(*env, (IloInt)con->_coefficients.size());
    IloNumVarArray vars(*env, (IloInt)con->_variables.size());
    
    for(unsigned int i = 0; i < con->_variables.size(); ++i){
        DBG("\tAdding variable to CPLEX\n%s", "");
        IloNumVar var_ptr;
        
        if(con->_variables[i]->_var == NULL){
            
            IloNumVar::Type type;
        
            if(con->_variables[i]->_continuous) type = IloNumVar::Float;
            // else if(con->_lower == 0 && con->_upper == 1) type = IloNumVar::Bool;
            else type = IloNumVar::Int;

            var_ptr = IloNumVar(getEnv(),
                                con->_variables[i]->_lower, // LB
                                con->_variables[i]->_upper, // UB
                                type);

            int *var_id = new int;
            *var_id = variables->getSize();
            variables->add(var_ptr);
            variableptrs.push_back(var_id);
            con->_variables[i]->_var = (void*) var_id;
            
            DBG("Created new variable with id %d. type:%c lb:%f ub:%f coef:%f\n", *var_id, type, con->_variables[i]->_lower, con->_variables[i]->_upper, coef);
        } else {
            var_ptr = (*variables)[*(int*)(con->_variables[i]->_var)];
        }

        vars[i] = (*variables)[*(int*)(con->_variables[i]->_var)];
        weights[i] = con->_coefficients[i];
    }

    IloNumExprArg lin_expr = IloScalProd(weights, vars);

    if(coef < -0.1){
        model->add(IloMinimize(*env, lin_expr));
    } else if(coef > 0.1){
        model->add(IloMaximize(*env, lin_expr));
    } else {
        if(con->_lhs > -INFINITY && con->_rhs < INFINITY){
            if(con->_lhs == con->_rhs) {
                model->add(lin_expr == con->_lhs);
            } else {
                model->add(IloRange(*env, con->_lhs, lin_expr, con->_rhs));
            }
        } else if(con->_lhs > -INFINITY) model->add(lin_expr >= con->_lhs);
        else if(con->_rhs < INFINITY) model->add(lin_expr <= con->_rhs);
    }
}


int CPLEXSolver::solve() {
    if(!has_been_added) initialise();

    double start_time = cplex->getCplexTime();
    cplex->extract(*model);
    cplex->setParam(IloCplex::MIPSearch, IloCplex::Traditional);
    cplex->solve();
    cplextime = (cplex->getCplexTime() - start_time);

    optimstatus = cplex->getStatus();
    return is_sat();
}

void CPLEXSolver::setNodeLimit(const int cutoff) {
    if(cutoff < 0) {
        std::cerr << "Warning: cannot specify a negative node limit, ignoring." << std::endl;
    } else {
        cplex->setParam(IloCplex::NodeLim, cutoff);
    }
}

void CPLEXSolver::setTimeLimit(const int cutoff) {
    if(cutoff < 0) {
        std::cerr << "Warning: cannot specify a negative node limit, ignoring." << std::endl;
    } else {
        cplex->setParam(IloCplex::TiLim, cutoff);
    }
}

void CPLEXSolver::setThreadCount(const int nr_threads){
    if(nr_threads < 0){
        std::cerr << "Warning: cannot specify a negative thread count, ignoring." << std::endl;
    } else {
        cplex->setParam(IloCplex::Threads, nr_threads);
    }
}

void CPLEXSolver::setVerbosity(const int degree) {
    // http://pic.dhe.ibm.com/infocenter/cosinfoc/v12r2/topic/ilog.odms.cplex.help/Content/Optimization/Documentation/CPLEX/_pubskel/CPLEX999.html
    if(degree >= 0 && degree <= 5)
        cplex->setParam(IloCplex::MIPDisplay, degree);
}

void CPLEXSolver::setRandomSeed(const int seed) {
    cplex->setParam(IloCplex::RandomSeed, seed);
}

void CPLEXSolver::setWorkMem(const int mb){
    cplex->setParam(IloCplex::WorkMem, mb);
}

int CPLEXSolver::getWorkMem(){
    return cplex->getParam(IloCplex::WorkMem);
}

bool CPLEXSolver::is_opt() {
    return optimstatus == IloAlgorithm::Optimal;
}

bool CPLEXSolver::is_sat() {
    return optimstatus == IloAlgorithm::Optimal || optimstatus == IloAlgorithm::Feasible;
}

bool CPLEXSolver::is_unsat() {
    return optimstatus == IloAlgorithm::Infeasible;
}

void CPLEXSolver::output_lp(const char *filename) {
    cplex->exportModel(filename);
}

void CPLEXSolver::output_mps(const char *filename) {
    cplex->exportModel(filename);
}

void CPLEXSolver::printStatistics() {
    std::cout << "\td Time: " << getTime() << "\tNodes:" << getNodes() << std::endl;
}

int CPLEXSolver::getNodes() {
    return cplex->getNnodes();
}

double CPLEXSolver::getTime() {
    return cplextime;
}

double CPLEXSolver::getOptimalityGap(){
    return cplex->getMIPRelativeGap();
}

double CPLEXSolver::get_value(void *ptr){
    int index = *(int*)(ptr);
    if(index >= 0 && index < variables->getSize()){
        IloNumVar v = (*variables)[index];
        return cplex->getValue(v, -1);
    }
    return 0.0;
}
