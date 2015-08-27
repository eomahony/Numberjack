
#ifndef GUROBI_H
#define GUROBI_H

#include <iostream>
#include <math.h>
#include <vector>
#include <gurobi_c++.h>
#include "MipWrapper.hpp"

/**
     The solver itself
*/
class GurobiSolver : public MipWrapperSolver{
private:

    GRBEnv *env;
    GRBModel *model;
    vector<GRBVar> *variables;
    int _verbosity, optimstatus;
    bool has_been_added;
    void add_in_constraint(LinearConstraint *con, double coef=0);
    void add_variables_from(LinearConstraint *con, double coef=0);

public:

    GurobiSolver();
    virtual ~GurobiSolver();
    
    // initialise the solver before solving (no more calls to add after this)
    void initialise(MipWrapperExpArray& arg);
    void initialise();

    // solving methods
    int solve();
    
    // parameter tuning methods
    void setTimeLimit(const int cutoff);
    void setNodeLimit(const int cutoff);
    void setThreadCount(const int nr_threads);
    void setVerbosity(const int degree);

    // statistics methods
    int getNodes();
    bool is_sat();
    bool is_unsat();
    bool is_opt();
    void output_lp(const char *filename);
    void output_mps(const char *filename);
    void printStatistics();
    double getTime();
    double getOptimalityGap();
    
    // Value stuff
    virtual double get_value(void *ptr);

};

#endif
