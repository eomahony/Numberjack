
#ifndef _CPLEX_H
#define _CPLEX_H

#include <iostream>
#include <vector>
#include <ilcplex/ilocplex.h>
#include "MipWrapper.hpp"

class CPLEXSolver : public MipWrapperSolver {
private:
    IloEnv *env;
    IloModel *model;
    IloCplex *cplex;
    IloNumVarArray *variables;
    IloAlgorithm::Status optimstatus;
    std::vector<int*> variableptrs;
    int _verbosity, var_counter;
    double cplextime;
    bool has_been_added;
    void add_in_constraint(LinearConstraint *con, double coef=0);

public:

    CPLEXSolver();
    virtual ~CPLEXSolver();

    inline IloEnv getEnv() {return *env;}

    void initialise(MipWrapperExpArray &arg); // used to initialise search on a given subset of variables
    void initialise();

    int solve();

    void setNodeLimit(const int cutoff);
    void setTimeLimit(const int cutoff);
    void setThreadCount(const int nr_threads);
    void setVerbosity(const int degree);
    void setRandomSeed(const int seed);
    void setWorkMem(const int mb);
    int getWorkMem();

    bool is_opt();
    bool is_sat();
    bool is_unsat();
    void output_lp(const char *filename);
    void output_mps(const char *filename);
    void printStatistics();
    int getNodes();
    double getTime();
    double getOptimalityGap();
    virtual double get_value(void *ptr);
};

#endif // _CPLEX_H
