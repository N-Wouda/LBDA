#ifndef DEQFORM_H
#define DEQFORM_H

#include "problem.h"

#include <gurobi_c++.h>
#include <iosfwd>

class DeqForm
{
    enum class status
    {
        SOLVED,
        UNSOLVED,
        INFEASIBLE
    };

public:
    size_t d_n1;
    GRBModel d_model;
    GRBVar *d_xVars;

    status d_status;

    double *d_xVals;

    double d_objVal;
    double d_objBound;

    bool d_MIP;  // is the problem a MIP?
    double d_MIPGap;
    double d_runTime;

    // using the problem data, d_model is initialized
    DeqForm(GRBEnv &env, Problem &problem);

    ~DeqForm();

    void init_fs(size_t n1,
                 size_t p1,
                 size_t m1,
                 size_t fs_leq,
                 size_t fs_geq,
                 double *lb,
                 double *ub,
                 double *c,
                 double *rhs,
                 std::vector<std::vector<double>> &Amat);

    void init_ss(size_t n1,
                 size_t n2,
                 size_t p2,
                 size_t m2,
                 size_t S,
                 size_t ss_leq,
                 size_t ss_geq,
                 double *lb,
                 double *ub,
                 double *probs,
                 double *q,
                 std::vector<std::vector<double>> &Tmat,
                 std::vector<std::vector<double>> &Wmat,
                 std::vector<std::vector<double>> &omega);

    // solves grbModel, updates status, updates xvals, updates MIPGap
    void solve(double time_limit = 1e20);
};

#endif