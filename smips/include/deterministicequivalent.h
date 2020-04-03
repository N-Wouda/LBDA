#ifndef DETERMINISTICEQUIVALENT_H
#define DETERMINISTICEQUIVALENT_H

#include "problem.h"

#include <armadillo>
#include <gurobi_c++.h>
#include <iosfwd>
#include <memory>

class DeterministicEquivalent
{
    enum class status
    {
        SOLVED,
        UNSOLVED,
        INFEASIBLE
    };

    Problem const &d_problem;
    GRBModel d_model;

    status d_status;

    double d_objVal = 0;  // TODO use these?
    double d_objBound = 0;
    double d_MIPGap = 0;
    double d_runTime = 0;

    GRBVar *d_xVars;

    void initFirstStage(size_t n1,
                        size_t p1,
                        size_t fs_leq,
                        size_t fs_geq,
                        double const *lb,
                        double const *ub,
                        double const *c,
                        double const *rhs,
                        arma::mat const &Amat);

    void initSecondStage(size_t n1,
                         size_t n2,
                         size_t p2,
                         size_t m2,
                         size_t S,
                         size_t ss_leq,
                         size_t ss_geq,
                         double const *lb,
                         double const *ub,
                         double const *probs,
                         double const *q,
                         arma::mat const &Tmat,
                         arma::mat const &Wmat,
                         arma::mat const &omega);

public:
    DeterministicEquivalent(GRBEnv &env, Problem const &problem);

    ~DeterministicEquivalent();

    std::unique_ptr<arma::vec> solve(double time_limit = 1e20);
};

#endif