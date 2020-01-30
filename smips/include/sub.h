#ifndef SUB_H
#define SUB_H

#include "problem.h"

#include <armadillo>
#include <gurobi_c++.h>
#include <iosfwd>

class Sub
{
    size_t d_m2;
    size_t d_n2;

    arma::vec const &d_q;

    GRBModel d_model;
    GRBConstr *d_constrs;
    GRBVar *d_vars;

public:
    Sub(GRBEnv &env, Problem const &problem);

    Sub(Sub const &other);

    ~Sub();

    void update(arma::vec &rhs);

    void solve();

    struct GomInfo  // TODO maybe combine with Multipliers?
    {
        arma::vec lambda;
        arma::Col<int> vBasis;
        arma::Col<int> cBasis;
    };

    struct Multipliers
    {
        arma::vec lambda;
        arma::vec pi_u;
    };

    Multipliers const multipliers();

    GomInfo const gomInfo();
};

#endif
