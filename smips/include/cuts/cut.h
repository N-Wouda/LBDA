#ifndef CUT_H
#define CUT_H

#include "problem.h"

#include <armadillo>
#include <gurobi_c++.h>


class Cut
{
protected:
    GRBEnv &d_env;
    GRBModel d_model;
    Problem const &d_problem;

    explicit Cut(GRBEnv &env, Problem const &problem);

public:
    struct CutResult
    {
        arma::vec beta;
        double gamma;
    };

    // TODO make protected
    virtual Cut::CutResult computeCut(arma::vec const &x) = 0;

    virtual double solve();
};

#endif  // CUT_H
