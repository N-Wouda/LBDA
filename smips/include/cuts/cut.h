#ifndef CUT_H
#define CUT_H

#include "problem.h"

#include <armadillo>
#include <gurobi_c++.h>


class Cut
{
protected:
    GRBEnv &d_env;

    Problem const &d_problem;

    struct CutResult
    {
        arma::vec beta;
        double gamma;
    };

    explicit Cut(GRBEnv &env, Problem const &problem);

    void computeTx(arma::vec const &x, arma::vec &Tx);

public:
    virtual CutResult computeCut(arma::vec const &x) = 0;
};

#endif  // CUT_H
