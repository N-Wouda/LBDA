#ifndef STRONGBENDERS_H
#define STRONGBENDERS_H

#include "cut.h"
#include "problem.h"

#include <armadillo>
#include <gurobi_c++.h>

class StrongBenders : public Cut
{
    GRBConstr *d_constrs;
    GRBVar *d_z_vars;

    void update(arma::vec &rhs, arma::vec &pi);

public:
    StrongBenders(GRBEnv &env, Problem const &problem);

    ~StrongBenders();

    Cut::CutResult computeCut(arma::vec const &x) override;
};

#endif  // STRONGBENDERS_H
