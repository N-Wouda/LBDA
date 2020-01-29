#ifndef BENDERS_H
#define BENDERS_H

#include "cuts/cut.h"
#include "decompositions/cglp.h"
#include "decompositions/gomory.h"
#include "decompositions/lagrangian.h"
#include "master.h"
#include "problem.h"
#include "sub.h"

#include <armadillo>

class Benders
{
    size_t d_nCuts;
    double d_runTime;

    GRBEnv &d_env;

    Problem &d_problem;  // contains problem data
    Master d_master;     // master problem

public:
    Benders(GRBEnv &env, GRBenv *c_env, Problem &problem);

    Benders(Benders const &other);

    std::unique_ptr<arma::vec> solve(Cut &cut, double tol = 1e-4);
};

#endif
