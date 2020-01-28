#ifndef BENDERS_H
#define BENDERS_H

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
    Lagrangian d_lr;     // lagrangian relaxation
    Gomory d_gomory;     // Gomory relaxation

    // for each scenario, we store the basis matrices that we
    // have visited (encoded by vBasis, cBasis)
    std::vector<std::vector<std::vector<double>>> d_visited;

    // for each visited basis matrix, we store the
    // corresponding gomory objective value
    std::vector<std::vector<double>> d_objectives;

    double computeGomory(size_t s, int *vBasis, int *cBasis, arma::vec &rhs);

    void computeTx(arma::vec const &x, arma::vec &Tx);

    void lbdaCut(arma::vec const &x,
                 arma::vec const &alpha,
                 arma::vec &beta,
                 double &gamma);

    void lpCut(arma::vec const &x, arma::vec &beta, double &gamma);

    void sbCut(arma::vec const &x, arma::vec &beta, double &gamma);

public:
    Benders(GRBEnv &env, GRBenv *c_env, Problem &problem);

    Benders(Benders const &other);

    std::unique_ptr<arma::vec> lpSolve(double tol = 1e-4);

    std::unique_ptr<arma::vec> strongBenders(double tol = 1e-4);

    std::unique_ptr<arma::vec> lbda(arma::vec const &alpha,
                                    double timeLimit = 1e6,
                                    double tol = 1e-4);

    // TODO this is the dream.
    //std::unique_ptr<arma::vec> solve(Cut &cut, double tol = 1e-4);
};

#endif
