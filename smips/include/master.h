#ifndef MASTER_H
#define MASTER_H

#include "cuts/cut.h"
#include "problem.h"

#include <armadillo>
#include <gurobi_c++.h>
#include <iosfwd>
#include <memory>

class Master
{
    // TODO is this useful?
    // GRBVar *d_xVars;
    // GRBVar d_theta;
    // GRBModel d_model;
    GRBmodel *d_cmodel;

    // internal storage: only valid for regular l-shaped and lbda cuts
    std::vector<std::vector<double>> d_xCoeffs;
    std::vector<double> d_cuts;

    // slack variable identities slack = kappa * theta - beta * x - gamma
    // TODO is this useful?
    // std::vector<double> d_kappa;
    // std::vector<std::vector<double>> d_beta;
    // std::vector<double> d_gamma;

    size_t d_n1;
    size_t d_nSlacks;

public:
    struct Solution
    {
        std::unique_ptr<arma::vec> x;
        double theta;
    };

    Master(GRBEnv &env, GRBenv *c_env, Problem &problem);

    Master(Master const &other);

    ~Master();

    /**
     * Determines if the proposed cut is violated by the current solution.
     */
    static bool isValidCut(Cut::CutResult const &cutResult,
                    Master::Solution const &sol,
                    double tol);

    /**
     * Adds cut <code>theta >= beta^T x + gamma</code>.
     */
    void addCut(Cut::CutResult &cutResult);

    [[nodiscard]] std::vector<double> const &cuts() const;

    [[nodiscard]] std::vector<std::vector<double>> const &xCoeffs() const;

    [[nodiscard]] size_t n1() const;

    Solution const solve();
};

#endif
