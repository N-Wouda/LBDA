#ifndef MASTERPROBLEM_H
#define MASTERPROBLEM_H

#include "decompositions/decomposition.h"
#include "problem.h"

#include <armadillo>
#include <gurobi_c++.h>
#include <iosfwd>
#include <memory>

class MasterProblem
{
    /**
     * Reference to the problem instance. Cannot be declared const, as
     * Gurobi expects non-const pointers to some of the Problem data.
     */
    Problem &d_problem;

    /**
     * Pointer to the Gurobi (C) model. The C API gives access to advanced
     * simplex routines. TODO do we need these?
     */
    GRBmodel *d_cmodel;

    /**
     * Number of slack variables in the Gurobi model instance.
     */
    size_t d_nSlacks;

public:
    struct Solution
    {
        std::unique_ptr<arma::vec> x;
        double theta;
    };

    MasterProblem(GRBenv *c_env, Problem &problem);

    MasterProblem(MasterProblem const &other);

    ~MasterProblem();

    /**
     * Determines if the proposed decomposition is violated by the current solution.
     */
    static bool isValidCut(Decomposition::Cut const &cut,
                           MasterProblem::Solution const &sol,
                           double tol);

    /**
     * Adds decomposition <code>theta >= beta^T x + gamma</code>.
     */
    void addCut(Decomposition::Cut &cut);

    /**
     * Solves the current instance of the master problem. This method returns
     * the first-stage decisions (x), and an estimate for the expected second-
     * stage costs (Q(x); theta).
     */
    Solution const solve();
};

#endif
