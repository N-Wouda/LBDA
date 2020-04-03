#ifndef PROBLEM_H
#define PROBLEM_H

#include "data.h"

#include <armadillo>
#include <gurobi_c++.h>
#include <iosfwd>


class Problem
{
    Data d_gen;  // used to generate (random) data

    // TODO this is not the place - should move to the actual solver stuff.
    GRBModel d_sub;  // subproblem, used for evaluation of cx + Q(x) (more
                     // precisely, of v(omega, x)).

    // useful to quickly update rhs of d_sub (heap allocated)
    GRBConstr *d_constrs;

    // is d_sub initialized. (destructor will only call delete[] on d_constrs
    // if yes)
    bool d_sub_initialized;

    size_t d_nFirstStageIntVars = 0;    // TODO make constant!
    size_t d_nSecondStageIntVars = 0;

    arma::mat d_Amat;
    arma::mat d_Tmat;
    arma::mat d_Wmat;

    // Each column in d_omegas correspond to a single scenario (omega).
    arma::mat d_omegas;

    void initSub();  // initializes the subproblem, and sets rhs = 0. Called by
                     // evaluate() when evaluate is called for the first time.

    void clearSub();  // should be called if problem data changes

public:
    // TODO make these members private
    double d_L;  // lb of Q - TODO do we need this?

    // number of >= and <= constraints in the first and second stage
    size_t d_fs_leq;
    size_t d_fs_geq;

    size_t d_ss_leq;
    size_t d_ss_geq;

    arma::vec d_l1;
    arma::vec d_u1;
    arma::vec d_l2;
    arma::vec d_u2;

    arma::vec d_c;
    arma::vec d_b;
    arma::vec d_q;

    arma::vec d_probabilities;

    Problem(Data &generator, GRBEnv &env);

    Problem(const Problem &other) = delete;

    ~Problem();

    void enforceCcr(double penalty);  // TODO: enforce CCR assumption

    void ssv95(size_t S,
               bool fs_continuous,
               bool ss_binary,
               bool standard_T = true);

    // evaluates cx + Q(x) (does not check feasibility)
    double evaluate(arma::vec const &x);

    [[nodiscard]] size_t nFirstStageIntVars() const;

    [[nodiscard]] size_t nSecondStageIntVars() const;

    [[nodiscard]] size_t nScenarios() const;

    [[nodiscard]] bool isMixedIntegerProblem() const;

    [[nodiscard]] arma::mat const &Amat() const;
    [[nodiscard]] arma::mat &Amat();

    [[nodiscard]] arma::mat const &Wmat() const;

    [[nodiscard]] arma::mat const &Tmat() const;

    [[nodiscard]] arma::mat const &scenarios() const;
};

inline size_t Problem::nFirstStageIntVars() const
{
    return d_nFirstStageIntVars;
};

inline size_t Problem::nSecondStageIntVars() const
{
    return d_nSecondStageIntVars;
};

inline size_t Problem::nScenarios() const
{
    return d_omegas.n_cols;
}

inline bool Problem::isMixedIntegerProblem() const
{
    // TODO should this not be or? Either makes the problem an integer problem,
    // right?
    return nFirstStageIntVars() != 0 and nSecondStageIntVars() != 0;
};

inline arma::mat const &Problem::Amat() const
{
    return d_Amat;
}

inline arma::mat &Problem::Amat()
{
    return d_Amat;
}

inline arma::mat const &Problem::Wmat() const
{
    return d_Wmat;
}

inline arma::mat const &Problem::Tmat() const
{
    return d_Tmat;
}

inline arma::mat const &Problem::scenarios() const
{
    return d_omegas;
}

#endif