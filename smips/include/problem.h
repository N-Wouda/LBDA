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

    void initSub();  // initializes the subproblem, and sets rhs = 0. Called by
                     // evaluate() when evaluate is called for the first time.

    void clearSub();  // should be called if problem data changes

public:
    // TODO make these members private
    // size parameters
    size_t d_m1;  // number of rows of A
    size_t d_m2;  // number of rows of W (and T)
    size_t d_n1;  // number of columns of A (and T)
    size_t d_n2;  // number of columns of W
    size_t d_p1;  // number of integer first-stage decision variables
    size_t d_p2;  // idem (second stage)
    size_t d_q1;  // TODO unused?
    size_t d_q2;  // TODO unused?
    size_t d_S;

    double d_L;  // lb of Q

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

    arma::mat d_Amat;
    arma::mat d_Tmat;
    arma::mat d_Wmat;

    // rows of d_omega correspond to scenarios
    arma::mat d_omega;
    arma::vec d_probs;

    Problem(Data &generator, GRBEnv &env);

    Problem(const Problem &other) = delete;

    ~Problem();

    // Initializes A, T, W, b, c, q (>= constraints)
    void randomInstance(int A_low = 1,
                        int A_high = 6,
                        int T_low = 1,
                        int T_high = 6,
                        int W_low = 1,
                        int W_high = 6,
                        int c_low = 1,
                        int c_high = 5,
                        int b_low = 1,
                        int b_high = 5,
                        int q_low = 5,
                        int q_high = 10);

    void enforceCcr(double penalty);  // TODO: enforce CCR assumption

    // initializes d_omega
    void setGaussianOmega(double mean, double sd);

    void setBounds(arma::vec &l1, arma::vec &u1, arma::vec &l2, arma::vec &u2);

    void sizes(size_t S);

    void ssv95(size_t S,
               bool fs_continuous,
               bool ss_binary,
               bool standard_T = true);

    void sslp(size_t nServers, size_t nClients, size_t S);

    void classic_ri();

    // evaluates cx + Q(x) (does not check feasibility)
    double evaluate(arma::vec const &x);
};

#endif