#include "benders.h"

void Benders::lpCut(double *x, double *beta, double &gamma)
{
    auto &omega = d_problem.d_omega;
    auto &Tmat = d_problem.d_Tmat;
    auto &probs = d_problem.d_probs;

    gamma = 0;

    arma::vec Tx(d_problem.d_m2);
    computeTx(x, Tx);

    double dual[d_problem.d_m2];
    std::fill(dual, dual + d_problem.d_m2, 0.0);

    auto sub = Sub(d_env, d_problem);

    for (size_t s = 0; s != d_problem.d_S; ++s)
    {
        arma::vec ws(omega[s]);
        arma::vec rhs = ws - Tx;

        sub.update(rhs);
        sub.solve();

        auto const info = sub.multipliers();

        double *lambda = info.lambda;
        double *pi_u = info.pi_u;
        double prob = probs[s];

        for (size_t row = 0; row != d_problem.d_m2; ++row)
        {
            dual[row] -= prob * lambda[row];
            gamma += prob * lambda[row] * ws[row];
        }

        for (size_t var = 0; var != d_problem.d_n2; ++var)
            gamma += prob * pi_u[var] * d_problem.d_u2[var];

        delete[] lambda;
        delete[] pi_u;
    }

    for (size_t col = 0; col != d_problem.d_n1; ++col)
    {
        beta[col] = 0.0;
        for (size_t row = 0; row != d_problem.d_m2; ++row)
            beta[col] += dual[row] * Tmat[row][col];
    }
}