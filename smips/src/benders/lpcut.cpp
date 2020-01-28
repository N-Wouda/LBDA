#include "benders.h"

void Benders::lpCut(arma::vec const &x, arma::vec &beta, double &gamma)
{
    arma::vec Tx(d_problem.d_m2);
    computeTx(x, Tx);

    arma::vec dual = arma::zeros(d_problem.d_m2);  // cut coefficients

    auto sub = Sub(d_env, d_problem);

    for (size_t s = 0; s != d_problem.d_S; ++s)
    {
        arma::vec ws(d_problem.d_omega[s]);
        arma::vec rhs = ws - Tx;

        sub.update(rhs);
        sub.solve();

        auto const info = sub.multipliers();

        double *lambda = info.lambda;
        double *pi_u = info.pi_u;
        double const prob = d_problem.d_probs[s];

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
        for (size_t row = 0; row != d_problem.d_m2; ++row)
            beta[col] += dual[row] * d_problem.d_Tmat[row][col];
    }
}