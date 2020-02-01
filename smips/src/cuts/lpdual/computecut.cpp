#include "cuts/lpdual.h"
#include "subproblem.h"


LpDual::CutResult LpDual::computeCut(arma::vec const &x)
{
    arma::vec Tx(d_problem.d_m2);
    computeTx(x, Tx);

    arma::vec beta = arma::zeros(d_problem.d_n1);
    double gamma = 0;

    arma::vec dual = arma::zeros(d_problem.d_m2);  // cut coefficients

    auto sub = SubProblem(d_env, d_problem);

    for (size_t s = 0; s != d_problem.d_S; ++s)
    {
        arma::vec ws(d_problem.d_omega[s]);
        arma::vec rhs = ws - Tx;

        sub.update(rhs);
        sub.solve();

        auto const info = sub.multipliers();
        double const prob = d_problem.d_probs[s];

        dual -= prob * info.lambda;
        gamma += prob * arma::dot(info.lambda, ws);
        gamma += prob * arma::dot(info.pi_u, d_problem.d_u2);
    }

    for (size_t col = 0; col != d_problem.d_n1; ++col)
    {
        for (size_t row = 0; row != d_problem.d_m2; ++row)
            beta[col] += dual[row] * d_problem.d_Tmat[row][col];
    }

    return CutResult{beta, gamma};
}
