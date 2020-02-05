#include "cuts/loosebenders.h"
#include "subproblem.h"


LooseBenders::CutResult LooseBenders::computeCut(arma::vec const &x)
{
    arma::vec Tx(d_problem.d_m2);
    computeTx(x, Tx);

    arma::vec beta = arma::zeros(d_problem.d_n1);
    double gamma = 0;

    arma::vec dual = arma::zeros(d_problem.d_m2);  // cut coefficients

    auto sub = SubProblem(d_env, d_problem);

    for (size_t s = 0; s != d_problem.d_S; ++s)
    {
        arma::vec omega(d_problem.d_omega[s]);

        sub.update(omega - Tx);
        sub.solve();

        auto const info = sub.gomInfo();
        double const prob = d_problem.d_probs[s];

        // Gomory is lambda^T (omega - alpha) + psi(omega - alpha), so we add
        // lambda^T alpha.
        arma::vec rhs = omega - d_alpha;
        gamma += prob * computeGomory(s, rhs, info.vBasis, info.cBasis);
        gamma += prob * arma::dot(info.lambda, d_alpha);

        dual -= prob * info.lambda;
    }

    for (size_t col = 0; col != d_problem.d_n1; ++col)
    {
        for (size_t row = 0; row != d_problem.d_m2; ++row)
            beta[col] += dual[row] * d_problem.d_Tmat[row][col];
    }

    return CutResult{beta, gamma};
}
