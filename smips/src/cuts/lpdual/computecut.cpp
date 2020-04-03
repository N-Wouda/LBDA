#include "cuts/lpdual.h"
#include "subproblem.h"


LpDual::CutResult LpDual::computeCut(arma::vec const &x)
{
    auto const &Tmat = d_problem.Tmat();

    arma::vec Tx = Tmat * x;
    arma::vec dual = arma::zeros(Tmat.n_cols);  // cut coeffs

    auto sub = SubProblem(d_env, d_problem);

    double gamma = 0;

    for (size_t scenario = 0; scenario != d_problem.nScenarios(); ++scenario)
    {
        arma::vec omega = d_problem.scenarios().col(scenario);

        sub.update(omega - Tx);
        sub.solve();

        auto const info = sub.multipliers();
        double const prob = d_problem.d_probabilities[scenario];

        gamma += prob * arma::dot(info.lambda, omega);
        gamma += prob * arma::dot(info.pi_u, d_problem.d_u2);

        dual -= prob * info.lambda;
    }

    return CutResult{Tmat * dual, gamma};
}
