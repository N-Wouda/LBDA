#include "decompositions/lpdual.h"

#include "subproblem.h"

LpDual::LpDual(GRBEnv &env, Problem const &problem) :
    Decomposition(env, problem)
{
}

LpDual::Cut LpDual::computeCut(arma::vec const &x)
{
    auto const &Tmat = d_problem.Tmat();

    arma::vec Tx = (x.t() * Tmat).t();          // TODO simplify
    arma::vec dual = arma::zeros(Tmat.n_cols);  // decomposition coeffs

    auto sub = SubProblem(d_env, d_problem);

    double gamma = 0;

    for (size_t scenario = 0; scenario != d_problem.nScenarios(); ++scenario)
    {
        arma::vec omega = d_problem.scenarios().col(scenario);

        sub.update(omega - Tx);
        sub.solve();

        auto const info = sub.multipliers();
        double const prob = d_problem.probability(scenario);

        gamma += prob * arma::dot(info.lambda, omega);

        // TODO discuss with NvdL (what happens when unconstrained?)
        if (d_problem.d_secondStageUpperBound.is_finite())
        {
            auto val = arma::dot(info.pi_u, d_problem.d_secondStageUpperBound);
            gamma += prob * val;
        }

        dual -= prob * info.lambda;
    }

    return Cut{Tmat * dual, gamma};
}
