#include "cuts/strongbenders.h"
#include "subproblem.h"


StrongBenders::CutResult StrongBenders::computeCut(arma::vec const &x)
{
    auto const &Tmat = d_problem.Tmat();

    arma::vec Tx = Tmat * x;
    arma::vec beta = arma::zeros(Tmat.n_rows);

    auto sub = SubProblem(d_env, d_problem);

    double gamma = 0;
    double lw = 0;

    for (size_t s = 0; s != d_problem.nScenarios(); ++s)
    {
        arma::vec omega = d_problem.scenarios().col(s);

        sub.update(omega - Tx);
        sub.solve();

        auto const info = sub.multipliers();
        double const prob = d_problem.d_probabilities[s];

        arma::vec pi = Tmat * info.lambda;

        d_lr.update(omega, pi);
        gamma += prob * d_lr.solve();

        beta -= prob * pi;

        lw += prob * arma::dot(info.lambda, omega);
    }

    std::cout << "lw = " << lw << ". L(pi, w) = " << gamma << '\n';
    return CutResult{beta, gamma};
}
