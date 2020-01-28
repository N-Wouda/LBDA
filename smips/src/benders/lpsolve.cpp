#include "benders.h"

#include <chrono>

std::unique_ptr<arma::vec> Benders::lpSolve(double tol)
{
    using seconds = std::chrono::seconds;
    using clock = std::chrono::high_resolution_clock;

    auto t1 = clock::now();
    size_t iter = 1;

    while (true)
    {
        ++iter;

        auto sol = d_master.solve();

        // derive cut
        arma::vec beta(d_problem.d_n1, arma::fill::zeros);
        double gamma = 0;
        lpCut(*sol.x, beta, gamma);

        // add the cut (conditional on it being violated by the current
        // solution)
        if (d_master.addCut(*sol.x, beta, gamma, sol.theta, tol))
        {
            auto t2 = clock::now();

            d_runTime += std::chrono::duration_cast<seconds>(t2 - t1).count();
            d_nCuts += iter;

            return std::move(sol.x);
        }
    }
}
