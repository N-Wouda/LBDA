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
        double *x = sol.x->memptr();

        // derive cut
        double beta[d_problem.d_n1];
        double gamma;
        lpCut(x, beta, gamma);

        // add the cut (conditional on it being violated by the current
        // solution)
        if (d_master.addCut(x, beta, gamma, sol.theta, tol))
        {
            auto t2 = clock::now();

            d_runTime += std::chrono::duration_cast<seconds>(t2 - t1).count();
            d_nCuts += iter;

            return std::move(sol.x);
        }
    }
}
