#include "benders.h"

#include <chrono>

std::unique_ptr<arma::vec> Benders::strongBenders(double tol)
{
    using seconds = std::chrono::seconds;
    using clock = std::chrono::high_resolution_clock;

    auto t1 = clock::now();
    size_t iterations = 1;

    while (true)
    {
        iterations++;

        auto sol = d_master.solve();

        // derive cut
        arma::vec beta = arma::zeros(d_problem.d_n1);
        double gamma = 0;
        sbCut(*sol.x, beta, gamma);

        // add the cut conditional on it being violated by the current solution
        if (d_master.addCut(*sol.x, beta, gamma, sol.theta, tol))
        {
            auto t2 = clock::now();

            d_runTime += std::chrono::duration_cast<seconds>(t2 - t1).count();
            d_nCuts += iterations;

            return std::move(sol.x);
        }
    }
}
