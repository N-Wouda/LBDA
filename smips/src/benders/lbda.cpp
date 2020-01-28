#include "benders.h"

#include <chrono>


std::unique_ptr<arma::vec> Benders::lbda(double *alpha,
                                         double timeLimit,
                                         double tol)
{
    using clock = std::chrono::high_resolution_clock;
    using seconds = std::chrono::seconds;

    d_gomory.setTimeLimit(timeLimit);

    auto t1 = clock::now();

    size_t iterations = 1;

    while (true)
    {
        ++iterations;

        auto sol = d_master.solve();

        // derive cut
        arma::vec beta(d_problem.d_n1, arma::fill::zeros);
        double gamma = 0;
        lbdaCut(*sol.x, alpha, beta, gamma);

        // add the cut (conditional on it being violated by the current
        // solution)
        if (d_master.addCut(*sol.x, beta, gamma, sol.theta, tol))
        {
            auto t2 = clock::now();

            d_runTime += std::chrono::duration_cast<seconds>(t2 - t1).count();
            d_nCuts += iterations;

            return std::move(sol.x);
        }
    }
}
