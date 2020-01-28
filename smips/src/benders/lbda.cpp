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

        // solve master problem, and collect x and theta
        auto sol = d_master.solve();

        double *x = sol.xVals->memptr();
        double theta = sol.thetaVal;

        // derive cut
        double beta[d_problem.d_n1];
        double gamma;
        lbdaCut(x, alpha, beta, gamma);

        // add the cut (conditional on it being violated by the current
        // solution)
        if (d_master.addCut(x, beta, gamma, theta, tol))
        {
            auto t2 = clock::now();

            d_runTime += std::chrono::duration_cast<seconds>(t2 - t1).count();
            d_nCuts += iterations;

            return std::move(sol.xVals);
        }
    }
}
