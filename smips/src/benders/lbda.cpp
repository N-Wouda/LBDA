#include "benders.h"

#include <chrono>


void Benders::lbda(double *alpha, double gomoryTimeLimit, double tol)
{
    using clock = std::chrono::high_resolution_clock;
    using seconds = std::chrono::seconds;

    d_gomory.setTimeLimit(gomoryTimeLimit);

    auto start = clock::now();

    bool stop = false;
    size_t iterations = 0;

    while (not stop)
    {
        ++iterations;

        // solve master problem, and collect x and theta
        Master::Solution sol = d_master.solve();

        double *x = sol.xVals;
        double theta = sol.thetaVal;
        double beta[d_n1];

        // derive cut
        double gamma = lbdaCut(x, beta, alpha);  // beta is RBA

        // add the cut (conditional on it being violated by the current
        // solution)
        stop = d_master.addCut(beta, gamma, x, theta, tol);

        if (stop)
            std::copy(x, x + d_n1, d_xvals);

        delete[] x;
    }

    auto final = clock::now();

    d_runTime += std::chrono::duration_cast<seconds>(final - start).count();
    d_nCuts += iterations - 1;
}
