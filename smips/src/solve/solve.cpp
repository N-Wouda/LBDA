#include "solve.h"


std::unique_ptr<arma::vec> solve(MasterProblem &master, Cut &cut, double tol)
{
    while (true)
    {
        auto sol = master.solve();
        auto cutResult = cut.computeCut(*sol.x);

        if (MasterProblem::isValidCut(cutResult, sol, tol))
            master.addCut(cutResult);
        else
            return std::move(sol.x);
    }
}
