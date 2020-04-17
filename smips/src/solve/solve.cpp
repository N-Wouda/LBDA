#include "solve.h"


std::unique_ptr<arma::vec> solve(MasterProblem &master,
                                 Decomposition &decomposition,
                                 double tol)
{
    while (true)
    {
        auto sol = master.solve();
        auto cut = decomposition.computeCut(*sol.x);

        if (MasterProblem::isValidCut(cut, sol, tol))
            master.addCut(cut);
        else
            return std::move(sol.x);
    }
}
