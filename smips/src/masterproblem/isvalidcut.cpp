#include "masterproblem.h"


bool MasterProblem::isValidCut(Decomposition::Cut const &cut,
                               MasterProblem::Solution const &sol,
                               double tol)
{
    return cut.gamma + arma::dot(*sol.x, cut.beta) >= sol.theta + tol;
}
