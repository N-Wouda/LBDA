#include "master.h"


bool Master::isValidCut(Cut::CutResult const &cutResult,
                        Master::Solution const &sol,
                        double tol)
{
    return cutResult.gamma + arma::dot(*sol.x, cutResult.beta)
           >= sol.theta + tol;
}
