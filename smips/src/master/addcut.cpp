#include "master.h"

bool Master::addCut(arma::vec const &x,
                    arma::vec const &beta,
                    double gamma,
                    double theta,
                    double tol)
{
    double betaxgamma = gamma + arma::dot(x, beta);

    if (betaxgamma <= theta + tol)  // betaxgamma >= theta, no cut added
        return true;

    ++d_nSlacks;

    // slack
    GRBaddvar(d_cmodel, 0, nullptr, nullptr, 0, 0, 1e20, GRB_CONTINUOUS, nullptr);

    // slack variable index (there are d_n1 + 1 + nSlacks variables in the
    // Gurobi model)
    size_t slackIdx = d_n1 + d_nSlacks;

    int cind[d_n1 + 2];
    std::iota(cind, cind + d_n1 + 1, 0);

    // refers to the last variable (i.e. the slack)
    cind[d_n1 + 1] = slackIdx;
    double cval[d_n1 + 1];

    for (size_t var = 0; var != d_n1; ++var)
        cval[var + 1] = -beta[var];

    cval[0] = 1;
    cval[d_n1 + 1] = -1;  // >= constraint, so slack features with -1

    GRBaddconstr(d_cmodel, d_n1 + 2, cind, cval, GRB_EQUAL, gamma, nullptr);

    // add cut to internal storage of master
    d_xCoeffs.emplace_back(beta.memptr(), beta.memptr() + d_n1);
    d_cuts.emplace_back(gamma);

    return false;
}