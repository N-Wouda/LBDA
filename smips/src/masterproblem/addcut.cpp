#include "masterproblem.h"

void MasterProblem::addCut(Cut::CutResult &cutResult)
{
    ++d_nSlacks;

    // slack
    GRBaddvar(d_cmodel,
              0,
              nullptr,
              nullptr,
              0,
              0,
              arma::datum::inf,
              GRB_CONTINUOUS,
              nullptr);

    size_t const n1 = d_problem.Amat().n_rows;

    // slack variable index (there are n1 + 1 + nSlacks variables in the
    // Gurobi model)
    size_t slackIdx = n1 + d_nSlacks;

    int cind[n1 + 2];
    std::iota(cind, cind + n1 + 1, 0);

    // refers to the last variable (i.e. the slack)
    cind[n1 + 1] = slackIdx;
    double cval[n1 + 1];

    for (size_t var = 0; var != n1; ++var)
        cval[var + 1] = -cutResult.beta[var];

    cval[0] = 1;
    cval[n1 + 1] = -1;  // >= constraint, so slack features with -1

    GRBaddconstr(d_cmodel,
                 n1 + 2,
                 cind,
                 cval,
                 GRB_EQUAL,
                 cutResult.gamma,
                 nullptr);

    // add cut to internal storage of master
    // TODO place cuts, not these.
    d_xCoeffs.emplace_back(cutResult.beta.memptr(),
                           cutResult.beta.memptr() + n1);

    d_cuts.emplace_back(cutResult.gamma);
}