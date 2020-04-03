#include "masterproblem.h"


MasterProblem::Solution const MasterProblem::solve()
{
    GRBoptimize(d_cmodel);

    size_t const n1 = d_problem.Amat().n_rows;

    auto xVals = std::make_unique<arma::vec>(n1);
    GRBgetdblattrarray(d_cmodel, "X", 1, n1, xVals->memptr());

    double theta;
    GRBgetdblattrelement(d_cmodel, "X", 0, &theta);

    return Solution{std::move(xVals), theta};
}
