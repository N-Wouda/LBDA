#include "benders.h"

void Benders::computeTx(arma::vec const &x, arma::vec &Tx)
{
    auto &Tmat = d_problem.d_Tmat;

    for (size_t zvar = 0; zvar != d_problem.d_m2; ++zvar)
    {
        Tx[zvar] = 0.0;
        for (size_t xvar = 0; xvar != d_problem.d_n1; ++xvar)
            Tx[zvar] += Tmat[zvar][xvar] * x[xvar];
    }
}
