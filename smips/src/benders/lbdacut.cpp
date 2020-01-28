#include "benders.h"

void Benders::lbdaCut(arma::vec const &x,
                      double *alpha,
                      arma::vec &beta,
                      double &gamma)
{
    arma::vec Tx(d_problem.d_m2);
    computeTx(x, Tx);

    // cut coefficients: initialize to zero
    double dual[d_problem.d_m2];
    std::fill(dual, dual + d_problem.d_m2, 0.0);

    auto sub = Sub(d_env, d_problem);

    for (size_t s = 0; s != d_problem.d_S; ++s)
    {
        arma::vec ws(d_problem.d_omega[s]);
        arma::vec rhs = ws - Tx;

        sub.update(rhs);
        sub.solve();

        auto const info = sub.gomInfo();

        double *lambda = info.lambda;  // lambda (for optimality cut)
        int *vBasis = info.vBasis;     // vBasis (to update gomory relaxation)
        int *cBasis = info.cBasis;     // cBasis (to update gomory relaxation)

        double const gom_obj = computeGomory(s,
                                             vBasis,
                                             cBasis,
                                             ws.memptr(),
                                             alpha);

        double const prob = d_problem.d_probs[s];

        gamma += prob * gom_obj;  // gom_obj = lambda^T (omega - alpha) +
                                  // psi(omega - alpha), thus, we add lambda^T
                                  // alpha in the following loop

        for (size_t row = 0; row != d_problem.d_m2; ++row)
        {
            dual[row] -= prob * lambda[row];
            gamma += prob * lambda[row] * alpha[row];
        }

        delete[] lambda;
        delete[] vBasis;
        delete[] cBasis;
    }

    for (size_t col = 0; col != d_problem.d_n1; ++col)
    {
        for (size_t row = 0; row != d_problem.d_m2; ++row)
            beta[col] += dual[row] * d_problem.d_Tmat[row][col];
    }
}
