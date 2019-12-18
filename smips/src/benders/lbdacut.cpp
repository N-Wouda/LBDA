#include "benders.h"

double Benders::lbdaCut(double *x, double *beta, double *alpha)
{
    auto &omega = d_problem.d_omega;
    auto &Tmat = d_problem.d_Tmat;
    auto &probs = d_problem.d_probs;

    double Tx[d_m2];
    computeTx(x, Tx);  // Tx is rba

    // cut coefficients: initialize to zero
    double gamma = 0.0;
    double dual[d_m2];
    std::fill(dual, dual + d_m2, 0.0);

    for (size_t s = 0; s != d_S; ++s)
    {
        double *ws = omega[s].data();  // scenario (c-style array pointer)
                                       // compute rhs, update subproblem
        double rhs[d_m2];  // rhs vector of subproblem (c-style array)

        for (size_t row = 0; row != d_m2; ++row)  // compute element-by-element
            rhs[row] = ws[row] - Tx[row];

        d_sub.update(rhs);
        Sub::GomInfo info = d_sub.solve2();  // solve subproblem

        double *lambda = info.lambda;  // extract lambda (for optimality cut)

        // extract vBasis (to update gomory relaxation)
        int *vBasis = info.vBasis;

        // extract vBasis (to update gomory relaxation)
        int *cBasis = info.cBasis;

        double gom_obj = compute_gomory(s, vBasis, cBasis, ws, alpha);
        double prob = probs[s];


        gamma += prob * gom_obj;  // gom_obj = lambda^T (omega - alpha) +
                                  // psi(omega - alpha), thus, we add lambda^T
                                  // alpha in the following loop

        for (size_t row = 0; row != d_m2; ++row)
        {
            dual[row] -= prob * lambda[row];
            gamma += prob * lambda[row] * alpha[row];
        }

        delete[] lambda;
        delete[] vBasis;
        delete[] cBasis;
    }

    for (size_t col = 0; col != d_n1; ++col)
    {
        beta[col] = 0.0;

        for (size_t row = 0; row != d_m2; ++row)
            beta[col] += dual[row] * Tmat[row][col];
    }

    return gamma;
}
