#include "benders.h"

void Benders::sbCut(double *x, double *beta, double &gamma)
{
    auto &omega = d_problem.d_omega;
    auto &Tmat = d_problem.d_Tmat;
    auto &probs = d_problem.d_probs;

    gamma = 0;

    double Tx[d_problem.d_m2];
    computeTx(x, Tx);

    double lw = 0;

    auto sub = Sub(d_env, d_problem);

    for (size_t s = 0; s != d_problem.d_S; ++s)
    {
        double *ws = omega[s].data();  // scenario (c-style array pointer)
        double prob = probs[s];

        double rhs[d_problem.d_m2];

        for (size_t row = 0; row != d_problem.d_m2; ++row)
            rhs[row] = ws[row] - Tx[row];

        sub.update(rhs);
        sub.solve();

        auto const info = sub.multipliers();

        double *lambda = info.lambda;
        double *pi_u = info.pi_u;

        double pi[d_problem.d_n1];  // pi = lambda T
        for (size_t var = 0; var != d_problem.d_n1; ++var)
        {
            pi[var] = 0.0;
            for (size_t row = 0; row != d_problem.d_m2; ++row)
                pi[var] += lambda[row] * Tmat[row][var];
            beta[var] += -prob * pi[var];
        }

        d_lr.update(ws, pi);
        double Lpiw = d_lr.solve();
        gamma += prob * Lpiw;

        for (size_t row = 0; row != d_problem.d_m2; ++row)
            lw += prob * lambda[row] * ws[row];

        delete[] lambda;
        delete[] pi_u;
    }

    std::cout << "lw =  " << lw << ". L(pi, w) = " << gamma << '\n';
}