#include "benders.h"

void Benders::sbCut(arma::vec const &x, arma::vec &beta, double &gamma)
{
    arma::vec Tx(d_problem.d_m2);
    computeTx(x, Tx);

    double lw = 0;

    auto sub = Sub(d_env, d_problem);

    for (size_t s = 0; s != d_problem.d_S; ++s)
    {
        arma::vec ws(d_problem.d_omega[s]);
        double const prob = d_problem.d_probs[s];

        arma::vec rhs = ws - Tx;

        sub.update(rhs);
        sub.solve();

        auto const info = sub.multipliers();

        arma::vec pi(d_problem.d_n1, arma::fill::zeros);

        for (size_t var = 0; var != d_problem.d_n1; ++var)
        {
            for (size_t row = 0; row != d_problem.d_m2; ++row)
                pi[var] += info.lambda[row] * d_problem.d_Tmat[row][var];

            beta[var] -= prob * pi[var];
        }

        d_lr.update(ws, pi);
        gamma += prob * d_lr.solve();

        for (size_t row = 0; row != d_problem.d_m2; ++row)
            lw += prob * info.lambda[row] * ws[row];

        delete[] info.lambda;
        delete[] info.pi_u;
    }

    std::cout << "lw = " << lw << ". L(pi, w) = " << gamma << '\n';
}