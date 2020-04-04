#include "subproblem.h"


SubProblem::Multipliers const SubProblem::multipliers()
{
    // TODO clean this up
    auto const &Wmat = d_problem.Wmat();

    // computing shadow prices of upper bounds of y variables
    double *pi_u;

    // No complete recourse. Dirty fix: ignore the part of the dual objective
    // corresponding to the upper bounds (so we underestimate QLP)
    if (d_model.get(GRB_IntAttr_Status) == GRB_INFEASIBLE)
        pi_u = new double[Wmat.n_rows]();
    else
    {
        // basis information
        int *vbasis = d_model.get(GRB_IntAttr_VBasis, d_vars, Wmat.n_rows);

        // largest value of obj coeff for which current basis remains optimal
        pi_u = d_model.get(GRB_DoubleAttr_SAObjUp, d_vars, Wmat.n_rows);

        for (size_t var = 0; var != Wmat.n_rows; ++var)
        {
            // if variable is not at the upper bound, the shadow price is zero
            if (vbasis[var] != -2)  // TODO magic number
                pi_u[var] = 0;
            else
                // else: equal to the reduced costs
                pi_u[var] = d_problem.d_secondStageCoeffs[var] - pi_u[var];
        }

        delete[] vbasis;
    }

    auto const *lambda = d_model.get(GRB_DoubleAttr_Pi, d_constrs, Wmat.n_cols);

    return Multipliers{arma::vec{lambda, Wmat.n_cols},
                       arma::vec{pi_u, Wmat.n_rows}};
}
