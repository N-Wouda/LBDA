#include "decompositions/gomory.h"

void Gomory::update(double *rhs, int const *vBasis, int const *cBasis)
{
    // <= constraints - relax if the constraint is non-binding.
    for (size_t con = 0; con != d_problem.d_ss_leq; ++con)
        if (cBasis[con] == 0)
            rhs[con] = arma::datum::inf;

    // >= constraints
    for (size_t con = d_problem.d_ss_leq;
         con != d_problem.d_ss_leq + d_problem.d_ss_geq;
         ++con)
        if (cBasis[con] == 0)
            rhs[con] = -arma::datum::inf;

    auto const &Wmat = d_problem.Wmat();

    d_model.set(GRB_DoubleAttr_RHS, d_constrs, rhs, Wmat.n_cols);

    double lb[Wmat.n_rows];  // relax appropriate variable bounds.
    double ub[Wmat.n_rows];

    for (size_t var = 0; var != Wmat.n_rows; ++var)
    {
        // lb = -infinity if if var not at lower bound
        lb[var] = (vBasis[var] == -1) ? d_problem.d_l2[var] : -arma::datum::inf;

        // ub =  infinity if x[var] not at upper bound
        ub[var] = (vBasis[var] == -2) ? d_problem.d_u2[var] : arma::datum::inf;
    }

    // Sets the bounds to their appropriate values.
    d_model.set(GRB_DoubleAttr_LB, d_vars, lb, Wmat.n_rows);
    d_model.set(GRB_DoubleAttr_UB, d_vars, ub, Wmat.n_rows);
}
