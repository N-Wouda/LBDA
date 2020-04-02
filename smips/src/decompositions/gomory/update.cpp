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

    d_model.set(GRB_DoubleAttr_RHS, d_constrs, rhs, d_problem.d_m2);

    double lb[d_problem.d_n2];  // relax appropriate variable bounds.
    double ub[d_problem.d_n2];

    for (size_t var = 0; var != d_problem.d_n2; ++var)
    {
        // lb = -infinity if if var not at lower bound
        lb[var] = (vBasis[var] == -1) ? d_problem.d_l2[var] : -arma::datum::inf;

        // ub =  infinity if x[var] not at upper bound
        ub[var] = (vBasis[var] == -2) ? d_problem.d_u2[var] : arma::datum::inf;
    }

    // Sets the bounds to their appropriate values.
    d_model.set(GRB_DoubleAttr_LB, d_vars, lb, d_problem.d_n2);
    d_model.set(GRB_DoubleAttr_UB, d_vars, ub, d_problem.d_n2);
}
