#include "decompositions/loosebenders.h"

void LooseBenders::update(arma::vec &rhs,
                          arma::Col<int> const &vBasis,
                          arma::Col<int> const &cBasis)
{
    // TODO clean this up
    // <= constraints - relax if the constraint is non-binding.
    for (size_t con = 0; con != d_problem.d_nSecondStageLeqConstraints; ++con)
        if (cBasis[con] == GRB_BASIC)
            rhs[con] = arma::datum::inf;

    // >= constraints
    for (size_t con = d_problem.d_nSecondStageLeqConstraints;
         con
         != d_problem.d_nSecondStageLeqConstraints
                + d_problem.d_nSecondStageGeqConstraints;
         ++con)
        if (cBasis[con] == GRB_BASIC)
            rhs[con] = -arma::datum::inf;

    auto const &Wmat = d_problem.Wmat();

    d_model.set(GRB_DoubleAttr_RHS, d_constrs, rhs.memptr(), Wmat.n_cols);

    arma::vec lb = d_problem.d_secondStageLowerBound;
    arma::vec ub = d_problem.d_secondStageUpperBound;

    // Relax appropriate variable bounds if the bound is not tight.
    lb.elem(arma::find(vBasis != GRB_NONBASIC_LOWER)).fill(-arma::datum::inf);
    ub.elem(arma::find(vBasis != GRB_NONBASIC_UPPER)).fill(arma::datum::inf);

    d_model.set(GRB_DoubleAttr_LB, d_vars, lb.memptr(), Wmat.n_rows);
    d_model.set(GRB_DoubleAttr_UB, d_vars, ub.memptr(), Wmat.n_rows);
}
