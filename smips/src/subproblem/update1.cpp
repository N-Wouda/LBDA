#include "subproblem.h"

void SubProblem::update(arma::vec &rhs)
{
    d_model.set(GRB_DoubleAttr_RHS, d_constrs, rhs.memptr(), d_problem.d_Wmat.n_cols);
}
