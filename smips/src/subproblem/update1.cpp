#include "subproblem.h"

void SubProblem::update(arma::vec &rhs)
{
    d_model.set(GRB_DoubleAttr_RHS, d_constrs, rhs.memptr(), rhs.n_elem);
}