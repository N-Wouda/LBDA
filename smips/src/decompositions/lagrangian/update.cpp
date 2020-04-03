#include "decompositions/lagrangian.h"


void Lagrangian::update(arma::vec &rhs, arma::vec &pi)
{
    d_model.set(GRB_DoubleAttr_RHS, d_constrs, rhs.memptr(), d_problem.d_Wmat.n_cols);
    d_model.set(GRB_DoubleAttr_Obj, d_z_vars, pi.memptr(), d_problem.d_Amat.n_rows);
}
