#include "sub.h"

void Sub::update(arma::vec &rhs)
{
    d_model.set(GRB_DoubleAttr_RHS, d_constrs, rhs.memptr(), d_m2);
}
