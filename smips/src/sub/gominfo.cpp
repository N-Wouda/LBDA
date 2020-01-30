#include "sub.h"


Sub::GomInfo const Sub::gomInfo()
{
    // clang-format off
    return GomInfo
    {
        arma::vec(d_model.get(GRB_DoubleAttr_Pi, d_constrs, d_m2), d_m2),
        arma::Col<int>(d_model.get(GRB_IntAttr_VBasis, d_vars, d_n2), d_n2),
        arma::Col<int>(d_model.get(GRB_IntAttr_CBasis, d_constrs, d_m2), d_m2)
    };
    // clang-format on
}
