#include "subproblem.h"


SubProblem::GomInfo const SubProblem::gomInfo()
{
    // clang-format off
    return GomInfo
    {
        arma::vec(d_model.get(GRB_DoubleAttr_Pi, d_constrs, d_problem.d_Wmat.n_cols),
                  d_problem.d_Wmat.n_cols),
        arma::Col<int>(d_model.get(GRB_IntAttr_VBasis, d_vars, d_problem.d_n2),
                       d_problem.d_n2),
        arma::Col<int>(d_model.get(GRB_IntAttr_CBasis, d_constrs, d_problem.d_Wmat.n_cols),
                       d_problem.d_Wmat.n_cols)
    };
    // clang-format on
}
