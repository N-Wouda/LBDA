#include "deterministicequivalent.h"


std::unique_ptr<arma::vec> DeterministicEquivalent::solve(double time_limit)
{
    d_model.set(GRB_DoubleParam_TimeLimit, time_limit);
    d_model.optimize();

    int status = d_model.get(GRB_IntAttr_Status);

    if (status == GRB_INFEASIBLE)
    {
        d_status = status::INFEASIBLE;
        // TODO check this.
        return std::make_unique<arma::vec>(nullptr);
    }

    d_status = status::SOLVED;

    if (d_problem.isMixedIntegerProblem())
        d_MIPGap = d_model.get(GRB_DoubleAttr_MIPGap);

    d_objVal = d_model.get(GRB_DoubleAttr_ObjVal);
    d_objBound = d_model.get(GRB_DoubleAttr_ObjBound);
    d_runTime = d_model.get(GRB_DoubleAttr_Runtime);

    auto const &Amat = d_problem.Amat();

    arma::vec xVals(d_model.get(GRB_DoubleAttr_X, d_xVars, Amat.n_rows),
                    Amat.n_rows);
    return std::make_unique<arma::vec>(xVals);
}
