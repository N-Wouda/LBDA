#include "subproblem.h"

SubProblem::SubProblem(GRBEnv &env, Problem const &problem) :
    d_model(env),
    d_problem(problem)
{
    auto const &Wmat = d_problem.Wmat();

    // variable types
    char vTypes[Wmat.n_rows];
    std::fill(vTypes, vTypes + Wmat.n_rows, GRB_CONTINUOUS);

    d_vars = d_model.addVars(problem.d_secondStageLowerBound.memptr(),
                             problem.d_secondStageUpperBound.memptr(),
                             problem.d_secondStageCoeffs.memptr(),
                             vTypes,
                             nullptr,
                             Wmat.n_rows);

    size_t const ss_leq = problem.d_nSecondStageLeqConstraints;
    size_t const ss_geq = problem.d_nSecondStageGeqConstraints;

    // constraint senses
    char senses[Wmat.n_cols];
    std::fill(senses, senses + ss_leq, GRB_LESS_EQUAL);
    std::fill(senses + ss_leq, senses + ss_leq + ss_geq, GRB_GREATER_EQUAL);
    std::fill(senses + ss_leq + ss_geq, senses + Wmat.n_cols, GRB_EQUAL);

    // constraint rhs
    double rhs[Wmat.n_cols];
    std::fill(rhs, rhs + Wmat.n_cols, 0.0);

    // constraint lhs
    GRBLinExpr Wy[Wmat.n_cols];
    for (size_t conIdx = 0; conIdx != Wmat.n_cols; ++conIdx)
        Wy[conIdx].addTerms(Wmat.colptr(conIdx), d_vars, Wmat.n_rows);

    d_constrs = d_model.addConstrs(Wy, senses, rhs, nullptr, Wmat.n_cols);
    d_model.update();
}
