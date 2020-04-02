#include "problem.h"

void Problem::initSub()
{
    char vTypes[d_n2];
    std::fill(vTypes, vTypes + nSecondStageIntVars(), GRB_INTEGER);
    std::fill(vTypes + nSecondStageIntVars(), vTypes + d_n2, GRB_CONTINUOUS);

    GRBVar *vars = d_sub.addVars(d_l2.memptr(),
                                 d_u2.memptr(),
                                 d_q.memptr(),
                                 vTypes,
                                 nullptr,
                                 d_n2);

    // constraint senses
    char senses[d_Wmat.n_cols];
    std::fill(senses, senses + d_ss_leq, GRB_LESS_EQUAL);
    std::fill(senses + d_ss_leq,
              senses + d_ss_leq + d_ss_geq,
              GRB_GREATER_EQUAL);
    std::fill(senses + d_ss_leq + d_ss_geq, senses + d_Wmat.n_cols, GRB_EQUAL);

    // constraint rhs
    double rhs[d_Wmat.n_cols];
    std::fill(rhs, rhs + d_Wmat.n_cols, 0.0);

    GRBLinExpr Wy[d_Wmat.n_cols];
    for (size_t conIdx = 0; conIdx != d_Wmat.n_cols; ++conIdx)
        Wy[conIdx].addTerms(d_Wmat.colptr(conIdx), vars, d_Wmat.n_rows);

    // add constraints
    d_constrs = d_sub.addConstrs(Wy, senses, rhs, nullptr, d_Wmat.n_cols);
    d_sub_initialized = true;

    delete[] vars;
}
