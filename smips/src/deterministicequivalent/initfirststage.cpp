#include "deterministicequivalent.h"

void DeterministicEquivalent::initFirstStage()
{
    auto const &Amat = d_problem.Amat();

    // variables
    char vTypes[Amat.n_rows];
    std::fill_n(vTypes, d_problem.nFirstStageIntVars(), GRB_INTEGER);
    std::fill_n(vTypes + d_problem.nFirstStageIntVars(),
                Amat.n_rows - d_problem.nFirstStageIntVars(),
                GRB_CONTINUOUS);
    d_xVars = d_model.addVars(d_problem.d_l1.memptr(),
                              d_problem.d_u1.memptr(),
                              d_problem.d_c.memptr(),
                              vTypes,
                              nullptr,
                              Amat.n_rows);

    // constraints
    GRBLinExpr lhsExprs[Amat.n_cols];
    for (size_t conIdx = 0; conIdx != Amat.n_cols; ++conIdx)
        lhsExprs[conIdx].addTerms(Amat.colptr(conIdx), d_xVars, Amat.n_rows);

    char senses[Amat.n_cols];
    std::fill(senses, senses + d_problem.d_fs_leq, GRB_LESS_EQUAL);
    std::fill(senses + d_problem.d_fs_leq,
              senses + d_problem.d_fs_leq + d_problem.d_fs_geq,
              GRB_GREATER_EQUAL);
    std::fill(senses + d_problem.d_fs_leq + d_problem.d_fs_geq,
              senses + Amat.n_cols,
              GRB_EQUAL);

    GRBConstr *constrs = d_model.addConstrs(lhsExprs,
                                            senses,
                                            d_problem.d_b.memptr(),
                                            nullptr,
                                            Amat.n_cols);

    delete[] constrs;
}
