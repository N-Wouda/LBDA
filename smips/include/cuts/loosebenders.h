#ifndef LOOSEBENDERS_H
#define LOOSEBENDERS_H

#include "cut.h"
#include "decompositions/gomory.h"


class LooseBenders : public Cut
{
    arma::vec const &d_alpha;

    Gomory d_gomory;

    // For each scenario, we store the basis matrices that we have visited
    // (encoded by vBasis, cBasis).
    std::vector<std::vector<std::vector<double>>> d_visited;

    // For each visited basis matrix, we store the corresponding Gomory
    // objective value.
    std::vector<std::vector<double>> d_objectives;

    double computeGomory(size_t s, int *vBasis, int *cBasis, arma::vec &rhs);

public:
    LooseBenders(GRBEnv &env,
                 Problem const &problem,
                 arma::vec const &alpha,
                 double timeLimit = 1e6);

    CutResult computeCut(arma::vec const &x) override;
};

#endif  // LOOSEBENDERS_H
