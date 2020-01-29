#ifndef LPDUAL_H
#define LPDUAL_H

#include "cut.h"


class LpDual : public Cut
{
public:
    LpDual(GRBEnv &env, Problem const &problem);

    CutResult computeCut(arma::vec const &x) override;
};

#endif  // LPDUAL_H
