#ifndef GOMORY_H
#define GOMORY_H

#include "problem.h"
#include "relaxation.h"

#include <gurobi_c++.h>


class Gomory : public Relaxation
{
    GRBConstr *d_constrs;
    GRBVar *d_vars;

public:
    Gomory(GRBEnv &env, Problem const &problem);

    Gomory(Gomory const &other);

    ~Gomory() override;

    void update(double *rhs, int const *vBasis, int const *cBasis);

    void setTimeLimit(double timeLimit);
};

#endif
