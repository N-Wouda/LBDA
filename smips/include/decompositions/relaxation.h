#ifndef RELAXATION_H
#define RELAXATION_H

#include "problem.h"

#include <gurobi_c++.h>

class Relaxation  // TODO think about name
{
protected:
    GRBModel d_model;
    Problem const &d_problem;

public:
    explicit Relaxation(GRBEnv &env, Problem const &problem);

    Relaxation(Relaxation const &other) = default;

    virtual ~Relaxation() = default;

    virtual double solve();
};

#endif  // RELAXATION_H
