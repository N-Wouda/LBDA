#include "decompositions/lagrangian.h"

#include <algorithm>

Lagrangian::Lagrangian(Lagrangian const &other) :
    Relaxation(other),
    d_constrs(d_model.getConstrs())
{
    GRBVar *vars = d_model.getVars();

    d_z_vars = new GRBVar[d_problem.d_n1];  // deallocated in destructor

    // We only need the first d_n1
    std::copy_n(vars, d_problem.d_n1, d_z_vars);

    delete[] vars;
}
