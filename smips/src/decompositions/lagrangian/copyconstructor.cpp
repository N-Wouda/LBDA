#include "decompositions/lagrangian.h"

#include <algorithm>

Lagrangian::Lagrangian(Lagrangian const &other) :
    Relaxation(other),
    d_constrs(d_model.getConstrs())
{
    GRBVar *vars = d_model.getVars();

    auto const &Amat = d_problem.Amat();

    // We only need the first n1, those associated with the first-stage problem.
    d_z_vars = new GRBVar[Amat.n_rows];
    std::copy_n(vars, Amat.n_rows, d_z_vars);

    delete[] vars;
}
