#ifndef SMIPS_SOLVE_H
#define SMIPS_SOLVE_H

#include "decompositions/decomposition.h"
#include "masterproblem.h"

#include <armadillo>

// TODO STATISTICS - number of cuts, run-time?

/**
 * TODO this should be a method on MasterProblem (merge with
 *   MasterProblem::solve?)
 * Solves the master problem using the decomposition strategy. The master
 * problem is solved with an optimality gap smaller than tol.
 *
 * @return The (near) optimal first-stage decisions (x).
 */
std::unique_ptr<arma::vec> solve(MasterProblem &master,
                                 Decomposition &decomposition,
                                 double tol = 1e-4);

#endif  // SMIPS_SOLVE_H
