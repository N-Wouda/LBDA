#ifndef SMIPS_SOLVE_H
#define SMIPS_SOLVE_H

#include "cuts/cut.h"
#include "masterproblem.h"

#include <armadillo>

// TODO STATISTICS - number of cuts, run-time?

std::unique_ptr<arma::vec> solve(MasterProblem &master,
                                 Cut &cut,
                                 double tol = 1e-4);

#endif  // SMIPS_SOLVE_H
