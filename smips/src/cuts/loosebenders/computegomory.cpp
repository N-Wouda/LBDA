#include "cuts/loosebenders.h"

#include <algorithm>

double LooseBenders::computeGomory(size_t s,
                                   arma::vec &rhs,
                                   arma::Col<int> const &vBasis,
                                   arma::Col<int> const &cBasis)
{
    // TODO overhaul this method
    std::vector<int> basis(d_problem.d_n2 + d_problem.d_Wmat.n_cols);
    std::copy(vBasis.memptr(), vBasis.memptr() + d_problem.d_n2, basis.begin());
    std::copy(cBasis.memptr(),
              cBasis.memptr() + d_problem.d_Wmat.n_cols,
              basis.begin() + d_problem.d_n2);

    auto &visited = d_visited[s];
    auto it = std::find(visited.begin(), visited.end(), basis);

    if (it != visited.end())
    {
        // find index and retrieve corresponding objective value
        size_t idx = std::distance(visited.begin(), it);
        return d_objectives[s][idx];
    }

    d_gomory.update(rhs.memptr(), vBasis.memptr(), cBasis.memptr());
    double gom_obj = d_gomory.solve();

    visited.emplace_back(basis);
    d_objectives[s].emplace_back(gom_obj);

    return gom_obj;
}
