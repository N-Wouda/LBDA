#include "data.h"

std::vector<std::vector<double>> Data::rand_gaus_mat(size_t nRows,
                                                     size_t nCols,
                                                     double mean,
                                                     double sd)
{
    std::normal_distribution<double> gaus(mean, sd);
    std::vector<std::vector<double>> mat(nRows);  // initialize
                                                  // and fill element-wise
    for (size_t row = 0; row != nRows; ++row)
        for (size_t col = 0; col != nCols; ++col)
            mat[row].push_back(gaus(d_engine));

    return mat;
}
