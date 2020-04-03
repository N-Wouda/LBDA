#include "masterproblem.h"

MasterProblem::MasterProblem(MasterProblem const &other) :
    d_problem(other.d_problem),
    d_xCoeffs(other.d_xCoeffs),
    d_cuts(other.d_cuts),
    d_nSlacks(other.d_nSlacks)
{
    GRBupdatemodel(other.d_cmodel);
    d_cmodel = GRBcopymodel(other.d_cmodel);
}
