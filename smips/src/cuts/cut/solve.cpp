#include "cuts/cut.h"


double Cut::solve()
{
    d_model.optimize();
    return d_model.get(GRB_DoubleAttr_ObjVal);
}
