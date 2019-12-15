#include "ald.h"

Ald::~Ald()
{
    GRBfreemodel(d_model);
}
