#include "benders.h"

Benders::Benders(Master &master) :
    d_nCuts(0),
    d_runTime(0),
    d_master(master)
{
}
