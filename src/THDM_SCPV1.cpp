#include "THDM_SCPV1.h"


using namespace std;

THDM_SCPV1::THDM_SCPV1():Basic_Model(3)
{
    // The basic DOF is 3
}

void THDM_SCPV1::FindLocalMinima()
{
    if (_Solved)
    {
        return;
    }

    Clear_Local_Cache();

    // * 1. phi_I = 0;
    
}