#include "cxSM_CP.h"

cxSM_CP::cxSM_CP():Basic_Model(3)
{

}

void cxSM_CP::Set_Potential_Parameters(double mu2, double lam, double del2, double b2, double d2, double del3, double b1, double d1, double d3)
{
    _mu2 = mu2;
    _lam = lam;
    _del2 = del2;
    _b2 = b2;
    _d2 = d2;
    _del3 = del3;
    _b1 = b1;
    _d1 = d1;
    _d3 = d3;
}