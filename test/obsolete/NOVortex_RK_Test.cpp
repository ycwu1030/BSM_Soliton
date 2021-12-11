#include "RungeKutta.h"
#include <iostream>
#include <cmath>
using namespace std;


const int n = 1;
const double e2 = 0.091725;
const double lam = 1;

// NO type Vortex:
VD ODE1(double x, VD y)
{
    VD res;
    res.push_back(-1/x*y[0]+n*n/x/x*y[1]*pow(y[3]-1,2)+lam/2*y[1]*(y[1]*y[1]-1));
    res.push_back(y[0]);
    res.push_back(1/x*y[2]+2*e2*y[1]*y[1]*(y[3]-1));
    res.push_back(y[2]);
    return res;
}

int main(int argc, char const *argv[])
{
    VD BOUND_end = {5e-7,1,1e-5,1};
    RungeKutta RK(4);
    RK.SetBound(10,0.001,BOUND_end);
    RK.SetODE(ODE1);
    RK.ODEINTEGRAL(0.001,1e-4);
    RK.DumpSolution("RK_Vortex_TEST.dat");
    return 0;
}
