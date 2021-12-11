#include "RungeKutta.h"
#include <iostream>
using namespace std;


// Test Case 1:
VD ODE1(double x, VD y)
{
    // y[0], y[1];
    // dy[0]/dx = y[0](2-y[1])
    // dy[1]/dx = y[1](y[0]-1)
    // y[0](0) = 1;
    // y[1](0) = 2.7;
    // From x=0 to x=10
    VD res;
    res.push_back(y[0]*(2-y[1]));
    res.push_back(y[1]*(y[0]-1));
    return res;
}

int main(int argc, char const *argv[])
{
    VD BOUND = {1,2.7};
    RungeKutta RK(2);
    RK.SetBound(0,10,BOUND);
    RK.SetODE(ODE1);
    RK.ODEINTEGRAL(0.01);
    // RK.PrintSolution();
    RK.DumpSolution("RK_TEST.dat");
    return 0;
}
