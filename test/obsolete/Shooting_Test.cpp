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
    VD BOUND_begin = {1,0};
    VD BOUND_end = {0,2.11533};
    VD BOUND_guess = {2.1};
    VD delta_Bound = {0.05};
    VD eps = {1e-5};
    vector<bool> bound_beg_Q = {true, false};
    vector<bool> bound_end_Q = {false, true};
    RungeKutta RK(2);
    RK.SetBound(0,10,BOUND_begin,BOUND_end,bound_beg_Q,bound_end_Q);
    RK.SetODE(ODE1);
    RK.SHOOTING(BOUND_guess, delta_Bound, eps);
    RK.DumpSolution("SHOOTING_TEST.dat");
    return 0;
}
