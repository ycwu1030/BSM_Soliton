/*
 * @Description  : 
 * @Author       : Yongcheng Wu
 * @Date         : 2020-01-27 17:56:36
 * @LastEditors  : Yongcheng Wu
 * @LastEditTime : 2020-02-06 00:27:17
 */

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
    VD BOUND_begin = {0,0,0,0};
    VD BOUND_end = {0,1,0,1};
    VD BOUND_guess = {1e-5,1e-5};
    VD delta_Bound = {1e-3,1e-3};
    VD eps = {1e-4,1e-4};
    vector<bool> bound_beg_Q = {false, true, false, true};
    vector<bool> bound_end_Q = {false, true, false, true};
    RungeKutta RK(4);
    RK.SetBound(10,0.001,BOUND_end,BOUND_begin,bound_end_Q,bound_beg_Q);
    RK.SetODE(ODE1);
    RK.SHOOTING(BOUND_guess, delta_Bound, eps);
    RK.DumpSolution("SHOOTING_Vortex_TEST.dat");
    return 0;
}
