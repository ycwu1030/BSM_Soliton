#include "SM_cxSM.h"
#include "SolitonSolver.h"
#include "RungeKutta.h"
#include <iostream>
using namespace std;


CXSM model;
VD ODE1(double x, VD y)
{
    VD res;
    double phiH = y[0]*model.vev;
    double phiS = y[1]*model.vev;
    res.push_back(model.dVdH(phiH,phiS,0)/pow(model.vev,3));
    res.push_back(y[0]);
    res.push_back(model.dVdS(phiH,phiS,0)/pow(model.vev,3));
    res.push_back(y[2]);
    return res;
}
int main(int argc, char const *argv[])
{
    model.SetInput(157.633,42.647,0.0,0.335465);
    model.PrintPotentialParameter();
    model.PrintLocalMinima();
    cout<<"Stability: "<<model.CheckStability()<<endl;
    cout<<"Unitarity: "<<model.CheckUnitarity()<<endl;
    cout<<"Global Min: "<<model.CheckGlobalMinimum()<<endl;
    
    VD BOUND_begin = {0,model.vev/model.vev,0,-model.VS/model.vev};
    VD BOUND_end = {0,model.vev/model.vev,0,model.VS/model.vev};
    VD BOUND_guess = {1e-5,1e-5};
    VD delta_Bound = {1e-2,1e-2};
    VD eps = {1e-4,1e-4};
    vector<bool> bound_beg_Q = {false, true, false, true};
    vector<bool> bound_end_Q = {false, true, false, true};
    RungeKutta RK(4);
    RK.SetBound(-25,25,BOUND_begin,BOUND_end,bound_beg_Q,bound_end_Q);
    RK.SetODE(ODE1);
    RK.SHOOTING(BOUND_guess, delta_Bound, eps);
    RK.DumpSolution("SHOOTING_Soliton_TEST.dat");
    return 0;
}
