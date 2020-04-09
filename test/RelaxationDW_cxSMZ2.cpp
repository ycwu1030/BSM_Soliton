#include "cxSM_Z2.h"
#include "DWSolver.h"
// #include "Relaxation.h"
#include <iostream>
#include <cmath>

using namespace std;


int main(int argc, char const *argv[])
{
    cxSM_Z2 model;
    // model.SetInput(157.633,42.647,0.335465);
    // model.Set_Physical_Parameters(157.633,0.335465,42.647);
    model.Set_Physical_Parameters(15.0,0.1,66.0);
    model.PrintParameters();
    VD left = {model.GetVEV(),-15.};
    VD right = {model.GetVEV(),15.};
    DWSolver solver(&model,left,right);
    VD X;
    VVD Y;
    bool good = solver.Solve(X,Y);
    solver.DumpSolution("cxSMZ2_DW_sol.dat");
    cout<<"Energy: "<<model.GetTotalEnergy(X,Y)<<endl;
    cout<<"Tension: "<<model.GetTension(X,Y)<<endl;
    return 0;
}

