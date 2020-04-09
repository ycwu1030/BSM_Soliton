#include "cxSM_CP_reduced.h"
#include "DWSolver.h"
// #include "Relaxation.h"
#include <iostream>
#include <cmath>

using namespace std;


int main(int argc, char const *argv[])
{
    cxSM_CP_reduced model;
    // model.SetInput(157.633,42.647,0.335465);
    // model.Set_Physical_Parameters(157.633,0.335465,42.647);
    // model.Set_Physical_Parameters(15.0,0.1,66.0);
    double vsr,vsi;
    double vs = 100000;
    double MHH = 10000;
    double MHA = 10100;
    double theta1 = 0.001;
    double theta3 = 0.0001;
    bool good = model.Set_Physical_Parameters_vs_theta(vs,MHH,MHA,theta1,theta3);
    model.PrintParameters();
    cout<<"Unitarity: "<<model.CheckUnitarity()<<endl;
    cout<<"Stability: "<<model.CheckStability()<<endl;
    cout<<"Global Min: "<<model.CheckGlobalMinimum()<<endl;
    DWSolver solver(&model);
    model.GetVS(vsr,vsi);
    VD left = {model.GetVEV(),vsr,-vsi};
    VD right = {model.GetVEV(),vsr,vsi};
    solver.SetBoundary(left,right);
    solver.SetOverallScale({model.GetVEV(),1,vsi});
    VD X;
    VVD Y;
    good = solver.Solve(X,Y);
    // solver.PrintSolution();
    solver.DumpSolution("cxSMCP_reduced_DW_sol.dat");
    model.DumpEnergyDensity(X,Y,"cxSMCP_reduced_Density.dat");
    cout<<"Energy: "<<model.GetTotalEnergy(X,Y)<<endl;
    cout<<"Tension: "<<model.GetTension(X,Y)<<endl;
    // RX.DumpSolution("Relaxation_DW_test.dat");
    return 0;
}

