#include "cxSM_CP_reduced.h"
#include "PathDeformation.h"
// #include "Relaxation.h"
#include <iostream>
#include <iomanip>
#include <fstream>
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
    model.GetVS(vsr,vsi);
    VD left = {model.GetVEV(),vsr,-vsi};
    VD right = {model.GetVEV(),vsr,vsi};
    VVD pts_init;
    pts_init.push_back(left);
    pts_init.push_back(right);
    VnD vtol = [&](VD field, double scale){
        return model.Vtotal(field,scale);
    };
    dVnD dvtol = [&](VD field, double scale){
        return model.dVtotal(field,scale);
    };
    KinknD sol = fullKink(pts_init,vtol,dvtol);
    // ofstream output("cxSMCP_reduced_PDDW_test.dat");
    // output<<"r\tphi0\tphi1\tphi2\tphi1D\tdphi1D"<<endl;
    // for (int i = 0; i < sol.R.size(); i++)
    // {
    //     output<<scientific<<setprecision(10)<<sol.R[i]<<"\t"<<sol.Phi[i]<<"\t"<<sol.Phi_1D[i]<<"\t"<<sol.dPhi_1D[i]<<endl;
    // }
    // model.DumpEnergyDensity(sol.R,sol.Phi,"cxSMCP_reduced_PDDW_density.dat");
    model.DumpFullSolution(sol.R,sol.Phi,"cxSMCP_reduced_PDDW_full.dat");
    cout<<"Energy: "<<model.GetTotalEnergy(sol.R,sol.Phi)<<endl;
    cout<<"Tension: "<<model.GetTension(sol.R,sol.Phi)<<endl;
    // RX.DumpSolution("Relaxation_DW_test.dat");
    return 0;
}

