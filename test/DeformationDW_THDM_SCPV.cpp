#include "THDM_SCPV1.h"
#include "PathDeformation.h"
// #include "Relaxation.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

using namespace std;


int main(int argc, char const *argv[])
{
    THDM_SCPV1 model;
    // model.SetInput(157.633,42.647,0.335465);
    // model.Set_Physical_Parameters(157.633,0.335465,42.647);
    // model.Set_Physical_Parameters(15.0,0.1,66.0);
    double beta = atan(2);
    double MHH = 400;
    double MHA = 410;
    double MHpm = 410;
    double alpha = beta-M_PI_2;
    double alphac = 1e-5;
    double theta;
    bool good = model.Set_Physical_Parameters(beta,MHH,MHA,MHpm,alpha,alphac);
    model.PrintParameters();
    cout<<"Unitarity: "<<model.CheckUnitarity()<<endl;
    cout<<"Stability: "<<model.CheckStability()<<endl;
    cout<<"Global Min: "<<model.CheckGlobalMinimum()<<endl;
    // model.GetVS(vsr,vsi);
    model.GetTheta(theta);
    double vev = model.GetVEV();
    double v1 = vev*cos(beta);
    double v2r = vev*sin(beta)*cos(theta);
    double v2i = vev*sin(beta)*sin(theta);
    VD left = {v1,v2r,-v2i};
    VD right = {v1,v2r,v2i};
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
    model.DumpFullSolution(sol.R,sol.Phi,"THDM_SCPV_PDDW_full.dat");
    cout<<"Energy: "<<model.GetTotalEnergy(sol.R,sol.Phi)<<endl;
    cout<<"Tension: "<<model.GetTension(sol.R,sol.Phi)<<endl;
    // RX.DumpSolution("Relaxation_DW_test.dat");
    return 0;
}

