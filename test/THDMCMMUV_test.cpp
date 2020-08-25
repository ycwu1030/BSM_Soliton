#include "THDM_CMMonopoleSolverUVReg.h"
#include <iostream>
#include <cmath>

using namespace std;

int main(int argc, char const *argv[])
{
    THDMCMMSolverUV solver;
    double sw=solver.GetSW();
    double g=solver.Getgweak();
    double rho0=solver.GetVEV();
    double mw=solver.GetMW();
    // cout<<"g: "<<solver.Getgweak()<<endl;
    double A0=mw/2.0;
    
    double beta = atan(1.0);//atan(2.0);
    double cba = 0.0;//0.1;
    double alpha = beta - acos(cba);
    double mhh = 400;
    double mha = 400;//450;
    double mpm = 400;//500;
    double m122 = 40000;//80000;

    solver.Set_Physical_Parameters(beta,alpha,mhh,mha,mpm,m122);

    VD left = {0,0,1.0/sw,0,-0.2*g};
    VD right = {cos(beta),sin(beta),0,A0/rho0,0};
    
    solver.SetBoundary(left,right);
    solver.SetUVRegular(0.0);
    solver.SetMeshPoints(200);
    solver.SetXRange(0.5,50);

    VD X;
    VVD Y;
    bool good;
    double E0,E1;

    solver.ExtendtoZero(false);
    good = solver.Solve(X,Y);
    solver.DumpSolution("THDMCMMUV_fa_sol_test.dat");

    solver.GetEnergy(E0,E1);
    cout<<"The Energy is: "<<endl;
    cout<<"E0: "<<E0<<endl;
    cout<<"E1: "<<E1<<endl;
    cout<<"Energy total: "<<E0+E1<<endl;

    solver.ExtendtoZero(true);
    good = solver.Solve(X,Y);
    solver.DumpSolution("THDMCMMUV_tr_sol_test.dat");

    solver.GetEnergy(E0,E1);
    cout<<"The Energy is: "<<endl;
    cout<<"E0: "<<E0<<endl;
    cout<<"E1: "<<E1<<endl;
    cout<<"Energy total: "<<E0+E1<<endl;


    return 0;
}
