#include "THDM_CMMonopoleSolver.h"
#include <iostream>
#include <cmath>

using namespace std;

int main(int argc, char const *argv[])
{
    THDMCMMSolver solver;
    double g=solver.Getgweak();
    double rho0=solver.GetVEV();
    double mw=solver.GetMW();
    // cout<<"g: "<<solver.Getgweak()<<endl;
    double A0=mw/2.0;
    
    double beta = atan(2.0);
    double cba = 0.1;
    double alpha = beta - acos(cba);
    double mhh = 400;
    double mha = 450;
    double mpm = 500;
    double m122 = 80000;

    solver.Set_Physical_Parameters(beta,alpha,mhh,mha,mpm,m122);

    VD left = {0,0,1,0,-A0/rho0};
    VD right = {cos(beta),sin(beta),0,A0/rho0,0};

    solver.SetBoundary(left,right);
    solver.SetMeshPoints(200);
    solver.SetXRange(0.5,50);

    VD X;
    VVD Y;
    bool good;

    good = solver.Solve(X,Y);
    solver.DumpSolution("THDMCMM_sol_test.dat");



    return 0;
}
