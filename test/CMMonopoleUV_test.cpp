#include "CMMonopoleSolverUVReg.h"
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;
double RunSolver(CMMonopoleSolverUV &sol, double mh, double xmin, double xmax, bool DumpSol=true, bool print_energy=false, bool b00=false)
{
    sol.SetXRange(xmin,xmax);
    sol.SetMHL(mh);
    VD X; VVD Y;
    bool good = sol.Solve(X,Y);
    char tmp[200];
    if (b00)
    {
        sprintf(tmp,"CMMonopole_UV_sol_mh%.1f_xmin%.1f_b00.dat",mh,xmin);
    }
    else
    {
        sprintf(tmp,"CMMonopole_UV_sol_mh%.1f_xmin%.1f.dat",mh,xmin);
    }
    if (DumpSol)
    {
        sol.DumpSolution(tmp);
    }
    double E0,E1;
    sol.GetEnergy(E0,E1);
    if (print_energy)
    {
        cout<<"The Energy for mh = "<<mh<<" x_min = "<<xmin<<" is: "<<endl;
        cout<<"KA: "<<E0<<endl;
        cout<<"KB: "<<E1<<endl;
        cout<<"Energy total: "<<E0+E1<<endl;
    }
    return E0+E1;
}

int main(int argc, char const *argv[])
{
    CMMonopoleSolverUV solver;
    double sw=solver.GetSW();
    double g=solver.Getgweak();
    double rho0=solver.GetVEV();
    double mw=solver.GetMW();
    cout<<"f0: "<<1.0/sw<<endl;
    double A0=mw/2.0;

    VD left = {0,1.0/sw,0,-0.2*g};
    VD right = {1,0,A0/rho0,0};
    double Eall;
    solver.SetBoundary(left,right);
    solver.SetUVRegular(0.0);
    solver.SetMeshPoints(200);

    Eall=RunSolver(solver,125,0.2,50,true,true);
    Eall=RunSolver(solver,125,0.3,50,true,true);
    Eall=RunSolver(solver,125,0.5,50,true,true);
    Eall=RunSolver(solver,125,0.8,50,true,true);

    double f0;
    double Mass;
    ofstream outmass("CMMUV_Mass_f0.dat");
    outmass<<"f0\tM"<<endl;
    for (int i = 0; i < 50; i++)
    {
        f0 = (0.5+i*(3.0-0.5)/(49.0));
        left[1] = f0;
        solver.SetBoundary(left,right);
        solver.SetUVRegular(0.0);
        Eall=RunSolver(solver,125,0.5,50,false);
        Mass = (Eall);
        outmass<<f0<<"\t"<<Mass<<endl;
    }
    outmass.close();

    return 0;
}
