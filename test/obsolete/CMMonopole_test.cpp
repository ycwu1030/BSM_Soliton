#include "CMMonopoleSolver.h"
#include <iostream>
#include <cmath>

using namespace std;
void RunSolver(CMMonopoleSolver &sol, double mh, double xmin, double xmax, bool print_energy=false, bool b00=false)
{
    sol.SetXRange(xmin,xmax);
    sol.SetMHL(mh);
    VD X; VVD Y;
    bool good = sol.Solve(X,Y);
    char tmp[200];
    if (b00)
    {
        sprintf(tmp,"CMMonopole_sol_mh%.1f_xmin%.1f_b00.dat",mh,xmin);
    }
    else
    {
        sprintf(tmp,"CMMonopole_sol_mh%.1f_xmin%.1f.dat",mh,xmin);
    }
    sol.DumpSolution(tmp);
    if (print_energy)
    {
        double KA,KB,KPhi,VPhi;
        sol.GetEnergy(KA,KPhi,VPhi,KB);
        cout<<"The Energy for mh = "<<mh<<" x_min = "<<xmin<<" is: "<<endl;
        cout<<"KA: "<<KA<<endl;
        cout<<"KB: "<<KB<<endl;
        cout<<"KPhi: "<<KPhi<<endl;
        cout<<"VPhi: "<<VPhi<<endl;
        cout<<"Energy total: "<<KA+KPhi+VPhi+KB<<endl;
    }
}

int main(int argc, char const *argv[])
{
    CMMonopoleSolver solver;
    double g=solver.Getgweak();
    double rho0=solver.GetVEV();
    double mw=solver.GetMW();
    // cout<<"g: "<<solver.Getgweak()<<endl;
    double A0=mw/2.0;

    VD left = {0,1,0,-0.2*g};
    VD right = {1,0,A0/rho0,0};

    solver.SetBoundary(left,right);
    solver.SetMeshPoints(200);

    RunSolver(solver,10,0.2,50);

    RunSolver(solver,125,0.2,50,true);
    RunSolver(solver,125,0.3,50,true);
    RunSolver(solver,125,0.5,50,true);
    RunSolver(solver,125,0.8,50,true);

    RunSolver(solver,1000,0.5,50);

    solver.SetMeshPoints(500);
    RunSolver(solver,10000,0.5,20);

    
    VD left0 = {0,1,0,0};
    VD right0 = {1,0,A0/rho0,0};
    solver.SetBoundary(left0,right0);
    solver.SetMeshPoints(400);


    RunSolver(solver,10,0.2,50,false,true);

    RunSolver(solver,125,0.2,50,true,true);
    RunSolver(solver,125,0.3,50,true,true);
    RunSolver(solver,125,0.5,50,true,true);
    RunSolver(solver,125,0.8,50,true,true);

    RunSolver(solver,1000,0.5,50,false,true);

    solver.SetMeshPoints(500);
    RunSolver(solver,10000,0.5,20,false,true);


    return 0;
}
