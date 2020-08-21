#include "CMMonopoleSolver.h"
#include <iostream>
#include <cmath>

using namespace std;

int main(int argc, char const *argv[])
{
    CMMonopoleSolver solver;
    double g=solver.Getgweak();
    double rho0=solver.GetVEV();
    double mw=solver.GetMW();
    // cout<<"g: "<<solver.Getgweak()<<endl;
    double A0=mw/2.0;

    VD left = {0,1,0,-A0/rho0};
    VD right = {1,0,A0/rho0,0};

    solver.SetBoundary(left,right);
    solver.SetMeshPoints(200);
    solver.SetXRange(0.2,50);

    VD X;
    VVD Y;
    bool good;
    
    solver.SetMHL(10.0);
    good = solver.Solve(X,Y);
    solver.DumpSolution("CMMonopole_sol_ms10.dat");

    solver.SetXRange(0.2,50);
    solver.SetMHL(125.0);
    good = solver.Solve(X,Y);
    solver.DumpSolution("CMMonopole_sol_ms125.dat");
    solver.DumpSolution("CMMonopole_sol_ms125_min0x2.dat");
    double KA,KB,KPhi,VPhi;
    solver.GetEnergy(KA,KPhi,VPhi,KB);
    cout<<"The Energy is: "<<endl;
    cout<<"KA: "<<KA<<endl;
    cout<<"KB: "<<KB<<endl;
    cout<<"KPhi: "<<KPhi<<endl;
    cout<<"VPhi: "<<VPhi<<endl;
    cout<<"Energy total: "<<KA+KPhi+VPhi+KB<<endl;

    solver.SetXRange(0.3,50);
    solver.SetMHL(125.0);
    good = solver.Solve(X,Y);
    solver.DumpSolution("CMMonopole_sol_ms125_min0x3.dat");
    // double KA,KB,KPhi,VPhi;
    solver.GetEnergy(KA,KPhi,VPhi,KB);
    cout<<"The Energy is: "<<endl;
    cout<<"KA: "<<KA<<endl;
    cout<<"KB: "<<KB<<endl;
    cout<<"KPhi: "<<KPhi<<endl;
    cout<<"VPhi: "<<VPhi<<endl;
    cout<<"Energy total: "<<KA+KPhi+VPhi+KB<<endl;

    solver.SetMHL(1000.0);
    good = solver.Solve(X,Y);
    solver.DumpSolution("CMMonopole_sol_ms1000.dat");

    solver.SetMHL(10000.0);
    solver.SetMeshPoints(500);
    solver.SetXRange(0.5,20);
    good = solver.Solve(X,Y);
    solver.DumpSolution("CMMonopole_sol_ms10000.dat");

    VD left0 = {0,1,0,0};
    VD right0 = {1,0,A0/rho0,0};

    solver.SetBoundary(left0,right0);
    solver.SetMeshPoints(500);
    solver.SetXRange(0.3,50);

    solver.SetMHL(10.0);
    good = solver.Solve(X,Y);
    solver.DumpSolution("CMMonopole_sol_ms10_b00.dat");

    solver.SetMHL(125.0);
    good = solver.Solve(X,Y);
    solver.DumpSolution("CMMonopole_sol_ms125_b00.dat");
    // solver.GetEnergy(KA,KPhi,VPhi,KB);
    // cout<<"The Energy is: "<<endl;
    // cout<<"KA: "<<KA<<endl;
    // cout<<"KB: "<<KB<<endl;
    // cout<<"KPhi: "<<KPhi<<endl;
    // cout<<"VPhi: "<<VPhi<<endl;
    // cout<<"Energy total: "<<KA+KPhi+VPhi+KB<<endl;

    solver.SetMHL(1000.0);
    good = solver.Solve(X,Y);
    solver.DumpSolution("CMMonopole_sol_ms1000_b00.dat");

    solver.SetMHL(10000.0);
    solver.SetMeshPoints(500);
    solver.SetXRange(0.5,20);
    good = solver.Solve(X,Y);
    solver.DumpSolution("CMMonopole_sol_ms10000_b00.dat");
    return 0;
}
