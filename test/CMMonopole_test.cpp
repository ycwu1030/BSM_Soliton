#include "CMMonopoleSolver.h"
#include <iostream>
#include <cmath>

using namespace std;

int main(int argc, char const *argv[])
{
    CMMonopoleSolver solver;
    cout<<"g: "<<solver.Getgweak()<<endl;
    VD left = {0,1,0,-0.21*solver.Getgweak()};
    VD right = {1,0,0.25*solver.Getgweak(),0};

    solver.SetBoundary(left,right);
    solver.SetMeshPoints(200);
    solver.SetXRange(0.5,50);

    VD X;
    VVD Y;
    bool good;
    
    solver.SetMHL(10.0);
    good = solver.Solve(X,Y);
    solver.DumpSolution("CMMonopole_sol_ms10.dat");

    solver.SetMHL(1000.0);
    good = solver.Solve(X,Y);
    solver.DumpSolution("CMMonopole_sol_ms1000.dat");

    solver.SetMHL(10000.0);
    solver.SetMeshPoints(500);
    solver.SetXRange(0.5,20);
    good = solver.Solve(X,Y);
    solver.DumpSolution("CMMonopole_sol_ms10000.dat");
    return 0;
}
