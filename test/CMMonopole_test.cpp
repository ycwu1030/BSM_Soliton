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
    bool good = solver.Solve(X,Y);
    solver.DumpSolution("CMMonopole_sol_paper.dat");
    return 0;
}
