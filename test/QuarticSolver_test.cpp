#include "QuarticSolver.h"
#include <iostream>
using namespace std;

int main(int argc, char const *argv[])
{
    QuarticSolver solver;
    solver.Solve();
    for (int i = 0; i < 4; i++)
    {
        cout<<solver.SOLUTIONS[i]<<endl;
    }
    vector<double> realsol = solver.GetRealSolution();
    if (realsol.size()>0)
    {
        for (int i = 0; i < realsol.size(); i++)
        {
            cout<<realsol[i]<<endl;
        }
    }
    else
    {
        cout<<"NO REAL SOLUTION"<<endl;
    }
    
    
    solver.Solve(1,3,4,2,-9);
    for (int i = 0; i < 4; i++)
    {
        cout<<solver.SOLUTIONS[i]<<endl;
    }
    realsol = solver.GetRealSolution();
    if (realsol.size()>0)
    {
        for (int i = 0; i < realsol.size(); i++)
        {
            cout<<realsol[i]<<endl;
        }
    }
    else
    {
        cout<<"NO REAL SOLUTION"<<endl;
    }

    solver.Solve(5,3,4,2,-9);
    for (int i = 0; i < 4; i++)
    {
        cout<<solver.SOLUTIONS[i]<<endl;
    }
    realsol = solver.GetRealSolution();
    if (realsol.size()>0)
    {
        for (int i = 0; i < realsol.size(); i++)
        {
            cout<<realsol[i]<<endl;
        }
    }
    else
    {
        cout<<"NO REAL SOLUTION"<<endl;
    }
    return 0;
}
