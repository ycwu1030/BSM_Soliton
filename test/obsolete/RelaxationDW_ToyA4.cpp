#include "ToyA4.h"
#include "DWSolver.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

using namespace std;

int main(int argc, char const *argv[])
{
    ToyA4 model;
    model.Set_Potential_Parameters(-1.0,3.0/2.0,4.0);
    model.PrintParameters();
    model.PrintLocalMinima();

    double v1 = model.GetV1();
    double v2 = model.GetV2();

    VD left = {v1,0.0};
    VD right = {-v1,0.0};
    DWSolver sol(&model,left,right);
    sol.SetZRange(20);
    VD X;
    VVD Y;
    bool good = sol.Solve(X,Y);
    model.DumpFullSolution(X,Y,"ToyA4_2.dat");
    cout<<"Energy: "<<model.GetTotalEnergy(X,Y)<<endl;
    cout<<"Tension: "<<model.GetTension(X,Y)<<endl;

    VD left1 = {v1,0.1};
    VD right1 = {-v1,0.1};
    DWSolver sol1(&model,left1,right1);
    sol1.SetZRange(20);
    VD X1;
    VVD Y1;
    good = sol1.Solve(X1,Y1);
    model.DumpFullSolution(X1,Y1,"ToyA4_3.dat");
    cout<<"Energy: "<<model.GetTotalEnergy(X1,Y1)<<endl;
    cout<<"Tension: "<<model.GetTension(X1,Y1)<<endl;

    return 0;
}
