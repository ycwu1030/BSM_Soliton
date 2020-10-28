#include "ToyD4.h"
#include "PathDeformation.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

using namespace std;

int main(int argc, char const *argv[])
{
    ToyD4 model;
    model.Set_Potential_Parameters(1.0);
    model.PrintParameters();
    model.PrintLocalMinima();

    double f1 = model.Getf1();
    double f2 = model.Getf2();

    
    VnD vtol = [&](VD field, double scale){
        return model.Vtotal(field,scale);
    };
    dVnD dvtol = [&](VD field, double scale){
        return model.dVtotal(field,scale);
    };
    VVD pts_init;
    VD left = {f1,0};
    VD right = {0,f2};
    pts_init.push_back(left);
    pts_init.push_back(right);
    KinknD sol = fullKink(pts_init,vtol,dvtol);

    model.DumpFullSolution(sol.R,sol.Phi,"ToyD4_1.dat");
    cout<<"Energy: "<<model.GetTotalEnergy(sol.R,sol.Phi)<<endl;
    cout<<"Tension: "<<model.GetTension(sol.R,sol.Phi)<<endl;

    return 0;
}
