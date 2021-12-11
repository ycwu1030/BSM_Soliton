#include "ToyA4.h"
#include "PathDeformation.h"
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

    
    VnD vtol = [&](VD field, double scale){
        return model.Vtotal(field,scale);
    };
    dVnD dvtol = [&](VD field, double scale){
        return model.dVtotal(field,scale);
    };
    VVD pts_init;
    VD left = {v1,0};
    VD right = {0,v2};
    pts_init.push_back(left);
    pts_init.push_back(right);
    KinknD sol = fullKink(pts_init,vtol,dvtol);

    model.DumpFullSolution(sol.R,sol.Phi,"ToyA4_1.dat");
    cout<<"Energy: "<<model.GetTotalEnergy(sol.R,sol.Phi)<<endl;
    cout<<"Tension: "<<model.GetTension(sol.R,sol.Phi)<<endl;

    VVD pts_init1;
    VD left1 = {v1,0};
    VD mid1 = {0,v2/3};
    VD right1 = {-v1,0};
    pts_init1.push_back(left1);
    pts_init1.push_back(mid1);
    pts_init1.push_back(right1);
    KinknD sol1 = fullKink(pts_init1,vtol,dvtol);

    model.DumpFullSolution(sol1.R,sol1.Phi,"ToyA4_2.dat");
    cout<<"Energy: "<<model.GetTotalEnergy(sol1.R,sol1.Phi)<<endl;
    cout<<"Tension: "<<model.GetTension(sol1.R,sol1.Phi)<<endl;

    // VVD pts_init2;
    // VD left2 = {v1,0.01};
    // VD right2 = {-v1,0.01};
    // pts_init2.push_back(left2);
    // pts_init2.push_back(right2);
    // KinknD sol2 = fullKink(pts_init2,vtol,dvtol);

    // model.DumpFullSolution(sol2.R,sol2.Phi,"ToyA4_3.dat");
    // cout<<"Energy: "<<model.GetTotalEnergy(sol2.R,sol2.Phi)<<endl;
    // cout<<"Tension: "<<model.GetTension(sol2.R,sol2.Phi)<<endl;

    return 0;
}
