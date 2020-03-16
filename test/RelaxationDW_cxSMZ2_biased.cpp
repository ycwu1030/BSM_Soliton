#include "cxSM_Z2_biased.h"
#include "DWSolver.h"
// #include "Relaxation.h"
#include <iostream>
#include <cmath>

using namespace std;


int main(int argc, char const *argv[])
{
    cxSM_Z2_biased model;
    DWSolver solver(&model);
    bool good;
    // model.SetInput(157.633,42.647,0.335465);
    // model.Set_Physical_Parameters(157.633,0.335465,42.647);
    model.Set_Physical_Parameters_a1_c1_c2(500,0.1,400,0,0,0,0);
    model.PrintLocalMinima();
    if (!model.CheckStability()||!model.CheckUnitarity()||!model.CheckGlobalMinimum())
    {
        cout<<"Theoretical not good"<<endl;
        cout<<"Stability: "<<model.CheckStability()<<endl;
        cout<<"Unitarity: "<<model.CheckUnitarity()<<endl;
        cout<<"Global Min: "<<model.CheckGlobalMinimum()<<endl;
        return -1;
    }
    VD left(2);
    VD right(2);
    good = model.GetBiasedMirrorMinimum(left);
    if (!good) return -1;
    right[0] = model.GetVEV();
    right[1] = 500;
    solver.SetBoundary(left,right);
    solver.SetOverallScale(model.GetVEV());
    solver.SetXRange();
    VD X;
    VVD Y;
    good = solver.Solve(X,Y);
    // solver.PrintSolution();
    solver.DumpSolution("cxSMZ2_biased_DW_sol.dat");
    cout<<"Energy: "<<model.GetTotalEnergy(X,Y)<<endl;
    cout<<"Tension: "<<model.GetTension(X,Y)<<endl;
    // RX.DumpSolution("Relaxation_DW_test.dat");
    return 0;
}

