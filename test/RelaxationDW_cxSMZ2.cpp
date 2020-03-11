#include "cxSM_Z2.h"
#include "DWSolver.h"
// #include "Relaxation.h"
#include <iostream>
#include <cmath>

using namespace std;


int main(int argc, char const *argv[])
{
    cxSM_Z2 model;
    // model.SetInput(157.633,42.647,0.335465);
    model.Set_Physical_Parameters(157.633,0.335465,42.647);
    VD left = {model.GetVEV(),-157.633};
    VD right = {model.GetVEV(),157.633};
    DWSolver solver(&model,left,right);
    solver.SetOverallScale(model.GetVEV());
    solver.SetXRange();
    VD X;
    VVD Y;
    bool good = solver.Solve(X,Y);
    for (size_t i = 0; i < X.size(); i++)
    {
        cout<<X[i]<<"\t";
        for (size_t j = 0; j < Y[i].size(); j++)
        {
            cout<<Y[i][j]<<"\t";
        }
        cout<<endl;
    }
    cout<<"Energy: "<<model.GetTotalEnergy(X,Y)<<endl;
    cout<<"Tension: "<<model.GetTension(X,Y)<<endl;
    // RX.DumpSolution("Relaxation_DW_test.dat");
    return 0;
}

