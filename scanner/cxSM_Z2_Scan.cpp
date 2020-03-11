#include "cxSM_Z2.h"
#include "DWSolver.h"
#include <iostream>
#include <fstream>
using namespace std;

double RandomReal(double min, double max)
{
    return (rand()/(double)RAND_MAX)*(max-min)+min;
}
int main(int argc, char const *argv[])
{
    cxSM_Z2 model;
    DWSolver solver(&model);
    solver.SetOverallScale(model.GetVEV());
    VD left(2);
    VD right(2);
    double vs;
    double m2;
    double theta;
    srand (time(NULL));
    int NGOT = 0;
    int NTRIED = 0;
    ofstream output(argv[1]);
    VD X;
    VVD Y;
    bool good;
    while (NGOT < 10000)
    {
        ++NTRIED;
        vs = RandomReal(0,200);
        m2 = RandomReal(0,500);
        theta = RandomReal(0,0.5);
        model.Set_Physical_Parameters(vs,theta,m2);
        if (!model.CheckStability()||!model.CheckUnitarity()||!model.CheckGlobalMinimum()) continue;
        left[0] = model.GetVEV();
        left[1] = -vs;
        right[0] = model.GetVEV();
        right[1] = vs;
        solver.SetBoundary(left,right);
        good = solver.Solve(X,Y);
        if (!good) continue;
        ++NGOT;
        output<<NGOT<<"\t"<<vs<<"\t"<<m2<<"\t"<<theta<<"\t"<<model.GetTotalEnergy(X,Y)<<"\t"<<model.GetTension(X,Y)<<endl;
    }
    // SS.PrintSolitonSolution();
    // SS.DumpSolitonSolution("solution_cxSM.dat");
    return 0;
}