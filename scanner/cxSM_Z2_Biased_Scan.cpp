#include "cxSM_Z2_biased.h"
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
    cxSM_Z2_biased model;
    DWSolver solver(&model);
    solver.SetOverallScale(model.GetVEV());
    VD left(2);
    VD right(2);
    double vs;
    double MHH;
    double MHA;
    double theta;
    double del1;
    double c1;
    double c2;
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
        MHH = RandomReal(0,500);
        MHA = RandomReal(0,200);
        theta = RandomReal(0,0.5);
        del1 = RandomReal(-0.1,0.1);
        c1 = RandomReal(-0.1,0.1);
        c2 = RandomReal(-0.1,0.1);
        model.Set_Physical_Parameters(vs,theta,MHH,MHA,del1,c1,c2);
        if (!model.CheckStability()||!model.CheckUnitarity()||!model.CheckGlobalMinimum()) continue;
        good = model.GetBiasedMirrorMinimum(left);
        if (!good) continue;
        right[0] = model.GetVEV();
        right[1] = vs;
        solver.SetBoundary(left,right);
        good = solver.Solve(X,Y);
        if (!good) continue;
        ++NGOT;
        output<<NGOT<<"\t"<<vs<<"\t"<<MHH<<"\t"<<MHA<<"\t"<<theta<<"\t"<<del1<<"\t"<<c1<<"\t"<<c2<<"\t"<<model.GetTotalEnergy(X,Y)<<"\t"<<model.GetTension(X,Y)<<endl;
    }
    // SS.PrintSolitonSolution();
    // SS.DumpSolitonSolution("solution_cxSM.dat");
    return 0;
}