#include "SM_cxSM.h"
#include "SolitonSolver.h"
#include <iostream>
using namespace std;

double V(VD Fields, void *param)
{
    Toy *model = (Toy*)param;
    double phiH = Fields[0];
    return model->Vtot(phiH);
}
VD dV(VD Fields, void *param)
{
    Toy *model = (Toy*)param;
    double phiH = Fields[0];
    VD res;
    res.push_back(model->dVdphi(phiH));
    return res;
}

int main(int argc, char const *argv[])
{
    Toy model(1,1);

    SolitonSolver SS(1,1000);
    SS.SetPotentials(V,dV);
    SS.SetParam(&model);
    VD LEFT = {-1};
    VD RIGHT = {+1};
    SS.SetBoundary(LEFT,RIGHT);
    SS.SetScale(1);
    SS.SetV0Global(model.Get_V0_global());
    SS.Solve();
    SS.PrintSolitonSolution();
    SS.DumpSolitonSolution("solution_toy.dat");
    return 0;
}
