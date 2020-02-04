/*
 * @Description  : 
 * @Author       : Yongcheng Wu
 * @Date         : 2020-01-27 17:56:36
 * @LastEditors  : Yongcheng Wu
 * @LastEditTime : 2020-02-04 12:15:44
 */
#include "SM_cxSM.h"
#include "SolitonSolver.h"
#include <iostream>
using namespace std;

double V(VD Fields, void *param)
{
    CXSM *model = (CXSM*)param;
    double phiH = Fields[0];
    double phiS = Fields[1];
    return model->Vtot(phiH,phiS,0);
}
VD dV(VD Fields, void *param)
{
    CXSM *model = (CXSM*)param;
    double phiH = Fields[0];
    double phiS = Fields[1];
    VD res;
    res.push_back(model->dVdH(phiH,phiS,0));
    res.push_back(model->dVdS(phiH,phiS,0));
    return res;
}

int main(int argc, char const *argv[])
{
    CXSM model;
    model.SetInput(157.633,42.647,0.0,0.335465);
    model.PrintPotentialParameter();
    model.PrintLocalMinima();
    cout<<"Stability: "<<model.CheckStability()<<endl;
    cout<<"Unitarity: "<<model.CheckUnitarity()<<endl;
    cout<<"Global Min: "<<model.CheckGlobalMinimum()<<endl;

    SolitonSolver SS(2,500);
    SS.SetPotentials(V,dV);
    SS.SetParam(&model);
    VD LEFT = {model.vev, -model.VS};
    VD RIGHT = {model.vev, model.VS};
    SS.SetBoundary(LEFT,RIGHT);
    SS.SetScale(model.vev);
    SS.SetV0Global(model.Get_V0_global());
    SS.Solve();
    SS.PrintSolitonSolution();
    SS.DumpSolitonSolution("solution_cxSM.dat");
    return 0;
}
