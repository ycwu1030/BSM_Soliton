/*
 * @Description  : 
 * @Author       : Yongcheng Wu
 * @Date         : 2020-01-30 21:21:49
 * @LastEditors  : Yongcheng Wu
 * @LastEditTime : 2020-01-30 22:07:52
 */
#include "SM_cxSM.h"
#include "SolitonSolver.h"
#include <iostream>
#include <fstream>
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
double RandomReal(double min, double max)
{
    return (rand()/(double)RAND_MAX)*(max-min)+min;
}
int main(int argc, char const *argv[])
{
    CXSM model;
    SolitonSolver SS(2,100);
    SS.SetPotentials(V,dV);
    SS.SetParam(&model);
    
    double vs;
    double m2;
    double theta;
    srand (time(NULL));
    int NGOT = 0;
    int NTRIED = 0;
    ofstream output(argv[1]);
    while (NGOT < 10000)
    {
        ++NTRIED;
        vs = RandomReal(0,200);
        m2 = RandomReal(0,500);
        theta = RandomReal(0,0.5);
        model.SetInput(vs,m2,0,theta);
        if (!model.CheckStability()||!model.CheckUnitarity()||!model.CheckGlobalMinimum()) continue;
        ++NGOT;
        VD LEFT = {model.vev, -model.VS};
        VD RIGHT = {model.vev, model.VS};
        SS.SetBoundary(LEFT,RIGHT);
        SS.SetScale(model.vev);
        SS.SetV0Global(model.Get_V0_global());
        SS.Solve();
        output<<NGOT<<"\t"<<vs<<"\t"<<m2<<"\t"<<theta<<"\t"<<SS.GetTotalEnergy()<<"\t"<<SS.GetTension()<<endl;
    }
    // SS.PrintSolitonSolution();
    // SS.DumpSolitonSolution("solution_cxSM.dat");
    return 0;
}