#include "cxSM_Z2_biased.h"
#include <iostream>
using namespace std;

int main(int argc, char const *argv[])
{
    cxSM_Z2_biased model;
    model.Set_Physical_Parameters_del1_c1_c2(50,0.1,200,10,0.1,0.1,0.1);
    model.PrintParameters();
    model.PrintLocalMinima();
    cout<<"Stability: "<<model.CheckStability()<<endl;
    cout<<"Unitarity: "<<model.CheckUnitarity()<<endl;
    cout<<"Global Min: "<<model.CheckGlobalMinimum()<<endl;
    VD x(2);
    x[0]=0;
    x[1]=1;
    VVD f(2,VD(2));
    VD ff(2);
    ff[0] = 246;
    ff[1] = 10;
    f.push_back(ff);
    ff[0] = 246;
    ff[1] = -10;
    f.push_back(ff);
    cout<<"Total Energy: "<<model.GetTotalEnergy(x,f)<<endl;
    return 0;
}
