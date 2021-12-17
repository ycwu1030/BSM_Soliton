#include "THDM_SCPV1.h"
#include <iostream>
using namespace std;

int main(int argc, char const *argv[])
{
    THDM_SCPV1 model;
    // model.Set_Physical_Parameters(10,10,200,66,0.01,0.02,0.03);
    bool good = model.Set_Physical_Parameters(atan(5.0),300,350,380,atan(5.0)-M_PI_2,0.1);
    if (!good)
    {
        cout<<"Bad Input!"<<endl;
        return -1;
    }
    model.PrintParameters();
    model.PrintLocalMinima();
    cout<<"Stability: "<<model.CheckStability()<<endl;
    cout<<"Unitarity: "<<model.CheckUnitarity()<<endl;
    cout<<"Global Min: "<<model.CheckGlobalMinimum()<<endl;
    return 0;
}
