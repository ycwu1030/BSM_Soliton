#include "cxSM_Z2.h"
#include <iostream>
using namespace std;

int main(int argc, char const *argv[])
{
    cxSM_Z2 model;
    model.Set_Physical_Parameters(10,0.1,200);
    model.PrintParameters();
    model.PrintLocalMinima();
    cout<<"Stability: "<<model.CheckStability()<<endl;
    cout<<"Unitarity: "<<model.CheckUnitarity()<<endl;
    cout<<"Global Min: "<<model.CheckGlobalMinimum()<<endl;
    return 0;
}
