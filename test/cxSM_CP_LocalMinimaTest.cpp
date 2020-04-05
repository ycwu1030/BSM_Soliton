#include "cxSM_CP.h"
#include <iostream>
using namespace std;

int main(int argc, char const *argv[])
{
    cxSM_CP model;
    // model.Set_Physical_Parameters(10,10,200,66,0.01,0.02,0.03);
    model.Set_Physical_Parameters(31.961,42.5982,410.83,456.265,1.24781,0.407,-0.487092);
    model.PrintParameters();
    model.PrintLocalMinima();
    cout<<"Stability: "<<model.CheckStability()<<endl;
    cout<<"Unitarity: "<<model.CheckUnitarity()<<endl;
    cout<<"Global Min: "<<model.CheckGlobalMinimum()<<endl;
    return 0;
}
