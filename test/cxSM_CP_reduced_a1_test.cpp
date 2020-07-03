#include "cxSM_CP_reduced_a1.h"
#include <random>
#include <iostream>
using namespace std;

double ScaleTo(double x, double rmin, double rmax)
{
    return rmin + x*(rmax-rmin);
}
int main(int argc, char const *argv[])
{
    cxSM_CP_reduced_a1 model;
    // bool notgot = true;
    bool good;
    double vs,alpha;
    double vsr,vsi,MHH,MHA,theta1,theta3;
    default_random_engine gen;
    uniform_real_distribution<double> uni(0,1);
    while (true)
    {
        vs = ScaleTo(uni(gen),0,100);
        MHH = ScaleTo(uni(gen),70,500);
        MHA = ScaleTo(uni(gen),0,500);
        theta1 = ScaleTo(uni(gen),0,M_PI_2);
        theta3 = ScaleTo(uni(gen),0,M_PI_2);
        good = model.Set_Physical_Parameters_vs_theta(vs,MHH,MHA,theta1,theta3);
        if (good)
        {
            break;
        }
    }
    model.PrintParameters();
    model.PrintLocalMinima();
    cout<<"Stability: "<<model.CheckStability()<<endl;
    cout<<"Unitarity: "<<model.CheckUnitarity()<<endl;
    cout<<"Global Min: "<<model.CheckGlobalMinimum()<<endl;


    while (true)
    {
        vs = ScaleTo(uni(gen),0,100);
        MHH = ScaleTo(uni(gen),70,500);
        MHA = ScaleTo(uni(gen),0,500);
        theta1 = ScaleTo(uni(gen),0,M_PI_2);
        theta3 = ScaleTo(uni(gen),0,M_PI_2);
        good = model.Set_Physical_Parameters_vsr_theta(vs,MHH,MHA,theta1,theta3);
        if (good)
        {
            break;
        }
    }
    model.PrintParameters();
    model.PrintLocalMinima();
    cout<<"Stability: "<<model.CheckStability()<<endl;
    cout<<"Unitarity: "<<model.CheckUnitarity()<<endl;
    cout<<"Global Min: "<<model.CheckGlobalMinimum()<<endl;
    return 0;
}
