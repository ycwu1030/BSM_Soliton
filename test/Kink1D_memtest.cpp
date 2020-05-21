#include "Kink1D.h"
// #include "Relaxation.h"
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

#define lam (0.5)
#define vev (20000.0)
#define vevleft (-179.0)
#define vevc (1.0)
#define vevright (180)

double pot_z2(double phi)
{
    return lam/4.0*pow(phi*phi-vev*vev,2);
}
double dpot_z2(double phi)
{
    return lam*phi*(phi*phi-vev*vev);
}
double d2pot_z2(double phi)
{
    return lam*(3*phi*phi-vev*vev);
}

double pot_bias(double phi)
{
    return lam/12.0*phi*(3*phi*phi*phi-4*phi*phi*(vevleft+vevc+vevright)+6*phi*(vevleft*(vevc+vevright)+vevc*vevright)-12*vevleft*vevc*vevright);
}
double dpot_bias(double phi)
{
    return lam*(phi-vevleft)*(phi-vevc)*(phi-vevright);
}
double d2pot_bias(double phi)
{
    return lam*((phi-vevleft)*(phi-vevc)+(phi-vevleft)*(phi-vevright)+(phi-vevc)*(phi-vevright));
}
int main(int argc, char const *argv[])
{
    // Kink1D sol(vevleft,vevright,pot_bias,dpot_bias,d2pot_bias,1e-8,true);
    int test = 0;
    while (test<=1000)
    {
        Kink1D sol(vev,-vev,pot_z2,dpot_z2,d2pot_z2,1e-8,false);
        VD R;
        VD phi;
        VD dphi;
        VD Rerr;
        tie(R,phi,dphi,Rerr)=sol.findProfile();
        ofstream output("Kink1D_z2_test.dat");
        output<<"r\tphi\tdphi"<<endl;
        for (int i = 0; i < R.size(); i++)
        {
            output<<R[i]<<"\t"<<phi[i]<<"\t"<<dphi[i]<<endl;
        }
        output.close();
        test++;
    }
    return 0;
}

