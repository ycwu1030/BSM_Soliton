#include "PathDeformation.h"
// #include "Relaxation.h"
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

#define lamh (0.3)
#define lams (0.4)
#define lamhs (0.25)
#define vh (246.221)
#define vs (10000.0)

double Potential(VD X, double scale)
{
    double h = X[0];
    double s = X[1];
    double pot = 0;
    pot += (-2*lamh*vh*vh - lamhs*vs*vs)*h*h;
    pot += (-2*lams*vs*vs - lamhs*vh*vh)*s*s;
    pot += lamh*h*h*h*h;
    pot += lams*s*s*s*s;
    pot += lamhs*h*h*s*s;
    return pot;
}
VD dPotential(VD X, double scale)
{
    double h = X[0];
    double s = X[1];
    VD res(2);
    res[0] = 4*lamh*(h*h-vh*vh)*h + 2*lamhs*(s*s-vs*vs)*h;
    res[1] = 4*lams*(s*s-vs*vs)*s + 2*lamhs*(h*h-vh*vh)*s;
    return res;
}
int main(int argc, char const *argv[])
{
    VVD pts_init(2);
    pts_init[0] = {vh,-vs};
    pts_init[1] = {vh, vs};
    KinknD sol = fullKink(pts_init,Potential,dPotential);
    ofstream output("KinknD_test.dat");
    output<<"r\tphi1d\tdphi1d\tphi0\tphi1"<<endl;
    for (int i = 0; i < sol.R.size(); i++)
    {
        output<<sol.R[i]<<"\t"<<sol.Phi_1D[i]<<"\t"<<sol.dPhi_1D[i]<<"\t"<<sol.Phi[i]<<endl;
    }
    output.close();
    return 0;
}

