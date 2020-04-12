#include <cmath>
#include "PathDeformation.h"
#include <iostream>
#include <fstream>

using namespace std;

#define N 150

#define lamh (0.3)
#define lams (0.4)
#define lamhs (0.25)
#define vh (246.221)
#define vs (1000.0)

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
    const int n = N;

    // ! Generating the path linking two local minima
    // ! (vh, vs) -- (vh, -vs)

    ofstream fdata("Spline_Potential_Data.dat");
    fdata<<"phi0\tphi1\ts\tV"<<endl;
    VD min1 = {vh,-vs};
    VD min2 = {vh, vs};
    VVD pts;
    double R = sqrt(vh*vh+vs*vs);
    double theta0 = atan(vs/vh);
    double dist;
    double phi0,phi1;
    double theta;
    for (int i = 0; i < n; i++)
    {
        theta = -theta0 + i*2*theta0/(n-1);
        phi0 = R*cos(theta);
        phi1 = R*sin(theta);
        if (i==0)
        {
            phi0 = vh;
            phi1 = -vs;
        }

        if (i== n-1)
        {
            phi0 = vh;
            phi1 = vs;
        }
        pts.push_back({phi0,phi1});
        dist = R*(theta+theta0);
        fdata << phi0<<"\t"<<phi1<<"\t"<<dist<<"\t"<<Potential({phi0,phi1},1)<<endl;  
    }
    fdata.close();

    SplinePath spt(pts,Potential);
    double V,dV,d2V;
    VD pts_pre;
    ofstream fout("Spline_Potential_Inter.dat");
    fout<<"s\tV\tdV\td2V\tphi0\tphi1"<<endl;
    for (double s = -0.1*dist; s < dist*1.1; s+=dist/200.0)
    {
        pts_pre = spt.pts_at_dist(s);
        V = spt.V(s);
        dV = spt.dV(s);
        d2V = spt.d2V(s);
        fout<<s<<"\t"<<V<<"\t"<<dV<<"\t"<<d2V<<"\t"<<pts_pre<<endl;
    }
    fout.close();

    return 0;
}
