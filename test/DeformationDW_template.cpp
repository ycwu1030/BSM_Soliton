#include "PathDeformation.h"
#include "Potential.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

using namespace std;

class SelfPotential: public Potential
{
public:
    SelfPotential(int field_dimension=2):Potential(field_dimension){}; // Please specify the field dimension
    ~SelfPotential(){};

    double Vtotal(VD field_values, double scale = 1) {
        // This is your potential
        // Following is just a random example, please change to your own potential
        double h = field_values[0];
        double s = field_values[1];
        double lamh = 0.3;
        double lams = 0.4;
        double lamhs = 0.25;
        double vh = 246.221;
        double vs = 1000.0;

        // Potential
        double pot = (-2*lamh*vh*vh - lamhs*vs*vs)*h*h + (-2*lams*vs*vs - lamhs*vh*vh)*s*s + lamh*pow(h,4) + lams*pow(s,4) + lamhs*pow(h*s,2);

        return pot/pow(scale,4);
    }
    double V0_global(double scale = 1) {
        // This is the potential value at the global minimum point
        double vh = 246.221;
        double vs = 1000.0;
        double V0 = Vtotal({vh,vs},scale);
        return V0;
    }
    VD dVtotal(VD field_values, double scale = 1) {
        // This is the derivative of the potential respect to fields
        double h = field_values[0];
        double s = field_values[1];
        double lamh = 0.3;
        double lams = 0.4;
        double lamhs = 0.25;
        double vh = 246.221;
        double vs = 1000.0;

        VD result(2);

        result[0] = 4*lamh*(h*h-vh*vh)*h + 2*lamhs*(s*s-vs*vs)*h; // This is  dV/dh
        result[1] = 4*lams*(s*s-vs*vs)*s + 2*lamhs*(h*h-vh*vh)*s; // This is  dV/ds
        return result/pow(scale,3);
    }

};

int main(int argc, char const *argv[])
{
    SelfPotential pot;
    double vh = 246.221;
    double vs = 1000.0;
    VD left = {vh,vs};  // This provides the left side boundary
    VD right = {vh,-vs}; // This provides the right side boundary
    VVD pts_init;
    pts_init.push_back(left);
    pts_init.push_back(right);
    VnD vtol = [&](VD field, double scale){
        return pot.Vtotal(field,scale);
    };
    dVnD dvtol = [&](VD field, double scale){
        return pot.dVtotal(field,scale);
    };
    KinknD sol = fullKink(pts_init,vtol,dvtol);
    pot.DumpFullSolution(sol.R,sol.Phi,"PDDW_Solution_full.dat");
    cout<<"Energy: "<<pot.GetTotalEnergy(sol.R,sol.Phi)<<endl;
    return 0;
}

