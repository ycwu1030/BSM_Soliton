#include "ToyA4.h"
#include "PathDeformation.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

using namespace std;

int main(int argc, char const *argv[])
{
    ToyA4 model;
    model.Set_Potential_Parameters(-1.0,3.0/2.0,4.0);
    // model.PrintParameters();
    // model.PrintLocalMinima();
    
    VnD vtol = [&](VD field, double scale){
        return model.Vtotal(field,scale);
    };
    dVnD dvtol = [&](VD field, double scale){
        return model.dVtotal(field,scale);
    };

    ofstream out("ToyD4_Scan.dat");
    out<<"g1 g2 tension"<<endl;
    for (double g1 = 0.1; g1 <= 5.0; g1+=0.1)
    {
        for (double g2 = 0.1; g2 <= 5.0; g2+=0.1)
        {
            model.Set_Potential_Parameters(-1.0,g1,g2);
            double v1 = model.GetV1();
            double v2 = model.GetV2();
            VVD pts_init;
            VD left = {v1,0};
            VD right = {0,v2};
            pts_init.push_back(left);
            pts_init.push_back(right);
            KinknD sol = fullKink(pts_init,vtol,dvtol);
            out<<g1<<" "<<g2<<" "<<model.GetTotalEnergy(sol.R,sol.Phi)<<endl;
            // cout<<"Energy: "<<model.GetTotalEnergy(sol.R,sol.Phi)<<endl;
        }
    }

    return 0;
}
