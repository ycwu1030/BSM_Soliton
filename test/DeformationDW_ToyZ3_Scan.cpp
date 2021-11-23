#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "PathDeformation.h"
#include "ToyZ3.h"

using namespace std;

int main(int argc, char const *argv[]) {
    ToyZ3 model;
    model.Set_Potential_Parameters(0.5);
    // model.PrintParameters();
    // model.PrintLocalMinima();

    VnD vtol = [&](VD field, double scale) { return model.Vtotal(field, scale); };
    dVnD dvtol = [&](VD field, double scale) { return model.dVtotal(field, scale); };

    ofstream out("ToyZ3_Scan.dat");
    out << "r tension width" << endl;
    for (int rid = 1; rid <= 50; rid += 1) {
        double r = rid / 10.0;
        model.Set_Potential_Parameters(r);
        VVD pts_init;
        VD left = model.GetLocalMinima(0);
        VD right = model.GetLocalMinima(1);
        double v1 = sqrt(r * r + 1) - r;
        double v2 = sqrt(r * r + 1) + r;
        pts_init.push_back({-v2, 0});
        pts_init.push_back({-v1 / 2, -v1 * sqrt(3) / 2});
        pts_init.push_back({v2 / 2, -v2 * sqrt(3) / 2});
        KinknD sol = fullKink(pts_init, vtol, dvtol);
        char fname[200];
        sprintf(fname, "ToyZ3_Kink/ToyZ3_%d.dat", rid);
        model.DumpFullSolution(sol.R, sol.Phi, fname);
        out << r << " " << model.GetTotalEnergy(sol.R, sol.Phi) << " " << model.GetWallWidth(sol.R, sol.Phi) << endl;
    }

    return 0;
}
