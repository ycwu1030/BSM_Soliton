#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "PathDeformation.h"
#include "ToyZ4.h"

using namespace std;

int main(int argc, char const *argv[]) {
    ToyTest::ToyZ4 model;
    model.Set_Potential_Parameters(0.5, 0.5);
    // model.PrintParameters();
    // model.PrintLocalMinima();

    VnD vtol = [&](VD field, double scale) { return model.Vtotal(field, scale); };
    dVnD dvtol = [&](VD field, double scale) { return model.dVtotal(field, scale); };

    {
        double beta = 1.0 / 4.0;
        double delta_beta = 3.0 / 4.0;
        model.Set_Potential_Parameters(beta, delta_beta);
        VVD pts_init;
        // VD left = model.GetLocalMinima(0);
        // VD right = model.GetLocalMinima(1);
        double v1 = 1.0 / sqrt(delta_beta);
        pts_init.push_back({v1, 0});
        pts_init.push_back({0, v1});
        KinknD sol_adj = fullKink(pts_init, vtol, dvtol);
        char fname[200];
        sprintf(fname, "ToyZ4_Alter_Kink/ToyZ4_1over4_adj.dat");
        model.DumpFullSolution(sol_adj.R, sol_adj.Phi, sol_adj.dPhi_1D, fname);
        pts_init.clear();
        pts_init.push_back({v1, 0});
        pts_init.push_back({-v1, 0});
        KinknD sol_ops = fullKink(pts_init, vtol, dvtol);
        sprintf(fname, "ToyZ4_Alter_Kink/ToyZ4_1over4_ops_straight.dat");
        model.DumpFullSolution(sol_ops.R, sol_ops.Phi, sol_ops.dPhi_1D, fname);

        pts_init.clear();
        pts_init.push_back({v1, 0});
        pts_init.push_back({0, v1 * 0.5});
        pts_init.push_back({-v1, 0});
        KinknD sol_ops_d = fullKink(pts_init, vtol, dvtol);
        sprintf(fname, "ToyZ4_Alter_Kink/ToyZ4_1over4_ops_detour.dat");
        model.DumpFullSolution(sol_ops_d.R, sol_ops_d.Phi, sol_ops_d.dPhi_1D, fname);
    }
    return 0;
}
