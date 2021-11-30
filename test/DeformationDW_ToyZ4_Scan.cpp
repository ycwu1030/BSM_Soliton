#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "PathDeformation.h"
#include "ToyZ4.h"

using namespace std;

int main(int argc, char const *argv[]) {
    ToyZ4 model;
    model.Set_Potential_Parameters(0.5);
    // model.PrintParameters();
    // model.PrintLocalMinima();

    VnD vtol = [&](VD field, double scale) { return model.Vtotal(field, scale); };
    dVnD dvtol = [&](VD field, double scale) { return model.dVtotal(field, scale); };

    ofstream out("ToyZ4_Scan.dat");
    out << "beta tension_adj width_adj tension_ops width_ops" << endl;
    for (int rid = -40; rid < 0; rid += 1) {
        double beta = pow(10, rid / 10.0);
        double r = 2 * beta / (1 - beta);
        model.Set_Potential_Parameters(beta);
        VVD pts_init;
        VD left = model.GetLocalMinima(0);
        VD right = model.GetLocalMinima(1);
        double v1 = 1;
        double v2 = 1.0 / sqrt(2 + 2 * r);
        pts_init.push_back({1, 0});
        pts_init.push_back({0, 1});
        KinknD sol_adj = fullKink(pts_init, vtol, dvtol);
        char fname[200];
        sprintf(fname, "ToyZ4_Kink/ToyZ4_%d_adj.dat", rid);
        model.DumpFullSolution(sol_adj.R, sol_adj.Phi, fname);
        pts_init.clear();
        pts_init.push_back({1, 0});
        pts_init.push_back({-1, 0});
        KinknD sol_ops = fullKink(pts_init, vtol, dvtol);
        sprintf(fname, "ToyZ4_Kink/ToyZ4_%d_ops.dat", rid);
        model.DumpFullSolution(sol_ops.R, sol_ops.Phi, fname);
        out << beta << " " << model.GetTotalEnergy(sol_adj.R, sol_adj.Phi) << " "
            << model.GetWallWidth(sol_adj.R, sol_adj.Phi) << " " << model.GetTotalEnergy(sol_ops.R, sol_ops.Phi) << " "
            << model.GetWallWidth(sol_ops.R, sol_ops.Phi) << endl;
    }
    ofstream out1("ToyZ4_Scan_close_to_1.dat");
    out1 << "beta tension_adj width_adj tension_ops width_ops" << endl;
    for (int rid = -90; rid < 0; rid += 1) {
        double beta = pow(10, rid / 1000.0);
        double r = 2 * beta / (1 - beta);
        model.Set_Potential_Parameters(beta);
        VVD pts_init;
        VD left = model.GetLocalMinima(0);
        VD right = model.GetLocalMinima(1);
        double v1 = 1;
        double v2 = 1.0 / sqrt(2 + 2 * r);
        pts_init.push_back({1, 0});
        pts_init.push_back({0, 1});
        KinknD sol_adj = fullKink(pts_init, vtol, dvtol);
        char fname[200];
        sprintf(fname, "ToyZ4_Kink/ToyZ4_%d_10_adj.dat", rid);
        model.DumpFullSolution(sol_adj.R, sol_adj.Phi, fname);
        pts_init.clear();
        pts_init.push_back({1, 0});
        pts_init.push_back({-1, 0});
        KinknD sol_ops = fullKink(pts_init, vtol, dvtol);
        sprintf(fname, "ToyZ4_Kink/ToyZ4_%d_10_ops.dat", rid);
        model.DumpFullSolution(sol_ops.R, sol_ops.Phi, fname);
        out1 << beta << " " << model.GetTotalEnergy(sol_adj.R, sol_adj.Phi) << " "
             << model.GetWallWidth(sol_adj.R, sol_adj.Phi) << " " << model.GetTotalEnergy(sol_ops.R, sol_ops.Phi) << " "
             << model.GetWallWidth(sol_ops.R, sol_ops.Phi) << endl;
    }
    double beta = 3.0 / 4.0;
    model.Set_Potential_Parameters(beta);
    VVD pts_init;
    // VD left = model.GetLocalMinima(0);
    // VD right = model.GetLocalMinima(1);
    pts_init.push_back({1, 0});
    pts_init.push_back({0, 1});
    KinknD sol = fullKink(pts_init, vtol, dvtol);
    char fname[200];
    sprintf(fname, "ToyZ4_Kink/ToyZ4_3over4_adj.dat");
    model.DumpFullSolution(sol.R, sol.Phi, fname);
    return 0;
}
