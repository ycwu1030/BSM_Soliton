#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "PathDeformation.h"
#include "ToyZ4.h"

using namespace std;

int main(int argc, char const *argv[]) {
    ToyZ4 model;
    model.Set_Potential_Parameters(0.5, 0.5);
    // model.PrintParameters();
    // model.PrintLocalMinima();

    VnD vtol = [&](VD field, double scale) { return model.Vtotal(field, scale); };
    dVnD dvtol = [&](VD field, double scale) { return model.dVtotal(field, scale); };

    ofstream out("ToyZ4_Scan.dat");
    out << "beta tension_adj width_adj tension_ops width_ops" << endl;
    out << std::scientific;
    out << std::showpos;
    out << std::setprecision(10);
    for (int rid = -40; rid < 0; rid += 1) {
        double beta = pow(10, rid / 10.0);
        double delta_beta = 1.0 - beta;
        double r = 2 * beta / delta_beta;
        model.Set_Potential_Parameters(beta, delta_beta);
        VVD pts_init;
        VD left = model.GetLocalMinima(0);
        VD right = model.GetLocalMinima(1);
        double v1 = 1;
        double v2 = 1.0 / sqrt(2 + 2 * r);
        pts_init.push_back({1, 0});
        pts_init.push_back({0, 1});
        cout << "Adj beta = " << beta << endl;
        KinknD sol_adj = fullKink(pts_init, vtol, dvtol);
        char fname[200];
        sprintf(fname, "ToyZ4_Kink/ToyZ4_%d_adj.dat", rid);
        model.DumpFullSolution(sol_adj.R, sol_adj.Phi, sol_adj.dPhi_1D, fname);
        pts_init.clear();
        pts_init.push_back({1, 0});
        pts_init.push_back({-1, 0});
        cout << "Ops beta = " << beta << endl;
        KinknD sol_ops = fullKink(pts_init, vtol, dvtol);
        sprintf(fname, "ToyZ4_Kink/ToyZ4_%d_ops.dat", rid);
        model.DumpFullSolution(sol_ops.R, sol_ops.Phi, sol_ops.dPhi_1D, fname);
        cout << " Dump energy and wall width" << endl;
        out << beta << " " << model.GetTotalEnergy(sol_adj.R, sol_adj.Phi, sol_adj.dPhi_1D) << " "
            << model.GetWallWidth(sol_adj.R, sol_adj.Phi, sol_adj.dPhi_1D) << " "
            << model.GetTotalEnergy(sol_ops.R, sol_ops.Phi, sol_ops.dPhi_1D) << " "
            << model.GetWallWidth(sol_ops.R, sol_ops.Phi, sol_ops.dPhi_1D) << endl;
    }
    ofstream out1("ToyZ4_Scan_close_to_1.dat");
    out1 << "beta tension_adj width_adj tension_ops width_ops" << endl;
    out1 << std::scientific;
    out1 << std::showpos;
    out1 << std::setprecision(10);
    for (int rid = -40; rid < 0; rid += 1) {
        double delta_beta = pow(10, rid / 10.0);
        double beta = 1 - delta_beta;
        double r = 2 * beta / (1 - beta);
        model.Set_Potential_Parameters(beta, delta_beta);
        VVD pts_init;
        VD left = model.GetLocalMinima(0);
        VD right = model.GetLocalMinima(1);
        double v1 = 1;
        double v2 = 1.0 / sqrt(2 + 2 * r);
        pts_init.push_back({1, 0});
        pts_init.push_back({0, 0});
        pts_init.push_back({0, 1});
        KinknD sol_adj = fullKink(pts_init, vtol, dvtol, 50, 0.03, false, 400);
        char fname[200];
        sprintf(fname, "ToyZ4_Kink/ToyZ4_Close_One_%d_adj.dat", rid);
        model.DumpFullSolution(sol_adj.R, sol_adj.Phi, sol_adj.dPhi_1D, fname);
        pts_init.clear();
        pts_init.push_back({1, 0});
        pts_init.push_back({-1, 0});
        KinknD sol_ops = fullKink(pts_init, vtol, dvtol, 50, 0.03, false, 400);
        sprintf(fname, "ToyZ4_Kink/ToyZ4_Close_One_%d_ops.dat", rid);
        model.DumpFullSolution(sol_ops.R, sol_ops.Phi, sol_ops.dPhi_1D, fname);
        out1 << delta_beta << " " << model.GetTotalEnergy(sol_adj.R, sol_adj.Phi, sol_adj.dPhi_1D) << " "
             << model.GetWallWidth(sol_adj.R, sol_adj.Phi, sol_adj.dPhi_1D) << " "
             << model.GetTotalEnergy(sol_ops.R, sol_ops.Phi, sol_ops.dPhi_1D) << " "
             << model.GetWallWidth(sol_ops.R, sol_ops.Phi, sol_ops.dPhi_1D) << endl;
    }
    {
        double beta = 3.0 / 4.0;
        double delta_beta = 1.0 / 4.0;
        model.Set_Potential_Parameters(beta, delta_beta);
        VVD pts_init;
        // VD left = model.GetLocalMinima(0);
        // VD right = model.GetLocalMinima(1);
        pts_init.push_back({1, 0});
        pts_init.push_back({0, 1});
        KinknD sol_adj = fullKink(pts_init, vtol, dvtol);
        char fname[200];
        sprintf(fname, "ToyZ4_Kink/ToyZ4_3over4_adj.dat");
        model.DumpFullSolution(sol_adj.R, sol_adj.Phi, sol_adj.dPhi_1D, fname);
        pts_init.clear();
        pts_init.push_back({1, 0});
        pts_init.push_back({-1, 0});
        KinknD sol_ops = fullKink(pts_init, vtol, dvtol);
        sprintf(fname, "ToyZ4_Kink/ToyZ4_3over4_ops.dat");
        model.DumpFullSolution(sol_ops.R, sol_ops.Phi, sol_ops.dPhi_1D, fname);
        out << beta << " " << model.GetTotalEnergy(sol_adj.R, sol_adj.Phi, sol_adj.dPhi_1D) << " "
            << model.GetWallWidth(sol_adj.R, sol_adj.Phi, sol_adj.dPhi_1D) << " "
            << model.GetTotalEnergy(sol_ops.R, sol_ops.Phi, sol_ops.dPhi_1D) << " "
            << model.GetWallWidth(sol_ops.R, sol_ops.Phi, sol_ops.dPhi_1D) << endl;
    }
    {
        double beta = 1.0 / 3.0;
        double delta_beta = 2.0 / 3.0;
        model.Set_Potential_Parameters(beta, delta_beta);
        VVD pts_init;
        // VD left = model.GetLocalMinima(0);
        // VD right = model.GetLocalMinima(1);
        pts_init.push_back({1, 0});
        pts_init.push_back({0, 1});
        KinknD sol_adj = fullKink(pts_init, vtol, dvtol);
        char fname[200];
        sprintf(fname, "ToyZ4_Kink/ToyZ4_1over3_adj.dat");
        model.DumpFullSolution(sol_adj.R, sol_adj.Phi, sol_adj.dPhi_1D, fname);
        pts_init.clear();
        pts_init.push_back({1, 0});
        pts_init.push_back({-1, 0});
        KinknD sol_ops = fullKink(pts_init, vtol, dvtol);
        sprintf(fname, "ToyZ4_Kink/ToyZ4_1over3_ops.dat");
        model.DumpFullSolution(sol_ops.R, sol_ops.Phi, sol_ops.dPhi_1D, fname);
        out << beta << " " << model.GetTotalEnergy(sol_adj.R, sol_adj.Phi, sol_adj.dPhi_1D) << " "
            << model.GetWallWidth(sol_adj.R, sol_adj.Phi, sol_adj.dPhi_1D) << " "
            << model.GetTotalEnergy(sol_ops.R, sol_ops.Phi, sol_ops.dPhi_1D) << " "
            << model.GetWallWidth(sol_ops.R, sol_ops.Phi, sol_ops.dPhi_1D) << endl;
    }
    return 0;
}
