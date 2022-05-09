#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "PathDeformation.h"
#include "ToyZ6.h"

using namespace std;

int main(int argc, char const *argv[]) {
    ToyTest::ToyZ6 model;

    VnD vtol = [&](VD field, double scale) { return model.Vtotal(field, scale); };
    dVnD dvtol = [&](VD field, double scale) { return model.dVtotal(field, scale); };

    ofstream out("ToyZ6_Scan.dat");
    out << "beta tension_adj width_adj tension_oth width_oth tension_ops width_ops" << endl;
    out << std::scientific;
    out << std::showpos;
    out << std::setprecision(10);
    for (int rid = -30; rid < -10; rid += 1) {
        double beta = pow(10, rid / 10.0);
        model.Set_potential_Parameters(-beta);
        // model.PrintLocalMinima();
        model.PrintMinimaInOrder();
        VVD pts_init;
        if (model.vev_order_by_angle.size() < 6) continue;
        VD p1 = model.vev_order_by_angle[0];
        VD p2 = model.vev_order_by_angle[1];
        VD p3 = model.vev_order_by_angle[2];
        VD p4 = model.vev_order_by_angle[3];

        // * Adjcent
        pts_init.push_back(p1);
        pts_init.push_back(p2);
        KinknD sol_adj = fullKink(pts_init, vtol, dvtol);
        char fname[200];
        sprintf(fname, "ToyZ6_Kink/ToyZ6_%d_adj.dat", rid);
        model.DumpFullSolution(sol_adj.R, sol_adj.Phi, sol_adj.dPhi_1D, fname);

        // * Every the other
        pts_init.clear();
        pts_init.push_back(p1);
        pts_init.push_back(p3);
        KinknD sol_oth = fullKink(pts_init, vtol, dvtol);
        sprintf(fname, "ToyZ6_Kink/ToyZ6_%d_oth.dat", rid);
        model.DumpFullSolution(sol_oth.R, sol_oth.Phi, sol_oth.dPhi_1D, fname);

        // * Opposite
        pts_init.clear();
        pts_init.push_back(p1);
        pts_init.push_back({0, abs(p1[0]) / 4.0});
        pts_init.push_back(p4);
        cout << "From p1=(" << p1[0] << "," << p1[1] << ") to p4=(" << p4[0] << "," << p4[1] << ")" << endl;
        KinknD sol_ops = fullKink(pts_init, vtol, dvtol);
        sprintf(fname, "ToyZ6_Kink/ToyZ6_%d_ops.dat", rid);
        model.DumpFullSolution(sol_ops.R, sol_ops.Phi, sol_ops.dPhi_1D, fname);

        out << beta << " " << model.GetTotalEnergy(sol_adj.R, sol_adj.Phi, sol_adj.dPhi_1D) << " "
            << model.GetWallWidth(sol_adj.R, sol_adj.Phi, sol_adj.dPhi_1D) << " "
            << model.GetTotalEnergy(sol_oth.R, sol_oth.Phi, sol_oth.dPhi_1D) << " "
            << model.GetWallWidth(sol_oth.R, sol_oth.Phi, sol_oth.dPhi_1D) << " "
            << model.GetTotalEnergy(sol_ops.R, sol_ops.Phi, sol_ops.dPhi_1D) << " "
            << model.GetWallWidth(sol_ops.R, sol_ops.Phi, sol_ops.dPhi_1D) << endl;
    }

    return 0;
}
