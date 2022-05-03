#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "PathDeformation.h"
#include "ToyZ5.h"

using namespace std;

int main(int argc, char const *argv[]) {
    ToyZ5 model;
    model.Set_Potential_Parameters(0.5);
    // model.PrintParameters();
    // model.PrintLocalMinima();

    VnD vtol = [&](VD field, double scale) { return model.Vtotal(field, scale); };
    dVnD dvtol = [&](VD field, double scale) { return model.dVtotal(field, scale); };

    ofstream out("ToyZ5_Scan.dat");
    out << "beta tension_adj width_adj tension_ops width_ops" << endl;
    out << std::scientific;
    out << std::showpos;
    out << std::setprecision(10);
    for (int rid = -30; rid < -10; rid += 1) {
        double beta = pow(10, rid / 10.0);
        double delta_beta = 1.0 - beta;
        model.Set_Potential_Parameters(-beta);
        model.PrintLocalMinima();
        VVD pts_init;
        VD p1 = model.vev_order_by_angle[0];
        VD p2 = model.vev_order_by_angle[1];
        VD p3 = model.vev_order_by_angle[2];
        pts_init.push_back(p1);
        pts_init.push_back(p2);
        cout << "Adj beta = " << beta << endl;
        cout << "From (" << p1[0] << "," << p1[1] << ") -> (" << p2[0] << "," << p2[1] << ")" << endl;
        KinknD sol_adj = fullKink(pts_init, vtol, dvtol);
        char fname[200];
        sprintf(fname, "ToyZ5_Kink/ToyZ5_%d_adj.dat", rid);
        model.DumpFullSolution(sol_adj.R, sol_adj.Phi, sol_adj.dPhi_1D, fname);
        pts_init.clear();
        pts_init.push_back(p1);
        pts_init.push_back(p3);
        cout << "Ops beta = " << beta << endl;
        cout << "From (" << p1[0] << "," << p1[1] << ") -> (" << p3[0] << "," << p3[1] << ")" << endl;
        KinknD sol_ops = fullKink(pts_init, vtol, dvtol);
        sprintf(fname, "ToyZ5_Kink/ToyZ5_%d_ops.dat", rid);
        model.DumpFullSolution(sol_ops.R, sol_ops.Phi, sol_ops.dPhi_1D, fname);
        cout << " Dump energy and wall width" << endl;
        out << beta << " " << model.GetTotalEnergy(sol_adj.R, sol_adj.Phi, sol_adj.dPhi_1D) << " "
            << model.GetWallWidth(sol_adj.R, sol_adj.Phi, sol_adj.dPhi_1D) << " "
            << model.GetTotalEnergy(sol_ops.R, sol_ops.Phi, sol_ops.dPhi_1D) << " "
            << model.GetWallWidth(sol_ops.R, sol_ops.Phi, sol_ops.dPhi_1D) << endl;
    }
    return 0;
}
