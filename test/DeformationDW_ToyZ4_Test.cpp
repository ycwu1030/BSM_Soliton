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

    ofstream out("ToyZ4_Kink_solution_test.dat");
    VVD pts_init;
    pts_init.push_back({1, 0});
    pts_init.push_back({0, 1});
    KinknD sol_adj = fullKink(pts_init, vtol, dvtol);
    out << std::scientific;
    out << std::showpos;
    out << std::setprecision(8);
    for (size_t i = 0; i < sol_adj.R.size(); i++) {
        out << sol_adj.R[i] << "\t" << sol_adj.Phi_1D[i] << "\t" << sol_adj.dPhi_1D[i] << "\t" << sol_adj.Phi[i]
            << endl;
    }

    return 0;
}
