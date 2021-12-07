#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "DWSolver.h"
#include "ToyZ4.h"

using namespace std;

int main(int argc, char const *argv[]) {
    ToyZ4 model;
    model.Set_Potential_Parameters(0.5, 0.5);

    ofstream out("ToyZ4_Scan_Relax.dat");
    out << "beta tension_adj width_adj tension_ops width_ops" << endl;
    for (int rid = -40; rid < 0; rid += 1) {
        double beta = pow(10, rid / 10.0);
        double delta_beta = 1.0 - beta;
        double r = 2 * beta / delta_beta;
        model.Set_Potential_Parameters(beta, delta_beta);
        VD left = {1, 0};
        VD right = {0, 1};
        DWSolver sol_adj(&model, left, right);
        sol_adj.SetZRange(20);
        char fname[200];
        VD X_adj;
        VVD Y_adj;
        sol_adj.Solve(X_adj, Y_adj);
        sprintf(fname, "ToyZ4_Kink/ToyZ4_Relax_%d_adj.dat", rid);
        model.DumpFullSolution(X_adj, Y_adj, fname);

        out << beta << " " << model.GetTotalEnergy(X_adj, Y_adj) << " " << model.GetWallWidth(X_adj, Y_adj) << endl;
    }
    ofstream out1("ToyZ4_Scan_close_to_1_Relax.dat");
    out1 << "delta_beta tension_adj width_adj tension_ops width_ops" << endl;
    for (int rid = -40; rid < 0; rid += 1) {
        double delta_beta = pow(10, rid / 10.0);
        double beta = 1 - delta_beta;
        model.Set_Potential_Parameters(beta, delta_beta);
        VD left = {1, 0};
        VD right = {0, 1};
        DWSolver sol_adj(&model, left, right);
        sol_adj.SetZRange(20);
        char fname[200];
        VD X_adj;
        VVD Y_adj;
        sol_adj.Solve(X_adj, Y_adj);
        sprintf(fname, "ToyZ4_Kink/ToyZ4_Relax_Close_One_%d_adj.dat", rid);
        model.DumpFullSolution(X_adj, Y_adj, fname);

        out1 << delta_beta << " " << model.GetTotalEnergy(X_adj, Y_adj) << " " << model.GetWallWidth(X_adj, Y_adj)
             << endl;
    }
    {
        double beta = 3.0 / 4.0;
        double delta_beta = 1.0 / 4.0;
        model.Set_Potential_Parameters(beta, delta_beta);
        VD left = {1, 0};
        VD right = {0, 1};
        DWSolver sol_adj(&model, left, right);
        sol_adj.SetZRange(20);
        VD X_adj;
        VVD Y_adj;
        sol_adj.Solve(X_adj, Y_adj);
        char fname[200];
        sprintf(fname, "ToyZ4_Kink/ToyZ4_Relax_3over4_adj.dat");
        model.DumpFullSolution(X_adj, Y_adj, fname);
    }
    {
        double beta = 1.0 / 3.0;
        double delta_beta = 2.0 / 3.0;
        model.Set_Potential_Parameters(beta, delta_beta);
        VD left = {1, 0};
        VD right = {0, 1};
        DWSolver sol_adj(&model, left, right);
        sol_adj.SetZRange(20);
        VD X_adj;
        VVD Y_adj;
        sol_adj.Solve(X_adj, Y_adj);
        char fname[200];
        sprintf(fname, "ToyZ4_Kink/ToyZ4_Relax_1over3_adj.dat");
        model.DumpFullSolution(X_adj, Y_adj, fname);
    }
    return 0;
}
