#include "ToyZ4.h"

#include <cmath>
#include <iostream>

#include "DWSolver.h"

using namespace std;

namespace ToyTest {
ToyZ4::ToyZ4() : BSM_Soliton::BaseModel(2) {}

void ToyZ4::Set_Potential_Parameters(double beta_, double delta_beta_) {
    beta = beta_;
    delta_beta = delta_beta_;
    r = 2 * beta / delta_beta;
    Solved = false;
    Find_Local_Extrema();

    _v1global = 0;
    _v2global = 0;
    _vr = 0;
    Input_Minimum_ID = -1;
    double Vtmp = 0;
    for (size_t i = 0; i < N_Local_Extrema; i++) {
        if (is_Local_Minima[i]) {
            double vtt = sqrt(pow(Local_Extrema[i][0], 2) + pow(Local_Extrema[i][1], 2));
            if (vtt > _vr) {
                _vr = vtt;
            }
        }
        if (Potential_at_Extrema[i] < Vtmp && is_Local_Minima[i]) {
            _v2global = Local_Extrema[i][1];
            _v1global = Local_Extrema[i][0];
            Input_Minimum_ID = i;
            Vtmp = Potential_at_Extrema[i];
        }
    }
}

double ToyZ4::V(VD field_values, double scale) {
    double p1 = field_values[0];
    double p2 = field_values[1];
    double pt = p1 * p1 + p2 * p2;
    double vtot = -pt / 2 + pt * pt / 4 + r * p1 * p1 * p2 * p2;
    return vtot / pow(scale, 4);
}

VD ToyZ4::dV(VD field_values, double scale) {
    VD res(2);
    double p1 = field_values[0];
    double p2 = field_values[1];
    res[0] = p1 * (p1 * p1 + p2 * p2) + 2 * r * p1 * p2 * p2 - p1;
    res[1] = p2 * (p1 * p1 + p2 * p2) + 2 * r * p1 * p1 * p2 - p2;
    return res / pow(scale, 3);
}

VVD ToyZ4::d2V(VD field_values, double scale) {
    VVD res(2, VD(2));
    double p1 = field_values[0];
    double p2 = field_values[1];
    res[0][0] = 3 * p1 * p1 + (1 + 2 * r) * p2 * p2 - 1;
    res[0][1] = 2 * (1 + 2 * r) * p1 * p2;
    res[1][0] = res[0][1];
    res[1][1] = (1 + 2 * r) * p1 * p1 + 3 * p2 * p2 - 1;

    return res / pow(scale, 2);
}

double ToyZ4::V_min(double scale) { return V({_v1global, _v2global}, scale); }

VD ToyZ4::Quartic_Couplings(VD field_values) { return {1.0 / 4.0, 1.0 / 4.0}; }

void ToyZ4::Calculate_Local_Extrema() {
    Add_Local_Extremum({0, 0});

    double r1 = 1.0;
    Add_Local_Extremum({r1, 0});
    Add_Local_Extremum({-r1, 0});
    Add_Local_Extremum({0, r1});
    Add_Local_Extremum({0, -r1});

    double t2_tmp = 2 + 2 * r;
    if (t2_tmp > 0) {
        double r2 = 1.0 / sqrt(2 + 2 * r);

        Add_Local_Extremum({r2, r2});
        Add_Local_Extremum({r2, -r2});
        Add_Local_Extremum({-r2, r2});
        Add_Local_Extremum({-r2, -r2});
    }
}

void ToyZ4::Print_Parameters() {
    cout << "Potential Parameter:" << endl;
    cout << "r:\t" << r << endl;
}
}  // namespace ToyTest

int main(int argc, char const *argv[]) {
    ToyTest::ToyZ4 model;
    model.Set_Potential_Parameters(0.99, 0.01);
    BSM_Soliton::DomainWallSolver DW(&model);
    DW.Solve({1, 0}, {0, 1});
    double tension = DW.Get_Tension();
    double dw_width = DW.Get_Wall_Width();
    DW.Dump_Solution("ToyZ4_Relaxation_test.dat");
    return 0;
}
