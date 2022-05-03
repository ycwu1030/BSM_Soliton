#include "ToyZ5.h"

#include <algorithm>
#include <cmath>

using namespace std;

double get_angle(vector<double> a) {
    double angle = atan2(a[1], a[0]);
    if (angle < 0) angle = angle + 2 * M_PI;
    return angle;
}

bool compare_in_angle(vector<double> a1, vector<double> a2) { return get_angle(a1) < get_angle(a2); }

ToyZ5::ToyZ5() : Basic_Model(2) {}

double ToyZ5::Vtotal(VD field_values, double scale) {
    double p1 = field_values[0];
    double p2 = field_values[1];
    double pt = p1 * p1 + p2 * p2;
    double vtot = -pt / 2 + pt * pt / 4 + beta_ * (pow(p1, 5) - 10 * pow(p1, 3) * pow(p2, 2) + 5 * p1 * pow(p2, 4));
    return vtot / pow(scale, 4);
}
void ToyZ5::Set_Potential_Parameters(double beta) {
    beta_ = beta;
    _Solved = false;
    FindLocalMinima();

    v1_global_ = 0;
    v2_global_ = 0;
    _IndexInput = -1;
    vev_order_by_angle.clear();
    double Vtmp = 0;
    for (size_t i = 0; i < _NLocalExtreme; i++) {
        if (_LocalMinimaQ[i]) {
            vev_order_by_angle.push_back({_localExtreme[i][0], _localExtreme[i][1]});
        }
        if (_Vtotal[i] < Vtmp && _LocalMinimaQ[i]) {
            v1_global_ = _localExtreme[i][0];
            v2_global_ = _localExtreme[i][1];
            _IndexInput = i;
            Vtmp = _Vtotal[i];
        }
    }
    sort(vev_order_by_angle.begin(), vev_order_by_angle.end(), compare_in_angle);
}

VD ToyZ5::dVtotal(VD field_values, double scale) {
    VD res(2);
    double p1 = field_values[0];
    double p2 = field_values[1];
    res[0] = p1 * (p1 * p1 + p2 * p2 - 1) + 5 * beta_ * (pow(p1, 4) - 6 * p1 * p1 * p2 * p2 + pow(p2, 4));
    res[1] = p2 * (p1 * p1 + p2 * p2 - 1) + 20 * beta_ * p1 * p2 * (p2 * p2 - p1 * p1);
    return res / pow(scale, 3);
}

VVD ToyZ5::d2Vtotal(VD field_values, double scale) {
    VVD res(2, VD(2));
    double p1 = field_values[0];
    double p2 = field_values[1];
    res[0][0] = 3 * p1 * p1 + p2 * p2 - 1 + 20 * beta_ * (pow(p1, 3) - 3 * p1 * p2 * p2);
    res[0][1] = 2 * p1 * p2 + 20 * beta_ * p2 * (p2 * p2 - 3 * p1 * p1);
    res[1][0] = res[0][1];
    res[1][1] = p1 * p1 + 3 * p2 * p2 - 1 - 20 * beta_ * p1 * (p1 * p1 - 3 * p2 * p2);
    return res / pow(scale, 2);
}

double ToyZ5::V0_global(double scale) { return Vtotal({v1_global_, v2_global_}, scale); }

VD ToyZ5::QuarticCoupling(VD field_values) { return {1.0 / 4.0, 1.0 / 4.0}; }

void ToyZ5::FindLocalMinima() {
    if (_Solved) {
        return;
    }

    Clear_Local_Cache();

    _localExtreme.push_back({0, 0});
    AppendLocalExtreme();

    // * We first solve the minima condition in polar system
    // * In the angluar direction
    // * We will always have following solution
    int N_ANGULAR = 5;
    double angulars_even[5] = {0 * M_PI / 5.0, 2 * M_PI / 5.0, 4 * M_PI / 5.0, 6 * M_PI / 5.0, 8 * M_PI / 5.0};
    sol_.Solve(5 * beta_, 1, 0, -1);
    auto res_even = sol_.GetRealSolution();
    for (auto &res_ : res_even) {
        // cout << "Even solution: " << res_ << endl;
        // if (res_ < 0) continue;
        for (int i = 0; i < N_ANGULAR; i++) {
            _localExtreme.push_back({res_ * cos(angulars_even[i]), res_ * sin(angulars_even[i])});
            AppendLocalExtreme();
        }
    }

    // double angulars_odd[5] = {1 * M_PI / 5.0, 3 * M_PI / 5.0, 5 * M_PI / 5.0, 7 * M_PI / 5.0, 9 * M_PI / 5.0};
    // sol_.Solve(-5 * beta_, 1, 0, -1);
    // auto res_odd = sol_.GetRealSolution();
    // for (auto &res_ : res_odd) {
    //     cout << "Odd solution: " << res_ << endl;
    //     for (int i = 0; i < N_ANGULAR; i++) {
    //         _localExtreme.push_back({res_ * cos(angulars_odd[i]), res_ * sin(angulars_odd[i])});
    //         AppendLocalExtreme();
    //     }
    // }
    _Solved = true;
}

void ToyZ5::PrintParameters() {
    cout << "Potential Parameter:" << endl;
    cout << "beta:\t" << beta_ << endl;
}
