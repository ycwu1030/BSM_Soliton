#include "ToyZ6.h"

#include <cmath>
#include <iostream>

using namespace std;

double get_angle(vector<double> a) {
    double angle = atan2(a[1], a[0]);
    if (angle < 0) angle = angle + 2 * M_PI;
    return angle;
}

bool compare_in_angle(vector<double> a1, vector<double> a2) { return get_angle(a1) < get_angle(a2); }

namespace ToyTest {
ToyZ6::ToyZ6() : Basic_Model(2) {}

void ToyZ6::Set_potential_Parameters(double beta) {
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

double ToyZ6::Vtotal(VD field_values, double scale) {
    double p1 = field_values[0];
    double p2 = field_values[1];
    double pt = p1 * p1 + p2 * p2;
    double vtot = -pt / 2 + pt * pt / 4 +
                  beta_ * (pow(p1, 6) - 15 * pow(p1, 4) * pow(p2, 2) + 15 * pow(p1, 2) * pow(p2, 4) - pow(p2, 6));
    return vtot / pow(scale, 4);
}

VD ToyZ6::dVtotal(VD field_values, double scale) {
    VD res(2);
    double p1 = field_values[0];
    double p2 = field_values[1];
    res[0] =
        p1 * (p1 * p1 + p2 * p2 - 1) + 6 * beta_ * (pow(p1, 5) - 10 * pow(p1, 3) * pow(p2, 2) + 5 * p1 * pow(p2, 4));
    res[1] =
        p2 * (p1 * p1 + p2 * p2 - 1) - 6 * beta_ * (pow(p2, 5) - 10 * pow(p1, 2) * pow(p2, 3) + 5 * pow(p1, 4) * p2);
    return res / pow(scale, 3);
}

VVD ToyZ6::d2Vtotal(VD field_values, double scale) {
    VVD res(2, VD(2));
    double p1 = field_values[0];
    double p2 = field_values[1];
    res[0][0] = 3 * p1 * p1 + p2 * p2 - 1 + 30 * beta_ * (pow(p1, 4) - 6 * p1 * p1 * p2 * p2 + pow(p2, 4));
    res[0][1] = 2 * p1 * p2 + 120 * beta_ * p1 * p2 * (p2 * p2 - p1 * p1);
    res[1][0] = res[0][1];
    res[1][1] = 3 * p2 * p2 + p1 * p1 - 1 - 30 * beta_ * (pow(p1, 4) - 6 * p1 * p1 * p2 * p2 + pow(p2, 4));
    return res / pow(scale, 2);
}

double ToyZ6::V0_global(double scale) { return Vtotal({v1_global_, v2_global_}, scale); }

VD ToyZ6::QuarticCoupling(VD field_values) { return {1.0 / 4.0, 1.0 / 4.0}; }

void ToyZ6::FindLocalMinima() {
    if (_Solved) {
        return;
    }

    Clear_Local_Cache();

    _localExtreme.push_back({0, 0});
    AppendLocalExtreme();

    static const int N_ANGULAR = 6;
    double angulars_even[N_ANGULAR] = {0 * M_PI / 6.0, 2 * M_PI / 6.0, 4 * M_PI / 6.0,
                                       6 * M_PI / 6.0, 8 * M_PI / 6.0, 10 * M_PI / 6.0};
    vector<double> res_even;
    if (24 * beta_ + 1 >= 0) {
        double tmp = sqrt(24 * beta_ + 1);
        if ((tmp - 1) / beta_ >= 0) {
            res_even.push_back(sqrt((tmp - 1) / beta_) / 2.0 / sqrt(3));
        }
        if (-(tmp + 1) / beta_ >= 0) {
            res_even.push_back(sqrt(-(tmp + 1) / beta_) / 2.0 / sqrt(3));
        }
    }
    for (auto &res_ : res_even) {
        for (int i = 0; i < N_ANGULAR; i++) {
            _localExtreme.push_back({res_ * cos(angulars_even[i]), res_ * sin(angulars_even[i])});
            AppendLocalExtreme();
        }
    }

    double angulars_odd[N_ANGULAR] = {1 * M_PI / 6.0, 3 * M_PI / 6.0, 5 * M_PI / 6.0,
                                      7 * M_PI / 6.0, 9 * M_PI / 6.0, 11 * M_PI / 6.0};
    vector<double> res_odd;
    if (1 - 24 * beta_ >= 0) {
        double tmp = sqrt(1 - 24 * beta_);
        if ((tmp + 1) / beta_ >= 0) {
            res_odd.push_back(sqrt((tmp + 1) / beta_) / 2.0 / sqrt(3));
        }
        if ((1 - tmp) / beta_ >= 0) {
            res_odd.push_back(sqrt((1 - tmp) / beta_) / 2.0 / sqrt(3));
        }
    }

    for (auto &res_ : res_odd) {
        for (int i = 0; i < N_ANGULAR; i++) {
            _localExtreme.push_back({res_ * cos(angulars_odd[i]), res_ * sin(angulars_odd[i])});
            AppendLocalExtreme();
        }
    }

    _Solved = true;
}

void ToyZ6::PrintParameters() {
    cout << "Potential Parameter:" << endl;
    cout << "beta:\t" << beta_ << endl;
}

void ToyZ6::PrintMinimaInOrder() {
    for (int i = 0; i < vev_order_by_angle.size(); i++) {
        cout << " p" << i << ": (" << vev_order_by_angle[i][0] << "," << vev_order_by_angle[i][1] << ")" << endl;
    }
}
}  // namespace ToyTest
