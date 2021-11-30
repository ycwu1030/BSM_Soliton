#include "ToyZ4.h"

#include <cmath>
#include <iostream>

using namespace std;

ToyZ4::ToyZ4() : Basic_Model(2) {}

void ToyZ4::Set_Potential_Parameters(double beta_) {
    beta = beta_;
    r = 2 * beta / (1 - beta);
    _Solved = false;
    FindLocalMinima();

    _v1global = 0;
    _v2global = 0;
    _vr = 0;
    _IndexInput = -1;
    double Vtmp = 0;
    for (size_t i = 0; i < _NLocalExtreme; i++) {
        if (_LocalMinimaQ[i]) {
            double vtt = sqrt(pow(_localExtreme[i][0], 2) + pow(_localExtreme[i][1], 2));
            if (vtt > _vr) {
                _vr = vtt;
            }
        }
        if (_Vtotal[i] < Vtmp && _LocalMinimaQ[i]) {
            _v2global = _localExtreme[i][1];
            _v1global = _localExtreme[i][0];
            _IndexInput = i;
            Vtmp = _Vtotal[i];
        }
    }
}

double ToyZ4::Vtotal(VD field_values, double scale) {
    double p1 = field_values[0];
    double p2 = field_values[1];
    double pt = p1 * p1 + p2 * p2;
    double vtot = -pt / 2 + pt * pt / 4 + r * p1 * p1 * p2 * p2;
    return vtot / pow(scale, 4);
}

VD ToyZ4::dVtotal(VD field_values, double scale) {
    VD res(2);
    double p1 = field_values[0];
    double p2 = field_values[1];
    res[0] = p1 * (p1 * p1 + p2 * p2) + 2 * r * p1 * p2 * p2 - p1;
    res[1] = p2 * (p1 * p1 + p2 * p2) + 2 * r * p1 * p1 * p2 - p2;
    return res / pow(scale, 3);
}

VVD ToyZ4::d2Vtotal(VD field_values, double scale) {
    VVD res(2, VD(2));
    double p1 = field_values[0];
    double p2 = field_values[1];
    res[0][0] = 3 * p1 * p1 + (1 + 2 * r) * p2 * p2 - 1;
    res[0][1] = 2 * (1 + 2 * r) * p1 * p2;
    res[1][0] = res[0][1];
    res[1][1] = (1 + 2 * r) * p1 * p1 + 3 * p2 * p2 - 1;

    return res / pow(scale, 2);
}

double ToyZ4::V0_global(double scale) { return Vtotal({_v1global, _v2global}, scale); }

VD ToyZ4::QuarticCoupling(VD field_values) { return {1.0 / 4.0, 1.0 / 4.0}; }

void ToyZ4::FindLocalMinima() {
    if (_Solved) {
        return;
    }

    Clear_Local_Cache();

    _localExtreme.push_back({0, 0});
    AppendLocalExtreme();

    double r1 = 1.0;
    _localExtreme.push_back({r1, 0});
    AppendLocalExtreme();

    _localExtreme.push_back({-r1, 0});
    AppendLocalExtreme();

    _localExtreme.push_back({0, r1});
    AppendLocalExtreme();

    _localExtreme.push_back({0, -r1});
    AppendLocalExtreme();

    double t2_tmp = 2 + 2 * r;
    if (t2_tmp > 0) {
        double r2 = 1.0 / sqrt(2 + 2 * r);

        _localExtreme.push_back({r2, r2});
        AppendLocalExtreme();

        _localExtreme.push_back({r2, -r2});
        AppendLocalExtreme();

        _localExtreme.push_back({-r2, r2});
        AppendLocalExtreme();

        _localExtreme.push_back({-r2, -r2});
        AppendLocalExtreme();
    }

    _Solved = true;
}

void ToyZ4::PrintParameters() {
    cout << "Potential Parameter:" << endl;
    cout << "r:\t" << r << endl;
}
