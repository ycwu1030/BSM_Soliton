#include "ToyZ3.h"

#include <cmath>
#include <iostream>

using namespace std;

ToyZ3::ToyZ3() : Basic_Model(2) {}

void ToyZ3::Set_Potential_Parameters(double r_) {
    r = r_;
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

double ToyZ3::Vtotal(VD field_values, double scale) {
    double p1 = field_values[0];
    double p2 = field_values[1];
    double pt = p1 * p1 + p2 * p2;
    double p3 = p1 * p1 * p1 - 3 * p1 * p2 * p2;
    double vtot = -pt / 2 + pt * pt / 4 - 2 * r / 3 * p3;
    return vtot / pow(scale, 4);
}

VD ToyZ3::dVtotal(VD field_values, double scale) {
    VD res(2);
    double p1 = field_values[0];
    double p2 = field_values[1];
    double pt = p1 * p1 + p2 * p2;
    res[0] = -2 * r * (p1 * p1 - p2 * p2) + p1 * pt - p1;
    res[1] = p2 * (pt + 4 * r * p1 - 1);
    return res / pow(scale, 3);
}

VVD ToyZ3::d2Vtotal(VD field_values, double scale) {
    VVD res(2, VD(2));
    double p1 = field_values[0];
    double p2 = field_values[1];
    res[0][0] = -4 * r * p1 + 3 * p1 * p1 + p2 * p2 - 1;
    res[0][1] = 2 * p2 * (p1 + 2 * r);
    res[1][0] = res[0][1];
    res[1][1] = 4 * r * p1 + p1 * p1 + 3 * p2 * p2 - 1;

    return res / pow(scale, 2);
}

double ToyZ3::V0_global(double scale) { return Vtotal({_v1global, _v2global}, scale); }

VD ToyZ3::QuarticCoupling(VD field_values) { return {1.0 / 4.0, 1.0 / 4.0}; }

void ToyZ3::FindLocalMinima() {
    if (_Solved) {
        return;
    }

    Clear_Local_Cache();

    _localExtreme.push_back({0, 0});
    AppendLocalExtreme();

    double r1 = sqrt(r * r + 1) + r;
    double r2 = -(sqrt(r * r + 1) - r);
    for (int i = 0; i < 3; i++) {
        double th = i * 2.0 * M_PI / 3.0;

        _localExtreme.push_back({r1 * cos(th), r1 * sin(th)});
        AppendLocalExtreme();

        _localExtreme.push_back({r2 * cos(th), r2 * sin(th)});
        AppendLocalExtreme();
    }

    _Solved = true;
}

void ToyZ3::PrintParameters() {
    cout << "Potential Parameter:" << endl;
    cout << "r:\t" << r << endl;
}
