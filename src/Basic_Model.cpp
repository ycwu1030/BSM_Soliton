#include "Basic_Model.h"

#include <cmath>
#include <iostream>

using namespace std;

namespace BSM_Soliton {
SM::SM() {
    alpha_EW = 1.0 / alpha1;
    ee = sqrt(4 * M_PI * alpha_EW);
    vev = pow(M_SQRT2 * GF, -0.5);
    double tmp = sqrt(M_PI * alpha_EW) * vev;
    theta_w = asin(2 * tmp / MZ) / 2.0;
    MW = tmp / sin(theta_w);
    MW2 = MW * MW;
    g_weak = ee / sin(theta_w);
    g_hyper = ee / cos(theta_w);
    yt = M_SQRT2 * MT / vev;
}
}  // namespace BSM_Soliton

SM::SM() {
    alpha = 1.0 / alpha1;
    ee = sqrt(4 * Pi * alpha);
    vev = pow(sqrt(2) * GF, -0.5);
    double A = sqrt(Pi * alpha) * vev;
    thetaW = asin(2 * A / MZ) / 2.0;
    MW = A / sin(thetaW);
    MW2 = MW * MW;
    g_weak = ee / sin(thetaW);
    gp_hyper = ee / cos(thetaW);
    yt = sqrt(2) * MT / vev;
#ifdef DEBUG
    cout << "A:  " << A << endl;
    cout << "vev: " << vev << endl;
    cout << "thetaW: " << thetaW << endl;
    cout << "MW: " << MW << endl;
    cout << "yt: " << yt << endl;
    cout << "g_weak: " << g_weak << endl;
    cout << "gp_hyper: " << gp_hyper << endl;
#endif
}

Basic_Model::Basic_Model() : Potential() {
    _Solved = false;
    _N_VEVs = 0;
    Clear_Local_Cache();
}
Basic_Model::Basic_Model(int Field_Dim) : Potential(Field_Dim) {
    _Solved = false;
    _N_VEVs = Field_Dim;
    Clear_Local_Cache();
}
void Basic_Model::Clear_Local_Cache() {
    _NLocalMinima = 0;
    _NLocalExtreme = 0;
    _IndexInput = -1;
    _localExtreme.clear();
    _MinimaIndex.clear();
    _localExtreme.clear();
    _LocalMinimaQ.clear();
    _Vtotal.clear();
}
bool Basic_Model::CheckGlobalMinimum() {
    FindLocalMinima();
    if (_IndexInput >= 0) {
        if (!_LocalMinimaQ[_IndexInput]) {
            // ! The input vacuum is not local minimum
            return false;
        }
        double VtotEW = _Vtotal[_IndexInput];
        bool goodEW = true;
        for (size_t i = 0; i < _NLocalMinima; i++) {
            if (i == _IndexInput) continue;
            goodEW *= (VtotEW <= (_Vtotal[_MinimaIndex[i]]) + 1e-5 * abs(_Vtotal[_MinimaIndex[i]]));
        }
        return goodEW;
    } else {
        cout << "Error in finding the input vacuum" << endl;
        return false;
    }
}
void Basic_Model::AppendLocalExtreme() {
    _Vtotal.push_back(Vtotal(_localExtreme[_NLocalExtreme]));
    _LocalMinimaQ.push_back(CheckHessian(_localExtreme[_NLocalExtreme]));
    if (_LocalMinimaQ[_NLocalExtreme]) {
        _MinimaIndex.push_back(_NLocalExtreme);
        _NLocalMinima++;
    }
    ++_NLocalExtreme;
}
void Basic_Model::PrintParameters() {
    cout << "No Parameter for Basic Model, please implement your own specific model from the basic one" << endl;
}
void Basic_Model::PrintLocalMinima() {
    FindLocalMinima();
    cout << "The Local Extreme points are: " << endl;
    cout << "ID\t";
    for (size_t i = 0; i < _N_VEVs; i++) {
        cout << "v_" << i << "\t";
    }
    cout << "LocalMinimumQ\tVtot" << endl;
    for (size_t i = 0; i < _NLocalExtreme; i++) {
        cout << i << "\t";
        cout << _localExtreme[i] << "\t";
        // for (size_t j = 0; j < _N_VEVs; j++)
        // {
        //     cout<<_localExtreme[i][j]<<"\t";
        // }
        cout << _LocalMinimaQ[i] << "\t" << _Vtotal[i] << endl;
    }
}
