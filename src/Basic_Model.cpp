#include "Basic_Model.h"

#include <cmath>
#include <iostream>

#include "Eigen/Dense"

using namespace std;
using namespace Eigen;

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

BaseModel::BaseModel(int FieldSpaceDimension)
    : Field_Space_Dimension(FieldSpaceDimension), Solved(false), Input_Minimum(FieldSpaceDimension, 0) {}

void BaseModel::Clear_Local_Cache() {
    Input_Minimum_ID = -1;
    N_Local_Extrema = 0;
    N_Local_Minima = 0;
    Minima_Index.clear();
    Local_Extrema.clear();
    is_Local_Minima.clear();
    Potential_at_Extrema.clear();
}

bool BaseModel::Check_Hessian_Matrix(VD field_values) {
    MatrixXd HM(Field_Space_Dimension, Field_Space_Dimension);
    VVD HM_VVD = d2V(field_values);
    for (size_t i = 0; i < Field_Space_Dimension; i++) {
        for (size_t j = 0; j < Field_Space_Dimension; j++) {
            HM(i, j) = HM_VVD[i][j];
        }
    }
    SelfAdjointEigenSolver<MatrixXd> eigensolver(HM);
    if (eigensolver.info() != Success) return false;
    return ((eigensolver.eigenvalues()).array() > 0).all();
}

void BaseModel::Find_Local_Extrema() {
    if (Solved == true) return;
    Clear_Local_Cache();
    Calculate_Local_Extrema();
    N_Local_Extrema = Local_Extrema.size();
    N_Local_Minima = Minima_Index.size();
    for (size_t i = 0; i < N_Local_Extrema; i++) {
        if (CloseQ(Input_Minimum, Local_Extrema[i])) {
            Input_Minimum_ID = i;
        }
    }
    Solved = true;
}

void BaseModel::Add_Local_Extremum(VD local_extremum) {
    Local_Extrema.push_back(local_extremum);
    Potential_at_Extrema.push_back(V(local_extremum));

    bool is_minimum = Check_Hessian_Matrix(local_extremum);
    is_Local_Minima.push_back(is_minimum);
    if (is_minimum) {
        Minima_Index.push_back(Local_Extrema.size() - 1);
    }
}

void BaseModel::Calculate_Local_Extrema() { Add_Local_Extremum({0, 0}); }

bool BaseModel::Check_Global_Minimum() {
    Find_Local_Extrema();
    if (Input_Minimum_ID < 0) {
        // * Didn't find the input minimum in the solutions
        return false;
    }
    if (!is_Local_Minima[Input_Minimum_ID]) {
        // * The input minimum is not actually a local minimum
        return false;
    }
    // * We find the input minimum and it is indeed a local minimum
    // * Then check whether it is the global minimum
    double V_at_input_minimum = Potential_at_Extrema[Input_Minimum_ID];
    bool good_input = true;
    for (int i = 0; i < N_Local_Minima; i++) {
        int index = Minima_Index[i];
        if (index == Input_Minimum_ID) continue;
        double V_test = Potential_at_Extrema[index];
        good_input *= (V_at_input_minimum <= (V_test + 1e-5 * abs(V_test)));
    }
    return good_input;
}

void BaseModel::Print_Local_Extrema() const {
    cout << "The Local Extreme points are: " << endl;
    cout << "ID\t";
    for (int i = 0; i < Field_Space_Dimension; i++) {
        cout << "v_" << i << "\t";
    }
    cout << "is_local_minimum\tV" << endl;
    for (int i = 0; i < N_Local_Extrema; i++) {
        cout << i << "\t" << Local_Extrema[i] << "\t" << is_Local_Minima[i] << "\t" << Potential_at_Extrema[i] << endl;
    }
}

VD BaseModel::Get_Local_Minimum(unsigned id) const {
    if (id >= N_Local_Minima) id = 0;
    return Local_Extrema[Minima_Index[id]];
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
