#include "Potential.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>

#include "Eigen/Dense"

using namespace std;
using namespace Eigen;

namespace BSM_Soliton {
Potential::Potential(int FieldSpaceDimension)
    : Field_Space_Dimension(FieldSpaceDimension), Solved(false), Input_Minimum(FieldSpaceDimension, 0) {}

void Potential::Clear_Local_Cache() {
    Input_Minimum_ID = -1;
    N_Local_Extrema = 0;
    N_Local_Minima = 0;
    Minima_Index.clear();
    Local_Extrema.clear();
    is_Local_Minima.clear();
    Potential_at_Extrema.clear();
}

bool Potential::Check_Hessian_Matrix(VD field_values) {
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

void Potential::Find_Local_Extrema() {
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

void Potential::Add_Local_Extremum(VD local_extremum) {
    Local_Extrema.push_back(local_extremum);
    Potential_at_Extrema.push_back(V(local_extremum));

    bool is_minimum = Check_Hessian_Matrix(local_extremum);
    is_Local_Minima.push_back(is_minimum);
    if (is_minimum) {
        Minima_Index.push_back(Local_Extrema.size() - 1);
    }
}

void Potential::Calculate_Local_Extrema() { Add_Local_Extremum({0, 0}); }

bool Potential::Check_Global_Minimum() {
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

void Potential::Print_Local_Extrema() const {
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

VD Potential::Get_Local_Minimum(unsigned id) const {
    if (id >= N_Local_Minima) id = 0;
    return Local_Extrema[Minima_Index[id]];
}
}  // namespace BSM_Soliton

double Potential::GetTotalEnergy(VD x, VVD fields, VD dfields) {
    double energy = 0;
    double V0 = V0_global();
    for (size_t i = 0; i < x.size(); i++) {
        double density = 0;
        double DeltaZ = 0;
        if (i == 0) {
            DeltaZ = (x[1] - x[0]) / 2;
        } else if (i == x.size() - 1) {
            DeltaZ = (x[i] - x[i - 1]) / 2;
        } else {
            DeltaZ = (x[i + 1] - x[i - 1]) / 2;
        }
        density += dfields[i] * dfields[i] / 2;
        density += Vtotal(fields[i]) - V0;
        energy += density * DeltaZ;
    }
    return energy;
}
double Potential::GetWallWidth(VD x, VVD fields, VD dfields, double criteria) {
    VD x_aver;
    VD accumulated_energy;
    double energy = 0;
    double V0 = V0_global();
    for (size_t i = 0; i < x.size(); i++) {
        double density = 0;
        double DeltaZ = 0;
        if (i == 0) {
            DeltaZ = (x[1] - x[0]) / 2;
        } else if (i == x.size() - 1) {
            DeltaZ = (x[i] - x[i - 1]) / 2;
        } else {
            DeltaZ = (x[i + 1] - x[i - 1]) / 2;
        }
        density += dfields[i] * dfields[i] / 2;
        density += Vtotal(fields[i]) - V0;
        energy += density * DeltaZ;
        x_aver.push_back(x[i]);
        accumulated_energy.push_back(energy);
    }
    double energy_total = accumulated_energy.back();
    double p1 = energy_total * (1.0 - criteria) / 2.0;
    double p2 = energy_total * (1.0 + criteria) / 2.0;
    auto p1_iter = std::lower_bound(accumulated_energy.begin(), accumulated_energy.end(), p1);
    auto p2_iter = std::lower_bound(accumulated_energy.begin(), accumulated_energy.end(), p2);
    size_t ip1 = std::distance(accumulated_energy.begin(), p1_iter);
    size_t ip2 = std::distance(accumulated_energy.begin(), p2_iter);
    ip1 = ip1 < 1 ? 1 : ip1;
    ip2 = ip2 < 1 ? 1 : ip2;
    double r1 = (p1 - accumulated_energy[ip1 - 1]) / (accumulated_energy[ip1] - accumulated_energy[ip1 - 1]);
    double r2 = (p2 - accumulated_energy[ip2 - 1]) / (accumulated_energy[ip2] - accumulated_energy[ip2 - 1]);
    double x1 = x_aver[ip1 - 1] + r1 * (x_aver[ip1] - x_aver[ip1 - 1]);
    double x2 = x_aver[ip2 - 1] + r2 * (x_aver[ip2] - x_aver[ip2 - 1]);
    return x2 - x1;
}
void Potential::DumpFullSolution(VD x, VVD fields, VD dfields, std::string filename) {
    double density;
    std::ofstream output(filename.c_str());
    double V0 = V0_global();
    output << "z";
    for (size_t i = 0; i < fields.front().size(); i++) {
        output << "\tphi" << i;
    }
    output << "\tdensity" << endl;
    output << std::scientific;
    output << std::showpos;
    output << std::setprecision(10);
    for (int i = 0; i < x.size(); i++) {
        double density = 0;
        double DeltaZ = 0;
        if (i == 0) {
            DeltaZ = (x[1] - x[0]) / 2;
        } else if (i == x.size() - 1) {
            DeltaZ = (x[i] - x[i - 1]) / 2;
        } else {
            DeltaZ = (x[i + 1] - x[i - 1]) / 2;
        }
        density += dfields[i] * dfields[i] / 2;
        density += Vtotal(fields[i]) - V0;
        output << x[i] << "\t" << fields[i] << "\t" << density << endl;
    }
}
void Potential::DumpEnergyDensity(VD x, VVD fields, VD dfields, std::string filename) {
    double density;
    std::ofstream output(filename.c_str());
    double V0 = V0_global();
    output << "x\tdensity" << std::endl;
    for (int i = 0; i < x.size(); i++) {
        double density = 0;
        double DeltaZ = 0;
        if (i == 0) {
            DeltaZ = (x[1] - x[0]) / 2;
        } else if (i == x.size() - 1) {
            DeltaZ = (x[i] - x[i - 1]) / 2;
        } else {
            DeltaZ = (x[i + 1] - x[i - 1]) / 2;
        }
        density += dfields[i] * dfields[i] / 2;
        density += Vtotal(fields[i]) - V0;
        output << x[i] << "\t" << density << std::endl;
    }
}
double Potential::GetTension(VD x, VVD fields, VD dfields) {
    double tension = 0;
    for (size_t i = 0; i < x.size() - 1; i++) {
        double density = 0;
        double DeltaZ = 0;
        if (i == 0) {
            DeltaZ = (x[1] - x[0]) / 2;
        } else if (i == x.size() - 1) {
            DeltaZ = (x[i] - x[i - 1]) / 2;
        } else {
            DeltaZ = (x[i + 1] - x[i - 1]) / 2;
        }
        density += dfields[i] * dfields[i] / 2;
        tension += density * DeltaZ;
    }
    return tension;
}
bool Potential::CheckHessian(VD field_values) {
    MatrixXd HM(_Field_Dim, _Field_Dim);
    VVD HM_VVD = d2Vtotal(field_values);
    for (size_t i = 0; i < _Field_Dim; i++) {
        for (size_t j = 0; j < _Field_Dim; j++) {
            HM(i, j) = HM_VVD[i][j];
        }
    }
    SelfAdjointEigenSolver<MatrixXd> eigensolver(HM);
    if (eigensolver.info() != Success) return false;
#if DEBUG
    cout << "Eigen Values: " << endl;
    cout << eigensolver.eigenvalues() << endl;
    cout << "--" << endl;
#endif
    return ((eigensolver.eigenvalues()).array() > 0).all();
}
