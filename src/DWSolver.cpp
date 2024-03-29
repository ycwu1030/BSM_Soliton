#include "DWSolver.h"

#include <algorithm>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
using namespace std;

namespace BSM_Soliton {
DomainWallSolverOriginal::DomainWallSolverOriginal(BaseModel *model)
    : mod(model),
      Mesh_Size(401),
      Field_Space_Dim(model->Get_Field_Space_Dimension()),
      RelaxationODE(model->Get_Field_Space_Dimension() * 2, model->Get_Field_Space_Dimension()),
      Fields_at_Left(model->Get_Field_Space_Dimension(), 0),
      Fields_at_Right(model->Get_Field_Space_Dimension(), 0) {
    solver = new Relaxation(this, 1.0);
}

void DomainWallSolverOriginal::dYdX(MeshPoint &point) {
    int n = Field_Space_Dim;
    VD &y = point.Y;
    double t = point.X;
    VD fields(y.begin(), y.begin() + n);
    VD dV = mod->dV(fields);
    VVD d2V = mod->d2V(fields);
    for (size_t i = 0; i < n; i++) {
        point.Result[i] = y[i + n] / cos(t) / cos(t);
        point.Result[i + n] = dV[i] / cos(t) / cos(t);
    }

    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            double delta_ij = i == j ? 1 : 0;
            point.dResultdY[i][j] = 0;
            point.dResultdY[i][j + n] = delta_ij / cos(t) / cos(t);
            point.dResultdY[i + n][j] = d2V[i][j] / cos(t) / cos(t);
            point.dResultdY[i + n][j + n] = 0;
        }
    }
}

void DomainWallSolverOriginal::Left_Boundary_Constraints(MeshPoint &point) {
    VD &y = point.Y;
    int n = Field_Space_Dim;
    for (size_t i = 0; i < n; i++) {
        point.Result[i] = y[i] - Fields_at_Left[i];
        for (size_t j = 0; j < n; j++) {
            double delta_ij = i == j ? 1 : 0;
            point.dResultdY[i][j] = delta_ij;
            point.dResultdY[i][j + n] = 0;
        }
    }
}

void DomainWallSolverOriginal::Right_Boundary_Constraints(MeshPoint &point) {
    VD &y = point.Y;
    int n = Field_Space_Dim;
    for (size_t i = 0; i < n; i++) {
        point.Result[i] = y[i] - Fields_at_Right[i];
        for (size_t j = 0; j < n; j++) {
            double delta_ij = i == j ? 1 : 0;
            point.dResultdY[i][j] = delta_ij;
            point.dResultdY[i][j + n] = 0;
        }
    }
}

bool DomainWallSolverOriginal::Solve(const VD &field_at_left, const VD &field_at_right) {
    Fields_at_Left = field_at_left;
    Fields_at_Right = field_at_right;
    double z_max_guess = Guess_Z_Range();
    VD qs = linspace(-M_PI_2, M_PI_2, Mesh_Size);
    VD left_boundary(DOF, 0);
    VD right_boundary(DOF, 0);
    int n = Field_Space_Dim;
    for (size_t i = 0; i < n; i++) {
        left_boundary[i] = Fields_at_Left[i];
        left_boundary[i + n] = 0;
        right_boundary[i] = Fields_at_Right[i];
        right_boundary[i + n] = 0;
    }

    VVD grids = linspace(left_boundary, right_boundary, Mesh_Size);
    bool solved = solver->Solve(qs, grids);
    result = solver->Get_Solution();
    return solved;
}

double DomainWallSolverOriginal::Guess_Z_Range() {
    double zmax = 10.0;
    VD couplings = mod->Quartic_Couplings(Fields_at_Left);
    for (size_t i = 0; i < Field_Space_Dim; i++) {
        if (CloseQ(Fields_at_Left[i], Fields_at_Right[i])) continue;
        double zi = 6.0 / sqrt(couplings[i] / 2);
        if (zi > zmax) {
            zmax = zi;
        }
    }
    return zmax;
}

double DomainWallSolverOriginal::Get_Tension() {
    double tension = 0;
    double V0 = mod->V_min();
    int n = Field_Space_Dim;
    for (size_t i = 0; i < result.size() - 1; i++) {
        VD &y_l = result[i].Y;
        double z_l = result[i].X;
        VD phi_l(y_l.begin(), y_l.begin() + n);
        VD dphidz_l(y_l.begin() + n, y_l.begin() + n + n);
        double density_l = dphidz_l * dphidz_l / 2.0;
        density_l += mod->V(phi_l) - V0;

        VD &y_h = result[i + 1].Y;
        double z_h = result[i + 1].X;
        VD phi_h(y_h.begin(), y_h.begin() + n);
        VD dphidz_h(y_h.begin() + n, y_h.begin() + n + n);
        double density_h = dphidz_h * dphidz_h / 2.0;
        density_h += mod->V(phi_h) - V0;

        tension += (density_l + density_h) * (z_h - z_l) / 2.0;
    }
    return tension;
}

double DomainWallSolverOriginal::Get_Wall_Width(double criteria) {
    VD z_aver;
    VD accumulated_tension;
    double tension = 0;
    double V0 = mod->V_min();
    int n = Field_Space_Dim;
    for (size_t i = 0; i < result.size() - 1; i++) {
        VD &y_l = result[i].Y;
        double z_l = result[i].X;
        VD phi_l(y_l.begin(), y_l.begin() + n);
        VD dphidz_l(y_l.begin() + n, y_l.begin() + n + n);
        double density_l = dphidz_l * dphidz_l / 2.0;
        density_l += mod->V(phi_l) - V0;

        VD &y_h = result[i + 1].Y;
        double z_h = result[i + 1].X;
        VD phi_h(y_h.begin(), y_h.begin() + n);
        VD dphidz_h(y_h.begin() + n, y_h.begin() + n + n);
        double density_h = dphidz_h * dphidz_h / 2.0;
        density_h += mod->V(phi_h) - V0;

        tension += (density_l + density_h) * (z_h - z_l) / 2.0;

        z_aver.push_back((z_l + z_h) / 2.0);
        accumulated_tension.push_back(tension);
    }
    double p1 = tension * (1.0 - criteria) / 2.0;
    double p2 = tension * (1.0 + criteria) / 2.0;
    auto p1_iter = std::lower_bound(accumulated_tension.begin(), accumulated_tension.end(), p1);
    auto p2_iter = std::lower_bound(accumulated_tension.begin(), accumulated_tension.end(), p2);
    size_t ip1 = std::distance(accumulated_tension.begin(), p1_iter);
    size_t ip2 = std::distance(accumulated_tension.begin(), p2_iter);
    ip1 = ip1 < 1 ? 1 : ip1;
    ip2 = ip2 < 1 ? 1 : ip2;
    double r1 = (p1 - accumulated_tension[ip1 - 1]) / (accumulated_tension[ip1] - accumulated_tension[ip1 - 1]);
    double r2 = (p2 - accumulated_tension[ip2 - 1]) / (accumulated_tension[ip2] - accumulated_tension[ip2 - 1]);
    double x1 = z_aver[ip1 - 1] + r1 * (z_aver[ip1] - z_aver[ip1 - 1]);
    double x2 = z_aver[ip2 - 1] + r2 * (z_aver[ip2] - z_aver[ip2 - 1]);
    return x2 - x1;
}

void DomainWallSolverOriginal::Dump_Solution(string filename) {
    ofstream output(filename.c_str());
    int n = Field_Space_Dim;
    output << "z";
    for (size_t i = 0; i < n; i++) {
        output << "\tphi" << i;
    }
    for (size_t i = 0; i < n; i++) {
        output << "\tdphi" << i << "dz";
    }
    output << "\tepsilon_z" << endl;
    double V0 = mod->V_min();
    output << std::scientific;
    output << std::showpos;
    output << std::setprecision(8);
    for (size_t i = 0; i < result.size(); i++) {
        VD &y_l = result[i].Y;
        double z_l = result[i].X;
        VD phi_l(y_l.begin(), y_l.begin() + n);
        VD dphidz_l(y_l.begin() + n, y_l.begin() + n + n);
        double density_l = dphidz_l * dphidz_l / 2.0;
        density_l += mod->V(phi_l) - V0;

        output << z_l;
        for (size_t j = 0; j < n; j++) {
            output << "\t" << y_l[j];
        }
        for (size_t j = 0; j < n; j++) {
            output << "\t" << y_l[j + n];
        }
        output << "\t" << density_l << endl;
    }
    output.close();
}

DomainWallSolver::DomainWallSolver(BaseModel *model)
    : mod(model),
      Mesh_Size(401),
      Field_Space_Dim(model->Get_Field_Space_Dimension()),
      F_ratio(1.0),
      RelaxationODE(model->Get_Field_Space_Dimension() * 2 + 3, model->Get_Field_Space_Dimension() + 1),
      dF_value(model->Get_Field_Space_Dimension() * 2 + 3, 0),
      Fields_at_Left(model->Get_Field_Space_Dimension(), 0),
      Fields_at_Right(model->Get_Field_Space_Dimension(), 0) {
    solver = new Relaxation(this, 1.0);
}

// void DomainWallSolver::Calc_F_Variables(MeshPoint &point) {
//     VD &y = point.Y;
//     double tmp = 0;
//     double c1 = 1.0 / (1.0 + F_ratio);
//     double c2 = F_ratio / (1.0 + F_ratio);
//     for (size_t i = 0; i < Field_Space_Dim; i++) {
//         tmp += pow(y[i + Field_Space_Dim + 1], 2) / (pow(y[i + 1], 2) + 1e-3);
//         dF_value[i + 1] = -2 * c2 * pow(y[i + Field_Space_Dim + 1], 2) / pow(pow(y[i + 1], 2) + 1e-3, 2) * (y[i +
//         1]); dF_value[i + Field_Space_Dim + 1] = 2 * c2 * y[i + Field_Space_Dim + 1] / (pow(y[i + 1], 2) + 1e-3);
//     }
//     F_value = c1 * Field_Space_Dim / 4.0 + c2 * tmp;
//     dF_value[0] = 0;
//     dF_value[2 * Field_Space_Dim + 1] = 0;
//     dF_value[2 * Field_Space_Dim + 2] = 0;
// }

void DomainWallSolver::Calc_F_Variables(MeshPoint &point) {
    VD &y = point.Y;
    int n = Field_Space_Dim;
    double tmp = 0;
    double c1 = 1.0 / (1.0 + F_ratio);
    double c2 = F_ratio / (1.0 + F_ratio);
    double epsilon = 1e-3;
    VD field(y.begin() + 1, y.begin() + 1 + n);
    VD dV = mod->dV(field);
    VVD d2V = mod->d2V(field);
    for (size_t i = 0; i < n; i++) {
        tmp += pow(y[i + 1 + n], 2) / (pow(y[i + 1], 2) + epsilon);
        tmp += pow(dV[i], 2) / (pow(y[i + 1 + n], 2) + epsilon);
    }
    F_value = c1 * Field_Space_Dim / 2.0 + c2 * tmp;

    dF_value[0] = 0;
    dF_value[2 * n + 1] = 0;
    dF_value[2 * n + 2] = 0;
    for (size_t i = 0; i < n; i++) {
        dF_value[i + 1] = -2 * c2 * pow(y[i + 1 + n], 2) / pow(pow(y[i + 1], 2) + epsilon, 2) * y[i + 1];
        for (size_t j = 0; j < n; j++) {
            dF_value[i + 1] += 2 * c2 * dV[j] * d2V[i][j] / (pow(y[j + 1 + n], 2) + epsilon);
        }
        dF_value[i + 1 + n] = 2 * c2 * y[i + 1 + n] / (pow(y[i + 1], 2) + epsilon) -
                              2 * c2 * y[i + 1 + n] * pow(dV[i], 2) / pow(pow(y[i + 1 + n], 2) + epsilon, 2);
    }
}

void DomainWallSolver::dYdX(MeshPoint &point) {
    Calc_F_Variables(point);
    int n = Field_Space_Dim;
    VD &y = point.Y;
    VD field(y.begin() + 1, y.begin() + 1 + n);
    VD dV = mod->dV(field);
    VVD d2V = mod->d2V(field);
    double tmp = y[2 * n + 2] / F_value;
    point.Result[0] = tmp;
    for (size_t i = 0; i < n; i++) {
        point.Result[i + 1] = y[i + 1 + n] * y[2 * n + 1] * tmp;
        point.Result[i + 1 + n] = dV[i] * y[2 * n + 1] * tmp;
    }
    point.Result[2 * n + 1] = 0;
    point.Result[2 * n + 2] = 0;

    point.dResultdY[0][0] = 0;
    for (int i = 0; i < n; i++) {
        point.dResultdY[0][i + 1] = -y[2 * n + 2] / F_value / F_value * dF_value[i + 1];
        point.dResultdY[0][i + 1 + n] = -y[2 * n + 2] / F_value / F_value * dF_value[i + 1 + n];
    }
    point.dResultdY[0][2 * n + 1] = 0;
    point.dResultdY[0][2 * n + 2] = 1.0 / F_value;

    for (int i = 0; i < n; i++) {
        point.dResultdY[i + 1][0] = 0;
        point.dResultdY[i + 1 + n][0] = 0;
        for (int j = 0; j < n; j++) {
            double delta_ij = i == j ? 1.0 : 0.0;
            double tmp = (y[2 * n + 1]) * y[2 * n + 2] / F_value;
            point.dResultdY[i + 1][j + 1] = -y[i + 1 + n] * tmp / F_value * dF_value[j + 1];
            point.dResultdY[i + 1][j + 1 + n] = tmp * (delta_ij - y[i + 1 + n] / F_value * dF_value[j + 1 + n]);
            point.dResultdY[i + 1 + n][j + 1] = tmp * (d2V[i][j] - dV[i] / F_value * dF_value[j + 1]);
            point.dResultdY[i + 1 + n][j + 1 + n] = -dV[i] * tmp / F_value * dF_value[j + 1 + n];
        }
        point.dResultdY[i + 1][2 * n + 1] = y[i + 1 + n] * y[2 * n + 2] / F_value;
        point.dResultdY[i + 1][2 * n + 2] = y[i + 1 + n] * y[2 * n + 1] / F_value;

        point.dResultdY[i + 1 + n][2 * n + 1] = dV[i] * y[2 * n + 2] / F_value;
        point.dResultdY[i + 1 + n][2 * n + 2] = dV[i] * y[2 * n + 1] / F_value;
    }

    for (int i = 0; i < 2 * n + 3; i++) {
        point.dResultdY[2 * n + 1][i] = 0;
        point.dResultdY[2 * n + 2][i] = 0;
    }
}

void DomainWallSolver::Left_Boundary_Constraints(MeshPoint &point) {
    int n = Field_Space_Dim;
    VD &y = point.Y;
    point.Result[0] = y[0];
    for (int i = 0; i < n; i++) {
        point.Result[i + 1] = y[i + 1] - Fields_at_Left[i];
    }

    point.dResultdY[0][0] = 1;
    for (int i = 0; i < n; i++) {
        point.dResultdY[0][i + 1] = 0;
        point.dResultdY[0][i + 1 + n] = 0;
    }
    point.dResultdY[0][2 * n + 1] = 0;
    point.dResultdY[0][2 * n + 2] = 0;

    for (int i = 0; i < n; i++) {
        point.dResultdY[i + 1][0] = 0;
        for (int j = 0; j < n; j++) {
            double delta_ij = i == j ? 1 : 0;
            point.dResultdY[i + 1][j + 1] = delta_ij;
            point.dResultdY[i + 1][j + 1 + n] = 0;
        }
        point.dResultdY[i + 1][2 * n + 1] = 0;
        point.dResultdY[i + 1][2 * n + 2] = 0;
    }
}

void DomainWallSolver::Right_Boundary_Constraints(MeshPoint &point) {
    int n = Field_Space_Dim;
    VD &y = point.Y;
    point.Result[0] = y[0] - 1;
    double tmp = 0;
    for (int i = 0; i < n; i++) {
        point.Result[i + 1] = y[i + 1] - Fields_at_Right[i];
        tmp += abs(y[i + 1 + n]);
    }
    point.Result[n + 1] = tmp;

    point.dResultdY[0][0] = 1;
    for (int i = 0; i < n; i++) {
        point.dResultdY[0][i + 1] = 0;
        point.dResultdY[0][i + 1 + n] = 0;
    }
    point.dResultdY[0][2 * n + 1] = 0;
    point.dResultdY[0][2 * n + 2] = 0;

    for (int i = 0; i < n; i++) {
        point.dResultdY[i + 1][0] = 0;
        for (int j = 0; j < n; j++) {
            double delta_ij = i == j ? 1 : 0;
            point.dResultdY[i + 1][j + 1] = delta_ij;
            point.dResultdY[i + 1][j + 1 + n] = 0;
        }
        point.dResultdY[i + 1][2 * n + 1] = 0;
        point.dResultdY[i + 1][2 * n + 2] = 0;
    }

    point.dResultdY[n + 1][0] = 0;
    for (int i = 0; i < n; i++) {
        point.dResultdY[n + 1][i + 1] = 0;
        point.dResultdY[n + 1][i + 1 + n] = y[i + 1 + n] >= 0 ? 1 : -1;
    }
    point.dResultdY[n + 1][2 * n + 1] = 0;
    point.dResultdY[n + 1][2 * n + 2] = 0;
}

bool DomainWallSolver::Solve(const VD &field_at_left, const VD &field_at_right) {
    Fields_at_Left = field_at_left;
    Fields_at_Right = field_at_right;
    VD qs = linspace(0.0, Mesh_Size - 1.0, Mesh_Size);
    VD left_boundary(DOF, 0);
    VD right_boundary(DOF, 0);
    double z_max_guess = Guess_Z_Range();
    int n = Field_Space_Dim;
    left_boundary[0] = 0;
    right_boundary[0] = 1;
    for (size_t i = 0; i < n; i++) {
        left_boundary[i + 1] = Fields_at_Left[i];
        left_boundary[i + 1 + n] = 0;
        right_boundary[i + 1] = Fields_at_Right[i];
        right_boundary[i + 1 + n] = 0;
    }
    left_boundary[2 * n + 1] = 2 * z_max_guess;
    left_boundary[2 * n + 2] = 1;
    right_boundary[2 * n + 1] = 2 * z_max_guess;
    right_boundary[2 * n + 2] = 1;

    VVD grids = linspace(left_boundary, right_boundary, Mesh_Size);
    bool solved = solver->Solve(qs, grids);
    result = solver->Get_Solution();
    return solved;
}

double DomainWallSolver::Guess_Z_Range() {
    double zmax = 10.0;
    VD couplings = mod->Quartic_Couplings(Fields_at_Left);
    for (size_t i = 0; i < Field_Space_Dim; i++) {
        if (CloseQ(Fields_at_Left[i], Fields_at_Right[i])) continue;
        double zi = 6.0 / sqrt(couplings[i] / 2);
        if (zi > zmax) {
            zmax = zi;
        }
    }
    return zmax;
}

double DomainWallSolver::Get_Tension() {
    double tension = 0;
    double V0 = mod->V_min();
    int n = Field_Space_Dim;
    double z_min = 0;
    double z_max = result[0].Y[2 * n + 1];
    for (size_t i = 0; i < result.size() - 1; i++) {
        VD &y_l = result[i].Y;
        double z_l = z_min + (z_max - z_min) * y_l[0];
        VD phi_l(y_l.begin() + 1, y_l.begin() + 1 + n);
        VD dphidz_l(y_l.begin() + 1 + n, y_l.begin() + 1 + n + n);
        double density_l = dphidz_l * dphidz_l / 2.0;
        density_l += mod->V(phi_l) - V0;

        VD &y_h = result[i + 1].Y;
        double z_h = z_min + (z_max - z_min) * y_h[0];
        VD phi_h(y_h.begin() + 1, y_h.begin() + 1 + n);
        VD dphidz_h(y_h.begin() + 1 + n, y_h.begin() + 1 + n + n);
        double density_h = dphidz_h * dphidz_h / 2.0;
        density_h += mod->V(phi_h) - V0;

        tension += (density_l + density_h) * (z_h - z_l) / 2.0;
    }
    return tension;
}

double DomainWallSolver::Get_Wall_Width(double criteria) {
    VD z_aver;
    VD accumulated_tension;
    double tension = 0;
    double V0 = mod->V_min();
    int n = Field_Space_Dim;
    double z_min = 0;
    double z_max = result[0].Y[2 * n + 1];
    for (size_t i = 0; i < result.size() - 1; i++) {
        VD &y_l = result[i].Y;
        double z_l = z_min + (z_max - z_min) * y_l[0];
        VD phi_l(y_l.begin() + 1, y_l.begin() + 1 + n);
        VD dphidz_l(y_l.begin() + 1 + n, y_l.begin() + 1 + n + n);
        double density_l = dphidz_l * dphidz_l / 2.0;
        density_l += mod->V(phi_l) - V0;

        VD &y_h = result[i + 1].Y;
        double z_h = z_min + (z_max - z_min) * y_h[0];
        VD phi_h(y_h.begin() + 1, y_h.begin() + 1 + n);
        VD dphidz_h(y_h.begin() + 1 + n, y_h.begin() + 1 + n + n);
        double density_h = dphidz_h * dphidz_h / 2.0;
        density_h += mod->V(phi_h) - V0;

        tension += (density_l + density_h) * (z_h - z_l) / 2.0;

        z_aver.push_back((z_l + z_h) / 2.0);
        accumulated_tension.push_back(tension);
    }
    double p1 = tension * (1.0 - criteria) / 2.0;
    double p2 = tension * (1.0 + criteria) / 2.0;
    auto p1_iter = std::lower_bound(accumulated_tension.begin(), accumulated_tension.end(), p1);
    auto p2_iter = std::lower_bound(accumulated_tension.begin(), accumulated_tension.end(), p2);
    size_t ip1 = std::distance(accumulated_tension.begin(), p1_iter);
    size_t ip2 = std::distance(accumulated_tension.begin(), p2_iter);
    ip1 = ip1 < 1 ? 1 : ip1;
    ip2 = ip2 < 1 ? 1 : ip2;
    double r1 = (p1 - accumulated_tension[ip1 - 1]) / (accumulated_tension[ip1] - accumulated_tension[ip1 - 1]);
    double r2 = (p2 - accumulated_tension[ip2 - 1]) / (accumulated_tension[ip2] - accumulated_tension[ip2 - 1]);
    double x1 = z_aver[ip1 - 1] + r1 * (z_aver[ip1] - z_aver[ip1 - 1]);
    double x2 = z_aver[ip2 - 1] + r2 * (z_aver[ip2] - z_aver[ip2 - 1]);
    return x2 - x1;
}

void DomainWallSolver::Dump_Solution(string filename) {
    ofstream output(filename.c_str());
    int n = Field_Space_Dim;
    output << "q\tt\tz";
    for (size_t i = 0; i < n; i++) {
        output << "\tphi" << i;
    }
    for (size_t i = 0; i < n; i++) {
        output << "\tdphi" << i << "dz";
    }
    output << "\tzmax\tepsilon_z" << endl;
    double V0 = mod->V_min();
    double z_max = result[0].Y[2 * n + 1];
    output << std::scientific;
    output << std::showpos;
    output << std::setprecision(8);
    for (size_t i = 0; i < result.size(); i++) {
        VD &y_l = result[i].Y;
        double z_l = (z_max)*y_l[0];
        VD phi_l(y_l.begin() + 1, y_l.begin() + 1 + n);
        VD dphidz_l(y_l.begin() + 1 + n, y_l.begin() + 1 + n + n);
        double density_l = dphidz_l * dphidz_l / 2.0;
        density_l += mod->V(phi_l) - V0;

        output << result[i].X << "\t" << y_l[0] << "\t" << z_l;
        for (size_t j = 0; j < n; j++) {
            output << "\t" << y_l[j + 1];
        }
        for (size_t j = 0; j < n; j++) {
            output << "\t" << y_l[j + 1 + n];
        }
        output << "\t" << y_l[2 * n + 1];
        output << "\t" << density_l << endl;
    }
    output.close();
}
}  // namespace BSM_Soliton

void DIFEQ_DW(const Relaxation_Param relax_param, void *param, VVD &S) {
    DWSolver *solver = (DWSolver *)param;
    solver->SetDWODE(relax_param, S);
}
DWSolver::DWSolver() {
    _mod = nullptr;
    _z_range = -1;
    SetMeshPoints();
    // SetZRange();
}
DWSolver::DWSolver(Potential *mod, int mesh_points) {
    _mod = mod;
    _N_Fields = _mod->GetFieldDimension();
    _ODE_DOF = 2 * _N_Fields;
    _N_Left_Bound = _N_Fields;
    _N_Right_Bound = _N_Fields;
    _mesh_points = mesh_points;
    _z_range = -1;
    // SetZRange();
}
DWSolver::DWSolver(Potential *mod, VD Left_Bound, VD Right_Bound, int mesh_points) {
    _mod = mod;
    _N_Fields = _mod->GetFieldDimension();
    _ODE_DOF = 2 * _N_Fields;
    _N_Left_Bound = _N_Fields;
    _N_Right_Bound = _N_Fields;

    _z_range = -1;
    SetBoundary(Left_Bound, Right_Bound);
    SetMeshPoints(mesh_points);
}
void DWSolver::SetZRange(double z_range) {
    _z_range = z_range;
    _x_half_range = atan(z_range);
}
void DWSolver::SetZRange() {
    double zmax = 0;
    VD couplings = _mod->QuarticCoupling(_Left_Bound);
    for (int i = 0; i < _N_Fields; i++) {
        if (_Left_Bound[i] * _Right_Bound[i] < 0) {
            double zi = 6.0 / sqrt(couplings[i] / 2);
            if (zi > zmax) {
                zmax = zi;
            }
        }
    }
    _z_range = zmax;
    _x_half_range = atan(_z_range);
}
void DWSolver::SetBoundary(VD Left_Bound, VD Right_Bound) {
    _Left_Bound = Left_Bound;
    _Right_Bound = Right_Bound;
    SetZRange();
    VD scales;
    for (int i = 0; i < _N_Fields; i++) {
        if (_Left_Bound[i] * _Right_Bound[i] < 0) {
            scales.push_back(max(abs(_Left_Bound[i]), abs(_Right_Bound[i])));
        } else {
            double center = (_Left_Bound[i] + _Right_Bound[i]) / 2;
            double left_left = _Left_Bound[i] - center;
            double right_left = _Right_Bound[i] - center;
            double scale = max(abs(left_left), abs(right_left));
            scale = scale < 1 ? 1 : scale;
            scales.push_back(scale);
        }
    }
    SetOverallScale(scales);
}
void DWSolver::SetOverallScale(VD overall_scale) {
    _overall_scale = abs(overall_scale);
    _z_scale = *max_element(_overall_scale.begin(), _overall_scale.end());
}
void DWSolver::SetOverallScale(double overall_scale) {
    _overall_scale = VD(_mod->GetFieldDimension(), overall_scale);
    _z_scale = overall_scale;
}
void DWSolver::SetInitial() {
    VD X_init;
    VVD Y_init(_ODE_DOF);
    double z_max = tan(_x_half_range);
    double dz = 2 * z_max / (_mesh_points - 1);
    for (int i = 0; i < _mesh_points; i++) {
        X_init.push_back(atan(-z_max + i * dz));
    }
    X_init[0] = -_x_half_range;
    X_init.back() = _x_half_range;
    // cout<<X_init<<endl;
    _Field_Basis.clear();
    for (int i = 0; i < _N_Fields; i++) {
        VD Yi;
        VD dYi;
        function<double(double)> f, df;
        double f_aver = (_Left_Bound[i] + _Right_Bound[i]) / 2.0;
        _Field_Basis.push_back(f_aver);
        double left_left = _Left_Bound[i] - f_aver;
        double right_left = _Right_Bound[i] - f_aver;
        if (_Left_Bound[i] * _Right_Bound[i] > 0) {
            f = [&](double x) {
                double r = tan(_x_half_range) / 5.0;
                return left_left / _overall_scale[i] +
                       (right_left - left_left) / _overall_scale[i] / (2 * _x_half_range) * (x + _x_half_range) +
                       1 / pow(cosh(tan(x) / r), 2);
            };
            df = [&](double x) {
                double r = tan(_x_half_range) / 5.0;
                return (right_left - left_left) / _overall_scale[i] / (2 * _x_half_range) -
                       2 * tanh(tan(x) / r) / r / pow(cosh(tan(x) / r), 2);
            };
        } else {
            f = [&](double x) {
                double r = tan(_x_half_range) / 5.0;
                return -left_left / _overall_scale[i] * tanh(tan(x) / r);
            };
            df = [&](double x) {
                double r = tan(_x_half_range) / 5.0;
                return -left_left / _overall_scale[i] / pow(cosh(tan(x) / r), 2);
            };
        }
        // f = [&](double x){double r = _x_half_range/5.0; return left_left/_overall_scale[i] + (right_left -
        // left_left)/_overall_scale[i]/(2*_x_half_range)*(x+_x_half_range);}; df = [&](double x){double r =
        // _x_half_range/5.0; return (right_left - left_left)/_overall_scale[i]/(2*_x_half_range);};
        for (int j = 0; j < _mesh_points; j++) {
            Yi.push_back(f(X_init[j]));
            dYi.push_back(df(X_init[j]));
        }
        int index_end = _mesh_points * 0.05;
        double Yiend = Yi[index_end];
        double Yibegin = left_left / _overall_scale[i];
        for (int j = 0; j <= index_end; j++) {
            Yi[j] = Yibegin + j * (Yiend - Yibegin) / index_end;
            dYi[j] = (Yiend - Yibegin) / (X_init[index_end] - X_init[0]);
        }
        int index_begin = _mesh_points - index_end;
        Yibegin = Yi[index_begin];
        Yiend = right_left / _overall_scale[i];
        for (int j = index_begin; j < _mesh_points; j++) {
            Yi[j] = Yibegin + (j - index_begin) * (Yiend - Yibegin) / (_mesh_points - index_begin - 1);
            dYi[j] = (Yiend - Yibegin) / (X_init[_mesh_points - 1] - X_init[index_begin]);
        }
        Y_init[i] = Yi;
        Y_init[i + _N_Fields] = dYi;
    }
    Y_init = transpose(Y_init);
    // cout<<Y_init.size()<<endl;
    _ODESolver.SetBoundary(X_init, Y_init);
    _X = X_init;
    _Y = Y_init;
    // DumpSolution("Initial_Guess.dat");
    _dV_replace = [&](VD y_aver) {
        VD field;
        for (int i = 0; i < _N_Fields; i++) {
            field.push_back(_Field_Basis[i] + y_aver[i] * _overall_scale[i]);
        }
        return _mod->dVtotal(field);
    };
    _d2V_replace = [&](VD y_aver) {
        VD field;
        for (int i = 0; i < _N_Fields; i++) {
            field.push_back(_Field_Basis[i] + y_aver[i] * _overall_scale[i]);
        }
        return _mod->d2Vtotal(field);
    };
}
bool DWSolver::Solve(VD &X, VVD &Y) {
    _ODESolver.SetDOF(_ODE_DOF, _N_Left_Bound, _mesh_points);

    SetInitial();
    _ODESolver.SetMaxIteration(10000);
    _ODESolver.SetConvergeCriterion(0.5, 1e-8);
    _ODESolver.SetODESystem(DIFEQ_DW, this);
    VD scales(_ODE_DOF, 1);
    _ODESolver.SetScales(scales);
    bool good = _ODESolver.SOLVDE();
    // if (!good) return good;
    VD X_Solved = _ODESolver.GetX();
    VVD Y_Solved = _ODESolver.GetY();
    // _ODESolver.SetBoundary(X_Solved,Y_Solved);
    // good = _ODESolver.SOLVDE();
    // X_Solved = _ODESolver.GetX();
    // Y_Solved = _ODESolver.GetY();
    // _X = X_Solved/_z_scale;
    _X.clear();
    _Y.clear();
    for (size_t i = 0; i < Y_Solved.size(); i++) {
        VD y;
        for (size_t j = 0; j < _N_Fields; j++) {
            y.push_back(_Field_Basis[j] + Y_Solved[i][j] * _overall_scale[j]);
        }
        _X.push_back(tan(X_Solved[i]) / _z_scale);
        _Y.push_back(y);
    }
    X = _X;
    Y = _Y;
    return good;
}
void DWSolver::SetDWODE(const Relaxation_Param relax_param, VVD &S) {
    // ! Clean S
    for (size_t i = 0; i < S.size(); i++) {
        for (size_t j = 0; j < S[i].size(); j++) {
            S[i][j] = 0;
        }
    }

    if (relax_param.k == relax_param.k_init) {
        SetDWODE_LeftBoundary(relax_param, S);
    } else if (relax_param.k > relax_param.k_final) {
        SetDWODE_RightBoundary(relax_param, S);
    } else {
        SetDWODE_Body(relax_param, S);
    }
}
void DWSolver::SetDWODE_LeftBoundary(const Relaxation_Param relax_param, VVD &S) {
    VD y1 = relax_param.y1;
    for (size_t i = 0; i < _N_Left_Bound; i++) {
        for (size_t j = 0; j < _N_Fields; j++) {
            if (i == j) {
                S[i + _N_Right_Bound][j + _ODE_DOF] = 1;
            }
        }
        S[i + _N_Right_Bound][relax_param.k_coeff] = y1[i] - (_Left_Bound[i] - _Field_Basis[i]) / _overall_scale[i];
    }
}
void DWSolver::SetDWODE_RightBoundary(const Relaxation_Param relax_param, VVD &S) {
    VD y1 = relax_param.y1;
    for (size_t i = 0; i < _N_Right_Bound; i++) {
        for (size_t j = 0; j < _N_Fields; j++) {
            if (i == j) {
                S[i][j + _ODE_DOF] = 1;
            }
        }
        S[i][relax_param.k_coeff] = y1[i] - (_Right_Bound[i] - _Field_Basis[i]) / _overall_scale[i];
    }
}
void DWSolver::SetDWODE_Body(const Relaxation_Param relax_param, VVD &S) {
    int k = relax_param.k;
    double x1 = relax_param.x1;
    double x2 = relax_param.x2;
    VD y1 = relax_param.y1;
    VD y2 = relax_param.y2;
    double dz = x2 - x1;
    double z_aver = (x1 + x2) / 2.0;
    double sz = sin(z_aver);
    double cz3 = pow(cos(z_aver), 3);
    double cz4 = pow(cos(z_aver), 4);
    VD y_aver = (y1 + y2) / 2.0;
    VD dy = y2 - y1;

    for (size_t i = 0; i < _N_Fields; i++) {
        for (size_t j = 0; j < _ODE_DOF; j++) {
            if (i == j) {
                S[i][j] = -1;
                S[i][j + _ODE_DOF] = 1;
            }
            if (i + _N_Fields == j) {
                S[i][j] = -dz / 2;
                S[i][j + _ODE_DOF] = -dz / 2;
            }
        }
        S[i][relax_param.k_coeff] = dy[i] - dz * y_aver[i + _N_Fields];
    }

    // VD field_aver;
    // for (size_t i = 0; i < _N_Fields; i++)
    // {
    //     field_aver.push_back(y_aver[i]*_overall_scale[i]);
    // }
    // VVD HM = _mod->d2Vtotal(field_aver);
    // VD dV = _mod->dVtotal(field_aver);
    VVD HM = _d2V_replace(y_aver);
    VD dV = _dV_replace(y_aver);

    for (size_t i = 0; i < _N_Fields; i++) {
        for (size_t j = 0; j < _N_Fields; j++) {
            S[i + _N_Fields][j] = -dz / 2 * HM[i][j] / _z_scale / _z_scale / _overall_scale[i] * _overall_scale[j];
            S[i + _N_Fields][j + _ODE_DOF] =
                -dz / 2 * HM[i][j] / _z_scale / _z_scale / _overall_scale[i] * _overall_scale[j];
            if (i == j) {
                S[i + _N_Fields][j + _N_Fields] = -cz4 - dz * sz * cz3;
                S[i + _N_Fields][j + _N_Fields + _ODE_DOF] = cz4 - dz * sz * cz3;
            }
        }
        S[i + _N_Fields][relax_param.k_coeff] = cz4 * dy[i + _N_Fields] - 2 * dz * sz * cz3 * y_aver[i + _N_Fields] -
                                                dz * dV[i] / _z_scale / _z_scale / _overall_scale[i];
    }
}

void DWSolver::PrintSolution() {
    cout << "The Solution is:" << endl;
    cout << "x\t";
    for (size_t i = 0; i < _N_Fields; i++) {
        cout << "y_" << i << "\t";
    }
    cout << endl;
    for (size_t i = 0; i < _X.size(); i++) {
        cout << _X[i] << "\t";
        for (size_t j = 0; j < _N_Fields; j++) {
            cout << _Y[i][j] << "\t";
        }
        cout << endl;
    }
}
void DWSolver::DumpSolution(string filename) {
    ofstream output(filename.c_str());
    // output<<"The Solution is:"<<endl;
    output << "x\t";
    for (size_t i = 0; i < _N_Fields; i++) {
        output << "y_" << i << "\t";
    }
    output << endl;
    output << scientific << setprecision(10);
    for (size_t i = 0; i < _X.size(); i++) {
        output << _X[i] << "\t";
        for (size_t j = 0; j < _N_Fields; j++) {
            output << _Y[i][j] << "\t";
        }
        output << endl;
    }
}
