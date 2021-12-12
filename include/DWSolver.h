#ifndef DWSOLVER_H
#define DWSOLVER_H

#include <functional>
#include <string>

#include "Basic_Model.h"
#include "Potential.h"
#include "Relaxation.h"

namespace BSM_Soliton {
class DomainWallSolver : public RelaxationODE {
    // * DOF convention:
    // * internally solved in terms of q from 0 to 1
    // * 0: t
    // * 1 to n: Field
    // * n+1 to 2n: dField/dz
    // * 2n+1: z_min
    // * 2n+2: z_max
    // * 2n+3: C
    // * z = z_min + t*(z_max-z_min)
    // * dt/dq = C/F
public:
    DomainWallSolver(BaseModel *mod);
    void Set_Allocation_Parameters(double ratio = 1.0) { F_ratio = ratio; }
    void Set_Mesh_Size(int MeshSize = 401) { Mesh_Size = MeshSize; }
    // void Set_Boundaries(const VD &field_at_left, const VD &field_at_right);

    bool Solve(const VD &field_at_left, const VD &field_at_right);
    bool Solve(const VVD &fields_grid);
    virtual void dYdX(MeshPoint &point) override final;
    virtual void Left_Boundary_Constraints(MeshPoint &point) override final;
    virtual void Right_Boundary_Constraints(MeshPoint &point) override final;

private:
    BaseModel *mod;
    Relaxation *solver;

    int Field_Space_Dim;
    int Mesh_Size;

    // F = 1/(1+F_ratio) + F_ratio/(1+F_ratio)*\sum_i (Y_{i+n}/Y_i)^2
    double F_ratio;

    double F_value;
    VD dF_value;
    void Calc_F_Variables(MeshPoint &point);

    VD Fields_at_Left;
    VD Fields_at_Right;
};
}  // namespace BSM_Soliton

class DWSolver {
private:
    int _N_Fields;
    int _ODE_DOF;
    int _N_Left_Bound;
    int _N_Right_Bound;
    int _mesh_points;
    VD _X;
    VVD _Y;

    Relaxation _ODESolver;
    Potential *_mod;
    VD _overall_scale;
    double _z_scale;
    double _z_range;
    double _x_half_range;
    VD _Left_Bound;
    VD _Right_Bound;
    VD _Field_Basis;
    void SetOverallScale(double overall_scale);

    void SetDWODE_LeftBoundary(const Relaxation_Param relax_param, VVD &S);
    void SetDWODE_RightBoundary(const Relaxation_Param relax_param, VVD &S);
    void SetDWODE_Body(const Relaxation_Param relax_param, VVD &S);
    void SetInitial();
    std::function<VD(VD)> _dV_replace;
    std::function<VVD(VD)> _d2V_replace;

public:
    DWSolver();
    DWSolver(Potential *mod, int mesh_points = 400);
    DWSolver(Potential *mod, VD Left_Bound, VD Right_Bound, int mesh_points = 400);
    ~DWSolver(){};

    void SetZRange();  // Set the z range automatically according to the potential
    void SetZRange(double z_range);
    void SetMeshPoints(int mesh_points = 400) { _mesh_points = mesh_points; };
    void SetBoundary(VD Left_Bound, VD Right_Bound);
    void SetOverallScale(VD overall_scale);

    bool Solve(VD &X, VVD &Y);

    void PrintSolution();
    void DumpSolution(std::string filename);

    void SetDWODE(const Relaxation_Param relax_param, VVD &S);  // Used by ODE solver, not for users
};

#endif  // DWSolver_H
