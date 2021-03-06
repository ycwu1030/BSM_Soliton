#ifndef CMMSOLVER_H
#define CMMSOLVER_H

#include "Relaxation.h"
#include "Potential.h"
#include "THDM_CPC.h"
#include <string>
#include <functional>

class THDMCMMSolver: public THDM_CPC
{
private:
    double g2;
    double gpp2;
    int _N_Fields;
    int _ODE_DOF;
    int _N_Left_Bound;
    int _N_Right_Bound;
    int _mesh_points;
    VD _X;
    VVD _Y;

    Relaxation _ODESolver;
    double _x_min;
    double _x_max;
    VD _Left_Bound;
    VD _Right_Bound;

    void SetODE_LeftBoundary(const Relaxation_Param relax_param, VVD &S);
    void SetODE_RightBoundary(const Relaxation_Param relax_param, VVD &S);
    void SetODE_Body(const Relaxation_Param relax_param, VVD &S);
    void SetInitial();

    // * Get Energy 
    VD GetKAIntegrand();
    VD GetKPhiIntegrand();
    VD GetVPhiIntegrand();

public:
    
    THDMCMMSolver(int mesh_points=400);
    THDMCMMSolver(VD Left_Bound, VD Right_Bound, int mesh_points = 400);
    ~THDMCMMSolver(){};

    void SetXRange(double x_min=0.1,double x_max=25);
    void SetMeshPoints(int mesh_points = 400) { _mesh_points = mesh_points;};
    void SetBoundary(VD Left_Bound, VD Right_Bound);

    bool Solve(VD &X, VVD &Y);

    void SetODE(const Relaxation_Param relax_param, VVD &S); // Used by ODE Solver, not for users.

    void GetEnergy(double &KA, double &KPhi, double &VPhi, double &KB);
    void PrintSolution();
    void DumpSolution(std::string filename);
};




#endif