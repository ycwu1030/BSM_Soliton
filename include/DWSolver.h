#ifndef DWSOLVER_H
#define DWSOLVER_H

#include "Relaxation.h"
#include "Potential.h"
#include <string>

class DWSolver
{
private:
    int _N_Fields;
    int _ODE_DOF;
    int _N_Left_Bound;
    int _N_Right_Bound;
    VD _X;
    VVD _Y;

    Relaxation _ODESolver;
    Potential *_mod;
    double _overall_scale;
    double _x_half_range;
    VD _Left_Bound;
    VD _Right_Bound;
    void SetDWODE_LeftBoundary(const Relaxation_Param relax_param, VVD &S);
    void SetDWODE_RightBoundary(const Relaxation_Param relax_param, VVD &S);
    void SetDWODE_Body(const Relaxation_Param relax_param, VVD &S);

public:
    DWSolver();
    DWSolver(Potential *mod);
    DWSolver(Potential *mod, VD Left_Bound, VD Right_Bound);
    ~DWSolver(){};

    void SetXRange(double x_range=25);
    void SetOverallScale(double overall_scale);
    void SetBoundary(VD Left_Bound, VD Right_Bound);
    void SetDWODE(const Relaxation_Param relax_param, VVD &S);

    bool Solve(VD &X, VVD &Y);

    void PrintSolution();
    void DumpSolution(std::string filename);
};



#endif //DWSOLVER_H