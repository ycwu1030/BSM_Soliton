#include "DWSolver.h"
#include <iostream>
#include <fstream>
using namespace std;

void DIFEQ_DW(const Relaxation_Param relax_param, void *param, VVD &S)
{
    DWSolver *solver = (DWSolver *)param;
    solver->SetDWODE(relax_param,S);
}
DWSolver::DWSolver()
{
    _mod = nullptr;
    SetXRange();
}
DWSolver::DWSolver(Potential *mod)
{
    _mod = mod;
    SetXRange();
}
DWSolver::DWSolver(Potential *mod, VD Left_Bound, VD Right_Bound)
{
    _mod = mod;
    _Left_Bound = Left_Bound;
    _Right_Bound = Right_Bound;
    SetXRange();
}
void DWSolver::SetXRange(double x_range)
{
    _x_half_range = x_range;
}
void DWSolver::SetBoundary(VD Left_Bound, VD Right_Bound)
{
    _Left_Bound = Left_Bound;
    _Right_Bound = Right_Bound;
}
void DWSolver::SetOverallScale(double overall_scale)
{
    _overall_scale = overall_scale;
}
bool DWSolver::Solve(VD &X, VVD &Y)
{
    _N_Fields = _mod->GetFieldDimension();
    _ODE_DOF = 2*_N_Fields;
    _N_Left_Bound = _N_Fields;
    _N_Right_Bound = _N_Fields;
    _ODESolver.SetDOF(_ODE_DOF,_N_Left_Bound,200);
    VD left_bound;
    VD right_bound;
    for (size_t i = 0; i < _N_Fields; i++)
    {
        left_bound.push_back(_Left_Bound[i]/_overall_scale);
        right_bound.push_back(_Right_Bound[i]/_overall_scale);
    }
    for (size_t i = 0; i < _N_Fields; i++)
    {
        left_bound.push_back(0);
        right_bound.push_back(0);
    }
    _ODESolver.SetBoundary(-_x_half_range,_x_half_range,left_bound,right_bound);
    _ODESolver.SetMaxIteration(10000);
    _ODESolver.SetConvergeCriterion(1.0,1e-8);
    _ODESolver.SetODESystem(DIFEQ_DW,this);
    VD scales(4,1);
    _ODESolver.SetScales(scales);
    bool good = _ODESolver.SOLVDE();
    if (!good) return good;
    VD X_Solved = _ODESolver.GetX();
    VVD Y_Solved = _ODESolver.GetY();
    _X = X_Solved/_overall_scale;
    _Y.clear();
    for (size_t i = 0; i < Y_Solved.size(); i++)
    {
        VD y;
        for (size_t j = 0; j < _N_Fields; j++)
        {
            y.push_back(Y_Solved[i][j]*_overall_scale);
        }
        _Y.push_back(y);
    }
    X = _X;
    Y = _Y;
    return good;
}
void DWSolver::SetDWODE(const Relaxation_Param relax_param, VVD &S)
{
    // ! Clean S
    for (size_t i = 0; i < S.size(); i++)
    {
        for (size_t j = 0; j < S[i].size(); j++)
        {
            S[i][j] = 0;
        }
    }

    if (relax_param.k == relax_param.k_init)
    {
        SetDWODE_LeftBoundary(relax_param, S);
    }
    else if (relax_param.k > relax_param.k_final)
    {
        SetDWODE_RightBoundary(relax_param,S);
    }
    else
    {
        SetDWODE_Body(relax_param,S);
    }
}
void DWSolver::SetDWODE_LeftBoundary(const Relaxation_Param relax_param, VVD &S)
{
    VD y1 = relax_param.y1;
    for (size_t i = 0; i < _N_Left_Bound; i++)
    {
        for (size_t j = 0; j < _N_Fields; j++)
        {
            if (i == j)
            {
                S[i+_N_Right_Bound][j+_ODE_DOF] = 1;
            }
        }
        S[i+_N_Right_Bound][relax_param.k_coeff] = y1[i] - _Left_Bound[i]/_overall_scale;
    }
}
void DWSolver::SetDWODE_RightBoundary(const Relaxation_Param relax_param, VVD &S)
{
    VD y1 = relax_param.y1;
    for (size_t i = 0; i < _N_Right_Bound; i++)
    {
        for (size_t j = 0; j < _N_Fields; j++)
        {
            if (i == j)
            {
                S[i][j+_ODE_DOF] = 1;
            }
        }
        S[i][relax_param.k_coeff] = y1[i] - _Right_Bound[i]/_overall_scale;
    }
}
void DWSolver::SetDWODE_Body(const Relaxation_Param relax_param, VVD &S)
{
    int k = relax_param.k;
    double x1 = relax_param.x1;
    double x2 = relax_param.x2;
    VD y1 = relax_param.y1;
    VD y2 = relax_param.y2;
    double dz = x2-x1;
    double z_aver = (x1+x2)/2.0;
    VD y_aver = (y1+y2)/2.0;
    VD dy = y2 - y1;

    for (size_t i = 0; i < _N_Fields; i++)
    {
        for (size_t j = 0; j < _ODE_DOF; j++)
        {
            if (i == j)
            {
                S[i][j] = -1;
                S[i][j+_ODE_DOF] = 1;
            }
            if (i + _N_Fields == j)
            {
                S[i][j] = -dz/2;
                S[i][j+_ODE_DOF] = -dz/2;
            }
        }
        S[i][relax_param.k_coeff] = dy[i] - dz*y_aver[i+_N_Fields]; 
    }
    
    VD field_aver;
    for (size_t i = 0; i < _N_Fields; i++)
    {
        field_aver.push_back(y_aver[i]*_overall_scale);
    }
    VVD HM = _mod->d2Vtotal(field_aver,_overall_scale);
    VD dV = _mod->dVtotal(field_aver,_overall_scale);
    
    for (size_t i = 0; i < _N_Fields; i++)
    {
        for (size_t j = 0; j < _N_Fields; j++)
        {
            S[i+_N_Fields][j] = -dz/2*HM[i][j];
            S[i+_N_Fields][j+_ODE_DOF] = -dz/2*HM[i][j];
            if (i == j)
            {
                S[i+_N_Fields][j+_N_Fields] = -1;
                S[i+_N_Fields][j+_N_Fields+_ODE_DOF] = 1;
            }
        }
        S[i+_N_Fields][relax_param.k_coeff] = dy[i+_N_Fields] - dz*dV[i];
    }
}

void DWSolver::PrintSolution()
{
    cout<<"The Solution is:"<<endl;
    cout<<"x\t";
    for (size_t i = 0; i < _N_Fields; i++)
    {
        cout<<"y_"<<i<<"\t";
    }
    cout<<endl;
    for (size_t i = 0; i < _X.size(); i++)
    {
        cout<<_X[i]<<"\t";
        for (size_t j = 0; j < _N_Fields; j++)
        {
            cout<<_Y[i][j]<<"\t";
        }
        cout<<endl;
    }
}
void DWSolver::DumpSolution(string filename)
{
    ofstream output(filename.c_str());
    // output<<"The Solution is:"<<endl;
    output<<"x\t";
    for (size_t i = 0; i < _N_Fields; i++)
    {
        output<<"y_"<<i<<"\t";
    }
    output<<endl;
    for (size_t i = 0; i < _X.size(); i++)
    {
        output<<_X[i]<<"\t";
        for (size_t j = 0; j < _N_Fields; j++)
        {
            output<<_Y[i][j]<<"\t";
        }
        output<<endl;
    }
}