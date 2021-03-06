#include "DWSolver.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <functional>
using namespace std;

void DIFEQ_DW(const Relaxation_Param relax_param, void *param, VVD &S)
{
    DWSolver *solver = (DWSolver *)param;
    solver->SetDWODE(relax_param,S);
}
DWSolver::DWSolver()
{
    _mod = nullptr;
    _z_range = -1;
    SetMeshPoints();
    // SetZRange();
}
DWSolver::DWSolver(Potential *mod, int mesh_points)
{
    _mod = mod;
    _N_Fields = _mod->GetFieldDimension();
    _ODE_DOF = 2*_N_Fields;
    _N_Left_Bound = _N_Fields;
    _N_Right_Bound = _N_Fields;
    _mesh_points = mesh_points;
    _z_range = -1;
    // SetZRange();
}
DWSolver::DWSolver(Potential *mod, VD Left_Bound, VD Right_Bound, int mesh_points)
{
    _mod = mod;
    _N_Fields = _mod->GetFieldDimension();
    _ODE_DOF = 2*_N_Fields;
    _N_Left_Bound = _N_Fields;
    _N_Right_Bound = _N_Fields;

    _z_range = -1;
    SetBoundary(Left_Bound,Right_Bound);
    SetMeshPoints(mesh_points);
}
void DWSolver::SetZRange(double z_range)
{
    _z_range = z_range;
    _x_half_range = atan(z_range);
}
void DWSolver::SetZRange()
{
    double zmax = 0;
    VD couplings = _mod->QuarticCoupling(_Left_Bound);
    for (int i = 0; i < _N_Fields; i++)
    {
        if (_Left_Bound[i]*_Right_Bound[i]<0)
        {
            double zi = 6.0/sqrt(couplings[i]/2);
            if (zi>zmax)
            {
                zmax = zi;
            }
        }
    }
    _z_range = zmax;
    _x_half_range = atan(_z_range);
}
void DWSolver::SetBoundary(VD Left_Bound, VD Right_Bound)
{
    _Left_Bound = Left_Bound;
    _Right_Bound = Right_Bound;
    SetZRange();
    VD scales;
    for (int i = 0; i < _N_Fields; i++)
    {
        if (_Left_Bound[i]*_Right_Bound[i]<0)
        {
            scales.push_back(max(abs(_Left_Bound[i]),abs(_Right_Bound[i])));
        }
        else
        {
            double center = (_Left_Bound[i] + _Right_Bound[i])/2;
            double left_left = _Left_Bound[i] - center;
            double right_left = _Right_Bound[i] - center;
            double scale = max(abs(left_left),abs(right_left));
            scale = scale<1?1:scale;
            scales.push_back(scale);
        }
    }
    SetOverallScale(scales);
}
void DWSolver::SetOverallScale(VD overall_scale)
{
    _overall_scale = abs(overall_scale);
    _z_scale = *max_element(_overall_scale.begin(),_overall_scale.end());
}
void DWSolver::SetOverallScale(double overall_scale)
{
    _overall_scale = VD(_mod->GetFieldDimension(),overall_scale);
    _z_scale = overall_scale;
}
void DWSolver::SetInitial()
{
    VD X_init;
    VVD Y_init(_ODE_DOF);
    double z_max = tan(_x_half_range);
    double dz = 2*z_max/(_mesh_points-1);
    for (int i = 0; i < _mesh_points; i++)
    {
        X_init.push_back(atan(-z_max+i*dz));
    }
    X_init[0] = -_x_half_range;
    X_init.back() = _x_half_range;
    // cout<<X_init<<endl;
    _Field_Basis.clear();
    for (int i = 0; i < _N_Fields; i++)
    {
        VD Yi;
        VD dYi;
        function<double(double)> f,df;
        double f_aver = (_Left_Bound[i] + _Right_Bound[i])/2.0;
        _Field_Basis.push_back(f_aver);
        double left_left = _Left_Bound[i] - f_aver;
        double right_left = _Right_Bound[i] - f_aver;
        if (_Left_Bound[i]*_Right_Bound[i]>0)
        {
            f = [&](double x){double r = tan(_x_half_range)/5.0; return left_left/_overall_scale[i] + (right_left - left_left)/_overall_scale[i]/(2*_x_half_range)*(x+_x_half_range) + 1/pow(cosh(tan(x)/r),2);};
            df = [&](double x){double r = tan(_x_half_range)/5.0; return (right_left - left_left)/_overall_scale[i]/(2*_x_half_range) - 2*tanh(tan(x)/r)/r/pow(cosh(tan(x)/r),2);};
        }
        else
        {
            f = [&](double x){double r = tan(_x_half_range)/5.0; return -left_left/_overall_scale[i]*tanh(tan(x)/r);};
            df = [&](double x){double r = tan(_x_half_range)/5.0; return -left_left/_overall_scale[i]/pow(cosh(tan(x)/r),2);};
        }
        // f = [&](double x){double r = _x_half_range/5.0; return left_left/_overall_scale[i] + (right_left - left_left)/_overall_scale[i]/(2*_x_half_range)*(x+_x_half_range);};
        // df = [&](double x){double r = _x_half_range/5.0; return (right_left - left_left)/_overall_scale[i]/(2*_x_half_range);};
        for (int j = 0; j < _mesh_points; j++)
        {
            Yi.push_back(f(X_init[j]));
            dYi.push_back(df(X_init[j]));
        }
        int index_end = _mesh_points*0.05;
        double Yiend = Yi[index_end];
        double Yibegin = left_left/_overall_scale[i];
        for (int j = 0; j <= index_end; j++)
        {
            Yi[j] =  Yibegin + j*(Yiend-Yibegin)/index_end;
            dYi[j] = (Yiend-Yibegin)/(X_init[index_end]-X_init[0]);
        }
        int index_begin = _mesh_points-index_end;
        Yibegin = Yi[index_begin];
        Yiend = right_left/_overall_scale[i];
        for (int j = index_begin; j < _mesh_points; j++)
        {
            Yi[j] = Yibegin + (j-index_begin)*(Yiend-Yibegin)/(_mesh_points-index_begin-1);
            dYi[j] = (Yiend-Yibegin)/(X_init[_mesh_points-1]-X_init[index_begin]);
        }
        Y_init[i] = Yi;
        Y_init[i+_N_Fields] = dYi;
    }
    Y_init = transpose(Y_init);
    // cout<<Y_init.size()<<endl;
    _ODESolver.SetBoundary(X_init,Y_init);
    _X = X_init;
    _Y = Y_init;
    // DumpSolution("Initial_Guess.dat");
    _dV_replace = [&](VD y_aver){
        VD field;
        for (int i = 0; i < _N_Fields; i++)
        {
            field.push_back(_Field_Basis[i] + y_aver[i]*_overall_scale[i]);
        }
        return _mod->dVtotal(field);
    };
    _d2V_replace = [&](VD y_aver){
        VD field;
        for (int i = 0; i < _N_Fields; i++)
        {
            field.push_back(_Field_Basis[i] + y_aver[i]*_overall_scale[i]);
        }
        return _mod->d2Vtotal(field);
    };
}
bool DWSolver::Solve(VD &X, VVD &Y)
{
    _ODESolver.SetDOF(_ODE_DOF,_N_Left_Bound,_mesh_points);
    
    SetInitial();
    _ODESolver.SetMaxIteration(10000);
    _ODESolver.SetConvergeCriterion(0.5,1e-8);
    _ODESolver.SetODESystem(DIFEQ_DW,this);
    VD scales(_ODE_DOF,1);
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
    for (size_t i = 0; i < Y_Solved.size(); i++)
    {
        VD y;
        for (size_t j = 0; j < _N_Fields; j++)
        {
            y.push_back(_Field_Basis[j]+Y_Solved[i][j]*_overall_scale[j]);
        }
        _X.push_back(tan(X_Solved[i])/_z_scale);
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
        S[i+_N_Right_Bound][relax_param.k_coeff] = y1[i] - (_Left_Bound[i]-_Field_Basis[i])/_overall_scale[i];
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
        S[i][relax_param.k_coeff] = y1[i] - (_Right_Bound[i]-_Field_Basis[i])/_overall_scale[i];
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
    double sz = sin(z_aver);
    double cz3 = pow(cos(z_aver),3);
    double cz4 = pow(cos(z_aver),4);
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
    
    // VD field_aver;
    // for (size_t i = 0; i < _N_Fields; i++)
    // {
    //     field_aver.push_back(y_aver[i]*_overall_scale[i]);
    // }
    // VVD HM = _mod->d2Vtotal(field_aver);
    // VD dV = _mod->dVtotal(field_aver);
    VVD HM = _d2V_replace(y_aver);
    VD dV = _dV_replace(y_aver);
    
    for (size_t i = 0; i < _N_Fields; i++)
    {
        for (size_t j = 0; j < _N_Fields; j++)
        {
            S[i+_N_Fields][j] = -dz/2*HM[i][j]/_z_scale/_z_scale/_overall_scale[i]*_overall_scale[j];
            S[i+_N_Fields][j+_ODE_DOF] = -dz/2*HM[i][j]/_z_scale/_z_scale/_overall_scale[i]*_overall_scale[j];
            if (i == j)
            {
                S[i+_N_Fields][j+_N_Fields] = -cz4 - dz*sz*cz3;
                S[i+_N_Fields][j+_N_Fields+_ODE_DOF] = cz4 - dz*sz*cz3;
            }
        }
        S[i+_N_Fields][relax_param.k_coeff] = cz4*dy[i+_N_Fields] - 2*dz*sz*cz3*y_aver[i+_N_Fields] - dz*dV[i]/_z_scale/_z_scale/_overall_scale[i];
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
    output<<scientific<<setprecision(10);
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