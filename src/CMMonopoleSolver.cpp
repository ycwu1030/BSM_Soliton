#include "CMMonopoleSolver.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <functional>
using namespace std;

void DIFEQ_CMMonopole(const Relaxation_Param relax_param, void *param, VVD &S)
{
    CMMonopoleSolver *solver = (CMMonopoleSolver *)param;
    solver->SetODE(relax_param,S);
}

CMMonopoleSolver::CMMonopoleSolver(int mesh_points)
{
    _N_Fields = 4;
    _ODE_DOF = 2*_N_Fields;
    _N_Left_Bound = _N_Fields;
    _N_Right_Bound = _N_Fields;
    _x_min = 0.01;
    _x_max = 25.01;
    SetMeshPoints(mesh_points);

    lamh = MHL2/2.0/vev/vev;
    g2 = pow(g_weak,2);
    gpp2 = pow(g_weak,2)+pow(gp_hyper,2);
}

CMMonopoleSolver::CMMonopoleSolver(VD Left_Bound, VD Right_Bound, int mesh_points)
{
    _N_Fields = 4;
    _ODE_DOF = 2*_N_Fields;
    _N_Left_Bound = _N_Fields;
    _N_Right_Bound = _N_Fields;
    _x_min = 0.01;
    _x_max = 25.01;
    SetMeshPoints(mesh_points);

    lamh = MHL2/2.0/vev/vev;
    g2 = pow(g_weak,2);
    gpp2 = pow(g_weak,2)+pow(gp_hyper,2);

    SetBoundary(Left_Bound,Right_Bound);
}

void CMMonopoleSolver::SetXRange(double xmin, double xmax)
{
    _x_min = xmin;
    _x_max = xmax;
}
void CMMonopoleSolver::SetBoundary(VD Left_Bound, VD Right_Bound)
{
    _Left_Bound = Left_Bound;
    _Right_Bound = Right_Bound;
}

void CMMonopoleSolver::SetInitial()
{
    VD X_init;
    VVD Y_init(_ODE_DOF);

    double dx = (_x_max-_x_min)/(_mesh_points-1);
    for (int i = 0; i < _mesh_points; i++)
    {
        X_init.push_back(_x_min + i*dx);
    }
    X_init[0] = _x_min;
    X_init.back() = _x_max;

    for (int i = 0; i < _N_Fields; i++)
    {
        VD Yi;
        VD dYi;
        double Ystep = (_Right_Bound[i]-_Left_Bound[i])/(_mesh_points-1);
        double dY_esitimate = (_Right_Bound[i]-_Left_Bound[i])/(_x_max-_x_min);
        for (int j = 0; j < _mesh_points; j++)
        {
            Yi.push_back(_Left_Bound[i]+j*Ystep);
            dYi.push_back(dY_esitimate);
        }
        Y_init[i] = Yi;
        Y_init[i+_N_Fields] = dYi;
    }
    Y_init = transpose(Y_init);
    _ODESolver.SetBoundary(X_init,Y_init);
    _X = X_init;
    _Y = Y_init;
}

bool CMMonopoleSolver::Solve(VD &X, VVD &Y)
{
    _ODESolver.SetDOF(_ODE_DOF,_N_Left_Bound,_mesh_points);

    SetInitial();
    
    _ODESolver.SetMaxIteration(20000);
    _ODESolver.SetConvergeCriterion(0.6,1e-9);
    _ODESolver.SetODESystem(DIFEQ_CMMonopole,this);
    VD scales(_ODE_DOF,1);
    _ODESolver.SetScales(scales);
    bool good = _ODESolver.SOLVDE();
    _X = _ODESolver.GetX();
    _Y = _ODESolver.GetY();
    
    X = _X;
    Y = _Y;
    return good;
}

void CMMonopoleSolver::SetODE(const Relaxation_Param relax_param, VVD &S)
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
        SetODE_LeftBoundary(relax_param, S);
    }
    else if (relax_param.k > relax_param.k_final)
    {
        SetODE_RightBoundary(relax_param,S);
    }
    else
    {
        SetODE_Body(relax_param,S);
    }
}

void CMMonopoleSolver::SetODE_LeftBoundary(const Relaxation_Param relax_param, VVD &S)
{
    VD y1 = relax_param.y1;
    for (int i = 0; i < _N_Left_Bound; i++)
    {
        for (int j = 0; j < _N_Fields; j++)
        {
            if (i==j)
            {
                S[i+_N_Right_Bound][j+_ODE_DOF] = 1;
            }
        }
        S[i+_N_Right_Bound][relax_param.k_coeff] = y1[i] - _Left_Bound[i];
    }
}

void CMMonopoleSolver::SetODE_RightBoundary(const Relaxation_Param relax_param, VVD &S)
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
        S[i][relax_param.k_coeff] = y1[i] - _Right_Bound[i];
    }
}

void CMMonopoleSolver::SetODE_Body(const Relaxation_Param relax_param, VVD &S)
{
    int k = relax_param.k;
    double x1 = relax_param.x1;
    double x2 = relax_param.x2;
    VD y1 = relax_param.y1;
    VD y2 = relax_param.y2;
    double h = x2-x1;
    double xa = (x1+x2)/2.0;
    VD ya = (y1+y2)/2.0;
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
                S[i][j] = -h/2;
                S[i][j+_ODE_DOF] = -h/2;
            }
        }
        S[i][relax_param.k_coeff] = dy[i] - h*ya[i+_N_Fields]; 
    }

    S[4][0] = (-2*ya[1]*ya[1]/xa/xa+ya[3]*ya[3]+lamh*(2-6*ya[0]*ya[0]))*h/8.;
    S[4][1] = -h*ya[0]*ya[1]/2/xa/xa;
    S[4][2] = 0;
    S[4][3] = h*ya[0]*ya[3]/4;
    S[4][4] = h/xa-1;
    S[4][5] = 0;
    S[4][6] = 0;
    S[4][7] = 0;
    S[4][8] = (-2*ya[1]*ya[1]/xa/xa+ya[3]*ya[3]+lamh*(2-6*ya[0]*ya[0]))*h/8.;
    S[4][9] = -h*ya[0]*ya[1]/2/xa/xa;
    S[4][10] = 0;
    S[4][11] = h*ya[0]*ya[3]/4;
    S[4][12] = h/xa+1;
    S[4][13] = 0;
    S[4][14] = 0;
    S[4][15] = 0;

    S[4][relax_param.k_coeff] = dy[4] + h*ya[0]*(2*lamh-2*lamh*ya[0]*ya[0]+ya[3]*ya[3])/4 + 2*h*ya[4]/xa - h*ya[0]*ya[1]*ya[1]/2/xa/xa;

    S[5][0] = -g2*h*ya[0]*ya[1]/4;
    S[5][1] = -(g2*xa*xa*ya[0]*ya[0]+12*ya[1]*ya[1]-4*xa*xa*ya[2]*ya[2]-4)*h/8/xa/xa;
    S[5][2] = h*ya[0]*ya[1];
    S[5][3] = 0;
    S[5][4] = 0;
    S[5][5] = -1;
    S[5][6] = 0;
    S[5][7] = 0;
    S[5][8] = -g2*h*ya[0]*ya[1]/4;
    S[5][9] = -(g2*xa*xa*ya[0]*ya[0]+12*ya[1]*ya[1]-4*xa*xa*ya[2]*ya[2]-4)*h/8/xa/xa;
    S[5][10] = h*ya[1]*ya[2];
    S[5][11] = 0;
    S[5][12] = 0;
    S[5][13] = 1;
    S[5][14] = 0;
    S[5][15] = 0;

    S[5][relax_param.k_coeff] = dy[5] - h*ya[1]*(g2*xa*xa*ya[0]*ya[0]-4*xa*xa*ya[2]*ya[2]+4*ya[1]*ya[1]-4)/4/xa/xa;

    S[6][0] = -g2*h*ya[0]*ya[3]/4;
    S[6][1] = -2*h*ya[1]*ya[2]/xa/xa;
    S[6][2] = -h*ya[1]*ya[1]/xa/xa;
    S[6][3] = -g2*h*ya[0]*ya[0]/8;
    S[6][4] = 0;
    S[6][5] = 0;
    S[6][6] = h/xa-1;
    S[6][7] = 0;
    S[6][8] = -g2*h*ya[0]*ya[3]/4;
    S[6][9] = -2*h*ya[1]*ya[2]/xa/xa;
    S[6][10] = -h*ya[1]*ya[1]/xa/xa;
    S[6][11] = -g2*h*ya[0]*ya[0]/8;
    S[6][12] = 0;
    S[6][13] = 0;
    S[6][14] = h/xa+1;
    S[6][15] = 0;

    S[6][relax_param.k_coeff] = dy[6] - (xa*(g2*xa*ya[0]*ya[0]*ya[3]-8*ya[6])+8*ya[2]*ya[1]*ya[1])*h/4/xa/xa;
    

    S[7][0] = -gpp2*h*ya[0]*ya[3]/4;
    S[7][1] = -2*h*ya[1]*ya[2]/xa/xa;
    S[7][2] = -h*ya[1]*ya[1]/xa/xa;
    S[7][3] = -gpp2*h*ya[0]*ya[0]/8;
    S[7][4] = 0;
    S[7][5] = 0;
    S[7][6] = 0;
    S[7][7] = h/xa-1;
    S[7][8] = -gpp2*h*ya[0]*ya[3]/4;
    S[7][9] = -2*h*ya[1]*ya[2]/xa/xa;
    S[7][10] = -h*ya[1]*ya[1]/xa/xa;
    S[7][11] = -gpp2*h*ya[0]*ya[0]/8;
    S[7][12] = 0;
    S[7][13] = 0;
    S[7][14] = 0;
    S[7][15] = h/xa+1;

    S[7][relax_param.k_coeff] = dy[7] - (xa*(gpp2*xa*ya[0]*ya[0]*ya[3]-8*ya[7])+8*ya[2]*ya[1]*ya[1])*h/4/xa/xa;
}

void CMMonopoleSolver::PrintSolution()
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
void CMMonopoleSolver::DumpSolution(string filename)
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
