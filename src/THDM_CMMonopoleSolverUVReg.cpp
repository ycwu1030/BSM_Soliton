#include "THDM_CMMonopoleSolverUVReg.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <functional>
using namespace std;

void DIFEQ_THDMCMMonopoleUV(const Relaxation_Param relax_param, void *param, VVD &S)
{
    THDMCMMSolverUV *solver = (THDMCMMSolverUV *)param;
    solver->SetODE(relax_param,S);
}

THDMCMMSolverUV::THDMCMMSolverUV(int mesh_points):THDM_CPC()
{
    _N_Fields = 5;
    _ODE_DOF = 2*_N_Fields;
    _N_Left_Bound = _N_Fields;
    _N_Right_Bound = _N_Fields;
    _x_min = 0.01;
    _x_max = 25.01;
    SetMeshPoints(mesh_points);

    Set_EW_Parameters();
    Set_Physical_Parameters();
    ExtendtoZero();
}

THDMCMMSolverUV::THDMCMMSolverUV(VD Left_Bound, VD Right_Bound, int mesh_points):THDM_CPC()
{
    _N_Fields = 4;
    _ODE_DOF = 2*_N_Fields;
    _N_Left_Bound = _N_Fields;
    _N_Right_Bound = _N_Fields;
    _x_min = 0.01;
    _x_max = 25.01;
    SetMeshPoints(mesh_points);

    Set_EW_Parameters();
    Set_Physical_Parameters();    
    ExtendtoZero();
    SetBoundary(Left_Bound,Right_Bound);
}
void THDMCMMSolverUV::Set_EW_Parameters()
{
    g2 = pow(g_weak,2);
    gp2 = pow(gp_hyper,2);
    gpp2 = pow(g_weak,2)+pow(gp_hyper,2);
    SW2 = pow(sin(thetaW),2);
}
void THDMCMMSolverUV::SetXRange(double xmin, double xmax)
{
    _x_min = xmin;
    _x_max = xmax;
}
void THDMCMMSolverUV::SetBoundary(VD Left_Bound, VD Right_Bound)
{
    _Left_Bound = Left_Bound;
    _Right_Bound = Right_Bound;
}

void THDMCMMSolverUV::SetInitial()
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

bool THDMCMMSolverUV::Solve(VD &X, VVD &Y)
{
    _ODESolver.SetDOF(_ODE_DOF,_N_Left_Bound,_mesh_points);

    SetInitial();
    
    _ODESolver.SetMaxIteration(20000);
    _ODESolver.SetConvergeCriterion(0.6,1e-9);
    _ODESolver.SetODESystem(DIFEQ_THDMCMMonopoleUV,this);
    VD scales(_ODE_DOF,1);
    _ODESolver.SetScales(scales);
    bool good = _ODESolver.SOLVDE();
    _X = _ODESolver.GetX();
    _Y = _ODESolver.GetY();
    
    X = _X;
    Y = _Y;
    return good;
}

void THDMCMMSolverUV::SetODE(const Relaxation_Param relax_param, VVD &S)
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

void THDMCMMSolverUV::SetODE_LeftBoundary(const Relaxation_Param relax_param, VVD &S)
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

void THDMCMMSolverUV::SetODE_RightBoundary(const Relaxation_Param relax_param, VVD &S)
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

void THDMCMMSolverUV::SetODE_Body(const Relaxation_Param relax_param, VVD &S)
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
    VD dV = dVtotal({ya[0]*_vev,ya[1]*_vev},_vev);
    VVD d2V = d2Vtotal({ya[0]*_vev,ya[1]*_vev},_vev);

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

    S[5][0] = (-2.0*(1.0+_gamma)*ya[2]*ya[2]/xa/xa+ya[4]*ya[4]-4.0*d2V[0][0])*h/8.0;
    S[5][1] = -h*d2V[0][1]/2.0;
    S[5][2] = -h*(1.0+_gamma)*ya[0]*ya[2]/2.0/xa/xa;
    S[5][3] = 0;
    S[5][4] = h*ya[0]*ya[4]/4.0;
    S[5][5] = h/xa-1.0;
    S[5][6] = 0;
    S[5][7] = 0;
    S[5][8] = 0;
    S[5][9] = 0;
    S[5][10] = (-2.0*(1.0+_gamma)*ya[2]*ya[2]/xa/xa+ya[4]*ya[4]-4.0*d2V[0][0])*h/8.0;
    S[5][11] = -h*d2V[0][1]/2.0;
    S[5][12] = -h*(1.0+_gamma)*ya[0]*ya[2]/2.0/xa/xa;
    S[5][13] = 0;
    S[5][14] = h*ya[0]*ya[4]/4.0;
    S[5][15] = h/xa+1.0;
    S[5][16] = 0;
    S[5][17] = 0;
    S[5][18] = 0;
    S[5][19] = 0;

    S[5][relax_param.k_coeff] = dy[5] - h*dV[0] - h*(1.0+_gamma)*ya[0]*ya[2]*ya[2]/2.0/xa/xa + h*ya[0]*ya[4]*ya[4]/4.0 + 2.0*h*ya[5]/xa;

    S[6][0] = -h*d2V[0][1]/2.0;
    S[6][1] = (-2.0*(1.0+_gamma)*ya[2]*ya[2]/xa/xa+ya[4]*ya[4]-4.0*d2V[1][1])*h/8.0;
    S[6][2] = -h*(1.0+_gamma)*ya[1]*ya[2]/2.0/xa/xa;
    S[6][3] = 0;
    S[6][4] = h*ya[1]*ya[4]/4.0;
    S[6][5] = 0;
    S[6][6] = h/xa-1.0;
    S[6][7] = 0;
    S[6][8] = 0;
    S[6][9] = 0;
    S[6][10] = -h*d2V[0][1]/2.0;
    S[6][11] = (-2.0*(1.0+_gamma)*ya[2]*ya[2]/xa/xa+ya[4]*ya[4]-4.0*d2V[1][1])*h/8.0;
    S[6][12] = -h*(1.0+_gamma)*ya[1]*ya[2]/2.0/xa/xa;
    S[6][13] = 0;
    S[6][14] = h*ya[1]*ya[4]/4.0;
    S[6][15] = 0;
    S[6][16] = h/xa+1.0;
    S[6][17] = 0;
    S[6][18] = 0;
    S[6][19] = 0;

    S[6][relax_param.k_coeff] = dy[6] - h*dV[1] - h*(1.0+_gamma)*ya[1]*ya[2]*ya[2]/2.0/xa/xa + h*ya[1]*ya[4]*ya[4]/4.0 + 2.0*h*ya[6]/xa;


    S[7][0] = -g2*(1.0+_gamma)*h*ya[0]*ya[2]/4.0;
    S[7][1] = -g2*(1.0+_gamma)*h*ya[1]*ya[2]/4.0;
    S[7][2] = h*(_f0*_f0*(4.0*ya[3]*ya[3]*xa*xa-g2*(1.0+_gamma)*(ya[0]*ya[0]+ya[1]*ya[1])*xa*xa+4.0+4.0*_alpha)-12.0*(1.0+_alpha)*ya[2]*ya[2])/_f0/_f0/xa/xa/8.0;
    S[7][3] = h*ya[2]*ya[3];
    S[7][4] = 0;
    S[7][5] = 0;
    S[7][6] = 0;
    S[7][7] = -1;
    S[7][8] = 0;
    S[7][9] = 0;
    S[7][10] = -g2*(1.0+_gamma)*h*ya[0]*ya[2]/4.0;
    S[7][11] = -g2*(1.0+_gamma)*h*ya[1]*ya[2]/4.0;
    S[7][12] = h*(_f0*_f0*(4.0*ya[3]*ya[3]*xa*xa-g2*(1.0+_gamma)*(ya[0]*ya[0]+ya[1]*ya[1])*xa*xa+4.0+4.0*_alpha)-12.0*(1.0+_alpha)*ya[2]*ya[2])/_f0/_f0/xa/xa/8.0;
    S[7][13] = h*ya[2]*ya[3];
    S[7][14] = 0;
    S[7][15] = 0;
    S[7][16] = 0;
    S[7][17] = 1;
    S[7][18] = 0;
    S[7][19] = 0;

    S[7][relax_param.k_coeff] = dy[7] - h*(1.0+_gamma)*g2*(ya[0]*ya[0]+ya[1]*ya[1])*ya[2]/4.0 + h*ya[2]*(1.0+_alpha)*(_f0*_f0-ya[2]*ya[2])/_f0/_f0/xa/xa + h*ya[2]*ya[3]*ya[3];

    S[8][0] = -g2*h*ya[0]*ya[4]/4.0;
    S[8][1] = -g2*h*ya[1]*ya[4]/4.0;
    S[8][2] = -2.0*h*ya[2]*ya[3]/xa/xa;
    S[8][3] = -h*ya[2]*ya[2]/xa/xa;
    S[8][4] = -g2*h*(ya[0]*ya[0]+ya[1]*ya[1])/8.0;
    S[8][5] = 0;
    S[8][6] = 0;
    S[8][7] = 0;
    S[8][8] = h/xa-1.0;
    S[8][9] = 0;
    S[8][10] = -g2*h*ya[0]*ya[4]/4.0;
    S[8][11] = -g2*h*ya[1]*ya[4]/4.0;
    S[8][12] = -2.0*h*ya[2]*ya[3]/xa/xa;
    S[8][13] = -h*ya[2]*ya[2]/xa/xa;
    S[8][14] = -g2*h*(ya[0]*ya[0]+ya[1]*ya[1])/8.0;
    S[8][15] = 0;
    S[8][16] = 0;
    S[8][17] = 0;
    S[8][18] = h/xa+1.0;
    S[8][19] = 0;

    S[8][relax_param.k_coeff] = dy[8] - (xa*(g2*xa*(ya[0]*ya[0]+ya[1]*ya[1])*ya[4]-8.0*ya[8])+8.0*ya[3]*ya[2]*ya[2])*h/4.0/xa/xa;
    

    S[9][0] = -gpp2*h*ya[0]*ya[4]/4.0;
    S[9][1] = -gpp2*h*ya[1]*ya[4]/4.0;
    S[9][2] = -2.0*h*ya[2]*ya[3]/xa/xa;
    S[9][3] = -h*ya[2]*ya[2]/xa/xa;
    S[9][4] = -gpp2*h*(ya[0]*ya[0]+ya[1]*ya[1])/8.0;
    S[9][5] = 0;
    S[9][6] = 0;
    S[9][7] = 0;
    S[9][8] = 0;
    S[9][9] = h/xa-1.0;
    S[9][10] = -gpp2*h*ya[0]*ya[4]/4.0;
    S[9][11] = -gpp2*h*ya[1]*ya[4]/4.0;
    S[9][12] = -2.0*h*ya[2]*ya[3]/xa/xa;
    S[9][13] = -h*ya[2]*ya[2]/xa/xa;
    S[9][14] = -gpp2*h*(ya[0]*ya[0]+ya[1]*ya[1])/8.0;
    S[9][15] = 0;
    S[9][16] = 0;
    S[9][17] = 0;
    S[9][18] = 0;
    S[9][19] = h/xa+1.0;

    S[9][relax_param.k_coeff] = dy[9] - (xa*(gpp2*xa*(ya[0]*ya[0]+ya[1]*ya[1])*ya[4]-8.0*ya[9])+8.0*ya[3]*ya[2]*ya[2])*h/4.0/xa/xa;
}

void THDMCMMSolverUV::PrintSolution()
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
void THDMCMMSolverUV::DumpSolution(string filename)
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

VD THDMCMMSolverUV::GetE0Integrand()
{
    // ! 0: rho1/rho0, 1: rho2/rho0, 2: f, 3: A/rho0, 4: Z/rho0
    // ! 5: drho1/dr/rho0^2, 6: drho2/dr/rho0^2, 7: df/dr/rho0, 8: dA/dr/rho0^2, 9: dZ/dr/rho0^2
    // ! X = r*rho0
    VD E0Integrand(_X.size());
    double rho0 = vev;
    double f, r;
    for (int i = 0; i < _X.size(); i++)
    {
        r = _X[i]/rho0;
        f = _Y[i][1];
        E0Integrand[i] = pow(f*f-_f0*_f0,2)/pow(_f0,4)/SW2/r/r;
    }
    
    return E0Integrand*2.0*Pi/g2;
}
VD THDMCMMSolverUV::GetE1Integrand()
{
    VD E1Integrand(_X.size());
    double rho0 = _vev;
    double rho1, rho2, f, A, Z, B, drho1, drho2, df, dA, dZ, dB, r;
    for (int i = 0; i < _X.size(); i++)
    {
        r = _X[i]/rho0;
        rho1 = _Y[i][0]*rho0;
        rho2 = _Y[i][1]*rho0;
        f = _Y[i][2];
        A = _Y[i][3]*rho0;
        Z = _Y[i][4]*rho0;
        B = A-Z;
        drho1 = _Y[i][5]*rho0*rho0;
        drho2 = _Y[i][6]*rho0*rho0;
        df = _Y[i][7]*rho0;
        dA = _Y[i][8]*rho0*rho0;
        dZ = _Y[i][9]*rho0*rho0;
        dB = dA - dZ;
        E1Integrand[i] = pow(r*drho1,2) + pow(r*drho2,2) + 2.0/g2*(df*df+f*f*A*A) + 1.0/g2*pow(r*dA,2) + 1.0/gp2*pow(r*dB,2) + (1.0+_gamma)*(rho1*rho1+rho2*rho2)*f*f/2.0 + r*r*Z*Z*(rho1*rho1+rho2*rho2)/4.0 + 2.0*r*r*Vtotal({rho1,rho2});
    }
    return E1Integrand*2.0*Pi;
}
void THDMCMMSolverUV::GetEnergy(double &E0, double &E1)
{
    double rho0 = _vev;
    E0 = Simpson(_X/rho0,GetE0Integrand());
    E1 = Simpson(_X/rho0,GetE1Integrand());
}