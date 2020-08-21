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
    _deltam = (sqrt(3)-1.0)/2.0;
    _deltap = (sqrt(3)+1.0)/2.0;
    SetMeshPoints(mesh_points);

    SetMHL();
}

CMMonopoleSolver::CMMonopoleSolver(VD Left_Bound, VD Right_Bound, int mesh_points)
{
    _N_Fields = 4;
    _ODE_DOF = 2*_N_Fields;
    _N_Left_Bound = _N_Fields;
    _N_Right_Bound = _N_Fields;
    _x_min = 0.01;
    _x_max = 25.01;
    _deltam = (sqrt(3)-1.0)/2.0;
    _deltap = (sqrt(3)+1.0)/2.0;
    SetMeshPoints(mesh_points);

    SetMHL();    

    SetBoundary(Left_Bound,Right_Bound);
}

void CMMonopoleSolver::SetMHL(double MS)
{
    mS = MS;
    lamh = mS*mS/vev/vev;
    g2 = pow(g_weak,2);
    gpp2 = pow(g_weak,2)+pow(gp_hyper,2);
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
    VD X_sol = _ODESolver.GetX();
    VVD Y_sol = _ODESolver.GetY();
    
    // * WE Extend the solution to x->0, using the asymptotic form:
    // * y0 = c0 x^(deltam)
    // * y1 = 1 + c1 x^2
    // * y2 = c2 x
    // * y3 = y30 + c2 x + c3 x^(2*deltap)
    // * y4 = deltam c0 x^(deltam-1)
    // * y5 = 2 c1 x
    // * y6 = c2
    // * y7 = c2 + 2 deltap c3 x^(2 deltap - 1)
    double c0 = Y_sol[0][0]/pow(X_sol[0],_deltam);
    double c1 = (Y_sol[0][1]-1.0)/pow(X_sol[0],2);
    double c2 = Y_sol[0][2]/X_sol[0];
    double c3 = (Y_sol[0][3]-Y_sol[0][2]-_Left_Bound[3])/pow(X_sol[0],2*_deltap);

    VD X_ext = linspace(1e-3,_x_min,50);
    VVD Y_ext;
    for (int i = 0; i < X_ext.size(); i++)
    {
        VD Ytmp(_ODE_DOF);
        Ytmp[0] = c0*pow(X_ext[i],_deltam);
        Ytmp[1] = 1.0 + c1*pow(X_ext[i],2);
        Ytmp[2] = c2*X_ext[i];
        Ytmp[3] = _Left_Bound[3] + c2*X_ext[i] + c3*pow(X_ext[i],2*_deltap);
        Ytmp[4] = _deltam*c0*pow(X_ext[i],_deltam-1.0);
        Ytmp[5] = 2.0*c1*X_ext[i];
        Ytmp[6] = c2;
        Ytmp[7] = c2 + 2.0*_deltap*c3*pow(X_ext[i],2*_deltap-1);
        Y_ext.push_back(Ytmp);
    }
    
    _X = VD(X_ext.begin(),X_ext.end()-1);
    _Y = VVD(Y_ext.begin(),Y_ext.end()-1);

    _X.insert(_X.end(),X_sol.begin(),X_sol.end());
    _Y.insert(_Y.end(),Y_sol.begin(),Y_sol.end());

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
    // ! The left boundaries are fixed using asymptotic form at x->0
    VD y1 = relax_param.y1;
    double x = relax_param.x1;
    // for (int i = 0; i < _N_Left_Bound; i++)
    // {
    //     for (int j = 0; j < _N_Fields; j++)
    //     {
    //         if (i==j)
    //         {
    //             S[i+_N_Right_Bound][j+_ODE_DOF] = 1;
    //         }
    //     }
    //     S[i+_N_Right_Bound][relax_param.k_coeff] = y1[i] - _Left_Bound[i];
    // }
    S[4][8] = _deltam;
    S[4][12] = -_x_min;

    S[4][relax_param.k_coeff] = _deltam*y1[0] - y1[4]*_x_min;
    
    S[5][9] = 2.0;
    S[5][13] = -_x_min;

    S[5][relax_param.k_coeff] = 2.0*y1[1] - 2.0 - y1[5]*_x_min;
    
    S[6][10] = 1.0;
    S[6][14] = -_x_min;

    S[6][relax_param.k_coeff] = y1[2] - y1[6]*_x_min;

    // S[7][10] = -1.0;
    // S[7][11] = 1.0;
    // S[7][14] = _x_min/2.0/_deltap;
    // S[7][15] = -_x_min/2.0/_deltap;

    // S[7][relax_param.k_coeff] = y1[3] - y1[2] - _Left_Bound[3] - (y1[7] - y1[6])*_x_min/2.0/_deltap;
    S[7][11] = 1.0;
    S[7][14] = -_x_min + _x_min/2.0/_deltap;
    S[7][15] = -_x_min/2.0/_deltap;

    S[7][relax_param.k_coeff] = y1[3] - _Left_Bound[3] - y1[6]*_x_min + y1[6]*_x_min/2.0/_deltap - y1[7]*_x_min/2.0/_deltap; 
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
                S[i][j+_ODE_DOF] = 1.0;
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
                S[i][j] = -1.0;
                S[i][j+_ODE_DOF] = 1.0;
            }
            if (i + _N_Fields == j)
            {
                S[i][j] = -h/2.0;
                S[i][j+_ODE_DOF] = -h/2.0;
            }
        }
        S[i][relax_param.k_coeff] = dy[i] - h*ya[i+_N_Fields]; 
    }

    S[4][0] = (-2.0*ya[1]*ya[1]/xa/xa+ya[3]*ya[3]+lamh*(2.0-6.0*ya[0]*ya[0]))*h/8.0;
    S[4][1] = -h*ya[0]*ya[1]/2.0/xa/xa;
    S[4][2] = 0;
    S[4][3] = h*ya[0]*ya[3]/4.0;
    S[4][4] = h/xa-1.0;
    S[4][5] = 0;
    S[4][6] = 0;
    S[4][7] = 0;
    S[4][8] = (-2.0*ya[1]*ya[1]/xa/xa+ya[3]*ya[3]+lamh*(2.0-6.0*ya[0]*ya[0]))*h/8.0;
    S[4][9] = -h*ya[0]*ya[1]/2.0/xa/xa;
    S[4][10] = 0;
    S[4][11] = h*ya[0]*ya[3]/4.0;
    S[4][12] = h/xa+1.0;
    S[4][13] = 0;
    S[4][14] = 0;
    S[4][15] = 0;

    S[4][relax_param.k_coeff] = dy[4] + h*ya[0]*(2.0*lamh-2.0*lamh*ya[0]*ya[0]+ya[3]*ya[3])/4.0 + 2.0*h*ya[4]/xa - h*ya[0]*ya[1]*ya[1]/2.0/xa/xa;

    S[5][0] = -g2*h*ya[0]*ya[1]/4.0;
    S[5][1] = -(g2*xa*xa*ya[0]*ya[0]+12.0*ya[1]*ya[1]-4.0*xa*xa*ya[2]*ya[2]-4.0)*h/8.0/xa/xa;
    S[5][2] = h*ya[1]*ya[2];
    S[5][3] = 0;
    S[5][4] = 0;
    S[5][5] = -1.0;
    S[5][6] = 0;
    S[5][7] = 0;
    S[5][8] = -g2*h*ya[0]*ya[1]/4.0;
    S[5][9] = -(g2*xa*xa*ya[0]*ya[0]+12.0*ya[1]*ya[1]-4.0*xa*xa*ya[2]*ya[2]-4.0)*h/8.0/xa/xa;
    S[5][10] = h*ya[1]*ya[2];
    S[5][11] = 0;
    S[5][12] = 0;
    S[5][13] = 1;
    S[5][14] = 0;
    S[5][15] = 0;

    S[5][relax_param.k_coeff] = dy[5] - h*ya[1]*(g2*xa*xa*ya[0]*ya[0]-4.0*xa*xa*ya[2]*ya[2]+4.0*ya[1]*ya[1]-4.0)/4.0/xa/xa;

    S[6][0] = -g2*h*ya[0]*ya[3]/4.0;
    S[6][1] = -2.0*h*ya[1]*ya[2]/xa/xa;
    S[6][2] = -h*ya[1]*ya[1]/xa/xa;
    S[6][3] = -g2*h*ya[0]*ya[0]/8.0;
    S[6][4] = 0;
    S[6][5] = 0;
    S[6][6] = h/xa-1.0;
    S[6][7] = 0;
    S[6][8] = -g2*h*ya[0]*ya[3]/4.0;
    S[6][9] = -2.0*h*ya[1]*ya[2]/xa/xa;
    S[6][10] = -h*ya[1]*ya[1]/xa/xa;
    S[6][11] = -g2*h*ya[0]*ya[0]/8.0;
    S[6][12] = 0;
    S[6][13] = 0;
    S[6][14] = h/xa+1.0;
    S[6][15] = 0;

    S[6][relax_param.k_coeff] = dy[6] - (xa*(g2*xa*ya[0]*ya[0]*ya[3]-8.0*ya[6])+8.0*ya[2]*ya[1]*ya[1])*h/4.0/xa/xa;
    

    S[7][0] = -gpp2*h*ya[0]*ya[3]/4.0;
    S[7][1] = -2.0*h*ya[1]*ya[2]/xa/xa;
    S[7][2] = -h*ya[1]*ya[1]/xa/xa;
    S[7][3] = -gpp2*h*ya[0]*ya[0]/8.0;
    S[7][4] = 0;
    S[7][5] = 0;
    S[7][6] = 0;
    S[7][7] = h/xa-1.0;
    S[7][8] = -gpp2*h*ya[0]*ya[3]/4.0;
    S[7][9] = -2.0*h*ya[1]*ya[2]/xa/xa;
    S[7][10] = -h*ya[1]*ya[1]/xa/xa;
    S[7][11] = -gpp2*h*ya[0]*ya[0]/8.0;
    S[7][12] = 0;
    S[7][13] = 0;
    S[7][14] = 0;
    S[7][15] = h/xa+1.0;

    S[7][relax_param.k_coeff] = dy[7] - (xa*(gpp2*xa*ya[0]*ya[0]*ya[3]-8.0*ya[7])+8.0*ya[2]*ya[1]*ya[1])*h/4.0/xa/xa;
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
    for (size_t i = 0; i < _N_Fields*2; i++)
    {
        output<<"y_"<<i<<"\t";
    }
    output<<endl;
    output<<scientific<<setprecision(10);
    for (size_t i = 0; i < _X.size(); i++)
    {
        output<<_X[i]<<"\t";
        for (size_t j = 0; j < _N_Fields*2; j++)
        {
            output<<_Y[i][j]<<"\t";
        }
        output<<endl;
    }
}
VD CMMonopoleSolver::GetKAIntegrand()
{
    // ! 0: rho/rho0, 1: f, 2: A/rho0, 3: Z/rho0
    // ! 4: drho/dr/rho0^2, 5: df/dr/rho0, 6: dA/dr/rho0^2, 7: dZ/dr/rho0^2
    // ! X = r*rho0
    VD KAIntegrand(_X.size());
    double rho0 = vev;
    double rho, f, A, drho, df, r;
    for (int i = 0; i < _X.size(); i++)
    {
        r = _X[i]/rho0;
        f = _Y[i][1];
        A = _Y[i][2]*rho0;
        df = _Y[i][5]*rho0;
        A = 0;
        KAIntegrand[i] = A*A*f*f + df*df + (f*f-1)*(f*f-1)/2/r/r;
    }
    
    return KAIntegrand*4.0*Pi/g2;
}
VD CMMonopoleSolver::GetKPhiIntegrand()
{
    VD KPhiIntegrand(_X.size());
    double rho0 = vev;
    double rho, f, A, drho, df, r;
    for (int i = 0; i < _X.size(); i++)
    {
        r = _X[i]/rho0;
        drho = _Y[i][4]*rho0*rho0;
        KPhiIntegrand[i] = pow(r*drho,2);
    }
    return KPhiIntegrand*2.0*Pi;
}
VD CMMonopoleSolver::GetVPhiIntegrand()
{
    VD VPhiIntegrand(_X.size());
    double rho0 = vev;
    double rho, f, A, drho, df, r;
    for (int i = 0; i < _X.size(); i++)
    {
        r = _X[i]/rho0;
        rho = _Y[i][0]*rho0;
        VPhiIntegrand[i] = r*r*lamh*pow(rho*rho-rho0*rho0,2);
    }
    return VPhiIntegrand*Pi/2.0;
}
void CMMonopoleSolver::GetEnergy(double &KA, double &KPhi, double &VPhi, double &KB)
{
    double rho0 = vev;
    KA = Simpson(_X/rho0,GetKAIntegrand());
    KPhi = Simpson(_X/rho0,GetKPhiIntegrand());
    VPhi = Simpson(_X/rho0,GetVPhiIntegrand());
    KB = KPhi + 3*VPhi - KA;
}
