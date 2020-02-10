/*
 * @Description  : 
 * @Author       : Yongcheng Wu
 * @Date         : 2020-01-27 14:19:34
 * @LastEditors  : Yongcheng Wu
 * @LastEditTime : 2020-02-09 23:37:02
 */
#include <iostream>
#include <cmath>
#include "SM_cxSM.h"

using namespace std;
using namespace Eigen;

bool CloseQ(double x1, double x2)
{
    return abs(x1-x2)<1e-2;
}

SM::SM()
{
    alpha = 1.0/alpha1;
    ee = sqrt(4*Pi*alpha);
    vev = pow(sqrt(2)*GF,-0.5);
    double A = sqrt(Pi*alpha)*vev;
    thetaW = asin(2*A/MZ)/2.0;
    MW = A/sin(thetaW);
    MW2 = MW*MW;
    g_weak = ee/sin(thetaW);
    gp_hyper = ee/cos(thetaW);
    yt = sqrt(2)*MT/vev;
#ifdef DEBUG
    cout<<"A:  "<<A<<endl;
    cout<<"vev: "<<vev<<endl;
    cout<<"thetaW: "<<thetaW<<endl;
    cout<<"MW: "<<MW<<endl;
    cout<<"yt: "<<yt<<endl;
    cout<<"g_weak: "<<g_weak<<endl;
    cout<<"gp_hyper: "<<gp_hyper<<endl;
#endif
}

CXSM::CXSM()
{
    SetInput();
}
CXSM::CXSM(double VSin, double MHHin, double MHAin, double thetain)
{
    SetInput(VSin, MHHin, MHAin, thetain);
}
void CXSM::SetInput(double VSin, double MHHin, double MHAin, double thetain)
{
    VS = VSin;
    MHH = MHHin;
    MHH2 = MHH*MHH;
    MHA = MHAin;
    MHA2 = MHA*MHA;
    theta = thetain;

    double cth = cos(theta);
    double sth = sin(theta);
    double cth2 = cth*cth;
    double sth2 = sth*sth;

    del1 = -4*sqrt(2)*MHA2*VS/vev/vev;
    lam = (MHL2*cth2+MHH2*sth2)/2/vev/vev;
    del2 = (4*(MHL2-MHH2)*sth*cth-sqrt(2)*del1*vev)/2/vev/VS;
    d2 = (8*MHL2*VS*sth2+8*MHH2*VS*cth2+sqrt(2)*del1*vev*vev)/4/pow(VS,3);
    mu2 = -(4*lam*vev*vev+del2*VS*VS+sqrt(2)*del1*VS)/4;
    b2 = -(2*d2*pow(VS,3)+vev*vev*(sqrt(2)*del1+2*del2*VS))/4/VS;

    Solved = false;
}

double CXSM::Vtot(double phiH, double phiS, double phiA)
{
    double vtot = 0;
    double h2 = phiH*phiH;
    double s2 = phiS*phiS;
    double A2 = phiA*phiA;
    double h4 = h2*h2;
    vtot += mu2/2*h2+b2/4*(A2+s2);
    vtot += del1/4/sqrt(2)*h2*phiS;
    vtot += lam/4*h4 + (A2+s2)*(d2*(A2+s2)+2*del2*h2)/16;
    return vtot;
}
double CXSM::Get_V0_global()
{
    return Vtot(vev,VS,0);
}
double CXSM::dVdH(double phiH, double phiS, double phiA)
{
    double dvtot = 0;
    double h2 = phiH*phiH;
    double s2 = phiS*phiS;
    double A2 = phiA*phiA;
    dvtot = phiH*(del2*A2+4*(h2*lam+mu2)+del2*s2+sqrt(2)*del1*phiS)/4.0;
    return dvtot;
}
double CXSM::dVdS(double phiH, double phiS, double phiA)
{
    double dvtot = 0;
    double h2 = phiH*phiH;
    double s2 = phiS*phiS;
    double A2 = phiA*phiA;
    dvtot = (2*d2*phiS*(A2+s2)+4*b2*phiS+h2*(sqrt(2)*del1+2*del2*phiS))/8.0;
    return dvtot;
}
double CXSM::dVdA(double phiH, double phiS, double phiA)
{
    double dvtot = 0;
    double h2 = phiH*phiH;
    double s2 = phiS*phiS;
    double A2 = phiA*phiA;
    dvtot = phiA*(A2*d2+2*b2+d2*s2+del2*h2)/4.0;
    return dvtot;
}

bool CXSM::CheckStability()
{
    return (lam>0)&&(d2>0)&&((lam*d2>del2*del2)||del2>0);
}

bool CXSM::CheckUnitarity(double MAX)
{
    double EigenA0[4];
    EigenA0[0] = d2/2;
    EigenA0[1] = 2*lam;
    EigenA0[2] = (d2+6*lam-sqrt(d2*d2-12*d2*lam+2*del2*del2+36*lam*lam))/2;
    EigenA0[3] = (d2+6*lam+sqrt(d2*d2-12*d2*lam+2*del2*del2+36*lam*lam))/2;
    bool good=true;
    for (int i = 0; i < 4; ++i)
    {
        good*=(abs(EigenA0[i])<MAX*16.0*Pi);
        if (!good)
        {
            return good;
        }
    }

    return good;
}

void CXSM::PrintPotentialParameter()
{
    cout<<"Potential Parameter: "<<endl;
    cout<<"mu2: "<<mu2<<endl;
    cout<<"b2: "<<b2<<endl;
    cout<<"lam: "<<lam<<endl;
    cout<<"del2: "<<del2<<endl;
    cout<<"d2: "<<d2<<endl;
    cout<<"del1: "<<del1<<endl;
}

void CXSM::SolveCubicEquation(double A[4], double *results, int &NSolution)
{
    NSolution = 0;
    Solver.Solve(A[3],A[2],A[1],A[0]);
    if (Solver.STATE==ONEREAL||Solver.STATE==THREEEQUALREAL)
    {
        results[NSolution]=Solver.SOLUTIONS[0];
        NSolution++;
    }
    else if (Solver.STATE==ONETWO)
    {
        for (int i = 0; i < 2; ++i)
        {
            results[NSolution]=Solver.SOLUTIONS[i];
            NSolution++;
        }
    }
    else if (Solver.STATE==THREEREAL)
    {
        for (int i = 0; i < 3; ++i)
        {
            results[NSolution]=Solver.SOLUTIONS[i];
            NSolution++;
        }
    }
    else
    {
#ifndef DEBUG
        std::cout<<"Warning in Cubic Solver 1"<<endl;
#endif
    }
}

void CXSM::FindLocalMinima()
{
    if (Solved)
    {
        return;
    }
    
    NLocalMinima=0; // Track we are getting which solution
    NLocalExtreme=0;
    IndexInput = -1;
// ! A. For phiH = 0

    // ! A-1. First solution would be phiH = 0, phiS = 0, phiA = 0
    localH[NLocalExtreme]=0;
    localS[NLocalExtreme]=0;
    localA[NLocalExtreme]=0;
    Vtotal[NLocalExtreme]=Vtot(localH[NLocalExtreme],localS[NLocalExtreme],localA[NLocalExtreme]);
    if ((LocalMinimaQ[NLocalExtreme]=CheckHessianMatrix(localH[NLocalExtreme],localS[NLocalExtreme],localA[NLocalExtreme])))
    {
        MinimaIndex[NLocalMinima++]=NLocalExtreme;
    }
    ++NLocalExtreme;

    // ! A-2. Second solution would be phiH = 0 and any point along phiS^2 + phiA^2 = -2*b2/d2; But the whole valley is degenerate, so we can just pick one point as phiS > 0 and phiA = 0
    if (b2/d2 < 0)
    {
        localH[NLocalExtreme]=0;
        localS[NLocalExtreme]=sqrt(-2*b2/d2);
        localA[NLocalExtreme]=0;
        Vtotal[NLocalExtreme]=Vtot(localH[NLocalExtreme],localS[NLocalExtreme],localA[NLocalExtreme]);
        if ((LocalMinimaQ[NLocalExtreme]=CheckHessianMatrix(localH[NLocalExtreme],localS[NLocalExtreme],localA[NLocalExtreme])))
        {
            MinimaIndex[NLocalMinima++]=NLocalExtreme;
        }
        ++NLocalExtreme;
    }

// ! B. For phiH != 0, then phiH^2 = - (phiA^2*del2+4*mu2+del2*phiS^2+sqrt(2)*del1*phiS)/4/lam;

    // ! B-1. For phiA = 0;
    // ! Then phiS satisfies:
    // ! -(del2^2-4*d2*lam)/16/lam phiS^3 - 3*sqrt(2)*del1*del2/32/lam phiS^2 - (del1^2+4*del2*mu2-8*b2*lam)/16/lam phiS - sqrt(2)*del1*mu2/8/lam = 0;
    int NCurrent=0;
    double results[3];
    double stemp;
    double h2temp;
    double AA[4] = {-sqrt(2)*del1*mu2/8/lam,-(del1*del1+4*del2*mu2-8*b2*lam)/16/lam,-3*sqrt(2)*del1*del2/32/lam,-(del2*del2-4*d2*lam)/16/lam};
    SolveCubicEquation(AA,results,NCurrent);
    for (int i = 0; i < NCurrent; ++i)
    {
        stemp = results[i];
        h2temp = -(4*mu2+del2*stemp*stemp+sqrt(2)*del1*stemp)/4/lam;
        if (h2temp < 0)
        {
            continue;
        }
        localH[NLocalExtreme]=sqrt(h2temp);
        localS[NLocalExtreme]=stemp;
        localA[NLocalExtreme]=0;
        Vtotal[NLocalExtreme]=Vtot(localH[NLocalExtreme],localS[NLocalExtreme],localA[NLocalExtreme]);
        if ((LocalMinimaQ[NLocalExtreme]=CheckHessianMatrix(localH[NLocalExtreme],localS[NLocalExtreme],localA[NLocalExtreme])))
        {
            MinimaIndex[NLocalMinima++]=NLocalExtreme;
        }
        if (CloseQ(sqrt(h2temp),vev)&&CloseQ(stemp,VS))
        {
            IndexInput = NLocalExtreme;
        }
        ++NLocalExtreme;
    }
    if (IndexInput < 0)
    {
        cout<<"Error in getting local extreme point: Didn't find local extreme point corresponding to the input."<<endl;
    }
    

    // ! B-2. phiA != 0;
    // ! This is the most non-trivial case, where phiH != 0, phiS != 0 and phiH != 0
    double A2temp;
    if (del1 != 0)
    {
        stemp = sqrt(2)*(b2*del2 - 2*d2*mu2)/del1/d2;
        A2temp = (-8*b2*lam+4*del2*mu2+stemp*stemp*(del2*del2-4*d2*lam)+sqrt(2)*del1*del2*stemp)/(4*d2*lam-del2*del2);
        if (A2temp >= 0)
        {
            h2temp = -(A2temp*del2+4*mu2+del2*stemp*stemp+sqrt(2)*del1*stemp)/4/lam;
            if (h2temp > 0)
            {
                localH[NLocalExtreme]=sqrt(h2temp);
                localS[NLocalExtreme]=stemp;
                localA[NLocalExtreme]=sqrt(A2temp);
                Vtotal[NLocalExtreme]=Vtot(localH[NLocalExtreme],localS[NLocalExtreme],localA[NLocalExtreme]);
                if ((LocalMinimaQ[NLocalExtreme]=CheckHessianMatrix(localH[NLocalExtreme],localS[NLocalExtreme],localA[NLocalExtreme])))
                {
                    MinimaIndex[NLocalMinima++]=NLocalExtreme;
                }
                ++NLocalExtreme;
            }         
        }
    }
    Solved = true;
}


void CXSM::PrintLocalMinima()
{
    FindLocalMinima();
    cout<<"The Local Extreme points are: "<<endl;
    cout<<"ID\tVH\tVS\tVA\tLocalMinimumQ\tVtot"<<endl;
    for (size_t i = 0; i < NLocalExtreme; i++)
    {
        cout<<i<<"\t"<<localH[i]<<"\t"<<localS[i]<<"\t"<<localA[i]<<"\t"<<LocalMinimaQ[i]<<"\t"<<Vtotal[i]<<endl;
    }
}


void CXSM::GetHessian(double phiH, double phiS, double phiA)
{
    HessianMatrix(0,0)=(del2*phiA*phiA+sqrt(2)*phiS*del1+del2*phiS*phiS+4*(3*lam*phiH*phiH+mu2))/4.0;
    HessianMatrix(0,1)=phiH*(sqrt(2)*del1+2*phiS*del2)/4.0;
    HessianMatrix(0,2)=phiA*phiH*del2/2.0;
    HessianMatrix(1,0)=HessianMatrix(0,1);
    HessianMatrix(1,1)=(d2*phiA*phiA+3*d2*phiS*phiS+2*b2+del2*phiH*phiH)/4.0;
    HessianMatrix(1,2)=phiA*phiS*d2/2.0;
    HessianMatrix(2,0)=HessianMatrix(0,2);
    HessianMatrix(2,1)=HessianMatrix(1,2);
    HessianMatrix(2,2)=(3*d2*phiA*phiA+d2*phiS*phiS+2*b2+del2*phiH*phiH)/4.0;
}

bool CXSM::CheckHessianMatrix(double phiH, double phiS, double phiA)
{
    GetHessian(phiH, phiS, phiA);
    SelfAdjointEigenSolver<Matrix3d> eigensolver(HessianMatrix);
    if (eigensolver.info() != Success) return false;
#if DEBUG
    cout<<"Eigen Values: "<<endl;
    cout<<eigensolver.eigenvalues()<<endl;
    cout<<"--"<<endl;
#endif
    return ((eigensolver.eigenvalues()).array()>=0).all();
}

bool CXSM::CheckGlobalMinimum()
{
    FindLocalMinima();
    if (IndexInput >= 0)
    {
        if (!LocalMinimaQ[IndexInput])
        {
            // ! The input vacuum is not local minimum
            return false;
        }
        double VtotEW = Vtotal[IndexInput];
        bool goodEW = true;
        for (size_t i = 0; i < NLocalMinima; i++)
        {
            goodEW *= (VtotEW <= (Vtotal[MinimaIndex[i]]+1e-3));
        }
        return goodEW;
    }
    else
    {
        cout<<"Error in finding the input vacuum"<<endl;
        return false;
    } 
}


Toy::Toy(double lambda, double eta)
{
    _lambda = lambda;
    _eta = eta;
}

double Toy::Vtot(double phi)
{
    return _lambda/4.0*pow(phi*phi-_eta*_eta,2);
}
double Toy::dVdphi(double phi)
{
    return _lambda*phi*(phi*phi-_eta*_eta);
}
double Toy::Get_V0_global()
{
    return Vtot(_eta);
}

double CXSM::GetTotalEnergy(VD x, VVD y)
{
    double energy = 0;
    double V0 = Vtot(vev,VS,0);
    for (size_t i = 0; i < x.size()-1; i++)
    {
        double density = 0;
        double DeltaZ = x[i+1]-x[i];
        VD yaver = (y[i]+y[i+1])/2;
        density += pow(vev,4)*(y[i+1]-y[i])*(y[i+1]-y[i])/pow(DeltaZ,2)/2;
        density += Vtot(yaver[0]*vev,yaver[1]*vev,0) - V0;
        energy += density*DeltaZ/vev;
    }
    return energy;
}
double CXSM::GetTension(VD x, VVD y)
{
    double tension = 0;
    for (size_t i = 0; i < x.size()-1; i++)
    {
        double _DeltaZ = x[i+1]-x[i];
        VD dfield_dZ = (y[i+1]-y[i])/_DeltaZ;
        tension += dfield_dZ*dfield_dZ * _DeltaZ;
    }
    return tension*pow(vev,3);
}