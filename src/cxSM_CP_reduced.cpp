#include "cxSM_CP_reduced.h"
// #include "QuarticSolver.h"
#include <functional>
#include <iostream>

using namespace std;
using namespace Eigen;

cxSM_CP_reduced::cxSM_CP_reduced():cxSM_CP()
{

}

void cxSM_CP_reduced::Set_Potential_Parameters(double mu2, double lam, double del2, double b2, double d2, double b1, double d1)
{
    cxSM_CP::Set_Potential_Parameters(mu2,lam,del2,b2,d2,0.0,b1,d1,0.0);
}

bool cxSM_CP_reduced::Set_Physical_Parameters_vsr_vsi_theta(double vsr, double vsi, double MHH, double MHA, double theta1)
{
    double theta2,theta3;
    bool good = _GetTheta2Theta3(vsr,vsi,MHH,MHA,theta1,theta2,theta3);
    if (!good)
    {
        return false;
    }
    Set_Physical_Parameters(vsr,vsi,MHH,MHA,theta1,theta2,theta3);
    return true;
}

bool cxSM_CP_reduced::_GetTheta2Theta3(const double vsr, const double vsi, const double MHH, const double MHA, const double theta1, double &theta2, double &theta3)
{
    double k = vsr/vsi;
    double k2 = k*k;
    double ta1 = tan(theta1);
    double dm132 = MHL2 - MHA*MHA;
    double dm232 = MHH*MHH - MHA*MHA;
    double m32 = MHA*MHA;
    double C1 = 1/2.0/ta1;
    double C2 = (1-k2)*m32/dm232/ta1/2.0;
    double C3 = dm232*ta1/dm132/2.0;
    double D11 = -((1 + pow(k,2))*(-1 + 4*pow(C3,2)*pow(ta1,2)));
    double D12 = 2*(-1 + pow(k,2))*(1 + 4*C1*C3*pow(ta1,2));
    double D13 = -((1 + pow(k,2))*(-1 + 4*pow(C1,2)*pow(ta1,2)));
    double D14 = 8*C2*C3*(-1 + pow(k,2))*pow(ta1,2) + 4*pow(k,2)*(-1 + pow(ta1,2));
    double D15 = -8*C1*C2*(1 + pow(k,2))*pow(ta1,2);
    double D16 = -4*pow(C2,2)*(1 + pow(k,2))*pow(ta1,2);
    double D21 = -(C3*(1 + pow(k,2))*(1 + 2*C3*ta1));
    double D22 = (-1 + pow(k,2))*(C1 - C3 + 4*C1*C3*ta1);
    double D23 = -(C1*(1 + pow(k,2))*(-1 + 2*C1*ta1));
    double D24 = 2*pow(k,2)*ta1 + C2*(-1 + pow(k,2))*(1 + 4*C3*ta1);
    double D25 = -(C2*(1 + pow(k,2))*(-1 + 4*C1*ta1));
    double D26 = -2*pow(C2,2)*(1 + pow(k,2))*ta1;


    // * With further simplification K4 = 0, K3 = 0
    // double K4 = pow(D13,2)*pow(D21,2) + D23*(pow(D12,2)*D21 - D11*D12*D22 + pow(D11,2)*D23) + D13*(-(D12*D21*D22) + D11*(pow(D22,2) - 2*D21*D23));
    // double K3 = pow(D12,2)*D21*D25 + D13*(2*D15*pow(D21,2) - D14*D21*D22 - D12*D21*D24 + 2*D11*D22*D24 - 2*D11*D21*D25) - D12*(D15*D21*D22 - 2*D14*D21*D23 + D11*D23*D24 + D11*D22*D25) + D11*(D15*pow(D22,2) - 2*D15*D21*D23 - D14*D22*D23 + 2*D11*D23*D25);
    double K2 = pow(D15,2)*pow(D21,2) - D12*D16*D21*D22 + D11*D16*pow(D22,2) + pow(D14,2)*D21*D23 - 2*D11*D16*D21*D23 - D11*D14*D23*D24 + 2*D12*D14*D21*D25 - D11*D14*D22*D25 - D11*D12*D24*D25 + pow(D11,2)*pow(D25,2) - D15*(D14*D21*D22 + D12*D21*D24 - 2*D11*D22*D24 + 2*D11*D21*D25) + pow(D12,2)*D21*D26 - D11*D12*D22*D26 + 2*pow(D11,2)*D23*D26 + D13*(2*D16*pow(D21,2) - D14*D21*D24 + D11*pow(D24,2) - 2*D11*D21*D26);
    double K1 = -(D12*D16*D21*D24) + 2*D11*D16*D22*D24 + pow(D14,2)*D21*D25 - 2*D11*D16*D21*D25 - D11*D12*D24*D26 + 2*pow(D11,2)*D25*D26 + D15*(2*D16*pow(D21,2) - D14*D21*D24 + D11*pow(D24,2) - 2*D11*D21*D26) - D14*(D16*D21*D22 + D11*D24*D25 - 2*D12*D21*D26 + D11*D22*D26);
    double K0 = pow(D16,2)*pow(D21,2) + D26*(pow(D14,2)*D21 - D11*D14*D24 + pow(D11,2)*D26) + D16*(-(D14*D21*D24) + D11*(pow(D24,2) - 2*D21*D26));

    double Delta = K1*K1-4*K2*K0;
    if (Delta<0)
    {
        return false;
    }
    VD y_sols = {(-K1+sqrt(Delta))/2.0/K2,(-K1-sqrt(Delta))/2.0/K2};
    // cout<<"K2\tK1\tK0"<<endl;
    // cout<<K2<<"\t"<<K1<<"\t"<<K0<<endl;
    // QuarticSolver solver(K4,K3,K2,K1,K0);
    // solver.Solve();
    // VD y_sols = solver.GetRealSolution();
    // for (int i = 0; i < 4; i++)
    // {
        // cout<<solver.SOLUTIONS[i]<<endl;
    // }
    
    // if (y_sols.size()==0)
    // {
    //     return false;
    // }
    double y;
    double x1,x2;
    double a = D11;
    double b,c;
    double delta;
    bool good = false;
    double Rm,Rp;
    double Rm2,Rp2;
    function<bool(double &, double &)> getmixing = [&](double &th2, double &th3){ 
        double R22 = (Rm+Rp)/2.0;
        if (abs(R22)>1)
        {
            return false;
        }
        double R32 = (Rp-Rm)/k/2.0;
        if (abs(R32)>1)
        {
            return false;
        }
        double R21 = C1*Rp+C2/Rm-C3*Rm;
        if (abs(R21)>1)
        {
            return false;
        }
        double R31 = (C1*Rp+C2/Rm+C3*Rm)/k;
        if (abs(R31)>1)
        {
            return false;
        }
        if (R21*R21+R22*R22>1)
        {
            return false;
        }
        if (R31*R31+R32*R32>1)
        {
            return false;
        }
        if (R21*R21+R31*R31>1)
        {
            return false;
        }
        if (R22*R22+R32*R32>1)
        {
            return false;
        }
        double cth1 = cos(theta1);
        double sth1 = sin(theta1);
        double cth3 = cth1*R22-sth1*R21;
        double sth3 = sth1*R31-cth1*R32;
        if (abs(cth3)>1 || abs(sth3)>1)
        {
            return false;
        }
        double sth2 = (cth1*cth3-R22)/sth1/sth3;
        double cth2 = sqrt(1-sth2*sth2);// negative value is equivalent
        th2 = atan2(sth2,cth2);
        th3 = atan2(sth3,cth3);
        return true;
    };
    for (int i = 0; i < y_sols.size(); i++)
    {
        y = y_sols[i];
        // cout<<"y"<<i<<" "<<y<<endl;
        if (abs(y) > 1.0)
        {
            continue;
        }
        b = D12*y+D14;
        c = D13*y*y+D15*y+D16;
        delta = b*b-4*a*c;
        if (delta<0)
        {
            continue;
        }
        x1 = (-b+sqrt(delta))/a/2.0;
        x2 = (-b-sqrt(delta))/a/2.0;
        if ((x1<0||x1>1) && (x2<0||x2>1))
        {
            continue;
        }
        if (x1 > 0)
        {
            Rm = sqrt(x1);// negative value is equivalent
            Rp = y/Rm;
            bool test = getmixing(theta2,theta3);
            if (test)
            {
                good = true;
                break;
            }
        } 
        if (x2 > 0)
        {
            Rm = sqrt(x2);
            Rp = y/Rm;
            bool test = getmixing(theta2,theta3);
            if (test)
            {
                good = true;
                break;
            }
        }        
    }
    return good;
    
}

bool cxSM_CP_reduced::Set_Physical_Parameters_vs_theta(double vs, double MHH, double MHA, double theta1, double theta3)
{
    double alpha,theta2;
    bool good = _GetAlphaTheta2(vs,MHH,MHA,theta1,theta3,alpha,theta2);
    if (!good)
    {
        return false;
    }
    double vsr = vs*cos(alpha);
    double vsi = vs*sin(alpha);
    Set_Physical_Parameters(vsr,vsi,MHH,MHA,theta1,theta2,theta3);
    return true;
}
bool cxSM_CP_reduced::_GetAlphaTheta2(const double vs, const double MHH, const double MHA, const double theta1, double theta3, double &alpha, double &theta2)
{
    // In this situation, we always have solution
    double m12 = MHL2;
    double m22 = MHH*MHH;
    double m32 = MHA*MHA;

    // * First solution
    theta2 = M_PI_2; // * Only consider positive one, negative one should be equivalent (! not varified)
    alpha = atan((m12+m22-(m22-m12)*cos(2*(theta1+theta3)))/(m12+m22+(m22-m12)*cos(2*(theta1+theta3))));

    double K2 = 4*(m12 + m22 - 2*m32 + (m12 - m22)*cos(2*theta1))*(2*m12*m22 - (m12 + m22)*m32 + (m12 - m22)*m32*cos(2*theta1))*cos(2*theta3);
    double K1 = -8*(m12 - m22)*(m12*m22 - pow(m32,2) + (m12 - m22)*m32*cos(2*theta1))*sin(2*theta1)*sin(2*theta3);
    double K0 = -4*pow(m12 - m22,2)*m32*cos(2*theta3)*pow(sin(2*theta1),2);

    double delta = K1*K1-4*K2*K0;
    if (delta < 0)
    {
        return true; // Always good, since we have theta2= pi/2 solution
    }
    double x1 = (-K1+sqrt(delta))/2/K2;
    double x2 = (-K1-sqrt(delta))/2/K2;
    double sth2_tmp = abs(x1)<abs(x2)?x1:x2;
    if (abs(sth2_tmp) > 1)
    {
        return true; // We have pi/2 solution
    }
    
    theta2 = asin(sth2_tmp); // * assume cth2 > 0, cth2 < 0 should be equivalent (!not varified)

    double ta = ((-m12 + m22)*cos(theta1)*cos(theta3)*sin(theta1) + sth2_tmp*(m32 - m12*pow(cos(theta1),2) - m22*pow(sin(theta1),2))*sin(theta3))/(sth2_tmp*cos(theta3)*(m32 - m12*pow(cos(theta1),2) - m22*pow(sin(theta1),2)) + (m12 - m22)*cos(theta1)*sin(theta1)*sin(theta3));
    
    alpha = atan(1/ta);
   
    return true;
}
