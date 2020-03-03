/*
 * @Description  : 
 * @Author       : Yongcheng Wu
 * @Date         : 2020-02-08 15:48:52
 * @LastEditors  : Yongcheng Wu
 * @LastEditTime : 2020-02-09 23:44:59
 */
#include "SM_cxSM.h"
#include "Relaxation.h"
#include <iostream>
#include <cmath>

using namespace std;

struct cxSMCP_Param
{
    double mu2;
    double b1;
    double b2;
    double lambda;
    double d1;
    double d2;
    double delta2;
    double vs;
    double vev;
};


void DIFEQ_ODD(const Relaxation_Param relax_param, void *param, VVD &S)
{
    for (size_t i = 0; i < S.size(); i++)
    {
        for (size_t j = 0; j < S[i].size(); j++)
        {
            S[i][j] = 0;
        }
    }
    
    int k = relax_param.k;
    // cout<<"Setting k: "<<k<<endl;
    double x1 = relax_param.x1;
    double x2 = relax_param.x2;
    VD y1 = relax_param.y1;
    VD y2 = relax_param.y2;
    double dr = x2-x1;
    double raver = (x1+x2)/2.0;
    VD yaver = (y1+y2)/2.0;
    VD dy = y2 - y1;
    cxSMCP_Param *spa = (cxSMCP_Param*)param;
    double vev = spa->vev;
    double vs = spa->vs;
    double mu2hat = spa->mu2/vev/vev;
    double b1hat = spa->b1/vev/vev;
    double b2hat = spa->b2/vev/vev;
    double lambda = spa->lambda;
    double d1 = spa->d1;
    double d2 = spa->d2;
    double delta2 = spa->delta2;
    if (k==relax_param.k_init)
    {
        S[3][6] = 1.0;
        S[3][7] = 0.0;
        S[3][8] = 0.0;
        S[3][9] = 0.0;
        S[3][10] = 0.0;
        S[3][11] = 0.0;

        S[4][6] = 0.0;
        S[4][7] = 1.0;
        S[4][8] = 0.0;
        S[4][9] = 0.0;
        S[4][10] = 0.0;
        S[4][11] = 0.0;

        S[5][6] = 0.0;
        S[5][7] = 0.0;
        S[5][8] = 1.0;
        S[5][9] = 0.0;
        S[5][10] = 0.0;
        S[5][11] = 0.0;

    // at -inf:
    // hhat = 1;
    // shat = -vs/v;
    // alpha = 1/2 acos(b1/-d1/vs/vs);

        S[3][relax_param.k_coeff] = y1[0] - 1;
        S[4][relax_param.k_coeff] = y1[1] + vs/vev;
        S[5][relax_param.k_coeff] = y1[2] - acos(-b1hat*vev*vev/d1/vs/vs)/2;

    // at -inf:
    // hhat = 1;
    // shat = vs/v;
    // alpha = -1/2 acos(b1/-d1/vs/vs);
        // S[3][relax_param.k_coeff] = y1[0] - 1;
        // S[4][relax_param.k_coeff] = y1[1] - vs/vev;
        // S[5][relax_param.k_coeff] = y1[2] + acos(-b1hat*vev*vev/d1/vs/vs)/2;
        return;
    }
    if (k > relax_param.k_final)
    {
        S[0][6] = 1.0;
        S[0][7] = 0.0;
        S[0][8] = 0.0;
        S[0][9] = 0.0;
        S[0][10] = 0.0;
        S[0][11] = 0.0;

        S[1][6] = 0.0;
        S[1][7] = 1.0;
        S[1][8] = 0.0;
        S[1][9] = 0.0;
        S[1][10] = 0.0;
        S[1][11] = 0.0;

        S[2][6] = 0.0;
        S[2][7] = 0.0;
        S[2][8] = 1.0;
        S[2][9] = 0.0;
        S[2][10] = 0.0;
        S[2][11] = 0.0;

    // at +inf:
    // hhat = 1;
    // shat = vs/v;
    // alpha = 1/2 acos(b1/-d1/vs/vs);

        S[0][relax_param.k_coeff] = y1[0] - 1;
        S[1][relax_param.k_coeff] = y1[1] - vs/vev;
        S[2][relax_param.k_coeff] = y1[2] - acos(-b1hat*vev*vev/d1/vs/vs)/2;

        return;
    }

    S[0][0] = -dr*mu2hat/2 - 3*dr*lambda*yaver[0]*yaver[0]/2 - dr*delta2*yaver[1]*yaver[1]/8;
    S[0][1] = -dr*delta2*yaver[0]*yaver[1]/4;
    S[0][2] = 0;
    S[0][3] = -1;
    S[0][4] = 0;
    S[0][5] = 0;
    S[0][6] = -dr*mu2hat/2 - 3*dr*lambda*yaver[0]*yaver[0]/2 - dr*delta2*yaver[1]*yaver[1]/8;
    S[0][7] = -dr*delta2*yaver[0]*yaver[1]/4;
    S[0][8] = 0;
    S[0][9] = 1;
    S[0][10] = 0;
    S[0][11] = 0;
    
    S[1][0] = -dr*delta2*yaver[0]*yaver[1]/4;
    S[1][1] = -dr*yaver[5]*yaver[5]/2 - dr*delta2*yaver[0]*yaver[0]/8-dr*b2hat/4 - 3*dr*d2*yaver[1]*yaver[1]/8 - dr*b1hat*cos(2*yaver[2])/4 - 3*dr*d1*yaver[1]*yaver[1]*cos(4*yaver[2])/8;
    S[1][2] = dr*b1hat*yaver[1]*sin(2*yaver[2])/2 + dr*d1*pow(yaver[1],3)*sin(4*yaver[2])/2;
    S[1][3] = 0;
    S[1][4] = -1;
    S[1][5] = -dr*yaver[1]*yaver[5];
    S[1][6] = -dr*delta2*yaver[0]*yaver[1]/4;
    S[1][7] = -dr*yaver[5]*yaver[5]/2 - dr*delta2*yaver[0]*yaver[0]/8-dr*b2hat/4 - 3*dr*d2*yaver[1]*yaver[1]/8 - dr*b1hat*cos(2*yaver[2])/4 - 3*dr*d1*yaver[1]*yaver[1]*cos(4*yaver[2])/8;
    S[1][8] = dr*b1hat*yaver[1]*sin(2*yaver[2])/2 + dr*d1*pow(yaver[1],3)*sin(4*yaver[2])/2;
    S[1][9] = 0;
    S[1][10] = 1;
    S[1][11] = -dr*yaver[1]*yaver[5];

    S[2][0] = 0;
    S[2][1] = dy[5]*yaver[1] + dr*yaver[4]*yaver[5] + dr*b1hat*yaver[1]*sin(2*yaver[2])/2 + dr*d1*pow(yaver[1],3)*sin(4*yaver[2])/2;
    S[2][2] = dr*b1hat*pow(yaver[1],2)*cos(2*yaver[2])/2 + dr*d1*pow(yaver[1],4)*cos(4*yaver[4])/2;
    S[2][3] = 0;
    S[2][4] = dr*yaver[1]*yaver[5];
    S[2][5] = -yaver[1]*yaver[1] + dr*yaver[1]*yaver[4];
    S[2][6] = 0;
    S[2][7] = dy[5]*yaver[1] + dr*yaver[4]*yaver[5] + dr*b1hat*yaver[1]*sin(2*yaver[2])/2 + dr*d1*pow(yaver[1],3)*sin(4*yaver[2])/2;
    S[2][8] = dr*b1hat*pow(yaver[1],2)*cos(2*yaver[2])/2 + dr*d1*pow(yaver[1],4)*cos(4*yaver[4])/2;
    S[2][9] = 0;
    S[2][10] = dr*yaver[1]*yaver[5];
    S[2][11] = yaver[1]*yaver[1] + dr*yaver[1]*yaver[4];

    S[3][0] = -1;
    S[3][1] = 0;
    S[3][2] = 0;
    S[3][3] = -dr/2;
    S[3][4] = 0;
    S[3][5] = 0;
    S[3][6] = 1;
    S[3][7] = 0;
    S[3][8] = 0;
    S[3][9] = -dr/2;
    S[3][10] = 0;
    S[3][11] = 0;

    S[4][0] = 0;
    S[4][1] = -1;
    S[4][2] = 0;
    S[4][3] = 0;
    S[4][4] = -dr/2;
    S[4][5] = 0;
    S[4][6] = 0;
    S[4][7] = 1;
    S[4][8] = 0;
    S[4][9] = 0;
    S[4][10] = -dr/2;
    S[4][11] = 0;

    S[5][0] = 0;
    S[5][1] = 0;
    S[5][2] = -1;
    S[5][3] = 0;
    S[5][4] = 0;
    S[5][5] = -dr/2;
    S[5][6] = 0;
    S[5][7] = 0;
    S[5][8] = 1;
    S[5][9] = 0;
    S[5][10] = 0;
    S[5][11] = -dr/2;

    S[0][relax_param.k_coeff] = dy[3] - dr*(mu2hat*yaver[0] + lambda*pow(yaver[0],3) + delta2*yaver[0]*yaver[1]*yaver[1]/4);
    S[1][relax_param.k_coeff] = dy[4] - dr*yaver[1]*yaver[5]*yaver[5] - dr*(delta2*yaver[0]*yaver[0]*yaver[1]/4+b2hat*yaver[1]/2+d2*pow(yaver[1],3)/4 + b1hat*yaver[1]*cos(2*yaver[2])/2 + d1*pow(yaver[1],3)*cos(4*yaver[2])/4);
    S[2][relax_param.k_coeff] = yaver[1]*yaver[1]*dy[5] + 2*dr*yaver[1]*yaver[4]*yaver[5] + dr*(b1hat*yaver[1]*yaver[1]*sin(2*yaver[2])/2 + d1*pow(yaver[1],4)*sin(4*yaver[2])/4);
    S[3][relax_param.k_coeff] = dy[0] - dr*yaver[3];
    S[4][relax_param.k_coeff] = dy[1] - dr*yaver[4];
    S[5][relax_param.k_coeff] = dy[2] - dr*yaver[5];
    
}

double GetTension(VD x, VVD y)
{
    double res = 0;
    for (int i = 0; i < x.size()-1; i++)
    {
        double dx = x[i+1]-x[i];
        VD dy = y[i+1] - y[i];
        VD yaver = (y[i+1]+y[i])/2;
        double dhdx = dy[0]/dx;
        double dsdx = dy[1]/dx;
        double dadx = dy[2]/dx;

        res += dx*(dhdx*dhdx + dsdx*dsdx + yaver[1]*yaver[1]*dadx*dadx)/2;
    }
    return res;
}

int main(int argc, char const *argv[])
{
    double mu2 = -0.8*246*246;
    double b1 = -0.9*50*50;
    double b2 = -0.75*50*50;
    double lambda = 1.4;
    double d1 = 1.1;
    double d2 = 0.6;
    double delta2 = 0.6;
    double vs = 50;
    double vev = 246;
    cxSMCP_Param spa = {mu2,b1,b2,lambda,d1,d2,delta2,vs,vev};
    double alp = acos(-b1/d1/vs/vs)/2;
    Relaxation RX(6,3,501);
    VD y_beg = {1,-vs/vev,alp,0,0,0};
    VD y_end = {1,vs/vev,alp,0,0,0};
    VD scales(6);
    scales[0] = 1.0;
    scales[1] = 1.0;
    scales[2] = 1.0;
    scales[3] = 1.0;
    scales[4] = 1.0;
    scales[5] = 1.0;
    RX.SetBoundary(-20,20,y_beg,y_end);
    RX.SetMaxIteration(10000);
    RX.SetConvergeCriterion(1.0,1e-10);
    RX.SetODESystem(DIFEQ_ODD,&spa);
    RX.SetScales(scales);
    // RX._INIT();
    // RX.PrintSolution();
    RX.SOLVDE();
    VD X = RX.GetX();
    VVD Y = RX.GetY();
    cout<<"Tension: "<<GetTension(X,Y)<<endl;
    RX.DumpSolution("Relaxation_CPDW_test_1.dat");
    return 0;
}

