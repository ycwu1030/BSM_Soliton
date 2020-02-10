/*
 * @Description  : 
 * @Author       : Yongcheng Wu
 * @Date         : 2020-02-08 15:48:52
 * @LastEditors  : Yongcheng Wu
 * @LastEditTime : 2020-02-09 21:25:33
 */
#include "SM_cxSM.h"
#include "Relaxation.h"
#include <iostream>
#include <cmath>

using namespace std;

struct cxSM_Param
{
    double lambda;
    double delta2;
    double d2;
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
    cxSM_Param *spa = (cxSM_Param*)param;
    double lambda = spa->lambda;
    double delta2 = spa->delta2;
    double d2 = spa->d2;
    double vs = spa->vs;
    double vev = spa->vev;
    if (k==relax_param.k_init)
    {
        S[2][4] = 1.0;
        S[2][5] = 0.0;
        S[2][6] = 0.0;
        S[2][7] = 0.0;

        S[3][4] = 0.0;
        S[3][5] = 1.0;
        S[3][6] = 0.0;
        S[3][7] = 0.0;

        S[2][relax_param.k_coeff] = y1[0] - 1;
        S[3][relax_param.k_coeff] = y1[1] + vs/vev;
        return;
    }
    if (k > relax_param.k_final)
    {
        S[0][4] = 1.0;
        S[0][5] = 0.0;
        S[0][6] = 0.0;
        S[0][7] = 0.0;

        S[1][4] = 0.0;
        S[1][5] = 1.0;
        S[1][6] = 0.0;
        S[1][7] = 0.0;

        S[0][relax_param.k_coeff] = y1[0]-1;
        S[1][relax_param.k_coeff] = y1[1]-vs/vev;

        return;
    }

    S[0][0] = -1.0;
    S[0][1] = 0.0;
    S[0][2] = -dr/2;
    S[0][3] = 0.0;
    S[0][4] = 1.0;
    S[0][5] = 0.0;
    S[0][6] = -dr/2;
    S[0][7] = 0.0;

    S[1][0] = 0.0;
    S[1][1] = -1.0;
    S[1][2] = 0.0;
    S[1][3] = -dr/2;
    S[1][4] = 0.0;
    S[1][5] = 1.0;
    S[1][6] = 0.0;
    S[1][7] = -dr/2;

    S[2][0] = dr*(lambda+delta2/4*vs*vs/vev/vev)/2.0-3*lambda*pow(yaver[0],2)*dr/2.0 - delta2/4*pow(yaver[1],2)*dr/2;
    S[2][1] = -dr*delta2/4.0*yaver[0]*yaver[1];
    S[2][2] = -1;
    S[2][3] = 0.0;
    S[2][4] = S[2][0];
    S[2][5] = S[2][1];
    S[2][6] = 1;
    S[2][7] = 0.0;

    S[3][0] = -delta2/2.0*yaver[0]*yaver[1]*dr/2;
    S[3][1] = (delta2+d2*vs*vs/vev/vev)*dr/8-delta2/4*pow(yaver[0],2)*dr/2-3*d2/4*pow(yaver[1],2)*dr/2;
    S[3][2] = 0.0;
    S[3][3] = -1;
    S[3][4] = S[3][0];
    S[3][5] = S[3][1];
    S[3][6] = 0.0;
    S[3][7] = 1;

    S[0][relax_param.k_coeff] = dy[0] - dr*yaver[2];
    S[1][relax_param.k_coeff] = dy[1] - dr*yaver[3];
    S[2][relax_param.k_coeff] = dy[2] + dr*(lambda+delta2/4*vs*vs/vev/vev)*yaver[0] - dr*lambda*pow(yaver[0],3) - dr*delta2/4*yaver[0]*pow(yaver[1],2);
    S[3][relax_param.k_coeff] = dy[3] +dr/4*(delta2+d2*vs*vs/vev/vev)*yaver[1] - dr*delta2/4*pow(yaver[0],2)*yaver[1] - dr*d2/4.0*pow(yaver[1],3);
}

int main(int argc, char const *argv[])
{
    CXSM model;
    model.SetInput(157.633,42.647,0.0,0.335465);
    cxSM_Param spa = {model.lam,model.del2,model.d2,model.VS,model.vev};
    Relaxation RX(4,2,101);
    VD y_beg = {1,-1,0.0,0.0};
    VD y_end = {1,1,0.0,0.0};
    VD scales(4);
    scales[0] = 1.0;
    scales[1] = 1.0;
    scales[2] = 1.0;
    scales[3] = 1.0;
    RX.SetBoundary(-100,100,y_beg,y_end);
    RX.SetMaxIteration(10000);
    RX.SetConvergeCriterion(1.0,1e-10);
    RX.SetODESystem(DIFEQ_ODD,&spa);
    RX.SetScales(scales);
    // RX._INIT();
    // RX.PrintSolution();
    RX._SOLVDE();
    RX.PrintSolution();
    RX.DumpSolution("Relaxation_DW_test.dat");
    return 0;
}

