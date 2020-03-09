#include "Relaxation.h"
#include <iostream>
#include <cmath>

using namespace std;

struct Vortex_Param
{
    int n2;
    double lambda;
    double e2;
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
    double sth = sin(raver);
    double cth = cos(raver);
    double sth2 = sth*sth;
    double cth2 = cth*cth;
    double cth4 = cth2*cth2;
    double th = sth/cth;
    VD yaver = (y1+y2)/2.0;
    VD dy = y2 - y1;
    Vortex_Param *spa = (Vortex_Param*)param;
    int n2 = spa->n2;
    double lambda = spa->lambda;
    double e2 = spa->e2;
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

        S[2][relax_param.k_coeff] = y1[0];
        S[3][relax_param.k_coeff] = y1[1];
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
        S[1][relax_param.k_coeff] = y1[1]-1;

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

    S[2][0] = -dr*n2*pow(yaver[1],2)/cth2/sth2/2.0 - dr*lambda/cth4*(3*pow(yaver[0],2)-1)/4.0;
    S[2][1] = -n2*dr/sth2/cth2*yaver[0]*(yaver[1]-1);
    S[2][2] = -1 + dr*(1/sth/cth-2*th)/2.0;
    S[2][3] = 0.0;
    S[2][4] = S[2][0];
    S[2][5] = S[2][1];
    S[2][6] = 1 + dr*(1/sth/cth-2*th)/2.0;
    S[2][7] = 0.0;

    S[3][0] = -dr*2*e2*yaver[0]*(yaver[1]-1)/cth4;
    S[3][1] = -dr*2*e2*pow(yaver[0],2)/cth4/2.0;
    S[3][2] = 0.0;
    S[3][3] = -1 - dr*(1/sth/cth+2*th)/2.0;
    S[3][4] = S[3][0];
    S[3][5] = S[3][1];
    S[3][6] = 0.0;
    S[3][7] = 1 - dr*(1/sth/cth+2*th)/2.0;

    S[0][relax_param.k_coeff] = dy[0] - dr*yaver[2];
    S[1][relax_param.k_coeff] = dy[1] - dr*yaver[3];
    S[2][relax_param.k_coeff] = dy[2] + dr*(1/sth/cth-2*th)*yaver[2]-dr*n2*yaver[0]*pow(yaver[1]-1,2)/sth2/cth2 - dr*lambda*yaver[0]*(pow(yaver[0],2)-1)/2/cth4;
    S[3][relax_param.k_coeff] = dy[3] - dr*(1/sth/cth+2*th)*yaver[3] - dr*2*e2*pow(yaver[0],2)*(yaver[1]-1)/cth4;
}

int main(int argc, char const *argv[])
{
    double PI = 3.14159265;
    Vortex_Param spa = {1,1,0.091725};
    Relaxation RX(4,2,101);
    VD y_beg = {0,0,0.0,0.0};
    VD y_end = {1,1,0.0,0.0};
    VD scales(4);
    scales[0] = 1.0;
    scales[1] = 1.0;
    scales[2] = 1.0;
    scales[3] = 1.0;
    RX.SetBoundary(0,PI/2.0,y_beg,y_end);
    RX.SetMaxIteration(5000000);
    RX.SetODESystem(DIFEQ_ODD,&spa);
    RX.SetScales(scales);
    // RX._INIT();
    // RX.PrintSolution();
    RX.SOLVDE();
    RX.PrintSolution();
    RX.DumpSolution("Relaxation_vortex_theta_test.dat");
    return 0;
}

