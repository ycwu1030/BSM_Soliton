#include "Relaxation.h"
#include <iostream>
#include <cmath>

using namespace std;
/*
 * Spheroidal Harmonics example for Relaxation
 * The Original ODE is:
 * d/dx((1-x^2)dS/dx) + (\lambda - c^2 x^2 - m^2/(1-x^2))S = 0
 * After remove the singular point and some other simplifications, 
 * we got following ODEs (More detail, please check 'Numerical Recipes')
 * DOF = 3
 * dy1/dx = y2
 * dy2/dx = 1/(1-x^2)(2x(m+1)y2-(y3-c^2x^2)y1)
 * dy3/dx = 0
 * The solution range is chose to be [0,1].
 * The boundary conditions are considered separately according to n-m
 * If n-m odd:
 *      At x=0:
 *          y1 = 0
 *      At x=1:
 *          y2 = (y3-c^2)y1/2/(m+1)
 *          y1 = (-1)^m(n+m)!/2^m/m!/(n-m)! 
*/
struct Spheroidal_Param
{
    int m;
    int n;
    double c2;
    double gamma;
};

double plgndr(int l, int m, double x)
{
    double fact, pll, pmm, pmmp1, somx2;

    pmm = 1.0;
    if (m>0)
    {
        somx2 = sqrt((1-x)*(1+x));
        fact = 1.0;
        for (size_t i = 0; i < m; i++)
        {
            pmm *= -fact*somx2;
            fact += 2.0;
        }
    }
    if (l == m)
    {
        return pmm;
    }
    else
    {
        pmmp1 = x*(2*m+1)*pmm;
        if (l==(m+1))
        {
            return pmmp1;
        }
        else
        {
            for (size_t ll = m+2; ll <= l; ll++)
            {
                pll = (x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m);
                pmm = pmmp1;
                pmmp1=pll;
            }
            return pll;
        }   
    }
}

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
    Spheroidal_Param *spa = (Spheroidal_Param*)param;
    if (k==relax_param.k_init)
    {
        S[2][3] = 1.0;
        S[2][4] = 0.0;
        S[2][5] = 0.0;
        S[2][relax_param.k_coeff] = y1[0];
        return;
    }
    if (k > relax_param.k_final)
    {
        S[0][3] = -(y1[2]-spa->c2)/(2*(spa->m+1));
        S[0][4] = 1.0;
        S[0][5] = -y1[0]/(2*(spa->m+1));
        S[0][relax_param.k_coeff] = y1[1]-(y1[2]-spa->c2)*y1[0]/(2*(spa->m+1));

        S[1][3] = 1.0;
        S[1][4] = 0.0;
        S[1][5] = 0.0;
        S[1][relax_param.k_coeff] = y1[0] - spa->gamma;

        return;
    }

    double h = x2-x1;
    S[0][0] = -1.0;
    S[0][1] = -h/2;
    S[0][2] = 0.0;
    S[0][3] = 1.0;
    S[0][4] = -h/2;
    S[0][5] = 0.0;

    double temp1 = x1 + x2;
    double temp = h/(1.0-temp1*temp1/4.0);
    double temp2 = (y1[2]+y2[2])/2.0 - spa->c2*temp1*temp1/4;
    S[1][0] = temp*temp2/2;
    S[1][1] = -1.0 - temp*(spa->m+1)*temp1/2;
    S[1][2] = temp*(y1[0]+y2[0])/4;
    S[1][3] = S[1][0];
    S[1][4] = 2.0 + S[1][1];
    S[1][5] = S[1][2];

    S[2][0] = 0.0;
    S[2][1] = 0.0;
    S[2][2] = -1.0;
    S[2][3] = 0.0;
    S[2][4] = 0.0;
    S[2][5] = 1.0;

    S[0][relax_param.k_coeff] = y2[0]-y1[0] - h*(y2[1]+y1[1])/2;
    S[1][relax_param.k_coeff] = y2[1]-y1[1] - temp*((x1+x2)*(spa->m+1)*(y1[1]+y2[1])/2-temp2*(y1[0]+y2[0])/2);
    S[2][relax_param.k_coeff] = y2[2] - y1[2];
}
double fac(int n)
{
    double res=1.0;
    for (size_t i = 1; i <= n; i++)
    {
        res*=i;
    }
    return res;
}

int main(int argc, char const *argv[])
{
    int m = 2;
    int n = 5;
    double c2 = 1.0;
    double gamma = pow(-1,m)*fac(n+m)/pow(2,m)/fac(m)/fac(n-m);
    Spheroidal_Param spa = {m,n,c2,gamma};
    Relaxation RX(3,1,41);
    VD y_beg = {0.1,0.1,20.0};
    VD y_end = {105,2.0/3.0*105,20.0};
    VD scales(3);
    scales[0] = abs(gamma);
    scales[1] = y_end[1]>abs(gamma)?y_end[1]:abs(gamma);
    scales[2] = y_end[2]>1?y_end[2]:1;
    VD X_Guess;
    VVD Y_Guess;
    VD _Y(3);
    double h = 1.0/40.0;
    double fac1,fac2,deriv;
    for (size_t k = 0; k < 40; k++)
    {
        X_Guess.push_back(k*h);
        fac1 = 1.0 - k*h*k*h;
        fac2 = exp((-m/2.0)*log(fac1));
        _Y[0] = plgndr(n,m,k*h)*fac2;
        deriv = -((n-m+1)*plgndr(n+1,m,k*h)-(n+1)*k*h*plgndr(n,m,k*h))/fac1;
        _Y[1] = m*k*h*_Y[0]/fac1 + deriv*fac2;
        _Y[2] = n*(n+1) - m*(m+1);
        Y_Guess.push_back(_Y);
    }
    X_Guess.push_back(1);
    _Y[0] = gamma;
    _Y[2] = n*(n+1) - m*(m+1);
    _Y[1] = (_Y[2]-c2)*_Y[0]/(2.0*(m+1));
    Y_Guess.push_back(_Y);
    RX.SetBoundary(X_Guess,Y_Guess); 
    // RX.SetBoundary(0,1,y_beg,y_end);
    RX.SetMaxIteration(100);
    RX.SetODESystem(DIFEQ_ODD,&spa);
    RX.SetScales(scales);
    // RX._INIT();
    // RX.PrintSolution();
    RX.SOLVDE();
    RX.PrintSolution();
    RX.DumpSolution("Relaxation_test.dat");
    return 0;
}

