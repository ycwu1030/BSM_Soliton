#include "QuarticSolver.h"
#include <cmath>

using namespace std;

// cdouble pow(const cdouble &x,double n)
// {
//     double r = abs(x);
//     if (r>0)
//     {
//         double angle = arg(x);
//         r = pow(r,n);
//         angle *= n;
//         return cdouble(r*cos(angle),r*sin(angle));
//     }
//     return cdouble();
// }
QuarticSolver::QuarticSolver()
{
    omg1 = cdouble(-0.5,sqrt(3)/2);
    omg2 = cdouble(-0.5,-sqrt(3)/2);
    SetUpEquation();
}
QuarticSolver::QuarticSolver(double a4, double a3, double a2, double a1, double a0)
{
    omg1 = cdouble(-0.5,sqrt(3)/2);
    omg2 = cdouble(-0.5,-sqrt(3)/2);
    SetUpEquation(a4,a3,a2,a1,a0);
}
void QuarticSolver::Solve(double a4, double a3, double a2, double a1, double a0)
{
    SetUpEquation(a4,a3,a2,a1,a0);
    Solve();
}
void QuarticSolver::SetUpEquation(double a4, double a3, double a2, double a1, double a0)
{
    a = cdouble(a4,0);
    b = cdouble(a3/a4,0);
    c = cdouble(a2/a4,0);
    d = cdouble(a1/a4,0);
    e = cdouble(a0/a4,0);

    P = (c*c + 12.0*e - 3.0*b*d)/9.0;
    Q = (27.0*d*d + 2.0*c*c*c + 27.0*b*b*e - 72.0*c*e - 9.0*b*c*d)/54.0;
    D = pow(Q*Q-P*P*P,0.5);
    u = Q + D;
    v = Q - D;
}
void QuarticSolver::Solve()
{
    if (abs(v) > abs(u))
    {
        u = pow(v,1.0/3.0);
    }
    else
    {
        u = pow(u,1.0/3.0);
    }
    if (abs(u) > 0.0)
    {
        v = P/u;
        cdouble &yMax = SOLUTIONS[0];
        double m2 = 0.0;
        double m2Max = 0.0;
        int iMax = -1.0;
        for (int i = 0; i < 3; i++)
        {
            y = u + v + c/3.0;
            u *= omg1;
            v *= omg2;
            a = b*b + 4.0*(y-c);
            m2 = norm(a);
            if (0 == i || m2Max < m2)
            {
                m2Max = m2;
                yMax = y;
                iMax = i;
            }
        }
        y = yMax;
    }
    else
    {
        y = c/3.0;
    }
    cdouble m = pow(b*b + 4.0*(y-c),0.5);
    if (norm(m) > __DBL_MIN__)
    {
        cdouble n = (b*y - 2.0*d)/m;
        a = pow((b+m)*(b+m)-8.0*(y+n),0.5);
        SOLUTIONS[0] = (-(b+m)+a)/4.0;
        SOLUTIONS[1] = (-(b+m)-a)/4.0;
        a = pow((b-m)*(b-m)-8.0*(y-n),0.5);
        SOLUTIONS[2] = (-(b-m)+a)/4.0;
        SOLUTIONS[3] = (-(b-m)-a)/4.0;
    }
    else
    {
        a = pow(b*b-8.0*y,0.5);
        SOLUTIONS[0] = (-b + a)/4.0;
        SOLUTIONS[1] = (-b + a)/4.0;
        SOLUTIONS[2] = (-b - a)/4.0;
        SOLUTIONS[3] = (-b - a)/4.0;
    }
}
vector<double> QuarticSolver::GetRealSolution()
{
    vector<double> res;
    for (int i = 0; i < 4; i++)
    {
        if (abs(SOLUTIONS[i].imag()) < 1e-8)
        {
            res.push_back(SOLUTIONS[i].real());
        }
    }
    return res;
}
