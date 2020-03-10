#include "CubicSolver.h"
#include "Constants.h"
#include <cmath>
// #include <bits/stdc++.h>

using namespace std;

CubicSolver::CubicSolver()
{
    SetUpEquation();
    // Solve();
}

CubicSolver::CubicSolver(double a3, double a2, double a1, double a0)
{
    SetUpEquation(a3,a2,a1,a0);
    // Solve();
}

void CubicSolver::SetUpEquation(double a3, double a2, double a1, double a0)
{
    A0 = a0;
    A1 = a1;
    A2 = a2;
    A3 = a3;

    b = A2/(3.0*A3);
    c = A1/(6.0*A3);
    d = A0/(2.0*A3);

    alpha = -pow(b,3) + 3.0*b*c - d;
    beta = pow(b,2) - 2.0*c;
    Delta = pow(alpha,2) - pow(beta,3);
}

void CubicSolver::Solve(double a3, double a2, double a1, double a0)
{
    SetUpEquation(a3, a2, a1, a0);
    Solve();
}

void CubicSolver::Solve()
{
    if (Delta > 0)
    {
        STATE = ONEREAL;
        R1 = cbrt(alpha+sqrt(Delta));
        R2 = cbrt(alpha-sqrt(Delta));
        SOLUTIONS[0] = - b + R1 + R2;
        SOLUTIONS[1] = -999;
        SOLUTIONS[2] = -999;
    }
    else if (Delta == 0)
    {
        R1=cbrt(alpha);
        R2=R1;
        SOLUTIONS[0] = -b + 2*R1;
        SOLUTIONS[1] = -b - R1;
        SOLUTIONS[2] = -b - R1;
        if (R1 == 0)
        {
            STATE = THREEEQUALREAL;
        }
        else
        {
            STATE = ONETWO;
        }
    }
    else
    {
        STATE = THREEREAL;
        theta = acos(alpha/sqrt(pow(beta,3)));
        SOLUTIONS[0] = -b + 2*sqrt(beta)*cos(theta/3.0);
        SOLUTIONS[1] = -b + 2*sqrt(beta)*cos(theta/3.0+2.0*Pi/3.0);
        SOLUTIONS[2] = -b + 2*sqrt(beta)*cos(theta/3.0-2.0*Pi/3.0); 
    }
}