#include <cmath>
#include <iostream>

#include "Relaxation.h"

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
struct Spheroidal_Param {
    int m;
    int n;
    double c2;
    double gamma;
};

double plgndr(int l, int m, double x) {
    double fact, pll, pmm, pmmp1, somx2;

    pmm = 1.0;
    if (m > 0) {
        somx2 = sqrt((1 - x) * (1 + x));
        fact = 1.0;
        for (size_t i = 0; i < m; i++) {
            pmm *= -fact * somx2;
            fact += 2.0;
        }
    }
    if (l == m) {
        return pmm;
    } else {
        pmmp1 = x * (2 * m + 1) * pmm;
        if (l == (m + 1)) {
            return pmmp1;
        } else {
            for (size_t ll = m + 2; ll <= l; ll++) {
                pll = (x * (2 * ll - 1) * pmmp1 - (ll + m - 1) * pmm) / (ll - m);
                pmm = pmmp1;
                pmmp1 = pll;
            }
            return pll;
        }
    }
}

double fac(int n) {
    double res = 1.0;
    for (size_t i = 1; i <= n; i++) {
        res *= i;
    }
    return res;
}

class SpheroidalHarmODE : public BSM_Soliton::RelaxationODE {
public:
    SpheroidalHarmODE() : BSM_Soliton::RelaxationODE(3, 1), m(2), n(5), c2(1), gamma(1){};
    ~SpheroidalHarmODE(){};

    void Set_Parameter(int m_, int n_, double c2_) {
        m = m_;
        n = n_;
        c2 = c2_;
        gamma = pow(-1, m) * fac(n + m) / pow(2, m) / fac(m) / fac(n - m);
    }
    virtual void dYdX(BSM_Soliton::MeshPoint &point) override final {
        double x = point.X;
        VD &y = point.Y;
        point.Result[0] = y[1];
        point.Result[1] = (2 * x * (m + 1) * y[1] - (y[2] - c2 * x * x) * y[0]) / (1.0 - x * x);
        point.Result[2] = 0;

        point.dResultdY[0][0] = 0;
        point.dResultdY[0][1] = 1;
        point.dResultdY[0][2] = 0;

        point.dResultdY[1][0] = -(y[2] - c2 * x * x) / (1.0 - x * x);
        point.dResultdY[1][1] = 2 * x * (m + 1) / (1.0 - x * x);
        point.dResultdY[1][2] = -y[0] / (1.0 - x * x);

        point.dResultdY[2][0] = 0;
        point.dResultdY[2][1] = 0;
        point.dResultdY[2][2] = 0;
    }
    virtual void Left_Boundary_Constraints(BSM_Soliton::MeshPoint &point) override final {
        point.Result[0] = point.Y[0];
        point.dResultdY[0][0] = 1;
        point.dResultdY[0][1] = 0;
        point.dResultdY[0][2] = 0;
    }
    virtual void Right_Boundary_Constraints(BSM_Soliton::MeshPoint &point) override final {
        double x = point.X;
        VD &y = point.Y;

        point.Result[0] = y[0] - gamma;
        point.Result[1] = y[1] - y[0] * (y[2] - c2) / (2.0 * (m + 1));

        point.dResultdY[0][0] = 1;
        point.dResultdY[0][1] = 0;
        point.dResultdY[0][2] = 0;

        point.dResultdY[1][0] = -(y[2] - c2) / (2.0 * (m + 1));
        point.dResultdY[1][1] = 1;
        point.dResultdY[1][2] = -y[0] / (2.0 * (m + 1));
    }

private:
    int m, n;
    double c2, gamma;
};

int main(int argc, char const *argv[]) {
    int m = 2;
    int n = 5;
    double c2 = 1.0;
    double gamma = pow(-1, m) * fac(n + m) / pow(2, m) / fac(m) / fac(n - m);
    SpheroidalHarmODE fode;
    fode.Set_Parameter(m, n, c2);
    BSM_Soliton::Relaxation RX(&fode);
    VD y_beg = {0.1, 0.1, 20.0};
    VD y_end = {105, 2.0 / 3.0 * 105, 20.0};
    VD X_Guess;
    VVD Y_Guess;
    VD _Y(3);
    double h = 1.0 / 40.0;
    double fac1, fac2, deriv;
    for (size_t k = 0; k < 40; k++) {
        X_Guess.push_back(k * h);
        fac1 = 1.0 - k * h * k * h;
        fac2 = exp((-m / 2.0) * log(fac1));
        _Y[0] = plgndr(n, m, k * h) * fac2;
        deriv = -((n - m + 1) * plgndr(n + 1, m, k * h) - (n + 1) * k * h * plgndr(n, m, k * h)) / fac1;
        _Y[1] = m * k * h * _Y[0] / fac1 + deriv * fac2;
        _Y[2] = n * (n + 1) - m * (m + 1);
        Y_Guess.push_back(_Y);
    }
    X_Guess.push_back(1);
    _Y[0] = gamma;
    _Y[2] = n * (n + 1) - m * (m + 1);
    _Y[1] = (_Y[2] - c2) * _Y[0] / (2.0 * (m + 1));
    Y_Guess.push_back(_Y);
    RX.Solve(X_Guess, Y_Guess);
    RX.DumpSolution("Relaxation_test.dat");
    return 0;
}
