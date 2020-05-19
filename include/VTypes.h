#ifndef VTypes_H
#define VTypes_H
#include <cmath>
#include <vector>
#include <iostream>
#include <functional>

typedef std::vector<bool> VB;
typedef std::vector<int> VI;
typedef std::vector<double> VD;
typedef std::vector<std::vector<double> > VVD;
typedef std::vector<std::vector<std::vector<double> > > VVVD;

typedef std::function<double(double)> V1D;
typedef std::function<double(VD,double)> VnD; // n-dimension potential
typedef std::function<VD(VD,double)> dVnD; // n-dimension potential derivative

VD abs(const VD &input);
VD pow(const VD &input, double power);

VD operator+(const VD &lhs, const VD &rhs);
VD operator+(const VD &lhs, const double &cons);
VD operator+(const double &cons, const VD &rhs);
VD operator+(const VD &lhs, const VD &rhs);


VD operator-(const VD &lhs, const VD &rhs);
VD operator-(const VD &lhs, const double &rhs);
VD operator-(const double &lhs, const VD &rhs);
VD operator-(const VD &rhs);

VD operator*(const VD &lhs, const double &s);
VD operator*(const double &s, const VD &rhs);
double operator*(const VD &lhs, const VD &rhs); // Scalar Product

VD operator/(const VD &lhs, const VD &rhs); // elementary-wise divide
VD operator/(const VD &lhs, const double &s);

std::ostream& operator<<(std::ostream& out, const VD& s);

VVD operator*(const VVD &lhs, const double &s);
VVD operator*(const double &s, const VVD &rhs);
VVD operator/(const VVD &lhs, const double &s);
VVD operator+(const VVD &lhs, const VVD &rhs);
VVD operator-(const VVD &lhs, const VVD &rhs);
VD operator*(const VVD &lhs, const VVD &rhs);
VVD operator*(const VVD &lhs, const VD &rhs);
VVD operator*(const VD &lhs, const VVD &rhs);

double Simpson(VD X, VD Y);

VD cumtrapz(VD pts, VD x = {}, double dx = 1.0, double initial=0.0);

VVD transpose(const VVD &mat);

VD linspace(double start, double end, int n);
VD linspace(double start, double end, double dx);

template<class T>
T deriv14_const_dx(T y, double dx = 1.0)
{
    int N = y.size();
    T dy;
    dy.push_back(-25.0*y[0]+48.0*y[1]-36.0*y[2]+16.0*y[3]-3*y[4]);
    dy.push_back(-3.0*y[0]-10.0*y[1]+18.0*y[2]-6.0*y[3]+y[4]);
    for (int i = 2; i < N-2; i++)
    {
        dy.push_back(y[i-2]-8.0*y[i-1]+8.0*y[i+1]-y[i+2]);
    }
    dy.push_back(3.0*y[N-1]+10.0*y[N-2]-18.0*y[N-3]+6.0*y[N-4]-y[N-5]);
    dy.push_back(25.0*y[N-1]-48.0*y[N-2]+36.0*y[N-3]-16.0*y[N-4]+3.0*y[N-5]);

    return dy/(12.0*dx);
}

bool CloseQ(double x1, double x2);

bool CloseQ(VD x1, VD x2);


double GetMag(double vre, double vim);
double GetArg(double vre, double vim);
double GetAngle(double sin_th, double cos_th);

#endif