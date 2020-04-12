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

VD abs(const VD &input);

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

VVD transpose(const VVD &mat);

bool CloseQ(double x1, double x2);

bool CloseQ(VD x1, VD x2);

double GetAngle(double sin_th, double cos_th);

#endif