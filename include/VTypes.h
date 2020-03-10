#ifndef VTypes_H
#define VTypes_H
#include <cmath>
#include <vector>

typedef std::vector<bool> VB;
typedef std::vector<double> VD;
typedef std::vector<std::vector<double> > VVD;
typedef std::vector<std::vector<std::vector<double> > > VVVD;

VD abs(const VD &input);

VD operator+(const VD &lhs, const VD &rhs);
VD operator+(const VD &lhs, const double &cons);
VD operator+(const double &cons, const VD &rhs);
VD operator+(const VD &lhs, const VD &rhs);

VD operator-(const VD &lhs, const VD &rhs);

VD operator*(const VD &lhs, const double &s);
VD operator*(const double &s, const VD &rhs);
double operator*(const VD &lhs, const VD &rhs); // Scalar Product

VD operator/(const VD &lhs, const VD &rhs); // elementary-wise divide
VD operator/(const VD &lhs, const double &s);

VVD operator*(const VVD &lhs, const double &s);
VVD operator*(const double &s, const VVD &rhs);

VVD operator/(const VVD &lhs, const double &s);

bool CloseQ(double x1, double x2);

#endif