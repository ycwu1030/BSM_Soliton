#ifndef VTypes_H
#define VTypes_H
#include <vector>

typedef std::vector<bool> VB;
typedef std::vector<double> VD;
typedef std::vector<std::vector<double> > VVD;
typedef std::vector<std::vector<std::vector<double> > > VVVD;

std::vector<double> abs(const std::vector<double> &input);

std::vector<double> operator+(const std::vector<double> &lhs, const std::vector<double> &rhs);
std::vector<double> operator+(const std::vector<double> &lhs, const double &cons);
std::vector<double> operator+(const double &cons, const std::vector<double> &rhs);
std::vector<double> operator+(const std::vector<double> &lhs, const std::vector<double> &rhs);

std::vector<double> operator-(const std::vector<double> &lhs, const std::vector<double> &rhs);

std::vector<double> operator*(const std::vector<double> &lhs, const double &s);
std::vector<double> operator*(const double &s, const std::vector<double> &rhs);
double operator*(const std::vector<double> &lhs, const std::vector<double> &rhs); // Scalar Product

std::vector<double> operator/(const std::vector<double> &lhs, const std::vector<double> &rhs); // elementary-wise divide
std::vector<double> operator/(const std::vector<double> &lhs, const double &s);

#endif