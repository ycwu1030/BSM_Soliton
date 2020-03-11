#include "VTypes.h"
#include <cmath>
using namespace std;

bool CloseQ(double x1, double x2)
{
    return abs(x1-x2)<1e-2;
}
bool CloseQ(VD x1, VD x2)
{
    VD difs = x1-x2;
    double distance = sqrt(difs*difs);
    return distance<1e-2;
}
vector<double> operator+(const vector<double> &lhs, const vector<double> &rhs){
    vector<double> res;
    for (int i = 0; i < lhs.size(); ++i)
    {
        res.push_back(lhs[i]+rhs[i]);
    }
    return res;
}

vector<double> operator+(const vector<double> &lhs, const double &cons){
    vector<double> res;
    for (int i = 0; i < lhs.size(); ++i)
    {
        res.push_back(lhs[i]+cons);
    }
    return res;
}

vector<double> operator+(const double &cons, const vector<double> &rhs){
    vector<double> res;
    for (int i = 0; i < rhs.size(); ++i)
    {
        res.push_back(cons+rhs[i]);
    }
    return res;
}


vector<double> operator-(const vector<double> &lhs, const vector<double> &rhs){
    vector<double> res;
    for (int i = 0; i < lhs.size(); ++i)
    {
        res.push_back(lhs[i]-rhs[i]);
    }
    return res;
}


vector<double> operator*(const vector<double> &lhs, const double &s)
{
    vector<double> res;
    for (int i = 0; i < lhs.size(); ++i)
    {
        res.push_back(lhs[i]*s);
    }
    return res;
}


vector<double> operator*(const double &s, const vector<double> &rhs)
{
    vector<double> res;
    for (int i = 0; i < rhs.size(); ++i)
    {
        res.push_back(rhs[i]*s);
    }
    return res;
}


double operator*(const vector<double> &lhs, const vector<double> &rhs)
{
    double res = 0;
    for (int i = 0; i < lhs.size(); ++i)
    {
        res += lhs[i]*rhs[i];
    }
    return res;
}

vector<double> operator/(const vector<double> &lhs, const double &s)
{
    vector<double> res;
    for (int i = 0; i < lhs.size(); ++i)
    {
        res.push_back(lhs[i]/s);
    }
    return res;
}

vector<double> operator/(const vector<double> &lhs, const vector<double> &rhs)
{
    vector<double> res;
    for (int i = 0; i < lhs.size(); ++i)
    {
        res.push_back(lhs[i]/rhs[i]);
    }
    return res;
}

vector<double> abs(const vector<double> &input)
{
    vector<double> res;
    for (int i = 0; i < input.size(); ++i)
    {
        res.push_back(abs(input[i]));
    }
    return res;
}

VVD operator*(const VVD &lhs, const double &s)
{
    VVD res;
    for (int i = 0; i < lhs.size(); ++i)
    {
        res.push_back(lhs[i]*s);
    }
    return res;
}
VVD operator*(const double &s, const VVD &rhs)
{
    VVD res;
    for (int i = 0; i < rhs.size(); ++i)
    {
        res.push_back(rhs[i]*s);
    }
    return res;
}
VVD operator/(const VVD &lhs, const double &s)
{
    VVD res;
    for (int i = 0; i < lhs.size(); ++i)
    {
        res.push_back(lhs[i]/s);
    }
    return res;
}

ostream& operator<<(ostream& out, const VD& s)
{
    for (size_t i = 0; i < s.size()-1; i++)
    {
        out<<s[i]<<"\t";
    }
    out<<s[s.size()-1];
    return out;
}