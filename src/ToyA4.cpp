#include "ToyA4.h"
#include <iostream>

using namespace std;

ToyA4::ToyA4():Basic_Model(2)
{

}

void ToyA4::Set_Potential_Parameters(double mu2, double g1, double g2)
{
    _mu2 = mu2;
    _g1 = g1;
    _g2 = g2;

    _Solved = false;
    FindLocalMinima();

    _v1global = 0;
    _v2global = 0;
    _v1 = 0;
    _v2 = 0;
    _IndexInput = -1;
    double Vtmp=0;
    for (size_t i = 0; i < _NLocalExtreme; i++)
    {
        if (_LocalMinimaQ[i])
        {
            if (abs(_localExtreme[i][0])>_v1)
            {
                _v1 = abs(_localExtreme[i][0]);
            }
            if (abs(_localExtreme[i][1])>_v2)
            {
                _v2 = abs(_localExtreme[i][1]);
            }
        }
        if (_Vtotal[i] < Vtmp && _LocalMinimaQ[i])
        {
            _v2global         = _localExtreme[i][1];
            _v1global        = _localExtreme[i][0];
            _IndexInput = i;
            Vtmp = _Vtotal[i];
        }
    }

}

double ToyA4::Vtotal(VD field_values, double scale)
{
    double chi1 = field_values[0];
    double chi2 = field_values[1];
    double I1 = chi1*chi1 + chi2*chi2;
    double I2 = chi1*chi1*chi2*chi2;

    double vtot = _mu2*I1/2.0 + _g1*I1*I1/4.0 + _g2*I2/4.0;
    return vtot/pow(scale,4);
}
VD ToyA4::QuarticCoupling(VD field_values)
{
    return {_g1, _g1};
}
VD ToyA4::dVtotal(VD field_values, double scale)
{
    double chi1 = field_values[0];
    double chi2 = field_values[1];

    VD res(2);
    res[0] = _mu2*chi1 + _g1*chi1*(chi1*chi1+chi2*chi2) + _g2*chi1*chi2*chi2/2.0;
    res[1] = _mu2*chi2 + _g1*chi2*(chi1*chi1+chi2*chi2) + _g2*chi1*chi1*chi2/2.0;

    return res/pow(scale,3);
}
VVD ToyA4::d2Vtotal(VD field_values, double scale)
{
    double chi1 = field_values[0];
    double chi2 = field_values[1];

    VVD res(2, VD(2));
    res[0][0] = _mu2 + 3.0*_g1*chi1*chi1 + (2.0*_g1+_g2)*chi2*chi2/2.0;
    res[0][1] = (2.0*_g1+_g2)*chi1*chi2;
    res[1][0] = res[0][1];
    res[1][1] = _mu2 + 3.0*_g1*chi2*chi2 + (2.0*_g1+_g2)*chi1*chi1/2.0;

    return res/pow(scale,2); 
}
double ToyA4::V0_global(double scale)
{
    return Vtotal({_v1global,_v2global},scale);
}
void ToyA4::FindLocalMinima()
{
    if (_Solved)
    {
        return;
    }
    
    Clear_Local_Cache();

    _localExtreme.push_back({0,0});
    AppendLocalExtreme();

    if (-_mu2/_g1 >= 0)
    {
        _localExtreme.push_back({0,sqrt(-_mu2/_g1)});
        AppendLocalExtreme();

        _localExtreme.push_back({0,-sqrt(-_mu2/_g1)});
        AppendLocalExtreme();

        _localExtreme.push_back({sqrt(-_mu2/_g1),0});
        AppendLocalExtreme();

        _localExtreme.push_back({-sqrt(-_mu2/_g1),0});
        AppendLocalExtreme();
    }

    if (_mu2/(-4*_g1-_g2) >= 0)
    {
        _localExtreme.push_back({sqrt(2*_mu2/(-4*_g1-_g2)),sqrt(2*_mu2/(-4*_g1-_g2))});
        AppendLocalExtreme();

        _localExtreme.push_back({sqrt(2*_mu2/(-4*_g1-_g2)),-sqrt(2*_mu2/(-4*_g1-_g2))});
        AppendLocalExtreme();

        _localExtreme.push_back({-sqrt(2*_mu2/(-4*_g1-_g2)),sqrt(2*_mu2/(-4*_g1-_g2))});
        AppendLocalExtreme();

        _localExtreme.push_back({-sqrt(2*_mu2/(-4*_g1-_g2)),-sqrt(2*_mu2/(-4*_g1-_g2))});
        AppendLocalExtreme();
    }
    

    _Solved = true;

}

void ToyA4::PrintParameters()
{
    cout<<"Potential Parameter:"<<endl;
    cout<<"mu2:\t"<<_mu2<<endl;
    cout<<"g1:\t"<<_g1<<endl;
    cout<<"g2:\t"<<_g2<<endl;

    cout<<"Physical Parameter:"<<endl;
    cout<<"v1:\t"<<_v1<<endl;
    cout<<"v2:\t"<<_v2<<endl;
}