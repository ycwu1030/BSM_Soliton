#include "ToyD4.h"
#include <iostream>

using namespace std;

ToyD4::ToyD4():Basic_Model(2)
{

}

void ToyD4::Set_Potential_Parameters(double g)
{

    _g = g;

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
            _v2global        = _localExtreme[i][1];
            _v1global        = _localExtreme[i][0];
            _IndexInput = i;
            Vtmp = _Vtotal[i];
        }
    }

}

double ToyD4::Vtotal(VD field_values, double scale)
{
    double f1 = field_values[0];
    double f2 = field_values[1];
    double f12 = f1*f1;
    double f22 = f2*f2;
    double f14 = f12*f12;
    double f24 = f22*f22;

    double vtot = -(f12+f22)/2.0 + (f14+f24)/4.0 + (1.0+_g)*f12*f22/2.0;
    return vtot/pow(scale,4);
}
VD ToyD4::QuarticCoupling(VD field_values)
{
    return {1.0, 1.0};
}
VD ToyD4::dVtotal(VD field_values, double scale)
{
    double f1 = field_values[0];
    double f2 = field_values[1];

    VD res(2);
    res[0] = (-1.0 + (f1*f1+f2*f2) + _g*f2*f2/2.0)*f1;
    res[1] = (-1.0 + (f1*f1+f2*f2) + _g*f1*f1/2.0)*f2;

    return res/pow(scale,3);
}
VVD ToyD4::d2Vtotal(VD field_values, double scale)
{
    double f1 = field_values[0];
    double f2 = field_values[1];

    VVD res(2, VD(2));
    res[0][0] = -1.0 + 3.0*f1*f1 + (1.0+_g)*f2*f2;
    res[0][1] = 2.0*(1.0+_g)*f1*f2;
    res[1][0] = res[0][1];
    res[1][1] = -1.0 + 3.0*f2*f2 + (1.0+_g)*f1*f1;

    return res/pow(scale,2); 
}
double ToyD4::V0_global(double scale)
{
    return Vtotal({_v1global,_v2global},scale);
}
void ToyD4::FindLocalMinima()
{
    if (_Solved)
    {
        return;
    }
    
    Clear_Local_Cache();

    _localExtreme.push_back({0,0});
    AppendLocalExtreme();

    _localExtreme.push_back({0,1});
    AppendLocalExtreme();

    _localExtreme.push_back({0,-1});
    AppendLocalExtreme();

    _localExtreme.push_back({1,0});
    AppendLocalExtreme();

    _localExtreme.push_back({-1,0});
    AppendLocalExtreme();

    if (2.0 + _g > 0)
    {
        double tmp = 1.0/sqrt(2.0+_g);
        _localExtreme.push_back({tmp,tmp});
        AppendLocalExtreme();

        _localExtreme.push_back({tmp,-tmp});
        AppendLocalExtreme();

        _localExtreme.push_back({-tmp,tmp});
        AppendLocalExtreme();

        _localExtreme.push_back({-tmp,-tmp});
        AppendLocalExtreme();
    }
    

    _Solved = true;

}

void ToyD4::PrintParameters()
{
    cout<<"Potential Parameter:"<<endl;
    cout<<"g:\t"<<_g<<endl;

    cout<<"Physical Parameter:"<<endl;
    cout<<"v1:\t"<<_v1<<endl;
    cout<<"v2:\t"<<_v2<<endl;
}