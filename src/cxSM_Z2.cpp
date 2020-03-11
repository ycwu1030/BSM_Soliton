#include "cxSM_Z2.h"
#include "Constants.h"
#include <iostream>

using namespace std;


cxSM_Z2::cxSM_Z2():Basic_Model(2)
{
    // Set_Physical_Parameters(10,0.1,200);
}

void cxSM_Z2::Set_Potential_Parameters(double mu2, double lam, double del2, double b2, double d2)
{
    _mu2  = mu2;
    _lam  = lam;
    _del2 = del2;
    _b2   = b2;
    _d2   = d2;

    if ((2*_b2*_del2-4*_d2*_mu2)/(4*_d2*_lam-_del2*_del2) >= 0 )
    {
        _vev = sqrt((2*_b2*_del2-4*_d2*_mu2)/(4*_d2*_lam-_del2*_del2));
    }
    else
    {
        _vev = sqrt(abs((2*_b2*_del2-4*_d2*_mu2)/(4*_d2*_lam-_del2*_del2)));
        cerr<<"[Error]: Wrong potential parameters given negative vev^2"<<endl;
    }
    
    if ((_del2*_mu2-2*_b2*_lam)/(4*_d2*_lam-_del2*_del2) >= 0)
    {
        _vs = 2*sqrt((_del2*_mu2-2*_b2*_lam)/(4*_d2*_lam-_del2*_del2));
    }
    else
    {
        _vs = 2*sqrt(abs((_del2*_mu2-2*_b2*_lam)/(4*_d2*_lam-_del2*_del2)));
        cerr<<"[Error]: Wrong potential parameters given negative vs^2"<<endl;
    }

    double M11 = 2*_lam*_vev*_vev;
    double M12 = _del2*_vs*_vev/2;
    double M22 = _d2*_vs*_vs/2;

    _MHL2 = (M11+M22 - sqrt(pow(M11-M22,2)+4*M22*M22))/2;
    _MHH2 = (M11+M22 + sqrt(pow(M11-M22,2)+4*M22*M22))/2;

    _MHL = sqrt(_MHL2);
    _MHH = sqrt(_MHH2);

    double S2A = 2*M12/(_MHL2 - _MHH2);
    double C2A = (M11-M22)/(_MHL2 - _MHH2);

    _theta = asin(S2A)/2.0;
    if (C2A < 0)
    {
        if (S2A < 0)
        {
            _theta = Pi/2.0 - _theta;
        }
        else
        {
            _theta = -Pi/2.0 - _theta;
        }
    }
    _Solved = false;
}

void cxSM_Z2::Set_Physical_Parameters(double vs, double theta, double MHH)
{
    _MHL = MHL;
    _MHH = MHH;

    _MHL2 = _MHL * _MHL;
    _MHH2 = _MHH * _MHH;

    _vev   = vev;
    _vs    = vs;
    _theta = theta;

    double cth = cos(_theta);
    double sth = sin(_theta);

    _lam  = (_MHL2*cth*cth + _MHH2*sth*sth)/2/_vev/_vev;
    _del2 = 2*sth*cth*(_MHL2-_MHH2)/_vs/_vev;
    _d2   = 2*(_MHL2*sth*sth+_MHH2*cth*cth)/_vs/_vs;
    _mu2  = -(4*_lam*_vev*_vev + _del2*_vs*_vs)/4;
    _b2   = -(_d2*_vs*_vs + _del2*_vev*_vev)/2;

    _Solved = false;
}

double cxSM_Z2::Vtotal(VD field_values, double scale)
{
    double vh  = field_values[0];
    double vs  = field_values[1];
    double vh2 = vh*vh;
    double vs2 = vs*vs;

    double vtot = _mu2*vh2/2 + _b2*vs2/4 + _lam*vh2*vh2/4 + _d2*vs2*vs2/16 + _del2*vh2*vs2/8;
    return vtot/pow(scale,4);
}
VD cxSM_Z2::dVtotal(VD field_values, double scale)
{
    double vh  = field_values[0];
    double vs  = field_values[1];
    double vh2 = vh*vh;
    double vs2 = vs*vs;

    VD res(2);
    res[0] = _mu2*vh + _lam*vh2*vh + _del2*vh*vs2/4;
    res[1] = _b2*vs/2 + _d2*vs2*vs/4 + _del2*vh2*vs/4;

    return res/pow(scale,3); 
}
VVD cxSM_Z2::d2Vtotal(VD field_values, double scale)
{
    double vh  = field_values[0];
    double vs  = field_values[1];
    double vh2 = vh*vh;
    double vs2 = vs*vs;

    VVD res(2,VD(2));

    res[0][0] = (_mu2 + 3*_lam*vh2 + _del2*vs2/4)/pow(scale,2);
    res[0][1] = (_del2*vh*vs/2)/pow(scale,2);
    res[1][0] = res[0][1];
    res[1][1] = (_b2/2 + 3*_d2*vs2/4 + _del2*vh2/4)/pow(scale,2);
    
    return res;
}
double cxSM_Z2::V0_global(double scale)
{
    return Vtotal({_vev,_vs},scale);
}

bool cxSM_Z2::CheckStability()
{
    return (_lam>0)&&(_d2>0)&&((_lam*_d2>_del2*_del2)||_del2>0);
}

bool cxSM_Z2::CheckUnitarity(double MAX)
{
    double EigenA0[4];
    EigenA0[0] = _d2/2;
    EigenA0[1] = 2*_lam;
    EigenA0[2] = (_d2+6*_lam-sqrt(_d2*_d2-12*_d2*_lam+2*_del2*_del2+36*_lam*_lam))/2;
    EigenA0[3] = (_d2+6*_lam+sqrt(_d2*_d2-12*_d2*_lam+2*_del2*_del2+36*_lam*_lam))/2;
    bool good=true;
    for (int i = 0; i < 4; ++i)
    {
        good*=(abs(EigenA0[i])<MAX*16.0*Pi);
        if (!good)
        {
            return good;
        }
    }

    return good;
}

void cxSM_Z2::FindLocalMinima()
{
    if (_Solved)
    {
        return;
    }
    
    _NLocalMinima=0; // Track we are getting which solution
    _NLocalExtreme=0;
    _IndexInput = -1;

// ! First solution is vh = 0, vs = 0
    _localExtreme.push_back({0,0});
    AppendLocalExtreme();

// ! Second solution is vh = sqrt(-mu2/lam), vs = 0, if possible;
    if (-_mu2/_lam >= 0)
    {
        _localExtreme.push_back({sqrt(-_mu2/_lam),0});
        AppendLocalExtreme();
    }

// ! Third solution is vh = 0, vs = sqrt(-2*b2/d2), if possible;
    if (-2*_b2/_d2 >= 0)
    {
        _localExtreme.push_back({0,sqrt(-2*_b2/_d2)});
        AppendLocalExtreme();
    }
    
// ! The last one is vh = sqrt((2*b2*del2-4*d2*mu2)/(4*d2*lam-del2*del2)), vs = 2*sqrt((del2*mu2-2*b2*lam)/(4*d2*lam-del2*del2))
    if ((2*_b2*_del2-4*_d2*_mu2)/(4*_d2*_lam-_del2*_del2) >= 0 && (_del2*_mu2-2*_b2*_lam)/(4*_d2*_lam-_del2*_del2) >= 0)
    {
        _localExtreme.push_back({sqrt((2*_b2*_del2-4*_d2*_mu2)/(4*_d2*_lam-_del2*_del2)),2*sqrt((_del2*_mu2-2*_b2*_lam)/(4*_d2*_lam-_del2*_del2))});
        if (CloseQ(_localExtreme[_NLocalExtreme],{_vev,_vs}))
        {
            _IndexInput = _NLocalExtreme;
        }
        AppendLocalExtreme();
    }
    _Solved=true;
}

void cxSM_Z2::PrintParameters()
{
    cout<<"Potential Parameter: "<<endl;
    cout<<"mu2:\t"<<_mu2<<endl;
    cout<<"b2:\t"<<_b2<<endl;
    cout<<"lam:\t"<<_lam<<endl;
    cout<<"del2:\t"<<_del2<<endl;
    cout<<"d2:\t"<<_d2<<endl;

    cout<<"Physical Parameter: "<<endl;
    cout<<"vev:\t"<<_vev<<endl;
    cout<<"vs:\t"<<_vs<<endl;
    cout<<"theta:\t"<<_theta<<endl;
    cout<<"MHL:\t"<<_MHL<<endl;
    cout<<"MHH:\t"<<_MHH<<endl;
}