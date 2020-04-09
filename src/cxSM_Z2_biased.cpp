#include "cxSM_Z2_biased.h"
#include "Constants.h"
#include <iostream>

using namespace std;


cxSM_Z2_biased::cxSM_Z2_biased()//:cxSM_Z2()
{
    // Set_Physical_Parameters(10,0.1,200);
    // Set_Physical_Parameters(10,0.1,200,250,0.1,0.1,0.1);
}

void cxSM_Z2_biased::Set_Potential_Parameters(double mu2, double lam, double del2, double b2, double d2, double del1, double a1, double c1, double c2)
{
    _mu2  = mu2;
    _lam  = lam;
    _del2 = del2;
    _b2   = b2;
    _d2   = d2;
    _del1 = del1;
    _a1   = a1;
    _c1   = c1;
    _c2   = c2;
    _Solved = false;

    FindLocalMinima();

    _vev = 0;
    _vs = 0;
    _IndexInput = -1;
    for (size_t i = 0; i < _NLocalExtreme; i++)
    {
        if (_localExtreme[i][1] > _vs && _LocalMinimaQ[i])
        {
            _vs         = _localExtreme[i][1];
            _vev        = _localExtreme[i][0];
            _IndexInput = i;
        }
    }

    double M11 = 2*_lam*_vev*_vev;
    double M12 = _del2*_vs*_vev/2 + _del1/2/sqrt(2)*_vev;
    double M22 = _d2*_vs*_vs/2 - sqrt(2)*_a1/_vs + (_c1 + _c2)/2/sqrt(2)*_vs - _del1/4/sqrt(2)*_vev*_vev/_vs;

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

    _MHA2 = -sqrt(2)*_a1/_vs - 3*_c1/2/sqrt(2)*_vs - _c2/6/sqrt(2)*_vs - _del1/4/sqrt(2)*_vev*_vev/_vs;
    _MHA = sqrt(_MHA2);
}

void cxSM_Z2_biased::Set_Physical_Parameters_del1_c1_c2(double vs, double theta, double MHH, double MHA, double del1, double c1, double c2)
{
    _MHL = MHL;
    _MHH = MHH;
    _MHA = MHA;

    _MHL2 = _MHL * _MHL;
    _MHH2 = _MHH * _MHH;
    _MHA2 = _MHA * _MHA;

    _vev   = vev;
    _vs    = vs;
    _theta = theta;

    _del1 = del1;
    _c1   = c1;
    _c2   = c2;

    double cth = cos(_theta);
    double sth = sin(_theta);

    _lam  = (_MHL2*cth*cth + _MHH2*sth*sth)/2/_vev/_vev;
    _del2 = (4*sth*cth*(_MHL2-_MHH2)-sqrt(2)*_del1*_vev)/2/_vs/_vev;
    _d2   = -2*(3*sqrt(2)*_c1*_vs+sqrt(2)*_c2*_vs-3*_MHL2*sth*sth-3*_MHH2*cth*cth+3*_MHA2)/3/_vs/_vs;
    _a1   = (-2*_vs*(_vs*(9*_c1+_c2)+6*sqrt(2)*_MHA2)-3*_del1*_vev*_vev)/24;
    _mu2  = -(4*_lam*_vev*_vev + _del2*_vs*_vs + sqrt(2)*_del1*_vs)/4;
    _b2   = -(8*sqrt(2)*_a1 + 2*sqrt(2)*(_c1+_c2)*_vs*_vs + 2*_d2*_vs*_vs*_vs + sqrt(2)*_del1*_vev*_vev + 2*_del2*_vev*_vev*_vs)/4/_vs;

    _Solved = false;
    FindLocalMinima();
    _IndexInput = -1;
    for (size_t i = 0; i < _NLocalExtreme; i++)
    {
        if (CloseQ(_localExtreme[i],{_vev,_vs}))
        {
            _IndexInput = i;
        }
    }
}

void cxSM_Z2_biased::Set_Physical_Parameters_a1_c1_c2(double vs, double theta, double MHH, double MHA, double a1, double c1, double c2)
{
    _MHL = MHL;
    _MHH = MHH;
    _MHA = MHA;

    _MHL2 = _MHL * _MHL;
    _MHH2 = _MHH * _MHH;
    _MHA2 = _MHA * _MHA;

    _vev   = vev;
    _vs    = vs;
    _theta = theta;

    _a1   = a1;
    _c1   = c1;
    _c2   = c2;

    double cth = cos(_theta);
    double sth = sin(_theta);

    _lam  = (_MHL2*cth*cth + _MHH2*sth*sth)/2/_vev/_vev;
    _del2 = (12*sqrt(2)*_a1+_vs*(sqrt(2)*_vs*(9*_c1+_c2)+12*_MHA2)+3*_vev*2*sth*cth*(_MHL2-_MHH2))/3/_vs/_vev/_vev;
    _d2   = -2*(3*sqrt(2)*_c1*_vs+sqrt(2)*_c2*_vs-3*_MHL2*sth*sth-3*_MHH2*cth*cth+3*_MHA2)/3/_vs/_vs;
    _del1   = -2*(12*_a1+_vs*(_vs*(9*_c1+_c2)+6*sqrt(2)*_MHA2))/3/_vev/_vev;
    _mu2  = -(4*_lam*_vev*_vev + _del2*_vs*_vs + sqrt(2)*_del1*_vs)/4;
    _b2   = -(8*sqrt(2)*_a1 + 2*sqrt(2)*(_c1+_c2)*_vs*_vs + 2*_d2*_vs*_vs*_vs + sqrt(2)*_del1*_vev*_vev + 2*_del2*_vev*_vev*_vs)/4/_vs;

    _Solved = false;
    FindLocalMinima();
    _IndexInput = -1;
    for (size_t i = 0; i < _NLocalExtreme; i++)
    {
        if (CloseQ(_localExtreme[i],{_vev,_vs}))
        {
            _IndexInput = i;
        }
    }
}

double cxSM_Z2_biased::Vtotal(VD field_values, double scale)
{
    double vh  = field_values[0];
    double vs  = field_values[1];
    double vh2 = vh*vh;
    double vs2 = vs*vs;

    double vtot = sqrt(2)*_a1*vs + _b2*vs2/4 + (_c1 + _c2)/6/sqrt(2)*vs2*vs + _d2*vs2*vs2/16 + _del1/4/sqrt(2)*vh2*vs + _lam*vh2*vh2/4 + _mu2*vh2/2 + _del2*vh2*vs2/8;
    return vtot/pow(scale,4);
}
VD cxSM_Z2_biased::QuarticCoupling(VD field_values)
{
    return {_lam, _d2/4.0};
}
VD cxSM_Z2_biased::dVtotal(VD field_values, double scale)
{
    double vh  = field_values[0];
    double vs  = field_values[1];
    double vh2 = vh*vh;
    double vs2 = vs*vs;

    VD res(2);
    res[0] = _mu2*vh + _lam*vh2*vh + _del2*vh*vs2/4 + _del1/2/sqrt(2)*vh*vs;
    res[1] = sqrt(2)*_a1 + _b2*vs/2 + _d2*vs2*vs/4 + _del2*vh2*vs/4 + (_c1 + _c2)/2/sqrt(2)*vs2 + _del1/4/sqrt(2)*vh2;

    return res/pow(scale,3); 
}
VVD cxSM_Z2_biased::d2Vtotal(VD field_values, double scale)
{
    double vh  = field_values[0];
    double vs  = field_values[1];
    double vh2 = vh*vh;
    double vs2 = vs*vs;

    VVD res(2,VD(2));

    res[0][0] = _mu2 + 3*_lam*vh2 + _del2*vs2/4 + _del1*vs/2/sqrt(2);
    res[0][1] = _del2*vh*vs/2 + _del1*vh/2/sqrt(2);
    res[1][0] = res[0][1];
    res[1][1] = _b2/2 + 3*_d2*vs2/4 + _del2*vh2/4 + (_c1 + _c2)/sqrt(2)*vs;
    
    return res/pow(scale,2);
}
double cxSM_Z2_biased::V0_global(double scale)
{
    return Vtotal({_vev,_vs},scale);
}

void cxSM_Z2_biased::SolveCubicEquation(double A[4], double *results, int &NSolution)
{
    NSolution = 0;
    _Solver.Solve(A[3],A[2],A[1],A[0]);
    if (_Solver.STATE==ONEREAL||_Solver.STATE==THREEEQUALREAL)
    {
        results[NSolution]=_Solver.SOLUTIONS[0];
        NSolution++;
    }
    else if (_Solver.STATE==ONETWO)
    {
        for (int i = 0; i < 2; ++i)
        {
            results[NSolution]=_Solver.SOLUTIONS[i];
            NSolution++;
        }
    }
    else if (_Solver.STATE==THREEREAL)
    {
        for (int i = 0; i < 3; ++i)
        {
            results[NSolution]=_Solver.SOLUTIONS[i];
            NSolution++;
        }
    }
    else
    {
#ifndef DEBUG
        std::cout<<"Warning in Cubic Solver of cxSM_Z2_biased"<<endl;
#endif
    }
}
void cxSM_Z2_biased::FindLocalMinima()
{
    if (_Solved)
    {
        return;
    }
    
    Clear_Local_Cache();

    int NCurrent;
    double results[3];
// ! First case is vh = 0, vs is solved from sqrt(2)*a1+b2/2*vs+(c1+c2)/2/sqrt(2)*vs^2+d2/4*vs^3 = 0
    double A[4] = {sqrt(2)*_a1,_b2/2,(_c1+_c2)/2/sqrt(2),_d2/4};
    SolveCubicEquation(A,results,NCurrent);
    for (size_t i = 0; i < NCurrent; i++)
    {
        _localExtreme.push_back({0,results[i]});
        AppendLocalExtreme();
    }
    
// ! Second case is vh != 0:
// ! vs is solved from:
// !    (8*a1*lam-del1*mu2)/4/sqrt(2)/lam - (-8*b2*lam+4*del2*mu2+del1*del1)/16/lam * vs + (8*lam*(c1+c2)-3*del2*del1)/16/sqrt(2)/lam * vs^2 - (del2*del2-4*d2*lam)/16/lam * vs^3 = 0
// ! vh is solve from:
// !    vh^2 = -(mu2+del1*vs/2/sqrt(2)+del2*vs^2/4)/lam
    A[0] = (8*_a1*_lam - _del1*_mu2)/4/sqrt(2)/_lam;
    A[1] = -(-8*_b2*_lam+4*_del2*_mu2+_del1*_del1)/16/_lam;
    A[2] = (8*_lam*(_c1+_c2)-3*_del2*_del1)/16/sqrt(2)/_lam;
    A[3] = -(_del2*_del2-4*_d2*_lam)/16/_lam;
    SolveCubicEquation(A,results,NCurrent);
    double vstemp;
    double vh2temp;
    for (size_t i = 0; i < NCurrent; i++)
    {
        vstemp = results[i];
        vh2temp = -(_mu2 + _del1/2/sqrt(2)*vstemp + _del2/4*vstemp*vstemp)/_lam;
        if (vh2temp < 0)
        {
            continue;
        }
        _localExtreme.push_back({sqrt(vh2temp),vstemp});
        AppendLocalExtreme();
    }
    
    _Solved=true;
}
bool cxSM_Z2_biased::GetBiasedMirrorMinimum(VD &mirror)
{
    VD current;
    VD localtemp;
    VD guessed = {_vev, -_vs};
    double dis_cur;
    double dis_tmp;
    bool good=false;
    int id;
    if (_NLocalMinima < 2) return good;
    for (size_t i = 0; i < _NLocalMinima; i++)
    {
        id = _MinimaIndex[i];
        if (_IndexInput == id)
        {
            continue;
        }
        if (!good)
        {
            current = _localExtreme[id];
            good = true;
        }
        else
        {
            localtemp = _localExtreme[id];
            dis_cur = sqrt((current-guessed)*(current-guessed));
            dis_tmp = sqrt((localtemp-guessed)*(localtemp-guessed));
            if (dis_tmp < dis_cur)
            {
                current = localtemp;
            }
        }
    }
    mirror = current;
    return good;
}
void cxSM_Z2_biased::PrintParameters()
{
    cout<<"Potential Parameter: "<<endl;
    cout<<"mu2:\t"<<_mu2<<endl;
    cout<<"b2:\t"<<_b2<<endl;
    cout<<"lam:\t"<<_lam<<endl;
    cout<<"del2:\t"<<_del2<<endl;
    cout<<"d2:\t"<<_d2<<endl;
    cout<<"del1:\t"<<_del1<<endl;
    cout<<"a1:\t"<<_a1<<endl;
    cout<<"c1:\t"<<_c1<<endl;
    cout<<"c2:\t"<<_c2<<endl;

    cout<<"Physical Parameter: "<<endl;
    cout<<"vev:\t"<<_vev<<endl;
    cout<<"vs:\t"<<_vs<<endl;
    cout<<"theta:\t"<<_theta<<endl;
    cout<<"MHL:\t"<<_MHL<<endl;
    cout<<"MHH:\t"<<_MHH<<endl;
    cout<<"MHA:\t"<<_MHA<<endl;
}