#include "THDM_CPC.h"
#include "QuarticSolver.h"


using namespace std;
using namespace Eigen;

THDM_CPC::THDM_CPC():Basic_Model(2)
{
    // The basic DOF is 2, rho1 and rho2
}

// bool THDM_CPC::Set_Potential_Parameters(double m112, double m222, double m122, double lam1, double lam2, double lam3, double lam4, double lam5)
// {
//     _m112 = m112;
//     _m222 = m222;
//     _m122 = m122;
//     _lam1 = lam1;
//     _lam2 = lam2;
//     _lam3 = lam3;
//     _lam4 = lam4;
//     _lam5 = lam5;
//     _Solved = false;

//     FindLocalMinima();

//     _v1 = 0;
//     _v2 = 0;
//     _IndexInput = -1;

//     for (size_t i = 0; i < _NLocalExtreme; i++)
//     {
//         if (_LocalMinimaQ[i] && _localExtreme[i][0] > _v1 && _localExtreme[i][1] > _v2)
//         {
//             _v1 = _localExtreme[i][0];
//             _v2 = _localExtreme[i][1];
//             _IndexInput = i;
//         }
//     }
//     _vev = GetMag(_v1,_v2);
//     _vev2 = _vev*_vev;
//     _beta = GetArg(_v1,_v2);

//     _tb = tan(_beta);
//     _sb = sin(_beta);
//     _cb = cos(_beta);


//     // ! Get the masses and alpha
//     blabla


//     if (_MHL2 < 0 || _MHH2 < 0 || _MHA2 < 0)
//     {
//         return false;
//     }
//     _MHL = sqrt(_MHL2);
//     _MHH = sqrt(_MHH2);
//     _MHA = sqrt(_MHA2);


//     return true;
// }

bool THDM_CPC::Set_Physical_Parameters(double beta, double alpha, double mhh, double mha, double mpm, double m122)
{
    _vev = vev;
    _vev2 = _vev*_vev;
    _beta = beta;
    _tb = tan(_beta);
    _sb = sin(_beta);
    _cb = cos(_beta);
    _v1 = _vev*_cb;
    _v2 = _vev*_sb;

    _MHL = MHL;
    _MHL2 = _MHL*_MHL;

    _MHH = mhh;
    _MHH2 = _MHH*_MHH;

    _MHA = mha;
    _MHA2 = _MHA*_MHA;
    
    _MHpm = mpm;
    _MHpm2 = _MHpm*_MHpm;

    _alpha = alpha;
    _ca = cos(_alpha);
    _sa = sin(_alpha);
    _cba = cos(_beta-_alpha);
    _sba = sin(_beta-_alpha);
    _m122 = m122;

    // * 1. Calculate lambda's
    _lam1 = (_cb*((_MHL2-_MHH2)*pow(_cba*_sb-_cb*_sba,2)+_MHH2)-_m122*_sb)/pow(_cb,3)/_vev2;
    _lam2 = (_sb*((_MHH2-_MHL2)*pow(_cba*_sb-_cb*_sba,2)+_MHL2)-_m122*_cb)/pow(_sb,3)/_vev2;
    _lam3 = ((_MHH2-_MHL2)*(_cba*_sb-_cb*_sba)*(_cb*_cba+_sb*_sba)+2.0*_cb*_sb*_MHpm2-_m122)/_cb/_sb/_vev2;
    _lam4 = (_cb*_sb*(_MHA2-2*_MHpm2)+_m122)/_cb/_sb/_vev2;
    _lam5 = (_m122-_cb*_sb*_MHA2)/_sb/_cb/_vev2;
    // * 2. Calculate m112 m222
    _m112 = _m122*_sb/_cb + _MHL2*_sba*(_cba*_sb-_cb*_sba)/2.0/_cb - _MHH2*_cba*(_cb*_cba+_sb*_sba)/2.0/_cb;
    _m222 = _m122*_cb/_sb - _MHL2*_sba*(_cb*_cba+_sb*_sba)/2.0/_sb + _MHH2*_cba*(_cb*_sba-_cba*_sb)/2.0/_sb;


    // cout<<_lam1<<"\t"<<_lam2<<"\t"<<_lam3<<"\t"<<_lam4<<"\t"<<_lam5<<endl;
    // cout<<_m112<<"\t"<<_m222<<"\t"<<_m122<<endl;

    _Solved = false;
    FindLocalMinima();
    _IndexInput = -1;
    for (size_t i = 0; i < _NLocalExtreme; i++)
    {
        if (CloseQ(_localExtreme[i],{_v1,_v2}))
        {
            _IndexInput = i;
        }
    }
    return true;
}


void THDM_CPC::FindLocalMinima()
{
    if (_Solved)
    {
        return;
    }

    Clear_Local_Cache();
    // cout<<"Cache Cleared"<<endl;
    // * 1. v1 = 0, v2 = 0
    _localExtreme.push_back({0,0});
    AppendLocalExtreme();
    // cout<<"Got First trivial solution"<<endl;

    // * 2. k = v2/v1
    // * k satisfy: lam2*m122*k^4 + (lam345*m222-lam2*m112)*k^3 + (lam1*m222-lam345*m112)*k - lam1*m122 = 0
    QuarticSolver _solver4(_lam2*_m122,(_lam3+_lam4+_lam5)*_m222-_lam2*_m112,0,_lam1*_m222-(_lam3+_lam4+_lam5)*_m112,-_lam1*_m122);
    _solver4.Solve();
    // cout<<"Solve the quartic equation"<<endl;
    VD k_solutions = _solver4.GetRealSolution();
    for (int i = 0; i < k_solutions.size(); i++)
    {
        double k = k_solutions[i];
        double _v12temp = (2*_m122*k-2*_m112)/((_lam3+_lam4+_lam5)*k*k+_lam1);
        if (_v12temp < 0) continue;
        _localExtreme.push_back({sqrt(_v12temp),k*sqrt(_v12temp)});
        AppendLocalExtreme();

        // ! temporarily ignore v1 < 0 solution, 
        // _localExtreme.push_back({-sqrt(_v12temp),-k*sqrt(_v12temp),0});
        // AppendLocalExtreme();
    }

    _Solved = true;
}

void THDM_CPC::PrintParameters()
{
    // Potential Parameter
    cout<<"Potential Parameters:"<<endl;
    cout<<"m112:\t"<<_m112<<endl;
    cout<<"m222:\t"<<_m222<<endl;
    cout<<"m122:\t"<<_m122<<endl;
    cout<<"lam1:\t"<<_lam1<<endl;
    cout<<"lam2:\t"<<_lam2<<endl;
    cout<<"lam3:\t"<<_lam3<<endl;
    cout<<"lam4:\t"<<_lam4<<endl;
    cout<<"lam5:\t"<<_lam5<<endl;
// Physical Parameter
    cout<<"Physical Parameters:"<<endl;
    cout<<"vev:\t"<<_vev<<endl;
    cout<<"beta:\t"<<_beta<<endl;
    cout<<"tb:\t"<<_tb<<endl;
    cout<<"v1:\t"<<_v1<<endl;
    cout<<"v2:\t"<<_v2<<endl;
    cout<<"m1:\t"<<_MHL<<endl;
    cout<<"m2:\t"<<_MHH<<endl;
    cout<<"m3:\t"<<_MHA<<endl;
    cout<<"m+-:\t"<<_MHpm<<endl;
    cout<<"alpha:\t"<<_alpha<<endl;
}

double THDM_CPC::Vtotal(VD field_values, double scale)
{
    // * This is the potential for CM Monopole
    double phi1 = field_values[0];
    double phi2 = field_values[1];
    double phi12 = phi1*phi1;
    double phi22 = phi2*phi2;
    double tmp = 0;
    
    tmp += _m122*pow(sqrt(_tb)*phi1-sqrt(1.0/_tb)*phi2,2)/2;
    tmp += _lam1*pow(phi12-_v1*_v1,2)/8;
    tmp += _lam2*pow(phi22-_v2*_v2,2)/8;
    tmp += (_lam3+_lam4+_lam5)*(phi12-_v1*_v1)*(phi22-_v2*_v2)/4;
    return tmp/pow(scale,4);
}
VD THDM_CPC::dVtotal(VD field_values, double scale)
{
    // * This is the potential for CM Monopole
    double phi1 = field_values[0];
    double phi2 = field_values[1];
    VD res(_N_VEVs);
    res[0] = _m122*(phi1*_sb/_cb-phi2) + _lam1*phi1*(phi1*phi1-_v1*_v1)/2 + (_lam3+_lam4+_lam5)*phi1*(phi2*phi2-_v2*_v2)/2;
    res[1] = _m122*(phi2*_cb/_sb-phi1) + _lam2*phi2*(phi2*phi2-_v2*_v2)/2 + (_lam3+_lam4+_lam5)*phi2*(phi1*phi1-_v1*_v1)/2;
    return res/pow(scale,3);
}
VVD THDM_CPC::d2Vtotal(VD field_values, double scale)
{
    // * This is the potential for CM Monopole
    double phi1 = field_values[0];
    double phi2 = field_values[1];
    VVD res(_N_VEVs,VD(_N_VEVs));
    res[0][0] = _m122*_sb/_cb + _lam1*(3*phi1*phi1-_v1*_v1)/2 + (_lam3+_lam4+_lam5)*(phi2*phi2-_v2*_v2)/2;
    res[0][1] = (_lam3+_lam4+_lam5)*phi1*phi2 - _m122;
    res[1][0] = res[0][1];
    res[1][1] = _m122*_cb/_sb + _lam2*(3*phi2*phi2-_v2*_v2)/2 + (_lam3+_lam4+_lam5)*(phi1*phi1-_v1*_v1)/2;
    return res/pow(scale,2);
}
double THDM_CPC::V0_global(double scale)
{
    return Vtotal({_v1,_v2},scale);
}
VD THDM_CPC::QuarticCoupling(VD field_values)
{
    return {_lam1/8,_lam2/8};
}
bool THDM_CPC::CheckUnitarity(double MAX)
{
    double a0Eigens[12];
    a0Eigens[0] = _lam3+_lam4;
    a0Eigens[1] = _lam3-_lam4;
    a0Eigens[2] = _lam3+_lam5;
    a0Eigens[3] = _lam3-_lam5;
    a0Eigens[4] = _lam3+2*_lam4+3*_lam5;
    a0Eigens[5] = _lam3+2*_lam4-3*_lam5;
    a0Eigens[6] = (_lam1+_lam2+sqrt(pow(_lam1-_lam2,2)+4*pow(_lam5,2)))/2;
    a0Eigens[7] = (_lam1+_lam2-sqrt(pow(_lam1-_lam2,2)+4*pow(_lam5,2)))/2;
    a0Eigens[8] = (_lam1+_lam2+sqrt(pow(_lam1-_lam2,2)+4*pow(_lam4,2)))/2;
    a0Eigens[9] = (_lam1+_lam2-sqrt(pow(_lam1-_lam2,2)+4*pow(_lam4,2)))/2;
    a0Eigens[10] = (3*(_lam1+_lam2)+sqrt(9*pow(_lam1-_lam2,2)+4*pow(2*_lam3+_lam4,2)))/2;
    a0Eigens[11] = (3*(_lam1+_lam2)-sqrt(9*pow(_lam1-_lam2,2)+4*pow(2*_lam3+_lam4,2)))/2;
    bool good = true;
    for (int i = 0; i < 12; i++)
    {
        good*=(abs(a0Eigens[i])<MAX*16.0*Pi);
        if (!good)
        {
            return good;
        }
    }
    return good;
}
bool THDM_CPC::CheckStability()
{
    if (_lam1 <= 0 || _lam2 <= 0)
    {
        return false;
    }
    
    if (_lam3 <= -sqrt(_lam1*_lam2) || _lam3+_lam4-_lam5 <= -sqrt(_lam1*_lam2))
    {
        return false;
    }

    return true;
}

string THDM_CPC::repr()
{
    return "vev\tbeta\talpha\tMHL\tMHH\tMHA\tMHpm\tm112\tm222\tm122\tlam1\tlam2\tlam3\tlam4\tlam5";
}
ostream &operator<<(std::ostream &os, const THDM_CPC &mod)
{
    os<<mod._vev<<"\t"<<mod._beta<<"\t"<<mod._alpha<<"\t"<<mod._MHL<<"\t"<<mod._MHH<<"\t"<<mod._MHA<<"\t"<<mod._MHpm<<"\t"<<mod._m112<<"\t"<<mod._m222<<"\t"<<mod._m122<<"\t"<<mod._lam1<<"\t"<<mod._lam2<<"\t"<<mod._lam3<<"\t"<<mod._lam4<<"\t"<<mod._lam5;
    return os;
}