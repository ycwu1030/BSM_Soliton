#include "THDM_SCPV1.h"
#include "QuarticSolver.h"


using namespace std;
using namespace Eigen;

THDM_SCPV1::THDM_SCPV1():Basic_Model(3)
{
    // The basic DOF is 3
}

bool THDM_SCPV1::Set_Potential_Parameters(double m112, double m222, double m122, double lam1, double lam2, double lam3, double lam4, double lam5)
{
    _m112 = m112;
    _m222 = m222;
    _m122 = m122;
    _lam1 = lam1;
    _lam2 = lam2;
    _lam3 = lam3;
    _lam4 = lam4;
    _lam5 = lam5;
    _Solved = false;

    FindLocalMinima();

    _v1 = 0;
    _v2r = 0;
    _v2i = 0;
    _IndexInput = -1;

    for (size_t i = 0; i < _NLocalExtreme; i++)
    {
        if (_LocalMinimaQ[i] && _localExtreme[i][0] > _v1 && _localExtreme[i][1] > _v2r)
        {
            _v1 = _localExtreme[i][0];
            _v2r = _localExtreme[i][1];
            _v2i = _localExtreme[i][2];
            _IndexInput = i;
        }
    }
    _v2 = GetMag(_v2r,_v2i);
    _theta = GetArg(_v2r,_v2i);
    _vev = GetMag(_v1,_v2);
    _vev2 = _vev*_vev;
    _beta = GetArg(_v1,_v2);

    _tb = tan(_beta);
    _sb = sin(_beta);
    _cb = cos(_beta);

    _sth = sin(_theta);
    _cth = cos(_theta);

    Matrix3d MM2;
    MM2(0,0) = _lam1*_cb*_cb + _lam5*_cth*_cth*_sb*_sb;
    MM2(0,1) = (_lam3+_lam4-_lam5*_sth*_sth)*_sb*_cb;
    MM2(0,2) = -_lam5*_sth*_cth*_sb;
    MM2(1,0) = MM2(0,1);
    MM2(1,1) = _lam2*_sb*_sb+_lam5*_cth*_cth*_cb*_cb;
    MM2(1,2) = -_lam5*_sth*_cth*_cb;
    MM2(2,0) = MM2(0,2);
    MM2(2,1) = MM2(1,2);
    MM2(2,2) = _lam5*_sth*_sth;

    SelfAdjointEigenSolver<Matrix3d> eigensolver(MM2);
    Vector3d masses = eigensolver.eigenvalues();
    _MHL2 = masses(0)*_vev2;
    _MHH2 = masses(1)*_vev2;
    _MHA2 = masses(2)*_vev2;

    if (_MHL2 < 0 || _MHH2 < 0 || _MHA2 < 0)
    {
        return false;
    }
    _MHL = sqrt(_MHL2);
    _MHH = sqrt(_MHH2);
    _MHA = sqrt(_MHA2);

    _R = eigensolver.eigenvectors();

    _GetAlphas();
    return true;
}

bool THDM_SCPV1::Set_Physical_Parameters(double beta, double m2, double m3, double mpm, double alpha, double alphac)
{
    _vev = vev;
    _vev2 = _vev*vev;
    _beta = beta;
    _tb = tan(_beta);
    _sb = sin(_beta);
    _cb = cos(_beta);
    _v1 = _vev*_cb;
    _v2 = _vev*_sb;

    _MHL = MHL;
    _MHL2 = _MHL*_MHL;

    _MHH = m2;
    _MHH2 = _MHH*_MHH;

    _MHA = m3;
    _MHA2 = _MHA*_MHA;
    
    _MHpm = mpm;
    _MHpm2 = _MHpm*_MHpm;

    _alpha = alpha;
    _alphac = alphac;

    // * 1. Get alphab
    double sb = sin(_alphac)*cos(_alphac)*(_MHH2-_MHA2)/((_MHL2-pow(sin(_alphac),2)*_MHH2-pow(cos(_alphac),2)*_MHA2)*tan(_alpha+_beta));
    if (abs(sb)>1.0) return false;
    _alphab = asin(sb); // Just keep one such solution is enough. Other solution will be degenerate
    // cout<<"alphab:\t"<<_alphab<<endl;
    // * 2. Get mixing matrix
    _GetR();
    // cout<<"Got R"<<endl;
    // * 3. Get theta
    _theta = atan(-_sb*(_MHL2*_R(0,2)*_R(0,2)+_MHH2*_R(1,2)*_R(1,2)+_MHA2*_R(2,2)*_R(2,2))/(_MHL2*_R(0,0)*_R(0,2)+_MHH2*_R(1,0)*_R(1,2)+_MHA2*_R(2,0)*_R(2,2)));
    _tth = tan(_theta);
    _sth = sin(_theta);
    _cth = cos(_theta);
    _v2r = _v2*_cth;
    _v2i = _v2*_sth;
    // cout<<"theta:\t"<<_theta<<endl;
    // * 4. Calculate lambda's
    _lam1 = (_MHL2*(_R(0,0)*_R(0,0)-_sb*_sb/_tth/_tth*_R(0,2)*_R(0,2))+_MHH2*(_R(1,0)*_R(1,0)-_sb*_sb/_tth/_tth*_R(1,2)*_R(1,2))+_MHA2*(_R(2,0)*_R(2,0)-_sb*_sb/_tth/_tth*_R(2,2)*_R(2,2)))/_cb/_cb/_vev2;
    _lam2 = (_MHL2*(_R(0,1)*_R(0,1)-_cb*_cb/_tth/_tth*_R(0,2)*_R(0,2))+_MHH2*(_R(1,1)*_R(1,1)-_cb*_cb/_tth/_tth*_R(1,2)*_R(1,2))+_MHA2*(_R(2,1)*_R(2,1)-_cb*_cb/_tth/_tth*_R(2,2)*_R(2,2)))/_sb/_sb/_vev2;
    _lam3 = 2*_MHpm2/_vev2 + (_MHL2*(_R(0,0)*_R(0,1)-_sb*_cb/_tth/_tth*_R(0,2)*_R(0,2))+_MHH2*(_R(1,0)*_R(1,1)-_sb*_cb/_tth/_tth*_R(1,2)*_R(1,2))+_MHA2*(_R(2,0)*_R(2,1)-_sb*_cb/_tth/_tth*_R(2,2)*_R(2,2)))/_sb/_cb/_vev2;
    _lam4 = -2*_MHpm2/_vev2 + (_MHL2*(_R(0,2)*_R(0,2))+_MHH2*(_R(1,2)*_R(1,2))+_MHA2*(_R(2,2)*_R(2,2)))/_sth/_sth/_vev2;
    _lam5 = (_MHL2*(_R(0,2)*_R(0,2))+_MHH2*(_R(1,2)*_R(1,2))+_MHA2*(_R(2,2)*_R(2,2)))/_sth/_sth/_vev2;
    // * 5. Calculate m112 m222 m122
    _m112 = -_lam1*_cb*_cb*_vev2/2.0 - (_lam3+_lam4-_lam5)*_sb*_sb*_vev2/2.0;
    _m222 = -_lam2*_sb*_sb*_vev2/2.0 - (_lam3+_lam4-_lam5)*_cb*_cb*_vev2/2.0;
    _m122 = _lam5*_sb*_cb*_cth*_vev2;

    // cout<<_lam1<<"\t"<<_lam2<<"\t"<<_lam3<<"\t"<<_lam4<<"\t"<<_lam5<<endl;
    // cout<<_m112<<"\t"<<_m222<<"\t"<<_m122<<endl;

    _Solved = false;
    FindLocalMinima();
    _IndexInput = -1;
    for (size_t i = 0; i < _NLocalExtreme; i++)
    {
        if (CloseQ(_localExtreme[i],{_v1,_v2r,_v2i}))
        {
            _IndexInput = i;
        }
    }
    return true;
}

void THDM_SCPV1::_GetAlphas()
{
    double sab = _R(0,2);
    double cab = sqrt(1-sab*sab);

    double sac = _R(1,2)/cab;
    double cac = _R(2,2)/cab;
    double sa = -_R(0,0)/cab;
    double ca = _R(0,1)/cab;

    _alpha = GetAngle(sa,ca);
    _alphab = GetAngle(sab,cab);
    _alphac = GetAngle(sac,cac);
}

void THDM_SCPV1::_GetR()
{
    double ca = cos(_alpha);
    double sa = sin(_alpha);
    double cb = cos(_alphab);
    double sb = sin(_alphab);
    double cc = cos(_alphac);
    double sc = sin(_alphac);

    _R(0,0) = -sa*cb;
    _R(0,1) = ca*cb;
    _R(0,2) = sb;

    _R(1,0) = ca*cc + sa*sb*sc;
    _R(1,1) = sa*cc - ca*sb*sc;
    _R(1,2) = cb*sc;

    _R(2,0) = -ca*sc + sa*sb*cc;
    _R(2,1) = -sa*sc - ca*sb*cc;
    _R(2,2) = cb*cc;
}

void THDM_SCPV1::FindLocalMinima()
{
    if (_Solved)
    {
        return;
    }

    Clear_Local_Cache();
    // cout<<"Cache Cleared"<<endl;
    // * 1. v2i = 0;
    // * 1.1 v1 = 0, v2r = 0, v2i = 0
    _localExtreme.push_back({0,0,0});
    AppendLocalExtreme();
    // cout<<"Got First trivial solution"<<endl;

    // * 1.2 k = v2r/v1
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
        _localExtreme.push_back({sqrt(_v12temp),k*sqrt(_v12temp),0});
        AppendLocalExtreme();

        // ! temporarily ignore v1 < 0 solution, 
        // _localExtreme.push_back({-sqrt(_v12temp),-k*sqrt(_v12temp),0});
        // AppendLocalExtreme();
    }
    

    // * 2. v2i != 0
    // * 2.1 trivial case, v1 = 0, v2r = 0
    {
        double _v2i2temp = -2*_m222/_lam2;
        if (_v2i2temp >= 0)
        {
            _localExtreme.push_back({0,0,sqrt(_v2i2temp)});
            AppendLocalExtreme();

            _localExtreme.push_back({0,0,-sqrt(_v2i2temp)});
            AppendLocalExtreme();
        }
        
    }
    // * 2.2 non-trivial case, v1 != 0, v2r != 0
    double lam345p = _lam3+_lam4-_lam5;
    double _v12temp = (2*_m222*lam345p - 2*_m112*_lam2)/(_lam1*_lam2-lam345p*lam345p);
    if (_v12temp >= 0)
    {
        double _v2r2temp = _m122*_m122/_lam5/_lam5/_v12temp;
        double _v2i2temp = -(2*_m222+lam345p*_v12temp + _lam2*_v2r2temp)/_lam2;
        double _v1temp;
        double _v2rtemp;
        double _v2itemp;
        if (_v2i2temp >= 0)
        {
            _v1temp = sqrt(_v12temp);
            _v2rtemp = _m122/_lam5/_v1temp;
            _v2itemp = sqrt(_v2i2temp);
            _localExtreme.push_back({_v1temp,_v2rtemp,_v2itemp});
            AppendLocalExtreme();

            _v1temp = sqrt(_v12temp);
            _v2rtemp = _m122/_lam5/_v1temp;
            _v2itemp = -sqrt(_v2i2temp);
            _localExtreme.push_back({_v1temp,_v2rtemp,_v2itemp});
            AppendLocalExtreme();

            // ! temporarily ignore v1 < 0 solution, 
            // _v1temp = -sqrt(_v12temp);
            // _v2rtemp = _m122/_lam5/_v1temp;
            // _v2itemp = sqrt(_v2i2temp);
            // _localExtreme.push_back({_v1temp,_v2rtemp,_v2itemp});
            // AppendLocalExtreme();

            // _v1temp = -sqrt(_v12temp);
            // _v2rtemp = _m122/_lam5/_v1temp;
            // _v2itemp = -sqrt(_v2i2temp);
            // _localExtreme.push_back({_v1temp,_v2rtemp,_v2itemp});
            // AppendLocalExtreme();
        } 
    }

    _Solved = true;
}

void THDM_SCPV1::PrintParameters()
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
    cout<<"theta:\t"<<_theta<<endl;
    cout<<"v1:\t"<<_v1<<endl;
    cout<<"v2:\t"<<_v2<<endl;
    cout<<"v2r:\t"<<_v2r<<endl;
    cout<<"v2i:\t"<<_v2i<<endl;
    cout<<"m1:\t"<<_MHL<<endl;
    cout<<"m2:\t"<<_MHH<<endl;
    cout<<"m3:\t"<<_MHA<<endl;
    cout<<"m+-:\t"<<_MHpm<<endl;
    cout<<"alpha:\t"<<_alpha<<endl;
    cout<<"alphab:\t"<<_alphab<<endl;
    cout<<"alphac:\t"<<_alphac<<endl;
}

double THDM_SCPV1::Vtotal(VD field_values, double scale)
{
    double phi1 = field_values[0];
    double phi2R = field_values[1];
    double phi2I = field_values[2];
    double phi22 = phi2R*phi2R + phi2I*phi2I;
    double tmp = 0;
    tmp += 4*_m112*phi1*phi1;
    tmp += -8*_m122*phi1*phi2R;
    tmp += 4*_m222*phi22;
    tmp += _lam1*pow(phi1,4);
    tmp += -2*_lam5*phi1*phi1*phi2I*phi2I;
    tmp += 2*_lam3*phi1*phi1*phi22;
    tmp += 2*_lam4*phi1*phi1*phi22;
    tmp += 2*_lam5*phi1*phi1*phi2R*phi2R;
    tmp += _lam2*phi22*phi22;
    return tmp/8.0/pow(scale,4);
}
VD THDM_SCPV1::dVtotal(VD field_values, double scale)
{
    double phi1 = field_values[0];
    double phi2R = field_values[1];
    double phi2I = field_values[2];
    double phi22 = phi2R*phi2R + phi2I*phi2I;
    VD res(_N_VEVs);
    res[0] = (2*_m112*phi1-2*_m122*phi2R+phi1*(_lam1*phi1*phi1 + (_lam3+_lam4-_lam5)*phi2I*phi2I + (_lam3+_lam4+_lam5)*phi2R*phi2R))/2;
    res[1] = (-2*_m122*phi1+phi2R*(2*_m222+phi1*phi1*(_lam3+_lam4+_lam5)+_lam2*phi22))/2;
    res[2] = phi2I*(2*_m222+phi1*phi1*(_lam3+_lam4-_lam5)+_lam2*phi22)/2;
    return res/pow(scale,3);
}
VVD THDM_SCPV1::d2Vtotal(VD field_values, double scale)
{
    double phi1 = field_values[0];
    double phi2R = field_values[1];
    double phi2I = field_values[2];
    double phi22 = phi2R*phi2R + phi2I*phi2I;
    VVD res(_N_VEVs,VD(_N_VEVs));
    res[0][0] = (3*_lam1*phi1*phi1+2*_m112+(_lam3+_lam4-_lam5)*phi2I*phi2I+(_lam3+_lam4+_lam5)*phi2R*phi2R)/2;
    res[0][1] = -_m122+(_lam3+_lam4+_lam5)*phi1*phi2R;
    res[0][2] = (_lam3+_lam4-_lam5)*phi1*phi2I;
    res[1][0] = res[0][1];
    res[1][1] = (2*_m222+(_lam3+_lam4+_lam5)*phi1*phi1+_lam2*phi2I*phi2I+3*_lam2*phi2R*phi2R)/2;
    res[1][2] = _lam2*phi2I*phi2R;
    res[2][0] = res[0][2];
    res[2][1] = res[1][2];
    res[2][2] = (2*_m222+(_lam3+_lam4+_lam5)*phi1*phi1+_lam2*phi2R*phi2R+3*_lam2*phi2I*phi2I)/2;

    return res/pow(scale,2);
}
double THDM_SCPV1::V0_global(double scale)
{
    return Vtotal({_v1,_v2r,_v2i},scale);
}
VD THDM_SCPV1::QuarticCoupling(VD field_values)
{
    return {_lam1/8,_lam2/8,_lam2/8};
}
bool THDM_SCPV1::CheckUnitarity(double MAX)
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
bool THDM_SCPV1::CheckStability()
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

string THDM_SCPV1::repr()
{
    return "vev\tbeta\ttheta\talpha\talphab\talphac\tm1\tm2\tm3\tmpm\tm112\tm222\tm122\tlam1\tlam2\tlam3\tlam4\tlam5";
}
ostream &operator<<(std::ostream &os, const THDM_SCPV1 &mod)
{
    os<<mod._vev<<"\t"<<mod._beta<<"\t"<<mod._theta<<"\t"<<mod._alpha<<"\t"<<mod._alphab<<"\t"<<mod._alphac<<"\t"<<mod._MHL<<"\t"<<mod._MHH<<"\t"<<mod._MHA<<"\t"<<mod._MHpm<<"\t"<<mod._m112<<"\t"<<mod._m222<<"\t"<<mod._m122<<"\t"<<mod._lam1<<"\t"<<mod._lam2<<"\t"<<mod._lam3<<"\t"<<mod._lam4<<"\t"<<mod._lam5;
    return os;
}