#include "cxSM_CP.h"

using namespace Eigen;

cxSM_CP::cxSM_CP():Basic_Model(3)
{

}

void cxSM_CP::Set_Potential_Parameters(double mu2, double lam, double del2, double b2, double d2, double del3, double b1, double d1, double d3)
{
    _mu2 = mu2;
    _lam = lam;
    _del2 = del2;
    _b2 = b2;
    _d2 = d2;
    _del3 = del3;
    _b1 = b1;
    _d1 = d1;
    _d3 = d3;
    _Solved = false;

    FindLocalMinima();

    _vev = 0;
    _vs = 0;
    _alpha = 0;
    _IndexInput = -1;

    for (size_t i = 0; i < _NLocalExtreme; i++)
    {
        if (_LocalMinimaQ[i] && _localExtreme[i][0] > _vev && _localExtreme[i][1] > _vs)
        {
            _vev = _localExtreme[i][0];
            _vs = _localExtreme[i][1];
            _alpha = _localExtreme[i][2];
            _IndexInput = i;
        }
    }
    Matrix3d MM2;
    MM2(0,0) = 2*_lam*_vev*_vev;
    MM2(0,1) = (_del2+_del3)*cos(_alpha)*_vev*_vs/2;
    MM2(0,2) = (_del2-_del3)*sin(_alpha)*_vev*_vs/2;
    MM2(1,0) = MM2(0,1);
    MM2(1,1) = (_d1+_d2+_d3)*pow(cos(_alpha)*_vs,2)/2;
    MM2(1,2) = -(3*_d1-_d2)*sin(2*_alpha)*_vs*_vs/4;
    MM2(2,0) = MM2(0,2);
    MM2(2,1) = MM2(1,2);
    MM2(2,2) = (_d1+_d2-_d3)*pow(sin(_alpha)*_vs,2)/2;

    SelfAdjointEigenSolver<Matrix3d> eigensolver(MM2);
    Vector3d masses = eigensolver.eigenvalues();
    _MHL2 = masses(0);
    _MHH2 = masses(1);
    _MHA2 = masses(2);

    if (_MHL2 < 0 || _MHH2 < 0 || _MHA2 < 0)
    {
        return;
    }
    _MHL = sqrt(_MHL2);
    _MHH = sqrt(_MHH2);
    _MHA = sqrt(_MHA2);

    _R = eigensolver.eigenvectors();

}


void cxSM_CP::FindLocalMinima()
{
    if (_Solved)
    {
        return;
    }

    Clear_Local_Cache();

    // ! v = 0, vsr = 0, vsi = 0
    _localExtreme.push_back({0,0,0});
    AppendLocalExtreme();

    // ! v2 = -mu2/lam, vsr = 0, vsi = 0
    if (-_mu2/_lam >= 0)
    {
        _localExtreme.push_back({sqrt(-_mu2/_lam),0,0});
        AppendLocalExtreme();
    }

    // ! v = 0, vsr2 = - 2(b1+b2)/(d1+d2+d3), vsi = 0
    if (-2*(_b1+_b2)/(_d1+_d2+_d3) >= 0 )
    {
        _localExtreme.push_back({0,sqrt(-2*(_b1+_b2)/(_d1+_d2+_d3)),0});
        AppendLocalExtreme();

        _localExtreme.push_back({0,sqrt(-2*(_b1+_b2)/(_d1+_d2+_d3)),M_PI});
        AppendLocalExtreme();
    }

    // ! v = 0, vsr = 0, vsi2 = -2(b2-b1)/(d1+d2-d3)
    if (-2*(_b2-_b1)/(_d1+_d2-_d3))
    {
        _localExtreme.push_back({0,sqrt(-2*(_b2-_b1)/(_d1+_d2-_d3)),M_PI_2});
        AppendLocalExtreme();

        _localExtreme.push_back({0,sqrt(-2*(_b2-_b1)/(_d1+_d2-_d3)),-M_PI_2});
        AppendLocalExtreme();
    }

    // ! v = 0, vsr2 != 0 vsi2 != 0
    double tmpa1 = 2*(_b2+_b1);
    double tmpa2 = 2*(_b2-_b1);
    double tmpc1 = _d1+_d2+_d3;
    double tmpc2 = _d2-3*_d1;
    double tmpd1 = _d2-3*_d1;
    double tmpd2 = _d1+_d2-_d3;
    double _vsr2 = -(tmpa2*tmpd1-tmpa1*tmpd2)/(tmpc2*tmpd1-tmpc1*tmpd2);
    double _vsi2 = -(tmpa1*tmpc2-tmpa2*tmpc1)/(tmpc2*tmpd1-tmpc1*tmpd2);
    double vsamp;
    double vsarg;
    if (_vsr2 >=0 && _vsi2 >= 0)
    {
        _vsr = sqrt(_vsr2);
        _vsi = sqrt(_vsi2);
        vsamp=_mag(_vsr,_vsi);
        vsarg=_arg(_vsr,_vsi);
        _localExtreme.push_back({0,vsamp,vsarg});
        AppendLocalExtreme();

        _localExtreme.push_back({0,vsamp,-vsarg});
        AppendLocalExtreme();

        _localExtreme.push_back({0,vsamp,M_PI-vsarg});
        AppendLocalExtreme();

        _localExtreme.push_back({0,vsamp,-M_PI+vsarg});
        AppendLocalExtreme();

    }

    // ! vsi = 0, v2 != 0, vsr2 != 0
    tmpa1 = 4*_mu2;
    tmpa2 = 2*(_b1+_b2);
    tmpc1 = 4*_lam;
    tmpc2 = _del2 + _del3;
    tmpd1 = _del2 - _del3;
    tmpd2 = _d1 + _d2 + _d3;

    double _vev2 = -(tmpa2*tmpd1-tmpa1*tmpd2)/(tmpc2*tmpd1-tmpc1*tmpd2);
    _vsr2 = -(tmpa1*tmpc2-tmpa2*tmpc1)/(tmpc2*tmpd1-tmpc1*tmpd2);
    
    if (_vev2 >= 0 && _vsr2 >= 0)
    {
        _vev = sqrt(_vev2);
        _vsr = sqrt(_vsr2);
        _vsi = 0;
        
        _localExtreme.push_back({_vev,_vsr,0});
        AppendLocalExtreme();

        _localExtreme.push_back({_vev,_vsr,M_PI});
        AppendLocalExtreme();

    }

    // ! vsr = 0, v2 != 0, vsi2 != 0
    tmpa1 = 4*_mu2;
    tmpa2 = 2*(_b2-_b1);
    tmpc1 = 4*_lam;
    tmpc2 = _del2 - _del3;
    tmpd1 = _del2 + _del3;
    tmpd2 = _d1 + _d2 - _d3;

    _vev2 = -(tmpa2*tmpd1-tmpa1*tmpd2)/(tmpc2*tmpd1-tmpc1*tmpd2);
    _vsi2 = -(tmpa1*tmpc2-tmpa2*tmpc1)/(tmpc2*tmpd1-tmpc1*tmpd2);

    if (_vev2 >=0 && _vsi2 >=0 )
    {
        _vev = sqrt(_vev2);
        _vsr = 0;
        _vsi = sqrt(_vsi2);

        _localExtreme.push_back({_vev,_vsi,M_PI_2});
        AppendLocalExtreme();

        _localExtreme.push_back({_vev,_vsi,-M_PI_2});
        AppendLocalExtreme();
    }

    // ! v!=0, vsr!=0, vsi!=0
    tmpa1 = 4*_mu2;
    tmpa2 = 2*(_b2+_b1);
    double tmpa3 = 2*(_b2-_b1);
    double tmpb1 = 4*_lam;
    double tmpb2 = _del2 + _del3;
    double tmpb3 = _del2 - _del3;
    tmpc1 = _del2 + _del3;
    tmpc2 = _d1 + _d2 + _d3;
    double tmpc3 = _d2 - 3*_d1;
    tmpd1 = _del2 - _del3;
    tmpd2 = _d2 - 3*_d1;
    double tmpd3 = _d1 + _d2 - _d3;

    _vev2 = -(-tmpa1*tmpc2*tmpd3+tmpa1*tmpc3*tmpd2+tmpa2*tmpc1*tmpd3-tmpa2*tmpc3*tmpd1-tmpa3*tmpc1*tmpd2+tmpa3*tmpc2*tmpd1)/(-tmpb1*tmpc2*tmpd3+tmpb1*tmpc3*tmpd2+tmpb2*tmpc1*tmpd3-tmpb2*tmpc3*tmpd1-tmpb3*tmpc1*tmpd2+tmpb3*tmpc2*tmpd1);
    _vsr2 = -(tmpa1*tmpb2*tmpd3-tmpa1*tmpb3*tmpd2-tmpa2*tmpb1*tmpd3+tmpa2*tmpb3*tmpd1+tmpa3*tmpb1*tmpd2-tmpa3*tmpb2*tmpd1)/(-tmpb1*tmpc2*tmpd3+tmpb1*tmpc3*tmpd2+tmpb2*tmpc1*tmpd3-tmpb2*tmpc3*tmpd1-tmpb3*tmpc1*tmpd2+tmpb3*tmpc2*tmpd1);
    _vsi2 = -(tmpa1*tmpb2*tmpc3-tmpa1*tmpb3*tmpc2-tmpa2*tmpb1*tmpc3+tmpa2*tmpb3*tmpc1+tmpa3*tmpb1*tmpc2-tmpa3*tmpb2*tmpc1)/(tmpb1*tmpc2*tmpd3-tmpb1*tmpc3*tmpd2-tmpb2*tmpc1*tmpd3+tmpb2*tmpc3*tmpd1+tmpb3*tmpc1*tmpd2-tmpb3*tmpc2*tmpd1);
    
    if (_vev2 >= 0 && _vsr2 >= 0 && _vsi2 >= 0)
    {
        _vev = sqrt(_vev2);
        _vsr = sqrt(_vsr2);
        _vsi = sqrt(_vsi2);

        vsamp = _mag(_vsr,_vsi);
        vsarg = _arg(_vsr,_vsi);

        _localExtreme.push_back({_vev,vsamp,vsarg});
        AppendLocalExtreme();

        _localExtreme.push_back({_vev,vsamp,-vsarg});
        AppendLocalExtreme();

        _localExtreme.push_back({_vev,vsamp,M_PI-vsarg});
        AppendLocalExtreme();

        _localExtreme.push_back({_vev,vsamp,-M_PI+vsarg});
        AppendLocalExtreme();
  
    }
    
    _Solved = true;    
    
    
}
double cxSM_CP::_mag(double vsr, double vsi)
{
    return sqrt(vsr*vsr+vsi*vsi);
}
double cxSM_CP::_arg(double vsr, double vsi)
{
    double amg = _mag(vsr,vsi);
    if (amg == 0)
    {
        return 0;
    }
    double sth = vsi/amg;
    double cth = vsr/amg;

    double angle = asin(sth);
    if (cth < 0)
    {
        if (sth >= 0)
        {
            angle = M_PI - angle;
        }
        else
        {
            angle = -M_PI - angle;
        }
    }
    return angle;
}