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
            _vsr = _localExtreme[i][1];
            _vsi = _localExtreme[i][2];
            _IndexInput = i;
        }
    }
    Matrix3d MM2;
    MM2(0,0) = 2*_lam*_vev*_vev;
    MM2(0,1) = (_del2+_del3)*_vev*_vsr/2;
    MM2(0,2) = (_del2-_del3)*_vev*_vsi/2;
    MM2(1,0) = MM2(0,1);
    MM2(1,1) = (_d1+_d2+_d3)*pow(_vsr,2)/2;
    MM2(1,2) = (_d2-3*_d1)*_vsi*_vsr/2;
    MM2(2,0) = MM2(0,2);
    MM2(2,1) = MM2(1,2);
    MM2(2,2) = (_d1+_d2-_d3)*pow(_vsi,2)/2;

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

    _GetThetas();

}

void cxSM_CP::Set_Physical_Parameters(double vsr, double vsi, double MHH, double MHA, double theta1, double theta2, double theta3)
{
    _vev = vev;
    _vsr = vsr;
    _vsi = vsi;

    _MHL = MHL;
    _MHH = MHH;
    _MHA = MHA;

    _MHL2 = MHL2;
    _MHH2 = MHH*MHH;
    _MHA2 = MHA*MHA;

    _theta1 = theta1;
    _theta2 = theta2;
    _theta3 = theta3;

    _GetR();

    _lam  = (_MHL2*pow(_R(0,0),2)+_MHH2*pow(_R(0,1),2)+_MHA2*pow(_R(0,2),2))/2/_vev/_vev;
    _del2 = (_MHA2*_R(0,2)*_R(1,2)*_vsi + _MHA2*_R(0,2)*_R(2,2)*_vsr + _MHH2*_R(0,1)*_R(1,1)*_vsi + _MHH2*_R(0,1)*_R(2,1)*_vsr + _MHL2*_R(0,0)*_R(1,0)*_vsi + _MHL2*_R(0,0)*_R(2,0)*_vsr)/_vev/_vsr/_vsi;
    _d2   = (3*_vsi*_vsi*(_MHA2*pow(_R(1,2),2)+_MHH2*pow(_R(1,1),2+_MHL2*pow(_R(1,0),2)))+2*_vsi*_vsr*(_MHA2*_R(1,2)*_R(2,2)+_MHH2*_R(1,1)*_R(2,1)+_MHL2*_R(1,0)*_R(2,0))+3*_vsr*_vsr*(_MHA2*pow(_R(2,2),2)+_MHH2*pow(_R(2,1),2)+_MHL2*pow(_R(2,0),2)))/4/_vsi/_vsi/_vsr/_vsr;
    _del3 = (_vsi*(_MHA2*_R(0,2)*_R(1,2)+_MHH2*_R(0,1)*_R(1,1)+_MHL2*_R(0,0)*_R(1,0))-_vsr*(_MHA2*_R(0,2)*_R(2,2)+_MHH2*_R(0,1)*_R(2,1)+_MHL2*_R(0,0)*_R(2,0)))/_vev/_vsi/_vsr;
    _d1   = (_MHA2*pow(_R(1,2)*_vsi-_R(2,2)*_vsr,2)+_MHH2*pow(_R(1,1)*_vsi-_R(2,1)*_vsr,2)+_MHL2*pow(_R(1,0)*_vsi-_R(2,0)*_vsr,2))/4/_vsr/_vsr/_vsi/_vsi;
    _d3   = (_MHA2*pow(_R(1,2),2)+_MHH2*pow(_R(1,1),2)+_MHL2*pow(_R(1,0),2))/_vsr/_vsr - (_MHA2*pow(_R(2,2),2)+_MHH2*pow(_R(2,1),2)+_MHL2*pow(_R(2,0),2))/_vsi/_vsi;

    _mu2 = (-4*_lam*_vev*_vev + (_del3-_del2)*_vsi*_vsi - (_del2 + _del3)*_vsr*_vsr)/4;
    _b1  = (4*_d1*(_vsi*_vsi - _vsr*_vsr) - _d3*(_vsi*_vsi+_vsr*_vsr)-2*_del3*_vev*_vev)/4;
    _b2  = (2*_d1*(_vsi*_vsi+_vsr*_vsr)-2*_d2*(_vsi*_vsi+_vsr*_vsr)+_d3*(_vsi*_vsi-_vsr*_vsr)-2*_del2*_vev*_vev)/4;

    _Solved = false;
    FindLocalMinima();
    _IndexInput = -1;
    for (size_t i = 0; i < _NLocalExtreme; i++)
    {
        if (CloseQ(_localExtreme[i],{_vev,_vsr,_vsi}))
        {
            _IndexInput = i;
        }
    }
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

        _localExtreme.push_back({0,-sqrt(-2*(_b1+_b2)/(_d1+_d2+_d3)),0});
        AppendLocalExtreme();
    }

    // ! v = 0, vsr = 0, vsi2 = -2(b2-b1)/(d1+d2-d3)
    if (-2*(_b2-_b1)/(_d1+_d2-_d3))
    {
        _localExtreme.push_back({0,0,sqrt(-2*(_b2-_b1)/(_d1+_d2-_d3))});
        AppendLocalExtreme();

        _localExtreme.push_back({0,0,-sqrt(-2*(_b2-_b1)/(_d1+_d2-_d3))});
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
    double _vsrtmp;
    double _vsitmp;
    double _vevtmp;
    // double vsamp;
    // double vsarg;
    if (_vsr2 >=0 && _vsi2 >= 0)
    {
        _vsrtmp = sqrt(_vsr2);
        _vsitmp = sqrt(_vsi2);
        // vsamp=_mag(_vsr,_vsi);
        // vsarg=_arg(_vsr,_vsi);
        _localExtreme.push_back({0,_vsrtmp,_vsitmp});
        AppendLocalExtreme();

        _localExtreme.push_back({0,_vsrtmp,-_vsitmp});
        AppendLocalExtreme();

        _localExtreme.push_back({0,-_vsrtmp,_vsitmp});
        AppendLocalExtreme();

        _localExtreme.push_back({0,-_vsrtmp,-_vsitmp});
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
        _vevtmp = sqrt(_vev2);
        _vsrtmp = sqrt(_vsr2);
        _vsitmp = 0;
        
        _localExtreme.push_back({_vevtmp,_vsrtmp,0});
        AppendLocalExtreme();

        _localExtreme.push_back({_vevtmp,-_vsrtmp,0});
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
        _vevtmp = sqrt(_vev2);
        _vsrtmp = 0;
        _vsitmp = sqrt(_vsi2);

        _localExtreme.push_back({_vevtmp,0,_vsitmp});
        AppendLocalExtreme();

        _localExtreme.push_back({_vevtmp,0,-_vsitmp});
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
        _vevtmp = sqrt(_vev2);
        _vsrtmp = sqrt(_vsr2);
        _vsitmp = sqrt(_vsi2);

        // vsamp = _mag(_vsr,_vsi);
        // vsarg = _arg(_vsr,_vsi);

        _localExtreme.push_back({_vevtmp,_vsrtmp,_vsitmp});
        AppendLocalExtreme();

        _localExtreme.push_back({_vevtmp,_vsrtmp,-_vsitmp});
        AppendLocalExtreme();

        _localExtreme.push_back({_vevtmp,-_vsrtmp,_vsitmp});
        AppendLocalExtreme();

        _localExtreme.push_back({_vevtmp,-_vsrtmp,-_vsitmp});
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

    // double angle = asin(sth);
    // if (cth < 0)
    // {
    //     if (sth >= 0)
    //     {
    //         angle = M_PI - angle;
    //     }
    //     else
    //     {
    //         angle = -M_PI - angle;
    //     }
    // }
    // return angle;
    return GetAngle(sth,cth);
}
void cxSM_CP::_GetThetas()
{
    double s2 = _R(0,2);
    double c2 = sqrt(1-s2*s2);

    double c1 = _R(0,0)/c2;
    double s1 = _R(0,1)/c2;
    double s3 = _R(1,2)/c2;
    double c3 = _R(2,2)/c2;

    _theta1 = GetAngle(s1,c1);
    _theta2 = GetAngle(s2,c2);
    _theta3 = GetAngle(s3,c3);
}
void cxSM_CP::_GetR()
{
    double c1 = cos(_theta1);
    double s1 = sin(_theta1);
    double c2 = cos(_theta2);
    double s2 = sin(_theta2);
    double c3 = cos(_theta3);
    double s3 = sin(_theta3);

    _R(0,0) = c1*c2;
    _R(0,1) = s1*c2;
    _R(0,2) = s2;
    _R(1,0) = -(c1*s2*s3+s1*c3);
    _R(1,1) = c1*c3-s1*s2*s3;
    _R(1,2) = c2*s3;
    _R(2,0) = -c1*s2*c3+s1*s3;
    _R(2,1) = -(c1*s2+s1*s2*c3);
    _R(2,2) = c2*c3;
}
double cxSM_CP::Vtotal(VD field_values, double scale = 1)
{
    double vh = field_values[0];
    double vsr = field_values[1];
    double vsi = field_values[2];
    double vh2 = vh*vh;
    double vsr2 = vsr*vsr;
    double vsi2 = vsi*vsi;

    return (-_b1*vsi2/4+_b1*vsr2/4+_b2*vsi2/4+_b2*vsr2/4+_d1*vsi2*vsi2/16-3*_d1*vsi2*vsr2/8+_d1*vsr2*vsr2/16+_d2*vsi2*vsi2/16+_d2*vsi2*vsr2/8+_d2*vsr2*vsr2/16-_d3*vsi2*vsi2/16+_d3*vsr2*vsr2/16+_lam*vh2*vh2/4+_mu2*vh2/2+_del2*vh2*vsi2/8-_del3*vh2*vsi2/8+_del2*vh2*vsr2/8+_del3*vh2*vsr2/8)/pow(scale,4);
}
VD cxSM_CP::dVtotal(VD field_values, double scale = 1)
{
    double vh = field_values[0];
    double vsr = field_values[1];
    double vsi = field_values[2];
    double vh2 = vh*vh;
    double vsr2 = vsr*vsr;
    double vsi2 = vsi*vsi;

    VD res(_N_VEVs);
    res[0] = vh*(4*(_mu2+_lam*vh2)+vsi2*(_del2-_del3)+vsr2*(_del2+_del3))/4;
    res[1] = vsr*(2*_b1+2*_b2+vsr2*(_d1+_d2+_d3)+vsi2*(_d2-3*_d1)+vh2*(_del2+_del3))/4;
    res[2] = vsi*(-2*_b1+2*_b2+vsi2*(_d1+_d2-_d3)+vsr2*(_d2-3*_d1)+vh2*(_del2-_del3))/4;

    return res/pow(scale,3);
}
VVD cxSM_CP::d2Vtotal(VD field_values, double scale = 1)
{
    double vh = field_values[0];
    double vsr = field_values[1];
    double vsi = field_values[2];
    double vh2 = vh*vh;
    double vsr2 = vsr*vsr;
    double vsi2 = vsi*vsi;

    VVD res(_N_VEVs,VD(_N_VEVs));
    res[0][0] = ((_del2-_del3)*vsi2+(_del2+_del3)*vsr2+4*(_mu2+3*_lam*vh2))/4;
    res[0][1] = (_del2+_del3)*vh*vsr/2;
    res[0][2] = (_del2-_del3)*vh*vsi/2;
    res[1][0] = res[0][1];
    res[1][1] = ((_del2+_del3)*vh2+(_d2-3*_d1)*vsi2+3*(_d1+_d2+_d3)*vsr2+2*_b1+2*_b2)/4;
    res[1][2] = (_d2-3*_d1)*vsi*vsr/2;
    res[2][0] = res[0][2];
    res[2][1] = res[1][2];
    res[2][2] = ((_del2-_del3)*vh2+3*(_d1+_d2-_d3)*vsi2+(_d2-3*_d1)*vsr2-2*_b1+2*_b2)/4;

    return res/pow(scale,2);
}
double cxSM_CP::V0_global(double scale = 1)
{
    return Vtotal({_vev,_vsr,_vsi},scale);
}