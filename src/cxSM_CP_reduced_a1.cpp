#include "cxSM_CP_reduced_a1.h"
#include "CubicSolver.h"
// #include "QuarticSolver.h"
#include <functional>
#include <iostream>

using namespace std;
using namespace Eigen;

cxSM_CP_reduced_a1::cxSM_CP_reduced_a1():cxSM_CP()
{

}

void cxSM_CP_reduced_a1::Set_Potential_Parameters(double mu2, double lam, double del2, double b2, double d2, double a1, double b1)
{
    cxSM_CP::Set_Potential_Parameters(mu2,lam,del2,b2,d2,0.0,b1,0.0,0.0);
    _a1 = a1;

    _Solved = false;

    FindLocalMinima();

    _vev = 0;
    _vsr = 0;
    _alpha = 0;
    _IndexInput = -1;

    for (size_t i = 0; i < _NLocalExtreme; i++)
    {
        if (_LocalMinimaQ[i] && _localExtreme[i][0] > _vev && _localExtreme[i][1] > _vsr)
        {
            _vev = _localExtreme[i][0];
            _vsr = _localExtreme[i][1];
            _vsi = _localExtreme[i][2];
            _IndexInput = i;
        }
    }
    _vs = GetMag(_vsr,_vsi);
    _alpha = GetArg(_vsr,_vsi);

    Matrix3d MM2;
    MM2(0,0) = 2*_lam*_vev*_vev;
    MM2(0,1) = (_del2)*_vev*_vsr/2;
    MM2(0,2) = (_del2)*_vev*_vsi/2;
    MM2(1,0) = MM2(0,1);
    MM2(1,1) = (_d2)*pow(_vsr,2)/2 + _b1;
    MM2(1,2) = (_d2)*_vsi*_vsr/2;
    MM2(2,0) = MM2(0,2);
    MM2(2,1) = MM2(1,2);
    MM2(2,2) = (_d2)*pow(_vsi,2)/2;

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

bool cxSM_CP_reduced_a1::Set_Physical_Parameters_vs_theta(double vs, double MHH, double MHA, double theta1, double theta3)
{
    double alpha,theta2;
    bool good = _GetAlphaTheta2(MHH,MHA,theta1,theta3,alpha,theta2);
    if (!good)
    {
        return false;
    }
    double vsr = vs*cos(alpha);
    double vsi = vs*sin(alpha);
    Set_Physical_Parameters(vsr,vsi,MHH,MHA,theta1,theta2,theta3);
    return true;
}
bool cxSM_CP_reduced_a1::Set_Physical_Parameters_vsr_theta(double vsr, double MHH, double MHA, double theta1, double theta3)
{
    double alpha,theta2;
    bool good = _GetAlphaTheta2(MHH,MHA,theta1,theta3,alpha,theta2);
    if (!good)
    {
        return false;
    }
    double vsi = tan(alpha)*vsr;
    Set_Physical_Parameters(vsr,vsi,MHH,MHA,theta1,theta2,theta3);
    return true;
}
bool cxSM_CP_reduced_a1::_GetAlphaTheta2(const double MHH, const double MHA, const double theta1, double theta3, double &alpha, double &theta2)
{
    // In this situation, we always have solution
    double m12 = MHL2;
    double m22 = MHH*MHH;
    double m32 = MHA*MHA;

    // * First solution
    theta2 = M_PI_2; // * Only consider positive one, negative one should be equivalent (! not varified)
    alpha = atan(((m12-m22)*cos(2*theta1+2*theta3)+m12+m22)/((m12-m22)*sin(2*theta1+2*theta3)));

    // * The other solution
    double sth2 = 2*m32*sin(theta1)*cos(theta1)*(m22-m12)*cos(theta3)/sin(theta3)/(m32*cos(2*theta1)*(m12-m22)-m32*(m12+m22)+2*m12*m22);
    if (abs(sth2)>1)
    {
        return true;
    }
    theta2 = asin(sth2);
    double s1 = sin(theta1);
    double s2 = sin(theta2);
    double s3 = sin(theta3);
    double c1 = cos(theta1);
    double c2 = cos(theta2);
    double c3 = cos(theta3);

    double ta = (m12*pow(c1*s2*c3-s1*s3,2)+m22*pow(s1*s2*c3+c1*s3,2)+m32*pow(c2*c3,2))/(c1*c1*s3*c3*(m12*s2*s2-m22)+s1*c1*s2*(c3*c3-s3*s3)*(m12-m22)+s3*c3*(s1*s1*(m22*s2*s2-m12)+m32*c2*c2));
    alpha = atan(ta);
    return true;
}

bool cxSM_CP_reduced_a1::Set_Physical_Parameters(double vsr, double vsi, double MHH, double MHA, double theta1, double theta2, double theta3)
{
    _vev = vev;
    _vsr = vsr;
    _vsi = vsi;

    _vs = GetMag(_vsr,_vsi);
    _alpha = GetArg(_vsr,_vsi);

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

    _lam = (_MHL2*pow(_R(0,0),2)+_MHH2*pow(_R(0,1),2)+_MHA2*pow(_R(0,2),2))/2/_vev/_vev;
    _del2 = 2*(_MHL2*_R(0,0)*_R(1,0)+_MHH2*_R(0,1)*_R(1,1)+_MHA2*_R(0,2)*_R(1,2))/_vev/_vsr;
    _d2 = 2*(_MHL2*pow(_R(2,0),2)+_MHH2*pow(_R(2,1),2)+_MHA2*pow(_R(2,2),2))/_vsi/_vsi;
    _b1 = _MHL2*(pow(_R(1,0),2)-pow(_R(2,0)*_vsr/_vsi,2))+_MHH2*(pow(_R(1,1),2)-pow(_R(2,1)*_vsr/_vsi,2))+_MHA2*(pow(_R(1,2),2)-pow(_R(2,2)*_vsr/_vsi,2));

    _mu2 = -(4*_lam*_vev*_vev+_del2*_vsi*_vsi+_del2*_vsr*_vsr)/4;
    _a1 = -_b1*_vsr/sqrt(2);
    _b2 = _b1 - _d2*(_vsi*_vsi+_vsr*_vsr)/2 - _del2*_vev*_vev/2;

    _del3 = 0;
    _d1 = 0;
    _d3 = 0;

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
void cxSM_CP_reduced_a1::_GetR()
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
    _R(2,1) = -(c1*s3+s1*s2*c3);
    _R(2,2) = c2*c3;
}
void cxSM_CP_reduced_a1::_GetThetas()
{
    // Please note that following treatment is not correct actually.
    // I keep it temporarily, I need to update the algorithm later.
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
void cxSM_CP_reduced_a1::FindLocalMinima()
{
    if (_Solved)
    {
        return;
    }

    Clear_Local_Cache();
    
    CubicSolver _solver;
    double a0,a1,a2,a3;
    // ! v = 0, vsr != 0, vsi = 0
    a0 = 4*sqrt(2)*_a1;
    a1 = 2*(_b1+_b2);
    a2 = 0;
    a3 = _d2;
    _solver.Solve(a3,a2,a1,a0);
    VD _sol = _solver.GetRealSolution();
    for (size_t i = 0; i < _sol.size(); i++)
    {
        _localExtreme.push_back({0,_sol[i],0});
        AppendLocalExtreme();
    }

    // ! v=0, vsr != 0, vsi != 0
    double vsrtmp = -sqrt(2)*_a1/_b1;
    double vs2tmp = 2*(_b1-_b2)/_d2;
    double vsi2tmp = vs2tmp - pow(vsrtmp,2);
    if (vsi2tmp >= 0)
    {
        _localExtreme.push_back({0,vsrtmp,sqrt(vsi2tmp)});
        AppendLocalExtreme();

        _localExtreme.push_back({0,vsrtmp,-sqrt(vsi2tmp)});
        AppendLocalExtreme();
    }

    // ! v != 0, vsr != 0, vsi = 0
    a0 = 4*sqrt(2)*_a1;
    a1 = 2*(_b1+_b2)-_del2*_mu2/_lam;
    a2 = 0;
    a3 = _d2 - _del2*_del2/4/_lam;
    _solver.Solve(a3,a2,a1,a0);
    _sol = _solver.GetRealSolution();
    double vh2tmp;
    for (size_t i = 0; i < _sol.size(); i++)
    {
        vh2tmp = -(4*_mu2+_del2*_sol[i]*_sol[i])/4/_lam;
        if (vh2tmp >= 0)
        {
            _localExtreme.push_back({sqrt(vh2tmp),_sol[i],0}); // We don't need to consider vh < 0;
            AppendLocalExtreme();
        }
    }

    // ! v != 0, vsr != 0, vsi != 0
    vsrtmp = -sqrt(2)*_a1/_b1;
    vh2tmp = 2*(_b1*_del2-_b2*_del2+2*_d2*_mu2)/(_del2*_del2-4*_d2*_lam);
    vs2tmp = -4*(2*_b1*_lam-2*_b2*_lam+_del2*_mu2)/(_del2*_del2-4*_d2*_lam);
    vsi2tmp = vs2tmp - vsrtmp*vsrtmp;
    if (vh2tmp >= 0 && vsi2tmp >= 0)
    {
        _localExtreme.push_back({sqrt(vh2tmp),vsrtmp,sqrt(vsi2tmp)});
        AppendLocalExtreme();

        _localExtreme.push_back({sqrt(vh2tmp),vsrtmp,-sqrt(vsi2tmp)});
        AppendLocalExtreme();
    }
    
    _Solved = true;
    
}

double cxSM_CP_reduced_a1::Vtotal(VD field_values, double scale)
{
    double vsr = field_values[1];
    return cxSM_CP::Vtotal(field_values,scale) + sqrt(2)*_a1*vsr/pow(scale,4);
}
VD cxSM_CP_reduced_a1::dVtotal(VD field_values, double scale)
{
    double vsr = field_values[1];
    VD res = cxSM_CP::dVtotal(field_values,scale);
    res[1] += sqrt(2)*_a1/pow(scale,3);
    return res;
}
double cxSM_CP_reduced_a1::V0_global(double scale)
{
    return Vtotal({_vev,_vsr,_vsi},scale);
}
void cxSM_CP_reduced_a1::PrintParameters()
{
    cout<<"Potential Parameter: "<<endl;
    cout<<"mu2:\t"<<_mu2<<endl;
    cout<<"b1:\t"<<_b1<<endl;
    cout<<"b2:\t"<<_b2<<endl;
    cout<<"lam:\t"<<_lam<<endl;
    cout<<"del2:\t"<<_del2<<endl;
    cout<<"del3:\t"<<_del3<<endl;
    cout<<"d1:\t"<<_d1<<endl;
    cout<<"d2:\t"<<_d2<<endl;
    cout<<"d3:\t"<<_d3<<endl;
    cout<<"a1:\t"<<_a1<<endl;

    cout<<"Physical Parameter: "<<endl;
    cout<<"vev:\t"<<_vev<<endl;
    cout<<"vsr:\t"<<_vsr<<endl;
    cout<<"vsi:\t"<<_vsi<<endl;
    cout<<"vs:\t"<<_vs<<endl;
    cout<<"alpha:\t"<<_alpha<<endl;
    cout<<"MHL:\t"<<_MHL<<endl;
    cout<<"MHH:\t"<<_MHH<<endl;
    cout<<"MHA:\t"<<_MHA<<endl;
    cout<<"theta1:\t"<<_theta1<<endl;
    cout<<"theta2:\t"<<_theta2<<endl;
    cout<<"theta3:\t"<<_theta3<<endl;
}