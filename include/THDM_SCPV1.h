#ifndef THDM_SCPV1_H
#define THDM_SCPV1_H

#include "Basic_Model.h"
#include <Eigen/Dense>

class THDM_SCPV1: public Basic_Model
{
protected:
// * In total the potential is controlled by 8 parameters;
// * In interaction basis, they are m112 m222 m122 lam1 lam2 lam3 lam4 lam5
// * In physical basis,  we choose vev beta mh mH mA mpm alpha alphc as free parameters, theta and alphab will be obtained from those free parameters
// Potential Parameter
    double _m112;
    double _m222;
    double _m122;
    double _lam1;
    double _lam2;
    double _lam3;
    double _lam4;
    double _lam5;
// Physical Parameter
    double _vev;
    double _beta;
    double _theta;
    double _v1;
    double _v2r;
    double _v2i;
    double _MHL;
    double _MHL2;
    double _MHH;
    double _MHH2;
    double _MHA;
    double _MHA2;
    double _MHpm;
    double _MHpm2;
    double _alpha;
    double _alphab;
    double _alphac;

    // The mixing matrix for neutral scalars
    Eigen::Matrix3d _R;

public:
    THDM_SCPV1();
    ~THDM_SCPV1(){};

    void Set_Potential_Parameters(double m112, double m222, double m122, double lam1, double lam2, double lam3, double lam4, double lam5);
    void Set_Physical_Parameters(double beta, double m2, double m3, double mpm, double alpha, double alphac);

    virtual double Vtotal(VD field_values, double scale = 1);
    virtual VD dVtotal(VD field_values, double scale = 1);
    virtual VVD d2Vtotal(VD field_values, double scale = 1);
    virtual double V0_global(double scale = 1);
    virtual VD QuarticCoupling(VD field_values);

    virtual void FindLocalMinima(); 
    virtual bool CheckStability();
    virtual bool CheckUnitarity(double MAX=0.5);

    virtual void PrintParameters();

};

#endif //THDM_SCPV1_H