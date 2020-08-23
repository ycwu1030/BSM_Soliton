#ifndef THDM_SCPV1_H
#define THDM_SCPV1_H

#include "Basic_Model.h"
#include <Eigen/Dense>

class THDM_CPC: public Basic_Model
{
protected:
// * In total the potential is controlled by 8 parameters;
// * In interaction basis, they are m112 m222 m122 lam1 lam2 lam3 lam4 lam5
// * In physical basis,  we choose vev beta alpha mh mH mA mpm m122 as free parameters,
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
    double _vev2;
    double _beta;
    double _v1;
    double _v2;
    double _MHL;
    double _MHL2;
    double _MHH;
    double _MHH2;
    double _MHA;
    double _MHA2;
    double _MHpm;
    double _MHpm2;
    double _alpha;

    // Auxillary
    double _tb;
    double _sb;
    double _cb;
    double _ca;
    double _sa;
    double _cba;
    double _sba;


private:


public:
    THDM_CPC();
    ~THDM_CPC(){};

    // bool Set_Potential_Parameters(double m112, double m222, double m122, double lam1, double lam2, double lam3, double lam4, double lam5);
    bool Set_Physical_Parameters(double beta=atan(2.0), double alpha=(atan(2.0)-Pi/2.0), double mhh=400.0, double mha=450.0, double mpm=500.0, double m122=80000.0);

    virtual double Vtotal(VD field_values, double scale = 1);
    virtual VD dVtotal(VD field_values, double scale = 1);
    virtual VVD d2Vtotal(VD field_values, double scale = 1);
    virtual double V0_global(double scale = 1);
    virtual VD QuarticCoupling(VD field_values);

    virtual void FindLocalMinima(); 
    virtual bool CheckStability();
    virtual bool CheckUnitarity(double MAX=0.5);

    virtual void PrintParameters();

    std::string repr();
    friend std::ostream &operator<<(std::ostream &os, const THDM_CPC &mod);
};

#endif //THDM_SCPV1_H