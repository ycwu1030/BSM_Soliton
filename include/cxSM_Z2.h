#ifndef CXSM_WITH_Z2
#define CXSM_WITH_Z2

#include "Basic_Model.h"
// #include "Potential.h"

class cxSM_Z2: public Basic_Model
{
protected:
// Potential Parameter
    double _mu2;
    double _lam;
    double _del2;
    double _b2;
    double _d2;
// Physical Parameter
    double _vs;
    double _theta;
    double _MHH;
    double _MHH2;

    double _vev; // This two will usually be the value from SM
    double _MHL; 
    double _MHL2;

public:
    cxSM_Z2();
    ~cxSM_Z2(){};

    void Set_Potential_Parameters(double mu2, double lam, double del2, double b2, double d2);
    void Set_Physical_Parameters(double vs, double theta, double MHH);

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

// class Toy: public Potential, public SM
// {
// // Toy model following Brawn's thesis Sec. E.2
// private:
//     double _lambda;
//     double _eta;
// public:
//     Toy(){};
//     Toy(double lambda, double eta);
//     ~Toy(){};

//     virtual double Vtotal(VD field_values, double scale = 1);
//     virtual VD dVtotal(VD field_values, double scale = 1);
//     virtual VVD d2Vtotal(VD field_values, double scale = 1);
//     virtual double V0_global(double scale = 1);

// };

#endif //CXSM_WITH_Z2