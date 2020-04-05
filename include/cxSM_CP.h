#ifndef CXSM_CP
#define CXSM_CP

#include "Basic_Model.h"
#include <Eigen/Dense>

class cxSM_CP: public Basic_Model
{
// ! Note that, in this version, the DW/soliton ODEs are written in terms of vev, Re(vs), Im(vs), not vev, |vs|, alpha !
protected:
// * V = mu2 |Phi|^2 + lam |Phi|^4 + del2/2 |Phi|^2|S|^2 + b2/2 |S|^2 + d2/4 |S|^4 
// * + ( del3/4 |Phi|^2 S^2 + b1/4 S^2 + d1/8 S^4 + d3/8 S^2 |S|^2 +c.c)
// Potential Parameter
    double _mu2; //- Z2
    double _lam; //- Z2
    double _del2; //- Z2
    double _b2; //- Z2
    double _d2; //- Z2
    double _del3;
    double _b1;
    double _d1;
    double _d3;
    
// Physical Parameter
    double _MHL;
    double _MHH;
    double _MHA;
    double _vev;
    double _vsr;
    double _vsi;
    double _theta1;
    double _theta2;
    double _theta3;

    double _MHL2;
    double _MHH2;
    double _MHA2;

// Another way as vs and alpha:
    double _vs;
    double _alpha;

// The mixing matrix elements
    Eigen::Matrix3d _R;


// Local parameter store local extreme points
    // From cxSM_Z2, no need to redefine;

private:
    double _mag(double vsr, double vsi);
    double _arg(double vsr, double vsi);
    void _GetThetas();
    void _GetR();

    // Used for stability 
    double _Vstab(double sth2, double c2phi);

public:
    cxSM_CP();
    ~cxSM_CP(){};

    void Set_Potential_Parameters(double mu2, double lam, double del2, double b2, double d2, double del3, double b1, double d1, double d3);
    void Set_Physical_Parameters(double vsr, double vsi, double MHH, double MHA, double theta1, double theta2, double theta3);

    virtual double Vtotal(VD field_values, double scale = 1);
    virtual VD dVtotal(VD field_values, double scale = 1);
    virtual VVD d2Vtotal(VD field_values, double scale = 1);
    virtual double V0_global(double scale = 1);

    virtual void FindLocalMinima(); 
    virtual bool CheckStability();
    virtual bool CheckUnitarity(double MAX=0.5);

    virtual void PrintParameters();

};

#endif //CXSM_CP