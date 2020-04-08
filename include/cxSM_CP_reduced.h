#ifndef CXSM_CP_REDUCED
#define CXSM_CP_REDUCED

#include "cxSM_CP.h"
#include <iostream>

class cxSM_CP_reduced: public cxSM_CP
{
// ! Note that, in this version, the DW/soliton ODEs are solved in terms of vev, Re(vs), Im(vs), not vev, |vs|, alpha !
protected:
// * V = mu2 |Phi|^2 + lam |Phi|^4 + del2/2 |Phi|^2|S|^2 + b2/2 |S|^2 + d2/4 |S|^4 
// * + ( b1/4 S^2 + d1/8 S^4 + c.c)
// Potential Parameter
    // * Inherit from cxSM_CP
    // double _mu2; //- Z2
    // double _lam; //- Z2
    // double _del2; //- Z2
    // double _b2; //- Z2
    // double _d2; //- Z2
    // double _del3; // Should be zero
    // double _b1;
    // double _d1;
    // double _d3; // Should be zero
    
// Physical Parameter
    // * Inherit from cxSM_CP
    // double _MHL;
    // double _MHH;
    // double _MHA;
    // double _vev;
    // double _vsr;
    // double _vsi;
    // double _theta1;
    // double _theta2;
    // double _theta3;

    // double _MHL2;
    // double _MHH2;
    // double _MHA2;

// Another way as vs and alpha:
    // * Inherit from cxSM_CP
    // double _vs;
    // double _alpha;

// The mixing matrix elements
    // * Inherit from cxSM_CP;
    // Eigen::Matrix3d _R;


// Local parameter store local extreme points
    // From cxSM_Z2, no need to redefine;

private:
    bool _GetTheta2Theta3(const double vsr, const double vsi, const double MHH, const double MHA, const double theta1, double &theta2, double &theta3);
    bool _GetAlphaTheta2(const double vs, const double MHH, const double MHA, const double theta1, double theta3, double &alpha, double &theta2);

    // * Inherit from cxSM_CP
    // double _mag(double vsr, double vsi);
    // double _arg(double vsr, double vsi);
    // void _GetThetas();
    // void _GetR();

    // // Used for stability 
    // double _Vstab(double sth2, double c2phi);

public:
    cxSM_CP_reduced();
    ~cxSM_CP_reduced(){};

    // void Set_Potential_Parameters(double mu2, double lam, double del2, double b2, double d2, double del3, double b1, double d1, double d3);
    // void Set_Physical_Parameters(double vsr, double vsi, double MHH, double MHA, double theta1, double theta2, double theta3);
    void Set_Potential_Parameters(double mu2, double lam, double del2, double b2, double d2, double b1, double d1);
    bool Set_Physical_Parameters_vsr_vsi_theta(double vsr, double vsi, double MHH, double MHA, double theta1);
    bool Set_Physical_Parameters_vs_theta(double vs, double MHH, double MHA, double theta1, double theta3);

    void GetVS(double &vs_real, double &vs_img){vs_real = _vsr; vs_img = _vsi;};
    void GetTheta(int id, double &theta)
    {
        switch (id)
        {
        case 1:
            theta = _theta1;
            return;
        case 2:
            theta = _theta2;
            return;
        case 3:
            theta = _theta3;
            return;
        default:
            theta = _theta2;
            return;
        }
    }
    friend std::ostream &operator<<(std::ostream &os, const cxSM_CP_reduced &mod){
        os << mod._MHL << "\t" << mod._MHH << "\t" << mod._MHA << "\t" << mod._vs;
        os << "\t" << mod._theta1 << "\t" << mod._theta2 << "\t" << mod._theta3;
        os << "\t" << mod._mu2 << "\t" << mod._b1 << "\t" << mod._b2;
        os << "\t" << mod._lam << "\t" << mod._del2 << "\t" << mod._d1;
        os << "\t" << mod._d2 << "\t" << mod._alpha;
        return os;
    };
    // * Inherit from cxSM_CP
    // virtual double Vtotal(VD field_values, double scale = 1);
    // virtual VD dVtotal(VD field_values, double scale = 1);
    // virtual VVD d2Vtotal(VD field_values, double scale = 1);
    // virtual double V0_global(double scale = 1);

    // virtual void FindLocalMinima(); 
    // virtual bool CheckStability();
    // virtual bool CheckUnitarity(double MAX=0.5);

    // virtual void PrintParameters();

};

#endif //CXSM_CP_REDUCED