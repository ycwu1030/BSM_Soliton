#ifndef CXSM_Z2_BIASED
#define CXSM_Z2_BIASED

#include "cxSM_Z2.h"
#include "CubicSolver.h"

class cxSM_Z2_biased: public cxSM_Z2
{
protected:
// ! Just define extra parameters based on cxSM_Z2
// Potential Parameter
    double _del1;
    double _a1;
    double _c1;
    double _c2;
    
// Physical Parameter
    double _MHA;
    double _MHA2;

// Local parameter store local extreme points
    // From cxSM_Z2, no need to redefine;

private:
    CubicSolver _Solver;
    void SolveCubicEquation(double A[4], double *results, int &NSolution);

public:
    cxSM_Z2_biased();
    ~cxSM_Z2_biased(){};

    void Set_Potential_Parameters(double mu2, double lam, double del2, double b2, double d2, double del1, double a1, double c1, double c2);
    void Set_Physical_Parameters(double vs, double theta, double MHH, double MHA, double del1, double c1, double c2);

    virtual double Vtotal(VD field_values, double scale = 1);
    virtual VD dVtotal(VD field_values, double scale = 1);
    virtual VVD d2Vtotal(VD field_values, double scale = 1);
    virtual double V0_global(double scale = 1);

    virtual void FindLocalMinima(); 
    // virtual bool CheckStability();
    // virtual bool CheckUnitarity(double MAX=0.5);
    // virtual bool CheckGlobalMinimum();

    virtual void PrintParameters();

    bool GetBiasedMirrorMinimum(VD &mirror);

};

#endif //CXSM_Z2_BIASED