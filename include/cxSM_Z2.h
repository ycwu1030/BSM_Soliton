#ifndef CXSM_WITH_Z2
#define CXSM_WITH_Z2

#include "SM.h"
#include "Potential.h"

class cxSM_Z2: public Potential, public SM
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

// Local parameter store local extreme points
    static const int _NExtremeMax = 6;
    int _NLocalExtreme;
    int _NLocalMinima;
    int _IndexInput;
    int _MinimaIndex[_NExtremeMax];
    double _localH[_NExtremeMax];
    double _localS[_NExtremeMax];
    bool _LocalMinimaQ[_NExtremeMax];
    double _Vtotal[_NExtremeMax];
    bool _Solved;

private:
    void SetLocalExtreme();

public:
    cxSM_Z2();
    ~cxSM_Z2(){};

    void Set_Potential_Parameters(double mu2, double lam, double del2, double b2, double d2);
    void Set_Physical_Parameters(double vs, double theta, double MHH);

    virtual double Vtotal(VD field_values, double scale = 1);
    virtual VD dVtotal(VD field_values, double scale = 1);
    virtual VVD d2Vtotal(VD field_values, double scale = 1);
    virtual double V0_global(double scale = 1);

    virtual void FindLocalMinima(); 
    virtual bool CheckStability();
    virtual bool CheckUnitarity(double MAX=0.5);
    virtual bool CheckGlobalMinimum();

    virtual void PrintParameters();
    virtual void PrintLocalMinima();

};

#endif //CXSM_WITH_Z2