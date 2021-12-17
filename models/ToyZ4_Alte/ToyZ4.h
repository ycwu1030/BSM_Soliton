#ifndef _TOY_Z4_H_
#define _TOY_Z4_H_

#include "Basic_Model.h"

namespace ToyTest {
class ToyZ4 : public Basic_Model {
protected:
    double beta;
    double delta_beta;  // 1-beta

    double _v1global;
    double _v2global;
    double _vr;

public:
    ToyZ4();
    ~ToyZ4(){};

    void Set_Potential_Parameters(double beta, double delta_betas);

    virtual double Vtotal(VD field_values, double scale = 1);
    virtual VD dVtotal(VD field_values, double scale = 1);
    virtual VVD d2Vtotal(VD field_values, double scale = 1);
    virtual double V0_global(double scale = 1);
    virtual VD QuarticCoupling(VD field_values);

    virtual void FindLocalMinima();
    virtual void PrintParameters();

    double Get_Vr() { return _vr; }
};
}  // namespace ToyTest

#endif  //_TOY_Z3_H_
