#ifndef _TOY_Z4_H_
#define _TOY_Z4_H_

#include "Basic_Model.h"

namespace ToyTest {
class ToyZ4 : public BSM_Soliton::BaseModel {
public:
    ToyZ4();
    ~ToyZ4(){};

    void Set_Potential_Parameters(double beta, double delta_betas);

    virtual double V(VD field_values, double scale = 1) override;
    virtual VD dV(VD field_values, double scale = 1) override;
    virtual VVD d2V(VD field_values, double scale = 1) override;
    virtual double V_min(double scale = 1) override;
    virtual VD Quartic_Couplings(VD field_values) override;

    double Get_Vr() { return _vr; }
    void Print_Parameters();

protected:
    double r;
    double beta;
    double delta_beta;  // 1-beta

    double _v1global;
    double _v2global;
    double _vr;

    virtual void Calculate_Local_Extrema() override;
};
}  // namespace ToyTest

#endif  //_TOY_Z3_H_
