#ifndef _TOY_Z4_H_
#define _TOY_Z4_H_

#include "Basic_Model.h"

class ToyZ4 : public Basic_Model {
protected:
    double r;

    double _v1global;
    double _v2global;
    double _vr;

public:
    ToyZ4();
    ~ToyZ4(){};

    void Set_Potential_Parameters(double r);

    virtual double Vtotal(VD field_values, double scale = 1);
    virtual VD dVtotal(VD field_values, double scale = 1);
    virtual VVD d2Vtotal(VD field_values, double scale = 1);
    virtual double V0_global(double scale = 1);
    virtual VD QuarticCoupling(VD field_values);

    virtual void FindLocalMinima();
    virtual void PrintParameters();

    double Get_Vr() { return _vr; }
};

#endif  //_TOY_Z3_H_
