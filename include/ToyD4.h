#ifndef TOY_D4_H
#define TOY_D4_H

#include "Basic_Model.h"

class ToyD4: public Basic_Model
{
protected:
    double _g;

    double _v1global;
    double _v2global;
    double _v1;
    double _v2;
    // double _v3;

public:
    ToyD4();
    ~ToyD4(){};

    void Set_Potential_Parameters(double g);

    virtual double Vtotal(VD field_values, double scale = 1);
    virtual VD dVtotal(VD field_values, double scale = 1);
    virtual VVD d2Vtotal(VD field_values, double scale = 1);
    virtual double V0_global(double scale = 1);
    virtual VD QuarticCoupling(VD field_values);

    virtual void FindLocalMinima(); 
    virtual void PrintParameters();

    double Getf1(){return _v1;}
    double Getf2(){return _v2;}
};




#endif // TOY_A4_H