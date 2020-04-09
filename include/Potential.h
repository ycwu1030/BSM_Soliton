#ifndef POTENTIAL_H
#define POTENTIAL_H
#include "VTypes.h"
#include <string>

class Potential
{
// When using any function in this class, suppose the input value of x and fields are with good scale
// Whether the output value is scaled depends on the input scale. output = true_value/pow(scale,n), where n is the energy dimension of each function.
private:
    int _Field_Dim; // The dimension of the field space

public:
    Potential(){_Field_Dim=0;}
    Potential(int Field_Dim){_Field_Dim = Field_Dim;}
    ~Potential(){};

    virtual double Vtotal(VD field_values, double scale = 1) {return 0;}
    virtual double V0_global(double scale = 1) {return 0;}
    virtual VD dVtotal(VD field_values, double scale = 1) {return VD(1,0);}
    virtual VVD d2Vtotal(VD field_values, double scale = 1) {return VVD(1,VD(1,0));}
    virtual VD QuarticCoupling(VD field_values){return VD(1,0);};

    int GetFieldDimension(){return _Field_Dim;}

    bool CheckHessian(VD field_values);

    double GetTotalEnergy(VD x, VVD fields);
    double GetTension(VD x, VVD fields);
    void DumpEnergyDensity(VD x, VVD fields, std::string filename);

};



#endif //POTENTIAL_H