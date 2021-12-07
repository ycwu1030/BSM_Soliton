#ifndef _POTENTIAL_H_
#define _POTENTIAL_H_
#include <string>

#include "VTypes.h"

namespace BSM_Soliton {
class Potential {
public:
    Potential(int FieldSpaceDimension = 0);
    virtual ~Potential(){};

    // * Several functions used to calculation the potential value, the derivatives etc.
    // * The scale an energy scale, which will be used to normalize the potential as well as the field value;
    virtual double V(VD field_values, double scale = 1) { return 0; }
    virtual double V_min(double scale = 1) { return 0; }
    virtual VD dV(VD field_values, double scale = 1) { return VD(Field_Space_Dimension, 0); }
    virtual VVD d2V(VD field_values, double scale = 1) {
        return VVD(Field_Space_Dimension, VD(Field_Space_Dimension, 0));
    }

    // * Used in relaxation method for soliton to set scale for position space.
    virtual VD Quartic_Couplings(VD field_values) { return VD(Field_Space_Dimension, 0); }

    // * Theoretical test on the potential
    virtual bool Check_Stability() { return true; }
    virtual bool Check_Unitarity(double max_eigen_value = 0.5) { return true; }
    virtual bool Check_Global_Minimum();

    void Print_Local_Extrema() const;

    int Get_Local_Minima_Size() const { return N_Local_Minima; }
    VD Get_Local_Minimum(unsigned id = 0) const;
    int Get_Field_Space_Dimension() const { return Field_Space_Dimension; }

protected:
    const int Field_Space_Dimension;
    VD Input_Minimum;
    bool Solved;

    int Input_Minimum_ID;
    int N_Local_Extrema;
    int N_Local_Minima;
    VI Minima_Index;
    VVD Local_Extrema;
    VB is_Local_Minima;
    VD Potential_at_Extrema;

    void Find_Local_Extrema();
    virtual void Calculate_Local_Extrema();
    void Add_Local_Extremum(VD local_extremum);
    void Clear_Local_Cache();

    bool Check_Hessian_Matrix(VD field_values);
};
}  // namespace BSM_Soliton

class Potential {
    // When using any function in this class, suppose the input value of x and fields are with good scale
    // Whether the output value is scaled depends on the input scale. output = true_value/pow(scale,n), where n is the
    // energy dimension of each function.
private:
    int _Field_Dim;  // The dimension of the field space

public:
    Potential() { _Field_Dim = 0; }
    Potential(int Field_Dim) { _Field_Dim = Field_Dim; }
    ~Potential(){};

    virtual double Vtotal(VD field_values, double scale = 1) { return 0; }
    virtual double V0_global(double scale = 1) { return 0; }
    virtual VD dVtotal(VD field_values, double scale = 1) { return VD(1, 0); }
    virtual VVD d2Vtotal(VD field_values, double scale = 1) { return VVD(1, VD(1, 0)); }
    virtual VD QuarticCoupling(VD field_values) {
        return VD(1, 0);
    };  // Used in relaxation method for soliton, used to set the scale in z.

    int GetFieldDimension() { return _Field_Dim; }

    bool CheckHessian(VD field_values);

    double GetTotalEnergy(VD x, VVD fields);
    double GetTension(VD x, VVD fields);
    VD GetEnergyDensity(VD x, VVD fields);
    double GetWallWidth(VD x, VVD fields, double criteria = 0.64);
    void DumpFullSolution(VD x, VVD fields, std::string filename);
    void DumpEnergyDensity(VD x, VVD fields, std::string filename);
};

#endif  // _POTENTIAL_H_
