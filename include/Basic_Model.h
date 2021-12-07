#ifndef _BASIC_MODEL_H_
#define _BASIC_MODEL_H_

#include "Constants.h"
#include "Potential.h"
#include "VTypes.h"

namespace BSM_Soliton {
class SM {
public:
    SM();
    virtual ~SM(){};

    double Get_MZ() const { return MZ; }
    double Get_MZ2() const { return MZ2; }
    double Get_MW() const { return MW; }
    double Get_MW2() const { return MW2; }
    double Get_MT() const { return MT; }

    double Get_theta_w() const { return theta_w; }
    double Get_sin_theta_w() const { return sin(theta_w); }

    double Get_g_weak() const { return g_weak; }
    double Get_g_hyper() const { return g_hyper; }

    double Get_vev() const { return vev; }

protected:
    double MW;
    double MW2;
    double theta_w;
    double alpha_EW;
    double vev;
    double ee;
    double g_weak;
    double g_hyper;
    double yt;
};

class BaseModel : public SM {
public:
    BaseModel(int FieldSpaceDimension = 0);
    virtual ~BaseModel(){};

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

class SM {
public:
    SM();
    ~SM(){};

    double GetMZ() { return MZ; }
    double GetMZ2() { return MZ2; }
    double GetMW() { return MW; }
    double GetMW2() { return MW2; }
    double GetMT() { return MT; }
    double Getgweak() { return g_weak; }
    double Getghyper() { return gp_hyper; }
    double GetVEV() { return vev; }
    double GetSW() { return sin(thetaW); }

protected:
    double MW;
    double MW2;
    double thetaW;
    double alpha;
    double vev;       // (sqrt(2)GF)^(-0.5)
    double ee;        // ee = sqrt(4*Pi*alpha)
    double g_weak;    // g = ee/sw
    double gp_hyper;  // g' = ee/cw
    double yt;        // sqrt(2)mt/vev;
};

class Basic_Model : public SM, public Potential {
protected:
    int _N_VEVs;

    int _NLocalExtreme;
    int _NLocalMinima;
    int _IndexInput;  // This is the index for the local extreme that is supposed to be the global minimum.
    VI _MinimaIndex;
    VVD _localExtreme;
    VB _LocalMinimaQ;
    VD _Vtotal;
    bool _Solved;

    // The functions to get the local extreme, and setting the local extreme
    // Both should be re-implemented accordingly in each model
    void Clear_Local_Cache();
    virtual void FindLocalMinima() {
        if (_Solved) return;
        _NLocalExtreme = 0;
        _NLocalMinima = 0;
        _IndexInput = -1;
        _Solved = true;
    };
    virtual void AppendLocalExtreme();

public:
    Basic_Model();
    Basic_Model(int Field_Dim);
    ~Basic_Model(){};

    // Theoretical checks
    // The Unitarity and Stability checks should be re-implemented in each derived class
    virtual bool CheckStability() { return true; }
    virtual bool CheckUnitarity(double MAX = 0.5) { return true; }
    // The global minimum check can be re-used, if possible.
    virtual bool CheckGlobalMinimum();

    // Print parameters
    //  PrintParameters should be re-implemented in each derived class
    virtual void PrintParameters();
    // PrintLocalMinima can be re-used, no change is needed, but users can change it accordingly.
    virtual void PrintLocalMinima();

    int GetNLocalMinima() { return _NLocalMinima; }
    VD GetLocalMinima(int id) {
        if (id >= _NLocalMinima) return _localExtreme[_MinimaIndex[0]];
        return _localExtreme[_MinimaIndex[id]];
    }
};

#endif  // BASIC_MODEL_H
