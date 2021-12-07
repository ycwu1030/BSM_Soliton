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
