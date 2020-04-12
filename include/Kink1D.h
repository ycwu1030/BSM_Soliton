#ifndef KINK_1D_H
#define KINK_1D_H

#include "VTypes.h"
#include "RungeKutta.h"
#include <tuple>

enum CONVERGENCETYPE
{
    UNDERSHOOT = -1,
    CONVERGED = 0,
    OVERSHOOT = 1,
    TOSMALLSTEP = 2,
    TOLARGEZ = 3,
    NONE      = -9
};

class Kink1D
{
private:
    double phi_left;
    double phi_right;

    // * Actually, for a stable DW, V(phi_left) = V(phi_right)
    // * But if we deal with some small biased term, we still need to 
    // * treat meta minimum and absolute minimum differently.
    bool match_left_meta; // True for phi_left = phi_meta, False otherwise
    double phi_meta; 
    double phi_abs;
    double phi_bar;
    double phi_bar_top;

    V1D V;
    V1D dV;
    V1D d2V;

    double phi_eps_rel;
    double phi_eps_abs;

    double zscale;
    double phiscale;

    RungeKutta _rk_calculator;

    void findBarrierLocation();
    void findScales();

    // * Start at phi_bar_top, where we set z = 0, integrate towards phi_abs/phi_meta to get the profile.
    // * Note, no matter for phi_abs or phi_meta, z is positive.
    // * So for one of them, we need to flip later.
    // * Only when from phi_bar_top to phi_abs integral converged, we can continue to consider the phi_meta part.
    std::tuple<double, VD, CONVERGENCETYPE> integrateProfile(VD y0, VD y_desire, double dr0, double phi_eps_rel_, double drmin, double rmax);

    // * integrate and also save the field profile
    std::tuple<VD, VD, VD, double> integrateAndSaveProfile(VD R, VD y0, double dr, double phi_eps_rel_, double drmin);

public:
    Kink1D(double phi_left_, double phi_right_, V1D V_, V1D dV_, V1D d2V_, double phi_eps_rel_ = 1e-3);
    ~Kink1D(){};

    VD equationOfMotion(double z, VD Y);
    std::tuple<VD,VD,VD,VD> findProfile(double dphi0_tol_rel=1e-12, double phi_tol_rel = 1e-8, int npoints = 200, double rmax=1e4);

    double VvalatX(double x){return V(x);};
};



#endif //KINK_1D_H