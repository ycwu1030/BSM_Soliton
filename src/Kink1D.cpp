#include "Kink1D.h"
#include "GSL_Wraper.h"
#include <cmath>
#include <iostream>
#include <iomanip>

using namespace std;

VD func_for_rkqc(double x, VD y, void *param)
{
    Kink1D *mod = (Kink1D*)param;
    return mod->equationOfMotion(x,y);
}
Kink1D::Kink1D(double phi_left_, double phi_right_, V1D V_, V1D dV_, V1D d2V_, double phi_eps_rel_, bool scaled_)
{
    V = V_;
    dV = dV_;
    d2V = d2V_;

    _rk_calculator.SetDOF(2); // * y0 = phi, y1 = dphi;
    _rk_calculator.SetODE(func_for_rkqc);
    _rk_calculator.SetParams(this);

    phi_left = phi_left_;
    phi_right = phi_right_;
    scaled = scaled_;
    findBarrierLocation();

    findScales();

    phi_eps_rel = phi_eps_rel_;
    phi_eps_abs = phi_eps_rel*abs(phi_left-phi_right);
}
VD Kink1D::equationOfMotion(double z, VD Y)
{
    VD res(2);
    res[0] = Y[1];
    res[1] = scaled?dV_scaled(Y[0])/phiscale*zscale*zscale:dV(Y[0]);
    return res;
}
double func_for_findBarrierTop(double x, void *params)
{
    Kink1D *mod = (Kink1D*) params;
    return -mod->VvalatX(x);
}
void Kink1D::findBarrierLocation()
{
    // * First find the meta and absolute minima location
    if (V(phi_right) <= V(phi_left))
    {
        phi_abs = phi_right;
        phi_meta = phi_left;
        match_left_meta = true;
    }
    else
    {
        phi_abs = phi_left;
        phi_meta = phi_right;
        match_left_meta = false;
    }

    // * Find the position where V(phi_bar) === V(phi_meta);
    // cout<<"abs(V(phi_abs)-V(phi_meta)): "<<abs(V(phi_abs)-V(phi_meta))<<endl;
    if (abs(V(phi_abs)-V(phi_meta))<1e-3)
    {
        // * If the V at abs/meta point are quite close, just set the barrier at phi_abs;
        phi_bar = phi_abs;
    }
    else
    {
        // * Otherwise, use binary search to find the barrier position
        double phi_tol = abs(phi_abs-phi_meta)*1e-12;
        double V_meta = V(phi_meta);
        double phiH = phi_meta;
        double phiL = phi_abs;
        double phiM = (phiH+phiL)/2;
        double V0;
        while (abs(phiH-phiL) > phi_tol)
        {
            V0 = V(phiM);
            if (V0 > V_meta)
            {
                phiH = phiM;
            }
            else
            {
                phiL = phiM;
            }
            phiM = (phiH + phiL)/2;
        }
        phi_bar = phiM;
    }
    
    // * Then between phi_meta and phi_bar, we try to find the top of the barrier
    double phi_tol = abs(phi_bar - phi_meta)*1e-8;
    double x1 = min(phi_bar,phi_meta);
    double x2 = max(phi_bar,phi_meta);
    // cout<<"x1: "<<x1<<" x2: "<<x2<<endl;
    phi_bar_top = find_min_arg_gsl_wraper(func_for_findBarrierTop,this,x2,x1,phi_tol);
    // cout<<"phi_bar_top: "<<phi_bar_top<<endl;
    if (phi_bar_top + phi_tol > x2 || phi_bar_top - phi_tol < x1)
    {
        cout<<"In findBarrierLocation: No Barrier for the potential, can't find the top position"<<endl;
    }
}
 
void Kink1D::findScales()
{
    double Vtop = V(phi_bar_top) - V(phi_abs);
    double xtop = phi_bar_top - phi_abs;
    if (Vtop <= 0)
    {
        cout<<"In findZScale: No barrier for the potential: non-positive barrier height."<<endl;
    }

    // * From phi_bar_top to phi_abs, assume that it is a quartic potential
    // * V(x) = Vtop/xtop^4 (x^2-xtop^2)^2; 
    // * The solution to d^2x/dt^2 = dV(x)/dx is 
    // * x(t) = xtop*tanh(sqrt(2*Vtop)/xtop*t)
    zscale = abs(xtop)/sqrt(abs(2*Vtop));
    phiscale = abs(phi_abs);

    V_scaled = [&](double field_hat){
        return V(field_hat*phiscale);
    };
    dV_scaled = [&](double field_hat){
        return dV(field_hat*phiscale);
    };
    d2V_scaled = [&](double field_hat){
        return d2V(field_hat*phiscale);
    };
    phi_abs_scaled = phi_abs/phiscale;
    phi_meta_scaled = phi_meta/phiscale;
    phi_bar_scaled = phi_bar/phiscale;
    phi_bar_top_scaled = phi_bar_top/phiscale;
}
struct cubic_param
{
    double y0;
    double dy0;
    double y1;
    double dy1;
    double diff;
};
double cubicInterpolation(double x, void *param)
{
    cubic_param* mod = (cubic_param*)param;
    double mt = 1-x;
    double c3 = mod->y1;
    double c2 = mod->y1 - mod->dy1/3.0;
    double c1 = mod->y0 + mod->dy0/3.0;
    double c0 = mod->y0;
    return c0*pow(mt,3) + 3*c1*mt*mt*x + 3*c2*mt*x*x + c3*pow(x,3) - mod->diff;
}
tuple<double, VD, CONVERGENCETYPE> Kink1D::integrateProfile(VD y0, VD y_desired, double dr0, double phi_eps_rel_, double drmin, double rmax)
{
    VD y_final_desired_value = y_desired;
    double phi_final = y_desired[0];
    VD y_tol_scale(2);// = {abs(phi_final-phi_bar_top), abs(phi_final-phi_bar_top)/zscale};
    if (scaled)
    {
        y_tol_scale[0] = abs(phi_final-phi_bar_top_scaled);
        y_tol_scale[1] = y_tol_scale[0];
    }
    else
    {
        y_tol_scale[0] = abs(phi_final-phi_bar_top);
        y_tol_scale[1] = abs(phi_final-phi_bar_top)/zscale;
    }
    
    
    // cout<<"desired: "<<y_desired<<endl;
    // cout<<"y-tol-scale: "<<y_tol_scale<<endl;
    VD y_diff;
    double dr_guess = dr0;
    double dr_did,dr_next;
    double r = 0;
    VD y = y0;
    VD dydr = equationOfMotion(r, y);
    double r_cache;
    VD y_cache;
    VD dydr_cache;
    VD y_scale;
    VD y_inter(2);
    int ysign = signbit(y[0] - phi_final)?-1:1;

    CONVERGENCETYPE convergQ = NONE;

    cubic_param inter_param;
    double x;
    // int steps = 0;
    // cout<<"Desired: "<<y_final_desired_value<<endl;
    // cout<<steps<<"\t"<<r<<"\t"<<y<<endl;
    while (true)
    {
        // y_scale = {abs(phi_final-phi_bar_top),abs(phi_final-phi_bar_top)/zscale};
        // steps++;
        y_scale = abs(y) + abs(dydr*dr_guess);
        r_cache = r;
        y_cache = y;
        dydr_cache = dydr;

        _rk_calculator._RKQC_SingleStep(r_cache,y_cache,dydr_cache,dr_guess,phi_eps_rel_,y_scale,dr_did,dr_next);
        dydr_cache = equationOfMotion(r_cache,y_cache);

        // cout<<"Step-"<<steps<<endl;
        // cout<<"\t"<<r_cache<<"\t"<<y_cache<<endl;
        y_diff = abs(y_cache - y_final_desired_value);

        if (r_cache > rmax)
        {
            convergQ = TOOLARGEZ;
            break;
        }
        
        if (dr_did < drmin)
        {
            convergQ = TOOSMALLSTEP;
            break;
        }
        
        // cout<<"\t"<<y_diff<<endl;
        // cout<<"\t"<<y_tol_scale<<"\t"<<phi_eps_rel_<<endl;
        if (y_diff[0] < y_tol_scale[0]*phi_eps_rel_ && y_diff[1] < y_tol_scale[1]*phi_eps_rel_)
        {
            // cout<<"desired: "<<y_desired<<endl;
            // cout<<y_diff<<endl;
            // cout<<y_tol_scale<<endl;
            convergQ = CONVERGED;
            r = r_cache;
            y = y_cache;
            break;
        }

        if (y_cache[1]*ysign > 0)
        {
            // * This means that the field is rolling back, will never reach phi_final;
            // * So we try to stop the integrate at y1=0 (dphi/dr = 0)
            convergQ = UNDERSHOOT;
            inter_param = {y[1], dydr[1]*dr_did, y_cache[1], dydr_cache[1]*dr_did,0};
            x = find_root_gsl_wraper(&cubicInterpolation,&inter_param,1,0);
            r += dr_did*x;
            y_inter[1] = cubicInterpolation(x,&inter_param);
            inter_param = {y[0], dydr[0]*dr_did, y_cache[0], dydr_cache[0]*dr_did,0};
            y_inter[0] = cubicInterpolation(x,&inter_param);
            y = y_inter;
            break;
        }

        if ((y_cache[0]-phi_final)*ysign<0)
        {
            // * This means that the field is already passing the desired ending point
            // * Then we try to stop the integrate at y0 = phi_final;
            convergQ = OVERSHOOT;
            inter_param = {y[0],dydr[0]*dr_did,y_cache[0],dydr_cache[0]*dr_did,phi_final};
            x = find_root_gsl_wraper(&cubicInterpolation,&inter_param,1,0);
            r += dr_did*x;
            inter_param = {y[1],dydr[1]*dr_did,y_cache[1],dydr_cache[1]*dr_did,0};
            y_inter[1] = cubicInterpolation(x,&inter_param);
            inter_param = {y[0],dydr[0]*dr_did,y_cache[0],dydr_cache[0]*dr_did,0};
            y_inter[0] = cubicInterpolation(x,&inter_param);
            y = y_inter;
            break;
        }
        r = r_cache;
        y = y_cache;
        dydr = dydr_cache;
        dr_guess = dr_next;
    }
    y_diff = abs(y-y_final_desired_value);
    if (y_diff[0] < y_tol_scale[0]*phi_eps_rel_ && y_diff[1] < y_tol_scale[1]*phi_eps_rel_)
    {
        convergQ = CONVERGED;
    }
    return make_tuple(r,y,convergQ);
}
tuple<VD, VD, VD, double> Kink1D::integrateAndSaveProfile(VD R, VD y0, double dr, double phi_eps_rel_, double drmin)
{
    int N = R.size();
    double r0 = R[0];
    VVD Yout(y0.size(),VD(N,0));
    Yout[0][0] = y0[0];
    Yout[1][0] = y0[1];
    VD dydr0 = equationOfMotion(r0,y0);
    double Rerr = NAN;

    int i = 1;
    double r = r0;
    VD y = y0;
    VD dydr = dydr0;
    double r_cache;
    VD y_cache;
    VD dydr_cache;
    double dr_guess = dr;
    double dr_did,dr_next;
    VD y_scale(2);
    cubic_param inter_param;
    while (i<N)
    {
        y_scale = abs(y)+abs(dydr*dr_guess);
        r_cache = r;
        y_cache = y;
        dydr_cache = dydr;
        _rk_calculator._RKQC_SingleStep(r_cache,y_cache,dydr_cache,dr_guess,phi_eps_rel_,y_scale,dr_did,dr_next);
        if (dr_did < drmin)
        {
            y_cache = y + (y_cache-y)*drmin/dr_did;
            dr_did = drmin;
            dr_next = drmin;
            r_cache = r + dr_did;
            if (!(std::isnan(Rerr)))
            {
                Rerr = r_cache;
            }
        }
        dydr_cache = equationOfMotion(r_cache,y_cache);
        if (r < R[i] && R[i] <= r_cache)
        {
            while (i < N && r < R[i] && R[i] <= r_cache)
            {
                double x = (R[i]-r)/dr_did;
                inter_param = {y[0], dr_did*dydr[0], y_cache[0], dr_did*dydr_cache[0], 0};
                Yout[0][i] = cubicInterpolation(x, &inter_param);
                inter_param = {y[1], dr_did*dydr[1], y_cache[1], dr_did*dydr_cache[1], 0};
                Yout[1][i] = cubicInterpolation(x, &inter_param);
                i += 1;
            }   
        }

        r = r_cache;
        y = y_cache;
        dydr = dydr_cache;
        dr_guess = dr_next;
    }
    
    return make_tuple(R,Yout[0],Yout[1],Rerr);
}
tuple<VD,VD,VD,VD> Kink1D::findProfile(double dphi0_tol_rel, double phi_tol_rel, int npoints, double rmax)
{
    double dphi0_initial_guess = sqrt(2*abs(V(phi_bar_top)-V(phi_abs)));
    if (scaled)
    {
        dphi0_initial_guess = dphi0_initial_guess*zscale/phiscale;
    }
    dphi0_initial_guess = phi_bar_top<phi_abs?dphi0_initial_guess:-dphi0_initial_guess;
    // double dphi0_tol_rel = dphi0_tol_rel;// * dphi0_initial_guess;
    double dphi0_min_rel = 0.9;//*dphi0_initial_guess;
    double dphi0_max_rel = 1.1;//*dphi0_initial_guess;

    double dr0 = scaled?5.0/npoints:5.0*zscale/npoints;
    double drmin = dr0*1e-2;
    if (!scaled)
    {
        rmax *= zscale;
    }

    double rf_abs = NAN;
    double phi0 = scaled?phi_bar_top_scaled:phi_bar_top;
    double dphi0_rel = (dphi0_min_rel+dphi0_max_rel)/2.0;
    double dphi0 = dphi0_rel*dphi0_initial_guess;
    VD y0_abs,y0_meta;
    VD y_desired;
    VD yf;
    CONVERGENCETYPE ctype;

    // * Shooting to reach absolute minimum
    y_desired = {phi_abs, 0};
    if (scaled)
    {
        y_desired = y_desired/phiscale;
    }
    while (true)
    {
        y0_abs = {phi0, dphi0};
        tie(rf_abs,yf,ctype) = integrateProfile(y0_abs,y_desired,dr0,phi_tol_rel,drmin,rmax);
        if (ctype == CONVERGED)
        {
            // cout<<"CONVERGED: "<<dphi0<<endl;
            break;
        }
        else if (ctype == UNDERSHOOT)
        {
            // cout<<setprecision(15)<<"UNDERSHOOT: ("<<dphi0_min_rel<<","<<dphi0_rel<<","<<dphi0_max_rel<<") -> (";
            dphi0_min_rel = dphi0_rel;
            dphi0_rel = (dphi0_min_rel+dphi0_max_rel)/2.0;
            // cout<<dphi0_min_rel<<","<<dphi0_rel<<","<<dphi0_max_rel<<")"<<endl;
            dphi0 = dphi0_rel * dphi0_initial_guess;
            
        }
        else if (ctype == OVERSHOOT)
        {
            // cout<<setprecision(15)<<"OVERSHOOT: ("<<dphi0_min_rel<<","<<dphi0_rel<<","<<dphi0_max_rel<<") -> (";
            dphi0_max_rel = dphi0_rel;
            dphi0_rel = (dphi0_min_rel+dphi0_max_rel)/2.0;
            // cout<<dphi0_min_rel<<","<<dphi0_rel<<","<<dphi0_max_rel<<")"<<endl;
            dphi0 = dphi0_rel * dphi0_initial_guess;
        }
        
        if (abs(dphi0_max_rel - dphi0_min_rel) < dphi0_tol_rel)
        {
            break;
        }
    }

    double dphi0_solution = dphi0;
    // * Find the r_max for the meta side
    double rf_meta = NAN;
    y_desired = {phi_meta, 0};
    y0_meta = {phi_bar_top, -dphi0_solution};
    if (scaled)
    {
        y_desired[0] = y_desired[0]/phiscale;
        y0_meta[0] = y0_meta[0]/phiscale;
    }
    
    tie(rf_meta,yf,ctype)=integrateProfile(y0_meta,y_desired,dr0,phi_tol_rel,drmin,rmax);

    if (ctype == UNDERSHOOT)
    {
        cout<<"In findProfile: converge in absolute minimum side but undershoot in meta-stable side. Something wrong in the potential!"<<endl;
    }
    
    // * Save the profile from phi_bar_top to phi_abs;
    VD R_abs(npoints);
    for (int i = 0; i < npoints; i++)
    {
        R_abs[i] = 0 + i*rf_abs/(npoints - 1.0);
    }
    VD phi_abs;
    VD dphi_abs;
    double Rerr_abs;
    y0_abs = {phi_bar_top,dphi0_solution};
    if (scaled)
    {
        y0_abs[0] = y0_abs[0]/phiscale;
    }
    tie(R_abs,phi_abs,dphi_abs,Rerr_abs) = integrateAndSaveProfile(R_abs,y0_abs,dr0,phi_tol_rel,drmin);

    VD R_meta(npoints);
    for (int i = 0; i < npoints; i++)
    {
        R_meta[i] = 0 + i*rf_meta/(npoints - 1.0);
    }
    VD phi_meta;
    VD dphi_meta;
    double Rerr_meta;
    y0_meta = {phi_bar_top,-dphi0_solution};
    if (scaled)
    {
        y0_meta[0] = y0_meta[0]/phiscale;
    }
    tie(R_meta,phi_meta,dphi_meta,Rerr_meta) = integrateAndSaveProfile(R_meta,y0_meta,dr0,phi_tol_rel,drmin);
    
    // * Combine the two side
    VD R_final;
    VD phi_final;
    VD dphi_final;
    VD Rerr_final;
    if (match_left_meta)
    {
        // * need to flip meta part;

        // * First insert the abs part as it is
        R_final = R_abs;
        phi_final = phi_abs;
        dphi_final = dphi_abs;
        Rerr_final.push_back(Rerr_abs);

        // * Then insert the reversed meta part
        R_meta = -R_meta;
        dphi_meta = -dphi_meta;
        R_final.insert(R_final.begin(),R_meta.rbegin(),R_meta.rend()-1);
        phi_final.insert(phi_final.begin(),phi_meta.rbegin(),phi_meta.rend()-1);
        dphi_final.insert(dphi_final.begin(),dphi_meta.rbegin(),dphi_meta.rend()-1);
        Rerr_final.insert(Rerr_final.begin(),-Rerr_meta);
    }
    else
    {
        // * need to flip abs part;

        // * First insert the meta part as it is
        R_final = R_meta;
        phi_final = phi_meta;
        dphi_final = dphi_meta;
        Rerr_final.push_back(Rerr_meta);

        // * Then insert the reversed abs part
        R_abs = -R_abs;
        dphi_abs = -dphi_abs;
        R_final.insert(R_final.begin(),R_abs.rbegin(),R_abs.rend()-1);
        phi_final.insert(phi_final.begin(),phi_abs.rbegin(),phi_abs.rend()-1);
        dphi_final.insert(dphi_final.begin(),dphi_abs.rbegin(),dphi_abs.rend()-1);
        Rerr_final.insert(Rerr_final.begin(),-Rerr_abs);
    }
    
    if (scaled)
    {
        R_final = R_final*zscale;
        phi_final = phi_final*phiscale;
        dphi_final = dphi_final*phiscale;
        Rerr_final = Rerr_final*zscale;
    }
    
    return make_tuple(R_final,phi_final,dphi_final,Rerr_final);
}
tuple<VD, VD> Kink1D::evenlySpacedPhi(VD phi, VD dphi, int npoint, int k)
{
    // phi.insert(phi.begin(),phi_left);
    // phi.insert(phi.end(),phi_right);
    // dphi.insert(dphi.begin(),0.0);
    // dphi.insert(dphi.end(),0.0);
    VVD fullPhi = transpose({phi,dphi});
    sort(fullPhi.begin(),fullPhi.end(),[](VD x1, VD x2){return x1[0]<x2[0];});
    VVD::iterator iter = unique(fullPhi.begin(),fullPhi.end(),[](VD x1, VD x2){return x1[0]==x2[0];});
    fullPhi.resize(distance(fullPhi.begin(),iter));
    fullPhi=transpose(fullPhi);

    GSL_Spline_Inter inter;
    inter.SetData(&fullPhi[1],&fullPhi[0]);
    VD p = linspace(phi.front(),phi.back(),npoint);

    return make_tuple(p,inter.valAt(p));
}