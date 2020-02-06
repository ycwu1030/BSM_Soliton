/*
 * @Description  : Following Numerical Recipes in C (1988 version), Section 15.1 and 15.2
 * @Author       : Yongcheng Wu
 * @Date         : 2020-02-04 14:23:09
 * @LastEditors  : Yongcheng Wu
 * @LastEditTime : 2020-02-06 14:37:03
 */
#ifndef RungeKutta_H
#define RungeKutta_H

#include "VTypes.h"
#include <string>

typedef VD (*F_ODE)(double x, VD y);

class RungeKutta
{
private:
    int _DOF;
    double _x_begin, _x_end;
    VD _X; // Position, dimension depends on whether we use the adaptive stepsize control.
    VVD _Y; // The vector (as the same dimension as _X) of vector (dimension _DOF) of Y
    VVD _dYdX; // The vector (as the same dimension as _X) of vector (dimension _DOF) of dY/dX
    VD _BOUND_CONDITION;  // The boundary condition at the starting point
    VD _Y_SCALE; // The scale 
    F_ODE _derivs; // The ODE functions, arguments are the position and the Y values 

    void _RESET(); // Whenever we change the DOF, boundary condition or the ODE function, we need to reset X,Y,dYdX and be prepared for re-do the RK iteration;
    void _RK4_SingleStep(double X_CUR, VD Y_CUR, VD dY_CUR, double step_size, VD &dY_NEXT); // This is the usual 4th-order Runge-Kutta method, take one step forward.
    void _RKQC_SingleStep(double &X, VD &Y, VD dY, double step_size_guess, double eps, VD Y_Scale, double &step_size_did, double &step_size_further); // This is the Runge-Kutta method with quality controlled, which will achieve 5-th order accuracy. (adaptive stepsize)

// Following are used for two-point boundary problem
    std::vector<bool> _BOUND_begin_Q;
    std::vector<bool> _BOUND_end_Q;
    VD _BOUND_begin;
    VD _BOUND_end;
    int _N_BOUND_begin;
    int _N_BOUND_end; // _N_BOUND_end = _DOF - _N_BOUND_begin;
    VD _Score(VD _Y_END); // _Y_END in dimension _DOF
    void _Load(VD BOUND_GUESS); // Loading the initial condition at one boundary; BOUND_GUESS in dimension _N_BOUND_end;
    void _SHOOTING(VD BOUND_GUESS, VD delta_Bound, VD &BOUND_FURTHER, VD &SCORES); // Shooting method for two point boundary problem; BOUND_GUESS/delta_Bound in dimension _N_BOUND_end;


public:
    RungeKutta();
    RungeKutta(int DOF);
    ~RungeKutta(){};

    void SetDOF(int DOF=1);
    void SetBound(double x_begin, double x_end, VD BOUND); // This is used for usual ODE problem where the boundary conditions are given only at x_begin;
    void SetBound(double x_begin, double x_end, VD BOUND_begin, VD BOUND_end, std::vector<bool> At_Begin_Q, std::vector<bool> At_End_Q);
    void SetODE(F_ODE derivs);

    void ODEINTEGRAL(double step_start, double eps=1e-6);

    void SHOOTING(VD BOUND_GUESS, VD delta_Bound, VD eps);

    void PrintSolution();
    void DumpSolution(std::string filename);

};

#endif