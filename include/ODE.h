#ifndef ODE_H
#define ODE_H

#include "VTypes.h"

class ODE
{
private:
    int _DOF;
    double _X_BEGIN;
    double _X_END;
    
    VD _X;
    VVD _Y;
    VVD _dYdX;
    
    VD _BOUND_CONDITION;
    VD _BOUND_BEGIN;
    VD _BOUND_END;
    VB _BOUND_AT_BEGIN_Q;
    VB _BOUND_AT_END_Q;
    int _N_BOUND_BEGIN;
    int _N_BOUND_END;

public:
    ODE();
    ODE(int DOF);
    ODE(int DOF, double x_begin, double x_end);
    ~ODE();

    void SetXBoundary(double x_begin, double x_end);
    void SetYBoundary(VD BOUNDARY);
    void SetYBoundary(VD BOUND_BEGIN, VD BOUND_END, VB BOUND_BEGIN_Q, VB BOUND_END_Q);
    void SetYBoundaryGuess(VD BOUND_GUESS);

};


#endif