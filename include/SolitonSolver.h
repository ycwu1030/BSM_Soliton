/*
 * @Description  : 
 * @Author       : Yongcheng Wu
 * @Date         : 2020-01-28 12:19:40
 * @LastEditors  : Yongcheng Wu
 * @LastEditTime : 2020-01-30 13:48:33
 */
#ifndef SolitonSolver_H
#define SolitonSolver_H

#include "VTypes.h"
#include <string>

typedef double (*ScalarFunction)(VD Field_Point, void *param);
typedef VD (*dScalarFunction)(VD Field_Point, void *param);

// ! For simplicity, we can scale the field and position such that we can deal with dimensionless variables.
// ! However, the final physical results should be dimensionful.
// ! In order to be consistent, I should follow:
// ! 1. Only GridPositions, FieldValue are stored dimensionless.
// ! 2. Any other quantities are stored/calculated in usual physical dimension.
class SolitonSolver
{
private:
    const int _NSpaceDim = 1;
    const double _energy_rel_error = 1e-20;
    const double _LEFT = -25;
    const double _RIGHT = 25;
    const double _del_t = 0.05;

    int _NFieldDim;
    int _MAXROUNDS;
    int _NGrid;
    int _NPoint; // _NGrid + 1;
    double _Scaling; // field_hat = field/_Scaling;  z_hat = v*z;
    VD _GridPositions_NoDim; // Dim: NPoint, dimensionless
    VVD _FieldValue_NoDim_Last; // Dim: NPoint x NFieldDim, dimensionless
    VVD _FieldValue_NoDim_Current; // Dim: NPoint x NFieldDim, dimensionless
    VD _bound_left; // Dimensional
    VD _bound_right; // Dimensional
    VD _EnergyDensity; // Dim: NGrid, Dimensional 
    double _TotalEnergy_Last;
    double _TotalEnergy_Current;
    double _DeltaZ;
    double _DeltaT;

    ScalarFunction _potential = nullptr;
    dScalarFunction _dpotential = nullptr;
    void *_param = nullptr;
    double _V0_global;

    void _Initiative();
    void _Iterating();
    double _Calc_EnergyDensity(VD point_beg, VD point_end);

public:
    SolitonSolver(){};
    SolitonSolver(int NFieldDim, int NGrid = 500, int MAXROUNDS = 1e7);
    ~SolitonSolver(){};

    void SetDimension(int NFieldDim);
    void SetGridPoints(int NGrid = 500, int MAXROUNDS = 1e7);

    void SetPotentials(ScalarFunction potential, dScalarFunction dpotential);
    void SetBoundary(VD left, VD right);
    void SetParam(void *param);
    void SetScale(double scaling = 246.0);
    void SetV0Global(double V0global);

    bool Solve();

    void PrintSolitonSolution();
    void DumpSolitonSolution(std::string filename);

};

#endif