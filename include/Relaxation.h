/*
 * @Description  : This class provides the APIs for Relaxation method in solving two-point boundary ODE
 * @Author       : Yongcheng Wu
 * @Date         : 2020-02-03 17:23:33
 * @LastEditors  : Yongcheng Wu
 * @LastEditTime : 2020-02-04 13:30:13
 */
#ifndef Relaxation_H
#define Relaxation_H

#include "VTypes.h"

class Relaxation
{
private:
    int _MeshPoints;// Number of mesh points
    int _NE; // Number of equations, please reduce the ODEs into the standard form with only first directives.
    int _NField; // The number of fields, _NE = 2*_NField
    int _NB_Left; // The number of boundary condition at the left hand side
    int _ITE_MAX; // Maximum number of iterations
    double _conv; // Convergence criterion
    VD _scales; // Scales for each variables.
    VVD _S;
    VVD _Y;
    VVVD _C;
    VD _indexV;

public:
    Relaxation();
    Relaxation(double NField, int NB_Left);
    ~Relaxation(){};

    void SetGrid(int MeshPoints=100);
    void SetMaxIteration(int itemax=100);
    void SetDOF(int NField, int NB_Left);
    void SetConvergeCriterion(double conv=1e-6);
    void SetBoundaryCondition();
    
};

#endif