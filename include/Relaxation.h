/*
 * @Description  : This class provides the APIs for Relaxation method in solving two-point boundary ODE
 * @Author       : Yongcheng Wu
 * @Date         : 2020-02-03 17:23:33
 * @LastEditors  : Yongcheng Wu
 * @LastEditTime : 2020-02-09 10:48:42
 */
#ifndef Relaxation_H
#define Relaxation_H

#include "VTypes.h"
#include <string>

/*
 * About the Relaxation method:
 * For more detail, please consult the 'Numerical Recipes'
 * In short, the whole space (from x_begin to x_end), is divided into mesh points (in total _N_MeshPoints, the first one is at x_begin, the last one is at x_end)
 * For a d.o.f _DOF ODE, we will have _DOF variables.
 * Under this context, solving the ODE means we want to find out all tabulated values of the variables at each mesh point,
 * Thus, in total we need to find out _N_MeshPoints x _DOF variables.
 * We achieve by iterating from initial guess. By each iteration, we improve the guess a little bit.
 * With the ODE and the boundary condition (some at x_begin, some at x_end), we then have a linear function like:
 * S * dV = F;
 * where S is a matrix determined by the ODE and boundary condition
 * F is the difference between current guess and the final prospects
 * dV is solved from above equation, and is used to update our guess.
 * Since, we convert ODE to FDE which links two neighbour points, so the form of matrix S is special:
 * X X X X X                                V  F
 * X X X X X  <- Initial Block              V  F
 * X X X X X                                V  F
 * X X X X X X X X X X                      V  F
 * X X X X X X X X X X                      V  F
 * X X X X X X X X X X                      V  F
 * X X X X X X X X X X                      V  F
 * X X X X X X X X X X                      V  F
 *           X X X X X X X X X X            V  F
 *           X X X X X X X X X X            V  F
 *           X X X X X X X X X X            V  F
 *           X X X X X X X X X X            V  F
 *           X X X X X X X X X X            V  F
 *                     X X X X X X X X X X  V  F
 *                     X X X X X X X X X X  V  F
 *                     X X X X X X X X X X  V  F
 *                     X X X X X X X X X X  V  F
 *                     X X X X X X X X X X  V  F
 *           Final Block ->      X X X X X  V  F
 *                               X X X X X  V  F
 * Above is an example for DOF = 5, N_MESH = 4, and 3 boundary condition at the beginning, 2 boundary condition at the end
 * What we want to do is reduce above form by Guassian Elimination into:
 * C * dV = C'
 * 1     X X                                V  F
 *   1   X X                                V  F
 *     1 X X                                V  F
 *       1         X X                      V  F
 *         1       X X                      V  F
 *           1     X X                      V  F
 *             1   X X                      V  F
 *               1 X X                      V  F
 *                 1         X X            V  F
 *                   1       X X            V  F
 *                     1     X X            V  F
 *                       1   X X            V  F
 *                         1 X X            V  F
 *                           1         X X  V  F
 *                             1       X X  V  F
 *                               1     X X  V  F
 *                                 1   X X  V  F
 *                                   1 X X  V  F
 *                                     1    V  F
 *                                       1  V  F
 * The detailed reduction procedure can be found in 'Numerical recipes'
 * Here we only list the way we store the S matrix, and the final form C
 * For S:
 * Beside the initial and final block, the sub-matrix in S is the same: DOF x 2DOF;
 * We also need to store the difference F, thus the matrix used to store S is DOF x (2*DOF+1)
 * Note that we do the reduction point by point, so we don't need a whole matrix 
 * But for the initial and final block, they don't span full DOF x 2 DOF matrix, 
 * but they still be put into such matrix by following convention:
 * Initial block:
 * - - - - - - - - - - -
 * - - - - - - - - - - -
 * - - - - - X X X X X F
 * - - - - - X X X X X F
 * - - - - - X X X X X F
 * Final block:
 * - - - - - X X X X X F
 * - - - - - X X X X X F
 * - - - - - - - - - - -
 * - - - - - - - - - - -
 * - - - - - - - - - - -
 * 
 * For C:
 * After reduction, only few of the coefficients need to be stored
 * According the above reduced form, we find that for each point, we need just a matrix of DOF x (DOF - NB)
 * with final difference, we need DOF x (DOF - NB + 1). Here NB is the number of boundary condition at the initial point
 * C C F
 * C C F
 * C C F
 * C C F
 * C C F
 * Again, for initial block
 * - - -
 * - - -
 * C C F
 * C C F
 * C C F
 * For final block
 * - - F
 * - - F
 * - - -
 * - - -
 * - - -
 * 
 * In total, the storage is
 * For S: just DOF x (2*DOF + 1)   (Just needed for one point)
 * For C: N_MeshPoints x DOF x (DOF - NB + 1)  (Need to store all coefficients for all points)
*/

struct Relaxation_Param
{
    double x1; 
    double x2;
    VD y1;
    VD y2;

    int k_init; // starting meshpoint,  k==k_init means we need first boundary condition
    int k_final; // ending meshpoint,  k > k_final means we need the second boundary condition

    int k; // Current meshpoint (if k = k_init, deal with the first boundary, if k > k_final, deal with the final boundary, otherwise, deal with the FDE linking k-1 and k)

    // Used for filling S matrix:
    int k_coeff; // This is the column index where the difference should be filled
    int index_row_begin; // This is the starting row index in S need to be filled, This will also be passed to reduce and pinvs
    int index_row_end; // This is the last row index in S need to be filled, This will also be passed to reduce and pinvs

    // Used for reducing leading part of S matrix to 0:
    // Reducing a matrix means from:
    // 1 0 0 S S 
    // 0 1 0 S S
    // 0 0 1 S S
    // Z Z Z D D D D D A A
    // Z Z Z D D D D D A A
    // Z Z Z D D D D D A A
    // Z Z Z D D D D D A A
    // Z Z Z D D D D D A A
    // to
    // 1 0 0 S S 
    // 0 1 0 S S
    // 0 0 1 S S
    // 0 0 0 d d D D D A A
    // 0 0 0 d d D D D A A
    // 0 0 0 d d D D D A A
    // 0 0 0 d d D D D A A
    // 0 0 0 d d D D D A A
    int index_zero_begin; // The first column need to be set to zero;
    int index_zero_end; // The last column need to be set to zero;
    int index_mod_begin; // The first column need to be modified in order to zero leading coeffs, this is also used by pinvs
    int index_mod_end; // The last column need to be modified in order to zero leading coeffs


    // Used for diagonalize the the sub-matrix of S
    int index_column_off; // The column index of C at which start to store the reduced coeffs; It should be DOF-NB for last block (the final boundary condition), but 0 for all other cases


    VD indexv;

};


typedef void (*DIFEQ)(const Relaxation_Param relax_param, void *param, VVD &S);

class Relaxation
{
private:
    int _N_MeshPoints;// Number of mesh points
    int _DOF; // Number of equations, please reduce the ODEs into the standard form with only first directives.
    int _N_BOUND_LEFT; // The number of boundary condition at the left hand side
    int _N_BOUND_RIGHT; // = _DOF - _N_BOUND_LEFT
    int _ITE_MAX; // Maximum number of iterations
    double _conv; // Convergence criterion
    double _slowc;
    
    double _X_BEGIN;
    double _X_END;
    VD _Y_BEGIN;
    VD _Y_END;

    VD _scales; // Scales for each variables.
    VD _X;
    VVD _Y;

    VD _X_Guess;
    VVD _Y_Guess;

    VVD _S; // _DOF * (2*_DOF+1)
    VVVD _C; // (_N_MeshPoints+1) * _DOF * (_DOF-_N_BOUND_LEFT+1);
    
    VD _indexV;

    Relaxation_Param _relax_param;
    DIFEQ _difeq;
    void *_param;

    void _reduce();
    void _pinvs();
    void _bksub();
    
    
    

public:
    Relaxation();
    Relaxation(int DOF, int NB_Left, int MeshPoints=100);
    ~Relaxation(){};

    void SetDOF(int DOF, int NB_Left, int MeshPoints=100);
    void SetBoundary(double x_begin, double x_end, VD y_begin, VD y_end);// The inputs help setting the inital value for relaxation
    void SetBoundary(VD X_Guess, VVD Y_Guess){_X_Guess = X_Guess; _Y_Guess=Y_Guess;} // Directly set the guess solution
    void SetMaxIteration(int itemax=100);
    void SetScales(VD scales){_scales = scales;}
    void SetConvergeCriterion(double slowc=1.0, double conv=1e-6);
    void SetODESystem(DIFEQ difeq, void *param);
    void _INIT();
    void _SOLVDE();

    void PrintS();
    void PrintC();
    void PrintC(int k);
    void PrintSolution();
    void DumpSolution(std::string filename);
};

#endif