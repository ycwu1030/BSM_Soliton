#ifndef Relaxation_H
#define Relaxation_H

#include <iostream>
#include <string>

#include "ODE.h"
#include "VTypes.h"

/*
 * About the Relaxation method:
 * For more detail, please consult the 'Numerical Recipes'
 * In short, the whole space (from x_begin to x_end), is divided into mesh points (in total _N_MeshPoints, the first one
 * is at x_begin, the last one is at x_end) For a d.o.f _DOF ODE, we will have _DOF variables. Under this context,
 * solving the ODE means we want to find out all tabulated values of the variables at each mesh point, Thus, in total we
 * need to find out _N_MeshPoints x _DOF variables. We achieve by iterating from initial guess. By each iteration, we
 * improve the guess a little bit. With the ODE and the boundary condition (some at x_begin, some at x_end), we then
 * have a linear function like: S * dV = F; where S is a matrix determined by the ODE and boundary condition F is the
 * difference between current guess and the final prospects dV is solved from above equation, and is used to update our
 * guess. Since, we convert ODE to FDE which links two neighbour points, so the form of matrix S is special:
 *
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
 * Above is an example for DOF = 5, N_MESH = 4, and 3 boundary condition at the beginning, 2 boundary condition at the
 * end What we want to do is reduce above form by Guassian Elimination into:
 *
 * C * dV = C'
 *
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
 *
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
 * For C: (N_MeshPoints + 1) x DOF x (DOF - NB + 1)  (Need to store all coefficients for all points)
 */

namespace BSM_Soliton {
struct MeshPoint {
    MeshPoint(int dof) : Y(dof, 0), Result(dof, 0), dResultdY(dof, VD(dof, 0)) {}
    double X;
    VD Y;

    VD Result;
    VVD dResultdY;
};

class RelaxationODE {
public:
    RelaxationODE(unsigned dof, unsigned left_boundary_size);
    virtual ~RelaxationODE();

    unsigned Get_DOF() const { return DOF; }

    virtual void dYdX(MeshPoint &point) = 0;

    unsigned Get_Left_Boundary_Size() const { return Left_Boundary_Size; }
    virtual void Left_Boundary_Constraints(MeshPoint &point) = 0;

    unsigned Get_Right_Boundary_Size() const { return Right_Boundary_Size; }
    virtual void Right_Boundary_Constraints(MeshPoint &point) = 0;

protected:
    const unsigned DOF;
    const unsigned Left_Boundary_Size;
    const unsigned Right_Boundary_Size;
};

// * Struct used to calculate and store current S matrix
struct RelaxationMatrix {
    RelaxationMatrix(int dof);
    VVD S;

    void Calc_S_at_Left_Boundary(RelaxationODE *ode, MeshPoint &p_left);
    void Calc_S_at_Middle(RelaxationODE *ode, MeshPoint &p_km1, MeshPoint &p_k);
    void Calc_S_at_Right_Boundary(RelaxationODE *ode, MeshPoint &p_right);
};

class Relaxation {
public:
    typedef std::vector<MeshPoint> Grid;
    explicit Relaxation(RelaxationODE *fode, double rel_error_threshold = 0.5, double converge_criteria = 1e-6,
                        int MeshSize = 400);

    void Set_Mesh_Size(int MeshSize = 400);
    bool Solve(const VD &x, const VVD &y);  // * Solve the ODE with relaxation method using x, y as initial guess.
    void DumpSolution(std::string filename);

private:
    RelaxationODE *ode;
    const int DOF;
    const int Left_Boundary_Size;
    int Mesh_Size;
    double rel_error_threshold;
    double converge_criteria;

    Grid mesh_grid;
    VVVD Relax_C;  // (Mesh_Size+1)*DOF*(DOF - left_boundary_size + 1);
    RelaxationMatrix Relax_S;

    void Reduce_to_Zero(int mesh_id);
    void Pivot_Elimination(int mesh_id);
    void Backsubstitution();
    void Relax();  // relax the solution by one step
    double Update_Grid();
};

}  // namespace BSM_Soliton

std::ostream &operator<<(std::ostream &os, const BSM_Soliton::MeshPoint &point);

struct Relaxation_Param {
    double x1;
    double x2;
    VD y1;
    VD y2;

    int k_init;   // starting meshpoint,  k==k_init means we need first boundary condition
    int k_final;  // ending meshpoint,  k > k_final means we need the second boundary condition

    int k;  // Current meshpoint (if k = k_init, deal with the first boundary, if k > k_final, deal with the final
            // boundary, otherwise, deal with the FDE linking k-1 and k)

    // Used for filling S matrix:
    int k_coeff;          // This is the column index where the difference should be filled
    int index_row_begin;  // This is the starting row index in S need to be filled, This will also be passed to reduce
                          // and pinvs
    int index_row_end;    // This is the last row index in S need to be filled, This will also be passed to reduce and
                          // pinvs

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
    int index_zero_begin;  // The first column need to be set to zero;
    int index_zero_end;    // The last column need to be set to zero;
    int index_mod_begin;   // The first column need to be modified in order to zero leading coeffs, this is also used by
                           // pinvs
    int index_mod_end;     // The last column need to be modified in order to zero leading coeffs

    // Used for diagonalize the the sub-matrix of S
    int index_column_off;  // The column index of C at which start to store the reduced coeffs; It should be DOF-NB for
                           // last block (the final boundary condition), but 0 for all other cases

    VD indexv;
};

typedef void (*DIFEQ)(const Relaxation_Param relax_param, void *param, VVD &S);

class Relaxation {
private:
    int _N_MeshPoints;  // Number of mesh points
    int _DOF;  // Number of equations, please reduce the ODEs into the standard form with only first directives.
    int _N_BOUND_LEFT;   // The number of boundary condition at the left hand side
    int _N_BOUND_RIGHT;  // = _DOF - _N_BOUND_LEFT
    int _ITE_MAX;        // Maximum number of iterations
    double _conv;        // Convergence criterion
    double _slowc;

    double _X_BEGIN;
    double _X_END;
    VD _Y_BEGIN;
    VD _Y_END;

    VD _scales;  // Scales for each variables.
    VD _X;
    VVD _Y;

    VD _X_Guess;
    VVD _Y_Guess;

    VVD _S;   // _DOF * (2*_DOF+1)
    VVVD _C;  // (_N_MeshPoints+1) * _DOF * (_DOF-_N_BOUND_LEFT+1);

    VD _indexV;

    Relaxation_Param _relax_param;
    DIFEQ _difeq;
    void *_param;

    void _reduce();
    void _pinvs();
    void _bksub();

    void _INIT();

public:
    Relaxation();
    Relaxation(int DOF, int NB_Left, int MeshPoints = 100);
    ~Relaxation(){};

    void SetDOF(int DOF, int NB_Left, int MeshPoints = 100);
    void SetBoundary(double x_begin, double x_end, VD y_begin,
                     VD y_end);  // The inputs help setting the inital value for relaxation
    void SetBoundary(VD X_Guess, VVD Y_Guess) {
        _X_Guess = X_Guess;
        _Y_Guess = Y_Guess;
    }  // Directly set the guess solution
    void SetMaxIteration(int itemax = 100);
    void SetScales(VD scales) { _scales = scales; }
    void SetConvergeCriterion(double slowc = 1.0, double conv = 1e-6);
    void SetODESystem(DIFEQ difeq, void *param);

    bool SOLVDE();

    void PrintS();
    void PrintC();
    void PrintC(int k);
    void PrintSolution();
    void DumpSolution(std::string filename);

    VVD GetY() { return _Y; }
    VD GetX() { return _X; }
};

#endif
