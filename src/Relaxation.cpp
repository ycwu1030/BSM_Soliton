#include "Relaxation.h"

#include <cmath>
#include <fstream>
#include <iostream>
using namespace std;

ostream &operator<<(ostream &os, const BSM_Soliton::MeshPoint &point) {
    os << point.X << "\t";
    os << point.Y;
    return os;
}

namespace BSM_Soliton {

void PrintS(VVD &S) {
    cout << "---------------------->" << endl;
    for (int i = 0; i < S.size(); i++) {
        for (int j = 0; j < S[i].size(); j++) {
            cout << S[i][j] << "\t";
        }
        cout << endl;
    }
    cout << "<----------------------" << endl;
}

RelaxationODE::RelaxationODE(unsigned dof, unsigned left_boundary_size)
    : DOF(dof),
      Left_Boundary_Size(left_boundary_size <= dof ? left_boundary_size : dof),
      Right_Boundary_Size(left_boundary_size <= dof ? dof - left_boundary_size : 0) {}

RelaxationODE::~RelaxationODE() {}

RelaxationMatrix::RelaxationMatrix(int dof) : DOF(dof), S(dof, VD(2 * dof + 1, 0)) {}

void RelaxationMatrix::Clean() {
    S.clear();
    S.resize(DOF, VD(2 * DOF + 1, 0));
}

void RelaxationMatrix::Calc_S_at_Left_Boundary(RelaxationODE *ode, MeshPoint &p_left) {
    /*
     * Convention for Left Boundary:
     * - - - - - - - - - - -
     * - - - - - - - - - - -
     * - - - - - X X X X X F
     * - - - - - X X X X X F
     * - - - - - X X X X X F
     */
    Clean();
    int dof = ode->Get_DOF();
    int left_boundary_size = ode->Get_Left_Boundary_Size();
    int offset = dof - left_boundary_size;
    ode->Left_Boundary_Constraints(p_left);
    for (int i = 0; i < left_boundary_size; i++) {
        for (int j = 0; j < dof; j++) {
            S[i + offset][j + dof] = p_left.dResultdY[i][j];
        }
        S[i + offset][2 * dof] = -p_left.Result[i];
    }
}

void RelaxationMatrix::Calc_S_at_Right_Boundary(RelaxationODE *ode, MeshPoint &p_right) {
    /*
     * Convention for Right Boundary:
     * - - - - - X X X X X F
     * - - - - - X X X X X F
     * - - - - - - - - - - -
     * - - - - - - - - - - -
     * - - - - - - - - - - -
     */
    Clean();
    int dof = ode->Get_DOF();
    int right_boundary_size = ode->Get_Right_Boundary_Size();
    ode->Right_Boundary_Constraints(p_right);
    for (int i = 0; i < right_boundary_size; i++) {
        for (int j = 0; j < dof; j++) {
            S[i][j + dof] = p_right.dResultdY[i][j];
        }
        S[i][2 * dof] = -p_right.Result[i];
    }
}

void RelaxationMatrix::Calc_S_at_Middle(RelaxationODE *ode, MeshPoint &p_km1, MeshPoint &p_k) {
    Clean();
    int dof = ode->Get_DOF();
    MeshPoint p_mid(dof);
    double xm1 = p_km1.X;
    double x = p_k.X;
    VD &ym1 = p_km1.Y;
    VD &y = p_k.Y;
    p_mid.X = (xm1 + x) / 2.0;
    p_mid.Y = (ym1 + y) / 2.0;
    ode->dYdX(p_mid);
    for (int i = 0; i < dof; i++) {
        for (int j = 0; j < dof; j++) {
            S[i][j] = -(x - xm1) * p_mid.dResultdY[i][j] / 2.0;
            S[i][j + dof] = S[i][j];
            if (i == j) {
                S[i][j] += -1.0;
                S[i][j + dof] += 1.0;
            }
        }
        S[i][2 * dof] = -(y[i] - ym1[i] - (x - xm1) * p_mid.Result[i]);
    }
}

Relaxation::Relaxation(RelaxationODE *fode, double threshold, double converge)
    : ode(fode),
      //   Mesh_Size(MeshSize),
      DOF(fode->Get_DOF()),
      Left_Boundary_Size(fode->Get_Left_Boundary_Size()),
      Relax_S(fode->Get_DOF()),
      ITER_MAX(1e4),
      //   Relax_C(MeshSize + 1, VVD(fode->Get_DOF(), VD(fode->Get_Right_Boundary_Size() + 1, 0))),
      //   mesh_grid(fode->Get_DOF(), MeshSize),
      rel_error_threshold(threshold),
      converge_criteria(converge) {}

void Relaxation::Set_Mesh_Size(int mesh_size) {
    Mesh_Size = mesh_size;
    mesh_grid.clear();
    mesh_grid.resize(Mesh_Size, MeshPoint(DOF));
    Relax_C.clear();
    Relax_C.resize(Mesh_Size + 1, VVD(DOF, VD(ode->Get_Right_Boundary_Size() + 1, 0)));
}

void Reduce_to_Zero(const unsigned s_row_dim, const unsigned offset, const unsigned s_row_max, const unsigned s_col_beg,
                    const VVD &c, VVD &s) {
    // * Reduce the left-upper most s_row_max*(s_row_dim-offset) sub-matrix of S to zero
    // * according to following convention;
    // * From
    // *
    // *     offset
    // *       ^
    // *     |   |
    // * 1 0 0 0 0 C C              V C  <-\
    // * 0 1 0 0 0 C C              V C  <-\
    // * 0 0 1 0 0 C C              V C  <-\
    // * 0 0 0 1 0 C C              V C  <--+-- Stored in Relax_C
    // * 0 0 0 0 1 C C              V C  <-/
    // * 0 0 Z Z Z D D D D D A A    V S  <-\
    // * 0 0 Z Z Z D D D D D A A    V S  <-\
    // * 0 0 Z Z Z D D D D D A A    V S  <--+-- Current Relax_S
    // * 0 0 Z Z Z D D D D D A A    V S  <-/
    // * 0 0 Z Z Z D D D D D A A    V S  <-/
    // * To
    // * 1 0 0 0 0 C C              V C  <-\
    // * 0 1 0 0 0 C C              V C  <-\
    // * 0 0 1 0 0 C C              V C  <-\
    // * 0 0 0 1 0 C C              V C  <--+-- Stored in Relax_C
    // * 0 0 0 0 1 C C              V C  <-/
    // * 0 0 0 0 0 E E D D D A A    V R  <-\
    // * 0 0 0 0 0 E E D D D A A    V R  <-\
    // * 0 0 0 0 0 E E D D D A A    V R  <--+-- Current Relax_S
    // * 0 0 0 0 0 E E D D D A A    V R  <-/
    // * 0 0 0 0 0 E E D D D A A    V R  <-/
    // *
    // * ir = row-index
    // * ic = col-index
    unsigned ir_s_beg = 0;
    unsigned ir_s_end = s_row_max > s_row_dim ? s_row_dim : s_row_max;
    unsigned ic_s_to_be_zero_beg = 0;
    unsigned ic_s_to_be_zero_end = offset;
    unsigned ir_c_offset = s_row_dim - offset;
    for (unsigned ir_s = ir_s_beg; ir_s < ir_s_end; ir_s++) {
        // * Go through the S matrix line by line;
        for (unsigned ic_s = ic_s_to_be_zero_beg; ic_s < ic_s_to_be_zero_end; ic_s++) {
            // * Go through the corresponding Z -> 0 position in current line
            double Z = s[ir_s][ic_s + s_col_beg];
            for (unsigned ic_c = 0; ic_c < s_row_dim - offset; ic_c++) {
                // * Go through the corrsponding D->E position
                double C = c[ic_s + ir_c_offset][ic_c];
                s[ir_s][ic_c + s_col_beg + ic_s_to_be_zero_end] -= Z * C;
            }
            double C = c[ic_s + ir_c_offset][s_row_dim - offset];
            s[ir_s][2 * s_row_dim] -= Z * C;
        }
    }
}

void Gaussian_Elimination_with_Partial_Pivot(const unsigned r_beg, const unsigned c_beg, const unsigned dim, VVD &s,
                                             VVD &c) {
    // * Reduce the submatrix of s of dimension dim and starting from (r_beg, c_beg) to identity matrix
    // * Other element of s matrix will change simultaneously
    // * Gaussian elimination is used with partial pivot
    // * S matrix is
    // *
    // * S S S S S S S S S S    S
    // * S S S S S S S S S S    S
    // * S S S S S S S S S S    S
    // * S S S S S S S S S S    S
    // * S S S S S S S S S S    S
    // *
    // * Reduce to
    // *
    // * R R R R R R R R R R    R
    // * R R R 1 0 0 R R R R    R  <- r_beg  --\
    // * R R R 0 1 0 R R R R    R               > dim
    // * R R R 0 0 1 R R R R    R            --/
    // * R R R R R R R R R R    R
    // *       ^
    // *       |
    // *       c_beg

    vector<int> row_permutation_index(dim, -1);  // Used to store the perumutation information for rows
    VD row_scale(dim, 0);                        // The scale in each row

    // Record the biggest element in each row to use as scale for each row;

    unsigned ir_beg = r_beg;
    unsigned ir_end = r_beg + dim;
    unsigned ic_beg = c_beg;
    unsigned ic_end = ic_beg + dim;
    unsigned i_coeff = s[0].size();
    for (size_t i = ir_beg; i < ir_end; i++) {
        double largest_in_row = 0;
        for (size_t j = ic_beg; j < ic_end; j++) {
            if (abs(s[i][j]) > largest_in_row) {
                largest_in_row = abs(s[i][j]);
            }
        }
        if (largest_in_row == 0) {
            cout << "Matrix s has one row with all zeros" << endl;
            return;
        }
        row_scale[i - ir_beg] = 1.0 / largest_in_row;
        row_permutation_index[i - ir_beg] = -1;
    }
    // cout << "scale for each row: " << row_scale << endl;

    double pivot;
    double pivot_inverse;
    double dum;
    int ic_largest;
    int ir_pivot, ic_pivot;

    for (size_t id = 0; id < dim; id++) {
        pivot = 0.0;

        // * Finding the current pivot (the largest element)
        for (size_t i = ir_beg; i < ir_end; i++) {
            // * Do not consider the row that already contains pivot
            if (row_permutation_index[i - ir_beg] >= 0) continue;
            double largest = 0.0;
            for (size_t j = ic_beg; j < ic_end; j++) {
                if (abs(s[i][j]) > largest) {
                    ic_largest = j;
                    largest = abs(s[i][j]);
                }
            }
            if (largest * row_scale[i - ir_beg] > pivot) {
                ir_pivot = i;
                ic_pivot = ic_largest;
                pivot = largest * row_scale[i - ir_beg];
            }
        }
        // PrintS(s);
        // cout << "pivot-" << id << ": " << s[ir_pivot][ic_pivot] << endl;
        // * Now current pivot is at (ir_pivot,ic_pivot)
        if (s[ir_pivot][ic_pivot] == 0.0) {
            cout << "The whole matrix is zero" << endl;
            return;
        }
        // * Then this row should become ic_pivot-th row,
        // * But the ir_pivot and ic_pivot are starting from ir_beg and ic_beg respectively
        // * We need to compensate that
        row_permutation_index[ir_pivot - ir_beg] = ic_pivot - ic_beg;

        // * Then rescale the row where the pivot is in, such that pivot will be 1;
        // * This is specific for relaxation application, so we won't scale the elements before the starting column
        pivot_inverse = 1.0 / s[ir_pivot][ic_pivot];
        for (size_t j = ic_beg; j < s[ir_pivot].size(); j++) {
            s[ir_pivot][j] *= pivot_inverse;
        }
        // * The pivot is set explicitly to 1
        s[ir_pivot][ic_pivot] = 1.0;

        // * Then we need to eliminate the elements in the same column as current pivot
        // * The other elements will change simultaneously
        for (size_t i = ir_beg; i < ir_end; i++) {
            // * Do not consider the row where current pivot is
            if (i == ir_pivot) continue;
            if (s[i][ic_pivot]) {
                dum = s[i][ic_pivot];
                for (size_t j = ic_beg; j < s[i].size(); j++) {
                    s[i][j] -= dum * s[ir_pivot][j];
                }
                s[i][ic_pivot] = 0.0;
            }
        }
    }

    // * Now store the part after the identity matrix into c
    int ic_s_coeff = s[0].size() - 1;
    int ic_c_coeff = c[0].size() - 1;
    for (size_t i = ir_beg; i < ir_end; i++) {
        // * the row index in c is determined from row_permutation_index,
        // * but will be shifted according to ir_beg.
        int ir_c = row_permutation_index[i - ir_beg] + ir_beg;
        for (size_t j = ic_end; j < ic_s_coeff; j++) {
            c[ir_c][j - ic_end] = s[i][j];
        }
        c[ir_c][ic_c_coeff] = s[i][ic_s_coeff];
    }
}

void Backsubstitution(const unsigned dimension, const unsigned offset, const unsigned mesh_size, VVVD &c) {
    /*
     * 1     X X                                V  F  -\
     *   1   X X                                V  F  -+  offset   - 0th block in c
     *     1 X X                                V  F  -/
     *       1         X X                      V  F   -\
     *         1       X X                      V  F   -\
     *           1     X X                      V  F   -+   dimension  - 1th block in c
     *             1   X X                      V  F   -/
     *               1 X X                      V  F   -/
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
     *                                     1    V  F  -\
     *                                       1  V  F  -/  Last block in c
     *
     * Find the solutions to those V's
     * only X's and F's are stored inside c
     */
    // * Backsubstitution to find the solution
    int ic_coeff = dimension - offset;
    for (int mesh_id = mesh_size - 1; mesh_id >= 0; --mesh_id) {
        // * The storage in first block in c does not start from 0;
        int ir_beg_in_c = mesh_id == 0 ? dimension - offset : 0;
        for (int ic_pre = 0; ic_pre < dimension - offset; ic_pre++) {
            double res = c[mesh_id + 1][ic_pre][ic_coeff];
            for (int ir_c = ir_beg_in_c; ir_c < dimension; ir_c++) {
                c[mesh_id][ir_c][ic_coeff] -= c[mesh_id][ir_c][ic_pre] * res;
            }
        }
    }

    // * Then move all the results into the first column
    // * And arrange according to the mesh point
    for (int mesh_id = 0; mesh_id < mesh_size; mesh_id++) {
        for (int ir = 0; ir < offset; ir++) {
            c[mesh_id][ir][0] = c[mesh_id][ir + dimension - offset][ic_coeff];
        }
        for (int ir = 0; ir < dimension - offset; ir++) {
            c[mesh_id][ir + offset][0] = c[mesh_id + 1][ir][ic_coeff];
        }
    }
    // * Now the results are stored in c mesh-point by mesh-point in the first column
}

void Relaxation::Reduce_to_Zero(int mesh_id) {
    if (mesh_id <= 0) return;  // * For the S matrix corresponding to left-boundary, we don't need to do this
    unsigned row_max = mesh_id >= Mesh_Size ? DOF - Left_Boundary_Size : DOF;
    unsigned col_beg = mesh_id >= Mesh_Size ? DOF : 0;
    BSM_Soliton::Reduce_to_Zero(DOF, Left_Boundary_Size, row_max, col_beg, Relax_C[mesh_id - 1], Relax_S.S);
}

void Relaxation::Pivot_Elimination(int mesh_id) {
    unsigned row_begin = mesh_id <= 0 ? DOF - Left_Boundary_Size : 0;
    unsigned col_begin = mesh_id <= 0 ? DOF : (mesh_id >= Mesh_Size ? DOF + Left_Boundary_Size : Left_Boundary_Size);
    unsigned dimension = mesh_id <= 0 ? Left_Boundary_Size : (mesh_id >= Mesh_Size ? DOF - Left_Boundary_Size : DOF);
    BSM_Soliton::Gaussian_Elimination_with_Partial_Pivot(row_begin, col_begin, dimension, Relax_S.S, Relax_C[mesh_id]);
}

void Relaxation::Backsubstitution() { BSM_Soliton::Backsubstitution(DOF, Left_Boundary_Size, Mesh_Size, Relax_C); }

void Relaxation::Relax() {
    // * Relax
    // * Left Boundary
    Relax_S.Calc_S_at_Left_Boundary(ode, mesh_grid[0]);
    Reduce_to_Zero(0);
    Pivot_Elimination(0);

    // * Mesh Point
    for (size_t mesh_id = 1; mesh_id < Mesh_Size; mesh_id++) {
        Relax_S.Calc_S_at_Middle(ode, mesh_grid[mesh_id - 1], mesh_grid[mesh_id]);
        Reduce_to_Zero(mesh_id);
        Pivot_Elimination(mesh_id);
    }
    // * Right Boundary
    Relax_S.Calc_S_at_Right_Boundary(ode, mesh_grid[Mesh_Size - 1]);
    Reduce_to_Zero(Mesh_Size);
    Pivot_Elimination(Mesh_Size);

    // * Backsubstitution
    Backsubstitution();
}

double Relaxation::Update_Grid() {
    // * Using the results in Relax_C to update the results in grid

    // * First determine the error
    double rel_error_total = 0;
    for (size_t index_y = 0; index_y < DOF; index_y++) {
        for (size_t mesh_id = 0; mesh_id < Mesh_Size; mesh_id++) {
            double dy_cur = abs(Relax_C[mesh_id][index_y][0]);
            double y_cur = abs(mesh_grid[mesh_id].Y[index_y]);
            if (y_cur < 1) y_cur = 1;
            rel_error_total += dy_cur / y_cur;
        }
    }
    double rel_error_average = rel_error_total / DOF / Mesh_Size;
    double frac = rel_error_average > rel_error_threshold ? rel_error_threshold / rel_error_average : 1;

    // * Update the grid
    for (size_t index_y = 0; index_y < DOF; index_y++) {
        for (size_t mesh_id = 0; mesh_id < Mesh_Size; mesh_id++) {
            mesh_grid[mesh_id].Y[index_y] += Relax_C[mesh_id][index_y][0] * frac;
        }
    }
    return rel_error_average;
}

bool Relaxation::Solve(const VD &x, const VVD &y) {
    if (x.size() != y.size()) {
        cout << "The input grid size does not match between x and y" << endl;
        cout << "x.size() = " << x.size() << " != "
             << " y.size() = " << y.size() << endl;
        return false;
    }
    if (y[0].size() != DOF) {
        cout << "The input DOF of y does not match with the ode" << endl;
        cout << " dof of y = " << y[0].size() << " != "
             << " dof of ode = " << DOF << endl;
    }
    Set_Mesh_Size(x.size());
    for (size_t mesh_id = 0; mesh_id < Mesh_Size; mesh_id++) {
        mesh_grid[mesh_id].X = x[mesh_id];
        mesh_grid[mesh_id].Y = y[mesh_id];
    }
    size_t iter = 0;
    // size_t ITER_MAX = 1000;
    for (; iter < ITER_MAX; iter++) {
        cout << "Iter-" << iter << endl;
        Relax();
        double err = Update_Grid();
        cout << "Error = " << err << endl;
        if (err < converge_criteria) {
            break;
        }
    }
    if (iter == ITER_MAX) return false;
    return true;
}

void Relaxation::DumpSolution(string filename) {
    ofstream output(filename.c_str());
    // output<<"The Solution is:"<<endl;
    output << "x\t";
    for (size_t i = 0; i < DOF; i++) {
        output << "y_" << i << "\t";
    }
    output << endl;
    for (size_t i = 0; i < Mesh_Size; i++) {
        output << mesh_grid[i] << endl;
    }
}

}  // namespace BSM_Soliton

Relaxation::Relaxation() {
    SetDOF(1, 1);
    SetMaxIteration();
    SetConvergeCriterion();
}
Relaxation::Relaxation(int DOF, int NB_Left, int MeshPoints) {
    SetDOF(DOF, NB_Left, MeshPoints);
    SetMaxIteration();
    SetConvergeCriterion();
}
void Relaxation::SetMaxIteration(int itemax) { _ITE_MAX = itemax; }
void Relaxation::SetDOF(int DOF, int NB_Left, int MeshPoints) {
    _DOF = DOF;
    _N_BOUND_LEFT = NB_Left;
    _N_BOUND_RIGHT = _DOF - _N_BOUND_LEFT;
    _N_MeshPoints = MeshPoints;
}
void Relaxation::SetConvergeCriterion(double slowc, double conv) {
    _slowc = slowc;
    _conv = conv;
}
void Relaxation::SetBoundary(double x_begin, double x_end, VD y_begin, VD y_end) {
    _X_BEGIN = x_begin;
    _X_END = x_end;
    _Y_BEGIN = y_begin;
    _Y_END = y_end;
}
void Relaxation::SetODESystem(DIFEQ difeq, void *param) {
    _difeq = difeq;
    _param = param;
}
void Relaxation::_INIT() {
    // ! Initiative the vectors
    _X.clear();
    _Y.clear();
    if (_X_Guess.size() != 0 && _Y_Guess.size() == _X_Guess.size()) {
        _X = _X_Guess;
        _Y = _Y_Guess;
    } else {
        double x_stepsize = (_X_END - _X_BEGIN) / (_N_MeshPoints - 1);
        VD y_stepsize = (_Y_END - _Y_BEGIN) / (_N_MeshPoints - 1);
        for (size_t i = 0; i < _N_MeshPoints; i++) {
            _X.push_back(_X_BEGIN + i * x_stepsize);
            _Y.push_back(_Y_BEGIN + i * y_stepsize);
        }
    }
    _S.clear();
    _C.clear();
    _S = VVD(_DOF, VD(2 * _DOF + 1, 0));
    _C = VVVD(_N_MeshPoints + 1, VVD(_DOF, VD(_DOF - _N_BOUND_LEFT + 1, 0)));
}
bool Relaxation::SOLVDE() {
    _INIT();
    _relax_param.k_init = 0;
    _relax_param.k_final = _N_MeshPoints - 1;

    int nvars = _DOF * _N_MeshPoints;
    vector<int> kmax(_DOF, 0);
    VD ermax(_DOF, 0);

    for (size_t iter = 0; iter < _ITE_MAX; iter++) {
        // if ((iter+1)%1000==0)
        // {
        //     cout<<"ITER: "<<iter+1<<endl;
        // }
        // The first boundary at k=k_init;
        _relax_param.k = _relax_param.k_init;

        _relax_param.x1 = _X[_relax_param.k];
        _relax_param.x2 = _X[_relax_param.k];
        _relax_param.y1 = _Y[_relax_param.k];
        _relax_param.y2 = _Y[_relax_param.k];

        _relax_param.k_coeff = 2 * _DOF;
        _relax_param.index_row_begin = _DOF - _N_BOUND_LEFT;
        _relax_param.index_row_end = _DOF - 1;
        _relax_param.index_mod_begin = _DOF;
        _relax_param.index_column_off = 0;

        _difeq(_relax_param, _param, _S);
        // cout<<"Iter: "<<iter<<"  "<<"S at left:"<<endl;
        // PrintS();
        _pinvs();

        // Intermidate points
        for (size_t k = _relax_param.k_init + 1; k <= _relax_param.k_final; k++) {
            _relax_param.k = k;

            _relax_param.x1 = _X[_relax_param.k - 1];
            _relax_param.x2 = _X[_relax_param.k];
            _relax_param.y1 = _Y[_relax_param.k - 1];
            _relax_param.y2 = _Y[_relax_param.k];

            _relax_param.index_row_begin = 0;
            _relax_param.index_row_end = _DOF - 1;
            _relax_param.index_zero_begin = 0;
            _relax_param.index_zero_end = _N_BOUND_LEFT - 1;
            _relax_param.index_mod_begin = _N_BOUND_LEFT;
            _relax_param.index_mod_end = _DOF - 1;
            _relax_param.k_coeff = 2 * _DOF;

            _relax_param.index_column_off = 0;

            _difeq(_relax_param, _param, _S);
            _reduce();
            _pinvs();
            // cout<<"Iter: "<<iter<<"  C at "<<k<<endl;
            // PrintC();
        }

        // The final boundary condition
        _relax_param.k = _relax_param.k_final + 1;

        _relax_param.x1 = _X[_relax_param.k - 1];
        _relax_param.x2 = _X[_relax_param.k - 1];
        _relax_param.y1 = _Y[_relax_param.k - 1];
        _relax_param.y2 = _Y[_relax_param.k - 1];

        _relax_param.index_row_begin = 0;
        _relax_param.index_row_end = _DOF - _N_BOUND_LEFT - 1;
        _relax_param.index_zero_begin = 0 + _DOF;
        _relax_param.index_zero_end = _DOF + _N_BOUND_LEFT - 1;
        _relax_param.index_mod_begin = _DOF + _N_BOUND_LEFT;
        _relax_param.index_mod_end = 2 * _DOF - 1;
        _relax_param.k_coeff = 2 * _DOF;

        _relax_param.index_column_off = _DOF - _N_BOUND_LEFT;

        _difeq(_relax_param, _param, _S);
        _reduce();
        _pinvs();

        // Back substitution
        _bksub();

        double err = 0.0;
        for (size_t j = 0; j < _DOF; j++) {
            double errj = 0, vmax = 0;
            double vz;
            int km = 0;
            for (size_t k = _relax_param.k_init; k <= _relax_param.k_final; k++) {
                vz = abs(_C[k][j][0]);
                if (isnormal(vz) && vz > vmax) {
                    vmax = vz;
                    km = k;
                }
                errj += vz;
            }
            err += errj / _scales[j];
            ermax[j] = _C[km][j][0] / _scales[j];
            kmax[j] = km;
        }
        err /= nvars;
        double fac = (err > _slowc ? _slowc / err : 1.0);
        for (size_t j = 0; j < _DOF; j++) {
            for (size_t k = _relax_param.k_init; k <= _relax_param.k_final; k++) {
                _Y[k][j] -= fac * _C[k][j][0];
            }
        }
        // cout<<"Iter.\tError\tFAC"<<endl;
        // cout<<iter<<"\t"<<err<<"\t"<<fac<<endl;
        // cout<<"Var\tKmax\tMax.Error"<<endl;
        // for (size_t id = 0; id < _DOF; id++)
        // {
        // cout<<id<<"\t"<<kmax[id]<<"\t"<<ermax[id]<<endl;
        // }
        // if (((iter+1)%100==0))
        // {
        //     DumpSolution("steps_caches/profile_step_"+to_string(iter+1)+".dat");
        // }

        if (err < _conv) {
            return true;
        }
    }
    cout << "TOO MANY TRIES" << endl;
    return false;
}

void Relaxation::_reduce() {
    //
    int k_pre = _relax_param.k - 1;  // The index of previous mesh point // kc

    int index_row_beg = _relax_param.index_row_begin;  // iz1
    int index_row_end = _relax_param.index_row_end;    // iz2

    int index_zero_beg = _relax_param.index_zero_begin;  // jz1
    int index_zero_end = _relax_param.index_zero_end;    // jz2
    int index_mod_beg = _relax_param.index_mod_begin;    // jm1
    int index_mod_end = _relax_param.index_mod_end;      // jm2

    int index_coeff = _relax_param.k_coeff;  // jmf

    int index_C_row_beg = _DOF - _N_BOUND_LEFT;       // ic1
    int index_C_column_beg = 0;                       // jc1
    int index_C_column_final = _DOF - _N_BOUND_LEFT;  // jcf

    int loff = index_C_column_beg - index_mod_beg;  // The align the _C and the _S
    int ic = index_C_row_beg;
    double vx;

    for (size_t j = index_zero_beg; j <= index_zero_end; j++) {
        for (size_t l = index_mod_beg; l <= index_mod_end; l++) {
            vx = _C[k_pre][ic][l + loff];
            for (size_t i = index_row_beg; i <= index_row_end; i++) {
                _S[i][l] -= _S[i][j] * vx;
            }
        }
        vx = _C[k_pre][ic][index_C_column_final];
        for (size_t i = index_row_beg; i <= index_row_end; i++) {
            _S[i][index_coeff] -= _S[i][j] * vx;
        }
        ic += 1;
    }
}

void Relaxation::_pinvs() {
    int k = _relax_param.k;

    int index_row_beg = _relax_param.index_row_begin;  // ie1
    int index_row_end = _relax_param.index_row_end;    // ie2

    int index_column_beg = _relax_param.index_mod_begin;  // je1
    int index_column_end = index_column_beg + (index_row_end - index_row_beg);
    int index_store_beg = index_column_end + 1;

    int index_coeff = _relax_param.k_coeff;  // jsf

    vector<int> indxr(index_row_end + 1, 0);
    VD pscl(index_row_end + 1, 0);

    // Record the biggest element in each row;
    double big;
    for (size_t i = index_row_beg; i <= index_row_end; i++) {
        big = 0;
        for (size_t j = index_column_beg; j <= index_column_end; j++) {
            if (abs(_S[i][j]) > big) {
                big = abs(_S[i][j]);
            }
        }
        if (big == 0) {
            cout << "Singular Matrix 1 in PINVS" << endl;
            return;
        }
        pscl[i] = 1.0 / big;
        indxr[i] = 0;
    }

    double piv;
    double pivinv;
    double dum;
    int jp;
    int ipiv, jpiv;

    for (size_t id = index_row_beg; id <= index_row_end; id++) {
        piv = 0.0;
        for (size_t i = index_row_beg; i <= index_row_end; i++) {
            if (indxr[i] == 0) {
                big = 0.0;
                for (size_t j = index_column_beg; j <= index_column_end; j++) {
                    if (abs(_S[i][j]) > big) {
                        jp = j;
                        big = abs(_S[i][j]);
                    }
                }
                if (big * pscl[i] > piv) {
                    ipiv = i;
                    jpiv = jp;
                    piv = big * pscl[i];
                }
            }
        }
        if (_S[ipiv][jpiv] == 0.0) {
            cout << "Singular matrix 2 in PINVS" << endl;
            return;
        }
        indxr[ipiv] = jpiv;
        pivinv = 1.0 / _S[ipiv][jpiv];
        for (size_t j = index_column_beg; j <= index_coeff; j++) {
            _S[ipiv][j] *= pivinv;
        }
        _S[ipiv][jpiv] = 1.0;
        for (size_t i = index_row_beg; i <= index_row_end; i++) {
            if (indxr[i] != jpiv) {
                if (_S[i][jpiv]) {
                    dum = _S[i][jpiv];
                    for (size_t j = index_column_beg; j <= index_coeff; j++) {
                        _S[i][j] -= dum * _S[ipiv][j];
                    }
                    _S[i][jpiv] = 0.0;
                }
            }
        }
    }
    // jc1
    int jcoff = _relax_param.index_column_off - index_store_beg;
    int icoff = index_row_beg - index_column_beg;
    int irow;
    for (size_t i = index_row_beg; i <= index_row_end; i++) {
        irow = indxr[i] + icoff;
        for (size_t j = index_store_beg; j <= index_coeff; j++) {
            _C[k][irow][j + jcoff] = _S[i][j];
        }
    }
}

void Relaxation::_bksub() {
    double xx;
    int im = 0;
    int k_pre;
    int index_column_coeff = _DOF - _N_BOUND_LEFT;
    for (int k = _relax_param.k_final; k >= _relax_param.k_init; k--) {
        if (k == _relax_param.k_init) {
            im = _N_BOUND_RIGHT;
        }
        k_pre = k + 1;
        for (int j = 0; j < _N_BOUND_RIGHT; j++) {
            xx = _C[k_pre][j][index_column_coeff];
            for (int i = im; i < _DOF; i++) {
                _C[k][i][index_column_coeff] -= _C[k][i][j] * xx;
            }
        }
    }

    for (int k = _relax_param.k_init; k <= _relax_param.k_final; k++) {
        k_pre = k + 1;
        for (int i = 0; i < _N_BOUND_LEFT; i++) {
            _C[k][i][0] = _C[k][i + _N_BOUND_RIGHT][index_column_coeff];
        }
        for (int i = 0; i < _N_BOUND_RIGHT; i++) {
            _C[k][i + _N_BOUND_LEFT][0] = _C[k_pre][i][index_column_coeff];
        }
    }
}

void Relaxation::PrintSolution() {
    cout << "The Solution is:" << endl;
    cout << "x\t";
    for (size_t i = 0; i < _DOF; i++) {
        cout << "y_" << i << "\t";
    }
    cout << endl;
    for (size_t i = 0; i < _X.size(); i++) {
        cout << _X[i] << "\t";
        for (size_t j = 0; j < _DOF; j++) {
            cout << _Y[i][j] << "\t";
        }
        cout << endl;
    }
}
void Relaxation::DumpSolution(string filename) {
    ofstream output(filename.c_str());
    // output<<"The Solution is:"<<endl;
    output << "x\t";
    for (size_t i = 0; i < _DOF; i++) {
        output << "y_" << i << "\t";
    }
    output << endl;
    for (size_t i = 0; i < _X.size(); i++) {
        output << _X[i] << "\t";
        for (size_t j = 0; j < _DOF; j++) {
            output << _Y[i][j] << "\t";
        }
        output << endl;
    }
}

void Relaxation::PrintS() {
    cout << "S: " << endl;
    for (int i = 0; i < _DOF; ++i) {
        for (int j = 0; j < 2 * _DOF + 1; ++j) {
            cout << _S[i][j] << "  ";
        }
        cout << endl;
    }
    cout << endl;
}
void Relaxation::PrintC() {
    cout << "C: " << endl;
    for (size_t k = 0; k < _N_MeshPoints + 1; k++) {
        cout << "At " << k << endl;
        for (size_t i = 0; i < _DOF; i++) {
            for (size_t j = 0; j < _DOF - _N_BOUND_LEFT + 1; j++) {
                cout << _C[k][i][j] << "\t";
            }
            cout << endl;
        }
        cout << endl;
    }
    cout << endl;
}

void Relaxation::PrintC(int k) {
    cout << "C at k=" << k << ": " << endl;
    for (size_t i = 0; i < _DOF; i++) {
        for (size_t j = 0; j < _DOF - _N_BOUND_LEFT + 1; j++) {
            cout << _C[k][i][j] << "\t";
        }
        cout << endl;
    }
    cout << endl;
}
