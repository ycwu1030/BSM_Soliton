#include <iostream>

#include "Relaxation.h"

using namespace std;
using namespace BSM_Soliton;

namespace BSM_Soliton {
void Reduce_to_Zero(const unsigned s_row_dim, const unsigned offset, const unsigned s_row_max, const unsigned s_col_beg,
                    const VVD &c, VVD &s);
void Gaussian_Elimination_with_Partial_Pivot(const unsigned r_beg, const unsigned c_beg, const unsigned dim, VVD &s,
                                             VVD &c);
}  // namespace BSM_Soliton

int main(int argc, char const *argv[]) {
    int dof = 5;
    int n_left = 3;
    {
        cout << "Test middle blocks" << endl;
        VVD C{{3, 4, 5}, {1, 2, 7}, {1, 2, 3}, {4, 6, 7}, {5, 8, 2}};
        VVD S{{1, 2, 3, 4, 5, 6, 7, 8, 9, 3, 4},
              {4, 5, 2, 2, 6, 7, 8, 4, 2, 5, 5},
              {3, 4, 5, 6, 7, 3, 2, 4, 5, 6, 7},
              {2, 3, 5, 6, 7, 8, 9, 0, 3, 2, 3},
              {3, 4, 5, 6, 3, 1, 3, 4, 5, 6, 2}};
        Reduce_to_Zero(dof, n_left, dof, 0, C, S);
        for (int i = 0; i < dof; i++) {
            for (int j = 0; j < dof * 2 + 1; j++) {
                cout << S[i][j] << " ";
            }
            cout << endl;
        }
        cout << endl;
        Gaussian_Elimination_with_Partial_Pivot(0, n_left, dof, S, C);
        for (int i = 0; i < dof; i++) {
            for (int j = 0; j < dof * 2 + 1; j++) {
                cout << S[i][j] << " ";
            }
            cout << endl;
        }

        cout << endl;
        for (int i = 0; i < dof; i++) {
            for (int j = 0; j < dof - n_left + 1; j++) {
                cout << C[i][j] << " ";
            }
            cout << endl;
        }
        cout << endl;
    }

    {
        cout << "Test initial block" << endl;
        VVD C{{0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
        VVD S{{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
              {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
              {0, 0, 0, 0, 0, 3, 2, 4, 5, 6, 7},
              {0, 0, 0, 0, 0, 8, 9, 7, 3, 2, 3},
              {0, 0, 0, 0, 0, 1, 3, 4, 5, 6, 2}};
        Gaussian_Elimination_with_Partial_Pivot(dof - n_left, dof, n_left, S, C);
        for (int i = 0; i < dof; i++) {
            for (int j = 0; j < dof * 2 + 1; j++) {
                cout << S[i][j] << " ";
            }
            cout << endl;
        }
        cout << endl;
        for (int i = 0; i < dof; i++) {
            for (int j = 0; j < dof - n_left + 1; j++) {
                cout << C[i][j] << " ";
            }
            cout << endl;
        }
        cout << endl;
    }

    {
        cout << "Test the final block " << endl;
        VVD C{{3, 4, 5}, {1, 2, 7}, {1, 2, 3}, {4, 6, 7}, {5, 8, 2}};
        VVD S{{0, 0, 0, 0, 0, 3, 2, 4, 5, 6, 7},
              {0, 0, 0, 0, 0, 8, 9, 7, 3, 2, 3},
              {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
              {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
              {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};
        Reduce_to_Zero(dof, n_left, dof - n_left, dof, C, S);
        for (int i = 0; i < dof; i++) {
            for (int j = 0; j < dof * 2 + 1; j++) {
                cout << S[i][j] << " ";
            }
            cout << endl;
        }
        cout << endl;
        Gaussian_Elimination_with_Partial_Pivot(0, n_left + dof, dof - n_left, S, C);
        for (int i = 0; i < dof; i++) {
            for (int j = 0; j < dof * 2 + 1; j++) {
                cout << S[i][j] << " ";
            }
            cout << endl;
        }
        cout << endl;
        for (int i = 0; i < dof; i++) {
            for (int j = 0; j < dof - n_left + 1; j++) {
                cout << C[i][j] << " ";
            }
            cout << endl;
        }
    }

    return 0;
}
