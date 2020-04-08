#include "Relaxation.h"
#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;
Relaxation::Relaxation()
{
    SetDOF(1,1);
    SetMaxIteration();
    SetConvergeCriterion();
}
Relaxation::Relaxation(int DOF, int NB_Left, int MeshPoints)
{
    SetDOF(DOF,NB_Left,MeshPoints);
    SetMaxIteration();
    SetConvergeCriterion();
}
void Relaxation::SetMaxIteration(int itemax)
{
    _ITE_MAX = itemax;
}
void Relaxation::SetDOF(int DOF, int NB_Left, int MeshPoints)
{
    _DOF = DOF;
    _N_BOUND_LEFT = NB_Left;
    _N_BOUND_RIGHT = _DOF - _N_BOUND_LEFT;
    _N_MeshPoints = MeshPoints;
}
void Relaxation::SetConvergeCriterion(double slowc, double conv)
{
    _slowc = slowc;
    _conv = conv;
}
void Relaxation::SetBoundary(double x_begin, double x_end, VD y_begin, VD y_end)
{
    _X_BEGIN = x_begin;
    _X_END = x_end;
    _Y_BEGIN = y_begin;
    _Y_END = y_end;
}
void Relaxation::SetODESystem(DIFEQ difeq, void *param)
{
    _difeq = difeq;
    _param = param;
}
void Relaxation::_INIT()
{
    // ! Initiative the vectors
    _X.clear();
    _Y.clear();
    if (_X_Guess.size()!=0 && _Y_Guess.size() == _X_Guess.size())
    {
        _X = _X_Guess;
        _Y = _Y_Guess;
    }
    else
    {
        double x_stepsize = (_X_END - _X_BEGIN)/(_N_MeshPoints-1);
        VD y_stepsize = (_Y_END - _Y_BEGIN)/(_N_MeshPoints-1);
        for (size_t i = 0; i < _N_MeshPoints; i++)
        {
            _X.push_back(_X_BEGIN+i*x_stepsize);
            _Y.push_back(_Y_BEGIN+i*y_stepsize);
        }   
    }
    _S.clear();
    _C.clear();
    _S = VVD(_DOF,VD(2*_DOF+1,0));
    _C = VVVD(_N_MeshPoints+1,VVD(_DOF,VD(_DOF-_N_BOUND_LEFT+1,0)));
}
bool Relaxation::SOLVDE()
{
    _INIT();
    _relax_param.k_init = 0;
    _relax_param.k_final = _N_MeshPoints-1;

    int nvars = _DOF*_N_MeshPoints;
    vector<int> kmax(_DOF,0);
    VD ermax(_DOF,0);

    for (size_t iter = 0; iter < _ITE_MAX; iter++)
    {
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

        _relax_param.k_coeff = 2*_DOF;
        _relax_param.index_row_begin = _DOF - _N_BOUND_LEFT;
        _relax_param.index_row_end = _DOF-1;
        _relax_param.index_mod_begin = _DOF;
        _relax_param.index_column_off = 0;

        _difeq(_relax_param,_param,_S);
        // cout<<"Iter: "<<iter<<"  "<<"S at left:"<<endl;
        // PrintS();
        _pinvs();
        
        // Intermidate points
        for (size_t k = _relax_param.k_init + 1; k <= _relax_param.k_final; k++)
        {
            _relax_param.k = k;

            _relax_param.x1 = _X[_relax_param.k-1];
            _relax_param.x2 = _X[_relax_param.k];
            _relax_param.y1 = _Y[_relax_param.k-1];
            _relax_param.y2 = _Y[_relax_param.k];

            _relax_param.index_row_begin = 0;
            _relax_param.index_row_end = _DOF - 1;
            _relax_param.index_zero_begin = 0;
            _relax_param.index_zero_end = _N_BOUND_LEFT - 1;
            _relax_param.index_mod_begin = _N_BOUND_LEFT;
            _relax_param.index_mod_end = _DOF - 1;
            _relax_param.k_coeff = 2*_DOF;

            _relax_param.index_column_off = 0;
            
            _difeq(_relax_param,_param,_S);
            _reduce();
            _pinvs();
            // cout<<"Iter: "<<iter<<"  C at "<<k<<endl;
            // PrintC();
        }

        // The final boundary condition
        _relax_param.k = _relax_param.k_final + 1;

        _relax_param.x1 = _X[_relax_param.k-1];
        _relax_param.x2 = _X[_relax_param.k-1];
        _relax_param.y1 = _Y[_relax_param.k-1];
        _relax_param.y2 = _Y[_relax_param.k-1];

        _relax_param.index_row_begin = 0;
        _relax_param.index_row_end = _DOF - _N_BOUND_LEFT - 1;
        _relax_param.index_zero_begin = 0 + _DOF;
        _relax_param.index_zero_end = _DOF + _N_BOUND_LEFT - 1;
        _relax_param.index_mod_begin = _DOF + _N_BOUND_LEFT;
        _relax_param.index_mod_end = 2*_DOF - 1;
        _relax_param.k_coeff = 2*_DOF;

        _relax_param.index_column_off = _DOF - _N_BOUND_LEFT;

        _difeq(_relax_param,_param,_S);
        _reduce();
        _pinvs();
    
        // Back substitution
        _bksub();
        
        double err = 0.0;
        for (size_t j = 0; j < _DOF; j++)
        {
            double errj=0,vmax=0;
            double vz;
            int km=0;
            for (size_t k = _relax_param.k_init; k <= _relax_param.k_final; k++)
            {
                vz = abs(_C[k][j][0]);
                if (isnormal(vz) && vz > vmax)
                {
                    vmax = vz;
                    km = k;
                }
                errj += vz;
            }
            err += errj/_scales[j];
            ermax[j] = _C[km][j][0]/_scales[j];
            kmax[j]=km;
        }
        err /= nvars;
        double fac = (err > _slowc? _slowc/err: 1.0);
        for (size_t j = 0; j < _DOF; j++)
        {
            for (size_t k = _relax_param.k_init; k <= _relax_param.k_final; k++)
            {
                _Y[k][j] -= fac*_C[k][j][0];
            }
        }
        // cout<<"Iter.\tError\tFAC"<<endl;
        // cout<<iter<<"\t"<<err<<"\t"<<fac<<endl;
        // cout<<"Var\tKmax\tMax.Error"<<endl;
        // for (size_t id = 0; id < _DOF; id++)
        // {
            // cout<<id<<"\t"<<kmax[id]<<"\t"<<ermax[id]<<endl;
        // }
        if (((iter+1)%100==0))
        {
            DumpSolution("steps_caches/profile_step_"+to_string(iter+1)+".dat");
        }
        
        if (err < _conv)
        {
            return true;
        }
    }
    cout<<"TOO MANY TRIES"<<endl;
    return false;
}

void Relaxation::_reduce()
{
    // 
    int k_pre = _relax_param.k - 1; // The index of previous mesh point // kc

    int index_row_beg = _relax_param.index_row_begin; // iz1
    int index_row_end = _relax_param.index_row_end; // iz2

    int index_zero_beg = _relax_param.index_zero_begin; // jz1
    int index_zero_end = _relax_param.index_zero_end; // jz2
    int index_mod_beg = _relax_param.index_mod_begin; // jm1
    int index_mod_end = _relax_param.index_mod_end; // jm2

    int index_coeff = _relax_param.k_coeff; // jmf

    int index_C_row_beg = _DOF - _N_BOUND_LEFT; // ic1
    int index_C_column_beg = 0; // jc1
    int index_C_column_final = _DOF - _N_BOUND_LEFT; // jcf

    int loff = index_C_column_beg - index_mod_beg; // The align the _C and the _S
    int ic = index_C_row_beg;
    double vx;

    for (size_t j = index_zero_beg; j <= index_zero_end; j++)
    {
        for (size_t l = index_mod_beg; l <= index_mod_end; l++)
        {
            vx = _C[k_pre][ic][l+loff];
            for (size_t i = index_row_beg; i <= index_row_end; i++)
            {
                _S[i][l] -= _S[i][j]*vx;
            }
        }
        vx = _C[k_pre][ic][index_C_column_final];
        for (size_t i = index_row_beg; i <= index_row_end; i++)
        {
            _S[i][index_coeff] -= _S[i][j]*vx;
        }
        ic += 1;
    }
}

void Relaxation::_pinvs()
{
    int k = _relax_param.k;
 
    int index_row_beg = _relax_param.index_row_begin; // ie1
    int index_row_end = _relax_param.index_row_end; // ie2

    int index_column_beg = _relax_param.index_mod_begin; // je1
    int index_column_end = index_column_beg + (index_row_end - index_row_beg);
    int index_store_beg = index_column_end + 1;

    int index_coeff = _relax_param.k_coeff; // jsf

    vector<int> indxr(index_row_end+1,0);
    VD pscl(index_row_end+1,0);


    // Record the biggest element in each row;
    double big;
    for (size_t i = index_row_beg; i <= index_row_end; i++)
    {
        big = 0;
        for (size_t j = index_column_beg; j <= index_column_end; j++)
        {
            if (abs(_S[i][j])>big)
            {
                big = abs(_S[i][j]);
            }   
        }
        if (big == 0)
        {
            cout<<"Singular Matrix 1 in PINVS"<<endl;
            return;
        }
        pscl[i]=1.0/big;
        indxr[i]=0;
    }
    
    
    double piv;
    double pivinv;
    double dum;
    int jp;
    int ipiv,jpiv;

    for (size_t id = index_row_beg; id <= index_row_end; id++)
    {
        piv = 0.0;
        for (size_t i = index_row_beg; i <= index_row_end; i++)
        {
            if (indxr[i]==0)
            {
                big = 0.0;
                for (size_t j = index_column_beg; j <= index_column_end; j++)
                {
                    if (abs(_S[i][j]) > big)
                    {
                        jp = j;
                        big = abs(_S[i][j]);
                    }
                }
                if (big*pscl[i] > piv)
                {
                    ipiv = i;
                    jpiv = jp;
                    piv = big*pscl[i];
                }
            }
        }
        if (_S[ipiv][jpiv] == 0.0)
        {
            cout<<"Singular matrix 2 in PINVS"<<endl;
            return;
        }
        indxr[ipiv] = jpiv;
        pivinv = 1.0/_S[ipiv][jpiv];
        for (size_t j = index_column_beg; j <= index_coeff; j++)
        {
            _S[ipiv][j] *= pivinv;
        }
        _S[ipiv][jpiv] = 1.0;
        for (size_t i = index_row_beg; i <= index_row_end; i++)
        {
            if (indxr[i] != jpiv)
            {
                if (_S[i][jpiv])
                {
                    dum = _S[i][jpiv];
                    for (size_t j = index_column_beg; j <= index_coeff; j++)
                    {
                        _S[i][j] -= dum*_S[ipiv][j];
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
    for (size_t i = index_row_beg; i <= index_row_end; i++)
    {
        irow = indxr[i] + icoff;
        for (size_t j = index_store_beg; j <= index_coeff; j++)
        {
            _C[k][irow][j+jcoff] = _S[i][j];
        }
    }
}

void Relaxation::_bksub()
{

    double xx;
    int im = 0;
    int k_pre;
    int index_column_coeff = _DOF - _N_BOUND_LEFT;
    for (int k = _relax_param.k_final; k >= _relax_param.k_init; k--)
    {
        if (k == _relax_param.k_init)
        {
            im = _N_BOUND_RIGHT;
        }
        k_pre = k+1;
        for (int j = 0; j < _N_BOUND_RIGHT; j++)
        {
            xx = _C[k_pre][j][index_column_coeff];
            for (int i = im; i < _DOF; i++)
            {
                _C[k][i][index_column_coeff] -= _C[k][i][j]*xx;
            }
        }
    }

    for (int k = _relax_param.k_init; k <= _relax_param.k_final; k++)
    {
        k_pre = k+1;
        for (int i = 0; i < _N_BOUND_LEFT; i++)
        {
            _C[k][i][0] = _C[k][i+_N_BOUND_RIGHT][index_column_coeff];
        }
        for (int i = 0; i < _N_BOUND_RIGHT; i++)
        {
            _C[k][i+_N_BOUND_LEFT][0] = _C[k_pre][i][index_column_coeff];
        }
    }
}

void Relaxation::PrintSolution()
{
    cout<<"The Solution is:"<<endl;
    cout<<"x\t";
    for (size_t i = 0; i < _DOF; i++)
    {
        cout<<"y_"<<i<<"\t";
    }
    cout<<endl;
    for (size_t i = 0; i < _X.size(); i++)
    {
        cout<<_X[i]<<"\t";
        for (size_t j = 0; j < _DOF; j++)
        {
            cout<<_Y[i][j]<<"\t";
        }
        cout<<endl;
    }
}
void Relaxation::DumpSolution(string filename)
{
    ofstream output(filename.c_str());
    // output<<"The Solution is:"<<endl;
    output<<"x\t";
    for (size_t i = 0; i < _DOF; i++)
    {
        output<<"y_"<<i<<"\t";
    }
    output<<endl;
    for (size_t i = 0; i < _X.size(); i++)
    {
        output<<_X[i]<<"\t";
        for (size_t j = 0; j < _DOF; j++)
        {
            output<<_Y[i][j]<<"\t";
        }
        output<<endl;
    }
}

void Relaxation::PrintS()
{
    cout<<"S: "<<endl;
    for (int i = 0; i < _DOF; ++i)
    {
        for (int j = 0; j < 2*_DOF+1; ++j)
        {
            cout<<_S[i][j]<<"  ";
        }
        cout<<endl;
    }
    cout<<endl;
}
void Relaxation::PrintC()
{
    cout<<"C: "<<endl;
    for (size_t k = 0; k < _N_MeshPoints + 1; k++)
    {
        cout<<"At "<<k<<endl;
        for (size_t i = 0; i < _DOF; i++)
        {
            for (size_t j = 0; j < _DOF-_N_BOUND_LEFT+1; j++)
            {
                cout<<_C[k][i][j]<<"\t";
            }
            cout<<endl;
        }
        cout<<endl;
    }
    cout<<endl;
}

void Relaxation::PrintC(int k)
{
    cout<<"C at k="<<k<<": "<<endl;
    for (size_t i = 0; i < _DOF; i++)
    {
        for (size_t j = 0; j < _DOF-_N_BOUND_LEFT+1; j++)
        {
            cout<<_C[k][i][j]<<"\t";
        }
        cout<<endl;
    }
    cout<<endl;
}