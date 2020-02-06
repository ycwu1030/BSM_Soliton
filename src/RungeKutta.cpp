/*
 * @Description  : 
 * @Author       : Yongcheng Wu
 * @Date         : 2020-02-04 14:24:16
 * @LastEditors  : Yongcheng Wu
 * @LastEditTime : 2020-02-06 10:49:36
 */
#include "RungeKutta.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

#define POW_GROW (-0.2)
#define POW_SHRINK (-0.25)
#define SAFETY (0.9)
#define MAXSTEPS (10000)
#define TINY (1e-30)

bool abs_compare(double i, double j) {return abs(i)<abs(j);}

RungeKutta::RungeKutta()
{
    SetDOF();
    _derivs = nullptr;
}
RungeKutta::RungeKutta(int DOF)
{
    SetDOF(DOF);
    _derivs = nullptr;
}
void RungeKutta::SetDOF(int DOF)
{
    _DOF = DOF;
}
void RungeKutta::SetBound(double x_begin, double x_end, VD BOUND)
{
    _x_begin = x_begin;
    _x_end = x_end;
    _BOUND_CONDITION = BOUND;
}
void RungeKutta::SetODE(F_ODE derivs)
{
    _derivs = derivs;
}
void RungeKutta::_RESET()
{
    _X.clear();
    _Y.clear();
    _dYdX.clear();
    _X.push_back(_x_begin);
    _Y.push_back(_BOUND_CONDITION);
}
void RungeKutta::_RK4_SingleStep(double X_CUR, VD Y_CUR, VD dY_CUR, double step_size, VD &Y_NEXT)
{
    double half_step = step_size/2;

// The first step, just use the current dY/dX 
    VD dY_Step1 = step_size*dY_CUR;

// The second step, half_step in x, half dY_Step1 in Y
    VD dY_Step2 = step_size*_derivs(X_CUR+half_step,Y_CUR+dY_Step1/2);

// The third step, half_step in x, half dY_Step2 in Y
    VD dY_Step3 = step_size*_derivs(X_CUR+half_step,Y_CUR+dY_Step2/2);

// The fourth step, full step in x, full dY_Step3 in Y
    VD dY_Step4 = step_size*_derivs(X_CUR+step_size,Y_CUR+dY_Step3);
    
    Y_NEXT = Y_CUR + dY_Step1/6 + dY_Step2/3 + dY_Step3/3 + dY_Step4/6;
}

void RungeKutta::_RKQC_SingleStep(double &X, VD &Y, VD dY, double step_size_guess, double eps, VD Y_Scale, double &step_size_did, double &step_size_further)
{
    double x_cache = X;
    double step_size = step_size_guess;
    double half_step_size;
    double error_max = 0;
    double min_step_size = 1e-5*step_size_guess; // Since we will adapt the step size according to the estimated error, we need to terminate such operate at some point, otherwise, we will stuck at here for some case.
    double max_step_size; // We will enlarge the step size, when we are within the precision. But we need to limit it within a reasonable range.
    double step_size_temp;

    VD y_cache = Y;
    VD dy_cache = dY;
    VD y_temp;
    VD Delta_y;
    VD error_temp;

    while (true)
    {
        // Take two half steps
        // cout<<"Step Size: "<<step_size<<endl;
        half_step_size = step_size/2.0;
        _RK4_SingleStep(x_cache,y_cache,dy_cache,half_step_size,y_temp);
        X = x_cache + half_step_size;
        dY = _derivs(X,y_temp);
        _RK4_SingleStep(X,y_temp,dY,half_step_size,Y);

        // Take one large step
        X = x_cache + step_size;
        _RK4_SingleStep(x_cache,y_cache,dy_cache,step_size,y_temp);
        Delta_y = Y - y_temp;
        error_temp = abs(Delta_y/Y_Scale);
        // error_temp = abs(y_temp);
        error_max = *max_element(error_temp.begin(),error_temp.end());
        error_max /= eps;
        if (error_max <= 1.0)
        {
            step_size_did = step_size;
            step_size_temp = SAFETY*step_size*exp(POW_GROW*log(error_max));
            max_step_size = 4*step_size_did;
            step_size_further = min(step_size_temp,max_step_size);
            break;
        }
        step_size_temp = step_size*SAFETY*exp(POW_SHRINK*log(error_max));
        if (abs(step_size_temp) < abs(min_step_size))
        {
            step_size_did = step_size;
            step_size_further = step_size;
            break;
        }
        step_size = step_size_temp;
    }
    Y = Y + Delta_y/15;
}

void RungeKutta::ODEINTEGRAL(double step_start,double eps)
{
    _RESET();
    double x = _X[0];
    VD y = _Y[0];
    VD dydx = _derivs(x,y);
    _dYdX.push_back(dydx);
    double step_size = (_x_end > _x_begin)? abs(step_start):-abs(step_start);
    double step_size_did;
    double step_size_next;
    for (size_t nstp = 0; nstp < MAXSTEPS; ++nstp)
    {
        _Y_SCALE = abs(y)+abs(dydx*step_size);
        for (size_t i = 0; i < _DOF; i++)
        {
            _Y_SCALE[i] = min(1.0,_Y_SCALE[i]);
        }
        // Check whether the step size is too big, and already pass the end point
        if ( (step_size > 0 && x+step_size > _x_end) || (step_size < 0 && x+step_size < _x_end) )
        {
            step_size = _x_end - x;
        }
        _RKQC_SingleStep(x,y,dydx,step_size,eps,_Y_SCALE,step_size_did,step_size_next);
        _X.push_back(x);
        _Y.push_back(y);
        dydx=_derivs(x,y);
        _dYdX.push_back(dydx);
        if ((x - _x_end)*(_x_end-_x_begin)>=0)
        {
            // We are finished;
            return;
        }
        // step_size = abs(step_size_next) > 1e-3*abs(_x_end-_x_begin)?step_size_next:1e-3*(_x_end-_x_begin);
        step_size = step_size_next;
        if (abs(step_size_next) > 1e-1*abs(_x_end-_x_begin))
        {
            step_size = 1e-1*(_x_end-_x_begin);
        }
        if (abs(step_size_next) < 1e-5*abs(_x_end-_x_begin))
        {
            step_size = 1e-5*(_x_end-_x_begin);
        }
    }
    cout<<"TAKE TOO MANY STEPS"<<endl;
}

void RungeKutta::PrintSolution()
{
    cout<<"The Solution is:"<<endl;
    cout<<"x\t";
    for (size_t i = 0; i < _DOF; i++)
    {
        cout<<"y_"<<i<<"\t";
    }
    for (size_t i = 0; i < _DOF; i++)
    {
        cout<<"dy_"<<i<<"/dx"<<"\t";
    }
    cout<<endl;
    for (size_t i = 0; i < _X.size(); i++)
    {
        cout<<_X[i]<<"\t";
        for (size_t j = 0; j < _DOF; j++)
        {
            cout<<_Y[i][j]<<"\t";
        }
        for (size_t j = 0; j < _DOF; j++)
        {
            cout<<_dYdX[i][j]<<"\t";
        }
        cout<<endl;
    }
}
void RungeKutta::DumpSolution(string filename)
{
    ofstream output(filename.c_str());
    output<<"The Solution is:"<<endl;
    output<<"x\t";
    for (size_t i = 0; i < _DOF; i++)
    {
        output<<"y_"<<i<<"\t";
    }
    for (size_t i = 0; i < _DOF; i++)
    {
        output<<"dy_"<<i<<"/dx"<<"\t";
    }
    output<<endl;
    for (size_t i = 0; i < _X.size(); i++)
    {
        output<<_X[i]<<"\t";
        for (size_t j = 0; j < _DOF; j++)
        {
            output<<_Y[i][j]<<"\t";
        }
        for (size_t j = 0; j < _DOF; j++)
        {
            output<<_dYdX[i][j]<<"\t";
        }
        output<<endl;
    }
}


void RungeKutta::SetBound(double x_begin, double x_end, VD BOUND_begin, VD BOUND_end, vector<bool> At_Begin_Q, vector<bool> At_End_Q)
{
    _x_begin = x_begin;
    _x_end = x_end;
    _BOUND_begin = BOUND_begin;
    _BOUND_end = BOUND_end;
    _BOUND_begin_Q = At_Begin_Q;
    _BOUND_end_Q = At_End_Q;
    _N_BOUND_begin = 0;
    _N_BOUND_end = 0;
    for (size_t i = 0; i < _DOF; i++)
    {
        if (_BOUND_begin_Q[i])
        {
            ++_N_BOUND_begin;
        }
        if (_BOUND_end_Q[i])
        {
            ++_N_BOUND_end;
        }
    }
    if (_N_BOUND_end != _DOF - _N_BOUND_begin) cout<<"Error in Boundary Condition"<<endl;
    // printf("DOF %d\nAt Begin: %d\nAt End: %d\n",_DOF,_N_BOUND_begin,_N_BOUND_end);
}
void RungeKutta::_Load(VD BOUND_GUESS)
{
    int j = 0;
    // cout<<"LOADING BOUNDARY"<<endl;
    _BOUND_CONDITION.clear();
    for (size_t i = 0; i < _DOF; i++)
    {
        if (_BOUND_begin_Q[i])
        {
            _BOUND_CONDITION.push_back(_BOUND_begin[i]);
        }
        else
        {
            _BOUND_CONDITION.push_back(BOUND_GUESS[j++]);
        }
        // cout<<"BOUND AT "<<i<<": "<<_BOUND_CONDITION[i]<<endl;
    }
}
VD RungeKutta::_Score(VD _Y_END)
{
    VD res;
    for (size_t i = 0; i < _DOF; i++)
    {
        if (_BOUND_end_Q[i])
        {
            res.push_back(_Y_END[i]-_BOUND_end[i]);
        }
    }
    return res;
}
void RungeKutta::_SHOOTING(VD Bound_Guess, VD delta_Bound, VD &Bound_further, VD &score_last)
{
    VectorXd F_score(_N_BOUND_end);
    MatrixXd dFdV(_N_BOUND_end,_N_BOUND_end);

    _Load(Bound_Guess);
    ODEINTEGRAL((_x_end-_x_begin)/100);
    VD score = _Score(_Y.back());
    VD score_temp;
    VD Bound_cache = Bound_Guess;
    double bound_cache;
    for (size_t i = 0; i < _N_BOUND_end; i++)
    {

        F_score(i) = -score[i];
        bound_cache = Bound_Guess[i];
        Bound_Guess[i] += delta_Bound[i];
        _Load(Bound_Guess);
        ODEINTEGRAL((_x_end-_x_begin)/100);
        score_temp = _Score(_Y.back());
        for (size_t j = 0; j < _N_BOUND_end; j++)
        {
            dFdV(j,i) = (score_temp[j]-score[j])/delta_Bound[i];
        }
        Bound_Guess[i] = bound_cache;
    }
    VectorXd dV = dFdV.colPivHouseholderQr().solve(F_score);
    Bound_further.clear();
    for (size_t i = 0; i < _N_BOUND_end; i++)
    {
        Bound_further.push_back(Bound_Guess[i]+dV[i]);
    }
    _Load(Bound_further);
    ODEINTEGRAL((_x_end-_x_begin)/100);
    score_last=_Score(_Y.back());
}
void RungeKutta::SHOOTING(VD BOUND_GUESS, VD delta_Bound, VD eps)
{
    bool Finished = false;
    VD BOUND_NEXT;
    VD scores;
    int TRIED = 0;
    while (TRIED < 1000)
    {
        ++TRIED;
        _SHOOTING(BOUND_GUESS,delta_Bound,BOUND_NEXT,scores);
        Finished = true;
        for (size_t i = 0; i < _N_BOUND_end; i++)
        {
            Finished *= (abs(scores[i])<abs(eps[i]));
        }
        if (Finished) break;
        BOUND_GUESS = BOUND_NEXT;
    }
    if (TRIED>=1000) cout<<"SHOOTING TOO MANY TIMES"<<endl;
}