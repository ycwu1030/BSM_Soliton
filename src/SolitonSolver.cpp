/*
 * @Description  : 
 * @Author       : Yongcheng Wu
 * @Date         : 2020-01-28 12:46:03
 * @LastEditors  : Yongcheng Wu
 * @LastEditTime : 2020-02-04 12:09:14
 */
#include <iostream>
#include <fstream>
#include <cmath>
#include "SolitonSolver.h"

using namespace std;
void SolitonSolver::SetDimension(int NFieldDim)
{
    _NFieldDim = NFieldDim;
}
void SolitonSolver::SetGridPoints(int NGrid, int MAXROUNDS)
{
    _NGrid = NGrid;
    _NPoint = _NGrid+1;
    _MAXROUNDS = MAXROUNDS;
}
void SolitonSolver::SetPotentials(ScalarFunction potential, dScalarFunction dpotential)
{
    _potential = potential;
    _dpotential = dpotential;
}
void SolitonSolver::SetParam(void *param)
{
    _param = param;
}
void SolitonSolver::SetBoundary(VD left, VD right)
{
    _bound_left = left;
    _bound_right = right;
}
void SolitonSolver::SetScale(double scaling)
{
    _Scaling = scaling;
}
void SolitonSolver::SetV0Global(double V0global)
{
    _V0_global = V0global;
}
SolitonSolver::SolitonSolver(int NFieldDim, int NGrid, int MAXROUNDS)
{
    SetDimension(NFieldDim);
    SetGridPoints(NGrid,MAXROUNDS);
    SetScale();
}
double SolitonSolver::_Calc_EnergyDensity(VD point_beg, VD point_end)
{
    double density = 0;
    density += pow(_Scaling,4)*(point_end-point_beg)*(point_end-point_beg)/pow(_DeltaZ,2)/2;
    density += _potential((point_beg+point_end)/2*_Scaling,_param);
    return density - _V0_global;
}
void SolitonSolver::_Initiative()
{
    _GridPositions_NoDim.clear();
    _FieldValue_NoDim_Current.clear();
    _FieldValue_NoDim_Last.clear();
    _EnergyDensity.clear();
    _TotalEnergy_Last = 0;
    _TotalEnergy_Current = 0;
    _DeltaZ = (_RIGHT - _LEFT)/_NGrid;
    _DeltaT = _del_t*pow(_DeltaZ,2);
    VD FieldSteps = (_bound_right - _bound_left)/_NGrid;
    for (size_t i = 0; i < _NPoint; i++)
    {
        _GridPositions_NoDim.push_back(_LEFT + _DeltaZ*i);
        _FieldValue_NoDim_Current.push_back((_bound_left + FieldSteps*i)/_Scaling);
        _FieldValue_NoDim_Last.push_back((_bound_left + FieldSteps*i)/_Scaling);
    }
    for (size_t i = 0; i < _NGrid; i++)
    {
        _EnergyDensity.push_back(0);
    }
    // cout<<"Initial: "<<endl;
    // PrintSolitonSolution();
    // cout<<_TotalEnergy_Last<<"\t"<<_TotalEnergy_Current<<endl;
}
void SolitonSolver::_Iterating()
{
    int SUBROUNDS = 1000;
    for (size_t i = 0; i < _MAXROUNDS/SUBROUNDS; i++)
    {
        #if VERBOSE == 1
        cout<<"ROUND: "<<i<<"*"<<SUBROUNDS<<endl;
        #endif
        for (size_t k = 0; k < SUBROUNDS; k++)
        {
            // Store the results from last step
            for (size_t j = 0; j < _NPoint; j++)
            {
                _FieldValue_NoDim_Last[j] = _FieldValue_NoDim_Current[j];
            }
            // Iterate again;
            for (size_t j = 1; j < _NPoint-1; j++)
            {
                _FieldValue_NoDim_Current[j] = _FieldValue_NoDim_Last[j] + _del_t*(_FieldValue_NoDim_Last[j+1]+_FieldValue_NoDim_Last[j-1]-2*_FieldValue_NoDim_Last[j]);
                _FieldValue_NoDim_Current[j] = _FieldValue_NoDim_Current[j] - _DeltaT*_dpotential(_FieldValue_NoDim_Last[j]*_Scaling,_param)/pow(_Scaling,3);
            }
        }
        // Calculate the Energy_density/Energy, checking the improvement, if too small, jump out of the loop
        _TotalEnergy_Last = _TotalEnergy_Current;
        _TotalEnergy_Current = 0;
        for (size_t j = 0; j < _NGrid; j++)
        {
            _TotalEnergy_Current += _DeltaZ/_Scaling*(_EnergyDensity[j] = _Calc_EnergyDensity(_FieldValue_NoDim_Current[j],_FieldValue_NoDim_Current[j+1]));
        }
        if (_TotalEnergy_Last != 0)
        {
            double rel = abs(_TotalEnergy_Current-_TotalEnergy_Last)/_TotalEnergy_Last;
            #if VERBOSE == 1
            cout<<"Last: "<<_TotalEnergy_Last<<"  Current: "<<_TotalEnergy_Current<<"  rel_error: "<<rel<<endl;
            #endif
            if (rel < _energy_rel_error)
            {
                break;
            }
        }
    }
}
bool SolitonSolver::Solve()
{
    if (!_potential || !_dpotential || !_param)
    {
        cout<<"At least one of the potential, gradient or parameter is not set!"<<endl;
        return false;
    }
    _Initiative();
    _Iterating();
}
void SolitonSolver::PrintSolitonSolution()
{
    cout<<"The Soliton Solution: "<<endl;
    cout<<"Position_ID\tPosition\t";
    for (size_t i = 0; i < _NFieldDim; i++)
    {
        cout<<"Field_"<<i+1<<"\t";
    }
    cout<<"Energy_Density"<<endl;
    for (size_t i = 0; i < 10; i++)
    {
        cout<<_GridPositions_NoDim[i]<<"\t";
        cout<<_GridPositions_NoDim[i]/_Scaling<<"\t";
        for (size_t j = 0; j < _NFieldDim; j++)
        {
            cout<<_FieldValue_NoDim_Current[i][j]*_Scaling<<"\t";
        }
        if (i < _NGrid)
        {
            cout<<_EnergyDensity[i]<<endl;
        }
        else
        {
            cout<<endl;
        }
    }
    cout<<"......"<<endl;
    for (size_t i = _NPoint-10; i < _NPoint; i++)
    {
        cout<<_GridPositions_NoDim[i]<<"\t";
        cout<<_GridPositions_NoDim[i]/_Scaling<<"\t";
        for (size_t j = 0; j < _NFieldDim; j++)
        {
            cout<<_FieldValue_NoDim_Current[i][j]*_Scaling<<"\t";
        }
        if (i < _NGrid)
        {
            cout<<_EnergyDensity[i]<<endl;
        }
        else
        {
            cout<<endl;
        }
        
    }
    cout<<"Total Energy is: "<<_TotalEnergy_Current<<endl;
    cout<<"The Tension is: "<<GetTension()<<endl;
}
void SolitonSolver::DumpSolitonSolution(string filename)
{
    ofstream output(filename.c_str());
    output<<"Position_ID\tPosition\t";
    for (size_t i = 0; i < _NFieldDim; i++)
    {
        output<<"Field_"<<i+1<<"\t";
    }
    output<<"Energy_Density"<<endl;
    for (size_t i = 0; i < _NPoint; i++)
    {
        output<<_GridPositions_NoDim[i]<<"\t";
        output<<_GridPositions_NoDim[i]/_Scaling<<"\t";
        for (size_t j = 0; j < _NFieldDim; j++)
        {
            output<<_FieldValue_NoDim_Current[i][j]*_Scaling<<"\t";
        }
        if (i < _NGrid)
        {
            output<<_EnergyDensity[i]<<endl;
        }
        else
        {
            output<<endl;
        }
    }
    output<<"Total Energy is: "<<_TotalEnergy_Current<<endl;
    output<<"The Tension is: "<<GetTension()<<endl;
    output.close();
}
double SolitonSolver::GetTension()
{
    double tension = 0;
    for (size_t i = 0; i < _NGrid; i++)
    {
        VD dfield_dZ = (_FieldValue_NoDim_Current[i+1]-_FieldValue_NoDim_Current[i])/_DeltaZ;
        tension += dfield_dZ*dfield_dZ * _DeltaZ;
    }
    return tension*pow(_Scaling,3);
}