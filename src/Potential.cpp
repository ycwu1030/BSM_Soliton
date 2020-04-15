#include "Potential.h"
#include <cmath>
#include <iomanip>
#include <fstream>
#include "Eigen/Dense"

using namespace std;
using namespace Eigen;

double Potential::GetTotalEnergy(VD x, VVD fields)
{
    double energy = 0;
    double V0 = V0_global();
    for (size_t i = 0; i < x.size()-1; i++)
    {
        double density = 0;
        double DeltaZ = x[i+1]-x[i];
        VD yaver = (fields[i+1] + fields[i])/2;
        VD ydif  = (fields[i+1] - fields[i]);
        density += (ydif/DeltaZ*ydif/DeltaZ)/2;//(fields[i+1]-fields[i])*(fields[i+1]-fields[i])/pow(DeltaZ,2)/2;
        density += Vtotal(yaver) - V0;
        energy += density*DeltaZ;
    }
    return energy;
}
void Potential::DumpFullSolution(VD x, VVD fields, std::string filename)
{
    double density;
    std::ofstream output(filename.c_str());
    double V0 = V0_global();
    output<<"z";
    for (size_t i = 0; i < fields.front().size(); i++)
    {
        output<<"\tphi"<<i;
    }
    output<<"\tdensity"<<endl;
    for (int i = 0; i < x.size()-1; i++)
    {
        double DeltaZ = x[i+1]-x[i];
        double xaver = (x[i+1]+x[i])/2;
        VD yaver = (fields[i+1]+fields[i])/2;
        VD ydif = (fields[i+1]-fields[i]);
        density=((ydif/DeltaZ*ydif/DeltaZ)/2 + Vtotal(yaver) - V0);
        output<<scientific<<setprecision(10)<<x[i]<<"\t"<<fields[i]<<"\t"<<density<<endl;
    }
    output<<scientific<<setprecision(10)<<x[x.size()-1]<<"\t"<<fields[fields.size()-1]<<"\t"<<density<<endl;
}
void Potential::DumpEnergyDensity(VD x, VVD fields,std::string filename)
{
    double density;
    std::ofstream output(filename.c_str());
    double V0 = V0_global();
    output<<"x\tdensity"<<std::endl;
    for (int i = 0; i < x.size()-1; i++)
    {
        double DeltaZ = x[i+1]-x[i];
        double xaver = (x[i+1]+x[i])/2;
        VD yaver = (fields[i+1]+fields[i])/2;
        VD ydif = (fields[i+1]-fields[i]);
        density=((ydif/DeltaZ*ydif/DeltaZ)/2 + Vtotal(yaver) - V0);
        output<<xaver<<"\t"<<density<<std::endl;
    }
}
double Potential::GetTension(VD x, VVD fields)
{
    double tension = 0;
    for (size_t i = 0; i < x.size()-1; i++)
    {
        double DeltaZ = x[i+1]-x[i];
        VD dfield_dZ = (fields[i+1]-fields[i])/DeltaZ;
        tension += dfield_dZ*dfield_dZ * DeltaZ;
    }
    return tension;
}
bool Potential::CheckHessian(VD field_values)
{
    MatrixXd HM(_Field_Dim,_Field_Dim);
    VVD HM_VVD = d2Vtotal(field_values);
    for (size_t i = 0; i < _Field_Dim; i++)
    {
        for (size_t j = 0; j < _Field_Dim; j++)
        {
            HM(i,j) = HM_VVD[i][j];
        }
    }
    SelfAdjointEigenSolver<MatrixXd> eigensolver(HM);
    if (eigensolver.info() != Success) return false;
#if DEBUG
    cout<<"Eigen Values: "<<endl;
    cout<<eigensolver.eigenvalues()<<endl;
    cout<<"--"<<endl;
#endif
    return ((eigensolver.eigenvalues()).array()>0).all();
}