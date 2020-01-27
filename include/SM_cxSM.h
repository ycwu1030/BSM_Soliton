/*
 * @Description  : 
 * @Author       : Yongcheng Wu
 * @Date         : 2020-01-27 14:19:33
 * @LastEditors  : Yongcheng Wu
 * @LastEditTime : 2020-01-27 18:02:30
 */
#ifndef CXSM_H
#define CXSM_H

#include "CubicSolver.h"
#include <Eigen/Dense>
#include <stdio.h>
#include <vector>

#ifndef alpha1
#define alpha1 127.9
#endif

#ifndef MZ
#define MZ 91.1876
#endif

#ifndef MZ2
#define MZ2 (MZ*MZ)
#endif

#ifndef GF
#define GF 0.0000116637
#endif

#ifndef MT
#define MT 173.5
#endif

#ifndef MT2
#define MT2 (MT*MT)
#endif

#ifndef MHL
#define MHL 125.0
#endif

#ifndef MHL2
#define MHL2 (MHL*MHL)
#endif

class SM
{
public:
    SM();
    ~SM(){};

    double MW;
    double thetaW;  
    double alpha;
    double vev;     // (sqrt(2)GF)^(-0.5)
    double ee;      // ee = sqrt(4*Pi*alpha)
    double g_weak;  // g = ee/sw
    double gp_hyper; // g' = ee/cw
    double yt; // sqrt(2)mt/vev;

    double MW2;
    
};


class CXSM:public SM
{
public:
    CXSM();
    CXSM(double VSin, double MHHin, double MHAin, double thetain);
    ~CXSM(){};

    void SetInput(double VSin=50, double MHHin=400, double MHAin=0, double thetain=0.1);

// The Potential
    double Vtot(double phiH, double phiS, double phiA);

//Hessian Matrix
    void GetHessian(double phiH, double phiS, double phiA);
    bool CheckHessianMatrix(double phiH, double phiS, double phiA);
    
//Leading Order Stability and Unitarity
    bool CheckStability();
    bool CheckUnitarity(double MAX=0.5);
    bool CheckGlobalMinimum();

//Find LO local minimum position, used as starting points for critical point finding.
    void SolveCubicEquation(double A[4],double *results, int &NSolution);
    void FindLocalMinima();
    void PrintLocalMinima();
    void PrintPotentialParameter();

// Potential Parameters:
    double mu2;
    double lam;
    double del2;
    double b2;
    double d2;
    double del1;

// Physical Parameters:
    double MHH, MHH2;
    double MHA, MHA2;
    double VS;
    double theta;
    // double MHL;
    // double vev;

// Hessian Matrix:
    Eigen::Matrix3d HessianMatrix;

// Local Extreme @ T=0
    static const int NExtremeMax = 6;
    CubicSolver Solver;
    int NLocalExtreme;
    int NLocalMinima;
    int IndexInput;
    int MinimaIndex[NExtremeMax];
    double localH[NExtremeMax];
    double localS[NExtremeMax];
    double localA[NExtremeMax];
    bool LocalMinimaQ[NExtremeMax];
    double Vtotal[NExtremeMax];
    bool Solved;
};


#endif
