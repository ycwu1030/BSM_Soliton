#include "cxSM_CP_reduced.h"
#include "DWSolver.h"
#include <iostream>
#include <fstream>
using namespace std;

double RandomReal(double min, double max)
{
    return (rand()/(double)RAND_MAX)*(max-min)+min;
}
int main(int argc, char const *argv[])
{
    cxSM_CP_reduced model;
    DWSolver solver(&model);
    solver.SetOverallScale(model.GetVEV());
    VD left(3); // For CP case, the DOF is 3
    VD right(3);
    double vs;
    double alpha;
    double vsr;
    double vsi;
    double MHH;
    double MHA;
    double theta1;
    double theta2;
    double theta3;
    srand (time(NULL));
    int NGOT = 0;
    int NTRIED = 0;
    ofstream output(argv[1]);
    VD X;
    VVD Y;
    bool good;
    VD vs_candi = {10,20,30,40,50,60,70,80,90,100};
    VD M2_candi = {10,20,30,40,50,60,70,80,90,100};
    for (int i = 0; i < 10; i++)
    {
        for (int j = 0; j < 10; j++)
        {
            vs = vs_candi[i]*1000.0; // 100 TeV;
            MHH = M2_candi[j]*1000.0;
            MHA = MHH + 100;
            theta1 = 0.001;
            theta3 = 0.0001; 
            good = model.Set_Physical_Parameters_vs_theta(vs,MHH,MHA,theta1,theta3);
            if (!good||!model.CheckStability()||!model.CheckUnitarity()||!model.CheckGlobalMinimum()) continue;
            model.GetVS(vsr,vsi);
            model.GetTheta(2,theta2);
            ++NGOT;
            output<<model<<"\t"<<model.CheckUnitarity()<<"\t"<<model.CheckStability()<<"\t"<<model.CheckGlobalMinimum()<<endl;
        }
    }

    // while (NGOT < 2*100)
    // {
    //     /* 
    //      * Second way to provide parameters:
    //     */
    //     ++NTRIED;
    //     // * Feel free to change the range of these parameters
    //     vs = RandomReal(0,200);
    //     alpha = RandomReal(-M_PI,M_PI);
    //     MHH = RandomReal(0,500);
    //     MHA = RandomReal(0,200);
    //     theta1 = RandomReal(0,0.5);
    //     vsr = vs*cos(alpha);
    //     vsi = vs*sin(alpha);
    //     good = model.Set_Physical_Parameters_vsr_vsi_theta(vsr,vsi,MHH,MHA,theta1);
    //     if (!good||!model.CheckStability()||!model.CheckUnitarity()||!model.CheckGlobalMinimum()) continue;
    //     // ! We have following three possible boundary at one end, choice any one of them!!!!!
    //     left = {model.GetVEV(),vsr,-vsi};
    //     // left = {model.GetVEV(),-vsr,vsi};
    //     // left = {model.GetVEV(),-vsr,-vsi};
    //     // ! At the other end, it is the vacuum we choose for EW vacuum.
    //     right = {model.GetVEV(),vsr,vsi};
    //     solver.SetBoundary(left,right);
    //     good = solver.Solve(X,Y);
    //     if (!good) continue;
    //     model.GetTheta(2,theta2);
    //     model.GetTheta(3,theta3);
    //     ++NGOT;
    //     output<<NGOT<<"\t"<<vsr<<"\t"<<vsi<<"\t"<<MHH<<"\t"<<MHA<<"\t"<<theta1<<"\t"<<theta2<<"\t"<<theta3<<"\t"<<model.GetTotalEnergy(X,Y)<<endl;
    // }

    return 0;
}