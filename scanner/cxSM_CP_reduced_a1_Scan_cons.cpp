#include "cxSM_CP_reduced_a1.h"
#include "PathDeformation.h"
#include <iostream>
#include <fstream>
using namespace std;

double RandomReal(double min, double max)
{
    return (rand()/(double)RAND_MAX)*(max-min)+min;
}
int main(int argc, char const *argv[])
{
    cxSM_CP_reduced_a1 model;
    VnD vtol = [&](VD field, double scale){
        return model.Vtotal(field,scale);
    };
    dVnD dvtol = [&](VD field, double scale){
        return model.dVtotal(field,scale);
    };

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
    output<<"ID\tMHL\tMHH\tMHA\tvs\ttheta1\ttheta2\ttheta3\tmu2\tb2\tlam\tdel2\td2\ta1\tb1\talpha\tUnitarity\tStability\tGlobal\tEnergy"<<endl;
    while ( NGOT < 10000 )
    {
            vs = RandomReal(10,100);
            MHH = RandomReal(10,100);
            vs *= 1000.0; // 100 TeV;
            MHH *= 1000.0;
            MHA = MHH + 100;
            theta1 = 0.0001;
            theta3 = 0.01; 
            good = model.Set_Physical_Parameters_vs_theta(vs,MHH,MHA,theta1,theta3);
            if (!good||!model.CheckStability()||!model.CheckUnitarity()||!model.CheckGlobalMinimum()) continue;
            model.GetVS(vsr,vsi);
            VD left = {model.GetVEV(),vsr,-abs(vsi)};
            VD right = {model.GetVEV(),vsr,abs(vsi)};
            VVD pts_init;
            pts_init.push_back(left);
            pts_init.push_back(right);
            KinknD sol = fullKink(pts_init,vtol,dvtol);
            ++NGOT;
            output<<NGOT<<"\t"<<model<<"\t"<<model.CheckUnitarity()<<"\t"<<model.CheckStability()<<"\t"<<model.CheckGlobalMinimum()<<"\t"<<model.GetTotalEnergy(sol.R,sol.Phi)<<endl;
    }

    return 0;
}