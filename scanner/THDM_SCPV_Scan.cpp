#include "THDM_SCPV1.h"
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
    THDM_SCPV1 model;
    VnD vtol = [&](VD field, double scale){
        return model.Vtotal(field,scale);
    };
    dVnD dvtol = [&](VD field, double scale){
        return model.dVtotal(field,scale);
    };

    double tb;
    double beta;
    double alpha;
    double MHH;
    double MHA;
    double MHpm;
    double alphab;
    double alphac;
    double theta;
    srand (time(NULL));
    int NGOT = 0;
    int NTRIED = 0;
    ofstream output(argv[1]);
    VD X;
    VVD Y;
    bool good;
    output<<"id\t"<<model.repr()<<"\tUnitarity\tStability\tGlobal\tSigma"<<endl;

    while (NGOT < 10000)
    {
        ++NTRIED;
        // cout<<"Try "<<NTRIED<<endl;
        MHH = RandomReal(150,1000);
        MHA = RandomReal(150,1000);
        MHpm = RandomReal(150,1000);
        tb = RandomReal(0.5,20);
        beta = atan(tb);
        alpha = beta-M_PI_2;
        alphac = pow(10,RandomReal(-10,-1));
        good = model.Set_Physical_Parameters(beta,MHH,MHA,MHpm,alpha,alphac);
        if (!good) {
            // cout<<NTRIED<<": Bad input"<<endl;
            continue;
        }
        // cout<<NTRIED<<": Good input"<<endl;
        if (!model.CheckStability()||!model.CheckUnitarity()||!model.CheckGlobalMinimum()) continue;
        cout<<NTRIED<<": Good Points"<<endl;
        cout<<"\t"<<beta<<"\t"<<MHH<<"\t"<<MHA<<"\t"<<MHpm<<"\t"<<alpha<<"\t"<<alphac<<endl;
        model.GetTheta(theta);
        double vev = model.GetVEV();
        double v1 = vev*cos(beta);
        double v2r = vev*sin(beta)*cos(theta);
        double v2i = vev*sin(beta)*sin(theta);
        VD left = {v1,v2r,-v2i};
        VD right = {v1,v2r,v2i};
        VVD pts_init;
        pts_init.push_back(left);
        pts_init.push_back(right);
        KinknD sol = fullKink(pts_init,vtol,dvtol);
        ++NGOT;
        output<<NGOT<<"\t"<<model<<"\t"<<model.CheckUnitarity()<<"\t"<<model.CheckStability()<<"\t"<<model.CheckGlobalMinimum()<<"\t"<<model.GetTotalEnergy(sol.R,sol.Phi)<<endl;
    }
    

    
    return 0;
}