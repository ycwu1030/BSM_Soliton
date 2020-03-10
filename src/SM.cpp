#include "SM.h"
#include <cmath>
#include "Constants.h"
#include <iostream>

using namespace std;

SM::SM()
{
    alpha = 1.0/alpha1;
    ee = sqrt(4*Pi*alpha);
    vev = pow(sqrt(2)*GF,-0.5);
    double A = sqrt(Pi*alpha)*vev;
    thetaW = asin(2*A/MZ)/2.0;
    MW = A/sin(thetaW);
    MW2 = MW*MW;
    g_weak = ee/sin(thetaW);
    gp_hyper = ee/cos(thetaW);
    yt = sqrt(2)*MT/vev;
#ifdef DEBUG
    cout<<"A:  "<<A<<endl;
    cout<<"vev: "<<vev<<endl;
    cout<<"thetaW: "<<thetaW<<endl;
    cout<<"MW: "<<MW<<endl;
    cout<<"yt: "<<yt<<endl;
    cout<<"g_weak: "<<g_weak<<endl;
    cout<<"gp_hyper: "<<gp_hyper<<endl;
#endif
}