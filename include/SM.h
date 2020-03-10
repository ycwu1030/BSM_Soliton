#ifndef SM_BASIC_H
#define SM_BASIC_H

class SM
{
public:
    SM();
    ~SM(){};

    double GetMW() {return MW;}
    double GetMW2() {return MW2;}
    double GetVEV() {return vev;}
    
protected:
    double MW;
    double MW2;
    double thetaW;  
    double alpha;
    double vev;     // (sqrt(2)GF)^(-0.5)
    double ee;      // ee = sqrt(4*Pi*alpha)
    double g_weak;  // g = ee/sw
    double gp_hyper; // g' = ee/cw
    double yt; // sqrt(2)mt/vev;
    
};
#endif //SM_BASIC_H
