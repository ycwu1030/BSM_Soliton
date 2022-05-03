#ifndef _TOY_Z5_H_
#define _TOY_Z5_H_

#include <utility>
#include <vector>

#include "Basic_Model.h"
#include "CubicSolver.h"

class ToyZ5 : public Basic_Model {
private:
    double beta_;
    double v1_global_;
    double v2_global_;
    CubicSolver sol_;

public:
    ToyZ5();
    ~ToyZ5(){};

    void Set_Potential_Parameters(double beta);
    virtual double Vtotal(VD field_values, double scale = 1);
    virtual VD dVtotal(VD field_values, double scale = 1);
    virtual VVD d2Vtotal(VD field_values, double scale = 1);
    virtual double V0_global(double scale = 1);
    virtual VD QuarticCoupling(VD field_values);

    virtual void FindLocalMinima();
    virtual void PrintParameters();

    std::vector<std::vector<double> > vev_order_by_angle;
};

#endif  //_TOY_Z5_H_
