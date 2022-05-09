#ifndef _TOY_Z6_H_
#define _TOY_Z6_H_

#include <vector>

#include "Basic_Model.h"

namespace ToyTest {
class ToyZ6 : public Basic_Model {
protected:
    double beta_;
    double v1_global_;
    double v2_global_;

public:
    ToyZ6();
    ~ToyZ6(){};

    void Set_potential_Parameters(double beta);
    virtual double Vtotal(VD field_values, double scale = 1);
    virtual VD dVtotal(VD field_values, double scale = 1);
    virtual VVD d2Vtotal(VD field_values, double scale = 1);
    virtual double V0_global(double scale = 1);
    virtual VD QuarticCoupling(VD field_values);

    virtual void FindLocalMinima();
    virtual void PrintParameters();

    void PrintMinimaInOrder();

    std::vector<std::vector<double> > vev_order_by_angle;
};

}  // namespace ToyTest

#endif  //_TOY_Z6_H_
